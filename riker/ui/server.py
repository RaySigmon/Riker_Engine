# Riker Engine - Condition-Agnostic Transcriptomics Pipeline
# Copyright (C) 2024-2026 Ray Sigmon
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""FastAPI web server for the Riker Engine UI.

Provides REST endpoints for file upload, dataset inspection, config validation,
pipeline execution, and result retrieval. WebSocket endpoint streams real-time
progress during pipeline runs.
"""

import asyncio
import csv
import gzip
import io
import json
import tempfile
from pathlib import Path
from queue import Empty

import yaml
from fastapi import FastAPI, File, UploadFile, WebSocket, WebSocketDisconnect
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from starlette.requests import Request
from jinja2 import Environment, FileSystemLoader

from riker.ui import runner

app = FastAPI(title="Riker Engine", version="0.3.0")

# Template directory
_TEMPLATE_DIR = Path(__file__).parent / "templates"
_jinja_env = Environment(loader=FileSystemLoader(str(_TEMPLATE_DIR)))

# Temp directory for uploads
_UPLOAD_DIR = Path(tempfile.mkdtemp(prefix="riker_uploads_"))


# ---------------------------------------------------------------------------
# Page
# ---------------------------------------------------------------------------

@app.get("/", response_class=HTMLResponse)
async def index():
    """Serve the single-page UI."""
    template = _jinja_env.get_template("index.html")
    return template.render()


# ---------------------------------------------------------------------------
# File Uploads
# ---------------------------------------------------------------------------

@app.post("/api/upload/seeds")
async def upload_seeds(file: UploadFile = File(...)):
    """Upload a seed gene CSV file."""
    dest = _UPLOAD_DIR / file.filename
    content = await file.read()
    dest.write_bytes(content)

    # Parse preview
    text = content.decode("utf-8", errors="replace")
    reader = csv.DictReader(io.StringIO(text))
    genes = []
    # Find the symbol column
    fieldnames = reader.fieldnames or []
    symbol_col = None
    for candidate in ["symbol", "gene_symbol", "Symbol", "Gene_Symbol", "gene", "Gene"]:
        if candidate in fieldnames:
            symbol_col = candidate
            break
    if symbol_col is None and fieldnames:
        symbol_col = fieldnames[0]

    for row in reader:
        val = row.get(symbol_col, "").strip()
        if val:
            genes.append(val)

    return {
        "path": str(dest),
        "gene_count": len(genes),
        "preview": genes[:10],
        "columns": fieldnames,
        "symbol_column": symbol_col,
    }


@app.post("/api/upload/dataset")
async def upload_dataset(file: UploadFile = File(...)):
    """Upload a GEO series matrix file (.txt.gz)."""
    dest = _UPLOAD_DIR / file.filename
    content = await file.read()
    dest.write_bytes(content)

    # Quick inspect for sample count
    info = _inspect_series_matrix(dest)
    return {
        "path": str(dest),
        "filename": file.filename,
        **info,
    }


@app.post("/api/upload/platform")
async def upload_platform(file: UploadFile = File(...)):
    """Upload a platform annotation file."""
    dest = _UPLOAD_DIR / file.filename
    content = await file.read()
    dest.write_bytes(content)
    return {"path": str(dest), "filename": file.filename}


# ---------------------------------------------------------------------------
# Dataset Inspection
# ---------------------------------------------------------------------------

@app.post("/api/inspect-dataset")
async def inspect_dataset(request: Request):
    """Inspect a dataset file and return metadata."""
    body = await request.json()
    path = Path(body["path"])
    if not path.exists():
        return JSONResponse({"error": "File not found"}, status_code=404)
    info = _inspect_series_matrix(path)
    return info


def _inspect_series_matrix(path: Path) -> dict:
    """Parse a series matrix for metadata."""
    try:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", errors="replace") as f:
            lines = f.readlines()
    except Exception as e:
        return {"error": str(e)}

    sample_count = 0
    phenotype_fields = []
    probe_format = ""
    value_range = ""

    # Collect unique characteristic rows
    char_data: dict[str, list[str]] = {}

    for line in lines:
        if line.startswith("!Sample_geo_accession"):
            sample_count = len(line.strip().split("\t")) - 1

        if line.startswith("!Sample_characteristics_ch"):
            parts = line.strip().split("\t")
            field_name = parts[0].strip()
            vals = [p.strip('"') for p in parts[1:]]
            unique_vals = sorted(set(vals))
            # Determine the key (e.g., "disease status" from "disease status: AD")
            if unique_vals:
                sample_key = unique_vals[0].split(":")[0].strip() if ":" in unique_vals[0] else field_name
                entry_key = f"{field_name} | {sample_key}"
                if entry_key not in char_data:
                    char_data[entry_key] = unique_vals

        if line.startswith("!Sample_source_name_ch1"):
            parts = line.strip().split("\t")
            vals = [p.strip('"') for p in parts[1:]]
            unique_vals = sorted(set(vals))
            char_data["Sample_source_name_ch1 | source"] = unique_vals

        if line.startswith("!Sample_title"):
            parts = line.strip().split("\t")
            vals = [p.strip('"') for p in parts[1:6]]
            unique_vals = sorted(set(vals))
            char_data["Sample_title | title"] = unique_vals[:10]

        if line.startswith('"ID_REF"') or (not line.startswith("!") and "\t" in line and not line.startswith('"ID_REF"')):
            if not probe_format and not line.startswith("!"):
                parts = line.strip().split("\t")
                probe_id = parts[0].strip('"')
                if probe_id and probe_id != "ID_REF":
                    if probe_id.startswith("ILMN_"):
                        probe_format = "Illumina (ILMN_*)"
                    elif probe_id.startswith("ENSG"):
                        probe_format = "Ensembl (ENSG*)"
                    elif "_at" in probe_id or "_s_at" in probe_id:
                        probe_format = "Affymetrix (*_at)"
                    elif probe_id.startswith("GI_"):
                        probe_format = "GenInfo (GI_*)"
                    elif probe_id.replace(".", "").isdigit():
                        probe_format = "Numeric"
                    else:
                        probe_format = f"Other ({probe_id[:15]})"
                    # Check value range
                    try:
                        vals_num = [float(v) for v in parts[1:6] if v.strip()]
                        if vals_num:
                            mn, mx = min(vals_num), max(vals_num)
                            if mx < 5:
                                value_range = f"Log-ratio ({mn:.2f} to {mx:.2f})"
                            elif mx < 20:
                                value_range = f"Log2-intensity ({mn:.1f} to {mx:.1f})"
                            elif mx > 100:
                                value_range = f"Raw intensity ({mn:.0f} to {mx:.0f})"
                            else:
                                value_range = f"{mn:.2f} to {mx:.2f}"
                    except (ValueError, IndexError):
                        pass

    # Build phenotype fields list
    for key, vals in char_data.items():
        parts = key.split(" | ")
        phenotype_fields.append({
            "field": parts[0],
            "label": parts[1] if len(parts) > 1 else parts[0],
            "values": vals[:30],  # cap at 30 unique values
        })

    return {
        "sample_count": sample_count,
        "probe_format": probe_format,
        "value_range": value_range,
        "phenotype_fields": phenotype_fields,
    }


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

@app.post("/api/validate")
async def validate_config(request: Request):
    """Validate a pipeline config."""
    body = await request.json()
    try:
        # Normalize UI keys to engine format, then write temp YAML
        normalized = runner._normalize_config(body, str(_UPLOAD_DIR / "validate_output"))
        tmp = _UPLOAD_DIR / "validate_config.yaml"
        with open(tmp, "w") as f:
            yaml.dump(normalized, f, default_flow_style=False)

        from riker.config import load_config
        config = load_config(str(tmp))
        return {
            "valid": True,
            "condition": config.condition,
            "n_datasets": len(config.datasets),
            "n_discovery": sum(1 for d in config.datasets if d.role == "discovery"),
            "n_replication": sum(1 for d in config.datasets if d.role == "replication"),
        }
    except Exception as e:
        return {"valid": False, "error": str(e)}


# ---------------------------------------------------------------------------
# Pipeline Execution
# ---------------------------------------------------------------------------

@app.post("/api/run")
async def start_pipeline(request: Request):
    """Start a pipeline run in the background."""
    body = await request.json()
    run_id = runner.create_run(body)
    runner.start_run(run_id)
    return {"run_id": run_id}


@app.websocket("/ws/progress/{run_id}")
async def ws_progress(websocket: WebSocket, run_id: str):
    """Stream pipeline progress via WebSocket."""
    await websocket.accept()
    state = runner.get_run(run_id)
    if not state:
        await websocket.send_text(json.dumps({"type": "error", "message": "Run not found"}))
        await websocket.close()
        return

    try:
        while True:
            # Drain the queue
            try:
                while True:
                    msg = state.log_queue.get_nowait()
                    await websocket.send_text(msg)
                    # Check if this was a terminal message
                    parsed = json.loads(msg)
                    if parsed.get("type") in ("complete", "error"):
                        await websocket.close()
                        return
            except Empty:
                pass

            # Check if run finished without queue message
            if state.status in ("complete", "error"):
                if state.status == "complete":
                    await websocket.send_text(json.dumps({"type": "complete", "summary": state.summary}))
                else:
                    await websocket.send_text(json.dumps({"type": "error", "message": state.error}))
                await websocket.close()
                return

            await asyncio.sleep(0.3)

    except WebSocketDisconnect:
        pass


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

@app.get("/api/results/{run_id}")
async def get_results(run_id: str):
    """Get full results for a completed run."""
    results = runner.get_results(run_id)
    if not results:
        return JSONResponse({"error": "Run not found or not complete"}, status_code=404)
    return results


@app.get("/api/results/{run_id}/files")
async def list_files(run_id: str):
    """List downloadable result files."""
    files = runner.list_result_files(run_id)
    return {"files": files}


@app.get("/api/results/{run_id}/download/{filename}")
async def download_file(run_id: str, filename: str):
    """Download a specific result file."""
    path = runner.get_result_file_path(run_id, filename)
    if not path:
        return JSONResponse({"error": "File not found"}, status_code=404)
    return FileResponse(path, filename=filename)


@app.get("/api/results/{run_id}/download-all")
async def download_all(run_id: str):
    """Download all results as a ZIP file."""
    zip_path = runner.create_results_zip(run_id)
    if not zip_path:
        return JSONResponse({"error": "Run not found"}, status_code=404)
    return FileResponse(zip_path, filename=f"riker_results_{run_id}.zip",
                        media_type="application/zip")


# ---------------------------------------------------------------------------
# Lifecycle
# ---------------------------------------------------------------------------

@app.on_event("startup")
async def startup():
    """Clean up old runs on startup."""
    runner.cleanup_old_runs()
