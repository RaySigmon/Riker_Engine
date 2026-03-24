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

"""Pipeline execution wrapper for the web UI.

Runs the Riker Engine pipeline in a background thread, capturing log output
and forwarding it through a queue for WebSocket streaming.
"""

import json
import logging
import shutil
import threading
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from queue import Queue

import yaml

RUNS_DIR = Path.home() / ".riker" / "runs"

# Suppress noisy warnings from UMAP/sklearn in log output
_FILTERED_PATTERNS = [
    "FutureWarning",
    "n_jobs value 1 overridden",
    "UserWarning: n_jobs",
    "numba.core.ssa",
    "list.append failed unexpectedly",
]


@dataclass
class RunState:
    """Tracks the state of a pipeline run."""

    run_id: str
    status: str = "pending"  # pending, running, complete, error
    current_phase: int = 0
    output_dir: Path = field(default_factory=Path)
    error: str = ""
    start_time: float = 0.0
    end_time: float = 0.0
    log_queue: Queue = field(default_factory=Queue)
    summary: dict = field(default_factory=dict)


# Global registry of runs
_runs: dict[str, RunState] = {}


class QueueLogHandler(logging.Handler):
    """Logging handler that forwards records to a queue for WebSocket streaming."""

    def __init__(self, queue: Queue):
        super().__init__()
        self.queue = queue

    def emit(self, record: logging.LogRecord) -> None:
        msg = self.format(record)
        # Filter noisy warnings
        if any(p in msg for p in _FILTERED_PATTERNS):
            return
        # Detect phase transitions
        if "PHASE " in msg and ":" in msg:
            for i in range(1, 7):
                if f"PHASE {i}:" in msg:
                    self.queue.put(json.dumps({"type": "phase", "phase": i, "name": msg.split(":", 1)[-1].strip()}))
                    break
        self.queue.put(json.dumps({"type": "log", "message": msg}))


def _normalize_config(config_dict: dict, output_dir: str) -> dict:
    """Translate UI JSON config into the engine's expected YAML format."""
    normalized = {
        "condition": config_dict.get("condition", ""),
        "seed_genes": config_dict.get("seed_genes", config_dict.get("seed_genes_path", "")),
        "hgnc_path": config_dict.get("hgnc_path", "auto"),
        "output_dir": output_dir,
    }

    # Datasets — translate UI keys to engine keys
    datasets = []
    for ds in config_dict.get("datasets", []):
        entry = {
            "id": ds.get("id", ds.get("datasetId", "")),
            "series_matrix": ds.get("series_matrix", ds.get("matrix_path", "")),
            "platform": ds.get("platform", ds.get("platform_path", "")),
            "role": ds.get("role", "discovery"),
            "tissue": ds.get("tissue", "brain"),
        }
        pf = ds.get("phenotype_field")
        if pf:
            entry["phenotype_field"] = pf
        cv = ds.get("case_values")
        if cv:
            entry["case_values"] = cv if isinstance(cv, list) else [cv]
        ctrl = ds.get("control_values")
        if ctrl:
            entry["control_values"] = ctrl if isinstance(ctrl, list) else [ctrl]
        datasets.append(entry)
    normalized["datasets"] = datasets

    # Phase parameters — accept both nested and flat formats
    params = config_dict.get("parameters", {})
    phase1 = config_dict.get("phase1", {})
    normalized["phase1"] = {
        "p_threshold": params.get("phase1_p_threshold", phase1.get("p_threshold", 0.05)),
        "min_datasets": params.get("phase1_min_datasets", phase1.get("min_datasets", 2)),
    }
    phase3 = config_dict.get("phase3", {})
    if phase3:
        normalized["phase3"] = phase3
    phase4 = config_dict.get("phase4", {})
    normalized["phase4"] = {
        "n_permutations": params.get("phase4_n_permutations", phase4.get("n_permutations", 10000)),
        "seed": phase4.get("seed", 42),
    }

    return normalized


def create_run(config_dict: dict) -> str:
    """Create a new run, write config, return run_id."""
    run_id = uuid.uuid4().hex[:12]
    run_dir = RUNS_DIR / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    output_dir = run_dir / "output"
    output_dir.mkdir(exist_ok=True)

    # Normalize UI config to engine format and write YAML
    normalized = _normalize_config(config_dict, str(output_dir))
    config_path = run_dir / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(normalized, f, default_flow_style=False)

    state = RunState(run_id=run_id, output_dir=output_dir)
    _runs[run_id] = state
    return run_id


def get_run(run_id: str) -> RunState | None:
    """Get run state by ID."""
    return _runs.get(run_id)


def start_run(run_id: str) -> None:
    """Start the pipeline in a background thread."""
    state = _runs.get(run_id)
    if not state:
        return

    thread = threading.Thread(target=_run_pipeline, args=(run_id,), daemon=True)
    thread.start()


def _run_pipeline(run_id: str) -> None:
    """Execute the pipeline (runs in background thread)."""
    state = _runs[run_id]
    state.status = "running"
    state.start_time = time.time()

    # Set up log capture
    handler = QueueLogHandler(state.log_queue)
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s", datefmt="%H:%M:%S"))
    handler.setLevel(logging.INFO)

    riker_logger = logging.getLogger("riker")
    riker_logger.addHandler(handler)
    riker_logger.setLevel(logging.INFO)

    config_path = RUNS_DIR / run_id / "config.yaml"

    try:
        from riker.config import load_config
        from riker.cli import cmd_run

        # Build a minimal args object that cmd_run expects
        class Args:
            config = str(config_path)
            verbose = False

        result = cmd_run(Args())

        state.end_time = time.time()

        if result == 0:
            state.status = "complete"
            # Load summary
            summary_path = state.output_dir / "pipeline_summary.json"
            if summary_path.exists():
                with open(summary_path) as f:
                    state.summary = json.load(f)
            state.log_queue.put(json.dumps({
                "type": "complete",
                "summary": state.summary,
            }))
        else:
            state.status = "error"
            state.error = "Pipeline returned non-zero exit code"
            state.log_queue.put(json.dumps({"type": "error", "message": state.error}))

    except Exception as e:
        state.end_time = time.time()
        state.status = "error"
        state.error = str(e)
        state.log_queue.put(json.dumps({"type": "error", "message": str(e)}))

    finally:
        riker_logger.removeHandler(handler)


def get_results(run_id: str) -> dict:
    """Load all result files for a completed run."""
    state = _runs.get(run_id)
    if not state or state.status != "complete":
        return {}

    results = {"summary": state.summary}
    output_dir = state.output_dir

    # Load CSV files as lists of dicts
    import csv
    for name in ["phase1_study_genes", "phase4_core_genes", "phase4_all_levels",
                  "phase5_verdicts", "phase6_meta_analysis"]:
        path = output_dir / f"{name}.csv"
        if path.exists():
            with open(path) as f:
                results[name] = list(csv.DictReader(f))

    # Load QC report
    qc_path = output_dir / "qc_report.json"
    if qc_path.exists():
        with open(qc_path) as f:
            results["qc_report"] = json.load(f)

    return results


def list_result_files(run_id: str) -> list[str]:
    """List downloadable result files."""
    state = _runs.get(run_id)
    if not state:
        return []
    return [f.name for f in state.output_dir.iterdir() if f.is_file()]


def get_result_file_path(run_id: str, filename: str) -> Path | None:
    """Get the full path to a result file."""
    state = _runs.get(run_id)
    if not state:
        return None
    path = state.output_dir / filename
    if path.exists() and path.is_file():
        return path
    return None


def create_results_zip(run_id: str) -> Path | None:
    """Create a ZIP of all results and return the path."""
    state = _runs.get(run_id)
    if not state:
        return None
    zip_path = RUNS_DIR / run_id / "results"
    shutil.make_archive(str(zip_path), "zip", str(state.output_dir))
    return Path(f"{zip_path}.zip")


def cleanup_old_runs(max_age_hours: int = 24) -> int:
    """Remove runs older than max_age_hours. Returns count removed."""
    if not RUNS_DIR.exists():
        return 0
    removed = 0
    cutoff = time.time() - (max_age_hours * 3600)
    for run_dir in RUNS_DIR.iterdir():
        if run_dir.is_dir() and run_dir.stat().st_mtime < cutoff:
            run_id = run_dir.name
            if run_id in _runs:
                del _runs[run_id]
            shutil.rmtree(run_dir, ignore_errors=True)
            removed += 1
    return removed
