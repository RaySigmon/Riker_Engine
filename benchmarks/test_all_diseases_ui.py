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

"""Automated UI backend test for all 5 validated diseases.

Starts the FastAPI server, uploads seed genes and datasets for each disease,
validates config, runs the pipeline, and verifies results.

Usage:
    python benchmarks/test_all_diseases_ui.py [--diseases ASD,T2D,IBD,AD,BRCA]
"""

import asyncio
import json
import sys
import time
import threading
import requests
import argparse

# ---- Disease configurations (from validated YAML configs) ----

DISEASES = {
    "ASD": {
        "condition": "ASD",
        "seed_genes": "/home/kai001/asd_full_run/data/seed_genes/sfari_asd_genes_clean.csv",
        "hgnc_path": "/home/kai001/.riker/hgnc_complete_set.txt",
        "datasets": [
            {
                "id": "GSE28521", "role": "discovery", "tissue": "brain",
                "series_matrix": "/home/kai001/asd_full_run/data/bulk_geo/GSE28521_series_matrix.txt.gz",
                "platform": "/home/kai001/asd_full_run/data/bulk_geo/GPL6883_clean.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["disease status: autism"],
                "control_values": ["disease status: control"],
            },
            {
                "id": "GSE28475", "role": "discovery", "tissue": "brain",
                "series_matrix": "/home/kai001/asd_full_run/data/bulk_geo/GSE28475-GPL6883_fixed_series_matrix.txt.gz",
                "platform": "/home/kai001/asd_full_run/data/bulk_geo/GPL6883_clean.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["diagnosis: autism"],
                "control_values": ["diagnosis: control"],
            },
            {
                "id": "GSE64018", "role": "discovery", "tissue": "brain",
                "series_matrix": "/home/kai001/asd_full_run/data/bulk_geo/GSE64018_reconstructed_series_matrix.txt.gz",
                "platform": "/home/kai001/asd_full_run/data/bulk_geo/GPL11154_ensembl.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["diagnosis: asd"],
                "control_values": ["diagnosis: ctl"],
            },
            {
                "id": "GSE102741", "role": "replication", "tissue": "brain",
                "series_matrix": "/home/kai001/asd_full_run/data/bulk_geo/GSE102741_reconstructed_series_matrix.txt.gz",
                "platform": "/home/kai001/asd_full_run/data/bulk_geo/GPL11154_ensembl.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["disease status: autism spectrum disorder"],
                "control_values": ["disease status: healthy control"],
            },
            {
                "id": "GSE18123", "role": "replication", "tissue": "blood",
                "series_matrix": "/home/kai001/asd_full_run/data/bulk_geo/GSE18123-GPL6244_series_matrix.txt.gz",
                "platform": "/home/kai001/asd_full_run/data/bulk_geo/GPL6244_clean.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["diagnosis: autism", "diagnosis: pdd-nos", "diagnosis: asperger"],
                "control_values": ["diagnosis: control"],
            },
            {
                "id": "GSE26415", "role": "replication", "tissue": "blood",
                "series_matrix": "/home/kai001/asd_full_run/data/bulk_geo/GSE26415_series_matrix.txt.gz",
                "platform": "/home/kai001/asd_full_run/data/bulk_geo/GPL6480_clean.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["disease: autism spectrum disorder"],
                "control_values": ["disease: nonaustistic control"],
            },
            {
                "id": "GSE42133", "role": "replication", "tissue": "blood",
                "series_matrix": "/home/kai001/asd_full_run/data/bulk_geo/GSE42133_series_matrix.txt.gz",
                "platform": "/home/kai001/asd_full_run/data/bulk_geo/GPL10558_clean.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["dx (diagnosis): asd"],
                "control_values": ["dx (diagnosis): control"],
            },
        ],
    },
    "T2D": {
        "condition": "T2D",
        "seed_genes": "/home/kai001/t2d_validation/data/seed_genes/t2d_seed_genes.csv",
        "hgnc_path": "/home/kai001/.riker/hgnc_complete_set.txt",
        "datasets": [
            {
                "id": "GSE41762", "role": "discovery", "tissue": "islet",
                "series_matrix": "/home/kai001/t2d_validation/data/bulk_geo/GSE41762_series_matrix.txt.gz",
                "platform": "/home/kai001/t2d_validation/data/bulk_geo/GPL6244.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["status: diabetic donor"],
                "control_values": ["status: non-diabetic donor"],
            },
            {
                "id": "GSE25724", "role": "discovery", "tissue": "islet",
                "series_matrix": "/home/kai001/t2d_validation/data/bulk_geo/GSE25724_series_matrix.txt.gz",
                "platform": "/home/kai001/t2d_validation/data/bulk_geo/GPL96.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["disease state: type 2 diabetes"],
                "control_values": ["disease state: non-diabetic"],
            },
            {
                "id": "GSE20966", "role": "discovery", "tissue": "islet",
                "series_matrix": "/home/kai001/t2d_validation/data/bulk_geo/GSE20966_series_matrix.txt.gz",
                "platform": "/home/kai001/t2d_validation/data/bulk_geo/GPL1352.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["disease: type 2 diabetes"],
                "control_values": ["disease: non-diabetic control"],
            },
            {
                "id": "GSE86468", "role": "replication", "tissue": "islet",
                "series_matrix": "/home/kai001/t2d_validation/data/bulk_geo/GSE86468_reconstructed_series_matrix.txt.gz",
                "platform": "/home/kai001/t2d_validation/data/bulk_geo/ensembl_to_symbol.txt",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["disease: type 2 diabetic"],
                "control_values": ["disease: non-diabetic"],
            },
        ],
    },
    "IBD": {
        "condition": "IBD",
        "seed_genes": "/home/kai001/ibd_validation/data/seed_genes/ibd_seed_genes.csv",
        "hgnc_path": "/home/kai001/.riker/hgnc_complete_set.txt",
        "datasets": [
            {
                "id": "GSE75214", "role": "discovery", "tissue": "colon",
                "series_matrix": "/home/kai001/ibd_validation/data/bulk_geo/GSE75214_series_matrix.txt.gz",
                "platform": "/home/kai001/ibd_validation/data/bulk_geo/GPL6244.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["disease: cd", "disease: crohn", "disease: ulcerative colitis"],
                "control_values": ["disease: control"],
            },
            {
                "id": "GSE16879", "role": "discovery", "tissue": "colon",
                "series_matrix": "/home/kai001/ibd_validation/data/bulk_geo/GSE16879_series_matrix.txt.gz",
                "platform": "/home/kai001/ibd_validation/data/bulk_geo/GPL570.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["disease: cd", "disease: uc"],
                "control_values": ["disease: control"],
            },
            {
                "id": "GSE59071", "role": "discovery", "tissue": "colon",
                "series_matrix": "/home/kai001/ibd_validation/data/bulk_geo/GSE59071_series_matrix.txt.gz",
                "platform": "/home/kai001/ibd_validation/data/bulk_geo/GPL6244.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["disease: crohn", "disease: ulcerative colitis"],
                "control_values": ["disease: control"],
            },
            {
                "id": "GSE87466", "role": "replication", "tissue": "colon",
                "series_matrix": "/home/kai001/ibd_validation/data/bulk_geo/GSE87466_series_matrix.txt.gz",
                "platform": "/home/kai001/ibd_validation/data/bulk_geo/GPL13158.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["disease: ulcerative colitis"],
                "control_values": ["disease: normal"],
            },
            {
                "id": "GSE38713", "role": "replication", "tissue": "colon",
                "series_matrix": "/home/kai001/ibd_validation/data/bulk_geo/GSE38713_series_matrix.txt.gz",
                "platform": "/home/kai001/ibd_validation/data/bulk_geo/GPL570.annot",
                "phenotype_field": "Sample_source_name_ch1",
                "case_values": ["uc patient"],
                "control_values": ["non-inflammatory control"],
            },
            {
                "id": "GSE36807", "role": "replication", "tissue": "colon",
                "series_matrix": "/home/kai001/ibd_validation/data/bulk_geo/GSE36807_series_matrix.txt.gz",
                "platform": "/home/kai001/ibd_validation/data/bulk_geo/GPL570.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["diagnosis: crohn", "diagnosis: ulcerative colitis"],
                "control_values": ["diagnosis: healthy control"],
            },
        ],
    },
    "AD": {
        "condition": "alzheimers_disease",
        "seed_genes": "/home/kai001/ad_validation/data/seed_genes/ad_genes.csv",
        "hgnc_path": "/home/kai001/.riker/hgnc_complete_set.txt",
        "datasets": [
            {
                "id": "GSE33000", "role": "discovery", "tissue": "brain",
                "series_matrix": "/home/kai001/ad_validation/data/bulk_geo/GSE33000_AD_only_series_matrix.txt.gz",
                "platform": "/home/kai001/ad_validation/data/bulk_geo/GPL4372_clean.annot",
                "phenotype_field": "Sample_characteristics_ch2",
                "case_values": ["disease status: alzheimer"],
                "control_values": ["disease status: non-demented"],
            },
            {
                "id": "GSE44770", "role": "discovery", "tissue": "brain",
                "series_matrix": "/home/kai001/ad_validation/data/bulk_geo/GSE44770_series_matrix.txt.gz",
                "platform": "/home/kai001/ad_validation/data/bulk_geo/GPL4372_clean.annot",
                "phenotype_field": "Sample_characteristics_ch2",
                "case_values": ["disease: a"],
                "control_values": ["disease: n"],
            },
            {
                "id": "GSE118553", "role": "discovery", "tissue": "brain",
                "series_matrix": "/home/kai001/ad_validation/data/bulk_geo/GSE118553_TC_series_matrix.txt.gz",
                "platform": "/home/kai001/ad_validation/data/bulk_geo/GPL10558_clean.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["disease state: ad"],
                "control_values": ["disease state: control"],
            },
            {
                "id": "GSE5281", "role": "replication", "tissue": "brain",
                "series_matrix": "/home/kai001/ad_validation/data/bulk_geo/GSE5281_SFG_series_matrix.txt.gz",
                "platform": "/home/kai001/ad_validation/data/bulk_geo/GPL570_clean.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["disease state: alzheimer"],
                "control_values": ["disease state: normal"],
            },
            {
                "id": "GSE15222", "role": "replication", "tissue": "brain",
                "series_matrix": "/home/kai001/ad_validation/data/bulk_geo/GSE15222_series_matrix.txt.gz",
                "platform": "/home/kai001/ad_validation/data/bulk_geo/GPL2700_clean.annot",
                "phenotype_field": "Sample_source_name_ch1",
                "case_values": ["alzheimer"],
                "control_values": ["normal"],
            },
        ],
    },
    "BRCA": {
        "condition": "breast_cancer",
        "seed_genes": "/home/kai001/breast_cancer_test/seed_genes.csv",
        "hgnc_path": "/home/kai001/.riker/hgnc_complete_set.txt",
        "datasets": [
            {
                "id": "GSE10810", "role": "discovery", "tissue": "breast",
                "series_matrix": "/home/kai001/breast_cancer_test/datasets/GSE10810_series_matrix.txt.gz",
                "platform": "/home/kai001/breast_cancer_test/platforms/GPL570_clean.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["tumor (t) vs healthy (s): t"],
                "control_values": ["tumor (t) vs healthy (s): s"],
            },
            {
                "id": "GSE42568", "role": "discovery", "tissue": "breast",
                "series_matrix": "/home/kai001/breast_cancer_test/datasets/GSE42568_series_matrix.txt.gz",
                "platform": "/home/kai001/breast_cancer_test/platforms/GPL570_clean.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["tissue: breast cancer"],
                "control_values": ["tissue: normal breast"],
            },
            {
                "id": "GSE15852", "role": "discovery", "tissue": "breast",
                "series_matrix": "/home/kai001/breast_cancer_test/datasets/GSE15852_series_matrix.txt.gz",
                "platform": "/home/kai001/breast_cancer_test/platforms/GPL96_clean.annot",
                "phenotype_field": "Sample_source_name_ch1",
                "case_values": ["breast tumor tissue"],
                "control_values": ["normal breast tissue"],
            },
            {
                "id": "GSE45827", "role": "replication", "tissue": "breast",
                "series_matrix": "/home/kai001/breast_cancer_test/datasets/GSE45827_series_matrix.txt.gz",
                "platform": "/home/kai001/breast_cancer_test/platforms/GPL570_clean.annot",
                "phenotype_field": "Sample_source_name_ch1",
                "case_values": ["tumor"],
                "control_values": ["normal"],
            },
            {
                "id": "GSE65194", "role": "replication", "tissue": "breast",
                "series_matrix": "/home/kai001/breast_cancer_test/datasets/GSE65194_series_matrix.txt.gz",
                "platform": "/home/kai001/breast_cancer_test/platforms/GPL570_clean.annot",
                "phenotype_field": "Sample_characteristics_ch1",
                "case_values": ["sample_group: tnbc", "sample_group: her2", "sample_group: luminal"],
                "control_values": ["sample_group: healthy"],
            },
        ],
    },
}


def upload_file(base_url, endpoint, filepath):
    """Upload a file and return the response JSON."""
    with open(filepath, "rb") as f:
        r = requests.post(f"{base_url}{endpoint}", files={"file": f}, timeout=120)
    r.raise_for_status()
    return r.json()


def run_disease(base_url, name, config):
    """Run a single disease through the UI backend. Returns result dict."""
    result = {"disease": name, "status": "FAIL", "error": "", "study_genes": 0,
              "core_genes": 0, "survived": 0, "meta_sig": 0, "qc": "", "time": 0}
    t0 = time.time()

    try:
        # Upload seeds
        seed_data = upload_file(base_url, "/api/upload/seeds", config["seed_genes"])
        seed_path = seed_data["path"]
        print(f"  Seeds: {seed_data['gene_count']} genes")

        # Upload datasets
        ds_entries = []
        platform_cache = {}
        for ds in config["datasets"]:
            # Upload series matrix
            mat_data = upload_file(base_url, "/api/upload/dataset", ds["series_matrix"])

            # Upload platform (cache to avoid re-uploading same file)
            plat_path = ds["platform"]
            if plat_path not in platform_cache:
                plat_data = upload_file(base_url, "/api/upload/platform", plat_path)
                platform_cache[plat_path] = plat_data["path"]

            ds_entries.append({
                "id": ds["id"],
                "matrix_path": mat_data["path"],
                "platform_path": platform_cache[plat_path],
                "role": ds["role"],
                "tissue": ds["tissue"],
                "phenotype_field": ds["phenotype_field"],
                "case_values": ds["case_values"],
                "control_values": ds["control_values"],
            })
            print(f"  Dataset {ds['id']}: {mat_data.get('sample_count', '?')} samples ({ds['role']})")

        # Build config
        ui_config = {
            "condition": config["condition"],
            "seed_genes_path": seed_path,
            "hgnc_path": config["hgnc_path"],
            "datasets": ds_entries,
        }

        # Validate
        vr = requests.post(f"{base_url}/api/validate", json=ui_config, timeout=30).json()
        if not vr.get("valid"):
            result["error"] = f"Validation failed: {vr.get('error', 'unknown')}"
            return result
        print(f"  Config valid: {vr.get('n_discovery')} discovery + {vr.get('n_replication')} replication")

        # Run
        run_resp = requests.post(f"{base_url}/api/run", json=ui_config,
                                  headers={"Content-Type": "application/json"}, timeout=30).json()
        run_id = run_resp["run_id"]
        print(f"  Run started: {run_id}")

        # Wait for completion via WebSocket
        import websockets

        async def wait_for_completion():
            uri = f"ws://127.0.0.1:8050/ws/progress/{run_id}"
            async with websockets.connect(uri) as ws:
                while True:
                    msg = await asyncio.wait_for(ws.recv(), timeout=1800)
                    data = json.loads(msg)
                    if data["type"] == "phase":
                        phase_name = data.get("name", "").split(":")[-1].strip()[:30]
                        print(f"    Phase {data['phase']}: {phase_name}")
                    elif data["type"] == "complete":
                        return data["summary"]
                    elif data["type"] == "error":
                        raise RuntimeError(data["message"])

        summary = asyncio.run(wait_for_completion())

        result["study_genes"] = summary.get("phase1_study_genes", 0)
        result["core_genes"] = summary.get("phase4_core_genes", 0)
        result["survived"] = summary.get("phase5_survived", 0)
        result["meta_sig"] = summary.get("phase6_significant_random", 0)
        result["qc"] = summary.get("qc_status", "UNKNOWN")
        result["time"] = round(time.time() - t0, 1)

        # Verify results
        res_data = requests.get(f"{base_url}/api/results/{run_id}", timeout=30).json()
        core_list = res_data.get("phase4_core_genes", [])

        checks = []
        checks.append(("core_genes > 0", result["core_genes"] > 0))
        checks.append(("QC passed", result["qc"] == "PASSED"))
        checks.append(("core gene data loaded", len(core_list) > 0))
        checks.append(("study_genes > 0", result["study_genes"] > 0))

        all_pass = all(ok for _, ok in checks)
        if all_pass:
            result["status"] = "PASS"
        else:
            failed = [name for name, ok in checks if not ok]
            result["error"] = f"Failed checks: {', '.join(failed)}"

    except Exception as e:
        result["error"] = str(e)[:200]
        result["time"] = round(time.time() - t0, 1)

    return result


def main():
    parser = argparse.ArgumentParser(description="Test Riker Engine UI with all diseases")
    parser.add_argument("--diseases", default="ASD,T2D,IBD,AD,BRCA",
                        help="Comma-separated disease list (default: all)")
    parser.add_argument("--port", type=int, default=8050, help="Server port")
    args = parser.parse_args()

    selected = [d.strip().upper() for d in args.diseases.split(",")]

    # Start server
    print("Starting Riker Engine UI server...")
    import uvicorn
    from riker.ui.server import app

    server_thread = threading.Thread(
        target=uvicorn.run,
        kwargs={"app": app, "host": "127.0.0.1", "port": args.port, "log_level": "warning"},
        daemon=True,
    )
    server_thread.start()
    time.sleep(2)

    base_url = f"http://127.0.0.1:{args.port}"

    # Verify server is up
    try:
        r = requests.get(base_url, timeout=5)
        assert r.status_code == 200
        print(f"Server running at {base_url}\n")
    except Exception as e:
        print(f"Server failed to start: {e}")
        sys.exit(1)

    # Run each disease
    results = []
    for name in selected:
        if name not in DISEASES:
            print(f"Unknown disease: {name}")
            continue
        print(f"\n{'='*60}")
        print(f"TESTING: {name}")
        print(f"{'='*60}")
        result = run_disease(base_url, name, DISEASES[name])
        results.append(result)
        status_icon = "PASS" if result["status"] == "PASS" else "FAIL"
        print(f"\n  Result: {status_icon}")
        if result["error"]:
            print(f"  Error: {result['error']}")

    # Summary table
    print(f"\n\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"{'Disease':<10} {'Status':<8} {'Study':<8} {'Core':<8} {'Survived':<10} {'Meta-sig':<10} {'Time':<8} {'QC':<8}")
    print("-" * 80)
    for r in results:
        print(f"{r['disease']:<10} {r['status']:<8} {r['study_genes']:<8} {r['core_genes']:<8} "
              f"{r['survived']:<10} {r['meta_sig']:<10} {r['time']:<8} {r['qc']:<8}")

    n_pass = sum(1 for r in results if r["status"] == "PASS")
    n_total = len(results)
    print(f"\n{n_pass}/{n_total} diseases passed")

    if any(r["error"] for r in results):
        print("\nErrors:")
        for r in results:
            if r["error"]:
                print(f"  {r['disease']}: {r['error']}")

    return 0 if n_pass == n_total else 1


if __name__ == "__main__":
    sys.exit(main())
