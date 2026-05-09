#!/usr/bin/env python3
"""Run N negative control trials with matched ASD-blind config.

Each trial uses a different random 500-gene seed file, run through the
full Riker pipeline with ASD-blind-matched settings (7 datasets, same
phase3/phase4 config). Produces a summary CSV with per-trial metrics.

Usage:
    python scripts/run_negative_control_trials.py \
        --config configs/negative_control_matched.yaml \
        --seed-dir data/seeds/negative_control_trials \
        --output-dir results/negative_control_matched \
        --n-trials 50
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import yaml


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Template config (matched to ASD blind)")
    parser.add_argument("--seed-dir", required=True, help="Directory with trial_001.csv ... trial_N.csv")
    parser.add_argument("--output-dir", required=True, help="Base output directory")
    parser.add_argument("--n-trials", type=int, default=50)
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    seed_dir = Path(args.seed_dir).resolve()
    output_base = Path(args.output_dir).resolve()
    output_base.mkdir(parents=True, exist_ok=True)

    # Load template config
    with open(config_path) as f:
        template = yaml.safe_load(f)

    # Resolve all relative paths against config directory
    config_dir = config_path.parent
    for field in ["hgnc_path"]:
        if field in template and not os.path.isabs(template[field]):
            template[field] = str((config_dir / template[field]).resolve())
    for ds in template.get("datasets", []):
        for field in ["series_matrix", "platform"]:
            if field in ds and not os.path.isabs(ds[field]):
                ds[field] = str((config_dir / ds[field]).resolve())

    summary_path = output_base / "negative_control_summary.csv"
    log_path = output_base / "runner.log"

    results = []

    with open(log_path, "w") as logf:
        for trial_num in range(1, args.n_trials + 1):
            seed_file = seed_dir / f"trial_{trial_num:03d}.csv"
            if not seed_file.exists():
                logf.write(f"SKIP trial {trial_num}: {seed_file} not found\n")
                logf.flush()
                continue

            trial_dir = output_base / f"trial_{trial_num:03d}"
            trial_dir.mkdir(exist_ok=True)

            # Check if already completed
            summary_file = trial_dir / "pipeline_summary.json"
            if summary_file.exists():
                logf.write(f"SKIP trial {trial_num}: already complete\n")
                logf.flush()
                try:
                    with open(summary_file) as sf:
                        d = json.load(sf)
                    results.append({
                        "trial": trial_num,
                        "seed_file": str(seed_file),
                        "phase1_study_genes": d.get("phase1_study_genes", 0),
                        "phase4_core_genes": d.get("phase4_core_genes", 0),
                        "phase4_significant_clusters": d.get("phase4_significant_clusters", 0),
                        "phase5_survived": d.get("phase5_survived", 0),
                        "phase6_significant": d.get("phase6_significant_random_effects", d.get("phase6_significant_random", 0)),
                        "qc_status": d.get("qc_status", "UNKNOWN"),
                        "wall_seconds": 0,
                        "status": "CACHED",
                    })
                except Exception:
                    pass
                continue

            # Build per-trial config
            trial_config = dict(template)
            trial_config["seed_genes"] = str(seed_file)
            trial_config["output_dir"] = str(trial_dir)

            tmp_config = trial_dir / "config.yaml"
            with open(tmp_config, "w") as cf:
                yaml.dump(trial_config, cf, default_flow_style=False)

            logf.write(f"START trial {trial_num}/{args.n_trials} at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            logf.flush()
            print(f"Trial {trial_num}/{args.n_trials} starting...")

            t0 = time.time()
            try:
                proc = subprocess.run(
                    ["riker", "run", str(tmp_config)],
                    capture_output=True, text=True, timeout=3600
                )
                wall = time.time() - t0
                status = "OK" if proc.returncode == 0 else f"EXIT_{proc.returncode}"

                # Write stdout/stderr to trial dir
                with open(trial_dir / "stdout.log", "w") as f:
                    f.write(proc.stdout)
                with open(trial_dir / "stderr.log", "w") as f:
                    f.write(proc.stderr)

            except subprocess.TimeoutExpired:
                wall = time.time() - t0
                status = "TIMEOUT"
            except Exception as e:
                wall = time.time() - t0
                status = f"ERROR: {e}"

            # Read pipeline_summary if it exists
            row = {
                "trial": trial_num,
                "seed_file": str(seed_file),
                "phase1_study_genes": 0,
                "phase4_core_genes": 0,
                "phase4_significant_clusters": 0,
                "phase5_survived": 0,
                "phase6_significant": 0,
                "qc_status": "NO_OUTPUT",
                "wall_seconds": round(wall, 1),
                "status": status,
            }

            if summary_file.exists():
                try:
                    with open(summary_file) as sf:
                        d = json.load(sf)
                    row["phase1_study_genes"] = d.get("phase1_study_genes", 0)
                    row["phase4_core_genes"] = d.get("phase4_core_genes", 0)
                    row["phase4_significant_clusters"] = d.get("phase4_significant_clusters", 0)
                    row["phase5_survived"] = d.get("phase5_survived", 0)
                    row["phase6_significant"] = d.get("phase6_significant_random_effects",
                                                       d.get("phase6_significant_random", 0))
                    row["qc_status"] = d.get("qc_status", "UNKNOWN")
                except Exception:
                    pass

            results.append(row)
            logf.write(f"DONE  trial {trial_num} in {wall:.0f}s — status={status}, "
                       f"p1={row['phase1_study_genes']}, p4_core={row['phase4_core_genes']}, "
                       f"p6_sig={row['phase6_significant']}\n")
            logf.flush()

            # Write running summary after each trial
            _write_summary(summary_path, results)

    _write_summary(summary_path, results)
    print(f"\nDone. {len(results)} trials. Summary: {summary_path}")

    # Print aggregate stats
    p1_counts = [r["phase1_study_genes"] for r in results if r["status"] in ("OK", "CACHED")]
    p6_counts = [r["phase6_significant"] for r in results if r["status"] in ("OK", "CACHED")]
    if p1_counts:
        print(f"Phase 1 study genes: mean={sum(p1_counts)/len(p1_counts):.1f}, "
              f"min={min(p1_counts)}, max={max(p1_counts)}")
    if p6_counts:
        print(f"Phase 6 significant: mean={sum(p6_counts)/len(p6_counts):.1f}, "
              f"min={min(p6_counts)}, max={max(p6_counts)}")
        nonzero = sum(1 for c in p6_counts if c > 0)
        print(f"Trials with >0 Phase 6 genes: {nonzero}/{len(p6_counts)}")


def _write_summary(path, results):
    if not results:
        return
    fieldnames = list(results[0].keys())
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    main()
