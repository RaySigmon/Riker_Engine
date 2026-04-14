#!/usr/bin/env python3
"""
Riker Engine — Stability Profiling

Runs the engine N times on the same config to quantify per-gene stochastic
variance from UMAP/HDBSCAN clustering. Produces a stability report classifying
each gene as iron-clad, borderline, or stochastic.

Usage:
    python scripts/stability_profiling.py configs/examples/asd_bulk.yaml -n 10
    python scripts/stability_profiling.py configs/examples/asd_blind.yaml -n 50 --output-dir stability_results/asd_blind

Output:
    <output_dir>/
        stability_scores.csv     — per-gene appearance count and stability class
        run_summary.csv          — per-run core gene count
        stability_report.txt     — human-readable summary
        runs/run_001/            — individual pipeline outputs
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

import yaml


def modify_config_output(config_path: str, new_output_dir: str) -> str:
    """Create a temporary config file with a modified output_dir."""
    with open(config_path) as f:
        config = yaml.safe_load(f)

    config["output_dir"] = new_output_dir

    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".yaml", delete=False, prefix="riker_stability_"
    )
    yaml.dump(config, tmp, default_flow_style=False)
    tmp.close()
    return tmp.name


def run_pipeline(config_path: str, run_dir: str, run_num: int, total: int) -> dict:
    """Run the Riker pipeline once and return results summary."""
    tmp_config = modify_config_output(config_path, run_dir)

    try:
        print(f"\n{'='*60}")
        print(f"RUN {run_num}/{total}")
        print(f"Output: {run_dir}")
        print(f"{'='*60}")

        start = time.time()
        # Find the riker entry point — check common locations
        riker_bin = shutil.which("riker")
        if not riker_bin:
            # Check ~/.local/bin (pip install --user)
            local_bin = os.path.expanduser("~/.local/bin/riker")
            if os.path.exists(local_bin):
                riker_bin = local_bin
        if not riker_bin:
            print("ERROR: 'riker' command not found. Install with: pip install .",
                  file=sys.stderr)
            sys.exit(1)

        result = subprocess.run(
            [riker_bin, "run", tmp_config],
            capture_output=True,
            text=True,
        )
        elapsed = time.time() - start

        # Check for phase4_core_genes.csv
        core_genes_path = os.path.join(run_dir, "phase4_core_genes.csv")
        if os.path.exists(core_genes_path):
            with open(core_genes_path) as f:
                reader = csv.DictReader(f)
                genes = [row["gene"] for row in reader]
            print(f"  Run {run_num}: {len(genes)} core genes in {elapsed:.1f}s")
            return {
                "run": run_num,
                "success": True,
                "core_gene_count": len(genes),
                "genes": genes,
                "elapsed": elapsed,
            }
        else:
            print(f"  Run {run_num}: FAILED — no phase4_core_genes.csv")
            if result.stderr:
                # Print last 5 lines of stderr for debugging
                for line in result.stderr.strip().split("\n")[-5:]:
                    print(f"    {line}")
            return {
                "run": run_num,
                "success": False,
                "core_gene_count": 0,
                "genes": [],
                "elapsed": elapsed,
            }
    finally:
        os.unlink(tmp_config)


def classify_stability(count: int, total_runs: int) -> str:
    """Classify a gene's stability based on appearance frequency."""
    rate = count / total_runs
    if rate >= 0.9:
        return "iron-clad"
    elif rate >= 0.5:
        return "borderline"
    else:
        return "stochastic"


def write_results(output_dir: str, run_results: list, config_path: str):
    """Write stability analysis results."""
    os.makedirs(output_dir, exist_ok=True)

    successful_runs = [r for r in run_results if r["success"]]
    n_runs = len(successful_runs)

    if n_runs == 0:
        print("ERROR: No successful runs. Cannot produce stability report.")
        return

    # Count gene appearances across all successful runs
    gene_counts = Counter()
    for r in successful_runs:
        for gene in r["genes"]:
            gene_counts[gene] += 1

    # Write stability_scores.csv
    scores_path = os.path.join(output_dir, "stability_scores.csv")
    with open(scores_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "gene", "appearance_count", "appearance_rate",
            "stability_class", "total_runs"
        ])
        for gene, count in sorted(gene_counts.items(), key=lambda x: -x[1]):
            rate = count / n_runs
            cls = classify_stability(count, n_runs)
            writer.writerow([gene, count, f"{rate:.4f}", cls, n_runs])

    # Write run_summary.csv
    summary_path = os.path.join(output_dir, "run_summary.csv")
    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["run", "success", "core_gene_count", "elapsed_seconds"])
        for r in run_results:
            writer.writerow([r["run"], r["success"], r["core_gene_count"],
                             f"{r['elapsed']:.1f}"])

    # Classify all genes
    iron_clad = [g for g, c in gene_counts.items()
                 if classify_stability(c, n_runs) == "iron-clad"]
    borderline = [g for g, c in gene_counts.items()
                  if classify_stability(c, n_runs) == "borderline"]
    stochastic = [g for g, c in gene_counts.items()
                  if classify_stability(c, n_runs) == "stochastic"]

    # Compute stats
    counts = [r["core_gene_count"] for r in successful_runs]
    mean_count = sum(counts) / len(counts)
    min_count = min(counts)
    max_count = max(counts)
    total_elapsed = sum(r["elapsed"] for r in run_results)

    # Write human-readable report
    report_path = os.path.join(output_dir, "stability_report.txt")
    with open(report_path, "w") as f:
        f.write("Riker Engine — Stability Profiling Report\n")
        f.write(f"{'='*60}\n\n")
        f.write(f"Config: {config_path}\n")
        f.write(f"Runs attempted: {len(run_results)}\n")
        f.write(f"Runs successful: {n_runs}\n")
        f.write(f"Total time: {total_elapsed:.1f}s "
                f"({total_elapsed/60:.1f}m)\n\n")

        f.write(f"Core gene counts across runs:\n")
        f.write(f"  Min: {min_count}\n")
        f.write(f"  Max: {max_count}\n")
        f.write(f"  Mean: {mean_count:.1f}\n")
        f.write(f"  Variance: ±{((max_count - min_count) / 2):.1f} genes "
                f"({((max_count - min_count) / mean_count * 100 / 2):.1f}%)\n\n")

        f.write(f"Unique genes seen: {len(gene_counts)}\n")
        f.write(f"  Iron-clad (>=90% of runs): {len(iron_clad)}\n")
        f.write(f"  Borderline (50-89%):       {len(borderline)}\n")
        f.write(f"  Stochastic (<50%):         {len(stochastic)}\n\n")

        f.write(f"Iron-clad genes ({len(iron_clad)}):\n")
        for g in sorted(iron_clad):
            c = gene_counts[g]
            f.write(f"  {g}: {c}/{n_runs} ({c/n_runs*100:.0f}%)\n")

        if borderline:
            f.write(f"\nBorderline genes ({len(borderline)}):\n")
            for g in sorted(borderline):
                c = gene_counts[g]
                f.write(f"  {g}: {c}/{n_runs} ({c/n_runs*100:.0f}%)\n")

        if stochastic:
            f.write(f"\nStochastic genes ({len(stochastic)}):\n")
            for g in sorted(stochastic):
                c = gene_counts[g]
                f.write(f"  {g}: {c}/{n_runs} ({c/n_runs*100:.0f}%)\n")

    print(f"\n{'='*60}")
    print("STABILITY PROFILING COMPLETE")
    print(f"{'='*60}")
    print(f"Successful runs: {n_runs}/{len(run_results)}")
    print(f"Core genes: {min_count}–{max_count} (mean {mean_count:.1f})")
    print(f"Unique genes seen: {len(gene_counts)}")
    print(f"  Iron-clad (>=90%): {len(iron_clad)}")
    print(f"  Borderline (50-89%): {len(borderline)}")
    print(f"  Stochastic (<50%): {len(stochastic)}")
    print(f"\nResults: {output_dir}/")
    print(f"  stability_scores.csv")
    print(f"  run_summary.csv")
    print(f"  stability_report.txt")


def main():
    parser = argparse.ArgumentParser(
        description="Run the Riker Engine N times to profile per-gene stability",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "config",
        help="Path to YAML config file",
    )
    parser.add_argument(
        "-n", "--num-runs",
        type=int,
        default=10,
        help="Number of runs (default: 10)",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory for stability results (default: stability_<condition>/)",
    )
    parser.add_argument(
        "--keep-runs",
        action="store_true",
        help="Keep individual run outputs in runs/ subdirectory (default: delete to save space)",
    )
    parser.add_argument(
        "--seed-list",
        default=None,
        help="Optional: path to a seed gene CSV to cross-reference against "
             "(e.g., SFARI genes for ASD). Adds 'in_seed_list' column to stability_scores.csv",
    )

    args = parser.parse_args()

    config_path = os.path.abspath(args.config)
    if not os.path.exists(config_path):
        print(f"Config file not found: {config_path}", file=sys.stderr)
        sys.exit(1)

    # Determine output directory
    with open(config_path) as f:
        config = yaml.safe_load(f)
    condition = config.get("condition", "unknown")

    output_dir = args.output_dir or f"stability_{condition}"
    output_dir = os.path.abspath(output_dir)
    runs_dir = os.path.join(output_dir, "runs")
    os.makedirs(runs_dir, exist_ok=True)

    # Load seed list for cross-reference if provided
    seed_genes = set()
    if args.seed_list:
        with open(args.seed_list) as f:
            reader = csv.DictReader(f)
            # Try common column names
            for row in reader:
                for col in ["symbol", "gene", "gene_symbol", "Symbol", "Gene"]:
                    if col in row:
                        seed_genes.add(row[col])
                        break
        print(f"Loaded {len(seed_genes)} seed genes for cross-reference")

    print(f"Stability Profiling: {condition}")
    print(f"Config: {config_path}")
    print(f"Runs: {args.num_runs}")
    print(f"Output: {output_dir}")

    # Run the pipeline N times
    run_results = []
    for i in range(1, args.num_runs + 1):
        run_dir = os.path.join(runs_dir, f"run_{i:03d}")
        result = run_pipeline(config_path, run_dir, i, args.num_runs)
        run_results.append(result)

    # Write results
    write_results(output_dir, run_results, config_path)

    # Add seed list cross-reference if provided
    if seed_genes:
        scores_path = os.path.join(output_dir, "stability_scores.csv")
        # Read existing scores
        rows = []
        with open(scores_path) as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames + ["in_seed_list"]
            for row in reader:
                row["in_seed_list"] = "yes" if row["gene"] in seed_genes else "no"
                rows.append(row)

        # Rewrite with seed column
        with open(scores_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        # Count novel discoveries
        successful = [r for r in run_results if r["success"]]
        n_runs = len(successful)
        novel_iron = [r for r in rows
                      if r["in_seed_list"] == "no"
                      and r["stability_class"] == "iron-clad"]
        novel_borderline = [r for r in rows
                           if r["in_seed_list"] == "no"
                           and r["stability_class"] == "borderline"]

        print(f"\nSeed list cross-reference:")
        print(f"  Novel iron-clad genes (not in seed list): {len(novel_iron)}")
        for r in novel_iron:
            print(f"    {r['gene']}: {r['appearance_count']}/{n_runs}")
        print(f"  Novel borderline genes: {len(novel_borderline)}")

    # Clean up individual runs if not keeping
    if not args.keep_runs:
        print(f"\nCleaning up individual run outputs...")
        shutil.rmtree(runs_dir, ignore_errors=True)
        print("Done. Use --keep-runs to preserve individual run outputs.")


if __name__ == "__main__":
    main()
