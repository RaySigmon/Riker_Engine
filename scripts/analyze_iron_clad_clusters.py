#!/usr/bin/env python3
"""
Analyze iron-clad gene cluster organization from a 50-run stability profile.

Reads per-run phase3_cluster_assignments.csv and phase4_core_genes.csv to
determine each iron-clad gene's modal cluster (most frequently assigned
cluster across runs).

Usage:
    python scripts/analyze_iron_clad_clusters.py \
        disease_days/2026-04-28_asd_v033/stability_50run

Output:
    iron_clad_cluster_analysis.csv  — per-gene modal cluster and assignment count
    iron_clad_cluster_summary.csv   — per-cluster member list and stability
"""

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_iron_clad_clusters.py <stability_output_dir>")
        sys.exit(1)

    stability_dir = Path(sys.argv[1])
    runs_dir = stability_dir / "runs"
    scores_path = stability_dir / "stability_scores.csv"

    if not scores_path.exists():
        print(f"stability_scores.csv not found at {scores_path}")
        sys.exit(1)

    # Load iron-clad genes
    iron_clad = set()
    with open(scores_path) as f:
        for row in csv.DictReader(f):
            if row["stability_class"] == "iron-clad":
                iron_clad.add(row["gene"])
    print(f"Iron-clad genes: {len(iron_clad)}")

    # Track cluster assignments per iron-clad gene across runs
    gene_cluster_counts = defaultdict(Counter)
    n_runs = 0

    for run_dir in sorted(runs_dir.iterdir()):
        if not run_dir.is_dir() or not run_dir.name.startswith("run_"):
            continue
        p3 = run_dir / "phase3_cluster_assignments.csv"
        p4 = run_dir / "phase4_core_genes.csv"
        if not p3.exists() or not p4.exists():
            continue

        n_runs += 1
        run_core = set()
        with open(p4) as f:
            for row in csv.DictReader(f):
                run_core.add(row["gene"])

        with open(p3) as f:
            for row in csv.DictReader(f):
                gene = row["gene"]
                if gene in iron_clad and gene in run_core:
                    gene_cluster_counts[gene][int(row["cluster_id"])] += 1

    print(f"Runs analyzed: {n_runs}")

    # Per-gene analysis
    out_gene = stability_dir / "iron_clad_cluster_analysis.csv"
    with open(out_gene, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene", "modal_cluster", "modal_count", "total_runs",
                         "modal_fraction", "n_distinct_clusters"])
        for gene in sorted(iron_clad):
            if gene in gene_cluster_counts:
                modal_cid, modal_cnt = gene_cluster_counts[gene].most_common(1)[0]
                n_distinct = len(gene_cluster_counts[gene])
            else:
                modal_cid, modal_cnt, n_distinct = -1, 0, 0
            writer.writerow([gene, modal_cid, modal_cnt, n_runs,
                             f"{modal_cnt/n_runs:.3f}" if n_runs > 0 else "0",
                             n_distinct])
    print(f"Wrote: {out_gene}")

    # Per-cluster summary
    cluster_members = defaultdict(list)
    cluster_modal_counts = defaultdict(list)
    for gene in sorted(iron_clad):
        if gene in gene_cluster_counts:
            modal_cid, modal_cnt = gene_cluster_counts[gene].most_common(1)[0]
            cluster_members[modal_cid].append(gene)
            cluster_modal_counts[modal_cid].append(modal_cnt)

    out_cluster = stability_dir / "iron_clad_cluster_summary.csv"
    with open(out_cluster, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cluster_id", "n_iron_clad_members", "members",
                         "min_modal_count", "max_modal_count", "mean_modal_count"])
        for cid in sorted(cluster_members.keys()):
            members = cluster_members[cid]
            counts = cluster_modal_counts[cid]
            writer.writerow([
                cid, len(members), "|".join(members),
                min(counts), max(counts),
                f"{sum(counts)/len(counts):.1f}",
            ])
    print(f"Wrote: {out_cluster}")

    # Print summary
    print(f"\nClusters with >= 5 iron-clad members:")
    for cid in sorted(cluster_members.keys()):
        if len(cluster_members[cid]) >= 5 and cid != -1:
            print(f"  Cluster {cid}: {len(cluster_members[cid])} genes")


if __name__ == "__main__":
    main()
