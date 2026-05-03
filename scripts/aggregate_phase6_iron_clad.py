#!/usr/bin/env python3
"""
Aggregate Phase 6 meta-analysis results across stability runs for iron-clad genes.

The 50 stability runs test the same meta-analysis hypothesis on the same
underlying data with stochastic variation only in clustering and permutation
seeds. We are not combining independent tests; we are estimating the variance
of a single meta-analysis result under stochastic re-clustering.

In this framing:
- Median p-value characterizes the typical significance estimate
- Std dev of effect size characterizes its variability across runs
- Direction concordance flags genes where the effect direction itself is
  unstable, not just the magnitude

Median is the default reporting metric because it is robust to outlier runs.
If a reviewer pushes back, harmonic mean p-value (HMP) is a recognized
alternative; re-aggregation with HMP can be done from the per-run data.

For direction: "modal direction" means whichever of 'up' or 'down' is more
common across runs. direction_concordance reports the fraction of runs
agreeing. Genes with low concordance (< 0.80) deserve flagging because
their effect direction is itself unstable.

Usage:
    python scripts/aggregate_phase6_iron_clad.py \
        disease_days/2026-04-28_asd_v033/stability_50run

Output:
    <stability_dir>/iron_clad_aggregated.csv
"""

import csv
import sys
from collections import Counter
from pathlib import Path

import numpy as np


def main():
    if len(sys.argv) < 2:
        print("Usage: python aggregate_phase6_iron_clad.py <stability_output_dir>")
        sys.exit(1)

    stability_dir = Path(sys.argv[1])
    runs_dir = stability_dir / "runs"
    scores_path = stability_dir / "stability_report.csv"

    if not scores_path.exists():
        print(f"stability_report.csv not found at {scores_path}")
        sys.exit(1)

    # Load iron-clad genes
    iron_clad = set()
    with open(scores_path) as f:
        for row in csv.DictReader(f):
            if row["stability_class"] == "iron-clad":
                iron_clad.add(row["gene"])
    print(f"Iron-clad genes: {len(iron_clad)}")

    # Read phase6_meta_analysis.csv from each run
    gene_effects = {}  # gene -> list of random_effect floats
    gene_pvalues = {}  # gene -> list of random_p floats
    gene_directions = {}  # gene -> list of direction strings

    for gene in iron_clad:
        gene_effects[gene] = []
        gene_pvalues[gene] = []
        gene_directions[gene] = []

    n_runs_read = 0
    for run_dir in sorted(runs_dir.iterdir()):
        if not run_dir.is_dir() or not run_dir.name.startswith("run_"):
            continue
        p6_path = run_dir / "phase6_meta_analysis.csv"
        if not p6_path.exists():
            continue
        n_runs_read += 1
        with open(p6_path) as f:
            for row in csv.DictReader(f):
                gene = row["gene"]
                if gene in iron_clad:
                    gene_effects[gene].append(float(row["random_effect"]))
                    gene_pvalues[gene].append(float(row["random_p"]))
                    gene_directions[gene].append(row["direction"])

    print(f"Runs read: {n_runs_read}")

    # Compute aggregated statistics per gene
    out_path = stability_dir / "iron_clad_aggregated.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "gene", "n_runs_in_phase6",
            "median_random_effect", "std_random_effect", "iqr_random_effect",
            "median_random_p", "direction", "direction_concordance",
        ])

        for gene in sorted(iron_clad):
            effects = gene_effects[gene]
            pvalues = gene_pvalues[gene]
            directions = gene_directions[gene]
            n_runs = len(effects)

            if n_runs == 0:
                writer.writerow([gene, 0, 0, 0, 0, 1.0, "none", 0.0])
                continue

            effects_arr = np.array(effects)
            pvalues_arr = np.array(pvalues)

            median_effect = float(np.median(effects_arr))
            std_effect = float(np.std(effects_arr))
            p25 = float(np.percentile(effects_arr, 25))
            p75 = float(np.percentile(effects_arr, 75))
            iqr_effect = p75 - p25
            median_p = float(np.median(pvalues_arr))

            # Modal direction and concordance
            dir_counts = Counter(directions)
            modal_dir = dir_counts.most_common(1)[0][0]
            concordance = dir_counts[modal_dir] / n_runs

            writer.writerow([
                gene, n_runs,
                f"{median_effect:.6f}", f"{std_effect:.6f}", f"{iqr_effect:.6f}",
                f"{median_p:.6e}", modal_dir, f"{concordance:.4f}",
            ])

    print(f"Wrote: {out_path}")
    print(f"Genes with 0 phase6 appearances: {sum(1 for g in iron_clad if len(gene_effects[g]) == 0)}")

    # Quick summary
    all_concordances = []
    for gene in iron_clad:
        dirs = gene_directions[gene]
        if dirs:
            dc = Counter(dirs)
            all_concordances.append(dc.most_common(1)[0][1] / len(dirs))
    if all_concordances:
        low_concordance = sum(1 for c in all_concordances if c < 0.80)
        print(f"Genes with direction concordance < 0.80: {low_concordance}")


if __name__ == "__main__":
    main()
