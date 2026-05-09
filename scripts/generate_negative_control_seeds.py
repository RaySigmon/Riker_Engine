#!/usr/bin/env python3
"""Generate N independent random 500-gene seed files for negative control trials.

Each trial draws 500 genes from the protein-coding genome, excluding any genes
in the ASD curated seed list. Seeds are deterministic from (master_seed + trial_num).

Usage:
    python scripts/generate_negative_control_seeds.py --n-trials 50 --output-dir data/seeds/negative_control_trials
"""

import argparse
import csv
import os
import random


def load_gene_list(path):
    with open(path) as f:
        reader = csv.DictReader(f)
        col = "symbol" if "symbol" in reader.fieldnames else reader.fieldnames[0]
        return [row[col].strip() for row in reader if row[col].strip()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-trials", type=int, default=50)
    parser.add_argument("--genes-per-trial", type=int, default=500)
    parser.add_argument("--protein-coding", default="data/seeds/all_protein_coding_genes.csv")
    parser.add_argument("--exclude", default="data/seeds/asd_sfari_genes.csv",
                        help="Gene list to exclude (ASD curated seeds)")
    parser.add_argument("--output-dir", default="data/seeds/negative_control_trials")
    parser.add_argument("--master-seed", type=int, default=42)
    args = parser.parse_args()

    all_genes = load_gene_list(args.protein_coding)
    exclude = set()
    if os.path.exists(args.exclude):
        exclude = set(load_gene_list(args.exclude))
        print(f"Excluding {len(exclude)} genes from {args.exclude}")

    eligible = [g for g in all_genes if g not in exclude]
    print(f"Eligible pool: {len(eligible)} genes (from {len(all_genes)} total, {len(exclude)} excluded)")

    os.makedirs(args.output_dir, exist_ok=True)

    for trial in range(1, args.n_trials + 1):
        rng = random.Random(args.master_seed + trial)
        sample = rng.sample(eligible, args.genes_per_trial)
        sample.sort()

        path = os.path.join(args.output_dir, f"trial_{trial:03d}.csv")
        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["symbol"])
            for gene in sample:
                writer.writerow([gene])

    print(f"Wrote {args.n_trials} seed files to {args.output_dir}")


if __name__ == "__main__":
    main()
