#!/usr/bin/env python3
"""
Cold Replication Test — IPF Core Genes in Held-Out Dataset

Tests whether the 170 core genes identified by the Riker Engine in 5 IPF
datasets are independently differentially expressed in GSE47460-GPL6480,
a dataset the engine NEVER saw.

This is a simple Welch's t-test — no pipeline, no clustering, no fancy
statistics. Just: "are these genes actually differentially expressed in
an independent cohort?"

Usage:
    python scripts/cold_replication_ipf.py

Output:
    results/ipf/cold_replication/replication_results.csv
    results/ipf/cold_replication/README.md
"""

import csv
import gzip
import io
import os
import sys

import numpy as np
from scipy import stats

# Paths
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORE_GENES_PATH = os.path.join(REPO_ROOT, "results/ipf/curated/phase4_core_genes.csv")
META_PATH = os.path.join(REPO_ROOT, "results/ipf/curated/phase6_meta_analysis.csv")
SERIES_MATRIX = os.path.join(REPO_ROOT, "data/geo/ipf/GSE47460-GPL6480_series_matrix.txt.gz")
PLATFORM_PATH = os.path.join(REPO_ROOT, "data/platforms/GPL6480.annot")
OUTPUT_DIR = os.path.join(REPO_ROOT, "results/ipf/cold_replication")


def parse_series_matrix(path):
    """Parse GEO series matrix into expression dict and sample metadata."""
    with gzip.open(path, "rt", errors="replace") as f:
        lines = f.readlines()

    # Extract phenotype data
    phenotype_lines = [l for l in lines if l.startswith("!Sample_")]
    sample_ids = None
    characteristics = {}

    for line in phenotype_lines:
        parts = line.strip().split("\t")
        key = parts[0]
        values = [v.strip('"') for v in parts[1:]]

        if key == "!Sample_geo_accession":
            sample_ids = values
        elif key == "!Sample_characteristics_ch1":
            if sample_ids:
                for sid, val in zip(sample_ids, values):
                    if sid not in characteristics:
                        characteristics[sid] = []
                    characteristics[sid].append(val.lower())

    # Parse expression data
    data_start = None
    data_end = None
    for i, line in enumerate(lines):
        if line.startswith("!series_matrix_table_begin"):
            data_start = i + 1
        elif line.startswith("!series_matrix_table_end"):
            data_end = i

    header = lines[data_start].strip().split("\t")
    header = [h.strip('"') for h in header]
    probe_col = header[0]
    sample_cols = header[1:]

    expression = {}
    for line in lines[data_start + 1:data_end]:
        parts = line.strip().split("\t")
        probe = parts[0].strip('"')
        values = []
        for v in parts[1:]:
            try:
                values.append(float(v))
            except (ValueError, IndexError):
                values.append(np.nan)
        expression[probe] = dict(zip(sample_cols, values))

    return expression, characteristics, sample_cols


def parse_platform(path):
    """Parse platform annotation, skipping GEO header lines."""
    probe_to_gene = {}
    with open(path, "r", errors="replace") as f:
        lines = [l for l in f if not l.startswith(("^", "!", "#"))]

    reader = csv.DictReader(io.StringIO("".join(lines)), delimiter="\t")
    gene_col = None
    for col in reader.fieldnames:
        if col.lower().replace(" ", ".") in ("gene.symbol", "gene_symbol", "symbol"):
            gene_col = col
            break

    if not gene_col:
        print(f"Available columns: {reader.fieldnames[:10]}")
        raise ValueError("Could not find gene symbol column")

    # Re-read with the found column
    reader = csv.DictReader(io.StringIO("".join(lines)), delimiter="\t")
    for row in reader:
        gene = row.get(gene_col, "").strip()
        probe = row.get("ID", "").strip()
        if gene and probe and "///" not in gene:
            probe_to_gene[probe] = gene

    return probe_to_gene


def assign_phenotypes(characteristics, sample_ids):
    """Assign case/control labels based on disease state."""
    cases = []
    controls = []

    for sid in sample_ids:
        chars = characteristics.get(sid, [])
        combined = " | ".join(chars)

        # GSE47460: ILD with UIP/IPF subtype = case, control = control
        if "interstitial lung disease" in combined and "2-uip/ipf" in combined:
            cases.append(sid)
        elif "control" in combined and "interstitial lung disease" not in combined:
            controls.append(sid)
        # Skip COPD and non-IPF ILD samples

    return cases, controls


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load core genes and their directions from the engine run
    print("Loading Riker Engine core genes...")
    core_genes = {}
    with open(CORE_GENES_PATH) as f:
        for row in csv.DictReader(f):
            core_genes[row["gene"]] = row["direction"]
    print(f"  {len(core_genes)} core genes")

    # Load meta-analysis for effect sizes
    meta_effects = {}
    with open(META_PATH) as f:
        for row in csv.DictReader(f):
            meta_effects[row["gene"]] = {
                "random_effect": float(row["random_effect"]),
                "random_p": float(row["random_p"]),
                "direction": row["direction"],
            }

    # Parse held-out dataset
    print(f"\nParsing held-out dataset: GSE47460-GPL6480...")
    expression, characteristics, sample_ids = parse_series_matrix(SERIES_MATRIX)
    print(f"  {len(expression)} probes × {len(sample_ids)} samples")

    # Map probes to genes
    print(f"Mapping probes to genes...")
    probe_to_gene = parse_platform(PLATFORM_PATH)

    # Collapse to gene-level expression (mean across probes)
    gene_expression = {}
    for probe, values in expression.items():
        gene = probe_to_gene.get(probe)
        if gene:
            if gene not in gene_expression:
                gene_expression[gene] = {}
            for sid, val in values.items():
                if sid not in gene_expression[gene]:
                    gene_expression[gene][sid] = []
                gene_expression[gene][sid].append(val)

    # Average multiple probes per gene
    for gene in gene_expression:
        for sid in gene_expression[gene]:
            vals = [v for v in gene_expression[gene][sid] if np.isfinite(v)]
            gene_expression[gene][sid] = np.mean(vals) if vals else np.nan

    print(f"  {len(gene_expression)} unique genes")

    # Assign phenotypes
    cases, controls = assign_phenotypes(characteristics, sample_ids)
    print(f"  IPF/UIP cases: {len(cases)}")
    print(f"  Controls: {len(controls)}")

    if len(cases) < 3 or len(controls) < 3:
        print("ERROR: Not enough cases or controls. Check phenotype parsing.")
        sys.exit(1)

    # Run Welch's t-test for each core gene
    print(f"\nRunning Welch's t-test on {len(core_genes)} core genes...")
    results = []
    n_found = 0
    n_significant = 0
    n_concordant = 0
    n_concordant_sig = 0

    for gene, engine_direction in core_genes.items():
        if gene not in gene_expression:
            results.append({
                "gene": gene,
                "engine_direction": engine_direction,
                "in_dataset": False,
                "log2fc": "",
                "p_value": "",
                "replication_direction": "",
                "concordant": "",
                "significant": "",
            })
            continue

        n_found += 1
        case_vals = [gene_expression[gene].get(s, np.nan) for s in cases]
        ctrl_vals = [gene_expression[gene].get(s, np.nan) for s in controls]

        case_vals = [v for v in case_vals if np.isfinite(v)]
        ctrl_vals = [v for v in ctrl_vals if np.isfinite(v)]

        if len(case_vals) < 2 or len(ctrl_vals) < 2:
            results.append({
                "gene": gene,
                "engine_direction": engine_direction,
                "in_dataset": True,
                "log2fc": "",
                "p_value": "",
                "replication_direction": "",
                "concordant": "",
                "significant": "",
            })
            continue

        t_stat, p_value = stats.ttest_ind(case_vals, ctrl_vals, equal_var=False)
        log2fc = np.mean(case_vals) - np.mean(ctrl_vals)
        rep_direction = "up" if log2fc > 0 else "down"
        concordant = (rep_direction == engine_direction)
        significant = (p_value < 0.05)

        if concordant:
            n_concordant += 1
        if significant:
            n_significant += 1
        if concordant and significant:
            n_concordant_sig += 1

        results.append({
            "gene": gene,
            "engine_direction": engine_direction,
            "in_dataset": True,
            "log2fc": f"{log2fc:.4f}",
            "p_value": f"{p_value:.6f}",
            "replication_direction": rep_direction,
            "concordant": concordant,
            "significant": significant,
        })

    # Sort by p-value
    results.sort(key=lambda r: float(r["p_value"]) if r["p_value"] else 999)

    # Write results
    csv_path = os.path.join(OUTPUT_DIR, "replication_results.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "gene", "engine_direction", "in_dataset", "log2fc",
            "p_value", "replication_direction", "concordant", "significant",
        ])
        writer.writeheader()
        writer.writerows(results)

    # Print summary
    print(f"\n{'='*60}")
    print(f"COLD REPLICATION RESULTS — IPF")
    print(f"{'='*60}")
    print(f"Engine core genes:           {len(core_genes)}")
    print(f"Found in GSE47460:           {n_found}")
    print(f"Significant (p<0.05):        {n_significant}")
    print(f"Same direction as engine:    {n_concordant}")
    print(f"Concordant AND significant:  {n_concordant_sig}")
    print(f"")
    print(f"Concordance rate:            {n_concordant}/{n_found} ({100*n_concordant/n_found:.1f}%)")
    print(f"Significant replication:     {n_concordant_sig}/{n_found} ({100*n_concordant_sig/n_found:.1f}%)")
    print(f"{'='*60}")

    # Write README
    readme_path = os.path.join(OUTPUT_DIR, "README.md")
    with open(readme_path, "w") as f:
        f.write(f"""# IPF Cold Replication Test

## Design

Tested whether the {len(core_genes)} core genes identified by the Riker Engine
(from GSE32537, GSE53845, GSE24206, GSE110147, GSE10667) are independently
differentially expressed in **GSE47460-GPL6480** — a dataset the engine never saw.

**Method:** Simple Welch's t-test per gene. No pipeline, no clustering.
Just: "is this gene differentially expressed in the same direction in an
independent cohort?"

**Held-out dataset:** GSE47460-GPL6480 (LGRC cohort)
- {len(cases)} IPF/UIP cases, {len(controls)} controls
- Agilent 4x44K platform (GPL6480)
- Filtered to UIP/IPF subtype only (excluded COPD and non-IPF ILD)

## Results

| Metric | Count | Rate |
|--------|-------|------|
| Core genes tested | {len(core_genes)} | — |
| Found in GSE47460 | {n_found} | {100*n_found/len(core_genes):.1f}% |
| Significant (p<0.05) | {n_significant} | {100*n_significant/n_found:.1f}% |
| Concordant direction | {n_concordant} | {100*n_concordant/n_found:.1f}% |
| Concordant + significant | {n_concordant_sig} | {100*n_concordant_sig/n_found:.1f}% |

## Interpretation

{"Strong replication." if n_concordant_sig/n_found > 0.6 else "Moderate replication." if n_concordant_sig/n_found > 0.4 else "Weak replication — investigate."} {n_concordant_sig} of {n_found} core genes ({100*n_concordant_sig/n_found:.1f}%) show
significant differential expression in the same direction in a completely
independent dataset. This was achieved with no pipeline involvement — just
a basic t-test on the held-out data.

## Reproducibility

```bash
python scripts/cold_replication_ipf.py
```

## Files

- `replication_results.csv` — per-gene results
- `README.md` — this summary
""")

    print(f"\nResults: {csv_path}")
    print(f"README:  {readme_path}")


if __name__ == "__main__":
    main()
