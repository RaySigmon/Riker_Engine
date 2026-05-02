# Negative Control Summary

**Config:** `configs/examples/negative_control.yaml`
**Seed file:** `data/seeds/negative_control_random_500.csv` (500 randomly selected protein-coding genes, SFARI ASD genes excluded)
**Engine version at time of run:** v0.3.2 (pre-SOP, exploratory)
**Datasets:** 3 discovery (brain cortex microarray), 4 replication (1 brain + 3 blood)

## Results

| Phase | Count | Rate (of 500 seeds) |
|-------|-------|---------------------|
| Phase 1 study genes | 38 | 7.6% |
| Phase 4 locked core genes | 5 | 1.0% |
| Phase 5 survivors | 4 | 0.8% |
| Phase 6 random-effects significant | 1 | 0.2% |

**QC status:** PASSED

**Core genes identified:** BBS7, EPB41L3, EXTL2, RETREG2, SNU13
**Phase 5 eliminated:** RETREG2 (1 gene)
**Phase 6 random-effects significant:** SNU13 (p ~ 0.00185)

## Interpretation

The negative control was **not empty**. It produced low-yield residual signal:
5/500 locked core genes (1.0%), 4/500 Phase 5 survivors (0.8%), and 1/500
random-effects meta-significant gene (0.2%).

Claims about the negative control must be framed as **low false-positive yield**,
not zero false positives. A 1.0% core-gene yield from random seeds is consistent
with the expected base rate from chance co-expression across independent datasets.

## Limitations of this negative control

This run used a **different configuration** from the ASD blind disease-day run:
- Missing `phase3:` and `phase4:` blocks (ran with defaults instead of the 15-config
  sweep and 10,000-permutation settings used in disease blind runs)
- Used 5 datasets (missing GSE64018 and GSE102741 RNA-seq datasets present in ASD blind)

Because of these configuration differences, this negative control is classified as
**exploratory evidence only**. A matched negative-control rerun using identical
dataset and parameter settings as the disease blind config is required before
publication-grade false-positive claims can be made.

## Source artifacts

All outputs are in `results/negative_control/`:
- `pipeline_summary.json` — machine-readable summary
- `phase1_study_genes.csv` — 38 genes passing Phase 1
- `phase4_core_genes.csv` — 5 locked core genes
- `phase4_all_levels.csv` — progressive sensitivity levels
- `phase5_verdicts.csv` — replication verdicts
- `phase6_meta_analysis.csv` — meta-analysis statistics
- `qc_report.json` — QC gate results
