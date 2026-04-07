# IPF Cold Replication Test

## Design

Tested whether the 190 core genes identified by the Riker Engine
(from GSE32537, GSE53845, GSE24206, GSE110147, GSE10667) are independently
differentially expressed in **GSE47460-GPL6480** — a dataset the engine never saw.

**Method:** Simple Welch's t-test per gene. No pipeline, no clustering.
Just: "is this gene differentially expressed in the same direction in an
independent cohort?"

**Held-out dataset:** GSE47460-GPL6480 (LGRC cohort)
- 38 IPF/UIP cases, 17 controls
- Agilent 4x44K platform (GPL6480)
- Filtered to UIP/IPF subtype only (excluded COPD and non-IPF ILD)

## Results

| Metric | Count | Rate |
|--------|-------|------|
| Core genes tested | 190 | — |
| Found in GSE47460 | 153 | 80.5% |
| Significant (p<0.05) | 134 | 87.6% |
| Concordant direction | 148 | 96.7% |
| Concordant + significant | 132 | 86.3% |

## Interpretation

Strong replication. 132 of 153 core genes (86.3%) show
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
