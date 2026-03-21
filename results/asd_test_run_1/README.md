# ASD Test Run #1 — Results

**Date:** March 20-21, 2026
**Engine Version:** v0.1.1 (262 tests)
**Condition:** Autism Spectrum Disorder
**Seed Genes:** 1,267 (SFARI Gene database)

## Bulk Tissue Analysis

**Datasets:** 7 GEO expression datasets (3 brain discovery, 1 brain replication, 3 blood/LCL replication)

| Metric | Value |
|---|---|
| Seed genes | 1,267 |
| Study genes (Phase 1) | 142 |
| Consensus clusters (Phase 3) | 13 |
| Core genes (Phase 4) | 35 |
| Survived replication (Phase 5) | 35/35 |
| Meta-significant (Phase 6) | 13 |
| QC status | PASSED (0 warnings) |

**Key genes:** ATP2B2 (log2FC=-0.161, down), SEZ6L2 (log2FC=-0.274, down)

## Cell-Type-Specific snRNA-seq Analysis

**Datasets:** Velmeshev 2019 (104,559 nuclei, pseudo-bulked) + Wamsley 2024 (published DE)

| Cell Type | Study Genes | Clusters | Core Genes |
|---|---|---|---|
| L2/3 Excitatory | 45 | 4 | 0 |
| SST Interneurons | 55 | 2 | **7** |
| Microglia | 22 | 2 | **3** |
| Oligodendrocytes | 24 | 2 | **3** |
| Astrocytes | 43 | 4 | 0 |

**13 cell-type-specific core genes:**
- SST: ABL2, CNTNAP4, SLC27A4, STK39, TAF1, TRAPPC9, UIMC1
- Microglia: CHD4, HEPACAM, SETD1A
- Oligodendrocytes: CD276, GRIK4, MAGEC3

## Files

- `bulk/` — Full pipeline output from 7-dataset bulk analysis
- `snrnaseq_dual/` — Per-cell-type output from dual-dataset snRNA-seq
- `asd_bulk_config.yaml` — Configuration file used for bulk run
