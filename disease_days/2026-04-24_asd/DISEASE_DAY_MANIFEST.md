# ASD Disease-Day Manifest

**Date:** 2026-04-24 (configs prepared) / 2026-04-26 (executed)
**Disease:** Autism Spectrum Disorder (ASD)
**Master seed:** 42 (canonical)
**Executor:** Kai on Ghost (Pi 5, 8GB RAM)
**SOP version:** v1.0.1
**Riker version:** 0.3.2
**Git commit at run time:** 2a7d28e

## Runs completed

### Run 1: Curated (SFARI seeds, 1,267 genes)
- **Config:** disease_days/2026-04-24_asd/curated/config.yaml
- **Wall clock:** ~2 min (08:10:43 → 08:12:39)
- Phase 1 study genes: 141
- Phase 4 core genes: **35**
- Phase 4 significant clusters: 3
- Phase 5 survived: 20 / eliminated: 15
- Phase 6 significant (random effects): 9
- QC: PASSED (4 passed, 0 warnings, 0 critical)
- **Reproduction status:** EXACT match to all historical curated runs

### Run 2: Protein-coding blind (19,296 HGNC protein-coding genes)
- **Config:** disease_days/2026-04-24_asd/blind_pc/config.yaml
- **Wall clock:** ~8 min (08:22:33 → 08:30:21)
- Phase 1 study genes: 1,793
- Phase 4 core genes: **401**
- Phase 4 significant clusters: 21
- Phase 5 survived: 220 / eliminated: 181
- Phase 6 significant (random effects): 135
- QC: PASSED (4 passed, 0 warnings, 0 critical)
- **Reproduction status:** Within expected stochastic range (historical: 401-403)

### Run 3: All-expressed blind (39,185 archived seed identifiers)
- **Config:** disease_days/2026-04-24_asd/allexpressed/config.yaml
- **Wall clock:** ~8 min (08:31:47 → 08:39:33)
- Phase 1 study genes: 1,813
- Phase 4 core genes: **415**
- Phase 4 significant clusters: 21
- Phase 5 survived: 225 / eliminated: 190
- Phase 6 significant (random effects): 140
- QC: PASSED (3 passed, 1 warning, 0 critical)
- **Note:** Seed file preserved at data/seeds/asd_all_expressed.csv (provenance in .provenance companion file). Initial run failed due to # comment headers in CSV; fixed by stripping comments to companion file.
- **Note:** 9 gene symbols remapped via HGNC (e.g., NCL → NUCLEOLIN, ZNF131 → ZBTB35)

## Cross-run comparison

| Metric | Curated | Blind PC | All-expressed |
|--------|---------|----------|---------------|
| Seed genes | 1,267 | 19,296 | 39,185 |
| Phase 1 study genes | 141 (11.1%) | 1,793 (9.3%) | 1,813 (4.6%) |
| Phase 4 core genes | 35 | 401 | 415 |
| Significant clusters | 3 | 21 | 21 |
| Phase 5 survived | 20 (57.1%) | 220 (54.9%) | 225 (54.2%) |
| Phase 6 significant | 9 | 135 | 140 |
| Wall clock | ~2 min | ~8 min | ~8 min |

**Observations:**
- Phase 5 survival rate is consistent across all three modes (~54-57%)
- All-expressed (39K seeds) produces only marginally more study genes than protein-coding (19K) — 1,813 vs 1,793 — because the additional non-coding identifiers mostly don't appear in the expression matrices
- Core gene counts: 415 (all-expressed) vs 401 (protein-coding) — 14 additional genes from the broader seed set

## Cohort inclusion (pre-specified)
- Discovery (3 brain): GSE28521, GSE28475, GSE64018
- Replication (4): GSE102741 (brain) + GSE18123, GSE26415, GSE42133 (blood)
- 332 ASD / 332 control / 664 samples total
- No cohorts excluded; all available ASD tissue-relevant cohorts used

## Flag 6 annotation (from blind_pc run)
- n_significant_clusters (Phase 4 permutation): 21
- n_clusters_contributing_core_genes: 74 (from April 23 audit of run_001 with identical config)
- Cluster-size filter (≥3 Level 2 survivors): excluded 74 Level 2 survivors from 50 clusters with <3 survivors

## Reproducibility verification
- Curated run: EXACT reproduction of historical 35-core-gene result (all phase counts identical)
- Blind PC run: 401 core genes (historical range: 401-403, within stochastic variance from UMAP/HDBSCAN)
- Master seed 42 confirmed deterministic for curated mode; blind mode shows expected ±2-gene stochastic range

## Execution notes
- First SOP-compliant disease-day execution
- SOP v1.0.1 patches applied (SOP_LESSONS.md Entry 001)
- All-expressed CSV required comment header stripping (provenance preserved in .provenance companion file)
- System state at execution: 76GB disk free, 4.0GB RAM available, load avg ~1.8

---

*This manifest follows docs/DISEASE_DAY_SOP.md §"Required provenance metadata".*
