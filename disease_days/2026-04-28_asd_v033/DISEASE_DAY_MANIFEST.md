# Disease Day Manifest — ASD

**Disease:** Autism Spectrum Disorder (ASD)
**Date:** 2026-04-28
**Git commit at run time:** 1b250ed5c7fb436b25292971911fe743a3ff0339
**Riker version:** 0.3.3
**Master seed:** 42
**Executor:** Kai (Claude Code CLI agent) under Cody Sigmon direction

**Note:** ASD ran under v0.3.3 (pre-.2 hotfix). The v0.3.3.2 fixes (symbol_column wiring, tissue requirement) do not affect ASD because (a) ASD's seed file uses `symbol` as column name (matching the silently-used default), and (b) ASD's config had explicit tissue on every dataset. No re-run needed.

---

## Seed file checksums

| Run | Seed file | SHA256 |
|-----|-----------|--------|
| Curated | `asd_sfari_genes.csv` (1,267 genes) | `b74143fcf038970c60f93d61b34d2a6c205164e4b9c08ee821dbb89fecbc8eb3` |
| Blind | `all_protein_coding_genes.csv` (19,296 genes) | `e9ab2790c45a3501545d60298a0fc51c7cd2791b040002387665f671c98639aa` |
| Stability | same as Blind | same as Blind |

---

## GEO datasets used

| Accession | Role | Tissue | Platform |
|-----------|------|--------|----------|
| GSE28521 | discovery | brain | GPL6883 |
| GSE28475 | discovery | brain | GPL6883 |
| GSE64018 | discovery | brain | Ensembl (RNA-seq) |
| GSE102741 | replication | brain | Ensembl (RNA-seq) |
| GSE18123 | replication | blood | GPL6244 |
| GSE26415 | replication | blood | GPL6480 |
| GSE42133 | replication | blood | GPL10558 |

Cohort selection: All available ASD brain tissue cohorts with clear case/control definitions used as discovery. Replication includes brain (GSE102741) and blood (GSE18123, GSE26415, GSE42133) cohorts.

---

## Run results

### Tier 1 — Curated SFARI run
- Wall clock: 2 min 4 sec
- Phase 1 study genes: 141
- Phase 4 core genes: 29
- Phase 5 survived: 29
- Phase 6 meta-significant (random effects): 10

### Tier 2 — Blind protein-coding run
- Wall clock: 7 min 14 sec
- Phase 1 study genes: 1,793
- Phase 4 core genes: 421
- Phase 5 survived: 414
- Phase 6 meta-significant (random effects): 233

### Tier 3 — 50-run stability profile
- Wall clock: 18,011 sec (5.0 hours)
- Runs completed: 50/50 (0 failures)
- Total unique genes seen: 457
- Iron-clad (>=90%): 394
- Borderline (50-89%): 35
- Stochastic (<50%): 28
- Iron-clad fraction: 86.21%
- Pairwise Jaccard median: 0.933 (p25=0.924, p75=0.942)
- Genes at perfect 50/50 reproducibility: 293

---

## Flag 6 annotation
- `n_significant_clusters` (Tier 2 blind): 15
- Cluster-size filter impact: documented in `phase4_all_levels.csv`

---

## Manifest status
This manifest was retrospectively generated on 2026-05-01 from existing artifacts (`pipeline_summary.json`, `config.yaml`, `stability_summary.json`, `DISEASE_REPORT_DRAFT.md`). All numbers verified against committed output files.
