# Disease Day Manifest — IPF

**Disease:** Idiopathic Pulmonary Fibrosis (IPF)
**Date:** 2026-04-29
**Git commit at run time:** 939663e0f1f39e01efd53f61253b73f084ee74fd
**Riker version:** 0.3.3.2
**Master seed:** 42
**Executor:** Kai (Claude Code CLI agent) under Cody Sigmon direction

---

## Seed file checksums

| Run | Seed file | SHA256 |
|-----|-----------|--------|
| Curated | `ipf_curated_genes.csv` (352 genes, 8 HGNC-remapped) | `c28c863b71363c5b43bce997c5bc7d302b3a2113890cad6dc5de75bb50fab1d7` |
| Blind | `all_protein_coding_genes.csv` (19,296 genes) | `e9ab2790c45a3501545d60298a0fc51c7cd2791b040002387665f671c98639aa` |
| Stability | same as Blind | same as Blind |

---

## GEO datasets used

| Accession | Role | Tissue | Platform |
|-----------|------|--------|----------|
| GSE32537 | discovery | lung | GPL6244 |
| GSE53845 | discovery | lung | GPL6480 |
| GSE24206 | discovery | lung | GPL570 |
| GSE110147 | replication | lung | GPL6244 |
| GSE10667 | replication | lung | GPL4133 |

Cohort selection: Three largest available IPF lung tissue cohorts with clear case/control definitions used as discovery. Two additional lung tissue cohorts used as replication. All tissue-matched (lung-to-lung); no cross-tissue tolerance required.

---

## Run results

### Tier 1 — Curated IPF run
- Wall clock: 1 min 36 sec
- Phase 1 study genes: 241
- Phase 4 core genes: 186
- Phase 5 survived: 166 (20 eliminated)
- Phase 6 meta-significant (random effects): 155

### Tier 2 — Blind protein-coding run
- Wall clock: 13 min 31 sec
- Phase 1 study genes: 4,449
- Phase 4 core genes: 2,451
- Phase 5 survived: 1,856 (595 eliminated)
- Phase 6 meta-significant (random effects): 1,584

### Tier 3 — 50-run stability profile
- Wall clock: 33,272 sec (9.2 hours)
- Runs completed: 50/50 (0 failures)
- Total unique genes seen: 2,545
- Iron-clad (>=90%): 2,309
- Borderline (50-89%): 180
- Stochastic (<50%): 56
- Iron-clad fraction: 90.73%
- Pairwise Jaccard median: 0.9547 (p25=0.9512, p75=0.9576)
- Genes at perfect 50/50 reproducibility: 1,946
- Phase-4-stable but Phase-5-unstable: 548 (23.7% of iron-clad)

---

## Flag 6 annotation
- `n_significant_clusters` (Tier 2 blind): documented in `pipeline_summary.json`
- 204 clusters with >=5 iron-clad members
- Cluster-size filter impact: documented in `phase4_all_levels.csv`

---

## Inaugural stability profile
This is the first-ever stability profile for IPF. No v0.3.2 stability baseline exists for comparison. v0.3.2 only had a curated single-run (190 core genes). Section 9 of the disease report uses inaugural-stability-profile framing.

## v0.3.2 comparison (curated only)
v0.3.2 curated: 190 core genes. v0.3.3.2 curated: 186 core genes. 186 stable across versions (97.9%). 4 dropped (FGF7, GREM1, LUM, PLOD2), 0 newly added.

---

## Provenance status
This disease-day predates SOP v1.0.3. `exact_commands.txt`, `environment.txt`,
and `pip_freeze.txt` were not captured at execution time. Available provenance
includes `config.yaml`, `run.log`, `pipeline_summary.json`, `qc_report.json`,
and all phase output files.

This manifest was retrospectively generated on 2026-05-01 from existing artifacts
(`pipeline_summary.json`, `config.yaml`, `stability_summary.json`,
`DISEASE_REPORT_DRAFT.md`). All numbers verified against committed output files.
