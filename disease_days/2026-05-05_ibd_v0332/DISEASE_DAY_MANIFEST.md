# Disease Day Manifest — Inflammatory Bowel Disease

**Disease:** Inflammatory Bowel Disease (IBD)
**Date:** 2026-05-05
**Git commit at run time:** 21d5c5117d311cff364d503a7b65b31d1515ff86
**Riker version:** 0.3.3.2
**Master seed:** 42
**Executor:** Kai (Claude Code CLI agent) under Cody Sigmon direction
**SOP version:** 1.0.3

---

## Seed file checksums

| Run | Seed file | SHA256 |
|-----|-----------|--------|
| Curated | `ibd_open_targets_genes.csv` (762 genes) | `2956e45d3a620ca41d76ebe8477b165d1f62e79350c21efab03013e027035668` |
| Blind | `all_protein_coding_genes.csv` (19,296 genes) | `e9ab2790c45a3501545d60298a0fc51c7cd2791b040002387665f671c98639aa` |

---

## GEO datasets used

| Accession | Role | Tissue | Platform | Description |
|-----------|------|--------|----------|-------------|
| GSE75214 | discovery | colon | GPL6244 | CD + UC vs control |
| GSE16879 | discovery | colon | GPL570 | CD + UC vs control |
| GSE59071 | discovery | colon | GPL6244 | Crohn + UC vs control |
| GSE87466 | replication | colon | GPL13158 | UC vs normal |
| GSE38713 | replication | colon | GPL570 | UC patient vs non-inflammatory control |
| GSE36807 | replication | colon | GPL570 | Crohn + UC vs healthy control |

3 discovery + 3 replication. All colon tissue. Mixed platforms (GPL6244, GPL570, GPL13158).

---

## Run results

### Tier 1 — Curated IBD run (publishable validation chain)
- Phase 1 study genes: 408
- Phase 4 core genes: 301
- Phase 4 significant clusters: 3
- Phase 5 survived: 288
- Phase 5 eliminated: 11
- Phase 6 meta-significant (random effects): 284
- QC status: PASSED

### v0.3.2 baseline comparison
- v0.3.2: 304 core / 302 survived
- v0.3.3.2: 301 core / 288 survived / 284 meta-sig
- Delta: -3 core genes (within ±2% threshold)
- Survival rate: v0.3.2 99.3% vs v0.3.3.2 95.7% — slight tightening from mean log2FC dilution fix

### Tier 2 — Blind protein-coding run (discovery mode)
- Phase 1 study genes: 8,432 (43.7% yield)
- Phase 4 core genes: 6,027
- Phase 4 significant clusters: 0
- Phase 5 survived: 5,664
- Phase 5 eliminated: 271
- Phase 6 meta-significant (random effects): 5,488
- QC status: PASSED
- **Regime: Global** — 43.7% yield with 0 significant clusters, consistent with CRC (42.9%) and Psoriasis (41.1%)

### Tier 3 — Stability profile (50-run blind)
- Total genes seen across 50 runs: 6,209
- Iron-clad genes (>=90% appearance): 5,729 (92.3%)
- Borderline genes (50-89%): 416
- Stochastic genes (<50%): 64
- Pairwise Jaccard similarity: median 0.9592 (range 0.9479-0.9716)
- Core gene count per run: min 5,968 / max 6,068 / mean 6,018
- Total wall time: 19.7 hours (50 runs × ~24 min each)
- Engine commit at stability run: fab1375ea59be4d8fa7bf1498a6a29abb3432de6

Note: Tier 3 was restarted clean after a power outage interrupted the first attempt at run 6/50. All 50 runs in this profile are from a single continuous execution session.

---

## Cross-disease context

IBD is the 7th of 8 diseases in the v0.3.3.2 validation cycle. Its blind-mode result (43.7% yield, 0 sig clusters) confirms the global regime pattern: diseases with Phase 1 yield >~40% produce 0 significant clusters and a global rewiring signal rather than localized gene modules.

The curated result (301/288/284) is the publishable IBD finding. The 3 significant clusters at curated tier indicate IBD's curated core gene set is dominated by one or two main programs rather than multiple distinct programs.

---

## Files committed

- `curated/pipeline_summary.json` — full curated pipeline output
- `blind_pc/pipeline_summary.json` — full blind pipeline output
- `blind_pc/config.yaml` — blind run configuration
- `stability_50run/stability_summary.json` — machine-readable stability metrics
- `stability_50run/stability_report.csv` — per-gene stability classification
- `stability_50run/stability_pairwise_jaccard.csv` — all 1,225 pairwise Jaccard values
- `stability_50run/run_summary.csv` — per-run core gene counts
- `stability_50run/stability_report.txt` — human-readable summary
- `stability_50run/profiler.log` — profiler execution log
- `DISEASE_DAY_MANIFEST.md` — this file

## Provenance note — Tier 3 engine_commit field

The `stability_summary.json` records `engine_commit: fab1375` (HEAD at summary write time). The engine code that executed all 50 stability runs is identical to commit `a397083` (the blind-tier run commit) — both are v0.3.3.2 with no engine code modifications between them. Only data files (T2D Tier 3 results) were added between these commits. Verify with: `git diff a397083..fab1375 -- riker/` (empty diff).

This cosmetic mismatch arises because the stability profiler records HEAD at write time, not at run start. Tracked as a v0.4.0 improvement.

---

Per-run outputs (`runs/run_001/` through `runs/run_050/`) are not committed. They are regenerable from the committed config and master seed using:
```bash
python3 scripts/stability_profiling.py disease_days/2026-05-05_ibd_v0332/blind_pc/config.yaml -n 50 --output-dir disease_days/2026-05-05_ibd_v0332/stability_50run
```
