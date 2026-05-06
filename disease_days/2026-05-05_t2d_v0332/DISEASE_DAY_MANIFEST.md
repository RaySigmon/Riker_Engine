# Disease Day Manifest — Type 2 Diabetes

**Disease:** Type 2 Diabetes (T2D)
**Date:** 2026-05-05
**Git commit at run time:** f3234ee7dac68ef76a05b5bc67cf4e89fb7194e8
**Riker version:** 0.3.3.2
**Master seed:** 42
**Executor:** Kai (Claude Code CLI agent) under Cody Sigmon direction
**SOP version:** 1.0.3

---

## Seed file checksums

| Run | Seed file | SHA256 |
|-----|-----------|--------|
| Curated | `t2d_open_targets_genes.csv` (443 genes) | `5bfe54717a2079cb97d6b614c947b39a9ea54f5f88d622f047ced2d8de565796` |
| Blind | `all_protein_coding_genes.csv` (19,296 genes) | `e9ab2790c45a3501545d60298a0fc51c7cd2791b040002387665f671c98639aa` |

---

## GEO datasets used

| Accession | Role | Tissue | Platform | Description |
|-----------|------|--------|----------|-------------|
| GSE41762 | discovery | islet | GPL6244 | Pancreatic islets, diabetic vs non-diabetic |
| GSE25724 | discovery | islet | GPL96 | Pancreatic islets, T2D vs non-diabetic |
| GSE20966 | discovery | islet | GPL1352 | Beta cells (LCM), T2D vs non-diabetic |
| GSE86468 | replication | islet | Ensembl (RNA-seq) | Islet RNA-seq, T2D vs non-diabetic |

3 discovery + 1 replication. All pancreatic islet tissue. Mixed platforms
(2 microarray + 1 LCM microarray + 1 RNA-seq replication).

---

## Run results

### Tier 1 — Curated T2D run (publishable validation chain)
- Wall clock: ~1 min
- Phase 1 study genes: 56
- Phase 4 core genes: 8
- Phase 4 significant clusters: 1
- Phase 5 survived: 8 (0 eliminated)
- Phase 6 meta-significant (random effects): 8

**Validation:** Exact match to v0.3.2 baseline (8/8/8). Core genes: ABCC8,
CACNA1D, CASR, LRFN2, OTULINL, RASGRP1, RNF6, SLC2A2. All 8 are 100%
meta-significant. The tightest, most precise result in the validation.

### Tier 2 — Blind protein-coding run (publishable — localized regime)
- Wall clock: ~6 min
- Phase 1 study genes: 1,545 (8.0% yield — lowest observed)
- Phase 4 core genes: 169
- **Phase 4 significant clusters: 16** (localized regime confirmed)
- Phase 5 survived: 162 (0 eliminated, 7 insufficient data)
- Phase 6 meta-significant (random effects): 78

**Key finding: IAPP (islet amyloid polypeptide) present in blind core genes.**
IAPP is not in the Open Targets curated seed list. The engine independently
identified the pathological hallmark of T2D from raw expression data across
4 pancreatic islet datasets without being told to look for it. IAPP:
- Forms islet amyloid deposits (T2D pathological hallmark)
- Co-secreted with insulin by beta cells
- Current drug target (pramlintide, FDA-approved IAPP analog)

This is the engine's strongest blind-mode novel rediscovery: a canonical
disease gene found without prior hypothesis.

### Tier 3 — 50-run stability profile
- Wall clock: 14,069 sec (3.9 hours)
- Runs completed: 50/50 (0 failures)
- Total unique genes seen: 205
- Iron-clad (>=90%): 138
- Borderline (50-89%): 30
- Stochastic (<50%): 37
- Iron-clad fraction: 67.3%
- Pairwise Jaccard median: 0.8514
- Core gene counts: min=156, max=179, mean=166.0

**IAPP stability: 49/50 runs (98%) — iron-clad.** The engine's headline
T2D blind rediscovery is rock-solid reproducible.

**Note on lower iron-clad fraction:** T2D's 67.3% iron-clad (vs 86-94% for
other diseases) reflects the weak underlying signal (8% Phase 1 yield, only
1,545 study genes). With fewer genes passing Phase 1, the consensus clustering
has less material to work with and more stochastic variance per run. This is
the expected behavior for the engine's weakest-signal disease — the method
is honest about where reproducibility degrades.

---

## Regime classification

T2D is **strongly localized** — 8.0% Phase 1 yield with 16 significant clusters.
The lowest yield and smallest core gene set in the validation, producing the
cleanest discriminating signal.

| Disease | Phase 1 yield | Sig clusters | Regime |
|---------|--------------|-------------|--------|
| **T2D** | **8.0%** | **16** | **Localized** |
| ASD | 9.3% | 15 | Localized |
| IPF | 23.2% | 24 | Localized |
| BrCa | 27.5% | 35 | Localized |
| Psoriasis | 41.1% | 0 | Global |
| CRC | 42.9% | 0 | Global |
| AD | 68.3% | OOM | Predicted global |

---

## v0.3.2 comparison
- v0.3.2 curated: 8 core genes, 8 survived
- v0.3.3.2 curated: 8 core genes, 8 survived
- Delta: 0 (exact match)

---

## Provenance (SOP v1.0.3)
- `exact_commands.txt` — captured at execution time
- `environment.txt` — system state at execution time
- `pip_freeze.txt` — exact package versions at execution time

## Manifest status
Generated on 2026-05-05. Curated and blind numbers verified against committed
output files. Tier 3 stability pending.
