# Disease Day Manifest — Breast Cancer

**Disease:** Breast Cancer (BrCa)
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
| Curated | `breast_cancer_curated_genes.csv` (653 genes) | `f21ac606583af5fee09c59da13f54814b4b6526246fe246ba282b9f44e6c210e` |
| Blind | `all_protein_coding_genes.csv` (19,296 genes) | `e9ab2790c45a3501545d60298a0fc51c7cd2791b040002387665f671c98639aa` |

---

## GEO datasets used

| Accession | Role | Tissue | Platform | Description |
|-----------|------|--------|----------|-------------|
| GSE10810 | discovery | breast | GPL570 | Breast tumor vs healthy |
| GSE42568 | discovery | breast | GPL570 | Breast cancer vs normal |
| GSE15852 | discovery | breast | GPL96 | Breast tumor vs normal tissue |
| GSE45827 | replication | breast | GPL570 | Breast cancer diagnosis |
| GSE65194 | replication | breast | GPL570 | TNBC + HER2 + luminal vs healthy |

3 discovery + 2 replication. All breast tissue. Mixed platforms (GPL570 + GPL96).

---

## Run results

### Tier 1 — Curated BrCa run (publishable validation chain)
- Wall clock: ~2 min
- Phase 1 study genes: 227
- Phase 4 core genes: 153
- Phase 4 significant clusters: 2
- Phase 5 survived: 141 (12 eliminated)
- Phase 6 meta-significant (random effects): 114

**Validation:** Close to v0.3.2 baseline (152/139). +1 core, +2 survived.

### Tier 2 — Blind protein-coding run (publishable — localized regime)
- Wall clock: ~16 min
- Phase 1 study genes: 5,315 (27.5% yield)
- Phase 4 core genes: 3,697
- **Phase 4 significant clusters: 35** (highest of any disease — localized regime)
- Phase 5 survived: 3,380 (317 eliminated, 8.6% elimination rate)
- Phase 6 meta-significant (random effects): 2,555

**Key finding:** BrCa blind at 27.5% Phase 1 yield produces **35 Bonferroni-significant
clusters** — more than any other disease tested, including those with higher yields.
This result establishes that Phase 1 yield alone does not determine the regime;
biological subtype heterogeneity creates real localized cluster structure even at
moderate-to-high yields.

### Tier 3 — 50-run stability profile
- Wall clock: 41,472 sec (11.5 hours)
- Runs completed: 50/50 (0 failures)
- Total unique genes seen: 3,789
- Iron-clad (>=90%): 3,502
- Borderline (50-89%): 252
- Stochastic (<50%): 35
- Iron-clad fraction: 92.43%
- Pairwise Jaccard median: 0.9606 (p25=0.9580, p75=0.9635)
- Core gene counts: min=3,644, max=3,715, mean=3,679.2
- Clusters with >=5 iron-clad members: 347

---

## Molecular subtype reconstruction

The engine independently separated known breast cancer molecular subtypes into
distinct iron-clad gene clusters **without being told subtypes exist**:

| Subtype program | Key markers found | Cluster distribution |
|----------------|-------------------|---------------------|
| ER+/Luminal | ESR1, FOXA1, GATA3, TFF1, TFF3, XBP1, AGR2 | 6 clusters |
| HER2 | ERBB2, PGAP3 | 2 clusters |
| Proliferation | TOP2A, AURKA, BIRC5, CDK1, CCNB1, BUB1B | 5 clusters |
| Basal/TNBC | EGFR, KRT14 | 2 clusters |
| Immune infiltrate | CXCL9, CXCL10 | 2 clusters |
| ECM/Stromal | FN1, FAP, MMP11 | 3 clusters |

Markers distribute across multiple consensus clusters (not one mega-cluster per
subtype), consistent with the SOP's binding terminology: biological annotations
span multiple computational clusters.

This reproduces the canonical molecular subtype classification framework from raw
expression data — the same classification that took the field decades to establish
through clinical assays, survival data, and immunohistochemistry.

---

## Regime classification

BrCa is **localized** — the strongest localized result in the validation, right
at the regime threshold edge:

| Disease | Phase 1 yield | Sig clusters | Regime |
|---------|--------------|-------------|--------|
| T2D | 8.0% | 16 | Localized |
| ASD | 9.3% | 15 | Localized |
| IPF | 23.2% | 24 | Localized |
| **BrCa** | **27.5%** | **35** | **Localized (strongest)** |
| Psoriasis | 41.1% | 0 | Global |
| CRC | 42.9% | 0 | Global |

---

## v0.3.2 comparison (curated only)
- v0.3.2 curated: 152 core genes, 139 survived
- v0.3.3.2 curated: 153 core genes, 141 survived
- Delta: +1 core, +2 survived

---

## Provenance (SOP v1.0.3)
- `exact_commands.txt` — captured at execution time
- `environment.txt` — system state at execution time
- `pip_freeze.txt` — exact package versions at execution time

## Manifest status
Generated on 2026-05-06 from pipeline outputs. All numbers verified against
committed output files.
