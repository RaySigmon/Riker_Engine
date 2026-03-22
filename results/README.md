# Riker Engine — Validation Results

**Engine version:** v0.2.0
**Date:** 2026-03-21
**Platform:** Raspberry Pi 5 (8GB), Python 3.13, Claude Code CLI

Three diseases tested with zero code modifications between runs.
Each condition includes a **curated seed run** (disease-specific gene list)
and a **blind all-genes run** (every expressed gene as seeds).

---

## Cross-Disease Summary

| | ASD | T2D | IBD |
|--|-----|-----|-----|
| **Tissue** | Brain cortex | Pancreatic islets | Intestinal mucosa |
| **Seed source** | SFARI (1,267) | Open Targets (443) | Open Targets (762) |
| **Discovery datasets** | 3 | 3 | 3 |
| **Replication datasets** | 4 | 1 | 3 |
| **Phase 1 study genes** | 141 (11.1%) | 56 (12.6%) | 408 (53.5%) |
| **Phase 4 core genes** | 35 | 8 | 304 |
| **Phase 5 survived** | 35/35 | 8/8 | 302/304 |
| **Phase 6 significant** | 13/35 | 8/8 | 296/302 |
| **QC status** | PASSED | PASSED | PASSED |
| **Fold changes valid** | Yes | Yes | Yes |

## All-Genes Blind Validation

| | ASD | T2D | IBD |
|--|-----|-----|-----|
| **All-genes seeds** | 39,185 | 26,800 | 24,651 |
| **Study genes** | 1,813 | 1,605 | 9,120 |
| **Core genes** | 415 | 177 | 6,458 |
| **Curated core recovery** | 27/35 (77.1%) | 5/8 (62.5%) | 297/304 (97.7%) |
| **All pass Phase 1?** | 35/35 (100%) | 8/8 (100%) | 304/304 (100%) |

Every curated core gene passes Phase 1 in the hypothesis-free run across
all three diseases. Genes that drop from core status do so because of
stricter FDR correction (larger denominator) or clustering noise at scale —
never because the biology fails.

---

## ASD (Autism Spectrum Disorder)

**Directory:** `asd_test_run_1/`

### Datasets
| ID | Samples | Platform | Role | Tissue |
|----|---------|----------|------|--------|
| GSE28521 | 79 | GPL6883 | Discovery | Brain cortex |
| GSE28475 | 123 | GPL6883 | Discovery | Brain cortex |
| GSE64018 | 24 | GPL11154 (RNA-seq) | Discovery | Brain cortex |
| GSE102741 | 52 | GPL11154 (RNA-seq) | Replication | Brain DLPFC |
| GSE18123 | 186 | GPL6244 | Replication | Blood (LCL) |
| GSE26415 | 84 | GPL6480 | Replication | Blood |
| GSE42133 | 147 | GPL10558 | Replication | Blood |

### Key Results
- 35 core genes across 8 clusters, all surviving replication
- Top genes: ZFHX3 (log2FC=0.582, up), TBR1 (-0.572, down), ATP2B2 (-0.161, down), SEZ6L2 (-0.274, down)
- 13/35 significant by random effects meta-analysis
- All-genes run recovered 27/35 (77.1%) — 8 dropped due to 31x stricter FDR denominator

### Notes
- RNA-seq datasets (GSE64018, GSE102741) required reconstruction from supplementary files
- GSE28475 had mixed tissue types; only diagnosis-labeled samples used
- GSE18123 matched only "AUTISM" diagnosis (PDD-NOS, Asperger's excluded)

---

## T2D (Type 2 Diabetes)

**Directory:** `t2d_test_run_1/`

### Datasets
| ID | Samples | Platform | Role | Tissue |
|----|---------|----------|------|--------|
| GSE41762 | 77 (20 T2D + 43 ND + 14 unlabeled) | GPL6244 | Discovery | Pancreatic islets |
| GSE25724 | 13 | GPL96 | Discovery | Pancreatic islets |
| GSE20966 | 20 | GPL1352 | Discovery | Pancreatic beta cells (LCM) |
| GSE86468 | 8 (Baseline only) | GPL18573 (RNA-seq) | Replication | Pancreatic islets |

### Key Results
- 8 core genes in 1 significant cluster, ALL downregulated in T2D islets
- **SLC2A2** (GLUT2): log2FC=-1.300 — primary beta cell glucose transporter
- **ABCC8** (sulfonylurea receptor): log2FC=-0.863 — known drug target
- **CACNA1D**: log2FC=-0.870 — voltage-gated calcium channel for insulin secretion
- **CASR**: log2FC=-0.392 — calcium-sensing receptor regulating insulin release
- All 8/8 significant by random effects meta-analysis
- All-genes run recovered 5/8 (62.5%)

### Notes
- GSE159984 returned 404 (not available on FTP)
- GSE86468 RNA-seq required reconstruction; only Baseline samples used (3 T2D, 5 ND)
- GWAS-top genes (TCF7L2, KCNJ11, GCK) did not pass Phase 1 — GWAS identifies risk variants, not DE genes
- Phase 5 tissue limitation: "islet" doesn't trigger brain/blood elimination logic

---

## IBD (Inflammatory Bowel Disease)

**Directory:** `ibd_test_run_1/`

### Datasets
| ID | Samples | Platform | Role | Tissue |
|----|---------|----------|------|--------|
| GSE75214 | 194 | GPL6244 | Discovery | Colonic + ileal mucosa |
| GSE16879 | 133 | GPL570 | Discovery | Colonic mucosa |
| GSE59071 | 116 | GPL6244 | Discovery | Colonic mucosa |
| GSE87466 | 108 | GPL13158 | Replication | Colonic mucosa |
| GSE38713 | 43 | GPL570 | Replication | Colonic mucosa |
| GSE36807 | 35 | GPL570 | Replication | Colonic tissue |

### Key Results
- 304 core genes across 50 clusters — strongest signal of any disease tested
- 53.5% of seed genes passed Phase 1 (vs 11-13% for ASD and T2D)
- 302/304 survived replication, 296 significant by random effects
- Top upregulated: MUC1 (2.014), CCL2 (2.212), CCL20 (1.639), IL7R (1.700), PDE4B (1.555)
- Top downregulated: ABCB1 (-2.388), SLC22A5 (-1.710), CDKN2B (-1.459), ALPI (-1.411)
- Known IBD genes confirmed: NOD2, IL10, JAK2, STAT3, TNFAIP3, TLR4, ICAM1, IFNG, HNF4A
- All-genes run recovered 297/304 (97.7%) — highest recovery of any disease

### Notes
- Both Crohn's Disease and Ulcerative Colitis treated as cases
- GSE16879 includes pre/post-infliximab samples; all CD/UC samples used as cases
- GSE38713 phenotype encoded in Sample_source_name_ch1 (not characteristics)
- IL23R, ATG16L1, IRGM, CARD9 did not pass Phase 1 — risk variants, not expression changes
- Phase 5 tissue limitation: "colon" doesn't trigger brain/blood elimination logic

---

## Directory Structure

```
results/
├── README.md
├── asd_test_run_1/
│   ├── curated/          # 1,267 SFARI seed genes
│   │   ├── config.yaml
│   │   ├── pipeline_summary.json
│   │   ├── qc_report.json
│   │   ├── phase1_study_genes.csv
│   │   ├── phase4_core_genes.csv
│   │   ├── phase4_all_levels.csv
│   │   ├── phase5_verdicts.csv
│   │   └── phase6_meta_analysis.csv
│   └── allgenes/         # 39,185 expressed genes as seeds
│       └── (same files)
├── t2d_test_run_1/
│   ├── curated/          # 443 Open Targets seed genes
│   └── allgenes/         # 26,800 expressed genes as seeds
└── ibd_test_run_1/
    ├── curated/          # 762 Open Targets seed genes
    └── allgenes/         # 24,651 expressed genes as seeds
```

---

## Reproducibility

All runs used:
- Riker Engine v0.2.0 (`git checkout v0.2.0`)
- `phase1.p_threshold: 0.05`, `phase1.min_datasets: 2`
- Default Phase 3 parameters: 3 n_neighbors x 5 seeds = 15 configurations
- Default Phase 4: 10,000 permutations, seed 42
- HGNC complete set for gene symbol resolution
- Platform annotation files from GEO FTP (`.annot` format)
- YAML configs included in each run directory

Hardware: Raspberry Pi 5 (8GB RAM), Debian Bookworm, Python 3.13.
Run times: 1-2 min (curated), 5-70 min (all-genes, depending on study gene count).
