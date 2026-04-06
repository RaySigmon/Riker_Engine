# Riker Engine — Data Dependencies Manifest

Generated: 2026-04-06
Source: Harvested from all configs across riker_project/, .riker/runs/, and repo configs

## ASD (Autism Spectrum Disorder)
**Seed genes:** SFARI Gene database (sfari_genes.csv, 1,267 genes)
**Tissue:** Brain cortex
**GEO Datasets:**
| GEO ID | Role | Platform |
|--------|------|----------|
| GSE28521 | discovery | Illumina HumanRef-8 v3.0 |
| GSE26415 | discovery | Illumina HumanHT-12 V3.0 |
| GSE18123-GPL570 | discovery | Affymetrix HG-U133 Plus 2.0 |
| GSE18123-GPL6244 | discovery | Affymetrix HuGene-1.0-ST |
| GSE64018 | discovery | RNA-seq (FPKM) |
| GSE28475 | replication | Affymetrix HG-U133 Plus 2.0 |
| GSE42133 | replication | Illumina HumanHT-12 V4.0 |
| GSE89594 | replication | Affymetrix HuGene-2.0-ST |
| GSE102741 | replication | RNA-seq (RPKM) |

## T2D (Type 2 Diabetes)
**Seed genes:** Open Targets platform (443 genes)
**Tissue:** Pancreatic islets
**GEO Datasets:** See `results/t2d_test_run_1/curated/config.yaml` for full dataset list

## IBD (Inflammatory Bowel Disease)
**Seed genes:** Open Targets platform (762 genes)
**Tissue:** Intestinal mucosa
**GEO Datasets:**
| GEO ID | Role | Platform |
|--------|------|----------|
| GSE75214 | discovery | Affymetrix HuGene-1.0-ST |
See `results/ibd_test_run_1/curated/config.yaml` for full list

## AD (Alzheimer's Disease)
**Seed genes:** Assembled from GVC tables + Open Targets (801 genes)
**Tissue:** Brain cortex
**GEO Datasets:** 5 brain cortex datasets (see results/ad/ for assembly details)

## Breast Cancer
**Seed genes:** Custom curated (653 genes)
**Tissue:** Breast tumor
**GEO Datasets:**
| GEO ID | Role | Platform |
|--------|------|----------|
| GSE10810 | discovery | Affymetrix HG-U133 Plus 2.0 (GPL570) |
| GSE42568 | discovery | Affymetrix HG-U133 Plus 2.0 (GPL570) |
| GSE15852 | discovery | Affymetrix HG-U133A (GPL96) |
| GSE45827 | replication | Affymetrix HG-U133 Plus 2.0 (GPL570) |
| GSE65194 | replication | Affymetrix HG-U133 Plus 2.0 (GPL570) |

## Shared Dependencies
- **HGNC complete set:** `hgnc_complete_set.txt` (~17 MB) from genenames.org
- **Pathway databases:** KEGG, Reactome, MSigDB Hallmark gene sets
