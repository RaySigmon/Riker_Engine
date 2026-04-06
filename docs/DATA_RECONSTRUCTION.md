# Data Reconstruction Guide

Some datasets used in the Riker Engine validation cannot be directly downloaded as
standard GEO series matrix files. They require either reconstruction from supplementary
RNA-seq data or filtering to specific sample subsets.

This document provides step-by-step instructions for each.

## RNA-seq Datasets (Reconstruction Required)

These datasets provide count or FPKM data in supplementary files rather than standard
series matrix format. The Riker Engine expects a genes-by-samples matrix in series
matrix format, so these must be reconstructed.

### GSE64018 — ASD brain cortex RNA-seq (Gupta et al.)

**What you need:** The supplementary file `GSE64018_adjfpkm.txt.gz` contains adjusted
FPKM values with Ensembl gene IDs.

1. Download: `wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE64nnn/GSE64018/suppl/GSE64018_adjfpkm.txt.gz"`
2. The file has Ensembl IDs as row names. Use `data/mappings/ensembl_to_symbol.txt` to
   convert to gene symbols.
3. Log2-transform: values are FPKM, apply log2(FPKM + 1).
4. Format as a series matrix: gene symbols as rows, sample IDs as columns.
5. Save as `GSE64018_reconstructed_series_matrix.txt.gz` in `data/geo/asd/`.
6. Use `data/mappings/ensembl_to_symbol.txt` as the platform file in your config.

**Phenotype:** `Sample_characteristics_ch1` — case: `"diagnosis: asd"`, control: `"diagnosis: ctl"`

### GSE102741 — ASD brain cortex RNA-seq

**What you need:** Supplementary file `GSE102741_log2RPKMcounts.xlsx`

1. Download: `wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102741/suppl/GSE102741_log2RPKMcounts.xlsx"`
2. Convert Excel to TSV. Values are already log2 RPKM.
3. Map Ensembl IDs to symbols using `data/mappings/ensembl_to_symbol.txt`.
4. Format as series matrix and save as `GSE102741_reconstructed_series_matrix.txt.gz`.
5. Use `data/mappings/ensembl_to_symbol.txt` as the platform file.

**Phenotype:** `Sample_characteristics_ch1` — case: `"disease status: autism spectrum disorder"`,
control: `"disease status: healthy control"`

### GSE86468 — T2D pancreatic islet RNA-seq (Lawlor et al.)

**What you need:** Supplementary count data.

1. Download the series matrix: the standard download works, but the expression data
   is not in the matrix. Download supplementary counts instead.
2. Filter to baseline samples only (exclude stimulated conditions).
3. Apply log2(counts + 1) normalization.
4. Map Ensembl IDs to symbols.
5. Save as `GSE86468_reconstructed_series_matrix.txt.gz` in `data/geo/t2d/`.

**Phenotype:** case: `"disease: type 2 diabetic"`, control: `"disease: non-diabetic"`

## Filtered Datasets

These datasets download normally from GEO but contain samples from multiple conditions
or brain regions. They must be filtered to the relevant subset before use.

### GSE33000 — AD prefrontal cortex (must remove Huntington's samples)

The full dataset contains Alzheimer's, Huntington's, and non-demented samples.

1. Download the full series matrix (the download script handles this).
2. Filter to keep only AD and non-demented samples. Remove all HD samples.
3. Save as `GSE33000_AD_only_series_matrix.txt.gz`.

**Phenotype field:** `Sample_characteristics_ch2` — case: `"disease status: alzheimer"`,
control: `"disease status: non-demented"`

**Platform:** GPL4372 (Rosetta/Merck log-ratio array — values centered near zero,
not log2 intensities)

### GSE118553 — AD temporal cortex (must extract TC subset)

The full dataset contains multiple brain regions and AsymAD samples.

1. Download the full series matrix.
2. Filter to temporal cortex samples only.
3. Exclude AsymAD (asymptomatic AD) samples — keep only definite AD and controls.
4. Save as `GSE118553_TC_series_matrix.txt.gz`.

**Phenotype field:** `Sample_characteristics_ch1` — case: `"disease state: ad"`,
control: `"disease state: control"`

**Platform:** GPL10558

### GSE5281 — AD superior frontal gyrus (must extract SFG subset)

The full dataset contains multiple brain regions.

1. Download the full series matrix.
2. Filter to superior frontal gyrus (SFG) samples only.
3. Save as `GSE5281_SFG_series_matrix.txt.gz`.

**Phenotype field:** `Sample_characteristics_ch1` — case: `"disease state: alzheimer"`,
control: `"disease state: normal"`

**Platform:** GPL570

## Verification

After reconstruction/filtering, verify each dataset:

1. Check that the matrix has the expected number of samples (compare with GEO listing).
2. Check expression value ranges: log2-transformed data should typically be in the
   range of 0–16 for microarray, 0–20 for RNA-seq log2(counts+1).
3. Confirm case/control labels parse correctly by running:
   ```bash
   riker validate-config your_config.yaml
   ```

## File Naming Convention

Place reconstructed files alongside the standard downloads:
```
data/geo/asd/
├── GSE28521_series_matrix.txt.gz          # standard download
├── GSE64018_reconstructed_series_matrix.txt.gz  # reconstructed
└── ...
```

Update your config YAML paths to point to these files.
