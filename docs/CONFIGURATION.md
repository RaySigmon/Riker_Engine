# Configuration Reference

Riker Engine is driven by a single YAML configuration file. To study a different
condition you only need to change the seed genes and datasets -- the pipeline
logic stays the same. This document describes every field the config loader
(`riker/config.py`) recognises.

---

## 1. Top-Level Fields

| Field | Required | Default | Description |
|-------|----------|---------|-------------|
| `condition` | **yes** | -- | Short label for the condition (e.g. `ASD`, `T2D`, `alzheimers_disease`). Used in output filenames and logs. |
| `seed_genes` | **yes** | -- | Path to a CSV file containing known disease genes. Must have a column with HGNC symbols. |
| `hgnc_path` | no | `auto` | Path to the HGNC complete set file. When set to `auto` the engine downloads it on first run and caches it under `~/.riker/`. |
| `output_dir` | no | `riker_output` | Directory where all pipeline outputs are written. Created if it does not exist. |
| `random_seed` | no | `42` | Base random seed used throughout the pipeline for reproducibility. |

---

## 2. Dataset Fields

Datasets are listed under the `datasets` key. At least one dataset with
`role: discovery` is required.

```yaml
datasets:
  - id: GSE28521
    series_matrix: /path/to/GSE28521_series_matrix.txt.gz
    platform: /path/to/GPL6883_clean.annot
    role: discovery
    tissue: brain
    phenotype_field: "Sample_characteristics_ch1"
    case_values: ["disease status: autism"]
    control_values: ["disease status: control"]
```

| Field | Required | Default | Description |
|-------|----------|---------|-------------|
| `id` | **yes** | -- | Unique identifier for the dataset, typically a GEO accession (e.g. `GSE28521`). |
| `series_matrix` | **yes** | -- | Path to the GEO series matrix file (`.txt` or `.txt.gz`). |
| `platform` | **yes** | -- | Path to the platform annotation file mapping probes to gene symbols. |
| `role` | no | `discovery` | Either `discovery` or `replication`. Discovery datasets are used in Phase 1 differential expression. Replication datasets provide independent validation. |
| `tissue` | no | `brain` | Tissue of origin (e.g. `brain`, `blood`, `islet`). Used for grouping in reports. |
| `phenotype_field` | no | auto-detect | Metadata field that encodes case/control status. See the Phenotype Override Guide below. |
| `case_values` | no | auto-detect | List of values in `phenotype_field` that identify **case** samples. |
| `control_values` | no | auto-detect | List of values in `phenotype_field` that identify **control** samples. |
| `gene_column` | no | auto-detect | Column name in the platform annotation file containing gene symbols. Use when auto-detection picks the wrong column. |
| `probe_column` | no | auto-detect | Column name in the platform annotation file containing probe IDs. Use when auto-detection picks the wrong column. |

---

## 3. Phenotype Override Guide

The engine tries to auto-detect which metadata field encodes disease status. When
auto-detection fails -- or when a dataset uses non-obvious labels -- provide
`phenotype_field`, `case_values`, and `control_values` explicitly.

### Which metadata field to use

GEO series matrices store sample metadata in rows prefixed with
`!Sample_characteristics_ch1`, `!Sample_characteristics_ch2`,
`!Sample_source_name_ch1`, and so on. Open the series matrix (or check GEO
online) and find the row that distinguishes cases from controls.

| Situation | `phenotype_field` value |
|-----------|------------------------|
| Disease status in channel 1 characteristics | `Sample_characteristics_ch1` |
| Disease status in channel 2 characteristics (dual-channel arrays) | `Sample_characteristics_ch2` |
| Disease status in sample source name | `Sample_source_name_ch1` |

### How matching works

Values are matched **case-insensitively** and by **substring**. You only need to
supply the distinguishing portion of the value. For example, given the raw
metadata value `"disease status: Alzheimer's Disease"`, a `case_values` entry
of `"disease status: alzheimer"` will match.

### Common encoding patterns

Different datasets encode case/control status in many ways. Below are patterns
seen across ASD, T2D, and Alzheimer's disease configs:

```yaml
# Explicit disease status
case_values: ["disease status: autism"]
control_values: ["disease status: control"]

# Short diagnostic codes
case_values: ["diagnosis: asd"]
control_values: ["diagnosis: ctl"]

# Donor status (T2D islet studies)
case_values: ["status: diabetic donor"]
control_values: ["status: non-diabetic donor"]

# Disease state with brief labels (AD)
case_values: ["disease: a"]
control_values: ["disease: n"]

# Multiple case subtypes matched to one group
case_values: ["diagnosis: autism", "diagnosis: pdd-nos", "diagnosis: asperger"]
control_values: ["diagnosis: control"]
```

### Tips

- Always quote values that contain colons.
- If a dataset has subtypes you want to include as cases (e.g. PDD-NOS and
  Asperger's alongside Autism), list them all in `case_values`.
- Samples that match neither `case_values` nor `control_values` are silently
  excluded.

---

## 4. Phase Parameters

Phase parameters are optional. Defaults work well for most conditions.

### phase1 -- Differential Expression

```yaml
phase1:
  p_threshold: 0.05
  min_datasets: 2
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `p_threshold` | `0.05` | Significance threshold for per-dataset differential expression. |
| `min_datasets` | `2` | A gene must be significant in at least this many discovery datasets to pass Phase 1. |

### phase3 -- UMAP + HDBSCAN Clustering

```yaml
phase3:
  n_neighbors: [10, 15, 30]
  seeds: [42, 123, 456, 789, 1024]
  min_cluster_size: 5
  min_samples: 3
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_neighbors` | `[10, 15, 30]` | List of UMAP neighbor counts. The pipeline runs clustering at each value and takes the consensus. |
| `seeds` | `[42, 123, 456, 789, 1024]` | Random seeds for UMAP. Multiple seeds stabilize the embedding. |
| `min_cluster_size` | `5` | HDBSCAN minimum cluster size. |
| `min_samples` | `3` | HDBSCAN minimum samples parameter. |

### phase4 -- Permutation Testing

```yaml
phase4:
  n_permutations: 10000
  seed: 42
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_permutations` | `10000` | Number of permutations for significance testing. |
| `seed` | `42` | Random seed for the permutation procedure (independent of `random_seed`). |

---

## 5. Full Annotated Example

A minimal two-dataset config for a hypothetical condition:

```yaml
# ---- General ----
condition: T2D
seed_genes: /home/user/t2d/data/seed_genes/t2d_seed_genes.csv
hgnc_path: /home/user/.riker/hgnc_complete_set.txt   # or "auto"
output_dir: /home/user/t2d/output
random_seed: 42

# ---- Datasets ----
datasets:
  # Discovery: pancreatic islets, microarray
  - id: GSE41762
    series_matrix: /home/user/t2d/data/GSE41762_series_matrix.txt.gz
    platform: /home/user/t2d/data/GPL6244.annot
    role: discovery
    tissue: islet
    phenotype_field: Sample_characteristics_ch1
    case_values: ["status: diabetic donor"]
    control_values: ["status: non-diabetic donor"]

  # Discovery: pancreatic islets, different platform
  - id: GSE25724
    series_matrix: /home/user/t2d/data/GSE25724_series_matrix.txt.gz
    platform: /home/user/t2d/data/GPL96.annot
    role: discovery
    tissue: islet
    phenotype_field: Sample_characteristics_ch1
    case_values: ["disease state: type 2 diabetes"]
    control_values: ["disease state: non-diabetic"]

  # Replication: pancreatic islets, RNA-seq
  - id: GSE86468
    series_matrix: /home/user/t2d/data/GSE86468_reconstructed_series_matrix.txt.gz
    platform: /home/user/t2d/data/ensembl_to_symbol.txt
    role: replication
    tissue: islet
    phenotype_field: Sample_characteristics_ch1
    case_values: ["disease: type 2 diabetic"]
    control_values: ["disease: non-diabetic"]

# ---- Phase Parameters (all optional, defaults shown) ----
phase1:
  p_threshold: 0.05    # per-dataset significance cutoff
  min_datasets: 2       # gene must be DE in >= N discovery datasets

phase3:
  n_neighbors: [10, 15, 30]            # UMAP neighbor sweep
  seeds: [42, 123, 456, 789, 1024]     # UMAP random seeds
  min_cluster_size: 5                   # HDBSCAN min cluster size
  min_samples: 3                        # HDBSCAN min samples

phase4:
  n_permutations: 10000   # permutation count for significance
  seed: 42                 # permutation RNG seed
```
