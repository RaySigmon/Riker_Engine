# Riker Engine

**A condition-agnostic transcriptomics pipeline for discovering replicated gene modules from public expression data.**

[![Tests](https://img.shields.io/badge/tests-280%2B%20passing-brightgreen)]()
[![Python](https://img.shields.io/badge/python-3.11%2B-blue)]()
[![License](https://img.shields.io/badge/license-MIT-green)]()

---

## What Is This?

The Riker Engine takes a list of candidate genes for any condition (autism, Alzheimer's, Parkinson's — anything with public transcriptomic data) and systematically determines which genes form reproducible, statistically defensible modules across multiple independent datasets.

It was born from [Project Riker](https://github.com/RaySigmon/Project-Riker), a proof-of-concept study that discovered replicated gene modules in autism spectrum disorder using a d20-based narrative framework. The engine generalizes that methodology into a reusable tool.

### The Problem It Solves

Candidate gene lists from GWAS, literature reviews, or database mining are noisy. Most genes on any list won't replicate across datasets. Traditional approaches test genes one at a time, miss coordinated patterns, and are vulnerable to p-hacking. The Riker Engine addresses all three problems through a six-phase progressive filtering pipeline with built-in anti-fishing protections.

### Key Principles

- **Agnostic Discovery**: Clusters are found by the data, not imposed by the researcher. No biological hypotheses during clustering.
- **Progressive Honesty**: Each phase applies stricter statistical filters. Genes that can't survive don't advance.
- **Pre-Specification Protocol**: The core gene list is locked *before* replication data is touched. No post-hoc additions or removals.
- **Full Seed Set FDR**: Benjamini-Hochberg correction uses the entire seed gene count as denominator, not just the study set. This prevents the inflated significance that killed the AD stress test.

---

## Pipeline Overview

```
Seed Genes + GEO Datasets
        │
        ▼
┌─────────────────────────┐
│  Phase 1: Cross-Ref     │  Welch's t-test per gene per dataset
│  p<0.05 in ≥2 datasets  │  → Study gene set
└─────────┬───────────────┘
          ▼
┌─────────────────────────┐
│  Phase 2: Pathways      │  KEGG + Reactome + Hallmark mapping
│  Feature matrix build   │  → Binary pathways + expression stats
└─────────┬───────────────┘
          ▼
┌─────────────────────────┐
│  Phase 3: Clustering    │  UMAP + HDBSCAN × 15 configurations
│  Consensus matrix       │  → Gene modules discovered
└─────────┬───────────────┘
          ▼
┌─────────────────────────┐
│  Phase 4: Robustness    │  Permutation + sensitivity + LOO
│  Core gene ID           │  → Pre-specified replication set
└─────────┬───────────────┘
          ▼
┌─────────────────────────┐
│  Phase 5: Replication   │  Independent directional testing
│  Elimination protocol   │  → Survivors + eliminated genes
└─────────┬───────────────┘
          ▼
┌─────────────────────────┐
│  Phase 6: Meta-Analysis │  IVW fixed + random effects
│  Forest plot data       │  → Final effect sizes + I²
└─────────┴───────────────┘
          │
          ▼
    Structured Output
    (CSV + JSON + QC Report)
```

---

## Installation

### Requirements

- Python 3.11 or higher
- A system with at least 2GB RAM (tested on Raspberry Pi 5)

### Install from Source

```bash
git clone https://github.com/RaySigmon/Riker_Engine.git
cd Riker_Engine

# Core install (Phases 1-2, 4-6, all stats)
pip install -e .

# Full install with clustering dependencies (Phase 3)
pip install -e ".[clustering]"
```

### Verify Installation

```bash
# Run the full test suite
python -m pytest tests/ -q

# Check the CLI works
riker --help
```

### Dependencies

**Core** (installed automatically):
- numpy, pandas, scipy — numerical computing
- scikit-learn — HDBSCAN clustering, utilities
- PyYAML — configuration parsing
- matplotlib — plotting (optional, for forest plots)

**Clustering extras** (`pip install -e ".[clustering]"`):
- umap-learn — UMAP dimensionality reduction
- hdbscan — standalone HDBSCAN (fallback; sklearn's built-in is preferred)

**snRNA-seq extras** (`pip install scanpy`):
- scanpy — AnnData/h5ad file support for single-nucleus data

---

## Quick Start

### 1. Prepare Your Data

You need three things:

**a) Seed gene list** — A CSV file with candidate genes for your condition:
```csv
gene_symbol,source,evidence
SHANK3,SFARI,syndromic
CHD8,SFARI,high_confidence
SCN2A,SFARI,high_confidence
DYRK1A,literature,replicated
...
```
The only required column is `gene_symbol`. Additional columns are preserved but not used by the engine.

**b) GEO series matrix files** — Download from [GEO](https://www.ncbi.nlm.nih.gov/geo/). Each dataset needs:
- The series matrix file (e.g., `GSE28521_series_matrix.txt.gz`)
- The platform annotation file (e.g., `GPL6244.txt`)

**c) YAML configuration file** — Tells the engine what to run:

```yaml
# riker_asd.yaml — Example ASD configuration
condition: ASD
seed_genes: data/asd_seed_genes.csv
hgnc_path: auto  # downloads HGNC automatically
output_dir: output/asd_run

datasets:
  # Discovery datasets (used for Phases 1-4)
  - id: GSE28521
    series_matrix: data/GSE28521_series_matrix.txt.gz
    platform: data/GPL6244.txt
    role: discovery

  - id: GSE38322
    series_matrix: data/GSE38322_series_matrix.txt.gz
    platform: data/GPL6244.txt
    role: discovery

  - id: GSE33000
    series_matrix: data/GSE33000_series_matrix.txt.gz
    platform: data/GPL4372.txt
    role: discovery
    # Override phenotype extraction for this dataset
    phenotype_field: Sample_characteristics_ch1
    case_values: ["alzheimer"]
    control_values: ["control"]

  # Replication datasets (used in Phase 5 only)
  - id: GSE64018
    series_matrix: data/GSE64018_series_matrix.txt.gz
    platform: data/GPL11154.txt
    role: replication
    tissue: brain

  - id: GSE18123
    series_matrix: data/GSE18123_series_matrix.txt.gz
    platform: data/GPL570.txt
    role: replication
    tissue: blood

# Phase parameters (all optional — defaults shown)
phase1:
  p_threshold: 0.05
  min_datasets: 2

phase3:
  n_neighbors: [10, 15, 30]
  seeds: [42, 123, 456, 789, 1024]
  min_cluster_size: 5
  min_samples: 3

phase4:
  n_permutations: 10000
  seed: 42
```

### Using snRNA-seq Data

For single-nucleus RNA-seq data, pseudo-bulk by donor before running the pipeline:

```python
from riker.ingestion.snrnaseq import pseudo_bulk_from_counts, pseudo_bulk_from_h5ad

# From an h5ad file (AnnData format)
result = pseudo_bulk_from_h5ad(
    "data/asd_snrnaseq.h5ad",
    donor_key="donor_id",
    condition_key="disease",
    cell_type_key="cell_type",
    target_cell_type="Excitatory",
    case_values=["ASD"],
    control_values=["Control"],
)

# result.expression is a genes × samples DataFrame
# result.phenotypes is a {sample_id: 'case'/'control'} dict
# Both plug directly into Phase 1 cross-referencing

# Or from a CSV count matrix
import pandas as pd
counts = pd.read_csv("data/counts_with_metadata.csv")
result = pseudo_bulk_from_counts(
    counts,
    donor_col="donor_id",
    condition_col="condition",
    cell_type_col="cell_type",
    target_cell_type="Excitatory",
    case_values=["ASD"],
    control_values=["Control"],
)
```

The pseudo-bulking module:
- Aggregates nuclei per donor using sum (recommended for DE analysis)
- Applies CPM + log2(CPM+1) normalization
- Filters donors with fewer than 10 nuclei (configurable)
- Filters nuclei with fewer than 200 detected genes (QC)
- Auto-detects case/control from condition labels
- Outputs a genes × samples matrix that feeds directly into Phase 1

### 2. Validate Your Config

```bash
riker validate riker_asd.yaml
```

This checks that the YAML is well-formed, all required fields are present, and at least one discovery dataset is specified. It does NOT check that data files exist (those are checked at runtime).

### 3. Run the Pipeline

```bash
riker run riker_asd.yaml
```

The engine will:
1. Load seed genes and resolve HGNC symbols
2. Parse each GEO dataset (series matrix + platform annotation)
3. Auto-detect and apply log2 normalization if needed
4. Run all 6 phases with QC gates between each
5. Write structured output to the configured directory

### 4. Check Your Output

```
output/asd_run/
├── pipeline_summary.json      # Overall results + locked gene list
├── qc_report.json             # QC checks with pass/warn/critical
├── phase1_study_genes.csv     # Genes passing Phase 1 filter
├── phase4_core_genes.csv      # Core genes (pre-specified set)
├── phase5_verdicts.csv        # Replication results per gene
└── phase6_meta_analysis.csv   # Forest plot data + effect sizes
```

---

## Configuration Reference

### Top-Level Fields

| Field | Required | Default | Description |
|---|---|---|---|
| `condition` | Yes | — | Condition name (e.g., 'ASD', 'AD') |
| `seed_genes` | Yes | — | Path to seed gene CSV |
| `hgnc_path` | No | `auto` | Path to HGNC file, or 'auto' to download |
| `output_dir` | No | `riker_output` | Output directory |
| `random_seed` | No | `42` | Base random seed |

### Dataset Fields

| Field | Required | Default | Description |
|---|---|---|---|
| `id` | Yes | — | Unique dataset identifier |
| `series_matrix` | Yes | — | Path to GEO series matrix file |
| `platform` | Yes | — | Path to platform annotation file |
| `role` | No | `discovery` | `discovery` or `replication` |
| `tissue` | No | `brain` | `brain` or `blood` (replication only) |
| `phenotype_field` | No | auto-detect | Metadata field for case/control |
| `case_values` | No | auto-detect | Values indicating case samples |
| `control_values` | No | auto-detect | Values indicating control samples |

### Per-Dataset Phenotype Overrides

GEO metadata varies wildly between datasets. The engine auto-detects case/control assignments by scanning for keywords like "autism", "alzheimer", "disease", "control", "normal" in the `Sample_characteristics_ch1` field. When auto-detection fails (and it will for some datasets), use explicit overrides:

```yaml
- id: GSE33000
  phenotype_field: Sample_characteristics_ch1
  case_values: ["alzheimer disease"]
  control_values: ["normal"]
```

The engine performs case-insensitive substring matching on the values.

### Phase Parameters

**Phase 1 (Cross-Referencing)**:
- `p_threshold`: Significance cutoff (default: 0.05). Intentionally lenient.
- `min_datasets`: Gene must be significant in at least this many datasets (default: 2).

**Phase 3 (Clustering)**:
- `n_neighbors`: UMAP neighborhood sizes to sweep (default: [10, 15, 30]).
- `seeds`: Random seeds for UMAP (default: [42, 123, 456, 789, 1024]).
- `min_cluster_size`: HDBSCAN minimum cluster size (default: 5).
- `min_samples`: HDBSCAN minimum samples (default: 3).
- Total configurations = len(n_neighbors) × len(seeds). Default: 15.

**Phase 4 (Robustness)**:
- `n_permutations`: Permutations per cluster significance test (default: 10000).
- `seed`: Base random seed for permutation tests (default: 42).

---

## How Each Phase Works

### Phase 1: Cross-Referencing

For each seed gene in each discovery dataset, computes case vs control log2 fold change using Welch's t-test (exact t-distribution, not normal approximation). Genes significant (p < 0.05) in ≥2 datasets advance to the study gene set.

**QC checks**: Fold change range validation catches impossible values (the "log2FC = 119 bug" from development — caused by running stats on raw intensities instead of log2-transformed data). The normalizer auto-detects and applies log2(x+1) when median expression ≥ 20.

### Phase 2: Pathway Mapping

Maps study genes to biological pathways and builds a feature matrix for clustering. Each gene gets a feature vector combining:
- Binary pathway memberships (one column per pathway)
- Average log2FC across datasets
- -log10(minimum p-value), capped at 10
- Number of significant datasets
- Directional consistency score
- Confidence tier score

All features are min-max normalized to [0, 1].

**Anti-circularity rule**: Individual pathway IDs are features. Pre-assigned biological category labels (e.g., "chromatin remodeling") are NEVER used as features — that would create tautological enrichment.

### Phase 3: Consensus Clustering

Runs UMAP (2D, min_dist=0.0) + HDBSCAN across all configurations (3 n_neighbors × 5 seeds = 15 by default). Builds a co-association consensus matrix where entry (i,j) is the fraction of configurations placing genes i and j in the same cluster (excluding noise assignments). Final clusters come from HDBSCAN with `metric='precomputed'` on the distance matrix (1 - consensus).

This makes results robust to UMAP/HDBSCAN parameter choices. Genes that only cluster together in one specific configuration won't survive.

### Phase 4: Robustness Testing

Four progressive filters:

1. **Permutation significance**: Is each cluster's mean |log2FC| more extreme than random gene sets of the same size? Bonferroni-corrected across all clusters.

2. **Sensitivity analysis**: Tests genes at four levels:
   - Level 1: p < 0.05 in ≥2 datasets (baseline)
   - Level 2: p < 0.01 in ≥2 datasets
   - Level 3: FDR q < 0.10 (BH correction over FULL seed set)
   - Level 4: FDR q < 0.05 (BH correction over FULL seed set)

3. **Leave-one-dataset-out (LOO)**: Removes each dataset in turn and checks if ≥80% of cluster genes still meet criteria. Catches dataset-dependent clusters.

4. **Core gene identification**: Level 2 survivors with ≥3 cluster-mates become the pre-specified set for replication. **This list is locked and recorded.**

### Phase 5: Independent Replication

Tests each core gene in held-out replication datasets for directional concordance.

**Elimination rules**:
- **ELIMINATE**: Significant (p < 0.05) in the OPPOSITE direction in a brain dataset
- **RETAIN**: Non-significant opposite effects (noted but not eliminated)
- **RETAIN**: Blood non-replication does NOT trigger elimination (brain-specific signals are expected to fail in blood)

This is the protocol that caught KANK1 in the ASD proof-of-concept — it passed all discovery filters but showed significant upregulation in brain replication (opposite to discovery direction of downregulation).

### Phase 6: Effect Size Meta-Analysis

Computes inverse-variance weighted meta-analysis for each surviving gene:
- **Fixed effects**: Assumes homogeneous true effect across datasets
- **Random effects** (DerSimonian-Laird): Allows between-study heterogeneity — this is the primary result
- **Heterogeneity**: Cochran's Q, I², τ² statistics
- **Expression scale check**: Flags datasets with max expression < 5.0 (log-ratio format like GSE33000)

Output includes forest plot data (per-dataset effect sizes + SEs + pooled estimates).

---

## Output Files

### `pipeline_summary.json`

Top-level results:
```json
{
  "condition": "ASD",
  "phase1_study_genes": 42,
  "phase4_core_genes": 18,
  "phase5_survived": 15,
  "phase5_eliminated": 3,
  "phase6_significant_random": 12,
  "qc_status": "PASSED",
  "locked_core_genes": ["ATP2B2", "SEZ6L2", "RELN", ...]
}
```

### `phase4_core_genes.csv`

The pre-specified gene list (locked before replication):
```
gene,cluster_id,max_level_survived,mean_log2fc,direction
ATP2B2,0,3,-0.82,down
SEZ6L2,0,4,-0.71,down
...
```

### `phase5_verdicts.csv`

Per-gene replication results:
```
gene,cluster_id,status,discovery_direction,n_brain_concordant,n_brain_discordant,reason
ATP2B2,0,survived,down,2,0,"Brain replication: 2/2 concordant..."
KANK1,0,eliminated,down,0,1,"Significant opposite direction in brain..."
```

### `phase6_meta_analysis.csv`

Forest plot data:
```
gene,random_effect,random_se,random_p,i_squared,direction
ATP2B2,-0.78,0.12,0.0001,15.3,down
SEZ6L2,-0.65,0.14,0.0003,8.7,down
```

### `qc_report.json`

QC checkpoint results with severity levels:
```json
{
  "pipeline_ok": true,
  "checks": [
    {"check_name": "study_gene_yield", "phase": "phase1", "passed": true, "severity": "info"},
    {"check_name": "cluster_count", "phase": "phase3", "passed": true, "severity": "info"},
    ...
  ]
}
```

---

## Key Design Decisions

These are hard-won lessons from the ASD and AD proof-of-concept runs:

| Decision | Why | What Happens If Violated |
|---|---|---|
| FDR uses full seed set | Study-set-only FDR inflates significance by excluding non-significant genes from the denominator | AD stress test: 420 genes pass at q<0.10 (false). With full seed set: 0 survive (correct). |
| Welch's uses exact t-distribution | Normal approximation underestimates p-values at small sample sizes | Marginal genes get false confidence |
| log2FC range check (±10) | Raw intensities produce log2FC=119 (biologically impossible) | Entire downstream analysis is garbage |
| Blood doesn't eliminate | Brain-specific signals naturally fail in blood | Real discoveries get killed by tissue mismatch |
| Core genes locked before replication | Prevents post-hoc fishing through replication data | Undermines the entire pre-specification framework |
| Random effects is primary | Real transcriptomic data has between-study heterogeneity | Fixed effects overestimates precision |
| Expression scale check | GSE33000 uses log-ratios (max ~2.0), not log2-intensities (max ~16.0) | SE recovery produces wrong values; meta-analysis weights are wrong |

---

## HGNC Gene Symbol Resolution

The engine automatically resolves outdated and alias gene symbols using the HGNC complete set:

- **Approved symbols**: Used as-is (e.g., `BRCA1` → `BRCA1`)
- **Previous symbols**: Mapped to current (e.g., `MRE11A` → `MRE11`)
- **Aliases**: Mapped to approved (e.g., `BRCC1` → `BRCA1`)

Resolution priority: approved > previous > alias. This ensures genes aren't missed due to naming inconsistencies across GEO platforms.

Set `hgnc_path: auto` in config to download the latest HGNC data automatically, or provide a local path for offline operation.

---

## Development

### Running Tests

```bash
# Full suite
python -m pytest tests/ -q

# Just stats layer
python -m pytest tests/test_stats.py -v

# Just ingestion layer
python -m pytest tests/test_ingestion.py -v

# Just pipeline phases
python -m pytest tests/test_phases.py -v

# With coverage
python -m pytest tests/ --cov=riker --cov-report=term-missing
```

### Architecture

```
riker/
├── stats/              # Statistical primitives (locked, 96 tests)
│   ├── welch.py        # Welch's t-test with exact t-distribution
│   ├── fdr.py          # Benjamini-Hochberg with scope enforcement
│   ├── meta.py         # IVW meta-analysis (fixed + random effects)
│   └── permutation.py  # Generic permutation framework
│
├── ingestion/          # Data loading (locked, 74 tests)
│   ├── gene_db.py      # Seed genes + HGNC resolution
│   ├── normalizer.py   # Log2 detection + fold change validation
│   ├── geo_parser.py   # GEO series matrix + platform annotation
│   └── snrnaseq.py     # snRNA-seq pseudo-bulking (h5ad + CSV)
│
├── phases/             # Pipeline phases (locked, 89 tests)
│   ├── phase1_crossref.py    # Cross-referencing
│   ├── phase2_pathways.py    # Pathway mapping + feature matrix
│   ├── phase3_clustering.py  # UMAP/HDBSCAN consensus clustering
│   ├── phase4_robustness.py  # Permutation, sensitivity, LOO
│   ├── phase5_replication.py # Pre-specification + elimination
│   └── phase6_meta.py        # Effect size meta-analysis
│
├── qc/                 # Quality control
│   └── checks.py       # Phase boundary QC gates
│
├── io/                 # Input/Output
│   └── outputs.py      # CSV + JSON writers
│
├── config.py           # YAML configuration loader
└── cli.py              # Command-line interface
```

### Adding a New Condition

1. Compile a seed gene CSV from relevant databases (SFARI for ASD, AlzGene for AD, etc.)
2. Download GEO series matrix files for 3+ discovery datasets and 1+ replication datasets
3. Write a YAML config (use the quick-start example as template)
4. Run `riker validate config.yaml` then `riker run config.yaml`

The engine handles everything else — normalization, symbol resolution, phenotype detection, statistical testing, clustering, and output generation.

---

## Statistical Methods

### Meta-Analysis
Phase 6 uses inverse-variance weighted (IVW) meta-analysis with both fixed-effects and random-effects models. The random-effects model is primary and uses **REML (Restricted Maximum Likelihood)** for tau-squared estimation, with DerSimonian-Laird as fallback when REML does not converge. REML produces more accurate confidence intervals than DL when the number of studies is small (<10), which is typical for transcriptomic meta-analysis.

### Differential Expression
Phase 1 uses Welch's t-test on log2-transformed expression values. This is a conservative choice — limma or DESeq2 may provide better power for RNA-seq count data by modeling mean-variance relationships. However, Welch's on log2 data is well-established for microarray analysis and produces valid (if conservative) p-values for RNA-seq when applied to normalized data.

### FDR Correction
Benjamini-Hochberg FDR uses the **full seed gene count** as denominator, not just the study gene set. This prevents significance inflation that occurs when FDR is computed only over genes that passed an initial filter.

### Consensus Clustering
Phase 3 sweeps 15 UMAP+HDBSCAN configurations (3 n_neighbors x 5 seeds) and builds a consensus co-association matrix. Optional PCA embedding (`embedding_methods=["pca"]` or `["umap", "pca"]`) validates that discovered clusters are not artifacts of the UMAP embedding.

---

## Known Limitations

1. **Welch's t-test for RNA-seq**: For raw RNA-seq count data, variance-stabilizing methods (limma-voom, DESeq2) are generally preferred over Welch's t-test. The engine mitigates this by requiring log2 transformation and cross-dataset replication, but users working primarily with RNA-seq count data should consider this limitation.

2. **Single embedding method by default**: The default consensus clustering uses only UMAP embeddings. While the 15-configuration sweep reduces stochastic artifacts, PCA validation (opt-in via `embedding_methods=["umap", "pca"]`) provides additional confidence.

3. **No built-in pathway database download**: Pathway databases (KEGG, Reactome, MSigDB Hallmarks) must be provided manually. Automatic download is planned for a future version.

4. **Platform annotation preprocessing**: GEO platform annotation files (.annot) now handle standard header formats natively. Non-standard formats may still require manual preprocessing.

---

## Validation Results

The engine has been validated across three diseases with zero code modifications:

| Condition | Tissue | Core Genes | Replication Survival | All-Genes Recovery |
|-----------|--------|------------|---------------------|-------------------|
| **ASD** | Brain cortex | 35 | 35/35 (100%) | 27/35 (77.1%) |
| **T2D** | Pancreatic islets | 8 | 8/8 (100%) | 5/8 (62.5%) |
| **IBD** | Intestinal mucosa | 304 | 302/304 (99.3%) | 297/304 (97.7%) |

Each condition was run twice: once with curated disease-specific seed genes, and once with every expressed gene as seeds (hypothesis-free). Every curated core gene passes Phase 1 in the blind run across all three diseases.

Full results, configs, and methodology: [`results/`](results/)

---

## Roadmap

- [x] snRNA-seq pseudo-bulking module (`riker/ingestion/snrnaseq.py`) for single-nucleus data
- [ ] KEGG REST API loader for automatic pathway downloads
- [ ] Reactome bulk download parser
- [ ] Forest plot visualization (`riker/io/plots.py`)
- [ ] HTML report generator (`riker/io/report.py`)
- [ ] PyPI release

---

## Citation

If you use the Riker Engine in your research, please cite:

> Sigmon, R. (2026). Riker Engine: A condition-agnostic transcriptomics pipeline for discovering replicated gene modules. GitHub: https://github.com/RaySigmon/Riker_Engine

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

*Built on a Raspberry Pi 5 named Ghost, using Claude Code CLI as the engineering layer.*
*Named after a boy who inspired a father to search across dimensions for answers.*
