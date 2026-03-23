# Riker Engine

**A condition-agnostic transcriptomics pipeline for discovering replicated gene modules across multiple independent datasets.**

[![Tests](https://img.shields.io/badge/tests-274%20passing-brightgreen)]()
[![Python](https://img.shields.io/badge/python-3.11%2B-blue)]()
[![License](https://img.shields.io/badge/license-AGPL--3.0-blue)]()
[![Version](https://img.shields.io/badge/version-v0.2.1-orange)]()

---

## What It Does

The Riker Engine takes a list of candidate genes for any disease and determines which ones form reproducible, statistically defensible modules across multiple independent expression datasets. It enforces pre-specification, full-seed-set FDR correction, consensus clustering, and directional replication — returning a minimal core gene set that survives all filters.

## Why It Matters

Existing tools like WGCNA analyze one dataset at a time and produce large gene modules (hundreds to thousands of genes) with no built-in cross-dataset replication. Riker Engine automates multi-dataset analysis, replication filtering, and meta-analysis in a single pipeline. Switching diseases requires changing only the seed gene list and dataset configuration — zero code modifications.

---

## Validation Results

Validated across four diseases with zero code changes:

| Metric | ASD | T2D | IBD | AD |
|---|---|---|---|---|
| **Tissue** | Brain cortex | Pancreatic islets | Intestinal mucosa | Brain cortex |
| **Seed genes** | 1,267 (SFARI) | 443 (Open Targets) | 762 (Open Targets) | 801 (Open Targets) |
| **Datasets** | 7 | 4 | 6 | 5 |
| **Phase 1 yield** | 11.1% | 12.6% | 53.5% | 54.7% |
| **Core genes** | 35 | 8 | 304 | 394 |
| **Survived replication** | 35 (100%) | 8 (100%) | 302 (99.3%) | 340 (86.3%) |
| **Meta-significant** | 13 | 8 | 296 | 312 |
| **Blind recovery (Phase 1)** | 100% | 100% | 100% | 100% |
| **Blind full recovery** | 77.1% | 62.5% | 97.7% | Phase 1 only* |

*AD blind run Phase 1 completed (100% recovery); Phases 3-6 require >8GB RAM for 14,442 study genes.

Each condition was run twice: once with curated disease-specific seeds, and once with every expressed gene as seeds (hypothesis-free). Every curated core gene passes Phase 1 in the blind run across all four diseases.

---

## WGCNA Benchmark

Head-to-head comparison on the same ASD brain cortex data, same hardware (Raspberry Pi 5, 8GB):

| | Riker Engine | WGCNA |
|---|---|---|
| **Runtime** | ~8 min (3 discovery + 4 replication, all 6 phases) | ~3h 50min (3 discovery only) |
| **Core output** | 35 genes | 1,427–9,995 genes per significant module |
| **Cross-dataset validation** | Built-in (Phases 4–5) | None |
| **Memory** | Ran full gene set without issue | OOM on 18K genes; required filtering to 10K |

Full benchmark: [`benchmarks/benchmark_report.md`](benchmarks/benchmark_report.md)

---

## Quick Start

```bash
git clone https://github.com/RaySigmon/Riker_Engine.git
cd Riker_Engine
pip install -e ".[clustering]"
riker run your_config.yaml
```

See [`docs/CONFIGURATION.md`](docs/CONFIGURATION.md) for how to write a config file, or use an example from [`configs/examples/`](configs/examples/).

---

## The Six-Phase Pipeline

```
Seed Genes + GEO Datasets
        |
        v
+-------------------------+
|  Phase 1: Cross-Ref     |  Welch's t-test per gene per dataset
|  p<0.05 in >=2 datasets |  -> Study gene set
+---------+---------------+
          v
+-------------------------+
|  Phase 2: Pathways      |  KEGG + Reactome + Hallmark mapping
|  Feature matrix build   |  -> Binary pathways + expression stats
+---------+---------------+
          v
+-------------------------+
|  Phase 3: Clustering    |  UMAP + HDBSCAN x 15 configurations
|  Consensus matrix       |  -> Gene modules discovered
+---------+---------------+
          v
+-------------------------+
|  Phase 4: Robustness    |  Permutation + sensitivity + LOO
|  Core gene ID           |  -> Pre-specified replication set
+---------+---------------+
          v
+-------------------------+
|  Phase 5: Replication   |  Independent directional testing
|  Elimination protocol   |  -> Survivors + eliminated genes
+---------+---------------+
          v
+-------------------------+
|  Phase 6: Meta-Analysis |  IVW fixed + REML random effects
|  Forest plot data       |  -> Final effect sizes + I-squared
+---------+---------------+
          |
          v
    Structured Output
    (CSV + JSON + QC Report)
```

1. **Phase 1 — Cross-Referencing**: Welch's t-test (exact t-distribution) per gene per dataset. Genes significant at p < 0.05 in 2+ datasets advance.
2. **Phase 2 — Pathway Mapping**: Builds a feature matrix combining pathway memberships and expression statistics for clustering.
3. **Phase 3 — Consensus Clustering**: Sweeps 15 UMAP + HDBSCAN configurations, builds a co-association consensus matrix. Genes must cluster together across multiple parameter settings.
4. **Phase 4 — Robustness Testing**: Bonferroni-corrected permutation tests, four progressive stringency levels (including full-seed-set FDR), and leave-one-dataset-out stability. Core genes are locked before replication.
5. **Phase 5 — Replication**: Tests core genes in held-out datasets for directional concordance. Genes with significant opposite-direction effects in brain tissue are eliminated.
6. **Phase 6 — Meta-Analysis**: Inverse-variance weighted meta-analysis with REML random effects (primary) and fixed effects. Reports effect sizes, confidence intervals, and heterogeneity statistics.

Detailed methodology: [`docs/PIPELINE.md`](docs/PIPELINE.md)

---

## Installation

### Requirements

- Python 3.11 or higher
- Minimum: 8GB RAM (tested on Raspberry Pi 5)
- Recommended: 16GB+ RAM for blind genome-wide runs with >10K study genes

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
python -m pytest tests/ -q   # 274 tests
riker --help
```

### Dependencies

**Core** (installed automatically): numpy, pandas, scipy, scikit-learn, PyYAML, matplotlib

**Clustering extras** (`pip install -e ".[clustering]"`): umap-learn, hdbscan

**snRNA-seq** (optional): scanpy for AnnData/h5ad support

---

## Configuration

Switching diseases requires only changing the seed gene file and dataset list. See [`docs/CONFIGURATION.md`](docs/CONFIGURATION.md) for the full reference.

Minimal example:
```yaml
condition: ASD
seed_genes: data/asd_seed_genes.csv
hgnc_path: auto
output_dir: output/asd_run

datasets:
  - id: GSE28521
    series_matrix: data/GSE28521_series_matrix.txt.gz
    platform: data/GPL6883.annot
    role: discovery
    tissue: brain
    case_values: ["disease status: autism"]
    control_values: ["disease status: control"]

  - id: GSE64018
    series_matrix: data/GSE64018_series_matrix.txt.gz
    platform: data/GPL11154.annot
    role: replication
    tissue: brain
```

Example configs for all four validated diseases: [`configs/examples/`](configs/examples/)

---

## Output Files

```
output/
├── pipeline_summary.json      # Overall results + locked gene list
├── qc_report.json             # QC checks with pass/warn/critical
├── phase1_study_genes.csv     # Genes passing Phase 1 filter
├── phase4_core_genes.csv      # Core genes (pre-specified set)
├── phase4_all_levels.csv      # Progressive confidence pyramid
├── phase5_verdicts.csv        # Replication results per gene
└── phase6_meta_analysis.csv   # Forest plot data + effect sizes
```

---

## Key Design Decisions

| Decision | Why |
|---|---|
| FDR uses full seed set | Study-set-only FDR inflates significance by excluding non-significant genes |
| Welch's uses exact t-distribution | Normal approximation underestimates p-values at small sample sizes |
| log2FC range check (±10) | Raw intensities produce biologically impossible fold changes |
| Blood doesn't eliminate | Brain-specific signals naturally fail in blood — tissue mismatch is not evidence against |
| Core genes locked before replication | Prevents post-hoc fishing through replication data |
| Random effects is primary | Real transcriptomic data has between-study heterogeneity |

---

## Architecture

```
riker/
├── stats/              # Statistical primitives (96 tests)
│   ├── welch.py        # Welch's t-test with exact t-distribution
│   ├── fdr.py          # Benjamini-Hochberg with scope enforcement
│   ├── meta.py         # IVW meta-analysis (fixed + random effects)
│   └── permutation.py  # Generic permutation framework
├── ingestion/          # Data loading (74 tests)
│   ├── gene_db.py      # Seed genes + HGNC resolution
│   ├── normalizer.py   # Log2 detection + fold change validation
│   ├── geo_parser.py   # GEO series matrix + platform annotation
│   └── snrnaseq.py     # snRNA-seq pseudo-bulking (h5ad + CSV)
├── phases/             # Pipeline phases (89 tests)
│   ├── phase1_crossref.py    # Cross-referencing
│   ├── phase2_pathways.py    # Pathway mapping + feature matrix
│   ├── phase3_clustering.py  # Consensus clustering
│   ├── phase4_robustness.py  # Permutation, sensitivity, LOO
│   ├── phase5_replication.py # Pre-specification + elimination
│   └── phase6_meta.py        # Effect size meta-analysis
├── qc/checks.py        # Phase boundary QC gates
├── io/outputs.py       # CSV + JSON writers
├── config.py           # YAML configuration loader
└── cli.py              # Command-line interface
```

---

## Known Limitations

1. **Welch's t-test for RNA-seq**: Variance-stabilizing methods (limma-voom, DESeq2) are generally preferred for raw count data. The engine mitigates this by requiring log2 transformation and cross-dataset replication.
2. **Memory for genome-wide runs**: Blind runs with >10K study genes require >8GB RAM for the consensus clustering matrix.
3. **No built-in pathway database download**: Pathway databases must be provided manually.
4. **Single embedding default**: Default consensus clustering uses UMAP only. PCA validation is opt-in via `embedding_methods`.

---

## Citation

Paper in preparation. Please cite the GitHub repository for now:

> Sigmon, R. (2026). Riker Engine: A condition-agnostic transcriptomics pipeline for discovering replicated gene modules. GitHub: https://github.com/RaySigmon/Riker_Engine

---

## License

Licensed under [AGPL-3.0](LICENSE). Free to use for research and academic purposes. Commercial use as a hosted service requires release of source code under AGPL-3.0. See [LICENSE](LICENSE) for details.

---

## About

Built by Ray Sigmon, independent computational researcher. Developed and validated entirely on a Raspberry Pi 5 named Ghost, using Claude Code CLI as the engineering layer.

Named after his son Riker.
