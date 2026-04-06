# Riker Engine

A condition-agnostic transcriptomics pipeline for discovering replicated gene modules across multiple independent datasets.

[![Tests](https://github.com/RaySigmon/Riker_Engine/actions/workflows/test.yml/badge.svg)](https://github.com/RaySigmon/Riker_Engine/actions/workflows/test.yml) ![Python](https://img.shields.io/badge/python-3.11%2B-blue) ![License](https://img.shields.io/badge/license-AGPL--3.0-blue) ![Version](https://img.shields.io/badge/version-0.3.0-orange)

## What It Does

The Riker Engine takes a list of candidate genes for any disease and determines which ones form reproducible, statistically defensible modules across multiple independent expression datasets. It enforces pre-specification, full-seed-set FDR correction, consensus clustering, and directional replication — returning a minimal core gene set that survives all filters.

## Why It Matters

Existing tools like WGCNA analyze one dataset at a time and produce large gene modules (hundreds to thousands of genes) with no built-in cross-dataset replication. Riker Engine automates multi-dataset analysis, replication filtering, and meta-analysis in a single pipeline. Switching diseases requires changing only the seed gene list and dataset configuration — zero code modifications.

Across all five validated diseases, known drug targets emerge naturally without the engine being directed to seek them: ABCC8 (sulfonylureas) in T2D, JAK2 (tofacitinib) in IBD, TOP2A (doxorubicin), ERBB2 (trastuzumab), and ESR1 (tamoxifen) in breast cancer. GWAS-identified risk genes (TCF7L2, PSEN1, IL23R, PIK3CA) behave differently — genetic variants don't necessarily change transcript abundance. The engine finds functional consequences at the expression level, not genetic causes.

## Validation Results

Validated across five diseases with zero code changes:

| Metric | ASD | T2D | IBD | AD | Breast Ca. |
|---|---|---|---|---|---|
| Tissue | Brain cortex | Pancreatic islets | Intestinal mucosa | Brain cortex | Breast tumor |
| Seed genes | 1,267 (SFARI) | 443 (Open Targets) | 762 (Open Targets) | 801 (Open Targets) | 653 (Open Targets) |
| Datasets | 7 | 4 | 6 | 5 | 5 |
| Phase 1 yield | 11.1% | 12.6% | 53.5% | 54.7% | 34.8% |
| Core genes | 35 | 8 | 304 | 394 | 152 |
| Survived replication | 35 (100%) | 8 (100%) | 302 (99.3%) | 340 (86.3%) | 152 (100%) |
| Meta-significant | 13 | 8 | 296 | 312 | 121 |
| Blind recovery (Phase 1) | 100% | 100% | 100% | 100% | N/R |
| Blind full recovery | 77.1% | 62.5% | 97.7% | Phase 1 only* | N/R |

\*AD blind run Phase 1 completed (100% recovery); Phases 3–6 require >8GB RAM for 14,442 study genes. N/R = not yet run.

Each condition was run twice: once with curated disease-specific seeds, and once with every expressed gene as seeds (hypothesis-free). Every curated core gene passes Phase 1 in the blind run across all tested diseases.

### Signal Strength Calibration

The pipeline is calibrated — not biased toward any disease. Weak transcriptomic signals (ASD 11%, T2D 13%) produce small, precise modules (8–35 core genes). Moderate signals (breast cancer 35%) produce medium-sized sets (152 genes). Strong signals (IBD 54%, AD 55%) produce comprehensive modules (304–394 genes). All thresholds are identical across diseases.

### Highlights

- **T2D blind run**: From 26,800 genes with zero prior hypothesis, the engine returned IAPP (islet amyloid polypeptide) as the #1 signal — the pathological hallmark of T2D — and reconstructed the complete beta cell failure cascade (PCSK1, ERO1B, MAFB, CHGB).
- **Breast cancer**: The engine independently separated ER biology (ESR1), HER2 biology (ERBB2), and proliferation biology (TOP2A/AURKA) into distinct modules — reconstructing the current clinical subtype classification from raw expression data without being told the subtypes exist.
- **Alzheimer's**: Found TREM2, APOE, APP, MAPT, CLU, BIN1, CD33. PSEN1 correctly absent (genetic variant, not expression change).

## WGCNA Benchmark

Head-to-head comparison on the same ASD brain cortex data, same hardware (Raspberry Pi 5, 8GB):

| | Riker Engine | WGCNA |
|---|---|---|
| Runtime | ~8 min (3 discovery + 4 replication, all 6 phases) | ~3h 50min (3 discovery only, crashed once) |
| Core output | 35 genes | 1,427–9,995 genes per significant module |
| Cross-dataset validation | Built-in (Phases 4–5) | None (must be done manually) |
| Memory | Ran full 18K gene set without issue | OOM on 18K genes; required filtering to 10K |
| Scale-free fit | N/A (not required) | Failed on GSE28475 (R² = 0.35 at power 20) |
| Gene overlap | N/A | 34/35 Riker core genes appear in WGCNA modules |

WGCNA is not wrong — 34 of 35 core genes appear within its significant modules. This comparison demonstrates the workflow advantage of integrated multi-dataset analysis over manual single-dataset analysis followed by cross-study comparison. WGCNA excels at comprehensive single-dataset co-expression network analysis; the Riker Engine is designed to extract minimal, replicated gene sets across multiple independent datasets. The output size difference (35 vs. 1,427–9,995 genes) reflects different design goals, not a quality difference in the underlying methods. The engine automates what currently requires a skilled bioinformatician weeks of manual work: running WGCNA per dataset, comparing module preservation across studies, filtering for replication, and performing meta-analysis.

Full benchmark: `benchmarks/benchmark_report.md`

## Quick Start

```bash
git clone https://github.com/RaySigmon/Riker_Engine.git
cd Riker_Engine
pip install -e ".[clustering]"

# Download data for a disease (e.g., IBD — all datasets auto-download)
python scripts/download_data.py ibd

# Run the pipeline
riker run configs/examples/ibd_bulk.yaml
```

To see all available diseases and data dependencies: `python scripts/download_data.py --list`

Seed gene files are included in the repo (`data/seeds/`). Some RNA-seq datasets require manual reconstruction — see `docs/DATA_RECONSTRUCTION.md`.

See `docs/CONFIGURATION.md` for how to write a config file, or use an example from `configs/examples/`.

### Web UI

```bash
riker ui
```

Opens a browser-based interface for running and exploring pipeline results. Requires optional UI dependencies (`pip install -e ".[ui]"`).

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

1. **Phase 1 — Cross-Referencing**: Welch's t-test (exact t-distribution) per gene per dataset. Genes significant at p < 0.05 in 2+ datasets advance. Intentionally lenient — this is dimensionality reduction, not discovery.
2. **Phase 2 — Pathway Mapping**: Builds a feature matrix combining KEGG, Reactome, and MSigDB Hallmark pathway memberships with five expression statistics. Anti-circularity rule: individual pathway IDs are features; pre-assigned category labels are never used.
3. **Phase 3 — Consensus Clustering**: Sweeps 15 UMAP + HDBSCAN configurations (3 n_neighbors × 5 random seeds), builds a co-association consensus matrix. Genes must cluster together across multiple parameter settings. Converts a parameter-sensitive method into a parameter-robust one.
4. **Phase 4 — Robustness Testing**: Bonferroni-corrected permutation tests (10,000 permutations), four progressive stringency levels (including full-seed-set FDR using the entire seed set as the denominator), and leave-one-dataset-out stability. Core genes are locked before replication.
5. **Phase 5 — Replication**: Tests core genes in held-out datasets for directional concordance. Genes with significant opposite-direction effects in same-tissue datasets are eliminated. Cross-tissue non-replication (e.g., brain signal absent in blood) is tolerated — tissue mismatch is not evidence against.
6. **Phase 6 — Meta-Analysis**: Inverse-variance weighted meta-analysis with REML random effects (primary) and fixed effects. REML selected over DerSimonian-Laird because DL underestimates tau-squared with few studies (<10). Reports effect sizes, confidence intervals, and heterogeneity statistics (Cochran's Q, I², τ²).

Detailed methodology: `docs/PIPELINE.md`

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

# Web UI
pip install -e ".[ui]"
```

### Verify Installation

```bash
python -m pytest tests/ -q   # 300 tests
riker --help
```

### Dependencies

**Core** (installed automatically): numpy, pandas, scipy, scikit-learn, PyYAML, matplotlib

**Clustering extras** (`pip install -e ".[clustering]"`): umap-learn, hdbscan

**UI extras** (`pip install -e ".[ui]"`): FastAPI, uvicorn, jinja2

**snRNA-seq** (optional): scanpy for AnnData/h5ad support — enables pseudo-bulking of single-nucleus RNA-seq data for use in the standard pipeline

## Configuration

Switching diseases requires only changing the seed gene file and dataset list. See `docs/CONFIGURATION.md` for the full reference.

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

Example configs for all five validated diseases: `configs/examples/`

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

## Key Design Decisions

| Decision | Why |
|---|---|
| FDR uses full seed set | Study-set-only FDR inflates significance by excluding non-significant genes from the denominator. Stress test: study-set FDR produced 420 false positives; full-seed-set FDR correctly returned 0. |
| Welch's uses exact t-distribution | Normal approximation underestimates p-values at small sample sizes (n=8–15 per group) |
| Welch's over limma/DESeq2 | Uniform across platforms (microarray + RNA-seq), conservative, independent per gene |
| log2FC range check (±10) | Raw intensities produce biologically impossible fold changes |
| Blood doesn't eliminate | Brain-specific signals naturally fail in blood — tissue mismatch is not evidence against |
| Core genes locked before replication | Prevents post-hoc fishing through replication data (same logic as clinical trial pre-registration) |
| Random effects is primary | Real transcriptomic data has between-study heterogeneity; REML over DL for accurate CI with few studies |
| 15 consensus configurations | 3 n_neighbors × 5 seeds captures structural and stochastic sensitivity; consensus matrix makes results robust to parameter choices |

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
├── cli.py              # Command-line interface
└── ui/                 # Web UI (FastAPI + Jinja2)
```

## Known Limitations

1. **Welch's t-test for RNA-seq**: Variance-stabilizing methods (limma-voom, DESeq2) are generally preferred for raw count data. The engine mitigates this by requiring log2 transformation and cross-dataset replication. Conservative behavior (fewer false positives) is accepted as a tradeoff for cross-platform uniformity.
2. **Memory for genome-wide runs**: Blind runs with >10K study genes require >8GB RAM for the consensus clustering matrix (scales O(n²)). AD blind run (14,442 study genes) exceeded 8GB. 16GB+ recommended.
3. **Pathway databases require manual setup**: KEGG and Reactome can be downloaded programmatically, but MSigDB Hallmark gene sets require registration at [gsea-msigdb.org](https://www.gsea-msigdb.org). See `docs/DATA_RECONSTRUCTION.md`.
4. **Single embedding default**: Default consensus clustering uses UMAP only. PCA validation is opt-in via `embedding_methods` configuration.
5. **Phase 5 tissue handling**: The elimination protocol uses same-tissue vs. cross-tissue logic. Pass `discovery_tissues` in the config to specify which tissue types were used in discovery; replication datasets matching those tissues can trigger elimination, while cross-tissue non-replication is tolerated.
6. **Single developer**: Bus-factor risk. The codebase is well-tested (300 tests) and reproducible (cold-start proven), but institutional adoption would benefit from independent code review.

## Citation

Dissertation and JOSS paper in preparation. Please cite the GitHub repository:

> Sigmon, R. (2026). Riker Engine: A condition-agnostic transcriptomics pipeline for discovering replicated gene modules. GitHub: [https://github.com/RaySigmon/Riker_Engine](https://github.com/RaySigmon/Riker_Engine)

## License

Licensed under [AGPL-3.0](LICENSE). Free to use for research and academic purposes. Commercial use as a hosted service requires release of source code under AGPL-3.0. Open-core model: the engine itself is and will remain open source.

---

Built by Ray Sigmon, independent computational researcher. Developed and validated entirely on a Raspberry Pi 5 named Ghost, using Claude Code CLI as the engineering layer.

Named after his son Riker.
