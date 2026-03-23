# Changelog

## v0.3.0 — Public release preparation (2026-03-22)

### Documentation
- Rewrote README for public-facing clarity: validation table (4 diseases), WGCNA benchmark summary, streamlined quick start
- Added `docs/PIPELINE.md` — detailed six-phase methodology reference
- Added `docs/CONFIGURATION.md` — full YAML config reference with examples
- Added `CONTRIBUTING.md`
- Added AD validation report (`docs/ad_validation_report.md`)
- Example configs for all four validated diseases (`configs/examples/`)

### Benchmarks
- WGCNA head-to-head comparison on ASD brain cortex data
- Runtime, memory, specificity, and cross-dataset consistency analysis
- Full results in `benchmarks/`

### Validation
- Alzheimer's Disease validation: 801 seeds, 394 core genes, 340 survived replication, 100% blind Phase 1 recovery
- Cross-disease table now covers ASD, T2D, IBD, and AD

## v0.2.1 — Relicensed to AGPL-3.0 (2026-03-22)

- Relicensed from MIT to AGPL-3.0
- Added AGPL-3.0 header to all source files
- Copyright: 2024-2026 Ray Sigmon

## v0.2.0 — Steel-manned for publication (2026-03-21)

### Statistical Improvements
- **REML tau-squared estimator** (primary) with DerSimonian-Laird fallback for meta-analysis. More accurate confidence intervals when number of studies is small (<10).
- **Normalizer handles background-subtracted negatives**: Raw intensity data with a few negative values from background subtraction is now correctly detected and log2-transformed (clamp negatives to 0 first). Previously, any negative values caused the normalizer to skip transformation.
- **Probe mapping rate validation**: Warning emitted when <10% of probes map to genes, indicating platform annotation mismatch.
- **CLI normalizer preserves DataFrame**: Normalization now returns a DataFrame when given a DataFrame, preserving gene index and sample columns.

### Parser Robustness
- **GEO annotation headers handled natively**: Platform .annot files with ^, !, and # header lines are now parsed correctly without manual preprocessing.
- **HGNC download URL**: Updated to Google Storage mirror as primary, with EBI FTP as fallback. Clearer error message if both fail.

### Clustering Rigor
- **PCA embedding option**: `run_consensus_clustering()` accepts `embedding_methods=["pca"]` or `["umap", "pca"]` for validation that clusters are not UMAP artifacts. Default remains UMAP only.

### Documentation
- Fixed clone URL (Project-Riker -> Riker_Engine)
- Added Statistical Methods section to README
- Added Known Limitations section to README
- This CHANGELOG

### Test Count
- 280+ tests passing (262 from v0.1.1 + new tests for all changes)

## v0.1.1 — Phase 17 progressive confidence output (2026-03-20)

- Added `phase4_all_levels.csv` output (progressive confidence pyramid for all study genes)
- Wired into CLI and snRNA-seq runner
- Moved PHASE*.md build instructions to docs/build_phases/
- Added .gitignore, removed build artifacts from tracking
- 262 tests passing

## v0.1.0 — Initial release (2026-03-20)

- Complete 6-phase pipeline: cross-referencing, pathways, consensus clustering, robustness, replication, meta-analysis
- GEO series matrix parser with auto-detection of phenotype metadata
- HGNC gene symbol resolution (approved, previous, alias symbols)
- snRNA-seq pseudo-bulking from h5ad and count matrices
- QC framework with pass/warn/critical gates
- 260 tests passing
