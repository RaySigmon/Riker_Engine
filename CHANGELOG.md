# Changelog

## v0.3.2 — Independent validation & PI outreach (2026-04-12)

### Independent Validation
- Psoriasis validated by Gemini CLI agent: 50 core genes, 50/50 survived, 100% replication, 94.0% blind recovery
- CRC validated by Claude Code CLI agent: 264 core genes, 245/264 survived (92.8%), 97.7% blind recovery
- Both validations performed on a separate Pi 5 with zero author involvement
- Validation table now covers 8 diseases across 6 tissue types

### Blind Recovery
- Added blind recovery rates for all diseases: ASD 77.1%, T2D 62.5%, IBD 97.7%, AD 88.2% (curated scope), Breast Ca 99.3%, Psoriasis 94.0%, CRC 97.7%
- Re-ran ASD and T2D with locked engine to confirm numbers
- Re-ran AD curated with locked engine, standardized seed column

### Documentation
- Added Psoriasis/CRC to preprint and dissertation validation tables
- Aligned all CRC numbers to committed results (264/245/219)
- Updated validation docs with blind recovery and corrected numbers

## v0.3.1 — Cold-start reproducibility (2026-04-06)

### New Disease: IPF
- 6th validated disease: 354 seeds → 190 core → 170 survived → 157 meta-significant
- Cold replication in held-out GSE47460: 86.3% concordant+significant
- Stability profiling: 52 iron-clad genes across leave-one-out runs

### Cold-Start Bug Fixes (9 total)
- GPL URL bucket calculation for 3-digit platform IDs
- Seed_genes null guard — clear error when config is missing seeds
- Reproducibility variance note added to README
- IPF config and negative control config moved to examples/
- GSE47460 held-out dataset added to IPF download script
- Non-editable install as default, developer setup separated
- Gene_column and probe_column exposed in YAML dataset config
- NEW_DISEASE_GUIDE.md for configuring novel diseases
- Tool framing section ("What This Tool Does") and Phase 2 pathway docs

### Compatibility
- Pinned numpy<2 and hdbscan<0.9 for ARM/piwheels ABI compatibility
- Seed list size guidance added to NEW_DISEASE_GUIDE.md

### Web UI
- Local web UI: `riker ui` command launches browser-based interface
- Phenotype auto-detection, field sanitization, dropdown filtering
- WebSocket reconnect, config key mapping, results display fixes
- Automated 5-disease UI backend test script

### Benchmarks & Evidence
- WGCNA benchmark outputs committed (34/35 multi-dataset vs 21/35 single)
- Blind run results for breast cancer and AD committed
- CITATION.cff added
- Breast cancer validation added to README
- RunPod setup downloads filtered AD datasets from GitHub Release

### Reproducibility
- Results, data scripts, CI, Phase 5 fix overhaul
- 300 tests passing

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
