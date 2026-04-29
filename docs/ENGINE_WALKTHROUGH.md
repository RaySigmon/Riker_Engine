# Riker Engine — End-to-End Walkthrough

**Engine version:** 0.3.3 (commit 1b250ed)
**Walkthrough date:** 2026-04-29
**Walkthrough by:** Kai (Claude Code agent)
**Reviewer pending:** Cody Sigmon, Claude (chat)

---

## Pipeline overview

```
GEO Series Matrix Files + Seed Gene CSV + HGNC Database
    |
    v
[Ingestion Layer]  geo_parser.py, gene_db.py, normalizer.py
    |  Produces: gene-level expression DataFrames + phenotype assignments
    v
[Phase 1]  Cross-Referencing (phase1_crossref.py)
    |  Welch's t-test per gene per dataset. Filter: p < 0.05 in >= 2 datasets
    |  Produces: study gene set (typically 5-15% of seeds)
    v
[Phase 2]  Feature Construction (phase2_pathways.py)
    |  Builds per-gene feature matrix (expression stats + optional pathways)
    |  Produces: min-max normalized feature matrix
    v
[Phase 3]  Consensus Clustering (phase3_clustering.py)  [STOCHASTIC]
    |  15x UMAP+HDBSCAN, co-association consensus, final HDBSCAN
    |  Produces: cluster assignments per gene
    v
[Phase 4]  Robustness Testing (phase4_robustness.py)  [STOCHASTIC]
    |  Permutation significance, sensitivity analysis, LOO stability
    |  Core gene identification: Level 2 survivors in clusters with >= 3 members
    |  Produces: locked core gene list
    v
[Phase 5]  Replication (phase5_replication.py)
    |  Directional concordance in held-out datasets
    |  Tissue-aware elimination protocol
    |  Produces: survived/eliminated verdicts
    v
[Phase 6]  Meta-Analysis (phase6_meta.py)
    |  IVW fixed + REML random effects across discovery datasets
    |  Produces: pooled effect sizes, heterogeneity stats
    v
[Output]  CSV + JSON files (io/outputs.py)
```

Orchestrated by `riker/cli.py`. Configured by YAML files parsed by `riker/config.py`.

---

## Phase-by-phase walkthrough

### Phase 1 — Cross-Referencing

**File(s):** `riker/phases/phase1_crossref.py` (331 lines)
**Dependencies:** `riker/stats/welch.py`, `riker/ingestion/normalizer.py`

**Input:** Gene-level expression DataFrames (genes x samples, log2-transformed, index = gene symbols) produced by the ingestion layer. One DataFrame per discovery dataset. Phenotype assignments (dict of sample_id -> 'case' or 'control') per dataset. List of resolved seed gene symbols from SeedGeneDB.

**Computation:** For each seed gene in each discovery dataset, computes Welch's unequal-variance t-test (two-sided) between case and control expression values. The effect size is the mean difference on the log2 scale (equivalent to log2 fold change for log2-transformed input). A gene passes Phase 1 if it reaches nominal significance (p < p_threshold) in at least min_datasets discovery datasets. Two mean log2FC values are computed: `mean_log2fc` (average across all detected datasets) and `mean_log2fc_sig` (average across significant datasets only, falling back to all-datasets mean when no datasets are significant).

**Output:** `Phase1Result` containing `study_genes` dict (gene -> GeneResult with per-dataset DE statistics) and `excluded_genes` dict. Written to `phase1_study_genes.csv` with columns: gene, n_datasets_detected, n_datasets_significant, mean_log2fc, mean_log2fc_sig, consistent_direction. Consumed by Phase 2 (feature matrix construction), Phase 4 (sensitivity analysis, LOO stability, permutation test, core gene identification), Phase 5 (replication uses same Welch's t-test), and Phase 6 (DE results provide effect sizes for meta-analysis).

**Stochasticity:** Deterministic. No random operations. Welch's t-test uses `scipy.stats.t` (exact t-distribution).

**Configuration parameters:**

| Config field | Default | Description |
|---|---|---|
| `phase1.p_threshold` | 0.05 | Nominal significance cutoff per gene per dataset |
| `phase1.min_datasets` | 2 | Minimum datasets where gene must reach p_threshold |

**Structural observations:**

- `riker/phases/phase1_crossref.py:41` — `WelchResult` is imported but not used in the module body. Only `welch_ttest` (the function) is called.
- `validate_fold_changes()` from `normalizer.py` is called on a dict keyed by `gene_datasetid` (line 273), producing per-gene-per-dataset QC. The fold change threshold of 10.0 is hardcoded in `normalizer.py:29` (`FOLD_CHANGE_MAX = 10.0`), not configurable from YAML.
- `symbol_column` appears in YAML config files (e.g., `symbol_column: symbol`) but is not read by `load_config()` in `config.py`. The CLI passes no `symbol_column` to `SeedGeneDB`, which defaults to `"symbol"`. The field is consumed only by the UI server (`riker/ui/server.py:97`).

### Phase 2 — Pathway Mapping and Feature Construction

**File(s):** `riker/phases/phase2_pathways.py` (413 lines)

**Input:** `Phase1Result.study_genes` dict from Phase 1. Optional `PathwayDatabase` with KEGG/Reactome/Hallmark pathway-to-gene mappings.

**Computation:** Builds a feature matrix for downstream clustering. For each study gene, computes 5 expression-based features: `avg_log2fc` (from `mean_log2fc_sig`), `neg_log10_min_p` (capped at 10.0), `n_sig_datasets`, `direction_consistency`, and `tier_score`. If a PathwayDatabase is provided, binary pathway membership columns are added (1 if gene is in pathway, 0 otherwise). Pathways are filtered by study gene overlap (min 3 study genes, max 500 total genes, min 2% study fraction, top 100 by overlap count). All features are min-max normalized to [0, 1].

**Output:** `Phase2Result` containing `feature_matrix` DataFrame (genes x features). Written to `phase2_feature_matrix.csv`. Consumed by Phase 3 (clustering input).

**Stochasticity:** Deterministic. Sorting and normalization are deterministic operations.

**Configuration parameters:**

| Config field | Default | Description |
|---|---|---|
| `phase2.min_study_genes` | 3 | Min study genes per pathway (filter) |
| `phase2.max_total_genes` | 500 | Max total genes per pathway (filter) |
| `phase2.min_study_fraction` | 0.02 | Min study gene fraction per pathway |
| `phase2.max_pathways` | 100 | Max pathways retained |

Note: Phase 2 filter parameters are defined in `filter_pathways()` function defaults but are NOT exposed in `load_config()` or in any current YAML config. They are passed as hardcoded defaults in `run_phase2()` unless overridden programmatically.

**Structural observations:**

- `riker/phases/phase2_pathways.py:34` — `warnings` is imported but not used.
- In all current disease configs, no pathway database is configured (`config.pathways` is empty). Phase 2 runs with 5 expression features only (0 pathway features). The entire PathwayDatabase infrastructure (classes, filtering logic, ~150 lines) is defined but not exercised by any current pipeline invocation.
- `gene_tiers` parameter in `build_feature_matrix()` (line 255) is never passed by `run_phase2()` — it always receives `None`. The `tier_score` feature column is always 0.0 for all genes in current usage.
- Phase 2 pathway filter parameters (`min_study_genes`, `max_total_genes`, etc.) have no corresponding fields in `PipelineConfig` dataclass or `load_config()`. A user cannot configure them via YAML without code changes.

### Phase 3 — Consensus Clustering

**File(s):** `riker/phases/phase3_clustering.py` (507 lines)
**Dependencies:** `umap-learn`, `hdbscan` (or `sklearn.cluster.HDBSCAN`)

**Input:** `Phase2Result.feature_matrix` DataFrame (genes x features, min-max normalized to [0,1]).

**Computation:** Runs UMAP dimensionality reduction to 2D followed by HDBSCAN density-based clustering across multiple parameter configurations (default: 3 n_neighbors values x 5 random seeds = 15 configurations). For each gene pair, computes the fraction of configurations in which both genes were assigned to the same cluster (excluding noise assignments where either gene has label -1). This produces an n x n co-association consensus matrix. Final cluster assignments are derived by converting the consensus matrix to a distance matrix (1 - consensus) and applying HDBSCAN with `metric='precomputed'`.

**Output:** `Phase3Result` containing `cluster_labels` dict (gene -> cluster_id, -1 for noise), `consensus_matrix` DataFrame, `cluster_info` dict, and `per_config_labels`. Written to `phase3_cluster_assignments.csv` with columns: gene, cluster_id. Consumed by Phase 4 (cluster significance, sensitivity analysis, LOO stability, core gene identification).

**Stochasticity:** Stochastic. UMAP `random_state` parameter varies across configurations. Configured via `phase3.seeds` in YAML (default: [42, 123, 456, 789, 1024]).

**Configuration parameters:**

| Config field | Default | Description |
|---|---|---|
| `phase3.n_neighbors` | [10, 15, 30] | UMAP neighbor counts to sweep |
| `phase3.seeds` | [42, 123, 456, 789, 1024] | Random seeds for UMAP |
| `phase3.min_cluster_size` | 5 | HDBSCAN min_cluster_size |
| `phase3.min_samples` | 3 | HDBSCAN min_samples |

**Structural observations:**

- `embedding_methods` parameter (line 289) accepts `["umap"]`, `["pca"]`, or `["umap", "pca"]` but is not exposed in `load_config()` or YAML config. Only reachable via programmatic API. PCA path (lines 377-410) is functional but not testable via CLI.
- `final_min_cluster_size` and `final_min_samples` parameters (lines 291-292) allow overriding HDBSCAN parameters for the final consensus clustering step. Not exposed in config. Default behavior: same parameters for per-config and final clustering.
- `build_consensus_matrix()` (line 181) uses triple-nested Python loops: outer over configurations, middle and inner over gene pairs. Complexity is O(n^2 * n_configs) where n = number of study genes. For 1,800 study genes x 15 configs, this is the dominant computation cost in Phase 3.
- HDBSCAN import (lines 42-52) tries `sklearn.cluster.HDBSCAN` first, then falls back to standalone `hdbscan` package. The `_HDBSCAN_SOURCE` variable tracks which was loaded.
- `n_neighbors` is clamped to `min(nn, n_genes - 1)` with a floor of 2 (lines 335-337). This prevents UMAP errors when the study gene count is very small.

### Phase 4 — Robustness Testing

**File(s):** `riker/phases/phase4_robustness.py` (558 lines)
**Dependencies:** `riker/stats/permutation.py`, `riker/stats/fdr.py`

**Input:** `Phase1Result.study_genes` (per-gene DE statistics), `Phase3Result.cluster_info` (cluster membership), seed gene count (total seeds for FDR scope), discovery dataset IDs (for LOO analysis).

**Computation:** Four independent assessments:

1. **Permutation significance** (`evaluate_cluster_significance`): For each cluster, tests whether the cluster's mean |log2FC| exceeds random gene groups of the same size. 10,000 permutations, Bonferroni correction across clusters. Significance threshold: corrected p < 0.05. Reported but does NOT gate core gene selection.
2. **Sensitivity analysis** (`sensitivity_analysis`): Tests each gene at 4 progressive thresholds — Level 1 (p<0.05 in >=2 datasets), Level 2 (p<0.01 in >=2), Level 3 (FDR q<0.10 in >=2), Level 4 (FDR q<0.05 in >=2). FDR uses full-seed-set scope enforcement (pads with p=1.0 up to seed count).
3. **LOO stability** (`loo_stability`): Removes each discovery dataset in turn, checks if >=80% of cluster genes still meet Phase 1 criterion with remaining datasets.
4. **Core gene identification** (`identify_core_genes`): Core genes = Level 2 survivors (p<0.01 in >=2 datasets) from clusters with >=3 Level 2 survivors. This list is locked before Phase 5.

**Output:** `Phase4Result` containing cluster significance, sensitivity results, LOO stability, and `core_genes` dict (gene -> CoreGene). Written to `phase4_core_genes.csv` (core genes) and `phase4_all_levels.csv` (all study genes with level classifications). Consumed by Phase 5 (core gene list) and Phase 6 (cluster assignments for reporting).

**Stochasticity:** Stochastic. Permutation test uses `np.random.default_rng(seed)` where seed = `phase4.seed` + cluster_id. Different clusters get different seeds.

**Configuration parameters:**

| Config field | Default | Description |
|---|---|---|
| `phase4.n_permutations` | 10000 | Permutations per cluster |
| `phase4.seed` | 42 | Base random seed for permutation test |

**Structural observations:**

- `riker/phases/phase4_robustness.py:33` — `warnings` is imported but not used.
- `riker/phases/phase4_robustness.py:38` — `PermutationResult` is imported but not referenced in the module body (the return type is used implicitly via `cluster_permutation_test`).
- `riker/phases/phase4_robustness.py:39` — `fdr_survivors` and `FDRResult` are imported but not used in the module body.
- `riker/phases/phase4_robustness.py:288` — FDR levels 3-4 check `q_val < threshold` but also require `n_sig >= min_datasets` using hardcoded `p < 0.05` for dataset counting (line 288). The hardcoded 0.05 is not the level's FDR threshold — it's a raw p-value gate for the "in how many datasets" criterion. This mixing of FDR and raw-p thresholds is documented in the code comment but not in CONFIGURATION.md or PIPELINE.md.
- `SENSITIVITY_LEVELS` dict (lines 46-49) defines the four thresholds as module-level constants, not configurable from YAML.
- The `min_genes_per_cluster` parameter in `identify_core_genes()` (line 420, default 3) is not configurable from YAML. It controls the anti-singleton guard threshold.
- The LOO stability threshold of 0.80 (line 338) and the Bonferroni significance threshold of 0.05 (line 204) are function parameter defaults, not YAML-configurable.

### Phase 5 — Independent Replication

**File(s):** `riker/phases/phase5_replication.py` (487 lines)
**Dependencies:** `riker/stats/welch.py`

**Input:** `Phase4Result.core_genes` dict (locked, immutable core gene list), replication dataset expression DataFrames, phenotype assignments, tissue types per dataset, discovery tissue types (for cross-tissue tolerance).

**Computation:** Tests each core gene in each replication dataset using Welch's t-test (same implementation as Phase 1). Compares the direction of effect to the discovery direction. Elimination rule: a gene is eliminated if it shows significant (p < 0.05) opposite-direction expression in a same-tissue replication dataset. Cross-tissue replication failures (e.g., brain discovery signal absent in blood) do not trigger elimination. When `discovery_tissues` is provided, replication datasets matching those tissues are treated as same-tissue; others are cross-tissue. Per-cluster verdicts are assigned: `replicated`, `tissue_specific`, `partially_replicated`, or `failed`.

**Output:** `Phase5Result` containing `gene_verdicts` dict (gene -> GeneVerdict with status, reason, per-dataset results), `cluster_verdicts` dict, and `locked_core_genes` list. Written to `phase5_verdicts.csv`. Consumed by Phase 6 (determines which genes enter meta-analysis).

**Stochasticity:** Deterministic. Uses the same Welch's t-test as Phase 1.

**Configuration parameters:**

Replication datasets and tissue types are configured per-dataset in the YAML `datasets` section (`role: replication`, `tissue: brain`). The elimination p-value threshold of 0.05 is hardcoded (line 225), not configurable from YAML. This is by design per the METHODS doc.

**Structural observations:**

- `riker/phases/phase5_replication.py:35` — `warnings` is imported but not used.
- `replicate_gene()` function (line 146) has `tissue: str = "brain"` as a default parameter. This default is never reached in practice because `run_elimination_protocol()` always passes the tissue from `dataset_tissues` dict.
- The elimination threshold `p_threshold: float = 0.05` (line 225) is a function parameter default in `run_elimination_protocol()`. It is not configurable from YAML and the CLI does not pass it — the default is always used.
- When no replication datasets are configured, the CLI (cli.py lines 320-339) creates a synthetic `Phase5Result` with all core genes marked as `survived` with reason "No replication datasets configured; retained by default."

### Phase 6 — Effect Size Meta-Analysis

**File(s):** `riker/phases/phase6_meta.py` (319 lines)
**Dependencies:** `riker/stats/meta.py`

**Input:** `Phase1Result.study_genes` (per-dataset DE statistics for effect sizes and SEs), `Phase5Result.gene_verdicts` (determines which genes are analyzed — only `survived` genes), optional `dataset_tissues` dict and `dataset_expression_ranges` for metadata and scale checking.

**Computation:** For each gene surviving Phase 5, collects per-dataset effect sizes (log2FC) and standard errors from Phase 1 DE results. If SE is invalid (zero or non-finite), attempts recovery from log2FC, p-value, and sample sizes using the inverse t-distribution. Runs both fixed-effects and random-effects (REML with DL fallback) inverse-variance weighted meta-analysis via `run_meta_analysis()`. Computes heterogeneity statistics (Cochran's Q, I-squared, tau-squared). Reports pooled effect sizes, confidence intervals (95%, z-based with 1.96 multiplier), and two-sided p-values.

**Output:** `Phase6Result` containing `gene_results` dict (gene -> GeneMetaResult with per-dataset forest plot data and pooled estimates). Written to `phase6_meta_analysis.csv` with columns: gene, cluster_id, random_effect, random_se, random_p, fixed_effect, fixed_se, fixed_p, i_squared, tau_squared, cochran_q, n_datasets, direction. Consumed by disease reports and downstream analysis.

**Stochasticity:** Deterministic. All meta-analysis computations are closed-form (IVW) or iterative with deterministic starting points (REML Fisher scoring). Verified empirically: sampled 3 genes across 5 stability runs, all Phase 6 values identical to 17-18 significant digits.

**Configuration parameters:** Phase 6 has no user-configurable parameters. The random-effects model is always primary. REML is always attempted before DL fallback.

**Structural observations:**

- `riker/phases/phase6_meta.py:35` — `pd` (pandas) is imported but not used in the module body.
- `riker/stats/meta.py:302-303, 413-414` — Confidence intervals use hardcoded 1.96 (z-score for 95% CI). Not configurable. This is standard practice for meta-analysis but differs from the t-distribution used in Phase 1 for individual tests.
- `check_expression_scale()` (called from `run_phase6` line 265) receives `np.array([max_val])` — a single-element array — rather than the full expression distribution. The scale check flags datasets with max expression below 5.0 as potential log-ratio format. The single-value check is a proxy, not a full distribution analysis.
- `compute_gene_meta()` iterates over `de_results` from Phase 1, which includes all discovery datasets where the gene was detected — not just significant ones. The effect sizes entering meta-analysis are from all detected datasets, not filtered by significance. This is standard meta-analysis practice (include all studies, weight by precision).

---

## Wrapper layers

### CLI (riker/cli.py)

**File:** `riker/cli.py` (486 lines)

**Function:** Entry point for the pipeline. Provides three subcommands: `run` (full pipeline), `validate` (config check only), `ui` (web interface). No subcommand defaults to `ui`. The `run` command orchestrates all six phases sequentially, loading data, running each phase, checking QC gates between phases, and writing outputs.

**Structural observations:**

- Imports are done inline within `cmd_run()` (e.g., `from riker.phases.phase1_crossref import run_phase1` inside the Phase 1 try block). This is consistent across all phases — a deliberate pattern that defers import cost until needed.
- `setup_logging()` is called only for `run` and `validate` commands, not for `ui`. UI logging is handled by uvicorn.
- The `--host` and `--port` arguments are defined both at the top-level parser (lines 439-442) and on the `ui` subparser (lines 458-461). The top-level versions are used when no subcommand is given (default UI mode).
- Phase 5 skip path (lines 320-339): when no replication datasets exist, creates a synthetic Phase5Result with all core genes marked as survived. This was a v0.3.3 fix (Round 2, #3).

### Configuration loading (riker/config.py)

**File:** `riker/config.py` (247 lines)

**Function:** Parses YAML config files into `PipelineConfig` and `DatasetConfig` dataclasses. Resolves relative paths against the config file's parent directory. Derives phase-specific seeds from `random_seed` when they aren't explicitly set.

**Structural observations:**

- `symbol_column` is present in many YAML config files but is not read by `load_config()`. The field is consumed only by the UI server. CLI-driven runs always use the default `"symbol"` column name.
- `_resolve_path()` handles the sentinel value `"auto"` for `hgnc_path`. No other sentinel values exist across current configs.
- `DatasetConfig.tissue` defaults to `"brain"` (line 59). For non-brain diseases (T2D islets, IBD colon, etc.), configs must explicitly set `tissue:` to the correct value. If omitted, the tissue is silently classified as brain.
- Phase 2 filter parameters (`min_study_genes`, `max_total_genes`, `min_study_fraction`, `max_pathways`) are not represented in `PipelineConfig` and cannot be configured via YAML.
- `embedding_methods` for Phase 3 is not represented in `PipelineConfig`.

### IO and outputs (riker/io/outputs.py)

**File:** `riker/io/outputs.py` (305 lines)

**Function:** Writes pipeline results to structured CSV and JSON files. Functions: `write_phase1_summary`, `write_phase2_feature_matrix`, `write_phase3_cluster_assignments`, `write_phase4_core_genes`, `write_phase4_all_levels`, `write_phase5_verdicts`, `write_phase6_meta`, `write_qc_report`, `write_pipeline_summary`. Also provides `_get_git_commit_hash()` for provenance.

**Structural observations:**

- `riker/io/outputs.py:30` — `asdict` from dataclasses is imported but not used.
- `riker/io/plots.py` and `riker/io/report.py` exist as files but contain only `pass` statements (19 lines each, all license header). Placeholder modules with no functionality.
- Output column names are not versioned or documented in a schema file. Column names are defined inline in each write function. Adding a column in a future version would not be detectable by consumers without reading the code or diff.

### Stability profiler (scripts/stability_profiling.py)

**File:** `scripts/stability_profiling.py` (530 lines)

**Function:** Runs the full pipeline N times (default 50) with deterministically varied seeds derived from a master seed. Each run gets unique UMAP seeds (Phase 3) and a unique permutation seed (Phase 4, derived via v0.3.3 random_seed wiring). Classifies genes by appearance frequency: iron-clad (>=90%), borderline (50-89%), stochastic (<50%). Computes pairwise Jaccard similarity between all run pairs. Writes stability_summary.json, stability_scores.csv, stability_pairwise_jaccard.csv, run_summary.csv, and stability_report.txt.

**Structural observations:**

- `modify_config_for_run()` (line 38) reads the source YAML directly via `yaml.safe_load()`, not via `load_config()`. Config validation and path resolution are handled separately — paths are resolved against the source config's directory using `_resolve_path()` before writing the temp config to `/tmp/`. This was a v0.3.3.1 hotfix after the path resolver in Round 8 broke the profiler's temp config pattern.
- `scripts/stability_profiling.py:61` — TODO comment: "Refactor profiler to use load_config() + serialize back, so config validation and defaults apply consistently." Added in v0.3.3.
- The profiler does not aggregate Phase 6 effect sizes or p-values across runs. The separate script `scripts/aggregate_phase6_iron_clad.py` was written for this purpose. The aggregation confirmed that Phase 6 values are identical across runs (deterministic given Phase 5 survivors).
- `--seed-list` flag (line 283) adds cross-reference against a seed gene CSV but this feature has no test coverage.

---

## Cross-cutting concerns

### Random seed handling

Seeds flow from config through two paths:

1. **Direct config:** `phase3.seeds` (list of 5 UMAP random states) and `phase4.seed` (permutation base seed) are read from YAML and passed to their respective phases.
2. **Derivation:** If `random_seed` is explicitly set in the YAML but phase-specific seeds are not, `load_config()` derives them: `phase3_seeds` = 5 values from `Random(random_seed)`, `phase4_permutation_seed` = 1 value from `Random(random_seed + 1000)`. If `random_seed` is absent from the YAML, historical defaults are preserved ([42, 123, 456, 789, 1024] and 42).
3. **Stability profiler:** `modify_config_for_run()` sets `random_seed` per run and generates explicit `phase3.seeds`. It does NOT set `phase4.seed`, so the phase4 seed is derived from the per-run `random_seed` via the derivation path. This produces unique permutation seeds per run — a v0.3.3 improvement over v0.3.2 where all runs used the same phase4 seed.

### Logging and progress reporting

All phases log to the `riker` logger hierarchy via Python's `logging` module. CLI configures logging at INFO level (DEBUG with `-v` flag). Each phase prints a `====` banner at start. Per-gene and per-dataset progress is logged at INFO level. Warnings use `warnings.warn()` for user-facing issues (e.g., NaN values, fold change range violations) and `logger.warning()` for operational issues (e.g., low probe mapping rate). The profiler writes to a separate log file via stdout redirect.

### Error handling and QC gates

Four QC checkpoints in the CLI, one after each major phase boundary:

1. **Phase 1:** Study gene yield check. Critical if <1% of seeds pass, warning if <5%.
2. **Phase 3:** Cluster count check. Critical if 0 clusters, warning if <2.
3. **Phase 4:** Core gene count check. Critical if 0 core genes, warning if <3.
4. **Phase 5:** Replication survival check. Critical if 0 survived.

Critical failures halt the pipeline and write a partial QC report. Warnings are logged but do not halt. Phase 2 and Phase 6 have no QC gates. Each phase wraps its work in a try/except that catches all exceptions, logs the traceback, and returns exit code 1.

### Output schema versioning

Output file schemas are implicit in the code. Column names are defined inline in each `write_*` function in `riker/io/outputs.py`. There is no schema definition file, no version field in CSV headers, and no mechanism for consumers to detect schema changes between engine versions. `pipeline_summary.json` includes `package_version` and `code_version` fields (added in v0.3.3), which can serve as indirect schema version indicators.

### Caching and intermediate state

The engine does not cache intermediate results between phases. Each phase receives its input as an in-memory Python object from the previous phase's return value. There is no checkpointing mechanism — if the engine crashes at Phase 4, Phases 1-3 must be re-run. Phase outputs are written to disk as CSV/JSON for audit purposes but are not read back during the same pipeline invocation. The HGNC database is cached at `~/.riker/hgnc_complete_set.txt` after first download.

---

## Engine-level structural observations

- `riker/io/plots.py` (19 lines) and `riker/io/report.py` (19 lines) are placeholder modules containing only license headers and `pass`. No functionality.
- All 5 `__init__.py` files in subpackages (`phases/`, `stats/`, `ingestion/`, `qc/`, `io/`) contain only the AGPL license header. No re-exports, no package-level imports.
- `riker/__init__.py` contains only `__version__ = "0.3.3"`. No public API surface defined at package level.
- Three phases independently compute mean log2FC from `de_results` before the v0.3.3 refactor centralized it on `GeneResult.mean_log2fc_sig`. Post-v0.3.3, Phase 2 reads from GeneResult, Phase 4 reads from GeneResult, but Phase 4 `evaluate_cluster_significance()` also reads from GeneResult. The centralization is complete.
- The 95% CI multiplier 1.96 in `riker/stats/meta.py` (lines 302-303, 413-414) is a z-score, not a t-distribution quantile. Standard for meta-analysis but differs from the t-distribution-based Welch's test in Phase 1. The CI is z-based because the pooled estimate is treated as normally distributed (valid when combining multiple studies).
- `DatasetConfig.tissue` defaults to `"brain"` across the engine. For multi-tissue or non-brain diseases, every dataset must explicitly set `tissue:` in the YAML. There is no validation that tissue values are consistent within discovery or replication groups.
- The engine has no mechanism for resuming a partial pipeline run. If Phase 4 fails, all prior phases must be re-executed. The ~8 minute wall clock for a full blind run on Pi 5 makes this a minor concern in practice.
- The stability profiler's `modify_config_for_run()` reads YAML directly via `yaml.safe_load()` and resolves paths manually, bypassing `load_config()` validation. Config validation (required fields, dataset role checks) is not applied to per-run temp configs.
- `riker/ingestion/gene_db.py` — `csv` (line 30) and `os` (line 32) are imported but not used in the module body.
- `riker/phases/phase6_meta.py` — `warnings` (line 31) and `MetaResult` (line 41, imported from stats/meta.py) are imported but not used. The module defines its own `GeneMetaResult` dataclass and uses `run_meta_analysis()` which returns `MetaResult` objects, but the type name is never referenced.
- `riker/ingestion/snrnaseq.py` — `pseudo_bulk_from_h5ad()` (line 303) is defined but never called from anywhere in the engine, tests, or scripts. It is the only function in the engine with zero call sites. It requires scanpy (optional dependency) and is designed for h5ad file input.
- The `pathways` top-level config field is parsed by `load_config()` (line 239) and consumed by `cli.py` (lines 186-194), but is not documented in `docs/CONFIGURATION.md`.
- `PipelineConfig.random_seed` (line 113) is used only during config loading to derive `phase3_seeds` and `phase4_permutation_seed`. After `load_config()` returns, the `random_seed` attribute is never read by any phase, cli.py, or script. It is a derivation intermediate, not a consumed runtime parameter.
- Unused imports total: 21 instances across 11 files (gene_db.py: 3, geo_parser.py: 1, snrnaseq.py: 2, outputs.py: 1, phase1: 1, phase2: 1, phase4: 4, phase5: 1, phase6: 3, fdr.py: 1, ui/runner.py: 1).

---

## Verified determinism profile

| Phase | Deterministic? | Stochastic element |
|---|---|---|
| Phase 1 (differential expression) | Yes | None |
| Phase 2 (feature matrix) | Yes | None |
| Phase 3 (consensus clustering) | **No** | UMAP random_state |
| Phase 4 (robustness testing) | **No** | Permutation seed |
| Phase 5 (replication) | Yes | None |
| Phase 6 (meta-analysis) | Yes | None |

Verified empirically on 2026-04-29: Phase 6 values (random_effect, random_se, random_p, fixed_effect, fixed_p, tau_squared) identical to 17-18 significant digits across sampled stability runs for BAG3, CCN1, SFMBT2. Code inspection of phase1_crossref.py, phase5_replication.py, phase6_meta.py, and stats/meta.py confirmed zero stochastic operations.

---

*Walkthrough generated by Kai (Claude Code agent), 2026-04-29. Pending review by Cody Sigmon and Claude (chat).*