# Pipeline Reference

Detailed description of the Riker Engine's six-phase progressive filtering pipeline for discovering replicated gene modules from public transcriptomic data.

---

## Overview

The Riker Engine processes candidate gene lists through six sequential phases, each applying progressively stricter statistical filters. Genes that fail at any phase do not advance. The pipeline enforces a strict separation between discovery (Phases 1--4) and replication (Phase 5), with a pre-specification protocol that locks the core gene list before replication data is accessed.

```
Seed genes + Discovery datasets
  --> Phase 1: Cross-Referencing      (differential expression filter)
  --> Phase 2: Pathway Mapping        (feature matrix construction)
  --> Phase 3: Consensus Clustering   (gene module discovery)
  --> Phase 4: Robustness Testing     (permutation, sensitivity, LOO)
  --> [Core gene list locked]
  --> Phase 5: Replication            (independent directional testing)
  --> Phase 6: Meta-Analysis          (pooled effect sizes)
  --> Structured output (CSV + JSON + QC report)
```

Source: `riker/phases/phase1_crossref.py` through `riker/phases/phase6_meta.py`.

---

## Phase 1: Cross-Referencing

### What it does

For each seed gene in each discovery dataset, computes a case-versus-control differential expression test. Genes reaching nominal significance in a minimum number of datasets advance to the study gene set. This threshold is intentionally lenient; stricter thresholds are applied in Phase 4.

### Statistical method

**Welch's t-test** (unequal variances) on log2-transformed expression values. P-values are computed from the exact t-distribution with Welch--Satterthwaite degrees of freedom, not the normal approximation. The test is two-sided.

The effect size is the mean difference on the log2 scale (i.e., log2 fold change), computed as `mean(cases) - mean(controls)` on log2-transformed data.

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `p_threshold` | `0.05` | Nominal significance cutoff per gene per dataset |
| `min_datasets` | `2` | Minimum number of datasets in which a gene must reach `p_threshold` |

### Inputs

- **Seed gene list**: Resolved HGNC gene symbols (from `seed_genes` CSV).
- **Discovery datasets**: Dict of `dataset_id -> DataFrame` (genes x samples, log2-transformed expression, index = gene symbols).
- **Phenotype assignments**: Dict of `dataset_id -> {sample_id: 'case' or 'control'}`.

### Outputs

`Phase1Result` containing:

- `study_genes`: Genes passing the filter, with per-dataset DE statistics (log2FC, p-value, t-statistic, degrees of freedom, standard error, sample sizes, direction).
- `excluded_genes`: Genes failing the filter, with the same statistics.
- `per_dataset_coverage`: Number of seed genes detected per dataset.
- `qc_warnings`: Fold change range validation results. Flags gene-dataset pairs with |log2FC| > 10, which typically indicates raw intensity data was not log2-transformed.

### QC gate

The normalizer auto-detects expression scale and applies `log2(x + 1)` when median expression >= 20. The fold change validator (`validate_fold_changes`) catches impossible values (e.g., log2FC = 119) that would corrupt all downstream analysis.

### Configurable options

```yaml
phase1:
  p_threshold: 0.05   # Any value in (0, 1)
  min_datasets: 2      # Integer >= 1; must not exceed number of discovery datasets
```

---

## Phase 2: Pathway Mapping

### What it does

Maps study genes from Phase 1 to biological pathways and constructs a mixed feature matrix combining binary pathway memberships and continuous expression statistics. This matrix is the input to unsupervised clustering in Phase 3.

### Feature construction

Each gene receives a feature vector with two groups of features:

**Binary pathway features** (one column per retained pathway):
- Value is 1 if the gene belongs to that pathway, 0 otherwise.
- Pathway sources: KEGG, Reactome, MSigDB Hallmark gene sets.

**Expression statistic features** (5 columns):

| Feature | Computation |
|---|---|
| `avg_log2fc` | Mean log2FC across all discovery datasets where the gene was detected |
| `neg_log10_min_p` | `-log10(min p-value)` across datasets, capped at 10.0 |
| `n_sig_datasets` | Number of datasets where the gene was significant at `p_threshold` |
| `direction_consistency` | Fraction of datasets agreeing with the majority direction |
| `tier_score` | External confidence tier score (e.g., SFARI tier), or 0.0 if not provided |

All features are **min-max normalized to [0, 1]** across genes. Constant columns are set to 0.0 (or 1.0 if all values are nonzero).

### Anti-circularity rule

Individual pathway IDs (e.g., `hsa04020`, `R-HSA-123456`) are used as binary features. Pre-assigned biological category labels (e.g., "chromatin remodeling", "synaptic signaling") are never used as clustering features. Using category labels would create tautological enrichment where genes cluster by their assigned labels rather than by data-driven similarity.

### Pathway filtering

Before feature matrix construction, pathways are filtered to remove uninformative or overly broad pathways:

| Parameter | Default | Description |
|---|---|---|
| `min_study_genes` | `3` | Minimum study genes that must belong to the pathway |
| `max_total_genes` | `500` | Maximum total genes in the pathway (removes overly broad pathways) |
| `min_study_fraction` | `0.02` | Minimum fraction of pathway genes that are study genes |
| `max_pathways` | `100` | Maximum pathways retained (top by study gene count) |

### Inputs

- `Phase1Result`: Study genes with per-dataset DE statistics.
- `PathwayDatabase`: Combined pathway-to-gene mappings from KEGG, Reactome, and/or Hallmark.
- `gene_tiers` (optional): Confidence tier scores per gene.

### Outputs

`Phase2Result` containing:

- `feature_matrix`: DataFrame (genes x features), min-max normalized to [0, 1].
- `pathway_features`: List of pathway IDs used as features.
- `pathway_info`: Metadata per pathway (ID, name, source, total genes, study gene overlap).
- `expression_feature_names`: Names of the 5 expression statistic features.

### Configurable options

```yaml
phase2:
  min_study_genes: 3
  max_total_genes: 500
  min_study_fraction: 0.02
  max_pathways: 100
```

---

## Phase 3: Consensus Clustering

### What it does

Discovers gene modules by running UMAP dimensionality reduction followed by HDBSCAN density-based clustering across multiple parameter configurations, then derives stable consensus clusters from a co-association matrix. No biological hypotheses are imposed during clustering.

### Statistical method

1. **Per-configuration embedding and clustering**: For each combination of `n_neighbors` and random `seed`, UMAP reduces the feature matrix to 2 dimensions, then HDBSCAN identifies density-based clusters. Genes not assigned to any cluster are labeled as noise (-1).

2. **Co-association consensus matrix**: For each pair of genes (i, j), the consensus score is the fraction of configurations in which both genes were assigned to the same cluster, excluding configurations where either gene was noise. This produces an n x n symmetric matrix with values in [0, 1].

3. **Final consensus clustering**: The consensus matrix is converted to a distance matrix (`1 - consensus`) and HDBSCAN is applied with `metric='precomputed'` to derive final cluster assignments.

### UMAP parameters (per configuration)

| Parameter | Default | Description |
|---|---|---|
| `n_components` | `2` | Fixed. Embedding dimensionality |
| `min_dist` | `0.0` | Fixed. Allows tightly packed clusters |
| `metric` | `euclidean` | Fixed. Distance metric on feature space |
| `n_neighbors` | Swept: `[10, 15, 30]` | Controls local vs global structure |
| `random_state` | Swept: `[42, 123, 456, 789, 1024]` | Controls stochastic initialization |

### HDBSCAN parameters

| Parameter | Default | Description |
|---|---|---|
| `min_cluster_size` | `5` | Minimum genes to form a cluster |
| `min_samples` | `3` | Controls cluster density threshold |
| `final_min_cluster_size` | Same as `min_cluster_size` | Override for final consensus clustering step |
| `final_min_samples` | Same as `min_samples` | Override for final consensus clustering step |

### Configuration sweep

Total configurations = `len(n_neighbors_list) x len(seeds)`. Default: 3 x 5 = **15 configurations**.

This ensemble approach makes results robust to individual UMAP/HDBSCAN parameter choices. Genes that only cluster together under one specific configuration will not appear in the consensus clusters.

### Optional PCA validation

Set `embedding_methods: ["umap", "pca"]` to add PCA-based configurations alongside UMAP. PCA embeddings are deterministic (no seed sweep), so this adds `len(n_neighbors_list)` additional configurations. This validates that discovered clusters are not artifacts of the UMAP embedding.

### Inputs

- `Phase2Result.feature_matrix`: DataFrame (genes x features), min-max normalized.

### Outputs

`Phase3Result` containing:

- `cluster_labels`: Dict of `gene_symbol -> cluster_id` (-1 = noise).
- `consensus_matrix`: DataFrame (genes x genes), co-clustering frequencies.
- `cluster_info`: Per-cluster metadata (gene list, size, mean within-cluster consensus score).
- `per_config_labels`: Full label assignments from each configuration (for diagnostics).
- `n_clusters`, `n_noise`: Summary counts.

### Configurable options

```yaml
phase3:
  n_neighbors: [10, 15, 30]
  seeds: [42, 123, 456, 789, 1024]
  min_cluster_size: 5
  min_samples: 3
  embedding_methods: ["umap"]   # or ["umap", "pca"] for PCA validation
```

---

## Phase 4: Robustness Testing

### What it does

Applies four progressive statistical filters to clusters from Phase 3. Clusters and genes that survive all four filters constitute the pre-specified core gene set for replication.

### Filter 1: Permutation significance

Tests whether each cluster's mean |log2FC| is more extreme than expected by chance.

- **Test statistic**: Mean absolute log2FC of genes in the cluster.
- **Null distribution**: Generated by randomly sampling gene sets of the same size from all study genes, `n_permutations` times per cluster.
- **P-value**: Fraction of permuted statistics >= observed statistic.
- **Correction**: Bonferroni correction across all clusters (`bonf_p = min(raw_p * n_clusters, 1.0)`).
- **Significance**: Bonferroni-corrected p < 0.05.

| Parameter | Default | Description |
|---|---|---|
| `n_permutations` | `10000` | Permutations per cluster |
| `seed` | `42` | Base random seed (incremented by cluster_id for independence) |

### Filter 2: Sensitivity analysis

Tests genes at four progressively stricter thresholds to identify those robust to threshold choice:

| Level | Criterion | FDR scope |
|---|---|---|
| 1 | p < 0.05 in >= 2 datasets | N/A (raw p-values) |
| 2 | p < 0.01 in >= 2 datasets | N/A (raw p-values) |
| 3 | FDR q < 0.10 in >= 2 datasets | **Full seed gene set** |
| 4 | FDR q < 0.05 in >= 2 datasets | **Full seed gene set** |

**Full-seed-set FDR**: At Levels 3 and 4, Benjamini-Hochberg correction uses the total number of seed genes as the denominator, not just the study genes that passed Phase 1. This prevents significance inflation that occurs when FDR is computed only over genes that already passed an initial filter. The function `apply_fdr_with_scope()` enforces this.

The **dissolution point** is the first level at which fewer than 3 genes in a cluster survive.

### Filter 3: Leave-one-dataset-out (LOO) stability

For each cluster, each discovery dataset is removed in turn. The cluster is **stable** if >= 80% of its genes still meet the Phase 1 criterion (p < 0.05 in >= 2 of the remaining datasets) after any single dataset removal.

| Parameter | Default | Description |
|---|---|---|
| `stability_threshold` | `0.80` | Minimum retention fraction for stability |
| `p_threshold` | `0.05` | P-value threshold for retention check |
| `min_datasets` | `2` | Minimum datasets after removal |

Datasets whose removal causes > 20% gene loss are flagged as **fragile**.

### Filter 4: Core gene identification

Core genes are **Level 2 survivors** (p < 0.01 in >= 2 datasets) from clusters that have at least `min_genes_per_cluster` (default: 3) Level 2 survivors. Each core gene record includes:

- Cluster assignment
- Maximum sensitivity level survived (1--4)
- Per-dataset p-values and log2FC values
- Mean log2FC and predominant direction

**This core gene list is locked and recorded before Phase 5 begins.** No genes may be added or removed based on replication results.

### Inputs

- `Phase1Result`: Study genes with per-dataset DE statistics.
- `Phase3Result`: Cluster assignments and cluster info.
- `seed_gene_count`: Total number of seed genes (for FDR scope enforcement).
- `dataset_ids`: Discovery dataset IDs (for LOO analysis).

### Outputs

`Phase4Result` containing:

- `cluster_significance`: Permutation test results per cluster (raw p, Bonferroni p, observed statistic, null mean).
- `sensitivity`: Per-cluster sensitivity analysis (level survivors, dissolution point, core genes).
- `loo_stability`: Per-cluster LOO results (per-dataset retention fractions, stability flag, fragile datasets).
- `core_genes`: Dict of `gene_symbol -> CoreGene` (the locked pre-specification set).

### Configurable options

```yaml
phase4:
  n_permutations: 10000
  seed: 42
```

---

## Phase 5: Replication

### What it does

Tests each pre-specified core gene in held-out replication datasets for directional concordance. The core gene list was locked in Phase 4 and is recorded in the output. No genes are added or removed from the list based on replication results; the only action is to mark genes as survived or eliminated.

### Pre-specification protocol

1. The core gene list is finalized at the end of Phase 4 and stored as `locked_core_genes`.
2. Replication datasets are not accessed until Phase 5 begins.
3. The locked list is included in the output for audit purposes.

### Statistical method

For each core gene in each replication dataset, Welch's t-test is applied (same implementation as Phase 1) to compute log2FC and p-value. The direction (up/down) is compared to the discovery direction from Phase 1.

### Elimination rules

The elimination protocol applies tissue-aware directional concordance rules:

| Condition | Outcome | Rationale |
|---|---|---|
| Significant (p < 0.05) opposite direction in a **brain** dataset | **ELIMINATE** | Active contradictory evidence in matching tissue |
| Non-significant opposite direction in a brain dataset | **RETAIN** (noted) | Lack of power is not contradictory evidence |
| Any result in a **blood** dataset (concordant or discordant) | **RETAIN** | Brain-specific signals are expected to fail in blood; blood non-replication does not invalidate a brain discovery |
| No replication data available | **INSUFFICIENT DATA** | Gene is neither survived nor eliminated |

### Cluster-level verdicts

After gene-level elimination, each cluster receives a verdict:

| Verdict | Criterion |
|---|---|
| `replicated` | All core genes survived, no blood-specific failures |
| `brain_specific` | All core genes survived, but blood replication failed |
| `partially_replicated` | >= 50% of core genes survived |
| `failed` | < 50% of core genes survived, or none survived |

### Inputs

- `core_genes`: From `Phase4Result.core_genes` (locked, immutable).
- `replication_datasets`: Dict of `dataset_id -> DataFrame` (genes x samples).
- `replication_phenotypes`: Dict of `dataset_id -> {sample_id: 'case' or 'control'}`.
- `dataset_tissues`: Dict of `dataset_id -> 'brain' or 'blood'`.

### Outputs

`Phase5Result` containing:

- `gene_verdicts`: Per-gene elimination results (status, reason, per-dataset replication statistics, concordance counts).
- `cluster_verdicts`: Per-cluster replication verdicts.
- `locked_core_genes`: The sorted pre-specified core gene list (immutable record).
- `n_survived`, `n_eliminated`, `n_insufficient`: Summary counts.

### Configurable options

Replication datasets and tissue assignments are specified in the YAML configuration:

```yaml
datasets:
  - id: GSE64018
    role: replication
    tissue: brain
  - id: GSE18123
    role: replication
    tissue: blood
```

The elimination threshold is fixed at p < 0.05 and is not user-configurable (by design).

---

## Phase 6: Meta-Analysis

### What it does

Computes inverse-variance weighted (IVW) meta-analysis for each gene that survived Phase 5, pooling effect sizes across all discovery datasets. Produces both fixed-effects and random-effects estimates, heterogeneity statistics, and forest plot data.

### Statistical method

**Fixed-effects model** (IVW):
- Assumes a single true effect size across all datasets.
- Weights each study by the inverse of its variance (`w_i = 1 / SE_i^2`).
- Pooled effect: `sum(w_i * effect_i) / sum(w_i)`.
- Pooled SE: `1 / sqrt(sum(w_i))`.
- Significance via Z-test (normal distribution).

**Random-effects model** (primary result):
- Allows between-study heterogeneity in the true effect size.
- Tau-squared (between-study variance) is estimated using **REML** (Restricted Maximum Likelihood), with **DerSimonian-Laird** as fallback when REML does not converge.
- Study weights are adjusted: `w_i* = 1 / (SE_i^2 + tau^2)`.
- REML produces more accurate confidence intervals than DL when the number of studies is small (< 10), which is typical for transcriptomic meta-analysis.

**Heterogeneity statistics**:

| Statistic | Description |
|---|---|
| Cochran's Q | Weighted sum of squared deviations from the pooled estimate. Tested against chi-squared with k-1 degrees of freedom. |
| I-squared (I^2) | Percentage of total variability due to between-study heterogeneity rather than sampling error. Range: 0--100%. |
| Tau-squared (tau^2) | Estimated between-study variance. |

### SE recovery

If the standard error from Phase 1 is invalid (zero or non-finite), the engine attempts recovery from the log2FC, p-value, and sample sizes using the inverse of the t-distribution CDF, consistent with the Welch's t-test implementation.

### Expression scale check

Before pooling, the engine checks each dataset's expression scale by examining the maximum expression value. Datasets with max expression < 5.0 (indicating log-ratio format, e.g., GSE33000) are flagged with a scale warning. Log-ratio data has a compressed range (~-2 to ~2) compared to log2-intensity data (~4 to ~16), which affects SE recovery and meta-analysis weights.

### Inputs

- `Phase1Result`: Study genes with per-dataset DE statistics (log2FC, SE, p-value, sample sizes).
- `Phase5Result`: Gene verdicts (only `survived` genes are analyzed).
- `dataset_expression_ranges` (optional): Dict of `dataset_id -> max expression value`.

### Outputs

`Phase6Result` containing:

- `gene_results`: Per-gene meta-analysis results, each including:
  - Per-dataset effect sizes and standard errors (forest plot rows).
  - Fixed-effect and random-effect pooled estimates, SEs, and p-values.
  - Cochran's Q, I-squared, and tau-squared.
  - Overall direction and number of datasets.
- `n_genes_analyzed`: Genes with valid meta-analysis (requires >= 2 datasets).
- `n_significant_random`: Genes with random-effects p < 0.05.
- `n_high_heterogeneity`: Genes with I^2 > 75%.
- `scale_warnings`: Dataset IDs flagged for expression scale issues.

### Configurable options

Phase 6 has no user-configurable parameters. The random-effects model is always the primary result, and REML is always attempted before DL fallback.

---

## Pipeline-Level Design Decisions

These design choices are enforced across the pipeline and are not configurable:

| Decision | Enforcement | Consequence if violated |
|---|---|---|
| FDR uses full seed gene set | `apply_fdr_with_scope(scope="full_seed_set")` in Phase 4 | Significance inflation. AD stress test: 420 false positives at q < 0.10 with study-set FDR; 0 with full-seed-set FDR. |
| Welch's t-test uses exact t-distribution | `scipy.stats.t.sf()` in `riker/stats/welch.py` | Normal approximation underestimates p-values at small sample sizes. |
| Core gene list locked before replication | `locked_core_genes` recorded in `Phase5Result` | Post-hoc fishing through replication data undermines pre-specification. |
| Random effects is the primary meta-analysis result | REML with DL fallback in `riker/stats/meta.py` | Fixed effects overestimates precision when between-study heterogeneity exists. |
| Blood non-replication does not eliminate | Tissue-aware elimination rules in Phase 5 | Brain-specific discoveries would be incorrectly discarded. |
| log2FC range check | `validate_fold_changes()` in Phase 1 QC | Raw intensities produce biologically impossible fold changes (e.g., log2FC = 119). |
| Anti-circularity in feature construction | Pathway IDs only (no category labels) as clustering features in Phase 2 | Category labels create tautological enrichment. |
