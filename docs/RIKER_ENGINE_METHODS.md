# Riker Engine — Methods

**Version:** 1.0
**Last updated:** April 24, 2026
**Authors:** Ray Sigmon (Alpha Research Labs, Independent)
**Corresponding email:** alphalabsincorp@gmail.com
**Repository:** https://github.com/RaySigmon/Riker_Engine
**License:** AGPL-3.0
**Zenodo archive:** 10.5281/zenodo.19623672

---

## Abstract

The Riker Engine is a condition-agnostic transcriptomic discovery pipeline that identifies convergent differential expression signals across independent cohorts using a combination of per-gene statistical rigor, consensus clustering, stability profiling, and replication filtering. It is designed to surface robust molecular signatures from heterogeneous publicly-available transcriptomic data, with an explicit priority toward minimizing false positives. The engine runs on commodity hardware (validated on a Raspberry Pi 5) and produces fully-reproducible results from a fixed master seed. This document describes the engine's architecture, statistical methodology, design rationale, output interpretation, and explicit limitations.

---

## Conceptual overview

Most transcriptomic differential expression analyses operate on a single cohort at a time, identify genes that differ between cases and controls, then interpret what those genes mean biologically. This works well when you have a large, clean cohort. It works less well when you want to identify genes that are *consistently* affected by a disease across independent studies from different labs, different platforms, and different patient populations.

The Riker Engine is built for the multi-cohort convergence question: "which genes show the same directional change across multiple independent studies of the same disease?" Answering this cleanly requires more than running differential expression on pooled data. It requires understanding which signals are reproducible across study designs and which are artifacts of any single cohort.

The engine addresses this through six sequential phases, each answering a specific question:

1. **Per-cohort differential expression** — for each seed gene, is it significantly different between cases and controls in each cohort independently?
2. **Feature construction for clustering** — summarize each gene's cross-cohort behavior in a feature matrix
3. **Consensus clustering** — group genes by similar cross-cohort behavior using a bootstrap-consensus approach over UMAP+HDBSCAN
4. **Robustness testing** — apply per-gene statistical thresholds, identify which genes co-occur with other robust genes, assess cluster-level significance
5. **Replication filtering** — test candidate genes in held-out cohorts, eliminate genes with significant opposite-direction expression in same-tissue replication
6. **Meta-analysis** — compute pooled effect sizes with heterogeneity statistics across discovery cohorts

The output is a set of **core genes** that have passed every filter. A separate process, **stability profiling**, repeats the pipeline 50 times with varying random seeds to measure how reproducibly each core gene emerges under stochastic variation — distinguishing **iron-clad genes** (appearing in ≥90% of runs) from genes that are artifacts of any particular random seed.

The design philosophy prioritizes specificity over sensitivity. It is better to miss a real gene than to include a false positive. Users who want maximum sensitivity should supplement the Riker Engine's output with less-conservative methods.

---

## Architecture

### Phase 1: Cross-Cohort Differential Expression

**Input:** Per-cohort expression matrices (pre-normalized, log2-scale) and case/control phenotype labels.

**Procedure:**
- For each seed gene in each discovery cohort, compute Welch's t-test (unequal-variance two-sample t-test) between cases and controls
- Compute mean difference in log2-scale expression (`mean(cases) - mean(controls)`), which under log2-normalized input is equivalent to log2 fold change
- Retain genes significant at nominal p < 0.05 in at least 2 discovery cohorts

**Output:** A study gene set — genes showing coordinated case-vs-control differences across multiple independent cohorts. Typically 5-20% of input seed gene count; much smaller than the original seed space.

**Statistical details:**
- Welch's t-test uses the exact t-distribution via SciPy (`scipy.stats.t`)
- NaN values and insufficient-sample groups (<2 samples per group) are skipped
- QC guardrail flags any |log2 fold change| > 10 as potential non-log2 data

### Phase 2: Feature Matrix Construction

**Input:** Phase 1 study genes with their per-cohort statistics.

**Procedure:** Build a per-gene feature matrix for downstream clustering. Current implementation uses five summary expression features per gene:
- Average log2 fold change across cohorts
- Negative log10 of the minimum p-value across cohorts
- Number of cohorts where the gene is significant
- Direction consistency (fraction of cohorts with concordant direction)
- Tier score (composite rank)

All features are min-max normalized to [0, 1].

A pathway-membership feature infrastructure exists in the code but is not currently used by any disease config. When populated, binary pathway-membership features (e.g., KEGG, Reactome) are filtered to pathways containing ≥3 study genes and ≤500 total genes, then added to the feature matrix.

**Output:** Feature matrix with rows = genes, columns = features. Feeds into Phase 3 clustering.

### Phase 3: Consensus Clustering

**Input:** Phase 2 feature matrix.

**Procedure:**
- Run UMAP dimensionality reduction + HDBSCAN density-based clustering across multiple parameter configurations (default: 3 `n_neighbors` values × 5 random seeds = 15 configurations)
- For each gene pair, compute the fraction of configurations in which both genes were assigned to the same cluster (ignoring noise assignments)
- Apply final HDBSCAN on the consensus distance matrix (1 - co-association fraction, using `metric="precomputed"`)

**Parameters (defaults):**
- UMAP: `n_components=2`, `min_dist=0.0`, `metric="euclidean"`
- `n_neighbors`: `[10, 15, 30]`
- UMAP random seeds: `[42, 123, 456, 789, 1024]`
- HDBSCAN: `min_cluster_size=5`, `min_samples=3`

**Output:** Cluster assignments per gene. Genes not confidently assigned are marked as noise (cluster -1).

**Rationale for consensus:** Single UMAP+HDBSCAN runs are stochastic and parameter-sensitive. The consensus approach identifies gene groupings that are robust across parameter variation, filtering out spurious structure from any single configuration.

### Phase 4: Robustness Testing

**Input:** Phase 3 cluster assignments + Phase 1 per-gene statistics.

Four independent sub-steps run here, each answering a different question:

**Sub-step 4.1 — Permutation significance per cluster:**
- For each cluster, compute mean of |log2 fold change| across its member genes
- Permute gene-cluster assignments 10,000 times, compute the null distribution of cluster mean |log2FC|
- Apply Bonferroni correction across clusters
- Report clusters passing p < 0.05 (Bonferroni-corrected) as **significant clusters**

**This test is reported as a diagnostic. It is NOT used to gate core gene selection.** See "Design decisions" below for rationale.

**Sub-step 4.2 — Sensitivity analysis:**
Per-gene progressive significance filter across four thresholds:
- Level 1: p < 0.05 in at least 2 datasets
- Level 2: p < 0.01 in at least 2 datasets (the per-gene criterion)
- Level 3: FDR q < 0.10 in at least 2 datasets
- Level 4: FDR q < 0.05 in at least 2 datasets

FDR correction uses Benjamini-Hochberg with scope-enforcement: p-values are padded with 1.0 up to the full seed gene count before correction, preventing inflation from selective testing.

**Sub-step 4.3 — Leave-one-out stability:**
- For each discovery dataset, remove it and re-test which genes still meet the per-gene significance criterion (p < 0.05 in ≥2 of the remaining datasets)
- Report genes retained at ≥80% across all leave-one-out configurations
- Note: this is a *statistical* LOO (recomputes significance thresholds), not a *clustering* LOO (does not re-run UMAP+HDBSCAN)

**Sub-step 4.4 — Core gene identification:**
- Core genes = Level 2 survivors (p < 0.01 in ≥2 datasets) from clusters containing at least 3 Level 2 survivors
- This is the output that feeds Phase 5

The ≥3-per-cluster threshold is an anti-singleton guard: it requires that a gene's statistical robustness be corroborated by at least 2 other genes in the same co-expression module. Genes passing individual significance but sitting in tiny clusters (size 1 or 2) are excluded because their robustness lacks corroboration from neighboring genes.

**Output:** The core gene set for this run.

### Phase 5: Replication Filtering

**Input:** Phase 4 core genes + held-out replication cohort expression matrices.

**Procedure:**
- For each core gene, compute Welch's t-test in each replication cohort
- Classify the gene's direction in each replication cohort (concordant with discovery, discordant, or insignificant)
- **Elimination criterion:** A gene is eliminated from the core set if it shows *significant* expression in an *opposite direction* in a *same-tissue* replication cohort. A single discordant-significant same-tissue cohort triggers elimination.
- Genes eliminated here are permanently excluded from the final reported core set for this run

**Tissue logic:** The code supports a cross-tissue tolerance mode — when replication is in a different tissue from discovery, the elimination criterion is relaxed. This reflects the biological reality that, e.g., brain and blood can show opposite directions for tissue-specific metabolic genes. **Current pipeline invocations run under the conservative default where all replication is treated as same-tissue**, making the elimination criterion maximally strict. This is a latent plumbing issue slated for a v0.3.3 fix; in the interim, all results have been produced under the stricter filter.

**Output:** Surviving core genes for final meta-analysis.

### Phase 6: Meta-Analysis

**Input:** Phase 5 survivors + their per-cohort effect sizes from discovery cohorts only (replication cohorts are deliberately excluded to avoid double-counting evidence).

**Procedure:**
- For each surviving gene, compute inverse-variance weighted meta-analysis across discovery cohorts
- Produce both fixed-effects and random-effects estimates
- Random-effects method: REML for k ≥ 3 cohorts (the typical case), DerSimonian-Laird fallback if REML fails to converge
- Compute Cochran's Q, I², and τ² for heterogeneity
- Apply Benjamini-Hochberg FDR correction across all tested genes

**Output:** Pooled effect sizes, confidence intervals, p-values (fixed and random), q-values, and heterogeneity statistics per gene.

### Stability Profiling

A wrapper process that runs the entire six-phase pipeline 50 times with deterministically varied seeds:

- Seeds are derived from `master_seed + run_num` via `random.Random(master_seed + run_num).randint(0, 2**31 - 1)` for each run's base seed, plus 5 derived seeds for Phase 3 UMAP
- Each run produces a `phase4_core_genes.csv` containing that run's core gene set
- After all 50 runs, genes are classified by frequency of appearance:
  - **Iron-clad:** appears in ≥90% of 50 runs
  - **Borderline:** appears in 50-89% of runs
  - **Stochastic:** appears in <50% of runs
- Note: stability profiling measures Phase 4 output (pre-replication). A gene can be iron-clad in stability but eliminated in Phase 5 replication.

Stability profiling is the engine's primary reproducibility test. Iron-clad genes are the genes the engine is most confident about.

---

## Statistical methodology

### Choice of Welch's t-test

Welch's unequal-variance t-test is used in Phases 1 and 5. Rationale: differences in variance between cases and controls are biologically expected (diseased tissue often has higher expression heterogeneity than normal tissue). Student's t-test assumes equal variances and is more powerful when that assumption holds, but produces inflated false-positive rates when it doesn't. Welch's is the safer default for case-control expression data.

### Choice of Bonferroni for cluster permutation

Bonferroni correction is used on the Phase 4 permutation test because the number of clusters is typically small (tens to ~200) and Bonferroni's conservatism is appropriate at that scale. FDR would be more powerful but also more permissive, and the engine's philosophy prioritizes specificity.

### Choice of REML for meta-analysis

Random-effects meta-analysis uses REML (Restricted Maximum Likelihood) as the primary estimator because it has better small-sample properties than DL at typical k=3 cohort counts. DL is used as a fallback only when REML fails to converge — which is rare on the kinds of effect-size data the engine produces.

### Scope-enforced FDR

The FDR correction in Phase 4 pads the p-value vector with 1.0 up to the full seed gene count before correction. This prevents an inflated rejection rate from selectively testing only the genes that passed Phase 1. Without this padding, BH correction would be applied to a pre-filtered subset, artificially inflating discovery rates. This guardrail is essential for the engine's FDR claims to be meaningful.

### Fixed master seed for reproducibility

Every run's stochastic behavior is controlled by `phase3.seeds` (UMAP random states) and `phase4.seed` (permutation base seed), both specified in the YAML config. The default Phase 3 seeds are `[42, 123, 456, 789, 1024]` and the default Phase 4 seed is `42`. Two runs with the same seed configuration, same input data, and same code should produce identical outputs. Stability profiling wraps this by programmatically varying these seeds (derived from `master_seed + run_number`) to explore the variance landscape rather than eliminating it.

The reproducibility claim is about deterministic reproduction of the same result. Stability profiling then quantifies *how much* the result varies under seed variation, giving users a concrete measurement rather than just assurance.

---

## Design decisions and their rationale

### Why the cluster permutation test does NOT gate core gene selection

Phase 4.1 computes per-cluster permutation significance and reports it. Phase 4.4 identifies core genes using per-gene Level 2 survivorship + a cluster-size threshold (≥3 Level 2 survivors per contributing cluster). Cluster-level permutation significance is NOT a filter on core gene selection.

This is a deliberate design choice, not an omission. The rationale:

- **The per-gene criterion is the statistical rigor layer.** A gene passing p<0.01 in ≥2 independent cohorts has demonstrated its own statistical robustness. Requiring additionally that its parent cluster show elevated mean |effect size| over random groups would be a composite filter testing two different things.
- **The cluster-size threshold is the corroboration layer.** A gene's individual statistical robustness should be corroborated by the presence of at least 2 other robust genes in the same co-expression module. This catches lone-wolf genes whose significance lacks neighboring support.
- **Cluster mean |effect size| is not a coherent proxy for "convergent biology."** Mean |log2FC| of a cluster depends on cluster composition, not just signal strength. A cluster of one huge-effect gene plus two modest-effect genes can fail permutation even if all three genes are individually meaningful.

The result is that the engine identifies **iron-clad gene sets** whose biological coherence is established through enrichment analysis against external pathway databases, not through single-cluster identification by the pipeline itself. When a biological theme emerges (e.g., mitochondrial metabolism, immune signaling), it is a post-hoc annotation of the gene set, not a computational cluster identified by the pipeline. See "Output interpretation" below for the precise terminology.

### Why the pipeline is fast enough to run on a Raspberry Pi

Published results have been produced in 5-10 minutes of wall clock per pipeline run on a Raspberry Pi 5 (8GB RAM, ARM64). This reflects several deliberate design decisions that collectively give the engine a smaller computational footprint than tools like WGCNA:

- **Phase 1 pre-filtering is aggressive.** Typical runs reduce 17,000-26,000 input genes to 1,500-5,000 study genes before any clustering. All downstream phases operate on this reduced set.
- **Consensus clustering is structural, not iterative.** WGCNA iteratively optimizes soft-threshold parameters for a full gene-gene correlation matrix, which is O(n²) in genes and requires many iterations. Riker's consensus runs UMAP+HDBSCAN 15 times (fixed) on a small feature matrix and combines results.
- **HDBSCAN scales with samples, not features.** For clustering, genes are the "samples." A few thousand samples is small by HDBSCAN standards.
- **No iterative model fitting.** Phase 6 meta-analysis uses closed-form IVW + REML, which is fast.

The engine is not fast because it cuts corners. It is fast because the algorithmic structure works in a smaller computational regime than the tools it is often compared against. Benchmarking against WGCNA and MEGENA is ongoing.

### Why all-expressed seed sets are no longer recommended

Early iterations of the engine used per-disease "all expressed" seed lists, generated by taking the union of gene identifiers detected across each disease's GEO microarray platforms. This included non-coding transcripts, pseudogenes, and platform-specific probe identifiers. The approach was data-faithful but methodologically messy — seed list composition varied by disease depending on which platforms each cohort used.

Current methodology uses a standardized protein-coding gene list (`all_protein_coding_genes.csv`, ~19,296 HGNC protein-coding genes) for blind runs across all diseases. This provides methodological consistency across disease studies and is the canonical seed set going forward. Historical all-expressed runs remain valid reproducible artifacts but are not directly comparable to current runs.

---

## Output interpretation

The engine's outputs use specific terminology that must be preserved accurately in any paper or presentation citing the engine. The following terms have fixed, binding meanings:

- **Iron-clad gene:** A gene appearing in at least 90% of 50 stability-profile runs. This is the engine's strongest reproducibility designation.
- **Borderline gene:** Appears in 50-89% of stability runs. Moderately reproducible.
- **Stochastic gene:** Appears in <50% of stability runs. Not reliably reproducible under seed variation.
- **Core gene (single run):** A gene surviving Phase 4 filters in a single pipeline invocation. Exists only in the context of one specific run.
- **Level 2 survivor:** A gene passing p<0.01 in ≥2 discovery cohorts. The per-gene statistical criterion.
- **Significant cluster:** A cluster passing Bonferroni-corrected permutation significance in Phase 4.1. Reported as a diagnostic, NOT a gate on core gene selection. The number of significant clusters will almost always be smaller than the number of clusters contributing core genes.
- **Biological annotation:** A post-hoc functional grouping of iron-clad genes using external pathway or function databases (e.g., KEGG, Reactome, GO). NOT the same as a computational cluster identified by the pipeline. Biological annotations are typically more coarse than computational clusters and may span many clusters.

**A finding described as "the pipeline identified an X-gene cluster of [biological function]" is a methodological error** unless all X genes actually originated from a single Phase 3 consensus cluster. Most real-world biological findings span multiple clusters and are unified by post-hoc annotation. Papers must describe such findings as "an iron-clad gene set enriched for [biological function]" or similar language that accurately reflects the process.

---

## Performance Characteristics

The Riker Engine completes a full six-phase pipeline run on a Raspberry Pi 5 (8GB RAM, ARM64)
in 2-8 minutes for typical diseases. This section documents measured timing from the first
SOP-compliant ASD disease-day (April 26, 2026) and explains the architectural decisions that
make this possible.

### Measured per-phase timing (ASD, 2026-04-26, Pi 5)

Wall-clock measurements extracted from per-phase log timestamps across three runs of the same
disease under different seed-set configurations. Source: `disease_days/2026-04-24_asd/{run}/run.log`.

| Phase | Curated (1,267 seeds) | Blind protein-coding (19,296 seeds) | All-expressed (39,185 seeds) |
|---|---|---|---|
| 1: Cross-referencing (Welch's t per gene per dataset) | 25s | 79s | 81s |
| 2: Feature matrix construction | <1s | <1s | <1s |
| 3: Consensus clustering (15× UMAP+HDBSCAN) | 35s | 190s | 191s |
| 4: Robustness testing (10K permutations × N clusters) | 11s | 163s | 159s |
| 5: Replication filtering (4 replication datasets) | 45s | 36s | 35s |
| 6: Meta-analysis (closed-form IVW + REML/DL) | <1s | <1s | <1s |
| **Total** | **1.9 min** | **7.8 min** | **7.8 min** |

**Phase 2 and Phase 6 take <1 second each** — confirmed by identical timestamps at phase
boundaries in the run logs. Phase 2 builds a 5-feature numeric matrix per gene (no pathway
database is configured in the current pipeline). Phase 6 applies closed-form inverse-variance
weighted meta-analysis — no MCMC, no iterative fitting.

**Phase 3 + Phase 4 dominate blind runs** (76% of wall clock). Phase 3 runs 15 UMAP+HDBSCAN
configurations on the study gene set; Phase 4 runs 10,000 permutation tests per cluster.

**Phase 5 dominates the curated run** (39% of wall clock) because replication dataset loading
is a fixed I/O cost regardless of core gene count (35 curated vs 401 blind).

### Why it's fast: the Phase 1 funnel

The single most important performance fact is that **Phase 1 eliminates 89-96% of input genes
before any clustering happens.** All downstream phases operate on the filtered study gene set,
not the original transcriptome.

| Run | Input seeds | → Phase 1 study genes | Reduction |
|---|---|---|---|
| Curated | 1,267 | 141 | 89% eliminated |
| Blind protein-coding | 19,296 | 1,793 | 91% eliminated |
| All-expressed | 39,185 | 1,813 | 95% eliminated |

The all-expressed run starts with twice as many seeds as the protein-coding blind run but
produces only 20 more study genes (1,813 vs 1,793). The additional ~20,000 non-protein-coding
identifiers in the all-expressed seed list are not detected in expression matrices and are
eliminated immediately. This is why the all-expressed run takes the same time as the
protein-coding blind run (7.8 min vs 7.8 min) — the effective working set is identical.

### Why it's fast: architectural comparison to WGCNA

WGCNA computes a full soft-thresholded gene-gene correlation matrix (O(n²) in genes),
iteratively optimizes topology parameters, and requires the full transcriptome to remain
in memory throughout. For the same ASD datasets:

- **WGCNA**: ~3h 50min on Pi 5, OOM-crashed on 18K genes (required filtering to 10K),
  produced 1,427-9,995 genes per significant module
- **Riker Engine**: 7.8 min on Pi 5, ran full 19K gene set without issue, produced 401
  core genes

The engine avoids the O(n²) bottleneck by never computing a full gene-gene correlation matrix.
Instead, it reduces genes to a small feature matrix (5 features per gene), runs density-based
clustering (HDBSCAN) on the 2D UMAP embedding of that matrix, and builds consensus from 15
independent UMAP configurations. The co-association matrix is n×n, but n is the study gene
count (~1,800), not the full transcriptome (~19,000).

### The gene funnel (selectivity)

Speed is not the goal — specificity is. Each phase is a statistical filter:

```
19,296 protein-coding seeds
  → Phase 1:  1,793 study genes    (p<0.05 in ≥2 of 3 discovery datasets)
  → Phase 4:    401 core genes     (p<0.01 in ≥2, with ≥3 per cluster)
  → Phase 5:    220 survived       (replicated in held-out cohorts)
  → Phase 6:    135 meta-significant (random-effects p<0.05)
```

Cumulative pass-through: 1.1% of seeds reach Phase 5 survival. 0.7% reach meta-analytical
significance. The engine prioritizes specificity over sensitivity — it is better to miss a
real gene than to include a false positive.

For full statistical methodology, see `docs/RIKER_ENGINE_METHODS.md`.


---

## What the engine does NOT claim

Clear scope boundaries, to prevent overinterpretation:

- The engine does not claim that iron-clad genes are *causal* in the disease. They are convergent in expression across independent cohorts. Causation requires experimental validation (knockout, knock-in, perturbation studies) that the engine does not perform.
- The engine does not claim to identify *all* disease-relevant genes. It identifies genes that are robust across the specific cohorts used. Genes present in only one cohort, or with effects smaller than the per-gene significance threshold, will be missed. False negatives are an accepted cost of the specificity-first design.
- The engine does not replace targeted differential expression analysis. For deep exploration of a single cohort, tools like DESeq2, limma-voom, and edgeR are appropriate. The engine addresses a different question — cross-cohort convergence — and complements these tools rather than replacing them.
- The engine does not claim its clustering step identifies biologically-interpretable modules. Phase 3 clustering is structural — it provides the substrate for the cluster-size filter and for per-cluster reporting. Biological interpretation is a post-hoc step using enrichment analysis against curated pathway databases.
- The engine does not support cross-tissue studies in its current state (v0.3.2.1). The cross-tissue tolerance logic exists in code but is not passed from the CLI. Users running cross-tissue studies should be aware that the current conservative-default behavior may over-eliminate genes whose tissue-specific expression patterns legitimately differ between discovery and replication.

---

## Validation status

As of April 24, 2026, the following validations have been completed:

| Disease | Curated run | Protein-coding blind run | All-expressed blind run (historical) | 50-run stability |
|---|---|---|---|---|
| ASD | ✅ 35 core | ✅ 403 core (run_001) | ✅ 39,185 seeds → historical result | ✅ 376 iron-clad |
| T2D | ✅ 8 core | (SOP-compliant run pending) | ✅ 26,800 seeds → 177 core (archived) | (pending) |
| IBD | ✅ 304 core | (SOP-compliant run pending) | ✅ 24,651 seeds → historical result | (pending) |
| AD | ✅ 394 core | ✅ (earlier blind run) | (not done) | (pending) |
| Breast Cancer | ✅ 152 core | ✅ (earlier blind run) | (not done) | (pending) |
| IPF | ✅ 190 core | (SOP-compliant run pending) | (not done) | (pending) |
| Psoriasis | ✅ 50 core | ✅ (earlier blind run) | (not done) | (pending) |
| CRC | ✅ 264 core | ✅ (earlier blind run) | (not done) | (pending) |

**Cross-disease remarks:**
- The ASD 50-run stability profile (376 iron-clad genes) is the only fully SOP-compliant disease to date. Per-disease SOP-compliant validation is underway in the April-May 2026 window.
- Earlier blind runs (AD, Breast Cancer, CRC, Psoriasis, ASD, T2D, IBD) were conducted at various stages of engine development and used inconsistent seed-set definitions. Archived configs and output are preserved at `~/riker-archive/old-iterations/`. These are treated as historical validation evidence, not canonical current-state results.
- Independent cold-start reproducibility has been verified for the v0.3.2 engine by four independent validators (Jim R1, Jim R2, Jim R3 scrubbed environment, Claude Code agent) across 6 diseases in February-March 2026. All four produced matching core gene sets on the ASD curated pipeline. The v0.3.2 commit `801693d` captures this validated state.

A full re-validation of all 8 diseases under the Disease-Day SOP is the current focus of the April-May 2026 work window.

---

## Known limitations and v0.3.3 roadmap

The following items are known issues scheduled for the v0.3.3 release. None of them affect the correctness of published results, but they affect usability, portability, and documentation clarity.

- **`config.random_seed` field is dead.** The top-level `random_seed` in YAML configs is parsed but not consumed by any phase. No `--random-seed` CLI flag exists either. The actual stochastic controls are `phase3.seeds` (UMAP random states) and `phase4.seed` (permutation base seed), both set explicitly in the YAML config. The SOP-compliant workaround is to set all three fields explicitly in per-run disease-day configs. Slated fix: wire the YAML `random_seed` field to derive `phase3.seeds` and `phase4.seed` defaults automatically.
- **`discovery_tissues` is not plumbed to Phase 5.** Cross-tissue tolerance logic exists but never activates. Current behavior is conservative (same-tissue by default). Slated fix: pass `discovery_tissues` from CLI through the pipeline.
- **Config portability.** Archived configs use absolute `/home/kai001/` paths. Slated fix: relative path resolution with environment-variable override support.
- **Phase output standardization.** Some phases' intermediate output files are inconsistent across versions. Slated fix: standardize the required output set per phase (see SOP §"Required outputs per run").
- **Code version in pipeline output.** `pipeline_summary.json` does not embed the git commit hash that produced it. Slated fix: inject `code_version` at run start.
- **All-expressed config templates.** The archived all-expressed runs use seed-file paths and dataset conventions that don't match current relative-path standards. Slated fix: runnable all-expressed config template per disease with current paths.

---

## Usage

For operational procedure on running a disease validation, see `docs/DISEASE_DAY_SOP.md`.

For quick start and installation, see `README.md`.

---

## Citation

If you use the Riker Engine in your research, please cite:

> Sigmon, R. (2026). Riker Engine: A stability-profiling transcriptomic convergence pipeline. Zenodo. https://doi.org/10.5281/zenodo.19623672

And the associated disease-specific paper(s) as appropriate.

Alpha Labs / Alpha Research Labs
Independent computational research
Contact: alphalabsincorp@gmail.com

---

*This document is authoritative for the Riker Engine methodology. Changes are tracked through git commits and should include written rationale in the commit message.*
