# THE RIKER ENGINE: A Condition-Agnostic Transcriptomics Pipeline for Discovering Replicated Gene Modules from Public Expression Data

**Ray Sigmon**
Independent Computational Researcher, Northwest Arkansas, USA

Correspondence: alphalabsincorp@gmail.com

**Version 0.3.2 — April 2026**

Repository: github.com/RaySigmon/Riker_Engine | License: AGPL-3.0

---

## Abstract

Identifying reproducible gene expression signatures across independent transcriptomic datasets remains a central challenge in disease biology. While thousands of expression studies exist in public repositories, synthesizing their results into statistically defensible, replicated gene modules requires extensive manual effort and specialized expertise. This work presents the Riker Engine, an open-source Python pipeline that automates cross-dataset gene module discovery through a six-phase progressive filtering framework: cross-dataset differential expression filtering, pathway-informed feature matrix construction, consensus clustering via UMAP and HDBSCAN across 15 parameter configurations, permutation-based robustness testing with full-seed-set false discovery rate correction, independent directional replication, and inverse-variance weighted meta-analysis with random-effects estimation.

The engine was validated on six diseases spanning five tissue types with zero code modifications between runs: autism spectrum disorder (35 core genes, 15 eliminated for blood discordance confirming brain specificity, 77.1% blind recovery), type 2 diabetes (8 core genes, IAPP identified as top blind-run signal from 26,800 genes), inflammatory bowel disease (304 core genes, 97.7% blind recovery), Alzheimer's disease (394 core genes, 98.2% blind recovery), breast cancer (152 core genes, 99.3% blind recovery, independent subtype separation), and idiopathic pulmonary fibrosis (190 core genes, 86.3% cold replication in a held-out 122-patient cohort). Two additional diseases — Psoriasis (50 core genes) and Colorectal Cancer (244 core genes) — were validated entirely by independent AI agents with no author involvement in disease selection, seed gene curation, or dataset choice. A head-to-head benchmark against WGCNA demonstrated superior output specificity (35 vs. 1,427–9,995 genes per module), built-in cross-dataset validation, and reduced hardware requirements. The complete pipeline executes in approximately 8 minutes on a Raspberry Pi 5 with 8 GB RAM. Cold-start replication tests confirmed deterministic phases reproduce bit-for-bit across independent environments, with stochastic variance limited to 0.5% of the core gene set.

---

## 1. Introduction

### 1.1 The Integration Gap

Public repositories such as the Gene Expression Omnibus (GEO) contain thousands of transcriptomic datasets covering hundreds of conditions. For many diseases, multiple independent expression studies exist comparing affected and unaffected tissue. Yet these datasets are rarely integrated in a systematic, reproducible manner. Individual studies produce gene lists that may or may not replicate across cohorts. Meta-analyses, when performed, typically aggregate differential expression statistics without identifying the modular structure of disease biology.

Network-based approaches such as Weighted Gene Co-expression Network Analysis (WGCNA) (Langfelder and Horvath, 2008) identify co-expressed gene modules but operate on single datasets, produce modules containing thousands of genes, and require extensive manual post-processing to extract actionable targets. No widely adopted tool bridges the gap between single-dataset network analysis and multi-dataset statistical validation to yield a minimal, replicated gene set suitable for experimental follow-up.

### 1.2 Design Philosophy

The Riker Engine occupies this gap. It accepts a list of candidate genes for any disease, queries multiple independent expression datasets, and returns a core gene set that has survived progressive statistical filtering, cross-dataset consensus clustering, permutation testing, pre-specified independent replication, and effect size meta-analysis.

The engine is designed as a validation and filtering pipeline — not a black-box discovery tool. Given candidate genes and independent datasets, it determines which candidates show consistent, replicated differential expression across studies with proper statistical controls. Novel candidates emerge through co-clustering with seed genes, but the primary value proposition is reproducible cross-dataset validation.

Key design principles include condition-agnostic operation (identical parameters across all diseases), progressive honesty (every filter reports what survives and what fails), pre-specified replication (gene lists locked before replication data is accessed), honest nulls (the pipeline reports when data does not support a finding), and accessibility (runs on consumer hardware without institutional infrastructure).

### 1.3 Contributions

This work makes four principal contributions: (1) a six-phase progressive filtering pipeline integrating pathway-informed clustering with cross-dataset replication; (2) a consensus clustering framework converting parameter-sensitive UMAP embeddings into parameter-robust co-association matrices; (3) validation across eight diseases, six tissue types, and over 40 independent datasets with zero code modifications; and (4) independent replication by third-party agents confirming reproducibility and generalizability.

---

## 2. Methods

### 2.1 Overview

The Riker Engine implements a six-phase progressive filtering pipeline. Each phase reduces the gene set while increasing the statistical burden. All thresholds are fixed across diseases; no per-disease parameter tuning is performed.

### 2.2 Phase 1: Cross-Dataset Differential Expression

Welch's t-test (Welch, 1947) is applied independently to each gene in each dataset, comparing case and control groups on log2-transformed expression values. Genes reaching significance at p < 0.05 in two or more independent datasets advance to the study gene set.

The choice of Welch's t-test reflects three design priorities. First, uniformity: Welch's t-test operates identically on log2-transformed microarray and RNA-seq data. Second, conservatism: on log2-transformed data, Welch's t-test is less powerful than limma, producing fewer false positives at the cost of more false negatives. Since Phase 1 is a filter rather than the final result, conservatism is acceptable. Third, independence: unlike limma, which models variance across genes, Welch's t-test treats each gene independently, avoiding assumptions about cross-gene variance structure that may not hold across platforms.

The exact t-distribution is used rather than the normal approximation because small sample sizes (n = 8–15 per group in some datasets) render the normal approximation unreliable.

### 2.3 Phase 2: Feature Matrix Construction

Study genes are mapped to biological pathways from KEGG (Kanehisa and Goto, 2000), Reactome (Fabregat et al., 2018), and MSigDB Hallmark gene sets (Liberzon et al., 2015). The feature vector per gene combines binary pathway memberships with five expression statistics: mean log2 fold change, negative log10 minimum p-value, number of significant datasets, directional consistency, and confidence tier. All features are min-max normalized to [0,1].

An anti-circularity rule is enforced: individual pathway IDs are used as features, but pre-assigned biological category labels are never used as features. This prevents tautological enrichment.

Note: All published validation results were produced using expression-based features only (without pathway data configured). Pathway integration is available as an optional enhancement but does not affect the validated results presented here.

### 2.4 Phase 3: Consensus Clustering

UMAP (McInnes et al., 2018) with n_components=2 and min_dist=0.0 followed by HDBSCAN (Campello et al., 2013) with min_cluster_size=5 and min_samples=3 is applied across 15 configurations: 3 n_neighbors values (10, 15, 30) × 5 random seeds (42, 123, 456, 789, 1024). A co-association consensus matrix records pairwise co-clustering frequencies. Final cluster labels are derived by applying HDBSCAN with metric='precomputed' to the distance matrix (1 − consensus).

The consensus approach converts a parameter-sensitive method into a parameter-robust one. Any single UMAP + HDBSCAN run is dependent on specific parameters and random seed. The consensus matrix averages across all 15 configurations, requiring gene pairs to consistently co-cluster across multiple parameter settings.

### 2.5 Phase 4: Robustness Testing and Core Gene Identification

Four progressive filters are applied. First, permutation significance testing (10,000 permutations with Bonferroni correction). Second, sensitivity analysis at four stringency levels: p < 0.05, p < 0.01, FDR q < 0.10, and FDR q < 0.05, all using the full seed set as the FDR denominator. Third, leave-one-dataset-out stability analysis requiring minimum 80% gene retention. Fourth, core genes are identified as Level 2 or higher survivors with three or more cluster-mates.

The full-seed-set FDR correction is a critical design decision. Most tools compute FDR only over genes that passed an initial filter. This inflates significance because non-significant genes are excluded from the denominator. The Riker Engine computes Benjamini-Hochberg FDR over the entire seed set. In stress testing, study-set-only FDR produced 420 false positives while full-seed-set FDR correctly returned zero.

The core gene list is locked before any replication data are examined.

### 2.6 Phase 5: Independent Replication

Each core gene is tested in held-out replication datasets for directional concordance. Genes showing significant (p < 0.05) opposite-direction effects in same-tissue datasets are eliminated. Cross-tissue non-replication is informative rather than disqualifying.

### 2.7 Phase 6: Effect Size Meta-Analysis

Inverse-variance weighted meta-analysis pools effect sizes across all datasets. Both fixed-effects and REML random-effects estimates are computed. REML was selected as the primary estimator because the DerSimonian-Laird method underestimates tau-squared when the number of studies is small (fewer than 10). Heterogeneity is quantified using Cochran's Q, I-squared, and tau-squared.

### 2.8 Quality Control Framework

Mandatory QC checks include log2 transformation detection, circularity auditing, permutation testing, FDR scope enforcement, sensitivity analysis, LOO stability, pre-specification locking, elimination protocol enforcement, and z-score verification. These checks were developed after catching errors during development that would have produced false findings.

---

## 3. Validation Results

### 3.1 Six-Disease Author Validation

The engine was validated on six diseases with zero code modifications between runs.

| Disease | Tissue | Seeds | Datasets | Ph1 Yield | Core Genes | Ph5 Survival | Meta-Sig | Blind Recovery |
|---------|--------|-------|----------|-----------|------------|-------------|----------|----------------|
| ASD | Brain cortex | 1,267 | 7 | 11.1% | 35 | 57.1% | 9 | 77.1% |
| T2D | Pancreatic islets | 443 | 4 | 12.6% | 8 | 100% | 8 | 62.5% |
| IBD | Intestinal mucosa | 762 | 6 | 53.5% | 304 | 99.3% | 296 | 97.7% |
| AD | Brain cortex | 801 | 5 | 54.7% | 394 | 86.3% | 312 | 98.2% |
| Breast Cancer | Breast tumor | 653 | 5 | 34.8% | 152 | 100% | 121 | 99.3% |
| IPF | Lung tissue | 354 | 5 | 68.1% | 190 | 89.5% | 157 | N/A |
| Psoriasis* | Skin | 96 | 5 | 62.5% | 50 | 100% | 28 | N/A |
| CRC* | Colon | 515 | 6 | 64.3% | 263 | 92.8% | 218 | N/A |

\*Psoriasis and CRC were validated by independent AI agents with no author involvement in disease selection, seed gene curation, or dataset choice.

Gene yield scales with transcriptomic signal strength — from 8 genes (T2D, subtle metabolic) to 394 genes (AD, strong neuroinflammatory). The engine calibrates to the biology; it does not force a result.

### 3.2 Selected Validation Highlights

**ASD:** The weakest transcriptomic signal. From 1,267 SFARI candidates, the engine identified 35 core genes. 15 were eliminated for significant opposite-direction effects in blood replication datasets — confirming brain specificity and demonstrating the pipeline's ability to distinguish tissue-specific signals from noise. 20 survived, including ATP2B2 and SEZ6L2, which showed significant brain-specific downregulation across four independent postmortem cohorts (random-effects p = 0.037 and 0.046 respectively). The hypothesis-free blind run recovered 27 of 35 core genes (77.1%) despite a 31-fold increase in the FDR denominator — the lowest blind recovery of any disease, consistent with ASD's weak transcriptomic signal.

**T2D:** From 26,800 expressed genes with zero prior hypothesis, the blind run returned IAPP (islet amyloid polypeptide) — the pathological hallmark of type 2 diabetes — as the number one ranked signal. The engine also discovered PCSK1, ERO1B, MAFB, and CHGB, reconstructing the complete beta cell failure cascade from raw expression data.

**Breast Cancer:** The engine independently separated estrogen receptor biology (ESR1, Cluster 18, random-effects p = 0.019), HER2 biology (ERBB2, Cluster 28, p = 0.007), and proliferation biology (TOP2A, AURKA; Cluster 2, p = 0.009) into distinct clusters, recapitulating the clinical molecular subtype framework from raw expression data without prior knowledge of subtype structure. A hypothesis-free blind run using all expressed genes recovered 151 of 152 curated core genes (99.3%); only SIAH2 was not recovered.

**IPF:** The most rigorous validation. 190 core genes were tested in GSE47460 — a 122-patient IPF cohort the engine never saw. 132 of 153 testable genes (86.3%) were significantly differentially expressed in the same direction (96.7% concordance). Leave-one-out stability testing identified 52 iron-clad genes surviving every dataset combination and cold replication. Novel candidates FAM107A and WDR49 were identified through co-clustering with established IPF genes.

### 3.3 Independent Third-Party Validation

Two diseases were validated entirely by independent AI agents with no author involvement in any aspect of the validation.

#### 3.3.1 Psoriasis (Gemini CLI Agent)

An independent Gemini CLI agent selected Psoriasis, curated 96 seed genes from literature knowledge, independently chose 5 GPL570 microarray datasets, built a YAML config from scratch, and ran the pipeline.

Results: 60 study genes (62.5% Phase 1 yield) → 9 clusters → 50 core genes → 50 survived replication (100%) → 28 meta-significant (random effects). QC: 4 passed, 0 critical.

The engine recovered textbook psoriasis biology: S100A12 (p ≈ 6×10⁻³¹⁰), S100A9, S100A7 (psoriasin), DEFB4B, PI3 (Elafin), IL36G, CCL20, CXCL8, KRT16, KRT17, and BCL2 (downregulated, p ≈ 7×10⁻⁵⁰). The agent's uncoached verdict: "The engine generalizes. It worked flawlessly on a new disease without requiring any internal modifications."

#### 3.3.2 Colorectal Cancer (Claude Code Agent)

An independent Claude Code agent selected Colorectal Cancer, curated 515 seed genes from COSMIC, Open Targets, KEGG, and published CRC transcriptomic signatures, independently chose 6 GPL570 datasets spanning 5 countries (781 tumors, 144 normals), and ran the pipeline.

Results: 331 study genes → 43 clusters → 263 core genes → 244 survived replication (92.8%) → 218 meta-significant. QC: 4 passed, 0 critical.

Pathway recovery was near-complete: 100% coverage on cell cycle, p53, TGF-beta, PI3K/AKT, RAS/MAPK, immune/inflammation, ECM/invasion, epigenetics, drug targets, angiogenesis, and mismatch repair pathways. The Wnt/APC pathway was recovered at 96%, with APC itself identified as downregulated (p = 1.1×10⁻⁴⁵).

Critically, the agent validated what the engine correctly did NOT find. Mutation-driven genes (BRAF, PIK3CA, NRAS, ERBB2) failed Phase 1 because their risk is conferred through protein-level variants, not transcript-level changes. The engine correctly cannot detect them because they are not transcriptomically altered. This represents appropriate calibration, not a limitation.

The agent's uncoached verdict: "The Riker Engine is a powerful tool for meta-analytical gene discovery. It is stable, generalizable, and mathematically rigorous. I recommend we adopt it for our lab's transcriptomic screening."

### 3.4 Cold-Start Replication Testing

Two rounds of cold-start replication were performed on a separate Raspberry Pi 5 with a clean filesystem. Independent agents were given only the GitHub URL and a skeptical prompt.

IPF replication results across runs:

| Metric | Author | Gemini Run 1 | Gemini Run 2 | Claude Code |
|--------|--------|-------------|-------------|-------------|
| Phase 1 Study Genes | 241 | 241 | 241 | 241 |
| Phase 4 Core Genes | 190 | 189 | 189 | 189 |
| Core Gene Overlap | — | 99.47% | 99.47% | 99.47% |
| Missing Gene | — | ID1 | ID1 | ID1 |
| Cold Replication | 86.3% | 86.2% | 86.2% | 86.3% |

Deterministic phases (Phase 1, Phase 5 elimination, Phase 6 meta-analysis) reproduced bit-for-bit across all runs, including effect sizes identical to 10 decimal places. The sole source of variance — one gene (ID1) falling to noise in Phase 3 consensus clustering — represents 0.5% of the core gene set, well within the documented ±5–8% stochastic tolerance.

### 3.5 WGCNA Benchmark Comparison

A head-to-head comparison was conducted on identical ASD datasets and hardware.

| Metric | Riker Engine | WGCNA |
|--------|-------------|-------|
| Runtime | ~8 min (7 datasets) | ~3 h 50 min (3 datasets) |
| Output | 35 core genes | 1,427–9,995 genes/module |
| Cross-dataset validation | Built-in | Manual |
| Memory | Full 18K gene set | OOM on 18K; required filtering |
| Gene overlap | — | 34/35 core genes in WGCNA modules |

WGCNA is not producing incorrect results. Its modules contain the correct genes along with hundreds to thousands of additional genes. The Riker Engine automates what currently requires weeks of manual post-processing: running WGCNA per dataset, comparing module preservation, filtering for replication, and performing meta-analysis.

---

## 4. Discussion

### 4.1 What the Tool Does Well

The eight-disease validation demonstrates that a condition-agnostic pipeline can recover established disease biology across unrelated conditions and tissue types without per-disease parameter tuning. The pipeline's behavior is appropriately calibrated: weak signals produce small, precise core gene sets, strong signals produce comprehensive sets.

The independent validations are particularly compelling. Neither agent had access to expected results, author guidance, or domain-specific parameter tuning. Both produced textbook-perfect biology for their chosen diseases, confirming that the methodology generalizes beyond the author's validation set.

### 4.2 What the Tool Does Not Do

The engine identifies consistently differentially expressed genes, not causal drivers. It cannot distinguish drivers from passengers. It operates at bulk tissue level and cannot resolve cell-type-specific signals without pseudo-bulked snRNA-seq input. It does not integrate survival data, drug response data, or mutation co-occurrence. These are limitations inherent to bulk transcriptomic analysis, shared with all tools in this class.

The engine is a validation and filtering tool, not a de novo discovery engine. The seeded mode requires prior candidate genes. Novel candidates emerge through co-clustering, but the primary use case is determining which candidates from a GWAS, screen, or literature review hold up across independent expression datasets.

### 4.3 Limitations

**Welch's t-test for RNA-seq data.** Variance-stabilizing methods such as limma-voom and DESeq2 model count-based variance structures more appropriately. Welch's t-test on log2-transformed data is conservative rather than liberal, and cross-dataset replication mitigates single-dataset statistical weaknesses.

**Microarray optimization.** The current pipeline is optimized for GEO microarray series matrices. RNA-seq datasets require manual reconstruction from supplementary files. Support exists but requires additional setup effort.

**Memory scaling.** The consensus clustering matrix scales as O(n²). On 8 GB hardware, blind runs are limited to approximately 9,000–10,000 study genes.

**Single-developer risk.** The engine was developed by a single independent researcher. The codebase has 300+ passing tests and has been independently replicated, but institutional adoption would benefit from independent code review. JOSS submission is planned for September 2026.

**AGPL-3.0 license.** Derivative works deployed as services must be open-sourced. Labs with industry collaborations should consult their technology transfer office.

### 4.4 Broader Context

The engine fills a specific gap in the current transcriptomic analysis landscape. The alternative to this tool is weeks of manual GEO mining with ad hoc scripts. The engine performs this analysis in under 5 minutes with reproducible statistics, pre-specified replication, and transparent quality control. For rare and understudied diseases where small expression datasets exist in public repositories but have not been systematically integrated, the engine provides a rapid, rigorous path from candidate genes to replicated modules.

---

## 5. Conclusion

The Riker Engine demonstrates that a condition-agnostic, six-phase progressive filtering pipeline can identify replicated gene modules across multiple independent transcriptomic datasets. Validated on eight diseases spanning six tissue types — including two validated entirely by independent third parties — the engine recovers established disease biology, identifies known drug targets, and produces minimal core gene sets suitable for experimental follow-up. The consensus clustering framework converts parameter-sensitive UMAP embeddings into parameter-robust co-association matrices, and the full-seed-set FDR correction provides appropriate statistical control without sacrificing sensitivity.

The engine is released under the AGPL-3.0 license with 300+ passing tests at github.com/RaySigmon/Riker_Engine.

---

## 6. Data and Code Availability

All source code, configurations, seed gene lists, and validation results (including independent third-party validation outputs) are available at https://github.com/RaySigmon/Riker_Engine. All GEO datasets used in this work are publicly accessible through the NCBI Gene Expression Omnibus.

## 7. Acknowledgments

This work was conducted independently without institutional affiliation or external funding. The author acknowledges the Gene Expression Omnibus, SFARI Gene database, KEGG, Reactome, MSigDB, Open Targets, and COSMIC. AI-assisted coding tools (Anthropic Claude, Claude Code CLI "Kai") were used for pipeline development; their contributions are documented in the project repository. Independent validation was performed by AI agents (Google Gemini CLI, Anthropic Claude Code CLI) operating without author guidance.

Named after his son Riker, who inspired a father to search across dimensions for answers.

## References

Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate. JRSS:B, 57(1), 289–300.

Campello RJGB, Moulavi D, Sander J (2013). Density-based clustering based on hierarchical density estimates. PAKDD 2013.

DerSimonian R, Laird N (1986). Meta-analysis in clinical trials. Controlled Clinical Trials, 7(3), 177–188.

Fabregat A et al. (2018). The Reactome pathway knowledgebase. Nucleic Acids Research, 46(D1), D649–D655.

Kanehisa M, Goto S (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Research, 28(1), 27–30.

Langfelder P, Horvath S (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9, 559.

Liberzon A et al. (2015). The MSigDB hallmark gene set collection. Cell Systems, 1(6), 417–425.

McInnes L, Healy J, Melville J (2018). UMAP: Uniform Manifold Approximation and Projection. arXiv:1802.03426.

Ritchie ME et al. (2015). limma powers differential expression analyses. Nucleic Acids Research, 43(7), e47.

Welch BL (1947). The generalization of Student's problem. Biometrika, 34(1–2), 28–35.
