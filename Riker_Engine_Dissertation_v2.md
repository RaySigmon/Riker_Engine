# THE RIKER ENGINE
## A Condition-Agnostic Framework for Identifying Coordinately Dysregulated Gene Modules from Public Transcriptomic Data

### Comprehensive Technical Dissertation

**Ray Sigmon**
Independent Computational Researcher, Northwest Arkansas, USA

**Version 0.3.2 — April 2026**

Repository: github.com/RaySigmon/Riker_Engine | License: AGPL-3.0

*Named after his son Riker, who inspired a father to search across dimensions for answers.*

---

## TABLE OF CONTENTS

1. Origin Story: From a Father's Question to a Scientific Pipeline
2. What the Riker Engine Is (and Is Not)
3. Engine Architecture: The Six-Phase Pipeline
4. Phase 1: Cross-Dataset Differential Expression
5. Phase 2: Pathway Mapping and Feature Matrix Construction
6. Phase 3: Consensus Clustering
7. Phase 4: Robustness Testing and Core Gene Identification
8. Phase 5: Independent Replication
9. Phase 6: Effect Size Meta-Analysis
10. Quality Control Framework
11. Validation History: Six Author-Validated Diseases
12. Independent Third-Party Validation
13. Cold-Start Replication Testing
14. WGCNA Benchmark Comparison
15. Cross-Disease Patterns and Calibration
16. The Rare Disease Opportunity
17. snRNA-seq Pseudo-Bulking Module
18. Critical Review and Honest Assessment
19. Known Limitations and Honest Boundaries
20. Technical Specifications and Infrastructure
21. Development Methodology: AI-Assisted Engineering
22. Licensing, Sustainability, and Open-Core Model
23. Future Work and Roadmap
24. Conclusion

Appendix A: Validation Summary Table (8 Diseases)
Appendix B: Cold-Start Replication Data
Appendix C: Independent Critical Review (9 Points) with Author Responses
Appendix D: Technical Specifications
Appendix E: References

---

## 1. Origin Story: From a Father's Question to a Scientific Pipeline

The Riker Engine began not as a software project but as a question from a father.

In early 2026, Ray Sigmon — a personal trainer and self-taught programmer in Northwest Arkansas — began searching for computational approaches to understanding autism spectrum disorder. His son Riker had inspired him to look for patterns in the genetic data that professional researchers might be missing, not because he knew more biology than they did, but because he was willing to ask a different kind of question: across all these candidate genes, do any of them move together?

The SFARI Gene database catalogs over 1,200 genes implicated in ASD. Most studies examine these genes one at a time. Sigmon's insight was to treat the entire gene set as a feature space and ask whether agnostic clustering — without any prior hypothesis about which genes should group together — could identify coordinately dysregulated modules that survived independent replication.

### 1.1 Project Riker: The ASD Proof of Concept

The first implementation — Project Riker — cross-referenced 1,267 SFARI genes against three GEO transcriptomic datasets. The progressive narrowing identified 243 study genes, 20 consensus clusters, one uniformly downregulated 11-gene module (permutation p = 0.0026), a 4-gene FDR-surviving core (ATP10A, ATP2B2, KANK1, SEZ6L2), and ultimately 2 meta-analytically confirmed genes showing brain-specific downregulation across four independent postmortem cohorts (total N = 228).

ATP2B2 (plasma membrane calcium ATPase) and SEZ6L2 (complement-regulating synaptic protein) survived every stage of progressive filtering. Their co-identification through agnostic clustering — discovered without prior hypothesis about either gene — suggested a convergent mechanism: neurons with simultaneously impaired calcium extrusion and reduced complement protection may be both hyperexcitable and vulnerable to synaptic loss.

Blood replication failed, confirming brain specificity. KANK1 was eliminated for significant upregulation in brain replication — the opposite direction from discovery. These honest negatives strengthened the positives: the pipeline was not simply confirming everything, but actively eliminating inconsistent signals.

The preprint was submitted to bioRxiv three times and rejected each time for lack of institutional affiliation — a structural barrier that ultimately led to the JOSS submission pathway.

### 1.2 The Alzheimer's Stress Test

The pipeline was subsequently applied to Alzheimer's disease using 800 seed genes. This stress test produced a critical finding: FDR correction over the full seed set yielded zero survivors at q < 0.10, despite 119 validated core genes showing 87.5% replication. The pipeline correctly identified that bulk tissue data lacked sufficient resolution for FDR-surviving discoveries in AD — the high replication rate was consistent with cell-type composition artifacts from pervasive neuronal loss, not genuine regulatory changes.

This null result was the pipeline working correctly. It demonstrated the difference between "consistent signal" and "meaningful signal," and motivated the development of the snRNA-seq pseudo-bulking module.

### 1.3 From Project to Engine

The critical realization was that the methodology — not the ASD-specific application — was the primary contribution. The same pipeline that identified ATP2B2 and SEZ6L2 in autism could be applied to any disease with candidate genes and public expression data. This insight drove the generalization from "Project Riker" (an ASD study) to "Riker Engine" (a condition-agnostic tool).

The engine was rebuilt as a pip-installable Python package with a CLI interface, YAML configuration, and automated data download. The architecture was designed so that switching diseases required only a new configuration file — zero code modifications.

### 1.4 Hardware and Infrastructure

The entire project was developed on a Raspberry Pi 5 (8GB RAM) named "Ghost," using Claude Code CLI ("Kai") as the primary autonomous builder. Claude.ai served as an independent auditor and technical advisor, deliberately separated from Kai's project-level context to prevent confirmation bias. This dual-agent architecture — one agent building, one agent reviewing — became a defining feature of the development methodology.

---

## 2. What the Riker Engine Is (and Is Not)

### What It IS

The Riker Engine is a cross-dataset validation and filtering pipeline for disease gene candidates. Given a set of candidate genes and multiple independent expression datasets, it identifies which candidates show consistent, replicated differential expression across studies — with proper statistical controls, pre-specification, and independent replication.

**Use it when you have:**
- A list of candidate genes from a GWAS, screen, collaborator's paper, or pathway database
- A question: "Which of these genes hold up across independent datasets?"
- Access to 3+ GEO expression datasets for your disease

**What it produces:**
- A filtered core gene set that survives consensus clustering and cross-dataset replication
- Meta-analysis effect sizes and p-values for each surviving gene
- Novel candidates that co-cluster with your seeds but were not in your original list

### What It Is NOT

- A black-box tool that produces "answers." Every output requires human interpretation.
- A de novo discovery engine. The seeded mode requires prior candidate genes. Blind mode exists but is experimental and memory-intensive.
- A replacement for domain expertise. The pipeline identifies statistical patterns; biological interpretation requires a biologist.
- A causal inference tool. It identifies consistently differentially expressed genes, not drivers versus passengers.
- A cell-type deconvolution tool. It operates at bulk tissue level (unless fed pseudo-bulked snRNA-seq data).
- Guaranteed to find something. The AD stress test demonstrated correct null reporting.

The alternative to this tool is weeks of manual GEO mining with ad hoc scripts. The engine performs this analysis in under 5 minutes with reproducible statistics.

---

## 3. Engine Architecture: The Six-Phase Pipeline

The Riker Engine consists of a data ingestion layer feeding into a six-phase analytical pipeline. The ingestion layer is the only component that changes between applications; Phases 1–6 are identical regardless of disease, tissue, or data type.

| Phase | Function | Input | Output |
|-------|----------|-------|--------|
| 1 | Cross-Dataset DE Filtering | Seeds + expression matrices | Study gene set |
| 2 | Feature Matrix Construction | Study genes + pathway DBs | Normalized feature matrix |
| 3 | Consensus Clustering | Feature matrix | Cluster assignments + consensus matrix |
| 4 | Robustness Testing | Clusters + DE stats | Core genes + dissolution narratives |
| 5 | Independent Replication | Core genes + held-out datasets | Validated genes + elimination log |
| 6 | Meta-Analysis | All datasets | Forest plots + pooled estimates |

Phases 1–4 constitute the discovery framework. Phases 5–6 constitute the validation framework. Both are required for any defensible finding.

---

## 4. Phase 1: Cross-Dataset Differential Expression

For each seed gene in each dataset, the engine computes log2FC (mean cases minus mean controls on the log2 scale) and Welch's t-test p-value using the exact t-distribution. Genes reaching nominal significance (p < 0.05) in at least two datasets are retained as the study gene set.

This threshold is intentionally lenient; stricter thresholds are applied in Phase 4. The engine logs which genes were detectable on each platform, providing per-platform coverage statistics.

### 4.1 Why Welch's t-test

Three design priorities drive this choice. First, uniformity: Welch's t-test operates identically on log2-transformed microarray and RNA-seq data, enabling cross-platform analysis without method switching. Second, conservatism: on log2-transformed data, Welch's t-test is less powerful than limma or DESeq2, producing fewer false positives at the cost of more false negatives. Since Phase 1 is a dimensionality reduction step, conservatism is acceptable. Third, independence: unlike limma, which models variance across genes, Welch's t-test treats each gene independently.

### 4.2 Expression Scale Handling

Some platforms produce log-ratio values centered near zero rather than log2-intensity values. The engine detects this by checking if the expression range maximum is below 5.0, and flags these datasets for special handling. Standard errors are recovered from the native scale; z-score distributions are verified across datasets to confirm comparability before pooling.

### 4.3 Log2 Transformation Detection

Microarray data may or may not be pre-log2-transformed. The engine automatically detects transformation status by checking the median of all expression values. If median >= 20, data is raw intensity and requires log2(x+1) transformation. This check was added after detecting biologically impossible fold-change values (log2FC = 119) caused by computing differences on raw intensity data.

---

## 5. Phase 2: Pathway Mapping and Feature Matrix Construction

Study genes are mapped to biological pathways from three sources: KEGG (approximately 350 human pathways), Reactome (approximately 2,200 human pathways), and MSigDB Hallmark gene sets (50 curated biological state gene sets). Pathways are filtered by minimum 3 study gene members, maximum 500 total genes, and minimum 2% study gene representation.

The feature vector per gene combines binary pathway memberships with five expression statistics: average log2FC across datasets, negative log10 of minimum p-value, number of significant datasets, directional consistency score, and confidence tier score. All features are min-max normalized to [0,1] before clustering.

### 5.1 Anti-Circularity Rule

Do NOT use pre-assigned biological category labels (e.g., "chromatin remodeling") as clustering features. Use only individual pathway IDs as binary features. Using categories as both features and enrichment targets is tautological and will produce inflated statistics. This error was caught and corrected during Project Riker development.

### 5.2 Expression-Only Mode

All published validation results were produced without pathway data configured, using expression-based features only (5 dimensions). This is the validated default. The log message "No pathway data configured. Using expression features only." is informational, not an error. Pathway integration via MSigDB is available as an optional enhancement.

---

## 6. Phase 3: Consensus Clustering

UMAP (n_components=2, min_dist=0.0, metric=euclidean) followed by HDBSCAN (min_cluster_size=5, min_samples=3) is applied across 15 configurations: 3 n_neighbors values (10, 15, 30) × 5 random seeds (42, 123, 456, 789, 1024).

A co-association consensus matrix records pairwise co-clustering frequencies. Final cluster labels are derived by applying HDBSCAN with metric='precomputed' to the distance matrix (1 − consensus/15).

### 6.1 Why Consensus

Single runs of UMAP/HDBSCAN are sensitive to random seed and hyperparameter choices. Consensus clustering identifies which gene groupings are stable across multiple configurations. The three n_neighbors values capture local (10), medium (15), and broad (30) neighborhood structure, while the five random seeds capture UMAP's inherent stochasticity.

### 6.2 Stochastic Variance

Because UMAP initialization and HDBSCAN tie-breaking involve random elements, reproduced gene counts may vary by approximately ±5–8% from published numbers. Different numpy/scipy versions may also contribute to minor variance. The core gene lists should overlap substantially (>90%) even when counts differ slightly. In practice, observed variance has been limited to 0.5% (1 gene out of 190 in IPF replication).

---

## 7. Phase 4: Robustness Testing and Core Gene Identification

### 7.1 Bonferroni Permutation Testing

Permutation p-values are evaluated against a Bonferroni-corrected threshold (0.05 / number of clusters). 10,000 permutations are performed.

### 7.2 Progressive Threshold Sensitivity Analysis

For each cluster of interest, member genes are tested at four progressively stricter thresholds:
- Level 1: p < 0.05 in 2+ datasets (baseline)
- Level 2: p < 0.01 in 2+ datasets
- Level 3: FDR q < 0.10 in 2+ datasets (BH correction within each dataset over the FULL seed gene set)
- Level 4: FDR q < 0.05 in 2+ datasets

### 7.3 Full-Seed-Set FDR

This is a critical design decision. Most tools compute FDR only over genes that passed an initial filter, excluding non-significant genes from the denominator. This inflates significance. The Riker Engine computes Benjamini-Hochberg FDR over the entire seed set. In stress testing, study-set-only FDR produced 420 genes at q < 0.10; full-seed-set FDR correctly returned zero.

### 7.4 Leave-One-Dataset-Out Stability

Each discovery dataset is removed in turn and the retention criterion is recomputed. A cluster is LOO-stable if 80% or more of its genes survive removal of any single dataset.

### 7.5 Core Gene Identification and Pre-Specification Lock

Genes surviving at least Level 2 with 3 or more cluster-mates are designated core genes. The core gene list is locked before any replication data is accessed. This pre-specification protocol prevents post-hoc selection bias, mirroring pre-registration logic from clinical trials.

---

## 8. Phase 5: Independent Replication

### 8.1 Protocol

The core gene list is finalized and documented BEFORE any replication dataset is examined. No genes may be added or removed based on replication results. No exploratory analysis is performed on replication datasets.

### 8.2 Elimination Protocol

Any gene significantly (p < 0.05) in the OPPOSITE direction in a same-tissue replication dataset is eliminated. Non-significant opposite effects are noted but not automatically disqualifying. Cross-tissue non-replication (e.g., brain genes failing in blood) is informative rather than disqualifying.

This rule was empirically validated when KANK1 was correctly eliminated in ASD: it passed all discovery filters but showed significant upregulation in brain replication, opposite to its discovery direction. Similarly, CTNNB1 was correctly eliminated in CRC — its protein is constitutively active via APC loss, but mRNA levels are inconsistent across datasets.

---

## 9. Phase 6: Effect Size Meta-Analysis

Inverse-variance weighted meta-analysis pools effect sizes across all datasets. Standard errors are recovered from log2FC and p-values via the t-distribution. Both fixed-effects and random-effects models are computed.

REML (restricted maximum likelihood) is the primary estimator for the between-study variance component, with DerSimonian-Laird as fallback when REML does not converge. REML was selected because DerSimonian-Laird underestimates tau-squared when the number of studies is small (fewer than 10), a condition typical of transcriptomic meta-analysis.

Heterogeneity is quantified using Cochran's Q, I-squared, and tau-squared. The random-effects estimate is the primary result.

---

## 10. Quality Control Framework

These checks were developed during Project Riker after catching errors that would have produced false findings.

| Check | What It Catches | When Applied |
|-------|----------------|-------------|
| Log2 detection | Impossible fold changes from raw intensity | Phase 1, every dataset |
| Circularity audit | Category labels used as both features and enrichment targets | Phase 2, before clustering |
| Permutation test | Clusters expected by chance | Phase 3, every cluster |
| FDR scope enforcement | Inflated survival from reduced denominator | Phase 4, sensitivity analysis |
| Sensitivity analysis | Findings dependent on lenient thresholds | Phase 4, every cluster |
| LOO stability | Single-dataset-dependent clusters | Phase 4, every cluster |
| Pre-specification lock | Post-hoc gene selection bias | Phase 5, before data access |
| Elimination protocol | Retaining genes with opposite-direction replication | Phase 5, after each dataset |
| Z-score verification | Cross-platform scale incompatibility | Phase 6, before pooling |

---

## 11. Validation History: Six Author-Validated Diseases

### 11.1 Autism Spectrum Disorder

**Tissue:** Brain cortex | **Seeds:** 1,267 (SFARI) | **Datasets:** 7 | **Phase 1 yield:** 11.1%

The weakest transcriptomic signal. 35 core genes identified through consensus clustering. 15 eliminated for significant opposite-direction effects in blood replication datasets — confirming brain specificity and demonstrating the pipeline's ability to distinguish tissue-specific signals from noise. 20 survived replication, including ATP2B2, SEZ6L2, SNAP25, GABRG2, TBR1. 9 meta-significant under random effects. Blind recovery: 27/35 (77.1%).

Key finding: ATP2B2 and SEZ6L2 showed significant brain-specific downregulation across four independent postmortem cohorts (random-effects p = 0.037 and 0.046). The 15 eliminated genes (including EFR3A, FBXO33, PPP3CA, SPTAN1, APBA2, ATP1A1, ITPR1, NAV2, ZFHX3, and others) all showed significant upregulation in blood datasets — the opposite direction from their brain discovery signal. This is the engine correctly identifying brain-specific biology that does not generalize to peripheral tissue.

### 11.2 Type 2 Diabetes

**Tissue:** Pancreatic islets | **Seeds:** 443 (Open Targets) | **Datasets:** 4 | **Phase 1 yield:** 12.6%

8 core genes forming a single coherent module of insulin secretory machinery: SLC2A2 (GLUT2), ABCC8 (sulfonylurea drug target), CACNA1D, CASR. All survived replication.

Blind run highlight: From 26,800 expressed genes with zero prior hypothesis, IAPP (islet amyloid polypeptide) emerged as the number one ranked signal. IAPP aggregation is the pathological hallmark of T2D.

### 11.3 Inflammatory Bowel Disease

**Tissue:** Intestinal mucosa | **Seeds:** 762 | **Datasets:** 6 (629 samples) | **Phase 1 yield:** 53.5%

The strongest transcriptomic signal. 304 core genes including drug targets NOD2, JAK2 (tofacitinib target), STAT3, TNFAIP3, TLR4, HNF4A, SMAD7. Two genes eliminated (99.3% survival). Blind recovery: 97.7%.

### 11.4 Alzheimer's Disease

**Tissue:** Brain cortex | **Seeds:** 801 | **Datasets:** 5 | **Phase 1 yield:** 54.7%

394 core genes including TREM2, APOE, APP, MAPT, CLU, BIN1, CD33. PSEN1 absent from core — confers risk through protein-level variants, not transcript changes. Highest elimination rate (13.5%, 53 genes). Blind recovery: 387/394 core genes (98.2%) recovered when run with all 19,296 protein-coding genes as seeds (RunPod, 64 GB RAM). Only 7 genes not recovered.

### 11.5 Breast Cancer

**Tissue:** Breast tumor | **Seeds:** 653 | **Datasets:** 5 | **Phase 1 yield:** 34.8%

152 core genes across 28 clusters. The engine independently separated estrogen receptor biology (ESR1, Cluster 18, random-effects p = 0.019), HER2 biology (ERBB2, Cluster 28, p = 0.007), and proliferation biology (TOP2A, AURKA; Cluster 2, p = 0.009) into distinct modules — recapitulating the clinical molecular subtype framework from raw expression data. Blind recovery: 151/152 (99.3%) — only SIAH2 not recovered.

### 11.6 Idiopathic Pulmonary Fibrosis

**Tissue:** Lung | **Seeds:** 354 | **Datasets:** 5 | **Phase 1 yield:** 68.1%

The most rigorous validation. 190 core genes tested in GSE47460 — a 122-patient cohort the engine never saw. 132 of 153 testable genes (86.3%) were significantly differentially expressed in the same direction (96.7% concordance). Leave-one-out stability testing identified 52 iron-clad genes surviving every dataset combination AND cold replication.

Novel candidates FAM107A and WDR49 were identified through co-clustering with established IPF genes — these were not in the original seed list.

---

## 12. Independent Third-Party Validation

Two diseases were validated entirely by independent AI agents with no author involvement in disease selection, seed gene curation, dataset choice, or configuration.

### 12.1 Psoriasis (Gemini CLI Agent — "Jim")

**Disease selection:** Agent's independent choice
**Seeds:** 96 genes curated by the agent from literature knowledge
**Datasets:** 5 GPL570 microarray datasets independently selected
**Results:** 60 study genes (62.5% Phase 1 yield) → 50 core genes → 100% replication survival → 28 meta-significant (random effects)
**QC:** 4 passed, 0 warnings, 0 critical

The engine recovered textbook psoriasis biology:
- **Antimicrobial peptides:** S100A12 (p ≈ 6×10⁻³¹⁰), S100A9, S100A7 (psoriasin), DEFB4B, PI3 (Elafin)
- **Inflammatory mediators:** IL36G, IL36RN, CCL20, CXCL8, CXCL1, TNF, IL1B, IL12B
- **Keratinocyte markers:** KRT16, KRT17, KRT6A, IVL, TGM1
- **Apoptosis dysregulation:** BCL2 (downregulated, p ≈ 7×10⁻⁵⁰)

Agent's uncoached verdict: "I recommend we adopt it for our lab's transcriptomic screening."

### 12.2 Colorectal Cancer (Claude Code Agent)

**Disease selection:** Agent's independent choice
**Seeds:** 515 genes curated from COSMIC, Open Targets, KEGG, published CRC signatures
**Datasets:** 6 GPL570 datasets spanning 5 countries (781 tumors, 144 normals)
**Results:** 263 core genes → 244 survived replication (92.8%) → 218 meta-significant
**QC:** 4 passed, 0 warnings, 0 critical

Pathway coverage was near-complete:

| Pathway | Coverage |
|---------|----------|
| Cell cycle / proliferation | 28/28 (100%) |
| p53 / DNA damage repair | 25/25 (100%) |
| TGF-beta / SMAD | 10/10 (100%) |
| PI3K/AKT/mTOR | 9/9 (100%) |
| RAS/MAPK/EGFR | 10/10 (100%) |
| Immune / inflammation | 32/32 (100%) |
| ECM / invasion / metastasis | 22/22 (100%) |
| Epigenetics / chromatin | 11/11 (100%) |
| Drug targets | 7/7 (100%) |
| Wnt/APC pathway | 22/23 (96%) |
| Angiogenesis | 6/6 (100%) |
| CRC tissue markers | 9/11 (82%) |
| Mismatch repair (MMR) | 3/3 (100%) |

Top hits included APC (downregulated, p = 1.1×10⁻⁴⁵ — the initiating event in ~80% of CRC), MMP7 (log2FC = +4.51), AURKA (active drug target), WNT2, CXCL12, BCL2, and DNMT1.

The agent also validated what the engine correctly did NOT find: BRAF, PIK3CA, NRAS, and ERBB2 failed Phase 1 because they are mutation-driven, not expression-driven. The engine correctly cannot detect them because they are not transcriptomically altered. This represents appropriate calibration.

Agent's uncoached verdict: "The Riker Engine is a powerful tool for meta-analytical gene discovery. It is stable, generalizable, and mathematically rigorous."

---

## 13. Cold-Start Replication Testing

### 13.1 Methodology

Two rounds of cold-start replication were performed on a separate Raspberry Pi 5 ("Ghost") with a clean filesystem. Independent AI agents were given only the GitHub URL and a skeptical PI-style prompt: "I'm skeptical but willing to look. [...] I want to know if this is reproducible science or a black box."

The agents were NOT told which disease to replicate, which datasets to use, or what results to expect.

### 13.2 Isolation Protocol

Before each test: all prior clones, output directories, config files, virtual environments, and leftover artifacts were deleted. The agent confirmed the correct commit hash (e6bf09e) before beginning.

### 13.3 Results

IPF replication was performed by three independent agents (two Gemini runs, one Claude Code run):

| Metric | Author | Gemini R1 | Gemini R2 | Claude Code |
|--------|--------|-----------|-----------|-------------|
| Phase 1 Study Genes | 241 | 241 | 241 | 241 |
| Phase 4 Core Genes | 190 | 189 | 189 | 189 |
| Core Gene Overlap | — | 99.47% | 99.47% | 99.47% |
| Missing Gene | — | ID1 | ID1 | ID1 |
| Ph5 Eliminated | 20 | 20 | 20 | 20 |
| Cold Replication | 86.3% | 86.2% | 86.2% | 86.3% |

### 13.4 Determinism Analysis

Deterministic phases reproduced bit-for-bit across all runs:
- Phase 1: All 241 genes, all effect sizes, all p-values — identical to 10 decimal places
- Phase 5: Same 20 genes eliminated for the same documented reasons
- Phase 6: All overlapping meta-analysis estimates identical to 10 decimal places, 100% direction agreement

The sole source of variance: one gene (ID1) falling to noise in Phase 3 consensus clustering. This represents 0.5% of the core gene set — well within the documented ±5–8% stochastic tolerance.

### 13.5 Bug Discovery and Resolution

The cold-start testing process identified and resolved 12 issues across two fix rounds:

**Round 1 (Fixes 1–6):**
1. GPL URL bucket calculation bug for 3-digit platform IDs
2. IPF missing from download script (confirmed already fixed in v0.3.1)
3. GPL96 probe mapping failure (resolved by Fix 1)
4. seed_genes: null crash (clear error message added)
5. Reproducibility variance note added to README
6. New disease configuration guide (docs/NEW_DISEASE_GUIDE.md)

**Round 2 (Fixes 7–12):**
7. GSE47460 held-out dataset added to IPF download script
8. ipf_curated.yaml moved to configs/examples/
9. Install docs changed from editable (-e) to non-editable mode
10. gene_column and probe_column exposed through YAML config
11. NEW_DISEASE_GUIDE.md written (297 lines)
12. Hardcoded hgnc_path confirmed already fixed

---

## 14. WGCNA Benchmark Comparison

A head-to-head comparison was conducted on identical ASD brain cortex datasets using the same hardware (Raspberry Pi 5, 8 GB RAM).

| Metric | Riker Engine | WGCNA |
|--------|-------------|-------|
| Runtime | ~8 min (7 datasets, all 6 phases) | ~3 h 50 min (3 datasets only) |
| Output | 35 core genes | 1,427–9,995 genes per module |
| Cross-dataset validation | Built-in (Phases 4–5) | None (manual) |
| Memory | Full 18K gene set | OOM on 18K; required filtering to 10K |
| Scale-free fit | N/A (not required) | Failed on GSE28475 (R² = 0.35) |
| Gene overlap | — | 34/35 Riker core genes in WGCNA modules |

The gene overlap result is important: WGCNA is not producing incorrect results. Its modules contain the correct genes along with hundreds to thousands of additional genes that the Riker Engine's progressive filtering eliminates. The engine automates what currently requires a skilled bioinformatician weeks of manual work.

---

## 15. Cross-Disease Patterns and Calibration

Five patterns emerged consistently across all validated diseases:

**1. Phase 1 blind recovery.** Every curated core gene passed Phase 1 when tested against the full transcriptome in blind runs (100% Phase 1 recovery), confirming genuine signals regardless of seed list composition.

**2. Signal-proportional calibration.** Weak signals (ASD 11%, T2D 13% Phase 1 yield) produce small, precise core gene sets (8–35 genes). Strong signals (IBD 54%, AD 55%) produce comprehensive sets (304–394 genes). The engine calibrates to the biology; it does not force a result.

**3. GWAS genes behave differently.** Top GWAS hits such as TCF7L2 (T2D), PSEN1 (AD), IL23R (IBD), and PIK3CA (breast cancer) do not pass Phase 1, because genetic risk variants do not necessarily alter transcript abundance. The engine identifies functional consequences at the expression level, not genetic causes.

**4. Known biology first, then extensions.** In every disease, the engine first recovers established biology and then discovers additional genes within the same functional modules.

**5. Drug targets emerge naturally.** ABCC8 (sulfonylureas) in T2D, JAK2 (tofacitinib) in IBD, TOP2A (doxorubicin), ERBB2 (trastuzumab), ESR1 (tamoxifen) in breast cancer, AURKA (alisertib) in CRC — all identified without being directed.

---

## 16. The Rare Disease Opportunity

The validation demonstrates that the engine reliably finds known biology with as few as 8 samples per group (T2D). Rare disease datasets with 10–15 samples per condition are therefore viable inputs.

Target conditions occupy a specific sweet spot: diseases uncommon enough to be under-researched but common enough to have generated 2–3 GEO expression datasets. Idiopathic pulmonary fibrosis (now validated), scleroderma, primary biliary cholangitis, Huntington's disease, Duchenne muscular dystrophy, amyotrophic lateral sclerosis, certain epilepsy syndromes, and pediatric cancers all have small expression studies in public databases that have not been systematically integrated.

Three converging trends make this application timely. The FDA's Plausible Mechanism Framework requires validated targets for individualized therapies. Gene therapy development costs are decreasing while targetable conditions increase. Public expression data in GEO grows continuously.

---

## 17. snRNA-seq Pseudo-Bulking Module

### 17.1 The Problem

Bulk tissue transcriptomics measures an average signal across all cell types. In neurodegenerative diseases, cell death creates pervasive correlated expression changes indistinguishable from true regulatory changes without deconvolution.

### 17.2 The Solution

Single-nucleus RNA-seq provides per-cell expression profiles. Pseudo-bulking aggregates these into one expression value per gene per cell type per donor, producing a format identical to bulk tissue data but with cell-type resolution. The engine then runs the standard Phase 1–6 pipeline separately for each cell type.

### 17.3 Implementation

Input: AnnData (.h5ad) or CSV count matrices with cell-type annotations and donor IDs. Processing: filter cells by QC metrics, sum raw counts per cell type per donor, normalize by library size (CPM), log2(x+1) transform. Output: one expression matrix per cell type formatted identically to a GEO series matrix.

---

## 18. Critical Review and Honest Assessment

An independent Claude Code agent conducted a detailed critical review after completing both replication and novel disease validation. Nine concerns were identified. None are dealbreakers; several have been addressed through documentation updates.

The full critical review with author responses is included in Appendix C.

The most important criticism — "it finds what you already suspect" — is partially correct. The seeded mode is a refinement and replication engine. However, the criticism undersells the tool's ability to discover novel candidates through co-clustering. FAM107A and WDR49 in IPF were not in the seed list; they emerged from clustering with established IPF genes.

The recommended framing is: lead with "validation and filtering" and let discovery be the pleasant surprise, not the promise.

---

## 19. Known Limitations and Honest Boundaries

**Bulk tissue resolution.** Without cell-type deconvolution or snRNA-seq pseudo-bulking, the engine cannot distinguish true regulatory changes from cell-type composition artifacts.

**Welch's t-test approximation.** Not optimal for count-based RNA-seq data. Mitigated by cross-dataset replication. Perfectly appropriate for microarray data.

**Microarray optimization.** RNA-seq datasets require manual reconstruction. Most new transcriptomic data is RNA-seq or single-cell. The tool is optimized for retrospective mining of public microarray data.

**Memory scaling.** The consensus matrix scales as O(n²). On 8 GB hardware, blind runs are limited to approximately 9,000–10,000 study genes.

**Single-developer risk (bus factor of 1).** The codebase has 300+ passing tests, has been independently replicated multiple times, and cold-start reproduction is confirmed. JOSS submission is planned for September 2026.

**AGPL-3.0 license.** Derivative works deployed as services must be open-sourced. Labs with industry collaborations should consult their technology transfer office.

**No causal inference.** The engine identifies correlated expression patterns, not causal mechanisms.

**No clinical context integration.** No survival data, drug response, or mutation co-occurrence integration. Future work.

---

## 20. Technical Specifications and Infrastructure

| Specification | Value |
|--------------|-------|
| Language | Python 3.11+ |
| License | AGPL-3.0 (OSI-approved) |
| Tests | 300+ (pytest) |
| Core dependencies | numpy, pandas, scipy, scikit-learn, PyYAML, matplotlib |
| Clustering dependencies | umap-learn, hdbscan (optional install) |
| UI dependencies | FastAPI, uvicorn, jinja2 (optional install) |
| snRNA-seq dependencies | scanpy (optional) |
| Minimum hardware | Raspberry Pi 5, 8 GB RAM |
| Recommended | 16 GB+ RAM for genome-wide blind runs |
| CLI entry point | riker run config.yaml |
| Config format | YAML |
| Output format | JSON (summary, QC) + CSV (gene lists, meta-analysis) |
| Repository | github.com/RaySigmon/Riker_Engine |
| Version | v0.3.2 |

---

## 21. Development Methodology: AI-Assisted Engineering

The Riker Engine was developed using a dual-agent architecture:

**Kai (Claude Code CLI):** The primary autonomous builder, running on the Raspberry Pi 5. Kai maintains project-level memory via CLAUDE.md files and executes all coding, testing, and deployment tasks.

**Claude.ai:** The independent auditor and technical advisor. Claude.ai's memories are entirely separate from Kai's project context. This separation is deliberately exploited to prevent confirmation bias — when validating results, the auditor has no knowledge of what the builder "expects" to find.

The methodology is the author's; the AI assisted with implementation. Key design decisions — full-seed FDR, consensus clustering, directional elimination, pre-specification, the six-phase progressive filtering architecture — are human contributions. The AI tools accelerated implementation of these designs and provided independent verification.

---

## 22. Licensing, Sustainability, and Open-Core Model

The engine is released under AGPL-3.0: the engine itself is and will remain open source. An open-core model is planned where a future web-based SaaS platform for non-CLI users may be offered as a commercial product.

The AGPL-3.0 license requires that derivative works deployed as services be open-sourced. For academic research use, this imposes no restrictions. For industry collaborations, institutional technology transfer offices should review the license terms.

---

## 23. Future Work and Roadmap

**Near-term (pre-JOSS submission, September 2026):**
- PI collaborator outreach for independent human validation
- Cancer validation as second paper (CRC data already available)
- Web UI improvements
- PyPI packaging for pip install

**Medium-term (post-publication):**
- RNA-seq first-class support (streamlined ingestion)
- Cell-type deconvolution integration
- Clinical context modules (survival, drug response)
- Additional rare disease validations through collaborator partnerships

**Long-term:**
- Spatial transcriptomics integration
- Multi-condition comparative analysis
- NSF SBIR Phase I application

---

## 24. Conclusion

The Riker Engine demonstrates that a single person, working from a Raspberry Pi in Northwest Arkansas, can build a tool that recovers established disease biology across eight unrelated conditions, survives independent replication with 99.47% gene overlap, and produces textbook-perfect results when run by agents who have never seen the tool before.

The engine is not the most sophisticated transcriptomic analysis tool ever built. It is the most accessible one that also enforces statistical rigor. It runs in under 5 minutes on consumer hardware, requires no institutional infrastructure, and produces results that are reproducible, transparent, and honestly reported.

The validation framework — progressive filtering, pre-specified replication, full-seed-set FDR, consensus clustering, directional elimination — is the primary contribution. The clustering algorithm is commodity software. The rigorous progressive filtering and honest reporting pipeline is what makes the results defensible.

Eight diseases. Three independent validators. Zero code modifications between diseases. Built on a Raspberry Pi named Ghost, using Claude Code CLI. Named after a boy who inspired a father to search across dimensions for answers.

---

## Appendix A: Validation Summary Table (8 Diseases)

| Disease | Tissue | Seeds | Datasets | Ph1 Yield | Core Genes | Ph5 Survival | Meta-Sig | Validated By |
|---------|--------|-------|----------|-----------|------------|-------------|----------|-------------|
| ASD | Brain | 1,267 | 7 | 11.1% | 35 | 57.1% | 9 | Author (77.1% blind) |
| T2D | Islets | 443 | 4 | 12.6% | 8 | 100% | 8 | Author (62.5% blind) |
| IBD | Mucosa | 762 | 6 | 53.5% | 304 | 99.3% | 296 | Author (97.7% blind) |
| AD | Brain | 801 | 5 | 54.7% | 394 | 86.3% | 312 | Author (98.2% blind) |
| Breast Ca. | Tumor | 653 | 5 | 34.8% | 152 | 100% | 121 | Author (99.3% blind) |
| IPF | Lung | 354 | 5 | 68.1% | 190 | 89.5% | 157 | Author + Gemini + Claude |
| Psoriasis | Skin | 96 | 5 | 62.5% | 50 | 100% | 28 | Gemini (94.0% blind) |
| CRC | Colon | 515 | 6 | 64.3% | 263 | 92.8% | 218 | Claude (97.7% blind) |

## Appendix B: Cold-Start Replication Data

See Section 13 for full details.

## Appendix C: Independent Critical Review (9 Points) with Author Responses

**1. Seed dependency:** The seeded mode is a refinement engine by design. Novel candidates emerge through co-clustering (e.g., FAM107A in IPF). Framing: "validation and filtering tool" with discovery as a bonus.

**2. Microarray-centric:** Valid. RNA-seq support exists but requires more setup. The microarray-first approach has merit for cross-study meta-analysis due to greater standardization.

**3. No cell-type resolution:** True. Same limitation as WGCNA and all bulk co-expression tools. The snRNA-seq module provides a path forward.

**4. Welch's t-test for RNA-seq:** Conservative rather than liberal. Cross-dataset replication mitigates weakness. Appropriate for microarray.

**5. Single developer, AGPL:** Bus factor of 1 is real. JOSS publication addresses peer review gap. AGPL implications documented for labs with industry collaborations.

**6. Phase 2 pathway data absent:** All validations ran without pathway data. Documented as optional enhancement. Impact is an untested variable — future work.

**7. Dependency fragility:** numpy<2 pin resolves the hdbscan ABI crash. Added to pyproject.toml.

**8. No causal inference:** Correct and not a flaw. No transcriptomic tool does causal inference. Documentation updated to be explicit.

**9. No clinical context:** Survival data, drug response, mutation co-occurrence are roadmap items. Not current gaps.

## Appendix D: Technical Specifications

See Section 20.

## Appendix E: References

Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate. JRSS:B, 57(1), 289–300.

Campello RJGB, Moulavi D, Sander J (2013). Density-based clustering based on hierarchical density estimates. PAKDD 2013.

DerSimonian R, Laird N (1986). Meta-analysis in clinical trials. Controlled Clinical Trials, 7(3), 177–188.

Fabregat A et al. (2018). The Reactome pathway knowledgebase. Nucleic Acids Research, 46(D1), D649–D655.

Kanehisa M, Goto S (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Research, 28(1), 27–30.

Langfelder P, Horvath S (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9, 559.

Liberzon A et al. (2015). The MSigDB hallmark gene set collection. Cell Systems, 1(6), 417–425.

Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550.

McInnes L, Healy J, Melville J (2018). UMAP: Uniform Manifold Approximation and Projection. arXiv:1802.03426.

Ritchie ME et al. (2015). limma powers differential expression analyses. Nucleic Acids Research, 43(7), e47.

Welch BL (1947). The generalization of Student's problem. Biometrika, 34(1–2), 28–35.

---

*Riker Engine v0.3.2 | github.com/RaySigmon/Riker_Engine | AGPL-3.0*

*Built on a Raspberry Pi 5 named Ghost, using Claude Code CLI.*

*Named after a boy who inspired a father to search across dimensions for answers.*
