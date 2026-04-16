# Genome-Wide Stability Profiling Identifies a Mitochondrial Oxidative Phosphorylation Cluster in Autism Spectrum Disorder Brain Transcriptomics

*Multi-Run Consensus Analysis Across Multiple Public ASD Cohorts*

Ray Sigmon
Independent Computational Researcher, Northwest Arkansas, USA
alphalabsincorp@gmail.com

Alpha Research Labs | April 2026 | Preprint v2 (revised)
Repository: github.com/RaySigmon/Riker_Engine | License: AGPL-3.0

---

## Abstract

Autism spectrum disorder (ASD) is genetically heterogeneous, with over 1,200 candidate genes cataloged in the SFARI database. This heterogeneity has made it difficult to identify convergent transcriptomic signatures that replicate across independent cohorts. We applied the Riker Engine — a condition-agnostic, six-phase progressive filtering pipeline — to perform a hypothesis-free, genome-wide stability analysis using three ASD postmortem brain discovery cohorts, one brain replication cohort, and three blood replication cohorts (332 ASD cases, 332 controls, 664 total samples). The pipeline was executed 50 times with varied stochastic parameters (unique UMAP seeds per run) to distinguish biologically stable gene co-expression modules from clustering artifacts. From approximately 26,800 expressed genes, 376 achieved ≥90% cross-run stability (iron-clad), including 352 genes not in the SFARI database. Among the novel findings, a 26-gene mitochondrial and energy metabolism core cluster reached iron-clad stability, encoding nearly the complete oxidative phosphorylation system: ATP synthase subunits and assembly factors, Complex III and IV components, TCA cycle enzymes, coenzyme Q biosynthesis enzymes, and mitochondrial solute carriers. When expanded to include glycolysis, V-ATPase, and additional energy metabolism genes at iron-clad stability, the full energy metabolism cluster comprises 41 genes (23 upregulated, 18 downregulated), consistent with a mixed metabolic state across cell types and CDR stages. These findings are consistent with the Cell Danger Response (CDR) hypothesis, which predicts mitochondrial metabolic dysfunction as central to ASD pathophysiology. Additional iron-clad novel candidates include PVALB (parvalbumin interneurons), C1QB and C1QC (complement-mediated synaptic pruning), and GABRA1 (GABAergic inhibitory signaling) — each corresponding to active ASD research directions. All data, code, and stability reports are publicly available.

## 1. Introduction

The SFARI Gene database catalogs over 1,200 genes implicated in ASD through genetic studies. Yet most transcriptomic analyses examine these genes individually or within single datasets, producing gene lists with limited cross-cohort reproducibility. The fundamental challenge is heterogeneity: different patients carry different genetic variants, and bulk brain tissue averages across cell types, diluting cell-type-specific signals.

Network-based approaches such as WGCNA identify co-expressed modules but operate on single datasets and produce modules containing thousands of genes. No widely adopted tool performs cross-dataset consensus clustering with built-in replication testing and meta-analysis to yield a minimal, replicated gene set suitable for experimental follow-up.

We applied the Riker Engine, a condition-agnostic pipeline validated on eight diseases across six tissue types, to ASD brain transcriptomics. Rather than reporting a single pipeline run, we executed 50 independent runs with varied stochastic parameters across the entire expressed genome to quantify per-gene stability and separate robust biology from clustering artifacts.

### 1.1 Connection to the Cell Danger Response

The Cell Danger Response (CDR) hypothesis, proposed by Naviaux (2014), elaborated as a cell biology framework (Naviaux, 2020), and formalized in a three-hit metabolic signaling model (Naviaux, 2025), posits that ASD arises from persistent activation of an evolutionarily conserved mitochondrial stress response. When the CDR becomes chronic, mitochondria shift from normal energy production to a defensive metabolic state, disrupting oxidative phosphorylation and releasing extracellular ATP as a danger signal. Naviaux's metabolomic work identifies this dysfunction from the systemic side — abnormal metabolite levels in blood reflecting mitochondrial impairment. The present study approaches from the complementary direction: which specific mitochondrial genes are consistently dysregulated in ASD brain tissue across multiple postmortem cohorts?

## 2. Methods

### 2.1 The Riker Engine Pipeline

The Riker Engine (v0.3.2, AGPL-3.0) implements a six-phase progressive filtering pipeline. Phase 1 applies Welch's t-test (Welch, 1947; exact t-distribution) independently per gene per dataset, retaining genes significant (p < 0.05) in two or more datasets. Phase 2 constructs a feature matrix from expression statistics. Phase 3 performs consensus clustering via UMAP (McInnes et al., 2018; n_components=2, min_dist=0.0) and HDBSCAN (Campello et al., 2013) across 15 configurations (3 n_neighbors values × 5 random seeds). Phase 4 applies 10,000-permutation significance testing with Bonferroni correction, progressive sensitivity analysis using full-seed-set FDR (Benjamini & Hochberg, 1995; correction over the entire seed gene count, not the reduced study set), and leave-one-dataset-out stability testing. Phase 5 tests core genes in held-out replication datasets for directional concordance. Phase 6 computes inverse-variance weighted meta-analysis with REML random-effects estimation (DerSimonian & Laird, 1986).

### 2.2 ASD Datasets

Seven ASD datasets from the Gene Expression Omnibus were used: three brain cortex discovery datasets, one brain cortex replication dataset, and three blood replication datasets. The core gene identification (Phases 1–4) is driven entirely by the three brain discovery datasets. Blood datasets enter only at Phase 5 for cross-tissue replication testing. SFARI Gene database candidate genes (accessed March 2026, 1,267 genes) were used for blind recovery cross-referencing but were not used as pipeline inputs in blind mode.

**Table 4. Cohort demographics.**

| Dataset | Role | Tissue | ASD | Controls | Total | Ages | Reference |
|---------|------|--------|-----|----------|-------|------|-----------|
| GSE28521 | Discovery | Brain (frontal cortex, temporal cortex, cerebellum) | 39 samples (19 donors) | 40 samples (17 donors) | 79 | 2–56 yrs | Voineagu et al., Nature 2011 |
| GSE28475 | Discovery | Prefrontal cortex | 52 | 61 | 113 | 2–56 yrs (young 2–14, adult 15–56), all male | Chow et al., PLoS Genet 2012 |
| GSE64018 | Discovery | Temporal cortex (BA41/42) | 12 | 12 | 24 | ASD 15–60, controls 16–56, age/sex-matched | Irimia et al., Cell 2014 |
| GSE102741 | Replication | Prefrontal cortex (DLPFC) | 13 | 39 | 52 | ASD 4–67, controls 2–69, mean 22, matched ±6 yrs | Wright et al., Transl Psychiatry 2017 |
| GSE18123 | Replication | Blood (LCL) | 104 | 82 | 186 | 2–22 yrs, mean 8.1 | Kong et al., PLoS ONE 2012 |
| GSE26415 | Replication | Blood (leukocytes) | 21 | 42† | 63 | Age-matched | Kuwano et al., PLoS ONE 2011 |
| GSE42133 | Replication | Blood (leukocytes) | 91 | 56 | 147 | 1–4 yrs (toddlers) | Kong et al., PLoS ONE 2012 |
| **Total** | | | **332** | **332** | **664** | | |

Discovery (brain): 103 ASD samples, 113 control samples across 3 cohorts. GSE28521 includes multiple brain regions per donor (32 frontal cortex, 26 temporal cortex, 21 cerebellum); the ASD transcriptomic signal in the original study was primarily cortical (Voineagu et al., 2011). Replication (brain): 13 ASD, 39 controls. Replication (blood): 216 ASD, 180 controls across 3 cohorts. †GSE26415 controls include 21 ASD-matched subjects and 21 healthy mothers of ASD children (ctrlMO); both groups are labeled "nonaustistic control" in GEO metadata (Kuwano et al., 2011). Note: GSE28521 and GSE28475 share overlapping donors from the Autism Tissue Program and Harvard Brain Bank.

### 2.3 Blind Genome-Wide Mode

Rather than seeding with SFARI candidate genes, the pipeline was run in blind mode using all 19,296 protein-coding genes as seeds. This eliminates any selection bias from the seed list and allows the pipeline to discover gene modules from the full expressed transcriptome without prior hypothesis.

### 2.4 Stability Profiling Protocol

The full pipeline was executed 50 times. Each run received a unique base random seed derived deterministically from a master seed (42) plus the run number, and five unique UMAP seeds were generated per run. This produces genuinely different UMAP embeddings and HDBSCAN clustering outcomes per run while maintaining full reproducibility (any researcher using the same master seed will obtain identical results). After all 50 runs, each gene was scored by the fraction of runs in which it appeared in the Phase 4 core gene set. Genes were classified as iron-clad (≥90%), borderline (50–89%), or stochastic (<50%).

## 3. Results

### 3.1 Stability Overview

All 50 runs completed successfully. Core gene counts ranged from 389 to 410 (mean 401, ±2.6%), confirming genuine stochastic variance from varied UMAP initializations. A total of 438 unique genes appeared in at least one run.

**Table 1. Stability classification across 50 independent runs with varied stochastic parameters.**

| Classification | Count | SFARI | Novel | At 100% |
|----------------|-------|-------|-------|---------|
| Iron-clad (≥90%) | 376 | 24 | 352 | 307 (289 novel + 18 SFARI) |
| Borderline (50–89%) | 33 | 4 | 29 | — |
| Stochastic (<50%) | 29 | 2 | 27 | — |
| **Total** | **438** | **30** | **408** | |

### 3.2 SFARI Gene Recovery

Of 1,267 SFARI candidate genes, 30 were recovered in the blind genome-wide analysis. 18 appeared in all 50 runs (100% stability): ACTL6B, APBA2, ATP1A1, ATP2B2, DPP10, EFR3A, GABRG2, ICA1, KIAA0232, LHX2, LIN7B, NPAS3, NR3C2, PRPF19, SLC12A5, SLC45A1, SNAP25, and ZNF385B. Six more appeared in 46–49/50 runs (ZFHX3, DEPDC5, PPP3CA, SLC25A12, ITPR1, SLC25A27). ATP2B2, originally identified in the curated ASD analysis as part of a calcium/complement co-expression module, achieved perfect stability in the blind genome-wide context, confirming its robustness independent of seed gene composition.

SEZ6L2, the second headline gene from the original curated analysis, appeared in only 14/50 blind runs (28% stability). This indicates that its co-clustering behavior is dependent on the seed gene context: it reliably co-clusters with SFARI genes in the curated mode but does not consistently co-cluster in the genome-wide setting. This is an honest finding about the gene's transcriptomic neighborhood, not a failure of the method.

### 3.3 Mitochondrial and Energy Metabolism Cluster

The most striking finding is a mitochondrial and energy metabolism cluster at iron-clad stability. The original 26-gene core (Table 2a) encodes nearly the complete oxidative phosphorylation system and was discovered without any hypothesis about mitochondrial involvement. When expanded to include glycolysis enzymes, V-ATPase subunits, and additional energy metabolism genes at iron-clad stability, the full cluster comprises 41 genes (Table 2b). Directions were verified as consistent across all 50 runs with zero flips.

**Table 2a. Original 26-gene mitochondrial/OxPhos core (25 iron-clad + 1 borderline).**

| Category | Gene | Stability | Direction | log2FC | Function |
|----------|------|-----------|-----------|--------|----------|
| ATP synthase | ATP5F1A | 50/50 | up | +0.19 | F1 subunit alpha |
| | ATP5F1B | 50/50 | up | +0.03 | F1 subunit beta |
| | ATPAF1 | 50/50 | up | +0.24 | Assembly factor 1 |
| Complex III | UQCRC1 | 50/50 | down | -0.11 | Core protein 1 |
| | UQCRC2 | 50/50 | up | +0.01 | Core protein 2 |
| | CYC1 | 33/50 | up | +0.06 | Cytochrome c1 (borderline) |
| Complex IV | COX7A1 | 47/50 | down | -0.15 | Cytochrome c oxidase 7A1 |
| Complex II / TCA | SDHA | 50/50 | down | -0.002 | Succinate dehydrogenase A |
| TCA cycle | IDH3A | 50/50 | down | -0.006 | Isocitrate DH 3 alpha |
| | IDH3B | 50/50 | up | +0.17 | Isocitrate DH 3 beta |
| | IDH3G | 50/50 | up | +0.06 | Isocitrate DH 3 gamma |
| | OGDHL | 50/50 | down | -0.25 | 2-oxoglutarate DH-like |
| | ME3 | 50/50 | up | +0.08 | Malic enzyme 3 (mito) |
| CoQ biosynthesis | COQ3 | 50/50 | up | +0.05 | Ubiquinone methyltransferase |
| | COQ6 | 50/50 | up | +0.23 | Ubiquinone monooxygenase |
| Mito. transport | SLC25A3 | 50/50 | up | +0.14 | Phosphate carrier |
| | SLC25A12 | 48/50 | down | -0.13 | Aspartate/glutamate carrier (SFARI) |
| | SLC25A27 | 46/50 | down | -0.21 | Uncoupling protein 4 (SFARI) |
| | SLC25A37 | 50/50 | up | +0.45 | Mitoferrin-1 (iron import) |
| | VDAC3 | 47/50 | down | -0.02 | Voltage-dependent anion channel 3 |
| Assembly/chaperone | CHCHD4 | 50/50 | up | +0.15 | Mito. disulfide relay |
| | CHCHD7 | 50/50 | up | +0.004 | Coiled-coil-helix domain |
| | HSPD1 | 50/50 | up | +0.64 | Mitochondrial chaperonin (Hsp60) |
| | NDUFAF5 | 45/50 | down | -0.005 | Complex I assembly factor |
| Other energy | CYBA | 50/50 | up | +0.87 | Cytochrome b-245 alpha (ROS) |
| | NME4 | 50/50 | up | +0.37 | Nucleoside diphosphate kinase D |

Only SLC25A12 and SLC25A27 are in the SFARI database. The remaining 24 have not been individually linked to ASD in the existing literature.

**Table 2b. Additional energy metabolism genes at iron-clad stability (expanded cluster).**

| Category | Gene | Stability | Direction | log2FC | Function |
|----------|------|-----------|-----------|--------|----------|
| TCA cycle | IDH1 | 50/50 | up | +0.61 | Isocitrate DH 1 (cytoplasmic) |
| | FH | 50/50 | up | +0.13 | Fumarase |
| Amino acid metabolism | GOT1 | 50/50 | down | -0.29 | Aspartate aminotransferase (cyto) |
| | GOT2 | 50/50 | up | +0.03 | Aspartate aminotransferase (mito) |
| V-ATPase (lysosomal) | ATP6V1A | 50/50 | down | -0.09 | V-type proton ATPase subunit A |
| | ATP6V1C1 | 50/50 | down | -0.04 | V-type proton ATPase subunit C1 |
| | ATP6V1E1 | 50/50 | down | -0.13 | V-type proton ATPase subunit E1 |
| Iron-sulfur cluster | NFS1 | 50/50 | down | -0.10 | Cysteine desulfurase |
| Glycolysis | HK2 | 49/50 | up | +0.62 | Hexokinase 2 (mito-bound) |
| | ALDOC | 50/50 | up | +0.10 | Aldolase C (brain) |
| | ENO2 | 50/50 | down | -0.04 | Enolase 2 (neuronal) |
| | GPI | 50/50 | down | -0.07 | Glucose-6-phosphate isomerase |
| | PFKM | 50/50 | up | +0.11 | Phosphofructokinase (muscle) |
| | PFKP | 46/50 | down | -0.02 | Phosphofructokinase (platelet) |
| | PYGB | 50/50 | down | -0.06 | Glycogen phosphorylase (brain) |

**Direction summary (full 41-gene cluster): 23 upregulated, 18 downregulated.** All directions verified as consistent across all 50 runs (zero genes changed direction between runs).

This cluster spans the core of the mitochondrial energy production system: the TCA cycle (IDH3A/B/G, OGDHL, SDHA, ME3, FH, IDH1), the electron transport chain (Complex I assembly: NDUFAF5; Complex II: SDHA; Complex III: UQCRC1/2, CYC1; Complex IV: COX7A1), ATP synthesis (ATP5F1A/B, ATPAF1), coenzyme Q biosynthesis (COQ3, COQ6), mitochondrial protein import and assembly (CHCHD4, HSPD1, NDUFAF5), solute transport (SLC25A3/12/27/37, VDAC3), ROS production (CYBA), lysosomal acidification (ATP6V1A/C1/E1), and glycolysis (HK2, ALDOC, ENO2, GPI, PFKM, PFKP, PYGB). Of the original 26-gene core, only two (SLC25A12 and SLC25A27) are in the SFARI database.

### 3.4 Additional Biologically Notable Findings

**Table 3. Biologically notable genes beyond the mitochondrial cluster.**

| Gene | Stability | ASD Literature | Biological Relevance |
|------|-----------|----------------|---------------------|
| PVALB | 50/50 | Active research area, not in SFARI | Parvalbumin interneuron marker. PV+ interneuron dysfunction is a leading ASD cellular hypothesis. |
| C1QB | 50/50 | 1–2 papers | Complement C1q subunit B. Complement-mediated synaptic pruning implicated in ASD. |
| C1QC | 50/50 | 1 paper | Complement C1q subunit C. Co-clusters with C1QB. |
| GABRA1 | 45/50 | Minimal | GABA-A receptor alpha-1. Excitation/inhibition imbalance is a core ASD theory. |
| GABRG2 | 50/50 | In SFARI | GABA-A receptor gamma-2. Confirms GABAergic involvement. |
| SNAP25 | 50/50 | In SFARI | Synaptic vesicle docking. Validates known synaptic biology. |
| ATP2B2 | 50/50 | In SFARI | Calcium extrusion pump. Original Riker Engine ASD finding, confirmed at 100% blind stability. |

## 4. Discussion

### 4.1 The Mitochondrial Cluster and the CDR

The central finding — a mitochondrial and energy metabolism cluster at iron-clad stability in blind ASD brain transcriptomics — is consistent with the Cell Danger Response hypothesis from a fundamentally different methodological angle. Naviaux's metabolomic work identified the CDR from blood metabolite profiles, measuring the systemic consequences of mitochondrial dysfunction. The present study identifies the specific mitochondrial genes that are consistently dysregulated in ASD brain tissue. These are complementary observations: one captures what the mitochondria are producing (metabolites), the other captures what the mitochondria are expressing (transcripts).

The mixed directionality pattern (23 upregulated, 18 downregulated) is consistent with the heterogeneous cellular composition of bulk brain tissue. Different cell types — neurons, interneurons, astroglia, microglia, oligodendroglia — may occupy different stages of the CDR (CDR1, CDR2, CDR3), each with distinct metabolic profiles. The upregulation of ATP synthase components (ATP5F1A/B, ATPAF1) alongside downregulation of electron transport chain components (UQCRC1, COX7A1, NDUFAF5) and V-ATPase subunits (ATP6V1A/C1/E1) suggests a complex compensatory response rather than uniform suppression. CYBA (+0.87, strongest effect) and HSPD1 (+0.64) — encoding ROS-generating NADPH oxidase and the mitochondrial stress chaperone, respectively — point toward active danger signaling consistent with CDR1.

A critical question is whether these genes would normalize in expression when the CDR is pharmacologically resolved. If suramin treatment (Naviaux et al., 2017) in mouse models produces corresponding changes in these specific transcripts, that would establish a mechanistic bridge between the metabolomic and transcriptomic levels of the CDR.

### 4.2 Convergence with Other ASD Research Directions

PVALB at 50/50 stability is notable because parvalbumin interneuron dysfunction is one of the most actively investigated cellular mechanisms in ASD. Multiple mouse models show reduced PV+ interneuron density and altered inhibitory circuit function. The engine found PVALB independently from raw expression data without being told to look for interneuron markers. C1QB and C1QC (complement cascade, 50/50) align with research on excessive synaptic pruning in ASD. GABRA1 (GABAergic receptor, 45/50) supports the excitation/inhibition imbalance hypothesis. Each of these was discovered through the same blind, genome-wide, multi-run stability protocol, suggesting they reflect genuine convergent biology rather than hypothesis-driven bias.

### 4.3 SEZ6L2 Context Dependence

SEZ6L2 appeared in only 14/50 blind runs despite being a core finding in the curated SFARI-seeded analysis. This indicates seed-context dependence: SEZ6L2 reliably co-clusters with other SFARI genes when the seed space is restricted to 1,267 candidates, but its co-expression neighborhood shifts in the full genome-wide context. This does not invalidate SEZ6L2 as an ASD-relevant gene — its brain-specific differential expression is statistically robust — but its transcriptomic module membership is not stable across the genome-wide landscape. We report this transparently as a methodological observation.

### 4.4 Limitations

This analysis uses bulk brain tissue transcriptomics, which cannot resolve cell-type-specific signals. The mitochondrial genes may reflect changes in specific neuronal or glial populations that are diluted by whole-tissue averaging. Single-nucleus RNA-seq pseudo-bulking could reveal whether the mitochondrial cluster is specific to particular cell types. The discovery datasets are microarray-based postmortem samples; RNA-seq and in vivo data may yield additional candidates. The pipeline identifies consistently differentially expressed genes, not causal drivers. Experimental validation in cell culture or animal models is required to establish functional roles. GSE28521 and GSE28475 share overlapping donors from the same brain banks; while they use different platforms and brain regions, they are not fully independent cohorts.

## 5. Pipeline Validation

The Riker Engine has been validated on eight diseases across six tissue types with zero code modifications: ASD (35 core genes curated, 77.1% blind recovery), T2D (8 genes, IAPP as #1 blind signal), IBD (304 genes, 97.7% blind), Alzheimer's disease (394 genes, 98.2% blind), breast cancer (152 core genes, 139 survived replication, independent subtype separation), IPF (190 genes, 86.3% cold replication), psoriasis (50 genes, validated independently by a Gemini CLI agent), and colorectal cancer (264 genes, validated independently by a Claude Code agent). The pipeline has been independently replicated by four separate testers who were given only the GitHub URL. It has 300 unit tests with CI/CD. Full details are available in the Riker Engine preprint (Sigmon, 2026) and at github.com/RaySigmon/Riker_Engine.

## 6. Conclusion

Genome-wide stability profiling across 50 independent pipeline runs identifies 376 iron-clad genes in ASD brain transcriptomics, including a 26-gene mitochondrial core cluster (expanded to 41 with glycolysis and V-ATPase) showing a mixed directional pattern (23 up, 18 down) consistent with heterogeneous CDR staging across cell types. This finding is consistent with the Cell Danger Response hypothesis and identifies specific molecular candidates for experimental investigation. The stability profiling methodology itself — running a stochastic pipeline multiple times with varied parameters and reporting per-gene stability scores — offers a general framework for distinguishing robust co-expression patterns from clustering artifacts in any transcriptomic analysis.

## 7. Data and Code Availability

All source code, configurations, stability reports, seed gene lists, and validation results are publicly available at https://github.com/RaySigmon/Riker_Engine. Key files for reproducing this analysis:

- **Stability profiling script:** `scripts/stability_profiling.py`
- **ASD blind configuration:** `configs/examples/asd_blind.yaml`
- **Per-gene stability scores:** `stability_ASD_blind/stability_scores.csv`
- **Per-run core gene counts:** `stability_ASD_blind/run_summary.csv`
- **Independent validation reports:** `results/INDEPENDENT_VALIDATION.md`, `results/REPLICATION_LOG.md`
- **SFARI seed gene list:** `data/seeds/asd_sfari_genes.csv`

All GEO datasets used are publicly accessible through the NCBI Gene Expression Omnibus. Results are reproducible using master seed 42.

## 8. Acknowledgments

This work was conducted independently without institutional affiliation or external funding. The author acknowledges the Gene Expression Omnibus, SFARI Gene database, and the researchers who generated and deposited the postmortem brain expression datasets used in this analysis. AI-assisted coding tools (Anthropic Claude, Claude Code CLI) were used for pipeline development; their contributions are documented in the project repository. The methodology, experimental design, and all analytical decisions are the author's.

*Named after his son Riker, who inspired a father to search across dimensions for answers.*

## References

Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate. JRSS:B, 57(1), 289–300.

Campello RJGB, Moulavi D, Sander J (2013). Density-based clustering based on hierarchical density estimates. PAKDD 2013.

Chow ML, Pramparo T, Winn ME, et al. (2012). Age-dependent brain gene expression and copy number anomalies in autism suggest distinct pathological processes at young versus mature ages. PLoS Genet, 8(3), e1002592.

DerSimonian R, Laird N (1986). Meta-analysis in clinical trials. Controlled Clinical Trials, 7(3), 177–188.

Irimia M, Weatheritt RJ, Ellis JD, et al. (2014). A highly conserved program of neuronal microexons is misregulated in autistic brains. Cell, 159(7), 1511–1523.

Kong SW, Collins CD, Shimizu-Motohashi Y, et al. (2012). Characteristics and predictive value of blood transcriptome signature in males with autism spectrum disorders. PLoS ONE, 7(12), e49475.

Kuwano Y, Kamio Y, Kawai T, et al. (2011). Autism-associated gene expression in peripheral leucocytes commonly observed between subjects with autism and healthy women having autistic children. PLoS ONE, 6(9), e24723.

McInnes L, Healy J, Melville J (2018). UMAP: Uniform Manifold Approximation and Projection. arXiv:1802.03426.

Naviaux RK (2014). Metabolic features of the cell danger response. Mitochondrion, 16, 7–17.

Naviaux RK (2020). Perspective: Cell danger response Biology. Mitochondrion, 51, 40–45.

Naviaux RK (2025). A 3-hit metabolic signaling model for the core symptoms of autism spectrum disorder. Mitochondrion, 87, 102096. ePub 2025 Nov 14.

Naviaux RK et al. (2017). Low-dose suramin in autism spectrum disorder: a small, phase I/II, randomized clinical trial. Ann Clin Transl Neurol.

Sigmon R (2026). Riker Engine: A condition-agnostic transcriptomics pipeline for discovering replicated gene modules. GitHub: github.com/RaySigmon/Riker_Engine.

Voineagu I, Wang X, Johnston P, et al. (2011). Transcriptomic analysis of autistic brain reveals convergent molecular pathology. Nature, 474, 380–384.

Welch BL (1947). The generalization of Student's problem. Biometrika, 34(1–2), 28–35.

Wright C, Shin JH, Rajpurohit A, et al. (2017). Altered expression of histamine signaling genes in autism spectrum disorder. Transl Psychiatry, 7(5), e1126.

---

Riker Engine v0.3.2 | github.com/RaySigmon/Riker_Engine | AGPL-3.0 | Built on a Raspberry Pi 5 named Ghost
