# Genome-Wide Stability Profiling Identifies a Mitochondrial Oxidative Phosphorylation Cluster in Autism Spectrum Disorder Brain Transcriptomics

Ray Sigmon

Alpha Research Labs, Northwest Arkansas, USA

Correspondence: alphalabsincorp@gmail.com

**Keywords:** autism spectrum disorder, mitochondrial dysfunction, oxidative phosphorylation, transcriptomics, stability profiling, cell danger response, gene co-expression, purinergic signaling, postmortem brain, meta-analysis

---

## Abstract

**Background:** Autism spectrum disorder (ASD) is genetically heterogeneous, with over 1,200 candidate genes in the SFARI database, making it difficult to identify convergent transcriptomic signatures that replicate across cohorts. The Cell Danger Response (CDR) hypothesis proposes persistent mitochondrial metabolic dysfunction as central to ASD pathophysiology [2–4]. We applied hypothesis-free, genome-wide stability analysis to test whether a mitochondrial transcriptomic signal replicates across independent postmortem brain cohorts.

**Methods:** The Riker Engine, a condition-agnostic six-phase pipeline [5], was applied to three ASD postmortem brain discovery cohorts, one brain replication cohort, and three blood replication cohorts (332 ASD samples, 332 control samples, 664 total). The pipeline was executed 50 times with varied stochastic parameters to distinguish stable co-expression modules from clustering artifacts. Each gene was scored by the fraction of runs in which it appeared in the core set. AI-assisted coding tools were used during development; the author made all analytical decisions.

**Results:** Of approximately 26,800 expressed genes, 376 achieved ≥90% cross-run stability (iron-clad), including 352 not in SFARI. A 26-gene mitochondrial core cluster spanning the oxidative phosphorylation system reached iron-clad stability, expanding to 41 genes with glycolysis and V-ATPase components (23 upregulated, 18 downregulated). Inverse-variance weighted meta-analysis across brain discovery cohorts showed significant pooled effects for CYBA, GPI, IDH3A, OGDHL, and SLC25A27, with high between-study heterogeneity (I² > 80%) for others. A targeted query of 27 purinergic signaling genes found no receptor at any stability tier, though platform probe limitations for P2RY12 and P2RX4 constrained detection in two of three cohorts.

**Limitations:** Bulk postmortem brain tissue cannot resolve cell-type-specific signals. Discovery cohorts (ages 2–60) fall outside the 18–36 month neurodevelopmental window emphasized by CDR. Two discovery cohorts share overlapping donors. Effect sizes show between-study heterogeneity for some genes; unweighted mean log2FC values should be interpreted directionally rather than as precise magnitudes.

**Conclusions:** Genome-wide stability profiling identifies a convergent mitochondrial and energy metabolism signature in ASD brain transcriptomics, consistent with the CDR hypothesis. Co-expression stability is distinct from effect size homogeneity; both are reported transparently.

---

## Background

The SFARI Gene database catalogs over 1,200 genes implicated in ASD through genetic studies. Yet most transcriptomic analyses examine these genes individually or within single datasets, producing gene lists with limited cross-cohort reproducibility. The fundamental challenge is heterogeneity: different patients carry different genetic variants, and bulk brain tissue averages across cell types, diluting cell-type-specific signals.

Network-based approaches such as weighted gene co-expression network analysis (WGCNA) identify co-expressed modules but operate on single datasets and produce modules containing thousands of genes [6]. No widely adopted tool performs cross-dataset consensus clustering with built-in replication testing and meta-analysis to yield a minimal, replicated gene set suitable for experimental follow-up.

We applied the Riker Engine, a condition-agnostic pipeline validated on eight diseases across six tissue types [5], to ASD brain transcriptomics. Rather than reporting a single pipeline run, we executed 50 independent runs with varied stochastic parameters across the entire expressed genome to quantify per-gene stability and separate robust biology from clustering artifacts.

### Connection to the Cell Danger Response

The Cell Danger Response (CDR) hypothesis, proposed by Naviaux [2], elaborated as a cell biology framework [3], and formalized in a three-hit metabolic signaling model [4], posits that ASD arises from persistent activation of an evolutionarily conserved mitochondrial stress response. When the CDR becomes chronic, mitochondria shift from normal energy production to a defensive metabolic state, disrupting oxidative phosphorylation and releasing extracellular ATP as a danger signal. Naviaux's metabolomic work identifies this dysfunction from the systemic side — abnormal metabolite levels in blood reflecting mitochondrial impairment. Prior transcriptomic evidence is also consistent: Ginsberg et al. [1] reported down-regulation of mitochondrial oxidative phosphorylation genes in postmortem ASD brain using 9 idiopathic ASD and 9 control samples from cerebellar hemisphere and Brodmann area 19 cortex, with complexes I, III, and ATP synthase particularly affected. Their WGCNA also identified a cerebellar gene module enriched for purinergic signaling and immune response ontology terms that correlated with ADI-R social interaction impairment scores (r = −0.94, p = 5 × 10⁻⁴). The present study tests whether this kind of mitochondrial signal replicates and extends across larger postmortem cohorts using an independent methodology: which specific mitochondrial genes are consistently dysregulated in ASD brain tissue, and do those genes show reproducible co-expression stability across 50 independent pipeline runs?

## Methods

### The Riker Engine pipeline

The Riker Engine (v0.3.2, AGPL-3.0) implements a six-phase progressive filtering pipeline [5]. Phase 1 applies Welch's t-test [7] (exact t-distribution) independently per gene per dataset, retaining genes significant (p < 0.05) in two or more datasets. Phase 2 constructs a feature matrix from expression statistics. Phase 3 performs consensus clustering via UMAP [8] (n_components=2, min_dist=0.0) and HDBSCAN [9] across 15 configurations (3 n_neighbors values × 5 random seeds). Phase 4 applies 10,000-permutation significance testing with Bonferroni correction, progressive sensitivity analysis using full-seed-set FDR [10] (correction over the entire seed gene count, not the reduced study set), and leave-one-dataset-out stability testing. Phase 5 tests core genes in held-out replication datasets for directional concordance. Phase 6 computes inverse-variance weighted meta-analysis with DerSimonian-Laird random-effects estimation [11] (REML fallback when k ≥ 3) across discovery datasets only.

### ASD datasets

Seven ASD datasets from the Gene Expression Omnibus were used (Table 1): three brain discovery datasets, one brain cortex replication dataset, and three blood replication datasets. Core gene identification (Phases 1–4) is driven entirely by the three brain discovery datasets. Blood datasets enter only at Phase 5 for cross-tissue replication testing. SFARI Gene database candidate genes (accessed March 2026, 1,267 genes) were used for blind recovery cross-referencing but were not used as pipeline inputs.

**Table 1. Cohort demographics.**

| Dataset | Role | Tissue | ASD | Controls | Total | Ages | Reference |
|---------|------|--------|-----|----------|-------|------|-----------|
| GSE28521 | Discovery | Brain (frontal, temporal, cerebellum) | 39 samples (19 donors) | 40 samples (17 donors) | 79 | 2–56 yrs | Voineagu et al. [12] |
| GSE28475 | Discovery | Prefrontal cortex | 52 | 61 | 113 | 2–56 yrs, all male | Chow et al. [13] |
| GSE64018 | Discovery | Temporal cortex (BA41/42) | 12 | 12 | 24 | 15–60 yrs | Irimia et al. [14] |
| GSE102741 | Replication | Prefrontal cortex (DLPFC) | 13 | 39 | 52 | 4–67 yrs | Wright et al. [15] |
| GSE18123 | Replication | Blood (LCL) | 104 | 82 | 186 | 2–22 yrs | Kong et al. [16] |
| GSE26415 | Replication | Blood (leukocytes) | 21 | 42† | 63 | Age-matched | Kuwano et al. [17] |
| GSE42133 | Replication | Blood (leukocytes) | 91 | 56 | 147 | 1–4 yrs | Kong et al. [16] |
| **Total** | | | **332** | **332** | **664** | | |

Discovery (brain): 103 ASD samples from 83 donors, 113 control samples from 90 donors across 3 cohorts. GSE28521 includes multiple brain regions per donor; all other cohorts have 1:1 sample-to-donor ratios. Replication (brain): 13 ASD, 39 controls. Replication (blood): 216 ASD, 180 controls. †GSE26415 controls include 21 ASD-matched subjects and 21 healthy mothers of ASD children. Note: GSE28521 and GSE28475 share overlapping donors from the Autism Tissue Program and Harvard Brain Bank.

### Blind genome-wide mode

Rather than seeding with SFARI candidate genes, the pipeline was run in blind mode using all 19,296 protein-coding genes as seeds. The merged expression matrix across discovery datasets contained approximately 26,800 genes; the 19,296 protein-coding subset was used as the seed list.

### Stability profiling protocol

The full pipeline was executed 50 times. Each run received a unique base random seed derived deterministically from a master seed (42) plus the run number, and five unique UMAP seeds were generated per run. After all 50 runs, each gene was scored by the fraction of runs in which it appeared in the Phase 4 core gene set. Genes were classified as iron-clad (≥90%), borderline (50–89%), or stochastic (<50%).

### Use of large language models

AI-assisted coding tools (Anthropic Claude, Claude Code CLI) were used during pipeline development and manuscript preparation. Their contributions are documented in the project repository [5]. The author reviewed all code, made all analytical decisions, and takes full responsibility for the content of this publication.

## Results

### Stability overview

All 50 runs completed successfully. Core gene counts ranged from 389 to 410 (mean 401, ±2.6%), confirming genuine stochastic variance from varied UMAP initializations. A total of 438 unique genes appeared in at least one run.

**Table 2. Stability classification across 50 independent runs.**

| Classification | Count | SFARI | Novel | At 100% |
|----------------|-------|-------|-------|---------|
| Iron-clad (≥90%) | 376 | 24 | 352 | 307 (289 novel + 18 SFARI) |
| Borderline (50–89%) | 33 | 4 | 29 | — |
| Stochastic (<50%) | 29 | 2 | 27 | — |
| **Total** | **438** | **30** | **408** | |

### SFARI gene recovery

Of 1,267 SFARI candidate genes, 30 were recovered in the blind genome-wide analysis. Eighteen appeared in all 50 runs (100% stability): ACTL6B, APBA2, ATP1A1, ATP2B2, DPP10, EFR3A, GABRG2, ICA1, KIAA0232, LHX2, LIN7B, NPAS3, NR3C2, PRPF19, SLC12A5, SLC45A1, SNAP25, and ZNF385B. Six more appeared in 46–49/50 runs (ZFHX3, DEPDC5, PPP3CA, SLC25A12, ITPR1, SLC25A27).

### Mitochondrial and energy metabolism cluster

The most striking finding is a mitochondrial and energy metabolism cluster at iron-clad stability. The original 26-gene core (Table 3) spans the oxidative phosphorylation system and adjacent mitochondrial metabolism and was discovered without any hypothesis about mitochondrial involvement. When expanded to include glycolysis enzymes, V-ATPase subunits, and additional energy metabolism genes at iron-clad stability, the full cluster comprises 41 genes (Table 4). Directions were verified as consistent across all 50 runs with zero flips.

**Table 3. 26-gene mitochondrial/OxPhos core (25 iron-clad + 1 borderline).**

| Category | Gene | Stability | Direction | Mean log2FC | Function |
|----------|------|-----------|-----------|-------------|----------|
| ATP synthase | ATP5F1A | 50/50 | up | +0.19 | F1 subunit alpha |
| | ATP5F1B | 50/50 | up | +0.03 | F1 subunit beta |
| | ATPAF1 | 50/50 | up | +0.24 | Assembly factor 1 |
| Complex III | UQCRC1 | 50/50 | down | −0.11 | Core protein 1 |
| | UQCRC2 | 50/50 | up | +0.01 | Core protein 2 |
| | CYC1 | 33/50 | up | +0.06 | Cytochrome c1 (borderline) |
| Complex IV | COX7A1 | 47/50 | down | −0.15 | Cytochrome c oxidase 7A1 |
| Complex II/TCA | SDHA | 50/50 | down | −0.002 | Succinate dehydrogenase A |
| TCA cycle | IDH3A | 50/50 | down | −0.006 | Isocitrate DH 3 alpha |
| | IDH3B | 50/50 | up | +0.17 | Isocitrate DH 3 beta |
| | IDH3G | 50/50 | up | +0.06 | Isocitrate DH 3 gamma |
| | OGDHL | 50/50 | down | −0.25 | 2-oxoglutarate DH-like |
| | ME3 | 50/50 | up | +0.08 | Malic enzyme 3 (mito) |
| CoQ biosynthesis | COQ3 | 50/50 | up | +0.05 | Ubiquinone methyltransferase |
| | COQ6 | 50/50 | up | +0.23 | Ubiquinone monooxygenase |
| Mito. transport | SLC25A3 | 50/50 | up | +0.14 | Phosphate carrier |
| | SLC25A12 | 48/50 | down | −0.13 | Aspartate/glutamate carrier (SFARI) |
| | SLC25A27 | 46/50 | down | −0.21 | Uncoupling protein 4 (SFARI) |
| | SLC25A37 | 50/50 | up | +0.45 | Mitoferrin-1 (iron import) |
| | VDAC3 | 47/50 | down | −0.02 | Voltage-dependent anion channel 3 |
| Assembly/chaperone | CHCHD4 | 50/50 | up | +0.15 | Mito. disulfide relay |
| | CHCHD7 | 50/50 | up | +0.004 | Coiled-coil-helix domain |
| | HSPD1 | 50/50 | up | +0.64 | Mitochondrial chaperonin (Hsp60) |
| | NDUFAF5 | 45/50 | down | −0.005 | Complex I assembly factor |
| Other energy | CYBA | 50/50 | up | +0.87 | Cytochrome b-245 alpha (ROS) |
| | NME4 | 50/50 | up | +0.37 | Nucleoside diphosphate kinase D |

Mean log2FC = unweighted cross-dataset average from Phase 1. Only SLC25A12 and SLC25A27 are in SFARI. Inverse-variance weighted meta-analysis with 95% CIs is reported in Additional file 2.

**Table 4. Additional energy metabolism genes at iron-clad stability.**

| Category | Gene | Stability | Direction | Mean log2FC | Function |
|----------|------|-----------|-----------|-------------|----------|
| TCA cycle | IDH1 | 50/50 | up | +0.61 | Isocitrate DH 1 (cytoplasmic) |
| | FH | 50/50 | up | +0.13 | Fumarase |
| Amino acid metab. | GOT1 | 50/50 | down | −0.29 | Aspartate aminotransferase (cyto) |
| | GOT2 | 50/50 | up | +0.03 | Aspartate aminotransferase (mito) |
| V-ATPase | ATP6V1A | 50/50 | down | −0.09 | V-type proton ATPase subunit A |
| | ATP6V1C1 | 50/50 | down | −0.04 | V-type proton ATPase subunit C1 |
| | ATP6V1E1 | 50/50 | down | −0.13 | V-type proton ATPase subunit E1 |
| Iron-sulfur cluster | NFS1 | 50/50 | down | −0.10 | Cysteine desulfurase |
| Glycolysis | HK2 | 49/50 | up | +0.62 | Hexokinase 2 (mito-bound) |
| | ALDOC | 50/50 | up | +0.10 | Aldolase C (brain) |
| | ENO2 | 50/50 | down | −0.04 | Enolase 2 (neuronal) |
| | GPI | 50/50 | down | −0.07 | Glucose-6-phosphate isomerase |
| | PFKM | 50/50 | up | +0.11 | Phosphofructokinase (muscle) |
| | PFKP | 46/50 | down | −0.02 | Phosphofructokinase (platelet) |
| | PYGB | 50/50 | down | −0.06 | Glycogen phosphorylase (brain) |

Mean log2FC = unweighted cross-dataset average from Phase 1.

**Direction summary (full 41-gene cluster): 23 upregulated, 18 downregulated.** All directions verified as consistent across all 50 runs (zero genes changed direction between runs).

### Meta-analysis of effect sizes

Inverse-variance weighted random-effects meta-analysis was computed across the three brain discovery cohorts for all genes surviving Phase 5 replication (Additional file 2). Of the 41 mitochondrial cluster genes, 19 survived Phase 5 and entered Phase 6 meta-analysis; the remaining 22 were eliminated at Phase 5 due to significant opposite-direction expression in blood replication cohorts (GSE18123 or GSE26415), consistent with expected brain-blood tissue divergence for brain-specific mitochondrial genes.

Among the 19 genes with meta-analysis data, CYBA (RE = 0.78, 95% CI [0.30, 1.26], p = 0.001), GPI (RE = −0.29, 95% CI [−0.43, −0.15], p = 3.4 × 10⁻⁵), IDH3A (RE = −0.35, 95% CI [−0.58, −0.12], p = 0.003), OGDHL (RE = −0.46, 95% CI [−0.81, −0.10], p = 0.011), and SLC25A27 (RE = −0.42, 95% CI [−0.63, −0.21], p = 7.3 × 10⁻⁵) showed significant pooled effects. Other genes showed high between-study heterogeneity (I² > 80% for ATP5F1A, PFKM, NFS1, ALDOC, FH, PYGB), reflecting differences in microarray platform, brain region, and age distribution across the three discovery cohorts. The stability finding (consistent co-clustering across 50 runs) is distinct from effect size homogeneity; both are reported in Additional file 2.

### Additional notable findings

**Table 5. Biologically notable genes beyond the mitochondrial cluster.**

| Gene | Stability | ASD literature | Biological relevance |
|------|-----------|----------------|---------------------|
| PVALB | 50/50 | Active research area | Parvalbumin interneuron marker; PV+ interneuron dysfunction is a leading ASD cellular hypothesis [18,19] |
| C1QB | 50/50 | 1–2 papers | Complement C1q subunit B; complement-mediated synaptic pruning implicated in ASD [20] |
| C1QC | 50/50 | 1 paper | Complement C1q subunit C; co-clusters with C1QB |
| GABRA1 | 45/50 | Minimal | GABA-A receptor alpha-1; excitation/inhibition imbalance is a core ASD theory [21] |
| GABRG2 | 50/50 | In SFARI | GABA-A receptor gamma-2; confirms GABAergic involvement |
| SNAP25 | 50/50 | In SFARI | Synaptic vesicle docking |
| ATP2B2 | 50/50 | In SFARI | Calcium extrusion pump |

### Purinergic receptor query

We queried all 27 canonical purinergic signaling genes across the three discovery cohorts (7 P2X receptors, 10 P2Y receptors, 4 adenosine/P1 receptors, 3 pannexins, and 3 ecto-enzymes). No purinergic receptor reached core gene status in any of the 50 stability runs. Two key brain-expressed receptors — P2RY12 (the canonical microglial purinergic receptor) and P2RX4 (expressed on microglia, astroglia, and oligodendroglia) — lack probes on the GPL6883 microarray used by two of three discovery cohorts, creating a platform-level detection gap. ADA (adenosine deaminase) was differentially expressed in all three discovery cohorts with consistent upregulation (mean log2FC = 0.76) but did not achieve core stability, remaining at the lowest clustering level. Full results including per-gene evidentiary tier classifications and platform probe availability are reported in Additional file 1.

## Discussion

### The mitochondrial cluster and the CDR

The central finding — a mitochondrial and energy metabolism cluster at iron-clad stability in blind ASD brain transcriptomics — is consistent with the Cell Danger Response hypothesis from a fundamentally different methodological angle. Naviaux's metabolomic work identified the CDR from blood metabolite profiles, measuring the systemic consequences of mitochondrial dysfunction [2–4]. The present study identifies specific mitochondrial genes that are consistently dysregulated in ASD brain tissue. These are complementary observations: one captures what the mitochondria are producing (metabolites), the other captures what the mitochondria are expressing (transcripts).

The mixed directionality pattern (23 upregulated, 18 downregulated) is consistent with a summed cell-specific purinergic response in bulk brain tissue. Different brain cell types express distinct subsets of purinergic receptors: microglia express P2Y12, P2X4, P2X7, and adenosine A1/A2A/A3 receptors; astroglia express a broader set including P2Y1/2/4/6/12/13/14, P2X1/2/4/7, and A1/A2A/A2B/A3; and oligodendroglia express a partially overlapping repertoire [4, Supplementary Tables 1–2]. Extracellular ATP signaling through these different receptor subsets drives cell-specific metabolic responses — some cells shift toward glycolysis, others toward oxidative phosphorylation — producing the mixed directional pattern observed in bulk tissue [4, Supplementary Tables 1–2]. The upregulation of ATP synthase components alongside downregulation of electron transport chain components and V-ATPase subunits is consistent with active danger signaling and a mixed metabolic state across cell types, as predicted by the CDR framework [2–4]. CYBA (+0.87, strongest effect) and HSPD1 (+0.64) — encoding ROS-generating NADPH oxidase and the mitochondrial stress chaperone, respectively — are consistent with ongoing CDR activation.

The 41-gene cluster represents a potential transcriptomic readout for antipurinergic therapy (APT) studies. If suramin treatment [22] in mouse models produces corresponding changes in these specific transcripts, that would establish a mechanistic bridge between the metabolomic and transcriptomic levels of the CDR.

### Convergence with other ASD research directions

PVALB at 50/50 stability is notable because parvalbumin interneuron dysfunction is one of the most actively investigated cellular mechanisms in ASD, with multiple mouse models showing reduced PV+ interneuron density and altered inhibitory circuit function [18,19]. C1QB and C1QC (complement cascade, 50/50) align with research on complement-mediated synaptic pruning in ASD [20]. GABRA1 (GABAergic receptor, 45/50) supports the excitation/inhibition imbalance hypothesis [21]. Each was discovered through the same blind, genome-wide, multi-run stability protocol, suggesting genuine convergent biology rather than hypothesis-driven bias.

### SEZ6L2 context dependence

SEZ6L2 appeared in only 14/50 blind runs despite being a core finding in the curated SFARI-seeded analysis. This indicates seed-context dependence: SEZ6L2 reliably co-clusters with other SFARI genes when the seed space is restricted to 1,267 candidates, but its co-expression neighborhood shifts in the full genome-wide context. This does not invalidate SEZ6L2 as an ASD-relevant gene but indicates its module membership is context-dependent.

## Limitations

This analysis uses bulk brain tissue transcriptomics, which cannot resolve cell-type-specific signals. The mitochondrial genes may reflect changes in specific neuronal or glial populations that are diluted by whole-tissue averaging. Single-nucleus RNA-seq pseudo-bulking could reveal whether the mitochondrial cluster is specific to particular cell types.

The CDR hypothesis as formulated by Naviaux [4] emphasizes a critical neurodevelopmental window from late first trimester through approximately 18–36 months of age. Our discovery cohorts range in age from 2 to 60 years, with mean ages substantially outside this window. The mitochondrial signature we detect therefore represents a postmortem transcriptomic state in chronic CDR rather than the developmental window in which CDR stacking is hypothesized to drive ASD emergence. Prospective or earlier-age samples would be needed to test whether this signature is present during the critical window.

The discovery datasets are microarray-based postmortem samples; RNA-seq and in vivo data may yield additional candidates. GSE28521 and GSE28475 share overlapping donors from the same brain banks; while they use different platforms and brain regions, they are not fully independent cohorts.

Effect size estimates vary across cohorts for some genes. The unweighted mean log2FC values reported in the main tables should be interpreted as directional indicators rather than precise magnitudes. Inverse-variance weighted meta-analysis with 95% confidence intervals (Additional file 2) shows significant pooled effects for a subset of genes (e.g., CYBA, GPI, IDH3A, OGDHL, SLC25A27) but high between-study heterogeneity (I² > 80%) for others (e.g., ATP5F1A, PFKM, NFS1). This heterogeneity reflects differences in microarray platform, brain region, and age distribution across the three discovery cohorts. The stability finding — consistent co-clustering across 50 independent runs — is methodologically distinct from effect size homogeneity and should not be conflated with it.

Twenty-two of 41 mitochondrial cluster genes were eliminated at Phase 5 replication due to significant opposite-direction expression in blood replication cohorts, consistent with expected brain-blood tissue divergence for brain-specific transcripts. Phase 6 meta-analysis and confidence intervals are therefore available for only 19 of 41 genes. The pipeline identifies consistently differentially expressed genes, not causal drivers; experimental validation is required to establish functional roles.

## Conclusions

Genome-wide stability profiling across 50 independent pipeline runs identifies 376 iron-clad genes in ASD brain transcriptomics, including a 26-gene mitochondrial core cluster (expanded to 41 with glycolysis and V-ATPase) showing a mixed directional pattern (23 up, 18 down) consistent with summed cell-specific metabolic responses across brain cell types, as predicted by the CDR framework. No purinergic receptor reached stability-supported status, though platform limitations constrained detection of key brain-expressed receptors. The stability profiling methodology — running a stochastic pipeline multiple times with varied parameters and reporting per-gene stability scores — offers a general framework for distinguishing robust co-expression patterns from clustering artifacts in any transcriptomic analysis.

## Declarations

### Ethics approval and consent to participate

Not applicable. This study is a secondary analysis of publicly available, de-identified gene expression data from the NCBI Gene Expression Omnibus.

### Consent for publication

Not applicable.

### Availability of data and materials

All source code, configurations, stability reports, seed gene lists, and validation results are publicly available.

The dataset(s) supporting the conclusions of this article are available in the GitHub repository [5] and archived in Zenodo (https://doi.org/10.5281/zenodo.19623672).

- **Project name:** Riker Engine
- **Project home page:** https://github.com/RaySigmon/Riker_Engine [5]
- **Archived version:** https://doi.org/10.5281/zenodo.19623672
- **Operating system:** Platform independent
- **Programming language:** Python
- **License:** AGPL-3.0
- **Restrictions:** None

All GEO datasets used are publicly accessible through the NCBI Gene Expression Omnibus under accession numbers GSE28521, GSE28475, GSE64018, GSE102741, GSE18123, GSE26415, and GSE42133. Results are reproducible using master seed 42.

### Competing interests

The author declares no competing interests.

### Funding

This work was conducted independently without institutional affiliation or external funding.

### Authors' contributions

RS conceived and designed the study, developed the computational pipeline, performed all analyses, and wrote the manuscript.

### Acknowledgements

The author acknowledges the Gene Expression Omnibus, the SFARI Gene database, and the researchers who generated and deposited the postmortem brain and blood expression datasets used in this analysis. AI-assisted coding tools (Anthropic Claude, Claude Code CLI) were used during pipeline development and manuscript preparation; their contributions are documented in the project repository.

## References

1. Ginsberg MR, Rubin RA, Falcone T, Ting AH, Natowicz MR. Brain transcriptional and epigenetic associations with autism. PLoS ONE. 2012;7(9):e44736.
2. Naviaux RK. Metabolic features of the cell danger response. Mitochondrion. 2014;16:7–17.
3. Naviaux RK. Perspective: Cell danger response biology. Mitochondrion. 2020;51:40–45.
4. Naviaux RK. A 3-hit metabolic signaling model for the core symptoms of autism spectrum disorder. Mitochondrion. 2025;87:102096.
5. Sigmon R. Riker Engine: A condition-agnostic transcriptomics pipeline for discovering replicated gene modules. GitHub. 2026. https://github.com/RaySigmon/Riker_Engine. Accessed 16 Apr 2026.
6. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics. 2008;9:559.
7. Welch BL. The generalization of Student's problem. Biometrika. 1947;34(1–2):28–35.
8. McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for dimension reduction. arXiv. 2018;1802.03426.
9. Campello RJGB, Moulavi D, Sander J. Density-based clustering based on hierarchical density estimates. In: PAKDD 2013. Berlin: Springer; 2013.
10. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B. 1995;57(1):289–300.
11. DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7(3):177–188.
12. Voineagu I, Wang X, Johnston P, Lowe JK, Tian Y, Horvath S, et al. Transcriptomic analysis of autistic brain reveals convergent molecular pathology. Nature. 2011;474:380–384.
13. Chow ML, Pramparo T, Winn ME, Barnes CC, Li HR, Weiss L, et al. Age-dependent brain gene expression and copy number anomalies in autism suggest distinct pathological processes at young versus mature ages. PLoS Genet. 2012;8(3):e1002592.
14. Irimia M, Weatheritt RJ, Ellis JD, Parikshak NN, Gonatopoulos-Pournatzis T, Babber M, et al. A highly conserved program of neuronal microexons is misregulated in autistic brains. Cell. 2014;159(7):1511–1523.
15. Wright C, Shin JH, Rajpurohit A, Deep-Soboslay A, Collado-Torres L, Brandon NJ, et al. Altered expression of histamine signaling genes in autism spectrum disorder. Transl Psychiatry. 2017;7(5):e1126.
16. Kong SW, Collins CD, Shimizu-Motohashi Y, Koze IA, Barber S, Bhatt CK, et al. Characteristics and predictive value of blood transcriptome signature in males with autism spectrum disorders. PLoS ONE. 2012;7(12):e49475.
17. Kuwano Y, Kamio Y, Kawai T, Katsuura S, Inada N, Takaki A, et al. Autism-associated gene expression in peripheral leucocytes commonly observed between subjects with autism and healthy women having autistic children. PLoS ONE. 2011;6(9):e24723.
18. Filice F, Vörckel KJ, Sungur AÖ, Wöhr M, Bhatt D, Bhatt PS. Reduction in parvalbumin expression not loss of the parvalbumin-expressing GABA interneuron subpopulation in genetic parvalbumin and shank mouse models of autism. Mol Brain. 2016;9:10.
19. Hashemi E, Ariza J, Rogers H, Noctor SC, Martínez-Cerdeño V. The number of parvalbumin-expressing interneurons is decreased in the prefrontal cortex in autism. Cereb Cortex. 2017;27(3):1931–1943.
20. Magdalon J, Mansur F, Teles E Silva AL, de Goes VA, Reber JF, Bhatt PS. Complement system in brain architecture and neurodevelopmental disorders. Front Neurosci. 2020;14:23.
21. Rubenstein JL, Merzenich MM. Model of autism: increased ratio of excitation/inhibition in key neural systems. Genes Brain Behav. 2003;2(5):255–267.
22. Naviaux RK, Curtis B, Li K, Naviaux JC, Leber AT, Haas RH, et al. Low-dose suramin in autism spectrum disorder: a small, phase I/II, randomized clinical trial. Ann Clin Transl Neurol. 2017;4(7):491–505.

---

## Additional files

**Additional file 1** — Purinergic receptor query results across three ASD brain discovery cohorts (Supplementary Table S1 from preprint v3). Per-gene evidentiary tier classifications and platform probe availability.

**Additional file 2** — Inverse-variance weighted random-effects meta-analysis across three brain discovery cohorts for 41 mitochondrial and energy metabolism cluster genes. Includes random-effects pooled estimate, standard error, 95% confidence interval, p-value, and I² heterogeneity. For genes eliminated at Phase 5 replication (opposite direction in blood cohorts), Phase 5 verdict is reported in lieu of meta-analysis.
