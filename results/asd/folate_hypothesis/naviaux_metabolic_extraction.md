# Naviaux CDR/Suramin Metabolic Data Extraction

**Date:** 2026-04-17
**Purpose:** Structured extraction of metabolic pathway data from Naviaux publications for folate hypothesis cross-reference

---

## Source Papers

1. **Naviaux et al. (2017)** "Low-dose suramin in autism spectrum disorder: a small, phase I/II, randomized clinical trial." *Ann Clin Transl Neurol* 4(7):491-505. PMID: 28695149
2. **Naviaux (2019)** "Metabolic features and regulation of the healing cycle—A new model for chronic disease pathogenesis and treatment." *Mitochondrion* 46:278-297. PMID: 30099222
3. **Naviaux (2014)** "Metabolic features of the cell danger response." *Mitochondrion* 16:7-17. PMID: 23981537
4. **Naviaux et al. (2014b)** "Reversal of autism-like behaviors and metabolism in adult mice with single-dose antipurinergic therapy." *Transl Psychiatry* 4:e400. PMID: 24937094
5. **Naviaux et al. (2021)** "Metabolic and behavioral features of acute hyperpurinergia and the MIA mouse model of ASD." *PLOS ONE* 16(3):e0248771. PMID: 33735311
6. **Lingampelly et al. (2024)** "Metabolic network analysis of pre-ASD newborns and 5-year-old children with ASD." *Commun Biol* 7:1458. (Naviaux senior author)

---

## SECTION 1: Metabolic Pathways Disturbed in ASD and Corrected by Suramin

### From the 2017 SAT-1 Human Clinical Trial

**Study parameters:** 431 of 610 targeted metabolites measurable in plasma; 63 biochemical pathways analyzed; 48 metabolites significant at 6 weeks post-suramin.

#### Table: Top 22 Pathways Changed by Suramin at 6 Weeks (Table 3 in paper)

| Rank | Pathway | Measured | Hits | Fold Enrichment | Direction |
|------|---------|----------|------|-----------------|-----------|
| 1 | Purine metabolism | 26 | 5 | 1.7 | Corrected |
| 2 | SAM/SAH/methionine/cysteine/glutathione | 15 | 5 | 3.0 | Corrected |
| 3 | Microbiome metabolism | 18 | 4 | 2.0 | Corrected |
| 4 | Branch-chain amino acid metabolism | 12 | 4 | 3.0 | Corrected |
| 5 | Bile acid metabolism | 6 | 3 | 4.5 | Corrected |
| 6 | Fatty acid oxidation/synthesis | -- | -- | -- | Decreased |
| 7 | Amino acid metabolism (alanine) | -- | -- | -- | Corrected |
| 8 | Krebs cycle | -- | -- | 2.0 | Corrected |
| 9 | Pyrimidine metabolism | -- | -- | 2.0 | Corrected |
| 10 | Sphingomyelin metabolism | -- | -- | 0.5 | Decreased |
| 11 | 1-carbon/folate/formate/glycine/serine | -- | -- | 3.6 | Corrected |
| 12 | GABA/glutamate/arginine/ornithine/proline | -- | -- | 3.0 | Corrected |
| 13 | Tyrosine/phenylalanine metabolism | -- | -- | 6.0 | Corrected |
| 14 | Cholesterol/cortisol/nongonadal steroid | -- | -- | 1.1 | Corrected |
| 15 | Gamma-glutamyl and other dipeptides | -- | -- | 4.5 | Corrected |
| 16 | Histidine/histamine/carnosine metabolism | -- | -- | 2.2 | Corrected |
| 17 | Nitric oxide/superoxide/peroxide metabolism | -- | -- | 4.5 | Corrected |
| 18 | Tryptophan/kynurenine/serotonin/melatonin | -- | -- | 1.5 | Corrected |
| 19 | Glycolysis/gluconeogenesis | -- | -- | 1.3 | Corrected |
| 20 | Vitamin C (ascorbate) metabolism | -- | -- | 4.5 | Corrected |
| 21 | Amino-sugar/hexose metabolism | -- | -- | 1.8 | Corrected |
| 22 | Phospholipid metabolism | -- | -- | 0.1 | Decreased |

**Key finding:** 75% (21/28) of pathways altered in ASD children were also changed by suramin in MIA and Fragile X mouse models.

#### Named Metabolites from SAT-1 (Human, 2017)

| Metabolite | Domain | Direction in ASD | Suramin Corrected? |
|------------|--------|-----------------|-------------------|
| AICAR | Purine | Decreased | Yes (increased) |
| 1-Methyladenine (1-MA) | Purine | Decreased | Yes (increased) |
| cAMP | Purine | Elevated | Yes (decreased) |
| dGDP | Purine | Elevated | Yes (decreased) |

Note: The paper states Figure 4 shows the "rank order of the top 35 of 48 significant metabolites" but these are only available as a figure image, not enumerated in text. Supplementary Tables S2-S4 contain the full individual metabolite data.

---

### From the 2014 Mouse MIA Model (Naviaux et al. 2014b)

**Study parameters:** 478 metabolites from 44 pathways measured. 48 discriminant metabolites. Single suramin dose corrected 17 of 18 disturbed pathways.

#### Pathway Impact Breakdown (Mouse MIA Model)

| Rank | Pathway | % of Impact | Metabolites | Corrected by Suramin |
|------|---------|-------------|-------------|---------------------|
| 1 | Purine metabolism | 23% | 11 (9 increased, ATP & allantoin decreased) | 9/11 (82%) |
| 2 | Microbiome metabolism | 15.1% | 6 | 6/6 (100%) |
| 3 | Phospholipid metabolism | 8.4% | 4 | 4/4 (100%) |
| 4 | Bile salt metabolism | 7.9% | 3 | **0/3 (NOT corrected)** |
| 5 | Sphingolipid metabolism | 7.1% | 4 | 4/4 (100%) |
| 6 | Cholesterol/steroid metabolism | 7.0% | 4 | 4/4 (100%) |
| 7 | Glycolysis/gluconeogenesis | 5.4% | 3 | 3/3 (100%) |
| 8 | Tryptophan metabolism | 3.5% | -- | Yes |

Additional pathways corrected: pyrimidine, Krebs cycle, oxalate/glyoxylate, vitamin B3 (NAD), vitamin B2, GABA/glutamate/arginine, amino-sugar/galactose, SAM/SAH/methionine/cysteine/glutathione, thyroxine, biopterin/neopterin.

**Only 1 pathway NOT corrected:** Bile salt metabolism.

---

## SECTION 2: Comprehensive Metabolite Directory by Domain

### Domain 1: PURINES

| Metabolite | ASD Status | Direction | Suramin Corrected | CDR Phase | Source |
|------------|-----------|-----------|-------------------|-----------|--------|
| ATP (extracellular) | Disturbed | Elevated (pericellular) | Target of suramin | CDR1 trigger | 2014, 2019 |
| ATP (plasma) | Disturbed | Decreased | Yes | -- | 2014b mouse |
| AICAR | Disturbed | Decreased in ASD | Yes (increased) | -- | 2017 human |
| 1-Methyladenine | Disturbed | Decreased | Yes (increased) | -- | 2017 human |
| 1-Methyladenosine | Disturbed | Elevated at 4h post-ATP | -- | CDR1 | 2021 mouse |
| cAMP | Disturbed | Elevated | Yes (decreased) | -- | 2017 human |
| dGDP | Disturbed | Elevated | Yes (decreased) | -- | 2017 human |
| Inosine | Disturbed | Elevated in MIA | Over-corrected | -- | 2014b mouse |
| Deoxyinosine | Disturbed | Elevated in MIA | Over-corrected | -- | 2014b mouse |
| Allantoin | Disturbed | Decreased in MIA | Yes | -- | 2014b mouse |
| Xanthosine | Disturbed | Elevated at 30min | -- | CDR1 acute | 2021 mouse |
| Hypoxanthine | Disturbed | Elevated at 30min | -- | CDR1 acute | 2021 mouse |
| Xanthine | Disturbed | Elevated at 30min | -- | CDR1 acute | 2021 mouse |
| Uric acid | Disturbed | Elevated at 30min | -- | CDR1 acute | 2021 mouse |
| 7-Methylguanine | Disturbed | Elevated at 30min | -- | CDR1 acute | 2021 mouse |
| 7-Methylguanosine | Disturbed | Elevated in ASD | -- | -- | 2024 network |
| Pseudouridine | Disturbed | Elevated at 4h | -- | CDR1 | 2021 mouse |
| GMP | Disturbed | Network dysregulated | -- | -- | 2024 network |
| AMP | Disturbed | Network dysregulated | -- | -- | 2024 network |
| ADP | Disturbed | Network dysregulated | -- | -- | 2024 network |
| IMP | Disturbed | Network dysregulated | -- | -- | 2024 network |
| XMP | Disturbed | Network dysregulated | -- | -- | 2024 network |
| dGMP | Disturbed | Decreased in newborns | -- | -- | 2024 network |
| Guanine | Disturbed | Decreased in newborns | -- | -- | 2024 network |

### Domain 2: PYRIMIDINES

| Metabolite | ASD Status | Direction | Suramin Corrected | Source |
|------------|-----------|-----------|-------------------|--------|
| Orotic acid | Disturbed | Elevated at 30min | -- | 2021 mouse |
| Beta-alanine | Disturbed (pyrimidine product) | Elevated at 30min | -- | 2021 mouse |
| Pyrimidine pathway (aggregate) | Disturbed | -- | Yes (2.0x enrichment) | 2017 human |

### Domain 3: SPHINGOLIPIDS

| Metabolite | ASD Status | Direction | Suramin Corrected | Source |
|------------|-----------|-----------|-------------------|--------|
| Ceramides (15 species) | Disturbed | Decreased at 30min post-ATP | -- | 2021 mouse |
| Ceramides (multiple) | Disturbed | **Elevated** in ASD children | -- | 2024 network |
| Sphingomyelins | Disturbed | Decreased at 30min | -- | 2021 mouse |
| Sphingomyelins (parent) | Disturbed | Decreased in ASD children | -- | 2024 network |
| 2-Hydroxy sphingomyelins | Disturbed | Decreased | -- | 2024 network |
| SM(d18:1/26:0 OH) | Disturbed | Decreased | -- | 2024 network |
| SM(d18:1/18:2 OH) | Disturbed | Decreased | -- | 2024 network |
| Sphingosine-1-phosphate (S1P) | Disturbed | **Elevated** in ASD | -- | 2024 network |
| Monohexosylceramides (MHCs) | Disturbed | Elevated | -- | 2024 network |
| Trihexosylceramides (THCs) | Disturbed | Elevated | -- | 2024 network |
| Sphingolipid pathway (aggregate) | Disturbed | 7.1% of impact | Yes (4/4 in mouse) | 2014b mouse |
| Sphingomyelin pathway (aggregate) | Disturbed | -- | Pathway changed (0.5x) | 2017 human |

**CDR role:** Sphingolipid and cholesterol-enriched membrane lipid rafts act as rheostats for tuning cellular sensitivity to purinergic signaling (Naviaux 2019, 2023).

### Domain 4: PHOSPHOLIPIDS

| Metabolite | ASD Status | Direction | Suramin Corrected | Source |
|------------|-----------|-----------|-------------------|--------|
| Phosphatidylethanolamine (PE) species | Disturbed | Decreased at 30min; Elevated in 5yr ASD | -- | 2021, 2024 |
| Phosphatidylinositol (PI) species | Disturbed | Decreased at 30min; Elevated in 5yr ASD | -- | 2021, 2024 |
| Phosphatidylcholine (PC) species | Disturbed | Decreased in ASD | -- | 2024 network |
| Phosphatidylserine (PS) species | Disturbed | Increased at 30min; Decreased in ASD | -- | 2021, 2024 |
| PS(18:0/20:4) | Disturbed | Decreased | -- | 2024 network |
| Choline | Disturbed | Elevated at 30min | -- | 2021 mouse |
| Phosphorylcholine | Disturbed | Elevated at 30min; Elevated in ASD | -- | 2021, 2024 |
| Ethanolamine | Disturbed | Elevated at 30min | -- | 2021 mouse |
| Myoinositol | Disturbed | Elevated at 30min | -- | 2021 mouse |
| Cardiolipins | Disturbed | Decreased in ASD | -- | 2024 network |
| Bis(monoacylglycero)phosphates (BMPs) | Disturbed | Decreased at 30min | -- | 2021 mouse |
| BMP(18:1/16:1) | Disturbed | New correlations with eicosanoids in ASD | -- | 2024 network |
| Phospholipid pathway (aggregate) | Disturbed | 8.4% of impact | Yes (4/4 in mouse) | 2014b mouse |

### Domain 5: GANGLIOSIDES

| Metabolite | ASD Status | Direction | Source |
|------------|-----------|-----------|--------|
| GM3(d18:1/18:1) | Disturbed | Network dysregulated | 2024 network |
| GM1 | Accumulates with suramin exposure in dorsal root ganglion cultures | -- | In vitro data |

**Note:** Gangliosides were not a primary focus of the suramin trial metabolomics but are mentioned in network analysis and CDR framework as glycosphingolipids.

### Domain 6: EICOSANOIDS / OXYLIPIDS

| Metabolite | ASD Status | Direction | Source |
|------------|-----------|-----------|--------|
| 5-HETE | Disturbed | Elevated at 30min post-ATP | 2021 mouse |
| 8-HETE | Disturbed | Network dysregulated in ASD | 2024 network |
| 9-HETE | Disturbed | Network dysregulated in ASD | 2024 network |
| 12-HETE | Disturbed | Network dysregulated in ASD | 2024 network |
| 15-HETE | Disturbed | Network dysregulated in ASD | 2024 network |
| 18-HETE | Disturbed | Network dysregulated in ASD | 2024 network |
| 13S-HODE | Disturbed | Elevated at 30min | 2021 mouse |
| Anandamide (endocannabinoid) | Disturbed | Elevated at 30min | 2021 mouse |
| Prostaglandins (aggregate) | Disturbed | Network dysregulated | 2024 network |
| Leukotrienes (aggregate) | Disturbed | Network dysregulated | 2024 network |

**Key finding from 2024 network study:** ASD eicosanoid network had 4x more positive and 3x more negative correlations vs. TD. New negative correlations found between eicosanoids and L-asparagine, pentose phosphate pathway, 1-methyladenosine, and CoQ10.

### Domain 7: MICROBIOME METABOLITES

| Metabolite | ASD Status | Direction | Suramin Corrected | Source |
|------------|-----------|-----------|-------------------|--------|
| O-acetylserine | Disturbed | Elevated at 30min | -- | 2021 mouse |
| Isopropylmalic acid | Disturbed | Elevated at 30min | -- | 2021 mouse |
| Trimethylamine-oxide (TMAO) | Disturbed | Elevated at 30min | -- | 2021 mouse |
| Imidazoleacetic acid | Disturbed | Elevated at 30min | -- | 2021 mouse |
| 4-Hydroxyphenyllactic acid | Disturbed | Elevated at 30min | -- | 2021 mouse |
| 4-Hydroxyphenylpyruvic acid | Disturbed | Elevated at 30min | -- | 2021 mouse |
| Indoxyl-3-sulfate | Disturbed | Elevated at 30min; Decreased in 5yr ASD | -- | 2021, 2024 |
| Butyrylcarnitine | Disturbed | Decreased at 30min | -- | 2021 mouse |
| Vitamin K2 (menaquinone) | Disturbed | Decreased at 30min | -- | 2021 mouse |
| Quinolinic acid | Disturbed | Decreased in MIA model | Yes | 2014b mouse |
| Microbiome pathway (aggregate) | Disturbed | 15.1% of impact; 2.0x enrichment | Yes (6/6 in mouse; corrected in human) | 2014b, 2017 |

### Domain 8: REDOX / ONE-CARBON / FOLATE

| Metabolite | ASD Status | Direction | Suramin Corrected | Source |
|------------|-----------|-----------|-------------------|--------|
| Oxidized glutathione (GSSG) | Disturbed | Elevated at 4h post-ATP | -- | 2021 mouse |
| GSH/GSSG ratio | Disturbed | Decreased in ASD | -- | 2024 network |
| Cystine (CysS-SCys) | Disturbed | Elevated at 4h | -- | 2021 mouse |
| 5-Methyltetrahydrofolic acid (5-mTHF) | Disturbed | Elevated at 4h; **Decreased in ASD children** | -- | 2021, 2024 |
| Serine | Disturbed | Decreased at 30min (Z = -5.0) | -- | 2021 mouse |
| Glycine | Disturbed | Decreased at 30min | -- | 2021 mouse |
| Trimethylglycine (betaine) | Disturbed | Decreased at 30min | -- | 2021 mouse |
| Dimethylglycine | Disturbed | Elevated at 4h | -- | 2021 mouse |
| S-Adenosylmethionine (SAM) | Disturbed | Decreased in ASD | -- | 2024 network |
| Dimethylarginine | Disturbed | Dysregulated (NO pathway) | -- | 2024 network |
| 3-Methylthiopropionate | Disturbed | Dysregulated (polyamine/methionine salvage) | -- | 2024 network |
| Carnosine | Disturbed | Decreased in ASD | -- | 2024 network |
| CoQ10 | Disturbed | Decreased in ASD | -- | 2024 network |
| SAM/SAH/Met/Cys/GSH pathway | Disturbed | 3.0x enrichment | Yes | 2017 human |
| 1-Carbon/folate/formate/Gly/Ser pathway | Disturbed | 3.6x enrichment | Yes | 2017 human |

---

## SECTION 3: Cell Danger Response Framework

### The 8 Functional Changes of the Acute CDR (Naviaux 2014)

1. **Metabolic shift** — from polymer to monomer synthesis (prevents pathogen hijacking of cellular resources)
2. **Membrane stiffening** — limits superinfection and pathogen egress
3. **Antimicrobial release** — antiviral and antimicrobial chemicals deployed into pericellular environment
4. **Autophagy upregulation** — increased mitochondrial fission and mitophagy for pathogen removal
5. **Epigenetic reprogramming** — DNA methylation and histone modifications alter gene expression
6. **Retroelement mobilization** — endogenous retroviruses and LINEs activated for genetic variation
7. **Danger signaling** — extracellular nucleotides, H2O2, eicosanoids, metabolites, and cytokines warn neighboring/distant cells
8. **Behavioral modification** — altered host behavior to prevent spread; sleep pattern changes to facilitate healing

### Metabolic Cascade Domains (Naviaux 2014)

The CDR produces changes in:
- Cellular electron flow
- Oxygen consumption
- Redox balance
- Membrane fluidity
- Lipid dynamics
- Bioenergetics
- Carbon and sulfur resource allocation
- Protein folding and aggregation
- Vitamin availability
- Metal homeostasis
- Indole metabolism
- Pterin metabolism
- 1-Carbon and polyamine metabolism
- Polymer formation

### CDR Phase Characteristics (Naviaux 2019)

#### CDR1: Inflammation / Damage Control
- **Mitochondria:** M1 (punctate, proinflammatory)
- **Energy metabolism:** Anaerobic glycolysis
- **eATP levels:** HIGH
- **ROS:** HIGH
- **Lactate:** HIGH
- **mtDNA copy number:** LOW
- **FAO:** For ROS/NLRP3
- **Gap junctions:** Lost
- **Key signaling:** NFkB, NOX2 (NADPH oxidase), P2X7-pannexin channels
- **Key genes/proteins:** NFkB, NOX2, CD33-related SIGLECs, P2X7

#### CDR2: Proliferation / Biomass Replacement
- **Mitochondria:** M0 (intermediate, multipotential)
- **Energy metabolism:** Aerobic glycolysis (Warburg)
- **eATP levels:** MODERATE
- **ROS:** Intermediate
- **Lactate:** Intermediate
- **mtDNA copy number:** Intermediate
- **PPP:** HIGH (NADPH for biosynthesis)
- **FAO:** Synthesis > oxidation
- **Key signaling:** HIF1alpha, stem cell recruitment
- **Key genes/proteins:** HIF1alpha, RECQL4 helicase (interacts with p53/POLG)
- **Pathology if blocked:** Replicative senescence; fibrosis; cancer

#### CDR3: Differentiation / Resolution
- **Mitochondria:** M2 (filamentous, anti-inflammatory)
- **Energy metabolism:** Oxidative phosphorylation (restored)
- **eATP levels:** LOW
- **ROS:** LOW
- **Lactate:** LOW
- **mtDNA copy number:** HIGH
- **FAO:** For oxphos
- **Gap junctions:** Re-established
- **NAD+/NADP+:** Restored
- **Key signaling:** TGFbeta, adaptive immune regulation, circadian reintegration
- **Key genes/proteins:** NRF2, mTOR
- **Pathology if blocked:** Autoimmune disease, ME/CFS, psychiatric conditions

### Table: Mitochondrial Phenotypes by CDR Phase

| Trait | M1 (CDR1) | M0 (CDR2) | M2 (CDR3) |
|-------|-----------|-----------|-----------|
| Energy metabolism | Glycolysis | Aerobic glycolysis | Oxidative phosphorylation |
| Morphology | Punctate | Intermediate | Filamentous |
| ROS production | High | Intermediate | Low |
| Lactate release | High | Intermediate | Low |
| PPP function | NOX support | Biosynthesis | Redox maintenance |
| FAO function | ROS/NLRP3 | Synthesis > oxidation | For oxphos |
| mtDNA copy number | Low | Intermediate | High |
| Inflammatory state | Proinflammatory | Neutral | Anti-inflammatory |

### Purinergic Receptor System (CDR Framework)

19 purinergic receptors cloned (suramin is a non-selective inhibitor of most):

**P2X receptors (7 subtypes, ionotropic):**
P2X1, P2X2, P2X3, P2X4, P2X5, P2X6, P2X7

**P2Y receptors (8 subtypes, metabotropic GPCRs):**
P2Y1, P2Y2, P2Y4, P2Y6, P2Y11, P2Y12, P2Y13, P2Y14

**P1/Adenosine receptors (4 subtypes):**
A1, A2A, A2B, A3

**Key CDR-specific receptor findings:**
- P2X7-pannexin channels: nucleotide release in CDR1
- P2Y2, P2X7: expression abnormalities corrected by suramin in MIA mouse
- A2A (A2AR): sleep induction via adenosine
- A1 (A1R): referenced in CDR3 resolution
- P2X6: age-related nuclear translocation affecting mRNA splicing
- CD39/CD73: convert eATP -> ADP -> AMP -> adenosine (extracellular degradation)

### Additional CDR Genes/Proteins

| Gene/Protein | CDR Role |
|-------------|----------|
| CD38 | Increases with age; synthesizes cADPR/NAADP from NAD+/NADP+ |
| PARP | Poly ADP ribose polymerase; NAD+ consumption |
| AMPK | Stimulates autophagy/mitophagy; activated by AICAR |
| DRP1 | Mitochondrial fission regulation |
| ERK1/2 | Signal transduction; corrected by suramin in MIA model |
| CAMKII | Signal transduction; corrected by suramin in MIA model |
| C1q (complement) | Corrected by suramin in Fragile X model |
| TDP43 | Corrected by suramin in Fragile X model |
| APP (amyloid precursor) | Corrected by suramin in Fragile X model |
| NLRP3 | Inflammasome; CDR1 activation |

---

## SECTION 4: The 14 Shared Pathways in Pre-ASD Newborns and 5yr ASD Children (2024)

From Lingampelly/Naviaux 2024 metabolic network analysis:

These 14 pathways account for 79-80% of metabolic impact:

1. Sphingolipids (~25% of total impact)
2. Phospholipids (20-26% of total impact)
3. Ceramide metabolism
4. Glycosphingolipid pathway
5. Bile acid metabolism
6. Purine metabolism
7. Eicosanoid pathway
8. Glycolysis
9. Fatty acid oxidation
10. One-carbon metabolism
11. Cholesterol metabolism
12. Neurotransmitter pathways
13. Carnitine/acyl-carnitine metabolism
14. (Lipid metabolism aggregate: 63-71% of total impact)

**Critical finding:** 80% of the metabolic impact was caused by decreased anti-inflammatory and antioxidant defenses, and increased physiologic stress molecules (lactate, glycerol, cholesterol, ceramides).

---

## SECTION 5: Summary Cross-Reference for Folate Hypothesis

### Metabolites Directly Relevant to Folate/One-Carbon Metabolism

| Metabolite | Status in ASD | CDR Connection |
|------------|--------------|----------------|
| 5-mTHF (active folate) | DECREASED in ASD children | One-carbon pathway; CDR depletes folate intermediates |
| SAM | DECREASED in ASD | Methyl donor; CDR diverts methionine |
| Serine | DECREASED acutely | Feeds into folate cycle; consumed in CDR |
| Glycine | DECREASED acutely | Feeds into folate cycle; consumed in CDR |
| Betaine (trimethylglycine) | DECREASED acutely | Alternative methyl donor |
| Dimethylglycine | ELEVATED at 4h | Betaine demethylation product |
| GSH/GSSG ratio | DECREASED in ASD | Glutathione redox; requires cysteine from transsulfuration |
| GSSG | ELEVATED | Oxidative stress marker |
| Cystine | ELEVATED | Oxidized cysteine dimer |
| CoQ10 | DECREASED in ASD | Mitochondrial electron transport |
| Carnosine | DECREASED in ASD | Antioxidant dipeptide |
| Homocysteine | Not directly measured in suramin trials | Implied disturbed via SAM/SAH pathway ranking |

### Key Pathway Convergence

The suramin trial shows the **SAM/SAH/methionine/cysteine/glutathione pathway** ranked #2 (3.0x enrichment) and **1-carbon/folate/formate/glycine/serine** ranked #11 (3.6x enrichment) among pathways corrected by suramin in ASD children. This directly connects the CDR/purinergic framework to folate metabolism disruption.

---

## SECTION 6: Vitamins and Cofactors

| Metabolite | Direction in CDR/ASD | Source |
|------------|---------------------|--------|
| Thiamine (B1) | Elevated acutely | 2021 mouse |
| Niacin (B3) | Elevated acutely; Decreased in ASD | 2021, 2024 |
| Pyridoxic acid (B6 metabolite) | Elevated acutely | 2021 mouse |
| FAD (B2/riboflavin) | Elevated at 4h; Decreased in ASD | 2021, 2024 |
| L-Carnitine | Elevated at 4h; Decreased in ASD | 2021, 2024 |
| Vitamin C (ascorbate) pathway | Disturbed | 2017 human (4.5x enrichment) |
| Vitamin D3 | Decreased at 30min | 2021 mouse |
| Vitamin K2 (menaquinone) | Decreased at 30min | 2021 mouse |

---

## Notes on Data Limitations

1. **Supplementary tables S2-S4** from the 2017 human trial contain the full individual metabolite data but are not freely accessible in web-crawlable form; they are embedded in the journal's supplementary file system.
2. The 2014 mouse study names pathways but does not enumerate all 48 discriminant metabolites individually in the main text.
3. Ganglioside data is sparse in the suramin metabolomics — GM3 appears only in the 2024 network study; GM1 only in an in vitro suramin exposure context.
4. Specific eicosanoid species (PGE2, PGD2, LTB4, LTC4 etc.) are not individually named — only HETE species and aggregate pathway references appear.
5. p-Cresol sulfate, indole-3-propionate, and hippurate are not individually named in these papers; the microbiome metabolites that are named are listed in Domain 7 above.
6. The CDR phase-to-specific-metabolite mapping at the individual metabolite level exists primarily in the 2019 and 2023 papers (behind paywalls at ScienceDirect); what is available from open-access sources and reviews is captured here.

---

## Sources

- [Naviaux et al. 2017 - PMC5497533](https://pmc.ncbi.nlm.nih.gov/articles/PMC5497533/)
- [Naviaux et al. 2014b Mouse Study - PMC4080315](https://pmc.ncbi.nlm.nih.gov/articles/PMC4080315/)
- [Naviaux et al. 2021 Hyperpurinergia - PLOS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0248771)
- [Lingampelly/Naviaux 2024 Network Analysis - PMC11549098](https://pmc.ncbi.nlm.nih.gov/articles/PMC11549098/)
- [Naviaux 2014 CDR - PubMed 23981537](https://pubmed.ncbi.nlm.nih.gov/23981537/)
- [Naviaux 2019 Healing Cycle - PubMed 30099222](https://pubmed.ncbi.nlm.nih.gov/30099222/)
- [Naviaux 2023 Salugenesis - PubMed 37120082](https://pubmed.ncbi.nlm.nih.gov/37120082/)
- [Incomplete Healing/Aging 2019 - PMC6627909](https://pmc.ncbi.nlm.nih.gov/articles/PMC6627909/)
