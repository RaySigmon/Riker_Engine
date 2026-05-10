# Experimental Proposal Outline — Naviaux Collaboration

**STATUS: DRAFT OUTLINE — DO NOT SEND**
**Last updated:** 2026-05-10

---

## 1. Executive Summary (1 page)

- Engine found 40/41-gene mito/energy cluster at iron-clad stability in ASD brain
- CDR-consistent mixed direction pattern (23 up / 18 down) — Naviaux confirmed
- Version update (v0.3.2 → v0.3.3.2) retained 98% of findings, improved CYC1
- Propose staged experimental validation starting with lowest-cost data sharing

## 2. Computational Findings Update (2-3 pages)

- v0.3.3.2 update: documented dilution fix, 40/41 cluster retention
- HSPD1 lost (single gene), CYC1 improved borderline → iron-clad
- All 8 FDR-significant genes retained as iron-clad
- Direction pattern unchanged (deterministic)
- 394 iron-clad (up from 376)
- Cross-disease validation: 7 additional diseases, regime model

## 3. Central Hypothesis (1 page)

**If the CDR model is correct, the 15 surviving mito-core genes should normalize in expression when the CDR is pharmacologically resolved (suramin treatment).**

Pre-specified directional predictions:
- Genes upregulated in ASD brain → should decrease toward control levels with suramin
- Genes downregulated in ASD brain → should increase toward control levels with suramin

This is the "mechanistic bridge" prediction stated in preprint v2, which Naviaux engaged with.

## 4. Staged Experimental Design (3-5 pages)

### Tier 1: Computational analysis of existing data (cost: zero)

**Ask:** Does the Naviaux lab have existing transcriptomic data (RNA-seq or microarray) from suramin-treated vs untreated mouse brain tissue (from the poly(I:C) MIA mouse model or similar)?

**If yes:** We run the Riker Engine on that data ourselves. Check whether:
- The 15 mito-core genes show direction-reversal with suramin treatment
- The 41-gene expanded cluster shows coordinated normalization
- The stability profiling produces iron-clad genes that overlap with the human ASD set

**Deliverable:** Computational analysis + manuscript draft within 4-6 weeks of data access.

**Publishable outcome regardless of result:**
- Positive (genes normalize): "Suramin resolves the transcriptomic CDR signature — mechanistic bridge between metabolomic and transcriptomic levels"
- Negative (genes don't change): "CDR metabolomic resolution does not propagate to the stable transcriptomic signature, suggesting post-transcriptional regulation or chronic epigenetic imprinting"
- Mixed (some normalize, some don't): "Partial transcriptomic normalization identifies suramin-responsive vs suramin-resistant CDR components"

### Tier 2: New RNA-seq on archived tissue (cost: ~$3-5k, 4-6 weeks)

**Only if Tier 1 data doesn't exist.**

**Ask:** RNA-seq on existing archived mouse brain tissue from suramin treatment experiments.

- Suramin-treated MIA mouse brain (n=5-8)
- Untreated MIA mouse brain (n=5-8)
- Control (non-MIA) mouse brain (n=5-8)
- Brain region: prefrontal cortex or whole cortex (to match human discovery cohorts)

**We provide:** Complete computational analysis pipeline, stability profiling, manuscript draft.

**Lab provides:** RNA-seq library prep + sequencing + basic QC.

### Tier 3: Targeted functional validation (contingent on Tier 1/2 results)

**Only if Tier 1 or 2 shows suramin-responsive genes.**

Options (ranked by feasibility):
1. qRT-PCR validation of top suramin-responsive genes in independent samples
2. Enzyme activity assays for OxPhos genes that normalize (NDUFAF5/Complex I, IDH3A/Krebs, HK2/glycolysis)
3. Seahorse XF bioenergetic profiling in patient-derived cells

Specific design depends on Tier 1/2 results.

## 5. Power Analysis (1-2 pages)

### For Tier 1 (if transcriptomic data exists):
- Standard RNA-seq differential expression analysis
- Primary test: are the 15 mito-core genes enriched among suramin-responsive genes? (hypergeometric/GSEA)
- Secondary test: directional concordance (binomial, p=0.03 for 5/5 match by chance)
- No new sample size needed — uses existing data

### For Tier 2 (new RNA-seq):
- n=5-8 per group gives ~80% power to detect log2FC ≥ 0.5 at p<0.05 (based on typical mouse brain RNA-seq variance)
- Three groups × 6 samples = 18 samples minimum
- Sequencing depth: 20-30M reads per sample (standard for differential expression)

## 6. Publication Strategy (1 page)

### Target journal: Molecular Autism (Naviaux's suggestion)

**Paper structure:** Combined computational + experimental
- Computational: Riker Engine identifies stable mito/OxPhos cluster in human ASD brain (updated from preprint v2 with v0.3.3.2 validation)
- Experimental: Suramin treatment normalizes/doesn't normalize the transcriptomic CDR signature in mouse model
- Together: first evidence that the computationally-identified transcriptomic CDR signature responds (or doesn't) to pharmacological CDR resolution

### Authorship framework
- Ray Sigmon: first author (computational pipeline, analysis, manuscript draft)
- Naviaux lab member(s): experimental contributor(s)
- Robert K. Naviaux: co-senior or co-corresponding author
- Exact arrangement to be discussed

### Timeline
- Tier 1 (if data exists): 2-3 months to manuscript
- Tier 2 (new RNA-seq): 4-6 months to manuscript
- Tier 3 (functional): 6-12 months, separate paper or revision

## 7. What We Offer (1 page)

- Complete, validated, open-source pipeline (AGPL-3.0)
- 8-disease validation including the ASD analysis he's already reviewed
- Regime model explaining when the engine discovers vs characterizes
- Full manuscript drafting and computational analysis
- Pre-specified hypothesis with falsifiable predictions

## 8. What We Ask (1 page)

- Tier 1: data sharing (existing suramin mouse transcriptomic data)
- Tier 2 (if needed): RNA-seq on ~18 archived tissue samples
- Scientific guidance on CDR-specific interpretation
- Co-authorship

## Appendices

- A: Full 41-gene cluster table with v0.3.3.2 stability scores
- B: Cross-version comparison (v0.3.2 → v0.3.3.2)
- C: Regime model summary
- D: Riker Engine repository link and reproduction instructions

---

## Research Questions Before Finalizing

- [x] Does Naviaux lab have published suramin mouse transcriptomic data?
  **PubMed search (2026-05-10): NO published transcriptomic data found.**
  Three suramin mouse papers (2013 poly(I:C), 2014 MIA adult, 2015 Fmr1 KO) — all metabolomics-based.
  No GEO-deposited expression datasets from suramin experiments.
  **Implication:** Tier 1 may not be possible with published data. Ask Naviaux if unpublished transcriptomic data exists, or if archived tissue is available for Tier 2 RNA-seq. The email should explicitly ask this question.
- [x] What mouse model does his lab use? **Both poly(I:C) MIA model AND Fragile X (Fmr1 KO).**
- [ ] What brain regions does his lab typically analyze?
- [ ] Is the suramin trial (NCT02508259) associated with any transcriptomic data?
- [ ] What's the typical turnaround for RNA-seq in his lab/core facility?
- [ ] Does UCSD have a genomics core that could process the RNA-seq?
