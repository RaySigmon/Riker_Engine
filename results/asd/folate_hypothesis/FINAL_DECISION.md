# FINAL DECISION: Suramin × Riker Engine Cross-Reference

**Date:** 2026-04-17
**Test:** Does the Riker Engine ASD iron-clad gene set independently detect the same CDR metabolic pathways that suramin corrected?

---

## Result: DROP

The CDR-specific cross-reference does not meet the threshold for proceeding.

---

## The Numbers

| Metric | Value |
|--------|-------|
| Iron-clad genes tested | 376 |
| Suramin-response genes | 121 (across 8 CDR domains) |
| Overlap | 7 genes |
| Expected by chance | 2.3 |
| Fold enrichment | 3.08× |
| Overall p-value | 0.008 |
| Domains enriched (p<0.05, ≥2 hits) | **1 of 8** |
| Direction consistency | 3/3 = 100% |

### Domain-by-Domain

| Domain | Hits/Total | p-value | Enriched? |
|--------|-----------|---------|-----------|
| **Purines** | **4/31** | **0.003** | **YES** |
| Redox/one-carbon | 1/17 | 0.276 | no |
| Phospholipids | 1/18 | 0.290 | no |
| Microbiome-related | 1/10 | 0.173 | no |
| Sphingolipids | 0/17 | 1.000 | no |
| Pyrimidines | 0/8 | 1.000 | no |
| Gangliosides | 0/8 | 1.000 | no |
| Eicosanoids | 0/12 | 1.000 | no |

### The 7 Overlapping Genes

| Gene | CDR Domain | Direction | log2FC | Match? |
|------|-----------|-----------|--------|--------|
| PNP | Purines | up | **+0.935** | MATCH |
| NME4 | Purines | up | +0.371 | undetermined |
| GMPR2 | Purines | up | +0.121 | MATCH |
| ENTPD6 | Purines | up | +0.045 | undetermined |
| MAT2B | Redox/one-carbon | up | +0.130 | MATCH |
| KYAT3 | Microbiome | down | -0.137 | undetermined |
| LPCAT4 | Phospholipids | down | -0.087 | undetermined |

---

## Applying the Rubric

| Criterion | Required | Actual | Met? |
|-----------|----------|--------|------|
| ≥3 CDR domains enriched | 3+ | **1** | **NO** |
| >65% direction match | >65% | 100% (3/3) | YES |

**Result: 1 domain enriched. Does not meet ≥3 threshold. → DROP.**

---

## What This Means

### The purine arm IS there
The purines domain is significantly enriched (p=0.003, 6.9× fold). PNP (purine nucleoside phosphorylase, log2FC=+0.935) is the strongest metabolic gene in the overlap — a large, reproducible effect at a key purine degradation enzyme. GMPR2 (purine salvage) and ENTPD6 (extracellular nucleotide metabolism) add to the purine picture.

This is consistent with disrupted purine metabolism in ASD, which IS part of the CDR model. But purine disruption alone is not CDR — it could be downstream of many causes (neuroinflammation, mitochondrial dysfunction, oxidative stress).

### The CDR signature as a whole is NOT there
Zero sphingolipid genes. Zero ganglioside genes. Zero eicosanoid genes. Zero pyrimidine genes. Naviaux's CDR model predicts coordinated disruption across ALL of these domains. The Riker Engine only finds the purine piece. That's not independent computational validation of CDR — it's one overlapping metabolic domain that could be explained by simpler biology.

### The direction consistency is perfect but meaningless at n=3
3 out of 3 determined matches = 100%. But with only 3 data points, this is not statistically powered to mean anything. A coin flip would give 100% in 12.5% of trials.

---

## What Survives for the ASD Paper

These findings are still independently valuable even without the CDR framing:

1. **PNP upregulation (log2FC = +0.935, iron-clad)** — a strong, reproducible signal at a purine metabolism gene. PNP deficiency causes severe T-cell immunodeficiency. PNP upregulation in ASD brain tissue is consistent with altered purine catabolism. Worth reporting as a finding.

2. **The 59-gene immune signature** — this was never dependent on CDR and remains the strongest finding.

3. **MAT2B (ASD-specific, stable, direction-consistent)** — small effect but real.

4. **The disease specificity** — these signals are not in AD, BC, CRC, or psoriasis.

---

## Closing the CDR Investigation

The folate → CDR hypothesis investigation is complete. Here's what it produced:

### Worth keeping:
- 59 ASD-specific immune genes (complement, IFN, antigen presentation) → ASD disease paper
- MAT2B as ASD-specific one-carbon metabolism finding → ASD disease paper
- PNP as strong purine metabolism signal → ASD disease paper
- Disease specificity data → ASD disease paper
- The alternative explanation analysis (all 3 failed) → ASD disease paper

### Worth archiving:
- Receptor competition model → interesting but the folate hypothesis didn't pan out
- B12 peripheral sink model → novel but based on unreplicated HR 17.6
- UMFA → FRAA investigation → plausible but untested, becomes a discussion paragraph
- CDR biomarker cross-reference → purine-only, not CDR-broad
- All attack reports → documentation of rigorous process

### Not worth pursuing:
- CDR-specific framing for Naviaux communication
- Hypothesis paper in Medical Hypotheses as originally planned
- The full folate chain as a primary thesis

---

## What to Tell Naviaux

If continuing the conversation (which is worth doing based on the immune signature alone):

> "We ran a comprehensive cross-reference of our ASD iron-clad gene set against CDR metabolic pathways from your suramin work. The purine metabolism domain was significantly enriched (p=0.003), with PNP and ENTPD6 both iron-clad. However, the other CDR domains (sphingolipid, ganglioside, eicosanoid) showed no signal. We can't claim CDR-specific detection — the overlap is limited to purines.
>
> What we CAN report is an ASD-specific neuroinflammatory signature (59 iron-clad genes, 95% upregulated) with a complement/interferon/antigen-presentation pattern consistent with autoantibody-driven pathology. This connects to the FRα autoantibody literature and may be relevant to your work."

That's honest and still interesting. It just isn't the CDR validation we were hoping for.

---

## Redirect

Per the plan: close the folate/CDR investigation. Redirect energy to:
1. **ASD disease findings paper** — immune signature + MAT2B + PNP + disease specificity + Cluster 21
2. **JOSS submission** — Riker Engine methodology (already in progress)
3. **Naviaux communication** — lead with immune signature, mention purine overlap, drop CDR claim

The overnight sessions produced 26 files of rigorous analysis. The original hypothesis didn't survive, but the process found real findings worth publishing. That's good science.

---

*"We tested it properly and got our answer. The cargo is the immune signature and the disease-specific metabolic findings. The vehicle (folate → CDR) gets parked."*
