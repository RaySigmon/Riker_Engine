# Attack 2: CDR Biomarker Cross-Reference

## Question
Does the Riker ASD iron-clad gene set show a specific Cell Danger Response (CDR) signature, or is the CDR framing an overreach based on a single purinergic gene (ENTPD6)?

## Method
Cross-referenced all 376 iron-clad genes (appearance rate >= 0.90, 50-run blind stability) against 9 CDR pathway categories derived from Naviaux 2014/2018.

---

## Results: CDR Gene Hits in Iron-Clad Set

### Category 1: Purinergic Signaling (CDR Core)
| Gene | Stability | Cluster | Direction | log2FC |
|------|-----------|---------|-----------|--------|
| ENTPD6 | 50/50 | 174 | up | +0.045 |

**1 hit out of 35 canonical genes.** No P2RX, P2RY, ADORA, PANX, NT5E, or ADA genes detected.

### Category 2: Sphingolipid/Ceramide Metabolism (CDR1)
| Gene | Stability | Cluster | Direction | log2FC |
|------|-----------|---------|-----------|--------|
| — | — | — | — | — |

**0 hits out of 28 canonical genes.**

Note: B4GALT6 (iron-clad, 48/50, cluster 156, down, log2fc=-0.317) is a lactosylceramide synthase functionally adjacent to this pathway but is NOT in Naviaux's canonical CDR sphingolipid gene list.

### Category 3: Phospholipid Remodeling (CDR1/2)
| Gene | Stability | Cluster | Direction | log2FC |
|------|-----------|---------|-----------|--------|
| LPCAT4 | 50/50 | 177 | down | -0.087 |

**1 hit out of 23 canonical genes.**

### Category 4: Ganglioside Metabolism (CDR2)
**0 hits out of 17 canonical genes.**

### Category 5: Mitochondrial Dynamics/Fission-Fusion (CDR)
**0 hits out of 15 canonical genes.** No DRP1, MFN1/2, OPA1, FIS1, PINK1, or Parkin.

### Category 6: DAMP/Innate Danger Signaling (CDR Activation)
| Gene | Stability | Cluster | Direction | log2FC |
|------|-----------|---------|-----------|--------|
| S100A8 | 50/50 | 34 | up | +1.087 |

**1 hit out of 19 canonical genes.** S100A8 is a strong hit (high effect size, calprotectin subunit, classic DAMP). However, no TLR2/4/9, no NLRP3, no cGAS-STING, no IL1B/IL18.

Note: S100A10 (iron-clad, 50/50, cluster 34, up, log2fc=+0.974) is in the same cluster but is NOT in Naviaux's canonical DAMP list. It is an annexin A2 ligand with distinct function.

### Category 7: Tryptophan/Kynurenine Pathway (CDR-Related)
| Gene | Stability | Cluster | Direction | log2FC |
|------|-----------|---------|-----------|--------|
| KYAT3 | 50/50 | 138 | down | -0.137 |

**1 hit out of 11 canonical genes.** KYAT3 (kynurenine aminotransferase 3) converts kynurenine to kynurenic acid (neuroprotective branch). Downregulation shifts balance toward neurotoxic metabolites. No IDO1/2, TDO2, or KMO detected.

### Category 8: Methyl/One-Carbon Metabolism (CDR1)
| Gene | Stability | Cluster | Direction | log2FC |
|------|-----------|---------|-----------|--------|
| MAT2B | 50/50 | 131 | up | +0.130 |
| ENTPD6 | 50/50 | 174 | up | +0.045 |

**2 hits** (already known from Phase 2 analysis).

### Category 9: Eicosanoid/Oxylipid Signaling (CDR)
**0 hits out of 12 canonical genes.** No COX1/2, no lipoxygenases.

---

## Summary Table

| CDR Category | Canonical Genes | Hits | Hit Rate |
|---|---|---|---|
| Purinergic signaling | 35 | 1 (ENTPD6) | 2.9% |
| Sphingolipid/ceramide | 28 | 0 | 0% |
| Phospholipid remodeling | 23 | 1 (LPCAT4) | 4.3% |
| Ganglioside metabolism | 17 | 0 | 0% |
| Mito dynamics/fission-fusion | 15 | 0 | 0% |
| DAMP/danger signaling | 19 | 1 (S100A8) | 5.3% |
| Tryptophan/kynurenine | 11 | 1 (KYAT3) | 9.1% |
| Methyl/one-carbon | counted separately | 2 (MAT2B, ENTPD6) | — |
| Eicosanoid/oxylipid | 12 | 0 | 0% |
| **TOTAL** | **~160 unique** | **5 unique genes** | **~3.1%** |

---

## CDR Phase Assessment

### Does this match CDR1 (oxidative shielding, innate activation)?
- S100A8 up: YES (DAMP alarm signal)
- MAT2B up: PARTIAL (methyl stress)
- KYAT3 down: PARTIAL (kynurenine shift)
- Missing: sphingolipid changes, ceramide accumulation, full inflammasome, eicosanoids

### Does this match CDR2 (membrane remodeling, proliferative)?
- LPCAT4 down: WEAK (one phospholipid gene, small effect)
- Missing: ganglioside changes, full phospholipid remodeling cascade

### Does this match CDR3 (neuroendocrine normalization)?
- No evidence

---

## Verdict

### Kill condition evaluation
The kill condition was: "If ENTPD6 is the ONLY CDR gene and no other categories show signal, CDR framing is overreach."

**ENTPD6 is NOT the only CDR gene.** Four additional canonical CDR genes are iron-clad: S100A8, LPCAT4, KYAT3, and MAT2B. This spans 4 of 9 CDR categories.

### Survive condition evaluation
The survive condition was: "Multiple CDR categories show hits, especially sphingolipid/DAMP/purinergic, with a pattern consistent with sustained CDR1."

**PARTIAL SURVIVE.** Multiple categories are hit (4/9), but:
- Sphingolipid category has ZERO canonical hits
- Only 1 DAMP gene (S100A8), no inflammasome or TLR genes
- Only 1 purinergic gene (ENTPD6), with tiny effect size
- 5 genes out of ~160 canonical CDR genes = 3.1% hit rate

### Final Assessment: CDR-ADJACENT, NOT CDR-SPECIFIC

The Riker ASD data shows **scattered CDR-adjacent signal** rather than a **coherent CDR signature**. The pattern is:

1. **S100A8 is real and strong** — this is a genuine danger signal (log2fc=+1.09, 50/50 stability). But S100A8 upregulation occurs in many inflammatory contexts (infection, autoimmunity, neuroinflammation). It is not CDR-specific.

2. **KYAT3 downregulation is interesting** — shifts kynurenine toward neurotoxic branch. This is consistent with neuroinflammation broadly, not CDR specifically.

3. **LPCAT4 and ENTPD6 are weak hits** — small effect sizes (log2fc < 0.1), and each is a single gene from a large pathway family. One LPCAT out of 4 and one ENTPD out of 8 does not constitute pathway activation.

4. **Critical absences destroy the CDR narrative:**
   - Zero sphingolipid/ceramide genes (CDR1 hallmark)
   - Zero mitochondrial fission/fusion genes (CDR hallmark)
   - Zero purinergic receptors (P2X/P2Y — CDR core mechanism)
   - Zero inflammasome components (NLRP3/CASP1/IL1B)
   - Zero eicosanoid genes

**The honest framing:** The Riker ASD data shows neuroinflammation (S100A8, immune clusters), mitochondrial dysfunction (OXPHOS genes, TCA cycle), and metabolic stress (MAT2B, KYAT3). These are features that OVERLAP with CDR but do not constitute a CDR-specific signature. Calling this "CDR activation" based on 5 scattered genes from a 160-gene framework is an overreach.

### Recommended language for the paper
- OK: "Several genes overlap with CDR-associated pathways (S100A8, KYAT3, ENTPD6), suggesting partial convergence with Naviaux's cell danger model"
- NOT OK: "The gene signature is consistent with CDR1 activation"
- NOT OK: "These findings support the CDR model of autism"
