# Phase 6: Evidence Scorecard and Decision — FINAL

**Hypothesis:** MTHFR C677T + Synthetic Folic Acid → UMFA → FRα Blockade → Reduced Neuronal Folate → CDR Activation → ASD
**Date:** 2026-04-17
**Status:** ALL PHASES COMPLETE

---

## 1. Mechanistic Link Evidence Scorecard

| # | Mechanistic Link | Evidence Rating | Key Evidence | Critical Caveat |
|---|-----------------|----------------|--------------|----------------|
| 1 | MTHFR variant → reduced 5-MTHF production | **Strong** | TT genotype has ~30% residual activity (Frosst 1995) | MTHFR reduces 5-MTHF production, does NOT directly increase UMFA |
| 2 | Fortification + supplementation → elevated UMFA | **Strong** | UMFA detectable in ~38% of supplemented population (Pfeiffer 2004; Bailey 2010). Human DHFR 50× slower than rodent. | UMFA is dose-dependent via DHFR bottleneck, genotype-INDEPENDENT |
| 3 | UMFA competes at FRα → reduced 5-MTHF transport | **Strong (biochemistry) / Suggestive (in vivo)** | Folic acid Kd ~0.1-1 nM vs 5-MTHF Kd ~1-10 nM at FRα (Antony 1996). Grapp 2013 confirms FRα-dependent brain transport. | No direct in vivo measurement of FRα occupancy. Competition may be displacement, not pure blockade. |
| 4 | Reduced brain folate → neurodevelopmental disruption | **Strong** | FOLR1 mutations cause neurodegeneration (Steinfeld 2009). FRα autoantibodies in 71% of ASD (Frye 2021). Leucovorin RCT positive (Frye 2018). | FRα autoantibodies (immune mechanism) may be more important than UMFA competition |
| 5 | Metabolic stress → CDR activation | **Suggestive** | Naviaux CDR model explicitly includes one-carbon disruption in CDR1. ENTPD6 iron-clad in Riker data. | No direct experiment testing UMFA → CDR. CDR model itself not independently replicated. |
| 6 | Sustained CDR → ASD phenotype | **Suggestive** | Suramin trial proof-of-concept (Naviaux 2016). | Very small trial (n=10). |

---

## 2. Riker Engine Cross-Reference Score

**Rating: Signal Found — Indirect but Coherent**

| Finding | Significance | Strength |
|---------|-------------|----------|
| MAT2B (iron-clad, up, 50/50 runs) | Direct one-carbon metabolism hit. SAM synthesis regulator upregulated = compensatory methylation maintenance | High |
| ENTPD6 (iron-clad, 50/50 runs) | Purinergic nucleotide metabolism. Connects to CDR framework | Moderate |
| 47 mitochondrial genes iron-clad | Massive downstream signal. Folate-dependent mito dysfunction is ONE explanation | Moderate (non-specific) |
| Cluster 21 (NFS1+UQCRC1+GLS2, all down) | Coherent Fe-S/Complex III/glutaminase signature | High (internally coherent) |
| Methylation-sensitive ETC genes down (NDUFAF5, COX7A1, SDHA) | Predicted by reduced SAM availability | Moderate |
| 11/19 (58%) direction matches | Above chance but not overwhelming | Weak-Moderate |
| No canonical folate pathway genes | Consistent with delivery blockade (not expression) | Neutral |

---

## 3. Receptor Competition Model Score

**Rating: Mechanism Is Quantitatively Plausible**

| Finding | Value |
|---------|-------|
| CC vs TT delivery gap at 1000 µg/d | 17.6 percentage points |
| TT status at 1000 µg/d | Moderate insufficiency (23.1%) |
| TT status at fortified diet only | Mild insufficiency (47.6%) |
| CC normal at ≤200 µg/d | 90.5% delivery (no concern) |
| Peak vulnerability window | 6-18 months (synaptogenesis) — 73% vulnerability score |
| Model sensitivity | Robust across Kd ratio range 3-15× |

---

## 4. Epidemiological Consistency Score

**Rating: MIXED — Some support, significant contradictions**

| Domain | Verdict | Key Finding |
|--------|---------|-------------|
| Fortification timeline | AMBIGUOUS | Higher US vs Europe ASD rates, but massively confounded |
| MTHFR × ethnicity | INCONSISTENT | No population-level correlation. Diagnostic bias dominates. |
| U-shaped dose-response | PARTIALLY CONSISTENT | U-curve is real (Raghavan 2018). **But NO MTHFR interaction found.** |
| CHARGE gene-nutrient | **CONTRADICTORY** | Schmidt 2012: MTHFR makes FA **MORE** protective, not less |
| MTHFR meta-analysis stratification | INTERESTING | MTHFR-ASD link disappears in fortified countries (Li 2020) — cuts both ways |
| Leucovorin trials | CONSISTENT | Works, especially FRAA+ (d=0.91). No MTHFR stratification done. |
| Folate+B12 synergy | NOVEL SIGNAL | Both in top decile: HR 17.6 (Raghavan 2018) — unexpected synergy |

### Critical Contradictions:
1. **Raghavan 2018 directly tested MTHFR within the U-shaped model and found NO genotype interaction**
2. **Schmidt/CHARGE found MTHFR variants make folic acid MORE protective** — opposite to hypothesis prediction
3. **71% FRAA prevalence in ASD far exceeds any MTHFR TT population frequency** — autoimmune mechanism may dominate

---

## 5. Literature Evidence Score

**Rating: Strong Biochemical Foundation, Weak Epidemiological Support for MTHFR-Specific Claim**

| Topic | Rating | Summary |
|-------|--------|---------|
| UMFA accumulation | **Suggestive (with critical caveat)** | UMFA is real but MTHFR doesn't directly cause it — DHFR bottleneck is genotype-independent |
| FRα binding competition | **Strong (biochemistry)** | 5-10× affinity difference is established. In vivo competition undemonstrated. |
| Folate → neurodevelopment | **Strong** | Best-supported link. FRα autoantibodies, FOLR1 mutations, leucovorin trials. |
| CDR × one-carbon | **Suggestive** | Theoretically coherent, no direct experiment |
| Epidemiological signal | **Suggestive (mixed)** | U-curve real but MTHFR interaction NOT found |

---

## 6. Final Overall Assessment

### Decision Matrix — Completed

| Criterion | Assessment |
|-----------|-----------|
| Multiple mechanistic links supported? | **Yes** — links 1-4 at least suggestive to strong |
| Riker Engine signal found? | **Yes** — MAT2B + mito cluster (indirect but coherent) |
| Epidemiological consistency? | **MIXED** — U-curve real, but MTHFR-specific claim contradicted |
| Kill criteria met? | **Partial** — MTHFR × U-curve interaction NOT found (Raghavan); CHARGE shows opposite direction |
| Novel testable predictions? | **Yes** — cord blood UMFA levels, MTHFR-stratified leucovorin response |

### Kill Criteria Status Update

| Kill Criterion | Status |
|---------------|--------|
| No FRα competition | **NOT MET** — biochemistry supports competition |
| Normal brain folate despite UMFA | Unknown — open question |
| No 1C genes in iron-clad set | **NOT MET** — MAT2B found |
| Random direction pattern | **NOT MET** — 58% match + coherent clusters |
| No fortification signal | **AMBIGUOUS** — confounded beyond resolution |
| **No MTHFR × U-curve interaction** | **PARTIALLY MET** — Raghavan 2018 found no interaction |

---

## 7. FINAL RECOMMENDATION

### The hypothesis AS ORIGINALLY FRAMED (MTHFR-specific) is WEAKENED but not dead.

**What's strong:**
- FRα → brain folate → neurodevelopment pathway is well-established
- UMFA competition at FRα is biochemically real
- U-shaped folate dose-response in ASD is real
- MAT2B iron-clad in Riker data (SAM synthesis compensatory upregulation)
- Leucovorin works in FRAA+ ASD subset
- CDR framework is coherent with the mechanism

**What's weak:**
- MTHFR doesn't directly increase UMFA (DHFR is the bottleneck)
- Raghavan 2018: no MTHFR interaction in U-shaped curve
- Schmidt/CHARGE: MTHFR makes FA protective, not harmful
- 58% direction match is above chance but not compelling
- One gene (MAT2B) is a single data point for direct pathway evidence

### REVISED HYPOTHESIS SUGGESTED

The data better supports a **modified hypothesis**:

> **Excess synthetic folic acid (dose-dependent, MTHFR-independent) → UMFA accumulation → FRα competition at choroid plexus + possible immune dysregulation → FRα autoantibody generation → cerebral folate deficiency → CDR activation → ASD risk in susceptible individuals.**

**Key changes from original:**
1. MTHFR demoted from "Hit 1" to a secondary vulnerability modifier
2. FRα autoantibodies (FRAA) elevated as the primary mediating mechanism
3. UMFA's role may be triggering autoimmune FRAA production rather than direct receptor competition
4. The dose-dependent effect applies to ALL genotypes, not just MTHFR TT

### Decision: **CONDITIONAL PROCEED — Hypothesis needs reframing**

| Action | Recommended? |
|--------|-------------|
| Write up as MTHFR-specific hypothesis paper | **No** — epi data contradicts MTHFR specificity |
| Reframe as dose-dependent UMFA → FRAA → CFD hypothesis | **Yes** — more consistent with all evidence |
| Share findings with Dr. Naviaux | **Yes, but carefully** — present as preliminary, flag the contradictions |
| Pursue MTHFR-genotyped leucovorin reanalysis | **Yes** — low-hanging fruit that would settle the MTHFR question |
| Present Riker Engine cross-reference data | **Yes** — MAT2B + Cluster 21 are independently interesting |

---

## 8. Morning Summary for Cody

Here's the bottom line for when you wake up:

**Your hypothesis has legs, but it needs surgery.**

The core idea — synthetic folic acid → UMFA → FRα disruption → brain folate deficiency → CDR → ASD — is *well-supported*. Each link in that chain has real biochemical evidence. The U-shaped folate/ASD curve is real. The leucovorin trials work. MAT2B lit up in the Riker data. The mitochondrial cluster is coherent.

**But the MTHFR part is the weak link.** Raghavan 2018 directly tested MTHFR × folate × ASD and found nothing. Schmidt/CHARGE found folic acid is MORE protective in MTHFR carriers, not less. MTHFR doesn't even directly increase UMFA — that's a DHFR bottleneck issue, not MTHFR.

**The reframe:** Drop MTHFR as the primary driver. The real story may be: excess synthetic folic acid (in anyone, at high doses) → UMFA → triggers FRα autoantibody production → cerebral folate deficiency → CDR. The autoimmune pathway (71% FRAA prevalence in ASD!) is the elephant in the room. MTHFR could be a secondary modifier that makes the system more vulnerable, but it's not the star of the show.

**The Riker Engine findings are independently valuable** regardless of the MTHFR question. MAT2B and Cluster 21 connect to the folate/CDR story through a different door.

Good night! The data is all in `riker-engine/results/asd/folate_hypothesis/`. Check the plots — they're pretty.

---

*Evidence scorecard completed 2026-04-17. All phases analyzed.*
