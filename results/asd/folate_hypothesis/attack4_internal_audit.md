# Attack 4: Internal Consistency Audit

**Date:** 2026-04-17
**Purpose:** Find internal contradictions, temporal plausibility issues, and effect size concerns.
**Stance:** ADVERSARIAL — try to break the hypothesis from inside.

---

## 1. Effect Size Audit — The Numbers Don't All Speak Equally

### KEY FINDING: The metabolic/folate genes have SMALL effect sizes; the immune genes have LARGE ones.

| Gene | Category | Direction | log2FC | Meta RE | Meta p | Assessment |
|------|----------|-----------|--------|---------|--------|------------|
| **C1QB** | Complement | up | **1.242** | N/A | N/A | **Very strong** — >2× fold change |
| **CD14** | Innate immune | up | **1.148** | +0.967 | 4.7e-5 | **Very strong** — large effect, highly significant |
| **NQO1** | Oxidative stress | up | **1.211** | +1.043 | 3.9e-4 | **Very strong** |
| **S100A8** | Alarmin/DAMP | up | **1.087** | N/A | N/A | **Very strong** |
| **IFITM3** | IFN response | up | **1.075** | +0.866 | 1.3e-5 | **Very strong** |
| **VIM** | DAMP/structural | up | **1.055** | +0.956 | 1.9e-2 | **Strong** |
| **C1QC** | Complement | up | **0.983** | N/A | N/A | **Strong** |
| **CXCL16** | Chemokine | up | **0.976** | +0.786 | 3.5e-5 | **Strong** |
| **LYN** | B-cell signaling | up | **0.918** | +0.673 | 1.7e-4 | **Strong** |
| **BST2** | IFN response | up | **0.862** | +0.783 | 1.0e-5 | **Strong** |
| **CYBA** | NADPH oxidase | up | **0.871** | +0.776 | 1.5e-3 | **Strong** |
| **TXNIP** | Inflammasome | up | **0.721** | +0.646 | 2.3e-2 | **Moderate-Strong** |
| **TRAF3IP2** | NF-kB/IL-17 | up | **0.741** | +0.645 | 1.2e-2 | **Moderate-Strong** |
| **B2M** | MHC-I | up | **0.707** | N/A | N/A | **Moderate** |
| **NFE2L2** | NRF2/OxStress | up | **0.592** | N/A | N/A | **Moderate** |
| **HSPD1** | Mito DAMP | up | **0.641** | N/A | N/A | **Moderate** |
| | | | | | | |
| **MAT2B** | One-carbon/SAM | up | **0.130** | N/A | N/A | **WEAK** — tiny fold change |
| **ENTPD6** | Purinergic/CDR | up/down? | **0.045** | -0.225 | 2.1e-3 | **PROBLEMATIC** — direction discrepancy |
| **GLS2** | Cluster 21 | down | **-0.213** | -0.425 | 2.1e-2 | **Moderate** (largest of Cluster 21) |
| **NFS1** | Cluster 21/Fe-S | down | **-0.099** | -0.214 | **0.20** | **WEAK** — meta-analysis NOT significant |
| **UQCRC1** | Cluster 21/CIII | down | **-0.107** | -0.351 | 4.1e-2 | **Weak-Moderate** |
| **SDHA** | Mito/CII | down | **-0.002** | N/A | N/A | **NEGLIGIBLE** — essentially zero |
| **NDUFAF5** | Mito/CI assembly | down | **-0.005** | -0.271 | 1.9e-2 | **NEGLIGIBLE** fold change; meta shows real effect |
| **COX7A1** | Mito/CIV | down | **-0.150** | -0.406 | 3.1e-2 | **Weak** |
| **IDH3A** | TCA cycle | down | **-0.006** | -0.346 | 3.3e-3 | **NEGLIGIBLE** fold change; meta shows effect |

### Critical Assessment

**The immune signal is robust.** 15+ genes with |log2FC| > 0.7, many with highly significant meta-analysis p-values. These are real, large, reproducible effects.

**The metabolic/folate-specific signal is fragile.** The genes that directly connect to the folate hypothesis have concerning effect sizes:

1. **MAT2B** (log2FC = 0.130): This is the centerpiece of the one-carbon metabolism argument. A ~9% increase in expression. While perfectly stable (50/50 runs), the biological significance of such a small change is questionable. Counter-argument: MAT2B is a regulatory subunit, not a catalytic enzyme — small expression changes can have amplified functional effects. But this needs to be flagged.

2. **ENTPD6** (log2FC = +0.045 core, RE = -0.225 meta): **Direction discrepancy.** Core genes say UP, meta-analysis says DOWN. The effect size is tiny either way. This is the ONLY purinergic gene in the iron-clad set and it can't even agree on its direction. This significantly weakens the CDR-purinergic argument.

3. **SDHA, NDUFAF5** (log2FC ≈ 0): These were cited as "methylation-sensitive ETC genes being silenced." But their core gene fold changes are essentially ZERO. The meta-analysis shows real effects for NDUFAF5 (-0.271), suggesting the effect is real but inconsistent across datasets. SDHA's fold change is -0.002 — this is noise, not biology.

4. **NFS1** (meta p = 0.20): Not statistically significant in meta-analysis. The Cluster 21 coherence argument relies partly on NFS1, but its signal is the weakest of the three.

### Verdict: MIXED
- The immune signature is bulletproof on effect sizes
- The metabolic/folate-specific claims are built on small effects and one directionally confused gene
- This doesn't kill the hypothesis but it REFRAMES the strength: the immune arm is the strong evidence, the folate-specific metabolic arm is suggestive at best

---

## 2. Direction Logic Audit — Do the 3 Downregulated Immune Genes Make Sense?

### STAT4 (down, log2FC = -0.758): Th1 transcription factor
- STAT4 drives Th1 differentiation (IFN-γ producing T cells)
- STAT4 down = Th1 suppression
- **Is FRAA-driven pathology Th1 or Th2?** Autoantibody production is classically Th2/Tfh-driven. STAT4 downregulation (less Th1) with simultaneous B-cell/antibody markers upregulated (LYN, PIK3AP1) is CONSISTENT with a Th2-skewed autoimmune response.
- **Verdict: CONSISTENT** — Th1 suppression with Th2/antibody activation fits autoimmune FRAA pathology.

### THEMIS (down, log2FC = -0.315): T-cell selection in thymus
- THEMIS is required for proper T-cell receptor signaling during thymic selection
- THEMIS down = impaired T-cell selection → potentially more autoreactive T cells escaping thymic deletion
- **Verdict: CONSISTENT** — impaired thymic selection could contribute to autoimmune susceptibility. Small effect size noted.

### PELI3 (down, log2FC = -0.389, RE = -0.484, p = 0.011): Pellino E3 ubiquitin ligase in TLR/IL-1R signaling
- PELI3 is a negative regulator of certain TLR signaling branches
- PELI3 down could mean: (a) feedback inhibition of TLR signaling reduced, leading to MORE inflammation, or (b) the specific TLR branch PELI3 regulates is less active
- **Verdict: AMBIGUOUS** — could be consistent (loss of negative regulation = more inflammation) or inconsistent. Not a strong contradiction.

### Overall Direction Audit
**No fatal direction contradictions found.** The 3 downregulated immune genes are either clearly consistent (STAT4) or ambiguous (PELI3). No gene that the hypothesis predicts should be UP is found DOWN, or vice versa, in a way that creates a logical impossibility.

---

## 3. Temporal Plausibility Audit

### The Problem
- Hypothesis claims: developmental insult during 0-36 months (synaptogenesis window)
- Transcriptomic data: from post-mortem brain tissue and blood samples of children/adults
- Question: Would a developmental insult leave a transcriptomic signature years later?

### Assessment
**This is addressed by the CDR model itself.** Naviaux explicitly argues CDR is SUSTAINED — it persists after the initial trigger is removed because the metabolic shift becomes self-reinforcing through:
1. Purinergic signaling feedback loops (extracellular ATP maintains CDR)
2. Epigenetic changes from methylation disruption (persist through cell divisions)
3. Complement-mediated synaptic pruning (structural changes are permanent)
4. Autoantibody persistence (FRAA can persist for years)

**Counter-argument:** If CDR is truly sustained, why does leucovorin supplementation improve symptoms? Frye's trial shows improvement with leucovorin, suggesting the system is still responsive to intervention — which means it's actively maintained, not just a scar.

**Verdict: PLAUSIBLE** — the sustained CDR + persistent FRAA + structural pruning changes argument holds. The leucovorin response data actually supports this (the system is actively maintained, so intervening works). But document the assumption explicitly in any paper.

---

## 4. The ENTPD6 Direction Problem — Detailed Analysis

This is the most concerning internal inconsistency.

**Core genes say:** UP (mean_log2fc = +0.045, cluster 174)
**Meta-analysis says:** DOWN (random_effect = -0.225, p = 0.002)

This happens when:
- Different datasets show different directions
- The core gene analysis (which uses mean across datasets) gets a slightly positive average
- The meta-analysis (which properly models heterogeneity) detects a significant downward trend

**I² for ENTPD6:** Would need to check, but the direction discrepancy suggests high heterogeneity.

**Impact on hypothesis:** ENTPD6 was cited as the CDR-purinergic connection. If its direction is uncertain, the specific CDR-purinergic claim weakens. However:
- The CDR argument now rests MORE on TXNIP (inflammasome bridge, log2FC = 0.721, clear direction) and DAMP signaling (S100A8, HSPD1) than on ENTPD6 alone
- Attack 2 (CDR biomarker cross-reference) will determine if there are other CDR genes beyond ENTPD6

**Verdict: PROBLEMATIC for ENTPD6 specifically, but not fatal for CDR argument if other CDR markers exist.**

---

## 5. NFE2L2 (NRF2) — The Oxidative Stress Elephant

NFE2L2 (NRF2) is the MASTER regulator of the cellular antioxidant response. It is iron-clad and upregulated (log2FC = 0.592). NQO1 (its canonical target) is also iron-clad and strongly upregulated (log2FC = 1.211).

**Problem for the hypothesis:** If NRF2 is upregulated, this could mean oxidative stress is the PRIMARY driver of the mitochondrial signature, with folate deprivation being unnecessary.

**Defense:** NRF2 activation is an expected DOWNSTREAM consequence of mitochondrial dysfunction from any cause, including folate deprivation. Impaired ETC → increased ROS → NRF2 activation is a well-established cascade. NRF2 being present doesn't identify the CAUSE of the oxidative stress — it just confirms oxidative stress IS present.

**But:** If NRF2 is the only "upstream" regulatory change, and folate pathway genes (MTHFR, DHFR, etc.) are absent, Occam's razor says oxidative stress is the simpler explanation. The folate hypothesis adds complexity.

**Verdict: NOT a kill shot, but a serious concern.** Attack 1 will assess this in detail. The hypothesis needs to demonstrate that the oxidative stress is SECONDARY to folate deprivation, not primary. The fact that NRF2 is a downstream responder (not a cause) helps, but the absence of upstream folate pathway genes weakens the causal chain from folate → oxidative stress.

---

## 6. Summary of Findings

### Issues Found (ranked by severity)

| # | Issue | Severity | Impact on Hypothesis |
|---|-------|----------|---------------------|
| 1 | MAT2B has tiny effect size (log2FC = 0.130) | **Serious** | The centerpiece one-carbon finding is biologically small |
| 2 | ENTPD6 has direction discrepancy (up in core, down in meta) | **Serious** | CDR-purinergic connection is unreliable |
| 3 | SDHA, NDUFAF5 have ~zero fold changes in core genes | **Moderate** | "Methylation-sensitive ETC silencing" claim weakened |
| 4 | NFS1 meta-analysis p = 0.20 (not significant) | **Moderate** | Cluster 21 coherence partially undermined |
| 5 | NFE2L2 (NRF2) is iron-clad and upregulated | **Moderate** | Oxidative stress may be primary, not folate |
| 6 | No canonical folate pathway genes in iron-clad set | **Moderate** | The folate hypothesis is entirely indirect |
| 7 | Direction logic for immune genes | **Minor** | All 3 downregulated genes are explainable |
| 8 | Temporal plausibility | **Minor** | CDR sustained model addresses this |

### Does the Hypothesis Survive Attack 4?

**PARTIALLY.** The internal audit reveals that:

1. **The immune arm is strong** — large effect sizes, consistent directions, no contradictions
2. **The metabolic/folate arm is weaker than previously claimed** — MAT2B is real but small, ENTPD6 is confused, several mito genes have negligible fold changes
3. **The correct framing is: "robust neuroinflammatory signature with a suggestive but small folate/methylation connection"** — not "folate deprivation drives everything"

### Recommendation
- Lead with the immune/autoimmune signature (strong effect sizes, unambiguous)
- Present MAT2B as "consistent with" the folate hypothesis rather than "proving" it
- Flag the ENTPD6 direction issue honestly
- Don't overweight SDHA/NDUFAF5 as "methylation-silenced genes" when their fold changes are essentially zero — the meta-analysis effects are real but the core gene data doesn't support dramatic claims
- The hypothesis SURVIVES this audit, but the framing needs to be more modest about the metabolic arm
