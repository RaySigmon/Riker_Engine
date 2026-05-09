# Regime Model — Phase 1 Yield Predicts Discovery vs Characterization Mode

**Date:** 2026-05-08
**Engine version:** 0.3.3.2
**Supporting data:** 7 blind-mode disease days on Pi 5 (AD predicted, untested)

---

## 1. Empirical Observation

When the Riker Engine runs in blind mode (all 19,296 protein-coding genes as input), two qualitatively different outcomes emerge depending on the fraction of genes that pass Phase 1 differential expression filtering:

| Disease | Phase 1 yield | Sig clusters | Iron-clad % | Jaccard | Regime |
|---------|---:|---:|---:|---:|--------|
| T2D | 8.0% | 16 | 67.3%* | 0.851* | Localized |
| ASD | 9.3% | 15 | 86.2% | 0.933 | Localized |
| IPF | 23.2% | 24 | 90.7% | 0.955 | Localized |
| BrCa | 27.5% | 35 | 92.4% | 0.961 | Localized |
| Psoriasis | 41.1% | 0 | 94.3% | 0.970 | Global |
| CRC | 43.0% | 0 | 92.8% | 0.963 | Global |
| IBD | 43.7% | 0 | 92.3% | 0.959 | Global |
| AD | 68.3% | OOM | — | — | Predicted global |

\* T2D metrics reflect small-denominator effect (205 total genes); see FINDINGS_SUMMARY Section 5.2.

**Pattern:** All diseases below ~30% Phase 1 yield produce significant clusters (localized regime). All diseases above ~40% produce zero significant clusters (global regime). The classification has zero misclassifications across 7 tested diseases.

---

## 2. Mechanism

Phase 4 tests each cluster's robustness by permutation: it randomly draws gene sets of the same size as the cluster from the Phase 1 study gene pool, clusters those random sets, and asks whether the real cluster is more coherent than what random gene sets produce.

At low yields, the Phase 1 study gene pool is small and enriched for specific biological programs. Random draws from this pool are unlikely to form coherent clusters — most random subsets span unrelated programs that don't co-cluster. Real biological modules (e.g., ER signaling, cell cycle, immune infiltrate) significantly outperform the permutation null because they have genuine co-expression structure that random subsets lack.

At high yields, the Phase 1 study gene pool is large and includes a broad swath of the transcriptome. Random draws from this pool are themselves enriched for multiple real biological programs — at 43% yield, any random 50-gene draw likely contains fragments of several active pathways. These random subsets form their own clusters with non-trivial coherence. The permutation null rises to meet the real clusters, and no individual cluster achieves significance above the elevated baseline.

This is not a failure of the engine. It is the engine correctly reporting that the disease's transcriptomic signature is diffuse rather than focal. When 40%+ of the protein-coding genome is differentially expressed, the disease is a transcriptome-wide perturbation, not a modular one.

---

## 3. Predictions and Confirmations

The regime model was developed iteratively across disease days. Below are the predictions made before each disease ran and their outcomes.

### 3.1 Confirmed predictions

**T2D (predicted localized, confirmed):** At 8% yield — the lowest in the set — T2D was predicted to produce a small number of tight clusters. Result: 16 significant clusters from 1,545 study genes. The engine identified IAPP (49/50 iron-clad), the canonical beta-cell amyloidosis gene, without any diabetes-specific priors.

**BrCa (predicted localized based on yield, confirmed):** At 27.5% yield, BrCa was predicted to produce localized structure. Result: 35 significant clusters — the highest cluster count of any disease, with independent reconstruction of PAM50 molecular subtypes (ER+, HER2, basal, proliferative, immune, stromal programs). BrCa sits at the upper edge of the localized regime.

**CRC (predicted global based on yield, confirmed):** CRC's 43% yield was predicted to produce 0 significant clusters based on the emerging model. Result: 0 significant clusters, global regime confirmed.

### 3.2 Falsified predictions

**Psoriasis (predicted localized, was global):** Initially predicted to be localized based on disease heterogeneity reasoning (psoriasis has distinct plaque, guttate, and pustular subtypes). At 41.1% yield, Psoriasis produced 0 significant clusters — global regime. This falsified the heterogeneity hypothesis as a primary predictor.

**IBD (predicted localized, was global):** Predicted to be localized (25–35% yield, 20–30 sig clusters) based on IBD's distinct CD/UC subtypes and mucosal heterogeneity. At 43.7% yield, IBD produced 0 significant clusters — global regime. This was the second falsification of the heterogeneity hypothesis.

### 3.3 Lesson from falsifications

Both falsified predictions used disease-type reasoning ("this disease has subtypes, so it should produce distinct modules") rather than yield-based reasoning. Both times, the data showed yield matters more than heterogeneity. **Phase 1 yield is the single best predictor.** Disease heterogeneity is irrelevant once yield is known.

### 3.4 Untested prediction

**AD (predicted global):** At 68.3% Phase 1 yield — the highest in the set — AD is predicted to produce 0 significant clusters with iron-clad fraction 93–95% and Jaccard 0.96–0.97. AD blind OOM'd on Pi 5 and will be tested on RunPod. If AD blind produces significant clusters despite 68% yield, the regime model needs revision.

---

## 4. The Transition Zone

The regime flips between 27.5% (BrCa, localized, 35 sig clusters) and 41.1% (Psoriasis, global, 0 sig clusters). No tested disease falls in this gap.

This is an honest limitation. The transition could be:
- **Sharp threshold** (~30–35%): a single critical yield above which the permutation null saturates
- **Gradual transition** (30–40%): a range where cluster count decreases from many to zero as yield increases
- **Disease-specific**: some diseases in the gap might be localized, others global, depending on the structure of their differential expression

Without data points between 27.5% and 41.1%, we cannot distinguish these hypotheses. The methods paper should state this as "a sharp transition observed between approximately 28% and 41% Phase 1 yield" without claiming a specific threshold value.

Future work: deliberately selecting diseases with expected yields in the 30–40% range would resolve the transition shape. Alternatively, adjusting the Phase 1 p-value threshold (e.g., p < 0.01 instead of p < 0.05) on existing diseases would shift yields downward and provide additional data points.

---

## 5. Implications for Use

### 5.1 Localized regime (yield < ~30%)

The blind tier produces a credible candidate gene set. Significant clusters correspond to distinct biological programs. Iron-clad genes from 50-run stability profiling are high-confidence candidates for follow-up.

**Appropriate claims:** "The engine identified N gene modules associated with [disease], including [specific rediscovery] at [stability] iron-clad status."

**What to report:** Curated-tier validation (does the engine reproduce known biology?), blind-tier discovery (what new candidates emerge?), stability profile (which candidates are robust?), cluster content analysis (what biological programs do they represent?).

### 5.2 Global regime (yield > ~40%)

The blind tier produces a transcriptomic-scale measurement, not a candidate list. Zero significant clusters means no individual gene module is more robust than what random gene subsets produce at this scale.

**Appropriate claims:** "The engine identified N genes as consistently differentially expressed across K independent cohorts, but the signal is transcriptome-wide rather than concentrated in specific modules."

**What to report:** Curated-tier validation remains informative (the curated seed list restricts to known biology and produces meaningful clusters). Blind-tier gene count characterizes the disease's transcriptomic breadth. The curated result is the publishable finding; the blind result contextualizes it.

### 5.3 Decision protocol for new diseases

1. Run curated tier first (validates known biology, fast)
2. Run blind tier (measures Phase 1 yield)
3. If yield < 30%: proceed to stability profiling and cluster content analysis
4. If yield > 40%: report curated tier as primary finding, blind tier as characterization
5. If yield 30–40%: inspect cluster count — may be either regime

---

## 6. Limitations and Open Questions

### 6.1 Hardware environment

All 7 data points come from a single hardware environment (Raspberry Pi 5, 8GB RAM). Cross-hardware reproducibility will be tested by re-running one disease on RunPod (cloud compute). If results differ, the regime model needs a hardware-dependence caveat.

### 6.2 AD prediction untested

AD's 68% yield is the highest in the set and has not been tested due to OOM constraints. The prediction (0 sig clusters, high stability) is confident but empirically unverified. Cloud compute will resolve this.

### 6.3 Transition zone uncharacterized

No data between 27.5% and 41.1% yield. The transition shape is unknown. Additional diseases in this yield range or Phase 1 threshold adjustments on existing diseases would resolve it.

### 6.4 Data type scope

All validation uses bulk transcriptomics (microarray + RNA-seq) from GEO. Whether the regime model generalizes to:
- Single-cell RNA-seq
- Proteomics
- Methylation arrays
- Other omics platforms

is untested. The mechanism (permutation null saturation at high yields) should apply to any feature-selection-then-clustering pipeline, but empirical verification is needed.

### 6.5 Phase 1 threshold sensitivity

All runs use p < 0.05, min 2 datasets. Different thresholds would produce different yields for the same disease. A stricter threshold (p < 0.01) would lower yields and potentially move global-regime diseases into the localized regime. Whether the regime model's biological conclusions change under different thresholds is an open question.

### 6.6 Regime as disease property vs analysis property

The current framing treats regime as a property of the disease ("IBD is a global-regime disease"). An alternative framing: regime is a property of the analysis configuration (disease × datasets × threshold × seed list). The same disease could be localized under one analysis and global under another. This distinction matters for methods paper framing — prefer "IBD under this analysis configuration falls in the global regime" over "IBD is inherently a global-regime disease."

---

## References

- FINDINGS_SUMMARY.md — full validation data (Section 3 for regime data, Section 5 for stability)
- Disease day manifests in `disease_days/2026-05-*/DISEASE_DAY_MANIFEST.md`
- Stability profiles in `disease_days/2026-05-*/stability_50run/stability_summary.json`
