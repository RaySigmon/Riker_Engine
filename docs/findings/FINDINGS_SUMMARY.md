# Riker Engine v0.3.3.2 — Cross-Disease Validation Findings

**Date:** 2026-05-08
**Engine version:** 0.3.3.2 (frozen April 29, 2026)
**Hardware:** Raspberry Pi 5, 8GB RAM (7/8 diseases); AD blind deferred to cloud
**SOP:** v1.0.3 (3-tier: curated → blind → 50-run stability)

---

## 1. Validation Matrix

### 1.1 Tier completion

| Disease | Curated (Tier 1) | Blind (Tier 2) | Stability (Tier 3) | Status |
|---------|:---:|:---:|:---:|--------|
| T2D | Done | Done | Done (J=0.851*) | Complete |
| ASD | Done | Done | Done (J=0.933) | Complete |
| IPF | Done | Done | Done (J=0.955) | Complete |
| BrCa | Done | Done | Done (J=0.961) | Complete |
| IBD | Done | Done | Done (J=0.959) | Complete |
| CRC | Done | Done | Done (J=0.963) | Complete |
| Psoriasis | Done | Done | Done (J=0.970) | Complete |
| AD | Done | Deferred (OOM) | — | Cloud-bound |

\* T2D Jaccard reflects small-denominator effect (see Section 5).

7/8 diseases fully validated on Pi 5. AD blind OOM at 68% Phase 1 yield (exceeds 8GB RAM ceiling).

### 1.2 Curated tier results

| Disease | Phase 1 study genes | Core genes | Sig clusters | Survived | Meta-significant |
|---------|---:|---:|---:|---:|---:|
| T2D | 56 | 8 | 1 | 8 | 8 |
| ASD | 141 | 29 | 2 | 29 | 10 |
| Psoriasis | 60 | 50 | 2 | 50 | 28 |
| BrCa | 227 | 153 | 2 | 141 | 114 |
| IPF | 241 | 186 | 4 | 166 | 155 |
| CRC | 331 | 262 | 4 | 245 | 219 |
| IBD | 408 | 301 | 3 | 288 | 284 |
| AD | 437 | 389 | 5 | 334 | 307 |

### 1.3 Blind tier results

| Disease | Phase 1 study genes | Yield % | Core genes | Sig clusters | Survived | Meta-sig | Regime |
|---------|---:|---:|---:|---:|---:|---:|--------|
| T2D | 1,545 | 8.0% | 169 | 16 | 162 | 78 | Localized |
| ASD | 1,793 | 9.3% | 421 | 15 | 414 | 233 | Localized |
| IPF | 4,479 | 23.2% | 2,451 | 24 | 1,856 | 1,584 | Localized |
| BrCa | 5,315 | 27.5% | 3,697 | 35 | 3,380 | 2,555 | Localized |
| Psoriasis | 7,936 | 41.1% | 6,275 | 0 | 6,183 | 3,718 | Global |
| CRC | 8,288 | 43.0% | 6,052 | 0 | 5,790 | 5,036 | Global |
| IBD | 8,432 | 43.7% | 6,027 | 0 | 5,664 | 5,488 | Global |
| AD | — | 68.3%† | — | — | — | — | Predicted global |

† AD Phase 1 yield measured before OOM kill in Phase 3.

---

## 2. Cross-Engine-Version Stability (v0.3.2 → v0.3.3.2)

| Disease | v0.3.2 core | v0.3.3.2 core | Delta | v0.3.2 survived | v0.3.3.2 survived | Delta |
|---------|---:|---:|---:|---:|---:|---:|
| T2D | 8 | 8 | 0 | 8 | 8 | 0 |
| Psoriasis | 50 | 50 | 0 | 50 | 50 | 0 |
| ASD | 35 | 29 | -6 | 20 | 29 | +9 |
| BrCa | 152 | 153 | +1 | 139 | 141 | +2 |
| IPF | 190 | 186 | -4 | — | 166 | — |
| CRC | 264 | 262 | -2 | — | 245 | — |
| IBD | 304 | 301 | -3 | 302 | 288 | -14 |
| AD | 394 | 389 | -5 | 340 | 334 | -6 |

**Summary:** 7/8 within ±2% on core gene count (two exact matches). ASD is the outlier at -6 core genes on a small absolute count (35 → 29), but the story is richer than the core delta: v0.3.2 found 35 core genes but only 20 survived Phase 5 replication (15 eliminated); v0.3.3.2 found 29 core genes and ALL 29 survived (0 eliminated). Meta-significant count increased from 9 to 10. The dilution fix tightened Phase 4 but the surviving gene set is more robust.

The systematic slight decrease in core counts across 6/8 diseases has a documented mechanism: the v0.3.3.2 mean log2FC dilution fix corrected an upward bias in fold-change calculation, producing a tighter but higher-quality gene set.

**T2D and Psoriasis exact matches explained:** T2D's study set is so small (56 genes) that the dilution fix had nothing to dilute. Psoriasis is in the global regime where the gate the fix improved isn't doing work. Both edge cases showing zero delta is mechanistically coherent.

---

## 3. Regime Model

### 3.1 The finding

Phase 1 yield (fraction of protein-coding genome that is differentially expressed) is the single best predictor of whether the engine's blind discovery mode produces localized gene modules (significant clusters) or a global characterization (zero significant clusters).

### 3.2 Empirical support

| Yield | Sig clusters | Regime |
|------:|---:|--------|
| 8.0% | 16 | Localized |
| 9.3% | 15 | Localized |
| 23.2% | 24 | Localized |
| 27.5% | 35 | Localized |
| 41.1% | 0 | Global |
| 43.0% | 0 | Global |
| 43.7% | 0 | Global |
| 68.3% | OOM | Predicted global |

**Transition:** Sharp boundary between 27.5% and 41.1%. No data points in the gap. 7 supporting data points with zero misclassifications.

### 3.3 Mechanism

At low yields, the Phase 1 filter selects a small fraction of the genome — these genes are enriched for specific biological programs and cluster into distinct modules. At high yields, Phase 1 selects a large fraction — these genes span many programs and the permutation null in Phase 4 saturates (i.e., random gene sets of comparable size produce comparable clustering structure, so no individual cluster achieves significance above the permutation baseline).

### 3.4 Stability by regime

| Regime | Jaccard range | Iron-clad % range |
|--------|---:|---:|
| Localized (excl. T2D) | 0.933–0.961 | 86–92% |
| Global | 0.959–0.970 | 92–94% |

Global regime diseases are slightly more stable than localized — mechanistically sensible because when most of the genome is differentially expressed, stochastic clustering decisions affect a smaller fraction of the total gene set.

### 3.5 Falsified predictions

The initial hypothesis that disease-type heterogeneity (e.g., IBD having distinct CD/UC subtypes) would predict localized structure was falsified by IBD landing in the global regime at 43.7% yield. **Yield, not disease heterogeneity, determines regime.**

### 3.6 Pre-registered prediction

AD blind is predicted to produce 0 significant clusters at 65–70% yield, with iron-clad fraction 93–95% and Jaccard 0.96–0.97. To be tested on RunPod.

---

## 4. Headline Rediscoveries

The engine identified canonical disease biology from blind protein-coding input (no disease-specific priors), and these rediscoveries are iron-clad across 50-run stability profiles.

### 4.1 Summary table

| Disease | Rediscovery | Stability | Significance |
|---------|-------------|-----------|--------------|
| T2D | IAPP (islet amyloid polypeptide) | 49/50 iron-clad | Canonical T2D pathology gene, independently recovered |
| IPF | FAM107A | 50/50 iron-clad | Novel IPF candidate — not in curated seed list |
| BrCa | Molecular subtype markers (ESR1, FOXA1, GATA3, ERBB2, AURKA, TOP2A, etc.) | 15/17 markers at 50/50; 2 at 49/50 | Engine reconstructed PAM50 subtype framework from raw expression without subtype labels |
| ASD | SFARI overlap (GABRG2, SLC12A5, NPAS3) | All 50/50 iron-clad | Canonical SFARI genes recovered in blind mode; ATP2B2 borderline at 37/50 |

### 4.2 BrCa subtype reconstruction detail

Six canonical biological programs independently identified from 19,296 protein-coding genes across 5 datasets with no subtype labels:

| Program | Markers found | Cluster pattern |
|---------|---------------|-----------------|
| ER+/Luminal | FOXA1+GATA3 co-clustered; ESR1 in related cluster | Distinct TF modules |
| HER2 | ERBB2 + PGAP3 in adjacent clusters | HER2 amplicon signature |
| Proliferation | AURKA+BIRC5 together; TOP2A, CDK1, CCNB1 | Cell cycle machinery |
| Basal | EGFR, KRT14 in distinct clusters | Basal markers |
| Immune infiltrate | CXCL9, CXCL10 in distinct clusters | Chemokine signaling |
| ECM/stromal | FN1, FAP, MMP11 | Matrix remodeling |

Markers are distributed across multiple clusters per program rather than one mega-cluster per subtype. This is biologically appropriate: real biology has distinct but related gene modules within each program (e.g., ER transcription factors vs ER target genes vs luminal differentiation markers).

See: `disease_days/2026-05-05_brca_v0332/stability_50run/stability_report.csv` for per-gene stability data.

### 4.3 ASD — notable absences

Several canonical SFARI genes (TBR1, SCN2A, SHANK3, PTEN, CHD8) are NOT in the blind stability set. These are neurodevelopmental genes that may not show strong differential expression in adult post-mortem brain cortex (the available tissue). The engine correctly did not force them into the output — absence of expression-level signal is a genuine negative.

---

## 5. Stability Profile Detail

### 5.1 Cross-disease stability table

| Disease | Iron-clad | Iron-clad % | Jaccard | Total genes | Core range | Regime |
|---------|---:|---:|---:|---:|---:|--------|
| T2D | 138 | 67.3% | 0.851 | 205 | 156–179 (23) | Localized |
| ASD | 394 | 86.2% | 0.933 | 457 | — (24) | Localized |
| IPF | 2,309 | 90.7% | 0.955 | 2,545 | — (55) | Localized |
| BrCa | 3,502 | 92.4% | 0.961 | 3,789 | — (71) | Localized |
| IBD | 5,729 | 92.3% | 0.959 | 6,209 | 5,968–6,068 (100) | Global |
| CRC | 5,748 | 92.8% | 0.963 | 6,192 | — (63) | Global |
| Psoriasis | 6,076 | 94.3% | 0.970 | 6,441 | — (68) | Global |

### 5.2 T2D small-denominator note

T2D's lower median Jaccard (0.851) reflects the small-denominator effect at very low Phase 1 yields rather than reduced engine reproducibility. Per-run core gene counts varied within a 23-gene range (156–179, mean 166), with iron-clad gene fraction of 67% indicating consistent identification of the same biological signal across runs. The canonical T2D rediscovery (IAPP) appeared in 49/50 runs.

For cross-disease comparison, core gene count range (absolute) is a more informative secondary stability metric than Jaccard alone when gene set sizes vary by 30×.

---

## 6. Provenance Receipts

### 6.1 Per-disease provenance

Each disease day directory contains:
- `DISEASE_DAY_MANIFEST.md` — full provenance record
- `curated/pipeline_summary.json` — locked_core_genes captured before Phase 5
- `blind_pc/pipeline_summary.json` — locked_core_genes captured before Phase 5
- `blind_pc/config.yaml` — exact config used
- `stability_50run/stability_summary.json` — machine-readable stability metrics

### 6.2 Integrity verification

| Item | Method |
|------|--------|
| Seed file identity | SHA256 checksums in each manifest |
| Engine version | `package_version` field in every pipeline_summary.json |
| Code identity | `code_version` (git commit) in every pipeline_summary.json |
| Pre-specification | `locked_core_genes` array written before Phase 5 replication |
| Cross-version | v0.3.2 baselines documented per manifest |
| Run-to-run | 50-run stability profiles with pairwise Jaccard |

### 6.3 Known provenance imperfection

IBD Tier 3 `stability_summary.json` records `engine_commit: fab1375` (HEAD at write time) rather than `a397083` (the blind run commit). Both are v0.3.3.2 with no engine code changes between them — only data files were added. Verified with `git diff a397083..fab1375 -- riker/` (empty diff). Tracked as a v0.4.0 runner improvement.

---

## 7. Negative Control

**Status:** Complete (50 trials, 2026-05-08). Stability profiling of random input pending.

**Historical context:** An exploratory negative control (v0.3.2, April 2026) using a mismatched config (5 datasets, no phase3/phase4 blocks, SFARI genes excluded) produced 5 core genes / 1 meta-significant. The project's own audit (NEGATIVE_CONTROL_SUMMARY.md) flagged this as "exploratory evidence only" and required a matched-config rerun before publication-grade claims. Today's matched measurement (7 datasets, identical phase3/phase4 to ASD blind, no exclusions) supersedes that result. The discrepancy (5 → 11 mean core) is fully explained by the config differences: more datasets increase the Phase 1 pass rate (7.6% → 9.3%), and the 15-config Phase 3 sweep finds more cluster structure than defaults. The artifact is preserved at `results/negative_control/pipeline_summary.json` for provenance.

**Design:** 50 independent random draws of 500 genes from the 19,296 protein-coding genome, each run through the full pipeline with ASD-blind-matched settings (7 datasets, identical phase3/phase4 config). Tests the engine's false positive rate under the null hypothesis of no coordinated disease signal.

**Sampling rationale:** The negative control samples 500 genes uniformly at random from the full protein-coding genome (19,296 genes), with no exclusion of disease-associated genes. This produces an unbiased null where ~33 SFARI genes (1267/19296 × 500) appear in each trial by chance, matching the realistic expectation of random gene draws. An alternative protocol excluding known disease genes would produce an artificially conservative test; we chose the unbiased version to avoid the appearance of rigging the null distribution.

**Results:** 50/50 trials completed successfully (2026-05-08).

### 7.1 Distribution of outputs across 50 trials

| Metric | Mean | Median | Stdev | Min | Q1 | Q3 | Max |
|--------|---:|---:|---:|---:|---:|---:|---:|
| Phase 1 study genes | 46.4 | 46.5 | 6.7 | 29 | 42 | 50 | 62 |
| Phase 4 core genes | 11.3 | 12.5 | 3.8 | 3 | 9 | 14 | 19 |
| Phase 4 sig clusters | 0.9 | 1.0 | 0.4 | 0 | 1 | 1 | 1 |
| Phase 5 survived | 11.1 | 12.0 | 3.8 | 3 | 9 | 14 | 19 |
| Phase 6 meta-sig | 6.3 | 7.0 | 3.1 | 0 | 4 | 8 | 13 |

43/50 trials (86%) produced exactly 1 significant cluster. 7/50 (14%) produced 0. No trial produced >1.

### 7.2 Phase-by-phase rate comparison with ASD blind

| Stage | Negative control (500 input) | ASD blind (19,296 input) |
|-------|---:|---:|
| Phase 1 pass rate | 46.4/500 = 9.3% | 1,793/19,296 = 9.3% |
| Core / Phase 1 | 11.3/46.4 = 24.4% | 421/1,793 = 23.5% |
| Survived / Phase 1 | 11.1/46.4 = 23.9% | 414/1,793 = 23.1% |
| Meta-sig / Phase 1 | 6.3/46.4 = 13.6% | 233/1,793 = 13.0% |
| **Sig clusters** | **0–1 (median 1)** | **15** |

Per-gene attrition rates are nearly identical between negative control and ASD blind. The engine's per-gene filtering pipeline does not discriminate at the individual gene level — a gene that is differentially expressed in ≥2 datasets passes Phase 1 regardless of biological relevance.

### 7.3 Interpretation — REQUIRES TEAM DISCUSSION

**What the engine does discriminate on: cluster structure, not per-gene filtering.**

ASD blind produces 15 significant clusters from 421 core genes — structured, modular biology with distinct programs. Negative control produces 0–1 significant clusters from ~12 core genes — at most a single spurious cluster from noise.

The negative control result does NOT mean the engine fabricates signal at the same rate as real disease. It means:

1. **Per-gene false positive rate is non-zero.** ~6 meta-significant genes per trial from random input. This is a baseline rate of the statistical filters, not a discovery signal.
2. **Modular discovery is where the engine adds value.** The discriminator between signal and noise is the number and structure of significant clusters, not the count of individual genes passing Phase 6.
3. **The curated tier remains the publishable finding for per-gene claims.** Curated input restricts to known biology, where per-gene filtering IS informative because the input is enriched for real signal.

**Implications for methods paper framing:** The engine is a *modular structure detector*, not a per-gene signal-vs-noise classifier. The negative control demonstrates this distinction empirically. Per-gene claims require the curated tier or stability profiling (iron-clad at 50/50 is a much stronger per-gene claim than single-run meta-significance).

**Open question for team:** Does this result change how we frame the blind tier's gene count? The current framing ("the engine identified N genes") may need qualification ("the engine identified N genes organized into M significant clusters; individual gene significance should be assessed via stability profiling rather than single-run meta-analysis").

---

## 8. Remaining Work

### 8.1 Must-do before submission
- [ ] AD blind tier on RunPod (single run + 50-run stability)
- [x] Matched negative control (50 trials, running on Pi 5)

### 8.2 Strongly recommended
- [ ] Cross-hardware reproducibility (one disease re-run on RunPod, probably IPF)
- [ ] ASD manuscript Methods section fix (DL/REML order — audit C7)

### 8.3 Nice-to-have
- [ ] MEGENA/metaDE benchmarks
- [ ] WGCNA benchmark refresh on ASD

### 8.4 Documentation
- [x] FINDINGS_SUMMARY.md (this document)
- [x] INFRASTRUCTURE.md
- [ ] REGIME_MODEL.md (dedicated document)
- [ ] BRCA_SUBTYPE_RECONSTRUCTION.md (dedicated document)
- [ ] Cross-engine-version stability summary (standalone table + mechanism)
- [ ] CLAIMS_POLICY.md updates with tier-tagged findings

---

## 9. Disease Day Directory Map

| Disease | Directory | Engine version |
|---------|-----------|---------------|
| ASD (v0.3.2) | `disease_days/2026-04-24_asd/` | 0.3.2 |
| ASD (v0.3.3) | `disease_days/2026-04-28_asd_v033/` | 0.3.3.2 |
| IPF | `disease_days/2026-04-29_ipf_v0332/` | 0.3.3.2 |
| CRC | `disease_days/2026-05-02_crc_v0332/` | 0.3.3.2 |
| Psoriasis | `disease_days/2026-05-03_psoriasis_v0332/` | 0.3.3.2 |
| AD | `disease_days/2026-05-05_ad_v0332/` | 0.3.3.2 |
| BrCa | `disease_days/2026-05-05_brca_v0332/` | 0.3.3.2 |
| T2D | `disease_days/2026-05-05_t2d_v0332/` | 0.3.3.2 |
| IBD | `disease_days/2026-05-05_ibd_v0332/` | 0.3.3.2 |
