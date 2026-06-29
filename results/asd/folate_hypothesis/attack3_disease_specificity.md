# Attack 3: Disease Specificity Test

**Question:** Is the folate/immune/CDR signal found in ASD actually specific to ASD, or does it appear in every disease the Riker Engine analyzes?

**Date:** 2026-04-17

## Critical Caveat

ASD data comes from 50-run stability profiling (376 iron-clad genes). All other diseases are single-run results (lower confidence). A gene absent from another disease may simply not have been detected in that single run. This asymmetry means we can confirm cross-disease presence but cannot definitively confirm absence.

---

## Diseases Compared

| Disease | Core Genes | Data Quality |
|---------|-----------|--------------|
| ASD | 438 (376 iron-clad) | 50-run stability (high) |
| Alzheimer's (AD) | 394 | Single run |
| Breast Cancer | 152 | Single run |
| Colorectal Cancer (CRC) | 264 | Single run |
| Psoriasis | 50 | Single run |
| T2D | 8 | Single run (too small for meaningful comparison) |

---

## Results by Gene Set

### Set 1: One-Carbon Metabolism (26 genes tested)

| Gene | ASD | AD | Breast | CRC | Psoriasis | T2D |
|------|-----|-----|--------|-----|-----------|-----|
| MAT2B | iron-clad | - | - | - | - | - |
| MTHFD1 | - | down (c10) | - | - | - | - |
| DNMT1 | - | - | - | up (c12) | - | - |
| DNMT3A | - | - | - | up (c17) | - | - |
| DNMT3B | - | - | - | up (c42) | - | - |

All other one-carbon genes (MTHFR, MTHFD2, MTR, MTRR, MAT1A, MAT2A, SHMT1/2, DHFR, TYMS, ALDH1L1/2, SLC19A1, SLC46A1, FOLR1/2, AHCY, BHMT, CBS, CTH): **not found in any disease**.

**Summary:** 1/26 in ASD (MAT2B), 1/26 in AD (MTHFD1), 3/26 in CRC (all DNMTs), 0 in Breast/Psoriasis/T2D.

**Interpretation:** The one-carbon signal is sparse everywhere. MAT2B is unique to ASD. CRC shows DNMT upregulation, which is a known cancer epigenetics signal, not folate metabolism per se. AD's single MTHFD1 hit is a different gene than ASD's MAT2B. No disease shares ASD's specific one-carbon finding.

### Set 2: Immune Signature (15 genes tested)

| Gene | ASD | AD | Breast | CRC | Psoriasis | T2D |
|------|-----|-----|--------|-----|-----------|-----|
| C1QB | iron-clad (1.00) | - | - | - | - | - |
| C1QC | iron-clad (1.00) | - | - | - | - | - |
| B2M | iron-clad (1.00) | - | - | **down** (c22) | - | - |
| CD14 | iron-clad (1.00) | - | - | - | - | - |
| IRF8 | iron-clad (1.00) | - | - | - | - | - |
| IFITM3 | iron-clad (1.00) | - | - | - | - | - |
| TRAF3IP2 | iron-clad (1.00) | - | - | - | **up** (c6) | - |
| TXNIP | iron-clad (1.00) | - | **down** (c16) | - | - | - |
| S100A8 | iron-clad (1.00) | - | - | - | - | - |
| NFE2L2 | iron-clad (1.00) | - | - | - | - | - |
| LYN | iron-clad (1.00) | - | - | - | - | - |
| NFKBIZ | iron-clad (0.98) | - | - | - | - | - |
| CXCL16 | iron-clad (1.00) | - | - | - | - | - |
| STAT4 | iron-clad (1.00) | - | - | - | - | - |
| THEMIS | iron-clad (1.00) | - | - | - | - | - |

**Summary:** ASD: 15/15. AD: 0/15. Breast: 1/15 (TXNIP, down). CRC: 1/15 (B2M, down). Psoriasis: 1/15 (TRAF3IP2, up). T2D: 0/15.

**Interpretation:** The 15-gene immune signature is overwhelmingly ASD-specific. Even the three individual hits in other diseases are singletons (not the coordinated multi-gene pattern seen in ASD), and notably:
- CRC's B2M is **downregulated** (immune evasion, a known cancer mechanism), not the immune activation pattern seen in ASD
- Breast cancer's TXNIP is **downregulated** (metabolic, not immune context)
- Psoriasis's TRAF3IP2 is expected (IL-17 signaling is the canonical psoriasis pathway) but appears alone, without the complement/innate immune genes that define the ASD pattern

The complement pair C1QB+C1QC, the microglial marker CD14, the interferon response (IFITM3, IRF8), and the adaptive immune components (STAT4, THEMIS) are entirely absent from all other diseases.

### Set 3: CDR (Cell Danger Response) Markers (17 genes tested)

| Gene | ASD | AD | Breast | CRC | Psoriasis | T2D |
|------|-----|-----|--------|-----|-----------|-----|
| ENTPD6 | iron-clad (1.00) | - | - | - | - | - |
| MAT2B | iron-clad (1.00) | - | - | - | - | - |
| KYAT3 | iron-clad (1.00) | - | - | - | - | - |
| ASAH1 | - | down (c38) | - | - | - | - |
| PTGS2 | - | - | - | up (c33) | - | - |

All other CDR genes (HMGB1, PANX1/2, SMPD1/3, CERS2/6, SPHK1, DNM1L, MFN1/2, ALOX5): **not found in any disease**.

**Summary:** ASD: 3/17. AD: 1/17 (ASAH1). CRC: 1/17 (PTGS2). Others: 0/17.

**Interpretation:** CDR markers are sparse in all diseases. ASD has the most (3), including the purinergic ENTPD6 and kynurenine pathway KYAT3, which are absent everywhere else. CRC's PTGS2 (COX-2) is a well-known CRC pathway gene, not a CDR signal. AD's ASAH1 is a ceramide enzyme, appearing alone without the broader sphingolipid pattern.

### Set 4: ASD-Specific Genes (MAT2B + Cluster 21)

| Gene | ASD | AD | Breast | CRC | Psoriasis | T2D |
|------|-----|-----|--------|-----|-----------|-----|
| MAT2B | iron-clad (1.00) | - | - | - | - | - |
| NFS1 | iron-clad (1.00) | - | - | - | - | - |
| UQCRC1 | iron-clad (1.00) | - | - | - | - | - |
| GLS2 | iron-clad (1.00) | - | - | - | - | - |

**Summary:** All 4 genes appear ONLY in ASD. Zero of these genes appear in any other disease.

**Interpretation:** The MAT2B finding and the mitochondrial cluster (NFS1/UQCRC1/GLS2) are completely ASD-specific across all diseases tested.

---

## Overall Gene Overlap

ASD iron-clad genes (376) overlapping with other disease core gene sets:

| Comparison | Overlap | % of ASD iron-clad |
|------------|---------|-------------------|
| ASD vs AD | 13 genes | 3.5% |
| ASD vs Breast Cancer | 7 genes | 1.9% |
| ASD vs CRC | 8 genes | 2.1% |
| ASD vs Psoriasis | 1 gene | 0.3% |
| ASD vs T2D | 0 genes | 0.0% |

The overlapping genes are general cellular machinery (ABCA1, SOX4, etc.), not the folate/immune/CDR signature genes.

---

## Answers to Key Questions

### Q1: Is MAT2B in any other disease core set?
**NO.** MAT2B appears only in ASD, where it is iron-clad (50/50 runs, 100% stability). It is absent from AD, breast cancer, CRC, psoriasis, and T2D.

### Q2: Do Cluster 21 genes (NFS1, UQCRC1, GLS2) co-occur in other diseases?
**NO.** None of the three genes appear in any other disease. The mitochondrial iron-sulfur/ETC/glutamine cluster is entirely ASD-specific.

### Q3: How many immune genes appear in AD vs non-neurological diseases?
Strikingly, **AD shows zero** of the 15 ASD immune genes despite being neurological and having a comparable core gene count (394 vs 438). This is notable because AD has its own well-characterized neuroinflammatory signature (TREM2, complement C4A, etc.) that is compositionally distinct from the ASD immune pattern. The ASD immune signature is not simply "neuroinflammation" -- it is a specific configuration.

### Q4: Does the immune signature composition differ between ASD and other diseases?
**Yes, dramatically.** ASD shows a coordinated 15-gene immune pattern spanning complement (C1QB/C1QC), innate immunity (CD14, S100A8, NFE2L2), interferon response (IFITM3, IRF8), NF-kB signaling (NFKBIZ, TRAF3IP2), and adaptive immune regulation (STAT4, THEMIS). No other disease shows more than 1 of these 15 genes, and even those singletons appear in different functional contexts (cancer immune evasion, psoriasis IL-17 pathway).

---

## Kill Condition Evaluation

**Kill condition:** MAT2B + immune signature + CDR markers appear identically across all diseases.

**Result: SURVIVED.** The kill condition is decisively not met:
- MAT2B is found in zero other diseases
- The 15-gene immune signature is found in zero other diseases as a group (max 1/15 anywhere)
- CDR markers ENTPD6 and KYAT3 are found in zero other diseases
- The NFS1/UQCRC1/GLS2 mitochondrial cluster is found in zero other diseases

---

## Verdict

**The folate/immune/CDR signal is highly specific to ASD.** It is not a generic tissue-stress artifact that appears in every Riker Engine analysis.

**Strength of evidence:** Strong, with caveats.

**Caveats that prevent this from being conclusive:**
1. **Asymmetric data quality:** ASD has 50-run stability; others have single runs. Some genes absent from other diseases might appear with stability profiling. However, the 15/15 vs 0-1/15 immune gap is too large to explain by sensitivity alone.
2. **Limited disease panel:** Five comparator diseases is useful but not exhaustive. Other neurodevelopmental conditions (schizophrenia, ADHD, intellectual disability) would be more informative comparators.
3. **Single-run CRC found 3 DNMTs:** If CRC stability profiling were done, it might reveal more one-carbon genes. This would not undermine ASD specificity (CRC epigenetic dysregulation is a known distinct mechanism) but would reduce the uniqueness of "one-carbon pathway involvement" as a talking point.

**What would strengthen this further:**
- Run stability profiling on AD (the closest neurological comparator)
- Test against schizophrenia, ADHD, or epilepsy datasets
- Confirm the ASD immune gene directions (need phase4 data, not just stability scores)
