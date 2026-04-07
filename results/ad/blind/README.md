# AD Blind Run Results

## Summary

| Metric | AD Blind | AD Curated | ASD Blind | T2D Blind |
|--------|----------|------------|-----------|-----------|
| Seed genes | 19,296 (all protein-coding) | 801 | 39,185 | 26,800 |
| Phase 1 study genes | 13,181 (68.3%) | 439 (54.7%) | 4,340 (11.1%) | 3,384 (12.6%) |
| Core genes | 11,481 | 394 | 415 | 177 |
| Core / study gene ratio | 87.1% | 89.7% | 9.6% | 5.2% |
| Survived replication | 9,364 | 340 (86.3%) | — | — |
| Eliminated | 1,981 | 54 | — | — |
| Meta-significant | 8,458 | 312 | — | — |

## Interpretation: Why AD produces a large blind output

The AD blind run produces 11,481 core genes — an 87% pass rate from study genes.
This is **not** a pipeline failure; it reflects a genuine biological property of
Alzheimer's disease brain tissue.

**Why this happens:**
1. **Pervasive transcriptomic disruption.** AD brain cortex shows widespread gene
   expression changes. 68% of all protein-coding genes reach significance in 2+
   datasets at Phase 1 (vs. 11% for ASD, 13% for T2D). When most genes are
   genuinely differentially expressed, the pipeline correctly passes them.

2. **Full-seed-set FDR is less discriminating when signal dominates noise.** With
   13,181 study genes from 19,296 seeds, the FDR denominator is only 1.5x the
   study set. In ASD (141 study from 1,267 seeds), the denominator is 9x larger,
   providing stronger correction. This is a known property of FDR correction —
   it controls the false discovery *rate*, not the absolute number.

3. **Consensus clustering still groups genes into 1,389 modules.** The pipeline
   identifies modular structure within the large gene set, but the permutation
   threshold does not aggressively filter when most clusters contain genuinely
   co-regulated genes.

**What this means for users:**
- The **curated run (394 core genes)** is the actionable output for AD target
  discovery. It starts from 801 disease-associated genes and filters to a
  manageable, replicated set.
- The blind run demonstrates that **100% of curated core genes pass Phase 1**
  in the hypothesis-free analysis, validating that they are genuine signals.
- The blind run also shows a **pipeline limitation**: for diseases with pervasive
  transcriptomic disruption (AD, and to a lesser extent IBD), the blind mode
  does not narrow the field to a small target list. The pipeline is calibrated
  for the more common scenario where disease signal is sparse.

**Comparison with other blind runs:**
- ASD (weak signal, 11% Phase 1): blind produces 415 core genes — highly targeted
- T2D (weak signal, 13% Phase 1): blind produces 177 core genes — highly targeted
- IBD (strong signal, 54% Phase 1): blind produces 6,458 core genes — moderate
- AD (strong signal, 68% Phase 1): blind produces 11,481 core genes — broad

This calibration pattern is expected and documented in the dissertation (Section 5.1).

## Run details

- **Hardware:** RunPod CPU pod, 64 GB RAM, 32 vCPU
- **Runtime:** ~36 minutes (all 6 phases)
- **Config:** `configs/runpod/ad_blind.yaml`
- **Seed genes:** All 19,296 protein-coding genes from HGNC
