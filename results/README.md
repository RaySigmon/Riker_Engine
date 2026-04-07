# Riker Engine — Validation Results

**Engine version:** v0.3.1
**Platform:** Raspberry Pi 5 (8GB) + RunPod (64GB for blind runs)

Six diseases validated with zero code modifications between runs.

## Results Directory

```
results/
├── asd_test_run_1/     # ASD: 35 core genes (curated + allgenes blind)
├── t2d_test_run_1/     # T2D: 8 core genes (curated + allgenes blind)
├── ibd_test_run_1/     # IBD: 304 core genes (curated + allgenes blind)
├── ad/                 # AD: 394 curated, 11,481 blind
│   ├── curated/
│   └── blind/
├── breast_cancer/      # Breast cancer: 152 curated, 3,675 blind
│   ├── curated/
│   └── blind/
├── ipf/                # IPF: 190 curated + cold replication + LOO stability
│   ├── curated/
│   ├── cold_replication/
│   └── loo_stability/
├── negative_control/   # 500 random genes → 5 core (1% FP rate)
└── drug_target_analysis/  # Systematic hit rates across all diseases
```

## Summary

| Disease | Tissue | Seeds | Datasets | Core genes | Survived | Meta-sig |
|---------|--------|-------|----------|-----------|----------|----------|
| ASD | Brain cortex | 1,267 | 7 | 35 | 35 (100%) | 13 |
| T2D | Pancreatic islets | 443 | 4 | 8 | 8 (100%) | 8 |
| IBD | Intestinal mucosa | 762 | 6 | 304 | 302 (99.3%) | 296 |
| AD | Brain cortex | 801 | 5 | 394 | 340 (86.3%) | 312 |
| Breast Ca. | Breast tumor | 653 | 5 | 152 | 152 (100%) | 121 |
| **IPF** | **Lung** | **354** | **5** | **190** | **170 (89.5%)** | **157** |

## IPF Cold Replication

The IPF validation includes independent verification in a held-out dataset
(GSE47460) the engine never saw:
- 132/153 genes concordant + significant (**86.3%**)
- 96.7% directional concordance
- 52 iron-clad genes survive all leave-one-out configurations + cold replication
- FAM107A identified as novel candidate (zero IPF literature)

See `results/ipf/cold_replication/` and `results/ipf/loo_stability/`.
