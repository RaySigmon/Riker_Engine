# Riker Engine — Validation Results

**Engine version:** v0.3.2
**Platform:** Raspberry Pi 5 (8GB) + RunPod (64GB for blind runs)

Eight diseases validated with zero code modifications between runs.

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
├── psoriasis/          # Psoriasis: 50 core genes (independent validation)
│   └── independent_validation/
├── crc/                # CRC: 264 core genes (independent validation)
│   └── independent_validation/
├── negative_control/   # 500 random genes → 5 core (1% FP rate)
├── drug_target_analysis/  # Systematic hit rates across all diseases
├── REPLICATION_LOG.md  # Third-party cold-start replication (6 diseases)
└── INDEPENDENT_VALIDATION.md  # Psoriasis + CRC by independent AI agents
```

## Summary

| Disease | Tissue | Seeds | Datasets | Core genes | Survived | Meta-sig |
|---------|--------|-------|----------|-----------|----------|----------|
| ASD | Brain cortex | 1,267 | 7 | 35 | 20 (57.1%) | 9 |
| T2D | Pancreatic islets | 443 | 4 | 8 | 8 (100%) | 8 |
| IBD | Intestinal mucosa | 762 | 6 | 304 | 302 (99.3%) | 296 |
| AD | Brain cortex | 801 | 5 | 394 | 340 (86.3%) | 312 |
| Breast Ca. | Breast tumor | 653 | 5 | 152 | 139 (91.4%) | 112 |
| **IPF** | **Lung** | **354** | **5** | **190** | **170 (89.5%)** | **157** |
| Psoriasis* | Skin | 96 | 5 | 50 | 50 (100%) | 28 |
| CRC* | Colon | 515 | 6 | 264 | 245 (92.8%) | 219 |

\*Validated by independent AI agents with no author involvement.

## IPF Cold Replication

The IPF validation includes independent verification in a held-out dataset
(GSE47460) the engine never saw:
- 132/153 genes concordant + significant (**86.3%**)
- 96.7% directional concordance
- 52 iron-clad genes survive all leave-one-out configurations + cold replication
- FAM107A identified as novel candidate (zero IPF literature)

See `results/ipf/cold_replication/` and `results/ipf/loo_stability/`.

## ASD Stability Profiling

50-run stability profiling with varied UMAP seeds (see `stability_ASD_blind/`):
- Core gene range: 389–410 (mean 401)
- 376 iron-clad genes (>=90% of runs), 33 borderline, 29 stochastic
- 289 novel non-SFARI genes at 50/50 stability
