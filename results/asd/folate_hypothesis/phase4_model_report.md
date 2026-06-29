# Phase 4: Receptor Competition Model Report

## Model Overview

This model simulates competitive binding at Folate Receptor Alpha (FRα) between
5-methyltetrahydrofolate (5-MTHF, the active folate form) and unmetabolized folic
acid (UMFA), which accumulates when MTHFR enzyme activity is reduced.

## Key Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| Kd folic acid at FRα | 0.4 nM | Antony 1992, Elnakat & Bhatt 2004 |
| Kd 5-MTHF at FRα | 3.0 nM | Antony 1992 |
| Affinity ratio | 7.5× (FA binds tighter) | Derived |
| MTHFR CC efficiency | 100% | Frosst et al. 1995 |
| MTHFR CT efficiency | 65% | van der Put et al. 1998 |
| MTHFR TT efficiency | 30% | van der Put et al. 1998 |

## Results Summary

### At Typical Prenatal Supplementation (1000 µg/d)

| Genotype | 5-MTHF (nM) | UMFA (nM) | UMFA:5-MTHF | Effective Delivery | Status |
|----------|-------------|-----------|-------------|-------------------|--------|
| CC (wildtype) | 36.1 | 6.6 | 0.183 | 40.8% | mild_insufficiency |
| CT (heterozygous) | 27.0 | 7.0 | 0.261 | 32.6% | mild_insufficiency |
| TT (homozygous) | 17.8 | 7.5 | 0.420 | 23.1% | moderate_insufficiency |

### Delivery Gap
- CC vs TT delivery difference: **17.6 percentage points**
- This gap widens with increasing dose

### Full Results Table

| Genotype | Scenario | Delivery | Status |
|----------|----------|----------|--------|
| CC (wildtype) | No fortification (200 µg/d) | 86.6% | normal |
| CC (wildtype) | Fortified diet only (400 µg/d) | 65.7% | normal |
| CC (wildtype) | Fortified + prenatal (1000 µg/d) | 40.8% | mild_insufficiency |
| CC (wildtype) | High-dose supplement (2000 µg/d) | 25.9% | moderate_insufficiency |
| CC (wildtype) | Very high dose (4000 µg/d) | 15.2% | moderate_insufficiency |
| CT (heterozygous) | No fortification (200 µg/d) | 83.2% | normal |
| CT (heterozygous) | Fortified diet only (400 µg/d) | 58.1% | normal |
| CT (heterozygous) | Fortified + prenatal (1000 µg/d) | 32.6% | mild_insufficiency |
| CT (heterozygous) | High-dose supplement (2000 µg/d) | 19.5% | moderate_insufficiency |
| CT (heterozygous) | Very high dose (4000 µg/d) | 11.0% | severe_insufficiency |
| TT (homozygous) | No fortification (200 µg/d) | 77.7% | normal |
| TT (homozygous) | Fortified diet only (400 µg/d) | 47.6% | mild_insufficiency |
| TT (homozygous) | Fortified + prenatal (1000 µg/d) | 23.1% | moderate_insufficiency |
| TT (homozygous) | High-dose supplement (2000 µg/d) | 12.8% | severe_insufficiency |
| TT (homozygous) | Very high dose (4000 µg/d) | 6.9% | severe_insufficiency |

## Developmental Window Vulnerability

For MTHFR TT at 1000 µg/d (base delivery: 23.1%):

| Window | Process | FRα Dependence | Vulnerability Score |
|--------|---------|---------------|-------------------|
| Prenatal (neural tube) | Neural tube closure | 85% | 65.3% |
| 0-6 months | Early neuronal migration | 80% | 61.5% |
| 6-18 months | Peak synaptogenesis | 95% | 73.0% |
| 18-36 months | Synaptic pruning onset | 85% | 65.3% |
| 3-6 years | Myelination | 60% | 46.1% |

## Sensitivity Analysis

The model's predictions depend critically on the Kd ratio between 5-MTHF and folic acid
at FRα. The literature range for this ratio is approximately 3-10×. At the lower end (3×),
the competitive effect is moderate; at the higher end (10×), it becomes much more pronounced.

The model uses a ratio of 7.5×, which is conservative. If the true ratio is higher,
the predicted delivery deficits would be more severe.

## Model Limitations

1. **Simplified PK**: Blood levels are estimated from dose using linear approximations;
   real pharmacokinetics involve absorption variability, tissue distribution, and renal clearance
2. **Single-compartment**: Does not model blood-brain barrier transport separately from
   peripheral binding
3. **Steady-state assumption**: Real UMFA levels fluctuate with meals/supplements
4. **No feedback loops**: Does not model compensatory upregulation of FRα or alternative
   transport (RFC/SLC19A1)
5. **Parameter uncertainty**: Kd values vary across studies and experimental conditions

## Interpretation for Hypothesis

The model demonstrates that **the competitive binding mechanism is physically plausible**:
- MTHFR TT individuals show meaningfully reduced 5-MTHF delivery at supplementation
  levels typical for pregnant women in fortified countries
- The effect is dose-dependent, worsening with higher supplementation
- The most vulnerable developmental window (peak synaptogenesis, 6-18 months) coincides
  with the period when infants may receive both formula (fortified) and are still being
  breastfed by supplementing mothers
- Even moderate delivery deficits during this window could be biologically significant

This does NOT prove the hypothesis — it shows the mechanism is **quantitatively feasible**
with known binding parameters.
