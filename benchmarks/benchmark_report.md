# WGCNA vs Riker Engine — ASD Brain Cortex Benchmark

Generated: 2026-03-22 17:05

## Runtime & Resource Comparison

| | Riker Engine | WGCNA |
|---|---|---|
| **Total wall time** | ~8 minutes (all 6 phases, 3 discovery + 4 replication datasets, 10K permutations) | ~3 hours 50 min (3 discovery datasets only, no replication) |
| **Memory** | Ran full gene set (~18K genes) without issue | OOM-killed on GSE28475 (18K genes); required filtering to top 10K most variable genes |
| **Preprocessing** | None required — handles raw series matrices directly | Required gene variance filtering to fit in 8GB RAM (standard WGCNA practice, but a forced compromise) |
| **Cross-dataset validation** | Built-in: Phases 4-5 pre-specify and replicate across held-out datasets | None: each dataset analyzed independently; cross-dataset comparison must be done manually |
| **Hardware** | Raspberry Pi 5, 8GB RAM | Same — but exceeded memory on one dataset |

Note: GSE28475 failed to reach R² > 0.80 for scale-free topology fit, even after gene filtering. WGCNA used power=20 (best available R²=0.35), which means the network for that dataset does not meet WGCNA's own quality criteria. Additionally, WGCNA found only 1 module (+ grey) for GSE28475 — effectively no useful clustering.

## Riker Engine Results (reference)

- **Core genes**: 35
- **Clusters**: 8
  - Cluster 1 (3 genes): EFR3A, FBXO33, SLC45A1
  - Cluster 2 (10 genes): ATP2B2, DEPDC5, GABRG2, NR3C2, PPP3CA, PRPF19, QRICH1, SLC12A5, SPTAN1, ZNF385B
  - Cluster 3 (4 genes): NAV2, NPAS3, TBR1, ZFHX3
  - Cluster 4 (5 genes): ACTL6B, APBA2, ICA1, SLC25A12, SLC25A27
  - Cluster 5 (3 genes): ATP1A1, DPP10, ITPR1
  - Cluster 7 (4 genes): LIN7B, NPTN, RGS7, RPH3A
  - Cluster 9 (3 genes): ELP2, KIAA0232, SNX14
  - Cluster 12 (3 genes): CADPS, SEZ6L2, SNAP25

## GSE28521 — WGCNA Results

- **Total modules found**: 12
- **Significant modules (p<0.05)**: 4
  - Module 5: r=-0.294, p=0.0085, 653 genes
  - Module 6: r=0.309, p=0.0056, 648 genes
  - Module 11: r=0.299, p=0.0074, 80 genes
  - Module 12: r=0.278, p=0.0131, 46 genes
- **Core genes found in dataset**: 35/35
- **Core genes in grey (unassigned)**: 0
- **Core genes assigned to modules**: 35
- **Core genes in significant modules**: 17/35

### Core gene module distribution in GSE28521

| WGCNA Module | # Core Genes | Core Genes | Significant? |
|---|---|---|---|
| 1 (1848 total) | 14 | CADPS, DEPDC5, DPP10, ELP2, FBXO33, KIAA0232, LIN7B, NPTN, PPP3CA, PRPF19, RGS7, SEZ6L2, SNX14, SPTAN1 | No |
| 2 (1516 total) | 1 | NAV2 | No |
| 3 (899 total) | 2 | NR3C2, SLC45A1 | No |
| 4 (783 total) | 1 | RPH3A | No |
| 5 (653 total) | 15 | ACTL6B, APBA2, ATP1A1, ATP2B2, EFR3A, GABRG2, ICA1, ITPR1, QRICH1, SLC12A5, SLC25A12, SLC25A27, SNAP25, TBR1, ZNF385B | Yes |
| 6 (648 total) | 2 | NPAS3, ZFHX3 | Yes |

## GSE28475 — WGCNA Results

- **Total modules found**: 1
- **Significant modules (p<0.05)**: 1
  - Module 0: r=0.230, p=0.0144, 5 genes
- **Core genes found in dataset**: 31/35
- **Core genes in grey (unassigned)**: 0
- **Core genes assigned to modules**: 31
- **Core genes in significant modules**: 0/35
- **Core genes absent from dataset**: DEPDC5, PRPF19, SLC45A1, ZNF385B

### Core gene module distribution in GSE28475

| WGCNA Module | # Core Genes | Core Genes | Significant? |
|---|---|---|---|
| 1 (9995 total) | 31 | ACTL6B, APBA2, ATP1A1, ATP2B2, CADPS, DPP10, EFR3A, ELP2, FBXO33, GABRG2, ICA1, ITPR1, KIAA0232, LIN7B, NAV2, NPAS3, NPTN, NR3C2, PPP3CA, QRICH1, RGS7, RPH3A, SEZ6L2, SLC12A5, SLC25A12, SLC25A27, SNAP25, SNX14, SPTAN1, TBR1, ZFHX3 | No |

## GSE64018 — WGCNA Results

- **Total modules found**: 23
- **Significant modules (p<0.05)**: 12
  - Module 1: r=-0.597, p=0.0021, 3075 genes
  - Module 5: r=0.672, p=0.0003, 943 genes
  - Module 6: r=-0.504, p=0.0120, 737 genes
  - Module 7: r=-0.546, p=0.0057, 663 genes
  - Module 9: r=-0.514, p=0.0102, 483 genes
  - Module 10: r=0.472, p=0.0198, 392 genes
  - Module 14: r=0.575, p=0.0033, 203 genes
  - Module 15: r=0.455, p=0.0254, 178 genes
  - Module 16: r=0.551, p=0.0052, 156 genes
  - Module 17: r=-0.588, p=0.0025, 152 genes
  - Module 18: r=-0.670, p=0.0003, 92 genes
  - Module 23: r=-0.487, p=0.0157, 34 genes
- **Core genes found in dataset**: 35/35
- **Core genes in grey (unassigned)**: 0
- **Core genes assigned to modules**: 35
- **Core genes in significant modules**: 33/35

### Core gene module distribution in GSE64018

| WGCNA Module | # Core Genes | Core Genes | Significant? |
|---|---|---|---|
| 1 (3075 total) | 14 | ACTL6B, ATP1A1, DPP10, EFR3A, ELP2, FBXO33, GABRG2, NPTN, PPP3CA, QRICH1, RGS7, SLC25A12, SLC25A27, SNX14 | Yes |
| 2 (2836 total) | 2 | NAV2, ZFHX3 | No |
| 6 (737 total) | 5 | ICA1, ITPR1, NR3C2, SPTAN1, ZNF385B | Yes |
| 7 (663 total) | 9 | APBA2, CADPS, LIN7B, PRPF19, RPH3A, SEZ6L2, SLC12A5, SNAP25, TBR1 | Yes |
| 9 (483 total) | 1 | KIAA0232 | Yes |
| 16 (156 total) | 1 | NPAS3 | Yes |
| 17 (152 total) | 2 | DEPDC5, SLC45A1 | Yes |
| 18 (92 total) | 1 | ATP2B2 | Yes |

## Cross-Dataset Consistency

### Core gene placement across datasets

| Gene | Riker Cluster | GSE28521 | GSE28475 | GSE64018 | Consistent? |
|---|---|---|---|---|---|
| ACTL6B | 4 | 5* | 1 | 1* | Yes |
| APBA2 | 4 | 5* | 1 | 7* | Yes |
| ATP1A1 | 5 | 5* | 1 | 1* | Yes |
| ATP2B2 | 2 | 5* | 1 | 18* | Yes |
| CADPS | 12 | 1 | 1 | 7* | Yes |
| DEPDC5 | 2 | 1 | — | 17* | No |
| DPP10 | 5 | 1 | 1 | 1* | Yes |
| EFR3A | 1 | 5* | 1 | 1* | Yes |
| ELP2 | 9 | 1 | 1 | 1* | Yes |
| FBXO33 | 1 | 1 | 1 | 1* | Yes |
| GABRG2 | 2 | 5* | 1 | 1* | Yes |
| ICA1 | 4 | 5* | 1 | 6* | Yes |
| ITPR1 | 5 | 5* | 1 | 6* | Yes |
| KIAA0232 | 9 | 1 | 1 | 9* | Yes |
| LIN7B | 7 | 1 | 1 | 7* | Yes |
| NAV2 | 3 | 2 | 1 | 2 | Yes |
| NPAS3 | 3 | 6* | 1 | 16* | Yes |
| NPTN | 7 | 1 | 1 | 1* | Yes |
| NR3C2 | 2 | 3 | 1 | 6* | Yes |
| PPP3CA | 2 | 1 | 1 | 1* | Yes |
| PRPF19 | 2 | 1 | — | 7* | No |
| QRICH1 | 2 | 5* | 1 | 1* | Yes |
| RGS7 | 7 | 1 | 1 | 1* | Yes |
| RPH3A | 7 | 4 | 1 | 7* | Yes |
| SEZ6L2 | 12 | 1 | 1 | 7* | Yes |
| SLC12A5 | 2 | 5* | 1 | 7* | Yes |
| SLC25A12 | 4 | 5* | 1 | 1* | Yes |
| SLC25A27 | 4 | 5* | 1 | 1* | Yes |
| SLC45A1 | 1 | 3 | — | 17* | No |
| SNAP25 | 12 | 5* | 1 | 7* | Yes |
| SNX14 | 9 | 1 | 1 | 1* | Yes |
| SPTAN1 | 2 | 1 | 1 | 6* | Yes |
| TBR1 | 3 | 5* | 1 | 7* | Yes |
| ZFHX3 | 3 | 6* | 1 | 2 | Yes |
| ZNF385B | 2 | 5* | — | 6* | No |

*Module numbers marked with * are significant (p<0.05)*

### Riker cluster cohesion in WGCNA

Do genes that Riker groups into the same cluster also land in the same WGCNA module?

**GSE28521**:
  - Riker cluster 1 (3 genes): **SCATTERED** — mod 1: FBXO33; mod 3: SLC45A1; mod 5: EFR3A
  - Riker cluster 2 (10 genes): **SCATTERED** — mod 1: DEPDC5, PPP3CA, PRPF19, SPTAN1; mod 3: NR3C2; mod 5: ATP2B2, GABRG2, QRICH1, SLC12A5, ZNF385B
  - Riker cluster 3 (4 genes): **SCATTERED** — mod 2: NAV2; mod 5: TBR1; mod 6: NPAS3, ZFHX3
  - Riker cluster 4 (5 genes): **COHESIVE** — mod 5: ACTL6B, APBA2, ICA1, SLC25A12, SLC25A27
  - Riker cluster 5 (3 genes): **SCATTERED** — mod 1: DPP10; mod 5: ATP1A1, ITPR1
  - Riker cluster 7 (4 genes): **SCATTERED** — mod 1: LIN7B, NPTN, RGS7; mod 4: RPH3A
  - Riker cluster 9 (3 genes): **COHESIVE** — mod 1: ELP2, KIAA0232, SNX14
  - Riker cluster 12 (3 genes): **SCATTERED** — mod 1: CADPS, SEZ6L2; mod 5: SNAP25

**GSE28475**:
  - Riker cluster 1 (3 genes): **SCATTERED** — mod 1: EFR3A, FBXO33; absent: SLC45A1
  - Riker cluster 2 (10 genes): **SCATTERED** — mod 1: ATP2B2, GABRG2, NR3C2, PPP3CA, QRICH1, SLC12A5, SPTAN1; absent: DEPDC5, PRPF19, ZNF385B
  - Riker cluster 3 (4 genes): **COHESIVE** — mod 1: NAV2, NPAS3, TBR1, ZFHX3
  - Riker cluster 4 (5 genes): **COHESIVE** — mod 1: ACTL6B, APBA2, ICA1, SLC25A12, SLC25A27
  - Riker cluster 5 (3 genes): **COHESIVE** — mod 1: ATP1A1, DPP10, ITPR1
  - Riker cluster 7 (4 genes): **COHESIVE** — mod 1: LIN7B, NPTN, RGS7, RPH3A
  - Riker cluster 9 (3 genes): **COHESIVE** — mod 1: ELP2, KIAA0232, SNX14
  - Riker cluster 12 (3 genes): **COHESIVE** — mod 1: CADPS, SEZ6L2, SNAP25

**GSE64018**:
  - Riker cluster 1 (3 genes): **SCATTERED** — mod 1: EFR3A, FBXO33; mod 17: SLC45A1
  - Riker cluster 2 (10 genes): **SCATTERED** — mod 1: GABRG2, PPP3CA, QRICH1; mod 17: DEPDC5; mod 18: ATP2B2; mod 6: NR3C2, SPTAN1, ZNF385B; mod 7: PRPF19, SLC12A5
  - Riker cluster 3 (4 genes): **SCATTERED** — mod 16: NPAS3; mod 2: NAV2, ZFHX3; mod 7: TBR1
  - Riker cluster 4 (5 genes): **SCATTERED** — mod 1: ACTL6B, SLC25A12, SLC25A27; mod 6: ICA1; mod 7: APBA2
  - Riker cluster 5 (3 genes): **SCATTERED** — mod 1: ATP1A1, DPP10; mod 6: ITPR1
  - Riker cluster 7 (4 genes): **SCATTERED** — mod 1: NPTN, RGS7; mod 7: LIN7B, RPH3A
  - Riker cluster 9 (3 genes): **SCATTERED** — mod 1: ELP2, SNX14; mod 9: KIAA0232
  - Riker cluster 12 (3 genes): **COHESIVE** — mod 7: CADPS, SEZ6L2, SNAP25

## Summary

### Recovery of Riker's 35 core genes by WGCNA

- Core genes in a significant WGCNA module in **at least one** dataset: **34/35** (97.1%)
- Core genes **never** in any significant WGCNA module: **1/35** (2.9%)
  - Genes WGCNA misses: NAV2
- Core genes assigned to a non-grey module in **all 3 datasets**: **31/35** (88.6%)

### Scale comparison
- GSE28521: WGCNA significant module genes = **1427** (vs Riker's 35 core genes)
- GSE28475: WGCNA significant module genes = **5** (vs Riker's 35 core genes)
- GSE64018: WGCNA significant module genes = **7108** (vs Riker's 35 core genes)

### Key takeaways

1. **Specificity**: Riker Engine identifies 35 core genes through progressive filtering. WGCNA's significant modules contain orders of magnitude more genes, making biological interpretation harder.
2. **Cross-dataset consistency**: Riker requires genes to survive across multiple datasets by design. WGCNA runs independently per dataset with no built-in cross-dataset validation.
3. **Cluster cohesion**: See per-dataset analysis above for whether Riker's gene groupings are preserved or scattered by WGCNA.

---

### Methodological note

This comparison demonstrates the **workflow advantage** of integrated multi-dataset
analysis over manual single-dataset analysis followed by cross-study comparison. It
is not a claim that the Riker Engine's clustering method is superior to WGCNA's
network construction — the tools answer different questions at different scales.

WGCNA is designed for comprehensive single-dataset network analysis and excels at
identifying the full co-expression structure within a cohort. The Riker Engine is
designed to extract a minimal, replicated gene set across multiple independent
datasets. Comparing output sizes (35 vs. 1,427–9,995 genes) reflects this
difference in design goals, not a quality difference in the underlying methods.

A fair assessment: WGCNA correctly identifies the disease biology (34/35 Riker core
genes appear in WGCNA significant modules) but embeds it within much larger modules
that require manual curation to extract actionable targets. The Riker Engine
automates this extraction by design.

---
*Benchmark run on Raspberry Pi 5 (Ghost)*