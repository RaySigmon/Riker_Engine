# IPF Disease Report — v0.3.3.2 (DRAFT)

---

## 1. Run summary

- **Disease:** IPF (Idiopathic Pulmonary Fibrosis)
- **Date:** 2026-04-29
- **Engine version:** 0.3.3.2
- **Engine commit:** 939663e0f1f39e01efd53f61253b73f084ee74fd
- **Master seed:** 42
- **Tier 1 wall clock:** 1 min 36 sec
- **Tier 2 wall clock:** 13 min 31 sec
- **Tier 3 wall clock:** 33272 sec (9.2 hours)
- Three tiers completed: curated single run, protein-coding blind single run, 50-run stability profile. All runs successful, all QC checks passed.

## 2. Tier 1 — Curated IPF run

- **Seed set:** 352 curated IPF genes (8 HGNC-remapped)
- **Phase 1 study genes:** 241
- **Phase 4 core genes:** 186
- **Phase 4 significant clusters:** 4
- **Phase 5 survived / eliminated:** 166 / 20
- **Phase 6 meta-significant (random effects):** 155

**Core genes (186, alphabetical):**

- ABCA13
- ABCA3
- ACADL
- ACVRL1
- ADAM17
- ADAMTS14
- AFF3
- AGER
- AKAP13
- ANXA3
- AOC3
- ASPN
- BCHE
- BMP4
- C12orf75
- CA2
- CA4
- CACNA2D2
- CCL21
- CDH13
- CDH2
- CDH3
- CDH5
- CFH
- CHI3L2
- CLCA2
- CLDN1
- CLEC4E
- COL14A1
- COL15A1
- COL17A1
- COL1A1
- COL1A2
- COL3A1
- COL5A1
- COL5A2
- COL6A3
- COMP
- CPA3
- CRTAC1
- CSF3R
- CTSK
- CXCL1
- CXCL12
- CXCL13
- CXCL14
- CXCL8
- DCLK1
- DIO2
- EDNRB
- ENG
- EPHA3
- FAM107A
- FAP
- FCN3
- FGF1
- FKBP5
- FMO5
- FN1
- FNDC1
- FPR1
- GEM
- GPNMB
- GPR4
- GPR87
- GPX3
- GREM2
- GSTA1
- HGF
- HHIP
- HHLA2
- HIF3A
- HMCN1
- HS6ST2
- HSD17B6
- HSPA4L
- ID1
- IFNG
- IGF1
- IGFL2
- IL13RA2
- IL18R1
- IL18RAP
- IL1R2
- IL1RL1
- IL33
- IL6
- ITGA11
- ITGA5
- ITGAV
- ITGB8
- IVD
- KDR
- KRT14
- KRT15
- KRT17
- KRT5
- LAMA3
- LIN7A
- LOXL2
- LPL
- LRRC17
- LRRN1
- LTBP1
- MARCKS
- MGAM
- MME
- MMP1
- MMP10
- MMP13
- MMP2
- MMP3
- MMP7
- MS4A2
- MSMB
- MUC5B
- MXRA5
- NEK11
- NELL2
- NPR1
- OGN
- P2RY1
- PAPSS2
- PDGFRA
- PEBP4
- PECAM1
- PGC
- PLA2G1B
- PLA2G4F
- PLLP
- POSTN
- PPP1R14A
- PROK2
- PRSS12
- PSD3
- PTPRB
- S100A12
- S100A8
- SACK1D
- SCN7A
- SERPINB3
- SERPIND1
- SERPINI2
- SFRP2
- SFRP4
- SFTPA2
- SLAMF7
- SLC27A2
- SLC39A8
- SLC6A4
- SLCO2A1
- SLCO4A1
- SLN
- SMAD6
- SMAD7
- SOX2
- SPARC
- SPP1
- STX11
- SULF1
- SYNPO2
- SYTL2
- TAMALIN
- TDO2
- TGFB2
- THBS1
- THBS2
- THY1
- TIMP2
- TIMP3
- TLR2
- TMEM100
- TMEM45A
- TMPRSS4
- TNNC1
- TP63
- TRIM2
- TSPAN1
- VCAM1
- VIPR1
- VSIG1
- WDR49
- WNT5A
- ZBTB16
- ZBTB7C
- ZKSCAN1

## 3. Tier 2 — Blind protein-coding run

- **Seed set:** 19,296 protein-coding genes
- **Phase 1 study genes:** 4479
- **Phase 4 core genes:** 2451
- **Phase 4 significant clusters:** 24
- **Phase 5 survived / eliminated:** 1856 / 595
- **Phase 6 genes analyzed:** 1856
- **Phase 6 significant (random effects):** 1584

**Top 20 genes by Phase 6 random-effects -log10(p):**

| Rank | Gene | Random effect | p-value | -log10(p) | Direction | Datasets |
|---|---|---|---|---|---|---|
| 1 | ASPN | 2.8061 | 1.73e-109 | 108.76 | up | 3 |
| 2 | SDR16C5 | -1.6775 | 2.12e-83 | 82.67 | down | 3 |
| 3 | BCAT2 | -0.6848 | 3.12e-72 | 71.51 | down | 3 |
| 4 | DENND3 | -0.9571 | 2.40e-70 | 69.62 | down | 3 |
| 5 | CD24 | 2.1520 | 6.43e-69 | 68.19 | up | 3 |
| 6 | LSS | -0.9154 | 2.49e-67 | 66.60 | down | 3 |
| 7 | FCN3 | -2.2836 | 8.88e-67 | 66.05 | down | 3 |
| 8 | HHAT | 0.7281 | 3.92e-66 | 65.41 | up | 3 |
| 9 | KRT5 | 2.4403 | 6.43e-63 | 62.19 | up | 3 |
| 10 | TRIM2 | 1.0712 | 1.42e-60 | 59.85 | up | 3 |
| 11 | RCAN2 | 0.6741 | 1.32e-59 | 58.88 | up | 3 |
| 12 | S100A8 | -1.6453 | 1.58e-57 | 56.80 | down | 3 |
| 13 | SNCAIP | 0.6612 | 3.83e-57 | 56.42 | up | 3 |
| 14 | S1PR1 | -1.0316 | 6.83e-57 | 56.17 | down | 3 |
| 15 | CRABP2 | 0.7957 | 1.55e-55 | 54.81 | up | 3 |
| 16 | CFH | 1.1680 | 1.82e-54 | 53.74 | up | 3 |
| 17 | FAM169A | 0.5970 | 5.47e-54 | 53.26 | up | 3 |
| 18 | PRTFDC1 | 0.8711 | 1.03e-53 | 52.99 | up | 3 |
| 19 | VSNL1 | 0.9305 | 1.64e-53 | 52.78 | up | 3 |
| 20 | PGM2L1 | 0.8669 | 1.98e-52 | 51.70 | up | 3 |

## 4. Tier 3 — 50-run stability profile

- **Runs completed:** 50 / 50 (0 failures)
- **Master seed:** 42
- **Total unique genes seen:** 2545

**Stability classifications:**

| Class | Count | Fraction | Threshold |
|---|---|---|---|
| Iron-clad | 2309 | 90.7% | >= 90% appearance |
| Borderline | 180 | — | 50-89% appearance |
| Stochastic | 56 | — | < 50% appearance |

**Per-run core gene counts:**

| Metric | Value |
|---|---|
| Min | 2416 |
| Max | 2471 |
| Mean | 2442.6 |
| Std dev | 11.42 |

**Pairwise Jaccard similarity:**

| Metric | Value |
|---|---|
| Median | 0.9547 |
| 25th percentile | 0.9512 |
| 75th percentile | 0.9576 |
| Min | 0.9354 |
| Max | 0.9718 |
| Pairs computed | 1225 |

**Total wall clock:** 33272 sec (9.2 hours)

## 5. Iron-clad gene set (n=2309)

Full list: `stability_50run/stability_scores.csv` (filtered by stability_class = iron-clad)

1946 of 2309 iron-clad genes appeared in all 50 runs.

## 6. Cluster organization of iron-clad genes

Analysis script: `scripts/analyze_iron_clad_clusters.py`
Output files: `stability_50run/iron_clad_cluster_analysis.csv`, `stability_50run/iron_clad_cluster_summary.csv`

**204 clusters contain >= 5 iron-clad members.**
Largest cluster: 44 iron-clad genes.
Full cluster details in `iron_clad_cluster_summary.csv`.

## 7. Cross-tier comparison

**Curated core genes (186) AND iron-clad: 179**

**Curated core genes NOT in iron-clad: 7**

- CHI3L2
- CXCL1
- FMO5
- FPR1
- SFTPA2
- TIMP2
- TNNC1

**Blind single core genes (2451) AND iron-clad: 2281**

**Iron-clad genes NOT in blind single run: 28**

## 8. Findings

Candidates identified by objective criteria applied mechanically to run outputs.

Effect size and p-value are deterministic given Phase 5 survivors — Phase 6 meta-analysis math is not stochastic, verified empirically. Cross-run variation appears in n_runs_in_phase6 (how often each iron-clad gene reaches the meta-analysis stage), not in the values themselves. See `iron_clad_aggregated.csv` for the full per-gene audit trail.

### A. Highest median effect size (top 10 iron-clad by median |random_effect|)

| Gene | Median effect | Direction | n_runs |
|---|---|---|---|
| MMP7 | 3.3706 | up | 50 |
| MMP1 | 3.3392 | up | 49 |
| SFRP2 | 3.2288 | up | 50 |
| ASPN | 2.8061 | up | 48 |
| CXCL14 | 2.7692 | up | 50 |
| CP | 2.5752 | up | 50 |
| S100A12 | -2.4679 | down | 50 |
| UBD | 2.4664 | up | 50 |
| SLCO4A1 | -2.4649 | down | 50 |
| BTNL9 | -2.4439 | down | 50 |

### B. Strongest statistical evidence (top 10 iron-clad by -log10(median p))

| Gene | -log10(median p) | Median p | Median effect |
|---|---|---|---|
| ASPN | 108.76 | 1.73e-109 | 2.8061 |
| SDR16C5 | 82.67 | 2.12e-83 | -1.6775 |
| BCAT2 | 71.51 | 3.12e-72 | -0.6848 |
| DENND3 | 69.62 | 2.40e-70 | -0.9571 |
| LSS | 66.60 | 2.49e-67 | -0.9154 |
| FCN3 | 66.05 | 8.88e-67 | -2.2836 |
| HHAT | 65.41 | 3.92e-66 | 0.7281 |
| KRT5 | 62.19 | 6.43e-63 | 2.4403 |
| TRIM2 | 59.85 | 1.42e-60 | 1.0712 |
| RCAN2 | 58.88 | 1.32e-59 | 0.6741 |

### C. Perfect reproducibility (50/50 runs): 1946 genes

1946 iron-clad genes appeared in all 50 runs.

### D. Cluster anchors: 2062 genes

2062 iron-clad genes reside in modal clusters with >= 5 other iron-clad members.

### E. Curated-seed-AND-iron-clad: 184 genes

Of the 352 a priori IPF candidate genes in ipf_curated_genes.csv, 184 appear in the iron-clad set from the 50-run blind stability profile. This is distinct from Section 7 (curated core genes that are also iron-clad). Section 8E captures iron-clad genes that were nominated as IPF candidates by the curated seed list but may not have survived the curated single-run pipeline — useful for identifying boundary-case candidates the engine recovers under blind stability but not under curated single-run conditions.

- ABCA13
- ABCA3
- ACADL
- ACVRL1
- ADAM17
- ADAMTS14
- AFF3
- AGER
- AKAP13
- ANXA3
- AOC3
- ASPN
- BCHE
- BMP4
- C12orf75
- CA2
- CA4
- CACNA2D2
- CCL21
- CDH13
- CDH2
- CDH3
- CDH5
- CFH
- CLCA2
- CLDN1
- CLEC4E
- COL14A1
- COL15A1
- COL17A1
- COL1A1
- COL1A2
- COL3A1
- COL5A1
- COL5A2
- COL6A3
- COMP
- CPA3
- CRTAC1
- CSF3R
- CTSK
- CXCL12
- CXCL13
- CXCL14
- CXCL8
- DCLK1
- DIO2
- EDNRB
- ENG
- EPHA3
- FAM107A
- FAP
- FCN3
- FGF1
- FGF7
- FKBP5
- FN1
- FNDC1
- GEM
- GPNMB
- GPR4
- GPR87
- GPX3
- GREM1
- GREM2
- GSTA1
- HGF
- HHIP
- HHLA2
- HIF3A
- HMCN1
- HS6ST2
- HSD17B6
- HSPA4L
- ID1
- IFNG
- IGF1
- IGFL2
- IL13RA2
- IL18R1
- IL18RAP
- IL1R2
- IL1RL1
- IL1RN
- IL33
- IL6
- ITGA11
- ITGA5
- ITGAV
- ITGB8
- IVD
- JAG2
- KDR
- KRT14
- KRT15
- KRT17
- KRT5
- LAMA3
- LIN7A
- LOXL2
- LPL
- LRRC17
- LRRN1
- LTBP1
- LUM
- MARCKS
- MGAM
- MME
- MMP1
- MMP10
- MMP13
- MMP2
- MMP3
- MMP7
- MS4A2
- MSMB
- MUC5B
- MXRA5
- NEK11
- NELL2
- NPR1
- OGN
- P2RY1
- PAPSS2
- PDGFRA
- PEBP4
- PECAM1
- PGC
- PLA2G1B
- PLA2G4F
- PLLP
- PLOD2
- POSTN
- PPP1R14A
- PROK2
- PRSS12
- PSD3
- PTPRB
- RRM2
- S100A12
- S100A8
- SCN7A
- SERPINB3
- SERPIND1
- SERPINI2
- SFRP2
- SFRP4
- SLAMF7
- SLC27A2
- SLC39A8
- SLC6A4
- SLCO2A1
- SLCO4A1
- SLN
- SMAD6
- SMAD7
- SOX2
- SPARC
- SPP1
- STX11
- SULF1
- SYNPO2
- SYTL2
- TDO2
- TGFB2
- THBS1
- THBS2
- THY1
- TIMP3
- TLR2
- TMEM100
- TMEM45A
- TMPRSS4
- TP63
- TRIM2
- TSPAN1
- VCAM1
- VIPR1
- VSIG1
- WDR49
- WNT5A
- ZBTB16
- ZBTB7C
- ZKSCAN1

### F. Phase-4-stable but Phase-5-unstable genes

548 iron-clad genes are consistently identified in Phase 4 clustering but consistently eliminated at Phase 5 (replication tier).

Full list (548 genes): see `iron_clad_aggregated.csv` where n_runs_in_phase6 = 0.

## 9. Comparison to historical v0.3.2

### Curated single-run comparison

v0.3.2 historical: 190 core genes (from /home/kai001/ipf_validation/, date unknown, single curated run)
v0.3.3.2 current: 186 core genes
Stable across versions: 186 (97.9%)
Dropped under v0.3.3.2: 4 (FGF7, GREM1, LUM, PLOD2)
Newly identified: 0

The v0.3.3.2 dilution fix (Round 4 of v0.3.3 release) removed 4 boundary-case genes from the curated single-run output. The other 186 curated core genes are stable across the methodology change. No genes were newly identified by v0.3.3.2 that weren't already in v0.3.2.

### Stability profile comparison

No v0.3.2 stability profile was produced for IPF. The 50-run blind stability profiling protocol was introduced during the v0.3.3 release cycle; the v0.3.3.2 IPF stability profile (2,309 iron-clad genes) is the inaugural stability characterization for this disease. Historical v0.3.2 data exists only for the curated single run (above), not for blind stability profiling.

## 10. Cluster-level findings

204 clusters with >= 5 iron-clad members. Full membership in `iron_clad_cluster_summary.csv`.

Top 10 largest clusters by iron-clad membership:

| Cluster ID | Iron-clad members | Min modal count | Max modal count |
|---|---|---|---|
| 0 | 44 | 2 | 4 |
| 90 | 40 | 2 | 3 |
| 437 | 25 | 2 | 3 |
| 369 | 24 | 1 | 2 |
| 187 | 23 | 2 | 3 |
| 341 | 23 | 2 | 2 |
| 8 | 22 | 2 | 4 |
| 140 | 22 | 2 | 3 |
| 233 | 22 | 2 | 3 |
| 275 | 22 | 2 | 3 |

## 11. Observations

- 1946 of 2309 iron-clad genes (84.3%) appeared in all 50 runs.
- Per-run core gene count ranged from 2416 to 2471 (std dev 11.42).
- Pairwise Jaccard minimum was 0.9354, meaning even the two most dissimilar runs shared 93.5% of their union.
- 7 of 186 curated core genes are not in the iron-clad set.
- 28 iron-clad genes are not in the Tier 2 blind single-run core gene set.
- 548 iron-clad genes (23.7%) have 0 Phase 6 appearances (eliminated at Phase 5 in every run they appeared in).
- 204 clusters contain >= 5 iron-clad members each.
- 0 iron-clad genes had direction concordance below 0.80.

---

*Draft generated by Kai (Claude Code agent), 2026-04-30. Pending review by Cody Sigmon.*