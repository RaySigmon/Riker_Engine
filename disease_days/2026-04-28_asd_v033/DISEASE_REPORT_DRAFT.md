# ASD Disease Report — v0.3.3 (DRAFT)

---

## 1. Run summary

- **Disease:** ASD (Autism Spectrum Disorder)
- **Date:** 2026-04-28
- **Engine version:** 0.3.3
- **Engine commit:** 1b250ed5c7fb436b25292971911fe743a3ff0339
- **Master seed:** 42
- **Tier 1 wall clock:** 2 min 4 sec
- **Tier 2 wall clock:** 7 min 14 sec
- **Tier 3 wall clock:** 18011 sec (5.0 hours)
- Three tiers completed: curated single run, protein-coding blind single run, 50-run stability profile. All runs successful, all QC checks passed.

## 2. Tier 1 — Curated SFARI run

- **Seed set:** 1,267 SFARI genes
- **Phase 1 study genes:** 141
- **Phase 4 core genes:** 29
- **Phase 4 significant clusters:** 2
- **Phase 5 survived / eliminated:** 29 / 0
- **Phase 6 meta-significant (random effects):** 10

**Core genes (29, alphabetical):**

- ACTL6B
- APBA2
- ATP1A1
- ATP2B2
- DPP10
- DPP6
- EFR3A
- ELP2
- FBXO33
- GABRG2
- ICA1
- ITPR1
- KIAA0232
- LIN7B
- NAV2
- NPAS3
- NPTN
- PPP3CA
- RGS7
- RPH3A
- SLC12A5
- SLC25A12
- SLC25A27
- SLC45A1
- SNX14
- SRGAP3
- TBR1
- UNC80
- ZFHX3

## 3. Tier 2 — Blind protein-coding run

- **Seed set:** 19,296 protein-coding genes
- **Phase 1 study genes:** 1793
- **Phase 4 core genes:** 421
- **Phase 4 significant clusters:** 15
- **Phase 5 survived / eliminated:** 414 / 7
- **Phase 6 genes analyzed:** 414
- **Phase 6 significant (random effects):** 233

**Top 20 genes by Phase 6 random-effects -log10(p):**

| Rank | Gene | Random effect | p-value | -log10(p) | Direction | Datasets |
|---|---|---|---|---|---|---|
| 1 | SFMBT2 | 0.5058 | 1.87e-11 | 10.73 | up | 3 |
| 2 | PLXDC2 | 0.4740 | 1.33e-10 | 9.88 | up | 3 |
| 3 | H2BC21 | 0.5707 | 3.67e-10 | 9.44 | up | 3 |
| 4 | DTNA | 0.5020 | 2.44e-08 | 7.61 | up | 3 |
| 5 | SERF2 | 0.2776 | 2.65e-08 | 7.58 | up | 3 |
| 6 | WWC1 | 0.4237 | 2.95e-08 | 7.53 | up | 3 |
| 7 | ATXN7 | 0.2251 | 3.94e-08 | 7.40 | up | 2 |
| 8 | RPL36 | 0.1914 | 5.24e-08 | 7.28 | up | 3 |
| 9 | SHD | -0.5873 | 5.58e-08 | 7.25 | down | 3 |
| 10 | EPB41L2 | 0.3064 | 6.17e-08 | 7.21 | up | 3 |
| 11 | CPNE3 | 0.3846 | 7.25e-08 | 7.14 | up | 3 |
| 12 | APBA2 | -0.2468 | 1.23e-07 | 6.91 | down | 3 |
| 13 | PITPNC1 | 0.3616 | 1.46e-07 | 6.83 | up | 3 |
| 14 | H2AC20 | 0.4056 | 3.03e-07 | 6.52 | up | 3 |
| 15 | GBP2 | 0.9275 | 4.20e-07 | 6.38 | up | 3 |
| 16 | CPLX2 | -0.4373 | 4.25e-07 | 6.37 | down | 3 |
| 17 | CTSH | 0.5061 | 6.04e-07 | 6.22 | up | 3 |
| 18 | QRICH1 | -0.1787 | 6.90e-07 | 6.16 | down | 3 |
| 19 | HEY2 | 0.4411 | 8.08e-07 | 6.09 | up | 3 |
| 20 | PLOD1 | 0.4091 | 9.77e-07 | 6.01 | up | 3 |

## 4. Tier 3 — 50-run stability profile

- **Runs completed:** 50 / 50 (0 failures)
- **Master seed:** 42
- **Total unique genes seen:** 457

**Stability classifications:**

| Class | Count | Fraction | Threshold |
|---|---|---|---|
| Iron-clad | 394 | 86.2% | >= 90% appearance |
| Borderline | 35 | — | 50-89% appearance |
| Stochastic | 28 | — | < 50% appearance |

**Per-run core gene counts:**

| Metric | Value |
|---|---|
| Min | 410 |
| Max | 434 |
| Mean | 420.5 |
| Std dev | 5.53 |

**Pairwise Jaccard similarity:**

| Metric | Value |
|---|---|
| Median | 0.933 |
| 25th percentile | 0.9238 |
| 75th percentile | 0.9421 |
| Min | 0.8826 |
| Max | 0.9725 |
| Pairs computed | 1225 |

**Total wall clock:** 18011 sec (5.0 hours)

## 5. Iron-clad gene set (n=394)

| Gene | Appearance count |
|---|---|
| AARS1 | 50/50 |
| ABCA1 | 50/50 |
| ABCA5 | 50/50 |
| ABCG4 | 50/50 |
| ACTL6B | 49/50 |
| ACTR1B | 48/50 |
| ACTR6 | 50/50 |
| ADAP2 | 50/50 |
| ADGRG1 | 50/50 |
| AEBP1 | 50/50 |
| AGGF1 | 49/50 |
| ALDOC | 50/50 |
| ANAPC5 | 50/50 |
| ANKRD29 | 50/50 |
| ANO6 | 47/50 |
| APBA2 | 48/50 |
| APOC1 | 45/50 |
| ARHGAP26 | 50/50 |
| ARL3 | 50/50 |
| ARMCX5 | 50/50 |
| ASB13 | 50/50 |
| ASXL2 | 50/50 |
| ATP1A1 | 49/50 |
| ATP5F1A | 50/50 |
| ATP5F1B | 50/50 |
| ATP6V1A | 50/50 |
| ATP6V1C1 | 50/50 |
| ATPAF1 | 50/50 |
| ATRNL1 | 50/50 |
| ATXN7 | 50/50 |
| AVPI1 | 50/50 |
| B2M | 48/50 |
| B4GALT6 | 50/50 |
| BAG3 | 50/50 |
| BIRC3 | 48/50 |
| BMP7 | 50/50 |
| BST2 | 50/50 |
| BTG3 | 50/50 |
| C1QB | 50/50 |
| CACHD1 | 50/50 |
| CADPS | 50/50 |
| CBLN2 | 49/50 |
| CBX7 | 50/50 |
| CCDC25 | 47/50 |
| CCL20 | 50/50 |
| CCN1 | 50/50 |
| CCNB1IP1 | 45/50 |
| CD14 | 50/50 |
| CD84 | 50/50 |
| CDCA5 | 49/50 |
| CDKN1A | 50/50 |
| CEBPD | 47/50 |
| CEND1 | 50/50 |
| CHCHD4 | 49/50 |
| CHCHD7 | 50/50 |
| CHST11 | 50/50 |
| CHSY1 | 50/50 |
| CIRBP | 50/50 |
| CNN3 | 50/50 |
| CNTNAP1 | 50/50 |
| COLGALT1 | 45/50 |
| COPG1 | 50/50 |
| COQ3 | 49/50 |
| COQ6 | 50/50 |
| CORO2A | 50/50 |
| COX7A1 | 50/50 |
| CPEB1 | 50/50 |
| CPLX2 | 50/50 |
| CPNE3 | 50/50 |
| CTSH | 50/50 |
| CXCL16 | 50/50 |
| CYBA | 50/50 |
| CYC1 | 50/50 |
| DDAH2 | 50/50 |
| DDX1 | 50/50 |
| DDX24 | 50/50 |
| DEPDC5 | 50/50 |
| DHX35 | 47/50 |
| DKK3 | 50/50 |
| DNAJB1 | 50/50 |
| DOK1 | 50/50 |
| DOP1A | 50/50 |
| DPP6 | 50/50 |
| DTD1 | 50/50 |
| DTNA | 50/50 |
| DTX2 | 50/50 |
| DYNLT1 | 50/50 |
| DYNLT5 | 50/50 |
| EFR3A | 50/50 |
| EHD4 | 50/50 |
| EID1 | 50/50 |
| EIF3K | 48/50 |
| ELF1 | 50/50 |
| ELP2 | 49/50 |
| EMC1 | 50/50 |
| ENO2 | 50/50 |
| EPB41L2 | 50/50 |
| EXOSC9 | 50/50 |
| EXTL2 | 50/50 |
| FAHD2B | 50/50 |
| FAM111A | 50/50 |
| FAM216A | 46/50 |
| FAM234B | 50/50 |
| FAM241B | 50/50 |
| FAM81A | 50/50 |
| FAR2 | 49/50 |
| FBLN7 | 50/50 |
| FERMT2 | 50/50 |
| FH | 49/50 |
| FOCAD | 50/50 |
| FXR2 | 50/50 |
| GABRA1 | 49/50 |
| GABRG2 | 50/50 |
| GALNT17 | 48/50 |
| GAREM1 | 50/50 |
| GBP2 | 50/50 |
| GLMN | 50/50 |
| GLS | 48/50 |
| GLS2 | 50/50 |
| GMPR2 | 50/50 |
| GNA12 | 48/50 |
| GNG5 | 50/50 |
| GOT1 | 50/50 |
| GOT2 | 50/50 |
| GPI | 49/50 |
| GPRASP2 | 50/50 |
| H2AC20 | 50/50 |
| H2BC21 | 49/50 |
| HACL1 | 50/50 |
| HAGH | 48/50 |
| HEATR5B | 49/50 |
| HEY2 | 50/50 |
| HK2 | 50/50 |
| HMCES | 50/50 |
| HMG20A | 48/50 |
| HMGB3 | 48/50 |
| HPS6 | 49/50 |
| HSPA1A | 50/50 |
| HSPA5 | 50/50 |
| HSPB1 | 50/50 |
| HYLS1 | 46/50 |
| ICA1 | 50/50 |
| IDH1 | 50/50 |
| IDH3A | 50/50 |
| IDH3B | 50/50 |
| IFITM3 | 50/50 |
| IFRD1 | 50/50 |
| INA | 50/50 |
| IRF8 | 45/50 |
| ITPR1 | 50/50 |
| JUN | 50/50 |
| KCNN3 | 50/50 |
| KCNS1 | 49/50 |
| KCNS2 | 50/50 |
| KCTD9 | 50/50 |
| KIFBP | 50/50 |
| KLHDC2 | 49/50 |
| KLHDC9 | 50/50 |
| KRT222 | 50/50 |
| KYAT3 | 50/50 |
| LAPTM5 | 50/50 |
| LCMT1 | 50/50 |
| LHX2 | 50/50 |
| LIMS1 | 50/50 |
| LIN7B | 47/50 |
| LPCAT4 | 49/50 |
| LPP | 49/50 |
| LYN | 48/50 |
| LYRM9 | 48/50 |
| MAGED2 | 50/50 |
| MAGEE1 | 50/50 |
| MANF | 50/50 |
| MAP3K9 | 50/50 |
| MAPKAPK3 | 50/50 |
| MARF1 | 50/50 |
| MAT2B | 50/50 |
| MCL1 | 50/50 |
| MDFIC | 49/50 |
| ME3 | 50/50 |
| MGST1 | 50/50 |
| MICU3 | 50/50 |
| MPP1 | 50/50 |
| MPZL1 | 50/50 |
| MRPL2 | 49/50 |
| MRPL22 | 45/50 |
| MSI2 | 50/50 |
| MSN | 50/50 |
| MTMR7 | 50/50 |
| MTRES1 | 48/50 |
| N4BP2L1 | 50/50 |
| NAE1 | 50/50 |
| NAPB | 50/50 |
| NARS1 | 49/50 |
| NCKIPSD | 50/50 |
| NDEL1 | 50/50 |
| NDRG3 | 48/50 |
| NDRG4 | 48/50 |
| NDUFAF5 | 50/50 |
| NDUFB8 | 50/50 |
| NECAP2 | 50/50 |
| NEFH | 50/50 |
| NET1 | 50/50 |
| NFE2L2 | 50/50 |
| NFKBIZ | 50/50 |
| NFS1 | 50/50 |
| NGFR | 49/50 |
| NIF3L1 | 50/50 |
| NME4 | 50/50 |
| NMNAT2 | 50/50 |
| NOMO1 | 49/50 |
| NPAS3 | 50/50 |
| NPTN | 48/50 |
| NQO1 | 50/50 |
| NR3C2 | 50/50 |
| NSG1 | 50/50 |
| NTAQ1 | 48/50 |
| NTHL1 | 50/50 |
| NTN4 | 50/50 |
| OBI1 | 50/50 |
| OGDHL | 50/50 |
| OLA1 | 50/50 |
| ORC2 | 50/50 |
| P3H2 | 50/50 |
| P4HA1 | 50/50 |
| PABPC4 | 50/50 |
| PALLD | 50/50 |
| PARL | 48/50 |
| PARN | 46/50 |
| PARP14 | 48/50 |
| PARP2 | 50/50 |
| PARP9 | 50/50 |
| PDIA5 | 50/50 |
| PDIA6 | 50/50 |
| PDYN | 50/50 |
| PFKM | 50/50 |
| PFKP | 50/50 |
| PIGZ | 49/50 |
| PIK3AP1 | 50/50 |
| PIK3R4 | 50/50 |
| PIP4P2 | 48/50 |
| PITPNC1 | 50/50 |
| PLIN2 | 50/50 |
| PLOD1 | 50/50 |
| PLPP6 | 48/50 |
| PLXDC2 | 49/50 |
| PNP | 50/50 |
| POMGNT2 | 50/50 |
| POU3F2 | 50/50 |
| PPP1R14B | 50/50 |
| PPP3CA | 50/50 |
| PRC1 | 48/50 |
| PRKCE | 50/50 |
| PRPF19 | 50/50 |
| PRR14 | 50/50 |
| PSME2 | 50/50 |
| PSMG3 | 50/50 |
| PTAFR | 47/50 |
| PTH2R | 50/50 |
| PTPN13 | 50/50 |
| PTTG1IP | 50/50 |
| PYGB | 50/50 |
| QRICH1 | 50/50 |
| RAB29 | 50/50 |
| RAB31 | 50/50 |
| RALGPS1 | 50/50 |
| RASGEF1B | 50/50 |
| RCAN2 | 50/50 |
| RCOR3 | 50/50 |
| RDH10 | 50/50 |
| RETREG2 | 49/50 |
| RFX5 | 50/50 |
| RHBDF2 | 50/50 |
| RHBDL3 | 50/50 |
| RHPN2 | 50/50 |
| RIOK1 | 50/50 |
| RMND1 | 50/50 |
| RNMT | 49/50 |
| RPL12 | 50/50 |
| RPL36 | 47/50 |
| RPL39 | 47/50 |
| RPLP1 | 50/50 |
| RPS16 | 50/50 |
| RPS2 | 50/50 |
| RWDD2A | 50/50 |
| RWDD2B | 50/50 |
| S100A10 | 49/50 |
| S100A8 | 50/50 |
| SAYSD1 | 50/50 |
| SBF2 | 50/50 |
| SCN1B | 50/50 |
| SCN2B | 50/50 |
| SCN4B | 49/50 |
| SDF2L1 | 50/50 |
| SDHA | 50/50 |
| SEC11A | 50/50 |
| SEC14L2 | 50/50 |
| SEC61A1 | 47/50 |
| SERF2 | 50/50 |
| SERGEF | 50/50 |
| SERP1 | 49/50 |
| SERPINH1 | 50/50 |
| SERTAD1 | 46/50 |
| SEZ6L2 | 50/50 |
| SFMBT2 | 50/50 |
| SH3GL2 | 49/50 |
| SH3GLB1 | 50/50 |
| SHLD1 | 50/50 |
| SIPA1 | 50/50 |
| SLA | 50/50 |
| SLC12A5 | 50/50 |
| SLC15A3 | 50/50 |
| SLC16A14 | 50/50 |
| SLC16A3 | 50/50 |
| SLC25A1 | 47/50 |
| SLC25A12 | 50/50 |
| SLC25A27 | 50/50 |
| SLC25A3 | 50/50 |
| SLC25A37 | 48/50 |
| SLC2A5 | 50/50 |
| SLC33A2 | 50/50 |
| SLC35B3 | 50/50 |
| SLC37A4 | 50/50 |
| SLC45A1 | 50/50 |
| SLC6A15 | 46/50 |
| SNAP25 | 50/50 |
| SNAP47 | 50/50 |
| SNU13 | 50/50 |
| SNX14 | 50/50 |
| SNX17 | 49/50 |
| SOX2 | 50/50 |
| SOX4 | 50/50 |
| SP110 | 50/50 |
| SPATA7 | 50/50 |
| SPIN2B | 50/50 |
| SPINT2 | 50/50 |
| SPOCK2 | 50/50 |
| SRGN | 50/50 |
| STAT4 | 47/50 |
| SUPV3L1 | 48/50 |
| SUV39H1 | 50/50 |
| SYNGR1 | 50/50 |
| TARBP1 | 50/50 |
| TBC1D4 | 50/50 |
| THEMIS | 49/50 |
| TIMP1 | 50/50 |
| TIPARP | 45/50 |
| TMEM268 | 50/50 |
| TMIGD3 | 48/50 |
| TMX3 | 48/50 |
| TMX4 | 45/50 |
| TOB1 | 49/50 |
| TOX | 50/50 |
| TRAF3IP2 | 46/50 |
| TRAM2 | 50/50 |
| TRIM44 | 50/50 |
| TRMT11 | 50/50 |
| TRPC1 | 50/50 |
| TSNAX | 50/50 |
| TTC1 | 45/50 |
| TTC13 | 50/50 |
| TTC19 | 50/50 |
| TUBB2B | 50/50 |
| TXNIP | 47/50 |
| UBE2L6 | 49/50 |
| UBE4A | 50/50 |
| UBLCP1 | 48/50 |
| UBTD1 | 48/50 |
| UCHL5 | 50/50 |
| UGDH | 50/50 |
| UQCRC1 | 45/50 |
| UQCRC2 | 50/50 |
| URGCP | 48/50 |
| VAMP1 | 46/50 |
| VDAC3 | 46/50 |
| VIM | 50/50 |
| VPS16 | 48/50 |
| VPS33B | 50/50 |
| VPS36 | 50/50 |
| VPS39 | 50/50 |
| VRK1 | 49/50 |
| VSNL1 | 49/50 |
| WDR54 | 50/50 |
| WSCD1 | 50/50 |
| WWC1 | 50/50 |
| YAP1 | 50/50 |
| YBX1 | 46/50 |
| YBX3 | 50/50 |
| YPEL4 | 49/50 |
| ZBTB24 | 50/50 |
| ZFAND3 | 50/50 |
| ZFHX3 | 50/50 |
| ZFYVE21 | 50/50 |
| ZNF25 | 46/50 |
| ZNF385B | 50/50 |

## 6. Cluster organization of iron-clad genes

Analysis script: `scripts/analyze_iron_clad_clusters.py`
Output files: `stability_50run/iron_clad_cluster_analysis.csv`, `stability_50run/iron_clad_cluster_summary.csv`

**Clusters with >= 5 iron-clad members (by modal cluster assignment):**

| Cluster ID | Iron-clad members | Min modal count | Max modal count | Mean modal count |
|---|---|---|---|---|
| 0 | 14 | 3 | 9 | 6.9 |
| 2 | 5 | 4 | 4 | 4.0 |
| 4 | 14 | 3 | 7 | 5.0 |
| 7 | 6 | 3 | 3 | 3.0 |
| 10 | 9 | 2 | 3 | 2.7 |
| 17 | 12 | 2 | 3 | 2.5 |
| 20 | 6 | 2 | 4 | 3.7 |
| 21 | 13 | 2 | 3 | 2.6 |
| 23 | 7 | 3 | 3 | 3.0 |
| 24 | 7 | 3 | 3 | 3.0 |
| 26 | 10 | 3 | 3 | 3.0 |
| 27 | 5 | 3 | 3 | 3.0 |
| 29 | 8 | 2 | 2 | 2.0 |
| 30 | 5 | 2 | 2 | 2.0 |
| 34 | 14 | 2 | 2 | 2.0 |
| 35 | 9 | 3 | 3 | 3.0 |
| 36 | 6 | 3 | 3 | 3.0 |
| 46 | 10 | 3 | 3 | 3.0 |
| 52 | 8 | 2 | 2 | 2.0 |
| 57 | 8 | 2 | 2 | 2.0 |
| 68 | 8 | 2 | 3 | 2.1 |
| 75 | 10 | 2 | 2 | 2.0 |
| 79 | 6 | 3 | 3 | 3.0 |
| 89 | 6 | 2 | 2 | 2.0 |
| 93 | 6 | 2 | 3 | 2.5 |
| 124 | 14 | 2 | 3 | 2.4 |
| 139 | 11 | 2 | 3 | 2.3 |
| 140 | 5 | 3 | 3 | 3.0 |
| 149 | 6 | 3 | 3 | 3.0 |
| 150 | 5 | 3 | 3 | 3.0 |
| 152 | 5 | 3 | 3 | 3.0 |
| 164 | 5 | 3 | 3 | 3.0 |
| 167 | 10 | 3 | 3 | 3.0 |
| 168 | 5 | 2 | 2 | 2.0 |
| 173 | 6 | 3 | 3 | 3.0 |

## 7. Cross-tier comparison

**Curated core genes (29) AND iron-clad: 19**

- ACTL6B
- APBA2
- ATP1A1
- DPP6
- EFR3A
- ELP2
- GABRG2
- ICA1
- ITPR1
- LIN7B
- NPAS3
- NPTN
- PPP3CA
- SLC12A5
- SLC25A12
- SLC25A27
- SLC45A1
- SNX14
- ZFHX3

**Curated core genes NOT in iron-clad: 10**

- ATP2B2
- DPP10
- FBXO33
- KIAA0232
- NAV2
- RGS7
- RPH3A
- SRGAP3
- TBR1
- UNC80

**Blind single core genes (421) AND iron-clad: 389**

**Iron-clad genes NOT in blind single run: 5**

- COLGALT1
- SERTAD1
- TMX4
- TRAF3IP2
- YBX1

## 8. Findings

Candidates identified by objective criteria applied mechanically to run outputs.

### A. Highest meta-analysis effect size (top 10 iron-clad by |random_effect|)

| Gene | Random effect | Direction |
|---|---|---|
| BAG3 | 1.4722 | up |
| CCN1 | 1.4007 | up |
| HSPB1 | 1.2985 | up |
| TIMP1 | 1.1720 | up |
| C1QB | 1.1328 | up |
| PDYN | 1.1264 | up |
| RDH10 | 1.0550 | up |
| S100A8 | 1.0438 | up |
| NQO1 | 1.0435 | up |
| YBX3 | 1.0298 | up |

### B. Strongest statistical confidence (top 10 iron-clad by -log10(p_random))

| Gene | -log10(p) | p-value |
|---|---|---|
| SFMBT2 | 10.73 | 1.87e-11 |
| PLXDC2 | 9.88 | 1.33e-10 |
| H2BC21 | 9.44 | 3.67e-10 |
| DTNA | 7.61 | 2.44e-08 |
| SERF2 | 7.58 | 2.65e-08 |
| WWC1 | 7.53 | 2.95e-08 |
| ATXN7 | 7.40 | 3.94e-08 |
| RPL36 | 7.28 | 5.24e-08 |
| EPB41L2 | 7.21 | 6.17e-08 |
| CPNE3 | 7.14 | 7.25e-08 |

### C. Perfect reproducibility (50/50 runs): 293 genes

293 iron-clad genes appeared in all 50 runs.

### D. Cluster anchors: 284 genes

284 iron-clad genes reside in modal clusters with >= 5 other iron-clad members.

### E. SFARI-AND-iron-clad: 28 genes

- ACTL6B
- APBA2
- ATP1A1
- CADPS
- DEPDC5
- DPP6
- EFR3A
- ELP2
- GABRG2
- ICA1
- ITPR1
- LHX2
- LIN7B
- NPAS3
- NPTN
- NR3C2
- PPP3CA
- PRPF19
- QRICH1
- SEZ6L2
- SLC12A5
- SLC25A12
- SLC25A27
- SLC45A1
- SNAP25
- SNX14
- ZFHX3
- ZNF385B

## 9. Comparison to historical v0.3.2

- **Source:** `stability_ASD_blind/stability_scores.csv`
- **v0.3.2 iron-clad count:** 376
- **v0.3.3 iron-clad count:** 394
- **Intersection (in both):** 331
- **v0.3.2 only:** 45
- **v0.3.3 only:** 63

**Genes in v0.3.2 iron-clad but NOT v0.3.3 iron-clad:**

- ACBD6
- ANXA5
- ARHGAP45
- ATP2B2
- ATP6V1E1
- C1QC
- C3orf62
- DEPP1
- DPP10
- DPP8
- ENTPD6
- ERBB2
- FAM131B
- FIBP
- GADD45B
- GNAI3
- GPAT3
- GPCPD1
- HSPD1
- IDH3G
- KAT7
- KIAA0232
- MAPKAPK2
- NAP1L5
- NBPF9
- NEK6
- NIBAN2
- PEA15
- PELI3
- PHYH
- PPL
- PPP1R3F
- PVALB
- RAMP1
- RFXANK
- RNF41
- RPL18A
- RPS19
- SERTAD4
- SHD
- SLC35B4
- SYT13
- TOPBP1
- VCPKMT
- ZFP36L1

**Genes in v0.3.3 iron-clad but NOT v0.3.2 iron-clad:**

- ADGRG1
- AGGF1
- ANKRD29
- APOC1
- BAG3
- CADPS
- CBLN2
- CBX7
- CCL20
- CCN1
- CCNB1IP1
- CNTNAP1
- CYC1
- DKK3
- DOP1A
- DPP6
- ELP2
- GAREM1
- GLS
- HMCES
- HPS6
- HSPA1A
- HSPB1
- KCNN3
- KCNS1
- MSI2
- NARS1
- NDEL1
- NDRG4
- NDUFB8
- NPTN
- OBI1
- PDIA5
- PIGZ
- PIP4P2
- PLPP6
- PRR14
- QRICH1
- RAB29
- RASGEF1B
- RCOR3
- RIOK1
- SEC61A1
- SERP1
- SERPINH1
- SEZ6L2
- SHLD1
- SLC25A1
- SLC2A5
- SLC33A2
- SLC37A4
- SNX14
- SNX17
- SUPV3L1
- TMIGD3
- TMX3
- TRPC1
- TTC1
- UBTD1
- VPS16
- WSCD1
- YPEL4
- ZNF25

## 10. Cluster-level findings

Clusters with >= 5 iron-clad members (by modal assignment across 50 runs):

### Cluster 0 (14 iron-clad members)

Min modal count: 3, Max: 9, Mean: 6.9

Members:
- ATXN7
- BMP7
- DDAH2
- DTNA
- EPB41L2
- GOT2
- MCL1
- PDIA6
- RPS2
- SLC16A14
- SNAP25
- WWC1
- ZFHX3
- ZFYVE21

### Cluster 2 (5 iron-clad members)

Min modal count: 4, Max: 4, Mean: 4.0

Members:
- CHST11
- HSPA5
- P3H2
- SBF2
- SERF2

### Cluster 4 (14 iron-clad members)

Min modal count: 3, Max: 7, Mean: 5.0

Members:
- AARS1
- C1QB
- CDKN1A
- FXR2
- GBP2
- IFITM3
- NDEL1
- NQO1
- QRICH1
- RCOR3
- SNAP47
- SPIN2B
- TIMP1
- YBX3

### Cluster 7 (6 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- IFRD1
- P4HA1
- RPLP1
- SEC11A
- SLA
- UGDH

### Cluster 10 (9 iron-clad members)

Min modal count: 2, Max: 3, Mean: 2.7

Members:
- ATP5F1A
- FAHD2B
- H2AC20
- NFE2L2
- PFKM
- PITPNC1
- PPP1R14B
- SEC14L2
- SOX2

### Cluster 17 (12 iron-clad members)

Min modal count: 2, Max: 3, Mean: 2.5

Members:
- ALDOC
- COLGALT1
- COQ6
- DEPDC5
- NGFR
- PIGZ
- RPL36
- RPL39
- SERP1
- SPOCK2
- TRIM44
- UBE2L6

### Cluster 20 (6 iron-clad members)

Min modal count: 2, Max: 4, Mean: 3.7

Members:
- FAM216A
- KCTD9
- ME3
- MTMR7
- TTC19
- WDR54

### Cluster 21 (13 iron-clad members)

Min modal count: 2, Max: 3, Mean: 2.6

Members:
- AEBP1
- CACHD1
- CHCHD7
- EHD4
- FERMT2
- MAGED2
- N4BP2L1
- NAE1
- RIOK1
- SLC15A3
- SLC33A2
- SNX14
- VPS39

### Cluster 23 (7 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- ABCA1
- BST2
- CD14
- LAPTM5
- MSN
- PDYN
- YAP1

### Cluster 24 (7 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- BAG3
- CNN3
- CYBA
- HSPA1A
- PLIN2
- PNP
- RDH10

### Cluster 26 (10 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- ADAP2
- CPNE3
- DYNLT1
- H2BC21
- IDH1
- LIMS1
- PARP9
- PLXDC2
- SRGN
- TRAM2

### Cluster 27 (5 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- CXCL16
- GNG5
- SDF2L1
- SOX4
- TUBB2B

### Cluster 29 (8 iron-clad members)

Min modal count: 2, Max: 2, Mean: 2.0

Members:
- CHSY1
- ELF1
- FAM111A
- HEY2
- MANF
- NECAP2
- PLOD1
- RAB29

### Cluster 30 (5 iron-clad members)

Min modal count: 2, Max: 2, Mean: 2.0

Members:
- CPEB1
- MRPL22
- SUV39H1
- UBE4A
- URGCP

### Cluster 34 (14 iron-clad members)

Min modal count: 2, Max: 2, Mean: 2.0

Members:
- ABCG4
- ACTL6B
- AVPI1
- CDCA5
- COPG1
- COX7A1
- DDX24
- GLMN
- KLHDC9
- MPP1
- SCN4B
- SLC25A27
- TBC1D4
- UCHL5

### Cluster 35 (9 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- COQ3
- GPRASP2
- MAP3K9
- NIF3L1
- PRPF19
- SHLD1
- TOX
- UBLCP1
- VPS36

### Cluster 36 (6 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- CBX7
- DPP6
- GALNT17
- LIN7B
- SNX17
- YPEL4

### Cluster 46 (10 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- ANO6
- CTSH
- GNA12
- HK2
- LHX2
- LYN
- MAPKAPK3
- PSME2
- TIPARP
- WSCD1

### Cluster 52 (8 iron-clad members)

Min modal count: 2, Max: 2, Mean: 2.0

Members:
- B2M
- BIRC3
- CCDC25
- DHX35
- HMG20A
- PARL
- PARN
- TMIGD3

### Cluster 57 (8 iron-clad members)

Min modal count: 2, Max: 2, Mean: 2.0

Members:
- ARHGAP26
- CEND1
- EXOSC9
- NR3C2
- RMND1
- RNMT
- SERGEF
- ZNF385B

### Cluster 68 (8 iron-clad members)

Min modal count: 2, Max: 3, Mean: 2.1

Members:
- ACTR1B
- HAGH
- IDH3A
- NEFH
- OGDHL
- OLA1
- PYGB
- TOB1

### Cluster 75 (10 iron-clad members)

Min modal count: 2, Max: 2, Mean: 2.0

Members:
- ARL3
- CORO2A
- ELP2
- FH
- NDUFB8
- SAYSD1
- SLC25A12
- STAT4
- TARBP1
- VAMP1

### Cluster 79 (6 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- CYC1
- DOP1A
- HACL1
- MARF1
- NDUFAF5
- PIK3R4

### Cluster 89 (6 iron-clad members)

Min modal count: 2, Max: 2, Mean: 2.0

Members:
- AGGF1
- APOC1
- CCNB1IP1
- CEBPD
- IRF8
- TTC1

### Cluster 93 (6 iron-clad members)

Min modal count: 2, Max: 3, Mean: 2.5

Members:
- CCN1
- LYRM9
- PIP4P2
- SERPINH1
- TRMT11
- VIM

### Cluster 124 (14 iron-clad members)

Min modal count: 2, Max: 3, Mean: 2.4

Members:
- ANKRD29
- B4GALT6
- EID1
- FAM234B
- INA
- KCNS2
- KYAT3
- NTN4
- PPP3CA
- PRKCE
- SLC6A15
- TSNAX
- UQCRC2
- ZNF25

### Cluster 139 (11 iron-clad members)

Min modal count: 2, Max: 3, Mean: 2.3

Members:
- ATP6V1A
- CBLN2
- EXTL2
- FAM241B
- GABRA1
- KCNS1
- KRT222
- MICU3
- NMNAT2
- RHBDL3
- TMX4

### Cluster 140 (5 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- FAM81A
- GLS2
- ICA1
- MAGEE1
- UQCRC1

### Cluster 149 (6 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- HMCES
- LPCAT4
- ORC2
- SLC37A4
- SNU13
- ZBTB24

### Cluster 150 (5 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- MPZL1
- NPAS3
- PABPC4
- POU3F2
- RPS16

### Cluster 152 (5 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- FBLN7
- GABRG2
- GOT1
- RCAN2
- SH3GL2

### Cluster 164 (5 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- APBA2
- GMPR2
- NFS1
- PRC1
- SDHA

### Cluster 167 (10 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- BTG3
- CD84
- CHCHD4
- GLS
- JUN
- KLHDC2
- MRPL2
- NPTN
- SFMBT2
- TMX3

### Cluster 168 (5 iron-clad members)

Min modal count: 2, Max: 2, Mean: 2.0

Members:
- ARMCX5
- ENO2
- IDH3B
- PARP2
- RFX5

### Cluster 173 (6 iron-clad members)

Min modal count: 3, Max: 3, Mean: 3.0

Members:
- ATP1A1
- ITPR1
- PTH2R
- SCN2B
- THEMIS
- VSNL1

## 11. Observations

- 293 of 394 iron-clad genes (74.4%) appeared in all 50 runs.
- 101 iron-clad genes appeared in 45-49 runs, placing them near the iron-clad threshold.
- Per-run core gene count ranged from 410 to 434 (std dev 5.53).
- Pairwise Jaccard minimum was 0.8826, meaning even the two most dissimilar runs shared 88.3% of their union.
- 10 of 29 curated core genes are not in the iron-clad set.
- 5 iron-clad genes are not in the Tier 2 blind single-run core gene set, indicating they appear in >=90% of stability runs but were absent from the specific single run used for Tier 2.
- 35 clusters contain >= 5 iron-clad members each.

---

*Draft generated by Kai (Claude Code agent), 2026-04-29. Pending review by Cody Sigmon.*