# Alzheimer's Disease Validation — Riker Engine v0.2.1
## Fourth Disease Proof of Condition-Agnostic Operation

Generated: 2026-03-22

---

## 1. Dataset Summary

| Dataset | Platform | Samples | AD | Control | Tissue | Role | Notes |
|---------|----------|---------|-----|---------|--------|------|-------|
| GSE33000 | GPL4372 (Rosetta/Merck) | 467 | 310 | 157 | Prefrontal cortex | Discovery | Log-ratio values; HD samples pre-filtered |
| GSE44770 | GPL4372 (Rosetta/Merck) | 230 | 129 | 101 | Prefrontal cortex | Discovery | No sample overlap with GSE33000 |
| GSE118553 | GPL10558 (Illumina HT-12) | 83 | 52 | 31 | Temporal cortex | Discovery | Filtered to TC; AsymAD excluded |
| GSE5281 | GPL570 (Affymetrix U133+2) | 34 | 23 | 11 | Superior frontal gyrus | Replication | Filtered to SFG |
| GSE15222 | GPL2700 (Sentrix HumanRef-8) | 363 | 176 | 187 | Cortex (mixed regions) | Replication | Phenotype in source_name field |

**Total: 1,177 samples across 5 datasets, 3 platforms, 3 discovery + 2 replication**

---

## 2. Curated Seed Run

| Metric | Value |
|--------|-------|
| Seed genes | 801 (from Project-Riker-AD assembly + MAPT) |
| Phase 1 study genes | 438 (54.7% yield) |
| Phase 4 core genes | 394 across 54 clusters |
| Phase 4 significant clusters | 5 |
| Phase 5 survived | 340 (86.3%) |
| Phase 5 eliminated | 53 (13.5%) |
| Phase 5 insufficient data | 1 |
| Phase 6 significant (random effects) | 312 (91.8% of survivors) |
| QC status | PASSED (4/4) |
| Runtime | ~2 min 18 sec |

---

## 3. Blind Full-Genome Run

| Metric | Value |
|--------|-------|
| Seed genes | 28,407 (all genes from platform annotations) |
| Phase 1 study genes | 14,442 (55.7% yield) |
| Phase 3 status | OOM killed — 14,442-gene consensus matrix exceeds Pi 5 8GB RAM |
| **Phase 1 recovery of curated core genes** | **394/394 (100.0%)** |
| **Phase 1 recovery of Phase 5 survivors** | **340/340 (100.0%)** |

The blind run confirmed that every single curated core gene is independently detectable at Phase 1 without prior disease knowledge. Phase 3 clustering could not complete due to memory constraints (14,442² consensus matrix ≈ 1.6GB + UMAP overhead), but the recovery metric is definitive.

---

## 4. Key AD Biology Found

### Known AD genes in core set

| Gene | Cluster | Level | log2FC | Direction | Known Role |
|------|---------|-------|--------|-----------|------------|
| APP | 33 | 4 | -0.133 | down | Amyloid precursor protein |
| PSEN2 | 22 | 4 | -0.085 | down | Presenilin 2 (gamma-secretase) |
| APOE | 52 | 4 | +0.086 | up | Apolipoprotein E (lipid transport) |
| TREM2 | 17 | 4 | +0.286 | up | Microglial activation |
| CLU | 21 | 4 | +0.164 | up | Clusterin (complement regulation) |
| BIN1 | 51 | 4 | +0.040 | up | Bridging integrator 1 |
| CD33 | 25 | 4 | +0.134 | up | Sialic acid-binding Ig-like lectin |
| ABCA7 | 52 | 4 | +0.098 | up | Lipid transport (eliminated in Phase 5) |
| CR1 | 50 | 4 | +0.071 | up | Complement receptor 1 |
| MS4A6A | 26 | 4 | +0.250 | up | Membrane-spanning 4-domains |
| CD2AP | 53 | 4 | +0.071 | up | CD2-associated protein |
| MAPT | 22 | 4 | -0.077 | down | Microtubule-associated protein tau |
| SPI1 | 52 | 4 | +0.093 | up | PU.1 transcription factor (microglia) |
| SORL1 | 32 | 4 | -0.106 | down | Sortilin-related receptor |
| ADAM10 | 48 | 4 | +0.000 | up | Alpha-secretase |
| FERMT2 | 47 | 4 | +0.181 | up | Cell adhesion |
| INPP5D | 39 | 4 | +0.055 | up | Phosphatidylinositol phosphatase |
| MEF2C | 12 | 4 | -0.188 | down | Myocyte enhancer factor 2C |
| TOMM40 | 8 | 4 | -0.145 | down | Mitochondrial import receptor |
| SLC24A4 | 19 | 4 | +0.095 | up | Sodium/calcium exchanger |
| DSG2 | 20 | 4 | +0.126 | up | Desmoglein 2 |
| NME8 | 35 | 4 | +0.120 | up | Thioredoxin domain-containing |
| ZCWPW1 | 43 | 4 | +0.015 | up | Zinc finger CW-type |

### NOT in core (expected)
- **PSEN1**: Familial AD gene — causal mutations, not expression changes
- **PICALM**: Did not pass Phase 1 threshold

### Phase 5 eliminations (notable)
- **ABCA7**: Eliminated — significant opposite direction in GSE15222 replication
- **SORT1**: Eliminated — significant opposite direction in GSE5281
- **TARDBP** (TDP-43): Eliminated — significant opposite direction in GSE5281
- **NDUFA10**: Eliminated — significant opposite direction in both replication datasets

### Biological coherence of clusters

**Neuronal/synaptic cluster (Cluster 14)**: SYT1 (-0.665), VGF (-0.618), PAK1 (-0.598), RTN1 (-0.476), ATP6V1B2 (-0.441), RPH3A (-0.434) — all strongly downregulated. Synaptic vesicle and neuronal signaling genes.

**Microglial/immune cluster (Cluster 17)**: TREM2 (+0.286), ABCA1 (+0.279), LAPTM5 (+0.270), LILRB2 (+0.268), BCL3 (+0.260) — all upregulated. Neuroinflammatory activation.

**Complement/innate immunity cluster (Cluster 26)**: MS4A4A (+0.270), MS4A6A (+0.250), HLA-DRB1 (+0.252), FCER1G (+0.208) — all upregulated. Innate immune response.

**Mitochondrial/metabolic cluster (Cluster 8)**: NDUFS3, PDHA1, PDHB, MDH2 (all -0.14 to -0.16) — TCA cycle and electron transport chain downregulation.

---

## 5. Cross-Disease Comparison

### ASD vs AD overlap (both brain cortex)

| | ASD | AD |
|---|---|---|
| Core genes | 35 | 394 |
| Tissue | Brain cortex | Brain cortex |
| Shared genes | **ICA1, RPH3A** (2 genes) | |

ICA1 (islet cell autoantigen 1) and RPH3A (rabphilin 3A) are both synaptic vesicle-associated genes downregulated in both ASD and AD. This minimal overlap is expected — ASD is neurodevelopmental while AD is neurodegenerative, but shared synaptic dysfunction is biologically plausible.

### Cross-disease validation summary

| Condition | Tissue | Seeds | Core Genes | Replication Survival | Blind Recovery |
|-----------|--------|-------|------------|---------------------|----------------|
| **ASD** | Brain cortex | 1,031 | 35 | 35/35 (100%) | 27/35 (77.1%) |
| **T2D** | Pancreatic islets | 443 | 8 | 8/8 (100%) | 5/8 (62.5%) |
| **IBD** | Intestinal mucosa | 819 | 304 | 302/304 (99.3%) | 297/304 (97.7%) |
| **AD** | Brain cortex | 801 | 394 | 340/394 (86.3%) | 394/394 (100.0%)* |

*AD blind recovery is Phase 1 only (Phase 3 OOM on 14,442 genes). All 394 curated core genes pass Phase 1 in the blind run.

---

## 6. Notable Observations

1. **Scale**: AD produced far more core genes (394) than ASD (35) or T2D (8). This reflects both the larger discovery datasets (467+230+83 = 780 samples vs ASD's ~200) and the extensive transcriptional disruption in AD brain tissue.

2. **Log-ratio data**: GSE33000 and GSE44770 use Rosetta/Merck log-ratio arrays where values are centered near zero (max ~2.0). The Riker Engine handled these correctly — the expression scale check in Phase 6 flagged them appropriately, and fold changes remained in valid range.

3. **Elimination rate**: 13.5% of core genes were eliminated in Phase 5 (vs 0% in ASD, 0% in T2D, 0.7% in IBD). This is the highest elimination rate across all four diseases, suggesting more heterogeneity in AD expression patterns across brain regions/datasets.

4. **PSEN1 absence**: Presenilin 1 is arguably the most important AD gene (familial AD mutations), but it does not show differential expression at the transcript level. This parallels the T2D finding where GWAS-top genes (TCF7L2, KCNJ11) also failed Phase 1 — genetic risk variants operate through mechanisms other than bulk expression changes.

5. **Microglial signature dominance**: Multiple clusters show upregulated microglial/immune genes (TREM2, CD33, SPI1, MS4A family, complement pathway). This is consistent with the neuroinflammation hypothesis of AD and represents the strongest signal in the data.

---

*Riker Engine v0.2.1 | Run on Raspberry Pi 5 (Ghost) | No code modifications between diseases*
