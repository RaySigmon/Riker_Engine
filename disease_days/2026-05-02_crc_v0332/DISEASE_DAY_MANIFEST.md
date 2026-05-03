# Disease Day Manifest — CRC

**Disease:** Colorectal Cancer (CRC)
**Date:** 2026-05-02
**Git commit at run time:** 11fb6921892018b6d116c6e034f9ef634517e1ff
**Riker version:** 0.3.3.2
**Master seed:** 42
**Executor:** Kai (Claude Code CLI agent) under Cody Sigmon direction
**SOP version:** 1.0.3

---

## Seed file checksums

| Run | Seed file | SHA256 |
|-----|-----------|--------|
| Curated | `crc_curated_genes.csv` (515 genes, 1 HGNC-remapped: H2AFX→H2AX) | `10f62164533d9f7388e2b561723e25235a44380d6b5c47abd4b560597e694197` |
| Blind | `all_protein_coding_genes.csv` (19,296 genes) | `e9ab2790c45a3501545d60298a0fc51c7cd2791b040002387665f671c98639aa` |
| Stability | same as Blind | same as Blind |

---

## GEO datasets used

| Accession | Role | Tissue | Platform | Description |
|-----------|------|--------|----------|-------------|
| GSE20916 | discovery | colon | GPL570 | Colon adenocarcinoma vs normal (~66 CRC + ~34 normal) |
| GSE9348 | discovery | colon | GPL570 | Early-stage CRC tumor vs healthy mucosa (70 + 12) |
| GSE23878 | discovery | colon | GPL570 | Colon tumour vs normal paired tissue (35 + 24) |
| GSE32323 | replication | colon | GPL570 | Paired CRC tumor/normal (17 + 17) |
| GSE37364 | replication | colon | GPL570 | CRC adenocarcinoma vs normal mucosa |
| GSE39582 | replication | colon | GPL570 | CRC molecular subtypes (566 tumors + 19 normal) |

All 6 datasets use GPL570 (Affymetrix HG-U133 Plus 2.0). All tissue-matched (colon-to-colon). No cross-tissue tolerance required.

---

## Run results

### Tier 1 — Curated CRC run
- Wall clock: ~2 min
- Phase 1 study genes: 331
- Phase 4 core genes: 262
- Phase 4 significant clusters: 4
- Phase 5 survived: 245 (17 eliminated)
- Phase 6 meta-significant (random effects): 219

### Tier 2 — Blind protein-coding run (transcriptomic characterization)
- Wall clock: ~34 min
- Phase 1 study genes: 8,288 (42.9% of protein-coding genome — global cancer rewiring)
- Phase 4 core genes: 6,052
- **Phase 4 significant clusters: 0** (no localized structure detected by permutation test)
- Phase 5 survived: 5,790 (262 eliminated, 4.3% elimination rate)
- Phase 6 meta-significant (random effects): 5,036

**Interpretation:** The blind run detected massive, genome-wide differential expression
between CRC tumor and normal colon (43% of the protein-coding genome passes Phase 1).
The Bonferroni-corrected permutation test found zero significant clusters — meaning no
cluster showed coherence distinguishable from random gene groups drawn from this highly
perturbed background. The 6,052 core genes reflect the unfiltered breadth of CRC
transcriptomic rewiring, not localized disease-specific biology.

**Contrast with IPF blind** (23% Phase 1 yield, 24 significant clusters): IPF's blind
run earned its 2,451 core genes through real statistical gates. CRC's blind run did not.

**The curated run (Tier 1) is the publishable CRC validation.** The blind Tier 2/3
results characterize the transcriptomic scale of CRC but should not be cited as a
CRC-specific gene set.

### Tier 3 — 50-run stability profile (engine reproducibility under global rewiring)
- Wall clock: 70,764 sec (19.7 hours)
- Runs completed: 50/50 (0 failures)
- Total unique genes seen: 6,192
- Iron-clad (>=90%): 5,748
- Borderline (50-89%): 393
- Stochastic (<50%): 51
- Iron-clad fraction: 92.83%
- Pairwise Jaccard median: 0.9627 (p25=0.9606, p75=0.9648)
- Core gene counts: min=5,994, max=6,057, mean=6,027.4
- Clusters with >=5 iron-clad members: 577

**Note:** The high iron-clad fraction (92.8%) reflects engine-level reproducibility
under stochastic variation, not biological specificity to CRC. The engine is stable
even when processing an unfiltered high-yield disease.

---

## Flag 6 annotation
- `n_significant_clusters` (Tier 2 blind): **0** — no localized structure
- `n_significant_clusters` (Tier 1 curated): **4** — real localized structure
- The curated run passes the permutation gate; the blind run does not
- This is the expected behavior when >40% of the genome is differentially expressed

---

## v0.3.2 comparison (curated only)
- v0.3.2 curated: 264 core genes (from `results/crc/curated/`)
- v0.3.3.2 curated: 262 core genes
- Delta: -2 genes

---

## Provenance (SOP v1.0.3)
- `exact_commands.txt` — captured at execution time
- `environment.txt` — system state at execution time
- `pip_freeze.txt` — exact package versions at execution time

## Manifest status
Generated on 2026-05-03 from pipeline outputs. All numbers verified against committed output files.
