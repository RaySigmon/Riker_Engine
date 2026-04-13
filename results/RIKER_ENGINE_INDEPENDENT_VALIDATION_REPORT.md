# Riker Engine — Independent Cold Replication & Novel Disease Validation Report
## April 12, 2026 | Ray Sigmon / Alpha Research Labs

---

## Executive Summary

The Riker Engine (v0.3.1) was subjected to two rounds of fully independent cold-start testing using an AI agent (Gemini CLI, designated "Jim") with zero prior knowledge of the engine, its methodology, or its expected results. The agent was given only the GitHub URL and a skeptical PI-style prompt asking for replication and generalization testing.

**Result:** The agent independently confirmed the engine as "reproducible science" and recommended adoption for lab use — without prompting, coaching, or access to expected outcomes.

---

## Test Environment

- **Hardware:** Raspberry Pi 5 (8GB RAM), clean filesystem with no prior Riker Engine artifacts
- **Agent:** Google Gemini CLI (codename "Jim") — chosen deliberately as a weaker instruction-follower than Claude to simulate a less-careful human user
- **Isolation protocol:** All prior clones, output directories, config files, and virtual environments were deleted before each test. Agent was instructed to clone from GitHub only.
- **Repository state:** Commit `e6bf09e` (latest at time of test), including all fixes from the day's two fix rounds (8 total commits)
- **Verification:** Agent confirmed correct commit hash before testing began

---

## Test 1: Direct Replication (IPF — Idiopathic Pulmonary Fibrosis)

### Protocol

The agent was given this prompt:

> "You are a computational biology postdoc in my lab. A colleague sent me a gene clustering pipeline called 'Riker Engine' and claims it produces strong results for identifying disease-relevant gene networks. I'm skeptical but willing to look. [...] I need you to attempt a replication. [...] I want to know if this is reproducible science or a black box."

The agent was NOT told which disease to replicate, which datasets to use, or what results to expect.

### Agent's Actions

1. Cloned repository from GitHub (fresh clone, verified commit `e6bf09e`)
2. Read README.md, identified IPF as the most rigorous validation
3. Located `configs/examples/ipf_curated.yaml` (found immediately — Fix 8 working)
4. Ran `python scripts/download_data.py ipf` — all 6 datasets downloaded including held-out GSE47460 (Fix 7 working), all platform annotations downloaded without 404 errors (Fix 1 working)
5. Created virtual environment, installed with `pip install ".[clustering]"`
6. Ran pipeline: `riker run configs/examples/ipf_curated.yaml`
7. Compared output against published results in `results/ipf/curated/pipeline_summary.json`
8. Ran cold replication script against held-out dataset GSE47460

### Results

| Metric | Published | Agent's Replication | Difference |
|--------|-----------|-------------------|------------|
| Phase 1 Study Genes | 241 | 241 | 0.0% |
| Phase 4 Core Genes | 190 | 189 | -0.5% |
| Phase 5 Survived Replication | 170 (89.5%) | 169 (89.4%) | -0.6% |
| Phase 6 Meta-Significant | 157 | 156 | -0.6% |
| Core Gene Overlap | — | 189/190 (99.47%) | — |
| Missing Gene | — | ID1 | — |
| Cold Replication (GSE47460) | 86.3% (132/153) | 86.2% (131/152) | -0.1% |
| Concordance Rate | 96.7% | 96.7% | 0.0% |
| QC Report | 4 passed, 0 critical | 4 passed, 0 critical | — |

### Agent's Uncoached Assessment

> "The Riker Engine is reproducible science, not a black box. The authors' claims for the IPF validation hold up under independent replication with negligible variance."

### Consistency Across Runs

This replication was performed twice (once pre-fix, once post-fix). Both runs produced identical results: 189/190 core genes, ID1 as the only variance, 86.2% cold replication. The consistency across independent runs confirms the engine's stability.

---

## Test 2: Novel Disease Generalization (Psoriasis)

### Protocol

After the successful replication, the agent was given this follow-up prompt:

> "Good work. I'm cautiously impressed but I need one more test. Pick a disease we haven't looked at — something in your wheelhouse — and run the engine against it from scratch. Choose your own GEO datasets, build your own config, and tell me if the results make biological sense. I want to know if this tool generalizes or if it only works on their cherry-picked examples."

The agent was also pointed to `docs/NEW_DISEASE_GUIDE.md`.

### Agent's Actions

1. Read `docs/NEW_DISEASE_GUIDE.md` for configuration guidance
2. Selected **Psoriasis** as the test disease (agent's independent choice)
3. Curated a seed gene list of ~96 genes from literature knowledge (IL-17/IL-23 axis, S100 proteins, keratins, antimicrobial peptides, inflammatory mediators)
4. Selected 5 GEO datasets, all GPL570 platform:
   - **Discovery:** GSE13355 (58 LS vs 58 NL), GSE30999 (85 LS vs 85 NL), GSE41662 (28 LS vs 28 NL)
   - **Replication:** GSE50790 (17 LS vs 17 NL), GSE54456 (20 LS vs 20 NL)
5. Inspected GEO metadata to determine correct phenotype fields and case/control labels
6. Caught that two initially selected datasets (GSE63979, GSE121212) were RNA-seq, not microarray — replaced them without modifying source code
7. Built YAML config from scratch following the guide
8. Ran full pipeline: `riker run configs/psoriasis_bulk.yaml`
9. Evaluated biological relevance of results

### Key Observations During Setup

- Agent read and followed NEW_DISEASE_GUIDE.md (Fix 11 confirmed effective)
- Agent did NOT modify any source code (Fix 9 — non-editable install — confirmed effective)
- Agent kept all datasets on the same platform (GPL570) per guide recommendations
- Agent correctly used the metadata inspection commands documented in the guide
- Agent independently caught RNA-seq vs microarray mismatch and resolved it

### Results

| Metric | Value |
|--------|-------|
| Seed Genes | 96 |
| Phase 1 Study Genes | 60 |
| Phase 3 Clusters | 9 |
| Phase 4 Core Genes | 50 |
| Phase 5 Survived Replication | 50 (100% in GSE50790) |
| Phase 6 Meta-Significant | 28 (random effects) |
| QC Report | 4 passed, 0 warnings, 0 critical |

### Top Replicated Genes (Biological Validation)

The engine identified genes with well-established roles in psoriasis pathology:

**Antimicrobial Peptides / S100 Proteins (hallmark of psoriasis):**
- S100A12 (p ≈ 6×10⁻³¹⁰) — calgranulin C, neutrophil marker
- S100A9 — calprotectin subunit
- S100A7 — psoriasin (named for its discovery in psoriatic skin)
- DEFB4B — beta-defensin 4B
- PI3 — Elafin, serine protease inhibitor

**Inflammatory Mediators (IL-17/IL-23 axis):**
- IL36G, IL36RN — IL-36 pathway (key psoriasis cytokines)
- CCL20 — Th17 cell chemoattractant
- CXCL8, CXCL1 — neutrophil chemoattractants
- TNF, IL1B, IL12B — pro-inflammatory cytokines
- IL19, IL20, IL24 — IL-20 subfamily (keratinocyte proliferation)

**Keratinocyte Hyperproliferation Markers:**
- KRT16, KRT17, KRT6A — keratins upregulated in hyperproliferative epidermis
- IVL — involucrin (cornification marker)
- TGM1, TGM3 — transglutaminases
- SERPINB3, SERPINB4 — serine protease inhibitors

**Apoptosis Dysregulation:**
- BCL2 (p ≈ 7×10⁻⁵⁰, downregulated) — consistent with altered apoptosis in psoriatic plaques

### Agent's Uncoached Assessment

> "The Riker Engine is a powerful tool for meta-analytical gene discovery. It is stable, generalizable, and mathematically rigorous. I recommend we adopt it for our lab's transcriptomic screening."

> "The docs/NEW_DISEASE_GUIDE.md is accurate. I was able to go from raw GEO accessions to a finished meta-analysis in under 30 minutes."

---

## Comparison: Previous Cold Test vs. Current Cold Test

The same novel disease test was attempted in an earlier session (pre-fixes). The comparison demonstrates the impact of the documentation and usability improvements:

| Aspect | Pre-Fix (MS attempt) | Post-Fix (Psoriasis) |
|--------|---------------------|---------------------|
| Disease Guide Available | No | Yes (NEW_DISEASE_GUIDE.md) |
| Source Code Modified | Yes — rewrote ProbeGeneMapper, config.py, cli.py | No modifications |
| Platform Consistency | Mixed (GPL570, GPL96, GPL13158, GPL10558) | Single platform (GPL570) |
| Pipeline Completion | HALTED at Phase 5 | COMPLETE — all phases passed |
| Core Genes Identified | 62 (from broken methodology) | 50 (clean methodology) |
| Biological Validation | Partial (some relevant genes, no replication) | Full (textbook psoriasis biology, replicated) |
| Time Spent on Metadata | ~80% of session | ~20% of session |
| Agent's Verdict | "High-integrity tool" (generous given the struggles) | "Recommend adoption for our lab" |

---

## Cumulative Validation Summary

Including independent tests, the Riker Engine has now been validated on **8 diseases**:

| Disease | Core Genes | Blind Recovery | Replication Rate | Validated By |
|---------|-----------|---------------|-----------------|-------------|
| ASD | 35 | 77.1% | 57.1% (15 eliminated, brain-specific) | Author |
| T2D | 8 | 62.5% | 100% | Author |
| IBD | 304 | 97.7% | 99.3% | Author + Independent Agent |
| Alzheimer's Disease | 394 | 98.2% | 86.3% | Author |
| Breast Cancer | 152 | 99.3% | 100% | Author |
| IPF | 190 | — | 89.5% (86.3% cold replication) | Author + Independent Agent |
| **Psoriasis** | **50** | **94.0%** | **100%** | **Independent Agent (Gemini)** |
| **CRC** | **264** | **97.7%** | **92.8% (245/264)** | **Independent Agent (Claude)** |

Psoriasis and CRC are the first diseases validated entirely by independent parties with no author involvement in dataset selection, seed gene curation, or configuration.

---

## Files Preserved

All output files from the independent validation are preserved at:

**IPF Replication:**
- `results/ipf/independent_replication_summary.json`

**Psoriasis (Novel Disease):**
- `results/psoriasis/independent_validation/pipeline_summary.json`
- `results/psoriasis/independent_validation/phase6_meta_analysis.csv`
- `results/psoriasis/independent_validation/phase4_core_genes.csv`
- `results/psoriasis/independent_validation/qc_report.json`
- `results/psoriasis/independent_validation/psoriasis_bulk.yaml` (agent's config)

---

## Methodology Note

The testing agent (Google Gemini CLI) was deliberately chosen as a less instruction-compliant model than Claude to simulate a careless or inexperienced human user. Despite this:

- The agent independently identified the correct disease for replication (IPF)
- The agent independently selected a biologically appropriate novel disease (Psoriasis)
- The agent independently curated seed genes from literature knowledge
- The agent caught and corrected its own mistakes (RNA-seq datasets) without modifying source code
- The agent arrived at a positive verdict without coaching or expected-result framing

A follow-up test using Claude Code CLI is planned to provide cross-model verification.

---

*Report prepared by Ray Sigmon / Alpha Research Labs*
*Repository: https://github.com/RaySigmon/Riker_Engine*
*License: AGPL-3.0*
