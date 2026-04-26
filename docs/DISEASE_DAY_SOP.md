# Riker Engine — Disease-Day Standard Operating Procedure

**Version:** 1.0
**Effective date:** April 24, 2026
**Owner:** Ray Sigmon / Alpha Research Labs
**Applies to:** All disease validation runs after this date
**Status:** Authoritative — all 8-disease validations from this date forward follow this procedure

---

## Purpose

Ensure every disease validation is performed consistently, produces comparable outputs, and yields results that can be clearly and accurately described in papers. Standardizes the protocol so that reviewers, collaborators, and future papers can trust that each of the 8 diseases was treated identically.

This SOP is binding for all disease runs going forward. Deviations must be documented explicitly in the run's output directory with a `DEVIATION.md` note stating what was changed and why.

---

## Scope

A "disease-day" is a structured validation effort for one disease, producing two or three independent result sets plus a stability profile, depending on what seed-file artifacts are available.

**Mandatory components:**

1. **Curated run** — single run using a biologically-motivated seed list (e.g., SFARI for ASD, Open Targets for T2D, disease-specific curated sets)
2. **Protein-coding blind run** — single run using `data/seeds/all_protein_coding_genes.csv` (19,296 curated HGNC protein-coding genes — the canonical going-forward standard)
3. **50-run stability profile** — 50 runs of the protein-coding blind config with deterministically varied seeds, producing per-gene stability classifications

**Optional component:**

4. **All-expressed blind run** — single run using the preserved archive file at `data/seeds/<disease>_all_expressed.csv` in the repo. These files are copied (not regenerated) from `/home/kai001/riker-archive/old-iterations/<disease>_validation/data/seed_genes/all_expressed_genes.csv` with a header comment documenting source, date captured, and originating procedure. If neither the repo copy nor the archived file exists for a given disease, the all-expressed run is SKIPPED. Do not regenerate all-expressed files — the historical procedure used raw platform-specific probe identifiers that cannot be cleanly reproduced. Preservation-by-copy is allowed and encouraged.

All outputs are preserved in a single `disease_days/YYYY-MM-DD_<disease>/` directory structure with full provenance metadata.

---

## Why multiple single runs plus one stability profile?

Each single run answers a different scientific question:

| Run type | Question answered | Contribution to the paper |
|---|---|---|
| Curated | Within known disease-associated biology, what is convergent across cohorts? | Validates field knowledge |
| Protein-coding blind | What does the standardized protein-coding gene space yield without hypothesis? | Canonical discovery result |
| All-expressed blind (optional) | For historical-comparison only: what did the broader platform-measured identifier set yield in earlier methodology? | Methodology evolution note |

The 50-run stability profile is done only on the protein-coding blind config. Stability profiling characterizes reproducibility under seed variation — a property of the statistical pipeline, not of the seed set. Running it once per disease on the canonical seed set is sufficient. Reviewers can see seed-set effect from the single-run comparisons; they can see stochastic reproducibility from the 50-run profile.

---

## Cohort selection pre-specification

Cohort inclusion is documented per-disease following these pre-specified criteria:

- **Minimum 3 discovery cohorts** of the disease-relevant tissue type
- **Minimum 1 replication cohort**, tissue-matched where possible
- **Sample size ≥ 15 case + 15 control per cohort** (targeting ≥ 30/30 where available)
- **Platform**: Any microarray or RNA-seq platform with sufficient annotation coverage
- **Quality**: Peer-reviewed publication associated with the cohort; clear case/control definitions
- **Exclusion rationale** for any cohort meeting the above but not used must be documented (e.g., "GSE####### excluded due to unclear disease stratification")

Cohort selection is pre-specified per disease-day. Cohorts are not added or removed after results are seen. If a cohort choice is reconsidered post-hoc, the entire disease-day is redone with the revised cohort set and the prior result is archived with an explanation.

Where a disease has many cohorts available (IPF, IBD, breast cancer, CRC), selection favors tissue-matched, larger cohorts but is not required to include all available data. Where a disease has limited cohorts (ASD brain tissue), all available tissue-relevant cohorts are used and documented as such.

---

## Pre-run checklist

Before starting a disease-day, the executor must verify:

- [ ] Disease has curated seed file in `data/seeds/<disease>_curated_genes.csv` or equivalent
- [ ] Protein-coding seed file exists at `data/seeds/all_protein_coding_genes.csv` (shared, 19,296 genes)
- [ ] Archived all-expressed file status checked — if exists, the optional run is included; if not, the optional run is skipped (no regeneration)
- [ ] Disease's GEO data files are present in `data/geo/<disease>/`
- [ ] All required configs exist and use relative paths (not absolute `/home/kai001/` paths)
- [ ] `riker` CLI is on PATH (`which riker` returns `/home/kai001/.local/bin/riker`)
- [ ] `riker.__version__` matches the intended version (should be `0.3.2` or later)
- [ ] Ghost has ≥20GB free disk and ≥3GB available RAM (empirically validated floor on Pi 5)
- [ ] No competing Ghost workload scheduled for the disease-day window

If any checkbox fails, halt and resolve before proceeding.

---

## Standard invocation

All runs use **master seed 42** as the canonical reference seed. This is the seed that reproduces the published ASD result. Deviations from master seed 42 must be documented explicitly.

### Curated run (mandatory)

Create a per-run disease-day config by copying the base disease config and setting
`output_dir`, `random_seed`, `phase3.seeds`, and `phase4.seed` in the YAML directly:

```bash
cd /home/kai001/riker-engine
DD="disease_days/$(date +%Y-%m-%d)_<disease>"
mkdir -p "$DD/curated"
cp configs/examples/<disease>_curated.yaml "$DD/curated/config.yaml"
# Edit config: set output_dir=$DD/curated, random_seed=42, phase3.seeds, phase4.seed
riker run "$DD/curated/config.yaml" 2>&1 | tee "$DD/curated/run.log"
```

### Protein-coding blind run (mandatory)

Same per-run config pattern as curated. Base config is `configs/examples/<disease>_blind.yaml`:

```bash
cd /home/kai001/riker-engine
DD="disease_days/$(date +%Y-%m-%d)_<disease>"
mkdir -p "$DD/blind_pc"
cp configs/examples/<disease>_blind.yaml "$DD/blind_pc/config.yaml"
# Edit config: set output_dir=$DD/blind_pc, random_seed=42, phase3.seeds, phase4.seed
riker run "$DD/blind_pc/config.yaml" 2>&1 | tee "$DD/blind_pc/run.log"
```

### All-expressed blind run (optional — only if archive seed file exists)

```bash
# Verify the preserved seed file exists in the repo:
ls data/seeds/<disease>_all_expressed.csv || echo "No all-expressed file — skipping"

cd /home/kai001/riker-engine
DD="disease_days/$(date +%Y-%m-%d)_<disease>"
mkdir -p "$DD/allexpressed"
# Start from blind_pc config and swap seed file to the all-expressed version
cp "$DD/blind_pc/config.yaml" "$DD/allexpressed/config.yaml"
# Edit config: set seed_genes=data/seeds/<disease>_all_expressed.csv, output_dir=$DD/allexpressed
riker run "$DD/allexpressed/config.yaml" 2>&1 | tee "$DD/allexpressed/run.log"
```

### 50-run stability profile (mandatory, overnight)

```bash
cd /home/kai001/riker-engine
mkdir -p disease_days/$(date +%Y-%m-%d)_<disease>/stability_50run
nohup python3 scripts/stability_profiling.py \
  configs/examples/<disease>_blind.yaml \
  -n 50 \
  --master-seed 42 \
  --keep-runs \
  --output-dir disease_days/$(date +%Y-%m-%d)_<disease>/stability_50run \
  > disease_days/$(date +%Y-%m-%d)_<disease>/stability_50run/profiler.log 2>&1 &
```

Use `tmux` or `screen` when executing interactively. `nohup` alone protects against SSH disconnect; `tmux` additionally protects against terminal close.

---

## Required outputs per run

### Single runs (curated, protein-coding blind, and all-expressed if run)

Each single run must produce the full intermediate output set, preserving phase-by-phase transparency:

- `pipeline_summary.json` — top-level results, gene counts per phase
- `qc_report.json` — QC pass/fail per check
- `phase1_study_genes.csv` — genes surviving Phase 1
- `phase2_feature_matrix.csv` — feature matrix produced for clustering (if phase writes it)
- `phase3_cluster_assignments.csv` — final cluster assignment per gene
- `phase4_all_levels.csv` — all genes with per-level survivorship (pre-filter view — enables full auditability)
- `phase4_core_genes.csv` — core genes (post ≥3-per-cluster filter)
- `phase5_verdicts.csv` — per-gene replication verdict
- `phase6_meta_analysis.csv` — pooled effect sizes and meta-analysis statistics
- `run.log` — stdout/stderr from the run
- `timings.csv` — per-phase wall clock (if produced by pipeline)

### 50-run stability profile

- Per-run subdirectories (`runs/run_001/` through `runs/run_050/`), each containing at minimum `phase4_core_genes.csv` and `pipeline_summary.json`
- `runs/run_001/` additionally contains the full intermediate output set (as specified for single runs above) — provides a representative single-run view alongside the aggregate
- `stability_report.csv` — per-gene stability classification (iron-clad / borderline / stochastic)
- `stability_summary.json` — summary statistics, per-run runtime, total wall time
- `profiler.log` — profiler-level log

---

## Required provenance metadata

Each disease-day directory must contain a `DISEASE_DAY_MANIFEST.md` at its root with:

- Disease name
- Date of run
- Git commit hash at time of run (`git rev-parse HEAD`)
- Riker version (`riker --version`)
- Master seed used (should be 42 unless otherwise documented)
- Seed file SHA256 for each run (curated, protein-coding, all-expressed-if-run)
- GEO dataset accessions used with file SHA256s
- Cohort inclusion criteria applied and any exclusion rationale
- Ghost system state at start (disk free, RAM available)
- Total wall clock per run
- Executor (Kai, Cody manual, RunPod, etc.)

This manifest is the reproducibility receipt. Any future question of "how was this result produced" must be answerable from this file alone.

For runs NOT meeting these manifest standards (existing archived runs), do not modify them. Instead, create a retrospective manifest in the archive noting what IS known from the preserved artifacts.

---

## Post-run validation

After every disease-day completes, execute these validation steps:

1. **Gene count consistency check:**
   - Curated core gene count matches historical record for already-validated diseases
   - Protein-coding blind count is sensible (expect 100-500 for most diseases; higher counts warrant Phase 4 filter review)
   - All-expressed count (if run) can be compared to archived result as methodology continuity check
   - Stability iron-clad count is reported in `stability_report.csv`

2. **Flag 6 annotation** — note in the manifest:
   - `n_significant_clusters` from Phase 4 permutation test (reported, not a gate)
   - `n_clusters_contributing_core_genes` (always ≥ n_significant_clusters under current design)
   - Cluster-size filter impact (how many Level 2 survivors excluded by the ≥3 threshold)

3. **Stability convergence check** — count direction flips across the 50-run profile. Zero direction flips means signal direction is stable. Non-zero flips require investigation.

4. **Runtime capture** — per-run wall clock should be 5-10 minutes on Ghost for typical diseases. Outliers (>15 min or <3 min) warrant investigation.

5. **Cluster provenance for headline findings** — for any finding that will be described in a paper, trace it back through `phase4_core_genes.csv` to identify which clusters contributed. Confirm whether the biological grouping is a single computational cluster (rare) or a post-hoc biological annotation spanning multiple clusters (typical).

---

## Precise terminology (binding)

To prevent framing errors like the ASD manuscript's cluster characterization, the following terms have fixed meanings and must be used consistently in all papers, presentations, and documentation:

| Term | Definition | Example |
|---|---|---|
| **Iron-clad gene** | A gene appearing in ≥90% of 50 stability runs | ASD iron-clad count: 376 |
| **Borderline gene** | A gene appearing in 50-89% of 50 stability runs | |
| **Stochastic gene** | A gene appearing in <50% of 50 stability runs | |
| **Core gene (single run)** | A gene surviving Phase 4 (Level 2 + cluster-size filter) in a single pipeline run | ASD run_001 core count: 403 |
| **Core gene (curated)** | Same definition, from a curated-seed run — a different gene set than a blind run | ASD curated core count: 35 |
| **Significant cluster** | A cluster passing Bonferroni-corrected permutation test in Phase 4 — REPORTED but NOT A GATE on core gene selection | |
| **Level 2 survivor** | A gene with p<0.01 in ≥2 discovery datasets — the per-gene statistical criterion | |
| **Biological annotation** | Post-hoc functional grouping of iron-clad genes using pathway/function databases | ASD mito/OxPhos 41-gene set = biological annotation across 27 consensus clusters, NOT a single computational cluster |

Binding rule for papers: Never describe a post-hoc biological grouping as a "cluster identified by the pipeline." Describe findings as "iron-clad gene sets with biological annotation" when the genes span multiple consensus clusters.

---

## Paper-ready criteria

A disease result is paper-ready when:

1. Full disease-day has been completed (curated + protein-coding blind + 50-run stability; all-expressed if archive exists)
2. All required outputs present and QC-passed
3. `DISEASE_DAY_MANIFEST.md` is complete
4. Cluster-provenance annotation documented for the headline finding
5. Biological annotations labeled as such (not framed as computational clusters)
6. Cross-run comparison tabulated in the manifest

Disease-days passing all six criteria can be cited in papers. In-progress disease-days can be referenced but not cited as paper-ready.

---

## Deferred items (not blocking SOP compliance)

Queued for v0.3.3 / v0.4.0 releases:

- **`config.random_seed` wiring** — currently dead; `--random-seed` CLI flag is used instead
- **`discovery_tissues` plumbing** — currently not passed to Phase 5; cross-tissue tolerance inactive (conservative default)
- **Config portability** — replace absolute `/home/kai001/` paths with relative/env-var paths
- **`code_version` in pipeline_summary.json** — embed git commit hash in output for automatic provenance
- **Phase output standardization** — ensure the required intermediate output set (above) is consistently produced by all phases

These are improvements, not blockers. The SOP can be followed without them.

---

## Revision history

- **v1.0 (2026-04-24):** Initial SOP.
- **v1.0.1 (2026-04-24):** Patched during first execution. See `SOP_LESSONS.md` Entry 001.
  Corrected CLI flag error (use per-run YAML configs, not non-existent `--output-dir`/`--random-seed`).
  Clarified all-expressed preservation path (copy to `data/seeds/<disease>_all_expressed.csv`).
  Lowered RAM threshold from 4GB to 3GB. Fixed stability profiler arg syntax. Established after April 23 internal audit. All-expressed blind run made optional after determining the archived generation procedure uses raw platform-specific probe identifiers that cannot be cleanly reproduced. Added cohort selection pre-specification criteria. First SOP-compliant disease-day: ASD three-run comparison + IPF 50-run stability (depending on execution sequence).

---

*This SOP lives at `docs/DISEASE_DAY_SOP.md` in the repository. All changes are tracked through git commits. Amendments require written rationale preserved in the commit message.*
