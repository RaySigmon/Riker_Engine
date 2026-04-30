# Session Handoff — 2026-04-30

## What is the Riker Engine

A condition-agnostic transcriptomics pipeline that identifies reproducible
gene modules across independent datasets. Six-phase pipeline: cross-referencing,
feature construction, consensus clustering, robustness testing, replication,
meta-analysis. Built by Ray Sigmon / Alpha Research Labs.

## Current engine version

**v0.3.3.2** — tagged and pushed to GitHub on 2026-04-29.

The engine is **frozen** for the 8-disease re-validation. No code changes
unless an emergency bug is found during validation. If something breaks,
stop and report — do not fix in place.

## What happened in this session (April 28-30, 2026)

### Engine work (completed)
- **v0.3.3 release** (8 rounds of fixes):
  - Round 1: tissue labels, dead random_seed, broken test, RNG migration
  - Round 2: geo_parser row drop, p-value floor, discovery_tissues plumbing, Phase 5 skip path
  - Round 3: phase2/phase3 CSV output standardization
  - Round 4: mean log2FC dilution fix (Option C — significant-only averaging)
  - Round 5: code_version in pipeline output, --version CLI flag
  - Round 6: documentation sweep (METHODS v1.1, SOP v1.0.2)
  - Round 7: stability profiler pairwise Jaccard + stability_summary.json
  - Round 8: config portability (path resolution relative to config file)
- **v0.3.3.2 hotfix**: symbol_column wiring from YAML, tissue field required per dataset
- **Profiler hotfix**: path resolution for temp configs written to /tmp/

### Validation (completed)
- **ASD disease day** — all 3 tiers complete, disease report written
- **IPF disease day** — all 3 tiers complete, disease report written

### Documentation (completed)
- ENGINE_WALKTHROUGH.md — full structural review
- CLAUDE.md — operational instructions for Kai sessions
- Updated SOP, METHODS, CONFIGURATION docs

## What needs to happen next

### Immediate: 6 remaining disease days

The 8-disease re-validation is 2/8 complete. Remaining diseases in
suggested order:

1. **IBD** — needs blind config created (same pattern as IPF blind)
2. **AD (Alzheimer's)** — needs blind config created
3. **Breast Cancer** — needs blind config created
4. **CRC** — blind config exists (`crc_blind.yaml`)
5. **Psoriasis** — blind config exists (`psoriasis_blind.yaml`)
6. **T2D** — needs blind config created

For each disease, the protocol is identical (SOP v1.0.2):
1. Copy curated config to `disease_days/YYYY-MM-DD_<disease>_v0332/curated/config.yaml`
2. Adjust paths (../../ → ../../../ for 3-level depth), set output_dir: .
3. Run Tier 1 (curated), report raw numbers, wait for approval
4. Copy blind config similarly, run Tier 2, report, wait for approval
5. Launch Tier 3 (50-run stability, ~7-13 hours depending on disease)
6. After Tier 3 completes: run analysis scripts, build disease report
7. Commit

### Key scripts for disease reports
```bash
# Cluster analysis
python3 scripts/analyze_iron_clad_clusters.py <stability_dir>

# Cross-run effect size aggregation
python3 scripts/aggregate_phase6_iron_clad.py <stability_dir>

# Disease report (currently manual — build_disease_report.py is ASD-specific)
```

### After all 8 diseases complete
- Clean up the computer (old test data scattered across multiple locations)
- Clean up the GitHub repo
- Update README with v0.3.3.2 validation numbers
- Engine paper preparation

## File locations that matter

### Current (v0.3.3.2)
- **Engine source**: `~/riker-engine/riker/`
- **Tests**: `~/riker-engine/tests/` (319 tests)
- **Configs**: `~/riker-engine/configs/examples/`
- **Disease day outputs**: `~/riker-engine/disease_days/`
- **Scripts**: `~/riker-engine/scripts/`
- **Docs**: `~/riker-engine/docs/`

### Historical (pre-v0.3.3, DO NOT USE for current work)
- `~/riker-engine/output/` — pre-SOP runs from April 12-14
- `~/riker-engine/results/` — committed v0.3.2 results
- `~/riker-archive/old-iterations/` — earliest runs from March
- `~/riker-engine/stability_ASD_blind/` — v0.3.2 ASD stability (376 iron-clad)
- `~/riker-engine/disease_days/2026-04-24_asd/` — v0.3.2 SOP run
- `~/riker-engine/disease_days/2026-04-28_asd_v033_verify/` — engineering verification artifacts
- `~/riker-engine/disease_days/2026-04-28_asd_v033_r4_verify/` — same
- `~/riker-engine/disease_days/2026-04-28_asd_v033_r8_verify/` — same

### Config path convention
Configs in `configs/examples/` are 2 levels deep from repo root. When
copying to `disease_days/YYYY-MM-DD_<disease>/curated/` (3 levels deep),
all `../../` paths must become `../../../`. Set `output_dir: .` for
disease day runs.

## Diseases with blind configs ready
- ASD: `configs/examples/asd_blind.yaml` ✓
- IPF: `configs/examples/ipf_blind.yaml` ✓
- CRC: `configs/examples/crc_blind.yaml` ✓
- Psoriasis: `configs/examples/psoriasis_blind.yaml` ✓

## Diseases needing blind configs created
- IBD: create from `ibd_bulk.yaml`, swap seed to `all_protein_coding_genes.csv`
- AD: create from `ad_bulk.yaml`, same swap
- Breast Cancer: create from `breast_cancer_bulk.yaml`, same swap
- T2D: create from `t2d_bulk.yaml`, same swap

Pattern: copy curated config, change condition name to `<DISEASE>_blind`,
change seed_genes to `../../data/seeds/all_protein_coding_genes.csv`,
change output_dir, add explicit phase3/phase4 blocks.

## Completed disease day results

### ASD (2026-04-28, engine v0.3.3)
- Curated: 29 core, 29 survived, 10 meta-significant
- Blind: 421 core, 414 survived, 233 meta-significant
- Stability: 394 iron-clad, Jaccard median 0.933
- Report: `disease_days/2026-04-28_asd_v033/DISEASE_REPORT_DRAFT.md`

### IPF (2026-04-29, engine v0.3.3.2)
- Curated: 186 core, 166 survived, 155 meta-significant
- Blind: 2451 core, 1856 survived, 1584 meta-significant
- Stability: 2309 iron-clad, Jaccard median 0.955
- Report: `disease_days/2026-04-29_ipf_v0332/DISEASE_REPORT_DRAFT.md`

## Standing rules (from this session)
- Report raw numbers only during validation. No interpretation.
- Flag judgment calls before executing — team discusses first.
- Architectural claims require empirical verification.
- Engine frozen during validation — stop and report bugs, don't fix.

## Unpushed commits
Several commits are local but not pushed to GitHub. Next session should
push after Cody provides a PAT, or Cody can push manually.

## v0.3.4 roadmap
Saved in chat from April 30. Key items: unused import cleanup, plots.py
buildout, WGCNA benchmark, sensitivity analysis figure, hardcoded
thresholds to config, stability profiler refactor.
