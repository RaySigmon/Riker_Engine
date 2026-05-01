# Riker Engine — Project Briefing for AI Agents

## What this is

The Riker Engine is a condition-agnostic transcriptomics pipeline that
identifies reproducible gene modules across multiple independent expression
datasets. It enforces pre-specification, full-seed-set FDR correction,
consensus clustering, and directional replication — returning a minimal core
gene set that survives all filters.

The engine operates as a six-phase progressive filtering pipeline. The phases
are locked. The methodology is what's being validated, with eight diseases as
test cases.

## Current phase (2026-05)

The engine is FROZEN at v0.3.3.2 for an 8-disease re-validation under a unified
SOP (Standard Operating Procedure for disease-day runs). No pipeline code
changes during validation. Two diseases are complete (ASD, IPF). Six remain.
After all eight complete, results consolidate into a methods paper
plus per-disease application papers.

## How this project is organized

- **Owner:** Cody Sigmon (publishing as Ray Sigmon / Alpha Research Labs)
- **Implementation:** Kai — Claude Code CLI on Ghost (Raspberry Pi 5)
- **QA / architect / audit:** Claude.ai conversations
- **External cold reviewers:** rotated, document-restricted to fight bias

## Operating norms (rules every AI agent should follow)

- **Don't infer when uncertain — ask Cody for files.** GitHub serves stale
  cache for this repo. The mirror at rikerengine.quickaffordablesites.com is the
  freshest publicly-readable surface, but Cody's local repo on Ghost is
  ground truth.
- **The six pipeline phases are LOCKED and never modified.** Any phase code
  change during the 8-disease validation is a regression, not an improvement.
- **Kai instructions must be clean technical directives only.** No
  interpretive commentary inside Kai instruction blocks.
- **"Everything checks out" should trigger more verification, not less.**
  Reassuring conclusions are when audits fail.
- **Always clean up stale artifacts when spotted.** Don't leave half-finished
  things floating.

## ASD result terminology (precision matters here)

- **376** = iron-clad genes across 50-run blind stability profile (>=90% appearance, v0.3.2-era)
- **394** = iron-clad genes from v0.3.3 ASD disease day, 50-run blind stability
- **403** = core genes from run_001 of the v0.3.3 ASD blind stability profile (single run only)
- **35** = core genes from the *curated* ASD run with SFARI seeds (different config entirely)

These are distinct numbers meaning distinct things. Do not conflate.

## Where things live in this mirror

- Engine code: `/raw/riker/`
- Pipeline phases (LOCKED): `/raw/riker/phases/`
- Statistical primitives: `/raw/riker/stats/`
- Documentation: `/raw/docs/`
- Disease day outputs: `/raw/disease_days/`
- Tests: `/raw/tests/`
- SOP: `/raw/docs/DISEASE_DAY_SOP.md`
- Methods paper: `/raw/docs/RIKER_ENGINE_METHODS.md`
- Engine walkthrough: `/raw/docs/ENGINE_WALKTHROUGH.md`

## Data files not served

Raw GEO datasets, platform annotation files, and HGNC databases are referenced
in configs but not served by this mirror. They are obtained directly from
GEO/HGNC per the recipe in `docs/DATA_RECONSTRUCTION.md`.

## Active blockers

See `/STATE.json` field `disease_days_completed` and the disease list there for
current status. As of last update:

- IBD disease day blocked by missing GEO data and platform file GPL13158.annot.
- All other pending diseases (CRC, Psoriasis, Breast Cancer, AD, T2D)
  have data and configs ready.

## Contact

Cody Sigmon — see /STATE.json for current project info.
