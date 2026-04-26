# Riker Engine SOP — Lessons Log

Running log of corrections and refinements to DISEASE_DAY_SOP.md discovered during execution.
Each entry documents: what the SOP said, what actually happened, what was changed, and why.

This log preserves institutional knowledge about why the SOP evolved. When future maintainers
ask "why does the SOP say X," they can trace back here.

---

## Entry 001 — 2026-04-24 — First SOP-compliant execution (ASD disease-day)

### Lesson 1.1: CLI flags in SOP did not exist

**SOP said:** `riker run <config> --output-dir <path> --random-seed <int>`

**Reality:** `riker run` only takes a config file path and optional `-v`. No `--output-dir`
or `--random-seed` flags exist in the current CLI.

**Correction:** SOP Standard Invocation rewritten. Output directory and random seed are now
set inside a per-run disease-day YAML config (copied from the base config for that disease,
modified for `output_dir`, `random_seed`, `phase3.seeds`, and `phase4.seed`). Each run has its
own config file preserved in its disease-day directory.

**Why this is better:** Every run's exact configuration is preserved as a file, not just in
shell history. More reproducible. Matches how the stability profiler already works.

### Lesson 1.2: All-expressed file portability contradiction

**SOP said two things that contradicted:**
- Pre-run checklist: "All required configs use relative paths (not absolute `/home/kai001/`)"
- Optional Component: "Do not regenerate all-expressed files"

**Reality:** The archived all-expressed files exist only at absolute paths in `riker-archive/`.
Following both rules made the all-expressed run impossible.

**Correction:** All-expressed files are copied (not regenerated) into the repo at
`data/seeds/<disease>_all_expressed.csv`. Each copy preserves the original file's content
verbatim and includes a header comment documenting the source archive path, date captured,
gene count, and originating procedure.

**Why this is better:** Preservation-by-copying is not the same as regeneration. The files
are now version-controlled, reproducible for anyone who clones the repo, and survive archive
reorganization. The "do not regenerate" rule was about not recomputing from raw probe data
(which introduces variance). Copying verbatim preserves the archived snapshot without that risk.

### Lesson 1.3: RAM threshold too strict

**SOP said:** "Ghost has ≥4GB available RAM"

**Reality:** Pi 5 typically shows 3.5-3.8GB available after OS/background services. ASD blind
runs have completed successfully at these levels.

**Correction:** SOP threshold lowered to ≥3GB with a note that this is empirically validated
on Pi 5.

### Lesson 1.4: RIKER_ENGINE_METHODS.md line 288 error

**Methods doc said:** "`--random-seed` CLI flag is the actual reproducibility control"

**Reality:** No `--random-seed` CLI flag exists. The YAML `random_seed` field is parsed but
not consumed by any phase. The actual stochastic control is `phase3.seeds` (UMAP seeds) and
`phase4.seed` (permutation seed), both set in the YAML config.

**Correction deferred:** Methods doc will be updated when v0.3.3 wires `config.random_seed`
to derive phase3/phase4 seeds. For now, the SOP correctly specifies setting all three fields
explicitly in per-run configs.

---

*Each subsequent disease-day should add entries here if anything surprising happens.*
