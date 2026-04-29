# CLAUDE.md — Operational Instructions for Riker Engine

This file is read by Claude Code (Kai) at the start of each session.

## Project overview

Riker Engine is a condition-agnostic transcriptomics pipeline for discovering
replicated gene modules across independent datasets. Repository owner: Ray
Sigmon / Alpha Research Labs.

## Key locations

- Engine source: `riker/`
- Tests: `tests/` (318 tests as of v0.3.3)
- Configs: `configs/examples/` (paths are relative to config file location)
- Disease day outputs: `disease_days/YYYY-MM-DD_<disease>/`
- Scripts: `scripts/` (stability_profiling.py, aggregate_phase6_iron_clad.py, etc.)
- Documentation: `docs/` (SOP, METHODS, CONFIGURATION, PIPELINE)

## Running the engine

```bash
riker run <config.yaml>        # Full pipeline
riker validate <config.yaml>   # Config check only
riker --version                # Print version
python -m pytest tests/ -q     # Run test suite
```

## Verification discipline

When reporting on engine behavior, distinguish between:

- **Empirical findings:** "the data shows X" (backed by run outputs)
- **Architectural claims:** "the engine is designed to do X" (backed by
  code inspection or reasoning)

Architectural claims require empirical verification before being treated
as conclusions. When in doubt, verify.

**Standing rule:** When offering an architectural explanation for unexpected
results, verify empirically before presenting it as a conclusion. Plausible
explanations that sound correct can still be wrong. Run the data check first,
then report both the claim and the verification.

## Reporting discipline

- Report raw numbers. Do not characterize results as "matches expected,"
  "looks reasonable," or "deviates from baseline."
- Flag judgment calls before executing them so the team can discuss.
- When comparing across runs or versions, report overlap counts and gene
  lists. Let Cody interpret.

## Engine determinism profile

| Phase | Deterministic? | Stochastic element |
|---|---|---|
| Phase 1 (differential expression) | Yes | None |
| Phase 2 (feature matrix) | Yes | None |
| Phase 3 (consensus clustering) | **No** | UMAP random_state |
| Phase 4 (robustness testing) | **No** | Permutation seed |
| Phase 5 (replication) | Yes | None |
| Phase 6 (meta-analysis) | Yes | None |

Verified empirically April 29, 2026: Phase 6 values identical to 17-18
significant digits across sampled stability runs. Code inspection confirmed
zero stochastic operations in Phases 1, 5, and 6.

## Config path resolution

Paths in YAML configs are resolved relative to the config file's parent
directory. Absolute paths pass through unchanged. The sentinel value "auto"
for hgnc_path is preserved. When copying configs to deeper directories
(e.g., disease_days), adjust relative path depth accordingly.

## Current version

v0.3.3 (tagged April 28, 2026). Engine is frozen for 8-disease re-validation.
No code changes during validation. If unexpected behavior is found, stop and
report — do not fix in place.
