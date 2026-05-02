# Riker Engine — Claims Policy

This file governs how the project describes its capabilities and results.
All documentation, README text, papers, presentations, and agent-facing
summaries must use language consistent with the tier assigned here.

---

## Claim tiers

### Tier 1 — Validated
Claims supported by completed SOP-compliant disease-day runs with
reproducibility receipts (manifest, config, pipeline outputs, stability
profile). Reviewers can audit the evidence chain.

### Tier 2 — Supported
Claims supported by internal validation but not yet independently
reproduced or SOP-compliant. May rely on pre-SOP single runs.

### Tier 3 — Experimental / hypothesis-generating
Features or results that generate hypotheses for further investigation.
Must be labeled as such in any context.

---

## Current claim assignments

| Claim | Tier | Evidence |
|-------|------|----------|
| Cross-dataset candidate-gene filtering | 1 — Validated | 8 diseases, SOP-compliant runs for ASD + IPF |
| Directional replication with tissue-aware elimination | 1 — Validated | Phase 5 implemented, verified in audit |
| Full-seed-set FDR correction | 1 — Validated | Stress test: 420 false positives prevented |
| Consensus clustering robustness | 1 — Validated | 15-config sweep, stability profiling |
| REML meta-analysis with heterogeneity stats | 1 — Validated | Code-audited, output verified |
| Iron-clad gene stability profiling (50-run) | 1 — Validated | ASD (394, J=0.933), IPF (2309, J=0.955) |
| Per-disease validation counts (curated single-run) | 2 — Supported | 8 diseases run under v0.3.2, 2 re-validated under v0.3.3+ |
| WGCNA benchmark comparison | 2 — Supported | Single-disease (ASD) comparison, not multi-disease |
| Cold replication (held-out dataset) | 2 — Supported | IPF GSE47460: 86.3% replication |
| Blind all-gene mode | 3 — Experimental | Works empirically but labeled hypothesis-generating |
| Novel candidate discovery | 3 — Experimental | "Novel" = survived filters from full protein-coding seed set |
| Drug target recovery | 2 — Supported | Known targets recovered, not new targets proposed |
| Clinical subtype separation | 3 — Experimental | Breast cancer ER/HER2 modules observed, not clinically validated |

---

## Prohibited wording (without specific qualification)

These words require explicit evidence backing when used in any project output:

- **"proves"** — the engine provides statistical evidence, not proof
- **"clinical-grade"** — no clinical validation has been performed
- **"diagnostic"** — the engine is a research tool, not a diagnostic
- **"causal"** — the engine identifies expression associations, not causal mechanisms
- **"drug target"** (unqualified) — distinguish "known target recovered" from "new target proposed"
- **"completely blind"** — blind mode uses all protein-coding genes as seed set, which is a defined universe, not truly unconstrained
- **"clinical classification from scratch"** — observed expression modules that correspond to known subtypes, not independent clinical classification
- **"zero author involvement"** — requires signed protocol or agent transcript as evidence

---

## Revision history

- v1.0 (2026-05-01): Initial policy established during project cleanup session.
