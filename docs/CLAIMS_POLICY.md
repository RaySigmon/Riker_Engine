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
| Iron-clad gene stability profiling (50-run) | 1 — Validated | 7 diseases, J=0.851–0.970, see FINDINGS_SUMMARY §5 |
| Per-disease validation counts (curated single-run) | 1 — Validated | 8/8 diseases within ±2% of v0.3.2 baseline (v0.3.3.2) |
| Cross-engine-version stability (v0.3.2 → v0.3.3.2) | 1 — Validated | 7/8 within ±2%, 2 exact matches, see FINDINGS_SUMMARY §2 |
| Regime model (yield → localized/global) | 1 — Validated | 7 diseases, 0 misclassifications, see REGIME_MODEL.md |
| BrCa subtype reconstruction (blind mode) | 1 — Validated | 6 canonical programs, 15/17 markers 50/50 iron-clad |
| Headline rediscoveries iron-clad | 1 — Validated | IAPP 49/50, FAM107A 50/50, SFARI overlap 50/50 |
| WGCNA benchmark comparison | 2 — Supported | Single-disease (ASD) comparison, not multi-disease |
| Cold replication (held-out dataset) | 2 — Supported | IPF GSE47460: 86.3% replication |
| Blind all-gene mode | 1 — Validated | 7 diseases with full 3-tier validation; regime model characterizes output |
| Novel candidate discovery (blind mode) | 1 — Validated | 4 localized-regime diseases with documented rediscoveries |
| Drug target recovery | 2 — Supported | Known targets recovered, not new targets proposed |
| Clinical subtype separation | 1 — Validated | BrCa: 6 canonical programs independently recovered, stability-verified |
| Negative control false-positive characterization | 1 — Validated | 50-trial matched config; supersedes exploratory v0.3.2 result (5 core/1.0%); see FINDINGS_SUMMARY §7 |

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
- v1.1 (2026-05-08): Updated with 8-disease validation results. Promoted blind mode, subtype separation, and stability profiling to Tier 1. Added regime model, headline rediscoveries, cross-version stability, and negative control entries. Retired exploratory negative control claim (5 core/1.0%); superseded by matched 50-trial measurement.
