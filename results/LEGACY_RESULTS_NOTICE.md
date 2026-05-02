# Legacy Results Notice

The `/results/` directory contains **pre-SOP and exploratory outputs** produced
during engine development and early validation (v0.3.0 through v0.3.2). These
runs predate the Disease-Day SOP (v1.0, established 2026-04-24) and do not
include the standardized provenance, three-tier validation structure, or
stability profiling required by the current methodology.

## Where to find current validation evidence

**SOP-compliant disease-day results** are in `/disease_days/`:

- `disease_days/2026-04-28_asd_v033/` — ASD v0.3.3, three-tier SOP-compliant
- `disease_days/2026-04-29_ipf_v0332/` — IPF v0.3.3.2, three-tier SOP-compliant

Each SOP-compliant disease day includes a `DISEASE_DAY_MANIFEST.md` with full
provenance (git commit, engine version, seed file checksums, run timings).

## What the legacy results contain

Legacy results may include curated and blind single-run outputs for various
diseases (ASD, T2D, IBD, AD, Breast Cancer, IPF, Psoriasis, CRC, negative
control). These were produced under earlier engine versions and configurations.
Numbers may differ from current SOP-compliant runs due to methodology fixes
applied in v0.3.3 (e.g., mean log2FC dilution fix, random seed wiring,
discovery_tissues plumbing).

## Evidence hierarchy

For any claim or citation:

1. **Use `disease_days/` outputs** when available (SOP-compliant, current engine)
2. **Use `results/` outputs** only for historical comparison or when no disease-day
   run exists yet for that disease
3. **Label the source version** whenever citing numbers from either location

The `/output/` directory (if present) contains intermediate pre-SOP runs from
April 12-14, 2026 and should not be cited.
