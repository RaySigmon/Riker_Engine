# Stability Profile Regeneration

Per-run outputs (`runs/run_001/` through `runs/run_050/`) are not committed to git.
They are regenerable from the committed config and master seed.

## To regenerate

```bash
cd /home/kai001/riker-engine
python3 scripts/stability_profiling.py \
  disease_days/2026-04-29_ipf_v0332/blind_pc/config.yaml \
  -n 50 \
  --master-seed 42 \
  --keep-runs \
  --output-dir disease_days/2026-04-29_ipf_v0332/stability_50run
```

## Requirements
- Engine version: v0.3.3.2 (commit 939663e)
- Master seed: 42
- Config: `blind_pc/config.yaml` (committed)
- Expected wall clock: ~9.2 hours on Pi 5

## Committed summary files
- `stability_summary.json` — aggregate statistics (iron-clad count, Jaccard, appearance distribution)
- `stability_report.csv` — per-gene appearance frequency and classification
- `stability_pairwise_jaccard.csv` — all 1,225 pairwise Jaccard values
- `stability_report.txt` — human-readable summary
- `run_summary.csv` — per-run core gene counts and timings
- `iron_clad_cluster_analysis.csv` — cluster provenance for iron-clad genes
- `iron_clad_cluster_summary.csv` — cluster-level aggregation
- `iron_clad_aggregated.csv` — cross-run Phase 6 effect size aggregation
