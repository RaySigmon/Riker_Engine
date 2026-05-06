# Stability Profile Regeneration

Per-run outputs (`runs/run_001/` through `runs/run_050/`) are not committed to git.

## To regenerate

```bash
cd /home/kai001/riker-engine
python3 scripts/stability_profiling.py \
  disease_days/2026-05-05_t2d_v0332/blind_pc/config.yaml \
  -n 50 --master-seed 42 --keep-runs \
  --output-dir disease_days/2026-05-05_t2d_v0332/stability_50run
```

- Engine version: v0.3.3.2
- Master seed: 42
- Expected wall clock: ~3.9 hours on Pi 5
