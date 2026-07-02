# RunPod Runbook — Finish AD (Alzheimer's) Blind Disease Day

Completes the deferred AD blind protein-coding tier on cloud compute so AD joins the
other 7 diseases as a clean 8th. Drafted 2026-06-29 (verified from source). Full
analysis: task output `wxr7uqd3o`.

## ⚠️ Engine-version pin (corrected)
Run on the **canonical engine on `origin/main`** — the `riker/` source tree is
**byte-identical from commit `11fb692` through HEAD** (`git diff 11fb692 HEAD -- riker/`
is empty), and that's the engine the other diseases ran on. So clone GitHub
`origin/main` and use its HEAD.
- **Do NOT** `git checkout v0.3.3.2` (tag = `939663e`, predates `11fb692`; older docstrings + missing the `phase6_significant_random_effects` summary key).
- **Do NOT** use local commit `be2d4f3` (that's the unpushed `wip` artifacts branch — not on the remote).
- The REML meta computation is identical across all these commits, so AD will be numerically consistent with the set regardless; pin for clean provenance.
- Cross-arch note: prior runs were ARM (Pi); RunPod is x86_64. Phases 1/2/5/6 deterministic; 3/4 stochastic — the ≥90%-of-50 iron-clad rule is designed to absorb minor BLAS/numba clustering differences.

## Why this exists
AD blind OOM'd on the 8 GB Pi at **Phase 3 consensus clustering** (UMAP+HDBSCAN over 13,181 study genes; killed >6.4 GB RSS). Phases 1–2 completed; AD **curated** completed fine — only **blind** OOM'd.

## Pod
- **CPU pod** (clustering is CPU-only; GPU not needed). Highest single-core clock available; 4–8 vCPU (UMAP/HDBSCAN single-threaded when `random_state` fixed — extra cores don't help a single run).
- **RAM ≥ 32 GB** (peak ~8–12 GB; 16 GB hard minimum).
- **Disk ≥ 20 GB** (repo+data ~0.4 GB; 50-run `--keep-runs` adds a few GB).

## Environment
```bash
cd /workspace && git clone <origin/main> riker-engine && cd riker-engine
git rev-parse HEAD            # confirm a commit >= 11fb692 on main
python3 -m venv .venv && source .venv/bin/activate   # Python >=3.11
pip install -U pip && pip install -r requirements-lock.txt && pip install -e ".[clustering]"
riker --version              # 0.3.3.2 (expected)
```
Pinned (must match the 8-disease validation env): `numpy==1.26.4` (NOT 2.x — ABI-crashes hdbscan), `scipy==1.17.1`, `pandas==2.3.3`, `scikit-learn==1.8.0`, `umap-learn==0.5.11`, `hdbscan==0.8.41`, `PyYAML==6.0.2`.

## Data to sync (~289 MB; same repo-relative layout → no config edits)
- `data/geo/ad/`: GSE33000_AD_only (~83M), GSE44770 (~42M), GSE118553_TC (~30M), GSE5281_SFG (~7M), GSE15222 (~16M)
- `data/platforms/`: GPL4372.annot (~12M), GPL10558.annot (~28M), GPL570.annot (~38M), GPL2700.annot (~17M)
- `data/hgnc/hgnc_complete_set.txt` (~16M)
- `data/seeds/all_protein_coding_genes.csv` (~140K, 19,296 genes, col `symbol`)
```bash
rsync -avz kai:/home/kai001/riker-engine/data/ /workspace/riker-engine/data/
riker validate disease_days/2026-05-05_ad_v0332/blind_pc/config.yaml   # sanity
```

## Tier 2 — single blind run (the part that OOM'd)
```bash
riker run disease_days/2026-05-05_ad_v0332/blind_pc/config.yaml \
  2>&1 | tee disease_days/2026-05-05_ad_v0332/blind_pc/run.log
```
Watch `free -g`. Phase 3 ("15 UMAP+HDBSCAN configs") must now clear past HDBSCAN. ~1–2 h.

## Tier 3 — 50-run stability + iron-clad aggregation (where the iron-clad count comes from)
```bash
nohup python3 scripts/stability_profiling.py \
  disease_days/2026-05-05_ad_v0332/blind_pc/config.yaml \
  -n 50 --master-seed 42 --keep-runs \
  --output-dir disease_days/2026-05-05_ad_v0332/stability_50run \
  > disease_days/2026-05-05_ad_v0332/stability_50run/profiler.log 2>&1 &
python3 scripts/aggregate_phase6_iron_clad.py disease_days/2026-05-05_ad_v0332/stability_50run
```
Sequential; ~24–48 h wall (AD ~13.2k genes vs CRC ~8.3k). Keep `--master-seed 42` (matches the set). Don't shard unless seed derivation `master_seed + run_number` is preserved.

## Success check
```bash
grep "PIPELINE COMPLETE" disease_days/2026-05-05_ad_v0332/blind_pc/run.log
grep -i "Successful runs: 50/50" disease_days/2026-05-05_ad_v0332/stability_50run/profiler.log
grep -i "Iron-clad" disease_days/2026-05-05_ad_v0332/stability_50run/profiler.log
```
Confirm every `runs/run_*/pipeline_summary.json` + the blind summary carry the on-main `code_version` and `qc_status: PASSED`.

## Provenance + pull back (do NOT commit/push from the pod)
```bash
riker --version > .../blind_pc/environment.cloud.txt; uname -a >> .../blind_pc/environment.cloud.txt
pip freeze > .../blind_pc/pip_freeze.cloud.txt; git rev-parse HEAD > .../blind_pc/git_commit.cloud.txt
```
`rsync` `blind_pc/` + `stability_50run/` back to the Pi; update `DISEASE_DAY_MANIFEST.md` Tier 2/3; let the normal Pi repo flow handle versioning.
