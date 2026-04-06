# RIKER ENGINE — FINAL ACTION PLAN
## For Claude Code CLI (Kai) — Sequential Execution
## Created: April 6, 2026

**Context:** The Riker Engine repo (github.com/RaySigmon/Riker_Engine) was reviewed cold from a PI perspective. The methodology and code structure passed. The repo's packaging, reproducibility, and evidence did not. This plan fixes everything needed before collaborator outreach begins. Work through these in order — later phases depend on earlier ones.

**Canonical repo:** `/home/kai001/riker-engine/`
**All changes go here.** No other directory is the source of truth.

---

## RULES FOR KAI (READ BEFORE STARTING)

1. **DO NOT modify any pipeline algorithm code** (anything in `riker/phases/`, `riker/stats/`, `riker/ingestion/`, `riker/qc/`, `riker/io/`) unless explicitly required by Task 3.7 (Phase 5 tissue logic). The pipeline is frozen.
2. **DO NOT delete any files without archiving first.** Move to `riker-archive/`, never `rm`.
3. **Run `pytest` as a gate before moving to the next phase.** All tests must still pass after every phase of work.
4. **If you discover that validation results for any disease truly don't exist anywhere on this machine, STOP and report immediately.** We need to know if a re-run is required before proceeding.
5. **If any task is ambiguous, ask.** Don't guess on scientific decisions.
6. **Commit frequently** with descriptive messages. A reviewer will read your git log.

---

## PHASE 0: SECURITY (DO THIS FIRST, BEFORE ANYTHING ELSE)

### 0.1 — Rotate exposed GitHub PATs
The `.git/config` files in these directories have Personal Access Tokens embedded in remote URLs:
- `/home/kai001/riker_project/.git/config`
- `/home/kai001/Project-Riker-AD/.git/config`

**Action:**
1. Log into GitHub → Settings → Developer Settings → Personal Access Tokens
2. Identify and **revoke** any token that matches what's in those config files
3. Generate a new token if needed for the canonical repo
4. Update remotes to use `https://github.com/...` (no embedded token) or switch to SSH
5. Confirm the canonical repo at `/home/kai001/riker-engine/` does NOT have a PAT in its `.git/config`

**This is a live security issue. Do not proceed to other tasks until confirmed resolved.**

---

## PHASE 1: HARVEST BEFORE ARCHIVE

Extract valuable artifacts from old directories BEFORE moving them to archive. Everything harvested goes into or informs the canonical repo at `/home/kai001/riker-engine/`.

### 1.1 — Harvest GEO filenames for download script
In `/home/kai001/riker_project/data/raw/` (~400 MB of GEO data):
- [ ] List every GEO series matrix filename and accession number (e.g., `GSE28521_series_matrix.txt.gz`)
- [ ] List any SFARI gene list files, HGNC database files, and pathway database files
- [ ] Record the exact filenames and directory structure
- [ ] Map these to the disease validations they were used for (ASD primarily)
- [ ] **Do NOT delete these files yet** — we need the filenames to build the download script in Phase 3

Also check these locations for additional data files used in non-ASD validations:
- Any YAML configs in `/home/kai001/.riker/runs/` (they reference data paths)
- `/home/kai001/Project_Riker_ Alzheimer/` — has GWAS catalog zip and gene lists
- `/home/kai001/RIKERENGINENOTES/` — may have validation configs

**Output:** A single manifest file listing every data dependency per disease: GEO accession numbers, URLs where known, seed gene sources, pathway databases, and the config YAML that referenced them.

### 1.2 — Harvest AD validation results
- [ ] Compare contents of `/home/kai001/Project-Riker-AD/` with `/home/kai001/riker_phase1/` through `riker_phase6/`
- [ ] Confirm the phase directories are duplicates (or identify unique content in each)
- [ ] Pull the most complete AD results into `riker-engine/results/ad/`
- [ ] These are the "missing AD results" the PI review flagged — the README claims them but the repo doesn't have them

### 1.3 — Harvest breast cancer results
- [ ] Check all old directories for breast cancer validation outputs
- [ ] Check `/home/kai001/.riker/runs/` — the 11 run directories may contain these
- [ ] If they exist anywhere on disk, pull into `riker-engine/results/breast_cancer/`
- [ ] **If they don't exist on disk, STOP and report** — they'll need to be re-run

### 1.4 — Harvest T2D and IBD results
- [ ] The repo has empty subdirs for these (`results/t2d_test_run_1/`, `results/ibd_test_run_1/`)
- [ ] Check `/home/kai001/.riker/runs/` — the 11 run directories may contain these outputs
- [ ] Check `/home/kai001/RIKERENGINENOTES/User RUn/` for run CSVs
- [ ] Check `/home/kai001/ibd_validation/` for IBD results (contains a zipped run result)
- [ ] Pull complete results into `riker-engine/results/t2d/` and `riker-engine/results/ibd/`

### 1.5 — Harvest useful documentation
- [ ] Scan `RIKERENGINENOTES/`, `Project Riker/`, `RE_Foundation/`, `Project_Riker_ Alzheimer/`
- [ ] Identify anything that adds value beyond what's already in `riker-engine/docs/`
- [ ] Pull in useful content; discard redundant copies

### 1.6 — Harvest config/data references
- [ ] From all old directories, collect any YAML configs, seed gene CSVs, pathway database files
- [ ] These inform what the download script needs to fetch
- [ ] Collect GEO accession numbers used for ALL 5 diseases across all directories

**For each disease (ASD, T2D, IBD, AD, breast cancer), confirm you have:**
- [ ] Phase 4 core gene list (the locked pre-replication list)
- [ ] Phase 5 replication results
- [ ] Phase 6 meta-analysis output
- [ ] Pipeline summary JSON
- [ ] QC report

---

## PHASE 2: FILESYSTEM CLEANUP

Only after Phase 1 harvesting is complete.

### 2.1 — Create archive structure
```
mkdir -p /home/kai001/riker-archive/old-iterations
mkdir -p /home/kai001/riker-archive/planning-docs
mkdir -p /home/kai001/riker-archive/downloads-deduped
```

### 2.2 — Archive old code iterations
Move (not delete) to `riker-archive/old-iterations/`:
- [ ] `riker-engine-v020-build-archive/` — frozen v0.2.0 snapshot
- [ ] `riker_project/` — original ASD research (1.5 GB, has raw GEO data)
- [ ] `Project-Riker-AD/` — AD standalone repo
- [ ] `riker_phase1/` through `riker_phase6/` — AD manual phase run artifacts
- [ ] `ibd_validation/` — after harvesting results

### 2.3 — Archive planning docs
Move to `riker-archive/planning-docs/`:
- [ ] `Project Riker/` (action plan, outreach blueprint, architecture docs)
- [ ] `Project_Riker_ Alzheimer/` (phase directives, GWAS catalog)
- [ ] `RE_Foundation/` (context transfer, master blueprint)
- [ ] `RIKERENGINENOTES/` + `RIKERENGINENOTES.md`

### 2.4 — Deduplicate Downloads
- [ ] Keep `Riker_Engine_Dissertation-1-1.pdf` (newest), remove older copies
- [ ] Keep one copy of each unique document (`Project_Riker_Final`, `PROJECT_RIKER_PROGRESS.md`, `Project_Riker_Preprint_Draft.docx`)
- [ ] Move surviving Riker files to `riker-archive/downloads-deduped/`
- [ ] **DO NOT TOUCH** `RIKER-CARE-PROTOCOL.pdf` (personal, about Riker the person)

### 2.5 — Clean up misc
- [ ] Archive `/home/kai001/Documents/start riker engine` (tiny text snippet)
- [ ] Leave `/home/kai001/.riker/` in place (runtime data, needed by engine)
- [ ] Leave `/home/kai001/.openclaw/workspace/memory/Riker/` in place (personal)

### 2.6 — Verify
- [ ] Confirm `riker-engine/` is unchanged and `pytest` still passes
- [ ] Confirm no symlinks or references from the canonical repo point to archived directories
- [ ] Home directory should now be clean: `riker-engine/`, `riker-archive/`, `.riker/`, and non-project files

---

## PHASE 3: REPO FIXES (Engineering)

All work in `/home/kai001/riker-engine/`. These are the issues that would make a PI's postdoc unable to clone and verify.

### 3.1 — Fix version mismatch
- [ ] `riker/__init__.py` says `__version__ = "0.1.0"` → change to `0.3.0`
- [ ] `pyproject.toml` says `version = "0.3.0"` → confirm correct
- [ ] Search for any other files referencing a version number and align them
- [ ] Two-minute fix. Do it first.

### 3.2 — Commit all 5 disease validation results

Using the outputs harvested in Phase 1, create this structure in the repo:

```
results/
├── asd/
│   ├── config.yaml           ← exact config used for this run
│   ├── pipeline_summary.json
│   ├── core_genes.csv
│   ├── meta_analysis.csv
│   └── qc_report.json
├── t2d/
│   └── (same structure)
├── ibd/
│   └── (same structure)
├── ad/
│   └── (same structure)
└── breast_cancer/
    └── (same structure)
```

- [ ] ASD: already has results committed — verify completeness against this structure
- [ ] T2D: commit harvested results
- [ ] IBD: commit harvested results
- [ ] AD: commit harvested results
- [ ] Breast cancer: commit harvested results (or flag if re-run needed)

**Critical:** Include the YAML config for each disease. A reviewer needs to see exactly what inputs produced these outputs. If configs reference local paths, update them to use relative paths or placeholder paths with comments explaining where to get the data.

### 3.3 — Build data download script

**This is the single most important task in this entire document.** Without it, nobody can use the tool. The README currently says pathway databases "must be provided manually. Automatic download planned." That stops a postdoc dead.

Using the manifest from Task 1.1, build a download system:

**Option A: Makefile targets**
```
make download-asd    # Downloads ASD GEO data + SFARI genes
make download-t2d    # Downloads T2D GEO data + Open Targets genes
make download-all    # Everything
```

**Option B: CLI subcommand**
```
riker download --condition asd --data-dir ./data/
```

**What the script must do per disease:**
1. Download all GEO series matrix files (use `wget` or `urllib` from NCBI FTP: `ftp://ftp.ncbi.nlm.nih.gov/geo/series/...`)
2. Download seed gene list from source (SFARI, Open Targets, etc.) or bundle a frozen copy in `data/seeds/` (with version/date noted)
3. Download pathway databases (KEGG, Reactome, MSigDB Hallmark) or document how to obtain them (MSigDB requires registration — note this clearly)
4. Download HGNC gene mapping data
5. Place everything in the correct directory structure matching the config YAML
6. Verify downloads (check file sizes, MD5 if available)

**Target end state:** `git clone && ./scripts/download_data.sh asd && riker run configs/asd.yaml` works on a fresh machine.

### 3.4 — Add GitHub Actions CI

Create `.github/workflows/test.yml`:
```yaml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.11", "3.12"]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}
      - run: pip install -e ".[dev]"
      - run: pytest --tb=short -q
```

- [ ] Update the README badge to link to actual Actions workflow status, not a static badge
- [ ] Consider adding a minimal integration test: small synthetic dataset through all 6 phases

### 3.5 — Fix test_qc.py
- [ ] Currently contains only `pass` — the QC module has zero test coverage
- [ ] Write actual tests for `qc/checks.py`: test that QC catches bad inputs, flags low-sample datasets, validates config structure, halts execution on failures
- [ ] If the QC module itself is a stub, either implement it or **remove test_qc.py entirely** and update the test count
- [ ] **Do not ship a test file that is `pass`.** It actively damages credibility.

### 3.6 — Add .gitignore and clean artifacts
The repo has no `.gitignore`. Build artifacts are getting committed.

- [ ] Create a proper Python `.gitignore` (`__pycache__/`, `*.egg-info/`, `.pytest_cache/`, `dist/`, `build/`, `*.pyc`)
- [ ] Clean committed build artifacts:
```bash
git rm -r --cached riker_engine.egg-info/
git rm -r --cached __pycache__/
# etc.
```
- [ ] Move the 16 `PHASE*.md` files from root into `docs/development/` or remove them. Root should contain: README, LICENSE, pyproject.toml, Makefile, .gitignore, and source directories. That's it.

### 3.7 — Generalize Phase 5 tissue logic

**Problem:** Phase 5 directional replication elimination is hardcoded for brain/blood tissue pairs only. For T2D (islets), IBD (mucosa), and breast cancer (tumor), Phase 5 is effectively a pass-through. This means 3 of 5 validated diseases don't actually experience the "six-phase progressive filtering" as described.

**Fix (implement both):**
- **(a) Generalize tissue matching** — Accept a `tissue_type` field in the config YAML. Compare replication datasets only against same-tissue discovery datasets. Any same-tissue opposite-direction significant effect triggers elimination, regardless of tissue type.
- **(b) Document the limitation** — In docs/PIPELINE.md: "Phase 5 elimination applies when same-tissue replication datasets are available. For conditions where all datasets are same-tissue (e.g., IBD mucosal biopsies), all replication comparisons are same-tissue and elimination applies normally. For conditions with only one tissue type and no held-out replication datasets, Phase 5 reports concordance statistics without elimination."

**After the fix:** Re-examine whether T2D, IBD, and breast cancer results would change. Document the answer in each disease's results directory.

### 3.8 — Pin dependencies
- [ ] `pyproject.toml` uses `>=` with no upper bounds — dangerous for reproducibility
- [ ] Add a `requirements-lock.txt` or use `poetry.lock` / `uv.lock`
- [ ] Especially critical: pin `umap-learn` — API changes across minor versions could break Phase 3
- [ ] Pin `hdbscan`, `scikit-learn`, `scipy` as well

### 3.9 — Fix README claims
- [ ] Update test badge to reflect actual passing test count after Task 3.5
- [ ] Fix install instructions: `Project-Riker.git` → `Riker_Engine.git`
- [ ] Fix citation section: references should point to `Riker_Engine`, not `Project-Riker`
- [ ] Verify license badge says AGPL-3.0 (check pyproject.toml license field — previously caused sidebar to show MIT)

### 3.10 — Dockerfile (lower priority)
If time permits after all above tasks:
```dockerfile
FROM python:3.11-slim
COPY . /app
WORKDIR /app
RUN pip install -e .
ENTRYPOINT ["riker"]
```
Not a blocker for v0.3.1 but removes one more "I couldn't install it" excuse.

---

## PHASE 4: SCIENTIFIC ADDITIONS

These address the PI review's methodological concerns. Some require running new analyses.

### 4.1 — Negative control experiment

**Why:** The engine has only been tested on diseases with known strong signals. A skeptic asks: "Does it return nothing when there's nothing to find?" Without a negative control, you can't answer.

**Run at least one of these:**
1. **Random gene set:** Take 500 randomly selected genes (not associated with any disease), run them through the pipeline against the ASD datasets. The engine should return zero or near-zero core genes. If it returns a large module, that's a problem.
2. **Shuffled labels:** Take the real ASD config but randomize the case/control labels. The engine should find no consistent signal.
3. **Unrelated tissue:** Run ASD seed genes against non-brain datasets (e.g., the T2D islet data). The engine should find minimal signal.

**Commit the results** to `results/negative_control/` with a README explaining what was tested and what was found.

### 4.2 — Systematic drug target hit rate

**Problem:** The dissertation and README cherry-pick individual drug targets (ABCC8, JAK2, ERBB2, etc.) without reporting systematic hit rates.

**For each of the 5 diseases:**
1. Compile a list of ALL known drug targets for that disease (use Open Targets, DrugBank, or similar — document which source)
2. Count how many are in the core gene set
3. Count how many passed Phase 1 but didn't make core
4. Count how many weren't expressed / didn't pass Phase 1
5. Report: `X of Y known drug targets found in core genes (Z%)`

**This turns anecdote into evidence.** Add a table to the README or a file in `results/drug_target_analysis/`.

### 4.3 — Commit blind run outputs
- [ ] Commit the full raw outputs from all completed blind (hypothesis-free) runs
- [ ] These substantiate the "100% Phase 1 blind recovery" and blind recovery rate claims
- [ ] Without them, the claims are unverifiable
- [ ] Place in `results/<disease>/blind_run/`

### 4.4 — AD blind run (requires 16+ GB RAM)
- [ ] The dissertation notes it exceeded 8 GB RAM at Phase 3 (14,442 study genes)
- [ ] Either: run on a machine with 16+ GB RAM, or document the memory ceiling clearly
- [ ] This is the only disease where the blind run didn't complete
- [ ] **If no 16 GB machine is available, this is not a blocker** — document it as a known limitation and move on

### 4.5 — Breast cancer blind run
- [ ] Listed as "N/R" (not run) in the validation table
- [ ] Run it and commit results, or document why it wasn't feasible
- [ ] If it requires >8 GB, same approach as AD — document the limitation

### 4.6 — WGCNA comparison reframing (docs only)
The PI review correctly noted the WGCNA comparison is a workflow comparison, not a method comparison. The dissertation partially acknowledges this but the framing still reads like a head-to-head.

**In `docs/PIPELINE.md` and the README, add a clarifying statement:**
"The comparison demonstrates the workflow advantage of integrated multi-dataset analysis over manual single-dataset analysis followed by cross-study comparison. It is not a claim that the Riker Engine's clustering method is superior to WGCNA's network construction; the tools answer different questions at different scales."

---

## PHASE 5: PUBLICATION READINESS

After Phases 0–4 are complete.

### 5.1 — Preprint
- [ ] The dissertation is essentially a preprint draft — needs formatting for submission
- [ ] Add: negative control results, systematic drug target analysis, complete blind run data
- [ ] bioRxiv rejected 3x for affiliation requirements — alternatives: arXiv q-bio (needs endorser), Zenodo (DOI but no peer review), Research Square, OSF Preprints
- [ ] This is the single most important credibility signal for PI outreach

### 5.2 — JOSS submission (September 2026 target)
- [ ] Requires 6 months of public development history (clock started when repo went public)
- [ ] Most requirements will exist after Phase 3: CI, documentation, community guidelines
- [ ] JOSS scrutinizes AI-assisted development — public issues, PRs, and iterative development history over the next 6 months strengthens the submission

### 5.3 — External cold-start reproduction
- [ ] Write a step-by-step guide: from `git clone` to reproduced ASD results
- [ ] Have someone else (NOT you, NOT Kai) follow it and report where they get stuck
- [ ] This is worth more than any amount of README polish
- [ ] Their feedback becomes your final cleanup pass before outreach

---

## PHASE 5: FINAL VERIFICATION

After all above phases are complete, before pushing:

### 5.1 — End-to-end test
1. Clone the repo fresh into `/tmp/riker-test/`
2. Run the download script for ASD
3. Run the ASD pipeline from the committed config
4. Verify outputs match committed results
5. Run `pytest` — all tests must pass
6. Delete `/tmp/riker-test/`

### 5.2 — Push
1. `git add -A && git status` — review everything being committed
2. Commit: `"v0.3.1: reproducibility overhaul — results, data scripts, CI, cleanup"`
3. `git tag v0.3.1`
4. `git push origin main --tags`
5. Verify GitHub Actions CI runs and passes
6. Verify the repo looks clean from a browser

---

## PRIORITY ORDER IF TIME-CONSTRAINED

If you can't do everything, this is the order that matters most:

1. **Phase 0** — Security / rotate PATs (10 minutes)
2. **Task 3.1** — Version mismatch (2 minutes)
3. **Task 3.2** — Commit all validation results (1–2 hours depending on harvest)
4. **Task 3.3** — Data download script (half day — the single biggest impact item)
5. **Task 3.5** — Fix test_qc.py (30 minutes)
6. **Task 3.6** — .gitignore + clean root artifacts (15 minutes)
7. **Task 3.9** — Fix README claims (15 minutes)
8. **Task 3.4** — GitHub Actions CI (30 minutes)
9. **Task 4.1** — Negative control experiment (1–2 hours for the run)
10. **Task 3.7** — Phase 5 tissue logic generalization (1–2 hours)
11. **Task 3.8** — Pin dependencies (30 minutes)
12. **Task 4.2** — Drug target hit rate analysis (1–2 hours)
13. **Task 4.6** — WGCNA reframing (15 minutes)
14. Everything else

---

## EXECUTION TIMELINE

```
Phase 0  ──→  Phase 1  ──→  Phase 2  ──→  Phase 3  ──→  Phase 4  ──→  Phase 5
(security)   (harvest)    (cleanup)     (repo fixes)   (science)     (publish)
 10 min       2-3 hrs      1 hr          multi-day      multi-day     weeks
```

Phases 0–2 are filesystem work (can do today).
Phase 3 is engineering work (this week).
Phase 4 requires running analyses (depends on compute time).
Phase 5 is writing/publishing (ongoing).

---

## FILES REFERENCE

| Item | Location |
|------|----------|
| Canonical repo | `/home/kai001/riker-engine/` |
| Runtime data | `/home/kai001/.riker/` |
| Dissertation (latest) | `/home/kai001/Downloads/Riker_Engine_Dissertation-1-1.pdf` |
| This plan | `/home/kai001/RIKER_ENGINE_ACTION_PLAN.md` |
| Archive (after Phase 2) | `/home/kai001/riker-archive/` |
