# INSTRUCTION SET FOR KAI — PHASE 16: README UPDATE + TEST COUNT RESTORATION

## TWO TASKS

### TASK A: Update README.md

The snRNA-seq pseudo-bulking module is now built and tested. The README
currently lists it as a roadmap TODO. Update the following sections:

1. **In the "Installation" section**, add scanpy as an optional dependency:
```
**snRNA-seq extras** (`pip install scanpy`):
- scanpy — AnnData/h5ad file support for single-nucleus data
```

2. **In the "Quick Start" section**, add a subsection after the YAML
   example showing snRNA-seq usage:

```markdown
### Using snRNA-seq Data

For single-nucleus RNA-seq data, pseudo-bulk by donor before running the pipeline:

\```python
from riker.ingestion.snrnaseq import pseudo_bulk_from_counts, pseudo_bulk_from_h5ad

# From an h5ad file (AnnData format)
result = pseudo_bulk_from_h5ad(
    "data/asd_snrnaseq.h5ad",
    donor_key="donor_id",
    condition_key="disease",
    cell_type_key="cell_type",
    target_cell_type="Excitatory",
    case_values=["ASD"],
    control_values=["Control"],
)

# result.expression is a genes × samples DataFrame
# result.phenotypes is a {sample_id: 'case'/'control'} dict
# Both plug directly into Phase 1 cross-referencing

# Or from a CSV count matrix
import pandas as pd
counts = pd.read_csv("data/counts_with_metadata.csv")
result = pseudo_bulk_from_counts(
    counts,
    donor_col="donor_id",
    condition_col="condition",
    cell_type_col="cell_type",
    target_cell_type="Excitatory",
    case_values=["ASD"],
    control_values=["Control"],
)
\```

The pseudo-bulking module:
- Aggregates nuclei per donor using sum (recommended for DE analysis)
- Applies CPM + log2(CPM+1) normalization
- Filters donors with fewer than 10 nuclei (configurable)
- Filters nuclei with fewer than 200 detected genes (QC)
- Auto-detects case/control from condition labels
- Outputs a genes × samples matrix that feeds directly into Phase 1
```

3. **In the "Roadmap" section**, move snRNA-seq from TODO to DONE:
   Change `- [ ] snRNA-seq pseudo-bulking module` to
   `- [x] snRNA-seq pseudo-bulking module (riker/ingestion/snrnaseq.py)`

4. **In the "Architecture" section**, add snrnaseq.py to the ingestion tree:
```
├── ingestion/          # Data loading
│   ├── gene_db.py      # Seed genes + HGNC resolution
│   ├── normalizer.py   # Log2 detection + fold change validation
│   ├── geo_parser.py   # GEO series matrix + platform annotation
│   └── snrnaseq.py     # snRNA-seq pseudo-bulking (h5ad + CSV)
```

5. **Update the test count badge** at the top from 248 to the actual
   final count after Task B is complete.

---

### TASK B: Restore Missing Phase Tests

The Phase 14 rewrite of tests/test_phases.py lost approximately 37 tests
during file compression. The current file has ~90 tests but should have
~130 (77 original phase tests + ~12 operational shell tests + the
Phase 13 meta-analysis tests).

**DO NOT rewrite the entire file.** Instead:

1. First, count the current tests:
```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_phases.py --collect-only -q 2>&1 | tail -5
```

2. Then compare against what each phase should have by running each
   test class individually and counting:
```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_phases.py --collect-only -q 2>&1 | grep "test_" | wc -l
```

3. Check which test classes exist vs what's expected:

**Expected test classes and approximate test counts:**
- TestCrossReferenceGene: 5 tests (Phase 1 single gene)
- TestRunPhase1: 7 tests (Phase 1 full pipeline)
- TestPathwayDatabase: 5 tests (Phase 2 DB operations)
- TestFilterPathways: 4 tests (Phase 2 filtering)
- TestBuildFeatureMatrix: 6 tests (Phase 2 feature matrix)
- TestRunPhase2: 2 tests (Phase 2 integration)
- TestBuildConsensusMatrix: 5 tests (Phase 3 consensus)
- TestRunConsensusClustering: 8 tests (Phase 3 full pipeline)
- TestClusterSignificance: 3 tests (Phase 4 permutation)
- TestSensitivityAnalysis: 4 tests (Phase 4 sensitivity)
- TestLOOStability: 3 tests (Phase 4 LOO)
- TestCoreGenes: 2 tests (Phase 4 core genes)
- TestRunPhase4: 1 test (Phase 4 integration)
- TestReplicateGene: 4 tests (Phase 5 single gene)
- TestEliminationProtocol: 4 tests (Phase 5 elimination)
- TestClusterVerdicts: 2 tests (Phase 5 verdicts)
- TestRunPhase5: 3 tests (Phase 5 integration)
- TestComputeGeneMeta: 5 tests (Phase 6 meta)
- TestRunPhase6: 4 tests (Phase 6 integration)
- TestConfig: 6 tests (config loading)
- TestQCReport: 4 tests (QC system)
- TestOutputWriters: 2 tests (IO)

**Total expected: ~87 tests in test_phases.py**

4. For any MISSING test classes or tests, add them back. Use `str_replace`
   to insert at the appropriate location — do NOT rewrite the entire file.

5. After restoration, run:
```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_phases.py --collect-only -q 2>&1 | tail -5
cd /home/kai001/riker-engine && python -m pytest tests/ -q 2>&1
```

Report both the collected test count and the full suite result.

---

## EXECUTION ORDER

1. Task B first (restore tests)
2. Run full suite to get final count
3. Task A (update README with correct test count)
4. Final verification:
```bash
cd /home/kai001/riker-engine && python -m pytest tests/ -q 2>&1
```

Report the final full-suite output.
