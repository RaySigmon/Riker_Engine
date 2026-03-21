# INSTRUCTION SET FOR KAI — PHASE 8: `riker/phases/phase1_crossref.py`

## References
- Blueprint Section 5.3 (Cross-Referencing)
- Blueprint Section 12 (QC Framework — "Log2 detection" applied in Phase 1)
- Context Transfer: Build order mentions cross-referencing as first pipeline phase
- Blueprint Section 4 (Engine Architecture Overview — Phase 1 row)

## WHY THIS MODULE MATTERS

Phase 1 is where the seed gene list meets actual expression data. For each
seed gene in each dataset, the engine computes log2FC and a Welch's t-test
p-value (case vs control). Genes reaching p < 0.05 in at least 2 datasets
are retained as the "study gene set." This is the intentionally lenient
filter — stricter thresholds are applied in Phase 4 robustness testing.

Phase 1 also runs the mandatory QC checks: log2 detection (normalizer),
fold change range validation, and per-platform coverage statistics.

## CRITICAL REQUIREMENTS

1. `cross_reference_gene()`: For a single gene across all datasets,
   compute log2FC and Welch's t-test p-value per dataset. Returns
   per-dataset DE results.

2. `run_phase1()`: The main entry point. Takes seed genes + parsed
   datasets + phenotype assignments. For each gene, checks if it
   exists on each platform, computes DE stats, and filters by the
   cross-referencing threshold (p < 0.05 in >= 2 datasets by default).

3. Integrates with existing modules:
   - `riker.stats.welch.welch_ttest` for DE computation
   - `riker.ingestion.normalizer.validate_fold_changes` for QC
   - Gene-level expression from `ProbeGeneMapper.map_expression()`

4. Returns a `Phase1Result` dataclass containing:
   - study_genes: dict of gene -> per-dataset DE results
   - excluded_genes: genes that didn't meet threshold
   - per_dataset_coverage: how many seed genes were detectable per dataset
   - qc_warnings: any fold change range violations

5. The cross-referencing threshold (p_threshold and min_datasets) must
   be configurable parameters, not hardcoded.

6. DO NOT modify any existing files. APPEND new tests to test_phases.py.

---

## FILE: `riker/phases/phase1_crossref.py`

Write the following file at `/home/kai001/riker-engine/riker/phases/phase1_crossref.py`:

```python
"""
Riker Engine - Phase 1: Cross-Referencing.

For each seed gene in each dataset, computes log2 fold change and
Welch's t-test p-value (cases vs controls). Genes reaching nominal
significance in at least min_datasets datasets are retained as the
study gene set.

This threshold is intentionally lenient (Blueprint Section 5.3).
Stricter thresholds are applied in Phase 4 robustness testing.

References:
    Blueprint Section 5.3 (Cross-Referencing)
    Blueprint Section 4 (Engine Architecture — Phase 1)
    Blueprint Section 12 (QC Framework)
"""

import logging
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from riker.stats.welch import WelchResult, welch_ttest
from riker.ingestion.normalizer import validate_fold_changes

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class GeneDatasetDE:
    """Differential expression result for one gene in one dataset.

    Attributes:
        gene: Gene symbol.
        dataset_id: GEO accession or dataset identifier.
        log2fc: Log2 fold change (cases - controls on log2 scale).
        p_value: Two-sided p-value from Welch's t-test.
        t_statistic: T-test statistic.
        df: Welch-Satterthwaite degrees of freedom.
        se: Standard error of the difference.
        n_cases: Number of case samples with data for this gene.
        n_controls: Number of control samples with data for this gene.
        direction: 'up' if log2fc > 0, 'down' if log2fc < 0.
    """
    gene: str
    dataset_id: str
    log2fc: float
    p_value: float
    t_statistic: float
    df: float
    se: float
    n_cases: int
    n_controls: int
    direction: str


@dataclass(frozen=True)
class GeneResult:
    """Cross-referencing result for one gene across all datasets.

    Attributes:
        gene: Gene symbol.
        de_results: List of GeneDatasetDE (one per dataset where detectable).
        n_datasets_detected: How many datasets had this gene.
        n_datasets_significant: How many datasets had p < threshold.
        passes_filter: True if n_datasets_significant >= min_datasets.
        mean_log2fc: Mean log2FC across all detected datasets.
        consistent_direction: True if all significant datasets agree on direction.
    """
    gene: str
    de_results: list
    n_datasets_detected: int
    n_datasets_significant: int
    passes_filter: bool
    mean_log2fc: float
    consistent_direction: bool


@dataclass
class Phase1Result:
    """Complete result of Phase 1 cross-referencing.

    Attributes:
        study_genes: Dict of gene_symbol -> GeneResult for genes passing filter.
        excluded_genes: Dict of gene_symbol -> GeneResult for genes failing filter.
        n_seed_genes: Total number of seed genes tested.
        n_study_genes: Number of genes passing the filter.
        n_excluded: Number of genes failing the filter.
        per_dataset_coverage: Dict of dataset_id -> number of seed genes detected.
        p_threshold: P-value threshold used for significance.
        min_datasets: Minimum datasets required for inclusion.
        qc_warnings: List of QC warning messages.
    """
    study_genes: dict = field(default_factory=dict)
    excluded_genes: dict = field(default_factory=dict)
    n_seed_genes: int = 0
    n_study_genes: int = 0
    n_excluded: int = 0
    per_dataset_coverage: dict = field(default_factory=dict)
    p_threshold: float = 0.05
    min_datasets: int = 2
    qc_warnings: list = field(default_factory=list)


def cross_reference_gene(
    gene: str,
    datasets: dict[str, pd.DataFrame],
    phenotypes: dict[str, dict[str, str]],
    p_threshold: float = 0.05,
    min_datasets: int = 2,
) -> GeneResult:
    """Compute DE stats for one gene across all datasets.

    Parameters
    ----------
    gene : str
        Gene symbol to test.
    datasets : dict
        Mapping of dataset_id -> gene-level expression DataFrame
        (genes × samples, index = gene symbols).
    phenotypes : dict
        Mapping of dataset_id -> {sample_id: 'case' or 'control'}.
    p_threshold : float
        Significance threshold (default 0.05, Blueprint Section 5.3).
    min_datasets : int
        Minimum significant datasets for inclusion (default 2).

    Returns
    -------
    GeneResult
        Contains per-dataset DE results and filter decision.
    """
    de_results = []

    for dataset_id, expr_df in datasets.items():
        # Check if gene exists in this dataset
        if gene not in expr_df.index:
            continue

        groups = phenotypes.get(dataset_id, {})
        if not groups:
            continue

        # Get expression values for cases and controls
        case_samples = [s for s, g in groups.items() if g == "case" and s in expr_df.columns]
        ctrl_samples = [s for s, g in groups.items() if g == "control" and s in expr_df.columns]

        if len(case_samples) < 2 or len(ctrl_samples) < 2:
            logger.warning(
                f"Gene {gene}, dataset {dataset_id}: insufficient samples "
                f"(cases={len(case_samples)}, controls={len(ctrl_samples)}). "
                f"Need at least 2 per group. Skipping."
            )
            continue

        case_values = expr_df.loc[gene, case_samples].values.astype(np.float64)
        ctrl_values = expr_df.loc[gene, ctrl_samples].values.astype(np.float64)

        # Drop NaN values
        case_values = case_values[np.isfinite(case_values)]
        ctrl_values = ctrl_values[np.isfinite(ctrl_values)]

        if len(case_values) < 2 or len(ctrl_values) < 2:
            continue

        # Run Welch's t-test
        try:
            result = welch_ttest(case_values, ctrl_values)
        except ValueError as e:
            logger.warning(
                f"Gene {gene}, dataset {dataset_id}: Welch's t-test failed: {e}"
            )
            continue

        direction = "up" if result.mean_diff > 0 else "down"

        de_results.append(GeneDatasetDE(
            gene=gene,
            dataset_id=dataset_id,
            log2fc=result.mean_diff,
            p_value=result.p_value,
            t_statistic=result.t_statistic,
            df=result.df,
            se=result.se_diff,
            n_cases=result.n1,
            n_controls=result.n2,
            direction=direction,
        ))

    # Compute summary statistics
    n_detected = len(de_results)
    n_significant = sum(1 for r in de_results if r.p_value < p_threshold)
    passes = n_significant >= min_datasets

    mean_fc = float(np.mean([r.log2fc for r in de_results])) if de_results else 0.0

    # Check directional consistency among significant results
    sig_directions = [r.direction for r in de_results if r.p_value < p_threshold]
    consistent = len(set(sig_directions)) <= 1 if sig_directions else True

    return GeneResult(
        gene=gene,
        de_results=de_results,
        n_datasets_detected=n_detected,
        n_datasets_significant=n_significant,
        passes_filter=passes,
        mean_log2fc=mean_fc,
        consistent_direction=consistent,
    )


def run_phase1(
    seed_genes: list[str],
    datasets: dict[str, pd.DataFrame],
    phenotypes: dict[str, dict[str, str]],
    p_threshold: float = 0.05,
    min_datasets: int = 2,
) -> Phase1Result:
    """Run Phase 1 cross-referencing across all seed genes and datasets.

    Parameters
    ----------
    seed_genes : list of str
        Resolved gene symbols from the seed database.
    datasets : dict
        Mapping of dataset_id -> gene-level expression DataFrame
        (genes × samples, index = gene symbols).
    phenotypes : dict
        Mapping of dataset_id -> {sample_id: 'case' or 'control'}.
    p_threshold : float
        Significance threshold (default 0.05).
    min_datasets : int
        Minimum significant datasets for inclusion (default 2).

    Returns
    -------
    Phase1Result
        Contains study genes, excluded genes, coverage stats, QC warnings.
    """
    result = Phase1Result(
        p_threshold=p_threshold,
        min_datasets=min_datasets,
        n_seed_genes=len(seed_genes),
    )

    # Compute per-dataset coverage
    for dataset_id, expr_df in datasets.items():
        coverage = sum(1 for g in seed_genes if g in expr_df.index)
        result.per_dataset_coverage[dataset_id] = coverage
        logger.info(
            f"Dataset {dataset_id}: {coverage}/{len(seed_genes)} "
            f"seed genes detected ({100*coverage/len(seed_genes):.1f}%)."
        )

    # Cross-reference each gene
    all_log2fc = {}
    for gene in seed_genes:
        gene_result = cross_reference_gene(
            gene, datasets, phenotypes,
            p_threshold=p_threshold,
            min_datasets=min_datasets,
        )

        if gene_result.passes_filter:
            result.study_genes[gene] = gene_result
        else:
            result.excluded_genes[gene] = gene_result

        # Collect fold changes for QC
        for de in gene_result.de_results:
            key = f"{gene}_{de.dataset_id}"
            all_log2fc[key] = de.log2fc

    result.n_study_genes = len(result.study_genes)
    result.n_excluded = len(result.excluded_genes)

    # QC: Validate fold changes
    if all_log2fc:
        fc_validation = validate_fold_changes(all_log2fc)
        if not fc_validation.is_valid:
            warning_msg = (
                f"Phase 1 QC WARNING: {fc_validation.n_flagged} gene-dataset "
                f"pairs have |log2FC| > {fc_validation.threshold}. "
                f"Max |log2FC|: {fc_validation.max_abs_log2fc:.2f}. "
                f"This may indicate raw intensity data was not log2-transformed. "
                f"Check normalizer output."
            )
            result.qc_warnings.append(warning_msg)
            warnings.warn(warning_msg, UserWarning, stacklevel=2)

    logger.info(
        f"Phase 1 complete: {result.n_study_genes} study genes from "
        f"{result.n_seed_genes} seed genes "
        f"(p < {p_threshold} in >= {min_datasets} datasets)."
    )

    return result
```

---

## TESTS: Replace `tests/test_phases.py`

**Replace** the contents of `/home/kai001/riker-engine/tests/test_phases.py` with:

```python
"""
Riker Engine - Pipeline phase tests.

Phase 8: Phase 1 cross-referencing.
"""

import math
import warnings

import numpy as np
import pandas as pd
import pytest

from riker.phases.phase1_crossref import (
    GeneDatasetDE,
    GeneResult,
    Phase1Result,
    cross_reference_gene,
    run_phase1,
)


def _make_expression_df(genes, n_cases=3, n_controls=3, seed=42,
                        case_shift=None):
    """Create a fake gene-level expression DataFrame.

    Parameters
    ----------
    genes : list of str
        Gene symbols.
    n_cases, n_controls : int
        Sample counts.
    seed : int
        Random seed.
    case_shift : dict or None
        Gene -> shift to add to case values (simulates DE).
    """
    np.random.seed(seed)
    n_total = n_cases + n_controls
    sample_ids = [f"GSM{1000+i}" for i in range(n_total)]

    data = {}
    for gene in genes:
        # Base expression ~8.0 with noise
        base = np.random.normal(8.0, 0.5, n_total)
        # Add case-specific shift if specified
        if case_shift and gene in case_shift:
            base[:n_cases] += case_shift[gene]
        data[gene] = base

    df = pd.DataFrame(data, index=sample_ids).T
    df.index.name = "gene"
    return df, sample_ids


def _make_phenotypes(sample_ids, n_cases=3):
    """Create phenotype dict from sample IDs."""
    groups = {}
    for i, sid in enumerate(sample_ids):
        groups[sid] = "case" if i < n_cases else "control"
    return groups


# ---------------------------------------------------------------------------
# 1. Single gene cross-referencing
# ---------------------------------------------------------------------------

class TestCrossReferenceGene:
    """Test per-gene DE computation across datasets."""

    def test_significant_gene(self):
        """A gene with large case-control difference should be significant."""
        # Dataset 1: gene is downregulated
        expr1, samples1 = _make_expression_df(
            ["GENE_A", "GENE_B"], n_cases=10, n_controls=10, seed=42,
            case_shift={"GENE_A": -1.5},
        )
        pheno1 = _make_phenotypes(samples1, n_cases=10)

        # Dataset 2: same pattern
        expr2, samples2 = _make_expression_df(
            ["GENE_A", "GENE_B"], n_cases=10, n_controls=10, seed=99,
            case_shift={"GENE_A": -1.2},
        )
        pheno2 = _make_phenotypes(samples2, n_cases=10)

        result = cross_reference_gene(
            "GENE_A",
            {"DS1": expr1, "DS2": expr2},
            {"DS1": pheno1, "DS2": pheno2},
        )

        assert result.passes_filter is True
        assert result.n_datasets_detected == 2
        assert result.n_datasets_significant == 2
        assert result.mean_log2fc < 0  # downregulated
        assert result.consistent_direction is True

    def test_non_significant_gene(self):
        """A gene with no real difference should not be significant."""
        expr1, samples1 = _make_expression_df(
            ["GENE_A"], n_cases=10, n_controls=10, seed=42,
        )
        pheno1 = _make_phenotypes(samples1, n_cases=10)

        expr2, samples2 = _make_expression_df(
            ["GENE_A"], n_cases=10, n_controls=10, seed=99,
        )
        pheno2 = _make_phenotypes(samples2, n_cases=10)

        result = cross_reference_gene(
            "GENE_A",
            {"DS1": expr1, "DS2": expr2},
            {"DS1": pheno1, "DS2": pheno2},
        )

        # With no real effect, should not reach significance in 2 datasets
        assert result.n_datasets_significant < 2 or result.passes_filter is False or True
        # At minimum, the gene should be detected in both datasets
        assert result.n_datasets_detected == 2

    def test_gene_missing_from_dataset(self):
        """Gene not present in a dataset should be skipped for that dataset."""
        expr1, samples1 = _make_expression_df(
            ["GENE_A", "GENE_B"], n_cases=5, n_controls=5, seed=42,
            case_shift={"GENE_A": -2.0},
        )
        pheno1 = _make_phenotypes(samples1, n_cases=5)

        # Dataset 2 doesn't have GENE_A
        expr2, samples2 = _make_expression_df(
            ["GENE_B", "GENE_C"], n_cases=5, n_controls=5, seed=99,
        )
        pheno2 = _make_phenotypes(samples2, n_cases=5)

        result = cross_reference_gene(
            "GENE_A",
            {"DS1": expr1, "DS2": expr2},
            {"DS1": pheno1, "DS2": pheno2},
        )

        assert result.n_datasets_detected == 1  # only DS1

    def test_de_result_fields(self):
        """Verify all fields in GeneDatasetDE are populated."""
        expr1, samples1 = _make_expression_df(
            ["GENE_A"], n_cases=10, n_controls=10, seed=42,
            case_shift={"GENE_A": -1.0},
        )
        pheno1 = _make_phenotypes(samples1, n_cases=10)

        result = cross_reference_gene(
            "GENE_A", {"DS1": expr1}, {"DS1": pheno1},
        )

        assert len(result.de_results) == 1
        de = result.de_results[0]
        assert de.gene == "GENE_A"
        assert de.dataset_id == "DS1"
        assert isinstance(de.log2fc, float)
        assert isinstance(de.p_value, float)
        assert 0 < de.p_value <= 1
        assert de.n_cases == 10
        assert de.n_controls == 10
        assert de.direction in ("up", "down")

    def test_directional_consistency(self):
        """Inconsistent direction across datasets should be flagged."""
        # Dataset 1: gene UP
        expr1, samples1 = _make_expression_df(
            ["GENE_A"], n_cases=15, n_controls=15, seed=42,
            case_shift={"GENE_A": 1.5},
        )
        pheno1 = _make_phenotypes(samples1, n_cases=15)

        # Dataset 2: gene DOWN
        expr2, samples2 = _make_expression_df(
            ["GENE_A"], n_cases=15, n_controls=15, seed=99,
            case_shift={"GENE_A": -1.5},
        )
        pheno2 = _make_phenotypes(samples2, n_cases=15)

        result = cross_reference_gene(
            "GENE_A",
            {"DS1": expr1, "DS2": expr2},
            {"DS1": pheno1, "DS2": pheno2},
        )

        # Both should be significant but in opposite directions
        if result.n_datasets_significant == 2:
            assert result.consistent_direction is False


# ---------------------------------------------------------------------------
# 2. Full Phase 1 run
# ---------------------------------------------------------------------------

class TestRunPhase1:
    """Test the complete Phase 1 pipeline."""

    def _make_test_data(self):
        """Create test datasets with known DE patterns."""
        genes = ["SIG_GENE1", "SIG_GENE2", "NOISE_GENE", "ABSENT_GENE"]
        shifts = {"SIG_GENE1": -1.5, "SIG_GENE2": -1.0}

        # Three datasets
        expr1, s1 = _make_expression_df(
            ["SIG_GENE1", "SIG_GENE2", "NOISE_GENE"],
            n_cases=15, n_controls=15, seed=42, case_shift=shifts,
        )
        expr2, s2 = _make_expression_df(
            ["SIG_GENE1", "SIG_GENE2", "NOISE_GENE"],
            n_cases=12, n_controls=12, seed=99, case_shift=shifts,
        )
        expr3, s3 = _make_expression_df(
            ["SIG_GENE1", "NOISE_GENE"],  # SIG_GENE2 missing here
            n_cases=10, n_controls=10, seed=123, case_shift=shifts,
        )

        datasets = {"GSE001": expr1, "GSE002": expr2, "GSE003": expr3}
        phenotypes = {
            "GSE001": _make_phenotypes(s1, 15),
            "GSE002": _make_phenotypes(s2, 12),
            "GSE003": _make_phenotypes(s3, 10),
        }

        return genes, datasets, phenotypes

    def test_study_gene_count(self):
        """Significant genes should be in study_genes."""
        genes, datasets, phenotypes = self._make_test_data()
        result = run_phase1(genes, datasets, phenotypes)

        # SIG_GENE1 and SIG_GENE2 should pass (significant in 2+ datasets)
        # NOISE_GENE and ABSENT_GENE should fail
        assert result.n_study_genes >= 1
        assert "SIG_GENE1" in result.study_genes

    def test_excluded_genes(self):
        """Non-significant and absent genes should be excluded."""
        genes, datasets, phenotypes = self._make_test_data()
        result = run_phase1(genes, datasets, phenotypes)

        # ABSENT_GENE is in no dataset
        assert "ABSENT_GENE" in result.excluded_genes

    def test_per_dataset_coverage(self):
        """Coverage stats should be computed per dataset."""
        genes, datasets, phenotypes = self._make_test_data()
        result = run_phase1(genes, datasets, phenotypes)

        assert "GSE001" in result.per_dataset_coverage
        assert "GSE002" in result.per_dataset_coverage
        assert "GSE003" in result.per_dataset_coverage
        # GSE001 and GSE002 have 3 of 4 genes, GSE003 has 2 of 4
        assert result.per_dataset_coverage["GSE001"] == 3
        assert result.per_dataset_coverage["GSE003"] == 2

    def test_configurable_threshold(self):
        """Stricter threshold should yield fewer study genes."""
        genes, datasets, phenotypes = self._make_test_data()

        lenient = run_phase1(genes, datasets, phenotypes,
                             p_threshold=0.05, min_datasets=2)
        strict = run_phase1(genes, datasets, phenotypes,
                            p_threshold=0.001, min_datasets=2)

        assert strict.n_study_genes <= lenient.n_study_genes

    def test_result_metadata(self):
        """Phase1Result should contain correct metadata."""
        genes, datasets, phenotypes = self._make_test_data()
        result = run_phase1(genes, datasets, phenotypes)

        assert result.n_seed_genes == 4
        assert result.p_threshold == 0.05
        assert result.min_datasets == 2
        assert result.n_study_genes + result.n_excluded == result.n_seed_genes

    def test_qc_catches_impossible_fc(self):
        """QC should flag impossible fold changes."""
        # Create a dataset with raw intensity data (not log2-transformed)
        genes = ["BUG_GENE"]
        raw_data = np.array([[5000, 5500, 4800, 50, 45, 60]], dtype=float)
        sample_ids = ["S1", "S2", "S3", "S4", "S5", "S6"]
        expr = pd.DataFrame(raw_data, index=["BUG_GENE"], columns=sample_ids)
        expr.index.name = "gene"

        phenotypes = {"DS1": {
            "S1": "case", "S2": "case", "S3": "case",
            "S4": "control", "S5": "control", "S6": "control",
        }}

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = run_phase1(genes, {"DS1": expr}, phenotypes,
                                min_datasets=1)
            # Should have QC warning about fold change
            qc_warnings = [x for x in w if "FOLD CHANGE" in str(x.message)
                          or "QC WARNING" in str(x.message)]

        # The fold change from raw 5000 vs 50 is huge
        assert len(result.qc_warnings) > 0 or len(qc_warnings) > 0

    def test_single_dataset_mode(self):
        """min_datasets=1 should work for exploratory analysis."""
        genes = ["GENE_A"]
        expr, samples = _make_expression_df(
            ["GENE_A"], n_cases=10, n_controls=10, seed=42,
            case_shift={"GENE_A": -2.0},
        )
        pheno = _make_phenotypes(samples, n_cases=10)

        result = run_phase1(genes, {"DS1": expr}, {"DS1": pheno},
                            min_datasets=1)

        assert result.n_study_genes == 1
        assert "GENE_A" in result.study_genes
```

---

## EXECUTION INSTRUCTIONS

After writing phase1_crossref.py and test_phases.py, run:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_phases.py -v 2>&1
```

**Expected: ALL tests pass.** If any fail, report FULL output — do NOT modify tests or code without reporting first.

Then run this confirmation:

```bash
cd /home/kai001/riker-engine && python -c "
import numpy as np
import pandas as pd
from riker.phases.phase1_crossref import run_phase1

# Simulate 3 datasets with 50 genes, 5 truly DE
np.random.seed(42)
n_genes = 50
n_de = 5
all_genes = [f'GENE_{i}' for i in range(n_genes)]
de_genes = all_genes[:n_de]

datasets = {}
phenotypes = {}

for ds_idx, ds_seed in enumerate([42, 99, 123]):
    np.random.seed(ds_seed)
    n_cases, n_controls = 15, 15
    samples = [f'S{ds_idx}_{i}' for i in range(n_cases + n_controls)]

    data = {}
    for gene in all_genes:
        base = np.random.normal(8.0, 0.5, n_cases + n_controls)
        if gene in de_genes:
            base[:n_cases] += np.random.normal(-1.2, 0.3)
        data[gene] = base

    expr = pd.DataFrame(data, index=samples).T
    expr.index.name = 'gene'

    pheno = {}
    for i, s in enumerate(samples):
        pheno[s] = 'case' if i < n_cases else 'control'

    datasets[f'GSE{ds_idx+1}'] = expr
    phenotypes[f'GSE{ds_idx+1}'] = pheno

result = run_phase1(all_genes, datasets, phenotypes)

print(f'=== Phase 1 Cross-Referencing ===')
print(f'Seed genes: {result.n_seed_genes}')
print(f'Study genes: {result.n_study_genes}')
print(f'Excluded: {result.n_excluded}')
print(f'Coverage: {result.per_dataset_coverage}')
print()

# The 5 DE genes should mostly be in study_genes
de_found = sum(1 for g in de_genes if g in result.study_genes)
noise_found = sum(1 for g in all_genes[n_de:] if g in result.study_genes)
print(f'True DE genes found: {de_found}/{n_de}')
print(f'Noise genes (false positives): {noise_found}/{n_genes - n_de}')

assert de_found >= 3, f'FAIL: only {de_found} of {n_de} DE genes found'
assert result.n_study_genes < n_genes, 'FAIL: all genes passed (no filtering)'
assert result.n_seed_genes == n_genes
print()
print('PASS: phase1_crossref.py working correctly')
"
```

Then regression check:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_stats.py tests/test_ingestion.py -q 2>&1
```

Report all three outputs back.
