# INSTRUCTION SET FOR KAI — PHASE 6: `riker/ingestion/normalizer.py`

## References
- Blueprint Section 5.2 (Transcriptomic Data Ingestion — log2 transformation detection)
- Blueprint Section 5.4 (Expression Scale Handling)
- Blueprint Section 12 (QC Framework — "Log2 detection" row)
- Context Transfer: "Log2 transformation detection is critical"
- Context Transfer: "failure to detect raw intensity data produced log2FC values of 119"

## WHY THIS MODULE MATTERS

During Project Riker development, failing to detect that expression data was in
raw intensity format (not log2-transformed) produced log2FC values of 119 — a
biologically impossible result that passed all automated tests at the time.
This module detects whether data needs log2 transformation, applies it when
needed, and includes mandatory range checks that would catch impossible fold
changes before they propagate.

The expression value 119 as a log2FC means one group has ~6.6 × 10^35 times
the expression of the other. This is obviously wrong, but without a range
check, the pipeline happily computed statistics on it.

## CRITICAL REQUIREMENTS

1. `detect_log2_status()`: Examine expression values and determine if they
   appear to be raw intensities (need log2) or already log2-transformed.
   Rule from Blueprint Section 5.2: if median >= 20, data is likely raw
   intensities and needs log2(x+1) transformation.

2. `normalize_expression()`: Apply log2(x+1) transformation if needed.
   Returns the transformed data and a flag indicating what was done.

3. `validate_fold_changes()`: Check computed log2FC values. ANY absolute
   log2FC > 10 triggers a warning. This is the range check from Blueprint
   Section 5.2 that would have caught the 119 bug.

4. `detect_expression_scale()`: Integration with meta.py's check_expression_scale.
   Flags datasets with expression range max < 5.0 as potential log-ratio
   (GSE33000 issue, Blueprint Section 5.4).

5. All functions work on numpy arrays and pandas DataFrames.

6. DO NOT modify any existing files. APPEND new tests to test_ingestion.py.

---

## FILE: `riker/ingestion/normalizer.py`

Write the following file at `/home/kai001/riker-engine/riker/ingestion/normalizer.py`:

```python
"""
Riker Engine - Expression data normalization and validation.

Detects whether expression data needs log2 transformation, applies it
when needed, and validates computed fold changes against biologically
plausible ranges.

During Project Riker development, failure to detect raw intensity data
produced log2FC values of 119, which passed all automated tests but
were biologically impossible. The range checks in this module prevent
that class of error.

References:
    Blueprint Section 5.2 (log2 transformation detection)
    Blueprint Section 5.4 (Expression scale handling)
    Blueprint Section 12 (QC Framework — log2 detection check)
"""

import logging
import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Thresholds from the Blueprint
LOG2_DETECTION_MEDIAN_THRESHOLD = 20.0  # Section 5.2: median >= 20 means raw
FOLD_CHANGE_MAX = 10.0  # Section 5.2: |log2FC| > 10 triggers warning
LOG_RATIO_RANGE_THRESHOLD = 5.0  # Section 5.4: range max < 5 means log-ratio


@dataclass(frozen=True)
class NormalizationResult:
    """Result of expression normalization.

    Attributes:
        data: The (possibly transformed) expression matrix.
        was_transformed: True if log2(x+1) was applied.
        original_median: Median of the original data.
        transformed_median: Median after transformation (same as original
            if no transformation was applied).
        n_negative_values: Count of negative values in original data
            (negative values cannot be log2-transformed).
        detection_reason: Human-readable explanation of the decision.
    """
    data: np.ndarray
    was_transformed: bool
    original_median: float
    transformed_median: float
    n_negative_values: int
    detection_reason: str


@dataclass(frozen=True)
class FoldChangeValidation:
    """Result of fold change range validation.

    Attributes:
        n_checked: Number of genes checked.
        n_flagged: Number of genes with |log2FC| > threshold.
        max_abs_log2fc: Maximum absolute log2FC observed.
        flagged_genes: List of (gene, log2fc) tuples for flagged genes.
        is_valid: True if no genes were flagged (all within range).
        threshold: The threshold used for flagging.
    """
    n_checked: int
    n_flagged: int
    max_abs_log2fc: float
    flagged_genes: list
    is_valid: bool
    threshold: float


def detect_log2_status(expression_values) -> tuple[bool, str]:
    """Detect whether expression data needs log2 transformation.

    Parameters
    ----------
    expression_values : array-like or DataFrame
        Expression values (genes × samples or flat array).

    Returns
    -------
    tuple of (bool, str)
        (needs_log2, reason)
        needs_log2 is True if data appears to be raw intensities.
        reason is a human-readable explanation.

    Detection Logic (Blueprint Section 5.2):
    - If median >= 20: likely raw intensities → needs log2(x+1)
    - If median < 20 and no negative values: likely already log2
    - If negative values present: already log2 (raw intensities are non-negative)
    """
    vals = _to_flat_array(expression_values)

    if len(vals) == 0:
        return False, "Empty data — no transformation applied."

    median_val = float(np.nanmedian(vals))
    has_negatives = bool(np.any(vals < 0))

    if has_negatives:
        return False, (
            f"Data contains negative values — already log-transformed. "
            f"Median: {median_val:.2f}."
        )

    if median_val >= LOG2_DETECTION_MEDIAN_THRESHOLD:
        return True, (
            f"Median ({median_val:.2f}) >= {LOG2_DETECTION_MEDIAN_THRESHOLD}. "
            f"Data appears to be raw intensities — log2(x+1) recommended."
        )

    return False, (
        f"Median ({median_val:.2f}) < {LOG2_DETECTION_MEDIAN_THRESHOLD}. "
        f"Data appears to be already log2-transformed."
    )


def normalize_expression(
    expression_data,
    force_log2: bool | None = None,
) -> NormalizationResult:
    """Normalize expression data, applying log2 if needed.

    Parameters
    ----------
    expression_data : array-like, DataFrame, or numpy array
        Expression matrix (genes × samples) or flat values.
    force_log2 : bool or None
        If True, force log2(x+1) transformation regardless of detection.
        If False, skip transformation regardless of detection.
        If None (default), auto-detect using detect_log2_status().

    Returns
    -------
    NormalizationResult
        Contains transformed data and metadata about what was done.
    """
    if isinstance(expression_data, pd.DataFrame):
        vals = expression_data.values.astype(np.float64)
        is_dataframe = True
    else:
        vals = np.asarray(expression_data, dtype=np.float64)
        is_dataframe = False

    original_median = float(np.nanmedian(vals))
    n_negative = int(np.sum(vals < 0))

    # Decide whether to transform
    if force_log2 is None:
        needs_log2, reason = detect_log2_status(vals)
    elif force_log2:
        needs_log2 = True
        reason = "Log2 transformation forced by user."
    else:
        needs_log2 = False
        reason = "Log2 transformation skipped by user."

    if needs_log2:
        if n_negative > 0:
            warnings.warn(
                f"Cannot apply log2(x+1): {n_negative} negative values found. "
                f"Negative values would produce NaN. Skipping transformation.",
                UserWarning,
                stacklevel=2,
            )
            transformed = vals.copy()
            reason += " WARNING: Skipped due to negative values."
            was_transformed = False
        else:
            # log2(x + 1) — the +1 handles zeros
            transformed = np.log2(vals + 1)
            was_transformed = True
            logger.info(
                f"Applied log2(x+1) transformation. "
                f"Median: {original_median:.2f} -> {float(np.nanmedian(transformed)):.2f}"
            )
    else:
        transformed = vals.copy()
        was_transformed = False

    transformed_median = float(np.nanmedian(transformed))

    return NormalizationResult(
        data=transformed,
        was_transformed=was_transformed,
        original_median=original_median,
        transformed_median=transformed_median,
        n_negative_values=n_negative,
        detection_reason=reason,
    )


def validate_fold_changes(
    gene_log2fc: dict[str, float],
    threshold: float = FOLD_CHANGE_MAX,
) -> FoldChangeValidation:
    """Validate computed log2FC values against biologically plausible range.

    This is the critical QC check that would have caught the log2FC=119 bug
    during Project Riker development (Blueprint Section 5.2, Section 12).

    Parameters
    ----------
    gene_log2fc : dict
        Mapping of gene_symbol -> log2FC value.
    threshold : float
        Maximum allowed absolute log2FC. Default 10.0 (Blueprint Section 5.2).
        |log2FC| > 10 means > 1024-fold change, which is biologically
        implausible for most transcriptomic comparisons.

    Returns
    -------
    FoldChangeValidation
        Contains flagged genes and validation status.
    """
    if len(gene_log2fc) == 0:
        return FoldChangeValidation(
            n_checked=0,
            n_flagged=0,
            max_abs_log2fc=0.0,
            flagged_genes=[],
            is_valid=True,
            threshold=threshold,
        )

    flagged = []
    max_abs = 0.0

    for gene, fc in gene_log2fc.items():
        abs_fc = abs(fc)
        if abs_fc > max_abs:
            max_abs = abs_fc
        if abs_fc > threshold:
            flagged.append((gene, fc))

    if flagged:
        flagged.sort(key=lambda x: abs(x[1]), reverse=True)
        warnings.warn(
            f"FOLD CHANGE RANGE CHECK FAILED: {len(flagged)} gene(s) have "
            f"|log2FC| > {threshold}. Maximum: {max_abs:.2f}. "
            f"This likely indicates raw intensity data was not log2-transformed. "
            f"Top offenders: {flagged[:5]}. "
            f"See Blueprint Section 5.2.",
            UserWarning,
            stacklevel=2,
        )

    return FoldChangeValidation(
        n_checked=len(gene_log2fc),
        n_flagged=len(flagged),
        max_abs_log2fc=max_abs,
        flagged_genes=flagged,
        is_valid=(len(flagged) == 0),
        threshold=threshold,
    )


def _to_flat_array(data) -> np.ndarray:
    """Convert various input types to a flat float64 array."""
    if isinstance(data, pd.DataFrame):
        arr = data.values.astype(np.float64).ravel()
    elif isinstance(data, pd.Series):
        arr = data.values.astype(np.float64).ravel()
    else:
        arr = np.asarray(data, dtype=np.float64).ravel()

    # Remove NaN for statistics
    return arr[np.isfinite(arr)]
```

---

## TESTS: APPEND to `tests/test_ingestion.py`

**APPEND** the following to the END of `/home/kai001/riker-engine/tests/test_ingestion.py`.
Do NOT delete or modify any existing test classes.

```python


# ===========================================================================
# PHASE 6: NORMALIZER TESTS
# ===========================================================================

import warnings

import numpy as np
import pandas as pd

from riker.ingestion.normalizer import (
    FoldChangeValidation,
    NormalizationResult,
    detect_log2_status,
    normalize_expression,
    validate_fold_changes,
)


# ---------------------------------------------------------------------------
# 3. Log2 detection
# ---------------------------------------------------------------------------

class TestLog2Detection:
    """Verify log2 transformation detection logic."""

    def test_raw_intensities_detected(self):
        """High median (>= 20) should be flagged as raw."""
        # Typical raw microarray: values in hundreds to thousands
        vals = np.random.uniform(100, 10000, size=1000)
        needs_log2, reason = detect_log2_status(vals)
        assert needs_log2 is True
        assert "raw intensities" in reason.lower()

    def test_log2_data_not_flagged(self):
        """Low median (< 20) without negatives should not be flagged."""
        # Typical log2-transformed: values 3-16
        vals = np.random.uniform(3.0, 16.0, size=1000)
        needs_log2, reason = detect_log2_status(vals)
        assert needs_log2 is False

    def test_negative_values_mean_already_log2(self):
        """Presence of negative values means already transformed."""
        vals = np.random.normal(0, 2, size=1000)  # centered near 0, has negatives
        needs_log2, reason = detect_log2_status(vals)
        assert needs_log2 is False
        assert "negative" in reason.lower()

    def test_empty_data(self):
        needs_log2, reason = detect_log2_status(np.array([]))
        assert needs_log2 is False

    def test_median_exactly_at_threshold(self):
        """Median exactly 20 should trigger transformation."""
        vals = np.array([20.0, 20.0, 20.0])
        needs_log2, _ = detect_log2_status(vals)
        assert needs_log2 is True

    def test_dataframe_input(self):
        """Should accept pandas DataFrame."""
        df = pd.DataFrame(np.random.uniform(100, 5000, size=(10, 5)))
        needs_log2, _ = detect_log2_status(df)
        assert needs_log2 is True


# ---------------------------------------------------------------------------
# 4. Expression normalization
# ---------------------------------------------------------------------------

class TestNormalization:
    """Verify expression normalization applies log2 correctly."""

    def test_raw_data_transformed(self):
        """Raw intensity data should be log2(x+1) transformed."""
        raw = np.array([0, 1, 100, 1000, 10000], dtype=float)
        result = normalize_expression(raw)
        assert result.was_transformed is True
        # log2(1001) ≈ 9.97, log2(10001) ≈ 13.29
        assert result.data[-1] == pytest.approx(np.log2(10001))

    def test_log2_data_untouched(self):
        """Already log2-transformed data should not be changed."""
        log2_data = np.random.uniform(3, 16, size=100)
        result = normalize_expression(log2_data)
        assert result.was_transformed is False
        np.testing.assert_array_equal(result.data, log2_data)

    def test_force_log2(self):
        """force_log2=True should transform regardless of detection."""
        data = np.array([5.0, 8.0, 12.0])  # looks like log2
        result = normalize_expression(data, force_log2=True)
        assert result.was_transformed is True
        assert result.data[0] == pytest.approx(np.log2(6.0))

    def test_force_no_log2(self):
        """force_log2=False should skip regardless of detection."""
        data = np.array([100.0, 500.0, 1000.0])  # looks like raw
        result = normalize_expression(data, force_log2=False)
        assert result.was_transformed is False

    def test_negative_values_prevent_transform(self):
        """Negative values should prevent log2 even if median is high."""
        data = np.array([100.0, 200.0, -5.0, 300.0])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = normalize_expression(data, force_log2=True)
            assert len(w) >= 1
            assert "negative" in str(w[0].message).lower()
        assert result.was_transformed is False

    def test_zeros_handled(self):
        """Zeros should be handled by log2(0+1) = 0."""
        data = np.array([0.0, 0.0, 100.0])
        result = normalize_expression(data)
        assert result.data[0] == 0.0  # log2(1) = 0
        assert result.data[1] == 0.0

    def test_dataframe_input(self):
        """Should accept and return numpy array from DataFrame."""
        df = pd.DataFrame(np.random.uniform(100, 5000, size=(5, 3)))
        result = normalize_expression(df)
        assert result.was_transformed is True
        assert isinstance(result.data, np.ndarray)

    def test_metadata_populated(self):
        data = np.array([100.0, 200.0, 300.0])
        result = normalize_expression(data)
        assert result.original_median == 200.0
        assert result.transformed_median < 200.0  # log2 compresses
        assert result.n_negative_values == 0


# ---------------------------------------------------------------------------
# 5. Fold change validation (the 119 bug catcher)
# ---------------------------------------------------------------------------

class TestFoldChangeValidation:
    """Verify the range check that catches impossible fold changes."""

    def test_normal_fold_changes_pass(self):
        """Typical log2FC values (-2 to 2) should all pass."""
        fc = {f"GENE{i}": np.random.uniform(-2, 2) for i in range(100)}
        result = validate_fold_changes(fc)
        assert result.is_valid is True
        assert result.n_flagged == 0

    def test_impossible_fold_change_caught(self):
        """The exact bug from development: log2FC of 119."""
        fc = {"NORMAL_GENE": -0.5, "BUG_GENE": 119.0, "ANOTHER_NORMAL": 1.2}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = validate_fold_changes(fc)
            assert len(w) == 1
            assert "FOLD CHANGE RANGE CHECK FAILED" in str(w[0].message)
        assert result.is_valid is False
        assert result.n_flagged == 1
        assert result.flagged_genes[0] == ("BUG_GENE", 119.0)
        assert result.max_abs_log2fc == 119.0

    def test_threshold_boundary(self):
        """log2FC of exactly 10 should NOT be flagged (boundary is >)."""
        fc = {"GENE1": 10.0}
        result = validate_fold_changes(fc)
        assert result.is_valid is True

    def test_threshold_exceeded(self):
        """log2FC of 10.1 should be flagged."""
        fc = {"GENE1": 10.1}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = validate_fold_changes(fc)
        assert result.is_valid is False

    def test_negative_extreme_caught(self):
        """Large negative fold changes should also be caught."""
        fc = {"GENE1": -15.0}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = validate_fold_changes(fc)
        assert result.is_valid is False
        assert result.flagged_genes[0] == ("GENE1", -15.0)

    def test_multiple_offenders_sorted(self):
        """Multiple flagged genes should be sorted by |log2FC| descending."""
        fc = {"A": 50.0, "B": -100.0, "C": 12.0, "D": 0.5}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = validate_fold_changes(fc)
        assert result.n_flagged == 3
        assert result.flagged_genes[0][0] == "B"  # -100 is largest absolute
        assert result.flagged_genes[1][0] == "A"  # 50 is next

    def test_custom_threshold(self):
        """Custom threshold should be respected."""
        fc = {"GENE1": 5.5}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = validate_fold_changes(fc, threshold=5.0)
        assert result.is_valid is False
        assert result.threshold == 5.0

    def test_empty_input(self):
        result = validate_fold_changes({})
        assert result.is_valid is True
        assert result.n_checked == 0
```

---

## EXECUTION INSTRUCTIONS

After writing normalizer.py and appending tests, run:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_ingestion.py -v 2>&1
```

**Expected: ALL tests pass (Phase 5 gene_db + Phase 6 normalizer).** If any fail, report FULL output — do NOT modify tests or code without reporting first.

Then run this confirmation:

```bash
cd /home/kai001/riker-engine && python -c "
import numpy as np
from riker.ingestion.normalizer import (
    detect_log2_status,
    normalize_expression,
    validate_fold_changes,
)

# Scenario 1: Raw intensity data (the bug scenario)
print('=== Scenario 1: Raw Intensity Detection ===')
raw_data = np.random.uniform(50, 15000, size=1000)
needs_log2, reason = detect_log2_status(raw_data)
print(f'Needs log2: {needs_log2}')
print(f'Reason: {reason}')
assert needs_log2, 'FAIL: did not detect raw intensities!'

# Scenario 2: Apply normalization
result = normalize_expression(raw_data)
print(f'Transformed: {result.was_transformed}')
print(f'Median before: {result.original_median:.1f}')
print(f'Median after:  {result.transformed_median:.2f}')

# Scenario 3: The 119 bug
print()
print('=== Scenario 2: The 119 Bug Detection ===')
import warnings
buggy_fc = {'GENE_A': -0.45, 'GENE_B': 119.0, 'GENE_C': 0.8}
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter('always')
    validation = validate_fold_changes(buggy_fc)
    assert len(w) == 1, 'FAIL: no warning emitted!'
    print(f'Warning: {str(w[0].message)[:80]}...')

assert not validation.is_valid, 'FAIL: validation should have failed!'
assert validation.n_flagged == 1
print(f'Flagged: {validation.flagged_genes}')
print(f'Max |log2FC|: {validation.max_abs_log2fc}')

# Scenario 4: Normal fold changes
print()
print('=== Scenario 3: Normal Fold Changes ===')
normal_fc = {f'G{i}': np.random.uniform(-2, 2) for i in range(500)}
validation = validate_fold_changes(normal_fc)
assert validation.is_valid, 'FAIL: normal FCs should pass!'
print(f'Checked: {validation.n_checked}, Flagged: {validation.n_flagged}')
print(f'Max |log2FC|: {validation.max_abs_log2fc:.4f}')

print()
print('PASS: normalizer.py working correctly')
"
```

Then regression check:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_stats.py -q 2>&1
```

Report all three outputs back.
