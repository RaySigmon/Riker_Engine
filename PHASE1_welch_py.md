# INSTRUCTION SET FOR KAI — PHASE 1: `riker/stats/welch.py`

## References
- Blueprint Section 5.3 (Cross-Referencing — Welch's t-test usage)
- Blueprint Section 15 (Known Limitations — normal approximation note)
- Context Transfer: "Welch's T-Test Uses Normal Approximation" decision
- Context Transfer: "Kai's custom welch_ttest uses math.erf (normal CDF) instead of the t-distribution"

## CRITICAL REQUIREMENTS (read before writing ANY code)

1. The PRIMARY implementation MUST use `scipy.stats.t.sf()` for the exact t-distribution p-value. This is NON-NEGOTIABLE.
2. A FALLBACK using the normal approximation via `math.erf` is provided ONLY for environments where scipy is unavailable.
3. The Welch-Satterthwaite degrees of freedom formula MUST be used — not pooled df.
4. The function MUST return a named result object (not a bare tuple) containing: t_statistic, p_value, df, mean_diff, se_diff, n1, n2, method (either "t-distribution" or "normal-approximation").
5. Two-sided test is the default. One-sided must be available via parameter.
6. Input validation: reject empty arrays, single-element arrays, zero-variance arrays. Raise ValueError with descriptive messages.
7. All inputs are converted to numpy arrays internally. Accept lists, tuples, pandas Series, or numpy arrays.

---

## FILE: `riker/stats/welch.py`

Write the following file at `/home/kai001/riker-engine/riker/stats/welch.py`:

```python
"""
Riker Engine - Welch's t-test with exact t-distribution.

Computes Welch's unequal variances t-test using the exact t-distribution
via scipy. Falls back to normal approximation only if scipy is unavailable.

References:
    Blueprint Section 5.3, Section 15
    Context Transfer: Welch's T-Test Uses Normal Approximation (decision to fix)
"""

import math
import warnings
from dataclasses import dataclass
from typing import Literal

import numpy as np

# Attempt scipy import; flag availability
try:
    from scipy.stats import t as t_dist
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False


@dataclass(frozen=True)
class WelchResult:
    """Result of a Welch's t-test.

    Attributes:
        t_statistic: The t-test statistic.
        p_value: The p-value (two-sided or one-sided depending on call).
        df: Welch-Satterthwaite degrees of freedom.
        mean_diff: Difference in means (group1 - group2).
        se_diff: Standard error of the difference in means.
        n1: Sample size of group 1.
        n2: Sample size of group 2.
        method: Which distribution was used for the p-value.
    """
    t_statistic: float
    p_value: float
    df: float
    mean_diff: float
    se_diff: float
    n1: int
    n2: int
    method: str


def _validate_group(group: np.ndarray, name: str) -> np.ndarray:
    """Validate and clean a single input group.

    - Converts to float64 numpy array
    - Strips NaN values (with warning if any found)
    - Rejects empty, single-element, or zero-variance arrays
    """
    arr = np.asarray(group, dtype=np.float64).ravel()

    # Strip NaN
    nan_mask = np.isnan(arr)
    if nan_mask.any():
        n_nan = int(nan_mask.sum())
        warnings.warn(
            f"{name}: dropped {n_nan} NaN value(s) from {len(arr)} observations.",
            stacklevel=3,
        )
        arr = arr[~nan_mask]

    if len(arr) == 0:
        raise ValueError(f"{name}: array is empty (after NaN removal).")
    if len(arr) < 2:
        raise ValueError(
            f"{name}: need at least 2 observations, got {len(arr)}."
        )
    if np.var(arr, ddof=1) == 0.0:
        raise ValueError(
            f"{name}: zero variance (all values identical: {arr[0]})."
        )

    return arr


def _welch_satterthwaite_df(var1: float, n1: int, var2: float, n2: int) -> float:
    """Compute Welch-Satterthwaite degrees of freedom.

    df = (s1^2/n1 + s2^2/n2)^2 / ( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    """
    v1 = var1 / n1
    v2 = var2 / n2
    numerator = (v1 + v2) ** 2
    denominator = (v1 ** 2) / (n1 - 1) + (v2 ** 2) / (n2 - 1)

    if denominator == 0.0:
        return float(n1 + n2 - 2)  # fallback to pooled df

    return numerator / denominator


def _p_value_normal(t_stat: float, alternative: str) -> float:
    """Compute p-value using normal approximation (fallback only).

    Uses the complementary error function for numerical stability.
    """
    if alternative == "two-sided":
        return math.erfc(abs(t_stat) / math.sqrt(2))
    elif alternative == "greater":
        return 0.5 * math.erfc(t_stat / math.sqrt(2))
    else:  # less
        return 0.5 * math.erfc(-t_stat / math.sqrt(2))


def _p_value_t(t_stat: float, df: float, alternative: str) -> float:
    """Compute p-value using the exact t-distribution (primary method)."""
    if alternative == "two-sided":
        return float(2.0 * t_dist.sf(abs(t_stat), df))
    elif alternative == "greater":
        return float(t_dist.sf(t_stat, df))
    else:  # less
        return float(t_dist.cdf(t_stat, df))


def welch_ttest(
    group1,
    group2,
    alternative: Literal["two-sided", "greater", "less"] = "two-sided",
) -> WelchResult:
    """Perform Welch's unequal variances t-test.

    Parameters
    ----------
    group1 : array-like
        Expression values for condition group (e.g., cases).
    group2 : array-like
        Expression values for control group.
    alternative : str
        'two-sided' (default), 'greater' (group1 > group2),
        or 'less' (group1 < group2).

    Returns
    -------
    WelchResult
        Dataclass with t_statistic, p_value, df, mean_diff, se_diff,
        n1, n2, and method.

    Raises
    ------
    ValueError
        If inputs are empty, single-element, zero-variance, or
        alternative is invalid.
    """
    if alternative not in ("two-sided", "greater", "less"):
        raise ValueError(
            f"alternative must be 'two-sided', 'greater', or 'less', "
            f"got '{alternative}'."
        )

    # Validate inputs
    g1 = _validate_group(group1, "group1")
    g2 = _validate_group(group2, "group2")

    n1 = len(g1)
    n2 = len(g2)
    mean1 = float(np.mean(g1))
    mean2 = float(np.mean(g2))
    var1 = float(np.var(g1, ddof=1))
    var2 = float(np.var(g2, ddof=1))

    # Standard error of the difference
    se = math.sqrt(var1 / n1 + var2 / n2)

    # T-statistic
    t_stat = (mean1 - mean2) / se

    # Welch-Satterthwaite degrees of freedom
    df = _welch_satterthwaite_df(var1, n1, var2, n2)

    # P-value: exact t-distribution (primary) or normal approximation (fallback)
    if _HAS_SCIPY:
        p_val = _p_value_t(t_stat, df, alternative)
        method = "t-distribution"
    else:
        warnings.warn(
            "scipy not available — using normal approximation for p-values. "
            "This is slightly conservative for small sample sizes (n < 40). "
            "Install scipy for exact t-distribution p-values.",
            UserWarning,
            stacklevel=2,
        )
        p_val = _p_value_normal(t_stat, alternative)
        method = "normal-approximation"

    return WelchResult(
        t_statistic=t_stat,
        p_value=p_val,
        df=df,
        mean_diff=mean1 - mean2,
        se_diff=se,
        n1=n1,
        n2=n2,
        method=method,
    )
```

---

## FILE: `tests/test_stats.py`

**Replace** the contents of `/home/kai001/riker-engine/tests/test_stats.py` with the following. These tests validate welch.py against scipy.stats.ttest_ind as ground truth:

```python
"""
Riker Engine - Statistical module tests.

Phase 1: Welch's t-test validation.
Tests validate against scipy.stats.ttest_ind as ground truth.
"""

import math
import warnings

import numpy as np
import pytest
from scipy.stats import ttest_ind

from riker.stats.welch import WelchResult, welch_ttest


# ---------------------------------------------------------------------------
# 1. Agreement with scipy ground truth
# ---------------------------------------------------------------------------

class TestWelchVsScipy:
    """Verify our implementation matches scipy.stats.ttest_ind(equal_var=False)."""

    def test_basic_two_groups(self):
        """Standard two-group comparison with known different means."""
        np.random.seed(42)
        g1 = np.random.normal(loc=5.0, scale=1.0, size=30)
        g2 = np.random.normal(loc=4.0, scale=1.5, size=25)

        result = welch_ttest(g1, g2)
        scipy_t, scipy_p = ttest_ind(g1, g2, equal_var=False)

        assert result.method == "t-distribution"
        assert math.isclose(result.t_statistic, scipy_t, rel_tol=1e-10)
        assert math.isclose(result.p_value, scipy_p, rel_tol=1e-10)

    def test_equal_groups(self):
        """Two groups drawn from the same distribution — p should be large."""
        np.random.seed(123)
        g1 = np.random.normal(loc=10.0, scale=2.0, size=50)
        g2 = np.random.normal(loc=10.0, scale=2.0, size=50)

        result = welch_ttest(g1, g2)
        scipy_t, scipy_p = ttest_ind(g1, g2, equal_var=False)

        assert math.isclose(result.t_statistic, scipy_t, rel_tol=1e-10)
        assert math.isclose(result.p_value, scipy_p, rel_tol=1e-10)
        assert result.p_value > 0.05  # should not reject null

    def test_very_different_variances(self):
        """Unequal variances — the whole point of Welch's test."""
        np.random.seed(456)
        g1 = np.random.normal(loc=0, scale=1.0, size=40)
        g2 = np.random.normal(loc=0.5, scale=5.0, size=20)

        result = welch_ttest(g1, g2)
        scipy_t, scipy_p = ttest_ind(g1, g2, equal_var=False)

        assert math.isclose(result.t_statistic, scipy_t, rel_tol=1e-10)
        assert math.isclose(result.p_value, scipy_p, rel_tol=1e-10)

    def test_small_samples(self):
        """Small n where normal approximation diverges from t-distribution.
        This is the exact scenario from the Context Transfer document."""
        np.random.seed(789)
        g1 = np.random.normal(loc=3.0, scale=0.5, size=8)
        g2 = np.random.normal(loc=2.5, scale=0.8, size=6)

        result = welch_ttest(g1, g2)
        scipy_t, scipy_p = ttest_ind(g1, g2, equal_var=False)

        assert result.method == "t-distribution"
        assert math.isclose(result.t_statistic, scipy_t, rel_tol=1e-10)
        assert math.isclose(result.p_value, scipy_p, rel_tol=1e-10)

    def test_large_samples(self):
        """Large n where normal and t should converge."""
        np.random.seed(1024)
        g1 = np.random.normal(loc=0, scale=1, size=500)
        g2 = np.random.normal(loc=0.2, scale=1, size=500)

        result = welch_ttest(g1, g2)
        scipy_t, scipy_p = ttest_ind(g1, g2, equal_var=False)

        assert math.isclose(result.t_statistic, scipy_t, rel_tol=1e-10)
        assert math.isclose(result.p_value, scipy_p, rel_tol=1e-10)


# ---------------------------------------------------------------------------
# 2. One-sided alternatives
# ---------------------------------------------------------------------------

class TestWelchAlternatives:
    """Verify one-sided tests agree with scipy."""

    def _get_data(self):
        np.random.seed(42)
        g1 = np.random.normal(loc=5.0, scale=1.0, size=30)
        g2 = np.random.normal(loc=4.0, scale=1.5, size=25)
        return g1, g2

    def test_greater(self):
        g1, g2 = self._get_data()
        result = welch_ttest(g1, g2, alternative="greater")
        scipy_t, scipy_p = ttest_ind(g1, g2, equal_var=False, alternative="greater")
        assert math.isclose(result.p_value, scipy_p, rel_tol=1e-10)

    def test_less(self):
        g1, g2 = self._get_data()
        result = welch_ttest(g1, g2, alternative="less")
        scipy_t, scipy_p = ttest_ind(g1, g2, equal_var=False, alternative="less")
        assert math.isclose(result.p_value, scipy_p, rel_tol=1e-10)

    def test_two_sided_equals_two_times_smaller_one_sided(self):
        """Two-sided p should equal 2x the smaller one-sided p."""
        g1, g2 = self._get_data()
        two = welch_ttest(g1, g2, alternative="two-sided")
        greater = welch_ttest(g1, g2, alternative="greater")
        less = welch_ttest(g1, g2, alternative="less")

        smaller_one_sided = min(greater.p_value, less.p_value)
        assert math.isclose(two.p_value, 2 * smaller_one_sided, rel_tol=1e-10)


# ---------------------------------------------------------------------------
# 3. Return value correctness
# ---------------------------------------------------------------------------

class TestWelchReturnValues:
    """Verify all fields in the WelchResult are correct."""

    def test_result_type(self):
        result = welch_ttest([1, 2, 3, 4, 5], [6, 7, 8, 9, 10])
        assert isinstance(result, WelchResult)

    def test_mean_diff(self):
        g1 = [10, 20, 30]
        g2 = [5, 15, 25]
        result = welch_ttest(g1, g2)
        assert math.isclose(result.mean_diff, 5.0)

    def test_sample_sizes(self):
        g1 = [1, 2, 3]
        g2 = [4, 5, 6, 7, 8]
        result = welch_ttest(g1, g2)
        assert result.n1 == 3
        assert result.n2 == 5

    def test_df_is_welch_satterthwaite(self):
        """Verify df is NOT pooled (n1+n2-2) but Welch-Satterthwaite."""
        np.random.seed(42)
        g1 = np.random.normal(loc=0, scale=1.0, size=20)
        g2 = np.random.normal(loc=0, scale=3.0, size=10)
        result = welch_ttest(g1, g2)

        pooled_df = 20 + 10 - 2
        # Welch-Satterthwaite df should be LESS than pooled for unequal var
        assert result.df < pooled_df
        assert result.df > 0

    def test_se_diff_positive(self):
        result = welch_ttest([1, 2, 3, 4], [5, 6, 7, 8])
        assert result.se_diff > 0


# ---------------------------------------------------------------------------
# 4. Input validation
# ---------------------------------------------------------------------------

class TestWelchInputValidation:
    """Ensure bad inputs are caught with clear errors."""

    def test_empty_group1(self):
        with pytest.raises(ValueError, match="empty"):
            welch_ttest([], [1, 2, 3])

    def test_empty_group2(self):
        with pytest.raises(ValueError, match="empty"):
            welch_ttest([1, 2, 3], [])

    def test_single_element(self):
        with pytest.raises(ValueError, match="at least 2"):
            welch_ttest([1], [2, 3, 4])

    def test_zero_variance(self):
        with pytest.raises(ValueError, match="zero variance"):
            welch_ttest([5, 5, 5], [1, 2, 3])

    def test_invalid_alternative(self):
        with pytest.raises(ValueError, match="alternative"):
            welch_ttest([1, 2, 3], [4, 5, 6], alternative="bigger")


# ---------------------------------------------------------------------------
# 5. NaN handling
# ---------------------------------------------------------------------------

class TestWelchNanHandling:
    """Verify NaN values are stripped with warnings."""

    def test_nan_stripped_with_warning(self):
        g1 = [1.0, 2.0, float("nan"), 4.0, 5.0]
        g2 = [10.0, 20.0, 30.0]

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = welch_ttest(g1, g2)
            assert len(w) == 1
            assert "NaN" in str(w[0].message) or "nan" in str(w[0].message).lower()

        # Should have used 4 observations, not 5
        assert result.n1 == 4

    def test_all_nan_raises(self):
        with pytest.raises(ValueError, match="empty"):
            welch_ttest([float("nan"), float("nan")], [1, 2, 3])


# ---------------------------------------------------------------------------
# 6. Input type flexibility
# ---------------------------------------------------------------------------

class TestWelchInputTypes:
    """Accept lists, tuples, numpy arrays, pandas Series."""

    def test_lists(self):
        result = welch_ttest([1, 2, 3, 4], [5, 6, 7, 8])
        assert isinstance(result, WelchResult)

    def test_tuples(self):
        result = welch_ttest((1, 2, 3, 4), (5, 6, 7, 8))
        assert isinstance(result, WelchResult)

    def test_numpy_arrays(self):
        result = welch_ttest(np.array([1, 2, 3, 4]), np.array([5, 6, 7, 8]))
        assert isinstance(result, WelchResult)

    def test_pandas_series(self):
        import pandas as pd
        result = welch_ttest(pd.Series([1, 2, 3, 4]), pd.Series([5, 6, 7, 8]))
        assert isinstance(result, WelchResult)

    def test_mixed_types(self):
        """List for group1, numpy array for group2."""
        result = welch_ttest([1, 2, 3, 4], np.array([5, 6, 7, 8]))
        assert isinstance(result, WelchResult)
```

---

## EXECUTION INSTRUCTIONS

After writing both files, run:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_stats.py -v 2>&1
```

**Expected result: ALL tests pass.** If any test fails, report the FULL pytest output — do NOT attempt to fix without reporting first.

Then run this confirmation check:

```bash
cd /home/kai001/riker-engine && python -c "
from riker.stats.welch import welch_ttest
result = welch_ttest([1,2,3,4,5], [6,7,8,9,10])
print(f't = {result.t_statistic:.6f}')
print(f'p = {result.p_value:.6f}')
print(f'df = {result.df:.4f}')
print(f'method = {result.method}')
print(f'mean_diff = {result.mean_diff:.2f}')
assert result.method == 't-distribution', 'FAIL: not using exact t-distribution!'
print('PASS: welch.py using exact t-distribution')
"
```

Report both outputs back to the architect for QA review before proceeding.
