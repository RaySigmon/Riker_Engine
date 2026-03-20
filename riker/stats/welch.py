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
