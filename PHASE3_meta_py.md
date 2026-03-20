# INSTRUCTION SET FOR KAI — PHASE 3: `riker/stats/meta.py`

## References
- Blueprint Section 10 (Phase 6: Meta-Analysis)
- Blueprint Section 5.4 (Expression Scale Handling)
- Blueprint Section 15 (Welch's t-test approximation note)
- Context Transfer: "GSE33000 Log-Ratio Scale Issue"
- Context Transfer: "Standard errors are recovered from log2FC and p-values via the t-distribution"

## WHY THIS MODULE MATTERS

Phase 6 meta-analysis pools effect sizes across all brain datasets to produce
the final verdict on each gene. SE recovery from log2FC + p-value + sample sizes
is the bridge between per-dataset DE results and pooled estimates. The GSE33000
scale issue (log-ratio values centered near zero, max ~2.0) means SEs must be
recovered on the native scale — and the engine must detect and flag datasets
with unusual expression ranges before pooling.

## CRITICAL REQUIREMENTS

1. SE recovery MUST use the t-distribution (via scipy.stats.t.isf), NOT the
   normal distribution. This is consistent with welch.py using exact t-distribution.
   Formula: SE = abs(log2FC) / abs(t_stat), where t_stat is recovered from
   the p-value and df using the inverse survival function.

2. Two meta-analysis models are REQUIRED:
   - Fixed-effects (inverse-variance weighted): weight_i = 1/SE_i^2
   - DerSimonian-Laird random-effects: incorporates between-study variance (tau^2)
   The RANDOM-EFFECTS model is the PRIMARY result (Blueprint Section 10).

3. Heterogeneity statistics are REQUIRED: Cochran's Q, I-squared, tau-squared.

4. Scale detection: flag any dataset where the expression value range maximum
   is below 5.0 as potential log-ratio format (Blueprint Section 5.4).

5. All inputs/outputs use dataclasses with clear field documentation.

6. DO NOT modify welch.py or fdr.py. DO NOT modify existing tests.
   APPEND new tests to the end of test_stats.py.

---

## FILE: `riker/stats/meta.py`

Write the following file at `/home/kai001/riker-engine/riker/stats/meta.py`:

```python
"""
Riker Engine - Inverse-variance weighted meta-analysis.

Implements fixed-effects and DerSimonian-Laird random-effects meta-analysis
for pooling gene-level effect sizes across multiple transcriptomic datasets.
The random-effects model is the primary result (Blueprint Section 10).

Standard errors are recovered from log2FC, p-values, and sample sizes using
the exact t-distribution, consistent with the Welch's t-test implementation.

References:
    Blueprint Section 10 (Phase 6: Meta-Analysis)
    Blueprint Section 5.4 (Expression Scale Handling)
    Context Transfer: "GSE33000 Log-Ratio Scale Issue"
"""

import math
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.stats import norm, t as t_dist


@dataclass(frozen=True)
class StudyEffect:
    """Per-study effect size with metadata for meta-analysis.

    Attributes:
        dataset_id: GEO accession or dataset identifier.
        log2fc: Log2 fold change (cases minus controls on log2 scale).
        se: Standard error of the log2FC estimate.
        n_cases: Number of case samples.
        n_controls: Number of control samples.
        p_value: Original p-value from Welch's t-test.
        tissue: Tissue type (e.g., 'brain', 'blood').
    """
    dataset_id: str
    log2fc: float
    se: float
    n_cases: int
    n_controls: int
    p_value: float
    tissue: str = "brain"


@dataclass(frozen=True)
class MetaResult:
    """Result of meta-analysis for a single gene.

    Attributes:
        gene: Gene symbol.
        pooled_log2fc: Pooled effect size estimate.
        pooled_se: Standard error of the pooled estimate.
        pooled_z: Z-statistic for the pooled estimate.
        pooled_p: P-value for the pooled estimate (two-sided).
        ci_lower: Lower bound of 95% confidence interval.
        ci_upper: Upper bound of 95% confidence interval.
        model: 'fixed' or 'random' (DerSimonian-Laird).
        cochran_q: Cochran's Q statistic for heterogeneity.
        q_p_value: P-value for Cochran's Q (chi-squared test).
        i_squared: I-squared heterogeneity percentage (0-100).
        tau_squared: Between-study variance estimate.
        n_studies: Number of studies included.
        study_weights: Per-study weights (normalized to sum to 1).
        study_ids: Dataset IDs in the order they were provided.
    """
    gene: str
    pooled_log2fc: float
    pooled_se: float
    pooled_z: float
    pooled_p: float
    ci_lower: float
    ci_upper: float
    model: str
    cochran_q: float
    q_p_value: float
    i_squared: float
    tau_squared: float
    n_studies: int
    study_weights: list
    study_ids: list


@dataclass(frozen=True)
class ScaleCheck:
    """Result of expression scale detection for a dataset.

    Attributes:
        dataset_id: GEO accession.
        expression_range_max: Maximum absolute expression value observed.
        is_log_ratio: True if range_max < 5.0 (potential log-ratio format).
        warning: Human-readable warning message if flagged.
    """
    dataset_id: str
    expression_range_max: float
    is_log_ratio: bool
    warning: str


def recover_se(
    log2fc: float,
    p_value: float,
    n_cases: int,
    n_controls: int,
) -> float:
    """Recover standard error from log2FC, p-value, and sample sizes.

    Uses the t-distribution with Welch-Satterthwaite approximated degrees
    of freedom (approximated as n_cases + n_controls - 2 for recovery
    purposes, since we don't have the per-group variances here).

    For very small p-values (< 1e-300), clips to 1e-300 to avoid
    numerical issues with the inverse survival function.

    Parameters
    ----------
    log2fc : float
        Log2 fold change.
    p_value : float
        Two-sided p-value from Welch's t-test.
    n_cases : int
        Number of case samples.
    n_controls : int
        Number of control samples.

    Returns
    -------
    float
        Recovered standard error.

    Raises
    ------
    ValueError
        If p_value is not in (0, 1], or sample sizes are < 2.
    """
    if not (0.0 < p_value <= 1.0):
        raise ValueError(
            f"p_value must be in (0, 1], got {p_value}."
        )
    if n_cases < 2 or n_controls < 2:
        raise ValueError(
            f"Need n >= 2 per group, got n_cases={n_cases}, "
            f"n_controls={n_controls}."
        )

    # Degrees of freedom (approximation for SE recovery)
    df = n_cases + n_controls - 2

    # Clip very small p-values to avoid inf from inverse survival
    p_clipped = max(p_value, 1e-300)

    # Recover |t| from two-sided p-value using t-distribution
    # Two-sided p = 2 * P(T > |t|), so P(T > |t|) = p/2
    abs_t = float(t_dist.isf(p_clipped / 2.0, df))

    if abs_t == 0.0:
        # p_value is 1.0 or very close — effect is zero
        # Return a large SE to downweight in meta-analysis
        return abs(log2fc) * 100.0 if log2fc != 0.0 else 1.0

    se = abs(log2fc) / abs_t

    # Guard against zero SE when log2fc is exactly zero
    if se == 0.0:
        # log2fc is zero — recover SE from t-distribution quantile alone
        # Use the minimum detectable effect at this sample size
        se = 1.0 / abs_t

    return se


def check_expression_scale(
    dataset_id: str,
    expression_values: np.ndarray,
    threshold: float = 5.0,
) -> ScaleCheck:
    """Detect potential log-ratio expression scale.

    Datasets using Rosetta/Merck arrays or similar platforms produce
    log-ratio values centered near zero (max ~2.0) rather than
    log2-intensity values (typical range 3-16). This affects
    interpretation of effect sizes in meta-analysis.

    Parameters
    ----------
    dataset_id : str
        GEO accession for the dataset.
    expression_values : array-like
        All expression values from the dataset.
    threshold : float
        Maximum absolute value below which the dataset is flagged.
        Default 5.0 (per Blueprint Section 5.4).

    Returns
    -------
    ScaleCheck
        Contains dataset_id, range_max, is_log_ratio flag, and warning.
    """
    vals = np.asarray(expression_values, dtype=np.float64).ravel()
    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        return ScaleCheck(
            dataset_id=dataset_id,
            expression_range_max=0.0,
            is_log_ratio=True,
            warning=f"{dataset_id}: no finite expression values found.",
        )

    range_max = float(np.max(np.abs(vals)))
    is_log_ratio = range_max < threshold

    if is_log_ratio:
        warning = (
            f"{dataset_id}: expression range max ({range_max:.2f}) is below "
            f"{threshold}. This dataset may use log-ratio format "
            f"(e.g., Rosetta/Merck arrays). Standard errors are recovered "
            f"on the native scale, but verify z-score distributions before "
            f"pooling. See Blueprint Section 5.4."
        )
        warnings.warn(warning, UserWarning, stacklevel=2)
    else:
        warning = ""

    return ScaleCheck(
        dataset_id=dataset_id,
        expression_range_max=range_max,
        is_log_ratio=is_log_ratio,
        warning=warning,
    )


def _validate_studies(studies: list[StudyEffect]) -> None:
    """Validate study inputs for meta-analysis."""
    if len(studies) < 2:
        raise ValueError(
            f"Meta-analysis requires at least 2 studies, got {len(studies)}."
        )
    for s in studies:
        if s.se <= 0:
            raise ValueError(
                f"Study '{s.dataset_id}' has non-positive SE ({s.se})."
            )


def fixed_effects_meta(
    gene: str,
    studies: list[StudyEffect],
) -> MetaResult:
    """Fixed-effects inverse-variance weighted meta-analysis.

    Parameters
    ----------
    gene : str
        Gene symbol being analyzed.
    studies : list of StudyEffect
        Per-study effect sizes and standard errors.

    Returns
    -------
    MetaResult
        Pooled estimate under the fixed-effects model.
    """
    _validate_studies(studies)

    effects = np.array([s.log2fc for s in studies])
    ses = np.array([s.se for s in studies])

    # Weights: w_i = 1 / SE_i^2
    weights = 1.0 / (ses ** 2)
    total_weight = np.sum(weights)

    # Pooled estimate
    pooled_effect = float(np.sum(weights * effects) / total_weight)
    pooled_se = float(math.sqrt(1.0 / total_weight))

    # Z-test for pooled estimate
    pooled_z = pooled_effect / pooled_se
    pooled_p = float(2.0 * norm.sf(abs(pooled_z)))

    # 95% CI
    ci_lower = pooled_effect - 1.96 * pooled_se
    ci_upper = pooled_effect + 1.96 * pooled_se

    # Heterogeneity: Cochran's Q
    cochran_q = float(np.sum(weights * (effects - pooled_effect) ** 2))
    k = len(studies)
    q_df = k - 1
    q_p_value = float(1.0 - _chi2_cdf(cochran_q, q_df)) if q_df > 0 else 1.0

    # I-squared
    i_squared = max(0.0, (cochran_q - q_df) / cochran_q * 100.0) if cochran_q > 0 else 0.0

    # Tau-squared (DL estimator, reported even in fixed model for reference)
    tau_sq = _estimate_tau_squared(cochran_q, weights, k)

    # Normalized weights for reporting
    norm_weights = [float(w / total_weight) for w in weights]

    return MetaResult(
        gene=gene,
        pooled_log2fc=pooled_effect,
        pooled_se=pooled_se,
        pooled_z=pooled_z,
        pooled_p=pooled_p,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        model="fixed",
        cochran_q=cochran_q,
        q_p_value=q_p_value,
        i_squared=i_squared,
        tau_squared=tau_sq,
        n_studies=k,
        study_weights=norm_weights,
        study_ids=[s.dataset_id for s in studies],
    )


def random_effects_meta(
    gene: str,
    studies: list[StudyEffect],
) -> MetaResult:
    """DerSimonian-Laird random-effects meta-analysis.

    This is the PRIMARY meta-analysis model for the Riker Engine
    (Blueprint Section 10). It accommodates between-study variability
    by adding tau-squared to each study's variance.

    Parameters
    ----------
    gene : str
        Gene symbol being analyzed.
    studies : list of StudyEffect
        Per-study effect sizes and standard errors.

    Returns
    -------
    MetaResult
        Pooled estimate under the random-effects model.
    """
    _validate_studies(studies)

    effects = np.array([s.log2fc for s in studies])
    ses = np.array([s.se for s in studies])
    k = len(studies)

    # Step 1: Fixed-effects weights for Q calculation
    fe_weights = 1.0 / (ses ** 2)
    fe_total = np.sum(fe_weights)
    fe_pooled = float(np.sum(fe_weights * effects) / fe_total)

    # Step 2: Cochran's Q
    cochran_q = float(np.sum(fe_weights * (effects - fe_pooled) ** 2))
    q_df = k - 1
    q_p_value = float(1.0 - _chi2_cdf(cochran_q, q_df)) if q_df > 0 else 1.0

    # Step 3: Tau-squared (DerSimonian-Laird estimator)
    tau_sq = _estimate_tau_squared(cochran_q, fe_weights, k)

    # Step 4: Random-effects weights
    re_weights = 1.0 / (ses ** 2 + tau_sq)
    re_total = np.sum(re_weights)

    # Step 5: Pooled estimate
    pooled_effect = float(np.sum(re_weights * effects) / re_total)
    pooled_se = float(math.sqrt(1.0 / re_total))

    # Z-test
    pooled_z = pooled_effect / pooled_se
    pooled_p = float(2.0 * norm.sf(abs(pooled_z)))

    # 95% CI
    ci_lower = pooled_effect - 1.96 * pooled_se
    ci_upper = pooled_effect + 1.96 * pooled_se

    # I-squared
    i_squared = max(0.0, (cochran_q - q_df) / cochran_q * 100.0) if cochran_q > 0 else 0.0

    # Normalized weights
    norm_weights = [float(w / re_total) for w in re_weights]

    return MetaResult(
        gene=gene,
        pooled_log2fc=pooled_effect,
        pooled_se=pooled_se,
        pooled_z=pooled_z,
        pooled_p=pooled_p,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        model="random",
        cochran_q=cochran_q,
        q_p_value=q_p_value,
        i_squared=i_squared,
        tau_squared=tau_sq,
        n_studies=k,
        study_weights=norm_weights,
        study_ids=[s.dataset_id for s in studies],
    )


def _estimate_tau_squared(
    cochran_q: float,
    fe_weights: np.ndarray,
    k: int,
) -> float:
    """DerSimonian-Laird estimator for between-study variance.

    tau^2 = max(0, (Q - (k-1)) / (sum(w) - sum(w^2)/sum(w)))

    Parameters
    ----------
    cochran_q : float
        Cochran's Q statistic.
    fe_weights : np.ndarray
        Fixed-effects weights (1/SE^2).
    k : int
        Number of studies.

    Returns
    -------
    float
        Estimated tau-squared (non-negative).
    """
    if k < 2:
        return 0.0

    sum_w = float(np.sum(fe_weights))
    sum_w2 = float(np.sum(fe_weights ** 2))

    denominator = sum_w - sum_w2 / sum_w

    if denominator <= 0:
        return 0.0

    tau_sq = (cochran_q - (k - 1)) / denominator
    return max(0.0, tau_sq)


def _chi2_cdf(x: float, df: int) -> float:
    """Chi-squared CDF using scipy.

    Used for Cochran's Q p-value calculation.
    """
    from scipy.stats import chi2
    if df <= 0 or x < 0:
        return 0.0
    return float(chi2.cdf(x, df))


def run_meta_analysis(
    gene: str,
    studies: list[StudyEffect],
) -> tuple[MetaResult, MetaResult]:
    """Run both fixed and random effects meta-analysis for a gene.

    This is the standard entry point for Phase 6. Returns both models;
    the random-effects result is the primary finding.

    Parameters
    ----------
    gene : str
        Gene symbol.
    studies : list of StudyEffect
        Per-study results. Should include only brain datasets
        (blood datasets are excluded from meta-analysis pooling
        per Blueprint Section 10).

    Returns
    -------
    tuple of (MetaResult, MetaResult)
        (fixed_result, random_result). Use random_result as primary.
    """
    fixed = fixed_effects_meta(gene, studies)
    random = random_effects_meta(gene, studies)
    return fixed, random
```

---

## TESTS: APPEND to `tests/test_stats.py`

**APPEND** the following to the END of `/home/kai001/riker-engine/tests/test_stats.py`.
Do NOT delete or modify any existing test classes (Phase 1 Welch + Phase 2 FDR).

```python


# ===========================================================================
# PHASE 3: META-ANALYSIS TESTS
# ===========================================================================

from riker.stats.meta import (
    MetaResult,
    ScaleCheck,
    StudyEffect,
    check_expression_scale,
    fixed_effects_meta,
    random_effects_meta,
    recover_se,
    run_meta_analysis,
)


# ---------------------------------------------------------------------------
# 11. SE recovery from log2FC and p-value
# ---------------------------------------------------------------------------

class TestSERecovery:
    """Verify SE recovery uses t-distribution and is mathematically correct."""

    def test_known_values(self):
        """Recover SE from a known t-test result.
        If log2FC = 1.0, n_cases=30, n_controls=30, df=58.
        With p=0.01, |t| = t.isf(0.005, 58) ≈ 2.663.
        SE = 1.0 / 2.663 ≈ 0.3755."""
        from scipy.stats import t as t_dist
        expected_t = float(t_dist.isf(0.005, 58))
        expected_se = 1.0 / expected_t

        se = recover_se(log2fc=1.0, p_value=0.01, n_cases=30, n_controls=30)
        assert math.isclose(se, expected_se, rel_tol=1e-6)

    def test_round_trip_with_welch(self):
        """Recover SE from a welch_ttest result, then verify consistency."""
        from riker.stats.welch import welch_ttest
        np.random.seed(42)
        g1 = np.random.normal(5.0, 1.0, 30)
        g2 = np.random.normal(4.0, 1.5, 25)

        result = welch_ttest(g1, g2)
        recovered = recover_se(
            log2fc=result.mean_diff,
            p_value=result.p_value,
            n_cases=result.n1,
            n_controls=result.n2,
        )
        # Recovered SE should be close to original (not exact due to
        # df approximation in recovery: pooled vs Welch-Satterthwaite)
        assert math.isclose(recovered, result.se_diff, rel_tol=0.15)

    def test_large_effect(self):
        """Large fold change should give small relative SE."""
        se = recover_se(log2fc=5.0, p_value=0.001, n_cases=50, n_controls=50)
        assert se > 0
        assert se < 5.0  # SE should be much less than effect

    def test_zero_log2fc(self):
        """Zero fold change with high p-value should still return valid SE."""
        se = recover_se(log2fc=0.0, p_value=0.95, n_cases=20, n_controls=20)
        assert se > 0

    def test_very_small_pvalue(self):
        """Very small p-value should not cause overflow."""
        se = recover_se(log2fc=2.0, p_value=1e-50, n_cases=100, n_controls=100)
        assert se > 0
        assert np.isfinite(se)

    def test_pvalue_one(self):
        """p=1.0 means t=0, should return large SE."""
        se = recover_se(log2fc=0.5, p_value=1.0, n_cases=20, n_controls=20)
        assert se > 1.0  # Large SE to downweight

    def test_invalid_pvalue(self):
        with pytest.raises(ValueError, match="p_value"):
            recover_se(log2fc=1.0, p_value=0.0, n_cases=10, n_controls=10)

    def test_invalid_sample_size(self):
        with pytest.raises(ValueError, match="n >= 2"):
            recover_se(log2fc=1.0, p_value=0.05, n_cases=1, n_controls=10)


# ---------------------------------------------------------------------------
# 12. Expression scale detection (GSE33000 issue)
# ---------------------------------------------------------------------------

class TestScaleDetection:
    """Verify detection of log-ratio vs log2-intensity scale."""

    def test_log2_intensity_not_flagged(self):
        """Typical microarray data (range 3-16) should not be flagged."""
        vals = np.random.uniform(3.0, 16.0, size=1000)
        result = check_expression_scale("GSE5281", vals)
        assert not result.is_log_ratio
        assert result.warning == ""

    def test_log_ratio_flagged(self):
        """Log-ratio data (range -2 to 2) should be flagged."""
        vals = np.random.uniform(-2.0, 2.0, size=1000)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = check_expression_scale("GSE33000", vals)
            assert len(w) == 1
            assert "log-ratio" in str(w[0].message).lower()
        assert result.is_log_ratio
        assert result.expression_range_max < 5.0

    def test_threshold_boundary(self):
        """Value exactly at threshold should not be flagged."""
        vals = np.array([0.0, 5.0])
        result = check_expression_scale("TEST", vals)
        assert not result.is_log_ratio

    def test_empty_data(self):
        vals = np.array([])
        result = check_expression_scale("EMPTY", vals)
        assert result.is_log_ratio  # flagged as problematic


# ---------------------------------------------------------------------------
# 13. Fixed-effects meta-analysis
# ---------------------------------------------------------------------------

class TestFixedEffects:
    """Verify fixed-effects IVW meta-analysis."""

    def _make_studies(self):
        """Create test studies with known properties."""
        return [
            StudyEffect("GSE1", log2fc=-0.5, se=0.15, n_cases=30,
                        n_controls=30, p_value=0.001, tissue="brain"),
            StudyEffect("GSE2", log2fc=-0.4, se=0.20, n_cases=25,
                        n_controls=25, p_value=0.05, tissue="brain"),
            StudyEffect("GSE3", log2fc=-0.6, se=0.18, n_cases=35,
                        n_controls=30, p_value=0.001, tissue="brain"),
            StudyEffect("GSE4", log2fc=-0.3, se=0.25, n_cases=20,
                        n_controls=20, p_value=0.10, tissue="brain"),
        ]

    def test_pooled_effect_direction(self):
        """Pooled effect should be negative (all studies negative)."""
        result = fixed_effects_meta("ATP2B2", self._make_studies())
        assert result.pooled_log2fc < 0

    def test_pooled_effect_manual(self):
        """Verify against manual IVW calculation."""
        studies = self._make_studies()
        effects = np.array([s.log2fc for s in studies])
        ses = np.array([s.se for s in studies])
        weights = 1.0 / ses ** 2
        expected = float(np.sum(weights * effects) / np.sum(weights))

        result = fixed_effects_meta("ATP2B2", studies)
        assert math.isclose(result.pooled_log2fc, expected, rel_tol=1e-10)

    def test_model_label(self):
        result = fixed_effects_meta("TEST", self._make_studies())
        assert result.model == "fixed"

    def test_ci_contains_estimate(self):
        result = fixed_effects_meta("TEST", self._make_studies())
        assert result.ci_lower < result.pooled_log2fc < result.ci_upper

    def test_weights_sum_to_one(self):
        result = fixed_effects_meta("TEST", self._make_studies())
        assert math.isclose(sum(result.study_weights), 1.0, rel_tol=1e-10)

    def test_n_studies(self):
        result = fixed_effects_meta("TEST", self._make_studies())
        assert result.n_studies == 4

    def test_heterogeneity_stats(self):
        result = fixed_effects_meta("TEST", self._make_studies())
        assert result.cochran_q >= 0
        assert 0 <= result.i_squared <= 100
        assert result.tau_squared >= 0
        assert 0 <= result.q_p_value <= 1

    def test_requires_two_studies(self):
        with pytest.raises(ValueError, match="at least 2"):
            fixed_effects_meta("TEST", [self._make_studies()[0]])


# ---------------------------------------------------------------------------
# 14. Random-effects meta-analysis (PRIMARY model)
# ---------------------------------------------------------------------------

class TestRandomEffects:
    """Verify DerSimonian-Laird random-effects meta-analysis."""

    def _make_studies(self):
        return [
            StudyEffect("GSE1", log2fc=-0.5, se=0.15, n_cases=30,
                        n_controls=30, p_value=0.001, tissue="brain"),
            StudyEffect("GSE2", log2fc=-0.4, se=0.20, n_cases=25,
                        n_controls=25, p_value=0.05, tissue="brain"),
            StudyEffect("GSE3", log2fc=-0.6, se=0.18, n_cases=35,
                        n_controls=30, p_value=0.001, tissue="brain"),
            StudyEffect("GSE4", log2fc=-0.3, se=0.25, n_cases=20,
                        n_controls=20, p_value=0.10, tissue="brain"),
        ]

    def test_model_label(self):
        result = random_effects_meta("TEST", self._make_studies())
        assert result.model == "random"

    def test_wider_ci_than_fixed(self):
        """Random effects CI should be >= fixed effects CI width."""
        studies = self._make_studies()
        fixed = fixed_effects_meta("TEST", studies)
        random = random_effects_meta("TEST", studies)
        fixed_width = fixed.ci_upper - fixed.ci_lower
        random_width = random.ci_upper - random.ci_lower
        assert random_width >= fixed_width - 1e-10

    def test_pooled_effect_direction(self):
        result = random_effects_meta("ATP2B2", self._make_studies())
        assert result.pooled_log2fc < 0

    def test_heterogeneous_studies(self):
        """When studies disagree, tau_squared should be positive."""
        studies = [
            StudyEffect("GSE1", log2fc=-1.0, se=0.1, n_cases=50,
                        n_controls=50, p_value=0.001, tissue="brain"),
            StudyEffect("GSE2", log2fc=0.5, se=0.1, n_cases=50,
                        n_controls=50, p_value=0.001, tissue="brain"),
        ]
        result = random_effects_meta("TEST", studies)
        assert result.tau_squared > 0
        assert result.i_squared > 50  # high heterogeneity

    def test_homogeneous_studies(self):
        """When studies agree perfectly, tau_squared should be 0."""
        studies = [
            StudyEffect("GSE1", log2fc=-0.5, se=0.2, n_cases=30,
                        n_controls=30, p_value=0.01, tissue="brain"),
            StudyEffect("GSE2", log2fc=-0.5, se=0.2, n_cases=30,
                        n_controls=30, p_value=0.01, tissue="brain"),
        ]
        result = random_effects_meta("TEST", studies)
        assert result.tau_squared == 0.0
        assert result.i_squared == 0.0

    def test_weights_sum_to_one(self):
        result = random_effects_meta("TEST", self._make_studies())
        assert math.isclose(sum(result.study_weights), 1.0, rel_tol=1e-10)


# ---------------------------------------------------------------------------
# 15. run_meta_analysis entry point
# ---------------------------------------------------------------------------

class TestRunMetaAnalysis:
    """Verify the combined entry point returns both models."""

    def test_returns_both_models(self):
        studies = [
            StudyEffect("GSE1", log2fc=-0.5, se=0.15, n_cases=30,
                        n_controls=30, p_value=0.001, tissue="brain"),
            StudyEffect("GSE2", log2fc=-0.4, se=0.20, n_cases=25,
                        n_controls=25, p_value=0.05, tissue="brain"),
        ]
        fixed, random = run_meta_analysis("ATP2B2", studies)
        assert fixed.model == "fixed"
        assert random.model == "random"
        assert fixed.gene == "ATP2B2"
        assert random.gene == "ATP2B2"

    def test_random_is_primary(self):
        """Random effects should be the second element (primary result)."""
        studies = [
            StudyEffect("GSE1", log2fc=-0.5, se=0.15, n_cases=30,
                        n_controls=30, p_value=0.001, tissue="brain"),
            StudyEffect("GSE2", log2fc=-0.6, se=0.18, n_cases=35,
                        n_controls=30, p_value=0.001, tissue="brain"),
        ]
        _, primary = run_meta_analysis("SEZ6L2", studies)
        assert primary.model == "random"
```

---

## EXECUTION INSTRUCTIONS

After writing meta.py and appending tests, run:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_stats.py -v 2>&1
```

**Expected: ALL tests pass (Phase 1 + Phase 2 + Phase 3).** If any fail, report FULL output — do NOT modify tests or code without reporting first.

Then run this meta-analysis confirmation:

```bash
cd /home/kai001/riker-engine && python -c "
from riker.stats.meta import StudyEffect, run_meta_analysis, recover_se

# Simulate ATP2B2 from the ASD proof-of-concept (4 brain datasets)
studies = [
    StudyEffect('GSE28521', log2fc=-0.45, se=recover_se(-0.45, 0.003, 15, 18),
                n_cases=15, n_controls=18, p_value=0.003, tissue='brain'),
    StudyEffect('GSE38322', log2fc=-0.38, se=recover_se(-0.38, 0.012, 18, 18),
                n_cases=18, n_controls=18, p_value=0.012, tissue='brain'),
    StudyEffect('GSE64018', log2fc=-0.52, se=recover_se(-0.52, 0.008, 12, 12),
                n_cases=12, n_controls=12, p_value=0.008, tissue='brain'),
    StudyEffect('GSE102741', log2fc=-0.30, se=recover_se(-0.30, 0.045, 20, 20),
                n_cases=20, n_controls=20, p_value=0.045, tissue='brain'),
]

fixed, random = run_meta_analysis('ATP2B2', studies)

print('=== ATP2B2 Meta-Analysis (simulated) ===')
print(f'Fixed:  pooled_log2fc={fixed.pooled_log2fc:.4f}, p={fixed.pooled_p:.6f}')
print(f'Random: pooled_log2fc={random.pooled_log2fc:.4f}, p={random.pooled_p:.6f}')
print(f'I-squared: {random.i_squared:.1f}%')
print(f'Tau-squared: {random.tau_squared:.6f}')
print(f'Studies: {random.n_studies}')
print(f'Model: {random.model}')
print()
assert random.model == 'random', 'FAIL: primary model is not random effects!'
assert random.pooled_log2fc < 0, 'FAIL: pooled effect should be negative!'
assert random.pooled_p < 0.05, 'FAIL: pooled effect should be significant!'
print('PASS: meta.py producing correct pooled estimates')
"
```

Report both outputs back to the architect for QA review. Do NOT modify anything if tests fail — report back.
