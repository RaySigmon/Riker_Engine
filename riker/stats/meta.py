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
