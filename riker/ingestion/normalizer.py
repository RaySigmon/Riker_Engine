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
    n_negative = int(np.sum(vals < 0))
    pct_negative = n_negative / len(vals) if len(vals) > 0 else 0.0

    if median_val >= LOG2_DETECTION_MEDIAN_THRESHOLD:
        # Raw intensity scale — even if a few negatives from background subtraction
        if n_negative > 0:
            return True, (
                f"Median ({median_val:.2f}) >= {LOG2_DETECTION_MEDIAN_THRESHOLD}. "
                f"Raw intensities with {n_negative} negative values ({pct_negative:.1%}) "
                f"from background subtraction. Will clamp to 0 and apply log2(x+1)."
            )
        return True, (
            f"Median ({median_val:.2f}) >= {LOG2_DETECTION_MEDIAN_THRESHOLD}. "
            f"Data appears to be raw intensities — log2(x+1) recommended."
        )

    if n_negative > 0:
        # Low median + negatives = already log-transformed
        return False, (
            f"Data contains {n_negative} negative values with median {median_val:.2f} "
            f"— already log-transformed."
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
        df_index = expression_data.index
        df_columns = expression_data.columns
    else:
        vals = np.asarray(expression_data, dtype=np.float64)
        is_dataframe = False
        df_index = None
        df_columns = None

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
        if n_negative > 0 and original_median >= LOG2_DETECTION_MEDIAN_THRESHOLD:
            # Background-subtracted raw intensities — clamp negatives, then log2
            logger.warning(
                f"Clamping {n_negative} negative values to 0 "
                f"(background subtraction artifacts in raw intensity data)."
            )
            vals = np.maximum(vals, 0)
            transformed = np.log2(vals + 1)
            was_transformed = True
            reason += f" Clamped {n_negative} negatives to 0."
            logger.info(
                f"Applied log2(x+1) transformation (after clamping). "
                f"Median: {original_median:.2f} -> {float(np.nanmedian(transformed)):.2f}"
            )
        elif n_negative > 0:
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

    # Preserve DataFrame type if input was DataFrame
    if is_dataframe:
        transformed = pd.DataFrame(transformed, index=df_index, columns=df_columns)

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
