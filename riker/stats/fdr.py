# Riker Engine - Condition-Agnostic Transcriptomics Pipeline
# Copyright (C) 2024-2026 Ray Sigmon
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
Riker Engine - Benjamini-Hochberg FDR with enforced scope control.

The FDR scope enforcement is the single most critical statistical guardrail
in the engine. During the AD stress test, reducing FDR scope from the full
seed gene set (800 genes) to the study set only (420 genes) inflated results
by reducing the multiple-testing burden. This module enforces full-seed-set
FDR as the default and flags any reduced-scope analysis as exploratory.

References:
    Blueprint Section 8.2 (FDR scope enforcement)
    Blueprint Section 15 (FDR sensitivity to seed set size)
    Context Transfer: "FDR Scope Is Critical"
"""

import warnings
from dataclasses import dataclass, field
from typing import Literal

import numpy as np


@dataclass(frozen=True)
class FDRResult:
    """Result of FDR correction with scope metadata.

    Attributes:
        q_values: Mapping of gene_symbol -> BH-adjusted q-value.
        n_tested: The n used for BH denominator (seed count or study count).
        n_significant: Count of genes with q < threshold.
        scope: 'full_seed_set' or 'exploratory_study_set_only'.
        threshold: The q-value threshold used for n_significant count.
    """
    q_values: dict
    n_tested: int
    n_significant: int
    scope: str
    threshold: float


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Standard Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    p_values : array-like of float
        Raw p-values. Must all be in [0, 1].

    Returns
    -------
    np.ndarray
        BH-adjusted q-values in the SAME ORDER as input p-values.

    Algorithm
    ---------
    1. Sort p-values ascending.
    2. For rank i (1-based) out of n total: q_i = p_i * n / i
    3. Enforce monotonicity: cumulative minimum from bottom up.
    4. Clip to [0, 1].
    5. Return q-values in original input order.
    """
    pvals = np.asarray(p_values, dtype=np.float64)

    if len(pvals) == 0:
        return np.array([], dtype=np.float64)

    n = len(pvals)

    # Sort ascending, track original indices
    sorted_indices = np.argsort(pvals)
    sorted_pvals = pvals[sorted_indices]

    # BH adjustment: q_i = p_i * n / rank (1-based)
    ranks = np.arange(1, n + 1, dtype=np.float64)
    adjusted = sorted_pvals * n / ranks

    # Enforce monotonicity: cumulative minimum from the right
    # (ensures q[i] <= q[i+1] in sorted order)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]

    # Clip to [0, 1]
    adjusted = np.clip(adjusted, 0.0, 1.0)

    # Restore original order
    q_values = np.empty(n, dtype=np.float64)
    q_values[sorted_indices] = adjusted

    return q_values


def _validate_p_values(p_values: dict) -> None:
    """Validate all p-values are in [0, 1]."""
    for gene, pval in p_values.items():
        if not (0.0 <= pval <= 1.0):
            raise ValueError(
                f"p-value for gene '{gene}' is {pval}, must be in [0, 1]."
            )


def apply_fdr_with_scope(
    p_values: dict,
    seed_gene_count: int,
    scope: Literal["full_seed_set", "study_set_only"] = "full_seed_set",
    threshold: float = 0.10,
) -> FDRResult:
    """Apply BH FDR correction with enforced scope control.

    This is the main entry point for FDR in the Riker Engine. Downstream
    phases (Phase 4 sensitivity analysis) call this function.

    Parameters
    ----------
    p_values : dict
        Mapping of gene_symbol -> p_value for genes being tested.
        These are the STUDY genes (genes that passed cross-referencing).
    seed_gene_count : int
        Total number of seed genes in the disease gene database.
        For ASD: 1267 (SFARI). For AD: 800 (ADSP+Agora+GWAS).
        This is the denominator for full-seed-set FDR.
    scope : str
        'full_seed_set' (default, enforced): Uses seed_gene_count as n.
        'study_set_only' (exploratory): Uses len(p_values) as n.
        WARNING: study_set_only inflates results by reducing the
        multiple-testing burden. It is logged as exploratory.
    threshold : float
        Q-value threshold for counting significant genes (default 0.10).

    Returns
    -------
    FDRResult
        Contains q_values dict, n_tested, n_significant, scope label,
        and threshold.

    Raises
    ------
    ValueError
        If p-values are outside [0, 1], seed_gene_count < len(p_values),
        or seed_gene_count <= 0.
    """
    # --- Input validation ---
    if seed_gene_count <= 0:
        raise ValueError(
            f"seed_gene_count must be > 0, got {seed_gene_count}."
        )
    if seed_gene_count < len(p_values):
        raise ValueError(
            f"seed_gene_count ({seed_gene_count}) must be >= number of "
            f"study genes ({len(p_values)}). The seed set is the superset."
        )
    if not (0.0 < threshold <= 1.0):
        raise ValueError(
            f"threshold must be in (0, 1], got {threshold}."
        )
    _validate_p_values(p_values)

    if len(p_values) == 0:
        return FDRResult(
            q_values={},
            n_tested=seed_gene_count,
            n_significant=0,
            scope="full_seed_set",
            threshold=threshold,
        )

    # --- Determine scope ---
    if scope == "study_set_only":
        warnings.warn(
            "FDR scope set to 'study_set_only'. This reduces the "
            "multiple-testing burden and INFLATES significance. "
            f"Using n={len(p_values)} instead of n={seed_gene_count}. "
            "Results will be flagged as EXPLORATORY. "
            "See Blueprint Section 8.2 and Context Transfer: "
            "'FDR Scope Is Critical'.",
            UserWarning,
            stacklevel=2,
        )
        scope_label = "exploratory_study_set_only"
        n_for_bh = len(p_values)
    elif scope == "full_seed_set":
        scope_label = "full_seed_set"
        n_for_bh = seed_gene_count
    else:
        raise ValueError(
            f"scope must be 'full_seed_set' or 'study_set_only', "
            f"got '{scope}'."
        )

    # --- Build the padded p-value array ---
    genes = list(p_values.keys())
    pvals = np.array([p_values[g] for g in genes], dtype=np.float64)

    if scope == "full_seed_set" and n_for_bh > len(pvals):
        # Pad with p=1.0 for untested seed genes
        n_pad = n_for_bh - len(pvals)
        padded_pvals = np.concatenate([pvals, np.ones(n_pad)])
    else:
        padded_pvals = pvals

    # --- Run BH on the full (possibly padded) array ---
    all_q_values = benjamini_hochberg(padded_pvals)

    # Extract q-values for the real genes only (first len(genes) entries)
    gene_q_values = {
        gene: float(all_q_values[i]) for i, gene in enumerate(genes)
    }

    # --- Count survivors ---
    n_sig = sum(1 for q in gene_q_values.values() if q < threshold)

    return FDRResult(
        q_values=gene_q_values,
        n_tested=n_for_bh,
        n_significant=n_sig,
        scope=scope_label,
        threshold=threshold,
    )


def fdr_survivors(
    fdr_result: FDRResult,
    threshold: float | None = None,
) -> list[tuple[str, float]]:
    """Extract genes surviving FDR correction.

    Parameters
    ----------
    fdr_result : FDRResult
        Output from apply_fdr_with_scope().
    threshold : float, optional
        Override the threshold from the FDRResult. If None, uses
        the threshold stored in the result.

    Returns
    -------
    list of (gene_symbol, q_value) tuples
        Sorted by q-value ascending.
    """
    thresh = threshold if threshold is not None else fdr_result.threshold
    survivors = [
        (gene, q) for gene, q in fdr_result.q_values.items() if q < thresh
    ]
    survivors.sort(key=lambda x: x[1])
    return survivors
