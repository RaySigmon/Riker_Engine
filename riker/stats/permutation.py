"""
Riker Engine - Permutation testing framework.

Tests whether an observed statistic for a gene group is more extreme
than expected by chance, using random gene groups of the same size
drawn from the full study gene set.

The default test statistic is mean absolute log2FC, used in Phase 3
to assess cluster significance after consensus clustering.

References:
    Blueprint Section 7 (Phase 3 — permutation testing)
    Blueprint Section 8.1 (Bonferroni correction for multiple clusters)
    Blueprint Section 12 (QC Framework)
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass(frozen=True)
class PermutationResult:
    """Result of a permutation test.

    Attributes:
        observed: The observed test statistic for the real gene group.
        p_value: Permutation p-value (proportion of null >= observed).
        n_permutations: Number of permutations performed.
        n_extreme: Count of null statistics >= observed.
        null_mean: Mean of the null distribution.
        null_std: Standard deviation of the null distribution.
        null_min: Minimum of the null distribution.
        null_max: Maximum of the null distribution.
    """
    observed: float
    p_value: float
    n_permutations: int
    n_extreme: int
    null_mean: float
    null_std: float
    null_min: float
    null_max: float


def permutation_test(
    observed_genes: list[str],
    all_genes: list[str],
    gene_values: dict[str, float],
    stat_func: Callable[[list[float]], float],
    n_permutations: int = 10_000,
    seed: int | None = 42,
    alternative: str = "greater",
) -> PermutationResult:
    """Generic permutation test for a gene group.

    Draws random gene groups of the same size as observed_genes from
    all_genes, computes the test statistic for each, and returns
    the permutation p-value.

    Parameters
    ----------
    observed_genes : list of str
        Gene symbols in the group being tested (e.g., a cluster).
    all_genes : list of str
        All gene symbols in the study set (the sampling pool).
    gene_values : dict
        Mapping of gene_symbol -> numeric value used by stat_func
        (e.g., gene -> mean log2FC across datasets).
    stat_func : callable
        Function that takes a list of float values and returns a
        single float statistic. Applied to both the observed group
        and each permuted group.
    n_permutations : int
        Number of random permutations (default 10,000).
    seed : int or None
        Random seed for reproducibility. None for no seeding.
    alternative : str
        'greater' (default): p = P(null >= observed).
        'less': p = P(null <= observed).
        'two-sided': p = P(|null - null_mean| >= |observed - null_mean|).

    Returns
    -------
    PermutationResult
        Contains observed stat, p-value, null distribution summary.

    Raises
    ------
    ValueError
        If observed_genes is empty, larger than all_genes, or if
        gene_values is missing entries for any gene in all_genes.
    """
    # --- Input validation ---
    if len(observed_genes) == 0:
        raise ValueError("observed_genes cannot be empty.")
    if len(observed_genes) > len(all_genes):
        raise ValueError(
            f"observed_genes ({len(observed_genes)}) cannot be larger "
            f"than all_genes ({len(all_genes)})."
        )
    if n_permutations < 1:
        raise ValueError(
            f"n_permutations must be >= 1, got {n_permutations}."
        )
    if alternative not in ("greater", "less", "two-sided"):
        raise ValueError(
            f"alternative must be 'greater', 'less', or 'two-sided', "
            f"got '{alternative}'."
        )

    # Check all genes have values
    missing = [g for g in all_genes if g not in gene_values]
    if missing:
        raise ValueError(
            f"{len(missing)} genes in all_genes missing from gene_values. "
            f"First 5: {missing[:5]}"
        )

    # --- Compute observed statistic ---
    observed_values = [gene_values[g] for g in observed_genes]
    observed_stat = stat_func(observed_values)

    # --- Build null distribution ---
    rng = np.random.RandomState(seed)
    all_genes_arr = np.array(all_genes)
    k = len(observed_genes)
    null_stats = np.empty(n_permutations, dtype=np.float64)

    for i in range(n_permutations):
        perm_indices = rng.choice(len(all_genes_arr), size=k, replace=False)
        perm_genes = all_genes_arr[perm_indices]
        perm_values = [gene_values[g] for g in perm_genes]
        null_stats[i] = stat_func(perm_values)

    # --- Compute p-value ---
    # Conservative formula: (n_extreme + 1) / (n_permutations + 1)
    if alternative == "greater":
        n_extreme = int(np.sum(null_stats >= observed_stat))
    elif alternative == "less":
        n_extreme = int(np.sum(null_stats <= observed_stat))
    else:  # two-sided
        null_mean = float(np.mean(null_stats))
        obs_dev = abs(observed_stat - null_mean)
        n_extreme = int(np.sum(np.abs(null_stats - null_mean) >= obs_dev))

    p_value = (n_extreme + 1) / (n_permutations + 1)

    return PermutationResult(
        observed=observed_stat,
        p_value=p_value,
        n_permutations=n_permutations,
        n_extreme=n_extreme,
        null_mean=float(np.mean(null_stats)),
        null_std=float(np.std(null_stats)),
        null_min=float(np.min(null_stats)),
        null_max=float(np.max(null_stats)),
    )


def mean_abs_log2fc(values: list[float]) -> float:
    """Compute mean absolute log2FC for a gene group.

    This is the default test statistic for cluster permutation testing
    in Phase 3 (Blueprint Section 7).

    Parameters
    ----------
    values : list of float
        log2FC values for genes in the group.

    Returns
    -------
    float
        Mean of absolute values.
    """
    if len(values) == 0:
        return 0.0
    return float(np.mean(np.abs(values)))


def cluster_permutation_test(
    cluster_genes: list[str],
    all_study_genes: list[str],
    gene_log2fc: dict[str, float],
    n_permutations: int = 10_000,
    seed: int | None = 42,
) -> PermutationResult:
    """Convenience function: permutation test for cluster mean |log2FC|.

    Tests whether a cluster's mean absolute fold change is greater
    than expected by chance from random gene groups of the same size.

    This is the standard test applied in Phase 3 after consensus
    clustering (Blueprint Section 7).

    Parameters
    ----------
    cluster_genes : list of str
        Gene symbols in the cluster.
    all_study_genes : list of str
        All study genes (the pool for random sampling).
    gene_log2fc : dict
        Mapping of gene_symbol -> average log2FC across datasets.
    n_permutations : int
        Number of permutations (default 10,000).
    seed : int or None
        Random seed for reproducibility.

    Returns
    -------
    PermutationResult
        Result with p-value for the cluster.
    """
    return permutation_test(
        observed_genes=cluster_genes,
        all_genes=all_study_genes,
        gene_values=gene_log2fc,
        stat_func=mean_abs_log2fc,
        n_permutations=n_permutations,
        seed=seed,
        alternative="greater",
    )
