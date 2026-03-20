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
