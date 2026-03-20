"""
Riker Engine - Phase 6: Effect Size Meta-Analysis.

Computes inverse-variance weighted meta-analysis for each surviving
gene across all discovery datasets. Produces forest plot data and
heterogeneity statistics.

References:
    Blueprint Section 10 (Phase 6: Effect Size Meta-Analysis)
    Blueprint Section 10.1 (IVW fixed/random effects)
    Context Transfer: GSE33000 scale check
"""

import logging
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from riker.stats.meta import (
    run_meta_analysis,
    recover_se,
    check_expression_scale,
    MetaResult,
    StudyEffect,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class GeneEffectSize:
    """Per-dataset effect size for one gene.

    Attributes:
        dataset_id: Dataset identifier.
        log2fc: Log2 fold change (effect size).
        se: Standard error of the effect size.
        p_value: P-value from Welch's t-test.
        n_cases: Number of case samples.
        n_controls: Number of control samples.
        scale_warning: True if dataset has expression scale issues.
    """
    dataset_id: str
    log2fc: float
    se: float
    p_value: float
    n_cases: int
    n_controls: int
    scale_warning: bool


@dataclass(frozen=True)
class GeneMetaResult:
    """Meta-analysis result for one gene.

    Attributes:
        gene: Gene symbol.
        cluster_id: From Phase 4/5 assignment.
        per_dataset: List of GeneEffectSize (forest plot rows).
        fixed_effect: Fixed-effect pooled estimate.
        fixed_se: Fixed-effect standard error.
        fixed_p: Fixed-effect p-value.
        random_effect: Random-effects pooled estimate (DerSimonian-Laird).
        random_se: Random-effects standard error.
        random_p: Random-effects p-value.
        cochran_q: Cochran's Q heterogeneity statistic.
        i_squared: I² heterogeneity percentage.
        tau_squared: Between-study variance (τ²).
        n_datasets: Number of datasets in meta-analysis.
        direction: Overall direction from random-effects estimate.
    """
    gene: str
    cluster_id: int
    per_dataset: list
    fixed_effect: float
    fixed_se: float
    fixed_p: float
    random_effect: float
    random_se: float
    random_p: float
    cochran_q: float
    i_squared: float
    tau_squared: float
    n_datasets: int
    direction: str


@dataclass
class Phase6Result:
    """Complete result of Phase 6 meta-analysis.

    Attributes:
        gene_results: Dict of gene_symbol -> GeneMetaResult.
        n_genes_analyzed: Number of genes with meta-analysis.
        n_significant_random: Genes significant (p<0.05) under random effects.
        n_high_heterogeneity: Genes with I² > 75%.
        scale_warnings: List of dataset IDs with expression scale issues.
    """
    gene_results: dict = field(default_factory=dict)
    n_genes_analyzed: int = 0
    n_significant_random: int = 0
    n_high_heterogeneity: int = 0
    scale_warnings: list = field(default_factory=list)


def compute_gene_meta(
    gene: str,
    cluster_id: int,
    de_results: list,
    dataset_scales: dict[str, bool] | None = None,
) -> GeneMetaResult | None:
    """Compute meta-analysis for one gene.

    Parameters
    ----------
    gene : str
        Gene symbol.
    cluster_id : int
        Cluster assignment.
    de_results : list
        List of GeneDatasetDE from Phase 1 (per-dataset stats).
    dataset_scales : dict or None
        Dataset_id -> True if scale warning applies.

    Returns
    -------
    GeneMetaResult or None
        None if fewer than 2 datasets available.
    """
    if len(de_results) < 2:
        logger.warning(
            f"Gene {gene}: only {len(de_results)} dataset(s), "
            f"need >= 2 for meta-analysis. Skipping."
        )
        return None

    study_effects = []
    per_dataset_info = []

    for de in de_results:
        effect = de.log2fc
        se = de.se

        # If SE looks invalid, try to recover from log2fc and p-value
        if se <= 0 or not np.isfinite(se):
            try:
                se = recover_se(
                    log2fc=de.log2fc,
                    p_value=de.p_value,
                    n_cases=de.n_cases,
                    n_controls=de.n_controls,
                )
            except (ValueError, ZeroDivisionError):
                logger.warning(
                    f"Gene {gene}, dataset {de.dataset_id}: "
                    f"could not recover SE. Skipping dataset."
                )
                continue

        scale_warn = False
        if dataset_scales and de.dataset_id in dataset_scales:
            scale_warn = dataset_scales[de.dataset_id]

        study_effects.append(StudyEffect(
            dataset_id=de.dataset_id,
            log2fc=effect,
            se=se,
            n_cases=de.n_cases,
            n_controls=de.n_controls,
            p_value=de.p_value,
            tissue="brain" # Standard assumption for meta-analysis pooling
        ))
        
        per_dataset_info.append(GeneEffectSize(
            dataset_id=de.dataset_id,
            log2fc=effect,
            se=se,
            p_value=de.p_value,
            n_cases=de.n_cases,
            n_controls=de.n_controls,
            scale_warning=scale_warn,
        ))

    if len(study_effects) < 2:
        return None

    # Run IVW meta-analysis
    fixed, random = run_meta_analysis(
        gene=gene,
        studies=study_effects,
    )

    direction = "down" if random.pooled_log2fc < 0 else "up"

    return GeneMetaResult(
        gene=gene,
        cluster_id=cluster_id,
        per_dataset=per_dataset_info,
        fixed_effect=fixed.pooled_log2fc,
        fixed_se=fixed.pooled_se,
        fixed_p=fixed.pooled_p,
        random_effect=random.pooled_log2fc,
        random_se=random.pooled_se,
        random_p=random.pooled_p,
        cochran_q=random.cochran_q,
        i_squared=random.i_squared,
        tau_squared=random.tau_squared,
        n_datasets=len(study_effects),
        direction=direction,
    )


def run_phase6(
    phase1_result,
    phase5_result,
    dataset_expression_ranges: dict[str, float] | None = None,
) -> Phase6Result:
    """Run Phase 6 meta-analysis for all surviving genes.

    Parameters
    ----------
    phase1_result : Phase1Result
        Contains study_genes with per-dataset DE stats.
    phase5_result : Phase5Result
        Contains gene_verdicts (survived/eliminated) and locked core genes.
    dataset_expression_ranges : dict or None
        Dataset_id -> max expression value (for scale check).
        If None, scale check is skipped.

    Returns
    -------
    Phase6Result
    """
    result = Phase6Result()

    # Check expression scales
    dataset_scales = {}
    if dataset_expression_ranges:
        for ds_id, max_val in dataset_expression_ranges.items():
            # Adaptive wrapper for check_expression_scale(dataset_id, expression_values)
            check = check_expression_scale(ds_id, np.array([max_val]))
            if check.is_log_ratio:
                dataset_scales[ds_id] = True
                result.scale_warnings.append(f"{ds_id}: {check.warning}")
                logger.warning(f"Scale check: {ds_id}: {check.warning}")
            else:
                dataset_scales[ds_id] = False

    # Get surviving genes from Phase 5
    surviving_genes = [
        gene for gene, verdict in phase5_result.gene_verdicts.items()
        if verdict.status == "survived"
    ]

    logger.info(
        f"Phase 6: Computing meta-analysis for {len(surviving_genes)} "
        f"surviving genes."
    )

    for gene in surviving_genes:
        # Get Phase 1 DE results for this gene
        if gene not in phase1_result.study_genes:
            continue

        gene_result = phase1_result.study_genes[gene]
        cluster_id = phase5_result.gene_verdicts[gene].cluster_id

        meta = compute_gene_meta(
            gene, cluster_id, gene_result.de_results,
            dataset_scales=dataset_scales or None,
        )

        if meta is not None:
            result.gene_results[gene] = meta

    result.n_genes_analyzed = len(result.gene_results)
    result.n_significant_random = sum(
        1 for m in result.gene_results.values() if m.random_p < 0.05
    )
    result.n_high_heterogeneity = sum(
        1 for m in result.gene_results.values() if m.i_squared > 75.0
    )

    logger.info(
        f"Phase 6 complete: {result.n_genes_analyzed} genes analyzed, "
        f"{result.n_significant_random} significant (random effects), "
        f"{result.n_high_heterogeneity} high heterogeneity."
    )

    return result
