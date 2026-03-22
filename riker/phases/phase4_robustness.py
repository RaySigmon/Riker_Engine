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
Riker Engine - Phase 4: Robustness Testing.

Progressive statistical filtering of clusters from Phase 3:
permutation testing, sensitivity analysis, leave-one-dataset-out
stability, and core gene identification.

FDR correction uses the FULL seed gene set as denominator
(Blueprint Section 8.2). This is enforced via apply_fdr_with_scope().

References:
    Blueprint Section 8 (Phase 4: Robustness Testing)
    Blueprint Section 12 (QC Framework)
"""

import logging
import warnings
from dataclasses import dataclass, field

import numpy as np

from riker.stats.permutation import cluster_permutation_test, PermutationResult
from riker.stats.fdr import apply_fdr_with_scope, fdr_survivors, FDRResult

logger = logging.getLogger(__name__)


# Sensitivity levels from Blueprint Section 8.2
SENSITIVITY_LEVELS = {
    1: {"description": "p < 0.05 in 2+ datasets", "threshold": 0.05, "use_fdr": False},
    2: {"description": "p < 0.01 in 2+ datasets", "threshold": 0.01, "use_fdr": False},
    3: {"description": "FDR q < 0.10 in 2+ datasets", "threshold": 0.10, "use_fdr": True},
    4: {"description": "FDR q < 0.05 in 2+ datasets", "threshold": 0.05, "use_fdr": True},
}


@dataclass(frozen=True)
class ClusterSignificance:
    """Permutation test result for one cluster.

    Attributes:
        cluster_id: Cluster label.
        permutation_p: Raw permutation p-value.
        bonferroni_p: Bonferroni-corrected p-value.
        is_significant: True if bonferroni_p < 0.05.
        n_permutations: Number of permutations used.
        observed_stat: Observed test statistic.
        null_mean: Mean of null distribution.
    """
    cluster_id: int
    permutation_p: float
    bonferroni_p: float
    is_significant: bool
    n_permutations: int
    observed_stat: float
    null_mean: float


@dataclass(frozen=True)
class SensitivityResult:
    """Sensitivity analysis result for one cluster.

    Attributes:
        cluster_id: Cluster label.
        level_survivors: Dict of level (1-4) -> list of surviving gene symbols.
        dissolution_point: Level at which < 3 genes remain (or None if all survive).
        core_genes: Genes surviving Level 2 (the pre-replication set).
    """
    cluster_id: int
    level_survivors: dict
    dissolution_point: int | None
    core_genes: list


@dataclass(frozen=True)
class LOOResult:
    """Leave-one-dataset-out stability result for one cluster.

    Attributes:
        cluster_id: Cluster label.
        is_stable: True if >= 80% genes survive removal of any single dataset.
        per_dataset_retention: Dict of dataset_id -> fraction of genes retained.
        min_retention: Minimum retention fraction across all datasets.
        fragile_datasets: Datasets whose removal causes > 20% gene loss.
    """
    cluster_id: int
    is_stable: bool
    per_dataset_retention: dict
    min_retention: float
    fragile_datasets: list


@dataclass(frozen=True)
class CoreGene:
    """A single core gene identified in Phase 4.

    Attributes:
        gene: Gene symbol.
        cluster_id: Which cluster this gene belongs to.
        max_level_survived: Highest sensitivity level survived (1-4).
        per_dataset_pvalues: Dict of dataset_id -> p-value.
        per_dataset_log2fc: Dict of dataset_id -> log2FC.
        mean_log2fc: Mean log2FC across all datasets.
        direction: Predominant direction ('up' or 'down').
    """
    gene: str
    cluster_id: int
    max_level_survived: int
    per_dataset_pvalues: dict
    per_dataset_log2fc: dict
    mean_log2fc: float
    direction: str


@dataclass
class Phase4Result:
    """Complete result of Phase 4 robustness testing.

    Attributes:
        cluster_significance: Dict of cluster_id -> ClusterSignificance.
        sensitivity: Dict of cluster_id -> SensitivityResult.
        loo_stability: Dict of cluster_id -> LOOResult.
        core_genes: Dict of gene_symbol -> CoreGene.
        n_core_genes: Total number of core genes across all clusters.
        n_clusters_significant: Clusters passing Bonferroni permutation test.
        n_clusters_stable: Clusters passing LOO stability.
    """
    cluster_significance: dict = field(default_factory=dict)
    sensitivity: dict = field(default_factory=dict)
    loo_stability: dict = field(default_factory=dict)
    core_genes: dict = field(default_factory=dict)
    n_core_genes: int = 0
    n_clusters_significant: int = 0
    n_clusters_stable: int = 0


def evaluate_cluster_significance(
    cluster_info: dict,
    study_genes: dict,
    all_study_gene_symbols: list[str],
    n_permutations: int = 10000,
    seed: int = 42,
) -> dict[int, ClusterSignificance]:
    """Run permutation tests for all clusters with Bonferroni correction.

    Parameters
    ----------
    cluster_info : dict
        From Phase3Result.cluster_info (cluster_id -> ClusterInfo).
    study_genes : dict
        From Phase1Result.study_genes (gene -> GeneResult).
    all_study_gene_symbols : list
        All study gene symbols (the permutation pool).
    n_permutations : int
        Number of permutations per cluster.
    seed : int
        Base random seed (incremented per cluster for independence).

    Returns
    -------
    dict of cluster_id -> ClusterSignificance
    """
    n_clusters = len(cluster_info)
    if n_clusters == 0:
        return {}

    # Build gene -> mean |log2FC| map for permutation stat
    gene_log2fc = {}
    for gene, result in study_genes.items():
        if result.de_results:
            gene_log2fc[gene] = float(
                np.mean([d.log2fc for d in result.de_results])
            )
        else:
            gene_log2fc[gene] = 0.0

    results = {}
    for cid, info in cluster_info.items():
        perm_result = cluster_permutation_test(
            cluster_genes=info.gene_symbols,
            all_study_genes=all_study_gene_symbols,
            gene_log2fc=gene_log2fc,
            n_permutations=n_permutations,
            seed=seed + cid,  # different seed per cluster
        )

        bonf_p = min(perm_result.p_value * n_clusters, 1.0)

        results[cid] = ClusterSignificance(
            cluster_id=cid,
            permutation_p=perm_result.p_value,
            bonferroni_p=bonf_p,
            is_significant=(bonf_p < 0.05),
            n_permutations=n_permutations,
            observed_stat=perm_result.observed,
            null_mean=perm_result.null_mean,
        )

        logger.info(
            f"Cluster {cid}: permutation p={perm_result.p_value:.4f}, "
            f"Bonferroni p={bonf_p:.4f} "
            f"({'SIGNIFICANT' if bonf_p < 0.05 else 'not significant'})"
        )

    return results


def sensitivity_analysis(
    cluster_info: dict,
    study_genes: dict,
    seed_gene_count: int,
    min_datasets: int = 2,
) -> dict[int, SensitivityResult]:
    """Progressive threshold sensitivity analysis for all clusters.

    Tests genes at four progressively stricter levels (Blueprint Section 8.2).
    FDR correction at Levels 3-4 uses the FULL seed gene count.

    Parameters
    ----------
    cluster_info : dict
        From Phase3Result.cluster_info.
    study_genes : dict
        From Phase1Result.study_genes.
    seed_gene_count : int
        Total seed genes for FDR scope enforcement.
    min_datasets : int
        Minimum datasets for a gene to count as significant.

    Returns
    -------
    dict of cluster_id -> SensitivityResult
    """
    # Pre-compute per-gene, per-dataset p-values and FDR q-values
    # Collect all p-values across all genes for FDR computation
    all_min_pvalues = {}
    gene_dataset_pvals = {}

    for gene, result in study_genes.items():
        per_ds = {}
        for de in result.de_results:
            per_ds[de.dataset_id] = de.p_value
        gene_dataset_pvals[gene] = per_ds

        if per_ds:
            all_min_pvalues[gene] = min(per_ds.values())

    # Compute FDR q-values over full seed set (for Levels 3-4)
    fdr_result = apply_fdr_with_scope(
        all_min_pvalues,
        seed_gene_count=seed_gene_count,
        scope="full_seed_set",
        threshold=0.10,
    )

    results = {}
    for cid, info in cluster_info.items():
        level_survivors = {}
        cluster_genes = set(info.gene_symbols)

        for level, config in SENSITIVITY_LEVELS.items():
            threshold = config["threshold"]
            use_fdr = config["use_fdr"]
            survivors = []

            for gene in info.gene_symbols:
                if gene not in study_genes:
                    continue

                if use_fdr:
                    # Use FDR q-value
                    q_val = fdr_result.q_values.get(gene, 1.0)
                    if q_val < threshold:
                        # Also check min_datasets criterion on raw p-values
                        n_sig = sum(
                            1 for p in gene_dataset_pvals.get(gene, {}).values()
                            if p < 0.05  # raw p for dataset counting
                        )
                        if n_sig >= min_datasets:
                            survivors.append(gene)
                else:
                    # Use raw p-value threshold
                    n_sig = sum(
                        1 for p in gene_dataset_pvals.get(gene, {}).values()
                        if p < threshold
                    )
                    if n_sig >= min_datasets:
                        survivors.append(gene)

            level_survivors[level] = survivors

        # Find dissolution point
        dissolution = None
        for level in sorted(SENSITIVITY_LEVELS.keys()):
            if len(level_survivors[level]) < 3:
                dissolution = level
                break

        # Core genes = Level 2 survivors
        core = level_survivors.get(2, [])

        results[cid] = SensitivityResult(
            cluster_id=cid,
            level_survivors=level_survivors,
            dissolution_point=dissolution,
            core_genes=core,
        )

        logger.info(
            f"Cluster {cid} sensitivity: "
            + ", ".join(
                f"L{lv}={len(level_survivors[lv])}"
                for lv in sorted(level_survivors)
            )
            + f" | dissolution={dissolution} | core={len(core)}"
        )

    return results


def loo_stability(
    cluster_info: dict,
    study_genes: dict,
    dataset_ids: list[str],
    p_threshold: float = 0.05,
    min_datasets: int = 2,
    stability_threshold: float = 0.80,
) -> dict[int, LOOResult]:
    """Leave-one-dataset-out stability analysis.

    For each cluster, removes each dataset in turn and checks what
    fraction of genes still meet the retention criterion.

    Parameters
    ----------
    cluster_info : dict
        From Phase3Result.cluster_info.
    study_genes : dict
        From Phase1Result.study_genes.
    dataset_ids : list
        All discovery dataset IDs.
    p_threshold : float
        P-value threshold for retention (default 0.05).
    min_datasets : int
        Minimum datasets after removal (default 2).
    stability_threshold : float
        Minimum retention fraction for stability (default 0.80).

    Returns
    -------
    dict of cluster_id -> LOOResult
    """
    results = {}

    for cid, info in cluster_info.items():
        per_ds_retention = {}
        cluster_genes = info.gene_symbols

        for removed_ds in dataset_ids:
            remaining_ds = [d for d in dataset_ids if d != removed_ds]

            if len(remaining_ds) < min_datasets:
                # Can't evaluate with fewer datasets than min required
                per_ds_retention[removed_ds] = 1.0
                continue

            n_surviving = 0
            for gene in cluster_genes:
                if gene not in study_genes:
                    continue
                gene_result = study_genes[gene]
                # Count significant datasets excluding the removed one
                n_sig = sum(
                    1 for de in gene_result.de_results
                    if de.dataset_id != removed_ds and de.p_value < p_threshold
                )
                if n_sig >= min_datasets:
                    n_surviving += 1

            retention = n_surviving / len(cluster_genes) if cluster_genes else 0.0
            per_ds_retention[removed_ds] = retention

        min_ret = min(per_ds_retention.values()) if per_ds_retention else 0.0
        is_stable = min_ret >= stability_threshold
        fragile = [ds for ds, ret in per_ds_retention.items()
                    if ret < stability_threshold]

        results[cid] = LOOResult(
            cluster_id=cid,
            is_stable=is_stable,
            per_dataset_retention=per_ds_retention,
            min_retention=min_ret,
            fragile_datasets=fragile,
        )

        logger.info(
            f"Cluster {cid} LOO: min_retention={min_ret:.2f} "
            f"({'STABLE' if is_stable else 'FRAGILE'})"
        )

    return results


def identify_core_genes(
    sensitivity_results: dict[int, SensitivityResult],
    study_genes: dict,
    min_genes_per_cluster: int = 3,
) -> dict[str, CoreGene]:
    """Identify core genes from sensitivity analysis.

    Core genes are those surviving Level 2 (p < 0.01 in 2+ datasets)
    in clusters with at least min_genes_per_cluster Level 2 survivors.

    Parameters
    ----------
    sensitivity_results : dict
        From sensitivity_analysis().
    study_genes : dict
        From Phase1Result.study_genes.
    min_genes_per_cluster : int
        Minimum Level 2 survivors for a cluster to contribute core genes.

    Returns
    -------
    dict of gene_symbol -> CoreGene
    """
    core_genes = {}

    for cid, sens in sensitivity_results.items():
        level2_genes = sens.core_genes

        if len(level2_genes) < min_genes_per_cluster:
            logger.info(
                f"Cluster {cid}: {len(level2_genes)} Level 2 survivors "
                f"< {min_genes_per_cluster} minimum. No core genes."
            )
            continue

        for gene in level2_genes:
            if gene not in study_genes:
                continue

            gene_result = study_genes[gene]

            # Determine max level survived
            max_level = 1
            for level in [2, 3, 4]:
                if gene in sens.level_survivors.get(level, []):
                    max_level = level

            # Collect per-dataset stats
            per_ds_p = {}
            per_ds_fc = {}
            for de in gene_result.de_results:
                per_ds_p[de.dataset_id] = de.p_value
                per_ds_fc[de.dataset_id] = de.log2fc

            mean_fc = gene_result.mean_log2fc
            direction = "down" if mean_fc < 0 else "up"

            core_genes[gene] = CoreGene(
                gene=gene,
                cluster_id=cid,
                max_level_survived=max_level,
                per_dataset_pvalues=per_ds_p,
                per_dataset_log2fc=per_ds_fc,
                mean_log2fc=mean_fc,
                direction=direction,
            )

    logger.info(f"Identified {len(core_genes)} core genes across all clusters.")
    return core_genes


def run_phase4(
    phase1_result,
    phase3_result,
    seed_gene_count: int,
    dataset_ids: list[str],
    n_permutations: int = 10000,
    permutation_seed: int = 42,
) -> Phase4Result:
    """Run complete Phase 4 robustness testing.

    Parameters
    ----------
    phase1_result : Phase1Result
        From Phase 1.
    phase3_result : Phase3Result
        From Phase 3.
    seed_gene_count : int
        Total seed genes for FDR scope enforcement.
    dataset_ids : list
        Discovery dataset IDs for LOO analysis.
    n_permutations : int
        Permutation count per cluster.
    permutation_seed : int
        Base seed for permutation tests.

    Returns
    -------
    Phase4Result
    """
    study_genes = phase1_result.study_genes
    all_symbols = list(study_genes.keys())
    cluster_info = phase3_result.cluster_info

    # 1. Permutation significance
    significance = evaluate_cluster_significance(
        cluster_info, study_genes, all_symbols,
        n_permutations=n_permutations,
        seed=permutation_seed,
    )

    # 2. Sensitivity analysis
    sens = sensitivity_analysis(
        cluster_info, study_genes,
        seed_gene_count=seed_gene_count,
    )

    # 3. LOO stability
    loo = loo_stability(
        cluster_info, study_genes, dataset_ids,
    )

    # 4. Core gene identification
    core = identify_core_genes(sens, study_genes)

    n_sig = sum(1 for s in significance.values() if s.is_significant)
    n_stable = sum(1 for l in loo.values() if l.is_stable)

    result = Phase4Result(
        cluster_significance=significance,
        sensitivity=sens,
        loo_stability=loo,
        core_genes=core,
        n_core_genes=len(core),
        n_clusters_significant=n_sig,
        n_clusters_stable=n_stable,
    )

    logger.info(
        f"Phase 4 complete: {n_sig} significant clusters, "
        f"{n_stable} LOO-stable clusters, {len(core)} core genes."
    )

    return result
