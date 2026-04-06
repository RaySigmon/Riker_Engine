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
Riker Engine - Phase 5: Independent Replication.

Tests core genes from Phase 4 in held-out replication datasets.
The core gene list is PRE-SPECIFIED and LOCKED before replication
data is accessed. No genes may be added or removed based on results.

The elimination protocol removes genes showing significant opposite
direction in same-tissue replication datasets. Cross-tissue non-replication
(e.g., brain signal absent in blood) is expected and does not trigger
elimination.

References:
    Blueprint Section 9 (Phase 5: Independent Replication)
    Blueprint Section 12 (QC Framework — pre-specification, elimination)
"""

import logging
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from riker.stats.welch import welch_ttest

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ReplicationResult:
    """Replication test result for one gene in one dataset.

    Attributes:
        gene: Gene symbol.
        dataset_id: Replication dataset identifier.
        tissue: 'brain' or 'blood'.
        log2fc: Log2 fold change in replication dataset.
        p_value: P-value from Welch's t-test.
        direction: 'up' or 'down'.
        discovery_direction: Direction from discovery (Phase 1).
        is_concordant: True if replication direction matches discovery.
        is_significant: True if p < 0.05.
        n_cases: Number of case samples.
        n_controls: Number of control samples.
    """
    gene: str
    dataset_id: str
    tissue: str
    log2fc: float
    p_value: float
    direction: str
    discovery_direction: str
    is_concordant: bool
    is_significant: bool
    n_cases: int
    n_controls: int


@dataclass(frozen=True)
class GeneVerdict:
    """Elimination verdict for one core gene.

    Attributes:
        gene: Gene symbol.
        cluster_id: From Phase 4 core gene assignment.
        status: 'survived', 'eliminated', or 'insufficient_data'.
        reason: Human-readable explanation.
        replication_results: List of ReplicationResult across all datasets.
        n_same_tissue_concordant: Same-tissue datasets with concordant direction.
        n_same_tissue_discordant: Same-tissue datasets with significant opposite direction.
        n_cross_tissue_tested: Cross-tissue datasets tested.
        n_cross_tissue_concordant: Cross-tissue datasets with concordant direction.
        discovery_direction: Direction from discovery.
    """
    gene: str
    cluster_id: int
    status: str
    reason: str
    replication_results: list
    n_same_tissue_concordant: int
    n_same_tissue_discordant: int
    n_cross_tissue_tested: int
    n_cross_tissue_concordant: int
    discovery_direction: str


@dataclass(frozen=True)
class ClusterVerdict:
    """Replication verdict for one cluster.

    Attributes:
        cluster_id: Cluster label.
        verdict: 'replicated', 'partially_replicated', 'brain_specific', 'failed'.
        n_core_genes: Total core genes in this cluster.
        n_survived: Core genes that survived elimination.
        n_eliminated: Core genes eliminated.
        survived_genes: List of surviving gene symbols.
        eliminated_genes: List of eliminated gene symbols.
    """
    cluster_id: int
    verdict: str
    n_core_genes: int
    n_survived: int
    n_eliminated: int
    survived_genes: list
    eliminated_genes: list


@dataclass
class Phase5Result:
    """Complete result of Phase 5 replication.

    Attributes:
        gene_verdicts: Dict of gene_symbol -> GeneVerdict.
        cluster_verdicts: Dict of cluster_id -> ClusterVerdict.
        n_survived: Total genes surviving elimination.
        n_eliminated: Total genes eliminated.
        n_insufficient: Genes with insufficient replication data.
        locked_core_genes: The pre-specified core gene list (immutable record).
    """
    gene_verdicts: dict = field(default_factory=dict)
    cluster_verdicts: dict = field(default_factory=dict)
    n_survived: int = 0
    n_eliminated: int = 0
    n_insufficient: int = 0
    locked_core_genes: list = field(default_factory=list)


def replicate_gene(
    gene: str,
    discovery_direction: str,
    expression: pd.DataFrame,
    phenotype: dict[str, str],
    dataset_id: str,
    tissue: str = "brain",
) -> ReplicationResult | None:
    """Test one gene in one replication dataset.

    Parameters
    ----------
    gene : str
        Gene symbol to test.
    discovery_direction : str
        Direction from discovery ('up' or 'down').
    expression : DataFrame
        Gene-level expression (genes × samples).
    phenotype : dict
        Sample_id -> 'case' or 'control'.
    dataset_id : str
        Replication dataset identifier.
    tissue : str
        'brain' or 'blood'.

    Returns
    -------
    ReplicationResult or None
        None if gene is not in the dataset or insufficient samples.
    """
    if gene not in expression.index:
        return None

    case_samples = [s for s, g in phenotype.items()
                    if g == "case" and s in expression.columns]
    ctrl_samples = [s for s, g in phenotype.items()
                    if g == "control" and s in expression.columns]

    if len(case_samples) < 2 or len(ctrl_samples) < 2:
        return None

    case_vals = expression.loc[gene, case_samples].values.astype(np.float64)
    ctrl_vals = expression.loc[gene, ctrl_samples].values.astype(np.float64)

    case_vals = case_vals[np.isfinite(case_vals)]
    ctrl_vals = ctrl_vals[np.isfinite(ctrl_vals)]

    if len(case_vals) < 2 or len(ctrl_vals) < 2:
        return None

    try:
        result = welch_ttest(case_vals, ctrl_vals)
    except ValueError:
        return None

    direction = "up" if result.mean_diff > 0 else "down"
    is_concordant = (direction == discovery_direction)

    return ReplicationResult(
        gene=gene,
        dataset_id=dataset_id,
        tissue=tissue,
        log2fc=result.mean_diff,
        p_value=result.p_value,
        direction=direction,
        discovery_direction=discovery_direction,
        is_concordant=is_concordant,
        is_significant=(result.p_value < 0.05),
        n_cases=result.n1,
        n_controls=result.n2,
    )


def run_elimination_protocol(
    core_genes: dict,
    replication_datasets: dict[str, pd.DataFrame],
    replication_phenotypes: dict[str, dict[str, str]],
    dataset_tissues: dict[str, str],
    discovery_tissues: set[str] | None = None,
    p_threshold: float = 0.05,
) -> dict[str, GeneVerdict]:
    """Apply elimination protocol to all core genes.

    Genes are eliminated when they show significant opposite-direction
    effects in same-tissue replication datasets. Cross-tissue
    non-replication is tolerated (e.g., brain signal absent in blood).

    Parameters
    ----------
    core_genes : dict
        From Phase4Result.core_genes (gene -> CoreGene).
    replication_datasets : dict
        Dataset_id -> gene-level expression DataFrame.
    replication_phenotypes : dict
        Dataset_id -> {sample_id: 'case' or 'control'}.
    dataset_tissues : dict
        Dataset_id -> tissue type string (e.g., 'brain', 'colon', 'islet').
    discovery_tissues : set of str, optional
        Tissue types used in discovery datasets. When provided, replication
        datasets with matching tissue are considered same-tissue. When None,
        all replication tissues are treated as same-tissue (conservative).
    p_threshold : float
        Significance threshold for elimination (default 0.05).

    Returns
    -------
    dict of gene_symbol -> GeneVerdict
    """
    verdicts = {}

    for gene, core_info in core_genes.items():
        discovery_dir = core_info.direction
        cluster_id = core_info.cluster_id
        rep_results = []

        # Test in each replication dataset
        for ds_id, expr_df in replication_datasets.items():
            tissue = dataset_tissues.get(ds_id, "other")
            pheno = replication_phenotypes.get(ds_id, {})

            rep = replicate_gene(
                gene, discovery_dir, expr_df, pheno, ds_id, tissue
            )
            if rep is not None:
                rep_results.append(rep)

        # Split results into same-tissue and cross-tissue
        if discovery_tissues is not None:
            same_tissue_results = [
                r for r in rep_results if r.tissue in discovery_tissues
            ]
            cross_tissue_results = [
                r for r in rep_results if r.tissue not in discovery_tissues
            ]
        else:
            # No discovery tissue info: treat all as same-tissue (conservative)
            same_tissue_results = rep_results
            cross_tissue_results = []

        n_same_concordant = sum(
            1 for r in same_tissue_results if r.is_concordant
        )
        n_same_discordant_sig = sum(
            1 for r in same_tissue_results
            if not r.is_concordant and r.is_significant
        )
        n_cross_concordant = sum(
            1 for r in cross_tissue_results if r.is_concordant
        )

        # Elimination decision
        if not rep_results:
            status = "insufficient_data"
            reason = "No replication data available for this gene."
        elif n_same_discordant_sig > 0:
            # ELIMINATE: significant opposite direction in same tissue
            discordant_ds = [
                r.dataset_id for r in same_tissue_results
                if not r.is_concordant and r.is_significant
            ]
            status = "eliminated"
            reason = (
                f"Significant opposite direction in same-tissue dataset(s): "
                f"{discordant_ds}. Discovery direction: {discovery_dir}."
            )
        else:
            status = "survived"
            if same_tissue_results:
                reason = (
                    f"Same-tissue replication: {n_same_concordant}/"
                    f"{len(same_tissue_results)} concordant. "
                    f"Cross-tissue: {n_cross_concordant}/"
                    f"{len(cross_tissue_results)} concordant."
                )
            else:
                reason = (
                    "No same-tissue replication datasets; retained by default. "
                    f"Cross-tissue: {n_cross_concordant}/"
                    f"{len(cross_tissue_results)} concordant."
                )

        verdicts[gene] = GeneVerdict(
            gene=gene,
            cluster_id=cluster_id,
            status=status,
            reason=reason,
            replication_results=rep_results,
            n_same_tissue_concordant=n_same_concordant,
            n_same_tissue_discordant=n_same_discordant_sig,
            n_cross_tissue_tested=len(cross_tissue_results),
            n_cross_tissue_concordant=n_cross_concordant,
            discovery_direction=discovery_dir,
        )

        logger.info(
            f"Gene {gene}: {status} — {reason}"
        )

    return verdicts


def assign_cluster_verdicts(
    gene_verdicts: dict[str, GeneVerdict],
    core_genes: dict,
) -> dict[int, ClusterVerdict]:
    """Assign replication verdicts per cluster.

    Parameters
    ----------
    gene_verdicts : dict
        From run_elimination_protocol().
    core_genes : dict
        From Phase4Result.core_genes.

    Returns
    -------
    dict of cluster_id -> ClusterVerdict
    """
    # Group genes by cluster
    cluster_genes: dict[int, list[str]] = {}
    for gene, core_info in core_genes.items():
        cid = core_info.cluster_id
        if cid not in cluster_genes:
            cluster_genes[cid] = []
        cluster_genes[cid].append(gene)

    verdicts = {}
    for cid, genes in cluster_genes.items():
        survived = []
        eliminated = []

        for gene in genes:
            if gene in gene_verdicts:
                gv = gene_verdicts[gene]
                if gv.status == "survived":
                    survived.append(gene)
                elif gv.status == "eliminated":
                    eliminated.append(gene)
                # insufficient_data genes are neither survived nor eliminated

        n_core = len(genes)
        n_surv = len(survived)
        n_elim = len(eliminated)

        # Determine cluster verdict
        if n_surv == 0:
            verdict = "failed"
        elif n_surv == n_core:
            # Check if blood failed (brain-specific pattern)
            has_blood_fail = False
            for gene in survived:
                gv = gene_verdicts[gene]
                if gv.n_cross_tissue_tested > 0 and gv.n_cross_tissue_concordant == 0:
                    has_blood_fail = True
                    break
            verdict = "brain_specific" if has_blood_fail else "replicated"
        elif n_surv >= n_core * 0.5:
            verdict = "partially_replicated"
        else:
            verdict = "failed"

        verdicts[cid] = ClusterVerdict(
            cluster_id=cid,
            verdict=verdict,
            n_core_genes=n_core,
            n_survived=n_surv,
            n_eliminated=n_elim,
            survived_genes=survived,
            eliminated_genes=eliminated,
        )

        logger.info(
            f"Cluster {cid}: {verdict} "
            f"({n_surv}/{n_core} survived, {n_elim} eliminated)"
        )

    return verdicts


def run_phase5(
    core_genes: dict,
    replication_datasets: dict[str, pd.DataFrame],
    replication_phenotypes: dict[str, dict[str, str]],
    dataset_tissues: dict[str, str],
    discovery_tissues: set[str] | None = None,
) -> Phase5Result:
    """Run Phase 5: Independent Replication.

    Parameters
    ----------
    core_genes : dict
        From Phase4Result.core_genes. This list is LOCKED —
        no modifications allowed based on replication results.
    replication_datasets : dict
        Dataset_id -> gene-level expression DataFrame.
    replication_phenotypes : dict
        Dataset_id -> {sample_id: 'case' or 'control'}.
    dataset_tissues : dict
        Dataset_id -> tissue type string (e.g., 'brain', 'colon', 'islet').
    discovery_tissues : set of str, optional
        Tissue types from discovery datasets. Same-tissue replication
        datasets can trigger elimination; cross-tissue cannot.

    Returns
    -------
    Phase5Result
    """
    # Record the locked gene list
    locked_list = sorted(core_genes.keys())

    logger.info(
        f"Phase 5: Testing {len(locked_list)} pre-specified core genes "
        f"in {len(replication_datasets)} replication datasets."
    )

    # Run elimination protocol
    gene_verdicts = run_elimination_protocol(
        core_genes, replication_datasets,
        replication_phenotypes, dataset_tissues,
        discovery_tissues=discovery_tissues,
    )

    # Assign cluster verdicts
    cluster_verdicts = assign_cluster_verdicts(gene_verdicts, core_genes)

    n_survived = sum(1 for v in gene_verdicts.values() if v.status == "survived")
    n_eliminated = sum(1 for v in gene_verdicts.values() if v.status == "eliminated")
    n_insufficient = sum(1 for v in gene_verdicts.values() if v.status == "insufficient_data")

    logger.info(
        f"Phase 5 complete: {n_survived} survived, {n_eliminated} eliminated, "
        f"{n_insufficient} insufficient data."
    )

    return Phase5Result(
        gene_verdicts=gene_verdicts,
        cluster_verdicts=cluster_verdicts,
        n_survived=n_survived,
        n_eliminated=n_eliminated,
        n_insufficient=n_insufficient,
        locked_core_genes=locked_list,
    )
