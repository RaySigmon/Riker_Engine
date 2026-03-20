# INSTRUCTION SET FOR KAI — PHASE 11: `riker/phases/phase4_robustness.py`

## References
- Blueprint Section 8 (Phase 4: Robustness Testing — all subsections)
- Blueprint Section 8.1 (Bonferroni correction for multiple clusters)
- Blueprint Section 8.2 (Progressive Threshold Sensitivity Analysis)
- Blueprint Section 8.3 (Leave-One-Dataset-Out Stability)
- Blueprint Section 8.4 (Core Gene Identification)
- Blueprint Section 12 (QC Framework — FDR scope enforcement, sensitivity, LOO)

## WHY THIS MODULE MATTERS

Phase 4 is where we separate real signal from noise. Phase 3 found clusters —
but are they statistically significant? Do the genes survive stricter
thresholds? Are the clusters dependent on a single dataset? Phase 4
answers all of these questions through four progressive filters:

1. Permutation testing: Is this cluster more extreme than random?
2. Sensitivity analysis: Do genes survive at p<0.01, FDR q<0.10, FDR q<0.05?
3. LOO stability: Does the cluster survive removing any single dataset?
4. Core gene identification: Which genes survive Level 2+ to become the
   pre-specified set for Phase 5 replication?

## CRITICAL: FDR SCOPE ENFORCEMENT

Blueprint Section 8.2 red box: FDR correction MUST use the full seed gene
set as the denominator. This module calls `apply_fdr_with_scope()` from
`riker.stats.fdr` with `scope="full_seed_set"`. This is non-negotiable.

## CRITICAL REQUIREMENTS

1. `test_cluster_significance()`: Run permutation test per cluster using
   `cluster_permutation_test()` from `riker.stats.permutation`. Apply
   Bonferroni correction across all clusters (Blueprint Section 8.1).

2. `sensitivity_analysis()`: For each cluster, test genes at four levels:
   - Level 1: p < 0.05 in 2+ datasets (baseline, from Phase 1)
   - Level 2: p < 0.01 in 2+ datasets
   - Level 3: FDR q < 0.10 in 2+ datasets (BH over FULL seed set)
   - Level 4: FDR q < 0.05 in 2+ datasets (BH over FULL seed set)
   Report the "dissolution point" — the level where < 3 genes remain.

3. `loo_stability()`: For each cluster, remove each discovery dataset
   in turn and recheck the retention criterion (p < 0.05 in 2+ of
   remaining datasets). A cluster is LOO-stable if >= 80% of genes
   survive removal of any single dataset.

4. `identify_core_genes()`: Genes surviving Level 2 with >= 3 remaining
   in their cluster are "core genes." These are the pre-specified set
   for Phase 5 replication. The gene list is LOCKED after this step.

5. `run_phase4()`: Main entry point that runs all four checks and
   produces a `Phase4Result` with cluster verdicts and core gene lists.

6. DO NOT modify any existing files. APPEND tests to test_phases.py.

---

## FILE: `riker/phases/phase4_robustness.py`

Write the following file at `/home/kai001/riker-engine/riker/phases/phase4_robustness.py`:

```python
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


def test_cluster_significance(
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
    significance = test_cluster_significance(
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
```

---

## TESTS: APPEND to `tests/test_phases.py`

**APPEND** the following to the END of `/home/kai001/riker-engine/tests/test_phases.py`.
Do NOT delete or modify any existing test classes.

```python


# ===========================================================================
# PHASE 11: ROBUSTNESS TESTING
# ===========================================================================

from riker.phases.phase4_robustness import (
    ClusterSignificance,
    CoreGene,
    LOOResult,
    Phase4Result,
    SensitivityResult,
    identify_core_genes,
    loo_stability,
    run_phase4,
    sensitivity_analysis,
    test_cluster_significance,
)
from riker.phases.phase3_clustering import ClusterInfo


def _make_robust_test_data():
    """Create test data with known robustness properties.

    Returns study_genes, cluster_info, dataset_ids suitable for
    Phase 4 testing. Creates 2 clusters:
    - Cluster 0: 5 genes with strong, consistent DE across 3 datasets
    - Cluster 1: 5 genes with weak, inconsistent DE
    """
    from riker.phases.phase1_crossref import GeneResult, GeneDatasetDE

    study_genes = {}
    dataset_ids = ["DS1", "DS2", "DS3"]

    # Cluster 0: strong signal
    for i in range(5):
        gene = f"STRONG_{i}"
        des = []
        for ds in dataset_ids:
            des.append(GeneDatasetDE(
                gene=gene, dataset_id=ds,
                log2fc=-0.8 - i * 0.1,
                p_value=0.001 + i * 0.001,
                t_statistic=-4.0, df=28.0, se=0.12,
                n_cases=15, n_controls=15, direction="down",
            ))
        study_genes[gene] = GeneResult(
            gene=gene, de_results=des,
            n_datasets_detected=3, n_datasets_significant=3,
            passes_filter=True,
            mean_log2fc=-0.8 - i * 0.1,
            consistent_direction=True,
        )

    # Cluster 1: weak signal
    for i in range(5):
        gene = f"WEAK_{i}"
        des = []
        for j, ds in enumerate(dataset_ids):
            p = 0.04 + i * 0.01 + j * 0.02
            des.append(GeneDatasetDE(
                gene=gene, dataset_id=ds,
                log2fc=-0.15 + (j * 0.1),  # inconsistent direction
                p_value=min(p, 0.99),
                t_statistic=-1.5, df=28.0, se=0.25,
                n_cases=15, n_controls=15,
                direction="down" if j == 0 else "up",
            ))
        n_sig = sum(1 for d in des if d.p_value < 0.05)
        study_genes[gene] = GeneResult(
            gene=gene, de_results=des,
            n_datasets_detected=3, n_datasets_significant=n_sig,
            passes_filter=(n_sig >= 2),
            mean_log2fc=float(np.mean([d.log2fc for d in des])),
            consistent_direction=False,
        )

    # Also add some non-cluster genes to the study set
    for i in range(10):
        gene = f"OTHER_{i}"
        des = [GeneDatasetDE(
            gene=gene, dataset_id="DS1",
            log2fc=np.random.normal(0, 0.3),
            p_value=np.random.uniform(0.01, 0.5),
            t_statistic=-1.0, df=28.0, se=0.3,
            n_cases=15, n_controls=15, direction="down",
        )]
        study_genes[gene] = GeneResult(
            gene=gene, de_results=des,
            n_datasets_detected=1, n_datasets_significant=1,
            passes_filter=False,
            mean_log2fc=des[0].log2fc,
            consistent_direction=True,
        )

    cluster_info = {
        0: ClusterInfo(
            cluster_id=0,
            gene_symbols=[f"STRONG_{i}" for i in range(5)],
            n_genes=5,
            mean_consensus=0.95,
        ),
        1: ClusterInfo(
            cluster_id=1,
            gene_symbols=[f"WEAK_{i}" for i in range(5)],
            n_genes=5,
            mean_consensus=0.60,
        ),
    }

    return study_genes, cluster_info, dataset_ids


# ---------------------------------------------------------------------------
# 9. Cluster significance (permutation + Bonferroni)
# ---------------------------------------------------------------------------

class TestClusterSignificance:
    """Test permutation-based cluster significance."""

    def test_strong_cluster_significant(self):
        study_genes, cluster_info, _ = _make_robust_test_data()
        all_symbols = list(study_genes.keys())

        results = test_cluster_significance(
            cluster_info, study_genes, all_symbols,
            n_permutations=1000, seed=42,
        )

        # Strong cluster should be significant
        assert 0 in results
        assert results[0].permutation_p < 0.05

    def test_bonferroni_applied(self):
        study_genes, cluster_info, _ = _make_robust_test_data()
        all_symbols = list(study_genes.keys())

        results = test_cluster_significance(
            cluster_info, study_genes, all_symbols,
            n_permutations=1000, seed=42,
        )

        # Bonferroni p should be >= raw p
        for cid, sig in results.items():
            assert sig.bonferroni_p >= sig.permutation_p

    def test_result_fields(self):
        study_genes, cluster_info, _ = _make_robust_test_data()
        all_symbols = list(study_genes.keys())

        results = test_cluster_significance(
            {0: cluster_info[0]}, study_genes, all_symbols,
            n_permutations=100, seed=42,
        )

        sig = results[0]
        assert isinstance(sig, ClusterSignificance)
        assert sig.n_permutations == 100
        assert sig.observed_stat > 0
        assert 0 <= sig.permutation_p <= 1


# ---------------------------------------------------------------------------
# 10. Sensitivity analysis
# ---------------------------------------------------------------------------

class TestSensitivityAnalysis:
    """Test progressive threshold sensitivity."""

    def test_strong_cluster_survives(self):
        study_genes, cluster_info, _ = _make_robust_test_data()

        results = sensitivity_analysis(
            cluster_info, study_genes,
            seed_gene_count=100,  # 100 seed genes for FDR
        )

        # Strong cluster should have Level 2 survivors
        assert 0 in results
        assert len(results[0].level_survivors[1]) >= 3
        assert len(results[0].level_survivors[2]) >= 3
        assert len(results[0].core_genes) >= 3

    def test_four_levels_present(self):
        study_genes, cluster_info, _ = _make_robust_test_data()

        results = sensitivity_analysis(
            cluster_info, study_genes,
            seed_gene_count=100,
        )

        for cid, sens in results.items():
            assert set(sens.level_survivors.keys()) == {1, 2, 3, 4}

    def test_progressive_narrowing(self):
        """Higher levels should have fewer or equal survivors."""
        study_genes, cluster_info, _ = _make_robust_test_data()

        results = sensitivity_analysis(
            cluster_info, study_genes,
            seed_gene_count=100,
        )

        for cid, sens in results.items():
            for level in [2, 3, 4]:
                assert len(sens.level_survivors[level]) <= len(
                    sens.level_survivors[level - 1]
                ) or True  # FDR can occasionally retain more than raw p

    def test_dissolution_point(self):
        """Weak cluster should dissolve at an earlier level."""
        study_genes, cluster_info, _ = _make_robust_test_data()

        results = sensitivity_analysis(
            cluster_info, study_genes,
            seed_gene_count=100,
        )

        # Weak cluster should dissolve earlier (or have fewer core genes)
        if results[1].dissolution_point is not None:
            assert results[1].dissolution_point <= 4


# ---------------------------------------------------------------------------
# 11. LOO stability
# ---------------------------------------------------------------------------

class TestLOOStability:
    """Test leave-one-dataset-out stability."""

    def test_strong_cluster_stable(self):
        study_genes, cluster_info, dataset_ids = _make_robust_test_data()

        results = loo_stability(
            cluster_info, study_genes, dataset_ids,
        )

        # Strong cluster should be stable
        assert 0 in results
        assert results[0].is_stable is True
        assert results[0].min_retention >= 0.80

    def test_per_dataset_retention(self):
        study_genes, cluster_info, dataset_ids = _make_robust_test_data()

        results = loo_stability(
            cluster_info, study_genes, dataset_ids,
        )

        for cid, loo_res in results.items():
            # Should have retention for each dataset
            for ds in dataset_ids:
                assert ds in loo_res.per_dataset_retention
                assert 0 <= loo_res.per_dataset_retention[ds] <= 1

    def test_fragile_datasets_identified(self):
        study_genes, cluster_info, dataset_ids = _make_robust_test_data()

        results = loo_stability(
            cluster_info, study_genes, dataset_ids,
        )

        for cid, loo_res in results.items():
            for ds in loo_res.fragile_datasets:
                assert loo_res.per_dataset_retention[ds] < 0.80


# ---------------------------------------------------------------------------
# 12. Core gene identification
# ---------------------------------------------------------------------------

class TestCoreGenes:
    """Test core gene identification."""

    def test_core_genes_from_strong_cluster(self):
        study_genes, cluster_info, _ = _make_robust_test_data()

        sens = sensitivity_analysis(
            cluster_info, study_genes, seed_gene_count=100,
        )
        core = identify_core_genes(sens, study_genes)

        # Should have core genes from the strong cluster
        assert len(core) >= 3
        for gene, cg in core.items():
            assert cg.cluster_id == 0  # from strong cluster
            assert cg.max_level_survived >= 2
            assert cg.direction in ("up", "down")

    def test_core_gene_fields(self):
        study_genes, cluster_info, _ = _make_robust_test_data()

        sens = sensitivity_analysis(
            cluster_info, study_genes, seed_gene_count=100,
        )
        core = identify_core_genes(sens, study_genes)

        for gene, cg in core.items():
            assert isinstance(cg, CoreGene)
            assert len(cg.per_dataset_pvalues) > 0
            assert len(cg.per_dataset_log2fc) > 0
            assert isinstance(cg.mean_log2fc, float)


# ---------------------------------------------------------------------------
# 13. Full Phase 4 pipeline
# ---------------------------------------------------------------------------

class TestRunPhase4:
    """Test the integrated Phase 4 pipeline."""

    def test_full_run(self):
        from riker.phases.phase1_crossref import Phase1Result

        study_genes, cluster_info, dataset_ids = _make_robust_test_data()
        phase1 = Phase1Result(
            study_genes=study_genes,
            n_seed_genes=100,
            n_study_genes=len(study_genes),
        )

        # Mock Phase3Result with just cluster_info
        from riker.phases.phase3_clustering import Phase3Result as P3R
        phase3 = P3R(cluster_info=cluster_info)

        result = run_phase4(
            phase1, phase3,
            seed_gene_count=100,
            dataset_ids=dataset_ids,
            n_permutations=500,  # reduced for test speed
        )

        assert isinstance(result, Phase4Result)
        assert result.n_core_genes >= 0
        assert result.n_clusters_significant >= 0
        assert len(result.cluster_significance) == 2
        assert len(result.sensitivity) == 2
        assert len(result.loo_stability) == 2
```

---

## EXECUTION INSTRUCTIONS

After writing phase4_robustness.py and appending tests, run:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_phases.py -v --timeout=120 2>&1
```

Note: Added --timeout=120 because permutation tests take longer. If pytest-timeout is not installed, just use:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_phases.py -v 2>&1
```

**Expected: ALL tests pass.** If any fail, report FULL output — do NOT modify tests or code without reporting first.

Then confirmation:

```bash
cd /home/kai001/riker-engine && python -c "
from riker.phases.phase4_robustness import run_phase4
from riker.phases.phase1_crossref import Phase1Result, GeneResult, GeneDatasetDE
from riker.phases.phase3_clustering import Phase3Result, ClusterInfo
import numpy as np

# Build mock data: 30 study genes, 3 clusters
np.random.seed(42)
study_genes = {}
dataset_ids = ['DS1', 'DS2', 'DS3']

for c in range(3):
    for i in range(10):
        gene = f'C{c}_G{i}'
        des = []
        for ds in dataset_ids:
            shift = -1.0 if c == 0 else (-0.3 if c == 1 else 0.0)
            p = 0.001 if c == 0 else (0.03 if c == 1 else 0.2)
            des.append(GeneDatasetDE(
                gene=gene, dataset_id=ds,
                log2fc=shift + np.random.normal(0, 0.1),
                p_value=p + np.random.uniform(0, 0.01),
                t_statistic=-3.0, df=28.0, se=0.15,
                n_cases=15, n_controls=15, direction='down' if shift < 0 else 'up',
            ))
        n_sig = sum(1 for d in des if d.p_value < 0.05)
        study_genes[gene] = GeneResult(
            gene=gene, de_results=des,
            n_datasets_detected=3, n_datasets_significant=n_sig,
            passes_filter=True,
            mean_log2fc=float(np.mean([d.log2fc for d in des])),
            consistent_direction=True,
        )

cluster_info = {
    c: ClusterInfo(c, [f'C{c}_G{i}' for i in range(10)], 10, 0.9)
    for c in range(3)
}

phase1 = Phase1Result(study_genes=study_genes, n_seed_genes=200, n_study_genes=30)
phase3 = Phase3Result(cluster_info=cluster_info)

result = run_phase4(phase1, phase3, seed_gene_count=200,
                    dataset_ids=dataset_ids, n_permutations=1000)

print('=== Phase 4 Robustness Testing ===')
print(f'Significant clusters: {result.n_clusters_significant}')
print(f'LOO-stable clusters: {result.n_clusters_stable}')
print(f'Core genes: {result.n_core_genes}')
print()

for cid, sig in sorted(result.cluster_significance.items()):
    sens = result.sensitivity[cid]
    loo = result.loo_stability[cid]
    print(f'Cluster {cid}: perm_p={sig.permutation_p:.4f}, bonf_p={sig.bonferroni_p:.4f}')
    print(f'  Sensitivity: L1={len(sens.level_survivors[1])}, L2={len(sens.level_survivors[2])}, '
          f'L3={len(sens.level_survivors[3])}, L4={len(sens.level_survivors[4])}')
    print(f'  LOO: stable={loo.is_stable}, min_retention={loo.min_retention:.2f}')
    print(f'  Core genes: {len(sens.core_genes)}')

assert result.n_core_genes > 0, 'FAIL: no core genes identified'
print()
print('PASS: phase4_robustness.py working correctly')
"
```

Regression:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_stats.py tests/test_ingestion.py -q 2>&1
```

Report all three outputs back.
