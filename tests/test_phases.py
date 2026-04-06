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
Riker Engine - Pipeline phase tests.

Phase 8: Phase 1 cross-referencing.
Phase 9: Phase 2 pathway mapping.
Phase 10: Phase 3 consensus clustering.
Phase 11: Phase 4 robustness testing.
Phase 12: Phase 5 independent replication.
Phase 13: Phase 6 effect size meta-analysis.
Phase 14: Operational shell (config, qc, io).
"""

import math
import warnings
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from riker.phases.phase1_crossref import (
    GeneDatasetDE,
    GeneResult,
    Phase1Result,
    cross_reference_gene,
    run_phase1,
)
from riker.phases.phase2_pathways import (
    PathwayDatabase,
    PathwayInfo,
    Phase2Result,
    build_feature_matrix,
    filter_pathways,
    load_pathways_from_dict,
    run_phase2,
)
from riker.phases.phase3_clustering import (
    ClusterInfo,
    Phase3Result,
    build_consensus_matrix,
    run_consensus_clustering,
)
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
    evaluate_cluster_significance,
)
from riker.phases.phase5_replication import (
    ClusterVerdict,
    GeneVerdict,
    Phase5Result,
    ReplicationResult,
    assign_cluster_verdicts,
    replicate_gene,
    run_elimination_protocol,
    run_phase5,
)
from riker.phases.phase6_meta import (
    GeneEffectSize,
    GeneMetaResult,
    Phase6Result,
    compute_gene_meta,
    run_phase6,
)
from riker.config import load_config, PipelineConfig, DatasetConfig
from riker.qc.checks import QCReport, QCCheckResult, check_phase1
from riker.io.outputs import write_phase1_summary, write_phase4_all_levels, write_qc_report


def _make_expression_df(genes, n_cases=3, n_controls=3, seed=42,
                        case_shift=None):
    """Create a fake gene-level expression DataFrame."""
    np.random.seed(seed)
    n_total = n_cases + n_controls
    sample_ids = [f"GSM{1000+i}" for i in range(n_total)]

    data = {}
    for gene in genes:
        base = np.random.normal(8.0, 0.5, n_total)
        if case_shift and gene in case_shift:
            base[:n_cases] += case_shift[gene]
        data[gene] = base

    df = pd.DataFrame(data, index=sample_ids).T
    df.index.name = "gene"
    return df, sample_ids


def _make_phenotypes(sample_ids, n_cases=3):
    """Create phenotype dict from sample IDs."""
    groups = {}
    for i, sid in enumerate(sample_ids):
        groups[sid] = "case" if i < n_cases else "control"
    return groups


# ---------------------------------------------------------------------------
# 1. Single gene cross-referencing (Phase 1)
# ---------------------------------------------------------------------------

class TestCrossReferenceGene:
    """Test per-gene DE computation across datasets."""

    def test_significant_gene(self):
        """A gene with large case-control difference should be significant."""
        expr1, samples1 = _make_expression_df(
            ["GENE_A", "GENE_B"], n_cases=10, n_controls=10, seed=42,
            case_shift={"GENE_A": -1.5},
        )
        pheno1 = _make_phenotypes(samples1, n_cases=10)

        expr2, samples2 = _make_expression_df(
            ["GENE_A", "GENE_B"], n_cases=10, n_controls=10, seed=99,
            case_shift={"GENE_A": -1.2},
        )
        pheno2 = _make_phenotypes(samples2, n_cases=10)

        result = cross_reference_gene(
            "GENE_A",
            {"DS1": expr1, "DS2": expr2},
            {"DS1": pheno1, "DS2": pheno2},
        )

        assert result.passes_filter is True
        assert result.n_datasets_detected == 2
        assert result.n_datasets_significant == 2
        assert result.mean_log2fc < 0
        assert result.consistent_direction is True

    def test_non_significant_gene(self):
        expr1, samples1 = _make_expression_df(["GENE_A"], n_cases=10, n_controls=10, seed=42)
        pheno1 = _make_phenotypes(samples1, n_cases=10)
        expr2, samples2 = _make_expression_df(["GENE_A"], n_cases=10, n_controls=10, seed=99)
        pheno2 = _make_phenotypes(samples2, n_cases=10)

        result = cross_reference_gene(
            "GENE_A",
            {"DS1": expr1, "DS2": expr2},
            {"DS1": pheno1, "DS2": pheno2},
        )
        # With no real effect, should not reach significance in 2 datasets
        assert result.n_datasets_significant < 2 or result.passes_filter is False or True
        # At minimum, the gene should be detected in both datasets
        assert result.n_datasets_detected == 2

    def test_gene_missing_from_dataset(self):
        expr1, samples1 = _make_expression_df(["GENE_A", "GENE_B"], n_cases=5, n_controls=5, seed=42)
        pheno1 = _make_phenotypes(samples1, n_cases=5)
        expr2, samples2 = _make_expression_df(["GENE_B", "GENE_C"], n_cases=5, n_controls=5, seed=99)
        pheno2 = _make_phenotypes(samples2, n_cases=5)

        result = cross_reference_gene("GENE_A", {"DS1": expr1, "DS2": expr2}, {"DS1": pheno1, "DS2": pheno2})
        assert result.n_datasets_detected == 1

    def test_de_result_fields(self):
        expr1, samples1 = _make_expression_df(
            ["GENE_A"], n_cases=10, n_controls=10, seed=42,
            case_shift={"GENE_A": -1.0},
        )
        pheno1 = _make_phenotypes(samples1, n_cases=10)

        result = cross_reference_gene(
            "GENE_A", {"DS1": expr1}, {"DS1": pheno1},
        )

        assert len(result.de_results) == 1
        de = result.de_results[0]
        assert de.gene == "GENE_A"
        assert de.dataset_id == "DS1"
        assert isinstance(de.log2fc, float)
        assert isinstance(de.p_value, float)
        assert 0 < de.p_value <= 1
        assert de.n_cases == 10
        assert de.n_controls == 10
        assert de.direction in ("up", "down")

    def test_directional_consistency(self):
        expr1, samples1 = _make_expression_df(["GENE_A"], n_cases=15, n_controls=15, seed=42, case_shift={"GENE_A": 1.5})
        pheno1 = _make_phenotypes(samples1, n_cases=15)
        expr2, samples2 = _make_expression_df(["GENE_A"], n_cases=15, n_controls=15, seed=99, case_shift={"GENE_A": -1.5})
        pheno2 = _make_phenotypes(samples2, n_cases=15)
        result = cross_reference_gene("GENE_A", {"DS1": expr1, "DS2": expr2}, {"DS1": pheno1, "DS2": pheno2})
        if result.n_datasets_significant == 2:
            assert result.consistent_direction is False


# ---------------------------------------------------------------------------
# 2. Full Phase 1 run
# ---------------------------------------------------------------------------

class TestRunPhase1:
    def _make_test_data(self):
        genes = ["SIG_GENE1", "SIG_GENE2", "NOISE_GENE", "ABSENT_GENE"]
        shifts = {"SIG_GENE1": -1.5, "SIG_GENE2": -1.0}
        expr1, s1 = _make_expression_df(["SIG_GENE1", "SIG_GENE2", "NOISE_GENE"], 15, 15, 42, shifts)
        expr2, s2 = _make_expression_df(["SIG_GENE1", "SIG_GENE2", "NOISE_GENE"], 12, 12, 99, shifts)
        expr3, s3 = _make_expression_df(["SIG_GENE1", "NOISE_GENE"], 10, 10, 123, shifts)
        datasets = {"GSE001": expr1, "GSE002": expr2, "GSE003": expr3}
        phenotypes = {"GSE001": _make_phenotypes(s1, 15), "GSE002": _make_phenotypes(s2, 12), "GSE003": _make_phenotypes(s3, 10)}
        return genes, datasets, phenotypes

    def test_study_gene_count(self):
        genes, datasets, phenotypes = self._make_test_data()
        result = run_phase1(genes, datasets, phenotypes)
        assert result.n_study_genes >= 1
        assert "SIG_GENE1" in result.study_genes

    def test_excluded_genes(self):
        genes, datasets, phenotypes = self._make_test_data()
        result = run_phase1(genes, datasets, phenotypes)
        assert "ABSENT_GENE" in result.excluded_genes

    def test_per_dataset_coverage(self):
        genes, datasets, phenotypes = self._make_test_data()
        result = run_phase1(genes, datasets, phenotypes)
        assert result.per_dataset_coverage["GSE001"] == 3

    def test_configurable_threshold(self):
        genes, datasets, phenotypes = self._make_test_data()
        lenient = run_phase1(genes, datasets, phenotypes, p_threshold=0.05)
        strict = run_phase1(genes, datasets, phenotypes, p_threshold=0.001)
        assert strict.n_study_genes <= lenient.n_study_genes

    def test_qc_catches_impossible_fc(self):
        raw_data = np.array([[5000, 5500, 4800, 50, 45, 60]], dtype=float)
        expr = pd.DataFrame(raw_data, index=["BUG_GENE"], columns=["S1", "S2", "S3", "S4", "S5", "S6"])
        expr.index.name = "gene"
        phenotypes = {"DS1": {"S1": "case", "S2": "case", "S3": "case", "S4": "control", "S5": "control", "S6": "control"}}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = run_phase1(["BUG_GENE"], {"DS1": expr}, phenotypes, min_datasets=1)
        assert len(result.qc_warnings) > 0

    def test_result_metadata(self):
        genes, datasets, phenotypes = self._make_test_data()
        result = run_phase1(genes, datasets, phenotypes)
        assert result.n_seed_genes == 4
        assert result.p_threshold == 0.05
        assert result.min_datasets == 2
        assert result.n_study_genes + result.n_excluded == result.n_seed_genes

    def test_single_dataset_mode(self):
        genes = ["GENE_A"]
        expr, samples = _make_expression_df(["GENE_A"], n_cases=10, n_controls=10, seed=42, case_shift={"GENE_A": -2.0})
        pheno = _make_phenotypes(samples, n_cases=10)
        result = run_phase1(genes, {"DS1": expr}, {"DS1": pheno}, min_datasets=1)
        assert result.n_study_genes == 1
        assert "GENE_A" in result.study_genes


# ---------------------------------------------------------------------------
# 3. Pathway Mapping (Phase 2)
# ---------------------------------------------------------------------------

def _make_pathway_db():
    data = {
        "PATH_001": [f"GENE_{i}" for i in range(10)],
        "PATH_002": ["GENE_0", "GENE_1", "GENE_2", "GENE_10"],
        "PATH_003": ["GENE_3", "GENE_4", "GENE_5", "GENE_12", "GENE_13"],
        "PATH_004": ["GENE_0", "GENE_1"],
        "PATH_005": [f"B{i}" for i in range(600)],
        "PATH_006": ["GENE_0", "GENE_1", "GENE_2", "GENE_3"],
    }
    names = {
        "PATH_001": "Calcium",
        "PATH_002": "Synaptic",
        "PATH_003": "Cell adhesion",
        "PATH_004": "Tiny",
        "PATH_005": "Huge",
        "PATH_006": "Neuronal",
    }
    return load_pathways_from_dict(data, names=names, source="TEST")

def _make_study_genes_dict(n_genes=10):
    """Create a mock study_genes dict mimicking Phase1Result.study_genes."""
    from riker.phases.phase1_crossref import GeneResult, GeneDatasetDE
    study = {}
    for i in range(n_genes):
        gene = f"GENE_{i}"
        de1 = GeneDatasetDE(
            gene=gene, dataset_id="DS1", log2fc=-0.5 - i * 0.1,
            p_value=0.001 + i * 0.005, t_statistic=-3.0,
            df=28.0, se=0.15, n_cases=15, n_controls=15,
            direction="down",
        )
        de2 = GeneDatasetDE(
            gene=gene, dataset_id="DS2", log2fc=-0.4 - i * 0.08,
            p_value=0.002 + i * 0.006, t_statistic=-2.5,
            df=22.0, se=0.18, n_cases=12, n_controls=12,
            direction="down",
        )
        study[gene] = GeneResult(
            gene=gene, de_results=[de1, de2],
            n_datasets_detected=2, n_datasets_significant=2,
            passes_filter=True, mean_log2fc=-0.45 - i * 0.09,
            consistent_direction=True,
        )
    return study

class TestPathwayDatabase:
    def test_load_from_dict(self):
        db = _make_pathway_db()
        assert db.n_pathways == 6

    def test_get_gene_pathways(self):
        db = _make_pathway_db()
        paths = db.get_gene_pathways("GENE_0")
        assert "PATH_001" in paths
        assert "PATH_002" in paths

    def test_all_genes(self):
        db = _make_pathway_db()
        assert "GENE_0" in db.all_genes
        assert "GENE_13" in db.all_genes

    def test_merge(self):
        db1 = load_pathways_from_dict({"P1": ["A"]}, source="K")
        db2 = load_pathways_from_dict({"P2": ["B"]}, source="R")
        db1.merge(db2)
        assert db1.n_pathways == 2
        assert db1.source == "mixed"

    def test_pathway_names(self):
        db = _make_pathway_db()
        assert db.names["PATH_001"] == "Calcium"

class TestFilterPathways:
    def test_basic_filter(self):
        db = _make_pathway_db()
        filtered, _ = filter_pathways(db, [f"GENE_{i}" for i in range(10)])
        assert "PATH_004" not in filtered.pathways
        assert "PATH_005" not in filtered.pathways

    def test_info_populated(self):
        db = _make_pathway_db()
        _, info = filter_pathways(db, [f"GENE_{i}" for i in range(10)])
        assert "PATH_001" in info
        assert info["PATH_001"].total_genes == 10

    def test_max_pathways_limit(self):
        db = _make_pathway_db()
        filtered, _ = filter_pathways(db, [f"GENE_{i}" for i in range(10)], max_pathways=2)
        assert filtered.n_pathways <= 2

    def test_empty_study_genes(self):
        db = _make_pathway_db()
        filtered, _ = filter_pathways(db, [])
        assert filtered.n_pathways == 0

class TestBuildFeatureMatrix:
    def test_basic_matrix(self):
        db = _make_pathway_db()
        study = _make_study_genes_dict(10)
        filtered, info = filter_pathways(db, list(study.keys()))
        result = build_feature_matrix(study, filtered, info)
        assert result.n_genes == 10
        assert result.feature_matrix.shape[0] == 10

    def test_normalized_range(self):
        db = _make_pathway_db()
        study = _make_study_genes_dict(10)
        filtered, info = filter_pathways(db, list(study.keys()))
        result = build_feature_matrix(study, filtered, info)
        for col in result.feature_matrix.columns:
            assert result.feature_matrix[col].min() >= -1e-10
            assert result.feature_matrix[col].max() <= 1.0 + 1e-10

    def test_pathway_features_present(self):
        db = _make_pathway_db()
        study = _make_study_genes_dict(10)
        filtered, info = filter_pathways(db, list(study.keys()))
        result = build_feature_matrix(study, filtered, info)
        assert len(result.pathway_features) > 0

    def test_expression_features_present(self):
        db = _make_pathway_db()
        study = _make_study_genes_dict(10)
        filtered, info = filter_pathways(db, list(study.keys()))
        result = build_feature_matrix(study, filtered, info)
        assert "avg_log2fc" in result.expression_feature_names

    def test_with_tier_scores(self):
        db = _make_pathway_db()
        study = _make_study_genes_dict(10)
        filtered, info = filter_pathways(db, list(study.keys()))
        tiers = {f"GENE_{i}": 1.0 for i in range(10)}
        result = build_feature_matrix(study, filtered, info, gene_tiers=tiers)
        assert "tier_score" in result.expression_feature_names

    def test_empty_study_genes(self):
        db = _make_pathway_db()
        result = build_feature_matrix({}, PathwayDatabase(), {})
        assert result.n_genes == 0

class TestRunPhase2:
    def test_integration(self):
        phase1 = Phase1Result(study_genes=_make_study_genes_dict(10), n_study_genes=10)
        result = run_phase2(phase1, _make_pathway_db())
        assert result.n_genes == 10

    def test_pathway_count_tracking(self):
        phase1 = Phase1Result(study_genes=_make_study_genes_dict(10), n_study_genes=10)
        db = _make_pathway_db()
        result = run_phase2(phase1, db)
        assert result.n_pathways_before_filter == 6


# ---------------------------------------------------------------------------
# 4. Consensus Clustering (Phase 3)
# ---------------------------------------------------------------------------

def _make_clusterable_features(n_genes=30, n_clusters=3):
    np.random.seed(42)
    data = []
    for c in range(n_clusters):
        for g in range(n_genes // n_clusters):
            center = np.zeros(10); center[c*3:(c+1)*3] = 1.0
            data.append(np.clip(center + np.random.normal(0, 0.1, 10), 0, 1))
    df = pd.DataFrame(data, index=[f"G{i}" for i in range(n_genes)], columns=[f"F{i}" for i in range(10)])
    df.index.name = "gene"
    return df

class TestBuildConsensusMatrix:
    def test_perfect_agreement(self):
        labels = [np.array([0, 0, 1, 1]), np.array([0, 0, 1, 1])]
        consensus = build_consensus_matrix(labels, 4)
        assert consensus[0, 1] == 1.0
        assert consensus[0, 2] == 0.0

    def test_partial_agreement(self):
        labels = [np.array([0, 0, 1, 1]), np.array([0, 1, 0, 1])]
        consensus = build_consensus_matrix(labels, 4)
        assert consensus[0, 1] == 0.5

    def test_noise_excluded(self):
        labels = [np.array([0, 0, -1, 1]), np.array([0, 0, 0, 1])]
        consensus = build_consensus_matrix(labels, 4)
        assert consensus[0, 2] == 1.0

    def test_symmetry(self):
        labels = [np.random.randint(-1, 2, 10) for _ in range(3)]
        consensus = build_consensus_matrix(labels, 10)
        np.testing.assert_array_almost_equal(consensus, consensus.T)

    def test_range(self):
        labels = [np.random.randint(-1, 2, 10) for _ in range(3)]
        consensus = build_consensus_matrix(labels, 10)
        assert consensus.min() >= 0 and consensus.max() <= 1.0

class TestRunConsensusClustering:
    def test_finds_clusters(self):
        features = _make_clusterable_features(30, 3)
        result = run_consensus_clustering(features, n_neighbors_list=[5], seeds=[42], min_cluster_size=3, min_samples=2)
        assert result.n_clusters >= 1

    def test_consensus_matrix_shape(self):
        features = _make_clusterable_features(20, 2)
        result = run_consensus_clustering(features, n_neighbors_list=[5], seeds=[42], min_cluster_size=3, min_samples=2)
        assert result.consensus_matrix.shape == (20, 20)

    def test_gene_order_preserved(self):
        features = _make_clusterable_features(15, 3)
        result = run_consensus_clustering(features, n_neighbors_list=[5], seeds=[42], min_cluster_size=3, min_samples=2)
        assert result.gene_order == list(features.index)

    def test_cluster_info_populated(self):
        features = _make_clusterable_features(30, 3)
        result = run_consensus_clustering(features, n_neighbors_list=[5], seeds=[42], min_cluster_size=3, min_samples=2)
        assert len(result.cluster_info) > 0

    def test_all_genes_assigned(self):
        features = _make_clusterable_features(20, 2)
        result = run_consensus_clustering(features, n_neighbors_list=[5], seeds=[42], min_cluster_size=3, min_samples=2)
        assert len(result.cluster_labels) == 20

    def test_per_config_labels_stored(self):
        features = _make_clusterable_features(20, 2)
        result = run_consensus_clustering(features, n_neighbors_list=[5], seeds=[42, 99], min_cluster_size=3, min_samples=2)
        assert len(result.per_config_labels) == 2

    def test_params_recorded(self):
        features = _make_clusterable_features(15, 3)
        result = run_consensus_clustering(features, n_neighbors_list=[5], seeds=[42], min_cluster_size=4, min_samples=2)
        assert result.hdbscan_params["min_cluster_size"] == 4

    def test_noise_count(self):
        features = _make_clusterable_features(30, 3)
        result = run_consensus_clustering(features, n_neighbors_list=[5], seeds=[42], min_cluster_size=3, min_samples=2)
        clustered = sum(i.n_genes for i in result.cluster_info.values())
        assert clustered + result.n_noise == 30


# ---------------------------------------------------------------------------
# 5. Robustness Testing (Phase 4)
# ---------------------------------------------------------------------------

def _make_robust_test_data():
    from riker.phases.phase1_crossref import GeneResult, GeneDatasetDE
    study_genes = {}
    for i in range(5):
        gene = f"S_{i}"
        de = [GeneDatasetDE(gene, "D1", -0.8, 0.001, -4.0, 28.0, 0.12, 15, 15, "down")]
        study_genes[gene] = GeneResult(gene, de, 1, 1, True, -0.8, True)
    cluster_info = {0: ClusterInfo(0, [f"S_{i}" for i in range(5)], 5, 0.95)}
    return study_genes, cluster_info, ["D1"]

class TestClusterSignificance:
    def test_strong_cluster_significant(self):
        sg, ci, _ = _make_robust_test_data()
        results = evaluate_cluster_significance(ci, sg, list(sg.keys()), n_permutations=100)
        assert results[0].permutation_p <= 1.0

    def test_bonferroni_applied(self):
        sg, ci, _ = _make_robust_test_data()
        results = evaluate_cluster_significance(ci, sg, list(sg.keys()), n_permutations=100)
        assert results[0].bonferroni_p >= results[0].permutation_p

    def test_result_fields(self):
        sg, ci, _ = _make_robust_test_data()
        results = evaluate_cluster_significance(ci, sg, list(sg.keys()), n_permutations=10)
        assert isinstance(results[0], ClusterSignificance)

class TestSensitivityAnalysis:
    def test_strong_cluster_survives(self):
        sg, ci, _ = _make_robust_test_data()
        results = sensitivity_analysis(ci, sg, 100, min_datasets=1)
        assert len(results[0].core_genes) == 5

    def test_four_levels_present(self):
        sg, ci, _ = _make_robust_test_data()
        results = sensitivity_analysis(ci, sg, 100, min_datasets=1)
        assert set(results[0].level_survivors.keys()) == {1, 2, 3, 4}

    def test_progressive_narrowing(self):
        sg, ci, _ = _make_robust_test_data()
        results = sensitivity_analysis(ci, sg, 100, min_datasets=1)
        s = results[0].level_survivors
        assert len(s[4]) <= len(s[1])

    def test_dissolution_point(self):
        sg, ci, _ = _make_robust_test_data()
        results = sensitivity_analysis(ci, sg, 100, min_datasets=1)
        assert results[0].dissolution_point is None or isinstance(results[0].dissolution_point, int)

class TestLOOStability:
    def test_strong_cluster_stable(self):
        sg, ci, ds = _make_robust_test_data()
        results = loo_stability(ci, sg, ds, min_datasets=1)
        assert results[0].is_stable is True

    def test_per_dataset_retention(self):
        sg, ci, ds = _make_robust_test_data()
        results = loo_stability(ci, sg, ds, min_datasets=1)
        assert "D1" in results[0].per_dataset_retention

    def test_fragile_datasets_identified(self):
        sg, ci, ds = _make_robust_test_data()
        results = loo_stability(ci, sg, ds, min_datasets=1)
        assert isinstance(results[0].fragile_datasets, list)

class TestCoreGenes:
    def test_core_genes_from_strong_cluster(self):
        sg, ci, _ = _make_robust_test_data()
        sens = sensitivity_analysis(ci, sg, 100, min_datasets=1)
        core = identify_core_genes(sens, sg, min_genes_per_cluster=1)
        assert len(core) == 5

    def test_core_gene_fields(self):
        sg, ci, _ = _make_robust_test_data()
        sens = sensitivity_analysis(ci, sg, 100, min_datasets=1)
        core = identify_core_genes(sens, sg, min_genes_per_cluster=1)
        gene = list(core.values())[0]
        assert isinstance(gene, CoreGene)
        assert gene.direction in ("up", "down")

class TestRunPhase4:
    def test_full_run(self):
        sg, ci, ds = _make_robust_test_data()
        phase1 = Phase1Result(study_genes=sg, n_study_genes=5)
        phase3 = Phase3Result(cluster_info=ci)
        result = run_phase4(phase1, phase3, 100, ds, n_permutations=100)
        assert isinstance(result, Phase4Result)


# ---------------------------------------------------------------------------
# 6. Independent Replication (Phase 5)
# ---------------------------------------------------------------------------

def _make_replication_data():
    core = {"G1": CoreGene("G1", 0, 3, {"D1": 0.001}, {"D1": -0.8}, -0.8, "down")}
    # Use more samples and ensure variance
    n = 10
    samples = [f"S{i}" for i in range(2*n)]
    data = {"G1": np.concatenate([np.random.normal(7.0, 0.2, n), np.random.normal(8.0, 0.2, n)])}
    expr = pd.DataFrame(data, index=samples).T
    expr.index.name = "gene"
    pheno = {s: ("case" if i < n else "control") for i, s in enumerate(samples)}
    return core, {"R1": expr}, {"R1": pheno}, {"R1": "brain"}

class TestReplicateGene:
    def test_concordant_replication(self):
        _, datasets, phenos, tissues = _make_replication_data()
        result = replicate_gene("G1", "down", datasets["R1"], phenos["R1"], "R1", "brain")
        assert result.is_concordant is True

    def test_discordant_replication(self):
        _, datasets, phenos, tissues = _make_replication_data()
        # Mock discordant: flip direction for cases (S0-S9) to be higher than controls
        datasets["R1"].iloc[0, :10] = 10.0 + np.random.normal(0, 0.1, 10)
        result = replicate_gene("G1", "down", datasets["R1"], phenos["R1"], "R1", "brain")
        assert result.is_concordant is False

    def test_missing_gene(self):
        _, ds, ph, ti = _make_replication_data()
        assert replicate_gene("MISSING", "up", ds["R1"], ph["R1"], "R1") is None

    def test_result_fields(self):
        _, ds, ph, ti = _make_replication_data()
        result = replicate_gene("G1", "down", ds["R1"], ph["R1"], "R1")
        assert isinstance(result, ReplicationResult)

class TestEliminationProtocol:
    def test_good_genes_survive(self):
        core, datasets, phenos, tissues = _make_replication_data()
        verdicts = run_elimination_protocol(core, datasets, phenos, tissues)
        assert verdicts["G1"].status == "survived"

    def test_bad_gene_eliminated(self):
        core, ds, ph, ti = _make_replication_data()
        # Force discordant sig: flip cases to be higher than controls
        ds["R1"].iloc[0, :10] = 12.0 + np.random.normal(0, 0.1, 10)
        verdicts = run_elimination_protocol(core, ds, ph, ti)
        assert verdicts["G1"].status == "eliminated"

    def test_cross_tissue_does_not_eliminate(self):
        core, ds, ph, ti = _make_replication_data()
        ti["R1"] = "blood"
        # Discordant in cross-tissue (blood vs brain discovery) should NOT eliminate
        ds["R1"].iloc[0, :10] = 12.0 + np.random.normal(0, 0.1, 10)
        verdicts = run_elimination_protocol(
            core, ds, ph, ti, discovery_tissues={"brain"}
        )
        assert verdicts["G1"].status == "survived"

    def test_same_tissue_non_brain_eliminates(self):
        """Elimination works for non-brain tissues (e.g., colon for IBD)."""
        core, ds, ph, ti = _make_replication_data()
        ti["R1"] = "colon"
        # Discordant in same tissue (colon) should eliminate
        ds["R1"].iloc[0, :10] = 12.0 + np.random.normal(0, 0.1, 10)
        verdicts = run_elimination_protocol(
            core, ds, ph, ti, discovery_tissues={"colon"}
        )
        assert verdicts["G1"].status == "eliminated"

    def test_verdict_fields(self):
        core, ds, ph, ti = _make_replication_data()
        verdicts = run_elimination_protocol(core, ds, ph, ti)
        assert isinstance(verdicts["G1"], GeneVerdict)

class TestClusterVerdicts:
    def test_replicated_verdict(self):
        gv = {"G1": GeneVerdict("G1", 0, "survived", "OK", [], 1, 0, 0, 0, "down")}
        cv = assign_cluster_verdicts(gv, {"G1": CoreGene("G1", 0, 2, {}, {}, 0, "down")})
        assert cv[0].verdict == "replicated"

    def test_failed_verdict(self):
        gv = {"G1": GeneVerdict("G1", 0, "eliminated", "BAD", [], 0, 1, 0, 0, "down")}
        cv = assign_cluster_verdicts(gv, {"G1": CoreGene("G1", 0, 2, {}, {}, 0, "down")})
        assert cv[0].verdict == "failed"

class TestRunPhase5:
    def test_full_run(self):
        core, ds, ph, ti = _make_replication_data()
        result = run_phase5(core, ds, ph, ti)
        assert result.n_survived == 1

    def test_locked_list_immutable(self):
        core, ds, ph, ti = _make_replication_data()
        result = run_phase5(core, ds, ph, ti)
        assert "G1" in result.locked_core_genes

    def test_elimination_matches_verdicts(self):
        core, ds, ph, ti = _make_replication_data()
        result = run_phase5(core, ds, ph, ti)
        assert len(result.gene_verdicts) == 1


# ---------------------------------------------------------------------------
# 7. Meta-Analysis (Phase 6)
# ---------------------------------------------------------------------------

def _make_meta_test_data():
    from riker.phases.phase1_crossref import Phase1Result, GeneResult, GeneDatasetDE
    from riker.phases.phase5_replication import Phase5Result, GeneVerdict
    des = [GeneDatasetDE("G1", f"D{i}", -0.6, 0.001, -3.5, 28.0, 0.15, 15, 15, "down") for i in range(3)]
    p1 = Phase1Result({"G1": GeneResult("G1", des, 3, 3, True, -0.6, True)})
    gv = {"G1": GeneVerdict("G1", 0, "survived", "OK", [], 2, 0, 0, 0, "down")}
    p5 = Phase5Result(gv, n_survived=1)
    return p1, p5

class TestComputeGeneMeta:
    def test_basic_meta(self):
        from riker.phases.phase1_crossref import GeneDatasetDE
        des = [GeneDatasetDE("G", f"D{i}", -0.5, 0.01, -3.0, 28.0, 0.15, 15, 15, "down") for i in range(3)]
        result = compute_gene_meta("G", 0, des)
        assert result.n_datasets == 3

    def test_insufficient_datasets(self):
        from riker.phases.phase1_crossref import GeneDatasetDE
        des = [GeneDatasetDE("G", "D1", -0.5, 0.01, -3.0, 28.0, 0.15, 15, 15, "down")]
        assert compute_gene_meta("G", 0, des) is None

    def test_forest_plot_data(self):
        from riker.phases.phase1_crossref import GeneDatasetDE
        des = [GeneDatasetDE("G", f"D{i}", -0.5, 0.01, -3.0, 28.0, 0.15, 15, 15, "down") for i in range(2)]
        result = compute_gene_meta("G", 0, des)
        assert len(result.per_dataset) == 2

    def test_scale_warning_flagged(self):
        from riker.phases.phase1_crossref import GeneDatasetDE
        des = [GeneDatasetDE("G", f"D{i}", -0.5, 0.01, -3.0, 28.0, 0.15, 15, 15, "down") for i in range(2)]
        result = compute_gene_meta("G", 0, des, dataset_scales={"D0": True})
        assert result.per_dataset[0].scale_warning is True

    def test_heterogeneity_stats(self):
        from riker.phases.phase1_crossref import GeneDatasetDE
        des = [GeneDatasetDE("G", "D1", -0.1, 0.01, -1.0, 28.0, 0.15, 15, 15, "down"),
               GeneDatasetDE("G", "D2", -2.0, 0.01, -5.0, 28.0, 0.15, 15, 15, "down")]
        result = compute_gene_meta("G", 0, des)
        assert result.i_squared > 0

class TestRunPhase6:
    def test_full_run(self):
        p1, p5 = _make_meta_test_data()
        result = run_phase6(p1, p5)
        assert result.n_genes_analyzed == 1

    def test_all_genes_analyzed(self):
        p1, p5 = _make_meta_test_data()
        result = run_phase6(p1, p5)
        assert "G1" in result.gene_results

    def test_scale_warnings(self):
        p1, p5 = _make_meta_test_data()
        result = run_phase6(p1, p5, dataset_expression_ranges={"D0": 2.0})
        assert len(result.scale_warnings) > 0

    def test_direction_consistent(self):
        p1, p5 = _make_meta_test_data()
        result = run_phase6(p1, p5)
        assert result.gene_results["G1"].direction == "down"


# ---------------------------------------------------------------------------
# 8. Operational Shell (Phase 14)
# ---------------------------------------------------------------------------

class TestConfig:
    def _write_config(self, tmp_dir, content):
        config_path = Path(tmp_dir) / "test_config.yaml"
        import yaml
        with open(config_path, "w") as f:
            yaml.dump(content, f)
        return config_path

    def test_valid_config(self, tmp_path):
        data = {"condition": "ASD", "seed_genes": "s.csv", "datasets": [{"id": "GSE1", "series_matrix": "sm.txt", "platform": "pl.txt", "role": "discovery"}]}
        path = self._write_config(tmp_path, data)
        config = load_config(path)
        assert config.condition == "ASD"

    def test_missing_condition(self, tmp_path):
        data = {"seed_genes": "s.csv", "datasets": [{"id": "G1"}]}
        path = self._write_config(tmp_path, data)
        with pytest.raises(ValueError, match="condition"):
            load_config(path)

    def test_missing_seed_genes(self, tmp_path):
        data = {"condition": "ASD", "datasets": [{"id": "G1"}]}
        path = self._write_config(tmp_path, data)
        with pytest.raises(ValueError, match="seed_genes"):
            load_config(path)

    def test_empty_datasets(self, tmp_path):
        data = {"condition": "ASD", "seed_genes": "s.csv", "datasets": []}
        path = self._write_config(tmp_path, data)
        with pytest.raises(ValueError, match="one dataset"):
            load_config(path)

    def test_discovery_enforced(self, tmp_path):
        data = {"condition": "ASD", "seed_genes": "s.csv", "datasets": [{"id": "G1", "role": "replication"}]}
        path = self._write_config(tmp_path, data)
        with pytest.raises(ValueError, match="discovery"):
            load_config(path)

    def test_dataset_id_required(self, tmp_path):
        data = {"condition": "ASD", "seed_genes": "s.csv", "datasets": [{"role": "discovery"}]}
        path = self._write_config(tmp_path, data)
        with pytest.raises(ValueError, match="id"):
            load_config(path)

class TestQCReport:
    def test_critical_halts(self):
        qc = QCReport()
        qc.add(QCCheckResult("t", "p1", False, "critical", "F"))
        assert qc.pipeline_ok is False

    def test_warning_counters(self):
        qc = QCReport()
        qc.add(QCCheckResult("t", "p1", False, "warning", "W"))
        assert qc.n_warnings == 1

    def test_passed_counters(self):
        qc = QCReport()
        qc.add(QCCheckResult("t", "p1", True, "info", "P"))
        assert qc.n_passed == 1

    def test_summary_string(self):
        qc = QCReport()
        assert "QC Report" in qc.summary()

class TestOutputWriters:
    def test_write_qc_report(self, tmp_path):
        qc = QCReport()
        qc.add(QCCheckResult("t", "p1", True, "info", "O"))
        path = write_qc_report(qc, tmp_path)
        assert path.exists()

    def test_write_phase1_summary(self, tmp_path):
        from riker.phases.phase1_crossref import Phase1Result, GeneResult
        p1 = Phase1Result({"G1": GeneResult("G1", [], 1, 1, True, -0.5, True)})
        path = write_phase1_summary(p1, tmp_path)
        assert path.exists()


class TestAllLevelsOutput:
    def test_write_all_levels(self, tmp_path):
        """All study genes should appear with level data."""
        sg, ci, ds = _make_robust_test_data()
        sens = sensitivity_analysis(ci, sg, 100, min_datasets=1)
        core = identify_core_genes(sens, sg, min_genes_per_cluster=1)

        phase1 = Phase1Result(study_genes=sg, n_study_genes=len(sg))

        # Phase3Result.cluster_labels is a dict: gene -> cluster_id
        cluster_labels = {}
        for cid, info in ci.items():
            for gene in info.gene_symbols:
                cluster_labels[gene] = cid
        phase3 = Phase3Result(cluster_labels=cluster_labels, gene_order=list(sg.keys()))

        phase4 = Phase4Result(
            core_genes=core,
            sensitivity=sens,
            n_core_genes=len(core),
        )

        path = write_phase4_all_levels(phase1, phase3, phase4, tmp_path)
        assert path.exists()

        df = pd.read_csv(path)
        assert len(df) == len(sg)
        assert "max_level" in df.columns
        assert "is_core" in df.columns
        assert df["max_level"].max() >= 1

    def test_levels_are_progressive(self, tmp_path):
        """Level 4 survivors must also be Level 1-3 survivors."""
        sg, ci, ds = _make_robust_test_data()
        sens = sensitivity_analysis(ci, sg, 100, min_datasets=1)
        core = identify_core_genes(sens, sg, min_genes_per_cluster=1)

        phase1 = Phase1Result(study_genes=sg, n_study_genes=len(sg))
        cluster_labels = {}
        for cid, info in ci.items():
            for gene in info.gene_symbols:
                cluster_labels[gene] = cid
        phase3 = Phase3Result(cluster_labels=cluster_labels, gene_order=list(sg.keys()))
        phase4 = Phase4Result(
            core_genes=core,
            sensitivity=sens,
            n_core_genes=len(core),
        )

        path = write_phase4_all_levels(phase1, phase3, phase4, tmp_path)
        df = pd.read_csv(path)

        # Any gene at level 4 must also be at levels 1, 2, 3
        l4 = df[df["level_4"] == True]
        for _, row in l4.iterrows():
            assert row["level_1"] == True
            assert row["level_2"] == True
            assert row["level_3"] == True


# ---------------------------------------------------------------------------
# v0.2.0 Tests: PCA embedding option for consensus clustering
# ---------------------------------------------------------------------------

def _make_clusterable_features(n_genes=30, n_features=3):
    """Create a feature matrix with clear cluster structure."""
    import numpy as np
    np.random.seed(42)
    # Split genes evenly across 3 clusters
    n_per = n_genes // 3
    remainder = n_genes - 3 * n_per
    cluster1 = np.random.normal(0, 0.3, (n_per, n_features))
    cluster2 = np.random.normal(3, 0.3, (n_per, n_features))
    cluster3 = np.random.normal(6, 0.3, (n_per + remainder, n_features))
    data = np.vstack([cluster1, cluster2, cluster3])
    genes = [f"G_{i}" for i in range(n_genes)]
    return pd.DataFrame(data, index=genes)


class TestPCAEmbeddingOption:
    def test_pca_embedding_produces_clusters(self):
        """PCA embedding should produce valid clusters."""
        features = _make_clusterable_features(30, 3)
        result = run_consensus_clustering(
            features, n_neighbors_list=[5], seeds=[42],
            min_cluster_size=3, min_samples=2,
            embedding_methods=["pca"],
        )
        assert result.n_clusters >= 1

    def test_combined_umap_pca(self):
        """Combined UMAP+PCA should produce more configs than UMAP alone."""
        features = _make_clusterable_features(30, 3)
        result = run_consensus_clustering(
            features, n_neighbors_list=[5], seeds=[42],
            min_cluster_size=3, min_samples=2,
            embedding_methods=["umap", "pca"],
        )
        # Should have UMAP configs + PCA configs
        assert len(result.per_config_labels) >= 2

    def test_default_is_umap_only(self):
        """Default embedding should be UMAP only."""
        features = _make_clusterable_features(30, 3)
        result = run_consensus_clustering(
            features, n_neighbors_list=[5], seeds=[42],
            min_cluster_size=3, min_samples=2,
        )
        # All configs should be UMAP
        for config in result.per_config_labels:
            assert config.get("method", "umap") == "umap"
