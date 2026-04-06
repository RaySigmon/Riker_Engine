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

"""Tests for the QC checkpoint module."""

from types import SimpleNamespace

from riker.qc.checks import (
    QCCheckResult,
    QCReport,
    check_phase1,
    check_phase3,
    check_phase4,
    check_phase5,
)


# --- QCCheckResult and QCReport dataclass tests ---


class TestQCReport:
    def test_empty_report(self):
        report = QCReport()
        assert report.n_passed == 0
        assert report.n_warnings == 0
        assert report.n_critical == 0
        assert report.pipeline_ok is True
        assert report.checks == []

    def test_add_passing_check(self):
        report = QCReport()
        report.add(QCCheckResult("test", "phase1", True, "info", "ok"))
        assert report.n_passed == 1
        assert report.pipeline_ok is True

    def test_add_warning_check(self):
        report = QCReport()
        report.add(QCCheckResult("test", "phase1", False, "warning", "low yield"))
        assert report.n_warnings == 1
        assert report.pipeline_ok is True  # warnings don't halt

    def test_add_critical_check_halts_pipeline(self):
        report = QCReport()
        report.add(QCCheckResult("test", "phase1", False, "critical", "no genes"))
        assert report.n_critical == 1
        assert report.pipeline_ok is False

    def test_multiple_checks_accumulate(self):
        report = QCReport()
        report.add(QCCheckResult("a", "phase1", True, "info", "ok"))
        report.add(QCCheckResult("b", "phase1", False, "warning", "low"))
        report.add(QCCheckResult("c", "phase3", True, "info", "ok"))
        assert report.n_passed == 2
        assert report.n_warnings == 1
        assert report.n_critical == 0
        assert report.pipeline_ok is True
        assert len(report.checks) == 3

    def test_critical_after_passing_still_halts(self):
        report = QCReport()
        report.add(QCCheckResult("a", "phase1", True, "info", "ok"))
        report.add(QCCheckResult("b", "phase4", False, "critical", "zero core genes"))
        assert report.pipeline_ok is False

    def test_summary_ok(self):
        report = QCReport()
        report.add(QCCheckResult("a", "phase1", True, "info", "ok"))
        assert "OK" in report.summary()
        assert "1 passed" in report.summary()

    def test_summary_halted(self):
        report = QCReport()
        report.add(QCCheckResult("a", "phase1", False, "critical", "fail"))
        assert "HALTED" in report.summary()


# --- Phase 1 QC ---


class TestCheckPhase1:
    def _make_phase1_result(self, n_study_genes, qc_warnings=None):
        return SimpleNamespace(
            n_study_genes=n_study_genes,
            qc_warnings=qc_warnings or [],
        )

    def test_healthy_yield(self):
        result = self._make_phase1_result(n_study_genes=100)
        checks = check_phase1(result, seed_gene_count=1000)
        yield_check = checks[0]
        assert yield_check.passed is True
        assert yield_check.severity == "info"

    def test_low_yield_warning(self):
        # 3% yield: below 5% threshold, above 1%
        result = self._make_phase1_result(n_study_genes=30)
        checks = check_phase1(result, seed_gene_count=1000)
        yield_check = checks[0]
        assert yield_check.passed is False
        assert yield_check.severity == "warning"

    def test_very_low_yield_critical(self):
        # 0.5% yield: below 1% threshold
        result = self._make_phase1_result(n_study_genes=5)
        checks = check_phase1(result, seed_gene_count=1000)
        yield_check = checks[0]
        assert yield_check.passed is False
        assert yield_check.severity == "critical"

    def test_zero_seed_genes(self):
        result = self._make_phase1_result(n_study_genes=0)
        checks = check_phase1(result, seed_gene_count=0)
        yield_check = checks[0]
        assert yield_check.passed is False

    def test_fold_change_warnings_reported(self):
        result = self._make_phase1_result(
            n_study_genes=100,
            qc_warnings=["gene1 fold change out of range"],
        )
        checks = check_phase1(result, seed_gene_count=1000)
        assert len(checks) == 2
        fc_check = checks[1]
        assert fc_check.check_name == "fold_change_range"
        assert fc_check.severity == "warning"

    def test_no_fold_change_warnings(self):
        result = self._make_phase1_result(n_study_genes=100)
        checks = check_phase1(result, seed_gene_count=1000)
        assert len(checks) == 1  # only yield check


# --- Phase 3 QC ---


class TestCheckPhase3:
    def _make_phase3_result(self, n_clusters, n_noise=0):
        return SimpleNamespace(n_clusters=n_clusters, n_noise=n_noise)

    def test_healthy_clusters(self):
        checks = check_phase3(self._make_phase3_result(5, 10))
        assert checks[0].passed is True

    def test_one_cluster_warning(self):
        checks = check_phase3(self._make_phase3_result(1, 0))
        assert checks[0].passed is False
        assert checks[0].severity == "warning"

    def test_zero_clusters_critical(self):
        checks = check_phase3(self._make_phase3_result(0, 50))
        assert checks[0].passed is False
        assert checks[0].severity == "critical"


# --- Phase 4 QC ---


class TestCheckPhase4:
    def _make_phase4_result(self, n_core_genes, n_clusters_significant=1):
        return SimpleNamespace(
            n_core_genes=n_core_genes,
            n_clusters_significant=n_clusters_significant,
        )

    def test_healthy_core_genes(self):
        checks = check_phase4(self._make_phase4_result(35, 3))
        assert checks[0].passed is True

    def test_few_core_genes_warning(self):
        checks = check_phase4(self._make_phase4_result(2, 1))
        assert checks[0].passed is False
        assert checks[0].severity == "warning"

    def test_zero_core_genes_critical(self):
        checks = check_phase4(self._make_phase4_result(0, 0))
        assert checks[0].passed is False
        assert checks[0].severity == "critical"


# --- Phase 5 QC ---


class TestCheckPhase5:
    def _make_phase5_result(self, n_survived, n_eliminated=0, n_insufficient=0):
        return SimpleNamespace(
            n_survived=n_survived,
            n_eliminated=n_eliminated,
            n_insufficient=n_insufficient,
        )

    def test_healthy_replication(self):
        checks = check_phase5(self._make_phase5_result(30, 2, 3))
        assert checks[0].passed is True
        assert checks[0].severity == "info"

    def test_zero_survived_critical(self):
        checks = check_phase5(self._make_phase5_result(0, 10, 0))
        assert checks[0].passed is False
        assert checks[0].severity == "critical"

    def test_all_insufficient(self):
        checks = check_phase5(self._make_phase5_result(0, 0, 35))
        assert checks[0].passed is False
        assert checks[0].severity == "critical"


# --- Integration: QC report with phase checks ---


class TestQCIntegration:
    def test_full_pipeline_passing(self):
        report = QCReport()
        ph1 = SimpleNamespace(n_study_genes=100, qc_warnings=[])
        for c in check_phase1(ph1, 1000):
            report.add(c)
        ph3 = SimpleNamespace(n_clusters=5, n_noise=10)
        for c in check_phase3(ph3):
            report.add(c)
        ph4 = SimpleNamespace(n_core_genes=35, n_clusters_significant=3)
        for c in check_phase4(ph4):
            report.add(c)
        ph5 = SimpleNamespace(n_survived=33, n_eliminated=2, n_insufficient=0)
        for c in check_phase5(ph5):
            report.add(c)
        assert report.pipeline_ok is True
        assert report.n_passed == 4

    def test_pipeline_halts_on_phase4_critical(self):
        report = QCReport()
        ph1 = SimpleNamespace(n_study_genes=100, qc_warnings=[])
        for c in check_phase1(ph1, 1000):
            report.add(c)
        ph3 = SimpleNamespace(n_clusters=5, n_noise=10)
        for c in check_phase3(ph3):
            report.add(c)
        ph4 = SimpleNamespace(n_core_genes=0, n_clusters_significant=0)
        for c in check_phase4(ph4):
            report.add(c)
        assert report.pipeline_ok is False
        assert report.n_critical == 1
