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
Riker Engine - QC Checkpoint Runner.

Runs quality control checks at each phase boundary. Critical failures
halt the pipeline; warnings are logged and collected.

References:
    Blueprint Section 12 (QC Framework)
"""

import logging
import warnings
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


@dataclass
class QCCheckResult:
    """Result of a single QC check.

    Attributes:
        check_name: Name of the check.
        phase: Which phase boundary this runs at.
        passed: True if check passed.
        severity: 'info', 'warning', or 'critical'.
        message: Human-readable description.
    """
    check_name: str
    phase: str
    passed: bool
    severity: str
    message: str


@dataclass
class QCReport:
    """Accumulated QC results across all phases.

    Attributes:
        checks: List of all QC check results.
        n_passed: Number of checks that passed.
        n_warnings: Number of warning-level issues.
        n_critical: Number of critical failures.
        pipeline_ok: True if no critical failures.
    """
    checks: list = field(default_factory=list)
    n_passed: int = 0
    n_warnings: int = 0
    n_critical: int = 0
    pipeline_ok: bool = True

    def add(self, result: QCCheckResult) -> None:
        """Add a check result and update counters."""
        self.checks.append(result)
        if result.passed:
            self.n_passed += 1
        elif result.severity == "critical":
            self.n_critical += 1
            self.pipeline_ok = False
            logger.error(f"QC CRITICAL [{result.phase}] {result.check_name}: {result.message}")
        elif result.severity == "warning":
            self.n_warnings += 1
            logger.warning(f"QC WARNING [{result.phase}] {result.check_name}: {result.message}")
        else:
            logger.info(f"QC INFO [{result.phase}] {result.check_name}: {result.message}")

    def summary(self) -> str:
        """Return a human-readable summary."""
        return (
            f"QC Report: {self.n_passed} passed, "
            f"{self.n_warnings} warnings, {self.n_critical} critical. "
            f"Pipeline {'OK' if self.pipeline_ok else 'HALTED'}."
        )


def check_phase1(phase1_result, seed_gene_count: int) -> list[QCCheckResult]:
    """QC checks after Phase 1."""
    results = []

    # Check study gene yield
    yield_pct = 100 * phase1_result.n_study_genes / seed_gene_count if seed_gene_count > 0 else 0
    results.append(QCCheckResult(
        check_name="study_gene_yield",
        phase="phase1",
        passed=yield_pct >= 5.0,
        severity="critical" if yield_pct < 1.0 else ("warning" if yield_pct < 5.0 else "info"),
        message=f"{phase1_result.n_study_genes}/{seed_gene_count} seed genes passed "
                f"({yield_pct:.1f}%). Minimum recommended: 5%.",
    ))

    # Check for QC warnings from Phase 1
    if phase1_result.qc_warnings:
        results.append(QCCheckResult(
            check_name="fold_change_range",
            phase="phase1",
            passed=False,
            severity="warning",
            message=f"{len(phase1_result.qc_warnings)} fold change warning(s). "
                    f"Check normalizer output.",
        ))

    return results


def check_phase3(phase3_result) -> list[QCCheckResult]:
    """QC checks after Phase 3."""
    results = []

    results.append(QCCheckResult(
        check_name="cluster_count",
        phase="phase3",
        passed=phase3_result.n_clusters >= 2,
        severity="critical" if phase3_result.n_clusters == 0 else (
            "warning" if phase3_result.n_clusters < 2 else "info"
        ),
        message=f"{phase3_result.n_clusters} clusters found, "
                f"{phase3_result.n_noise} noise genes.",
    ))

    return results


def check_phase4(phase4_result) -> list[QCCheckResult]:
    """QC checks after Phase 4."""
    results = []

    results.append(QCCheckResult(
        check_name="core_gene_count",
        phase="phase4",
        passed=phase4_result.n_core_genes >= 3,
        severity="critical" if phase4_result.n_core_genes == 0 else (
            "warning" if phase4_result.n_core_genes < 3 else "info"
        ),
        message=f"{phase4_result.n_core_genes} core genes identified. "
                f"{phase4_result.n_clusters_significant} significant clusters.",
    ))

    return results


def check_phase5(phase5_result) -> list[QCCheckResult]:
    """QC checks after Phase 5."""
    results = []

    total = phase5_result.n_survived + phase5_result.n_eliminated + phase5_result.n_insufficient
    results.append(QCCheckResult(
        check_name="replication_survival",
        phase="phase5",
        passed=phase5_result.n_survived >= 1,
        severity="critical" if phase5_result.n_survived == 0 else "info",
        message=f"{phase5_result.n_survived} survived, "
                f"{phase5_result.n_eliminated} eliminated, "
                f"{phase5_result.n_insufficient} insufficient data "
                f"(of {total} total).",
    ))

    return results
