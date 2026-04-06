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
Riker Engine - Output Writer.

Writes pipeline results to structured JSON and CSV files.
All outputs go to the configured output directory.

References:
    Blueprint Section 15 (Output Specification)
"""

import json
import logging
from dataclasses import asdict
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def write_phase1_summary(phase1_result, output_dir: Path) -> Path:
    """Write Phase 1 study gene summary to CSV."""
    rows = []
    for gene, gr in phase1_result.study_genes.items():
        rows.append({
            "gene": gene,
            "n_datasets_detected": gr.n_datasets_detected,
            "n_datasets_significant": gr.n_datasets_significant,
            "mean_log2fc": gr.mean_log2fc,
            "consistent_direction": gr.consistent_direction,
        })

    df = pd.DataFrame(rows)
    path = output_dir / "phase1_study_genes.csv"
    df.to_csv(path, index=False)
    logger.info(f"Phase 1 summary: {path}")
    return path


def write_phase4_core_genes(phase4_result, output_dir: Path) -> Path:
    """Write Phase 4 core gene list to CSV."""
    rows = []
    for gene, cg in phase4_result.core_genes.items():
        rows.append({
            "gene": gene,
            "cluster_id": cg.cluster_id,
            "max_level_survived": cg.max_level_survived,
            "mean_log2fc": cg.mean_log2fc,
            "direction": cg.direction,
        })

    df = pd.DataFrame(rows)
    path = output_dir / "phase4_core_genes.csv"
    df.to_csv(path, index=False)
    logger.info(f"Phase 4 core genes: {path}")
    return path


def write_phase4_all_levels(
    phase1_result,
    phase3_result,
    phase4_result,
    output_dir: Path,
) -> Path:
    """Write progressive confidence levels for all study genes.

    Shows every gene from Phase 1, its cluster assignment, and
    which sensitivity levels it survived. Gives researchers a
    ranked confidence list beyond the binary core/not-core output.
    """
    rows = []

    # Collect all genes that survived each level across all clusters
    level_genes = {1: set(), 2: set(), 3: set(), 4: set()}
    for cluster_id, cluster_sens in phase4_result.sensitivity.items():
        for level, survivors in cluster_sens.level_survivors.items():
            level_genes[level].update(survivors)

    # Get core gene set
    core_set = set(phase4_result.core_genes.keys())

    # Build rows for ALL study genes
    for gene, gene_result in phase1_result.study_genes.items():
        cluster_id = phase3_result.cluster_labels.get(gene, -1)
        l1 = gene in level_genes[1]
        l2 = gene in level_genes[2]
        l3 = gene in level_genes[3]
        l4 = gene in level_genes[4]

        if l4:
            max_level = 4
        elif l3:
            max_level = 3
        elif l2:
            max_level = 2
        elif l1:
            max_level = 1
        else:
            max_level = 0

        rows.append({
            "gene": gene,
            "cluster_id": cluster_id,
            "level_1": l1,
            "level_2": l2,
            "level_3": l3,
            "level_4": l4,
            "max_level": max_level,
            "is_core": gene in core_set,
            "mean_log2fc": gene_result.mean_log2fc,
            "direction": "down" if gene_result.mean_log2fc < 0 else "up",
        })

    df = pd.DataFrame(rows)
    df = df.sort_values(
        ["max_level", "is_core"],
        ascending=[False, False],
    )

    path = output_dir / "phase4_all_levels.csv"
    df.to_csv(path, index=False)
    logger.info(f"Phase 4 all levels: {path} ({len(df)} genes)")
    return path


def write_phase5_verdicts(phase5_result, output_dir: Path) -> Path:
    """Write Phase 5 replication verdicts to CSV."""
    rows = []
    for gene, gv in phase5_result.gene_verdicts.items():
        rows.append({
            "gene": gene,
            "cluster_id": gv.cluster_id,
            "status": gv.status,
            "discovery_direction": gv.discovery_direction,
            "n_same_tissue_concordant": gv.n_same_tissue_concordant,
            "n_same_tissue_discordant": gv.n_same_tissue_discordant,
            "n_cross_tissue_concordant": gv.n_cross_tissue_concordant,
            "reason": gv.reason,
        })

    df = pd.DataFrame(rows)
    path = output_dir / "phase5_verdicts.csv"
    df.to_csv(path, index=False)
    logger.info(f"Phase 5 verdicts: {path}")
    return path


def write_phase6_meta(phase6_result, output_dir: Path) -> Path:
    """Write Phase 6 meta-analysis results to CSV."""
    rows = []
    for gene, meta in phase6_result.gene_results.items():
        rows.append({
            "gene": gene,
            "cluster_id": meta.cluster_id,
            "random_effect": meta.random_effect,
            "random_se": meta.random_se,
            "random_p": meta.random_p,
            "fixed_effect": meta.fixed_effect,
            "fixed_se": meta.fixed_se,
            "fixed_p": meta.fixed_p,
            "i_squared": meta.i_squared,
            "tau_squared": meta.tau_squared,
            "cochran_q": meta.cochran_q,
            "n_datasets": meta.n_datasets,
            "direction": meta.direction,
        })

    df = pd.DataFrame(rows)
    path = output_dir / "phase6_meta_analysis.csv"
    df.to_csv(path, index=False)
    logger.info(f"Phase 6 meta-analysis: {path}")
    return path


def write_qc_report(qc_report, output_dir: Path) -> Path:
    """Write QC report to JSON."""
    data = {
        "summary": qc_report.summary(),
        "n_passed": qc_report.n_passed,
        "n_warnings": qc_report.n_warnings,
        "n_critical": qc_report.n_critical,
        "pipeline_ok": qc_report.pipeline_ok,
        "checks": [
            {
                "check_name": c.check_name,
                "phase": c.phase,
                "passed": c.passed,
                "severity": c.severity,
                "message": c.message,
            }
            for c in qc_report.checks
        ],
    }

    path = output_dir / "qc_report.json"
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info(f"QC report: {path}")
    return path


def write_pipeline_summary(
    config,
    phase1_result,
    phase4_result,
    phase5_result,
    phase6_result,
    qc_report,
    output_dir: Path,
) -> Path:
    """Write overall pipeline summary to JSON."""
    data = {
        "condition": config.condition,
        "seed_genes": config.seed_genes_path,
        "n_datasets_discovery": len([d for d in config.datasets if d.role == "discovery"]),
        "n_datasets_replication": len([d for d in config.datasets if d.role == "replication"]),
        "phase1_study_genes": phase1_result.n_study_genes,
        "phase4_core_genes": phase4_result.n_core_genes,
        "phase4_significant_clusters": phase4_result.n_clusters_significant,
        "phase5_survived": phase5_result.n_survived,
        "phase5_eliminated": phase5_result.n_eliminated,
        "phase6_genes_analyzed": phase6_result.n_genes_analyzed,
        "phase6_significant_random": phase6_result.n_significant_random,
        "qc_status": "PASSED" if qc_report.pipeline_ok else "FAILED",
        "locked_core_genes": phase5_result.locked_core_genes,
    }

    path = output_dir / "pipeline_summary.json"
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    logger.info(f"Pipeline summary: {path}")
    return path
