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


def write_phase5_verdicts(phase5_result, output_dir: Path) -> Path:
    """Write Phase 5 replication verdicts to CSV."""
    rows = []
    for gene, gv in phase5_result.gene_verdicts.items():
        rows.append({
            "gene": gene,
            "cluster_id": gv.cluster_id,
            "status": gv.status,
            "discovery_direction": gv.discovery_direction,
            "n_brain_concordant": gv.n_brain_concordant,
            "n_brain_discordant": gv.n_brain_discordant,
            "n_blood_concordant": gv.n_blood_concordant,
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
