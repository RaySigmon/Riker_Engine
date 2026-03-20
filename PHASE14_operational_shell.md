# INSTRUCTION SET FOR KAI — PHASE 14: OPERATIONAL SHELL
# config.py + cli.py + qc/checks.py + io/outputs.py

## References
- Blueprint Section 13 (YAML Configuration)
- Blueprint Section 12 (QC Framework)
- Blueprint Section 15 (Output Specification)
- Blueprint Section 4 (Engine Architecture — pipeline flow)

## WHY THIS MATTERS

This is the operational layer that makes the engine runnable. After this,
a user can write a YAML config file and run `riker run config.yaml` from
the command line. The engine reads the config, loads data, runs all 6
phases in sequence, and writes structured output.

## BUILD ORDER (4 files)

1. `riker/config.py` — YAML config loader and validation
2. `riker/qc/checks.py` — QC checkpoint runner
3. `riker/io/outputs.py` — Output writer (JSON + CSV)
4. `riker/cli.py` — CLI entry point (`riker run config.yaml`)

## CRITICAL REQUIREMENTS

- Config must support per-dataset overrides for phenotype extraction
- Config must separate discovery vs replication datasets with tissue labels
- QC checks must run at each phase boundary and halt on critical failures
- Output must be machine-readable (JSON) with human-readable summaries (CSV)
- CLI must handle errors gracefully and report which phase failed

---

## FILE 1: `riker/config.py`

Write at `/home/kai001/riker-engine/riker/config.py`:

```python
"""
Riker Engine - YAML Configuration Loader.

Parses and validates the pipeline configuration file.
Supports per-dataset phenotype overrides and separate
discovery/replication dataset specifications.

References:
    Blueprint Section 13 (YAML Configuration)
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)


@dataclass
class DatasetConfig:
    """Configuration for a single dataset.

    Attributes:
        dataset_id: Unique identifier (e.g., 'GSE28521').
        series_matrix_path: Path to series matrix file.
        platform_path: Path to platform annotation file.
        role: 'discovery' or 'replication'.
        tissue: 'brain' or 'blood' (for replication datasets).
        phenotype_field: Override metadata field for case/control.
        case_values: Override values indicating case samples.
        control_values: Override values indicating control samples.
    """
    dataset_id: str = ""
    series_matrix_path: str = ""
    platform_path: str = ""
    role: str = "discovery"
    tissue: str = "brain"
    phenotype_field: str | None = None
    case_values: list | None = None
    control_values: list | None = None


@dataclass
class PipelineConfig:
    """Complete pipeline configuration.

    Attributes:
        condition: Condition being studied (e.g., 'ASD', 'Alzheimers').
        seed_genes_path: Path to seed gene CSV file.
        hgnc_path: Path to HGNC complete set (or 'auto' for download).
        datasets: List of DatasetConfig.
        output_dir: Directory for output files.
        pathways: Dict with pathway database paths/settings.
        phase1: Phase 1 parameters.
        phase3: Phase 3 parameters.
        phase4: Phase 4 parameters.
        n_permutations: Number of permutations for significance tests.
        random_seed: Base random seed for reproducibility.
    """
    condition: str = "unknown"
    seed_genes_path: str = ""
    hgnc_path: str = "auto"
    datasets: list = field(default_factory=list)
    output_dir: str = "riker_output"
    pathways: dict = field(default_factory=dict)

    # Phase parameters
    phase1_p_threshold: float = 0.05
    phase1_min_datasets: int = 2
    phase3_n_neighbors: list = field(default_factory=lambda: [10, 15, 30])
    phase3_seeds: list = field(default_factory=lambda: [42, 123, 456, 789, 1024])
    phase3_min_cluster_size: int = 5
    phase3_min_samples: int = 3
    phase4_n_permutations: int = 10000
    phase4_permutation_seed: int = 42
    random_seed: int = 42


def load_config(path: str | Path) -> PipelineConfig:
    """Load and validate a YAML configuration file.

    Parameters
    ----------
    path : str or Path
        Path to the YAML config file.

    Returns
    -------
    PipelineConfig

    Raises
    ------
    FileNotFoundError
        If config file doesn't exist.
    ValueError
        If required fields are missing or invalid.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with open(path, "r") as f:
        raw = yaml.safe_load(f)

    if not isinstance(raw, dict):
        raise ValueError("Config file must be a YAML dictionary.")

    config = PipelineConfig()

    # Required fields
    if "condition" not in raw:
        raise ValueError("Config missing required field: 'condition'")
    config.condition = raw["condition"]

    if "seed_genes" not in raw:
        raise ValueError("Config missing required field: 'seed_genes'")
    config.seed_genes_path = raw["seed_genes"]

    config.hgnc_path = raw.get("hgnc_path", "auto")
    config.output_dir = raw.get("output_dir", "riker_output")
    config.random_seed = raw.get("random_seed", 42)

    # Datasets
    if "datasets" not in raw or not raw["datasets"]:
        raise ValueError("Config must specify at least one dataset.")

    for ds_raw in raw["datasets"]:
        ds = DatasetConfig(
            dataset_id=ds_raw.get("id", ""),
            series_matrix_path=ds_raw.get("series_matrix", ""),
            platform_path=ds_raw.get("platform", ""),
            role=ds_raw.get("role", "discovery"),
            tissue=ds_raw.get("tissue", "brain"),
            phenotype_field=ds_raw.get("phenotype_field"),
            case_values=ds_raw.get("case_values"),
            control_values=ds_raw.get("control_values"),
        )
        if not ds.dataset_id:
            raise ValueError("Each dataset must have an 'id' field.")
        config.datasets.append(ds)

    # Validate at least one discovery dataset
    discovery = [d for d in config.datasets if d.role == "discovery"]
    if not discovery:
        raise ValueError("Config must include at least one discovery dataset.")

    # Phase parameters
    p1 = raw.get("phase1", {})
    config.phase1_p_threshold = p1.get("p_threshold", 0.05)
    config.phase1_min_datasets = p1.get("min_datasets", 2)

    p3 = raw.get("phase3", {})
    config.phase3_n_neighbors = p3.get("n_neighbors", [10, 15, 30])
    config.phase3_seeds = p3.get("seeds", [42, 123, 456, 789, 1024])
    config.phase3_min_cluster_size = p3.get("min_cluster_size", 5)
    config.phase3_min_samples = p3.get("min_samples", 3)

    p4 = raw.get("phase4", {})
    config.phase4_n_permutations = p4.get("n_permutations", 10000)
    config.phase4_permutation_seed = p4.get("seed", 42)

    # Pathways
    config.pathways = raw.get("pathways", {})

    logger.info(
        f"Config loaded: condition={config.condition}, "
        f"{len(discovery)} discovery + "
        f"{len(config.datasets) - len(discovery)} replication datasets."
    )

    return config
```

---

## FILE 2: `riker/qc/checks.py`

Write at `/home/kai001/riker-engine/riker/qc/checks.py`:

```python
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
```

---

## FILE 3: `riker/io/outputs.py`

Write at `/home/kai001/riker-engine/riker/io/outputs.py`:

```python
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
```

---

## FILE 4: `riker/cli.py`

Write at `/home/kai001/riker-engine/riker/cli.py`:

```python
"""
Riker Engine - Command Line Interface.

Usage:
    riker run config.yaml
    riker validate config.yaml

References:
    Blueprint Section 4 (Engine Architecture — pipeline flow)
"""

import argparse
import logging
import sys
from pathlib import Path

from riker.config import load_config

logger = logging.getLogger("riker")


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for the pipeline."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def cmd_validate(args) -> int:
    """Validate a config file without running the pipeline."""
    try:
        config = load_config(args.config)
        print(f"Config valid: {config.condition}")
        print(f"  Seed genes: {config.seed_genes_path}")
        print(f"  Datasets: {len(config.datasets)}")
        for ds in config.datasets:
            print(f"    {ds.dataset_id} ({ds.role}, {ds.tissue})")
        print(f"  Output dir: {config.output_dir}")
        return 0
    except (FileNotFoundError, ValueError) as e:
        print(f"Config error: {e}", file=sys.stderr)
        return 1


def cmd_run(args) -> int:
    """Run the full pipeline."""
    try:
        config = load_config(args.config)
    except (FileNotFoundError, ValueError) as e:
        logger.error(f"Config error: {e}")
        return 1

    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    from riker.qc.checks import QCReport, check_phase1, check_phase3, check_phase4, check_phase5

    qc = QCReport()

    # ---- Phase 1: Cross-Referencing ----
    logger.info("=" * 60)
    logger.info("PHASE 1: Cross-Referencing")
    logger.info("=" * 60)

    try:
        from riker.ingestion.gene_db import SeedGeneDB, HGNCResolver
        from riker.ingestion.geo_parser import GEOSeriesMatrix, PhenotypeExtractor, ProbeGeneMapper
        from riker.ingestion.normalizer import normalize_expression
        from riker.phases.phase1_crossref import run_phase1

        # Load seed genes
        if config.hgnc_path == "auto":
            resolver = HGNCResolver()
        else:
            resolver = HGNCResolver(hgnc_path=config.hgnc_path)

        seed_db = SeedGeneDB(
            csv_path=config.seed_genes_path,
            resolver=resolver,
        )
        seed_genes = [g.resolved_symbol for g in seed_db.genes]
        seed_gene_count = len(seed_genes)
        logger.info(f"Loaded {seed_gene_count} seed genes.")

        # Parse discovery datasets
        discovery_datasets = {}
        discovery_phenotypes = {}
        dataset_ids = []

        for ds_config in config.datasets:
            if ds_config.role != "discovery":
                continue

            ds_id = ds_config.dataset_id
            dataset_ids.append(ds_id)
            logger.info(f"Loading dataset: {ds_id}")

            # Parse series matrix
            geo = GEOSeriesMatrix(ds_config.series_matrix_path)
            geo_result = geo.get_result()

            # Extract phenotypes
            if ds_config.phenotype_field:
                extractor = PhenotypeExtractor(
                    override_field=ds_config.phenotype_field,
                    override_case_values=ds_config.case_values,
                    override_control_values=ds_config.control_values,
                )
            else:
                extractor = PhenotypeExtractor()

            pheno = extractor.extract(
                geo_result.sample_metadata, geo_result.sample_ids
            )
            discovery_phenotypes[ds_id] = pheno.groups

            # Map probes to genes
            mapper = ProbeGeneMapper(ds_config.platform_path, resolver=resolver)
            gene_expr = mapper.map_expression(geo_result.expression)

            # Normalize
            norm = normalize_expression(gene_expr)
            if norm.was_transformed:
                gene_expr = norm.data

            discovery_datasets[ds_id] = gene_expr

        # Run Phase 1
        phase1_result = run_phase1(
            seed_genes, discovery_datasets, discovery_phenotypes,
            p_threshold=config.phase1_p_threshold,
            min_datasets=config.phase1_min_datasets,
        )

        # QC
        for check in check_phase1(phase1_result, seed_gene_count):
            qc.add(check)
        if not qc.pipeline_ok:
            logger.error("Pipeline halted after Phase 1 QC.")
            from riker.io.outputs import write_qc_report
            write_qc_report(qc, output_dir)
            return 1

        from riker.io.outputs import write_phase1_summary
        write_phase1_summary(phase1_result, output_dir)

    except Exception as e:
        logger.error(f"Phase 1 failed: {e}", exc_info=True)
        return 1

    # ---- Phase 2: Pathway Mapping ----
    logger.info("=" * 60)
    logger.info("PHASE 2: Pathway Mapping")
    logger.info("=" * 60)

    try:
        from riker.phases.phase2_pathways import run_phase2, load_pathways_from_dict

        # Load pathways (from config or empty)
        if config.pathways and "data" in config.pathways:
            pathway_db = load_pathways_from_dict(
                config.pathways["data"],
                names=config.pathways.get("names"),
                source=config.pathways.get("source", "custom"),
            )
        else:
            # Empty pathway DB — features will be expression-only
            pathway_db = load_pathways_from_dict({}, source="none")
            logger.warning("No pathway data configured. Using expression features only.")

        phase2_result = run_phase2(phase1_result, pathway_db)

    except Exception as e:
        logger.error(f"Phase 2 failed: {e}", exc_info=True)
        return 1

    # ---- Phase 3: Consensus Clustering ----
    logger.info("=" * 60)
    logger.info("PHASE 3: Consensus Clustering")
    logger.info("=" * 60)

    try:
        from riker.phases.phase3_clustering import run_consensus_clustering

        phase3_result = run_consensus_clustering(
            phase2_result.feature_matrix,
            n_neighbors_list=config.phase3_n_neighbors,
            seeds=config.phase3_seeds,
            min_cluster_size=config.phase3_min_cluster_size,
            min_samples=config.phase3_min_samples,
        )

        for check in check_phase3(phase3_result):
            qc.add(check)
        if not qc.pipeline_ok:
            logger.error("Pipeline halted after Phase 3 QC.")
            from riker.io.outputs import write_qc_report
            write_qc_report(qc, output_dir)
            return 1

    except Exception as e:
        logger.error(f"Phase 3 failed: {e}", exc_info=True)
        return 1

    # ---- Phase 4: Robustness Testing ----
    logger.info("=" * 60)
    logger.info("PHASE 4: Robustness Testing")
    logger.info("=" * 60)

    try:
        from riker.phases.phase4_robustness import run_phase4

        phase4_result = run_phase4(
            phase1_result, phase3_result,
            seed_gene_count=seed_gene_count,
            dataset_ids=dataset_ids,
            n_permutations=config.phase4_n_permutations,
            permutation_seed=config.phase4_permutation_seed,
        )

        for check in check_phase4(phase4_result):
            qc.add(check)
        if not qc.pipeline_ok:
            logger.error("Pipeline halted after Phase 4 QC.")
            from riker.io.outputs import write_qc_report
            write_qc_report(qc, output_dir)
            return 1

        from riker.io.outputs import write_phase4_core_genes
        write_phase4_core_genes(phase4_result, output_dir)

    except Exception as e:
        logger.error(f"Phase 4 failed: {e}", exc_info=True)
        return 1

    # ---- Phase 5: Independent Replication ----
    logger.info("=" * 60)
    logger.info("PHASE 5: Independent Replication")
    logger.info("=" * 60)

    try:
        from riker.phases.phase5_replication import run_phase5

        # Load replication datasets
        replication_datasets = {}
        replication_phenotypes = {}
        dataset_tissues = {}

        for ds_config in config.datasets:
            if ds_config.role != "replication":
                continue

            ds_id = ds_config.dataset_id
            logger.info(f"Loading replication dataset: {ds_id}")

            geo = GEOSeriesMatrix(ds_config.series_matrix_path)
            geo_result = geo.get_result()

            if ds_config.phenotype_field:
                extractor = PhenotypeExtractor(
                    override_field=ds_config.phenotype_field,
                    override_case_values=ds_config.case_values,
                    override_control_values=ds_config.control_values,
                )
            else:
                extractor = PhenotypeExtractor()

            pheno = extractor.extract(
                geo_result.sample_metadata, geo_result.sample_ids
            )

            mapper = ProbeGeneMapper(ds_config.platform_path, resolver=resolver)
            gene_expr = mapper.map_expression(geo_result.expression)

            replication_datasets[ds_id] = gene_expr
            replication_phenotypes[ds_id] = pheno.groups
            dataset_tissues[ds_id] = ds_config.tissue

        if replication_datasets:
            phase5_result = run_phase5(
                phase4_result.core_genes,
                replication_datasets,
                replication_phenotypes,
                dataset_tissues,
            )
        else:
            logger.warning("No replication datasets configured. Skipping Phase 5.")
            from riker.phases.phase5_replication import Phase5Result
            phase5_result = Phase5Result(
                locked_core_genes=sorted(phase4_result.core_genes.keys()),
                n_survived=len(phase4_result.core_genes),
            )

        for check in check_phase5(phase5_result):
            qc.add(check)

        from riker.io.outputs import write_phase5_verdicts
        write_phase5_verdicts(phase5_result, output_dir)

    except Exception as e:
        logger.error(f"Phase 5 failed: {e}", exc_info=True)
        return 1

    # ---- Phase 6: Meta-Analysis ----
    logger.info("=" * 60)
    logger.info("PHASE 6: Effect Size Meta-Analysis")
    logger.info("=" * 60)

    try:
        from riker.phases.phase6_meta import run_phase6

        phase6_result = run_phase6(phase1_result, phase5_result)

        from riker.io.outputs import write_phase6_meta
        write_phase6_meta(phase6_result, output_dir)

    except Exception as e:
        logger.error(f"Phase 6 failed: {e}", exc_info=True)
        return 1

    # ---- Write Final Outputs ----
    logger.info("=" * 60)
    logger.info("WRITING FINAL OUTPUTS")
    logger.info("=" * 60)

    from riker.io.outputs import write_qc_report, write_pipeline_summary
    write_qc_report(qc, output_dir)
    write_pipeline_summary(
        config, phase1_result, phase4_result,
        phase5_result, phase6_result, qc, output_dir,
    )

    logger.info("=" * 60)
    logger.info(f"PIPELINE COMPLETE: {config.condition}")
    logger.info(qc.summary())
    logger.info(f"Output: {output_dir}")
    logger.info("=" * 60)

    return 0


def main() -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="riker",
        description="Riker Engine — Condition-agnostic transcriptomics pipeline",
    )
    subparsers = parser.add_subparsers(dest="command")

    # run command
    run_parser = subparsers.add_parser("run", help="Run the full pipeline")
    run_parser.add_argument("config", help="Path to YAML config file")
    run_parser.add_argument("-v", "--verbose", action="store_true",
                            help="Enable debug logging")

    # validate command
    val_parser = subparsers.add_parser("validate", help="Validate config file")
    val_parser.add_argument("config", help="Path to YAML config file")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return 1

    setup_logging(verbose=getattr(args, "verbose", False))

    if args.command == "validate":
        return cmd_validate(args)
    elif args.command == "run":
        return cmd_run(args)

    return 1


if __name__ == "__main__":
    sys.exit(main())
```

---

## TESTS: APPEND to `tests/test_phases.py`

**APPEND** the following to the END of `tests/test_phases.py`:

```python


# ===========================================================================
# PHASE 14: OPERATIONAL SHELL TESTS
# ===========================================================================

import tempfile
import json

from riker.config import load_config, PipelineConfig, DatasetConfig
from riker.qc.checks import QCReport, QCCheckResult, check_phase1
from riker.io.outputs import write_phase1_summary, write_qc_report


# ---------------------------------------------------------------------------
# 20. Config loading
# ---------------------------------------------------------------------------

class TestConfig:
    """Test YAML config loading and validation."""

    def _write_config(self, tmp_dir, content):
        """Write a YAML config to a temp file."""
        config_path = Path(tmp_dir) / "test_config.yaml"
        import yaml
        with open(config_path, "w") as f:
            yaml.dump(content, f)
        return config_path

    def test_valid_config(self, tmp_path):
        config_data = {
            "condition": "ASD",
            "seed_genes": "/path/to/seeds.csv",
            "datasets": [
                {"id": "GSE28521", "series_matrix": "/path/to/sm.txt",
                 "platform": "/path/to/pl.txt", "role": "discovery"},
            ],
            "output_dir": str(tmp_path / "output"),
        }
        config_path = self._write_config(tmp_path, config_data)
        config = load_config(config_path)

        assert config.condition == "ASD"
        assert len(config.datasets) == 1
        assert config.datasets[0].dataset_id == "GSE28521"

    def test_missing_condition(self, tmp_path):
        config_data = {
            "seed_genes": "/path/to/seeds.csv",
            "datasets": [{"id": "GSE1", "role": "discovery"}],
        }
        config_path = self._write_config(tmp_path, config_data)
        with pytest.raises(ValueError, match="condition"):
            load_config(config_path)

    def test_no_discovery_dataset(self, tmp_path):
        config_data = {
            "condition": "ASD",
            "seed_genes": "/path/to/seeds.csv",
            "datasets": [
                {"id": "REP1", "role": "replication", "tissue": "brain"},
            ],
        }
        config_path = self._write_config(tmp_path, config_data)
        with pytest.raises(ValueError, match="discovery"):
            load_config(config_path)

    def test_phenotype_overrides(self, tmp_path):
        config_data = {
            "condition": "AD",
            "seed_genes": "/path/to/seeds.csv",
            "datasets": [
                {"id": "GSE33000", "series_matrix": "/path/sm.txt",
                 "platform": "/path/pl.txt", "role": "discovery",
                 "phenotype_field": "Sample_characteristics_ch1",
                 "case_values": ["alzheimer"], "control_values": ["control"]},
            ],
        }
        config_path = self._write_config(tmp_path, config_data)
        config = load_config(config_path)
        ds = config.datasets[0]
        assert ds.phenotype_field == "Sample_characteristics_ch1"
        assert ds.case_values == ["alzheimer"]

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_config("/nonexistent/config.yaml")

    def test_default_parameters(self, tmp_path):
        config_data = {
            "condition": "TEST",
            "seed_genes": "/path/seeds.csv",
            "datasets": [{"id": "DS1", "role": "discovery"}],
        }
        config_path = self._write_config(tmp_path, config_data)
        config = load_config(config_path)

        assert config.phase1_p_threshold == 0.05
        assert config.phase1_min_datasets == 2
        assert config.phase3_min_cluster_size == 5
        assert config.phase4_n_permutations == 10000


# ---------------------------------------------------------------------------
# 21. QC checks
# ---------------------------------------------------------------------------

class TestQCReport:
    """Test QC checkpoint system."""

    def test_add_passed(self):
        qc = QCReport()
        qc.add(QCCheckResult("test", "p1", True, "info", "OK"))
        assert qc.n_passed == 1
        assert qc.pipeline_ok is True

    def test_critical_halts(self):
        qc = QCReport()
        qc.add(QCCheckResult("test", "p1", False, "critical", "FAIL"))
        assert qc.n_critical == 1
        assert qc.pipeline_ok is False

    def test_warning_does_not_halt(self):
        qc = QCReport()
        qc.add(QCCheckResult("test", "p1", False, "warning", "WARN"))
        assert qc.n_warnings == 1
        assert qc.pipeline_ok is True

    def test_summary(self):
        qc = QCReport()
        qc.add(QCCheckResult("a", "p1", True, "info", "OK"))
        qc.add(QCCheckResult("b", "p1", False, "warning", "W"))
        summary = qc.summary()
        assert "1 passed" in summary
        assert "1 warnings" in summary


# ---------------------------------------------------------------------------
# 22. Output writers
# ---------------------------------------------------------------------------

class TestOutputWriters:
    """Test output file writing."""

    def test_write_qc_report(self, tmp_path):
        qc = QCReport()
        qc.add(QCCheckResult("test_check", "phase1", True, "info", "OK"))
        path = write_qc_report(qc, tmp_path)

        assert path.exists()
        with open(path) as f:
            data = json.load(f)
        assert data["n_passed"] == 1
        assert data["pipeline_ok"] is True

    def test_write_phase1_summary(self, tmp_path):
        from riker.phases.phase1_crossref import Phase1Result, GeneResult
        phase1 = Phase1Result(
            study_genes={
                "GENE_A": GeneResult(
                    gene="GENE_A", de_results=[],
                    n_datasets_detected=3, n_datasets_significant=3,
                    passes_filter=True, mean_log2fc=-0.5,
                    consistent_direction=True,
                ),
            },
            n_study_genes=1,
        )
        path = write_phase1_summary(phase1, tmp_path)
        assert path.exists()
        df = pd.read_csv(path)
        assert len(df) == 1
        assert df.iloc[0]["gene"] == "GENE_A"
```

---

## ALSO: Ensure `riker/qc/__init__.py` and `riker/io/__init__.py` exist

Create empty init files:

```bash
touch /home/kai001/riker-engine/riker/qc/__init__.py
touch /home/kai001/riker-engine/riker/io/__init__.py
```

---

## EXECUTION INSTRUCTIONS

After writing all 4 files, creating init files, and appending tests, run:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_phases.py -v 2>&1
```

Then confirmation — validate a test config:

```bash
cd /home/kai001/riker-engine && python -c "
import tempfile, yaml, json
from pathlib import Path
from riker.config import load_config
from riker.cli import cmd_validate

# Create a test config
config_data = {
    'condition': 'ASD',
    'seed_genes': '/path/to/asd_seed_genes.csv',
    'hgnc_path': 'auto',
    'datasets': [
        {'id': 'GSE28521', 'series_matrix': '/data/GSE28521_sm.txt',
         'platform': '/data/GPL570.txt', 'role': 'discovery'},
        {'id': 'GSE38322', 'series_matrix': '/data/GSE38322_sm.txt',
         'platform': '/data/GPL6244.txt', 'role': 'discovery'},
        {'id': 'GSE64018', 'series_matrix': '/data/GSE64018_sm.txt',
         'platform': '/data/GPL11154.txt', 'role': 'replication',
         'tissue': 'brain'},
    ],
    'output_dir': '/tmp/riker_test_output',
    'phase1': {'p_threshold': 0.05, 'min_datasets': 2},
    'phase3': {'n_neighbors': [10, 15, 30], 'seeds': [42, 123, 456]},
    'phase4': {'n_permutations': 5000},
}

with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
    yaml.dump(config_data, f)
    config_path = f.name

config = load_config(config_path)
print(f'=== Config Validation ===')
print(f'Condition: {config.condition}')
print(f'Seed genes: {config.seed_genes_path}')
print(f'Datasets: {len(config.datasets)}')
for ds in config.datasets:
    print(f'  {ds.dataset_id}: {ds.role} ({ds.tissue})')
print(f'Phase 1 threshold: {config.phase1_p_threshold}')
print(f'Phase 3 neighbors: {config.phase3_n_neighbors}')
print(f'Phase 4 permutations: {config.phase4_n_permutations}')
print(f'Output: {config.output_dir}')
print()
print('PASS: operational shell working correctly')
"
```

Regression:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_stats.py tests/test_ingestion.py -q 2>&1
```

Report all three outputs back.
