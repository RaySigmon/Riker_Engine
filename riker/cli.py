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

        from riker.io.outputs import write_phase4_core_genes, write_phase4_all_levels
        write_phase4_core_genes(phase4_result, output_dir)
        write_phase4_all_levels(phase1_result, phase3_result, phase4_result, output_dir)

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


def cmd_ui(args) -> int:
    """Launch the web UI."""
    try:
        import uvicorn
        from riker.ui.server import app
    except ImportError:
        print(
            "UI dependencies not installed.\n"
            "Run: pip install riker-engine[ui]\n"
            "  or: pip install fastapi uvicorn jinja2 python-multipart",
            file=sys.stderr,
        )
        return 1

    host = args.host
    port = args.port

    if not args.no_browser:
        import webbrowser
        import threading
        threading.Timer(1.5, lambda: webbrowser.open(f"http://{host}:{port}")).start()

    print(f"Starting Riker Engine UI at http://{host}:{port}")
    uvicorn.run(app, host=host, port=port, log_level="warning")
    return 0


def main() -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="riker",
        description="Riker Engine — Condition-agnostic transcriptomics pipeline. "
                    "Run with no arguments to launch the web UI.",
    )
    parser.add_argument("--no-browser", action="store_true",
                        help="(UI mode) Do not auto-open browser")
    parser.add_argument("--host", default="127.0.0.1",
                        help="(UI mode) Host to bind to")
    parser.add_argument("--port", type=int, default=8000,
                        help="(UI mode) Port to bind to")
    subparsers = parser.add_subparsers(dest="command")

    # run command
    run_parser = subparsers.add_parser("run", help="Run the full pipeline")
    run_parser.add_argument("config", help="Path to YAML config file")
    run_parser.add_argument("-v", "--verbose", action="store_true",
                            help="Enable debug logging")

    # validate command
    val_parser = subparsers.add_parser("validate", help="Validate config file")
    val_parser.add_argument("config", help="Path to YAML config file")

    # ui command
    ui_parser = subparsers.add_parser("ui", help="Launch the web interface")
    ui_parser.add_argument("--host", default="127.0.0.1", help="Host to bind to")
    ui_parser.add_argument("--port", type=int, default=8000, help="Port to bind to")
    ui_parser.add_argument("--no-browser", action="store_true",
                           help="Do not auto-open browser")

    args = parser.parse_args()

    # No subcommand = launch the UI (the default experience)
    if args.command is None:
        args.host = "127.0.0.1"
        args.port = 8000
        args.no_browser = False
        return cmd_ui(args)

    setup_logging(verbose=getattr(args, "verbose", False))

    if args.command == "validate":
        return cmd_validate(args)
    elif args.command == "run":
        return cmd_run(args)
    elif args.command == "ui":
        return cmd_ui(args)

    return 1


if __name__ == "__main__":
    sys.exit(main())
