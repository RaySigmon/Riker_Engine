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
    if not config.seed_genes_path:
        raise ValueError(
            "seed_genes must be a path to a seed gene CSV file (got null/empty). "
            "Blind discovery mode (no seed genes) is not yet supported."
        )

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
