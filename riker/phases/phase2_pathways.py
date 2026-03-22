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
Riker Engine - Phase 2: Pathway Mapping and Feature Construction.

Maps study genes to biological pathways (KEGG, Reactome, MSigDB Hallmark)
and constructs a feature matrix for unsupervised clustering in Phase 3.

Anti-circularity rule (Blueprint Section 7): Individual pathway IDs are
used as binary features. Pre-assigned biological category labels must
NEVER be used as features — this would create tautological enrichment.

References:
    Blueprint Section 6 (Phase 2: Pathway Mapping)
    Blueprint Section 7 (Feature vector construction)
    Blueprint Section 12 (QC Framework — circularity audit)
"""

import logging
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PathwayInfo:
    """Metadata for a single pathway.

    Attributes:
        pathway_id: Unique identifier (e.g., 'hsa04020', 'R-HSA-123456').
        name: Human-readable pathway name.
        source: 'KEGG', 'Reactome', or 'Hallmark'.
        total_genes: Total number of genes in this pathway.
        study_genes: Number of study genes mapped to this pathway.
        study_gene_list: List of study gene symbols in this pathway.
    """
    pathway_id: str
    name: str
    source: str
    total_genes: int
    study_genes: int
    study_gene_list: list


@dataclass
class Phase2Result:
    """Complete result of Phase 2 pathway mapping.

    Attributes:
        feature_matrix: DataFrame (genes × features), min-max normalized.
        pathway_features: List of pathway IDs used as features.
        pathway_info: Dict of pathway_id -> PathwayInfo.
        n_pathways_before_filter: Total pathways loaded.
        n_pathways_after_filter: Pathways surviving filter.
        n_genes: Number of study genes in the feature matrix.
        n_features: Total number of features (pathways + expression stats).
        expression_feature_names: Names of non-pathway features.
    """
    feature_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    pathway_features: list = field(default_factory=list)
    pathway_info: dict = field(default_factory=dict)
    n_pathways_before_filter: int = 0
    n_pathways_after_filter: int = 0
    n_genes: int = 0
    n_features: int = 0
    expression_feature_names: list = field(default_factory=list)


class PathwayDatabase:
    """Holds pathway-to-gene mappings from one or more sources.

    Parameters
    ----------
    pathways : dict
        Mapping of pathway_id -> set of gene symbols.
    names : dict or None
        Mapping of pathway_id -> human-readable name.
    source : str
        Source label ('KEGG', 'Reactome', 'Hallmark', or 'mixed').
    """

    def __init__(
        self,
        pathways: dict[str, set[str]] | None = None,
        names: dict[str, str] | None = None,
        source: str = "unknown",
    ):
        self.pathways: dict[str, set[str]] = pathways or {}
        self.names: dict[str, str] = names or {}
        self.source = source

    def add_pathway(self, pathway_id: str, genes: set[str],
                    name: str = "") -> None:
        """Add a single pathway."""
        self.pathways[pathway_id] = genes
        if name:
            self.names[pathway_id] = name

    def merge(self, other: "PathwayDatabase") -> None:
        """Merge another PathwayDatabase into this one."""
        self.pathways.update(other.pathways)
        self.names.update(other.names)
        if self.source != other.source:
            self.source = "mixed"

    def get_gene_pathways(self, gene: str) -> list[str]:
        """Get all pathway IDs that contain a given gene."""
        return [pid for pid, genes in self.pathways.items() if gene in genes]

    @property
    def n_pathways(self) -> int:
        return len(self.pathways)

    @property
    def all_genes(self) -> set[str]:
        """All unique genes across all pathways."""
        result = set()
        for genes in self.pathways.values():
            result.update(genes)
        return result


def load_pathways_from_dict(
    data: dict[str, list[str] | set[str]],
    names: dict[str, str] | None = None,
    source: str = "custom",
) -> PathwayDatabase:
    """Build a PathwayDatabase from a plain dict.

    Parameters
    ----------
    data : dict
        Mapping of pathway_id -> list or set of gene symbols.
    names : dict or None
        Mapping of pathway_id -> name.
    source : str
        Source label.

    Returns
    -------
    PathwayDatabase
    """
    pathways = {pid: set(genes) for pid, genes in data.items()}
    return PathwayDatabase(pathways=pathways, names=names or {}, source=source)


def filter_pathways(
    db: PathwayDatabase,
    study_genes: list[str],
    min_study_genes: int = 3,
    max_total_genes: int = 500,
    min_study_fraction: float = 0.02,
    max_pathways: int = 100,
) -> tuple[PathwayDatabase, dict[str, PathwayInfo]]:
    """Filter pathways by study gene membership criteria.

    Parameters
    ----------
    db : PathwayDatabase
        Unfiltered pathway database.
    study_genes : list of str
        Genes from Phase 1 study gene set.
    min_study_genes : int
        Minimum study gene members per pathway (default 3).
    max_total_genes : int
        Maximum total genes per pathway (default 500).
    min_study_fraction : float
        Minimum fraction of pathway genes that are study genes (default 0.02).
    max_pathways : int
        Maximum pathways to retain after filtering (top by study gene count).

    Returns
    -------
    tuple of (PathwayDatabase, dict)
        Filtered database and pathway info dict.
    """
    study_set = set(study_genes)
    filtered = PathwayDatabase(source=db.source)
    info = {}

    # Score each pathway
    candidates = []
    for pid, genes in db.pathways.items():
        total = len(genes)
        if total > max_total_genes:
            continue
        if total == 0:
            continue

        study_overlap = genes & study_set
        n_study = len(study_overlap)

        if n_study < min_study_genes:
            continue

        fraction = n_study / total
        if fraction < min_study_fraction:
            continue

        candidates.append((pid, genes, n_study, total, study_overlap))

    # Sort by study gene count descending, take top max_pathways
    candidates.sort(key=lambda x: x[2], reverse=True)
    candidates = candidates[:max_pathways]

    for pid, genes, n_study, total, overlap in candidates:
        filtered.add_pathway(pid, genes, db.names.get(pid, ""))
        info[pid] = PathwayInfo(
            pathway_id=pid,
            name=db.names.get(pid, ""),
            source=db.source,
            total_genes=total,
            study_genes=n_study,
            study_gene_list=sorted(overlap),
        )

    logger.info(
        f"Pathway filter: {db.n_pathways} -> {filtered.n_pathways} pathways "
        f"(min_study={min_study_genes}, max_total={max_total_genes}, "
        f"min_frac={min_study_fraction})."
    )

    return filtered, info


def build_feature_matrix(
    study_genes: dict,
    filtered_db: PathwayDatabase,
    pathway_info: dict[str, PathwayInfo],
    gene_tiers: dict[str, float] | None = None,
) -> Phase2Result:
    """Build the feature matrix for clustering.

    Parameters
    ----------
    study_genes : dict
        Mapping of gene_symbol -> GeneResult from Phase 1.
        Each GeneResult has de_results with per-dataset stats.
    filtered_db : PathwayDatabase
        Filtered pathway database from filter_pathways().
    pathway_info : dict
        PathwayInfo dict from filter_pathways().
    gene_tiers : dict or None
        Mapping of gene_symbol -> confidence tier score (numeric).
        Higher = more confident. If None, tier feature is set to 0.

    Returns
    -------
    Phase2Result
        Contains the feature matrix and metadata.
    """
    gene_symbols = sorted(study_genes.keys())
    pathway_ids = sorted(filtered_db.pathways.keys())

    if not gene_symbols:
        return Phase2Result()

    # --- Binary pathway features ---
    pathway_data = {}
    for pid in pathway_ids:
        members = filtered_db.pathways[pid]
        pathway_data[pid] = [1.0 if g in members else 0.0 for g in gene_symbols]

    # --- Expression statistic features ---
    avg_log2fc = []
    neg_log10_min_p = []
    n_sig_datasets = []
    direction_score = []
    tier_score = []

    for gene in gene_symbols:
        gene_result = study_genes[gene]
        de_list = gene_result.de_results

        # Average log2FC
        if de_list:
            avg_fc = float(np.mean([d.log2fc for d in de_list]))
        else:
            avg_fc = 0.0
        avg_log2fc.append(avg_fc)

        # -log10(min p-value), capped at 10
        if de_list:
            min_p = min(d.p_value for d in de_list)
            neg_log10 = min(-np.log10(max(min_p, 1e-10)), 10.0)
        else:
            neg_log10 = 0.0
        neg_log10_min_p.append(neg_log10)

        # Number of significant datasets
        n_sig = gene_result.n_datasets_significant
        n_sig_datasets.append(float(n_sig))

        # Directional consistency: fraction of datasets with majority direction
        if de_list:
            directions = [d.direction for d in de_list]
            n_up = directions.count("up")
            n_down = directions.count("down")
            consistency = max(n_up, n_down) / len(directions)
        else:
            consistency = 0.0
        direction_score.append(consistency)

        # Tier score (inverted: tier 1 = highest confidence = highest score)
        if gene_tiers and gene in gene_tiers:
            tier_score.append(gene_tiers[gene])
        else:
            tier_score.append(0.0)

    # --- Assemble feature matrix ---
    expr_features = {
        "avg_log2fc": avg_log2fc,
        "neg_log10_min_p": neg_log10_min_p,
        "n_sig_datasets": n_sig_datasets,
        "direction_consistency": direction_score,
        "tier_score": tier_score,
    }

    all_features = {}
    all_features.update(pathway_data)
    all_features.update(expr_features)

    feature_df = pd.DataFrame(all_features, index=gene_symbols)
    feature_df.index.name = "gene"

    # --- Min-max normalize to [0, 1] ---
    for col in feature_df.columns:
        col_min = feature_df[col].min()
        col_max = feature_df[col].max()
        if col_max > col_min:
            feature_df[col] = (feature_df[col] - col_min) / (col_max - col_min)
        else:
            # Constant column — set to 0 (or 1 if all nonzero)
            feature_df[col] = 0.0 if col_min == 0 else 1.0

    expr_feature_names = list(expr_features.keys())

    return Phase2Result(
        feature_matrix=feature_df,
        pathway_features=pathway_ids,
        pathway_info=pathway_info,
        n_pathways_before_filter=0,  # set by caller if needed
        n_pathways_after_filter=len(pathway_ids),
        n_genes=len(gene_symbols),
        n_features=len(feature_df.columns),
        expression_feature_names=expr_feature_names,
    )


def run_phase2(
    phase1_result,
    pathway_db: PathwayDatabase,
    gene_tiers: dict[str, float] | None = None,
    min_study_genes: int = 3,
    max_total_genes: int = 500,
    min_study_fraction: float = 0.02,
    max_pathways: int = 100,
) -> Phase2Result:
    """Run Phase 2: filter pathways and build feature matrix.

    Parameters
    ----------
    phase1_result : Phase1Result
        Output from Phase 1 cross-referencing.
    pathway_db : PathwayDatabase
        Combined pathway database (KEGG + Reactome + Hallmark).
    gene_tiers : dict or None
        Confidence tier scores per gene.
    min_study_genes, max_total_genes, min_study_fraction, max_pathways :
        Pathway filter parameters (see filter_pathways).

    Returns
    -------
    Phase2Result
    """
    study_gene_list = list(phase1_result.study_genes.keys())

    n_before = pathway_db.n_pathways
    filtered_db, pathway_info = filter_pathways(
        pathway_db, study_gene_list,
        min_study_genes=min_study_genes,
        max_total_genes=max_total_genes,
        min_study_fraction=min_study_fraction,
        max_pathways=max_pathways,
    )

    result = build_feature_matrix(
        phase1_result.study_genes, filtered_db, pathway_info,
        gene_tiers=gene_tiers,
    )
    result.n_pathways_before_filter = n_before

    logger.info(
        f"Phase 2 complete: {result.n_genes} genes × {result.n_features} features "
        f"({result.n_pathways_after_filter} pathway + "
        f"{len(result.expression_feature_names)} expression features)."
    )

    return result
