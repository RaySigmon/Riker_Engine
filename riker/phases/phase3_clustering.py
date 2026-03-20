"""
Riker Engine - Phase 3: Consensus Clustering.

Applies UMAP dimensionality reduction followed by HDBSCAN density clustering
across multiple configurations, then derives consensus clusters from the
co-association matrix.

No biological hypotheses are imposed during clustering. This is agnostic
discovery (Blueprint Section 3). Post-hoc pathway enrichment labels are
descriptive only and must NOT be reported as independent evidence.

References:
    Blueprint Section 7 (Phase 3: Consensus Clustering)
    Blueprint Section 8.1 (Bonferroni correction)
"""

import logging
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Import UMAP
try:
    from umap import UMAP
    _HAS_UMAP = True
except ImportError:
    _HAS_UMAP = False

# Import HDBSCAN: prefer sklearn (core dep), fallback to standalone
try:
    from sklearn.cluster import HDBSCAN as _HDBSCAN
    _HDBSCAN_SOURCE = "sklearn"
except ImportError:
    try:
        from hdbscan import HDBSCAN as _HDBSCAN
        _HDBSCAN_SOURCE = "hdbscan"
    except ImportError:
        _HDBSCAN = None
        _HDBSCAN_SOURCE = None

# Default configuration from Blueprint Section 7
DEFAULT_N_NEIGHBORS = [10, 15, 30]
DEFAULT_SEEDS = [42, 123, 456, 789, 1024]
DEFAULT_MIN_CLUSTER_SIZE = 5
DEFAULT_MIN_SAMPLES = 3


@dataclass(frozen=True)
class ClusterInfo:
    """Information about a single consensus cluster.

    Attributes:
        cluster_id: Integer cluster label (0-indexed).
        gene_symbols: List of gene symbols in this cluster.
        n_genes: Number of genes.
        mean_consensus: Mean pairwise consensus score within cluster.
    """
    cluster_id: int
    gene_symbols: list
    n_genes: int
    mean_consensus: float


@dataclass
class Phase3Result:
    """Complete result of Phase 3 consensus clustering.

    Attributes:
        cluster_labels: Dict of gene_symbol -> cluster_id (-1 = noise).
        consensus_matrix: DataFrame (genes × genes) with co-clustering freq.
        n_clusters: Number of clusters (excluding noise).
        n_noise: Number of unclustered genes (label = -1).
        cluster_info: Dict of cluster_id -> ClusterInfo.
        gene_order: List of gene symbols in matrix order.
        per_config_labels: List of dicts, one per configuration run.
        n_configurations: Number of UMAP/HDBSCAN configs used.
        umap_params: Dict of UMAP parameters used.
        hdbscan_params: Dict of HDBSCAN parameters used.
    """
    cluster_labels: dict = field(default_factory=dict)
    consensus_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    n_clusters: int = 0
    n_noise: int = 0
    cluster_info: dict = field(default_factory=dict)
    gene_order: list = field(default_factory=list)
    per_config_labels: list = field(default_factory=list)
    n_configurations: int = 0
    umap_params: dict = field(default_factory=dict)
    hdbscan_params: dict = field(default_factory=dict)


def _check_dependencies() -> None:
    """Verify clustering dependencies are available."""
    if not _HAS_UMAP:
        raise ImportError(
            "UMAP not available. Install with: pip install umap-learn"
        )
    if _HDBSCAN is None:
        raise ImportError(
            "HDBSCAN not available. Install with: "
            "pip install hdbscan or pip install scikit-learn>=1.3"
        )


def _run_single_config(
    feature_matrix: np.ndarray,
    n_neighbors: int,
    random_state: int,
    min_cluster_size: int,
    min_samples: int,
) -> np.ndarray:
    """Run one UMAP + HDBSCAN configuration.

    Parameters
    ----------
    feature_matrix : np.ndarray
        Genes × features matrix (already normalized).
    n_neighbors : int
        UMAP n_neighbors parameter.
    random_state : int
        Random seed for UMAP.
    min_cluster_size : int
        HDBSCAN min_cluster_size.
    min_samples : int
        HDBSCAN min_samples.

    Returns
    -------
    np.ndarray
        Cluster labels per gene (-1 = noise).
    """
    # UMAP reduction to 2D
    reducer = UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=0.0,
        metric="euclidean",
        random_state=random_state,
    )

    embedding = reducer.fit_transform(feature_matrix)

    # HDBSCAN clustering on the 2D embedding
    if _HDBSCAN_SOURCE == "sklearn":
        clusterer = _HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
        )
        labels = clusterer.fit_predict(embedding)
    else:
        clusterer = _HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
        )
        clusterer.fit(embedding)
        labels = clusterer.labels_

    return labels


def build_consensus_matrix(
    all_labels: list[np.ndarray],
    n_genes: int,
) -> np.ndarray:
    """Build the co-association consensus matrix.

    For each pair of genes (i, j), the consensus score is the fraction
    of configurations in which both genes were assigned to the SAME
    cluster (excluding configurations where either gene was noise = -1).

    Parameters
    ----------
    all_labels : list of np.ndarray
        Cluster labels from each configuration.
    n_genes : int
        Number of genes.

    Returns
    -------
    np.ndarray
        n × n symmetric consensus matrix with values in [0, 1].
    """
    consensus = np.zeros((n_genes, n_genes), dtype=np.float64)
    counts = np.zeros((n_genes, n_genes), dtype=np.float64)

    for labels in all_labels:
        for i in range(n_genes):
            if labels[i] == -1:
                continue
            for j in range(i, n_genes):
                if labels[j] == -1:
                    continue
                counts[i, j] += 1
                counts[j, i] += 1
                if labels[i] == labels[j]:
                    consensus[i, j] += 1
                    consensus[j, i] += 1

    # Normalize: consensus / counts (where counts > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        consensus = np.divide(consensus, counts,
                              out=np.zeros_like(consensus),
                              where=counts > 0)

    # Diagonal = 1 (gene always co-clusters with itself)
    np.fill_diagonal(consensus, 1.0)

    return consensus


def _cluster_from_consensus(
    consensus: np.ndarray,
    min_cluster_size: int,
    min_samples: int,
) -> np.ndarray:
    """Derive final clusters from the consensus matrix.

    Converts the consensus matrix to a distance matrix (1 - consensus)
    and applies HDBSCAN with metric='precomputed'.

    Parameters
    ----------
    consensus : np.ndarray
        n × n consensus matrix.
    min_cluster_size : int
        HDBSCAN min_cluster_size for final clustering.
    min_samples : int
        HDBSCAN min_samples for final clustering.

    Returns
    -------
    np.ndarray
        Final cluster labels.
    """
    distance = 1.0 - consensus
    # Ensure diagonal is 0 and matrix is symmetric
    np.fill_diagonal(distance, 0.0)
    distance = (distance + distance.T) / 2.0
    # Clip to valid range
    distance = np.clip(distance, 0.0, 1.0)

    if _HDBSCAN_SOURCE == "sklearn":
        clusterer = _HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            metric="precomputed",
        )
        labels = clusterer.fit_predict(distance)
    else:
        clusterer = _HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            metric="precomputed",
        )
        clusterer.fit(distance)
        labels = clusterer.labels_

    return labels


def run_consensus_clustering(
    feature_matrix: pd.DataFrame,
    n_neighbors_list: list[int] | None = None,
    seeds: list[int] | None = None,
    min_cluster_size: int = DEFAULT_MIN_CLUSTER_SIZE,
    min_samples: int = DEFAULT_MIN_SAMPLES,
    final_min_cluster_size: int | None = None,
    final_min_samples: int | None = None,
) -> Phase3Result:
    """Run Phase 3 consensus clustering.

    Parameters
    ----------
    feature_matrix : DataFrame
        Genes × features from Phase 2 (min-max normalized, index = gene symbols).
    n_neighbors_list : list of int
        UMAP n_neighbors values. Default: [10, 15, 30].
    seeds : list of int
        Random seeds. Default: [42, 123, 456, 789, 1024].
    min_cluster_size : int
        HDBSCAN min_cluster_size for per-config runs. Default: 5.
    min_samples : int
        HDBSCAN min_samples for per-config runs. Default: 3.
    final_min_cluster_size : int or None
        Override min_cluster_size for final consensus clustering.
        If None, uses same as per-config.
    final_min_samples : int or None
        Override min_samples for final consensus clustering.

    Returns
    -------
    Phase3Result
    """
    _check_dependencies()

    if n_neighbors_list is None:
        n_neighbors_list = DEFAULT_N_NEIGHBORS
    if seeds is None:
        seeds = DEFAULT_SEEDS

    gene_order = list(feature_matrix.index)
    n_genes = len(gene_order)
    X = feature_matrix.values.astype(np.float64)

    if n_genes < min_cluster_size:
        warnings.warn(
            f"Only {n_genes} genes — fewer than min_cluster_size "
            f"({min_cluster_size}). Clustering may produce all noise.",
            UserWarning,
            stacklevel=2,
        )

    # Run all configurations
    all_labels = []
    per_config_labels = []
    n_configs = len(n_neighbors_list) * len(seeds)

    logger.info(
        f"Running {n_configs} UMAP+HDBSCAN configurations "
        f"({len(n_neighbors_list)} n_neighbors × {len(seeds)} seeds)..."
    )

    for nn in n_neighbors_list:
        # Cap n_neighbors at n_genes - 1
        effective_nn = min(nn, n_genes - 1)
        if effective_nn < 2:
            effective_nn = 2

        for seed in seeds:
            try:
                labels = _run_single_config(
                    X, effective_nn, seed, min_cluster_size, min_samples
                )
            except Exception as e:
                logger.warning(
                    f"Config n_neighbors={nn}, seed={seed} failed: {e}. Skipping."
                )
                continue

            all_labels.append(labels)
            config_dict = {
                gene: int(labels[i]) for i, gene in enumerate(gene_order)
            }
            per_config_labels.append({
                "n_neighbors": nn,
                "seed": seed,
                "labels": config_dict,
                "n_clusters": len(set(labels)) - (1 if -1 in labels else 0),
            })

    if not all_labels:
        warnings.warn(
            "All clustering configurations failed. Returning empty result.",
            UserWarning,
            stacklevel=2,
        )
        return Phase3Result(gene_order=gene_order, n_configurations=0)

    # Build consensus matrix
    consensus = build_consensus_matrix(all_labels, n_genes)

    # Final clustering from consensus
    final_mcs = final_min_cluster_size or min_cluster_size
    final_ms = final_min_samples or min_samples
    final_labels = _cluster_from_consensus(consensus, final_mcs, final_ms)

    # Build result
    cluster_labels = {
        gene: int(final_labels[i]) for i, gene in enumerate(gene_order)
    }

    # Compute cluster info
    unique_labels = set(final_labels)
    unique_labels.discard(-1)
    cluster_info = {}

    for cid in sorted(unique_labels):
        mask = final_labels == cid
        genes_in_cluster = [gene_order[i] for i in range(n_genes) if mask[i]]

        # Mean within-cluster consensus
        indices = np.where(mask)[0]
        if len(indices) > 1:
            pairs = []
            for ii in range(len(indices)):
                for jj in range(ii + 1, len(indices)):
                    pairs.append(consensus[indices[ii], indices[jj]])
            mean_cons = float(np.mean(pairs)) if pairs else 1.0
        else:
            mean_cons = 1.0

        cluster_info[cid] = ClusterInfo(
            cluster_id=cid,
            gene_symbols=genes_in_cluster,
            n_genes=len(genes_in_cluster),
            mean_consensus=mean_cons,
        )

    n_clusters = len(unique_labels)
    n_noise = sum(1 for l in final_labels if l == -1)

    consensus_df = pd.DataFrame(
        consensus, index=gene_order, columns=gene_order
    )

    logger.info(
        f"Phase 3 complete: {n_clusters} clusters, {n_noise} noise genes "
        f"from {len(all_labels)} configurations."
    )

    return Phase3Result(
        cluster_labels=cluster_labels,
        consensus_matrix=consensus_df,
        n_clusters=n_clusters,
        n_noise=n_noise,
        cluster_info=cluster_info,
        gene_order=gene_order,
        per_config_labels=per_config_labels,
        n_configurations=len(all_labels),
        umap_params={
            "n_components": 2,
            "min_dist": 0.0,
            "metric": "euclidean",
            "n_neighbors_list": n_neighbors_list,
            "seeds": seeds,
        },
        hdbscan_params={
            "min_cluster_size": min_cluster_size,
            "min_samples": min_samples,
            "final_min_cluster_size": final_mcs,
            "final_min_samples": final_ms,
            "source": _HDBSCAN_SOURCE,
        },
    )
