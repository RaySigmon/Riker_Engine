# INSTRUCTION SET FOR KAI — PHASE 10: `riker/phases/phase3_clustering.py`

## References
- Blueprint Section 7 (Phase 3: Consensus Clustering)
- Blueprint Section 8.1 (Bonferroni correction for multiple clusters)
- Blueprint Section 12 (QC Framework — permutation test row)

## WHY THIS MODULE MATTERS

Phase 3 is where the engine discovers gene modules. UMAP reduces the
feature matrix from Phase 2 into 2D, then HDBSCAN identifies density
clusters. Running this across 15 configurations (3 n_neighbors × 5 seeds)
and building a consensus matrix makes the results robust to parameter
choices. This is the "agnostic discovery" principle from Blueprint
Section 3 — no biological hypotheses are imposed during clustering.

## BLUEPRINT SECTION 7 SPECIFICATION

- UMAP: n_components=2, min_dist=0.0, metric='euclidean'
- 3 n_neighbors values: [10, 15, 30]
- 5 random seeds: [42, 123, 456, 789, 1024]
- = 15 total configurations
- HDBSCAN: min_cluster_size=5, min_samples=3
- Co-association consensus matrix: pairwise co-clustering frequencies
- Final clusters: HDBSCAN(metric='precomputed') on distance matrix (1 - consensus)
- Post-hoc: Fisher's exact test for pathway enrichment (labels ONLY, not evidence)

## CRITICAL REQUIREMENTS

1. `run_consensus_clustering()`: Main function that:
   - Takes the feature matrix from Phase 2
   - Runs UMAP + HDBSCAN across all 15 configurations
   - Builds the co-association consensus matrix
   - Derives final cluster labels from the consensus
   - Returns cluster assignments + consensus matrix

2. Use `sklearn.cluster.HDBSCAN` as primary (sklearn 1.8.0 confirmed on Ghost).
   Fall back to `hdbscan.HDBSCAN` if sklearn version doesn't have it.

3. The consensus matrix is an n×n symmetric matrix where entry (i,j) is the
   fraction of configurations in which genes i and j were placed in the same
   cluster (excluding noise assignments where label = -1).

4. Final clustering applies HDBSCAN with metric='precomputed' to the
   DISTANCE matrix (1 - consensus_matrix).

5. Parameters must be configurable (n_neighbors list, seeds, HDBSCAN params).

6. Return a `Phase3Result` dataclass with cluster assignments, consensus
   matrix, per-configuration labels, and cluster statistics.

7. Since UMAP can be slow on Pi 5, keep the test data SMALL (< 50 genes)
   and use fewer configurations for tests (e.g., 2 neighbors × 2 seeds = 4).

8. DO NOT modify any existing files. APPEND tests to test_phases.py.

---

## FILE: `riker/phases/phase3_clustering.py`

Write the following file at `/home/kai001/riker-engine/riker/phases/phase3_clustering.py`:

```python
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
```

---

## TESTS: APPEND to `tests/test_phases.py`

**APPEND** the following to the END of `/home/kai001/riker-engine/tests/test_phases.py`.
Do NOT delete or modify any existing test classes.

These tests use SMALL data and FEW configurations to keep runtime reasonable on Pi 5.

```python


# ===========================================================================
# PHASE 10: CONSENSUS CLUSTERING TESTS
# ===========================================================================

from riker.phases.phase3_clustering import (
    ClusterInfo,
    Phase3Result,
    build_consensus_matrix,
    run_consensus_clustering,
)


def _make_clusterable_features(n_genes=30, n_features=10, n_clusters=3, seed=42):
    """Create a feature matrix with clear cluster structure.

    Generates n_clusters groups of genes with distinct feature patterns,
    making them easily separable by UMAP + HDBSCAN.
    """
    np.random.seed(seed)
    genes_per_cluster = n_genes // n_clusters
    gene_names = []
    data = []

    for c in range(n_clusters):
        for g in range(genes_per_cluster):
            gene_name = f"C{c}_G{g}"
            gene_names.append(gene_name)
            # Base features: cluster-specific center + noise
            center = np.zeros(n_features)
            center[c * (n_features // n_clusters):(c + 1) * (n_features // n_clusters)] = 1.0
            features = center + np.random.normal(0, 0.1, n_features)
            features = np.clip(features, 0, 1)
            data.append(features)

    df = pd.DataFrame(data, index=gene_names,
                      columns=[f"F{i}" for i in range(n_features)])
    df.index.name = "gene"
    return df


# ---------------------------------------------------------------------------
# 7. Consensus matrix construction
# ---------------------------------------------------------------------------

class TestBuildConsensusMatrix:
    """Test the co-association consensus matrix."""

    def test_perfect_agreement(self):
        """If all configs agree, consensus should be 1 within clusters."""
        # 6 genes, 2 clusters of 3 each, all configs agree
        labels1 = np.array([0, 0, 0, 1, 1, 1])
        labels2 = np.array([0, 0, 0, 1, 1, 1])
        consensus = build_consensus_matrix([labels1, labels2], 6)

        # Within cluster 0: should be 1.0
        assert consensus[0, 1] == 1.0
        assert consensus[0, 2] == 1.0
        # Between clusters: should be 0.0
        assert consensus[0, 3] == 0.0
        # Diagonal should be 1.0
        assert consensus[0, 0] == 1.0

    def test_partial_agreement(self):
        """Mixed agreement should give intermediate values."""
        labels1 = np.array([0, 0, 1, 1])
        labels2 = np.array([0, 1, 0, 1])  # different grouping
        consensus = build_consensus_matrix([labels1, labels2], 4)

        # Gene 0 and Gene 1: together in config 1, apart in config 2 = 0.5
        assert consensus[0, 1] == 0.5
        # Gene 0 and Gene 3: apart in both = 0.0
        assert consensus[0, 3] == 0.0

    def test_noise_excluded(self):
        """Noise labels (-1) should be excluded from consensus counting."""
        labels1 = np.array([0, 0, -1, 1])
        labels2 = np.array([0, 0, 0, 1])
        consensus = build_consensus_matrix([labels1, labels2], 4)

        # Gene 2 was noise in config 1: only 1 valid config for pairs with gene 2
        # Gene 0 and Gene 2: same cluster in config 2 only (config 1 excluded)
        assert consensus[0, 2] == 1.0  # 1/1 valid configs

    def test_symmetry(self):
        """Consensus matrix should be symmetric."""
        np.random.seed(42)
        labels = [np.random.randint(-1, 3, size=10) for _ in range(5)]
        consensus = build_consensus_matrix(labels, 10)
        np.testing.assert_array_almost_equal(consensus, consensus.T)

    def test_range(self):
        """All values should be in [0, 1]."""
        np.random.seed(42)
        labels = [np.random.randint(-1, 3, size=10) for _ in range(5)]
        consensus = build_consensus_matrix(labels, 10)
        assert consensus.min() >= 0.0
        assert consensus.max() <= 1.0


# ---------------------------------------------------------------------------
# 8. Full consensus clustering pipeline
# ---------------------------------------------------------------------------

class TestRunConsensusClustering:
    """Test the full Phase 3 pipeline.

    Uses small data and few configs to keep Pi 5 runtime reasonable.
    """

    def test_finds_clusters(self):
        """Should find clusters in clearly separable data."""
        features = _make_clusterable_features(n_genes=30, n_clusters=3)
        result = run_consensus_clustering(
            features,
            n_neighbors_list=[5, 10],
            seeds=[42, 99],
            min_cluster_size=3,
            min_samples=2,
        )

        # Should find at least 2 clusters
        assert result.n_clusters >= 2
        assert result.n_configurations == 4  # 2 neighbors × 2 seeds

    def test_consensus_matrix_shape(self):
        """Consensus matrix should be n_genes × n_genes."""
        features = _make_clusterable_features(n_genes=20, n_clusters=2)
        result = run_consensus_clustering(
            features,
            n_neighbors_list=[5],
            seeds=[42],
            min_cluster_size=3,
            min_samples=2,
        )

        assert result.consensus_matrix.shape == (20, 20)

    def test_gene_order_preserved(self):
        """Gene order should match feature matrix index."""
        features = _make_clusterable_features(n_genes=15, n_clusters=3)
        result = run_consensus_clustering(
            features,
            n_neighbors_list=[5],
            seeds=[42],
            min_cluster_size=3,
            min_samples=2,
        )

        assert result.gene_order == list(features.index)

    def test_cluster_info_populated(self):
        """ClusterInfo should be created for each cluster."""
        features = _make_clusterable_features(n_genes=30, n_clusters=3)
        result = run_consensus_clustering(
            features,
            n_neighbors_list=[5, 10],
            seeds=[42, 99],
            min_cluster_size=3,
            min_samples=2,
        )

        for cid, info in result.cluster_info.items():
            assert info.n_genes == len(info.gene_symbols)
            assert info.n_genes > 0
            assert 0 <= info.mean_consensus <= 1.0

    def test_all_genes_assigned(self):
        """Every gene should have a label (even if -1 for noise)."""
        features = _make_clusterable_features(n_genes=20, n_clusters=2)
        result = run_consensus_clustering(
            features,
            n_neighbors_list=[5],
            seeds=[42],
            min_cluster_size=3,
            min_samples=2,
        )

        assert len(result.cluster_labels) == 20
        for gene in features.index:
            assert gene in result.cluster_labels

    def test_per_config_labels_stored(self):
        """Should store labels from each configuration."""
        features = _make_clusterable_features(n_genes=20, n_clusters=2)
        result = run_consensus_clustering(
            features,
            n_neighbors_list=[5, 10],
            seeds=[42, 99],
            min_cluster_size=3,
            min_samples=2,
        )

        assert len(result.per_config_labels) == result.n_configurations
        for cfg in result.per_config_labels:
            assert "n_neighbors" in cfg
            assert "seed" in cfg
            assert "labels" in cfg
            assert "n_clusters" in cfg

    def test_params_recorded(self):
        """UMAP and HDBSCAN params should be recorded in result."""
        features = _make_clusterable_features(n_genes=15, n_clusters=3)
        result = run_consensus_clustering(
            features,
            n_neighbors_list=[5],
            seeds=[42],
            min_cluster_size=4,
            min_samples=2,
        )

        assert result.umap_params["n_components"] == 2
        assert result.umap_params["min_dist"] == 0.0
        assert result.hdbscan_params["min_cluster_size"] == 4

    def test_noise_count(self):
        """n_clusters + n_noise labels should cover all genes."""
        features = _make_clusterable_features(n_genes=30, n_clusters=3)
        result = run_consensus_clustering(
            features,
            n_neighbors_list=[5, 10],
            seeds=[42, 99],
            min_cluster_size=3,
            min_samples=2,
        )

        clustered = sum(info.n_genes for info in result.cluster_info.values())
        assert clustered + result.n_noise == len(features)
```

---

## EXECUTION INSTRUCTIONS

After writing phase3_clustering.py and appending tests, run:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_phases.py -v 2>&1
```

**Expected: ALL tests pass (Phase 8 + 9 + 10).** If any fail, report FULL output — do NOT modify tests or code without reporting first.

Note: UMAP tests will be slower than previous tests due to dimensionality reduction. Expect 10-30 seconds total on Pi 5.

Then confirmation:

```bash
cd /home/kai001/riker-engine && python -c "
import numpy as np
import pandas as pd
from riker.phases.phase3_clustering import run_consensus_clustering

# Create 40 genes with 4 clear clusters
np.random.seed(42)
n_per = 10
n_features = 12
genes = []
data = []

for c in range(4):
    for g in range(n_per):
        genes.append(f'CLUSTER{c}_GENE{g}')
        center = np.zeros(n_features)
        center[c*3:(c+1)*3] = 1.0
        feat = center + np.random.normal(0, 0.08, n_features)
        feat = np.clip(feat, 0, 1)
        data.append(feat)

features = pd.DataFrame(data, index=genes,
                         columns=[f'F{i}' for i in range(n_features)])
features.index.name = 'gene'

# Run with reduced configs for speed
result = run_consensus_clustering(
    features,
    n_neighbors_list=[5, 10],
    seeds=[42, 99, 123],
    min_cluster_size=3,
    min_samples=2,
)

print(f'=== Phase 3 Consensus Clustering ===')
print(f'Configurations: {result.n_configurations}')
print(f'Clusters found: {result.n_clusters}')
print(f'Noise genes: {result.n_noise}')
print(f'Consensus matrix: {result.consensus_matrix.shape}')
print()

for cid, info in sorted(result.cluster_info.items()):
    print(f'  Cluster {cid}: {info.n_genes} genes, '
          f'consensus={info.mean_consensus:.3f}')
    # Check if cluster genes share a prefix (ground truth)
    prefixes = set(g.split('_')[0] for g in info.gene_symbols)
    print(f'    Prefixes: {prefixes}')

# Verify basic properties
assert result.n_clusters >= 2, f'FAIL: only {result.n_clusters} clusters'
assert result.consensus_matrix.shape == (40, 40)
clustered = sum(i.n_genes for i in result.cluster_info.values())
assert clustered + result.n_noise == 40
print()
print('PASS: phase3_clustering.py working correctly')
"
```

Regression:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_stats.py tests/test_ingestion.py -q 2>&1
```

Report all three outputs back.
