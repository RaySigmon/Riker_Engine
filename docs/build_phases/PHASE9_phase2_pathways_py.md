# INSTRUCTION SET FOR KAI — PHASE 9: `riker/phases/phase2_pathways.py`

## References
- Blueprint Section 6 (Phase 2: Pathway Mapping)
- Blueprint Section 7 (Phase 3 — feature vector description)
- Blueprint Section 12 (QC Framework — "Circularity audit" row)
- Blueprint Section 14 (Technology Stack — KEGG REST API, Reactome bulk download)

## WHY THIS MODULE MATTERS

Phase 2 transforms the study gene set from Phase 1 into a feature matrix
suitable for unsupervised clustering. Each gene gets a feature vector
combining binary pathway memberships, expression statistics, and
confidence tier scores. This feature matrix is what Phase 3's
UMAP/HDBSCAN operates on.

The pathway sources are KEGG (~350 human pathways), Reactome (~2200),
and MSigDB Hallmark (50 curated sets). Pathways are filtered by gene
membership criteria to remove overly broad or sparsely populated sets.

## CRITICAL REQUIREMENT: ANTI-CIRCULARITY

Blueprint Section 7 red box: "Do NOT use pre-assigned biological category
labels (e.g., 'chromatin remodeling') as clustering features. Use only
individual pathway IDs as binary features. Using categories as both
features and enrichment targets is tautological."

## DESIGN FOR OFFLINE/TESTABLE OPERATION

Since Ghost (Raspberry Pi 5) may not always have network access, and
since tests MUST run without network, the module must:
1. Accept pre-loaded pathway data (dicts/DataFrames) — no API calls in core logic
2. Provide separate loader functions that CAN make API calls
3. Tests use small fixture pathway data, never real KEGG/Reactome

## CRITICAL REQUIREMENTS

1. `PathwayDatabase` class that holds pathway-to-gene mappings from any source.
   Simple structure: dict of pathway_id -> set of gene symbols.

2. `load_pathways_from_dict()`: Build PathwayDatabase from a plain dict.
   This is what tests and offline mode use.

3. `filter_pathways()`: Apply Blueprint Section 6 filters:
   - Minimum 3 study gene members per pathway
   - Maximum 500 total genes per pathway (excludes overly broad)
   - Minimum 2% study gene representation
   Returns filtered PathwayDatabase.

4. `build_feature_matrix()`: Construct the feature matrix per Blueprint Section 7:
   - Binary pathway memberships (one column per pathway)
   - Average log2FC across datasets
   - Negative log10 of minimum p-value (capped at 10)
   - Number of significant datasets
   - Directional consistency score (fraction of datasets agreeing)
   - Confidence tier score (inverted so higher = more confident)
   - All features min-max normalized to [0, 1]

5. Return a `Phase2Result` dataclass with the feature matrix, pathway
   annotations, and filtering statistics.

6. DO NOT modify any existing files. APPEND new tests to test_phases.py.

---

## FILE: `riker/phases/phase2_pathways.py`

Write the following file at `/home/kai001/riker-engine/riker/phases/phase2_pathways.py`:

```python
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
```

---

## TESTS: APPEND to `tests/test_phases.py`

**APPEND** the following to the END of `/home/kai001/riker-engine/tests/test_phases.py`.
Do NOT delete or modify any existing test classes.

```python


# ===========================================================================
# PHASE 9: PATHWAY MAPPING TESTS
# ===========================================================================

from riker.phases.phase2_pathways import (
    PathwayDatabase,
    PathwayInfo,
    Phase2Result,
    build_feature_matrix,
    filter_pathways,
    load_pathways_from_dict,
    run_phase2,
)


def _make_pathway_db():
    """Create a test pathway database with known structure."""
    data = {
        "PATH_001": ["GENE_0", "GENE_1", "GENE_2", "GENE_3", "GENE_4",
                      "GENE_5", "GENE_6", "GENE_7", "GENE_8", "GENE_9"],
        "PATH_002": ["GENE_0", "GENE_1", "GENE_2", "GENE_10", "GENE_11"],
        "PATH_003": ["GENE_3", "GENE_4", "GENE_5", "GENE_12", "GENE_13"],
        "PATH_004": ["GENE_0", "GENE_1"],  # too few study genes (< 3)
        "PATH_005": [f"BIG_{i}" for i in range(600)],  # too big (> 500)
        "PATH_006": ["GENE_0", "GENE_1", "GENE_2", "GENE_3"],
    }
    names = {
        "PATH_001": "Calcium signaling",
        "PATH_002": "Synaptic transmission",
        "PATH_003": "Cell adhesion",
        "PATH_004": "Tiny pathway",
        "PATH_005": "Huge pathway",
        "PATH_006": "Neuronal development",
    }
    return load_pathways_from_dict(data, names=names, source="TEST")


def _make_study_genes_dict(n_genes=10):
    """Create a mock study_genes dict mimicking Phase1Result.study_genes."""
    from riker.phases.phase1_crossref import GeneResult, GeneDatasetDE
    study = {}
    for i in range(n_genes):
        gene = f"GENE_{i}"
        de1 = GeneDatasetDE(
            gene=gene, dataset_id="DS1", log2fc=-0.5 - i * 0.1,
            p_value=0.001 + i * 0.005, t_statistic=-3.0,
            df=28.0, se=0.15, n_cases=15, n_controls=15,
            direction="down",
        )
        de2 = GeneDatasetDE(
            gene=gene, dataset_id="DS2", log2fc=-0.4 - i * 0.08,
            p_value=0.002 + i * 0.006, t_statistic=-2.5,
            df=22.0, se=0.18, n_cases=12, n_controls=12,
            direction="down",
        )
        study[gene] = GeneResult(
            gene=gene, de_results=[de1, de2],
            n_datasets_detected=2, n_datasets_significant=2,
            passes_filter=True, mean_log2fc=-0.45 - i * 0.09,
            consistent_direction=True,
        )
    return study


# ---------------------------------------------------------------------------
# 3. PathwayDatabase operations
# ---------------------------------------------------------------------------

class TestPathwayDatabase:
    """Test PathwayDatabase construction and operations."""

    def test_load_from_dict(self):
        db = _make_pathway_db()
        assert db.n_pathways == 6
        assert db.source == "TEST"

    def test_get_gene_pathways(self):
        db = _make_pathway_db()
        paths = db.get_gene_pathways("GENE_0")
        assert "PATH_001" in paths
        assert "PATH_002" in paths
        assert "PATH_004" in paths

    def test_all_genes(self):
        db = _make_pathway_db()
        all_g = db.all_genes
        assert "GENE_0" in all_g
        assert "GENE_13" in all_g

    def test_merge(self):
        db1 = load_pathways_from_dict({"P1": ["A", "B"]}, source="KEGG")
        db2 = load_pathways_from_dict({"P2": ["C", "D"]}, source="Reactome")
        db1.merge(db2)
        assert db1.n_pathways == 2
        assert db1.source == "mixed"

    def test_pathway_names(self):
        db = _make_pathway_db()
        assert db.names["PATH_001"] == "Calcium signaling"


# ---------------------------------------------------------------------------
# 4. Pathway filtering
# ---------------------------------------------------------------------------

class TestFilterPathways:
    """Test pathway filtering by study gene criteria."""

    def test_basic_filter(self):
        db = _make_pathway_db()
        study_genes = [f"GENE_{i}" for i in range(10)]
        filtered, info = filter_pathways(db, study_genes)

        # PATH_004 should be filtered (only 2 study genes < 3 minimum)
        assert "PATH_004" not in filtered.pathways
        # PATH_005 should be filtered (600 genes > 500 max)
        assert "PATH_005" not in filtered.pathways
        # PATH_001, PATH_002, PATH_003, PATH_006 should survive
        assert "PATH_001" in filtered.pathways

    def test_info_populated(self):
        db = _make_pathway_db()
        study_genes = [f"GENE_{i}" for i in range(10)]
        _, info = filter_pathways(db, study_genes)

        assert "PATH_001" in info
        p1 = info["PATH_001"]
        assert p1.source == "TEST"
        assert p1.total_genes == 10
        assert p1.study_genes > 0
        assert len(p1.study_gene_list) == p1.study_genes

    def test_max_pathways_limit(self):
        db = _make_pathway_db()
        study_genes = [f"GENE_{i}" for i in range(10)]
        filtered, _ = filter_pathways(db, study_genes, max_pathways=2)
        assert filtered.n_pathways <= 2

    def test_empty_study_genes(self):
        db = _make_pathway_db()
        filtered, info = filter_pathways(db, [])
        assert filtered.n_pathways == 0


# ---------------------------------------------------------------------------
# 5. Feature matrix construction
# ---------------------------------------------------------------------------

class TestBuildFeatureMatrix:
    """Test feature matrix construction for clustering."""

    def test_basic_matrix(self):
        db = _make_pathway_db()
        study = _make_study_genes_dict(10)
        study_genes = list(study.keys())
        filtered, info = filter_pathways(db, study_genes)

        result = build_feature_matrix(study, filtered, info)

        assert result.n_genes == 10
        assert result.n_features > 0
        assert result.feature_matrix.shape[0] == 10

    def test_normalized_range(self):
        """All features should be in [0, 1] after min-max normalization."""
        db = _make_pathway_db()
        study = _make_study_genes_dict(10)
        study_genes = list(study.keys())
        filtered, info = filter_pathways(db, study_genes)

        result = build_feature_matrix(study, filtered, info)

        for col in result.feature_matrix.columns:
            assert result.feature_matrix[col].min() >= -1e-10
            assert result.feature_matrix[col].max() <= 1.0 + 1e-10

    def test_pathway_features_binary_before_normalization(self):
        """Pathway columns should be binary (0 or 1) before normalization.
        After normalization they may still be 0/1 or rescaled."""
        db = _make_pathway_db()
        study = _make_study_genes_dict(10)
        study_genes = list(study.keys())
        filtered, info = filter_pathways(db, study_genes)

        result = build_feature_matrix(study, filtered, info)

        # Check that pathway features exist
        assert len(result.pathway_features) > 0
        for pid in result.pathway_features:
            if pid in result.feature_matrix.columns:
                vals = result.feature_matrix[pid].unique()
                # After normalization, should have at most 2 unique values
                assert len(vals) <= 2

    def test_expression_features_present(self):
        db = _make_pathway_db()
        study = _make_study_genes_dict(10)
        study_genes = list(study.keys())
        filtered, info = filter_pathways(db, study_genes)

        result = build_feature_matrix(study, filtered, info)

        assert "avg_log2fc" in result.expression_feature_names
        assert "neg_log10_min_p" in result.expression_feature_names
        assert "n_sig_datasets" in result.expression_feature_names
        assert "direction_consistency" in result.expression_feature_names
        assert "tier_score" in result.expression_feature_names

    def test_with_tier_scores(self):
        db = _make_pathway_db()
        study = _make_study_genes_dict(10)
        study_genes = list(study.keys())
        filtered, info = filter_pathways(db, study_genes)

        tiers = {f"GENE_{i}": float(3 - (i % 3)) for i in range(10)}
        result = build_feature_matrix(study, filtered, info, gene_tiers=tiers)

        # Tier scores should vary (not all zero)
        tier_col = result.feature_matrix["tier_score"]
        assert tier_col.max() > tier_col.min()

    def test_empty_study_genes(self):
        db = _make_pathway_db()
        filtered, info = filter_pathways(db, [])
        result = build_feature_matrix({}, filtered, info)
        assert result.n_genes == 0


# ---------------------------------------------------------------------------
# 6. Full Phase 2 run
# ---------------------------------------------------------------------------

class TestRunPhase2:
    """Test the integrated Phase 2 pipeline."""

    def test_integration(self):
        """Run Phase 2 with Phase 1 output and pathway DB."""
        from riker.phases.phase1_crossref import Phase1Result

        # Build a mock Phase1Result
        study = _make_study_genes_dict(10)
        phase1 = Phase1Result(
            study_genes=study,
            n_seed_genes=50,
            n_study_genes=10,
        )

        db = _make_pathway_db()
        result = run_phase2(phase1, db)

        assert result.n_genes == 10
        assert result.n_pathways_after_filter > 0
        assert result.n_features > 0
        assert result.feature_matrix.shape[0] == 10

    def test_pathway_count_tracking(self):
        from riker.phases.phase1_crossref import Phase1Result

        study = _make_study_genes_dict(10)
        phase1 = Phase1Result(study_genes=study, n_study_genes=10)

        db = _make_pathway_db()
        result = run_phase2(phase1, db)

        assert result.n_pathways_before_filter == 6  # our test DB has 6
        assert result.n_pathways_after_filter <= 6
```

---

## EXECUTION INSTRUCTIONS

After writing phase2_pathways.py and appending tests, run:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_phases.py -v 2>&1
```

**Expected: ALL tests pass (Phase 8 + Phase 9).** If any fail, report FULL output — do NOT modify tests or code without reporting first.

Then confirmation:

```bash
cd /home/kai001/riker-engine && python -c "
from riker.phases.phase1_crossref import Phase1Result, GeneResult, GeneDatasetDE
from riker.phases.phase2_pathways import load_pathways_from_dict, run_phase2
import numpy as np

# Build mock Phase 1 output (20 study genes)
study = {}
for i in range(20):
    gene = f'GENE_{i}'
    de1 = GeneDatasetDE(gene=gene, dataset_id='DS1', log2fc=-0.5-i*0.05,
                        p_value=0.001+i*0.002, t_statistic=-3.0, df=28.0,
                        se=0.15, n_cases=15, n_controls=15, direction='down')
    de2 = GeneDatasetDE(gene=gene, dataset_id='DS2', log2fc=-0.4-i*0.04,
                        p_value=0.002+i*0.003, t_statistic=-2.5, df=22.0,
                        se=0.18, n_cases=12, n_controls=12, direction='down')
    study[gene] = GeneResult(gene=gene, de_results=[de1, de2],
                             n_datasets_detected=2, n_datasets_significant=2,
                             passes_filter=True, mean_log2fc=-0.45-i*0.045,
                             consistent_direction=True)

phase1 = Phase1Result(study_genes=study, n_seed_genes=100, n_study_genes=20)

# Build pathway DB with 10 pathways
pathways = {}
for p in range(10):
    genes = [f'GENE_{g}' for g in range(p, min(p+8, 20))]
    pathways[f'hsa{p:05d}'] = genes

db = load_pathways_from_dict(pathways, source='KEGG')

# Run Phase 2
result = run_phase2(phase1, db)

print(f'=== Phase 2 Pathway Mapping ===')
print(f'Genes: {result.n_genes}')
print(f'Features: {result.n_features}')
print(f'Pathways before filter: {result.n_pathways_before_filter}')
print(f'Pathways after filter: {result.n_pathways_after_filter}')
print(f'Expression features: {result.expression_feature_names}')
print(f'Feature matrix shape: {result.feature_matrix.shape}')
print(f'Value range: [{result.feature_matrix.min().min():.4f}, {result.feature_matrix.max().max():.4f}]')

assert result.n_genes == 20
assert result.n_features > 5  # at least expression features
assert result.feature_matrix.min().min() >= -1e-10
assert result.feature_matrix.max().max() <= 1.0 + 1e-10
print()
print('PASS: phase2_pathways.py working correctly')
"
```

Regression:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_stats.py tests/test_ingestion.py -q 2>&1
```

Report all three outputs back.
