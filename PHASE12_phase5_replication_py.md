# INSTRUCTION SET FOR KAI — PHASE 12: `riker/phases/phase5_replication.py`

## References
- Blueprint Section 9 (Phase 5: Independent Replication)
- Blueprint Section 9.1 (Pre-Specification Protocol)
- Blueprint Section 9.2 (Replication Strategy)
- Blueprint Section 9.3 (Elimination Protocol)
- Blueprint Section 12 (QC Framework — "Pre-specification lock" and "Elimination protocol")

## WHY THIS MODULE MATTERS

Phase 5 is the anti-fishing protection. The core gene list from Phase 4
is LOCKED before any replication data is touched. No genes may be added
or removed based on replication results. Each core gene is then tested
in held-out datasets for directional concordance with discovery.

The elimination protocol removes any gene that is significantly (p < 0.05)
in the OPPOSITE direction in a brain replication dataset. This is what
caught KANK1 in the ASD proof-of-concept — it was eliminated for
significant discordance in brain replication (Blueprint Section 3.1).

Blood non-replication is EXPECTED for brain-specific signals and does
NOT trigger elimination (Blueprint Section 9.3).

## CRITICAL: PRE-SPECIFICATION PROTOCOL

Blueprint Section 9.1: "The core gene list is finalized and documented
BEFORE any replication dataset is downloaded or examined. No genes may
be added or removed based on replication results. No exploratory analysis
is performed on replication datasets."

The module enforces this by taking the locked core gene list as input
and only performing pre-specified directional tests.

## NAMING NOTE

Do NOT name any public function starting with "test_" to avoid
pytest collection conflicts (lesson from Phase 11).

## CRITICAL REQUIREMENTS

1. `replicate_gene()`: For one core gene in one replication dataset,
   compute log2FC, p-value, direction, and concordance with discovery.

2. `run_elimination_protocol()`: For each core gene across all replication
   datasets, apply the elimination rules:
   - ELIMINATE if significantly (p < 0.05) in OPPOSITE direction in brain
   - NOTE but don't eliminate non-significant opposite effects
   - Blood non-replication does NOT trigger elimination
   - Return elimination verdict per gene with full justification

3. `assign_cluster_verdicts()`: For each cluster, assign a verdict:
   - 'replicated': majority of core genes survive with concordant direction
   - 'partially_replicated': some genes survive
   - 'brain_specific': brain replicates but blood fails (expected)
   - 'failed': no genes survive or direction reversed

4. `run_phase5()`: Main entry point. Takes locked core genes + replication
   datasets + phenotypes. Returns Phase5Result with per-gene verdicts,
   eliminations, and cluster-level summaries.

5. DO NOT modify any existing files. APPEND tests to test_phases.py.

---

## FILE: `riker/phases/phase5_replication.py`

Write the following file at `/home/kai001/riker-engine/riker/phases/phase5_replication.py`:

```python
"""
Riker Engine - Phase 5: Independent Replication.

Tests core genes from Phase 4 in held-out replication datasets.
The core gene list is PRE-SPECIFIED and LOCKED before replication
data is accessed. No genes may be added or removed based on results.

The elimination protocol removes genes showing significant opposite
direction in brain replication datasets. Blood non-replication is
expected for brain-specific signals and does not trigger elimination.

References:
    Blueprint Section 9 (Phase 5: Independent Replication)
    Blueprint Section 12 (QC Framework — pre-specification, elimination)
"""

import logging
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from riker.stats.welch import welch_ttest

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ReplicationResult:
    """Replication test result for one gene in one dataset.

    Attributes:
        gene: Gene symbol.
        dataset_id: Replication dataset identifier.
        tissue: 'brain' or 'blood'.
        log2fc: Log2 fold change in replication dataset.
        p_value: P-value from Welch's t-test.
        direction: 'up' or 'down'.
        discovery_direction: Direction from discovery (Phase 1).
        is_concordant: True if replication direction matches discovery.
        is_significant: True if p < 0.05.
        n_cases: Number of case samples.
        n_controls: Number of control samples.
    """
    gene: str
    dataset_id: str
    tissue: str
    log2fc: float
    p_value: float
    direction: str
    discovery_direction: str
    is_concordant: bool
    is_significant: bool
    n_cases: int
    n_controls: int


@dataclass(frozen=True)
class GeneVerdict:
    """Elimination verdict for one core gene.

    Attributes:
        gene: Gene symbol.
        cluster_id: From Phase 4 core gene assignment.
        status: 'survived', 'eliminated', or 'insufficient_data'.
        reason: Human-readable explanation.
        replication_results: List of ReplicationResult across all datasets.
        n_brain_concordant: Brain datasets with concordant direction.
        n_brain_discordant: Brain datasets with significant opposite direction.
        n_blood_tested: Blood datasets tested.
        n_blood_concordant: Blood datasets with concordant direction.
        discovery_direction: Direction from discovery.
    """
    gene: str
    cluster_id: int
    status: str
    reason: str
    replication_results: list
    n_brain_concordant: int
    n_brain_discordant: int
    n_blood_tested: int
    n_blood_concordant: int
    discovery_direction: str


@dataclass(frozen=True)
class ClusterVerdict:
    """Replication verdict for one cluster.

    Attributes:
        cluster_id: Cluster label.
        verdict: 'replicated', 'partially_replicated', 'brain_specific', 'failed'.
        n_core_genes: Total core genes in this cluster.
        n_survived: Core genes that survived elimination.
        n_eliminated: Core genes eliminated.
        survived_genes: List of surviving gene symbols.
        eliminated_genes: List of eliminated gene symbols.
    """
    cluster_id: int
    verdict: str
    n_core_genes: int
    n_survived: int
    n_eliminated: int
    survived_genes: list
    eliminated_genes: list


@dataclass
class Phase5Result:
    """Complete result of Phase 5 replication.

    Attributes:
        gene_verdicts: Dict of gene_symbol -> GeneVerdict.
        cluster_verdicts: Dict of cluster_id -> ClusterVerdict.
        n_survived: Total genes surviving elimination.
        n_eliminated: Total genes eliminated.
        n_insufficient: Genes with insufficient replication data.
        locked_core_genes: The pre-specified core gene list (immutable record).
    """
    gene_verdicts: dict = field(default_factory=dict)
    cluster_verdicts: dict = field(default_factory=dict)
    n_survived: int = 0
    n_eliminated: int = 0
    n_insufficient: int = 0
    locked_core_genes: list = field(default_factory=list)


def replicate_gene(
    gene: str,
    discovery_direction: str,
    expression: pd.DataFrame,
    phenotype: dict[str, str],
    dataset_id: str,
    tissue: str = "brain",
) -> ReplicationResult | None:
    """Test one gene in one replication dataset.

    Parameters
    ----------
    gene : str
        Gene symbol to test.
    discovery_direction : str
        Direction from discovery ('up' or 'down').
    expression : DataFrame
        Gene-level expression (genes × samples).
    phenotype : dict
        Sample_id -> 'case' or 'control'.
    dataset_id : str
        Replication dataset identifier.
    tissue : str
        'brain' or 'blood'.

    Returns
    -------
    ReplicationResult or None
        None if gene is not in the dataset or insufficient samples.
    """
    if gene not in expression.index:
        return None

    case_samples = [s for s, g in phenotype.items()
                    if g == "case" and s in expression.columns]
    ctrl_samples = [s for s, g in phenotype.items()
                    if g == "control" and s in expression.columns]

    if len(case_samples) < 2 or len(ctrl_samples) < 2:
        return None

    case_vals = expression.loc[gene, case_samples].values.astype(np.float64)
    ctrl_vals = expression.loc[gene, ctrl_samples].values.astype(np.float64)

    case_vals = case_vals[np.isfinite(case_vals)]
    ctrl_vals = ctrl_vals[np.isfinite(ctrl_vals)]

    if len(case_vals) < 2 or len(ctrl_vals) < 2:
        return None

    try:
        result = welch_ttest(case_vals, ctrl_vals)
    except ValueError:
        return None

    direction = "up" if result.mean_diff > 0 else "down"
    is_concordant = (direction == discovery_direction)

    return ReplicationResult(
        gene=gene,
        dataset_id=dataset_id,
        tissue=tissue,
        log2fc=result.mean_diff,
        p_value=result.p_value,
        direction=direction,
        discovery_direction=discovery_direction,
        is_concordant=is_concordant,
        is_significant=(result.p_value < 0.05),
        n_cases=result.n1,
        n_controls=result.n2,
    )


def run_elimination_protocol(
    core_genes: dict,
    replication_datasets: dict[str, pd.DataFrame],
    replication_phenotypes: dict[str, dict[str, str]],
    dataset_tissues: dict[str, str],
    p_threshold: float = 0.05,
) -> dict[str, GeneVerdict]:
    """Apply elimination protocol to all core genes.

    Parameters
    ----------
    core_genes : dict
        From Phase4Result.core_genes (gene -> CoreGene).
    replication_datasets : dict
        Dataset_id -> gene-level expression DataFrame.
    replication_phenotypes : dict
        Dataset_id -> {sample_id: 'case' or 'control'}.
    dataset_tissues : dict
        Dataset_id -> 'brain' or 'blood'.
    p_threshold : float
        Significance threshold for elimination (default 0.05).

    Returns
    -------
    dict of gene_symbol -> GeneVerdict
    """
    verdicts = {}

    for gene, core_info in core_genes.items():
        discovery_dir = core_info.direction
        cluster_id = core_info.cluster_id
        rep_results = []

        # Test in each replication dataset
        for ds_id, expr_df in replication_datasets.items():
            tissue = dataset_tissues.get(ds_id, "brain")
            pheno = replication_phenotypes.get(ds_id, {})

            rep = replicate_gene(
                gene, discovery_dir, expr_df, pheno, ds_id, tissue
            )
            if rep is not None:
                rep_results.append(rep)

        # Apply elimination rules
        brain_results = [r for r in rep_results if r.tissue == "brain"]
        blood_results = [r for r in rep_results if r.tissue == "blood"]

        n_brain_concordant = sum(
            1 for r in brain_results if r.is_concordant
        )
        n_brain_discordant_sig = sum(
            1 for r in brain_results
            if not r.is_concordant and r.is_significant
        )
        n_blood_concordant = sum(
            1 for r in blood_results if r.is_concordant
        )

        # Elimination decision
        if not rep_results:
            status = "insufficient_data"
            reason = "No replication data available for this gene."
        elif n_brain_discordant_sig > 0:
            # ELIMINATE: significant opposite direction in brain
            discordant_ds = [
                r.dataset_id for r in brain_results
                if not r.is_concordant and r.is_significant
            ]
            status = "eliminated"
            reason = (
                f"Significant opposite direction in brain dataset(s): "
                f"{discordant_ds}. Discovery direction: {discovery_dir}. "
                f"See Blueprint Section 9.3."
            )
        else:
            status = "survived"
            if brain_results:
                reason = (
                    f"Brain replication: {n_brain_concordant}/{len(brain_results)} "
                    f"concordant. Blood: {n_blood_concordant}/{len(blood_results)} "
                    f"concordant."
                )
            else:
                reason = "No brain replication datasets; retained by default."

        verdicts[gene] = GeneVerdict(
            gene=gene,
            cluster_id=cluster_id,
            status=status,
            reason=reason,
            replication_results=rep_results,
            n_brain_concordant=n_brain_concordant,
            n_brain_discordant=n_brain_discordant_sig,
            n_blood_tested=len(blood_results),
            n_blood_concordant=n_blood_concordant,
            discovery_direction=discovery_dir,
        )

        logger.info(
            f"Gene {gene}: {status} — {reason}"
        )

    return verdicts


def assign_cluster_verdicts(
    gene_verdicts: dict[str, GeneVerdict],
    core_genes: dict,
) -> dict[int, ClusterVerdict]:
    """Assign replication verdicts per cluster.

    Parameters
    ----------
    gene_verdicts : dict
        From run_elimination_protocol().
    core_genes : dict
        From Phase4Result.core_genes.

    Returns
    -------
    dict of cluster_id -> ClusterVerdict
    """
    # Group genes by cluster
    cluster_genes: dict[int, list[str]] = {}
    for gene, core_info in core_genes.items():
        cid = core_info.cluster_id
        if cid not in cluster_genes:
            cluster_genes[cid] = []
        cluster_genes[cid].append(gene)

    verdicts = {}
    for cid, genes in cluster_genes.items():
        survived = []
        eliminated = []

        for gene in genes:
            if gene in gene_verdicts:
                gv = gene_verdicts[gene]
                if gv.status == "survived":
                    survived.append(gene)
                elif gv.status == "eliminated":
                    eliminated.append(gene)
                # insufficient_data genes are neither survived nor eliminated

        n_core = len(genes)
        n_surv = len(survived)
        n_elim = len(eliminated)

        # Determine cluster verdict
        if n_surv == 0:
            verdict = "failed"
        elif n_surv == n_core:
            # Check if blood failed (brain-specific pattern)
            has_blood_fail = False
            for gene in survived:
                gv = gene_verdicts[gene]
                if gv.n_blood_tested > 0 and gv.n_blood_concordant == 0:
                    has_blood_fail = True
                    break
            verdict = "brain_specific" if has_blood_fail else "replicated"
        elif n_surv >= n_core * 0.5:
            verdict = "partially_replicated"
        else:
            verdict = "failed"

        verdicts[cid] = ClusterVerdict(
            cluster_id=cid,
            verdict=verdict,
            n_core_genes=n_core,
            n_survived=n_surv,
            n_eliminated=n_elim,
            survived_genes=survived,
            eliminated_genes=eliminated,
        )

        logger.info(
            f"Cluster {cid}: {verdict} "
            f"({n_surv}/{n_core} survived, {n_elim} eliminated)"
        )

    return verdicts


def run_phase5(
    core_genes: dict,
    replication_datasets: dict[str, pd.DataFrame],
    replication_phenotypes: dict[str, dict[str, str]],
    dataset_tissues: dict[str, str],
) -> Phase5Result:
    """Run Phase 5: Independent Replication.

    Parameters
    ----------
    core_genes : dict
        From Phase4Result.core_genes. This list is LOCKED —
        no modifications allowed based on replication results.
    replication_datasets : dict
        Dataset_id -> gene-level expression DataFrame.
    replication_phenotypes : dict
        Dataset_id -> {sample_id: 'case' or 'control'}.
    dataset_tissues : dict
        Dataset_id -> 'brain' or 'blood'.

    Returns
    -------
    Phase5Result
    """
    # Record the locked gene list
    locked_list = sorted(core_genes.keys())

    logger.info(
        f"Phase 5: Testing {len(locked_list)} pre-specified core genes "
        f"in {len(replication_datasets)} replication datasets."
    )

    # Run elimination protocol
    gene_verdicts = run_elimination_protocol(
        core_genes, replication_datasets,
        replication_phenotypes, dataset_tissues,
    )

    # Assign cluster verdicts
    cluster_verdicts = assign_cluster_verdicts(gene_verdicts, core_genes)

    n_survived = sum(1 for v in gene_verdicts.values() if v.status == "survived")
    n_eliminated = sum(1 for v in gene_verdicts.values() if v.status == "eliminated")
    n_insufficient = sum(1 for v in gene_verdicts.values() if v.status == "insufficient_data")

    logger.info(
        f"Phase 5 complete: {n_survived} survived, {n_eliminated} eliminated, "
        f"{n_insufficient} insufficient data."
    )

    return Phase5Result(
        gene_verdicts=gene_verdicts,
        cluster_verdicts=cluster_verdicts,
        n_survived=n_survived,
        n_eliminated=n_eliminated,
        n_insufficient=n_insufficient,
        locked_core_genes=locked_list,
    )
```

---

## TESTS: APPEND to `tests/test_phases.py`

**APPEND** the following to the END of `/home/kai001/riker-engine/tests/test_phases.py`.
Do NOT delete or modify any existing test classes.

```python


# ===========================================================================
# PHASE 12: REPLICATION TESTS
# ===========================================================================

from riker.phases.phase5_replication import (
    ClusterVerdict,
    GeneVerdict,
    Phase5Result,
    ReplicationResult,
    assign_cluster_verdicts,
    replicate_gene,
    run_elimination_protocol,
    run_phase5,
)
from riker.phases.phase4_robustness import CoreGene


def _make_replication_data():
    """Create replication test data with known patterns.

    Returns core_genes, replication datasets, phenotypes, tissues.

    Core genes (from "discovery"):
    - GOOD_GENE_0..2: discovery direction "down", should replicate
    - BAD_GENE: discovery direction "down", replicates in OPPOSITE direction
    - BLOOD_FAIL_GENE: replicates in brain but fails in blood
    """
    # Core genes from Phase 4
    core_genes = {}
    for i in range(3):
        core_genes[f"GOOD_GENE_{i}"] = CoreGene(
            gene=f"GOOD_GENE_{i}", cluster_id=0,
            max_level_survived=3,
            per_dataset_pvalues={"DS1": 0.001, "DS2": 0.003},
            per_dataset_log2fc={"DS1": -0.8, "DS2": -0.7},
            mean_log2fc=-0.75, direction="down",
        )

    core_genes["BAD_GENE"] = CoreGene(
        gene="BAD_GENE", cluster_id=0,
        max_level_survived=2,
        per_dataset_pvalues={"DS1": 0.005},
        per_dataset_log2fc={"DS1": -0.6},
        mean_log2fc=-0.6, direction="down",
    )

    core_genes["BLOOD_FAIL_GENE"] = CoreGene(
        gene="BLOOD_FAIL_GENE", cluster_id=1,
        max_level_survived=2,
        per_dataset_pvalues={"DS1": 0.002},
        per_dataset_log2fc={"DS1": -0.9},
        mean_log2fc=-0.9, direction="down",
    )

    # Replication datasets
    np.random.seed(42)
    n_cases, n_controls = 12, 12
    n_total = n_cases + n_controls

    # Brain replication dataset
    brain_samples = [f"BS{i}" for i in range(n_total)]
    brain_data = {}
    for gene in ["GOOD_GENE_0", "GOOD_GENE_1", "GOOD_GENE_2"]:
        vals = np.random.normal(8.0, 0.5, n_total)
        vals[:n_cases] -= 0.8  # concordant downregulation
        brain_data[gene] = vals
    # BAD_GENE: opposite direction in brain (UP instead of down)
    vals = np.random.normal(8.0, 0.5, n_total)
    vals[:n_cases] += 1.2  # OPPOSITE direction
    brain_data["BAD_GENE"] = vals
    # BLOOD_FAIL_GENE: concordant in brain
    vals = np.random.normal(8.0, 0.5, n_total)
    vals[:n_cases] -= 0.7
    brain_data["BLOOD_FAIL_GENE"] = vals

    brain_expr = pd.DataFrame(brain_data, index=brain_samples).T
    brain_expr.index.name = "gene"
    brain_pheno = {s: ("case" if i < n_cases else "control")
                   for i, s in enumerate(brain_samples)}

    # Blood replication dataset
    blood_samples = [f"BL{i}" for i in range(n_total)]
    blood_data = {}
    for gene in ["GOOD_GENE_0", "GOOD_GENE_1", "GOOD_GENE_2", "BAD_GENE"]:
        blood_data[gene] = np.random.normal(8.0, 0.5, n_total)  # no effect
    # BLOOD_FAIL_GENE: no effect in blood (expected for brain-specific)
    blood_data["BLOOD_FAIL_GENE"] = np.random.normal(8.0, 0.5, n_total)

    blood_expr = pd.DataFrame(blood_data, index=blood_samples).T
    blood_expr.index.name = "gene"
    blood_pheno = {s: ("case" if i < n_cases else "control")
                   for i, s in enumerate(blood_samples)}

    replication_datasets = {
        "REP_BRAIN": brain_expr,
        "REP_BLOOD": blood_expr,
    }
    replication_phenotypes = {
        "REP_BRAIN": brain_pheno,
        "REP_BLOOD": blood_pheno,
    }
    dataset_tissues = {
        "REP_BRAIN": "brain",
        "REP_BLOOD": "blood",
    }

    return core_genes, replication_datasets, replication_phenotypes, dataset_tissues


# ---------------------------------------------------------------------------
# 14. Single gene replication
# ---------------------------------------------------------------------------

class TestReplicateGene:
    """Test per-gene replication testing."""

    def test_concordant_replication(self):
        core, datasets, phenos, tissues = _make_replication_data()
        result = replicate_gene(
            "GOOD_GENE_0", "down",
            datasets["REP_BRAIN"], phenos["REP_BRAIN"],
            "REP_BRAIN", "brain",
        )
        assert result is not None
        assert result.is_concordant is True
        assert result.direction == "down"
        assert result.tissue == "brain"

    def test_discordant_replication(self):
        core, datasets, phenos, tissues = _make_replication_data()
        result = replicate_gene(
            "BAD_GENE", "down",
            datasets["REP_BRAIN"], phenos["REP_BRAIN"],
            "REP_BRAIN", "brain",
        )
        assert result is not None
        assert result.is_concordant is False
        assert result.direction == "up"

    def test_missing_gene(self):
        core, datasets, phenos, tissues = _make_replication_data()
        result = replicate_gene(
            "NONEXISTENT_GENE", "down",
            datasets["REP_BRAIN"], phenos["REP_BRAIN"],
            "REP_BRAIN", "brain",
        )
        assert result is None

    def test_result_fields(self):
        core, datasets, phenos, tissues = _make_replication_data()
        result = replicate_gene(
            "GOOD_GENE_0", "down",
            datasets["REP_BRAIN"], phenos["REP_BRAIN"],
            "REP_BRAIN", "brain",
        )
        assert isinstance(result, ReplicationResult)
        assert result.gene == "GOOD_GENE_0"
        assert result.dataset_id == "REP_BRAIN"
        assert 0 < result.p_value <= 1
        assert result.n_cases > 0
        assert result.n_controls > 0


# ---------------------------------------------------------------------------
# 15. Elimination protocol
# ---------------------------------------------------------------------------

class TestEliminationProtocol:
    """Test the elimination protocol."""

    def test_bad_gene_eliminated(self):
        core, datasets, phenos, tissues = _make_replication_data()
        verdicts = run_elimination_protocol(
            core, datasets, phenos, tissues,
        )
        assert "BAD_GENE" in verdicts
        assert verdicts["BAD_GENE"].status == "eliminated"
        assert verdicts["BAD_GENE"].n_brain_discordant > 0

    def test_good_genes_survive(self):
        core, datasets, phenos, tissues = _make_replication_data()
        verdicts = run_elimination_protocol(
            core, datasets, phenos, tissues,
        )
        for i in range(3):
            gene = f"GOOD_GENE_{i}"
            assert gene in verdicts
            assert verdicts[gene].status == "survived"

    def test_blood_does_not_eliminate(self):
        """Blood non-replication should NOT trigger elimination."""
        core, datasets, phenos, tissues = _make_replication_data()
        verdicts = run_elimination_protocol(
            core, datasets, phenos, tissues,
        )
        # BLOOD_FAIL_GENE should survive (brain concordant, blood irrelevant)
        assert verdicts["BLOOD_FAIL_GENE"].status == "survived"

    def test_verdict_fields(self):
        core, datasets, phenos, tissues = _make_replication_data()
        verdicts = run_elimination_protocol(
            core, datasets, phenos, tissues,
        )
        for gene, v in verdicts.items():
            assert isinstance(v, GeneVerdict)
            assert v.discovery_direction in ("up", "down")
            assert v.status in ("survived", "eliminated", "insufficient_data")
            assert len(v.reason) > 0


# ---------------------------------------------------------------------------
# 16. Cluster verdicts
# ---------------------------------------------------------------------------

class TestClusterVerdicts:
    """Test cluster-level replication verdicts."""

    def test_cluster_with_elimination(self):
        core, datasets, phenos, tissues = _make_replication_data()
        gene_verdicts = run_elimination_protocol(
            core, datasets, phenos, tissues,
        )
        cluster_verdicts = assign_cluster_verdicts(gene_verdicts, core)

        # Cluster 0 has 3 good genes + 1 bad gene
        assert 0 in cluster_verdicts
        cv = cluster_verdicts[0]
        assert cv.n_eliminated >= 1
        assert "BAD_GENE" in cv.eliminated_genes

    def test_verdict_types(self):
        core, datasets, phenos, tissues = _make_replication_data()
        gene_verdicts = run_elimination_protocol(
            core, datasets, phenos, tissues,
        )
        cluster_verdicts = assign_cluster_verdicts(gene_verdicts, core)

        for cid, cv in cluster_verdicts.items():
            assert cv.verdict in (
                "replicated", "partially_replicated",
                "brain_specific", "failed"
            )
            assert cv.n_survived + cv.n_eliminated <= cv.n_core_genes


# ---------------------------------------------------------------------------
# 17. Full Phase 5 pipeline
# ---------------------------------------------------------------------------

class TestRunPhase5:
    """Test the integrated Phase 5 pipeline."""

    def test_full_run(self):
        core, datasets, phenos, tissues = _make_replication_data()
        result = run_phase5(core, datasets, phenos, tissues)

        assert isinstance(result, Phase5Result)
        assert result.n_survived >= 3  # at least the 3 good genes
        assert result.n_eliminated >= 1  # BAD_GENE
        assert len(result.locked_core_genes) == 5  # all 5 core genes recorded

    def test_locked_list_immutable(self):
        """The locked core gene list should be a record of what was tested."""
        core, datasets, phenos, tissues = _make_replication_data()
        result = run_phase5(core, datasets, phenos, tissues)

        # Locked list should contain all original core genes
        for gene in core:
            assert gene in result.locked_core_genes

    def test_elimination_matches_verdicts(self):
        core, datasets, phenos, tissues = _make_replication_data()
        result = run_phase5(core, datasets, phenos, tissues)

        n_surv = sum(1 for v in result.gene_verdicts.values() if v.status == "survived")
        n_elim = sum(1 for v in result.gene_verdicts.values() if v.status == "eliminated")
        assert n_surv == result.n_survived
        assert n_elim == result.n_eliminated
```

---

## EXECUTION INSTRUCTIONS

After writing phase5_replication.py and appending tests, run:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_phases.py -v 2>&1
```

**Expected: ALL tests pass.** If any fail, report FULL output — do NOT modify tests or code without reporting first.

Then confirmation:

```bash
cd /home/kai001/riker-engine && python -c "
from riker.phases.phase5_replication import run_phase5
from riker.phases.phase4_robustness import CoreGene
import numpy as np
import pandas as pd

# Simulate: 6 core genes, 2 brain + 1 blood replication datasets
np.random.seed(42)
core_genes = {}
for i in range(4):
    core_genes[f'SURVIVOR_{i}'] = CoreGene(
        gene=f'SURVIVOR_{i}', cluster_id=0, max_level_survived=3,
        per_dataset_pvalues={'D1': 0.001}, per_dataset_log2fc={'D1': -0.8},
        mean_log2fc=-0.8, direction='down',
    )
core_genes['ELIMINATED_1'] = CoreGene(
    gene='ELIMINATED_1', cluster_id=0, max_level_survived=2,
    per_dataset_pvalues={'D1': 0.005}, per_dataset_log2fc={'D1': -0.5},
    mean_log2fc=-0.5, direction='down',
)
core_genes['BRAIN_ONLY'] = CoreGene(
    gene='BRAIN_ONLY', cluster_id=1, max_level_survived=2,
    per_dataset_pvalues={'D1': 0.002}, per_dataset_log2fc={'D1': -0.9},
    mean_log2fc=-0.9, direction='down',
)

# Build replication expression data
n = 15
samples_b1 = [f'B1_{i}' for i in range(2*n)]
samples_b2 = [f'B2_{i}' for i in range(2*n)]
samples_bl = [f'BL_{i}' for i in range(2*n)]

def make_expr(samples, genes_shifts):
    data = {}
    for gene, shift in genes_shifts.items():
        vals = np.random.normal(8.0, 0.4, len(samples))
        vals[:n] += shift
        data[gene] = vals
    df = pd.DataFrame(data, index=samples).T
    df.index.name = 'gene'
    return df

brain1 = make_expr(samples_b1, {
    **{f'SURVIVOR_{i}': -0.7 for i in range(4)},
    'ELIMINATED_1': 1.0,  # OPPOSITE direction
    'BRAIN_ONLY': -0.6,
})
brain2 = make_expr(samples_b2, {
    **{f'SURVIVOR_{i}': -0.5 for i in range(4)},
    'ELIMINATED_1': 0.8,  # OPPOSITE direction again
    'BRAIN_ONLY': -0.5,
})
blood = make_expr(samples_bl, {
    **{f'SURVIVOR_{i}': 0.0 for i in range(4)},
    'ELIMINATED_1': 0.0,
    'BRAIN_ONLY': 0.0,  # no effect in blood (brain-specific)
})

pheno_b1 = {s: ('case' if i < n else 'control') for i, s in enumerate(samples_b1)}
pheno_b2 = {s: ('case' if i < n else 'control') for i, s in enumerate(samples_b2)}
pheno_bl = {s: ('case' if i < n else 'control') for i, s in enumerate(samples_bl)}

result = run_phase5(
    core_genes,
    {'BRAIN1': brain1, 'BRAIN2': brain2, 'BLOOD1': blood},
    {'BRAIN1': pheno_b1, 'BRAIN2': pheno_b2, 'BLOOD1': pheno_bl},
    {'BRAIN1': 'brain', 'BRAIN2': 'brain', 'BLOOD1': 'blood'},
)

print('=== Phase 5 Replication ===')
print(f'Locked core genes: {result.locked_core_genes}')
print(f'Survived: {result.n_survived}')
print(f'Eliminated: {result.n_eliminated}')
print()
for gene, v in sorted(result.gene_verdicts.items()):
    print(f'  {gene}: {v.status} (brain_concordant={v.n_brain_concordant}, '
          f'brain_discordant={v.n_brain_discordant}, blood={v.n_blood_concordant}/{v.n_blood_tested})')
print()
for cid, cv in sorted(result.cluster_verdicts.items()):
    print(f'  Cluster {cid}: {cv.verdict} ({cv.n_survived}/{cv.n_core_genes} survived)')

assert result.n_eliminated >= 1, 'FAIL: ELIMINATED_1 should be eliminated'
assert result.n_survived >= 4, 'FAIL: survivors too few'
assert result.gene_verdicts['ELIMINATED_1'].status == 'eliminated'
for i in range(4):
    assert result.gene_verdicts[f'SURVIVOR_{i}'].status == 'survived'
print()
print('PASS: phase5_replication.py working correctly')
"
```

Regression:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_stats.py tests/test_ingestion.py -q 2>&1
```

Report all three outputs back.
