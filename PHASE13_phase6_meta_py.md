# INSTRUCTION SET FOR KAI — PHASE 13: `riker/phases/phase6_meta.py`

## References
- Blueprint Section 10 (Phase 6: Effect Size Meta-Analysis)
- Blueprint Section 10.1 (Inverse-Variance Weighted fixed/random effects)
- Blueprint Section 12 (QC Framework — "Expression scale check" row)
- Context Transfer: "GSE33000 log-ratio scale" note

## WHY THIS MODULE MATTERS

Phase 6 computes formal effect size estimates for each surviving gene
using inverse-variance weighted meta-analysis across all discovery
datasets. This produces the forest plot data and heterogeneity statistics
that appear in the final report.

The module wires the existing `riker.stats.meta` functions (already built
and locked in Phase 3 of the stats layer) into the pipeline phase structure.

## CRITICAL: EXPRESSION SCALE CHECK

Blueprint Section 10.1 and Context Transfer note: GSE33000 expression
values are log-ratios centered near zero (max ~2.0), not log2-intensities.
Effect sizes from this dataset are on a DIFFERENT scale. The module must
call `check_expression_scale()` from `riker.stats.meta` before computing
meta-analysis, and flag datasets that need SE adjustment.

## CRITICAL REQUIREMENTS

1. `compute_gene_meta()`: For one surviving gene, collect per-dataset
   effect sizes (log2FC) and standard errors, then run IVW meta-analysis
   (both fixed and random effects). Returns forest plot data.

2. `run_phase6()`: Main entry point. For each gene surviving Phase 5,
   compute meta-analysis. Returns Phase6Result with per-gene synthesis,
   heterogeneity stats, and forest plot data.

3. SE recovery: When SE is not directly available, recover from t-statistic
   and sample sizes using `riker.stats.meta.recover_se()`.

4. Expression scale check: Call `check_expression_scale()` per dataset
   and warn if any dataset appears to be on a different scale.

5. DO NOT modify any existing files. APPEND tests to test_phases.py.

---

## FILE: `riker/phases/phase6_meta.py`

Write the following file at `/home/kai001/riker-engine/riker/phases/phase6_meta.py`:

```python
"""
Riker Engine - Phase 6: Effect Size Meta-Analysis.

Computes inverse-variance weighted meta-analysis for each surviving
gene across all discovery datasets. Produces forest plot data and
heterogeneity statistics.

References:
    Blueprint Section 10 (Phase 6: Effect Size Meta-Analysis)
    Blueprint Section 10.1 (IVW fixed/random effects)
    Context Transfer: GSE33000 scale check
"""

import logging
import warnings
from dataclasses import dataclass, field

import numpy as np

from riker.stats.meta import (
    ivw_meta_analysis,
    recover_se,
    check_expression_scale,
    MetaAnalysisResult,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class GeneEffectSize:
    """Per-dataset effect size for one gene.

    Attributes:
        dataset_id: Dataset identifier.
        log2fc: Log2 fold change (effect size).
        se: Standard error of the effect size.
        p_value: P-value from Welch's t-test.
        n_cases: Number of case samples.
        n_controls: Number of control samples.
        scale_warning: True if dataset has expression scale issues.
    """
    dataset_id: str
    log2fc: float
    se: float
    p_value: float
    n_cases: int
    n_controls: int
    scale_warning: bool


@dataclass(frozen=True)
class GeneMetaResult:
    """Meta-analysis result for one gene.

    Attributes:
        gene: Gene symbol.
        cluster_id: From Phase 4/5 assignment.
        per_dataset: List of GeneEffectSize (forest plot rows).
        fixed_effect: Fixed-effect pooled estimate.
        fixed_se: Fixed-effect standard error.
        fixed_p: Fixed-effect p-value.
        random_effect: Random-effects pooled estimate (DerSimonian-Laird).
        random_se: Random-effects standard error.
        random_p: Random-effects p-value.
        cochran_q: Cochran's Q heterogeneity statistic.
        i_squared: I² heterogeneity percentage.
        tau_squared: Between-study variance (τ²).
        n_datasets: Number of datasets in meta-analysis.
        direction: Overall direction from random-effects estimate.
    """
    gene: str
    cluster_id: int
    per_dataset: list
    fixed_effect: float
    fixed_se: float
    fixed_p: float
    random_effect: float
    random_se: float
    random_p: float
    cochran_q: float
    i_squared: float
    tau_squared: float
    n_datasets: int
    direction: str


@dataclass
class Phase6Result:
    """Complete result of Phase 6 meta-analysis.

    Attributes:
        gene_results: Dict of gene_symbol -> GeneMetaResult.
        n_genes_analyzed: Number of genes with meta-analysis.
        n_significant_random: Genes significant (p<0.05) under random effects.
        n_high_heterogeneity: Genes with I² > 75%.
        scale_warnings: List of dataset IDs with expression scale issues.
    """
    gene_results: dict = field(default_factory=dict)
    n_genes_analyzed: int = 0
    n_significant_random: int = 0
    n_high_heterogeneity: int = 0
    scale_warnings: list = field(default_factory=list)


def compute_gene_meta(
    gene: str,
    cluster_id: int,
    de_results: list,
    dataset_scales: dict[str, bool] | None = None,
) -> GeneMetaResult | None:
    """Compute meta-analysis for one gene.

    Parameters
    ----------
    gene : str
        Gene symbol.
    cluster_id : int
        Cluster assignment.
    de_results : list
        List of GeneDatasetDE from Phase 1 (per-dataset stats).
    dataset_scales : dict or None
        Dataset_id -> True if scale warning applies.

    Returns
    -------
    GeneMetaResult or None
        None if fewer than 2 datasets available.
    """
    if len(de_results) < 2:
        logger.warning(
            f"Gene {gene}: only {len(de_results)} dataset(s), "
            f"need >= 2 for meta-analysis. Skipping."
        )
        return None

    effects = []
    ses = []
    per_dataset = []

    for de in de_results:
        effect = de.log2fc
        se = de.se

        # If SE looks invalid, try to recover from t-stat and sample sizes
        if se <= 0 or not np.isfinite(se):
            try:
                se = recover_se(
                    t_statistic=de.t_statistic,
                    n1=de.n_cases,
                    n2=de.n_controls,
                )
            except (ValueError, ZeroDivisionError):
                logger.warning(
                    f"Gene {gene}, dataset {de.dataset_id}: "
                    f"could not recover SE. Skipping dataset."
                )
                continue

        scale_warn = False
        if dataset_scales and de.dataset_id in dataset_scales:
            scale_warn = dataset_scales[de.dataset_id]

        effects.append(effect)
        ses.append(se)
        per_dataset.append(GeneEffectSize(
            dataset_id=de.dataset_id,
            log2fc=effect,
            se=se,
            p_value=de.p_value,
            n_cases=de.n_cases,
            n_controls=de.n_controls,
            scale_warning=scale_warn,
        ))

    if len(effects) < 2:
        return None

    # Run IVW meta-analysis
    meta = ivw_meta_analysis(
        effects=effects,
        standard_errors=ses,
    )

    direction = "down" if meta.random_effect < 0 else "up"

    return GeneMetaResult(
        gene=gene,
        cluster_id=cluster_id,
        per_dataset=per_dataset,
        fixed_effect=meta.fixed_effect,
        fixed_se=meta.fixed_se,
        fixed_p=meta.fixed_p,
        random_effect=meta.random_effect,
        random_se=meta.random_se,
        random_p=meta.random_p,
        cochran_q=meta.cochran_q,
        i_squared=meta.i_squared,
        tau_squared=meta.tau_squared,
        n_datasets=len(effects),
        direction=direction,
    )


def run_phase6(
    phase1_result,
    phase5_result,
    dataset_expression_ranges: dict[str, float] | None = None,
) -> Phase6Result:
    """Run Phase 6 meta-analysis for all surviving genes.

    Parameters
    ----------
    phase1_result : Phase1Result
        Contains study_genes with per-dataset DE stats.
    phase5_result : Phase5Result
        Contains gene_verdicts (survived/eliminated) and locked core genes.
    dataset_expression_ranges : dict or None
        Dataset_id -> max expression value (for scale check).
        If None, scale check is skipped.

    Returns
    -------
    Phase6Result
    """
    result = Phase6Result()

    # Check expression scales
    dataset_scales = {}
    if dataset_expression_ranges:
        for ds_id, max_val in dataset_expression_ranges.items():
            is_ok, msg = check_expression_scale(max_val)
            if not is_ok:
                dataset_scales[ds_id] = True
                result.scale_warnings.append(f"{ds_id}: {msg}")
                logger.warning(f"Scale check: {ds_id}: {msg}")
            else:
                dataset_scales[ds_id] = False

    # Get surviving genes from Phase 5
    surviving_genes = [
        gene for gene, verdict in phase5_result.gene_verdicts.items()
        if verdict.status == "survived"
    ]

    logger.info(
        f"Phase 6: Computing meta-analysis for {len(surviving_genes)} "
        f"surviving genes."
    )

    for gene in surviving_genes:
        # Get Phase 1 DE results for this gene
        if gene not in phase1_result.study_genes:
            continue

        gene_result = phase1_result.study_genes[gene]
        cluster_id = phase5_result.gene_verdicts[gene].cluster_id

        meta = compute_gene_meta(
            gene, cluster_id, gene_result.de_results,
            dataset_scales=dataset_scales or None,
        )

        if meta is not None:
            result.gene_results[gene] = meta

    result.n_genes_analyzed = len(result.gene_results)
    result.n_significant_random = sum(
        1 for m in result.gene_results.values() if m.random_p < 0.05
    )
    result.n_high_heterogeneity = sum(
        1 for m in result.gene_results.values() if m.i_squared > 75.0
    )

    logger.info(
        f"Phase 6 complete: {result.n_genes_analyzed} genes analyzed, "
        f"{result.n_significant_random} significant (random effects), "
        f"{result.n_high_heterogeneity} high heterogeneity."
    )

    return result
```

---

## TESTS: APPEND to `tests/test_phases.py`

**APPEND** the following to the END of `/home/kai001/riker-engine/tests/test_phases.py`.
Do NOT delete or modify any existing test classes.

```python


# ===========================================================================
# PHASE 13: META-ANALYSIS PIPELINE TESTS
# ===========================================================================

from riker.phases.phase6_meta import (
    GeneEffectSize,
    GeneMetaResult,
    Phase6Result,
    compute_gene_meta,
    run_phase6,
)


def _make_meta_test_data():
    """Create test data for Phase 6 meta-analysis.

    Returns phase1_result, phase5_result with surviving genes
    that have multi-dataset DE results.
    """
    from riker.phases.phase1_crossref import Phase1Result, GeneResult, GeneDatasetDE
    from riker.phases.phase5_replication import Phase5Result, GeneVerdict

    # Phase 1 study genes with per-dataset DE
    study_genes = {}
    for i in range(5):
        gene = f"META_GENE_{i}"
        des = []
        for j, ds in enumerate(["DS1", "DS2", "DS3"]):
            des.append(GeneDatasetDE(
                gene=gene, dataset_id=ds,
                log2fc=-0.6 - i * 0.1 + j * 0.05,
                p_value=0.001 + i * 0.002,
                t_statistic=-3.5 + j * 0.3,
                df=28.0, se=0.15 + j * 0.02,
                n_cases=15, n_controls=15,
                direction="down",
            ))
        study_genes[gene] = GeneResult(
            gene=gene, de_results=des,
            n_datasets_detected=3, n_datasets_significant=3,
            passes_filter=True,
            mean_log2fc=float(np.mean([d.log2fc for d in des])),
            consistent_direction=True,
        )

    phase1 = Phase1Result(
        study_genes=study_genes,
        n_seed_genes=100,
        n_study_genes=5,
    )

    # Phase 5 verdicts: all survived
    gene_verdicts = {}
    for i in range(5):
        gene = f"META_GENE_{i}"
        gene_verdicts[gene] = GeneVerdict(
            gene=gene, cluster_id=0,
            status="survived",
            reason="Brain concordant.",
            replication_results=[],
            n_brain_concordant=2, n_brain_discordant=0,
            n_blood_tested=1, n_blood_concordant=0,
            discovery_direction="down",
        )

    phase5 = Phase5Result(
        gene_verdicts=gene_verdicts,
        n_survived=5,
        locked_core_genes=[f"META_GENE_{i}" for i in range(5)],
    )

    return phase1, phase5


# ---------------------------------------------------------------------------
# 18. Single gene meta-analysis
# ---------------------------------------------------------------------------

class TestComputeGeneMeta:
    """Test per-gene meta-analysis computation."""

    def test_basic_meta(self):
        from riker.phases.phase1_crossref import GeneDatasetDE
        des = [
            GeneDatasetDE(
                gene="G", dataset_id=f"DS{i}",
                log2fc=-0.5 - i * 0.1, p_value=0.01,
                t_statistic=-3.0, df=28.0, se=0.15,
                n_cases=15, n_controls=15, direction="down",
            )
            for i in range(3)
        ]
        result = compute_gene_meta("G", 0, des)
        assert result is not None
        assert result.n_datasets == 3
        assert result.random_effect < 0  # downregulated
        assert result.direction == "down"
        assert 0 <= result.i_squared <= 100

    def test_insufficient_datasets(self):
        from riker.phases.phase1_crossref import GeneDatasetDE
        des = [GeneDatasetDE(
            gene="G", dataset_id="DS1",
            log2fc=-0.5, p_value=0.01,
            t_statistic=-3.0, df=28.0, se=0.15,
            n_cases=15, n_controls=15, direction="down",
        )]
        result = compute_gene_meta("G", 0, des)
        assert result is None

    def test_forest_plot_data(self):
        from riker.phases.phase1_crossref import GeneDatasetDE
        des = [
            GeneDatasetDE(
                gene="G", dataset_id=f"DS{i}",
                log2fc=-0.5, p_value=0.01,
                t_statistic=-3.0, df=28.0, se=0.15,
                n_cases=15, n_controls=15, direction="down",
            )
            for i in range(3)
        ]
        result = compute_gene_meta("G", 0, des)
        assert len(result.per_dataset) == 3
        for row in result.per_dataset:
            assert isinstance(row, GeneEffectSize)
            assert row.se > 0

    def test_scale_warning_flagged(self):
        from riker.phases.phase1_crossref import GeneDatasetDE
        des = [
            GeneDatasetDE(
                gene="G", dataset_id=f"DS{i}",
                log2fc=-0.5, p_value=0.01,
                t_statistic=-3.0, df=28.0, se=0.15,
                n_cases=15, n_controls=15, direction="down",
            )
            for i in range(3)
        ]
        scales = {"DS0": True, "DS1": False, "DS2": False}
        result = compute_gene_meta("G", 0, des, dataset_scales=scales)
        # DS0 should have scale_warning=True
        assert result.per_dataset[0].scale_warning is True
        assert result.per_dataset[1].scale_warning is False

    def test_heterogeneity_stats(self):
        from riker.phases.phase1_crossref import GeneDatasetDE
        des = [
            GeneDatasetDE(
                gene="G", dataset_id=f"DS{i}",
                log2fc=-0.5 + i * 0.8,  # very different effects
                p_value=0.01,
                t_statistic=-3.0, df=28.0, se=0.15,
                n_cases=15, n_controls=15, direction="down",
            )
            for i in range(4)
        ]
        result = compute_gene_meta("G", 0, des)
        # Heterogeneous effects should produce high I²
        assert result.cochran_q > 0
        assert result.i_squared > 0


# ---------------------------------------------------------------------------
# 19. Full Phase 6 pipeline
# ---------------------------------------------------------------------------

class TestRunPhase6:
    """Test the integrated Phase 6 pipeline."""

    def test_full_run(self):
        phase1, phase5 = _make_meta_test_data()
        result = run_phase6(phase1, phase5)

        assert isinstance(result, Phase6Result)
        assert result.n_genes_analyzed == 5
        assert result.n_significant_random >= 0

    def test_all_genes_analyzed(self):
        phase1, phase5 = _make_meta_test_data()
        result = run_phase6(phase1, phase5)

        for i in range(5):
            gene = f"META_GENE_{i}"
            assert gene in result.gene_results
            meta = result.gene_results[gene]
            assert meta.n_datasets == 3

    def test_scale_warnings(self):
        phase1, phase5 = _make_meta_test_data()
        # DS1 has max expression < 5.0 (log-ratio scale)
        result = run_phase6(
            phase1, phase5,
            dataset_expression_ranges={"DS1": 2.0, "DS2": 14.0, "DS3": 15.0},
        )
        assert len(result.scale_warnings) >= 1
        assert any("DS1" in w for w in result.scale_warnings)

    def test_direction_consistent(self):
        phase1, phase5 = _make_meta_test_data()
        result = run_phase6(phase1, phase5)

        for gene, meta in result.gene_results.items():
            # All test genes are downregulated
            assert meta.direction == "down"
            assert meta.random_effect < 0
```

---

## EXECUTION INSTRUCTIONS

After writing phase6_meta.py and appending tests, run:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_phases.py -v 2>&1
```

**Expected: ALL tests pass.** If any fail, report FULL output — do NOT modify tests or code without reporting first.

Then confirmation:

```bash
cd /home/kai001/riker-engine && python -c "
from riker.phases.phase6_meta import run_phase6
from riker.phases.phase1_crossref import Phase1Result, GeneResult, GeneDatasetDE
from riker.phases.phase5_replication import Phase5Result, GeneVerdict
import numpy as np

# Build mock data: 8 genes across 4 datasets
np.random.seed(42)
study_genes = {}
for i in range(8):
    gene = f'GENE_{i}'
    des = []
    for j, ds in enumerate(['DS1', 'DS2', 'DS3', 'DS4']):
        fc = -0.5 - i * 0.05 + np.random.normal(0, 0.1)
        des.append(GeneDatasetDE(
            gene=gene, dataset_id=ds,
            log2fc=fc, p_value=0.001 + i * 0.003,
            t_statistic=-3.0 + np.random.normal(0, 0.3),
            df=28.0, se=0.12 + j * 0.01,
            n_cases=15, n_controls=15, direction='down',
        ))
    study_genes[gene] = GeneResult(
        gene=gene, de_results=des,
        n_datasets_detected=4, n_datasets_significant=4,
        passes_filter=True,
        mean_log2fc=float(np.mean([d.log2fc for d in des])),
        consistent_direction=True,
    )

phase1 = Phase1Result(study_genes=study_genes, n_seed_genes=200, n_study_genes=8)

# Phase 5: all survived
gene_verdicts = {}
for i in range(8):
    gene = f'GENE_{i}'
    gene_verdicts[gene] = GeneVerdict(
        gene=gene, cluster_id=0, status='survived',
        reason='OK', replication_results=[],
        n_brain_concordant=2, n_brain_discordant=0,
        n_blood_tested=0, n_blood_concordant=0,
        discovery_direction='down',
    )
phase5 = Phase5Result(gene_verdicts=gene_verdicts, n_survived=8,
                      locked_core_genes=[f'GENE_{i}' for i in range(8)])

# Run with scale check
result = run_phase6(phase1, phase5,
                    dataset_expression_ranges={'DS1': 14.0, 'DS2': 2.0,
                                               'DS3': 15.0, 'DS4': 13.0})

print('=== Phase 6 Meta-Analysis ===')
print(f'Genes analyzed: {result.n_genes_analyzed}')
print(f'Significant (random): {result.n_significant_random}')
print(f'High heterogeneity: {result.n_high_heterogeneity}')
print(f'Scale warnings: {result.scale_warnings}')
print()
for gene, meta in sorted(result.gene_results.items()):
    print(f'  {gene}: RE={meta.random_effect:.3f} (p={meta.random_p:.4f}), '
          f'I²={meta.i_squared:.1f}%, direction={meta.direction}')

assert result.n_genes_analyzed == 8
assert all(m.direction == 'down' for m in result.gene_results.values())
assert len(result.scale_warnings) >= 1  # DS2 has max=2.0
print()
print('PASS: phase6_meta.py working correctly')
"
```

Regression:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_stats.py tests/test_ingestion.py -q 2>&1
```

Report all three outputs back.
