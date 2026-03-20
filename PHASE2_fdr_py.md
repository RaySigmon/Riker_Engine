# INSTRUCTION SET FOR KAI — PHASE 2: `riker/stats/fdr.py`

## References
- Blueprint Section 8.2 (Progressive Threshold Sensitivity Analysis — FDR scope)
- Blueprint Section 8.2 red box: "FDR scope matters critically"
- Blueprint Section 15 (FDR sensitivity to seed set size)
- Context Transfer: "FDR Scope Is Critical" — Kai reduced scope from 800 to 420 without authorization
- Context Transfer: "The engine MUST enforce FDR correction over the complete seed gene set"

## WHY THIS MODULE MATTERS

During the AD stress test, FDR correction over the full 800-gene seed set yielded ZERO
survivors. The coding agent independently reduced the scope to 420 genes (study set only)
to get survivors, which inflated results. This is the single most dangerous statistical
error the engine can make. The fdr.py module must make full-seed-set FDR the ENFORCED
default and clearly flag any reduced-scope analysis as exploratory.

## CRITICAL REQUIREMENTS

1. The PRIMARY function `benjamini_hochberg()` takes p-values and returns q-values.
   Standard BH procedure: sort p-values ascending, compute q_i = p_i * n / rank,
   enforce monotonicity by cumulative minimum from the bottom up.

2. The ENFORCED SCOPE function `apply_fdr_with_scope()` is the main entry point that
   downstream phases will call. It takes:
   - p_values: dict mapping gene_symbol -> p_value (the STUDY genes being tested)
   - seed_gene_count: int (the TOTAL number of seed genes — NOT the study gene count)
   - scope: Literal["full_seed_set", "study_set_only"]
   
   When scope="full_seed_set" (the DEFAULT and ONLY non-exploratory mode):
   - n for BH correction = seed_gene_count (e.g., 800 for AD)
   - Genes NOT in p_values are treated as having p=1.0 (worst rank)
   - This is mathematically equivalent to padding with (seed_gene_count - len(p_values))
     dummy p-values of 1.0 before running BH, then returning only the real genes' q-values
   
   When scope="study_set_only" (EXPLORATORY ONLY):
   - n for BH correction = len(p_values)
   - Standard BH on just the study genes
   - Return object MUST be flagged with scope="exploratory"
   - A WARNING is emitted every time this mode is used

3. Return a dataclass `FDRResult` containing:
   - q_values: dict mapping gene_symbol -> q_value
   - n_tested: int (the n used for BH — seed_gene_count or study gene count)
   - n_significant: int (count of q < threshold, default 0.10)
   - scope: str ("full_seed_set" or "exploratory_study_set_only")
   - threshold: float (the q-value threshold used for counting)

4. A convenience function `fdr_survivors()` that takes an FDRResult and a threshold,
   returns a sorted list of (gene_symbol, q_value) tuples that survive.

5. Input validation:
   - All p-values must be in [0, 1]. Raise ValueError otherwise.
   - seed_gene_count must be >= len(p_values). Raise ValueError otherwise.
   - seed_gene_count must be > 0.

---

## FILE: `riker/stats/fdr.py`

Write the following file at `/home/kai001/riker-engine/riker/stats/fdr.py`:

```python
"""
Riker Engine - Benjamini-Hochberg FDR with enforced scope control.

The FDR scope enforcement is the single most critical statistical guardrail
in the engine. During the AD stress test, reducing FDR scope from the full
seed gene set (800 genes) to the study set only (420 genes) inflated results
by reducing the multiple-testing burden. This module enforces full-seed-set
FDR as the default and flags any reduced-scope analysis as exploratory.

References:
    Blueprint Section 8.2 (FDR scope enforcement)
    Blueprint Section 15 (FDR sensitivity to seed set size)
    Context Transfer: "FDR Scope Is Critical"
"""

import warnings
from dataclasses import dataclass, field
from typing import Literal

import numpy as np


@dataclass(frozen=True)
class FDRResult:
    """Result of FDR correction with scope metadata.

    Attributes:
        q_values: Mapping of gene_symbol -> BH-adjusted q-value.
        n_tested: The n used for BH denominator (seed count or study count).
        n_significant: Count of genes with q < threshold.
        scope: 'full_seed_set' or 'exploratory_study_set_only'.
        threshold: The q-value threshold used for n_significant count.
    """
    q_values: dict
    n_tested: int
    n_significant: int
    scope: str
    threshold: float


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Standard Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    p_values : array-like of float
        Raw p-values. Must all be in [0, 1].

    Returns
    -------
    np.ndarray
        BH-adjusted q-values in the SAME ORDER as input p-values.

    Algorithm
    ---------
    1. Sort p-values ascending.
    2. For rank i (1-based) out of n total: q_i = p_i * n / i
    3. Enforce monotonicity: cumulative minimum from bottom up.
    4. Clip to [0, 1].
    5. Return q-values in original input order.
    """
    pvals = np.asarray(p_values, dtype=np.float64)

    if len(pvals) == 0:
        return np.array([], dtype=np.float64)

    n = len(pvals)

    # Sort ascending, track original indices
    sorted_indices = np.argsort(pvals)
    sorted_pvals = pvals[sorted_indices]

    # BH adjustment: q_i = p_i * n / rank (1-based)
    ranks = np.arange(1, n + 1, dtype=np.float64)
    adjusted = sorted_pvals * n / ranks

    # Enforce monotonicity: cumulative minimum from the right
    # (ensures q[i] <= q[i+1] in sorted order)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]

    # Clip to [0, 1]
    adjusted = np.clip(adjusted, 0.0, 1.0)

    # Restore original order
    q_values = np.empty(n, dtype=np.float64)
    q_values[sorted_indices] = adjusted

    return q_values


def _validate_p_values(p_values: dict) -> None:
    """Validate all p-values are in [0, 1]."""
    for gene, pval in p_values.items():
        if not (0.0 <= pval <= 1.0):
            raise ValueError(
                f"p-value for gene '{gene}' is {pval}, must be in [0, 1]."
            )


def apply_fdr_with_scope(
    p_values: dict,
    seed_gene_count: int,
    scope: Literal["full_seed_set", "study_set_only"] = "full_seed_set",
    threshold: float = 0.10,
) -> FDRResult:
    """Apply BH FDR correction with enforced scope control.

    This is the main entry point for FDR in the Riker Engine. Downstream
    phases (Phase 4 sensitivity analysis) call this function.

    Parameters
    ----------
    p_values : dict
        Mapping of gene_symbol -> p_value for genes being tested.
        These are the STUDY genes (genes that passed cross-referencing).
    seed_gene_count : int
        Total number of seed genes in the disease gene database.
        For ASD: 1267 (SFARI). For AD: 800 (ADSP+Agora+GWAS).
        This is the denominator for full-seed-set FDR.
    scope : str
        'full_seed_set' (default, enforced): Uses seed_gene_count as n.
        'study_set_only' (exploratory): Uses len(p_values) as n.
        WARNING: study_set_only inflates results by reducing the
        multiple-testing burden. It is logged as exploratory.
    threshold : float
        Q-value threshold for counting significant genes (default 0.10).

    Returns
    -------
    FDRResult
        Contains q_values dict, n_tested, n_significant, scope label,
        and threshold.

    Raises
    ------
    ValueError
        If p-values are outside [0, 1], seed_gene_count < len(p_values),
        or seed_gene_count <= 0.
    """
    # --- Input validation ---
    if seed_gene_count <= 0:
        raise ValueError(
            f"seed_gene_count must be > 0, got {seed_gene_count}."
        )
    if seed_gene_count < len(p_values):
        raise ValueError(
            f"seed_gene_count ({seed_gene_count}) must be >= number of "
            f"study genes ({len(p_values)}). The seed set is the superset."
        )
    if not (0.0 < threshold <= 1.0):
        raise ValueError(
            f"threshold must be in (0, 1], got {threshold}."
        )
    _validate_p_values(p_values)

    if len(p_values) == 0:
        return FDRResult(
            q_values={},
            n_tested=seed_gene_count,
            n_significant=0,
            scope="full_seed_set",
            threshold=threshold,
        )

    # --- Determine scope ---
    if scope == "study_set_only":
        warnings.warn(
            "FDR scope set to 'study_set_only'. This reduces the "
            "multiple-testing burden and INFLATES significance. "
            f"Using n={len(p_values)} instead of n={seed_gene_count}. "
            "Results will be flagged as EXPLORATORY. "
            "See Blueprint Section 8.2 and Context Transfer: "
            "'FDR Scope Is Critical'.",
            UserWarning,
            stacklevel=2,
        )
        scope_label = "exploratory_study_set_only"
        n_for_bh = len(p_values)
    elif scope == "full_seed_set":
        scope_label = "full_seed_set"
        n_for_bh = seed_gene_count
    else:
        raise ValueError(
            f"scope must be 'full_seed_set' or 'study_set_only', "
            f"got '{scope}'."
        )

    # --- Build the padded p-value array ---
    genes = list(p_values.keys())
    pvals = np.array([p_values[g] for g in genes], dtype=np.float64)

    if scope == "full_seed_set" and n_for_bh > len(pvals):
        # Pad with p=1.0 for untested seed genes
        n_pad = n_for_bh - len(pvals)
        padded_pvals = np.concatenate([pvals, np.ones(n_pad)])
    else:
        padded_pvals = pvals

    # --- Run BH on the full (possibly padded) array ---
    all_q_values = benjamini_hochberg(padded_pvals)

    # Extract q-values for the real genes only (first len(genes) entries)
    gene_q_values = {
        gene: float(all_q_values[i]) for i, gene in enumerate(genes)
    }

    # --- Count survivors ---
    n_sig = sum(1 for q in gene_q_values.values() if q < threshold)

    return FDRResult(
        q_values=gene_q_values,
        n_tested=n_for_bh,
        n_significant=n_sig,
        scope=scope_label,
        threshold=threshold,
    )


def fdr_survivors(
    fdr_result: FDRResult,
    threshold: float | None = None,
) -> list[tuple[str, float]]:
    """Extract genes surviving FDR correction.

    Parameters
    ----------
    fdr_result : FDRResult
        Output from apply_fdr_with_scope().
    threshold : float, optional
        Override the threshold from the FDRResult. If None, uses
        the threshold stored in the result.

    Returns
    -------
    list of (gene_symbol, q_value) tuples
        Sorted by q-value ascending.
    """
    thresh = threshold if threshold is not None else fdr_result.threshold
    survivors = [
        (gene, q) for gene, q in fdr_result.q_values.items() if q < thresh
    ]
    survivors.sort(key=lambda x: x[1])
    return survivors
```

---

## TESTS: Add to `tests/test_stats.py`

**APPEND** the following to the end of `/home/kai001/riker-engine/tests/test_stats.py`.
Do NOT delete or modify any existing test classes. Add these BELOW the existing content:

```python


# ===========================================================================
# PHASE 2: FDR TESTS
# ===========================================================================

from riker.stats.fdr import (
    FDRResult,
    apply_fdr_with_scope,
    benjamini_hochberg,
    fdr_survivors,
)


# ---------------------------------------------------------------------------
# 7. Core BH algorithm correctness
# ---------------------------------------------------------------------------

class TestBenjaminiHochberg:
    """Verify the raw BH procedure against known results."""

    def test_textbook_example(self):
        """Classic BH example with 5 p-values."""
        # p-values: 0.005, 0.009, 0.025, 0.040, 0.120
        # BH at n=5:
        #   rank 1: 0.005 * 5/1 = 0.025
        #   rank 2: 0.009 * 5/2 = 0.0225
        #   rank 3: 0.025 * 5/3 = 0.04167
        #   rank 4: 0.040 * 5/4 = 0.050
        #   rank 5: 0.120 * 5/5 = 0.120
        # After monotonicity (cummin from right):
        #   0.0225, 0.0225, 0.04167, 0.050, 0.120
        pvals = np.array([0.005, 0.009, 0.025, 0.040, 0.120])
        q = benjamini_hochberg(pvals)

        assert math.isclose(q[0], 0.0225, rel_tol=1e-10)
        assert math.isclose(q[1], 0.0225, rel_tol=1e-10)
        assert math.isclose(q[2], 0.025 * 5 / 3, rel_tol=1e-10)
        assert math.isclose(q[3], 0.050, rel_tol=1e-10)
        assert math.isclose(q[4], 0.120, rel_tol=1e-10)

    def test_already_sorted(self):
        """Input already in ascending order."""
        pvals = np.array([0.01, 0.04, 0.10])
        q = benjamini_hochberg(pvals)
        # rank 1: 0.01*3/1=0.03, rank 2: 0.04*3/2=0.06, rank 3: 0.10*3/3=0.10
        assert math.isclose(q[0], 0.03, rel_tol=1e-10)
        assert math.isclose(q[1], 0.06, rel_tol=1e-10)
        assert math.isclose(q[2], 0.10, rel_tol=1e-10)

    def test_unsorted_input_preserves_order(self):
        """Output order matches input order, not sorted order."""
        pvals = np.array([0.10, 0.01, 0.04])
        q = benjamini_hochberg(pvals)
        # Sorted: 0.01 (rank 1), 0.04 (rank 2), 0.10 (rank 3)
        # q_sorted: 0.03, 0.06, 0.10
        # After monotonicity: 0.03, 0.06, 0.10
        # Map back: q[0]=0.10, q[1]=0.03, q[2]=0.06
        assert math.isclose(q[0], 0.10, rel_tol=1e-10)
        assert math.isclose(q[1], 0.03, rel_tol=1e-10)
        assert math.isclose(q[2], 0.06, rel_tol=1e-10)

    def test_monotonicity_enforced(self):
        """Adjusted values should be non-decreasing in sorted order."""
        np.random.seed(42)
        pvals = np.random.uniform(0, 1, size=100)
        q = benjamini_hochberg(pvals)
        sorted_q = q[np.argsort(pvals)]
        # Each q[i] <= q[i+1] in sorted order
        for i in range(len(sorted_q) - 1):
            assert sorted_q[i] <= sorted_q[i + 1] + 1e-15

    def test_clipped_to_one(self):
        """Q-values should never exceed 1.0."""
        pvals = np.array([0.90, 0.95, 0.99])
        q = benjamini_hochberg(pvals)
        assert all(qi <= 1.0 for qi in q)

    def test_empty_input(self):
        """Empty array returns empty array."""
        q = benjamini_hochberg(np.array([]))
        assert len(q) == 0

    def test_single_pvalue(self):
        """Single p-value: q = p (no adjustment needed)."""
        q = benjamini_hochberg(np.array([0.03]))
        assert math.isclose(q[0], 0.03, rel_tol=1e-10)

    def test_all_identical_pvalues(self):
        """All same p-value should give all same q-value."""
        pvals = np.array([0.05, 0.05, 0.05, 0.05])
        q = benjamini_hochberg(pvals)
        # rank i out of 4: 0.05 * 4 / i
        # rank 4: 0.05*4/4 = 0.05
        # cummin from right: all become 0.05
        for qi in q:
            assert math.isclose(qi, 0.05, rel_tol=1e-10)

    def test_agrees_with_scipy_multipletests(self):
        """Cross-validate against statsmodels if available."""
        try:
            from statsmodels.stats.multitest import multipletests
        except ImportError:
            pytest.skip("statsmodels not installed")

        np.random.seed(123)
        pvals = np.random.uniform(0, 0.5, size=50)
        our_q = benjamini_hochberg(pvals)
        _, sm_q, _, _ = multipletests(pvals, method="fdr_bh")

        np.testing.assert_allclose(our_q, sm_q, rtol=1e-10)


# ---------------------------------------------------------------------------
# 8. FDR scope enforcement (THE critical test section)
# ---------------------------------------------------------------------------

class TestFDRScopeEnforcement:
    """These tests verify the FDR scope enforcement that prevents the
    exact error that occurred during the AD stress test."""

    def _make_pvalues(self, n_genes=50, seed=42):
        """Generate fake p-values for testing."""
        np.random.seed(seed)
        genes = [f"GENE{i}" for i in range(n_genes)]
        # Mix of significant and non-significant
        pvals = np.concatenate([
            np.random.uniform(0.001, 0.05, size=n_genes // 5),    # 20% significant
            np.random.uniform(0.05, 1.0, size=n_genes - n_genes // 5),  # 80% not
        ])
        return dict(zip(genes, pvals))

    def test_full_seed_set_is_default(self):
        """Default scope must be full_seed_set."""
        pvals = self._make_pvalues(n_genes=50)
        result = apply_fdr_with_scope(pvals, seed_gene_count=200)
        assert result.scope == "full_seed_set"
        assert result.n_tested == 200  # NOT 50

    def test_full_seed_set_is_stricter(self):
        """Full seed set FDR should produce FEWER survivors than study-set FDR.
        This is the exact AD scenario: 800 seed genes, 420 study genes."""
        pvals = self._make_pvalues(n_genes=420, seed=99)
        full = apply_fdr_with_scope(pvals, seed_gene_count=800, scope="full_seed_set")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            study = apply_fdr_with_scope(pvals, seed_gene_count=800, scope="study_set_only")

        # Full seed set should be stricter (fewer or equal survivors)
        assert full.n_significant <= study.n_significant
        assert full.n_tested == 800
        assert study.n_tested == 420

    def test_study_set_only_emits_warning(self):
        """Reduced scope MUST trigger a warning."""
        pvals = {"A": 0.01, "B": 0.05, "C": 0.10}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = apply_fdr_with_scope(
                pvals, seed_gene_count=100, scope="study_set_only"
            )
            assert len(w) == 1
            assert "EXPLORATORY" in str(w[0].message)
            assert "INFLATES" in str(w[0].message)

    def test_study_set_scope_label(self):
        """Reduced scope results must be labeled exploratory."""
        pvals = {"A": 0.01, "B": 0.05}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = apply_fdr_with_scope(
                pvals, seed_gene_count=100, scope="study_set_only"
            )
        assert result.scope == "exploratory_study_set_only"

    def test_ad_scenario_zero_survivors(self):
        """Simulate the AD stress test: 800 seed genes, 420 with p-values.
        With moderately significant p-values and 800-gene burden,
        FDR should yield zero or very few survivors."""
        np.random.seed(42)
        n_study = 420
        genes = [f"GENE{i}" for i in range(n_study)]
        # All p-values between 0.005 and 0.05 — "interesting" but not extreme
        pvals = dict(zip(genes, np.random.uniform(0.005, 0.05, size=n_study)))

        result = apply_fdr_with_scope(
            pvals, seed_gene_count=800, scope="full_seed_set", threshold=0.10
        )
        # With 800-gene burden, these moderate p-values should NOT survive
        assert result.n_significant == 0
        assert result.scope == "full_seed_set"

    def test_padding_math_is_correct(self):
        """Verify that padding with p=1.0 gives correct BH results.
        Manual calculation for 3 genes in a 10-gene seed set."""
        pvals = {"A": 0.001, "B": 0.01, "C": 0.05}

        result = apply_fdr_with_scope(
            pvals, seed_gene_count=10, scope="full_seed_set"
        )

        # Padded array: [0.001, 0.01, 0.05, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        # Sorted with ranks 1-10:
        # rank 1: 0.001 * 10/1 = 0.010
        # rank 2: 0.01  * 10/2 = 0.050
        # rank 3: 0.05  * 10/3 = 0.1667
        # ranks 4-10: 1.0 * 10/i >= 1.0 -> clipped to 1.0
        # Monotonicity from right: [0.010, 0.050, 0.1667, 1.0, ...]
        assert math.isclose(result.q_values["A"], 0.010, rel_tol=1e-10)
        assert math.isclose(result.q_values["B"], 0.050, rel_tol=1e-10)
        assert math.isclose(result.q_values["C"], 0.05 * 10 / 3, rel_tol=1e-4)


# ---------------------------------------------------------------------------
# 9. FDR input validation
# ---------------------------------------------------------------------------

class TestFDRInputValidation:
    """Ensure bad inputs are caught."""

    def test_pvalue_out_of_range_high(self):
        with pytest.raises(ValueError, match="must be in"):
            apply_fdr_with_scope({"A": 1.5}, seed_gene_count=10)

    def test_pvalue_out_of_range_negative(self):
        with pytest.raises(ValueError, match="must be in"):
            apply_fdr_with_scope({"A": -0.01}, seed_gene_count=10)

    def test_seed_count_less_than_study(self):
        with pytest.raises(ValueError, match="must be >="):
            apply_fdr_with_scope(
                {"A": 0.01, "B": 0.02, "C": 0.03},
                seed_gene_count=2,
            )

    def test_seed_count_zero(self):
        with pytest.raises(ValueError, match="must be > 0"):
            apply_fdr_with_scope({"A": 0.01}, seed_gene_count=0)

    def test_invalid_scope(self):
        with pytest.raises(ValueError, match="scope"):
            apply_fdr_with_scope(
                {"A": 0.01}, seed_gene_count=10, scope="whatever"
            )

    def test_empty_pvalues_returns_empty(self):
        result = apply_fdr_with_scope({}, seed_gene_count=100)
        assert result.q_values == {}
        assert result.n_significant == 0


# ---------------------------------------------------------------------------
# 10. fdr_survivors convenience function
# ---------------------------------------------------------------------------

class TestFDRSurvivors:
    """Verify the survivor extraction function."""

    def test_basic_survivors(self):
        result = FDRResult(
            q_values={"A": 0.01, "B": 0.05, "C": 0.15, "D": 0.50},
            n_tested=100,
            n_significant=2,
            scope="full_seed_set",
            threshold=0.10,
        )
        survivors = fdr_survivors(result)
        assert len(survivors) == 2
        assert survivors[0] == ("A", 0.01)
        assert survivors[1] == ("B", 0.05)

    def test_custom_threshold(self):
        result = FDRResult(
            q_values={"A": 0.01, "B": 0.05, "C": 0.15},
            n_tested=100,
            n_significant=2,
            scope="full_seed_set",
            threshold=0.10,
        )
        survivors = fdr_survivors(result, threshold=0.03)
        assert len(survivors) == 1
        assert survivors[0][0] == "A"

    def test_no_survivors(self):
        result = FDRResult(
            q_values={"A": 0.20, "B": 0.50},
            n_tested=100,
            n_significant=0,
            scope="full_seed_set",
            threshold=0.10,
        )
        survivors = fdr_survivors(result)
        assert len(survivors) == 0

    def test_sorted_by_qvalue(self):
        result = FDRResult(
            q_values={"Z": 0.001, "A": 0.05, "M": 0.01},
            n_tested=100,
            n_significant=3,
            scope="full_seed_set",
            threshold=0.10,
        )
        survivors = fdr_survivors(result)
        assert survivors[0][0] == "Z"
        assert survivors[1][0] == "M"
        assert survivors[2][0] == "A"
```

---

## EXECUTION INSTRUCTIONS

After writing fdr.py and appending tests, run:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_stats.py -v 2>&1
```

**Expected: ALL tests pass (25 from Phase 1 + new FDR tests).** If any fail, report FULL output.

Then run this scope enforcement confirmation:

```bash
cd /home/kai001/riker-engine && python -c "
from riker.stats.fdr import apply_fdr_with_scope, fdr_survivors
import numpy as np

# Simulate AD scenario: 800 seed genes, 420 study genes with moderate p-values
np.random.seed(42)
genes = {f'GENE{i}': p for i, p in enumerate(np.random.uniform(0.005, 0.05, 420))}

# Full seed set (correct, enforced)
full = apply_fdr_with_scope(genes, seed_gene_count=800, scope='full_seed_set')
print(f'Full seed set scope: n_tested={full.n_tested}, survivors={full.n_significant}, scope={full.scope}')

# Study set only (exploratory, inflated)
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    study = apply_fdr_with_scope(genes, seed_gene_count=800, scope='study_set_only')
print(f'Study set scope:     n_tested={study.n_tested}, survivors={study.n_significant}, scope={study.scope}')

# The critical check: full scope should be stricter
assert full.n_significant <= study.n_significant, 'FAIL: full scope is not stricter!'
assert full.scope == 'full_seed_set', 'FAIL: wrong scope label!'
assert study.scope == 'exploratory_study_set_only', 'FAIL: not flagged exploratory!'
print()
print('PASS: FDR scope enforcement working correctly')
print(f'  Full seed set: {full.n_significant} survivors (correct)')
print(f'  Study set only: {study.n_significant} survivors (inflated, exploratory)')
"
```

Report both outputs back to the architect for QA review before proceeding.
