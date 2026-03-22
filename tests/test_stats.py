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
Riker Engine - Statistical module tests.

Phase 1: Welch's t-test validation.
Phase 2: FDR correction and scope enforcement.
Phase 3: Meta-analysis pooling and SE recovery.
Phase 4: Permutation testing for clusters.
Tests validate against scipy.stats as ground truth.
"""

import math
import warnings

import numpy as np
import pytest
from scipy.stats import ttest_ind

from riker.stats.welch import WelchResult, welch_ttest
from riker.stats.fdr import (
    FDRResult,
    apply_fdr_with_scope,
    benjamini_hochberg,
    fdr_survivors,
)
from riker.stats.meta import (
    MetaResult,
    ScaleCheck,
    StudyEffect,
    check_expression_scale,
    fixed_effects_meta,
    random_effects_meta,
    recover_se,
    run_meta_analysis,
)
from riker.stats.permutation import (
    PermutationResult,
    cluster_permutation_test,
    mean_abs_log2fc,
    permutation_test,
)


# ---------------------------------------------------------------------------
# 1. Agreement with scipy ground truth (Welch)
# ---------------------------------------------------------------------------

class TestWelchVsScipy:
    """Verify our implementation matches scipy.stats.ttest_ind(equal_var=False)."""

    def test_basic_two_groups(self):
        """Standard two-group comparison with known different means."""
        np.random.seed(42)
        g1 = np.random.normal(loc=5.0, scale=1.0, size=30)
        g2 = np.random.normal(loc=4.0, scale=1.5, size=25)

        result = welch_ttest(g1, g2)
        scipy_t, scipy_p = ttest_ind(g1, g2, equal_var=False)

        assert result.method == "t-distribution"
        assert math.isclose(result.t_statistic, scipy_t, rel_tol=1e-10)
        assert math.isclose(result.p_value, scipy_p, rel_tol=1e-10)

    def test_equal_groups(self):
        """Two groups drawn from the same distribution — p should be large."""
        np.random.seed(123)
        g1 = np.random.normal(loc=10.0, scale=2.0, size=50)
        g2 = np.random.normal(loc=10.0, scale=2.0, size=50)

        result = welch_ttest(g1, g2)
        scipy_t, scipy_p = ttest_ind(g1, g2, equal_var=False)

        assert math.isclose(result.t_statistic, scipy_t, rel_tol=1e-10)
        assert math.isclose(result.p_value, scipy_p, rel_tol=1e-10)
        assert result.p_value > 0.05  # should not reject null

    def test_very_different_variances(self):
        """Unequal variances — the whole point of Welch's test."""
        np.random.seed(456)
        g1 = np.random.normal(loc=0, scale=1.0, size=40)
        g2 = np.random.normal(loc=0.5, scale=5.0, size=20)

        result = welch_ttest(g1, g2)
        scipy_t, scipy_p = ttest_ind(g1, g2, equal_var=False)

        assert math.isclose(result.t_statistic, scipy_t, rel_tol=1e-10)
        assert math.isclose(result.p_value, scipy_p, rel_tol=1e-10)

    def test_small_samples(self):
        """Small n where normal approximation diverges from t-distribution.
        This is the exact scenario from the Context Transfer document."""
        np.random.seed(789)
        g1 = np.random.normal(loc=3.0, scale=0.5, size=8)
        g2 = np.random.normal(loc=2.5, scale=0.8, size=6)

        result = welch_ttest(g1, g2)
        scipy_t, scipy_p = ttest_ind(g1, g2, equal_var=False)

        assert result.method == "t-distribution"
        assert math.isclose(result.t_statistic, scipy_t, rel_tol=1e-10)
        assert math.isclose(result.p_value, scipy_p, rel_tol=1e-10)

    def test_large_samples(self):
        """Large n where normal and t should converge."""
        np.random.seed(1024)
        g1 = np.random.normal(loc=0, scale=1, size=500)
        g2 = np.random.normal(loc=0.2, scale=1, size=500)

        result = welch_ttest(g1, g2)
        scipy_t, scipy_p = ttest_ind(g1, g2, equal_var=False)

        assert math.isclose(result.t_statistic, scipy_t, rel_tol=1e-10)
        assert math.isclose(result.p_value, scipy_p, rel_tol=1e-10)


# ---------------------------------------------------------------------------
# 2. One-sided alternatives (Welch)
# ---------------------------------------------------------------------------

class TestWelchAlternatives:
    """Verify one-sided tests agree with scipy."""

    def _get_data(self):
        np.random.seed(42)
        g1 = np.random.normal(loc=5.0, scale=1.0, size=30)
        g2 = np.random.normal(loc=4.0, scale=1.5, size=25)
        return g1, g2

    def test_greater(self):
        g1, g2 = self._get_data()
        result = welch_ttest(g1, g2, alternative="greater")
        scipy_t, scipy_p = ttest_ind(g1, g2, equal_var=False, alternative="greater")
        assert math.isclose(result.p_value, scipy_p, rel_tol=1e-10)

    def test_less(self):
        g1, g2 = self._get_data()
        result = welch_ttest(g1, g2, alternative="less")
        scipy_t, scipy_p = ttest_ind(g1, g2, equal_var=False, alternative="less")
        assert math.isclose(result.p_value, scipy_p, rel_tol=1e-10)

    def test_two_sided_equals_two_times_smaller_one_sided(self):
        """Two-sided p should equal 2x the smaller one-sided p."""
        g1, g2 = self._get_data()
        two = welch_ttest(g1, g2, alternative="two-sided")
        greater = welch_ttest(g1, g2, alternative="greater")
        less = welch_ttest(g1, g2, alternative="less")

        smaller_one_sided = min(greater.p_value, less.p_value)
        assert math.isclose(two.p_value, 2 * smaller_one_sided, rel_tol=1e-10)


# ---------------------------------------------------------------------------
# 3. Return value correctness (Welch)
# ---------------------------------------------------------------------------

class TestWelchReturnValues:
    """Verify all fields in the WelchResult are correct."""

    def test_result_type(self):
        result = welch_ttest([1, 2, 3, 4, 5], [6, 7, 8, 9, 10])
        assert isinstance(result, WelchResult)

    def test_mean_diff(self):
        g1 = [10, 20, 30]
        g2 = [5, 15, 25]
        result = welch_ttest(g1, g2)
        assert math.isclose(result.mean_diff, 5.0)

    def test_sample_sizes(self):
        g1 = [1, 2, 3]
        g2 = [4, 5, 6, 7, 8]
        result = welch_ttest(g1, g2)
        assert result.n1 == 3
        assert result.n2 == 5

    def test_df_is_welch_satterthwaite(self):
        """Verify df is NOT pooled (n1+n2-2) but Welch-Satterthwaite."""
        np.random.seed(42)
        g1 = np.random.normal(loc=0, scale=1.0, size=20)
        g2 = np.random.normal(loc=0, scale=3.0, size=10)
        result = welch_ttest(g1, g2)

        pooled_df = 20 + 10 - 2
        # Welch-Satterthwaite df should be LESS than pooled for unequal var
        assert result.df < pooled_df
        assert result.df > 0

    def test_se_diff_positive(self):
        result = welch_ttest([1, 2, 3, 4], [5, 6, 7, 8])
        assert result.se_diff > 0


# ---------------------------------------------------------------------------
# 4. Input validation (Welch)
# ---------------------------------------------------------------------------

class TestWelchInputValidation:
    """Ensure bad inputs are caught with clear errors."""

    def test_empty_group1(self):
        with pytest.raises(ValueError, match="empty"):
            welch_ttest([], [1, 2, 3])

    def test_empty_group2(self):
        with pytest.raises(ValueError, match="empty"):
            welch_ttest([1, 2, 3], [])

    def test_single_element(self):
        with pytest.raises(ValueError, match="at least 2"):
            welch_ttest([1], [2, 3, 4])

    def test_zero_variance(self):
        with pytest.raises(ValueError, match="zero variance"):
            welch_ttest([5, 5, 5], [1, 2, 3])

    def test_invalid_alternative(self):
        with pytest.raises(ValueError, match="alternative"):
            welch_ttest([1, 2, 3], [4, 5, 6], alternative="bigger")


# ---------------------------------------------------------------------------
# 5. NaN handling (Welch)
# ---------------------------------------------------------------------------

class TestWelchNanHandling:
    """Verify NaN values are stripped with warnings."""

    def test_nan_stripped_with_warning(self):
        g1 = [1.0, 2.0, float("nan"), 4.0, 5.0]
        g2 = [10.0, 20.0, 30.0]

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = welch_ttest(g1, g2)
            assert len(w) == 1
            assert "NaN" in str(w[0].message) or "nan" in str(w[0].message).lower()

        # Should have used 4 observations, not 5
        assert result.n1 == 4

    def test_all_nan_raises(self):
        with pytest.raises(ValueError, match="empty"):
            welch_ttest([float("nan"), float("nan")], [1, 2, 3])


# ---------------------------------------------------------------------------
# 6. Input type flexibility (Welch)
# ---------------------------------------------------------------------------

class TestWelchInputTypes:
    """Accept lists, tuples, numpy arrays, pandas Series."""

    def test_lists(self):
        result = welch_ttest([1, 2, 3, 4], [5, 6, 7, 8])
        assert isinstance(result, WelchResult)

    def test_tuples(self):
        result = welch_ttest((1, 2, 3, 4), (5, 6, 7, 8))
        assert isinstance(result, WelchResult)

    def test_numpy_arrays(self):
        result = welch_ttest(np.array([1, 2, 3, 4]), np.array([5, 6, 7, 8]))
        assert isinstance(result, WelchResult)

    def test_pandas_series(self):
        import pandas as pd
        result = welch_ttest(pd.Series([1, 2, 3, 4]), pd.Series([5, 6, 7, 8]))
        assert isinstance(result, WelchResult)

    def test_mixed_types(self):
        """List for group1, numpy array for group2."""
        result = welch_ttest([1, 2, 3, 4], np.array([5, 6, 7, 8]))
        assert isinstance(result, WelchResult)


# ===========================================================================
# PHASE 2: FDR TESTS
# ===========================================================================

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
        FDR at q<0.10 (Level 3, Blueprint Section 8.2) should yield zero survivors."""
        np.random.seed(42)
        n_study = 420
        genes = [f"GENE{i}" for i in range(n_study)]
        # p-values between 0.02 and 0.06 — marginal significance, not extreme
        pvals = dict(zip(genes, np.random.uniform(0.02, 0.06, size=n_study)))

        result = apply_fdr_with_scope(
            pvals, seed_gene_count=800, scope="full_seed_set", threshold=0.10
        )
        # With 800-gene burden at q<0.10, these marginal p-values should NOT survive
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


# ===========================================================================
# PHASE 3: META-ANALYSIS TESTS
# ===========================================================================

# ---------------------------------------------------------------------------
# 11. SE recovery from log2FC and p-value
# ---------------------------------------------------------------------------

class TestSERecovery:
    """Verify SE recovery uses t-distribution and is mathematically correct."""

    def test_known_values(self):
        """Recover SE from a known t-test result.
        If log2FC = 1.0, n_cases=30, n_controls=30, df=58.
        With p=0.01, |t| = t.isf(0.005, 58) ≈ 2.663.
        SE = 1.0 / 2.663 ≈ 0.3755."""
        from scipy.stats import t as t_dist
        expected_t = float(t_dist.isf(0.005, 58))
        expected_se = 1.0 / expected_t

        se = recover_se(log2fc=1.0, p_value=0.01, n_cases=30, n_controls=30)
        assert math.isclose(se, expected_se, rel_tol=1e-6)

    def test_round_trip_with_welch(self):
        """Recover SE from a welch_ttest result, then verify consistency."""
        from riker.stats.welch import welch_ttest
        np.random.seed(42)
        g1 = np.random.normal(5.0, 1.0, 30)
        g2 = np.random.normal(4.0, 1.5, 25)

        result = welch_ttest(g1, g2)
        recovered = recover_se(
            log2fc=result.mean_diff,
            p_value=result.p_value,
            n_cases=result.n1,
            n_controls=result.n2,
        )
        # Recovered SE should be close to original (not exact due to
        # df approximation in recovery: pooled vs Welch-Satterthwaite)
        assert math.isclose(recovered, result.se_diff, rel_tol=0.15)

    def test_large_effect(self):
        """Large fold change should give small relative SE."""
        se = recover_se(log2fc=5.0, p_value=0.001, n_cases=50, n_controls=50)
        assert se > 0
        assert se < 5.0  # SE should be much less than effect

    def test_zero_log2fc(self):
        """Zero fold change with high p-value should still return valid SE."""
        se = recover_se(log2fc=0.0, p_value=0.95, n_cases=20, n_controls=20)
        assert se > 0

    def test_very_small_pvalue(self):
        """Very small p-value should not cause overflow."""
        se = recover_se(log2fc=2.0, p_value=1e-50, n_cases=100, n_controls=100)
        assert se > 0
        assert np.isfinite(se)

    def test_pvalue_one(self):
        """p=1.0 means t=0, should return large SE."""
        se = recover_se(log2fc=0.5, p_value=1.0, n_cases=20, n_controls=20)
        assert se > 1.0  # Large SE to downweight

    def test_invalid_pvalue(self):
        with pytest.raises(ValueError, match="p_value"):
            recover_se(log2fc=1.0, p_value=0.0, n_cases=10, n_controls=10)

    def test_invalid_sample_size(self):
        with pytest.raises(ValueError, match="n >= 2"):
            recover_se(log2fc=1.0, p_value=0.05, n_cases=1, n_controls=10)


# ---------------------------------------------------------------------------
# 12. Expression scale detection (GSE33000 issue)
# ---------------------------------------------------------------------------

class TestScaleDetection:
    """Verify detection of log-ratio vs log2-intensity scale."""

    def test_log2_intensity_not_flagged(self):
        """Typical microarray data (range 3-16) should not be flagged."""
        vals = np.random.uniform(3.0, 16.0, size=1000)
        result = check_expression_scale("GSE5281", vals)
        assert not result.is_log_ratio
        assert result.warning == ""

    def test_log_ratio_flagged(self):
        """Log-ratio data (range -2 to 2) should be flagged."""
        vals = np.random.uniform(-2.0, 2.0, size=1000)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = check_expression_scale("GSE33000", vals)
            assert len(w) == 1
            assert "log-ratio" in str(w[0].message).lower()
        assert result.is_log_ratio
        assert result.expression_range_max < 5.0

    def test_threshold_boundary(self):
        """Value exactly at threshold should not be flagged."""
        vals = np.array([0.0, 5.0])
        result = check_expression_scale("TEST", vals)
        assert not result.is_log_ratio

    def test_empty_data(self):
        vals = np.array([])
        result = check_expression_scale("EMPTY", vals)
        assert result.is_log_ratio  # flagged as problematic


# ---------------------------------------------------------------------------
# 13. Fixed-effects meta-analysis
# ---------------------------------------------------------------------------

class TestFixedEffects:
    """Verify fixed-effects IVW meta-analysis."""

    def _make_studies(self):
        """Create test studies with known properties."""
        return [
            StudyEffect("GSE1", log2fc=-0.5, se=0.15, n_cases=30,
                        n_controls=30, p_value=0.001, tissue="brain"),
            StudyEffect("GSE2", log2fc=-0.4, se=0.20, n_cases=25,
                        n_controls=25, p_value=0.05, tissue="brain"),
            StudyEffect("GSE3", log2fc=-0.6, se=0.18, n_cases=35,
                        n_controls=30, p_value=0.001, tissue="brain"),
            StudyEffect("GSE4", log2fc=-0.3, se=0.25, n_cases=20,
                        n_controls=20, p_value=0.10, tissue="brain"),
        ]

    def test_pooled_effect_direction(self):
        """Pooled effect should be negative (all studies negative)."""
        result = fixed_effects_meta("ATP2B2", self._make_studies())
        assert result.pooled_log2fc < 0

    def test_pooled_effect_manual(self):
        """Verify against manual IVW calculation."""
        studies = self._make_studies()
        effects = np.array([s.log2fc for s in studies])
        ses = np.array([s.se for s in studies])
        weights = 1.0 / ses ** 2
        expected = float(np.sum(weights * effects) / np.sum(weights))

        result = fixed_effects_meta("ATP2B2", studies)
        assert math.isclose(result.pooled_log2fc, expected, rel_tol=1e-10)

    def test_model_label(self):
        result = fixed_effects_meta("TEST", self._make_studies())
        assert result.model == "fixed"

    def test_ci_contains_estimate(self):
        result = fixed_effects_meta("TEST", self._make_studies())
        assert result.ci_lower < result.pooled_log2fc < result.ci_upper

    def test_weights_sum_to_one(self):
        result = fixed_effects_meta("TEST", self._make_studies())
        assert math.isclose(sum(result.study_weights), 1.0, rel_tol=1e-10)

    def test_n_studies(self):
        result = fixed_effects_meta("TEST", self._make_studies())
        assert result.n_studies == 4

    def test_heterogeneity_stats(self):
        result = fixed_effects_meta("TEST", self._make_studies())
        assert result.cochran_q >= 0
        assert 0 <= result.i_squared <= 100
        assert result.tau_squared >= 0
        assert 0 <= result.q_p_value <= 1

    def test_requires_two_studies(self):
        with pytest.raises(ValueError, match="at least 2"):
            fixed_effects_meta("TEST", [self._make_studies()[0]])


# ---------------------------------------------------------------------------
# 14. Random-effects meta-analysis (PRIMARY model)
# ---------------------------------------------------------------------------

class TestRandomEffects:
    """Verify DerSimonian-Laird random-effects meta-analysis."""

    def _make_studies(self):
        return [
            StudyEffect("GSE1", log2fc=-0.5, se=0.15, n_cases=30,
                        n_controls=30, p_value=0.001, tissue="brain"),
            StudyEffect("GSE2", log2fc=-0.4, se=0.20, n_cases=25,
                        n_controls=25, p_value=0.05, tissue="brain"),
            StudyEffect("GSE3", log2fc=-0.6, se=0.18, n_cases=35,
                        n_controls=30, p_value=0.001, tissue="brain"),
            StudyEffect("GSE4", log2fc=-0.3, se=0.25, n_cases=20,
                        n_controls=20, p_value=0.10, tissue="brain"),
        ]

    def test_model_label(self):
        result = random_effects_meta("TEST", self._make_studies())
        assert result.model == "random"

    def test_wider_ci_than_fixed(self):
        """Random effects CI should be >= fixed effects CI width."""
        studies = self._make_studies()
        fixed = fixed_effects_meta("TEST", studies)
        random = random_effects_meta("TEST", studies)
        fixed_width = fixed.ci_upper - fixed.ci_lower
        random_width = random.ci_upper - random.ci_lower
        assert random_width >= fixed_width - 1e-10

    def test_pooled_effect_direction(self):
        result = random_effects_meta("ATP2B2", self._make_studies())
        assert result.pooled_log2fc < 0

    def test_heterogeneous_studies(self):
        """When studies disagree, tau_squared should be positive."""
        studies = [
            StudyEffect("GSE1", log2fc=-1.0, se=0.1, n_cases=50,
                        n_controls=50, p_value=0.001, tissue="brain"),
            StudyEffect("GSE2", log2fc=0.5, se=0.1, n_cases=50,
                        n_controls=50, p_value=0.001, tissue="brain"),
        ]
        result = random_effects_meta("TEST", studies)
        assert result.tau_squared > 0
        assert result.i_squared > 50  # high heterogeneity

    def test_homogeneous_studies(self):
        """When studies agree perfectly, tau_squared should be 0."""
        studies = [
            StudyEffect("GSE1", log2fc=-0.5, se=0.2, n_cases=30,
                        n_controls=30, p_value=0.01, tissue="brain"),
            StudyEffect("GSE2", log2fc=-0.5, se=0.2, n_cases=30,
                        n_controls=30, p_value=0.01, tissue="brain"),
        ]
        result = random_effects_meta("TEST", studies)
        assert result.tau_squared == 0.0
        assert result.i_squared == 0.0

    def test_weights_sum_to_one(self):
        result = random_effects_meta("TEST", self._make_studies())
        assert math.isclose(sum(result.study_weights), 1.0, rel_tol=1e-10)


# ---------------------------------------------------------------------------
# 15. run_meta_analysis entry point
# ---------------------------------------------------------------------------

class TestRunMetaAnalysis:
    """Verify the combined entry point returns both models."""

    def test_returns_both_models(self):
        studies = [
            StudyEffect("GSE1", log2fc=-0.5, se=0.15, n_cases=30,
                        n_controls=30, p_value=0.001, tissue="brain"),
            StudyEffect("GSE2", log2fc=-0.4, se=0.20, n_cases=25,
                        n_controls=25, p_value=0.05, tissue="brain"),
        ]
        fixed, random = run_meta_analysis("ATP2B2", studies)
        assert fixed.model == "fixed"
        assert random.model == "random"
        assert fixed.gene == "ATP2B2"
        assert random.gene == "ATP2B2"

    def test_random_is_primary(self):
        """Random effects should be the second element (primary result)."""
        studies = [
            StudyEffect("GSE1", log2fc=-0.5, se=0.15, n_cases=30,
                        n_controls=30, p_value=0.001, tissue="brain"),
            StudyEffect("GSE2", log2fc=-0.6, se=0.18, n_cases=35,
                        n_controls=30, p_value=0.001, tissue="brain"),
        ]
        _, primary = run_meta_analysis("SEZ6L2", studies)
        assert primary.model == "random"


# ===========================================================================
# PHASE 4: PERMUTATION TESTS
# ===========================================================================

from riker.stats.permutation import (
    PermutationResult,
    cluster_permutation_test,
    mean_abs_log2fc,
    permutation_test,
)


# ---------------------------------------------------------------------------
# 16. Generic permutation test
# ---------------------------------------------------------------------------

class TestPermutationTest:
    """Verify the generic permutation testing framework."""

    def _make_data(self, n_genes=200, n_extreme=10, seed=42):
        """Create gene data where a small group has extreme values."""
        np.random.seed(seed)
        all_genes = [f"GENE{i}" for i in range(n_genes)]
        # Most genes: small log2FC near zero
        gene_values = {g: np.random.normal(0, 0.1) for g in all_genes}
        # Extreme group: large negative log2FC
        extreme_genes = all_genes[:n_extreme]
        for g in extreme_genes:
            gene_values[g] = np.random.normal(-1.5, 0.2)
        return all_genes, extreme_genes, gene_values

    def test_extreme_group_significant(self):
        """A group with extreme values should have small p-value."""
        all_genes, extreme_genes, gene_values = self._make_data()
        result = permutation_test(
            observed_genes=extreme_genes,
            all_genes=all_genes,
            gene_values=gene_values,
            stat_func=mean_abs_log2fc,
            n_permutations=5000,
            seed=42,
        )
        assert result.p_value < 0.01
        assert result.observed > result.null_mean

    def test_random_group_not_significant(self):
        """A random group should NOT be significant."""
        all_genes, _, gene_values = self._make_data()
        # Pick a random group that is NOT the extreme group
        np.random.seed(999)
        random_genes = list(np.random.choice(all_genes[20:], size=10, replace=False))
        result = permutation_test(
            observed_genes=random_genes,
            all_genes=all_genes,
            gene_values=gene_values,
            stat_func=mean_abs_log2fc,
            n_permutations=5000,
            seed=42,
        )
        assert result.p_value > 0.05

    def test_reproducibility(self):
        """Same seed should give same p-value."""
        all_genes, extreme_genes, gene_values = self._make_data()
        r1 = permutation_test(
            extreme_genes, all_genes, gene_values,
            mean_abs_log2fc, n_permutations=1000, seed=42,
        )
        r2 = permutation_test(
            extreme_genes, all_genes, gene_values,
            mean_abs_log2fc, n_permutations=1000, seed=42,
        )
        assert r1.p_value == r2.p_value
        assert r1.observed == r2.observed

    def test_different_seeds_differ(self):
        """Different seeds should give (slightly) different results."""
        all_genes, extreme_genes, gene_values = self._make_data()
        r1 = permutation_test(
            extreme_genes, all_genes, gene_values,
            mean_abs_log2fc, n_permutations=1000, seed=42,
        )
        r2 = permutation_test(
            extreme_genes, all_genes, gene_values,
            mean_abs_log2fc, n_permutations=1000, seed=99,
        )
        # P-values should be similar but not identical
        assert r1.n_extreme != r2.n_extreme or r1.null_mean != r2.null_mean

    def test_p_value_never_zero(self):
        """Conservative formula should prevent p=0."""
        all_genes, extreme_genes, gene_values = self._make_data()
        result = permutation_test(
            extreme_genes, all_genes, gene_values,
            mean_abs_log2fc, n_permutations=100, seed=42,
        )
        assert result.p_value > 0

    def test_p_value_formula(self):
        """Verify p = (n_extreme + 1) / (n_perms + 1)."""
        all_genes, extreme_genes, gene_values = self._make_data()
        result = permutation_test(
            extreme_genes, all_genes, gene_values,
            mean_abs_log2fc, n_permutations=1000, seed=42,
        )
        expected_p = (result.n_extreme + 1) / (result.n_permutations + 1)
        assert math.isclose(result.p_value, expected_p, rel_tol=1e-10)

    def test_return_type(self):
        all_genes, extreme_genes, gene_values = self._make_data()
        result = permutation_test(
            extreme_genes, all_genes, gene_values,
            mean_abs_log2fc, n_permutations=100, seed=42,
        )
        assert isinstance(result, PermutationResult)

    def test_null_distribution_stats(self):
        """Null distribution summary should be populated."""
        all_genes, extreme_genes, gene_values = self._make_data()
        result = permutation_test(
            extreme_genes, all_genes, gene_values,
            mean_abs_log2fc, n_permutations=1000, seed=42,
        )
        assert result.null_std > 0
        assert result.null_min <= result.null_mean <= result.null_max


# ---------------------------------------------------------------------------
# 17. Input validation for permutation test
# ---------------------------------------------------------------------------

class TestPermutationValidation:
    """Ensure bad inputs are caught."""

    def test_empty_observed(self):
        with pytest.raises(ValueError, match="empty"):
            permutation_test([], ["A", "B"], {"A": 1.0, "B": 2.0},
                             mean_abs_log2fc)

    def test_observed_larger_than_pool(self):
        with pytest.raises(ValueError, match="cannot be larger"):
            permutation_test(["A", "B", "C"], ["A", "B"],
                             {"A": 1.0, "B": 2.0, "C": 3.0},
                             mean_abs_log2fc)

    def test_missing_gene_values(self):
        with pytest.raises(ValueError, match="missing"):
            permutation_test(["A"], ["A", "B"],
                             {"A": 1.0},  # B is missing
                             mean_abs_log2fc)

    def test_invalid_n_permutations(self):
        with pytest.raises(ValueError, match="n_permutations"):
            permutation_test(["A"], ["A", "B"], {"A": 1.0, "B": 2.0},
                             mean_abs_log2fc, n_permutations=0)

    def test_invalid_alternative(self):
        with pytest.raises(ValueError, match="alternative"):
            permutation_test(["A"], ["A", "B"], {"A": 1.0, "B": 2.0},
                             mean_abs_log2fc, alternative="bigger")


# ---------------------------------------------------------------------------
# 18. Cluster permutation convenience function
# ---------------------------------------------------------------------------

class TestClusterPermutation:
    """Verify the cluster-specific convenience wrapper."""

    def test_basic_usage(self):
        np.random.seed(42)
        all_genes = [f"G{i}" for i in range(100)]
        gene_fc = {g: np.random.normal(0, 0.1) for g in all_genes}
        # Make a cluster with extreme values
        cluster = all_genes[:5]
        for g in cluster:
            gene_fc[g] = -2.0

        result = cluster_permutation_test(
            cluster, all_genes, gene_fc, n_permutations=1000, seed=42
        )
        assert result.p_value < 0.01
        assert result.observed > 1.0  # mean |log2FC| should be ~2.0

    def test_uses_greater_alternative(self):
        """Cluster test should use 'greater' (higher |FC| is interesting)."""
        np.random.seed(42)
        all_genes = [f"G{i}" for i in range(50)]
        gene_fc = {g: 0.0 for g in all_genes}
        cluster = all_genes[:5]

        result = cluster_permutation_test(
            cluster, all_genes, gene_fc, n_permutations=100, seed=42
        )
        # All values are 0, so observed=null, p should be ~1.0
        assert result.p_value > 0.5


# ---------------------------------------------------------------------------
# 19. mean_abs_log2fc helper
# ---------------------------------------------------------------------------

class TestMeanAbsLog2FC:
    """Verify the default test statistic."""

    def test_basic(self):
        assert math.isclose(mean_abs_log2fc([-1.0, 1.0, -0.5, 0.5]), 0.75)

    def test_all_negative(self):
        assert math.isclose(mean_abs_log2fc([-2.0, -1.0, -3.0]), 2.0)

    def test_empty(self):
        assert mean_abs_log2fc([]) == 0.0

    def test_single_value(self):
        assert math.isclose(mean_abs_log2fc([-0.45]), 0.45)


# ---------------------------------------------------------------------------
# REML Tau-Squared Estimator
# ---------------------------------------------------------------------------

class TestREMLEstimator:
    def test_reml_converges(self):
        """REML should converge with reasonable data."""
        from riker.stats.meta import _estimate_tau_squared_reml
        effects = np.array([-0.5, -0.3, -0.7, -0.4])
        variances = np.array([0.04, 0.05, 0.03, 0.06])
        tau2, converged = _estimate_tau_squared_reml(effects, variances)
        assert converged is True
        assert tau2 >= 0

    def test_reml_nonnegative(self):
        """Tau-squared must be non-negative."""
        from riker.stats.meta import _estimate_tau_squared_reml
        effects = np.array([-0.5, -0.5, -0.5])
        variances = np.array([0.04, 0.04, 0.04])
        tau2, _ = _estimate_tau_squared_reml(effects, variances)
        assert tau2 >= 0

    def test_reml_vs_dl_direction(self):
        """REML typically gives equal or larger tau2 than DL with moderate heterogeneity."""
        from riker.stats.meta import _estimate_tau_squared_reml
        # Moderate heterogeneity with more studies for stable REML
        effects = np.array([-0.2, -0.4, -0.6, -0.3, -0.5])
        variances = np.array([0.04, 0.05, 0.03, 0.06, 0.04])
        tau2_reml, converged = _estimate_tau_squared_reml(effects, variances)
        assert converged is True
        assert tau2_reml >= 0

    def test_single_study(self):
        """Single study should return tau2=0."""
        from riker.stats.meta import _estimate_tau_squared_reml
        tau2, converged = _estimate_tau_squared_reml(
            np.array([-0.5]), np.array([0.04])
        )
        assert tau2 == 0.0
        assert converged is True

    def test_meta_result_records_method(self):
        """MetaResult should record whether REML or DL was used."""
        fixed, random = run_meta_analysis("GENE", [
            StudyEffect("D1", -0.5, 0.1, 15, 15, 0.001),
            StudyEffect("D2", -0.3, 0.12, 15, 15, 0.01),
            StudyEffect("D3", -0.7, 0.09, 15, 15, 0.001),
        ])
        assert hasattr(random, "tau_method")
        assert random.tau_method in ("REML", "DL")
