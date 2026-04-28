# Riker Engine - Condition-Agnostic Transcriptomics Pipeline
# Copyright (C) 2024-2026 Ray Sigmon
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""Tests for the stability profiling pairwise Jaccard computation."""

import sys
from pathlib import Path

import pytest

# Add scripts/ to path so we can import the profiler module
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
from stability_profiling import compute_pairwise_jaccard


class TestPairwiseJaccard:
    """Verify pairwise Jaccard similarity with hand-computable cases."""

    def test_identical_sets(self):
        """All runs produce the same genes -> Jaccard = 1.0 for all pairs."""
        sets = [{"A", "B", "C"}, {"A", "B", "C"}, {"A", "B", "C"}]
        pairs = compute_pairwise_jaccard(sets)
        assert len(pairs) == 3  # C(3,2) = 3 pairs
        for p in pairs:
            assert p["jaccard"] == 1.0
            assert p["n_intersection"] == 3
            assert p["n_union"] == 3

    def test_disjoint_sets(self):
        """No overlap -> Jaccard = 0.0."""
        sets = [{"A", "B"}, {"C", "D"}]
        pairs = compute_pairwise_jaccard(sets)
        assert len(pairs) == 1
        assert pairs[0]["jaccard"] == 0.0
        assert pairs[0]["n_intersection"] == 0
        assert pairs[0]["n_union"] == 4

    def test_partial_overlap(self):
        """Known overlap: {A,B,C} vs {B,C,D} -> Jaccard = 2/4 = 0.5."""
        sets = [{"A", "B", "C"}, {"B", "C", "D"}]
        pairs = compute_pairwise_jaccard(sets)
        assert len(pairs) == 1
        assert pairs[0]["jaccard"] == pytest.approx(0.5)
        assert pairs[0]["n_intersection"] == 2
        assert pairs[0]["n_union"] == 4

    def test_five_runs_pair_count(self):
        """5 runs -> C(5,2) = 10 pairs."""
        sets = [{"A"}, {"A", "B"}, {"B", "C"}, {"A", "C"}, {"A", "B", "C"}]
        pairs = compute_pairwise_jaccard(sets)
        assert len(pairs) == 10

    def test_five_runs_specific_pair(self):
        """Verify a specific pair from a 5-run set."""
        sets = [{"A"}, {"A", "B"}, {"B", "C"}, {"A", "C"}, {"A", "B", "C"}]
        pairs = compute_pairwise_jaccard(sets)
        # Pair (0,1): {A} vs {A,B} -> intersection={A}, union={A,B} -> 1/2=0.5
        pair_01 = [p for p in pairs if p["run_i"] == 1 and p["run_j"] == 2][0]
        assert pair_01["jaccard"] == pytest.approx(0.5)
        # Pair (2,4): {B,C} vs {A,B,C} -> intersection={B,C}, union={A,B,C} -> 2/3
        pair_24 = [p for p in pairs if p["run_i"] == 3 and p["run_j"] == 5][0]
        assert pair_24["jaccard"] == pytest.approx(2.0 / 3.0)

    def test_empty_sets(self):
        """Empty gene sets -> Jaccard = 0.0."""
        sets = [set(), set()]
        pairs = compute_pairwise_jaccard(sets)
        assert len(pairs) == 1
        assert pairs[0]["jaccard"] == 0.0

    def test_single_run_no_pairs(self):
        """One run -> no pairs."""
        sets = [{"A", "B"}]
        pairs = compute_pairwise_jaccard(sets)
        assert len(pairs) == 0

    def test_run_indices_are_1_based(self):
        """Run indices in output should be 1-based (matching run_001, etc)."""
        sets = [{"A"}, {"B"}, {"C"}]
        pairs = compute_pairwise_jaccard(sets)
        run_is = [p["run_i"] for p in pairs]
        run_js = [p["run_j"] for p in pairs]
        assert min(run_is) == 1
        assert max(run_js) == 3
