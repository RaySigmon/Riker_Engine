"""
Riker Engine - Ingestion module tests.

Phase 5: Gene database loading and HGNC symbol resolution.
Phase 6: Expression data normalization and validation.
Phase 7: GEO series matrix parser and probe-to-gene mapper.
Uses mini fixture files (no network access required).
"""

import math
from pathlib import Path
import warnings

import pytest
import numpy as np
import pandas as pd

from riker.ingestion.gene_db import GeneEntry, HGNCResolver, SeedGeneDB
from riker.ingestion.normalizer import (
    FoldChangeValidation,
    NormalizationResult,
    detect_log2_status,
    normalize_expression,
    validate_fold_changes,
)
from riker.ingestion.geo_parser import (
    GEOSeriesMatrix,
    PhenotypeExtractor,
    ProbeGeneMapper,
)

# Path to test fixtures
FIXTURES_DIR = Path(__file__).parent / "fixtures"
MINI_HGNC = FIXTURES_DIR / "mini_hgnc.txt"
MINI_SEEDS = FIXTURES_DIR / "mini_seed_genes.csv"
MINI_MATRIX = FIXTURES_DIR / "mini_series_matrix.txt"
MINI_PLATFORM = FIXTURES_DIR / "mini_platform.txt"


# ---------------------------------------------------------------------------
# 1. HGNCResolver with mini fixture
# ---------------------------------------------------------------------------

class TestHGNCResolver:
    """Test HGNC symbol resolution with mini fixture data."""

    @pytest.fixture
    def resolver(self):
        return HGNCResolver(hgnc_path=MINI_HGNC)

    def test_loads_approved_symbols(self, resolver):
        """Approved symbols should resolve to themselves."""
        assert resolver.resolve("ATP2B2") == "ATP2B2"
        assert resolver.resolve("SEZ6L2") == "SEZ6L2"
        assert resolver.resolve("TP53") == "TP53"
        assert resolver.resolve("BRCA1") == "BRCA1"

    def test_resolves_previous_symbols(self, resolver):
        """Previous symbols should map to current approved symbol."""
        assert resolver.resolve("PMCA2") == "ATP2B2"
        assert resolver.resolve("MRE11A") == "MRE11"
        assert resolver.resolve("BRCC1") == "BRCA1"
        assert resolver.resolve("ANKRD15") == "KANK1"
        assert resolver.resolve("MMAC1") == "PTEN"

    def test_resolves_aliases(self, resolver):
        """Alias symbols should map to current approved symbol."""
        assert resolver.resolve("BSRP-A") == "SEZ6L2"
        assert resolver.resolve("RNF53") == "BRCA1"
        assert resolver.resolve("FANCS") == "BRCA1"
        assert resolver.resolve("LFS1") == "TP53"
        assert resolver.resolve("TRP53") == "TP53"
        assert resolver.resolve("KANK") == "KANK1"

    def test_case_insensitive(self, resolver):
        """Resolution should be case-insensitive."""
        assert resolver.resolve("atp2b2") == "ATP2B2"
        assert resolver.resolve("Pmca2") == "ATP2B2"
        assert resolver.resolve("mre11a") == "MRE11"

    def test_unknown_symbol_unchanged(self, resolver):
        """Unknown symbols should pass through unchanged."""
        assert resolver.resolve("FAKE_GENE") == "FAKE_GENE"
        assert resolver.resolve("ZZZZZ999") == "ZZZZZ999"

    def test_empty_string(self, resolver):
        """Empty input should return empty."""
        assert resolver.resolve("") == ""
        assert resolver.resolve("  ") == ""

    def test_approved_priority_over_alias(self, resolver):
        """If a symbol is both approved AND an alias for another gene,
        the approved identity should win."""
        # MMAC1 is both a previous symbol for PTEN and listed under prev_symbol
        # The approved symbol PTEN should resolve to itself
        assert resolver.resolve("PTEN") == "PTEN"

    def test_is_approved(self, resolver):
        assert resolver.is_approved("ATP2B2")
        assert resolver.is_approved("TP53")
        assert not resolver.is_approved("PMCA2")  # previous, not approved
        assert not resolver.is_approved("FAKE_GENE")

    def test_batch_resolve(self, resolver):
        result = resolver.resolve_batch(["ATP2B2", "PMCA2", "FAKE_GENE"])
        assert result["ATP2B2"] == "ATP2B2"
        assert result["PMCA2"] == "ATP2B2"
        assert result["FAKE_GENE"] == "FAKE_GENE"

    def test_n_approved(self, resolver):
        """Mini fixture has 8 approved symbols."""
        assert resolver.n_approved == 8

    def test_total_mappings_greater_than_approved(self, resolver):
        """Total mappings should include aliases and previous symbols."""
        assert resolver.n_total_mappings > resolver.n_approved

    def test_pipe_delimited_previous_symbols(self, resolver):
        """KANK1 has pipe-delimited previous symbols: ANKRD15|CPSQ2."""
        assert resolver.resolve("ANKRD15") == "KANK1"
        assert resolver.resolve("CPSQ2") == "KANK1"

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            HGNCResolver(hgnc_path="/nonexistent/path.txt")


# ---------------------------------------------------------------------------
# 2. SeedGeneDB loading and resolution
# ---------------------------------------------------------------------------

class TestSeedGeneDB:
    """Test seed gene database loading with HGNC resolution."""

    @pytest.fixture
    def resolver(self):
        return HGNCResolver(hgnc_path=MINI_HGNC)

    def test_basic_load(self, resolver):
        """Load seed genes and verify count."""
        db = SeedGeneDB(
            MINI_SEEDS,
            symbol_column="symbol",
            tier_column="tier",
            chromosome_column="chromosome",
            source_name="TEST",
            resolver=resolver,
        )
        # 8 rows in CSV, but PMCA2 resolves to ATP2B2 (duplicate) -> 7 unique
        # and BRCC1 resolves to BRCA1
        # ATP2B2, SEZ6L2, KANK1, MRE11 (from MRE11A), ATP2B2 (from PMCA2 - dedup),
        # BRCA1 (from BRCC1), FAKE_GENE, TP53
        # = 7 unique after dedup
        assert db.n_genes == 7

    def test_remapped_count(self, resolver):
        """MRE11A->MRE11, PMCA2->ATP2B2 (deduped), BRCC1->BRCA1 = 2 remapped that survive dedup."""
        db = SeedGeneDB(
            MINI_SEEDS, symbol_column="symbol",
            resolver=resolver,
        )
        # MRE11A -> MRE11 (remapped, kept)
        # BRCC1 -> BRCA1 (remapped, kept)
        # PMCA2 -> ATP2B2 (remapped, but deduped because ATP2B2 already loaded)
        assert db.n_remapped == 2

    def test_resolved_symbols(self, resolver):
        db = SeedGeneDB(
            MINI_SEEDS, symbol_column="symbol",
            resolver=resolver,
        )
        symbols = db.symbols
        assert "ATP2B2" in symbols
        assert "SEZ6L2" in symbols
        assert "KANK1" in symbols
        assert "MRE11" in symbols  # resolved from MRE11A
        assert "BRCA1" in symbols  # resolved from BRCC1
        assert "TP53" in symbols
        assert "FAKE_GENE" in symbols  # unknown, kept as-is

        # Previous symbols should NOT appear as separate entries
        assert "MRE11A" not in symbols
        assert "PMCA2" not in symbols
        assert "BRCC1" not in symbols

    def test_contains(self, resolver):
        db = SeedGeneDB(
            MINI_SEEDS, symbol_column="symbol",
            resolver=resolver,
        )
        assert db.contains("ATP2B2")
        assert db.contains("atp2b2")  # case-insensitive
        assert db.contains("MRE11")
        assert not db.contains("MRE11A")  # original symbol, not in resolved set
        assert not db.contains("NONEXISTENT")

    def test_get_entry(self, resolver):
        db = SeedGeneDB(
            MINI_SEEDS, symbol_column="symbol",
            tier_column="tier",
            resolver=resolver,
        )
        entry = db.get_entry("MRE11")
        assert entry is not None
        assert entry.original_symbol == "MRE11A"
        assert entry.resolved_symbol == "MRE11"
        assert entry.was_remapped is True

    def test_tier_loaded(self, resolver):
        db = SeedGeneDB(
            MINI_SEEDS, symbol_column="symbol",
            tier_column="tier",
            resolver=resolver,
        )
        entry = db.get_entry("ATP2B2")
        assert entry is not None
        assert entry.tier == "1"

    def test_load_without_resolver(self):
        """Loading without resolver should keep original symbols."""
        db = SeedGeneDB(MINI_SEEDS, symbol_column="symbol")
        assert db.n_genes == 8  # no dedup since no resolution
        assert "MRE11A" in db.symbols
        assert "MRE11" not in db.symbols

    def test_file_not_found(self, resolver):
        with pytest.raises(FileNotFoundError):
            SeedGeneDB("/nonexistent/file.csv", resolver=resolver)

    def test_bad_column_name(self, resolver):
        with pytest.raises(ValueError, match="not found"):
            SeedGeneDB(
                MINI_SEEDS, symbol_column="wrong_column",
                resolver=resolver,
            )

    def test_source_name(self, resolver):
        db = SeedGeneDB(
            MINI_SEEDS, symbol_column="symbol",
            source_name="SFARI",
            resolver=resolver,
        )
        for gene in db.genes:
            assert gene.source == "SFARI"

    def test_get_remapped_genes(self, resolver):
        db = SeedGeneDB(
            MINI_SEEDS, symbol_column="symbol",
            resolver=resolver,
        )
        remapped = db.get_remapped_genes()
        original_symbols = [g.original_symbol for g in remapped]
        assert "MRE11A" in original_symbols
        assert "BRCC1" in original_symbols


# ===========================================================================
# PHASE 6: NORMALIZER TESTS
# ===========================================================================

# ---------------------------------------------------------------------------
# 3. Log2 detection
# ---------------------------------------------------------------------------

class TestLog2Detection:
    """Verify log2 transformation detection logic."""

    def test_raw_intensities_detected(self):
        """High median (>= 20) should be flagged as raw."""
        # Typical raw microarray: values in hundreds to thousands
        vals = np.random.uniform(100, 10000, size=1000)
        needs_log2, reason = detect_log2_status(vals)
        assert needs_log2 is True
        assert "raw intensities" in reason.lower()

    def test_log2_data_not_flagged(self):
        """Low median (< 20) without negatives should not be flagged."""
        # Typical log2-transformed: values 3-16
        vals = np.random.uniform(3.0, 16.0, size=1000)
        needs_log2, reason = detect_log2_status(vals)
        assert needs_log2 is False

    def test_negative_values_mean_already_log2(self):
        """Presence of negative values means already transformed."""
        vals = np.random.normal(0, 2, size=1000)  # centered near 0, has negatives
        needs_log2, reason = detect_log2_status(vals)
        assert needs_log2 is False
        assert "negative" in reason.lower()

    def test_empty_data(self):
        needs_log2, reason = detect_log2_status(np.array([]))
        assert needs_log2 is False

    def test_median_exactly_at_threshold(self):
        """Median exactly 20 should trigger transformation."""
        vals = np.array([20.0, 20.0, 20.0])
        needs_log2, _ = detect_log2_status(vals)
        assert needs_log2 is True

    def test_dataframe_input(self):
        """Should accept pandas DataFrame."""
        df = pd.DataFrame(np.random.uniform(100, 5000, size=(10, 5)))
        needs_log2, _ = detect_log2_status(df)
        assert needs_log2 is True


# ---------------------------------------------------------------------------
# 4. Expression normalization
# ---------------------------------------------------------------------------

class TestNormalization:
    """Verify expression normalization applies log2 correctly."""

    def test_raw_data_transformed(self):
        """Raw intensity data should be log2(x+1) transformed."""
        raw = np.array([0, 1, 100, 1000, 10000], dtype=float)
        result = normalize_expression(raw)
        assert result.was_transformed is True
        # log2(1001) ≈ 9.97, log2(10001) ≈ 13.29
        assert result.data[-1] == pytest.approx(np.log2(10001))

    def test_log2_data_untouched(self):
        """Already log2-transformed data should not be changed."""
        log2_data = np.random.uniform(3, 16, size=100)
        result = normalize_expression(log2_data)
        assert result.was_transformed is False
        np.testing.assert_array_equal(result.data, log2_data)

    def test_force_log2(self):
        """force_log2=True should transform regardless of detection."""
        data = np.array([5.0, 8.0, 12.0])  # looks like log2
        result = normalize_expression(data, force_log2=True)
        assert result.was_transformed is True
        assert result.data[0] == pytest.approx(np.log2(6.0))

    def test_force_no_log2(self):
        """force_log2=False should skip regardless of detection."""
        data = np.array([100.0, 500.0, 1000.0])  # looks like raw
        result = normalize_expression(data, force_log2=False)
        assert result.was_transformed is False

    def test_negative_background_subtracted_clamp_and_transform(self):
        """Raw intensities with negatives should clamp to 0 and log2-transform."""
        data = np.array([100.0, 200.0, -5.0, 300.0])
        result = normalize_expression(data)
        assert result.was_transformed is True
        assert result.data.min() >= 0  # negatives clamped

    def test_zeros_handled(self):
        """Zeros should be handled by log2(0+1) = 0."""
        data = np.array([0.0, 0.0, 100.0])
        result = normalize_expression(data)
        assert result.data[0] == 0.0  # log2(1) = 0
        assert result.data[1] == 0.0

    def test_dataframe_input(self):
        """Should accept DataFrame and return DataFrame with preserved index/columns."""
        df = pd.DataFrame(
            np.random.uniform(100, 5000, size=(5, 3)),
            index=[f"GENE_{i}" for i in range(5)],
            columns=[f"GSM{i}" for i in range(3)],
        )
        result = normalize_expression(df)
        assert result.was_transformed is True
        assert isinstance(result.data, pd.DataFrame)
        assert list(result.data.index) == list(df.index)
        assert list(result.data.columns) == list(df.columns)

    def test_metadata_populated(self):
        data = np.array([100.0, 200.0, 300.0])
        result = normalize_expression(data)
        assert result.original_median == 200.0
        assert result.transformed_median < 200.0  # log2 compresses
        assert result.n_negative_values == 0


# ---------------------------------------------------------------------------
# 5. Fold change validation (the 119 bug catcher)
# ---------------------------------------------------------------------------

class TestFoldChangeValidation:
    """Verify the range check that catches impossible fold changes."""

    def test_normal_fold_changes_pass(self):
        """Typical log2FC values (-2 to 2) should all pass."""
        fc = {f"GENE{i}": np.random.uniform(-2, 2) for i in range(100)}
        result = validate_fold_changes(fc)
        assert result.is_valid is True
        assert result.n_flagged == 0

    def test_impossible_fold_change_caught(self):
        """The exact bug from development: log2FC of 119."""
        fc = {"NORMAL_GENE": -0.5, "BUG_GENE": 119.0, "ANOTHER_NORMAL": 1.2}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = validate_fold_changes(fc)
            assert len(w) == 1
            assert "FOLD CHANGE RANGE CHECK FAILED" in str(w[0].message)
        assert result.is_valid is False
        assert result.n_flagged == 1
        assert result.flagged_genes[0] == ("BUG_GENE", 119.0)
        assert result.max_abs_log2fc == 119.0

    def test_threshold_boundary(self):
        """log2FC of exactly 10 should NOT be flagged (boundary is >)."""
        fc = {"GENE1": 10.0}
        result = validate_fold_changes(fc)
        assert result.is_valid is True

    def test_threshold_exceeded(self):
        """log2FC of 10.1 should be flagged."""
        fc = {"GENE1": 10.1}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = validate_fold_changes(fc)
        assert result.is_valid is False

    def test_negative_extreme_caught(self):
        """Large negative fold changes should also be caught."""
        fc = {"GENE1": -15.0}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = validate_fold_changes(fc)
        assert result.is_valid is False
        assert result.flagged_genes[0] == ("GENE1", -15.0)

    def test_multiple_offenders_sorted(self):
        """Multiple flagged genes should be sorted by |log2FC| descending."""
        fc = {"A": 50.0, "B": -100.0, "C": 12.0, "D": 0.5}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = validate_fold_changes(fc)
        assert result.n_flagged == 3
        assert result.flagged_genes[0][0] == "B"  # -100 is largest absolute
        assert result.flagged_genes[1][0] == "A"  # 50 is next

    def test_custom_threshold(self):
        """Custom threshold should be respected."""
        fc = {"GENE1": 5.5}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = validate_fold_changes(fc, threshold=5.0)
        assert result.is_valid is False
        assert result.threshold == 5.0

    def test_empty_input(self):
        result = validate_fold_changes({})
        assert result.is_valid is True
        assert result.n_checked == 0


# ===========================================================================
# PHASE 7: GEO PARSER TESTS
# ===========================================================================

# ---------------------------------------------------------------------------
# 6. GEO Series Matrix Parsing
# ---------------------------------------------------------------------------

class TestGEOSeriesMatrix:
    """Test GEO series matrix file parsing."""

    def test_basic_parse(self):
        parser = GEOSeriesMatrix(MINI_MATRIX)
        result = parser.get_result()
        assert result.accession == "GSE99999"
        assert result.platform == "GPL99999"
        assert result.n_samples == 6
        assert result.n_probes == 5

    def test_sample_ids(self):
        parser = GEOSeriesMatrix(MINI_MATRIX)
        result = parser.get_result()
        assert len(result.sample_ids) == 6
        assert "GSM100001" in result.sample_ids
        assert "GSM100006" in result.sample_ids

    def test_expression_matrix_shape(self):
        parser = GEOSeriesMatrix(MINI_MATRIX)
        result = parser.get_result()
        assert result.expression.shape == (5, 6)

    def test_expression_values(self):
        parser = GEOSeriesMatrix(MINI_MATRIX)
        result = parser.get_result()
        # PROBE_001, GSM100001 should be 8.5
        assert result.expression.loc["PROBE_001", "GSM100001"] == 8.5

    def test_metadata_extracted(self):
        parser = GEOSeriesMatrix(MINI_MATRIX)
        result = parser.get_result()
        assert "Sample_characteristics_ch1" in result.sample_metadata
        ch1 = result.sample_metadata["Sample_characteristics_ch1"]
        assert "GSM100001" in ch1
        assert "autism" in ch1["GSM100001"].lower()

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            GEOSeriesMatrix("/nonexistent/file.txt")


# ---------------------------------------------------------------------------
# 7. Phenotype Extraction
# ---------------------------------------------------------------------------

class TestPhenotypeExtractor:
    """Test phenotype (case/control) extraction from metadata."""

    def _get_metadata(self):
        parser = GEOSeriesMatrix(MINI_MATRIX)
        result = parser.get_result()
        return result.sample_metadata, result.sample_ids

    def test_auto_detection(self):
        """Auto-detect case/control from 'disease status' field."""
        metadata, sample_ids = self._get_metadata()
        extractor = PhenotypeExtractor()
        assignment = extractor.extract(metadata, sample_ids)

        assert assignment.n_cases == 3
        assert assignment.n_controls == 3
        assert assignment.method == "auto"
        assert assignment.groups["GSM100001"] == "case"
        assert assignment.groups["GSM100004"] == "control"

    def test_override_extraction(self):
        """Override with explicit field and value mappings."""
        metadata, sample_ids = self._get_metadata()
        extractor = PhenotypeExtractor(
            override_field="Sample_characteristics_ch1",
            override_case_values=["autism"],
            override_control_values=["control"],
        )
        assignment = extractor.extract(metadata, sample_ids)

        assert assignment.n_cases == 3
        assert assignment.n_controls == 3
        assert assignment.method == "override"
        assert assignment.field_used == "Sample_characteristics_ch1"

    def test_override_missing_field(self):
        """Override with non-existent field should raise."""
        metadata, sample_ids = self._get_metadata()
        extractor = PhenotypeExtractor(
            override_field="NONEXISTENT_FIELD",
            override_case_values=["case"],
            override_control_values=["control"],
        )
        with pytest.raises(ValueError, match="not found"):
            extractor.extract(metadata, sample_ids)

    def test_override_without_values_raises(self):
        """Override field without value mappings should raise."""
        metadata, sample_ids = self._get_metadata()
        extractor = PhenotypeExtractor(
            override_field="Sample_characteristics_ch1",
        )
        with pytest.raises(ValueError, match="override_case_values"):
            extractor.extract(metadata, sample_ids)


# ---------------------------------------------------------------------------
# 8. Probe-to-Gene Mapping
# ---------------------------------------------------------------------------

class TestProbeGeneMapper:
    """Test probe-to-gene mapping with platform annotation."""

    def test_basic_mapping(self):
        mapper = ProbeGeneMapper(MINI_PLATFORM)
        result = mapper.get_result()
        assert result.n_mapped == 5
        assert result.probe_to_gene["PROBE_001"] == "ATP2B2"
        assert result.probe_to_gene["PROBE_002"] == "SEZ6L2"

    def test_multi_gene_handled(self):
        """PROBE_004 has 'BRCA1 /// TP53' — should use first (BRCA1)."""
        mapper = ProbeGeneMapper(MINI_PLATFORM)
        result = mapper.get_result()
        assert result.probe_to_gene["PROBE_004"] == "BRCA1"
        assert result.n_multi_gene >= 1

    def test_hgnc_resolution(self):
        """MRE11A (PROBE_005) should resolve to MRE11 via HGNC."""
        resolver = HGNCResolver(hgnc_path=MINI_HGNC)
        mapper = ProbeGeneMapper(MINI_PLATFORM, resolver=resolver)
        result = mapper.get_result()
        assert result.probe_to_gene["PROBE_005"] == "MRE11"

    def test_without_resolver(self):
        """Without HGNC resolver, MRE11A stays as MRE11A."""
        mapper = ProbeGeneMapper(MINI_PLATFORM)
        result = mapper.get_result()
        assert result.probe_to_gene["PROBE_005"] == "MRE11A"

    def test_auto_detect_columns(self):
        """Should auto-detect 'ID' and 'Gene Symbol' columns."""
        mapper = ProbeGeneMapper(MINI_PLATFORM)
        result = mapper.get_result()
        assert result.gene_column_used == "Gene Symbol"

    def test_map_expression(self):
        """Collapse probe-level expression to gene-level."""
        parser = GEOSeriesMatrix(MINI_MATRIX)
        geo_result = parser.get_result()

        mapper = ProbeGeneMapper(MINI_PLATFORM)
        gene_expr = mapper.map_expression(geo_result.expression)

        assert "ATP2B2" in gene_expr.index
        assert "SEZ6L2" in gene_expr.index
        assert len(gene_expr) <= 5  # at most 5 genes (could be less with multi-probe)

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            ProbeGeneMapper("/nonexistent/annotation.txt")


# ===========================================================================
# snRNA-seq PSEUDO-BULKING TESTS
# ===========================================================================

from riker.ingestion.snrnaseq import (
    PseudoBulkResult,
    pseudo_bulk_from_counts,
)


def _make_snrnaseq_data(n_donors=6, n_nuclei_per_donor=50, n_genes=250,
                         seed=42):
    """Create simulated snRNA-seq count data.

    3 case donors, 3 control donors, 2 cell types.
    """
    np.random.seed(seed)
    rows = []

    donors = [f"DONOR_{i}" for i in range(n_donors)]
    conditions = ["ASD"] * (n_donors // 2) + ["Control"] * (n_donors // 2)
    cell_types = ["Excitatory", "Microglia"]
    gene_names = [f"GENE_{i}" for i in range(n_genes)]

    for donor, condition in zip(donors, conditions):
        for ct in cell_types:
            for nuc in range(n_nuclei_per_donor):
                row = {
                    "donor_id": donor,
                    "condition": condition,
                    "cell_type": ct,
                }
                # Generate count data with condition-specific shift for some genes
                for g_idx, gene in enumerate(gene_names):
                    base_count = np.random.poisson(5)
                    # First 10 genes are DE in excitatory neurons
                    if g_idx < 10 and ct == "Excitatory" and condition == "ASD":
                        base_count += np.random.poisson(3)
                    row[gene] = base_count
                rows.append(row)

    return pd.DataFrame(rows), gene_names


class TestPseudoBulkFromCounts:
    """Test pseudo-bulking from count matrices."""

    def test_basic_pseudo_bulk(self):
        counts, genes = _make_snrnaseq_data()
        result = pseudo_bulk_from_counts(
            counts, donor_col="donor_id", condition_col="condition",
            case_values=["ASD"], control_values=["Control"],
        )

        assert isinstance(result, PseudoBulkResult)
        assert result.n_donors == 6  # 6 donors survive
        assert result.cell_type == "all"
        assert result.expression.shape[1] == 6  # 6 samples

    def test_cell_type_filter(self):
        counts, genes = _make_snrnaseq_data()
        result = pseudo_bulk_from_counts(
            counts, donor_col="donor_id", condition_col="condition",
            cell_type_col="cell_type", target_cell_type="Excitatory",
            case_values=["ASD"], control_values=["Control"],
        )

        assert result.cell_type == "Excitatory"
        assert result.n_donors == 6

    def test_phenotype_assignment(self):
        counts, genes = _make_snrnaseq_data()
        result = pseudo_bulk_from_counts(
            counts, donor_col="donor_id", condition_col="condition",
            case_values=["ASD"], control_values=["Control"],
        )

        n_cases = sum(1 for v in result.phenotypes.values() if v == "case")
        n_controls = sum(1 for v in result.phenotypes.values() if v == "control")
        assert n_cases == 3
        assert n_controls == 3

    def test_normalized_output(self):
        counts, genes = _make_snrnaseq_data()
        result = pseudo_bulk_from_counts(
            counts, donor_col="donor_id", condition_col="condition",
            case_values=["ASD"], control_values=["Control"],
            normalize=True,
        )

        # After log2(CPM+1), values should be moderate (not raw counts)
        assert result.expression.max().max() < 25  # log2 scale
        assert result.expression.min().min() >= 0  # log2(x+1) >= 0

    def test_unnormalized_output(self):
        counts, genes = _make_snrnaseq_data()
        result = pseudo_bulk_from_counts(
            counts, donor_col="donor_id", condition_col="condition",
            case_values=["ASD"], control_values=["Control"],
            normalize=False,
        )

        # Raw sums should be larger than log2 values
        assert result.expression.max().max() > 25

    def test_min_nuclei_filter(self):
        counts, genes = _make_snrnaseq_data(n_nuclei_per_donor=5)
        # 5 nuclei per donor * 2 cell types = 10 total per donor
        # min_nuclei=20 should filter all
        with pytest.raises(ValueError, match="No valid pseudo-bulk samples produced"):
            pseudo_bulk_from_counts(
                counts, donor_col="donor_id", condition_col="condition",
                case_values=["ASD"], control_values=["Control"],
                min_nuclei=20,
            )

    def test_expression_matrix_format(self):
        """Output should be genes × samples (matching pipeline expectation)."""
        counts, genes = _make_snrnaseq_data()
        result = pseudo_bulk_from_counts(
            counts, donor_col="donor_id", condition_col="condition",
            case_values=["ASD"], control_values=["Control"],
        )

        # Rows = genes, columns = samples
        assert result.expression.index.name == "gene"
        assert all("DONOR" in col for col in result.expression.columns)

    def test_auto_detect_case_control(self):
        counts, genes = _make_snrnaseq_data()
        # Don't provide case/control values — should auto-detect
        result = pseudo_bulk_from_counts(
            counts, donor_col="donor_id", condition_col="condition",
        )

        n_cases = sum(1 for v in result.phenotypes.values() if v == "case")
        n_controls = sum(1 for v in result.phenotypes.values() if v == "control")
        assert n_cases == 3  # ASD
        assert n_controls == 3  # Control

    def test_missing_column_raises(self):
        counts, genes = _make_snrnaseq_data()
        with pytest.raises(ValueError, match="not found"):
            pseudo_bulk_from_counts(
                counts, donor_col="NONEXISTENT", condition_col="condition",
            )

    def test_sample_metadata_populated(self):
        counts, genes = _make_snrnaseq_data()
        result = pseudo_bulk_from_counts(
            counts, donor_col="donor_id", condition_col="condition",
            case_values=["ASD"], control_values=["Control"],
        )

        for sample_id, meta in result.sample_metadata.items():
            assert "donor" in meta
            assert "condition" in meta
            assert "n_nuclei" in meta
            assert meta["n_nuclei"] > 0

    def test_pipeline_integration(self):
        """Pseudo-bulk output should be directly usable by Phase 1."""
        from riker.phases.phase1_crossref import cross_reference_gene

        counts, genes = _make_snrnaseq_data(n_genes=250)
        result = pseudo_bulk_from_counts(
            counts, donor_col="donor_id", condition_col="condition",
            cell_type_col="cell_type", target_cell_type="Excitatory",
            case_values=["ASD"], control_values=["Control"],
        )

        # Feed directly into cross_reference_gene
        gene_result = cross_reference_gene(
            "GENE_0",  # one of the DE genes
            {"snRNA_DS": result.expression},
            {"snRNA_DS": result.phenotypes},
            min_datasets=1,
        )

        assert gene_result.n_datasets_detected == 1
        assert len(gene_result.de_results) == 1


# ---------------------------------------------------------------------------
# v0.2.0 Tests: Normalizer negative handling, probe validation, HGNC URL
# ---------------------------------------------------------------------------

class TestNormalizerNegativeValues:
    def test_background_subtracted_raw(self):
        """Raw intensities with a few negative values should be log2-transformed."""
        data = np.random.normal(200, 50, (100, 10))
        data[0, 0] = -5.0
        data[1, 1] = -22.0
        result = normalize_expression(data)
        assert result.was_transformed is True
        assert result.data.max() < 20
        assert result.data.min() >= 0

    def test_genuine_log2_with_negatives(self):
        """Already log2 data with negative fold changes should NOT be transformed."""
        data = np.random.normal(8.0, 2.0, (100, 10))
        data[0, 0] = -1.5
        result = normalize_expression(data)
        assert result.was_transformed is False

    def test_normalize_preserves_dataframe(self):
        """Normalization should preserve DataFrame index and columns."""
        data = pd.DataFrame(
            np.random.uniform(100, 50000, (10, 5)),
            index=[f"GENE_{i}" for i in range(10)],
            columns=[f"GSM{i}" for i in range(5)],
        )
        result = normalize_expression(data)
        assert result.was_transformed is True
        assert isinstance(result.data, pd.DataFrame)
        assert list(result.data.index) == list(data.index)
        assert list(result.data.columns) == list(data.columns)


class TestHGNCURLFallback:
    def test_hgnc_urls_has_multiple(self):
        """HGNC_URLS should contain multiple mirror URLs."""
        from riker.ingestion.gene_db import HGNC_URLS
        assert len(HGNC_URLS) >= 2
        assert "storage.googleapis.com" in HGNC_URLS[0]
