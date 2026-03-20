# INSTRUCTION SET FOR KAI — PHASE 5: `riker/ingestion/gene_db.py`

## References
- Blueprint Section 5.1 (Gene Database Ingestion)
- Blueprint Section 14 (Technology Stack — HGNC complete set TSV)
- Context Transfer: "Probe-to-Gene Mapping Requires HGNC Alias Resolution"
- Context Transfer: "Without alias resolution via the HGNC complete set, 10-15% of seed genes go undetected"

## WHY THIS MODULE MATTERS

Gene symbols change over time. SFARI might list a gene as "MRE11A" but the
GEO platform annotation file calls it "MRE11". Without alias resolution,
that gene is silently missed during cross-referencing. The HGNC (HUGO Gene
Nomenclature Committee) maintains a complete set of approved symbols,
previous symbols, and aliases. This module downloads that set once, caches
it, and builds a lookup table that maps ANY known symbol (current, previous,
or alias) to the current approved symbol.

## CRITICAL REQUIREMENTS

1. `HGNCResolver` class that:
   - Downloads the HGNC complete set TSV from genenames.org (or loads from cache)
   - Builds a dict mapping: previous_symbol -> current_symbol, alias_symbol -> current_symbol
   - The `resolve(symbol)` method returns the current approved symbol, or the
     original symbol if no mapping is found (do NOT drop unmapped genes)
   - Handles case-insensitive matching (gene symbols can appear in mixed case)
   - Caches the downloaded TSV to a configurable path (default: ~/.riker/hgnc_complete_set.txt)

2. `SeedGeneDB` class that:
   - Loads seed genes from a CSV file with configurable column names
   - Applies HGNC resolution to all loaded symbols
   - Tracks: original_symbol, resolved_symbol, confidence_tier (optional), chromosome (optional)
   - Reports how many symbols were remapped and which ones

3. Both classes must work OFFLINE after first download. The HGNC file is ~15MB
   and should be downloaded once, not on every run.

4. For testing, create a SMALL fixture file (not the real 15MB HGNC file).
   Tests must run without network access.

5. DO NOT modify any files in riker/stats/ or existing tests.

---

## FILE: `riker/ingestion/gene_db.py`

Write the following file at `/home/kai001/riker-engine/riker/ingestion/gene_db.py`:

```python
"""
Riker Engine - Seed gene database loading and HGNC symbol resolution.

Downloads and caches the HGNC complete set to build a symbol lookup table
that maps previous symbols and aliases to current approved symbols. This
resolves the 10-15% gene detection loss identified during Project Riker
development.

References:
    Blueprint Section 5.1 (Gene Database Ingestion)
    Context Transfer: "Probe-to-Gene Mapping Requires HGNC Alias Resolution"
"""

import csv
import logging
import os
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# HGNC complete set download URL
HGNC_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/"
    "hgnc_complete_set.txt"
)

# Default cache directory
DEFAULT_CACHE_DIR = Path.home() / ".riker"


@dataclass
class GeneEntry:
    """A single gene in the seed database.

    Attributes:
        original_symbol: The symbol as it appeared in the source file.
        resolved_symbol: The HGNC-resolved current approved symbol.
        was_remapped: True if original != resolved.
        tier: Confidence tier from the source database (optional).
        chromosome: Chromosomal location (optional).
        source: Which database this gene came from (e.g., 'SFARI', 'ADSP').
    """
    original_symbol: str
    resolved_symbol: str
    was_remapped: bool = False
    tier: str = ""
    chromosome: str = ""
    source: str = ""


class HGNCResolver:
    """Resolves gene symbols to current HGNC approved symbols.

    Downloads the HGNC complete set once and caches it locally.
    Builds a lookup: {any_known_symbol_upper -> approved_symbol}.

    Priority order for conflicts:
    1. Approved symbol (identity mapping, highest priority)
    2. Previous symbol
    3. Alias symbol

    Parameters
    ----------
    cache_dir : Path or str
        Directory to cache the HGNC file. Default: ~/.riker/
    hgnc_path : Path or str or None
        Explicit path to a pre-downloaded HGNC TSV file. If provided,
        skips download entirely. Use this for testing or air-gapped
        environments.
    """

    def __init__(
        self,
        cache_dir: Path | str = DEFAULT_CACHE_DIR,
        hgnc_path: Path | str | None = None,
    ):
        self._lookup: dict[str, str] = {}
        self._approved: set[str] = set()

        if hgnc_path is not None:
            self._load_from_file(Path(hgnc_path))
        else:
            cache_dir = Path(cache_dir)
            cached_file = cache_dir / "hgnc_complete_set.txt"
            if cached_file.exists():
                self._load_from_file(cached_file)
            else:
                self._download_and_cache(cache_dir, cached_file)

    def _download_and_cache(self, cache_dir: Path, cached_file: Path) -> None:
        """Download HGNC complete set and cache locally."""
        import requests

        cache_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Downloading HGNC complete set from {HGNC_URL}...")

        try:
            resp = requests.get(HGNC_URL, timeout=120)
            resp.raise_for_status()
        except Exception as e:
            raise RuntimeError(
                f"Failed to download HGNC complete set: {e}. "
                f"You can manually download from {HGNC_URL} and pass "
                f"the path via hgnc_path parameter."
            ) from e

        cached_file.write_text(resp.text, encoding="utf-8")
        logger.info(f"HGNC complete set cached at {cached_file}")
        self._load_from_file(cached_file)

    def _load_from_file(self, path: Path) -> None:
        """Parse HGNC TSV and build the lookup table."""
        if not path.exists():
            raise FileNotFoundError(f"HGNC file not found: {path}")

        logger.info(f"Loading HGNC data from {path}...")

        df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)

        # Required column: symbol (approved symbol)
        if "symbol" not in df.columns:
            # Try alternate column name
            if "Approved symbol" in df.columns:
                df = df.rename(columns={"Approved symbol": "symbol"})
            else:
                raise ValueError(
                    f"HGNC file missing 'symbol' column. "
                    f"Found columns: {list(df.columns)[:10]}"
                )

        # Step 1: Register all approved symbols (identity mapping, highest priority)
        for _, row in df.iterrows():
            approved = str(row.get("symbol", "")).strip()
            if approved:
                key = approved.upper()
                self._lookup[key] = approved
                self._approved.add(key)

        # Step 2: Register previous symbols (lower priority than approved)
        prev_col = None
        for col_name in ["prev_symbol", "Previous symbols", "Previous symbol"]:
            if col_name in df.columns:
                prev_col = col_name
                break

        if prev_col:
            for _, row in df.iterrows():
                approved = str(row.get("symbol", "")).strip()
                prev_raw = str(row.get(prev_col, "")).strip()
                if not approved or not prev_raw or prev_raw == "nan":
                    continue
                for prev in self._split_symbols(prev_raw):
                    key = prev.upper()
                    # Don't overwrite approved symbols
                    if key not in self._approved:
                        self._lookup[key] = approved

        # Step 3: Register alias symbols (lowest priority)
        alias_col = None
        for col_name in ["alias_symbol", "Alias symbols", "Alias symbol"]:
            if col_name in df.columns:
                alias_col = col_name
                break

        if alias_col:
            for _, row in df.iterrows():
                approved = str(row.get("symbol", "")).strip()
                alias_raw = str(row.get(alias_col, "")).strip()
                if not approved or not alias_raw or alias_raw == "nan":
                    continue
                for alias in self._split_symbols(alias_raw):
                    key = alias.upper()
                    # Don't overwrite approved or previous symbols
                    if key not in self._lookup:
                        self._lookup[key] = approved

        n_approved = len(self._approved)
        n_total = len(self._lookup)
        logger.info(
            f"HGNC resolver loaded: {n_approved} approved symbols, "
            f"{n_total - n_approved} aliases/previous symbols mapped."
        )

    @staticmethod
    def _split_symbols(raw: str) -> list[str]:
        """Split a pipe-delimited or comma-delimited symbol string."""
        # HGNC uses pipe-delimited in TSV downloads
        if "|" in raw:
            return [s.strip().strip('"') for s in raw.split("|") if s.strip()]
        elif "," in raw:
            return [s.strip().strip('"') for s in raw.split(",") if s.strip()]
        else:
            return [raw.strip().strip('"')]

    def resolve(self, symbol: str) -> str:
        """Resolve a gene symbol to current HGNC approved symbol.

        Parameters
        ----------
        symbol : str
            Any gene symbol (current, previous, or alias).

        Returns
        -------
        str
            The current approved HGNC symbol if found, otherwise
            the original symbol unchanged.
        """
        if not symbol or not symbol.strip():
            return symbol
        key = symbol.strip().upper()
        return self._lookup.get(key, symbol.strip())

    def is_approved(self, symbol: str) -> bool:
        """Check if a symbol is a current approved HGNC symbol."""
        return symbol.strip().upper() in self._approved

    def resolve_batch(self, symbols: list[str]) -> dict[str, str]:
        """Resolve multiple symbols at once.

        Returns
        -------
        dict
            Mapping of original_symbol -> resolved_symbol.
        """
        return {s: self.resolve(s) for s in symbols}

    @property
    def n_approved(self) -> int:
        """Number of approved symbols in the lookup."""
        return len(self._approved)

    @property
    def n_total_mappings(self) -> int:
        """Total number of symbol mappings (approved + aliases + previous)."""
        return len(self._lookup)


class SeedGeneDB:
    """Load and standardize a seed gene database.

    Reads a CSV file of disease-associated genes, resolves all symbols
    via HGNC, and provides the standardized gene list for downstream
    pipeline phases.

    Parameters
    ----------
    csv_path : str or Path
        Path to the seed gene CSV file.
    symbol_column : str
        Column name containing gene symbols. Default: 'symbol'.
    tier_column : str or None
        Column name for confidence tiers. Default: None.
    chromosome_column : str or None
        Column name for chromosomal location. Default: None.
    source_name : str
        Name of the gene database (e.g., 'SFARI', 'ADSP'). Default: ''.
    resolver : HGNCResolver or None
        Pre-built resolver. If None, symbols are used as-is.
    """

    def __init__(
        self,
        csv_path: str | Path,
        symbol_column: str = "symbol",
        tier_column: str | None = None,
        chromosome_column: str | None = None,
        source_name: str = "",
        resolver: HGNCResolver | None = None,
    ):
        self.csv_path = Path(csv_path)
        self.source_name = source_name
        self.genes: list[GeneEntry] = []
        self._symbol_to_entry: dict[str, GeneEntry] = {}

        self._load(
            symbol_column=symbol_column,
            tier_column=tier_column,
            chromosome_column=chromosome_column,
            resolver=resolver,
        )

    def _load(
        self,
        symbol_column: str,
        tier_column: str | None,
        chromosome_column: str | None,
        resolver: HGNCResolver | None,
    ) -> None:
        """Load and resolve seed genes from CSV."""
        if not self.csv_path.exists():
            raise FileNotFoundError(
                f"Seed gene file not found: {self.csv_path}"
            )

        df = pd.read_csv(self.csv_path, dtype=str)

        if symbol_column not in df.columns:
            raise ValueError(
                f"Column '{symbol_column}' not found in {self.csv_path}. "
                f"Available columns: {list(df.columns)}"
            )

        n_remapped = 0
        seen_symbols: set[str] = set()

        for _, row in df.iterrows():
            original = str(row[symbol_column]).strip()
            if not original or original == "nan":
                continue

            # Resolve symbol
            if resolver is not None:
                resolved = resolver.resolve(original)
            else:
                resolved = original

            # Deduplicate by resolved symbol
            resolved_upper = resolved.upper()
            if resolved_upper in seen_symbols:
                continue
            seen_symbols.add(resolved_upper)

            was_remapped = (resolved.upper() != original.upper())
            if was_remapped:
                n_remapped += 1

            tier = ""
            if tier_column and tier_column in df.columns:
                tier = str(row.get(tier_column, "")).strip()
                if tier == "nan":
                    tier = ""

            chrom = ""
            if chromosome_column and chromosome_column in df.columns:
                chrom = str(row.get(chromosome_column, "")).strip()
                if chrom == "nan":
                    chrom = ""

            entry = GeneEntry(
                original_symbol=original,
                resolved_symbol=resolved,
                was_remapped=was_remapped,
                tier=tier,
                chromosome=chrom,
                source=self.source_name,
            )
            self.genes.append(entry)
            self._symbol_to_entry[resolved_upper] = entry

        logger.info(
            f"Loaded {len(self.genes)} seed genes from {self.csv_path}. "
            f"{n_remapped} symbols remapped via HGNC."
        )
        if n_remapped > 0:
            remapped = [g for g in self.genes if g.was_remapped]
            for g in remapped[:10]:
                logger.info(
                    f"  Remapped: {g.original_symbol} -> {g.resolved_symbol}"
                )
            if len(remapped) > 10:
                logger.info(f"  ... and {len(remapped) - 10} more.")

    @property
    def symbols(self) -> list[str]:
        """List of resolved gene symbols."""
        return [g.resolved_symbol for g in self.genes]

    @property
    def n_genes(self) -> int:
        """Total number of unique seed genes."""
        return len(self.genes)

    @property
    def n_remapped(self) -> int:
        """Number of genes that were remapped via HGNC."""
        return sum(1 for g in self.genes if g.was_remapped)

    def contains(self, symbol: str) -> bool:
        """Check if a gene symbol is in the seed database (case-insensitive)."""
        return symbol.strip().upper() in self._symbol_to_entry

    def get_entry(self, symbol: str) -> GeneEntry | None:
        """Get the GeneEntry for a symbol (case-insensitive)."""
        return self._symbol_to_entry.get(symbol.strip().upper())

    def get_remapped_genes(self) -> list[GeneEntry]:
        """Get all genes that were remapped via HGNC."""
        return [g for g in self.genes if g.was_remapped]
```

---

## TEST FIXTURES

Before writing tests, create the fixture files that tests will use.

### Fixture 1: Mini HGNC file

Create file at `/home/kai001/riker-engine/tests/fixtures/mini_hgnc.txt` with this exact content:

```
hgnc_id	symbol	name	prev_symbol	alias_symbol
HGNC:1	A1BG	alpha-1-B glycoprotein		
HGNC:100	ATP2B2	ATPase plasma membrane Ca2+ transporting 2	PMCA2	ATP2B2A
HGNC:200	SEZ6L2	seizure related 6 homolog like 2		BSRP-A
HGNC:300	KANK1	KN motif and ankyrin repeat domains 1	ANKRD15|CPSQ2	KANK
HGNC:400	BRCA1	BRCA1 DNA repair associated	BRCC1	RNF53|FANCS
HGNC:500	TP53	tumor protein p53		LFS1|TRP53
HGNC:600	MRE11	MRE11 homologous recombination 11	MRE11A	HNGS1
HGNC:700	PTEN	phosphatase and tensin homolog	MMAC1	TEP1|MMAC1
```

### Fixture 2: Mini seed gene CSV

Create file at `/home/kai001/riker-engine/tests/fixtures/mini_seed_genes.csv` with this exact content:

```
symbol,tier,chromosome
ATP2B2,1,3p25.3
SEZ6L2,1,16p11.2
KANK1,2,9p24.3
MRE11A,2,11q21
PMCA2,1,3p25.3
BRCC1,3,17q21.31
FAKE_GENE,3,1p36.33
TP53,1,17p13.1
```

Note the test data: MRE11A is a previous symbol for MRE11, PMCA2 is a previous symbol for ATP2B2,
BRCC1 is a previous symbol for BRCA1, and FAKE_GENE has no HGNC mapping.

---

## TESTS: Create `tests/test_ingestion.py`

**Replace** the contents of `/home/kai001/riker-engine/tests/test_ingestion.py` with:

```python
"""
Riker Engine - Ingestion module tests.

Phase 5: Gene database loading and HGNC symbol resolution.
Uses mini fixture files (no network access required).
"""

import math
from pathlib import Path

import pytest

from riker.ingestion.gene_db import GeneEntry, HGNCResolver, SeedGeneDB

# Path to test fixtures
FIXTURES_DIR = Path(__file__).parent / "fixtures"
MINI_HGNC = FIXTURES_DIR / "mini_hgnc.txt"
MINI_SEEDS = FIXTURES_DIR / "mini_seed_genes.csv"


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
```

---

## EXECUTION INSTRUCTIONS

1. First, create the fixtures directory and files:

```bash
mkdir -p /home/kai001/riker-engine/tests/fixtures
```

2. Then write gene_db.py, the fixture files, and test_ingestion.py as specified above.

3. Run:

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_ingestion.py -v 2>&1
```

**Expected: ALL tests pass.** If any fail, report FULL output — do NOT modify tests or code without reporting first.

4. Then run this confirmation:

```bash
cd /home/kai001/riker-engine && python -c "
from pathlib import Path
from riker.ingestion.gene_db import HGNCResolver, SeedGeneDB

fixtures = Path('tests/fixtures')
resolver = HGNCResolver(hgnc_path=fixtures / 'mini_hgnc.txt')

print(f'HGNC resolver: {resolver.n_approved} approved, {resolver.n_total_mappings} total mappings')
print()

# Test the key remappings from both ASD and AD runs
test_cases = [
    ('ATP2B2', 'ATP2B2'),   # approved -> itself
    ('PMCA2', 'ATP2B2'),    # previous -> approved
    ('MRE11A', 'MRE11'),    # previous -> approved
    ('BRCC1', 'BRCA1'),     # previous -> approved
    ('BSRP-A', 'SEZ6L2'),   # alias -> approved
    ('FAKE', 'FAKE'),        # unknown -> unchanged
]

all_pass = True
for input_sym, expected in test_cases:
    result = resolver.resolve(input_sym)
    status = 'OK' if result == expected else f'FAIL (got {result})'
    if result != expected:
        all_pass = False
    print(f'  {input_sym:12s} -> {result:12s} {status}')

print()

# Load seed genes
db = SeedGeneDB(
    fixtures / 'mini_seed_genes.csv',
    symbol_column='symbol',
    tier_column='tier',
    source_name='TEST',
    resolver=resolver,
)
print(f'Seed DB: {db.n_genes} genes loaded, {db.n_remapped} remapped')
print(f'Symbols: {db.symbols}')

assert all_pass, 'FAIL: some resolution tests failed!'
assert db.n_genes == 7, f'FAIL: expected 7 unique genes, got {db.n_genes}'
print()
print('PASS: gene_db.py working correctly')
"
```

5. Also verify the stats tests still pass (regression check):

```bash
cd /home/kai001/riker-engine && python -m pytest tests/test_stats.py -q 2>&1
```

Report all three outputs back to the architect for QA review.
