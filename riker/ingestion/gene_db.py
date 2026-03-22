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

# HGNC complete set download URLs (try in order)
HGNC_URLS = [
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
    "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",
]

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
        """Download HGNC complete set and cache locally.

        Tries multiple mirror URLs in order. Falls back gracefully
        if the primary mirror is unavailable.
        """
        import requests

        cache_dir.mkdir(parents=True, exist_ok=True)

        for url in HGNC_URLS:
            try:
                logger.info(f"Downloading HGNC complete set from {url}...")
                resp = requests.get(url, timeout=120)
                resp.raise_for_status()
                cached_file.write_text(resp.text, encoding="utf-8")
                logger.info(f"HGNC complete set cached at {cached_file}")
                self._load_from_file(cached_file)
                return
            except Exception as e:
                logger.warning(f"HGNC download failed from {url}: {e}")
                continue

        raise RuntimeError(
            "Could not download HGNC data from any mirror. "
            "Download manually from https://www.genenames.org/download/archive/ "
            "and set hgnc_path in your config."
        )

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
            return ""
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
