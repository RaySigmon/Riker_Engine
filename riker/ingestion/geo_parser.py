"""
Riker Engine - GEO series matrix parser and probe-to-gene mapper.

Parses GEO series matrix files to extract expression data and sample
metadata. Handles the many format variations encountered across GEO
datasets, with per-dataset override support for metadata extraction.

References:
    Blueprint Section 5.2 (Transcriptomic Data Ingestion)
    Context Transfer: "GEO Metadata Parsing Is Fragile"
    Context Transfer: "Platform Annotation Files Vary Wildly"
"""

import logging
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ParsedGEOMatrix:
    """Result of parsing a GEO series matrix file.

    Attributes:
        accession: GEO series accession (e.g., 'GSE28521').
        expression: Expression DataFrame (probes × samples).
        sample_metadata: Dict of metadata_field -> {sample_id: value}.
        sample_ids: List of sample IDs (GSM accessions) in column order.
        platform: Platform accession (e.g., 'GPL570').
        n_probes: Number of probes in expression matrix.
        n_samples: Number of samples.
    """
    accession: str
    expression: pd.DataFrame
    sample_metadata: dict
    sample_ids: list
    platform: str
    n_probes: int
    n_samples: int


@dataclass(frozen=True)
class PhenotypeAssignment:
    """Result of phenotype extraction.

    Attributes:
        groups: Dict mapping sample_id -> 'case' or 'control'.
        n_cases: Number of case samples.
        n_controls: Number of control samples.
        field_used: Which metadata field was used for assignment.
        method: 'auto' or 'override'.
    """
    groups: dict
    n_cases: int
    n_controls: int
    field_used: str
    method: str


@dataclass(frozen=True)
class ProbeMapping:
    """Result of probe-to-gene mapping.

    Attributes:
        probe_to_gene: Dict mapping probe_id -> gene_symbol.
        n_mapped: Number of probes successfully mapped.
        n_unmapped: Number of probes with no gene symbol.
        n_multi_gene: Number of probes mapping to multiple genes.
        gene_column_used: Which column was used for gene symbols.
    """
    probe_to_gene: dict
    n_mapped: int
    n_unmapped: int
    n_multi_gene: int
    gene_column_used: str


# ---------------------------------------------------------------------------
# GEO Series Matrix Parser
# ---------------------------------------------------------------------------

class GEOSeriesMatrix:
    """Parse a GEO series matrix file.

    GEO series matrix files have a header section (lines starting with !)
    containing sample metadata, followed by an expression matrix between
    !series_matrix_table_begin and !series_matrix_table_end markers.

    Parameters
    ----------
    filepath : str or Path
        Path to the series matrix file (.txt or .txt.gz).
    """

    def __init__(self, filepath: str | Path):
        self.filepath = Path(filepath)
        if not self.filepath.exists():
            raise FileNotFoundError(
                f"Series matrix file not found: {self.filepath}"
            )

        self._metadata_lines: list[str] = []
        self._expression_df: pd.DataFrame | None = None
        self._sample_ids: list[str] = []
        self._accession: str = ""
        self._platform: str = ""
        self._sample_metadata: dict[str, dict[str, str]] = {}

        self._parse()

    def _parse(self) -> None:
        """Parse the series matrix file into metadata + expression."""
        import gzip

        open_func = gzip.open if str(self.filepath).endswith(".gz") else open
        open_kwargs = {"mode": "rt", "encoding": "utf-8"} if str(self.filepath).endswith(".gz") else {"mode": "r", "encoding": "utf-8"}

        header_lines = []
        data_lines = []
        in_table = False

        with open_func(self.filepath, **open_kwargs) as f:
            for line in f:
                line = line.rstrip("\n\r")
                if line.startswith("!series_matrix_table_begin"):
                    in_table = True
                    continue
                elif line.startswith("!series_matrix_table_end"):
                    in_table = False
                    continue

                if in_table:
                    data_lines.append(line)
                elif line.startswith("!"):
                    header_lines.append(line)

        # Parse metadata from header
        self._parse_header(header_lines)

        # Parse expression matrix from data
        if data_lines:
            self._parse_expression(data_lines)

    def _parse_header(self, lines: list[str]) -> None:
        """Extract sample metadata from header lines."""
        for line in lines:
            # Extract accession
            if line.startswith("!Series_geo_accession"):
                match = re.search(r'"(GSE\d+)"', line)
                if match:
                    self._accession = match.group(1)

            # Extract platform
            if line.startswith("!Series_platform_id") or line.startswith("!Platform_geo_accession"):
                match = re.search(r'"(GPL\d+)"', line)
                if match:
                    self._platform = match.group(1)

            # Extract sample IDs
            if line.startswith("!Sample_geo_accession"):
                self._sample_ids = self._extract_quoted_values(line)

            # Store all Sample_ metadata fields
            if line.startswith("!Sample_"):
                field_name = line.split("\t")[0].lstrip("!")
                values = self._extract_quoted_values(line)
                if values:
                    if field_name not in self._sample_metadata:
                        self._sample_metadata[field_name] = {}

                    # Handle multi-row fields (characteristics can repeat)
                    if self._sample_ids and len(values) == len(self._sample_ids):
                        existing = self._sample_metadata[field_name]
                        for sid, val in zip(self._sample_ids, values):
                            if sid in existing:
                                existing[sid] = existing[sid] + " | " + val
                            else:
                                existing[sid] = val

    def _parse_expression(self, lines: list[str]) -> None:
        """Parse the expression matrix from data lines."""
        if not lines:
            return

        # First line is header (probe ID + sample IDs)
        header = lines[0].split("\t")
        probe_col = header[0].strip('"')
        sample_cols = [c.strip('"') for c in header[1:]]

        # If we don't have sample IDs from metadata, use expression header
        if not self._sample_ids:
            self._sample_ids = sample_cols

        # Parse data rows
        rows = []
        index_vals = []
        for line in lines[1:]:
            if not line.strip():
                continue
            parts = line.split("\t")
            probe_id = parts[0].strip('"')
            try:
                values = [float(v.strip('"')) if v.strip('"') not in ("", "null", "NA", "NaN") else np.nan for v in parts[1:]]
            except (ValueError, IndexError):
                continue

            if len(values) == len(sample_cols):
                index_vals.append(probe_id)
                rows.append(values)

        if rows:
            self._expression_df = pd.DataFrame(
                rows,
                index=index_vals,
                columns=sample_cols,
                dtype=np.float64,
            )
            self._expression_df.index.name = "ID_REF"

    @staticmethod
    def _extract_quoted_values(line: str) -> list[str]:
        """Extract tab-separated quoted values from a metadata line."""
        parts = line.split("\t")
        if len(parts) < 2:
            return []
        return [p.strip('"').strip() for p in parts[1:] if p.strip()]

    def get_result(self) -> ParsedGEOMatrix:
        """Return the parsed result as a dataclass."""
        expr = self._expression_df if self._expression_df is not None else pd.DataFrame()
        return ParsedGEOMatrix(
            accession=self._accession,
            expression=expr,
            sample_metadata=self._sample_metadata,
            sample_ids=self._sample_ids,
            platform=self._platform,
            n_probes=len(expr),
            n_samples=len(self._sample_ids),
        )


# ---------------------------------------------------------------------------
# Phenotype Extractor
# ---------------------------------------------------------------------------

# Common field names where disease status appears
_PHENOTYPE_FIELDS = [
    "Sample_characteristics_ch1",
    "Sample_characteristics_ch2",
    "Sample_source_name_ch1",
    "Sample_title",
    "Sample_description",
]

# Common disease-related keywords to search for
_DISEASE_KEYWORDS = [
    "disease", "diagnosis", "status", "condition", "group",
    "phenotype", "tissue", "state",
]


class PhenotypeExtractor:
    """Extract case/control assignments from GEO sample metadata.

    Supports automatic detection of phenotype fields and per-dataset
    override configuration.

    Parameters
    ----------
    override_field : str or None
        Exact metadata field name to use. Bypasses auto-detection.
    override_case_values : list of str or None
        Values that indicate case samples (case-insensitive).
    override_control_values : list of str or None
        Values that indicate control samples (case-insensitive).
    """

    def __init__(
        self,
        override_field: str | None = None,
        override_case_values: list[str] | None = None,
        override_control_values: list[str] | None = None,
    ):
        self.override_field = override_field
        self.override_case_values = [v.lower() for v in (override_case_values or [])]
        self.override_control_values = [v.lower() for v in (override_control_values or [])]

    def extract(
        self,
        sample_metadata: dict[str, dict[str, str]],
        sample_ids: list[str],
    ) -> PhenotypeAssignment:
        """Assign case/control labels to samples.

        Parameters
        ----------
        sample_metadata : dict
            From ParsedGEOMatrix.sample_metadata.
            Structure: {field_name: {sample_id: value}}.
        sample_ids : list
            All sample IDs to assign.

        Returns
        -------
        PhenotypeAssignment
            Contains group assignments and metadata.
        """
        if self.override_field:
            return self._extract_with_override(sample_metadata, sample_ids)
        else:
            return self._extract_auto(sample_metadata, sample_ids)

    def _extract_with_override(
        self,
        sample_metadata: dict,
        sample_ids: list,
    ) -> PhenotypeAssignment:
        """Extract using explicit field and value overrides."""
        if self.override_field not in sample_metadata:
            raise ValueError(
                f"Override field '{self.override_field}' not found in metadata. "
                f"Available fields: {list(sample_metadata.keys())}"
            )

        if not self.override_case_values or not self.override_control_values:
            raise ValueError(
                "When using override_field, both override_case_values and "
                "override_control_values must be provided."
            )

        field_data = sample_metadata[self.override_field]
        groups = {}

        for sid in sample_ids:
            val = field_data.get(sid, "").lower().strip()
            if any(cv in val for cv in self.override_case_values):
                groups[sid] = "case"
            elif any(cv in val for cv in self.override_control_values):
                groups[sid] = "control"
            else:
                warnings.warn(
                    f"Sample {sid} value '{val}' did not match any "
                    f"case ({self.override_case_values}) or "
                    f"control ({self.override_control_values}) values.",
                    UserWarning,
                    stacklevel=3,
                )

        n_cases = sum(1 for g in groups.values() if g == "case")
        n_controls = sum(1 for g in groups.values() if g == "control")

        return PhenotypeAssignment(
            groups=groups,
            n_cases=n_cases,
            n_controls=n_controls,
            field_used=self.override_field,
            method="override",
        )

    def _extract_auto(
        self,
        sample_metadata: dict,
        sample_ids: list,
    ) -> PhenotypeAssignment:
        """Auto-detect phenotype field and assign groups."""
        # Try each candidate field
        for field_name in _PHENOTYPE_FIELDS:
            if field_name not in sample_metadata:
                continue

            field_data = sample_metadata[field_name]
            if not field_data:
                continue

            # Check if this field contains disease-related keywords
            all_values = [str(v).lower() for v in field_data.values()]
            combined = " ".join(all_values)

            has_keyword = any(kw in combined for kw in _DISEASE_KEYWORDS)
            if not has_keyword:
                continue

            # Try to identify case/control groups from unique values
            unique_vals = set(all_values)
            groups = self._try_assign_groups(field_data, sample_ids, unique_vals)

            if groups is not None:
                n_cases = sum(1 for g in groups.values() if g == "case")
                n_controls = sum(1 for g in groups.values() if g == "control")

                if n_cases > 0 and n_controls > 0:
                    logger.info(
                        f"Auto-detected phenotype from '{field_name}': "
                        f"{n_cases} cases, {n_controls} controls."
                    )
                    return PhenotypeAssignment(
                        groups=groups,
                        n_cases=n_cases,
                        n_controls=n_controls,
                        field_used=field_name,
                        method="auto",
                    )

        # Auto-detection failed
        raise ValueError(
            "Could not auto-detect case/control assignments. "
            "Available metadata fields: "
            f"{list(sample_metadata.keys())}. "
            "Use override_field and override_case_values/override_control_values "
            "to specify the exact field and value mappings. "
            "See Context Transfer: 'GEO Metadata Parsing Is Fragile'."
        )

    @staticmethod
    def _try_assign_groups(
        field_data: dict[str, str],
        sample_ids: list[str],
        unique_vals: set[str],
    ) -> dict[str, str] | None:
        """Try to assign case/control from field values.

        Common patterns:
        - "disease status: autism" vs "disease status: control"
        - "diagnosis: AD" vs "diagnosis: normal"
        - "group: case" vs "group: control"
        """
        # Known control keywords
        control_kw = {
            "control", "normal", "healthy", "unaffected", "baseline",
            "reference", "wt", "wild type", "wildtype", "non-disease",
            "neurotypical", "td",
        }

        groups = {}
        for sid in sample_ids:
            val = field_data.get(sid, "").lower().strip()

            # Extract value after colon if present (e.g., "disease status: autism")
            if ":" in val:
                val_part = val.split(":", 1)[1].strip()
            else:
                val_part = val

            # Check for control keywords
            is_control = any(kw in val_part for kw in control_kw)

            if is_control:
                groups[sid] = "control"
            elif val_part:  # non-empty and not control = case
                groups[sid] = "case"

        if len(groups) == len(sample_ids) and len(set(groups.values())) == 2:
            return groups

        return None


# ---------------------------------------------------------------------------
# Probe-to-Gene Mapper
# ---------------------------------------------------------------------------

# Column name variants for gene symbol in platform annotations
_GENE_SYMBOL_COLUMNS = [
    "Gene Symbol",
    "Gene symbol",
    "gene_symbol",
    "GENE_SYMBOL",
    "Symbol",
    "symbol",
    "gene_assignment",
    "GENE",
    "Gene",
    "ORF",
    "ILMN_Gene",
]

# Column name variants for probe ID
_PROBE_ID_COLUMNS = [
    "ID",
    "ID_REF",
    "SPOT_ID",
    "probe_id",
    "Probe_Id",
    "PROBE_ID",
]


class ProbeGeneMapper:
    """Map probe IDs to gene symbols using platform annotation.

    Handles the many variations in platform annotation file format
    encountered across GEO platforms.

    Parameters
    ----------
    annotation_path : str or Path
        Path to the platform annotation file (TSV/CSV).
    resolver : HGNCResolver or None
        If provided, resolves gene symbols via HGNC.
    gene_column : str or None
        Explicit column name for gene symbols. If None, auto-detects.
    probe_column : str or None
        Explicit column name for probe IDs. If None, auto-detects.
    """

    def __init__(
        self,
        annotation_path: str | Path,
        resolver=None,
        gene_column: str | None = None,
        probe_column: str | None = None,
    ):
        self.annotation_path = Path(annotation_path)
        self.resolver = resolver
        self._probe_to_gene: dict[str, str] = {}
        self._n_multi_gene = 0
        self._gene_column_used = ""

        if not self.annotation_path.exists():
            raise FileNotFoundError(
                f"Platform annotation file not found: {self.annotation_path}"
            )

        self._load(gene_column, probe_column)

    def _load(self, gene_column: str | None, probe_column: str | None) -> None:
        """Load and parse the platform annotation file."""
        # Try tab-separated first, then comma
        try:
            df = pd.read_csv(
                self.annotation_path, sep="\t", dtype=str,
                low_memory=False, comment="#",
            )
        except Exception:
            df = pd.read_csv(
                self.annotation_path, sep=",", dtype=str,
                low_memory=False, comment="#",
            )

        # Find probe ID column
        if probe_column:
            if probe_column not in df.columns:
                raise ValueError(
                    f"Probe column '{probe_column}' not found. "
                    f"Available: {list(df.columns)[:10]}"
                )
            pid_col = probe_column
        else:
            pid_col = self._find_column(df, _PROBE_ID_COLUMNS)
            if pid_col is None:
                # Use first column as fallback
                pid_col = df.columns[0]
                logger.warning(
                    f"No standard probe ID column found. "
                    f"Using first column: '{pid_col}'."
                )

        # Find gene symbol column
        if gene_column:
            if gene_column not in df.columns:
                raise ValueError(
                    f"Gene column '{gene_column}' not found. "
                    f"Available: {list(df.columns)[:10]}"
                )
            sym_col = gene_column
        else:
            sym_col = self._find_column(df, _GENE_SYMBOL_COLUMNS)
            if sym_col is None:
                raise ValueError(
                    f"No gene symbol column found in annotation file. "
                    f"Tried: {_GENE_SYMBOL_COLUMNS}. "
                    f"Available columns: {list(df.columns)}. "
                    f"Use gene_column parameter to specify."
                )

        self._gene_column_used = sym_col

        # Build probe -> gene mapping
        for _, row in df.iterrows():
            probe_id = str(row[pid_col]).strip()
            gene_raw = str(row[sym_col]).strip()

            if not probe_id or probe_id == "nan":
                continue
            if not gene_raw or gene_raw == "nan" or gene_raw == "---":
                continue

            # Handle multi-gene annotations: "GENE1 /// GENE2"
            if "///" in gene_raw:
                symbols = [s.strip() for s in gene_raw.split("///") if s.strip()]
                self._n_multi_gene += 1
            elif "//" in gene_raw:
                # gene_assignment format: "NM_001234 // GENE1 // ..."
                parts = [s.strip() for s in gene_raw.split("//")]
                symbols = [p for p in parts if p and not p.startswith("NM_")
                           and not p.startswith("NR_") and not p.startswith("XM_")
                           and len(p) < 20 and p.upper() == p.replace(" ", "")]
                if not symbols:
                    symbols = [gene_raw]
                self._n_multi_gene += 1
            else:
                symbols = [gene_raw]

            # Use first valid symbol
            for sym in symbols:
                sym = sym.strip()
                if not sym or sym == "---":
                    continue

                # Apply HGNC resolution if available
                if self.resolver is not None:
                    sym = self.resolver.resolve(sym)

                self._probe_to_gene[probe_id] = sym
                break  # use first valid symbol only

    @staticmethod
    def _find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
        """Find the first matching column name from candidates."""
        for col in candidates:
            if col in df.columns:
                return col
        return None

    def get_result(self) -> ProbeMapping:
        """Return the mapping result as a dataclass."""
        return ProbeMapping(
            probe_to_gene=self._probe_to_gene,
            n_mapped=len(self._probe_to_gene),
            n_unmapped=0,  # computed at usage time
            n_multi_gene=self._n_multi_gene,
            gene_column_used=self._gene_column_used,
        )

    def map_expression(
        self,
        expression: pd.DataFrame,
    ) -> pd.DataFrame:
        """Collapse probe-level expression to gene-level.

        For genes with multiple probes, takes the probe with the
        highest mean expression (standard practice).

        Parameters
        ----------
        expression : DataFrame
            Probes × samples expression matrix (index = probe IDs).

        Returns
        -------
        DataFrame
            Genes × samples expression matrix.
        """
        # Map probe IDs to gene symbols
        mapped_probes = []
        for probe_id in expression.index:
            gene = self._probe_to_gene.get(str(probe_id))
            if gene:
                mapped_probes.append((probe_id, gene))

        if not mapped_probes:
            warnings.warn(
                "No probes could be mapped to genes. Check platform "
                "annotation file and column settings.",
                UserWarning,
                stacklevel=2,
            )
            return pd.DataFrame()

        # Create mapped expression with gene column
        probe_ids, genes = zip(*mapped_probes)
        mapped_expr = expression.loc[list(probe_ids)].copy()
        mapped_expr["_gene"] = genes

        # For multi-probe genes, keep the probe with highest mean expression
        mapped_expr["_mean"] = mapped_expr.drop(columns=["_gene"]).mean(axis=1)
        gene_expr = (
            mapped_expr
            .sort_values("_mean", ascending=False)
            .drop_duplicates(subset="_gene", keep="first")
            .drop(columns=["_mean"])
            .set_index("_gene")
        )
        gene_expr.index.name = "gene"

        n_unmapped = len(expression) - len(mapped_probes)
        logger.info(
            f"Probe-to-gene mapping: {len(mapped_probes)} mapped, "
            f"{n_unmapped} unmapped, {len(gene_expr)} unique genes."
        )

        return gene_expr
