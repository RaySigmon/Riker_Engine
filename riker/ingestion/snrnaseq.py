"""
Riker Engine - snRNA-seq Pseudo-Bulking.

Converts single-nucleus RNA-seq data into pseudo-bulk expression
matrices suitable for the Riker pipeline. Aggregates nuclei within
each donor (or donor × cell type) to produce genes × samples matrices.

Supports multiple input formats:
- AnnData (.h5ad) via scanpy (if available)
- CSV/TSV count matrices with metadata
- Pre-computed pseudo-bulk matrices (passthrough)

References:
    Blueprint Section 5.5 (snRNA-seq Pseudo-Bulking)
"""

import logging
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Try to import scanpy for h5ad support
try:
    import scanpy as sc
    _HAS_SCANPY = True
except ImportError:
    _HAS_SCANPY = False


@dataclass(frozen=True)
class PseudoBulkResult:
    """Result of pseudo-bulking snRNA-seq data.

    Attributes:
        expression: DataFrame (genes × samples) of pseudo-bulk expression.
        sample_metadata: Dict of sample_id -> metadata dict.
        phenotypes: Dict of sample_id -> 'case' or 'control'.
        cell_type: Cell type used for pseudo-bulking (or 'all').
        n_donors: Number of unique donors.
        n_nuclei_total: Total nuclei before aggregation.
        n_nuclei_per_sample: Dict of sample_id -> nucleus count.
        min_nuclei_per_sample: Minimum nuclei in any sample.
        genes_detected: Number of genes with nonzero expression.
    """
    expression: pd.DataFrame
    sample_metadata: dict
    phenotypes: dict
    cell_type: str
    n_donors: int
    n_nuclei_total: int
    n_nuclei_per_sample: dict
    min_nuclei_per_sample: int
    genes_detected: int


def pseudo_bulk_from_counts(
    counts: pd.DataFrame,
    donor_col: str,
    condition_col: str,
    cell_type_col: str | None = None,
    target_cell_type: str | None = None,
    case_values: list[str] | None = None,
    control_values: list[str] | None = None,
    min_nuclei: int = 10,
    min_genes: int = 200,
    aggregation: str = "sum",
    normalize: bool = True,
) -> PseudoBulkResult:
    """Pseudo-bulk snRNA-seq count data by donor.

    Parameters
    ----------
    counts : DataFrame
        Nuclei × genes count matrix. Index = nucleus barcodes,
        columns = gene symbols. Can also include metadata columns.
    donor_col : str
        Column name for donor/sample identity.
    condition_col : str
        Column name for case/control condition.
    cell_type_col : str or None
        Column name for cell type annotations. If None, use all nuclei.
    target_cell_type : str or None
        Specific cell type to pseudo-bulk. If None and cell_type_col
        is provided, pseudo-bulk all cell types together.
    case_values : list or None
        Values in condition_col indicating case (e.g., ['ASD', 'AD']).
        If None, auto-detect.
    control_values : list or None
        Values in condition_col indicating control.
        If None, auto-detect.
    min_nuclei : int
        Minimum nuclei per donor to include (default 10).
    min_genes : int
        Minimum detected genes per nucleus for QC (default 200).
    aggregation : str
        'sum' (recommended for DE) or 'mean'.
    normalize : bool
        If True, apply CPM normalization then log2(CPM + 1).

    Returns
    -------
    PseudoBulkResult

    Raises
    ------
    ValueError
        If required columns are missing or no valid samples remain.
    """
    # Validate required columns exist
    if donor_col not in counts.columns:
        raise ValueError(f"Donor column '{donor_col}' not found in data.")
    if condition_col not in counts.columns:
        raise ValueError(f"Condition column '{condition_col}' not found in data.")

    # Separate metadata from gene counts
    meta_cols = [donor_col, condition_col]
    if cell_type_col and cell_type_col in counts.columns:
        meta_cols.append(cell_type_col)

    gene_cols = [c for c in counts.columns if c not in meta_cols]
    metadata = counts[meta_cols].copy()
    gene_counts = counts[gene_cols].copy()

    # Ensure numeric
    gene_counts = gene_counts.apply(pd.to_numeric, errors="coerce").fillna(0)

    n_nuclei_total = len(counts)
    logger.info(f"Input: {n_nuclei_total} nuclei, {len(gene_cols)} genes.")

    # Filter by cell type if specified
    if cell_type_col and target_cell_type:
        if cell_type_col not in metadata.columns:
            raise ValueError(f"Cell type column '{cell_type_col}' not found.")
        mask = metadata[cell_type_col].str.lower() == target_cell_type.lower()
        metadata = metadata[mask]
        gene_counts = gene_counts.loc[mask.values]
        cell_type_label = target_cell_type
        logger.info(
            f"Filtered to {len(gene_counts)} nuclei "
            f"of type '{target_cell_type}'."
        )
    else:
        cell_type_label = "all"

    if len(gene_counts) == 0:
        raise ValueError("No nuclei remaining after cell type filter.")

    # Basic QC: filter nuclei with too few detected genes
    genes_per_nucleus = (gene_counts > 0).sum(axis=1)
    qc_mask = genes_per_nucleus >= min_genes
    n_filtered = (~qc_mask).sum()
    if n_filtered > 0:
        logger.info(
            f"QC: removed {n_filtered} nuclei with < {min_genes} genes detected."
        )
        gene_counts = gene_counts[qc_mask]
        metadata = metadata[qc_mask]

    if len(gene_counts) == 0:
        raise ValueError("No nuclei remaining after QC filter.")

    # Auto-detect case/control if not provided
    if case_values is None or control_values is None:
        unique_conditions = metadata[condition_col].unique()
        if len(unique_conditions) != 2:
            raise ValueError(
                f"Expected 2 conditions, found {len(unique_conditions)}: "
                f"{unique_conditions}. Please specify case_values and "
                f"control_values explicitly."
            )
        # Try to detect which is case vs control
        case_values, control_values = _detect_case_control(unique_conditions)

    # Assign case/control
    conditions = metadata[condition_col].astype(str).str.lower()
    case_mask = conditions.isin([v.lower() for v in case_values])
    ctrl_mask = conditions.isin([v.lower() for v in control_values])

    if case_mask.sum() == 0:
        raise ValueError(f"No nuclei matched case values: {case_values}")
    if ctrl_mask.sum() == 0:
        raise ValueError(f"No nuclei matched control values: {control_values}")

    # Aggregate by donor
    donors = metadata[donor_col].values
    unique_donors = sorted(set(donors))

    pseudo_bulk_data = {}
    sample_metadata = {}
    phenotypes = {}
    n_nuclei_per_sample = {}

    for donor in unique_donors:
        donor_mask = donors == donor
        donor_counts = gene_counts[donor_mask]
        n_nuc = len(donor_counts)

        if n_nuc < min_nuclei:
            logger.warning(
                f"Donor {donor}: {n_nuc} nuclei < {min_nuclei} minimum. Skipping."
            )
            continue

        # Aggregate
        if aggregation == "sum":
            agg = donor_counts.sum(axis=0).values.astype(np.float64)
        elif aggregation == "mean":
            agg = donor_counts.mean(axis=0).values.astype(np.float64)
        else:
            raise ValueError(f"Unknown aggregation method: {aggregation}")

        # Determine condition for this donor
        donor_conditions = metadata.loc[donor_mask, condition_col]
        donor_condition = donor_conditions.mode().iloc[0]
        condition_lower = str(donor_condition).lower()

        if condition_lower in [v.lower() for v in case_values]:
            group = "case"
        elif condition_lower in [v.lower() for v in control_values]:
            group = "control"
        else:
            logger.warning(f"Donor {donor}: unknown condition '{donor_condition}'. Skipping.")
            continue

        sample_id = f"{donor}_{cell_type_label}"
        pseudo_bulk_data[sample_id] = agg
        phenotypes[sample_id] = group
        n_nuclei_per_sample[sample_id] = n_nuc
        sample_metadata[sample_id] = {
            "donor": donor,
            "condition": str(donor_condition),
            "cell_type": cell_type_label,
            "n_nuclei": n_nuc,
            "group": group,
        }

    if not pseudo_bulk_data:
        raise ValueError("No valid pseudo-bulk samples produced.")

    # Build expression matrix (genes × samples)
    expr_df = pd.DataFrame(pseudo_bulk_data, index=gene_cols)
    expr_df.index.name = "gene"

    # Normalize if requested
    if normalize:
        # CPM normalization per sample
        col_sums = expr_df.sum(axis=0)
        col_sums = col_sums.replace(0, 1)  # avoid division by zero
        expr_df = expr_df.div(col_sums, axis=1) * 1e6
        # Log2(CPM + 1)
        expr_df = np.log2(expr_df + 1)
        logger.info("Applied CPM + log2(CPM+1) normalization.")

    # Filter genes with zero expression across all samples
    nonzero_mask = expr_df.sum(axis=1) > 0
    genes_detected = nonzero_mask.sum()
    expr_df = expr_df[nonzero_mask]
    logger.info(f"Genes with nonzero expression: {genes_detected}")

    min_nuc = min(n_nuclei_per_sample.values()) if n_nuclei_per_sample else 0

    logger.info(
        f"Pseudo-bulk complete: {len(pseudo_bulk_data)} samples "
        f"({sum(1 for v in phenotypes.values() if v == 'case')} cases, "
        f"{sum(1 for v in phenotypes.values() if v == 'control')} controls), "
        f"{genes_detected} genes."
    )

    return PseudoBulkResult(
        expression=expr_df,
        sample_metadata=sample_metadata,
        phenotypes=phenotypes,
        cell_type=cell_type_label,
        n_donors=len(pseudo_bulk_data),
        n_nuclei_total=n_nuclei_total,
        n_nuclei_per_sample=n_nuclei_per_sample,
        min_nuclei_per_sample=min_nuc,
        genes_detected=int(genes_detected),
    )


def pseudo_bulk_from_h5ad(
    h5ad_path: str | Path,
    donor_key: str = "donor_id",
    condition_key: str = "disease",
    cell_type_key: str = "cell_type",
    target_cell_type: str | None = None,
    case_values: list[str] | None = None,
    control_values: list[str] | None = None,
    min_nuclei: int = 10,
    min_genes: int = 200,
    layer: str | None = None,
) -> PseudoBulkResult:
    """Pseudo-bulk from an AnnData h5ad file.

    Parameters
    ----------
    h5ad_path : str or Path
        Path to .h5ad file (AnnData format from scanpy).
    donor_key : str
        Key in adata.obs for donor identity.
    condition_key : str
        Key in adata.obs for case/control condition.
    cell_type_key : str
        Key in adata.obs for cell type annotations.
    target_cell_type : str or None
        Specific cell type to pseudo-bulk.
    case_values, control_values : list or None
        Condition values for case/control groups.
    min_nuclei : int
        Minimum nuclei per donor.
    min_genes : int
        Minimum detected genes per nucleus.
    layer : str or None
        AnnData layer to use. None = adata.X.

    Returns
    -------
    PseudoBulkResult
    """
    if not _HAS_SCANPY:
        raise ImportError(
            "scanpy is required for h5ad support. "
            "Install with: pip install scanpy"
        )

    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        raise FileNotFoundError(f"h5ad file not found: {h5ad_path}")

    logger.info(f"Loading h5ad: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)

    # Extract count matrix
    if layer and layer in adata.layers:
        X = adata.layers[layer]
    else:
        X = adata.X

    # Convert sparse to dense if needed
    if hasattr(X, "toarray"):
        X = X.toarray()

    # Build a DataFrame with metadata columns
    gene_names = list(adata.var_names)
    counts_df = pd.DataFrame(X, columns=gene_names, index=adata.obs_names)

    # Add metadata columns
    counts_df[donor_key] = adata.obs[donor_key].values
    counts_df[condition_key] = adata.obs[condition_key].values
    if cell_type_key in adata.obs.columns:
        counts_df[cell_type_key] = adata.obs[cell_type_key].values
        ct_col = cell_type_key
    else:
        ct_col = None

    return pseudo_bulk_from_counts(
        counts=counts_df,
        donor_col=donor_key,
        condition_col=condition_key,
        cell_type_col=ct_col,
        target_cell_type=target_cell_type,
        case_values=case_values,
        control_values=control_values,
        min_nuclei=min_nuclei,
        min_genes=min_genes,
        aggregation="sum",
        normalize=True,
    )


def _detect_case_control(conditions) -> tuple[list[str], list[str]]:
    """Auto-detect which condition values are case vs control.

    Uses keyword matching similar to PhenotypeExtractor.
    """
    control_keywords = [
        "control", "ctrl", "normal", "healthy", "unaffected",
        "neurotypical", "nt", "wt", "wild"
    ]

    case_vals = []
    ctrl_vals = []

    for cond in conditions:
        cond_lower = str(cond).lower()
        if any(kw in cond_lower for kw in control_keywords):
            ctrl_vals.append(str(cond))
        else:
            case_vals.append(str(cond))

    if not case_vals or not ctrl_vals:
        raise ValueError(
            f"Could not auto-detect case/control from conditions: "
            f"{list(conditions)}. Please specify explicitly."
        )

    return case_vals, ctrl_vals
