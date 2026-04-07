# Riker Engine — Data Dependencies

For the complete list of datasets, platforms, and seed genes per disease, run:

```bash
python scripts/download_data.py --list
```

This shows all 6 validated diseases (ASD, T2D, IBD, AD, breast cancer, IPF),
their GEO dataset accessions, platform annotations, and seed gene sources.

## Quick Start

```bash
# Download data for a specific disease
python scripts/download_data.py ipf

# Download everything
python scripts/download_data.py all
```

## What's included in the repo

- `data/seeds/` — Curated seed gene lists for all 6 diseases (committed)
- `data/mappings/` — Ensembl-to-symbol gene mapping (committed)

## What's downloaded by the script

- `data/geo/<disease>/` — GEO series matrix files (gitignored, ~5-80 MB each)
- `data/platforms/` — Platform annotation files (gitignored, ~15-40 MB each)
- `data/hgnc/` — HGNC complete gene set (gitignored, ~17 MB)

## RNA-seq datasets

Two ASD datasets (GSE64018, GSE102741) and one T2D dataset (GSE86468) are
RNA-seq and require manual reconstruction from supplementary data.
See `docs/DATA_RECONSTRUCTION.md` for instructions.

All other datasets (IBD, AD, breast cancer, IPF) are microarray and
download automatically.
