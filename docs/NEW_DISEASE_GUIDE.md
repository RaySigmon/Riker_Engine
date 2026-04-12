# Adding a New Disease to the Riker Engine

This guide walks through the full process of configuring the Riker Engine for a disease it has never seen. You will need:

- A list of candidate genes (seed genes) for your disease
- 3+ GEO microarray datasets with case/control samples
- The platform annotation files for those datasets

**Time estimate:** 15–30 minutes for an experienced user, longer for your first attempt.

---

## Seed List Guidelines

The seed gene list is your starting hypothesis — the candidates the engine will validate across datasets. Seed list size affects pipeline behavior:

- **200–500 genes** is the sweet spot for most diseases. This gives consensus clustering enough signal to find meaningful modules while keeping the analysis focused.
- **50–100 genes** works but may produce fewer clusters and less granular pathway separation. Use this range when your candidate list is naturally small (e.g., a single GWAS locus).
- **500+ genes** works well. More seeds give the pipeline more signal for cross-referencing. The engine filters aggressively — it's better to be inclusive and let the pipeline do the filtering than to be too conservative upfront.

**Where to get seed genes:**
- [Open Targets](https://platform.opentargets.org/) — broad disease-gene associations (recommended starting point)
- GWAS Catalog — genetic associations
- [DisGeNET](https://www.disgenet.org/) — curated disease-gene databases
- COSMIC — for cancer-specific genes
- Published differential expression gene lists from individual studies
- KEGG pathway members relevant to your disease biology

The validated diseases in this repository used seed lists ranging from 354 (IPF) to 1,267 (ASD), all producing robust results.

---

## Step 1: Find GEO Datasets

Search [GEO DataSets](https://www.ncbi.nlm.nih.gov/gds/) for your disease. Look for:

- **Microarray** studies (not RNA-seq — those require reconstruction; see `docs/DATA_RECONSTRUCTION.md`)
- **Case/control** design (not time-series or drug-response)
- **Standard Affymetrix or Illumina platforms** (GPL570, GPL6244, GPL96, GPL6480, GPL10558 are well-tested)
- **Minimum 5 samples per group** (more is better)

You need at least 3 discovery datasets and ideally 2+ replication datasets. Use different studies (not subseries of the same experiment) for true independence.

---

## Step 2: Identify the Platform

Every GEO series matrix encodes which microarray platform was used. To check:

```bash
zcat GSExxxxx_series_matrix.txt.gz | grep "!Series_platform_id"
```

This returns something like `"GPL570"`. Download the platform annotation:

```bash
python scripts/download_data.py --platforms
# Or download a specific one manually from GEO
```

**Critical:** Keep discovery and replication datasets on the same platform family when possible (e.g., all Affymetrix HG-U133). Mixing platforms can cause gene symbol mismatches at Phase 5 that look like replication failures but are actually probe resolution differences.

---

## Step 3: Inspect Phenotype Labels

This is the hardest part. GEO metadata is not standardized — every dataset encodes disease status differently. You need to find the exact field and values that distinguish cases from controls.

### Find the phenotype field

```bash
zcat GSExxxxx_series_matrix.txt.gz | grep "^!Sample_characteristics_ch1"
zcat GSExxxxx_series_matrix.txt.gz | grep "^!Sample_source_name_ch1"
zcat GSExxxxx_series_matrix.txt.gz | grep "^!Sample_title"
```

Look for lines containing disease-related terms (diagnosis, disease state, group, phenotype, status).

### Understand multi-line characteristics

`Sample_characteristics_ch1` can span multiple lines in a series matrix. The parser joins them with ` | ` and uses **substring matching**. So if the metadata looks like:

```
!Sample_characteristics_ch1  "gender: female"  "gender: male"  ...
!Sample_characteristics_ch1  "disease: IPF"    "disease: control"  ...
```

The parser sees each sample as `"gender: female | disease: IPF"`. Your `case_values` and `control_values` match as substrings within that combined string. So `case_values: ["disease: ipf"]` will match correctly (matching is case-insensitive).

### Verify your values match

Before running the full pipeline, confirm your case/control values actually match:

```bash
zcat GSExxxxx_series_matrix.txt.gz | grep "^!Sample_characteristics_ch1" | \
  tr '\t' '\n' | sort | uniq -c | sort -rn | head -20
```

This shows you the unique values and their frequencies. Your `case_values` and `control_values` must be substrings of these exact values.

---

## Step 4: Verify Probe Mapping

Before running the full pipeline, confirm the platform annotation maps probes to genes correctly:

```python
from riker.ingestion.geo_parser import ProbeGeneMapper

mapper = ProbeGeneMapper("data/platforms/GPLxxxx.annot")
result = mapper.get_result()
print(f"Probes mapped: {result.n_mapped}")
print(f"Gene column used: {result.gene_column_used}")

# Spot-check a few mappings
items = list(result.probe_to_gene.items())[:5]
for probe, gene in items:
    print(f"  {probe} -> {gene}")
```

If `n_mapped` is 0 or very low, the annotation file may use non-standard column names. Use `gene_column` and `probe_column` in your YAML config to override auto-detection:

```yaml
datasets:
  - id: GSExxxxx
    platform: data/platforms/GPLxxxx.annot
    gene_column: "Gene Symbol"    # exact column name from the annotation file
    probe_column: "ID"            # exact column name for probe IDs
```

---

## Step 5: Prepare Seed Genes

Create a CSV file with one column of HGNC gene symbols:

```csv
symbol
BRCA1
TP53
ERBB2
```

Sources for seed genes:
- [Open Targets](https://platform.opentargets.org/) — search your disease, export associated genes
- [SFARI Gene](https://gene.sfari.org/) — for autism
- [DisGeNET](https://www.disgenet.org/) — broad disease-gene associations
- Literature curation — compile from published GWAS/transcriptomic studies

Place the file in `data/seeds/your_disease_genes.csv`.

---

## Step 6: Write the Config YAML

```yaml
condition: your_disease
seed_genes: data/seeds/your_disease_genes.csv
hgnc_path: data/hgnc/hgnc_complete_set.txt
output_dir: output/your_disease

datasets:
  # Discovery (minimum 3)
  - id: GSE11111
    series_matrix: data/geo/your_disease/GSE11111_series_matrix.txt.gz
    platform: data/platforms/GPL570.annot
    role: discovery
    tissue: your_tissue
    phenotype_field: Sample_characteristics_ch1
    case_values: ["disease state: your disease"]
    control_values: ["disease state: control"]

  - id: GSE22222
    series_matrix: data/geo/your_disease/GSE22222_series_matrix.txt.gz
    platform: data/platforms/GPL570.annot
    role: discovery
    tissue: your_tissue
    phenotype_field: Sample_characteristics_ch1
    case_values: ["diagnosis: your disease"]
    control_values: ["diagnosis: healthy"]

  - id: GSE33333
    series_matrix: data/geo/your_disease/GSE33333_series_matrix.txt.gz
    platform: data/platforms/GPL570.annot
    role: discovery
    tissue: your_tissue
    phenotype_field: Sample_source_name_ch1
    case_values: ["disease tissue"]
    control_values: ["normal tissue"]

  # Replication (minimum 2)
  - id: GSE44444
    series_matrix: data/geo/your_disease/GSE44444_series_matrix.txt.gz
    platform: data/platforms/GPL570.annot
    role: replication
    tissue: your_tissue
    phenotype_field: Sample_characteristics_ch1
    case_values: ["disease: your disease"]
    control_values: ["disease: control"]
```

Note: every dataset can use a different `phenotype_field` and different `case_values`/`control_values`. This is by design — GEO metadata is not standardized.

---

## Step 7: Run the Pipeline

```bash
# Download HGNC data if you haven't already
python scripts/download_data.py --hgnc

# Run
riker run configs/your_disease_config.yaml
```

Check `output/your_disease/qc_report.json` first. If the pipeline halts, the QC report explains why.

---

## Common Pitfalls

### Exon array probe IDs

Some platforms (e.g., GPL6244 HuGene-1.0-ST) use numeric probe IDs like `7892501` while others (e.g., GPL570 HG-U133 Plus 2.0) use Affymetrix-style IDs like `1007_s_at`. Both formats work fine — the probe-to-gene mapper handles both. But do not use an annotation file from one platform with expression data from another.

### Disease labels are not what you expect

GEO submitters use inconsistent terminology. "Multiple sclerosis" might appear as:
- `disease: MS`
- `diagnosis: CIS` (clinically isolated syndrome — early MS)
- `disease state: relapsing remitting multiple sclerosis`
- Buried in `Sample_title` instead of `Sample_characteristics_ch1`

Always inspect the raw metadata before writing your config.

### Too few genes pass Phase 1

If fewer than ~20 genes pass Phase 1 (cross-referencing), the clustering in Phase 3 won't have enough signal. Possible causes:
- Seed gene list too small or too noisy
- Datasets have very small sample sizes (n < 5 per group)
- Case/control labels are wrong (swapped or misassigned)
- Platform mismatch (expression data and annotation file from different platforms)

### Phase 5 fails to find core genes

If replication datasets can't find core genes, check:
1. Are the replication datasets from the same tissue as discovery? Cross-tissue non-replication is expected and tolerated.
2. Is the platform annotation mapping correctly? Verify with the probe mapping check in Step 4.
3. Are the case/control labels correct for the replication datasets?

### Quoted vs unquoted values

The parser strips quotes from GEO metadata. Use the unquoted form in your config:
- Correct: `case_values: ["disease: IPF"]`
- Wrong: `case_values: ["\"disease: IPF\""]`

---

## Worked Example: Parkinson's Disease

Here's a complete walkthrough for adding Parkinson's Disease.

### 1. Find datasets on GEO

Search GEO for "Parkinson's disease" microarray case-control studies in substantia nigra or brain tissue. Good candidates:

- **GSE20141** — Substantia nigra, GPL570, 10 PD vs 8 controls
- **GSE20163** — Substantia nigra, GPL96, 8 PD vs 9 controls
- **GSE20164** — Substantia nigra, GPL96, 6 PD vs 5 controls
- **GSE49036** — Frontal cortex + cerebellum, GPL6244, 16 PD vs 12 controls
- **GSE8397** — Substantia nigra + frontal cortex, GPL96, 15 PD vs 12 controls

### 2. Inspect metadata

```bash
zcat GSE20141_series_matrix.txt.gz | grep "^!Sample_characteristics_ch1"
# Output: "disease state: Parkinson's disease"  "disease state: control"
```

```bash
zcat GSE20163_series_matrix.txt.gz | grep "^!Sample_source_name_ch1"
# May need to check which field has the labels
```

### 3. Get seed genes

Download Parkinson's Disease associated genes from [Open Targets](https://platform.opentargets.org/disease/MONDO_0005180/associations). Export the gene list, filter to high-confidence associations, save as `data/seeds/pd_open_targets_genes.csv`.

### 4. Write the config

```yaml
condition: PD
seed_genes: data/seeds/pd_open_targets_genes.csv
hgnc_path: data/hgnc/hgnc_complete_set.txt
output_dir: output/pd_curated

datasets:
  - id: GSE20141
    series_matrix: data/geo/pd/GSE20141_series_matrix.txt.gz
    platform: data/platforms/GPL570.annot
    role: discovery
    tissue: brain
    phenotype_field: Sample_characteristics_ch1
    case_values: ["disease state: parkinson"]
    control_values: ["disease state: control"]

  # ... additional datasets following the same pattern
```

### 5. Run

```bash
python scripts/download_data.py --hgnc
riker run configs/pd_curated.yaml
```

The pipeline handles everything from here — cross-referencing, clustering, robustness testing, replication, and meta-analysis. Check `output/pd_curated/pipeline_summary.json` for final results.
