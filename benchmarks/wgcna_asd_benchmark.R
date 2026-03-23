#!/usr/bin/env Rscript
# WGCNA benchmark on ASD brain cortex discovery datasets
# Matches Riker Engine's discovery data exactly: GSE28521, GSE28475, GSE64018

.libPaths("~/R/library")

suppressPackageStartupMessages({
  library(WGCNA)
})

# Allow multi-threading on Pi 5
allowWGCNAThreads(nThreads = 4)

options(stringsAsFactors = FALSE)

cat("=== WGCNA ASD Benchmark ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# --- Helper functions ---

parse_series_matrix <- function(gz_path) {
  # Read all lines
  con <- gzfile(gz_path, "rt")
  lines <- readLines(con)
  close(con)

  # Extract metadata lines
  meta_lines <- lines[grepl("^!", lines)]

  # Find data start (first line that starts with "ID_REF")
  data_start <- which(grepl("^\"ID_REF\"", lines))
  # Find data end (line with "!series_matrix_table_end" or end of file)
  data_end_marker <- which(grepl("^!series_matrix_table_end", lines))
  if (length(data_end_marker) > 0) {
    data_end <- data_end_marker[1] - 1
  } else {
    data_end <- length(lines)
  }

  # Parse expression data
  data_lines <- lines[data_start:data_end]
  # Write to temp file for fast reading
  tmp <- tempfile(fileext = ".tsv")
  writeLines(data_lines, tmp)
  expr <- read.delim(tmp, header = TRUE, row.names = 1, check.names = FALSE,
                     quote = "\"")
  unlink(tmp)

  # Parse sample characteristics
  char_lines <- meta_lines[grepl("^!Sample_characteristics_ch1", meta_lines)]
  # Parse sample geo accessions for column matching
  geo_line <- meta_lines[grepl("^!Sample_geo_accession", meta_lines)][1]
  geo_ids <- unlist(strsplit(geo_line, "\t"))[-1]
  geo_ids <- gsub("\"", "", geo_ids)

  list(expr = expr, char_lines = char_lines, geo_ids = geo_ids, meta = meta_lines)
}

extract_phenotypes <- function(char_lines, geo_ids, case_pattern, control_pattern) {
  # Each char_line is tab-separated: field_name \t val1 \t val2 ...
  # Find the line that contains case/control info
  pheno <- rep(NA, length(geo_ids))

  for (cl in char_lines) {
    parts <- unlist(strsplit(cl, "\t"))
    vals <- gsub("\"", "", parts[-1])

    case_match <- grepl(case_pattern, vals, ignore.case = TRUE)
    ctrl_match <- grepl(control_pattern, vals, ignore.case = TRUE)

    if (any(case_match) && any(ctrl_match)) {
      pheno[case_match] <- "ASD"
      pheno[ctrl_match] <- "Control"
      break
    }
  }

  names(pheno) <- geo_ids
  pheno
}

map_probes_to_genes <- function(expr, annot_path) {
  annot <- read.delim(annot_path, header = TRUE, check.names = FALSE,
                      comment.char = "#")

  # Find gene symbol column
  gene_col <- NULL
  for (candidate in c("Gene symbol", "Gene Symbol", "Gene.symbol", "Gene.Symbol")) {
    if (candidate %in% colnames(annot)) {
      gene_col <- candidate
      break
    }
  }
  if (is.null(gene_col)) {
    stop("Cannot find gene symbol column in annotation file: ", annot_path)
  }

  id_col <- colnames(annot)[1]

  # Filter annotation to probes in expression data
  annot_sub <- annot[annot[[id_col]] %in% rownames(expr), ]
  annot_sub <- annot_sub[annot_sub[[gene_col]] != "" & !is.na(annot_sub[[gene_col]]), ]

  # Some probes map to multiple genes (///), take first
  annot_sub$gene <- sapply(strsplit(as.character(annot_sub[[gene_col]]), "///"), function(x) trimws(x[1]))
  annot_sub <- annot_sub[annot_sub$gene != "", ]

  # Merge
  expr_annot <- expr[as.character(annot_sub[[id_col]]), , drop = FALSE]
  expr_annot$gene <- annot_sub$gene

  # Collapse probes to genes by taking the one with highest variance
  genes <- unique(expr_annot$gene)
  gene_expr <- matrix(NA, nrow = length(genes), ncol = ncol(expr))
  rownames(gene_expr) <- genes
  colnames(gene_expr) <- colnames(expr)

  for (i in seq_along(genes)) {
    g <- genes[i]
    rows <- which(expr_annot$gene == g)
    if (length(rows) == 1) {
      gene_expr[i, ] <- as.numeric(expr_annot[rows, 1:ncol(expr)])
    } else {
      sub <- as.matrix(expr_annot[rows, 1:ncol(expr)])
      vars <- apply(sub, 1, var, na.rm = TRUE)
      gene_expr[i, ] <- as.numeric(sub[which.max(vars), ])
    }
  }

  as.data.frame(gene_expr)
}


run_wgcna <- function(expr_genes, phenotypes, dataset_id, out_dir) {
  cat("\n========================================\n")
  cat("Processing", dataset_id, "\n")
  cat("Expression matrix:", nrow(expr_genes), "genes x", ncol(expr_genes), "samples\n")
  cat("Phenotype counts:", table(phenotypes), "\n")
  cat("========================================\n")

  # Transpose for WGCNA: samples in rows, genes in columns
  datExpr <- as.data.frame(t(expr_genes))

  # Good samples/genes check
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    cat("Removing", sum(!gsg$goodGenes), "genes and", sum(!gsg$goodSamples), "samples\n")
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
    # Update phenotypes to match
    phenotypes <- phenotypes[gsg$goodSamples]
  }

  cat("After QC:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")

  # --- Pick soft threshold ---
  cat("Picking soft threshold power...\n")
  powers <- c(1:10, seq(12, 20, 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 2,
                           networkType = "signed")

  # Pick power: first to reach R^2 >= 0.80, or if none, use the one with highest R^2
  r2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  above_threshold <- which(r2 > 0.80)
  if (length(above_threshold) > 0) {
    power <- powers[above_threshold[1]]
  } else {
    power <- powers[which.max(r2)]
    cat("WARNING: No power reached R^2 > 0.80. Using best available:", power,
        "with R^2 =", max(r2), "\n")
  }
  cat("Selected soft threshold power:", power, "\n")

  # Save soft threshold data
  write.csv(sft$fitIndices, file.path(out_dir, paste0(dataset_id, "_sft.csv")),
            row.names = FALSE)

  # --- Build network and detect modules ---
  cat("Building TOM (this may take a while)...\n")
  t_start <- Sys.time()

  net <- blockwiseModules(
    datExpr,
    power = power,
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 3,
    maxBlockSize = ncol(datExpr)  # single block
  )

  t_end <- Sys.time()
  cat("TOM + clustering took:", format(t_end - t_start), "\n")

  moduleLabels <- net$colors
  moduleColors <- labels2colors(moduleLabels)
  n_modules <- length(unique(moduleLabels)) - 1  # exclude grey/0
  cat("Found", n_modules, "modules (+1 grey/unassigned)\n")

  # --- Module-trait correlation ---
  cat("Computing module-trait correlations...\n")
  trait <- ifelse(phenotypes == "ASD", 1, 0)
  names(trait) <- names(phenotypes)
  # Keep only samples in datExpr
  trait <- trait[rownames(datExpr)]

  MEs <- net$MEs
  # Ensure ordered by module number
  MEs <- orderMEs(MEs)

  # Correlation of each ME with trait
  modTraitCor <- cor(MEs, trait, use = "p")
  modTraitPval <- corPvalueStudent(modTraitCor, nrow(datExpr))

  module_results <- data.frame(
    module = colnames(MEs),
    correlation = as.numeric(modTraitCor),
    pvalue = as.numeric(modTraitPval),
    n_genes = as.numeric(table(moduleLabels)[as.character(
      as.integer(gsub("ME", "", colnames(MEs)))
    )])
  )
  module_results$significant <- module_results$pvalue < 0.05
  module_results$fdr <- p.adjust(module_results$pvalue, method = "BH")

  cat("\nModule-trait correlations:\n")
  print(module_results[order(module_results$pvalue), ], row.names = FALSE)

  # --- Gene-module assignments ---
  gene_modules <- data.frame(
    gene = colnames(datExpr),
    module_number = moduleLabels,
    module_color = moduleColors
  )

  # Save results
  write.csv(gene_modules,
            file.path(out_dir, paste0(dataset_id, "_gene_modules.csv")),
            row.names = FALSE)
  write.csv(module_results,
            file.path(out_dir, paste0(dataset_id, "_module_trait.csv")),
            row.names = FALSE)

  cat("Results saved for", dataset_id, "\n")

  list(
    gene_modules = gene_modules,
    module_results = module_results,
    power = power,
    n_modules = n_modules
  )
}


# ============================================================
# DATASET 1: GSE28521 — Voineagu et al., brain cortex microarray
# ============================================================
cat("\n--- Loading GSE28521 ---\n")
d1 <- parse_series_matrix("/home/kai001/asd_full_run/data/bulk_geo/GSE28521_series_matrix.txt.gz")
p1 <- extract_phenotypes(d1$char_lines, d1$geo_ids,
                         "disease status: autism", "disease status: control")
cat("Phenotypes:", table(p1), "\n")

# Map probes to genes
cat("Mapping probes to gene symbols...\n")
expr1 <- map_probes_to_genes(d1$expr,
                              "/home/kai001/asd_full_run/data/bulk_geo/GPL6883_clean.annot")
cat("Mapped to", nrow(expr1), "unique genes\n")

# Filter to samples with known phenotype
keep1 <- !is.na(p1) & (names(p1) %in% colnames(expr1))
expr1 <- expr1[, names(p1)[keep1]]
p1 <- p1[keep1]

out_dir <- "/home/kai001/riker-engine/benchmarks/results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

res1 <- run_wgcna(expr1, p1, "GSE28521", out_dir)


# ============================================================
# DATASET 2: GSE28475 — brain cortex microarray (GPL6883)
# ============================================================
cat("\n--- Loading GSE28475 ---\n")
d2 <- parse_series_matrix("/home/kai001/asd_full_run/data/bulk_geo/GSE28475-GPL6883_fixed_series_matrix.txt.gz")
p2 <- extract_phenotypes(d2$char_lines, d2$geo_ids,
                         "diagnosis: autism", "diagnosis: control")
cat("Phenotypes:", table(p2), "\n")

expr2 <- map_probes_to_genes(d2$expr,
                              "/home/kai001/asd_full_run/data/bulk_geo/GPL6883_clean.annot")
cat("Mapped to", nrow(expr2), "unique genes\n")

keep2 <- !is.na(p2) & (names(p2) %in% colnames(expr2))
expr2 <- expr2[, names(p2)[keep2]]
p2 <- p2[keep2]

res2 <- run_wgcna(expr2, p2, "GSE28475", out_dir)


# ============================================================
# DATASET 3: GSE64018 — Gupta et al., brain cortex RNA-seq
# Already gene-level (ENSG IDs), need to map to symbols
# ============================================================
cat("\n--- Loading GSE64018 ---\n")
d3 <- parse_series_matrix("/home/kai001/asd_full_run/data/bulk_geo/GSE64018_reconstructed_series_matrix.txt.gz")
p3 <- extract_phenotypes(d3$char_lines, d3$geo_ids,
                         "diagnosis: asd", "diagnosis: ctl")
cat("Phenotypes:", table(p3), "\n")

# Map ENSG to gene symbols
annot3 <- read.delim("/home/kai001/asd_full_run/data/bulk_geo/GPL11154_ensembl.annot",
                     header = TRUE, check.names = FALSE)
# annot3 has columns: ID, Gene Symbol
ensg_to_gene <- annot3[annot3[["Gene Symbol"]] != "" & !is.na(annot3[["Gene Symbol"]]), ]
rownames(ensg_to_gene) <- ensg_to_gene$ID

# Filter expression to annotated genes
common_ids <- intersect(rownames(d3$expr), ensg_to_gene$ID)
expr3 <- d3$expr[common_ids, ]
rownames(expr3) <- ensg_to_gene[common_ids, "Gene Symbol"]

# Remove duplicate gene symbols (keep highest variance)
if (any(duplicated(rownames(expr3)))) {
  vars <- apply(expr3, 1, var, na.rm = TRUE)
  expr3$gene <- rownames(expr3)
  expr3$var <- vars
  expr3 <- expr3[order(-expr3$var), ]
  expr3 <- expr3[!duplicated(expr3$gene), ]
  rownames(expr3) <- expr3$gene
  expr3$gene <- NULL
  expr3$var <- NULL
}

cat("Mapped to", nrow(expr3), "unique genes\n")

keep3 <- !is.na(p3) & (names(p3) %in% colnames(expr3))
expr3 <- expr3[, names(p3)[keep3]]
p3 <- p3[keep3]

res3 <- run_wgcna(expr3, p3, "GSE64018", out_dir)


# ============================================================
# Summary
# ============================================================
cat("\n\n=== WGCNA BENCHMARK SUMMARY ===\n")
cat("Dataset       | Power | Modules | Sig Modules (p<0.05)\n")
cat("------------- | ----- | ------- | --------------------\n")
for (r in list(list("GSE28521", res1), list("GSE28475", res2), list("GSE64018", res3))) {
  id <- r[[1]]
  res <- r[[2]]
  n_sig <- sum(res$module_results$significant, na.rm = TRUE)
  cat(sprintf("%-13s | %5d | %7d | %d\n", id, res$power, res$n_modules, n_sig))
}

cat("\nEnd time:", format(Sys.time()), "\n")
cat("All results saved to:", out_dir, "\n")
