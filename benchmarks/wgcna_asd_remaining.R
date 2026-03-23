#!/usr/bin/env Rscript
# WGCNA benchmark — remaining datasets (GSE28475 + GSE64018)
# GSE28521 already completed successfully.
# GSE28475 requires top-10K variable gene filter to fit in Pi 5 memory.

.libPaths("~/R/library")

suppressPackageStartupMessages({
  library(WGCNA)
})

allowWGCNAThreads(nThreads = 4)
options(stringsAsFactors = FALSE)

cat("=== WGCNA ASD Benchmark (remaining datasets) ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# --- Helper functions (same as main script) ---

parse_series_matrix <- function(gz_path) {
  con <- gzfile(gz_path, "rt")
  lines <- readLines(con)
  close(con)

  meta_lines <- lines[grepl("^!", lines)]
  data_start <- which(grepl("^\"ID_REF\"", lines))
  data_end_marker <- which(grepl("^!series_matrix_table_end", lines))
  if (length(data_end_marker) > 0) {
    data_end <- data_end_marker[1] - 1
  } else {
    data_end <- length(lines)
  }

  data_lines <- lines[data_start:data_end]
  tmp <- tempfile(fileext = ".tsv")
  writeLines(data_lines, tmp)
  expr <- read.delim(tmp, header = TRUE, row.names = 1, check.names = FALSE, quote = "\"")
  unlink(tmp)

  char_lines <- meta_lines[grepl("^!Sample_characteristics_ch1", meta_lines)]
  geo_line <- meta_lines[grepl("^!Sample_geo_accession", meta_lines)][1]
  geo_ids <- unlist(strsplit(geo_line, "\t"))[-1]
  geo_ids <- gsub("\"", "", geo_ids)

  list(expr = expr, char_lines = char_lines, geo_ids = geo_ids, meta = meta_lines)
}

extract_phenotypes <- function(char_lines, geo_ids, case_pattern, control_pattern) {
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
  annot <- read.delim(annot_path, header = TRUE, check.names = FALSE, comment.char = "#")
  gene_col <- NULL
  for (candidate in c("Gene symbol", "Gene Symbol", "Gene.symbol", "Gene.Symbol")) {
    if (candidate %in% colnames(annot)) { gene_col <- candidate; break }
  }
  if (is.null(gene_col)) stop("Cannot find gene symbol column in: ", annot_path)
  id_col <- colnames(annot)[1]

  annot_sub <- annot[annot[[id_col]] %in% rownames(expr), ]
  annot_sub <- annot_sub[annot_sub[[gene_col]] != "" & !is.na(annot_sub[[gene_col]]), ]
  annot_sub$gene <- sapply(strsplit(as.character(annot_sub[[gene_col]]), "///"), function(x) trimws(x[1]))
  annot_sub <- annot_sub[annot_sub$gene != "", ]

  expr_annot <- expr[as.character(annot_sub[[id_col]]), , drop = FALSE]
  expr_annot$gene <- annot_sub$gene

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

  datExpr <- as.data.frame(t(expr_genes))

  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK) {
    cat("Removing", sum(!gsg$goodGenes), "genes and", sum(!gsg$goodSamples), "samples\n")
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
    phenotypes <- phenotypes[gsg$goodSamples]
  }
  cat("After QC:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")

  cat("Picking soft threshold power...\n")
  powers <- c(1:10, seq(12, 20, 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 2, networkType = "signed")

  r2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  above_threshold <- which(r2 > 0.80)
  if (length(above_threshold) > 0) {
    power <- powers[above_threshold[1]]
  } else {
    power <- powers[which.max(r2)]
    cat("WARNING: No power reached R^2 > 0.80. Using best:", power, "R^2 =", max(r2), "\n")
  }
  cat("Selected soft threshold power:", power, "\n")

  write.csv(sft$fitIndices, file.path(out_dir, paste0(dataset_id, "_sft.csv")), row.names = FALSE)

  cat("Building TOM (this may take a while)...\n")
  t_start <- Sys.time()

  net <- blockwiseModules(
    datExpr, power = power, networkType = "signed", TOMType = "signed",
    minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE,
    verbose = 3, maxBlockSize = ncol(datExpr)
  )

  t_end <- Sys.time()
  cat("TOM + clustering took:", format(t_end - t_start), "\n")

  moduleLabels <- net$colors
  moduleColors <- labels2colors(moduleLabels)
  n_modules <- length(unique(moduleLabels)) - 1
  cat("Found", n_modules, "modules (+1 grey/unassigned)\n")

  trait <- ifelse(phenotypes == "ASD", 1, 0)
  names(trait) <- names(phenotypes)
  trait <- trait[rownames(datExpr)]

  MEs <- orderMEs(net$MEs)
  modTraitCor <- cor(MEs, trait, use = "p")
  modTraitPval <- corPvalueStudent(modTraitCor, nrow(datExpr))

  module_results <- data.frame(
    module = colnames(MEs),
    correlation = as.numeric(modTraitCor),
    pvalue = as.numeric(modTraitPval),
    n_genes = as.numeric(table(moduleLabels)[as.character(as.integer(gsub("ME", "", colnames(MEs))))])
  )
  module_results$significant <- module_results$pvalue < 0.05
  module_results$fdr <- p.adjust(module_results$pvalue, method = "BH")

  cat("\nModule-trait correlations:\n")
  print(module_results[order(module_results$pvalue), ], row.names = FALSE)

  gene_modules <- data.frame(gene = colnames(datExpr), module_number = moduleLabels, module_color = moduleColors)

  write.csv(gene_modules, file.path(out_dir, paste0(dataset_id, "_gene_modules.csv")), row.names = FALSE)
  write.csv(module_results, file.path(out_dir, paste0(dataset_id, "_module_trait.csv")), row.names = FALSE)

  cat("Results saved for", dataset_id, "\n")
  list(gene_modules = gene_modules, module_results = module_results, power = power, n_modules = n_modules)
}


out_dir <- "/home/kai001/riker-engine/benchmarks/results"

# ============================================================
# DATASET 2: GSE28475 — with top-10K variable gene filter
# ============================================================
cat("\n--- Loading GSE28475 ---\n")
d2 <- parse_series_matrix("/home/kai001/asd_full_run/data/bulk_geo/GSE28475-GPL6883_fixed_series_matrix.txt.gz")
p2 <- extract_phenotypes(d2$char_lines, d2$geo_ids, "diagnosis: autism", "diagnosis: control")
cat("Phenotypes:", table(p2), "\n")

expr2 <- map_probes_to_genes(d2$expr, "/home/kai001/asd_full_run/data/bulk_geo/GPL6883_clean.annot")
cat("Mapped to", nrow(expr2), "unique genes\n")

keep2 <- !is.na(p2) & (names(p2) %in% colnames(expr2))
expr2 <- expr2[, names(p2)[keep2]]
p2 <- p2[keep2]

# Filter to top 10K most variable genes (standard WGCNA practice)
# Required because full 18K gene TOM matrix exceeds Pi 5 8GB RAM
gene_vars <- apply(expr2, 1, var, na.rm = TRUE)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(10000, nrow(expr2))]
n_before <- nrow(expr2)
expr2 <- expr2[top_genes, ]
cat("Filtered from", n_before, "to", nrow(expr2), "most variable genes (memory constraint)\n")

# Force garbage collection before heavy computation
gc()

res2 <- run_wgcna(expr2, p2, "GSE28475", out_dir)

# Clean up memory before next dataset
rm(d2, expr2, p2, res2)
gc()


# ============================================================
# DATASET 3: GSE64018 — small dataset, no filtering needed
# ============================================================
cat("\n--- Loading GSE64018 ---\n")
d3 <- parse_series_matrix("/home/kai001/asd_full_run/data/bulk_geo/GSE64018_reconstructed_series_matrix.txt.gz")
p3 <- extract_phenotypes(d3$char_lines, d3$geo_ids, "diagnosis: asd", "diagnosis: ctl")
cat("Phenotypes:", table(p3), "\n")

annot3 <- read.delim("/home/kai001/asd_full_run/data/bulk_geo/GPL11154_ensembl.annot",
                     header = TRUE, check.names = FALSE)
ensg_to_gene <- annot3[annot3[["Gene Symbol"]] != "" & !is.na(annot3[["Gene Symbol"]]), ]
rownames(ensg_to_gene) <- ensg_to_gene$ID

common_ids <- intersect(rownames(d3$expr), ensg_to_gene$ID)
expr3 <- d3$expr[common_ids, ]
rownames(expr3) <- ensg_to_gene[common_ids, "Gene Symbol"]

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

cat("\n\nAll remaining datasets complete.\n")
cat("End time:", format(Sys.time()), "\n")
