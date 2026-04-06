# MEGENA Benchmark — ASD Brain Cortex
# Compares MEGENA (Multiscale Embedded Gene Co-expression Network Analysis)
# against Riker Engine on the same ASD discovery datasets.
#
# MEGENA builds a planar filtered network and identifies multiscale modules
# using a hierarchical clustering approach. It is designed for single-dataset
# analysis, like WGCNA, but uses a different network construction method.
#
# Usage: Rscript benchmarks/megena_asd_benchmark.R
# Requires: MEGENA package from Bioconductor
# Output: benchmarks/results/megena_*

library(MEGENA)

cat("============================================\n")
cat("MEGENA Benchmark — ASD Brain Cortex\n")
cat("============================================\n")

# Load the same ASD discovery dataset used for WGCNA benchmark
# GSE28521 — Voineagu et al., brain cortex
cat("\nLoading GSE28521 expression data...\n")

# Read the series matrix (same data as WGCNA benchmark)
# Adjust path based on where data was downloaded
data_path <- "data/geo/asd/GSE28521_series_matrix.txt.gz"
if (!file.exists(data_path)) {
    cat("ERROR: Data file not found. Run: python scripts/download_data.py asd\n")
    quit(status = 1)
}

# Parse GEO series matrix
lines <- readLines(gzfile(data_path))
data_start <- which(grepl("^!series_matrix_table_begin", lines)) + 1
data_end <- which(grepl("^!series_matrix_table_end", lines)) - 1
expr_data <- read.delim(textConnection(lines[data_start:data_end]),
                         row.names = 1, check.names = FALSE)

cat(sprintf("Expression matrix: %d probes x %d samples\n",
            nrow(expr_data), ncol(expr_data)))

# Map probes to genes using platform annotation
platform_path <- "data/platforms/GPL6883.annot"
if (!file.exists(platform_path)) {
    cat("ERROR: Platform file not found. Run: python scripts/download_data.py --platforms\n")
    quit(status = 1)
}

platform <- read.delim(platform_path, comment.char = "#", stringsAsFactors = FALSE)
probe_to_gene <- platform[, c("ID", "Gene.symbol")]
probe_to_gene <- probe_to_gene[probe_to_gene$Gene.symbol != "" &
                                 !is.na(probe_to_gene$Gene.symbol), ]

# Merge and collapse to gene-level (mean per gene)
expr_data$probe <- rownames(expr_data)
merged <- merge(expr_data, probe_to_gene, by.x = "probe", by.y = "ID")
merged$probe <- NULL
gene_expr <- aggregate(. ~ Gene.symbol, data = merged, FUN = mean)
rownames(gene_expr) <- gene_expr$Gene.symbol
gene_expr$Gene.symbol <- NULL

cat(sprintf("Gene-level matrix: %d genes x %d samples\n",
            nrow(gene_expr), ncol(gene_expr)))

# Run MEGENA
cat("\nRunning MEGENA pipeline...\n")
start_time <- proc.time()

# Step 1: Calculate correlation
cat("  Step 1: Calculating correlations...\n")
cor_mat <- cor(t(gene_expr), method = "pearson")

# Step 2: Planar filtered network
cat("  Step 2: Building planar filtered network...\n")
pfn <- calculate.PFN(cor_mat, doPar = 1, num.cores = 1)

# Step 3: Multiscale clustering
cat("  Step 3: Multiscale module detection...\n")
megena_result <- do.MEGENA(pfn,
                            mod.pval = 0.05,
                            hub.pval = 0.05,
                            remove.unsig = TRUE,
                            min.size = 10,
                            max.size = nrow(gene_expr) / 2,
                            doPar = FALSE)

elapsed <- proc.time() - start_time
cat(sprintf("\nMEGENA completed in %.1f minutes\n", elapsed["elapsed"] / 60))

# Extract results
modules <- megena_result$module.output
if (!is.null(modules)) {
    n_modules <- length(modules$modules)
    module_sizes <- sapply(modules$modules, length)

    cat(sprintf("\nResults:\n"))
    cat(sprintf("  Total modules: %d\n", n_modules))
    cat(sprintf("  Module sizes: %d - %d (median %d)\n",
                min(module_sizes), max(module_sizes), median(module_sizes)))
    cat(sprintf("  Total genes in modules: %d\n", sum(module_sizes)))

    # Save module assignments
    module_df <- data.frame(gene = character(), module = integer(),
                            stringsAsFactors = FALSE)
    for (i in seq_along(modules$modules)) {
        genes <- modules$modules[[i]]
        module_df <- rbind(module_df,
                           data.frame(gene = genes, module = i,
                                      stringsAsFactors = FALSE))
    }
    write.csv(module_df, "benchmarks/results/megena_GSE28521_modules.csv",
              row.names = FALSE)

    # Check overlap with Riker core genes
    riker_core <- read.csv("results/asd_test_run_1/curated/phase4_core_genes.csv")
    riker_genes <- riker_core$gene
    in_megena <- riker_genes[riker_genes %in% module_df$gene]
    cat(sprintf("\n  Riker core genes in MEGENA modules: %d/%d\n",
                length(in_megena), length(riker_genes)))
} else {
    cat("\nNo significant modules found by MEGENA.\n")
}

# Save timing
timing_df <- data.frame(
    tool = "MEGENA",
    dataset = "GSE28521",
    n_genes = nrow(gene_expr),
    n_samples = ncol(gene_expr),
    n_modules = ifelse(!is.null(modules), length(modules$modules), 0),
    runtime_seconds = elapsed["elapsed"],
    timestamp = Sys.time()
)
write.csv(timing_df, "benchmarks/results/megena_GSE28521_timing.csv",
          row.names = FALSE)

cat("\nOutputs saved to benchmarks/results/megena_*\n")
