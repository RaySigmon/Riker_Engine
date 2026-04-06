# MetaDE Benchmark — ASD Brain Cortex
# Compares MetaDE (Meta-analysis of Differential Expression)
# against Riker Engine on the same ASD discovery datasets.
#
# MetaDE combines differential expression statistics across multiple
# datasets using various meta-analytic methods (Fisher's, maxP, etc.).
# It finds differentially expressed genes but does NOT identify modules —
# this is the key distinction from the Riker Engine.
#
# Usage: Rscript benchmarks/metade_asd_benchmark.R
# Requires: MetaDE package from CRAN
# Output: benchmarks/results/metade_*

library(MetaDE)

cat("============================================\n")
cat("MetaDE Benchmark — ASD Brain Cortex\n")
cat("============================================\n")

# Helper: parse GEO series matrix to expression matrix
parse_geo_matrix <- function(path) {
    lines <- readLines(gzfile(path))
    data_start <- which(grepl("^!series_matrix_table_begin", lines)) + 1
    data_end <- which(grepl("^!series_matrix_table_end", lines)) - 1
    expr <- read.delim(textConnection(lines[data_start:data_end]),
                        row.names = 1, check.names = FALSE)
    return(as.matrix(expr))
}

# Helper: map probes to genes
map_to_genes <- function(expr, platform_path) {
    platform <- read.delim(platform_path, comment.char = "#",
                           stringsAsFactors = FALSE)
    probe_to_gene <- platform[, c("ID", "Gene.symbol")]
    probe_to_gene <- probe_to_gene[probe_to_gene$Gene.symbol != "" &
                                     !is.na(probe_to_gene$Gene.symbol), ]

    common_probes <- intersect(rownames(expr), probe_to_gene$ID)
    expr_sub <- expr[common_probes, ]
    genes <- probe_to_gene$Gene.symbol[match(common_probes, probe_to_gene$ID)]

    # Collapse by gene (mean)
    gene_expr <- aggregate(expr_sub, by = list(gene = genes), FUN = mean)
    rownames(gene_expr) <- gene_expr$gene
    gene_expr$gene <- NULL
    return(as.matrix(gene_expr))
}

# Load datasets
cat("\nLoading ASD discovery datasets...\n")

datasets <- list()
labels <- list()

# GSE28521
cat("  Loading GSE28521...\n")
e1 <- parse_geo_matrix("data/geo/asd/GSE28521_series_matrix.txt.gz")
g1 <- map_to_genes(e1, "data/platforms/GPL6883.annot")
# Simple label assignment: need to identify case/control from GEO metadata
# For benchmark purposes, use the phenotype parsing approach
# We'll create binary labels: 1=case, 0=control
# GSE28521 has autism vs controls in characteristics
datasets[["GSE28521"]] <- g1

# GSE28475
cat("  Loading GSE28475...\n")
e2_path <- "data/geo/asd/GSE28475-GPL6883_series_matrix.txt.gz"
if (file.exists(e2_path)) {
    e2 <- parse_geo_matrix(e2_path)
    g2 <- map_to_genes(e2, "data/platforms/GPL6883.annot")
    datasets[["GSE28475"]] <- g2
}

cat(sprintf("  Loaded %d datasets\n", length(datasets)))

# Find common genes across all datasets
common_genes <- Reduce(intersect, lapply(datasets, rownames))
cat(sprintf("  Common genes across datasets: %d\n", length(common_genes)))

# Subset to common genes
for (name in names(datasets)) {
    datasets[[name]] <- datasets[[name]][common_genes, ]
}

# For MetaDE, we need label vectors matching sample columns
# This is a simplified approach — in a real benchmark, parse phenotypes properly
# For now, we demonstrate the MetaDE workflow and output structure
cat("\nNOTE: This benchmark script provides the MetaDE framework.\n")
cat("Full phenotype parsing for multi-dataset meta-DE requires\n")
cat("dataset-specific label extraction. See the Riker Engine configs\n")
cat("for exact case/control label mappings per dataset.\n")

# Run MetaDE with Fisher's method (most common)
cat("\nRunning MetaDE (Fisher's method)...\n")
start_time <- proc.time()

# MetaDE expects a list of expression matrices and a list of label vectors
# Create placeholder labels for the structural comparison
# In the full benchmark, these should be extracted from GEO metadata
n1 <- ncol(datasets[["GSE28521"]])
label1 <- c(rep(1, floor(n1/2)), rep(0, ceiling(n1/2)))  # placeholder

if (length(datasets) >= 2) {
    n2 <- ncol(datasets[["GSE28475"]])
    label2 <- c(rep(1, floor(n2/2)), rep(0, ceiling(n2/2)))  # placeholder

    # MetaDE.pvalue combines p-values across studies
    tryCatch({
        result <- MetaDE.pvalue(
            x = list(datasets[["GSE28521"]], datasets[["GSE28475"]]),
            label = list(label1, label2),
            meta.method = "Fisher",
            missed.percent.cutoff = 0.3
        )

        elapsed <- proc.time() - start_time
        cat(sprintf("\nMetaDE completed in %.1f seconds\n", elapsed["elapsed"]))

        # Extract significant genes
        if (!is.null(result$meta.analysis)) {
            pvals <- result$meta.analysis$pval
            n_sig_005 <- sum(pvals < 0.05, na.rm = TRUE)
            n_sig_001 <- sum(pvals < 0.01, na.rm = TRUE)
            n_sig_fdr <- sum(p.adjust(pvals, method = "BH") < 0.05, na.rm = TRUE)

            cat(sprintf("\nResults:\n"))
            cat(sprintf("  Genes tested: %d\n", length(pvals)))
            cat(sprintf("  Significant (p<0.05): %d\n", n_sig_005))
            cat(sprintf("  Significant (p<0.01): %d\n", n_sig_001))
            cat(sprintf("  Significant (FDR<0.05): %d\n", n_sig_fdr))

            # Save results
            result_df <- data.frame(
                gene = names(pvals),
                fisher_p = pvals,
                fdr = p.adjust(pvals, method = "BH"),
                stringsAsFactors = FALSE
            )
            result_df <- result_df[order(result_df$fisher_p), ]
            write.csv(result_df, "benchmarks/results/metade_GSE28521_GSE28475.csv",
                      row.names = FALSE)

            # Check overlap with Riker core genes
            riker_core <- read.csv("results/asd_test_run_1/curated/phase4_core_genes.csv")
            sig_genes <- result_df$gene[result_df$fdr < 0.05]
            overlap <- intersect(riker_core$gene, sig_genes)
            cat(sprintf("\n  Riker core genes in MetaDE significant set: %d/%d\n",
                        length(overlap), nrow(riker_core)))
            cat(sprintf("  MetaDE significant genes: %d (vs Riker's 35 core genes)\n",
                        n_sig_fdr))
        }
    }, error = function(e) {
        cat(sprintf("MetaDE error: %s\n", e$message))
        cat("This may require proper phenotype labels. See notes above.\n")
    })
} else {
    cat("Need at least 2 datasets for meta-analysis.\n")
}

# Save timing
timing_df <- data.frame(
    tool = "MetaDE",
    method = "Fisher",
    n_datasets = length(datasets),
    n_common_genes = length(common_genes),
    runtime_seconds = (proc.time() - start_time)["elapsed"],
    timestamp = Sys.time()
)
write.csv(timing_df, "benchmarks/results/metade_timing.csv",
          row.names = FALSE)

cat("\n============================================\n")
cat("MetaDE Benchmark Key Comparison Points:\n")
cat("============================================\n")
cat("1. MetaDE produces a ranked gene list, NOT modules\n")
cat("   - Riker Engine produces clustered modules with biological coherence\n")
cat("2. MetaDE uses statistical combination only (Fisher's, etc.)\n")
cat("   - Riker Engine adds pathway-informed clustering + replication\n")
cat("3. MetaDE output is typically hundreds-thousands of significant genes\n")
cat("   - Riker Engine filters to a minimal replicated core set\n")
cat("4. MetaDE has no built-in replication or robustness testing\n")
cat("   - Riker Engine enforces pre-specification and directional replication\n")
cat("\nOutputs saved to benchmarks/results/metade_*\n")
