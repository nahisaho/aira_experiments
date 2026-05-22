#!/usr/bin/env Rscript
# =============================================================================
# Module 2: RNA–Protein Expression Discordance & Translational Control
# =============================================================================
# Input:  RNA-seq TPM matrix, MaxQuant proteinGroups LFQ intensities
# Output: Discordance scores, translational efficiency, pathway enrichment
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(yaml)
})

config <- yaml::read_yaml("config/pipeline_config.yaml")
cfg    <- config$rna_protein_discordance

cat("=== Module 2: RNA–Protein Discordance Analysis ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load and harmonize RNA-seq & Proteomics data
# -----------------------------------------------------------------------------
load_and_harmonize <- function(rna_path, prot_path, min_samples = 10) {
  rna  <- read.delim(rna_path, row.names = 1, check.names = FALSE)
  prot <- read.delim(prot_path, stringsAsFactors = FALSE)

  # Extract LFQ intensities from MaxQuant proteinGroups
  lfq_cols <- grep("^LFQ\\.intensity\\.", colnames(prot), value = TRUE)
  prot_mat <- as.matrix(prot[, lfq_cols])
  rownames(prot_mat) <- prot$Gene.names
  colnames(prot_mat) <- gsub("^LFQ\\.intensity\\.", "", colnames(prot_mat))

  # log2 transform (replace 0 with NA)
  prot_mat[prot_mat == 0] <- NA
  prot_mat <- log2(prot_mat)

  # Match gene symbols
  rna_genes  <- rownames(rna)
  prot_genes <- rownames(prot_mat)
  # Take first gene name if multiple separated by ";"
  prot_genes_clean <- sapply(strsplit(prot_genes, ";"), `[`, 1)
  rownames(prot_mat) <- prot_genes_clean

  common_genes   <- intersect(rownames(rna), prot_genes_clean)
  common_samples <- intersect(colnames(rna), colnames(prot_mat))

  cat(sprintf("  Genes: RNA=%d, Protein=%d, Common=%d\n",
    nrow(rna), length(prot_genes_clean), length(common_genes)))
  cat(sprintf("  Samples: RNA=%d, Protein=%d, Common=%d\n",
    ncol(rna), ncol(prot_mat), length(common_samples)))

  rna_matched  <- as.matrix(rna[common_genes, common_samples])
  prot_matched <- prot_mat[common_genes, common_samples]

  # Filter: require detection in ≥ min_samples
  rna_detected  <- rowSums(!is.na(rna_matched) & rna_matched > 0) >= min_samples
  prot_detected <- rowSums(!is.na(prot_matched)) >= min_samples
  keep <- rna_detected & prot_detected

  cat(sprintf("  After detection filter: %d genes\n", sum(keep)))

  list(
    rna  = log2(rna_matched[keep, ] + 1),
    prot = prot_matched[keep, ],
    genes = common_genes[keep],
    samples = common_samples
  )
}

# -----------------------------------------------------------------------------
# Step 2: Compute per-gene RNA–protein correlation
# -----------------------------------------------------------------------------
compute_rna_protein_correlation <- function(rna_mat, prot_mat, method = "spearman") {
  n_genes <- nrow(rna_mat)
  results <- data.frame(
    gene      = rownames(rna_mat),
    rho       = numeric(n_genes),
    p_value   = numeric(n_genes),
    rna_mean  = numeric(n_genes),
    prot_mean = numeric(n_genes),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_genes)) {
    r <- rna_mat[i, ]
    p <- prot_mat[i, ]
    valid <- !is.na(r) & !is.na(p)
    if (sum(valid) >= 5) {
      ct <- cor.test(r[valid], p[valid], method = method)
      results$rho[i]     <- ct$estimate
      results$p_value[i] <- ct$p.value
    } else {
      results$rho[i]     <- NA
      results$p_value[i] <- NA
    }
    results$rna_mean[i]  <- mean(r, na.rm = TRUE)
    results$prot_mean[i] <- mean(p, na.rm = TRUE)
  }

  results$padj <- p.adjust(results$p_value, method = "BH")
  results <- results[order(results$rho, na.last = TRUE), ]

  cat(sprintf("  Median Spearman ρ: %.3f\n", median(results$rho, na.rm = TRUE)))
  cat(sprintf("  Positive correlation (ρ>0.3): %d genes\n", sum(results$rho > 0.3, na.rm = TRUE)))
  cat(sprintf("  Negative/weak (ρ<0.1): %d genes\n", sum(results$rho < 0.1, na.rm = TRUE)))

  return(results)
}

# -----------------------------------------------------------------------------
# Step 3: Translational efficiency estimation (regression residuals)
# -----------------------------------------------------------------------------
estimate_translational_efficiency <- function(rna_mat, prot_mat) {
  # For each sample, fit lm(protein ~ rna) per gene across samples
  # Residuals = protein levels unexplained by RNA → proxy for translational control
  n_genes <- nrow(rna_mat)
  te_matrix <- matrix(NA, nrow = n_genes, ncol = ncol(rna_mat))
  rownames(te_matrix) <- rownames(rna_mat)
  colnames(te_matrix) <- colnames(rna_mat)

  for (i in seq_len(n_genes)) {
    r <- rna_mat[i, ]
    p <- prot_mat[i, ]
    valid <- !is.na(r) & !is.na(p)
    if (sum(valid) >= 5) {
      fit <- lm(p[valid] ~ r[valid])
      te_matrix[i, valid] <- residuals(fit)
    }
  }

  # Summarize: mean absolute TE per gene → translational control strength
  te_summary <- data.frame(
    gene = rownames(te_matrix),
    te_mean_abs = rowMeans(abs(te_matrix), na.rm = TRUE),
    te_sd       = apply(te_matrix, 1, sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  te_summary <- te_summary[order(-te_summary$te_mean_abs), ]

  cat(sprintf("  Top translationally controlled genes:\n"))
  head_genes <- head(te_summary, 10)
  for (j in seq_len(nrow(head_genes))) {
    cat(sprintf("    %s: |TE| = %.3f\n", head_genes$gene[j], head_genes$te_mean_abs[j]))
  }

  return(list(matrix = te_matrix, summary = te_summary))
}

# -----------------------------------------------------------------------------
# Step 4: Identify discordant genes (RNA↑/Protein↓ or vice versa)
# -----------------------------------------------------------------------------
find_discordant_genes <- function(rna_mat, prot_mat, threshold = 1.5) {
  # Compute median z-score per gene for each layer
  rna_z  <- t(scale(t(rna_mat)))
  prot_z <- t(scale(t(prot_mat)))

  rna_median  <- apply(rna_z, 1, median, na.rm = TRUE)
  prot_median <- apply(prot_z, 1, median, na.rm = TRUE)

  discordance <- data.frame(
    gene         = rownames(rna_mat),
    rna_zscore   = rna_median,
    prot_zscore  = prot_median,
    delta        = prot_median - rna_median,
    stringsAsFactors = FALSE
  )

  discordance$category <- ifelse(
    abs(discordance$delta) > threshold,
    ifelse(discordance$delta > 0, "Protein_High", "RNA_High"),
    "Concordant"
  )

  cat(sprintf("  Discordant genes (|Δ| > %.1f): %d\n",
    threshold, sum(discordance$category != "Concordant")))
  cat(sprintf("    RNA > Protein: %d\n", sum(discordance$category == "RNA_High")))
  cat(sprintf("    Protein > RNA: %d\n", sum(discordance$category == "Protein_High")))

  return(discordance)
}

# -----------------------------------------------------------------------------
# Step 5: Pathway enrichment of discordant genes
# -----------------------------------------------------------------------------
enrich_discordant_pathways <- function(discordance_df, output_prefix = "results/discordance") {
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

  rna_high_genes <- discordance_df$gene[discordance_df$category == "RNA_High"]
  prot_high_genes <- discordance_df$gene[discordance_df$category == "Protein_High"]
  bg_genes <- discordance_df$gene

  run_enrichment <- function(genes, label) {
    if (length(genes) < 5) return(NULL)
    ego <- enrichGO(
      gene     = genes,
      universe = bg_genes,
      OrgDb    = org.Hs.eg.db,
      keyType  = "SYMBOL",
      ont      = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.05
    )
    if (!is.null(ego) && nrow(ego) > 0) {
      write.csv(as.data.frame(ego),
        paste0(output_prefix, "_", label, "_GO_BP.csv"), row.names = FALSE)
      cat(sprintf("  %s enriched GO terms: %d\n", label, nrow(ego)))
    }
    return(ego)
  }

  rna_enrich  <- run_enrichment(rna_high_genes, "RNA_high")
  prot_enrich <- run_enrichment(prot_high_genes, "Prot_high")

  return(list(rna_high = rna_enrich, prot_high = prot_enrich))
}

# -----------------------------------------------------------------------------
# Step 6: Visualization
# -----------------------------------------------------------------------------
plot_rna_protein <- function(cor_results, discordance_df, output_dir = "figures/") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # (A) Correlation distribution
  p1 <- ggplot(cor_results, aes(x = rho)) +
    geom_histogram(bins = 50, fill = "#2166ac", alpha = 0.7) +
    geom_vline(xintercept = median(cor_results$rho, na.rm = TRUE),
               linetype = "dashed", color = "red") +
    labs(x = "Spearman correlation (RNA vs Protein)",
         y = "Number of genes",
         title = "RNA–Protein Correlation Distribution") +
    theme_minimal(base_size = 14) +
    annotate("text", x = 0.6, y = Inf, vjust = 2,
      label = sprintf("Median ρ = %.3f", median(cor_results$rho, na.rm = TRUE)))

  ggsave(file.path(output_dir, "rna_protein_correlation_hist.pdf"),
    p1, width = 8, height = 5)

  # (B) Discordance scatter
  p2 <- ggplot(discordance_df, aes(x = rna_zscore, y = prot_zscore, color = category)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c(
      "Concordant"   = "grey70",
      "RNA_High"     = "#d6604d",
      "Protein_High" = "#4393c3"
    )) +
    labs(x = "RNA Z-score (median)", y = "Protein Z-score (median)",
         title = "RNA–Protein Expression Discordance") +
    theme_minimal(base_size = 14)

  ggsave(file.path(output_dir, "rna_protein_discordance_scatter.pdf"),
    p2, width = 8, height = 7)

  cat("  Figures saved to", output_dir, "\n")
}

# =============================================================================
# Main
# =============================================================================
main <- function() {
  cat("\n--- Loading Data ---\n")
  dat <- load_and_harmonize(cfg$rnaseq_counts, cfg$protein_intensities,
    min_samples = cfg$min_samples_detected)

  cat("\n--- RNA–Protein Correlation ---\n")
  cor_res <- compute_rna_protein_correlation(dat$rna, dat$prot,
    method = cfg$correlation_method)
  write.csv(cor_res, "results/rna_protein_correlations.csv", row.names = FALSE)

  cat("\n--- Translational Efficiency ---\n")
  te <- estimate_translational_efficiency(dat$rna, dat$prot)
  write.csv(te$summary, "results/translational_efficiency_scores.csv", row.names = FALSE)

  cat("\n--- Discordance Analysis ---\n")
  disc <- find_discordant_genes(dat$rna, dat$prot,
    threshold = cfg$discordance_threshold)
  write.csv(disc, "results/discordant_genes.csv", row.names = FALSE)

  cat("\n--- Pathway Enrichment ---\n")
  enrich <- enrich_discordant_pathways(disc)

  cat("\n--- Visualization ---\n")
  plot_rna_protein(cor_res, disc)

  cat("\n=== Module 2 complete ===\n")
}

if (!interactive()) main()
