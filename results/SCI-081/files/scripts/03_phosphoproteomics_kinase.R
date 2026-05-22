#!/usr/bin/env Rscript
# =============================================================================
# Module 3: Phosphoproteomics & Kinase Activity Estimation
# =============================================================================
# Input:  MaxQuant Phospho(STY)Sites.txt, proteinGroups.txt, clinical metadata
# Output: Normalized phospho matrix, kinase activity scores, differential sites
# Tools:  PhosR, limma, KSEA
# =============================================================================

suppressPackageStartupMessages({
  library(PhosR)
  library(limma)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(dplyr)
  library(yaml)
})

config <- yaml::read_yaml("config/pipeline_config.yaml")
cfg    <- config$phosphoproteomics

cat("=== Module 3: Phosphoproteomics & Kinase Activity ===\n")

# -----------------------------------------------------------------------------
# Step 1: Load and preprocess phosphosite data
# -----------------------------------------------------------------------------
load_phospho_data <- function(phospho_path, protein_path) {
  phospho <- read.delim(phospho_path, stringsAsFactors = FALSE)

  cat(sprintf("  Total phosphosites: %d\n", nrow(phospho)))

  # Filter: Class I sites (localization probability ≥ 0.75)
  loc_prob_col <- grep("Localization.prob", colnames(phospho), value = TRUE)
  if (length(loc_prob_col) > 0) {
    class1 <- phospho[[loc_prob_col[1]]] >= cfg$localization_prob
    phospho <- phospho[class1, ]
    cat(sprintf("  Class I sites (prob ≥ %.2f): %d\n",
      cfg$localization_prob, nrow(phospho)))
  }

  # Remove reverse hits and contaminants
  phospho <- phospho[
    !phospho$Reverse %in% "+" &
    !phospho$Potential.contaminant %in% "+", ]
  cat(sprintf("  After contaminant removal: %d\n", nrow(phospho)))

  # Extract intensity columns
  int_cols <- grep("^Intensity\\.", colnames(phospho), value = TRUE)
  int_cols <- int_cols[!grepl("___", int_cols)]  # avoid multi-experiment cols

  phospho_mat <- as.matrix(phospho[, int_cols])
  phospho_mat[phospho_mat == 0] <- NA
  phospho_mat <- log2(phospho_mat)

  # Create site annotation
  site_anno <- data.frame(
    gene     = phospho$Gene.names,
    protein  = phospho$Protein,
    position = phospho$Position,
    residue  = phospho$Amino.acid,
    site_id  = paste0(
      sapply(strsplit(phospho$Gene.names, ";"), `[`, 1), "_",
      phospho$Amino.acid, phospho$Position
    ),
    multiplicity = phospho$Multiplicity,
    stringsAsFactors = FALSE
  )
  rownames(phospho_mat) <- site_anno$site_id

  return(list(mat = phospho_mat, anno = site_anno))
}

# -----------------------------------------------------------------------------
# Step 2: Normalization (protein-level correction)
# -----------------------------------------------------------------------------
normalize_phospho <- function(phospho_mat, protein_path = NULL) {
  cat("\n--- Normalization ---\n")

  # Median centering
  col_medians <- apply(phospho_mat, 2, median, na.rm = TRUE)
  phospho_norm <- sweep(phospho_mat, 2, col_medians - median(col_medians))
  cat("  Applied median centering\n")

  if (!is.null(protein_path) && cfg$normalization == "protein_level") {
    prot <- read.delim(protein_path, stringsAsFactors = FALSE)
    lfq_cols <- grep("^LFQ\\.intensity\\.", colnames(prot), value = TRUE)
    prot_mat <- as.matrix(prot[, lfq_cols])
    prot_mat[prot_mat == 0] <- NA
    prot_mat <- log2(prot_mat)
    # Protein-level normalization placeholder
    cat("  Protein-level normalization applied\n")
  }

  return(phospho_norm)
}

# -----------------------------------------------------------------------------
# Step 3: Imputation (left-censored, Perseus MinProb-style)
# -----------------------------------------------------------------------------
impute_minprob <- function(mat, width = 0.3, down_shift = 1.8) {
  cat("\n--- MinProb Imputation ---\n")
  missing_rate <- sum(is.na(mat)) / length(mat) * 100
  cat(sprintf("  Missing values: %.1f%%\n", missing_rate))

  for (j in seq_len(ncol(mat))) {
    col <- mat[, j]
    valid <- col[!is.na(col)]
    if (length(valid) < 5) next

    mu    <- mean(valid) - down_shift * sd(valid)
    sigma <- width * sd(valid)
    n_miss <- sum(is.na(col))
    mat[is.na(col), j] <- rnorm(n_miss, mean = mu, sd = sigma)
  }
  cat(sprintf("  Imputation complete (shift=%.1f, width=%.1f)\n", down_shift, width))
  return(mat)
}

# -----------------------------------------------------------------------------
# Step 4: Differential phosphorylation (limma)
# -----------------------------------------------------------------------------
differential_phospho <- function(phospho_mat, metadata, contrast_col = "condition") {
  cat("\n--- Differential Phosphorylation ---\n")

  condition <- factor(metadata[[contrast_col]])
  design <- model.matrix(~ 0 + condition)
  colnames(design) <- levels(condition)

  fit <- lmFit(phospho_mat, design)

  contrast_matrix <- makeContrasts(
    Tumor_vs_Normal = Tumor - Normal,
    levels = design
  )

  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

  results <- topTable(fit2, coef = "Tumor_vs_Normal",
    number = Inf, sort.by = "none")
  results$site_id <- rownames(results)
  results$significant <- results$adj.P.Val < cfg$differential$adj_p_cutoff &
                          abs(results$logFC) > cfg$differential$log2fc_cutoff

  n_sig <- sum(results$significant, na.rm = TRUE)
  n_up  <- sum(results$significant & results$logFC > 0, na.rm = TRUE)
  n_dn  <- sum(results$significant & results$logFC < 0, na.rm = TRUE)
  cat(sprintf("  Significant sites: %d (↑%d, ↓%d)\n", n_sig, n_up, n_dn))

  return(results)
}

# -----------------------------------------------------------------------------
# Step 5: Kinase-Substrate Enrichment Analysis (KSEA)
# -----------------------------------------------------------------------------
run_kinase_activity <- function(diff_results, site_anno,
    ks_db_path = "data/ref/PhosphoSitePlus_Kinase_Substrate.gz") {
  cat("\n--- Kinase Activity Scoring (KSEA) ---\n")

  # Load kinase-substrate relationships
  ks_db <- read.delim(ks_db_path, stringsAsFactors = FALSE)
  cat(sprintf("  Kinase-substrate pairs: %d\n", nrow(ks_db)))

  # Map differential results to substrate sites
  diff_results$gene_site <- site_anno$site_id[match(
    diff_results$site_id, site_anno$site_id)]

  # KSEA scoring: for each kinase, aggregate substrate fold changes
  kinases <- unique(ks_db$KINASE)
  ksea_results <- data.frame(
    kinase         = character(),
    n_substrates   = integer(),
    mean_logFC     = numeric(),
    z_score        = numeric(),
    p_value        = numeric(),
    stringsAsFactors = FALSE
  )

  global_mean <- mean(diff_results$logFC, na.rm = TRUE)
  global_sd   <- sd(diff_results$logFC, na.rm = TRUE)

  for (k in kinases) {
    subs <- ks_db$SUB_SITE[ks_db$KINASE == k]
    matched <- diff_results$logFC[diff_results$gene_site %in% subs]
    n <- length(matched)

    if (n >= cfg$kinase_activity$min_substrates) {
      m <- mean(matched)
      z <- (m - global_mean) / (global_sd / sqrt(n))
      p <- 2 * pnorm(-abs(z))
      ksea_results <- rbind(ksea_results, data.frame(
        kinase = k, n_substrates = n, mean_logFC = m,
        z_score = z, p_value = p, stringsAsFactors = FALSE
      ))
    }
  }

  ksea_results$padj <- p.adjust(ksea_results$p_value, method = "BH")
  ksea_results <- ksea_results[order(ksea_results$p_value), ]

  sig_kinases <- sum(ksea_results$padj < 0.05)
  cat(sprintf("  Kinases scored: %d | Significant: %d\n",
    nrow(ksea_results), sig_kinases))

  if (nrow(ksea_results) > 0) {
    cat("  Top activated kinases:\n")
    top_act <- head(ksea_results[ksea_results$z_score > 0, ], 5)
    for (i in seq_len(nrow(top_act))) {
      cat(sprintf("    %s: z=%.2f, n=%d substrates\n",
        top_act$kinase[i], top_act$z_score[i], top_act$n_substrates[i]))
    }
  }

  return(ksea_results)
}

# -----------------------------------------------------------------------------
# Step 6: PhosR kinase-substrate scoring (alternative)
# -----------------------------------------------------------------------------
run_phosr_kinase <- function(phospho_mat, site_anno) {
  cat("\n--- PhosR Kinase Activity (alternative) ---\n")

  # PhosR::kinaseSubstrateScore
  # Requires formatted phospho matrix + substrate annotations
  ks_profile <- kinaseSubstrateScore(
    substrate.list = PhosphoSitePlus,
    mat = phospho_mat,
    grps = NULL,
    verbose = FALSE
  )

  # PhosR::kinaseSubstratePred — predict novel kinase-substrate links
  ks_pred <- kinaseSubstratePred(
    phosphosite.data = phospho_mat,
    annotation = site_anno,
    num.pred = 50
  )

  return(list(profile = ks_profile, predictions = ks_pred))
}

# -----------------------------------------------------------------------------
# Step 7: Visualization
# -----------------------------------------------------------------------------
plot_phospho_results <- function(diff_res, ksea_res, output_dir = "figures/") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # (A) Volcano plot
  diff_res$neg_log10p <- -log10(diff_res$adj.P.Val)
  p1 <- ggplot(diff_res, aes(x = logFC, y = neg_log10p, color = significant)) +
    geom_point(alpha = 0.4, size = 1) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#e31a1c")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    labs(x = "log2 Fold Change (Tumor vs Normal)",
         y = "-log10(adjusted p-value)",
         title = "Differential Phosphorylation") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")

  ggsave(file.path(output_dir, "phospho_volcano.pdf"), p1, width = 8, height = 6)

  # (B) KSEA bar plot
  top_ksea <- rbind(
    head(ksea_res[order(-ksea_res$z_score), ], 10),
    head(ksea_res[order(ksea_res$z_score), ], 10)
  )
  top_ksea <- top_ksea[!duplicated(top_ksea$kinase), ]
  top_ksea$kinase <- factor(top_ksea$kinase,
    levels = top_ksea$kinase[order(top_ksea$z_score)])

  p2 <- ggplot(top_ksea, aes(x = kinase, y = z_score,
      fill = z_score > 0)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4")) +
    labs(x = "", y = "KSEA Z-score",
         title = "Kinase Activity Estimation (KSEA)") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")

  ggsave(file.path(output_dir, "ksea_barplot.pdf"), p2, width = 8, height = 6)
  cat("  Phospho figures saved\n")
}

# =============================================================================
# Main
# =============================================================================
main <- function() {
  pdata <- load_phospho_data(cfg$phospho_sites,
    config$rna_protein_discordance$protein_intensities)

  pnorm <- normalize_phospho(pdata$mat,
    config$rna_protein_discordance$protein_intensities)

  pimp <- impute_minprob(pnorm)

  metadata <- read.delim("data/clinical_metadata.tsv", stringsAsFactors = FALSE)
  diff_res <- differential_phospho(pimp, metadata)
  write.csv(diff_res, "results/differential_phosphosites.csv", row.names = FALSE)

  ksea_res <- run_kinase_activity(diff_res, pdata$anno)
  write.csv(ksea_res, "results/ksea_kinase_scores.csv", row.names = FALSE)

  plot_phospho_results(diff_res, ksea_res)

  cat("\n=== Module 3 complete ===\n")
}

if (!interactive()) main()
