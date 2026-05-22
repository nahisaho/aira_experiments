#!/usr/bin/env Rscript

set.seed(42)

run_cibersortx_deconvolution <- function(base_dir = NULL) {
  message("[02] Starting immune deconvolution analysis...")

  suppressPackageStartupMessages({
    pkgs <- c("immunedeconv", "ggplot2", "tidyverse", "reshape2", "ComplexHeatmap", "RColorBrewer")
    missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      stop(sprintf("Missing required packages: %s", paste(missing_pkgs, collapse = ", ")))
    }
    invisible(lapply(pkgs, library, character.only = TRUE))
  })
  message("[02] CIBERSORTx signal is simulated in this script; no external CIBERSORT binary is required.")

  resolve_base_dir <- function(path_hint = NULL) {
    if (!is.null(path_hint)) return(normalizePath(path_hint, mustWork = TRUE))
    wd <- normalizePath(getwd(), mustWork = TRUE)
    if (basename(wd) == "scripts") dirname(wd) else wd
  }

  base_dir <- resolve_base_dir(base_dir)
  fig_dir <- file.path(base_dir, "figures")
  res_dir <- file.path(base_dir, "results")
  data_dir <- file.path(base_dir, "data")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  lm22_types <- c(
    "B_cells_naive", "B_cells_memory", "Plasma_cells", "T_cells_CD8", "T_cells_CD4_naive",
    "T_cells_CD4_memory_resting", "T_cells_CD4_memory_activated", "Tregs", "T_cells_gamma_delta",
    "NK_cells_resting", "NK_cells_activated", "Monocytes", "Macrophages_M0", "Macrophages_M1",
    "Macrophages_M2", "Dendritic_cells_resting", "Dendritic_cells_activated", "Mast_cells_resting",
    "Mast_cells_activated", "Eosinophils", "Neutrophils", "Follicular_helper_T"
  )
  group_levels <- c("RA_active", "RA_remission", "SLE", "HC")
  sample_ids <- sprintf("S%03d", seq_len(80))
  metadata <- data.frame(
    sample_id = sample_ids,
    group = factor(rep(group_levels, each = 20), levels = group_levels),
    DAS28 = c(rnorm(20, 5.9, 0.6), rnorm(20, 2.8, 0.4), rnorm(20, 4.6, 0.7), rnorm(20, 1.6, 0.3)),
    CRP = c(rlnorm(20, 1.5, 0.35), rlnorm(20, 0.6, 0.25), rlnorm(20, 1.2, 0.30), rlnorm(20, -0.1, 0.20)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample_id

  base_profile <- c(
    0.05, 0.04, 0.02, 0.08, 0.04,
    0.06, 0.04, 0.04, 0.02, 0.05,
    0.03, 0.08, 0.06, 0.06, 0.07,
    0.05, 0.03, 0.03, 0.01, 0.01,
    0.07, 0.04
  )
  names(base_profile) <- lm22_types

  adjust_profile <- function(group_name) {
    alpha <- base_profile
    if (group_name == "RA_active") {
      alpha["T_cells_CD8"] <- alpha["T_cells_CD8"] * 2.0
      alpha["Tregs"] <- alpha["Tregs"] * 0.45
      alpha["NK_cells_resting"] <- alpha["NK_cells_resting"] * 0.60
      alpha["NK_cells_activated"] <- alpha["NK_cells_activated"] * 0.65
      alpha["Macrophages_M1"] <- alpha["Macrophages_M1"] * 2.2
      alpha["Macrophages_M2"] <- alpha["Macrophages_M2"] * 0.65
      alpha["Monocytes"] <- alpha["Monocytes"] * 1.30
      alpha["Neutrophils"] <- alpha["Neutrophils"] * 1.35
    }
    if (group_name == "RA_remission") {
      alpha["Tregs"] <- alpha["Tregs"] * 0.70
      alpha["T_cells_CD8"] <- alpha["T_cells_CD8"] * 1.25
      alpha["Macrophages_M1"] <- alpha["Macrophages_M1"] * 1.25
      alpha["Macrophages_M2"] <- alpha["Macrophages_M2"] * 0.85
      alpha["NK_cells_resting"] <- alpha["NK_cells_resting"] * 0.85
    }
    if (group_name == "SLE") {
      alpha["Tregs"] <- alpha["Tregs"] * 0.55
      alpha["Plasma_cells"] <- alpha["Plasma_cells"] * 1.8
      alpha["B_cells_memory"] <- alpha["B_cells_memory"] * 1.4
      alpha["NK_cells_resting"] <- alpha["NK_cells_resting"] * 0.70
      alpha["Dendritic_cells_activated"] <- alpha["Dendritic_cells_activated"] * 1.5
    }
    alpha
  }

  rdirichlet_simple <- function(alpha) {
    draw <- rgamma(length(alpha), shape = alpha * 80, rate = 1)
    draw / sum(draw)
  }

  proportion_matrix <- sapply(metadata$group, function(g) rdirichlet_simple(adjust_profile(as.character(g))))
  colnames(proportion_matrix) <- metadata$sample_id
  rownames(proportion_matrix) <- lm22_types
  write.csv(t(proportion_matrix), file.path(data_dir, "synthetic_cibersortx_proportions.csv"), row.names = TRUE)
  write.csv(cbind(sample_id = metadata$sample_id, t(proportion_matrix)), file.path(res_dir, "cibersortx_sample_proportions.csv"), row.names = FALSE)

  long_df <- reshape2::melt(proportion_matrix)
  colnames(long_df) <- c("cell_type", "sample_id", "proportion")
  long_df <- dplyr::left_join(long_df, metadata, by = "sample_id")

  dunn_posthoc <- function(df_cell) {
    if (requireNamespace("FSA", quietly = TRUE)) {
      out <- FSA::dunnTest(proportion ~ group, data = df_cell, method = "bh")$res
      out$comparison <- out$Comparison
      out$adj_p <- out$P.adj
      return(out[, c("comparison", "Z", "P.unadj", "adj_p")])
    }
    pw <- pairwise.wilcox.test(df_cell$proportion, df_cell$group, p.adjust.method = "BH")
    out <- reshape2::melt(pw$p.value, varnames = c("group1", "group2"), value.name = "adj_p")
    out <- out[!is.na(out$adj_p), ]
    out$comparison <- paste(out$group1, out$group2, sep = " vs ")
    out$Z <- NA_real_
    out$P.unadj <- NA_real_
    out[, c("comparison", "Z", "P.unadj", "adj_p")]
  }

  result_list <- lapply(split(long_df, long_df$cell_type), function(df_cell) {
    kw <- kruskal.test(proportion ~ group, data = df_cell)
    n <- nrow(df_cell)
    k <- length(unique(df_cell$group))
    epsilon_sq <- max((kw$statistic - k + 1) / (n - k), 0)
    posthoc <- dunn_posthoc(df_cell)
    if (nrow(posthoc) == 0) {
      posthoc <- data.frame(comparison = NA_character_, Z = NA_real_, P.unadj = NA_real_, adj_p = NA_real_)
    }
    transform(posthoc,
      cell_type = unique(df_cell$cell_type),
      kruskal_p = kw$p.value,
      epsilon_sq = unname(epsilon_sq),
      ra_active_median = median(df_cell$proportion[df_cell$group == "RA_active"]),
      hc_median = median(df_cell$proportion[df_cell$group == "HC"]),
      median_diff = median(df_cell$proportion[df_cell$group == "RA_active"]) - median(df_cell$proportion[df_cell$group == "HC"])
    )
  })
  results_df <- dplyr::bind_rows(result_list) %>%
    dplyr::relocate(cell_type, comparison, kruskal_p, epsilon_sq)
  write.csv(results_df, file.path(res_dir, "cibersortx_results.csv"), row.names = FALSE)

  mean_prop <- long_df %>%
    dplyr::group_by(group, cell_type) %>%
    dplyr::summarise(mean_prop = mean(proportion), .groups = "drop")
  stack_plot <- ggplot2::ggplot(mean_prop, ggplot2::aes(x = group, y = mean_prop, fill = cell_type)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(22)) +
    ggplot2::labs(
      title = "Mean immune cell composition by group",
      x = "Group",
      y = "Mean proportion",
      fill = "LM22 cell type"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "right")
  ggplot2::ggsave(file.path(fig_dir, "deconvolution_stacked_barplot.pdf"), stack_plot, width = 11, height = 7)

  key_cells <- c("Tregs", "T_cells_CD8", "NK_cells_resting", "Macrophages_M1", "Macrophages_M2")
  violin_plot <- long_df %>%
    dplyr::filter(cell_type %in% key_cells) %>%
    ggplot2::ggplot(ggplot2::aes(x = group, y = proportion, fill = group)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.85) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
    ggplot2::facet_wrap(~ cell_type, scales = "free_y") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(
      title = "Key immune cell differences across groups",
      x = "Group",
      y = "Cell proportion",
      fill = "Group"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(file.path(fig_dir, "deconvolution_violin_key_cells.pdf"), violin_plot, width = 12, height = 7)

  corr_mat <- cor(t(proportion_matrix), method = "spearman")
  write.csv(corr_mat, file.path(res_dir, "cibersortx_celltype_correlation.csv"))
  pdf(file.path(fig_dir, "deconvolution_celltype_correlation_heatmap.pdf"), width = 10, height = 9)
  ComplexHeatmap::Heatmap(
    corr_mat,
    name = "rho",
    col = circlize::colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_title = "Cell type correlation matrix",
    row_title = "Cell type"
  )
  dev.off()

  heatmap_matrix <- proportion_matrix[, metadata$sample_id]
  column_ha <- ComplexHeatmap::HeatmapAnnotation(
    Group = metadata$group,
    DAS28 = metadata$DAS28,
    CRP = metadata$CRP,
    col = list(
      Group = c(RA_active = "#D73027", RA_remission = "#FDAE61", SLE = "#4575B4", HC = "#66BD63"),
      DAS28 = circlize::colorRamp2(c(min(metadata$DAS28), max(metadata$DAS28)), c("grey95", "black")),
      CRP = circlize::colorRamp2(c(min(metadata$CRP), max(metadata$CRP)), c("white", "darkred"))
    )
  )
  pdf(file.path(fig_dir, "deconvolution_abundance_heatmap.pdf"), width = 12, height = 8)
  ComplexHeatmap::Heatmap(
    heatmap_matrix,
    name = "Proportion",
    top_annotation = column_ha,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = FALSE,
    column_title = "Immune cell abundance by sample",
    row_title = "LM22 cell type"
  )
  dev.off()

  message("[02] Immune deconvolution analysis completed.")
  invisible(list(metadata = metadata, proportions = proportion_matrix, results = results_df))
}

if (sys.nframe() == 0L) {
  run_cibersortx_deconvolution()
}
