#!/usr/bin/env Rscript

set.seed(42)

run_integration_pipeline <- function(base_dir = NULL) {
  message("[07] Starting master integration pipeline...")

  suppressPackageStartupMessages({
    pkgs <- c("ggplot2", "dplyr", "tidyr", "gridExtra", "reshape2")
    missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      stop(sprintf("Missing required packages: %s", paste(missing_pkgs, collapse = ", ")))
    }
    invisible(lapply(pkgs, library, character.only = TRUE))
  })

  resolve_base_dir <- function(path_hint = NULL) {
    if (!is.null(path_hint)) return(normalizePath(path_hint, mustWork = TRUE))
    wd <- normalizePath(getwd(), mustWork = TRUE)
    if (basename(wd) == "scripts") dirname(wd) else wd
  }

  base_dir <- resolve_base_dir(base_dir)
  fig_dir <- file.path(base_dir, "figures")
  res_dir <- file.path(base_dir, "results")
  script_dir <- file.path(base_dir, "scripts")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

  source(file.path(script_dir, "01_multiomics_integration.R"))
  source(file.path(script_dir, "02_cibersortx_deconvolution.R"))
  source(file.path(script_dir, "03_cytokine_ode_model.R"))
  source(file.path(script_dir, "04_singlecell_checkpoint.R"))
  source(file.path(script_dir, "05_treatment_response_prediction.R"))
  source(file.path(script_dir, "06_tolerance_insilico.R"))

  multi_res <- run_multiomics_integration(base_dir)
  deconv_res <- run_cibersortx_deconvolution(base_dir)
  ode_res <- run_cytokine_ode_model(base_dir)
  sc_res <- run_singlecell_checkpoint(base_dir)
  pred_res <- run_treatment_response_prediction(base_dir)
  tol_res <- run_tolerance_insilico(base_dir)

  variance_df <- read.csv(file.path(res_dir, "mofa_variance_explained.csv"))
  cibersort_df <- read.csv(file.path(res_dir, "cibersortx_sample_proportions.csv"))
  ode_df <- read.csv(file.path(res_dir, "ode_steady_states.csv"))
  checkpoint_freq <- read.csv(file.path(res_dir, "checkpoint_exhausted_tcell_frequency.csv"))
  model_perf <- read.csv(file.path(res_dir, "model_performance.csv"))
  tol_eff <- read.csv(file.path(res_dir, "tolerance_efficacy.csv"))
  checkpoint_summary <- read.csv(file.path(res_dir, "checkpoint_expression_summary.csv"))

  p1 <- variance_df %>%
    dplyr::filter(factor %in% paste0("Factor", 1:5)) %>%
    ggplot2::ggplot(ggplot2::aes(x = factor, y = r2, fill = view)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::labs(title = "MOFA variance", x = "Factor", y = "R2", fill = "View")

  mean_cib <- cibersort_df %>%
    tidyr::pivot_longer(cols = -sample_id, names_to = "cell_type", values_to = "prop") %>%
    dplyr::mutate(group = rep(rep(c("RA_active", "RA_remission", "SLE", "HC"), each = 20), each = ncol(cibersort_df) - 1)) %>%
    dplyr::filter(cell_type %in% c("Tregs", "T_cells_CD8", "NK_cells_resting", "Macrophages_M1")) %>%
    dplyr::group_by(group, cell_type) %>%
    dplyr::summarise(prop = mean(prop), .groups = "drop")
  p2 <- ggplot2::ggplot(mean_cib, ggplot2::aes(x = group, y = prop, fill = cell_type)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::labs(title = "Deconvolution signatures", x = "Group", y = "Mean proportion", fill = "Cell type")

  p3 <- ggplot2::ggplot(ode_df, ggplot2::aes(x = inflammatory_score, y = regulatory_score, color = scenario)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::labs(title = "Cytokine steady states", x = "Inflammatory score", y = "Regulatory score", color = "Scenario")

  p4 <- ggplot2::ggplot(checkpoint_freq, ggplot2::aes(x = cell_type, y = exhausted_frequency, fill = condition)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::labs(title = "Exhausted T-cell frequency", x = "Cell type", y = "Frequency", fill = "Condition")

  p5 <- model_perf %>%
    dplyr::filter(subset == "Overall") %>%
    ggplot2::ggplot(ggplot2::aes(x = model, y = AUC, fill = model)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::coord_cartesian(ylim = c(0.5, 1.0)) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::labs(title = "Treatment response AUC", x = "Model", y = "AUC")

  p6 <- ggplot2::ggplot(tol_eff, ggplot2::aes(x = reorder(strategy, composite_rank_score), y = composite_rank_score, fill = tolerance_score)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::labs(title = "Tolerance strategies", x = "Strategy", y = "Composite score", fill = "Tolerance")

  pdf(file.path(fig_dir, "master_integration_figure.pdf"), width = 16, height = 10)
  gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3)
  dev.off()

  auc_by_arm <- model_perf %>% dplyr::filter(model == "XGBoost", subset != "Overall") %>% dplyr::select(treatment_arm = subset, AUC)
  cytokine_map <- data.frame(
    treatment_arm = c("MTX", "Anti-TNF", "JAK inhibitor", "IL-6R inhibitor"),
    scenario = c("RA_active", "Anti_TNF_treatment", "JAK_inhibitor", "RA_active")
  )
  treatment_vs_cytokine <- dplyr::left_join(auc_by_arm, cytokine_map, by = "treatment_arm") %>%
    dplyr::left_join(ode_df[, c("scenario", "inflammatory_score", "regulatory_score", "inflammatory_ratio")], by = "scenario")
  cor1 <- cor(treatment_vs_cytokine$AUC, treatment_vs_cytokine$inflammatory_ratio, method = "spearman")

  deconv_group <- cibersort_df %>%
    tidyr::pivot_longer(cols = -sample_id, names_to = "cell_type", values_to = "prop") %>%
    dplyr::mutate(group = rep(rep(c("RA_active", "RA_remission", "SLE", "HC"), each = 20), each = ncol(cibersort_df) - 1)) %>%
    dplyr::group_by(group, cell_type) %>%
    dplyr::summarise(prop = mean(prop), .groups = "drop")
  checkpoint_group <- checkpoint_summary %>%
    dplyr::mutate(group = ifelse(condition == "RA_active", "RA_active", "HC"))
  tregs_foxp3 <- dplyr::left_join(
    deconv_group %>% dplyr::filter(cell_type == "Tregs") %>% dplyr::select(group, Treg_prop = prop),
    checkpoint_group %>% dplyr::filter(cell_type == "CD4 T") %>% dplyr::select(group, FOXP3_mean),
    by = "group"
  )
  cd8_pd1 <- dplyr::left_join(
    deconv_group %>% dplyr::filter(cell_type == "T_cells_CD8") %>% dplyr::select(group, CD8_prop = prop),
    checkpoint_group %>% dplyr::filter(cell_type == "CD8 T") %>% dplyr::select(group, PDCD1_mean),
    by = "group"
  )
  cor2 <- cor(cd8_pd1$CD8_prop, cd8_pd1$PDCD1_mean, method = "spearman")

  correlation_df <- data.frame(
    analysis_pair = c("Treatment AUC vs cytokine inflammatory ratio", "CD8 proportion vs PDCD1 expression", "Treg proportion vs FOXP3 expression"),
    metric = c(cor1, cor2, cor(tregs_foxp3$Treg_prop, tregs_foxp3$FOXP3_mean, method = "spearman")),
    stringsAsFactors = FALSE
  )
  write.csv(correlation_df, file.path(res_dir, "integration_correlations.csv"), row.names = FALSE)

  system_summary <- data.frame(
    metric = c(
      "Number of MOFA factors", "Best MOFA factor transcriptome R2", "RA-active exhausted T-cell frequency",
      "Best overall predictive AUC", "Most effective tolerance strategy score", "Spearman correlation AUC vs inflammatory ratio"
    ),
    value = c(
      10,
      max(variance_df$r2[variance_df$view == "transcriptome"]),
      max(checkpoint_freq$exhausted_frequency[checkpoint_freq$condition == "RA_active"]),
      max(model_perf$AUC[model_perf$subset == "Overall"], na.rm = TRUE),
      max(tol_eff$composite_rank_score),
      cor1
    ),
    stringsAsFactors = FALSE
  )
  write.csv(system_summary, file.path(res_dir, "system_summary.csv"), row.names = FALSE)

  message("[07] Master integration pipeline completed.")
  invisible(list(
    multiomics = multi_res,
    deconvolution = deconv_res,
    ode = ode_res,
    scrna = sc_res,
    prediction = pred_res,
    tolerance = tol_res,
    summary = system_summary
  ))
}

if (sys.nframe() == 0L) {
  run_integration_pipeline()
}
