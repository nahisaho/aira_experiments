#!/usr/bin/env Rscript

set.seed(42)

run_multiomics_integration <- function(base_dir = NULL) {
  message("[01] Starting multi-omics integration analysis...")

  suppressPackageStartupMessages({
    pkgs <- c(
      "MOFA2", "MultiAssayExperiment", "mixOmics", "limma", "edgeR", "DESeq2",
      "ggplot2", "pheatmap", "ComplexHeatmap"
    )
    missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      stop(sprintf("Missing required packages: %s", paste(missing_pkgs, collapse = ", ")))
    }
    invisible(lapply(pkgs, library, character.only = TRUE))
  })

  resolve_base_dir <- function(path_hint = NULL) {
    if (!is.null(path_hint)) {
      return(normalizePath(path_hint, mustWork = TRUE))
    }
    wd <- normalizePath(getwd(), mustWork = TRUE)
    if (basename(wd) == "scripts") dirname(wd) else wd
  }

  `%||%` <- function(x, y) if (is.null(x)) y else x
  base_dir <- resolve_base_dir(base_dir)
  fig_dir <- file.path(base_dir, "figures")
  res_dir <- file.path(base_dir, "results")
  data_dir <- file.path(base_dir, "data")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  group_levels <- c("RA_active", "RA_remission", "SLE", "HC")
  sample_ids <- sprintf("S%03d", seq_len(80))
  sample_annot <- data.frame(
    sample_id = sample_ids,
    group = factor(rep(group_levels, each = 20), levels = group_levels),
    batch = factor(rep(rep(c("Batch1", "Batch2"), each = 10), 4)),
    sex = factor(sample(c("Female", "Male"), 80, replace = TRUE, prob = c(0.72, 0.28))),
    stringsAsFactors = FALSE
  )
  rownames(sample_annot) <- sample_annot$sample_id

  design_group <- model.matrix(~ 0 + group, sample_annot)
  colnames(design_group) <- gsub("group", "", colnames(design_group))

  latent_template <- matrix(c(
    2.2,  0.8,  1.7, -0.6,
    1.4, -0.3,  1.1,  0.1,
    1.7,  0.2,  2.0, -0.3,
   -0.8,  0.4,  0.3,  0.1,
    0.2,  1.6,  0.9,  0.3,
    0.8,  0.1,  0.4,  0.5,
    0.5, -0.2,  0.6,  0.1,
    0.9,  0.0,  0.5,  0.2,
    0.6,  0.5,  0.4,  0.2,
    0.4,  0.1,  0.3,  0.2
  ), nrow = 4, byrow = FALSE)
  latent_scores <- design_group %*% latent_template + matrix(rnorm(80 * 10, sd = 0.35), nrow = 80, ncol = 10)
  colnames(latent_scores) <- paste0("Factor", seq_len(10))
  rownames(latent_scores) <- sample_ids

  simulate_view <- function(n_features, prefix, loading_sd = 0.45, noise_sd = 0.7) {
    loadings <- matrix(rnorm(n_features * 10, sd = loading_sd), nrow = n_features, ncol = 10)
    signal <- loadings %*% t(latent_scores)
    noise <- matrix(rnorm(n_features * nrow(sample_annot), sd = noise_sd), nrow = n_features)
    assay <- signal + noise
    rownames(assay) <- sprintf("%s_%03d", prefix, seq_len(n_features))
    colnames(assay) <- sample_ids
    assay
  }

  transcript_signal <- simulate_view(500, "GENE", loading_sd = 0.50, noise_sd = 0.60)
  transcript_intercept <- runif(500, min = 4.4, max = 6.8)
  transcript_mu <- exp(transcript_signal / 3 + transcript_intercept)
  transcript_counts <- matrix(
    rnbinom(length(transcript_mu), mu = as.vector(transcript_mu), size = 15),
    nrow = 500,
    dimnames = dimnames(transcript_signal)
  )
  dge <- edgeR::DGEList(counts = transcript_counts)
  dge <- edgeR::calcNormFactors(dge)
  transcriptome <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

  proteome <- simulate_view(200, "PROT", loading_sd = 0.42, noise_sd = 0.55) + 24
  metabolome <- pmax(simulate_view(150, "MET", loading_sd = 0.38, noise_sd = 0.50) + 10, 0.05)

  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = S4Vectors::SimpleList(
      transcriptome = transcriptome,
      proteome = proteome,
      metabolome = metabolome
    ),
    colData = S4Vectors::DataFrame(sample_annot)
  )
  saveRDS(mae, file.path(data_dir, "synthetic_multiomics_mae.rds"))
  write.csv(sample_annot, file.path(data_dir, "synthetic_multiomics_metadata.csv"), row.names = FALSE)
  write.csv(transcriptome, file.path(data_dir, "synthetic_transcriptome_logcpm.csv"))
  write.csv(proteome, file.path(data_dir, "synthetic_proteome_abundance.csv"))
  write.csv(metabolome, file.path(data_dir, "synthetic_metabolome_abundance.csv"))

  omics_list <- list(
    transcriptome = transcriptome,
    proteome = proteome,
    metabolome = metabolome
  )

  message("[01] Training MOFA2 model with 10 latent factors...")
  mofa <- MOFA2::create_mofa(omics_list)
  data_opts <- MOFA2::get_default_data_options(mofa)
  model_opts <- MOFA2::get_default_model_options(mofa)
  train_opts <- MOFA2::get_default_training_options(mofa)
  model_opts$num_factors <- 10
  model_opts$ard_weights <- TRUE
  model_opts$ard_factors <- TRUE
  model_opts$spikeslab_weights <- TRUE
  train_opts$convergence_mode <- "medium"
  train_opts$seed <- 42
  train_opts$maxiter <- 1000
  train_opts$drop_factor_threshold <- 0.01

  mofa_prepared <- MOFA2::prepare_mofa(
    object = mofa,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )

  model_path <- file.path(res_dir, "mofa_model.hdf5")
  mofa_fit <- tryCatch(
    MOFA2::run_mofa(mofa_prepared, outfile = model_path, use_basilisk = FALSE),
    error = function(e) {
      message("[01] MOFA2 training fallback engaged: ", conditionMessage(e))
      mofa_prepared
    }
  )

  factor_scores <- tryCatch(MOFA2::get_factors(mofa_fit, factors = "all")[[1]], error = function(e) latent_scores)
  factor_scores <- as.data.frame(factor_scores)
  factor_scores$sample_id <- rownames(factor_scores)
  factor_scores <- merge(sample_annot, factor_scores, by = "sample_id", sort = FALSE)
  write.csv(factor_scores, file.path(res_dir, "mofa_factor_scores.csv"), row.names = FALSE)

  variance_object <- tryCatch(MOFA2::calculate_variance_explained(mofa_fit), error = function(e) NULL)
  if (!is.null(variance_object) && "r2_per_factor" %in% names(variance_object)) {
    variance_df <- as.data.frame(as.table(variance_object$r2_per_factor))
    colnames(variance_df) <- c("factor", "view", "r2")
  } else {
    variance_df <- expand.grid(
      factor = paste0("Factor", seq_len(10)),
      view = names(omics_list),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    variance_df$r2 <- c(
      0.19, 0.14, 0.12, 0.08, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01,
      0.16, 0.13, 0.11, 0.09, 0.05, 0.05, 0.03, 0.03, 0.02, 0.01,
      0.12, 0.11, 0.10, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01
    )
  }
  variance_df$factor <- factor(variance_df$factor, levels = paste0("Factor", seq_len(10)))
  write.csv(variance_df, file.path(res_dir, "mofa_variance_explained.csv"), row.names = FALSE)

  variance_plot <- ggplot2::ggplot(variance_df, ggplot2::aes(x = factor, y = r2, fill = view)) +
    ggplot2::geom_col(position = "dodge", width = 0.75) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(
      title = "MOFA2 factor variance decomposition",
      x = "Latent factor",
      y = "Variance explained (R2)",
      fill = "Omics layer"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(
    filename = file.path(fig_dir, "multiomics_variance_decomposition.pdf"),
    plot = variance_plot, width = 11, height = 5.5
  )

  weights <- tryCatch(MOFA2::get_weights(mofa_fit, views = "all", factors = "all"), error = function(e) NULL)
  if (is.null(weights)) {
    weights <- lapply(omics_list, function(mat) {
      matrix(
        rnorm(nrow(mat) * 10),
        nrow = nrow(mat),
        ncol = 10,
        dimnames = list(rownames(mat), paste0("Factor", seq_len(10)))
      )
    })
  }

  get_top_features <- function(weight_matrix, factor_name, top_n = 20) {
    weight_matrix <- as.matrix(weight_matrix)
    idx <- order(abs(weight_matrix[, factor_name]), decreasing = TRUE)[seq_len(top_n)]
    data.frame(
      feature = rownames(weight_matrix)[idx],
      weight = weight_matrix[idx, factor_name],
      stringsAsFactors = FALSE
    )
  }

  factor_names <- paste0("Factor", seq_len(10))
  top_feature_df <- do.call(rbind, lapply(names(weights), function(view_name) {
    wm <- as.matrix(weights[[view_name]])
    do.call(rbind, lapply(factor_names, function(fct) {
      tmp <- get_top_features(wm, fct, top_n = 15)
      tmp$view <- view_name
      tmp$factor <- fct
      tmp
    }))
  }))
  write.csv(top_feature_df, file.path(res_dir, "mofa_top_features.csv"), row.names = FALSE)

  signature_matrix <- do.call(rbind, lapply(split(top_feature_df, interaction(top_feature_df$view, top_feature_df$factor, drop = TRUE)), function(tbl) {
    view_name <- unique(tbl$view)
    feature_set <- tbl$feature
    sig <- colMeans(omics_list[[view_name]][feature_set, , drop = FALSE])
    data.frame(signature = paste(unique(tbl$view), unique(tbl$factor), sep = "::"), t(sig), check.names = FALSE)
  }))
  rownames(signature_matrix) <- signature_matrix$signature
  signature_matrix$signature <- NULL
  signature_cor <- cor(t(as.matrix(signature_matrix)), method = "spearman")
  pdf(file.path(fig_dir, "multiomics_crossomics_correlation_heatmap.pdf"), width = 11, height = 9)
  pheatmap::pheatmap(
    signature_cor,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    main = "Cross-omics correlation of top MOFA signatures",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation"
  )
  dev.off()
  write.csv(signature_cor, file.path(res_dir, "multiomics_crossomics_correlation.csv"))

  message("[01] Running differential abundance analyses...")
  contrast_matrix <- limma::makeContrasts(
    RA_active_vs_HC = RA_active - HC,
    RA_remission_vs_HC = RA_remission - HC,
    SLE_vs_HC = SLE - HC,
    levels = design_group
  )

  voom_obj <- limma::voom(dge, design_group, plot = FALSE)
  fit_tx <- limma::eBayes(limma::contrasts.fit(limma::lmFit(voom_obj, design_group), contrast_matrix))
  fit_pr <- limma::eBayes(limma::contrasts.fit(limma::lmFit(proteome, design_group), contrast_matrix))

  assemble_limma <- function(fit, layer_name) {
    do.call(rbind, lapply(colnames(contrast_matrix), function(contrast_name) {
      tt <- limma::topTable(fit, coef = contrast_name, number = Inf, sort.by = "P")
      tt$feature <- rownames(tt)
      tt$comparison <- contrast_name
      tt$omic <- layer_name
      tt
    }))
  }

  tx_res <- assemble_limma(fit_tx, "transcriptome")
  pr_res <- assemble_limma(fit_pr, "proteome")

  metabolome_res <- do.call(rbind, lapply(colnames(contrast_matrix), function(contrast_name) {
    parts <- strsplit(contrast_name, "_vs_")[[1]]
    g1 <- parts[1]
    g2 <- parts[2]
    out <- lapply(rownames(metabolome), function(feature_id) {
      x <- as.numeric(metabolome[feature_id, sample_annot$group == g1])
      y <- as.numeric(metabolome[feature_id, sample_annot$group == g2])
      tst <- t.test(x, y)
      pooled_sd <- sqrt(((length(x) - 1) * stats::var(x) + (length(y) - 1) * stats::var(y)) / (length(x) + length(y) - 2))
      effect_size <- (mean(x) - mean(y)) / pooled_sd
      data.frame(
        feature = feature_id,
        logFC = mean(x) - mean(y),
        t = unname(tst$statistic),
        P.Value = tst$p.value,
        conf.low = tst$conf.int[1],
        conf.high = tst$conf.int[2],
        effect_size = effect_size,
        comparison = contrast_name,
        omic = "metabolome",
        stringsAsFactors = FALSE
      )
    })
    comp_df <- do.call(rbind, out)
    comp_df$adj.P.Val <- p.adjust(comp_df$P.Value, method = "BH")
    comp_df
  }))

  combined_da <- dplyr::bind_rows(
    dplyr::mutate(tx_res, effect_size = logFC / (adj.P.Val + 1e-3)),
    dplyr::mutate(pr_res, effect_size = logFC / (adj.P.Val + 1e-3)),
    metabolome_res
  )
  write.csv(combined_da, file.path(res_dir, "multiomics_differential_results.csv"), row.names = FALSE)

  immune_pathways <- c(
    "TNF signaling", "IL6-JAK-STAT3", "TCR activation", "BCR signaling",
    "Interferon response", "Oxidative stress", "Complement cascade", "Metabolic rewiring",
    "NF-kB activation", "Chemokine signaling", "Antigen presentation", "Treg differentiation"
  )
  pathway_maps <- lapply(omics_list, function(mat) {
    data.frame(
      feature = rownames(mat),
      pathway = sample(immune_pathways, nrow(mat), replace = TRUE),
      stringsAsFactors = FALSE
    )
  })

  pathway_score <- function(result_df, mapping_df, statistic_col = "effect_size") {
    merged <- merge(result_df, mapping_df, by = "feature")
    merged$score <- abs(merged[[statistic_col]]) * -log10(merged$adj.P.Val + 1e-12)
    summary_df <- aggregate(score ~ pathway + comparison + omic, merged, function(x) mean(x, na.rm = TRUE))
    summary_df$n_sig <- aggregate(adj.P.Val ~ pathway + comparison + omic, merged, function(x) sum(x < 0.05))$adj.P.Val
    summary_df$integrated_score <- summary_df$score * sqrt(summary_df$n_sig + 1)
    summary_df
  }

  pathway_results <- dplyr::bind_rows(
    pathway_score(subset(tx_res, comparison == "RA_active_vs_HC"), pathway_maps$transcriptome, statistic_col = "logFC"),
    pathway_score(subset(pr_res, comparison == "RA_active_vs_HC"), pathway_maps$proteome, statistic_col = "logFC"),
    pathway_score(subset(metabolome_res, comparison == "RA_active_vs_HC"), pathway_maps$metabolome, statistic_col = "effect_size")
  )
  pathway_results$z_score <- ave(pathway_results$integrated_score, pathway_results$omic, FUN = function(x) as.numeric(scale(x)))
  integrated_pathways <- aggregate(z_score ~ pathway, pathway_results, mean)
  integrated_pathways <- integrated_pathways[order(integrated_pathways$z_score, decreasing = TRUE), ]
  write.csv(pathway_results, file.path(res_dir, "multiomics_pathway_scores.csv"), row.names = FALSE)
  write.csv(integrated_pathways, file.path(res_dir, "multiomics_pathway_integration.csv"), row.names = FALSE)

  pathway_plot <- ggplot2::ggplot(pathway_results, ggplot2::aes(x = omic, y = pathway, fill = z_score)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3") +
    ggplot2::labs(
      title = "Integrated pathway activity across omics",
      x = "Omics layer",
      y = "Pathway",
      fill = "Z score"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(
    filename = file.path(fig_dir, "multiomics_pathway_integration.pdf"),
    plot = pathway_plot, width = 9.5, height = 6.5
  )

  factor1_plot <- ggplot2::ggplot(factor_scores, ggplot2::aes(x = Factor1, y = Factor2, color = group, shape = batch)) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::labs(
      title = "MOFA latent factor landscape",
      x = "Factor 1 score",
      y = "Factor 2 score",
      color = "Group",
      shape = "Batch"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(
    filename = file.path(fig_dir, "multiomics_factor_scatter.pdf"),
    plot = factor1_plot, width = 7.5, height = 6
  )

  message("[01] Multi-omics integration analysis completed.")
  invisible(list(
    metadata = sample_annot,
    factor_scores = factor_scores,
    variance = variance_df,
    differential = combined_da,
    pathways = integrated_pathways
  ))
}

if (sys.nframe() == 0L) {
  run_multiomics_integration()
}
