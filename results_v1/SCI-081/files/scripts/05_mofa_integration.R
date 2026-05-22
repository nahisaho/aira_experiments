#!/usr/bin/env Rscript
# =============================================================================
# Module 5: MOFA+ Multi-Omics Factor Analysis & Patient Stratification
# =============================================================================
# Input:  mRNA, Proteome, Phosphoproteome, CNA matrices + clinical metadata
# Output: MOFA model, factor loadings, patient clusters, survival analysis
# =============================================================================

suppressPackageStartupMessages({
  library(MOFA2)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(survival)
  library(survminer)
  library(pheatmap)
  library(yaml)
})

config <- yaml::read_yaml("config/pipeline_config.yaml")
cfg    <- config$mofa_plus

cat("=== Module 5: MOFA+ Multi-Omics Factor Analysis ===\n")

# -----------------------------------------------------------------------------
# Step 1: Prepare multi-omics data views
# -----------------------------------------------------------------------------
prepare_mofa_input <- function() {
  views <- list()

  for (v in cfg$views) {
    cat(sprintf("  Loading view: %s (%s)\n", v$name, v$file))

    mat <- read.delim(v$file, row.names = 1, check.names = FALSE)
    mat <- as.matrix(mat)

    # Handle MaxQuant proteinGroups format
    if (v$name %in% c("Proteome", "Phosphoproteome")) {
      raw <- read.delim(v$file, stringsAsFactors = FALSE)
      if ("LFQ.intensity." %in% substr(colnames(raw), 1, 15)[1]) {
        int_cols <- grep("^LFQ\\.intensity\\.", colnames(raw), value = TRUE)
        mat <- as.matrix(raw[, int_cols])
        rownames(mat) <- make.unique(
          sapply(strsplit(raw$Gene.names, ";"), `[`, 1))
        colnames(mat) <- gsub("^LFQ\\.intensity\\.", "", colnames(mat))
      }
      mat[mat == 0] <- NA
      mat <- log2(mat)
    }

    # Select top variable features
    n_feat <- v$features
    if (is.numeric(n_feat) && n_feat < nrow(mat)) {
      vars <- apply(mat, 1, var, na.rm = TRUE)
      top_idx <- order(vars, decreasing = TRUE)[seq_len(n_feat)]
      mat <- mat[top_idx, ]
    }

    # Z-score normalization
    mat <- t(scale(t(mat)))

    views[[v$name]] <- mat
    cat(sprintf("    Dimensions: %d features × %d samples\n", nrow(mat), ncol(mat)))
  }

  # Harmonize sample names across views
  all_samples <- Reduce(intersect, lapply(views, colnames))
  cat(sprintf("  Common samples across all views: %d\n", length(all_samples)))

  for (nm in names(views)) {
    views[[nm]] <- views[[nm]][, all_samples, drop = FALSE]
  }

  return(list(views = views, samples = all_samples))
}

# -----------------------------------------------------------------------------
# Step 2: Create and train MOFA model
# -----------------------------------------------------------------------------
train_mofa <- function(views, n_factors = 15, seed = 42) {
  cat("\n--- Training MOFA+ ---\n")

  mofa_obj <- create_mofa(views)

  # Set data options
  data_opts <- get_default_data_options(mofa_obj)
  data_opts$scale_views <- TRUE
  data_opts$scale_groups <- FALSE

  # Set model options
  model_opts <- get_default_model_options(mofa_obj)
  model_opts$num_factors <- n_factors
  model_opts$likelihoods <- rep("gaussian", length(views))
  names(model_opts$likelihoods) <- names(views)

  # Set training options
  train_opts <- get_default_training_options(mofa_obj)
  train_opts$convergence_mode <- cfg$convergence_mode
  train_opts$drop_factor_threshold <- cfg$drop_factor_threshold
  train_opts$seed <- seed
  train_opts$gpu_mode <- cfg$gpu
  train_opts$verbose <- TRUE

  mofa_obj <- prepare_mofa(mofa_obj,
    data_options    = data_opts,
    model_options   = model_opts,
    training_options = train_opts
  )

  mofa_trained <- run_mofa(mofa_obj,
    outfile = cfg$output_model,
    use_basilisk = TRUE
  )

  cat(sprintf("  Trained model with %d factors\n",
    get_dimensions(mofa_trained)$K))
  cat(sprintf("  Model saved: %s\n", cfg$output_model))

  return(mofa_trained)
}

# -----------------------------------------------------------------------------
# Step 3: Variance decomposition analysis
# -----------------------------------------------------------------------------
analyze_variance <- function(mofa_model, output_dir = "results/") {
  cat("\n--- Variance Decomposition ---\n")

  r2 <- calculate_variance_explained(mofa_model)

  # Per-view variance explained
  r2_total <- r2$r2_total
  for (view in names(r2_total[[1]])) {
    cat(sprintf("  %s: total R² = %.1f%%\n",
      view, r2_total[[1]][[view]] * 100))
  }

  # Per-factor variance explained
  r2_per_factor <- r2$r2_per_factor[[1]]

  # Save
  write.csv(r2_per_factor,
    file.path(output_dir, "mofa_variance_per_factor.csv"))

  return(r2)
}

# -----------------------------------------------------------------------------
# Step 4: Factor interpretation & patient clustering
# -----------------------------------------------------------------------------
cluster_patients <- function(mofa_model, metadata, n_clusters = 3) {
  cat("\n--- Patient Clustering ---\n")

  # Extract factor values
  Z <- get_factors(mofa_model)[[1]]
  cat(sprintf("  Factor matrix: %d patients × %d factors\n", nrow(Z), ncol(Z)))

  # Consensus clustering on factor values
  set.seed(42)
  km <- kmeans(Z, centers = n_clusters, nstart = 25, iter.max = 100)
  clusters <- km$cluster

  cluster_df <- data.frame(
    sample  = names(clusters),
    cluster = paste0("C", clusters),
    stringsAsFactors = FALSE
  )

  # Merge with clinical metadata
  if (!is.null(metadata)) {
    cluster_df <- merge(cluster_df, metadata, by.x = "sample", by.y = "sample_id",
      all.x = TRUE)
  }

  cat(sprintf("  Cluster sizes: %s\n",
    paste(table(cluster_df$cluster), collapse = ", ")))

  write.csv(cluster_df, "results/mofa_patient_clusters.csv", row.names = FALSE)
  return(cluster_df)
}

# -----------------------------------------------------------------------------
# Step 5: Survival analysis by cluster
# -----------------------------------------------------------------------------
survival_analysis <- function(cluster_df, output_dir = "figures/") {
  cat("\n--- Survival Analysis ---\n")

  if (!all(c("os_time", "os_status") %in% colnames(cluster_df))) {
    cat("  Missing survival data — skipping\n")
    return(NULL)
  }

  surv_obj <- Surv(time = cluster_df$os_time, event = cluster_df$os_status)

  # Log-rank test
  diff <- survdiff(surv_obj ~ cluster, data = cluster_df)
  p_val <- 1 - pchisq(diff$chisq, length(unique(cluster_df$cluster)) - 1)
  cat(sprintf("  Log-rank p-value: %.2e\n", p_val))

  # Kaplan-Meier plot
  fit <- survfit(surv_obj ~ cluster, data = cluster_df)
  p <- ggsurvplot(fit, data = cluster_df,
    pval = TRUE,
    risk.table = TRUE,
    palette = c("#2166ac", "#b2182b", "#1b7837"),
    title = "Overall Survival by MOFA Cluster",
    xlab = "Time (months)",
    ylab = "Survival probability",
    legend.labs = paste0("Cluster ", sort(unique(cluster_df$cluster)))
  )

  pdf(file.path(output_dir, "mofa_survival_km.pdf"), width = 8, height = 7)
  print(p)
  dev.off()
  cat("  KM plot saved\n")

  return(list(fit = fit, pval = p_val))
}

# -----------------------------------------------------------------------------
# Step 6: Factor–clinical association
# -----------------------------------------------------------------------------
factor_clinical_association <- function(mofa_model, metadata, output_dir = "results/") {
  cat("\n--- Factor–Clinical Association ---\n")

  Z <- get_factors(mofa_model)[[1]]

  clinical_vars <- c("molecular_subtype", "tumor_grade", "stage")
  available_vars <- intersect(clinical_vars, colnames(metadata))

  assoc_results <- data.frame()
  for (cvar in available_vars) {
    for (f in seq_len(ncol(Z))) {
      factor_vals <- Z[, f]
      clin_vals   <- metadata[[cvar]][match(rownames(Z), metadata$sample_id)]

      if (is.numeric(clin_vals)) {
        ct <- cor.test(factor_vals, clin_vals, method = "spearman")
        assoc_results <- rbind(assoc_results, data.frame(
          factor = paste0("Factor", f), variable = cvar,
          statistic = ct$estimate, p_value = ct$p.value,
          test = "spearman"))
      } else {
        kw <- kruskal.test(factor_vals ~ as.factor(clin_vals))
        assoc_results <- rbind(assoc_results, data.frame(
          factor = paste0("Factor", f), variable = cvar,
          statistic = kw$statistic, p_value = kw$p.value,
          test = "kruskal"))
      }
    }
  }

  assoc_results$padj <- p.adjust(assoc_results$p_value, method = "BH")
  write.csv(assoc_results, file.path(output_dir, "mofa_factor_clinical_assoc.csv"),
    row.names = FALSE)

  sig <- sum(assoc_results$padj < 0.05)
  cat(sprintf("  Significant factor-clinical associations: %d\n", sig))
  return(assoc_results)
}

# -----------------------------------------------------------------------------
# Step 7: Visualization
# -----------------------------------------------------------------------------
plot_mofa <- function(mofa_model, r2, cluster_df, output_dir = "figures/") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # (A) Variance explained heatmap
  pdf(file.path(output_dir, "mofa_variance_heatmap.pdf"), width = 10, height = 6)
  plot_variance_explained(mofa_model, plot_total = TRUE)
  dev.off()

  # (B) Factor scatter (Factor1 vs Factor2)
  Z <- get_factors(mofa_model)[[1]]
  z_df <- as.data.frame(Z[, 1:min(ncol(Z), 3)])
  z_df$sample <- rownames(Z)
  z_df <- merge(z_df, cluster_df[, c("sample", "cluster")],
    by = "sample", all.x = TRUE)

  p <- ggplot(z_df, aes(x = Factor1, y = Factor2, color = cluster)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = c("#2166ac", "#b2182b", "#1b7837")) +
    labs(title = "MOFA Factor Space — Patient Clusters",
         x = "Factor 1", y = "Factor 2") +
    theme_minimal(base_size = 14)

  ggsave(file.path(output_dir, "mofa_factor_scatter.pdf"), p, width = 8, height = 6)

  # (C) Top weights per factor
  pdf(file.path(output_dir, "mofa_top_weights.pdf"), width = 10, height = 12)
  plot_top_weights(mofa_model, view = "Proteome", factor = 1, nfeatures = 20)
  dev.off()

  cat("  MOFA figures saved\n")
}

# =============================================================================
# Main
# =============================================================================
main <- function() {
  input <- prepare_mofa_input()

  mofa <- train_mofa(input$views, n_factors = cfg$n_factors, seed = cfg$seed)

  r2 <- analyze_variance(mofa)

  metadata <- read.delim(cfg$groups, stringsAsFactors = FALSE)
  clusters <- cluster_patients(mofa, metadata)

  surv <- survival_analysis(clusters)

  assoc <- factor_clinical_association(mofa, metadata)

  plot_mofa(mofa, r2, clusters)

  cat("\n=== Module 5 complete ===\n")
}

if (!interactive()) main()
