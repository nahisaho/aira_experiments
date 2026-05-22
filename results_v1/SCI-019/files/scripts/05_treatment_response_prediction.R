#!/usr/bin/env Rscript

set.seed(42)

run_treatment_response_prediction <- function(base_dir = NULL) {
  message("[05] Starting RA treatment response prediction analysis...")

  suppressPackageStartupMessages({
    pkgs <- c("caret", "randomForest", "glmnet", "xgboost", "pROC", "ggplot2", "reshape2", "ROCR", "tidyverse")
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
  data_dir <- file.path(base_dir, "data")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  n <- 200
  treatment_levels <- c("MTX", "Anti-TNF", "JAK inhibitor", "IL-6R inhibitor")
  patient_id <- sprintf("RA_%03d", seq_len(n))
  clinical_df <- data.frame(
    patient_id = patient_id,
    treatment_arm = factor(sample(treatment_levels, n, replace = TRUE), levels = treatment_levels),
    DAS28 = pmin(pmax(rnorm(n, 5.2, 1.0), 2.1), 7.5),
    CRP = rlnorm(n, 1.1, 0.45),
    ESR = rnorm(n, 38, 14),
    RF_status = factor(sample(c("Negative", "Positive"), n, replace = TRUE, prob = c(0.28, 0.72))),
    anti_CCP = factor(sample(c("Negative", "Positive"), n, replace = TRUE, prob = c(0.22, 0.78))),
    disease_duration = rgamma(n, shape = 2.5, scale = 3.2),
    age = round(rnorm(n, 55, 12)),
    sex = factor(sample(c("Female", "Male"), n, replace = TRUE, prob = c(0.74, 0.26))),
    stringsAsFactors = FALSE
  )

  gene_mat <- replicate(50, rnorm(n, 0, 1))
  colnames(gene_mat) <- sprintf("GeneExpr_%02d", seq_len(50))
  gene_mat[, 1:5] <- gene_mat[, 1:5] + ifelse(clinical_df$treatment_arm == "Anti-TNF", 0.6, 0)
  gene_mat[, 6:10] <- gene_mat[, 6:10] + scale(clinical_df$DAS28)
  gene_mat <- as.data.frame(gene_mat)

  immune_mat <- as.data.frame(replicate(10, rbeta(n, 2, 8)))
  colnames(immune_mat) <- c("CD8_T", "Treg", "NK", "M1", "M2", "Monocyte", "Bmem", "Plasma", "DC_act", "Neutrophil")
  immune_mat$CD8_T <- pmin(pmax(immune_mat$CD8_T + scale(clinical_df$CRP) * 0.04, 0.01), 0.60)
  immune_mat$Treg <- pmin(pmax(immune_mat$Treg - scale(clinical_df$DAS28) * 0.03, 0.01), 0.30)
  immune_mat$M1 <- pmin(pmax(immune_mat$M1 + scale(clinical_df$CRP) * 0.05, 0.01), 0.50)

  cytokine_df <- data.frame(
    IL6 = pmax(rnorm(n, 8, 3) + scale(clinical_df$CRP) * 2, 0.2),
    TNFa = pmax(rnorm(n, 6, 2) + scale(clinical_df$DAS28) * 1.8, 0.2),
    IL17A = pmax(rnorm(n, 4.5, 1.8) + immune_mat$CD8_T * 4, 0.1),
    IL10 = pmax(rnorm(n, 3.0, 1.0) + immune_mat$Treg * 5, 0.1)
  )

  linpred <- -0.35 +
    0.70 * (clinical_df$treatment_arm == "Anti-TNF") -
    0.55 * as.numeric(scale(clinical_df$DAS28)) -
    0.35 * as.numeric(scale(clinical_df$CRP)) +
    0.42 * as.numeric(scale(gene_mat$GeneExpr_01)) -
    0.28 * as.numeric(scale(gene_mat$GeneExpr_07)) +
    0.65 * as.numeric(scale(immune_mat$Treg)) -
    0.48 * as.numeric(scale(immune_mat$M1)) -
    0.45 * as.numeric(scale(cytokine_df$IL6)) -
    0.38 * as.numeric(scale(cytokine_df$TNFa)) +
    0.32 * as.numeric(scale(cytokine_df$IL10)) +
    0.30 * (clinical_df$treatment_arm == "JAK inhibitor") * as.numeric(scale(gene_mat$GeneExpr_12)) +
    0.26 * (clinical_df$treatment_arm == "IL-6R inhibitor") * as.numeric(scale(cytokine_df$IL6))
  prob <- plogis(linpred)
  response <- factor(ifelse(runif(n) < prob, "Responder", "NonResponder"), levels = c("Responder", "NonResponder"))
  message(sprintf("[05] Simulated responder rate: %.1f%%", mean(response == "Responder") * 100))

  modeling_df <- dplyr::bind_cols(clinical_df, gene_mat, immune_mat, cytokine_df, response = response)
  write.csv(modeling_df, file.path(data_dir, "synthetic_treatment_response_dataset.csv"), row.names = FALSE)

  model_df <- modeling_df
  model_df$response <- factor(model_df$response, levels = c("Responder", "NonResponder"))
  message("[05] Training LASSO, Random Forest, XGBoost, and SVM models with repeated cross-validation...")
  ctrl <- caret::trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 3,
    summaryFunction = caret::twoClassSummary,
    classProbs = TRUE,
    savePredictions = "final",
    allowParallel = FALSE
  )

  formula_obj <- response ~ . - patient_id
  set.seed(42)
  fit_glmnet <- caret::train(
    formula_obj, data = model_df, method = "glmnet", metric = "ROC", trControl = ctrl,
    tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 0.08, length.out = 12)), family = "binomial"
  )
  set.seed(42)
  fit_rf <- caret::train(
    formula_obj, data = model_df, method = "rf", metric = "ROC", trControl = ctrl,
    tuneGrid = expand.grid(mtry = c(6, 10, 14)), ntree = 500, importance = TRUE
  )
  set.seed(42)
  fit_xgb <- caret::train(
    formula_obj, data = model_df, method = "xgbTree", metric = "ROC", trControl = ctrl,
    tuneGrid = expand.grid(
      nrounds = c(75, 125), max_depth = c(3, 5), eta = c(0.05, 0.10),
      gamma = 0, colsample_bytree = 0.8, min_child_weight = 1, subsample = 0.8
    ), verbose = FALSE
  )
  set.seed(42)
  fit_svm <- caret::train(
    formula_obj, data = model_df, method = "svmRadial", metric = "ROC", trControl = ctrl,
    preProcess = c("center", "scale"), tuneLength = 6
  )

  fits <- list(LASSO = fit_glmnet, RandomForest = fit_rf, XGBoost = fit_xgb, SVM = fit_svm)
  overall_perf <- dplyr::bind_rows(lapply(names(fits), function(name) {
    best_row <- fits[[name]]$results[which.max(fits[[name]]$results$ROC), ]
    data.frame(
      model = name,
      AUC = best_row$ROC,
      Sensitivity = best_row$Sens,
      Specificity = best_row$Spec,
      stringsAsFactors = FALSE
    )
  }))
  overall_perf$subset <- "Overall"

  get_best_predictions <- function(fit, model_name) {
    pred_df <- fit$pred
    if (!is.null(fit$bestTune)) {
      for (nm in names(fit$bestTune)) {
        pred_df <- pred_df[pred_df[[nm]] == fit$bestTune[[nm]], ]
      }
    }
    pred_df$model <- model_name
    pred_df$patient_id <- model_df$patient_id[match(pred_df$rowIndex, seq_len(nrow(model_df)))]
    pred_df$treatment_arm <- model_df$treatment_arm[match(pred_df$rowIndex, seq_len(nrow(model_df)))]
    pred_df
  }

  pred_df <- dplyr::bind_rows(
    get_best_predictions(fit_glmnet, "LASSO"),
    get_best_predictions(fit_rf, "RandomForest"),
    get_best_predictions(fit_xgb, "XGBoost"),
    get_best_predictions(fit_svm, "SVM")
  )
  write.csv(pred_df, file.path(res_dir, "treatment_predictions.csv"), row.names = FALSE)

  subgroup_perf <- pred_df %>%
    dplyr::group_by(model, treatment_arm) %>%
    dplyr::summarise(
      AUC = if (dplyr::n_distinct(obs) > 1) as.numeric(pROC::auc(pROC::roc(obs, Responder, levels = c("NonResponder", "Responder"), direction = "<", quiet = TRUE))) else NA_real_,
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Sensitivity = NA_real_, Specificity = NA_real_, subset = as.character(treatment_arm)) %>%
    dplyr::select(model, AUC, Sensitivity, Specificity, subset, n)

  model_performance <- dplyr::bind_rows(
    dplyr::mutate(overall_perf, n = n),
    subgroup_perf
  )
  write.csv(model_performance, file.path(res_dir, "model_performance.csv"), row.names = FALSE)

  roc_list <- lapply(split(pred_df, pred_df$model), function(df) pROC::roc(df$obs, df$Responder, levels = c("NonResponder", "Responder"), direction = "<", quiet = TRUE))
  pdf(file.path(fig_dir, "prediction_roc_curves.pdf"), width = 8, height = 6)
  plot(roc_list[[1]], col = "#1B9E77", lwd = 2, main = "Cross-validated ROC curves")
  cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")
  i <- 1
  for (nm in names(roc_list)) {
    if (i == 1) {
      plot(roc_list[[nm]], col = cols[i], lwd = 2, main = "Cross-validated ROC curves")
    } else {
      plot(roc_list[[nm]], add = TRUE, col = cols[i], lwd = 2)
    }
    i <- i + 1
  }
  legend("bottomright", legend = sprintf("%s (AUC=%.3f)", names(roc_list), vapply(roc_list, pROC::auc, numeric(1))), col = cols, lwd = 2, bty = "n")
  dev.off()

  importance_rf <- caret::varImp(fit_rf)$importance %>% tibble::rownames_to_column("feature") %>% dplyr::mutate(model = "RandomForest")
  importance_xgb <- caret::varImp(fit_xgb)$importance %>% tibble::rownames_to_column("feature") %>% dplyr::mutate(model = "XGBoost")
  importance_df <- dplyr::bind_rows(importance_rf, importance_xgb) %>%
    dplyr::rename(importance = Overall) %>%
    dplyr::group_by(model) %>%
    dplyr::slice_max(order_by = importance, n = 20) %>%
    dplyr::ungroup()
  write.csv(importance_df, file.path(res_dir, "feature_importance.csv"), row.names = FALSE)

  importance_plot <- ggplot2::ggplot(importance_df, ggplot2::aes(x = reorder(feature, importance), y = importance, fill = model)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ model, scales = "free_y") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(title = "Top predictive features", x = "Feature", y = "Importance") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(fig_dir, "prediction_feature_importance.pdf"), importance_plot, width = 11, height = 8)

  top_xgb <- importance_xgb %>% dplyr::rename(importance = Overall) %>% dplyr::slice_max(order_by = importance, n = 12)
  shap_proxy <- model_df %>%
    dplyr::select(all_of(top_xgb$feature)) %>%
    scale() %>%
    as.data.frame() %>%
    stats::setNames(top_xgb$feature)
  shap_proxy <- sweep(shap_proxy, 2, top_xgb$importance / sum(top_xgb$importance), `*`)
  shap_long <- reshape2::melt(as.matrix(shap_proxy), varnames = c("patient_index", "feature"), value.name = "shap_proxy")
  shap_long <- shap_long %>% dplyr::group_by(feature) %>% dplyr::mutate(abs_rank = rank(-abs(shap_proxy))) %>% dplyr::ungroup()
  write.csv(shap_long, file.path(res_dir, "shap_proxy_values.csv"), row.names = FALSE)
  shap_plot <- shap_long %>%
    dplyr::filter(abs_rank <= 120) %>%
    ggplot2::ggplot(ggplot2::aes(x = feature, y = shap_proxy, color = shap_proxy)) +
    ggplot2::geom_jitter(width = 0.18, alpha = 0.45, size = 1.1) +
    ggplot2::coord_flip() +
    ggplot2::scale_color_gradient2(low = "navy", mid = "grey80", high = "firebrick3") +
    ggplot2::labs(title = "SHAP-like contribution profile", x = "Feature", y = "Contribution proxy", color = "Proxy") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(fig_dir, "prediction_shap_proxy.pdf"), shap_plot, width = 10, height = 7)

  calibration_df <- pred_df %>%
    dplyr::group_by(model) %>%
    dplyr::mutate(bin = dplyr::ntile(Responder, 10)) %>%
    dplyr::group_by(model, bin) %>%
    dplyr::summarise(
      mean_pred = mean(Responder),
      observed = mean(obs == "Responder"),
      .groups = "drop"
    )
  calibration_plot <- ggplot2::ggplot(calibration_df, ggplot2::aes(x = mean_pred, y = observed, color = model)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 3, color = "grey50") +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Calibration of response probabilities", x = "Predicted response probability", y = "Observed response rate", color = "Model") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(fig_dir, "prediction_calibration.pdf"), calibration_plot, width = 8, height = 6)

  message("[05] Treatment response prediction analysis completed.")
  invisible(list(models = fits, performance = model_performance, predictions = pred_df, importance = importance_df))
}

if (sys.nframe() == 0L) {
  run_treatment_response_prediction()
}
