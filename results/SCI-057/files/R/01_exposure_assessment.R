#!/usr/bin/env Rscript
# =============================================================================
# 01_exposure_assessment.R
# 暴露評価モデル: LUR, 化学輸送モデル(CTM)統合, 衛星データ融合
# =============================================================================

library(tidyverse)
library(sf)
library(mgcv)

# --- 1. Land Use Regression (LUR) モデル ---
build_lur_model <- function(monitor_data, predictors) {
  formula_str <- paste("pm25_annual ~", paste(predictors, collapse = " + "))
  
  lur_full <- lm(as.formula(formula_str), data = monitor_data)
  lur_step <- step(lur_full, direction = "both", k = log(nrow(monitor_data)))
  
  # Leave-one-out cross-validation
  cv_results <- sapply(1:nrow(monitor_data), function(i) {
    fit <- lm(formula(lur_step), data = monitor_data[-i, ])
    predict(fit, newdata = monitor_data[i, ])
  })
  
  loocv_r2 <- cor(cv_results, monitor_data$pm25_annual)^2
  loocv_rmse <- sqrt(mean((cv_results - monitor_data$pm25_annual)^2))
  
  list(model = lur_step, loocv_r2 = loocv_r2, loocv_rmse = loocv_rmse)
}

# --- 2. 化学輸送モデル(CTM)校正 ---
calibrate_ctm <- function(ctm_grid, monitor_obs) {
  merged <- inner_join(monitor_obs, ctm_grid, by = c("lon", "lat", "date"))
  bias_model <- lm(pm25_obs ~ pm25_ctm + temp + rh + pbl_height, data = merged)
  merged$residual <- residuals(bias_model)
  list(bias_model = bias_model, r2 = summary(bias_model)$r.squared)
}

# --- 3. 衛星AODデータ融合 ---
satellite_data_fusion <- function(aod_data, monitor_data, met_data) {
  merged <- aod_data %>%
    inner_join(monitor_data, by = c("date", "grid_id")) %>%
    inner_join(met_data, by = c("date", "grid_id"))
  
  library(lme4)
  fusion_model <- lmer(
    pm25 ~ aod + temp + rh + pbl_height + wind_speed +
      (1 + aod | date) + (1 | grid_id),
    data = merged
  )
  
  list(model = fusion_model, r2 = MuMIn::r.squaredGLMM(fusion_model))
}

# --- 4. 暴露推定パイプライン統合 ---
estimate_exposure <- function(subject_addresses, exposure_surfaces) {
  exposures <- subject_addresses %>%
    rowwise() %>%
    mutate(
      pm25_mean = extract_exposure(lon, lat, date_start, date_end,
                                    exposure_surfaces, "pm25"),
      o3_mean   = extract_exposure(lon, lat, date_start, date_end,
                                    exposure_surfaces, "o3")
    )
  exposures
}

cat("01_exposure_assessment.R loaded successfully\n")
