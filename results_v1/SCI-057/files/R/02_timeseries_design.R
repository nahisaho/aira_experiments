#!/usr/bin/env Rscript
# =============================================================================
# 02_timeseries_design.R
# 時系列研究デザイン: ケースクロスオーバー & DLNM
# =============================================================================

library(tidyverse)
library(dlnm)
library(gnm)
library(splines)
library(survival)

# --- 1. ケースクロスオーバー研究 ---
run_case_crossover <- function(mortality_data, exposure_data, confounders) {
  analysis_data <- mortality_data %>%
    mutate(
      year  = year(date),
      month = month(date),
      dow   = wday(date),
      stratum = paste(year, month, dow, sep = "-")
    ) %>%
    left_join(exposure_data, by = "date") %>%
    left_join(confounders, by = "date")
  
  model_pm25 <- clogit(
    death ~ pm25_lag01 + temp_lag01 + rh_lag01 + strata(stratum),
    data = analysis_data
  )
  
  model_o3 <- clogit(
    death ~ o3_8h_lag01 + temp_lag01 + rh_lag01 + strata(stratum),
    data = analysis_data
  )
  
  or_pm25 <- exp(coef(model_pm25)["pm25_lag01"] * 10)
  ci_pm25 <- exp(confint(model_pm25)["pm25_lag01", ] * 10)
  
  list(model_pm25 = model_pm25, model_o3 = model_o3,
       or_pm25_per10 = or_pm25, ci_pm25_per10 = ci_pm25)
}

# --- 2. DLNM (Distributed Lag Non-linear Model) ---
run_dlnm_analysis <- function(ts_data, max_lag = 21) {
  cb_pm25 <- crossbasis(
    ts_data$pm25, lag = max_lag,
    argvar = list(fun = "ns", df = 4),
    arglag = list(fun = "ns", df = 3)
  )
  
  cb_o3 <- crossbasis(
    ts_data$o3_8h, lag = max_lag,
    argvar = list(fun = "ns", df = 4),
    arglag = list(fun = "ns", df = 3)
  )
  
  cb_temp <- crossbasis(
    ts_data$temp, lag = max_lag,
    argvar = list(fun = "ns", df = 5),
    arglag = list(fun = "ns", df = 4)
  )
  
  model_pm25 <- glm(
    deaths ~ cb_pm25 + cb_temp + ns(rh, df = 3) +
      ns(date_num, df = 7 * n_years) + factor(dow),
    family = quasipoisson(), data = ts_data
  )
  
  model_o3 <- glm(
    deaths ~ cb_o3 + cb_temp + ns(rh, df = 3) +
      ns(date_num, df = 7 * n_years) + factor(dow),
    family = quasipoisson(), data = ts_data
  )
  
  pred_pm25 <- crosspred(cb_pm25, model_pm25,
    at = seq(0, 100, by = 1), cen = median(ts_data$pm25))
  pred_o3 <- crosspred(cb_o3, model_o3,
    at = seq(0, 120, by = 1), cen = median(ts_data$o3_8h))
  
  list(model_pm25 = model_pm25, model_o3 = model_o3,
       pred_pm25 = pred_pm25, pred_o3 = pred_o3)
}

# --- 3. DLNM感度分析 ---
dlnm_sensitivity <- function(ts_data, max_lags = c(14, 21, 28),
                              df_vars = c(3, 4, 5), df_lags = c(3, 4, 5)) {
  results <- expand.grid(max_lag = max_lags, df_var = df_vars, df_lag = df_lags)
  results$rr_per10 <- NA; results$ci_low <- NA; results$ci_high <- NA
  
  for (i in 1:nrow(results)) {
    cb <- crossbasis(ts_data$pm25, lag = results$max_lag[i],
      argvar = list(fun = "ns", df = results$df_var[i]),
      arglag = list(fun = "ns", df = results$df_lag[i]))
    model <- glm(deaths ~ cb + ns(date_num, df = 7 * n_years) + factor(dow),
      family = quasipoisson(), data = ts_data)
    pred <- crosspred(cb, model,
      at = c(median(ts_data$pm25), median(ts_data$pm25) + 10),
      cen = median(ts_data$pm25), cumul = TRUE)
    results$rr_per10[i] <- pred$allRRfit["10"]
    results$ci_low[i]   <- pred$allRRlow["10"]
    results$ci_high[i]  <- pred$allRRhigh["10"]
  }
  results
}

cat("02_timeseries_design.R loaded successfully\n")
