#!/usr/bin/env Rscript
# =============================================================================
# 04_nonlinear_erf.R
# 暴露反応関数の非線形モデリング (GAM / スプライン)
# =============================================================================

library(tidyverse)
library(mgcv)
library(splines)

# --- 1. GAMによる暴露-反応関数推定 ---
fit_gam_erf <- function(ts_data) {
  gam_pm25 <- gam(
    deaths ~ s(pm25, bs = "cr", k = 10) + s(temp, bs = "cr", k = 10) +
      s(rh, bs = "cr", k = 5) + s(date_num, bs = "cr", k = 100) + factor(dow),
    family = quasipoisson(), data = ts_data, method = "REML"
  )
  gam_o3 <- gam(
    deaths ~ s(o3_8h, bs = "cr", k = 10) + s(temp, bs = "cr", k = 10) +
      s(rh, bs = "cr", k = 5) + s(date_num, bs = "cr", k = 100) + factor(dow),
    family = quasipoisson(), data = ts_data, method = "REML"
  )
  list(gam_pm25 = gam_pm25, gam_o3 = gam_o3)
}

# --- 2. 線形 vs 非線形モデル比較 ---
compare_linearity <- function(ts_data) {
  m_linear <- glm(deaths ~ pm25 + ns(temp, df = 5) +
    ns(date_num, df = 7 * n_years) + factor(dow),
    family = quasipoisson(), data = ts_data)
  m_ns <- glm(deaths ~ ns(pm25, df = 4) + ns(temp, df = 5) +
    ns(date_num, df = 7 * n_years) + factor(dow),
    family = quasipoisson(), data = ts_data)
  m_gam <- gam(deaths ~ s(pm25, bs = "cr", k = 10) +
    s(temp, bs = "cr", k = 10) + s(date_num, bs = "cr", k = 100) + factor(dow),
    family = quasipoisson(), data = ts_data, method = "REML")
  
  tibble(
    model = c("Linear", "Natural Spline (4df)", "GAM (REML)"),
    aic = c(AIC(m_linear), AIC(m_ns), AIC(m_gam)),
    deviance = c(deviance(m_linear), deviance(m_ns), deviance(m_gam))
  )
}

# --- 3. Shape-Constrained Additive Model (SCAM) ---
fit_scam_erf <- function(ts_data) {
  library(scam)
  scam_pm25 <- scam(
    deaths ~ s(pm25, bs = "mpi", k = 10) + s(temp, bs = "cr", k = 10) +
      s(date_num, bs = "cr", k = 100) + factor(dow),
    family = quasipoisson(), data = ts_data
  )
  scam_pm25
}

cat("04_nonlinear_erf.R loaded successfully\n")
