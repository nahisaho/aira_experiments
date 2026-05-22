#!/usr/bin/env Rscript
# =============================================================================
# 05_sensitivity_evalue.R
# 感度分析: 未測定交絡のE-value計算 & その他感度分析
# =============================================================================

library(tidyverse)
library(EValue)

# --- 1. E-value計算 ---
compute_evalues <- function(results_table) {
  results_table %>%
    rowwise() %>%
    mutate(
      evalue_point = evalues.HR(HR, rare = FALSE)["E-value", "point"],
      evalue_ci    = evalues.HR(HR, rare = FALSE)["E-value",
                       ifelse(HR > 1, "lower", "upper")]
    ) %>%
    ungroup()
}

# --- 2. Quantitative Bias Analysis (QBA) ---
run_qba <- function(hr_crude, n_sims = 10000,
                    alpha_e = 4, beta_e = 6,
                    alpha_u = 2, beta_u = 8,
                    log_rr_mean = log(1.5), log_rr_sd = 0.3) {
  p_conf_exposed   <- rbeta(n_sims, alpha_e, beta_e)
  p_conf_unexposed <- rbeta(n_sims, alpha_u, beta_u)
  rr_confounder    <- rlnorm(n_sims, log_rr_mean, log_rr_sd)
  
  bias_factor <- (p_conf_exposed * (rr_confounder - 1) + 1) /
                 (p_conf_unexposed * (rr_confounder - 1) + 1)
  hr_adjusted <- hr_crude / bias_factor
  
  list(
    hr_adjusted_median = median(hr_adjusted),
    hr_adjusted_ci = quantile(hr_adjusted, c(0.025, 0.975)),
    bias_factor_median = median(bias_factor)
  )
}

# --- 3. 二汚染物質モデル ---
two_pollutant_sensitivity <- function(ts_data) {
  m1 <- glm(deaths ~ pm25 + ns(temp, df = 5) +
    ns(date_num, df = 7 * n_years) + factor(dow),
    family = quasipoisson(), data = ts_data)
  m2 <- glm(deaths ~ o3_8h + ns(temp, df = 5) +
    ns(date_num, df = 7 * n_years) + factor(dow),
    family = quasipoisson(), data = ts_data)
  m3 <- glm(deaths ~ pm25 + o3_8h + ns(temp, df = 5) +
    ns(date_num, df = 7 * n_years) + factor(dow),
    family = quasipoisson(), data = ts_data)
  
  extract_rr <- function(model, var) {
    rr <- exp(coef(model)[var] * 10)
    ci <- exp(confint.default(model)[var, ] * 10)
    c(rr = rr, ci_low = ci[1], ci_high = ci[2])
  }
  
  tibble(
    model = c("PM2.5 only", "O3 only", "PM2.5 (adj O3)", "O3 (adj PM2.5)"),
    pollutant = c("PM2.5", "O3", "PM2.5", "O3"),
    rr_per10 = c(extract_rr(m1,"pm25")[1], extract_rr(m2,"o3_8h")[1],
                 extract_rr(m3,"pm25")[1], extract_rr(m3,"o3_8h")[1]),
    ci_low = c(extract_rr(m1,"pm25")[2], extract_rr(m2,"o3_8h")[2],
               extract_rr(m3,"pm25")[2], extract_rr(m3,"o3_8h")[2]),
    ci_high = c(extract_rr(m1,"pm25")[3], extract_rr(m2,"o3_8h")[3],
                extract_rr(m3,"pm25")[3], extract_rr(m3,"o3_8h")[3])
  )
}

# --- 4. Negative Control Outcomes ---
run_negative_controls <- function(ts_data, exposure_var = "pm25") {
  negative_outcomes <- c("fractures", "appendicitis", "dental_caries", "hernia")
  map_dfr(negative_outcomes, function(outcome) {
    formula <- as.formula(paste(outcome, "~", exposure_var,
      "+ ns(temp, df = 5) + ns(date_num, df = 7 * n_years) + factor(dow)"))
    model <- glm(formula, family = quasipoisson(), data = ts_data)
    rr <- exp(coef(model)[exposure_var] * 10)
    ci <- exp(confint.default(model)[exposure_var, ] * 10)
    tibble(outcome = outcome, rr_per10 = rr, ci_low = ci[1], ci_high = ci[2])
  })
}

cat("05_sensitivity_evalue.R loaded successfully\n")
