#!/usr/bin/env Rscript
# =============================================================================
# 03_cohort_confounding.R
# 長期コホート研究の交絡調整戦略
# =============================================================================

library(tidyverse)
library(survival)
library(mgcv)
library(lme4)

# --- 1. Cox比例ハザードモデル (段階的交絡調整) ---
run_cohort_cox <- function(cohort_data) {
  cox_m1 <- coxph(Surv(follow_up_years, death) ~ pm25_mean + age + sex,
                  data = cohort_data)
  cox_m2 <- coxph(Surv(follow_up_years, death) ~ pm25_mean + age + sex +
    bmi + smoking_status + alcohol + education + income +
    physical_activity + diet_score, data = cohort_data)
  cox_m3 <- coxph(Surv(follow_up_years, death) ~ pm25_mean + age + sex +
    bmi + smoking_status + alcohol + education + income +
    physical_activity + diet_score +
    area_deprivation + urbanicity + greenspace_pct, data = cohort_data)
  cox_m4 <- coxph(Surv(follow_up_years, death) ~ pm25_mean + o3_mean +
    age + sex + bmi + smoking_status + alcohol + education + income +
    physical_activity + diet_score +
    area_deprivation + urbanicity + greenspace_pct, data = cohort_data)
  
  extract_hr <- function(model, var = "pm25_mean") {
    hr <- exp(coef(model)[var] * 10)
    ci <- exp(confint(model)[var, ] * 10)
    c(HR = hr, CI_low = ci[1], CI_high = ci[2])
  }
  
  tibble(
    model = paste0("Model ", 1:4),
    description = c("Age + Sex", "+ Individual confounders",
                    "+ Area-level confounders", "+ O3 co-adjustment"),
    HR_per10 = sapply(list(cox_m1, cox_m2, cox_m3, cox_m4),
                      function(m) extract_hr(m)[1]),
    CI_low = sapply(list(cox_m1, cox_m2, cox_m3, cox_m4),
                    function(m) extract_hr(m)[2]),
    CI_high = sapply(list(cox_m1, cox_m2, cox_m3, cox_m4),
                     function(m) extract_hr(m)[3])
  )
}

# --- 2. 傾向スコア重み付け (IPW) ---
run_ipw_analysis <- function(cohort_data) {
  cohort_data <- cohort_data %>%
    mutate(pm25_high = ifelse(pm25_mean > median(pm25_mean), 1, 0))
  
  ps_model <- glm(pm25_high ~ age + sex + bmi + smoking_status + alcohol +
    education + income + physical_activity + diet_score +
    area_deprivation + urbanicity, family = binomial(), data = cohort_data)
  
  cohort_data$ps <- predict(ps_model, type = "response")
  cohort_data <- cohort_data %>%
    mutate(
      iptw_stab = ifelse(pm25_high == 1, mean(pm25_high)/ps,
                         (1-mean(pm25_high))/(1-ps)),
      iptw_trim = pmin(pmax(iptw_stab, quantile(iptw_stab, 0.01)),
                       quantile(iptw_stab, 0.99))
    )
  
  cox_ipw <- coxph(Surv(follow_up_years, death) ~ pm25_high,
                   weights = iptw_trim, data = cohort_data)
  list(ps_model = ps_model, cox_ipw = cox_ipw)
}

cat("03_cohort_confounding.R loaded successfully\n")
