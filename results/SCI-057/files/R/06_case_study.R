#!/usr/bin/env Rscript
# =============================================================================
# 06_case_study.R
# PM2.5/O3 全死亡・心血管疾患リスク評価 ケーススタディ
# 実行: Rscript R/06_case_study.R
# =============================================================================

library(tidyverse)
library(dlnm)
library(mgcv)
library(survival)
library(EValue)
library(splines)

source("R/02_timeseries_design.R")
source("R/03_cohort_confounding.R")
source("R/04_nonlinear_erf.R")
source("R/05_sensitivity_evalue.R")

# =============================================================================
# シミュレーションデータ生成
# =============================================================================
set.seed(42)
n_days <- 365 * 10
dates <- seq(as.Date("2013-01-01"), by = "day", length.out = n_days)

ts_data <- tibble(
  date = dates, date_num = as.numeric(dates - min(dates)),
  dow = wday(dates), year = year(dates), month = month(dates),
  temp = 15 + 10*sin(2*pi*date_num/365) + rnorm(n_days, 0, 3),
  rh = 60 + 15*sin(2*pi*(date_num+90)/365) + rnorm(n_days, 0, 8),
  pm25 = pmax(5, 25 + 15*sin(2*pi*(date_num+180)/365) + 0.3*temp + rnorm(n_days,0,10)),
  o3_8h = pmax(10, 40 + 20*sin(2*pi*date_num/365) + 0.5*temp - 0.2*pm25 + rnorm(n_days,0,12)),
  log_mu = log(50) + 0.0008*pm25 + 0.0005*o3_8h + 0.001*(temp-20)^2/100 -
    0.0002*date_num/365 + 0.02*sin(2*pi*date_num/365),
  deaths = rpois(n_days, exp(log_mu))
) %>% mutate(
  pm25_lag01 = (pm25 + lag(pm25))/2, o3_8h_lag01 = (o3_8h + lag(o3_8h))/2,
  temp_lag01 = (temp + lag(temp))/2, rh_lag01 = (rh + lag(rh))/2
) %>% filter(!is.na(pm25_lag01))

n_years <- length(unique(ts_data$year))

# コホートデータ
set.seed(123)
n_subjects <- 50000
cohort_data <- tibble(
  subject_id = 1:n_subjects, area_id = sample(1:100, n_subjects, replace = TRUE),
  age = rnorm(n_subjects, 55, 12), sex = rbinom(n_subjects, 1, 0.48),
  bmi = rnorm(n_subjects, 25, 4),
  smoking_status = sample(0:2, n_subjects, replace=TRUE, prob=c(0.5,0.2,0.3)),
  alcohol = rbinom(n_subjects, 1, 0.4), education = sample(1:4, n_subjects, replace=TRUE),
  income = rlnorm(n_subjects, log(40000), 0.5),
  physical_activity = rbinom(n_subjects, 1, 0.35),
  diet_score = rnorm(n_subjects, 50, 15),
  area_deprivation = rnorm(n_subjects, 0, 1),
  urbanicity = rbinom(n_subjects, 1, 0.7),
  greenspace_pct = rbeta(n_subjects, 3, 7)*100,
  pm25_mean = pmax(5, 12 + 0.5*area_deprivation - 0.1*greenspace_pct + rnorm(n_subjects,0,4)),
  o3_mean = pmax(20, 35 - 0.3*area_deprivation + rnorm(n_subjects,0,8))
) %>% mutate(
  log_hazard = log(0.008) + 0.006*pm25_mean + 0.003*o3_mean +
    0.03*(age-55)/10 - 0.15*sex + 0.02*(bmi-25) +
    0.3*(smoking_status==2) + 0.1*(smoking_status==1),
  hazard = exp(log_hazard),
  follow_up_years = pmin(rexp(n_subjects, hazard), 15),
  death = ifelse(follow_up_years < 15, 1, 0)
)

# =============================================================================
# 解析実行
# =============================================================================
cat("=== Running DLNM ===\n")
dlnm_res <- run_dlnm_analysis(ts_data, max_lag = 21)

cat("=== Running GAM ERF ===\n")
gam_res <- fit_gam_erf(ts_data)

cat("=== Running Cox Models ===\n")
cox_res <- run_cohort_cox(cohort_data)
print(cox_res)

cat("=== Computing E-values ===\n")
evals <- compute_evalues(cox_res %>% rename(HR = HR_per10))
print(evals)

cat("=== Two-pollutant sensitivity ===\n")
two_poll <- two_pollutant_sensitivity(ts_data)
print(two_poll)

# 結果保存
write_csv(cox_res, "results/cox_model_results.csv")
write_csv(two_poll, "results/two_pollutant_results.csv")

cat("\n=== Case Study Complete ===\n")
