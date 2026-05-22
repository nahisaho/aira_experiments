#!/usr/bin/env Rscript
# ============================================================================
# 01_simulate_data.R — リアルワールドデータのシミュレーション
# TND研究デザインに対応した合成データ生成
# ============================================================================
source("R/00_setup.R")
library(dplyr)
library(lubridate)

set.seed(20260523)
N <- 50000

simulate_tnd_data <- function(N) {
  df <- data.frame(
    id                = 1:N,
    age               = pmin(pmax(round(rnorm(N, 55, 18)), 18), 95),
    sex               = rbinom(N, 1, 0.48),
    comorbidity       = rbinom(N, 1, 0.35),
    healthcare_worker = rbinom(N, 1, 0.08),
    prior_infection   = rbinom(N, 1, 0.15),
    calendar_week     = sample(1:52, N, replace = TRUE)
  )

  # Healthy vaccinee bias: health-seeking → higher vax uptake
  df$health_score <- with(df,
    -0.3 * (age > 65) + 0.4 * healthcare_worker -
    0.2 * comorbidity + rnorm(N, 0, 0.5))
  df$p_vax <- plogis(0.2 + 0.5 * df$health_score +
    0.3 * (df$age > 50) + 0.2 * df$healthcare_worker)
  df$vaccinated <- rbinom(N, 1, df$p_vax)

  df$days_since_vax <- ifelse(df$vaccinated == 1,
    pmax(7, round(rexp(N, 1/120))), NA)
  df$dose <- ifelse(df$vaccinated == 0, 0,
    ifelse(df$days_since_vax < 21, 1, ifelse(runif(N) < 0.6, 3, 2)))

  # Variant by calendar period
  df$variant <- case_when(
    df$calendar_week <= 13 ~ "Alpha",
    df$calendar_week <= 26 ~ "Delta",
    df$calendar_week <= 39 ~ "Omicron_BA1",
    TRUE                   ~ "Omicron_BA5")

  # True VE with waning
  base_ve <- case_when(
    df$variant == "Alpha"       ~ 0.88,
    df$variant == "Delta"       ~ 0.82,
    df$variant == "Omicron_BA1" ~ 0.55,
    df$variant == "Omicron_BA5" ~ 0.45)
  waning <- ifelse(is.na(df$days_since_vax), 1,
    exp(-0.005 * pmax(0, df$days_since_vax - 14)))
  booster <- ifelse(df$dose == 3, 0.15, 0)
  true_ve <- pmin(0.95, base_ve * waning + booster)
  true_ve <- ifelse(df$vaccinated == 0, 0, true_ve)
  df$true_ve <- true_ve

  # TND outcomes
  df$p_target <- 0.25 * (1 - true_ve) *
    (1 + 0.3 * df$comorbidity + 0.002 * (df$age - 50))
  df$p_target <- pmin(pmax(df$p_target, 0.01), 0.60)
  df$p_nontarget <- 0.15 * (1 + 0.1 * df$comorbidity)
  df$target_positive <- rbinom(N, 1, df$p_target)
  df$nontarget_positive <- ifelse(df$target_positive == 0,
    rbinom(N, 1, df$p_nontarget / (1 - df$p_target)), 0)
  df$tnd_case <- df$target_positive

  # Hospitalization
  hosp_ve <- pmin(0.98, true_ve * 1.15 + 0.05)
  hosp_ve <- ifelse(df$vaccinated == 0, 0, hosp_ve)
  df$p_hospitalization <- ifelse(df$tnd_case == 1,
    0.08 * (1 - hosp_ve) * (1 + 0.5 * df$comorbidity + 0.01 * pmax(0, df$age - 60)),
    0.02 * (1 + 0.3 * df$comorbidity))
  df$p_hospitalization <- pmin(df$p_hospitalization, 0.50)
  df$hospitalized <- rbinom(N, 1, df$p_hospitalization)

  # Time categories
  df$time_since_vax_cat <- cut(df$days_since_vax,
    breaks = c(0, 13, 59, 119, 179, Inf),
    labels = c("0-13d", "14-59d", "60-119d", "120-179d", "180d+"), right = TRUE)
  df$vax_status <- ifelse(df$vaccinated == 0, "Unvaccinated",
    paste0("Dose", df$dose, "_", as.character(df$time_since_vax_cat)))

  return(df)
}

tnd_data <- simulate_tnd_data(N)
saveRDS(tnd_data, "data/tnd_simulated.rds")
write.csv(tnd_data, "data/tnd_simulated.csv", row.names = FALSE)

cat(sprintf("✓ Simulated TND data: N=%d, Cases=%d, Controls=%d\n",
  nrow(tnd_data), sum(tnd_data$tnd_case), sum(tnd_data$tnd_case == 0)))
cat(sprintf("  Vaccinated: %d (%.1f%%)\n",
  sum(tnd_data$vaccinated), 100 * mean(tnd_data$vaccinated)))
