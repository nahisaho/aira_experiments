#!/usr/bin/env Rscript
# ============================================================================
# 02_tnd_analysis.R — Test-Negative Design (TND) の統計解析
# ============================================================================
source("R/00_setup.R")
library(dplyr); library(gnm); library(survival)
library(sandwich); library(lmtest); library(broom); library(ggplot2)

tnd_data <- readRDS("data/tnd_simulated.rds")

# --- 1. Conditional logistic regression (clogit) ---
cat("=== TND: Conditional Logistic Regression ===\n")
m_clogit <- clogit(
  tnd_case ~ vaccinated + age + sex + comorbidity +
    prior_infection + strata(calendar_week),
  data = tnd_data)
print(summary(m_clogit))

or_vax <- exp(coef(m_clogit)["vaccinated"])
ci_vax <- exp(confint(m_clogit)["vaccinated", ])
ve_overall <- 1 - or_vax
ve_ci <- 1 - rev(ci_vax)
cat(sprintf("\nOverall VE: %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
  ve_overall * 100, ve_ci[1] * 100, ve_ci[2] * 100))

# --- 2. gnm for large strata ---
cat("\n=== TND: gnm with eliminable strata ===\n")
m_gnm <- gnm(
  tnd_case ~ vaccinated + age + sex + comorbidity + prior_infection,
  family = binomial(link = "logit"),
  eliminate = factor(calendar_week),
  data = tnd_data)
print(summary(m_gnm))

# --- 3. Falsification test ---
cat("\n=== Falsification Test ===\n")
df_falsi <- tnd_data %>% filter(target_positive == 0)
m_falsi <- glm(
  nontarget_positive ~ vaccinated + age + sex + comorbidity +
    prior_infection + factor(calendar_week),
  family = binomial, data = df_falsi)
or_falsi <- exp(coef(m_falsi)["vaccinated"])
ci_falsi <- exp(confint.default(m_falsi)["vaccinated", ])
ve_falsi <- 1 - or_falsi
cat(sprintf("Falsification VE (non-target): %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
  ve_falsi * 100, (1 - ci_falsi[2]) * 100, (1 - ci_falsi[1]) * 100))

# --- 4. Robust SE ---
robust_se <- coeftest(m_gnm, vcov = vcovHC(m_gnm, type = "HC1"))
print(robust_se)

# --- Save ---
tnd_results <- data.frame(
  model = c("clogit", "gnm"), method = c("survival::clogit", "gnm::gnm"),
  ve_percent = c(ve_overall * 100, (1 - exp(coef(m_gnm)["vaccinated"])) * 100),
  ve_lower = c(ve_ci[1] * 100, NA), ve_upper = c(ve_ci[2] * 100, NA))
write.csv(tnd_results, "results/tnd_ve_estimates.csv", row.names = FALSE)
write.csv(data.frame(test = "Non-target VE", ve = ve_falsi * 100,
  ci_lower = (1 - ci_falsi[2]) * 100, ci_upper = (1 - ci_falsi[1]) * 100,
  pass = abs(ve_falsi) < 0.10), "results/falsification_test.csv", row.names = FALSE)
cat("\n✓ TND analysis complete.\n")
