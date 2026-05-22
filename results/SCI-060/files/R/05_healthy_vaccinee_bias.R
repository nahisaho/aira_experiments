#!/usr/bin/env Rscript
# ============================================================================
# 05_healthy_vaccinee_bias.R — 健康バイアス (Healthy Vaccinee Bias) の補正
# ============================================================================
source("R/00_setup.R")
library(dplyr); library(MatchIt); library(WeightIt); library(cobalt)
library(sandwich); library(lmtest); library(ggplot2); library(survival)

tnd_data <- readRDS("data/tnd_simulated.rds")

# ============================================================================
# 1. バイアスの検出: 共変量分布の不均衡
# ============================================================================
cat("=== Detecting Healthy Vaccinee Bias ===\n")

# ワクチン接種群 vs 非接種群の共変量バランス (調整前)
bal_unadj <- tnd_data %>%
  group_by(vaccinated) %>%
  summarise(
    mean_age       = mean(age),
    prop_female    = mean(sex),
    prop_comorbid  = mean(comorbidity),
    prop_hcw       = mean(healthcare_worker),
    prop_prior_inf = mean(prior_infection),
    mean_health    = mean(health_score),
    .groups = "drop")
cat("\n--- Unadjusted Balance ---\n")
print(bal_unadj)

# Standardized Mean Differences (SMD)
calc_smd <- function(x, g) {
  m1 <- mean(x[g == 1]); m0 <- mean(x[g == 0])
  s1 <- sd(x[g == 1]);   s0 <- sd(x[g == 0])
  (m1 - m0) / sqrt((s1^2 + s0^2) / 2)
}

covs <- c("age", "sex", "comorbidity", "healthcare_worker", "prior_infection")
smd_unadj <- sapply(covs, function(v) calc_smd(tnd_data[[v]], tnd_data$vaccinated))
cat("\nSMD (unadjusted):\n")
print(round(smd_unadj, 3))

# ============================================================================
# 2. 傾向スコア推定
# ============================================================================
cat("\n=== Propensity Score Estimation ===\n")

ps_model <- glm(
  vaccinated ~ age + I(age^2) + sex + comorbidity +
    healthcare_worker + prior_infection + factor(calendar_week),
  family = binomial, data = tnd_data)

tnd_data$ps <- predict(ps_model, type = "response")

# PS distribution
p_ps <- ggplot(tnd_data, aes(x = ps, fill = factor(vaccinated))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("0" = "#CB181D", "1" = "#2171B5"),
    labels = c("Unvaccinated", "Vaccinated")) +
  labs(title = "Propensity Score Distribution",
    x = "Propensity Score", y = "Density", fill = "") +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")
ggsave("figures/propensity_score_dist.png", p_ps, width = 8, height = 5, dpi = 300)

# ============================================================================
# 3. Method A: Inverse Probability of Treatment Weighting (IPTW)
# ============================================================================
cat("\n=== IPTW Adjustment ===\n")

W <- weightit(
  vaccinated ~ age + I(age^2) + sex + comorbidity +
    healthcare_worker + prior_infection,
  data = tnd_data, method = "ps", estimand = "ATE")

tnd_data$iptw <- W$weights

# Stabilized weights
tnd_data$iptw_stab <- ifelse(tnd_data$vaccinated == 1,
  mean(tnd_data$vaccinated) / tnd_data$ps,
  (1 - mean(tnd_data$vaccinated)) / (1 - tnd_data$ps))

# Truncation at 1st/99th percentiles
q_low <- quantile(tnd_data$iptw_stab, 0.01)
q_high <- quantile(tnd_data$iptw_stab, 0.99)
tnd_data$iptw_trunc <- pmin(pmax(tnd_data$iptw_stab, q_low), q_high)

cat(sprintf("Weight range (stabilized, truncated): [%.2f, %.2f]\n",
  min(tnd_data$iptw_trunc), max(tnd_data$iptw_trunc)))

# IPTW-weighted TND model
m_iptw <- glm(
  tnd_case ~ vaccinated + age + sex + comorbidity + prior_infection,
  family = binomial, data = tnd_data, weights = iptw_trunc)

# Robust SE with sandwich estimator
coef_iptw <- coeftest(m_iptw, vcov = vcovHC(m_iptw, type = "HC0"))
or_iptw <- exp(coef_iptw["vaccinated", "Estimate"])
ci_iptw <- exp(coef_iptw["vaccinated", "Estimate"] +
  c(-1, 1) * 1.96 * coef_iptw["vaccinated", "Std. Error"])
ve_iptw <- (1 - or_iptw) * 100
cat(sprintf("IPTW VE: %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
  ve_iptw, (1 - ci_iptw[2]) * 100, (1 - ci_iptw[1]) * 100))

# ============================================================================
# 4. Method B: Propensity Score Matching
# ============================================================================
cat("\n=== PS Matching ===\n")

m_match <- matchit(
  vaccinated ~ age + sex + comorbidity + healthcare_worker + prior_infection,
  data = tnd_data, method = "nearest", ratio = 1, caliper = 0.1)

matched_data <- match.data(m_match)
cat(sprintf("Matched sample: N = %d\n", nrow(matched_data)))

# Balance after matching
bal_matched <- bal.tab(m_match, stats = c("m", "v"), thresholds = c(m = 0.1))
print(bal_matched)

# VE in matched sample
m_matched <- clogit(
  tnd_case ~ vaccinated + strata(subclass),
  data = matched_data)
or_match <- exp(coef(m_matched)["vaccinated"])
ci_match <- exp(confint(m_matched)["vaccinated", ])
ve_match <- (1 - or_match) * 100
cat(sprintf("Matched VE: %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
  ve_match, (1 - ci_match[2]) * 100, (1 - ci_match[1]) * 100))

# ============================================================================
# 5. Method C: Negative Control Outcome Adjustment
# ============================================================================
cat("\n=== Negative Control Outcome ===\n")

# 非COVID入院をnegative control outcomeとして使用
# ワクチンが真に影響しないアウトカムでの「見かけの効果」をバイアスの指標とする
df_controls <- tnd_data %>% filter(tnd_case == 0)
m_nco <- glm(
  nontarget_positive ~ vaccinated + age + sex + comorbidity +
    prior_infection + factor(calendar_week),
  family = binomial, data = df_controls)
or_nco <- exp(coef(m_nco)["vaccinated"])
bias_factor <- or_nco
cat(sprintf("Negative control OR: %.3f (bias factor)\n", bias_factor))

# Bias-adjusted VE
m_naive <- glm(
  tnd_case ~ vaccinated + age + sex + comorbidity +
    prior_infection + factor(calendar_week),
  family = binomial, data = tnd_data)
or_naive <- exp(coef(m_naive)["vaccinated"])
or_adjusted <- or_naive / bias_factor
ve_adjusted <- (1 - or_adjusted) * 100
cat(sprintf("Bias-adjusted VE: %.1f%%\n", ve_adjusted))

# ============================================================================
# 6. E-value for unmeasured confounding
# ============================================================================
cat("\n=== E-value Analysis ===\n")

# E-value: 観察された効果をnullにするために必要な未測定交絡の強さ
# VE = 1 - OR → OR = 1 - VE/100
or_for_eval <- or_naive
if (or_for_eval < 1) {
  rr_approx <- or_for_eval  # approximate RR ≈ OR when outcome is rare
  eval_point <- rr_approx + sqrt(rr_approx * (rr_approx - 1))
  # For protective exposure, transform
  eval_point <- 1 / rr_approx + sqrt(1 / rr_approx * (1 / rr_approx - 1))
  cat(sprintf("E-value (point estimate): %.2f\n", eval_point))
  cat("Interpretation: unmeasured confounder would need to be associated with\n")
  cat(sprintf("both vaccination and outcome by a factor of ≥%.2f to explain away VE\n", eval_point))
}

# ============================================================================
# 7. 結果の比較とバランスプロット
# ============================================================================

comparison <- data.frame(
  method = c("Naive (no bias correction)", "IPTW", "PS Matching",
    "Negative Control Adjusted"),
  ve = c((1 - or_naive) * 100, ve_iptw, ve_match, ve_adjusted))

cat("\n--- VE Comparison ---\n")
print(comparison)

p_compare <- ggplot(comparison, aes(x = reorder(method, ve), y = ve)) +
  geom_col(fill = "#2171B5", alpha = 0.7, width = 0.6) +
  geom_hline(yintercept = 0, color = "gray50") +
  coord_flip() +
  labs(title = "VE Estimates: Impact of Healthy Vaccinee Bias Correction",
    x = "", y = "Vaccine Effectiveness (%)") +
  theme_minimal(base_size = 13)
ggsave("figures/bias_correction_comparison.png", p_compare, width = 10, height = 5, dpi = 300)

# Love plot for covariate balance
p_love <- love.plot(m_match, stats = "mean.diffs", abs = TRUE,
  thresholds = c(m = 0.1), colors = c("#CB181D", "#2171B5"),
  sample.names = c("Before Matching", "After Matching"),
  title = "Covariate Balance: Before vs After PS Matching")
ggsave("figures/love_plot_balance.png", p_love, width = 9, height = 6, dpi = 300)

write.csv(comparison, "results/bias_correction_comparison.csv", row.names = FALSE)
cat("\n✓ Healthy vaccinee bias analysis complete.\n")
