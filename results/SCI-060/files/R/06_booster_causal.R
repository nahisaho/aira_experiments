#!/usr/bin/env Rscript
# ============================================================================
# 06_booster_causal.R — ブースター接種の因果推定
# ============================================================================
source("R/00_setup.R")
library(dplyr); library(survival); library(WeightIt); library(cobalt)
library(sandwich); library(lmtest); library(ggplot2)

tnd_data <- readRDS("data/tnd_simulated.rds")

# ブースター接種者 vs 2回接種者の比較 (接種者のみ)
df_vax <- tnd_data %>%
  filter(vaccinated == 1, dose %in% c(2, 3)) %>%
  mutate(booster = as.integer(dose == 3))

cat(sprintf("Analysis population: N=%d (Dose2=%d, Booster=%d)\n",
  nrow(df_vax), sum(df_vax$booster == 0), sum(df_vax$booster == 1)))

# ============================================================================
# 1. Target Trial Emulation Framework
# ============================================================================
cat("\n=== Target Trial Emulation ===\n")
cat("
Target Trial Protocol:
  Eligibility: Adults ≥18 who completed primary series (2 doses)
  Treatment:   Booster dose vs no booster
  Assignment:  Emulated via observational data + IPW
  Follow-up:   From booster eligibility date
  Outcome:     COVID-19 (test-positive in TND)
  Causal contrast: Per-protocol effect
\n")

# ============================================================================
# 2. IPTW for Booster Effect
# ============================================================================
cat("=== IPTW for Booster ===\n")

W_boost <- weightit(
  booster ~ age + I(age^2) + sex + comorbidity + healthcare_worker +
    prior_infection + days_since_vax + factor(variant),
  data = df_vax, method = "ps", estimand = "ATE")

df_vax$w_boost <- W_boost$weights

# Stabilized & truncated
df_vax$w_stab <- ifelse(df_vax$booster == 1,
  mean(df_vax$booster) / W_boost$ps,
  (1 - mean(df_vax$booster)) / (1 - W_boost$ps))
q01 <- quantile(df_vax$w_stab, 0.01)
q99 <- quantile(df_vax$w_stab, 0.99)
df_vax$w_trunc <- pmin(pmax(df_vax$w_stab, q01), q99)

# Weighted outcome model
m_boost_iptw <- glm(
  tnd_case ~ booster + age + sex + comorbidity +
    prior_infection + factor(variant),
  family = binomial, data = df_vax, weights = w_trunc)

coef_boost <- coeftest(m_boost_iptw, vcov = vcovHC(m_boost_iptw, type = "HC0"))
or_boost <- exp(coef_boost["booster", "Estimate"])
ci_boost <- exp(coef_boost["booster", "Estimate"] +
  c(-1, 1) * 1.96 * coef_boost["booster", "Std. Error"])
rve_boost <- (1 - or_boost) * 100
cat(sprintf("Relative VE (Booster vs Dose2): %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
  rve_boost, (1 - ci_boost[2]) * 100, (1 - ci_boost[1]) * 100))

# ============================================================================
# 3. Doubly Robust Estimation (AIPW)
# ============================================================================
cat("\n=== Augmented IPW (Doubly Robust) ===\n")

# Outcome model (for augmentation)
m_out <- glm(
  tnd_case ~ booster + age + sex + comorbidity +
    prior_infection + days_since_vax + factor(variant),
  family = binomial, data = df_vax)

# Predicted potential outcomes
df_vax$mu1 <- predict(m_out, newdata = transform(df_vax, booster = 1), type = "response")
df_vax$mu0 <- predict(m_out, newdata = transform(df_vax, booster = 0), type = "response")
df_vax$ps_boost <- W_boost$ps

# AIPW estimator
n <- nrow(df_vax)
aipw_1 <- with(df_vax,
  mean(mu1 + booster / ps_boost * (tnd_case - mu1)))
aipw_0 <- with(df_vax,
  mean(mu0 + (1 - booster) / (1 - ps_boost) * (tnd_case - mu0)))

ate_aipw <- aipw_1 - aipw_0
rr_aipw  <- aipw_1 / aipw_0
rve_aipw <- (1 - rr_aipw) * 100

cat(sprintf("AIPW Risk Difference: %.4f\n", ate_aipw))
cat(sprintf("AIPW Risk Ratio:      %.3f\n", rr_aipw))
cat(sprintf("AIPW Relative VE:     %.1f%%\n", rve_aipw))

# Bootstrap CI for AIPW
set.seed(42)
n_boot <- 500
boot_rve <- numeric(n_boot)
for (b in 1:n_boot) {
  idx <- sample(n, n, replace = TRUE)
  df_b <- df_vax[idx, ]
  
  ps_b <- glm(booster ~ age + sex + comorbidity + healthcare_worker +
    prior_infection + days_since_vax, family = binomial, data = df_b)$fitted
  out_b <- glm(tnd_case ~ booster + age + sex + comorbidity +
    prior_infection + days_since_vax, family = binomial, data = df_b)
  mu1_b <- predict(out_b, newdata = transform(df_b, booster = 1), type = "response")
  mu0_b <- predict(out_b, newdata = transform(df_b, booster = 0), type = "response")
  
  a1 <- mean(mu1_b + df_b$booster / ps_b * (df_b$tnd_case - mu1_b))
  a0 <- mean(mu0_b + (1 - df_b$booster) / (1 - ps_b) * (df_b$tnd_case - mu0_b))
  boot_rve[b] <- (1 - a1 / a0) * 100
}
ci_aipw <- quantile(boot_rve, c(0.025, 0.975))
cat(sprintf("AIPW 95%% CI: (%.1f%%, %.1f%%)\n", ci_aipw[1], ci_aipw[2]))

# ============================================================================
# 4. Variant-specific booster effect
# ============================================================================
cat("\n=== Variant-Specific Booster Effect ===\n")

boost_by_var <- list()
for (v in unique(df_vax$variant)) {
  df_bv <- df_vax %>% filter(variant == v)
  if (sum(df_bv$booster) < 30 || sum(1 - df_bv$booster) < 30) next
  m_bv <- glm(tnd_case ~ booster + age + sex + comorbidity + prior_infection,
    family = binomial, data = df_bv)
  or_bv <- exp(coef(m_bv)["booster"])
  ci_bv <- exp(confint.default(m_bv)["booster", ])
  boost_by_var[[v]] <- data.frame(variant = v,
    rve = (1 - or_bv) * 100,
    ci_lo = (1 - ci_bv[2]) * 100, ci_hi = (1 - ci_bv[1]) * 100)
  cat(sprintf("  %s: rVE=%.1f%% (%.1f%%-%.1f%%)\n",
    v, (1 - or_bv) * 100, (1 - ci_bv[2]) * 100, (1 - ci_bv[1]) * 100))
}
boost_var_df <- bind_rows(boost_by_var)

# ============================================================================
# 5. Visualizations
# ============================================================================
results_df <- data.frame(
  method = c("IPTW", "AIPW (Doubly Robust)"),
  rve = c(rve_boost, rve_aipw),
  ci_lo = c((1 - ci_boost[2]) * 100, ci_aipw[1]),
  ci_hi = c((1 - ci_boost[1]) * 100, ci_aipw[2]))

p_boost <- ggplot(results_df, aes(x = method, y = rve)) +
  geom_point(size = 4, color = "#2171B5") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, color = "#2171B5") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Booster Dose Effectiveness (vs 2-dose primary series)",
    subtitle = "Relative VE estimated via IPTW and AIPW",
    x = "", y = "Relative VE (%)") +
  theme_minimal(base_size = 13) + coord_cartesian(ylim = c(-20, 80))
ggsave("figures/booster_effectiveness.png", p_boost, width = 8, height = 5, dpi = 300)

write.csv(results_df, "results/booster_ve_estimates.csv", row.names = FALSE)
write.csv(boost_var_df, "results/booster_ve_by_variant.csv", row.names = FALSE)
cat("\n✓ Booster causal analysis complete.\n")
