#!/usr/bin/env Rscript
# ============================================================================
# 03_waning_model.R — ワクチン効果減衰 (Waning) の推定
# ============================================================================
source("R/00_setup.R")
library(dplyr); library(survival); library(mgcv)
library(splines); library(ggplot2); library(patchwork)

tnd_data <- readRDS("data/tnd_simulated.rds")
df_vax <- tnd_data %>% filter(vaccinated == 1, !is.na(days_since_vax))

# --- 1. Piecewise constant VE ---
cat("=== Piecewise VE ===\n")
m_piece <- clogit(
  tnd_case ~ time_since_vax_cat + age + sex + comorbidity +
    prior_infection + strata(calendar_week), data = df_vax)
coefs_piece <- summary(m_piece)$coefficients
ve_piece <- data.frame(
  interval = gsub("time_since_vax_cat", "", rownames(coefs_piece)[1:4]),
  log_or = coefs_piece[1:4, "coef"], se = coefs_piece[1:4, "se(coef)"])
ve_piece$ve     <- (1 - exp(ve_piece$log_or)) * 100
ve_piece$ci_low <- (1 - exp(ve_piece$log_or + 1.96 * ve_piece$se)) * 100
ve_piece$ci_up  <- (1 - exp(ve_piece$log_or - 1.96 * ve_piece$se)) * 100
print(ve_piece[, c("interval", "ve", "ci_low", "ci_up")])

# --- 2. Natural spline ---
m_spline <- glm(
  tnd_case ~ ns(days_since_vax, df = 4) + age + sex + comorbidity +
    prior_infection + factor(calendar_week), family = binomial, data = df_vax)
pred_grid <- data.frame(days_since_vax = seq(14, 365), age = median(df_vax$age),
  sex = 0, comorbidity = 0, prior_infection = 0, calendar_week = 26)
pred_lp <- predict(m_spline, newdata = pred_grid, type = "link", se.fit = TRUE)
ref_lp <- pred_lp$fit[1]
pred_grid$ve <- (1 - exp(pred_lp$fit - ref_lp)) * 100

# --- 3. GAM ---
cat("\n=== GAM Waning ===\n")
m_gam <- gam(
  tnd_case ~ s(days_since_vax, bs = "cr", k = 10) +
    s(age, bs = "cr", k = 5) + sex + comorbidity + prior_infection +
    s(calendar_week, bs = "cc", k = 12),
  family = binomial, data = df_vax, method = "REML")
print(summary(m_gam))
pred_gam <- predict(m_gam, newdata = pred_grid, type = "link", se.fit = TRUE)
ref_gam <- pred_gam$fit[1]
pred_grid$ve_gam    <- (1 - exp(pred_gam$fit - ref_gam)) * 100
pred_grid$ve_gam_lo <- (1 - exp(pred_gam$fit - ref_gam + 1.96 * pred_gam$se.fit)) * 100
pred_grid$ve_gam_hi <- (1 - exp(pred_gam$fit - ref_gam - 1.96 * pred_gam$se.fit)) * 100

# --- 4. Exponential decay half-life ---
m_exp <- glm(
  tnd_case ~ days_since_vax + age + sex + comorbidity +
    prior_infection + factor(calendar_week), family = binomial, data = df_vax)
lambda <- coef(m_exp)["days_since_vax"]
half_life <- log(2) / abs(lambda)
cat(sprintf("Half-life: %.1f days (%.1f months)\n", half_life, half_life / 30))

# --- Plots ---
p_waning <- ggplot(pred_grid, aes(x = days_since_vax)) +
  geom_ribbon(aes(ymin = ve_gam_lo, ymax = ve_gam_hi), fill = "#2171B5", alpha = 0.2) +
  geom_line(aes(y = ve_gam, color = "GAM"), linewidth = 1) +
  geom_line(aes(y = ve, color = "Natural Spline"), linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = c("GAM" = "#2171B5", "Natural Spline" = "#CB181D")) +
  labs(title = "Vaccine Effectiveness Waning Over Time",
    x = "Days Since Vaccination", y = "Relative VE Change (%)", color = "Method") +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")
ggsave("figures/waning_curve.png", p_waning, width = 10, height = 6, dpi = 300)
ggsave("figures/waning_curve.svg", p_waning, width = 10, height = 6)

p_piece <- ggplot(ve_piece, aes(x = interval, y = ve)) +
  geom_col(fill = "#2171B5", alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_up), width = 0.2) +
  labs(title = "VE by Time Interval", x = "Time Since Vaccination", y = "VE (%)") +
  theme_minimal(base_size = 13) + coord_cartesian(ylim = c(0, 100))
ggsave("figures/ve_piecewise.png", p_piece, width = 8, height = 5, dpi = 300)

write.csv(ve_piece, "results/waning_piecewise.csv", row.names = FALSE)
write.csv(pred_grid, "results/waning_continuous.csv", row.names = FALSE)
write.csv(data.frame(lambda = lambda, half_life_days = half_life,
  half_life_months = half_life / 30), "results/waning_halflife.csv", row.names = FALSE)
cat("\n✓ Waning analysis complete.\n")
