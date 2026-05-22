#!/usr/bin/env Rscript
# ============================================================================
# 07_hospitalization_case_study.R — mRNAワクチンの入院予防効果評価
# ============================================================================
source("R/00_setup.R")
library(dplyr); library(survival); library(mgcv); library(splines)
library(sandwich); library(lmtest); library(ggplot2); library(patchwork)

tnd_data <- readRDS("data/tnd_simulated.rds")

# ============================================================================
# 1. 入院予防VE (全体)
# ============================================================================
cat("=== mRNA Vaccine Effectiveness Against Hospitalization ===\n")

# 入院を対象としたTND: COVID入院 vs non-COVID入院
df_hosp <- tnd_data %>% filter(hospitalized == 1)
cat(sprintf("Hospitalized individuals: N=%d\n", nrow(df_hosp)))

m_hosp <- glm(
  tnd_case ~ vaccinated + age + I(age^2) + sex + comorbidity +
    prior_infection + factor(calendar_week),
  family = binomial, data = df_hosp)

or_hosp <- exp(coef(m_hosp)["vaccinated"])
ci_hosp <- exp(confint.default(m_hosp)["vaccinated", ])
ve_hosp <- (1 - or_hosp) * 100
cat(sprintf("VE against hospitalization: %.1f%% (95%% CI: %.1f%% - %.1f%%)\n",
  ve_hosp, (1 - ci_hosp[2]) * 100, (1 - ci_hosp[1]) * 100))

# ============================================================================
# 2. 接種回数別の入院VE
# ============================================================================
cat("\n=== VE by Dose (Hospitalization) ===\n")

df_hosp$dose_f <- factor(df_hosp$dose, levels = c(0, 1, 2, 3),
  labels = c("Unvaccinated", "Dose 1", "Dose 2", "Booster"))

m_hosp_dose <- glm(
  tnd_case ~ dose_f + age + I(age^2) + sex + comorbidity +
    prior_infection + factor(calendar_week),
  family = binomial, data = df_hosp)

dose_coefs <- summary(m_hosp_dose)$coefficients
dose_labels <- c("Dose 1", "Dose 2", "Booster")
ve_dose <- data.frame(
  dose = dose_labels,
  or = exp(dose_coefs[2:4, "Estimate"]),
  ve = (1 - exp(dose_coefs[2:4, "Estimate"])) * 100,
  ci_lo = (1 - exp(dose_coefs[2:4, "Estimate"] + 1.96 * dose_coefs[2:4, "Std. Error"])) * 100,
  ci_hi = (1 - exp(dose_coefs[2:4, "Estimate"] - 1.96 * dose_coefs[2:4, "Std. Error"])) * 100)
print(ve_dose)

# ============================================================================
# 3. 変異株別の入院VE
# ============================================================================
cat("\n=== VE Against Hospitalization by Variant ===\n")

ve_hosp_var <- list()
for (v in unique(df_hosp$variant)) {
  df_hv <- df_hosp %>% filter(variant == v)
  if (nrow(df_hv) < 50) next
  m_hv <- glm(tnd_case ~ vaccinated + age + sex + comorbidity + prior_infection,
    family = binomial, data = df_hv)
  or_hv <- exp(coef(m_hv)["vaccinated"])
  ci_hv <- exp(confint.default(m_hv)["vaccinated", ])
  ve_hosp_var[[v]] <- data.frame(variant = v,
    ve = (1 - or_hv) * 100,
    ci_lo = (1 - ci_hv[2]) * 100, ci_hi = (1 - ci_hv[1]) * 100,
    n = nrow(df_hv))
  cat(sprintf("  %s: VE=%.1f%% (%.1f%%-%.1f%%), N=%d\n",
    v, (1 - or_hv) * 100, (1 - ci_hv[2]) * 100, (1 - ci_hv[1]) * 100, nrow(df_hv)))
}
ve_hosp_var_df <- bind_rows(ve_hosp_var)

# ============================================================================
# 4. 入院VEのWaning
# ============================================================================
cat("\n=== Hospitalization VE Waning ===\n")

df_hosp_vax <- df_hosp %>% filter(vaccinated == 1, !is.na(days_since_vax))
m_hosp_wan <- gam(
  tnd_case ~ s(days_since_vax, bs = "cr", k = 8) +
    s(age, bs = "cr", k = 5) + sex + comorbidity + prior_infection,
  family = binomial, data = df_hosp_vax, method = "REML")

pred_h <- data.frame(days_since_vax = seq(14, 300), age = 60,
  sex = 0, comorbidity = 0, prior_infection = 0)
pred_hp <- predict(m_hosp_wan, newdata = pred_h, type = "link", se.fit = TRUE)
ref_hp <- pred_hp$fit[1]
pred_h$ve <- (1 - exp(pred_hp$fit - ref_hp)) * 100
pred_h$ve_lo <- (1 - exp(pred_hp$fit - ref_hp + 1.96 * pred_hp$se.fit)) * 100
pred_h$ve_hi <- (1 - exp(pred_hp$fit - ref_hp - 1.96 * pred_hp$se.fit)) * 100

# ============================================================================
# 5. 年齢層別の入院VE
# ============================================================================
cat("\n=== VE by Age Group (Hospitalization) ===\n")

df_hosp$age_group <- cut(df_hosp$age,
  breaks = c(17, 39, 59, 74, 96),
  labels = c("18-39", "40-59", "60-74", "75+"))

ve_by_age <- list()
for (ag in levels(df_hosp$age_group)) {
  df_ag <- df_hosp %>% filter(age_group == ag)
  if (nrow(df_ag) < 30) next
  m_ag <- glm(tnd_case ~ vaccinated + sex + comorbidity + prior_infection,
    family = binomial, data = df_ag)
  or_ag <- exp(coef(m_ag)["vaccinated"])
  ci_ag <- exp(confint.default(m_ag)["vaccinated", ])
  ve_by_age[[ag]] <- data.frame(age_group = ag,
    ve = (1 - or_ag) * 100,
    ci_lo = (1 - ci_ag[2]) * 100, ci_hi = (1 - ci_ag[1]) * 100,
    n = nrow(df_ag))
}
ve_age_df <- bind_rows(ve_by_age)
print(ve_age_df)

# ============================================================================
# 6. Sensitivity Analyses
# ============================================================================
cat("\n=== Sensitivity Analyses ===\n")

# 6a. Excluding 0-13 days post-vax (partially vaccinated)
df_sens1 <- df_hosp %>% filter(is.na(days_since_vax) | days_since_vax >= 14)
m_sens1 <- glm(tnd_case ~ vaccinated + age + sex + comorbidity + prior_infection,
  family = binomial, data = df_sens1)
ve_sens1 <- (1 - exp(coef(m_sens1)["vaccinated"])) * 100
cat(sprintf("  Excluding <14d: VE = %.1f%%\n", ve_sens1))

# 6b. Restricting to ≥65 years
df_sens2 <- df_hosp %>% filter(age >= 65)
m_sens2 <- glm(tnd_case ~ vaccinated + age + sex + comorbidity + prior_infection,
  family = binomial, data = df_sens2)
ve_sens2 <- (1 - exp(coef(m_sens2)["vaccinated"])) * 100
cat(sprintf("  Age ≥65 only:    VE = %.1f%%\n", ve_sens2))

# 6c. Restricting to comorbid patients
df_sens3 <- df_hosp %>% filter(comorbidity == 1)
m_sens3 <- glm(tnd_case ~ vaccinated + age + sex + prior_infection,
  family = binomial, data = df_sens3)
ve_sens3 <- (1 - exp(coef(m_sens3)["vaccinated"])) * 100
cat(sprintf("  Comorbid only:   VE = %.1f%%\n", ve_sens3))

sensitivity <- data.frame(
  analysis = c("Primary", "Exclude <14d", "Age ≥65", "Comorbid only"),
  ve = c(ve_hosp, ve_sens1, ve_sens2, ve_sens3))

# ============================================================================
# 7. Comprehensive Figures
# ============================================================================

# Panel A: VE by dose
p_dose <- ggplot(ve_dose, aes(x = dose, y = ve)) +
  geom_point(size = 4, color = "#2171B5") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, color = "#2171B5") +
  labs(title = "A) VE Against Hospitalization by Dose", x = "", y = "VE (%)") +
  theme_minimal(base_size = 12) + coord_cartesian(ylim = c(0, 100))

# Panel B: VE waning
p_wan <- ggplot(pred_h, aes(x = days_since_vax)) +
  geom_ribbon(aes(ymin = ve_lo, ymax = ve_hi), fill = "#2171B5", alpha = 0.2) +
  geom_line(aes(y = ve), color = "#2171B5", linewidth = 1) +
  labs(title = "B) Hospitalization VE Waning", x = "Days Since Vaccination",
    y = "Relative VE Change (%)") +
  theme_minimal(base_size = 12)

# Panel C: VE by age group
p_age <- ggplot(ve_age_df, aes(x = age_group, y = ve)) +
  geom_point(size = 4, color = "#CB181D") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, color = "#CB181D") +
  labs(title = "C) VE by Age Group", x = "Age Group", y = "VE (%)") +
  theme_minimal(base_size = 12) + coord_cartesian(ylim = c(0, 100))

# Panel D: Sensitivity
p_sens <- ggplot(sensitivity, aes(x = reorder(analysis, ve), y = ve)) +
  geom_point(size = 4, color = "#238B45") +
  coord_flip() +
  labs(title = "D) Sensitivity Analyses", x = "", y = "VE (%)") +
  theme_minimal(base_size = 12) + coord_flip(ylim = c(0, 100))

p_combined <- (p_dose | p_wan) / (p_age | p_sens) +
  plot_annotation(
    title = "mRNA Vaccine Effectiveness Against COVID-19 Hospitalization",
    subtitle = "Test-Negative Design, Simulated Real-World Data (N=50,000)",
    theme = theme(plot.title = element_text(size = 16, face = "bold")))

ggsave("figures/hospitalization_case_study.png", p_combined,
  width = 14, height = 10, dpi = 300)
ggsave("figures/hospitalization_case_study.svg", p_combined,
  width = 14, height = 10)

# Save all results
write.csv(ve_dose, "results/hosp_ve_by_dose.csv", row.names = FALSE)
write.csv(ve_hosp_var_df, "results/hosp_ve_by_variant.csv", row.names = FALSE)
write.csv(ve_age_df, "results/hosp_ve_by_age.csv", row.names = FALSE)
write.csv(sensitivity, "results/hosp_sensitivity.csv", row.names = FALSE)
write.csv(pred_h, "results/hosp_waning_curve.csv", row.names = FALSE)

cat("\n✓ Hospitalization case study complete.\n")
