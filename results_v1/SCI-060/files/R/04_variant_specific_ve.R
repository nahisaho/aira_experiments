#!/usr/bin/env Rscript
# ============================================================================
# 04_variant_specific_ve.R — 変異株特異的VE推定
# ============================================================================
source("R/00_setup.R")
library(dplyr); library(survival); library(gnm); library(ggplot2); library(broom)
library(splines)

tnd_data <- readRDS("data/tnd_simulated.rds")

# --- 1. Stratified VE by variant ---
cat("=== Variant-Specific VE ===\n")
variants <- unique(tnd_data$variant)
ve_by_variant <- list()
for (v in variants) {
  df_v <- tnd_data %>% filter(variant == v)
  m_v <- glm(tnd_case ~ vaccinated + age + sex + comorbidity + prior_infection,
    family = binomial, data = df_v)
  or <- exp(coef(m_v)["vaccinated"])
  ci <- exp(confint.default(m_v)["vaccinated", ])
  ve <- (1 - or) * 100
  ve_by_variant[[v]] <- data.frame(variant = v, n = nrow(df_v), n_cases = sum(df_v$tnd_case),
    or = or, ve = ve, ci_lower = (1 - ci[2]) * 100, ci_upper = (1 - ci[1]) * 100)
  cat(sprintf("  %s: VE=%.1f%% (%.1f%%-%.1f%%), N=%d\n",
    v, ve, (1 - ci[2]) * 100, (1 - ci[1]) * 100, nrow(df_v)))
}
ve_variant_df <- bind_rows(ve_by_variant)

# --- 2. Interaction model ---
cat("\n=== Interaction: Vaccine × Variant ===\n")
tnd_data$variant_f <- relevel(factor(tnd_data$variant), ref = "Alpha")
m_interact <- glm(
  tnd_case ~ vaccinated * variant_f + age + sex + comorbidity +
    prior_infection + factor(calendar_week), family = binomial, data = tnd_data)
m_no_interact <- glm(
  tnd_case ~ vaccinated + variant_f + age + sex + comorbidity +
    prior_infection + factor(calendar_week), family = binomial, data = tnd_data)
lr_test <- anova(m_no_interact, m_interact, test = "LRT")
cat(sprintf("LR test: χ²=%.2f, p=%.2e\n", lr_test$Deviance[2], lr_test$`Pr(>Chi)`[2]))

# --- 3. Combined variant × waning ---
df_vax_v <- tnd_data %>% filter(vaccinated == 1, !is.na(days_since_vax))
m_combined <- glm(
  tnd_case ~ variant_f * ns(days_since_vax, df = 3) +
    age + sex + comorbidity + prior_infection,
  family = binomial, data = df_vax_v)
cat("AIC (combined):", AIC(m_combined), "\n")

# --- Plot ---
ve_variant_df$variant <- factor(ve_variant_df$variant,
  levels = c("Alpha", "Delta", "Omicron_BA1", "Omicron_BA5"))
p <- ggplot(ve_variant_df, aes(x = variant, y = ve)) +
  geom_point(size = 4, color = "#2171B5") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "#2171B5") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Vaccine Effectiveness by SARS-CoV-2 Variant", x = "Variant", y = "VE (%)") +
  theme_minimal(base_size = 13) + coord_cartesian(ylim = c(-10, 100))
ggsave("figures/ve_by_variant.png", p, width = 9, height = 6, dpi = 300)
ggsave("figures/ve_by_variant.svg", p, width = 9, height = 6)

write.csv(ve_variant_df, "results/ve_by_variant.csv", row.names = FALSE)
cat("\n✓ Variant-specific VE analysis complete.\n")
