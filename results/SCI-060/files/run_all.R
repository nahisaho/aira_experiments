#!/usr/bin/env Rscript
# ============================================================================
# run_all.R — 全パイプライン実行スクリプト
# ============================================================================
# Usage: Rscript run_all.R
#
# 前提: R 4.x + required packages (00_setup.R で自動インストール)
# 出力: data/, results/, figures/, logs/ に全成果物を保存
# ============================================================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  VE Estimation Framework — Real-World Data Pipeline         ║\n")
cat("║  Version 1.0 | 2026-05-23                                   ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

start_time <- Sys.time()

scripts <- c(
  "R/01_simulate_data.R",
  "R/02_tnd_analysis.R",
  "R/03_waning_model.R",
  "R/04_variant_specific_ve.R",
  "R/05_healthy_vaccinee_bias.R",
  "R/06_booster_causal.R",
  "R/07_hospitalization_case_study.R"
)

for (s in scripts) {
  cat(sprintf("\n{'='*60}\n▶ Running %s\n{'='*60}\n", s))
  tryCatch(
    source(s, local = new.env()),
    error = function(e) {
      cat(sprintf("✗ ERROR in %s: %s\n", s, e$message))
    }
  )
}

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("\n✓ Pipeline complete in %.1f minutes.\n", as.numeric(elapsed)))

# Process log
log_entry <- sprintf(
  '{"timestamp":"%s","phase":"pipeline","event_type":"run_completed","elapsed_min":%.1f}\n',
  format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), as.numeric(elapsed))
cat(log_entry, file = "logs/process-log.jsonl", append = TRUE)
