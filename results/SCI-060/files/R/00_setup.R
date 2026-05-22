#!/usr/bin/env Rscript
# ============================================================================
# 00_setup.R — パッケージインストールと環境設定
# VE Estimation Framework for Real-World Data
# ============================================================================

required_packages <- c(
  "survival", "gnm", "lme4", "mgcv", "splines",
  "MatchIt", "WeightIt", "cobalt", "sandwich", "lmtest",
  "dplyr", "tidyr", "ggplot2", "patchwork", "scales",
  "lubridate", "broom",
  "EValue", "timereg",
  "knitr", "rmarkdown", "gt"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org", quiet = TRUE)
  }
}
invisible(lapply(required_packages, install_if_missing))

cat("✓ All packages installed/verified.\n")

set.seed(20260523)
options(mc.cores = parallel::detectCores() - 1, warn = 1, scipen = 999, digits = 4)
