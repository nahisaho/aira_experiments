#!/usr/bin/env Rscript

set.seed(42)

run_tolerance_insilico <- function(base_dir = NULL) {
  message("[06] Starting in silico tolerance restoration analysis...")

  suppressPackageStartupMessages({
    pkgs <- c("deSolve", "ggplot2", "gridExtra", "reshape2", "tidyr", "dplyr", "viridis")
    missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      stop(sprintf("Missing required packages: %s", paste(missing_pkgs, collapse = ", ")))
    }
    invisible(lapply(pkgs, library, character.only = TRUE))
  })

  resolve_base_dir <- function(path_hint = NULL) {
    if (!is.null(path_hint)) return(normalizePath(path_hint, mustWork = TRUE))
    wd <- normalizePath(getwd(), mustWork = TRUE)
    if (basename(wd) == "scripts") dirname(wd) else wd
  }

  base_dir <- resolve_base_dir(base_dir)
  fig_dir <- file.path(base_dir, "figures")
  res_dir <- file.path(base_dir, "results")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

  hill_act <- function(x, k, n = 2) (x^n) / (k^n + x^n)
  hill_inh <- function(x, k, n = 2) 1 / (1 + (x / k)^n)

  base_params <- c(
    p_Teff = 0.40, p_Treg = 0.22, p_APCact = 0.36, p_APCtol = 0.18, p_IL10 = 0.25,
    p_TGFb = 0.23, p_IL6 = 0.44, p_TNFa = 0.42, p_IDO1 = 0.16, p_PDL1 = 0.20,
    d_Teff = 0.12, d_Treg = 0.08, d_APCact = 0.10, d_APCtol = 0.08, d_IL10 = 0.18,
    d_TGFb = 0.15, d_IL6 = 0.25, d_TNFa = 0.24, d_IDO1 = 0.12, d_PDL1 = 0.12,
    k_activation = 1.0, k_reg = 0.8, n_hill = 2,
    s_treg_supp = 1.10, s_tol_reprogram = 0.95, s_checkpoint = 0.90, s_ido = 0.85
  )

  tolerance_model <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      eff_drive <- p_Teff * (1 + 1.10 * hill_act(APC_active + IL6 + TNFa, k_activation, n_hill)) *
        hill_inh(Treg + IL10 + TGFb + PDL1 + IDO1, k_reg, n_hill)
      treg_drive <- p_Treg * (1 + 0.95 * hill_act(APC_tolerogenic + IL10 + TGFb, k_reg, n_hill))
      apc_active_drive <- p_APCact * (1 + 0.75 * hill_act(TNFa + IL6, k_activation, n_hill)) * hill_inh(IL10 + TGFb + IDO1 + PDL1, k_reg, n_hill)
      apc_tol_drive <- p_APCtol * (1 + s_tol_reprogram * hill_act(IL10 + TGFb + IDO1 + PDL1, k_reg, n_hill))
      il10_drive <- p_IL10 * (1 + 0.7 * hill_act(Treg + APC_tolerogenic, k_reg, n_hill))
      tgfb_drive <- p_TGFb * (1 + 0.6 * hill_act(Treg + APC_tolerogenic, k_reg, n_hill))
      il6_drive <- p_IL6 * (1 + hill_act(Teff + APC_active + TNFa, k_activation, n_hill)) * hill_inh(IL10 + TGFb + PDL1, k_reg, n_hill)
      tnf_drive <- p_TNFa * (1 + hill_act(Teff + APC_active + IL6, k_activation, n_hill)) * hill_inh(IL10 + TGFb + IDO1, k_reg, n_hill)
      ido_drive <- p_IDO1 * (1 + s_ido * hill_act(APC_tolerogenic + IL10 + TGFb, k_reg, n_hill))
      pdl1_drive <- p_PDL1 * (1 + s_checkpoint * hill_act(APC_tolerogenic + IL10 + IDO1, k_reg, n_hill))

      dTeff <- eff_drive - d_Teff * Teff
      dTreg <- treg_drive - d_Treg * Treg
      dAPC_active <- apc_active_drive - d_APCact * APC_active - 0.10 * hill_act(IL10 + TGFb + IDO1 + PDL1, k_reg, n_hill) * APC_active
      dAPC_tolerogenic <- apc_tol_drive - d_APCtol * APC_tolerogenic + 0.10 * hill_act(IL10 + TGFb + IDO1 + PDL1, k_reg, n_hill) * APC_active
      dIL10 <- il10_drive - d_IL10 * IL10
      dTGFb <- tgfb_drive - d_TGFb * TGFb
      dIL6 <- il6_drive - d_IL6 * IL6
      dTNFa <- tnf_drive - d_TNFa * TNFa
      dIDO1 <- ido_drive - d_IDO1 * IDO1
      dPDL1 <- pdl1_drive - d_PDL1 * PDL1

      list(c(dTeff, dTreg, dAPC_active, dAPC_tolerogenic, dIL10, dTGFb, dIL6, dTNFa, dIDO1, dPDL1))
    })
  }

  baseline_state <- c(
    Teff = 1.8, Treg = 0.6, APC_active = 1.3, APC_tolerogenic = 0.4,
    IL10 = 0.5, TGFb = 0.55, IL6 = 1.4, TNFa = 1.3, IDO1 = 0.35, PDL1 = 0.40
  )

  strategies <- list(
    `Treg adoptive transfer` = list(init = c(Treg = 1.5), pars = c(p_Treg = 0.28)),
    `Anti-inflammatory cytokines` = list(init = c(IL10 = 1.2, TGFb = 1.3), pars = c(p_IL10 = 0.32, p_TGFb = 0.30)),
    `Checkpoint agonist` = list(init = c(PDL1 = 1.1), pars = c(p_PDL1 = 0.36, s_checkpoint = 1.30)),
    `IDO1 induction` = list(init = c(IDO1 = 1.0), pars = c(p_IDO1 = 0.32, s_ido = 1.25)),
    `APC tolerogenic reprogramming` = list(init = c(APC_tolerogenic = 1.0, APC_active = 0.9), pars = c(p_APCtol = 0.30, s_tol_reprogram = 1.20)),
    `Combination therapy` = list(init = c(Treg = 1.4, IL10 = 1.2, TGFb = 1.1, PDL1 = 1.0), pars = c(p_Treg = 0.30, p_IL10 = 0.32, p_TGFb = 0.30, p_PDL1 = 0.34, s_checkpoint = 1.20))
  )

  simulate_strategy <- function(name, strategy, pars_override = NULL, times = seq(0, 30, by = 0.25)) {
    init <- baseline_state
    init[names(strategy$init)] <- strategy$init
    pars <- as.list(base_params)
    pars[names(strategy$pars)] <- strategy$pars
    if (!is.null(pars_override)) pars[names(pars_override)] <- pars_override
    out <- as.data.frame(deSolve::ode(y = init, times = times, func = tolerance_model, parms = pars, method = "lsoda"))
    out$strategy <- name
    out
  }

  traj_df <- dplyr::bind_rows(lapply(names(strategies), function(st) simulate_strategy(st, strategies[[st]])))
  write.csv(traj_df, file.path(res_dir, "tolerance_strategy_timecourse.csv"), row.names = FALSE)

  efficacy_df <- traj_df %>%
    dplyr::group_by(strategy) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Teff_Treg_ratio = Teff / Treg,
      inflammatory_score = (IL6 + TNFa + APC_active) / 3,
      tolerance_score = (Treg + APC_tolerogenic + IL10 + TGFb + IDO1 + PDL1) / 6 - inflammatory_score,
      composite_rank_score = scales::rescale(-Teff_Treg_ratio) + scales::rescale(-inflammatory_score) + scales::rescale(tolerance_score)
    ) %>%
    dplyr::arrange(dplyr::desc(composite_rank_score))
  write.csv(efficacy_df, file.path(res_dir, "tolerance_efficacy.csv"), row.names = FALSE)

  traj_long <- traj_df %>% tidyr::pivot_longer(cols = c(Teff, Treg, APC_active, APC_tolerogenic, IL10, TGFb, IL6, TNFa, IDO1, PDL1), names_to = "state", values_to = "value")
  traj_plot <- ggplot2::ggplot(traj_long %>% dplyr::filter(state %in% c("Teff", "Treg", "IL6", "TNFa", "IDO1", "PDL1")), ggplot2::aes(x = time, y = value, color = strategy)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~ state, scales = "free_y", ncol = 3) +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::labs(title = "Tolerance restoration trajectories", x = "Time (days)", y = "Abundance (a.u.)", color = "Strategy") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(fig_dir, "tolerance_dynamics.pdf"), traj_plot, width = 12, height = 8)

  ranking_plot <- ggplot2::ggplot(efficacy_df, ggplot2::aes(x = reorder(strategy, composite_rank_score), y = composite_rank_score, fill = tolerance_score)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_viridis_c(option = "plasma") +
    ggplot2::labs(title = "Composite efficacy of tolerance restoration strategies", x = "Strategy", y = "Composite efficacy score", fill = "Tolerance score") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(fig_dir, "tolerance_strategy_ranking.pdf"), ranking_plot, width = 10, height = 6)

  perturb_families <- c("Inflammatory production", "Regulatory production", "Checkpoint strength", "Cell turnover")
  sensitivity_runs <- lapply(seq_len(1000), function(run_id) {
    family <- sample(perturb_families, 1)
    pars_override <- numeric(0)
    if (family == "Inflammatory production") {
      scale_factor <- runif(1, 0.7, 1.3)
      pars_override <- c(p_IL6 = base_params["p_IL6"] * scale_factor, p_TNFa = base_params["p_TNFa"] * scale_factor, p_Teff = base_params["p_Teff"] * scale_factor)
    }
    if (family == "Regulatory production") {
      scale_factor <- runif(1, 0.7, 1.3)
      pars_override <- c(p_Treg = base_params["p_Treg"] * scale_factor, p_IL10 = base_params["p_IL10"] * scale_factor, p_TGFb = base_params["p_TGFb"] * scale_factor)
    }
    if (family == "Checkpoint strength") {
      scale_factor <- runif(1, 0.7, 1.3)
      pars_override <- c(p_PDL1 = base_params["p_PDL1"] * scale_factor, p_IDO1 = base_params["p_IDO1"] * scale_factor, s_checkpoint = base_params["s_checkpoint"] * scale_factor)
    }
    if (family == "Cell turnover") {
      scale_factor <- runif(1, 0.7, 1.3)
      pars_override <- c(d_Teff = base_params["d_Teff"] * scale_factor, d_Treg = base_params["d_Treg"] / scale_factor, d_APCact = base_params["d_APCact"] * scale_factor)
    }

    dplyr::bind_rows(lapply(names(strategies), function(st) {
      out <- simulate_strategy(st, strategies[[st]], pars_override = pars_override, times = seq(0, 30, by = 0.5))
      ss <- out[nrow(out), ]
      data.frame(
        run_id = run_id,
        perturbation_family = family,
        strategy = st,
        Teff_Treg_ratio = ss$Teff / ss$Treg,
        inflammatory_score = (ss$IL6 + ss$TNFa + ss$APC_active) / 3,
        tolerance_score = (ss$Treg + ss$APC_tolerogenic + ss$IL10 + ss$TGFb + ss$IDO1 + ss$PDL1) / 6 - ((ss$IL6 + ss$TNFa + ss$APC_active) / 3),
        stringsAsFactors = FALSE
      )
    }))
  })
  sensitivity_df <- dplyr::bind_rows(sensitivity_runs)
  write.csv(sensitivity_df, file.path(res_dir, "sensitivity_analysis.csv"), row.names = FALSE)

  heatmap_df <- sensitivity_df %>%
    dplyr::group_by(strategy, perturbation_family) %>%
    dplyr::summarise(mean_tolerance = mean(tolerance_score), robustness = mean(tolerance_score > 0), .groups = "drop")
  heatmap_plot <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = perturbation_family, y = strategy, fill = mean_tolerance)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_viridis_c(option = "magma") +
    ggplot2::labs(title = "Robustness of tolerance strategies across perturbations", x = "Parameter perturbation", y = "Strategy", fill = "Mean tolerance") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  ggplot2::ggsave(file.path(fig_dir, "tolerance_sensitivity_heatmap.pdf"), heatmap_plot, width = 10, height = 6)

  message("[06] In silico tolerance restoration analysis completed.")
  invisible(list(trajectories = traj_df, efficacy = efficacy_df, sensitivity = sensitivity_df))
}

if (sys.nframe() == 0L) {
  run_tolerance_insilico()
}
