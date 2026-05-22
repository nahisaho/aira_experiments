#!/usr/bin/env Rscript

set.seed(42)

run_cytokine_ode_model <- function(base_dir = NULL) {
  message("[03] Starting cytokine ODE network analysis...")

  suppressPackageStartupMessages({
    pkgs <- c("deSolve", "ggplot2", "gridExtra", "tidyr", "dplyr", "RColorBrewer")
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

  params <- c(
    p_IL6 = 0.45, p_TNF = 0.40, p_IL17 = 0.36, p_IL10 = 0.28, p_TGFb = 0.26, p_IFNg = 0.34,
    d_IL6 = 0.16, d_TNF = 0.20, d_IL17 = 0.18, d_IL10 = 0.12, d_TGFb = 0.10, d_IFNg = 0.17,
    a_TNF_IL6 = 1.20, a_TNF_IL17 = 1.10, a_IFNg_TNF = 0.95, a_IFNg_IL6 = 0.72,
    a_IL10_block = 0.90, a_TGFb_block = 0.88, a_TGFb_IL10 = 0.55, a_IL10_Treg = 0.42,
    a_IL17_TNF = 0.60, a_cross_Th1_Th17 = 0.48,
    K_infl = 1.0, K_reg = 0.8, n_hill = 2
  )

  cytokine_model <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      stim_tnf <- p_TNF * (1 + a_IFNg_TNF * hill_act(IFNg, K_infl, n_hill) + a_IL17_TNF * hill_act(IL17A, K_infl, n_hill)) *
        hill_inh(IL10, K_reg, n_hill) * hill_inh(TGFb, K_reg, n_hill)
      stim_il6 <- p_IL6 * (1 + a_TNF_IL6 * hill_act(TNFa, K_infl, n_hill) + a_IFNg_IL6 * hill_act(IFNg, K_infl, n_hill)) *
        hill_inh(IL10, K_reg, n_hill) * hill_inh(TGFb, K_reg, n_hill)
      stim_il17 <- p_IL17 * (1 + a_TNF_IL17 * hill_act(TNFa, K_infl, n_hill) + a_cross_Th1_Th17 * hill_act(IFNg, K_infl, n_hill)) *
        hill_inh(IL10, K_reg, n_hill) * hill_inh(TGFb, K_reg, n_hill)
      stim_il10 <- p_IL10 * (1 + a_TGFb_IL10 * hill_act(TGFb, K_reg, n_hill) + 0.25 * hill_act(IL6, K_infl, n_hill))
      stim_tgfb <- p_TGFb * (1 + a_IL10_Treg * hill_act(IL10, K_reg, n_hill) + 0.18 * hill_inh(TNFa, K_infl, n_hill))
      stim_ifng <- p_IFNg * (1 + 0.70 * hill_act(TNFa, K_infl, n_hill) + 0.42 * hill_act(IL6, K_infl, n_hill)) *
        hill_inh(TGFb, K_reg, n_hill)

      dIL6 <- stim_il6 - d_IL6 * IL6
      dTNFa <- stim_tnf - d_TNF * TNFa
      dIL17A <- stim_il17 - d_IL17 * IL17A
      dIL10 <- stim_il10 - d_IL10 * IL10
      dTGFb <- stim_tgfb - d_TGFb * TGFb
      dIFNg <- stim_ifng - d_IFNg * IFNg

      list(c(dIL6, dTNFa, dIL17A, dIL10, dTGFb, dIFNg))
    })
  }

  scenarios <- list(
    HC_baseline = list(
      init = c(IL6 = 0.35, TNFa = 0.30, IL17A = 0.28, IL10 = 0.75, TGFb = 0.80, IFNg = 0.40),
      pars = c(p_IL6 = 0.40, p_TNF = 0.35, p_IL17 = 0.30, p_IL10 = 0.30, p_TGFb = 0.28, p_IFNg = 0.30)
    ),
    RA_active = list(
      init = c(IL6 = 1.40, TNFa = 1.25, IL17A = 1.10, IL10 = 0.45, TGFb = 0.50, IFNg = 0.95),
      pars = c(p_IL6 = 0.65, p_TNF = 0.70, p_IL17 = 0.60, p_IL10 = 0.24, p_TGFb = 0.22, p_IFNg = 0.48)
    ),
    Anti_TNF_treatment = list(
      init = c(IL6 = 1.20, TNFa = 1.00, IL17A = 0.95, IL10 = 0.52, TGFb = 0.60, IFNg = 0.85),
      pars = c(p_IL6 = 0.52, p_TNF = 0.30, p_IL17 = 0.48, p_IL10 = 0.27, p_TGFb = 0.24, p_IFNg = 0.42)
    ),
    JAK_inhibitor = list(
      init = c(IL6 = 1.10, TNFa = 1.05, IL17A = 0.88, IL10 = 0.55, TGFb = 0.62, IFNg = 0.70),
      pars = c(p_IL6 = 0.36, p_TNF = 0.48, p_IL17 = 0.42, p_IL10 = 0.28, p_TGFb = 0.25, p_IFNg = 0.20)
    )
  )

  simulate_scenario <- function(name, scenario, times = seq(0, 72, by = 0.25)) {
    pars <- utils::modifyList(as.list(params), as.list(scenario$pars))
    out <- as.data.frame(deSolve::ode(y = scenario$init, times = times, func = cytokine_model, parms = pars, method = "lsoda"))
    out$scenario <- name
    out
  }

  sim_list <- lapply(names(scenarios), function(sc) simulate_scenario(sc, scenarios[[sc]]))
  sim_df <- dplyr::bind_rows(sim_list)
  write.csv(sim_df, file.path(res_dir, "cytokine_timecourse.csv"), row.names = FALSE)

  steady_states <- sim_df %>%
    dplyr::group_by(scenario) %>%
    dplyr::slice_max(order_by = time, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      inflammatory_score = (IL6 + TNFa + IL17A + IFNg) / 4,
      regulatory_score = (IL10 + TGFb) / 2,
      inflammatory_ratio = inflammatory_score / regulatory_score
    )
  write.csv(steady_states, file.path(res_dir, "ode_steady_states.csv"), row.names = FALSE)

  numerical_jacobian <- function(state_vec, pars, eps = 1e-6) {
    state_names <- names(state_vec)
    base <- unlist(cytokine_model(0, state_vec, pars)[[1]])
    jac <- matrix(0, nrow = length(state_vec), ncol = length(state_vec), dimnames = list(state_names, state_names))
    for (i in seq_along(state_vec)) {
      perturbed <- state_vec
      perturbed[i] <- perturbed[i] + eps
      deriv <- unlist(cytokine_model(0, perturbed, pars)[[1]])
      jac[, i] <- (deriv - base) / eps
    }
    jac
  }

  eigen_df <- dplyr::bind_rows(lapply(names(scenarios), function(sc) {
    pars <- utils::modifyList(as.list(params), as.list(scenarios[[sc]]$pars))
    state <- as.numeric(steady_states[steady_states$scenario == sc, c("IL6", "TNFa", "IL17A", "IL10", "TGFb", "IFNg")])
    names(state) <- c("IL6", "TNFa", "IL17A", "IL10", "TGFb", "IFNg")
    jac <- numerical_jacobian(state, pars)
    eig <- eigen(jac)$values
    data.frame(
      scenario = sc,
      eigenvalue = seq_along(eig),
      real = Re(eig),
      imaginary = Im(eig),
      stable = max(Re(eig)) < 0,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(eigen_df, file.path(res_dir, "ode_eigenvalues.csv"), row.names = FALSE)

  traj_long <- sim_df %>%
    tidyr::pivot_longer(cols = c(IL6, TNFa, IL17A, IL10, TGFb, IFNg), names_to = "cytokine", values_to = "value")
  timecourse_plot <- ggplot2::ggplot(traj_long, ggplot2::aes(x = time, y = value, color = scenario)) +
    ggplot2::geom_line(size = 0.8) +
    ggplot2::facet_wrap(~ cytokine, scales = "free_y", ncol = 3) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(
      title = "Cytokine network time course across scenarios",
      x = "Time (hours)",
      y = "Concentration (a.u.)",
      color = "Scenario"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(fig_dir, "cytokine_ode_timecourse.pdf"), timecourse_plot, width = 12, height = 8)

  phase_plot <- ggplot2::ggplot(sim_df, ggplot2::aes(x = IL6, y = IL17A, color = scenario)) +
    ggplot2::geom_path(linewidth = 1) +
    ggplot2::geom_point(data = steady_states, size = 2.8) +
    ggplot2::scale_color_brewer(palette = "Dark2") +
    ggplot2::labs(
      title = "Phase portrait of IL-6 and IL-17A dynamics",
      x = "IL-6",
      y = "IL-17A",
      color = "Scenario"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(fig_dir, "cytokine_ode_phase_portrait.pdf"), phase_plot, width = 7.5, height = 6)

  tnf_scan <- seq(0.2, 1.1, length.out = 45)
  bifurcation_df <- dplyr::bind_rows(lapply(tnf_scan, function(ptnf) {
    pars <- as.list(params)
    pars$p_TNF <- ptnf
    out <- as.data.frame(deSolve::ode(
      y = scenarios$RA_active$init,
      times = seq(0, 240, by = 1),
      func = cytokine_model,
      parms = pars,
      method = "lsoda"
    ))
    ss <- out[nrow(out), c("IL6", "TNFa", "IL17A", "IL10", "TGFb", "IFNg")]
    data.frame(p_TNF = ptnf, cytokine = names(ss), steady_state = as.numeric(ss), stringsAsFactors = FALSE)
  }))
  write.csv(bifurcation_df, file.path(res_dir, "ode_bifurcation_scan.csv"), row.names = FALSE)
  bifurcation_plot <- ggplot2::ggplot(bifurcation_df, ggplot2::aes(x = p_TNF, y = steady_state, color = cytokine)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_color_brewer(palette = "Paired") +
    ggplot2::labs(
      title = "Bifurcation scan across TNF production rates",
      x = "TNF production rate",
      y = "Steady-state abundance",
      color = "Cytokine"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(fig_dir, "cytokine_ode_bifurcation.pdf"), bifurcation_plot, width = 10, height = 6)

  stability_plot <- ggplot2::ggplot(eigen_df, ggplot2::aes(x = real, y = imaginary, color = scenario)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 3, color = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = 3, color = "grey50") +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_color_brewer(palette = "Set2") +
    ggplot2::labs(
      title = "Eigenvalue stability analysis",
      x = "Real component",
      y = "Imaginary component",
      color = "Scenario"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(fig_dir, "cytokine_ode_stability.pdf"), stability_plot, width = 7.5, height = 6)

  message("[03] Cytokine ODE network analysis completed.")
  invisible(list(simulation = sim_df, steady_states = steady_states, eigen = eigen_df, bifurcation = bifurcation_df))
}

if (sys.nframe() == 0L) {
  run_cytokine_ode_model()
}
