#!/usr/bin/env Rscript

# R-INLA / SPDE workflow template for Bayesian disease mapping
# ------------------------------------------------------------
# This template shows a complete end-to-end workflow for fitting a
# spatial Poisson disease model using the SPDE approach in R-INLA.
# It is intentionally written as a well-commented template: adapt the
# file paths, response variable names, and priors to your own study.

# 1) Load libraries ----------------------------------------------------------
library(INLA)     # Core INLA / SPDE functionality
library(sp)       # Spatial classes and coordinate handling
library(spdep)    # Neighbour structures / optional spatial diagnostics
library(ggplot2)  # Plotting

# 2) Prepare disease mapping data -------------------------------------------
# Example data frame structure. Replace this with your observed data import.
# disease_df <- read.csv("data/disease_mapping_input.csv")
set.seed(2025)
n_obs <- 180

disease_df <- data.frame(
  area_id = seq_len(n_obs),
  x = runif(n_obs, min = 0, max = 100),
  y = runif(n_obs, min = 0, max = 100),
  deprivation = rnorm(n_obs, mean = 0, sd = 1),
  expected = runif(n_obs, min = 20, max = 120)
)

# Simulated observed counts for template purposes only.
# In practice, `cases` would come from surveillance or registry data.
linear_predictor <- -2.1 + 0.25 * disease_df$deprivation
risk <- exp(linear_predictor)
disease_df$cases <- rpois(n_obs, lambda = disease_df$expected * risk)

coordinates(disease_df) <- ~ x + y
proj4string(disease_df) <- CRS("+proj=utm +zone=33 +datum=WGS84 +units=km +no_defs")

# Optional neighbour structure for exploratory diagnostics.
coords_matrix <- coordinates(disease_df)
knn_graph <- knearneigh(coords_matrix, k = 4)
nb_object <- knn2nb(knn_graph)
listw_object <- nb2listw(nb_object, style = "W")
# moran.test(disease_df$cases / disease_df$expected, listw_object)

# 3) Build the SPDE mesh ----------------------------------------------------
# max.edge and cutoff should reflect the scale of the coordinate system.
# Smaller values -> finer mesh, higher computational cost.
mesh <- inla.mesh.2d(
  loc = coords_matrix,
  max.edge = c(5, 12),
  cutoff = 1.5,
  offset = c(5, 10)
)

# plot(mesh)
# points(coords_matrix, col = "red", pch = 16, cex = 0.6)

# 4) Specify the Matérn SPDE model -----------------------------------------
# Penalised complexity priors: adjust to disease scale and scientific context.
# prior.range = c(r0, p0) means P(range < r0) = p0
# prior.sigma = c(s0, p0) means P(sigma > s0) = p0
spde_model <- inla.spde2.matern(
  mesh = mesh,
  alpha = 2,
  prior.range = c(15, 0.50),
  prior.sigma = c(1.0, 0.01)
)

# Create the latent field index for stacking.
field_index <- inla.spde.make.index(name = "spatial_field", n.spde = spde_model$n.spde)

# 5) Construct estimation and prediction stacks -----------------------------
# A matrix for observation locations
A_est <- inla.spde.make.A(mesh = mesh, loc = coords_matrix)

# Estimation stack (observed counts)
stack_est <- inla.stack(
  data = list(cases = disease_df$cases, expected = disease_df$expected),
  A = list(1, A_est),
  effects = list(
    data.frame(intercept = 1, deprivation = disease_df$deprivation),
    spatial_field = field_index
  ),
  tag = "estimation"
)

# Prediction grid across the study region
prediction_grid <- expand.grid(
  x = seq(min(coords_matrix[, 1]), max(coords_matrix[, 1]), length.out = 120),
  y = seq(min(coords_matrix[, 2]), max(coords_matrix[, 2]), length.out = 120)
)
prediction_grid$deprivation <- 0
prediction_grid$expected <- 100
A_pred <- inla.spde.make.A(mesh = mesh, loc = as.matrix(prediction_grid[, c("x", "y")]))

stack_pred <- inla.stack(
  data = list(cases = NA, expected = prediction_grid$expected),
  A = list(1, A_pred),
  effects = list(
    data.frame(intercept = 1, deprivation = prediction_grid$deprivation),
    spatial_field = field_index
  ),
  tag = "prediction"
)

# Combine stacks for one INLA call
stack_full <- inla.stack(stack_est, stack_pred)

# 6) Define the model formula -----------------------------------------------
# Offset(log(expected)) converts the count model into disease risk mapping.
# The `f()` term introduces the latent spatial SPDE field.
model_formula <- cases ~ 0 + intercept + deprivation +
  f(spatial_field, model = spde_model)

# 7) Fit the INLA model -----------------------------------------------------
# family = "poisson" is common for disease counts with expected counts E.
# control.compute requests model assessment measures and posterior config.
fit_spde <- inla(
  formula = model_formula,
  family = "poisson",
  data = inla.stack.data(stack_full),
  E = inla.stack.data(stack_full)$expected,
  control.predictor = list(
    A = inla.stack.A(stack_full),
    compute = TRUE,
    link = 1
  ),
  control.compute = list(
    dic = TRUE,
    waic = TRUE,
    cpo = TRUE,
    config = TRUE
  ),
  control.inla = list(strategy = "adaptive", int.strategy = "eb")
)

# 8) Extract posterior summaries -------------------------------------------
summary_fixed <- fit_spde$summary.fixed
summary_hyperpar <- fit_spde$summary.hyperpar

print(summary_fixed)
print(summary_hyperpar)

# Indices for observed and prediction components
idx_est <- inla.stack.index(stack_full, tag = "estimation")$data
idx_pred <- inla.stack.index(stack_full, tag = "prediction")$data

# Relative risk posterior summaries on the link scale and risk scale
prediction_grid$eta_mean <- fit_spde$summary.linear.predictor[idx_pred, "mean"]
prediction_grid$eta_sd <- fit_spde$summary.linear.predictor[idx_pred, "sd"]
prediction_grid$risk_mean <- exp(prediction_grid$eta_mean)
prediction_grid$risk_lower <- exp(prediction_grid$eta_mean - 1.96 * prediction_grid$eta_sd)
prediction_grid$risk_upper <- exp(prediction_grid$eta_mean + 1.96 * prediction_grid$eta_sd)

# 9) Visualise the predicted risk surface -----------------------------------
# Project the latent field onto a regular grid using the INLA mesh projector.
projector <- inla.mesh.projector(
  mesh,
  xlim = range(prediction_grid$x),
  ylim = range(prediction_grid$y),
  dims = c(200, 200)
)

spatial_mean_projected <- inla.mesh.project(
  projector,
  field = fit_spde$summary.random$spatial_field$mean
)

plot_df <- data.frame(
  expand.grid(x = projector$x, y = projector$y),
  spatial_mean = as.vector(spatial_mean_projected)
)

risk_plot <- ggplot(prediction_grid, aes(x = x, y = y, fill = risk_mean)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c(option = "viridis", name = "Risk") +
  labs(
    title = "Posterior mean disease risk",
    x = "x coordinate",
    y = "y coordinate"
  ) +
  theme_minimal(base_size = 12)

uncertainty_plot <- ggplot(prediction_grid, aes(x = x, y = y, fill = risk_upper - risk_lower)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c(option = "cividis", name = "95% CI width") +
  labs(
    title = "Posterior uncertainty in disease risk",
    x = "x coordinate",
    y = "y coordinate"
  ) +
  theme_minimal(base_size = 12)

print(risk_plot)
print(uncertainty_plot)

# 10) Compare alternative models --------------------------------------------
# Example baseline model without the spatial random effect.
fit_nonspatial <- inla(
  formula = cases ~ 0 + intercept + deprivation,
  family = "poisson",
  data = as.data.frame(disease_df),
  E = disease_df$expected,
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

comparison_table <- data.frame(
  model = c("Spatial SPDE", "Non-spatial"),
  DIC = c(fit_spde$dic$dic, fit_nonspatial$dic$dic),
  WAIC = c(fit_spde$waic$waic, fit_nonspatial$waic$waic),
  mean_neg_log_CPO = c(
    mean(-log(fit_spde$cpo$cpo), na.rm = TRUE),
    mean(-log(fit_nonspatial$cpo$cpo), na.rm = TRUE)
  )
)

print(comparison_table)

# Suggested next steps -------------------------------------------------------
# - Replace the simulated data block with observed disease counts and covariates.
# - Tune mesh density, priors, and boundary extension for the study region.
# - Add covariates such as age structure, pollution, or healthcare access.
# - Inspect residuals, PIT/CPO diagnostics, and sensitivity to prior choices.
# - Export the risk surfaces with ggsave() or sf/stars for GIS workflows.
