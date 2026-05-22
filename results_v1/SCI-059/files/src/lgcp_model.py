import json
import random
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from scipy.spatial.distance import cdist, pdist
from scipy.special import gamma, kv
from scipy.stats import gaussian_kde


SEED = 2025
WINDOW = (0.0, 1.0, 0.0, 1.0)
GRID_SIZE = 18
MATERN_NU = 1.5
MATERN_RHO = 0.18
MATERN_SIGMA = 0.9
LOG_BASE_INTENSITY = np.log(120.0)
JITTER = 1e-6


def matern_covariance(coords, sigma=1.0, rho=0.2, nu=1.5):
    distances = cdist(coords, coords)
    scaled = np.sqrt(2.0 * nu) * distances / rho
    covariance = np.empty_like(distances)
    zero_mask = distances == 0.0
    covariance[zero_mask] = sigma ** 2
    nonzero = ~zero_mask
    scaled_nonzero = scaled[nonzero]
    factor = (2.0 ** (1.0 - nu)) / gamma(nu)
    covariance[nonzero] = (
        sigma ** 2
        * factor
        * (scaled_nonzero ** nu)
        * kv(nu, scaled_nonzero)
    )
    return covariance


def build_grid(grid_size=GRID_SIZE, window=WINDOW):
    xmin, xmax, ymin, ymax = window
    x_edges = np.linspace(xmin, xmax, grid_size + 1)
    y_edges = np.linspace(ymin, ymax, grid_size + 1)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    xx, yy = np.meshgrid(x_centers, y_centers, indexing="xy")
    coords = np.column_stack([xx.ravel(), yy.ravel()])
    cell_area = (x_edges[1] - x_edges[0]) * (y_edges[1] - y_edges[0])
    return x_edges, y_edges, x_centers, y_centers, xx, yy, coords, cell_area


def simulate_lgcp(rng, coords, shape, cell_area):
    covariance = matern_covariance(
        coords,
        sigma=MATERN_SIGMA,
        rho=MATERN_RHO,
        nu=MATERN_NU,
    )
    covariance = covariance + JITTER * np.eye(covariance.shape[0])
    latent = rng.multivariate_normal(mean=np.zeros(coords.shape[0]), cov=covariance)
    latent_field = latent.reshape(shape)
    true_intensity = np.exp(LOG_BASE_INTENSITY + latent_field)
    counts = rng.poisson(true_intensity * cell_area)
    return latent_field, true_intensity, counts, covariance


def sample_points_from_counts(rng, counts, x_edges, y_edges):
    points = []
    ny, nx = counts.shape
    for iy in range(ny):
        for ix in range(nx):
            count = int(counts[iy, ix])
            if count <= 0:
                continue
            xs = rng.uniform(x_edges[ix], x_edges[ix + 1], size=count)
            ys = rng.uniform(y_edges[iy], y_edges[iy + 1], size=count)
            points.append(np.column_stack([xs, ys]))
    if not points:
        return np.empty((0, 2))
    return np.vstack(points)


def kde_intensity(points, xx, yy, total_points):
    if len(points) < 2:
        return np.full_like(xx, fill_value=total_points, dtype=float)
    kde = gaussian_kde(points.T, bw_method="scott")
    density = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
    return density * total_points


def lgcp_objective(params, counts_flat, cell_area, prior_factor):
    intercept = params[0]
    latent = params[1:]
    eta = intercept + latent
    mean_counts = cell_area * np.exp(eta)
    prior_term = cho_solve(prior_factor, latent)
    objective = np.sum(mean_counts - counts_flat * eta) + 0.5 * latent.dot(prior_term)
    gradient = np.empty_like(params)
    poisson_grad = mean_counts - counts_flat
    gradient[0] = np.sum(poisson_grad)
    gradient[1:] = poisson_grad + prior_term
    return objective, gradient


def fit_lgcp(counts, covariance, cell_area):
    counts_flat = counts.ravel().astype(float)
    prior_factor = cho_factor(covariance, lower=True, check_finite=False)
    total_points = counts_flat.sum()
    intercept0 = np.log((total_points + 1e-6) / (counts_flat.size * cell_area))
    x0 = np.concatenate([[intercept0], np.zeros_like(counts_flat)])

    result = minimize(
        lambda params: lgcp_objective(params, counts_flat, cell_area, prior_factor),
        x0,
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": 250, "maxfun": 400},
    )

    intercept = float(result.x[0])
    latent = result.x[1:].reshape(counts.shape)
    fitted_intensity = np.exp(intercept + latent)
    return intercept, latent, fitted_intensity, result


def rmse(predicted, observed):
    return float(np.sqrt(np.mean((predicted - observed) ** 2)))


def compute_k_function(points, area=1.0, max_r=0.35, n_steps=30):
    radii = np.linspace(max_r / n_steps, max_r, n_steps)
    if len(points) < 2:
        return radii, np.zeros_like(radii)
    dists = np.sort(pdist(points))
    pair_counts = np.searchsorted(dists, radii, side="right")
    k_values = area * 2.0 * pair_counts / (len(points) * (len(points) - 1))
    return radii, k_values


def compute_pair_correlation(points, area=1.0, max_r=0.35, n_bins=24):
    bins = np.linspace(0.0, max_r, n_bins + 1)
    radii = 0.5 * (bins[1:] + bins[:-1])
    if len(points) < 2:
        return radii, np.zeros_like(radii)
    dists = pdist(points)
    counts, _ = np.histogram(dists, bins=bins)
    annulus_areas = np.pi * (bins[1:] ** 2 - bins[:-1] ** 2)
    g_values = area * 2.0 * counts / (len(points) * (len(points) - 1) * annulus_areas)
    return radii, g_values


def plot_intensity_surface(path, x_edges, y_edges, true_intensity, fitted_intensity):
    vmax = max(true_intensity.max(), fitted_intensity.max())
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    images = []
    for ax, surface, title in zip(
        axes,
        [true_intensity, fitted_intensity],
        ["True LGCP intensity", "Estimated LGCP intensity"],
    ):
        image = ax.imshow(
            surface,
            origin="lower",
            extent=(x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]),
            cmap="viridis",
            aspect="auto",
            vmin=0.0,
            vmax=vmax,
        )
        images.append(image)
        ax.set_title(title)
        ax.set_xlabel("X coordinate")
        ax.set_ylabel("Y coordinate")
    cbar = fig.colorbar(images[-1], ax=axes, shrink=0.9)
    cbar.set_label("Intensity")
    fig.savefig(path, dpi=300)
    plt.close(fig)


def plot_point_pattern(path, x_edges, y_edges, true_intensity, points):
    fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)
    image = ax.imshow(
        true_intensity,
        origin="lower",
        extent=(x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]),
        cmap="cividis",
        aspect="auto",
    )
    if len(points) > 0:
        ax.scatter(points[:, 0], points[:, 1], s=18, c="white", edgecolors="black", linewidths=0.3)
    ax.set_title("Simulated disease case locations")
    ax.set_xlabel("X coordinate")
    ax.set_ylabel("Y coordinate")
    cbar = fig.colorbar(image, ax=ax, shrink=0.9)
    cbar.set_label("True intensity")
    fig.savefig(path, dpi=300)
    plt.close(fig)


def plot_kde_comparison(path, x_edges, y_edges, true_intensity, fitted_intensity, kde_surface, lgcp_rmse, kde_rmse):
    vmax = max(true_intensity.max(), fitted_intensity.max(), kde_surface.max())
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    titles = [
        "True intensity",
        f"LGCP estimate (RMSE={lgcp_rmse:.2f})",
        f"KDE estimate (RMSE={kde_rmse:.2f})",
    ]
    images = []
    for ax, surface, title in zip(axes, [true_intensity, fitted_intensity, kde_surface], titles):
        image = ax.imshow(
            surface,
            origin="lower",
            extent=(x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]),
            cmap="viridis",
            aspect="auto",
            vmin=0.0,
            vmax=vmax,
        )
        images.append(image)
        ax.set_title(title)
        ax.set_xlabel("X coordinate")
        ax.set_ylabel("Y coordinate")
    cbar = fig.colorbar(images[-1], ax=axes, shrink=0.9)
    cbar.set_label("Intensity")
    fig.savefig(path, dpi=300)
    plt.close(fig)


def save_results(path, results):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)


def main():
    random.seed(SEED)
    np.random.seed(SEED)
    rng = np.random.default_rng(SEED)

    base_dir = Path(__file__).resolve().parents[1]
    figures_dir = base_dir / "figures"
    results_dir = base_dir / "results"
    figures_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    x_edges, y_edges, x_centers, y_centers, xx, yy, coords, cell_area = build_grid()
    grid_shape = (len(y_centers), len(x_centers))

    latent_field, true_intensity, counts, covariance = simulate_lgcp(rng, coords, grid_shape, cell_area)
    points = sample_points_from_counts(rng, counts, x_edges, y_edges)

    intercept, latent_estimate, fitted_intensity, optimization = fit_lgcp(counts, covariance, cell_area)
    kde_surface = kde_intensity(points, xx, yy, total_points=len(points))

    lgcp_rmse = rmse(fitted_intensity, true_intensity)
    kde_rmse = rmse(kde_surface, true_intensity)
    lgcp_mae = float(np.mean(np.abs(fitted_intensity - true_intensity)))
    kde_mae = float(np.mean(np.abs(kde_surface - true_intensity)))

    radii_k, k_values = compute_k_function(points)
    radii_g, g_values = compute_pair_correlation(points)
    poisson_k = np.pi * radii_k ** 2

    plot_intensity_surface(figures_dir / "lgcp_intensity_surface.png", x_edges, y_edges, true_intensity, fitted_intensity)
    plot_point_pattern(figures_dir / "lgcp_point_pattern.png", x_edges, y_edges, true_intensity, points)
    plot_kde_comparison(
        figures_dir / "lgcp_kde_comparison.png",
        x_edges,
        y_edges,
        true_intensity,
        fitted_intensity,
        kde_surface,
        lgcp_rmse,
        kde_rmse,
    )

    results = {
        "seed": SEED,
        "grid_size": GRID_SIZE,
        "window": {"xmin": WINDOW[0], "xmax": WINDOW[1], "ymin": WINDOW[2], "ymax": WINDOW[3]},
        "matern_parameters": {
            "nu": MATERN_NU,
            "rho": MATERN_RHO,
            "sigma": MATERN_SIGMA,
            "log_base_intensity": float(LOG_BASE_INTENSITY),
        },
        "simulation_summary": {
            "total_cases": int(len(points)),
            "total_cells": int(counts.size),
            "mean_true_intensity": float(np.mean(true_intensity)),
            "max_true_intensity": float(np.max(true_intensity)),
            "mean_observed_count_per_cell": float(np.mean(counts)),
        },
        "lgcp_fit": {
            "optimization_success": bool(optimization.success),
            "optimizer_message": optimization.message,
            "iterations": int(getattr(optimization, "nit", 0)),
            "estimated_intercept": intercept,
            "rmse_vs_true_intensity": lgcp_rmse,
            "mae_vs_true_intensity": lgcp_mae,
        },
        "kde_fit": {
            "rmse_vs_true_intensity": kde_rmse,
            "mae_vs_true_intensity": kde_mae,
        },
        "summary_statistics": {
            "k_function_radii": radii_k.tolist(),
            "k_function_values": k_values.tolist(),
            "poisson_k_reference": poisson_k.tolist(),
            "pair_correlation_radii": radii_g.tolist(),
            "pair_correlation_values": g_values.tolist(),
        },
        "output_files": {
            "figures": [
                str(figures_dir / "lgcp_intensity_surface.png"),
                str(figures_dir / "lgcp_point_pattern.png"),
                str(figures_dir / "lgcp_kde_comparison.png"),
            ],
            "results": str(results_dir / "lgcp_results.json"),
        },
    }
    save_results(results_dir / "lgcp_results.json", results)

    print("LGCP disease risk analysis complete")
    print(f"Simulated cases: {len(points)}")
    print(f"Grid size: {GRID_SIZE} x {GRID_SIZE}")
    print(f"LGCP optimization success: {optimization.success}")
    print(f"LGCP RMSE vs. true intensity: {lgcp_rmse:.3f}")
    print(f"KDE RMSE vs. true intensity: {kde_rmse:.3f}")
    print(f"Results written to: {results_dir / 'lgcp_results.json'}")


if __name__ == "__main__":
    main()
