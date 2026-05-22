#!/usr/bin/env python3
"""Demonstration of INLA/SPDE ideas for Bayesian spatial disease mapping.

The script:
1. Generates synthetic disease count data on an irregular Delaunay mesh.
2. Builds finite-element mass (C) and stiffness (G) matrices.
3. Constructs the SPDE precision matrix
   Q = tau^2 (kappa^4 C + 2 kappa^2 G + G C^{-1} G)
   for a Matérn field with alpha = 2.
4. Fits a latent Gaussian Poisson model with a simplified Laplace approximation.
5. Predicts disease risk on a fine grid and saves posterior mean / uncertainty maps.
"""

from __future__ import annotations

import json
import random
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
import scipy.sparse.linalg as spla
from scipy.spatial import Delaunay

SEED = 2025


def ensure_output_dirs(root: Path) -> None:
    (root / "figures").mkdir(exist_ok=True)
    (root / "results").mkdir(exist_ok=True)


def generate_irregular_mesh(rng: np.random.Generator) -> tuple[np.ndarray, Delaunay]:
    boundary_axis = np.linspace(0.0, 1.0, 6)
    boundary_points = []
    for value in boundary_axis:
        boundary_points.extend(
            [(value, 0.0), (value, 1.0), (0.0, value), (1.0, value)]
        )
    boundary = np.unique(np.asarray(boundary_points, dtype=float), axis=0)
    interior = rng.uniform(0.03, 0.97, size=(55, 2))
    points = np.vstack([boundary, interior])
    triangulation = Delaunay(points)
    return points, triangulation


def local_fem_matrices(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    x1, y1 = coords[0]
    x2, y2 = coords[1]
    x3, y3 = coords[2]
    area = 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
    if area <= 1e-12:
        raise ValueError("Degenerate triangle encountered in mesh.")

    mass = (area / 12.0) * np.array(
        [[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]], dtype=float
    )

    b = np.array([y2 - y3, y3 - y1, y1 - y2], dtype=float)
    c = np.array([x3 - x2, x1 - x3, x2 - x1], dtype=float)
    stiffness = (np.outer(b, b) + np.outer(c, c)) / (4.0 * area)
    return mass, stiffness, area


def assemble_fem(points: np.ndarray, simplices: np.ndarray) -> tuple[sparse.csc_matrix, sparse.csc_matrix]:
    n_vertices = len(points)
    rows_c, cols_c, data_c = [], [], []
    rows_g, cols_g, data_g = [], [], []

    for simplex in simplices:
        local_c, local_g, _ = local_fem_matrices(points[simplex])
        for i_local, i_global in enumerate(simplex):
            for j_local, j_global in enumerate(simplex):
                rows_c.append(i_global)
                cols_c.append(j_global)
                data_c.append(local_c[i_local, j_local])
                rows_g.append(i_global)
                cols_g.append(j_global)
                data_g.append(local_g[i_local, j_local])

    mass = sparse.coo_matrix((data_c, (rows_c, cols_c)), shape=(n_vertices, n_vertices)).tocsc()
    stiffness = sparse.coo_matrix((data_g, (rows_g, cols_g)), shape=(n_vertices, n_vertices)).tocsc()
    return mass, stiffness


def build_precision(
    mass: sparse.csc_matrix,
    stiffness: sparse.csc_matrix,
    kappa: float,
    tau: float,
    regularizer: float = 1e-6,
) -> sparse.csc_matrix:
    n_vertices = mass.shape[0]
    mass_reg = mass + regularizer * sparse.eye(n_vertices, format="csc")
    c_inv_g = spla.spsolve(mass_reg, stiffness.toarray())
    g_c_inv_g = np.asarray(stiffness @ c_inv_g)
    q_dense = tau**2 * (
        (kappa**4) * mass.toarray()
        + 2.0 * (kappa**2) * stiffness.toarray()
        + g_c_inv_g
    )
    q_dense = 0.5 * (q_dense + q_dense.T) + regularizer * np.eye(n_vertices)
    return sparse.csc_matrix(q_dense)


def sample_from_precision(
    precision: sparse.csc_matrix, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray]:
    covariance = np.linalg.inv(precision.toarray())
    sample = rng.multivariate_normal(np.zeros(precision.shape[0]), covariance)
    return sample, covariance


def generate_synthetic_counts(
    n_vertices: int, precision: sparse.csc_matrix, rng: np.random.Generator
) -> dict[str, np.ndarray | float]:
    latent_field, _ = sample_from_precision(precision, rng)
    latent_field = latent_field - latent_field.mean()
    beta_true = -2.35
    exposure = rng.integers(60, 160, size=n_vertices).astype(float)
    eta = beta_true + latent_field
    risk = np.exp(eta)
    counts = rng.poisson(exposure * risk)
    return {
        "latent_field": latent_field,
        "beta_true": beta_true,
        "exposure": exposure,
        "eta": eta,
        "risk": risk,
        "counts": counts,
    }


def posterior_terms(
    theta: np.ndarray,
    counts: np.ndarray,
    exposure: np.ndarray,
    q_dense: np.ndarray,
    sigma_beta: float,
) -> tuple[float, np.ndarray, np.ndarray]:
    beta = theta[0]
    field = theta[1:]
    eta = beta + field
    mu = exposure * np.exp(eta)

    nlp = np.sum(mu - counts * eta)
    nlp += 0.5 * field @ (q_dense @ field)
    nlp += 0.5 * (beta / sigma_beta) ** 2

    residual = mu - counts
    grad = np.empty(theta.size)
    grad[0] = residual.sum() + beta / sigma_beta**2
    grad[1:] = residual + q_dense @ field

    hessian = np.zeros((theta.size, theta.size), dtype=float)
    hessian[0, 0] = mu.sum() + 1.0 / sigma_beta**2
    hessian[0, 1:] = mu
    hessian[1:, 0] = mu
    hessian[1:, 1:] = q_dense + np.diag(mu)
    return float(nlp), grad, hessian


def fit_laplace(
    precision: sparse.csc_matrix,
    counts: np.ndarray,
    exposure: np.ndarray,
    sigma_beta: float = 5.0,
    init_theta: np.ndarray | None = None,
    max_iter: int = 30,
    tol: float = 1e-6,
) -> dict[str, np.ndarray | float | int | list[dict[str, float]]]:
    n_vertices = counts.size
    q_dense = precision.toarray()
    if init_theta is None:
        theta = np.zeros(n_vertices + 1, dtype=float)
        theta[0] = np.log((counts.sum() + 0.5) / (exposure.sum() + 0.5))
    else:
        theta = init_theta.astype(float).copy()

    history: list[dict[str, float]] = []
    for iteration in range(1, max_iter + 1):
        nlp, grad, hessian = posterior_terms(theta, counts, exposure, q_dense, sigma_beta)
        grad_norm = float(np.linalg.norm(grad))
        history.append(
            {
                "iteration": float(iteration),
                "negative_log_posterior": float(nlp),
                "gradient_norm": grad_norm,
            }
        )
        if grad_norm < tol:
            break

        try:
            step = np.linalg.solve(hessian, -grad)
        except np.linalg.LinAlgError:
            step = np.linalg.solve(hessian + 1e-6 * np.eye(hessian.shape[0]), -grad)

        step_scale = 1.0
        accepted = False
        while step_scale >= 1e-4:
            candidate = theta + step_scale * step
            candidate_nlp, _, _ = posterior_terms(
                candidate, counts, exposure, q_dense, sigma_beta
            )
            if np.isfinite(candidate_nlp) and candidate_nlp < nlp:
                theta = candidate
                accepted = True
                break
            step_scale *= 0.5

        if not accepted:
            theta = theta + 0.1 * step

        if np.linalg.norm(step_scale * step, ord=np.inf) < tol:
            break

    nlp, grad, hessian = posterior_terms(theta, counts, exposure, q_dense, sigma_beta)
    covariance = np.linalg.inv(hessian)
    return {
        "theta_mode": theta,
        "posterior_cov": covariance,
        "hessian": hessian,
        "negative_log_posterior": float(nlp),
        "gradient_norm": float(np.linalg.norm(grad)),
        "iterations": int(len(history)),
        "history": history,
    }


def laplace_log_evidence(fit: dict[str, np.ndarray | float | int]) -> float:
    hessian = np.asarray(fit["hessian"])
    sign, logdet = np.linalg.slogdet(hessian)
    if sign <= 0:
        return -np.inf
    n_dim = hessian.shape[0]
    return float(
        -float(fit["negative_log_posterior"])
        + 0.5 * n_dim * np.log(2.0 * np.pi)
        - 0.5 * logdet
    )


def select_hyperparameters(
    mass: sparse.csc_matrix,
    stiffness: sparse.csc_matrix,
    counts: np.ndarray,
    exposure: np.ndarray,
) -> tuple[dict[str, object], list[dict[str, float]]]:
    kappa_grid = [3.5, 5.0, 7.0]
    tau_grid = [0.7, 0.9, 1.3]
    evaluations: list[dict[str, float]] = []
    best_score = -np.inf
    best_result: dict[str, object] | None = None
    warm_start: np.ndarray | None = None

    for kappa in kappa_grid:
        for tau in tau_grid:
            precision = build_precision(mass, stiffness, kappa=kappa, tau=tau)
            fit = fit_laplace(precision, counts, exposure, init_theta=warm_start)
            score = laplace_log_evidence(fit)
            warm_start = np.asarray(fit["theta_mode"])
            evaluations.append(
                {
                    "kappa": float(kappa),
                    "tau": float(tau),
                    "laplace_log_evidence": float(score),
                    "iterations": float(fit["iterations"]),
                }
            )
            if score > best_score:
                best_score = score
                best_result = {
                    "kappa": float(kappa),
                    "tau": float(tau),
                    "precision": precision,
                    "fit": fit,
                    "laplace_log_evidence": float(score),
                }

    if best_result is None:
        raise RuntimeError("Hyperparameter search failed.")
    return best_result, evaluations


def predict_on_grid(
    triangulation: Delaunay,
    theta_mode: np.ndarray,
    posterior_cov: np.ndarray,
    grid_size: int = 140,
) -> dict[str, np.ndarray]:
    x_coords = np.linspace(0.0, 1.0, grid_size)
    y_coords = np.linspace(0.0, 1.0, grid_size)
    xx, yy = np.meshgrid(x_coords, y_coords)
    grid_points = np.column_stack([xx.ravel(), yy.ravel()])

    simplex_ids = triangulation.find_simplex(grid_points)
    posterior_mean = np.full(grid_points.shape[0], np.nan)
    lower = np.full(grid_points.shape[0], np.nan)
    upper = np.full(grid_points.shape[0], np.nan)
    ci_width = np.full(grid_points.shape[0], np.nan)

    beta_mean = theta_mode[0]
    field_mean = theta_mode[1:]
    var_beta = posterior_cov[0, 0]
    cov_beta_field = posterior_cov[0, 1:]
    cov_field = posterior_cov[1:, 1:]

    for idx, simplex in enumerate(simplex_ids):
        if simplex < 0:
            continue

        transform = triangulation.transform[simplex]
        barycentric = transform[:2].dot(grid_points[idx] - transform[2])
        barycentric = np.append(barycentric, 1.0 - barycentric.sum())
        vertices = triangulation.simplices[simplex]

        eta_mean = beta_mean + barycentric @ field_mean[vertices]
        local_cov = cov_field[np.ix_(vertices, vertices)]
        eta_var = float(var_beta + barycentric @ local_cov @ barycentric)
        eta_var += float(2.0 * (barycentric @ cov_beta_field[vertices]))
        eta_var = max(eta_var, 0.0)
        eta_sd = np.sqrt(eta_var)

        posterior_mean[idx] = np.exp(eta_mean + 0.5 * eta_var)
        lower[idx] = np.exp(eta_mean - 1.96 * eta_sd)
        upper[idx] = np.exp(eta_mean + 1.96 * eta_sd)
        ci_width[idx] = upper[idx] - lower[idx]

    return {
        "x_coords": x_coords,
        "y_coords": y_coords,
        "posterior_mean": posterior_mean.reshape(grid_size, grid_size),
        "lower": lower.reshape(grid_size, grid_size),
        "upper": upper.reshape(grid_size, grid_size),
        "ci_width": ci_width.reshape(grid_size, grid_size),
    }


def plot_mesh(
    points: np.ndarray,
    triangulation: Delaunay,
    counts: np.ndarray,
    exposure: np.ndarray,
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.triplot(points[:, 0], points[:, 1], triangulation.simplices, color="0.60", linewidth=0.7)
    observed_rate = counts / exposure
    scatter = ax.scatter(
        points[:, 0],
        points[:, 1],
        c=observed_rate,
        cmap="viridis",
        s=36,
        edgecolor="black",
        linewidth=0.25,
    )
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("Observed rate")
    ax.set_title("Irregular SPDE mesh and observed disease rates")
    ax.set_xlabel("x coordinate")
    ax.set_ylabel("y coordinate")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_grid_map(
    field: np.ndarray,
    x_coords: np.ndarray,
    y_coords: np.ndarray,
    output_path: Path,
    title: str,
    colorbar_label: str,
    cmap: str,
) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    masked = np.ma.masked_invalid(field)
    image = ax.imshow(
        masked,
        origin="lower",
        extent=(x_coords.min(), x_coords.max(), y_coords.min(), y_coords.max()),
        cmap=cmap,
        interpolation="nearest",
        aspect="equal",
    )
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label(colorbar_label)
    ax.set_title(title)
    ax.set_xlabel("x coordinate")
    ax.set_ylabel("y coordinate")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_results(root: Path, payload: dict[str, object]) -> Path:
    output_path = root / "results" / "bayesian_spatial_results.json"
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return output_path


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    ensure_output_dirs(root)
    random.seed(SEED)
    rng = np.random.default_rng(SEED)

    points, triangulation = generate_irregular_mesh(rng)
    mass, stiffness = assemble_fem(points, triangulation.simplices)

    true_kappa = 5.0
    true_tau = 0.9
    q_true = build_precision(mass, stiffness, kappa=true_kappa, tau=true_tau)
    synthetic = generate_synthetic_counts(len(points), q_true, rng)
    counts = np.asarray(synthetic["counts"])
    exposure = np.asarray(synthetic["exposure"])

    best_model, hyperparameter_grid = select_hyperparameters(mass, stiffness, counts, exposure)
    best_fit = best_model["fit"]
    theta_mode = np.asarray(best_fit["theta_mode"])
    posterior_cov = np.asarray(best_fit["posterior_cov"])

    prediction = predict_on_grid(triangulation, theta_mode, posterior_cov)

    mesh_path = root / "figures" / "spde_mesh.png"
    mean_path = root / "figures" / "spde_posterior_mean.png"
    uncertainty_path = root / "figures" / "spde_posterior_uncertainty.png"

    plot_mesh(points, triangulation, counts, exposure, mesh_path)
    plot_grid_map(
        prediction["posterior_mean"],
        prediction["x_coords"],
        prediction["y_coords"],
        mean_path,
        title="Posterior mean disease risk",
        colorbar_label="Posterior mean risk",
        cmap="viridis",
    )
    plot_grid_map(
        prediction["ci_width"],
        prediction["x_coords"],
        prediction["y_coords"],
        uncertainty_path,
        title="Posterior uncertainty (95% credible interval width)",
        colorbar_label="95% CI width",
        cmap="cividis",
    )

    vertex_sd = np.sqrt(np.clip(np.diag(posterior_cov)[1:], 0.0, None))
    posterior_mean_grid = prediction["posterior_mean"]
    ci_width_grid = prediction["ci_width"]
    valid_mean = posterior_mean_grid[np.isfinite(posterior_mean_grid)]
    valid_width = ci_width_grid[np.isfinite(ci_width_grid)]

    results = {
        "seed": SEED,
        "mesh": {
            "n_vertices": int(points.shape[0]),
            "n_triangles": int(triangulation.simplices.shape[0]),
        },
        "spde": {
            "alpha": 2,
            "precision_formula": "Q = tau^2 (kappa^4 C + 2 kappa^2 G + G C^{-1} G)",
            "true_kappa": float(true_kappa),
            "true_tau": float(true_tau),
            "selected_kappa": float(best_model["kappa"]),
            "selected_tau": float(best_model["tau"]),
            "laplace_log_evidence": float(best_model["laplace_log_evidence"]),
        },
        "data_summary": {
            "total_counts": int(counts.sum()),
            "total_exposure": float(exposure.sum()),
            "mean_observed_rate": float(np.mean(counts / exposure)),
            "max_count": int(counts.max()),
        },
        "posterior_summary": {
            "beta_mode": float(theta_mode[0]),
            "beta_sd": float(np.sqrt(max(posterior_cov[0, 0], 0.0))),
            "mean_latent_sd": float(vertex_sd.mean()),
            "negative_log_posterior": float(best_fit["negative_log_posterior"]),
            "gradient_norm": float(best_fit["gradient_norm"]),
            "iterations": int(best_fit["iterations"]),
        },
        "prediction_summary": {
            "mean_risk": float(valid_mean.mean()),
            "min_risk": float(valid_mean.min()),
            "max_risk": float(valid_mean.max()),
            "mean_ci_width": float(valid_width.mean()),
            "max_ci_width": float(valid_width.max()),
        },
        "grid_search": hyperparameter_grid,
        "artifacts": {
            "mesh_figure": str(mesh_path.relative_to(root)),
            "posterior_mean_figure": str(mean_path.relative_to(root)),
            "posterior_uncertainty_figure": str(uncertainty_path.relative_to(root)),
            "results_json": "results/bayesian_spatial_results.json",
        },
    }
    save_results(root, results)

    print("Bayesian spatial disease mapping demonstration complete.")
    print(f"Mesh vertices: {results['mesh']['n_vertices']}, triangles: {results['mesh']['n_triangles']}")
    print(
        "Selected hyperparameters: "
        f"kappa={results['spde']['selected_kappa']:.2f}, tau={results['spde']['selected_tau']:.2f}"
    )
    print(
        "Posterior intercept (mode ± SD): "
        f"{results['posterior_summary']['beta_mode']:.3f} ± {results['posterior_summary']['beta_sd']:.3f}"
    )
    print(
        "Predicted risk range: "
        f"{results['prediction_summary']['min_risk']:.4f} to {results['prediction_summary']['max_risk']:.4f}"
    )
    print("Saved figures: figures/spde_mesh.png, figures/spde_posterior_mean.png, figures/spde_posterior_uncertainty.png")
    print("Saved results: results/bayesian_spatial_results.json")


if __name__ == "__main__":
    main()
