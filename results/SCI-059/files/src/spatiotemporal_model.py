import json
import random
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import SplineTransformer, StandardScaler


SEED = 123
N_MONTHS = 60
FUTURE_MONTHS = 12


def ensure_directories(root: Path) -> None:
    for name in ("figures", "results", "data", "logs"):
        (root / name).mkdir(parents=True, exist_ok=True)


def thin_plate_kernel(distance: np.ndarray) -> np.ndarray:
    distance = np.maximum(distance, 1e-8)
    return (distance ** 2) * np.log(distance)


def simulate_spatiotemporal_data(seed: int = SEED) -> dict:
    rng = np.random.default_rng(seed)
    random.seed(seed)

    n_locations = 90
    coords = rng.uniform(0.0, 1.0, size=(n_locations, 2))
    population = rng.integers(1800, 4200, size=n_locations)
    times = np.arange(N_MONTHS)

    def latent_log_rate(query_coords: np.ndarray, query_times: np.ndarray) -> np.ndarray:
        x = query_coords[:, 0][:, None]
        y = query_coords[:, 1][:, None]
        tt = query_times[None, :]
        bump_north = np.exp(-((x - 0.25) ** 2 + (y - 0.75) ** 2) / 0.03)
        bump_south = np.exp(-((x - 0.72) ** 2 + (y - 0.28) ** 2) / 0.05)
        ridge = 0.4 * np.sin(np.pi * x) * np.cos(np.pi * y)
        temporal = -3.4 + 0.012 * tt + 0.23 * np.sin(2 * np.pi * tt / 12) + 0.08 * np.cos(2 * np.pi * tt / 24)
        interaction = 0.30 * bump_north * np.sin(2 * np.pi * (tt - 4) / 18) - 0.18 * bump_south * np.cos(2 * np.pi * tt / 30)
        return temporal + 0.85 * bump_north - 0.60 * bump_south + ridge + interaction

    true_log_rate = latent_log_rate(coords, times)
    noise = rng.normal(scale=0.05, size=true_log_rate.shape)
    noisy_log_rate = true_log_rate + noise
    expected_cases = population[:, None] * np.exp(noisy_log_rate)
    cases = rng.poisson(expected_cases)

    location_id = np.repeat(np.arange(n_locations), N_MONTHS)
    time_id = np.tile(times, n_locations)
    response = np.log((cases.reshape(-1) + 0.5) / np.repeat(population, N_MONTHS))

    future_times = np.arange(N_MONTHS, N_MONTHS + FUTURE_MONTHS)
    future_true_log_rate = latent_log_rate(coords, future_times)

    return {
        "coords": coords,
        "population": population,
        "times": times,
        "future_times": future_times,
        "cases": cases,
        "response": response,
        "location_id": location_id,
        "time_id": time_id,
        "true_log_rate": true_log_rate,
        "future_true_log_rate": future_true_log_rate,
    }


def choose_knots(coords: np.ndarray, n_knots: int, seed: int) -> np.ndarray:
    model = KMeans(n_clusters=n_knots, n_init=20, random_state=seed)
    model.fit(coords)
    return model.cluster_centers_


def build_design_components(coords: np.ndarray, times: np.ndarray, knot_count: int, seed: int, spline_transformer=None) -> dict:
    knots = choose_knots(coords, knot_count, seed)
    radial = thin_plate_kernel(cdist(coords, knots))
    spatial_main = np.column_stack([coords, radial])

    if spline_transformer is None:
        spline_transformer = SplineTransformer(
            n_knots=8,
            degree=3,
            include_bias=False,
            extrapolation="continue",
        )
        spline_transformer.fit(times[:, None])
    temporal_basis = spline_transformer.transform(times[:, None])

    return {
        "knots": knots,
        "spatial_main": spatial_main,
        "radial": radial,
        "temporal_basis": temporal_basis,
        "spline_transformer": spline_transformer,
    }


def assemble_design(spatial_main: np.ndarray, radial: np.ndarray, temporal_basis: np.ndarray, location_idx: np.ndarray, time_idx: np.ndarray) -> np.ndarray:
    spatial_obs = spatial_main[location_idx]
    temporal_obs = temporal_basis[time_idx]
    interaction = radial[location_idx, :, None] * temporal_obs[:, None, :]
    return np.hstack([spatial_obs, temporal_obs, interaction.reshape(len(location_idx), -1)])


def fit_ridge_model(X: np.ndarray, y: np.ndarray, alpha: float) -> tuple:
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    model = Ridge(alpha=alpha, fit_intercept=True)
    model.fit(X_scaled, y)
    return model, scaler


def cross_validate_knots(data: dict, knot_candidates: list[int], alpha_grid: list[float], seed: int) -> tuple:
    coords = data["coords"]
    times = data["times"]
    y = data["response"]
    location_id = data["location_id"]
    time_id = data["time_id"]

    base_transformer = SplineTransformer(
        n_knots=8,
        degree=3,
        include_bias=False,
        extrapolation="continue",
    )
    base_transformer.fit(times[:, None])
    splitter = GroupKFold(n_splits=5)
    cv_results = []
    best_choice = None

    for knot_count in knot_candidates:
        components = build_design_components(coords, times, knot_count, seed, spline_transformer=base_transformer)
        X = assemble_design(components["spatial_main"], components["radial"], components["temporal_basis"], location_id, time_id)
        for alpha in alpha_grid:
            fold_rmse = []
            for train_idx, test_idx in splitter.split(X, y, groups=location_id):
                model, scaler = fit_ridge_model(X[train_idx], y[train_idx], alpha)
                pred = model.predict(scaler.transform(X[test_idx]))
                fold_rmse.append(float(np.sqrt(mean_squared_error(y[test_idx], pred))))
            mean_rmse = float(np.mean(fold_rmse))
            se_rmse = float(np.std(fold_rmse, ddof=1) / np.sqrt(len(fold_rmse)))
            cv_results.append({
                "knot_count": knot_count,
                "alpha": alpha,
                "rmse_mean": mean_rmse,
                "rmse_se": se_rmse,
            })
            if best_choice is None or mean_rmse < best_choice["rmse_mean"]:
                best_choice = {
                    "knot_count": knot_count,
                    "alpha": alpha,
                    "rmse_mean": mean_rmse,
                }
    return cv_results, best_choice, base_transformer


def ridge_prediction_variance(X_scaled: np.ndarray, model: Ridge, alpha: float, residual_var: float) -> np.ndarray:
    xtx = X_scaled.T @ X_scaled
    identity = np.eye(xtx.shape[0])
    inv_term = np.linalg.pinv(xtx + alpha * identity)
    cov_beta = residual_var * inv_term @ xtx @ inv_term
    return cov_beta


def save_knot_selection_figure(cv_results: list[dict], output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.1), constrained_layout=True)
    alphas = sorted({row["alpha"] for row in cv_results})
    cmap = plt.get_cmap("cividis")
    for idx, alpha in enumerate(alphas):
        rows = [row for row in cv_results if row["alpha"] == alpha]
        rows = sorted(rows, key=lambda row: row["knot_count"])
        ax.errorbar(
            [row["knot_count"] for row in rows],
            [row["rmse_mean"] for row in rows],
            yerr=[row["rmse_se"] for row in rows],
            marker="o",
            color=cmap(0.2 + 0.6 * idx / max(len(alphas) - 1, 1)),
            linewidth=2,
            capsize=4,
            label=f"alpha={alpha:g}",
        )
    ax.set_xlabel("Number of spatial knots")
    ax.set_ylabel("Cross-validated RMSE")
    ax.set_title("Cross-validation for knot selection")
    ax.legend(loc="best")
    fig.savefig(output_path, dpi=320)
    plt.close(fig)


def save_temporal_trend_figure(times: np.ndarray, future_times: np.ndarray, observed_mean: np.ndarray, trend: np.ndarray, ci_low: np.ndarray, ci_high: np.ndarray, output_path: Path) -> None:
    all_times = np.concatenate([times, future_times])
    fig, ax = plt.subplots(figsize=(10, 5.4), constrained_layout=True)
    ax.plot(times, observed_mean, color=plt.get_cmap("viridis")(0.35), linewidth=1.8, label="Observed mean log incidence")
    ax.plot(all_times, trend, color=plt.get_cmap("cividis")(0.85), linewidth=2.3, label="Predicted mean log incidence")
    ax.fill_between(all_times, ci_low, ci_high, color=plt.get_cmap("cividis")(0.55), alpha=0.25, label="Approx. 95% CI")
    ax.axvline(times[-1], color="black", linestyle="--", linewidth=1.2, label="Forecast start")
    ax.set_xlabel("Month index")
    ax.set_ylabel("Log incidence rate")
    ax.set_title("Spatiotemporal trend and forecast")
    ax.legend(loc="best")
    fig.savefig(output_path, dpi=320)
    plt.close(fig)


def save_risk_maps_figure(coords: np.ndarray, map_predictions: list[tuple[int, np.ndarray]], output_path: Path) -> None:
    grid_x, grid_y = np.meshgrid(np.linspace(0, 1, 60), np.linspace(0, 1, 60))
    grid_coords = np.column_stack([grid_x.ravel(), grid_y.ravel()])

    fig, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
    axes = axes.ravel()
    vmin = min(pred.min() for _, pred in map_predictions)
    vmax = max(pred.max() for _, pred in map_predictions)

    for ax, (time_point, prediction) in zip(axes, map_predictions):
        contour = ax.tricontourf(
            coords[:, 0],
            coords[:, 1],
            prediction,
            levels=15,
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
        )
        ax.scatter(coords[:, 0], coords[:, 1], s=12, c="black", alpha=0.35)
        ax.set_title(f"Month {time_point}")
        ax.set_xlabel("X coordinate")
        ax.set_ylabel("Y coordinate")
    fig.colorbar(contour, ax=axes.tolist(), label="Predicted log incidence")
    fig.savefig(output_path, dpi=320)
    plt.close(fig)


def summarize_to_stdout(results: dict) -> None:
    print("Spatiotemporal spline model summary")
    print(f"Selected spatial knots: {results['model_selection']['selected_knot_count']}")
    print(f"Selected ridge alpha: {results['model_selection']['selected_alpha']}")
    print(f"Training RMSE: {results['model_fit']['training_rmse']:.4f}")
    print(f"Forecast mean log incidence (next 12 months): {results['forecast']['mean_future_log_incidence']:.4f}")
    print(f"Peak forecast month: {results['forecast']['peak_month_index']}")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    ensure_directories(root)

    data = simulate_spatiotemporal_data(seed=SEED)
    knot_candidates = [4, 6, 8, 10, 12]
    alpha_grid = [0.1, 1.0, 10.0]
    cv_results, best_choice, base_transformer = cross_validate_knots(data, knot_candidates, alpha_grid, seed=SEED)

    components = build_design_components(
        data["coords"],
        data["times"],
        best_choice["knot_count"],
        seed=SEED,
        spline_transformer=base_transformer,
    )

    X = assemble_design(
        components["spatial_main"],
        components["radial"],
        components["temporal_basis"],
        data["location_id"],
        data["time_id"],
    )
    model, scaler = fit_ridge_model(X, data["response"], best_choice["alpha"])
    X_scaled = scaler.transform(X)
    fitted = model.predict(X_scaled)
    residuals = data["response"] - fitted
    training_rmse = float(np.sqrt(mean_squared_error(data["response"], fitted)))
    residual_var = float(np.var(residuals, ddof=X_scaled.shape[1] + 1))
    cov_beta = ridge_prediction_variance(X_scaled, model, best_choice["alpha"], residual_var)

    observed_mean = data["response"].reshape(len(data["coords"]), N_MONTHS).mean(axis=0)

    all_times = np.concatenate([data["times"], data["future_times"]])
    temporal_all = components["spline_transformer"].transform(all_times[:, None])
    location_sequence = np.repeat(np.arange(len(data["coords"])), len(all_times))
    time_sequence = np.tile(np.arange(len(all_times)), len(data["coords"]))
    X_all = assemble_design(components["spatial_main"], components["radial"], temporal_all, location_sequence, time_sequence)
    X_all_scaled = scaler.transform(X_all)
    pred_all = model.predict(X_all_scaled).reshape(len(data["coords"]), len(all_times))
    mean_trend = pred_all.mean(axis=0)

    mean_design = np.vstack([X_all_scaled[np.arange(i, len(X_all_scaled), len(all_times))].mean(axis=0) for i in range(len(all_times))])
    trend_se = np.sqrt(np.clip(np.sum((mean_design @ cov_beta) * mean_design, axis=1), 0.0, None) + residual_var / len(data["coords"]))
    ci_low = mean_trend - 1.96 * trend_se
    ci_high = mean_trend + 1.96 * trend_se

    forecast = pred_all[:, -FUTURE_MONTHS:]
    map_time_points = [0, 24, 48, 71]
    map_predictions = [(time_point, pred_all[:, time_point]) for time_point in map_time_points]

    save_knot_selection_figure(cv_results, root / "figures" / "spatiotemporal_knot_selection.png")
    save_temporal_trend_figure(data["times"], data["future_times"], observed_mean, mean_trend, ci_low, ci_high, root / "figures" / "spatiotemporal_trend.png")
    save_risk_maps_figure(data["coords"], map_predictions, root / "figures" / "spatiotemporal_risk_maps.png")

    results = {
        "simulation": {
            "seed": SEED,
            "n_locations": int(len(data["coords"])),
            "n_months": N_MONTHS,
            "future_months": FUTURE_MONTHS,
        },
        "model_selection": {
            "selected_knot_count": int(best_choice["knot_count"]),
            "selected_alpha": float(best_choice["alpha"]),
            "cv_results": cv_results,
        },
        "model_fit": {
            "training_rmse": training_rmse,
            "residual_variance": residual_var,
            "design_matrix_columns": int(X.shape[1]),
        },
        "forecast": {
            "mean_future_log_incidence": float(forecast.mean()),
            "peak_month_index": int(data["future_times"][np.argmax(forecast.mean(axis=0))]),
            "monthly_mean_log_incidence": [float(x) for x in forecast.mean(axis=0)],
            "monthly_ci": [[float(lo), float(hi)] for lo, hi in zip(ci_low[-FUTURE_MONTHS:], ci_high[-FUTURE_MONTHS:])],
        },
        "temporal_trend": {
            "month_index": [int(x) for x in all_times],
            "predicted_mean_log_incidence": [float(x) for x in mean_trend],
            "ci_low": [float(x) for x in ci_low],
            "ci_high": [float(x) for x in ci_high],
        },
        "risk_map_time_points": map_time_points,
    }

    with open(root / "results" / "spatiotemporal_results.json", "w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)

    summarize_to_stdout(results)


if __name__ == "__main__":
    main()
