from __future__ import annotations

import csv
import json
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from scipy import optimize, special
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import pdist

try:
    from libpysal.weights import DistanceBand
except Exception:
    DistanceBand = None

try:
    from esda import Geary as PySALGeary
    from esda import Moran as PySALMoran
    from esda import Moran_Local as PySALMoranLocal
except Exception:
    PySALGeary = None
    PySALMoran = None
    PySALMoranLocal = None


@dataclass
class SpatialWeights:
    matrix: np.ndarray
    threshold: float
    method: str
    pysal_w: Any | None = None


COLOR_MAP = {
    "Not significant": "#bdbdbd",
    "HH": "#E69F00",
    "HL": "#D55E00",
    "LH": "#56B4E9",
    "LL": "#0072B2",
}


def ensure_directories(root: Path) -> dict[str, Path]:
    directories = {
        "figures": root / "figures",
        "results": root / "results",
        "data": root / "data",
        "logs": root / "logs",
    }
    for path in directories.values():
        path.mkdir(parents=True, exist_ok=True)
    return directories


def row_standardize(weights: np.ndarray) -> np.ndarray:
    weights = weights.astype(float, copy=True)
    row_sums = weights.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0.0] = 1.0
    return weights / row_sums


def build_distance_weights(coords: np.ndarray, threshold: float | None = None) -> SpatialWeights:
    if threshold is None:
        threshold = math.sqrt(2.0) + 0.05

    if DistanceBand is not None:
        pysal_w = DistanceBand(coords, threshold=threshold, binary=True, silence_warnings=True)
        pysal_w.transform = "r"
        weights = np.asarray(pysal_w.full()[0], dtype=float)
        method = "distance_band_pysal"
        return SpatialWeights(matrix=weights, threshold=threshold, method=method, pysal_w=pysal_w)

    pairwise = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=2))
    weights = (pairwise <= threshold).astype(float)
    np.fill_diagonal(weights, 0.0)
    weights = row_standardize(weights)
    return SpatialWeights(matrix=weights, threshold=threshold, method="distance_band_manual", pysal_w=None)


def fdr_bh(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    n = p_values.size
    order = np.argsort(p_values)
    ranks = np.arange(1, n + 1)
    sorted_p = p_values[order]
    adjusted = sorted_p * n / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result


def permutation_p_value(observed: float, simulated: np.ndarray) -> float:
    extreme = np.sum(np.abs(simulated) >= abs(observed))
    return float((extreme + 1.0) / (simulated.size + 1.0))


def permutation_interval(simulated: np.ndarray, level: float = 0.95) -> list[float]:
    alpha = 1.0 - level
    low = float(np.quantile(simulated, alpha / 2.0))
    high = float(np.quantile(simulated, 1.0 - alpha / 2.0))
    return [low, high]


def global_morans_i(values: np.ndarray, weights: np.ndarray, permutations: int, rng: np.random.Generator) -> dict[str, Any]:
    x = np.asarray(values, dtype=float)
    z = x - x.mean()
    n = x.size
    s0 = float(weights.sum())
    denom = float(np.dot(z, z))
    observed = float((n / s0) * (z @ weights @ z) / denom)

    sims = np.empty(permutations, dtype=float)
    for i in range(permutations):
        zp = rng.permutation(z)
        sims[i] = float((n / s0) * (zp @ weights @ zp) / denom)

    return {
        "I": observed,
        "expected": float(-1.0 / (n - 1.0)),
        "p_value": permutation_p_value(observed, sims),
        "permutation_ci_95": permutation_interval(sims),
        "permutations": permutations,
        "permutation_mean": float(sims.mean()),
        "permutation_sd": float(sims.std(ddof=1)),
    }


def gearys_c(values: np.ndarray, weights: np.ndarray, permutations: int, rng: np.random.Generator) -> dict[str, Any]:
    x = np.asarray(values, dtype=float)
    n = x.size
    s0 = float(weights.sum())
    denom = float(np.sum((x - x.mean()) ** 2))
    diff = x[:, None] - x[None, :]
    observed = float(((n - 1.0) / (2.0 * s0)) * np.sum(weights * diff**2) / denom)

    sims = np.empty(permutations, dtype=float)
    for i in range(permutations):
        xp = rng.permutation(x)
        diffp = xp[:, None] - xp[None, :]
        sims[i] = float(((n - 1.0) / (2.0 * s0)) * np.sum(weights * diffp**2) / denom)

    return {
        "C": observed,
        "expected": 1.0,
        "p_value": permutation_p_value(observed - 1.0, sims - 1.0),
        "permutation_ci_95": permutation_interval(sims),
        "permutations": permutations,
        "permutation_mean": float(sims.mean()),
        "permutation_sd": float(sims.std(ddof=1)),
    }


def manual_local_morans_i(
    values: np.ndarray,
    weights: np.ndarray,
    permutations: int,
    alpha: float,
    rng: np.random.Generator,
) -> dict[str, Any]:
    x = np.asarray(values, dtype=float)
    z = (x - x.mean()) / x.std(ddof=1)
    lag_z = weights @ z
    observed = z * lag_z

    sims = np.empty((permutations, x.size), dtype=float)
    for i in range(permutations):
        zp = rng.permutation(z)
        sims[i, :] = zp * (weights @ zp)

    p_values = (np.sum(np.abs(sims) >= np.abs(observed), axis=0) + 1.0) / (permutations + 1.0)
    adjusted = fdr_bh(p_values)
    significant = adjusted < alpha

    quadrant = np.full(x.size, "Not significant", dtype=object)
    quadrant[(z >= 0.0) & (lag_z >= 0.0) & significant] = "HH"
    quadrant[(z < 0.0) & (lag_z < 0.0) & significant] = "LL"
    quadrant[(z >= 0.0) & (lag_z < 0.0) & significant] = "HL"
    quadrant[(z < 0.0) & (lag_z >= 0.0) & significant] = "LH"

    return {
        "I": observed,
        "lag_z": lag_z,
        "z": z,
        "p_values": p_values,
        "p_adjusted": adjusted,
        "quadrant": quadrant,
        "significant": significant,
        "permutations": permutations,
    }


def local_morans_i(
    values: np.ndarray,
    weights: np.ndarray,
    pysal_w: Any | None,
    permutations: int,
    alpha: float,
    seed: int,
    rng: np.random.Generator,
) -> dict[str, Any]:
    manual = manual_local_morans_i(values, weights, permutations, alpha, rng)

    if PySALMoranLocal is None or pysal_w is None:
        manual["implementation"] = "manual"
        return manual

    try:
        local = PySALMoranLocal(values, pysal_w, permutations=permutations, seed=seed)
        adjusted = fdr_bh(local.p_sim)
        significant = adjusted < alpha
        x = np.asarray(values, dtype=float)
        z = (x - x.mean()) / x.std(ddof=1)
        lag_z = weights @ z
        quadrant = np.full(x.size, "Not significant", dtype=object)
        quadrant[(local.q == 1) & significant] = "HH"
        quadrant[(local.q == 2) & significant] = "LH"
        quadrant[(local.q == 3) & significant] = "LL"
        quadrant[(local.q == 4) & significant] = "HL"
        return {
            "I": np.asarray(local.Is, dtype=float),
            "lag_z": lag_z,
            "z": z,
            "p_values": np.asarray(local.p_sim, dtype=float),
            "p_adjusted": adjusted,
            "quadrant": quadrant,
            "significant": significant,
            "permutations": permutations,
            "implementation": "esda",
        }
    except Exception:
        manual["implementation"] = "manual"
        return manual


def spherical_variogram(h: np.ndarray, nugget: float, partial_sill: float, range_param: float) -> np.ndarray:
    hr = np.asarray(h, dtype=float) / max(range_param, 1e-8)
    gamma = np.where(
        hr <= 1.0,
        nugget + partial_sill * (1.5 * hr - 0.5 * hr**3),
        nugget + partial_sill,
    )
    return gamma


def exponential_variogram(h: np.ndarray, nugget: float, partial_sill: float, range_param: float) -> np.ndarray:
    h = np.asarray(h, dtype=float)
    return nugget + partial_sill * (1.0 - np.exp(-h / max(range_param, 1e-8)))


def matern_variogram(h: np.ndarray, nugget: float, partial_sill: float, range_param: float, nu: float) -> np.ndarray:
    h = np.asarray(h, dtype=float)
    scaled = np.maximum(h / max(range_param, 1e-8), 1e-10)
    coefficient = 1.0 / (2.0 ** (nu - 1.0) * special.gamma(nu))
    corr = coefficient * (scaled**nu) * special.kv(nu, scaled)
    corr = np.where(h == 0.0, 1.0, corr)
    corr = np.clip(corr, 0.0, 1.0)
    return nugget + partial_sill * (1.0 - corr)


def empirical_variogram(coords: np.ndarray, values: np.ndarray, n_lags: int = 12) -> dict[str, np.ndarray]:
    distances = pdist(coords)
    semivariance = 0.5 * pdist(values[:, None], metric="sqeuclidean")

    max_dist = float(distances.max() * 0.6)
    bins = np.linspace(0.0, max_dist, n_lags + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])

    gamma = np.full(n_lags, np.nan)
    counts = np.zeros(n_lags, dtype=int)
    for i in range(n_lags):
        mask = (distances >= bins[i]) & (distances < bins[i + 1])
        counts[i] = int(mask.sum())
        if counts[i] > 0:
            gamma[i] = float(semivariance[mask].mean())

    valid = (~np.isnan(gamma)) & (counts > 0)
    return {
        "distances": centers[valid],
        "semivariance": gamma[valid],
        "pairs": counts[valid],
        "distance_pairs_all": distances,
        "semivariance_all": semivariance,
    }


def fit_variogram_models(empirical: dict[str, np.ndarray]) -> dict[str, dict[str, Any]]:
    h = empirical["distances"]
    gamma_hat = empirical["semivariance"]
    weights = np.sqrt(empirical["pairs"])

    nugget0 = float(max(0.0, np.nanmin(gamma_hat) * 0.25))
    sill0 = float(max(np.nanmax(gamma_hat) - nugget0, 1e-6))
    range0 = float(max(np.nanmedian(h), 1e-6))
    sigma = 1.0 / np.maximum(weights, 1.0)

    fitted: dict[str, dict[str, Any]] = {}
    model_specs = {
        "spherical": (spherical_variogram, [nugget0, sill0, range0], ([0.0, 1e-8, 1e-8], [np.inf, np.inf, np.inf])),
        "exponential": (exponential_variogram, [nugget0, sill0, range0], ([0.0, 1e-8, 1e-8], [np.inf, np.inf, np.inf])),
        "matern": (matern_variogram, [nugget0, sill0, range0, 1.0], ([0.0, 1e-8, 1e-8, 0.2], [np.inf, np.inf, np.inf, 5.0])),
    }

    for name, (model_func, initial, bounds) in model_specs.items():
        try:
            params, _ = optimize.curve_fit(
                model_func,
                h,
                gamma_hat,
                p0=initial,
                bounds=bounds,
                sigma=sigma,
                maxfev=20000,
            )
            fitted_values = model_func(h, *params)
            rmse = float(np.sqrt(np.mean((gamma_hat - fitted_values) ** 2)))
            param_names = ["nugget", "partial_sill", "range"] if name != "matern" else ["nugget", "partial_sill", "range", "nu"]
            fitted[name] = {
                "params": {key: float(value) for key, value in zip(param_names, params)},
                "rmse": rmse,
                "fitted": fitted_values,
            }
        except Exception as exc:
            fallback_values = model_func(h, *initial)
            rmse = float(np.sqrt(np.mean((gamma_hat - fallback_values) ** 2)))
            param_names = ["nugget", "partial_sill", "range"] if name != "matern" else ["nugget", "partial_sill", "range", "nu"]
            fitted[name] = {
                "params": {key: float(value) for key, value in zip(param_names, initial)},
                "rmse": rmse,
                "fitted": fallback_values,
                "warning": f"Curve fit failed: {exc}",
            }

    best = min(fitted.items(), key=lambda item: item[1]["rmse"])[0]
    fitted["best_model"] = {"name": best}
    return fitted


def generate_synthetic_data(rows: int, cols: int, seed: int) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    y_idx, x_idx = np.indices((rows, cols))
    smooth_noise = gaussian_filter(rng.normal(size=(rows, cols)), sigma=1.5, mode="reflect")
    smooth_noise = (smooth_noise - smooth_noise.mean()) / smooth_noise.std(ddof=1)
    trend = 0.18 * (x_idx / max(cols - 1, 1)) - 0.10 * (y_idx / max(rows - 1, 1))
    latent = smooth_noise + trend

    relative_risk = np.exp(0.45 * latent)
    population = rng.integers(120, 260, size=(rows, cols))
    baseline_rate = 0.06
    expected = population * baseline_rate
    cases = rng.poisson(expected * relative_risk)
    risk = (cases + 0.5) / (expected + 0.5)
    log_risk = np.log(risk)

    coords = np.column_stack([x_idx.ravel(), y_idx.ravel()]).astype(float)
    return {
        "x": x_idx.ravel(),
        "y": y_idx.ravel(),
        "coords": coords,
        "population": population.ravel(),
        "expected": expected.ravel(),
        "cases": cases.ravel(),
        "relative_risk": relative_risk.ravel(),
        "risk": risk.ravel(),
        "log_risk": log_risk.ravel(),
        "grid_log_risk": log_risk,
    }


def save_dataset(data: dict[str, np.ndarray], output_path: Path) -> None:
    fieldnames = ["x", "y", "population", "expected_cases", "observed_cases", "relative_risk", "smoothed_risk", "log_risk"]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for i in range(data["coords"].shape[0]):
            writer.writerow(
                {
                    "x": int(data["x"][i]),
                    "y": int(data["y"][i]),
                    "population": int(data["population"][i]),
                    "expected_cases": float(data["expected"][i]),
                    "observed_cases": int(data["cases"][i]),
                    "relative_risk": float(data["relative_risk"][i]),
                    "smoothed_risk": float(data["risk"][i]),
                    "log_risk": float(data["log_risk"][i]),
                }
            )


def plot_moran_scatter(
    z: np.ndarray,
    lag_z: np.ndarray,
    quadrant: np.ndarray,
    output_path: Path,
) -> None:
    colors = [COLOR_MAP[label] for label in quadrant]
    coeffs = np.polyfit(z, lag_z, 1)
    x_line = np.linspace(z.min() - 0.2, z.max() + 0.2, 200)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(z, lag_z, c=colors, s=40, edgecolor="white", linewidth=0.5, alpha=0.85)
    ax.plot(x_line, coeffs[0] * x_line + coeffs[1], color="#000000", linewidth=1.5, label="OLS fit")
    ax.axhline(0.0, color="#666666", linestyle="--", linewidth=1.0)
    ax.axvline(0.0, color="#666666", linestyle="--", linewidth=1.0)
    ax.set_xlabel("Standardized disease risk")
    ax.set_ylabel("Spatial lag of standardized risk")
    ax.set_title("Moran Scatter Plot")
    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_MAP[label], label=label, markersize=8)
        for label in ["HH", "HL", "LH", "LL", "Not significant"]
    ]
    legend_handles.append(Line2D([0], [0], color="#000000", lw=1.5, label="OLS fit"))
    ax.legend(handles=legend_handles, loc="best", frameon=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_variogram(
    empirical: dict[str, np.ndarray],
    fitted_models: dict[str, dict[str, Any]],
    output_path: Path,
) -> None:
    distances = empirical["distances"]
    semivariance = empirical["semivariance"]
    dense_x = np.linspace(0.0, float(distances.max()) * 1.05, 300)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(distances, semivariance, color="#000000", s=45, label="Empirical variogram")

    model_drawers = {
        "spherical": (spherical_variogram, "#0072B2"),
        "exponential": (exponential_variogram, "#009E73"),
        "matern": (matern_variogram, "#CC79A7"),
    }
    for name, (model_func, color) in model_drawers.items():
        params = fitted_models[name]["params"]
        if name == "matern":
            y = model_func(dense_x, params["nugget"], params["partial_sill"], params["range"], params["nu"])
        else:
            y = model_func(dense_x, params["nugget"], params["partial_sill"], params["range"])
        label = f"{name.title()} (RMSE={fitted_models[name]['rmse']:.3f})"
        ax.plot(dense_x, y, color=color, linewidth=2.0, label=label)

    ax.set_xlabel("Lag distance")
    ax.set_ylabel("Semivariance")
    ax.set_title("Empirical Variogram and Fitted Models")
    ax.legend(loc="best", frameon=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_lisa_map(cluster_grid: np.ndarray, output_path: Path) -> None:
    categories = ["Not significant", "HH", "HL", "LH", "LL"]
    code_map = {label: index for index, label in enumerate(categories)}
    grid_codes = np.vectorize(code_map.get)(cluster_grid)
    cmap = ListedColormap([COLOR_MAP[label] for label in categories])

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(grid_codes, origin="lower", cmap=cmap, interpolation="nearest")
    ax.set_xlabel("X coordinate")
    ax.set_ylabel("Y coordinate")
    ax.set_title("LISA Cluster Map")
    ax.set_xticks(np.arange(cluster_grid.shape[1]))
    ax.set_yticks(np.arange(cluster_grid.shape[0]))
    ax.grid(which="major", color="white", linestyle="-", linewidth=0.35)
    ax.set_xlim(-0.5, cluster_grid.shape[1] - 0.5)
    ax.set_ylim(-0.5, cluster_grid.shape[0] - 0.5)
    legend_handles = [Patch(facecolor=COLOR_MAP[label], edgecolor="none", label=label) for label in categories]
    ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_results_payload(
    data: dict[str, np.ndarray],
    weights: SpatialWeights,
    moran: dict[str, Any],
    local: dict[str, Any],
    geary: dict[str, Any],
    empirical: dict[str, np.ndarray],
    fitted_models: dict[str, dict[str, Any]],
    alpha: float,
) -> dict[str, Any]:
    significant_cells = []
    for idx, label in enumerate(local["quadrant"]):
        if label == "Not significant":
            continue
        significant_cells.append(
            {
                "index": int(idx),
                "x": int(data["x"][idx]),
                "y": int(data["y"][idx]),
                "cluster": str(label),
                "local_I": float(local["I"][idx]),
                "p_value": float(local["p_values"][idx]),
                "p_adjusted_fdr": float(local["p_adjusted"][idx]),
            }
        )

    cluster_counts = {label: int(np.sum(local["quadrant"] == label)) for label in ["HH", "HL", "LH", "LL", "Not significant"]}
    pysal_used = bool(weights.pysal_w is not None and (PySALMoran is not None or PySALMoranLocal is not None or PySALGeary is not None))

    payload = {
        "metadata": {
            "grid_shape": [int(data["grid_log_risk"].shape[0]), int(data["grid_log_risk"].shape[1])],
            "weights_method": weights.method,
            "distance_threshold": float(weights.threshold),
            "alpha": float(alpha),
            "pysal_available": pysal_used,
            "variable_analyzed": "log_risk",
        },
        "global_moran": moran,
        "geary_c": geary,
        "local_moran": {
            "implementation": local["implementation"],
            "cluster_counts": cluster_counts,
            "significant_cells": significant_cells,
        },
        "variogram": {
            "empirical": {
                "lag_distance": [float(x) for x in empirical["distances"]],
                "semivariance": [float(x) for x in empirical["semivariance"]],
                "pair_counts": [int(x) for x in empirical["pairs"]],
            },
            "models": {
                name: {
                    key: value if key != "fitted" else [float(v) for v in value]
                    for key, value in details.items()
                }
                for name, details in fitted_models.items()
                if name != "best_model"
            },
            "best_model": fitted_models["best_model"]["name"],
        },
        "summary": {
            "mean_observed_cases": float(np.mean(data["cases"])),
            "mean_relative_risk": float(np.mean(data["relative_risk"])),
            "mean_log_risk": float(np.mean(data["log_risk"])),
        },
    }

    if weights.pysal_w is not None and PySALMoran is not None:
        try:
            moran_obj = PySALMoran(data["log_risk"], weights.pysal_w, permutations=moran["permutations"])
            payload["global_moran"]["pysal_reference"] = {
                "I": float(moran_obj.I),
                "p_sim": float(moran_obj.p_sim),
            }
        except Exception:
            pass

    if weights.pysal_w is not None and PySALGeary is not None:
        try:
            geary_obj = PySALGeary(data["log_risk"], weights.pysal_w, permutations=geary["permutations"])
            payload["geary_c"]["pysal_reference"] = {
                "C": float(geary_obj.C),
                "p_sim": float(geary_obj.p_sim),
            }
        except Exception:
            pass

    return payload


def write_json(data: dict[str, Any], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)


def write_preprocessing_log(output_path: Path, rows: int, cols: int, seed: int) -> None:
    text = f"""# Preprocessing Log

- Dataset type: Synthetic spatial disease-risk grid
- Coordinate system: Cartesian grid coordinates (unit spacing)
- Grid size: {rows} x {cols}
- Random seeds: numpy={seed}, random={seed}
- Spatial autocorrelation generation: Gaussian-smoothed latent field plus mild linear trend
- Disease process: Poisson cases simulated from expected counts multiplied by latent relative risk
- Analysis variable: log-transformed risk ratio using continuity correction
- Missing data: none
- Spatial weights: distance band approximating Queen adjacency with threshold sqrt(2)+0.05
- Local inference correction: Benjamini-Hochberg FDR
- Variogram assumption: isotropic semivariance model fitted by weighted least squares
"""
    output_path.write_text(text, encoding="utf-8")


def main() -> None:
    seed = 42
    rows, cols = 12, 12
    global_permutations = 999
    local_permutations = 499
    alpha = 0.05

    np.random.seed(seed)
    random.seed(seed)
    rng = np.random.default_rng(seed)

    root = Path(__file__).resolve().parents[1]
    paths = ensure_directories(root)

    data = generate_synthetic_data(rows=rows, cols=cols, seed=seed)
    weights = build_distance_weights(data["coords"])

    moran = global_morans_i(data["log_risk"], weights.matrix, global_permutations, rng)
    geary = gearys_c(data["log_risk"], weights.matrix, global_permutations, rng)
    local = local_morans_i(
        data["log_risk"],
        weights.matrix,
        weights.pysal_w,
        local_permutations,
        alpha,
        seed,
        rng,
    )

    empirical = empirical_variogram(data["coords"], data["log_risk"], n_lags=12)
    fitted_models = fit_variogram_models(empirical)

    dataset_path = paths["data"] / "spatial_disease_data.csv"
    preprocessing_log_path = paths["data"] / "preprocessing-log.md"
    save_dataset(data, dataset_path)
    write_preprocessing_log(preprocessing_log_path, rows, cols, seed)

    moran_plot_path = paths["figures"] / "spatial_morans_i.png"
    variogram_plot_path = paths["figures"] / "spatial_variogram.png"
    lisa_plot_path = paths["figures"] / "spatial_lisa_map.png"

    plot_moran_scatter(local["z"], local["lag_z"], local["quadrant"], moran_plot_path)
    plot_variogram(empirical, fitted_models, variogram_plot_path)
    plot_lisa_map(local["quadrant"].reshape(rows, cols), lisa_plot_path)

    results = build_results_payload(data, weights, moran, local, geary, empirical, fitted_models, alpha)
    results_path = paths["results"] / "spatial_autocorrelation_results.json"
    write_json(results, results_path)

    print("Spatial autocorrelation analysis completed.")
    print(f"Global Moran's I: {moran['I']:.4f} (p={moran['p_value']:.4f})")
    print(f"Geary's C: {geary['C']:.4f} (p={geary['p_value']:.4f})")
    print(f"Best variogram model: {fitted_models['best_model']['name']}")
    print("Significant LISA clusters:")
    for label in ["HH", "HL", "LH", "LL"]:
        print(f"  {label}: {np.sum(local['quadrant'] == label)}")
    print(f"Results saved to: {results_path}")
    print(f"Figures saved to: {paths['figures']}")


if __name__ == "__main__":
    main()
