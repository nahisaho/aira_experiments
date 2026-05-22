#!/usr/bin/env python3
"""Brain organoid bioreactor shear stress and maturation modeling.

This script constructs phenomenological relationships between shear stress and
organoid maturation, quantifies an optimal operating window, and computes a
Pareto frontier balancing maturation gains against cell damage.

Outputs
-------
- results/shear_maturation_data.csv
- results/pareto_frontier.csv
- figures/shear_maturation_response.png
- figures/pareto_frontier.png
- figures/maturation_heatmap.png
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

np.random.seed(42)

BASE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = BASE_DIR / "figures"

plt.style.use("seaborn-v0_8-whitegrid")

SHEAR_RANGE = np.logspace(-3, 0, 320)
OPTIMAL_WINDOW = (0.01, 0.10)
TIME_DAYS = np.linspace(0.0, 60.0, 121)


def logistic(x: np.ndarray, x0: float, k: float) -> np.ndarray:
    """Standard logistic response."""
    return 1.0 / (1.0 + np.exp(-k * (x - x0)))


def log_gaussian(x: np.ndarray, center: float, width: float, amplitude: float = 1.0, baseline: float = 0.0) -> np.ndarray:
    """Gaussian response on a log10 shear axis."""
    z = (np.log10(x) - np.log10(center)) / width
    return baseline + amplitude * np.exp(-0.5 * z**2)


def compute_metrics(shear: np.ndarray) -> dict[str, np.ndarray]:
    """Generate maturation and damage metrics across shear stress."""
    map2 = np.clip(log_gaussian(shear, center=0.045, width=0.42, amplitude=0.72, baseline=0.25), 0.0, 1.0)
    tuj1 = np.clip(log_gaussian(shear, center=0.040, width=0.50, amplitude=0.68, baseline=0.28), 0.0, 1.0)
    gfap = np.clip(log_gaussian(shear, center=0.060, width=0.58, amplitude=0.50, baseline=0.20), 0.0, 1.0)

    viability_penalty = logistic(shear, x0=0.50, k=18.0)
    viability = np.clip(0.97 - 0.70 * viability_penalty, 0.15, 1.0)

    growth = np.clip(log_gaussian(shear, center=0.035, width=0.46, amplitude=0.62, baseline=0.30), 0.0, 1.0)
    electrophysiology = np.clip(
        0.18 + 0.82 * logistic(np.log10(shear), x0=np.log10(0.012), k=4.5) * (1.0 - logistic(shear, x0=0.35, k=16.0)),
        0.0,
        1.0,
    )

    marker_mean = (map2 + tuj1 + gfap) / 3.0
    maturation_score = np.clip(0.40 * marker_mean + 0.25 * growth + 0.35 * electrophysiology, 0.0, 1.0)
    cell_damage = np.clip(1.0 - viability, 0.0, 1.0)
    utility = maturation_score * viability

    return {
        "map2": map2,
        "tuj1": tuj1,
        "gfap": gfap,
        "viability": viability,
        "growth": growth,
        "electrophysiology": electrophysiology,
        "marker_mean": marker_mean,
        "maturation_score": maturation_score,
        "cell_damage": cell_damage,
        "utility": utility,
    }


def pareto_frontier(maturation: np.ndarray, damage: np.ndarray, shear: np.ndarray) -> np.ndarray:
    """Return non-dominated points maximizing maturation and minimizing damage."""
    keep = np.ones(shear.size, dtype=bool)
    for i in range(shear.size):
        dominates = (
            (maturation >= maturation[i])
            & (damage <= damage[i])
            & ((maturation > maturation[i]) | (damage < damage[i]))
        )
        dominates[i] = False
        if np.any(dominates):
            keep[i] = False
    frontier = np.column_stack([shear[keep], maturation[keep], damage[keep]])
    order = np.argsort(frontier[:, 2])
    return frontier[order]


def save_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    """Write rows to CSV."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def make_response_figure(shear: np.ndarray, metrics: dict[str, np.ndarray]) -> None:
    """Create a multi-panel response figure."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharex=True)
    cmap = plt.get_cmap("viridis")

    axes[0].semilogx(shear, metrics["map2"], color=cmap(0.15), linewidth=2.2, label="MAP2")
    axes[0].semilogx(shear, metrics["tuj1"], color=cmap(0.45), linewidth=2.2, label="TUJ1")
    axes[0].semilogx(shear, metrics["gfap"], color=cmap(0.75), linewidth=2.2, label="GFAP")
    axes[0].axvspan(*OPTIMAL_WINDOW, color="lightgray", alpha=0.3)
    axes[0].set_ylabel("Normalized expression")
    axes[0].set_title("Neural markers")
    axes[0].legend(frameon=True)

    axes[1].semilogx(shear, metrics["viability"], color=plt.get_cmap("cividis")(0.35), linewidth=2.2, label="Viability")
    axes[1].semilogx(shear, metrics["electrophysiology"], color=plt.get_cmap("cividis")(0.75), linewidth=2.2, label="Ephys maturity")
    axes[1].axvline(0.5, linestyle="--", color="black", linewidth=1.3, label="Damage threshold")
    axes[1].axvspan(*OPTIMAL_WINDOW, color="lightgray", alpha=0.3)
    axes[1].set_ylabel("Normalized response")
    axes[1].set_title("Viability and function")
    axes[1].legend(frameon=True)

    axes[2].semilogx(shear, metrics["growth"], color=cmap(0.30), linewidth=2.2, label="Growth")
    axes[2].semilogx(shear, metrics["maturation_score"], color=cmap(0.85), linewidth=2.2, label="Composite maturation")
    axes[2].axvspan(*OPTIMAL_WINDOW, color="lightgray", alpha=0.3)
    axes[2].set_ylabel("Normalized score")
    axes[2].set_title("Growth and composite score")
    axes[2].legend(frameon=True)

    for ax in axes:
        ax.set_xlabel("Shear stress (Pa)")
        ax.set_ylim(0.0, 1.05)

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "shear_maturation_response.png", dpi=300)
    plt.close(fig)


def make_pareto_figure(shear: np.ndarray, metrics: dict[str, np.ndarray], frontier: np.ndarray) -> None:
    """Plot Pareto-optimal operating conditions."""
    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    scatter = ax.scatter(metrics["cell_damage"], metrics["maturation_score"], c=np.log10(shear), cmap="viridis", s=24, alpha=0.8)
    ax.plot(frontier[:, 2], frontier[:, 1], color="black", linewidth=2.0, label="Pareto frontier")
    optimal_mask = (shear >= OPTIMAL_WINDOW[0]) & (shear <= OPTIMAL_WINDOW[1])
    ax.scatter(metrics["cell_damage"][optimal_mask], metrics["maturation_score"][optimal_mask], color="tomato", s=28, label="Literature-guided window")
    ax.set_xlabel("Cell damage index")
    ax.set_ylabel("Composite maturation score")
    ax.set_title("Pareto-optimal shear conditions")
    ax.legend(frameon=True)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("log10 shear stress (Pa)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "pareto_frontier.png", dpi=300)
    plt.close(fig)


def make_heatmap_figure(shear: np.ndarray, time_days: np.ndarray, metrics: dict[str, np.ndarray]) -> None:
    """Plot maturation score across shear stress and maturation time."""
    time_progress = 1.0 / (1.0 + np.exp(-(time_days - 24.0) / 6.5))
    maturation_surface = np.outer(metrics["maturation_score"], time_progress)

    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    mesh = ax.pcolormesh(time_days, shear, maturation_surface, shading="auto", cmap="cividis")
    ax.set_yscale("log")
    ax.set_xlabel("Culture time (days)")
    ax.set_ylabel("Shear stress (Pa)")
    ax.set_title("Maturation score landscape")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Composite maturation score")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "maturation_heatmap.png", dpi=300)
    plt.close(fig)


def main() -> None:
    """Generate shear-response tables and figures."""
    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)

    metrics = compute_metrics(SHEAR_RANGE)
    frontier = pareto_frontier(metrics["maturation_score"], metrics["cell_damage"], SHEAR_RANGE)

    in_window = (SHEAR_RANGE >= OPTIMAL_WINDOW[0]) & (SHEAR_RANGE <= OPTIMAL_WINDOW[1])
    best_idx = np.argmax(np.where(in_window, metrics["utility"], -np.inf))
    best_shear = SHEAR_RANGE[best_idx]

    data_rows = []
    for idx, shear in enumerate(SHEAR_RANGE):
        data_rows.append(
            {
                "shear_pa": f"{shear:.8f}",
                "map2": f"{metrics['map2'][idx]:.6f}",
                "tuj1": f"{metrics['tuj1'][idx]:.6f}",
                "gfap": f"{metrics['gfap'][idx]:.6f}",
                "viability": f"{metrics['viability'][idx]:.6f}",
                "growth": f"{metrics['growth'][idx]:.6f}",
                "electrophysiology": f"{metrics['electrophysiology'][idx]:.6f}",
                "maturation_score": f"{metrics['maturation_score'][idx]:.6f}",
                "cell_damage": f"{metrics['cell_damage'][idx]:.6f}",
                "utility": f"{metrics['utility'][idx]:.6f}",
                "in_literature_window": int(in_window[idx]),
                "is_best_window_point": int(idx == best_idx),
            }
        )

    frontier_rows = []
    for shear, maturation, damage in frontier:
        frontier_rows.append(
            {
                "shear_pa": f"{shear:.8f}",
                "maturation_score": f"{maturation:.6f}",
                "cell_damage": f"{damage:.6f}",
                "viability": f"{1.0 - damage:.6f}",
                "utility": f"{maturation * (1.0 - damage):.6f}",
                "in_literature_window": int(OPTIMAL_WINDOW[0] <= shear <= OPTIMAL_WINDOW[1]),
            }
        )

    save_csv(
        RESULTS_DIR / "shear_maturation_data.csv",
        data_rows,
        [
            "shear_pa",
            "map2",
            "tuj1",
            "gfap",
            "viability",
            "growth",
            "electrophysiology",
            "maturation_score",
            "cell_damage",
            "utility",
            "in_literature_window",
            "is_best_window_point",
        ],
    )
    save_csv(
        RESULTS_DIR / "pareto_frontier.csv",
        frontier_rows,
        ["shear_pa", "maturation_score", "cell_damage", "viability", "utility", "in_literature_window"],
    )

    make_response_figure(SHEAR_RANGE, metrics)
    make_pareto_figure(SHEAR_RANGE, metrics, frontier)
    make_heatmap_figure(SHEAR_RANGE, TIME_DAYS, metrics)

    print(f"Best literature-window shear: {best_shear:.4f} Pa")
    print("Saved shear stress maturation outputs.")


if __name__ == "__main__":
    main()
