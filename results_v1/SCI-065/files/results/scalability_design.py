#!/usr/bin/env python3
"""Compare scale-up strategies for brain organoid bioreactors."""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = BASE_DIR / "figures"
OUTPUT_CSV = RESULTS_DIR / "scalability_comparison.csv"


@dataclass
class ReactorMode:
    """Bioreactor operating mode definition."""

    name: str
    working_volume_l: float
    configuration: str
    organoid_density_per_ml: float
    recovery_efficiency: float
    medium_consumption_l_day: float
    supplement_cost_usd_l: float
    consumables_usd_batch: float
    size_cv_percent: float
    viability_percent: float
    impeller_diameter_m: float
    impeller_speed_s: float
    mean_velocity_m_s: float
    power_number: float
    limitation: str


MODES: List[ReactorMode] = [
    ReactorMode(
        name="Batch",
        working_volume_l=0.05,
        configuration="50 mL static culture, 12-well format",
        organoid_density_per_ml=2.6,
        recovery_efficiency=0.84,
        medium_consumption_l_day=0.020,
        supplement_cost_usd_l=290.0,
        consumables_usd_batch=34.0,
        size_cv_percent=19.0,
        viability_percent=83.0,
        impeller_diameter_m=0.015,
        impeller_speed_s=0.08,
        mean_velocity_m_s=1.2e-4,
        power_number=0.18,
        limitation="Labor-intensive handling and limited oxygen transfer",
    ),
    ReactorMode(
        name="Perfusion",
        working_volume_l=0.50,
        configuration="500 mL stirred tank with perfusion loop",
        organoid_density_per_ml=5.8,
        recovery_efficiency=0.88,
        medium_consumption_l_day=0.24,
        supplement_cost_usd_l=255.0,
        consumables_usd_batch=135.0,
        size_cv_percent=11.5,
        viability_percent=91.0,
        impeller_diameter_m=0.050,
        impeller_speed_s=1.20,
        mean_velocity_m_s=3.9e-3,
        power_number=1.15,
        limitation="Perfusion fouling risk and tubing dead volume",
    ),
    ReactorMode(
        name="Continuous",
        working_volume_l=5.0,
        configuration="5 L continuous perfusion with automated feeding",
        organoid_density_per_ml=7.2,
        recovery_efficiency=0.92,
        medium_consumption_l_day=1.80,
        supplement_cost_usd_l=232.0,
        consumables_usd_batch=780.0,
        size_cv_percent=8.5,
        viability_percent=94.0,
        impeller_diameter_m=0.120,
        impeller_speed_s=1.55,
        mean_velocity_m_s=7.8e-3,
        power_number=1.45,
        limitation="Automation complexity and sterility assurance burden",
    ),
]

RHO = 1000.0
MU = 1.0e-3
BASE_DIFFUSIVITY = 1.0e-9
REACTION_RATE = 6.0e-4
CHAR_LENGTH = 5.0e-4
MEDIUM_COST_BASE = 82.0
PROCESS_DURATION_D = 90.0


def effective_diffusivity(reynolds: float) -> float:
    """Approximate mixed diffusion enhancement."""
    return BASE_DIFFUSIVITY * (1.0 + 0.015 * np.sqrt(max(reynolds, 0.0)))


def compute_metrics(mode: ReactorMode, reference_yield: float) -> dict:
    """Calculate productivity, cost, and transport metrics."""
    yield_per_batch = mode.working_volume_l * 1000.0 * mode.organoid_density_per_ml * mode.recovery_efficiency
    total_medium_l = mode.medium_consumption_l_day * PROCESS_DURATION_D
    total_medium_cost = total_medium_l * (MEDIUM_COST_BASE + mode.supplement_cost_usd_l)
    total_cost = total_medium_cost + mode.consumables_usd_batch
    cost_per_organoid = total_cost / yield_per_batch
    reynolds = RHO * mode.impeller_speed_s * mode.impeller_diameter_m**2 / MU
    d_eff = effective_diffusivity(reynolds)
    damkohler = REACTION_RATE * CHAR_LENGTH**2 / d_eff
    peclet = mode.mean_velocity_m_s * CHAR_LENGTH / d_eff
    power_per_volume = mode.power_number * RHO * mode.impeller_speed_s**3 * mode.impeller_diameter_m**5 / mode.working_volume_l
    return {
        "mode": mode.name,
        "configuration": mode.configuration,
        "organoid_yield_per_batch": yield_per_batch,
        "medium_consumption_l_day": mode.medium_consumption_l_day,
        "cost_per_organoid_usd": cost_per_organoid,
        "size_uniformity_cv_percent": mode.size_cv_percent,
        "viability_percent": mode.viability_percent,
        "scale_up_factor": yield_per_batch / reference_yield,
        "reynolds_number": reynolds,
        "damkohler_number": damkohler,
        "peclet_number": peclet,
        "power_per_volume_w_m3": power_per_volume,
        "limitation": mode.limitation,
    }


def save_csv(rows: List[dict]) -> None:
    """Save scalability comparison to CSV."""
    header = list(rows[0].keys())
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        for row in rows:
            formatted = {
                key: (f"{value:.6f}" if isinstance(value, float) else value)
                for key, value in row.items()
            }
            writer.writerow(formatted)


def plot_scalability(rows: List[dict]) -> None:
    """Plot bar-chart comparison across production modes."""
    labels = [row["mode"] for row in rows]
    yield_values = np.array([row["organoid_yield_per_batch"] for row in rows])
    cost_values = np.array([row["cost_per_organoid_usd"] for row in rows])
    cv_values = np.array([row["size_uniformity_cv_percent"] for row in rows])
    viability_values = np.array([row["viability_percent"] for row in rows])

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    palettes = [plt.get_cmap("viridis"), plt.get_cmap("cividis"), plt.get_cmap("viridis"), plt.get_cmap("cividis")]
    datasets = [yield_values, cost_values, cv_values, viability_values]
    titles = ["Organoid yield per batch", "Cost per organoid", "Size uniformity (CV%)", "Viability"]
    ylabels = ["Organoids", "USD", "CV (%)", "Viability (%)"]
    for ax, cmap, data, title, ylabel in zip(axes.flatten(), palettes, datasets, titles, ylabels):
        colors = [cmap(v) for v in np.linspace(0.25, 0.85, len(labels))]
        ax.bar(labels, data, color=colors)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.25)
    fig.suptitle("Scalability comparison across brain organoid bioreactor modes", fontsize=14, weight="bold")
    fig.savefig(FIGURES_DIR / "scalability_comparison.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_dimensionless(rows: List[dict]) -> None:
    """Plot dimensionless transport and power metrics."""
    labels = [row["mode"] for row in rows]
    metrics = ["reynolds_number", "damkohler_number", "peclet_number", "power_per_volume_w_m3"]
    pretty = ["Reynolds", "Damkohler", "Peclet", "Power per volume"]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for idx, (ax, metric, title) in enumerate(zip(axes.flatten(), metrics, pretty)):
        values = np.array([row[metric] for row in rows])
        ax.plot(labels, values, marker="o", linewidth=2.5, color=plt.get_cmap("viridis")(0.2 + 0.2 * idx))
        ax.set_yscale("log")
        ax.set_title(title)
        ax.set_ylabel("Log scale")
        ax.grid(alpha=0.25, which="both")
    fig.suptitle("Dimensionless and energetic scale-up analysis", fontsize=14, weight="bold")
    fig.savefig(FIGURES_DIR / "dimensionless_analysis.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_cost_scaling(rows: List[dict]) -> None:
    """Plot cost-per-organoid scaling with production volume."""
    production = np.logspace(1.8, 5.0, 120)
    thresholds = [rows[0]["organoid_yield_per_batch"], rows[1]["organoid_yield_per_batch"]]
    costs = np.empty_like(production)
    for idx, count in enumerate(production):
        if count <= thresholds[0] * 2.0:
            costs[idx] = rows[0]["cost_per_organoid_usd"] * (count / thresholds[0]) ** -0.10
        elif count <= thresholds[1] * 2.5:
            costs[idx] = rows[1]["cost_per_organoid_usd"] * (count / thresholds[1]) ** -0.12
        else:
            costs[idx] = rows[2]["cost_per_organoid_usd"] * (count / rows[2]["organoid_yield_per_batch"]) ** -0.16
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.plot(production, costs, color="#2f4b7c", linewidth=2.8)
    ax.scatter(
        [row["organoid_yield_per_batch"] for row in rows],
        [row["cost_per_organoid_usd"] for row in rows],
        c=[0.25, 0.55, 0.85],
        cmap="viridis",
        s=70,
        edgecolor="black",
    )
    for row in rows:
        ax.annotate(row["mode"], (row["organoid_yield_per_batch"], row["cost_per_organoid_usd"]), xytext=(6, 6), textcoords="offset points")
    ax.set_xscale("log")
    ax.set_xlabel("Production volume (organoids per batch)")
    ax.set_ylabel("Cost per organoid (USD)")
    ax.set_title("Cost scaling with production volume")
    ax.grid(alpha=0.25, which="both")
    fig.savefig(FIGURES_DIR / "cost_scaling.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    """Run scalability analysis and save outputs."""
    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)
    reference_yield = MODES[0].working_volume_l * 1000.0 * MODES[0].organoid_density_per_ml * MODES[0].recovery_efficiency
    rows = [compute_metrics(mode, reference_yield) for mode in MODES]
    save_csv(rows)
    plot_scalability(rows)
    plot_dimensionless(rows)
    plot_cost_scaling(rows)
    print("Scalability analysis completed.")


if __name__ == "__main__":
    main()
