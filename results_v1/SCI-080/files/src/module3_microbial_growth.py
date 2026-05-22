from __future__ import annotations

import json
import math
import random
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

DPI = 300
SAFETY_THRESHOLD = 6.0
TEMPERATURES = np.array([5, 10, 15, 20, 25, 30, 37], dtype=float)

BASE_DIR = Path(__file__).resolve().parents[1]
FIGURES_DIR = BASE_DIR / "figures"
RESULTS_DIR = BASE_DIR / "results"
DATA_DIR = BASE_DIR / "data"
LOGS_DIR = BASE_DIR / "logs"
REPORT_PATH = BASE_DIR / "report.md"
PROCESS_LOG_PATH = LOGS_DIR / "process-log.jsonl"
PREPROCESSING_LOG_PATH = DATA_DIR / "preprocessing-log.md"
STAT_SUMMARY_PATH = RESULTS_DIR / "statistical-summary.md"
METRICS_PATH = RESULTS_DIR / "module3_metrics.json"
CURVES_PATH = DATA_DIR / "growth_curves.csv"


@dataclass(frozen=True)
class OrganismConfig:
    key: str
    display_name: str
    y0: float
    ymax: float
    h0: float
    pH_min: float
    pH_ref: float
    aw_min: float
    ratkowsky_guess: Tuple[float, float, float, float]
    mu_reference: Tuple[float, ...]


ORGANISMS: Dict[str, OrganismConfig] = {
    "salmonella": OrganismConfig(
        key="salmonella",
        display_name="Salmonella",
        y0=2.0,
        ymax=8.6,
        h0=2.2,
        pH_min=4.0,
        pH_ref=7.0,
        aw_min=0.94,
        ratkowsky_guess=(0.012, 4.0, 0.080, 47.0),
        mu_reference=(0.0, 0.015, 0.045, 0.110, 0.220, 0.380, 0.620),
    ),
    "e_coli": OrganismConfig(
        key="e_coli",
        display_name="E. coli",
        y0=2.0,
        ymax=8.8,
        h0=1.8,
        pH_min=4.3,
        pH_ref=7.0,
        aw_min=0.95,
        ratkowsky_guess=(0.012, 5.0, 0.070, 49.0),
        mu_reference=(0.0, 0.010, 0.040, 0.100, 0.240, 0.420, 0.720),
    ),
    "listeria": OrganismConfig(
        key="listeria",
        display_name="Listeria monocytogenes",
        y0=2.0,
        ymax=8.4,
        h0=1.6,
        pH_min=4.4,
        pH_ref=7.0,
        aw_min=0.92,
        ratkowsky_guess=(0.011, -2.0, 0.090, 45.0),
        mu_reference=(0.006, 0.018, 0.050, 0.120, 0.240, 0.360, 0.480),
    ),
    "s_aureus": OrganismConfig(
        key="s_aureus",
        display_name="S. aureus",
        y0=2.0,
        ymax=8.3,
        h0=2.0,
        pH_min=4.2,
        pH_ref=7.0,
        aw_min=0.86,
        ratkowsky_guess=(0.010, 6.0, 0.080, 46.0),
        mu_reference=(0.0, 0.008, 0.030, 0.090, 0.180, 0.280, 0.400),
    ),
}


FIGURE_FILES = {
    "baranyi": FIGURES_DIR / "fig3_baranyi_curves.png",
    "comparison": FIGURES_DIR / "fig3b_organism_comparison.png",
    "coldchain": FIGURES_DIR / "fig3c_coldchain_break.png",
    "montecarlo": FIGURES_DIR / "fig3d_monte_carlo.png",
    "boundary": FIGURES_DIR / "fig3e_growth_nogrowth.png",
}


class ProcessLogger:
    def __init__(self) -> None:
        self.events: List[dict] = []

    def record(
        self,
        phase: str,
        event_type: str,
        skill_or_tool: str,
        handoff_in: dict | None = None,
        handoff_out: dict | None = None,
        files_written: Iterable[str] | None = None,
        status: str = "ok",
    ) -> None:
        self.events.append(
            {
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "phase": phase,
                "event_type": event_type,
                "actor": "co-scientist",
                "skill_or_tool": skill_or_tool,
                "handoff_in": handoff_in or {},
                "handoff_out": handoff_out or {},
                "files_written": list(files_written or []),
                "status": status,
            }
        )

    def write(self, path: Path) -> None:
        with path.open("w", encoding="utf-8") as handle:
            for event in self.events:
                handle.write(json.dumps(event, ensure_ascii=False) + "\n")


LOGGER = ProcessLogger()


def ensure_directories() -> None:
    for directory in (FIGURES_DIR, RESULTS_DIR, DATA_DIR, LOGS_DIR):
        directory.mkdir(parents=True, exist_ok=True)


def ratkowsky_sqrt_model(temp_c: np.ndarray | float, b: float, t_min: float, c: float, t_max: float) -> np.ndarray:
    temp = np.asarray(temp_c, dtype=float)
    sqrt_mu = b * (temp - t_min) * (1.0 - np.exp(c * (temp - t_max)))
    sqrt_mu = np.where((temp <= t_min) | (temp >= t_max), 0.0, sqrt_mu)
    return np.clip(sqrt_mu, 0.0, None)


def ratkowsky_mu(temp_c: np.ndarray | float, params: Dict[str, float]) -> np.ndarray:
    sqrt_mu = ratkowsky_sqrt_model(temp_c, params["b"], params["Tmin"], params["c"], params["Tmax"])
    return np.square(sqrt_mu)


def baranyi_adjustment(t_h: np.ndarray | float, mu_max: float, h0: float) -> np.ndarray:
    time_h = np.asarray(t_h, dtype=float)
    if mu_max <= 0:
        return np.zeros_like(time_h)
    term = np.exp(-mu_max * time_h) + np.exp(-h0) - np.exp(-mu_max * time_h - h0)
    term = np.maximum(term, 1e-12)
    return time_h + np.log(term) / mu_max


def baranyi_growth(t_h: np.ndarray | float, y0: float, ymax: float, mu_max: float, h0: float) -> np.ndarray:
    time_h = np.asarray(t_h, dtype=float)
    if mu_max <= 0:
        return np.full_like(time_h, y0, dtype=float)
    adjustment = baranyi_adjustment(time_h, mu_max, h0)
    exp_term = np.exp(mu_max * adjustment)
    carrying = np.exp(ymax - y0)
    return y0 + mu_max * adjustment - np.log1p((exp_term - 1.0) / carrying)


def fit_secondary_model(config: OrganismConfig) -> Dict[str, object]:
    observed_mu = np.asarray(config.mu_reference, dtype=float)
    sqrt_obs = np.sqrt(np.clip(observed_mu, 0.0, None))
    bounds = ([0.001, -10.0, 0.001, 35.0], [0.200, 15.0, 1.000, 60.0])
    popt, pcov = curve_fit(
        ratkowsky_sqrt_model,
        TEMPERATURES,
        sqrt_obs,
        p0=config.ratkowsky_guess,
        bounds=bounds,
        maxfev=50000,
    )
    fitted_sqrt = ratkowsky_sqrt_model(TEMPERATURES, *popt)
    fitted_mu = fitted_sqrt**2
    residuals = observed_mu - fitted_mu
    rmse = float(np.sqrt(np.mean(np.square(residuals))))
    return {
        "params": {
            "b": float(popt[0]),
            "Tmin": float(popt[1]),
            "c": float(popt[2]),
            "Tmax": float(popt[3]),
        },
        "covariance": pcov,
        "observed_mu": observed_mu,
        "fitted_mu": fitted_mu,
        "rmse": rmse,
    }


def temperature_for_time(time_h: float, profile: List[Tuple[float, float, float]]) -> float:
    for start_h, end_h, temp_c in profile:
        if start_h <= time_h < end_h:
            return temp_c
    return profile[-1][2]


def simulate_temperature_profile(
    config: OrganismConfig,
    params: Dict[str, float],
    profile: List[Tuple[float, float, float]],
    dt_h: float = 0.1,
    y0_override: float | None = None,
    h0_override: float | None = None,
    ymax_override: float | None = None,
) -> pd.DataFrame:
    total_time_h = profile[-1][1]
    times = np.arange(0.0, total_time_h + dt_h, dt_h)
    counts = np.zeros_like(times)
    temperatures = np.zeros_like(times)
    counts[0] = config.y0 if y0_override is None else y0_override
    h_remaining = config.h0 if h0_override is None else h0_override
    ymax = config.ymax if ymax_override is None else ymax_override
    temperatures[0] = temperature_for_time(0.0, profile)

    for index in range(1, len(times)):
        current_time = times[index - 1]
        current_temp = temperature_for_time(current_time, profile)
        temperatures[index] = current_temp
        mu_max = float(ratkowsky_mu(current_temp, params))
        counts[index] = baranyi_growth(np.array([dt_h]), counts[index - 1], ymax, mu_max, h_remaining)[0]
        h_remaining = max(h_remaining - mu_max * dt_h, 0.0)

    return pd.DataFrame(
        {
            "time_h": times,
            "temperature_c": temperatures,
            "predicted_log_cfu_g": counts,
        }
    )


def generate_constant_growth_curves(fit_results: Dict[str, Dict[str, object]]) -> pd.DataFrame:
    records: List[pd.DataFrame] = []
    time_grid = np.linspace(0.0, 48.0, 241)
    for organism_key, config in ORGANISMS.items():
        params = fit_results[organism_key]["params"]
        for temp_c in TEMPERATURES:
            mu_max = float(ratkowsky_mu(temp_c, params))
            curve = pd.DataFrame(
                {
                    "organism": config.display_name,
                    "scenario": f"constant_{int(temp_c)}C",
                    "time_h": time_grid,
                    "temperature_c": temp_c,
                    "predicted_log_cfu_g": baranyi_growth(time_grid, config.y0, config.ymax, mu_max, config.h0),
                }
            )
            records.append(curve)
    return pd.concat(records, ignore_index=True)


def generate_cold_chain_scenarios(salmonella_params: Dict[str, float]) -> pd.DataFrame:
    config = ORGANISMS["salmonella"]
    scenarios = {
        "No abuse (4C)": [(0.0, 72.0, 4.0)],
        "Single break": [(0.0, 24.0, 4.0), (24.0, 48.0, 25.0), (48.0, 72.0, 8.0)],
        "Repeated breaks": [(0.0, 12.0, 4.0), (12.0, 24.0, 15.0), (24.0, 36.0, 4.0), (36.0, 54.0, 25.0), (54.0, 72.0, 10.0)],
    }
    outputs: List[pd.DataFrame] = []
    for scenario, profile in scenarios.items():
        curve = simulate_temperature_profile(config, salmonella_params, profile)
        curve.insert(0, "scenario", scenario)
        curve.insert(0, "organism", config.display_name)
        outputs.append(curve)
    return pd.concat(outputs, ignore_index=True)


def wilson_interval(successes: int, trials: int, z: float = 1.96) -> Tuple[float, float]:
    if trials == 0:
        return (0.0, 0.0)
    phat = successes / trials
    denominator = 1.0 + z**2 / trials
    center = (phat + z**2 / (2.0 * trials)) / denominator
    margin = z * math.sqrt((phat * (1.0 - phat) + z**2 / (4.0 * trials)) / trials) / denominator
    return (max(0.0, center - margin), min(1.0, center + margin))


def nearest_psd(matrix: np.ndarray) -> np.ndarray:
    symmetric = (matrix + matrix.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    eigenvalues = np.clip(eigenvalues, 1e-12, None)
    return eigenvectors @ np.diag(eigenvalues) @ eigenvectors.T


def run_monte_carlo(
    config: OrganismConfig,
    fit_result: Dict[str, object],
    profile: List[Tuple[float, float, float]],
    n_iter: int = 1000,
) -> Dict[str, object]:
    rng = np.random.default_rng(SEED)
    mean = np.array(
        [
            fit_result["params"]["b"],
            fit_result["params"]["Tmin"],
            fit_result["params"]["c"],
            fit_result["params"]["Tmax"],
        ],
        dtype=float,
    )
    covariance = nearest_psd(np.asarray(fit_result["covariance"], dtype=float))

    final_counts: List[float] = []
    for _ in range(n_iter):
        accepted = False
        while not accepted:
            sampled = rng.multivariate_normal(mean, covariance)
            b, t_min, c, t_max = sampled
            accepted = b > 0 and c > 0.001 and (t_max - t_min) > 8.0
        params = {"b": float(b), "Tmin": float(t_min), "c": float(c), "Tmax": float(t_max)}
        y0 = float(np.clip(rng.normal(config.y0, 0.15), 0.5, 4.0))
        ymax = float(np.clip(rng.normal(config.ymax, 0.18), y0 + 2.5, 9.5))
        h0 = float(np.clip(rng.normal(config.h0, 0.30), 0.1, 5.0))
        curve = simulate_temperature_profile(
            config,
            params,
            profile,
            dt_h=0.1,
            y0_override=y0,
            h0_override=h0,
            ymax_override=ymax,
        )
        final_counts.append(float(curve["predicted_log_cfu_g"].iloc[-1]))

    final_counts_array = np.asarray(final_counts)
    exceedances = int(np.sum(final_counts_array >= SAFETY_THRESHOLD))
    probability = float(exceedances / n_iter)
    ci_low, ci_high = wilson_interval(exceedances, n_iter)
    return {
        "iterations": n_iter,
        "final_counts": final_counts_array,
        "exceedance_probability": probability,
        "ci95": [float(ci_low), float(ci_high)],
        "summary": {
            "mean": float(np.mean(final_counts_array)),
            "sd": float(np.std(final_counts_array, ddof=1)),
            "p05": float(np.percentile(final_counts_array, 5)),
            "p50": float(np.percentile(final_counts_array, 50)),
            "p95": float(np.percentile(final_counts_array, 95)),
        },
    }


def growth_no_growth_boundary(config: OrganismConfig, fit_params: Dict[str, float], pH_values: np.ndarray) -> np.ndarray:
    boundary = fit_params["Tmin"] + 5.5 * np.power(np.maximum(config.pH_ref - pH_values, 0.0), 1.35)
    boundary = np.where(pH_values < config.pH_min, np.nan, boundary)
    return boundary


def save_figures(
    constant_curves: pd.DataFrame,
    cold_chain_df: pd.DataFrame,
    fit_results: Dict[str, Dict[str, object]],
    mc_results: Dict[str, object],
) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")

    salmonella_curves = constant_curves[constant_curves["organism"] == "Salmonella"]
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(TEMPERATURES)))
    for color, temp_c in zip(colors, TEMPERATURES):
        subset = salmonella_curves[salmonella_curves["scenario"] == f"constant_{int(temp_c)}C"]
        ax.plot(subset["time_h"], subset["predicted_log_cfu_g"], color=color, linewidth=2, label=f"{int(temp_c)}°C")
    ax.set_title("Baranyi growth curves for Salmonella")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Predicted load (log CFU/g)")
    ax.legend(title="Temperature", ncol=2, frameon=True)
    fig.tight_layout()
    fig.savefig(FIGURE_FILES["baranyi"], dpi=DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    compare_df = constant_curves[constant_curves["scenario"] == "constant_25C"]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, compare_df["organism"].nunique()))
    for color, organism in zip(colors, compare_df["organism"].drop_duplicates()):
        subset = compare_df[compare_df["organism"] == organism]
        ax.plot(subset["time_h"], subset["predicted_log_cfu_g"], color=color, linewidth=2, label=organism)
    ax.set_title("Organism comparison at 25°C")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Predicted load (log CFU/g)")
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(FIGURE_FILES["comparison"], dpi=DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, cold_chain_df["scenario"].nunique()))
    for color, scenario in zip(colors, cold_chain_df["scenario"].drop_duplicates()):
        subset = cold_chain_df[cold_chain_df["scenario"] == scenario]
        ax.plot(subset["time_h"], subset["predicted_log_cfu_g"], color=color, linewidth=2, label=scenario)
    ax.axhline(SAFETY_THRESHOLD, color="crimson", linestyle="--", linewidth=1.5, label="Safety threshold")
    ax.set_title("Cold chain break impact on Salmonella")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("Predicted load (log CFU/g)")
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(FIGURE_FILES["coldchain"], dpi=DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    counts = np.asarray(mc_results["final_counts"])
    ax.hist(counts, bins=30, color=plt.cm.viridis(0.55), edgecolor="white")
    ax.axvline(SAFETY_THRESHOLD, color="crimson", linestyle="--", linewidth=2, label="Safety threshold")
    ax.set_title("Monte Carlo final load distribution")
    ax.set_xlabel("Final load (log CFU/g)")
    ax.set_ylabel("Frequency")
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(FIGURE_FILES["montecarlo"], dpi=DPI)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    pH_values = np.linspace(4.0, 7.2, 160)
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(ORGANISMS)))
    for color, (organism_key, config) in zip(colors, ORGANISMS.items()):
        boundary = growth_no_growth_boundary(config, fit_results[organism_key]["params"], pH_values)
        ax.plot(pH_values, boundary, color=color, linewidth=2, label=config.display_name)
    ax.set_title("Growth/no-growth boundary")
    ax.set_xlabel("pH")
    ax.set_ylabel("Boundary temperature (°C)")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(FIGURE_FILES["boundary"], dpi=DPI)
    plt.close(fig)


def write_preprocessing_log() -> None:
    PREPROCESSING_LOG_PATH.write_text(
        "# Preprocessing log\n\n"
        "1. Defined ComBase-inspired reference growth-rate points for four organisms across 5-37°C.\n"
        "2. Fitted Ratkowsky square-root parameters with scipy.optimize.curve_fit.\n"
        "3. Simulated Baranyi growth curves under constant and dynamic temperature profiles.\n"
        "4. Quantified uncertainty with 1000 Monte Carlo iterations (seed=42).\n"
        "5. Exported figures, JSON metrics, CSV trajectories, and report artifacts.\n",
        encoding="utf-8",
    )


def write_statistical_summary(fit_results: Dict[str, Dict[str, object]], mc_results: Dict[str, object]) -> None:
    lines = [
        "# Statistical summary\n",
        "## Secondary model fit quality\n",
        "| Organism | RMSE (μmax) | Tmin (°C) | Tmax (°C) |",
        "|---|---:|---:|---:|",
    ]
    for organism_key, config in ORGANISMS.items():
        params = fit_results[organism_key]["params"]
        lines.append(
            f"| {config.display_name} | {fit_results[organism_key]['rmse']:.4f} | {params['Tmin']:.2f} | {params['Tmax']:.2f} |"
        )

    ci_low, ci_high = mc_results["ci95"]
    summary = mc_results["summary"]
    lines.extend(
        [
            "\n## Monte Carlo uncertainty\n",
            f"- Iterations: {mc_results['iterations']}\n",
            f"- Probability of exceeding {SAFETY_THRESHOLD:.1f} log CFU/g: {mc_results['exceedance_probability']:.3f} (95% CI {ci_low:.3f}-{ci_high:.3f})\n",
            f"- Mean final count: {summary['mean']:.2f} log CFU/g\n",
            f"- Standard deviation: {summary['sd']:.2f} log CFU/g\n",
            f"- 5th/50th/95th percentiles: {summary['p05']:.2f} / {summary['p50']:.2f} / {summary['p95']:.2f} log CFU/g\n",
        ]
    )
    STAT_SUMMARY_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(metrics: Dict[str, object]) -> None:
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        "# DRAFT — NOT FOR DISTRIBUTION\n",
        "## Module 3 microbial growth prediction\n",
        f"Generated: {timestamp}\n",
        "### Methods\n",
        "- Implemented the Baranyi primary growth model exactly in its logarithmic form.\n",
        "- Estimated Ratkowsky secondary-model parameters from ComBase-inspired reference μmax points using scipy.optimize.curve_fit.\n",
        "- Simulated constant-temperature growth, cold-chain abuse scenarios, Monte Carlo uncertainty, and a temperature-pH growth/no-growth boundary.\n",
        "\n### Results\n",
        f"- Safety threshold: {SAFETY_THRESHOLD:.1f} log CFU/g.\n",
        f"- Salmonella threshold exceedance probability under the single-break profile: {metrics['monte_carlo']['single_break_probability']:.3f}.\n",
        f"- Monte Carlo 95% CI: {metrics['monte_carlo']['single_break_probability_ci95'][0]:.3f}-{metrics['monte_carlo']['single_break_probability_ci95'][1]:.3f}.\n",
        f"- Highest fitted growth rate at 37°C: {metrics['organisms']['e_coli']['reference_fitted_mu_37C']:.3f} log CFU/g/h for E. coli.\n",
        "\n### Discussion\n",
        "- This script uses literature-inspired, ComBase-style reference points rather than directly downloading proprietary ComBase records.\n",
        "- Dynamic temperature abuse accelerates growth relative to uninterrupted refrigeration, while organism-specific Tmin values shift the feasible growth boundary.\n",
        "- The Monte Carlo module provides risk-oriented uncertainty estimates suitable for hazard screening, not regulatory validation.\n",
        "\n### File inventory\n",
        "- `src/module3_microbial_growth.py`: executable model script.\n",
        "- `figures/fig3_baranyi_curves.png`: Salmonella growth curves across temperatures.\n",
        "- `figures/fig3b_organism_comparison.png`: organism comparison at 25°C.\n",
        "- `figures/fig3c_coldchain_break.png`: cold-chain abuse impact.\n",
        "- `figures/fig3d_monte_carlo.png`: Monte Carlo final-count distribution.\n",
        "- `figures/fig3e_growth_nogrowth.png`: temperature-pH growth boundary.\n",
        "- `results/module3_metrics.json`: fitted parameters and risk metrics.\n",
        "- `results/statistical-summary.md`: fit quality and Monte Carlo summary.\n",
        "- `data/growth_curves.csv`: exported prediction trajectories.\n",
        "- `data/preprocessing-log.md`: preprocessing trace.\n",
        "- `logs/process-log.jsonl`: execution log.\n",
    ]
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_metrics(
    fit_results: Dict[str, Dict[str, object]],
    cold_chain_df: pd.DataFrame,
    mc_results: Dict[str, object],
) -> Dict[str, object]:
    organisms_metrics = {}
    for organism_key, config in ORGANISMS.items():
        params = fit_results[organism_key]["params"]
        organisms_metrics[organism_key] = {
            **asdict(config),
            "ratkowsky_fit": params,
            "fit_rmse": float(fit_results[organism_key]["rmse"]),
            "reference_fitted_mu_37C": float(ratkowsky_mu(37.0, params)),
            "boundary_temperature_at_pH55": float(growth_no_growth_boundary(config, params, np.array([5.5]))[0]),
        }
        organisms_metrics[organism_key].pop("ratkowsky_guess", None)
        organisms_metrics[organism_key]["mu_reference"] = list(config.mu_reference)

    cold_chain_summary = {
        scenario: float(group["predicted_log_cfu_g"].iloc[-1])
        for scenario, group in cold_chain_df.groupby("scenario")
    }
    return {
        "seed": SEED,
        "safety_threshold_log_cfu_g": SAFETY_THRESHOLD,
        "organisms": organisms_metrics,
        "cold_chain_final_counts": cold_chain_summary,
        "monte_carlo": {
            "iterations": mc_results["iterations"],
            "single_break_probability": float(mc_results["exceedance_probability"]),
            "single_break_probability_ci95": mc_results["ci95"],
            "summary": mc_results["summary"],
        },
    }


def main() -> None:
    ensure_directories()
    LOGGER.record(
        phase="plan",
        event_type="run_started",
        skill_or_tool="co-scientist-data-analysis",
        handoff_in={"seed": SEED, "temperatures": TEMPERATURES.tolist(), "iterations": 1000},
    )
    LOGGER.record(
        phase="plan",
        event_type="prompt_received",
        skill_or_tool="user_request",
        handoff_in={"task": "module3 microbial growth modeling"},
    )
    LOGGER.record(
        phase="plan",
        event_type="skill_selected",
        skill_or_tool="co-scientist-data-analysis",
        handoff_out={"reason": "Baranyi modeling, simulation, and visualization task"},
    )

    LOGGER.record(phase="execute", event_type="handoff_started", skill_or_tool="curve_fit")
    fit_results = {organism_key: fit_secondary_model(config) for organism_key, config in ORGANISMS.items()}
    LOGGER.record(
        phase="execute",
        event_type="handoff_completed",
        skill_or_tool="curve_fit",
        handoff_out={
            organism_key: fit_results[organism_key]["params"] for organism_key in ORGANISMS
        },
    )

    constant_curves = generate_constant_growth_curves(fit_results)
    salmonella_params = fit_results["salmonella"]["params"]
    cold_chain_df = generate_cold_chain_scenarios(salmonella_params)
    combined_curves = pd.concat([constant_curves, cold_chain_df], ignore_index=True)
    combined_curves.to_csv(CURVES_PATH, index=False)
    LOGGER.record(
        phase="execute",
        event_type="file_written",
        skill_or_tool="pandas.to_csv",
        files_written=[str(CURVES_PATH)],
    )

    single_break_profile = [(0.0, 24.0, 4.0), (24.0, 48.0, 25.0), (48.0, 72.0, 8.0)]
    mc_results = run_monte_carlo(ORGANISMS["salmonella"], fit_results["salmonella"], single_break_profile, n_iter=1000)

    save_figures(constant_curves, cold_chain_df, fit_results, mc_results)
    LOGGER.record(
        phase="execute",
        event_type="file_written",
        skill_or_tool="matplotlib",
        files_written=[str(path) for path in FIGURE_FILES.values()],
    )

    metrics = build_metrics(fit_results, cold_chain_df, mc_results)
    METRICS_PATH.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    LOGGER.record(
        phase="verify",
        event_type="file_written",
        skill_or_tool="json.dump",
        files_written=[str(METRICS_PATH)],
    )

    write_preprocessing_log()
    LOGGER.record(
        phase="verify",
        event_type="file_written",
        skill_or_tool="preprocessing-log",
        files_written=[str(PREPROCESSING_LOG_PATH)],
    )

    write_statistical_summary(fit_results, mc_results)
    LOGGER.record(
        phase="verify",
        event_type="file_written",
        skill_or_tool="statistical-summary",
        files_written=[str(STAT_SUMMARY_PATH)],
    )

    write_report(metrics)
    LOGGER.record(
        phase="report",
        event_type="report_finalized",
        skill_or_tool="report-writer",
        files_written=[str(REPORT_PATH)],
    )

    LOGGER.record(
        phase="log",
        event_type="run_completed",
        skill_or_tool="module3_microbial_growth.py",
        handoff_out={
            "metrics_file": str(METRICS_PATH),
            "csv_file": str(CURVES_PATH),
            "figure_count": len(FIGURE_FILES),
        },
    )
    LOGGER.write(PROCESS_LOG_PATH)


if __name__ == "__main__":
    main()
