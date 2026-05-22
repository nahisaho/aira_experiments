from __future__ import annotations

import json
import math
import random
import traceback
import warnings
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import minimize
from scipy.stats import binom, poisson, qmc, spearmanr

SEED = 42
np.random.seed(SEED)
random.seed(SEED)
warnings.filterwarnings("ignore")
sns.set_theme(style="whitegrid")
plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.dpi"] = 150
plt.rcParams["axes.prop_cycle"] = plt.cycler(color=sns.color_palette("viridis", 8))

BASE_DIR = Path(__file__).resolve().parent
FIG_DIR = BASE_DIR / "figures"
RESULTS_DIR = BASE_DIR / "results"
DATA_DIR = BASE_DIR / "data"
LOG_DIR = BASE_DIR / "logs"
LOG_PATH = LOG_DIR / "process-log.jsonl"
SUMMARY_PATH = RESULTS_DIR / "summary_metrics.json"

for directory in [FIG_DIR, RESULTS_DIR, DATA_DIR, LOG_DIR]:
    directory.mkdir(parents=True, exist_ok=True)


class ProcessLogger:
    def __init__(self, path: Path):
        self.path = path

    def log(
        self,
        phase: str,
        event_type: str,
        skill_or_tool: str,
        handoff_in: dict | None = None,
        handoff_out: dict | None = None,
        files_written: list[str] | None = None,
        status: str = "ok",
        message: str | None = None,
    ) -> None:
        payload = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "phase": phase,
            "event_type": event_type,
            "actor": "co-scientist",
            "skill_or_tool": skill_or_tool,
            "handoff_in": handoff_in or {},
            "handoff_out": handoff_out or {},
            "files_written": files_written or [],
            "status": status,
        }
        if message:
            payload["message"] = message
        with self.path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(payload, ensure_ascii=False) + "\n")


LOGGER = ProcessLogger(LOG_PATH)
LOGGER.log(
    phase="PLAN",
    event_type="run_started",
    skill_or_tool="adc_platform.py",
    handoff_in={"seed": SEED, "cwd": str(BASE_DIR)},
    handoff_out={"required_outputs": ["figures", "results", "data", "logs/process-log.jsonl"]},
)
LOGGER.log(
    phase="PLAN",
    event_type="prompt_received",
    skill_or_tool="adc_platform.py",
    handoff_in={
        "objective": "ADC payload-linker optimization platform",
        "modules": [
            "DAR distribution",
            "linker cleavage",
            "bystander diffusion",
            "optimization",
            "PKPD",
            "Monte Carlo sensitivity",
            "HER2 case study",
        ],
    },
)
LOGGER.log(
    phase="PLAN",
    event_type="skill_selected",
    skill_or_tool="co-scientist-data-analysis",
    handoff_out={"route": "simulation, optimization, visualization"},
)


def save_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False)
    LOGGER.log(
        phase="EXECUTE",
        event_type="file_written",
        skill_or_tool="pandas.to_csv",
        files_written=[str(path.relative_to(BASE_DIR))],
        handoff_out={"rows": int(len(df)), "columns": list(df.columns)},
    )



def save_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")
    LOGGER.log(
        phase="EXECUTE",
        event_type="file_written",
        skill_or_tool="write_text",
        files_written=[str(path.relative_to(BASE_DIR))],
    )



def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    LOGGER.log(
        phase="EXECUTE",
        event_type="file_written",
        skill_or_tool="matplotlib.savefig",
        files_written=[str(path.relative_to(BASE_DIR))],
    )



def bootstrap_ci(values: np.ndarray, n_boot: int = 1000) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return (float("nan"), float("nan"))
    boot = []
    for _ in range(n_boot):
        sample = rng.choice(values, size=values.size, replace=True)
        boot.append(float(np.mean(sample)))
    return tuple(np.percentile(boot, [2.5, 97.5]))



def safe_ratio(num: float, den: float) -> float:
    return float(num / den) if abs(den) > 1e-12 else float("nan")


summary_metrics: dict[str, dict] = {}
section_failures: list[str] = []



def run_section(name: str, func):
    print(f"[INFO] Running section: {name}")
    LOGGER.log(phase="EXECUTE", event_type="handoff_started", skill_or_tool=name)
    try:
        result = func()
        LOGGER.log(
            phase="VERIFY",
            event_type="handoff_completed",
            skill_or_tool=name,
            handoff_out={"status": "completed"},
        )
        print(f"[INFO] Completed section: {name}")
        return result
    except Exception as exc:  # noqa: BLE001
        section_failures.append(f"{name}: {exc}")
        LOGGER.log(
            phase="VERIFY",
            event_type="handoff_completed",
            skill_or_tool=name,
            status="error",
            message=str(exc),
            handoff_out={"traceback": traceback.format_exc()[-4000:]},
        )
        print(f"[ERROR] Section failed but execution continues: {name} -> {exc}")
        return None



def therapeutic_index_profile(dar_values: np.ndarray) -> np.ndarray:
    dar_values = np.asarray(dar_values, dtype=float)
    return 1.0 + 2.6 * np.exp(-0.5 * ((dar_values - 3.5) / 1.0) ** 2) - 0.12 * np.maximum(dar_values - 6.0, 0)



def section_dar_distribution() -> dict:
    species = np.arange(0, 9)
    n_sites = 8
    p_eff = 0.78
    rng = np.random.default_rng(SEED)
    binom_probs = binom.pmf(species, n_sites, p_eff)
    poisson_raw = poisson.pmf(species, mu=n_sites * p_eff)
    poisson_probs = poisson_raw / poisson_raw.sum()
    samples = rng.binomial(n_sites, p_eff, size=10_000)
    counts = np.bincount(samples, minlength=9)[:9]

    clearance_rate = 0.15 + 0.035 * species + 0.010 * np.maximum(species - 4, 0) ** 2
    hydrophobicity_penalty = 0.04 * species + 0.035 * np.maximum(species - 4, 0) ** 2
    therapeutic_index = therapeutic_index_profile(species)

    dar_df = pd.DataFrame(
        {
            "DAR": species,
            "binomial_probability": binom_probs,
            "poisson_probability": poisson_probs,
            "monte_carlo_count": counts,
            "monte_carlo_fraction": counts / counts.sum(),
            "clearance_rate_per_day": clearance_rate,
            "therapeutic_index_score": therapeutic_index,
            "hydrophobicity_penalty": hydrophobicity_penalty,
        }
    )
    save_csv(dar_df, RESULTS_DIR / "dar_analysis.csv")

    sample_df = pd.DataFrame({"sample_id": np.arange(1, len(samples) + 1), "DAR": samples})
    save_csv(sample_df, DATA_DIR / "dar_monte_carlo_samples.csv")

    fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
    width = 0.36
    axes[0].bar(species - width / 2, binom_probs, width=width, label="Binomial model", color=sns.color_palette("viridis", 5)[2])
    axes[0].bar(species + width / 2, poisson_probs, width=width, label="Poisson approximation", color=sns.color_palette("viridis", 5)[4], alpha=0.8)
    axes[0].set_ylabel("Probability")
    axes[0].set_title("DAR species distribution models")
    axes[0].legend()

    axes[1].hist(samples, bins=np.arange(-0.5, 9.5, 1), color=sns.color_palette("viridis", 6)[3], edgecolor="black", alpha=0.9)
    axes[1].axvspan(2.5, 4.5, color="orange", alpha=0.2, label="Therapeutic window (DAR 3-4)")
    axes[1].set_xlabel("Drug-to-antibody ratio (DAR)")
    axes[1].set_ylabel("Count")
    axes[1].set_title("Monte Carlo DAR distribution (n=10000)")
    ax2 = axes[1].twinx()
    ax2.plot(species, therapeutic_index, color="crimson", marker="o", label="Therapeutic index")
    ax2.set_ylabel("Therapeutic index score")
    axes[1].legend(loc="upper left")
    ax2.legend(loc="upper right")
    save_figure(fig, FIG_DIR / "01_dar_distribution.png")

    optimal_fraction = float(((samples >= 3) & (samples <= 4)).mean())
    mean_dar = float(np.mean(samples))
    summary_metrics["dar"] = {
        "mean_dar": mean_dar,
        "std_dar": float(np.std(samples, ddof=1)),
        "fraction_in_therapeutic_window": optimal_fraction,
        "target_dar": 8,
    }
    return summary_metrics["dar"]



def simulate_first_order_release(k_hr: float, t_hours: np.ndarray) -> np.ndarray:
    y0 = [1.0]

    def rhs(_, y):
        return [-k_hr * max(y[0], 0.0)]

    sol = solve_ivp(rhs, (float(t_hours[0]), float(t_hours[-1])), y0, t_eval=t_hours, method="RK45")
    return 1.0 - np.clip(sol.y[0], 0.0, 1.0)



def simulate_mm_release(vmax_uM_per_min: float, km_uM: float, scale: float, t_hours: np.ndarray) -> np.ndarray:
    s0 = 100.0
    y0 = [s0]
    vmax_hr = vmax_uM_per_min * 60.0 * scale

    def rhs(_, y):
        substrate = max(y[0], 0.0)
        rate = vmax_hr * substrate / (km_uM + substrate + 1e-12)
        return [-rate]

    sol = solve_ivp(rhs, (float(t_hours[0]), float(t_hours[-1])), y0, t_eval=t_hours, method="RK45", max_step=0.2)
    remaining = np.clip(sol.y[0], 0.0, s0)
    return 1.0 - remaining / s0



def section_linker_cleavage() -> dict:
    t_hours = np.linspace(0, 24, 241)
    records: list[dict] = []

    acid_params = {"k_max": 0.30, "pKa": 6.0, "sigma": 0.25, "plasma_pH": 7.4, "tumor_pH": 5.5}
    for compartment, pH in [("Plasma", acid_params["plasma_pH"]), ("Tumor/Lysosome", acid_params["tumor_pH"])]:
        k_hr = acid_params["k_max"] / (1.0 + math.exp(-(acid_params["pKa"] - pH) / acid_params["sigma"]))
        released = simulate_first_order_release(k_hr, t_hours)
        for t, rel in zip(t_hours, released, strict=False):
            records.append(
                {
                    "time_h": t,
                    "mechanism": "Acid-sensitive",
                    "compartment": compartment,
                    "released_fraction": rel,
                    "effective_rate_per_h": k_hr,
                }
            )

    enzyme_params = {"Vmax": 0.5, "Km": 50.0, "plasma_scale": 0.03, "tumor_scale": 1.0}
    for compartment, scale in [("Plasma", enzyme_params["plasma_scale"]), ("Tumor/Lysosome", enzyme_params["tumor_scale"])]:
        released = simulate_mm_release(enzyme_params["Vmax"], enzyme_params["Km"], scale, t_hours)
        effective_rate = enzyme_params["Vmax"] * scale / enzyme_params["Km"]
        for t, rel in zip(t_hours, released, strict=False):
            records.append(
                {
                    "time_h": t,
                    "mechanism": "Cathepsin B",
                    "compartment": compartment,
                    "released_fraction": rel,
                    "effective_rate_per_h": effective_rate,
                }
            )

    red_params = {"k_base": 0.35, "K_GSH": 0.5, "plasma_GSH_mM": 0.01, "tumor_GSH_mM": 5.0}
    for compartment, gsh in [("Plasma", red_params["plasma_GSH_mM"]), ("Tumor/Lysosome", red_params["tumor_GSH_mM"])]:
        k_hr = red_params["k_base"] * gsh / (red_params["K_GSH"] + gsh)
        released = simulate_first_order_release(k_hr, t_hours)
        for t, rel in zip(t_hours, released, strict=False):
            records.append(
                {
                    "time_h": t,
                    "mechanism": "Reductive disulfide",
                    "compartment": compartment,
                    "released_fraction": rel,
                    "effective_rate_per_h": k_hr,
                }
            )

    kinetics_df = pd.DataFrame(records)
    save_csv(kinetics_df, RESULTS_DIR / "linker_kinetics.csv")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    mechanisms = ["Acid-sensitive", "Cathepsin B", "Reductive disulfide"]
    for ax, mechanism in zip(axes, mechanisms, strict=False):
        subset = kinetics_df[kinetics_df["mechanism"] == mechanism]
        sns.lineplot(data=subset, x="time_h", y="released_fraction", hue="compartment", ax=ax, palette="viridis", linewidth=2)
        ax.set_title(mechanism)
        ax.set_xlabel("Time (h)")
        ax.set_ylabel("Released fraction")
        ax.set_ylim(0, 1.02)
        if ax is not axes[0]:
            ax.get_legend().remove()
    axes[0].legend(title="Compartment")
    fig.suptitle("Linker cleavage kinetics across compartments", y=1.02)
    save_figure(fig, FIG_DIR / "02_linker_cleavage_kinetics.png")

    end_release = (
        kinetics_df[kinetics_df["time_h"] == kinetics_df["time_h"].max()]
        .pivot(index="mechanism", columns="compartment", values="released_fraction")
        .reset_index()
    )
    save_csv(end_release, DATA_DIR / "linker_terminal_release_summary.csv")

    summary_metrics["linker"] = {
        "acid_tumor_release_24h": float(end_release.loc[end_release["mechanism"] == "Acid-sensitive", "Tumor/Lysosome"].iloc[0]),
        "enzyme_selectivity_ratio": float(
            end_release.loc[end_release["mechanism"] == "Cathepsin B", "Tumor/Lysosome"].iloc[0]
            / max(end_release.loc[end_release["mechanism"] == "Cathepsin B", "Plasma"].iloc[0], 1e-9)
        ),
        "reductive_plasma_release_24h": float(end_release.loc[end_release["mechanism"] == "Reductive disulfide", "Plasma"].iloc[0]),
    }
    return summary_metrics["linker"]



def simulate_diffusion_profile(source_scale: float = 1.0, diffusion_scale: float = 1.0, exposure_threshold: float = 5e-4) -> tuple[np.ndarray, dict[float, np.ndarray], float]:
    radius_cm = 0.1
    n_points = 100
    x = np.linspace(0.0, radius_cm, n_points)
    dx = x[1] - x[0]
    D = 1e-7 * diffusion_scale
    k_loss = (0.01 + 0.005) / 60.0
    dt = min(4.0, 0.45 * dx**2 / D)
    total_time_s = 24 * 3600
    n_steps = int(total_time_s / dt)
    concentration = np.zeros_like(x)
    snapshots: dict[float, np.ndarray] = {0.0: concentration.copy()}
    snapshot_hours = [1.0, 4.0, 12.0, 24.0]
    snapshot_seconds = {int(h * 3600): h for h in snapshot_hours}
    source_strength = 0.0025 * source_scale

    for step in range(1, n_steps + 1):
        t_s = int(step * dt)
        lap = np.zeros_like(concentration)
        lap[1:-1] = (concentration[2:] - 2 * concentration[1:-1] + concentration[:-2]) / dx**2
        lap[0] = 2 * (concentration[1] - concentration[0]) / dx**2
        lap[-1] = 2 * (concentration[-2] - concentration[-1]) / dx**2
        source = np.zeros_like(concentration)
        source[0] = source_strength * math.exp(-t_s / (6 * 3600))
        concentration = concentration + dt * (D * lap - k_loss * concentration + source)
        concentration = np.clip(concentration, 0.0, None)
        if t_s in snapshot_seconds:
            snapshots[snapshot_seconds[t_s]] = concentration.copy()

    final_profile = snapshots.get(24.0, concentration.copy())
    threshold = exposure_threshold
    bystander_radius_mm = 0.0
    above = np.where(final_profile >= threshold)[0]
    if len(above):
        bystander_radius_mm = float(x[above[-1]] * 10.0)
    return x, snapshots, bystander_radius_mm



def section_bystander_diffusion() -> dict:
    x, snapshots, bystander_radius_mm = simulate_diffusion_profile(source_scale=1.0, diffusion_scale=1.0)
    records = []
    for hour, profile in snapshots.items():
        for position, conc in zip(x, profile, strict=False):
            records.append(
                {
                    "time_h": hour,
                    "position_mm": position * 10.0,
                    "free_drug_concentration_au": conc,
                }
            )
    diffusion_df = pd.DataFrame(records)
    save_csv(diffusion_df, DATA_DIR / "bystander_diffusion_profiles.csv")

    fig, ax = plt.subplots(figsize=(10, 7))
    for hour in [0.0, 1.0, 4.0, 12.0, 24.0]:
        profile = snapshots.get(hour)
        if profile is not None:
            ax.plot(x * 10.0, profile, linewidth=2, label=f"t={hour:.0f} h")
    ax.set_xlabel("Tumor radial position (mm)")
    ax.set_ylabel("Free drug concentration (a.u.)")
    ax.set_title("Bystander diffusion profiles in tumor tissue")
    ax.legend()
    save_figure(fig, FIG_DIR / "03_bystander_diffusion.png")

    summary_metrics["bystander"] = {
        "radius_mm_at_24h": bystander_radius_mm,
        "center_concentration_24h": float(snapshots[24.0][0]),
        "edge_concentration_24h": float(snapshots[24.0][-1]),
    }
    return summary_metrics["bystander"]



def section_optimization() -> dict:
    k_plasma_vals = np.logspace(-3, -1, 55)
    k_tumor_vals = np.logspace(-1, 1, 80)
    horizon_days = 1.0
    landscape = []
    for k_p in k_plasma_vals:
        for k_t in k_tumor_vals:
            efficacy = 1.0 - math.exp(-k_t * horizon_days)
            toxicity = 1.0 - math.exp(-k_p * horizon_days)
            ratio_ok = (k_t / k_p) >= 50.0
            objective = efficacy - toxicity if ratio_ok else np.nan
            landscape.append(
                {
                    "k_cleavage_plasma_per_day": k_p,
                    "k_cleavage_tumor_per_day": k_t,
                    "efficacy_score": efficacy,
                    "toxicity_score": toxicity,
                    "objective": objective,
                    "selectivity_ratio": k_t / k_p,
                    "is_feasible": ratio_ok,
                }
            )
    landscape_df = pd.DataFrame(landscape)
    save_csv(landscape_df, DATA_DIR / "optimization_landscape.csv")

    feasible_df = landscape_df[landscape_df["is_feasible"]].copy()
    feasible_df = feasible_df.sort_values(["toxicity_score", "efficacy_score"], ascending=[True, False])
    pareto_rows = []
    best_efficacy = -np.inf
    for row in feasible_df.itertuples(index=False):
        if row.efficacy_score > best_efficacy:
            pareto_rows.append(row)
            best_efficacy = row.efficacy_score
    pareto_df = pd.DataFrame(pareto_rows)
    save_csv(pareto_df, RESULTS_DIR / "optimization_results.csv")

    best_row = feasible_df.loc[feasible_df["objective"].idxmax()]

    bounds = [(k_plasma_vals.min(), k_plasma_vals.max()), (max(50 * k_plasma_vals.min(), k_tumor_vals.min()), k_tumor_vals.max())]

    def objective_to_minimize(x):
        k_p, k_t = x
        efficacy = 1.0 - math.exp(-k_t * horizon_days)
        toxicity = 1.0 - math.exp(-k_p * horizon_days)
        penalty = max(0.0, 50.0 - (k_t / max(k_p, 1e-9))) * 10.0
        return -(efficacy - toxicity) + penalty

    opt = minimize(objective_to_minimize, x0=[best_row["k_cleavage_plasma_per_day"], best_row["k_cleavage_tumor_per_day"]], bounds=bounds, method="L-BFGS-B")

    grid = landscape_df.pivot(index="k_cleavage_tumor_per_day", columns="k_cleavage_plasma_per_day", values="objective")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sns.heatmap(grid, cmap="viridis", ax=axes[0], cbar_kws={"label": "Objective J"})
    axes[0].set_title("Plasma stability vs tumor release landscape")
    axes[0].set_xlabel("k_cleavage_plasma (grid index, log-scaled values)")
    axes[0].set_ylabel("k_cleavage_tumor (grid index, log-scaled values)")

    axes[1].scatter(feasible_df["toxicity_score"], feasible_df["efficacy_score"], s=18, alpha=0.3, color=sns.color_palette("viridis", 6)[2], label="Feasible grid")
    axes[1].plot(pareto_df["toxicity_score"], pareto_df["efficacy_score"], color="crimson", linewidth=2, label="Pareto front")
    axes[1].scatter(best_row["toxicity_score"], best_row["efficacy_score"], color="black", s=60, label="Best grid point")
    axes[1].set_xlabel("Toxicity score")
    axes[1].set_ylabel("Efficacy score")
    axes[1].set_title("Pareto front")
    axes[1].legend()
    save_figure(fig, FIG_DIR / "04_optimization_landscape.png")

    summary_metrics["optimization"] = {
        "best_grid_k_plasma": float(best_row["k_cleavage_plasma_per_day"]),
        "best_grid_k_tumor": float(best_row["k_cleavage_tumor_per_day"]),
        "best_objective": float(best_row["objective"]),
        "minimize_solution_fun": float(-opt.fun),
    }
    return summary_metrics["optimization"]


BASE_PKPD_PARAMS = {
    "CL": 0.5,
    "Vc": 3.0,
    "Vp": 6.0,
    "Q": 1.0,
    "Q_drug": 1.8,
    "CL_drug": 6.0,
    "CL_drug_p": 0.4,
    "k_release_plasma": 0.02,
    "k_release_tumor": 0.5,
    "k_internalization": 0.03,
    "k_internalization_tumor": 0.25,
    "k_distribution": 0.08,
    "k_dist_p": 0.05,
    "k_clearance_tumor": 0.18,
    "kon": 1.0,
    "koff": 0.1,
    "ksyn": 0.15,
    "kdeg": 0.05,
    "kdeg_complex": 0.20,
    "Emax": 0.8,
    "EC50": 0.1,
    "hill": 2.0,
    "kg": 0.1,
    "Kmax": 1.5,
    "dose_mg_per_kg": 6.4,
    "body_weight_kg": 70.0,
    "DAR": 8.0,
}



def initial_state_from_params(params: dict) -> np.ndarray:
    dose_mg = params["dose_mg_per_kg"] * params["body_weight_kg"]
    adc0 = dose_mg / params["Vc"]
    return np.array([adc0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0], dtype=float)



def enrich_params(params: dict) -> dict:
    p = deepcopy(params)
    dar_factor = p.get("DAR", 8.0) / 8.0
    p["k_release_plasma_eff"] = p["k_release_plasma"] * dar_factor
    p["k_release_tumor_eff"] = p["k_release_tumor"] * dar_factor
    return p



def pkpd_rhs(t: float, y: np.ndarray, params: dict) -> np.ndarray:
    adc_c, drug_c, adc_p, drug_p, adc_t, drug_t, target, drug_target, cell = y
    p = params
    drug_t_nonneg = max(drug_t, 0.0)
    target_nonneg = max(target, 0.0)
    cell_nonneg = max(cell, 0.0)

    d_adc_c = -(p["CL"] / p["Vc"]) * adc_c - (p["Q"] / p["Vc"]) * (adc_c - adc_p) - p["k_internalization"] * adc_c
    d_drug_c = p["k_release_plasma_eff"] * adc_c - (p["CL_drug"] / p["Vc"]) * drug_c - (p["Q_drug"] / p["Vc"]) * (drug_c - drug_p)
    d_adc_p = (p["Q"] / p["Vp"]) * (adc_c - adc_p) - p["k_dist_p"] * adc_p
    d_drug_p = (p["Q_drug"] / p["Vp"]) * (drug_c - drug_p) - p["CL_drug_p"] * drug_p
    d_adc_t = p["k_distribution"] * adc_c - p["k_internalization_tumor"] * adc_t
    d_drug_t = (
        p["k_release_tumor_eff"] * adc_t
        - p["k_clearance_tumor"] * drug_t
        - p["kon"] * drug_t_nonneg * target_nonneg
        + p["koff"] * drug_target
    )
    d_target = p["ksyn"] - p["kdeg"] * target_nonneg - p["kon"] * drug_t_nonneg * target_nonneg + p["koff"] * drug_target
    d_drug_target = p["kon"] * drug_t_nonneg * target_nonneg - p["koff"] * drug_target - p["kdeg_complex"] * drug_target
    effect = p["Emax"] * (drug_t_nonneg ** p["hill"]) / (p["EC50"] ** p["hill"] + drug_t_nonneg ** p["hill"] + 1e-12)
    d_cell = p["kg"] * cell_nonneg * (1.0 - cell_nonneg / p["Kmax"]) - effect * cell_nonneg
    return np.array([d_adc_c, d_drug_c, d_adc_p, d_drug_p, d_adc_t, d_drug_t, d_target, d_drug_target, d_cell], dtype=float)



def simulate_pkpd(params: dict, t_end: float = 21.0, n_points: int = 421) -> pd.DataFrame:
    p = enrich_params(params)
    y0 = initial_state_from_params(p)
    t_eval = np.linspace(0.0, t_end, n_points)
    try:
        sol = solve_ivp(lambda t, y: pkpd_rhs(t, y, p), (0.0, t_end), y0, t_eval=t_eval, method="LSODA")
    except Exception:
        sol = solve_ivp(lambda t, y: pkpd_rhs(t, y, p), (0.0, t_end), y0, t_eval=t_eval, method="RK45")
    if not sol.success:
        raise RuntimeError(sol.message)
    data = pd.DataFrame(
        {
            "time_day": sol.t,
            "ADC_central": sol.y[0],
            "Drug_central": sol.y[1],
            "ADC_peripheral": sol.y[2],
            "Drug_peripheral": sol.y[3],
            "ADC_tumor": sol.y[4],
            "Drug_tumor": sol.y[5],
            "Target": sol.y[6],
            "Drug_Target": sol.y[7],
            "Cell_viable": sol.y[8],
        }
    )
    return data.clip(lower=0.0)



def simulate_pkpd_euler(params: dict, t_end: float = 21.0, dt: float = 0.1) -> pd.DataFrame:
    p = enrich_params(params)
    times = np.arange(0.0, t_end + dt, dt)
    y = initial_state_from_params(p)
    records = []
    for t in times:
        records.append([t, *y.tolist()])
        y = y + dt * pkpd_rhs(t, y, p)
        y = np.clip(y, 0.0, None)
    cols = [
        "time_day",
        "ADC_central",
        "Drug_central",
        "ADC_peripheral",
        "Drug_peripheral",
        "ADC_tumor",
        "Drug_tumor",
        "Target",
        "Drug_Target",
        "Cell_viable",
    ]
    return pd.DataFrame(records, columns=cols)



def summarize_pkpd(df: pd.DataFrame) -> dict:
    return {
        "plasma_adc_auc": float(trapezoid(df["ADC_central"], df["time_day"])),
        "plasma_drug_auc": float(trapezoid(df["Drug_central"], df["time_day"])),
        "tumor_drug_auc": float(trapezoid(df["Drug_tumor"], df["time_day"])),
        "day21_cell_fraction": float(df["Cell_viable"].iloc[-1]),
    }



def section_pkpd() -> dict:
    timecourse = simulate_pkpd(BASE_PKPD_PARAMS)
    save_csv(timecourse, RESULTS_DIR / "pkpd_timecourse.csv")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes = axes.ravel()
    axes[0].plot(timecourse["time_day"], timecourse["ADC_central"], linewidth=2)
    axes[0].set_title("ADC plasma PK")
    axes[0].set_xlabel("Time (day)")
    axes[0].set_ylabel("ADC concentration (ug/mL)")

    axes[1].plot(timecourse["time_day"], timecourse["Drug_central"], linewidth=2)
    axes[1].set_title("Free drug in plasma")
    axes[1].set_xlabel("Time (day)")
    axes[1].set_ylabel("Drug concentration (ug/mL)")

    axes[2].plot(timecourse["time_day"], timecourse["Drug_tumor"], linewidth=2)
    axes[2].set_title("Tumor free drug exposure")
    axes[2].set_xlabel("Time (day)")
    axes[2].set_ylabel("Drug concentration (ug/mL)")

    axes[3].plot(timecourse["time_day"], timecourse["Cell_viable"], linewidth=2)
    axes[3].set_title("Viable tumor cell fraction")
    axes[3].set_xlabel("Time (day)")
    axes[3].set_ylabel("Cell fraction")
    save_figure(fig, FIG_DIR / "05_pkpd_simulation.png")

    summary_metrics["pkpd"] = summarize_pkpd(timecourse)
    return summary_metrics["pkpd"]



def section_monte_carlo() -> dict:
    sampler = qmc.LatinHypercube(d=5, seed=SEED)
    unit_samples = sampler.random(1000)
    ranges = {
        "CL": (BASE_PKPD_PARAMS["CL"] * 0.7, BASE_PKPD_PARAMS["CL"] * 1.3),
        "Vc": (BASE_PKPD_PARAMS["Vc"] * 0.7, BASE_PKPD_PARAMS["Vc"] * 1.3),
        "k_release_tumor": (BASE_PKPD_PARAMS["k_release_tumor"] * 0.7, BASE_PKPD_PARAMS["k_release_tumor"] * 1.3),
        "EC50": (BASE_PKPD_PARAMS["EC50"] * 0.7, BASE_PKPD_PARAMS["EC50"] * 1.3),
        "DAR": (BASE_PKPD_PARAMS["DAR"] * 0.7, BASE_PKPD_PARAMS["DAR"] * 1.3),
    }
    params_names = list(ranges.keys())
    lower = np.array([ranges[k][0] for k in params_names])
    upper = np.array([ranges[k][1] for k in params_names])
    scaled = qmc.scale(unit_samples, lower, upper)

    outputs = []
    for idx, values in enumerate(scaled, start=1):
        if idx % 100 == 0:
            print(f"[INFO] Monte Carlo progress: {idx}/1000")
        sample_params = deepcopy(BASE_PKPD_PARAMS)
        for name, value in zip(params_names, values, strict=False):
            sample_params[name] = float(value)
        try:
            df = simulate_pkpd_euler(sample_params, t_end=21.0, dt=0.1)
            tumor_auc = float(trapezoid(df["Drug_tumor"], df["time_day"]))
            cell_fraction = float(df["Cell_viable"].iloc[-1])
            cell_kill = float(max(0.0, 1.0 - cell_fraction))
            outputs.append(
                {
                    **{name: float(value) for name, value in zip(params_names, values, strict=False)},
                    "tumor_auc": tumor_auc,
                    "cell_fraction_day21": cell_fraction,
                    "cell_kill_day21": cell_kill,
                }
            )
        except Exception as exc:  # noqa: BLE001
            outputs.append(
                {
                    **{name: float(value) for name, value in zip(params_names, values, strict=False)},
                    "tumor_auc": np.nan,
                    "cell_fraction_day21": np.nan,
                    "cell_kill_day21": np.nan,
                    "error": str(exc),
                }
            )

    mc_df = pd.DataFrame(outputs)
    save_csv(mc_df, RESULTS_DIR / "monte_carlo_results.csv")

    clean_df = mc_df.dropna(subset=["tumor_auc", "cell_fraction_day21", "cell_kill_day21"]).copy()
    sensitivity_rows = []
    for pname in params_names:
        rho_auc, _ = spearmanr(clean_df[pname], clean_df["tumor_auc"])
        rho_kill, _ = spearmanr(clean_df[pname], clean_df["cell_kill_day21"])
        sensitivity_rows.append({"parameter": pname, "rho_tumor_auc": rho_auc, "rho_cell_kill": rho_kill})
    sensitivity_df = pd.DataFrame(sensitivity_rows).sort_values("rho_cell_kill", key=np.abs, ascending=True)
    save_csv(sensitivity_df, DATA_DIR / "monte_carlo_sensitivity_coefficients.csv")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    y_pos = np.arange(len(sensitivity_df))
    axes[0].barh(y_pos - 0.18, sensitivity_df["rho_tumor_auc"], height=0.35, color=sns.color_palette("viridis", 6)[2], label="Tumor AUC")
    axes[0].barh(y_pos + 0.18, sensitivity_df["rho_cell_kill"], height=0.35, color=sns.color_palette("viridis", 6)[4], label="Cell kill day 21")
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels(sensitivity_df["parameter"])
    axes[0].set_xlabel("Spearman correlation coefficient")
    axes[0].set_title("Monte Carlo sensitivity tornado chart")
    axes[0].legend()

    axes[1].hist(clean_df["cell_fraction_day21"], bins=30, color=sns.color_palette("viridis", 6)[3], edgecolor="black")
    axes[1].set_xlabel("Tumor viable cell fraction at day 21")
    axes[1].set_ylabel("Count")
    axes[1].set_title("Distribution of Monte Carlo outcomes")
    save_figure(fig, FIG_DIR / "06_monte_carlo_sensitivity.png")

    ci_low, ci_high = bootstrap_ci(clean_df["cell_fraction_day21"].to_numpy())
    summary_metrics["monte_carlo"] = {
        "n_success": int(len(clean_df)),
        "mean_tumor_auc": float(clean_df["tumor_auc"].mean()),
        "mean_cell_fraction_day21": float(clean_df["cell_fraction_day21"].mean()),
        "cell_fraction_day21_ci95": [float(ci_low), float(ci_high)],
    }

    stats_md = f"""# Statistical Summary\n\n- Random seed: {SEED}\n- Monte Carlo successful runs: {len(clean_df)} / 1000\n- Mean tumor AUC: {clean_df['tumor_auc'].mean():.3f}\n- Mean viable cell fraction at day 21: {clean_df['cell_fraction_day21'].mean():.3f}\n- 95% bootstrap CI for mean viable cell fraction: [{ci_low:.3f}, {ci_high:.3f}]\n\n## Sensitivity (Spearman correlation)\n\n{sensitivity_df.to_markdown(index=False)}\n"""
    save_text(RESULTS_DIR / "statistical-summary.md", stats_md)
    return summary_metrics["monte_carlo"]



def scenario_params(dar: float, cleavable: bool) -> dict:
    params = deepcopy(BASE_PKPD_PARAMS)
    params["DAR"] = dar
    if cleavable:
        params["k_release_tumor"] = 0.25 + 0.25 * (dar / 8.0)
        params["k_release_plasma"] = 0.02
    else:
        params["k_release_tumor"] = 0.08
        params["k_release_plasma"] = 0.005
    kd_nM = 0.1
    params["kon"] = 1.0
    params["koff"] = kd_nM * params["kon"]
    return params



def section_case_study() -> dict:
    scenarios = {
        "DAR4_cleavable": {"dar": 4.0, "cleavable": True, "label": "DAR=4 cleavable"},
        "DAR8_cleavable": {"dar": 8.0, "cleavable": True, "label": "DAR=8 cleavable"},
        "DAR8_noncleavable": {"dar": 8.0, "cleavable": False, "label": "DAR=8 non-cleavable"},
    }
    rows = []
    for key, config in scenarios.items():
        params = scenario_params(config["dar"], config["cleavable"])
        df = simulate_pkpd(params)
        pk_summary = summarize_pkpd(df)
        source_scale = (config["dar"] / 8.0) * (1.0 if config["cleavable"] else 0.35)
        _, _, bystander_radius = simulate_diffusion_profile(source_scale=source_scale, diffusion_scale=1.0)
        dar_ti = float(therapeutic_index_profile(np.array([config["dar"]]))[0])
        therapeutic_index = dar_ti * safe_ratio(pk_summary["tumor_drug_auc"], pk_summary["plasma_drug_auc"] + 1e-9)
        rows.append(
            {
                "scenario": config["label"],
                "DAR": config["dar"],
                "cleavable_linker": config["cleavable"],
                "therapeutic_index": therapeutic_index,
                "tumor_auc": pk_summary["tumor_drug_auc"],
                "plasma_drug_exposure": pk_summary["plasma_drug_auc"],
                "bystander_radius_mm": bystander_radius,
                "day21_viable_cell_fraction": pk_summary["day21_cell_fraction"],
            }
        )
    case_df = pd.DataFrame(rows)
    save_csv(case_df, RESULTS_DIR / "case_study_summary.csv")

    metrics = ["therapeutic_index", "tumor_auc", "plasma_drug_exposure", "bystander_radius_mm"]
    normalized = case_df.copy()
    for metric in metrics:
        max_val = normalized[metric].max()
        normalized[metric] = 100.0 * normalized[metric] / max_val if max_val > 0 else 0.0

    fig = plt.figure(figsize=(14, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0])
    ax_bar = fig.add_subplot(gs[0, 0])
    plot_df = normalized.melt(id_vars="scenario", value_vars=metrics, var_name="metric", value_name="normalized_value")
    sns.barplot(data=plot_df, x="metric", y="normalized_value", hue="scenario", palette="viridis", ax=ax_bar)
    ax_bar.set_ylabel("Normalized performance (%)")
    ax_bar.set_xlabel("Metric")
    ax_bar.set_title("HER2 ADC scenario comparison")
    ax_bar.tick_params(axis="x", rotation=20)

    ax_radar = fig.add_subplot(gs[0, 1], polar=True)
    radar_metrics = metrics
    angles = np.linspace(0, 2 * np.pi, len(radar_metrics), endpoint=False).tolist()
    angles += angles[:1]
    for _, row in normalized.iterrows():
        values = [row[m] for m in radar_metrics]
        values += values[:1]
        ax_radar.plot(angles, values, linewidth=2, label=row["scenario"])
        ax_radar.fill(angles, values, alpha=0.10)
    ax_radar.set_xticks(angles[:-1])
    ax_radar.set_xticklabels(radar_metrics)
    ax_radar.set_yticklabels([])
    ax_radar.set_title("Multi-parameter radar profile")
    ax_radar.legend(loc="upper right", bbox_to_anchor=(1.25, 1.15))
    save_figure(fig, FIG_DIR / "07_case_study_comparison.png")

    best_scenario = case_df.loc[case_df["therapeutic_index"].idxmax(), "scenario"]
    summary_metrics["case_study"] = {
        "best_scenario_by_ti": str(best_scenario),
        "max_therapeutic_index": float(case_df["therapeutic_index"].max()),
        "max_bystander_radius_mm": float(case_df["bystander_radius_mm"].max()),
    }
    return summary_metrics["case_study"]



def write_preprocessing_log() -> None:
    text = f"""# Preprocessing Log\n\n- Random seed fixed at {SEED} for numpy and random.\n- All outputs were generated from deterministic simulation settings unless Monte Carlo sampling was explicitly requested.\n- DAR sampling used binomial and Poisson models with support truncated to DAR 0-8.\n- Linker kinetics were normalized to released fraction for cross-mechanism comparison.\n- Reaction-diffusion PDE used a 1D finite-difference explicit solver with no-flux boundaries.\n- Optimization landscape stored complete grid in `data/optimization_landscape.csv`; Pareto front stored in `results/optimization_results.csv`.\n- PK/PD nominal dose assumed 6.4 mg/kg IV for a 70 kg patient.\n- Monte Carlo sensitivity used Latin Hypercube Sampling over ±30% parameter ranges.\n- Case study metrics were normalized only for visualization; raw values remain in `results/case_study_summary.csv`.\n"""
    save_text(DATA_DIR / "preprocessing-log.md", text)


results = {
    "dar": run_section("DAR distribution", section_dar_distribution),
    "linker": run_section("Linker cleavage", section_linker_cleavage),
    "bystander": run_section("Bystander diffusion", section_bystander_diffusion),
    "optimization": run_section("Optimization", section_optimization),
    "pkpd": run_section("PK/PD", section_pkpd),
    "monte_carlo": run_section("Monte Carlo sensitivity", section_monte_carlo),
    "case_study": run_section("Case study", section_case_study),
}

write_preprocessing_log()

summary_metrics["failures"] = {"messages": section_failures}
SUMMARY_PATH.write_text(json.dumps(summary_metrics, indent=2, ensure_ascii=False), encoding="utf-8")
LOGGER.log(
    phase="REPORT",
    event_type="file_written",
    skill_or_tool="json.dump",
    files_written=[str(SUMMARY_PATH.relative_to(BASE_DIR))],
    handoff_out={"sections": list(summary_metrics.keys())},
)

LOGGER.log(
    phase="REPORT",
    event_type="report_finalized",
    skill_or_tool="adc_platform.py",
    handoff_out={"summary_file": str(SUMMARY_PATH.relative_to(BASE_DIR)), "failures": section_failures},
    status="ok" if not section_failures else "warning",
)
LOGGER.log(
    phase="LOG",
    event_type="run_completed",
    skill_or_tool="adc_platform.py",
    handoff_out={"completed_sections": [k for k, v in results.items() if v is not None], "failed_sections": section_failures},
    status="ok" if not section_failures else "warning",
)

print("[INFO] ADC platform execution finished.")
if section_failures:
    print("[WARN] Some sections encountered issues:")
    for msg in section_failures:
        print(f" - {msg}")
