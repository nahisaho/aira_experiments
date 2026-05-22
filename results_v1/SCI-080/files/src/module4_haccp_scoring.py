from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
plt.style.use("seaborn-v0_8-whitegrid")

BASE_DIR = Path(__file__).resolve().parents[1]
FIGURES_DIR = BASE_DIR / "figures"
RESULTS_DIR = BASE_DIR / "results"
DATA_DIR = BASE_DIR / "data"
LOGS_DIR = BASE_DIR / "logs"
REPORT_PATH = BASE_DIR / "report.md"
PREPROCESS_LOG_PATH = DATA_DIR / "preprocessing-log.md"
STAT_SUMMARY_PATH = RESULTS_DIR / "statistical-summary.md"
METRICS_PATH = RESULTS_DIR / "module4_metrics.json"
DATA_PATH = DATA_DIR / "haccp_monitoring.csv"
PROCESS_LOG_PATH = LOGS_DIR / "process-log.jsonl"
DPI = 300

for directory in (FIGURES_DIR, RESULTS_DIR, DATA_DIR, LOGS_DIR):
    directory.mkdir(parents=True, exist_ok=True)


@dataclass(frozen=True)
class CCPConfig:
    ccp_id: str
    name: str
    parameter: str
    unit: str
    hazard: str
    readings_per_day: int
    direction: str
    target: float
    warning_low: float | None
    warning_high: float | None
    critical_low: float | None
    critical_high: float | None
    severity: int
    detectability: int
    normal_mean: float
    normal_sd: float
    seasonal_amplitude: float
    drift_after_day: int
    drift_per_day: float
    base_near_miss: float
    base_deviation: float
    corrective_action: str


CCP_CONFIGS = [
    CCPConfig(
        ccp_id="CCP1",
        name="Receiving",
        parameter="Product temperature",
        unit="°C",
        hazard="Pathogen growth during receiving",
        readings_per_day=4,
        direction="max",
        target=2.0,
        warning_low=None,
        warning_high=4.0,
        critical_low=None,
        critical_high=5.0,
        severity=8,
        detectability=3,
        normal_mean=2.3,
        normal_sd=0.6,
        seasonal_amplitude=0.7,
        drift_after_day=260,
        drift_per_day=0.002,
        base_near_miss=0.030,
        base_deviation=0.010,
        corrective_action="Reject or hold the lot, verify truck refrigeration records, and intensify receiving inspections.",
    ),
    CCPConfig(
        ccp_id="CCP2",
        name="Cold storage",
        parameter="Chiller temperature",
        unit="°C",
        hazard="Pathogen growth during chilled storage",
        readings_per_day=8,
        direction="max",
        target=2.0,
        warning_low=None,
        warning_high=3.5,
        critical_low=None,
        critical_high=4.0,
        severity=9,
        detectability=4,
        normal_mean=2.1,
        normal_sd=0.5,
        seasonal_amplitude=0.8,
        drift_after_day=210,
        drift_per_day=0.003,
        base_near_miss=0.035,
        base_deviation=0.012,
        corrective_action="Quarantine exposed product, restore refrigeration, and review door opening frequency and alarm response.",
    ),
    CCPConfig(
        ccp_id="CCP3",
        name="Cooking/thermal processing",
        parameter="Core temperature",
        unit="°C",
        hazard="Survival of Salmonella and Campylobacter",
        readings_per_day=6,
        direction="min",
        target=77.0,
        warning_low=75.0,
        warning_high=None,
        critical_low=74.0,
        critical_high=None,
        severity=10,
        detectability=3,
        normal_mean=77.8,
        normal_sd=0.9,
        seasonal_amplitude=0.2,
        drift_after_day=240,
        drift_per_day=-0.0015,
        base_near_miss=0.025,
        base_deviation=0.009,
        corrective_action="Stop the line, reprocess affected lots, verify cook settings, and recalibrate temperature probes.",
    ),
    CCPConfig(
        ccp_id="CCP4",
        name="Cooling",
        parameter="Hours to reach ≤5°C",
        unit="hours",
        hazard="Clostridium perfringens growth during cooling",
        readings_per_day=4,
        direction="max",
        target=4.0,
        warning_low=None,
        warning_high=5.0,
        critical_low=None,
        critical_high=6.0,
        severity=9,
        detectability=5,
        normal_mean=4.2,
        normal_sd=0.55,
        seasonal_amplitude=0.8,
        drift_after_day=220,
        drift_per_day=0.003,
        base_near_miss=0.040,
        base_deviation=0.015,
        corrective_action="Segregate product, accelerate cooling, inspect airflow/load depth, and document disposition decisions.",
    ),
    CCPConfig(
        ccp_id="CCP5",
        name="Metal detection",
        parameter="Detector challenge response",
        unit="%",
        hazard="Physical contamination from ferrous/non-ferrous metal",
        readings_per_day=3,
        direction="min",
        target=99.0,
        warning_low=94.0,
        warning_high=None,
        critical_low=90.0,
        critical_high=None,
        severity=9,
        detectability=2,
        normal_mean=98.6,
        normal_sd=0.9,
        seasonal_amplitude=0.1,
        drift_after_day=250,
        drift_per_day=-0.003,
        base_near_miss=0.018,
        base_deviation=0.007,
        corrective_action="Place product on hold, recalibrate detector, retest challenge packs, and inspect upstream equipment wear.",
    ),
    CCPConfig(
        ccp_id="CCP6",
        name="Packaging integrity",
        parameter="Seal strength",
        unit="N",
        hazard="Post-lethality contamination from package failure",
        readings_per_day=5,
        direction="min",
        target=20.0,
        warning_low=17.0,
        warning_high=None,
        critical_low=15.0,
        critical_high=None,
        severity=7,
        detectability=6,
        normal_mean=20.4,
        normal_sd=1.2,
        seasonal_amplitude=0.4,
        drift_after_day=200,
        drift_per_day=-0.0035,
        base_near_miss=0.035,
        base_deviation=0.013,
        corrective_action="Stop packaging, inspect seal jaws and film, rework or discard compromised packs, and increase leak checks.",
    ),
    CCPConfig(
        ccp_id="CCP7",
        name="Shipping/distribution",
        parameter="Trailer product temperature",
        unit="°C",
        hazard="Cold chain failure during distribution",
        readings_per_day=4,
        direction="max",
        target=2.5,
        warning_low=None,
        warning_high=3.5,
        critical_low=None,
        critical_high=4.0,
        severity=8,
        detectability=5,
        normal_mean=2.7,
        normal_sd=0.6,
        seasonal_amplitude=0.9,
        drift_after_day=230,
        drift_per_day=0.0025,
        base_near_miss=0.032,
        base_deviation=0.012,
        corrective_action="Hold shipment, correct reefer setpoint, review route dwell times, and verify customer receiving temperature.",
    ),
]


def log_event(event_type: str, skill_or_tool: str, handoff_in: dict[str, Any], handoff_out: dict[str, Any], files_written: list[str], status: str = "ok", phase: str = "module4") -> None:
    payload = {
        "timestamp": datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill_or_tool,
        "handoff_in": handoff_in,
        "handoff_out": handoff_out,
        "files_written": files_written,
        "status": status,
    }
    with PROCESS_LOG_PATH.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload) + "\n")


def safe_float(value: float | None) -> float | None:
    if value is None:
        return None
    return float(value)


def classify_state(cfg: CCPConfig, value: float) -> tuple[str, int, int]:
    if cfg.direction == "max":
        if value > float(cfg.critical_high):
            return "deviation", 1, 0
        if value >= float(cfg.warning_high):
            return "near_miss", 0, 1
        return "normal", 0, 0

    if value < float(cfg.critical_low):
        return "deviation", 1, 0
    if value <= float(cfg.warning_low):
        return "near_miss", 0, 1
    return "normal", 0, 0


def generate_value(cfg: CCPConfig, selected_state: str, baseline: float, rng: np.random.Generator) -> float:
    if cfg.direction == "max":
        if selected_state == "normal":
            value = rng.normal(baseline, cfg.normal_sd)
            return min(value, float(cfg.warning_high) - 0.05)
        if selected_state == "near_miss":
            low = float(cfg.warning_high) + 0.02
            high = float(cfg.critical_high) - 0.02
            return float(np.clip(rng.normal((low + high) / 2, 0.12), low, high))
        return float(cfg.critical_high + abs(rng.normal(0.5, 0.35)))

    if selected_state == "normal":
        value = rng.normal(baseline, cfg.normal_sd)
        return max(value, float(cfg.warning_low) + 0.05)
    if selected_state == "near_miss":
        low = float(cfg.critical_low) + 0.05
        high = float(cfg.warning_low) - 0.02
        return float(np.clip(rng.normal((low + high) / 2, 0.18), low, high))
    return float(cfg.critical_low - abs(rng.normal(1.0, 0.45)))


def generate_dataset(configs: list[CCPConfig]) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    all_rows: list[dict[str, Any]] = []
    dates = pd.date_range("2024-01-01", periods=365, freq="D")

    for cfg in configs:
        time_slots = np.linspace(6, 21, cfg.readings_per_day)
        for day_idx, day in enumerate(dates):
            season_signal = math.sin(2 * math.pi * day_idx / 365)
            drift_days = max(day_idx - cfg.drift_after_day, 0)
            drift = drift_days * cfg.drift_per_day
            multiplier = 1 + max(season_signal, 0) * 0.8 + (0.35 if drift_days > 0 else 0.0)
            near_prob = min(cfg.base_near_miss * multiplier, 0.14)
            deviation_prob = min(cfg.base_deviation * multiplier, 0.08)
            normal_prob = max(1 - near_prob - deviation_prob, 0.75)
            probabilities = np.array([normal_prob, near_prob, deviation_prob])
            probabilities = probabilities / probabilities.sum()

            for slot_idx, hour in enumerate(time_slots):
                minute = (slot_idx * 13) % 60
                timestamp = day + pd.Timedelta(hours=float(hour), minutes=int(minute))
                selected_state = rng.choice(["normal", "near_miss", "deviation"], p=probabilities)
                baseline = cfg.normal_mean + cfg.seasonal_amplitude * max(season_signal, 0) + drift
                value = generate_value(cfg, selected_state, baseline, rng)
                state, deviation_flag, near_miss_flag = classify_state(cfg, value)
                action = (
                    "No action required; continue routine verification."
                    if state == "normal"
                    else "Increase monitoring frequency, inspect the process setting, and prevent escalation."
                    if state == "near_miss"
                    else cfg.corrective_action
                )
                all_rows.append(
                    {
                        "timestamp": timestamp,
                        "date": timestamp.date().isoformat(),
                        "ccp_id": cfg.ccp_id,
                        "ccp_name": cfg.name,
                        "hazard": cfg.hazard,
                        "parameter": cfg.parameter,
                        "unit": cfg.unit,
                        "measurement": round(value, 3),
                        "target_value": cfg.target,
                        "warning_low": safe_float(cfg.warning_low),
                        "warning_high": safe_float(cfg.warning_high),
                        "critical_limit_low": safe_float(cfg.critical_low),
                        "critical_limit_high": safe_float(cfg.critical_high),
                        "severity": cfg.severity,
                        "detectability": cfg.detectability,
                        "state": state,
                        "near_miss_flag": near_miss_flag,
                        "deviation_flag": deviation_flag,
                        "corrective_action": action,
                    }
                )

    return pd.DataFrame(all_rows).sort_values(["ccp_id", "timestamp"]).reset_index(drop=True)


def compute_trend_label(series: pd.Series) -> tuple[str, float]:
    tail = series.tail(min(60, len(series)))
    if len(tail) < 3:
        return "stable", 0.0
    x = np.arange(len(tail), dtype=float)
    slope = float(np.polyfit(x, tail.to_numpy(dtype=float), 1)[0])
    if slope > 0.03:
        return "increasing", slope
    if slope < -0.03:
        return "decreasing", slope
    return "stable", slope


def compute_ewma(series: pd.Series, lam: float = 0.2) -> tuple[pd.Series, pd.Series, pd.Series]:
    values = series.to_numpy(dtype=float)
    center = float(np.mean(values))
    sigma = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
    ewma = np.zeros_like(values)
    ewma[0] = center
    for idx in range(1, len(values)):
        ewma[idx] = lam * values[idx] + (1 - lam) * ewma[idx - 1]
    t = np.arange(1, len(values) + 1)
    factor = np.sqrt((lam / (2 - lam)) * (1 - (1 - lam) ** (2 * t)))
    limit = 3 * sigma * factor
    upper = center + limit
    lower = center - limit
    return pd.Series(ewma, index=series.index), pd.Series(upper, index=series.index), pd.Series(lower, index=series.index)


def compute_cusum(series: pd.Series, k: float = 0.5, h: float = 5.0) -> tuple[pd.Series, pd.Series, pd.Series]:
    values = series.to_numpy(dtype=float)
    center = float(np.mean(values))
    sigma = float(np.std(values, ddof=1)) if len(values) > 1 else 1.0
    sigma = sigma if sigma > 0 else 1.0
    z = (values - center) / sigma
    cpos = np.zeros_like(z)
    cneg = np.zeros_like(z)
    signals = np.zeros_like(z)
    for idx in range(1, len(z)):
        cpos[idx] = max(0.0, cpos[idx - 1] + z[idx] - k)
        cneg[idx] = min(0.0, cneg[idx - 1] + z[idx] + k)
        if cpos[idx] > h or abs(cneg[idx]) > h:
            signals[idx] = 1
            cpos[idx] = 0.0
            cneg[idx] = 0.0
    return pd.Series(cpos, index=series.index), pd.Series(cneg, index=series.index), pd.Series(signals, index=series.index)


def beta_interval(alpha: float, beta: float, seed_offset: int) -> tuple[float, float]:
    rng = np.random.default_rng(SEED + seed_offset)
    samples = rng.beta(alpha, beta, size=12000)
    lower, upper = np.quantile(samples, [0.025, 0.975])
    return float(lower), float(upper)


def bootstrap_mean_ci(values: np.ndarray, seed_offset: int) -> tuple[float, float]:
    rng = np.random.default_rng(SEED + seed_offset)
    boot = []
    for _ in range(3000):
        sample = rng.choice(values, size=len(values), replace=True)
        boot.append(sample.mean())
    lower, upper = np.quantile(boot, [0.025, 0.975])
    return float(lower), float(upper)


def compute_risk_metrics(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any], dict[str, dict[str, Any]]]:
    frames: list[pd.DataFrame] = []
    summary: dict[str, Any] = {}
    diagnostics: dict[str, dict[str, Any]] = {}
    prior_alpha, prior_beta = 2.0, 48.0

    for idx, cfg in enumerate(CCP_CONFIGS, start=1):
        group = df[df["ccp_id"] == cfg.ccp_id].copy().sort_values("timestamp")
        observations = np.arange(1, len(group) + 1)
        cumulative_dev = group["deviation_flag"].cumsum().astype(float)
        alpha = prior_alpha + cumulative_dev
        beta = prior_beta + observations - cumulative_dev
        posterior_prob = alpha / (alpha + beta)
        rolling_rate = group["deviation_flag"].rolling(window=50, min_periods=10).mean()
        rolling_rate = rolling_rate.fillna(group["deviation_flag"].expanding().mean()).clip(lower=0)
        dynamic_rate = 0.55 * rolling_rate + 0.45 * posterior_prob
        likelihood = 1 + 9 * np.clip(dynamic_rate / 0.08, 0, 1)
        bayes_likelihood = 1 + 9 * np.clip(posterior_prob / 0.08, 0, 1)
        group["posterior_deviation_probability"] = posterior_prob.round(6)
        group["dynamic_likelihood"] = likelihood.round(3)
        group["rpn"] = (cfg.severity * cfg.detectability * likelihood).round(3)
        group["bayesian_likelihood"] = bayes_likelihood.round(3)
        group["bayesian_rpn"] = (cfg.severity * cfg.detectability * bayes_likelihood).round(3)

        center = group["measurement"].mean()
        sigma = group["measurement"].std(ddof=1) if len(group) > 1 else 0.0
        ucl = center + 3 * sigma
        lcl = center - 3 * sigma
        group["spc_center"] = center
        group["spc_ucl"] = ucl
        group["spc_lcl"] = lcl
        group["spc_signal"] = ((group["measurement"] > ucl) | (group["measurement"] < lcl)).astype(int)

        ewma, ewma_ucl, ewma_lcl = compute_ewma(group["measurement"])
        cusum_pos, cusum_neg, cusum_signal = compute_cusum(group["measurement"])
        group["ewma"] = ewma.round(3)
        group["ewma_ucl"] = ewma_ucl.round(3)
        group["ewma_lcl"] = ewma_lcl.round(3)
        group["ewma_signal"] = ((ewma > ewma_ucl) | (ewma < ewma_lcl)).astype(int)
        group["cusum_pos"] = cusum_pos.round(3)
        group["cusum_neg"] = cusum_neg.round(3)
        group["cusum_signal"] = cusum_signal.astype(int)

        trend_label, trend_slope = compute_trend_label(group["bayesian_rpn"])
        ci_low, ci_high = beta_interval(float(alpha.iloc[-1]), float(beta.iloc[-1]), idx)
        first_half = group.iloc[: len(group) // 2]["bayesian_rpn"].to_numpy(dtype=float)
        second_half = group.iloc[len(group) // 2 :]["bayesian_rpn"].to_numpy(dtype=float)
        mean_delta = float(second_half.mean() - first_half.mean())
        pooled_sd = float(
            np.sqrt(
                max(
                    ((len(first_half) - 1) * np.var(first_half, ddof=1) + (len(second_half) - 1) * np.var(second_half, ddof=1))
                    / (len(first_half) + len(second_half) - 2),
                    1e-9,
                )
            )
        )
        cohens_d = mean_delta / pooled_sd
        pair_count = min(len(first_half), len(second_half))
        paired_delta = second_half[:pair_count] - first_half[:pair_count]
        delta_ci = bootstrap_mean_ci(paired_delta, idx + 100)

        summary[cfg.ccp_id] = {
            "ccp_name": cfg.name,
            "hazard": cfg.hazard,
            "records": int(len(group)),
            "near_misses": int(group["near_miss_flag"].sum()),
            "deviations": int(group["deviation_flag"].sum()),
            "deviation_rate": round(float(group["deviation_flag"].mean()), 4),
            "posterior_deviation_probability": round(float(posterior_prob.iloc[-1]), 4),
            "posterior_95pct_interval": [round(ci_low, 4), round(ci_high, 4)],
            "final_dynamic_likelihood": round(float(group["dynamic_likelihood"].iloc[-1]), 3),
            "mean_rpn": round(float(group["rpn"].mean()), 3),
            "final_bayesian_rpn": round(float(group["bayesian_rpn"].iloc[-1]), 3),
            "trend": trend_label,
            "trend_slope": round(trend_slope, 4),
            "spc_signals": int(group["spc_signal"].sum()),
            "ewma_signals": int(group["ewma_signal"].sum()),
            "cusum_signals": int(group["cusum_signal"].sum()),
            "cohens_d_second_vs_first_half": round(float(cohens_d), 3),
            "mean_change_ci": [round(delta_ci[0], 3), round(delta_ci[1], 3)],
        }
        diagnostics[cfg.ccp_id] = {
            "center": float(center),
            "ucl": float(ucl),
            "lcl": float(lcl),
        }
        frames.append(group)

    combined = pd.concat(frames, ignore_index=True)
    monthly = (
        combined.assign(month=pd.to_datetime(combined["timestamp"]).dt.to_period("M").dt.to_timestamp())
        .groupby(["ccp_id", "month"], as_index=False)[["bayesian_rpn", "deviation_flag"]]
        .mean()
    )
    summary["overall"] = {
        "records": int(len(combined)),
        "ccp_count": len(CCP_CONFIGS),
        "overall_deviation_rate": round(float(combined["deviation_flag"].mean()), 4),
        "overall_mean_bayesian_rpn": round(float(combined["bayesian_rpn"].mean()), 3),
        "months_observed": int(monthly["month"].nunique()),
        "highest_final_risk_ccp": max(
            [key for key in summary.keys() if key.startswith("CCP")],
            key=lambda key: summary[key]["final_bayesian_rpn"],
        ),
    }
    return combined, summary, diagnostics


def create_heatmap_figure(df: pd.DataFrame) -> Path:
    weekly = (
        df.assign(week=pd.to_datetime(df["timestamp"]).dt.to_period("W").dt.start_time)
        .groupby(["ccp_name", "week"], as_index=False)["bayesian_rpn"]
        .mean()
    )
    pivot = weekly.pivot(index="ccp_name", columns="week", values="bayesian_rpn").reindex([cfg.name for cfg in CCP_CONFIGS])
    fig, ax = plt.subplots(figsize=(15, 6), constrained_layout=True)
    image = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="cividis")
    ax.set_title("Weekly Bayesian RPN heatmap")
    ax.set_ylabel("Critical control point")
    ax.set_xlabel("Week")
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    tick_positions = np.arange(0, len(pivot.columns), 4)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([pd.Timestamp(pivot.columns[i]).strftime("%b") for i in tick_positions], rotation=0)
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("Bayesian RPN")
    path = FIGURES_DIR / "fig4_risk_scores.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def create_spc_figure(df: pd.DataFrame, diagnostics: dict[str, dict[str, Any]]) -> Path:
    key_ids = ["CCP1", "CCP3", "CCP4", "CCP7"]
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharex=False, constrained_layout=True)
    for ax, ccp_id in zip(axes.ravel(), key_ids):
        group = df[df["ccp_id"] == ccp_id].copy()
        cfg = next(item for item in CCP_CONFIGS if item.ccp_id == ccp_id)
        sample = group.iloc[:: max(len(group) // 180, 1)].copy()
        ax.plot(sample["timestamp"], sample["measurement"], color="#1b9e77", linewidth=1.0, label="Measurement")
        ax.plot(sample["timestamp"], sample["ewma"], color="#7570b3", linewidth=1.2, label="EWMA")
        ax.axhline(diagnostics[ccp_id]["center"], color="#4d4d4d", linestyle="--", linewidth=1, label="Center")
        ax.axhline(diagnostics[ccp_id]["ucl"], color="#d95f02", linestyle=":", linewidth=1.1, label="UCL/LCL")
        ax.axhline(diagnostics[ccp_id]["lcl"], color="#d95f02", linestyle=":", linewidth=1.1)
        if cfg.critical_high is not None:
            ax.axhline(cfg.critical_high, color="#e7298a", linestyle="-.", linewidth=1, label="Critical limit")
        if cfg.critical_low is not None:
            ax.axhline(cfg.critical_low, color="#e7298a", linestyle="-.", linewidth=1, label="Critical limit")
        ewma_alerts = sample[sample["ewma_signal"] == 1]
        cusum_alerts = sample[sample["cusum_signal"] == 1]
        if not ewma_alerts.empty:
            ax.scatter(ewma_alerts["timestamp"], ewma_alerts["measurement"], color="#66a61e", s=16, label="EWMA alert")
        if not cusum_alerts.empty:
            ax.scatter(cusum_alerts["timestamp"], cusum_alerts["measurement"], color="#e6ab02", marker="x", s=24, label="CUSUM alert")
        ax.set_title(cfg.name)
        ax.set_ylabel(f"{cfg.parameter} ({cfg.unit})")
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=5, frameon=False)
    path = FIGURES_DIR / "fig4b_spc_charts.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def create_rpn_distribution_figure(df: pd.DataFrame) -> Path:
    totals = (
        df.groupby("ccp_name", as_index=False)["bayesian_rpn"].sum().sort_values("bayesian_rpn", ascending=False)
    )
    totals["cumulative_pct"] = 100 * totals["bayesian_rpn"].cumsum() / totals["bayesian_rpn"].sum()
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    for cfg, color in zip(CCP_CONFIGS, plt.cm.viridis(np.linspace(0.1, 0.9, len(CCP_CONFIGS)))):
        subset = df[df["ccp_id"] == cfg.ccp_id]
        axes[0].hist(subset["bayesian_rpn"], bins=18, alpha=0.45, label=cfg.name, color=color)
    axes[0].set_title("Bayesian RPN distribution by CCP")
    axes[0].set_xlabel("Bayesian RPN")
    axes[0].set_ylabel("Frequency")
    axes[0].legend(fontsize=7, frameon=False)

    axes[1].bar(totals["ccp_name"], totals["bayesian_rpn"], color=plt.cm.cividis(np.linspace(0.15, 0.85, len(totals))))
    axes[1].set_ylabel("Cumulative risk burden")
    axes[1].set_title("Pareto ranking of CCP risk burden")
    axes[1].tick_params(axis="x", rotation=45)
    ax2 = axes[1].twinx()
    ax2.plot(totals["ccp_name"], totals["cumulative_pct"], color="#d95f02", marker="o")
    ax2.set_ylabel("Cumulative contribution (%)")
    ax2.set_ylim(0, 110)
    ax2.axhline(80, color="#7570b3", linestyle="--", linewidth=1)
    path = FIGURES_DIR / "fig4c_rpn_distribution.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def create_bayesian_update_figure(df: pd.DataFrame) -> Path:
    weekly = (
        df.assign(week=pd.to_datetime(df["timestamp"]).dt.to_period("W").dt.start_time)
        .groupby(["ccp_name", "week"], as_index=False)[["posterior_deviation_probability", "bayesian_rpn"]]
        .last()
    )
    fig, axes = plt.subplots(2, 1, figsize=(15, 9), sharex=True, constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(CCP_CONFIGS)))
    for cfg, color in zip(CCP_CONFIGS, colors):
        subset = weekly[weekly["ccp_name"] == cfg.name]
        axes[0].plot(subset["week"], subset["posterior_deviation_probability"], label=cfg.name, color=color, linewidth=1.8)
        axes[1].plot(subset["week"], subset["bayesian_rpn"], label=cfg.name, color=color, linewidth=1.8)
    axes[0].set_title("Bayesian deviation probability update")
    axes[0].set_ylabel("Posterior probability")
    axes[1].set_title("Bayesian RPN update over time")
    axes[1].set_ylabel("Bayesian RPN")
    axes[1].set_xlabel("Week")
    axes[1].xaxis.set_major_locator(mdates.MonthLocator(interval=1))
    axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    fig.legend(loc="upper center", ncol=4, frameon=False)
    path = FIGURES_DIR / "fig4d_bayesian_update.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def write_preprocessing_log(df: pd.DataFrame) -> Path:
    content = f"""# HACCP preprocessing log\n\n- Seed: {SEED}\n- Generated records: {len(df)}\n- Date range: {df['timestamp'].min()} to {df['timestamp'].max()}\n- Synthetic facility: poultry processing facility with 7 CCPs\n- Monitoring cadence: 3 to 8 readings per day depending on CCP\n- Event classes: normal, near_miss, deviation\n- Bayesian prior: Beta(alpha=2, beta=48)\n- SPC methods: Individuals chart, CUSUM (k=0.5, h=5), EWMA (lambda=0.2)\n- Output backend: matplotlib Agg, DPI={DPI}\n\n## Transformation steps\n1. Generated timestamped synthetic readings for each CCP using seasonal pressure, late-year drift, and seeded randomness.\n2. Applied CCP-specific warning and critical limits to classify normal, near-miss, and deviation states.\n3. Calculated dynamic likelihood from rolling deviation rate and posterior deviation probability.\n4. Computed RPN = Severity × Likelihood × Detectability and Bayesian-updated RPN.\n5. Aggregated weekly and monthly views for dashboard-style figures and summary metrics.\n"""
    PREPROCESS_LOG_PATH.write_text(content, encoding="utf-8")
    return PREPROCESS_LOG_PATH


def write_statistical_summary(summary: dict[str, Any]) -> Path:
    lines = [
        "# Statistical summary",
        "",
        "This module uses effect sizes, bootstrap confidence intervals, and Bayesian credible intervals rather than null-hypothesis testing because the dataset is synthetic by design.",
        "",
        "| CCP | Deviation rate | Posterior deviation probability (95% interval) | Mean Bayesian RPN | Cohen's d (H2 vs H1) | Mean change 95% CI | Trend |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for cfg in CCP_CONFIGS:
        row = summary[cfg.ccp_id]
        lines.append(
            f"| {row['ccp_name']} | {row['deviation_rate']:.3f} | {row['posterior_deviation_probability']:.3f} "
            f"({row['posterior_95pct_interval'][0]:.3f}, {row['posterior_95pct_interval'][1]:.3f}) | {row['mean_rpn']:.2f} | "
            f"{row['cohens_d_second_vs_first_half']:.2f} | ({row['mean_change_ci'][0]:.2f}, {row['mean_change_ci'][1]:.2f}) | {row['trend']} |"
        )
    STAT_SUMMARY_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return STAT_SUMMARY_PATH


def write_report(summary: dict[str, Any], figure_paths: list[Path]) -> Path:
    highest_ccp = summary[summary["overall"]["highest_final_risk_ccp"]]["ccp_name"]
    lines = [
        "# DRAFT — NOT FOR DISTRIBUTION",
        "",
        "## HACCP critical control point risk scoring automation",
        "",
        f"- Timestamp: {datetime.now(UTC).isoformat(timespec='seconds').replace('+00:00', 'Z')}",
        f"- Synthetic dataset size: {summary['overall']['records']} monitoring records across {summary['overall']['ccp_count']} CCPs",
        f"- Highest final Bayesian risk: {highest_ccp}",
        "",
        "## Assumptions",
        "- The facility is a poultry processing plant with seven defined critical control points.",
        "- Monitoring data are synthetic but structured to mimic routine operations, near-misses, and true deviations over one calendar year.",
        "- Bayesian risk updating uses a Beta prior that assumes low baseline deviation frequency and is refreshed after each observation.",
        "",
        "## Methods",
        "- RPN was calculated as Severity × Likelihood × Detectability.",
        "- Likelihood was dynamically updated from rolling deviation frequency and posterior deviation probability.",
        "- SPC analytics included Individuals control limits, EWMA smoothing, and two-sided CUSUM change detection.",
        "- Effect sizes and confidence intervals were summarized in `results/statistical-summary.md`.",
        "",
        "## Key results",
    ]
    ranked = sorted([summary[cfg.ccp_id] for cfg in CCP_CONFIGS], key=lambda item: item["final_bayesian_rpn"], reverse=True)
    for item in ranked[:4]:
        lines.append(
            f"- {item['ccp_name']}: final Bayesian RPN {item['final_bayesian_rpn']:.1f}, deviation rate {item['deviation_rate']:.3f}, trend {item['trend']}, SPC/EWMA/CUSUM signals = {item['spc_signals']}/{item['ewma_signals']}/{item['cusum_signals']}."
        )
    lines.extend(
        [
            "",
            "## Figure inventory",
            *[f"- `{path.relative_to(BASE_DIR)}`" for path in figure_paths],
            "",
            "## File inventory",
            "- `src/module4_haccp_scoring.py`",
            "- `data/haccp_monitoring.csv`",
            "- `data/preprocessing-log.md`",
            "- `results/module4_metrics.json`",
            "- `results/statistical-summary.md`",
            "- `logs/process-log.jsonl`",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return REPORT_PATH


def write_metrics(summary: dict[str, Any], figure_paths: list[Path]) -> Path:
    payload = {
        "seed": SEED,
        "generated_at": datetime.now(UTC).isoformat(timespec="seconds").replace("+00:00", "Z"),
        "artifacts": {
            "data": str(DATA_PATH.relative_to(BASE_DIR)),
            "figures": [str(path.relative_to(BASE_DIR)) for path in figure_paths],
            "report": str(REPORT_PATH.relative_to(BASE_DIR)),
            "preprocessing_log": str(PREPROCESS_LOG_PATH.relative_to(BASE_DIR)),
            "statistical_summary": str(STAT_SUMMARY_PATH.relative_to(BASE_DIR)),
        },
        "ccp_summary": {cfg.ccp_id: summary[cfg.ccp_id] for cfg in CCP_CONFIGS},
        "overall": summary["overall"],
    }
    METRICS_PATH.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return METRICS_PATH


def main() -> None:
    if PROCESS_LOG_PATH.exists():
        PROCESS_LOG_PATH.unlink()

    log_event(
        "run_started",
        "module4_haccp_scoring",
        {"seed": SEED, "requested_output": str(Path(__file__).resolve())},
        {"base_dir": str(BASE_DIR)},
        [],
    )
    log_event(
        "prompt_received",
        "co-scientist-data-analysis",
        {"task": "Create HACCP CCP risk scoring automation for poultry processing"},
        {"ccp_count": len(CCP_CONFIGS), "figure_count": 4},
        [],
    )
    log_event(
        "skill_selected",
        "co-scientist-data-analysis",
        {"candidate_skill": "co-scientist-data-analysis"},
        {"selection_reason": "Synthetic monitoring analysis, visualization, and risk scoring"},
        [],
    )
    log_event(
        "handoff_started",
        "python",
        {"phase": "execute"},
        {"outputs": [str(DATA_PATH), str(METRICS_PATH)]},
        [],
    )

    df = generate_dataset(CCP_CONFIGS)
    df, summary, diagnostics = compute_risk_metrics(df)
    df.to_csv(DATA_PATH, index=False)
    log_event("file_written", "pandas.to_csv", {}, {"rows": len(df)}, [str(DATA_PATH.relative_to(BASE_DIR))])

    figure_paths = [
        create_heatmap_figure(df),
        create_spc_figure(df, diagnostics),
        create_rpn_distribution_figure(df),
        create_bayesian_update_figure(df),
    ]
    for figure_path in figure_paths:
        log_event("file_written", "matplotlib", {}, {"dpi": DPI}, [str(figure_path.relative_to(BASE_DIR))])

    preprocess_path = write_preprocessing_log(df)
    log_event("file_written", "python", {}, {"type": "preprocessing_log"}, [str(preprocess_path.relative_to(BASE_DIR))])

    stats_path = write_statistical_summary(summary)
    log_event("file_written", "python", {}, {"type": "statistical_summary"}, [str(stats_path.relative_to(BASE_DIR))])

    report_path = write_report(summary, figure_paths)
    log_event("report_finalized", "python", {}, {"type": "report"}, [str(report_path.relative_to(BASE_DIR))])

    metrics_path = write_metrics(summary, figure_paths)
    log_event("file_written", "json", {}, {"type": "metrics"}, [str(metrics_path.relative_to(BASE_DIR))])
    log_event(
        "handoff_completed",
        "python",
        {"phase": "verify"},
        {"records": len(df), "highest_final_risk_ccp": summary['overall']['highest_final_risk_ccp']},
        [
            str(DATA_PATH.relative_to(BASE_DIR)),
            str(METRICS_PATH.relative_to(BASE_DIR)),
            *(str(path.relative_to(BASE_DIR)) for path in figure_paths),
        ],
    )
    log_event(
        "run_completed",
        "module4_haccp_scoring",
        {"status": "success"},
        {"report": str(REPORT_PATH.relative_to(BASE_DIR))},
        [str(REPORT_PATH.relative_to(BASE_DIR))],
    )

    print("HACCP monitoring dataset, risk metrics, figures, and report generated successfully.")
    print(f"Records: {len(df)}")
    print(f"Metrics: {METRICS_PATH}")


if __name__ == "__main__":
    main()
