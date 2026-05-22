#!/usr/bin/env python3
"""Generate a biomarker monitoring strategy for brain organoid bioreactors."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches

SEED = 9303
RNG = np.random.default_rng(SEED)
BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = BASE_DIR / "figures"
DATA_CSV = DATA_DIR / "biomarker_timecourse.csv"
PROTOCOL_MD = RESULTS_DIR / "monitoring_protocol.md"


def logistic(x: np.ndarray, midpoint: float, slope: float, low: float, high: float) -> np.ndarray:
    """Smooth monotonic transition."""
    return low + (high - low) / (1.0 + np.exp(-(x - midpoint) / slope))


def gaussian(x: np.ndarray, center: float, width: float, amplitude: float, baseline: float = 0.0) -> np.ndarray:
    """Bell-shaped biomarker profile."""
    return baseline + amplitude * np.exp(-0.5 * ((x - center) / width) ** 2)


def add_noise(values: np.ndarray, sigma: float, floor: float = 0.0) -> np.ndarray:
    """Add Gaussian measurement noise while keeping non-negative signals."""
    return np.clip(values + RNG.normal(0.0, sigma, size=len(values)), floor, None)


def shewhart_limits(values: np.ndarray, baseline_mask: np.ndarray) -> tuple[float, float, float]:
    """Compute mean and three-sigma control limits."""
    mean = float(values[baseline_mask].mean())
    sigma = float(values[baseline_mask].std(ddof=1))
    return mean, mean - 3.0 * sigma, mean + 3.0 * sigma


def cusum(values: np.ndarray, target: float, sigma: float, k: float = 0.5, h: float = 5.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """One-sided CUSUM for upward and downward shifts."""
    if sigma < 1e-9:
        sigma = 1e-9
    pos = np.zeros_like(values)
    neg = np.zeros_like(values)
    flags = np.zeros_like(values, dtype=bool)
    for i in range(1, len(values)):
        pos[i] = max(0.0, pos[i - 1] + (values[i] - target) / sigma - k)
        neg[i] = max(0.0, neg[i - 1] - (values[i] - target) / sigma - k)
        flags[i] = (pos[i] > h) or (neg[i] > h)
    return pos, neg, flags


def build_dataset() -> Dict[str, np.ndarray]:
    """Generate synthetic online, at-line, and offline monitoring data."""
    day = np.arange(0.0, 120.5, 0.5)
    weekly_mask = np.isclose(day % 7.0, 0.0)
    daily_mask = np.isclose(day % 1.0, 0.0)
    baseline_mask = day <= 20.0

    anomaly_window = ((day >= 62.0) & (day <= 68.0)).astype(float)
    drift_window = np.clip((day - 92.0) / 16.0, 0.0, 1.0)

    data = {
        "day": day,
        "OCT4": np.where(weekly_mask, add_noise(logistic(day, 8.0, -3.2, 0.1, 1.0), 0.03), np.nan),
        "NANOG": np.where(weekly_mask, add_noise(logistic(day, 10.0, -3.8, 0.08, 0.95), 0.03), np.nan),
        "PAX6": np.where(weekly_mask, add_noise(gaussian(day, 18.0, 8.0, 0.85, 0.08), 0.03), np.nan),
        "SOX2": np.where(weekly_mask, add_noise(gaussian(day, 20.0, 9.0, 0.82, 0.10), 0.03), np.nan),
        "NESTIN": np.where(weekly_mask, add_noise(gaussian(day, 23.0, 10.0, 0.88, 0.12), 0.03), np.nan),
        "TBR1": np.where(weekly_mask, add_noise(logistic(day, 38.0, 6.5, 0.05, 0.92), 0.03), np.nan),
        "CTIP2": np.where(weekly_mask, add_noise(logistic(day, 45.0, 7.0, 0.06, 0.90), 0.03), np.nan),
        "SATB2": np.where(weekly_mask, add_noise(logistic(day, 55.0, 8.0, 0.05, 0.88), 0.03), np.nan),
        "MAP2": np.where(weekly_mask, add_noise(logistic(day, 48.0, 7.0, 0.08, 1.15), 0.04), np.nan),
        "SYN1": np.where(weekly_mask, add_noise(logistic(day, 56.0, 8.0, 0.05, 1.05), 0.04), np.nan),
        "GFAP": np.where(weekly_mask, add_noise(logistic(day, 63.0, 9.0, 0.05, 0.95), 0.04), np.nan),
    }

    glucose_consumption = add_noise(0.14 + 0.18 * logistic(day, 18.0, 5.5, 0.0, 1.0) + 0.12 * logistic(day, 48.0, 8.0, 0.0, 1.0), 0.01)
    lactate_production = add_noise(0.10 + 0.14 * logistic(day, 24.0, 6.5, 0.0, 1.0) + 0.10 * logistic(day, 52.0, 8.5, 0.0, 1.0), 0.01)
    ldh = add_noise(180.0 + 40.0 * logistic(day, 80.0, 8.0, 0.0, 1.0) + 75.0 * anomaly_window, 8.0)
    cytokine_panel = add_noise(1.2 + 0.9 * logistic(day, 45.0, 9.0, 0.0, 1.0) + 0.25 * anomaly_window, 0.08)

    data["glucose_consumption_rate"] = glucose_consumption
    data["lactate_production_rate"] = lactate_production
    data["LDH"] = np.where(daily_mask, ldh, np.nan)
    data["cytokine_panel"] = np.where(daily_mask, cytokine_panel, np.nan)

    ph = add_noise(7.34 - 0.04 * logistic(day, 44.0, 10.0, 0.0, 1.0) - 0.12 * anomaly_window - 0.04 * drift_window, 0.015)
    do = add_noise(66.0 - 6.0 * logistic(day, 48.0, 8.5, 0.0, 1.0) - 12.0 * anomaly_window - 3.0 * drift_window, 1.1)
    glucose = add_noise(17.5 - 4.6 * logistic(day, 30.0, 7.5, 0.0, 1.0) - 2.8 * anomaly_window - 1.4 * drift_window, 0.45)
    lactate = add_noise(1.1 + 5.8 * logistic(day, 36.0, 8.5, 0.0, 1.0) + 1.6 * anomaly_window + 0.9 * drift_window, 0.35)

    data["pH"] = ph
    data["DO_percent"] = do
    data["glucose_mM"] = glucose
    data["lactate_mM"] = lactate

    standardized = np.column_stack(
        [
            (ph - ph[baseline_mask].mean()) / ph[baseline_mask].std(ddof=1),
            (do - do[baseline_mask].mean()) / do[baseline_mask].std(ddof=1),
            (glucose - glucose[baseline_mask].mean()) / glucose[baseline_mask].std(ddof=1),
            (lactate - lactate[baseline_mask].mean()) / lactate[baseline_mask].std(ddof=1),
        ]
    )
    anomaly_score = np.sqrt(np.sum(standardized**2, axis=1))
    anomaly_flag = anomaly_score > 4.5
    data["anomaly_score"] = anomaly_score
    data["anomaly_flag"] = anomaly_flag.astype(int)
    return data


def save_csv(data: Dict[str, np.ndarray]) -> None:
    """Save monitoring dataset to CSV."""
    columns = list(data.keys())
    matrix = np.column_stack([data[col] for col in columns])
    header = ",".join(columns)
    np.savetxt(DATA_CSV, matrix, delimiter=",", header=header, comments="", fmt="%.6f")


def save_protocol() -> None:
    """Write markdown monitoring protocol."""
    protocol = """# Brain Organoid Bioreactor Monitoring Protocol\n\n## Scope\nThis protocol integrates online, at-line, and offline biomarker monitoring for maturation control in brain organoid bioreactors.\n\n## Monitoring cadence\n- **Online (continuous, 30 min acquisition):** pH, dissolved oxygen (DO), glucose, lactate\n- **At-line (daily):** LDH release, multiplex cytokine panel\n- **Offline (weekly):** qPCR and immunostaining for OCT4, NANOG, PAX6, SOX2, NESTIN, TBR1, CTIP2, SATB2, MAP2, SYN1, GFAP; weekly electrophysiology benchmarking\n\n## Release criteria by phase\n1. **Neural induction (day 0-6):** OCT4 and NANOG falling, pH 7.25-7.40, DO > 55%\n2. **Patterning (day 6-25):** PAX6/SOX2/NESTIN peak trajectory, glucose consumption increasing without LDH surge\n3. **Cortical differentiation (day 25-50):** TBR1 and CTIP2 increasing, lactate rise matched by stable DO\n4. **Maturation (day 50+):** MAP2/SYN1/GFAP increase, low LDH, sustained electrophysiology score\n\n## Process-control actions\n- Shewhart rule violation: verify sensor health, then inspect perfusion and aeration hardware within 2 h\n- CUSUM trigger without Shewhart breach: check slow drift in feed composition, calibration, and gas blending\n- Multivariate anomaly score > 4.5: hold automated feed escalation and perform targeted microscopy/qPCR confirmation\n\n## Escalation logic\n- DO low + lactate high -> increase gas transfer and inspect aggregate density\n- Glucose low + LDH high -> reduce residence time, refresh medium, and assess necrotic cores\n- Cytokine elevation without sensor deviation -> inspect contamination or inflammatory glial overgrowth\n"""
    PROTOCOL_MD.write_text(protocol, encoding="utf-8")


def plot_biomarkers(data: Dict[str, np.ndarray]) -> None:
    """Plot multi-panel temporal biomarker categories."""
    day = data["day"]
    groups = {
        "Pluripotency": ["OCT4", "NANOG"],
        "Neural progenitor": ["PAX6", "SOX2", "NESTIN"],
        "Cortical": ["TBR1", "CTIP2", "SATB2"],
        "Maturation": ["MAP2", "SYN1", "GFAP"],
        "Metabolic": ["glucose_consumption_rate", "lactate_production_rate", "LDH"],
    }
    fig, axes = plt.subplots(3, 2, figsize=(13, 11), constrained_layout=True)
    axes = axes.flatten()
    cmaps = [plt.get_cmap("viridis"), plt.get_cmap("cividis"), plt.get_cmap("viridis"), plt.get_cmap("cividis"), plt.get_cmap("viridis")]
    for ax, (title, markers), cmap in zip(axes, groups.items(), cmaps):
        for idx, marker in enumerate(markers):
            ax.plot(day, data[marker], marker="o", markersize=2.5, linewidth=1.8, label=marker, color=cmap(0.2 + 0.25 * idx))
        ax.set_title(title)
        ax.set_xlabel("Culture day")
        ax.set_ylabel("Signal (a.u.)")
        ax.grid(alpha=0.25)
        ax.legend(frameon=False, fontsize=8)
    axes[-1].axis("off")
    axes[-1].text(0.03, 0.70, "Weekly markers are sampled intermittently;\nat-line LDH is daily; metabolic rates are modeled continuously.", fontsize=11)
    axes[-1].text(0.03, 0.38, "Synthetic anomaly introduced at days 62-68\nto exercise control-chart and anomaly-detection logic.", fontsize=11, weight="bold")
    fig.suptitle("Temporal biomarker trajectories for brain organoid maturation", fontsize=15, weight="bold")
    fig.savefig(FIGURES_DIR / "biomarker_timecourse.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_control_charts(data: Dict[str, np.ndarray]) -> None:
    """Plot Shewhart control charts and annotate CUSUM alarms."""
    day = data["day"]
    baseline_mask = day <= 20.0
    sensors = [("pH", "pH"), ("DO_percent", "DO (%)"), ("glucose_mM", "Glucose (mM)"), ("lactate_mM", "Lactate (mM)")]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for ax, (key, ylabel), color in zip(axes.flatten(), sensors, ["#003f5c", "#58508d", "#bc5090", "#ff6361"]):
        values = data[key]
        mean, lcl, ucl = shewhart_limits(values, baseline_mask)
        pos, neg, flags = cusum(values, mean, values[baseline_mask].std(ddof=1))
        shewhart = (values < lcl) | (values > ucl)
        ax.plot(day, values, color=color, linewidth=1.8)
        ax.axhline(mean, color="black", linestyle="--", linewidth=1.0, label="Center line")
        ax.axhline(lcl, color="#2f4b7c", linestyle=":", linewidth=1.0, label="Control limits")
        ax.axhline(ucl, color="#2f4b7c", linestyle=":", linewidth=1.0)
        ax.scatter(day[shewhart], values[shewhart], color="#ffa600", s=20, label="Shewhart alert")
        ax.scatter(day[flags], values[flags], facecolors="none", edgecolors="#003f5c", s=45, label="CUSUM alert")
        ax.set_title(ylabel)
        ax.set_xlabel("Culture day")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False)
    fig.suptitle("Online sensor control charts for process monitoring", fontsize=15, weight="bold")
    fig.savefig(FIGURES_DIR / "control_charts.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def draw_box(ax, xy, text, width=0.23, height=0.12, fc="#d6e5fa"):
    """Draw a flowchart box."""
    rect = patches.FancyBboxPatch(xy, width, height, boxstyle="round,pad=0.02", facecolor=fc, edgecolor="#2f4b7c", linewidth=1.5)
    ax.add_patch(rect)
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center", fontsize=10)
    return rect


def arrow(ax, start, end):
    """Draw an arrow between flowchart nodes."""
    ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", color="#2f4b7c", linewidth=1.6))


def plot_decision_tree() -> None:
    """Create monitoring decision flowchart with matplotlib."""
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    draw_box(ax, (0.38, 0.86), "Sensor update")
    draw_box(ax, (0.38, 0.68), "Shewhart or CUSUM breach?", fc="#fbe7c6")
    draw_box(ax, (0.08, 0.47), "No -> continue scheduled sampling", fc="#d7f0d1")
    draw_box(ax, (0.38, 0.47), "Verify calibration\nand hardware state", fc="#fbe7c6")
    draw_box(ax, (0.68, 0.47), "Multivariate anomaly\nscore > 4.5?", fc="#fbe7c6")
    draw_box(ax, (0.08, 0.23), "Adjust aeration / perfusion\nand recheck in 2 h", fc="#cfe8ff")
    draw_box(ax, (0.38, 0.23), "Collect LDH, cytokines,\nmicroscopy, qPCR", fc="#cfe8ff")
    draw_box(ax, (0.68, 0.23), "Escalate to process hold\nand root-cause review", fc="#f8d7da")
    arrow(ax, (0.50, 0.86), (0.50, 0.80))
    arrow(ax, (0.50, 0.68), (0.20, 0.59))
    arrow(ax, (0.50, 0.68), (0.50, 0.59))
    arrow(ax, (0.50, 0.68), (0.80, 0.59))
    arrow(ax, (0.20, 0.47), (0.20, 0.35))
    arrow(ax, (0.50, 0.47), (0.50, 0.35))
    arrow(ax, (0.80, 0.47), (0.80, 0.35))
    ax.text(0.26, 0.62, "No", fontsize=10, weight="bold")
    ax.text(0.53, 0.62, "Yes", fontsize=10, weight="bold")
    ax.text(0.74, 0.62, "Check anomaly", fontsize=10, weight="bold")
    fig.suptitle("Decision tree for bioreactor process control", fontsize=15, weight="bold")
    fig.savefig(FIGURES_DIR / "monitoring_decision_tree.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    """Run monitoring strategy generation and save outputs."""
    DATA_DIR.mkdir(exist_ok=True)
    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)
    data = build_dataset()
    save_csv(data)
    save_protocol()
    plot_biomarkers(data)
    plot_control_charts(data)
    plot_decision_tree()
    print("Biomarker monitoring assets generated.")


if __name__ == "__main__":
    main()
