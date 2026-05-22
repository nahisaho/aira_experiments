"""Generate publication-quality figures from the InSAR analysis outputs."""

from __future__ import annotations

import csv
import json
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
nmat = matplotlib
nmat.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

FIG_DIR = ROOT / "figures"
RESULTS_DIR = ROOT / "results"
LOG_PATH = ROOT / "logs" / "process-log.jsonl"



def log_event(phase: str, event_type: str, files_written: list[str], status: str = "ok", skill_or_tool: str = "generate_figures") -> None:
    record = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill_or_tool,
        "files_written": files_written,
        "status": status,
    }
    with LOG_PATH.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, ensure_ascii=False) + "\n")



def load_time_series() -> dict[str, np.ndarray]:
    path = RESULTS_DIR / "time_series_decomposition.csv"
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    return {key: np.array([float(row[key]) for row in rows]) for key in rows[0].keys()}



def save_figure(fig: plt.Figure, name: str) -> str:
    path = FIG_DIR / name
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return f"figures/{name}"



def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)
    log_event("report", "run_started", [])

    psc_mask = np.load(RESULTS_DIR / "psc_mask.npy")
    ps_velocity = np.load(RESULTS_DIR / "ps_velocity.npy")
    sbas_velocity = np.load(RESULTS_DIR / "sbas_velocity.npy")
    transient_map = np.load(RESULTS_DIR / "transient_map.npy")
    coupling_map = np.load(RESULTS_DIR / "coupling_map.npy")
    east = np.load(RESULTS_DIR / "east_final.npy")
    north = np.load(RESULTS_DIR / "north_final.npy")
    up = np.load(RESULTS_DIR / "up_final.npy")
    clusters = np.load(RESULTS_DIR / "cluster_labels.npy")
    ts = load_time_series()
    with (RESULTS_DIR / "analysis_summary.json").open("r", encoding="utf-8") as handle:
        summary = json.load(handle)

    files_written: list[str] = []
    cmap_main = "cividis"
    cmap_div = "coolwarm"

    fig1, ax1 = plt.subplots(figsize=(6, 5))
    image1 = ax1.imshow(psc_mask, cmap="viridis")
    ax1.set_title("Figure 1. PS density map")
    ax1.set_xlabel("Column")
    ax1.set_ylabel("Row")
    plt.colorbar(image1, ax=ax1, label="PSC")
    files_written.append(save_figure(fig1, "figure1_ps_density.png"))

    fig2, ax2 = plt.subplots(figsize=(6, 5))
    image2 = ax2.imshow(ps_velocity, cmap=cmap_div)
    ax2.contour(sbas_velocity, levels=8, colors="k", linewidths=0.35, alpha=0.45)
    ax2.set_title("Figure 2. Velocity field map")
    ax2.set_xlabel("Column")
    ax2.set_ylabel("Row")
    plt.colorbar(image2, ax=ax2, label="Velocity (m/yr)")
    files_written.append(save_figure(fig2, "figure2_velocity_field.png"))

    fig3, ax3 = plt.subplots(figsize=(7, 4.5))
    ax3.plot(ts["years"], ts["observed_m"], label="Observed", color="#1b9e77")
    ax3.plot(ts["years"], ts["linear_m"], label="Linear", color="#d95f02")
    ax3.plot(ts["years"], ts["seasonal_m"], label="Seasonal", color="#7570b3")
    ax3.plot(ts["years"], ts["transient_m"], label="Transient", color="#e7298a")
    ax3.set_title("Figure 3. Time-series decomposition")
    ax3.set_xlabel("Time (years)")
    ax3.set_ylabel("Displacement (m)")
    ax3.legend(frameon=False, ncol=2)
    files_written.append(save_figure(fig3, "figure3_decomposition.png"))

    fig4, (ax41, ax42) = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    ax41.plot(ts["years"], ts["strain_rate_m_per_yr"], color="#1b9e77", label="Strain rate")
    ax41.plot(ts["years"], ts["cusum_positive"], color="#d95f02", label="CUSUM+")
    ax41.plot(ts["years"], ts["cusum_negative"], color="#7570b3", label="CUSUM-")
    ax41.legend(frameon=False)
    ax41.set_ylabel("Rate / CUSUM")
    ax41.set_title(f"Figure 4. Precursor detection ({summary['precursor_alert']['level']})")
    ax42.plot(ts["years"], ts["acceleration_zscore"], color="#66a61e", label="Acceleration z-score")
    ax42.plot(ts["years"], ts["alert_score"], color="#e7298a", label="Alert score")
    ax42.set_xlabel("Time (years)")
    ax42.set_ylabel("Score")
    ax42.legend(frameon=False)
    files_written.append(save_figure(fig4, "figure4_precursor_alerts.png"))

    fig5, axes5 = plt.subplots(1, 3, figsize=(13, 4.2), sharex=True, sharey=True)
    for ax, field, title in zip(axes5, (east, north, up), ("East", "North", "Up")):
        image = ax.imshow(field, cmap=cmap_div)
        ax.set_title(title)
        ax.set_xlabel("Column")
        plt.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="Displacement (m)")
    axes5[0].set_ylabel("Row")
    fig5.suptitle("Figure 5. 3D displacement field (ENU)")
    files_written.append(save_figure(fig5, "figure5_enu_displacement.png"))

    fig6, (ax61, ax62) = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True, sharey=True)
    image61 = ax61.imshow(coupling_map, cmap=cmap_main, vmin=0.0, vmax=1.0)
    ax61.set_title("Figure 6a. Nankai coupling map")
    ax61.set_xlabel("Column")
    ax61.set_ylabel("Row")
    plt.colorbar(image61, ax=ax61, label="Coupling coefficient")
    image62 = ax62.imshow(transient_map, cmap=cmap_div)
    ax62.contour(clusters >= 0, levels=[0.5], colors="yellow", linewidths=0.8)
    ax62.set_title(f"Figure 6b. Transient clusters ({summary['nankai_alert']['level']})")
    ax62.set_xlabel("Column")
    plt.colorbar(image62, ax=ax62, label="Transient (m)")
    files_written.append(save_figure(fig6, "figure6_nankai_coupling.png"))

    log_event("report", "file_written", files_written)
    log_event("report", "run_completed", files_written)
    print("Generated figures:")
    for item in files_written:
        print(f"- {item}")


if __name__ == "__main__":
    main()
