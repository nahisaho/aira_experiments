"""Generate synthetic InSAR data for Nankai Trough deformation monitoring."""

from __future__ import annotations

import csv
import json
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
LOG_PATH = ROOT / "logs" / "process-log.jsonl"



def log_event(phase: str, event_type: str, files_written: list[str], status: str = "ok", skill_or_tool: str = "generate_synthetic_data") -> None:
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



def wrap_phase(values: np.ndarray) -> np.ndarray:
    return (values + np.pi) % (2.0 * np.pi) - np.pi



def main() -> None:
    rng = np.random.default_rng(42)
    DATA_DIR.mkdir(exist_ok=True)
    (ROOT / "logs").mkdir(exist_ok=True)
    log_event("plan", "run_started", [])
    log_event("plan", "prompt_received", [])

    rows = cols = 100
    n_times = 50
    wavelength = 0.056
    start_date = datetime(2017, 1, 1)
    acquisition_days = np.arange(n_times, dtype=float) * 24.0
    acquisition_years = acquisition_days / 365.25
    dates = np.array([(start_date + timedelta(days=float(day))).strftime("%Y-%m-%d") for day in acquisition_days])
    perp_baselines = 110.0 * np.sin(np.linspace(0.0, 3.0 * np.pi, n_times)) + rng.normal(0.0, 22.0, n_times)

    yy, xx = np.indices((rows, cols))
    x_norm = (xx - cols / 2) / (cols / 2)
    y_norm = (yy - rows / 2) / (rows / 2)
    elevation = 120.0 + 2.5 * yy + 1.3 * xx + 45.0 * np.sin(xx / 14.0) + 25.0 * np.cos(yy / 17.0)

    sse_spatial = np.exp(-(((xx - 58.0) / 13.0) ** 2 + ((yy - 62.0) / 10.0) ** 2) / 2.0)
    sse_mask = sse_spatial > 0.45
    coastal_stable = ((yy < 35) & (xx > 12) & (xx < 88)) | ((xx < 28) & (yy > 42)) | sse_mask

    amplitude_mean = np.where(coastal_stable, 1200.0, 900.0)
    amplitude_sigma = np.where(coastal_stable, 90.0, 320.0)
    amplitude_stack = amplitude_mean[None, ...] + rng.normal(0.0, 1.0, size=(n_times, rows, cols)) * amplitude_sigma[None, ...]
    amplitude_stack = np.clip(amplitude_stack, 50.0, None)

    interseismic_up_rate = -0.010 - 0.004 * np.exp(-(((xx - 65.0) ** 2 + (yy - 32.0) ** 2) / 900.0))
    interseismic_east_rate = 0.004 * x_norm * (1.0 + 0.3 * np.cos(yy / 21.0))
    interseismic_north_rate = 0.0015 * np.sin(yy / 18.0)

    slow_slip = 0.022 * (1.0 / (1.0 + np.exp(-(acquisition_days - 620.0) / 22.0)) - 1.0 / (1.0 + np.exp(-(acquisition_days - 820.0) / 26.0)))
    seasonal_annual = 0.0045 * np.sin(2.0 * np.pi * acquisition_years - 0.25)
    seasonal_semi = 0.0018 * np.cos(4.0 * np.pi * acquisition_years + 0.4)
    seasonal_spatial = 0.65 + 0.35 * np.cos(xx / 20.0) * np.sin(yy / 18.0)

    east = acquisition_years[:, None, None] * interseismic_east_rate[None, ...] + 0.25 * slow_slip[:, None, None] * sse_spatial[None, ...]
    north = acquisition_years[:, None, None] * interseismic_north_rate[None, ...] + 0.10 * slow_slip[:, None, None] * sse_spatial[None, ...]
    up = (
        acquisition_years[:, None, None] * interseismic_up_rate[None, ...]
        + slow_slip[:, None, None] * sse_spatial[None, ...]
        + (seasonal_annual + seasonal_semi)[:, None, None] * seasonal_spatial[None, ...]
    )

    asc_incidence = np.deg2rad(34.0)
    desc_incidence = np.deg2rad(37.0)
    asc_azimuth = np.deg2rad(100.0)
    desc_azimuth = np.deg2rad(280.0)
    asc_vector = np.array([-np.sin(asc_incidence) * np.sin(asc_azimuth), np.sin(asc_incidence) * np.cos(asc_azimuth), np.cos(asc_incidence)])
    desc_vector = np.array([-np.sin(desc_incidence) * np.sin(desc_azimuth), np.sin(desc_incidence) * np.cos(desc_azimuth), np.cos(desc_incidence)])

    los_asc_true = east * asc_vector[0] + north * asc_vector[1] + up * asc_vector[2]
    los_desc_true = east * desc_vector[0] + north * desc_vector[1] + up * desc_vector[2]

    stratified = ((elevation - np.mean(elevation)) / 1000.0)[None, ...]
    time_wave = np.sin(np.linspace(0.0, 2.0 * np.pi, n_times, endpoint=False))[:, None, None]
    turbulent = gaussian_filter(rng.normal(0.0, 0.003, size=(n_times, rows, cols)), sigma=(0.0, 7.0, 7.0))
    aps = 0.006 * stratified * time_wave + turbulent
    noise_asc = rng.normal(0.0, 0.0018, size=(n_times, rows, cols))
    noise_desc = rng.normal(0.0, 0.0018, size=(n_times, rows, cols))

    los_asc_observed = los_asc_true + aps + noise_asc
    los_desc_observed = los_desc_true + 0.92 * aps + noise_desc
    wrapped_phase_asc = wrap_phase(4.0 * np.pi * los_asc_observed / wavelength)
    wrapped_phase_desc = wrap_phase(4.0 * np.pi * los_desc_observed / wavelength)

    gps_sites = np.array([[20, 18], [38, 74], [56, 58], [72, 30], [84, 82]], dtype=float)
    gps_points = []
    for row, col in gps_sites:
        r, c = int(row), int(col)
        gps_points.append(
            [
                r,
                c,
                east[-1, r, c] + rng.normal(0.0, 0.0010),
                north[-1, r, c] + rng.normal(0.0, 0.0010),
                up[-1, r, c] + rng.normal(0.0, 0.0012),
            ]
        )
    gps_points = np.asarray(gps_points, dtype=float)

    files_written = []
    arrays = {
        "amplitude_stack.npy": amplitude_stack,
        "elevation.npy": elevation,
        "times_days.npy": acquisition_days,
        "times_years.npy": acquisition_years,
        "dates.npy": dates,
        "perp_baselines.npy": perp_baselines,
        "true_east.npy": east,
        "true_north.npy": north,
        "true_up.npy": up,
        "los_asc_true.npy": los_asc_true,
        "los_desc_true.npy": los_desc_true,
        "los_asc_observed.npy": los_asc_observed,
        "los_desc_observed.npy": los_desc_observed,
        "wrapped_phase_asc.npy": wrapped_phase_asc,
        "wrapped_phase_desc.npy": wrapped_phase_desc,
        "sse_mask.npy": sse_mask.astype(np.uint8),
        "gps_points.npy": gps_points,
    }
    for name, array in arrays.items():
        np.save(DATA_DIR / name, array)
        files_written.append(f"data/{name}")

    metadata_path = DATA_DIR / "metadata.csv"
    with metadata_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_index", "date", "days_since_start", "decimal_years", "perpendicular_baseline_m", "asc_incidence_deg", "desc_incidence_deg"])
        for idx in range(n_times):
            writer.writerow([idx, dates[idx], acquisition_days[idx], acquisition_years[idx], perp_baselines[idx], np.rad2deg(asc_incidence), np.rad2deg(desc_incidence)])
    files_written.append("data/metadata.csv")

    gps_csv = DATA_DIR / "gps_points.csv"
    with gps_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["row", "col", "east_m", "north_m", "up_m"])
        writer.writerows(gps_points.tolist())
    files_written.append("data/gps_points.csv")

    preprocessing_log = DATA_DIR / "preprocessing-log.md"
    preprocessing_log.write_text(
        "# Preprocessing Log\n\n"
        "- Random seeds fixed with numpy default_rng(42).\n"
        "- 100x100 grid and 50 acquisitions with 24-day repeat cycle.\n"
        "- Synthetic signals: interseismic trend, slow slip event, annual/semi-annual deformation, atmospheric delays, and sensor noise.\n"
        "- LOS geometries represent ascending and descending Sentinel-1-like tracks.\n"
        "- GPS control points were sampled from the true ENU field with small Gaussian perturbations.\n",
        encoding="utf-8",
    )
    files_written.append("data/preprocessing-log.md")

    log_event("execute", "file_written", files_written)
    log_event("execute", "run_completed", files_written)
    print(f"Synthetic dataset written to {DATA_DIR}")
    print(f"Grid: {rows}x{cols}, acquisitions: {n_times}")
    print(f"SSE pixels: {int(sse_mask.sum())}, GPS stations: {len(gps_points)}")


if __name__ == "__main__":
    main()
