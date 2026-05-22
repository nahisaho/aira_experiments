"""Nankai Trough-specific deformation monitoring utilities."""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter

NANKAI_PARAMETERS = {
    "plate_convergence_rate_m_per_yr": 0.06,
    "back_arc_reference_stress": 1.0,
    "sse_catalog": [
        {"name": "Bungo Channel 2010", "magnitude": 6.5},
        {"name": "Tokai 2013", "magnitude": 6.2},
        {"name": "Kii Channel 2018", "magnitude": 6.4},
    ],
}



def estimate_coupling(velocity_field: np.ndarray, plate_rate: float = NANKAI_PARAMETERS["plate_convergence_rate_m_per_yr"]) -> np.ndarray:
    """Estimate interseismic coupling coefficients from velocity magnitude."""
    vel = np.asarray(velocity_field, dtype=float)
    coupling = np.clip(np.abs(vel) / max(plate_rate, 1.0e-6), 0.0, 1.0)
    return gaussian_filter(coupling, sigma=2.0)



def detect_sse(time_series: np.ndarray, threshold: float = 0.006, min_duration: int = 3) -> list[dict[str, float | int]]:
    """Detect slow slip events as sustained transient excursions."""
    series = np.asarray(time_series, dtype=float)
    mask = np.abs(series) >= threshold
    events: list[dict[str, float | int]] = []
    start = None
    for idx, flag in enumerate(mask):
        if flag and start is None:
            start = idx
        elif not flag and start is not None:
            if idx - start >= min_duration:
                segment = series[start:idx]
                events.append({"start": start, "end": idx - 1, "duration": idx - start, "peak": float(segment[np.argmax(np.abs(segment))])})
            start = None
    if start is not None and len(series) - start >= min_duration:
        segment = series[start:]
        events.append({"start": start, "end": len(series) - 1, "duration": len(series) - start, "peak": float(segment[np.argmax(np.abs(segment))])})
    return events



def monitor_stress_accumulation(coupling_map: np.ndarray) -> dict[str, np.ndarray | dict[str, float]]:
    """Monitor back-arc stress accumulation using coupling and spatial gradients."""
    coupling = np.asarray(coupling_map, dtype=float)
    gy, gx = np.gradient(coupling)
    stress = gaussian_filter(coupling, sigma=1.5) * (1.0 + np.hypot(gx, gy))
    summary = {"mean_stress": float(np.nanmean(stress)), "max_stress": float(np.nanmax(stress))}
    return {"stress": stress, "summary": summary}



def generate_alert(coupling_map: np.ndarray, sse_events: list[dict[str, float | int]], precursor_level: str) -> dict[str, float | str | int]:
    """Generate a Nankai regional alert from coupling, SSE, and precursor information."""
    mean_coupling = float(np.nanmean(coupling_map))
    max_coupling = float(np.nanmax(coupling_map))
    level_weight = {"GREEN": 0, "YELLOW": 1, "ORANGE": 2, "RED": 3}[precursor_level]
    score = 2.0 * mean_coupling + 1.5 * max_coupling + 0.4 * len(sse_events) + level_weight
    if score < 2.0:
        level = "GREEN"
    elif score < 3.2:
        level = "YELLOW"
    elif score < 4.5:
        level = "ORANGE"
    else:
        level = "RED"
    return {
        "level": level,
        "score": float(score),
        "mean_coupling": mean_coupling,
        "max_coupling": max_coupling,
        "sse_count": len(sse_events),
    }
