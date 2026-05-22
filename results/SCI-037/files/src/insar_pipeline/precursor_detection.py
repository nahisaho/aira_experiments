"""Earthquake precursor detection methods for deformation time series."""

from __future__ import annotations

import numpy as np
from sklearn.cluster import DBSCAN



def detect_strain_anomaly(series: np.ndarray, threshold_std: float = 3.0, drift: float = 0.1) -> dict[str, np.ndarray | float | int | None]:
    """Detect strain-rate anomalies with two-sided CUSUM."""
    y = np.asarray(series, dtype=float)
    baseline = y[: max(5, y.size // 3)]
    mean = float(np.nanmean(baseline))
    std = float(np.nanstd(baseline) + 1.0e-6)
    pos = np.zeros_like(y)
    neg = np.zeros_like(y)
    flags = np.zeros_like(y, dtype=bool)
    limit = threshold_std * std
    for idx in range(1, y.size):
        pos[idx] = max(0.0, pos[idx - 1] + y[idx] - mean - drift * std)
        neg[idx] = max(0.0, neg[idx - 1] - (y[idx] - mean) - drift * std)
        flags[idx] = max(pos[idx], neg[idx]) > limit
    onset = int(np.argmax(flags)) if np.any(flags) else None
    severity = float(max(np.nanmax(pos), np.nanmax(neg)) / (limit + 1.0e-6))
    return {"cusum_positive": pos, "cusum_negative": neg, "flags": flags, "onset_index": onset, "severity": severity}



def detect_acceleration(series: np.ndarray, time_years: np.ndarray, z_threshold: float = 2.5) -> dict[str, np.ndarray | float]:
    """Detect deformation acceleration using robust z-scores of the second derivative."""
    y = np.asarray(series, dtype=float)
    t = np.asarray(time_years, dtype=float)
    velocity = np.gradient(y, t)
    acceleration = np.gradient(velocity, t)
    median = np.nanmedian(acceleration)
    mad = 1.4826 * np.nanmedian(np.abs(acceleration - median)) + 1.0e-6
    zscore = (acceleration - median) / mad
    flags = np.abs(zscore) >= z_threshold
    return {"acceleration": acceleration, "zscore": zscore, "flags": flags, "max_abs_z": float(np.nanmax(np.abs(zscore)))}



def spatial_clustering(anomaly_mask: np.ndarray, eps: float = 3.5, min_samples: int = 6) -> dict[str, np.ndarray | int | list[int]]:
    """Cluster spatially coherent anomalous pixels using DBSCAN."""
    mask = np.asarray(anomaly_mask, dtype=bool)
    coords = np.argwhere(mask)
    labels_image = np.full(mask.shape, -1, dtype=int)
    if len(coords) == 0:
        return {"labels": labels_image, "n_clusters": 0, "cluster_sizes": []}
    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
    labels = clustering.labels_
    labels_image[coords[:, 0], coords[:, 1]] = labels
    valid = labels[labels >= 0]
    unique, counts = np.unique(valid, return_counts=True) if valid.size else (np.array([], dtype=int), np.array([], dtype=int))
    return {"labels": labels_image, "n_clusters": int(unique.size), "cluster_sizes": counts.tolist()}



def classify_alert_level(
    cusum_result: dict[str, np.ndarray | float | int | None],
    acceleration_result: dict[str, np.ndarray | float],
    cluster_result: dict[str, np.ndarray | int | list[int]],
) -> dict[str, float | str]:
    """Classify precursor alert level from temporal and spatial anomaly metrics."""
    cluster_sizes = cluster_result.get("cluster_sizes", [])
    largest_cluster = max(cluster_sizes) if cluster_sizes else 0
    score = (
        float(cusum_result.get("severity", 0.0))
        + 0.5 * float(acceleration_result.get("max_abs_z", 0.0))
        + 0.08 * largest_cluster
        + 0.3 * int(cluster_result.get("n_clusters", 0))
    )
    if score < 2.0:
        level = "GREEN"
    elif score < 4.0:
        level = "YELLOW"
    elif score < 6.5:
        level = "ORANGE"
    else:
        level = "RED"
    return {"level": level, "score": float(score)}
