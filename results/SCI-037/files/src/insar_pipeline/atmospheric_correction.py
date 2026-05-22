"""Atmospheric phase delay correction utilities."""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import pdist



def _as_cube(data: np.ndarray) -> np.ndarray:
    arr = np.asarray(data, dtype=float)
    if arr.ndim == 2:
        return arr[None, ...]
    if arr.ndim != 3:
        raise ValueError("Input must be 2-D or 3-D")
    return arr



def era5_correction(los_stack: np.ndarray, elevation: np.ndarray, incidence_angle: float | np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Apply a simplified ERA5-style stratified tropospheric correction.

    The delay is modeled as an elevation-dependent zenith delay mapped to the LOS.
    """
    cube = _as_cube(los_stack)
    elevation = np.asarray(elevation, dtype=float)
    incidence = np.asarray(incidence_angle, dtype=float)
    if elevation.shape != cube.shape[1:]:
        raise ValueError("elevation shape must match spatial dimensions of los_stack")
    cos_inc = np.cos(incidence)
    cos_inc = np.where(np.abs(cos_inc) < 1.0e-6, 1.0e-6, cos_inc)
    elev_km = (elevation - np.nanmean(elevation)) / 1000.0
    humidity_scale = 0.0025 + 0.0008 * np.sin(np.linspace(0, 2.0 * np.pi, cube.shape[0], endpoint=False))
    delay = humidity_scale[:, None, None] * elev_km[None, ...] / cos_inc
    corrected = cube - delay
    return (corrected[0] if los_stack.ndim == 2 else corrected, delay[0] if los_stack.ndim == 2 else delay)



def gacos_correction(los_stack: np.ndarray, elevation: np.ndarray, exponent: float = 1.2) -> tuple[np.ndarray, np.ndarray]:
    """Apply a GACOS-style power-law tropospheric correction."""
    cube = _as_cube(los_stack)
    elevation = np.asarray(elevation, dtype=float)
    if elevation.shape != cube.shape[1:]:
        raise ValueError("elevation shape must match spatial dimensions of los_stack")
    basis = np.sign(elevation - np.nanmedian(elevation)) * np.abs((elevation - np.nanmedian(elevation)) / 1000.0) ** exponent
    basis -= np.nanmean(basis)
    denom = np.nansum(basis**2)
    if denom <= 0:
        correction = np.zeros_like(cube)
    else:
        alpha = np.nansum(cube * basis[None, ...], axis=(1, 2)) / denom
        correction = alpha[:, None, None] * basis[None, ...]
    corrected = cube - correction
    return (corrected[0] if los_stack.ndim == 2 else corrected, correction[0] if los_stack.ndim == 2 else correction)



def spatial_filter(los_stack: np.ndarray, sigma_low: float = 12.0, sigma_high: float = 2.0) -> np.ndarray:
    """Band-pass filter the deformation cube to suppress long-wavelength APS."""
    cube = _as_cube(los_stack)
    filtered = np.empty_like(cube)
    for idx, slice_ in enumerate(cube):
        low = gaussian_filter(slice_, sigma=sigma_low, mode="nearest")
        high = gaussian_filter(slice_, sigma=sigma_high, mode="nearest")
        filtered[idx] = slice_ - low + high
    return filtered[0] if los_stack.ndim == 2 else filtered



def estimate_aps(los_stack: np.ndarray, sample_size: int = 400) -> tuple[np.ndarray, dict[str, float | list[float]]]:
    """Estimate atmospheric phase screens using a variogram-derived correlation range."""
    cube = _as_cube(los_stack)
    mean_field = np.nanmean(cube, axis=0)
    coords = np.argwhere(np.isfinite(mean_field))
    if coords.size == 0:
        aps = np.zeros_like(cube)
        return aps[0] if los_stack.ndim == 2 else aps, {"range_pixels": 1.0, "lags": [], "semivariance": []}
    if len(coords) > sample_size:
        step = max(len(coords) // sample_size, 1)
        coords = coords[::step]
    values = mean_field[coords[:, 0], coords[:, 1]]
    if len(coords) > 1:
        distances = pdist(coords.astype(float))
        semivariance = 0.5 * pdist(values[:, None], metric="sqeuclidean")
        bins = np.linspace(0, max(distances.max(), 1.0), 8)
        lag_means, semi_means = [], []
        for start, stop in zip(bins[:-1], bins[1:]):
            mask = (distances >= start) & (distances < stop)
            if np.any(mask):
                lag_means.append(float(distances[mask].mean()))
                semi_means.append(float(semivariance[mask].mean()))
        sill = float(np.nanmax(semi_means)) if semi_means else 1.0
        threshold = 0.95 * sill
        range_pixels = next((lag for lag, semi in zip(lag_means, semi_means) if semi >= threshold), 12.0)
    else:
        lag_means, semi_means, range_pixels = [], [], 12.0
    sigma = max(range_pixels / 4.0, 1.0)
    aps = np.stack([gaussian_filter(slice_ - np.nanmedian(slice_), sigma=sigma, mode="nearest") for slice_ in cube])
    meta = {"range_pixels": float(range_pixels), "lags": lag_means, "semivariance": semi_means}
    return (aps[0] if los_stack.ndim == 2 else aps, meta)
