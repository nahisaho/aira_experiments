"""Trend decomposition tools for crustal deformation time series."""

from __future__ import annotations

import numpy as np
from scipy import stats



def fit_linear_trend(time_years: np.ndarray, series: np.ndarray) -> dict[str, np.ndarray | float]:
    """Fit a linear trend and return slope, intercept, residuals, and 95% CI."""
    t = np.asarray(time_years, dtype=float)
    y = np.asarray(series, dtype=float)
    if t.ndim != 1 or y.ndim != 1 or t.size != y.size:
        raise ValueError("time_years and series must be matching 1-D arrays")
    design = np.column_stack([t, np.ones_like(t)])
    coeffs, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ coeffs
    residual = y - fitted
    dof = max(t.size - 2, 1)
    sigma2 = np.sum(residual**2) / dof
    sxx = np.sum((t - t.mean()) ** 2)
    slope_se = np.sqrt(sigma2 / sxx) if sxx > 0 else np.nan
    ci = coeffs[0] + np.array([-1.96, 1.96]) * slope_se
    return {
        "slope": float(coeffs[0]),
        "intercept": float(coeffs[1]),
        "fitted": fitted,
        "residual": residual,
        "slope_ci": ci,
    }



def fit_seasonal(time_years: np.ndarray, series: np.ndarray) -> dict[str, np.ndarray | float]:
    """Fit annual and semi-annual harmonics to the input series."""
    t = np.asarray(time_years, dtype=float)
    y = np.asarray(series, dtype=float)
    design = np.column_stack(
        [
            np.sin(2.0 * np.pi * t),
            np.cos(2.0 * np.pi * t),
            np.sin(4.0 * np.pi * t),
            np.cos(4.0 * np.pi * t),
            np.ones_like(t),
        ]
    )
    coeffs, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ coeffs
    annual_amp = float(np.hypot(coeffs[0], coeffs[1]))
    semi_amp = float(np.hypot(coeffs[2], coeffs[3]))
    return {
        "seasonal": fitted,
        "coefficients": coeffs,
        "annual_amplitude": annual_amp,
        "semi_annual_amplitude": semi_amp,
    }



def extract_transient(time_years: np.ndarray, series: np.ndarray, delta: float = 1.5, max_iter: int = 30) -> dict[str, np.ndarray | float]:
    """Extract transient deformation using Huber-robust harmonic regression."""
    t = np.asarray(time_years, dtype=float)
    y = np.asarray(series, dtype=float)
    design = np.column_stack(
        [
            t,
            np.ones_like(t),
            np.sin(2.0 * np.pi * t),
            np.cos(2.0 * np.pi * t),
            np.sin(4.0 * np.pi * t),
            np.cos(4.0 * np.pi * t),
        ]
    )
    weights = np.ones_like(y)
    beta = np.zeros(design.shape[1])
    for _ in range(max_iter):
        wdesign = design * weights[:, None]
        wy = y * weights
        beta_new, _, _, _ = np.linalg.lstsq(wdesign, wy, rcond=None)
        residual = y - design @ beta_new
        scale = 1.4826 * np.median(np.abs(residual - np.median(residual))) + 1.0e-6
        scaled = np.abs(residual) / scale
        new_weights = np.where(scaled <= delta, 1.0, delta / scaled)
        if np.allclose(beta, beta_new, atol=1.0e-8):
            beta = beta_new
            break
        beta, weights = beta_new, new_weights
    linear = beta[0] * t + beta[1]
    seasonal = design[:, 2:] @ beta[2:]
    transient = y - (linear + seasonal)
    return {"linear": linear, "seasonal": seasonal, "transient": transient, "coefficients": beta}



def kalman_filter_decompose(time_years: np.ndarray, series: np.ndarray, process_var: float = 1.0e-4, obs_var: float | None = None) -> dict[str, np.ndarray]:
    """Track a local linear trend using a 2-state Kalman filter."""
    t = np.asarray(time_years, dtype=float)
    y = np.asarray(series, dtype=float)
    dt = np.diff(t, prepend=t[0])
    if obs_var is None:
        obs_var = float(np.nanvar(np.diff(y))) if y.size > 2 else 1.0e-4
        obs_var = max(obs_var, 1.0e-6)
    state = np.array([y[0], 0.0])
    cov = np.eye(2)
    levels = np.zeros_like(y)
    slopes = np.zeros_like(y)
    for idx, step in enumerate(dt):
        transition = np.array([[1.0, step], [0.0, 1.0]])
        process = process_var * np.array([[step**4 / 4.0, step**3 / 2.0], [step**3 / 2.0, step**2]])
        state = transition @ state
        cov = transition @ cov @ transition.T + process
        innovation = y[idx] - state[0]
        s = cov[0, 0] + obs_var
        gain = cov[:, 0] / s
        state = state + gain * innovation
        cov = cov - np.outer(gain, cov[0, :])
        levels[idx], slopes[idx] = state
    return {"level": levels, "slope": slopes}
