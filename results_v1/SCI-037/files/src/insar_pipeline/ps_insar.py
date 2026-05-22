"""Persistent Scatterer InSAR processing routines."""

from __future__ import annotations

import numpy as np
from scipy.fft import dctn, idctn
from scipy.spatial import Delaunay

TWO_PI = 2.0 * np.pi


def _wrap_phase(values: np.ndarray) -> np.ndarray:
    return (values + np.pi) % TWO_PI - np.pi



def select_psc(amplitude_stack: np.ndarray, threshold: float = 0.25) -> tuple[np.ndarray, np.ndarray]:
    """Select Persistent Scatterer Candidates using amplitude dispersion.

    Parameters
    ----------
    amplitude_stack : np.ndarray
        Array of shape ``(time, rows, cols)`` containing SAR amplitudes.
    threshold : float, optional
        PSC threshold for the amplitude dispersion index ``DA = sigma / mu``.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Boolean PSC mask and amplitude dispersion map.
    """
    amplitude_stack = np.asarray(amplitude_stack, dtype=float)
    if amplitude_stack.ndim != 3:
        raise ValueError("amplitude_stack must be a 3-D array [time, rows, cols]")
    mean_amp = np.nanmean(amplitude_stack, axis=0)
    std_amp = np.nanstd(amplitude_stack, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        da = np.where(mean_amp > 0, std_amp / mean_amp, np.inf)
    mask = np.isfinite(da) & (da < threshold)
    return mask, da



def estimate_ps_coherence(phase_stack: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
    """Estimate temporal coherence from wrapped phase observations."""
    phase_stack = np.asarray(phase_stack, dtype=float)
    if phase_stack.ndim != 3:
        raise ValueError("phase_stack must be a 3-D array [time, rows, cols]")
    coherence = np.abs(np.nanmean(np.exp(1j * phase_stack), axis=0))
    if mask is not None:
        if mask.shape != coherence.shape:
            raise ValueError("mask shape must match the spatial dimensions of phase_stack")
        coherence = np.where(mask, coherence, np.nan)
    return coherence



def build_ps_network(ps_mask: np.ndarray) -> np.ndarray:
    """Construct a Delaunay PS network from PSC coordinates.

    Returns an ``(n_edges, 2)`` array of node indices into ``np.argwhere(ps_mask)``.
    """
    ps_mask = np.asarray(ps_mask, dtype=bool)
    points = np.argwhere(ps_mask)
    if len(points) < 3:
        return np.empty((0, 2), dtype=int)
    tri = Delaunay(points)
    edges: set[tuple[int, int]] = set()
    for simplex in tri.simplices:
        for i in range(3):
            a = int(simplex[i])
            b = int(simplex[(i + 1) % 3])
            edges.add(tuple(sorted((a, b))))
    return np.asarray(sorted(edges), dtype=int)



def _poisson_unwrap_2d(wrapped: np.ndarray) -> np.ndarray:
    gx = _wrap_phase(np.diff(wrapped, axis=1))
    gy = _wrap_phase(np.diff(wrapped, axis=0))

    div = np.zeros_like(wrapped)
    div[:, :-1] -= gx
    div[:, 1:] += gx
    div[:-1, :] -= gy
    div[1:, :] += gy

    rhs = dctn(div, type=2, norm="ortho")
    rows, cols = wrapped.shape
    yy, xx = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    denom = 2.0 * (np.cos(np.pi * xx / cols) - 1.0) + 2.0 * (np.cos(np.pi * yy / rows) - 1.0)
    denom[0, 0] = 1.0
    sol = rhs / denom
    sol[0, 0] = 0.0
    unwrapped = idctn(sol, type=2, norm="ortho")
    phase_offset = np.nanmedian(wrapped - _wrap_phase(unwrapped))
    return unwrapped + phase_offset



def unwrap_phase(wrapped_phase: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
    """Unwrap wrapped phase using a weighted least-squares minimum-cost surrogate.

    The implementation follows a SNAPHU-like strategy by minimizing the cost of
    wrapped phase gradients through a Poisson solver.
    """
    wrapped_phase = np.asarray(wrapped_phase, dtype=float)
    if wrapped_phase.ndim not in (2, 3):
        raise ValueError("wrapped_phase must be 2-D or 3-D")
    phase_cube = wrapped_phase[None, ...] if wrapped_phase.ndim == 2 else wrapped_phase
    unwrapped = np.empty_like(phase_cube)
    for index in range(phase_cube.shape[0]):
        slice_ = phase_cube[index]
        if mask is not None:
            if mask.shape != slice_.shape:
                raise ValueError("mask shape must match phase slice")
            fill = np.nanmedian(slice_[mask]) if np.any(mask) else 0.0
            slice_ = np.where(mask, slice_, fill)
        unwrapped[index] = _poisson_unwrap_2d(slice_)
        if mask is not None:
            unwrapped[index] = np.where(mask, unwrapped[index], np.nan)
    return unwrapped[0] if wrapped_phase.ndim == 2 else unwrapped



def estimate_velocity(
    unwrapped_phase: np.ndarray,
    time_years: np.ndarray,
    wavelength: float = 0.056,
    dem_sensitivity: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Estimate LOS velocity and DEM error using least squares.

    Parameters
    ----------
    unwrapped_phase : np.ndarray
        Unwrapped phase cube of shape ``(time, rows, cols)``.
    time_years : np.ndarray
        Acquisition times in decimal years relative to the first acquisition.
    wavelength : float, optional
        Radar wavelength in meters.
    dem_sensitivity : np.ndarray, optional
        1-D temporal sensitivity vector linked to perpendicular baselines.
    """
    unwrapped_phase = np.asarray(unwrapped_phase, dtype=float)
    time_years = np.asarray(time_years, dtype=float)
    if unwrapped_phase.ndim != 3:
        raise ValueError("unwrapped_phase must be 3-D")
    if unwrapped_phase.shape[0] != time_years.size:
        raise ValueError("time_years length must match the time dimension")

    displacement = unwrapped_phase * wavelength / (4.0 * np.pi)
    pixels = displacement.reshape(displacement.shape[0], -1)
    if dem_sensitivity is None:
        dem_sensitivity = np.zeros_like(time_years)
    dem_sensitivity = np.asarray(dem_sensitivity, dtype=float)
    if dem_sensitivity.shape != time_years.shape:
        raise ValueError("dem_sensitivity must match time_years shape")
    design = np.column_stack([time_years, dem_sensitivity, np.ones_like(time_years)])
    valid = np.all(np.isfinite(pixels), axis=0)
    coeffs = np.full((3, pixels.shape[1]), np.nan)
    if np.any(valid):
        pinv = np.linalg.pinv(design)
        coeffs[:, valid] = pinv @ pixels[:, valid]
    model = (design @ np.nan_to_num(coeffs)).reshape(unwrapped_phase.shape)
    residual = displacement - model
    velocity = coeffs[0].reshape(unwrapped_phase.shape[1:])
    dem_error = coeffs[1].reshape(unwrapped_phase.shape[1:])
    intercept = coeffs[2].reshape(unwrapped_phase.shape[1:])
    rms = np.sqrt(np.nanmean(residual**2, axis=0))
    return {
        "velocity": velocity,
        "dem_error": dem_error,
        "intercept": intercept,
        "residual_rms": rms,
        "modeled_displacement": model,
    }
