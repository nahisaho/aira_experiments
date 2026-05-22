"""3-D displacement decomposition utilities."""

from __future__ import annotations

import numpy as np



def _los_vector(incidence: float, azimuth: float) -> np.ndarray:
    east = -np.sin(incidence) * np.sin(azimuth)
    north = np.sin(incidence) * np.cos(azimuth)
    up = np.cos(incidence)
    return np.array([east, north, up], dtype=float)



def decompose_enu(
    los_measurements: np.ndarray,
    incidence_angles: np.ndarray,
    azimuth_angles: np.ndarray,
    north_prior: np.ndarray | None = None,
    north_sigma: float = 0.01,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Solve for ENU displacement components by weighted least squares."""
    obs = np.asarray(los_measurements, dtype=float)
    if obs.ndim < 2:
        raise ValueError("los_measurements must be at least 2-D with observation axis first")
    n_obs = obs.shape[0]
    incidence = np.asarray(incidence_angles, dtype=float).reshape(n_obs)
    azimuth = np.asarray(azimuth_angles, dtype=float).reshape(n_obs)
    shape = obs.shape[1:]
    A = np.vstack([_los_vector(inc, az) for inc, az in zip(incidence, azimuth)])
    data = obs.reshape(n_obs, -1)
    if north_prior is not None:
        north_prior = np.asarray(north_prior, dtype=float).reshape(1, -1)
        A = np.vstack([A, np.array([0.0, 1.0, 0.0])])
        data = np.vstack([data, north_prior])
        weights = np.diag(np.append(np.ones(n_obs), 1.0 / max(north_sigma, 1.0e-6) ** 2))
    else:
        weights = np.eye(n_obs)
    normal = A.T @ weights @ A
    rhs = A.T @ weights @ data
    solution = np.linalg.pinv(normal) @ rhs
    east, north, up = solution.reshape(3, *shape)
    return east, north, up



def los_to_3d(
    ascending_los: np.ndarray,
    descending_los: np.ndarray,
    asc_incidence: float,
    desc_incidence: float,
    asc_azimuth: float,
    desc_azimuth: float,
    north_prior: np.ndarray | None = None,
    north_sigma: float = 0.01,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Decompose ascending and descending LOS displacement into ENU components."""
    los = np.stack([ascending_los, descending_los])
    incidence = np.array([asc_incidence, desc_incidence], dtype=float)
    azimuth = np.array([asc_azimuth, desc_azimuth], dtype=float)
    return decompose_enu(los, incidence, azimuth, north_prior=north_prior, north_sigma=north_sigma)



def propagate_errors(
    obs_sigma: float | np.ndarray,
    incidence_angles: np.ndarray,
    azimuth_angles: np.ndarray,
    north_sigma: float | None = None,
) -> dict[str, np.ndarray]:
    """Propagate LOS errors into ENU standard deviations."""
    incidence = np.asarray(incidence_angles, dtype=float)
    azimuth = np.asarray(azimuth_angles, dtype=float)
    A = np.vstack([_los_vector(inc, az) for inc, az in zip(incidence, azimuth)])
    sigma = np.asarray(obs_sigma, dtype=float)
    sigma = np.repeat(sigma, len(incidence)) if sigma.ndim == 0 else sigma
    W = np.diag(1.0 / np.maximum(sigma, 1.0e-6) ** 2)
    if north_sigma is not None:
        A = np.vstack([A, np.array([0.0, 1.0, 0.0])])
        W = np.pad(W, ((0, 1), (0, 1)), mode="constant")
        W[-1, -1] = 1.0 / max(north_sigma, 1.0e-6) ** 2
    covariance = np.linalg.pinv(A.T @ W @ A)
    return {"covariance": covariance, "std_enu": np.sqrt(np.diag(covariance))}



def integrate_gps(
    enu_fields: tuple[np.ndarray, np.ndarray, np.ndarray],
    gps_points: np.ndarray,
) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray], dict[str, float]]:
    """Reference ENU fields to GPS using inverse-distance weighted residuals."""
    east, north, up = [np.asarray(component, dtype=float) for component in enu_fields]
    gps = np.asarray(gps_points, dtype=float)
    if gps.ndim != 2 or gps.shape[1] != 5:
        raise ValueError("gps_points must have columns [row, col, east, north, up]")
    rows, cols = east.shape
    yy, xx = np.indices((rows, cols))
    offsets = []
    for field, column in zip((east, north, up), (2, 3, 4)):
        model = np.array([field[int(r), int(c)] for r, c in gps[:, :2]])
        residual = gps[:, column] - model
        dist = np.sqrt((yy[..., None] - gps[:, 0]) ** 2 + (xx[..., None] - gps[:, 1]) ** 2)
        weights = 1.0 / np.maximum(dist, 1.0) ** 2
        offset = np.sum(weights * residual[None, None, :], axis=2) / np.sum(weights, axis=2)
        offsets.append(offset)
    corrected = tuple(field + offset for field, offset in zip((east, north, up), offsets))
    residuals = []
    for field, column in zip(corrected, (2, 3, 4)):
        pred = np.array([field[int(r), int(c)] for r, c in gps[:, :2]])
        residuals.append(gps[:, column] - pred)
    rms = float(np.sqrt(np.mean(np.concatenate(residuals) ** 2)))
    return corrected, {"gps_rms": rms}
