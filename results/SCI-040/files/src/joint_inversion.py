"""
Joint inversion of GNSS + InSAR + gravity data.

Handles:
  - Multi-type data weighting (VCE: Variance Component Estimation)
  - InSAR orbital ramp correction
  - Spatial covariance for InSAR
  - Cross-validation for weight optimization

References:
  Sambridge (1999), GJI
  Fukuda & Johnson (2008), GJI
"""

import numpy as np
from typing import Dict, Optional, Tuple, List
from dataclasses import dataclass

from .source_models import (
    MogiSource, mogi_displacement, mogi_gravity,
    SpheroidSource, spheroid_displacement
)
from .bayesian_inversion import InversionData


@dataclass
class JointInversionConfig:
    """Configuration for joint inversion."""
    # Relative weights for each data type
    w_gnss: float = 1.0
    w_insar: float = 1.0
    w_gravity: float = 1.0
    # InSAR covariance parameters
    insar_cov_model: str = "exponential"  # exponential / gaussian
    insar_cov_range: float = 5000.0       # correlation length [m]
    insar_cov_sill: float = 1e-6          # variance [m^2]
    # Orbital ramp correction
    insar_remove_ramp: bool = True
    ramp_order: int = 1                   # 1=linear, 2=quadratic
    # VCE iterations
    vce_iterations: int = 5
    vce_method: str = "helmert"  # helmert / minque


def build_covariance_matrix(
    x: np.ndarray,
    y: np.ndarray,
    model: str = "exponential",
    sill: float = 1e-6,
    range_param: float = 5000.0,
    nugget: float = 1e-8
) -> np.ndarray:
    """
    Build spatial covariance matrix for InSAR data.

    Parameters
    ----------
    x, y : (N,) coordinates
    model : covariance model type
    sill, range_param, nugget : model parameters

    Returns
    -------
    C : (N, N) covariance matrix
    """
    N = len(x)
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    dist = np.sqrt(dx**2 + dy**2)

    if model == "exponential":
        C = sill * np.exp(-dist / range_param)
    elif model == "gaussian":
        C = sill * np.exp(-0.5 * (dist / range_param)**2)
    elif model == "spherical":
        h = dist / range_param
        C = sill * np.where(h <= 1, 1 - 1.5*h + 0.5*h**3, 0)
    else:
        raise ValueError(f"Unknown covariance model: {model}")

    C += nugget * np.eye(N)
    return C


def remove_orbital_ramp(
    x: np.ndarray,
    y: np.ndarray,
    los: np.ndarray,
    order: int = 1
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Remove orbital ramp from InSAR LOS displacements.

    Parameters
    ----------
    x, y : coordinates
    los : LOS displacements
    order : polynomial order (1=linear, 2=quadratic)

    Returns
    -------
    los_corrected : ramp-removed LOS
    ramp_coeffs : fitted ramp coefficients
    """
    if order == 1:
        G = np.column_stack([np.ones_like(x), x, y])
    elif order == 2:
        G = np.column_stack([np.ones_like(x), x, y, x**2, x*y, y**2])
    else:
        raise ValueError("order must be 1 or 2")

    coeffs, _, _, _ = np.linalg.lstsq(G, los, rcond=None)
    ramp = G @ coeffs
    return los - ramp, coeffs


def variance_component_estimation(
    residuals_list: List[np.ndarray],
    sigma_list: List[np.ndarray],
    weights: List[float],
    n_params: int,
    method: str = "helmert",
    n_iter: int = 5
) -> List[float]:
    """
    Variance Component Estimation (VCE) to optimally weight
    multiple data types.

    Parameters
    ----------
    residuals_list : list of residual vectors per data type
    sigma_list : list of a priori uncertainty vectors
    weights : initial weights
    n_params : number of estimated parameters
    method : VCE method
    n_iter : number of iterations

    Returns
    -------
    updated_weights : optimized relative weights
    """
    n_types = len(residuals_list)
    w = np.array(weights, dtype=float)

    for iteration in range(n_iter):
        for i in range(n_types):
            v = residuals_list[i]
            sig = sigma_list[i]
            n_i = len(v)

            if method == "helmert":
                # Helmert VCE
                chi2_i = np.sum((v / sig)**2)
                redundancy_i = max(n_i - n_params / n_types, 1)
                w[i] = redundancy_i / chi2_i
            elif method == "minque":
                # MINQUE
                P_i = np.diag(1.0 / sig**2)
                w[i] = (v @ P_i @ v) / np.trace(P_i)

        # Normalize
        w = w / np.sum(w) * n_types

    return w.tolist()


def build_joint_design_matrix(
    data: InversionData,
    source_type: str = "mogi",
    params: Dict = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build the joint design matrix (Jacobian) for linearized inversion.

    Uses finite differences to compute partial derivatives.

    Parameters
    ----------
    data : observation data
    source_type : "mogi" or "spheroid"
    params : current parameter estimates

    Returns
    -------
    G : (N_total, N_params) design matrix
    d : (N_total,) data vector
    W : (N_total, N_total) weight matrix
    """
    if params is None:
        params = {'x': 0, 'y': 0, 'd': 5000, 'dV': 1e6}

    param_names = list(params.keys())
    n_params = len(param_names)

    # Collect data and predictions
    d_vec = []
    w_vec = []
    pred_func_list = []

    def predict_all(p):
        """Predict all data types given parameter dict."""
        preds = []
        src = MogiSource(x=p['x'], y=p['y'], d=p['d'], dV=p['dV'])

        if data.gnss_disp is not None:
            disp = mogi_displacement(
                data.obs_x[data.gnss_idx],
                data.obs_y[data.gnss_idx], src
            )
            preds.append(disp.flatten())

        if data.insar_los is not None:
            disp = mogi_displacement(
                data.obs_x[data.insar_idx],
                data.obs_y[data.insar_idx], src
            )
            los_pred = np.sum(disp * data.insar_look, axis=1)
            preds.append(los_pred)

        if data.gravity is not None:
            grav = mogi_gravity(
                data.obs_x[data.gravity_idx],
                data.obs_y[data.gravity_idx], src
            )
            preds.append(grav)

        return np.concatenate(preds)

    # Data vector
    if data.gnss_disp is not None:
        d_vec.append(data.gnss_disp.flatten())
        w_vec.append(1.0 / data.gnss_sigma.flatten()**2)
    if data.insar_los is not None:
        d_vec.append(data.insar_los)
        w_vec.append(1.0 / data.insar_sigma**2)
    if data.gravity is not None:
        d_vec.append(data.gravity)
        w_vec.append(1.0 / data.gravity_sigma**2)

    d_obs = np.concatenate(d_vec)
    w_diag = np.concatenate(w_vec)
    W = np.diag(w_diag)

    # Jacobian via finite differences
    n_data = len(d_obs)
    G = np.zeros((n_data, n_params))
    pred_0 = predict_all(params)

    for j, pname in enumerate(param_names):
        dp = dict(params)
        h = max(abs(params[pname]) * 1e-4, 1.0)
        dp[pname] = params[pname] + h
        pred_h = predict_all(dp)
        G[:, j] = (pred_h - pred_0) / h

    return G, d_obs, W


def joint_least_squares(
    G: np.ndarray,
    d: np.ndarray,
    W: np.ndarray,
    regularization: float = 0.0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Weighted least-squares solution with optional Tikhonov regularization.

    Returns
    -------
    m : parameter estimates
    Cm : posterior covariance matrix
    """
    GtW = G.T @ W
    N = GtW @ G + max(regularization, 1e-10) * np.eye(G.shape[1])
    rhs = GtW @ d

    m = np.linalg.solve(N, rhs)
    try:
        Cm = np.linalg.inv(N)
    except np.linalg.LinAlgError:
        Cm = np.linalg.pinv(N)

    return m, Cm


def iterative_joint_inversion(
    data: InversionData,
    config: JointInversionConfig,
    initial_params: Dict,
    max_iter: int = 20,
    convergence_tol: float = 1e-6
) -> Dict:
    """
    Iterative linearized joint inversion with VCE weight updates.

    Parameters
    ----------
    data : multi-type observation data
    config : inversion configuration
    initial_params : starting parameter estimates
    max_iter : maximum iterations
    convergence_tol : relative parameter change for convergence

    Returns
    -------
    result : dict with 'params', 'covariance', 'residuals', 'weights', 'iterations'
    """
    params = dict(initial_params)
    weights = [config.w_gnss, config.w_insar, config.w_gravity]

    history = []

    for it in range(max_iter):
        G, d_obs, W = build_joint_design_matrix(data, params=params)
        m, Cm = joint_least_squares(G, d_obs, W)

        # Update parameters
        param_names = list(params.keys())
        new_params = {}
        for j, pname in enumerate(param_names):
            new_params[pname] = params[pname] + m[j]

        # Check convergence
        rel_change = np.max([
            abs(new_params[k] - params[k]) / max(abs(params[k]), 1e-10)
            for k in param_names
        ])

        history.append({
            'iteration': it,
            'params': dict(new_params),
            'rel_change': rel_change
        })

        params = new_params

        if rel_change < convergence_tol:
            break

        # VCE weight update (every 3 iterations)
        if (it + 1) % 3 == 0 and config.vce_iterations > 0:
            pred = G @ m
            residuals_by_type = _split_residuals(d_obs - pred, data)
            sigmas_by_type = _split_sigmas(data)
            n_p = len(param_names)
            weights = variance_component_estimation(
                residuals_by_type, sigmas_by_type, weights, n_p,
                method=config.vce_method, n_iter=config.vce_iterations
            )

    return {
        'params': params,
        'covariance': Cm,
        'weights': weights,
        'iterations': len(history),
        'history': history,
        'converged': rel_change < convergence_tol
    }


def _split_residuals(residuals: np.ndarray, data: InversionData) -> List[np.ndarray]:
    """Split combined residual vector by data type."""
    parts = []
    idx = 0
    if data.gnss_disp is not None:
        n = data.gnss_disp.size
        parts.append(residuals[idx:idx+n])
        idx += n
    if data.insar_los is not None:
        n = len(data.insar_los)
        parts.append(residuals[idx:idx+n])
        idx += n
    if data.gravity is not None:
        n = len(data.gravity)
        parts.append(residuals[idx:idx+n])
        idx += n
    return parts


def _split_sigmas(data: InversionData) -> List[np.ndarray]:
    """Extract sigma arrays per data type."""
    parts = []
    if data.gnss_sigma is not None:
        parts.append(data.gnss_sigma.flatten())
    if data.insar_sigma is not None:
        parts.append(data.insar_sigma)
    if data.gravity_sigma is not None:
        parts.append(data.gravity_sigma)
    return parts
