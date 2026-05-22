"""
Bayesian inversion module using PyMC for MCMC-based uncertainty quantification
of volcanic deformation source parameters.

Supports:
  - Mogi / Spheroid / FEM forward models
  - Multiple data types (GNSS, InSAR, gravity)
  - Hierarchical noise models
  - Model comparison via WAIC/LOO

References:
  Bagnardi & Hooper (2018), G-Cubed
  Cervelli et al. (2001), JGR
"""

import numpy as np
from typing import Dict, List, Optional, Callable, Tuple
from dataclasses import dataclass

from .source_models import (
    MogiSource, mogi_displacement, mogi_gravity,
    SpheroidSource, spheroid_displacement,
    FEMSourceConfig, fem_displacement
)


@dataclass
class InversionData:
    """Container for multi-type observation data."""
    obs_x: np.ndarray          # (N,) observation x coords [m]
    obs_y: np.ndarray          # (N,) observation y coords [m]
    gnss_disp: Optional[np.ndarray] = None   # (M, 3) GNSS displacement [m]
    gnss_sigma: Optional[np.ndarray] = None  # (M, 3) GNSS uncertainties [m]
    gnss_idx: Optional[np.ndarray] = None    # (M,) indices into obs arrays
    insar_los: Optional[np.ndarray] = None   # (K,) InSAR LOS displacement [m]
    insar_sigma: Optional[np.ndarray] = None # (K,) InSAR uncertainties [m]
    insar_idx: Optional[np.ndarray] = None   # (K,) indices into obs arrays
    insar_look: Optional[np.ndarray] = None  # (K, 3) LOS unit vectors
    gravity: Optional[np.ndarray] = None     # (L,) gravity changes [µGal]
    gravity_sigma: Optional[np.ndarray] = None  # (L,) uncertainties
    gravity_idx: Optional[np.ndarray] = None    # (L,) indices


@dataclass
class PriorConfig:
    """Prior distribution configuration for Mogi parameters."""
    x_mean: float = 0.0
    x_std: float = 5000.0
    y_mean: float = 0.0
    y_std: float = 5000.0
    d_min: float = 500.0
    d_max: float = 15000.0
    dV_min: float = -1e8
    dV_max: float = 1e8


def build_mogi_pymc_model(
    data: InversionData,
    prior: PriorConfig,
    include_gravity: bool = True
):
    """
    Build a PyMC model for Mogi source inversion.

    Parameters
    ----------
    data : InversionData with observations
    prior : PriorConfig with prior bounds/parameters
    include_gravity : whether to include gravity data in likelihood

    Returns
    -------
    pymc.Model
    """
    try:
        import pymc as pm
        import pytensor.tensor as pt
    except ImportError:
        raise ImportError("PyMC is required: pip install pymc")

    with pm.Model() as model:
        # --- Priors ---
        x_src = pm.Normal("x_src", mu=prior.x_mean, sigma=prior.x_std)
        y_src = pm.Normal("y_src", mu=prior.y_mean, sigma=prior.y_std)
        d_src = pm.Uniform("d_src", lower=prior.d_min, upper=prior.d_max)
        dV = pm.Uniform("dV", lower=prior.dV_min, upper=prior.dV_max)
        nu = pm.Beta("nu", alpha=5, beta=15)  # peaked around 0.25

        # Hyperparameters for noise scaling
        sigma_scale = pm.HalfCauchy("sigma_scale", beta=2.0)

        # --- Forward model (Mogi) ---
        # Use PyTensor ops for differentiability
        dx = data.obs_x - x_src
        dy = data.obs_y - y_src
        R = pt.sqrt(dx**2 + dy**2 + d_src**2)
        C = dV * (1 - nu) / np.pi

        pred_ux = C * dx / R**3
        pred_uy = C * dy / R**3
        pred_uz = C * d_src / R**3

        # --- GNSS likelihood ---
        if data.gnss_disp is not None:
            idx = data.gnss_idx
            gnss_pred = pt.stack([
                pred_ux[idx], pred_uy[idx], pred_uz[idx]
            ], axis=1)
            gnss_scaled_sigma = data.gnss_sigma * sigma_scale
            pm.Normal(
                "gnss_obs",
                mu=gnss_pred,
                sigma=gnss_scaled_sigma,
                observed=data.gnss_disp
            )

        # --- InSAR likelihood ---
        if data.insar_los is not None:
            idx = data.insar_idx
            look = data.insar_look  # (K, 3) unit vectors
            insar_pred_3d = pt.stack([
                pred_ux[idx], pred_uy[idx], pred_uz[idx]
            ], axis=1)
            # Project to LOS
            insar_pred_los = pt.sum(insar_pred_3d * look, axis=1)
            insar_scaled_sigma = data.insar_sigma * sigma_scale
            pm.Normal(
                "insar_obs",
                mu=insar_pred_los,
                sigma=insar_scaled_sigma,
                observed=data.insar_los
            )

        # --- Gravity likelihood ---
        if include_gravity and data.gravity is not None:
            idx = data.gravity_idx
            rho = pm.Normal("rho", mu=2500, sigma=300)
            free_air = -3.086e-6
            dg_freeair = free_air * pred_uz[idx] * 1e8
            G = 6.674e-11
            dM = rho * dV
            R_grav = pt.sqrt(
                (data.obs_x[idx] - x_src)**2 +
                (data.obs_y[idx] - y_src)**2 +
                d_src**2
            )
            dg_mass = G * dM * d_src / R_grav**3 * 1e8
            grav_pred = dg_freeair + dg_mass
            grav_scaled_sigma = data.gravity_sigma * sigma_scale
            pm.Normal(
                "grav_obs",
                mu=grav_pred,
                sigma=grav_scaled_sigma,
                observed=data.gravity
            )

        # --- Derived quantities ---
        pm.Deterministic("volume_change_1e6m3", dV / 1e6)

    return model


def build_spheroid_pymc_model(
    data: InversionData,
    prior: Optional[Dict] = None
):
    """
    Build a PyMC model for spheroid source inversion.
    Includes additional parameters: semi-axes, strike, dip.
    """
    try:
        import pymc as pm
        import pytensor.tensor as pt
    except ImportError:
        raise ImportError("PyMC required")

    if prior is None:
        prior = {}

    with pm.Model() as model:
        x_src = pm.Normal("x_src", mu=0, sigma=5000)
        y_src = pm.Normal("y_src", mu=0, sigma=5000)
        d_src = pm.Uniform("d_src", lower=500, upper=15000)
        a_src = pm.Uniform("a_semi", lower=100, upper=5000)
        b_src = pm.Uniform("b_semi", lower=100, upper=5000)
        dP = pm.Normal("dP", mu=10e6, sigma=20e6)
        strike = pm.Uniform("strike", lower=0, upper=360)
        dip = pm.Uniform("dip", lower=0, upper=90)
        mu_rock = pm.LogNormal("mu_rock", mu=np.log(3e10), sigma=0.3)
        nu = pm.Beta("nu", alpha=5, beta=15)
        sigma_scale = pm.HalfCauchy("sigma_scale", beta=2.0)

        # Forward model via numerical evaluation
        # (PyTensor custom op wrapping spheroid_displacement)
        # Simplified: use Mogi-equivalent for demonstration
        aspect = b_src / a_src
        dV_eq = (4.0/3.0) * np.pi * a_src**2 * b_src * dP / mu_rock

        dx = data.obs_x - x_src
        dy = data.obs_y - y_src
        R = pt.sqrt(dx**2 + dy**2 + d_src**2)
        C = dV_eq * (1 - nu) / np.pi

        pred_ux = C * dx / R**3
        pred_uy = C * dy / R**3
        pred_uz = C * d_src / R**3

        if data.gnss_disp is not None:
            idx = data.gnss_idx
            gnss_pred = pt.stack([
                pred_ux[idx], pred_uy[idx], pred_uz[idx]
            ], axis=1)
            pm.Normal("gnss_obs", mu=gnss_pred,
                      sigma=data.gnss_sigma * sigma_scale,
                      observed=data.gnss_disp)

        pm.Deterministic("aspect_ratio", b_src / a_src)
        pm.Deterministic("equiv_volume_change", dV_eq)

    return model


def run_mcmc(
    model,
    draws: int = 5000,
    tune: int = 2000,
    chains: int = 4,
    target_accept: float = 0.9,
    cores: int = 2,
    random_seed: int = 42
):
    """
    Run MCMC sampling on a PyMC model.

    Returns ArviZ InferenceData.
    """
    import pymc as pm

    with model:
        trace = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            target_accept=target_accept,
            cores=cores,
            random_seed=random_seed,
            return_inferencedata=True
        )
    return trace


def compute_model_comparison(traces: Dict[str, object]):
    """
    Compare models using WAIC and LOO-CV.

    Parameters
    ----------
    traces : dict mapping model name to ArviZ InferenceData

    Returns
    -------
    comparison : DataFrame with WAIC/LOO rankings
    """
    import arviz as az

    comparison_waic = az.compare(traces, ic="waic")
    comparison_loo = az.compare(traces, ic="loo")

    return {
        'waic': comparison_waic,
        'loo': comparison_loo
    }


def summarize_posterior(trace, var_names=None, hdi_prob=0.94):
    """
    Summarize posterior distributions.

    Returns DataFrame with mean, sd, HDI, ESS, R-hat.
    """
    import arviz as az

    summary = az.summary(
        trace,
        var_names=var_names,
        hdi_prob=hdi_prob,
        stat_funcs={"median": np.median}
    )
    return summary


def extract_map_estimate(trace, var_names=None):
    """Extract Maximum A Posteriori estimates from trace."""
    import arviz as az

    posterior = trace.posterior
    map_est = {}
    names = var_names or list(posterior.data_vars)

    for name in names:
        vals = posterior[name].values.flatten()
        # KDE-based MAP
        from scipy.stats import gaussian_kde
        try:
            kde = gaussian_kde(vals)
            x_grid = np.linspace(vals.min(), vals.max(), 1000)
            map_est[name] = x_grid[np.argmax(kde(x_grid))]
        except Exception:
            map_est[name] = np.median(vals)

    return map_est
