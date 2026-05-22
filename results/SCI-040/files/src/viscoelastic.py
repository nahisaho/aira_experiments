"""
Viscoelastic crustal response correction for volcanic deformation.

Implements:
  - Standard Linear Solid (SLS) rheology
  - Maxwell rheology
  - Burgers body rheology
  - Time-dependent Green's functions for viscoelastic half-space
  - Correction factors for elastic inversion

References:
  Segall (2010), "Earthquake and Volcano Deformation"
  Del Negro et al. (2009), GJI
  Newman et al. (2006), JGR
"""

import numpy as np
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass

from .source_models import MogiSource, mogi_displacement


@dataclass
class RheologyParams:
    """Viscoelastic rheology parameters."""
    model: str = "sls"        # sls / maxwell / burgers
    mu_elastic: float = 3e10  # Elastic shear modulus [Pa]
    mu_kelvin: float = 1e10   # Kelvin element modulus [Pa] (SLS/Burgers)
    eta_maxwell: float = 1e18 # Maxwell viscosity [Pa·s]
    eta_kelvin: float = 1e17  # Kelvin viscosity [Pa·s] (SLS/Burgers)
    nu: float = 0.25          # Poisson's ratio

    @property
    def tau_maxwell(self) -> float:
        """Maxwell relaxation time [s]."""
        return self.eta_maxwell / self.mu_elastic

    @property
    def tau_kelvin(self) -> float:
        """Kelvin relaxation time [s]."""
        return self.eta_kelvin / self.mu_kelvin

    @property
    def tau_maxwell_days(self) -> float:
        return self.tau_maxwell / 86400.0

    @property
    def tau_kelvin_days(self) -> float:
        return self.tau_kelvin / 86400.0


def sls_relaxation_function(t: np.ndarray, rheology: RheologyParams) -> np.ndarray:
    """
    Standard Linear Solid relaxation function.

    G(t) = mu_elastic * [1 + (mu_e/mu_k) * exp(-t/tau_k)]

    Parameters
    ----------
    t : time array [seconds]
    rheology : rheology parameters

    Returns
    -------
    G : relaxation modulus as a function of time
    """
    mu_e = rheology.mu_elastic
    mu_k = rheology.mu_kelvin
    tau_k = rheology.tau_kelvin

    return mu_e * (1 + (mu_e / mu_k) * np.exp(-t / tau_k))


def maxwell_relaxation_function(t: np.ndarray, rheology: RheologyParams) -> np.ndarray:
    """
    Maxwell body relaxation function.

    G(t) = mu * exp(-t/tau_m)
    """
    return rheology.mu_elastic * np.exp(-t / rheology.tau_maxwell)


def burgers_relaxation_function(t: np.ndarray, rheology: RheologyParams) -> np.ndarray:
    """
    Burgers body relaxation function (Maxwell + Kelvin in series).

    J(t) = 1/mu_e + t/eta_m + (1/mu_k)*(1 - exp(-t/tau_k))
    G(t) ≈ 1/J(t)   (approximate inverse of creep function)
    """
    mu_e = rheology.mu_elastic
    mu_k = rheology.mu_kelvin
    eta_m = rheology.eta_maxwell
    tau_k = rheology.tau_kelvin

    J = 1.0/mu_e + t/eta_m + (1.0/mu_k) * (1 - np.exp(-t/tau_k))
    return 1.0 / J


def viscoelastic_correction_factor(
    t: np.ndarray,
    rheology: RheologyParams
) -> np.ndarray:
    """
    Compute time-dependent correction factor for converting
    elastic Mogi solution to viscoelastic response.

    The correction factor C(t) satisfies:
      u_ve(t) = C(t) * u_elastic

    For a Mogi source in a viscoelastic half-space:
      C(t) = 1 + A * (1 - exp(-t/tau))  (simplified)

    where A depends on the rheology model.

    Parameters
    ----------
    t : time array [seconds]
    rheology : rheology parameters

    Returns
    -------
    C : (len(t),) correction factors
    """
    if rheology.model == "maxwell":
        # Maxwell: amplification grows with time
        tau = rheology.tau_maxwell
        C = 1 + (1 - rheology.nu) * t / tau
        # Cap at reasonable value
        C = np.minimum(C, 5.0)
        return C

    elif rheology.model == "sls":
        # SLS: bounded amplification
        tau_k = rheology.tau_kelvin
        mu_ratio = rheology.mu_elastic / rheology.mu_kelvin
        A = mu_ratio / (1 + mu_ratio)
        C = 1 + A * (1 - np.exp(-t / tau_k))
        return C

    elif rheology.model == "burgers":
        # Burgers: combination of transient and steady-state
        tau_m = rheology.tau_maxwell
        tau_k = rheology.tau_kelvin
        mu_ratio = rheology.mu_elastic / rheology.mu_kelvin
        A_k = mu_ratio / (1 + mu_ratio)
        C = 1 + A_k * (1 - np.exp(-t / tau_k)) + (1 - rheology.nu) * t / tau_m
        C = np.minimum(C, 10.0)
        return C

    else:
        raise ValueError(f"Unknown rheology model: {rheology.model}")


def correct_elastic_displacement(
    disp_elastic: np.ndarray,
    t: float,
    rheology: RheologyParams
) -> np.ndarray:
    """
    Apply viscoelastic correction to elastic displacement predictions.

    Parameters
    ----------
    disp_elastic : (N, 3) elastic displacement
    t : time since pressure onset [seconds]
    rheology : rheology parameters

    Returns
    -------
    disp_ve : (N, 3) viscoelastically corrected displacement
    """
    C = viscoelastic_correction_factor(np.array([t]), rheology)[0]
    return disp_elastic * C


def mogi_viscoelastic_displacement(
    obs_x: np.ndarray,
    obs_y: np.ndarray,
    source: MogiSource,
    t: float,
    rheology: RheologyParams
) -> np.ndarray:
    """
    Mogi displacement with viscoelastic correction.

    Parameters
    ----------
    obs_x, obs_y : observation coordinates
    source : Mogi source parameters
    t : time since pressurization [seconds]
    rheology : rheology parameters

    Returns
    -------
    disp : (N, 3) viscoelastic displacement
    """
    disp_el = mogi_displacement(obs_x, obs_y, source)
    return correct_elastic_displacement(disp_el, t, rheology)


def compute_viscoelastic_timeseries(
    obs_x: np.ndarray,
    obs_y: np.ndarray,
    source: MogiSource,
    times: np.ndarray,
    rheology: RheologyParams,
    pressure_history: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Compute displacement time series with viscoelastic effects.

    For time-varying pressure, uses convolution with relaxation function.

    Parameters
    ----------
    obs_x, obs_y : observation coordinates
    source : Mogi source (reference state)
    times : time array [seconds]
    rheology : rheology parameters
    pressure_history : (len(times),) normalized pressure history
                      (default: step function at t=0)

    Returns
    -------
    disp_ts : (len(times), N, 3) displacement time series
    """
    N_obs = len(obs_x)
    N_t = len(times)
    disp_ts = np.zeros((N_t, N_obs, 3))

    if pressure_history is None:
        # Step function
        for i, t in enumerate(times):
            if t >= 0:
                disp_ts[i] = mogi_viscoelastic_displacement(
                    obs_x, obs_y, source, t, rheology
                )
    else:
        # Convolution with pressure history
        disp_el = mogi_displacement(obs_x, obs_y, source)

        for i, t in enumerate(times):
            # Discrete convolution
            C_total = 0.0
            for j in range(i + 1):
                dt = times[i] - times[j]
                if dt >= 0:
                    dP = pressure_history[j]
                    if j > 0:
                        dP -= pressure_history[j-1]
                    C_j = viscoelastic_correction_factor(
                        np.array([dt]), rheology
                    )[0]
                    C_total += dP * C_j

            disp_ts[i] = disp_el * C_total

    return disp_ts


def estimate_rheology_from_postseismic(
    obs_disp_ts: np.ndarray,
    times: np.ndarray,
    source: MogiSource,
    obs_x: np.ndarray,
    obs_y: np.ndarray,
    model: str = "sls"
) -> Dict:
    """
    Estimate rheological parameters from post-eruption relaxation data.

    Uses nonlinear least squares to fit relaxation curve.

    Parameters
    ----------
    obs_disp_ts : (N_t, N_obs, 3) observed displacement time series
    times : time array [seconds]
    source : assumed Mogi source
    obs_x, obs_y : station coordinates
    model : rheology model

    Returns
    -------
    result : dict with fitted parameters and residuals
    """
    from scipy.optimize import minimize

    disp_el = mogi_displacement(obs_x, obs_y, source)

    def misfit(params):
        if model == "sls":
            log_tau_k, log_mu_ratio = params
            tau_k = 10**log_tau_k
            mu_ratio = 10**log_mu_ratio
            rheology = RheologyParams(
                model="sls",
                mu_kelvin=rheology.mu_elastic / mu_ratio,
                eta_kelvin=rheology.mu_elastic / mu_ratio * tau_k
            )
        elif model == "maxwell":
            log_tau_m = params[0]
            tau_m = 10**log_tau_m
            rheology = RheologyParams(
                model="maxwell",
                eta_maxwell=3e10 * tau_m
            )

        C = viscoelastic_correction_factor(times, rheology)
        pred = np.array([disp_el * c for c in C])
        return np.sum((pred - obs_disp_ts)**2)

    if model == "sls":
        x0 = [np.log10(86400 * 365), 0]  # 1 year, ratio=1
        bounds = [(np.log10(86400), np.log10(86400*365*100)),
                  (-2, 2)]
    elif model == "maxwell":
        x0 = [np.log10(86400 * 365)]
        bounds = [(np.log10(86400), np.log10(86400*365*100))]

    result = minimize(misfit, x0, method='L-BFGS-B', bounds=bounds)

    return {
        'params': result.x,
        'success': result.success,
        'misfit': result.fun,
        'model': model
    }
