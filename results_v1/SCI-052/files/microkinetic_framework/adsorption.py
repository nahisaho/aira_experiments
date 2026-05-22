"""
Module 2: Adsorption Isotherm Models
=====================================
Langmuir, Temkin, and Fractal surface adsorption isotherms.
"""

import numpy as np
from dataclasses import dataclass
from typing import Callable


@dataclass
class AdsorptionParameters:
    """Parameters for adsorption models."""
    species: str
    delta_H_ads: float    # Adsorption enthalpy [eV]
    delta_S_ads: float    # Adsorption entropy [eV/K]
    site_type: str = "top"
    sticking_coefficient: float = 1.0


# --- Physical Constants ---
KB_EV = 8.617333e-5  # Boltzmann constant [eV/K]
KB = 1.380649e-23
H_PLANCK = 6.62607015e-34
R_GAS = 8.314462618
EV_TO_J = 1.602176634e-19


def langmuir_isotherm(P: float, T: float, params: AdsorptionParameters) -> float:
    """
    Langmuir adsorption isotherm.

    θ = K·P / (1 + K·P)
    K = exp(-ΔG_ads / k_B T)

    Parameters
    ----------
    P : float
        Partial pressure [bar].
    T : float
        Temperature [K].
    params : AdsorptionParameters

    Returns
    -------
    float
        Surface coverage θ (0 to 1).
    """
    delta_G = params.delta_H_ads - T * params.delta_S_ads
    K = np.exp(-delta_G / (KB_EV * T))
    theta = K * P / (1.0 + K * P)
    return np.clip(theta, 0, 1)


def competitive_langmuir(pressures: dict, T: float,
                         params_dict: dict) -> dict:
    """
    Competitive Langmuir isotherm for multiple species.

    θ_i = K_i·P_i / (1 + Σ K_j·P_j)

    Parameters
    ----------
    pressures : dict
        {species_name: partial_pressure [bar]}
    T : float
        Temperature [K].
    params_dict : dict
        {species_name: AdsorptionParameters}

    Returns
    -------
    dict
        {species_name: coverage}
    """
    K_values = {}
    for sp, par in params_dict.items():
        delta_G = par.delta_H_ads - T * par.delta_S_ads
        K_values[sp] = np.exp(-delta_G / (KB_EV * T))

    denom = 1.0 + sum(K_values[sp] * pressures.get(sp, 0) for sp in K_values)
    coverages = {sp: K_values[sp] * pressures.get(sp, 0) / denom
                 for sp in K_values}
    return coverages


def temkin_isotherm(P: float, T: float, params: AdsorptionParameters,
                    alpha: float = 0.5, delta_E: float = 0.3) -> float:
    """
    Temkin adsorption isotherm (linear decrease in adsorption energy with coverage).

    E_ads(θ) = E_ads_0 - α·δE·θ

    Parameters
    ----------
    P : float
        Partial pressure [bar].
    T : float
        Temperature [K].
    params : AdsorptionParameters
    alpha : float
        Temkin parameter (0 to 1).
    delta_E : float
        Range of adsorption energy variation [eV].

    Returns
    -------
    float
        Surface coverage θ.
    """
    # Self-consistent solution via iteration
    theta = 0.5  # initial guess
    for _ in range(100):
        E_eff = params.delta_H_ads + alpha * delta_E * theta
        delta_G = E_eff - T * params.delta_S_ads
        K = np.exp(-delta_G / (KB_EV * T))
        theta_new = K * P / (1.0 + K * P)
        if abs(theta_new - theta) < 1e-10:
            break
        theta = 0.5 * theta + 0.5 * theta_new
    return np.clip(theta, 0, 1)


def fractal_isotherm(P: float, T: float, params: AdsorptionParameters,
                     D_f: float = 2.5, E_min: float = -1.5,
                     E_max: float = -0.5, n_sites: int = 100) -> float:
    """
    Fractal surface adsorption isotherm.

    Integrates Langmuir isotherm over a distribution of adsorption energies
    determined by the fractal dimension D_f.

    f(E) ∝ E^(D_f - 3)

    Parameters
    ----------
    P : float
        Partial pressure [bar].
    T : float
        Temperature [K].
    params : AdsorptionParameters
    D_f : float
        Fractal dimension (2 ≤ D_f ≤ 3).
    E_min, E_max : float
        Range of adsorption energies [eV].
    n_sites : int
        Number of integration points.

    Returns
    -------
    float
        Average surface coverage θ.
    """
    D_f = np.clip(D_f, 2.0, 3.0)
    energies = np.linspace(E_min, E_max, n_sites)

    # Energy distribution for fractal surface
    weights = np.abs(energies) ** (D_f - 3.0)
    weights /= np.sum(weights)

    theta_avg = 0.0
    for E, w in zip(energies, weights):
        delta_G = E - T * params.delta_S_ads
        K = np.exp(-delta_G / (KB_EV * T))
        theta_local = K * P / (1.0 + K * P)
        theta_avg += w * theta_local

    return np.clip(theta_avg, 0, 1)


def sticking_rate(P: float, T: float, mass_amu: float,
                  S0: float = 1.0, site_density: float = 1.5e19) -> float:
    """
    Hertz-Knudsen impingement rate with sticking coefficient.

    r_ads = S0 * P / √(2π m k_B T) / n_sites

    Parameters
    ----------
    P : float
        Partial pressure [Pa].
    T : float
        Temperature [K].
    mass_amu : float
        Molecular mass [amu].
    S0 : float
        Sticking coefficient.
    site_density : float
        Surface site density [sites/m^2].

    Returns
    -------
    float
        Adsorption rate [1/(site·s)].
    """
    mass_kg = mass_amu * 1.66054e-27
    flux = P / np.sqrt(2.0 * np.pi * mass_kg * KB * T)
    return S0 * flux / site_density
