"""
Module 1: DFT-Based Rate Constant Calculation
==============================================
Transition State Theory (TST) with tunneling corrections.
Supports Wigner and Eckart tunneling approximations.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional

# Physical constants
KB = 1.380649e-23       # Boltzmann constant [J/K]
H_PLANCK = 6.62607015e-34  # Planck constant [J·s]
R_GAS = 8.314462618     # Gas constant [J/(mol·K)]
EV_TO_J = 1.602176634e-19  # eV to Joule
KCAL_TO_J = 4184.0      # kcal/mol to J/mol


@dataclass
class TransitionState:
    """Represents a transition state from DFT calculations."""
    label: str
    E_activation: float          # Activation energy [eV]
    E_reaction: float            # Reaction energy [eV]
    frequencies_real: list        # Real vibrational frequencies [cm^-1]
    frequency_imaginary: float    # Imaginary frequency at TS [cm^-1]
    symmetry_number: int = 1
    spin_multiplicity: int = 1


@dataclass
class RateConstant:
    """Calculated rate constant with metadata."""
    label: str
    k_forward: float              # Forward rate constant [1/s]
    k_reverse: float              # Reverse rate constant [1/s]
    K_eq: float                   # Equilibrium constant
    E_act_forward: float          # Forward activation energy [eV]
    E_act_reverse: float          # Reverse activation energy [eV]
    tunneling_correction: float   # Tunneling correction factor
    temperature: float            # Temperature [K]


def partition_function_vibration(frequencies_cm: list, T: float) -> float:
    """
    Calculate vibrational partition function (quantum harmonic oscillator).

    Parameters
    ----------
    frequencies_cm : list of float
        Vibrational frequencies in cm^-1.
    T : float
        Temperature in K.

    Returns
    -------
    float
        Vibrational partition function.
    """
    q_vib = 1.0
    for nu in frequencies_cm:
        if nu <= 0:
            continue
        # Convert cm^-1 to energy in J
        hv = H_PLANCK * nu * 2.998e10  # h * c * nu_tilde
        x = hv / (KB * T)
        q_vib *= 1.0 / (1.0 - np.exp(-x))
    return q_vib


def zero_point_energy(frequencies_cm: list) -> float:
    """Calculate zero-point energy from vibrational frequencies [eV]."""
    zpe = 0.0
    for nu in frequencies_cm:
        if nu <= 0:
            continue
        hv = H_PLANCK * nu * 2.998e10
        zpe += 0.5 * hv
    return zpe / EV_TO_J


def wigner_tunneling(frequency_imaginary: float, T: float) -> float:
    """
    Wigner tunneling correction factor.

    κ = 1 + (1/24) * (hν‡ / k_B T)^2

    Parameters
    ----------
    frequency_imaginary : float
        Imaginary frequency at the transition state [cm^-1].
    T : float
        Temperature [K].

    Returns
    -------
    float
        Tunneling correction factor κ.
    """
    hv = H_PLANCK * abs(frequency_imaginary) * 2.998e10
    x = hv / (KB * T)
    kappa = 1.0 + (1.0 / 24.0) * x**2
    return kappa


def eckart_tunneling(E_forward: float, E_reverse: float,
                     frequency_imaginary: float, T: float,
                     n_points: int = 1000) -> float:
    """
    Eckart tunneling correction using numerical integration.

    Parameters
    ----------
    E_forward : float
        Forward barrier height [eV].
    E_reverse : float
        Reverse barrier height [eV].
    frequency_imaginary : float
        Imaginary frequency [cm^-1].
    T : float
        Temperature [K].
    n_points : int
        Number of integration points.

    Returns
    -------
    float
        Eckart tunneling correction factor κ.
    """
    V1 = E_forward * EV_TO_J * 6.022e23  # J/mol
    V2 = E_reverse * EV_TO_J * 6.022e23

    if V1 <= 0 or V2 <= 0:
        return 1.0

    nu_im = abs(frequency_imaginary) * 2.998e10  # Hz
    if nu_im == 0:
        return 1.0

    # Eckart barrier parameters
    alpha1 = 2 * np.pi * V1 / (H_PLANCK * nu_im * 6.022e23)
    alpha2 = 2 * np.pi * V2 / (H_PLANCK * nu_im * 6.022e23)

    # Numerical integration of transmission probability
    E_max = max(V1, V2) * 3
    E_arr = np.linspace(0.001 * V1, E_max, n_points)

    beta = 1.0 / (R_GAS * T)
    integrand = np.zeros(n_points)

    for i, E in enumerate(E_arr):
        xi = E / V1
        # Simplified Eckart transmission coefficient
        a = np.sqrt(alpha1 * xi)
        b = np.sqrt(alpha2 * (xi - 1 + V2 / V1)) if xi > (1 - V2 / V1) else 0.0

        if a + b > 50:
            T_E = 1.0
        elif a > 0:
            cosh_sum = np.cosh(2 * (a + b)) if (a + b) < 350 else np.exp(2 * (a + b)) / 2
            cosh_diff = np.cosh(2 * (a - b)) if abs(a - b) < 350 else np.exp(2 * abs(a - b)) / 2
            d = np.cosh(2 * np.pi * np.sqrt(max(alpha1 * alpha2 - 0.25, 0)))
            T_E = 1.0 - (cosh_sum - cosh_diff) / (cosh_sum + d) if (cosh_sum + d) > 0 else 0.0
        else:
            T_E = 0.0

        T_E = np.clip(T_E, 0, 1)
        integrand[i] = T_E * np.exp(-beta * E)

    # Classical transmission
    classical = np.exp(-beta * V1) / beta if beta * V1 < 700 else 1e-300

    kappa = np.trapz(integrand, E_arr) * beta / max(classical * beta, 1e-300)
    return max(kappa, 1.0)


def calculate_tst_rate(ts: TransitionState, T: float,
                       tunneling: str = "wigner") -> RateConstant:
    """
    Calculate rate constant using Transition State Theory.

    k_TST = κ * (k_B T / h) * (Q‡ / Q_R) * exp(-E_a / k_B T)

    Parameters
    ----------
    ts : TransitionState
        Transition state information from DFT.
    T : float
        Temperature [K].
    tunneling : str
        Tunneling correction method: "none", "wigner", or "eckart".

    Returns
    -------
    RateConstant
        Calculated forward and reverse rate constants.
    """
    # Prefactor: k_B T / h
    prefactor = KB * T / H_PLANCK

    # Forward activation energy
    E_act_fwd = ts.E_activation  # eV

    # Reverse activation energy
    E_act_rev = ts.E_activation - ts.E_reaction  # eV

    # Tunneling correction
    if tunneling == "wigner":
        kappa = wigner_tunneling(ts.frequency_imaginary, T)
    elif tunneling == "eckart":
        kappa = eckart_tunneling(E_act_fwd, E_act_rev, ts.frequency_imaginary, T)
    else:
        kappa = 1.0

    # Vibrational partition function ratio (simplified)
    q_ts = partition_function_vibration(ts.frequencies_real, T)

    # TST rate constant (surface reaction, per site per second)
    k_fwd = (kappa * prefactor / ts.symmetry_number *
             np.exp(-E_act_fwd * EV_TO_J / (KB * T)))

    k_rev = (kappa * prefactor / ts.symmetry_number *
             np.exp(-E_act_rev * EV_TO_J / (KB * T)))

    # Equilibrium constant
    K_eq = k_fwd / k_rev if k_rev > 0 else np.inf

    return RateConstant(
        label=ts.label,
        k_forward=k_fwd,
        k_reverse=k_rev,
        K_eq=K_eq,
        E_act_forward=E_act_fwd,
        E_act_reverse=E_act_rev,
        tunneling_correction=kappa,
        temperature=T
    )


def arrhenius_parameters(ts: TransitionState,
                         T_range: tuple = (300, 1000),
                         n_points: int = 50,
                         tunneling: str = "wigner") -> dict:
    """
    Extract Arrhenius parameters (A, E_a) from TST rates over temperature range.

    Returns
    -------
    dict
        'A_forward', 'Ea_forward', 'A_reverse', 'Ea_reverse',
        'temperatures', 'k_forward', 'k_reverse'
    """
    temps = np.linspace(T_range[0], T_range[1], n_points)
    k_fwd = np.zeros(n_points)
    k_rev = np.zeros(n_points)

    for i, T in enumerate(temps):
        rc = calculate_tst_rate(ts, T, tunneling)
        k_fwd[i] = rc.k_forward
        k_rev[i] = rc.k_reverse

    # Linear fit: ln(k) = ln(A) - Ea/(R*T)
    inv_T = 1.0 / temps
    ln_k_fwd = np.log(np.maximum(k_fwd, 1e-300))
    ln_k_rev = np.log(np.maximum(k_rev, 1e-300))

    coeff_fwd = np.polyfit(inv_T, ln_k_fwd, 1)
    coeff_rev = np.polyfit(inv_T, ln_k_rev, 1)

    Ea_fwd = -coeff_fwd[0] * R_GAS / (EV_TO_J * 6.022e23)  # eV
    A_fwd = np.exp(coeff_fwd[1])
    Ea_rev = -coeff_rev[0] * R_GAS / (EV_TO_J * 6.022e23)
    A_rev = np.exp(coeff_rev[1])

    return {
        'A_forward': A_fwd, 'Ea_forward': Ea_fwd,
        'A_reverse': A_rev, 'Ea_reverse': Ea_rev,
        'temperatures': temps, 'k_forward': k_fwd, 'k_reverse': k_rev
    }
