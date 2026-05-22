"""
Module 3: Rate-Determining Step (RDS) Identification
=====================================================
Automatic identification via:
  - Campbell's degree of rate control (X_RC)
  - Sensitivity analysis
  - Apparent activation energy decomposition
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Tuple
from .rate_constants import RateConstant


@dataclass
class RDSResult:
    """Result of RDS analysis."""
    step_labels: List[str]
    X_RC: np.ndarray                # Degree of rate control for each step
    X_TRC: np.ndarray               # Thermodynamic rate control
    rds_label: str                  # Identified RDS
    rds_index: int                  # Index of RDS
    sensitivity_coefficients: Dict  # Sensitivity analysis results


def degree_of_rate_control(rate_constants: List[RateConstant],
                           rate_function,
                           coverages: np.ndarray,
                           conditions: dict,
                           delta: float = 0.01) -> RDSResult:
    """
    Calculate Campbell's degree of rate control.

    X_RC,i = (k_i / r) * (∂r / ∂k_i)_{K_eq,i, k_j≠i}

    Parameters
    ----------
    rate_constants : list of RateConstant
        All elementary step rate constants.
    rate_function : callable
        Function(rate_constants, coverages, conditions) -> overall_rate.
    coverages : np.ndarray
        Current surface coverages.
    conditions : dict
        Reaction conditions (T, P, etc.).
    delta : float
        Fractional perturbation for numerical derivative.

    Returns
    -------
    RDSResult
        Degree of rate control for each step.
    """
    n_steps = len(rate_constants)

    # Base rate
    r_base = rate_function(rate_constants, coverages, conditions)

    X_RC = np.zeros(n_steps)
    X_TRC = np.zeros(n_steps)
    sensitivity = {}

    for i in range(n_steps):
        # Perturb forward rate constant (keeping K_eq constant)
        rc_perturbed = [rc for rc in rate_constants]  # shallow copy
        original = rate_constants[i]

        # Forward perturbation: multiply k_fwd and k_rev by (1+delta) to keep K_eq
        perturbed_fwd = RateConstant(
            label=original.label,
            k_forward=original.k_forward * (1 + delta),
            k_reverse=original.k_reverse * (1 + delta),
            K_eq=original.K_eq,
            E_act_forward=original.E_act_forward,
            E_act_reverse=original.E_act_reverse,
            tunneling_correction=original.tunneling_correction,
            temperature=original.temperature
        )
        rc_perturbed_list = list(rate_constants)
        rc_perturbed_list[i] = perturbed_fwd

        r_perturbed = rate_function(rc_perturbed_list, coverages, conditions)

        # X_RC = (k/r) * (dr/dk)
        if abs(r_base) > 1e-30:
            X_RC[i] = (original.k_forward / r_base) * \
                       (r_perturbed - r_base) / (original.k_forward * delta)
        else:
            X_RC[i] = 0.0

        # Thermodynamic rate control (perturb K_eq only)
        perturbed_keq = RateConstant(
            label=original.label,
            k_forward=original.k_forward * (1 + delta),
            k_reverse=original.k_reverse,
            K_eq=original.K_eq * (1 + delta),
            E_act_forward=original.E_act_forward,
            E_act_reverse=original.E_act_reverse,
            tunneling_correction=original.tunneling_correction,
            temperature=original.temperature
        )
        rc_trc_list = list(rate_constants)
        rc_trc_list[i] = perturbed_keq
        r_trc = rate_function(rc_trc_list, coverages, conditions)

        if abs(r_base) > 1e-30:
            X_TRC[i] = (original.K_eq / r_base) * \
                        (r_trc - r_base) / (original.K_eq * delta)
        else:
            X_TRC[i] = 0.0

        sensitivity[rate_constants[i].label] = {
            'X_RC': X_RC[i],
            'X_TRC': X_TRC[i],
            'k_forward': original.k_forward,
            'E_act': original.E_act_forward
        }

    rds_idx = int(np.argmax(np.abs(X_RC)))
    labels = [rc.label for rc in rate_constants]

    return RDSResult(
        step_labels=labels,
        X_RC=X_RC,
        X_TRC=X_TRC,
        rds_label=labels[rds_idx],
        rds_index=rds_idx,
        sensitivity_coefficients=sensitivity
    )


def apparent_activation_energy(rate_function, rate_constants_func,
                               coverages: np.ndarray, conditions: dict,
                               T_center: float, delta_T: float = 5.0) -> float:
    """
    Calculate apparent activation energy from temperature derivative.

    E_app = -R * d(ln r) / d(1/T)

    Parameters
    ----------
    rate_function : callable
        Function(rate_constants, coverages, conditions) -> rate.
    rate_constants_func : callable
        Function(T) -> list of RateConstant.
    coverages : np.ndarray
        Current coverages.
    conditions : dict
        Reaction conditions.
    T_center : float
        Central temperature [K].
    delta_T : float
        Temperature perturbation [K].

    Returns
    -------
    float
        Apparent activation energy [eV].
    """
    EV_TO_J = 1.602176634e-19
    R = 8.314462618

    T1 = T_center - delta_T
    T2 = T_center + delta_T

    cond1 = dict(conditions, T=T1)
    cond2 = dict(conditions, T=T2)

    rc1 = rate_constants_func(T1)
    rc2 = rate_constants_func(T2)

    r1 = rate_function(rc1, coverages, cond1)
    r2 = rate_function(rc2, coverages, cond2)

    if r1 <= 0 or r2 <= 0:
        return 0.0

    ln_r1 = np.log(r1)
    ln_r2 = np.log(r2)

    E_app = -R * (ln_r2 - ln_r1) / (1.0 / T2 - 1.0 / T1)  # J/mol
    return E_app / (EV_TO_J * 6.022e23)  # eV


def energy_span_analysis(rate_constants: List[RateConstant],
                         intermediate_energies: List[float]) -> Dict:
    """
    Kozuch-Shaik Energy Span Model for catalytic cycle analysis.

    δE = T_TDTS - I_TDI + ΔG_r (if TDTS appears after TDI)
    δE = T_TDTS - I_TDI       (if TDTS appears before TDI)

    Parameters
    ----------
    rate_constants : list of RateConstant
        Elementary step rate constants.
    intermediate_energies : list of float
        Free energies of intermediates [eV].

    Returns
    -------
    dict
        Energy span analysis results.
    """
    n = len(rate_constants)
    ts_energies = []
    for rc in rate_constants:
        # TS energy = intermediate energy + activation energy
        idx = min(rate_constants.index(rc), len(intermediate_energies) - 1)
        ts_energies.append(intermediate_energies[idx] + rc.E_act_forward)

    # Find TDI (most stable intermediate) and TDTS (highest TS)
    tdi_idx = int(np.argmin(intermediate_energies))
    tdts_idx = int(np.argmax(ts_energies))

    tdi_energy = intermediate_energies[tdi_idx]
    tdts_energy = ts_energies[tdts_idx]

    # Reaction energy (overall)
    delta_G_r = intermediate_energies[-1] - intermediate_energies[0]

    # Energy span
    if tdts_idx >= tdi_idx:
        energy_span = tdts_energy - tdi_energy
    else:
        energy_span = tdts_energy - tdi_energy + delta_G_r

    return {
        'energy_span': energy_span,
        'TDTS_index': tdts_idx,
        'TDTS_label': rate_constants[tdts_idx].label,
        'TDTS_energy': tdts_energy,
        'TDI_index': tdi_idx,
        'TDI_energy': tdi_energy,
        'delta_G_reaction': delta_G_r,
        'TOF_estimate': (1.380649e-23 * 500 / 6.626e-34) *
                        np.exp(-energy_span * 1.602e-19 / (1.380649e-23 * 500))
    }
