"""
Module 4: Coverage-Dependent Lateral Interactions
==================================================
Models for adsorbate-adsorbate interactions on catalyst surfaces.
Supports:
  - Mean-field lateral interaction model
  - Pairwise interaction (nearest-neighbor)
  - Cluster expansion approach
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional


@dataclass
class LateralInteractionParams:
    """Parameters for lateral interaction model."""
    species_pair: Tuple[str, str]
    epsilon_nn: float     # Nearest-neighbor interaction energy [eV] (negative = repulsive)
    epsilon_nnn: float = 0.0   # Next-nearest-neighbor [eV]
    z_nn: int = 4         # Nearest-neighbor coordination number
    z_nnn: int = 4        # Next-nearest-neighbor coordination number


@dataclass
class CoverageState:
    """Represents the current surface coverage state."""
    species: List[str]
    coverages: np.ndarray      # Coverage for each species
    vacant_fraction: float     # Fraction of vacant sites
    total_energy_shift: float  # Total energy shift from lateral interactions [eV]


def mean_field_interaction_energy(coverages: Dict[str, float],
                                  interactions: List[LateralInteractionParams]) -> Dict[str, float]:
    """
    Calculate coverage-dependent energy shifts using mean-field approximation.

    ΔE_ads,i(θ) = Σ_j ε_{ij} * z * θ_j

    Parameters
    ----------
    coverages : dict
        {species: coverage}
    interactions : list of LateralInteractionParams
        Pairwise interaction parameters.

    Returns
    -------
    dict
        {species: energy_shift [eV]}
    """
    energy_shifts = {sp: 0.0 for sp in coverages}

    for inter in interactions:
        sp1, sp2 = inter.species_pair
        if sp1 in coverages and sp2 in coverages:
            # NN contribution
            energy_shifts[sp1] += inter.epsilon_nn * inter.z_nn * coverages[sp2]
            if sp1 != sp2:
                energy_shifts[sp2] += inter.epsilon_nn * inter.z_nn * coverages[sp1]
            # NNN contribution
            energy_shifts[sp1] += inter.epsilon_nnn * inter.z_nnn * coverages[sp2]
            if sp1 != sp2:
                energy_shifts[sp2] += inter.epsilon_nnn * inter.z_nnn * coverages[sp1]

    return energy_shifts


def modified_rate_constants(k_forward: float, k_reverse: float,
                            delta_E_reactant: float, delta_E_ts: float,
                            delta_E_product: float, T: float,
                            alpha_bep: float = 0.5) -> Tuple[float, float]:
    """
    Modify rate constants for coverage-dependent lateral interactions.

    Uses BEP (Brønsted-Evans-Polanyi) relation for TS energy shift:
    ΔE_TS = α * ΔE_product + (1-α) * ΔE_reactant

    Parameters
    ----------
    k_forward, k_reverse : float
        Original rate constants.
    delta_E_reactant, delta_E_ts, delta_E_product : float
        Energy shifts due to lateral interactions [eV].
    T : float
        Temperature [K].
    alpha_bep : float
        BEP coefficient (0 to 1).

    Returns
    -------
    tuple
        (k_forward_modified, k_reverse_modified)
    """
    KB_EV = 8.617333e-5

    # TS energy shift from BEP if not directly provided
    if delta_E_ts == 0.0:
        delta_E_ts = alpha_bep * delta_E_product + (1 - alpha_bep) * delta_E_reactant

    # Forward: change in activation barrier
    delta_Ea_fwd = delta_E_ts - delta_E_reactant
    # Reverse: change in activation barrier
    delta_Ea_rev = delta_E_ts - delta_E_product

    k_fwd_mod = k_forward * np.exp(-delta_Ea_fwd / (KB_EV * T))
    k_rev_mod = k_reverse * np.exp(-delta_Ea_rev / (KB_EV * T))

    return k_fwd_mod, k_rev_mod


def quasi_chemical_approximation(theta: float, epsilon: float,
                                 z: int, T: float) -> float:
    """
    Quasi-chemical approximation for lateral interactions.

    Gives a better estimate than mean-field for strong interactions.

    Parameters
    ----------
    theta : float
        Surface coverage.
    epsilon : float
        Pairwise interaction energy [eV].
    z : int
        Coordination number.
    T : float
        Temperature [K].

    Returns
    -------
    float
        Effective chemical potential shift [eV].
    """
    KB_EV = 8.617333e-5
    beta = 1.0 / (KB_EV * T)

    if abs(epsilon) < 1e-10:
        return 0.0

    # QCA parameter
    eta = np.exp(-beta * epsilon) - 1.0
    discriminant = (1.0 - 2.0 * theta)**2 + 4.0 * theta * (1.0 - theta) * np.exp(-beta * epsilon)

    if discriminant < 0:
        discriminant = 0.0

    sqrt_disc = np.sqrt(discriminant)
    if abs(1.0 - 2.0 * theta + sqrt_disc) < 1e-30:
        return z * epsilon * theta  # fall back to mean-field

    mu_shift = (z / 2.0) * KB_EV * T * np.log(
        (1.0 - 2.0 * theta + sqrt_disc) / (2.0 * (1.0 - theta))
    )
    return mu_shift


def solve_coverage_self_consistent(pressures: Dict[str, float],
                                   T: float,
                                   adsorption_energies: Dict[str, float],
                                   entropy_contributions: Dict[str, float],
                                   interactions: List[LateralInteractionParams],
                                   max_iter: int = 500,
                                   tol: float = 1e-8) -> Dict[str, float]:
    """
    Self-consistent solution of coverages with lateral interactions.

    Iteratively solves for coverages considering mutual interaction effects.

    Parameters
    ----------
    pressures : dict
        {species: partial pressure [bar]}
    T : float
        Temperature [K].
    adsorption_energies : dict
        {species: E_ads [eV]}
    entropy_contributions : dict
        {species: S_ads [eV/K]}
    interactions : list of LateralInteractionParams
    max_iter : int
        Maximum iterations.
    tol : float
        Convergence tolerance.

    Returns
    -------
    dict
        {species: converged coverage}
    """
    KB_EV = 8.617333e-5
    species = list(pressures.keys())
    n = len(species)

    # Initial guess
    coverages = {sp: 0.1 for sp in species}

    for iteration in range(max_iter):
        old_coverages = dict(coverages)

        # Calculate lateral interaction energy shifts
        energy_shifts = mean_field_interaction_energy(coverages, interactions)

        # Update coverages
        denom = 1.0
        K_values = {}
        for sp in species:
            E_eff = adsorption_energies[sp] + energy_shifts.get(sp, 0.0)
            delta_G = E_eff - T * entropy_contributions.get(sp, 0.0)
            K = np.exp(-delta_G / (KB_EV * T))
            K_values[sp] = K
            denom += K * pressures[sp]

        for sp in species:
            coverages[sp] = K_values[sp] * pressures[sp] / denom

        # Check convergence
        max_change = max(abs(coverages[sp] - old_coverages[sp]) for sp in species)
        if max_change < tol:
            break

        # Damping for stability
        damping = 0.3
        for sp in species:
            coverages[sp] = damping * coverages[sp] + (1 - damping) * old_coverages[sp]

    return coverages
