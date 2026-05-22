"""
HEA Descriptor Engineering
Composition-structure-property descriptor calculation for High Entropy Alloys.
Includes: atomic radius mismatch, VEC, mixing entropy/enthalpy, Omega parameter,
electronegativity difference, and CALPHAD-derived phase stability indicators.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
from itertools import combinations


# -----------------------------------------------------------------------
# Element property database (CrMnFeCoNi + common HEA elements)
# Sources: Kittel, Pauling, Materials Project
# -----------------------------------------------------------------------
ELEMENT_PROPS = {
    # elem: {atomic_radius(pm), VEC, electronegativity, melting_point(K),
    #         density(g/cm3), bulk_modulus(GPa), shear_modulus(GPa)}
    "Cr": {"r": 128, "VEC": 6, "chi": 1.66, "Tm": 2180, "rho": 7.19, "B": 190, "G": 115},
    "Mn": {"r": 127, "VEC": 7, "chi": 1.55, "Tm": 1519, "rho": 7.43, "B": 120, "G":  79},
    "Fe": {"r": 126, "VEC": 8, "chi": 1.83, "Tm": 1811, "rho": 7.87, "B": 170, "G":  82},
    "Co": {"r": 125, "VEC": 9, "chi": 1.88, "Tm": 1768, "rho": 8.90, "B": 180, "G":  75},
    "Ni": {"r": 124, "VEC":10, "chi": 1.91, "Tm": 1728, "rho": 8.91, "B": 200, "G":  76},
    "Al": {"r": 143, "VEC": 3, "chi": 1.61, "Tm":  933, "rho": 2.70, "B":  76, "G":  26},
    "Ti": {"r": 147, "VEC": 4, "chi": 1.54, "Tm": 1941, "rho": 4.51, "B": 110, "G":  44},
    "V":  {"r": 134, "VEC": 5, "chi": 1.63, "Tm": 2183, "rho": 6.11, "B": 160, "G":  47},
    "Mo": {"r": 139, "VEC": 6, "chi": 2.16, "Tm": 2896, "rho":10.22, "B": 230, "G": 120},
    "W":  {"r": 139, "VEC": 6, "chi": 2.36, "Tm": 3695, "rho":19.25, "B": 310, "G": 161},
    "Cu": {"r": 128, "VEC":11, "chi": 1.90, "Tm": 1358, "rho": 8.96, "B": 140, "G":  48},
    "Zr": {"r": 160, "VEC": 4, "chi": 1.33, "Tm": 2128, "rho": 6.52, "B":  94, "G":  36},
    "Nb": {"r": 146, "VEC": 5, "chi": 1.60, "Tm": 2750, "rho": 8.57, "B": 170, "G":  38},
    "Hf": {"r": 159, "VEC": 4, "chi": 1.30, "Tm": 2506, "rho":13.31, "B": 110, "G":  56},
    "Ta": {"r": 146, "VEC": 5, "chi": 1.50, "Tm": 3290, "rho":16.65, "B": 200, "G":  69},
}

R = 8.314  # J/(mol·K), gas constant


def _mean_prop(composition: Dict[str, float], prop: str) -> float:
    """Composition-weighted average of element property."""
    return sum(x * ELEMENT_PROPS[el][prop] for el, x in composition.items())


def atomic_radius_mismatch(composition: Dict[str, float]) -> float:
    """
    δ (delta) = sqrt(Σ c_i(1 - r_i/r_bar)^2)
    Quantifies lattice distortion; δ > 6.5% → solid solution destabilization.
    """
    r_bar = _mean_prop(composition, "r")
    delta_sq = sum(x * (1 - ELEMENT_PROPS[el]["r"] / r_bar) ** 2
                   for el, x in composition.items())
    return 100 * np.sqrt(delta_sq)  # percent


def valence_electron_concentration(composition: Dict[str, float]) -> float:
    """
    VEC = Σ c_i * VEC_i
    VEC < 6.87 → BCC phase preferred; VEC > 8.0 → FCC preferred.
    """
    return _mean_prop(composition, "VEC")


def mixing_entropy(composition: Dict[str, float]) -> float:
    """
    ΔS_mix = -R * Σ c_i * ln(c_i)   [J/(mol·K)]
    Maximum for equimolar; drives solid solution stability.
    """
    return -R * sum(x * np.log(x + 1e-12) for x in composition.values())


def mixing_enthalpy(composition: Dict[str, float]) -> float:
    """
    ΔH_mix = Σ_{i≠j} 4 * ΔH_AB^{mix} * c_i * c_j   [kJ/mol]
    Uses Miedema-derived binary interaction parameters.
    Reference: Takeuchi & Inoue (2005), ISIJ Int. 45, 1537.
    """
    # Miedema-derived binary interaction parameters (kJ/mol) - selected pairs
    OMEGA = {
        frozenset(["Cr", "Mn"]): -4.0,  frozenset(["Cr", "Fe"]): -1.0,
        frozenset(["Cr", "Co"]): -4.0,  frozenset(["Cr", "Ni"]): -7.0,
        frozenset(["Mn", "Fe"]): -1.0,  frozenset(["Mn", "Co"]): -5.0,
        frozenset(["Mn", "Ni"]): -8.0,  frozenset(["Fe", "Co"]): -1.0,
        frozenset(["Fe", "Ni"]): -2.0,  frozenset(["Co", "Ni"]): 0.0,
        frozenset(["Al", "Cr"]): -10.0, frozenset(["Al", "Fe"]): -11.0,
        frozenset(["Al", "Co"]): -19.0, frozenset(["Al", "Ni"]): -22.0,
        frozenset(["Al", "Ti"]): -30.0, frozenset(["Ti", "Cr"]): -7.0,
        frozenset(["Ti", "Fe"]): -17.0, frozenset(["Ti", "Ni"]): -35.0,
        frozenset(["Mo", "Ni"]): -7.0,  frozenset(["Mo", "Co"]): -3.0,
        frozenset(["Mo", "Fe"]): -2.0,  frozenset(["Mo", "Cr"]): 0.0,
        frozenset(["W",  "Ni"]): -7.0,  frozenset(["W",  "Fe"]): -7.0,
        frozenset(["V",  "Ni"]): -18.0, frozenset(["V",  "Fe"]): -7.0,
        frozenset(["Nb", "Ni"]): -30.0, frozenset(["Ta", "Ni"]): -29.0,
    }
    elements = list(composition.keys())
    h_mix = 0.0
    for i, j in combinations(range(len(elements)), 2):
        ei, ej = elements[i], elements[j]
        ci, cj = composition[ei], composition[ej]
        key = frozenset([ei, ej])
        omega_ij = OMEGA.get(key, 0.0)  # 0 if not tabulated
        h_mix += 4 * omega_ij * ci * cj
    return h_mix  # kJ/mol


def omega_parameter(composition: Dict[str, float], T: float = 1000.0) -> float:
    """
    Ω = T_m * ΔS_mix / |ΔH_mix|
    Ω > 1.1 → solid solution formation likely.
    """
    Tm_mean = _mean_prop(composition, "Tm")
    dS = mixing_entropy(composition)
    dH = abs(mixing_enthalpy(composition)) + 1e-6  # avoid div by zero
    return Tm_mean * dS / (abs(dH) * 1000)  # convert kJ to J


def electronegativity_mismatch(composition: Dict[str, float]) -> float:
    """Δχ = sqrt(Σ c_i(χ_i - χ_bar)^2)  — Pauling scale."""
    chi_bar = _mean_prop(composition, "chi")
    return np.sqrt(sum(x * (ELEMENT_PROPS[el]["chi"] - chi_bar) ** 2
                       for el, x in composition.items()))


def melting_point_mean(composition: Dict[str, float]) -> float:
    """T_m^avg = Σ c_i * T_m,i  [K]"""
    return _mean_prop(composition, "Tm")


def density_estimate(composition: Dict[str, float]) -> float:
    """Linear mixture rule estimate of density [g/cm³]."""
    return _mean_prop(composition, "rho")


def rule_of_mixture_modulus(composition: Dict[str, float]) -> Tuple[float, float]:
    """Voigt (upper bound) bulk and shear modulus [GPa]."""
    B = _mean_prop(composition, "B")
    G = _mean_prop(composition, "G")
    return B, G


def phase_stability_indicator(composition: Dict[str, float]) -> str:
    """
    Empirical phase selection rules (Zhang et al., 2012; Guo et al., 2011):
    - δ < 6.5 AND |ΔH_mix| < 15 kJ/mol AND -15 < ΔH_mix < 5 → solid solution
    - δ > 6.5 OR |ΔH_mix| > 15 → intermetallic / amorphous tendency
    - VEC < 6.87 → BCC; 6.87 < VEC < 8.0 → mixed; VEC > 8.0 → FCC
    """
    delta = atomic_radius_mismatch(composition)
    dH = mixing_enthalpy(composition)
    vec = valence_electron_concentration(composition)
    omega = omega_parameter(composition)

    phase = []
    if delta < 6.5 and abs(dH) < 15 and omega > 1.1:
        phase.append("SS")  # solid solution
    elif abs(dH) > 20 or delta > 9:
        phase.append("IM")  # intermetallic
    else:
        phase.append("SS+IM")  # mixed

    if vec < 6.87:
        phase.append("BCC")
    elif vec > 8.0:
        phase.append("FCC")
    else:
        phase.append("BCC+FCC")

    return "/".join(phase)


def compute_all_descriptors(composition: Dict[str, float], T: float = 1000.0) -> Dict:
    """
    Compute full descriptor vector for a given composition dict.
    Returns a flat dict suitable for DataFrame row construction.
    """
    B, G = rule_of_mixture_modulus(composition)
    E = 9 * B * G / (3 * B + G)  # Young's modulus (Voigt)
    nu = (3 * B - 2 * G) / (2 * (3 * B + G))  # Poisson's ratio

    desc = {
        "delta_r":     atomic_radius_mismatch(composition),
        "VEC":         valence_electron_concentration(composition),
        "dS_mix":      mixing_entropy(composition),
        "dH_mix":      mixing_enthalpy(composition),
        "Omega":       omega_parameter(composition, T),
        "delta_chi":   electronegativity_mismatch(composition),
        "Tm_mean":     melting_point_mean(composition),
        "density":     density_estimate(composition),
        "B_Voigt":     B,
        "G_Voigt":     G,
        "E_Voigt":     E,
        "nu_Voigt":    nu,
        "n_elements":  len(composition),
        "phase":       phase_stability_indicator(composition),
    }
    # Append individual element fractions
    for el, x in composition.items():
        desc[f"x_{el}"] = x
    return desc


def descriptors_dataframe(compositions: List[Dict[str, float]],
                          T: float = 1000.0) -> pd.DataFrame:
    """Batch computation of descriptors for a list of composition dicts."""
    rows = [compute_all_descriptors(c, T) for c in compositions]
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------
# CALPHAD-inspired Gibbs free energy (simplified regular solution model)
# -----------------------------------------------------------------------
def gibbs_free_energy(composition: Dict[str, float], T: float) -> float:
    """
    G_mix = ΔH_mix - T * ΔS_mix  [kJ/mol]
    Negative G_mix → thermodynamically stable solid solution at temperature T.
    """
    dH = mixing_enthalpy(composition)
    dS = mixing_entropy(composition) / 1000  # J → kJ
    return dH - T * dS


def calphad_phase_diagram_1d(el1: str, el2: str,
                              T_range: np.ndarray,
                              x_range: np.ndarray) -> np.ndarray:
    """
    Simplified binary G_mix landscape over composition and temperature.
    Returns 2D array [T, x] of G_mix values (kJ/mol).
    """
    G = np.zeros((len(T_range), len(x_range)))
    for i, T in enumerate(T_range):
        for j, x in enumerate(x_range):
            comp = {el1: x, el2: 1 - x}
            G[i, j] = gibbs_free_energy(comp, T)
    return G
