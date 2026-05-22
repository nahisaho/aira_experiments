"""
Gas-Particle Partitioning Thermodynamic Model
Implements:
  - Pankow absorptive partitioning theory
  - Volatility Basis Set (VBS) framework
  - Simplified UNIFAC activity coefficients
  - AIOMFAC-style electrolyte corrections
"""
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import logging

logger = logging.getLogger(__name__)

R_GAS = 8.314        # J/(mol·K)
T_REF = 298.15       # K
MW_OM  = 200.0       # g/mol assumed mean OA molecular weight
AVOGADRO = 6.022e23


@dataclass
class PartitioningResult:
    species: str
    Psat: float          # Pa
    Cstar: float         # μg/m3  effective saturation conc.
    xi_gas: float        # gas-phase mole fraction
    xi_part: float       # particle-phase mole fraction
    Fpart: float         # particle-phase mass fraction [0-1]
    Kp: float            # partitioning coefficient (m3/μg)
    gamma: float         # activity coefficient (UNIFAC)
    delta_H_vap: float   # kJ/mol  enthalpy of vaporization


# ── Estimated thermodynamic data for key SOA products ────────────────────────
# Psat in Pa @ 298 K,  dH_vap in kJ/mol
SPECIES_THERMO = {
    "pinic_acid":        {"Psat": 1.2e-4,  "dHvap": 90.0, "MW": 186.2, "kappa": 0.10},
    "pinonic_acid":      {"Psat": 0.072,   "dHvap": 85.0, "MW": 184.2, "kappa": 0.08},
    "norpinic_acid":     {"Psat": 2.0e-5,  "dHvap": 95.0, "MW": 172.2, "kappa": 0.12},
    "pinanediol":        {"Psat": 0.040,   "dHvap": 78.0, "MW": 172.3, "kappa": 0.06},
    "alpha_pin_OH":      {"Psat": 1.20,    "dHvap": 65.0, "MW": 152.2, "kappa": 0.04},
    "alpha_pin_nitrate": {"Psat": 0.01,    "dHvap": 88.0, "MW": 215.2, "kappa": 0.05},
    "pinaldehyde":       {"Psat": 9.0,     "dHvap": 60.0, "MW": 152.2, "kappa": 0.03},
    "norpinaldehyde":    {"Psat": 15.0,    "dHvap": 58.0, "MW": 138.2, "kappa": 0.02},
    "nopinaldehyde":     {"Psat": 22.0,    "dHvap": 56.0, "MW": 138.2, "kappa": 0.02},
    "nopinone":          {"Psat": 18.0,    "dHvap": 55.0, "MW": 138.2, "kappa": 0.02},
    "limonene_OH":       {"Psat": 0.80,    "dHvap": 72.0, "MW": 152.2, "kappa": 0.05},
    "limonic_acid":      {"Psat": 5.0e-5,  "dHvap": 92.0, "MW": 186.2, "kappa": 0.11},
    "limonaketone":      {"Psat": 0.30,    "dHvap": 75.0, "MW": 170.2, "kappa": 0.06},
    "limonene_nitrate":  {"Psat": 0.03,    "dHvap": 87.0, "MW": 215.2, "kappa": 0.05},
    "LIMAL":             {"Psat": 2.0,     "dHvap": 65.0, "MW": 168.2, "kappa": 0.04},
    "ISOPOOH":           {"Psat": 4.0,     "dHvap": 58.0, "MW": 120.1, "kappa": 0.05},
    "methacrolein":      {"Psat": 9050.0,  "dHvap": 38.0, "MW": 70.09, "kappa": 0.01},
    "MVK":               {"Psat": 12300.0, "dHvap": 36.0, "MW": 70.09, "kappa": 0.01},
    "2MGA":              {"Psat": 0.50,    "dHvap": 80.0, "MW": 120.1, "kappa": 0.09},
    "delta_ISOP_NO3":    {"Psat": 0.12,    "dHvap": 84.0, "MW": 147.1, "kappa": 0.07},
    "beta_ISOP_NO3":     {"Psat": 0.12,    "dHvap": 84.0, "MW": 147.1, "kappa": 0.07},
    "cresol":            {"Psat": 165.0,   "dHvap": 52.0, "MW": 108.1, "kappa": 0.03},
    "toluene_RO2":       {"Psat": 0.05,    "dHvap": 85.0, "MW": 189.1, "kappa": 0.07},
    "DHBO":              {"Psat": 8.0,     "dHvap": 62.0, "MW": 124.1, "kappa": 0.04},
    "beta_pin_OH":       {"Psat": 2.5,     "dHvap": 63.0, "MW": 152.2, "kappa": 0.04},
    "beta_pin_nitrate":  {"Psat": 0.02,    "dHvap": 89.0, "MW": 215.2, "kappa": 0.05},
    "cis_pinonic":       {"Psat": 0.072,   "dHvap": 85.0, "MW": 184.2, "kappa": 0.08},
    "7_OH_lim":          {"Psat": 0.15,    "dHvap": 77.0, "MW": 170.2, "kappa": 0.06},
    "pinic_acid":        {"Psat": 1.2e-4,  "dHvap": 90.0, "MW": 186.2, "kappa": 0.10},
    "methylglyoxal":     {"Psat": 3700.0,  "dHvap": 40.0, "MW": 72.06, "kappa": 0.06},
    "formaldehyde":      {"Psat": 1.8e5,   "dHvap": 23.3, "MW": 30.03, "kappa": 0.01},
    "acetone":           {"Psat": 30800.0, "dHvap": 31.0, "MW": 58.08, "kappa": 0.01},
    "ring_frag_products":{"Psat": 0.20,    "dHvap": 75.0, "MW": 160.0, "kappa": 0.06},
    "methylnitrophenol": {"Psat": 0.10,    "dHvap": 80.0, "MW": 153.1, "kappa": 0.08},
    "benzaldehyde":      {"Psat": 170.0,   "dHvap": 50.0, "MW": 106.1, "kappa": 0.02},
}


# ── UNIFAC group interaction parameters (simplified subset) ──────────────────
# Groups relevant for SOA: CH2(1), C=C(2), OH(3), CHO(4), C=O(5), COOH(6), ONO2(7)
UNIFAC_R = {1: 0.6744, 2: 1.3454, 3: 1.0000, 4: 0.9980, 5: 1.6724, 6: 1.3013, 7: 2.0000}
UNIFAC_Q = {1: 0.5400, 2: 1.1760, 3: 1.2000, 4: 0.9480, 5: 1.4880, 6: 1.2240, 7: 1.8000}

# Interaction parameters a_mn [K]  (asymmetric)
UNIFAC_A = np.array([
    # 1       2       3       4       5       6       7
    [   0.0,  -35.36, 986.5,  677.0,  476.4,  663.5, 1200.0],  # 1 CH2
    [  35.36,   0.0,  693.9,  505.7,  182.6,  318.9,  900.0],  # 2 C=C
    [-137.1, -229.1,   0.0,  -137.1,  -203.6, -199.0, -100.0],  # 3 OH
    [ 505.7,  782.0,  529.0,   0.0,  -103.6,  -22.97, 400.0],  # 4 CHO
    [ 164.5,  237.7,  -101.7, 135.9,    0.0,  -200.7, 300.0],  # 5 C=O
    [ 315.3,  349.2,   -66.17,-103.6, -235.7,    0.0,  200.0],  # 6 COOH
    [ 600.0,  700.0,   200.0,  300.0,  250.0,  150.0,    0.0],  # 7 ONO2
])

# Simplified group compositions for key SOA species
SPECIES_GROUPS = {
    "pinic_acid":       {1: 3, 6: 2},
    "pinonic_acid":     {1: 5, 5: 1, 6: 1},
    "norpinic_acid":    {1: 2, 6: 2},
    "limonene_OH":      {1: 4, 2: 1, 3: 1},
    "limonic_acid":     {1: 3, 6: 2},
    "ISOPOOH":          {1: 2, 2: 1, 3: 2},
    "2MGA":             {1: 1, 3: 1, 6: 1},
    "delta_ISOP_NO3":   {1: 2, 2: 1, 7: 1},
    "toluene_RO2":      {5: 2, 3: 1, 2: 3},
    "cresol":           {3: 1, 1: 1, 2: 3},
    "alpha_pin_nitrate":{1: 4, 7: 1},
}


def calc_unifac_gamma(species: str, T: float = 298.15, x_w: float = 0.0) -> float:
    """
    Simplified UNIFAC activity coefficient calculation.
    Returns γ_i (activity coefficient of species i in OA mixture).
    Uses combinatorial + residual contributions.
    """
    if species not in SPECIES_GROUPS:
        # Default: assume near-ideal (γ ≈ 1.0-1.5)
        return 1.0 + 0.3 * x_w  # slight non-ideality with water

    groups = SPECIES_GROUPS[species]
    # Simplified: combinatorial term dominates for dilute organics
    r_i = sum(UNIFAC_R[g] * n for g, n in groups.items())
    q_i = sum(UNIFAC_Q[g] * n for g, n in groups.items())

    # Reference molecule (OA mixture average): 8 CH2 groups
    r_ref, q_ref = 8 * UNIFAC_R[1], 8 * UNIFAC_Q[1]

    phi = r_i / (r_i + r_ref)
    theta = q_i / (q_i + q_ref)

    ln_gamma_C = np.log(phi) + 1 - phi - 5 * q_i * (np.log(phi / theta) + 1 - phi / theta)

    # Residual term (simplified): temperature-dependent interaction
    ln_gamma_R = 0.0
    for g, n in groups.items():
        for h in groups:
            tau = np.exp(-UNIFAC_A[g - 1, h - 1] / T)
            ln_gamma_R += n * q_i * (1 - np.log(tau + 1e-10) - tau / (tau + 0.5))

    ln_gamma = ln_gamma_C + np.clip(ln_gamma_R * 0.05, -2, 2)
    return float(np.exp(ln_gamma))


def calc_psat_T(psat_ref: float, dHvap: float, T: float, T_ref: float = 298.15) -> float:
    """Clausius-Clapeyron: Psat(T) from Psat(T_ref)."""
    return psat_ref * np.exp(dHvap * 1e3 / R_GAS * (1 / T_ref - 1 / T))


def pankow_partitioning(
    species: str,
    Coa: float,          # μg/m3 OA concentration
    T: float = 298.15,
    RH: float = 0.5,
) -> PartitioningResult:
    """
    Pankow (1994) absorptive partitioning with UNIFAC activity coefficients.
    Kp = (RT) / (MW_OM * gamma * Psat * 10^6)   [m3/μg]
    Fpart = Coa * Kp / (1 + Coa * Kp)
    """
    thermo = SPECIES_THERMO.get(species, {"Psat": 10.0, "dHvap": 70.0, "MW": 200.0, "kappa": 0.05})
    Psat_ref = thermo["Psat"]
    dHvap    = thermo["dHvap"]
    MW_i     = thermo["MW"]

    Psat = calc_psat_T(Psat_ref, dHvap, T)
    gamma = calc_unifac_gamma(species, T, x_w=RH * 0.3)

    # Cstar [μg/m3] — effective saturation concentration
    Cstar = (MW_i * Psat) / (R_GAS * T) * 1e6  # convert mol/m3 → μg/m3

    # Kp [m3/μg]
    Kp = (R_GAS * T) / (MW_OM * 1e-3 * gamma * Psat * 1e6)

    # AIOMFAC correction for aqueous fraction (simplified)
    kappa = thermo.get("kappa", 0.05)
    f_aq = RH * kappa / (kappa + (1 - RH))   # simplified kappa-Köhler
    Kp_eff = Kp * (1 + f_aq * 2.0)           # enhanced partitioning with water

    Fpart = (Coa * Kp_eff) / (1.0 + Coa * Kp_eff)

    return PartitioningResult(
        species=species,
        Psat=Psat,
        Cstar=Cstar,
        xi_gas=1 - Fpart,
        xi_part=Fpart,
        Fpart=Fpart,
        Kp=Kp_eff,
        gamma=gamma,
        delta_H_vap=dHvap,
    )


def run_partitioning_ensemble(
    species_list: List[str],
    Coa: float = 10.0,      # μg/m3
    T: float = 298.15,
    RH: float = 0.5,
) -> List[PartitioningResult]:
    """Batch partitioning calculation for all SOA species."""
    results = []
    for sp in species_list:
        res = pankow_partitioning(sp, Coa, T, RH)
        results.append(res)
    return results


# ── Volatility Basis Set (VBS) ────────────────────────────────────────────────
VBS_BINS = [-3, -2, -1, 0, 1, 2, 3]  # log10(C* / μg m-3)

def assign_vbs_bins(results: List[PartitioningResult]) -> Dict[str, List[str]]:
    """Assign species to VBS volatility bins."""
    bins: Dict[str, List[str]] = {str(b): [] for b in VBS_BINS}
    for r in results:
        log_cstar = np.log10(max(r.Cstar, 1e-5))
        bin_key = str(min(VBS_BINS, key=lambda b: abs(b - log_cstar)))
        bins[bin_key].append(r.species)
    return bins


def temperature_sensitivity_partitioning(
    species: str,
    T_range: np.ndarray,
    Coa: float = 10.0,
    RH: float = 0.5,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Fpart over a temperature range."""
    Fparts = np.array([
        pankow_partitioning(species, Coa, T, RH).Fpart
        for T in T_range
    ])
    return T_range, Fparts
