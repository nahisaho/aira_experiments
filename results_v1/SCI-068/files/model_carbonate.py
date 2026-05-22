"""
Module 1: Seawater CO2 Chemical Equilibrium (Carbonate System)
CO2SYS-based numerical calculation of the marine carbonate system.
"""

import numpy as np
import json
from scipy.optimize import fsolve

# === Constants (Lueker et al. 2000, Dickson & Millero 1987) ===
# Dissociation constants for carbonic acid in seawater (mol/kg-SW)

def compute_K1(T_K, S):
    """First dissociation constant of carbonic acid (Lueker et al. 2000)."""
    pK1 = (3633.86 / T_K - 61.2172 + 9.6777 * np.log(T_K)
            - 0.011555 * S + 0.0001152 * S**2)
    return 10**(-pK1)

def compute_K2(T_K, S):
    """Second dissociation constant of carbonic acid (Lueker et al. 2000)."""
    pK2 = (471.78 / T_K + 25.929 - 3.16967 * np.log(T_K)
            - 0.01781 * S + 0.0001122 * S**2)
    return 10**(-pK2)

def compute_Ksp_aragonite(T_K, S):
    """Solubility product of aragonite (Mucci 1983)."""
    log_Ksp = (-171.945 - 0.077993 * T_K + 2903.293 / T_K
               + 71.595 * np.log10(T_K)
               + (-0.068393 + 0.0017276 * T_K + 88.135 / T_K) * S**0.5
               - 0.10018 * S + 0.0059415 * S**1.5)
    return 10**log_Ksp

def compute_Kw(T_K, S):
    """Ion product of water in seawater."""
    lnKw = (148.9652 - 13847.26 / T_K - 23.6521 * np.log(T_K)
            + (-5.977 + 118.67 / T_K + 1.0495 * np.log(T_K)) * S**0.5
            - 0.01615 * S)
    return np.exp(lnKw)

def compute_KB(T_K, S):
    """Dissociation constant of boric acid (Dickson 1990)."""
    lnKB = ((-8966.90 - 2890.53 * S**0.5 - 77.942 * S
             + 1.728 * S**1.5 - 0.0996 * S**2) / T_K
            + 148.0248 + 137.1942 * S**0.5 + 1.62142 * S
            - (24.4344 + 25.085 * S**0.5 + 0.2474 * S) * np.log(T_K)
            + 0.053105 * S**0.5 * T_K)
    return np.exp(lnKB)

def carbonate_system(T_C, S, DIC, TA, pCO2_atm=None):
    """
    Solve the marine carbonate system given T, S, DIC, and TA.
    
    Parameters
    ----------
    T_C : float - Temperature in Celsius
    S : float - Salinity (PSU)
    DIC : float - Dissolved Inorganic Carbon (µmol/kg)
    TA : float - Total Alkalinity (µmol/kg)
    pCO2_atm : float, optional - Atmospheric pCO2 (µatm)
    
    Returns
    -------
    dict with pH, pCO2, [CO2], [HCO3-], [CO3 2-], Omega_aragonite
    """
    T_K = T_C + 273.15
    DIC_mol = DIC * 1e-6  # convert to mol/kg
    TA_mol = TA * 1e-6

    K1 = compute_K1(T_K, S)
    K2 = compute_K2(T_K, S)
    Kw = compute_Kw(T_K, S)
    KB = compute_KB(T_K, S)
    Ksp_arag = compute_Ksp_aragonite(T_K, S)

    # Total boron (Uppström 1974)
    BT = 0.000416 * S / 35.0  # mol/kg

    # Ca2+ concentration (Riley & Tongudai 1967)
    Ca = 0.01028 * S / 35.0  # mol/kg

    # Solve for [H+] from TA and DIC
    def alkalinity_eq(H):
        H = H[0]
        if H <= 0:
            return [1e10]
        CO2aq = DIC_mol / (1 + K1/H + K1*K2/H**2)
        HCO3 = DIC_mol * K1/H / (1 + K1/H + K1*K2/H**2)
        CO3 = DIC_mol * K1*K2/H**2 / (1 + K1/H + K1*K2/H**2)
        BOH4 = BT * KB / (KB + H)
        OH = Kw / H
        TA_calc = HCO3 + 2*CO3 + BOH4 + OH - H
        return [TA_calc - TA_mol]

    H_init = 10**(-8.1)
    H_solution = fsolve(alkalinity_eq, [H_init], full_output=False)[0]
    pH = -np.log10(H_solution)

    CO2aq = DIC_mol / (1 + K1/H_solution + K1*K2/H_solution**2)
    HCO3 = DIC_mol * K1/H_solution / (1 + K1/H_solution + K1*K2/H_solution**2)
    CO3 = DIC_mol * K1*K2/H_solution**2 / (1 + K1/H_solution + K1*K2/H_solution**2)

    # Henry's law constant (Weiss 1974)
    K0 = np.exp(93.4517 * (100/T_K) - 60.2409 + 23.3585 * np.log(T_K/100)
                + S * (0.023517 - 0.023656 * (T_K/100) + 0.0047036 * (T_K/100)**2))
    pCO2_calc = CO2aq / K0 * 1e6  # µatm

    Omega_arag = Ca * CO3 / Ksp_arag

    return {
        'pH': round(pH, 4),
        'pCO2_uatm': round(pCO2_calc, 1),
        'CO2aq_umol_kg': round(CO2aq * 1e6, 2),
        'HCO3_umol_kg': round(HCO3 * 1e6, 2),
        'CO3_umol_kg': round(CO3 * 1e6, 2),
        'Omega_aragonite': round(Omega_arag, 4),
        'T_C': T_C,
        'S': S,
        'DIC_umol_kg': DIC,
        'TA_umol_kg': TA
    }


def project_carbonate_scenarios():
    """
    Project carbonate chemistry under RCP scenarios (2020-2100).
    Uses atmospheric pCO2 trajectories to compute DIC changes.
    """
    # RCP atmospheric pCO2 trajectories (µatm)
    years = np.arange(2020, 2101, 5)
    scenarios = {
        'RCP2.6': 410 + 10 * np.sin(np.linspace(0, np.pi, len(years))) * 0.8,
        'RCP4.5': 410 + np.linspace(0, 140, len(years)),
        'RCP8.5': 410 + np.linspace(0, 550, len(years))
    }
    # Adjust RCP2.6 to peak and decline
    scenarios['RCP2.6'] = np.concatenate([
        np.linspace(410, 440, len(years)//2),
        np.linspace(440, 420, len(years) - len(years)//2)
    ])

    T_base = 25.0  # GBR baseline SST
    S = 35.0
    TA = 2300.0  # µmol/kg (approximately constant)

    # Temperature projections
    temp_increase = {
        'RCP2.6': np.linspace(0, 0.8, len(years)),
        'RCP4.5': np.linspace(0, 1.8, len(years)),
        'RCP8.5': np.linspace(0, 3.7, len(years))
    }

    results = {}
    for scenario in scenarios:
        scenario_results = []
        for i, year in enumerate(years):
            pCO2 = scenarios[scenario][i]
            T = T_base + temp_increase[scenario][i]
            T_K = T + 273.15
            K0 = np.exp(93.4517 * (100/T_K) - 60.2409 + 23.3585 * np.log(T_K/100)
                        + S * (0.023517 - 0.023656 * (T_K/100) + 0.0047036 * (T_K/100)**2))
            # Estimate DIC from pCO2 and TA
            CO2aq = K0 * pCO2 * 1e-6
            K1 = compute_K1(T_K, S)
            K2 = compute_K2(T_K, S)
            # Iteratively estimate DIC
            H = 10**(-8.1)
            for _ in range(50):
                alpha0 = 1 / (1 + K1/H + K1*K2/H**2)
                DIC_est = CO2aq / alpha0
                HCO3 = DIC_est * K1/H / (1 + K1/H + K1*K2/H**2)
                CO3 = DIC_est * K1*K2/H**2 / (1 + K1/H + K1*K2/H**2)
                KB = compute_KB(T_K, S)
                BT = 0.000416 * S / 35.0
                BOH4 = BT * KB / (KB + H)
                Kw = compute_Kw(T_K, S)
                OH = Kw / H
                TA_calc = HCO3 + 2*CO3 + BOH4 + OH - H
                TA_diff = TA * 1e-6 - TA_calc
                H = H * np.exp(-TA_diff / (HCO3 + 4*CO3))

            DIC = DIC_est * 1e6
            result = carbonate_system(T, S, DIC, TA)
            result['year'] = int(year)
            result['scenario'] = scenario
            result['T_C'] = round(T, 2)
            scenario_results.append(result)
        results[scenario] = scenario_results

    return results, years


if __name__ == '__main__':
    # Present-day GBR conditions
    print("=== Present-day GBR Carbonate Chemistry ===")
    present = carbonate_system(T_C=25.0, S=35.0, DIC=2050, TA=2300)
    for k, v in present.items():
        print(f"  {k}: {v}")

    # Run projections
    print("\n=== Projected Carbonate Chemistry ===")
    results, years = project_carbonate_scenarios()

    # Save results
    import json
    with open('results/carbonate_chemistry.json', 'w') as f:
        json.dump(results, f, indent=2)
    print("Results saved to results/carbonate_chemistry.json")
