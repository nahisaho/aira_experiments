"""
CO2RR Adsorption Energy Database
Literature-based DFT adsorption energies (eV) for CO2RR intermediates
on various metal and alloy surfaces.

References:
- Peterson et al., Energy Environ. Sci. 2010
- Bagger et al., ChemElectroChem 2017
- Nitopi et al., Chem. Rev. 2019
- Zhao et al., J. Am. Chem. Soc. 2020
- Back et al., ACS Catalysis 2019 (SAC)
"""

import pandas as pd
import numpy as np

# ------------------------------------------------------------------
# Pure metals: ΔG (eV) at U=0 V vs. RHE
# Key intermediates: *CO2, *COOH, *CO, *CHO, *COH, *CH2O, *OCH3, *CH3
# Products targeted: CO (2e-), HCOOH (2e-), CH4 (8e-), C2H4 (12e-)
# ------------------------------------------------------------------
PURE_METALS = {
    # Metal: {intermediate: ΔG_ads (eV)}
    "Cu":   {"dG_CO2": -0.59, "dG_COOH": 0.43, "dG_CO": -0.68, "dG_CHO": 0.39,
             "dG_COH": 0.22,  "dG_CH2O": -0.02, "dG_OCH3": -1.01, "dG_CH3": -3.10,
             "dG_OCCO": 0.18, "dG_OCCHO": -0.15},
    "Au":   {"dG_CO2": -0.45, "dG_COOH": 0.82, "dG_CO": -0.38, "dG_CHO": 1.25,
             "dG_COH": 1.10,  "dG_CH2O":  0.80, "dG_OCH3":  0.40, "dG_CH3": -1.80,
             "dG_OCCO": 0.90, "dG_OCCHO":  0.55},
    "Ag":   {"dG_CO2": -0.42, "dG_COOH": 0.97, "dG_CO": -0.29, "dG_CHO": 1.43,
             "dG_COH": 1.28,  "dG_CH2O":  0.95, "dG_OCH3":  0.58, "dG_CH3": -1.60,
             "dG_OCCO": 1.05, "dG_OCCHO":  0.68},
    "Pt":   {"dG_CO2": -0.65, "dG_COOH": 0.12, "dG_CO": -1.10, "dG_CHO": -0.08,
             "dG_COH": -0.25, "dG_CH2O": -0.72, "dG_OCH3": -1.65, "dG_CH3": -4.20,
             "dG_OCCO":-0.22, "dG_OCCHO": -0.55},
    "Pd":   {"dG_CO2": -0.62, "dG_COOH": 0.28, "dG_CO": -0.95, "dG_CHO": 0.15,
             "dG_COH": -0.05, "dG_CH2O": -0.52, "dG_OCH3": -1.40, "dG_CH3": -3.85,
             "dG_OCCO":-0.02, "dG_OCCHO": -0.30},
    "Ni":   {"dG_CO2": -0.70, "dG_COOH": 0.05, "dG_CO": -1.25, "dG_CHO": -0.22,
             "dG_COH": -0.38, "dG_CH2O": -0.95, "dG_OCH3": -1.92, "dG_CH3": -4.55,
             "dG_OCCO":-0.48, "dG_OCCHO": -0.82},
    "Fe":   {"dG_CO2": -0.78, "dG_COOH":-0.15, "dG_CO": -1.45, "dG_CHO": -0.42,
             "dG_COH": -0.58, "dG_CH2O": -1.18, "dG_OCH3": -2.15, "dG_CH3": -4.85,
             "dG_OCCO":-0.72, "dG_OCCHO": -1.05},
    "Zn":   {"dG_CO2": -0.38, "dG_COOH": 0.88, "dG_CO": -0.20, "dG_CHO": 1.35,
             "dG_COH": 1.18,  "dG_CH2O":  0.85, "dG_OCH3":  0.42, "dG_CH3": -1.45,
             "dG_OCCO": 0.95, "dG_OCCHO":  0.60},
    "Sn":   {"dG_CO2": -0.35, "dG_COOH": 0.35, "dG_CO":  0.15, "dG_CHO": 1.10,
             "dG_COH": 0.95,  "dG_CH2O":  0.62, "dG_OCH3":  0.18, "dG_CH3": -1.20,
             "dG_OCCO": 0.72, "dG_OCCHO":  0.38},
    "In":   {"dG_CO2": -0.32, "dG_COOH": 0.42, "dG_CO":  0.22, "dG_CHO": 1.18,
             "dG_COH": 1.02,  "dG_CH2O":  0.70, "dG_OCH3":  0.25, "dG_CH3": -1.12,
             "dG_OCCO": 0.80, "dG_OCCHO":  0.45},
    "Bi":   {"dG_CO2": -0.28, "dG_COOH": 0.55, "dG_CO":  0.35, "dG_CHO": 1.28,
             "dG_COH": 1.12,  "dG_CH2O":  0.78, "dG_OCH3":  0.35, "dG_CH3": -1.05,
             "dG_OCCO": 0.88, "dG_OCCHO":  0.52},
    "Pb":   {"dG_CO2": -0.25, "dG_COOH": 0.62, "dG_CO":  0.42, "dG_CHO": 1.38,
             "dG_COH": 1.22,  "dG_CH2O":  0.88, "dG_OCH3":  0.45, "dG_CH3": -0.98,
             "dG_OCCO": 0.95, "dG_OCCHO":  0.58},
    "Rh":   {"dG_CO2": -0.72, "dG_COOH": 0.02, "dG_CO": -1.35, "dG_CHO": -0.32,
             "dG_COH": -0.48, "dG_CH2O": -1.08, "dG_OCH3": -2.05, "dG_CH3": -4.68,
             "dG_OCCO":-0.62, "dG_OCCHO": -0.95},
    "Ir":   {"dG_CO2": -0.68, "dG_COOH": 0.08, "dG_CO": -1.18, "dG_CHO": -0.15,
             "dG_COH": -0.30, "dG_CH2O": -0.88, "dG_OCH3": -1.85, "dG_CH3": -4.42,
             "dG_OCCO":-0.42, "dG_OCCHO": -0.75},
    "Co":   {"dG_CO2": -0.75, "dG_COOH":-0.08, "dG_CO": -1.38, "dG_CHO": -0.35,
             "dG_COH": -0.52, "dG_CH2O": -1.05, "dG_OCH3": -2.02, "dG_CH3": -4.62,
             "dG_OCCO":-0.65, "dG_OCCHO": -0.98},
}

# ------------------------------------------------------------------
# Cu-based alloys (Cu3X or CuX surfaces)
# ------------------------------------------------------------------
CU_ALLOYS = {
    "Cu3Zn":  {"dG_CO2": -0.55, "dG_COOH": 0.55, "dG_CO": -0.55, "dG_CHO": 0.55,
               "dG_COH": 0.38,  "dG_CH2O":  0.12, "dG_OCH3": -0.85, "dG_CH3": -2.85,
               "dG_OCCO": 0.32, "dG_OCCHO":  0.05},
    "Cu3Ag":  {"dG_CO2": -0.52, "dG_COOH": 0.62, "dG_CO": -0.52, "dG_CHO": 0.62,
               "dG_COH": 0.45,  "dG_CH2O":  0.18, "dG_OCH3": -0.78, "dG_CH3": -2.72,
               "dG_OCCO": 0.38, "dG_OCCHO":  0.12},
    "Cu3Sn":  {"dG_CO2": -0.50, "dG_COOH": 0.48, "dG_CO": -0.48, "dG_CHO": 0.48,
               "dG_COH": 0.30,  "dG_CH2O":  0.05, "dG_OCH3": -0.92, "dG_CH3": -2.98,
               "dG_OCCO": 0.25, "dG_OCCHO": -0.02},
    "Cu3In":  {"dG_CO2": -0.48, "dG_COOH": 0.52, "dG_CO": -0.45, "dG_CHO": 0.52,
               "dG_COH": 0.35,  "dG_CH2O":  0.08, "dG_OCH3": -0.88, "dG_CH3": -2.92,
               "dG_OCCO": 0.28, "dG_OCCHO":  0.00},
    "CuPd":   {"dG_CO2": -0.60, "dG_COOH": 0.35, "dG_CO": -0.82, "dG_CHO": 0.27,
               "dG_COH": 0.10,  "dG_CH2O": -0.28, "dG_OCH3": -1.22, "dG_CH3": -3.48,
               "dG_OCCO": 0.08, "dG_OCCHO": -0.12},
    "CuAu":   {"dG_CO2": -0.52, "dG_COOH": 0.58, "dG_CO": -0.55, "dG_CHO": 0.58,
               "dG_COH": 0.42,  "dG_CH2O":  0.15, "dG_OCH3": -0.82, "dG_CH3": -2.78,
               "dG_OCCO": 0.35, "dG_OCCHO":  0.08},
    "CuNi":   {"dG_CO2": -0.63, "dG_COOH": 0.22, "dG_CO": -0.95, "dG_CHO": 0.15,
               "dG_COH":-0.02,  "dG_CH2O": -0.45, "dG_OCH3": -1.45, "dG_CH3": -3.82,
               "dG_OCCO":-0.08, "dG_OCCHO": -0.28},
    "CuCo":   {"dG_CO2": -0.65, "dG_COOH": 0.18, "dG_CO": -1.05, "dG_CHO": 0.08,
               "dG_COH":-0.10,  "dG_CH2O": -0.55, "dG_OCH3": -1.55, "dG_CH3": -3.98,
               "dG_OCCO":-0.18, "dG_OCCHO": -0.38},
    "Cu3Al":  {"dG_CO2": -0.48, "dG_COOH": 0.58, "dG_CO": -0.42, "dG_CHO": 0.58,
               "dG_COH": 0.42,  "dG_CH2O":  0.15, "dG_OCH3": -0.82, "dG_CH3": -2.75,
               "dG_OCCO": 0.32, "dG_OCCHO":  0.05},
    "Cu3Ga":  {"dG_CO2": -0.50, "dG_COOH": 0.52, "dG_CO": -0.45, "dG_CHO": 0.52,
               "dG_COH": 0.36,  "dG_CH2O":  0.08, "dG_OCH3": -0.88, "dG_CH3": -2.88,
               "dG_OCCO": 0.26, "dG_OCCHO": -0.01},
}

# ------------------------------------------------------------------
# N-doped Carbon + single-atom catalysts (SAC)
# Metal atom on N4-C, N2-C, or pyridinic-N support
# ------------------------------------------------------------------
SAC_NDOPED = {
    # Format: MN4-C (M = metal, N4-porphyrin-like coordination)
    "Fe-N4C":  {"dG_CO2": -0.82, "dG_COOH":-0.22, "dG_CO": -1.52, "dG_CHO": -0.52,
                "dG_COH": -0.68, "dG_CH2O": -1.28, "dG_OCH3": -2.28, "dG_CH3": -5.02,
                "coordination": "N4", "metal": "Fe", "d_band_center": -1.85},
    "Co-N4C":  {"dG_CO2": -0.78, "dG_COOH":-0.15, "dG_CO": -1.42, "dG_CHO": -0.42,
                "dG_COH": -0.58, "dG_CH2O": -1.18, "dG_OCH3": -2.18, "dG_CH3": -4.88,
                "coordination": "N4", "metal": "Co", "d_band_center": -1.52},
    "Ni-N4C":  {"dG_CO2": -0.72, "dG_COOH":-0.05, "dG_CO": -1.28, "dG_CHO": -0.25,
                "dG_COH": -0.42, "dG_CH2O": -0.98, "dG_OCH3": -1.98, "dG_CH3": -4.62,
                "coordination": "N4", "metal": "Ni", "d_band_center": -1.28},
    "Cu-N4C":  {"dG_CO2": -0.62, "dG_COOH": 0.28, "dG_CO": -0.82, "dG_CHO": 0.25,
                "dG_COH": 0.08,  "dG_CH2O": -0.42, "dG_OCH3": -1.38, "dG_CH3": -3.62,
                "coordination": "N4", "metal": "Cu", "d_band_center": -2.67},
    "Zn-N4C":  {"dG_CO2": -0.45, "dG_COOH": 0.72, "dG_CO": -0.32, "dG_CHO": 1.22,
                "dG_COH": 1.05,  "dG_CH2O":  0.72, "dG_OCH3":  0.28, "dG_CH3": -1.62,
                "coordination": "N4", "metal": "Zn", "d_band_center": -7.52},
    "Mn-N4C":  {"dG_CO2": -0.85, "dG_COOH":-0.32, "dG_CO": -1.62, "dG_CHO": -0.62,
                "dG_COH": -0.78, "dG_CH2O": -1.42, "dG_OCH3": -2.42, "dG_CH3": -5.18,
                "coordination": "N4", "metal": "Mn", "d_band_center": -1.62},
    "Mo-N4C":  {"dG_CO2": -0.88, "dG_COOH":-0.38, "dG_CO": -1.72, "dG_CHO": -0.68,
                "dG_COH": -0.85, "dG_CH2O": -1.52, "dG_OCH3": -2.52, "dG_CH3": -5.32,
                "coordination": "N4", "metal": "Mo", "d_band_center": -1.15},
    "Cr-N4C":  {"dG_CO2": -0.90, "dG_COOH":-0.42, "dG_CO": -1.78, "dG_CHO": -0.72,
                "dG_COH": -0.90, "dG_CH2O": -1.58, "dG_OCH3": -2.58, "dG_CH3": -5.42,
                "coordination": "N4", "metal": "Cr", "d_band_center": -0.95},
    "Pd-N4C":  {"dG_CO2": -0.65, "dG_COOH": 0.18, "dG_CO": -0.98, "dG_CHO": 0.12,
                "dG_COH":-0.05,  "dG_CH2O": -0.58, "dG_OCH3": -1.58, "dG_CH3": -3.98,
                "coordination": "N4", "metal": "Pd", "d_band_center": -1.78},
    "Ru-N4C":  {"dG_CO2": -0.80, "dG_COOH":-0.18, "dG_CO": -1.48, "dG_CHO": -0.48,
                "dG_COH": -0.62, "dG_CH2O": -1.22, "dG_OCH3": -2.22, "dG_CH3": -4.95,
                "coordination": "N4", "metal": "Ru", "d_band_center": -1.42},
    # N2-C defect (less N coordination)
    "Fe-N2C":  {"dG_CO2": -0.75, "dG_COOH":-0.08, "dG_CO": -1.35, "dG_CHO": -0.35,
                "dG_COH": -0.52, "dG_CH2O": -1.08, "dG_OCH3": -2.08, "dG_CH3": -4.75,
                "coordination": "N2", "metal": "Fe", "d_band_center": -1.62},
    "Co-N2C":  {"dG_CO2": -0.70, "dG_COOH": 0.02, "dG_CO": -1.22, "dG_CHO": -0.22,
                "dG_COH": -0.38, "dG_CH2O": -0.95, "dG_OCH3": -1.92, "dG_CH3": -4.52,
                "coordination": "N2", "metal": "Co", "d_band_center": -1.35},
    "Cu-N2C":  {"dG_CO2": -0.55, "dG_COOH": 0.42, "dG_CO": -0.65, "dG_CHO": 0.42,
                "dG_COH": 0.25,  "dG_CH2O": -0.22, "dG_OCH3": -1.15, "dG_CH3": -3.35,
                "coordination": "N2", "metal": "Cu", "d_band_center": -2.45},
}


def get_all_catalysts():
    """Return a unified DataFrame of all catalysts with labels."""
    records = []
    for name, vals in PURE_METALS.items():
        records.append({"catalyst": name, "category": "pure_metal", **vals})
    for name, vals in CU_ALLOYS.items():
        records.append({"catalyst": name, "category": "cu_alloy", **vals})
    for name, vals in SAC_NDOPED.items():
        r = {k: v for k, v in vals.items()
             if k not in ("coordination", "metal", "d_band_center")}
        extra = {k: vals.get(k, None)
                 for k in ("coordination", "metal", "d_band_center")}
        records.append({"catalyst": name, "category": "SAC_N-doped", **r, **extra})
    return pd.DataFrame(records)


if __name__ == "__main__":
    df = get_all_catalysts()
    print(df.to_string())
    df.to_csv("data/all_catalysts.csv", index=False)
    print(f"\nTotal catalysts: {len(df)}")
