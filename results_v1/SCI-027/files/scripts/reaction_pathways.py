"""
CO2RR Reaction Pathway Analysis
Implements the Computational Hydrogen Electrode (CHE) model (Norskov et al., 2004)
to calculate free energy diagrams for CO2 reduction pathways.

Pathways:
  Path A: CO2 → COOH* → CO* → (desorb as CO)
  Path B: CO2 → COOH* → CO* → CHO* → CH2O* → OCH3* → CH4
  Path C: CO2 → HCOO* → HCOOH (formic acid pathway)
  Path D: CO* + CO* → OC-CO* → ... → C2H4 (C-C coupling)
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import json
from datetime import datetime

# -----------------------------------------------------------------
# CHE constants
# -----------------------------------------------------------------
kB    = 8.617e-5   # eV/K
T     = 298.15     # K
kBT   = kB * T

# Reference chemical potentials (computational hydrogen electrode)
# At U=0 V vs. RHE: μ(H+ + e-) = 0.5 * μ(H2)
# G(H2O) - G(H2) - ... handled by CHE

# Zero-point energy and entropy corrections for CO2RR intermediates (eV)
ZPE_CORR = {
    "CO2(g)":  0.31,
    "COOH*":   0.52,
    "CO*":     0.18,
    "CHO*":    0.42,
    "COH*":    0.38,
    "CH2O*":   0.60,
    "OCH3*":   0.78,
    "CH3*":    1.05,
    "CH4(g)":  1.20,
    "HCOO*":   0.55,
    "HCOOH(g)":0.68,
    "OCCO*":   0.45,
    "OCCHO*":  0.62,
    "C2H4(g)": 1.52,
    "H2O(l)":  0.56,
    "CO(g)":   0.13,
}

TS_CORR = {   # T*S at 298 K (eV)  — vibrational contribution
    "CO2(g)":  0.66,
    "COOH*":   0.12,
    "CO*":     0.06,
    "CHO*":    0.10,
    "COH*":    0.09,
    "CH2O*":   0.12,
    "OCH3*":   0.18,
    "CH3*":    0.22,
    "CH4(g)":  0.62,
    "HCOO*":   0.13,
    "HCOOH(g)":0.55,
    "OCCO*":   0.11,
    "OCCHO*":  0.15,
    "C2H4(g)": 0.78,
    "H2O(l)":  0.58,
    "CO(g)":   0.61,
}


def zpe_ts(species: str) -> float:
    """Return ZPE - TS correction (eV)."""
    return ZPE_CORR.get(species, 0.0) - TS_CORR.get(species, 0.0)


def che_free_energy(dG_ads: float, U: float, n_electrons: int,
                    species: str = "") -> float:
    """
    Apply CHE: ΔG(U) = ΔG(0) - n*e*U + ZPE - TS
    Args:
        dG_ads:     DFT adsorption energy (eV)
        U:          electrode potential vs. RHE (V)
        n_electrons:number of electrons transferred to reach this state
        species:    species name for ZPE/TS lookup
    """
    return dG_ads - n_electrons * U + zpe_ts(species)


# -----------------------------------------------------------------
# Reference free energies (computational hydrogen electrode, eV)
# CO2(g) → reference = 0
# All free energies relative to CO2(g) + H2O(l) reference
# -----------------------------------------------------------------
G_REF = {
    "CO2(g)":   0.000,
    "CO(g)":   -0.517,   # ΔG_rxn = -0.52 eV at standard conditions
    "HCOOH(g)":-0.250,
    "CH4(g)":  -1.060,
    "C2H4(g)": -1.320,
}


def compute_pathway_A(dG_COOH: float, dG_CO: float, U: float,
                       solvation_corr: float = 0.0) -> dict:
    """
    CO2 → CO pathway (2e- reduction, CO2RR → CO product)
    Steps:
      CO2(g) + H+ + e- → COOH*          (step 1)
      COOH* + H+ + e-  → CO* + H2O      (step 2)
      CO*              → CO(g)           (desorption, no electron)
    """
    # Free energy at each state (relative to CO2(g) = 0)
    G0 = 0.0                                              # CO2(g)
    G1 = dG_COOH + zpe_ts("COOH*") - U + solvation_corr  # COOH* (1e-)
    G2 = dG_CO   + zpe_ts("CO*")   - 2*U + 0.27          # CO* (2e-), +0.27 from H2O
    G3 = G_REF["CO(g)"] - 2*U                             # CO(g) desorbed

    dG_step1 = G1 - G0
    dG_step2 = G2 - G1
    dG_desorb = G3 - G2

    return {
        "pathway": "CO2→CO",
        "states":  ["CO2(g)", "COOH*", "CO*", "CO(g)"],
        "G":       [G0, G1, G2, G3],
        "dG_steps":[dG_step1, dG_step2, dG_desorb],
        "PDS_index": int(np.argmax([dG_step1, dG_step2, max(dG_desorb, 0)])),
        "limiting_potential": -max(dG_step1, dG_step2, max(dG_desorb, 0)),
    }


def compute_pathway_B(dG_CO: float, dG_CHO: float, dG_CH2O: float,
                      dG_OCH3: float, U: float,
                      solvation_corr: float = 0.0) -> dict:
    """
    CO → CH4 pathway (8e- total from CO2; here 6e- from CO*)
    CO* → CHO* → CH2O* → OCH3* → CH3* → CH4(g)
    """
    base = dG_CO + zpe_ts("CO*") - 2*U
    G_CO   = base
    G_CHO  = G_CO  + (dG_CHO  - dG_CO)  + zpe_ts("CHO*")  - zpe_ts("CO*")  - U
    G_CH2O = G_CHO + (dG_CH2O - dG_CHO) + zpe_ts("CH2O*") - zpe_ts("CHO*") - U
    G_OCH3 = G_CH2O+ (dG_OCH3 - dG_CH2O)+ zpe_ts("OCH3*") - zpe_ts("CH2O*")- U
    G_CH4  = G_REF["CH4(g)"] - 8*U   # 8 electrons total from CO2

    dG_steps = [
        G_CHO - G_CO,
        G_CH2O - G_CHO,
        G_OCH3 - G_CH2O,
        G_CH4  - G_OCH3,
    ]

    return {
        "pathway":  "CO*→CH4",
        "states":   ["CO*", "CHO*", "CH2O*", "OCH3*", "CH4(g)"],
        "G":        [G_CO, G_CHO, G_CH2O, G_OCH3, G_CH4],
        "dG_steps": dG_steps,
        "PDS_index": int(np.argmax(dG_steps)),
        "limiting_potential": -max(dG_steps),
    }


def compute_pathway_D(dG_CO: float, dG_OCCO: float, dG_OCCHO: float,
                      U: float) -> dict:
    """
    C-C coupling pathway: CO* + CO* → OC-CO* → OC-CHO* → … → C2H4
    (12e- total from 2×CO2)
    """
    G_CO    = dG_CO   + zpe_ts("CO*")   - 2*U
    G_OCCO  = dG_OCCO + zpe_ts("OCCO*") - 4*U     # coupling: 2×CO*
    G_OCCHO = dG_OCCHO+ zpe_ts("OCCHO*")- 5*U
    G_C2H4  = G_REF["C2H4(g)"] - 12*U

    dG_steps = [
        G_OCCO  - 2*G_CO,    # C-C coupling
        G_OCCHO - G_OCCO,    # first hydrogenation of OCCO
        G_C2H4  - G_OCCHO,   # remaining steps lumped
    ]

    return {
        "pathway":  "2CO*→C2H4",
        "states":   ["2×CO*", "OC-CO*", "OC-CHO*", "C2H4(g)"],
        "G":        [2*G_CO, G_OCCO, G_OCCHO, G_C2H4],
        "dG_steps": dG_steps,
        "PDS_index": int(np.argmax(dG_steps)),
        "limiting_potential": -max(dG_steps),
    }


def analyze_catalyst(catalyst_name: str, data: dict, U: float = 0.0,
                     solvation: float = 0.0) -> dict:
    """Run all pathways for a given catalyst."""
    results = {"catalyst": catalyst_name, "U": U}

    # Pathway A: CO2 → CO
    pA = compute_pathway_A(data["dG_COOH"], data["dG_CO"], U, solvation)
    results["CO_limiting_U"]  = pA["limiting_potential"]
    results["CO_PDS"]         = pA["PDS_index"]

    # Pathway B: CO → CH4
    pB = compute_pathway_B(data["dG_CO"], data["dG_CHO"],
                            data["dG_CH2O"], data["dG_OCH3"], U, solvation)
    results["CH4_limiting_U"] = pB["limiting_potential"]
    results["CH4_PDS"]        = pB["PDS_index"]

    # Pathway D: C-C coupling
    if "dG_OCCO" in data and data["dG_OCCO"] is not None:
        pD = compute_pathway_D(data["dG_CO"], data["dG_OCCO"],
                                data.get("dG_OCCHO", data["dG_CO"]+0.3), U)
        results["C2H4_limiting_U"] = pD["limiting_potential"]
        results["C2H4_PDS"]        = pD["PDS_index"]

    results["_pathways"] = {"A": pA, "B": pB}
    return results


def plot_free_energy_diagrams(cu_data: dict, save_path: str,
                              U_values: list = [0.0, -0.5, -0.8]) -> None:
    """
    Plot free energy diagrams for Cu at multiple potentials.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    colors = ["#2196F3", "#FF9800", "#F44336"]

    for ax, U, col in zip(axes, U_values, colors):
        pA = compute_pathway_A(cu_data["dG_COOH"], cu_data["dG_CO"], U)
        ax.set_title(f"U = {U:.1f} V vs. RHE", fontsize=11, fontweight="bold")
        Gs = pA["G"]
        steps = list(range(len(Gs)))
        # Draw horizontal lines for each state
        for i, (step, g) in enumerate(zip(steps, Gs)):
            ax.hlines(g, step - 0.3, step + 0.3, colors=col, linewidth=2.5)
            if i < len(Gs) - 1:
                ax.plot([step + 0.3, step + 0.7], [g, Gs[i+1]],
                        color=col, linestyle="--", linewidth=1, alpha=0.6)
            ax.text(step, g + 0.04, f"{g:.2f}", ha="center",
                    fontsize=8, color=col)

        ax.set_xticks(steps)
        ax.set_xticklabels(pA["states"], fontsize=9, rotation=15)
        ax.set_ylabel("Free Energy (eV)" if ax == axes[0] else "")
        ax.axhline(0, color="gray", linewidth=0.5, linestyle=":")
        ax.set_xlabel("Reaction Coordinate")
        lp = pA["limiting_potential"]
        ax.annotate(f"U_L = {lp:.2f} V",
                    xy=(0.05, 0.92), xycoords="axes fraction",
                    fontsize=9, color=col,
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))

    fig.suptitle("CO2 → CO Free Energy Diagrams on Cu(111)", fontsize=13,
                 fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {save_path}")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ".")
    from data.adsorption_energies import PURE_METALS

    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)

    cu = PURE_METALS["Cu"]
    plot_free_energy_diagrams(cu, "figures/free_energy_diagram_Cu.png")

    results = []
    for name, data in PURE_METALS.items():
        r = analyze_catalyst(name, data)
        results.append({k: v for k, v in r.items() if k != "_pathways"})

    df = pd.DataFrame(results)
    df.to_csv("results/pathway_analysis.csv", index=False)
    print(df[["catalyst","CO_limiting_U","CH4_limiting_U"]].to_string())
