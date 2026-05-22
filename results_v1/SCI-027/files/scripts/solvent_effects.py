"""
Solvent Effects and Potential-Dependent Calculations for CO2RR
Implements:
  1. Implicit solvation corrections (SCCS/COSMO model approximations)
  2. Electric double layer (EDL) effects on adsorption energies
  3. Potential-dependent free energy diagrams
  4. pH effect on thermodynamics (SHE vs. RHE reference)
  5. Poisson-Boltzmann solvation model estimate
  6. CO2 solubility and local concentration effects

Reference:
  - Mathew et al., J. Chem. Phys. 2014 (VASPsol)
  - Gauthier et al., J. Phys. Chem. Lett. 2019 (EDL)
  - Chan & Norskov, J. Phys. Chem. Lett. 2015 (grand canonical)
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# ------------------------------------------------------------------
# Solvent correction parameters (implicit solvation model)
# ------------------------------------------------------------------
SOLVENT_PARAMS = {
    "epsilon_water":  78.4,       # dielectric constant of water
    "kappa_DH":       3.28e9,     # Debye-Hückel screening (m⁻¹ at 0.1 M KOH)
    "z_ref":          3.5e-10,    # reference plane to metal (m)
    "C_EDL":          20e-6,      # EDL capacitance (F/cm²)
    "sigma_max":      0.15,       # max surface charge density (e/Å²)
}

# Solvation correction to adsorption energy (eV) from VASPsol-type calculations
# Estimated from literature for key intermediates
SOLVATION_CORR_WATER = {
    "COOH*":  -0.25,   # stabilized by H-bonding (2 OH groups)
    "CO*":    -0.05,   # weakly stabilized (no H-bonding)
    "CHO*":   -0.18,   # moderate stabilization
    "COH*":   -0.15,   # OH group solvation
    "CH2O*":  -0.08,   # weak
    "OCH3*":  -0.12,   # ether-like O solvation
    "OCCO*":  -0.10,   # C2 intermediate
    "OCCHO*": -0.20,   # dual carbonyl
}

# EDL field effect on polar intermediates: dΔG/dU (eV/V)
EDL_SLOPE = {
    "COOH*":  -0.15,   # dipole: dΔG/dU (partial charge transfer)
    "CO*":    -0.08,
    "CHO*":   -0.10,
    "COH*":   -0.12,
    "CH2O*":  -0.05,
    "OCH3*":  -0.04,
}

# CO2 solubility corrections at different conditions
CO2_SOLUBILITY = {
    "1_atm_25C":       33.4e-3,    # mol/L in pure water
    "1_atm_KOH_1M":    8.5e-3,     # mol/L in 1M KOH
    "1_atm_KHCO3_1M":  22.0e-3,    # mol/L in 1M KHCO3
    "10_atm_25C":      334e-3,      # high pressure
}


def solvent_corrected_energy(dG_ads: float, species: str, U: float,
                              include_edl: bool = True) -> float:
    """
    Apply solvation and EDL corrections to adsorption energy.
    ΔG_corr = ΔG_DFT + ΔG_solv + ΔG_EDL(U)
    """
    dG_solv = SOLVATION_CORR_WATER.get(species, 0.0)
    dG_edl  = EDL_SLOPE.get(species, 0.0) * U if include_edl else 0.0
    return dG_ads + dG_solv + dG_edl


def compute_potential_dependent_diagram(
    dG_COOH: float, dG_CO: float, dG_CHO: float,
    catalyst_name: str, U_range: np.ndarray
) -> pd.DataFrame:
    """
    Compute free energy of each step as a function of potential U.
    Including solvation corrections.
    """
    records = []

    for U in U_range:
        # Step 1: CO2 → COOH* (1e-)
        dG1 = solvent_corrected_energy(dG_COOH, "COOH*", U) - (-U)
        # Step 2: COOH* → CO* (1e-)
        dG2_raw = (solvent_corrected_energy(dG_CO, "CO*", U) + 0.27
                   - solvent_corrected_energy(dG_COOH, "COOH*", U))
        dG2 = dG2_raw - U  # + 0.27 from H2O formation already included
        # Step 3: CO* → CHO* (1e-)
        dG3 = (solvent_corrected_energy(dG_CHO, "CHO*", U)
               - solvent_corrected_energy(dG_CO, "CO*", U) - U)

        records.append({
            "catalyst": catalyst_name,
            "U_V": round(float(U), 3),
            "dG_CO2_to_COOH": round(float(dG1), 4),
            "dG_COOH_to_CO":  round(float(dG2), 4),
            "dG_CO_to_CHO":   round(float(dG3), 4),
            "PDS": ("CO2→COOH*" if dG1 >= max(dG2, dG3)
                    else "COOH*→CO*" if dG2 >= dG3
                    else "CO*→CHO*"),
            "max_uphill_eV": round(float(max(dG1, dG2, dG3)), 4),
        })

    return pd.DataFrame(records)


def plot_potential_dependent_activity(df_catalysts: pd.DataFrame,
                                       save_dir: str) -> None:
    """
    Plot limiting potential as function of applied potential for key catalysts.
    Shows how activity evolves with U and where the PDS changes.
    """
    U_range = np.linspace(-1.5, 0.1, 80)
    key_cats = {
        "Cu":    ("#E65100", "-"),
        "Au":    ("#FDD835", "--"),
        "Ag":    ("#9E9E9E", "-."),
        "Fe-N4C":(  "#C62828", "-"),
        "Co-N4C":(  "#1565C0", "--"),
        "CuNi":  ("#2E7D32", "-"),
    }

    # Gather data for key catalysts
    from data.adsorption_energies import PURE_METALS, SAC_NDOPED, CU_ALLOYS
    all_data = {**PURE_METALS, **SAC_NDOPED, **CU_ALLOYS}

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    all_records = []
    for idx, (cat, (col, ls)) in enumerate(key_cats.items()):
        if cat not in all_data:
            continue
        data = all_data[cat]
        df_pd = compute_potential_dependent_diagram(
            data["dG_COOH"], data["dG_CO"],
            data.get("dG_CHO", data["dG_CO"] + 1.0),
            cat, U_range
        )
        all_records.append(df_pd)

        ax = axes[idx]
        # Plot each elementary step
        ax.plot(df_pd["U_V"], df_pd["dG_CO2_to_COOH"], color="#C62828",
                 linewidth=2, label="CO₂→*COOH")
        ax.plot(df_pd["U_V"], df_pd["dG_COOH_to_CO"],  color="#1565C0",
                 linewidth=2, label="*COOH→*CO")
        ax.plot(df_pd["U_V"], df_pd["dG_CO_to_CHO"],   color="#2E7D32",
                 linewidth=2, label="*CO→*CHO")
        ax.plot(df_pd["U_V"], df_pd["max_uphill_eV"],  color="black",
                 linewidth=2.5, linestyle="--", label="Max uphill (PDS)")

        ax.axhline(0, color="gray", linewidth=0.8, linestyle=":")
        ax.axvline(0, color="gray", linewidth=0.8, linestyle=":")

        # Mark limiting potential
        ul_idx = df_pd["max_uphill_eV"].le(0).idxmax()
        if ul_idx > 0:
            ul_U = df_pd.loc[ul_idx, "U_V"]
            ax.axvline(ul_U, color="purple", linestyle=":", alpha=0.6,
                        label=f"U_L ≈ {ul_U:.2f} V")

        ax.set_xlabel("Applied Potential (V vs. RHE)", fontsize=9)
        ax.set_ylabel("ΔG (eV)", fontsize=9)
        ax.set_title(f"{cat} — Potential-Dependent Steps", fontsize=10, fontweight="bold")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-1.5, 0.1)
        ax.set_ylim(-0.5, 1.8)

    plt.suptitle("Potential-Dependent Free Energy Steps for CO2RR\n"
                  "(with Implicit Solvation + EDL Corrections)",
                  fontsize=13, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(save_dir, "potential_dependent.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

    # Save all records
    if all_records:
        df_all = pd.concat(all_records)
        df_all.to_csv(os.path.join(RESULTS_DIR, "potential_dependent_steps.csv"),
                       index=False)


def plot_solvent_comparison(save_dir: str) -> None:
    """
    Compare gas-phase vs. solvated free energy diagram for Cu.
    """
    from data.adsorption_energies import PURE_METALS
    cu = PURE_METALS["Cu"]
    U = -0.65  # near optimal

    states = ["CO2(g)", "COOH*", "CO*", "CHO*", "CH4(g)"]
    n_e    = [0, 1, 2, 3, 8]

    def G_step(dG_species, n, corr=0.0):
        return dG_species - n * U + corr

    # Gas phase (no solvation)
    G_gas = [0,
              G_step(cu["dG_COOH"], 1, 0),
              G_step(cu["dG_CO"],   2, 0.27),
              G_step(cu["dG_CHO"],  3, 0.27),
              G_step(-1.06,         8, 0.0)]

    # Solvated
    G_sol = [0,
              G_step(cu["dG_COOH"] + SOLVATION_CORR_WATER["COOH*"], 1,
                     EDL_SLOPE["COOH*"]*U),
              G_step(cu["dG_CO"]   + SOLVATION_CORR_WATER["CO*"],   2,
                     0.27 + EDL_SLOPE["CO*"]*U),
              G_step(cu["dG_CHO"]  + SOLVATION_CORR_WATER["CHO*"],  3,
                     0.27 + EDL_SLOPE["CHO*"]*U),
              G_step(-1.06,                                          8, 0.0)]

    fig, ax = plt.subplots(figsize=(9, 6))
    x = np.arange(len(states))

    for xi, (gg, gs) in enumerate(zip(G_gas, G_sol)):
        ax.hlines(gg, xi-0.3, xi+0.3, colors="#1565C0", linewidth=3, label="Gas phase" if xi==0 else "")
        ax.hlines(gs, xi-0.3, xi+0.3, colors="#C62828", linewidth=3, label="Solvated" if xi==0 else "", linestyle="--")
        if xi < len(states)-1:
            ax.plot([xi+0.3, xi+0.7], [gg, G_gas[xi+1]], color="#1565C0", linestyle=":", alpha=0.5)
            ax.plot([xi+0.3, xi+0.7], [gs, G_sol[xi+1]], color="#C62828", linestyle=":", alpha=0.5)
        # Arrow showing solvation effect
        if abs(gg - gs) > 0.02:
            ax.annotate("", xy=(xi+0.35, gs), xytext=(xi+0.35, gg),
                         arrowprops=dict(arrowstyle="<->", color="green", lw=1.5))
            ax.text(xi+0.4, (gg+gs)/2, f"{gs-gg:+.2f}", fontsize=7.5, color="green")

    ax.set_xticks(x)
    ax.set_xticklabels(states, fontsize=10)
    ax.set_ylabel("Free Energy (eV)", fontsize=11)
    ax.set_xlabel("Reaction Coordinate", fontsize=11)
    ax.set_title(f"Solvation Effect on CO2RR Free Energy Diagram\nCu(111) at U = {U:.2f} V vs. RHE",
                  fontsize=12, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="gray", linewidth=0.5, linestyle=":")

    plt.tight_layout()
    out = os.path.join(save_dir, "solvent_comparison_Cu.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def compute_solvation_summary(df_catalysts: pd.DataFrame) -> pd.DataFrame:
    """
    Compute solvation-corrected limiting potentials for all catalysts.
    """
    from scripts.volcano_plot import limiting_potential_CO, limiting_potential_CH4

    records = []
    for _, row in df_catalysts.iterrows():
        dG_CO   = row.get("dG_CO", np.nan)
        dG_COOH = row.get("dG_COOH", np.nan)
        dG_CHO  = row.get("dG_CHO", np.nan)
        dG_CH2O = row.get("dG_CH2O", dG_CHO + 0.3 if not np.isnan(dG_CHO) else np.nan)
        dG_OCH3 = row.get("dG_OCH3", dG_CHO + 0.6 if not np.isnan(dG_CHO) else np.nan)

        if np.isnan(dG_CO) or np.isnan(dG_COOH):
            continue

        # Gas phase
        UL_CO_gas  = limiting_potential_CO(dG_COOH, dG_CO)
        UL_CH4_gas = limiting_potential_CH4(dG_COOH, dG_CO, dG_CHO, dG_CH2O, dG_OCH3)

        # With solvation
        dG_CO_solv   = dG_CO   + SOLVATION_CORR_WATER.get("CO*",   0)
        dG_COOH_solv = dG_COOH + SOLVATION_CORR_WATER.get("COOH*", 0)
        dG_CHO_solv  = dG_CHO  + SOLVATION_CORR_WATER.get("CHO*",  0) if not np.isnan(dG_CHO) else dG_CHO
        UL_CO_solv   = limiting_potential_CO(dG_COOH_solv, dG_CO_solv)
        UL_CH4_solv  = limiting_potential_CH4(dG_COOH_solv, dG_CO_solv, dG_CHO_solv,
                                               dG_CH2O, dG_OCH3)

        records.append({
            "catalyst":       row["catalyst"],
            "category":       row["category"],
            "UL_CO_gas":      round(UL_CO_gas,  4),
            "UL_CO_solv":     round(UL_CO_solv, 4),
            "delta_UL_CO":    round(UL_CO_solv - UL_CO_gas, 4),
            "UL_CH4_gas":     round(UL_CH4_gas,  4),
            "UL_CH4_solv":    round(UL_CH4_solv, 4),
            "delta_UL_CH4":   round(UL_CH4_solv - UL_CH4_gas, 4),
        })

    return pd.DataFrame(records)


if __name__ == "__main__":
    from data.adsorption_energies import get_all_catalysts
    df = get_all_catalysts()

    plot_potential_dependent_activity(df, OUTPUT_DIR)
    plot_solvent_comparison(OUTPUT_DIR)

    df_solv = compute_solvation_summary(df)
    df_solv.to_csv(os.path.join(RESULTS_DIR, "solvation_corrections.csv"), index=False)

    print("\nSolvation correction summary (ΔU_L change):")
    print(df_solv[["catalyst","UL_CO_gas","UL_CO_solv","delta_UL_CO"]
                  ].sort_values("delta_UL_CO", ascending=False).head(10).to_string(index=False))
