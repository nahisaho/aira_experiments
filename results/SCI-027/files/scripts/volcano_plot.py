"""
CO2RR Volcano Plot Generator
Uses the computational hydrogen electrode (CHE) model to generate
volcano plots for CO (2e-) and CH4 (8e-) production activity.

Volcano plot concept (Norskov/Bagger framework):
  - x-axis: ΔG(*CO) binding energy descriptor
  - y-axis: limiting potential U_L (the most negative step ΔG/e)
  - Optimal catalyst: minimizes the largest uphill step
  - Ideal *CO binding: neither too strong (poisoning) nor too weak (poor activation)
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# ZPE-TS corrections (condensed from reaction_pathways.py)
ZPE_TS = {"COOH*": 0.40, "CO*": 0.12, "CHO*": 0.32, "CH2O*": 0.48, "OCH3*": 0.60}

# Equilibrium potentials (V vs. RHE)
U_EQ = {
    "CO":   -0.106,   # CO2 + 2H+ + 2e- → CO + H2O
    "CH4":  +0.169,   # CO2 + 8H+ + 8e- → CH4 + 2H2O
    "C2H4": +0.064,   # 2CO2+12H++12e- → C2H4+4H2O
    "HCOOH":-0.250,   # CO2 + 2H+ + 2e- → HCOOH
}


def limiting_potential_CO(dG_COOH: float, dG_CO: float) -> float:
    """
    Limiting potential for CO production.
    Steps: CO2→COOH* (ΔG1), COOH*→CO* (ΔG2), CO*→CO(g) (ΔG3)
    """
    dG1 = dG_COOH + ZPE_TS["COOH*"]
    dG2 = (dG_CO + ZPE_TS["CO*"] + 0.27) - (dG_COOH + ZPE_TS["COOH*"])
    dG3 = -0.517 - (dG_CO + ZPE_TS["CO*"])   # CO(g) ref
    return -max(dG1, dG2, max(dG3, 0))


def limiting_potential_CH4(dG_COOH: float, dG_CO: float,
                            dG_CHO: float, dG_CH2O: float,
                            dG_OCH3: float) -> float:
    """
    Limiting potential for CH4 production (8e- pathway via *CO → *CHO → *CH2O → *OCH3 → CH4).
    """
    # CO2 → CO* : 2 electrons
    dG1 = dG_COOH + ZPE_TS["COOH*"]
    dG2 = (dG_CO + ZPE_TS["CO*"] + 0.27) - (dG_COOH + ZPE_TS["COOH*"])
    # CO* → CHO* → CH2O* → OCH3* → CH4 : 6 more electrons
    dG3 = (dG_CHO  + ZPE_TS["CHO*"])  - (dG_CO  + ZPE_TS["CO*"])
    dG4 = (dG_CH2O + ZPE_TS["CH2O*"]) - (dG_CHO + ZPE_TS["CHO*"])
    dG5 = (dG_OCH3 + ZPE_TS["OCH3*"]) - (dG_CH2O + ZPE_TS["CH2O*"])
    dG6 = -1.060 - 8*U_EQ["CH4"] - (dG_OCH3 + ZPE_TS["OCH3*"])   # CH4(g)
    return -max(dG1, dG2, dG3, dG4, dG5, max(dG6, 0))


def limiting_potential_C2H4(dG_CO: float, dG_OCCO: float,
                              dG_OCCHO: float) -> float:
    """
    Limiting potential for C2H4 production via CO* dimerization.
    """
    if dG_OCCO is None or np.isnan(dG_OCCO):
        return np.nan
    dG_coupling = dG_OCCO - 2*dG_CO           # C-C coupling step
    dG_hydro    = dG_OCCHO - dG_OCCO          # first hydrogenation
    dG_final    = -1.320 - 12*U_EQ["C2H4"] - dG_OCCHO
    return -max(dG_coupling, dG_hydro, max(dG_final, 0))


# ------------------------------------------------------------------
# Analytical volcano curve (for smooth background)
# ------------------------------------------------------------------
def volcano_CO_analytical(dG_CO_range: np.ndarray) -> np.ndarray:
    """
    Analytical volcano using *COOH scaling:  ΔG(*COOH) ≈ 0.87*ΔG(*CO) + 0.63
    Left leg  (strong binding, *CO poisoning): U_L = -ΔG(*CO) - ZPE_TS["CO*"]
    Right leg (weak binding, *COOH limit):    U_L = -(ΔG(*COOH) + ZPE_TS["COOH*"])
    """
    dG_COOH_scaled = 0.87 * dG_CO_range + 0.63
    U_left  = -(abs(dG_CO_range) + ZPE_TS["CO*"])       # poisoning side
    U_right = -(dG_COOH_scaled + ZPE_TS["COOH*"])        # activation side
    return np.maximum(U_left, U_right)   # take the LESS negative (higher activity)


def volcano_CH4_analytical(dG_CO_range: np.ndarray) -> np.ndarray:
    """
    Analytical volcano for CH4: branching at *CHO step.
    *CHO scaling: ΔG(*CHO) ≈ 0.76*ΔG(*CO) + 1.09
    Left leg:  -abs(ΔG(*CO))
    Right leg: -(ΔG(*CHO))
    """
    dG_CHO_scaled = 0.76 * dG_CO_range + 1.09
    U_left  = -(abs(dG_CO_range) + ZPE_TS["CO*"])
    U_right = -(dG_CHO_scaled + ZPE_TS["CHO*"])
    return np.maximum(U_left, U_right)


def build_volcano_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Compute limiting potentials for all catalysts."""
    records = []
    for _, row in df.iterrows():
        cat = row["catalyst"]
        dG_CO    = row.get("dG_CO",    np.nan)
        dG_COOH  = row.get("dG_COOH",  np.nan)
        dG_CHO   = row.get("dG_CHO",   np.nan)
        dG_CH2O  = row.get("dG_CH2O",  np.nan)
        dG_OCH3  = row.get("dG_OCH3",  np.nan)
        dG_OCCO  = row.get("dG_OCCO",  np.nan)
        dG_OCCHO = row.get("dG_OCCHO", np.nan)

        UL_CO   = limiting_potential_CO(dG_COOH, dG_CO)
        UL_CH4  = limiting_potential_CH4(dG_COOH, dG_CO, dG_CHO, dG_CH2O, dG_OCH3)
        UL_C2H4 = limiting_potential_C2H4(dG_CO, dG_OCCO, dG_OCCHO)

        records.append({
            "catalyst":     cat,
            "category":     row["category"],
            "dG_CO":        dG_CO,
            "dG_COOH":      dG_COOH,
            "dG_CHO":       dG_CHO,
            "UL_CO_V":      round(UL_CO, 4),
            "UL_CH4_V":     round(UL_CH4, 4),
            "UL_C2H4_V":    round(float(UL_C2H4) if not np.isnan(UL_C2H4) else np.nan, 4),
            "overpotential_CO":   round(abs(UL_CO)  - abs(U_EQ["CO"]),  4),
            "overpotential_CH4":  round(abs(UL_CH4) - abs(U_EQ["CH4"]), 4),
        })

    return pd.DataFrame(records).sort_values("UL_CO_V", ascending=False)


def plot_volcano(df_vol: pd.DataFrame, product: str,
                 UL_col: str, save_path: str) -> None:
    """
    Generate a volcano plot for a given product.
    """
    fig, ax = plt.subplots(figsize=(10, 7))

    x_range = np.linspace(-2.0, 1.5, 400)
    if product == "CO":
        y_volcano = volcano_CO_analytical(x_range)
        eq_U = U_EQ["CO"]
        color_volcano = "#1976D2"
    elif product == "CH4":
        y_volcano = volcano_CH4_analytical(x_range)
        eq_U = U_EQ["CH4"]
        color_volcano = "#C62828"
    else:
        y_volcano = volcano_CO_analytical(x_range)
        eq_U = U_EQ.get(product, 0)
        color_volcano = "#388E3C"

    # Volcano background
    ax.fill_between(x_range, y_volcano, y_volcano.min() - 0.5,
                     alpha=0.07, color=color_volcano)
    ax.plot(x_range, y_volcano, "-", color=color_volcano,
             linewidth=2.5, alpha=0.8, label=f"Analytical {product} volcano")

    # Equilibrium potential line
    ax.axhline(eq_U, color="green", linestyle=":", linewidth=1.5,
                label=f"U_eq = {eq_U:.3f} V")

    cat_colors  = {"pure_metal": "#1565C0", "cu_alloy": "#E65100", "SAC_N-doped": "#2E7D32"}
    cat_markers = {"pure_metal": "o",       "cu_alloy": "s",       "SAC_N-doped": "^"}
    cat_sizes   = {"pure_metal": 90,        "cu_alloy": 90,        "SAC_N-doped": 100}

    sub = df_vol.dropna(subset=["dG_CO", UL_col])
    for cat_type, grp in sub.groupby("category"):
        col = cat_colors.get(cat_type, "gray")
        mk  = cat_markers.get(cat_type, "o")
        sz  = cat_sizes.get(cat_type, 80)
        ax.scatter(grp["dG_CO"], grp[UL_col],
                   c=col, marker=mk, s=sz, alpha=0.9,
                   edgecolors="white", linewidth=0.8,
                   label=cat_type, zorder=5)

        # Label points near volcano peak (top performers)
        peak_x = x_range[np.argmax(y_volcano)]
        for _, row in grp.iterrows():
            dist_to_peak = abs(row["dG_CO"] - peak_x)
            ul_val = row[UL_col]
            # Always label well-known or top-performing catalysts
            always_label = {"Cu", "Au", "Ag", "Ni", "Fe", "Fe-N4C", "Co-N4C",
                             "Cu-N4C", "CuNi", "CuPd", "CuCo", "Cu3Zn"}
            if dist_to_peak < 0.6 or row["catalyst"] in always_label:
                ax.annotate(row["catalyst"],
                            (row["dG_CO"], ul_val),
                            fontsize=8, xytext=(5, 4),
                            textcoords="offset points", color=col,
                            path_effects=[pe.withStroke(linewidth=2.5,
                                                        foreground="white")],
                            zorder=6)

    # Highlight optimal region
    peak_idx = np.argmax(y_volcano)
    ax.axvspan(x_range[peak_idx]-0.25, x_range[peak_idx]+0.25,
                alpha=0.1, color="gold", label="Optimal *CO binding window")

    ax.set_xlabel("ΔG(*CO) [eV]", fontsize=12)
    ax.set_ylabel("Limiting Potential U_L [V vs. RHE]", fontsize=12)
    ax.set_title(f"CO2RR Volcano Plot — {product} Production",
                  fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-2.0, 1.5)

    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {save_path}")


def plot_combined_volcano(df_vol: pd.DataFrame, save_path: str) -> None:
    """
    Side-by-side CO and CH4 volcano plots.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    products = [
        ("CO",  "UL_CO_V",  "#1976D2", axes[0]),
        ("CH4", "UL_CH4_V", "#C62828", axes[1]),
    ]

    cat_colors  = {"pure_metal": "#1565C0", "cu_alloy": "#E65100", "SAC_N-doped": "#2E7D32"}
    cat_markers = {"pure_metal": "o",       "cu_alloy": "s",       "SAC_N-doped": "^"}

    for prod, ul_col, vcol, ax in products:
        x_range = np.linspace(-2.0, 1.5, 400)
        y_vol   = volcano_CO_analytical(x_range) if prod == "CO" \
                  else volcano_CH4_analytical(x_range)

        ax.fill_between(x_range, y_vol, y_vol.min()-0.3,
                         alpha=0.07, color=vcol)
        ax.plot(x_range, y_vol, "-", color=vcol, linewidth=2.5,
                 label=f"{prod} volcano")
        ax.axhline(U_EQ[prod], color="green", linestyle=":",
                    linewidth=1.5, label=f"U_eq = {U_EQ[prod]:.2f} V")

        sub = df_vol.dropna(subset=["dG_CO", ul_col])
        for ct, grp in sub.groupby("category"):
            ax.scatter(grp["dG_CO"], grp[ul_col],
                       c=cat_colors.get(ct, "gray"),
                       marker=cat_markers.get(ct, "o"),
                       s=75, alpha=0.9, edgecolors="white",
                       linewidth=0.8, label=ct)
            for _, row in grp.iterrows():
                if row["catalyst"] in {"Cu", "Au", "Ag", "Fe-N4C", "Co-N4C",
                                        "CuNi", "Cu3Zn", "CuPd"}:
                    ax.annotate(row["catalyst"],
                                (row["dG_CO"], row[ul_col]),
                                fontsize=8, xytext=(5, 4),
                                textcoords="offset points",
                                color=cat_colors.get(ct, "gray"),
                                path_effects=[pe.withStroke(linewidth=2,
                                                            foreground="white")])

        ax.set_xlabel("ΔG(*CO) [eV]", fontsize=11)
        ax.set_ylabel("U_L [V vs. RHE]", fontsize=11)
        ax.set_title(f"{prod} Production Volcano", fontsize=12, fontweight="bold")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-2.0, 1.5)

    plt.suptitle("CO2RR Volcano Plots: CO vs. CH4 Production Activity",
                  fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {save_path}")


if __name__ == "__main__":
    from data.adsorption_energies import get_all_catalysts
    df = get_all_catalysts()
    df_vol = build_volcano_dataframe(df)
    df_vol.to_csv(os.path.join(RESULTS_DIR, "volcano_data.csv"), index=False)

    plot_combined_volcano(df_vol,
        os.path.join(OUTPUT_DIR, "volcano_combined.png"))
    plot_volcano(df_vol, "CO",  "UL_CO_V",
        os.path.join(OUTPUT_DIR, "volcano_CO.png"))
    plot_volcano(df_vol, "CH4", "UL_CH4_V",
        os.path.join(OUTPUT_DIR, "volcano_CH4.png"))

    print("\nTop 5 catalysts for CO production (highest U_L):")
    print(df_vol.nlargest(5, "UL_CO_V")[
        ["catalyst","category","dG_CO","UL_CO_V","overpotential_CO"]
    ].to_string(index=False))
    print("\nTop 5 catalysts for CH4 production:")
    print(df_vol.nlargest(5, "UL_CH4_V")[
        ["catalyst","category","dG_CO","UL_CH4_V","overpotential_CH4"]
    ].to_string(index=False))
