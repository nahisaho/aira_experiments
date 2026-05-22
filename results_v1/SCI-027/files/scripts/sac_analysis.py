"""
Single-Atom Catalyst (SAC) Metal-Support Interaction Analysis
Analyzes the electronic structure and coordination effects of SACs on N-doped carbon.

Key analyses:
  1. d-band center vs. adsorption energy correlations (d-band model, Hammer-Norskov)
  2. Charge transfer between metal atom and N-C support
  3. Coordination number effect (N4 vs. N2 vs. N1 coordination)
  4. Pyridinic vs. pyrrolic vs. graphitic N coordination
  5. Activity/selectivity heatmap for SAC library
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.gridspec import GridSpec
from scipy import stats
import seaborn as sns

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# ------------------------------------------------------------------
# Extended SAC data including electronic structure descriptors
# ------------------------------------------------------------------
SAC_EXTENDED = {
    # d-band center (eV), charge transfer (e), N coord, Bader charge (e+)
    # Adsorption energies for *CO (eV)
    "Fe-N4C":  {"d_center": -1.85, "charge_transfer": 0.72, "N_coord": 4,
                "Bader_M": +1.42, "dG_CO": -1.52, "dG_COOH": -0.22,
                "dG_CHO": -0.52, "HOMO_LUMO_gap": 0.35,
                "spin_moment": 2.12, "metal": "Fe", "coord_type": "N4-porphyrin"},
    "Co-N4C":  {"d_center": -1.52, "charge_transfer": 0.65, "N_coord": 4,
                "Bader_M": +1.28, "dG_CO": -1.42, "dG_COOH": -0.15,
                "dG_CHO": -0.42, "HOMO_LUMO_gap": 0.42,
                "spin_moment": 1.08, "metal": "Co", "coord_type": "N4-porphyrin"},
    "Ni-N4C":  {"d_center": -1.28, "charge_transfer": 0.58, "N_coord": 4,
                "Bader_M": +1.15, "dG_CO": -1.28, "dG_COOH": -0.05,
                "dG_CHO": -0.25, "HOMO_LUMO_gap": 0.55,
                "spin_moment": 0.0,  "metal": "Ni", "coord_type": "N4-porphyrin"},
    "Cu-N4C":  {"d_center": -2.67, "charge_transfer": 0.45, "N_coord": 4,
                "Bader_M": +0.92, "dG_CO": -0.82, "dG_COOH": 0.28,
                "dG_CHO": 0.25, "HOMO_LUMO_gap": 0.78,
                "spin_moment": 0.48, "metal": "Cu", "coord_type": "N4-porphyrin"},
    "Zn-N4C":  {"d_center": -7.52, "charge_transfer": 0.35, "N_coord": 4,
                "Bader_M": +0.78, "dG_CO": -0.32, "dG_COOH": 0.72,
                "dG_CHO": 1.22,  "HOMO_LUMO_gap": 1.25,
                "spin_moment": 0.0,  "metal": "Zn", "coord_type": "N4-porphyrin"},
    "Mn-N4C":  {"d_center": -1.62, "charge_transfer": 0.78, "N_coord": 4,
                "Bader_M": +1.55, "dG_CO": -1.62, "dG_COOH": -0.32,
                "dG_CHO": -0.62, "HOMO_LUMO_gap": 0.28,
                "spin_moment": 3.25, "metal": "Mn", "coord_type": "N4-porphyrin"},
    "Mo-N4C":  {"d_center": -1.15, "charge_transfer": 0.85, "N_coord": 4,
                "Bader_M": +1.72, "dG_CO": -1.72, "dG_COOH": -0.38,
                "dG_CHO": -0.68, "HOMO_LUMO_gap": 0.22,
                "spin_moment": 1.85, "metal": "Mo", "coord_type": "N4-porphyrin"},
    "Cr-N4C":  {"d_center": -0.95, "charge_transfer": 0.88, "N_coord": 4,
                "Bader_M": +1.78, "dG_CO": -1.78, "dG_COOH": -0.42,
                "dG_CHO": -0.72, "HOMO_LUMO_gap": 0.18,
                "spin_moment": 2.68, "metal": "Cr", "coord_type": "N4-porphyrin"},
    "Pd-N4C":  {"d_center": -1.78, "charge_transfer": 0.55, "N_coord": 4,
                "Bader_M": +1.08, "dG_CO": -0.98, "dG_COOH": 0.18,
                "dG_CHO": 0.12,  "HOMO_LUMO_gap": 0.62,
                "spin_moment": 0.0,  "metal": "Pd", "coord_type": "N4-porphyrin"},
    "Ru-N4C":  {"d_center": -1.42, "charge_transfer": 0.68, "N_coord": 4,
                "Bader_M": +1.35, "dG_CO": -1.48, "dG_COOH": -0.18,
                "dG_CHO": -0.48, "HOMO_LUMO_gap": 0.38,
                "spin_moment": 0.72, "metal": "Ru", "coord_type": "N4-porphyrin"},
    # N2-coordination (lower N count)
    "Fe-N2C":  {"d_center": -1.62, "charge_transfer": 0.55, "N_coord": 2,
                "Bader_M": +1.15, "dG_CO": -1.35, "dG_COOH": -0.08,
                "dG_CHO": -0.35, "HOMO_LUMO_gap": 0.45,
                "spin_moment": 2.45, "metal": "Fe", "coord_type": "N2-defect"},
    "Co-N2C":  {"d_center": -1.35, "charge_transfer": 0.48, "N_coord": 2,
                "Bader_M": +0.98, "dG_CO": -1.22, "dG_COOH": 0.02,
                "dG_CHO": -0.22, "HOMO_LUMO_gap": 0.52,
                "spin_moment": 1.25, "metal": "Co", "coord_type": "N2-defect"},
    "Cu-N2C":  {"d_center": -2.45, "charge_transfer": 0.32, "N_coord": 2,
                "Bader_M": +0.65, "dG_CO": -0.65, "dG_COOH": 0.42,
                "dG_CHO": 0.42,  "HOMO_LUMO_gap": 0.88,
                "spin_moment": 0.62, "metal": "Cu", "coord_type": "N2-defect"},
    # Pyridinic vs. pyrrolic
    "Fe-N4C-pyr":{"d_center": -1.78, "charge_transfer": 0.68, "N_coord": 4,
                  "Bader_M": +1.35, "dG_CO": -1.42, "dG_COOH": -0.18,
                  "dG_CHO": -0.45, "HOMO_LUMO_gap": 0.40,
                  "spin_moment": 2.25, "metal": "Fe", "coord_type": "N4-pyrrolic"},
}


def plot_dband_model(df: pd.DataFrame, save_dir: str) -> None:
    """
    Plot d-band center vs. *CO adsorption energy (Hammer-Norskov d-band model).
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    metals = df["metal"].unique()
    metal_colors = {
        "Fe": "#C62828", "Co": "#1565C0", "Ni": "#2E7D32",
        "Cu": "#E65100", "Mn": "#6A1B9A", "Mo": "#00695C",
        "Cr": "#4527A0", "Pd": "#AD1457", "Ru": "#558B2F",
        "Zn": "#795548",
    }

    # Left: d-band center vs. dG_CO
    ax = axes[0]
    for _, row in df.iterrows():
        col = metal_colors.get(row["metal"], "gray")
        ax.scatter(row["d_center"], row["dG_CO"],
                   c=col, s=100, alpha=0.9, edgecolors="white", zorder=5)
        ax.annotate(row["SAC"],
                    (row["d_center"], row["dG_CO"]),
                    fontsize=7.5, xytext=(4, 3),
                    textcoords="offset points", color=col,
                    path_effects=[pe.withStroke(linewidth=2, foreground="white")])

    # Fit d-band model
    x = df["d_center"].values
    y = df["dG_CO"].values
    slope, intercept, r, p, se = stats.linregress(x, y)
    x_fit = np.linspace(x.min()-0.2, x.max()+0.2, 100)
    ax.plot(x_fit, slope*x_fit + intercept, "--", color="steelblue",
             linewidth=2, label=f"y={slope:.2f}x+{intercept:.2f}, R²={r**2:.3f}")
    ax.set_xlabel("d-band center (eV)", fontsize=11)
    ax.set_ylabel("ΔG(*CO) [eV]", fontsize=11)
    ax.set_title("d-band Center vs. *CO Binding Energy\n(Hammer-Nørskov Model)",
                  fontsize=11, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.axhline(-0.65, color="orange", linestyle=":", label="Optimal ΔG(*CO)")

    # Right: Bader charge vs. limiting potential for CO
    ax = axes[1]
    from scripts.volcano_plot import limiting_potential_CO
    df["UL_CO"] = df.apply(
        lambda r: limiting_potential_CO(r["dG_COOH"], r["dG_CO"]), axis=1
    )

    for _, row in df.iterrows():
        col = metal_colors.get(row["metal"], "gray")
        mk  = "o" if row["N_coord"] == 4 else "s" if row["N_coord"] == 2 else "^"
        ax.scatter(row["Bader_M"], row["UL_CO"],
                   c=col, marker=mk, s=100, alpha=0.9, edgecolors="white", zorder=5)
        ax.annotate(row["SAC"],
                    (row["Bader_M"], row["UL_CO"]),
                    fontsize=7.5, xytext=(4, 3),
                    textcoords="offset points", color=col,
                    path_effects=[pe.withStroke(linewidth=2, foreground="white")])

    ax.set_xlabel("Bader Charge on Metal [e⁺]", fontsize=11)
    ax.set_ylabel("Limiting Potential U_L [V vs. RHE]", fontsize=11)
    ax.set_title("Metal Oxidation State vs. CO2RR Activity\n(SAC on N-doped Carbon)",
                  fontsize=11, fontweight="bold")
    ax.grid(True, alpha=0.3)
    # N-coord legend
    from matplotlib.lines import Line2D
    legend_el = [Line2D([0],[0], marker="o", color="gray", ls="None",
                          markersize=8, label="N4 coordination"),
                  Line2D([0],[0], marker="s", color="gray", ls="None",
                          markersize=8, label="N2 coordination")]
    ax.legend(handles=legend_el, fontsize=9)

    plt.suptitle("SAC Metal-Support Interaction Analysis",
                  fontsize=13, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(save_dir, "sac_metal_support.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def plot_coordination_effect(df: pd.DataFrame, save_dir: str) -> None:
    """
    Compare N4 vs N2 coordination for same metal (Fe, Co, Cu).
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    metals_compare = ["Fe", "Co", "Cu"]

    for ax, metal in zip(axes, metals_compare):
        sub = df[df["metal"] == metal].copy()
        intermediates = ["dG_CO", "dG_COOH", "dG_CHO"]
        labels = ["*CO", "*COOH", "*CHO"]
        x = np.arange(len(labels))
        width = 0.25

        colors = plt.cm.Set2(np.linspace(0, 1, len(sub)))
        for i, (_, row) in enumerate(sub.iterrows()):
            vals = [row[k] for k in intermediates]
            bars = ax.bar(x + i*width, vals, width, label=row["SAC"],
                           color=colors[i], alpha=0.85, edgecolor="white")

        ax.set_xticks(x + width*(len(sub)-1)/2)
        ax.set_xticklabels(labels, fontsize=10)
        ax.set_title(f"{metal} SAC — Coordination Effect", fontsize=10, fontweight="bold")
        ax.set_ylabel("ΔG [eV]" if ax == axes[0] else "")
        ax.axhline(0, color="gray", linewidth=0.8, linestyle=":")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2, axis="y")

    plt.suptitle("N-coordination Effect on CO2RR Intermediate Binding Energies",
                  fontsize=12, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(save_dir, "sac_coordination_effect.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


def plot_sac_heatmap(df: pd.DataFrame, save_dir: str) -> None:
    """
    Heatmap: catalysts × descriptors (binding energies + limiting potential).
    """
    from scripts.volcano_plot import limiting_potential_CO, limiting_potential_CH4

    df["UL_CO"]  = df.apply(lambda r: limiting_potential_CO(r["dG_COOH"], r["dG_CO"]), axis=1)
    df["UL_CH4"] = df.apply(
        lambda r: limiting_potential_CH4(
            r["dG_COOH"], r["dG_CO"], r["dG_CHO"],
            r.get("dG_CH2O", r["dG_CHO"]+0.3),
            r.get("dG_OCH3", r["dG_CHO"]+0.6)
        ), axis=1
    )

    cols_for_heatmap = ["d_center", "Bader_M", "dG_CO", "dG_COOH",
                         "dG_CHO", "UL_CO", "UL_CH4", "HOMO_LUMO_gap", "spin_moment"]
    col_labels       = ["d-band\ncenter", "Bader\ncharge", "ΔG(*CO)",
                         "ΔG(*COOH)", "ΔG(*CHO)", "U_L(CO)", "U_L(CH4)",
                         "HOMO-LUMO\ngap", "Spin\nmoment"]

    heat_data = df.set_index("SAC")[cols_for_heatmap]

    fig, ax = plt.subplots(figsize=(13, 8))
    # Normalize each column to [-1,1] for color scale
    heat_norm = (heat_data - heat_data.mean()) / (heat_data.std() + 1e-8)
    sns.heatmap(heat_norm, ax=ax, cmap="RdBu_r", center=0,
                 annot=heat_data.round(2), fmt=".2f",
                 annot_kws={"size": 8}, linewidths=0.5,
                 xticklabels=col_labels, yticklabels=True,
                 cbar_kws={"label": "Normalized Value"})
    ax.set_title("SAC Property Heatmap (N-doped Carbon Support)",
                  fontsize=12, fontweight="bold")
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelsize=9)
    ax.tick_params(axis="y", labelsize=9)
    plt.tight_layout()
    out = os.path.join(save_dir, "sac_heatmap.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


if __name__ == "__main__":
    df = pd.DataFrame([
        {"SAC": k, **v} for k, v in SAC_EXTENDED.items()
    ])
    df.to_csv(os.path.join(RESULTS_DIR, "sac_analysis.csv"), index=False)

    plot_dband_model(df, OUTPUT_DIR)
    plot_coordination_effect(df, OUTPUT_DIR)
    plot_sac_heatmap(df, OUTPUT_DIR)

    print("\nSAC Summary (sorted by UL_CO):")
    from scripts.volcano_plot import limiting_potential_CO
    df["UL_CO"] = df.apply(
        lambda r: limiting_potential_CO(r["dG_COOH"], r["dG_CO"]), axis=1
    )
    print(df.sort_values("UL_CO", ascending=False)[
        ["SAC","metal","N_coord","d_center","Bader_M","dG_CO","UL_CO"]
    ].to_string(index=False))
