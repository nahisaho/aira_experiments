"""
CO2RR Scaling Relations Analysis
Establishes linear scaling relationships (LSRs) between adsorption energies
of CO2RR intermediates on metal surfaces.

Key scaling relations:
  ΔG(*COOH) = a * ΔG(*CO) + b   (primary CO2→CO descriptor)
  ΔG(*CHO)  = a * ΔG(*CO) + b   (CO→CH4 selectivity)
  ΔG(*COH)  vs ΔG(*CHO)         (mechanistic branching)
  Brønsted–Evans–Polanyi (BEP) relations for key steps
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/..")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy import stats
from sklearn.metrics import r2_score

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)


def fit_scaling_relation(x: np.ndarray, y: np.ndarray,
                          name: str = "") -> dict:
    """Fit y = a*x + b using ordinary least squares."""
    slope, intercept, r, p, se = stats.linregress(x, y)
    y_pred = slope * x + intercept
    mae  = np.mean(np.abs(y - y_pred))
    rmse = np.sqrt(np.mean((y - y_pred) ** 2))
    return {
        "relation": name,
        "slope":     round(float(slope), 4),
        "intercept": round(float(intercept), 4),
        "R2":        round(float(r**2), 4),
        "p_value":   float(p),
        "MAE_eV":    round(float(mae), 4),
        "RMSE_eV":   round(float(rmse), 4),
    }


def plot_scaling_relations(df: pd.DataFrame, save_dir: str) -> list:
    """
    Generate scaling relation plots for CO2RR intermediates.
    Returns list of fitted relation dicts.
    """
    # Color map by category
    cat_colors = {
        "pure_metal":  "#1565C0",
        "cu_alloy":    "#E65100",
        "SAC_N-doped": "#2E7D32",
    }
    cat_markers = {
        "pure_metal":  "o",
        "cu_alloy":    "s",
        "SAC_N-doped": "^",
    }

    relations_to_fit = [
        ("dG_CO", "dG_COOH",  "*COOH vs. *CO",       "ΔG(*CO) [eV]", "ΔG(*COOH) [eV]"),
        ("dG_CO", "dG_CHO",   "*CHO vs. *CO",        "ΔG(*CO) [eV]", "ΔG(*CHO) [eV]"),
        ("dG_CO", "dG_COH",   "*COH vs. *CO",        "ΔG(*CO) [eV]", "ΔG(*COH) [eV]"),
        ("dG_CO", "dG_OCCO",  "*OCCO vs. *CO",       "ΔG(*CO) [eV]", "ΔG(*OCCO) [eV]"),
        ("dG_CHO","dG_CH2O",  "*CH2O vs. *CHO",      "ΔG(*CHO) [eV]","ΔG(*CH2O) [eV]"),
        ("dG_COOH","dG_CO",   "*CO vs. *COOH",       "ΔG(*COOH) [eV]","ΔG(*CO) [eV]"),
    ]

    fitted = []
    n_cols = 3
    n_rows = 2
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 10))
    axes = axes.flatten()

    for idx, (xcol, ycol, title, xlabel, ylabel) in enumerate(relations_to_fit):
        ax = axes[idx]
        sub = df.dropna(subset=[xcol, ycol])

        # Plot all points colored by category
        for cat, grp in sub.groupby("category"):
            col = cat_colors.get(cat, "gray")
            mk  = cat_markers.get(cat, "o")
            ax.scatter(grp[xcol], grp[ycol], c=col, marker=mk,
                       s=60, alpha=0.85, label=cat, edgecolors="white", linewidth=0.5)
            # Annotate notable catalysts
            for _, row in grp.iterrows():
                if row["catalyst"] in ("Cu", "Au", "Ag", "Fe-N4C", "Co-N4C",
                                        "CuNi", "Ni", "CuPd", "Cu-N4C"):
                    ax.annotate(row["catalyst"],
                                (row[xcol], row[ycol]),
                                fontsize=7, xytext=(4, 3),
                                textcoords="offset points", color=col,
                                path_effects=[pe.withStroke(linewidth=2,
                                                            foreground="white")])

        # Fit and plot regression line
        x_arr = sub[xcol].values
        y_arr = sub[ycol].values
        fit = fit_scaling_relation(x_arr, y_arr, title)
        fitted.append(fit)

        x_range = np.linspace(x_arr.min() - 0.1, x_arr.max() + 0.1, 100)
        y_range = fit["slope"] * x_range + fit["intercept"]

        # 95% confidence band
        n  = len(x_arr)
        se = np.sqrt(np.sum((y_arr - (fit["slope"]*x_arr + fit["intercept"]))**2) / (n-2))
        x_mean = x_arr.mean()
        ci = 1.96 * se * np.sqrt(1/n + (x_range - x_mean)**2 / np.sum((x_arr-x_mean)**2))
        ax.fill_between(x_range, y_range - ci, y_range + ci,
                         alpha=0.12, color="steelblue")
        ax.plot(x_range, y_range, "--", color="steelblue", linewidth=2,
                label=f"y={fit['slope']:.2f}x+{fit['intercept']:.2f} (R²={fit['R2']:.3f})")

        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.legend(fontsize=7, loc="best")
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color="gray", linewidth=0.5, linestyle=":")
        ax.axvline(0, color="gray", linewidth=0.5, linestyle=":")

    plt.suptitle("CO2RR Linear Scaling Relations", fontsize=14, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(save_dir, "scaling_relations.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")
    return fitted


def plot_selectivity_map(df: pd.DataFrame, save_dir: str) -> None:
    """
    ΔG(*CHO) vs. ΔG(*CO) selectivity map: CO vs CH4/C2H4
    Regions:
      ΔG(*CO) > 0        → CO desorbs (CO product)
      ΔG(*CHO) < 0.4 eV → further reduction (CH4/C2)
    """
    fig, ax = plt.subplots(figsize=(9, 7))
    cat_colors = {
        "pure_metal":  "#1565C0",
        "cu_alloy":    "#E65100",
        "SAC_N-doped": "#2E7D32",
    }
    cat_markers = {
        "pure_metal":  "o",
        "cu_alloy":    "s",
        "SAC_N-doped": "^",
    }

    # Background regions
    ax.axhspan(-3, 0.4,  alpha=0.08, color="#FF9800",
                label="C2+/CH4 selective region")
    ax.axvspan(-0.1, 1.5, alpha=0.08, color="#4CAF50",
                label="CO selective region")

    for cat, grp in df.groupby("category"):
        col = cat_colors.get(cat, "gray")
        mk  = cat_markers.get(cat, "o")
        ax.scatter(grp["dG_CO"], grp["dG_CHO"], c=col, marker=mk,
                   s=80, alpha=0.9, label=cat, edgecolors="white", linewidth=0.5)
        for _, row in grp.iterrows():
            if not pd.isna(row["dG_CO"]) and not pd.isna(row["dG_CHO"]):
                ax.annotate(row["catalyst"],
                            (row["dG_CO"], row["dG_CHO"]),
                            fontsize=7.5, xytext=(4, 3),
                            textcoords="offset points", color=col,
                            path_effects=[pe.withStroke(linewidth=2, foreground="white")])

    # Ideal catalyst region
    ax.add_patch(plt.Circle((-0.65, 0.35), 0.15, color="gold",
                              fill=True, alpha=0.4, linewidth=2,
                              label="Ideal CO2RR catalyst region"))
    ax.text(-0.65, 0.35, "★", ha="center", va="center",
             fontsize=14, color="darkgoldenrod")

    ax.axhline(0.4,  color="darkorange", linestyle="--", linewidth=1.2,
               label="ΔG(*CHO) = 0.4 eV threshold")
    ax.axvline(-0.1, color="green",      linestyle="--", linewidth=1.2,
               label="ΔG(*CO) = -0.1 eV threshold")

    ax.set_xlabel("ΔG(*CO) [eV]", fontsize=12)
    ax.set_ylabel("ΔG(*CHO) [eV]", fontsize=12)
    ax.set_title("CO2RR Selectivity Map: *CO vs. *CHO Descriptor",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)

    out = os.path.join(save_dir, "selectivity_map.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")


if __name__ == "__main__":
    from data.adsorption_energies import get_all_catalysts
    df = get_all_catalysts()

    fitted = plot_scaling_relations(df, OUTPUT_DIR)
    plot_selectivity_map(df, OUTPUT_DIR)

    df_fit = pd.DataFrame(fitted)
    df_fit.to_csv(os.path.join(RESULTS_DIR, "scaling_relations.csv"), index=False)
    print("\nScaling relations summary:")
    print(df_fit[["relation","slope","intercept","R2","MAE_eV"]].to_string(index=False))
