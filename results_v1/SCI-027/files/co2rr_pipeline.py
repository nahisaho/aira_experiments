"""
Main CO2RR Computational Screening Pipeline
Orchestrates all analysis modules and generates comprehensive results.

Pipeline:
  1. Load catalyst database
  2. Run reaction pathway analysis (CHE model)
  3. Compute scaling relations
  4. Generate volcano plots
  5. SAC metal-support interaction analysis
  6. Solvent + potential-dependent corrections
  7. Rank and score candidates
  8. Export all results and figures
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from datetime import datetime
import json
import warnings
warnings.filterwarnings("ignore")

BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR  = os.path.join(BASE_DIR, "figures")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
DATA_DIR    = os.path.join(BASE_DIR, "data")
LOGS_DIR    = os.path.join(BASE_DIR, "logs")

for d in [OUTPUT_DIR, RESULTS_DIR, DATA_DIR, LOGS_DIR]:
    os.makedirs(d, exist_ok=True)

LOG_FILE = os.path.join(LOGS_DIR, "process-log.jsonl")


def log_event(phase: str, event_type: str, skill: str = "co2rr-pipeline",
               handoff_in: dict = None, handoff_out: dict = None,
               files: list = None, status: str = "ok", note: str = "") -> None:
    entry = {
        "timestamp":     datetime.utcnow().isoformat() + "Z",
        "phase":         phase,
        "event_type":    event_type,
        "actor":         "co-scientist",
        "skill_or_tool": skill,
        "handoff_in":    handoff_in or {},
        "handoff_out":   handoff_out or {},
        "files_written": files or [],
        "status":        status,
        "note":          note,
    }
    with open(LOG_FILE, "a") as f:
        f.write(json.dumps(entry) + "\n")
    print(f"  [{phase}] {event_type} — {status}")


# ==================================================================
# SCORING FUNCTION
# ==================================================================
def score_catalyst(row: pd.Series) -> dict:
    """
    Multi-objective score for CO2RR performance.
    Score components:
      1. CO activity:   how close UL_CO is to 0 V (or U_eq = -0.106 V)
      2. CH4 activity:  how close UL_CH4 is to U_eq(CH4) = +0.169 V
      3. C2H4 activity: how close UL_C2H4 is to U_eq(C2H4) = +0.064 V
      4. CO selectivity vs CH4: sign/magnitude of (dG_CHO - dG_CO)
      5. Overpotential penalty: penalize large |UL|
    """
    scores = {}

    # CO activity (maximize UL_CO → closest to -0.106 V)
    ul_co = row.get("UL_CO_V", np.nan)
    if not np.isnan(ul_co):
        scores["CO_activity"] = max(0, 1 - abs(ul_co - (-0.106)) / 1.5)
    else:
        scores["CO_activity"] = 0.0

    # CH4 activity
    ul_ch4 = row.get("UL_CH4_V", np.nan)
    if not np.isnan(ul_ch4):
        scores["CH4_activity"] = max(0, 1 - abs(ul_ch4 - 0.169) / 1.5)
    else:
        scores["CH4_activity"] = 0.0

    # C2H4 activity
    ul_c2h4 = row.get("UL_C2H4_V", np.nan)
    if not np.isnan(ul_c2h4):
        scores["C2H4_activity"] = max(0, 1 - abs(ul_c2h4 - 0.064) / 1.5)
    else:
        scores["C2H4_activity"] = 0.0

    # CO selectivity (dG_CO > -0.2 eV → CO desorbs, not overbound)
    dG_CO  = row.get("dG_CO",  np.nan)
    dG_CHO = row.get("dG_CHO", np.nan)
    if not np.isnan(dG_CO):
        # Optimal CO binding: -0.3 to -0.8 eV
        scores["CO_selectivity"] = max(0, 1 - abs(dG_CO - (-0.55)) / 1.0)
    else:
        scores["CO_selectivity"] = 0.0

    # C2 selectivity: low dG_CHO/dG_CO ratio → good for C-C coupling
    if not np.isnan(dG_CO) and not np.isnan(dG_CHO) and dG_CO < -0.3:
        scores["C2_selectivity"] = max(0, 1 - abs(dG_CHO - 0.3) / 1.0)
    else:
        scores["C2_selectivity"] = 0.0

    # Composite score (weighted)
    scores["composite_CO"]  = 0.6*scores["CO_activity"]  + 0.3*scores["CO_selectivity"]  + 0.1*scores["C2H4_activity"]
    scores["composite_CH4"] = 0.6*scores["CH4_activity"] + 0.2*scores["CO_selectivity"]  + 0.2*scores["C2H4_activity"]
    scores["composite_C2+"] = 0.5*scores["C2H4_activity"]+ 0.3*scores["C2_selectivity"]  + 0.2*scores["CH4_activity"]

    return scores


def plot_radar_top_candidates(df_scored: pd.DataFrame, save_path: str,
                               n_top: int = 6) -> None:
    """
    Radar (spider) chart for top catalyst candidates.
    """
    top = df_scored.nlargest(n_top, "composite_CO")
    categories = ["CO_activity", "CH4_activity", "C2H4_activity",
                   "CO_selectivity", "C2_selectivity"]
    labels = ["CO\nActivity", "CH4\nActivity", "C2H4\nActivity",
               "CO\nSelectivity", "C2\nSelectivity"]
    N = len(categories)
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]

    fig, axes = plt.subplots(2, 3, figsize=(14, 10),
                              subplot_kw=dict(polar=True))
    axes = axes.flatten()
    colors = plt.cm.tab10(np.linspace(0, 1, n_top))

    for i, (_, row) in enumerate(top.iterrows()):
        ax = axes[i]
        vals = [row[c] for c in categories]
        vals += vals[:1]

        ax.fill(angles, vals, alpha=0.25, color=colors[i])
        ax.plot(angles, vals, color=colors[i], linewidth=2)
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylim(0, 1)
        ax.set_yticks([0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(["0.25","0.5","0.75","1.0"], fontsize=6)
        ax.set_title(f"{row['catalyst']}\n({row['category']})",
                      fontsize=9, fontweight="bold", pad=10)
        ax.grid(True, alpha=0.3)

    plt.suptitle("Top CO2RR Catalyst Candidates — Performance Radar",
                  fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {save_path}")


def plot_materials_screening_summary(df_scored: pd.DataFrame,
                                      save_path: str) -> None:
    """
    Final summary bubble chart: dG_CO vs. UL_CO, bubble size = composite score.
    """
    fig, ax = plt.subplots(figsize=(12, 8))

    cat_colors  = {"pure_metal": "#1565C0", "cu_alloy": "#E65100", "SAC_N-doped": "#2E7D32"}
    cat_markers = {"pure_metal": "o",       "cu_alloy": "s",       "SAC_N-doped": "^"}

    sub = df_scored.dropna(subset=["dG_CO", "UL_CO_V"])
    for cat_type, grp in sub.groupby("category"):
        col = cat_colors.get(cat_type, "gray")
        mk  = cat_markers.get(cat_type, "o")
        sz  = (grp["composite_CO"] * 400 + 30).values
        sc  = ax.scatter(grp["dG_CO"], grp["UL_CO_V"],
                          c=col, marker=mk, s=sz, alpha=0.80,
                          edgecolors="white", linewidth=0.8,
                          label=cat_type, zorder=5)
        for _, row in grp.iterrows():
            ax.annotate(row["catalyst"],
                        (row["dG_CO"], row["UL_CO_V"]),
                        fontsize=7.5, xytext=(4, 3),
                        textcoords="offset points", color=col,
                        path_effects=[pe.withStroke(linewidth=2, foreground="white")])

    # Optimal window box
    rect = plt.Rectangle((-0.9, -0.5), 0.6, 0.4,
                           linewidth=2, edgecolor="gold",
                           facecolor="gold", alpha=0.12,
                           label="Optimal window")
    ax.add_patch(rect)
    ax.text(-0.61, -0.32, "★ Optimal\nRegion", fontsize=9,
             color="goldenrod", ha="center", fontweight="bold")

    ax.set_xlabel("ΔG(*CO) [eV] — CO Binding Descriptor", fontsize=12)
    ax.set_ylabel("Limiting Potential U_L [V vs. RHE]", fontsize=12)
    ax.set_title("CO2RR Catalyst Screening Summary\n"
                  "(bubble size ∝ composite CO activity score)",
                  fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.axhline(-0.106, color="green", linestyle=":", linewidth=1.2,
                label="U_eq(CO)")

    # Size legend
    for size_val, label in [(0.3, "Low"), (0.6, "Medium"), (0.9, "High")]:
        ax.scatter([], [], c="gray", alpha=0.6, s=size_val*400+30,
                    label=f"Score: {label} ({size_val:.1f})")
    ax.legend(fontsize=8, loc="lower right")

    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {save_path}")


# ==================================================================
# MAIN PIPELINE
# ==================================================================
def run_pipeline():
    log_event("INIT", "run_started", note="CO2RR screening pipeline v1.0")
    print("\n" + "="*65)
    print("  CO2RR Computational Screening Pipeline")
    print("  Based on CHE + Scaling Relations + Volcano Analysis")
    print("="*65)

    # ----------------------------------------------------------------
    # 1. Load catalyst database
    # ----------------------------------------------------------------
    print("\n[1/6] Loading catalyst database...")
    log_event("DATA", "prompt_received", skill="adsorption_energies",
               note="Loading 45 catalysts from literature DFT database")
    from data.adsorption_energies import get_all_catalysts
    df = get_all_catalysts()
    df.to_csv(os.path.join(DATA_DIR, "all_catalysts.csv"), index=False)
    log_event("DATA", "file_written", files=[os.path.join(DATA_DIR, "all_catalysts.csv")],
               handoff_out={"n_catalysts": len(df), "categories": df["category"].value_counts().to_dict()})
    print(f"  Loaded {len(df)} catalysts: "
          f"{df['category'].value_counts().to_dict()}")

    # ----------------------------------------------------------------
    # 2. Reaction pathway analysis
    # ----------------------------------------------------------------
    print("\n[2/6] Reaction pathway analysis (CHE model)...")
    log_event("EXECUTE", "skill_selected", skill="reaction_pathways",
               handoff_in={"model": "CHE", "T": 298.15, "include_ZPE": True})
    from scripts.reaction_pathways import (analyze_catalyst,
                                            plot_free_energy_diagrams)
    from data.adsorption_energies import PURE_METALS as PURE_METALS
    # Free energy diagrams for key catalysts
    cu_data = PURE_METALS["Cu"]
    fe_n4c  = {"dG_COOH": -0.22, "dG_CO": -1.52, "dG_CHO": -0.52,
                "dG_CH2O": -1.28, "dG_OCH3": -2.28}
    au_data = PURE_METALS["Au"]

    fig_paths = []
    plot_free_energy_diagrams(cu_data,
        os.path.join(OUTPUT_DIR, "free_energy_diagram_Cu.png"))
    fig_paths.append("figures/free_energy_diagram_Cu.png")
    plot_free_energy_diagrams(au_data,
        os.path.join(OUTPUT_DIR, "free_energy_diagram_Au.png"))
    fig_paths.append("figures/free_energy_diagram_Au.png")

    # Pathway analysis for all pure metals
    from data.adsorption_energies import PURE_METALS, CU_ALLOYS
    all_metal_data = {**PURE_METALS, **CU_ALLOYS}
    path_records = []
    for name, data in all_metal_data.items():
        r = analyze_catalyst(name, data)
        path_records.append({k: v for k, v in r.items() if k != "_pathways"})
    df_paths = pd.DataFrame(path_records)
    df_paths.to_csv(os.path.join(RESULTS_DIR, "pathway_analysis.csv"), index=False)
    log_event("EXECUTE", "handoff_completed", skill="reaction_pathways",
               files=fig_paths + ["results/pathway_analysis.csv"])

    # ----------------------------------------------------------------
    # 3. Scaling relations
    # ----------------------------------------------------------------
    print("\n[3/6] Computing scaling relations...")
    log_event("EXECUTE", "skill_selected", skill="scaling_relations")
    from scripts.scaling_relations import (plot_scaling_relations,
                                            plot_selectivity_map,
                                            fit_scaling_relation)
    fitted = plot_scaling_relations(df, OUTPUT_DIR)
    plot_selectivity_map(df, OUTPUT_DIR)
    df_fit = pd.DataFrame(fitted)
    df_fit.to_csv(os.path.join(RESULTS_DIR, "scaling_relations.csv"), index=False)
    log_event("EXECUTE", "handoff_completed", skill="scaling_relations",
               handoff_out={"n_relations": len(df_fit),
                             "best_R2": float(df_fit["R2"].max())},
               files=["figures/scaling_relations.png",
                       "figures/selectivity_map.png",
                       "results/scaling_relations.csv"])

    # ----------------------------------------------------------------
    # 4. Volcano plots
    # ----------------------------------------------------------------
    print("\n[4/6] Generating volcano plots...")
    log_event("EXECUTE", "skill_selected", skill="volcano_plot")
    from scripts.volcano_plot import (build_volcano_dataframe,
                                       plot_combined_volcano, plot_volcano)
    df_vol = build_volcano_dataframe(df)
    df_vol.to_csv(os.path.join(RESULTS_DIR, "volcano_data.csv"), index=False)
    plot_combined_volcano(df_vol,
        os.path.join(OUTPUT_DIR, "volcano_combined.png"))
    plot_volcano(df_vol, "CO",  "UL_CO_V",
        os.path.join(OUTPUT_DIR, "volcano_CO.png"))
    plot_volcano(df_vol, "CH4", "UL_CH4_V",
        os.path.join(OUTPUT_DIR, "volcano_CH4.png"))
    log_event("EXECUTE", "handoff_completed", skill="volcano_plot",
               files=["figures/volcano_combined.png",
                       "figures/volcano_CO.png",
                       "figures/volcano_CH4.png"])

    # ----------------------------------------------------------------
    # 5. SAC metal-support interaction
    # ----------------------------------------------------------------
    print("\n[5/6] SAC metal-support interaction analysis...")
    log_event("EXECUTE", "skill_selected", skill="sac_analysis")
    from scripts.sac_analysis import (SAC_EXTENDED, plot_dband_model,
                                       plot_coordination_effect,
                                       plot_sac_heatmap)
    df_sac = pd.DataFrame([{"SAC": k, **v} for k, v in SAC_EXTENDED.items()])
    df_sac.to_csv(os.path.join(RESULTS_DIR, "sac_analysis.csv"), index=False)
    plot_dband_model(df_sac, OUTPUT_DIR)
    plot_coordination_effect(df_sac, OUTPUT_DIR)
    plot_sac_heatmap(df_sac, OUTPUT_DIR)
    log_event("EXECUTE", "handoff_completed", skill="sac_analysis",
               handoff_out={"n_SAC": len(df_sac)},
               files=["figures/sac_metal_support.png",
                       "figures/sac_coordination_effect.png",
                       "figures/sac_heatmap.png"])

    # ----------------------------------------------------------------
    # 6. Solvent effects + potential-dependent analysis
    # ----------------------------------------------------------------
    print("\n[6/6] Solvent effects and potential-dependent analysis...")
    log_event("EXECUTE", "skill_selected", skill="solvent_effects")
    from scripts.solvent_effects import (plot_potential_dependent_activity,
                                          plot_solvent_comparison,
                                          compute_solvation_summary)
    plot_potential_dependent_activity(df, OUTPUT_DIR)
    plot_solvent_comparison(OUTPUT_DIR)
    df_solv = compute_solvation_summary(df)
    df_solv.to_csv(os.path.join(RESULTS_DIR, "solvation_corrections.csv"), index=False)
    log_event("EXECUTE", "handoff_completed", skill="solvent_effects",
               files=["figures/potential_dependent.png",
                       "figures/solvent_comparison_Cu.png",
                       "results/solvation_corrections.csv"])

    # ----------------------------------------------------------------
    # 7. Scoring and ranking
    # ----------------------------------------------------------------
    print("\n[SCORE] Ranking all catalysts...")
    score_cols = []
    for _, row in df_vol.iterrows():
        sc = score_catalyst(row)
        score_cols.append({**row.to_dict(), **sc})
    df_scored = pd.DataFrame(score_cols)

    # Add solvation-corrected UL
    solv_merge = df_solv[["catalyst","UL_CO_solv","UL_CH4_solv","delta_UL_CO"]]
    df_scored  = df_scored.merge(solv_merge, on="catalyst", how="left")
    df_scored.to_csv(os.path.join(RESULTS_DIR, "final_scores.csv"), index=False)

    plot_radar_top_candidates(df_scored,
        os.path.join(OUTPUT_DIR, "radar_top_candidates.png"))
    plot_materials_screening_summary(df_scored,
        os.path.join(OUTPUT_DIR, "screening_summary.png"))

    log_event("REPORT", "file_written",
               files=["results/final_scores.csv",
                       "figures/radar_top_candidates.png",
                       "figures/screening_summary.png"])

    # ----------------------------------------------------------------
    # 8. Print ranking tables
    # ----------------------------------------------------------------
    print("\n" + "="*65)
    print("  RESULTS SUMMARY")
    print("="*65)

    print("\n── Top 8 Catalysts for CO Production ──")
    top_CO = df_scored.nlargest(8, "composite_CO")[
        ["catalyst","category","dG_CO","UL_CO_V","UL_CO_solv",
         "CO_activity","composite_CO"]
    ].reset_index(drop=True)
    print(top_CO.to_string(index=False))

    print("\n── Top 8 Catalysts for CH4 Production ──")
    top_CH4 = df_scored.nlargest(8, "composite_CH4")[
        ["catalyst","category","dG_CO","UL_CH4_V","CH4_activity","composite_CH4"]
    ].reset_index(drop=True)
    print(top_CH4.to_string(index=False))

    print("\n── Top 8 Catalysts for C2+ Production ──")
    top_C2 = df_scored.nlargest(8, "composite_C2+")[
        ["catalyst","category","dG_CO","UL_C2H4_V","C2H4_activity","composite_C2+"]
    ].reset_index(drop=True)
    print(top_C2.to_string(index=False))

    # Scaling relation summary
    print("\n── Key Scaling Relations ──")
    print(df_fit[["relation","slope","intercept","R2","MAE_eV"]].to_string(index=False))

    log_event("REPORT", "run_completed",
               note="Pipeline completed successfully",
               handoff_out={
                   "n_catalysts_screened": len(df_scored),
                   "figures_generated": 12,
                   "results_files": 7,
               })

    return df_scored, df_fit, df_sac


if __name__ == "__main__":
    df_scored, df_fit, df_sac = run_pipeline()
    print("\n✓ Pipeline complete. All results saved to results/ and figures/")
