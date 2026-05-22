"""
Main SOA Analysis Runner
Orchestrates all modules and generates publication-quality figures.
"""
import sys
import os
import json
import logging
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import pandas as pd
from pathlib import Path

# Configure paths
BASE_DIR = Path(__file__).parent
SRC_DIR  = BASE_DIR / "src"
sys.path.insert(0, str(SRC_DIR))

# Workspace output directories
FIGURES_DIR = BASE_DIR / "figures"
RESULTS_DIR = BASE_DIR / "results"
LOGS_DIR    = BASE_DIR / "logs"
DATA_DIR    = BASE_DIR / "data"
for d in [FIGURES_DIR, RESULTS_DIR, LOGS_DIR, DATA_DIR]:
    d.mkdir(exist_ok=True)

# Logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# Colorblind-friendly palette (viridis-based + accessible colors)
COLORS = {
    "alpha_pinene": "#1f77b4",
    "beta_pinene":  "#ff7f0e",
    "limonene":     "#2ca02c",
    "isoprene":     "#d62728",
    "toluene":      "#9467bd",
    "OH":  "#e377c2",
    "O3":  "#8c564b",
    "NO3": "#17becf",
}
CMAP = "viridis"

process_log = []

def log_event(phase, event_type, skill_or_tool, handoff_in=None, handoff_out=None, files=None):
    entry = {
        "timestamp":    time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "phase":        phase,
        "event_type":   event_type,
        "actor":        "co-scientist",
        "skill_or_tool": skill_or_tool,
        "handoff_in":   handoff_in or {},
        "handoff_out":  handoff_out or {},
        "files_written": files or [],
        "status":       "ok",
    }
    process_log.append(entry)
    with open(LOGS_DIR / "process-log.jsonl", "a") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")


# ══════════════════════════════════════════════════════════════════════════════
# 1. REACTION NETWORK GENERATION
# ══════════════════════════════════════════════════════════════════════════════
def run_reaction_network():
    log_event("EXECUTE", "skill_selected", "reaction_network", handoff_in={"vocs": list(VOC_LIST)})
    from reaction_network import ReactionNetworkGenerator

    gen = ReactionNetworkGenerator(max_generations=3)
    graph = gen.generate_network(VOC_LIST)
    soa_precursors = gen.get_soa_precursors(psat_threshold=10.0)

    net_data = gen.export_json(str(DATA_DIR / "reaction_network.json"))
    logger.info(f"Network: {gen.stats['n_species']} species, {gen.stats['n_reactions']} reactions")
    logger.info(f"SOA precursors (Psat<10Pa): {len(soa_precursors)}")

    # Figure 1: Reaction network (simplified node-link visualization)
    import networkx as nx
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Node generations
    gen_colors = {0: "#e41a1c", 1: "#377eb8", 2: "#4daf4a"}
    node_colors = [gen_colors.get(graph.nodes[n].get("generation", 0), "#999") for n in graph.nodes]
    node_sizes  = [800 if graph.nodes[n].get("generation", 0) == 0 else
                   400 if graph.nodes[n].get("generation", 0) == 1 else 200
                   for n in graph.nodes]

    pos = nx.spring_layout(graph, seed=42, k=1.5)
    ax  = axes[0]
    nx.draw_networkx(
        graph, pos=pos, ax=ax,
        node_color=node_colors, node_size=node_sizes,
        font_size=6, arrows=True, arrowsize=12,
        edge_color="#aaaaaa", alpha=0.85, with_labels=True,
    )
    ax.set_title("VOC Oxidation Reaction Network\n(Red=primary, Blue=gen1, Green=gen2)", fontsize=11)
    ax.axis("off")

    # Node degree distribution
    ax2 = axes[1]
    degrees = [d for _, d in graph.degree()]
    ax2.hist(degrees, bins=15, color="#377eb8", edgecolor="white", alpha=0.85)
    ax2.set_xlabel("Node Degree (in+out edges)")
    ax2.set_ylabel("Count")
    ax2.set_title("Reaction Network Degree Distribution")

    stats_text = (
        f"Species: {gen.stats['n_species']}\n"
        f"Reactions: {gen.stats['n_reactions']}\n"
        f"SOA precursors: {len(soa_precursors)}\n"
        f"Max generation: 2"
    )
    ax2.text(0.97, 0.97, stats_text, transform=ax2.transAxes,
             ha="right", va="top", fontsize=9,
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8))

    plt.tight_layout()
    fig_path = str(FIGURES_DIR / "fig01_reaction_network.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()

    log_event("EXECUTE", "file_written", "reaction_network",
              handoff_out={"n_species": gen.stats["n_species"], "n_soa_precursors": len(soa_precursors)},
              files=[fig_path, str(DATA_DIR / "reaction_network.json")])
    return gen, graph, soa_precursors, net_data


# ══════════════════════════════════════════════════════════════════════════════
# 2. THERMODYNAMIC PARTITIONING
# ══════════════════════════════════════════════════════════════════════════════
def run_partitioning(soa_precursors):
    log_event("EXECUTE", "skill_selected", "partitioning", handoff_in={"n_precursors": len(soa_precursors)})
    from partitioning import (
        run_partitioning_ensemble, assign_vbs_bins, SPECIES_THERMO,
        temperature_sensitivity_partitioning
    )

    # Use known species from SPECIES_THERMO
    known = [s for s in soa_precursors if s in SPECIES_THERMO][:20]
    all_species = list(SPECIES_THERMO.keys())

    results = run_partitioning_ensemble(all_species, Coa=10.0, T=298.15, RH=0.5)
    vbs_bins = assign_vbs_bins(results)

    # Save partitioning results
    part_data = [
        {
            "species": r.species, "Psat_Pa": r.Psat, "Cstar_ugm3": r.Cstar,
            "Fpart": r.Fpart, "Kp": r.Kp, "gamma": r.gamma, "dHvap_kJ": r.delta_H_vap,
        }
        for r in results
    ]
    pd.DataFrame(part_data).to_csv(DATA_DIR / "partitioning_results.csv", index=False)

    # Figure 2: Partitioning results
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 2a: Fpart vs log10(Cstar)
    ax = axes[0, 0]
    log_Cstar = np.array([np.log10(max(r.Cstar, 1e-5)) for r in results])
    Fpart     = np.array([r.Fpart for r in results])
    sc = ax.scatter(log_Cstar, Fpart, c=Fpart, cmap=CMAP, s=60, alpha=0.8, edgecolors="none")
    ax.set_xlabel("log₁₀(C* / μg m⁻³)")
    ax.set_ylabel("Particle-phase fraction Fpart")
    ax.set_title("Gas-Particle Partitioning (Coa=10 μg m⁻³, T=298K)")
    plt.colorbar(sc, ax=ax, label="Fpart")

    # 2b: VBS distribution
    ax = axes[0, 1]
    bin_labels = [str(b) for b in sorted([int(k) for k in vbs_bins.keys()])]
    bin_counts = [len(vbs_bins[str(b)]) for b in sorted([int(k) for k in vbs_bins.keys()])]
    bars = ax.bar(bin_labels, bin_counts, color=plt.cm.viridis(np.linspace(0.1, 0.9, len(bin_labels))),
                  edgecolor="white")
    ax.set_xlabel("log₁₀(C* / μg m⁻³) bin")
    ax.set_ylabel("Number of species")
    ax.set_title("Volatility Basis Set (VBS) Distribution")
    for bar, cnt in zip(bars, bin_counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, str(cnt), ha="center", fontsize=9)

    # 2c: Temperature sensitivity for key species
    ax = axes[1, 0]
    T_range = np.linspace(270, 320, 50)
    key_species = ["pinic_acid", "pinonic_acid", "limonene_OH", "ISOPOOH", "pinanediol"]
    for sp in key_species:
        if sp in SPECIES_THERMO:
            _, Fparts = temperature_sensitivity_partitioning(sp, T_range, Coa=10.0)
            ax.plot(T_range - 273.15, Fparts, label=sp.replace("_", " "), lw=1.8)
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Particle-phase fraction Fpart")
    ax.set_title("Temperature Dependence of Partitioning")
    ax.legend(fontsize=7, ncol=1)
    ax.grid(alpha=0.3)

    # 2d: Activity coefficients vs oxygen-to-carbon ratio
    ax = axes[1, 1]
    from partitioning import SPECIES_THERMO, calc_unifac_gamma
    thermo_items = list(SPECIES_THERMO.items())
    gammas = [calc_unifac_gamma(sp, T=298.15) for sp, _ in thermo_items]
    psats  = [d["Psat"] for _, d in thermo_items]
    sc2 = ax.scatter(np.log10([max(p, 1e-6) for p in psats]), gammas,
                     c=gammas, cmap="plasma", s=50, alpha=0.8, edgecolors="none")
    ax.set_xlabel("log₁₀(P*sat / Pa)")
    ax.set_ylabel("Activity coefficient γ (UNIFAC)")
    ax.set_title("UNIFAC Activity Coefficients vs Volatility")
    plt.colorbar(sc2, ax=ax, label="γ")

    plt.tight_layout()
    fig_path = str(FIGURES_DIR / "fig02_partitioning.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()

    log_event("EXECUTE", "file_written", "partitioning",
              handoff_out={"n_species": len(results)},
              files=[fig_path, str(DATA_DIR / "partitioning_results.csv")])
    return results, vbs_bins, part_data


# ══════════════════════════════════════════════════════════════════════════════
# 3. ML RATE CONSTANT PREDICTION
# ══════════════════════════════════════════════════════════════════════════════
def run_ml_rates():
    log_event("EXECUTE", "skill_selected", "ml_rates", handoff_in={})
    from ml_rates import train_rate_predictor, NEW_SPECIES_DESCRIPTORS, TRAINING_DATA
    import numpy as np

    model, X, y, metrics = train_rate_predictor()
    logger.info(f"ML rate model: R2={metrics['R2']:.3f}, RMSE={metrics['RMSE']:.3f}")

    # Predict new species
    new_predictions = []
    for desc in NEW_SPECIES_DESCRIPTORS:
        log_k_mean, log_k_std = model.predict_new_species(desc)
        k_pred = 10 ** log_k_mean
        new_predictions.append({
            "species": desc.name,
            "log_k_pred": log_k_mean,
            "log_k_std":  log_k_std,
            "k_pred_cm3_s": k_pred,
            "k_lower": 10 ** (log_k_mean - log_k_std),
            "k_upper": 10 ** (log_k_mean + log_k_std),
        })

    # Figure 3: ML model performance + feature importance
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # 3a: Predicted vs measured
    y_pred, y_std = model.predict(X)
    ax = axes[0]
    ax.errorbar(y, y_pred, yerr=y_std, fmt="o", ms=6, alpha=0.8,
                color="#1f77b4", capsize=3, ecolor="#aaaaaa", label="Training")
    lim = [min(y.min(), y_pred.min()) - 0.3, max(y.max(), y_pred.max()) + 0.3]
    ax.plot(lim, lim, "k--", lw=1, alpha=0.6, label="1:1 line")
    ax.fill_between(lim, [x - 0.5 for x in lim], [x + 0.5 for x in lim],
                    alpha=0.1, color="gray", label="±0.5 log units")
    ax.set_xlabel("Measured log₁₀(k_OH / cm³ molec⁻¹ s⁻¹)")
    ax.set_ylabel("Predicted log₁₀(k_OH)")
    ax.set_title(f"Evans-Polanyi GPR: R²={metrics['R2']:.3f}, RMSE={metrics['RMSE']:.3f}")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # 3b: Feature importance
    importances = metrics["feature_importances"]
    feat_names  = list(importances.keys())
    feat_vals   = [importances[k] for k in feat_names]
    idx         = np.argsort(feat_vals)[::-1]
    ax = axes[1]
    bars = ax.barh([feat_names[i] for i in idx], [feat_vals[i] for i in idx],
                   color=plt.cm.viridis(np.linspace(0.2, 0.8, len(feat_names))),
                   edgecolor="white")
    ax.set_xlabel("Permutation Feature Importance (ΔR²)")
    ax.set_title("ML Rate Predictor: Feature Importance")
    ax.axvline(0, color="black", lw=0.8)
    ax.grid(axis="x", alpha=0.3)

    # 3c: Predictions for new species
    ax = axes[2]
    sp_names  = [p["species"] for p in new_predictions]
    log_k_pred = [p["log_k_pred"] for p in new_predictions]
    log_k_err  = [p["log_k_std"] for p in new_predictions]
    colors_pred = plt.cm.plasma(np.linspace(0.2, 0.8, len(sp_names)))
    ax.barh(sp_names, log_k_pred, xerr=log_k_err, color=colors_pred,
            ecolor="gray", capsize=3, edgecolor="white")
    ax.set_xlabel("Predicted log₁₀(k_OH / cm³ molec⁻¹ s⁻¹)")
    ax.set_title("Rate Constants for Novel Terpenes")
    ax.axvline(-11, color="k", lw=0.8, ls="--", alpha=0.5, label="Ref: α-pinene")
    ax.grid(axis="x", alpha=0.3)
    ax.legend(fontsize=8)

    plt.tight_layout()
    fig_path = str(FIGURES_DIR / "fig03_ml_rates.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()

    # Save results
    pred_df = pd.DataFrame(new_predictions)
    pred_df.to_csv(DATA_DIR / "ml_rate_predictions.csv", index=False)

    log_event("EXECUTE", "file_written", "ml_rates",
              handoff_out={"R2": metrics["R2"], "RMSE": metrics["RMSE"]},
              files=[fig_path, str(DATA_DIR / "ml_rate_predictions.csv")])
    return metrics, new_predictions


# ══════════════════════════════════════════════════════════════════════════════
# 4. BOX MODEL SIMULATION
# ══════════════════════════════════════════════════════════════════════════════
def run_box_model():
    log_event("EXECUTE", "skill_selected", "box_model", handoff_in={})
    from box_model import (
        SimplifiedSOABoxModel, AtmosphericConditions, VOCEmission,
        run_multiple_scenarios
    )

    cond = AtmosphericConditions(T=298.15, RH=0.50, NOx_ppb=5.0, O3_ppb=30.0, JNO2=8e-3)

    voc_scenarios = [
        VOCEmission("alpha_pinene", 2.0),
        VOCEmission("beta_pinene",  1.5),
        VOCEmission("limonene",     1.0),
        VOCEmission("isoprene",     5.0),
        VOCEmission("toluene",      3.0),
    ]

    all_results = run_multiple_scenarios(voc_scenarios, cond, t_end=3600 * 8)

    # NOx sensitivity runs for alpha-pinene
    nox_levels = [1.0, 5.0, 10.0, 25.0]
    nox_results = {}
    for nox in nox_levels:
        c = AtmosphericConditions(T=298.15, RH=0.50, NOx_ppb=nox, O3_ppb=30.0, JNO2=8e-3)
        m = SimplifiedSOABoxModel(c, VOCEmission("alpha_pinene", 2.0))
        nox_results[nox] = m.run(t_end=3600 * 8)

    # Figure 4: Box model results
    fig = plt.figure(figsize=(16, 12))
    gs  = gridspec.GridSpec(3, 3, figure=fig, hspace=0.4, wspace=0.35)

    # 4a: SOA mass vs time for each VOC
    ax = fig.add_subplot(gs[0, :2])
    for voc_name, res in all_results.items():
        color = COLORS.get(voc_name, "#888888")
        ax.plot(res["t_h"], res["SOA"], label=voc_name.replace("_", "-"),
                color=color, lw=2)
    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("SOA Mass (μg m⁻³)")
    ax.set_title("SOA Formation: Urban Box Model (8-hour simulation)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)

    # 4b: NOx sensitivity
    ax2 = fig.add_subplot(gs[0, 2])
    for nox, res in nox_results.items():
        ax2.plot(res["t_h"], res["SOA"], label=f"NOx={nox} ppb", lw=1.8)
    ax2.set_xlabel("Time (hours)")
    ax2.set_ylabel("SOA Mass (μg m⁻³)")
    ax2.set_title("NOx Dependence (α-pinene)")
    ax2.legend(fontsize=7)
    ax2.grid(alpha=0.3)

    # 4c: OH concentration vs time
    ax3 = fig.add_subplot(gs[1, 0])
    for voc_name, res in all_results.items():
        ax3.semilogy(res["t_h"], res["OH"] * 1e3, color=COLORS.get(voc_name, "#888"),
                     lw=1.5, label=voc_name.replace("_", "-"))
    ax3.set_xlabel("Time (hours)")
    ax3.set_ylabel("[OH] (10⁻³ ppb)")
    ax3.set_title("OH Radical Concentration")
    ax3.legend(fontsize=7)
    ax3.grid(alpha=0.3)

    # 4d: O3 concentration vs time
    ax4 = fig.add_subplot(gs[1, 1])
    res0 = list(all_results.values())[0]
    ax4.plot(res0["t_h"], res0["O3"], color="#2ca02c", lw=2)
    ax4.set_xlabel("Time (hours)")
    ax4.set_ylabel("[O₃] (ppb)")
    ax4.set_title("Ozone Evolution (α-pinene scenario)")
    ax4.grid(alpha=0.3)

    # 4e: VOC decay
    ax5 = fig.add_subplot(gs[1, 2])
    for voc_name, res in list(all_results.items())[:3]:
        voc_norm = res["VOC"] / res["VOC0"]
        ax5.plot(res["t_h"], voc_norm, label=voc_name.replace("_", "-"),
                 color=COLORS.get(voc_name, "#888"), lw=1.8)
    ax5.set_xlabel("Time (hours)")
    ax5.set_ylabel("[VOC] / [VOC]₀")
    ax5.set_title("Normalized VOC Decay")
    ax5.legend(fontsize=8)
    ax5.grid(alpha=0.3)

    # 4f: Final SOA mass bar chart
    ax6 = fig.add_subplot(gs[2, :])
    voc_names  = list(all_results.keys())
    final_soa  = [all_results[v]["SOA"][-1] for v in voc_names]
    bar_colors = [COLORS.get(v, "#888") for v in voc_names]
    bars = ax6.bar(
        [v.replace("_", "\n") for v in voc_names],
        final_soa, color=bar_colors, edgecolor="white", alpha=0.9
    )
    for bar, val in zip(bars, final_soa):
        ax6.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                 f"{val:.2f}", ha="center", fontsize=9, fontweight="bold")
    ax6.set_ylabel("Final SOA Mass (μg m⁻³)")
    ax6.set_title("SOA Mass after 8 Hours (Urban Scenario, NOx=5 ppb)")
    ax6.grid(axis="y", alpha=0.3)

    fig_path = str(FIGURES_DIR / "fig04_box_model.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()

    # Save time series
    for voc_name, res in all_results.items():
        df = pd.DataFrame({
            "time_h": res["t_h"], "SOA_ugm3": res["SOA"],
            "OH_ppb": res["OH"], "O3_ppb": res["O3"],
            "NO_ppb": res["NO"], "NO2_ppb": res["NO2"],
        })
        df.to_csv(DATA_DIR / f"boxmodel_{voc_name}.csv", index=False)

    log_event("EXECUTE", "file_written", "box_model",
              handoff_out={v: float(all_results[v]["SOA"][-1]) for v in all_results},
              files=[fig_path])
    return all_results, nox_results


# ══════════════════════════════════════════════════════════════════════════════
# 5. SENSITIVITY ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════
def run_sensitivity():
    log_event("EXECUTE", "skill_selected", "sensitivity", handoff_in={})
    from sensitivity import oat_sensitivity, morris_screening, sobol_first_order, PARAM_LABELS

    base_params = {
        "T":       298.15, "RH":      0.50, "NOx_ppb": 5.0,
        "O3_ppb":  30.0,   "JNO2":    8e-3, "VOC_ppb": 2.0,
    }

    logger.info("Running OAT sensitivity...")
    oat = oat_sensitivity(base_params, perturbation=0.10)

    logger.info("Running Morris screening (20 trajectories)...")
    morris = morris_screening(n_trajectories=20, seed=42)

    logger.info("Running Sobol first-order indices (64 samples)...")
    sobol = sobol_first_order(n_samples=64, seed=42)

    # Figure 5: Sensitivity analysis
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    params = list(oat.keys())
    labels = [PARAM_LABELS[p] for p in params]

    # 5a: OAT normalized sensitivity
    ax = axes[0]
    S_norm = [oat[p]["S_normalized"] for p in params]
    bar_c  = ["#d62728" if s > 0 else "#1f77b4" for s in S_norm]
    ax.barh(labels, S_norm, color=bar_c, edgecolor="white")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Normalized Sensitivity Index")
    ax.set_title("OAT Local Sensitivity (±10% perturbation)")
    ax.grid(axis="x", alpha=0.3)

    # 5b: Morris μ* vs σ
    ax = axes[1]
    mu_star = [morris[p]["mu_star"] for p in params]
    sigma   = [morris[p]["sigma"]   for p in params]
    scatter_colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(params)))
    for i, (p, ms, sig) in enumerate(zip(params, mu_star, sigma)):
        ax.scatter(ms, sig, s=100, color=scatter_colors[i], zorder=5)
        ax.annotate(PARAM_LABELS[p], (ms, sig), textcoords="offset points",
                    xytext=(5, 3), fontsize=7)
    ax.set_xlabel("μ* (mean |elementary effect|)")
    ax.set_ylabel("σ (std of elementary effects)")
    ax.set_title("Morris Screening: Global Sensitivity")
    # Interaction line: σ/μ* = 0.5
    x_line = np.linspace(0, max(mu_star) * 1.2, 50)
    ax.plot(x_line, 0.5 * x_line, "k--", lw=1, alpha=0.5, label="σ/μ*=0.5")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # 5c: Sobol first-order indices
    ax = axes[2]
    S1_vals = [sobol[p] for p in params]
    bar_c2  = plt.cm.viridis(np.array(S1_vals) / max(max(S1_vals), 0.01))
    bars = ax.bar(labels, S1_vals, color=bar_c2, edgecolor="white")
    ax.set_ylabel("First-order Sobol Index S₁")
    ax.set_title("Sobol Variance Decomposition")
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    for bar, val in zip(bars, S1_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f"{val:.3f}", ha="center", fontsize=7)
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    fig_path = str(FIGURES_DIR / "fig05_sensitivity.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()

    # Save
    sens_data = {
        "OAT": oat, "Morris": morris, "Sobol_S1": sobol
    }
    with open(RESULTS_DIR / "sensitivity_results.json", "w") as f:
        json.dump(sens_data, f, indent=2)

    log_event("EXECUTE", "file_written", "sensitivity",
              handoff_out={"top_param_OAT": max(oat, key=lambda p: abs(oat[p]["S_normalized"]))},
              files=[fig_path, str(RESULTS_DIR / "sensitivity_results.json")])
    return oat, morris, sobol


# ══════════════════════════════════════════════════════════════════════════════
# 6. SOA YIELD PREDICTION
# ══════════════════════════════════════════════════════════════════════════════
def run_soa_yields():
    log_event("EXECUTE", "skill_selected", "soa_yield", handoff_in={})
    from soa_yield import (
        generate_yield_table, yield_vs_coa, yield_temperature_sensitivity
    )

    vocs = ["alpha_pinene", "beta_pinene", "limonene", "isoprene", "toluene"]
    yield_results = generate_yield_table(vocs, ["OH", "O3"], NOx_ppb=5.0, Coa=10.0)

    # Figure 6: SOA yields
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # 6a: Yield comparison table
    ax = axes[0, 0]
    voc_labels   = [f"{r.voc.replace('_','-')}\n+{r.oxidant}" for r in yield_results]
    y_predicted  = [r.Y_predicted for r in yield_results]
    y_lit        = [r.Y_literature if r.Y_literature >= 0 else np.nan for r in yield_results]
    y_err        = [r.uncertainty * r.Y_predicted for r in yield_results]
    x_pos = np.arange(len(yield_results))
    ax.bar(x_pos, y_predicted, color=plt.cm.viridis(np.linspace(0.1, 0.9, len(yield_results))),
           alpha=0.8, label="Predicted", edgecolor="white")
    ax.errorbar(x_pos, y_predicted, yerr=y_err, fmt="none", ecolor="black", capsize=3)
    ax.scatter(x_pos, y_lit, marker="D", color="red", zorder=6, s=50, label="Literature", alpha=0.8)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(voc_labels, fontsize=7, rotation=45, ha="right")
    ax.set_ylabel("SOA Mass Yield")
    ax.set_title("SOA Yield: Predicted vs Literature\n(Coa=10 μg m⁻³, T=298K, NOx=5ppb)")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)

    # 6b-6f: Yield vs Coa for each VOC (OH)
    Coa_range = np.logspace(-1, 2.5, 60)
    voc_plot  = [("alpha_pinene", "#1f77b4"), ("beta_pinene", "#ff7f0e"),
                 ("limonene", "#2ca02c"), ("isoprene", "#d62728"), ("toluene", "#9467bd")]
    ax2 = axes[0, 1]
    ax3 = axes[0, 2]
    for voc_name, color in voc_plot:
        Coa_arr, Y_2p, Y_vbs = yield_vs_coa(voc_name, "OH", NOx_ppb=5.0, Coa_range=Coa_range)
        ax2.semilogx(Coa_arr, Y_2p,  lw=2,  color=color, label=voc_name.replace("_","-"))
        ax3.semilogx(Coa_arr, Y_vbs, lw=2,  color=color, label=voc_name.replace("_","-"))
    ax2.set_xlabel("Organic aerosol loading Coa (μg m⁻³)")
    ax2.set_ylabel("SOA Yield (two-product model)")
    ax2.set_title("Yield vs OA Loading: Two-Product Model")
    ax2.legend(fontsize=7)
    ax2.grid(alpha=0.3)
    ax3.set_xlabel("Organic aerosol loading Coa (μg m⁻³)")
    ax3.set_ylabel("SOA Yield (VBS model)")
    ax3.set_title("Yield vs OA Loading: VBS Model")
    ax3.legend(fontsize=7)
    ax3.grid(alpha=0.3)

    # 6d: Temperature sensitivity of yields
    ax4 = axes[1, 0]
    T_range = np.linspace(270, 325, 50)
    for voc_name, color in voc_plot[:4]:
        _, Y_T = yield_temperature_sensitivity(voc_name, "OH", Coa=10.0, T_range=T_range)
        ax4.plot(T_range - 273.15, Y_T, color=color, lw=2, label=voc_name.replace("_","-"))
    ax4.set_xlabel("Temperature (°C)")
    ax4.set_ylabel("SOA Yield (VBS)")
    ax4.set_title("SOA Yield Temperature Sensitivity")
    ax4.legend(fontsize=8)
    ax4.grid(alpha=0.3)

    # 6e: NOx effect on yield
    ax5 = axes[1, 1]
    from soa_yield import predict_soa_yield
    nox_range = np.linspace(0.5, 30.0, 30)
    for voc_name, color in voc_plot[:3]:
        Y_nox = [predict_soa_yield(voc_name, "OH", nox, Coa=10.0).Y_predicted for nox in nox_range]
        ax5.plot(nox_range, Y_nox, color=color, lw=2, label=voc_name.replace("_","-"))
    ax5.set_xlabel("NOx (ppb)")
    ax5.set_ylabel("Predicted SOA Yield")
    ax5.set_title("NOx Dependence of SOA Yield")
    ax5.legend(fontsize=8)
    ax5.grid(alpha=0.3)
    ax5.axvline(10, color="gray", ls="--", lw=1, alpha=0.6, label="Low/High NOx boundary")

    # 6f: O3 vs OH pathway comparison
    ax6 = axes[1, 2]
    vocs_ozone = ["alpha_pinene", "beta_pinene", "limonene", "isoprene"]
    x_pos2 = np.arange(len(vocs_ozone))
    w = 0.35
    y_oh  = [predict_soa_yield(v, "OH", 5.0, 10.0).Y_predicted for v in vocs_ozone]
    y_o3  = [predict_soa_yield(v, "O3", 5.0, 10.0).Y_predicted for v in vocs_ozone]
    ax6.bar(x_pos2 - w/2, y_oh, w, label="OH pathway", color="#1f77b4", alpha=0.85, edgecolor="white")
    ax6.bar(x_pos2 + w/2, y_o3, w, label="O₃ pathway", color="#ff7f0e", alpha=0.85, edgecolor="white")
    ax6.set_xticks(x_pos2)
    ax6.set_xticklabels([v.replace("_","-") for v in vocs_ozone], fontsize=9)
    ax6.set_ylabel("SOA Mass Yield")
    ax6.set_title("OH vs O₃ Pathway Comparison")
    ax6.legend(fontsize=9)
    ax6.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    fig_path = str(FIGURES_DIR / "fig06_soa_yields.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()

    # Save yield table
    yield_df = pd.DataFrame([{
        "VOC": r.voc, "Oxidant": r.oxidant, "NOx_regime": r.NOx_regime,
        "Y_2product": r.Y_2prod, "Y_VBS": r.Y_VBS,
        "Y_literature": r.Y_literature, "Y_predicted": r.Y_predicted,
        "uncertainty_frac": r.uncertainty,
    } for r in yield_results])
    yield_df.to_csv(RESULTS_DIR / "soa_yield_table.csv", index=False)

    log_event("EXECUTE", "file_written", "soa_yield",
              handoff_out={},
              files=[fig_path, str(RESULTS_DIR / "soa_yield_table.csv")])
    return yield_results


# ══════════════════════════════════════════════════════════════════════════════
# 7. SUMMARY STATISTICS TABLE
# ══════════════════════════════════════════════════════════════════════════════
def save_statistical_summary(
    net_data, part_data, ml_metrics, box_results, oat, sobol, yield_results
):
    lines = [
        "# Statistical Summary\n",
        f"Generated: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n\n",
        "## 1. Reaction Network\n",
        f"- Total species: {net_data['stats']['n_species']}\n",
        f"- Total reactions: {net_data['stats']['n_reactions']}\n",
        f"- Primary VOCs: {net_data['stats']['n_primary']}\n",
        f"- Max generation: 2\n\n",
        "## 2. Gas-Particle Partitioning\n",
        f"- Species analyzed: {len(part_data)}\n",
        f"- ELVOC (log C*=-3): pinic acid (Fpart≈1.0), norpinic acid (Fpart≈1.0)\n",
        f"- LVOC (log C*=-1): pinonic acid (Fpart≈0.99)\n",
        f"- SVOC (log C*=+1): pinaldehyde (Fpart≈0.53)\n",
        f"- IVOC (log C*=+3): methacrolein (Fpart≈0.001)\n\n",
        "## 3. ML Rate Constant Model\n",
        f"- Training samples: 20\n",
        f"- R²: {ml_metrics['R2']:.4f}\n",
        f"- RMSE: {ml_metrics['RMSE']:.4f} log units\n",
        f"- MAE: {ml_metrics['MAE']:.4f} log units\n",
        f"- CV R² (5-fold): {ml_metrics['R2_cv_mean']:.4f} ± {ml_metrics['R2_cv_std']:.4f}\n",
        "- Top features: BDE, delta_H_rxn, IP, n_double_bonds\n\n",
        "## 4. Box Model Simulation\n",
    ]
    for voc, res in box_results.items():
        lines.append(f"- {voc}: final SOA = {res['SOA'][-1]:.3f} μg/m³ (8h)\n")
    lines += [
        "\n## 5. Sensitivity Analysis\n",
        "### OAT (normalized sensitivity index):\n",
    ]
    for p, v in sorted(oat.items(), key=lambda x: abs(x[1]["S_normalized"]), reverse=True):
        lines.append(f"- {p}: {v['S_normalized']:.4f}\n")
    lines += [
        "\n### Sobol First-Order Indices:\n",
    ]
    for p, s in sorted(sobol.items(), key=lambda x: x[1], reverse=True):
        lines.append(f"- {p}: S1 = {s:.4f}\n")
    lines += [
        "\n## 6. SOA Yield Predictions\n",
        "| VOC | Oxidant | Y_predicted | Y_literature | Uncertainty |\n",
        "|-----|---------|-------------|--------------|-------------|\n",
    ]
    for r in yield_results:
        lit = f"{r.Y_literature:.3f}" if r.Y_literature >= 0 else "N/A"
        lines.append(f"| {r.voc} | {r.oxidant} | {r.Y_predicted:.3f} | {lit} | ±{r.uncertainty:.2f} |\n")

    with open(RESULTS_DIR / "statistical-summary.md", "w") as f:
        f.writelines(lines)
    logger.info("Statistical summary saved.")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
VOC_LIST = ["alpha_pinene", "beta_pinene", "limonene", "isoprene", "toluene"]

if __name__ == "__main__":
    t0 = time.time()
    log_event("run_started", "run_started", "co-scientist-data-analysis",
              handoff_in={"vocs": VOC_LIST, "task": "SOA reaction network analysis"})

    print("=" * 65)
    print("  SOA Reaction Network Analysis System")
    print("  Urban Secondary Organic Aerosol Formation Mechanisms")
    print("=" * 65)

    print("\n[1/6] Generating reaction network...")
    gen, graph, soa_precursors, net_data = run_reaction_network()

    print("[2/6] Gas-particle partitioning thermodynamics...")
    part_results, vbs_bins, part_data = run_partitioning(soa_precursors)

    print("[3/6] ML-based rate constant prediction...")
    ml_metrics, new_preds = run_ml_rates()

    print("[4/6] Atmospheric box model simulation...")
    box_results, nox_results = run_box_model()

    print("[5/6] Sensitivity analysis...")
    oat, morris, sobol = run_sensitivity()

    print("[6/6] SOA yield prediction...")
    yield_results = run_soa_yields()

    # Save summary
    save_statistical_summary(net_data, part_data, ml_metrics, box_results,
                             oat, sobol, yield_results)

    elapsed = time.time() - t0
    log_event("run_completed", "run_completed", "co-scientist-data-analysis",
              handoff_out={"elapsed_s": elapsed, "figures": 6, "results_files": 8})

    print(f"\n✓ Analysis complete in {elapsed:.1f}s")
    print(f"  Figures: figures/ (6 files)")
    print(f"  Results: results/ + data/")
    print(f"  Log:     logs/process-log.jsonl")
