"""
CrMnFeCoNi HEA Case Study — Ultra-High-Temperature Design
Full pipeline: descriptor calc → surrogate training → Bayesian optimization
→ active learning → Pareto front analysis → experimental recommendations.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
import json
import warnings
from datetime import datetime
warnings.filterwarnings("ignore")

from src.hea_descriptors import (
    compute_all_descriptors, descriptors_dataframe,
    gibbs_free_energy, calphad_phase_diagram_1d,
    ELEMENT_PROPS
)
from src.surrogate_models import HEAPropertySimulator, HEASurrogateModel, DESCRIPTOR_COLS, PROPERTY_COLS
from src.bayesian_opt import (
    MultiObjectiveBayesianOptimizer, ActiveLearningSelector,
    sample_compositions, latin_hypercube_compositions, ELEMENTS_5
)
from src.dft_and_databases import DFTSimulator, MaterialsDatabaseClient


# -----------------------------------------------------------------------
# Logging helper
# -----------------------------------------------------------------------
LOG_PATH = "logs/process-log.jsonl"
os.makedirs("logs", exist_ok=True)
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("data", exist_ok=True)


def log_event(phase: str, event_type: str, skill: str = "co-scientist-computational-materials",
              handoff_in: dict = None, handoff_out: dict = None,
              files_written: list = None, status: str = "ok"):
    entry = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill,
        "handoff_in": handoff_in or {},
        "handoff_out": handoff_out or {},
        "files_written": files_written or [],
        "status": status,
    }
    with open(LOG_PATH, "a") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")


# -----------------------------------------------------------------------
# STEP 1: CALPHAD-inspired thermodynamic landscape
# -----------------------------------------------------------------------
def step1_calphad(elements=ELEMENTS_5):
    print("\n" + "="*60)
    print("STEP 1: CALPHAD thermodynamic landscape")
    print("="*60)

    # Scan Gibbs free energy over temperature for equimolar
    equimolar = {el: 1/5 for el in elements}
    T_range = np.linspace(300, 1800, 50)
    G_vals = [gibbs_free_energy(equimolar, T) for T in T_range]

    # Binary phase diagram (Cr-Ni representative)
    T_scan = np.array([600, 800, 1000, 1200, 1400])
    x_scan = np.linspace(0.01, 0.99, 50)
    G_binary = calphad_phase_diagram_1d("Cr", "Ni", T_scan, x_scan)

    # Save data
    df_gibbs = pd.DataFrame({"T_K": T_range, "G_mix_kJmol": G_vals})
    df_gibbs.to_csv("data/calphad_equimolar_gibbs.csv", index=False)

    df_binary = pd.DataFrame(G_binary, index=[f"T={t}" for t in T_scan],
                             columns=[f"x={x:.2f}" for x in x_scan])
    df_binary.to_csv("data/calphad_CrNi_binary.csv")

    # Compute descriptors for equimolar
    desc = compute_all_descriptors(equimolar)
    print(f"  Equimolar CrMnFeCoNi descriptors:")
    for k, v in desc.items():
        if k not in ("phase",) and not k.startswith("x_"):
            print(f"    {k:15s} = {v:.4f}")
    print(f"    phase          = {desc['phase']}")

    log_event("CALPHAD", "file_written",
              files_written=["data/calphad_equimolar_gibbs.csv", "data/calphad_CrNi_binary.csv"],
              handoff_out={"G_min_kJ": min(G_vals), "T_stable": T_range[np.argmin(G_vals)]})
    return df_gibbs, G_binary, x_scan, T_scan


# -----------------------------------------------------------------------
# STEP 2: Generate training dataset
# -----------------------------------------------------------------------
def step2_generate_data(elements=ELEMENTS_5, n_train=200, n_candidates=1000, seed=42):
    print("\n" + "="*60)
    print("STEP 2: Dataset generation (LHS + DFT simulator)")
    print("="*60)

    # Extended composition space: allow Al, Ti substitution
    elements_ext = ["Cr", "Mn", "Fe", "Co", "Ni", "Al", "Ti"]

    # Training set: LHS sampling
    train_comps_5 = latin_hypercube_compositions(n_train // 2, elements, x_min=0.05)
    train_comps_7 = latin_hypercube_compositions(n_train // 2, elements_ext, x_min=0.02)
    all_comps = train_comps_5 + train_comps_7

    # Compute descriptors
    desc_df = descriptors_dataframe(all_comps, T=1000.0)

    # Simulate properties
    simulator = HEAPropertySimulator(noise_level=0.05, random_state=seed)
    df_train = simulator.simulate_dataset(desc_df)

    # DFT simulator
    dft = DFTSimulator(random_state=seed)
    df_dft = dft.generate_dataset(all_comps)
    df_dft.to_csv("data/dft_simulated.csv", index=False)

    # Candidate pool for active learning
    cand_comps = sample_compositions(n_candidates, elements, x_min=0.05, seed=seed+1)
    df_candidates = descriptors_dataframe(cand_comps, T=1000.0)
    df_candidates.to_csv("data/candidate_pool.csv", index=False)

    df_train.to_csv("data/training_dataset.csv", index=False)

    print(f"  Training samples: {len(df_train)}")
    print(f"  Candidate pool: {len(df_candidates)}")
    print(f"  Descriptor columns: {list(DESCRIPTOR_COLS)}")
    print(f"  Property stats:")
    for p in PROPERTY_COLS:
        print(f"    {p}: mean={df_train[p].mean():.2f}, std={df_train[p].std():.2f}")

    log_event("DataGeneration", "file_written",
              files_written=["data/training_dataset.csv", "data/dft_simulated.csv",
                             "data/candidate_pool.csv"],
              handoff_out={"n_train": len(df_train), "n_cand": len(df_candidates)})
    return df_train, df_candidates, all_comps, cand_comps


# -----------------------------------------------------------------------
# STEP 3: Surrogate model training + cross-validation
# -----------------------------------------------------------------------
def step3_train_surrogate(df_train: pd.DataFrame) -> HEASurrogateModel:
    print("\n" + "="*60)
    print("STEP 3: Surrogate GP model training")
    print("="*60)

    surrogate = HEASurrogateModel(features=DESCRIPTOR_COLS, random_state=42)

    # Cross-validate
    cv_scores = surrogate.cross_validate(df_train, df_train)
    print("  Cross-validation R² (5-fold):")
    for prop, r2 in cv_scores.items():
        print(f"    {prop}: R² = {r2:.4f}")

    # Full fit
    surrogate.fit(df_train, df_train)

    # Save CV results
    df_cv = pd.DataFrame({"property": list(cv_scores.keys()),
                           "R2_cv": list(cv_scores.values())})
    df_cv.to_csv("results/surrogate_cv_scores.csv", index=False)

    log_event("SurrogateTraining", "handoff_completed",
              files_written=["results/surrogate_cv_scores.csv"],
              handoff_out=cv_scores)
    return surrogate


# -----------------------------------------------------------------------
# STEP 4: Multi-objective Bayesian optimization
# -----------------------------------------------------------------------
def step4_bayesian_optimization(df_train: pd.DataFrame,
                                 all_comps: list,
                                 n_bo_iter: int = 15,
                                 batch_size: int = 4) -> MultiObjectiveBayesianOptimizer:
    print("\n" + "="*60)
    print("STEP 4: Multi-objective Bayesian optimization (qEHVI)")
    print("="*60)

    # Initial training data for BO
    n_init = min(60, len(df_train))
    init_comps = all_comps[:n_init]
    init_obj = df_train[PROPERTY_COLS].values[:n_init]

    optimizer = MultiObjectiveBayesianOptimizer(
        elements=ELEMENTS_5, x_min=0.05, x_max=0.55
    )
    optimizer.initialize(init_comps, init_obj)
    hv_history = [optimizer.hypervolume()]

    simulator = HEAPropertySimulator(noise_level=0.03, random_state=99)

    print(f"  Initial hypervolume: {hv_history[0]:.2f}")

    for iteration in range(n_bo_iter):
        # Suggest next batch
        try:
            suggested_comps = optimizer.suggest_next(
                batch_size=batch_size, n_candidates=300
            )
        except Exception as e:
            print(f"    BO iteration {iteration+1} failed: {e}")
            break

        # Evaluate on surrogate (simulating experiment/DFT)
        new_descs = descriptors_dataframe(suggested_comps, T=1000.0)
        new_obj = np.column_stack([
            [simulator.yield_strength(r.to_dict()) for _, r in new_descs.iterrows()],
            [simulator.elongation(r.to_dict()) for _, r in new_descs.iterrows()],
            [simulator.pitting_potential(r.to_dict()) for _, r in new_descs.iterrows()],
        ])

        optimizer.update(suggested_comps, new_obj)
        hv = optimizer.hypervolume()
        hv_history.append(hv)

        if (iteration + 1) % 5 == 0:
            print(f"  Iter {iteration+1:3d}: HV = {hv:.2f}")

    # Extract Pareto front
    pareto_Y, pareto_comps = optimizer.pareto_front()
    print(f"\n  Pareto-optimal solutions: {len(pareto_Y)}")
    print(f"  Best strength: {pareto_Y[:,0].max():.1f} MPa")
    print(f"  Best elongation: {pareto_Y[:,1].max():.1f} %")
    print(f"  Best E_pit: {pareto_Y[:,2].max():.3f} V_SCE")

    # Save results
    df_pareto = pd.DataFrame(pareto_Y, columns=["yield_strength_MPa",
                                                   "elongation_pct",
                                                   "pitting_potential_V"])
    for i, comp in enumerate(pareto_comps):
        for el, x in comp.items():
            df_pareto.loc[i, f"x_{el}"] = x
    df_pareto.to_csv("results/pareto_front.csv", index=False)

    df_hv = pd.DataFrame({"iteration": range(len(hv_history)),
                           "hypervolume": hv_history})
    df_hv.to_csv("results/bo_hypervolume_history.csv", index=False)

    log_event("BayesianOptimization", "handoff_completed",
              files_written=["results/pareto_front.csv", "results/bo_hypervolume_history.csv"],
              handoff_out={"n_pareto": len(pareto_Y), "final_HV": hv_history[-1]})
    return optimizer, hv_history, df_pareto


# -----------------------------------------------------------------------
# STEP 5: Active learning loop
# -----------------------------------------------------------------------
def step5_active_learning(surrogate: HEASurrogateModel,
                           df_train: pd.DataFrame,
                           df_candidates: pd.DataFrame,
                           n_al_iter: int = 5,
                           n_select: int = 5):
    print("\n" + "="*60)
    print("STEP 5: Active learning loop (uncertainty sampling)")
    print("="*60)

    selector = ActiveLearningSelector(strategy="hybrid")
    simulator = HEAPropertySimulator(noise_level=0.04, random_state=77)

    al_log = []
    current_train = df_train.copy()
    current_surrogate = surrogate

    for al_iter in range(n_al_iter):
        # Select most informative candidates
        sel_desc, sel_scores = selector.select(
            df_candidates, current_surrogate,
            n_select=n_select,
            existing_descriptors=current_train
        )

        # Simulate experiment on selected candidates
        new_rows = simulator.simulate_dataset(sel_desc)
        current_train = pd.concat([current_train, new_rows], ignore_index=True)

        # Retrain surrogate
        current_surrogate = HEASurrogateModel(features=DESCRIPTOR_COLS, random_state=42)
        cv = current_surrogate.cross_validate(current_train, current_train)
        current_surrogate.fit(current_train, current_train)

        al_log.append({
            "al_iteration": al_iter + 1,
            "n_train": len(current_train),
            "mean_score": float(sel_scores.mean()),
            **{f"R2_{p}": float(r2) for p, r2 in cv.items()}
        })
        print(f"  AL iter {al_iter+1}: N_train={len(current_train)}, "
              f"R²_strength={cv[PROPERTY_COLS[0]]:.4f}")

    df_al = pd.DataFrame(al_log)
    df_al.to_csv("results/active_learning_log.csv", index=False)

    log_event("ActiveLearning", "handoff_completed",
              files_written=["results/active_learning_log.csv"],
              handoff_out={"final_n_train": len(current_train),
                           "final_R2": al_log[-1]})
    return current_surrogate, df_al


# -----------------------------------------------------------------------
# STEP 6: Literature + DFT data summary
# -----------------------------------------------------------------------
def step6_database_summary():
    print("\n" + "="*60)
    print("STEP 6: Literature & database summary")
    print("="*60)

    db = MaterialsDatabaseClient()
    df_lit = db.get_curated_cantor_data()
    df_rhea = db.get_refractory_hea_data()

    df_lit.to_csv("data/literature_cantor_hea.csv", index=False)
    df_rhea.to_csv("data/literature_refractory_hea.csv", index=False)

    print(f"  Cantor-family literature records: {len(df_lit)}")
    print(f"  Refractory HEA records: {len(df_rhea)}")
    print(f"\n  Cantor data preview:")
    print(df_lit[["alloy", "sigma_y", "elong", "E_pit"]].to_string(index=False))

    log_event("DatabaseQuery", "file_written",
              files_written=["data/literature_cantor_hea.csv", "data/literature_refractory_hea.csv"])
    return df_lit, df_rhea


# -----------------------------------------------------------------------
# STEP 7: Visualization
# -----------------------------------------------------------------------
def step7_visualize(df_gibbs, G_binary, x_scan, T_scan,
                    df_train, df_pareto, hv_history, df_al, df_lit):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.gridspec import GridSpec
    import seaborn as sns

    print("\n" + "="*60)
    print("STEP 7: Generating figures")
    print("="*60)

    plt.rcParams.update({"font.size": 11, "axes.labelsize": 12, "figure.dpi": 150})
    colors = cm.viridis(np.linspace(0.15, 0.85, 6))

    # ---- Figure 1: CALPHAD Gibbs energy ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    axes[0].plot(df_gibbs["T_K"], df_gibbs["G_mix_kJmol"], color=colors[0], lw=2)
    axes[0].axhline(0, color="gray", lw=0.8, ls="--")
    axes[0].set_xlabel("Temperature (K)")
    axes[0].set_ylabel("G_mix (kJ/mol)")
    axes[0].set_title("Gibbs Free Energy — Equimolar CrMnFeCoNi")
    axes[0].fill_between(df_gibbs["T_K"], df_gibbs["G_mix_kJmol"], 0,
                         where=df_gibbs["G_mix_kJmol"] < 0,
                         alpha=0.25, color=colors[0], label="Stable region")
    axes[0].legend()

    im = axes[1].contourf(x_scan, T_scan, G_binary, levels=20, cmap="viridis")
    plt.colorbar(im, ax=axes[1], label="G_mix (kJ/mol)")
    axes[1].set_xlabel("x_Cr (Cr-Ni binary)")
    axes[1].set_ylabel("Temperature (K)")
    axes[1].set_title("Binary Phase Diagram — Cr-Ni System")

    plt.tight_layout()
    plt.savefig("figures/fig1_calphad_thermodynamics.png", dpi=150, bbox_inches="tight")
    plt.savefig("figures/fig1_calphad_thermodynamics.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig1_calphad_thermodynamics")

    # ---- Figure 2: Descriptor correlations ----
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    desc_pairs = [
        ("delta_r", "VEC", "yield_strength_MPa"),
        ("dH_mix", "dS_mix", "elongation_pct"),
        ("Omega", "Tm_mean", "yield_strength_MPa"),
        ("delta_r", "dH_mix", "pitting_potential_V"),
        ("VEC", "dS_mix", "elongation_pct"),
        ("Tm_mean", "B_Voigt", "yield_strength_MPa"),
    ]

    for ax, (xc, yc, color_prop) in zip(axes.flat, desc_pairs):
        sc = ax.scatter(df_train[xc], df_train[yc],
                        c=df_train[color_prop], cmap="viridis", s=15, alpha=0.7)
        plt.colorbar(sc, ax=ax, label=color_prop.replace("_", " "))
        ax.set_xlabel(xc)
        ax.set_ylabel(yc)
        ax.set_title(f"{xc} vs {yc}")

    plt.suptitle("HEA Descriptor Space — Property Colored", fontsize=13)
    plt.tight_layout()
    plt.savefig("figures/fig2_descriptor_correlations.png", dpi=150, bbox_inches="tight")
    plt.savefig("figures/fig2_descriptor_correlations.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig2_descriptor_correlations")

    # ---- Figure 3: Pareto front (3D objectives) ----
    fig = plt.figure(figsize=(14, 5))
    gs = GridSpec(1, 3, figure=fig)

    obj_pairs = [
        ("yield_strength_MPa", "elongation_pct"),
        ("yield_strength_MPa", "pitting_potential_V"),
        ("elongation_pct", "pitting_potential_V"),
    ]
    for i, (ax_i, (ox, oy)) in enumerate(zip(gs, obj_pairs)):
        ax = fig.add_subplot(ax_i)
        # All BO points (use training data as proxy)
        ax.scatter(df_train[ox].values, df_train[oy].values,
                   c="lightblue", s=12, alpha=0.4, label="Search space")
        ax.scatter(df_pareto[ox], df_pareto[oy],
                   c="crimson", s=50, zorder=5, label="Pareto front", marker="D")
        ax.set_xlabel(ox.replace("_", " "))
        ax.set_ylabel(oy.replace("_", " "))
        ax.set_title(f"Pareto Front: {ox.split('_')[0]} vs {oy.split('_')[0]}")
        ax.legend(fontsize=8)

    plt.suptitle("Multi-Objective Pareto Front — HEA Design", fontsize=13)
    plt.tight_layout()
    plt.savefig("figures/fig3_pareto_front.png", dpi=150, bbox_inches="tight")
    plt.savefig("figures/fig3_pareto_front.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig3_pareto_front")

    # ---- Figure 4: BO convergence + AL improvement ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    axes[0].plot(range(len(hv_history)), hv_history, "o-", color=colors[2], lw=2, ms=5)
    axes[0].set_xlabel("BO Iteration")
    axes[0].set_ylabel("Hypervolume Indicator")
    axes[0].set_title("Bayesian Optimization Convergence (qEHVI)")
    axes[0].fill_between(range(len(hv_history)), hv_history, hv_history[0],
                          alpha=0.2, color=colors[2])

    r2_cols = [c for c in df_al.columns if c.startswith("R2_")]
    for col in r2_cols:
        label = col.replace("R2_", "").replace("_pct", " (%)").replace("_MPa", " (MPa)").replace("_V", " (V)")
        axes[1].plot(df_al["al_iteration"], df_al[col], "o-", label=label, lw=2, ms=5)
    axes[1].set_xlabel("Active Learning Iteration")
    axes[1].set_ylabel("Cross-Validation R²")
    axes[1].set_title("Surrogate Model Improvement via Active Learning")
    axes[1].legend(fontsize=8)
    axes[1].set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig("figures/fig4_convergence.png", dpi=150, bbox_inches="tight")
    plt.savefig("figures/fig4_convergence.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig4_convergence")

    # ---- Figure 5: Literature validation + descriptor heatmap ----
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Scatter: literature σ_y vs VEC
    axes[0].scatter(df_lit["VEC"], df_lit["sigma_y"],
                    c=df_lit["E_pit"], cmap="RdYlGn", s=90, zorder=5, edgecolors="k", lw=0.5)
    im = axes[0].scatter(df_lit["VEC"], df_lit["sigma_y"],
                          c=df_lit["E_pit"], cmap="RdYlGn", s=90)
    plt.colorbar(im, ax=axes[0], label="E_pit (V_SCE)")
    axes[0].set_xlabel("VEC")
    axes[0].set_ylabel("Yield Strength (MPa)")
    axes[0].set_title("Literature CrMnFeCoNi-family HEAs\n(color = pitting potential)")

    # Heatmap of descriptor correlations
    feat_sub = ["delta_r", "VEC", "dH_mix", "dS_mix", "Omega", "Tm_mean"]
    corr = df_train[feat_sub + PROPERTY_COLS].corr()
    mask = np.ones_like(corr, dtype=bool)
    mask[len(feat_sub):, :len(feat_sub)] = False
    mask[:len(feat_sub), len(feat_sub):] = False
    sub_corr = corr.loc[feat_sub, PROPERTY_COLS]

    sns.heatmap(sub_corr, ax=axes[1], annot=True, fmt=".2f",
                cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                linewidths=0.5, cbar_kws={"label": "Pearson r"})
    axes[1].set_title("Descriptor–Property Correlations")
    axes[1].set_xticklabels(
        ["Yield Str.", "Elongation", "E_pit"], rotation=30, ha="right"
    )

    plt.tight_layout()
    plt.savefig("figures/fig5_literature_validation.png", dpi=150, bbox_inches="tight")
    plt.savefig("figures/fig5_literature_validation.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig5_literature_validation")

    # ---- Figure 6: Recommended compositions ----
    fig, ax = plt.subplots(figsize=(10, 5))
    top5 = df_pareto.nlargest(5, "yield_strength_MPa").reset_index(drop=True)
    el_cols = [c for c in top5.columns if c.startswith("x_")]
    labels = [f"Alloy {i+1}" for i in range(len(top5))]
    x = np.arange(len(labels))
    width = 0.12
    for j, ec in enumerate(el_cols):
        el = ec.replace("x_", "")
        offset = (j - len(el_cols)/2) * width
        ax.bar(x + offset, top5[ec], width, label=el)
    ax.set_xticks(x)
    ax.set_xticklabels([
        f"{labels[i]}\n({top5.loc[i,'yield_strength_MPa']:.0f} MPa)"
        for i in range(len(top5))
    ])
    ax.set_ylabel("Composition (mole fraction)")
    ax.set_title("Top-5 Pareto-Optimal HEA Compositions (by Yield Strength)")
    ax.legend(title="Element", fontsize=9, loc="upper right")
    ax.set_ylim(0, 0.7)

    plt.tight_layout()
    plt.savefig("figures/fig6_recommended_compositions.png", dpi=150, bbox_inches="tight")
    plt.savefig("figures/fig6_recommended_compositions.pdf", bbox_inches="tight")
    plt.close()
    print("  Saved fig6_recommended_compositions")

    log_event("Visualization", "file_written",
              files_written=[f"figures/fig{i}_{n}.png" for i, n in [
                  (1,"calphad_thermodynamics"), (2,"descriptor_correlations"),
                  (3,"pareto_front"), (4,"convergence"),
                  (5,"literature_validation"), (6,"recommended_compositions")
              ]])
    print("  All figures saved to figures/")


# -----------------------------------------------------------------------
# STEP 8: Experimental recommendations
# -----------------------------------------------------------------------
def step8_recommendations(df_pareto: pd.DataFrame) -> pd.DataFrame:
    print("\n" + "="*60)
    print("STEP 8: Experimental recommendations")
    print("="*60)

    # Score by weighted sum (normalized objectives)
    p_cols = ["yield_strength_MPa", "elongation_pct", "pitting_potential_V"]
    weights = [0.45, 0.30, 0.25]  # strength priority for high-temp HEA

    df = df_pareto.copy()
    for pc in p_cols:
        df[f"norm_{pc}"] = (df[pc] - df[pc].min()) / (df[pc].max() - df[pc].min() + 1e-9)

    df["composite_score"] = sum(w * df[f"norm_{pc}"] for w, pc in zip(weights, p_cols))
    df.sort_values("composite_score", ascending=False, inplace=True)

    top_n = min(10, len(df))
    df_rec = df.head(top_n).reset_index(drop=True)
    df_rec.to_csv("results/top_recommended_compositions.csv", index=False)

    print("  Top-5 recommended compositions:")
    el_cols = [c for c in df_rec.columns if c.startswith("x_")]
    for i, row in df_rec.head(5).iterrows():
        comp_str = " ".join([f"{c.replace('x_','')}:{row[c]:.3f}" for c in el_cols if row[c] > 0.01])
        print(f"  #{i+1}: {comp_str}")
        print(f"       σ_y={row['yield_strength_MPa']:.1f} MPa  "
              f"ε_f={row['elongation_pct']:.1f}%  "
              f"E_pit={row['pitting_potential_V']:.3f} V  "
              f"Score={row['composite_score']:.3f}")

    log_event("Recommendations", "file_written",
              files_written=["results/top_recommended_compositions.csv"],
              handoff_out={"n_recommendations": len(df_rec)})
    return df_rec


# -----------------------------------------------------------------------
# Main orchestration
# -----------------------------------------------------------------------
def main():
    log_event("Pipeline", "run_started",
              handoff_in={"task": "HEA composition optimization via ML",
                           "elements": ELEMENTS_5})

    # Step 1: CALPHAD
    df_gibbs, G_binary, x_scan, T_scan = step1_calphad()

    # Step 2: Data generation
    df_train, df_candidates, all_comps, cand_comps = step2_generate_data()

    # Step 3: Surrogate training
    surrogate = step3_train_surrogate(df_train)

    # Step 4: Bayesian optimization
    optimizer, hv_history, df_pareto = step4_bayesian_optimization(
        df_train, all_comps, n_bo_iter=8, batch_size=4
    )

    # Step 5: Active learning
    final_surrogate, df_al = step5_active_learning(surrogate, df_train, df_candidates)

    # Step 6: Database
    df_lit, df_rhea = step6_database_summary()

    # Step 7: Visualization
    step7_visualize(df_gibbs, G_binary, x_scan, T_scan,
                    df_train, df_pareto, hv_history, df_al, df_lit)

    # Step 8: Recommendations
    df_rec = step8_recommendations(df_pareto)

    log_event("Pipeline", "run_completed",
              files_written=["report.md"],
              handoff_out={"status": "complete", "n_pareto": len(df_pareto)})

    return {
        "df_train": df_train, "df_pareto": df_pareto,
        "hv_history": hv_history, "df_al": df_al,
        "df_rec": df_rec, "df_lit": df_lit, "df_rhea": df_rhea,
        "surrogate": final_surrogate,
    }


if __name__ == "__main__":
    results = main()
