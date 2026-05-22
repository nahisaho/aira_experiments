"""
Main Pipeline Runner: Lead-Free Perovskite High-Throughput Screening
=====================================================================
Executes the full screening pipeline and generates all results/figures.
"""

import sys
import os
import json
import time
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import FancyArrowPatch
from datetime import datetime

warnings.filterwarnings("ignore")

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from perovskite_screener.materials_database import get_all_candidates, KNOWN_BANDGAPS
from perovskite_screener.tolerance_factor import analyze_perovskite, mixed_halide_tolerance
from perovskite_screener.bandgap_ml import BandGapPredictor, compute_absorption_coefficient, slme
from perovskite_screener.defect_analysis import analyze_defects, defect_tolerance_classification
from perovskite_screener.neb_migration import run_neb, neb_temperature_dependence
from perovskite_screener.scaps_interface import run_scaps_simulation, get_perovskite_layer_params
from perovskite_screener.workflow import build_screening_workflow, save_aiida_workchain, generate_slurm_script
from perovskite_screener.ranking import rank_candidates, DEFAULT_WEIGHTS

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("results/scaps", exist_ok=True)
os.makedirs("data", exist_ok=True)
os.makedirs("logs", exist_ok=True)

LOG_FILE = "logs/process-log.jsonl"
START_TIME = datetime.now().isoformat()


def log_event(phase, event_type, skill, handoff_in=None, handoff_out=None,
              files=None, status="ok"):
    entry = {
        "timestamp": datetime.now().isoformat(),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill,
        "handoff_in": handoff_in or {},
        "handoff_out": handoff_out or {},
        "files_written": files or [],
        "status": status,
    }
    with open(LOG_FILE, "a") as f:
        f.write(json.dumps(entry, default=str) + "\n")


def save_colorblind_fig(fig, fname):
    fig.savefig(fname, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log_event("EXECUTE", "file_written", "matplotlib", files=[fname])


# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 1 – TOLERANCE FACTOR SCREENING
# ═══════════════════════════════════════════════════════════════════════════════

def run_phase1_tolerance(candidates):
    log_event("PLAN", "phase_start", "tolerance_factor",
              handoff_in={"n_candidates": len(candidates)})
    print(f"\n{'='*60}")
    print(f"PHASE 1: Tolerance Factor Screening ({len(candidates)} candidates)")
    print('='*60)

    results = []
    for c in candidates:
        try:
            r = analyze_perovskite(c["A"], c["B"], c["X"], c.get("B_ox", 2))
            results.append({
                **c,
                "goldschmidt_t": r.goldschmidt_t,
                "octahedral_mu": r.octahedral_mu,
                "bartel_tau": r.bartel_tau,
                "stability_class": r.stability_class,
                "distortion": r.distortion,
                "decomposition_risk": r.decomposition_risk,
                "phase1_pass": r.stability_class == "perovskite",
            })
        except Exception as e:
            results.append({**c, "phase1_pass": False, "error": str(e)})

    df = pd.DataFrame(results)
    passed = df[df["phase1_pass"] == True]
    print(f"  Passed: {len(passed)}/{len(candidates)} compositions")
    print(f"  Perovskite stability breakdown:")
    for sc, cnt in df["stability_class"].value_counts().items():
        print(f"    {sc:15s}: {cnt:3d}")

    df.to_csv("data/phase1_tolerance.csv", index=False)

    # Figure: t vs τ scatter
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    colors = {"Sn": "#0077BB", "Ge": "#33BBEE", "Bi": "#EE7733"}
    markers = {"I": "o", "Br": "s", "Cl": "^"}

    ax = axes[0]
    for B in ["Sn", "Ge", "Bi"]:
        sub = df[df["B"] == B]
        if sub.empty:
            continue
        for X in ["I", "Br", "Cl"]:
            s2 = sub[sub["X"] == X]
            if s2.empty:
                continue
            ax.scatter(s2["goldschmidt_t"], s2["bartel_tau"],
                       c=colors[B], marker=markers[X], s=60, alpha=0.7,
                       label=f"{B}/{X}" if X == "I" else None)

    ax.axvspan(0.80, 1.05, alpha=0.1, color="green", label="Stable t window")
    ax.axhline(4.18, color="red", ls="--", lw=1.5, label="τ = 4.18 (Bartel)")
    ax.set_xlabel("Goldschmidt Tolerance Factor $t$")
    ax.set_ylabel("Bartel $\\tau$")
    ax.set_title("Structural Stability Map")
    ax.legend(fontsize=7, ncol=2)
    ax.set_xlim(0.65, 1.15)
    ax.set_ylim(2.5, 9.0)

    ax2 = axes[1]
    stable = df[df["stability_class"] == "perovskite"]
    distortion_counts = stable["distortion"].value_counts()
    system_counts = df.groupby("B")["phase1_pass"].sum()
    wedge_colors = ["#0077BB", "#33BBEE", "#EE7733", "#BBBBBB"]
    ax2.pie(system_counts, labels=[f"{k}\n({int(v)})" for k, v in system_counts.items()],
            colors=wedge_colors[:len(system_counts)], autopct="%1.0f%%",
            startangle=90, pctdistance=0.75)
    ax2.set_title("Stable Perovskites by B-site System")

    plt.tight_layout()
    save_colorblind_fig(fig, "figures/phase1_tolerance_map.png")
    log_event("EXECUTE", "phase_complete", "tolerance_factor",
              handoff_out={"n_passed": int(len(passed))},
              files=["data/phase1_tolerance.csv", "figures/phase1_tolerance_map.png"])
    return passed.to_dict("records")


# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 2 – ML BAND GAP & ABSORPTION
# ═══════════════════════════════════════════════════════════════════════════════

def run_phase2_bandgap(candidates):
    log_event("PLAN", "phase_start", "bandgap_ml",
              handoff_in={"n_candidates": len(candidates)})
    print(f"\n{'='*60}")
    print(f"PHASE 2: DFT+ML Band Gap Screening ({len(candidates)} candidates)")
    print('='*60)

    predictor = BandGapPredictor()
    predictor.fit(verbose=True)

    print(f"  Model trained. LOO-CV MAE={predictor.cv_mae:.3f} eV, R²={predictor.cv_r2:.3f}")

    results = []
    for c in candidates:
        pred = predictor.predict(c["A"], c["B"], c["X"], c.get("B_ox", 2))
        Eg   = pred["Eg_predicted_eV"]

        # Absorption spectrum
        E_arr = np.linspace(0.3, 5.0, 500)
        direct = c["B"] != "Bi"   # Bi forms layered indirect-gap
        alpha  = compute_absorption_coefficient(Eg, E_arr, direct_gap=direct)
        slme_v = slme(Eg, alpha, E_arr, L_nm=500)

        # Window filter: 0.9–2.5 eV for useful solar absorption
        in_window = 0.9 <= Eg <= 2.5

        results.append({
            **c,
            "Eg_eV": Eg,
            "Eg_uncertainty": pred["uncertainty_eV"],
            "soc_correction": pred["soc_correction_eV"],
            "rashba_splitting": pred["rashba_splitting"],
            "direct_gap": direct,
            "slme_pct": slme_v,
            "alpha_at_1eV_above_gap": float(alpha[np.argmin(np.abs(E_arr - (Eg + 1.0)))]),
            "phase2_pass": in_window,
        })

    df = pd.DataFrame(results)
    passed = df[df["phase2_pass"] == True]
    print(f"  In 0.9–2.5 eV window: {len(passed)}/{len(candidates)}")
    print(f"\n  Top band gaps:")
    for _, row in df.nsmallest(5, "Eg_eV").iterrows():
        print(f"    {row['formula']:12s}: Eg={row['Eg_eV']:.3f} eV, SLME={row['slme_pct']:.1f}%")

    df.to_csv("data/phase2_bandgap.csv", index=False)

    # Figure: band gap distribution + parity plot
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    colors = {"Sn": "#0077BB", "Ge": "#33BBEE", "Bi": "#EE7733"}

    ax = axes[0]
    for B in ["Sn", "Ge", "Bi"]:
        sub = df[df["B"] == B]["Eg_eV"].dropna()
        ax.hist(sub, bins=12, alpha=0.65, label=B, color=colors[B], edgecolor="white")
    ax.axvline(1.34, color="red", ls="--", lw=1.5, label="SQ optimal (1.34 eV)")
    ax.axvspan(0.9, 2.5, alpha=0.08, color="green", label="Solar window")
    ax.set_xlabel("Predicted Band Gap (eV)")
    ax.set_ylabel("Count")
    ax.set_title("Band Gap Distribution by B-site")
    ax.legend(fontsize=8)

    # Parity plot (predicted vs known)
    ax2 = axes[1]
    known_X, known_y, known_labels = [], [], []
    B_ox_map = {"Pb": 2, "Sn": 2, "Ge": 2, "Bi": 3, "Sb": 3}
    for (A, B, X), data in KNOWN_BANDGAPS.items():
        try:
            pred = predictor.predict(A, B, X, B_ox_map.get(B, 2))
            known_X.append(pred["Eg_predicted_eV"])
            known_y.append(data["Eg"])
            known_labels.append(f"{A}{B}{X}")
        except Exception:
            pass
    known_X = np.array(known_X)
    known_y = np.array(known_y)
    mae = np.mean(np.abs(known_X - known_y))
    ax2.scatter(known_y, known_X, c="#0077BB", s=55, alpha=0.8, zorder=3)
    lim = [0.8, 3.2]
    ax2.plot(lim, lim, "k--", lw=1.5, label="Parity")
    ax2.fill_between(lim, [l - 0.3 for l in lim], [l + 0.3 for l in lim],
                     alpha=0.1, color="orange", label="±0.3 eV")
    for x, y, lab in zip(known_y, known_X, known_labels):
        if abs(x - y) > 0.25:
            ax2.annotate(lab, (x, y), fontsize=6, alpha=0.7)
    ax2.set_xlabel("Experimental / DFT-HSE06 Band Gap (eV)")
    ax2.set_ylabel("ML Predicted Band Gap (eV)")
    ax2.set_title(f"Parity Plot (MAE = {mae:.3f} eV)")
    ax2.legend(fontsize=8)
    ax2.set_xlim(lim); ax2.set_ylim(lim)

    # Feature importance
    ax3 = axes[2]
    fi = predictor.get_feature_importance()
    top_fi = sorted(fi.items(), key=lambda x: -x[1])[:10]
    names, vals = zip(*top_fi)
    bars = ax3.barh(range(len(names)), vals, color=cm.viridis(np.linspace(0.3, 0.9, len(names))))
    ax3.set_yticks(range(len(names)))
    ax3.set_yticklabels(names, fontsize=9)
    ax3.set_xlabel("Feature Importance")
    ax3.set_title("Top-10 ML Feature Importances")
    ax3.invert_yaxis()

    plt.tight_layout()
    save_colorblind_fig(fig, "figures/phase2_bandgap.png")

    # Absorption coefficient plot for top 5 materials
    fig2, ax = plt.subplots(figsize=(8, 5))
    top5 = df.nlargest(5, "slme_pct")
    E_plot = np.linspace(0.5, 4.0, 500)
    for _, row in top5.iterrows():
        Eg_i = row["Eg_eV"]
        direct_i = row["direct_gap"]
        alpha_i  = compute_absorption_coefficient(Eg_i, E_plot, direct_gap=direct_i)
        ax.semilogy(E_plot, alpha_i + 1, label=row["formula"])
    ax.axvline(1.34, color="gray", ls=":", lw=1)
    ax.set_xlabel("Photon Energy (eV)")
    ax.set_ylabel(r"Absorption Coefficient $\alpha$ (cm$^{-1}$)")
    ax.set_title("Optical Absorption — Top 5 Candidates")
    ax.legend(fontsize=9)
    ax.set_ylim(1e2, 2e5)
    plt.tight_layout()
    save_colorblind_fig(fig2, "figures/phase2_absorption.png")

    log_event("EXECUTE", "phase_complete", "bandgap_ml",
              handoff_out={"n_passed": int(len(passed)), "cv_mae": predictor.cv_mae},
              files=["data/phase2_bandgap.csv", "figures/phase2_bandgap.png",
                     "figures/phase2_absorption.png"])
    return passed.to_dict("records")


# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 3 – DEFECT ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

def run_phase3_defects(candidates):
    log_event("PLAN", "phase_start", "defect_analysis",
              handoff_in={"n_candidates": len(candidates)})
    print(f"\n{'='*60}")
    print(f"PHASE 3: Defect Formation Energy & Recombination ({len(candidates)})")
    print('='*60)

    results, defect_rows = [], []
    for c in candidates:
        Eg   = c.get("Eg_eV", 1.5)
        summ = analyze_defects(c["A"], c["B"], c["X"], Eg, c.get("B_ox", 2))
        tol  = defect_tolerance_classification(summ)

        results.append({
            **c,
            "defect_tolerance_score": summ.defect_tolerance_score,
            "n_deep_traps": summ.n_deep_traps,
            "Voc_nr_loss_mV": summ.Voc_nr_loss_mV,
            "dominant_defect": summ.dominant_defect,
            "defect_tolerance_class": tol,
            "phase3_pass": tol in ["defect-tolerant", "moderate"],
        })
        for d in summ.defects:
            defect_rows.append({
                "formula": c["formula"], "B": c["B"], "X": c["X"],
                "defect": d.defect_name,
                "type": d.defect_type,
                "formation_eV": d.formation_energy_eV,
                "deep_trap": d.is_deep_trap,
                "concentration": d.concentration_cm3,
                "SRH_rate": d.srh_rate_relative,
            })

    df   = pd.DataFrame(results)
    df_d = pd.DataFrame(defect_rows)
    passed = df[df["phase3_pass"] == True]
    print(f"  Defect-tolerant + moderate: {len(passed)}/{len(candidates)}")
    print(f"  Tolerance class distribution:")
    for tc, cnt in df["defect_tolerance_class"].value_counts().items():
        print(f"    {tc:20s}: {cnt:3d}")

    df.to_csv("data/phase3_defects.csv", index=False)
    df_d.to_csv("data/phase3_defect_details.csv", index=False)

    # Figure: defect formation energy heatmap
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax = axes[0]
    pivot = df.pivot_table(values="Voc_nr_loss_mV", index="B", columns="X", aggfunc="mean")
    im = ax.imshow(pivot.values, cmap="RdYlGn_r", aspect="auto", vmin=0, vmax=250)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    plt.colorbar(im, ax=ax, label="ΔVoc,nr (mV)")
    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):
            val = pivot.values[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.0f}", ha="center", va="center",
                        fontsize=9, color="black")
    ax.set_title("Non-Radiative Voc Loss (mV)")
    ax.set_xlabel("Halide (X)")
    ax.set_ylabel("B-site")

    ax2 = axes[1]
    for B, grp in df_d.groupby("B"):
        clr = {"Sn": "#0077BB", "Ge": "#33BBEE", "Bi": "#EE7733"}.get(B, "gray")
        deep = grp[grp["deep_trap"] == True]
        shallow = grp[grp["deep_trap"] == False]
        ax2.scatter(shallow["formation_eV"], np.log10(shallow["concentration"] + 1),
                    c=clr, marker="o", s=40, alpha=0.6)
        ax2.scatter(deep["formation_eV"], np.log10(deep["concentration"] + 1),
                    c=clr, marker="X", s=80, alpha=0.9, label=f"{B} deep")
    ax2.set_xlabel("Defect Formation Energy (eV)")
    ax2.set_ylabel("log₁₀ Concentration (cm⁻³)")
    ax2.set_title("Defect Formation Energy vs Concentration")
    ax2.legend(fontsize=8)
    ax2.axhline(15, color="red", ls="--", lw=1, label="N_trap = 10¹⁵")

    plt.tight_layout()
    save_colorblind_fig(fig, "figures/phase3_defects.png")
    log_event("EXECUTE", "phase_complete", "defect_analysis",
              handoff_out={"n_passed": int(len(passed))},
              files=["data/phase3_defects.csv", "figures/phase3_defects.png"])
    return passed.to_dict("records")


# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 4 – NEB ION MIGRATION
# ═══════════════════════════════════════════════════════════════════════════════

def run_phase4_neb(candidates):
    log_event("PLAN", "phase_start", "neb_migration",
              handoff_in={"n_candidates": len(candidates)})
    print(f"\n{'='*60}")
    print(f"PHASE 4: NEB Ion Migration ({len(candidates)} candidates)")
    print('='*60)

    results = []
    neb_paths = {}
    for c in candidates:
        neb = run_neb(c["A"], c["B"], c["X"], B_ox=c.get("B_ox", 2))
        results.append({
            **c,
            "barrier_eV": neb.barrier_eV,
            "D_cm2s": neb.diffusion_coeff_cm2s,
            "hop_rate_s": neb.hop_rate_s,
            "ion_mobility_risk": neb.ion_mobility_risk,
            "lit_Ea": neb.literature_Ea,
            "neb_conv_rms": neb.convergence_rms,
            "phase4_pass": neb.barrier_eV >= 0.10,
        })
        neb_paths[c["formula"]] = {
            "x": [img.position[0] for img in neb.images],
            "E": [img.energy for img in neb.images],
            "Ea": neb.barrier_eV,
        }

    df = pd.DataFrame(results)
    passed = df[df["phase4_pass"] == True]
    print(f"  Migration barrier ≥ 0.10 eV: {len(passed)}/{len(candidates)}")
    print(f"  Risk distribution:")
    for risk, cnt in df["ion_mobility_risk"].value_counts().items():
        print(f"    {risk:10s}: {cnt:3d}")

    df.to_csv("data/phase4_neb.csv", index=False)
    with open("data/neb_paths.json", "w") as f:
        json.dump(neb_paths, f, indent=2)

    # Figure: NEB profiles for key materials + barrier comparison
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    key_formulae = ["MASnI3", "CsSnI3", "MAGeI3", "CsGeI3", "MABiI3", "CsBiI3"]
    colors_neb   = plt.cm.viridis(np.linspace(0.1, 0.9, len(key_formulae)))

    ax = axes[0]
    for formula, color in zip(key_formulae, colors_neb):
        if formula in neb_paths:
            p = neb_paths[formula]
            ax.plot(p["x"], p["E"], "-o", color=color, ms=5,
                    label=f"{formula} (Ea={p['Ea']:.2f} eV)")
    ax.set_xlabel("Reaction Coordinate")
    ax.set_ylabel("Energy (eV)")
    ax.set_title("NEB Migration Profiles")
    ax.legend(fontsize=7, loc="upper right")
    ax.axhline(0, color="gray", ls=":", lw=0.8)

    ax2 = axes[1]
    T_range = np.linspace(250, 450, 100)
    for _, row in df.nlargest(3, "barrier_eV").iterrows():
        T_data = neb_temperature_dependence(row["barrier_eV"], T_range)
        ax2.semilogy(T_data["T_K"], T_data["D_cm2s"],
                     label=f"{row['formula']} (Ea={row['barrier_eV']:.2f} eV)")
    for _, row in df.nsmallest(3, "barrier_eV").iterrows():
        T_data = neb_temperature_dependence(row["barrier_eV"], T_range)
        ax2.semilogy(T_data["T_K"], T_data["D_cm2s"], "--",
                     label=f"{row['formula']} (Ea={row['barrier_eV']:.2f} eV)")
    ax2.axvline(300, color="gray", ls=":", lw=1, label="T=300K")
    ax2.set_xlabel("Temperature (K)")
    ax2.set_ylabel("Diffusion Coefficient D (cm²/s)")
    ax2.set_title("Ion Diffusion vs Temperature")
    ax2.legend(fontsize=7)

    plt.tight_layout()
    save_colorblind_fig(fig, "figures/phase4_neb.png")
    log_event("EXECUTE", "phase_complete", "neb_migration",
              handoff_out={"n_passed": int(len(passed))},
              files=["data/phase4_neb.csv", "data/neb_paths.json", "figures/phase4_neb.png"])
    return passed.to_dict("records")


# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 5 – SCAPS DEVICE SIMULATION
# ═══════════════════════════════════════════════════════════════════════════════

def run_phase5_scaps(candidates):
    log_event("PLAN", "phase_start", "scaps_interface",
              handoff_in={"n_candidates": len(candidates)})
    print(f"\n{'='*60}")
    print(f"PHASE 5: SCAPS-1D Device Simulation ({len(candidates)} candidates)")
    print('='*60)

    results = []
    for c in candidates:
        Eg     = c.get("Eg_eV", 1.5)
        dconc  = 1e15 * (2.0 if c["B"] == "Sn" else (1.5 if c["B"] == "Ge" else 0.5))
        dev    = run_scaps_simulation(c["A"], c["B"], c["X"], Eg,
                                       defect_conc=dconc,
                                       output_dir="results/scaps")
        results.append({
            **c,
            "Voc_V": dev.Voc_V,
            "Jsc_mAcm2": dev.Jsc_mAcm2,
            "FF": dev.FF,
            "device_pce_pct": dev.PCE_percent,
            "EQE_peak": dev.EQE_peak,
            "J0_mAcm2": dev.J0_mAcm2,
            "ideality_n": dev.ideality_n,
            "scaps_notes": dev.notes,
        })

    df = pd.DataFrame(results)
    print(f"  PCE range: {df['device_pce_pct'].min():.1f}% – {df['device_pce_pct'].max():.1f}%")
    print(f"  Top 5 simulated PCE:")
    for _, row in df.nlargest(5, "device_pce_pct").iterrows():
        print(f"    {row['formula']:12s}: PCE={row['device_pce_pct']:.2f}%,"
              f" Voc={row['Voc_V']:.3f}V, Jsc={row['Jsc_mAcm2']:.2f} mA/cm²,"
              f" FF={row['FF']:.3f}")

    df.to_csv("data/phase5_scaps.csv", index=False)

    # Figure: J-V curve comparison + PCE heatmap
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    top_devs = df.nlargest(6, "device_pce_pct")
    colors_dev = plt.cm.cividis(np.linspace(0.1, 0.9, len(top_devs)))

    ax = axes[0]
    for (_, row), color in zip(top_devs.iterrows(), colors_dev):
        Eg   = row["Eg_eV"]
        Voc  = row["Voc_V"]
        Jsc  = row["Jsc_mAcm2"]
        FF   = row["FF"]
        V    = np.linspace(0, Voc * 1.02, 200)
        # Ideal diode model approximation
        J0   = row["J0_mAcm2"]
        kT   = 8.617e-5 * 300
        n    = row["ideality_n"]
        J    = Jsc - J0 * (np.exp(V / (n * kT)) - 1)
        J    = np.maximum(J, 0)
        ax.plot(V, J, color=color, lw=1.8,
                label=f"{row['formula']} ({row['device_pce_pct']:.1f}%)")

    ax.set_xlabel("Voltage (V)")
    ax.set_ylabel("Current Density (mA/cm²)")
    ax.set_title("J-V Curves — Top 6 Materials")
    ax.legend(fontsize=7)
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)

    ax2 = axes[1]
    pivot_pce = df.pivot_table(values="device_pce_pct", index="B", columns="X", aggfunc="max")
    im2 = ax2.imshow(pivot_pce.values, cmap="viridis", aspect="auto", vmin=0, vmax=20)
    ax2.set_xticks(range(len(pivot_pce.columns)))
    ax2.set_xticklabels(pivot_pce.columns)
    ax2.set_yticks(range(len(pivot_pce.index)))
    ax2.set_yticklabels(pivot_pce.index)
    plt.colorbar(im2, ax=ax2, label="PCE (%)")
    for i in range(len(pivot_pce.index)):
        for j in range(len(pivot_pce.columns)):
            val = pivot_pce.values[i, j]
            if not np.isnan(val):
                ax2.text(j, i, f"{val:.1f}%", ha="center", va="center",
                         fontsize=8, color="white" if val > 10 else "black")
    ax2.set_title("Max PCE by B-site × Halide")
    ax2.set_xlabel("Halide (X)")
    ax2.set_ylabel("B-site")

    plt.tight_layout()
    save_colorblind_fig(fig, "figures/phase5_scaps.png")
    log_event("EXECUTE", "phase_complete", "scaps_interface",
              handoff_out={"max_pce": float(df["device_pce_pct"].max())},
              files=["data/phase5_scaps.csv", "figures/phase5_scaps.png"])
    return df.to_dict("records")


# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 6 – RANKING & WORKFLOW
# ═══════════════════════════════════════════════════════════════════════════════

def run_phase6_ranking(all_results_merged):
    log_event("PLAN", "phase_start", "ranking",
              handoff_in={"n_candidates": len(all_results_merged)})
    print(f"\n{'='*60}")
    print(f"PHASE 6: Multi-Criteria Ranking ({len(all_results_merged)} candidates)")
    print('='*60)

    ranked = rank_candidates(all_results_merged, top_n=None)

    print(f"\n  Top 10 Lead-Free Perovskite Candidates:")
    print(f"  {'Rank':>4} {'Formula':12} {'Score':>6} {'Eg':>5} "
          f"{'PCE%':>6} {'Voc':>5} {'Ea':>5} {'DT':>5} {'Pareto':>6}")
    print(f"  {'-'*70}")
    for c in ranked[:10]:
        pareto = "⭐" if c.pareto_optimal else " "
        print(f"  {c.rank:>4} {c.formula:12} {c.composite_score:>6.3f} "
              f"{c.Eg_eV:>5.2f} {c.device_pce_pct:>6.1f}%"
              f" {c.Voc_V:>5.3f} {c.neb_barrier_eV:>5.2f}"
              f" {c.defect_tolerance:>5.2f} {pareto:>6}")

    # Save ranking results
    rank_data = []
    for c in ranked:
        rank_data.append({
            "rank": c.rank,
            "formula": c.formula,
            "A": c.A, "B": c.B, "X": c.X,
            "system": c.system,
            "composite_score": c.composite_score,
            "Eg_eV": c.Eg_eV,
            "goldschmidt_t": c.goldschmidt_t,
            "bartel_tau": c.bartel_tau,
            "stability_class": c.stability_class,
            "defect_tolerance": c.defect_tolerance,
            "Voc_nr_loss_mV": c.voc_nr_loss_mV,
            "barrier_eV": c.neb_barrier_eV,
            "slme_pct": c.slme_pct,
            "device_pce_pct": c.device_pce_pct,
            "Voc_V": c.Voc_V,
            "Jsc_mAcm2": c.Jsc_mAcm2,
            "FF": c.FF,
            "pareto_optimal": c.pareto_optimal,
            "recommendation": c.recommendation,
            **{k: getattr(c, f"{k}_score" if not k.endswith("score") else k, 0)
               for k in ["band_gap_score", "stability_score", "ion_migration_score",
                         "slme_score", "device_pce_score", "voc_score"]},
        })

    df_rank = pd.DataFrame(rank_data)
    df_rank.to_csv("results/candidate_ranking.csv", index=False)

    # Figure: Comprehensive ranking visualization
    fig = plt.figure(figsize=(18, 11))
    top10 = ranked[:10]

    # 1. Composite score bar chart
    ax = fig.add_subplot(2, 3, 1)
    scores = [c.composite_score for c in top10]
    labels = [c.formula for c in top10]
    colors_bar = [{"Sn": "#0077BB", "Ge": "#33BBEE", "Bi": "#EE7733"}.get(c.B, "gray") for c in top10]
    bars = ax.barh(range(len(top10)), scores, color=colors_bar, edgecolor="white")
    ax.set_yticks(range(len(top10)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Composite Score")
    ax.set_title("Top-10 Ranking")
    ax.invert_yaxis()
    for bar, score in zip(bars, scores):
        ax.text(score + 0.005, bar.get_y() + bar.get_height()/2,
                f"{score:.3f}", va="center", fontsize=7)

    # 2. Radar chart for top 5 (polar projection)
    ax2 = fig.add_subplot(2, 3, 2, projection="polar")
    categories = ["Band Gap", "Stability", "Defect Tol.", "Ion Migr.", "SLME", "PCE"]
    N = len(categories)
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]
    ax2.set_theta_offset(np.pi / 2)
    ax2.set_theta_direction(-1)
    ax2.set_xticks(angles[:-1])
    ax2.set_xticklabels(categories, fontsize=7)
    ax2.set_ylim(0, 1)
    for c_i, color in zip(top10[:5], plt.cm.viridis(np.linspace(0.1, 0.9, 5))):
        vals = [c_i.band_gap_score, c_i.stability_score,
                c_i.defect_tolerance_score_n, c_i.ion_migration_score,
                c_i.slme_score, c_i.device_pce_score]
        vals += vals[:1]
        ax2.plot(angles, vals, "o-", lw=1.5, color=color, label=c_i.formula)
        ax2.fill(angles, vals, alpha=0.08, color=color)
    ax2.legend(fontsize=7, loc="lower right", bbox_to_anchor=(1.3, -0.1))
    ax2.set_title("Multi-Criteria Radar (Top 5)", pad=15)

    # 3. PCE vs Voc scatter
    ax3 = fig.add_subplot(2, 3, 3)
    for B, color in [("Sn", "#0077BB"), ("Ge", "#33BBEE"), ("Bi", "#EE7733")]:
        sub = [c for c in ranked if c.B == B]
        ax3.scatter([c.Voc_V for c in sub], [c.device_pce_pct for c in sub],
                    c=color, label=B, alpha=0.7, s=40)
    for c in top10[:5]:
        ax3.annotate(c.formula, (c.Voc_V, c.device_pce_pct),
                     fontsize=6, alpha=0.8)
    ax3.set_xlabel("Voc (V)")
    ax3.set_ylabel("PCE (%)")
    ax3.set_title("PCE vs Voc")
    ax3.legend()

    # 4. Band gap vs NEB barrier
    ax4 = fig.add_subplot(2, 3, 4)
    for B, color in [("Sn", "#0077BB"), ("Ge", "#33BBEE"), ("Bi", "#EE7733")]:
        sub = [c for c in ranked if c.B == B]
        sc = ax4.scatter([c.Eg_eV for c in sub], [c.neb_barrier_eV for c in sub],
                         c=[c.device_pce_pct for c in sub],
                         cmap="viridis", marker={"Sn": "o", "Ge": "s", "Bi": "^"}[B],
                         s=60, alpha=0.8, label=B, vmin=0, vmax=18)
    plt.colorbar(sc, ax=ax4, label="PCE (%)")
    ax4.axvspan(0.9, 2.5, alpha=0.06, color="green")
    ax4.axhline(0.22, color="red", ls="--", lw=1, label="MAPbI3 ref (0.22 eV)")
    ax4.set_xlabel("Band Gap (eV)")
    ax4.set_ylabel("NEB Migration Barrier (eV)")
    ax4.set_title("Band Gap vs Ion Migration Barrier")
    ax4.legend(fontsize=8)

    # 5. Defect tolerance vs Voc,nr loss
    ax5 = fig.add_subplot(2, 3, 5)
    for B, color in [("Sn", "#0077BB"), ("Ge", "#33BBEE"), ("Bi", "#EE7733")]:
        sub = [c for c in ranked if c.B == B]
        ax5.scatter([c.defect_tolerance for c in sub],
                    [c.voc_nr_loss_mV for c in sub],
                    c=color, label=B, alpha=0.7, s=40)
    ax5.axhline(80, color="orange", ls="--", lw=1.2, label="80 mV threshold")
    ax5.set_xlabel("Defect Tolerance Score")
    ax5.set_ylabel("Non-Radiative Voc Loss (mV)")
    ax5.set_title("Defect Tolerance vs Recombination Loss")
    ax5.legend()
    ax5.invert_yaxis()

    # 6. System comparison (Sn vs Ge vs Bi)
    ax6 = fig.add_subplot(2, 3, 6)
    systems = ["Sn", "Ge", "Bi"]
    metrics = ["band_gap_score", "stability_score", "defect_tolerance_score_n",
               "ion_migration_score", "device_pce_score"]
    metric_labels = ["Band Gap", "Stability", "Defect Tol.", "Ion Migr.", "PCE"]
    x_pos = np.arange(len(metrics))
    width = 0.25
    for i, B in enumerate(systems):
        sub = [c for c in ranked if c.B == B]
        means = [np.mean([getattr(c, m) for c in sub]) for m in metrics]
        color = {"Sn": "#0077BB", "Ge": "#33BBEE", "Bi": "#EE7733"}[B]
        ax6.bar(x_pos + i * width, means, width=width, label=B,
                color=color, alpha=0.8, edgecolor="white")
    ax6.set_xticks(x_pos + width)
    ax6.set_xticklabels(metric_labels, fontsize=8)
    ax6.set_ylabel("Mean Score (0–1)")
    ax6.set_title("System Comparison: Sn vs Ge vs Bi")
    ax6.legend()
    ax6.set_ylim(0, 1.05)

    plt.tight_layout()
    save_colorblind_fig(fig, "figures/phase6_ranking.png")
    log_event("EXECUTE", "phase_complete", "ranking",
              handoff_out={"top1": ranked[0].formula, "top1_score": ranked[0].composite_score},
              files=["results/candidate_ranking.csv", "figures/phase6_ranking.png"])
    return ranked


# ═══════════════════════════════════════════════════════════════════════════════
# WORKFLOW DESIGN
# ═══════════════════════════════════════════════════════════════════════════════

def run_workflow_design(candidates):
    print(f"\n{'='*60}")
    print("WORKFLOW DESIGN: AiiDA/FireWorks Pipeline")
    print('='*60)

    # Build FireWorks workflow
    wf = build_screening_workflow(candidates[:10])  # subset for demo
    wf.save_json("results/workflow_definition.json")
    print(f"  FireWorks workflow: {len(wf.fireworks)} Fireworks saved → results/workflow_definition.json")

    # AiiDA WorkChain
    wc_path = save_aiida_workchain("src/perovskite_screener")
    print(f"  AiiDA WorkChain saved → {wc_path}")

    # SLURM script
    slurm_path = generate_slurm_script("perovskite_screening", n_cores=32,
                                        memory_gb=64, walltime_h=48)
    print(f"  SLURM script saved → {slurm_path}")

    # Workflow diagram (simple DAG visualization)
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_xlim(0, 14); ax.set_ylim(0, 7)
    ax.axis("off")

    stages = [
        (1, 3.5, "Candidate\nGeneration\n(54 comps)", "#4477AA"),
        (3, 5.5, "Tolerance\nFilter\n(Phase 1)", "#228833"),
        (3, 2.0, "ML Band\nGap Screen\n(Phase 2)", "#228833"),
        (5, 3.5, "DFT\nRelaxation\n(Phase 3a)", "#AA3377"),
        (7, 5.5, "HSE06+SOC\nBands/DOS\n(Phase 3b)", "#AA3377"),
        (7, 2.0, "Defect\nCalcs\n(Phase 3c)", "#CCBB44"),
        (9, 5.5, "NEB Ion\nMigration\n(Phase 4)", "#CCBB44"),
        (9, 2.0, "SCAPS-1D\nDevice Sim\n(Phase 5)", "#66CCEE"),
        (11, 3.5, "Multi-Criteria\nRanking\n(Phase 6)", "#EE6677"),
        (13, 3.5, "Top-10\nCandidates\n+Report", "#AA3377"),
    ]

    for x, y, label, color in stages:
        rect = plt.Rectangle((x - 0.7, y - 0.55), 1.4, 1.1,
                              facecolor=color, alpha=0.75, edgecolor="white",
                              linewidth=1.5, zorder=3)
        ax.add_patch(rect)
        ax.text(x, y, label, ha="center", va="center", fontsize=7,
                color="white", fontweight="bold", zorder=4)

    arrows = [
        (1, 3.5, 2.3, 5.5), (1, 3.5, 2.3, 2.0),
        (3, 5.5, 4.3, 3.5), (3, 2.0, 4.3, 3.5),
        (5, 3.5, 6.3, 5.5), (5, 3.5, 6.3, 2.0),
        (7, 5.5, 8.3, 5.5), (7, 2.0, 8.3, 2.0),
        (9, 5.5, 10.3, 3.5), (9, 2.0, 10.3, 3.5),
        (11, 3.5, 12.3, 3.5),
    ]
    for x1, y1, x2, y2 in arrows:
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="->", color="gray", lw=1.5))

    # AiiDA/FireWorks label
    ax.text(7, 6.5, "AiiDA WorkChain / FireWorks DAG",
            ha="center", fontsize=10, style="italic", color="#333333")
    ax.text(7, 0.3, "SLURM/PBS Queue   |   VASP 6.4.1 + HSE06+SOC   |   Python post-processing",
            ha="center", fontsize=8, color="#555555")

    plt.tight_layout()
    save_colorblind_fig(fig, "figures/workflow_diagram.png")
    print("  Workflow diagram saved → figures/workflow_diagram.png")

    log_event("EXECUTE", "workflow_designed", "aiida_workflow",
              files=["results/workflow_definition.json",
                     "src/perovskite_screener/aiida_workchain.py",
                     "figures/workflow_diagram.png"])


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN EXECUTION
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    log_event("PLAN", "run_started", "pipeline",
              handoff_in={"timestamp": START_TIME})
    print("=" * 60)
    print(" Lead-Free Perovskite High-Throughput Screening")
    print(" Target: Sn/Ge/Bi-based ABX3 candidates")
    print("=" * 60)

    # Generate candidate pool
    candidates = get_all_candidates()
    print(f"\n  Total candidate pool: {len(candidates)} compositions")
    print(f"  Systems: Sn ({sum(1 for c in candidates if c['B']=='Sn')}), "
          f"Ge ({sum(1 for c in candidates if c['B']=='Ge')}), "
          f"Bi ({sum(1 for c in candidates if c['B']=='Bi')})")

    log_event("EXECUTE", "prompt_received", "pipeline",
              handoff_in={"n_total": len(candidates)})

    # Run phases
    p1 = run_phase1_tolerance(candidates)
    p2 = run_phase2_bandgap(p1)
    p3 = run_phase3_defects(p2)
    p4 = run_phase4_neb(p3)
    p5 = run_phase5_scaps(p4)

    # p5 already contains all upstream columns (via **c expansion in each phase)
    merged = pd.DataFrame(p5)
    # Drop suffixed duplicates if any exist from incremental merges
    dup_cols = [c for c in merged.columns if c.endswith("_x") or c.endswith("_y")]
    merged = merged.drop(columns=dup_cols, errors="ignore")
    merged.to_csv("results/all_candidates_merged.csv", index=False)

    # Ranking
    ranked = run_phase6_ranking(merged.to_dict("records"))

    # Workflow design
    run_workflow_design(candidates)

    # Summary statistics
    print(f"\n{'='*60}")
    print("SCREENING SUMMARY")
    print('='*60)
    print(f"  Total candidates:           {len(candidates):4d}")
    print(f"  After Phase 1 (tolerance):  {len(p1):4d}")
    print(f"  After Phase 2 (band gap):   {len(p2):4d}")
    print(f"  After Phase 3 (defects):    {len(p3):4d}")
    print(f"  After Phase 4 (NEB):        {len(p4):4d}")
    print(f"  After Phase 5 (device):     {len(p5):4d}")
    print(f"\n  #1 Candidate: {ranked[0].formula}")
    print(f"    Score:     {ranked[0].composite_score:.4f}")
    print(f"    Band gap:  {ranked[0].Eg_eV:.3f} eV")
    print(f"    PCE (sim): {ranked[0].device_pce_pct:.1f}%")
    print(f"    Voc:       {ranked[0].Voc_V:.3f} V")
    print(f"    NEB Ea:    {ranked[0].neb_barrier_eV:.3f} eV")
    print(f"    {ranked[0].recommendation}")

    log_event("REPORT", "run_completed", "pipeline",
              handoff_out={
                  "top1": ranked[0].formula,
                  "top1_pce": ranked[0].device_pce_pct,
                  "n_pareto": sum(1 for c in ranked if c.pareto_optimal),
              },
              files=[
                  "results/candidate_ranking.csv",
                  "results/all_candidates_merged.csv",
                  "results/workflow_definition.json",
              ])

    return ranked, merged


if __name__ == "__main__":
    ranked, merged = main()
