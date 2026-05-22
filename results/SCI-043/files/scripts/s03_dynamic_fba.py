#!/usr/bin/env python3
"""
Module 3: 動的FBA (dFBA) — SOAベース
"""
import cobra
from cobra.io import load_model
import pandas as pd
import numpy as np
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def dfba_simulation(model, initial, dt=0.1, t_end=10.0):
    Km_glc, Vmax_glc = 0.5, 10.0
    Km_o2, Vmax_o2 = 0.003, 20.0
    kLa, O2_sat = 100.0, 0.21
    X = initial.get("biomass", 0.05)
    Glc = initial.get("glucose", 20.0)
    Ac = initial.get("acetate", 0.0)
    O2 = initial.get("oxygen", 0.21)
    t, traj = 0.0, []

    while t <= t_end:
        traj.append({"time":round(t,4),"biomass":round(X,6),"glucose":round(max(Glc,0),6),
                      "acetate":round(max(Ac,0),6),"oxygen":round(max(O2,0),6)})
        if X < 1e-10 or Glc < 1e-10:
            while t <= t_end:
                t += dt
                traj.append({"time":round(t,4),"biomass":round(X,6),"glucose":round(max(Glc,0),6),
                              "acetate":round(max(Ac,0),6),"oxygen":round(max(O2,0),6)})
            break
        v_glc = Vmax_glc * Glc / (Km_glc + Glc)
        v_o2 = Vmax_o2 * O2 / (Km_o2 + O2)
        with model:
            model.reactions.get_by_id("EX_glc__D_e").lower_bound = -v_glc
            model.reactions.get_by_id("EX_o2_e").lower_bound = -v_o2
            if Glc < 0.1 and Ac > 0.1:
                try: model.reactions.get_by_id("EX_ac_e").lower_bound = -5.0*Ac/(0.5+Ac)
                except: pass
            sol = model.optimize()
            if sol.status != "optimal":
                t += dt; continue
            mu = sol.objective_value
            glc_f = sol.fluxes.get("EX_glc__D_e", 0)
            ac_f = sol.fluxes.get("EX_ac_e", 0)
            o2_f = sol.fluxes.get("EX_o2_e", 0)
        X += mu * X * dt
        Glc = max(0, Glc + glc_f * X * dt)
        Ac = max(0, Ac + ac_f * X * dt)
        O2 = max(0, O2 + (o2_f * X + kLa * (O2_sat - O2)) * dt)
        t += dt
    return pd.DataFrame(traj)

def run_dynamic_fba():
    print("=" * 60)
    print("Module 3: Dynamic FBA (dFBA)")
    print("=" * 60)
    model = load_model("textbook")
    ic = {"biomass": 0.05, "glucose": 20.0, "acetate": 0.0, "oxygen": 0.21}
    df = dfba_simulation(model, ic, dt=0.05, t_end=12.0)
    print(f"[3.1] Simulation: {len(df)} timepoints")

    # Phase detection
    glc_dep = df.loc[df["glucose"]<0.1,"time"].min() if (df["glucose"]<0.1).any() else -1
    results = {
        "summary": {
            "max_biomass": round(df["biomass"].max(), 4),
            "glucose_depletion_time": round(glc_dep, 2) if glc_dep > 0 else "N/A",
            "max_acetate": round(df["acetate"].max(), 4),
            "final_biomass": round(df["biomass"].iloc[-1], 4),
        },
        "phases": [
            {"phase": "Exponential Growth", "end_time": round(glc_dep, 2) if glc_dep > 0 else "N/A"},
            {"phase": "Diauxic Shift", "start_time": round(glc_dep, 2) if glc_dep > 0 else "N/A"},
        ]
    }
    print(f"[3.2] Max biomass: {results['summary']['max_biomass']} gDW/L")
    print(f"     Glucose depleted: {results['summary']['glucose_depletion_time']} h")

    # Multiple glucose levels
    glc_levels = [5, 10, 20, 40]
    multi = {}
    for g in glc_levels:
        d = dfba_simulation(model, {"biomass":0.05,"glucose":float(g),"acetate":0,"oxygen":0.21}, dt=0.05, t_end=15.0)
        multi[g] = d
        dep = d.loc[d["glucose"]<0.1,"time"].min() if (d["glucose"]<0.1).any() else -1
        results[f"glucose_{g}mM"] = {"max_biomass": round(d["biomass"].max(),4), "depletion_time": round(dep,2) if dep>0 else "N/A"}

    # Figures
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax1 = axes[0,0]; ax1.plot(df["time"], df["biomass"], 'b-', lw=2, label="Biomass")
    ax1.set_xlabel("Time (h)"); ax1.set_ylabel("Biomass (gDW/L)", color='b')
    ax2 = ax1.twinx(); ax2.plot(df["time"], df["glucose"], 'r--', lw=2, label="Glucose")
    ax2.set_ylabel("Glucose (mM)", color='r')
    ax1.set_title("dFBA: Biomass & Glucose"); ax1.grid(True,alpha=0.3)
    h1,l1=ax1.get_legend_handles_labels(); h2,l2=ax2.get_legend_handles_labels()
    ax1.legend(h1+h2,l1+l2,loc='center right')

    axes[0,1].plot(df["time"], df["acetate"], 'g-', lw=2, label="Acetate")
    axes[0,1].plot(df["time"], df["oxygen"]*100, 'm--', lw=2, label="O₂ (×100)")
    axes[0,1].set_xlabel("Time (h)"); axes[0,1].set_title("Acetate & O₂"); axes[0,1].legend(); axes[0,1].grid(True,alpha=0.3)

    colors = plt.cm.viridis(np.linspace(0.2,0.8,len(glc_levels)))
    for i,(g,d) in enumerate(multi.items()):
        axes[1,0].plot(d["time"], d["biomass"], color=colors[i], lw=2, label=f"Glc₀={g}")
        axes[1,1].plot(d["time"], d["glucose"], color=colors[i], lw=2, label=f"Glc₀={g}")
    axes[1,0].set_xlabel("Time (h)"); axes[1,0].set_ylabel("Biomass (gDW/L)")
    axes[1,0].set_title("Effect of Initial Glucose"); axes[1,0].legend(); axes[1,0].grid(True,alpha=0.3)
    axes[1,1].set_xlabel("Time (h)"); axes[1,1].set_ylabel("Glucose (mM)")
    axes[1,1].set_title("Glucose Consumption"); axes[1,1].legend(); axes[1,1].grid(True,alpha=0.3)

    plt.tight_layout()
    plt.savefig("figures/03_dynamic_fba.png", dpi=300, bbox_inches="tight")
    plt.savefig("figures/03_dynamic_fba.svg", bbox_inches="tight"); plt.close()
    df.to_csv("results/03_dfba_trajectory.csv", index=False)
    with open("results/03_dfba_results.json","w") as f: json.dump(results,f,indent=2,default=str)
    print("[Saved] figures/03_*, results/03_*")
    return results

if __name__ == "__main__":
    run_dynamic_fba()
