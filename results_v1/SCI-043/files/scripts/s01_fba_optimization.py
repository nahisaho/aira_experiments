#!/usr/bin/env python3
"""
Module 1: FBA制約条件設定最適化
- 標準FBA / pFBA / FVA
- グルコース・酸素感度解析
"""
import cobra
from cobra.io import load_model
import pandas as pd
import numpy as np
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def run_fba_optimization():
    print("=" * 60)
    print("Module 1: FBA Constraint Optimization")
    print("=" * 60)
    model = load_model("textbook")
    print(f"Model: {model.id}, Reactions: {len(model.reactions)}, "
          f"Metabolites: {len(model.metabolites)}, Genes: {len(model.genes)}")
    results = {}

    # 1.1 Standard FBA
    sol_std = model.optimize()
    results["standard_fba"] = {"objective_value": round(sol_std.objective_value, 6), "status": sol_std.status}
    print(f"\n[1.1] Standard FBA growth rate: {sol_std.objective_value:.6f} h⁻¹")

    # 1.2 pFBA
    sol_pfba = cobra.flux_analysis.pfba(model)
    total_std = sol_std.fluxes.abs().sum()
    total_pfba = sol_pfba.fluxes.abs().sum()
    results["pfba"] = {
        "objective_value": round(sol_pfba.objective_value, 6),
        "total_flux_standard": round(total_std, 4),
        "total_flux_pfba": round(total_pfba, 4),
        "flux_reduction_pct": round((1 - total_pfba / total_std) * 100, 2),
    }
    print(f"[1.2] pFBA total flux: {total_pfba:.4f} (vs std: {total_std:.4f}, "
          f"reduction: {results['pfba']['flux_reduction_pct']:.2f}%)")

    # 1.3 FVA
    fva_result = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=0.9)
    fva_result["range"] = fva_result["maximum"] - fva_result["minimum"]
    top_variable = fva_result.nlargest(10, "range")
    results["fva"] = {
        "fraction_of_optimum": 0.9, "mean_range": round(fva_result["range"].mean(), 4),
        "max_range": round(fva_result["range"].max(), 4),
        "top_variable_reactions": top_variable.index.tolist(),
    }
    print(f"[1.3] FVA: mean range = {results['fva']['mean_range']:.4f}")

    # 1.4 Glucose sensitivity
    glc_rxn = model.reactions.get_by_id("EX_glc__D_e")
    orig_lb = glc_rxn.lower_bound
    uptake_rates = np.linspace(-20, -1, 20)
    growth_rates = []
    for ur in uptake_rates:
        with model:
            glc_rxn.lower_bound = ur
            sol = model.optimize()
            growth_rates.append(sol.objective_value if sol.status == "optimal" else 0.0)
    glc_rxn.lower_bound = orig_lb

    # 1.5 O2 sensitivity
    o2_rxn = model.reactions.get_by_id("EX_o2_e")
    orig_o2 = o2_rxn.lower_bound
    o2_rates = np.linspace(-30, 0, 31)
    growth_o2 = []
    for o2r in o2_rates:
        with model:
            o2_rxn.lower_bound = o2r
            sol = model.optimize()
            growth_o2.append(sol.objective_value if sol.status == "optimal" else 0.0)
    o2_rxn.lower_bound = orig_o2

    results["glucose_sensitivity"] = {"max_growth": round(max(growth_rates), 6)}
    results["oxygen_sensitivity"] = {"aerobic": round(growth_o2[0], 6), "anaerobic": round(growth_o2[-1], 6)}

    # Figures
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes[0,0].plot(-np.array(uptake_rates), growth_rates, 'b-o', ms=4)
    axes[0,0].set_xlabel("Glucose Uptake Rate (mmol/gDW/h)"); axes[0,0].set_ylabel("Growth Rate (h⁻¹)")
    axes[0,0].set_title("Glucose Uptake vs Growth Rate"); axes[0,0].grid(True, alpha=0.3)

    axes[0,1].plot(-np.array(o2_rates), growth_o2, 'r-o', ms=4)
    axes[0,1].set_xlabel("O₂ Uptake Rate (mmol/gDW/h)"); axes[0,1].set_ylabel("Growth Rate (h⁻¹)")
    axes[0,1].set_title("Oxygen Uptake vs Growth Rate"); axes[0,1].grid(True, alpha=0.3)

    axes[1,0].barh(top_variable.index, top_variable["range"], color="steelblue")
    axes[1,0].set_xlabel("Flux Range (mmol/gDW/h)"); axes[1,0].set_title("Top 10 Variable Reactions (FVA)")
    axes[1,0].invert_yaxis()

    active = pd.DataFrame({"Standard": sol_std.fluxes.abs(), "pFBA": sol_pfba.fluxes.abs()})
    active = active[(active > 1e-6).any(axis=1)].nlargest(20, "Standard")
    x = np.arange(len(active))
    axes[1,1].bar(x-0.2, active["Standard"], 0.4, label="Standard FBA", color="coral")
    axes[1,1].bar(x+0.2, active["pFBA"], 0.4, label="pFBA", color="steelblue")
    axes[1,1].set_xticks(x); axes[1,1].set_xticklabels(active.index, rotation=90, fontsize=7)
    axes[1,1].set_ylabel("|Flux|"); axes[1,1].set_title("Standard FBA vs pFBA"); axes[1,1].legend()

    plt.tight_layout()
    plt.savefig("figures/01_fba_optimization.png", dpi=300, bbox_inches="tight")
    plt.savefig("figures/01_fba_optimization.svg", bbox_inches="tight"); plt.close()
    fva_result.to_csv("results/01_fva_results.csv")
    with open("results/01_fba_optimization.json", "w") as f: json.dump(results, f, indent=2)
    print("[Saved] figures/01_*, results/01_*")
    return results

if __name__ == "__main__":
    run_fba_optimization()
