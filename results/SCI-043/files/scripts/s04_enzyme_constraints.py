#!/usr/bin/env python3
"""
Module 4: 酵素容量制約（GECKO/sMOMENT）
"""
import cobra
from cobra.io import load_model
import pandas as pd
import numpy as np
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
np.random.seed(42)

KCAT_DATA = {
    "PFK":{"kcat":110,"mw":139,"name":"Phosphofructokinase"},
    "PYK":{"kcat":230,"mw":200,"name":"Pyruvate kinase"},
    "PGK":{"kcat":640,"mw":40.9,"name":"Phosphoglycerate kinase"},
    "ENO":{"kcat":87,"mw":93,"name":"Enolase"},
    "GAPD":{"kcat":245,"mw":145,"name":"G3P dehydrogenase"},
    "TPI":{"kcat":8600,"mw":53.9,"name":"Triose-P isomerase"},
    "FBA":{"kcat":17,"mw":156,"name":"F-bisP aldolase"},
    "PGI":{"kcat":1220,"mw":124,"name":"G6P isomerase"},
    "CS":{"kcat":65,"mw":192,"name":"Citrate synthase"},
    "ACONTa":{"kcat":18,"mw":189,"name":"Aconitase a"},
    "ACONTb":{"kcat":18,"mw":189,"name":"Aconitase b"},
    "ICDHyr":{"kcat":80,"mw":182,"name":"Isocitrate DH"},
    "AKGDH":{"kcat":48,"mw":400,"name":"α-KG DH"},
    "SUCOAS":{"kcat":89,"mw":140,"name":"Succinyl-CoA syn"},
    "SUCDi":{"kcat":580,"mw":260,"name":"Succinate DH"},
    "FUM":{"kcat":310,"mw":200,"name":"Fumarase"},
    "MDH":{"kcat":560,"mw":128,"name":"Malate DH"},
    "G6PDH2r":{"kcat":320,"mw":107,"name":"G6P DH"},
    "GND":{"kcat":120,"mw":156,"name":"6PG DH"},
    "PPC":{"kcat":69,"mw":396,"name":"PEP carboxylase"},
}
TOTAL_ENZYME = 0.096  # g/gDW — calibrated to full-growth enzyme requirement of core model

def solve_ec_fba(model, budget_frac=1.0):
    budget = TOTAL_ENZYME * budget_frac  # g/gDW
    sol_unc = model.optimize()
    mu_max, mu_min = sol_unc.objective_value, 0.0
    best_growth, best_alloc = 0, {}

    # cost: g_enzyme per (mmol/gDW/h) of flux
    # E (g/gDW) = v (mmol/gDW/h) * MW (g/mol) / (kcat (h⁻¹) * 1000 (mmol/mol))
    costs = {}
    for rid, p in KCAT_DATA.items():
        if rid in [r.id for r in model.reactions]:
            kcat_h = p["kcat"] * 3600  # s⁻¹ → h⁻¹
            mw_gpermol = p["mw"] * 1000  # kDa → g/mol
            cost = mw_gpermol / (kcat_h * 1000)  # g/(mmol/gDW/h)
            costs[rid] = {"cost": cost, "name": p["name"]}

    for _ in range(40):
        mu = (mu_max + mu_min) / 2
        with model:
            bm = model.reactions.get_by_id("Biomass_Ecoli_core")
            bm.lower_bound = mu
            bm.upper_bound = mu
            from cobra.flux_analysis import pfba
            try:
                sol = pfba(model)
            except Exception:
                mu_max = mu; continue
            total = sum(abs(sol.fluxes.get(r,0))*c["cost"] for r,c in costs.items())
            if total <= budget:
                mu_min = mu; best_growth = mu
                best_alloc = {r: {"flux":round(abs(sol.fluxes.get(r,0)),4),
                                   "enzyme_mg":round(abs(sol.fluxes.get(r,0))*c["cost"]*1000,4)}
                              for r,c in costs.items()}
            else:
                mu_max = mu
    return best_growth, best_alloc

def run_enzyme_constraints():
    print("=" * 60)
    print("Module 4: Enzyme Capacity Constraints (GECKO/sMOMENT)")
    print("=" * 60)
    model = load_model("textbook")
    sol_std = model.optimize()
    print(f"[4.1] Unconstrained growth: {sol_std.objective_value:.6f} h⁻¹")

    fracs = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5]
    budget_results = []
    for f in fracs:
        g, a = solve_ec_fba(model, f)
        red = round((1-g/sol_std.objective_value)*100, 2)
        budget_results.append({"fraction":f, "budget_mg":round(TOTAL_ENZYME*f,1),
                               "growth":round(g,6), "reduction_pct":red})
        print(f"     Budget {f}x: growth={g:.6f} ({red}% reduction)")

    g1, alloc1 = solve_ec_fba(model, 1.0)
    top = sorted(alloc1.items(), key=lambda x: x[1]["enzyme_mg"], reverse=True)[:10]

    results = {
        "standard_fba": {"growth_rate": round(sol_std.objective_value, 6)},
        "budget_analysis": budget_results,
        "enzyme_allocation_1x": {"growth_rate": round(g1, 6)},
        "top_enzyme_consumers": [{"reaction":k, **v} for k,v in top],
    }

    # Figures
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes[0,0].plot([r["budget_mg"] for r in budget_results],[r["growth"] for r in budget_results],'b-o',lw=2,ms=8)
    axes[0,0].axhline(sol_std.objective_value,color='r',ls='--',label=f"Unconstrained: {sol_std.objective_value:.4f}")
    axes[0,0].set_xlabel("Enzyme Budget (mg/gDW)"); axes[0,0].set_ylabel("Growth Rate (h⁻¹)")
    axes[0,0].set_title("Growth Rate vs Enzyme Budget"); axes[0,0].legend(); axes[0,0].grid(True,alpha=0.3)

    labels = [f"{k}" for k,v in top[:8]]; sizes = [v["enzyme_mg"] for k,v in top[:8]]
    other = sum(v["enzyme_mg"] for k,v in top[8:]); 
    if other>0: labels.append("Others"); sizes.append(other)
    if sum(sizes) > 0:
        axes[0,1].pie(sizes,labels=labels,autopct='%1.1f%%',textprops={'fontsize':8})
    axes[0,1].set_title("Enzyme Mass Allocation (1x Budget)")

    kcat_vals = [v["kcat"] for v in KCAT_DATA.values()]
    axes[1,0].hist(np.log10(kcat_vals),bins=12,color='teal',edgecolor='black',alpha=0.7)
    axes[1,0].set_xlabel("log₁₀(kcat) [s⁻¹]"); axes[1,0].set_title("kcat Distribution"); axes[1,0].grid(True,alpha=0.3)

    rxns = list(KCAT_DATA.keys())
    axes[1,1].bar(np.arange(len(rxns))-0.2,[abs(sol_std.fluxes.get(r,0)) for r in rxns],0.4,label="Unconstrained",color="coral")
    axes[1,1].set_xticks(np.arange(len(rxns))); axes[1,1].set_xticklabels(rxns,rotation=90,fontsize=7)
    axes[1,1].set_ylabel("|Flux|"); axes[1,1].set_title("Flux Distribution"); axes[1,1].legend()

    plt.tight_layout()
    plt.savefig("figures/04_enzyme_constraints.png", dpi=300, bbox_inches="tight")
    plt.savefig("figures/04_enzyme_constraints.svg", bbox_inches="tight"); plt.close()
    with open("results/04_enzyme_constraints.json","w") as f: json.dump(results,f,indent=2,default=str)
    print("[Saved] figures/04_*, results/04_*")
    return results

if __name__ == "__main__":
    run_enzyme_constraints()
