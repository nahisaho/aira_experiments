#!/usr/bin/env python3
"""
Module 5: 条件特異的モデル構築（RNA-seq統合 / GIMME）
"""
import cobra
from cobra.io import load_model
from cobra.flux_analysis import pfba
import pandas as pd
import numpy as np
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
np.random.seed(42)

def gen_rnaseq(model, cond="aerobic"):
    expr = {g.id: np.random.lognormal(3, 1.5) for g in model.genes}
    if cond == "aerobic":
        for r in model.reactions:
            if r.id in ["CS","ACONTa","ACONTb","ICDHyr","AKGDH","SUCOAS","SUCDi","FUM","MDH","CYTBD","NADH16","ATPS4r"]:
                for g in r.genes: expr[g.id] *= 5
    elif cond == "anaerobic":
        for r in model.reactions:
            if r.id in ["CS","ACONTa","ACONTb","ICDHyr","AKGDH","CYTBD","NADH16"]:
                for g in r.genes: expr[g.id] *= 0.1
            if r.id in ["LDH_D","PFL","ACALD","ALCD2x","PTAr","ACKr"]:
                for g in r.genes: expr[g.id] *= 8
    elif cond == "stress":
        for g in expr: expr[g] *= 0.3
        for r in model.reactions:
            if r.id in ["G6PDH2r","GND","TKT1","TKT2","TALA"]:
                for g in r.genes: expr[g.id] *= 10
    return expr

def gimme(model, expr, thresh_pct=25, min_growth_frac=0.1):
    thresh = np.percentile(list(expr.values()), thresh_pct)
    rxn_expr = {}
    for r in model.reactions:
        if r.genes:
            ge = [expr.get(g.id, thresh) for g in r.genes]
            rxn_expr[r.id] = max(ge) if ge else thresh
        else:
            rxn_expr[r.id] = thresh * 2
    with model:
        std_g = model.optimize().objective_value
        for r in model.reactions:
            if "Biomass" in r.id: r.lower_bound = std_g * min_growth_frac
        for rid, e in rxn_expr.items():
            if e < thresh:
                r = model.reactions.get_by_id(rid)
                s = e / thresh
                if r.lower_bound < 0: r.lower_bound *= s
                if r.upper_bound > 0: r.upper_bound *= s
        sol = pfba(model)
    active = {rid for rid in sol.fluxes.index if abs(sol.fluxes[rid]) > 1e-6}
    # Get actual growth rate from biomass reaction flux
    growth = sol.fluxes.get("Biomass_Ecoli_core", 0)
    return {"solution": sol, "active": active, "threshold": thresh,
            "growth_rate": growth, "n_active": len(active), "n_total": len(model.reactions)}

def run_condition_specific():
    print("=" * 60)
    print("Module 5: Condition-Specific Model (GIMME + RNA-seq)")
    print("=" * 60)
    model = load_model("textbook")
    conditions = ["aerobic", "anaerobic", "stress"]
    cond_res = {}; results = {}

    for c in conditions:
        print(f"\n[5.{conditions.index(c)+1}] Condition: {c}")
        e = gen_rnaseq(model, c)
        with model:
            if c == "anaerobic":
                model.reactions.get_by_id("EX_o2_e").lower_bound = 0
                model.reactions.get_by_id("EX_o2_e").upper_bound = 0
            gr = gimme(model, e)
        cond_res[c] = gr
        results[c] = {"growth_rate": round(gr["growth_rate"],6), "active_reactions": gr["n_active"],
                       "total_reactions": gr["n_total"]}
        print(f"     Growth: {gr['growth_rate']:.6f}, Active: {gr['n_active']}/{gr['n_total']}")

    # Jaccard matrix
    ov = np.zeros((3,3))
    for i,c1 in enumerate(conditions):
        for j,c2 in enumerate(conditions):
            s1,s2 = cond_res[c1]["active"],cond_res[c2]["active"]
            ov[i,j] = len(s1&s2)/len(s1|s2) if s1|s2 else 0
    results["jaccard"] = {f"{c1}_vs_{c2}": round(ov[i,j],4)
                          for i,c1 in enumerate(conditions) for j,c2 in enumerate(conditions) if i<=j}

    for c in conditions:
        sp = cond_res[c]["active"].copy()
        for o in conditions:
            if o!=c: sp -= cond_res[o]["active"]
        results[f"{c}_specific"] = list(sp)[:15]
        print(f"     {c}-specific reactions: {len(sp)}")

    # Figures
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    cols = ["#2196F3","#FF5722","#FF9800"]
    axes[0,0].bar(conditions, [results[c]["growth_rate"] for c in conditions], color=cols, edgecolor='black')
    axes[0,0].set_ylabel("Growth Rate (h⁻¹)"); axes[0,0].set_title("Growth Rate by Condition"); axes[0,0].grid(True,alpha=0.3,axis='y')

    ac = [results[c]["active_reactions"] for c in conditions]
    ic = [results[c]["total_reactions"]-results[c]["active_reactions"] for c in conditions]
    x = np.arange(3)
    axes[0,1].bar(x,ac,label="Active",color="steelblue"); axes[0,1].bar(x,ic,bottom=ac,label="Inactive",color="lightcoral")
    axes[0,1].set_xticks(x); axes[0,1].set_xticklabels(conditions); axes[0,1].set_title("Active vs Inactive"); axes[0,1].legend()

    im = axes[1,0].imshow(ov,cmap="YlOrRd",vmin=0,vmax=1)
    axes[1,0].set_xticks(range(3)); axes[1,0].set_xticklabels(conditions)
    axes[1,0].set_yticks(range(3)); axes[1,0].set_yticklabels(conditions)
    for i in range(3):
        for j in range(3): axes[1,0].text(j,i,f"{ov[i,j]:.2f}",ha="center",va="center",fontsize=12)
    axes[1,0].set_title("Jaccard Similarity"); plt.colorbar(im,ax=axes[1,0])

    kr = ["PFK","CS","G6PDH2r","PPC","ATPS4r","PYK"]
    fm = np.array([[abs(cond_res[c]["solution"].fluxes.get(r,0)) for c in conditions] for r in kr])
    w = 0.25; xp = np.arange(len(kr))
    for j,c in enumerate(conditions):
        axes[1,1].bar(xp+j*w, fm[:,j], w, label=c, color=cols[j])
    axes[1,1].set_xticks(xp+w); axes[1,1].set_xticklabels(kr)
    axes[1,1].set_title("Key Pathway Fluxes"); axes[1,1].legend(); axes[1,1].grid(True,alpha=0.3,axis='y')

    plt.tight_layout()
    plt.savefig("figures/05_condition_specific.png", dpi=300, bbox_inches="tight")
    plt.savefig("figures/05_condition_specific.svg", bbox_inches="tight"); plt.close()
    for c in conditions:
        pd.Series(gen_rnaseq(model,c)).to_csv(f"data/05_rnaseq_{c}.csv")
    with open("results/05_condition_specific.json","w") as f: json.dump(results,f,indent=2,default=str)
    print("[Saved] figures/05_*, results/05_*, data/05_*")
    return results

if __name__ == "__main__":
    run_condition_specific()
