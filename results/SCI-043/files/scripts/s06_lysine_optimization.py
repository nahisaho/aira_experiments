#!/usr/bin/env python3
"""
Module 6: 大腸菌リシン生産最適化ケーススタディ
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

def add_lysine_pathway(model):
    lys_c = cobra.Metabolite("lys__L_c", name="L-Lysine", compartment="c")
    lys_e = cobra.Metabolite("lys__L_e", name="L-Lysine", compartment="e")

    lys_synth = cobra.Reaction("LYS_SYNTH")
    lys_synth.name = "Lysine biosynthesis (DAP pathway)"
    lys_synth.bounds = (0, 1000)
    lys_synth.add_metabolites({
        model.metabolites.get_by_id("oaa_c"): -1,
        model.metabolites.get_by_id("pyr_c"): -1,
        model.metabolites.get_by_id("nadph_c"): -2,
        model.metabolites.get_by_id("glu__L_c"): -1,
        model.metabolites.get_by_id("succoa_c"): -1,
        model.metabolites.get_by_id("h2o_c"): -1,
        lys_c: 1,
        model.metabolites.get_by_id("co2_c"): 1,
        model.metabolites.get_by_id("nadp_c"): 2,
        model.metabolites.get_by_id("akg_c"): 1,
        model.metabolites.get_by_id("coa_c"): 1,
        model.metabolites.get_by_id("succ_c"): 1,
        model.metabolites.get_by_id("h_c"): 1,
    })

    lys_t = cobra.Reaction("LYSt"); lys_t.bounds = (0,1000)
    lys_t.add_metabolites({lys_c: -1, lys_e: 1})
    lys_ex = cobra.Reaction("EX_lys__L_e"); lys_ex.bounds = (0,1000)
    lys_ex.add_metabolites({lys_e: -1})
    model.add_reactions([lys_synth, lys_t, lys_ex])
    return model

def production_envelope(model, target, n=20):
    sol = model.optimize(); mu_max = sol.objective_value
    mus = np.linspace(0, mu_max*0.99, n)
    bm = [r for r in model.reactions if "Biomass" in r.id][0]
    maxp, minp = [], []
    for mu in mus:
        with model:
            bm.lower_bound = mu; bm.upper_bound = mu
            model.objective = target
            s = model.optimize()
            maxp.append(s.objective_value if s.status=="optimal" else 0)
            model.objective = {model.reactions.get_by_id(target): -1}
            s2 = model.optimize()
            minp.append(-s2.objective_value if s2.status=="optimal" else 0)
    return {"growth_rates": mus.tolist(), "max_prod": [max(0,v) for v in maxp], "min_prod": [max(0,v) for v in minp]}

def optknock_search(model, target, max_ko=2):
    with model:
        model.objective = target; s = model.optimize()
        wt_prod = s.objective_value if s.status=="optimal" else 0
    model.objective = "Biomass_Ecoli_core"
    sol_wt = model.optimize()
    hits = []

    # Single KOs
    for gene in model.genes:
        with model:
            gene.knock_out()
            sol = model.optimize()
            if sol.status=="optimal" and sol.objective_value > 0.01:
                g = sol.objective_value
                with model:
                    gene.knock_out()
                    [r for r in model.reactions if "Biomass" in r.id][0].lower_bound = g*0.99
                    model.objective = target
                    sp = model.optimize()
                    p = sp.objective_value if sp.status=="optimal" else 0
                if p > wt_prod * 0.1:
                    hits.append({"ko":[gene.id],"growth":round(g,6),"production":round(p,6)})

    # Double KOs
    if max_ko >= 2:
        gl = list(model.genes)
        for i in range(len(gl)):
            for j in range(i+1, len(gl)):
                with model:
                    gl[i].knock_out(); gl[j].knock_out()
                    sol = model.optimize()
                    if sol.status=="optimal" and sol.objective_value > 0.01:
                        g = sol.objective_value
                        with model:
                            gl[i].knock_out(); gl[j].knock_out()
                            [r for r in model.reactions if "Biomass" in r.id][0].lower_bound = g*0.99
                            model.objective = target
                            sp = model.optimize()
                            p = sp.objective_value if sp.status=="optimal" else 0
                        if p > wt_prod * 0.5:
                            hits.append({"ko":[gl[i].id,gl[j].id],"growth":round(g,6),"production":round(p,6)})
    hits.sort(key=lambda x: x["production"], reverse=True)
    return hits

def run_lysine_optimization():
    print("=" * 60)
    print("Module 6: Lysine Production Optimization")
    print("=" * 60)
    model = load_model("textbook")
    model = add_lysine_pathway(model)
    print(f"[6.1] Model: {len(model.reactions)} rxns, {len(model.metabolites)} mets")

    sol_wt = model.optimize()
    print(f"[6.2] WT growth: {sol_wt.objective_value:.6f} h⁻¹")

    with model:
        model.objective = "EX_lys__L_e"
        s = model.optimize(); max_lys = s.objective_value
    print(f"[6.3] Max theoretical lysine: {max_lys:.6f} mmol/gDW/h")

    with model:
        [r for r in model.reactions if "Biomass" in r.id][0].lower_bound = sol_wt.objective_value * 0.1
        model.objective = "EX_lys__L_e"
        s = model.optimize(); coupled_lys = s.objective_value if s.status=="optimal" else 0
    print(f"[6.4] Growth-coupled lysine: {coupled_lys:.6f}")

    print("[6.5] Production envelope...")
    env = production_envelope(model, "EX_lys__L_e")

    print("[6.6] OptKnock search...")
    ko_hits = optknock_search(model, "EX_lys__L_e")
    print(f"     Found {len(ko_hits)} strategies")
    for i,h in enumerate(ko_hits[:5]):
        print(f"     #{i+1}: KO={h['ko']}, growth={h['growth']:.4f}, prod={h['production']:.4f}")

    # Overexpression targets
    oe_targets = []
    for rid in ["LYS_SYNTH","CS","PPC","AKGDH","G6PDH2r","GND"]:
        try:
            with model:
                model.reactions.get_by_id(rid).upper_bound *= 2
                [r for r in model.reactions if "Biomass" in r.id][0].lower_bound = sol_wt.objective_value*0.1
                model.objective = "EX_lys__L_e"
                s = model.optimize()
                p = s.objective_value if s.status=="optimal" else 0
                oe_targets.append({"reaction":rid,"production":round(p,4),"increase":round(p-coupled_lys,4)})
        except: pass
    oe_targets.sort(key=lambda x: x["increase"], reverse=True)

    glc_up = abs(sol_wt.fluxes.get("EX_glc__D_e",-10))
    c_yield = (max_lys*6)/(glc_up*6) if glc_up>0 else 0

    results = {
        "wild_type": {"growth_rate": round(sol_wt.objective_value,6)},
        "max_theoretical_lysine": round(max_lys,6),
        "growth_coupled_lysine": round(coupled_lys,6),
        "carbon_yield": round(c_yield,4),
        "optknock_strategies": len(ko_hits),
        "top_strategies": ko_hits[:10],
        "overexpression_targets": oe_targets,
        "production_envelope": {
            "max_at_zero_growth": round(env["max_prod"][0],6),
            "max_at_50pct_growth": round(env["max_prod"][len(env["max_prod"])//2],6),
        }
    }

    # Figures
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes[0,0].fill_between(env["growth_rates"],env["min_prod"],env["max_prod"],alpha=0.3,color="steelblue")
    axes[0,0].plot(env["growth_rates"],env["max_prod"],'b-',lw=2,label="Max Lysine")
    axes[0,0].set_xlabel("Growth Rate (h⁻¹)"); axes[0,0].set_ylabel("Lysine (mmol/gDW/h)")
    axes[0,0].set_title("Production Envelope"); axes[0,0].legend(); axes[0,0].grid(True,alpha=0.3)

    if ko_hits:
        top = ko_hits[:10]
        axes[0,1].barh([",".join(k["ko"]) for k in top],[k["production"] for k in top],color="teal",edgecolor="black")
        axes[0,1].set_xlabel("Lysine (mmol/gDW/h)"); axes[0,1].set_title("Top Knockout Strategies"); axes[0,1].invert_yaxis()

    if oe_targets:
        axes[1,0].bar([t["reaction"] for t in oe_targets],[t["increase"] for t in oe_targets],
                      color=['green' if t["increase"]>0 else 'red' for t in oe_targets],edgecolor='black')
        axes[1,0].axhline(0,color='k',lw=0.5); axes[1,0].set_title("Overexpression Targets")
        axes[1,0].tick_params(axis='x',rotation=45); axes[1,0].grid(True,alpha=0.3,axis='y')

    with model:
        model.objective = "EX_lys__L_e"
        [r for r in model.reactions if "Biomass" in r.id][0].lower_bound = sol_wt.objective_value*0.1
        sol_opt = pfba(model)
    pr = ["EX_glc__D_e","PFK","PYK","CS","PPC","G6PDH2r","LYS_SYNTH","EX_lys__L_e","Biomass_Ecoli_core"]
    pf = [sol_opt.fluxes.get(r,0) for r in pr]
    axes[1,1].barh(pr,[abs(f) for f in pf],color=['red' if f<0 else 'steelblue' for f in pf],edgecolor='black')
    axes[1,1].set_xlabel("|Flux| (mmol/gDW/h)"); axes[1,1].set_title("Key Fluxes (Lys-Optimized)"); axes[1,1].invert_yaxis()

    plt.tight_layout()
    plt.savefig("figures/06_lysine_optimization.png", dpi=300, bbox_inches="tight")
    plt.savefig("figures/06_lysine_optimization.svg", bbox_inches="tight"); plt.close()
    with open("results/06_lysine_optimization.json","w") as f: json.dump(results,f,indent=2,default=str)
    print("[Saved] figures/06_*, results/06_*")
    return results

if __name__ == "__main__":
    run_lysine_optimization()
