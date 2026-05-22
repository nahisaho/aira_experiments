#!/usr/bin/env python3
"""
Module 2: 13C-MFA統合フレームワーク
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

CENTRAL_RXNS = {
    "PFK": "Glycolysis", "PYK": "Glycolysis", "CS": "TCA", "ICDHyr": "TCA",
    "AKGDH": "TCA", "G6PDH2r": "PPP", "GND": "PPP", "TKT1": "PPP",
    "ME1": "Anaplerotic", "ME2": "Anaplerotic", "PPC": "Anaplerotic", "PPCK": "Anaplerotic",
}

def run_13c_mfa_integration():
    print("=" * 60)
    print("Module 2: 13C-MFA Integration with FBA")
    print("=" * 60)
    model = load_model("textbook")
    results = {}

    sol_std = model.optimize()
    print(f"[2.1] Standard FBA growth: {sol_std.objective_value:.6f} h⁻¹")

    # Simulate 13C measurements
    flux_data = {}
    for rxn_id, pathway in CENTRAL_RXNS.items():
        if rxn_id in [r.id for r in model.reactions]:
            flux_data[rxn_id] = {"pathway": pathway, "fba_flux": round(sol_std.fluxes[rxn_id], 6)}

    measured = {}
    for rxn_id, info in flux_data.items():
        true_flux = abs(info["fba_flux"])
        if true_flux > 1e-6:
            noise = np.random.normal(0, 0.05 * max(true_flux, 0.1))
            measured[rxn_id] = {
                "measured": round(max(0, true_flux + noise), 4),
                "uncertainty": round(0.05 * max(true_flux, 0.1), 4),
                "true": round(true_flux, 4),
            }
    print(f"[2.2] 13C measurements: {len(measured)} reactions")

    # Constrain model
    with model:
        for rxn_id, m in measured.items():
            rxn = model.reactions.get_by_id(rxn_id)
            sigma = m["uncertainty"] * 2
            if rxn.lower_bound >= 0:
                rxn.lower_bound = max(0, m["measured"] - sigma)
            rxn.upper_bound = min(rxn.upper_bound, m["measured"] + sigma)
        sol_13c = model.optimize()
        growth_13c = sol_13c.objective_value if sol_13c.status == "optimal" else 0

    print(f"[2.3] 13C-constrained growth: {growth_13c:.6f} h⁻¹")

    # Split ratios
    def split_ratios(sol):
        f = sol.fluxes
        r = {}
        pfk, g6p = abs(f.get("PFK",0)), abs(f.get("G6PDH2r",0))
        if pfk+g6p > 1e-6:
            r["G6P→Glycolysis"] = round(pfk/(pfk+g6p), 4)
            r["G6P→PPP"] = round(g6p/(pfk+g6p), 4)
        pyk, ppc = abs(f.get("PYK",0)), abs(f.get("PPC",0))
        if pyk+ppc > 1e-6:
            r["PEP→PYR"] = round(pyk/(pyk+ppc), 4)
            r["PEP→OAA"] = round(ppc/(pyk+ppc), 4)
        return r

    ratios_std = split_ratios(sol_std)
    ratios_13c = split_ratios(sol_13c)

    # Chi-square
    chi2, dof = 0.0, 0
    for rxn_id, m in measured.items():
        pred = abs(sol_std.fluxes.get(rxn_id, 0))
        chi2 += ((pred - m["measured"]) / max(m["uncertainty"], 0.01)) ** 2
        dof += 1
    results = {
        "standard_fba": {"growth_rate": round(sol_std.objective_value, 6)},
        "13c_constrained": {"growth_rate": round(growth_13c, 6), "n_constrained": len(measured)},
        "flux_split_ratios": {"standard": ratios_std, "13c_constrained": ratios_13c},
        "chi_square": {"value": round(chi2, 4), "dof": dof, "reduced": round(chi2/max(dof,1), 4)},
    }
    print(f"[2.4] χ²={chi2:.4f}, dof={dof}, reduced χ²={chi2/max(dof,1):.4f}")

    # Figures
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    rxn_ids = list(measured.keys())
    meas_v = [measured[r]["measured"] for r in rxn_ids]
    true_v = [measured[r]["true"] for r in rxn_ids]
    err_v = [measured[r]["uncertainty"] for r in rxn_ids]

    axes[0,0].errorbar(true_v, meas_v, yerr=err_v, fmt='o', color='steelblue', capsize=3, ms=6)
    mx = max(max(true_v), max(meas_v))*1.1
    axes[0,0].plot([0,mx],[0,mx],'k--',alpha=0.5)
    axes[0,0].set_xlabel("FBA Predicted Flux"); axes[0,0].set_ylabel("13C-MFA Measured Flux")
    axes[0,0].set_title("FBA vs 13C-MFA Flux Comparison"); axes[0,0].grid(True,alpha=0.3)

    rk = list(ratios_std.keys())
    x = np.arange(len(rk))
    axes[0,1].bar(x-0.2, [ratios_std.get(k,0) for k in rk], 0.4, label="Standard", color="coral")
    axes[0,1].bar(x+0.2, [ratios_13c.get(k,0) for k in rk], 0.4, label="13C-constrained", color="steelblue")
    axes[0,1].set_xticks(x); axes[0,1].set_xticklabels(rk, rotation=45, ha='right', fontsize=8)
    axes[0,1].set_title("Flux Split Ratios"); axes[0,1].legend(); axes[0,1].grid(True,alpha=0.3)

    central = list(flux_data.keys())
    axes[1,0].barh(central, [abs(sol_std.fluxes.get(r,0)) for r in central], color="teal")
    axes[1,0].set_xlabel("Abs Flux (mmol/gDW/h)"); axes[1,0].set_title("Central Carbon Fluxes"); axes[1,0].invert_yaxis()

    residuals = [(abs(sol_std.fluxes.get(r,0))-measured[r]["measured"])/max(measured[r]["uncertainty"],0.01) for r in rxn_ids]
    axes[1,1].bar(rxn_ids, residuals, color="mediumpurple")
    axes[1,1].axhline(0,color='k',lw=0.5); axes[1,1].axhline(2,color='r',ls='--',alpha=0.5); axes[1,1].axhline(-2,color='r',ls='--',alpha=0.5)
    axes[1,1].set_title("Residual Analysis"); axes[1,1].tick_params(axis='x',rotation=45); axes[1,1].grid(True,alpha=0.3)

    plt.tight_layout()
    plt.savefig("figures/02_13c_mfa_integration.png", dpi=300, bbox_inches="tight")
    plt.savefig("figures/02_13c_mfa_integration.svg", bbox_inches="tight"); plt.close()
    pd.DataFrame(measured).T.to_csv("results/02_13c_measured_data.csv")
    with open("results/02_13c_mfa_results.json","w") as f: json.dump(results,f,indent=2,default=str)
    print("[Saved] figures/02_*, results/02_*")
    return results

if __name__ == "__main__":
    run_13c_mfa_integration()
