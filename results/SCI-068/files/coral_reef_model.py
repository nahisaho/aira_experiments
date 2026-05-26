#!/usr/bin/env python3
"""
Integrated Ocean Acidification – Coral Reef Ecosystem Model
============================================================
Modules:
  1. Seawater CO2 carbonate equilibrium (CO2SYS-style)
  2. Coral calcification rate (pH / Ω_arag dependence)
  3. Species interaction network (predation, competition, symbiosis)
  4. Temperature–pH compound stress (synergistic effects)
  5. Population genetics / evolutionary adaptation
  6. Great Barrier Reef 2100 projection scenarios
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import networkx as nx
import os, json

OUT = "figures"
os.makedirs(OUT, exist_ok=True)

# ── Global constants ──────────────────────────────────────────
R = 8.314462          # J mol⁻¹ K⁻¹
F_CONST = 96485.0     # C mol⁻¹

# ══════════════════════════════════════════════════════════════
# MODULE 1 – Seawater CO2 carbonate chemistry
# ══════════════════════════════════════════════════════════════

def K1_Lueker(T_K, S):
    """First dissociation constant of carbonic acid (Lueker 2000)."""
    pK1 = (3633.86/T_K - 61.2172 + 9.67770*np.log(T_K)
            - 0.011555*S + 1.0956e-4*S**2)
    return 10**(-pK1)

def K2_Lueker(T_K, S):
    """Second dissociation constant."""
    pK2 = (471.78/T_K + 25.9290 - 3.16967*np.log(T_K)
            - 0.01781*S + 1.122e-4*S**2)
    return 10**(-pK2)

def Ksp_aragonite(T_K, S):
    """Solubility product of aragonite (Mucci 1983)."""
    logKsp = (-171.945 - 0.077993*T_K + 2903.293/T_K
              + 71.595*np.log10(T_K)
              + (-0.068393 + 0.0017276*T_K + 88.135/T_K)*S**0.5
              - 0.10018*S + 0.0059415*S**1.5)
    return 10**(logKsp)

def carbonate_system(T_C, S, DIC, TA):
    """
    Solve carbonate system given T(°C), salinity, DIC, TA (all mol/kg).
    Returns dict with pH, pCO2, [CO3²⁻], Ω_arag, etc.
    """
    T_K = T_C + 273.15
    K1 = K1_Lueker(T_K, S)
    K2 = K2_Lueker(T_K, S)
    Ksp = Ksp_aragonite(T_K, S)
    Ca = 0.01028 * S / 35.0  # [Ca²⁺] mol/kg

    # Solve for [H⁺] from TA = [HCO3⁻] + 2[CO3²⁻] + [OH⁻] - [H⁺]  (simplified)
    Kw = np.exp(148.9802 - 13847.26/T_K - 23.6521*np.log(T_K)
                + (-5.977 + 118.67/T_K + 1.0495*np.log(T_K))*S**0.5
                - 0.01615*S)

    def residual(pH):
        H = 10**(-pH)
        alpha1 = K1*H / (H**2 + K1*H + K1*K2)
        alpha2 = K1*K2 / (H**2 + K1*H + K1*K2)
        HCO3 = DIC * alpha1
        CO3  = DIC * alpha2
        OH   = Kw / H
        return HCO3 + 2*CO3 + OH - H - TA

    # Find bracket where residual changes sign
    try:
        pH = brentq(residual, 6.0, 9.5)
    except ValueError:
        # Fallback: scan for sign change
        pH_scan = np.linspace(6.0, 9.5, 200)
        res_scan = [residual(p) for p in pH_scan]
        found = False
        for ii in range(len(res_scan)-1):
            if res_scan[ii]*res_scan[ii+1] < 0:
                pH = brentq(residual, pH_scan[ii], pH_scan[ii+1])
                found = True
                break
        if not found:
            pH = pH_scan[np.argmin(np.abs(res_scan))]
    H = 10**(-pH)
    alpha1 = K1*H / (H**2 + K1*H + K1*K2)
    alpha2 = K1*K2 / (H**2 + K1*H + K1*K2)
    CO2aq = DIC * (H**2 / (H**2 + K1*H + K1*K2))
    HCO3  = DIC * alpha1
    CO3   = DIC * alpha2
    omega_arag = Ca * CO3 / Ksp

    # Henry's law for pCO2
    K0 = np.exp(-60.2409 + 9345.17/T_K + 23.3585*np.log(T_K/100)
                + S*(0.023517 - 0.023656*(T_K/100) + 0.0047036*(T_K/100)**2))
    pCO2 = CO2aq / K0

    return {"pH": pH, "pCO2_uatm": pCO2*1e6,
            "CO3": CO3, "HCO3": HCO3, "CO2aq": CO2aq,
            "omega_arag": omega_arag, "DIC": DIC, "TA": TA}


# ══════════════════════════════════════════════════════════════
# MODULE 2 – Coral calcification rate model
# ══════════════════════════════════════════════════════════════

def calcification_rate(omega, T_C, omega_ref=3.5, G_max=10.0, n=1.4,
                       omega_crit=1.0, T_opt=27.0, T_sigma=3.0):
    """
    Calcification rate G (kg CaCO3 m⁻² yr⁻¹) as a function of Ω_arag and T.
    G = G_max * f(Ω) * f(T)
    f(Ω) = max(0, ((Ω - Ω_crit)/(Ω_ref - Ω_crit))^n)
    f(T) = exp(-((T - T_opt)/T_sigma)²)
    """
    raw = (omega - omega_crit)/(omega_ref - omega_crit)
    f_omega = np.where(raw > 0, raw**n, 0.0)
    f_T = np.exp(-((T_C - T_opt)/T_sigma)**2)
    return G_max * f_omega * f_T


# ══════════════════════════════════════════════════════════════
# MODULE 3 – Species interaction network
# ══════════════════════════════════════════════════════════════

SPECIES = ["HardCoral", "Zooxanthellae", "Macroalgae",
           "HerbivorousFish", "CarnivorousFish",
           "CorallivoreStarfish", "SeaUrchin"]

INTERACTIONS = [
    ("HardCoral", "Zooxanthellae", "mutualism", 0.3),
    ("Zooxanthellae", "HardCoral", "mutualism", 0.3),
    ("Macroalgae", "HardCoral", "competition", -0.2),
    ("HardCoral", "Macroalgae", "competition", -0.15),
    ("HerbivorousFish", "Macroalgae", "predation", -0.4),
    ("CarnivorousFish", "HerbivorousFish", "predation", -0.15),
    ("CorallivoreStarfish", "HardCoral", "predation", -0.25),
    ("SeaUrchin", "Macroalgae", "predation", -0.3),
]

def build_interaction_network():
    G = nx.DiGraph()
    for sp in SPECIES:
        G.add_node(sp)
    for src, tgt, itype, strength in INTERACTIONS:
        G.add_edge(src, tgt, interaction=itype, weight=strength)
    return G

def community_dynamics(t, y, params):
    """
    Lotka-Volterra community dynamics with OA modifiers.
    y = [HardCoral, Zooxanthellae, Macroalgae, HerbivorousFish,
         CarnivorousFish, CorallivoreStarfish, SeaUrchin]
    """
    r = params["r"]         # intrinsic growth rates
    K = params["K"]         # carrying capacities
    A = params["A"]         # interaction matrix 7x7
    oa_mod = params["oa_mod"]  # OA modifier per species

    dy = np.zeros(7)
    for i in range(7):
        interaction_sum = sum(A[i, j] * y[j] for j in range(7))
        dy[i] = r[i] * y[i] * (1 - y[i]/K[i]) + interaction_sum * y[i]
        dy[i] *= oa_mod[i]
        dy[i] = max(dy[i], -y[i]*0.5)  # prevent extreme negative
    return dy

def run_community_model(omega_arag, T_C, years=50):
    """Run community model with given environmental conditions."""
    # Intrinsic growth rates (yr⁻¹)
    r = np.array([0.05, 0.10, 0.15, 0.08, 0.06, 0.04, 0.07])
    K = np.array([0.6, 0.5, 0.8, 0.4, 0.3, 0.1, 0.2])

    # Interaction matrix
    A = np.zeros((7, 7))
    sp_idx = {s: i for i, s in enumerate(SPECIES)}
    for src, tgt, itype, strength in INTERACTIONS:
        A[sp_idx[tgt], sp_idx[src]] = strength * 0.1

    # OA modifiers
    oa_mod = np.ones(7)
    # Hard coral penalized by low Ω
    oa_mod[0] = np.clip(0.3 + 0.7 * (omega_arag - 1.0) / 2.5, 0.1, 1.0)
    # Zooxanthellae penalized by temperature
    oa_mod[1] = np.exp(-max(0, T_C - 29)**2 / 9)
    # Macroalgae benefit from low coral cover (already in dynamics)
    oa_mod[2] = 1.0 + 0.1 * max(0, 3.5 - omega_arag)

    y0 = np.array([0.4, 0.35, 0.15, 0.25, 0.15, 0.05, 0.12])
    params = {"r": r, "K": K, "A": A, "oa_mod": oa_mod}

    sol = solve_ivp(community_dynamics, [0, years], y0,
                    args=(params,), t_eval=np.linspace(0, years, years*4),
                    method='RK45', max_step=0.5)
    return sol


# ══════════════════════════════════════════════════════════════
# MODULE 4 – Temperature–pH compound stress
# ══════════════════════════════════════════════════════════════

def compound_stress(T_C, pH, T_opt=27.0, pH_opt=8.1):
    """
    Synergistic stress index S ∈ [0,1] where 0 = no stress, 1 = lethal.
    S = 1 - exp(-( α·ΔT² + β·ΔpH² + γ·ΔT·ΔpH ))
    """
    dT = max(0, T_C - T_opt)
    dpH = max(0, pH_opt - pH)
    alpha, beta, gamma = 0.02, 8.0, 0.5
    exponent = alpha * dT**2 + beta * dpH**2 + gamma * dT * dpH
    return 1 - np.exp(-exponent)

def bleaching_probability(T_C, DHW):
    """Bleaching probability from Degree Heating Weeks."""
    if DHW < 4:
        return 0.05
    elif DHW < 8:
        return 0.05 + 0.1 * (DHW - 4)
    else:
        return min(0.95, 0.45 + 0.0625 * (DHW - 8))

def mortality_from_stress(stress_idx, bleach_prob):
    """Annual mortality fraction combining OA stress and bleaching."""
    base_mortality = 0.01
    oa_mortality = 0.08 * stress_idx
    bleach_mortality = 0.20 * bleach_prob
    synergy = 0.05 * stress_idx * bleach_prob
    return min(0.95, base_mortality + oa_mortality + bleach_mortality + synergy)


# ══════════════════════════════════════════════════════════════
# MODULE 5 – Population genetics / evolutionary adaptation
# ══════════════════════════════════════════════════════════════

def evolutionary_response(generations, N_e, h2, sigma_p, delta_opt,
                          selection_strength=0.05):
    """
    Breeder's equation + drift: Δz̄ = h² · S + drift
    Track mean trait (thermal/pH tolerance) over generations.
    """
    z_mean = np.zeros(generations)
    z_var  = np.zeros(generations)
    z_mean[0] = 0.0
    z_var[0] = sigma_p**2
    optimal = np.linspace(0, delta_opt, generations)

    for g in range(1, generations):
        # Selection differential
        S = selection_strength * (optimal[g] - z_mean[g-1])
        # Response to selection
        response = h2 * S
        # Genetic drift
        drift = np.random.normal(0, np.sqrt(z_var[g-1] / (2*N_e)))
        z_mean[g] = z_mean[g-1] + response + drift
        # Variance erosion
        z_var[g] = z_var[g-1] * (1 - 1/(2*N_e)) + 0.001  # mutation adds variance

    fitness = np.exp(-0.5 * ((optimal - z_mean) / sigma_p)**2)
    return z_mean, z_var, optimal, fitness


# ══════════════════════════════════════════════════════════════
# MODULE 6 – GBR 2100 projection scenarios
# ══════════════════════════════════════════════════════════════

def rcp_scenario(scenario, year):
    """
    Return (T_anomaly, pCO2_atm) for a given RCP scenario and year.
    Baseline: 2020, T_base=27°C, pCO2_base=410 ppm.
    """
    t = max(0, year - 2020)
    if scenario == "RCP2.6":
        dT = 0.8 * (1 - np.exp(-t/60))
        pCO2 = 410 + 30 * (1 - np.exp(-t/40)) - 0.3*t*(t>40)
    elif scenario == "RCP4.5":
        dT = 1.5 * (1 - np.exp(-t/50))
        pCO2 = 410 + 150 * (1 - np.exp(-t/60))
    elif scenario == "RCP6.0":
        dT = 2.2 * (1 - np.exp(-t/55))
        pCO2 = 410 + 300 * (1 - np.exp(-t/65))
    else:  # RCP8.5
        dT = 3.5 * (1 - np.exp(-t/60))
        pCO2 = 410 + 3.0*t + 0.04*t**2
    pCO2 = max(280, pCO2)
    return dT, pCO2

def dic_from_pco2(pCO2_uatm, T_C, S=35.0, TA=2300e-6):
    """Estimate DIC from atmospheric pCO2 (assuming equilibrium)."""
    T_K = T_C + 273.15
    K0 = np.exp(-60.2409 + 9345.17/T_K + 23.3585*np.log(T_K/100)
                + S*(0.023517 - 0.023656*(T_K/100) + 0.0047036*(T_K/100)**2))
    CO2aq_target = K0 * pCO2_uatm * 1e-6
    K1 = K1_Lueker(T_K, S)
    K2 = K2_Lueker(T_K, S)

    # Use bisection on DIC to match target pCO2
    DIC_lo, DIC_hi = 1500e-6, 3000e-6
    for _ in range(60):
        DIC_mid = (DIC_lo + DIC_hi) / 2
        try:
            cs = carbonate_system(T_C, S, DIC_mid, TA)
            if cs["pCO2_uatm"] < pCO2_uatm:
                DIC_lo = DIC_mid
            else:
                DIC_hi = DIC_mid
        except:
            DIC_hi = DIC_mid
        if abs(DIC_hi - DIC_lo) < 1e-9:
            break
    return (DIC_lo + DIC_hi) / 2

def run_gbr_projection():
    """Run 2020–2100 projection for all RCP scenarios."""
    years = np.arange(2020, 2101)
    scenarios = ["RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"]
    S = 35.0
    TA = 2300e-6
    T_base = 27.0

    results = {}
    for sc in scenarios:
        res = {"year": years, "T": [], "pH": [], "omega": [],
               "calcification": [], "coral_cover": [],
               "stress": [], "mortality": []}
        coral_cover = 0.30  # initial 30%

        for yr in years:
            dT, pCO2 = rcp_scenario(sc, yr)
            T = T_base + dT
            DIC = dic_from_pco2(pCO2, T, S, TA)
            cs = carbonate_system(T, S, DIC, TA)

            G = calcification_rate(cs["omega_arag"], T)
            stress = compound_stress(T, cs["pH"])
            DHW = max(0, (T - 28.0)) * 4  # simplified DHW proxy
            bp = bleaching_probability(T, DHW)
            mort = mortality_from_stress(stress, bp)

            # Coral cover dynamics
            growth = 0.05 * coral_cover * (1 - coral_cover/0.6) * (G/10)
            loss = mort * coral_cover
            coral_cover = np.clip(coral_cover + growth - loss, 0.01, 0.60)

            res["T"].append(T)
            res["pH"].append(cs["pH"])
            res["omega"].append(cs["omega_arag"])
            res["calcification"].append(G)
            res["coral_cover"].append(coral_cover)
            res["stress"].append(stress)
            res["mortality"].append(mort)

        for k in res:
            if k != "year":
                res[k] = np.array(res[k])
        results[sc] = res

    return results


# ══════════════════════════════════════════════════════════════
# PLOTTING
# ══════════════════════════════════════════════════════════════

def plot_carbonate_chemistry():
    """Fig 1: Carbonate system sensitivity to DIC increase."""
    T_C, S, TA = 25.0, 35.0, 2300e-6
    DIC_range = np.linspace(1800e-6, 2400e-6, 80)
    pH_vals, omega_vals, pCO2_vals = [], [], []
    for d in DIC_range:
        cs = carbonate_system(T_C, S, d, TA)
        pH_vals.append(cs["pH"])
        omega_vals.append(cs["omega_arag"])
        pCO2_vals.append(cs["pCO2_uatm"])

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].plot(DIC_range*1e6, pH_vals, 'b-', lw=2)
    axes[0].set_xlabel("DIC (μmol/kg)"); axes[0].set_ylabel("pH")
    axes[0].set_title("pH vs DIC"); axes[0].grid(True, alpha=0.3)

    axes[1].plot(DIC_range*1e6, omega_vals, 'r-', lw=2)
    axes[1].axhline(y=1, color='k', ls='--', alpha=0.5, label='Ω=1 (undersaturation)')
    axes[1].axhline(y=3.3, color='g', ls='--', alpha=0.5, label='Pre-industrial Ω')
    axes[1].set_xlabel("DIC (μmol/kg)"); axes[1].set_ylabel("Ω_arag")
    axes[1].set_title("Aragonite Saturation vs DIC"); axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(DIC_range*1e6, pCO2_vals, 'g-', lw=2)
    axes[2].set_xlabel("DIC (μmol/kg)"); axes[2].set_ylabel("pCO₂ (μatm)")
    axes[2].set_title("pCO₂ vs DIC"); axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{OUT}/fig1_carbonate_chemistry.png", dpi=150)
    plt.close()
    print("  → fig1_carbonate_chemistry.png")

def plot_calcification():
    """Fig 2: Calcification rate as function of Ω and temperature."""
    omega_range = np.linspace(0.5, 5.0, 100)
    temps = [24, 26, 27, 28, 30, 32]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for T in temps:
        G = calcification_rate(omega_range, T)
        axes[0].plot(omega_range, G, label=f"T={T}°C", lw=1.5)
    axes[0].set_xlabel("Ω_arag"); axes[0].set_ylabel("G (kg CaCO₃ m⁻² yr⁻¹)")
    axes[0].set_title("Calcification Rate vs Ω_arag")
    axes[0].legend(fontsize=8); axes[0].grid(True, alpha=0.3)

    T_range = np.linspace(20, 35, 100)
    for omega in [2.0, 2.5, 3.0, 3.5, 4.0]:
        G = calcification_rate(omega, T_range)
        axes[1].plot(T_range, G, label=f"Ω={omega}", lw=1.5)
    axes[1].set_xlabel("Temperature (°C)"); axes[1].set_ylabel("G (kg CaCO₃ m⁻² yr⁻¹)")
    axes[1].set_title("Calcification Rate vs Temperature")
    axes[1].legend(fontsize=8); axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{OUT}/fig2_calcification.png", dpi=150)
    plt.close()
    print("  → fig2_calcification.png")

def plot_interaction_network():
    """Fig 3: Species interaction network diagram."""
    G = build_interaction_network()
    fig, ax = plt.subplots(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=42, k=2.0)

    edge_colors = []
    edge_styles = []
    for u, v, d in G.edges(data=True):
        if d["interaction"] == "mutualism":
            edge_colors.append("green")
        elif d["interaction"] == "predation":
            edge_colors.append("red")
        else:
            edge_colors.append("orange")

    node_colors = {"HardCoral": "#FF6B6B", "Zooxanthellae": "#FFD93D",
                   "Macroalgae": "#6BCB77", "HerbivorousFish": "#4D96FF",
                   "CarnivorousFish": "#1A1A2E", "CorallivoreStarfish": "#9B59B6",
                   "SeaUrchin": "#E67E22"}
    nc = [node_colors.get(n, "#999") for n in G.nodes()]

    nx.draw_networkx_nodes(G, pos, node_color=nc, node_size=2000, alpha=0.85, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=7, font_weight="bold", ax=ax)
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=2,
                           arrows=True, arrowsize=20, ax=ax,
                           connectionstyle="arc3,rad=0.1")

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0],[0], color='green', lw=2, label='Mutualism'),
        Line2D([0],[0], color='red', lw=2, label='Predation'),
        Line2D([0],[0], color='orange', lw=2, label='Competition'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10)
    ax.set_title("Coral Reef Species Interaction Network", fontsize=14)
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(f"{OUT}/fig3_interaction_network.png", dpi=150)
    plt.close()
    print("  → fig3_interaction_network.png")

def plot_community_dynamics():
    """Fig 4: Community dynamics under different OA scenarios."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    conditions = [
        ("Pre-industrial\nΩ=3.5, T=27°C", 3.5, 27),
        ("Moderate OA\nΩ=2.5, T=28.5°C", 2.5, 28.5),
        ("Severe OA\nΩ=1.5, T=30°C", 1.5, 30),
    ]

    for idx, (label, omega, T) in enumerate(conditions):
        sol = run_community_model(omega, T, years=50)
        for i, sp in enumerate(SPECIES):
            axes[idx].plot(sol.t, sol.y[i], label=sp, lw=1.5)
        axes[idx].set_title(label, fontsize=11)
        axes[idx].set_xlabel("Years"); axes[idx].set_ylabel("Relative Abundance")
        axes[idx].set_ylim(0, 0.8)
        axes[idx].grid(True, alpha=0.3)
        if idx == 0:
            axes[idx].legend(fontsize=6, loc='upper right')

    plt.tight_layout()
    plt.savefig(f"{OUT}/fig4_community_dynamics.png", dpi=150)
    plt.close()
    print("  → fig4_community_dynamics.png")

def plot_compound_stress():
    """Fig 5: Compound stress heatmap (T × pH)."""
    T_range = np.linspace(25, 34, 100)
    pH_range = np.linspace(7.4, 8.2, 100)
    T_grid, pH_grid = np.meshgrid(T_range, pH_range)
    stress_grid = np.vectorize(compound_stress)(T_grid, pH_grid)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    im = axes[0].contourf(T_grid, pH_grid, stress_grid, levels=20, cmap='YlOrRd')
    axes[0].set_xlabel("Temperature (°C)"); axes[0].set_ylabel("pH")
    axes[0].set_title("Compound Stress Index (T × pH)")
    plt.colorbar(im, ax=axes[0], label="Stress Index")
    # Mark RCP endpoints
    axes[0].plot(27.8, 8.05, 'ws', ms=10, label='RCP2.6 (2100)')
    axes[0].plot(28.5, 7.95, 'w^', ms=10, label='RCP4.5 (2100)')
    axes[0].plot(29.2, 7.85, 'wo', ms=10, label='RCP6.0 (2100)')
    axes[0].plot(30.5, 7.70, 'wD', ms=10, label='RCP8.5 (2100)')
    axes[0].legend(fontsize=7, loc='lower left')

    # Mortality surface
    stress_vals = np.linspace(0, 1, 50)
    bleach_vals = np.linspace(0, 1, 50)
    SG, BG = np.meshgrid(stress_vals, bleach_vals)
    mort_grid = np.vectorize(mortality_from_stress)(SG, BG)
    im2 = axes[1].contourf(SG, BG, mort_grid, levels=20, cmap='Reds')
    axes[1].set_xlabel("Stress Index"); axes[1].set_ylabel("Bleaching Probability")
    axes[1].set_title("Annual Mortality Rate")
    plt.colorbar(im2, ax=axes[1], label="Mortality Rate")

    plt.tight_layout()
    plt.savefig(f"{OUT}/fig5_compound_stress.png", dpi=150)
    plt.close()
    print("  → fig5_compound_stress.png")

def plot_evolutionary_response():
    """Fig 6: Evolutionary adaptation trajectories."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    np.random.seed(42)

    # Different Ne scenarios
    Ne_vals = [500, 5000, 50000]
    colors = ['#e74c3c', '#3498db', '#2ecc71']

    for i, Ne in enumerate(Ne_vals):
        z_mean, z_var, opt, fit = evolutionary_response(
            generations=80, N_e=Ne, h2=0.35, sigma_p=1.0, delta_opt=3.0)
        axes[0].plot(range(80), z_mean, color=colors[i], lw=2, label=f'Nₑ={Ne}')
        axes[1].plot(range(80), fit, color=colors[i], lw=2, label=f'Nₑ={Ne}')
        axes[2].plot(range(80), np.sqrt(z_var), color=colors[i], lw=2, label=f'Nₑ={Ne}')

    axes[0].plot(range(80), np.linspace(0, 3, 80), 'k--', alpha=0.5, label='Optimum shift')
    axes[0].set_xlabel("Generation"); axes[0].set_ylabel("Mean Trait Value")
    axes[0].set_title("Evolutionary Response"); axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.3)

    axes[1].set_xlabel("Generation"); axes[1].set_ylabel("Population Fitness")
    axes[1].set_title("Fitness Over Time"); axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.3)

    axes[2].set_xlabel("Generation"); axes[2].set_ylabel("Trait Std Dev")
    axes[2].set_title("Genetic Variation"); axes[2].legend(fontsize=8)
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{OUT}/fig6_evolutionary_response.png", dpi=150)
    plt.close()
    print("  → fig6_evolutionary_response.png")

def plot_gbr_projections():
    """Fig 7: GBR 2100 projections under RCP scenarios."""
    results = run_gbr_projection()
    scenarios = ["RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"]
    colors = {"RCP2.6": "#2ecc71", "RCP4.5": "#3498db",
              "RCP6.0": "#f39c12", "RCP8.5": "#e74c3c"}

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    titles = ["Temperature (°C)", "pH", "Ω_arag",
              "Calcification (kg m⁻² yr⁻¹)", "Coral Cover (%)", "Stress Index"]
    keys = ["T", "pH", "omega", "calcification", "coral_cover", "stress"]
    multipliers = [1, 1, 1, 1, 100, 1]

    for idx, (title, key, mult) in enumerate(zip(titles, keys, multipliers)):
        ax = axes[idx//3, idx%3]
        for sc in scenarios:
            r = results[sc]
            ax.plot(r["year"], r[key]*mult, color=colors[sc], lw=2, label=sc)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Year")
        ax.grid(True, alpha=0.3)
        if idx == 0:
            ax.legend(fontsize=8)

    plt.suptitle("Great Barrier Reef 2100 Projections", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig7_gbr_projections.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  → fig7_gbr_projections.png")

def plot_summary_dashboard():
    """Fig 8: Summary dashboard."""
    results = run_gbr_projection()

    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.3)

    # Panel A: Coral cover trajectories
    ax1 = fig.add_subplot(gs[0, 0])
    colors = {"RCP2.6": "#2ecc71", "RCP4.5": "#3498db",
              "RCP6.0": "#f39c12", "RCP8.5": "#e74c3c"}
    for sc in ["RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"]:
        r = results[sc]
        ax1.fill_between(r["year"], r["coral_cover"]*100*0.85,
                        r["coral_cover"]*100*1.15, alpha=0.15, color=colors[sc])
        ax1.plot(r["year"], r["coral_cover"]*100, color=colors[sc], lw=2, label=sc)
    ax1.set_ylabel("Coral Cover (%)"); ax1.set_xlabel("Year")
    ax1.set_title("A) Coral Cover Trajectories"); ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Panel B: 2100 bar chart
    ax2 = fig.add_subplot(gs[0, 1])
    sc_names = ["RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"]
    cover_2100 = [results[sc]["coral_cover"][-1]*100 for sc in sc_names]
    bars = ax2.bar(sc_names, cover_2100, color=[colors[sc] for sc in sc_names], alpha=0.8)
    ax2.set_ylabel("Coral Cover in 2100 (%)")
    ax2.set_title("B) Final Coral Cover by Scenario")
    for bar, val in zip(bars, cover_2100):
        ax2.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.5,
                f'{val:.1f}%', ha='center', fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')

    # Panel C: pH vs Ω phase space
    ax3 = fig.add_subplot(gs[1, 0])
    for sc in sc_names:
        r = results[sc]
        ax3.scatter(r["pH"], r["omega"], c=r["year"], cmap='viridis',
                   s=10, alpha=0.6, label=sc)
        ax3.annotate(f'{sc}\n2100', (r["pH"][-1], r["omega"][-1]),
                    fontsize=7, ha='center')
    ax3.set_xlabel("pH"); ax3.set_ylabel("Ω_arag")
    ax3.set_title("C) pH–Ω Phase Space (color = year)")
    ax3.axhline(y=1, color='r', ls='--', alpha=0.5, label='Dissolution threshold')
    ax3.grid(True, alpha=0.3)

    # Panel D: Cumulative mortality
    ax4 = fig.add_subplot(gs[1, 1])
    for sc in sc_names:
        r = results[sc]
        cum_mort = np.cumsum(r["mortality"])
        ax4.plot(r["year"], cum_mort, color=colors[sc], lw=2, label=sc)
    ax4.set_xlabel("Year"); ax4.set_ylabel("Cumulative Mortality Index")
    ax4.set_title("D) Cumulative Mortality Burden")
    ax4.legend(fontsize=8); ax4.grid(True, alpha=0.3)

    plt.savefig(f"{OUT}/fig8_summary_dashboard.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  → fig8_summary_dashboard.png")


# ══════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════

def collect_key_results():
    """Collect key numerical results for the report."""
    results = run_gbr_projection()
    kv = {}
    for sc in ["RCP2.6", "RCP4.5", "RCP6.0", "RCP8.5"]:
        r = results[sc]
        kv[sc] = {
            "T_2100": float(r["T"][-1]),
            "pH_2100": float(r["pH"][-1]),
            "omega_2100": float(r["omega"][-1]),
            "calcification_2100": float(r["calcification"][-1]),
            "coral_cover_2100_pct": float(r["coral_cover"][-1]*100),
            "stress_2100": float(r["stress"][-1]),
        }
    return kv

if __name__ == "__main__":
    print("=" * 60)
    print("Integrated Ocean Acidification – Coral Reef Model")
    print("=" * 60)

    print("\n[1/8] Carbonate chemistry sensitivity...")
    plot_carbonate_chemistry()

    print("[2/8] Calcification rate model...")
    plot_calcification()

    print("[3/8] Species interaction network...")
    plot_interaction_network()

    print("[4/8] Community dynamics...")
    plot_community_dynamics()

    print("[5/8] Compound stress analysis...")
    plot_compound_stress()

    print("[6/8] Evolutionary response...")
    plot_evolutionary_response()

    print("[7/8] GBR 2100 projections...")
    plot_gbr_projections()

    print("[8/8] Summary dashboard...")
    plot_summary_dashboard()

    print("\nCollecting key results...")
    kv = collect_key_results()
    with open("key_results.json", "w") as f:
        json.dump(kv, f, indent=2)

    print("\n" + "=" * 60)
    print("KEY RESULTS (2100 projections):")
    print("=" * 60)
    for sc, vals in kv.items():
        print(f"\n  {sc}:")
        print(f"    Temperature: {vals['T_2100']:.1f} °C")
        print(f"    pH:          {vals['pH_2100']:.3f}")
        print(f"    Ω_arag:      {vals['omega_2100']:.2f}")
        print(f"    Calcif.:     {vals['calcification_2100']:.2f} kg m⁻² yr⁻¹")
        print(f"    Coral cover: {vals['coral_cover_2100_pct']:.1f}%")
        print(f"    Stress idx:  {vals['stress_2100']:.3f}")

    print("\nAll figures saved to figures/")
    print("Done.")
