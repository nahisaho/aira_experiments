#!/usr/bin/env python3
"""
Systems Biology Framework for Predicting Diet-Gut Microbiome Interactions
=========================================================================
Integrates:
1. SHIME-inspired digestion/absorption kinetics
2. Generalized Lotka-Volterra (gLV) microbial community dynamics
3. SCFA flux prediction
4. Long-term dietary pattern simulation
5. Probiotic/prebiotic effect prediction
6. Fermented food case study
7. MICOM/gapseq-inspired community metabolic modeling
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

sns.set_theme(style="whitegrid", font_scale=1.1)
FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)

np.random.seed(42)

# ============================================================
# Species definitions for the gut microbiome model
# ============================================================
SPECIES = [
    "Bacteroides",       # 0
    "Faecalibacterium",  # 1 (butyrate producer)
    "Bifidobacterium",   # 2 (probiotic)
    "Roseburia",         # 3 (butyrate producer)
    "Lactobacillus",     # 4 (probiotic, fermented foods)
    "Prevotella",        # 5 (fiber degrader)
    "Clostridium",       # 6 (pathobiont)
    "Akkermansia",       # 7 (mucin degrader)
]
N_SPECIES = len(SPECIES)

# ============================================================
# 1. SHIME-inspired Digestion/Absorption Kinetics
# ============================================================
def shime_digestion(diet_composition, time_hours=48, dt=0.1):
    """
    Simulate nutrient digestion and absorption through GI compartments.
    Compartments: Stomach -> Small Intestine -> Ascending Colon -> Transverse Colon -> Descending Colon
    """
    t = np.arange(0, time_hours, dt)
    n_steps = len(t)
    
    # Diet components: [starch, fiber, protein, fat, simple_sugars]
    starch, fiber, protein, fat, sugars = diet_composition
    
    # Compartment concentrations [stomach, SI, AC, TC, DC]
    compartments = ['Stomach', 'Small Intestine', 'Ascending Colon', 
                    'Transverse Colon', 'Descending Colon']
    n_comp = 5
    
    # Transit rates (per hour)
    k_transit = np.array([0.5, 0.3, 0.15, 0.1, 0.08])
    # Digestion rates
    k_digest_starch = np.array([0.1, 0.8, 0.05, 0.02, 0.01])
    k_digest_protein = np.array([0.2, 0.6, 0.1, 0.05, 0.02])
    k_digest_fat = np.array([0.05, 0.7, 0.05, 0.02, 0.01])
    
    # Fiber is not digested in stomach/SI, only fermented in colon
    k_ferment_fiber = np.array([0.0, 0.0, 0.3, 0.2, 0.1])
    
    # Track substrate in each compartment
    results = {}
    for substrate_name, initial, k_dig in [
        ('Starch', starch, k_digest_starch),
        ('Fiber', fiber, k_ferment_fiber),
        ('Protein', protein, k_digest_protein),
        ('Fat', fat, k_digest_fat),
        ('Sugars', sugars, k_digest_starch * 1.5),
    ]:
        conc = np.zeros((n_steps, n_comp))
        conc[0, 0] = initial
        
        for i in range(1, n_steps):
            for j in range(n_comp):
                # Inflow from previous compartment
                inflow = k_transit[j-1] * conc[i-1, j-1] * dt if j > 0 else 0
                # Outflow to next compartment
                outflow = k_transit[j] * conc[i-1, j] * dt if j < n_comp - 1 else 0
                # Digestion/fermentation loss
                digestion = k_dig[j] * conc[i-1, j] * dt
                conc[i, j] = conc[i-1, j] + inflow - outflow - digestion
                conc[i, j] = max(0, conc[i, j])
        
        results[substrate_name] = conc
    
    return t, results, compartments


def plot_shime_digestion(t, results, compartments):
    """Plot SHIME digestion dynamics."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.flatten()
    
    colors = sns.color_palette("husl", 5)
    
    for idx, (substrate, conc) in enumerate(results.items()):
        ax = axes[idx]
        for j, comp in enumerate(compartments):
            ax.plot(t, conc[:, j], label=comp, linewidth=2, color=colors[j])
        ax.set_title(f'{substrate} Dynamics', fontsize=13, fontweight='bold')
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Concentration (g/L)')
        ax.legend(fontsize=8)
        ax.set_xlim(0, 48)
    
    axes[5].axis('off')
    plt.suptitle('SHIME-Inspired GI Tract Digestion Kinetics', fontsize=15, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/shime_digestion.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> Saved shime_digestion.png")


# ============================================================
# 2. Generalized Lotka-Volterra (gLV) Community Model
# ============================================================
def build_glv_parameters():
    """Build gLV interaction matrix and growth rates."""
    # Intrinsic growth rates (per day)
    mu = np.array([0.8, 0.6, 0.7, 0.55, 0.65, 0.75, 0.4, 0.5])
    
    # Interaction matrix (N x N)
    # Diagonal: self-limitation (carrying capacity)
    A = np.zeros((N_SPECIES, N_SPECIES))
    np.fill_diagonal(A, -np.array([0.01, 0.012, 0.011, 0.013, 0.012, 0.01, 0.015, 0.014]))
    
    # Cross-feeding and competition
    # Bacteroides-Faecalibacterium: cross-feeding (acetate -> butyrate)
    A[1, 0] = 0.003   # Faecal benefits from Bacteroides
    A[0, 1] = 0.001   # Slight benefit
    
    # Bifidobacterium-Roseburia: cross-feeding
    A[3, 2] = 0.004   # Roseburia benefits from Bifido (acetate/lactate)
    A[2, 3] = 0.001
    
    # Lactobacillus-Bifidobacterium: cooperation
    A[2, 4] = 0.002
    A[4, 2] = 0.002
    
    # Prevotella-Bacteroides: competition
    A[0, 5] = -0.005
    A[5, 0] = -0.005
    
    # Clostridium: inhibited by SCFA producers
    A[6, 1] = -0.004
    A[6, 2] = -0.003
    A[6, 3] = -0.004
    
    # Akkermansia: relatively independent
    A[7, 2] = 0.001  # Slight benefit from Bifido
    
    return mu, A


def glv_ode(t, x, mu, A, nutrient_modifier=None):
    """gLV ODE: dx_i/dt = x_i * (mu_i + sum_j(A_ij * x_j))"""
    x = np.maximum(x, 0)
    mu_eff = mu.copy()
    if nutrient_modifier is not None:
        mu_eff = mu_eff * nutrient_modifier
    dxdt = x * (mu_eff + A @ x)
    return dxdt


def simulate_glv(x0, mu, A, t_span, nutrient_modifier=None, t_eval=None):
    """Simulate gLV dynamics."""
    if t_eval is None:
        t_eval = np.linspace(t_span[0], t_span[1], 500)
    
    sol = solve_ivp(
        glv_ode, t_span, x0, args=(mu, A, nutrient_modifier),
        t_eval=t_eval, method='RK45', max_step=0.1,
        rtol=1e-8, atol=1e-10
    )
    return sol


def plot_glv_dynamics(sol, title_suffix="", filename="glv_dynamics.png"):
    """Plot gLV community dynamics."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    colors = sns.color_palette("Set2", N_SPECIES)
    
    for i, sp in enumerate(SPECIES):
        ax1.plot(sol.t, sol.y[i], label=sp, linewidth=2, color=colors[i])
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Abundance (cells/mL × 10⁸)')
    ax1.set_title(f'Community Dynamics {title_suffix}', fontweight='bold')
    ax1.legend(fontsize=9)
    ax1.set_ylim(bottom=0)
    
    # Relative abundance
    total = sol.y.sum(axis=0)
    total[total == 0] = 1
    for i, sp in enumerate(SPECIES):
        rel = sol.y[i] / total
        ax2.fill_between(sol.t, 
                         np.sum(sol.y[:i] / total[np.newaxis, :], axis=0) if i > 0 else np.zeros_like(sol.t),
                         np.sum(sol.y[:i+1] / total[np.newaxis, :], axis=0),
                         label=sp, color=colors[i], alpha=0.8)
    ax2.set_xlabel('Time (days)')
    ax2.set_ylabel('Relative Abundance')
    ax2.set_title(f'Relative Composition {title_suffix}', fontweight='bold')
    ax2.legend(fontsize=9, loc='center left', bbox_to_anchor=(1, 0.5))
    ax2.set_ylim(0, 1)
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/{filename}', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  -> Saved {filename}")


# ============================================================
# 3. SCFA Flux Prediction
# ============================================================
def predict_scfa(abundances, fiber_availability=1.0):
    """
    Predict SCFA production from species abundances.
    SCFA yield coefficients (mmol/cell·day) based on metabolic models.
    """
    # SCFA yield matrix: rows=species, cols=[acetate, propionate, butyrate]
    Y_scfa = np.array([
        [0.8, 0.4, 0.05],   # Bacteroides: acetate+propionate
        [0.3, 0.1, 0.9],    # Faecalibacterium: butyrate producer
        [0.7, 0.1, 0.1],    # Bifidobacterium: acetate
        [0.2, 0.05, 0.85],  # Roseburia: butyrate
        [0.5, 0.1, 0.05],   # Lactobacillus: acetate+lactate
        [0.6, 0.5, 0.1],    # Prevotella: acetate+propionate
        [0.3, 0.2, 0.3],    # Clostridium: mixed
        [0.4, 0.3, 0.15],   # Akkermansia: propionate
    ])
    
    # Scale by fiber availability
    Y_scaled = Y_scfa * fiber_availability
    
    # SCFA = abundances @ Y_scaled
    if abundances.ndim == 1:
        scfa = abundances @ Y_scaled
    else:
        scfa = abundances.T @ Y_scaled  # time x 3
    
    return scfa


def plot_scfa_dynamics(sol, filename="scfa_dynamics.png"):
    """Plot SCFA production over time."""
    scfa = predict_scfa(sol.y)  # (time, 3)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    scfa_names = ['Acetate', 'Propionate', 'Butyrate']
    colors = ['#2ecc71', '#3498db', '#e74c3c']
    
    for i, (name, color) in enumerate(zip(scfa_names, colors)):
        ax1.plot(sol.t, scfa[:, i], label=name, linewidth=2.5, color=color)
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('SCFA Production (mmol/L)')
    ax1.set_title('SCFA Production Dynamics', fontweight='bold')
    ax1.legend(fontsize=11)
    
    # Ratio plot
    total_scfa = scfa.sum(axis=1)
    total_scfa[total_scfa == 0] = 1
    for i, (name, color) in enumerate(zip(scfa_names, colors)):
        ax2.plot(sol.t, scfa[:, i] / total_scfa * 100, label=name, linewidth=2, color=color)
    ax2.set_xlabel('Time (days)')
    ax2.set_ylabel('SCFA Proportion (%)')
    ax2.set_title('SCFA Molar Ratios', fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.set_ylim(0, 100)
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/{filename}', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  -> Saved {filename}")
    return scfa


# ============================================================
# 4. Long-term Dietary Pattern Simulation
# ============================================================
def simulate_dietary_patterns(mu, A, x0, days=90):
    """Simulate different dietary patterns over 90 days."""
    diets = {
        'Western Diet': np.array([1.0, 0.8, 1.2, 1.3, 0.9, 0.7, 1.1, 0.8]),
        'High-Fiber Diet': np.array([1.0, 1.4, 1.3, 1.3, 1.1, 1.5, 0.7, 1.1]),
        'Mediterranean Diet': np.array([1.1, 1.3, 1.2, 1.1, 1.2, 1.3, 0.8, 1.2]),
        'Low-FODMAP Diet': np.array([1.0, 0.9, 0.8, 0.9, 0.9, 0.7, 1.0, 1.0]),
    }
    
    results = {}
    for diet_name, modifier in diets.items():
        sol = simulate_glv(x0, mu, A, (0, days), nutrient_modifier=modifier)
        results[diet_name] = sol
    
    return results


def plot_dietary_patterns(results, filename="dietary_patterns.png"):
    """Plot long-term dietary pattern effects."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    colors = sns.color_palette("Set2", N_SPECIES)
    
    for idx, (diet_name, sol) in enumerate(results.items()):
        ax = axes[idx // 2][idx % 2]
        total = sol.y.sum(axis=0)
        total[total == 0] = 1
        
        bottom = np.zeros_like(sol.t)
        for i, sp in enumerate(SPECIES):
            rel = sol.y[i] / total
            ax.fill_between(sol.t, bottom, bottom + rel, 
                          label=sp, color=colors[i], alpha=0.85)
            bottom += rel
        
        ax.set_title(diet_name, fontsize=13, fontweight='bold')
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Relative Abundance')
        ax.set_ylim(0, 1)
        if idx == 0:
            ax.legend(fontsize=8, loc='upper right')
    
    plt.suptitle('Long-term Dietary Pattern Effects on Gut Microbiota', 
                 fontsize=15, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/{filename}', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  -> Saved {filename}")


def plot_diet_scfa_comparison(results, filename="diet_scfa_comparison.png"):
    """Compare SCFA production across diets."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    scfa_names = ['Acetate', 'Propionate', 'Butyrate']
    diet_names = list(results.keys())
    x = np.arange(len(diet_names))
    width = 0.25
    colors = ['#2ecc71', '#3498db', '#e74c3c']
    
    for i, (scfa_name, color) in enumerate(zip(scfa_names, colors)):
        vals = []
        for diet_name, sol in results.items():
            scfa = predict_scfa(sol.y)
            vals.append(scfa[-1, i])  # Final SCFA value
        ax.bar(x + i * width, vals, width, label=scfa_name, color=color, alpha=0.85)
    
    ax.set_xticks(x + width)
    ax.set_xticklabels(diet_names, rotation=15)
    ax.set_ylabel('SCFA Production (mmol/L)')
    ax.set_title('SCFA Production by Dietary Pattern (Day 90)', fontweight='bold')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/{filename}', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  -> Saved {filename}")


# ============================================================
# 5. Probiotic/Prebiotic Effect Prediction
# ============================================================
def simulate_probiotic_prebiotic(mu, A, x0, days=60):
    """Simulate probiotic and prebiotic interventions."""
    scenarios = {}
    
    # Control
    sol_ctrl = simulate_glv(x0, mu, A, (0, days))
    scenarios['Control'] = sol_ctrl
    
    # Probiotic: add Bifidobacterium + Lactobacillus boost
    x0_prob = x0.copy()
    x0_prob[2] += 15  # Bifidobacterium boost
    x0_prob[4] += 10  # Lactobacillus boost
    sol_prob = simulate_glv(x0_prob, mu, A, (0, days))
    scenarios['Probiotic'] = sol_prob
    
    # Prebiotic (inulin/FOS): enhance fiber-degraders
    prebiotic_mod = np.array([1.0, 1.3, 1.4, 1.3, 1.1, 1.4, 0.8, 1.1])
    sol_preb = simulate_glv(x0, mu, A, (0, days), nutrient_modifier=prebiotic_mod)
    scenarios['Prebiotic (Inulin)'] = sol_preb
    
    # Synbiotic: probiotic + prebiotic
    sol_syn = simulate_glv(x0_prob, mu, A, (0, days), nutrient_modifier=prebiotic_mod)
    scenarios['Synbiotic'] = sol_syn
    
    return scenarios


def plot_probiotic_prebiotic(scenarios, filename="probiotic_prebiotic.png"):
    """Plot probiotic/prebiotic intervention results."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    colors = sns.color_palette("Set2", N_SPECIES)
    
    for idx, (name, sol) in enumerate(scenarios.items()):
        ax = axes[idx // 2][idx % 2]
        for i, sp in enumerate(SPECIES):
            ax.plot(sol.t, sol.y[i], label=sp, linewidth=2, color=colors[i])
        ax.set_title(name, fontsize=13, fontweight='bold')
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Abundance (×10⁸ cells/mL)')
        ax.set_ylim(bottom=0)
        if idx == 0:
            ax.legend(fontsize=8)
    
    plt.suptitle('Probiotic/Prebiotic Intervention Effects', fontsize=15, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/{filename}', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  -> Saved {filename}")
    
    # Shannon diversity comparison
    fig2, ax2 = plt.subplots(figsize=(10, 6))
    for name, sol in scenarios.items():
        total = sol.y.sum(axis=0)
        total[total == 0] = 1
        rel = sol.y / total[np.newaxis, :]
        rel = np.clip(rel, 1e-10, 1)
        shannon = -np.sum(rel * np.log(rel), axis=0)
        ax2.plot(sol.t, shannon, label=name, linewidth=2.5)
    
    ax2.set_xlabel('Time (days)')
    ax2.set_ylabel('Shannon Diversity Index')
    ax2.set_title('Diversity Changes Under Interventions', fontweight='bold')
    ax2.legend(fontsize=11)
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/diversity_interventions.png', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  -> Saved diversity_interventions.png")


# ============================================================
# 6. Fermented Food Case Study
# ============================================================
def simulate_fermented_food(mu, A, x0, days=90):
    """Simulate fermented food intake effects."""
    # Phase 1: baseline (0-30 days)
    # Phase 2: fermented food introduction (30-60 days)
    # Phase 3: post-intervention (60-90 days)
    
    t_eval_p1 = np.linspace(0, 30, 200)
    t_eval_p2 = np.linspace(30, 60, 200)
    t_eval_p3 = np.linspace(60, 90, 200)
    
    # Phase 1: baseline
    sol1 = simulate_glv(x0, mu, A, (0, 30), t_eval=t_eval_p1)
    
    # Phase 2: fermented food (yogurt/kimchi/kefir)
    # Boosts Lactobacillus and Bifidobacterium, increases diversity
    fermented_mod = np.array([1.0, 1.1, 1.3, 1.1, 1.4, 1.0, 0.85, 1.1])
    x0_p2 = sol1.y[:, -1].copy()
    x0_p2[4] += 8   # Lactobacillus from fermented food
    x0_p2[2] += 5   # Bifidobacterium
    sol2 = simulate_glv(x0_p2, mu, A, (30, 60), nutrient_modifier=fermented_mod, t_eval=t_eval_p2)
    
    # Phase 3: post-intervention (return to normal)
    x0_p3 = sol2.y[:, -1].copy()
    sol3 = simulate_glv(x0_p3, mu, A, (60, 90), t_eval=t_eval_p3)
    
    # Combine
    t_all = np.concatenate([sol1.t, sol2.t, sol3.t])
    y_all = np.concatenate([sol1.y, sol2.y, sol3.y], axis=1)
    
    return t_all, y_all


def plot_fermented_food(t_all, y_all, filename="fermented_food.png"):
    """Plot fermented food case study results."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    colors = sns.color_palette("Set2", N_SPECIES)
    
    # Abundance dynamics
    ax = axes[0, 0]
    for i, sp in enumerate(SPECIES):
        ax.plot(t_all, y_all[i], label=sp, linewidth=2, color=colors[i])
    ax.axvspan(30, 60, alpha=0.15, color='green', label='Fermented Food Period')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Abundance (×10⁸ cells/mL)')
    ax.set_title('Species Dynamics', fontweight='bold')
    ax.legend(fontsize=8)
    
    # Relative abundance
    ax = axes[0, 1]
    total = y_all.sum(axis=0)
    total[total == 0] = 1
    bottom = np.zeros_like(t_all)
    for i, sp in enumerate(SPECIES):
        rel = y_all[i] / total
        ax.fill_between(t_all, bottom, bottom + rel, label=sp, color=colors[i], alpha=0.85)
        bottom += rel
    ax.axvline(30, color='black', linestyle='--', linewidth=1)
    ax.axvline(60, color='black', linestyle='--', linewidth=1)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Relative Abundance')
    ax.set_title('Community Composition', fontweight='bold')
    ax.set_ylim(0, 1)
    
    # Shannon diversity
    ax = axes[1, 0]
    rel = y_all / total[np.newaxis, :]
    rel = np.clip(rel, 1e-10, 1)
    shannon = -np.sum(rel * np.log(rel), axis=0)
    ax.plot(t_all, shannon, linewidth=2.5, color='#8e44ad')
    ax.axvspan(30, 60, alpha=0.15, color='green')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Shannon Diversity Index')
    ax.set_title('α-Diversity Over Time', fontweight='bold')
    
    # SCFA dynamics
    ax = axes[1, 1]
    scfa = predict_scfa(y_all)
    scfa_names = ['Acetate', 'Propionate', 'Butyrate']
    scfa_colors = ['#2ecc71', '#3498db', '#e74c3c']
    for i, (name, color) in enumerate(zip(scfa_names, scfa_colors)):
        ax.plot(t_all, scfa[:, i], label=name, linewidth=2.5, color=color)
    ax.axvspan(30, 60, alpha=0.15, color='green')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('SCFA (mmol/L)')
    ax.set_title('SCFA Production', fontweight='bold')
    ax.legend(fontsize=11)
    
    plt.suptitle('Fermented Food Intervention Case Study', fontsize=15, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/{filename}', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  -> Saved {filename}")


# ============================================================
# 7. MICOM/gapseq-inspired Community Metabolic Modeling
# ============================================================
def community_metabolic_model():
    """
    Simplified MICOM/gapseq-inspired community metabolic model.
    Simulates metabolic exchanges and cooperative tradeoff.
    """
    # Species-specific metabolic capabilities (flux capacities)
    # Rows: species, Cols: [glucose_uptake, fiber_ferment, acetate_prod, 
    #   propionate_prod, butyrate_prod, lactate_prod, h2_prod, co2_prod]
    metabolic_matrix = np.array([
        [0.8, 0.5, 0.7, 0.5, 0.05, 0.1, 0.2, 0.3],   # Bacteroides
        [0.4, 0.6, 0.2, 0.05, 0.8, 0.05, 0.3, 0.2],   # Faecalibacterium
        [0.7, 0.3, 0.6, 0.05, 0.1, 0.5, 0.1, 0.2],    # Bifidobacterium
        [0.3, 0.5, 0.15, 0.03, 0.75, 0.05, 0.25, 0.15], # Roseburia
        [0.6, 0.2, 0.4, 0.05, 0.05, 0.7, 0.05, 0.15],  # Lactobacillus
        [0.5, 0.8, 0.5, 0.4, 0.1, 0.1, 0.2, 0.25],     # Prevotella
        [0.6, 0.3, 0.3, 0.2, 0.25, 0.1, 0.4, 0.3],     # Clostridium
        [0.2, 0.1, 0.3, 0.35, 0.1, 0.05, 0.1, 0.1],    # Akkermansia
    ])
    
    # Cooperative tradeoff simulation
    # Community growth rate optimization under different dietary conditions
    diets = {
        'Standard': np.array([10, 5, 0, 0, 0, 0, 0, 0]),
        'High-Fiber': np.array([5, 15, 0, 0, 0, 0, 0, 0]),
        'High-Protein': np.array([8, 3, 0, 0, 0, 0, 0, 0]),
        'Ketogenic': np.array([2, 2, 0, 0, 0, 0, 0, 0]),
    }
    
    results = {}
    for diet_name, diet_input in diets.items():
        # Compute growth rates as function of resource availability
        resource_availability = diet_input[:2]  # glucose, fiber
        
        growth_rates = np.zeros(N_SPECIES)
        metabolite_fluxes = np.zeros((N_SPECIES, 8))
        
        for i in range(N_SPECIES):
            # Growth proportional to resource utilization capacity
            growth = (metabolic_matrix[i, 0] * resource_availability[0] + 
                     metabolic_matrix[i, 1] * resource_availability[1])
            growth_rates[i] = growth
            metabolite_fluxes[i] = metabolic_matrix[i] * growth / 10.0
        
        # Cross-feeding: acetate and lactate as inputs for butyrate producers
        acetate_pool = metabolite_fluxes[:, 2].sum()
        lactate_pool = metabolite_fluxes[:, 5].sum()
        
        # Butyrate producers get extra growth from cross-feeding
        for i in [1, 3]:  # Faecalibacterium, Roseburia
            cross_feed_boost = (acetate_pool * 0.1 + lactate_pool * 0.05)
            growth_rates[i] += cross_feed_boost
            metabolite_fluxes[i, 4] += cross_feed_boost * 0.5
        
        results[diet_name] = {
            'growth_rates': growth_rates,
            'metabolite_fluxes': metabolite_fluxes,
            'total_scfa': metabolite_fluxes[:, 2:5].sum(axis=0),
        }
    
    return results


def plot_community_metabolic(results, filename="community_metabolic.png"):
    """Plot community metabolic modeling results."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Growth rates comparison
    ax = axes[0]
    diet_names = list(results.keys())
    x = np.arange(N_SPECIES)
    width = 0.2
    colors_diet = sns.color_palette("Set1", len(diet_names))
    
    for idx, (diet_name, res) in enumerate(results.items()):
        ax.bar(x + idx * width, res['growth_rates'], width, 
               label=diet_name, color=colors_diet[idx], alpha=0.85)
    
    ax.set_xticks(x + width * 1.5)
    ax.set_xticklabels(SPECIES, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Growth Rate (1/day)')
    ax.set_title('Species Growth Rates by Diet', fontweight='bold')
    ax.legend(fontsize=9)
    
    # Total SCFA by diet
    ax = axes[1]
    scfa_names = ['Acetate', 'Propionate', 'Butyrate']
    scfa_colors = ['#2ecc71', '#3498db', '#e74c3c']
    x = np.arange(len(diet_names))
    width = 0.25
    
    for i, (scfa_name, color) in enumerate(zip(scfa_names, scfa_colors)):
        vals = [results[d]['total_scfa'][i] for d in diet_names]
        ax.bar(x + i * width, vals, width, label=scfa_name, color=color, alpha=0.85)
    
    ax.set_xticks(x + width)
    ax.set_xticklabels(diet_names, rotation=15)
    ax.set_ylabel('Predicted SCFA Flux (mmol/day)')
    ax.set_title('Community SCFA Production', fontweight='bold')
    ax.legend()
    
    # Metabolic exchange network (heatmap)
    ax = axes[2]
    metabolite_labels = ['Glc↓', 'Fiber↓', 'Acetate↑', 'Prop↑', 'Butyrate↑', 
                         'Lactate↑', 'H₂↑', 'CO₂↑']
    
    hf_fluxes = results['High-Fiber']['metabolite_fluxes']
    sns.heatmap(hf_fluxes, ax=ax, cmap='YlOrRd', annot=True, fmt='.2f',
                xticklabels=metabolite_labels, yticklabels=SPECIES,
                cbar_kws={'label': 'Flux (mmol/day)'})
    ax.set_title('Metabolic Flux (High-Fiber)', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/{filename}', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  -> Saved {filename}")


# ============================================================
# 8. Interaction Network Heatmap
# ============================================================
def plot_interaction_matrix(A, filename="interaction_matrix.png"):
    """Plot the gLV interaction matrix."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    mask = np.zeros_like(A)
    # Only show significant interactions
    
    sns.heatmap(A, ax=ax, cmap='RdBu_r', center=0, annot=True, fmt='.4f',
                xticklabels=SPECIES, yticklabels=SPECIES,
                cbar_kws={'label': 'Interaction Coefficient'},
                linewidths=0.5)
    ax.set_title('gLV Species Interaction Matrix', fontsize=14, fontweight='bold')
    ax.set_xlabel('Effect of Species (column)')
    ax.set_ylabel('On Species (row)')
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/{filename}', dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  -> Saved {filename}")


# ============================================================
# 9. Summary Statistics
# ============================================================
def compute_summary_stats(results_diet, scenarios_prob, t_ferm, y_ferm, metab_results):
    """Compute and print summary statistics."""
    stats = {}
    
    # Diet comparison final abundances
    for diet_name, sol in results_diet.items():
        total = sol.y[:, -1].sum()
        rel = sol.y[:, -1] / total
        scfa = predict_scfa(sol.y[:, -1])
        
        # Shannon diversity
        rel_clip = np.clip(rel, 1e-10, 1)
        shannon = -np.sum(rel_clip * np.log(rel_clip))
        
        stats[diet_name] = {
            'dominant': SPECIES[np.argmax(rel)],
            'dominant_pct': np.max(rel) * 100,
            'shannon': shannon,
            'total_scfa': scfa.sum(),
            'butyrate': scfa[2],
            'acetate': scfa[0],
            'propionate': scfa[1],
        }
    
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    
    print("\n--- Dietary Pattern Effects (Day 90) ---")
    for diet_name, s in stats.items():
        print(f"\n  {diet_name}:")
        print(f"    Dominant species: {s['dominant']} ({s['dominant_pct']:.1f}%)")
        print(f"    Shannon Diversity: {s['shannon']:.3f}")
        print(f"    Total SCFA: {s['total_scfa']:.2f} mmol/L")
        print(f"    Butyrate: {s['butyrate']:.2f}, Acetate: {s['acetate']:.2f}, "
              f"Propionate: {s['propionate']:.2f}")
    
    # Probiotic/prebiotic comparison
    print("\n--- Probiotic/Prebiotic Intervention (Day 60) ---")
    for name, sol in scenarios_prob.items():
        total = sol.y[:, -1].sum()
        rel = sol.y[:, -1] / total
        rel_clip = np.clip(rel, 1e-10, 1)
        shannon = -np.sum(rel_clip * np.log(rel_clip))
        scfa = predict_scfa(sol.y[:, -1])
        print(f"  {name}: Shannon={shannon:.3f}, Butyrate={scfa[2]:.2f}")
    
    # Fermented food
    total_f = y_ferm[:, -1].sum()
    rel_f = y_ferm[:, -1] / total_f
    rel_f_clip = np.clip(rel_f, 1e-10, 1)
    shannon_f = -np.sum(rel_f_clip * np.log(rel_f_clip))
    
    total_f0 = y_ferm[:, 0].sum()
    rel_f0 = y_ferm[:, 0] / total_f0
    rel_f0_clip = np.clip(rel_f0, 1e-10, 1)
    shannon_f0 = -np.sum(rel_f0_clip * np.log(rel_f0_clip))
    
    print(f"\n--- Fermented Food Case Study ---")
    print(f"  Shannon Diversity: Baseline={shannon_f0:.3f} -> Post-intervention={shannon_f:.3f}")
    
    return stats


# ============================================================
# MAIN EXECUTION
# ============================================================
if __name__ == "__main__":
    print("="*80)
    print("Systems Biology Framework: Diet-Gut Microbiome Interactions")
    print("="*80)
    
    # Initial abundances (×10⁸ cells/mL)
    x0 = np.array([30, 20, 15, 10, 8, 12, 5, 7])
    
    # Build model parameters
    mu, A = build_glv_parameters()
    
    # --- 1. SHIME Digestion ---
    print("\n[1] SHIME-Inspired Digestion Kinetics...")
    diet_comp = [20, 15, 10, 8, 5]  # Standard diet
    t_dig, dig_results, compartments = shime_digestion(diet_comp)
    plot_shime_digestion(t_dig, dig_results, compartments)
    
    # --- 2. gLV Dynamics ---
    print("\n[2] gLV Community Dynamics...")
    sol_baseline = simulate_glv(x0, mu, A, (0, 60))
    plot_glv_dynamics(sol_baseline, "(Baseline)", "glv_dynamics.png")
    plot_interaction_matrix(A)
    
    # --- 3. SCFA Prediction ---
    print("\n[3] SCFA Flux Prediction...")
    scfa_baseline = plot_scfa_dynamics(sol_baseline, "scfa_dynamics.png")
    
    # --- 4. Long-term Dietary Patterns ---
    print("\n[4] Long-term Dietary Pattern Simulation...")
    results_diet = simulate_dietary_patterns(mu, A, x0, days=90)
    plot_dietary_patterns(results_diet)
    plot_diet_scfa_comparison(results_diet)
    
    # --- 5. Probiotic/Prebiotic ---
    print("\n[5] Probiotic/Prebiotic Effect Prediction...")
    scenarios_prob = simulate_probiotic_prebiotic(mu, A, x0, days=60)
    plot_probiotic_prebiotic(scenarios_prob)
    
    # --- 6. Fermented Food Case Study ---
    print("\n[6] Fermented Food Case Study...")
    t_ferm, y_ferm = simulate_fermented_food(mu, A, x0, days=90)
    plot_fermented_food(t_ferm, y_ferm)
    
    # --- 7. Community Metabolic Modeling ---
    print("\n[7] MICOM/gapseq-Inspired Community Metabolic Modeling...")
    metab_results = community_metabolic_model()
    plot_community_metabolic(metab_results)
    
    # --- Summary ---
    stats = compute_summary_stats(results_diet, scenarios_prob, t_ferm, y_ferm, metab_results)
    
    print("\n" + "="*80)
    print("All simulations complete. Figures saved to ./figures/")
    print("="*80)
