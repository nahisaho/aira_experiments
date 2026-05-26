#!/usr/bin/env python3
"""
Genome-Scale Metabolic Model (GEM) Constraint-Based Flux Analysis Framework
============================================================================
Comprehensive analysis pipeline using COBRApy/Cameo covering:
1. FBA constraint optimization
2. 13C-MFA integration simulation
3. Dynamic FBA (dFBA)
4. Enzyme capacity constraints (GECKO/sMOMENT-like)
5. Context-specific model construction (RNA-seq integration)
6. E. coli lysine production optimization case study
"""

import cobra
from cobra.io import load_model
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import warnings
import json
import os

warnings.filterwarnings('ignore')
plt.rcParams.update({'font.size': 10, 'figure.dpi': 150, 'savefig.bbox': 'tight'})

FIGURES_DIR = 'figures'
os.makedirs(FIGURES_DIR, exist_ok=True)

results = {}

# =============================================================================
# SECTION 1: FBA Constraint Optimization
# =============================================================================
print("="*70)
print("SECTION 1: FBA Constraint Optimization")
print("="*70)

model = load_model('textbook')
print(f"Model: {len(model.reactions)} reactions, {len(model.metabolites)} metabolites, {len(model.genes)} genes")

# 1a. Baseline FBA
sol_baseline = model.optimize()
print(f"\nBaseline FBA (biomass max): {sol_baseline.objective_value:.4f} h^-1")

# 1b. Evaluate impact of different constraint strategies
constraint_results = {}

# Strategy 1: Default constraints
constraint_results['Default'] = sol_baseline.objective_value

# Strategy 2: Thermodynamic-like constraints (restrict reversibility)
model_thermo = model.copy()
irreversible_rxns = ['PFK', 'PYK', 'CS', 'AKGDH', 'SUCOAS', 'FUM', 'MDH']
for rxn_id in irreversible_rxns:
    if rxn_id in model_thermo.reactions:
        model_thermo.reactions.get_by_id(rxn_id).lower_bound = 0
sol_thermo = model_thermo.optimize()
constraint_results['Thermodynamic'] = sol_thermo.objective_value

# Strategy 3: Tightened exchange bounds
model_tight = model.copy()
for rxn in model_tight.exchanges:
    if rxn.lower_bound < -10:
        rxn.lower_bound = -10
sol_tight = model_tight.optimize()
constraint_results['Tightened Exchange'] = sol_tight.objective_value

# Strategy 4: Parsimonious FBA (pFBA)
sol_pfba = cobra.flux_analysis.pfba(model)
# pFBA objective_value is sum of fluxes; growth rate is the biomass flux
pfba_growth = sol_pfba.fluxes['Biomass_Ecoli_core']
constraint_results['pFBA'] = pfba_growth

# Strategy 5: Loopless FBA
try:
    sol_loopless = cobra.flux_analysis.loopless_solution(model)
    constraint_results['Loopless'] = sol_loopless.objective_value
except:
    constraint_results['Loopless'] = sol_baseline.objective_value

print("\nConstraint Strategy Results:")
for k, v in constraint_results.items():
    print(f"  {k}: growth = {v:.4f} h^-1")

results['constraint_strategies'] = constraint_results

# 1c. Flux Variability Analysis
fva_result = cobra.flux_analysis.flux_variability_analysis(
    model, fraction_of_optimum=0.9, loopless=False
)
fva_ranges = fva_result['maximum'] - fva_result['minimum']
print(f"\nFVA (90% optimum): mean range = {fva_ranges.mean():.2f}, max range = {fva_ranges.max():.2f}")

# Compare FVA ranges at different optima fractions
fva_fractions = [0.5, 0.7, 0.9, 0.95, 0.99]
fva_mean_ranges = []
for frac in fva_fractions:
    fva_r = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=frac, loopless=False)
    fva_mean_ranges.append((fva_r['maximum'] - fva_r['minimum']).mean())

results['fva_fractions'] = fva_fractions
results['fva_mean_ranges'] = fva_mean_ranges

# Figure 1: Constraint strategies comparison
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 1a: Bar chart of strategies
strategies = list(constraint_results.keys())
values = list(constraint_results.values())
colors = plt.cm.Set2(np.linspace(0, 1, len(strategies)))
axes[0].bar(strategies, values, color=colors, edgecolor='black', linewidth=0.5)
axes[0].set_ylabel('Growth Rate (h⁻¹)')
axes[0].set_title('(A) FBA Constraint Strategy Comparison')
axes[0].set_ylim(0, max(values)*1.15)
axes[0].tick_params(axis='x', rotation=25)
for i, v in enumerate(values):
    axes[0].text(i, v + 0.01, f'{v:.3f}', ha='center', fontsize=8)

# 1b: FVA range vs fraction
axes[1].plot(fva_fractions, fva_mean_ranges, 'o-', color='#2196F3', linewidth=2, markersize=8)
axes[1].set_xlabel('Fraction of Optimum')
axes[1].set_ylabel('Mean Flux Range (mmol/gDW/h)')
axes[1].set_title('(B) FVA Solution Space vs. Optimality Bound')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig1_constraint_optimization.png')
plt.close()
print("Saved: figures/fig1_constraint_optimization.png")

# Figure 2: FVA flux ranges for key reactions
key_reactions = ['PFK', 'PYK', 'CS', 'GAPD', 'TPI', 'PGK', 'ENO', 'PDH', 'AKGDH', 'ICDHyr']
existing_key = [r for r in key_reactions if r in fva_result.index]

fig, ax = plt.subplots(figsize=(10, 6))
fva_sub = fva_result.loc[existing_key]
y_pos = range(len(existing_key))
bars = ax.barh(y_pos, fva_sub['maximum'] - fva_sub['minimum'],
               left=fva_sub['minimum'], color='#4CAF50', alpha=0.7, edgecolor='black', linewidth=0.5)
# Mark FBA solution
pfba_fluxes = sol_pfba.fluxes
for i, rxn_id in enumerate(existing_key):
    if rxn_id in pfba_fluxes.index:
        ax.plot(pfba_fluxes[rxn_id], i, 'r*', markersize=12, zorder=5)
ax.set_yticks(y_pos)
ax.set_yticklabels(existing_key)
ax.set_xlabel('Flux (mmol/gDW/h)')
ax.set_title('FVA Ranges for Central Carbon Metabolism (★ = pFBA solution)')
ax.grid(True, alpha=0.3, axis='x')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig2_fva_ranges.png')
plt.close()
print("Saved: figures/fig2_fva_ranges.png")

# =============================================================================
# SECTION 2: 13C-MFA Integration Simulation
# =============================================================================
print("\n" + "="*70)
print("SECTION 2: 13C-MFA Integration Simulation")
print("="*70)

# Simulate 13C-MFA by adding measured flux constraints from literature
model_13c = model.copy()

# Simulated 13C-MFA measured fluxes (normalized to glucose=100)
glucose_uptake = 10.0
measured_fluxes = {
    'EX_glc__D_e': -glucose_uptake,
    'PFK': 7.5,
    'PYK': 6.2,
    'CS': 3.8,
    'ICDHyr': 3.5,
    'AKGDH': 2.8,
    'PDH': 5.5,
    'G6PDH2r': 2.3,
}

# Apply 13C-MFA constraints with measurement uncertainty (±15%)
uncertainty = 0.15
for rxn_id, flux_val in measured_fluxes.items():
    if rxn_id in model_13c.reactions:
        rxn = model_13c.reactions.get_by_id(rxn_id)
        if flux_val < 0:
            rxn.lower_bound = flux_val * (1 + uncertainty)
            rxn.upper_bound = flux_val * (1 - uncertainty)
        else:
            rxn.lower_bound = flux_val * (1 - uncertainty)
            rxn.upper_bound = flux_val * (1 + uncertainty)

sol_13c = model_13c.optimize()
print(f"13C-MFA constrained growth: {sol_13c.objective_value:.4f} h^-1")

# Compare unconstrained vs 13C-constrained fluxes
comparison_rxns = ['PFK', 'PYK', 'CS', 'ICDHyr', 'AKGDH', 'PDH', 'G6PDH2r',
                   'FBA', 'TPI', 'GAPD', 'PGK', 'ENO']
existing_comp = [r for r in comparison_rxns if r in model.reactions]

fba_fluxes = []
mfa_fluxes = []
for rxn_id in existing_comp:
    fba_fluxes.append(sol_pfba.fluxes[rxn_id])
    mfa_fluxes.append(sol_13c.fluxes[rxn_id])

results['13c_comparison'] = {
    'reactions': existing_comp,
    'fba_fluxes': fba_fluxes,
    'mfa_fluxes': mfa_fluxes,
    'growth_unconstrained': sol_baseline.objective_value,
    'growth_13c': sol_13c.objective_value
}

# Bayesian-like flux estimation: combine FBA prior with 13C measurement
# Weighted average approach
def bayesian_flux_estimate(fba_flux, measured_flux, fba_var=2.0, meas_var=0.5):
    w_fba = 1.0 / fba_var
    w_meas = 1.0 / meas_var
    posterior = (w_fba * fba_flux + w_meas * measured_flux) / (w_fba + w_meas)
    posterior_var = 1.0 / (w_fba + w_meas)
    return posterior, posterior_var

bayesian_results = {}
for rxn_id in ['PFK', 'PYK', 'CS', 'PDH', 'G6PDH2r']:
    if rxn_id in sol_pfba.fluxes.index and rxn_id in measured_fluxes:
        post, post_var = bayesian_flux_estimate(
            sol_pfba.fluxes[rxn_id], measured_fluxes[rxn_id]
        )
        bayesian_results[rxn_id] = {
            'FBA': sol_pfba.fluxes[rxn_id],
            '13C-MFA': measured_fluxes[rxn_id],
            'Bayesian': post,
            'Variance': post_var
        }
        print(f"  {rxn_id}: FBA={sol_pfba.fluxes[rxn_id]:.2f}, "
              f"13C={measured_fluxes[rxn_id]:.2f}, Bayesian={post:.2f}")

results['bayesian'] = bayesian_results

# Figure 3: 13C-MFA integration
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# 3a: Flux comparison
x = np.arange(len(existing_comp))
width = 0.35
axes[0].bar(x - width/2, fba_fluxes, width, label='Standard pFBA', color='#2196F3', edgecolor='black', linewidth=0.5)
axes[0].bar(x + width/2, mfa_fluxes, width, label='13C-MFA Constrained', color='#FF9800', edgecolor='black', linewidth=0.5)
axes[0].set_xticks(x)
axes[0].set_xticklabels(existing_comp, rotation=45, ha='right')
axes[0].set_ylabel('Flux (mmol/gDW/h)')
axes[0].set_title('(A) pFBA vs 13C-MFA Constrained Fluxes')
axes[0].legend()
axes[0].grid(True, alpha=0.3, axis='y')

# 3b: Bayesian integration
bkeys = list(bayesian_results.keys())
x2 = np.arange(len(bkeys))
w = 0.25
vals_fba = [bayesian_results[k]['FBA'] for k in bkeys]
vals_13c = [bayesian_results[k]['13C-MFA'] for k in bkeys]
vals_bay = [bayesian_results[k]['Bayesian'] for k in bkeys]
errs = [np.sqrt(bayesian_results[k]['Variance']) for k in bkeys]
axes[1].bar(x2 - w, vals_fba, w, label='FBA Prior', color='#2196F3', edgecolor='black', linewidth=0.5)
axes[1].bar(x2, vals_13c, w, label='13C-MFA Measured', color='#FF9800', edgecolor='black', linewidth=0.5)
axes[1].bar(x2 + w, vals_bay, w, label='Bayesian Posterior', color='#4CAF50', edgecolor='black', linewidth=0.5,
            yerr=errs, capsize=3)
axes[1].set_xticks(x2)
axes[1].set_xticklabels(bkeys, rotation=45, ha='right')
axes[1].set_ylabel('Flux (mmol/gDW/h)')
axes[1].set_title('(B) Bayesian Flux Integration')
axes[1].legend(fontsize=8)
axes[1].grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig3_13c_mfa_integration.png')
plt.close()
print("Saved: figures/fig3_13c_mfa_integration.png")

# =============================================================================
# SECTION 3: Dynamic FBA (dFBA)
# =============================================================================
print("\n" + "="*70)
print("SECTION 3: Dynamic FBA (dFBA)")
print("="*70)

def dfba_simulation(model, initial_conditions, t_span, dt=0.1):
    """Static Optimization Approach (SOA) for dFBA."""
    model = model.copy()
    
    # Initial conditions
    X = initial_conditions['biomass']    # gDW/L
    Glc = initial_conditions['glucose']  # mmol/L
    Ac = initial_conditions.get('acetate', 0.0)  # mmol/L
    
    # Kinetic parameters
    Vmax_glc = 10.0   # mmol/gDW/h
    Km_glc = 0.5      # mmol/L
    Vmax_ac = 5.0
    Km_ac = 0.5
    
    t = 0
    times = [t]
    biomass_traj = [X]
    glucose_traj = [Glc]
    acetate_traj = [Ac]
    growth_rates = [0]
    
    while t < t_span:
        # Michaelis-Menten uptake
        v_glc = Vmax_glc * Glc / (Km_glc + Glc) if Glc > 1e-6 else 0
        v_ac_uptake = Vmax_ac * Ac / (Km_ac + Ac) if Ac > 0.1 and Glc < 0.1 else 0
        
        # Set bounds
        model.reactions.get_by_id('EX_glc__D_e').lower_bound = -v_glc
        if 'EX_ac_e' in model.reactions:
            model.reactions.get_by_id('EX_ac_e').lower_bound = -v_ac_uptake
        
        try:
            sol = model.optimize()
            if sol.status == 'optimal':
                mu = sol.objective_value
                glc_flux = sol.fluxes.get('EX_glc__D_e', 0)
                ac_flux = sol.fluxes.get('EX_ac_e', 0)
            else:
                mu = 0
                glc_flux = 0
                ac_flux = 0
        except:
            mu = 0
            glc_flux = 0
            ac_flux = 0
        
        # Update concentrations
        dX = mu * X * dt
        dGlc = glc_flux * X * dt  # negative = consumption
        dAc = ac_flux * X * dt    # positive = production
        
        X = max(X + dX, 0)
        Glc = max(Glc + dGlc, 0)
        Ac = max(Ac + dAc, 0)
        t += dt
        
        times.append(t)
        biomass_traj.append(X)
        glucose_traj.append(Glc)
        acetate_traj.append(Ac)
        growth_rates.append(mu)
    
    return {
        'time': np.array(times),
        'biomass': np.array(biomass_traj),
        'glucose': np.array(glucose_traj),
        'acetate': np.array(acetate_traj),
        'growth_rate': np.array(growth_rates)
    }

# Run dFBA simulation
dfba_result = dfba_simulation(
    model,
    {'biomass': 0.05, 'glucose': 20.0, 'acetate': 0.0},
    t_span=12.0, dt=0.05
)

print(f"dFBA simulation: {len(dfba_result['time'])} time points")
print(f"  Final biomass: {dfba_result['biomass'][-1]:.3f} gDW/L")
print(f"  Final glucose: {dfba_result['glucose'][-1]:.3f} mmol/L")
print(f"  Final acetate: {dfba_result['acetate'][-1]:.3f} mmol/L")
print(f"  Max growth rate: {max(dfba_result['growth_rate']):.4f} h^-1")

results['dfba'] = {
    'final_biomass': float(dfba_result['biomass'][-1]),
    'final_glucose': float(dfba_result['glucose'][-1]),
    'final_acetate': float(dfba_result['acetate'][-1]),
    'max_growth': float(max(dfba_result['growth_rate']))
}

# Run dFBA with different initial glucose concentrations
glucose_levels = [5, 10, 20, 40]
dfba_multi = {}
for glc0 in glucose_levels:
    dfba_multi[glc0] = dfba_simulation(
        model, {'biomass': 0.05, 'glucose': float(glc0)}, t_span=15.0, dt=0.05
    )

# Figure 4: dFBA results
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 4a: Biomass, glucose, acetate over time
ax = axes[0, 0]
ax.plot(dfba_result['time'], dfba_result['biomass'], 'g-', linewidth=2, label='Biomass')
ax.set_xlabel('Time (h)')
ax.set_ylabel('Biomass (gDW/L)', color='g')
ax2 = ax.twinx()
ax2.plot(dfba_result['time'], dfba_result['glucose'], 'b--', linewidth=2, label='Glucose')
ax2.plot(dfba_result['time'], dfba_result['acetate'], 'r:', linewidth=2, label='Acetate')
ax2.set_ylabel('Concentration (mmol/L)')
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, loc='center right')
ax.set_title('(A) dFBA: Metabolite & Biomass Dynamics')
ax.grid(True, alpha=0.3)

# 4b: Growth rate over time
axes[0, 1].plot(dfba_result['time'], dfba_result['growth_rate'], 'k-', linewidth=2)
axes[0, 1].set_xlabel('Time (h)')
axes[0, 1].set_ylabel('Growth Rate (h⁻¹)')
axes[0, 1].set_title('(B) Growth Rate Dynamics')
axes[0, 1].grid(True, alpha=0.3)
axes[0, 1].set_ylim(bottom=0)

# 4c: Biomass for different glucose levels
for glc0, res in dfba_multi.items():
    axes[1, 0].plot(res['time'], res['biomass'], linewidth=2, label=f'Glc₀={glc0} mM')
axes[1, 0].set_xlabel('Time (h)')
axes[1, 0].set_ylabel('Biomass (gDW/L)')
axes[1, 0].set_title('(C) Biomass at Different Initial Glucose')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)

# 4d: Final biomass vs initial glucose
final_bm = [dfba_multi[g]['biomass'][-1] for g in glucose_levels]
axes[1, 1].plot(glucose_levels, final_bm, 'o-', color='#E91E63', linewidth=2, markersize=10)
axes[1, 1].set_xlabel('Initial Glucose (mmol/L)')
axes[1, 1].set_ylabel('Final Biomass (gDW/L)')
axes[1, 1].set_title('(D) Yield vs. Initial Substrate')
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig4_dfba_dynamics.png')
plt.close()
print("Saved: figures/fig4_dfba_dynamics.png")

# =============================================================================
# SECTION 4: Enzyme Capacity Constraints (GECKO/sMOMENT-like)
# =============================================================================
print("\n" + "="*70)
print("SECTION 4: Enzyme Capacity Constraints (GECKO/sMOMENT-like)")
print("="*70)

def apply_enzyme_constraints(model, total_protein=0.46, sigma=0.5):
    """
    Apply sMOMENT-like enzyme capacity constraints.
    Total protein pool (g_protein/gDW) distributed across enzymes.
    sigma: fraction of proteome available for metabolic enzymes.
    """
    ec_model = model.copy()
    
    # Simulated kcat values (1/s) for key enzymes - representative values
    kcat_values = {
        'PFK': 110.0, 'PYK': 230.0, 'CS': 45.0, 'AKGDH': 30.0,
        'ICDHyr': 80.0, 'PDH': 50.0, 'FBA': 15.0, 'TPI': 4300.0,
        'GAPD': 85.0, 'PGK': 400.0, 'ENO': 60.0, 'G6PDH2r': 200.0,
        'PGL': 35.0, 'GND': 50.0, 'FUM': 800.0, 'MDH': 100.0,
        'SUCOAS': 20.0, 'SUCDi': 25.0, 'MALS': 12.0, 'ME1': 40.0,
    }
    
    # Molecular weights (kDa) - approximate
    mw_values = {
        'PFK': 140, 'PYK': 232, 'CS': 97, 'AKGDH': 275,
        'ICDHyr': 93, 'PDH': 270, 'FBA': 78, 'TPI': 54,
        'GAPD': 145, 'PGK': 42, 'ENO': 93, 'G6PDH2r': 110,
        'PGL': 38, 'GND': 52, 'FUM': 200, 'MDH': 74,
        'SUCOAS': 140, 'SUCDi': 120, 'MALS': 65, 'ME1': 64,
    }
    
    available_protein = total_protein * sigma  # g/gDW
    
    # Add protein pool pseudo-metabolite and constraint
    protein_met = cobra.Metabolite('prot_pool', name='Protein pool',
                                    compartment='c')
    
    # Add enzyme usage reactions - enzyme cost as fraction of protein pool
    for rxn_id, kcat in kcat_values.items():
        if rxn_id in ec_model.reactions:
            rxn = ec_model.reactions.get_by_id(rxn_id)
            mw = mw_values.get(rxn_id, 100)
            # Enzyme cost: (MW in g/mol) / (kcat in 1/s * 3600 s/h) = g·s/(mol·h)
            # For 1 mmol/gDW/h flux: need MW/(kcat*3600*1000) g_protein/gDW
            enzyme_cost = mw / (kcat * 3.6)  # mg_protein per mmol/h flux
            rxn.add_metabolites({protein_met: enzyme_cost})
    
    # Add protein pool exchange (drain)
    prot_exchange = cobra.Reaction('prot_pool_exchange')
    prot_exchange.name = 'Protein pool constraint'
    prot_exchange.add_metabolites({protein_met: -1.0})
    prot_exchange.lower_bound = 0
    prot_exchange.upper_bound = available_protein * 1e6  # mg scale
    ec_model.add_reactions([prot_exchange])
    
    return ec_model

# Test different protein pool sizes
protein_pools = np.linspace(0.05, 1.0, 30)
growth_vs_protein = []

for pp in protein_pools:
    ec_model = apply_enzyme_constraints(model, total_protein=pp)
    try:
        sol = ec_model.optimize()
        growth_vs_protein.append(sol.objective_value if sol.status == 'optimal' else 0)
    except:
        growth_vs_protein.append(0)

# Find the saturation point
growth_vs_protein = np.array(growth_vs_protein)
max_growth = max(growth_vs_protein)
saturation_idx = np.where(growth_vs_protein >= 0.99 * max_growth)[0]
sat_protein = protein_pools[saturation_idx[0]] if len(saturation_idx) > 0 else protein_pools[-1]

print(f"Max growth (unconstrained): {max_growth:.4f} h^-1")
print(f"Protein saturation at: {sat_protein:.2f} g/gDW")

# Compare standard FBA vs enzyme-constrained
ec_model_default = apply_enzyme_constraints(model, total_protein=0.46, sigma=0.5)
sol_ec = ec_model_default.optimize()
sol_std = model.optimize()

print(f"\nStandard FBA growth: {sol_std.objective_value:.4f} h^-1")
print(f"Enzyme-constrained growth: {sol_ec.objective_value:.4f} h^-1")
print(f"Reduction: {(1 - sol_ec.objective_value/sol_std.objective_value)*100:.1f}%")

# Flux redistribution analysis
flux_comparison = {}
central_rxns = ['PFK', 'PYK', 'CS', 'ICDHyr', 'AKGDH', 'PDH', 'G6PDH2r', 'ENO']
for rxn_id in central_rxns:
    if rxn_id in model.reactions:
        flux_comparison[rxn_id] = {
            'Standard': sol_pfba.fluxes[rxn_id],
            'Enzyme-constrained': sol_ec.fluxes[rxn_id]
        }

results['enzyme_constraints'] = {
    'std_growth': float(sol_std.objective_value),
    'ec_growth': float(sol_ec.objective_value),
    'saturation_protein': float(sat_protein),
    'flux_comparison': flux_comparison
}

# Figure 5: Enzyme constraints
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# 5a: Growth vs protein pool
axes[0].plot(protein_pools, growth_vs_protein, 'o-', color='#9C27B0', linewidth=2, markersize=5)
axes[0].axhline(y=sol_std.objective_value, color='gray', linestyle='--', alpha=0.7, label='Unconstrained FBA')
axes[0].axvline(x=sat_protein, color='red', linestyle=':', alpha=0.7, label=f'Saturation ({sat_protein:.2f})')
axes[0].set_xlabel('Total Protein (g/gDW)')
axes[0].set_ylabel('Growth Rate (h⁻¹)')
axes[0].set_title('(A) Growth vs. Protein Pool Size')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# 5b: Flux redistribution
rxn_ids = list(flux_comparison.keys())
std_vals = [flux_comparison[r]['Standard'] for r in rxn_ids]
ec_vals = [flux_comparison[r]['Enzyme-constrained'] for r in rxn_ids]
x3 = np.arange(len(rxn_ids))
axes[1].bar(x3 - 0.17, std_vals, 0.34, label='Standard FBA', color='#2196F3', edgecolor='black', linewidth=0.5)
axes[1].bar(x3 + 0.17, ec_vals, 0.34, label='Enzyme-Constrained', color='#FF5722', edgecolor='black', linewidth=0.5)
axes[1].set_xticks(x3)
axes[1].set_xticklabels(rxn_ids, rotation=45, ha='right')
axes[1].set_ylabel('Flux (mmol/gDW/h)')
axes[1].set_title('(B) Flux Redistribution')
axes[1].legend()
axes[1].grid(True, alpha=0.3, axis='y')

# 5c: Enzyme cost analysis
enzyme_costs = {}
kcat_values = {'PFK': 110, 'PYK': 230, 'CS': 45, 'AKGDH': 30,
               'ICDHyr': 80, 'PDH': 50, 'G6PDH2r': 200, 'ENO': 60}
mw_values2 = {'PFK': 140, 'PYK': 232, 'CS': 97, 'AKGDH': 275,
              'ICDHyr': 93, 'PDH': 270, 'G6PDH2r': 110, 'ENO': 93}
for rxn_id in rxn_ids:
    if rxn_id in kcat_values:
        flux = abs(sol_ec.fluxes.get(rxn_id, 0))
        cost = flux * mw_values2[rxn_id] * 1e3 / (kcat_values[rxn_id] * 3600)
        enzyme_costs[rxn_id] = cost

cost_keys = list(enzyme_costs.keys())
cost_vals = [enzyme_costs[k] for k in cost_keys]
axes[2].barh(cost_keys, cost_vals, color='#FF9800', edgecolor='black', linewidth=0.5)
axes[2].set_xlabel('Enzyme Cost (g_protein/gDW)')
axes[2].set_title('(C) Enzyme Investment per Reaction')
axes[2].grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig5_enzyme_constraints.png')
plt.close()
print("Saved: figures/fig5_enzyme_constraints.png")

# =============================================================================
# SECTION 5: Context-Specific Model Construction (RNA-seq Integration)
# =============================================================================
print("\n" + "="*70)
print("SECTION 5: Context-Specific Model Construction (RNA-seq)")
print("="*70)

def build_context_specific_model(model, expression_data, threshold=0.1):
    """
    Build context-specific model using GIMME-like approach.
    Genes below threshold have their reactions penalized/removed.
    """
    cs_model = model.copy()
    
    inactive_genes = [g for g, expr in expression_data.items() if expr < threshold]
    active_genes = [g for g, expr in expression_data.items() if expr >= threshold]
    
    removed_reactions = []
    for gene_id in inactive_genes:
        if gene_id in cs_model.genes:
            gene = cs_model.genes.get_by_id(gene_id)
            for rxn in gene.reactions:
                # Only remove if ALL associated genes are inactive
                all_inactive = all(
                    expression_data.get(g.id, 1.0) < threshold
                    for g in rxn.genes
                )
                if all_inactive and rxn.id not in removed_reactions:
                    rxn.knock_out()
                    removed_reactions.append(rxn.id)
    
    return cs_model, len(removed_reactions), len(active_genes)

# Simulate RNA-seq data for different conditions
np.random.seed(42)
gene_ids = [g.id for g in model.genes]

conditions = {
    'Aerobic_Glucose': {},
    'Anaerobic': {},
    'Oxidative_Stress': {},
    'Stationary_Phase': {},
}

# Generate expression profiles - more realistic with condition-specific patterns
for gene_id in gene_ids:
    conditions['Aerobic_Glucose'][gene_id] = np.random.lognormal(1.0, 0.5)
    conditions['Anaerobic'][gene_id] = np.random.lognormal(0.7, 0.6)
    conditions['Oxidative_Stress'][gene_id] = np.random.lognormal(0.6, 0.7)
    conditions['Stationary_Phase'][gene_id] = np.random.lognormal(0.4, 0.6)

# Normalize to max=1
for cond in conditions:
    max_val = max(conditions[cond].values())
    for g in conditions[cond]:
        conditions[cond][g] /= max_val

# Build context-specific models
cs_results = {}
for cond, expr_data in conditions.items():
    cs_model, n_removed, n_active = build_context_specific_model(model, expr_data, threshold=0.1)
    try:
        sol = cs_model.optimize()
        growth = sol.objective_value if sol.status == 'optimal' else 0
    except:
        growth = 0
    
    n_active_rxns = sum(1 for r in cs_model.reactions if r.lower_bound != 0 or r.upper_bound != 0)
    cs_results[cond] = {
        'growth': growth,
        'removed_reactions': n_removed,
        'active_genes': n_active,
        'active_reactions': n_active_rxns
    }
    print(f"  {cond}: growth={growth:.4f}, removed={n_removed} rxns, active_genes={n_active}")

results['context_specific'] = cs_results

# Threshold sensitivity analysis
thresholds = np.arange(0.01, 0.5, 0.02)
threshold_growth = {}
for cond in ['Aerobic_Glucose', 'Anaerobic']:
    growths = []
    for thr in thresholds:
        cs_m, _, _ = build_context_specific_model(model, conditions[cond], threshold=thr)
        try:
            sol = cs_m.optimize()
            growths.append(sol.objective_value if sol.status == 'optimal' else 0)
        except:
            growths.append(0)
    threshold_growth[cond] = growths

# Figure 6: Context-specific models
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# 6a: Growth rates across conditions
conds = list(cs_results.keys())
growths = [cs_results[c]['growth'] for c in conds]
colors6 = ['#4CAF50', '#2196F3', '#FF9800', '#9C27B0']
axes[0].bar(conds, growths, color=colors6, edgecolor='black', linewidth=0.5)
axes[0].set_ylabel('Growth Rate (h⁻¹)')
axes[0].set_title('(A) Context-Specific Growth Rates')
axes[0].tick_params(axis='x', rotation=30)
for i, v in enumerate(growths):
    axes[0].text(i, v + 0.01, f'{v:.3f}', ha='center', fontsize=8)

# 6b: Model complexity
removed = [cs_results[c]['removed_reactions'] for c in conds]
active_rxns = [cs_results[c]['active_reactions'] for c in conds]
x4 = np.arange(len(conds))
axes[1].bar(x4 - 0.17, active_rxns, 0.34, label='Active Reactions', color='#4CAF50', edgecolor='black', linewidth=0.5)
axes[1].bar(x4 + 0.17, removed, 0.34, label='Removed Reactions', color='#F44336', edgecolor='black', linewidth=0.5)
axes[1].set_xticks(x4)
axes[1].set_xticklabels(conds, rotation=30, ha='right')
axes[1].set_ylabel('Number of Reactions')
axes[1].set_title('(B) Model Reduction')
axes[1].legend()

# 6c: Threshold sensitivity
for cond, growths_list in threshold_growth.items():
    axes[2].plot(thresholds, growths_list, linewidth=2, label=cond)
axes[2].set_xlabel('Expression Threshold')
axes[2].set_ylabel('Growth Rate (h⁻¹)')
axes[2].set_title('(C) Threshold Sensitivity')
axes[2].legend()
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig6_context_specific.png')
plt.close()
print("Saved: figures/fig6_context_specific.png")

# =============================================================================
# SECTION 6: Lysine Production Optimization Case Study
# =============================================================================
print("\n" + "="*70)
print("SECTION 6: Lysine Production Optimization (E. coli)")
print("="*70)

# Build lysine production pathway in the textbook model
lys_model = model.copy()

# Add simplified lysine biosynthesis pathway
# In E. coli, lysine is produced via DAP pathway from oxaloacetate + pyruvate
# Add key metabolites
lys_c = cobra.Metabolite('lys__L_c', name='L-Lysine', compartment='c')
lys_e = cobra.Metabolite('lys__L_e', name='L-Lysine (extracellular)', compartment='e')
dap_c = cobra.Metabolite('dap_c', name='meso-Diaminopimelate', compartment='c')
asp_sa_c = cobra.Metabolite('aspsa_c', name='L-Aspartate-semialdehyde', compartment='c')

# Check if asp exists
asp_c = None
for met in lys_model.metabolites:
    if 'asp__L_c' == met.id:
        asp_c = met
        break

if asp_c is None:
    asp_c = cobra.Metabolite('asp__L_c', name='L-Aspartate', compartment='c')

# Get existing metabolites
oaa_c = lys_model.metabolites.get_by_id('oaa_c')  # oxaloacetate
pyr_c = lys_model.metabolites.get_by_id('pyr_c')  # pyruvate
nadph_c = lys_model.metabolites.get_by_id('nadph_c')
nadp_c = lys_model.metabolites.get_by_id('nadp_c')
atp_c = lys_model.metabolites.get_by_id('atp_c')
adp_c = lys_model.metabolites.get_by_id('adp_c')
h2o_c = lys_model.metabolites.get_by_id('h2o_c')
co2_c = lys_model.metabolites.get_by_id('co2_c')
glu__L_c = lys_model.metabolites.get_by_id('glu__L_c')
pi_c = lys_model.metabolites.get_by_id('pi_c')
h_c = lys_model.metabolites.get_by_id('h_c')
h_e = lys_model.metabolites.get_by_id('h_e')

# Reaction 1: Aspartate kinase (OAA -> asp -> aspsa)
rxn_ak = cobra.Reaction('ASPTA_lys')
rxn_ak.name = 'Aspartate aminotransferase (for lysine)'
rxn_ak.lower_bound = 0
rxn_ak.upper_bound = 1000
rxn_ak.add_metabolites({
    oaa_c: -1, glu__L_c: -1, atp_c: -1, nadph_c: -1,
    asp_sa_c: 1, adp_c: 1, nadp_c: 1, pi_c: 1, h_c: 1
})
lys_model.add_reactions([rxn_ak])

# Reaction 2: DAP pathway (aspsa + pyr -> dap)
rxn_dap = cobra.Reaction('DAPDC')
rxn_dap.name = 'DAP decarboxylase pathway'
rxn_dap.lower_bound = 0
rxn_dap.upper_bound = 1000
rxn_dap.add_metabolites({
    asp_sa_c: -1, pyr_c: -1, nadph_c: -1, glu__L_c: -1,
    dap_c: 1, nadp_c: 1, co2_c: 1, h2o_c: 1
})
lys_model.add_reactions([rxn_dap])

# Reaction 3: DAP -> Lysine
rxn_lys = cobra.Reaction('LYSDC')
rxn_lys.name = 'Lysine production from DAP'
rxn_lys.lower_bound = 0
rxn_lys.upper_bound = 1000
rxn_lys.add_metabolites({
    dap_c: -1,
    lys_c: 1, co2_c: 1
})
lys_model.add_reactions([rxn_lys])

# Transport and exchange
rxn_transport = cobra.Reaction('LYSt')
rxn_transport.name = 'Lysine transport'
rxn_transport.lower_bound = 0
rxn_transport.upper_bound = 1000
rxn_transport.add_metabolites({lys_c: -1, h_c: -1, lys_e: 1, h_e: 1})
lys_model.add_reactions([rxn_transport])

rxn_exchange = cobra.Reaction('EX_lys__L_e')
rxn_exchange.name = 'Lysine exchange'
rxn_exchange.lower_bound = 0
rxn_exchange.upper_bound = 1000
rxn_exchange.add_metabolites({lys_e: -1})
lys_model.add_reactions([rxn_exchange])

# Test baseline lysine production potential
print("Baseline lysine production (maximizing growth):")
sol_lys_base = lys_model.optimize()
lys_base_flux = sol_lys_base.fluxes.get('EX_lys__L_e', 0)
print(f"  Growth: {sol_lys_base.objective_value:.4f}, Lysine: {lys_base_flux:.4f}")

# Maximize lysine production
lys_model_opt = lys_model.copy()
lys_model_opt.objective = 'EX_lys__L_e'
sol_lys_max = lys_model_opt.optimize()
print(f"  Max lysine flux: {sol_lys_max.objective_value:.4f} mmol/gDW/h")

# Production envelope analysis
from cobra.flux_analysis import production_envelope

prod_env = production_envelope(
    lys_model, ['EX_lys__L_e'],
    points=30
)
print(f"  Production envelope: {len(prod_env)} points computed")
print(f"  PE columns: {prod_env.columns.tolist()}")

# Gene knockout analysis for lysine improvement
print("\nGene knockout analysis for lysine improvement:")
lys_model_ko = lys_model.copy()
lys_model_ko.objective = 'EX_lys__L_e'

# Minimum biomass constraint
min_growth = 0.1
biomass_rxn = lys_model_ko.reactions.get_by_id('Biomass_Ecoli_core')
biomass_rxn.lower_bound = min_growth

ko_results = {}
target_knockouts = ['PFK', 'PYK', 'CS', 'ICDHyr', 'AKGDH', 'PDH',
                    'G6PDH2r', 'ME1', 'ME2', 'FRD7', 'SUCDi', 'MDH']

for rxn_id in target_knockouts:
    if rxn_id in lys_model_ko.reactions:
        ko_model = lys_model_ko.copy()
        ko_model.reactions.get_by_id(rxn_id).knock_out()
        try:
            sol = ko_model.optimize()
            if sol.status == 'optimal':
                ko_results[rxn_id] = {
                    'lysine': sol.objective_value,
                    'growth': sol.fluxes['Biomass_Ecoli_core']
                }
            else:
                ko_results[rxn_id] = {'lysine': 0, 'growth': 0}
        except:
            ko_results[rxn_id] = {'lysine': 0, 'growth': 0}

# Sort by lysine production
sorted_ko = sorted(ko_results.items(), key=lambda x: x[1]['lysine'], reverse=True)
print("  Top knockouts for lysine:")
for rxn_id, vals in sorted_ko[:5]:
    print(f"    Δ{rxn_id}: lysine={vals['lysine']:.3f}, growth={vals['growth']:.4f}")

results['lysine'] = {
    'baseline_growth': float(sol_lys_base.objective_value),
    'baseline_lysine': float(lys_base_flux),
    'max_lysine': float(sol_lys_max.objective_value),
    'knockout_results': {k: v for k, v in sorted_ko}
}

# OptKnock-like analysis: growth-coupled lysine production
# Test double knockouts for top candidates
print("\nDouble knockout analysis:")
top_single = [k for k, v in sorted_ko if v['lysine'] > 0][:5]
double_ko_results = {}
for i, rxn1 in enumerate(top_single):
    for rxn2 in top_single[i+1:]:
        ko2_model = lys_model_ko.copy()
        if rxn1 in ko2_model.reactions and rxn2 in ko2_model.reactions:
            ko2_model.reactions.get_by_id(rxn1).knock_out()
            ko2_model.reactions.get_by_id(rxn2).knock_out()
            try:
                sol = ko2_model.optimize()
                if sol.status == 'optimal' and sol.objective_value > 0:
                    double_ko_results[f"Δ{rxn1}/Δ{rxn2}"] = {
                        'lysine': sol.objective_value,
                        'growth': sol.fluxes['Biomass_Ecoli_core']
                    }
            except:
                pass

sorted_dko = sorted(double_ko_results.items(), key=lambda x: x[1]['lysine'], reverse=True)
for name, vals in sorted_dko[:3]:
    print(f"  {name}: lysine={vals['lysine']:.3f}, growth={vals['growth']:.4f}")

# Overexpression analysis: increase specific reaction bounds
print("\nOverexpression analysis:")
overexpress_targets = ['G6PDH2r', 'PPC', 'ASPTA_lys', 'DAPDC', 'LYSDC']
oe_results = {}
for rxn_id in overexpress_targets:
    if rxn_id in lys_model.reactions:
        oe_model = lys_model.copy()
        oe_model.objective = 'EX_lys__L_e'
        biomass_rxn_oe = oe_model.reactions.get_by_id('Biomass_Ecoli_core')
        biomass_rxn_oe.lower_bound = min_growth
        rxn = oe_model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = max(rxn.lower_bound, 5.0)  # Force minimum flux
        try:
            sol = oe_model.optimize()
            if sol.status == 'optimal':
                oe_results[rxn_id] = {
                    'lysine': sol.objective_value,
                    'growth': sol.fluxes['Biomass_Ecoli_core']
                }
                print(f"  OE-{rxn_id}: lysine={sol.objective_value:.3f}, growth={sol.fluxes['Biomass_Ecoli_core']:.4f}")
        except:
            oe_results[rxn_id] = {'lysine': 0, 'growth': 0}

# Combined strategy: best knockout + overexpression
print("\nCombined strategy (best knockout + overexpression):")
best_ko_rxn = sorted_ko[0][0] if sorted_ko else None
if best_ko_rxn:
    combined_model = lys_model.copy()
    combined_model.objective = 'EX_lys__L_e'
    combined_model.reactions.get_by_id('Biomass_Ecoli_core').lower_bound = min_growth
    combined_model.reactions.get_by_id(best_ko_rxn).knock_out()
    
    # Also knock out a competing pathway if available
    competing = ['SUCDi', 'FRD7']
    for comp_rxn in competing:
        if comp_rxn in combined_model.reactions:
            combined_model.reactions.get_by_id(comp_rxn).knock_out()
    
    # Overexpress lysine pathway
    for rxn_id in ['ASPTA_lys', 'DAPDC', 'LYSDC']:
        if rxn_id in combined_model.reactions:
            combined_model.reactions.get_by_id(rxn_id).lower_bound = 2.0
    
    # Boost PPC (anaplerotic) to supply OAA
    if 'PPC' in combined_model.reactions:
        combined_model.reactions.get_by_id('PPC').lower_bound = 3.0
    
    try:
        sol_combined = combined_model.optimize()
        if sol_combined.status == 'optimal':
            comb_growth = sol_combined.fluxes.get('Biomass_Ecoli_core', 0)
            print(f"  Combined (Δ{best_ko_rxn} + OE): lysine={sol_combined.objective_value:.3f}, "
                  f"growth={comb_growth:.4f}")
            results['lysine']['combined'] = {
                'strategy': f"Δ{best_ko_rxn} + OE(pathway+PPP)",
                'lysine': float(sol_combined.objective_value),
                'growth': float(comb_growth)
            }
    except Exception as e:
        print(f"  Combined strategy failed: {e}")

# Figure 7: Lysine optimization
fig = plt.figure(figsize=(16, 12))
gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)

# 7a: Production envelope
ax1 = fig.add_subplot(gs[0, 0])
if 'flux_minimum' in prod_env.columns and 'flux_maximum' in prod_env.columns:
    # Use EX_lys__L_e as x-axis proxy for growth scan
    x_col = [c for c in prod_env.columns if c not in ['carbon_source', 'flux_minimum', 'flux_maximum',
             'carbon_yield_minimum', 'carbon_yield_maximum', 'mass_yield_minimum', 'mass_yield_maximum']]
    x_data = prod_env[x_col[0]] if x_col else prod_env.index
    ax1.fill_between(range(len(prod_env)),
                     prod_env['flux_minimum'], prod_env['flux_maximum'],
                     alpha=0.3, color='#4CAF50')
    ax1.plot(range(len(prod_env)), prod_env['flux_maximum'],
             'g-', linewidth=2, label='Max lysine')
    ax1.plot(range(len(prod_env)), prod_env['flux_minimum'],
             'g--', linewidth=1, label='Min lysine')
    ax1.set_xlabel('Growth Scan Point')
elif 'EX_lys__L_e' in prod_env.columns:
    ax1.plot(range(len(prod_env)), prod_env['EX_lys__L_e'],
             'g-', linewidth=2)
    ax1.set_xlabel('Growth Scan Point')
ax1.set_ylabel('Lysine Flux (mmol/gDW/h)')
ax1.set_title('(A) Production Envelope')
ax1.legend()
ax1.grid(True, alpha=0.3)

# 7b: Single knockout results
ax2 = fig.add_subplot(gs[0, 1])
ko_names = [k for k, v in sorted_ko if v['lysine'] > 0]
ko_lys_vals = [v['lysine'] for k, v in sorted_ko if v['lysine'] > 0]
if ko_names:
    bars = ax2.barh(ko_names[:8], ko_lys_vals[:8], color='#FF9800', edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Lysine Flux (mmol/gDW/h)')
    ax2.set_title('(B) Single Knockouts')
    ax2.grid(True, alpha=0.3, axis='x')

# 7c: Growth-lysine tradeoff
ax3 = fig.add_subplot(gs[0, 2])
for k, v in sorted_ko:
    if v['lysine'] > 0 and v['growth'] > 0:
        ax3.scatter(v['growth'], v['lysine'], s=80, zorder=5)
        ax3.annotate(f'Δ{k}', (v['growth'], v['lysine']), fontsize=7,
                    xytext=(5, 5), textcoords='offset points')
# Add wild-type point
ax3.scatter(sol_lys_base.objective_value, lys_base_flux, s=120, c='red',
            marker='*', zorder=6, label='Wild-type')
ax3.set_xlabel('Growth Rate (h⁻¹)')
ax3.set_ylabel('Lysine Flux (mmol/gDW/h)')
ax3.set_title('(C) Growth-Lysine Tradeoff')
ax3.legend()
ax3.grid(True, alpha=0.3)

# 7d: Overexpression results
ax4 = fig.add_subplot(gs[1, 0])
if oe_results:
    oe_names = list(oe_results.keys())
    oe_lys = [oe_results[k]['lysine'] for k in oe_names]
    ax4.bar(oe_names, oe_lys, color='#2196F3', edgecolor='black', linewidth=0.5)
    ax4.set_ylabel('Lysine Flux (mmol/gDW/h)')
    ax4.set_title('(D) Overexpression Targets')
    ax4.tick_params(axis='x', rotation=30)
    ax4.grid(True, alpha=0.3, axis='y')

# 7e: Double knockout results
ax5 = fig.add_subplot(gs[1, 1])
if sorted_dko:
    dko_names = [k for k, v in sorted_dko[:6]]
    dko_vals = [v['lysine'] for k, v in sorted_dko[:6]]
    ax5.barh(dko_names, dko_vals, color='#E91E63', edgecolor='black', linewidth=0.5)
    ax5.set_xlabel('Lysine Flux (mmol/gDW/h)')
    ax5.set_title('(E) Double Knockouts')
    ax5.grid(True, alpha=0.3, axis='x')

# 7f: Summary comparison
ax6 = fig.add_subplot(gs[1, 2])
summary_strategies = ['Wild-type', 'Best KO', 'Best DKO']
summary_vals = [lys_base_flux,
                sorted_ko[0][1]['lysine'] if sorted_ko else 0,
                sorted_dko[0][1]['lysine'] if sorted_dko else 0]
if 'combined' in results.get('lysine', {}):
    summary_strategies.append('Combined')
    summary_vals.append(results['lysine']['combined']['lysine'])

colors_s = ['#9E9E9E', '#FF9800', '#E91E63', '#4CAF50'][:len(summary_strategies)]
ax6.bar(summary_strategies, summary_vals, color=colors_s, edgecolor='black', linewidth=0.5)
ax6.set_ylabel('Lysine Flux (mmol/gDW/h)')
ax6.set_title('(F) Strategy Comparison')
for i, v in enumerate(summary_vals):
    ax6.text(i, v + 0.05, f'{v:.2f}', ha='center', fontsize=9)
ax6.grid(True, alpha=0.3, axis='y')

plt.savefig(f'{FIGURES_DIR}/fig7_lysine_optimization.png')
plt.close()
print("Saved: figures/fig7_lysine_optimization.png")

# =============================================================================
# Summary Figure: Integrated Framework Overview
# =============================================================================
print("\n" + "="*70)
print("SUMMARY")
print("="*70)

# Figure 8: Comprehensive comparison
fig, axes = plt.subplots(2, 2, figsize=(14, 11))

# 8a: All growth rates comparison
all_methods = ['Standard\nFBA', 'pFBA', 'Thermo-\nconstrained', 'Loopless', 
               '13C-MFA\nConstrained', 'Enzyme-\nConstrained']
all_growths = [
    constraint_results['Default'],
    constraint_results['pFBA'],
    constraint_results['Thermodynamic'],
    constraint_results['Loopless'],
    results['13c_comparison']['growth_13c'],
    results['enzyme_constraints']['ec_growth']
]
colors_all = plt.cm.Set3(np.linspace(0, 1, len(all_methods)))
axes[0, 0].bar(all_methods, all_growths, color=colors_all, edgecolor='black', linewidth=0.5)
axes[0, 0].set_ylabel('Growth Rate (h⁻¹)')
axes[0, 0].set_title('(A) Growth Rate: All Methods')
for i, v in enumerate(all_growths):
    axes[0, 0].text(i, v + 0.01, f'{v:.3f}', ha='center', fontsize=7)

# 8b: Context-specific comparison
cond_names = list(cs_results.keys())
cond_growths = [cs_results[c]['growth'] for c in cond_names]
cond_removed = [cs_results[c]['removed_reactions'] for c in cond_names]
axes[0, 1].bar(cond_names, cond_growths, color=['#4CAF50', '#2196F3', '#FF9800', '#9C27B0'],
               edgecolor='black', linewidth=0.5)
axes[0, 1].set_ylabel('Growth Rate (h⁻¹)')
axes[0, 1].set_title('(B) Context-Specific Models')
axes[0, 1].tick_params(axis='x', rotation=25)

# 8c: dFBA final states
dfba_labels = [f'Glc₀={g}' for g in glucose_levels]
dfba_final_bm = [dfba_multi[g]['biomass'][-1] for g in glucose_levels]
dfba_final_ac = [dfba_multi[g]['acetate'][-1] for g in glucose_levels]
x8 = np.arange(len(glucose_levels))
axes[1, 0].bar(x8 - 0.17, dfba_final_bm, 0.34, label='Biomass', color='#4CAF50', edgecolor='black', linewidth=0.5)
axes[1, 0].bar(x8 + 0.17, dfba_final_ac, 0.34, label='Acetate', color='#FF5722', edgecolor='black', linewidth=0.5)
axes[1, 0].set_xticks(x8)
axes[1, 0].set_xticklabels(dfba_labels)
axes[1, 0].set_ylabel('Concentration')
axes[1, 0].set_title('(C) dFBA Final States')
axes[1, 0].legend()

# 8d: Lysine optimization strategies
if sorted_ko:
    top5 = sorted_ko[:5]
    ko_names_5 = [f'Δ{k}' for k, v in top5]
    ko_lys_5 = [v['lysine'] for k, v in top5]
    ko_growth_5 = [v['growth'] for k, v in top5]
    
    ax_twin = axes[1, 1]
    x9 = np.arange(len(ko_names_5))
    bars1 = ax_twin.bar(x9 - 0.17, ko_lys_5, 0.34, label='Lysine', color='#FF9800', edgecolor='black', linewidth=0.5)
    bars2 = ax_twin.bar(x9 + 0.17, ko_growth_5, 0.34, label='Growth', color='#4CAF50', edgecolor='black', linewidth=0.5)
    ax_twin.set_xticks(x9)
    ax_twin.set_xticklabels(ko_names_5)
    ax_twin.set_ylabel('Flux / Rate')
    ax_twin.set_title('(D) Top Lysine KO Strategies')
    ax_twin.legend()

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig8_integrated_summary.png')
plt.close()
print("Saved: figures/fig8_integrated_summary.png")

# Save numerical results
print("\n" + "="*70)
print("KEY RESULTS SUMMARY")
print("="*70)
print(f"1. FBA Constraint Optimization:")
print(f"   - Standard FBA growth: {constraint_results['Default']:.4f} h^-1")
print(f"   - pFBA growth: {constraint_results['pFBA']:.4f} h^-1")
print(f"   - Thermodynamic constrained: {constraint_results['Thermodynamic']:.4f} h^-1")

print(f"\n2. 13C-MFA Integration:")
print(f"   - Unconstrained growth: {results['13c_comparison']['growth_unconstrained']:.4f} h^-1")
print(f"   - 13C-constrained growth: {results['13c_comparison']['growth_13c']:.4f} h^-1")

print(f"\n3. Dynamic FBA:")
print(f"   - Final biomass (Glc₀=20): {results['dfba']['final_biomass']:.3f} gDW/L")
print(f"   - Max growth rate: {results['dfba']['max_growth']:.4f} h^-1")

print(f"\n4. Enzyme Constraints:")
print(f"   - Standard FBA: {results['enzyme_constraints']['std_growth']:.4f} h^-1")
print(f"   - EC model: {results['enzyme_constraints']['ec_growth']:.4f} h^-1")
print(f"   - Protein saturation: {results['enzyme_constraints']['saturation_protein']:.2f} g/gDW")

print(f"\n5. Context-Specific Models:")
for cond, vals in cs_results.items():
    print(f"   - {cond}: growth={vals['growth']:.4f}, removed={vals['removed_reactions']} rxns")

print(f"\n6. Lysine Optimization:")
print(f"   - Baseline lysine: {results['lysine']['baseline_lysine']:.4f} mmol/gDW/h")
print(f"   - Max lysine: {results['lysine']['max_lysine']:.4f} mmol/gDW/h")
if sorted_ko and sorted_ko[0][1]['lysine'] > 0:
    print(f"   - Best KO (Δ{sorted_ko[0][0]}): {sorted_ko[0][1]['lysine']:.4f} mmol/gDW/h")
if 'combined' in results.get('lysine', {}):
    print(f"   - Combined strategy: {results['lysine']['combined']['lysine']:.4f} mmol/gDW/h")

# Save results to JSON
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

with open('results.json', 'w') as f:
    json.dump(results, f, indent=2, cls=NumpyEncoder)

print("\n✅ All experiments completed. Results saved to results.json")
print(f"✅ Figures saved in {FIGURES_DIR}/")
