#!/usr/bin/env python3
"""
Molecular Design Framework for Biodegradable Polymers with Controlled Degradation
==================================================================================
Comprehensive computational experiment covering:
1. Hydrolysis rate prediction model
2. Mechanical-degradation tradeoff optimization
3. Michaelis-Menten enzymatic degradation modeling
4. Marine environment degradation simulation
5. Combinatorial copolymer design exploration
6. PLA/PHA/PBS modification case studies
7. ML-based structure-degradability relationship model
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize, differential_evolution
from scipy.integrate import solve_ivp
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.inspection import permutation_importance
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
sns.set_theme(style="whitegrid", font_scale=1.1)

FIGURES_DIR = "figures"

# =============================================================================
# 1. HYDROLYSIS RATE PREDICTION MODEL
# =============================================================================

def hydrolysis_rate_model(bond_type_factor, crystallinity, Mw, T=310, pH=7.4):
    """
    Predict hydrolysis rate constant k_h (day^-1).
    k_h = A * bond_factor * exp(-Ea/RT) * (1 - X_c)^alpha * Mw^(-beta) * f(pH)
    """
    A = 1.0e8          # pre-exponential factor
    Ea = 80000          # activation energy J/mol
    R = 8.314           # gas constant
    alpha = 1.5         # crystallinity exponent
    beta = 0.3          # molecular weight exponent
    pH_ref = 7.4
    
    k_h = (A * bond_type_factor * np.exp(-Ea / (R * T)) *
           (1 - crystallinity) ** alpha *
           (Mw / 1e5) ** (-beta) *
           (1 + 0.3 * np.abs(pH - pH_ref)))
    return k_h

bond_types = {
    'ester': 1.0,
    'amide': 0.3,
    'anhydride': 3.0,
    'carbonate': 1.5,
    'urethane': 0.5,
    'orthoester': 5.0,
}

crystallinities = np.linspace(0.0, 0.8, 50)
Mws = np.logspace(3.5, 6, 50)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# (a) Bond type comparison
for name, factor in bond_types.items():
    kh = [hydrolysis_rate_model(factor, 0.3, 5e4) for _ in range(1)]
    axes[0].barh(name, kh[0], color=plt.cm.Set2(list(bond_types.keys()).index(name) / len(bond_types)))
axes[0].set_xlabel('Hydrolysis Rate $k_h$ (day$^{-1}$)')
axes[0].set_title('(a) Bond Type Dependence')

# (b) Crystallinity dependence
for name, factor in list(bond_types.items())[:3]:
    kh_vals = [hydrolysis_rate_model(factor, xc, 5e4) for xc in crystallinities]
    axes[1].plot(crystallinities * 100, kh_vals, label=name, linewidth=2)
axes[1].set_xlabel('Crystallinity (%)')
axes[1].set_ylabel('$k_h$ (day$^{-1}$)')
axes[1].set_title('(b) Crystallinity Dependence')
axes[1].legend()

# (c) Molecular weight dependence
for name, factor in list(bond_types.items())[:3]:
    kh_vals = [hydrolysis_rate_model(factor, 0.3, mw) for mw in Mws]
    axes[2].loglog(Mws, kh_vals, label=name, linewidth=2)
axes[2].set_xlabel('Molecular Weight (g/mol)')
axes[2].set_ylabel('$k_h$ (day$^{-1}$)')
axes[2].set_title('(c) Molecular Weight Dependence')
axes[2].legend()

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig1_hydrolysis_rate.png', dpi=150, bbox_inches='tight')
plt.close()
print("Figure 1: Hydrolysis rate model - DONE")

# =============================================================================
# 2. MECHANICAL-DEGRADATION TRADEOFF OPTIMIZATION
# =============================================================================

def tensile_strength(crystallinity, Mw, crosslink_density):
    """Estimate tensile strength (MPa)."""
    sigma_a = 20   # amorphous baseline MPa
    sigma_c = 80   # crystalline contribution
    Mw_ref = 1e5
    return (sigma_a + sigma_c * crystallinity) * (1 - np.exp(-Mw / Mw_ref)) * (1 + 50 * crosslink_density)

def elastic_modulus(crystallinity, Mw, crosslink_density):
    """Estimate elastic modulus (GPa)."""
    E_a = 0.5
    E_c = 3.0
    Mw_ref = 1e5
    return (E_a + E_c * crystallinity) * (1 - np.exp(-Mw / Mw_ref)) * (1 + 30 * crosslink_density)

def degradation_rate(crystallinity, Mw, crosslink_density, bond_factor=1.0):
    """Overall degradation rate (day^-1)."""
    kh = hydrolysis_rate_model(bond_factor, crystallinity, Mw)
    crosslink_penalty = np.exp(-100 * crosslink_density)
    return kh * crosslink_penalty

n_samples_trade = 500
xc_samples = np.random.uniform(0, 0.7, n_samples_trade)
mw_samples = 10 ** np.random.uniform(4, 5.5, n_samples_trade)
cl_samples = np.random.uniform(0, 0.05, n_samples_trade)

ts_vals = [tensile_strength(xc, mw, cl) for xc, mw, cl in zip(xc_samples, mw_samples, cl_samples)]
em_vals = [elastic_modulus(xc, mw, cl) for xc, mw, cl in zip(xc_samples, mw_samples, cl_samples)]
dr_vals = [degradation_rate(xc, mw, cl) for xc, mw, cl in zip(xc_samples, mw_samples, cl_samples)]

ts_arr = np.array(ts_vals)
em_arr = np.array(em_vals)
dr_arr = np.array(dr_vals)

# Pareto optimization
def pareto_front(obj1, obj2, maximize_both=True):
    """Find Pareto-optimal points."""
    if maximize_both:
        s1, s2 = obj1, obj2
    else:
        s1, s2 = obj1, -obj2
    pareto_mask = np.ones(len(s1), dtype=bool)
    for i in range(len(s1)):
        for j in range(len(s1)):
            if i != j:
                if s1[j] >= s1[i] and s2[j] >= s2[i] and (s1[j] > s1[i] or s2[j] > s2[i]):
                    pareto_mask[i] = False
                    break
    return pareto_mask

pareto_mask = pareto_front(ts_arr, dr_arr, maximize_both=True)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

sc = axes[0].scatter(ts_arr, dr_arr, c=xc_samples * 100, cmap='viridis', alpha=0.5, s=20)
axes[0].scatter(ts_arr[pareto_mask], dr_arr[pareto_mask], c='red', s=60, marker='*',
               label='Pareto optimal', zorder=5)
axes[0].set_xlabel('Tensile Strength (MPa)')
axes[0].set_ylabel('Degradation Rate (day$^{-1}$)')
axes[0].set_title('(a) Strength-Degradation Tradeoff')
axes[0].legend()
plt.colorbar(sc, ax=axes[0], label='Crystallinity (%)')

sc2 = axes[1].scatter(em_arr, dr_arr, c=mw_samples, cmap='plasma', alpha=0.5, s=20, norm=matplotlib.colors.LogNorm())
axes[1].set_xlabel('Elastic Modulus (GPa)')
axes[1].set_ylabel('Degradation Rate (day$^{-1}$)')
axes[1].set_title('(b) Modulus-Degradation Tradeoff')
plt.colorbar(sc2, ax=axes[1], label='Molecular Weight (g/mol)')

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig2_tradeoff.png', dpi=150, bbox_inches='tight')
plt.close()
print("Figure 2: Tradeoff optimization - DONE")

# =============================================================================
# 3. MICHAELIS-MENTEN ENZYMATIC DEGRADATION
# =============================================================================

def michaelis_menten_degradation(t, S, Vmax, Km, E0, kd=0.01):
    """
    ODE system for enzymatic degradation with enzyme deactivation.
    dS/dt = -Vmax * E * S / (Km + S)
    dE/dt = -kd * E
    """
    E = S[1]
    substrate = S[0]
    dSdt = -Vmax * E * substrate / (Km + substrate)
    dEdt = -kd * E
    return [dSdt, dEdt]

enzymes = {
    'Proteinase K (PLA)':    {'Vmax': 0.5,  'Km': 2.0, 'E0': 1.0, 'kd': 0.005},
    'PHA depolymerase':      {'Vmax': 0.8,  'Km': 1.5, 'E0': 1.0, 'kd': 0.01},
    'Lipase (PBS)':           {'Vmax': 0.3,  'Km': 3.0, 'E0': 1.0, 'kd': 0.008},
    'Cutinase (PLA)':         {'Vmax': 0.6,  'Km': 1.8, 'E0': 1.0, 'kd': 0.007},
}

t_span = (0, 100)
t_eval = np.linspace(0, 100, 500)
S0_init = 10.0  # initial substrate concentration (mg/mL)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for idx, (name, params) in enumerate(enzymes.items()):
    ax = axes[idx // 2, idx % 2]
    sol = solve_ivp(michaelis_menten_degradation, t_span, [S0_init, params['E0']],
                    args=(params['Vmax'], params['Km'], params['E0'], params['kd']),
                    t_eval=t_eval, method='RK45')
    ax.plot(sol.t, sol.y[0], 'b-', linewidth=2, label='Substrate')
    ax.plot(sol.t, sol.y[1] * S0_init, 'r--', linewidth=2, label='Enzyme (scaled)')
    
    # Michaelis-Menten rate at each time point
    rate = params['Vmax'] * sol.y[1] * sol.y[0] / (params['Km'] + sol.y[0])
    ax2 = ax.twinx()
    ax2.plot(sol.t, rate, 'g:', linewidth=1.5, label='Rate')
    ax2.set_ylabel('Rate (mg/mL/day)', color='green')
    
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Concentration (mg/mL)')
    ax.set_title(name)
    ax.legend(loc='upper right')

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig3_michaelis_menten.png', dpi=150, bbox_inches='tight')
plt.close()
print("Figure 3: Michaelis-Menten enzymatic degradation - DONE")

# =============================================================================
# 4. MARINE ENVIRONMENT DEGRADATION SIMULATION
# =============================================================================

def marine_degradation_model(t, state, params):
    """
    Marine degradation simulation.
    state = [M_polymer, M_oligomer, M_monomer, B_biomass]
    """
    M, O, Mon, B = state
    T, pH, k_h0, k_bio, mu_max, Ks, Y, kd = params
    
    # Temperature-dependent hydrolysis (Arrhenius)
    Ea = 75000
    R = 8.314
    T_ref = 298.15
    T_K = T + 273.15
    k_h = k_h0 * np.exp(-Ea / R * (1/T_K - 1/T_ref))
    
    # pH effect
    pH_opt = 8.1
    pH_effect = np.exp(-0.5 * ((pH - pH_opt) / 0.5) ** 2)
    
    # Abiotic hydrolysis: polymer -> oligomer
    dM = -k_h * pH_effect * M
    
    # Biotic degradation: oligomer -> monomer (Monod kinetics)
    mu = mu_max * O / (Ks + O) * (T_K / T_ref)
    dO = k_h * pH_effect * M - mu * B / Y
    
    # Monomer production and biomass uptake
    dMon = mu * B / Y - k_bio * Mon * B
    
    # Biomass growth
    dB = Y * mu * B - kd * B
    
    return [dM, dO, dMon, dB]

marine_conditions = {
    'Tropical Surface (28°C, pH 8.1)':   {'T': 28, 'pH': 8.1},
    'Temperate Surface (15°C, pH 8.1)':  {'T': 15, 'pH': 8.1},
    'Deep Sea (4°C, pH 7.8)':            {'T': 4,  'pH': 7.8},
    'Acidified Ocean (15°C, pH 7.6)':    {'T': 15, 'pH': 7.6},
}

base_params = {
    'k_h0': 0.005,    # base hydrolysis rate
    'k_bio': 0.02,    # biotic monomer uptake
    'mu_max': 0.1,    # max microbial growth rate
    'Ks': 0.5,        # half-saturation
    'Y': 0.4,         # yield coefficient
    'kd': 0.01,       # microbial death rate
}

t_marine = np.linspace(0, 365, 1000)  # 1 year

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

for idx, (name, cond) in enumerate(marine_conditions.items()):
    ax = axes[idx // 2, idx % 2]
    params = (cond['T'], cond['pH'], base_params['k_h0'], base_params['k_bio'],
              base_params['mu_max'], base_params['Ks'], base_params['Y'], base_params['kd'])
    
    sol = solve_ivp(marine_degradation_model, (0, 365), [100, 0, 0, 0.1],
                    args=(params,), t_eval=t_marine, method='RK45', max_step=1.0)
    
    ax.plot(sol.t, sol.y[0], '-', linewidth=2, label='Polymer', color=colors[0])
    ax.plot(sol.t, sol.y[1], '--', linewidth=2, label='Oligomer', color=colors[1])
    ax.plot(sol.t, sol.y[2], ':', linewidth=2, label='Monomer', color=colors[2])
    ax.plot(sol.t, sol.y[3], '-.', linewidth=2, label='Biomass', color=colors[3])
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Mass (mg)')
    ax.set_title(name)
    ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig4_marine_degradation.png', dpi=150, bbox_inches='tight')
plt.close()
print("Figure 4: Marine degradation simulation - DONE")

# Summary: half-life comparison
fig, ax = plt.subplots(figsize=(8, 5))
half_lives = []
for name, cond in marine_conditions.items():
    params = (cond['T'], cond['pH'], base_params['k_h0'], base_params['k_bio'],
              base_params['mu_max'], base_params['Ks'], base_params['Y'], base_params['kd'])
    sol = solve_ivp(marine_degradation_model, (0, 365*3), [100, 0, 0, 0.1],
                    args=(params,), t_eval=np.linspace(0, 365*3, 3000), method='RK45', max_step=1.0)
    # Find half-life
    idx_half = np.argmin(np.abs(sol.y[0] - 50))
    half_lives.append(sol.t[idx_half])

ax.barh(list(marine_conditions.keys()), half_lives, color=sns.color_palette('coolwarm', 4))
ax.set_xlabel('Half-life (days)')
ax.set_title('Polymer Half-life in Marine Environments')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig5_marine_halflife.png', dpi=150, bbox_inches='tight')
plt.close()
print("Figure 5: Marine half-life comparison - DONE")

# =============================================================================
# 5. COMBINATORIAL COPOLYMER DESIGN
# =============================================================================

monomers = {
    'Lactide':         {'bond_factor': 1.0, 'Tg': 60,  'Tm': 175, 'strength_contrib': 60, 'modulus_contrib': 2.5, 'cost': 3.0},
    'Glycolide':       {'bond_factor': 1.5, 'Tg': 36,  'Tm': 225, 'strength_contrib': 70, 'modulus_contrib': 3.0, 'cost': 4.0},
    'Caprolactone':    {'bond_factor': 0.6, 'Tg': -60, 'Tm': 60,  'strength_contrib': 25, 'modulus_contrib': 0.4, 'cost': 2.5},
    'Hydroxybutyrate': {'bond_factor': 0.8, 'Tg': 5,   'Tm': 175, 'strength_contrib': 40, 'modulus_contrib': 1.5, 'cost': 5.0},
    'Butylene Succinate': {'bond_factor': 0.7, 'Tg': -32, 'Tm': 115, 'strength_contrib': 35, 'modulus_contrib': 0.7, 'cost': 2.0},
    'Hydroxyvalerate': {'bond_factor': 0.9, 'Tg': -10, 'Tm': 108, 'strength_contrib': 30, 'modulus_contrib': 1.0, 'cost': 6.0},
}

def evaluate_copolymer(composition):
    """Evaluate copolymer properties from monomer fractions."""
    names = list(monomers.keys())
    props = {
        'degradation_rate': 0,
        'tensile_strength': 0,
        'elastic_modulus': 0,
        'Tg': 0,
        'cost': 0,
    }
    
    for i, name in enumerate(names):
        f = composition[i]
        m = monomers[name]
        props['degradation_rate'] += f * m['bond_factor']
        props['tensile_strength'] += f * m['strength_contrib']
        props['elastic_modulus'] += f * m['modulus_contrib']
        props['Tg'] += f * m['Tg']
        props['cost'] += f * m['cost']
    
    # Non-linear mixing effects
    heterogeneity = 1 - np.sum(np.array(composition) ** 2)
    props['degradation_rate'] *= (1 + 0.3 * heterogeneity)
    props['tensile_strength'] *= (1 - 0.15 * heterogeneity)
    
    return props

# Generate combinatorial library
n_copolymers = 2000
copolymer_data = []
monomer_names = list(monomers.keys())

for _ in range(n_copolymers):
    # Random Dirichlet composition
    comp = np.random.dirichlet(np.ones(len(monomers)))
    props = evaluate_copolymer(comp)
    entry = {f'f_{name}': comp[i] for i, name in enumerate(monomer_names)}
    entry.update(props)
    copolymer_data.append(entry)

df_copoly = pd.DataFrame(copolymer_data)

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# (a) Degradation vs Strength colored by composition
sc = axes[0].scatter(df_copoly['tensile_strength'], df_copoly['degradation_rate'],
                     c=df_copoly['f_Lactide'], cmap='YlOrRd', alpha=0.4, s=10)
axes[0].set_xlabel('Tensile Strength (MPa)')
axes[0].set_ylabel('Degradation Rate (rel.)')
axes[0].set_title('(a) Copolymer Design Space')
plt.colorbar(sc, ax=axes[0], label='Lactide Fraction')

# (b) Ternary-like: top 3 monomers
top3 = ['f_Lactide', 'f_Glycolide', 'f_Caprolactone']
for i, name in enumerate(top3):
    axes[1].scatter(df_copoly[name], df_copoly['degradation_rate'], alpha=0.3, s=10, label=name.replace('f_', ''))
axes[1].set_xlabel('Monomer Fraction')
axes[1].set_ylabel('Degradation Rate (rel.)')
axes[1].set_title('(b) Monomer Fraction vs Degradation')
axes[1].legend()

# (c) Cost-performance
axes[2].scatter(df_copoly['cost'], df_copoly['tensile_strength'], c=df_copoly['degradation_rate'],
               cmap='viridis', alpha=0.4, s=10)
axes[2].set_xlabel('Cost Index')
axes[2].set_ylabel('Tensile Strength (MPa)')
axes[2].set_title('(c) Cost-Performance Map')

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig6_combinatorial.png', dpi=150, bbox_inches='tight')
plt.close()
print("Figure 6: Combinatorial copolymer design - DONE")

# =============================================================================
# 6. PLA/PHA/PBS CASE STUDIES
# =============================================================================

polymer_systems = {
    'PLA': {
        'base': {'xc': 0.35, 'Mw': 1.5e5, 'cl': 0.0, 'bond_factor': 1.0},
        'modifications': {
            'PLA + 10% PEG': {'xc': 0.25, 'Mw': 1.2e5, 'cl': 0.0, 'bond_factor': 1.2},
            'PLA + Nanoclay': {'xc': 0.40, 'Mw': 1.5e5, 'cl': 0.005, 'bond_factor': 0.9},
            'PLA-co-GA (90:10)': {'xc': 0.20, 'Mw': 1.0e5, 'cl': 0.0, 'bond_factor': 1.3},
            'Stereocomplex PLA': {'xc': 0.55, 'Mw': 2.0e5, 'cl': 0.0, 'bond_factor': 0.7},
        }
    },
    'PHA': {
        'base': {'xc': 0.50, 'Mw': 3e5, 'cl': 0.0, 'bond_factor': 0.8},
        'modifications': {
            'P(HB-co-HV) 80:20': {'xc': 0.30, 'Mw': 2.5e5, 'cl': 0.0, 'bond_factor': 0.9},
            'PHA + Chain Extender': {'xc': 0.45, 'Mw': 5e5, 'cl': 0.002, 'bond_factor': 0.75},
            'PHA-g-MA': {'xc': 0.35, 'Mw': 2.8e5, 'cl': 0.01, 'bond_factor': 0.85},
            'PHA/Cellulose NF': {'xc': 0.55, 'Mw': 3e5, 'cl': 0.008, 'bond_factor': 0.7},
        }
    },
    'PBS': {
        'base': {'xc': 0.40, 'Mw': 8e4, 'cl': 0.0, 'bond_factor': 0.7},
        'modifications': {
            'PBS-co-BA (80:20)': {'xc': 0.25, 'Mw': 7e4, 'cl': 0.0, 'bond_factor': 0.85},
            'PBS + TiO2': {'xc': 0.42, 'Mw': 8e4, 'cl': 0.003, 'bond_factor': 0.65},
            'PBS + PLA Blend': {'xc': 0.30, 'Mw': 1e5, 'cl': 0.0, 'bond_factor': 0.9},
            'PBS-co-BF (70:30)': {'xc': 0.20, 'Mw': 6e4, 'cl': 0.0, 'bond_factor': 0.95},
        }
    }
}

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
all_case_results = {}

for pidx, (poly_name, system) in enumerate(polymer_systems.items()):
    ax = axes[pidx]
    base = system['base']
    
    base_ts = tensile_strength(base['xc'], base['Mw'], base['cl'])
    base_dr = degradation_rate(base['xc'], base['Mw'], base['cl'], base['bond_factor'])
    base_em = elastic_modulus(base['xc'], base['Mw'], base['cl'])
    
    names = [f'{poly_name} (base)']
    ts_list = [base_ts]
    dr_list = [base_dr]
    em_list = [base_em]
    
    all_case_results[f'{poly_name}_base'] = {'TS': base_ts, 'EM': base_em, 'DR': base_dr}
    
    for mod_name, mod_params in system['modifications'].items():
        mod_ts = tensile_strength(mod_params['xc'], mod_params['Mw'], mod_params['cl'])
        mod_dr = degradation_rate(mod_params['xc'], mod_params['Mw'], mod_params['cl'], mod_params['bond_factor'])
        mod_em = elastic_modulus(mod_params['xc'], mod_params['Mw'], mod_params['cl'])
        names.append(mod_name)
        ts_list.append(mod_ts)
        dr_list.append(mod_dr)
        em_list.append(mod_em)
        all_case_results[mod_name] = {'TS': mod_ts, 'EM': mod_em, 'DR': mod_dr}
    
    x = np.arange(len(names))
    width = 0.35
    ax.bar(x - width/2, [ts/base_ts for ts in ts_list], width, label='Rel. Strength', color='steelblue')
    ax.bar(x + width/2, [dr/base_dr for dr in dr_list], width, label='Rel. Degradation', color='coral')
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Relative Value')
    ax.set_title(f'{poly_name} Modification Design')
    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
    ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig7_case_studies.png', dpi=150, bbox_inches='tight')
plt.close()
print("Figure 7: PLA/PHA/PBS case studies - DONE")

# =============================================================================
# 7. ML-BASED STRUCTURE-DEGRADABILITY RELATIONSHIP MODEL
# =============================================================================

def generate_polymer_dataset(n=1000):
    """Generate synthetic polymer dataset with molecular descriptors."""
    data = []
    for _ in range(n):
        # Molecular descriptors
        mw = 10 ** np.random.uniform(3.5, 6)
        crystallinity = np.random.uniform(0, 0.7)
        bond_factor = np.random.choice([0.3, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 3.0, 5.0])
        hydrophilicity = np.random.uniform(0, 1)
        crosslink = np.random.uniform(0, 0.05)
        branching = np.random.uniform(0, 0.3)
        Tg = np.random.uniform(-60, 80)
        surface_area = np.random.uniform(0.5, 5.0)
        porosity = np.random.uniform(0, 0.5)
        
        # Calculate "true" degradation rate with noise
        kh = hydrolysis_rate_model(bond_factor, crystallinity, mw)
        bio_factor = 1 + 0.5 * hydrophilicity + 0.3 * surface_area - 0.2 * branching
        cl_factor = np.exp(-100 * crosslink)
        poro_factor = 1 + 2 * porosity
        
        true_rate = kh * bio_factor * cl_factor * poro_factor
        noise = np.random.lognormal(0, 0.15)
        observed_rate = true_rate * noise
        
        data.append({
            'Mw': mw,
            'log_Mw': np.log10(mw),
            'crystallinity': crystallinity,
            'bond_factor': bond_factor,
            'hydrophilicity': hydrophilicity,
            'crosslink_density': crosslink,
            'branching_degree': branching,
            'Tg': Tg,
            'surface_area': surface_area,
            'porosity': porosity,
            'degradation_rate': observed_rate,
        })
    
    return pd.DataFrame(data)

df_ml = generate_polymer_dataset(1500)

feature_cols = ['log_Mw', 'crystallinity', 'bond_factor', 'hydrophilicity',
                'crosslink_density', 'branching_degree', 'Tg', 'surface_area', 'porosity']
target = 'degradation_rate'

X = df_ml[feature_cols].values
y = np.log10(df_ml[target].values)  # log-transform target

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

scaler = StandardScaler()
X_train_s = scaler.fit_transform(X_train)
X_test_s = scaler.transform(X_test)

# Train models
models = {
    'Random Forest': RandomForestRegressor(n_estimators=200, max_depth=15, random_state=42),
    'Gradient Boosting': GradientBoostingRegressor(n_estimators=200, max_depth=6, random_state=42),
}

results = {}
for name, model in models.items():
    model.fit(X_train_s, y_train)
    y_pred = model.predict(X_test_s)
    
    r2 = r2_score(y_test, y_pred)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    mae = mean_absolute_error(y_test, y_pred)
    
    cv_scores = cross_val_score(model, X_train_s, y_train, cv=5, scoring='r2')
    
    results[name] = {
        'model': model,
        'y_pred': y_pred,
        'R2': r2,
        'RMSE': rmse,
        'MAE': mae,
        'CV_R2_mean': cv_scores.mean(),
        'CV_R2_std': cv_scores.std(),
    }
    print(f"{name}: R²={r2:.4f}, RMSE={rmse:.4f}, MAE={mae:.4f}, CV-R²={cv_scores.mean():.4f}±{cv_scores.std():.4f}")

# Figure 8: ML results
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# (a) Predicted vs Actual - RF
y_pred_rf = results['Random Forest']['y_pred']
axes[0, 0].scatter(y_test, y_pred_rf, alpha=0.4, s=15, c='steelblue')
lims = [min(y_test.min(), y_pred_rf.min()), max(y_test.max(), y_pred_rf.max())]
axes[0, 0].plot(lims, lims, 'r--', linewidth=2)
axes[0, 0].set_xlabel('Actual log$_{10}$(k)')
axes[0, 0].set_ylabel('Predicted log$_{10}$(k)')
axes[0, 0].set_title(f'(a) Random Forest (R²={results["Random Forest"]["R2"]:.3f})')

# (b) Predicted vs Actual - GBR
y_pred_gb = results['Gradient Boosting']['y_pred']
axes[0, 1].scatter(y_test, y_pred_gb, alpha=0.4, s=15, c='coral')
axes[0, 1].plot(lims, lims, 'r--', linewidth=2)
axes[0, 1].set_xlabel('Actual log$_{10}$(k)')
axes[0, 1].set_ylabel('Predicted log$_{10}$(k)')
axes[0, 1].set_title(f'(b) Gradient Boosting (R²={results["Gradient Boosting"]["R2"]:.3f})')

# (c) Feature importance
rf_model = results['Random Forest']['model']
imp = rf_model.feature_importances_
sorted_idx = np.argsort(imp)
axes[1, 0].barh(np.array(feature_cols)[sorted_idx], imp[sorted_idx], color='steelblue')
axes[1, 0].set_xlabel('Feature Importance')
axes[1, 0].set_title('(c) RF Feature Importance')

# (d) Permutation importance
perm_imp = permutation_importance(rf_model, X_test_s, y_test, n_repeats=10, random_state=42)
sorted_idx_perm = np.argsort(perm_imp.importances_mean)
axes[1, 1].barh(np.array(feature_cols)[sorted_idx_perm], perm_imp.importances_mean[sorted_idx_perm],
               xerr=perm_imp.importances_std[sorted_idx_perm], color='coral')
axes[1, 1].set_xlabel('Permutation Importance')
axes[1, 1].set_title('(d) Permutation Importance')

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig8_ml_model.png', dpi=150, bbox_inches='tight')
plt.close()
print("Figure 8: ML model results - DONE")

# Residual analysis
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for idx, (name, res) in enumerate(results.items()):
    residuals = y_test - res['y_pred']
    axes[idx].hist(residuals, bins=30, color='steelblue' if idx == 0 else 'coral', alpha=0.7, edgecolor='black')
    axes[idx].set_xlabel('Residual (log$_{10}$ scale)')
    axes[idx].set_ylabel('Count')
    axes[idx].set_title(f'{name} Residuals')
    axes[idx].axvline(x=0, color='red', linestyle='--')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig9_residuals.png', dpi=150, bbox_inches='tight')
plt.close()
print("Figure 9: Residual analysis - DONE")

# =============================================================================
# 8. COMPREHENSIVE HEATMAP: Structure-Property Relationships
# =============================================================================

corr_cols = ['log_Mw', 'crystallinity', 'bond_factor', 'hydrophilicity',
             'crosslink_density', 'branching_degree', 'Tg', 'surface_area', 'porosity',
             'degradation_rate']
df_corr = df_ml[corr_cols].copy()
df_corr['degradation_rate'] = np.log10(df_corr['degradation_rate'])

fig, ax = plt.subplots(figsize=(10, 8))
corr_matrix = df_corr.corr()
mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
sns.heatmap(corr_matrix, mask=mask, annot=True, fmt='.2f', cmap='RdBu_r',
            center=0, square=True, ax=ax, vmin=-1, vmax=1)
ax.set_title('Structure-Property Correlation Matrix')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig10_correlation.png', dpi=150, bbox_inches='tight')
plt.close()
print("Figure 10: Correlation matrix - DONE")

# =============================================================================
# SAVE RESULTS SUMMARY
# =============================================================================

print("\n" + "="*60)
print("EXPERIMENT COMPLETE - Summary of Numerical Results")
print("="*60)

print("\n--- ML Model Performance ---")
for name, res in results.items():
    print(f"{name}:")
    print(f"  R² = {res['R2']:.4f}")
    print(f"  RMSE = {res['RMSE']:.4f}")
    print(f"  MAE = {res['MAE']:.4f}")
    print(f"  CV-R² = {res['CV_R2_mean']:.4f} ± {res['CV_R2_std']:.4f}")

print("\n--- Marine Half-lives (days) ---")
for cond_name, hl in zip(marine_conditions.keys(), half_lives):
    print(f"  {cond_name}: {hl:.1f} days")

print("\n--- Case Study Results (Relative to Base) ---")
for sys_name, system in polymer_systems.items():
    base_key = f'{sys_name}_base'
    base_r = all_case_results[base_key]
    print(f"\n  {sys_name} Base: TS={base_r['TS']:.1f} MPa, EM={base_r['EM']:.2f} GPa, DR={base_r['DR']:.6f} day⁻¹")
    for mod_name in system['modifications'].keys():
        mod_r = all_case_results[mod_name]
        print(f"    {mod_name}: TS={mod_r['TS']:.1f} MPa ({mod_r['TS']/base_r['TS']:.2f}x), "
              f"DR={mod_r['DR']:.6f} ({mod_r['DR']/base_r['DR']:.2f}x)")

print("\n--- Top Feature Importances (RF) ---")
for feat, imp_val in sorted(zip(feature_cols, rf_model.feature_importances_), key=lambda x: -x[1]):
    print(f"  {feat}: {imp_val:.4f}")

# Save numerical results to CSV
df_results_summary = pd.DataFrame({
    'Model': list(results.keys()),
    'R2': [r['R2'] for r in results.values()],
    'RMSE': [r['RMSE'] for r in results.values()],
    'MAE': [r['MAE'] for r in results.values()],
    'CV_R2_mean': [r['CV_R2_mean'] for r in results.values()],
    'CV_R2_std': [r['CV_R2_std'] for r in results.values()],
})
df_results_summary.to_csv('ml_results.csv', index=False)

# Save case study results
pd.DataFrame(all_case_results).T.to_csv('case_study_results.csv')

print("\nAll figures saved to figures/ directory.")
print("Numerical results saved to ml_results.csv and case_study_results.csv")
