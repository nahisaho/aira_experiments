"""
Systems Biology Framework for Predicting Diet-Gut Microbiota Interactions
Implements: SHIME digestion dynamics, gLV competition, SCFA flux prediction,
            long-term diet simulation, probiotic/prebiotic effects, fermented food case study
"""

import sys, os, random
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.integrate import solve_ivp, odeint
from scipy.stats import pearsonr, spearmanr, ttest_ind
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Reproducibility
# ============================================================
np.random.seed(42)
random.seed(42)

FIGURES_DIR = '/app/projects/dd3351e2-1353-4208-8d78-052156f6d26a/workspace/figures'
DATA_DIR = '/app/projects/dd3351e2-1353-4208-8d78-052156f6d26a/workspace/data/raw'
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

print("=" * 60)
print("SYSTEMS BIOLOGY FRAMEWORK: DIET-GUT MICROBIOTA INTERACTIONS")
print("=" * 60)
print(f"Python: {sys.version.split()[0]}")
print(f"NumPy: {np.__version__}, Pandas: {pd.__version__}")

# ============================================================
# MODULE 1: SHIME Digestion Dynamics Model
# ============================================================
print("\n[MODULE 1] SHIME-inspired digestion dynamics model")

def shime_digestion_model(t, y, params):
    """
    Simplified SHIME (Simulator of Human Intestinal Microbial Ecosystem)
    Compartments: Stomach (S), Small Intestine (SI), Large Intestine proximal (LI_p), 
                  Large Intestine distal (LI_d)
    y = [substrate_S, substrate_SI, substrate_LI_p, substrate_LI_d,
         acetate, propionate, butyrate, pH_LI]
    """
    S, SI, LI_p, LI_d, acetate, propionate, butyrate, pH = y
    
    k_gastric  = params['k_gastric']   # gastric emptying rate (h^-1)
    k_SI_abs   = params['k_SI_abs']    # small intestinal absorption (h^-1)
    k_transit  = params['k_transit']   # intestinal transit rate (h^-1)
    k_ferm_p   = params['k_ferm_p']    # fermentation rate proximal LI (h^-1)
    k_ferm_d   = params['k_ferm_d']    # fermentation rate distal LI (h^-1)
    Y_acetate  = params['Y_acetate']   # acetate yield (mol/mol substrate)
    Y_propionate = params['Y_propionate']
    Y_butyrate = params['Y_butyrate']
    
    # Gastric emptying (first-order)
    dS_dt  = -k_gastric * S
    dSI_dt = k_gastric * S - k_SI_abs * SI - k_transit * SI
    
    # Large intestine – Michaelis-Menten fermentation
    Km = 5.0   # mM, half-saturation for fermentation
    ferm_rate_p = k_ferm_p * LI_p / (Km + LI_p)
    ferm_rate_d = k_ferm_d * LI_d / (Km + LI_d)
    
    dLI_p_dt = k_transit * SI - ferm_rate_p - k_transit * LI_p
    dLI_d_dt = k_transit * LI_p - ferm_rate_d
    
    # SCFA production
    total_ferm = ferm_rate_p + ferm_rate_d
    d_acetate_dt    = Y_acetate    * total_ferm - 0.05 * acetate
    d_propionate_dt = Y_propionate * total_ferm - 0.05 * propionate
    d_butyrate_dt   = Y_butyrate   * total_ferm - 0.05 * butyrate
    
    # pH dynamics (buffered by bicarbonate and SCFAs)
    total_SCFA = acetate + propionate + butyrate
    dpH_dt = -0.1 * (total_SCFA - 80) / 80  # approach equilibrium pH ~6.5
    
    return [dS_dt, dSI_dt, dLI_p_dt, dLI_d_dt, 
            d_acetate_dt, d_propionate_dt, d_butyrate_dt, dpH_dt]

# Define diet scenarios
diet_scenarios = {
    'High Fiber': {
        'k_gastric': 0.8, 'k_SI_abs': 0.3, 'k_transit': 0.5,
        'k_ferm_p': 2.5, 'k_ferm_d': 1.5,
        'Y_acetate': 0.50, 'Y_propionate': 0.25, 'Y_butyrate': 0.25,
        'initial_substrate': 80.0
    },
    'Western Diet': {
        'k_gastric': 1.2, 'k_SI_abs': 0.7, 'k_transit': 0.4,
        'k_ferm_p': 1.0, 'k_ferm_d': 0.5,
        'Y_acetate': 0.60, 'Y_propionate': 0.20, 'Y_butyrate': 0.20,
        'initial_substrate': 40.0
    },
    'Mediterranean': {
        'k_gastric': 0.9, 'k_SI_abs': 0.4, 'k_transit': 0.5,
        'k_ferm_p': 2.0, 'k_ferm_d': 1.2,
        'Y_acetate': 0.48, 'Y_propionate': 0.27, 'Y_butyrate': 0.25,
        'initial_substrate': 65.0
    },
    'Low Carb': {
        'k_gastric': 1.1, 'k_SI_abs': 0.8, 'k_transit': 0.35,
        'k_ferm_p': 0.6, 'k_ferm_d': 0.3,
        'Y_acetate': 0.55, 'Y_propionate': 0.25, 'Y_butyrate': 0.20,
        'initial_substrate': 25.0
    }
}

t_span = (0, 24)
t_eval = np.linspace(0, 24, 200)

shime_results = {}
for diet, params in diet_scenarios.items():
    substrate0 = params['initial_substrate']
    y0 = [substrate0, 0.0, 0.0, 0.0, 2.0, 1.0, 0.5, 7.0]  # initial pH 7
    sol = solve_ivp(shime_digestion_model, t_span, y0, t_eval=t_eval,
                    args=(params,), method='RK45', rtol=1e-6, atol=1e-8)
    shime_results[diet] = sol

print("  SHIME simulation complete.")

# Extract final SCFA concentrations (at t=24h)
shime_summary = {}
for diet, sol in shime_results.items():
    shime_summary[diet] = {
        'acetate_mM':    sol.y[4, -1],
        'propionate_mM': sol.y[5, -1],
        'butyrate_mM':   sol.y[6, -1],
        'total_SCFA_mM': sol.y[4, -1] + sol.y[5, -1] + sol.y[6, -1],
        'final_pH':      sol.y[7, -1]
    }
    print(f"  {diet}: Acetate={shime_summary[diet]['acetate_mM']:.2f}, "
          f"Propionate={shime_summary[diet]['propionate_mM']:.2f}, "
          f"Butyrate={shime_summary[diet]['butyrate_mM']:.2f} mM, "
          f"pH={shime_summary[diet]['final_pH']:.2f}")

# ============================================================
# MODULE 2: Generalized Lotka-Volterra (gLV) Community Model
# ============================================================
print("\n[MODULE 2] Generalized Lotka-Volterra community dynamics")

def gLV_model(t, N, r, A, K):
    """
    dN_i/dt = N_i * (r_i + sum_j(A_ij * N_j / K_j))
    r: intrinsic growth rates
    A: interaction matrix (competition/facilitation)
    K: carrying capacities
    """
    N = np.maximum(N, 0)  # non-negative abundances
    dN_dt = N * (r + A @ (N / K))
    return dN_dt

# 8 representative gut bacterial genera
species = ['Bacteroides', 'Prevotella', 'Ruminococcus', 'Faecalibacterium',
           'Bifidobacterium', 'Lactobacillus', 'Akkermansia', 'Clostridium']
n_species = len(species)

# Intrinsic growth rates (per day)
r_rates = np.array([0.8, 0.7, 0.6, 0.7, 0.9, 0.85, 0.65, 0.5])

# Carrying capacities (relative abundance, sum = 1)
K_cap = np.ones(n_species)

# Interaction matrix (negative = competition, positive = facilitation)
# Based on known ecological interactions from literature
np.random.seed(42)
A_base = -0.3 * np.eye(n_species)  # self-limitation

# Known positive interactions (cross-feeding)
cross_feeding = {
    (4, 2): 0.15,   # Bifidobacterium -> Ruminococcus (acetate/lactate cross-feeding)
    (5, 4): 0.10,   # Lactobacillus -> Bifidobacterium
    (6, 3): 0.12,   # Akkermansia -> Faecalibacterium
    (3, 2): 0.08,   # Faecalibacterium -> Ruminococcus
}
for (i, j), v in cross_feeding.items():
    A_base[i, j] = v
    A_base[j, i] = v * 0.5  # asymmetric facilitation

# Competition for resources
competition_pairs = [(0, 1), (0, 2), (1, 3), (2, 7), (6, 0)]
for (i, j) in competition_pairs:
    A_base[i, j] -= 0.10
    A_base[j, i] -= 0.10

# Diet-modulated resource availability (affects growth rates)
def diet_modulated_growth(diet_type):
    """Return growth rates modulated by diet composition"""
    base = r_rates.copy()
    if diet_type == 'High Fiber':
        base[2] *= 1.5   # Ruminococcus (fiber degrader) boosted
        base[3] *= 1.4   # Faecalibacterium (butyrate producer) boosted
        base[4] *= 1.3   # Bifidobacterium boosted
        base[6] *= 1.2   # Akkermansia boosted
        base[0] *= 0.9   # Bacteroides slightly reduced
    elif diet_type == 'Western Diet':
        base[0] *= 1.3   # Bacteroides boosted
        base[7] *= 1.2   # Clostridium boosted
        base[3] *= 0.7   # Faecalibacterium reduced
        base[6] *= 0.8   # Akkermansia reduced
    elif diet_type == 'Mediterranean':
        base[3] *= 1.3   # Faecalibacterium boosted
        base[4] *= 1.2   # Bifidobacterium boosted
        base[6] *= 1.15  # Akkermansia boosted
        base[1] *= 1.1   # Prevotella slightly boosted (plant foods)
    elif diet_type == 'Low Carb':
        base[0] *= 1.2   # Bacteroides boosted
        base[5] *= 0.9   # Lactobacillus reduced
        base[3] *= 0.85  # Faecalibacterium reduced
    return base

# Simulate for 60 days
t_glv = (0, 60)
t_eval_glv = np.linspace(0, 60, 600)

# Initial conditions: approximating healthy microbiome
N0_healthy = np.array([0.25, 0.15, 0.12, 0.10, 0.08, 0.07, 0.05, 0.18])

gLV_results = {}
for diet in diet_scenarios.keys():
    r_mod = diet_modulated_growth(diet)
    sol_glv = solve_ivp(gLV_model, t_glv, N0_healthy, t_eval=t_eval_glv,
                        args=(r_mod, A_base, K_cap), method='RK45', 
                        rtol=1e-6, atol=1e-9)
    # Normalize to relative abundance
    total = sol_glv.y.sum(axis=0)
    total = np.where(total > 0, total, 1)
    rel_abund = sol_glv.y / total[np.newaxis, :]
    gLV_results[diet] = {'time': t_eval_glv, 'rel_abund': rel_abund, 'sol': sol_glv}

# Compute diversity metrics (Shannon entropy)
def shannon_entropy(p):
    p = p[p > 1e-10]
    return -np.sum(p * np.log(p))

def simpson_diversity(p):
    p = p[p > 1e-10]
    return 1 - np.sum(p**2)

diversity_results = {}
for diet, res in gLV_results.items():
    shannon_final = shannon_entropy(res['rel_abund'][:, -1])
    simpson_final = simpson_diversity(res['rel_abund'][:, -1])
    diversity_results[diet] = {
        'shannon': shannon_final,
        'simpson': simpson_final,
        'final_composition': dict(zip(species, res['rel_abund'][:, -1]))
    }
    print(f"  {diet}: Shannon={shannon_final:.3f}, Simpson={simpson_final:.3f}")

print("  gLV simulation complete.")

# ============================================================
# MODULE 3: SCFA Flux Prediction via Flux Balance Analysis (FBA)
# ============================================================
print("\n[MODULE 3] SCFA flux prediction (FBA-inspired)")

# FBA-inspired approach: predict SCFA production based on 
# community composition and substrate availability
# Using stoichiometric coefficients from literature (Flint et al. 2012)

# Substrate utilization per species (g substrate per g biomass per day)
substrate_preference = {
    'Bacteroides':    {'inulin': 0.7, 'pectin': 0.8, 'starch': 0.9, 'arabinoxylan': 0.6},
    'Prevotella':     {'inulin': 0.5, 'pectin': 0.3, 'starch': 0.4, 'arabinoxylan': 0.8},
    'Ruminococcus':   {'inulin': 0.3, 'pectin': 0.2, 'starch': 0.8, 'arabinoxylan': 0.9},
    'Faecalibacterium':{'inulin':0.6, 'pectin': 0.4, 'starch': 0.3, 'arabinoxylan': 0.5},
    'Bifidobacterium':{'inulin': 0.9, 'pectin': 0.3, 'starch': 0.4, 'arabinoxylan': 0.7},
    'Lactobacillus':  {'inulin': 0.8, 'pectin': 0.2, 'starch': 0.5, 'arabinoxylan': 0.4},
    'Akkermansia':    {'inulin': 0.3, 'pectin': 0.9, 'starch': 0.2, 'arabinoxylan': 0.3},
    'Clostridium':    {'inulin': 0.4, 'pectin': 0.3, 'starch': 0.6, 'arabinoxylan': 0.4},
}

# SCFA yield per unit substrate consumed per species
# (acetate:propionate:butyrate ratios from metabolic pathways)
scfa_yield = {
    'Bacteroides':    {'acetate': 0.55, 'propionate': 0.30, 'butyrate': 0.15, 'succinate': 0.20},
    'Prevotella':     {'acetate': 0.45, 'propionate': 0.40, 'butyrate': 0.10, 'succinate': 0.25},
    'Ruminococcus':   {'acetate': 0.50, 'propionate': 0.15, 'butyrate': 0.35, 'succinate': 0.10},
    'Faecalibacterium':{'acetate':0.35, 'propionate': 0.15, 'butyrate': 0.50, 'succinate': 0.10},
    'Bifidobacterium':{'acetate': 0.65, 'propionate': 0.10, 'butyrate': 0.05, 'succinate': 0.10},
    'Lactobacillus':  {'acetate': 0.55, 'propionate': 0.15, 'butyrate': 0.10, 'succinate': 0.15},
    'Akkermansia':    {'acetate': 0.45, 'propionate': 0.35, 'butyrate': 0.15, 'succinate': 0.30},
    'Clostridium':    {'acetate': 0.45, 'propionate': 0.15, 'butyrate': 0.30, 'succinate': 0.15},
}

# Diet compositions (relative substrate availability)
diet_compositions = {
    'High Fiber':    {'inulin': 0.35, 'pectin': 0.25, 'starch': 0.20, 'arabinoxylan': 0.20},
    'Western Diet':  {'inulin': 0.05, 'pectin': 0.05, 'starch': 0.80, 'arabinoxylan': 0.10},
    'Mediterranean': {'inulin': 0.20, 'pectin': 0.25, 'starch': 0.30, 'arabinoxylan': 0.25},
    'Low Carb':      {'inulin': 0.10, 'pectin': 0.10, 'starch': 0.20, 'arabinoxylan': 0.60},
}

def predict_scfa_flux(composition, diet_type, total_substrate=50.0):
    """Predict community-level SCFA fluxes based on FBA-inspired calculation"""
    diet_comp = diet_compositions[diet_type]
    fluxes = {'acetate': 0, 'propionate': 0, 'butyrate': 0, 'succinate': 0}
    
    for sp, rel_abund in composition.items():
        for substrate, sub_avail in diet_comp.items():
            # Substrate consumed by this species
            consumed = rel_abund * substrate_preference[sp][substrate] * sub_avail * total_substrate
            # SCFA produced
            for scfa in fluxes:
                fluxes[scfa] += consumed * scfa_yield[sp][scfa]
    
    return fluxes

# Calculate for each diet at steady state
fba_results = {}
for diet in diet_scenarios.keys():
    final_comp = diversity_results[diet]['final_composition']
    fluxes = predict_scfa_flux(final_comp, diet)
    total = sum(fluxes[k] for k in ['acetate', 'propionate', 'butyrate'])
    fba_results[diet] = {**fluxes, 'total_SCFA': total}
    print(f"  {diet}: Acetate={fluxes['acetate']:.2f}, "
          f"Propionate={fluxes['propionate']:.2f}, "
          f"Butyrate={fluxes['butyrate']:.2f} mmol/day, "
          f"Total={total:.2f}")

print("  FBA-inspired SCFA flux prediction complete.")

# ============================================================
# MODULE 4: Long-term Diet-Microbiome Dynamics (180 days)
# ============================================================
print("\n[MODULE 4] Long-term diet transition dynamics (180 days)")

def simulate_diet_transition(diet_from, diet_to, transition_day, total_days=180):
    """Simulate microbiome response to diet transition"""
    t_span_full = (0, total_days)
    t_eval_full = np.linspace(0, total_days, total_days * 5)
    
    # Phase 1: initial diet
    r1 = diet_modulated_growth(diet_from)
    sol1 = solve_ivp(gLV_model, (0, transition_day), N0_healthy, 
                     t_eval=np.linspace(0, transition_day, transition_day * 5),
                     args=(r1, A_base, K_cap), method='RK45', rtol=1e-6, atol=1e-9)
    
    # Phase 2: new diet (starting from final state of phase 1)
    N_at_transition = sol1.y[:, -1]
    r2 = diet_modulated_growth(diet_to)
    sol2 = solve_ivp(gLV_model, (transition_day, total_days), N_at_transition,
                     t_eval=np.linspace(transition_day, total_days, (total_days - transition_day) * 5),
                     args=(r2, A_base, K_cap), method='RK45', rtol=1e-6, atol=1e-9)
    
    # Combine
    t_combined = np.concatenate([sol1.t, sol2.t[1:]])
    N_combined = np.concatenate([sol1.y, sol2.y[:, 1:]], axis=1)
    total_N = N_combined.sum(axis=0)
    total_N = np.where(total_N > 0, total_N, 1)
    rel_abund_combined = N_combined / total_N[np.newaxis, :]
    
    return t_combined, rel_abund_combined

# Simulate two key transitions
t_WtoHF, comp_WtoHF = simulate_diet_transition('Western Diet', 'High Fiber', 30, 120)
t_HFtoW, comp_HFtoW = simulate_diet_transition('High Fiber', 'Western Diet', 30, 120)

# Calculate Shannon diversity over time
shannon_WtoHF = [shannon_entropy(comp_WtoHF[:, i]) for i in range(comp_WtoHF.shape[1])]
shannon_HFtoW = [shannon_entropy(comp_HFtoW[:, i]) for i in range(comp_HFtoW.shape[1])]

# Recovery time (days to reach 90% of new steady state Shannon)
target_div_HF = diversity_results['High Fiber']['shannon']
target_div_W  = diversity_results['Western Diet']['shannon']

# Find recovery time for Western->High Fiber
recovery_time_WtoHF = None
for i, (t, sh) in enumerate(zip(t_WtoHF, shannon_WtoHF)):
    if t > 30 and sh >= 0.9 * target_div_HF:
        recovery_time_WtoHF = t - 30
        break

print(f"  Western → High Fiber transition: recovery ~{recovery_time_WtoHF:.1f} days")
print(f"  Target Shannon (High Fiber): {target_div_HF:.3f}")
print(f"  Final Shannon after transition: {shannon_WtoHF[-1]:.3f}")

print("  Long-term dynamics simulation complete.")

# ============================================================
# MODULE 5: Probiotic and Prebiotic Effect Prediction
# ============================================================
print("\n[MODULE 5] Probiotic/prebiotic effect prediction")

def simulate_probiotic_intervention(base_diet, probiotic_species_idx, 
                                     dose=0.05, duration=30, total_days=60):
    """
    Simulate probiotic (Lactobacillus=5 or Bifidobacterium=4) supplementation.
    probiotic_species_idx: index in species list
    dose: initial abundance boost
    """
    N0_with_probiotic = N0_healthy.copy()
    N0_with_probiotic[probiotic_species_idx] += dose
    # Normalize
    N0_with_probiotic /= N0_with_probiotic.sum()
    
    r_mod = diet_modulated_growth(base_diet)
    
    # Phase 1: intervention
    t_eval_p1 = np.linspace(0, duration, duration * 5)
    sol_p1 = solve_ivp(gLV_model, (0, duration), N0_with_probiotic,
                       t_eval=t_eval_p1, args=(r_mod, A_base, K_cap),
                       method='RK45', rtol=1e-6, atol=1e-9)
    
    # Phase 2: post-intervention (no more supplementation)
    N_post = sol_p1.y[:, -1]
    t_eval_p2 = np.linspace(duration, total_days, (total_days - duration) * 5)
    sol_p2 = solve_ivp(gLV_model, (duration, total_days), N_post,
                       t_eval=t_eval_p2, args=(r_mod, A_base, K_cap),
                       method='RK45', rtol=1e-6, atol=1e-9)
    
    t_comb = np.concatenate([sol_p1.t, sol_p2.t[1:]])
    N_comb = np.concatenate([sol_p1.y, sol_p2.y[:, 1:]], axis=1)
    total_N = N_comb.sum(axis=0)
    total_N = np.where(total_N > 0, total_N, 1)
    return t_comb, N_comb / total_N[np.newaxis, :]

def simulate_prebiotic(base_diet, prebiotic_boost, total_days=60):
    """Simulate prebiotic (inulin) supplementation - modifies growth rates"""
    r_mod = diet_modulated_growth(base_diet).copy()
    # Inulin boosts Bifidobacterium (idx=4), Lactobacillus (idx=5), Ruminococcus (idx=2)
    r_mod[4] *= (1 + prebiotic_boost)
    r_mod[5] *= (1 + prebiotic_boost * 0.7)
    r_mod[2] *= (1 + prebiotic_boost * 0.5)
    
    t_eval_pb = np.linspace(0, total_days, total_days * 5)
    sol_pb = solve_ivp(gLV_model, (0, total_days), N0_healthy,
                       t_eval=t_eval_pb, args=(r_mod, A_base, K_cap),
                       method='RK45', rtol=1e-6, atol=1e-9)
    total_N = sol_pb.y.sum(axis=0)
    total_N = np.where(total_N > 0, total_N, 1)
    return sol_pb.t, sol_pb.y / total_N[np.newaxis, :]

# Run probiotic/prebiotic simulations
probiotic_results = {}
# Lactobacillus probiotic (idx=5)
t_lb, comp_lb = simulate_probiotic_intervention('Western Diet', 5, dose=0.08)
probiotic_results['Lactobacillus probiotic'] = {'t': t_lb, 'comp': comp_lb}
# Bifidobacterium probiotic (idx=4)
t_bif, comp_bif = simulate_probiotic_intervention('Western Diet', 4, dose=0.08)
probiotic_results['Bifidobacterium probiotic'] = {'t': t_bif, 'comp': comp_bif}

# Inulin prebiotic (25% growth rate boost)
t_pre, comp_pre = simulate_prebiotic('Western Diet', prebiotic_boost=0.25)
probiotic_results['Inulin prebiotic'] = {'t': t_pre, 'comp': comp_pre}

# FOS prebiotic (combined boost)
t_fos, comp_fos = simulate_prebiotic('Western Diet', prebiotic_boost=0.40)
probiotic_results['FOS prebiotic'] = {'t': t_fos, 'comp': comp_fos}

# Compare Shannon diversity at end of interventions
print("  Intervention effects on Shannon diversity:")
for name, res in probiotic_results.items():
    sh = shannon_entropy(res['comp'][:, -1])
    print(f"    {name}: Shannon={sh:.3f}")
print("  Probiotic/prebiotic simulation complete.")

# ============================================================
# MODULE 6: Fermented Food Case Study
# ============================================================
print("\n[MODULE 6] Fermented food case study")

# Fermented food provides live cultures + bioactive metabolites
# Model: combination of probiotic effect + prebiotic fiber effect
# Based on: Wastyk et al. 2021 (Cell) fermented food diet study

fermented_foods = {
    'Yogurt':      {'lactobacillus_dose': 0.06, 'bifidobacterium_dose': 0.04, 'prebiotic': 0.10},
    'Kefir':       {'lactobacillus_dose': 0.08, 'bifidobacterium_dose': 0.03, 'prebiotic': 0.08},
    'Kimchi':      {'lactobacillus_dose': 0.07, 'bifidobacterium_dose': 0.01, 'prebiotic': 0.20},
    'Sauerkraut':  {'lactobacillus_dose': 0.05, 'bifidobacterium_dose': 0.01, 'prebiotic': 0.15},
    'Kombucha':    {'lactobacillus_dose': 0.03, 'bifidobacterium_dose': 0.02, 'prebiotic': 0.05},
}

def simulate_fermented_food(food_name, food_params, base_diet='Western Diet', total_days=60):
    """Combined probiotic + prebiotic effect from fermented food"""
    N0_ff = N0_healthy.copy()
    N0_ff[5] += food_params['lactobacillus_dose']  # Lactobacillus
    N0_ff[4] += food_params['bifidobacterium_dose'] # Bifidobacterium
    N0_ff /= N0_ff.sum()
    
    r_mod = diet_modulated_growth(base_diet).copy()
    r_mod[4] *= (1 + food_params['prebiotic'])
    r_mod[5] *= (1 + food_params['prebiotic'] * 0.8)
    r_mod[2] *= (1 + food_params['prebiotic'] * 0.4)  # Ruminococcus
    r_mod[3] *= (1 + food_params['prebiotic'] * 0.3)  # Faecalibacterium
    
    t_eval_ff = np.linspace(0, total_days, total_days * 5)
    sol_ff = solve_ivp(gLV_model, (0, total_days), N0_ff, t_eval=t_eval_ff,
                       args=(r_mod, A_base, K_cap), method='RK45', 
                       rtol=1e-6, atol=1e-9)
    total_N = sol_ff.y.sum(axis=0)
    total_N = np.where(total_N > 0, total_N, 1)
    return sol_ff.t, sol_ff.y / total_N[np.newaxis, :]

ff_results = {}
for food, params in fermented_foods.items():
    t_ff, comp_ff = simulate_fermented_food(food, params)
    sh_ff = shannon_entropy(comp_ff[:, -1])
    sim_ff = simpson_diversity(comp_ff[:, -1])
    # Predict SCFA with enriched composition
    final_comp_ff = dict(zip(species, comp_ff[:, -1]))
    scfa_ff = predict_scfa_flux(final_comp_ff, 'Western Diet', total_substrate=45.0)
    ff_results[food] = {
        't': t_ff, 'comp': comp_ff,
        'shannon': sh_ff, 'simpson': sim_ff,
        'scfa': scfa_ff
    }
    print(f"  {food}: Shannon={sh_ff:.3f}, Butyrate={scfa_ff['butyrate']:.2f} mmol/day")

# Baseline Western Diet (no fermented food)
_, comp_baseline = simulate_fermented_food('Baseline', 
    {'lactobacillus_dose': 0, 'bifidobacterium_dose': 0, 'prebiotic': 0.0})
sh_baseline = shannon_entropy(comp_baseline[:, -1])
base_comp = dict(zip(species, comp_baseline[:, -1]))
base_scfa = predict_scfa_flux(base_comp, 'Western Diet', total_substrate=45.0)
ff_results['Baseline (No FF)'] = {
    'shannon': sh_baseline, 'simpson': simpson_diversity(comp_baseline[:, -1]),
    'scfa': base_scfa
}
print(f"  Baseline (no fermented food): Shannon={sh_baseline:.3f}, "
      f"Butyrate={base_scfa['butyrate']:.2f} mmol/day")
print("  Fermented food case study complete.")

# ============================================================
# MODULE 7: ML-based SCFA Predictor (Random Forest + GBM)
# ============================================================
print("\n[MODULE 7] Machine learning SCFA prediction model")

# Generate synthetic dataset of microbiome profiles -> SCFA outputs
n_samples = 500
np.random.seed(42)

# Generate diverse microbiome compositions using Dirichlet distribution
# Different concentration parameters to simulate diverse populations
alpha_high_fiber = np.array([2.0, 1.5, 3.0, 3.5, 2.5, 2.0, 2.0, 1.0])
alpha_western    = np.array([4.0, 2.0, 1.5, 1.0, 1.5, 1.5, 0.8, 3.0])
alpha_healthy    = np.array([3.0, 2.0, 2.5, 2.5, 2.0, 1.8, 1.5, 2.0])

n_each = n_samples // 4
compositions = np.concatenate([
    np.random.dirichlet(alpha_high_fiber, n_each),
    np.random.dirichlet(alpha_western, n_each),
    np.random.dirichlet(alpha_healthy, n_each),
    np.random.dirichlet((alpha_high_fiber + alpha_western) / 2, n_samples - 3 * n_each)
])

diet_labels = (
    ['High Fiber'] * n_each + 
    ['Western Diet'] * n_each + 
    ['Mediterranean'] * n_each + 
    ['Low Carb'] * (n_samples - 3 * n_each)
)

# Compute SCFA for each sample
X_data = []
y_acetate, y_propionate, y_butyrate, y_total = [], [], [], []

for i, (comp, diet) in enumerate(zip(compositions, diet_labels)):
    comp_dict = dict(zip(species, comp))
    fluxes = predict_scfa_flux(comp_dict, diet)
    
    # Feature vector: composition + diet encoding
    diet_enc = {'High Fiber': [1,0,0,0], 'Western Diet': [0,1,0,0], 
                'Mediterranean': [0,0,1,0], 'Low Carb': [0,0,0,1]}
    features = list(comp) + diet_enc[diet]
    
    X_data.append(features)
    y_acetate.append(fluxes['acetate'])
    y_propionate.append(fluxes['propionate'])
    y_butyrate.append(fluxes['butyrate'])
    y_total.append(fluxes['acetate'] + fluxes['propionate'] + fluxes['butyrate'])

X_data = np.array(X_data)
y_butyrate = np.array(y_butyrate)
y_total = np.array(y_total)

# Add realistic noise (~10%)
noise_level = 0.10
y_butyrate_noisy = y_butyrate * (1 + np.random.normal(0, noise_level, n_samples))
y_total_noisy    = y_total    * (1 + np.random.normal(0, noise_level, n_samples))

# Feature names
feature_names = species + ['High Fiber', 'Western Diet', 'Mediterranean', 'Low Carb']

# Cross-validation for butyrate prediction
from sklearn.model_selection import KFold
kf = KFold(n_splits=5, shuffle=True, random_state=42)

rf_model = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1)
gbm_model = GradientBoostingRegressor(n_estimators=100, random_state=42, learning_rate=0.1)

rf_scores  = cross_val_score(rf_model, X_data, y_butyrate_noisy, cv=kf, scoring='r2')
gbm_scores = cross_val_score(gbm_model, X_data, y_butyrate_noisy, cv=kf, scoring='r2')

print(f"  RandomForest R² (5-fold CV): {rf_scores.mean():.4f} ± {rf_scores.std():.4f}")
print(f"  GBM R² (5-fold CV):          {gbm_scores.mean():.4f} ± {gbm_scores.std():.4f}")

# Fit final model and get feature importances
rf_model.fit(X_data, y_butyrate_noisy)
feature_importances = pd.Series(rf_model.feature_importances_, index=feature_names).sort_values(ascending=False)
print(f"  Top features: {', '.join(feature_importances.head(5).index.tolist())}")

# Total SCFA prediction
rf_total_scores = cross_val_score(rf_model, X_data, y_total_noisy, cv=kf, scoring='r2')
print(f"  RandomForest R² (Total SCFA): {rf_total_scores.mean():.4f} ± {rf_total_scores.std():.4f}")

print("  ML prediction models complete.")

# ============================================================
# FIGURE GENERATION
# ============================================================
print("\n[FIGURES] Generating publication-quality figures...")

# Set style
plt.style.use('default')
sns.set_palette("husl")
FIGSIZE = (12, 8)
DPI = 150

# --- Figure 1: SHIME Dynamics ---
fig1, axes = plt.subplots(2, 2, figsize=(14, 10))
fig1.suptitle('Figure 1: SHIME-inspired Digestion Dynamics and SCFA Production',
               fontsize=14, fontweight='bold')

colors_diet = {'High Fiber': '#2ecc71', 'Western Diet': '#e74c3c', 
               'Mediterranean': '#3498db', 'Low Carb': '#f39c12'}

for diet, sol in shime_results.items():
    color = colors_diet[diet]
    axes[0, 0].plot(sol.t, sol.y[0], label=diet, color=color, linewidth=2)
    axes[0, 1].plot(sol.t, sol.y[4], label=f'{diet} (Ac)', color=color, linewidth=2, linestyle='-')
    axes[0, 1].plot(sol.t, sol.y[5], color=color, linewidth=1.5, linestyle='--')
    axes[0, 1].plot(sol.t, sol.y[6], color=color, linewidth=1.5, linestyle=':')
    axes[1, 0].plot(sol.t, sol.y[2], label=diet, color=color, linewidth=2)
    axes[1, 1].plot(sol.t, sol.y[7], label=diet, color=color, linewidth=2)

axes[0, 0].set_xlabel('Time (h)'); axes[0, 0].set_ylabel('Substrate (mM)')
axes[0, 0].set_title('Gastric Substrate Kinetics'); axes[0, 0].legend(fontsize=8)
axes[0, 1].set_xlabel('Time (h)'); axes[0, 1].set_ylabel('Concentration (mM)')
axes[0, 1].set_title('SCFA Production (solid=Acetate, dashed=Propionate, dot=Butyrate)')
axes[1, 0].set_xlabel('Time (h)'); axes[1, 0].set_ylabel('Substrate (mM)')
axes[1, 0].set_title('Proximal Colon Substrate'); axes[1, 0].legend(fontsize=8)
axes[1, 1].set_xlabel('Time (h)'); axes[1, 1].set_ylabel('pH')
axes[1, 1].set_title('Luminal pH Dynamics'); axes[1, 1].legend(fontsize=8)
axes[1, 1].set_ylim(5.5, 7.5)

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig1_shime_dynamics.png', dpi=DPI, bbox_inches='tight')
plt.close()
print("  Saved: fig1_shime_dynamics.png")

# --- Figure 2: gLV Community Dynamics ---
fig2, axes2 = plt.subplots(2, 2, figsize=(14, 10))
fig2.suptitle('Figure 2: gLV Microbiota Community Dynamics (60 days)',
               fontsize=14, fontweight='bold')

diet_list = list(gLV_results.keys())
for idx, diet in enumerate(diet_list):
    ax = axes2[idx // 2, idx % 2]
    res = gLV_results[diet]
    for i, sp in enumerate(species):
        ax.plot(res['time'], res['rel_abund'][i, :], label=sp, linewidth=1.5)
    ax.set_title(diet)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Relative Abundance')
    ax.legend(fontsize=6, loc='right')
    ax.set_ylim(0, 0.55)

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig2_glv_dynamics.png', dpi=DPI, bbox_inches='tight')
plt.close()
print("  Saved: fig2_glv_dynamics.png")

# --- Figure 3: SCFA Flux Comparison ---
fig3, axes3 = plt.subplots(1, 2, figsize=(14, 6))
fig3.suptitle('Figure 3: Community-level SCFA Flux Prediction', fontsize=14, fontweight='bold')

diets_list = list(fba_results.keys())
scfa_types = ['acetate', 'propionate', 'butyrate']
scfa_colors = ['#3498db', '#e67e22', '#2ecc71']

x_pos = np.arange(len(diets_list))
width = 0.25

for j, (scfa, color) in enumerate(zip(scfa_types, scfa_colors)):
    values = [fba_results[d][scfa] for d in diets_list]
    axes3[0].bar(x_pos + j * width, values, width, label=scfa.capitalize(), color=color, alpha=0.8)

axes3[0].set_xlabel('Diet Type')
axes3[0].set_ylabel('Flux (mmol/day)')
axes3[0].set_title('SCFA Production by Diet')
axes3[0].set_xticks(x_pos + width)
axes3[0].set_xticklabels(diets_list, rotation=15, ha='right')
axes3[0].legend()

# Pie charts for SCFA composition
for idx_d, diet in enumerate(diets_list[:4]):
    vals = [fba_results[diet][s] for s in scfa_types]
    axes3[1].plot([], [])

# Stacked bar
vals_stacked = np.array([[fba_results[d][s] for s in scfa_types] for d in diets_list])
bars1 = axes3[1].bar(x_pos, vals_stacked[:, 0], color=scfa_colors[0], alpha=0.8, label='Acetate')
bars2 = axes3[1].bar(x_pos, vals_stacked[:, 1], bottom=vals_stacked[:, 0], 
                      color=scfa_colors[1], alpha=0.8, label='Propionate')
bars3 = axes3[1].bar(x_pos, vals_stacked[:, 2], 
                      bottom=vals_stacked[:, 0] + vals_stacked[:, 1],
                      color=scfa_colors[2], alpha=0.8, label='Butyrate')
axes3[1].set_xticks(x_pos)
axes3[1].set_xticklabels(diets_list, rotation=15, ha='right')
axes3[1].set_ylabel('Total SCFA Flux (mmol/day)')
axes3[1].set_title('SCFA Composition by Diet')
axes3[1].legend()

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig3_scfa_flux.png', dpi=DPI, bbox_inches='tight')
plt.close()
print("  Saved: fig3_scfa_flux.png")

# --- Figure 4: Diet Transition Dynamics ---
fig4, axes4 = plt.subplots(1, 2, figsize=(14, 6))
fig4.suptitle('Figure 4: Long-term Diet Transition Dynamics', fontsize=14, fontweight='bold')

axes4[0].plot(t_WtoHF, shannon_WtoHF, 'g-', linewidth=2, label='Western → High Fiber')
axes4[0].plot(t_HFtoW, shannon_HFtoW, 'r-', linewidth=2, label='High Fiber → Western')
axes4[0].axvline(x=30, color='k', linestyle='--', alpha=0.5, label='Diet Change (day 30)')
axes4[0].set_xlabel('Time (days)')
axes4[0].set_ylabel("Shannon Diversity (H')")
axes4[0].set_title('Diversity Dynamics During Diet Transition')
axes4[0].legend()
axes4[0].grid(True, alpha=0.3)

# Species composition at key timepoints for WtoHF
n_timepoints = len(t_WtoHF)
tidx_0  = 0
tidx_30 = np.searchsorted(t_WtoHF, 30)
tidx_60 = np.searchsorted(t_WtoHF, 60)
tidx_end = -1

timepoint_data = np.array([
    comp_WtoHF[:, tidx_0],
    comp_WtoHF[:, tidx_30],
    comp_WtoHF[:, tidx_60],
    comp_WtoHF[:, tidx_end]
])

im = axes4[1].imshow(timepoint_data.T, aspect='auto', cmap='YlOrRd', 
                      vmin=0, vmax=0.4)
axes4[1].set_xticks([0, 1, 2, 3])
axes4[1].set_xticklabels(['Day 0', 'Day 30\n(switch)', 'Day 60', 'Day 120'])
axes4[1].set_yticks(range(n_species))
axes4[1].set_yticklabels(species, fontsize=9)
axes4[1].set_title('Community Composition Over Time\n(Western → High Fiber)')
plt.colorbar(im, ax=axes4[1], label='Relative Abundance')

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig4_transition_dynamics.png', dpi=DPI, bbox_inches='tight')
plt.close()
print("  Saved: fig4_transition_dynamics.png")

# --- Figure 5: Probiotic/Prebiotic Effects ---
fig5, axes5 = plt.subplots(2, 2, figsize=(14, 10))
fig5.suptitle('Figure 5: Probiotic and Prebiotic Intervention Effects', fontsize=14, fontweight='bold')

# Panel A: Shannon diversity over time
for name, res in probiotic_results.items():
    sh_time = [shannon_entropy(res['comp'][:, i]) for i in range(res['comp'].shape[1])]
    axes5[0, 0].plot(res['t'], sh_time, linewidth=2, label=name)
baseline_sh = [shannon_entropy(comp_baseline[:, i]) for i in range(comp_baseline.shape[1])]
axes5[0, 0].plot(t_pre, baseline_sh[:len(t_pre)], 'k--', linewidth=2, label='No intervention')
axes5[0, 0].set_xlabel('Days')
axes5[0, 0].set_ylabel("Shannon Diversity (H')")
axes5[0, 0].set_title('Diversity Response to Interventions')
axes5[0, 0].legend(fontsize=8)
axes5[0, 0].axvline(x=30, color='gray', linestyle=':', alpha=0.5)

# Panel B: Final composition comparison
interv_names = list(probiotic_results.keys()) + ['No intervention']
bifidobact_abund = []
lactob_abund = []
faecalib_abund = []
for name in interv_names[:-1]:
    bifidobact_abund.append(probiotic_results[name]['comp'][4, -1])
    lactob_abund.append(probiotic_results[name]['comp'][5, -1])
    faecalib_abund.append(probiotic_results[name]['comp'][3, -1])
bifidobact_abund.append(comp_baseline[4, -1])
lactob_abund.append(comp_baseline[5, -1])
faecalib_abund.append(comp_baseline[3, -1])

x_interv = np.arange(len(interv_names))
w = 0.25
axes5[0, 1].bar(x_interv - w, bifidobact_abund, w, label='Bifidobacterium', color='#9b59b6')
axes5[0, 1].bar(x_interv,     lactob_abund,     w, label='Lactobacillus',   color='#3498db')
axes5[0, 1].bar(x_interv + w, faecalib_abund,   w, label='Faecalibacterium', color='#2ecc71')
axes5[0, 1].set_xticks(x_interv)
axes5[0, 1].set_xticklabels(interv_names, rotation=20, ha='right', fontsize=8)
axes5[0, 1].set_ylabel('Relative Abundance')
axes5[0, 1].set_title('Key Taxa After Intervention (Day 60)')
axes5[0, 1].legend(fontsize=8)

# Panel C: SCFA production comparison
int_names_short = list(probiotic_results.keys()) + ['No intervention']
but_vals = []
for name in int_names_short[:-1]:
    comp_final = dict(zip(species, probiotic_results[name]['comp'][:, -1]))
    fluxes_int = predict_scfa_flux(comp_final, 'Western Diet', 45.0)
    but_vals.append(fluxes_int['butyrate'])
comp_base_final = dict(zip(species, comp_baseline[:, -1]))
fluxes_base = predict_scfa_flux(comp_base_final, 'Western Diet', 45.0)
but_vals.append(fluxes_base['butyrate'])

axes5[1, 0].bar(x_interv, but_vals, 
               color=['#e74c3c', '#e74c3c', '#27ae60', '#27ae60', '#95a5a6'],
               alpha=0.8)
axes5[1, 0].set_xticks(x_interv)
axes5[1, 0].set_xticklabels(int_names_short, rotation=20, ha='right', fontsize=8)
axes5[1, 0].set_ylabel('Butyrate Flux (mmol/day)')
axes5[1, 0].set_title('Butyrate Production by Intervention')
axes5[1, 0].axhline(y=but_vals[-1], color='k', linestyle='--', alpha=0.5, label='Baseline')

# Panel D: Feature importance
feat_imp_top = feature_importances.head(10)
axes5[1, 1].barh(feat_imp_top.index[::-1], feat_imp_top.values[::-1], 
                 color='steelblue', alpha=0.8)
axes5[1, 1].set_xlabel('Feature Importance')
axes5[1, 1].set_title('RF Feature Importance for Butyrate Prediction')

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig5_probiotic_prebiotic.png', dpi=DPI, bbox_inches='tight')
plt.close()
print("  Saved: fig5_probiotic_prebiotic.png")

# --- Figure 6: Fermented Food Case Study ---
fig6, axes6 = plt.subplots(1, 3, figsize=(16, 6))
fig6.suptitle('Figure 6: Fermented Food Intake – Microbiome Diversity and SCFA Effects',
               fontsize=14, fontweight='bold')

ff_names  = list(ff_results.keys())
sh_values = [ff_results[f]['shannon'] for f in ff_names]
but_values = [ff_results[f]['scfa']['butyrate'] for f in ff_names]
total_scfa = [ff_results[f]['scfa']['acetate'] + ff_results[f]['scfa']['propionate'] + 
              ff_results[f]['scfa']['butyrate'] for f in ff_names]

bar_colors = ['#e74c3c' if 'Baseline' in f else '#27ae60' for f in ff_names]

axes6[0].bar(ff_names, sh_values, color=bar_colors, alpha=0.85)
axes6[0].set_xticklabels(ff_names, rotation=30, ha='right', fontsize=9)
axes6[0].set_ylabel("Shannon Diversity (H')")
axes6[0].set_title('Microbiome Diversity')
axes6[0].axhline(y=sh_values[-1], color='r', linestyle='--', alpha=0.5)

axes6[1].bar(ff_names, but_values, color=bar_colors, alpha=0.85)
axes6[1].set_xticklabels(ff_names, rotation=30, ha='right', fontsize=9)
axes6[1].set_ylabel('Butyrate Flux (mmol/day)')
axes6[1].set_title('Butyrate Production')
axes6[1].axhline(y=but_values[-1], color='r', linestyle='--', alpha=0.5)

# Composition radar-ish for fermented foods
ff_names_ferm = [f for f in ff_names if 'Baseline' not in f]
ff_compositions = {f: dict(zip(species, ff_results[f]['comp'][:, -1])) 
                   for f in ff_names_ferm 
                   if 'comp' in ff_results[f]}
# Heatmap of composition
comp_matrix = np.array([[ff_results[f]['comp'][i, -1] 
                         for i in range(n_species)] 
                        for f in ff_names_ferm])
im6 = axes6[2].imshow(comp_matrix, aspect='auto', cmap='Blues', vmin=0, vmax=0.35)
axes6[2].set_xticks(range(n_species))
axes6[2].set_xticklabels(species, rotation=45, ha='right', fontsize=8)
axes6[2].set_yticks(range(len(ff_names_ferm)))
axes6[2].set_yticklabels(ff_names_ferm)
axes6[2].set_title('Microbiome Composition\n(Fermented Food Consumers)')
plt.colorbar(im6, ax=axes6[2], label='Relative Abundance')

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig6_fermented_food.png', dpi=DPI, bbox_inches='tight')
plt.close()
print("  Saved: fig6_fermented_food.png")

# ============================================================
# SAVE RESULTS DATA
# ============================================================
print("\n[DATA] Saving tabulated results...")

# Table 1: SHIME Summary
df_shime = pd.DataFrame(shime_summary).T.reset_index().rename(columns={'index': 'Diet'})
df_shime = df_shime.round(3)
df_shime.to_csv(f'{DATA_DIR}/shime_summary.csv', index=False)
print("  Saved: shime_summary.csv")

# Table 2: Diversity Metrics
div_df = pd.DataFrame({
    'Diet': list(diversity_results.keys()),
    'Shannon': [diversity_results[d]['shannon'] for d in diversity_results],
    'Simpson': [diversity_results[d]['simpson'] for d in diversity_results]
}).round(4)
div_df.to_csv(f'{DATA_DIR}/diversity_metrics.csv', index=False)
print("  Saved: diversity_metrics.csv")

# Table 3: SCFA Flux
scfa_df = pd.DataFrame(fba_results).T.reset_index().rename(columns={'index': 'Diet'})
scfa_df = scfa_df.round(3)
scfa_df.to_csv(f'{DATA_DIR}/scfa_flux.csv', index=False)
print("  Saved: scfa_flux.csv")

# Table 4: ML Results
ml_df = pd.DataFrame({
    'Model': ['Random Forest', 'Gradient Boosting'],
    'R2_mean': [rf_scores.mean(), gbm_scores.mean()],
    'R2_std':  [rf_scores.std(),  gbm_scores.std()],
    'Metric':  ['5-fold CV R²', '5-fold CV R²']
}).round(4)
ml_df.to_csv(f'{DATA_DIR}/ml_results.csv', index=False)

# Table 5: Fermented food
ff_df = pd.DataFrame({
    'Food': ff_names,
    'Shannon': sh_values,
    'Butyrate_mmol_day': but_values,
    'Total_SCFA_mmol_day': total_scfa
}).round(4)
ff_df.to_csv(f'{DATA_DIR}/fermented_food_results.csv', index=False)
print("  Saved: fermented_food_results.csv")
print("  Saved: ml_results.csv")

# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("QUANTITATIVE RESULTS SUMMARY")
print("=" * 60)
print("\n[SHIME - SCFA at 24h]")
for diet, vals in shime_summary.items():
    print(f"  {diet}: Acetate={vals['acetate_mM']:.2f}, "
          f"Propionate={vals['propionate_mM']:.2f}, "
          f"Butyrate={vals['butyrate_mM']:.2f} mM, pH={vals['final_pH']:.2f}")

print("\n[gLV Diversity - 60 day steady state]")
for diet, vals in diversity_results.items():
    print(f"  {diet}: Shannon={vals['shannon']:.4f}, Simpson={vals['simpson']:.4f}")

print("\n[FBA SCFA Flux - mmol/day]")
for diet, vals in fba_results.items():
    print(f"  {diet}: Butyrate={vals['butyrate']:.3f}, Total={vals['total_SCFA']:.3f}")

print("\n[ML Models - 5-fold CV R²]")
print(f"  RandomForest:     {rf_scores.mean():.4f} ± {rf_scores.std():.4f}")
print(f"  GradientBoosting: {gbm_scores.mean():.4f} ± {gbm_scores.std():.4f}")

print("\n[Fermented Food - Shannon & Butyrate]")
for food in ff_names:
    print(f"  {food}: Shannon={ff_results[food]['shannon']:.4f}, "
          f"Butyrate={ff_results[food]['scfa']['butyrate']:.3f}")

print("\n[Diet Transition]")
print(f"  Western Diet Shannon:    {diversity_results['Western Diet']['shannon']:.4f}")
print(f"  High Fiber Shannon:      {diversity_results['High Fiber']['shannon']:.4f}")
print(f"  Recovery time (W→HF):   {recovery_time_WtoHF:.1f} days")
print(f"  Final Shannon (W→HF):   {shannon_WtoHF[-1]:.4f}")

print("\n[Top ML Features for Butyrate Prediction]")
for feat, imp in feature_importances.head(5).items():
    print(f"  {feat}: {imp:.4f}")

print("\n[Cross-validation scores detail]")
print(f"  RF folds: {rf_scores}")
print(f"  GBM folds: {gbm_scores}")

print("\nAll figures and data saved successfully!")
