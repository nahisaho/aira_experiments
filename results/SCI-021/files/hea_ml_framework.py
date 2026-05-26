#!/usr/bin/env python3
"""
Machine Learning Framework for High-Entropy Alloy (HEA) Composition Optimization
- CALPHAD-inspired thermodynamic database
- Descriptor design (atomic radius difference, VEC, mixing entropy)
- Multi-objective Bayesian optimization (strength, ductility, corrosion resistance)
- DFT-inspired data generation
- Active learning loop
- CrMnFeCoNi case study
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, RBF, ConstantKernel
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, r2_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import combinations
import warnings
import json
import os

warnings.filterwarnings('ignore')
np.random.seed(42)

ELEMENTS = ['Cr', 'Mn', 'Fe', 'Co', 'Ni']

# Elemental properties database
ELEMENT_PROPS = {
    'Cr': {'r': 1.249, 'VEC': 6, 'chi': 1.66, 'Tm': 2180, 'E': 279, 'mass': 52.0},
    'Mn': {'r': 1.350, 'VEC': 7, 'chi': 1.55, 'Tm': 1519, 'E': 198, 'mass': 54.9},
    'Fe': {'r': 1.241, 'VEC': 8, 'chi': 1.83, 'Tm': 1811, 'E': 211, 'mass': 55.8},
    'Co': {'r': 1.251, 'VEC': 9, 'chi': 1.88, 'Tm': 1768, 'E': 209, 'mass': 58.9},
    'Ni': {'r': 1.246, 'VEC': 10, 'chi': 1.91, 'Tm': 1728, 'E': 200, 'mass': 58.7},
}

# Binary mixing enthalpy (kJ/mol) - simplified Miedema model
BINARY_ENTHALPY = {
    ('Cr','Mn'): 2, ('Cr','Fe'): -1, ('Cr','Co'): -4, ('Cr','Ni'): -7,
    ('Mn','Fe'): 0, ('Mn','Co'): -5, ('Mn','Ni'): -8,
    ('Fe','Co'): -1, ('Fe','Ni'): -2,
    ('Co','Ni'): 0,
}

R_GAS = 8.314  # J/(mol·K)

# ============================================================================
# 1. CALPHAD-inspired Thermodynamic Calculations
# ============================================================================

def calc_mixing_entropy(x):
    """Configurational mixing entropy ΔS_mix = -R Σ x_i ln(x_i)"""
    x = np.array(x)
    x = x[x > 1e-10]
    return -R_GAS * np.sum(x * np.log(x))

def calc_mixing_enthalpy(x):
    """Mixing enthalpy ΔH_mix = Σ 4·ΔH_ij·x_i·x_j"""
    h_mix = 0
    for i in range(len(ELEMENTS)):
        for j in range(i+1, len(ELEMENTS)):
            key = (ELEMENTS[i], ELEMENTS[j])
            rev_key = (ELEMENTS[j], ELEMENTS[i])
            h_ij = BINARY_ENTHALPY.get(key, BINARY_ENTHALPY.get(rev_key, 0))
            h_mix += 4 * h_ij * x[i] * x[j]
    return h_mix

def calc_atomic_size_diff(x):
    """Atomic size difference δ = √(Σ x_i(1 - r_i/r_avg)²)"""
    r = np.array([ELEMENT_PROPS[e]['r'] for e in ELEMENTS])
    r_avg = np.sum(x * r)
    delta = np.sqrt(np.sum(x * (1 - r / r_avg)**2))
    return delta * 100  # percentage

def calc_vec(x):
    """Valence Electron Concentration VEC = Σ x_i·VEC_i"""
    vec = np.array([ELEMENT_PROPS[e]['VEC'] for e in ELEMENTS])
    return np.sum(x * vec)

def calc_electronegativity_diff(x):
    """Electronegativity difference Δχ = √(Σ x_i(χ_i - χ_avg)²)"""
    chi = np.array([ELEMENT_PROPS[e]['chi'] for e in ELEMENTS])
    chi_avg = np.sum(x * chi)
    return np.sqrt(np.sum(x * (chi - chi_avg)**2))

def calc_omega(x, T=1500):
    """Ω parameter: T·ΔS_mix / |ΔH_mix|"""
    s_mix = calc_mixing_entropy(x)
    h_mix = calc_mixing_enthalpy(x)
    if abs(h_mix) < 1e-10:
        return 100.0
    return T * s_mix / (abs(h_mix) * 1000)

def calc_all_descriptors(x):
    """Calculate all composition-based descriptors"""
    return {
        'delta_S_mix': calc_mixing_entropy(x),
        'delta_H_mix': calc_mixing_enthalpy(x),
        'delta_r': calc_atomic_size_diff(x),
        'VEC': calc_vec(x),
        'delta_chi': calc_electronegativity_diff(x),
        'Omega': calc_omega(x),
        'T_avg': np.sum(x * np.array([ELEMENT_PROPS[e]['Tm'] for e in ELEMENTS])),
        'E_avg': np.sum(x * np.array([ELEMENT_PROPS[e]['E'] for e in ELEMENTS])),
    }

# ============================================================================
# 2. CALPHAD Phase Diagram Calculation (Simplified)
# ============================================================================

def gibbs_energy_fcc(x, T):
    """Simplified Gibbs energy for FCC phase"""
    h_mix = calc_mixing_enthalpy(x) * 1000
    s_mix = calc_mixing_entropy(x)
    g_ref = np.sum(x * np.array([-20000, -15000, -18000, -17000, -19000]))
    return g_ref + h_mix - T * s_mix

def gibbs_energy_bcc(x, T):
    """Simplified Gibbs energy for BCC phase"""
    h_mix = calc_mixing_enthalpy(x) * 1000
    s_mix = calc_mixing_entropy(x)
    g_ref = np.sum(x * np.array([-18000, -14000, -17500, -18500, -20000]))
    return g_ref + h_mix * 1.1 - T * s_mix

def predict_phase(x, T=1500):
    """Predict stable phase based on Gibbs energy"""
    g_fcc = gibbs_energy_fcc(x, T)
    g_bcc = gibbs_energy_bcc(x, T)
    vec = calc_vec(x)
    if vec >= 8.0:
        return 'FCC', g_fcc
    elif vec <= 6.87:
        return 'BCC', g_bcc
    else:
        return 'FCC+BCC', min(g_fcc, g_bcc)

def calphad_phase_diagram(elem1_idx=0, elem2_idx=4, T_range=(800, 2200)):
    """Generate pseudo-binary phase diagram section"""
    T_vals = np.linspace(T_range[0], T_range[1], 100)
    x_vals = np.linspace(0.05, 0.45, 50)
    phases = np.zeros((len(T_vals), len(x_vals)))
    
    for i, T in enumerate(T_vals):
        for j, x1 in enumerate(x_vals):
            comp = np.array([0.2]*5)
            comp[elem1_idx] = x1
            remaining = 1.0 - x1
            for k in range(5):
                if k != elem1_idx:
                    comp[k] = remaining / 4
            phase, _ = predict_phase(comp, T)
            if phase == 'FCC':
                phases[i, j] = 0
            elif phase == 'BCC':
                phases[i, j] = 1
            else:
                phases[i, j] = 0.5
    return x_vals, T_vals, phases

# ============================================================================
# 3. DFT-inspired Data Generation
# ============================================================================

def simulate_dft_properties(x):
    """
    Simulate DFT-calculated properties with physically motivated models.
    Returns formation energy, lattice parameter, bulk modulus.
    """
    desc = calc_all_descriptors(x)
    
    # Formation energy (eV/atom) - negative = stable
    e_form = (desc['delta_H_mix'] * 0.01 
              + 0.5 * desc['delta_r']**2 
              - 0.02 * desc['delta_S_mix']
              + np.random.normal(0, 0.005))
    
    # Lattice parameter (Å) - weighted average with strain correction
    a0 = 2 * np.sum(x * np.array([ELEMENT_PROPS[e]['r'] for e in ELEMENTS]))
    a0 += desc['delta_r'] * 0.01 + np.random.normal(0, 0.002)
    
    # Bulk modulus (GPa)
    B = desc['E_avg'] * 0.8 + desc['VEC'] * 5 - desc['delta_r'] * 50
    B += np.random.normal(0, 2)
    
    return {
        'E_form': e_form,
        'a0': a0,
        'B': max(B, 50),
    }

def generate_property_targets(x):
    """
    Generate target properties for optimization:
    - Yield strength (MPa)
    - Elongation (%)
    - Corrosion resistance index
    """
    desc = calc_all_descriptors(x)
    dft = simulate_dft_properties(x)
    
    # Yield strength model (solid solution strengthening)
    sigma_ss = 200 + 1500 * desc['delta_r'] + 30 * abs(desc['delta_H_mix'])
    sigma_gb = 100 * (desc['VEC'] - 7)**2
    strength = sigma_ss + sigma_gb + np.random.normal(0, 15)
    strength = max(strength, 150)
    
    # Elongation model (inversely related to strength, VEC dependent)
    elong = 60 - 0.03 * strength + 5 * (desc['VEC'] - 8) + 10 * desc['delta_S_mix'] / R_GAS
    elong += np.random.normal(0, 3)
    elong = np.clip(elong, 2, 80)
    
    # Corrosion resistance (Cr content and passivation)
    cr_idx = ELEMENTS.index('Cr')
    corr = 50 * x[cr_idx] + 20 * desc['delta_chi'] + 10 * (desc['Omega'] if desc['Omega'] < 10 else 10)
    corr += np.random.normal(0, 2)
    corr = np.clip(corr, 0, 100)
    
    return {
        'strength': strength,
        'elongation': elong,
        'corrosion_resistance': corr,
    }

# ============================================================================
# 4. Dataset Generation
# ============================================================================

def generate_random_composition(n_elements=5, min_frac=0.05):
    """Generate random HEA composition with minimum fraction constraint"""
    while True:
        x = np.random.dirichlet(np.ones(n_elements))
        if np.all(x >= min_frac):
            return x

def build_dataset(n_samples=500):
    """Build comprehensive HEA dataset"""
    data = []
    for _ in range(n_samples):
        x = generate_random_composition()
        desc = calc_all_descriptors(x)
        props = generate_property_targets(x)
        dft = simulate_dft_properties(x)
        phase, g = predict_phase(x)
        
        row = {}
        for i, elem in enumerate(ELEMENTS):
            row[f'x_{elem}'] = x[i]
        row.update(desc)
        row.update(props)
        row.update(dft)
        row['phase'] = phase
        row['G'] = g
        data.append(row)
    
    return pd.DataFrame(data)

# ============================================================================
# 5. ML Model Training
# ============================================================================

def train_surrogate_models(df):
    """Train GP and RF surrogate models for each target property"""
    feature_cols = ['delta_S_mix', 'delta_H_mix', 'delta_r', 'VEC',
                    'delta_chi', 'Omega', 'T_avg', 'E_avg']
    targets = ['strength', 'elongation', 'corrosion_resistance']
    
    X = df[feature_cols].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    models = {}
    scores = {}
    
    for target in targets:
        y = df[target].values
        
        # GP model
        kernel = ConstantKernel(1.0) * Matern(nu=2.5)
        gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5, alpha=0.1)
        gp.fit(X_scaled, y)
        
        # RF model
        rf = RandomForestRegressor(n_estimators=200, max_depth=15, random_state=42)
        rf.fit(X_scaled, y)
        
        # GB model
        gb = GradientBoostingRegressor(n_estimators=200, max_depth=5, random_state=42)
        gb.fit(X_scaled, y)
        
        # Cross-validation scores
        kf = KFold(n_splits=5, shuffle=True, random_state=42)
        cv_rf = cross_val_score(rf, X_scaled, y, cv=kf, scoring='r2')
        cv_gb = cross_val_score(gb, X_scaled, y, cv=kf, scoring='r2')
        
        models[target] = {'gp': gp, 'rf': rf, 'gb': gb}
        scores[target] = {
            'rf_r2': cv_rf.mean(),
            'rf_r2_std': cv_rf.std(),
            'gb_r2': cv_gb.mean(),
            'gb_r2_std': cv_gb.std(),
        }
    
    return models, scaler, scores

# ============================================================================
# 6. Multi-Objective Bayesian Optimization
# ============================================================================

def expected_improvement(X, model, y_best, xi=0.01):
    """Calculate Expected Improvement acquisition function"""
    mu, sigma = model.predict(X.reshape(1, -1), return_std=True)
    sigma = max(sigma[0], 1e-9)
    imp = mu[0] - y_best - xi
    Z = imp / sigma
    ei = imp * norm.cdf(Z) + sigma * norm.pdf(Z)
    return ei

def multi_objective_acquisition(x_comp, models, scaler, y_bests, weights=(0.4, 0.3, 0.3)):
    """Multi-objective acquisition combining EI for all targets"""
    desc = calc_all_descriptors(x_comp)
    feature_cols = ['delta_S_mix', 'delta_H_mix', 'delta_r', 'VEC',
                    'delta_chi', 'Omega', 'T_avg', 'E_avg']
    X = np.array([desc[f] for f in feature_cols])
    X_scaled = scaler.transform(X.reshape(1, -1))
    
    targets = ['strength', 'elongation', 'corrosion_resistance']
    total_ei = 0
    for i, target in enumerate(targets):
        ei = expected_improvement(X_scaled[0], models[target]['gp'], y_bests[target])
        total_ei += weights[i] * ei
    return total_ei

def bayesian_optimization_loop(df, models, scaler, n_iterations=30):
    """Run multi-objective Bayesian optimization"""
    targets = ['strength', 'elongation', 'corrosion_resistance']
    y_bests = {t: df[t].max() for t in targets}
    
    history = []
    best_compositions = []
    pareto_front = []
    
    for iteration in range(n_iterations):
        best_ei = -np.inf
        best_comp = None
        
        # Random search for best acquisition value
        for _ in range(200):
            x_cand = generate_random_composition()
            ei = multi_objective_acquisition(x_cand, models, scaler, y_bests)
            if ei > best_ei:
                best_ei = ei
                best_comp = x_cand
        
        # Evaluate candidate
        props = generate_property_targets(best_comp)
        desc = calc_all_descriptors(best_comp)
        
        # Update bests
        for t in targets:
            if props[t] > y_bests[t]:
                y_bests[t] = props[t]
        
        history.append({
            'iteration': iteration,
            'composition': best_comp.tolist(),
            'strength': props['strength'],
            'elongation': props['elongation'],
            'corrosion_resistance': props['corrosion_resistance'],
            'EI': best_ei,
        })
        
        # Check Pareto optimality
        dominated = False
        for h in pareto_front:
            if (h['strength'] >= props['strength'] and 
                h['elongation'] >= props['elongation'] and
                h['corrosion_resistance'] >= props['corrosion_resistance']):
                dominated = True
                break
        if not dominated:
            pareto_front = [h for h in pareto_front if not (
                props['strength'] >= h['strength'] and
                props['elongation'] >= h['elongation'] and
                props['corrosion_resistance'] >= h['corrosion_resistance']
            )]
            pareto_front.append({
                'composition': best_comp.tolist(),
                **props
            })
    
    return history, pareto_front

# ============================================================================
# 7. Active Learning Loop
# ============================================================================

def uncertainty_sampling(models, scaler, n_candidates=500):
    """Select next experiment based on model uncertainty"""
    candidates = []
    uncertainties = []
    
    for _ in range(n_candidates):
        x = generate_random_composition()
        desc = calc_all_descriptors(x)
        feature_cols = ['delta_S_mix', 'delta_H_mix', 'delta_r', 'VEC',
                        'delta_chi', 'Omega', 'T_avg', 'E_avg']
        X = np.array([desc[f] for f in feature_cols])
        X_scaled = scaler.transform(X.reshape(1, -1))
        
        total_unc = 0
        for target in ['strength', 'elongation', 'corrosion_resistance']:
            _, sigma = models[target]['gp'].predict(X_scaled, return_std=True)
            total_unc += sigma[0]
        
        candidates.append(x)
        uncertainties.append(total_unc)
    
    # Return top-k most uncertain
    top_idx = np.argsort(uncertainties)[-10:]
    return [candidates[i] for i in top_idx], [uncertainties[i] for i in top_idx]

def active_learning_loop(initial_df, models, scaler, n_rounds=5, samples_per_round=10):
    """Run active learning to efficiently expand the dataset"""
    df = initial_df.copy()
    learning_curves = {'round': [], 'n_samples': [], 'mae_strength': [],
                       'mae_elongation': [], 'mae_corrosion': []}
    
    feature_cols = ['delta_S_mix', 'delta_H_mix', 'delta_r', 'VEC',
                    'delta_chi', 'Omega', 'T_avg', 'E_avg']
    
    for round_i in range(n_rounds):
        # Evaluate current model performance
        X = df[feature_cols].values
        X_scaled = scaler.fit_transform(X)
        
        for target in ['strength', 'elongation', 'corrosion_resistance']:
            models[target]['rf'].fit(X_scaled, df[target].values)
        
        # Record MAE via cross-validation
        kf = KFold(n_splits=min(5, len(df)//10), shuffle=True, random_state=42)
        maes = {}
        for target in ['strength', 'elongation', 'corrosion_resistance']:
            cv = cross_val_score(models[target]['rf'], X_scaled, df[target].values,
                                cv=kf, scoring='neg_mean_absolute_error')
            maes[target] = -cv.mean()
        
        learning_curves['round'].append(round_i)
        learning_curves['n_samples'].append(len(df))
        learning_curves['mae_strength'].append(maes['strength'])
        learning_curves['mae_elongation'].append(maes['elongation'])
        learning_curves['mae_corrosion'].append(maes['corrosion_resistance'])
        
        # Select informative samples
        uncertain_comps, _ = uncertainty_sampling(models, scaler, n_candidates=300)
        
        # "Run experiments" on selected compositions
        new_rows = []
        for x in uncertain_comps[:samples_per_round]:
            desc = calc_all_descriptors(x)
            props = generate_property_targets(x)
            dft = simulate_dft_properties(x)
            phase, g = predict_phase(x)
            row = {}
            for i, elem in enumerate(ELEMENTS):
                row[f'x_{elem}'] = x[i]
            row.update(desc)
            row.update(props)
            row.update(dft)
            row['phase'] = phase
            row['G'] = g
            new_rows.append(row)
        
        df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)
    
    return df, learning_curves

# ============================================================================
# 8. Visualization Functions
# ============================================================================

def plot_phase_diagram(save_path='figures/phase_diagram.png'):
    """Plot pseudo-binary phase diagram"""
    x_vals, T_vals, phases = calphad_phase_diagram()
    
    fig, ax = plt.subplots(figsize=(8, 6))
    c = ax.pcolormesh(x_vals, T_vals, phases, cmap='RdYlBu_r', shading='auto')
    ax.set_xlabel(f'{ELEMENTS[0]} fraction (balance: equal others)', fontsize=12)
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('Pseudo-binary Phase Diagram: Cr-CrMnFeCoNi', fontsize=14)
    cbar = fig.colorbar(c, ax=ax, ticks=[0, 0.5, 1])
    cbar.ax.set_yticklabels(['FCC', 'FCC+BCC', 'BCC'])
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()

def plot_descriptor_distributions(df, save_path='figures/descriptor_distributions.png'):
    """Plot distributions of key descriptors"""
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    descriptors = ['delta_S_mix', 'delta_H_mix', 'delta_r', 'VEC', 'delta_chi', 'Omega']
    labels = ['ΔS_mix (J/mol·K)', 'ΔH_mix (kJ/mol)', 'δ (%)', 'VEC', 'Δχ', 'Ω']
    
    for ax, desc, label in zip(axes.flat, descriptors, labels):
        ax.hist(df[desc], bins=30, color='steelblue', alpha=0.7, edgecolor='black')
        ax.set_xlabel(label, fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.axvline(df[desc].mean(), color='red', linestyle='--', label=f'Mean={df[desc].mean():.2f}')
        ax.legend(fontsize=9)
    
    fig.suptitle('Distribution of Composition Descriptors', fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()

def plot_property_correlations(df, save_path='figures/property_correlations.png'):
    """Plot property correlations with descriptors"""
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    
    plots = [
        ('VEC', 'strength', 'VEC vs Yield Strength'),
        ('delta_r', 'strength', 'δ vs Yield Strength'),
        ('delta_S_mix', 'elongation', 'ΔS_mix vs Elongation'),
        ('VEC', 'elongation', 'VEC vs Elongation'),
        ('delta_chi', 'corrosion_resistance', 'Δχ vs Corrosion Resistance'),
        ('strength', 'elongation', 'Strength-Ductility Trade-off'),
    ]
    
    for ax, (x_col, y_col, title) in zip(axes.flat, plots):
        colors = df['corrosion_resistance'] if x_col != 'strength' else df['VEC']
        sc = ax.scatter(df[x_col], df[y_col], c=colors, cmap='viridis', alpha=0.6, s=15)
        ax.set_xlabel(x_col, fontsize=11)
        ax.set_ylabel(y_col, fontsize=11)
        ax.set_title(title, fontsize=11)
        plt.colorbar(sc, ax=ax)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()

def plot_model_performance(df, models, scaler, save_path='figures/model_performance.png'):
    """Plot ML model prediction vs actual"""
    feature_cols = ['delta_S_mix', 'delta_H_mix', 'delta_r', 'VEC',
                    'delta_chi', 'Omega', 'T_avg', 'E_avg']
    X = scaler.transform(df[feature_cols].values)
    targets = ['strength', 'elongation', 'corrosion_resistance']
    labels = ['Yield Strength (MPa)', 'Elongation (%)', 'Corrosion Resistance']
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    for ax, target, label in zip(axes, targets, labels):
        y_true = df[target].values
        y_pred_rf = models[target]['rf'].predict(X)
        y_pred_gb = models[target]['gb'].predict(X)
        
        r2_rf = r2_score(y_true, y_pred_rf)
        r2_gb = r2_score(y_true, y_pred_gb)
        
        ax.scatter(y_true, y_pred_rf, alpha=0.4, s=15, label=f'RF (R²={r2_rf:.3f})', c='blue')
        ax.scatter(y_true, y_pred_gb, alpha=0.4, s=15, label=f'GB (R²={r2_gb:.3f})', c='red')
        lims = [min(y_true.min(), y_pred_rf.min()), max(y_true.max(), y_pred_rf.max())]
        ax.plot(lims, lims, 'k--', alpha=0.5)
        ax.set_xlabel(f'Actual {label}', fontsize=11)
        ax.set_ylabel(f'Predicted {label}', fontsize=11)
        ax.legend(fontsize=9)
        ax.set_title(f'{label}', fontsize=12)
    
    plt.suptitle('ML Model Performance: Predicted vs Actual', fontsize=14)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()

def plot_bayesian_optimization(history, save_path='figures/bayesian_optimization.png'):
    """Plot Bayesian optimization convergence"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    iters = [h['iteration'] for h in history]
    
    # Strength convergence
    strengths = [h['strength'] for h in history]
    best_s = np.maximum.accumulate(strengths)
    axes[0,0].plot(iters, strengths, 'o-', alpha=0.5, label='Sampled')
    axes[0,0].plot(iters, best_s, 'r-', linewidth=2, label='Best found')
    axes[0,0].set_ylabel('Yield Strength (MPa)')
    axes[0,0].set_title('Strength Optimization')
    axes[0,0].legend()
    
    # Elongation convergence
    elongs = [h['elongation'] for h in history]
    best_e = np.maximum.accumulate(elongs)
    axes[0,1].plot(iters, elongs, 'o-', alpha=0.5, label='Sampled')
    axes[0,1].plot(iters, best_e, 'r-', linewidth=2, label='Best found')
    axes[0,1].set_ylabel('Elongation (%)')
    axes[0,1].set_title('Ductility Optimization')
    axes[0,1].legend()
    
    # Corrosion convergence
    corrs = [h['corrosion_resistance'] for h in history]
    best_c = np.maximum.accumulate(corrs)
    axes[1,0].plot(iters, corrs, 'o-', alpha=0.5, label='Sampled')
    axes[1,0].plot(iters, best_c, 'r-', linewidth=2, label='Best found')
    axes[1,0].set_ylabel('Corrosion Resistance')
    axes[1,0].set_xlabel('Iteration')
    axes[1,0].set_title('Corrosion Resistance Optimization')
    axes[1,0].legend()
    
    # EI convergence
    eis = [h['EI'] for h in history]
    axes[1,1].plot(iters, eis, 'g-o', alpha=0.7)
    axes[1,1].set_ylabel('Acquisition Value (EI)')
    axes[1,1].set_xlabel('Iteration')
    axes[1,1].set_title('Expected Improvement')
    
    plt.suptitle('Multi-Objective Bayesian Optimization Progress', fontsize=14)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()

def plot_pareto_front(pareto_front, save_path='figures/pareto_front.png'):
    """Plot 3D Pareto front"""
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    s = [p['strength'] for p in pareto_front]
    e = [p['elongation'] for p in pareto_front]
    c = [p['corrosion_resistance'] for p in pareto_front]
    
    sc = ax.scatter(s, e, c, c=c, cmap='plasma', s=80, edgecolors='black')
    ax.set_xlabel('Yield Strength (MPa)', fontsize=11)
    ax.set_ylabel('Elongation (%)', fontsize=11)
    ax.set_zlabel('Corrosion Resistance', fontsize=11)
    ax.set_title('Pareto Front: Multi-Objective Optimization', fontsize=13)
    fig.colorbar(sc, ax=ax, label='Corrosion Resistance', shrink=0.6)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()

def plot_active_learning(learning_curves, save_path='figures/active_learning.png'):
    """Plot active learning curves"""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    rounds = learning_curves['n_samples']
    ax.plot(rounds, learning_curves['mae_strength'], 's-', label='Strength MAE', color='blue')
    ax.plot(rounds, learning_curves['mae_elongation'], 'o-', label='Elongation MAE', color='green')
    ax.plot(rounds, learning_curves['mae_corrosion'], '^-', label='Corrosion Resistance MAE', color='red')
    
    ax.set_xlabel('Number of Training Samples', fontsize=12)
    ax.set_ylabel('Mean Absolute Error', fontsize=12)
    ax.set_title('Active Learning: Model Improvement with Data', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()

def plot_feature_importance(df, models, scaler, save_path='figures/feature_importance.png'):
    """Plot feature importance from RF models"""
    feature_cols = ['delta_S_mix', 'delta_H_mix', 'delta_r', 'VEC',
                    'delta_chi', 'Omega', 'T_avg', 'E_avg']
    feature_labels = ['ΔS_mix', 'ΔH_mix', 'δ', 'VEC', 'Δχ', 'Ω', 'T_avg', 'E_avg']
    targets = ['strength', 'elongation', 'corrosion_resistance']
    target_labels = ['Yield Strength', 'Elongation', 'Corrosion Resistance']
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    for ax, target, tlabel in zip(axes, targets, target_labels):
        importances = models[target]['rf'].feature_importances_
        idx = np.argsort(importances)
        ax.barh(np.array(feature_labels)[idx], importances[idx], color='steelblue')
        ax.set_xlabel('Feature Importance', fontsize=11)
        ax.set_title(tlabel, fontsize=12)
    
    plt.suptitle('Random Forest Feature Importance', fontsize=14)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()

def plot_composition_heatmap(pareto_front, save_path='figures/composition_heatmap.png'):
    """Plot optimal compositions as heatmap"""
    comps = np.array([p['composition'] for p in pareto_front])
    
    fig, ax = plt.subplots(figsize=(10, max(4, len(pareto_front)*0.4)))
    im = ax.imshow(comps, cmap='YlOrRd', aspect='auto', vmin=0, vmax=0.5)
    
    ax.set_xticks(range(len(ELEMENTS)))
    ax.set_xticklabels(ELEMENTS, fontsize=12)
    ax.set_ylabel('Pareto Solution Index', fontsize=12)
    ax.set_title('Optimal Compositions on Pareto Front', fontsize=14)
    
    for i in range(comps.shape[0]):
        for j in range(comps.shape[1]):
            ax.text(j, i, f'{comps[i,j]:.2f}', ha='center', va='center', fontsize=9)
    
    fig.colorbar(im, ax=ax, label='Atomic Fraction')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()

def plot_case_study(save_path='figures/case_study_crmnfeconi.png'):
    """CrMnFeCoNi case study: vary compositions around equiatomic"""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Vary each element around equiatomic
    fracs = np.linspace(0.05, 0.45, 40)
    
    for idx, (elem, ax_row) in enumerate(zip(ELEMENTS[:3], [axes[0,0], axes[0,1], axes[0,2]])):
        strengths, elongs, corrs = [], [], []
        for f in fracs:
            x = np.array([0.2]*5)
            x[idx] = f
            remaining = 1.0 - f
            for k in range(5):
                if k != idx:
                    x[k] = remaining / 4
            props = generate_property_targets(x)
            strengths.append(props['strength'])
            elongs.append(props['elongation'])
            corrs.append(props['corrosion_resistance'])
        
        ax_row.plot(fracs, strengths, 'b-', label='Strength')
        ax2 = ax_row.twinx()
        ax2.plot(fracs, elongs, 'r--', label='Elongation')
        ax_row.set_xlabel(f'{elem} fraction')
        ax_row.set_ylabel('Strength (MPa)', color='b')
        ax2.set_ylabel('Elongation (%)', color='r')
        ax_row.set_title(f'Effect of {elem} content')
        ax_row.axvline(0.2, color='gray', linestyle=':', alpha=0.5)
    
    for idx, (elem, ax_row) in enumerate(zip(ELEMENTS[3:], [axes[1,0], axes[1,1]])):
        strengths, elongs, corrs = [], [], []
        for f in fracs:
            x = np.array([0.2]*5)
            x[idx+3] = f
            remaining = 1.0 - f
            for k in range(5):
                if k != idx+3:
                    x[k] = remaining / 4
            props = generate_property_targets(x)
            strengths.append(props['strength'])
            elongs.append(props['elongation'])
            corrs.append(props['corrosion_resistance'])
        
        ax_row.plot(fracs, strengths, 'b-', label='Strength')
        ax2 = ax_row.twinx()
        ax2.plot(fracs, corrs, 'g--', label='Corrosion Res.')
        ax_row.set_xlabel(f'{elem} fraction')
        ax_row.set_ylabel('Strength (MPa)', color='b')
        ax2.set_ylabel('Corrosion Resistance', color='g')
        ax_row.set_title(f'Effect of {elem} content')
        ax_row.axvline(0.2, color='gray', linestyle=':', alpha=0.5)
    
    # Phase stability map
    ax = axes[1,2]
    vec_vals, delta_vals, phases_str = [], [], []
    for _ in range(300):
        x = generate_random_composition()
        vec_vals.append(calc_vec(x))
        delta_vals.append(calc_atomic_size_diff(x))
        phase, _ = predict_phase(x)
        phases_str.append(0 if phase=='FCC' else (1 if phase=='BCC' else 0.5))
    
    sc = ax.scatter(vec_vals, delta_vals, c=phases_str, cmap='RdYlBu_r', s=20, alpha=0.7)
    ax.set_xlabel('VEC')
    ax.set_ylabel('δ (%)')
    ax.set_title('Phase Stability Map')
    plt.colorbar(sc, ax=ax, ticks=[0, 0.5, 1]).ax.set_yticklabels(['FCC', 'FCC+BCC', 'BCC'])
    
    plt.suptitle('CrMnFeCoNi System: Composition-Property Relationships', fontsize=14)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()

# ============================================================================
# 9. Main Pipeline
# ============================================================================

def main():
    print("=" * 70)
    print("HEA-ML: Machine Learning Framework for High-Entropy Alloy Design")
    print("=" * 70)
    
    # Step 1: Generate dataset
    print("\n[1/7] Generating HEA dataset with CALPHAD-inspired thermodynamics...")
    df = build_dataset(n_samples=500)
    df.to_csv('hea_dataset.csv', index=False)
    print(f"  Dataset: {len(df)} compositions, {len(df.columns)} features")
    print(f"  Phases: {df['phase'].value_counts().to_dict()}")
    
    # Step 2: Plot phase diagram
    print("\n[2/7] Computing phase diagram...")
    plot_phase_diagram()
    plot_descriptor_distributions(df)
    print("  Phase diagram and descriptor distributions saved.")
    
    # Step 3: Train ML models
    print("\n[3/7] Training surrogate models (GP, RF, GB)...")
    models, scaler, scores = train_surrogate_models(df)
    print("  Cross-validation R² scores:")
    for target, s in scores.items():
        print(f"    {target}: RF={s['rf_r2']:.3f}±{s['rf_r2_std']:.3f}, "
              f"GB={s['gb_r2']:.3f}±{s['gb_r2_std']:.3f}")
    
    # Step 4: Plot model performance and feature importance
    print("\n[4/7] Evaluating model performance...")
    plot_property_correlations(df)
    plot_model_performance(df, models, scaler)
    plot_feature_importance(df, models, scaler)
    print("  Model evaluation plots saved.")
    
    # Step 5: Bayesian optimization
    print("\n[5/7] Running multi-objective Bayesian optimization...")
    history, pareto_front = bayesian_optimization_loop(df, models, scaler, n_iterations=30)
    plot_bayesian_optimization(history)
    plot_pareto_front(pareto_front)
    plot_composition_heatmap(pareto_front)
    print(f"  Pareto front: {len(pareto_front)} non-dominated solutions")
    
    # Print top solutions
    print("\n  Top Pareto-optimal compositions:")
    for i, p in enumerate(pareto_front[:5]):
        comp_str = ', '.join([f"{ELEMENTS[j]}={p['composition'][j]:.3f}" for j in range(5)])
        print(f"    #{i+1}: [{comp_str}]")
        print(f"         σ={p['strength']:.0f} MPa, ε={p['elongation']:.1f}%, CR={p['corrosion_resistance']:.1f}")
    
    # Step 6: Active learning
    print("\n[6/7] Running active learning loop...")
    df_expanded, learning_curves = active_learning_loop(df, models, scaler, n_rounds=5, samples_per_round=10)
    plot_active_learning(learning_curves)
    df_expanded.to_csv('hea_dataset_expanded.csv', index=False)
    print(f"  Dataset expanded: {len(df)} → {len(df_expanded)} samples")
    print("  Learning curve saved.")
    
    # Step 7: Case study
    print("\n[7/7] Running CrMnFeCoNi case study...")
    plot_case_study()
    print("  Case study plots saved.")
    
    # Save results summary
    results = {
        'dataset_size': len(df),
        'expanded_dataset_size': len(df_expanded),
        'model_scores': scores,
        'n_pareto_solutions': len(pareto_front),
        'pareto_front': pareto_front,
        'bo_history': history,
        'learning_curves': learning_curves,
    }
    
    with open('experiment_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print("\n" + "=" * 70)
    print("All experiments completed successfully!")
    print("Generated files:")
    print("  - hea_dataset.csv")
    print("  - hea_dataset_expanded.csv")
    print("  - experiment_results.json")
    print("  - figures/phase_diagram.png")
    print("  - figures/descriptor_distributions.png")
    print("  - figures/property_correlations.png")
    print("  - figures/model_performance.png")
    print("  - figures/feature_importance.png")
    print("  - figures/bayesian_optimization.png")
    print("  - figures/pareto_front.png")
    print("  - figures/composition_heatmap.png")
    print("  - figures/active_learning.png")
    print("  - figures/case_study_crmnfeconi.png")
    print("=" * 70)
    
    return results

if __name__ == '__main__':
    results = main()
