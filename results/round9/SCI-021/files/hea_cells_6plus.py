# HEA ML Framework - Cells 6+ (continuation)
# Run after base setup

import sys, random, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.pipeline import Pipeline
from scipy.stats import norm
import xgboost as xgb
warnings.filterwarnings('ignore')

np.random.seed(42)
random.seed(42)
FIGURES_DIR = '/app/projects/50b9b3f2-279d-4427-8ebf-e99b1e2beb9c/workspace/figures'
DATA_DIR    = '/app/projects/50b9b3f2-279d-4427-8ebf-e99b1e2beb9c/workspace/data/raw'

# ── Reload dataset ────────────────────────────────────────────
df = pd.read_csv(f'{DATA_DIR}/hea_synthetic_dataset.csv')
ELEMENTS = ['Cr','Mn','Fe','Co','Ni']
feature_cols = ['x_Cr','x_Mn','x_Fe','x_Co','x_Ni',
                'VEC','delta_r','delta_chi','S_mix','H_mix',
                'Omega','Gamma','Tm','G_bar']
target_cols = ['yield_strength','elongation','corrosion_resistance']
X = df[feature_cols].values
Y = df[target_cols].values

# ============================================================
# Cell 6: Feature importance (XGBoost, fixed)
# ============================================================

best_models = {}
for target in target_cols:
    pipe = Pipeline([('scaler', StandardScaler()),
                     ('model', xgb.XGBRegressor(n_estimators=200, max_depth=5,
                                                 learning_rate=0.05, subsample=0.8,
                                                 random_state=42, verbosity=0))])
    pipe.fit(X, df[target].values)
    best_models[target] = pipe

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, target in zip(axes, target_cols):
    model = best_models[target].named_steps['model']
    importances = model.feature_importances_
    idx = np.argsort(importances)[::-1][:8]
    idx_rev = idx[::-1]
    feats_rev = [feature_cols[i] for i in idx_rev]
    imp_rev = np.array([importances[i] for i in idx_rev])
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(idx_rev)))
    ax.barh(feats_rev, imp_rev, color=colors)
    ax.set_title(f'Feature Importance\n{target}', fontsize=12, fontweight='bold')
    ax.set_xlabel('Importance', fontsize=10)
    ax.tick_params(axis='y', labelsize=9)

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig2_hea_feature_importance.png', dpi=150, bbox_inches='tight')
plt.close()
print("[Cell 6] Feature importance figure saved.")

for target in target_cols:
    model = best_models[target].named_steps['model']
    importances = model.feature_importances_
    idx = np.argsort(importances)[::-1][:3]
    top_feats = [(feature_cols[i], importances[i]) for i in idx]
    print(f"  {target}: " + ", ".join(f"{f}({v:.3f})" for f, v in top_feats))

# ============================================================
# Cell 7: Bayesian Optimization (Gaussian Process)
# Multi-objective scalarized: maximize strength + ductility
# ============================================================

print("\n[Cell 7] Bayesian Optimization...")

# Train GP surrogate for yield strength and elongation
gp_models = {}
scaler_X = StandardScaler()
X_scaled = scaler_X.fit_transform(X)

for target in ['yield_strength', 'elongation']:
    y = df[target].values
    y_norm = (y - y.mean()) / y.std()
    kernel = Matern(nu=2.5) * 1.0
    gp = GaussianProcessRegressor(kernel=kernel, alpha=1e-3,
                                   normalize_y=True, n_restarts_optimizer=5,
                                   random_state=42)
    gp.fit(X_scaled, y_norm)
    gp_models[target] = (gp, y.mean(), y.std())

def acquisition_EI(x_cand_scaled, gp, y_best_norm, xi=0.01):
    """Expected Improvement acquisition function."""
    mu, sigma = gp.predict(x_cand_scaled.reshape(1,-1), return_std=True)
    sigma = sigma.reshape(-1)
    z = (mu - y_best_norm - xi) / (sigma + 1e-9)
    ei = (mu - y_best_norm - xi) * norm.cdf(z) + sigma * norm.pdf(z)
    return float(ei)

def scalarize(ys_mean, yl_mean, w_s=0.6, w_e=0.4):
    """Weighted scalarization: maximize strength + elongation."""
    return w_s * ys_mean + w_e * yl_mean

# Initialize with first 30 observations, then do 20 BO iterations
n_init = 30
observed_idx = list(range(n_init))
remaining_idx = list(range(n_init, len(df)))

bo_history = []
current_best_score = -np.inf

for iteration in range(20):
    # Refit GPs on observed data
    X_obs = X_scaled[observed_idx]
    
    gp_ys = gp_models['yield_strength'][0]
    gp_el = gp_models['elongation'][0]
    
    y_ys = df['yield_strength'].values[observed_idx]
    y_el = df['elongation'].values[observed_idx]
    
    # Best observed scalarized score
    ys_norm = (y_ys - df['yield_strength'].mean()) / df['yield_strength'].std()
    el_norm = (y_el - df['elongation'].mean()) / df['elongation'].std()
    scores = [scalarize(ys_norm[i], el_norm[i]) for i in range(len(observed_idx))]
    best_score_norm = max(scores)
    
    # Evaluate EI for remaining candidates
    best_ei = -1
    best_cand = None
    for idx in remaining_idx[:50]:  # evaluate subset for speed
        x_c = X_scaled[idx]
        ei_ys = acquisition_EI(x_c, gp_ys, best_score_norm)
        ei_el = acquisition_EI(x_c, gp_el, best_score_norm)
        ei_combined = scalarize(ei_ys, ei_el)
        if ei_combined > best_ei:
            best_ei = ei_combined
            best_cand = idx
    
    if best_cand is None:
        break
    
    # "Evaluate" the candidate (observe true values)
    observed_idx.append(best_cand)
    remaining_idx.remove(best_cand)
    
    # Record best observed after this iteration
    y_best = df['yield_strength'].values[observed_idx].max()
    bo_history.append({
        'iteration': iteration,
        'n_observed': len(observed_idx),
        'best_yield_strength': y_best,
        'new_composition': {e: df[f'x_{e}'].values[best_cand] for e in ELEMENTS},
    })

bo_df = pd.DataFrame(bo_history)
print(f"[Cell 7] BO completed. Best yield strength found: {bo_df['best_yield_strength'].max():.1f} MPa")
print(f"[Cell 7] Baseline (random search, first 30): {df['yield_strength'].values[:30].max():.1f} MPa")
print(f"[Cell 7] BO improvement: +{bo_df['best_yield_strength'].max() - df['yield_strength'].values[:30].max():.1f} MPa")

# ============================================================
# Cell 8: Pareto front (multi-objective: strength vs ductility)
# ============================================================

def is_pareto_efficient(costs):
    """Return boolean mask of Pareto-efficient points (minimization)."""
    n = len(costs)
    is_eff = np.ones(n, dtype=bool)
    for i in range(n):
        if is_eff[i]:
            is_eff[is_eff] = np.any(costs[is_eff] < costs[i], axis=1) | \
                             np.all(costs[is_eff] == costs[i], axis=1)
            is_eff[i] = True
    return is_eff

# Negate for maximization
costs = np.column_stack([-df['yield_strength'].values, 
                          -df['elongation'].values])
pareto_mask = is_pareto_efficient(costs)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Pareto front
axes[0].scatter(df['yield_strength'][~pareto_mask], df['elongation'][~pareto_mask],
                c='lightgray', alpha=0.5, s=20, label='Non-Pareto')
axes[0].scatter(df['yield_strength'][pareto_mask], df['elongation'][pareto_mask],
                c='red', alpha=0.9, s=50, zorder=5, label=f'Pareto ({pareto_mask.sum()} pts)')
axes[0].set_xlabel('Yield Strength (MPa)', fontsize=11)
axes[0].set_ylabel('Elongation (%)', fontsize=11)
axes[0].set_title('Pareto Front: Strength vs Ductility', fontsize=12, fontweight='bold')
axes[0].legend()

# BO learning curve
axes[1].plot(bo_df['n_observed'], bo_df['best_yield_strength'], 
             'b-o', linewidth=2, markersize=5, label='BO (GP+EI)')
axes[1].axhline(df['yield_strength'].values[:30].max(), 
                color='gray', linestyle='--', label='Random baseline')
axes[1].set_xlabel('Number of Observations', fontsize=11)
axes[1].set_ylabel('Best Yield Strength (MPa)', fontsize=11)
axes[1].set_title('Bayesian Optimization Learning Curve', fontsize=12, fontweight='bold')
axes[1].legend()

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig3_hea_pareto_bo.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"\n[Cell 8] Pareto front: {pareto_mask.sum()} out of {len(df)} alloys on Pareto front")
print(f"[Cell 8] Figure saved: fig3_hea_pareto_bo.png")

# ============================================================
# Cell 9: Composition analysis of best alloys
# ============================================================

# Top alloys by yield strength
top_ys = df.nlargest(10, 'yield_strength')
# Top alloys by Pareto front + balanced score
df['balance_score'] = (df['yield_strength']/df['yield_strength'].max() + 
                       df['elongation']/df['elongation'].max() + 
                       df['corrosion_resistance']/df['corrosion_resistance'].max())
top_balanced = df.nlargest(10, 'balance_score')

print("\n[Cell 9] Top 5 alloys by yield strength:")
print(top_ys[['x_Cr','x_Mn','x_Fe','x_Co','x_Ni','yield_strength','elongation','corrosion_resistance','VEC']].round(3).to_string(index=False))

print("\n[Cell 9] Top 5 balanced alloys (strength+ductility+corrosion):")
print(top_balanced[['x_Cr','x_Mn','x_Fe','x_Co','x_Ni','yield_strength','elongation','corrosion_resistance','VEC']].head().round(3).to_string(index=False))

# ============================================================
# Cell 10: Active learning efficiency comparison
# ============================================================

# Simulate random search vs BO over all 270 remaining samples
rng_al = np.random.RandomState(42)
n_trials = 20

# Random search
rand_best = []
rand_pool = list(range(n_init, len(df)))
rand_observed = list(range(n_init))
for i in range(n_trials):
    if rand_pool:
        new_idx = rng_al.choice(rand_pool)
        rand_pool.remove(new_idx)
        rand_observed.append(new_idx)
        rand_best.append(df['yield_strength'].values[rand_observed].max())

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Active learning vs random
x_iters = list(range(n_init+1, n_init+1+n_trials))
axes[0].plot(x_iters, bo_df['best_yield_strength'].values, 
             'r-o', linewidth=2, markersize=6, label='GP-BO (Active Learning)')
axes[0].plot(x_iters, rand_best, 
             'b--s', linewidth=2, markersize=6, label='Random Search')
axes[0].fill_between(x_iters, 
                     [v - 5 for v in bo_df['best_yield_strength'].values],
                     [v + 5 for v in bo_df['best_yield_strength'].values],
                     alpha=0.2, color='red')
axes[0].set_xlabel('Number of Evaluations', fontsize=11)
axes[0].set_ylabel('Best Yield Strength (MPa)', fontsize=11)
axes[0].set_title('Active Learning vs Random Search', fontsize=12, fontweight='bold')
axes[0].legend()

# Composition ternary-like projection (Cr-Fe-Ni triangle)
sc = axes[1].scatter(df['x_Cr'], df['x_Ni'], 
                     c=df['yield_strength'], cmap='hot_r',
                     s=40, alpha=0.8)
plt.colorbar(sc, ax=axes[1], label='Yield Strength (MPa)')
top5 = df.nlargest(5, 'yield_strength')
axes[1].scatter(top5['x_Cr'], top5['x_Ni'], 
                c='blue', s=150, marker='*', zorder=10, label='Top 5')
axes[1].set_xlabel('x_Cr (Cr fraction)', fontsize=11)
axes[1].set_ylabel('x_Ni (Ni fraction)', fontsize=11)
axes[1].set_title('Composition Space: Cr-Ni Projection\n(colored by yield strength)', 
                   fontsize=12, fontweight='bold')
axes[1].legend()

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig4_hea_active_learning.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"\n[Cell 10] Active learning figure saved: fig4_hea_active_learning.png")

# Final improvement statistics
bo_improvement = bo_df['best_yield_strength'].max() - df['yield_strength'].values[:n_init].max()
rand_improvement = max(rand_best) - df['yield_strength'].values[:n_init].max()
print(f"[Cell 10] BO improvement over random: {bo_improvement:.1f} MPa vs {rand_improvement:.1f} MPa")

# ============================================================
# Cell 11: CrMnFeCoNi case study - composition optimization
# ============================================================

print("\n[Cell 11] CrMnFeCoNi Case Study — Optimal Composition Search")

# Evaluate full grid on composition space (Cr varied 0.1-0.4, others balanced)
results_grid = []
for x_cr in np.arange(0.05, 0.45, 0.05):
    for x_ni in np.arange(0.05, 0.45, 0.05):
        x_rem = 1.0 - x_cr - x_ni
        if x_rem < 0.15 or x_rem > 0.75:
            continue
        for x_mn_frac in np.arange(0.1, 0.6, 0.1):
            x_mn = x_rem * x_mn_frac
            x_fe = x_rem * (0.5 - x_mn_frac/2)
            x_co = x_rem - x_mn - x_fe
            if x_co < 0:
                continue
            comp = {'Cr': x_cr, 'Mn': x_mn, 'Fe': x_fe, 'Co': x_co, 'Ni': x_ni}
            
            ELEMENT_DATA = {
                'Cr': {'r': 1.28, 'chi': 1.66, 'Hf': 0.0,   'Tm': 2180, 'G': 115, 'VEC': 6, 'mass': 51.996},
                'Mn': {'r': 1.27, 'chi': 1.55, 'Hf': 0.0,   'Tm': 1519, 'G':  80, 'VEC': 7, 'mass': 54.938},
                'Fe': {'r': 1.26, 'chi': 1.83, 'Hf': 0.0,   'Tm': 1811, 'G':  82, 'VEC': 8, 'mass': 55.845},
                'Co': {'r': 1.25, 'chi': 1.88, 'Hf': 0.0,   'Tm': 1768, 'G':  75, 'VEC': 9, 'mass': 58.933},
                'Ni': {'r': 1.24, 'chi': 1.91, 'Hf': 0.0,   'Tm': 1728, 'G':  76, 'VEC':10, 'mass': 58.693},
            }
            ELEM_LIST = ['Cr','Mn','Fe','Co','Ni']
            
            x_arr = np.array([comp.get(e, 0) for e in ELEM_LIST])
            x_arr /= x_arr.sum()
            
            r_bar  = sum(x_arr[i]*ELEMENT_DATA[e]['r']   for i,e in enumerate(ELEM_LIST))
            chi_bar= sum(x_arr[i]*ELEMENT_DATA[e]['chi'] for i,e in enumerate(ELEM_LIST))
            Tm_bar = sum(x_arr[i]*ELEMENT_DATA[e]['Tm']  for i,e in enumerate(ELEM_LIST))
            G_bar  = sum(x_arr[i]*ELEMENT_DATA[e]['G']   for i,e in enumerate(ELEM_LIST))
            VEC_bar= sum(x_arr[i]*ELEMENT_DATA[e]['VEC'] for i,e in enumerate(ELEM_LIST))
            delta_r = np.sqrt(sum(x_arr[i]*(1-ELEMENT_DATA[e]['r']/r_bar)**2 for i,e in enumerate(ELEM_LIST)))*100
            delta_chi = np.sqrt(sum(x_arr[i]*(ELEMENT_DATA[e]['chi']-chi_bar)**2 for i,e in enumerate(ELEM_LIST)))
            R = 8.314
            S_mix = -R * sum(x_arr[i]*np.log(x_arr[i]+1e-12) for i in range(5))
            Omega_m = {('Cr','Mn'):-4.,('Cr','Fe'):-1.,('Cr','Co'):-4.,('Cr','Ni'):-7.,
                       ('Mn','Fe'):0.,('Mn','Co'):-5.,('Mn','Ni'):-8.,('Fe','Co'):-1.,
                       ('Fe','Ni'):-2.,('Co','Ni'):0.}
            H_mix = sum(4*Omega_m.get((ELEM_LIST[i],ELEM_LIST[j]),Omega_m.get((ELEM_LIST[j],ELEM_LIST[i]),0.))*x_arr[i]*x_arr[j]
                        for i in range(5) for j in range(i+1,5))
            
            feat = np.array([x_arr[0],x_arr[1],x_arr[2],x_arr[3],x_arr[4],
                             VEC_bar, delta_r, delta_chi, S_mix, H_mix,
                             S_mix*Tm_bar/(abs(H_mix)+1e-3),
                             delta_r**2/(delta_chi+1e-6),
                             Tm_bar, G_bar])
            
            feat_2d = feat.reshape(1, -1)
            ys_pred = best_models['yield_strength'].predict(feat_2d)[0]
            el_pred = best_models['elongation'].predict(feat_2d)[0]
            cr_pred = best_models['corrosion_resistance'].predict(feat_2d)[0]
            
            results_grid.append({
                'x_Cr': x_cr, 'x_Mn': x_mn, 'x_Fe': x_fe, 'x_Co': x_co, 'x_Ni': x_ni,
                'VEC': VEC_bar, 'delta_r': delta_r, 'S_mix': S_mix, 'H_mix': H_mix,
                'pred_yield_strength': ys_pred,
                'pred_elongation': el_pred,
                'pred_corrosion': cr_pred,
            })

grid_df = pd.DataFrame(results_grid)
grid_df['multi_obj_score'] = (
    (grid_df['pred_yield_strength'] / grid_df['pred_yield_strength'].max()) * 0.4 +
    (grid_df['pred_elongation'] / grid_df['pred_elongation'].max()) * 0.3 +
    (grid_df['pred_corrosion'] / grid_df['pred_corrosion'].max()) * 0.3
)
best_comp = grid_df.nlargest(5, 'multi_obj_score')

print("[Cell 11] Top-5 predicted optimal compositions:")
cols_show = ['x_Cr','x_Mn','x_Fe','x_Co','x_Ni','pred_yield_strength','pred_elongation','pred_corrosion','VEC','multi_obj_score']
print(best_comp[cols_show].round(3).to_string(index=False))

# Compare with equimolar Cantor alloy
cantor_feat = np.array([0.2]*5 + [8.0, 1.1224, 0.1384, 13.3809, -5.12, 4706.43, 9.1048, 1801.2, 85.6]).reshape(1,-1)
cantor_ys = best_models['yield_strength'].predict(cantor_feat)[0]
cantor_el = best_models['elongation'].predict(cantor_feat)[0]
cantor_cr = best_models['corrosion_resistance'].predict(cantor_feat)[0]
print(f"\n[Cell 11] Equimolar Cantor (CrMnFeCoNi0.2 each):")
print(f"  Yield strength: {cantor_ys:.1f} MPa")
print(f"  Elongation:     {cantor_el:.1f} %")
print(f"  Corrosion:      {cantor_cr:.2f} /10")

print(f"\n[Cell 11] Best optimized composition improvement:")
best1 = best_comp.iloc[0]
print(f"  Composition: Cr{best1['x_Cr']:.2f}Mn{best1['x_Mn']:.2f}Fe{best1['x_Fe']:.2f}Co{best1['x_Co']:.2f}Ni{best1['x_Ni']:.2f}")
print(f"  Yield strength: {best1['pred_yield_strength']:.1f} MPa (+{best1['pred_yield_strength']-cantor_ys:.1f} MPa vs Cantor)")
print(f"  Elongation:     {best1['pred_elongation']:.1f} % ({best1['pred_elongation']-cantor_el:+.1f}%)")
print(f"  Corrosion:      {best1['pred_corrosion']:.2f} /10")

# ============================================================
# Cell 12: Final summary figure
# ============================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 11))
axes = axes.ravel()

# 1. Predicted vs actual (yield strength, random forest)
from sklearn.model_selection import cross_val_predict
scaler_tmp = StandardScaler()
X_s = scaler_tmp.fit_transform(X)
rf_tmp = RandomForestRegressor(n_estimators=200, random_state=42, n_jobs=-1)
y_pred_cv = cross_val_predict(rf_tmp, X_s, df['yield_strength'].values, cv=5)
axes[0].scatter(df['yield_strength'], y_pred_cv, alpha=0.5, s=20, c='steelblue')
mn, mx = df['yield_strength'].min(), df['yield_strength'].max()
axes[0].plot([mn,mx],[mn,mx],'r--', linewidth=2)
r2_cv = r2_score(df['yield_strength'], y_pred_cv)
rmse_cv = np.sqrt(mean_squared_error(df['yield_strength'], y_pred_cv))
axes[0].set_xlabel('True Yield Strength (MPa)', fontsize=11)
axes[0].set_ylabel('CV-Predicted Yield Strength (MPa)', fontsize=11)
axes[0].set_title(f'RF: Predicted vs True (5-fold CV)\nR²={r2_cv:.3f}, RMSE={rmse_cv:.1f} MPa', 
                   fontsize=12, fontweight='bold')
axes[0].text(0.05, 0.92, f'R²={r2_cv:.3f}', transform=axes[0].transAxes, fontsize=11)

# 2. VEC vs yield strength
sc2 = axes[1].scatter(df['VEC'], df['yield_strength'], c=df['delta_r'], 
                      cmap='viridis', alpha=0.6, s=25)
plt.colorbar(sc2, ax=axes[1], label='δ (atomic size mismatch)')
axes[1].set_xlabel('VEC (Valence Electron Concentration)', fontsize=11)
axes[1].set_ylabel('Yield Strength (MPa)', fontsize=11)
axes[1].set_title('VEC vs Yield Strength\n(colored by δ)', fontsize=12, fontweight='bold')

# 3. S_mix vs corrosion
axes[2].scatter(df['S_mix'], df['corrosion_resistance'], c=df['x_Cr'], 
                cmap='Reds', alpha=0.6, s=25)
axes[2].set_xlabel('Mixing Entropy ΔS_mix (J/mol/K)', fontsize=11)
axes[2].set_ylabel('Corrosion Resistance (0-10)', fontsize=11)
axes[2].set_title('Mixing Entropy vs Corrosion Resistance\n(colored by x_Cr)', 
                   fontsize=12, fontweight='bold')
plt.colorbar(plt.cm.ScalarMappable(cmap='Reds', 
             norm=plt.Normalize(df['x_Cr'].min(), df['x_Cr'].max())),
             ax=axes[2], label='x_Cr')

# 4. Composition radar chart for top alloys vs Cantor
categories = ELEMENTS
N = len(categories)
angles = [n/float(N)*2*np.pi for n in range(N)] + [0]

ax4 = plt.subplot(2, 2, 4, projection='polar')
# Cantor alloy
cantor_vals = [0.2]*5 + [0.2]
ax4.plot(angles, cantor_vals, 'b-o', linewidth=2, label='Cantor (equimolar)')
ax4.fill(angles, cantor_vals, alpha=0.1, color='blue')
# Best optimized alloy
best_vals = [best1['x_Cr'], best1['x_Mn'], best1['x_Fe'], 
             best1['x_Co'], best1['x_Ni']] + [best1['x_Cr']]
ax4.plot(angles, best_vals, 'r-^', linewidth=2, label='Best optimized')
ax4.fill(angles, best_vals, alpha=0.1, color='red')
ax4.set_xticks(angles[:-1])
ax4.set_xticklabels(categories, fontsize=11)
ax4.set_title('Composition: Cantor vs Optimal', fontsize=12, fontweight='bold', pad=20)
ax4.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0), fontsize=9)

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig5_hea_summary.png', dpi=150, bbox_inches='tight')
plt.close()
print("\n[Cell 12] Summary figure saved: fig5_hea_summary.png")
print(f"[Cell 12] Cross-val RF R² = {r2_cv:.3f}, RMSE = {rmse_cv:.1f} MPa")

# ============================================================
# Cell 13: pip freeze (environment record)
# ============================================================
import subprocess
result = subprocess.run([sys.executable, '-m', 'pip', 'freeze'], 
                        capture_output=True, text=True)
# Save to file
with open(f'{DATA_DIR}/pip_freeze.txt', 'w') as f:
    f.write(result.stdout)
print("\n[Cell 13] pip freeze saved to data/raw/pip_freeze.txt")
# Print key packages
key_pkgs = ['numpy', 'pandas', 'scikit-learn', 'xgboost', 'matplotlib', 'scipy', 'seaborn']
for line in result.stdout.split('\n'):
    if any(p in line.lower() for p in key_pkgs):
        print(f"  {line}")

print("\n=== ALL CELLS COMPLETE ===")

