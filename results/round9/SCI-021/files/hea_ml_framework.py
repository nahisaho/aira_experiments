# ============================================================
# HEA Multi-objective ML Optimization Framework
# Seed: 42 | Python 3.11
# ============================================================

import sys, random, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy.stats import norm, spearmanr, pearsonr
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, Matern
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.decomposition import PCA
import xgboost as xgb
warnings.filterwarnings('ignore')

np.random.seed(42)
random.seed(42)
RANDOM_STATE = 42
FIGURES_DIR = '/app/projects/50b9b3f2-279d-4427-8ebf-e99b1e2beb9c/workspace/figures'
DATA_DIR = '/app/projects/50b9b3f2-279d-4427-8ebf-e99b1e2beb9c/workspace/data/raw'

# ============================================================
# Cell 1: HEA Physical Property Database (CrMnFeCoNi system)
# Reference atomic properties for descriptor calculation
# ============================================================

# Atomic properties for Cr, Mn, Fe, Co, Ni
# Sources: Goldschmidt radii, electronegativity (Pauling), 
# melting points, bulk moduli from literature
ELEMENT_DATA = {
    #       r(Å)   χ(Paul)  ΔHf(kJ/mol) Tm(K)  G(GPa) VEC
    'Cr': {'r': 1.28, 'chi': 1.66, 'Hf': 0.0,   'Tm': 2180, 'G': 115, 'VEC': 6, 'mass': 51.996},
    'Mn': {'r': 1.27, 'chi': 1.55, 'Hf': 0.0,   'Tm': 1519, 'G':  80, 'VEC': 7, 'mass': 54.938},
    'Fe': {'r': 1.26, 'chi': 1.83, 'Hf': 0.0,   'Tm': 1811, 'G':  82, 'VEC': 8, 'mass': 55.845},
    'Co': {'r': 1.25, 'chi': 1.88, 'Hf': 0.0,   'Tm': 1768, 'G':  75, 'VEC': 9, 'mass': 58.933},
    'Ni': {'r': 1.24, 'chi': 1.91, 'Hf': 0.0,   'Tm': 1728, 'G':  76, 'VEC':10, 'mass': 58.693},
}
ELEMENTS = ['Cr', 'Mn', 'Fe', 'Co', 'Ni']

print("[Cell 1] Element database loaded:", ELEMENTS)

# ============================================================
# Cell 2: Descriptor calculation functions
# (atomic size mismatch, VEC, mixing entropy, etc.)
# ============================================================

def calc_descriptors(comp_dict):
    """
    Compute HEA descriptors from a composition dictionary.
    comp_dict: {'Cr': 0.2, 'Mn': 0.2, ...} (must sum to 1)
    Returns a dict of descriptors.
    """
    elems = [e for e in ELEMENTS if e in comp_dict]
    x = np.array([comp_dict.get(e, 0.0) for e in ELEMENTS])
    x = x / x.sum()  # normalize
    
    props = ELEMENT_DATA
    
    # Mean values (compositionally weighted)
    r_bar  = sum(x[i]*props[e]['r']   for i,e in enumerate(ELEMENTS))
    chi_bar= sum(x[i]*props[e]['chi'] for i,e in enumerate(ELEMENTS))
    Tm_bar = sum(x[i]*props[e]['Tm']  for i,e in enumerate(ELEMENTS))
    G_bar  = sum(x[i]*props[e]['G']   for i,e in enumerate(ELEMENTS))
    VEC_bar= sum(x[i]*props[e]['VEC'] for i,e in enumerate(ELEMENTS))
    
    # Atomic size mismatch δ (Ye et al. formula)
    delta_r = np.sqrt(sum(x[i]*(1 - props[e]['r']/r_bar)**2 for i,e in enumerate(ELEMENTS))) * 100
    
    # Electronegativity difference Δχ
    delta_chi = np.sqrt(sum(x[i]*(props[e]['chi'] - chi_bar)**2 for i,e in enumerate(ELEMENTS)))
    
    # Mixing entropy ΔS_mix (ideal, J/mol/K)
    R = 8.314
    S_mix = -R * sum(x[i]*np.log(x[i]+1e-12) for i in range(len(ELEMENTS)))
    
    # Mixing enthalpy ΔH_mix (simplified regular solution)
    # Using Miedema-style interaction parameters (approximate)
    Omega_matrix = {
        ('Cr','Mn'): -4.0, ('Cr','Fe'): -1.0, ('Cr','Co'): -4.0, ('Cr','Ni'): -7.0,
        ('Mn','Fe'):  0.0, ('Mn','Co'): -5.0, ('Mn','Ni'): -8.0,
        ('Fe','Co'):  -1.0, ('Fe','Ni'): -2.0,
        ('Co','Ni'):  0.0,
    }
    H_mix = 0.0
    for i, ei in enumerate(ELEMENTS):
        for j, ej in enumerate(ELEMENTS):
            if i < j:
                key = (ei, ej) if (ei, ej) in Omega_matrix else (ej, ei)
                omega = Omega_matrix.get(key, 0.0)
                H_mix += 4 * omega * x[i] * x[j]
    
    # Ω parameter (stability criterion)
    T_melt_ref = 1000  # reference temperature K
    Omega_param = S_mix * Tm_bar / (abs(H_mix) + 1e-3)
    
    # Γ parameter (Yang & Zhang)
    Gamma = delta_r**2 / (delta_chi + 1e-6)
    
    return {
        'x_Cr': x[0], 'x_Mn': x[1], 'x_Fe': x[2], 'x_Co': x[3], 'x_Ni': x[4],
        'VEC': VEC_bar,
        'delta_r': delta_r,
        'delta_chi': delta_chi,
        'S_mix': S_mix,
        'H_mix': H_mix,
        'Omega': Omega_param,
        'Gamma': Gamma,
        'Tm': Tm_bar,
        'G_bar': G_bar,
    }

# Test with equimolar Cantor alloy
cantor = {'Cr': 0.2, 'Mn': 0.2, 'Fe': 0.2, 'Co': 0.2, 'Ni': 0.2}
desc_cantor = calc_descriptors(cantor)
print("[Cell 2] Cantor alloy descriptors:")
for k, v in desc_cantor.items():
    print(f"  {k:12s} = {v:.4f}")


# ============================================================
# Cell 3: Synthetic dataset generation
# Physics-informed mock dataset for CrMnFeCoNi space
# ============================================================

def generate_hea_dataset(n_samples=300, seed=42):
    """
    Generate a physics-informed synthetic dataset of CrMnFeCoNi alloys.
    Properties are modeled with physically motivated relationships plus noise.
    """
    rng = np.random.RandomState(seed)
    
    # Random compositions in CrMnFeCoNi simplex
    # Use Dirichlet distribution for uniform sampling on simplex
    concentrations = rng.dirichlet(np.ones(5), size=n_samples)
    
    rows = []
    for comp in concentrations:
        comp_dict = {e: comp[i] for i, e in enumerate(ELEMENTS)}
        d = calc_descriptors(comp_dict)
        
        # ── Yield Strength (MPa) ─────────────────────────────────
        # Based on solid solution strengthening model:
        # σ_y ~ f(δ, G_bar, VEC)
        # Higher δ and G → higher strength; FCC alloys (VEC~8) have moderate strength
        sigma_y = (
            150
            + 80 * d['delta_r']
            + 0.8 * d['G_bar']
            + 10 * abs(d['H_mix'])
            - 15 * abs(d['VEC'] - 8.0)   # penalty away from VEC=8 (FCC stability)
            + rng.normal(0, 15)
        )
        sigma_y = max(100, sigma_y)
        
        # ── Elongation (%) ───────────────────────────────────────
        # FCC structures (VEC 8-10) have higher ductility
        # Low δ favors ductility (less lattice distortion)
        elong = (
            40
            - 4 * d['delta_r']
            - 0.05 * d['G_bar']
            + 5 * max(0, 3 - abs(d['VEC'] - 8.5))
            + 2 * d['x_Ni']
            - 3 * d['x_Mn']
            + rng.normal(0, 3)
        )
        elong = max(2, min(70, elong))
        
        # ── Corrosion resistance (scale 0-10) ────────────────────
        # Cr content >15% strongly improves corrosion resistance
        # Ni also beneficial; Mn slightly detrimental
        corr = (
            4.0
            + 20 * d['x_Cr']
            + 8 * d['x_Ni']
            - 5 * d['x_Mn']
            - 0.2 * d['delta_chi']
            + rng.normal(0, 0.5)
        )
        corr = max(0, min(10, corr))
        
        # ── Phase stability (1=FCC, 0=mixed/BCC) ────────────────
        # FCC stable when VEC ≥ 8, Ω > 1.5, δ < 6.6
        fcc_prob = 1 / (1 + np.exp(-(d['VEC'] - 8.0)*2 - (d['Omega']/1000 - 1.5)))
        phase = int(rng.random() < fcc_prob)
        
        row = {**d, 'yield_strength': sigma_y, 'elongation': elong, 
               'corrosion_resistance': corr, 'is_fcc': phase}
        rows.append(row)
    
    df = pd.DataFrame(rows)
    return df

df = generate_hea_dataset(300, seed=42)
print(f"[Cell 3] Dataset shape: {df.shape}")
print(df[['yield_strength','elongation','corrosion_resistance','is_fcc','VEC','delta_r','S_mix']].describe().round(3))

# Save raw data
df.to_csv(f'{DATA_DIR}/hea_synthetic_dataset.csv', index=False)
print(f"[Cell 3] Saved to {DATA_DIR}/hea_synthetic_dataset.csv")


# ============================================================
# Cell 4: Feature correlation and PCA analysis
# ============================================================

feature_cols = ['x_Cr','x_Mn','x_Fe','x_Co','x_Ni',
                'VEC','delta_r','delta_chi','S_mix','H_mix',
                'Omega','Gamma','Tm','G_bar']
target_cols = ['yield_strength','elongation','corrosion_resistance']

X = df[feature_cols].values
Y = df[target_cols].values

# Correlation matrix between descriptors and targets
corr_data = df[feature_cols + target_cols]
corr_matrix = corr_data.corr()

fig, axes = plt.subplots(1, 2, figsize=(18, 7))

# Heatmap of descriptor correlations with targets
target_corr = corr_matrix.loc[feature_cols, target_cols]
sns.heatmap(target_corr, annot=True, fmt='.2f', cmap='RdBu_r', 
            center=0, ax=axes[0], vmin=-1, vmax=1,
            linewidths=0.5)
axes[0].set_title('Descriptor–Property Correlation Matrix\n(CrMnFeCoNi System)', 
                   fontsize=13, fontweight='bold')
axes[0].set_ylabel('Descriptors', fontsize=11)
axes[0].set_xlabel('Target Properties', fontsize=11)

# PCA biplot
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
pca = PCA(n_components=2, random_state=42)
X_pca = pca.fit_transform(X_scaled)
explained = pca.explained_variance_ratio_

scatter = axes[1].scatter(X_pca[:, 0], X_pca[:, 1],
                          c=df['yield_strength'], cmap='plasma',
                          s=30, alpha=0.7)
plt.colorbar(scatter, ax=axes[1], label='Yield Strength (MPa)')
axes[1].set_xlabel(f'PC1 ({explained[0]*100:.1f}% var.)', fontsize=11)
axes[1].set_ylabel(f'PC2 ({explained[1]*100:.1f}% var.)', fontsize=11)
axes[1].set_title('PCA Projection (colored by Yield Strength)', 
                   fontsize=13, fontweight='bold')

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig1_hea_descriptors_pca.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"[Cell 4] PCA variance explained: PC1={explained[0]*100:.1f}%, PC2={explained[1]*100:.1f}%")
print(f"[Cell 4] Total (2 PCs): {sum(explained)*100:.1f}%")
print(f"[Cell 4] Figure saved: fig1_hea_descriptors_pca.png")

# Top correlations with yield strength
ys_corr = [(f, abs(corr_matrix.loc[f, 'yield_strength'])) for f in feature_cols]
ys_corr_sorted = sorted(ys_corr, key=lambda x: -x[1])
print("[Cell 4] Top descriptor correlations with yield strength:")
for f, c in ys_corr_sorted[:5]:
    print(f"  {f:12s}: |r| = {c:.3f}")


# ============================================================
# Cell 5: Multi-target ML model training (RF + XGBoost)
# with 5-fold cross-validation
# ============================================================

from sklearn.model_selection import KFold, cross_validate
from sklearn.pipeline import Pipeline

kf = KFold(n_splits=5, shuffle=True, random_state=42)

# Models to evaluate
models = {
    'RandomForest': RandomForestRegressor(n_estimators=200, max_depth=8, 
                                          random_state=42, n_jobs=-1),
    'XGBoost': xgb.XGBRegressor(n_estimators=200, max_depth=5, learning_rate=0.05,
                                  subsample=0.8, random_state=42, verbosity=0),
}

results_cv = {}
for target in target_cols:
    y = df[target].values
    results_cv[target] = {}
    for model_name, model in models.items():
        pipe = Pipeline([('scaler', StandardScaler()), ('model', model)])
        cv_res = cross_validate(pipe, X, y, cv=kf, 
                                scoring=['r2', 'neg_root_mean_squared_error'],
                                return_train_score=True)
        results_cv[target][model_name] = {
            'R2_val_mean':  cv_res['test_r2'].mean(),
            'R2_val_std':   cv_res['test_r2'].std(),
            'RMSE_val_mean': -cv_res['test_neg_root_mean_squared_error'].mean(),
            'RMSE_val_std':  cv_res['test_neg_root_mean_squared_error'].std(),
            'R2_train_mean': cv_res['train_r2'].mean(),
        }

print("[Cell 5] 5-fold CV Results (mean ± std):")
print(f"{'Target':22s} {'Model':14s} {'R² val':>12s} {'RMSE val':>14s} {'R² train':>12s}")
print("-" * 80)
for target in target_cols:
    for model_name, res in results_cv[target].items():
        print(f"{target:22s} {model_name:14s} "
              f"{res['R2_val_mean']:>8.3f}±{res['R2_val_std']:.3f}  "
              f"{res['RMSE_val_mean']:>8.3f}±{res['RMSE_val_std']:.3f}  "
              f"{res['R2_train_mean']:>10.3f}")


# ============================================================
# Cell 6: Feature importance analysis (XGBoost)
# ============================================================

from sklearn.pipeline import Pipeline

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
    ax.barh([feature_cols[i] for i in reversed(idx)],
            importances[reversed(idx)],
            color=plt.cm.viridis(np.linspace(0.2, 0.9, 8)))
    ax.set_title(f'Feature Importance\n{target}', fontsize=12, fontweight='bold')
    ax.set_xlabel('Importance', fontsize=10)
    ax.tick_params(axis='y', labelsize=9)

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig2_hea_feature_importance.png', dpi=150, bbox_inches='tight')
plt.close()
print("[Cell 6] Feature importance figure saved: fig2_hea_feature_importance.png")

# Top features per target
print("[Cell 6] Top-3 features per target:")
for target in target_cols:
    model = best_models[target].named_steps['model']
    importances = model.feature_importances_
    idx = np.argsort(importances)[::-1][:3]
    top_feats = [(feature_cols[i], importances[i]) for i in idx]
    print(f"  {target}: " + ", ".join(f"{f}({v:.3f})" for f, v in top_feats))


# ============================================================
# Cell 6 (fixed): Feature importance analysis (XGBoost)
# ============================================================

from sklearn.pipeline import Pipeline

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
    # Reverse for horizontal bar
    idx_rev = idx[::-1]
    feats_rev = [feature_cols[i] for i in idx_rev]
    imp_rev = importances[idx_rev]
    ax.barh(feats_rev, imp_rev, color=plt.cm.viridis(np.linspace(0.2, 0.9, len(idx_rev))))
    ax.set_title(f'Feature Importance\n{target}', fontsize=12, fontweight='bold')
    ax.set_xlabel('Importance', fontsize=10)
    ax.tick_params(axis='y', labelsize=9)

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fig2_hea_feature_importance.png', dpi=150, bbox_inches='tight')
plt.close()
print("[Cell 6] Feature importance figure saved: fig2_hea_feature_importance.png")

# Top features per target
print("[Cell 6] Top-3 features per target:")
for target in target_cols:
    model = best_models[target].named_steps['model']
    importances = model.feature_importances_
    idx = np.argsort(importances)[::-1][:3]
    top_feats = [(feature_cols[i], importances[i]) for i in idx]
    print(f"  {target}: " + ", ".join(f"{f}({v:.3f})" for f, v in top_feats))

