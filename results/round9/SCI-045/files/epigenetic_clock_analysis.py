"""
Epigenetic Clock Analysis: Improved Biological Age Estimation
Cell-by-cell reference: [cell:N] format for paper citations

Environment:
- Python 3.x
- numpy, pandas, scikit-learn, xgboost, lightgbm, matplotlib, seaborn, scipy
- Random seed: 42 (all models)
"""

# ============================================================
# [cell:0] Environment setup and seed fixing
# ============================================================
import numpy as np
import pandas as pd
import random
import sys
import os

# Fix all random seeds for reproducibility
np.random.seed(42)
random.seed(42)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')

from sklearn.linear_model import ElasticNet, Lasso, Ridge
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, cross_val_score, train_test_split
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error
from sklearn.pipeline import Pipeline
import xgboost as xgb
import lightgbm as lgb

# Output directories
os.makedirs('figures', exist_ok=True)
os.makedirs('data/raw', exist_ok=True)

print("=" * 60)
print("EPIGENETIC CLOCK ANALYSIS")
print(f"Python: {sys.version.split()[0]}")
print(f"NumPy: {np.__version__}")
print(f"Pandas: {pd.__version__}")
print(f"sklearn: {__import__('sklearn').__version__}")
print(f"XGBoost: {xgb.__version__}")
print(f"LightGBM: {lgb.__version__}")
print("Random seed: 42")
print("=" * 60)


# ============================================================
# [cell:1] Simulate realistic DNA methylation dataset
# ============================================================
"""
Simulated dataset design:
- N=800 samples, 1000 CpG probes
- 5 tissue types: blood, brain, liver, lung, breast
- Age range: 18-90 years
- Age-correlated CpGs: 200 (100 hypermethylated, 100 hypomethylated)
- Noise: Gaussian sigma=0.05
- Intervention groups: control, exercise, diet, drug
- Data provenance: synthetically generated, saved to data/raw/
"""

np.random.seed(42)

N_SAMPLES = 800
N_CPGS = 1000
N_INFORMATIVE = 200  # CpGs correlated with age
TISSUES = ['blood', 'brain', 'liver', 'lung', 'breast']
TISSUE_EFFECT = {'blood': 0.0, 'brain': 2.5, 'liver': 1.2, 'lung': -1.0, 'breast': -2.0}

# Generate chronological ages
ages = np.random.uniform(18, 90, N_SAMPLES)

# Tissue assignment
tissue_ids = np.random.choice(TISSUES, N_SAMPLES)

# Generate methylation data
X = np.random.uniform(0.1, 0.9, (N_SAMPLES, N_CPGS))

# Add age-correlated signal
hypermeth_idx = np.arange(0, 100)   # CpGs that gain methylation with age
hypometh_idx = np.arange(100, 200)  # CpGs that lose methylation with age

for i, age in enumerate(ages):
    age_norm = (age - 18) / 72  # normalize 0-1
    # Hyper-methylation with age
    X[i, hypermeth_idx] += age_norm * 0.35 + np.random.normal(0, 0.05, 100)
    X[i, hypermeth_idx] = np.clip(X[i, hypermeth_idx], 0, 1)
    # Hypo-methylation with age
    X[i, hypometh_idx] -= age_norm * 0.25 + np.random.normal(0, 0.05, 100)
    X[i, hypometh_idx] = np.clip(X[i, hypometh_idx], 0, 1)

# Add tissue-specific effects
for i, tissue in enumerate(tissue_ids):
    tissue_noise = np.random.normal(TISSUE_EFFECT[tissue], 0.3, N_INFORMATIVE)
    X[i, :N_INFORMATIVE] += tissue_noise * 0.02
    X[i, :N_INFORMATIVE] = np.clip(X[i, :N_INFORMATIVE], 0, 1)

# Generate biological age (true) = chronological age + lifestyle effects + noise
lifestyle_effect = np.random.normal(0, 3, N_SAMPLES)
biological_age = ages + lifestyle_effect + np.random.normal(0, 2, N_SAMPLES)
biological_age = np.clip(biological_age, 15, 95)

# Age acceleration = biological - chronological
age_acceleration = biological_age - ages

# Intervention assignment
interventions = np.random.choice(['control', 'exercise', 'diet', 'drug'], N_SAMPLES, 
                                   p=[0.4, 0.2, 0.2, 0.2])

# Apply intervention effects (simulating biological age reversal)
intervention_effect = {
    'control': 0,
    'exercise': -1.8,  # 1.8 years younger
    'diet': -1.5,
    'drug': -2.5
}
for i, intervention in enumerate(interventions):
    biological_age[i] += intervention_effect[intervention] + np.random.normal(0, 1.0)
biological_age = np.clip(biological_age, 15, 95)

# Create DataFrame
cpg_names = [f'cg{i:07d}' for i in range(N_CPGS)]
df = pd.DataFrame(X, columns=cpg_names)
df['chronological_age'] = ages
df['biological_age'] = biological_age
df['age_acceleration'] = biological_age - ages
df['tissue'] = tissue_ids
df['intervention'] = interventions

# Save raw data
df.to_csv('data/raw/methylation_data.csv', index=False)

print(f"\n[cell:1] Dataset Summary:")
print(f"  N samples: {N_SAMPLES}")
print(f"  N CpG probes: {N_CPGS}")
print(f"  N informative CpGs: {N_INFORMATIVE}")
print(f"  Age range: {ages.min():.1f} - {ages.max():.1f} years")
print(f"  Mean age: {ages.mean():.1f} ± {ages.std():.1f}")
print(f"  Mean age acceleration: {age_acceleration.mean():.2f} ± {age_acceleration.std():.2f} years")
print(f"  Tissue distribution: {pd.Series(tissue_ids).value_counts().to_dict()}")
print(f"  Intervention distribution: {pd.Series(interventions).value_counts().to_dict()}")


# ============================================================
# [cell:2] Baseline Horvath-like clock (ElasticNet)
# ============================================================
print("\n" + "=" * 60)
print("[cell:2] Baseline Models (Horvath-like ElasticNet)")
print("=" * 60)

X_features = df[cpg_names].values
y_chrono = df['chronological_age'].values
y_bio = df['biological_age'].values

# Standard 5-fold cross-validation
kf = KFold(n_splits=5, shuffle=True, random_state=42)

# Model definitions
models = {
    'ElasticNet (Horvath-like)': Pipeline([
        ('scaler', StandardScaler()),
        ('model', ElasticNet(alpha=0.1, l1_ratio=0.5, random_state=42, max_iter=5000))
    ]),
    'Lasso': Pipeline([
        ('scaler', StandardScaler()),
        ('model', Lasso(alpha=0.05, random_state=42, max_iter=5000))
    ]),
    'Ridge': Pipeline([
        ('scaler', StandardScaler()),
        ('model', Ridge(alpha=1.0))
    ]),
    'Random Forest': Pipeline([
        ('model', RandomForestRegressor(n_estimators=100, max_features=0.3, 
                                        min_samples_leaf=5, random_state=42, n_jobs=-1))
    ]),
    'XGBoost': Pipeline([
        ('model', xgb.XGBRegressor(n_estimators=100, learning_rate=0.05, 
                                    max_depth=4, subsample=0.8, 
                                    colsample_bytree=0.3, random_state=42,
                                    verbosity=0))
    ]),
    'LightGBM': Pipeline([
        ('model', lgb.LGBMRegressor(n_estimators=100, learning_rate=0.05,
                                     max_depth=4, subsample=0.8,
                                     colsample_bytree=0.3, random_state=42,
                                     verbose=-1))
    ]),
    'Neural Network (MLP)': Pipeline([
        ('scaler', StandardScaler()),
        ('model', MLPRegressor(hidden_layer_sizes=(256, 128, 64), 
                                activation='relu', 
                                alpha=0.001,
                                batch_size=64,
                                learning_rate='adaptive',
                                max_iter=200, 
                                random_state=42))
    ])
}

results = {}

for model_name, model in models.items():
    # Predict chronological age
    mae_scores = []
    r2_scores = []
    
    y_pred_all = np.zeros_like(y_chrono, dtype=float)
    
    for fold, (train_idx, test_idx) in enumerate(kf.split(X_features)):
        X_train, X_test = X_features[train_idx], X_features[test_idx]
        y_train, y_test = y_chrono[train_idx], y_chrono[test_idx]
        
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        y_pred_all[test_idx] = y_pred
        
        mae_scores.append(mean_absolute_error(y_test, y_pred))
        r2_scores.append(r2_score(y_test, y_pred))
    
    corr, pval = pearsonr(y_chrono, y_pred_all)
    
    results[model_name] = {
        'mae_mean': np.mean(mae_scores),
        'mae_std': np.std(mae_scores),
        'r2_mean': np.mean(r2_scores),
        'r2_std': np.std(r2_scores),
        'pearson_r': corr,
        'pearson_p': pval,
        'y_pred': y_pred_all.copy()
    }
    
    print(f"  {model_name:35s}: MAE={np.mean(mae_scores):.2f}±{np.std(mae_scores):.2f}  "
          f"R²={np.mean(r2_scores):.3f}±{np.std(r2_scores):.3f}  "
          f"r={corr:.3f} (p={pval:.2e})")


# ============================================================
# [cell:3] Tissue-specific clock analysis
# ============================================================
print("\n" + "=" * 60)
print("[cell:3] Tissue-Specific Clock Analysis")
print("=" * 60)

tissue_results = {}
best_model_name = 'ElasticNet (Horvath-like)'

for tissue in TISSUES:
    mask = df['tissue'] == tissue
    X_t = X_features[mask]
    y_t = y_chrono[mask]
    
    if len(y_t) < 30:
        print(f"  {tissue}: N={len(y_t)} (insufficient samples)")
        continue
    
    # 5-fold CV per tissue
    model_t = Pipeline([
        ('scaler', StandardScaler()),
        ('model', ElasticNet(alpha=0.1, l1_ratio=0.5, random_state=42, max_iter=5000))
    ])
    
    kf_t = KFold(n_splits=min(5, len(y_t)//10), shuffle=True, random_state=42)
    mae_t = []
    r2_t = []
    y_pred_t = np.zeros_like(y_t, dtype=float)
    
    for train_i, test_i in kf_t.split(X_t):
        model_t.fit(X_t[train_i], y_t[train_i])
        pred = model_t.predict(X_t[test_i])
        y_pred_t[test_i] = pred
        mae_t.append(mean_absolute_error(y_t[test_i], pred))
        r2_t.append(r2_score(y_t[test_i], pred))
    
    corr_t, _ = pearsonr(y_t, y_pred_t)
    tissue_results[tissue] = {
        'n': len(y_t),
        'mae': np.mean(mae_t),
        'mae_std': np.std(mae_t),
        'r2': np.mean(r2_t),
        'r2_std': np.std(r2_t),
        'pearson_r': corr_t
    }
    print(f"  {tissue:10s} (N={len(y_t):3d}): MAE={np.mean(mae_t):.2f}±{np.std(mae_t):.2f}  "
          f"R²={np.mean(r2_t):.3f}±{np.std(r2_t):.3f}  r={corr_t:.3f}")

# Pan-tissue model
print(f"\n  Pan-tissue model performance:")
pan_mae = []
pan_r2 = []
y_pred_pan = np.zeros_like(y_chrono, dtype=float)

# Include tissue as a categorical encoding
tissue_encoded = pd.get_dummies(df['tissue'], prefix='tissue').values
X_pan = np.hstack([X_features, tissue_encoded])

model_pan = Pipeline([
    ('scaler', StandardScaler()),
    ('model', ElasticNet(alpha=0.1, l1_ratio=0.5, random_state=42, max_iter=5000))
])

for train_i, test_i in kf.split(X_pan):
    model_pan.fit(X_pan[train_i], y_chrono[train_i])
    pred = model_pan.predict(X_pan[test_i])
    y_pred_pan[test_i] = pred
    pan_mae.append(mean_absolute_error(y_chrono[test_i], pred))
    pan_r2.append(r2_score(y_chrono[test_i], pred))

corr_pan, p_pan = pearsonr(y_chrono, y_pred_pan)
print(f"  Pan-tissue with encoding: MAE={np.mean(pan_mae):.2f}±{np.std(pan_mae):.2f}  "
      f"R²={np.mean(pan_r2):.3f}±{np.std(pan_r2):.3f}  r={corr_pan:.3f}")


# ============================================================
# [cell:4] Age Acceleration Analysis
# ============================================================
print("\n" + "=" * 60)
print("[cell:4] Age Acceleration Biomarker Analysis")
print("=" * 60)

# Use ElasticNet predictions for age acceleration calculation
y_pred_best = results['ElasticNet (Horvath-like)']['y_pred']
computed_age_acc = y_bio - y_pred_best  # biological - predicted chronological

# Split by intervention
intervention_accel = {}
for intervention in ['control', 'exercise', 'diet', 'drug']:
    mask = df['intervention'] == intervention
    acc_vals = computed_age_acc[mask]
    intervention_accel[intervention] = {
        'mean': acc_vals.mean(),
        'std': acc_vals.std(),
        'n': mask.sum()
    }
    print(f"  {intervention:10s}: Age Accel = {acc_vals.mean():.2f} ± {acc_vals.std():.2f} "
          f"(N={mask.sum()})")

# ANOVA test for intervention effects
from scipy.stats import f_oneway
groups = [computed_age_acc[df['intervention'] == intv] for intv in ['control', 'exercise', 'diet', 'drug']]
f_stat, p_anova = f_oneway(*groups)
print(f"\n  ANOVA: F={f_stat:.3f}, p={p_anova:.4f}")

# Pairwise t-tests: control vs each intervention
from scipy.stats import ttest_ind
ctrl_acc = computed_age_acc[df['intervention'] == 'control']
print("\n  Pairwise t-tests (vs control):")
ttest_results = {}
for intv in ['exercise', 'diet', 'drug']:
    intv_acc = computed_age_acc[df['intervention'] == intv]
    t_stat, t_p = ttest_ind(ctrl_acc, intv_acc)
    effect_size = (ctrl_acc.mean() - intv_acc.mean()) / np.sqrt((ctrl_acc.std()**2 + intv_acc.std()**2)/2)
    ttest_results[intv] = {'t': t_stat, 'p': t_p, 'cohen_d': effect_size}
    print(f"    control vs {intv:10s}: t={t_stat:.3f}, p={t_p:.4f}, Cohen's d={effect_size:.3f}")


# ============================================================
# [cell:5] Longevity cohort validation strategy
# ============================================================
print("\n" + "=" * 60)
print("[cell:5] Longevity Cohort Validation")
print("=" * 60)

# Simulate centenarian cohort (ages 90-105)
np.random.seed(42)
N_LONG = 150
ages_long = np.random.uniform(90, 105, N_LONG)
X_long = np.random.uniform(0.1, 0.9, (N_LONG, N_CPGS))

for i, age in enumerate(ages_long):
    age_norm = (age - 18) / 87  # extend to centenarian range
    X_long[i, hypermeth_idx] += age_norm * 0.35 + np.random.normal(0, 0.04, 100)
    X_long[i, hypermeth_idx] = np.clip(X_long[i, hypermeth_idx], 0, 1)
    X_long[i, hypometh_idx] -= age_norm * 0.20 + np.random.normal(0, 0.04, 100)
    X_long[i, hypometh_idx] = np.clip(X_long[i, hypometh_idx], 0, 1)

# Biological age of long-lived individuals: systematically younger than predicted
bio_age_long = ages_long - np.random.exponential(4, N_LONG)  # longevity advantage ~4 yrs
bio_age_long = np.clip(bio_age_long, 80, 110)

# Train on full training cohort, test on longevity cohort
model_final = Pipeline([
    ('scaler', StandardScaler()),
    ('model', ElasticNet(alpha=0.1, l1_ratio=0.5, random_state=42, max_iter=5000))
])
model_final.fit(X_features, y_chrono)
pred_long = model_final.predict(X_long)

accel_long = bio_age_long - pred_long
mae_long = mean_absolute_error(ages_long, pred_long)
r2_long = r2_score(ages_long, pred_long)
corr_long, p_long = pearsonr(ages_long, pred_long)

print(f"  Longevity cohort (N={N_LONG}): MAE={mae_long:.2f}, R²={r2_long:.3f}, r={corr_long:.3f}")
print(f"  Mean age acceleration (longevity): {accel_long.mean():.2f} ± {accel_long.std():.2f} years")
print(f"  Expected: longevity cohort should show NEGATIVE age acceleration")

# Compare with training cohort
accel_train = y_bio - y_pred_best
print(f"  Mean age acceleration (training cohort): {accel_train.mean():.2f} ± {accel_train.std():.2f}")

# Statistical comparison
t_long, p_long_t = ttest_ind(accel_long, accel_train)
print(f"  t-test longevity vs training: t={t_long:.3f}, p={p_long_t:.4f}")


# ============================================================
# [cell:6] Neural network deep clock performance
# ============================================================
print("\n" + "=" * 60)
print("[cell:6] Neural Network (Deep Clock) Performance")
print("=" * 60)

nn_results = {}
X_train_nn, X_test_nn, y_train_nn, y_test_nn = train_test_split(
    X_features, y_chrono, test_size=0.2, random_state=42)

# Multi-layer architectures
nn_archs = {
    'MLP-Small': (64, 32),
    'MLP-Medium': (256, 128),
    'MLP-Large': (512, 256, 128),
    'MLP-Deep': (256, 128, 64, 32)
}

for arch_name, hidden_layers in nn_archs.items():
    nn_mae = []
    nn_r2 = []
    y_pred_nn = np.zeros_like(y_chrono, dtype=float)
    
    for train_i, test_i in kf.split(X_features):
        nn = Pipeline([
            ('scaler', StandardScaler()),
            ('model', MLPRegressor(hidden_layer_sizes=hidden_layers,
                                    activation='relu', alpha=0.001,
                                    batch_size=64, learning_rate='adaptive',
                                    max_iter=300, random_state=42))
        ])
        nn.fit(X_features[train_i], y_chrono[train_i])
        pred = nn.predict(X_features[test_i])
        y_pred_nn[test_i] = pred
        nn_mae.append(mean_absolute_error(y_chrono[test_i], pred))
        nn_r2.append(r2_score(y_chrono[test_i], pred))
    
    corr_nn, _ = pearsonr(y_chrono, y_pred_nn)
    nn_results[arch_name] = {
        'mae_mean': np.mean(nn_mae),
        'mae_std': np.std(nn_mae),
        'r2_mean': np.mean(nn_r2),
        'r2_std': np.std(nn_r2),
        'pearson_r': corr_nn,
        'y_pred': y_pred_nn.copy()
    }
    print(f"  {arch_name:20s}: MAE={np.mean(nn_mae):.2f}±{np.std(nn_mae):.2f}  "
          f"R²={np.mean(nn_r2):.3f}±{np.std(nn_r2):.3f}  r={corr_nn:.3f}")


# ============================================================
# [cell:7] Feature importance analysis (CpG selection)
# ============================================================
print("\n" + "=" * 60)
print("[cell:7] CpG Feature Importance Analysis")
print("=" * 60)

# Train Random Forest for feature importance
rf = RandomForestRegressor(n_estimators=200, max_features=0.3, random_state=42, n_jobs=-1)
rf.fit(X_features, y_chrono)
importances = rf.feature_importances_
top_k = 20
top_indices = np.argsort(importances)[-top_k:][::-1]
top_cpgs = [cpg_names[i] for i in top_indices]
top_importances = importances[top_indices]

print(f"  Top {top_k} CpGs by Random Forest importance:")
for i, (cpg, imp) in enumerate(zip(top_cpgs[:10], top_importances[:10])):
    in_informative = int(cpg.replace('cg', '')) < N_INFORMATIVE
    print(f"    {i+1:2d}. {cpg}: {imp:.4f} {'[informative]' if in_informative else ''}")

informative_in_top20 = sum(1 for cpg in top_cpgs if int(cpg.replace('cg', '')) < N_INFORMATIVE)
print(f"\n  Informative CpGs in top 20: {informative_in_top20}/20 ({informative_in_top20/20*100:.0f}%)")


# ============================================================
# [cell:8] FIGURE 1: Age prediction scatter plots
# ============================================================
print("\n[cell:8] Generating Figure 1: Model Comparison Scatter Plots")

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
axes = axes.flatten()

selected_models = ['ElasticNet (Horvath-like)', 'Lasso', 'Ridge', 'Random Forest', 
                   'XGBoost', 'LightGBM', 'Neural Network (MLP)']

for idx, model_name in enumerate(selected_models):
    ax = axes[idx]
    y_pred = results[model_name]['y_pred']
    mae = results[model_name]['mae_mean']
    r2 = results[model_name]['r2_mean']
    r = results[model_name]['pearson_r']
    
    ax.scatter(y_chrono, y_pred, alpha=0.4, s=15, c='steelblue', edgecolors='none')
    ax.plot([18, 90], [18, 90], 'r--', linewidth=2, label='Perfect prediction')
    
    # Add regression line
    z = np.polyfit(y_chrono, y_pred, 1)
    p = np.poly1d(z)
    x_line = np.linspace(18, 90, 100)
    ax.plot(x_line, p(x_line), 'k-', alpha=0.7, linewidth=1.5)
    
    ax.set_xlabel('Chronological Age (years)', fontsize=10)
    ax.set_ylabel('Predicted Age (years)', fontsize=10)
    ax.set_title(f'{model_name}\nMAE={mae:.2f}±{results[model_name]["mae_std"]:.2f}  '
                 f'R²={r2:.3f}  r={r:.3f}', fontsize=9)
    ax.legend(fontsize=8)
    ax.set_xlim(15, 95)
    ax.set_ylim(15, 95)

# Remove unused subplot
axes[-1].remove()

plt.suptitle('DNA Methylation-Based Age Prediction: Model Comparison\n(5-fold CV, N=800)', 
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/fig01_age_prediction_scatter.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/fig01_age_prediction_scatter.png")


# ============================================================
# [cell:9] FIGURE 2: Model performance comparison bar chart
# ============================================================
print("[cell:9] Generating Figure 2: Model Performance Bar Chart")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

model_names_short = [
    'ElasticNet', 'Lasso', 'Ridge', 'RF', 'XGBoost', 'LightGBM', 'MLP'
]
full_names = list(results.keys())
mae_means = [results[m]['mae_mean'] for m in full_names]
mae_stds = [results[m]['mae_std'] for m in full_names]
r2_means = [results[m]['r2_mean'] for m in full_names]
r2_stds = [results[m]['r2_std'] for m in full_names]

colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(full_names)))

# MAE comparison
bars = axes[0].bar(model_names_short, mae_means, yerr=mae_stds, 
                    capsize=5, color=colors, edgecolor='black', linewidth=0.5)
axes[0].set_xlabel('Model', fontsize=12)
axes[0].set_ylabel('Mean Absolute Error (years)', fontsize=12)
axes[0].set_title('Age Prediction MAE\n(5-fold CV ± SD)', fontsize=12)
axes[0].tick_params(axis='x', rotation=30)
for bar, mae in zip(bars, mae_means):
    axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
                  f'{mae:.2f}', ha='center', va='bottom', fontsize=8)

# R² comparison
bars2 = axes[1].bar(model_names_short, r2_means, yerr=r2_stds,
                     capsize=5, color=colors, edgecolor='black', linewidth=0.5)
axes[1].set_xlabel('Model', fontsize=12)
axes[1].set_ylabel('R² Score', fontsize=12)
axes[1].set_title('Age Prediction R²\n(5-fold CV ± SD)', fontsize=12)
axes[1].tick_params(axis='x', rotation=30)
axes[1].set_ylim(0, 1.1)
for bar, r2 in zip(bars2, r2_means):
    axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                  f'{r2:.3f}', ha='center', va='bottom', fontsize=8)

plt.suptitle('Epigenetic Clock Model Comparison', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/fig02_model_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/fig02_model_comparison.png")


# ============================================================
# [cell:10] FIGURE 3: Tissue-specific analysis
# ============================================================
print("[cell:10] Generating Figure 3: Tissue-Specific Analysis")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Tissue MAE comparison
tissue_names = list(tissue_results.keys())
tissue_maes = [tissue_results[t]['mae'] for t in tissue_names]
tissue_stds = [tissue_results[t]['mae_std'] for t in tissue_names]
tissue_r2s = [tissue_results[t]['r2'] for t in tissue_names]

colors_tissue = ['#2196F3', '#4CAF50', '#F44336', '#FF9800', '#9C27B0']
bars = axes[0].bar(tissue_names, tissue_maes, yerr=tissue_stds, capsize=5,
                    color=colors_tissue, edgecolor='black')
axes[0].axhline(y=np.mean(mae_means[:3]), color='red', linestyle='--', 
                 label=f'Pan-tissue ElasticNet ({np.mean(mae_means[:3]):.2f})')
axes[0].set_xlabel('Tissue Type', fontsize=12)
axes[0].set_ylabel('MAE (years)', fontsize=12)
axes[0].set_title('Tissue-Specific Clock Performance\n(MAE ± SD)', fontsize=12)
axes[0].legend()
for bar, mae in zip(bars, tissue_maes):
    axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                  f'{mae:.2f}', ha='center', va='bottom', fontsize=9)

# Age acceleration by tissue
tissue_accel = []
for tissue in TISSUES:
    mask = df['tissue'] == tissue
    accel_vals = computed_age_acc[mask]
    tissue_accel.append(accel_vals)

axes[1].boxplot(tissue_accel, labels=TISSUES, patch_artist=True,
                 medianprops={'linewidth': 2, 'color': 'black'})
for patch, color in zip(axes[1].patches, colors_tissue):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

axes[1].axhline(y=0, color='red', linestyle='--', label='No acceleration')
axes[1].set_xlabel('Tissue Type', fontsize=12)
axes[1].set_ylabel('Age Acceleration (years)', fontsize=12)
axes[1].set_title('Age Acceleration Distribution by Tissue', fontsize=12)
axes[1].legend()

plt.suptitle('Tissue-Specific Epigenetic Clock Analysis', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/fig03_tissue_specific.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/fig03_tissue_specific.png")


# ============================================================
# [cell:11] FIGURE 4: Age Acceleration by Intervention
# ============================================================
print("[cell:11] Generating Figure 4: Intervention Effects on Age Acceleration")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Box plots by intervention
intv_groups = [computed_age_acc[df['intervention'] == intv] 
               for intv in ['control', 'exercise', 'diet', 'drug']]
intv_labels = ['Control', 'Exercise', 'Diet', 'Drug']
colors_intv = ['#9E9E9E', '#4CAF50', '#2196F3', '#F44336']

bp = axes[0].boxplot(intv_groups, labels=intv_labels, patch_artist=True,
                      medianprops={'linewidth': 2, 'color': 'black'})
for patch, color in zip(bp['boxes'], colors_intv):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

axes[0].axhline(y=0, color='black', linestyle='-', linewidth=1)
axes[0].set_xlabel('Intervention Group', fontsize=12)
axes[0].set_ylabel('Epigenetic Age Acceleration (years)', fontsize=12)
axes[0].set_title(f'Intervention Effects on Age Acceleration\n(ANOVA: F={f_stat:.3f}, p={p_anova:.4f})',
                   fontsize=11)

# Add significance annotations
y_max = max([g.max() for g in intv_groups])
for i, (intv, intv_acc) in enumerate(zip(['exercise', 'diet', 'drug'], 
                                           [intv_groups[1], intv_groups[2], intv_groups[3]])):
    t, p = ttest_results[intv]['t'], ttest_results[intv]['p']
    sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'ns'))
    axes[0].text(i+2, y_max + 1.5, sig, ha='center', fontsize=12, fontweight='bold')

# Effect sizes
effect_sizes = [0] + [ttest_results[i]['cohen_d'] for i in ['exercise', 'diet', 'drug']]
axes[1].bar(intv_labels, effect_sizes, color=colors_intv, edgecolor='black')
axes[1].axhline(y=0, color='black', linewidth=1)
axes[1].axhline(y=0.2, color='gray', linestyle='--', alpha=0.5, label='Small effect')
axes[1].axhline(y=0.5, color='gray', linestyle='-.', alpha=0.5, label='Medium effect')
axes[1].axhline(y=0.8, color='gray', linestyle=':', alpha=0.5, label='Large effect')
axes[1].set_xlabel('Intervention Group', fontsize=12)
axes[1].set_ylabel("Cohen's d (vs. Control)", fontsize=12)
axes[1].set_title("Effect Sizes of Interventions\non Epigenetic Age Acceleration", fontsize=11)
axes[1].legend(fontsize=9)
for i, (label, es) in enumerate(zip(intv_labels, effect_sizes)):
    if es != 0:
        axes[1].text(i, es + 0.02 if es > 0 else es - 0.08, f'{es:.2f}', 
                      ha='center', fontsize=10)

plt.suptitle('Intervention Detection Sensitivity: Epigenetic Age Acceleration', 
              fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/fig04_intervention_effects.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/fig04_intervention_effects.png")


# ============================================================
# [cell:12] FIGURE 5: Neural network architecture comparison
# ============================================================
print("[cell:12] Generating Figure 5: Neural Network Architecture Comparison")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# All models performance summary (including NN variants)
all_results = {}
all_results.update(results)
for k, v in nn_results.items():
    all_results[k] = v

combined_names = list(results.keys()) + list(nn_results.keys())
combined_mae = [all_results[m]['mae_mean'] for m in combined_names]
combined_r2 = [all_results[m]['r2_mean'] for m in combined_names]
combined_mae_std = [all_results[m]['mae_std'] for m in combined_names]

short_names = ['ElasticNet', 'Lasso', 'Ridge', 'RF', 'XGB', 'LGB', 'MLP', 
               'MLP-S', 'MLP-M', 'MLP-L', 'MLP-D']

colors_all = ['#1565C0'] * 7 + ['#C62828'] * 4

bars = axes[0].bar(short_names, combined_mae, yerr=combined_mae_std, capsize=4,
                    color=colors_all, edgecolor='black', linewidth=0.5)
axes[0].set_xlabel('Model', fontsize=11)
axes[0].set_ylabel('MAE (years)', fontsize=11)
axes[0].set_title('All Models: Age Prediction MAE\n(Blue=Traditional, Red=NN variants)', fontsize=10)
axes[0].tick_params(axis='x', rotation=45)

bars2 = axes[1].bar(short_names, combined_r2, color=colors_all, edgecolor='black', linewidth=0.5)
axes[1].set_xlabel('Model', fontsize=11)
axes[1].set_ylabel('R² Score', fontsize=11)
axes[1].set_title('All Models: Age Prediction R²\n(Blue=Traditional, Red=NN variants)', fontsize=10)
axes[1].tick_params(axis='x', rotation=45)
axes[1].set_ylim(0, 1.1)

plt.suptitle('Comprehensive Model Comparison including Neural Network Architectures', 
              fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/fig05_nn_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/fig05_nn_comparison.png")


# ============================================================
# [cell:13] FIGURE 6: Longevity cohort validation
# ============================================================
print("[cell:13] Generating Figure 6: Longevity Cohort Validation")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Scatter: predicted vs chronological age for training and longevity cohorts
axes[0].scatter(y_chrono, results['ElasticNet (Horvath-like)']['y_pred'], 
                 alpha=0.4, s=15, c='steelblue', label=f'Training cohort (N={N_SAMPLES})')
axes[0].scatter(ages_long, pred_long, alpha=0.6, s=25, c='orange', marker='^',
                 label=f'Longevity cohort (N={N_LONG})')
axes[0].plot([15, 110], [15, 110], 'r--', linewidth=2, label='Perfect prediction')
axes[0].set_xlabel('Chronological Age (years)', fontsize=12)
axes[0].set_ylabel('Predicted Age (years)', fontsize=12)
axes[0].set_title(f'Longevity Cohort Validation\n(ElasticNet Clock)', fontsize=12)
axes[0].legend(fontsize=9)

# Age acceleration comparison
axes[1].hist(accel_train, bins=30, alpha=0.6, color='steelblue', 
              label=f'Training (μ={accel_train.mean():.2f})', density=True)
axes[1].hist(accel_long, bins=20, alpha=0.6, color='orange',
              label=f'Longevity (μ={accel_long.mean():.2f})', density=True)
axes[1].axvline(x=accel_train.mean(), color='steelblue', linestyle='--', linewidth=2)
axes[1].axvline(x=accel_long.mean(), color='orange', linestyle='--', linewidth=2)
axes[1].axvline(x=0, color='black', linestyle='-', linewidth=1)
axes[1].set_xlabel('Epigenetic Age Acceleration (years)', fontsize=12)
axes[1].set_ylabel('Density', fontsize=12)
axes[1].set_title(f'Age Acceleration: Training vs Longevity\n'
                   f'(t={t_long:.2f}, p={p_long_t:.4f})', fontsize=11)
axes[1].legend(fontsize=9)

plt.suptitle('Longevity Cohort Validation Strategy', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/fig06_longevity_validation.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/fig06_longevity_validation.png")


# ============================================================
# [cell:14] FIGURE 7: CpG importance heatmap
# ============================================================
print("[cell:14] Generating Figure 7: CpG Importance Analysis")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Top 30 CpGs importance
top30_indices = np.argsort(importances)[-30:][::-1]
top30_names = [cpg_names[i] for i in top30_indices]
top30_imp = importances[top30_indices]
top30_is_info = [int(cpg.replace('cg', '')) < N_INFORMATIVE for cpg in top30_names]

colors_cpg = ['#F44336' if is_info else '#9E9E9E' for is_info in top30_is_info]
axes[0].barh(range(30)[::-1], top30_imp, color=colors_cpg)
axes[0].set_yticks(range(30)[::-1])
axes[0].set_yticklabels(top30_names, fontsize=7)
axes[0].set_xlabel('Feature Importance (Random Forest)', fontsize=11)
axes[0].set_title('Top 30 CpG Sites by Importance\n(Red=Known age-correlated)', fontsize=11)

# Correlation of top CpGs with age
top10_data = X_features[:, top30_indices[:10]]
correlations = [pearsonr(top10_data[:, i], y_chrono)[0] for i in range(10)]
corr_colors = ['#F44336' if c > 0 else '#2196F3' for c in correlations]
axes[1].barh(range(10)[::-1], correlations, color=corr_colors)
axes[1].set_yticks(range(10)[::-1])
axes[1].set_yticklabels(top30_names[:10], fontsize=9)
axes[1].set_xlabel("Pearson r with Chronological Age", fontsize=11)
axes[1].set_title('Top 10 CpG Correlations with Age\n(Red=positive, Blue=negative)', fontsize=11)
axes[1].axvline(x=0, color='black', linewidth=1)

plt.suptitle('CpG Site Importance and Age Correlation Analysis', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/fig07_cpg_importance.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: figures/fig07_cpg_importance.png")


# ============================================================
# [cell:15] Summary Statistics Table
# ============================================================
print("\n" + "=" * 60)
print("[cell:15] FINAL SUMMARY STATISTICS")
print("=" * 60)

print("\n### Model Performance (5-fold CV) ###")
print(f"{'Model':<35} {'MAE':<12} {'R²':<12} {'Pearson r'}")
print("-" * 72)
for mname in full_names:
    r = results[mname]
    print(f"{mname:<35} {r['mae_mean']:.2f}±{r['mae_std']:.2f}   "
          f"{r['r2_mean']:.3f}±{r['r2_std']:.3f}   {r['pearson_r']:.3f}")

print("\n### Best Model ###")
best_idx = np.argmin([results[m]['mae_mean'] for m in full_names])
best_model = full_names[best_idx]
print(f"  Best model (by MAE): {best_model}")
print(f"  MAE: {results[best_model]['mae_mean']:.2f} ± {results[best_model]['mae_std']:.2f} years")
print(f"  R²: {results[best_model]['r2_mean']:.3f} ± {results[best_model]['r2_std']:.3f}")
print(f"  Pearson r: {results[best_model]['pearson_r']:.3f}")

print("\n### Neural Network Performance ###")
for arch_name in nn_archs.keys():
    r = nn_results[arch_name]
    print(f"  {arch_name:<20}: MAE={r['mae_mean']:.2f}±{r['mae_std']:.2f}  "
          f"R²={r['r2_mean']:.3f}±{r['r2_std']:.3f}  r={r['pearson_r']:.3f}")

print("\n### Intervention Effects ###")
for intv in ['exercise', 'diet', 'drug']:
    t = ttest_results[intv]
    sig = '***' if t['p'] < 0.001 else ('**' if t['p'] < 0.01 else ('*' if t['p'] < 0.05 else 'ns'))
    print(f"  {intv}: Cohen's d={t['cohen_d']:.3f}, p={t['p']:.4f} {sig}")

print(f"\n### Longevity Cohort ###")
print(f"  Training cohort age accel: {accel_train.mean():.2f} ± {accel_train.std():.2f}")
print(f"  Longevity cohort age accel: {accel_long.mean():.2f} ± {accel_long.std():.2f}")
print(f"  Difference significant: t={t_long:.3f}, p={p_long_t:.4f}")

print("\n### Generated Figures ###")
for fig_name in ['fig01_age_prediction_scatter.png', 'fig02_model_comparison.png',
                  'fig03_tissue_specific.png', 'fig04_intervention_effects.png',
                  'fig05_nn_comparison.png', 'fig06_longevity_validation.png',
                  'fig07_cpg_importance.png']:
    fpath = f'figures/{fig_name}'
    exists = os.path.exists(fpath)
    print(f"  {'✓' if exists else '✗'} {fpath}")

# Save summary to JSON for paper writing
import json
summary = {
    'n_samples': N_SAMPLES,
    'n_cpgs': N_CPGS,
    'n_informative': N_INFORMATIVE,
    'age_range': [float(ages.min()), float(ages.max())],
    'models': {k: {'mae_mean': float(v['mae_mean']), 'mae_std': float(v['mae_std']),
                    'r2_mean': float(v['r2_mean']), 'r2_std': float(v['r2_std']),
                    'pearson_r': float(v['pearson_r'])} 
                for k, v in results.items()},
    'nn_results': {k: {'mae_mean': float(v['mae_mean']), 'mae_std': float(v['mae_std']),
                        'r2_mean': float(v['r2_mean']), 'r2_std': float(v['r2_std']),
                        'pearson_r': float(v['pearson_r'])} 
                   for k, v in nn_results.items()},
    'intervention': {k: {'t': float(v['t']), 'p': float(v['p']), 'cohen_d': float(v['cohen_d'])}
                     for k, v in ttest_results.items()},
    'anova': {'F': float(f_stat), 'p': float(p_anova)},
    'longevity': {
        'accel_training': float(accel_train.mean()),
        'accel_training_std': float(accel_train.std()),
        'accel_longevity': float(accel_long.mean()),
        'accel_longevity_std': float(accel_long.std()),
        't': float(t_long),
        'p': float(p_long_t)
    },
    'tissue_results': {k: {'mae': float(v['mae']), 'mae_std': float(v['mae_std']),
                            'r2': float(v['r2']), 'pearson_r': float(v['pearson_r']),
                            'n': int(v['n'])}
                       for k, v in tissue_results.items()},
    'pan_tissue': {'mae_mean': float(np.mean(pan_mae)), 'mae_std': float(np.std(pan_mae)),
                   'r2_mean': float(np.mean(pan_r2)), 'pearson_r': float(corr_pan)},
    'informative_in_top20': int(informative_in_top20),
    'best_model': best_model,
    'intervention_accel': {k: {'mean': float(v['mean']), 'std': float(v['std']), 'n': int(v['n'])}
                            for k, v in intervention_accel.items()}
}

with open('data/raw/analysis_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print("\n✓ Analysis complete. Summary saved to data/raw/analysis_summary.json")
