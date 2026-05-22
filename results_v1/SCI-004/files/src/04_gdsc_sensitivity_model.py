"""
Module 4: Anticancer Drug Sensitivity Prediction (GDSC/CCLE-style)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.model_selection import cross_val_score, KFold, train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_regression
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ─────────────────────────────────────────────────────────────
# Simulate GDSC-like data
# n_cell_lines=200, n_drugs=20, n_genomic_features=500
# ─────────────────────────────────────────────────────────────
n_cell = 200
n_drugs = 20
n_genes_expr = 300
n_cnv = 100
n_mut = 50
n_features_total = n_genes_expr + n_cnv + n_mut

# Cell line genomic profiles
cancer_types = np.random.choice(
    ['BRCA','LUAD','COAD','SKCM','PAAD','GBM','PRAD','KIRC','OV','BLCA'],
    n_cell, p=[0.12,0.12,0.10,0.10,0.08,0.08,0.10,0.10,0.10,0.10]
)

gene_expression = np.random.randn(n_cell, n_genes_expr)
cnv_profiles    = np.random.choice([-2,-1,0,1,2], (n_cell, n_cnv),
                                    p=[0.05,0.10,0.70,0.10,0.05])
mutation_matrix = np.random.binomial(1, 0.05, (n_cell, n_mut))

# Known oncogene/driver gene effects
# EGFR expression → EGFR inhibitor sensitivity
egfr_idx = 0; kras_idx = 1; tp53_idx = 0; brca1_idx = 1
gene_expression[:, egfr_idx] += np.random.choice([0, 2, 4], n_cell,
                                                   p=[0.6, 0.25, 0.15])
mutation_matrix[:, tp53_idx] = np.random.binomial(1, 0.35, n_cell)
mutation_matrix[:, brca1_idx] = np.random.binomial(1, 0.12, n_cell)

X_genomic = np.hstack([gene_expression, cnv_profiles, mutation_matrix])

# Drug names and their primary targets
drugs = {
    'Erlotinib':    {'target': 'EGFR',   'pathway': 'RTK'},
    'Gefitinib':    {'target': 'EGFR',   'pathway': 'RTK'},
    'Vemurafenib':  {'target': 'BRAF',   'pathway': 'MAPK'},
    'Selumetinib':  {'target': 'MEK1/2', 'pathway': 'MAPK'},
    'Trametinib':   {'target': 'MEK1/2', 'pathway': 'MAPK'},
    'Olaparib':     {'target': 'PARP',   'pathway': 'DDR'},
    'Niraparib':    {'target': 'PARP',   'pathway': 'DDR'},
    'Gemcitabine':  {'target': 'RRM',    'pathway': 'Chemo'},
    'Paclitaxel':   {'target': 'Tubulin','pathway': 'Chemo'},
    'Doxorubicin':  {'target': 'TopoII', 'pathway': 'Chemo'},
    'Imatinib':     {'target': 'BCR-ABL','pathway': 'RTK'},
    'Dasatinib':    {'target': 'BCR-ABL','pathway': 'RTK'},
    'Palbociclib':  {'target': 'CDK4/6', 'pathway': 'CellCycle'},
    'Ribociclib':   {'target': 'CDK4/6', 'pathway': 'CellCycle'},
    'Everolimus':   {'target': 'mTOR',   'pathway': 'PI3K'},
    'Alpelisib':    {'target': 'PI3Kα',  'pathway': 'PI3K'},
    'Venetoclax':   {'target': 'BCL-2',  'pathway': 'Apoptosis'},
    'Navitoclax':   {'target': 'BCL-2',  'pathway': 'Apoptosis'},
    'Oxaliplatin':  {'target': 'DNA',    'pathway': 'Chemo'},
    'Carboplatin':  {'target': 'DNA',    'pathway': 'Chemo'},
}
drug_names = list(drugs.keys())

# Generate IC50 (log scale) for each drug × cell line
ln_ic50_matrix = np.zeros((n_cell, n_drugs))
for d_idx, drug_name in enumerate(drug_names):
    baseline_ic50 = np.random.normal(0, 1, n_cell)
    drug_info = drugs[drug_name]

    if drug_info['pathway'] == 'RTK' and 'EGFR' in drug_info['target']:
        # High EGFR expression → lower IC50 (sensitive)
        baseline_ic50 -= 0.8 * gene_expression[:, egfr_idx]
    if drug_info['pathway'] == 'DDR':
        # BRCA1 mutation → lower IC50
        baseline_ic50 -= 2.0 * mutation_matrix[:, brca1_idx]
    if drug_info['pathway'] == 'Chemo':
        # TP53 mutation → mild resistance
        baseline_ic50 += 0.5 * mutation_matrix[:, tp53_idx]
    # Cancer type effects
    brca_mask = (cancer_types == 'BRCA')
    if drug_info['pathway'] == 'DDR':
        baseline_ic50[brca_mask] -= 1.0

    ln_ic50_matrix[:, d_idx] = baseline_ic50 + np.random.randn(n_cell) * 0.5

df_ic50 = pd.DataFrame(ln_ic50_matrix, columns=drug_names)
df_ic50.insert(0, 'cell_line', [f'CL{i:03d}' for i in range(n_cell)])
df_ic50.insert(1, 'cancer_type', cancer_types)
df_ic50.to_csv('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/data/gdsc_ic50_synthetic.csv', index=False)

# ─────────────────────────────────────────────────────────────
# ML models for drug sensitivity prediction (per-drug)
# ─────────────────────────────────────────────────────────────
models_drug = {
    'Gradient Boosting': GradientBoostingRegressor(n_estimators=100, max_depth=3, random_state=42),
    'Random Forest':     RandomForestRegressor(n_estimators=100, random_state=42),
    'ElasticNet':        ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=2000),
    'Ridge':             Ridge(alpha=1.0),
}

drug_ml_results = {}
r2_summary = {m: [] for m in models_drug}

scaler_drug = StandardScaler()
X_scaled_drug = scaler_drug.fit_transform(X_genomic)

# Feature selection: top 50 features
selector = SelectKBest(f_regression, k=50)

target_drugs_eval = ['Erlotinib', 'Olaparib', 'Vemurafenib', 'Gemcitabine', 'Venetoclax']
for drug_name in target_drugs_eval:
    d_idx = drug_names.index(drug_name)
    y_drug = ln_ic50_matrix[:, d_idx]
    X_sel = selector.fit_transform(X_scaled_drug, y_drug)

    drug_ml_results[drug_name] = {}
    cv_kf = KFold(n_splits=5, shuffle=True, random_state=42)
    for m_name, model in models_drug.items():
        cv_r2  = cross_val_score(model, X_sel, y_drug, cv=cv_kf, scoring='r2')
        cv_mse = cross_val_score(model, X_sel, y_drug, cv=cv_kf, scoring='neg_mean_squared_error')
        drug_ml_results[drug_name][m_name] = {
            'r2_mean': float(cv_r2.mean()),
            'r2_std':  float(cv_r2.std()),
            'rmse_mean': float(np.sqrt(-cv_mse.mean())),
        }
        r2_summary[m_name].append(cv_r2.mean())

# ─────────────────────────────────────────────────────────────
# Figure 7: IC50 heatmap across cancer types
# ─────────────────────────────────────────────────────────────
ic50_by_cancer = df_ic50.groupby('cancer_type')[drug_names].mean()
fig, ax = plt.subplots(figsize=(14, 6))
sns.heatmap(ic50_by_cancer.T, cmap='RdBu_r', center=0, linewidths=0.3,
            ax=ax, cbar_kws={'label': 'ln(IC50)'}, annot=False)
ax.set_title('Drug Sensitivity (ln IC50) by Cancer Type\n(lower = more sensitive)')
ax.set_ylabel('Drug'); ax.set_xlabel('Cancer Type')
plt.tight_layout()
plt.savefig('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/figures/fig7_ic50_heatmap.png',
            dpi=150, bbox_inches='tight')
plt.close()

# ─────────────────────────────────────────────────────────────
# Figure 8: ML model R² comparison for selected drugs
# ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(11, 5))
x = np.arange(len(target_drugs_eval))
width = 0.2
bar_colors_m = ['#d73027','#4575b4','#1a9850','#762a83']

for i, (m_name, color) in enumerate(zip(models_drug.keys(), bar_colors_m)):
    r2_vals = [drug_ml_results[d][m_name]['r2_mean'] for d in target_drugs_eval]
    ax.bar(x + i*width, r2_vals, width, label=m_name, color=color, alpha=0.8)

ax.set_xticks(x + width*1.5)
ax.set_xticklabels(target_drugs_eval, rotation=20, ha='right')
ax.set_ylabel('Cross-validated R² Score')
ax.set_title('Drug Sensitivity Prediction: ML Model Comparison (5-fold CV)')
ax.legend(loc='upper right')
ax.axhline(0, color='gray', lw=0.5, linestyle='--')
plt.tight_layout()
plt.savefig('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/figures/fig8_drug_sensitivity_model.png',
            dpi=150, bbox_inches='tight')
plt.close()

# ─────────────────────────────────────────────────────────────
# Figure 9: Olaparib sensitivity vs BRCA1 mutation
# ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
brca1_status = mutation_matrix[:, brca1_idx].astype(str)
brca1_status = np.where(brca1_status == '1', 'BRCA1 Mutant', 'BRCA1 Wild-type')
olaparib_ic50 = ln_ic50_matrix[:, drug_names.index('Olaparib')]

df_plot = pd.DataFrame({'BRCA1_status': brca1_status, 'Olaparib_lnIC50': olaparib_ic50,
                         'cancer_type': cancer_types})
bp = axes[0].boxplot(
    [olaparib_ic50[brca1_status == 'BRCA1 Mutant'],
     olaparib_ic50[brca1_status == 'BRCA1 Wild-type']],
    labels=['BRCA1 Mutant', 'BRCA1 WT'], patch_artist=True
)
bp['boxes'][0].set_facecolor('#d73027'); bp['boxes'][0].set_alpha(0.7)
bp['boxes'][1].set_facecolor('#4575b4'); bp['boxes'][1].set_alpha(0.7)
axes[0].set_title('Olaparib Sensitivity: BRCA1 Mutation Status')
axes[0].set_ylabel('ln(IC50) — lower = more sensitive')

egfr_expr = gene_expression[:, egfr_idx]
erlotinib_ic50 = ln_ic50_matrix[:, drug_names.index('Erlotinib')]
axes[1].scatter(egfr_expr, erlotinib_ic50, alpha=0.4, c='#762a83', s=30)
z_fit = np.polyfit(egfr_expr, erlotinib_ic50, 1)
p_fit = np.poly1d(z_fit)
x_line = np.linspace(egfr_expr.min(), egfr_expr.max(), 100)
axes[1].plot(x_line, p_fit(x_line), 'r-', lw=2)
r_val = np.corrcoef(egfr_expr, erlotinib_ic50)[0,1]
axes[1].set_title(f'EGFR Expression vs Erlotinib Sensitivity (r={r_val:.3f})')
axes[1].set_xlabel('EGFR Expression Level')
axes[1].set_ylabel('ln(IC50)')

plt.tight_layout()
plt.savefig('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/figures/fig9_drug_biomarker.png',
            dpi=150, bbox_inches='tight')
plt.close()

results_gdsc = {
    'module': 'Anticancer Drug Sensitivity Prediction',
    'n_cell_lines': n_cell,
    'n_drugs': n_drugs,
    'n_genomic_features': n_features_total,
    'n_selected_features': 50,
    'drug_ml_results': drug_ml_results,
    'mean_r2_by_model': {m: float(np.mean(v)) for m, v in r2_summary.items() if v},
}
with open('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/results/gdsc_sensitivity_results.json', 'w') as f:
    json.dump(results_gdsc, f, indent=2, ensure_ascii=False)

print("[GDSC Module] Done")
for drug in target_drugs_eval:
    best = max(drug_ml_results[drug].items(), key=lambda x: x[1]['r2_mean'])
    print(f"  {drug}: best={best[0]} R²={best[1]['r2_mean']:.3f}")
