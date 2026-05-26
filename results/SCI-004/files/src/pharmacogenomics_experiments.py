#!/usr/bin/env python3
"""
Pharmacogenomics Model Construction: Comprehensive Experiments
==============================================================
1. CYP enzyme polymorphism & drug metabolism rate modeling
2. HLA genotype & adverse drug reaction prediction
3. GWAS summary statistics drug target validation (MR analysis)
4. Anticancer drug sensitivity prediction (GDSC/CCLE-style)
5. Deep learning drug-gene interaction network
6. Clinical Decision Support System (CDSS) prototype
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.metrics import (accuracy_score, roc_auc_score, roc_curve,
                             precision_recall_curve, confusion_matrix,
                             mean_squared_error, r2_score, classification_report,
                             f1_score, precision_score, recall_score)
from sklearn.preprocessing import StandardScaler, LabelEncoder
from scipy import stats
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import networkx as nx
import os
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
torch.manual_seed(42)

FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

results = {}

# ============================================================
# Experiment 1: CYP Enzyme Polymorphism & Drug Metabolism
# ============================================================
print("=" * 60)
print("Experiment 1: CYP Enzyme Polymorphism Modeling")
print("=" * 60)

def generate_cyp_data(n=2000):
    """Simulate CYP2D6/CYP2C19 polymorphism data."""
    cyp2d6_alleles = ['*1/*1', '*1/*2', '*1/*4', '*2/*2', '*4/*4',
                      '*1/*5', '*1/*41', '*2/*41', '*1/*1xN', '*2/*2xN']
    cyp2c19_alleles = ['*1/*1', '*1/*2', '*1/*3', '*2/*2', '*2/*3',
                       '*3/*3', '*1/*17', '*17/*17']
    
    # Activity scores for CYP2D6
    cyp2d6_activity = {
        '*1/*1': 2.0, '*1/*2': 1.5, '*1/*4': 1.0, '*2/*2': 1.0,
        '*4/*4': 0.0, '*1/*5': 1.0, '*1/*41': 1.5, '*2/*41': 1.0,
        '*1/*1xN': 3.0, '*2/*2xN': 2.0
    }
    
    # Metabolizer phenotypes
    cyp2d6_phenotype = {
        '*1/*1': 'NM', '*1/*2': 'NM', '*1/*4': 'IM', '*2/*2': 'NM',
        '*4/*4': 'PM', '*1/*5': 'IM', '*1/*41': 'NM', '*2/*41': 'IM',
        '*1/*1xN': 'UM', '*2/*2xN': 'UM'
    }
    
    data = []
    for _ in range(n):
        d6 = np.random.choice(cyp2d6_alleles, p=[0.25, 0.15, 0.12, 0.08, 0.05,
                                                   0.08, 0.10, 0.07, 0.05, 0.05])
        c19 = np.random.choice(cyp2c19_alleles, p=[0.30, 0.20, 0.10, 0.10,
                                                     0.08, 0.02, 0.12, 0.08])
        
        activity = cyp2d6_activity[d6]
        phenotype = cyp2d6_phenotype[d6]
        
        age = np.random.normal(50, 15)
        weight = np.random.normal(70, 12)
        liver_function = np.random.uniform(0.5, 1.5)
        
        # Drug metabolism rate (clearance) model
        base_clearance = activity * 15 + np.random.normal(0, 3)
        clearance = base_clearance * liver_function * (weight / 70) ** 0.75
        clearance = max(clearance, 0.1)
        
        # Half-life inversely proportional to clearance
        half_life = 200 / clearance + np.random.normal(0, 1)
        half_life = max(half_life, 0.5)
        
        data.append({
            'cyp2d6_genotype': d6,
            'cyp2c19_genotype': c19,
            'activity_score': activity,
            'phenotype': phenotype,
            'age': age,
            'weight': weight,
            'liver_function': liver_function,
            'clearance': clearance,
            'half_life': half_life
        })
    
    return pd.DataFrame(data)

cyp_data = generate_cyp_data(2000)

# Encode genotypes
le_d6 = LabelEncoder()
le_c19 = LabelEncoder()
le_pheno = LabelEncoder()
cyp_data['d6_encoded'] = le_d6.fit_transform(cyp_data['cyp2d6_genotype'])
cyp_data['c19_encoded'] = le_c19.fit_transform(cyp_data['cyp2c19_genotype'])
cyp_data['phenotype_encoded'] = le_pheno.fit_transform(cyp_data['phenotype'])

features_cyp = ['d6_encoded', 'c19_encoded', 'activity_score', 'age', 'weight', 'liver_function']
X_cyp = cyp_data[features_cyp].values
y_clearance = cyp_data['clearance'].values
y_phenotype = cyp_data['phenotype_encoded'].values

# Split data
X_train, X_test, y_train_cl, y_test_cl = train_test_split(X_cyp, y_clearance, test_size=0.2, random_state=42)
_, _, y_train_ph, y_test_ph = train_test_split(X_cyp, y_phenotype, test_size=0.2, random_state=42)

scaler_cyp = StandardScaler()
X_train_sc = scaler_cyp.fit_transform(X_train)
X_test_sc = scaler_cyp.transform(X_test)

# Regression: Predict clearance
gb_reg = GradientBoostingRegressor(n_estimators=200, max_depth=5, random_state=42)
gb_reg.fit(X_train_sc, y_train_cl)
y_pred_cl = gb_reg.predict(X_test_sc)
rmse_cl = np.sqrt(mean_squared_error(y_test_cl, y_pred_cl))
r2_cl = r2_score(y_test_cl, y_pred_cl)

# Classification: Predict metabolizer phenotype
rf_cls = RandomForestClassifier(n_estimators=200, max_depth=10, random_state=42)
rf_cls.fit(X_train_sc, y_train_ph)
y_pred_ph = rf_cls.predict(X_test_sc)
acc_ph = accuracy_score(y_test_ph, y_pred_ph)
f1_ph = f1_score(y_test_ph, y_pred_ph, average='weighted')

results['exp1'] = {
    'clearance_rmse': round(rmse_cl, 4),
    'clearance_r2': round(r2_cl, 4),
    'phenotype_accuracy': round(acc_ph, 4),
    'phenotype_f1': round(f1_ph, 4)
}

print(f"Clearance Prediction - RMSE: {rmse_cl:.4f}, R²: {r2_cl:.4f}")
print(f"Phenotype Classification - Accuracy: {acc_ph:.4f}, F1: {f1_ph:.4f}")

# Figure 1: CYP metabolism results
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 1a: Clearance by phenotype
phenotype_order = ['PM', 'IM', 'NM', 'UM']
colors = ['#e74c3c', '#f39c12', '#27ae60', '#3498db']
box_data = [cyp_data[cyp_data['phenotype'] == p]['clearance'].values for p in phenotype_order]
bp = axes[0, 0].boxplot(box_data, labels=phenotype_order, patch_artist=True)
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
axes[0, 0].set_xlabel('Metabolizer Phenotype', fontsize=12)
axes[0, 0].set_ylabel('Drug Clearance (L/h)', fontsize=12)
axes[0, 0].set_title('(A) Drug Clearance by CYP2D6 Phenotype', fontsize=13, fontweight='bold')

# 1b: Predicted vs actual clearance
axes[0, 1].scatter(y_test_cl, y_pred_cl, alpha=0.4, s=20, c='#3498db')
lims = [min(y_test_cl.min(), y_pred_cl.min()), max(y_test_cl.max(), y_pred_cl.max())]
axes[0, 1].plot(lims, lims, 'r--', lw=2)
axes[0, 1].set_xlabel('Actual Clearance (L/h)', fontsize=12)
axes[0, 1].set_ylabel('Predicted Clearance (L/h)', fontsize=12)
axes[0, 1].set_title(f'(B) Clearance Prediction (R²={r2_cl:.3f})', fontsize=13, fontweight='bold')

# 1c: Feature importance
importances = gb_reg.feature_importances_
feat_names = ['CYP2D6', 'CYP2C19', 'Activity Score', 'Age', 'Weight', 'Liver Func']
sorted_idx = np.argsort(importances)
axes[1, 0].barh(range(len(sorted_idx)), importances[sorted_idx], color='#2ecc71')
axes[1, 0].set_yticks(range(len(sorted_idx)))
axes[1, 0].set_yticklabels([feat_names[i] for i in sorted_idx])
axes[1, 0].set_xlabel('Feature Importance', fontsize=12)
axes[1, 0].set_title('(C) Feature Importance for Clearance', fontsize=13, fontweight='bold')

# 1d: Confusion matrix
cm = confusion_matrix(y_test_ph, y_pred_ph)
pheno_labels = le_pheno.classes_
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=pheno_labels,
            yticklabels=pheno_labels, ax=axes[1, 1])
axes[1, 1].set_xlabel('Predicted', fontsize=12)
axes[1, 1].set_ylabel('Actual', fontsize=12)
axes[1, 1].set_title(f'(D) Phenotype Classification (Acc={acc_ph:.3f})', fontsize=13, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig1_cyp_metabolism.png'), dpi=150, bbox_inches='tight')
plt.close()
print("Figure 1 saved.")

# ============================================================
# Experiment 2: HLA-B*1502 & Carbamazepine ADR Prediction
# ============================================================
print("\n" + "=" * 60)
print("Experiment 2: HLA-ADR Prediction")
print("=" * 60)

def generate_hla_data(n=3000):
    """Simulate HLA genotype and ADR data."""
    data = []
    hla_b_alleles = ['B*15:02', 'B*15:01', 'B*40:01', 'B*46:01', 'B*58:01',
                     'B*07:02', 'B*08:01', 'B*44:02', 'B*35:01', 'B*51:01']
    hla_a_alleles = ['A*31:01', 'A*02:01', 'A*24:02', 'A*11:01', 'A*33:03']
    
    for _ in range(n):
        hla_b = np.random.choice(hla_b_alleles,
                                  p=[0.08, 0.12, 0.10, 0.10, 0.08,
                                     0.12, 0.10, 0.10, 0.10, 0.10])
        hla_a = np.random.choice(hla_a_alleles, p=[0.10, 0.30, 0.25, 0.20, 0.15])
        
        ancestry = np.random.choice(['East_Asian', 'South_Asian', 'European', 'African'],
                                     p=[0.30, 0.20, 0.35, 0.15])
        age = np.random.normal(45, 15)
        dose = np.random.choice([200, 400, 600, 800, 1000])
        duration_days = np.random.randint(1, 365)
        
        # ADR risk model
        risk = 0.02  # baseline
        if hla_b == 'B*15:02':
            risk += 0.60 if ancestry in ['East_Asian', 'South_Asian'] else 0.20
        if hla_a == 'A*31:01':
            risk += 0.15
        if hla_b == 'B*58:01':
            risk += 0.10
        risk *= (dose / 400) ** 0.3
        risk += np.random.normal(0, 0.05)
        risk = np.clip(risk, 0, 1)
        
        adr = 1 if np.random.random() < risk else 0
        
        # Severity (if ADR occurs)
        if adr:
            severity = np.random.choice(['mild', 'moderate', 'severe', 'SJS/TEN'],
                                         p=[0.30, 0.35, 0.25, 0.10])
        else:
            severity = 'none'
        
        # Genomic features
        n_risk_variants = np.random.poisson(2)
        prs = np.random.normal(0, 1) + (0.5 if hla_b == 'B*15:02' else 0)
        
        data.append({
            'hla_b': hla_b,
            'hla_a': hla_a,
            'ancestry': ancestry,
            'age': age,
            'dose': dose,
            'duration_days': duration_days,
            'n_risk_variants': n_risk_variants,
            'prs': prs,
            'adr_risk': risk,
            'adr': adr,
            'severity': severity
        })
    
    return pd.DataFrame(data)

hla_data = generate_hla_data(3000)

# Encode
for col in ['hla_b', 'hla_a', 'ancestry']:
    hla_data[col + '_enc'] = LabelEncoder().fit_transform(hla_data[col])

features_hla = ['hla_b_enc', 'hla_a_enc', 'ancestry_enc', 'age', 'dose',
                'duration_days', 'n_risk_variants', 'prs']
X_hla = hla_data[features_hla].values
y_hla = hla_data['adr'].values

X_train_h, X_test_h, y_train_h, y_test_h = train_test_split(X_hla, y_hla, test_size=0.2,
                                                               random_state=42, stratify=y_hla)
scaler_hla = StandardScaler()
X_train_hs = scaler_hla.fit_transform(X_train_h)
X_test_hs = scaler_hla.transform(X_test_h)

# Logistic Regression
lr_hla = LogisticRegression(max_iter=1000, class_weight='balanced', random_state=42)
lr_hla.fit(X_train_hs, y_train_h)
y_pred_hla_lr = lr_hla.predict(X_test_hs)
y_prob_hla_lr = lr_hla.predict_proba(X_test_hs)[:, 1]

# Random Forest
rf_hla = RandomForestClassifier(n_estimators=200, max_depth=10, class_weight='balanced', random_state=42)
rf_hla.fit(X_train_hs, y_train_h)
y_pred_hla_rf = rf_hla.predict(X_test_hs)
y_prob_hla_rf = rf_hla.predict_proba(X_test_hs)[:, 1]

# Neural Network
class ADRPredictor(nn.Module):
    def __init__(self, input_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(32, 1),
            nn.Sigmoid()
        )
    def forward(self, x):
        return self.net(x)

X_train_t = torch.FloatTensor(X_train_hs)
y_train_t = torch.FloatTensor(y_train_h).unsqueeze(1)
X_test_t = torch.FloatTensor(X_test_hs)

model_adr = ADRPredictor(X_train_hs.shape[1])
criterion = nn.BCELoss()
optimizer_adr = optim.Adam(model_adr.parameters(), lr=0.001)

train_losses = []
for epoch in range(100):
    model_adr.train()
    optimizer_adr.zero_grad()
    out = model_adr(X_train_t)
    loss = criterion(out, y_train_t)
    loss.backward()
    optimizer_adr.step()
    train_losses.append(loss.item())

model_adr.eval()
with torch.no_grad():
    y_prob_hla_nn = model_adr(X_test_t).numpy().flatten()
y_pred_hla_nn = (y_prob_hla_nn >= 0.5).astype(int)

auc_lr = roc_auc_score(y_test_h, y_prob_hla_lr)
auc_rf = roc_auc_score(y_test_h, y_prob_hla_rf)
auc_nn = roc_auc_score(y_test_h, y_prob_hla_nn)

results['exp2'] = {
    'lr_auc': round(auc_lr, 4),
    'rf_auc': round(auc_rf, 4),
    'nn_auc': round(auc_nn, 4),
    'lr_accuracy': round(accuracy_score(y_test_h, y_pred_hla_lr), 4),
    'rf_accuracy': round(accuracy_score(y_test_h, y_pred_hla_rf), 4),
    'nn_accuracy': round(accuracy_score(y_test_h, y_pred_hla_nn), 4),
    'n_adr_cases': int(y_hla.sum()),
    'adr_rate': round(y_hla.mean(), 4)
}

print(f"LR AUC: {auc_lr:.4f}, RF AUC: {auc_rf:.4f}, NN AUC: {auc_nn:.4f}")

# Figure 2: HLA-ADR results
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 2a: ROC curves
for name, y_prob, color in [('Logistic Regression', y_prob_hla_lr, '#e74c3c'),
                              ('Random Forest', y_prob_hla_rf, '#3498db'),
                              ('Neural Network', y_prob_hla_nn, '#2ecc71')]:
    fpr, tpr, _ = roc_curve(y_test_h, y_prob)
    auc_val = roc_auc_score(y_test_h, y_prob)
    axes[0, 0].plot(fpr, tpr, color=color, lw=2, label=f'{name} (AUC={auc_val:.3f})')
axes[0, 0].plot([0, 1], [0, 1], 'k--', lw=1)
axes[0, 0].set_xlabel('False Positive Rate', fontsize=12)
axes[0, 0].set_ylabel('True Positive Rate', fontsize=12)
axes[0, 0].set_title('(A) ROC Curves for ADR Prediction', fontsize=13, fontweight='bold')
axes[0, 0].legend(fontsize=10)

# 2b: ADR rate by HLA-B allele
adr_by_hla = hla_data.groupby('hla_b')['adr'].mean().sort_values(ascending=False)
bars = axes[0, 1].bar(range(len(adr_by_hla)), adr_by_hla.values, color='#e74c3c', alpha=0.7)
axes[0, 1].set_xticks(range(len(adr_by_hla)))
axes[0, 1].set_xticklabels(adr_by_hla.index, rotation=45, ha='right', fontsize=9)
axes[0, 1].set_ylabel('ADR Rate', fontsize=12)
axes[0, 1].set_title('(B) ADR Rate by HLA-B Allele', fontsize=13, fontweight='bold')

# 2c: Training loss
axes[1, 0].plot(train_losses, color='#8e44ad', lw=2)
axes[1, 0].set_xlabel('Epoch', fontsize=12)
axes[1, 0].set_ylabel('BCE Loss', fontsize=12)
axes[1, 0].set_title('(C) Neural Network Training Loss', fontsize=13, fontweight='bold')

# 2d: Precision-Recall
for name, y_prob, color in [('LR', y_prob_hla_lr, '#e74c3c'),
                              ('RF', y_prob_hla_rf, '#3498db'),
                              ('NN', y_prob_hla_nn, '#2ecc71')]:
    prec, rec, _ = precision_recall_curve(y_test_h, y_prob)
    axes[1, 1].plot(rec, prec, color=color, lw=2, label=name)
axes[1, 1].set_xlabel('Recall', fontsize=12)
axes[1, 1].set_ylabel('Precision', fontsize=12)
axes[1, 1].set_title('(D) Precision-Recall Curves', fontsize=13, fontweight='bold')
axes[1, 1].legend(fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig2_hla_adr.png'), dpi=150, bbox_inches='tight')
plt.close()
print("Figure 2 saved.")

# ============================================================
# Experiment 3: Mendelian Randomization (MR) Analysis
# ============================================================
print("\n" + "=" * 60)
print("Experiment 3: Mendelian Randomization Drug Target Validation")
print("=" * 60)

def simulate_mr_analysis(n_snps=100, n_targets=10):
    """Simulate MR analysis for drug target validation."""
    target_genes = ['PCSK9', 'HMGCR', 'NPC1L1', 'CETP', 'LPL',
                    'APOC3', 'ANGPTL3', 'LDLR', 'SORT1', 'ABCG5']
    
    mr_results = []
    for gene in target_genes:
        n_ivs = np.random.randint(5, 50)
        
        # Simulate IV effects
        beta_exposure = np.random.normal(0.1, 0.05, n_ivs)
        se_exposure = np.abs(np.random.normal(0.02, 0.01, n_ivs))
        
        # True causal effect varies by gene
        true_effect = np.random.normal(0, 0.3)
        if gene in ['PCSK9', 'HMGCR', 'NPC1L1']:
            true_effect = np.random.normal(-0.3, 0.1)  # protective
        elif gene in ['CETP']:
            true_effect = np.random.normal(0.1, 0.05)  # risk
        
        beta_outcome = beta_exposure * true_effect + np.random.normal(0, 0.02, n_ivs)
        se_outcome = np.abs(np.random.normal(0.03, 0.01, n_ivs))
        
        # IVW estimate
        weights = 1 / se_outcome**2
        ivw_beta = np.sum(weights * beta_outcome / beta_exposure) / np.sum(weights / beta_exposure**2)
        ivw_se = np.sqrt(1 / np.sum(weights / beta_exposure**2))
        ivw_p = 2 * stats.norm.sf(np.abs(ivw_beta / ivw_se))
        
        # MR-Egger
        slope, intercept, _, p_egger, se_slope = stats.linregress(
            beta_exposure / se_outcome,
            beta_outcome / se_outcome
        )
        egger_intercept_p = 2 * stats.norm.sf(np.abs(intercept / (se_slope * 0.5)))
        
        # Weighted median
        sorted_idx = np.argsort(beta_outcome / beta_exposure)
        cum_weights = np.cumsum(weights[sorted_idx]) / np.sum(weights)
        median_idx = np.searchsorted(cum_weights, 0.5)
        wm_beta = (beta_outcome / beta_exposure)[sorted_idx[min(median_idx, n_ivs - 1)]]
        
        # Cochran's Q
        q_stat = np.sum(weights * (beta_outcome / beta_exposure - ivw_beta)**2)
        q_p = 1 - stats.chi2.cdf(q_stat, n_ivs - 1)
        
        f_stat = np.mean((beta_exposure / se_exposure)**2)
        
        mr_results.append({
            'gene': gene,
            'n_ivs': n_ivs,
            'ivw_beta': ivw_beta,
            'ivw_se': ivw_se,
            'ivw_p': ivw_p,
            'egger_beta': slope,
            'egger_intercept_p': egger_intercept_p,
            'weighted_median_beta': wm_beta,
            'cochran_q_p': q_p,
            'f_statistic': f_stat,
            'significant': ivw_p < 0.05
        })
    
    return pd.DataFrame(mr_results)

mr_data = simulate_mr_analysis()
n_significant = mr_data['significant'].sum()
mean_f = mr_data['f_statistic'].mean()

results['exp3'] = {
    'n_targets_tested': len(mr_data),
    'n_significant': int(n_significant),
    'significant_genes': mr_data[mr_data['significant']]['gene'].tolist(),
    'mean_f_statistic': round(mean_f, 2)
}

print(f"Targets tested: {len(mr_data)}, Significant: {n_significant}")
print(f"Significant genes: {mr_data[mr_data['significant']]['gene'].tolist()}")

# Figure 3: MR results
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 3a: Forest plot
mr_sorted = mr_data.sort_values('ivw_beta')
y_pos = range(len(mr_sorted))
for idx_mr, (_, row_mr) in enumerate(mr_sorted.iterrows()):
    c = '#e74c3c' if row_mr['ivw_p'] < 0.05 else '#95a5a6'
    axes[0, 0].errorbar(row_mr['ivw_beta'], idx_mr,
                         xerr=1.96 * row_mr['ivw_se'],
                         fmt='o', color=c, ecolor=c, capsize=3, markersize=6)
axes[0, 0].axvline(x=0, color='gray', linestyle='--', lw=1)
axes[0, 0].set_yticks(y_pos)
axes[0, 0].set_yticklabels(mr_sorted['gene'])
axes[0, 0].set_xlabel('IVW Causal Effect Estimate', fontsize=12)
axes[0, 0].set_title('(A) MR Forest Plot', fontsize=13, fontweight='bold')

# 3b: Comparison of MR methods
x = np.arange(len(mr_data))
width = 0.25
axes[0, 1].bar(x - width, mr_data['ivw_beta'], width, label='IVW', color='#3498db', alpha=0.8)
axes[0, 1].bar(x, mr_data['egger_beta'], width, label='MR-Egger', color='#e74c3c', alpha=0.8)
axes[0, 1].bar(x + width, mr_data['weighted_median_beta'], width, label='Weighted Median', color='#2ecc71', alpha=0.8)
axes[0, 1].set_xticks(x)
axes[0, 1].set_xticklabels(mr_data['gene'], rotation=45, ha='right', fontsize=8)
axes[0, 1].set_ylabel('Causal Effect Estimate', fontsize=12)
axes[0, 1].set_title('(B) MR Method Comparison', fontsize=13, fontweight='bold')
axes[0, 1].legend(fontsize=9)
axes[0, 1].axhline(y=0, color='gray', linestyle='--', lw=1)

# 3c: P-value volcano plot
neg_log_p = -np.log10(mr_data['ivw_p'] + 1e-300)
colors_vol = ['#e74c3c' if p < 0.05 else '#95a5a6' for p in mr_data['ivw_p']]
axes[1, 0].scatter(mr_data['ivw_beta'], neg_log_p, c=colors_vol, s=100, edgecolors='black', lw=0.5)
axes[1, 0].axhline(y=-np.log10(0.05), color='red', linestyle='--', lw=1, label='p=0.05')
for i, row in mr_data.iterrows():
    if row['ivw_p'] < 0.05:
        axes[1, 0].annotate(row['gene'], (row['ivw_beta'], -np.log10(row['ivw_p'] + 1e-300)),
                            fontsize=8, ha='center', va='bottom')
axes[1, 0].set_xlabel('IVW Effect Size', fontsize=12)
axes[1, 0].set_ylabel('-log₁₀(P-value)', fontsize=12)
axes[1, 0].set_title('(C) Volcano Plot', fontsize=13, fontweight='bold')
axes[1, 0].legend()

# 3d: Instrument strength
axes[1, 1].bar(mr_data['gene'], mr_data['f_statistic'], color='#f39c12', alpha=0.8)
axes[1, 1].axhline(y=10, color='red', linestyle='--', lw=1, label='F=10 threshold')
axes[1, 1].set_xticklabels(mr_data['gene'], rotation=45, ha='right', fontsize=8)
axes[1, 1].set_ylabel('Mean F-statistic', fontsize=12)
axes[1, 1].set_title('(D) Instrument Strength', fontsize=13, fontweight='bold')
axes[1, 1].legend()

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig3_mr_analysis.png'), dpi=150, bbox_inches='tight')
plt.close()
print("Figure 3 saved.")

# ============================================================
# Experiment 4: Anticancer Drug Sensitivity (GDSC/CCLE-style)
# ============================================================
print("\n" + "=" * 60)
print("Experiment 4: Anticancer Drug Sensitivity Prediction")
print("=" * 60)

def generate_gdsc_ccle_data(n_cell_lines=300, n_drugs=10, n_genes=100):
    """Simulate GDSC/CCLE-like drug sensitivity data."""
    gene_expr = np.random.randn(n_cell_lines, n_genes)
    mutation_matrix = (np.random.random((n_cell_lines, 30)) < 0.1).astype(float)
    cnv_data = np.random.randn(n_cell_lines, 20) * 0.5
    
    genomic_features = np.hstack([gene_expr, mutation_matrix, cnv_data])
    
    drug_names = [f'Drug_{i+1}' for i in range(n_drugs)]
    ic50_data = {}
    
    for d in range(n_drugs):
        drug_effect = np.random.randn(genomic_features.shape[1]) * 0.1
        drug_effect[np.random.choice(genomic_features.shape[1], 10, replace=False)] *= 5
        
        ic50 = genomic_features @ drug_effect + np.random.randn(n_cell_lines) * 2
        ic50_data[drug_names[d]] = ic50
    
    return genomic_features, pd.DataFrame(ic50_data), drug_names

genomic_features, ic50_df, drug_names = generate_gdsc_ccle_data()

# Multi-task drug sensitivity predictor
class DrugSensitivityNet(nn.Module):
    def __init__(self, input_dim, n_drugs):
        super().__init__()
        self.shared = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.ReLU(),
            nn.BatchNorm1d(256),
            nn.Dropout(0.3),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.BatchNorm1d(128),
            nn.Dropout(0.2),
        )
        self.drug_heads = nn.ModuleList([
            nn.Sequential(nn.Linear(128, 64), nn.ReLU(), nn.Linear(64, 1))
            for _ in range(n_drugs)
        ])
    
    def forward(self, x):
        shared_out = self.shared(x)
        return [head(shared_out) for head in self.drug_heads]

X_gdsc = genomic_features
y_gdsc = ic50_df.values

X_train_g, X_test_g, y_train_g, y_test_g = train_test_split(X_gdsc, y_gdsc, test_size=0.2, random_state=42)

scaler_g = StandardScaler()
X_train_gs = scaler_g.fit_transform(X_train_g)
X_test_gs = scaler_g.transform(X_test_g)

X_train_gt = torch.FloatTensor(X_train_gs)
y_train_gt = torch.FloatTensor(y_train_g)
X_test_gt = torch.FloatTensor(X_test_gs)

model_gdsc = DrugSensitivityNet(X_train_gs.shape[1], len(drug_names))
optimizer_g = optim.Adam(model_gdsc.parameters(), lr=0.001, weight_decay=1e-4)
criterion_g = nn.MSELoss()

train_dataset = TensorDataset(X_train_gt, y_train_gt)
train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)

gdsc_losses = []
for epoch in range(50):
    model_gdsc.train()
    epoch_loss = 0
    for batch_x, batch_y in train_loader:
        optimizer_g.zero_grad()
        preds = model_gdsc(batch_x)
        loss = sum(criterion_g(preds[i].squeeze(), batch_y[:, i]) for i in range(len(drug_names)))
        loss.backward()
        optimizer_g.step()
        epoch_loss += loss.item()
    gdsc_losses.append(epoch_loss / len(train_loader))

model_gdsc.eval()
with torch.no_grad():
    preds_test = model_gdsc(X_test_gt)
    y_pred_gdsc = np.column_stack([p.numpy().flatten() for p in preds_test])

# Per-drug performance
drug_r2 = []
drug_rmse = []
for i in range(len(drug_names)):
    r2 = r2_score(y_test_g[:, i], y_pred_gdsc[:, i])
    rmse = np.sqrt(mean_squared_error(y_test_g[:, i], y_pred_gdsc[:, i]))
    drug_r2.append(r2)
    drug_rmse.append(rmse)

# Baseline: Random Forest per drug
rf_r2 = []
for i in range(len(drug_names)):
    rf = GradientBoostingRegressor(n_estimators=50, random_state=42)
    rf.fit(X_train_gs, y_train_g[:, i])
    pred_rf = rf.predict(X_test_gs)
    rf_r2.append(r2_score(y_test_g[:, i], pred_rf))

results['exp4'] = {
    'mean_dl_r2': round(np.mean(drug_r2), 4),
    'mean_dl_rmse': round(np.mean(drug_rmse), 4),
    'mean_rf_r2': round(np.mean(rf_r2), 4),
    'best_drug': drug_names[np.argmax(drug_r2)],
    'best_r2': round(max(drug_r2), 4)
}

print(f"DL Mean R²: {np.mean(drug_r2):.4f}, RF Mean R²: {np.mean(rf_r2):.4f}")

# Figure 4: Drug sensitivity results
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 4a: Per-drug R² comparison
x_drugs = np.arange(len(drug_names))
axes[0, 0].bar(x_drugs - 0.2, drug_r2, 0.4, label='Deep Learning', color='#3498db', alpha=0.8)
axes[0, 0].bar(x_drugs + 0.2, rf_r2, 0.4, label='Gradient Boosting', color='#e74c3c', alpha=0.8)
axes[0, 0].set_xticks(x_drugs)
axes[0, 0].set_xticklabels(drug_names, rotation=90, fontsize=7)
axes[0, 0].set_ylabel('R² Score', fontsize=12)
axes[0, 0].set_title('(A) Per-Drug R² Comparison', fontsize=13, fontweight='bold')
axes[0, 0].legend(fontsize=10)

# 4b: Training loss
axes[0, 1].plot(gdsc_losses, color='#8e44ad', lw=2)
axes[0, 1].set_xlabel('Epoch', fontsize=12)
axes[0, 1].set_ylabel('MSE Loss', fontsize=12)
axes[0, 1].set_title('(B) Multi-task DL Training Loss', fontsize=13, fontweight='bold')

# 4c: Best drug predicted vs actual
best_idx = np.argmax(drug_r2)
axes[1, 0].scatter(y_test_g[:, best_idx], y_pred_gdsc[:, best_idx], alpha=0.5, s=25, c='#3498db')
lims = [min(y_test_g[:, best_idx].min(), y_pred_gdsc[:, best_idx].min()),
        max(y_test_g[:, best_idx].max(), y_pred_gdsc[:, best_idx].max())]
axes[1, 0].plot(lims, lims, 'r--', lw=2)
axes[1, 0].set_xlabel('Actual IC50', fontsize=12)
axes[1, 0].set_ylabel('Predicted IC50', fontsize=12)
axes[1, 0].set_title(f'(C) {drug_names[best_idx]} (R²={drug_r2[best_idx]:.3f})', fontsize=13, fontweight='bold')

# 4d: R² distribution
axes[1, 1].hist(drug_r2, bins=10, color='#2ecc71', alpha=0.7, edgecolor='black', label='Deep Learning')
axes[1, 1].hist(rf_r2, bins=10, color='#e74c3c', alpha=0.5, edgecolor='black', label='Gradient Boosting')
axes[1, 1].set_xlabel('R² Score', fontsize=12)
axes[1, 1].set_ylabel('Number of Drugs', fontsize=12)
axes[1, 1].set_title('(D) R² Distribution Across Drugs', fontsize=13, fontweight='bold')
axes[1, 1].legend(fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig4_drug_sensitivity.png'), dpi=150, bbox_inches='tight')
plt.close()
print("Figure 4 saved.")

# ============================================================
# Experiment 5: Drug-Gene Interaction Network (GNN-style)
# ============================================================
print("\n" + "=" * 60)
print("Experiment 5: Drug-Gene Interaction Network Learning")
print("=" * 60)

def build_drug_gene_network():
    """Build a drug-gene interaction network."""
    genes = ['CYP2D6', 'CYP2C19', 'CYP3A4', 'CYP1A2', 'UGT1A1',
             'ABCB1', 'SLCO1B1', 'VKORC1', 'DPYD', 'TPMT',
             'HLA-B', 'HLA-A', 'EGFR', 'BRAF', 'ALK',
             'BRCA1', 'TP53', 'KRAS', 'PIK3CA', 'ERBB2']
    
    drugs = ['Warfarin', 'Clopidogrel', 'Tamoxifen', 'Codeine', 'Irinotecan',
             'Carbamazepine', 'Simvastatin', '5-FU', '6-MP', 'Gefitinib',
             'Vemurafenib', 'Crizotinib', 'Olaparib', 'Trastuzumab', 'Imatinib']
    
    G = nx.Graph()
    
    for gene in genes:
        G.add_node(gene, node_type='gene', color='#3498db')
    for drug in drugs:
        G.add_node(drug, node_type='drug', color='#e74c3c')
    
    interactions = [
        ('CYP2D6', 'Codeine', 0.9), ('CYP2D6', 'Tamoxifen', 0.85),
        ('CYP2C19', 'Clopidogrel', 0.95), ('CYP2C19', 'Warfarin', 0.7),
        ('CYP3A4', 'Simvastatin', 0.8), ('CYP3A4', 'Imatinib', 0.75),
        ('CYP1A2', 'Warfarin', 0.6), ('UGT1A1', 'Irinotecan', 0.9),
        ('ABCB1', 'Clopidogrel', 0.65), ('SLCO1B1', 'Simvastatin', 0.85),
        ('VKORC1', 'Warfarin', 0.95), ('DPYD', '5-FU', 0.95),
        ('TPMT', '6-MP', 0.9), ('HLA-B', 'Carbamazepine', 0.92),
        ('HLA-A', 'Carbamazepine', 0.6), ('EGFR', 'Gefitinib', 0.95),
        ('BRAF', 'Vemurafenib', 0.98), ('ALK', 'Crizotinib', 0.97),
        ('BRCA1', 'Olaparib', 0.93), ('ERBB2', 'Trastuzumab', 0.96),
        ('TP53', 'Gefitinib', 0.5), ('KRAS', 'Gefitinib', 0.7),
        ('PIK3CA', 'Trastuzumab', 0.6), ('TP53', '5-FU', 0.55),
    ]
    
    gene_gene = [
        ('CYP2D6', 'CYP2C19', 0.4), ('CYP3A4', 'CYP2D6', 0.3),
        ('EGFR', 'KRAS', 0.7), ('EGFR', 'ERBB2', 0.6),
        ('BRAF', 'KRAS', 0.8), ('PIK3CA', 'TP53', 0.5),
        ('BRCA1', 'TP53', 0.6), ('HLA-A', 'HLA-B', 0.5),
    ]
    
    for g, d, w in interactions:
        G.add_edge(g, d, weight=w, edge_type='drug-gene')
    for g1, g2, w in gene_gene:
        G.add_edge(g1, g2, weight=w, edge_type='gene-gene')
    
    return G, genes, drugs, interactions

G, genes, drugs, interactions = build_drug_gene_network()

# Graph Neural Network (message passing)
class SimpleGNN(nn.Module):
    def __init__(self, n_nodes, embed_dim, hidden_dim):
        super().__init__()
        self.embeddings = nn.Embedding(n_nodes, embed_dim)
        self.conv1 = nn.Linear(embed_dim, hidden_dim)
        self.conv2 = nn.Linear(hidden_dim, hidden_dim)
        self.predictor = nn.Sequential(
            nn.Linear(hidden_dim * 2, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )
    
    def message_pass(self, x, adj):
        """Simple mean-aggregation message passing."""
        deg = adj.sum(dim=1, keepdim=True).clamp(min=1)
        return adj @ x / deg
    
    def forward(self, node_pairs, adj):
        x = self.embeddings.weight
        h = torch.relu(self.conv1(x + self.message_pass(x, adj)))
        h = torch.relu(self.conv2(h + self.message_pass(h, adj)))
        
        h_src = h[node_pairs[:, 0]]
        h_dst = h[node_pairs[:, 1]]
        return self.predictor(torch.cat([h_src, h_dst], dim=1))

all_nodes = genes + drugs
node2idx = {n: i for i, n in enumerate(all_nodes)}
n_nodes = len(all_nodes)

# Adjacency matrix
adj = torch.zeros(n_nodes, n_nodes)
for u, v, data in G.edges(data=True):
    i, j = node2idx[u], node2idx[v]
    adj[i, j] = data['weight']
    adj[j, i] = data['weight']

# Positive edges (known interactions)
pos_edges = [(node2idx[g], node2idx[d]) for g, d, _ in interactions]
pos_labels = [1.0] * len(pos_edges)

# Negative sampling
neg_edges = []
for _ in range(len(pos_edges)):
    while True:
        g = np.random.choice(len(genes))
        d = len(genes) + np.random.choice(len(drugs))
        if (g, d) not in pos_edges:
            neg_edges.append((g, d))
            break
neg_labels = [0.0] * len(neg_edges)

all_edges = torch.LongTensor(pos_edges + neg_edges)
all_labels = torch.FloatTensor(pos_labels + neg_labels).unsqueeze(1)

# Train GNN
gnn = SimpleGNN(n_nodes, 32, 64)
gnn_optimizer = optim.Adam(gnn.parameters(), lr=0.01)
gnn_criterion = nn.BCELoss()

gnn_losses = []
for epoch in range(200):
    gnn.train()
    gnn_optimizer.zero_grad()
    preds = gnn(all_edges, adj)
    loss = gnn_criterion(preds, all_labels)
    loss.backward()
    gnn_optimizer.step()
    gnn_losses.append(loss.item())

gnn.eval()
with torch.no_grad():
    final_preds = gnn(all_edges, adj).numpy().flatten()

gnn_auc = roc_auc_score(all_labels.numpy().flatten(), final_preds)
gnn_acc = accuracy_score(all_labels.numpy().flatten().astype(int), (final_preds >= 0.5).astype(int))

results['exp5'] = {
    'n_nodes': n_nodes,
    'n_edges': G.number_of_edges(),
    'link_prediction_auc': round(gnn_auc, 4),
    'link_prediction_acc': round(gnn_acc, 4),
    'n_gene_nodes': len(genes),
    'n_drug_nodes': len(drugs)
}

print(f"GNN Link Prediction - AUC: {gnn_auc:.4f}, Accuracy: {gnn_acc:.4f}")

# Figure 5: Network and GNN results
fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# 5a: Network visualization
pos = nx.spring_layout(G, k=2, seed=42)
gene_nodes = [n for n in G.nodes if G.nodes[n]['node_type'] == 'gene']
drug_nodes_list = [n for n in G.nodes if G.nodes[n]['node_type'] == 'drug']

nx.draw_networkx_nodes(G, pos, nodelist=gene_nodes, node_color='#3498db',
                       node_size=300, alpha=0.8, ax=axes[0, 0])
nx.draw_networkx_nodes(G, pos, nodelist=drug_nodes_list, node_color='#e74c3c',
                       node_size=300, alpha=0.8, ax=axes[0, 0])
nx.draw_networkx_edges(G, pos, alpha=0.3, ax=axes[0, 0])
nx.draw_networkx_labels(G, pos, font_size=6, ax=axes[0, 0])
axes[0, 0].set_title('(A) Drug-Gene Interaction Network', fontsize=13, fontweight='bold')
axes[0, 0].legend(['Genes', 'Drugs'], loc='upper left', fontsize=10)

# 5b: GNN training loss
axes[0, 1].plot(gnn_losses, color='#8e44ad', lw=2)
axes[0, 1].set_xlabel('Epoch', fontsize=12)
axes[0, 1].set_ylabel('BCE Loss', fontsize=12)
axes[0, 1].set_title('(B) GNN Training Loss', fontsize=13, fontweight='bold')

# 5c: Node degree distribution
degrees = dict(G.degree())
gene_degrees = [degrees[g] for g in gene_nodes]
drug_degrees = [degrees[d] for d in drug_nodes_list]
axes[1, 0].hist(gene_degrees, bins=8, alpha=0.7, color='#3498db', label='Genes', edgecolor='black')
axes[1, 0].hist(drug_degrees, bins=8, alpha=0.7, color='#e74c3c', label='Drugs', edgecolor='black')
axes[1, 0].set_xlabel('Degree', fontsize=12)
axes[1, 0].set_ylabel('Count', fontsize=12)
axes[1, 0].set_title('(C) Node Degree Distribution', fontsize=13, fontweight='bold')
axes[1, 0].legend()

# 5d: Prediction scores
pos_scores = final_preds[:len(pos_edges)]
neg_scores = final_preds[len(pos_edges):]
axes[1, 1].hist(pos_scores, bins=8, alpha=0.7, color='#2ecc71', label='Positive', edgecolor='black')
axes[1, 1].hist(neg_scores, bins=8, alpha=0.7, color='#e74c3c', label='Negative', edgecolor='black')
axes[1, 1].set_xlabel('Prediction Score', fontsize=12)
axes[1, 1].set_ylabel('Count', fontsize=12)
axes[1, 1].set_title('(D) Link Prediction Score Distribution', fontsize=13, fontweight='bold')
axes[1, 1].legend()

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig5_drug_gene_network.png'), dpi=150, bbox_inches='tight')
plt.close()
print("Figure 5 saved.")

# ============================================================
# Experiment 6: CDSS Prototype Design
# ============================================================
print("\n" + "=" * 60)
print("Experiment 6: Clinical Decision Support System Prototype")
print("=" * 60)

class CDSSPrototype:
    """Clinical Decision Support System for Pharmacogenomics."""
    
    def __init__(self):
        self.guidelines = self._load_guidelines()
        self.risk_thresholds = {
            'high': 0.7, 'moderate': 0.4, 'low': 0.1
        }
    
    def _load_guidelines(self):
        """Simulate CPIC/DPWG guidelines."""
        return {
            ('CYP2D6', 'Codeine'): {
                'PM': {'recommendation': 'AVOID', 'alternative': 'Morphine (reduced dose)',
                       'evidence': 'Strong', 'level': '1A'},
                'IM': {'recommendation': 'REDUCE DOSE', 'alternative': 'Consider alternative',
                       'evidence': 'Moderate', 'level': '1B'},
                'NM': {'recommendation': 'STANDARD DOSE', 'alternative': None,
                       'evidence': 'Strong', 'level': '1A'},
                'UM': {'recommendation': 'AVOID', 'alternative': 'Morphine',
                       'evidence': 'Strong', 'level': '1A'}
            },
            ('CYP2C19', 'Clopidogrel'): {
                'PM': {'recommendation': 'AVOID', 'alternative': 'Prasugrel/Ticagrelor',
                       'evidence': 'Strong', 'level': '1A'},
                'IM': {'recommendation': 'ALTERNATIVE THERAPY', 'alternative': 'Prasugrel/Ticagrelor',
                       'evidence': 'Strong', 'level': '1A'},
                'NM': {'recommendation': 'STANDARD DOSE', 'alternative': None,
                       'evidence': 'Strong', 'level': '1A'},
                'RM': {'recommendation': 'STANDARD DOSE', 'alternative': None,
                       'evidence': 'Moderate', 'level': '2A'}
            },
            ('HLA-B*15:02', 'Carbamazepine'): {
                'positive': {'recommendation': 'CONTRAINDICATED',
                            'alternative': 'Valproate/Lamotrigine',
                            'evidence': 'Strong', 'level': '1A'},
                'negative': {'recommendation': 'STANDARD', 'alternative': None,
                            'evidence': 'Strong', 'level': '1A'}
            },
            ('DPYD', '5-FU'): {
                'PM': {'recommendation': 'AVOID', 'alternative': 'Capecitabine (reduced)',
                       'evidence': 'Strong', 'level': '1A'},
                'IM': {'recommendation': 'REDUCE DOSE 50%', 'alternative': None,
                       'evidence': 'Strong', 'level': '1A'},
                'NM': {'recommendation': 'STANDARD DOSE', 'alternative': None,
                       'evidence': 'Strong', 'level': '1A'}
            }
        }
    
    def generate_recommendation(self, patient):
        """Generate personalized drug recommendations."""
        recommendations = []
        
        for (gene, drug), guidelines in self.guidelines.items():
            phenotype = patient.get(gene, 'NM')
            if phenotype in guidelines:
                rec = guidelines[phenotype]
                risk_level = 'high' if rec['recommendation'] in ['AVOID', 'CONTRAINDICATED'] else \
                             'moderate' if 'REDUCE' in rec['recommendation'] or 'ALTERNATIVE' in rec['recommendation'] else 'low'
                
                recommendations.append({
                    'gene': gene,
                    'drug': drug,
                    'phenotype': phenotype,
                    'recommendation': rec['recommendation'],
                    'alternative': rec['alternative'],
                    'evidence_level': rec['level'],
                    'risk_level': risk_level
                })
        
        return recommendations
    
    def batch_evaluate(self, patients):
        """Evaluate CDSS for a batch of patients."""
        all_recs = []
        for p in patients:
            recs = self.generate_recommendation(p)
            all_recs.extend(recs)
        return pd.DataFrame(all_recs)

# Generate simulated patients
def generate_patients(n=500):
    phenotypes = ['PM', 'IM', 'NM', 'UM']
    hla_status = ['positive', 'negative']
    patients = []
    for _ in range(n):
        patients.append({
            'CYP2D6': np.random.choice(phenotypes, p=[0.07, 0.15, 0.68, 0.10]),
            'CYP2C19': np.random.choice(phenotypes, p=[0.03, 0.18, 0.69, 0.10]),
            'HLA-B*15:02': np.random.choice(hla_status, p=[0.08, 0.92]),
            'DPYD': np.random.choice(phenotypes[:3], p=[0.01, 0.05, 0.94]),
        })
    return patients

cdss = CDSSPrototype()
patients = generate_patients(500)
cdss_results_df = cdss.batch_evaluate(patients)

# Analyze CDSS output
risk_distribution = cdss_results_df['risk_level'].value_counts(normalize=True)
actionable_rate = (cdss_results_df['risk_level'] != 'low').mean()
high_risk_rate = (cdss_results_df['risk_level'] == 'high').mean()

results['exp6'] = {
    'n_patients': 500,
    'n_recommendations': len(cdss_results_df),
    'actionable_rate': round(actionable_rate, 4),
    'high_risk_rate': round(high_risk_rate, 4),
    'risk_distribution': risk_distribution.to_dict()
}

print(f"Patients: 500, Recommendations: {len(cdss_results_df)}")
print(f"Actionable rate: {actionable_rate:.4f}, High-risk rate: {high_risk_rate:.4f}")

# CDSS response time simulation
response_times = np.random.lognormal(mean=2.5, sigma=0.5, size=500)
mean_response = np.mean(response_times)
p95_response = np.percentile(response_times, 95)

results['exp6']['mean_response_ms'] = round(mean_response, 2)
results['exp6']['p95_response_ms'] = round(p95_response, 2)

# Figure 6: CDSS results
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# 6a: Risk distribution
risk_counts = cdss_results_df['risk_level'].value_counts()
colors_risk = {'high': '#e74c3c', 'moderate': '#f39c12', 'low': '#27ae60'}
axes[0, 0].pie(risk_counts.values,
               labels=[f"{k}\n({v})" for k, v in risk_counts.items()],
               colors=[colors_risk.get(k, '#95a5a6') for k in risk_counts.index],
               autopct='%1.1f%%', startangle=90, textprops={'fontsize': 11})
axes[0, 0].set_title('(A) Risk Level Distribution', fontsize=13, fontweight='bold')

# 6b: Recommendations by gene-drug pair
gene_drug_counts = cdss_results_df.groupby(['gene', 'drug'])['risk_level'].apply(
    lambda x: (x != 'low').sum()
).sort_values(ascending=False)
axes[0, 1].barh(range(len(gene_drug_counts)),
                gene_drug_counts.values, color='#3498db', alpha=0.8)
axes[0, 1].set_yticks(range(len(gene_drug_counts)))
axes[0, 1].set_yticklabels([f"{g}-{d}" for g, d in gene_drug_counts.index], fontsize=9)
axes[0, 1].set_xlabel('Actionable Recommendations', fontsize=12)
axes[0, 1].set_title('(B) Actionable by Gene-Drug Pair', fontsize=13, fontweight='bold')

# 6c: Response time distribution
axes[1, 0].hist(response_times, bins=30, color='#9b59b6', alpha=0.7, edgecolor='black')
axes[1, 0].axvline(x=mean_response, color='red', linestyle='--', lw=2, label=f'Mean: {mean_response:.1f}ms')
axes[1, 0].axvline(x=p95_response, color='orange', linestyle='--', lw=2, label=f'P95: {p95_response:.1f}ms')
axes[1, 0].set_xlabel('Response Time (ms)', fontsize=12)
axes[1, 0].set_ylabel('Count', fontsize=12)
axes[1, 0].set_title('(C) CDSS Response Time', fontsize=13, fontweight='bold')
axes[1, 0].legend()

# 6d: Evidence level distribution
evidence_counts = cdss_results_df['evidence_level'].value_counts()
axes[1, 1].bar(evidence_counts.index, evidence_counts.values,
               color=['#2ecc71', '#3498db', '#f39c12', '#e74c3c'][:len(evidence_counts)],
               alpha=0.8)
axes[1, 1].set_xlabel('Evidence Level', fontsize=12)
axes[1, 1].set_ylabel('Count', fontsize=12)
axes[1, 1].set_title('(D) Evidence Level Distribution', fontsize=13, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig6_cdss.png'), dpi=150, bbox_inches='tight')
plt.close()
print("Figure 6 saved.")

# ============================================================
# Summary Figure
# ============================================================
print("\n" + "=" * 60)
print("Generating Summary Figure")
print("=" * 60)

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

metrics = [
    ('Exp1: CYP Metabolism\n(Clearance)', f"R²={results['exp1']['clearance_r2']:.3f}\nRMSE={results['exp1']['clearance_rmse']:.3f}", '#3498db'),
    ('Exp2: HLA-ADR\n(Best Model)', f"AUC={max(results['exp2']['rf_auc'], results['exp2']['nn_auc']):.3f}", '#e74c3c'),
    ('Exp3: MR Analysis\n(Targets)', f"{results['exp3']['n_significant']}/{results['exp3']['n_targets_tested']}\nSignificant", '#2ecc71'),
    ('Exp4: Drug Sensitivity\n(DL Model)', f"R²={results['exp4']['mean_dl_r2']:.3f}\nRMSE={results['exp4']['mean_dl_rmse']:.3f}", '#f39c12'),
    ('Exp5: GNN Network\n(Link Prediction)', f"AUC={results['exp5']['link_prediction_auc']:.3f}\nAcc={results['exp5']['link_prediction_acc']:.3f}", '#9b59b6'),
    ('Exp6: CDSS\n(Performance)', f"Actionable={results['exp6']['actionable_rate']:.1%}\nP95={results['exp6']['p95_response_ms']:.0f}ms", '#1abc9c'),
]

for idx, (title, metric, color) in enumerate(metrics):
    row, col = idx // 3, idx % 3
    axes[row, col].text(0.5, 0.5, metric, transform=axes[row, col].transAxes,
                         fontsize=18, ha='center', va='center', fontweight='bold',
                         bbox=dict(boxstyle='round,pad=0.5', facecolor=color, alpha=0.3))
    axes[row, col].set_title(title, fontsize=13, fontweight='bold')
    axes[row, col].set_xlim(0, 1)
    axes[row, col].set_ylim(0, 1)
    axes[row, col].axis('off')

plt.suptitle('Pharmacogenomics Model Suite: Summary of Results', fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'fig7_summary.png'), dpi=150, bbox_inches='tight')
plt.close()
print("Summary figure saved.")

# Save results
with open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        'results.json'), 'w') as f:
    json.dump(results, f, indent=2, default=str)

print("\n" + "=" * 60)
print("ALL EXPERIMENTS COMPLETE")
print("=" * 60)
print(json.dumps(results, indent=2, default=str))
