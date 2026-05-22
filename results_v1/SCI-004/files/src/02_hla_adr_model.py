"""
Module 2: HLA Genotype and Drug Adverse Reactions (Carbamazepine/HLA-B*1502)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold, train_test_split
from sklearn.metrics import (roc_auc_score, roc_curve, precision_recall_curve,
                              classification_report, confusion_matrix)
from sklearn.preprocessing import StandardScaler
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ─────────────────────────────────────────────────────────────
# Known HLA-drug associations (literature-based)
# ─────────────────────────────────────────────────────────────
HLA_DRUG_ASSOCIATIONS = {
    'HLA-B*15:02': {
        'drug': 'Carbamazepine',
        'reaction': 'SJS/TEN',
        'OR': 80.1,
        'sensitivity': 0.98,
        'specificity': 0.94,
        'ethnicity_risk': {'Asian': 0.06, 'European': 0.001},
    },
    'HLA-B*58:01': {
        'drug': 'Allopurinol',
        'reaction': 'SJS/TEN',
        'OR': 580.0,
        'sensitivity': 0.97,
        'specificity': 0.92,
        'ethnicity_risk': {'Asian': 0.07, 'European': 0.002},
    },
    'HLA-B*57:01': {
        'drug': 'Abacavir',
        'reaction': 'Hypersensitivity',
        'OR': 117.0,
        'sensitivity': 0.48,
        'specificity': 0.99,
        'ethnicity_risk': {'European': 0.055, 'African': 0.003},
    },
    'HLA-A*31:01': {
        'drug': 'Carbamazepine',
        'reaction': 'DRESS/MPE',
        'OR': 9.0,
        'sensitivity': 0.26,
        'specificity': 0.94,
        'ethnicity_risk': {'Asian': 0.05, 'European': 0.025},
    },
}

# ─────────────────────────────────────────────────────────────
# Synthetic cohort (n=2000): focus on CBZ/HLA-B*1502
# ─────────────────────────────────────────────────────────────
n = 2000
ethnicities = np.random.choice(['Asian', 'European', 'African', 'Other'],
                                n, p=[0.30, 0.45, 0.15, 0.10])

# HLA-B*15:02 prevalence by ethnicity
hla_b1502 = np.array([
    np.random.binomial(1, 0.06 if e == 'Asian' else
                          0.001 if e == 'European' else
                          0.005 if e == 'African' else 0.002)
    for e in ethnicities
])
hla_a3101 = np.array([
    np.random.binomial(1, 0.05 if e == 'Asian' else
                          0.025 if e == 'European' else 0.01)
    for e in ethnicities
])
hla_b5701 = np.array([
    np.random.binomial(1, 0.055 if e == 'European' else
                          0.003 if e == 'African' else 0.01)
    for e in ethnicities
])

# Carbamazepine SJS/TEN risk
# P(SJS/TEN | HLA-B*15:02) ~ 4-5%, P(SJS/TEN | no HLA-B*15:02) ~ 0.01%
p_sjs = np.where(hla_b1502 == 1, 0.05, 0.0005)
p_sjs = np.where((hla_b1502 == 1) & (hla_a3101 == 1), 0.08, p_sjs)

cbz_sjs   = np.array([np.random.binomial(1, p) for p in p_sjs])
cbz_dress = np.array([np.random.binomial(1, 0.03 if a == 1 else 0.001)
                       for a in hla_a3101])

# Additional covariates
dose_mg   = np.random.normal(600, 150, n).clip(200, 1200)
age       = np.random.randint(5, 80, n)
sex       = np.random.choice([0, 1], n)
prior_adr = np.random.binomial(1, 0.05, n)
renal_fn  = np.random.normal(90, 20, n).clip(10, 130)

df_hla = pd.DataFrame({
    'patient_id': [f'P{i+1:04d}' for i in range(n)],
    'ethnicity': ethnicities,
    'HLA_B1502': hla_b1502,
    'HLA_A3101': hla_a3101,
    'HLA_B5701': hla_b5701,
    'cbz_dose_mg': dose_mg,
    'age': age, 'sex': sex,
    'prior_adr': prior_adr,
    'renal_function_eGFR': renal_fn,
    'cbz_SJS_TEN': cbz_sjs,
    'cbz_DRESS': cbz_dress,
})
df_hla['any_severe_adr'] = ((df_hla['cbz_SJS_TEN'] | df_hla['cbz_DRESS'])).astype(int)
df_hla.to_csv('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/data/hla_drug_cohort.csv', index=False)

# ─────────────────────────────────────────────────────────────
# Model: predict SJS/TEN from HLA + covariates
# ─────────────────────────────────────────────────────────────
features_hla = ['HLA_B1502', 'HLA_A3101', 'HLA_B5701',
                 'cbz_dose_mg', 'age', 'sex', 'prior_adr', 'renal_function_eGFR']
X_hla = df_hla[features_hla].values
y_hla = df_hla['cbz_SJS_TEN'].values

scaler    = StandardScaler()
X_scaled  = scaler.fit_transform(X_hla)
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y_hla,
                                                      test_size=0.2, stratify=y_hla,
                                                      random_state=42)

models = {
    'Logistic Regression': LogisticRegression(class_weight='balanced', max_iter=1000),
    'Gradient Boosting':   GradientBoostingClassifier(n_estimators=100, random_state=42),
    'Random Forest':       RandomForestClassifier(n_estimators=100, class_weight='balanced', random_state=42),
}

model_results = {}
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

colors_roc = ['#1b7837', '#762a83', '#c51b7d']
for (name, model), color in zip(models.items(), colors_roc):
    model.fit(X_train, y_train)
    y_prob = model.predict_proba(X_test)[:, 1]
    auc    = roc_auc_score(y_test, y_prob)
    fpr, tpr, _ = roc_curve(y_test, y_prob)
    axes[0].plot(fpr, tpr, color=color, lw=2, label=f'{name} (AUC={auc:.3f})')
    model_results[name] = {'auc': auc}

axes[0].plot([0,1],[0,1],'k--', alpha=0.3)
axes[0].set_xlabel('False Positive Rate')
axes[0].set_ylabel('True Positive Rate')
axes[0].set_title('ROC Curves: CBZ SJS/TEN Prediction\n(HLA-B*15:02 + Covariates)')
axes[0].legend(loc='lower right')
axes[0].set_xlim([0, 1]); axes[0].set_ylim([0, 1.02])

# Sensitivity/Specificity by HLA carrier status (clinical test performance)
hla_pos  = df_hla[df_hla['HLA_B1502'] == 1]['cbz_SJS_TEN']
hla_neg  = df_hla[df_hla['HLA_B1502'] == 0]['cbz_SJS_TEN']
tp = hla_pos.sum()
fn = len(hla_pos) - tp
tn_rate = 1 - hla_neg.mean()
fp_rate = hla_neg.mean()

bar_data = {
    'HLA-B*15:02\nPositive': hla_pos.mean() * 100,
    'HLA-B*15:02\nNegative': hla_neg.mean() * 100,
}
eth_risk = df_hla.groupby('ethnicity')['cbz_SJS_TEN'].mean() * 100
for eth, risk in eth_risk.items():
    bar_data[f'{eth}'] = risk

bar_labels = list(bar_data.keys())
bar_vals   = list(bar_data.values())
bar_colors = ['#d73027' if 'Positive' in l else
               '#1a9850' if 'Negative' in l else '#4575b4' for l in bar_labels]
axes[1].bar(bar_labels, bar_vals, color=bar_colors, alpha=0.8)
axes[1].set_title('SJS/TEN Incidence Rate by HLA Status & Ethnicity (%)')
axes[1].set_ylabel('Incidence Rate (%)')
axes[1].set_xlabel('Group')
plt.setp(axes[1].get_xticklabels(), rotation=30, ha='right')

plt.tight_layout()
plt.savefig('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/figures/fig3_hla_drug_reaction.png',
            dpi=150, bbox_inches='tight')
plt.close()

# ─────────────────────────────────────────────────────────────
# Figure 4: HLA allele frequency by ethnicity heatmap
# ─────────────────────────────────────────────────────────────
hla_freq = df_hla.groupby('ethnicity')[['HLA_B1502','HLA_A3101','HLA_B5701']].mean() * 100
hla_freq.columns = ['HLA-B*15:02', 'HLA-A*31:01', 'HLA-B*57:01']

fig, ax = plt.subplots(figsize=(8, 4))
sns.heatmap(hla_freq, annot=True, fmt='.2f', cmap='YlOrRd',
            linewidths=0.5, ax=ax,
            cbar_kws={'label': 'Allele Frequency (%)'})
ax.set_title('HLA Allele Frequency by Ethnicity (%)')
ax.set_ylabel('Ethnicity'); ax.set_xlabel('HLA Allele')
plt.tight_layout()
plt.savefig('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/figures/fig4_hla_allele_frequency.png',
            dpi=150, bbox_inches='tight')
plt.close()

# Compute clinical metrics
sensitivity_hlab1502 = tp / (tp + fn) if (tp + fn) > 0 else 0
tn  = df_hla[(df_hla['HLA_B1502']==0) & (df_hla['cbz_SJS_TEN']==0)].shape[0]
fp2 = df_hla[(df_hla['HLA_B1502']==1) & (df_hla['cbz_SJS_TEN']==0)].shape[0]
specificity_hlab1502 = tn / (tn + fp2) if (tn + fp2) > 0 else 0
nnt = 1 / (hla_pos.mean() - hla_neg.mean()) if (hla_pos.mean() - hla_neg.mean()) > 0 else float('inf')

results_hla = {
    'module': 'HLA Drug Adverse Reaction Prediction',
    'n_patients': n,
    'cbz_SJS_TEN_prevalence': float(y_hla.mean()),
    'HLA_B1502_prevalence': float(hla_b1502.mean()),
    'clinical_test_performance': {
        'sensitivity': float(sensitivity_hlab1502),
        'specificity': float(specificity_hlab1502),
        'SJS_rate_HLA_positive': float(hla_pos.mean()),
        'SJS_rate_HLA_negative': float(hla_neg.mean()),
        'NNS_to_prevent_1_SJS': float(nnt),
    },
    'ml_model_aucs': model_results,
    'hla_drug_reference': HLA_DRUG_ASSOCIATIONS,
}
with open('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/results/hla_adr_results.json', 'w') as f:
    json.dump(results_hla, f, indent=2, ensure_ascii=False)

print("[HLA Module] Done")
print(f"  SJS/TEN prevalence: {y_hla.mean()*100:.2f}%")
print(f"  HLA-B*15:02 sensitivity: {sensitivity_hlab1502:.3f}, specificity: {specificity_hlab1502:.3f}")
for name, res in model_results.items():
    print(f"  {name} AUC: {res['auc']:.3f}")
