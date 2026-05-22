"""
Module 1: CYP2D6/CYP2C19 Polymorphism and Drug Metabolism Rate Modeling
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
from sklearn.preprocessing import LabelEncoder
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ─────────────────────────────────────────────────────────────
# Star allele activity scores (CPIC guidelines)
# ─────────────────────────────────────────────────────────────
CYP2D6_ACTIVITY_SCORES = {
    '*1': 1.0, '*2': 1.0, '*3': 0.0, '*4': 0.0,
    '*5': 0.0, '*6': 0.0, '*10': 0.25, '*17': 0.5,
    '*29': 0.5, '*41': 0.5, '*2xN': 2.0, '*1xN': 2.0,
}
CYP2C19_ACTIVITY_SCORES = {
    '*1': 1.0, '*2': 0.0, '*3': 0.0, '*4': 0.0,
    '*17': 1.5, '*1/*17': 1.25, '*2/*17': 0.5,
}

def classify_metabolizer(activity_score: float, gene: str = 'CYP2D6') -> str:
    if gene == 'CYP2D6':
        if activity_score == 0:      return 'Poor Metabolizer (PM)'
        elif activity_score < 1.25:  return 'Intermediate Metabolizer (IM)'
        elif activity_score <= 2.25: return 'Normal Metabolizer (NM)'
        else:                        return 'Ultrarapid Metabolizer (UM)'
    else:  # CYP2C19
        if activity_score == 0:      return 'Poor Metabolizer (PM)'
        elif activity_score < 1.25:  return 'Intermediate Metabolizer (IM)'
        elif activity_score <= 1.5:  return 'Normal Metabolizer (NM)'
        else:                        return 'Rapid/Ultrarapid Metabolizer (RM/UM)'

# ─────────────────────────────────────────────────────────────
# Synthetic patient cohort (n=1000)
# ─────────────────────────────────────────────────────────────
allele_pairs_2d6 = [
    ('*1','*1'), ('*1','*4'), ('*4','*4'), ('*1','*2xN'),
    ('*1','*10'), ('*4','*10'), ('*1','*41'), ('*2','*41')
]
allele_pairs_2c19 = [
    ('*1','*1'), ('*1','*2'), ('*2','*2'), ('*1','*17'),
    ('*2','*17'), ('*17','*17'), ('*1','*3'), ('*2','*3')
]
probs_2d6  = [0.35, 0.20, 0.07, 0.05, 0.10, 0.08, 0.08, 0.07]
probs_2c19 = [0.38, 0.20, 0.03, 0.22, 0.10, 0.04, 0.02, 0.01]

n_patients = 1000
idx_2d6  = np.random.choice(len(allele_pairs_2d6),  n_patients, p=probs_2d6)
idx_2c19 = np.random.choice(len(allele_pairs_2c19), n_patients, p=probs_2c19)

records = []
for i in range(n_patients):
    a1_2d6, a2_2d6   = allele_pairs_2d6[idx_2d6[i]]
    a1_2c19, a2_2c19 = allele_pairs_2c19[idx_2c19[i]]

    as_2d6  = CYP2D6_ACTIVITY_SCORES.get(a1_2d6, 1.0) + CYP2D6_ACTIVITY_SCORES.get(a2_2d6, 1.0)
    as_2c19 = CYP2C19_ACTIVITY_SCORES.get(a1_2c19, 1.0) + CYP2C19_ACTIVITY_SCORES.get(a2_2c19, 1.0)

    pheno_2d6  = classify_metabolizer(as_2d6,  'CYP2D6')
    pheno_2c19 = classify_metabolizer(as_2c19, 'CYP2C19')

    # Codeine (CYP2D6 substrate): simulate plasma concentration
    codeine_auc_base = 200.0
    if 'Poor' in pheno_2d6:
        codeine_auc = codeine_auc_base * np.random.normal(0.3, 0.05)
        toxicity    = np.random.binomial(1, 0.05)
    elif 'Ultrarapid' in pheno_2d6:
        codeine_auc = codeine_auc_base * np.random.normal(2.8, 0.3)
        toxicity    = np.random.binomial(1, 0.45)
    elif 'Intermediate' in pheno_2d6:
        codeine_auc = codeine_auc_base * np.random.normal(0.7, 0.1)
        toxicity    = np.random.binomial(1, 0.08)
    else:
        codeine_auc = codeine_auc_base * np.random.normal(1.0, 0.15)
        toxicity    = np.random.binomial(1, 0.10)

    # Clopidogrel (CYP2C19): simulate active metabolite
    clopi_active_base = 100.0
    if 'Poor' in pheno_2c19:
        clopi_active = clopi_active_base * np.random.normal(0.25, 0.05)
        clopi_efficacy = 0
    elif 'Rapid' in pheno_2c19 or 'Ultrarapid' in pheno_2c19:
        clopi_active = clopi_active_base * np.random.normal(1.6, 0.2)
        clopi_efficacy = 1
    elif 'Intermediate' in pheno_2c19:
        clopi_active = clopi_active_base * np.random.normal(0.6, 0.1)
        clopi_efficacy = np.random.binomial(1, 0.5)
    else:
        clopi_active = clopi_active_base * np.random.normal(1.0, 0.15)
        clopi_efficacy = np.random.binomial(1, 0.75)

    records.append({
        'patient_id': f'P{i+1:04d}',
        'CYP2D6_allele1': a1_2d6, 'CYP2D6_allele2': a2_2d6,
        'CYP2D6_activity_score': as_2d6, 'CYP2D6_phenotype': pheno_2d6,
        'CYP2C19_allele1': a1_2c19, 'CYP2C19_allele2': a2_2c19,
        'CYP2C19_activity_score': as_2c19, 'CYP2C19_phenotype': pheno_2c19,
        'codeine_AUC': max(0, codeine_auc),
        'codeine_toxicity': toxicity,
        'clopidogrel_active_metabolite': max(0, clopi_active),
        'clopidogrel_efficacy': clopi_efficacy,
        'age': np.random.randint(18, 80),
        'weight_kg': np.random.normal(70, 15),
        'sex': np.random.choice(['M', 'F']),
    })

df = pd.DataFrame(records)
df.to_csv('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/data/cyp_patient_cohort.csv', index=False)

# ─────────────────────────────────────────────────────────────
# Figure 1: Phenotype distribution
# ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
order_2d6  = ['Poor Metabolizer (PM)', 'Intermediate Metabolizer (IM)',
               'Normal Metabolizer (NM)', 'Ultrarapid Metabolizer (UM)']
order_2c19 = ['Poor Metabolizer (PM)', 'Intermediate Metabolizer (IM)',
               'Normal Metabolizer (NM)', 'Rapid/Ultrarapid Metabolizer (RM/UM)']
colors = ['#d73027', '#fc8d59', '#91bfdb', '#4575b4']

for ax, gene, order, col_name in zip(
    axes,
    ['CYP2D6', 'CYP2C19'],
    [order_2d6, order_2c19],
    ['CYP2D6_phenotype', 'CYP2C19_phenotype']
):
    counts = df[col_name].value_counts()
    vals = [counts.get(o, 0) for o in order]
    bars = ax.bar(range(len(order)), vals, color=colors)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(['PM','IM','NM','UM'], fontsize=11)
    ax.set_title(f'{gene} Metabolizer Phenotype Distribution (n={n_patients})', fontsize=12)
    ax.set_ylabel('Patient Count')
    ax.set_xlabel('Metabolizer Phenotype')
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                f'{val}\n({val/n_patients*100:.1f}%)', ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.savefig('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/figures/fig1_cyp_phenotype_distribution.png', dpi=150, bbox_inches='tight')
plt.close()

# ─────────────────────────────────────────────────────────────
# Figure 2: AUC by phenotype (Codeine)
# ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
pheno_order = ['Poor Metabolizer (PM)', 'Intermediate Metabolizer (IM)',
               'Normal Metabolizer (NM)', 'Ultrarapid Metabolizer (UM)']
pheno_labels = ['PM', 'IM', 'NM', 'UM']
data_by_pheno = [df[df['CYP2D6_phenotype'] == p]['codeine_AUC'].values for p in pheno_order]
bp = ax.boxplot(data_by_pheno, labels=pheno_labels, patch_artist=True)
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_title('Codeine AUC by CYP2D6 Phenotype', fontsize=12)
ax.set_xlabel('Metabolizer Phenotype')
ax.set_ylabel('AUC (ng·h/mL)')
ax.axhline(200, color='gray', linestyle='--', alpha=0.5, label='Reference AUC')
ax.legend()

ax = axes[1]
pheno_order_c19 = ['Poor Metabolizer (PM)', 'Intermediate Metabolizer (IM)',
                    'Normal Metabolizer (NM)', 'Rapid/Ultrarapid Metabolizer (RM/UM)']
pheno_labels_c19 = ['PM', 'IM', 'NM', 'RM/UM']
data_clopi = [df[df['CYP2C19_phenotype'] == p]['clopidogrel_active_metabolite'].values
               for p in pheno_order_c19]
bp2 = ax.boxplot(data_clopi, labels=pheno_labels_c19, patch_artist=True)
for patch, color in zip(bp2['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_title('Clopidogrel Active Metabolite by CYP2C19 Phenotype', fontsize=12)
ax.set_xlabel('Metabolizer Phenotype')
ax.set_ylabel('Active Metabolite Concentration (ng/mL)')

plt.tight_layout()
plt.savefig('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/figures/fig2_drug_auc_by_phenotype.png', dpi=150, bbox_inches='tight')
plt.close()

# ─────────────────────────────────────────────────────────────
# ML model: predict codeine toxicity from CYP2D6 genotype
# ─────────────────────────────────────────────────────────────
le_allele = LabelEncoder()
df['allele1_enc'] = le_allele.fit_transform(df['CYP2D6_allele1'])
le2 = LabelEncoder()
df['allele2_enc'] = le2.fit_transform(df['CYP2D6_allele2'])
le_sex = LabelEncoder()
df['sex_enc'] = le_sex.fit_transform(df['sex'])

features = ['CYP2D6_activity_score', 'allele1_enc', 'allele2_enc',
            'age', 'weight_kg', 'sex_enc']
X = df[features].values
y = df['codeine_toxicity'].values

clf = RandomForestClassifier(n_estimators=100, random_state=42)
cv   = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
aucs = cross_val_score(clf, X, y, cv=cv, scoring='roc_auc')
accs = cross_val_score(clf, X, y, cv=cv, scoring='accuracy')

clf.fit(X, y)
importances = pd.Series(clf.feature_importances_, index=features).sort_values(ascending=False)

results_cyp = {
    'module': 'CYP Metabolism Modeling',
    'n_patients': n_patients,
    'CYP2D6_phenotype_distribution': df['CYP2D6_phenotype'].value_counts().to_dict(),
    'CYP2C19_phenotype_distribution': df['CYP2C19_phenotype'].value_counts().to_dict(),
    'codeine_toxicity_rate_by_phenotype': df.groupby('CYP2D6_phenotype')['codeine_toxicity'].mean().to_dict(),
    'clopidogrel_efficacy_rate_by_phenotype': df.groupby('CYP2C19_phenotype')['clopidogrel_efficacy'].mean().to_dict(),
    'rf_toxicity_prediction': {
        'cv_roc_auc_mean': float(aucs.mean()),
        'cv_roc_auc_std': float(aucs.std()),
        'cv_accuracy_mean': float(accs.mean()),
        'feature_importances': importances.to_dict(),
    }
}
with open('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/results/cyp_metabolism_results.json', 'w') as f:
    json.dump(results_cyp, f, indent=2, ensure_ascii=False)

print("[CYP Module] Done")
print(f"  RF Toxicity AUC: {aucs.mean():.3f} ± {aucs.std():.3f}")
print(f"  RF Accuracy:     {accs.mean():.3f} ± {accs.std():.3f}")
print(f"  CYP2D6 phenotypes: {df['CYP2D6_phenotype'].value_counts().to_dict()}")
