"""
Ethical AI Evaluation Framework (EthicAI-Bench)
Integrates fairness, explainability, privacy, robustness, and environmental metrics.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.datasets import fetch_openml
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score, confusion_matrix
import shap
import json
import time
import os
import warnings
warnings.filterwarnings('ignore')

from fairlearn.metrics import (
    demographic_parity_difference,
    equalized_odds_difference,
    MetricFrame,
)

RESULTS = {}
FIGURES_DIR = 'figures'
os.makedirs(FIGURES_DIR, exist_ok=True)

np.random.seed(42)

# ============================================================
# 1. Data Preparation — Synthetic Medical Diagnosis Dataset
# ============================================================
def create_medical_dataset(n=3000):
    """Create synthetic medical diagnosis dataset with demographic attributes."""
    np.random.seed(42)
    age = np.random.normal(55, 15, n).clip(18, 90)
    gender = np.random.choice([0, 1], n, p=[0.48, 0.52])  # 0=Female, 1=Male
    ethnicity = np.random.choice([0, 1, 2], n, p=[0.6, 0.25, 0.15])
    bmi = np.random.normal(27, 5, n).clip(15, 50)
    blood_pressure = np.random.normal(130, 20, n).clip(80, 200)
    cholesterol = np.random.normal(200, 40, n).clip(100, 350)
    glucose = np.random.normal(100, 30, n).clip(50, 300)
    heart_rate = np.random.normal(75, 12, n).clip(40, 150)
    smoking = np.random.choice([0, 1], n, p=[0.7, 0.3])
    exercise = np.random.choice([0, 1, 2], n, p=[0.3, 0.4, 0.3])
    family_history = np.random.choice([0, 1], n, p=[0.65, 0.35])
    
    # Disease risk with intentional bias on gender/ethnicity
    logit = (
        0.03 * (age - 50)
        + 0.5 * gender  # intentional bias
        + 0.3 * (ethnicity == 2).astype(float)  # intentional bias
        + 0.04 * (bmi - 25)
        + 0.02 * (blood_pressure - 120)
        + 0.01 * (cholesterol - 180)
        + 0.015 * (glucose - 90)
        + 0.5 * smoking
        - 0.3 * (exercise == 2).astype(float)
        + 0.4 * family_history
        + np.random.normal(0, 0.5, n)
    )
    prob = 1 / (1 + np.exp(-logit + 1.5))
    diagnosis = (np.random.random(n) < prob).astype(int)
    
    df = pd.DataFrame({
        'age': age, 'gender': gender, 'ethnicity': ethnicity,
        'bmi': bmi, 'blood_pressure': blood_pressure,
        'cholesterol': cholesterol, 'glucose': glucose,
        'heart_rate': heart_rate, 'smoking': smoking,
        'exercise': exercise, 'family_history': family_history,
        'diagnosis': diagnosis
    })
    return df

print("=" * 60)
print("EthicAI-Bench: Ethical AI Evaluation Framework")
print("=" * 60)

df = create_medical_dataset()
feature_cols = [c for c in df.columns if c != 'diagnosis']
X = df[feature_cols].values
y = df['diagnosis'].values
sensitive_gender = df['gender'].values
sensitive_ethnicity = df['ethnicity'].values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42, stratify=y)
gender_train = sensitive_gender[df.index.isin(pd.Series(range(len(df))).sample(frac=0.7, random_state=42).values)]
gender_test = sensitive_gender[~df.index.isin(pd.Series(range(len(df))).sample(frac=0.7, random_state=42).values)]

# Proper split for sensitive attrs
idx = np.arange(len(df))
train_idx, test_idx = train_test_split(idx, test_size=0.3, random_state=42, stratify=y)
gender_train = sensitive_gender[train_idx]
gender_test = sensitive_gender[test_idx]
ethnicity_train = sensitive_ethnicity[train_idx]
ethnicity_test = sensitive_ethnicity[test_idx]

scaler = StandardScaler()
X_train_s = scaler.fit_transform(X_train)
X_test_s = scaler.transform(X_test)

# ============================================================
# 2. Train Models
# ============================================================
models = {
    'LogisticRegression': LogisticRegression(max_iter=1000, random_state=42),
    'RandomForest': RandomForestClassifier(n_estimators=100, random_state=42),
    'GradientBoosting': GradientBoostingClassifier(n_estimators=100, random_state=42),
}

trained = {}
predictions = {}
probabilities = {}

for name, model in models.items():
    print(f"\nTraining {name}...")
    t0 = time.time()
    model.fit(X_train_s, y_train)
    train_time = time.time() - t0
    y_pred = model.predict(X_test_s)
    y_prob = model.predict_proba(X_test_s)[:, 1]
    acc = accuracy_score(y_test, y_pred)
    auc = roc_auc_score(y_test, y_prob)
    f1 = f1_score(y_test, y_pred)
    trained[name] = model
    predictions[name] = y_pred
    probabilities[name] = y_prob
    print(f"  Accuracy={acc:.4f}, AUC={auc:.4f}, F1={f1:.4f}, Time={train_time:.2f}s")
    RESULTS.setdefault('performance', {})[name] = {
        'accuracy': round(acc, 4), 'auc': round(auc, 4),
        'f1': round(f1, 4), 'train_time': round(train_time, 3)
    }

# ============================================================
# 3. Fairness Evaluation (Statistical Parity, EO, Calibration)
# ============================================================
print("\n" + "=" * 60)
print("Module 1: Fairness Evaluation")
print("=" * 60)

fairness_results = {}

for name in models:
    y_pred = predictions[name]
    y_prob = probabilities[name]
    
    # Statistical Parity Difference (gender)
    spd_gender = demographic_parity_difference(y_test, y_pred, sensitive_features=gender_test)
    
    # Equalized Odds Difference (gender)
    eod_gender = equalized_odds_difference(y_test, y_pred, sensitive_features=gender_test)
    
    # Statistical Parity (ethnicity)
    spd_ethnicity = demographic_parity_difference(y_test, y_pred, sensitive_features=ethnicity_test)
    eod_ethnicity = equalized_odds_difference(y_test, y_pred, sensitive_features=ethnicity_test)
    
    # Calibration difference: |E[Y|S=a, score∈bin] - E[Y|S=b, score∈bin]| averaged
    def calibration_diff(y_true, y_prob, sensitive):
        groups = np.unique(sensitive)
        bins = np.linspace(0, 1, 11)
        cal_diffs = []
        for i in range(len(bins)-1):
            mask = (y_prob >= bins[i]) & (y_prob < bins[i+1])
            group_means = []
            for g in groups:
                gm = mask & (sensitive == g)
                if gm.sum() > 5:
                    group_means.append(y_true[gm].mean())
            if len(group_means) >= 2:
                cal_diffs.append(max(group_means) - min(group_means))
        return np.mean(cal_diffs) if cal_diffs else 0.0
    
    cal_gender = calibration_diff(y_test, y_prob, gender_test)
    cal_ethnicity = calibration_diff(y_test, y_prob, ethnicity_test)
    
    # Integrated Fairness Score (IFS) = 1 - (|SPD| + |EOD| + CalDiff) / 3
    ifs_gender = 1 - (abs(spd_gender) + abs(eod_gender) + cal_gender) / 3
    ifs_ethnicity = 1 - (abs(spd_ethnicity) + abs(eod_ethnicity) + cal_ethnicity) / 3
    
    fairness_results[name] = {
        'SPD_gender': round(spd_gender, 4), 'EOD_gender': round(eod_gender, 4),
        'Cal_gender': round(cal_gender, 4), 'IFS_gender': round(ifs_gender, 4),
        'SPD_ethnicity': round(spd_ethnicity, 4), 'EOD_ethnicity': round(eod_ethnicity, 4),
        'Cal_ethnicity': round(cal_ethnicity, 4), 'IFS_ethnicity': round(ifs_ethnicity, 4),
    }
    print(f"\n{name}:")
    print(f"  Gender  — SPD={spd_gender:.4f}, EOD={eod_gender:.4f}, Cal={cal_gender:.4f}, IFS={ifs_gender:.4f}")
    print(f"  Ethnicity — SPD={spd_ethnicity:.4f}, EOD={eod_ethnicity:.4f}, Cal={cal_ethnicity:.4f}, IFS={ifs_ethnicity:.4f}")

RESULTS['fairness'] = fairness_results

# ============================================================
# 4. Explainability — SHAP Consistency & Stability
# ============================================================
print("\n" + "=" * 60)
print("Module 2: Explainability Evaluation (SHAP)")
print("=" * 60)

explainability_results = {}

for name in models:
    model = trained[name]
    print(f"\nComputing SHAP for {name}...")
    
    # Use a subsample for speed
    X_sample = X_test_s[:200]
    
    if name == 'LogisticRegression':
        explainer = shap.LinearExplainer(model, X_train_s[:500])
    else:
        explainer = shap.TreeExplainer(model)
    
    shap_values = explainer.shap_values(X_sample)
    if isinstance(shap_values, list):
        shap_values = shap_values[1]
    
    # SHAP consistency: run multiple times (bootstrap) and measure rank correlation
    n_bootstrap = 10
    rank_correlations = []
    for _ in range(n_bootstrap):
        boot_idx = np.random.choice(len(X_sample), len(X_sample), replace=True)
        sv_boot = explainer.shap_values(X_sample[boot_idx])
        if isinstance(sv_boot, list):
            sv_boot = sv_boot[1]
        
        mean_orig = np.abs(shap_values).mean(axis=0)
        mean_boot = np.abs(sv_boot).mean(axis=0)
        
        from scipy.stats import spearmanr
        corr, _ = spearmanr(mean_orig, mean_boot)
        rank_correlations.append(corr)
    
    shap_consistency = np.mean(rank_correlations)
    shap_consistency_std = np.std(rank_correlations)
    
    # SHAP stability: variance of SHAP values across samples
    shap_stability = 1 / (1 + np.mean(np.var(shap_values, axis=0)))
    
    # Top-k feature agreement
    top_k = 5
    mean_abs = np.abs(shap_values).mean(axis=0)
    if mean_abs.ndim > 1:
        mean_abs = mean_abs.mean(axis=-1)
    top_features = list(np.argsort(-mean_abs)[:top_k])
    
    agreements = []
    for _ in range(n_bootstrap):
        boot_idx = np.random.choice(len(X_sample), len(X_sample), replace=True)
        sv_b = explainer.shap_values(X_sample[boot_idx])
        if isinstance(sv_b, list):
            sv_b = sv_b[1]
        m_b = np.abs(sv_b).mean(axis=0)
        if m_b.ndim > 1:
            m_b = m_b.mean(axis=-1)
        top_b = list(np.argsort(-m_b)[:top_k])
        agreements.append(len(set(top_features) & set(top_b)) / top_k)
    
    top_k_agreement = np.mean(agreements)
    
    explainability_results[name] = {
        'shap_consistency': round(shap_consistency, 4),
        'shap_consistency_std': round(shap_consistency_std, 4),
        'shap_stability': round(shap_stability, 4),
        'top_k_agreement': round(top_k_agreement, 4),
    }
    print(f"  Consistency (Spearman)={shap_consistency:.4f}±{shap_consistency_std:.4f}")
    print(f"  Stability={shap_stability:.4f}, Top-{top_k} Agreement={top_k_agreement:.4f}")

RESULTS['explainability'] = explainability_results

# Save SHAP summary plot for best model
best_model_name = 'GradientBoosting'
model = trained[best_model_name]
explainer = shap.TreeExplainer(model)
shap_vals = explainer.shap_values(X_test_s[:200])
if isinstance(shap_vals, list):
    shap_vals = shap_vals[1]

plt.figure(figsize=(10, 6))
shap.summary_plot(shap_vals, X_test_s[:200], feature_names=feature_cols, show=False)
plt.title(f'SHAP Feature Importance — {best_model_name}')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/shap_summary.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"\nSaved: {FIGURES_DIR}/shap_summary.png")

# ============================================================
# 5. Privacy Risk — Membership Inference Attack Resistance
# ============================================================
print("\n" + "=" * 60)
print("Module 3: Privacy Risk Evaluation (MIA)")
print("=" * 60)

privacy_results = {}

for name in models:
    model = trained[name]
    
    # Membership inference: compare train vs test confidence distributions
    train_probs = model.predict_proba(X_train_s)
    test_probs = model.predict_proba(X_test_s)
    
    train_conf = np.max(train_probs, axis=1)
    test_conf = np.max(test_probs, axis=1)
    
    # Simple threshold-based MIA
    thresholds = np.linspace(0.5, 1.0, 50)
    best_mia_acc = 0.5
    for t in thresholds:
        train_in = (train_conf >= t).mean()
        test_in = (test_conf >= t).mean()
        mia_acc = 0.5 * (train_in + (1 - test_in))
        best_mia_acc = max(best_mia_acc, mia_acc)
    
    # Privacy risk score: higher MIA accuracy = higher risk
    # Score = 1 - 2*(MIA_acc - 0.5) (0=high risk, 1=no risk)
    privacy_score = max(0, 1 - 2 * (best_mia_acc - 0.5))
    
    # Overfitting gap as privacy proxy
    train_acc = accuracy_score(y_train, model.predict(X_train_s))
    test_acc = accuracy_score(y_test, predictions[name])
    overfit_gap = train_acc - test_acc
    
    # Per-sample vulnerability
    train_loss = -np.log(np.clip(train_probs[np.arange(len(y_train)), y_train], 1e-10, 1))
    test_loss = -np.log(np.clip(test_probs[np.arange(len(y_test)), y_test], 1e-10, 1))
    vuln_ratio = (train_loss.mean()) / (test_loss.mean() + 1e-10)
    
    privacy_results[name] = {
        'mia_accuracy': round(best_mia_acc, 4),
        'privacy_score': round(privacy_score, 4),
        'overfit_gap': round(overfit_gap, 4),
        'loss_ratio': round(vuln_ratio, 4),
    }
    print(f"\n{name}:")
    print(f"  MIA Accuracy={best_mia_acc:.4f}, Privacy Score={privacy_score:.4f}")
    print(f"  Overfit Gap={overfit_gap:.4f}, Loss Ratio={vuln_ratio:.4f}")

RESULTS['privacy'] = privacy_results

# ============================================================
# 6. Robustness — Adversarial Perturbation & Distribution Shift
# ============================================================
print("\n" + "=" * 60)
print("Module 4: Robustness Evaluation")
print("=" * 60)

robustness_results = {}

for name in models:
    model = trained[name]
    y_pred_orig = predictions[name]
    base_acc = accuracy_score(y_test, y_pred_orig)
    
    # Adversarial perturbation (FGSM-like for feature space)
    epsilons = [0.01, 0.05, 0.1, 0.2, 0.5]
    adv_accs = []
    for eps in epsilons:
        X_perturbed = X_test_s + np.random.uniform(-eps, eps, X_test_s.shape)
        y_pred_adv = model.predict(X_perturbed)
        adv_acc = accuracy_score(y_test, y_pred_adv)
        adv_accs.append(adv_acc)
    
    # Robustness score: AUC of accuracy vs epsilon
    robustness_adv = np.trapz(adv_accs, epsilons) / (epsilons[-1] - epsilons[0])
    
    # Distribution shift: simulate covariate shift
    shift_results = {}
    for shift_mag in [0.5, 1.0, 2.0]:
        X_shifted = X_test_s + np.random.normal(0, shift_mag, X_test_s.shape) * 0.1
        y_pred_shift = model.predict(X_shifted)
        shift_acc = accuracy_score(y_test, y_pred_shift)
        shift_results[f'shift_{shift_mag}'] = round(shift_acc, 4)
    
    # Prediction consistency under noise
    n_noisy = 20
    consistent = np.zeros(len(X_test_s))
    for _ in range(n_noisy):
        X_noisy = X_test_s + np.random.normal(0, 0.05, X_test_s.shape)
        y_noisy = model.predict(X_noisy)
        consistent += (y_noisy == y_pred_orig)
    consistency_score = (consistent / n_noisy).mean()
    
    robustness_results[name] = {
        'base_accuracy': round(base_acc, 4),
        'adversarial_robustness': round(robustness_adv, 4),
        'consistency_score': round(consistency_score, 4),
        **shift_results,
        'adversarial_accs': {str(e): round(a, 4) for e, a in zip(epsilons, adv_accs)},
    }
    print(f"\n{name}:")
    print(f"  Base Acc={base_acc:.4f}, Adv Robustness={robustness_adv:.4f}")
    print(f"  Consistency={consistency_score:.4f}")
    for k, v in shift_results.items():
        print(f"  {k}: {v}")

RESULTS['robustness'] = robustness_results

# ============================================================
# 7. Environmental Impact — CO2 Estimation
# ============================================================
print("\n" + "=" * 60)
print("Module 5: Environmental Impact (CO2)")
print("=" * 60)

environmental_results = {}
TDP_WATTS = 65  # Typical CPU TDP
PUE = 1.58  # Power Usage Effectiveness
CARBON_INTENSITY = 0.475  # kgCO2/kWh (global average)

for name in models:
    model_cls = type(trained[name])
    
    # Measure training time
    t0 = time.time()
    m = model_cls(**trained[name].get_params())
    m.fit(X_train_s, y_train)
    train_time = time.time() - t0
    
    # Measure inference time
    t0 = time.time()
    for _ in range(100):
        m.predict(X_test_s)
    inference_time = (time.time() - t0) / 100
    
    # Energy = Power * Time
    train_energy_kwh = (TDP_WATTS * train_time) / (3600 * 1000)
    inference_energy_kwh = (TDP_WATTS * inference_time) / (3600 * 1000)
    
    # CO2 = Energy * PUE * Carbon Intensity
    train_co2_g = train_energy_kwh * PUE * CARBON_INTENSITY * 1000
    inference_co2_g = inference_energy_kwh * PUE * CARBON_INTENSITY * 1000
    
    # Annual projection (1M inferences/year)
    annual_co2_kg = inference_co2_g * 1_000_000 / 1000
    
    environmental_results[name] = {
        'train_time_s': round(train_time, 4),
        'inference_time_ms': round(inference_time * 1000, 4),
        'train_co2_g': round(train_co2_g, 6),
        'inference_co2_mg': round(inference_co2_g * 1000, 6),
        'annual_co2_kg': round(annual_co2_kg, 4),
    }
    print(f"\n{name}:")
    print(f"  Train: {train_time:.4f}s, CO2={train_co2_g:.6f}g")
    print(f"  Inference: {inference_time*1000:.4f}ms, CO2={inference_co2_g*1000:.6f}mg")
    print(f"  Annual (1M inferences): {annual_co2_kg:.4f} kgCO2")

RESULTS['environmental'] = environmental_results

# ============================================================
# 8. Integrated Ethics Score (EthicAI Score)
# ============================================================
print("\n" + "=" * 60)
print("Integrated EthicAI Score")
print("=" * 60)

ethics_scores = {}
weights = {'fairness': 0.25, 'explainability': 0.20, 'privacy': 0.20, 'robustness': 0.20, 'environmental': 0.15}

for name in models:
    f = fairness_results[name]
    e = explainability_results[name]
    p = privacy_results[name]
    r = robustness_results[name]
    env = environmental_results[name]
    
    fairness_score = (f['IFS_gender'] + f['IFS_ethnicity']) / 2
    explain_score = (e['shap_consistency'] + e['shap_stability'] + e['top_k_agreement']) / 3
    privacy_s = p['privacy_score']
    robust_score = (r['adversarial_robustness'] + r['consistency_score']) / 2
    # Environmental: normalize (lower is better, score inversely)
    max_annual = max(environmental_results[n]['annual_co2_kg'] for n in models)
    env_score = 1 - (env['annual_co2_kg'] / (max_annual + 1e-10)) if max_annual > 0 else 1.0
    
    ethicai = (
        weights['fairness'] * fairness_score +
        weights['explainability'] * explain_score +
        weights['privacy'] * privacy_s +
        weights['robustness'] * robust_score +
        weights['environmental'] * env_score
    )
    
    ethics_scores[name] = {
        'fairness_score': round(fairness_score, 4),
        'explainability_score': round(explain_score, 4),
        'privacy_score': round(privacy_s, 4),
        'robustness_score': round(robust_score, 4),
        'environmental_score': round(env_score, 4),
        'ethicai_score': round(ethicai, 4),
    }
    print(f"\n{name}: EthicAI = {ethicai:.4f}")
    print(f"  F={fairness_score:.4f} E={explain_score:.4f} P={privacy_s:.4f} R={robust_score:.4f} Env={env_score:.4f}")

RESULTS['ethics_scores'] = ethics_scores

# ============================================================
# 9. Generate Figures
# ============================================================
print("\n" + "=" * 60)
print("Generating Figures...")
print("=" * 60)

plt.rcParams.update({'font.size': 11})

# --- Figure 1: Fairness Metrics Comparison ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
model_names = list(models.keys())
x = np.arange(len(model_names))
width = 0.25

for ax, attr, title in [(axes[0], 'gender', 'Gender'), (axes[1], 'ethnicity', 'Ethnicity')]:
    spd = [abs(fairness_results[n][f'SPD_{attr}']) for n in model_names]
    eod = [abs(fairness_results[n][f'EOD_{attr}']) for n in model_names]
    cal = [fairness_results[n][f'Cal_{attr}'] for n in model_names]
    
    ax.bar(x - width, spd, width, label='|SPD|', color='#e74c3c', alpha=0.8)
    ax.bar(x, eod, width, label='|EOD|', color='#3498db', alpha=0.8)
    ax.bar(x + width, cal, width, label='CalDiff', color='#2ecc71', alpha=0.8)
    ax.axhline(y=0.1, color='gray', linestyle='--', alpha=0.5, label='Threshold (0.1)')
    ax.set_xlabel('Model')
    ax.set_ylabel('Disparity')
    ax.set_title(f'Fairness Metrics — {title}')
    ax.set_xticks(x)
    ax.set_xticklabels([n[:10] for n in model_names], rotation=15)
    ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/fairness_comparison.png', dpi=150, bbox_inches='tight')
plt.close()

# --- Figure 2: Radar Chart (EthicAI Scores) ---
categories = ['Fairness', 'Explainability', 'Privacy', 'Robustness', 'Environmental']
fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()
angles += angles[:1]
colors = ['#e74c3c', '#3498db', '#2ecc71']

for i, name in enumerate(model_names):
    scores = [
        ethics_scores[name]['fairness_score'],
        ethics_scores[name]['explainability_score'],
        ethics_scores[name]['privacy_score'],
        ethics_scores[name]['robustness_score'],
        ethics_scores[name]['environmental_score'],
    ]
    scores += scores[:1]
    ax.plot(angles, scores, 'o-', linewidth=2, label=name, color=colors[i])
    ax.fill(angles, scores, alpha=0.1, color=colors[i])

ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories)
ax.set_ylim(0, 1)
ax.set_title('EthicAI Radar — Multi-dimensional Ethics Evaluation', pad=20)
ax.legend(loc='lower right', bbox_to_anchor=(1.3, 0))
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/ethicai_radar.png', dpi=150, bbox_inches='tight')
plt.close()

# --- Figure 3: Robustness under Adversarial Perturbation ---
fig, ax = plt.subplots(figsize=(10, 6))
epsilons = [0.01, 0.05, 0.1, 0.2, 0.5]
for i, name in enumerate(model_names):
    accs = [robustness_results[name]['adversarial_accs'][str(e)] for e in epsilons]
    ax.plot(epsilons, accs, 'o-', linewidth=2, markersize=8, label=name, color=colors[i])

ax.set_xlabel('Perturbation Magnitude (ε)')
ax.set_ylabel('Accuracy')
ax.set_title('Model Robustness under Adversarial Perturbation')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/robustness_adversarial.png', dpi=150, bbox_inches='tight')
plt.close()

# --- Figure 4: Privacy Risk Comparison ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
mia_accs = [privacy_results[n]['mia_accuracy'] for n in model_names]
privacy_scores = [privacy_results[n]['privacy_score'] for n in model_names]
overfit = [privacy_results[n]['overfit_gap'] for n in model_names]

axes[0].bar(model_names, mia_accs, color=['#e74c3c', '#3498db', '#2ecc71'], alpha=0.8)
axes[0].axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='Random (0.5)')
axes[0].set_ylabel('MIA Accuracy')
axes[0].set_title('Membership Inference Attack Accuracy')
axes[0].legend()

axes[1].bar(model_names, overfit, color=['#e74c3c', '#3498db', '#2ecc71'], alpha=0.8)
axes[1].set_ylabel('Accuracy Gap')
axes[1].set_title('Overfitting Gap (Train - Test)')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/privacy_risk.png', dpi=150, bbox_inches='tight')
plt.close()

# --- Figure 5: Environmental Impact ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
train_co2 = [environmental_results[n]['train_co2_g'] for n in model_names]
annual_co2 = [environmental_results[n]['annual_co2_kg'] for n in model_names]

axes[0].bar(model_names, [c * 1000 for c in train_co2], color=['#e74c3c', '#3498db', '#2ecc71'], alpha=0.8)
axes[0].set_ylabel('Training CO₂ (mg)')
axes[0].set_title('Training Carbon Emission')

axes[1].bar(model_names, annual_co2, color=['#e74c3c', '#3498db', '#2ecc71'], alpha=0.8)
axes[1].set_ylabel('Annual CO₂ (kg)')
axes[1].set_title('Projected Annual CO₂ (1M inferences)')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/environmental_impact.png', dpi=150, bbox_inches='tight')
plt.close()

# --- Figure 6: Integrated EthicAI Score Summary ---
fig, ax = plt.subplots(figsize=(10, 6))
components = ['fairness_score', 'explainability_score', 'privacy_score', 'robustness_score', 'environmental_score']
comp_labels = ['Fairness', 'Explainability', 'Privacy', 'Robustness', 'Environmental']
bottom = np.zeros(len(model_names))
comp_colors = ['#e74c3c', '#f39c12', '#3498db', '#2ecc71', '#9b59b6']

for j, (comp, label) in enumerate(zip(components, comp_labels)):
    vals = [ethics_scores[n][comp] * list(weights.values())[j] for n in model_names]
    ax.bar(model_names, vals, bottom=bottom, label=label, color=comp_colors[j], alpha=0.85)
    bottom += vals

ax.set_ylabel('Weighted Score')
ax.set_title('EthicAI Integrated Score Breakdown')
ax.legend(loc='upper right')
ax.set_ylim(0, 1.1)

# Add total score labels
for i, name in enumerate(model_names):
    total = ethics_scores[name]['ethicai_score']
    ax.text(i, bottom[i] + 0.02, f'{total:.3f}', ha='center', fontweight='bold')

plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/ethicai_integrated.png', dpi=150, bbox_inches='tight')
plt.close()

# --- Figure 7: Medical Case Study Heatmap ---
fig, ax = plt.subplots(figsize=(10, 7))
matrix = []
row_labels = model_names
col_labels = ['Fair(G)', 'Fair(E)', 'Explain', 'Privacy', 'Robust', 'Environ', 'EthicAI']
for name in model_names:
    row = [
        ethics_scores[name]['fairness_score'],
        fairness_results[name]['IFS_ethnicity'],
        ethics_scores[name]['explainability_score'],
        ethics_scores[name]['privacy_score'],
        ethics_scores[name]['robustness_score'],
        ethics_scores[name]['environmental_score'],
        ethics_scores[name]['ethicai_score'],
    ]
    matrix.append(row)

matrix = np.array(matrix)
im = ax.imshow(matrix, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)
ax.set_xticks(np.arange(len(col_labels)))
ax.set_yticks(np.arange(len(row_labels)))
ax.set_xticklabels(col_labels, rotation=30, ha='right')
ax.set_yticklabels(row_labels)

for i in range(len(row_labels)):
    for j in range(len(col_labels)):
        ax.text(j, i, f'{matrix[i, j]:.3f}', ha='center', va='center', fontsize=11, fontweight='bold')

plt.colorbar(im, ax=ax, label='Score')
ax.set_title('Medical AI Ethics Audit — Comprehensive Evaluation Matrix')
plt.tight_layout()
plt.savefig(f'{FIGURES_DIR}/medical_audit_heatmap.png', dpi=150, bbox_inches='tight')
plt.close()

# Save results JSON
with open('results.json', 'w') as f:
    json.dump(RESULTS, f, indent=2, default=str)

print("\n" + "=" * 60)
print("All experiments completed. Results saved.")
print("=" * 60)
print(f"\nGenerated figures:")
for fig_file in sorted(os.listdir(FIGURES_DIR)):
    print(f"  - figures/{fig_file}")
