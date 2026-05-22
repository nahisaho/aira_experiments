#!/usr/bin/env python3
import json
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

np.random.seed(42)
random.seed(42)

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / 'results'
FIGURES_DIR = ROOT / 'figures'
RESULTS_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)

sns.set_theme(style='whitegrid')
COLORS = ['#0072B2', '#E69F00', '#009E73', '#CC79A7']
CLASS_MAP = {0: 'Healthy', 1: 'Parkinson', 2: 'ALS', 3: 'Cognitive'}


def simulate_group(label, n=40):
    gait_means = {
        0: [0.08, 0.07, 0.20, 0.85, 120, 112, 1.30],
        1: [0.18, 0.18, 0.62, 0.48, 78, 90, 0.94],
        2: [0.10, 0.09, 0.28, 0.72, 104, 108, 1.18],
        3: [0.11, 0.10, 0.26, 0.68, 98, 104, 1.12],
    }
    voice_means = {
        0: [0.36, 11.0, 0.18, 1.3, 23.5, 205, 8.0, 4.8],
        1: [0.42, 13.5, 0.22, 1.5, 21.5, 198, 9.5, 4.5],
        2: [1.15, 33.0, 0.70, 5.8, 13.5, 178, 18.0, 3.4],
        3: [0.55, 16.0, 0.28, 2.1, 20.0, 193, 10.5, 4.1],
    }
    touch_means = {
        0: [0.42, 0.10, 2.5, 520, 720, 0.48, 0.08, 0.04, 0.06, 2.8],
        1: [0.50, 0.13, 3.0, 485, 640, 0.52, 0.09, 0.05, 0.07, 2.7],
        2: [0.48, 0.12, 3.3, 470, 620, 0.50, 0.10, 0.06, 0.08, 2.6],
        3: [0.72, 0.24, 6.5, 355, 430, 0.76, 0.22, 0.16, 0.20, 1.9],
    }
    gait = np.random.normal(gait_means[label], [0.03, 0.04, 0.08, 0.08, 10, 8, 0.08], size=(n, 7))
    voice = np.random.normal(voice_means[label], [0.10, 3.5, 0.06, 0.5, 1.5, 10, 2.5, 0.3], size=(n, 8))
    touch = np.random.normal(touch_means[label], [0.06, 0.03, 0.8, 40, 65, 0.07, 0.03, 0.03, 0.04, 0.2], size=(n, 10))
    labels = np.full(n, label)
    return gait, voice, touch, labels


def late_fusion_probabilities(Xg, Xv, Xt, y, cv):
    base_models = [
        RandomForestClassifier(n_estimators=120, random_state=42),
        Pipeline([('scaler', StandardScaler()), ('mlp', MLPClassifier(hidden_layer_sizes=(20,), max_iter=500, random_state=42))]),
        RandomForestClassifier(n_estimators=120, random_state=42),
    ]
    probas = []
    for X, model in zip([Xg, Xv, Xt], base_models):
        probas.append(cross_val_predict(model, X, y, cv=cv, method='predict_proba'))
    weights = np.array([0.35, 0.30, 0.35])
    fused = sum(w * p for w, p in zip(weights, probas))
    return fused, probas


def main():
    generated = [simulate_group(label, 40) for label in range(4)]
    Xg = np.vstack([item[0] for item in generated])
    Xv = np.vstack([item[1] for item in generated])
    Xt = np.vstack([item[2] for item in generated])
    y = np.concatenate([item[3] for item in generated])

    X_early = np.hstack([Xg, Xv, Xt])
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    early_model = RandomForestClassifier(n_estimators=150, random_state=42)
    early_pred = cross_val_predict(early_model, X_early, y, cv=cv, method='predict')
    early_acc = accuracy_score(y, early_pred)

    late_proba, modality_probas = late_fusion_probabilities(Xg, Xv, Xt, y, cv)
    late_pred = late_proba.argmax(axis=1)
    late_acc = accuracy_score(y, late_pred)

    att_model = Pipeline([
        ('scaler', StandardScaler()),
        ('mlp', MLPClassifier(hidden_layer_sizes=(32, 16), max_iter=700, random_state=42)),
    ])
    att_pred = cross_val_predict(att_model, X_early, y, cv=cv, method='predict')
    att_proba = cross_val_predict(att_model, X_early, y, cv=cv, method='predict_proba')
    att_acc = accuracy_score(y, att_pred)

    disease_probability = 1 - att_proba[:, 0]
    severity_component = np.choose(att_pred, [0.10, 0.58, 0.72, 0.64])
    ndd_score = np.clip((0.65 * disease_probability + 0.35 * severity_component) * 100, 0, 100)

    comp_df = pd.DataFrame({
        'group': [CLASS_MAP[i] for i in y],
        'NDD_Score': ndd_score,
    })

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    method_df = pd.DataFrame({
        'Method': ['Early fusion', 'Late fusion', 'Attention fusion'],
        'Accuracy': [early_acc, late_acc, att_acc],
    })
    sns.barplot(data=method_df, x='Method', y='Accuracy', palette='viridis', ax=axes[0])
    axes[0].set_title('Fusion method comparison')
    axes[0].set_xlabel('Method')
    axes[0].set_ylabel('Accuracy')
    axes[0].set_ylim(0, 1.05)

    sns.violinplot(data=comp_df, x='group', y='NDD_Score', palette=COLORS, ax=axes[1])
    axes[1].set_title('Composite NDD-Score distribution')
    axes[1].set_xlabel('Group')
    axes[1].set_ylabel('NDD-Score (0-100)')
    axes[1].tick_params(axis='x', rotation=15)

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'multimodal_fusion.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    score_ranges = {'healthy_range': [0, 25], 'mild_risk': [25, 50], 'moderate_risk': [50, 75], 'high_risk': [75, 100]}
    group_stats = comp_df.groupby('group')['NDD_Score'].agg(['mean', 'std', 'min', 'max']).round(4).to_dict('index')

    results = {
        'n_samples': int(len(y)),
        'assumption': 'Four-class problem interpreted as 40 samples for each of Healthy, Parkinson, ALS, and Cognitive groups (160 total).',
        'fusion_accuracy': {
            'early_fusion': round(float(early_acc), 4),
            'late_fusion': round(float(late_acc), 4),
            'attention_fusion': round(float(att_acc), 4),
        },
        'composite_score_ranges': score_ranges,
        'composite_score_distribution': group_stats,
        'mean_modality_attention_proxy': {
            'gait_weight_proxy': round(float(np.mean(np.max(modality_probas[0], axis=1))), 4),
            'voice_weight_proxy': round(float(np.mean(np.max(modality_probas[1], axis=1))), 4),
            'touch_weight_proxy': round(float(np.mean(np.max(modality_probas[2], axis=1))), 4),
        },
    }

    with open(RESULTS_DIR / 'composite_score_results.json', 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)

    print(json.dumps(results['fusion_accuracy'], ensure_ascii=False))


if __name__ == '__main__':
    main()
