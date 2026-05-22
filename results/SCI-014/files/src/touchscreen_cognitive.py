#!/usr/bin/env python3
import json
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import accuracy_score, recall_score, roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import label_binarize

np.random.seed(42)
random.seed(42)

try:
    from xgboost import XGBClassifier
    XGB_AVAILABLE = True
except Exception:
    XGB_AVAILABLE = False

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / 'results'
FIGURES_DIR = ROOT / 'figures'
RESULTS_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)

sns.set_theme(style='whitegrid')
COLORS = ['#0072B2', '#E69F00', '#009E73']
CLASS_NAMES = ['Control', 'MCI', 'Mild dementia']


def simulate_sample(label):
    params = {
        0: {'iti': 0.42, 'iti_cv': 0.10, 'accuracy': 2.5, 'swipe_v': 520, 'swipe_a': 720, 'press': 0.48, 'dtap': 0.08, 'error': 0.04, 'corr': 0.06, 'entropy': 2.8},
        1: {'iti': 0.58, 'iti_cv': 0.18, 'accuracy': 4.8, 'swipe_v': 430, 'swipe_a': 560, 'press': 0.62, 'dtap': 0.15, 'error': 0.10, 'corr': 0.14, 'entropy': 2.3},
        2: {'iti': 0.76, 'iti_cv': 0.26, 'accuracy': 7.2, 'swipe_v': 340, 'swipe_a': 420, 'press': 0.81, 'dtap': 0.23, 'error': 0.18, 'corr': 0.23, 'entropy': 1.8},
    }[label]
    return {
        'inter_tap_interval_mean': np.random.normal(params['iti'], 0.06),
        'inter_tap_interval_cv': np.abs(np.random.normal(params['iti_cv'], 0.025)),
        'tap_accuracy': np.abs(np.random.normal(params['accuracy'], 0.8)),
        'swipe_velocity': np.random.normal(params['swipe_v'], 55),
        'swipe_acceleration': np.random.normal(params['swipe_a'], 80),
        'long_press_duration': np.random.normal(params['press'], 0.08),
        'double_tap_variability': np.abs(np.random.normal(params['dtap'], 0.03)),
        'error_rate': np.clip(np.random.normal(params['error'], 0.03), 0, 1),
        'correction_frequency': np.clip(np.random.normal(params['corr'], 0.04), 0, 1),
        'typing_rhythm_entropy': np.random.normal(params['entropy'], 0.22),
        'label': label,
    }


def build_model():
    if XGB_AVAILABLE:
        return XGBClassifier(
            n_estimators=100,
            max_depth=4,
            random_state=42,
            eval_metric='mlogloss',
            learning_rate=0.08,
            subsample=0.9,
            colsample_bytree=0.9,
        )
    return GradientBoostingClassifier(random_state=42)


def main():
    rows = [simulate_sample(1) for _ in range(60)] + [simulate_sample(2) for _ in range(60)] + [simulate_sample(0) for _ in range(60)]
    df = pd.DataFrame(rows)
    feature_cols = [c for c in df.columns if c != 'label']
    X = df[feature_cols].values
    y = df['label'].values

    model = build_model()
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    pred = cross_val_predict(model, X, y, cv=cv, method='predict')
    proba = cross_val_predict(model, X, y, cv=cv, method='predict_proba')

    accuracy = accuracy_score(y, pred)
    y_bin = label_binarize(y, classes=[0, 1, 2])
    auc = roc_auc_score(y_bin, proba, multi_class='ovr', average='macro')
    recalls = recall_score(y, pred, average=None)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    plot_df = df.copy()
    plot_df['class'] = plot_df['label'].map({0: 'Control', 1: 'MCI', 2: 'Mild dementia'})
    sns.violinplot(data=plot_df, x='class', y='inter_tap_interval_mean', palette=COLORS, ax=axes[0])
    axes[0].set_title('Inter-tap interval by class')
    axes[0].set_xlabel('Class')
    axes[0].set_ylabel('Inter-tap interval mean (s)')

    for idx, class_name in enumerate(CLASS_NAMES):
        fpr, tpr, _ = roc_curve(y_bin[:, idx], proba[:, idx])
        axes[1].plot(fpr, tpr, linewidth=2, color=COLORS[idx], label=class_name)
    axes[1].plot([0, 1], [0, 1], '--', color='gray')
    axes[1].set_title('One-vs-rest ROC curves')
    axes[1].set_xlabel('False positive rate')
    axes[1].set_ylabel('True positive rate')
    axes[1].legend(loc='lower right')

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'touchscreen_cognitive.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    results = {
        'n_samples': int(len(df)),
        'model': 'XGBoost' if XGB_AVAILABLE else 'GradientBoostingClassifier',
        'metrics': {
            'accuracy': round(float(accuracy), 4),
            'auc_roc_ovr': round(float(auc), 4),
            'per_class_recall': {
                'control': round(float(recalls[0]), 4),
                'mci': round(float(recalls[1]), 4),
                'mild_dementia': round(float(recalls[2]), 4),
            },
        },
        'feature_summary': df.groupby('label')[feature_cols].mean().round(4).rename(index={0: 'control', 1: 'mci', 2: 'mild_dementia'}).to_dict(),
    }

    with open(RESULTS_DIR / 'cognitive_decline_results.json', 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)

    print(json.dumps(results['metrics'], ensure_ascii=False))


if __name__ == '__main__':
    main()
