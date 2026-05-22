#!/usr/bin/env python3
import json
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import f1_score, mean_squared_error, r2_score
from sklearn.model_selection import StratifiedKFold, KFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.svm import SVC, SVR

np.random.seed(42)
random.seed(42)

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / 'results'
FIGURES_DIR = ROOT / 'figures'
RESULTS_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)

sns.set_theme(style='whitegrid')
COLORS = ['#0072B2', '#E69F00', '#009E73']


def simulate_sample(label):
    params = {
        0: {'jitter': 0.35, 'shimmer': 0.18, 'hnr': 23.0, 'f0': 205, 'f0_std': 8, 'rate': 4.8, 'score': 45},
        1: {'jitter': 0.75, 'shimmer': 0.42, 'hnr': 18.0, 'f0': 190, 'f0_std': 14, 'rate': 4.0, 'score': 33},
        2: {'jitter': 1.30, 'shimmer': 0.78, 'hnr': 12.5, 'f0': 175, 'f0_std': 20, 'rate': 3.2, 'score': 22},
    }[label]
    mfcc_means = np.linspace(-2, 2, 13) + np.random.normal(label * 0.45, 0.45, 13)
    mfcc_stds = np.abs(np.random.normal(0.9 + label * 0.25, 0.18, 13))
    jitter_local = np.random.normal(params['jitter'], 0.09)
    jitter_absolute = np.random.normal(params['jitter'] * 30, 2.5)
    shimmer_local = np.random.normal(params['shimmer'], 0.05)
    shimmer_db = np.random.normal(params['shimmer'] * 8, 0.3)
    hnr = np.random.normal(params['hnr'], 1.5)
    f0_mean = np.random.normal(params['f0'], 12)
    f0_std = np.abs(np.random.normal(params['f0_std'], 2.5))
    speaking_rate = np.random.normal(params['rate'], 0.35)
    progression = np.clip(np.random.normal(params['score'], 4.5), 5, 48)

    row = {
        'jitter_local': jitter_local,
        'jitter_absolute': jitter_absolute,
        'shimmer_local': shimmer_local,
        'shimmer_db': shimmer_db,
        'hnr': hnr,
        'f0_mean': f0_mean,
        'f0_std': f0_std,
        'speaking_rate': speaking_rate,
        'alsfrs_proxy': progression,
        'label': label,
    }
    for i, val in enumerate(mfcc_means, start=1):
        row[f'mfcc_{i}_mean'] = val
    for i, val in enumerate(mfcc_stds, start=1):
        row[f'mfcc_{i}_std'] = val
    return row


def rmse(y_true, y_pred):
    return float(np.sqrt(mean_squared_error(y_true, y_pred)))


def main():
    rows = [simulate_sample(1) for _ in range(50)] + [simulate_sample(2) for _ in range(50)] + [simulate_sample(0) for _ in range(50)]
    df = pd.DataFrame(rows)
    df = df[['label', 'alsfrs_proxy'] + [c for c in df.columns if c not in {'label', 'alsfrs_proxy'}]]
    feature_cols = [c for c in df.columns if c not in {'label', 'alsfrs_proxy'}]
    X = df[feature_cols].values
    y = df['label'].values
    y_reg = df['alsfrs_proxy'].values

    clf = Pipeline([
        ('scaler', StandardScaler()),
        ('svc', SVC(kernel='rbf', C=10, gamma='scale', probability=True, random_state=42)),
    ])
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    pred = cross_val_predict(clf, X, y, cv=cv, method='predict')
    accuracy = float((pred == y).mean())
    per_class_f1 = f1_score(y, pred, average=None)

    reg = Pipeline([
        ('scaler', StandardScaler()),
        ('svr', SVR(kernel='rbf', C=10, gamma='scale')),
    ])
    reg_pred = cross_val_predict(reg, X, y_reg, cv=KFold(n_splits=5, shuffle=True, random_state=42))
    r2 = r2_score(y_reg, reg_pred)
    reg_rmse = rmse(y_reg, reg_pred)

    mfcc_cols = [f'mfcc_{i}_mean' for i in range(1, 14)]
    heatmap_data = df.groupby('label')[mfcc_cols].mean().rename(index={0: 'Healthy', 1: 'Early ALS', 2: 'Moderate ALS'})

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    sns.heatmap(heatmap_data, cmap='viridis', ax=axes[0], cbar_kws={'label': 'Mean MFCC'})
    axes[0].set_title('MFCC profile by class')
    axes[0].set_xlabel('MFCC coefficient')
    axes[0].set_ylabel('Class')

    scatter_df = df.copy()
    scatter_df['class'] = scatter_df['label'].map({0: 'Healthy', 1: 'Early ALS', 2: 'Moderate ALS'})
    sns.scatterplot(data=scatter_df, x='jitter_local', y='shimmer_local', hue='class', palette=COLORS, ax=axes[1])
    axes[1].set_title('Jitter and shimmer separation')
    axes[1].set_xlabel('Jitter local (%)')
    axes[1].set_ylabel('Shimmer local (%)')

    axes[2].scatter(y_reg, reg_pred, color=COLORS[0], alpha=0.75)
    lims = [min(y_reg.min(), reg_pred.min()) - 1, max(y_reg.max(), reg_pred.max()) + 1]
    axes[2].plot(lims, lims, '--', color='gray')
    axes[2].set_xlim(lims)
    axes[2].set_ylim(lims)
    axes[2].set_title('ALS progression regression')
    axes[2].set_xlabel('Observed ALSFRS-R proxy')
    axes[2].set_ylabel('Predicted ALSFRS-R proxy')

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'als_voice_features.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    results = {
        'n_samples': int(len(df)),
        'metrics': {
            'multiclass_accuracy': round(accuracy, 4),
            'per_class_f1': {
                'healthy': round(float(per_class_f1[0]), 4),
                'early_als': round(float(per_class_f1[1]), 4),
                'moderate_als': round(float(per_class_f1[2]), 4),
            },
            'regression_r2': round(float(r2), 4),
            'regression_rmse': round(float(reg_rmse), 4),
        },
        'class_feature_means': df.groupby('label')[feature_cols + ['alsfrs_proxy']].mean().round(4).rename(index={0: 'healthy', 1: 'early_als', 2: 'moderate_als'}).to_dict(),
    }

    with open(RESULTS_DIR / 'als_monitoring_results.json', 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)

    print(json.dumps(results['metrics'], ensure_ascii=False))


if __name__ == '__main__':
    main()
