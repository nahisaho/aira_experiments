#!/usr/bin/env python3
import json
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (accuracy_score, confusion_matrix, recall_score, roc_auc_score,
                             roc_curve)
from sklearn.model_selection import StratifiedKFold, cross_val_predict

np.random.seed(42)
random.seed(42)

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / 'results'
FIGURES_DIR = ROOT / 'figures'
RESULTS_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)

sns.set_theme(style='whitegrid')
PALETTE = sns.color_palette('viridis', 6)


def bandpower_ratio(signal, fs, low_band=(0.5, 3.0), high_band=(3.0, 8.0)):
    freqs = np.fft.rfftfreq(signal.size, d=1 / fs)
    psd = np.abs(np.fft.rfft(signal)) ** 2
    low_mask = (freqs >= low_band[0]) & (freqs < low_band[1])
    high_mask = (freqs >= high_band[0]) & (freqs < high_band[1])
    low_power = np.trapz(psd[low_mask], freqs[low_mask]) if low_mask.any() else 0.0
    high_power = np.trapz(psd[high_mask], freqs[high_mask]) if high_mask.any() else 0.0
    return float(high_power / (low_power + 1e-6))


def autocorr_peak(signal):
    signal = signal - signal.mean()
    ac = np.correlate(signal, signal, mode='full')[signal.size - 1:]
    ac /= ac[0] + 1e-6
    search_window = ac[20:120]
    return float(np.max(search_window)) if search_window.size else 0.0


def peak_count(signal, threshold):
    peaks = []
    for i in range(1, len(signal) - 1):
        if signal[i] > threshold and signal[i] > signal[i - 1] and signal[i] > signal[i + 1]:
            peaks.append(i)
    return peaks


def simulate_subject(label):
    fs = 50
    duration = 12
    t = np.arange(0, duration, 1 / fs)
    is_pd = label == 1

    mean_stride = np.random.normal(1.18 if is_pd else 1.0, 0.06)
    stride_sd = np.random.uniform(0.13, 0.22) if is_pd else np.random.uniform(0.03, 0.08)
    n_intervals = int(duration / max(mean_stride, 0.6)) + 5
    stride_intervals = np.clip(np.random.normal(mean_stride, stride_sd, n_intervals), 0.5, 2.0)
    left_intervals = stride_intervals[::2]
    right_intervals = stride_intervals[1::2]
    if right_intervals.size == 0:
        right_intervals = left_intervals.copy()

    cadence = 60 / np.mean(stride_intervals)
    walking_speed = np.random.normal(0.92 if is_pd else 1.28, 0.08)
    step_freq = cadence / 60

    base_signal = np.sin(2 * np.pi * step_freq * t) + 0.35 * np.sin(4 * np.pi * step_freq * t)
    fog_component = (0.55 if is_pd else 0.15) * np.sin(2 * np.pi * np.random.uniform(3.5, 6.5) * t)
    noise = np.random.normal(0, 0.35 if is_pd else 0.18, size=t.size)
    accel = base_signal + fog_component + noise

    turn_base = np.random.normal(70 if is_pd else 110, 15)
    gyro = np.random.normal(0, 8 if is_pd else 6, size=t.size)
    turn_centers = np.random.choice(np.arange(75, t.size - 75), size=3, replace=False)
    for center in turn_centers:
        width = np.random.randint(8, 18)
        amplitude = np.random.normal(turn_base, 8)
        gyro += amplitude * np.exp(-0.5 * ((np.arange(t.size) - center) / width) ** 2)

    stride_cv = float(np.std(stride_intervals) / (np.mean(stride_intervals) + 1e-6))
    asymmetry = float(np.mean(np.abs(left_intervals[: right_intervals.size] - right_intervals[: left_intervals.size])) /
                      (0.5 * (np.mean(left_intervals) + np.mean(right_intervals)) + 1e-6) * 100)
    fog_ratio = bandpower_ratio(accel, fs)
    regularity = autocorr_peak(accel)
    turn_peaks = peak_count(gyro, np.percentile(gyro, 90))
    turning_speed = float(np.mean(gyro[turn_peaks])) if turn_peaks else float(np.max(gyro))

    return {
        'stride_length_variability': stride_cv,
        'gait_asymmetry_index': asymmetry,
        'fog_power_ratio': fog_ratio,
        'step_regularity': regularity,
        'turning_speed': turning_speed,
        'cadence': float(cadence),
        'walking_speed': float(walking_speed),
        'label': label,
    }


def main():
    rows = [simulate_subject(1) for _ in range(100)] + [simulate_subject(0) for _ in range(100)]
    df = pd.DataFrame(rows)
    feature_cols = [c for c in df.columns if c != 'label']
    X = df[feature_cols].values
    y = df['label'].values

    model = RandomForestClassifier(n_estimators=100, random_state=42)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    pred = cross_val_predict(model, X, y, cv=cv, method='predict')
    proba = cross_val_predict(model, X, y, cv=cv, method='predict_proba')[:, 1]

    accuracy = accuracy_score(y, pred)
    auc = roc_auc_score(y, proba)
    sensitivity = recall_score(y, pred, pos_label=1)
    specificity = recall_score(y, pred, pos_label=0)
    cm = confusion_matrix(y, pred)

    model.fit(X, y)
    importances = dict(sorted(zip(feature_cols, model.feature_importances_), key=lambda x: x[1], reverse=True))
    fpr, tpr, _ = roc_curve(y, proba)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    long_df = df.melt(id_vars='label', value_vars=['stride_length_variability', 'fog_power_ratio', 'walking_speed'],
                      var_name='feature', value_name='value')
    long_df['group'] = long_df['label'].map({0: 'Control', 1: 'Parkinson'})
    sns.boxplot(data=long_df, x='feature', y='value', hue='group', palette=[PALETTE[1], PALETTE[4]], ax=axes[0])
    axes[0].set_title('Key gait feature distributions')
    axes[0].set_xlabel('Feature')
    axes[0].set_ylabel('Value')
    axes[0].tick_params(axis='x', rotation=20)

    axes[1].plot(fpr, tpr, color=PALETTE[4], linewidth=2, label=f'ROC AUC = {auc:.3f}')
    axes[1].plot([0, 1], [0, 1], linestyle='--', color='gray')
    axes[1].set_title('ROC curve')
    axes[1].set_xlabel('False positive rate')
    axes[1].set_ylabel('True positive rate')
    axes[1].legend(loc='lower right')

    imp_series = pd.Series(importances)
    sns.barplot(x=imp_series.values, y=imp_series.index, palette='viridis', ax=axes[2])
    axes[2].set_title('Feature importance')
    axes[2].set_xlabel('Importance')
    axes[2].set_ylabel('Feature')

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'parkinson_gait_analysis.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    results = {
        'n_samples': int(len(df)),
        'metrics': {
            'accuracy': round(float(accuracy), 4),
            'auc_roc': round(float(auc), 4),
            'sensitivity': round(float(sensitivity), 4),
            'specificity': round(float(specificity), 4),
        },
        'feature_importances': {k: round(float(v), 4) for k, v in importances.items()},
        'confusion_matrix': cm.tolist(),
        'feature_summary': df.groupby('label')[feature_cols].mean().round(4).rename(index={0: 'control', 1: 'parkinson'}).to_dict(),
    }

    with open(RESULTS_DIR / 'parkinson_screening_results.json', 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)

    print(json.dumps(results['metrics'], ensure_ascii=False))


if __name__ == '__main__':
    main()
