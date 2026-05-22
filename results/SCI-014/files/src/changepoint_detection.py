#!/usr/bin/env python3
import json
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

np.random.seed(42)
random.seed(42)

try:
    import ruptures as rpt
    RUPTURES_AVAILABLE = True
except Exception:
    RUPTURES_AVAILABLE = False

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / 'results'
FIGURES_DIR = ROOT / 'figures'
RESULTS_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)

sns.set_theme(style='whitegrid')
COLORS = ['#0072B2', '#E69F00', '#009E73', '#CC79A7']


def simulate_patient(change_month=None, baseline=40, delta=0):
    months = np.arange(1, 25)
    trend = np.random.normal(0.12, 0.15, size=months.size).cumsum()
    signal = baseline + trend + np.random.normal(0, 1.2, size=months.size)
    if change_month is not None:
        signal[change_month - 1:] += delta
    return months, signal


def detect_cusum(signal):
    baseline = np.mean(signal[:6])
    centered = signal - baseline
    pos = np.cumsum(centered - 0.4)
    neg = np.cumsum(-centered - 0.4)
    stat = np.maximum(np.abs(pos), np.abs(neg))
    idx = int(np.argmax(stat))
    if stat[idx] > 6.0:
        return [idx + 1]
    return []


def custom_pelt(signal, min_size=4):
    n = len(signal)
    candidates = []
    for cp in range(min_size, n - min_size):
        left = signal[:cp]
        right = signal[cp:]
        sse = ((left - left.mean()) ** 2).sum() + ((right - right.mean()) ** 2).sum()
        gain = ((signal - signal.mean()) ** 2).sum() - sse
        candidates.append((gain, cp))
    best_gain, best_cp = max(candidates, key=lambda x: x[0])
    return [best_cp + 1] if best_gain > 12 else []


def detect_pelt(signal):
    if RUPTURES_AVAILABLE:
        model = rpt.Pelt(model='l2', min_size=4).fit(signal)
        result = model.predict(pen=8)
        return [cp for cp in result[:-1]]
    return custom_pelt(signal)


def detect_bocpd(signal):
    baseline_mean = np.mean(signal[:5])
    baseline_std = np.std(signal[:5]) + 1e-6
    scores = []
    for t in range(5, len(signal)):
        window = signal[max(0, t - 3): t + 1]
        z = abs(window.mean() - baseline_mean) / baseline_std
        scores.append((z, t + 1))
    best_score, month = max(scores, key=lambda x: x[0])
    return [month] if best_score > 2.8 else []


def ensemble_vote(method_outputs, tolerance=1):
    all_points = [pt for pts in method_outputs.values() for pt in pts]
    if not all_points:
        return []
    voted = []
    for pt in sorted(all_points):
        support = sum(any(abs(pt - cand) <= tolerance for cand in pts) for pts in method_outputs.values())
        if support >= 2:
            voted.append(pt)
    if not voted:
        return []
    clusters = []
    for pt in voted:
        if not clusters or pt - clusters[-1][-1] > tolerance:
            clusters.append([pt])
        else:
            clusters[-1].append(pt)
    return [int(round(np.mean(cluster))) for cluster in clusters]


def evaluate(patient_truth, patient_detected):
    tolerance = 1
    outcomes = []
    lags = []
    false_positives = 0
    for truth, detected in zip(patient_truth, patient_detected):
        if truth is None:
            if detected:
                false_positives += len(detected)
                outcomes.append('miss')
            else:
                outcomes.append('hit')
        else:
            matches = [cp for cp in detected if abs(cp - truth) <= tolerance]
            if matches:
                cp = matches[0]
                outcomes.append('hit')
                lags.append(cp - truth)
                false_positives += max(0, len(detected) - 1)
            else:
                outcomes.append('miss')
                false_positives += len(detected)
    return outcomes, false_positives / len(patient_truth), float(np.mean(lags)) if lags else None


def main():
    patients = [
        ('Patient 1', None, 38, 0),
        ('Patient 2', 12, 43, 9),
        ('Patient 3', 8, 46, 11),
    ]

    patient_outputs = []
    truths = []
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

    for ax, (name, truth, baseline, delta) in zip(axes, patients):
        months, signal = simulate_patient(truth, baseline, delta)
        outputs = {
            'cusum': detect_cusum(signal),
            'pelt': detect_pelt(signal),
            'bocpd': detect_bocpd(signal),
        }
        ensemble = ensemble_vote(outputs)
        patient_outputs.append({'name': name, 'truth': truth, 'signal': signal.tolist(), 'methods': outputs, 'ensemble': ensemble})
        truths.append(truth)

        ax.plot(months, signal, marker='o', color=COLORS[0], label='Composite biomarker score')
        if truth is not None:
            ax.axvline(truth, color=COLORS[3], linestyle='--', linewidth=2, label='True change point')
        for cp in ensemble:
            ax.axvline(cp, color=COLORS[1], linestyle='-', linewidth=2, alpha=0.8, label='Ensemble detection')
        ax.set_title(name)
        ax.set_ylabel('Score')
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        ax.legend(unique.values(), unique.keys(), loc='upper left')

    axes[-1].set_xlabel('Month')
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'changepoint_detection.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    outcomes, false_positive_rate, mean_lag = evaluate(truths, [p['ensemble'] for p in patient_outputs])
    detection_accuracy = sum(o == 'hit' for o in outcomes) / len(outcomes)

    results = {
        'n_patients': len(patient_outputs),
        'ruptures_available': RUPTURES_AVAILABLE,
        'metrics': {
            'detection_accuracy': round(float(detection_accuracy), 4),
            'false_positive_rate': round(float(false_positive_rate), 4),
            'mean_detection_lag_months': None if mean_lag is None else round(float(mean_lag), 4),
            'patient_outcomes': outcomes,
        },
        'patients': patient_outputs,
    }

    with open(RESULTS_DIR / 'changepoint_results.json', 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)

    print(json.dumps(results['metrics'], ensure_ascii=False))


if __name__ == '__main__':
    main()
