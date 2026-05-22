#!/usr/bin/env python3
import json
import math
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats

np.random.seed(42)
random.seed(42)

ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / 'results'
FIGURES_DIR = ROOT / 'figures'
RESULTS_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)

sns.set_theme(style='whitegrid')
COLORS = ['#0072B2', '#E69F00', '#009E73', '#CC79A7']


def bootstrap_corr_ci(x, y, n_boot=1000):
    corrs = []
    n = len(x)
    for _ in range(n_boot):
        idx = np.random.randint(0, n, n)
        corrs.append(np.corrcoef(x[idx], y[idx])[0, 1])
    lower, upper = np.percentile(corrs, [2.5, 97.5])
    return float(np.corrcoef(x, y)[0, 1]), float(lower), float(upper)


def concordance_corrcoef(x, y):
    x_mean, y_mean = np.mean(x), np.mean(y)
    cov = np.cov(x, y, ddof=1)[0, 1]
    x_var = np.var(x, ddof=1)
    y_var = np.var(y, ddof=1)
    return float((2 * cov) / (x_var + y_var + (x_mean - y_mean) ** 2 + 1e-12))


def standardized_response_mean(baseline, follow_up):
    change = follow_up - baseline
    return float(np.mean(change) / (np.std(change, ddof=1) + 1e-12))


def required_sample_size(effect_size, alpha=0.05, power=0.8):
    z_alpha = stats.norm.ppf(1 - alpha / 2)
    z_beta = stats.norm.ppf(power)
    return int(math.ceil(2 * ((z_alpha + z_beta) / effect_size) ** 2))


def main():
    n_each = 40
    pd_digital = np.random.normal(58, 12, n_each)
    pd_clinical = pd_digital * 0.92 + np.random.normal(0, 6, n_each)

    als_digital = np.random.normal(52, 11, n_each)
    als_impairment = als_digital * 0.88 + np.random.normal(0, 7, n_each)
    als_clinical = 48 - np.clip(als_impairment / 2.3, 0, 40)
    als_clinical_scaled = 48 - als_clinical

    cog_digital = np.random.normal(49, 13, n_each)
    cog_impairment = cog_digital * 0.90 + np.random.normal(0, 6, n_each)
    moca = 30 - np.clip(cog_impairment / 2.5, 0, 24)
    cog_clinical_scaled = 30 - moca

    modalities = {
        'Parkinson_vs_UPDRSIII': (pd_digital, pd_clinical),
        'ALS_vs_ALSFRSR_impairment': (als_digital, als_clinical_scaled),
        'Cognitive_vs_MoCA_impairment': (cog_digital, cog_clinical_scaled),
    }

    correlation_table = {}
    ccc_table = {}
    srm_table = {}
    power_table = {}

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    scatter_axes = axes.flatten()[:3]
    bland_ax = axes.flatten()[3]

    bland_x = []
    bland_y = []

    for idx, (ax, (name, (digital, clinical))) in enumerate(zip(scatter_axes, modalities.items())):
        r, ci_low, ci_high = bootstrap_corr_ci(digital, clinical, n_boot=1000)
        ccc = concordance_corrcoef(digital, clinical)
        baseline = digital + np.random.normal(0, 2.5, len(digital))
        follow_up = baseline + np.random.normal(0.2 * np.mean(baseline), 5.0, len(digital))
        srm = standardized_response_mean(baseline, follow_up)
        effect_size = (0.2 * np.mean(baseline)) / (np.std(baseline, ddof=1) + 1e-12)
        sample_size = required_sample_size(effect_size)

        correlation_table[name] = {
            'pearson_r': round(r, 4),
            'ci_95': [round(ci_low, 4), round(ci_high, 4)],
        }
        ccc_table[name] = round(ccc, 4)
        srm_table[name] = round(srm, 4)
        power_table[name] = sample_size

        ax.scatter(digital, clinical, alpha=0.75, color=COLORS[idx])
        coef = np.polyfit(digital, clinical, 1)
        xs = np.linspace(digital.min(), digital.max(), 100)
        ax.plot(xs, coef[0] * xs + coef[1], color='black', linewidth=1.5)
        ax.set_title(name.replace('_', ' '))
        ax.set_xlabel('Digital score')
        ax.set_ylabel('Clinical score')

        bland_x.extend(((digital + clinical) / 2).tolist())
        bland_y.extend((digital - clinical).tolist())

    bland_x = np.array(bland_x)
    bland_y = np.array(bland_y)
    bias = bland_y.mean()
    loa = 1.96 * bland_y.std(ddof=1)
    bland_ax.scatter(bland_x, bland_y, alpha=0.6, color=COLORS[3])
    bland_ax.axhline(bias, color='black', linestyle='-')
    bland_ax.axhline(bias + loa, color='gray', linestyle='--')
    bland_ax.axhline(bias - loa, color='gray', linestyle='--')
    bland_ax.set_title('Bland-Altman agreement')
    bland_ax.set_xlabel('Mean score')
    bland_ax.set_ylabel('Digital - clinical')

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / 'validation_strategy.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    results = {
        'n_patients': 120,
        'correlations': correlation_table,
        'ccc': ccc_table,
        'srm': srm_table,
        'power_analysis_required_sample_size': power_table,
        'bland_altman': {
            'bias': round(float(bias), 4),
            'upper_limit_of_agreement': round(float(bias + loa), 4),
            'lower_limit_of_agreement': round(float(bias - loa), 4),
        },
    }

    with open(RESULTS_DIR / 'validation_results.json', 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)

    print(json.dumps({'ccc': results['ccc'], 'power': power_table}, ensure_ascii=False))


if __name__ == '__main__':
    main()
