"""
Module 3: Mendelian Randomization (MR) Analysis for Drug Target Validation
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ─────────────────────────────────────────────────────────────
# Mendelian Randomization: IVW, Egger, Weighted Median
# ─────────────────────────────────────────────────────────────
def mr_ivw(beta_exp, beta_out, se_out):
    """Inverse-variance weighted MR estimator."""
    w = 1 / se_out**2
    beta_ivw = np.sum(w * beta_out * beta_exp) / np.sum(w * beta_exp**2)
    se_ivw   = np.sqrt(1 / np.sum(w * beta_exp**2))
    z         = beta_ivw / se_ivw
    p_val     = 2 * (1 - stats.norm.cdf(abs(z)))
    ci_lo     = beta_ivw - 1.96 * se_ivw
    ci_hi     = beta_ivw + 1.96 * se_ivw
    return {'estimate': beta_ivw, 'se': se_ivw, 'z': z, 'p': p_val,
            'ci_lower': ci_lo, 'ci_upper': ci_hi}

def mr_egger(beta_exp, beta_out, se_out):
    """MR-Egger with intercept test for directional pleiotropy."""
    w       = 1 / se_out**2
    bx_w    = beta_exp * np.sqrt(w)
    by_w    = beta_out * np.sqrt(w)
    ones_w  = np.sqrt(w)
    X       = np.column_stack([ones_w, bx_w])
    res     = np.linalg.lstsq(X, by_w, rcond=None)
    coef    = res[0]
    intercept, slope = coef[0], coef[1]
    # residual SE
    y_pred  = X @ coef
    resid   = by_w - y_pred
    rss     = np.sum(resid**2)
    n       = len(beta_exp)
    sigma2  = rss / (n - 2)
    XtX_inv = np.linalg.inv(X.T @ X)
    se_coef = np.sqrt(np.diag(sigma2 * XtX_inv))
    se_int, se_slope = se_coef[0], se_coef[1]
    p_slope = 2 * (1 - stats.t.cdf(abs(slope / se_slope), df=n-2))
    p_int   = 2 * (1 - stats.t.cdf(abs(intercept / se_int), df=n-2))
    return {
        'estimate': slope, 'se': se_slope, 'p': p_slope,
        'ci_lower': slope - 1.96*se_slope, 'ci_upper': slope + 1.96*se_slope,
        'intercept': intercept, 'se_intercept': se_int, 'p_intercept': p_int,
    }

def mr_weighted_median(beta_exp, beta_out, se_out, n_boot=1000):
    """Weighted median estimator (robust to up to 50% invalid IVs)."""
    ratio    = beta_out / beta_exp
    w        = (se_out / beta_exp)**(-2)
    w        = w / w.sum()
    idx_sort = np.argsort(ratio)
    ratio_s  = ratio[idx_sort]; w_s = w[idx_sort]
    cum_w    = np.cumsum(w_s)
    median_r = ratio_s[np.searchsorted(cum_w, 0.5)]
    # Bootstrap SE
    boot_est = []
    for _ in range(n_boot):
        samp = np.random.choice(len(ratio), len(ratio), replace=True)
        br = beta_out[samp] / beta_exp[samp]
        bw = (se_out[samp] / beta_exp[samp])**(-2)
        bw = bw / bw.sum()
        si = np.argsort(br); br_s = br[si]; bw_s = bw[si]
        boot_est.append(br_s[np.searchsorted(np.cumsum(bw_s), 0.5)])
    se = np.std(boot_est)
    z  = median_r / se
    p  = 2 * (1 - stats.norm.cdf(abs(z)))
    return {'estimate': median_r, 'se': se, 'z': z, 'p': p,
            'ci_lower': median_r - 1.96*se, 'ci_upper': median_r + 1.96*se}

# ─────────────────────────────────────────────────────────────
# Simulate GWAS summary statistics for drug targets
# ─────────────────────────────────────────────────────────────
drug_targets = {
    'PCSK9_LDL_CAD': {
        'description': 'PCSK9 inhibitors: LDL reduction → CAD risk',
        'n_ivs': 12, 'true_effect': -0.36,
        'pleiotropy': False,
    },
    'IL6R_CRP_CAD': {
        'description': 'IL-6 receptor: CRP → CAD risk',
        'n_ivs': 8, 'true_effect': -0.28,
        'pleiotropy': False,
    },
    'HMGCR_LDL_T2D': {
        'description': 'HMGCR (statin target): LDL → T2D risk',
        'n_ivs': 10, 'true_effect': 0.12,
        'pleiotropy': True,
    },
    'GLP1R_BMI_T2D': {
        'description': 'GLP1R (GLP-1 agonists): BMI → T2D risk',
        'n_ivs': 15, 'true_effect': -0.45,
        'pleiotropy': False,
    },
}

all_results = {}
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

for ax_idx, (target_name, target_info) in enumerate(drug_targets.items()):
    n_ivs  = target_info['n_ivs']
    beta_exp = np.random.uniform(0.01, 0.05, n_ivs)
    se_exp   = beta_exp * np.random.uniform(0.05, 0.15, n_ivs)
    f_stats  = (beta_exp / se_exp)**2
    true_eff = target_info['true_effect']

    beta_out = true_eff * beta_exp + np.random.normal(0, 0.005, n_ivs)
    if target_info['pleiotropy']:
        pleiotropic_ivs = np.random.choice(n_ivs, n_ivs // 3, replace=False)
        beta_out[pleiotropic_ivs] += np.random.normal(0.02, 0.01, len(pleiotropic_ivs))

    se_out = np.abs(beta_out) * np.random.uniform(0.1, 0.3, n_ivs) + 0.002

    res_ivw = mr_ivw(beta_exp, beta_out, se_out)
    res_egger = mr_egger(beta_exp, beta_out, se_out)
    res_wm = mr_weighted_median(beta_exp, beta_out, se_out)

    all_results[target_name] = {
        'description': target_info['description'],
        'n_instruments': n_ivs,
        'mean_F_statistic': float(f_stats.mean()),
        'IVW': {k: float(v) for k, v in res_ivw.items()},
        'MR_Egger': {k: float(v) for k, v in res_egger.items()},
        'Weighted_Median': {k: float(v) for k, v in res_wm.items()},
        'directional_pleiotropy': target_info['pleiotropy'],
    }

    ax = axes[ax_idx]
    # Scatter plot: IV beta_exp vs beta_out
    ax.scatter(beta_exp, beta_out, color='steelblue', s=60, zorder=3,
               label='Genetic instruments')
    x_range = np.linspace(0, beta_exp.max()*1.1, 100)
    ax.plot(x_range, res_ivw['estimate'] * x_range,
            'r-', lw=2, label=f"IVW: β={res_ivw['estimate']:.3f} (p={res_ivw['p']:.3e})")
    ax.plot(x_range, res_egger['intercept'] + res_egger['estimate'] * x_range,
            'g--', lw=2, label=f"Egger: β={res_egger['estimate']:.3f}")
    ax.plot(x_range, res_wm['estimate'] * x_range,
            'm:', lw=2, label=f"WM: β={res_wm['estimate']:.3f}")
    ax.axhline(0, color='gray', lw=0.5, alpha=0.5)
    ax.axvline(0, color='gray', lw=0.5, alpha=0.5)
    ax.set_xlabel('SNP-Exposure Effect (β_exp)')
    ax.set_ylabel('SNP-Outcome Effect (β_out)')
    ax.set_title(f'{target_name.replace("_"," ")}\n{target_info["description"][:50]}...',
                  fontsize=9)
    ax.legend(fontsize=7, loc='upper left')

plt.suptitle('Mendelian Randomization: Drug Target Validation', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/figures/fig5_mr_analysis.png',
            dpi=150, bbox_inches='tight')
plt.close()

# ─────────────────────────────────────────────────────────────
# Forest plot: MR estimates comparison
# ─────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))
methods = ['IVW', 'MR_Egger', 'Weighted_Median']
method_colors = {'IVW': '#d73027', 'MR_Egger': '#4575b4', 'Weighted_Median': '#1a9850'}
method_offsets = {'IVW': -0.15, 'MR_Egger': 0, 'Weighted_Median': 0.15}

y_ticks = []; y_labels = []
for i, (target_name, res) in enumerate(all_results.items()):
    y_base = i * 1.2
    for method in methods:
        y = y_base + method_offsets[method]
        est = res[method]['estimate']
        lo  = res[method]['ci_lower']
        hi  = res[method]['ci_upper']
        ax.plot([lo, hi], [y, y], color=method_colors[method], lw=2)
        ax.scatter([est], [y], color=method_colors[method], s=60, zorder=3)
    y_ticks.append(y_base)
    y_labels.append(target_name.replace('_', ' '))

ax.axvline(0, color='black', lw=1, linestyle='--', alpha=0.5)
ax.set_yticks(y_ticks)
ax.set_yticklabels(y_labels, fontsize=9)
ax.set_xlabel('MR Effect Estimate (β)')
ax.set_title('MR Analysis: Drug Target Validation – Forest Plot\n(IVW, MR-Egger, Weighted Median)')
patches = [mpatches.Patch(color=c, label=m) for m, c in method_colors.items()]
ax.legend(handles=patches, loc='lower right')
plt.tight_layout()
plt.savefig('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/figures/fig6_mr_forest_plot.png',
            dpi=150, bbox_inches='tight')
plt.close()

with open('/app/projects/fc776f29-54af-4cdf-a794-8bdabd570a60/workspace/results/mr_analysis_results.json', 'w') as f:
    json.dump(all_results, f, indent=2, ensure_ascii=False)

print("[MR Module] Done")
for name, res in all_results.items():
    sig = "✓" if res['IVW']['p'] < 0.05 else "✗"
    print(f"  {name}: IVW β={res['IVW']['estimate']:.3f} (p={res['IVW']['p']:.3e}) {sig}")
