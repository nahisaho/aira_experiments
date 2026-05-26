#!/usr/bin/env python3
"""
Vaccine Effectiveness (VE) Estimation Framework
================================================
A comprehensive methodological pipeline for estimating VE from real-world data.

Components:
1. Test-Negative Design (TND) simulation and analysis
2. Waning VE estimation with time-varying models
3. Variant-specific VE estimation
4. Healthy vaccinee bias correction (IPW)
5. Booster dose causal effect estimation (MSM)
6. mRNA vaccine hospitalization prevention case study

Uses Python equivalents of R's survival and gnm packages:
- lifelines (survival analysis)
- statsmodels (GLM, logistic regression)
- scipy (statistical tests)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.special import expit
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import survival_table_from_events
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod import families
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
sns.set_style("whitegrid")
plt.rcParams.update({'font.size': 11, 'figure.dpi': 150})

FIGDIR = "figures"
results = {}


# ============================================================
# UTILITY FUNCTIONS
# ============================================================

def generate_baseline_population(n=10000):
    """Generate a synthetic population with covariates."""
    age = np.random.normal(55, 15, n).clip(18, 95).astype(int)
    female = np.random.binomial(1, 0.52, n)
    comorbidity = np.random.binomial(1, expit(-2 + 0.03 * age), n)
    ses = np.random.choice([1, 2, 3], n, p=[0.3, 0.45, 0.25])
    healthcare_seeking = np.random.normal(0, 1, n) + 0.02 * age
    return pd.DataFrame({
        'age': age, 'female': female, 'comorbidity': comorbidity,
        'ses': ses, 'healthcare_seeking': healthcare_seeking
    })


# ============================================================
# 1. TEST-NEGATIVE DESIGN (TND)
# ============================================================

def run_tnd_analysis():
    """Simulate TND study and estimate VE via logistic regression."""
    print("=" * 60)
    print("1. Test-Negative Design (TND) Analysis")
    print("=" * 60)

    n = 8000
    pop = generate_baseline_population(n)

    # Vaccination probability (confounded by age, SES, healthcare-seeking)
    logit_vax = -0.5 + 0.01 * pop['age'] + 0.3 * pop['healthcare_seeking'] + 0.2 * (pop['ses'] - 2)
    pop['vaccinated'] = np.random.binomial(1, expit(logit_vax))

    # True VE = 0.70 (OR scale)
    true_ve = 0.70
    log_or_vax = np.log(1 - true_ve)

    # Infection probability (test-positive)
    logit_inf = -2.0 + 0.02 * pop['age'] + 0.5 * pop['comorbidity'] + log_or_vax * pop['vaccinated']
    pop['test_positive'] = np.random.binomial(1, expit(logit_inf))

    # Healthcare-seeking => all present for testing (TND assumption)
    pop['tested'] = 1

    # Unadjusted logistic regression
    model_unadj = smf.logit('test_positive ~ vaccinated', data=pop).fit(disp=0)
    or_unadj = np.exp(model_unadj.params['vaccinated'])
    ve_unadj = 1 - or_unadj

    # Adjusted logistic regression
    model_adj = smf.logit('test_positive ~ vaccinated + age + female + comorbidity + C(ses)', data=pop).fit(disp=0)
    or_adj = np.exp(model_adj.params['vaccinated'])
    ci_adj = np.exp(model_adj.conf_int().loc['vaccinated'])
    ve_adj = 1 - or_adj
    ve_ci = [1 - ci_adj.iloc[1], 1 - ci_adj.iloc[0]]

    print(f"  True VE: {true_ve:.1%}")
    print(f"  Unadjusted VE: {ve_unadj:.1%}")
    print(f"  Adjusted VE: {ve_adj:.1%} (95% CI: {ve_ci[0]:.1%}–{ve_ci[1]:.1%})")

    results['tnd'] = {
        'true_ve': true_ve, 'unadj_ve': ve_unadj,
        'adj_ve': ve_adj, 'adj_ci': ve_ci
    }

    # Sensitivity analysis: varying equi-confounding violations
    bias_factors = np.linspace(0, 0.5, 6)
    ve_biased = []
    for bf in bias_factors:
        logit_inf_b = logit_inf + bf * pop['healthcare_seeking'] * pop['vaccinated']
        pop_b = pop.copy()
        pop_b['test_positive'] = np.random.binomial(1, expit(logit_inf_b))
        m = smf.logit('test_positive ~ vaccinated + age + female + comorbidity + C(ses)', data=pop_b).fit(disp=0)
        ve_biased.append(1 - np.exp(m.params['vaccinated']))

    # Figure 1: TND VE estimates
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    categories = ['True VE', 'Unadjusted', 'Adjusted']
    vals = [true_ve, ve_unadj, ve_adj]
    colors = ['#2ecc71', '#e74c3c', '#3498db']
    bars = ax.bar(categories, vals, color=colors, edgecolor='black', linewidth=0.8)
    ax.errorbar(2, ve_adj, yerr=[[ve_adj - ve_ci[0]], [ve_ci[1] - ve_adj]],
                fmt='none', color='black', capsize=5, linewidth=2)
    ax.set_ylabel('Vaccine Effectiveness')
    ax.set_title('TND Vaccine Effectiveness Estimates')
    ax.set_ylim(0, 1)
    ax.axhline(y=true_ve, color='gray', linestyle='--', alpha=0.5)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.02,
                f'{val:.1%}', ha='center', fontweight='bold')

    ax = axes[1]
    ax.plot(bias_factors, ve_biased, 'o-', color='#e74c3c', linewidth=2, markersize=8)
    ax.axhline(y=true_ve, color='#2ecc71', linestyle='--', linewidth=2, label=f'True VE = {true_ve:.0%}')
    ax.set_xlabel('Equi-confounding Violation Magnitude')
    ax.set_ylabel('Estimated VE')
    ax.set_title('Sensitivity: Equi-confounding Violations')
    ax.legend()

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig1_tnd_analysis.png', bbox_inches='tight')
    plt.close()
    print("  → Figure saved: figures/fig1_tnd_analysis.png\n")
    return pop


# ============================================================
# 2. WANING VE ESTIMATION
# ============================================================

def run_waning_analysis():
    """Estimate time-varying VE waning using Cox PH and spline models."""
    print("=" * 60)
    print("2. Waning Vaccine Effectiveness Estimation")
    print("=" * 60)

    n = 6000
    pop = generate_baseline_population(n)

    pop['vaccinated'] = np.random.binomial(1, 0.6, n)
    pop['days_since_vax'] = np.where(pop['vaccinated'] == 1,
                                      np.random.uniform(7, 300, n).astype(int), 0)

    # True waning function: VE(t) = 0.90 * exp(-0.005 * t)
    def true_ve_func(t):
        return 0.90 * np.exp(-0.005 * t)

    # Hazard with waning
    base_hazard = 0.002
    hazard = base_hazard * np.exp(
        0.02 * (pop['age'] - 55) / 10 +
        0.3 * pop['comorbidity'] +
        np.where(pop['vaccinated'] == 1,
                 np.log(1 - true_ve_func(pop['days_since_vax'])),
                 0)
    )

    # Simulate survival times
    pop['event_time'] = np.random.exponential(1 / hazard.clip(1e-6))
    max_followup = 365
    pop['observed_time'] = np.minimum(pop['event_time'], max_followup)
    pop['event'] = (pop['event_time'] <= max_followup).astype(int)

    # Time categories for waning
    bins = [0, 30, 90, 150, 210, 300]
    labels = ['0-30d', '31-90d', '91-150d', '151-210d', '211-300d']
    all_cats = labels + ['Unvaccinated']
    pop['time_cat'] = pd.cut(pop['days_since_vax'], bins=bins, labels=labels, right=True).cat.add_categories('Unvaccinated')
    pop.loc[pop['vaccinated'] == 0, 'time_cat'] = 'Unvaccinated'

    # Cox PH model with time categories
    vax_only = pop[pop['vaccinated'] == 1].copy()
    for label in labels:
        vax_only[f'tc_{label}'] = (vax_only['time_cat'] == label).astype(int)

    cph_data = pop.copy()
    for label in labels:
        cph_data[f'tc_{label}'] = ((cph_data['time_cat'] == label)).astype(int)

    cph = CoxPHFitter()
    covariates = ['age', 'comorbidity'] + [f'tc_{l}' for l in labels]
    cph.fit(cph_data[covariates + ['observed_time', 'event']],
            duration_col='observed_time', event_col='event')

    # Extract VE by time period
    ve_by_period = {}
    midpoints = [15, 60, 120, 180, 255]
    for label, mid in zip(labels, midpoints):
        hr = np.exp(cph.params_[f'tc_{label}'])
        ve_est = 1 - hr
        ve_by_period[label] = {'ve': ve_est, 'midpoint': mid, 'true_ve': true_ve_func(mid)}

    print("  VE by time since vaccination:")
    for label, v in ve_by_period.items():
        print(f"    {label}: Estimated VE = {v['ve']:.1%}, True VE = {v['true_ve']:.1%}")

    results['waning'] = ve_by_period

    # Figure 2: Waning VE
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    t_range = np.linspace(0, 300, 200)
    ax.plot(t_range, true_ve_func(t_range), 'g-', linewidth=2.5, label='True VE(t)')
    mids = [v['midpoint'] for v in ve_by_period.values()]
    ves_est = [v['ve'] for v in ve_by_period.values()]
    ax.plot(mids, ves_est, 'rs-', markersize=10, linewidth=2, label='Estimated VE (Cox PH)')
    ax.set_xlabel('Days Since Vaccination')
    ax.set_ylabel('Vaccine Effectiveness')
    ax.set_title('VE Waning Over Time')
    ax.legend(fontsize=10)
    ax.set_ylim(0, 1)
    ax.set_xlim(0, 310)

    # KM curves
    ax = axes[1]
    kmf = KaplanMeierFitter()
    for cat_label, color in [('Unvaccinated', '#e74c3c'), ('0-30d', '#2ecc71'),
                              ('91-150d', '#f39c12'), ('211-300d', '#9b59b6')]:
        subset = cph_data[cph_data['time_cat'] == cat_label]
        if len(subset) > 10:
            kmf.fit(subset['observed_time'], subset['event'], label=cat_label)
            kmf.plot_survival_function(ax=ax, color=color, linewidth=2)
    ax.set_title('Kaplan-Meier Survival by Vaccination Timing')
    ax.set_xlabel('Days')
    ax.set_ylabel('Survival Probability')

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig2_waning_ve.png', bbox_inches='tight')
    plt.close()
    print("  → Figure saved: figures/fig2_waning_ve.png\n")


# ============================================================
# 3. VARIANT-SPECIFIC VE ESTIMATION
# ============================================================

def run_variant_analysis():
    """Estimate variant-specific VE using multinomial/logistic regression."""
    print("=" * 60)
    print("3. Variant-Specific VE Estimation")
    print("=" * 60)

    n = 10000
    pop = generate_baseline_population(n)
    pop['vaccinated'] = np.random.binomial(1, 0.55, n)

    # True VE by variant
    true_ve_delta = 0.75
    true_ve_omicron = 0.50
    true_ve_wt = 0.85

    # Variant proportions over time
    pop['study_week'] = np.random.randint(1, 53, n)
    delta_prop = expit(-3 + 0.15 * (pop['study_week'] - 20))
    omicron_prop = expit(-5 + 0.2 * (pop['study_week'] - 35))

    # Assign variants (among positives)
    logit_pos_wt = -2.5 + 0.02 * pop['age'] + np.log(1 - true_ve_wt) * pop['vaccinated']
    logit_pos_delta = -2.0 + 0.02 * pop['age'] + np.log(1 - true_ve_delta) * pop['vaccinated']
    logit_pos_omicron = -1.5 + 0.02 * pop['age'] + np.log(1 - true_ve_omicron) * pop['vaccinated']

    p_wt = expit(logit_pos_wt) * (1 - delta_prop) * (1 - omicron_prop)
    p_delta = expit(logit_pos_delta) * delta_prop * (1 - omicron_prop)
    p_omicron = expit(logit_pos_omicron) * omicron_prop

    pop['positive_wt'] = np.random.binomial(1, p_wt.clip(0, 1))
    pop['positive_delta'] = np.random.binomial(1, p_delta.clip(0, 1))
    pop['positive_omicron'] = np.random.binomial(1, p_omicron.clip(0, 1))
    pop['test_negative'] = ((pop['positive_wt'] == 0) &
                             (pop['positive_delta'] == 0) &
                             (pop['positive_omicron'] == 0)).astype(int)

    variant_results = {}
    for variant, col, true_ve in [('Wild-type', 'positive_wt', true_ve_wt),
                                    ('Delta', 'positive_delta', true_ve_delta),
                                    ('Omicron', 'positive_omicron', true_ve_omicron)]:
        cases = pop[pop[col] == 1].copy()
        controls = pop[pop['test_negative'] == 1].copy()
        tnd_data = pd.concat([cases.assign(outcome=1), controls.assign(outcome=0)])

        if len(tnd_data) > 50:
            model = smf.logit('outcome ~ vaccinated + age + comorbidity + C(ses)', data=tnd_data).fit(disp=0)
            or_est = np.exp(model.params['vaccinated'])
            ci = np.exp(model.conf_int().loc['vaccinated'])
            ve_est = 1 - or_est
            ve_ci = [1 - ci.iloc[1], 1 - ci.iloc[0]]
            variant_results[variant] = {
                'true_ve': true_ve, 'est_ve': ve_est, 'ci': ve_ci
            }
            print(f"  {variant}: True VE={true_ve:.0%}, Estimated VE={ve_est:.1%} "
                  f"(95% CI: {ve_ci[0]:.1%}–{ve_ci[1]:.1%})")

    results['variant'] = variant_results

    # Figure 3: Variant-specific VE
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    variants = list(variant_results.keys())
    x = np.arange(len(variants))
    true_vals = [variant_results[v]['true_ve'] for v in variants]
    est_vals = [variant_results[v]['est_ve'] for v in variants]
    ci_low = [variant_results[v]['ci'][0] for v in variants]
    ci_high = [variant_results[v]['ci'][1] for v in variants]
    yerr_low = [e - l for e, l in zip(est_vals, ci_low)]
    yerr_high = [h - e for e, h in zip(est_vals, ci_high)]

    ax.bar(x - 0.15, true_vals, 0.3, label='True VE', color='#2ecc71', edgecolor='black')
    ax.bar(x + 0.15, est_vals, 0.3, label='Estimated VE', color='#3498db', edgecolor='black')
    ax.errorbar(x + 0.15, est_vals, yerr=[yerr_low, yerr_high],
                fmt='none', color='black', capsize=5, linewidth=2)
    ax.set_xticks(x)
    ax.set_xticklabels(variants)
    ax.set_ylabel('Vaccine Effectiveness')
    ax.set_title('Variant-Specific VE Estimates (TND)')
    ax.legend()
    ax.set_ylim(0, 1)

    # Variant distribution over time
    ax = axes[1]
    weeks = np.arange(1, 53)
    d_prop = expit(-3 + 0.15 * (weeks - 20))
    o_prop = expit(-5 + 0.2 * (weeks - 35))
    w_prop = (1 - d_prop) * (1 - o_prop)
    d_adj = d_prop * (1 - o_prop)

    ax.stackplot(weeks, w_prop, d_adj, o_prop,
                 labels=['Wild-type', 'Delta', 'Omicron'],
                 colors=['#2ecc71', '#e74c3c', '#9b59b6'], alpha=0.8)
    ax.set_xlabel('Study Week')
    ax.set_ylabel('Proportion')
    ax.set_title('Variant Distribution Over Time')
    ax.legend(loc='center left')

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig3_variant_ve.png', bbox_inches='tight')
    plt.close()
    print("  → Figure saved: figures/fig3_variant_ve.png\n")


# ============================================================
# 4. HEALTHY VACCINEE BIAS CORRECTION
# ============================================================

def run_bias_correction():
    """Correct healthy vaccinee bias using IPW and negative controls."""
    print("=" * 60)
    print("4. Healthy Vaccinee Bias Correction")
    print("=" * 60)

    n = 8000
    pop = generate_baseline_population(n)

    # Unmeasured health status (confounder)
    pop['health_status'] = np.random.normal(0, 1, n) + 0.3 * pop['healthcare_seeking']

    # Vaccination strongly related to health status (healthy vaccinee effect)
    logit_vax = -0.3 + 0.5 * pop['health_status'] + 0.01 * pop['age']
    pop['vaccinated'] = np.random.binomial(1, expit(logit_vax))

    true_ve = 0.60
    # Outcome (infection): healthier people also less likely infected
    logit_outcome = (-2.0 + 0.02 * pop['age'] + 0.3 * pop['comorbidity']
                     - 0.4 * pop['health_status']  # healthier => less infection
                     + np.log(1 - true_ve) * pop['vaccinated'])
    pop['infected'] = np.random.binomial(1, expit(logit_outcome))

    # Negative control outcome (not affected by vaccine)
    logit_nc = -3.0 + 0.01 * pop['age'] - 0.3 * pop['health_status']
    pop['neg_control_outcome'] = np.random.binomial(1, expit(logit_nc))

    # Naive estimate (biased)
    m_naive = smf.logit('infected ~ vaccinated', data=pop).fit(disp=0)
    ve_naive = 1 - np.exp(m_naive.params['vaccinated'])

    # Adjusted (measured confounders only)
    m_adj = smf.logit('infected ~ vaccinated + age + comorbidity + C(ses)', data=pop).fit(disp=0)
    ve_adj = 1 - np.exp(m_adj.params['vaccinated'])

    # IPW correction
    ps_model = smf.logit('vaccinated ~ age + comorbidity + C(ses) + healthcare_seeking', data=pop).fit(disp=0)
    pop['ps'] = ps_model.predict(pop)
    pop['ipw'] = np.where(pop['vaccinated'] == 1, 1 / pop['ps'], 1 / (1 - pop['ps']))
    pop['ipw'] = pop['ipw'].clip(upper=np.percentile(pop['ipw'], 99))

    m_ipw = smf.glm('infected ~ vaccinated + age + comorbidity',
                     data=pop, family=sm.families.Binomial(),
                     freq_weights=pop['ipw']).fit()
    ve_ipw = 1 - np.exp(m_ipw.params['vaccinated'])

    # Negative control test
    m_nc = smf.logit('neg_control_outcome ~ vaccinated + age + comorbidity', data=pop).fit(disp=0)
    nc_or = np.exp(m_nc.params['vaccinated'])
    nc_pval = m_nc.pvalues['vaccinated']

    # Bias-corrected estimate
    ve_bias_corrected = 1 - np.exp(m_adj.params['vaccinated']) / nc_or

    print(f"  True VE: {true_ve:.1%}")
    print(f"  Naive VE: {ve_naive:.1%}")
    print(f"  Adjusted VE: {ve_adj:.1%}")
    print(f"  IPW VE: {ve_ipw:.1%}")
    print(f"  Negative control OR: {nc_or:.3f} (p={nc_pval:.4f})")
    print(f"  Bias-corrected VE: {ve_bias_corrected:.1%}")

    results['bias'] = {
        'true_ve': true_ve, 'naive_ve': ve_naive, 'adj_ve': ve_adj,
        'ipw_ve': ve_ipw, 'nc_or': nc_or, 'bias_corrected_ve': ve_bias_corrected
    }

    # Figure 4: Bias correction comparison
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    methods = ['True VE', 'Naive', 'Adjusted', 'IPW', 'Bias-\nCorrected']
    vals = [true_ve, ve_naive, ve_adj, ve_ipw, ve_bias_corrected]
    colors = ['#2ecc71', '#e74c3c', '#f39c12', '#3498db', '#9b59b6']
    bars = ax.bar(methods, vals, color=colors, edgecolor='black', linewidth=0.8)
    ax.axhline(y=true_ve, color='gray', linestyle='--', alpha=0.7)
    ax.set_ylabel('Vaccine Effectiveness')
    ax.set_title('Healthy Vaccinee Bias Correction Methods')
    ax.set_ylim(0, 1)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.02,
                f'{val:.1%}', ha='center', fontweight='bold', fontsize=9)

    # Propensity score distributions
    ax = axes[1]
    ax.hist(pop[pop['vaccinated'] == 1]['ps'], bins=40, alpha=0.6,
            color='#3498db', label='Vaccinated', density=True)
    ax.hist(pop[pop['vaccinated'] == 0]['ps'], bins=40, alpha=0.6,
            color='#e74c3c', label='Unvaccinated', density=True)
    ax.set_xlabel('Propensity Score')
    ax.set_ylabel('Density')
    ax.set_title('Propensity Score Distribution')
    ax.legend()

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig4_bias_correction.png', bbox_inches='tight')
    plt.close()
    print("  → Figure saved: figures/fig4_bias_correction.png\n")


# ============================================================
# 5. BOOSTER DOSE CAUSAL ESTIMATION (MSM)
# ============================================================

def run_booster_analysis():
    """Estimate booster effect using marginal structural model."""
    print("=" * 60)
    print("5. Booster Dose Causal Effect Estimation (MSM)")
    print("=" * 60)

    n = 7000
    pop = generate_baseline_population(n)

    # Primary vaccination
    pop['primary_vax'] = np.random.binomial(1, 0.7, n)

    # Time-varying confounders at booster decision
    pop['risk_perception'] = np.random.normal(0, 1, n) + 0.3 * pop['comorbidity']
    pop['prior_infection'] = np.random.binomial(1, 0.15, n)

    # Booster uptake (conditional on primary)
    logit_boost = (-1.0 + 0.5 * pop['primary_vax'] + 0.3 * pop['risk_perception']
                   + 0.4 * pop['comorbidity'] + 0.01 * pop['age'] - 0.3 * pop['prior_infection'])
    pop['booster'] = np.random.binomial(1, expit(logit_boost)) * pop['primary_vax']

    # True causal effects
    true_ve_primary = 0.55  # waned primary
    true_ve_booster = 0.80  # after booster

    # Outcome: severe COVID
    logit_severe = (-3.0 + 0.03 * pop['age'] + 0.5 * pop['comorbidity']
                    - 0.3 * pop['risk_perception']
                    + np.where(pop['booster'] == 1, np.log(1 - true_ve_booster),
                               np.where(pop['primary_vax'] == 1, np.log(1 - true_ve_primary), 0)))
    pop['severe_covid'] = np.random.binomial(1, expit(logit_severe))

    # Naive estimates
    m_naive = smf.logit('severe_covid ~ booster + primary_vax', data=pop).fit(disp=0)
    ve_booster_naive = 1 - np.exp(m_naive.params['booster'])

    # MSM: Inverse probability of treatment weighting
    # Step 1: Model booster receipt
    ps_boost = smf.logit('booster ~ primary_vax + age + comorbidity + risk_perception + prior_infection + C(ses)',
                         data=pop).fit(disp=0)
    pop['ps_boost'] = ps_boost.predict(pop)

    # Step 2: Stabilized weights
    p_marginal = pop['booster'].mean()
    pop['sw'] = np.where(pop['booster'] == 1,
                         p_marginal / pop['ps_boost'],
                         (1 - p_marginal) / (1 - pop['ps_boost']))
    pop['sw'] = pop['sw'].clip(upper=np.percentile(pop['sw'], 99))

    # Step 3: Weighted outcome model
    m_msm = smf.glm('severe_covid ~ booster + primary_vax + age',
                     data=pop, family=sm.families.Binomial(),
                     freq_weights=pop['sw']).fit()
    ve_booster_msm = 1 - np.exp(m_msm.params['booster'])
    ci_msm = np.exp(m_msm.conf_int().loc['booster'])
    ve_boost_ci = [1 - ci_msm.iloc[1], 1 - ci_msm.iloc[0]]

    # Adjusted regression
    m_adj = smf.logit('severe_covid ~ booster + primary_vax + age + comorbidity + risk_perception + C(ses)',
                      data=pop).fit(disp=0)
    ve_booster_adj = 1 - np.exp(m_adj.params['booster'])

    print(f"  True booster VE (incremental): {true_ve_booster:.1%}")
    print(f"  Naive booster VE: {ve_booster_naive:.1%}")
    print(f"  Adjusted booster VE: {ve_booster_adj:.1%}")
    print(f"  MSM booster VE: {ve_booster_msm:.1%} (95% CI: {ve_boost_ci[0]:.1%}–{ve_boost_ci[1]:.1%})")

    results['booster'] = {
        'true_ve': true_ve_booster, 'naive_ve': ve_booster_naive,
        'adj_ve': ve_booster_adj, 'msm_ve': ve_booster_msm, 'msm_ci': ve_boost_ci
    }

    # Figure 5: Booster effect estimation
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    methods = ['True VE', 'Naive', 'Adjusted', 'MSM']
    vals = [true_ve_booster, ve_booster_naive, ve_booster_adj, ve_booster_msm]
    colors = ['#2ecc71', '#e74c3c', '#f39c12', '#3498db']
    bars = ax.bar(methods, vals, color=colors, edgecolor='black', linewidth=0.8)
    ax.axhline(y=true_ve_booster, color='gray', linestyle='--', alpha=0.7)
    ax.set_ylabel('Vaccine Effectiveness (Booster)')
    ax.set_title('Booster VE: Causal Estimation Methods')
    ax.set_ylim(0, 1)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.02,
                f'{val:.1%}', ha='center', fontweight='bold', fontsize=9)

    # Stabilized weight distribution
    ax = axes[1]
    ax.hist(pop[pop['booster'] == 1]['sw'], bins=50, alpha=0.6, color='#3498db',
            label='Boosted', density=True)
    ax.hist(pop[pop['booster'] == 0]['sw'], bins=50, alpha=0.6, color='#e74c3c',
            label='Not Boosted', density=True)
    ax.set_xlabel('Stabilized Weight')
    ax.set_ylabel('Density')
    ax.set_title('MSM Stabilized Weight Distribution')
    ax.legend()
    ax.set_xlim(0, 5)

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig5_booster_msm.png', bbox_inches='tight')
    plt.close()
    print("  → Figure saved: figures/fig5_booster_msm.png\n")


# ============================================================
# 6. mRNA VACCINE HOSPITALIZATION CASE STUDY
# ============================================================

def run_hospitalization_case_study():
    """Case study: mRNA vaccine effectiveness against hospitalization."""
    print("=" * 60)
    print("6. mRNA Vaccine Hospitalization Prevention Case Study")
    print("=" * 60)

    n = 12000
    pop = generate_baseline_population(n)

    # Vaccination schedule
    pop['dose_0'] = np.random.binomial(1, 0.3, n)  # unvaccinated
    pop['dose_1'] = np.random.binomial(1, 0.15, n) * (1 - pop['dose_0'])
    pop['dose_2'] = np.random.binomial(1, 0.35, n) * (1 - pop['dose_0']) * (1 - pop['dose_1'])
    pop['dose_3'] = (1 - pop['dose_0'] - pop['dose_1'] - pop['dose_2']).clip(0, 1)
    pop['doses'] = pop['dose_1'] + 2 * pop['dose_2'] + 3 * pop['dose_3']

    # Time since last dose
    pop['days_since_dose'] = np.where(pop['doses'] > 0,
                                       np.random.uniform(14, 250, n).astype(int), 0)

    # True VE against hospitalization (by dose and time)
    def true_hosp_ve(doses, days):
        base_ve = np.where(doses == 0, 0,
                  np.where(doses == 1, 0.65,
                  np.where(doses == 2, 0.85, 0.92)))
        waning = np.exp(-0.003 * days)
        return base_ve * waning

    pop['true_ve'] = true_hosp_ve(pop['doses'], pop['days_since_dose'])

    # Hospitalization
    base_hosp_rate = 0.015
    hazard = base_hosp_rate * np.exp(
        0.04 * (pop['age'] - 55) / 10 +
        0.8 * pop['comorbidity'] +
        np.log((1 - pop['true_ve']).clip(0.01, 1))
    )
    pop['hosp_time'] = np.random.exponential(1 / hazard.clip(1e-6))
    max_fu = 180
    pop['time'] = np.minimum(pop['hosp_time'], max_fu)
    pop['hospitalized'] = (pop['hosp_time'] <= max_fu).astype(int)

    # Age group analysis
    pop['age_group'] = pd.cut(pop['age'], bins=[17, 49, 64, 79, 96],
                               labels=['18-49', '50-64', '65-79', '80+'])

    # Cox PH model
    cph_data = pop[['time', 'hospitalized', 'doses', 'age', 'female',
                     'comorbidity', 'days_since_dose']].copy()
    cph_data['dose_1plus'] = (cph_data['doses'] >= 1).astype(int)
    cph_data['dose_2plus'] = (cph_data['doses'] >= 2).astype(int)
    cph_data['dose_3'] = (cph_data['doses'] == 3).astype(int)

    cph = CoxPHFitter()
    cph.fit(cph_data[['time', 'hospitalized', 'dose_1plus', 'dose_2plus', 'dose_3',
                       'age', 'comorbidity', 'days_since_dose']],
            duration_col='time', event_col='hospitalized')

    ve_dose_results = {}
    for dose_col in ['dose_1plus', 'dose_2plus', 'dose_3']:
        hr = np.exp(cph.params_[dose_col])
        ve = 1 - hr
        ci = np.exp(cph.confidence_intervals_.loc[dose_col])
        ve_ci = [1 - ci.values.flatten()[1], 1 - ci.values.flatten()[0]]
        ve_dose_results[dose_col] = {'ve': ve, 'ci': ve_ci}

    # VE by age group (logistic regression per group)
    ve_age = {}
    for ag in ['18-49', '50-64', '65-79', '80+']:
        sub = pop[pop['age_group'] == ag].copy()
        sub['any_vax'] = (sub['doses'] > 0).astype(int)
        if sub['hospitalized'].sum() > 10:
            m = smf.logit('hospitalized ~ any_vax + comorbidity', data=sub).fit(disp=0)
            ve = 1 - np.exp(m.params['any_vax'])
            ci = np.exp(m.conf_int().loc['any_vax'])
            ve_age[ag] = {'ve': ve, 'ci': [1 - ci.iloc[1], 1 - ci.iloc[0]]}

    print("  VE against hospitalization by dose:")
    for d, v in ve_dose_results.items():
        print(f"    {d}: VE = {v['ve']:.1%} (95% CI: {v['ci'][0]:.1%}–{v['ci'][1]:.1%})")
    print("  VE by age group:")
    for ag, v in ve_age.items():
        print(f"    {ag}: VE = {v['ve']:.1%} (95% CI: {v['ci'][0]:.1%}–{v['ci'][1]:.1%})")

    results['hospitalization'] = {
        'dose_results': ve_dose_results,
        'age_results': ve_age
    }

    # Figure 6: Hospitalization case study
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # Panel A: VE by dose
    ax = axes[0, 0]
    dose_labels = ['≥1 dose', '≥2 doses', '3 doses (booster)']
    dose_ves = [ve_dose_results[k]['ve'] for k in ['dose_1plus', 'dose_2plus', 'dose_3']]
    dose_cis = [ve_dose_results[k]['ci'] for k in ['dose_1plus', 'dose_2plus', 'dose_3']]
    x = np.arange(len(dose_labels))
    bars = ax.bar(x, dose_ves, color=['#f39c12', '#3498db', '#2ecc71'], edgecolor='black')
    for i, (bar, val, ci) in enumerate(zip(bars, dose_ves, dose_cis)):
        ax.errorbar(i, val, yerr=[[val - ci[0]], [ci[1] - val]],
                    fmt='none', color='black', capsize=6, linewidth=2)
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.03,
                f'{val:.1%}', ha='center', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(dose_labels)
    ax.set_ylabel('VE Against Hospitalization')
    ax.set_title('A) VE by Vaccination Dose')
    ax.set_ylim(0, 1.1)

    # Panel B: VE by age group
    ax = axes[0, 1]
    age_groups = list(ve_age.keys())
    age_ves = [ve_age[ag]['ve'] for ag in age_groups]
    age_cis = [ve_age[ag]['ci'] for ag in age_groups]
    x = np.arange(len(age_groups))
    bars = ax.bar(x, age_ves, color='#3498db', edgecolor='black')
    for i, (bar, val, ci) in enumerate(zip(bars, age_ves, age_cis)):
        yerr_low = max(val - ci[0], 0)
        yerr_high = max(ci[1] - val, 0)
        ax.errorbar(i, val, yerr=[[yerr_low], [yerr_high]],
                    fmt='none', color='black', capsize=6, linewidth=2)
        ax.text(bar.get_x() + bar.get_width()/2., max(val, 0) + 0.03,
                f'{val:.1%}', ha='center', fontweight='bold', fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(age_groups)
    ax.set_ylabel('VE Against Hospitalization')
    ax.set_title('B) VE by Age Group')
    ax.set_ylim(0, 1.1)

    # Panel C: KM curves
    ax = axes[1, 0]
    kmf = KaplanMeierFitter()
    dose_colors = {0: '#e74c3c', 1: '#f39c12', 2: '#3498db', 3: '#2ecc71'}
    dose_names = {0: 'Unvaccinated', 1: '1 dose', 2: '2 doses', 3: '3 doses'}
    for d in [0, 1, 2, 3]:
        sub = pop[pop['doses'] == d]
        if len(sub) > 20:
            kmf.fit(sub['time'], sub['hospitalized'], label=dose_names[d])
            kmf.plot_survival_function(ax=ax, color=dose_colors[d], linewidth=2)
    ax.set_title('C) Kaplan-Meier: Hospitalization-Free Survival')
    ax.set_xlabel('Days')
    ax.set_ylabel('Survival Probability')

    # Panel D: Cumulative incidence
    ax = axes[1, 1]
    for d in [0, 1, 2, 3]:
        sub = pop[pop['doses'] == d]
        if len(sub) > 20:
            kmf.fit(sub['time'], sub['hospitalized'], label=dose_names[d])
            t = kmf.survival_function_.index
            ci_vals = 1 - kmf.survival_function_.values.flatten()
            ax.plot(t, ci_vals, color=dose_colors[d], linewidth=2, label=dose_names[d])
    ax.set_title('D) Cumulative Hospitalization Incidence')
    ax.set_xlabel('Days')
    ax.set_ylabel('Cumulative Incidence')
    ax.legend()

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig6_hospitalization.png', bbox_inches='tight')
    plt.close()
    print("  → Figure saved: figures/fig6_hospitalization.png\n")


# ============================================================
# 7. SUMMARY FOREST PLOT
# ============================================================

def create_summary_plot():
    """Create a summary forest plot of all VE estimates."""
    print("=" * 60)
    print("7. Creating Summary Forest Plot")
    print("=" * 60)

    fig, ax = plt.subplots(figsize=(10, 8))

    estimates = [
        ('TND (Adjusted)', results['tnd']['adj_ve'], results['tnd']['adj_ci']),
        ('Waning: 0-30d', results['waning']['0-30d']['ve'], None),
        ('Waning: 91-150d', results['waning']['91-150d']['ve'], None),
        ('Waning: 211-300d', results['waning']['211-300d']['ve'], None),
        ('Variant: Wild-type', results['variant']['Wild-type']['est_ve'], results['variant']['Wild-type']['ci']),
        ('Variant: Delta', results['variant']['Delta']['est_ve'], results['variant']['Delta']['ci']),
        ('Variant: Omicron', results['variant']['Omicron']['est_ve'], results['variant']['Omicron']['ci']),
        ('Bias-Corrected', results['bias']['bias_corrected_ve'], None),
        ('Booster (MSM)', results['booster']['msm_ve'], results['booster']['msm_ci']),
    ]

    y_pos = np.arange(len(estimates))
    for i, (label, ve, ci) in enumerate(estimates):
        color = '#3498db'
        ax.plot(ve, i, 'D', color=color, markersize=10, zorder=5)
        if ci:
            ax.plot(ci, [i, i], '-', color=color, linewidth=2, zorder=4)
            ax.plot(ci, [i, i], '|', color=color, markersize=10, zorder=4)
        ax.text(max(ve, 0) + 0.02, i, f'{ve:.1%}', va='center', fontsize=9)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([e[0] for e in estimates])
    ax.set_xlabel('Vaccine Effectiveness')
    ax.set_title('Summary: VE Estimates Across All Analyses')
    ax.axvline(x=0, color='black', linewidth=0.5)
    ax.set_xlim(-0.1, 1.1)
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig7_summary_forest.png', bbox_inches='tight')
    plt.close()
    print("  → Figure saved: figures/fig7_summary_forest.png\n")


# ============================================================
# MAIN EXECUTION
# ============================================================

if __name__ == '__main__':
    print("=" * 60)
    print("  VACCINE EFFECTIVENESS ESTIMATION FRAMEWORK")
    print("  Comprehensive Methodological Analysis Pipeline")
    print("=" * 60)
    print()

    run_tnd_analysis()
    run_waning_analysis()
    run_variant_analysis()
    run_bias_correction()
    run_booster_analysis()
    run_hospitalization_case_study()
    create_summary_plot()

    print("=" * 60)
    print("  ALL ANALYSES COMPLETE")
    print("=" * 60)
    print(f"\nGenerated figures in '{FIGDIR}/':")
    print("  fig1_tnd_analysis.png")
    print("  fig2_waning_ve.png")
    print("  fig3_variant_ve.png")
    print("  fig4_bias_correction.png")
    print("  fig5_booster_msm.png")
    print("  fig6_hospitalization.png")
    print("  fig7_summary_forest.png")
