"""
Systematic Comparison of Causal Inference Methods for Observational Data
=========================================================================
Implements: PSM, IV, DID, DML, Causal Forest with pharmacoepidemiology case study.
Uses DoWhy/EconML framework.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import warnings
warnings.filterwarnings('ignore')

from sklearn.linear_model import LogisticRegression, LinearRegression, Lasso
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier, RandomForestRegressor
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_predict

import dowhy
from dowhy import CausalModel

from econml.dml import LinearDML, CausalForestDML
try:
    from econml.cate_interpreter import SingleTreeCateInterpreter
except ImportError:
    SingleTreeCateInterpreter = None

np.random.seed(42)

# ============================================================
# 1. DATA GENERATION: Pharmacoepidemiology Simulation
# ============================================================

def generate_pharma_data(n=5000, true_ate=-2.5, heterogeneous=True):
    """
    Simulate RWD for a cardiovascular drug study.
    Outcome: blood pressure reduction (continuous).
    Treatment: new antihypertensive drug vs standard care.
    Confounders: age, BMI, baseline BP, comorbidity score, smoking status.
    """
    age = np.random.normal(60, 12, n).clip(30, 90)
    bmi = np.random.normal(27, 5, n).clip(18, 45)
    baseline_bp = np.random.normal(150, 20, n).clip(110, 200)
    comorbidity = np.random.poisson(2, n).clip(0, 8)
    smoking = np.random.binomial(1, 0.3, n)

    # Treatment assignment (confounded)
    propensity_logit = (
        -2.0
        + 0.01 * age
        + 0.02 * bmi
        + 0.005 * baseline_bp
        - 0.15 * comorbidity
        + 0.4 * smoking
    )
    propensity = 1 / (1 + np.exp(-propensity_logit))
    treatment = np.random.binomial(1, propensity, n)

    # Heterogeneous treatment effect
    if heterogeneous:
        tau = true_ate - 0.05 * (age - 60) + 0.1 * (bmi - 27) - 0.3 * comorbidity
    else:
        tau = np.full(n, true_ate)

    # Outcome: post-treatment blood pressure
    noise = np.random.normal(0, 5, n)
    outcome = (
        baseline_bp
        - 0.1 * age
        + 0.3 * bmi
        + 0.5 * comorbidity
        + 2.0 * smoking
        + tau * treatment
        + noise
    )

    # Instrument: physician preference (correlated with treatment, not outcome)
    physician_pref = 0.6 * treatment + np.random.normal(0, 0.5, n)
    physician_pref = (physician_pref > 0.5).astype(int)

    df = pd.DataFrame({
        'age': age, 'bmi': bmi, 'baseline_bp': baseline_bp,
        'comorbidity': comorbidity, 'smoking': smoking,
        'treatment': treatment, 'outcome': outcome,
        'physician_pref': physician_pref, 'true_effect': tau,
        'propensity': propensity
    })
    return df


def generate_did_data(n_units=500, n_periods=10, treat_period=5, true_effect=-3.0):
    """Generate panel data for DID analysis."""
    units = np.repeat(np.arange(n_units), n_periods)
    periods = np.tile(np.arange(n_periods), n_units)
    treated_group = np.repeat(np.random.binomial(1, 0.5, n_units), n_periods)
    post = (periods >= treat_period).astype(int)

    unit_fe = np.repeat(np.random.normal(0, 2, n_units), n_periods)
    time_fe = np.tile(np.linspace(0, 5, n_periods), n_units)
    noise = np.random.normal(0, 1, n_units * n_periods)

    # Parallel trends hold by construction
    outcome = 10 + unit_fe + time_fe + true_effect * treated_group * post + noise

    return pd.DataFrame({
        'unit': units, 'period': periods, 'treated_group': treated_group,
        'post': post, 'outcome': outcome
    })


# ============================================================
# 2. CAUSAL INFERENCE METHODS
# ============================================================

def propensity_score_matching(df, n_neighbors=1):
    """PSM with nearest-neighbor matching."""
    covariates = ['age', 'bmi', 'baseline_bp', 'comorbidity', 'smoking']
    X = df[covariates].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Estimate propensity scores
    ps_model = LogisticRegression(max_iter=1000, C=1.0)
    ps_model.fit(X_scaled, df['treatment'].values)
    ps = ps_model.predict_proba(X_scaled)[:, 1]

    treated_idx = np.where(df['treatment'] == 1)[0]
    control_idx = np.where(df['treatment'] == 0)[0]

    # Match on propensity scores
    nn = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean')
    nn.fit(ps[control_idx].reshape(-1, 1))
    distances, indices = nn.kneighbors(ps[treated_idx].reshape(-1, 1))

    matched_control_idx = control_idx[indices.flatten()]
    matched_treated_idx = np.repeat(treated_idx, n_neighbors)

    ate = df['outcome'].iloc[matched_treated_idx].values.mean() - \
          df['outcome'].iloc[matched_control_idx].values.mean()

    # Bootstrap SE
    boot_ates = []
    for _ in range(200):
        boot_idx = np.random.choice(len(matched_treated_idx), len(matched_treated_idx), replace=True)
        boot_ate = df['outcome'].iloc[matched_treated_idx[boot_idx]].values.mean() - \
                   df['outcome'].iloc[matched_control_idx[boot_idx]].values.mean()
        boot_ates.append(boot_ate)
    se = np.std(boot_ates)

    return {'method': 'PSM', 'ate': ate, 'se': se, 'ps': ps,
            'matched_treated': matched_treated_idx, 'matched_control': matched_control_idx}


def ipw_estimator(df):
    """Inverse Probability Weighting as PSM alternative."""
    covariates = ['age', 'bmi', 'baseline_bp', 'comorbidity', 'smoking']
    X = StandardScaler().fit_transform(df[covariates].values)
    ps_model = LogisticRegression(max_iter=1000, C=1.0)
    ps_model.fit(X, df['treatment'].values)
    ps = ps_model.predict_proba(X)[:, 1].clip(0.05, 0.95)

    T = df['treatment'].values
    Y = df['outcome'].values
    w1 = T / ps
    w0 = (1 - T) / (1 - ps)
    ate = (w1 * Y).sum() / w1.sum() - (w0 * Y).sum() / w0.sum()

    # Bootstrap SE
    boot_ates = []
    n = len(df)
    for _ in range(200):
        idx = np.random.choice(n, n, replace=True)
        w1b = T[idx] / ps[idx]
        w0b = (1 - T[idx]) / (1 - ps[idx])
        boot_ate = (w1b * Y[idx]).sum() / w1b.sum() - (w0b * Y[idx]).sum() / w0b.sum()
        boot_ates.append(boot_ate)
    se = np.std(boot_ates)

    return {'method': 'IPW', 'ate': ate, 'se': se}


def iv_estimation(df):
    """Two-stage least squares with physician preference as instrument."""
    covariates = ['age', 'bmi', 'baseline_bp', 'comorbidity', 'smoking']
    X = df[covariates].values
    Z = df['physician_pref'].values
    T = df['treatment'].values
    Y = df['outcome'].values

    # First stage
    X_first = np.column_stack([X, Z])
    first_stage = LinearRegression()
    first_stage.fit(X_first, T)
    T_hat = first_stage.predict(X_first)
    f_stat = np.var(T_hat) / np.var(T - T_hat) * (len(T) - X_first.shape[1] - 1)

    # Second stage
    X_second = np.column_stack([X, T_hat])
    second_stage = LinearRegression()
    second_stage.fit(X_second, Y)
    ate = second_stage.coef_[-1]

    # Bootstrap SE
    boot_ates = []
    n = len(df)
    for _ in range(200):
        idx = np.random.choice(n, n, replace=True)
        fs = LinearRegression().fit(X_first[idx], T[idx])
        th = fs.predict(X_first[idx])
        ss = LinearRegression().fit(np.column_stack([X[idx], th]), Y[idx])
        boot_ates.append(ss.coef_[-1])
    se = np.std(boot_ates)

    return {'method': 'IV-2SLS', 'ate': ate, 'se': se, 'f_stat': f_stat}


def iv_weak_instrument_analysis(df):
    """Analyze weak instrument problem with varying instrument strength."""
    covariates = ['age', 'bmi', 'baseline_bp', 'comorbidity', 'smoking']
    X = df[covariates].values
    T = df['treatment'].values
    Y = df['outcome'].values

    strengths = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
    results = []
    for s in strengths:
        Z = s * T + np.random.normal(0, 1, len(T))
        Z = (Z > 0.5).astype(int)
        X_first = np.column_stack([X, Z])
        fs = LinearRegression().fit(X_first, T)
        T_hat = fs.predict(X_first)
        f_stat = np.var(T_hat) / np.var(T - T_hat) * (len(T) - X_first.shape[1] - 1)
        X_second = np.column_stack([X, T_hat])
        ss = LinearRegression().fit(X_second, Y)
        results.append({'strength': s, 'ate': ss.coef_[-1], 'f_stat': f_stat})
    return pd.DataFrame(results)


def did_analysis(df_did, true_effect=-3.0):
    """Difference-in-Differences with parallel trends test."""
    # Standard DID
    did_model = LinearRegression()
    X_did = df_did[['treated_group', 'post']].values
    interaction = (df_did['treated_group'] * df_did['post']).values
    X_full = np.column_stack([X_did, interaction])
    did_model.fit(X_full, df_did['outcome'].values)
    ate = did_model.coef_[2]

    # Bootstrap SE
    units = df_did['unit'].unique()
    boot_ates = []
    for _ in range(200):
        boot_units = np.random.choice(units, len(units), replace=True)
        boot_df = pd.concat([df_did[df_did['unit'] == u] for u in boot_units])
        X_b = np.column_stack([boot_df[['treated_group', 'post']].values,
                               (boot_df['treated_group'] * boot_df['post']).values])
        m = LinearRegression().fit(X_b, boot_df['outcome'].values)
        boot_ates.append(m.coef_[2])
    se = np.std(boot_ates)

    # Parallel trends test (pre-treatment period)
    pre_df = df_did[df_did['post'] == 0].copy()
    pre_trends = pre_df.groupby(['period', 'treated_group'])['outcome'].mean().unstack()
    pre_trends.columns = ['control', 'treated']

    return {'method': 'DID', 'ate': ate, 'se': se, 'pre_trends': pre_trends}


def dml_estimation(df):
    """Double/Debiased Machine Learning using EconML."""
    covariates = ['age', 'bmi', 'baseline_bp', 'comorbidity', 'smoking']
    X = df[covariates].values
    T = df['treatment'].values.reshape(-1, 1)
    Y = df['outcome'].values

    dml = LinearDML(
        model_y=GradientBoostingRegressor(n_estimators=100, max_depth=3, random_state=42),
        model_t=GradientBoostingClassifier(n_estimators=100, max_depth=3, random_state=42),
        discrete_treatment=True,
        cv=5, random_state=42
    )
    dml.fit(Y, T, X=X)
    ate = dml.ate(X)
    ci = dml.ate_interval(X, alpha=0.05)
    cate = dml.effect(X)

    return {'method': 'DML', 'ate': ate, 'se': (ci[1] - ci[0]) / (2 * 1.96),
            'ci': ci, 'cate': cate, 'model': dml}


def causal_forest_estimation(df):
    """Causal Forest for heterogeneous treatment effects using EconML."""
    covariates = ['age', 'bmi', 'baseline_bp', 'comorbidity', 'smoking']
    X = df[covariates].values
    T = df['treatment'].values.reshape(-1, 1)
    Y = df['outcome'].values

    cf = CausalForestDML(
        model_y=GradientBoostingRegressor(n_estimators=100, max_depth=3, random_state=42),
        model_t=GradientBoostingClassifier(n_estimators=100, max_depth=3, random_state=42),
        discrete_treatment=True,
        n_estimators=200, min_samples_leaf=20,
        cv=5, random_state=42
    )
    cf.fit(Y, T, X=X)
    ate = cf.ate(X)
    ci = cf.ate_interval(X, alpha=0.05)
    cate = cf.effect(X)

    return {'method': 'Causal Forest', 'ate': ate, 'se': (ci[1] - ci[0]) / (2 * 1.96),
            'ci': ci, 'cate': cate, 'model': cf, 'X': X,
            'feature_names': covariates}


def dowhy_workflow(df):
    """DoWhy causal inference workflow with refutation tests."""
    model = CausalModel(
        data=df,
        treatment='treatment',
        outcome='outcome',
        common_causes=['age', 'bmi', 'baseline_bp', 'comorbidity', 'smoking'],
        instruments=['physician_pref']
    )
    identified = model.identify_effect(proceed_when_unidentifiable=True)

    # Backdoor estimation
    estimate_bd = model.estimate_effect(
        identified,
        method_name="backdoor.linear_regression"
    )

    # Refutation: placebo treatment
    refute_placebo = model.refute_estimate(
        identified, estimate_bd,
        method_name="placebo_treatment_refuter",
        placebo_type="permute", num_simulations=50
    )

    # Refutation: random common cause
    refute_rcc = model.refute_estimate(
        identified, estimate_bd,
        method_name="random_common_cause",
        num_simulations=50
    )

    # Refutation: data subset
    refute_subset = model.refute_estimate(
        identified, estimate_bd,
        method_name="data_subset_refuter",
        subset_fraction=0.8, num_simulations=50
    )

    return {
        'estimate': estimate_bd.value,
        'refute_placebo': refute_placebo,
        'refute_rcc': refute_rcc,
        'refute_subset': refute_subset
    }


# ============================================================
# 3. VISUALIZATION
# ============================================================

def plot_method_comparison(results, true_ate, save_path='figures/method_comparison.png'):
    """Bar chart comparing all methods' ATE estimates."""
    fig, ax = plt.subplots(figsize=(10, 6))
    methods = [r['method'] for r in results]
    ates = [r['ate'] for r in results]
    ses = [r['se'] for r in results]

    colors = ['#2196F3', '#4CAF50', '#FF9800', '#9C27B0', '#E91E63', '#00BCD4']
    bars = ax.bar(range(len(methods)), ates, yerr=[1.96*s for s in ses],
                  capsize=5, color=colors[:len(methods)], alpha=0.8, edgecolor='black')
    ax.axhline(y=true_ate, color='red', linestyle='--', linewidth=2, label=f'True ATE = {true_ate}')
    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(methods, rotation=30, ha='right', fontsize=11)
    ax.set_ylabel('Average Treatment Effect (ATE)', fontsize=12)
    ax.set_title('Comparison of Causal Inference Methods\n(with 95% CI)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(axis='y', alpha=0.3)

    for i, (a, s) in enumerate(zip(ates, ses)):
        ax.text(i, a + 1.96*s + 0.3, f'{a:.2f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_propensity_distribution(df, ps, save_path='figures/propensity_scores.png'):
    """Propensity score distribution by treatment group."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Before matching
    axes[0].hist(ps[df['treatment'] == 0], bins=40, alpha=0.6, label='Control', color='#2196F3', density=True)
    axes[0].hist(ps[df['treatment'] == 1], bins=40, alpha=0.6, label='Treated', color='#E91E63', density=True)
    axes[0].set_xlabel('Propensity Score', fontsize=12)
    axes[0].set_ylabel('Density', fontsize=12)
    axes[0].set_title('Propensity Score Distribution\n(Before Matching)', fontsize=13)
    axes[0].legend(fontsize=11)

    # True vs estimated propensity
    axes[1].scatter(df['propensity'], ps, alpha=0.3, s=10, c=df['treatment'].map({0: '#2196F3', 1: '#E91E63'}))
    axes[1].plot([0, 1], [0, 1], 'k--', alpha=0.5)
    axes[1].set_xlabel('True Propensity', fontsize=12)
    axes[1].set_ylabel('Estimated Propensity', fontsize=12)
    axes[1].set_title('True vs Estimated Propensity Scores', fontsize=13)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_weak_instrument(iv_weak_df, save_path='figures/weak_instrument.png'):
    """Weak instrument analysis visualization."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].plot(iv_weak_df['strength'], iv_weak_df['ate'], 'o-', color='#9C27B0', linewidth=2, markersize=8)
    axes[0].axhline(y=-2.5, color='red', linestyle='--', linewidth=2, label='True ATE = -2.5')
    axes[0].set_xlabel('Instrument Strength', fontsize=12)
    axes[0].set_ylabel('Estimated ATE', fontsize=12)
    axes[0].set_title('IV Estimate vs Instrument Strength', fontsize=13)
    axes[0].legend(fontsize=11)
    axes[0].grid(alpha=0.3)

    axes[1].bar(range(len(iv_weak_df)), iv_weak_df['f_stat'], color='#FF9800', alpha=0.8, edgecolor='black')
    axes[1].axhline(y=10, color='red', linestyle='--', linewidth=2, label='Stock-Yogo threshold (F=10)')
    axes[1].set_xticks(range(len(iv_weak_df)))
    axes[1].set_xticklabels([f'{s:.1f}' for s in iv_weak_df['strength']], fontsize=9)
    axes[1].set_xlabel('Instrument Strength', fontsize=12)
    axes[1].set_ylabel('First-Stage F-statistic', fontsize=12)
    axes[1].set_title('Weak Instrument Diagnostics', fontsize=13)
    axes[1].legend(fontsize=11)
    axes[1].grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_did_trends(pre_trends, did_ate, save_path='figures/did_parallel_trends.png'):
    """DID parallel trends visualization."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Pre-trends
    axes[0].plot(pre_trends.index, pre_trends['control'], 'o-', label='Control', color='#2196F3', linewidth=2)
    axes[0].plot(pre_trends.index, pre_trends['treated'], 's-', label='Treated', color='#E91E63', linewidth=2)
    axes[0].axvline(x=4.5, color='gray', linestyle=':', linewidth=2, alpha=0.7, label='Treatment onset')
    axes[0].set_xlabel('Period', fontsize=12)
    axes[0].set_ylabel('Mean Outcome', fontsize=12)
    axes[0].set_title('Parallel Trends Assessment\n(Pre-Treatment Periods)', fontsize=13)
    axes[0].legend(fontsize=11)
    axes[0].grid(alpha=0.3)

    # DID effect
    diff = pre_trends['treated'] - pre_trends['control']
    axes[1].plot(pre_trends.index, diff, 'D-', color='#4CAF50', linewidth=2, markersize=8)
    axes[1].axhline(y=diff.mean(), color='red', linestyle='--', alpha=0.7, label=f'Mean diff = {diff.mean():.2f}')
    axes[1].set_xlabel('Period', fontsize=12)
    axes[1].set_ylabel('Treated - Control Difference', fontsize=12)
    axes[1].set_title('Pre-Treatment Group Differences', fontsize=13)
    axes[1].legend(fontsize=11)
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_heterogeneous_effects(df, cate_dml, cate_cf, save_path='figures/heterogeneous_effects.png'):
    """Heterogeneous treatment effects visualization."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    covariates = ['age', 'bmi', 'baseline_bp', 'comorbidity', 'smoking']
    titles = ['Age', 'BMI', 'Baseline BP', 'Comorbidity Score', 'Smoking Status']

    for i, (cov, title) in enumerate(zip(covariates, titles)):
        row, col = divmod(i, 3)
        ax = axes[row, col]
        ax.scatter(df[cov], cate_cf.flatten(), alpha=0.2, s=8, color='#9C27B0', label='Causal Forest')
        ax.scatter(df[cov], df['true_effect'], alpha=0.1, s=8, color='red', label='True CATE')

        # Binned means
        bins = pd.qcut(df[cov], 10, duplicates='drop')
        binned_cf = pd.DataFrame({'x': df[cov], 'cate': cate_cf.flatten(), 'true': df['true_effect']}).groupby(bins).mean()
        ax.plot(binned_cf['x'], binned_cf['cate'], 'o-', color='#9C27B0', linewidth=2, markersize=6)
        ax.plot(binned_cf['x'], binned_cf['true'], 's--', color='red', linewidth=2, markersize=6)

        ax.set_xlabel(title, fontsize=11)
        ax.set_ylabel('CATE', fontsize=11)
        ax.set_title(f'CATE by {title}', fontsize=12)
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=9)

    # Last subplot: DML vs CF distribution
    ax = axes[1, 2]
    ax.hist(cate_dml.flatten(), bins=40, alpha=0.6, label='DML', color='#4CAF50', density=True)
    ax.hist(cate_cf.flatten(), bins=40, alpha=0.6, label='Causal Forest', color='#9C27B0', density=True)
    ax.hist(df['true_effect'], bins=40, alpha=0.4, label='True CATE', color='red', density=True)
    ax.set_xlabel('CATE', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title('CATE Distribution Comparison', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    plt.suptitle('Heterogeneous Treatment Effect Analysis', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_bias_rmse(results, true_ate, save_path='figures/bias_rmse.png'):
    """Bias and RMSE comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    methods = [r['method'] for r in results]
    biases = [r['ate'] - true_ate for r in results]
    abs_biases = [abs(b) for b in biases]
    rmses = [np.sqrt(b**2 + r['se']**2) for b, r in zip(biases, results)]
    colors = ['#2196F3', '#4CAF50', '#FF9800', '#9C27B0', '#E91E63', '#00BCD4']

    bars1 = axes[0].bar(range(len(methods)), abs_biases, color=colors[:len(methods)], alpha=0.8, edgecolor='black')
    axes[0].set_xticks(range(len(methods)))
    axes[0].set_xticklabels(methods, rotation=30, ha='right', fontsize=10)
    axes[0].set_ylabel('|Bias|', fontsize=12)
    axes[0].set_title('Absolute Bias by Method', fontsize=13, fontweight='bold')
    axes[0].grid(axis='y', alpha=0.3)
    for i, b in enumerate(abs_biases):
        axes[0].text(i, b + 0.02, f'{b:.3f}', ha='center', fontsize=9)

    bars2 = axes[1].bar(range(len(methods)), rmses, color=colors[:len(methods)], alpha=0.8, edgecolor='black')
    axes[1].set_xticks(range(len(methods)))
    axes[1].set_xticklabels(methods, rotation=30, ha='right', fontsize=10)
    axes[1].set_ylabel('RMSE', fontsize=12)
    axes[1].set_title('Root Mean Squared Error by Method', fontsize=13, fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3)
    for i, r in enumerate(rmses):
        axes[1].text(i, r + 0.02, f'{r:.3f}', ha='center', fontsize=9)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_dowhy_refutation(dowhy_results, save_path='figures/dowhy_refutation.png'):
    """DoWhy refutation test results."""
    fig, ax = plt.subplots(figsize=(8, 5))

    original = dowhy_results['estimate']
    tests = ['Original\nEstimate', 'Placebo\nTreatment', 'Random\nCommon Cause', 'Data\nSubset']
    values = [
        original,
        dowhy_results['refute_placebo'].new_effect,
        dowhy_results['refute_rcc'].new_effect,
        dowhy_results['refute_subset'].new_effect,
    ]
    colors_r = ['#2196F3', '#E91E63', '#FF9800', '#4CAF50']

    bars = ax.bar(range(len(tests)), values, color=colors_r, alpha=0.8, edgecolor='black')
    ax.set_xticks(range(len(tests)))
    ax.set_xticklabels(tests, fontsize=11)
    ax.set_ylabel('Estimated Effect', fontsize=12)
    ax.set_title('DoWhy Refutation Tests', fontsize=14, fontweight='bold')
    ax.axhline(y=original, color='blue', linestyle='--', alpha=0.5, label=f'Original = {original:.2f}')
    ax.legend(fontsize=11)
    ax.grid(axis='y', alpha=0.3)

    for i, v in enumerate(values):
        ax.text(i, v + 0.1 if v > 0 else v - 0.3, f'{v:.2f}', ha='center', fontsize=10, fontweight='bold')

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_causal_dag(save_path='figures/causal_dag.png'):
    """Draw a causal DAG for the study design."""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.axis('off')

    nodes = {
        'Confounders\n(Age, BMI, BP,\nComorbidity, Smoking)': (5, 6),
        'Treatment\n(Drug)': (2, 3.5),
        'Outcome\n(Blood Pressure)': (8, 3.5),
        'Instrument\n(Physician Pref.)': (0.5, 5.5),
    }
    for label, (x, y) in nodes.items():
        bbox = FancyBboxPatch((x-1.2, y-0.6), 2.4, 1.2, boxstyle="round,pad=0.1",
                              facecolor='#E3F2FD', edgecolor='#1565C0', linewidth=2)
        ax.add_patch(bbox)
        ax.text(x, y, label, ha='center', va='center', fontsize=9, fontweight='bold')

    arrows = [
        ((5, 5.4), (2.8, 4.1)),    # Confounders -> Treatment
        ((5, 5.4), (7.2, 4.1)),    # Confounders -> Outcome
        ((3.2, 3.5), (6.8, 3.5)),  # Treatment -> Outcome
        ((1.3, 5.0), (1.8, 4.1)),  # Instrument -> Treatment
    ]
    for start, end in arrows:
        ax.annotate('', xy=end, xytext=start,
                    arrowprops=dict(arrowstyle='->', color='#1565C0', lw=2))

    ax.set_title('Causal DAG: Pharmacoepidemiology Study Design', fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_cate_heatmap(df, cate_cf, save_path='figures/cate_heatmap.png'):
    """CATE heatmap by age and BMI groups."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    df_temp = df.copy()
    df_temp['cate_cf'] = cate_cf.flatten()
    df_temp['true_cate'] = df['true_effect']
    df_temp['age_bin'] = pd.cut(df_temp['age'], bins=5)
    df_temp['bmi_bin'] = pd.cut(df_temp['bmi'], bins=5)

    pivot_cf = df_temp.pivot_table(values='cate_cf', index='bmi_bin', columns='age_bin', aggfunc='mean')
    pivot_true = df_temp.pivot_table(values='true_cate', index='bmi_bin', columns='age_bin', aggfunc='mean')

    im1 = axes[0].imshow(pivot_cf.values, cmap='RdBu_r', aspect='auto', interpolation='nearest')
    axes[0].set_title('Estimated CATE (Causal Forest)', fontsize=12, fontweight='bold')
    axes[0].set_xlabel('Age Group', fontsize=11)
    axes[0].set_ylabel('BMI Group', fontsize=11)
    axes[0].set_xticks(range(len(pivot_cf.columns)))
    axes[0].set_xticklabels([str(c) for c in pivot_cf.columns], rotation=45, ha='right', fontsize=8)
    axes[0].set_yticks(range(len(pivot_cf.index)))
    axes[0].set_yticklabels([str(i) for i in pivot_cf.index], fontsize=8)
    plt.colorbar(im1, ax=axes[0], shrink=0.8)

    im2 = axes[1].imshow(pivot_true.values, cmap='RdBu_r', aspect='auto', interpolation='nearest')
    axes[1].set_title('True CATE', fontsize=12, fontweight='bold')
    axes[1].set_xlabel('Age Group', fontsize=11)
    axes[1].set_ylabel('BMI Group', fontsize=11)
    axes[1].set_xticks(range(len(pivot_true.columns)))
    axes[1].set_xticklabels([str(c) for c in pivot_true.columns], rotation=45, ha='right', fontsize=8)
    axes[1].set_yticks(range(len(pivot_true.index)))
    axes[1].set_yticklabels([str(i) for i in pivot_true.index], fontsize=8)
    plt.colorbar(im2, ax=axes[1], shrink=0.8)

    plt.suptitle('CATE Heatmap: Age × BMI Subgroups', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


# ============================================================
# MAIN EXECUTION
# ============================================================

if __name__ == '__main__':
    print("=" * 70)
    print("Causal Inference Methods: Systematic Comparison Framework")
    print("=" * 70)

    TRUE_ATE = -2.5

    # Generate data
    print("\n[1] Generating pharmacoepidemiology simulation data...")
    df = generate_pharma_data(n=5000, true_ate=TRUE_ATE)
    df_did = generate_did_data(n_units=500, n_periods=10, true_effect=-3.0)
    print(f"    Cross-sectional: {len(df)} observations")
    print(f"    Panel data: {len(df_did)} observations")
    print(f"    Treatment prevalence: {df['treatment'].mean():.3f}")
    print(f"    True ATE: {TRUE_ATE}")
    print(f"    True CATE range: [{df['true_effect'].min():.2f}, {df['true_effect'].max():.2f}]")

    # Plot causal DAG
    print("\n[2] Drawing causal DAG...")
    plot_causal_dag()

    # Run all methods
    print("\n[3] Running Propensity Score Matching...")
    psm_result = propensity_score_matching(df)
    print(f"    PSM ATE: {psm_result['ate']:.4f} (SE: {psm_result['se']:.4f})")

    print("\n[4] Running Inverse Probability Weighting...")
    ipw_result = ipw_estimator(df)
    print(f"    IPW ATE: {ipw_result['ate']:.4f} (SE: {ipw_result['se']:.4f})")

    print("\n[5] Running IV-2SLS...")
    iv_result = iv_estimation(df)
    print(f"    IV ATE: {iv_result['ate']:.4f} (SE: {iv_result['se']:.4f})")
    print(f"    First-stage F-stat: {iv_result['f_stat']:.2f}")

    print("\n[6] Weak instrument analysis...")
    iv_weak_df = iv_weak_instrument_analysis(df)

    print("\n[7] Running DID analysis...")
    did_result = did_analysis(df_did)
    print(f"    DID ATE: {did_result['ate']:.4f} (SE: {did_result['se']:.4f})")

    print("\n[8] Running Double/Debiased ML...")
    dml_result = dml_estimation(df)
    print(f"    DML ATE: {dml_result['ate']:.4f} (SE: {dml_result['se']:.4f})")
    print(f"    DML 95% CI: [{dml_result['ci'][0]:.4f}, {dml_result['ci'][1]:.4f}]")

    print("\n[9] Running Causal Forest...")
    cf_result = causal_forest_estimation(df)
    print(f"    CF ATE: {cf_result['ate']:.4f} (SE: {cf_result['se']:.4f})")
    print(f"    CF 95% CI: [{cf_result['ci'][0]:.4f}, {cf_result['ci'][1]:.4f}]")

    print("\n[10] Running DoWhy workflow with refutation tests...")
    dowhy_result = dowhy_workflow(df)
    print(f"    DoWhy Backdoor estimate: {dowhy_result['estimate']:.4f}")

    # Collect results
    all_results = [psm_result, ipw_result, iv_result, dml_result, cf_result]

    # Generate all plots
    print("\n[11] Generating visualizations...")
    plot_method_comparison(all_results, TRUE_ATE)
    plot_propensity_distribution(df, psm_result['ps'])
    plot_weak_instrument(iv_weak_df)
    plot_did_trends(did_result['pre_trends'], did_result['ate'])
    plot_heterogeneous_effects(df, dml_result['cate'], cf_result['cate'])
    plot_bias_rmse(all_results, TRUE_ATE)
    plot_dowhy_refutation(dowhy_result)
    plot_cate_heatmap(df, cf_result['cate'])

    # Summary table
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"{'Method':<20} {'ATE':>10} {'SE':>10} {'Bias':>10} {'RMSE':>10}")
    print("-" * 60)
    for r in all_results:
        bias = r['ate'] - TRUE_ATE
        rmse = np.sqrt(bias**2 + r['se']**2)
        print(f"{r['method']:<20} {r['ate']:>10.4f} {r['se']:>10.4f} {bias:>10.4f} {rmse:>10.4f}")
    print(f"{'DID':<20} {did_result['ate']:>10.4f} {did_result['se']:>10.4f} {did_result['ate']-(-3.0):>10.4f} {np.sqrt((did_result['ate']-(-3.0))**2 + did_result['se']**2):>10.4f}")
    print(f"\nTrue ATE: {TRUE_ATE}")
    print(f"DoWhy Backdoor: {dowhy_result['estimate']:.4f}")

    # CATE correlation
    cate_corr = np.corrcoef(df['true_effect'], cf_result['cate'].flatten())[0, 1]
    cate_rmse = np.sqrt(np.mean((df['true_effect'] - cf_result['cate'].flatten())**2))
    print(f"\nCATE Recovery (Causal Forest):")
    print(f"  Correlation with true CATE: {cate_corr:.4f}")
    print(f"  RMSE of CATE: {cate_rmse:.4f}")

    dml_cate_corr = np.corrcoef(df['true_effect'], dml_result['cate'].flatten())[0, 1]
    dml_cate_rmse = np.sqrt(np.mean((df['true_effect'] - dml_result['cate'].flatten())**2))
    print(f"\nCATE Recovery (DML):")
    print(f"  Correlation with true CATE: {dml_cate_corr:.4f}")
    print(f"  RMSE of CATE: {dml_cate_rmse:.4f}")

    # Save results to CSV
    results_df = pd.DataFrame([{
        'Method': r['method'], 'ATE': r['ate'], 'SE': r['se'],
        'Bias': r['ate'] - TRUE_ATE,
        'RMSE': np.sqrt((r['ate'] - TRUE_ATE)**2 + r['se']**2)
    } for r in all_results])
    results_df.to_csv('results_summary.csv', index=False)
    print("\nResults saved to results_summary.csv")
    print("All figures saved to figures/ directory")
    print("=" * 70)
