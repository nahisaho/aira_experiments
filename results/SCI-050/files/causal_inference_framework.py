"""
Systematic Comparison Framework for Causal Effect Estimation from Observational Data
=====================================================================================
DoWhy/EconML-based causal inference workflow covering:
1. Propensity Score Matching (PSM) — limitations & alternatives
2. Instrumental Variables (IV) — weak instrument diagnostics
3. Difference-in-Differences (DID) — parallel trends verification
4. Double/Debiased Machine Learning (DML)
5. Causal Forest — heterogeneous treatment effects
6. Pharmacoepidemiology case study (RWD)

Random seeds fixed for reproducibility.
"""

import json
import warnings
import datetime
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
from scipy import stats

warnings.filterwarnings("ignore")
np.random.seed(42)

# ── Output directories ──
BASE = os.path.dirname(os.path.abspath(__file__))
for d in ["figures", "results", "data", "logs"]:
    os.makedirs(os.path.join(BASE, d), exist_ok=True)

LOG_ENTRIES = []

def log_event(phase, event_type, skill, handoff_in=None, handoff_out=None, files=None, status="ok"):
    LOG_ENTRIES.append({
        "timestamp": datetime.datetime.now().isoformat(),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill,
        "handoff_in": handoff_in or {},
        "handoff_out": handoff_out or {},
        "files_written": files or [],
        "status": status,
    })

log_event("setup", "run_started", "co-scientist-causal-inference")

# ============================================================
# DATA GENERATION — Pharmacoepidemiology RWD Simulation
# ============================================================
# Scenario: Estimating the effect of a new antihypertensive drug
# on cardiovascular event risk reduction using real-world data.
# True ATE = -5.0 (mmHg systolic BP reduction)

def generate_rwd_data(n=3000, true_ate=-5.0, seed=42):
    """Simulate pharmacoepidemiology RWD with confounding."""
    rng = np.random.RandomState(seed)

    age = rng.normal(60, 10, n)
    sex = rng.binomial(1, 0.45, n)  # 1=female
    bmi = rng.normal(28, 5, n)
    comorbidity_score = rng.poisson(2, n)
    baseline_bp = 140 + 0.3 * age + 2 * comorbidity_score + rng.normal(0, 8, n)

    # Instrument: physician preference (quasi-random)
    physician_preference = rng.binomial(1, 0.5, n)

    # Treatment assignment with confounding
    propensity_logit = (
        -2.0
        + 0.03 * age
        - 0.5 * sex
        + 0.05 * bmi
        + 0.3 * comorbidity_score
        - 0.02 * baseline_bp
        + 1.5 * physician_preference
    )
    propensity = 1 / (1 + np.exp(-propensity_logit))
    treatment = rng.binomial(1, propensity, n)

    # Heterogeneous treatment effect: larger for older, higher comorbidity
    ite = true_ate - 0.1 * (age - 60) - 0.5 * comorbidity_score
    noise = rng.normal(0, 5, n)

    # Outcome: post-treatment systolic BP
    outcome = (
        baseline_bp
        + treatment * ite
        - 0.1 * age
        + 1.5 * comorbidity_score
        + noise
    )

    # Time variable for DID
    time_period = rng.choice([0, 1], n, p=[0.5, 0.5])
    # Pre-post outcome with parallel trend + treatment effect in post
    outcome_did = (
        baseline_bp
        + 3 * time_period  # time trend
        + treatment * time_period * true_ate
        - 0.1 * age
        + noise
    )

    df = pd.DataFrame({
        "age": age, "sex": sex, "bmi": bmi,
        "comorbidity_score": comorbidity_score,
        "baseline_bp": baseline_bp,
        "physician_preference": physician_preference,
        "treatment": treatment,
        "outcome": outcome,
        "ite_true": ite,
        "propensity_true": propensity,
        "time_period": time_period,
        "outcome_did": outcome_did,
    })
    return df

print("=" * 70)
print("CAUSAL INFERENCE FRAMEWORK — Systematic Comparison")
print("=" * 70)

df = generate_rwd_data()
df.to_csv(os.path.join(BASE, "data", "simulated_rwd.csv"), index=False)
log_event("data_generation", "file_written", "data-simulation",
          files=["data/simulated_rwd.csv"])

print(f"\nDataset: n={len(df)}, treatment rate={df['treatment'].mean():.2%}")
print(f"True ATE = -5.0 mmHg")

RESULTS = {"true_ate": -5.0, "methods": {}}

# ============================================================
# 1. PROPENSITY SCORE MATCHING (PSM) — Limitations & Alternatives
# ============================================================
print("\n" + "─" * 70)
print("1. Propensity Score Matching (PSM)")
print("─" * 70)

covariates = ["age", "sex", "bmi", "comorbidity_score", "baseline_bp"]
X = df[covariates]
T = df["treatment"]
Y = df["outcome"]

# Estimate propensity score
ps_model = LogisticRegression(max_iter=1000, random_state=42)
ps_model.fit(X, T)
ps = ps_model.predict_proba(X)[:, 1]
df["ps_estimated"] = ps

# AUC for propensity model
auc = roc_auc_score(T, ps)
print(f"  Propensity model AUC: {auc:.3f}")

# 1-to-1 nearest neighbor matching
treated_idx = df[df["treatment"] == 1].index
control_idx = df[df["treatment"] == 0].index

nn = NearestNeighbors(n_neighbors=1, metric="euclidean")
nn.fit(df.loc[control_idx, ["ps_estimated"]])
distances, indices = nn.kneighbors(df.loc[treated_idx, ["ps_estimated"]])
matched_control_idx = control_idx[indices.flatten()]

ate_psm = (
    df.loc[treated_idx, "outcome"].values - df.loc[matched_control_idx, "outcome"].values
).mean()
ate_psm_se = (
    df.loc[treated_idx, "outcome"].values - df.loc[matched_control_idx, "outcome"].values
).std() / np.sqrt(len(treated_idx))

print(f"  PSM ATE estimate: {ate_psm:.3f} (SE: {ate_psm_se:.3f})")
print(f"  PSM 95% CI: [{ate_psm - 1.96*ate_psm_se:.3f}, {ate_psm + 1.96*ate_psm_se:.3f}]")

# Covariate balance check (SMD)
smd_before = {}
smd_after = {}
for cov in covariates:
    t_mean = df.loc[treated_idx, cov].mean()
    c_mean_before = df.loc[control_idx, cov].mean()
    c_mean_after = df.loc[matched_control_idx, cov].mean()
    pooled_std = np.sqrt(
        (df.loc[treated_idx, cov].std()**2 + df.loc[control_idx, cov].std()**2) / 2
    )
    smd_before[cov] = abs(t_mean - c_mean_before) / pooled_std
    smd_after[cov] = abs(t_mean - c_mean_after) / pooled_std

print(f"\n  Covariate Balance (SMD before/after matching):")
for cov in covariates:
    print(f"    {cov:20s}: {smd_before[cov]:.3f} → {smd_after[cov]:.3f}")

# IPW as alternative
weights = T / ps + (1 - T) / (1 - ps)
ate_ipw = (weights * Y * (2*T - 1)).sum() / len(df)
# Trimmed IPW (trim extreme weights)
trim_mask = (ps > 0.05) & (ps < 0.95)
df_trimmed = df[trim_mask]
ps_tr = df_trimmed["ps_estimated"]
T_tr = df_trimmed["treatment"]
Y_tr = df_trimmed["outcome"]
w_tr = T_tr / ps_tr + (1 - T_tr) / (1 - ps_tr)
ate_ipw_trimmed = (w_tr * Y_tr * (2*T_tr - 1)).sum() / len(df_trimmed)

# AIPW (Augmented IPW / Doubly Robust)
from sklearn.ensemble import GradientBoostingRegressor as GBR
mu1_model = GBR(n_estimators=100, random_state=42).fit(X[T==1], Y[T==1])
mu0_model = GBR(n_estimators=100, random_state=42).fit(X[T==0], Y[T==0])
mu1_hat = mu1_model.predict(X)
mu0_hat = mu0_model.predict(X)

aipw_scores = (
    mu1_hat - mu0_hat
    + T * (Y - mu1_hat) / ps
    - (1 - T) * (Y - mu0_hat) / (1 - ps)
)
ate_aipw = aipw_scores.mean()
ate_aipw_se = aipw_scores.std() / np.sqrt(len(df))

print(f"\n  Alternative estimators:")
print(f"    IPW ATE:          {ate_ipw:.3f}")
print(f"    Trimmed IPW ATE:  {ate_ipw_trimmed:.3f}")
print(f"    AIPW (DR) ATE:    {ate_aipw:.3f} (SE: {ate_aipw_se:.3f})")

RESULTS["methods"]["PSM"] = {
    "ate": round(ate_psm, 3), "se": round(ate_psm_se, 3),
    "auc": round(auc, 3), "n_matched": int(len(treated_idx)),
}
RESULTS["methods"]["IPW"] = {"ate": round(ate_ipw, 3)}
RESULTS["methods"]["IPW_trimmed"] = {"ate": round(ate_ipw_trimmed, 3)}
RESULTS["methods"]["AIPW"] = {
    "ate": round(ate_aipw, 3), "se": round(ate_aipw_se, 3),
}

# ── Figure: PS distribution & balance ──
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

axes[0].hist(ps[T==1], bins=40, alpha=0.6, label="Treated", density=True, color="#2196F3")
axes[0].hist(ps[T==0], bins=40, alpha=0.6, label="Control", density=True, color="#FF9800")
axes[0].set_xlabel("Propensity Score")
axes[0].set_ylabel("Density")
axes[0].set_title("Propensity Score Distribution")
axes[0].legend()

cov_labels = [c[:12] for c in covariates]
y_pos = np.arange(len(covariates))
axes[1].barh(y_pos - 0.15, [smd_before[c] for c in covariates], 0.3,
             label="Before Matching", color="#e57373")
axes[1].barh(y_pos + 0.15, [smd_after[c] for c in covariates], 0.3,
             label="After Matching", color="#81c784")
axes[1].axvline(x=0.1, color="red", linestyle="--", label="SMD=0.1 threshold")
axes[1].set_yticks(y_pos)
axes[1].set_yticklabels(cov_labels)
axes[1].set_xlabel("Standardized Mean Difference")
axes[1].set_title("Covariate Balance (SMD)")
axes[1].legend(fontsize=8)

methods_psm = ["PSM", "IPW", "IPW\n(trimmed)", "AIPW"]
ates_psm = [ate_psm, ate_ipw, ate_ipw_trimmed, ate_aipw]
colors_psm = ["#2196F3", "#FF9800", "#4CAF50", "#9C27B0"]
axes[2].bar(methods_psm, ates_psm, color=colors_psm, alpha=0.8)
axes[2].axhline(y=-5.0, color="red", linestyle="--", label="True ATE (-5.0)")
axes[2].set_ylabel("Estimated ATE (mmHg)")
axes[2].set_title("PSM vs Alternative Estimators")
axes[2].legend()

plt.tight_layout()
plt.savefig(os.path.join(BASE, "figures", "01_psm_analysis.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(BASE, "figures", "01_psm_analysis.svg"), bbox_inches="tight")
plt.close()
print("  → Saved figures/01_psm_analysis.png")

# ============================================================
# 2. INSTRUMENTAL VARIABLES (IV) — Weak Instrument Diagnostics
# ============================================================
print("\n" + "─" * 70)
print("2. Instrumental Variables (IV)")
print("─" * 70)

Z = df["physician_preference"]

# First stage: T ~ Z + covariates
X_iv = sm.add_constant(pd.concat([Z, df[covariates]], axis=1))
first_stage = sm.OLS(T, X_iv).fit()
f_stat_instrument = first_stage.fvalue
print(f"  First-stage F-statistic: {f_stat_instrument:.2f}")
print(f"  Instrument coefficient: {first_stage.params.iloc[1]:.4f} (p={first_stage.pvalues.iloc[1]:.4e})")

# Stock-Yogo weak instrument test
is_weak = f_stat_instrument < 10
print(f"  Weak instrument (F < 10)? {'YES ⚠️' if is_weak else 'NO ✓'}")

# 2SLS
T_hat = first_stage.predict(X_iv)
X_second = sm.add_constant(pd.concat([pd.Series(T_hat, name="treatment_hat"),
                                       df[covariates]], axis=1))
second_stage = sm.OLS(Y, X_second).fit()
ate_iv = second_stage.params.iloc[1]
ate_iv_se = second_stage.bse.iloc[1]

print(f"  2SLS ATE estimate: {ate_iv:.3f} (SE: {ate_iv_se:.3f})")
print(f"  2SLS 95% CI: [{ate_iv - 1.96*ate_iv_se:.3f}, {ate_iv + 1.96*ate_iv_se:.3f}]")

# Anderson-Rubin confidence set (robust to weak instruments)
# Simplified: test H0: beta = beta0 for a grid
beta_grid = np.linspace(-15, 5, 200)
ar_pvals = []
for b0 in beta_grid:
    resid = Y - b0 * T
    X_ar = sm.add_constant(pd.concat([Z, df[covariates]], axis=1))
    ar_reg = sm.OLS(resid, X_ar).fit()
    ar_f = ar_reg.f_test("physician_preference = 0")
    ar_pvals.append(float(ar_f.pvalue))

ar_pvals = np.array(ar_pvals)
ar_ci_mask = ar_pvals > 0.05
if ar_ci_mask.any():
    ar_ci_low = beta_grid[ar_ci_mask].min()
    ar_ci_high = beta_grid[ar_ci_mask].max()
else:
    ar_ci_low, ar_ci_high = np.nan, np.nan
print(f"  Anderson-Rubin 95% CI: [{ar_ci_low:.3f}, {ar_ci_high:.3f}]")

# Sensitivity to instrument strength — simulate varying strengths
iv_strengths = [0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0]
iv_ates = []
iv_fstats = []
for gamma in iv_strengths:
    rng_s = np.random.RandomState(42)
    propensity_s = 1 / (1 + np.exp(-(
        -2.0 + 0.03*df["age"] - 0.5*df["sex"] + 0.05*df["bmi"]
        + 0.3*df["comorbidity_score"] - 0.02*df["baseline_bp"]
        + gamma * df["physician_preference"]
    )))
    T_s = rng_s.binomial(1, propensity_s)
    X_iv_s = sm.add_constant(pd.concat([Z, df[covariates]], axis=1))
    fs = sm.OLS(T_s, X_iv_s).fit()
    iv_fstats.append(fs.fvalue)
    T_hat_s = fs.predict(X_iv_s)
    X_ss = sm.add_constant(pd.concat([pd.Series(T_hat_s, name="t_hat"),
                                       df[covariates]], axis=1))
    ite_s = df["ite_true"].values
    Y_s = df["baseline_bp"] + T_s * ite_s - 0.1*df["age"] + np.random.RandomState(42).normal(0,5,len(df))
    ss = sm.OLS(Y_s, X_ss).fit()
    iv_ates.append(ss.params.iloc[1])

RESULTS["methods"]["IV_2SLS"] = {
    "ate": round(ate_iv, 3), "se": round(ate_iv_se, 3),
    "first_stage_F": round(f_stat_instrument, 2),
    "weak_instrument": is_weak,
    "AR_CI": [round(ar_ci_low, 3), round(ar_ci_high, 3)],
}

# ── Figure: IV diagnostics ──
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

axes[0].scatter(Z + np.random.normal(0, 0.05, len(Z)), T + np.random.normal(0, 0.05, len(T)),
                alpha=0.1, s=5, color="#2196F3")
axes[0].set_xlabel("Physician Preference (Instrument)")
axes[0].set_ylabel("Treatment Assignment")
axes[0].set_title(f"First Stage (F={f_stat_instrument:.1f})")

axes[1].plot(beta_grid, ar_pvals, color="#9C27B0", linewidth=1.5)
axes[1].axhline(y=0.05, color="red", linestyle="--", label="α=0.05")
axes[1].axvline(x=-5.0, color="green", linestyle="--", label="True ATE")
axes[1].fill_between(beta_grid, 0, 1, where=ar_ci_mask, alpha=0.15, color="#9C27B0",
                      label="AR 95% CS")
axes[1].set_xlabel("Causal Effect (β)")
axes[1].set_ylabel("p-value")
axes[1].set_title("Anderson-Rubin Confidence Set")
axes[1].legend(fontsize=8)

axes[2].plot(iv_fstats, iv_ates, "o-", color="#FF5722", markersize=6)
axes[2].axhline(y=-5.0, color="green", linestyle="--", label="True ATE")
axes[2].axvline(x=10, color="red", linestyle="--", alpha=0.6, label="F=10 threshold")
for i, g in enumerate(iv_strengths):
    axes[2].annotate(f"γ={g}", (iv_fstats[i], iv_ates[i]), fontsize=7,
                      textcoords="offset points", xytext=(5, 5))
axes[2].set_xlabel("First-Stage F-statistic")
axes[2].set_ylabel("2SLS ATE Estimate")
axes[2].set_title("IV Sensitivity to Instrument Strength")
axes[2].legend(fontsize=8)

plt.tight_layout()
plt.savefig(os.path.join(BASE, "figures", "02_iv_analysis.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(BASE, "figures", "02_iv_analysis.svg"), bbox_inches="tight")
plt.close()
print("  → Saved figures/02_iv_analysis.png")

# ============================================================
# 3. DIFFERENCE-IN-DIFFERENCES (DID) — Parallel Trends
# ============================================================
print("\n" + "─" * 70)
print("3. Difference-in-Differences (DID)")
print("─" * 70)

# Generate multi-period data for parallel trends test
n_periods = 8
rng_did = np.random.RandomState(42)
treatment_period = 5  # treatment starts at period 5

did_records = []
n_units = 500
unit_fe = rng_did.normal(0, 3, n_units)
treatment_units = rng_did.binomial(1, 0.4, n_units)

for t in range(n_periods):
    for i in range(n_units):
        post = int(t >= treatment_period)
        treated = treatment_units[i]
        y = (
            100 + unit_fe[i]
            + 2.0 * t  # common trend
            + (-5.0) * treated * post  # treatment effect
            + rng_did.normal(0, 3)
        )
        did_records.append({
            "unit": i, "period": t, "treatment": treated,
            "post": post, "outcome": y,
        })

df_did = pd.DataFrame(did_records)
df_did.to_csv(os.path.join(BASE, "data", "did_panel.csv"), index=False)

# DID regression
df_did["treat_post"] = df_did["treatment"] * df_did["post"]
X_did = sm.add_constant(df_did[["treatment", "post", "treat_post"]])
did_model = sm.OLS(df_did["outcome"], X_did).fit(cov_type="cluster",
                                                    cov_kwds={"groups": df_did["unit"]})
ate_did = did_model.params["treat_post"]
ate_did_se = did_model.bse["treat_post"]
print(f"  DID ATE estimate: {ate_did:.3f} (SE: {ate_did_se:.3f})")
print(f"  DID 95% CI: [{ate_did - 1.96*ate_did_se:.3f}, {ate_did + 1.96*ate_did_se:.3f}]")

# Parallel trends test: event study
event_study_coefs = []
event_study_se = []
event_study_periods = []
for t in range(n_periods):
    if t == treatment_period - 1:  # reference period
        event_study_coefs.append(0)
        event_study_se.append(0)
        event_study_periods.append(t - treatment_period)
        continue
    df_did[f"d_{t}"] = ((df_did["period"] == t) & (df_did["treatment"] == 1)).astype(int)
    event_study_periods.append(t - treatment_period)

# Event study regression
event_cols = [f"d_{t}" for t in range(n_periods) if t != treatment_period - 1]
X_event = sm.add_constant(df_did[["treatment", "post"] + event_cols])
event_model = sm.OLS(df_did["outcome"], X_event).fit(cov_type="cluster",
                                                       cov_kwds={"groups": df_did["unit"]})

# Rebuild with regression coefficients
event_study_coefs = []
event_study_se = []
event_study_periods = []
for t in range(n_periods):
    rel = t - treatment_period
    event_study_periods.append(rel)
    if t == treatment_period - 1:
        event_study_coefs.append(0.0)
        event_study_se.append(0.0)
    else:
        col = f"d_{t}"
        event_study_coefs.append(event_model.params[col])
        event_study_se.append(event_model.bse[col])

# Pre-trend test: joint F-test on pre-treatment dummies
pre_cols = [f"d_{t}" for t in range(treatment_period - 1)]
if pre_cols:
    f_test_str = " = ".join([f"{c} = 0" for c in pre_cols])
    try:
        pre_f_test = event_model.f_test(" = ".join(pre_cols) + " = 0")
        pre_f_pval = float(pre_f_test.pvalue)
    except:
        # manual test
        pre_coefs = [event_model.params[c] for c in pre_cols]
        pre_f_pval = 1 - stats.chi2.cdf(sum(c**2/se**2 for c, se in
                                              zip(pre_coefs, [event_model.bse[c] for c in pre_cols])),
                                          df=len(pre_cols))
    print(f"  Pre-trend joint test p-value: {pre_f_pval:.4f}")
    print(f"  Parallel trends {'supported ✓' if pre_f_pval > 0.05 else 'violated ⚠️'}")

RESULTS["methods"]["DID"] = {
    "ate": round(ate_did, 3), "se": round(ate_did_se, 3),
    "pre_trend_pvalue": round(pre_f_pval, 4),
    "parallel_trends_supported": pre_f_pval > 0.05,
}

# ── Figure: DID & Event Study ──
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Mean outcome by group and period
for grp, label, color in [(1, "Treated", "#2196F3"), (0, "Control", "#FF9800")]:
    means = df_did[df_did["treatment"] == grp].groupby("period")["outcome"].mean()
    axes[0].plot(means.index, means.values, "o-", label=label, color=color, linewidth=2)
axes[0].axvline(x=treatment_period - 0.5, color="red", linestyle="--", label="Treatment Start")
axes[0].set_xlabel("Time Period")
axes[0].set_ylabel("Mean Outcome")
axes[0].set_title("DID: Group Trends Over Time")
axes[0].legend()

# Event study plot
es_coefs = np.array(event_study_coefs)
es_se = np.array(event_study_se)
es_periods = np.array(event_study_periods)
axes[1].errorbar(es_periods, es_coefs, yerr=1.96*es_se, fmt="o-", color="#9C27B0",
                  capsize=4, linewidth=2, markersize=6)
axes[1].axhline(y=0, color="gray", linestyle="-", alpha=0.5)
axes[1].axvline(x=-0.5, color="red", linestyle="--", alpha=0.6, label="Treatment Start")
axes[1].fill_between(es_periods[es_periods < 0], es_coefs[es_periods < 0] - 1.96*es_se[es_periods < 0],
                      es_coefs[es_periods < 0] + 1.96*es_se[es_periods < 0], alpha=0.15, color="green",
                      label="Pre-treatment")
axes[1].set_xlabel("Relative Time Period")
axes[1].set_ylabel("Coefficient (Treatment Effect)")
axes[1].set_title("Event Study Plot")
axes[1].legend(fontsize=8)

# Placebo test: random treatment assignment
n_placebo = 500
placebo_ates = []
for _ in range(n_placebo):
    fake_treat = np.random.permutation(treatment_units)
    fake_effects = []
    for t in range(n_periods):
        for i in range(n_units):
            post = int(t >= treatment_period)
            did_records_fake = fake_treat[i] * post
    # Simplified: reshuffle treatment and re-estimate
    df_placebo = df_did.copy()
    fake_map = dict(zip(range(n_units), np.random.permutation(treatment_units)))
    df_placebo["treatment"] = df_placebo["unit"].map(fake_map)
    df_placebo["treat_post"] = df_placebo["treatment"] * df_placebo["post"]
    X_pl = sm.add_constant(df_placebo[["treatment", "post", "treat_post"]])
    pl_model = sm.OLS(df_placebo["outcome"], X_pl).fit()
    placebo_ates.append(pl_model.params["treat_post"])

axes[2].hist(placebo_ates, bins=40, density=True, alpha=0.7, color="#78909C",
              label="Placebo Distribution")
axes[2].axvline(x=ate_did, color="red", linewidth=2, label=f"Actual DID={ate_did:.2f}")
axes[2].axvline(x=-5.0, color="green", linestyle="--", label="True ATE")
p_val_placebo = np.mean(np.array(placebo_ates) <= ate_did)
axes[2].set_xlabel("Placebo ATE Estimates")
axes[2].set_ylabel("Density")
axes[2].set_title(f"Placebo Test (p={p_val_placebo:.3f})")
axes[2].legend(fontsize=8)

plt.tight_layout()
plt.savefig(os.path.join(BASE, "figures", "03_did_analysis.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(BASE, "figures", "03_did_analysis.svg"), bbox_inches="tight")
plt.close()
print("  → Saved figures/03_did_analysis.png")

# ============================================================
# 4. DOUBLE/DEBIASED MACHINE LEARNING (DML)
# ============================================================
print("\n" + "─" * 70)
print("4. Double/Debiased Machine Learning (DML)")
print("─" * 70)

from econml.dml import LinearDML, CausalForestDML
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier

X_dml = df[covariates].values
T_dml = df["treatment"].values
Y_dml = df["outcome"].values

# LinearDML
ldml = LinearDML(
    model_y=GradientBoostingRegressor(n_estimators=100, max_depth=4, random_state=42),
    model_t=GradientBoostingClassifier(n_estimators=100, max_depth=4, random_state=42),
    discrete_treatment=True,
    random_state=42,
    cv=5,
)
ldml.fit(Y_dml, T_dml, X=X_dml)

ate_ldml = ldml.ate(X=X_dml)
ate_ldml_ci = ldml.ate_interval(X=X_dml, alpha=0.05)

print(f"  LinearDML ATE: {ate_ldml:.3f}")
print(f"  LinearDML 95% CI: [{ate_ldml_ci[0]:.3f}, {ate_ldml_ci[1]:.3f}]")

# CausalForestDML for HTE
cfdml = CausalForestDML(
    model_y=GradientBoostingRegressor(n_estimators=100, max_depth=4, random_state=42),
    model_t=GradientBoostingClassifier(n_estimators=100, max_depth=4, random_state=42),
    discrete_treatment=True,
    n_estimators=200,
    random_state=42,
    cv=5,
)
cfdml.fit(Y_dml, T_dml, X=X_dml)

ate_cfdml = cfdml.ate(X=X_dml)
ate_cfdml_ci = cfdml.ate_interval(X=X_dml, alpha=0.05)

print(f"  CausalForestDML ATE: {ate_cfdml:.3f}")
print(f"  CausalForestDML 95% CI: [{ate_cfdml_ci[0]:.3f}, {ate_cfdml_ci[1]:.3f}]")

# CATE estimates
cate_dml = cfdml.effect(X=X_dml).flatten()
cate_ci = cfdml.effect_interval(X=X_dml, alpha=0.05)

RESULTS["methods"]["LinearDML"] = {
    "ate": round(float(ate_ldml), 3),
    "ci": [round(float(ate_ldml_ci[0]), 3), round(float(ate_ldml_ci[1]), 3)],
}
RESULTS["methods"]["CausalForestDML"] = {
    "ate": round(float(ate_cfdml), 3),
    "ci": [round(float(ate_cfdml_ci[0]), 3), round(float(ate_cfdml_ci[1]), 3)],
}

# ── Figure: DML results ──
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# CATE distribution
axes[0].hist(cate_dml, bins=50, density=True, alpha=0.7, color="#2196F3", label="Estimated CATE")
axes[0].hist(df["ite_true"], bins=50, density=True, alpha=0.5, color="#4CAF50", label="True ITE")
axes[0].axvline(x=-5.0, color="red", linestyle="--", label="True ATE")
axes[0].set_xlabel("Treatment Effect (mmHg)")
axes[0].set_ylabel("Density")
axes[0].set_title("DML: CATE Distribution")
axes[0].legend(fontsize=8)

# CATE vs true ITE
axes[1].scatter(df["ite_true"], cate_dml, alpha=0.15, s=8, color="#9C27B0")
min_v = min(df["ite_true"].min(), cate_dml.min())
max_v = max(df["ite_true"].max(), cate_dml.max())
axes[1].plot([min_v, max_v], [min_v, max_v], "r--", label="Perfect calibration")
corr = np.corrcoef(df["ite_true"], cate_dml)[0, 1]
axes[1].set_xlabel("True ITE")
axes[1].set_ylabel("Estimated CATE")
axes[1].set_title(f"CATE Calibration (r={corr:.3f})")
axes[1].legend()

# Feature importance for HTE
try:
    importances = cfdml.feature_importances_
except:
    importances = np.abs(np.corrcoef(X_dml.T, cate_dml)[-1, :-1])
sort_idx = np.argsort(importances)
axes[2].barh(np.array(covariates)[sort_idx], importances[sort_idx], color="#FF5722", alpha=0.8)
axes[2].set_xlabel("Feature Importance")
axes[2].set_title("HTE Drivers (CausalForestDML)")

plt.tight_layout()
plt.savefig(os.path.join(BASE, "figures", "04_dml_analysis.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(BASE, "figures", "04_dml_analysis.svg"), bbox_inches="tight")
plt.close()
print("  → Saved figures/04_dml_analysis.png")

# ============================================================
# 5. CAUSAL FOREST — Heterogeneous Treatment Effects
# ============================================================
print("\n" + "─" * 70)
print("5. Causal Forest — Heterogeneous Treatment Effects")
print("─" * 70)

from econml.grf import CausalForest

cf = CausalForest(
    n_estimators=500,
    max_depth=None,
    min_samples_leaf=20,
    random_state=42,
)
cf.fit(X_dml, T_dml.flatten(), Y_dml)

cate_cf = cf.predict(X_dml).flatten()
ate_cf = cate_cf.mean()

# Bootstrap CI for ATE
n_boot = 200
boot_ates = []
rng_boot = np.random.RandomState(42)
for _ in range(n_boot):
    idx_b = rng_boot.choice(len(X_dml), len(X_dml), replace=True)
    boot_ates.append(cate_cf[idx_b].mean())
ate_cf_ci = np.percentile(boot_ates, [2.5, 97.5])

print(f"  Causal Forest ATE: {ate_cf:.3f}")
print(f"  Bootstrap 95% CI: [{ate_cf_ci[0]:.3f}, {ate_cf_ci[1]:.3f}]")
print(f"  CATE range: [{cate_cf.min():.3f}, {cate_cf.max():.3f}]")

# Subgroup analysis
df["cate_cf"] = cate_cf
df["cate_dml"] = cate_dml

# Quintile analysis
df["cate_quintile"] = pd.qcut(cate_cf, 5, labels=["Q1\n(least)", "Q2", "Q3", "Q4", "Q5\n(most)"])
quintile_summary = df.groupby("cate_quintile").agg(
    mean_cate=("cate_cf", "mean"),
    mean_age=("age", "mean"),
    mean_comorbidity=("comorbidity_score", "mean"),
    mean_bmi=("bmi", "mean"),
    n=("cate_cf", "count"),
).reset_index()
print(f"\n  Subgroup Analysis by CATE Quintile:")
print(quintile_summary.to_string(index=False))

RESULTS["methods"]["CausalForest"] = {
    "ate": round(ate_cf, 3),
    "ci": [round(ate_cf_ci[0], 3), round(ate_cf_ci[1], 3)],
    "cate_range": [round(cate_cf.min(), 3), round(cate_cf.max(), 3)],
    "cate_ite_correlation": round(float(np.corrcoef(df["ite_true"], cate_cf)[0, 1]), 3),
}

# ── Figure: Causal Forest HTE ──
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# CATE by age
age_bins = pd.cut(df["age"], bins=10)
age_cate = df.groupby(age_bins).agg(
    mean_cate=("cate_cf", "mean"),
    mean_ite=("ite_true", "mean"),
    se_cate=("cate_cf", "sem"),
).reset_index()
x_pos = range(len(age_cate))
axes[0, 0].errorbar(x_pos, age_cate["mean_cate"], yerr=1.96*age_cate["se_cate"],
                     fmt="o-", color="#2196F3", capsize=3, label="Estimated CATE")
axes[0, 0].plot(x_pos, age_cate["mean_ite"], "s--", color="#4CAF50", label="True ITE")
axes[0, 0].set_xticks(x_pos)
axes[0, 0].set_xticklabels([str(b)[:8] for b in age_cate["age"]], rotation=45, fontsize=7)
axes[0, 0].set_xlabel("Age Bin")
axes[0, 0].set_ylabel("Treatment Effect (mmHg)")
axes[0, 0].set_title("HTE by Age")
axes[0, 0].legend(fontsize=8)

# CATE by comorbidity
comorb_cate = df.groupby("comorbidity_score").agg(
    mean_cate=("cate_cf", "mean"),
    mean_ite=("ite_true", "mean"),
    se_cate=("cate_cf", "sem"),
    n=("cate_cf", "count"),
).reset_index()
comorb_cate = comorb_cate[comorb_cate["n"] >= 30]
axes[0, 1].errorbar(comorb_cate["comorbidity_score"], comorb_cate["mean_cate"],
                     yerr=1.96*comorb_cate["se_cate"], fmt="o-", color="#FF5722",
                     capsize=3, label="Estimated CATE")
axes[0, 1].plot(comorb_cate["comorbidity_score"], comorb_cate["mean_ite"],
                 "s--", color="#4CAF50", label="True ITE")
axes[0, 1].set_xlabel("Comorbidity Score")
axes[0, 1].set_ylabel("Treatment Effect (mmHg)")
axes[0, 1].set_title("HTE by Comorbidity Score")
axes[0, 1].legend(fontsize=8)

# Quintile analysis
axes[1, 0].bar(quintile_summary["cate_quintile"], quintile_summary["mean_cate"],
                color=plt.cm.viridis(np.linspace(0.2, 0.8, 5)), alpha=0.8)
axes[1, 0].axhline(y=-5.0, color="red", linestyle="--", label="True ATE")
axes[1, 0].set_xlabel("CATE Quintile")
axes[1, 0].set_ylabel("Mean CATE (mmHg)")
axes[1, 0].set_title("Treatment Effect by Quintile")
axes[1, 0].legend()

# CATE heatmap: age × comorbidity
hte_pivot = df.pivot_table(values="cate_cf",
                            index=pd.cut(df["age"], bins=8),
                            columns=pd.cut(df["comorbidity_score"], bins=5),
                            aggfunc="mean")
im = axes[1, 1].imshow(hte_pivot.values, cmap="RdYlBu_r", aspect="auto")
axes[1, 1].set_xticks(range(hte_pivot.shape[1]))
axes[1, 1].set_xticklabels([str(c)[:6] for c in hte_pivot.columns], rotation=45, fontsize=7)
axes[1, 1].set_yticks(range(hte_pivot.shape[0]))
axes[1, 1].set_yticklabels([str(r)[:8] for r in hte_pivot.index], fontsize=7)
axes[1, 1].set_xlabel("Comorbidity Score")
axes[1, 1].set_ylabel("Age")
axes[1, 1].set_title("CATE Heatmap: Age × Comorbidity")
plt.colorbar(im, ax=axes[1, 1], label="CATE (mmHg)")

plt.tight_layout()
plt.savefig(os.path.join(BASE, "figures", "05_causal_forest.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(BASE, "figures", "05_causal_forest.svg"), bbox_inches="tight")
plt.close()
print("  → Saved figures/05_causal_forest.png")

# ============================================================
# 6. DoWhy WORKFLOW & COMPREHENSIVE COMPARISON
# ============================================================
print("\n" + "─" * 70)
print("6. DoWhy Causal Model & Method Comparison")
print("─" * 70)

import dowhy
from dowhy import CausalModel

# Build causal model with DAG
causal_graph = """
digraph {
    age -> treatment;
    age -> outcome;
    sex -> treatment;
    sex -> outcome;
    bmi -> treatment;
    bmi -> outcome;
    comorbidity_score -> treatment;
    comorbidity_score -> outcome;
    baseline_bp -> treatment;
    baseline_bp -> outcome;
    physician_preference -> treatment;
    treatment -> outcome;
}
"""

model = CausalModel(
    data=df,
    treatment="treatment",
    outcome="outcome",
    graph=causal_graph,
    instruments=["physician_preference"],
)

print("  DoWhy identified estimands:")
identified = model.identify_effect(proceed_when_unidentifiable=True)
print(f"    Backdoor variables: {identified.get_backdoor_variables()}")
print(f"    IV variables: {identified.get_instrumental_variables()}")

# Estimate with multiple methods
dowhy_results = {}

# Backdoor (linear regression)
est_lr = model.estimate_effect(
    identified,
    method_name="backdoor.linear_regression",
)
dowhy_results["DoWhy_LinearRegression"] = round(est_lr.value, 3)
print(f"  DoWhy Linear Regression: {est_lr.value:.3f}")

# Backdoor (propensity score matching)
est_psm = model.estimate_effect(
    identified,
    method_name="backdoor.propensity_score_matching",
)
dowhy_results["DoWhy_PSM"] = round(est_psm.value, 3)
print(f"  DoWhy PSM: {est_psm.value:.3f}")

# IV
est_iv = model.estimate_effect(
    identified,
    method_name="iv.instrumental_variable",
    method_params={"iv_instrument_name": "physician_preference"},
)
dowhy_results["DoWhy_IV"] = round(est_iv.value, 3)
print(f"  DoWhy IV: {est_iv.value:.3f}")

# Refutation tests
print("\n  Refutation Tests:")
refute_random = model.refute_estimate(
    identified, est_lr,
    method_name="random_common_cause",
    random_seed=42,
)
print(f"    Random common cause: effect={refute_random.new_effect:.3f}, p={refute_random.refutation_result['p_value']:.4f}")

refute_placebo_effect = "N/A (incompatible estimator)"
refute_placebo_p = "N/A"
try:
    refute_placebo = model.refute_estimate(
        identified, est_lr,
        method_name="placebo_treatment_refuter",
        placebo_type="permute",
        random_seed=42,
    )
    refute_placebo_effect = f"{refute_placebo.new_effect:.3f}"
    refute_placebo_p = f"{refute_placebo.refutation_result['p_value']:.4f}"
except Exception as e:
    pass
print(f"    Placebo treatment: effect={refute_placebo_effect}, p={refute_placebo_p}")

refute_subset = model.refute_estimate(
    identified, est_lr,
    method_name="data_subset_refuter",
    subset_fraction=0.8,
    random_seed=42,
)
print(f"    Data subset: effect={refute_subset.new_effect:.3f}, p={refute_subset.refutation_result['p_value']:.4f}")

RESULTS["methods"]["DoWhy"] = dowhy_results
RESULTS["refutation"] = {
    "random_common_cause": {
        "effect": round(refute_random.new_effect, 3),
        "p_value": round(refute_random.refutation_result['p_value'], 4),
    },
    "placebo_treatment": {
        "effect": refute_placebo_effect,
        "p_value": refute_placebo_p,
    },
    "data_subset": {
        "effect": round(refute_subset.new_effect, 3),
        "p_value": round(refute_subset.refutation_result['p_value'], 4),
    },
}

# ── Comprehensive Comparison Figure ──
all_methods = {
    "PSM": (RESULTS["methods"]["PSM"]["ate"], RESULTS["methods"]["PSM"]["se"]),
    "IPW": (RESULTS["methods"]["IPW"]["ate"], None),
    "IPW\n(trimmed)": (RESULTS["methods"]["IPW_trimmed"]["ate"], None),
    "AIPW": (RESULTS["methods"]["AIPW"]["ate"], RESULTS["methods"]["AIPW"]["se"]),
    "2SLS-IV": (RESULTS["methods"]["IV_2SLS"]["ate"], RESULTS["methods"]["IV_2SLS"]["se"]),
    "DID": (RESULTS["methods"]["DID"]["ate"], RESULTS["methods"]["DID"]["se"]),
    "LinearDML": (RESULTS["methods"]["LinearDML"]["ate"], None),
    "CF-DML": (RESULTS["methods"]["CausalForestDML"]["ate"], None),
    "Causal\nForest": (RESULTS["methods"]["CausalForest"]["ate"], None),
    "DoWhy\nLR": (dowhy_results["DoWhy_LinearRegression"], None),
    "DoWhy\nPSM": (dowhy_results["DoWhy_PSM"], None),
    "DoWhy\nIV": (dowhy_results["DoWhy_IV"], None),
}

fig, axes = plt.subplots(1, 2, figsize=(18, 7))

method_names = list(all_methods.keys())
method_ates = [all_methods[m][0] for m in method_names]
method_ses = [all_methods[m][1] if all_methods[m][1] else 0 for m in method_names]

colors_map = plt.cm.Set3(np.linspace(0, 1, len(method_names)))
bars = axes[0].barh(method_names, method_ates, xerr=[1.96*s for s in method_ses],
                     color=colors_map, capsize=3, alpha=0.85)
axes[0].axvline(x=-5.0, color="red", linewidth=2, linestyle="--", label="True ATE = -5.0")
axes[0].set_xlabel("Estimated ATE (mmHg)", fontsize=12)
axes[0].set_title("Comprehensive Method Comparison", fontsize=14)
axes[0].legend(fontsize=10)

# Bias comparison
biases = [abs(ate - (-5.0)) for ate in method_ates]
sort_idx = np.argsort(biases)
sorted_names = [method_names[i] for i in sort_idx]
sorted_biases = [biases[i] for i in sort_idx]
colors_bias = plt.cm.RdYlGn_r(np.linspace(0.1, 0.9, len(sorted_names)))
axes[1].barh(sorted_names, sorted_biases, color=colors_bias, alpha=0.85)
axes[1].set_xlabel("|Bias| = |Estimated ATE - True ATE|", fontsize=12)
axes[1].set_title("Estimation Bias (sorted)", fontsize=14)
axes[1].axvline(x=0.5, color="green", linestyle="--", alpha=0.6, label="|Bias| = 0.5")
axes[1].legend(fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(BASE, "figures", "06_comprehensive_comparison.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(BASE, "figures", "06_comprehensive_comparison.svg"), bbox_inches="tight")
plt.close()
print("  → Saved figures/06_comprehensive_comparison.png")

# ============================================================
# SAVE RESULTS
# ============================================================
with open(os.path.join(BASE, "results", "causal_estimates.json"), "w") as f:
    json.dump(RESULTS, f, indent=2, default=str)

quintile_summary.to_csv(os.path.join(BASE, "results", "hte_quintile_summary.csv"), index=False)

df[["age", "sex", "bmi", "comorbidity_score", "baseline_bp", "treatment",
    "outcome", "ite_true", "cate_cf", "cate_dml"]].to_csv(
    os.path.join(BASE, "results", "individual_cate_estimates.csv"), index=False
)

# Comparison summary table
comparison_df = pd.DataFrame([
    {"Method": name, "ATE": vals[0],
     "SE": vals[1] if vals[1] else "—",
     "Bias": round(abs(vals[0] - (-5.0)), 3)}
    for name, vals in all_methods.items()
])
comparison_df.to_csv(os.path.join(BASE, "results", "method_comparison.csv"), index=False)
print("\n" + "═" * 70)
print("METHOD COMPARISON SUMMARY")
print("═" * 70)
print(comparison_df.to_string(index=False))

# Save process log
log_event("completion", "run_completed", "co-scientist-causal-inference",
          files=["results/causal_estimates.json", "results/method_comparison.csv",
                 "results/hte_quintile_summary.csv", "results/individual_cate_estimates.csv"])

with open(os.path.join(BASE, "logs", "process-log.jsonl"), "w") as f:
    for entry in LOG_ENTRIES:
        f.write(json.dumps(entry, default=str) + "\n")

print("\n✅ All analyses complete. Results saved to results/ and figures/")
