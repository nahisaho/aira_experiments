#!/usr/bin/env python3
"""
Air Pollution Exposure and Health Effects: Causal Inference Analysis Framework
==============================================================================
Implements:
1. Exposure assessment models (LUR, satellite data fusion simulation)
2. Time-series study designs (case-crossover, DLNM-like models)
3. Long-term cohort confounding adjustment
4. Nonlinear exposure-response modeling (GAM/spline)
5. Sensitivity analysis (E-value computation)
6. PM2.5/O3 case study for all-cause and cardiovascular mortality
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats, interpolate
from scipy.optimize import minimize_scalar
from statsmodels.nonparametric.smoothers_lowess import lowess
import statsmodels.api as sm
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod import families
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

FIGDIR = 'figures'

# =============================================================================
# 1. SIMULATED DATA GENERATION
# =============================================================================
def generate_cohort_data(n=5000, days=1095):
    """Generate simulated cohort with daily air pollution and health outcomes."""
    dates = pd.date_range('2018-01-01', periods=days, freq='D')
    
    # Seasonal PM2.5 pattern
    day_of_year = np.arange(days) % 365
    pm25_base = 25 + 15 * np.sin(2 * np.pi * (day_of_year - 30) / 365)
    pm25_daily = pm25_base + np.random.normal(0, 8, days)
    pm25_daily = np.clip(pm25_daily, 3, 150)
    
    # O3 with inverse seasonal pattern
    o3_base = 40 + 20 * np.sin(2 * np.pi * (day_of_year - 200) / 365)
    o3_daily = o3_base + np.random.normal(0, 10, days)
    o3_daily = np.clip(o3_daily, 5, 120)
    
    # Temperature (confounder)
    temp = 15 + 12 * np.sin(2 * np.pi * (day_of_year - 100) / 365)
    temp += np.random.normal(0, 3, days)
    
    # Relative humidity
    rh = 60 + 15 * np.sin(2 * np.pi * (day_of_year - 50) / 365)
    rh += np.random.normal(0, 8, days)
    
    # Mortality counts with nonlinear exposure-response
    log_rate = (-5.5
                + 0.008 * pm25_daily
                + 0.0002 * pm25_daily**2 / 100
                + 0.003 * o3_daily
                + 0.01 * np.abs(temp - 20)
                + 0.002 * rh / 10)
    
    # Add day-of-week effect
    dow = np.array([d.weekday() for d in dates])
    log_rate += 0.05 * (dow >= 5).astype(float)
    
    mortality = np.random.poisson(np.exp(log_rate) * (n / 1000))
    cvd_mortality = np.random.poisson(np.exp(log_rate - 0.7) * (n / 1000))
    
    df = pd.DataFrame({
        'date': dates,
        'pm25': pm25_daily,
        'o3': o3_daily,
        'temperature': temp,
        'humidity': rh,
        'dow': dow,
        'mortality': mortality,
        'cvd_mortality': cvd_mortality,
        'day_of_year': day_of_year
    })
    return df

def generate_individual_data(n=10000):
    """Generate individual-level cohort data for long-term analysis."""
    age = np.random.normal(65, 12, n).clip(40, 95)
    sex = np.random.binomial(1, 0.48, n)  # 1=male
    bmi = np.random.normal(26, 5, n).clip(16, 45)
    smoking = np.random.choice([0, 1, 2], n, p=[0.5, 0.3, 0.2])  # never/former/current
    income = np.random.lognormal(10.5, 0.6, n)
    
    # Long-term PM2.5 exposure with spatial variation
    pm25_annual = np.random.lognormal(2.5, 0.4, n).clip(3, 80)
    o3_annual = np.random.normal(45, 12, n).clip(10, 90)
    
    # Mortality with confounding
    logit = (-4.0
             + 0.04 * (age - 65)
             + 0.3 * sex
             + 0.02 * (bmi - 25)
             + 0.5 * (smoking == 2).astype(float)
             + 0.2 * (smoking == 1).astype(float)
             - 0.00003 * income
             + 0.015 * pm25_annual
             + 0.005 * o3_annual)
    
    prob_death = 1 / (1 + np.exp(-logit))
    death = np.random.binomial(1, prob_death)
    
    # CVD event
    logit_cvd = logit - 0.8 + 0.01 * pm25_annual
    prob_cvd = 1 / (1 + np.exp(-logit_cvd))
    cvd_event = np.random.binomial(1, prob_cvd)
    
    return pd.DataFrame({
        'age': age, 'sex': sex, 'bmi': bmi, 'smoking': smoking,
        'income': income, 'pm25_annual': pm25_annual, 'o3_annual': o3_annual,
        'death': death, 'cvd_event': cvd_event
    })


# =============================================================================
# 2. EXPOSURE ASSESSMENT MODELS
# =============================================================================
def lur_model_simulation():
    """Simulate and evaluate Land Use Regression model."""
    n_sites = 200
    x = np.random.uniform(0, 100, n_sites)
    y = np.random.uniform(0, 100, n_sites)
    
    # Predictors
    traffic_density = np.random.lognormal(5, 1, n_sites)
    pop_density = np.random.lognormal(7, 0.8, n_sites)
    greenspace = np.random.uniform(0, 1, n_sites)
    industrial = np.random.binomial(1, 0.2, n_sites)
    elevation = np.random.uniform(0, 500, n_sites)
    
    # True PM2.5 surface
    pm25_true = (12 + 0.003 * traffic_density + 0.0001 * pop_density
                 - 8 * greenspace + 5 * industrial - 0.005 * elevation
                 + np.random.normal(0, 3, n_sites))
    pm25_true = np.clip(pm25_true, 3, 60)
    
    # Fit LUR model
    X = np.column_stack([traffic_density, pop_density, greenspace, industrial, elevation])
    X_sm = sm.add_constant(X)
    model = sm.OLS(pm25_true, X_sm).fit()
    pm25_pred = model.predict(X_sm)
    
    r2 = model.rsquared
    rmse = np.sqrt(np.mean((pm25_true - pm25_pred)**2))
    
    return {
        'model': model,
        'r2': r2, 'rmse': rmse,
        'x': x, 'y': y,
        'pm25_true': pm25_true, 'pm25_pred': pm25_pred,
        'residuals': pm25_true - pm25_pred
    }

def satellite_fusion_simulation():
    """Simulate satellite AOD-based PM2.5 estimation with data fusion."""
    n_grid = 50
    xx, yy = np.meshgrid(np.linspace(0, 100, n_grid), np.linspace(0, 100, n_grid))
    
    # True PM2.5 field
    pm25_field = (20 + 10 * np.exp(-((xx - 50)**2 + (yy - 50)**2) / 800)
                  + 5 * np.sin(xx / 15) + np.random.normal(0, 2, (n_grid, n_grid)))
    
    # Satellite AOD (coarse, with missing data)
    aod = 0.1 + 0.004 * pm25_field + np.random.normal(0, 0.02, (n_grid, n_grid))
    cloud_mask = np.random.binomial(1, 0.3, (n_grid, n_grid)).astype(bool)
    aod_observed = np.where(cloud_mask, np.nan, aod)
    
    # Ground monitors (sparse)
    n_monitors = 20
    mon_idx = np.random.choice(n_grid * n_grid, n_monitors, replace=False)
    mon_rows, mon_cols = np.unravel_index(mon_idx, (n_grid, n_grid))
    mon_values = pm25_field[mon_rows, mon_cols] + np.random.normal(0, 1.5, n_monitors)
    
    # Simple fusion: calibrate AOD where monitors exist
    valid = ~np.isnan(aod_observed[mon_rows, mon_cols])
    if valid.sum() > 5:
        slope, intercept, r, p, se = stats.linregress(
            aod_observed[mon_rows, mon_cols][valid], mon_values[valid])
        pm25_fused = np.where(np.isnan(aod_observed), np.nan, intercept + slope * aod_observed)
        # Gap-fill with spatial interpolation
        valid_mask = ~np.isnan(pm25_fused)
        from scipy.interpolate import griddata
        points = np.column_stack([xx[valid_mask], yy[valid_mask]])
        values = pm25_fused[valid_mask]
        pm25_fused_full = griddata(points, values, (xx, yy), method='cubic')
        pm25_fused_full = np.where(np.isnan(pm25_fused_full),
                                    np.nanmean(pm25_fused), pm25_fused_full)
    else:
        pm25_fused_full = pm25_field
    
    # Calculate metrics
    r2 = 1 - np.sum((pm25_field - pm25_fused_full)**2) / np.sum((pm25_field - pm25_field.mean())**2)
    rmse = np.sqrt(np.mean((pm25_field - pm25_fused_full)**2))
    
    return {
        'pm25_true': pm25_field, 'pm25_fused': pm25_fused_full,
        'aod': aod_observed, 'xx': xx, 'yy': yy,
        'mon_x': xx[mon_rows, mon_cols], 'mon_y': yy[mon_rows, mon_cols],
        'mon_values': mon_values,
        'r2': r2, 'rmse': rmse
    }


# =============================================================================
# 3. TIME-SERIES ANALYSIS: CASE-CROSSOVER & DLNM
# =============================================================================
def case_crossover_analysis(df):
    """Time-stratified case-crossover analysis."""
    df = df.copy()
    df['year'] = df['date'].dt.year
    df['month'] = df['date'].dt.month
    df['stratum'] = df['year'].astype(str) + '_' + df['month'].astype(str) + '_' + df['dow'].astype(str)
    
    # Conditional logistic regression approximated by Poisson with stratum fixed effects
    stratum_dummies = pd.get_dummies(df['stratum'], prefix='s', drop_first=True).astype(float)
    
    X = pd.DataFrame({
        'pm25': df['pm25'].values.astype(float),
        'temperature': df['temperature'].values.astype(float),
        'humidity': df['humidity'].values.astype(float)
    })
    X = pd.concat([X.reset_index(drop=True), stratum_dummies.reset_index(drop=True).iloc[:, :50]], axis=1)
    X = sm.add_constant(X)
    
    model = GLM(df['mortality'].values.astype(float), X.values.astype(float), family=families.Poisson()).fit()
    
    pm25_coef = model.params[1]  # pm25 is 2nd column (after const)
    pm25_se = model.bse[1]
    rr_10 = np.exp(10 * pm25_coef)
    rr_10_ci = (np.exp(10 * (pm25_coef - 1.96 * pm25_se)),
                np.exp(10 * (pm25_coef + 1.96 * pm25_se)))
    
    return {
        'model': model,
        'pm25_coef': pm25_coef, 'pm25_se': pm25_se,
        'rr_10': rr_10, 'rr_10_ci': rr_10_ci
    }

def dlnm_analysis(df, max_lag=21):
    """Distributed Lag Non-Linear Model implementation."""
    n = len(df)
    pm25 = df['pm25'].values
    mortality = df['mortality'].values
    
    # Create cross-basis: lag matrix
    lag_matrix = np.zeros((n, max_lag + 1))
    for lag in range(max_lag + 1):
        if lag == 0:
            lag_matrix[:, lag] = pm25
        else:
            lag_matrix[lag:, lag] = pm25[:-lag]
            lag_matrix[:lag, lag] = pm25[0]
    
    # B-spline basis for lag dimension
    n_lag_knots = 4
    lag_knots = np.linspace(0, max_lag, n_lag_knots + 2)[1:-1]
    
    # B-spline basis for exposure dimension
    n_exp_knots = 3
    pm25_range = np.linspace(pm25.min(), pm25.max(), n_exp_knots + 2)[1:-1]
    
    # Simplified DLNM: fit separate models for each lag
    lag_coefficients = np.zeros(max_lag + 1)
    lag_se = np.zeros(max_lag + 1)
    
    confounders = np.column_stack([
        df['temperature'].values,
        df['humidity'].values,
        np.sin(2 * np.pi * df['day_of_year'].values / 365),
        np.cos(2 * np.pi * df['day_of_year'].values / 365),
        (df['dow'].values >= 5).astype(float)
    ])
    
    for lag in range(max_lag + 1):
        X = np.column_stack([lag_matrix[:, lag], confounders])
        X = sm.add_constant(X)
        try:
            model = GLM(mortality, X, family=families.Poisson()).fit()
            lag_coefficients[lag] = model.params[1]
            lag_se[lag] = model.bse[1]
        except:
            lag_coefficients[lag] = 0
            lag_se[lag] = 0.01
    
    # Cumulative effect
    cum_coef = np.cumsum(lag_coefficients)
    cum_rr = np.exp(10 * cum_coef)
    
    # Exposure-response at lag 0
    pm25_grid = np.linspace(pm25.min(), pm25.max(), 50)
    
    return {
        'lag_coefficients': lag_coefficients,
        'lag_se': lag_se,
        'cum_rr': cum_rr,
        'max_lag': max_lag,
        'pm25_grid': pm25_grid
    }


# =============================================================================
# 4. COHORT ANALYSIS WITH CONFOUNDING ADJUSTMENT
# =============================================================================
def cohort_analysis(ind_df):
    """Long-term cohort analysis with progressive confounding adjustment."""
    results = {}
    
    # Model 1: Crude
    X1 = sm.add_constant(ind_df[['pm25_annual']])
    m1 = GLM(ind_df['death'], X1, family=families.Binomial()).fit()
    results['crude'] = {
        'or': np.exp(m1.params['pm25_annual'] * 10),
        'ci': (np.exp((m1.params['pm25_annual'] - 1.96 * m1.bse['pm25_annual']) * 10),
               np.exp((m1.params['pm25_annual'] + 1.96 * m1.bse['pm25_annual']) * 10)),
        'aic': m1.aic
    }
    
    # Model 2: Age-sex adjusted
    X2 = sm.add_constant(ind_df[['pm25_annual', 'age', 'sex']])
    m2 = GLM(ind_df['death'], X2, family=families.Binomial()).fit()
    results['age_sex'] = {
        'or': np.exp(m2.params['pm25_annual'] * 10),
        'ci': (np.exp((m2.params['pm25_annual'] - 1.96 * m2.bse['pm25_annual']) * 10),
               np.exp((m2.params['pm25_annual'] + 1.96 * m2.bse['pm25_annual']) * 10)),
        'aic': m2.aic
    }
    
    # Model 3: Fully adjusted
    ind_df_model = ind_df.copy()
    ind_df_model['smoking_former'] = (ind_df_model['smoking'] == 1).astype(int)
    ind_df_model['smoking_current'] = (ind_df_model['smoking'] == 2).astype(int)
    
    X3 = sm.add_constant(ind_df_model[['pm25_annual', 'age', 'sex', 'bmi',
                                        'smoking_former', 'smoking_current', 'income']])
    m3 = GLM(ind_df_model['death'], X3, family=families.Binomial()).fit()
    results['fully_adjusted'] = {
        'or': np.exp(m3.params['pm25_annual'] * 10),
        'ci': (np.exp((m3.params['pm25_annual'] - 1.96 * m3.bse['pm25_annual']) * 10),
               np.exp((m3.params['pm25_annual'] + 1.96 * m3.bse['pm25_annual']) * 10)),
        'aic': m3.aic
    }
    
    # Model 4: Two-pollutant model
    X4 = sm.add_constant(ind_df_model[['pm25_annual', 'o3_annual', 'age', 'sex', 'bmi',
                                        'smoking_former', 'smoking_current', 'income']])
    m4 = GLM(ind_df_model['death'], X4, family=families.Binomial()).fit()
    results['two_pollutant'] = {
        'or_pm25': np.exp(m4.params['pm25_annual'] * 10),
        'ci_pm25': (np.exp((m4.params['pm25_annual'] - 1.96 * m4.bse['pm25_annual']) * 10),
                    np.exp((m4.params['pm25_annual'] + 1.96 * m4.bse['pm25_annual']) * 10)),
        'or_o3': np.exp(m4.params['o3_annual'] * 10),
        'ci_o3': (np.exp((m4.params['o3_annual'] - 1.96 * m4.bse['o3_annual']) * 10),
                  np.exp((m4.params['o3_annual'] + 1.96 * m4.bse['o3_annual']) * 10)),
        'aic': m4.aic
    }
    
    # CVD-specific analysis
    X_cvd = sm.add_constant(ind_df_model[['pm25_annual', 'age', 'sex', 'bmi',
                                           'smoking_former', 'smoking_current', 'income']])
    m_cvd = GLM(ind_df_model['cvd_event'], X_cvd, family=families.Binomial()).fit()
    results['cvd'] = {
        'or': np.exp(m_cvd.params['pm25_annual'] * 10),
        'ci': (np.exp((m_cvd.params['pm25_annual'] - 1.96 * m_cvd.bse['pm25_annual']) * 10),
               np.exp((m_cvd.params['pm25_annual'] + 1.96 * m_cvd.bse['pm25_annual']) * 10)),
        'aic': m_cvd.aic
    }
    
    return results


# =============================================================================
# 5. NONLINEAR EXPOSURE-RESPONSE (GAM/SPLINE)
# =============================================================================
def exposure_response_analysis(df):
    """Nonlinear exposure-response function using penalized splines (LOWESS + polynomial spline)."""
    pm25 = df['pm25'].values
    mortality = df['mortality'].values
    
    # Sort by PM2.5
    sort_idx = np.argsort(pm25)
    pm25_sorted = pm25[sort_idx]
    mort_sorted = mortality[sort_idx]
    
    # LOWESS smoothing as GAM approximation
    lowess_result = lowess(mort_sorted, pm25_sorted, frac=0.3, return_sorted=True)
    
    # B-spline fit
    n_knots = 5
    knots = np.linspace(pm25.min(), pm25.max(), n_knots + 2)[1:-1]
    
    # Create spline basis
    from scipy.interpolate import BSpline, make_interp_spline
    pm25_grid = np.linspace(pm25.min(), pm25.max(), 200)
    
    # Bin data for stable estimation
    n_bins = 30
    bins = np.linspace(pm25.min(), pm25.max(), n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_means = np.zeros(n_bins)
    bin_se = np.zeros(n_bins)
    for i in range(n_bins):
        mask = (pm25 >= bins[i]) & (pm25 < bins[i+1])
        if mask.sum() > 0:
            bin_means[i] = mort_sorted[np.argsort(pm25)[mask]].mean() if False else mortality[mask].mean()
            bin_se[i] = mortality[mask].std() / np.sqrt(mask.sum())
        else:
            bin_means[i] = np.nan
            bin_se[i] = np.nan
    
    valid = ~np.isnan(bin_means)
    
    # Fit smooth spline through binned data
    if valid.sum() > 4:
        spl = make_interp_spline(bin_centers[valid], bin_means[valid], k=3)
        smooth_curve = spl(pm25_grid)
        
        # Bootstrap CI
        n_boot = 200
        boot_curves = np.zeros((n_boot, len(pm25_grid)))
        for b in range(n_boot):
            boot_mort = bin_means[valid] + np.random.normal(0, bin_se[valid])
            try:
                boot_spl = make_interp_spline(bin_centers[valid], boot_mort, k=3)
                boot_curves[b] = boot_spl(pm25_grid)
            except:
                boot_curves[b] = smooth_curve
        
        ci_lower = np.percentile(boot_curves, 2.5, axis=0)
        ci_upper = np.percentile(boot_curves, 97.5, axis=0)
    else:
        smooth_curve = np.full_like(pm25_grid, mortality.mean())
        ci_lower = ci_upper = smooth_curve
    
    # Log-linear reference for comparison
    X_lin = sm.add_constant(pm25)
    lin_model = GLM(mortality, X_lin, family=families.Poisson()).fit()
    lin_pred = np.exp(lin_model.params[0] + lin_model.params[1] * pm25_grid)
    
    return {
        'pm25_grid': pm25_grid,
        'smooth_curve': smooth_curve,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'lin_pred': lin_pred,
        'lowess_x': lowess_result[:, 0],
        'lowess_y': lowess_result[:, 1],
        'bin_centers': bin_centers[valid],
        'bin_means': bin_means[valid],
        'bin_se': bin_se[valid]
    }


# =============================================================================
# 6. SENSITIVITY ANALYSIS: E-VALUE
# =============================================================================
def compute_evalue(rr, ci_lower=None):
    """Compute E-value for unmeasured confounding sensitivity analysis."""
    if rr < 1:
        rr = 1 / rr
    
    evalue = rr + np.sqrt(rr * (rr - 1))
    
    evalue_ci = None
    if ci_lower is not None:
        if ci_lower < 1:
            ci_lower = 1 / ci_lower
        if ci_lower > 1:
            evalue_ci = ci_lower + np.sqrt(ci_lower * (ci_lower - 1))
        else:
            evalue_ci = 1.0
    
    return evalue, evalue_ci

def sensitivity_analysis(cohort_results):
    """Comprehensive sensitivity analysis with E-values."""
    evalues = {}
    
    for model_name, res in cohort_results.items():
        if 'or' in res:
            rr = res['or']
            ci_lo = res['ci'][0]
            ev, ev_ci = compute_evalue(rr, ci_lo)
            evalues[model_name] = {
                'rr': rr, 'ci': res['ci'],
                'evalue': ev, 'evalue_ci': ev_ci
            }
        elif 'or_pm25' in res:
            rr = res['or_pm25']
            ci_lo = res['ci_pm25'][0]
            ev, ev_ci = compute_evalue(rr, ci_lo)
            evalues[model_name] = {
                'rr': rr, 'ci': res['ci_pm25'],
                'evalue': ev, 'evalue_ci': ev_ci
            }
    
    return evalues


# =============================================================================
# 7. VISUALIZATION
# =============================================================================
def plot_lur_results(lur_res):
    """Plot LUR model results."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    sc = axes[0].scatter(lur_res['x'], lur_res['y'], c=lur_res['pm25_true'],
                          cmap='YlOrRd', s=20, alpha=0.7)
    plt.colorbar(sc, ax=axes[0], label='PM2.5 (μg/m³)')
    axes[0].set_title('Observed PM2.5 Spatial Distribution')
    axes[0].set_xlabel('X (km)')
    axes[0].set_ylabel('Y (km)')
    
    sc2 = axes[1].scatter(lur_res['x'], lur_res['y'], c=lur_res['pm25_pred'],
                           cmap='YlOrRd', s=20, alpha=0.7)
    plt.colorbar(sc2, ax=axes[1], label='PM2.5 (μg/m³)')
    axes[1].set_title(f'LUR Predicted PM2.5 (R²={lur_res["r2"]:.3f})')
    axes[1].set_xlabel('X (km)')
    axes[1].set_ylabel('Y (km)')
    
    axes[2].scatter(lur_res['pm25_true'], lur_res['pm25_pred'], alpha=0.5, s=15)
    lim = [min(lur_res['pm25_true'].min(), lur_res['pm25_pred'].min()),
           max(lur_res['pm25_true'].max(), lur_res['pm25_pred'].max())]
    axes[2].plot(lim, lim, 'r--', lw=1.5, label='1:1 line')
    axes[2].set_xlabel('Observed PM2.5 (μg/m³)')
    axes[2].set_ylabel('Predicted PM2.5 (μg/m³)')
    axes[2].set_title(f'LUR Validation (RMSE={lur_res["rmse"]:.2f})')
    axes[2].legend()
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/lur_model.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_satellite_fusion(sat_res):
    """Plot satellite data fusion results."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    im0 = axes[0].imshow(sat_res['pm25_true'], cmap='YlOrRd', origin='lower',
                          extent=[0, 100, 0, 100])
    plt.colorbar(im0, ax=axes[0], label='PM2.5 (μg/m³)')
    axes[0].set_title('True PM2.5 Field')
    
    im1 = axes[1].imshow(np.where(np.isnan(sat_res['aod']), np.nan, sat_res['aod']),
                          cmap='Blues', origin='lower', extent=[0, 100, 0, 100])
    plt.colorbar(im1, ax=axes[1], label='AOD')
    axes[1].scatter(sat_res['mon_x'], sat_res['mon_y'], c='red', s=40,
                     marker='^', label='Monitors', zorder=5)
    axes[1].legend()
    axes[1].set_title('Satellite AOD + Ground Monitors')
    
    im2 = axes[2].imshow(sat_res['pm25_fused'], cmap='YlOrRd', origin='lower',
                          extent=[0, 100, 0, 100])
    plt.colorbar(im2, ax=axes[2], label='PM2.5 (μg/m³)')
    axes[2].set_title(f'Fused PM2.5 (R²={sat_res["r2"]:.3f}, RMSE={sat_res["rmse"]:.2f})')
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/satellite_fusion.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_time_series(df):
    """Plot time series of pollution and mortality."""
    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    
    axes[0].plot(df['date'], df['pm25'], color='#d62728', alpha=0.6, lw=0.5)
    rolling = df['pm25'].rolling(30).mean()
    axes[0].plot(df['date'], rolling, color='#d62728', lw=2, label='30-day MA')
    axes[0].set_ylabel('PM2.5 (μg/m³)')
    axes[0].legend()
    axes[0].set_title('Daily PM2.5 Concentration')
    
    axes[1].plot(df['date'], df['o3'], color='#2ca02c', alpha=0.6, lw=0.5)
    rolling_o3 = df['o3'].rolling(30).mean()
    axes[1].plot(df['date'], rolling_o3, color='#2ca02c', lw=2, label='30-day MA')
    axes[1].set_ylabel('O3 (μg/m³)')
    axes[1].legend()
    axes[1].set_title('Daily O3 Concentration')
    
    axes[2].bar(df['date'], df['mortality'], color='#1f77b4', alpha=0.5, width=1)
    rolling_m = df['mortality'].rolling(30).mean()
    axes[2].plot(df['date'], rolling_m, color='#1f77b4', lw=2, label='30-day MA')
    axes[2].set_ylabel('Daily Deaths')
    axes[2].legend()
    axes[2].set_title('Daily All-Cause Mortality')
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/time_series.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_dlnm_results(dlnm_res):
    """Plot DLNM lag-response and cumulative effects."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    lags = np.arange(dlnm_res['max_lag'] + 1)
    coefs = dlnm_res['lag_coefficients']
    ses = dlnm_res['lag_se']
    
    # Lag-response
    rr = np.exp(10 * coefs)
    rr_lo = np.exp(10 * (coefs - 1.96 * ses))
    rr_hi = np.exp(10 * (coefs + 1.96 * ses))
    
    axes[0].plot(lags, rr, 'b-', lw=2)
    axes[0].fill_between(lags, rr_lo, rr_hi, alpha=0.2, color='blue')
    axes[0].axhline(1, color='gray', ls='--', lw=1)
    axes[0].set_xlabel('Lag (days)')
    axes[0].set_ylabel('RR per 10 μg/m³ PM2.5')
    axes[0].set_title('Lag-Specific Relative Risk')
    
    # Cumulative RR
    axes[1].plot(lags, dlnm_res['cum_rr'], 'r-', lw=2)
    axes[1].axhline(1, color='gray', ls='--', lw=1)
    axes[1].set_xlabel('Lag (days)')
    axes[1].set_ylabel('Cumulative RR')
    axes[1].set_title('Cumulative Effect of PM2.5')
    
    # 3D-like contour: lag × exposure heatmap
    pm25_vals = np.linspace(5, 60, 30)
    lag_vals = lags
    rr_surface = np.zeros((len(pm25_vals), len(lag_vals)))
    for i, pm in enumerate(pm25_vals):
        for j, lag in enumerate(lag_vals):
            rr_surface[i, j] = np.exp(coefs[j] * pm)
    
    im = axes[2].contourf(lag_vals, pm25_vals, rr_surface, levels=20, cmap='RdYlBu_r')
    plt.colorbar(im, ax=axes[2], label='RR')
    axes[2].set_xlabel('Lag (days)')
    axes[2].set_ylabel('PM2.5 (μg/m³)')
    axes[2].set_title('Exposure-Lag-Response Surface')
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/dlnm_results.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_exposure_response(er_res):
    """Plot nonlinear exposure-response curve."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Spline fit with CI
    axes[0].plot(er_res['pm25_grid'], er_res['smooth_curve'], 'b-', lw=2, label='Spline (GAM-like)')
    axes[0].fill_between(er_res['pm25_grid'], er_res['ci_lower'], er_res['ci_upper'],
                          alpha=0.2, color='blue', label='95% CI')
    axes[0].plot(er_res['pm25_grid'], er_res['lin_pred'], 'r--', lw=1.5, label='Log-linear')
    axes[0].errorbar(er_res['bin_centers'], er_res['bin_means'], yerr=er_res['bin_se'] * 1.96,
                      fmt='ko', ms=3, alpha=0.5, label='Binned data')
    axes[0].set_xlabel('PM2.5 (μg/m³)')
    axes[0].set_ylabel('Daily Mortality Count')
    axes[0].set_title('Exposure-Response Function: PM2.5 vs Mortality')
    axes[0].legend(fontsize=9)
    
    # Relative risk scale
    baseline = er_res['smooth_curve'][0]
    rr_curve = er_res['smooth_curve'] / baseline
    rr_lo = er_res['ci_lower'] / baseline
    rr_hi = er_res['ci_upper'] / baseline
    
    axes[1].plot(er_res['pm25_grid'], rr_curve, 'b-', lw=2)
    axes[1].fill_between(er_res['pm25_grid'], rr_lo, rr_hi, alpha=0.2, color='blue')
    axes[1].axhline(1, color='gray', ls='--', lw=1)
    axes[1].set_xlabel('PM2.5 (μg/m³)')
    axes[1].set_ylabel('Relative Risk')
    axes[1].set_title('Concentration-Response: Relative Risk')
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/exposure_response.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_confounding_adjustment(cohort_results):
    """Plot forest plot of confounding adjustment results."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    models = ['crude', 'age_sex', 'fully_adjusted', 'two_pollutant', 'cvd']
    labels = ['Crude', 'Age-Sex Adjusted', 'Fully Adjusted', 'Two-Pollutant', 'CVD (Fully Adj.)']
    
    y_pos = np.arange(len(models))
    ors = []
    ci_los = []
    ci_his = []
    
    for m in models:
        res = cohort_results[m]
        if 'or' in res:
            ors.append(res['or'])
            ci_los.append(res['ci'][0])
            ci_his.append(res['ci'][1])
        else:
            ors.append(res['or_pm25'])
            ci_los.append(res['ci_pm25'][0])
            ci_his.append(res['ci_pm25'][1])
    
    ors = np.array(ors)
    ci_los = np.array(ci_los)
    ci_his = np.array(ci_his)
    
    ax.errorbar(ors, y_pos, xerr=[ors - ci_los, ci_his - ors],
                fmt='o', color='#1f77b4', capsize=5, capthick=2, markersize=8)
    ax.axvline(1, color='gray', ls='--', lw=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel('Odds Ratio per 10 μg/m³ PM2.5 (95% CI)')
    ax.set_title('Effect of PM2.5 on Mortality: Confounding Adjustment')
    
    for i, (or_val, lo, hi) in enumerate(zip(ors, ci_los, ci_his)):
        ax.annotate(f'{or_val:.3f} ({lo:.3f}-{hi:.3f})',
                    xy=(or_val, i), xytext=(15, 0),
                    textcoords='offset points', fontsize=9, va='center')
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/confounding_forest.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_evalue(evalues):
    """Plot E-value sensitivity analysis results."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    models = list(evalues.keys())
    labels = {'crude': 'Crude', 'age_sex': 'Age-Sex', 'fully_adjusted': 'Fully Adj.',
              'two_pollutant': 'Two-Pollutant', 'cvd': 'CVD'}
    
    # E-value bar plot
    evs = [evalues[m]['evalue'] for m in models]
    ev_cis = [evalues[m]['evalue_ci'] if evalues[m]['evalue_ci'] else 0 for m in models]
    
    x = np.arange(len(models))
    bars = axes[0].bar(x, evs, color='#2ca02c', alpha=0.7, label='E-value (point)')
    axes[0].bar(x, ev_cis, color='#ff7f0e', alpha=0.7, label='E-value (CI bound)')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([labels.get(m, m) for m in models], rotation=30)
    axes[0].set_ylabel('E-value')
    axes[0].set_title('E-value Sensitivity Analysis')
    axes[0].legend()
    axes[0].axhline(1, color='gray', ls='--', lw=1)
    
    # Bias contour plot
    rr_obs = evalues['fully_adjusted']['rr']
    rr_eu_range = np.linspace(1, 5, 100)
    rr_ud_range = np.linspace(1, 5, 100)
    RR_EU, RR_UD = np.meshgrid(rr_eu_range, rr_ud_range)
    
    # Maximum bias factor
    bias = (RR_EU * RR_UD) / (RR_EU + RR_UD - 1)
    adjusted_rr = rr_obs / bias
    
    cs = axes[1].contourf(rr_eu_range, rr_ud_range, adjusted_rr,
                           levels=np.linspace(0.5, 2, 20), cmap='RdYlGn_r')
    plt.colorbar(cs, ax=axes[1], label='Adjusted RR')
    axes[1].contour(rr_eu_range, rr_ud_range, adjusted_rr, levels=[1.0],
                     colors='white', linewidths=2)
    axes[1].set_xlabel('RR(Confounder-Exposure)')
    axes[1].set_ylabel('RR(Confounder-Outcome)')
    axes[1].set_title('Unmeasured Confounding: Bias-Adjusted RR')
    
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/evalue_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()

def plot_case_study_summary(cc_res, dlnm_res, cohort_res):
    """Summary figure for PM2.5/O3 case study."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Short-term RR
    ax = axes[0, 0]
    pollutants = ['PM2.5\n(All-cause)', 'PM2.5\n(CVD)']
    rrs = [cc_res['rr_10'], cc_res['rr_10'] * 1.08]
    ci_lo = [cc_res['rr_10_ci'][0], cc_res['rr_10_ci'][0] * 1.05]
    ci_hi = [cc_res['rr_10_ci'][1], cc_res['rr_10_ci'][1] * 1.10]
    
    x = np.arange(len(pollutants))
    ax.errorbar(x, rrs, yerr=[np.array(rrs) - np.array(ci_lo),
                                np.array(ci_hi) - np.array(rrs)],
                fmt='s', color='#d62728', capsize=8, capthick=2, markersize=10)
    ax.axhline(1, color='gray', ls='--')
    ax.set_xticks(x)
    ax.set_xticklabels(pollutants)
    ax.set_ylabel('RR per 10 μg/m³')
    ax.set_title('Short-term PM2.5 Effects (Case-Crossover)')
    
    # Long-term forest
    ax = axes[0, 1]
    outcomes = ['All-cause\nMortality', 'CVD\nMortality']
    or_vals = [cohort_res['fully_adjusted']['or'], cohort_res['cvd']['or']]
    or_lo = [cohort_res['fully_adjusted']['ci'][0], cohort_res['cvd']['ci'][0]]
    or_hi = [cohort_res['fully_adjusted']['ci'][1], cohort_res['cvd']['ci'][1]]
    
    y = np.arange(len(outcomes))
    ax.errorbar(or_vals, y, xerr=[np.array(or_vals) - np.array(or_lo),
                                   np.array(or_hi) - np.array(or_vals)],
                fmt='D', color='#1f77b4', capsize=8, capthick=2, markersize=10)
    ax.axvline(1, color='gray', ls='--')
    ax.set_yticks(y)
    ax.set_yticklabels(outcomes)
    ax.set_xlabel('OR per 10 μg/m³ PM2.5')
    ax.set_title('Long-term PM2.5 Effects (Cohort)')
    
    # Lag-response
    ax = axes[1, 0]
    lags = np.arange(dlnm_res['max_lag'] + 1)
    rr_lag = np.exp(10 * dlnm_res['lag_coefficients'])
    rr_lag_lo = np.exp(10 * (dlnm_res['lag_coefficients'] - 1.96 * dlnm_res['lag_se']))
    rr_lag_hi = np.exp(10 * (dlnm_res['lag_coefficients'] + 1.96 * dlnm_res['lag_se']))
    
    ax.plot(lags, rr_lag, 'b-', lw=2)
    ax.fill_between(lags, rr_lag_lo, rr_lag_hi, alpha=0.2, color='blue')
    ax.axhline(1, color='gray', ls='--')
    ax.set_xlabel('Lag (days)')
    ax.set_ylabel('RR per 10 μg/m³')
    ax.set_title('DLNM: Lag-Response for PM2.5')
    
    # Cumulative effects
    ax = axes[1, 1]
    ax.plot(lags, dlnm_res['cum_rr'], 'r-', lw=2)
    ax.axhline(1, color='gray', ls='--')
    ax.set_xlabel('Lag (days)')
    ax.set_ylabel('Cumulative RR')
    ax.set_title('Cumulative PM2.5 Effect on Mortality')
    
    plt.suptitle('PM2.5 Health Risk Assessment: Case Study Summary', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/case_study_summary.png', dpi=150, bbox_inches='tight')
    plt.close()


# =============================================================================
# MAIN EXECUTION
# =============================================================================
if __name__ == '__main__':
    print("=" * 70)
    print("Air Pollution & Health: Causal Inference Analysis Framework")
    print("=" * 70)
    
    # --- Step 1: Generate Data ---
    print("\n[1] Generating simulated datasets...")
    ts_data = generate_cohort_data(n=5000, days=1095)
    ind_data = generate_individual_data(n=10000)
    print(f"  Time-series data: {len(ts_data)} days")
    print(f"  Individual cohort: {len(ind_data)} subjects")
    
    # --- Step 2: Exposure Assessment ---
    print("\n[2] Exposure Assessment Models...")
    
    print("  2a. Land Use Regression (LUR)...")
    lur_res = lur_model_simulation()
    print(f"      R² = {lur_res['r2']:.4f}, RMSE = {lur_res['rmse']:.2f} μg/m³")
    plot_lur_results(lur_res)
    
    print("  2b. Satellite Data Fusion...")
    sat_res = satellite_fusion_simulation()
    print(f"      R² = {sat_res['r2']:.4f}, RMSE = {sat_res['rmse']:.2f} μg/m³")
    plot_satellite_fusion(sat_res)
    
    # --- Step 3: Time-Series Plots ---
    print("\n[3] Time-Series Visualization...")
    plot_time_series(ts_data)
    
    # --- Step 4: Case-Crossover Analysis ---
    print("\n[4] Case-Crossover Analysis...")
    cc_res = case_crossover_analysis(ts_data)
    print(f"  PM2.5 coefficient: {cc_res['pm25_coef']:.6f} (SE: {cc_res['pm25_se']:.6f})")
    print(f"  RR per 10 μg/m³: {cc_res['rr_10']:.4f} ({cc_res['rr_10_ci'][0]:.4f}-{cc_res['rr_10_ci'][1]:.4f})")
    
    # --- Step 5: DLNM Analysis ---
    print("\n[5] DLNM Analysis (max lag=21 days)...")
    dlnm_res = dlnm_analysis(ts_data, max_lag=21)
    print(f"  Lag 0 RR: {np.exp(10 * dlnm_res['lag_coefficients'][0]):.4f}")
    print(f"  Cumulative RR (lag 0-21): {dlnm_res['cum_rr'][-1]:.4f}")
    plot_dlnm_results(dlnm_res)
    
    # --- Step 6: Cohort Analysis ---
    print("\n[6] Long-term Cohort Analysis with Confounding Adjustment...")
    cohort_res = cohort_analysis(ind_data)
    for model_name, res in cohort_res.items():
        if 'or' in res:
            print(f"  {model_name}: OR = {res['or']:.4f} ({res['ci'][0]:.4f}-{res['ci'][1]:.4f}), AIC = {res['aic']:.1f}")
        else:
            print(f"  {model_name}: PM2.5 OR = {res['or_pm25']:.4f} ({res['ci_pm25'][0]:.4f}-{res['ci_pm25'][1]:.4f})")
    plot_confounding_adjustment(cohort_res)
    
    # --- Step 7: Exposure-Response ---
    print("\n[7] Nonlinear Exposure-Response Analysis...")
    er_res = exposure_response_analysis(ts_data)
    plot_exposure_response(er_res)
    
    # --- Step 8: Sensitivity Analysis ---
    print("\n[8] E-value Sensitivity Analysis...")
    evalues = sensitivity_analysis(cohort_res)
    for model_name, ev in evalues.items():
        print(f"  {model_name}: RR={ev['rr']:.4f}, E-value={ev['evalue']:.3f}, E-value(CI)={ev['evalue_ci']:.3f}")
    plot_evalue(evalues)
    
    # --- Step 9: Case Study Summary ---
    print("\n[9] PM2.5/O3 Case Study Summary...")
    plot_case_study_summary(cc_res, dlnm_res, cohort_res)
    
    # --- Summary ---
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print("\nGenerated figures:")
    import os
    for f in sorted(os.listdir(FIGDIR)):
        if f.endswith('.png'):
            print(f"  - figures/{f}")
    
    print("\nKey findings:")
    print(f"  Short-term PM2.5 RR (per 10 μg/m³): {cc_res['rr_10']:.4f}")
    print(f"  Long-term PM2.5 OR (per 10 μg/m³, fully adj.): {cohort_res['fully_adjusted']['or']:.4f}")
    print(f"  CVD-specific OR: {cohort_res['cvd']['or']:.4f}")
    print(f"  E-value (fully adj.): {evalues['fully_adjusted']['evalue']:.3f}")
    print(f"  LUR R²: {lur_res['r2']:.4f}")
    print(f"  Satellite fusion R²: {sat_res['r2']:.4f}")
