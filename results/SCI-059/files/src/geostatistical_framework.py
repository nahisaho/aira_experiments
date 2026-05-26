#!/usr/bin/env python3
"""
Geostatistical Framework for Disease Risk Spatial Pattern Analysis and Prediction.

Implements:
1. Log-Gaussian Cox Process (LGCP) for spatial point patterns
2. Bayesian spatial modeling (INLA/SPDE-like approach via Python)
3. Spatial autocorrelation tests (Moran's I, variogram)
4. Ecological study confounding bias mitigation
5. Spatiotemporal knot-based spline prediction
6. Malaria/Dengue risk mapping case study
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
from scipy import stats, spatial, optimize, interpolate
from scipy.spatial.distance import pdist, squareform
from scipy.linalg import cholesky, cho_solve, cho_factor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, RBF, WhiteKernel
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from statsmodels.stats.diagnostic import het_breuschpagan
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
FIGURES_DIR = "figures"

# ============================================================
# PART 1: Data Generation — Synthetic Disease Risk Landscape
# ============================================================

def generate_spatial_domain(n_regions=100, grid_size=50):
    """Generate a synthetic spatial domain mimicking a geographic study area."""
    x = np.linspace(0, 10, grid_size)
    y = np.linspace(0, 10, grid_size)
    xx, yy = np.meshgrid(x, y)
    coords_grid = np.column_stack([xx.ravel(), yy.ravel()])

    # Random observation locations
    obs_x = np.random.uniform(0.5, 9.5, n_regions)
    obs_y = np.random.uniform(0.5, 9.5, n_regions)
    obs_coords = np.column_stack([obs_x, obs_y])

    return xx, yy, coords_grid, obs_coords


def generate_covariates(obs_coords):
    """Generate synthetic environmental covariates (temperature, rainfall, elevation, urbanization)."""
    n = len(obs_coords)
    x, y = obs_coords[:, 0], obs_coords[:, 1]

    temperature = 25 + 3 * np.sin(x * 0.5) + np.random.normal(0, 0.5, n)
    rainfall = 1200 + 200 * np.cos(y * 0.3) + 100 * np.sin(x * 0.4) + np.random.normal(0, 30, n)
    elevation = 500 - 50 * x + 30 * y + np.random.normal(0, 20, n)
    urbanization = 1 / (1 + np.exp(-(x - 5))) + np.random.normal(0, 0.05, n)
    urbanization = np.clip(urbanization, 0, 1)

    return pd.DataFrame({
        'temperature': temperature,
        'rainfall': rainfall,
        'elevation': elevation,
        'urbanization': urbanization
    })


def generate_true_risk_surface(xx, yy):
    """Generate a true underlying risk surface with spatial structure."""
    risk = (
        0.3 * np.exp(-((xx - 3)**2 + (yy - 7)**2) / 4) +
        0.5 * np.exp(-((xx - 7)**2 + (yy - 3)**2) / 3) +
        0.2 * np.exp(-((xx - 5)**2 + (yy - 5)**2) / 6) +
        0.1 * np.sin(xx * 0.8) * np.cos(yy * 0.6)
    )
    risk = (risk - risk.min()) / (risk.max() - risk.min())
    return risk


def generate_disease_cases(obs_coords, covariates, true_risk_grid, xx, yy):
    """Generate disease case counts based on risk surface and covariates."""
    from scipy.interpolate import griddata
    risk_at_obs = griddata(
        np.column_stack([xx.ravel(), yy.ravel()]),
        true_risk_grid.ravel(),
        obs_coords, method='cubic'
    )
    risk_at_obs = np.clip(risk_at_obs, 0.01, 0.99)

    log_rate = (
        -1.0 +
        0.05 * covariates['temperature'].values +
        0.001 * covariates['rainfall'].values -
        0.002 * covariates['elevation'].values +
        1.5 * covariates['urbanization'].values +
        2.0 * risk_at_obs
    )
    rate = np.exp(log_rate)
    population = np.random.randint(5000, 50000, len(obs_coords))
    cases = np.random.poisson(rate * population / 1000)

    return cases, population, risk_at_obs


# ============================================================
# PART 2: Log-Gaussian Cox Process (LGCP)
# ============================================================

def matern_covariance(d, sigma2, kappa, nu=1.5):
    """Matérn covariance function."""
    from scipy.special import kv, gamma as gamma_func
    d = np.asarray(d, dtype=float)
    result = np.zeros_like(d)
    mask = d > 0
    scaled = kappa * d[mask]
    coeff = sigma2 * (2**(1 - nu)) / gamma_func(nu)
    result[mask] = coeff * (scaled**nu) * kv(nu, scaled)
    result[~mask] = sigma2
    return result


def simulate_lgcp(domain_size=10, grid_n=50, sigma2=1.0, kappa=1.5, mu=-2.0):
    """Simulate a Log-Gaussian Cox Process on a grid."""
    x = np.linspace(0, domain_size, grid_n)
    y = np.linspace(0, domain_size, grid_n)
    xx, yy = np.meshgrid(x, y)
    coords = np.column_stack([xx.ravel(), yy.ravel()])

    dist_matrix = squareform(pdist(coords))
    cov_matrix = matern_covariance(dist_matrix, sigma2, kappa)
    cov_matrix += 1e-6 * np.eye(len(coords))

    L = cholesky(cov_matrix, lower=True)
    z = np.random.normal(0, 1, len(coords))
    log_intensity = mu + L @ z

    intensity = np.exp(log_intensity)
    cell_area = (domain_size / grid_n) ** 2
    counts = np.random.poisson(intensity * cell_area)

    point_x, point_y = [], []
    dx = domain_size / grid_n
    for i, (cx, cy) in enumerate(coords):
        for _ in range(counts[i]):
            px = cx + np.random.uniform(-dx/2, dx/2)
            py = cy + np.random.uniform(-dx/2, dx/2)
            point_x.append(px)
            point_y.append(py)

    return {
        'xx': xx, 'yy': yy,
        'log_intensity': log_intensity.reshape(grid_n, grid_n),
        'intensity': intensity.reshape(grid_n, grid_n),
        'counts': counts.reshape(grid_n, grid_n),
        'points_x': np.array(point_x),
        'points_y': np.array(point_y),
        'n_points': len(point_x)
    }


def plot_lgcp(lgcp_result):
    """Plot LGCP simulation results."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    im0 = axes[0].contourf(lgcp_result['xx'], lgcp_result['yy'],
                            lgcp_result['log_intensity'], levels=20, cmap='RdYlBu_r')
    axes[0].set_title('Log-Intensity Surface S(x)', fontsize=13)
    plt.colorbar(im0, ax=axes[0], shrink=0.8)

    im1 = axes[1].contourf(lgcp_result['xx'], lgcp_result['yy'],
                            lgcp_result['intensity'], levels=20, cmap='YlOrRd')
    axes[1].set_title('Intensity λ(x) = exp(μ + S(x))', fontsize=13)
    plt.colorbar(im1, ax=axes[1], shrink=0.8)

    axes[2].scatter(lgcp_result['points_x'], lgcp_result['points_y'],
                    s=2, alpha=0.5, color='darkred')
    axes[2].set_title(f'Simulated Point Pattern (n={lgcp_result["n_points"]})', fontsize=13)
    axes[2].set_xlim(0, 10); axes[2].set_ylim(0, 10)
    axes[2].set_aspect('equal')

    for ax in axes:
        ax.set_xlabel('Easting'); ax.set_ylabel('Northing')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/lgcp_simulation.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("[LGCP] Figure saved: figures/lgcp_simulation.png")


# ============================================================
# PART 3: Bayesian Spatial Model (INLA/SPDE-like via GP)
# ============================================================

def build_spde_mesh(obs_coords, n_knots=200):
    """Build a mesh approximation (simplified SPDE approach using knots)."""
    from scipy.spatial import Delaunay
    x_range = [obs_coords[:, 0].min() - 0.5, obs_coords[:, 0].max() + 0.5]
    y_range = [obs_coords[:, 1].min() - 0.5, obs_coords[:, 1].max() + 0.5]

    knot_x = np.random.uniform(x_range[0], x_range[1], n_knots)
    knot_y = np.random.uniform(y_range[0], y_range[1], n_knots)
    all_points = np.column_stack([
        np.concatenate([obs_coords[:, 0], knot_x]),
        np.concatenate([obs_coords[:, 1], knot_y])
    ])
    tri = Delaunay(all_points)
    return tri, all_points


def bayesian_spatial_model(obs_coords, cases, population, covariates):
    """
    Bayesian spatial model using Gaussian Process regression as
    a Python-accessible approximation of the INLA/SPDE approach.
    """
    sir = cases / (population / 1000)
    log_sir = np.log(sir + 0.01)

    X = covariates[['temperature', 'rainfall', 'elevation', 'urbanization']].values
    X_std = (X - X.mean(axis=0)) / X.std(axis=0)

    kernel = Matern(length_scale=2.0, nu=1.5) + WhiteKernel(noise_level=0.1)
    gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5, alpha=0.01)

    gp.fit(np.hstack([obs_coords, X_std]), log_sir)

    pred_mean, pred_std = gp.predict(np.hstack([obs_coords, X_std]), return_std=True)

    # Cross-validation
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = []
    for train_idx, test_idx in kf.split(obs_coords):
        gp_cv = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=3, alpha=0.01)
        X_train = np.hstack([obs_coords[train_idx], X_std[train_idx]])
        X_test = np.hstack([obs_coords[test_idx], X_std[test_idx]])
        gp_cv.fit(X_train, log_sir[train_idx])
        pred_cv = gp_cv.predict(X_test)
        cv_scores.append(r2_score(log_sir[test_idx], pred_cv))

    results = {
        'model': gp,
        'pred_mean': pred_mean,
        'pred_std': pred_std,
        'log_sir': log_sir,
        'cv_r2_mean': np.mean(cv_scores),
        'cv_r2_std': np.std(cv_scores),
        'kernel_params': gp.kernel_,
        'log_marginal_likelihood': gp.log_marginal_likelihood_value_
    }

    print(f"[Bayesian GP] CV R²: {results['cv_r2_mean']:.3f} ± {results['cv_r2_std']:.3f}")
    print(f"[Bayesian GP] Optimized kernel: {gp.kernel_}")
    print(f"[Bayesian GP] Log marginal likelihood: {gp.log_marginal_likelihood_value_:.2f}")

    return results


def plot_bayesian_model(obs_coords, bayes_results, xx, yy, true_risk):
    """Plot Bayesian spatial model results."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # True risk
    im0 = axes[0, 0].contourf(xx, yy, true_risk, levels=20, cmap='YlOrRd')
    axes[0, 0].scatter(obs_coords[:, 0], obs_coords[:, 1], c='black', s=10, alpha=0.5)
    axes[0, 0].set_title('True Risk Surface', fontsize=13)
    plt.colorbar(im0, ax=axes[0, 0], shrink=0.8)

    # Predicted mean
    sc1 = axes[0, 1].scatter(obs_coords[:, 0], obs_coords[:, 1],
                              c=bayes_results['pred_mean'], cmap='YlOrRd', s=40, edgecolors='k', linewidth=0.3)
    axes[0, 1].set_title('Predicted log(SIR) — GP Mean', fontsize=13)
    plt.colorbar(sc1, ax=axes[0, 1], shrink=0.8)

    # Posterior uncertainty
    sc2 = axes[1, 0].scatter(obs_coords[:, 0], obs_coords[:, 1],
                              c=bayes_results['pred_std'], cmap='viridis', s=40, edgecolors='k', linewidth=0.3)
    axes[1, 0].set_title('Posterior Uncertainty (Std)', fontsize=13)
    plt.colorbar(sc2, ax=axes[1, 0], shrink=0.8)

    # Observed vs Predicted
    axes[1, 1].scatter(bayes_results['log_sir'], bayes_results['pred_mean'], alpha=0.6, edgecolors='k', linewidth=0.3)
    lims = [min(bayes_results['log_sir'].min(), bayes_results['pred_mean'].min()),
            max(bayes_results['log_sir'].max(), bayes_results['pred_mean'].max())]
    axes[1, 1].plot(lims, lims, 'r--', linewidth=2)
    axes[1, 1].set_xlabel('Observed log(SIR)'); axes[1, 1].set_ylabel('Predicted log(SIR)')
    axes[1, 1].set_title(f'Obs vs Pred (CV R²={bayes_results["cv_r2_mean"]:.3f})', fontsize=13)

    for ax in axes.flat:
        ax.set_xlabel('Easting'); ax.set_ylabel('Northing')
    axes[1, 1].set_xlabel('Observed'); axes[1, 1].set_ylabel('Predicted')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/bayesian_spatial_model.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("[Bayesian GP] Figure saved: figures/bayesian_spatial_model.png")


# ============================================================
# PART 4: Spatial Autocorrelation (Moran's I, Variogram)
# ============================================================

def compute_morans_i(values, coords, k=8):
    """Compute Global Moran's I statistic."""
    from scipy.spatial import cKDTree
    n = len(values)
    tree = cKDTree(coords)
    z = values - values.mean()

    W = np.zeros((n, n))
    for i in range(n):
        dists, indices = tree.query(coords[i], k=k+1)
        for j_idx, j in enumerate(indices[1:]):
            W[i, j] = 1.0 / max(dists[j_idx + 1], 1e-10)

    W_sum = W.sum()
    numerator = n * np.sum(W * np.outer(z, z))
    denominator = W_sum * np.sum(z**2)
    I = numerator / denominator

    # Expected value and variance under null
    E_I = -1.0 / (n - 1)
    S1 = 0.5 * np.sum((W + W.T)**2)
    S2 = np.sum((W.sum(axis=1) + W.sum(axis=0))**2)
    S0 = W_sum
    n2 = n * n

    V_I = (n2 * S1 - n * S2 + 3 * S0**2) / (S0**2 * (n2 - 1)) - E_I**2
    Z_I = (I - E_I) / np.sqrt(max(V_I, 1e-10))
    p_value = 2 * (1 - stats.norm.cdf(abs(Z_I)))

    return {'I': I, 'E_I': E_I, 'V_I': V_I, 'Z': Z_I, 'p_value': p_value}


def compute_empirical_variogram(values, coords, n_bins=15, max_dist=None):
    """Compute empirical semivariogram."""
    n = len(values)
    if max_dist is None:
        max_dist = 0.5 * np.sqrt((coords[:, 0].max() - coords[:, 0].min())**2 +
                                  (coords[:, 1].max() - coords[:, 1].min())**2)

    dists = squareform(pdist(coords))
    bins = np.linspace(0, max_dist, n_bins + 1)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    semivariance = np.zeros(n_bins)
    counts = np.zeros(n_bins)

    for i in range(n):
        for j in range(i + 1, n):
            d = dists[i, j]
            for b in range(n_bins):
                if bins[b] <= d < bins[b + 1]:
                    semivariance[b] += 0.5 * (values[i] - values[j])**2
                    counts[b] += 1
                    break

    mask = counts > 0
    semivariance[mask] /= counts[mask]

    return bin_centers[mask], semivariance[mask], counts[mask]


def fit_variogram_model(h, gamma, model='spherical'):
    """Fit theoretical variogram model to empirical variogram."""
    def spherical_model(h, nugget, sill, range_param):
        result = np.where(
            h <= range_param,
            nugget + (sill - nugget) * (1.5 * h / range_param - 0.5 * (h / range_param)**3),
            sill
        )
        result[h == 0] = 0
        return result

    def exponential_model(h, nugget, sill, range_param):
        result = nugget + (sill - nugget) * (1 - np.exp(-3 * h / range_param))
        result[h == 0] = 0
        return result

    def gaussian_model(h, nugget, sill, range_param):
        result = nugget + (sill - nugget) * (1 - np.exp(-3 * (h / range_param)**2))
        result[h == 0] = 0
        return result

    models = {
        'spherical': spherical_model,
        'exponential': exponential_model,
        'gaussian': gaussian_model
    }

    best_model = None
    best_params = None
    best_sse = np.inf

    for name, func in models.items():
        try:
            p0 = [gamma.min(), gamma.max(), h.max() / 2]
            popt, _ = optimize.curve_fit(func, h, gamma, p0=p0, maxfev=5000,
                                          bounds=([0, 0, 0.01], [np.inf, np.inf, np.inf]))
            pred = func(h, *popt)
            sse = np.sum((gamma - pred)**2)
            if sse < best_sse:
                best_sse = sse
                best_model = name
                best_params = {'nugget': popt[0], 'sill': popt[1], 'range': popt[2]}
        except Exception:
            continue

    return best_model, best_params, best_sse


def plot_spatial_autocorrelation(obs_coords, sir, morans_result, h, gamma, vario_model, vario_params):
    """Plot spatial autocorrelation analysis."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Moran's I scatter
    z = sir - sir.mean()
    from scipy.spatial import cKDTree
    tree = cKDTree(obs_coords)
    spatial_lag = np.zeros(len(sir))
    for i in range(len(sir)):
        _, indices = tree.query(obs_coords[i], k=9)
        neighbors = indices[1:]
        spatial_lag[i] = z[neighbors].mean()

    axes[0].scatter(z, spatial_lag, alpha=0.6, edgecolors='k', linewidth=0.3)
    slope, intercept = np.polyfit(z, spatial_lag, 1)
    x_line = np.linspace(z.min(), z.max(), 100)
    axes[0].plot(x_line, slope * x_line + intercept, 'r-', linewidth=2,
                 label=f"Moran's I = {morans_result['I']:.3f}")
    axes[0].axhline(0, color='gray', linestyle='--', alpha=0.5)
    axes[0].axvline(0, color='gray', linestyle='--', alpha=0.5)
    axes[0].set_xlabel('Standardized SIR')
    axes[0].set_ylabel('Spatial Lag of SIR')
    axes[0].set_title(f"Moran's I Scatter Plot\n(I={morans_result['I']:.3f}, p={morans_result['p_value']:.4f})")
    axes[0].legend()

    # Empirical variogram with fitted model
    axes[1].scatter(h, gamma, s=60, c='steelblue', edgecolors='k', zorder=5, label='Empirical')
    h_fine = np.linspace(0, h.max(), 200)
    if vario_model == 'spherical':
        def model_func(h, n, s, r):
            result = np.where(h <= r, n + (s - n) * (1.5 * h / r - 0.5 * (h / r)**3), s)
            result[h == 0] = 0
            return result
    elif vario_model == 'exponential':
        def model_func(h, n, s, r):
            result = n + (s - n) * (1 - np.exp(-3 * h / r))
            result[h == 0] = 0
            return result
    else:
        def model_func(h, n, s, r):
            result = n + (s - n) * (1 - np.exp(-3 * (h / r)**2))
            result[h == 0] = 0
            return result

    gamma_fit = model_func(h_fine, vario_params['nugget'], vario_params['sill'], vario_params['range'])
    axes[1].plot(h_fine, gamma_fit, 'r-', linewidth=2,
                 label=f'{vario_model.capitalize()} model')
    axes[1].axhline(vario_params['sill'], color='gray', linestyle=':', alpha=0.5, label=f"Sill={vario_params['sill']:.2f}")
    axes[1].axvline(vario_params['range'], color='gray', linestyle='--', alpha=0.5, label=f"Range={vario_params['range']:.2f}")
    axes[1].set_xlabel('Distance (h)')
    axes[1].set_ylabel('Semivariance γ(h)')
    axes[1].set_title('Empirical & Fitted Variogram')
    axes[1].legend(fontsize=9)

    # Spatial distribution of SIR
    sc = axes[2].scatter(obs_coords[:, 0], obs_coords[:, 1], c=sir, cmap='YlOrRd',
                          s=50, edgecolors='k', linewidth=0.3)
    axes[2].set_title('Spatial Distribution of SIR')
    axes[2].set_xlabel('Easting'); axes[2].set_ylabel('Northing')
    plt.colorbar(sc, ax=axes[2], shrink=0.8, label='SIR')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/spatial_autocorrelation.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("[Spatial Autocorrelation] Figure saved: figures/spatial_autocorrelation.png")


# ============================================================
# PART 5: Ecological Study Confounding Bias
# ============================================================

def ecological_bias_analysis(obs_coords, cases, population, covariates):
    """Analyze and mitigate ecological confounding bias."""
    sir = cases / (population / 1000)
    log_sir = np.log(sir + 0.01)

    X = covariates[['temperature', 'rainfall', 'elevation', 'urbanization']].values
    X_with_const = sm.add_constant(X)
    ols_naive = sm.OLS(log_sir, X_with_const).fit()

    # Add spatial covariates to mitigate confounding
    X_spatial = np.column_stack([X, obs_coords[:, 0], obs_coords[:, 1],
                                  obs_coords[:, 0]**2, obs_coords[:, 1]**2,
                                  obs_coords[:, 0] * obs_coords[:, 1]])
    X_spatial_const = sm.add_constant(X_spatial)
    ols_spatial = sm.OLS(log_sir, X_spatial_const).fit()

    # Propensity score stratification for urbanization
    median_urban = covariates['urbanization'].median()
    high_urban = covariates['urbanization'] >= median_urban
    low_urban = ~high_urban

    results = {
        'naive_model': ols_naive,
        'spatial_model': ols_spatial,
        'naive_aic': ols_naive.aic,
        'spatial_aic': ols_spatial.aic,
        'naive_r2': ols_naive.rsquared,
        'spatial_r2': ols_spatial.rsquared,
        'coef_change': {},
    }

    var_names = ['temperature', 'rainfall', 'elevation', 'urbanization']
    for i, name in enumerate(var_names):
        naive_coef = ols_naive.params[i + 1]
        spatial_coef = ols_spatial.params[i + 1]
        change_pct = 100 * (spatial_coef - naive_coef) / abs(naive_coef) if naive_coef != 0 else 0
        results['coef_change'][name] = {
            'naive': naive_coef,
            'spatial': spatial_coef,
            'change_pct': change_pct
        }

    print("[Ecological Bias] Naive OLS R²:", f"{results['naive_r2']:.3f}, AIC: {results['naive_aic']:.1f}")
    print("[Ecological Bias] Spatial OLS R²:", f"{results['spatial_r2']:.3f}, AIC: {results['spatial_aic']:.1f}")
    for name, vals in results['coef_change'].items():
        print(f"  {name}: naive={vals['naive']:.4f}, spatial={vals['spatial']:.4f}, change={vals['change_pct']:.1f}%")

    return results


def plot_ecological_bias(eco_results):
    """Plot ecological bias analysis results."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    var_names = list(eco_results['coef_change'].keys())
    naive_coefs = [eco_results['coef_change'][v]['naive'] for v in var_names]
    spatial_coefs = [eco_results['coef_change'][v]['spatial'] for v in var_names]

    x_pos = np.arange(len(var_names))
    width = 0.35
    axes[0].bar(x_pos - width/2, naive_coefs, width, label='Naive OLS', color='salmon', edgecolor='k')
    axes[0].bar(x_pos + width/2, spatial_coefs, width, label='Spatial OLS', color='steelblue', edgecolor='k')
    axes[0].set_xticks(x_pos)
    axes[0].set_xticklabels(var_names, rotation=30, ha='right')
    axes[0].set_ylabel('Coefficient')
    axes[0].set_title('Covariate Effects: Naive vs Spatial Adjustment')
    axes[0].legend()
    axes[0].axhline(0, color='k', linestyle='-', linewidth=0.5)

    # AIC comparison
    models = ['Naive OLS', 'Spatial OLS']
    aics = [eco_results['naive_aic'], eco_results['spatial_aic']]
    colors = ['salmon', 'steelblue']
    axes[1].bar(models, aics, color=colors, edgecolor='k')
    axes[1].set_ylabel('AIC')
    axes[1].set_title('Model Comparison (AIC)')

    # R² comparison
    r2s = [eco_results['naive_r2'], eco_results['spatial_r2']]
    axes[2].bar(models, r2s, color=colors, edgecolor='k')
    axes[2].set_ylabel('R²')
    axes[2].set_title('Model Comparison (R²)')
    axes[2].set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/ecological_bias_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("[Ecological Bias] Figure saved: figures/ecological_bias_analysis.png")


# ============================================================
# PART 6: Spatiotemporal Knot-Based Spline Model
# ============================================================

def generate_spatiotemporal_data(obs_coords, n_times=12):
    """Generate spatiotemporal disease data (monthly over 1 year)."""
    n_locations = len(obs_coords)
    times = np.arange(n_times)
    records = []

    for t in times:
        seasonal = 0.5 * np.sin(2 * np.pi * t / 12 - np.pi / 3)
        for i in range(n_locations):
            x, y = obs_coords[i]
            spatial_effect = (
                0.3 * np.exp(-((x - 3)**2 + (y - 7)**2) / 4) +
                0.5 * np.exp(-((x - 7)**2 + (y - 3)**2) / 3)
            )
            trend = 0.02 * t
            log_rate = -1.5 + 2.0 * spatial_effect + seasonal + trend + np.random.normal(0, 0.3)
            rate = np.exp(log_rate)
            pop = np.random.randint(5000, 30000)
            cases = np.random.poisson(rate * pop / 1000)
            records.append({
                'x': x, 'y': y, 'time': t,
                'cases': cases, 'population': pop,
                'rate': cases / (pop / 1000)
            })

    return pd.DataFrame(records)


def knot_based_spline_model(st_data, n_spatial_knots=15, n_temporal_knots=6):
    """Fit a spatiotemporal knot-based spline model."""
    x = st_data['x'].values
    y = st_data['y'].values
    t = st_data['time'].values
    log_rate = np.log(st_data['rate'].values + 0.01)

    # Spatial knots
    x_knots = np.linspace(x.min(), x.max(), int(np.sqrt(n_spatial_knots)))
    y_knots = np.linspace(y.min(), y.max(), int(np.sqrt(n_spatial_knots)))
    xx_k, yy_k = np.meshgrid(x_knots, y_knots)
    spatial_knots = np.column_stack([xx_k.ravel(), yy_k.ravel()])

    # Temporal knots
    t_knots = np.linspace(t.min(), t.max(), n_temporal_knots)

    # Radial basis functions for spatial component
    def rbf(d, epsilon=1.0):
        return np.exp(-(epsilon * d)**2)

    n = len(x)
    n_sk = len(spatial_knots)
    n_tk = len(t_knots)

    # Spatial basis
    spatial_basis = np.zeros((n, n_sk))
    for j in range(n_sk):
        d = np.sqrt((x - spatial_knots[j, 0])**2 + (y - spatial_knots[j, 1])**2)
        spatial_basis[:, j] = rbf(d, epsilon=0.5)

    # Temporal basis (B-spline-like using Gaussian RBF)
    temporal_basis = np.zeros((n, n_tk))
    for j in range(n_tk):
        temporal_basis[:, j] = rbf(t - t_knots[j], epsilon=0.3)

    # Interaction terms (subset)
    n_interactions = min(n_sk * n_tk, 50)
    interaction_basis = np.zeros((n, n_interactions))
    idx = 0
    for j in range(n_sk):
        for k in range(n_tk):
            if idx >= n_interactions:
                break
            interaction_basis[:, idx] = spatial_basis[:, j] * temporal_basis[:, k]
            idx += 1
        if idx >= n_interactions:
            break

    # Full design matrix
    design = np.column_stack([np.ones(n), spatial_basis, temporal_basis, interaction_basis])

    # Ridge regression for regularization
    from sklearn.linear_model import Ridge
    ridge = Ridge(alpha=1.0)
    ridge.fit(design, log_rate)
    pred = ridge.predict(design)

    # Train/test split by time
    train_mask = t <= 9
    test_mask = t > 9
    ridge_train = Ridge(alpha=1.0)
    ridge_train.fit(design[train_mask], log_rate[train_mask])
    pred_test = ridge_train.predict(design[test_mask])

    rmse_train = np.sqrt(mean_squared_error(log_rate[train_mask], ridge.predict(design[train_mask])))
    rmse_test = np.sqrt(mean_squared_error(log_rate[test_mask], pred_test))
    r2_train = r2_score(log_rate[train_mask], ridge.predict(design[train_mask]))
    r2_test = r2_score(log_rate[test_mask], pred_test)
    mae_test = mean_absolute_error(log_rate[test_mask], pred_test)

    results = {
        'model': ridge,
        'pred': pred,
        'pred_test': pred_test,
        'log_rate': log_rate,
        'rmse_train': rmse_train,
        'rmse_test': rmse_test,
        'r2_train': r2_train,
        'r2_test': r2_test,
        'mae_test': mae_test,
        'spatial_knots': spatial_knots,
        'temporal_knots': t_knots,
        'train_mask': train_mask,
        'test_mask': test_mask,
    }

    print(f"[Spline Model] Train RMSE: {rmse_train:.3f}, R²: {r2_train:.3f}")
    print(f"[Spline Model] Test  RMSE: {rmse_test:.3f}, R²: {r2_test:.3f}, MAE: {mae_test:.3f}")

    return results


def plot_spatiotemporal_model(st_data, spline_results):
    """Plot spatiotemporal model results."""
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(2, 3, figure=fig)

    # Time series aggregated
    ax0 = fig.add_subplot(gs[0, 0])
    monthly_obs = st_data.groupby('time')['rate'].mean()
    monthly_pred = pd.Series(np.exp(spline_results['pred']), index=st_data.index).groupby(st_data['time']).mean()
    ax0.plot(monthly_obs.index, monthly_obs.values, 'bo-', label='Observed', linewidth=2)
    ax0.plot(monthly_pred.index, monthly_pred.values, 'r^--', label='Predicted', linewidth=2)
    ax0.axvline(9.5, color='gray', linestyle=':', label='Train/Test split')
    ax0.set_xlabel('Month'); ax0.set_ylabel('Mean Rate')
    ax0.set_title('Temporal Trend: Observed vs Predicted')
    ax0.legend()

    # Spatial knots
    ax1 = fig.add_subplot(gs[0, 1])
    ax1.scatter(st_data['x'], st_data['y'], c=st_data['rate'], cmap='YlOrRd', s=5, alpha=0.3)
    ax1.scatter(spline_results['spatial_knots'][:, 0], spline_results['spatial_knots'][:, 1],
                c='blue', marker='x', s=100, linewidths=2, label='Spatial knots')
    ax1.set_title('Spatial Knot Placement')
    ax1.set_xlabel('Easting'); ax1.set_ylabel('Northing')
    ax1.legend()

    # Observed vs predicted
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.scatter(spline_results['log_rate'], spline_results['pred'], alpha=0.3, s=10, c='steelblue')
    lims = [spline_results['log_rate'].min(), spline_results['log_rate'].max()]
    ax2.plot(lims, lims, 'r--', linewidth=2)
    ax2.set_xlabel('Observed log(rate)'); ax2.set_ylabel('Predicted log(rate)')
    ax2.set_title(f"Obs vs Pred (R²={spline_results['r2_train']:.3f})")

    # Risk maps for different months
    for i, month in enumerate([0, 5, 11]):
        ax = fig.add_subplot(gs[1, i])
        mask = st_data['time'] == month
        sub = st_data[mask]
        sc = ax.scatter(sub['x'], sub['y'],
                        c=np.exp(spline_results['pred'][mask.values]),
                        cmap='YlOrRd', s=50, edgecolors='k', linewidth=0.3)
        ax.set_title(f'Predicted Risk — Month {month + 1}')
        ax.set_xlabel('Easting'); ax.set_ylabel('Northing')
        plt.colorbar(sc, ax=ax, shrink=0.8)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/spatiotemporal_spline.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("[Spline Model] Figure saved: figures/spatiotemporal_spline.png")


# ============================================================
# PART 7: Malaria/Dengue Risk Mapping Case Study
# ============================================================

def disease_case_study():
    """Simulated malaria and dengue risk mapping case study."""
    np.random.seed(123)
    n_districts = 80
    coords = np.random.uniform(0, 10, (n_districts, 2))

    # Malaria risk — higher in rural, warm, wet areas
    malaria_risk = (
        0.4 * np.exp(-((coords[:, 0] - 2)**2 + (coords[:, 1] - 8)**2) / 5) +
        0.3 * np.exp(-((coords[:, 0] - 8)**2 + (coords[:, 1] - 2)**2) / 4) +
        np.random.normal(0, 0.02, n_districts)
    )
    malaria_risk = np.clip(malaria_risk, 0, 1)

    # Dengue risk — higher in urban, moderate-temperature areas
    dengue_risk = (
        0.5 * np.exp(-((coords[:, 0] - 5)**2 + (coords[:, 1] - 5)**2) / 3) +
        0.2 * np.exp(-((coords[:, 0] - 7)**2 + (coords[:, 1] - 7)**2) / 4) +
        np.random.normal(0, 0.02, n_districts)
    )
    dengue_risk = np.clip(dengue_risk, 0, 1)

    pop = np.random.randint(10000, 100000, n_districts)
    malaria_cases = np.random.poisson(malaria_risk * pop / 500)
    dengue_cases = np.random.poisson(dengue_risk * pop / 600)

    malaria_sir = malaria_cases / (pop / 1000)
    dengue_sir = dengue_cases / (pop / 1000)

    # GP risk surface prediction
    from sklearn.gaussian_process.kernels import Matern, WhiteKernel

    grid_x = np.linspace(0, 10, 50)
    grid_y = np.linspace(0, 10, 50)
    gx, gy = np.meshgrid(grid_x, grid_y)
    grid_coords = np.column_stack([gx.ravel(), gy.ravel()])

    results = {}
    for disease, sir_vals in [('Malaria', malaria_sir), ('Dengue', dengue_sir)]:
        kernel = Matern(length_scale=2.0, nu=1.5) + WhiteKernel(noise_level=0.1)
        gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=3, alpha=0.01)
        gp.fit(coords, np.log(sir_vals + 0.01))
        pred_mean, pred_std = gp.predict(grid_coords, return_std=True)

        # Exceedance probability P(risk > threshold)
        threshold = np.log(np.median(sir_vals) * 1.5)
        exceed_prob = 1 - stats.norm.cdf(threshold, loc=pred_mean, scale=pred_std)

        # Moran's I
        morans = compute_morans_i(sir_vals, coords)

        results[disease] = {
            'sir': sir_vals, 'cases': malaria_cases if disease == 'Malaria' else dengue_cases,
            'pred_mean': pred_mean.reshape(50, 50),
            'pred_std': pred_std.reshape(50, 50),
            'exceed_prob': exceed_prob.reshape(50, 50),
            'morans_I': morans['I'],
            'morans_p': morans['p_value'],
            'gp_kernel': str(gp.kernel_)
        }
        print(f"[{disease}] Moran's I: {morans['I']:.3f}, p-value: {morans['p_value']:.4f}")

    return coords, pop, gx, gy, results


def plot_disease_case_study(coords, pop, gx, gy, results):
    """Plot comprehensive disease risk mapping."""
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))

    for row, disease in enumerate(['Malaria', 'Dengue']):
        r = results[disease]

        # Observed SIR
        sc0 = axes[row, 0].scatter(coords[:, 0], coords[:, 1], c=r['sir'],
                                    cmap='YlOrRd', s=60, edgecolors='k', linewidth=0.3)
        axes[row, 0].set_title(f'{disease} — Observed SIR')
        plt.colorbar(sc0, ax=axes[row, 0], shrink=0.7)

        # Predicted risk surface
        im1 = axes[row, 1].contourf(gx, gy, r['pred_mean'], levels=20, cmap='YlOrRd')
        axes[row, 1].set_title(f'{disease} — Predicted log(SIR)')
        plt.colorbar(im1, ax=axes[row, 1], shrink=0.7)

        # Uncertainty
        im2 = axes[row, 2].contourf(gx, gy, r['pred_std'], levels=20, cmap='viridis')
        axes[row, 2].set_title(f'{disease} — Uncertainty (Std)')
        plt.colorbar(im2, ax=axes[row, 2], shrink=0.7)

        # Exceedance probability
        im3 = axes[row, 3].contourf(gx, gy, r['exceed_prob'], levels=20, cmap='RdYlGn_r')
        axes[row, 3].set_title(f'{disease} — P(Risk > threshold)')
        plt.colorbar(im3, ax=axes[row, 3], shrink=0.7)

        for ax in axes[row]:
            ax.set_xlabel('Easting'); ax.set_ylabel('Northing')

    plt.suptitle('Disease Risk Mapping Case Study: Malaria & Dengue', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/disease_risk_mapping.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("[Case Study] Figure saved: figures/disease_risk_mapping.png")


# ============================================================
# PART 8: SPDE Mesh Visualization
# ============================================================

def plot_spde_mesh(obs_coords):
    """Visualize SPDE triangulation mesh."""
    tri, all_points = build_spde_mesh(obs_coords)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    axes[0].triplot(all_points[:, 0], all_points[:, 1], tri.simplices, 'b-', alpha=0.3, linewidth=0.5)
    axes[0].scatter(obs_coords[:, 0], obs_coords[:, 1], c='red', s=30, zorder=5, label='Observation points')
    axes[0].scatter(all_points[len(obs_coords):, 0], all_points[len(obs_coords):, 1],
                    c='blue', s=5, alpha=0.5, label='Mesh vertices')
    axes[0].set_title('SPDE Triangulation Mesh')
    axes[0].set_xlabel('Easting'); axes[0].set_ylabel('Northing')
    axes[0].legend()

    # Basis function example
    from scipy.interpolate import LinearNDInterpolator
    center_idx = len(obs_coords)  # first knot
    basis = np.zeros(len(all_points))
    basis[center_idx] = 1.0
    interp = LinearNDInterpolator(all_points, basis)
    gx = np.linspace(all_points[:, 0].min(), all_points[:, 0].max(), 100)
    gy = np.linspace(all_points[:, 1].min(), all_points[:, 1].max(), 100)
    GX, GY = np.meshgrid(gx, gy)
    GZ = interp(GX, GY)
    GZ = np.nan_to_num(GZ, 0)
    axes[1].contourf(GX, GY, GZ, levels=20, cmap='Blues')
    axes[1].set_title('Example SPDE Basis Function (Single Vertex)')
    axes[1].set_xlabel('Easting'); axes[1].set_ylabel('Northing')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/spde_mesh.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("[SPDE Mesh] Figure saved: figures/spde_mesh.png")


# ============================================================
# PART 9: Model Comparison Summary
# ============================================================

def model_comparison_plot(bayes_results, spline_results, eco_results):
    """Create comprehensive model comparison visualization."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # R² comparison
    models = ['GP (Bayesian)', 'Spline (Train)', 'Spline (Test)', 'Naive OLS', 'Spatial OLS']
    r2_vals = [
        bayes_results['cv_r2_mean'],
        spline_results['r2_train'],
        spline_results['r2_test'],
        eco_results['naive_r2'],
        eco_results['spatial_r2']
    ]
    colors = ['steelblue', 'forestgreen', 'darkgreen', 'salmon', 'orange']
    axes[0].barh(models, r2_vals, color=colors, edgecolor='k')
    axes[0].set_xlabel('R²')
    axes[0].set_title('Model Performance Comparison (R²)')
    axes[0].set_xlim(0, 1)
    for i, v in enumerate(r2_vals):
        axes[0].text(v + 0.01, i, f'{v:.3f}', va='center')

    # RMSE comparison for spline
    metrics = ['RMSE (Train)', 'RMSE (Test)', 'MAE (Test)']
    vals = [spline_results['rmse_train'], spline_results['rmse_test'], spline_results['mae_test']]
    axes[1].bar(metrics, vals, color=['steelblue', 'salmon', 'orange'], edgecolor='k')
    axes[1].set_ylabel('Error')
    axes[1].set_title('Spatiotemporal Spline Model Errors')
    for i, v in enumerate(vals):
        axes[1].text(i, v + 0.01, f'{v:.3f}', ha='center')

    # AIC comparison
    aic_models = ['Naive OLS', 'Spatial OLS']
    aic_vals = [eco_results['naive_aic'], eco_results['spatial_aic']]
    axes[2].bar(aic_models, aic_vals, color=['salmon', 'steelblue'], edgecolor='k')
    axes[2].set_ylabel('AIC')
    axes[2].set_title('Model Selection (AIC — lower is better)')

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/model_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("[Comparison] Figure saved: figures/model_comparison.png")


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    print("=" * 70)
    print("GEOSTATISTICAL FRAMEWORK FOR DISEASE RISK ANALYSIS")
    print("=" * 70)

    # --- Data Generation ---
    print("\n[1] Generating synthetic spatial data...")
    xx, yy, coords_grid, obs_coords = generate_spatial_domain(n_regions=120)
    covariates = generate_covariates(obs_coords)
    true_risk = generate_true_risk_surface(xx, yy)
    cases, population, risk_at_obs = generate_disease_cases(obs_coords, covariates, true_risk, xx, yy)
    sir = cases / (population / 1000)
    print(f"  Generated {len(obs_coords)} observation locations, {cases.sum()} total cases")

    # --- LGCP ---
    print("\n[2] Log-Gaussian Cox Process simulation...")
    lgcp_result = simulate_lgcp(grid_n=30, sigma2=1.5, kappa=1.2, mu=-1.5)
    plot_lgcp(lgcp_result)

    # --- Bayesian Spatial Model ---
    print("\n[3] Bayesian spatial model (GP-SPDE approximation)...")
    bayes_results = bayesian_spatial_model(obs_coords, cases, population, covariates)
    plot_bayesian_model(obs_coords, bayes_results, xx, yy, true_risk)

    # --- SPDE Mesh ---
    print("\n[4] SPDE mesh construction...")
    plot_spde_mesh(obs_coords)

    # --- Spatial Autocorrelation ---
    print("\n[5] Spatial autocorrelation analysis...")
    morans_result = compute_morans_i(sir, obs_coords)
    h, gamma, counts = compute_empirical_variogram(sir, obs_coords)
    vario_model, vario_params, vario_sse = fit_variogram_model(h, gamma)
    print(f"  Moran's I: {morans_result['I']:.4f}, Z: {morans_result['Z']:.3f}, p: {morans_result['p_value']:.4f}")
    print(f"  Best variogram: {vario_model} (nugget={vario_params['nugget']:.3f}, "
          f"sill={vario_params['sill']:.3f}, range={vario_params['range']:.3f})")
    plot_spatial_autocorrelation(obs_coords, sir, morans_result, h, gamma, vario_model, vario_params)

    # --- Ecological Bias ---
    print("\n[6] Ecological confounding bias analysis...")
    eco_results = ecological_bias_analysis(obs_coords, cases, population, covariates)
    plot_ecological_bias(eco_results)

    # --- Spatiotemporal Spline ---
    print("\n[7] Spatiotemporal knot-based spline model...")
    st_data = generate_spatiotemporal_data(obs_coords, n_times=12)
    spline_results = knot_based_spline_model(st_data)
    plot_spatiotemporal_model(st_data, spline_results)

    # --- Disease Case Study ---
    print("\n[8] Malaria/Dengue risk mapping case study...")
    cs_coords, cs_pop, gx, gy, cs_results = disease_case_study()
    plot_disease_case_study(cs_coords, cs_pop, gx, gy, cs_results)

    # --- Model Comparison ---
    print("\n[9] Model comparison summary...")
    model_comparison_plot(bayes_results, spline_results, eco_results)

    # --- Summary Statistics ---
    print("\n" + "=" * 70)
    print("EXPERIMENT SUMMARY")
    print("=" * 70)
    summary = {
        'LGCP_n_points': lgcp_result['n_points'],
        'Bayesian_GP_CV_R2': f"{bayes_results['cv_r2_mean']:.3f} ± {bayes_results['cv_r2_std']:.3f}",
        'Morans_I': f"{morans_result['I']:.4f} (p={morans_result['p_value']:.4f})",
        'Variogram_model': f"{vario_model} (range={vario_params['range']:.2f})",
        'Naive_OLS_R2': f"{eco_results['naive_r2']:.3f}",
        'Spatial_OLS_R2': f"{eco_results['spatial_r2']:.3f}",
        'Spline_Train_R2': f"{spline_results['r2_train']:.3f}",
        'Spline_Test_R2': f"{spline_results['r2_test']:.3f}",
        'Spline_Test_RMSE': f"{spline_results['rmse_test']:.3f}",
        'Malaria_Morans_I': f"{cs_results['Malaria']['morans_I']:.3f} (p={cs_results['Malaria']['morans_p']:.4f})",
        'Dengue_Morans_I': f"{cs_results['Dengue']['morans_I']:.3f} (p={cs_results['Dengue']['morans_p']:.4f})",
    }
    for k, v in summary.items():
        print(f"  {k}: {v}")

    print("\nAll figures saved to figures/ directory.")
    print("=" * 70)

    return summary


if __name__ == '__main__':
    main()
