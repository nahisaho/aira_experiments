import json
import math
import random
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.linalg import solve
from scipy.special import gammaln
from scipy.stats import norm

SEED = 20250314


def ensure_directories(base_dir: Path) -> dict:
    paths = {
        'base': base_dir,
        'src': base_dir / 'src',
        'figures': base_dir / 'figures',
        'results': base_dir / 'results',
        'data': base_dir / 'data',
        'logs': base_dir / 'logs',
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def build_grid(n_rows: int = 14, n_cols: int = 14):
    area_ids = np.arange(n_rows * n_cols)
    rows, cols = np.divmod(area_ids, n_cols)
    neighbors = [[] for _ in area_ids]
    adjacency = np.zeros((len(area_ids), len(area_ids)), dtype=float)
    for idx, r, c in zip(area_ids, rows, cols):
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < n_rows and 0 <= nc < n_cols:
                j = nr * n_cols + nc
                adjacency[idx, j] = 1.0
                neighbors[idx].append(j)
    degree = np.diag(adjacency.sum(axis=1))
    precision = degree - adjacency
    coords = pd.DataFrame({'area_id': area_ids, 'row': rows, 'col': cols})
    return coords, adjacency, precision, neighbors


def simulate_spatial_field(precision: np.ndarray, rng: np.random.Generator, scale: float = 1.0, jitter: float = 0.15):
    n = precision.shape[0]
    cov = np.linalg.inv(precision + np.eye(n) * jitter)
    field = rng.multivariate_normal(mean=np.zeros(n), cov=cov)
    field = (field - field.mean()) / (field.std() + 1e-8)
    return scale * field


def logistic(x):
    return 1.0 / (1.0 + np.exp(-x))


def simulate_disease_data(coords: pd.DataFrame, precision: np.ndarray, n_months: int = 24):
    rng = np.random.default_rng(SEED)
    random.seed(SEED)

    n = len(coords)
    row_norm = (coords['row'].to_numpy() - coords['row'].mean()) / coords['row'].std()
    col_norm = (coords['col'].to_numpy() - coords['col'].mean()) / coords['col'].std()

    urban_hub_1 = np.exp(-(((row_norm + 0.8) ** 2) + ((col_norm - 0.6) ** 2)) / 0.55)
    urban_hub_2 = np.exp(-(((row_norm - 0.5) ** 2) + ((col_norm + 0.9) ** 2)) / 0.65)
    coastal_band = np.exp(-((row_norm + 1.1) ** 2) / 0.75)
    terrain_gradient = 1.4 * row_norm - 0.9 * col_norm

    elevation = 220 + 280 * (terrain_gradient - terrain_gradient.min()) / (terrain_gradient.max() - terrain_gradient.min())
    elevation += rng.normal(0, 25, size=n)
    elevation = np.clip(elevation, 30, None)

    pop_density = 120 + 740 * urban_hub_1 + 590 * urban_hub_2 + 240 * coastal_band
    pop_density += 50 * rng.lognormal(mean=0.0, sigma=0.45, size=n)
    pop_density = np.clip(pop_density, 60, None)

    urbanization = logistic(-1.1 + 0.0036 * pop_density - 0.0015 * elevation + rng.normal(0, 0.25, size=n))
    population = (2800 + 11 * pop_density + 1400 * urbanization + rng.normal(0, 260, size=n)).round().astype(int)
    population = np.clip(population, 1800, None)

    shared_field = simulate_spatial_field(precision, rng, scale=0.50)
    malaria_field = 0.7 * shared_field + simulate_spatial_field(precision, rng, scale=0.45)
    dengue_field = 0.6 * shared_field + simulate_spatial_field(precision, rng, scale=0.55)

    rows = []
    month_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    for month_index in range(n_months):
        month = month_index % 12 + 1
        year = 2023 + month_index // 12
        theta = 2 * np.pi * month_index / 12.0
        temp = 27.2 + 1.9 * np.sin(theta - 0.55) - 0.0065 * elevation + 0.75 * coastal_band
        temp += rng.normal(0, 0.45, size=n)
        precip = 155 + 85 * np.cos(theta - 1.2) + 40 * coastal_band + 25 * (row_norm + 0.6)
        precip += 10 * rng.normal(0, 1.0, size=n)
        precip = np.clip(precip, 30, None)

        temp_z = (temp - temp.mean()) / temp.std()
        precip_z = (precip - precip.mean()) / precip.std()
        elev_z = (elevation - elevation.mean()) / elevation.std()
        log_density = np.log(pop_density)
        density_z = (log_density - log_density.mean()) / log_density.std()
        urban_z = (urbanization - urbanization.mean()) / urbanization.std()

        malaria_lp = (
            -7.70
            + 0.24 * temp_z
            + 0.29 * precip_z
            - 0.17 * elev_z
            - 0.08 * density_z
            - 0.15 * urban_z
            + 0.28 * np.sin(theta - 0.3)
            + malaria_field
        )
        dengue_lp = (
            -8.00
            + 0.31 * temp_z
            + 0.12 * precip_z
            - 0.06 * elev_z
            + 0.14 * density_z
            + 0.24 * urban_z
            + 0.36 * np.cos(theta - 0.7)
            + dengue_field
        )

        malaria_mu = population * np.exp(malaria_lp)
        dengue_mu = population * np.exp(dengue_lp)
        malaria_cases = rng.poisson(np.clip(malaria_mu, 1e-5, None))
        dengue_cases = rng.poisson(np.clip(dengue_mu, 1e-5, None))

        frame = pd.DataFrame({
            'area_id': coords['area_id'],
            'row': coords['row'],
            'col': coords['col'],
            'year': year,
            'month': month,
            'month_label': f'{year}-{month:02d}',
            'month_name': month_names[month - 1],
            'temperature': np.round(temp, 3),
            'precipitation': np.round(precip, 3),
            'elevation': np.round(elevation, 3),
            'population_density': np.round(pop_density, 3),
            'urbanization_index': np.round(urbanization, 4),
            'population': population,
            'malaria_cases': malaria_cases,
            'dengue_cases': dengue_cases,
        })
        rows.append(frame)

    df = pd.concat(rows, ignore_index=True)
    return df


def aggregate_area_data(df: pd.DataFrame):
    area_df = (
        df.groupby(['area_id', 'row', 'col'], as_index=False)
        .agg({
            'population': 'first',
            'temperature': 'mean',
            'precipitation': 'mean',
            'elevation': 'first',
            'population_density': 'first',
            'urbanization_index': 'first',
            'malaria_cases': 'sum',
            'dengue_cases': 'sum',
        })
        .rename(columns={'temperature': 'temperature_mean', 'precipitation': 'precipitation_mean'})
    )
    for disease in ['malaria', 'dengue']:
        cases_col = f'{disease}_cases'
        total_cases = area_df[cases_col].sum()
        expected = total_cases * area_df['population'] / area_df['population'].sum()
        area_df[f'{disease}_expected'] = expected
        area_df[f'{disease}_smr'] = area_df[cases_col] / np.clip(expected, 1e-8, None)
        area_df[f'{disease}_rate_per_100k'] = area_df[cases_col] / area_df['population'] * 100000 / 2.0
    return area_df


def make_design_matrix(area_df: pd.DataFrame):
    covariate_columns = [
        'temperature_mean',
        'precipitation_mean',
        'elevation',
        'population_density',
        'urbanization_index',
    ]
    raw = area_df[covariate_columns].copy()
    raw['population_density'] = np.log(raw['population_density'])
    means = raw.mean()
    stds = raw.std(ddof=0).replace(0, 1.0)
    standardized = (raw - means) / stds
    standardized.insert(0, 'Intercept', 1.0)
    return standardized.to_numpy(), ['Intercept'] + covariate_columns, means.to_dict(), stds.to_dict()


def negative_log_posterior(
    y,
    E,
    X,
    u,
    v,
    beta,
    tau_s,
    tau_u,
    Q,
    prior_var=100.0,
    ridge=1e-3,
    sum_to_zero=5.0,
):
    eta = X @ beta + u + v
    mu = np.clip(E * np.exp(eta), 1e-10, None)
    log_like = np.sum(y * np.log(mu) - mu - gammaln(y + 1.0))
    rank_q = Q.shape[0] - 1
    penalty = (
        0.5 * tau_s * (u @ (Q @ u))
        + 0.5 * tau_u * np.dot(v, v)
        + 0.5 * ridge * np.dot(u, u)
        + 0.5 * sum_to_zero * (u.sum() ** 2)
    )
    beta_penalty = 0.5 * np.dot(beta, beta) / prior_var
    hyper_prior = 0.5 * (tau_s + tau_u) - 0.5 * (rank_q * np.log(tau_s) + len(v) * np.log(tau_u))
    return -(log_like) + penalty + beta_penalty + hyper_prior


def build_hessian_and_gradient(
    y,
    E,
    X,
    beta,
    u,
    v,
    tau_s,
    tau_u,
    Q,
    prior_var=100.0,
    ridge=1e-3,
    sum_to_zero=5.0,
):
    eta = X @ beta + u + v
    mu = np.clip(E * np.exp(eta), 1e-10, None)
    resid = mu - y
    w = mu
    n = len(y)
    p = X.shape[1]
    prior_precision = np.eye(p) / prior_var

    g_beta = X.T @ resid + prior_precision @ beta
    g_u = resid + tau_s * (Q @ u) + ridge * u + sum_to_zero * u.sum()
    g_v = resid + tau_u * v

    xw = X * w[:, None]
    h_beta = X.T @ xw + prior_precision
    cross = X.T * w
    diag_w = np.diag(w)
    h_u = diag_w + tau_s * Q + ridge * np.eye(n) + sum_to_zero * np.ones((n, n))
    h_v = diag_w + tau_u * np.eye(n)
    h_uv = diag_w

    hessian = np.block([
        [h_beta, cross, cross],
        [cross.T, h_u, h_uv],
        [cross.T, h_uv, h_v],
    ])
    gradient = np.concatenate([g_beta, g_u, g_v])
    return gradient, hessian, eta, mu


def fit_bym_laplace(y, E, X, Q, disease_name, tau_s_grid=None, tau_u_grid=None, sum_to_zero=5.0):
    if tau_s_grid is None:
        tau_s_grid = [0.4, 0.9, 1.8, 3.6]
    if tau_u_grid is None:
        tau_u_grid = [0.4, 0.9, 1.8, 3.6]

    n = len(y)
    p = X.shape[1]
    best = None
    beta_start = np.zeros(p)
    beta_start[0] = math.log((y.sum() + 1e-6) / (E.sum() + 1e-6))
    u_start = np.zeros(n)
    v_start = np.zeros(n)

    for tau_s in tau_s_grid:
        for tau_u in tau_u_grid:
            beta = beta_start.copy()
            u = u_start.copy()
            v = v_start.copy()
            current_obj = negative_log_posterior(y, E, X, u, v, beta, tau_s, tau_u, Q, sum_to_zero=sum_to_zero)

            for _ in range(40):
                gradient, hessian, _, _ = build_hessian_and_gradient(
                    y, E, X, beta, u, v, tau_s, tau_u, Q, sum_to_zero=sum_to_zero
                )
                try:
                    step = solve(hessian, gradient, assume_a='sym')
                except Exception:
                    step = np.linalg.solve(hessian, gradient)
                step_scale = 1.0
                improved = False
                theta = np.concatenate([beta, u, v])
                while step_scale > 1e-4:
                    proposal = theta - step_scale * step
                    beta_new = proposal[:p]
                    u_new = proposal[p:p + n]
                    v_new = proposal[p + n:]
                    u_new = u_new - u_new.mean()
                    proposal_obj = negative_log_posterior(
                        y, E, X, u_new, v_new, beta_new, tau_s, tau_u, Q, sum_to_zero=sum_to_zero
                    )
                    if np.isfinite(proposal_obj) and proposal_obj < current_obj:
                        beta, u, v = beta_new, u_new, v_new
                        current_obj = proposal_obj
                        improved = True
                        break
                    step_scale *= 0.5
                if not improved or np.max(np.abs(step_scale * step)) < 1e-5:
                    break

            gradient, hessian, eta, mu = build_hessian_and_gradient(
                y, E, X, beta, u, v, tau_s, tau_u, Q, sum_to_zero=sum_to_zero
            )
            sign, logdet = np.linalg.slogdet(hessian)
            if sign <= 0:
                continue
            laplace_score = current_obj + 0.5 * logdet
            if best is None or laplace_score < best['laplace_score']:
                covariance = np.linalg.inv(hessian)
                best = {
                    'disease': disease_name,
                    'tau_structured': tau_s,
                    'tau_unstructured': tau_u,
                    'beta': beta,
                    'u': u,
                    'v': v,
                    'eta': eta,
                    'mu': mu,
                    'gradient_norm': float(np.linalg.norm(gradient)),
                    'hessian': hessian,
                    'covariance': covariance,
                    'laplace_score': float(laplace_score),
                    'objective': float(current_obj),
                }

    if best is None:
        raise RuntimeError(f'Failed to fit BYM model for {disease_name}.')
    return best


def summarize_posterior(area_df: pd.DataFrame, X: np.ndarray, feature_names, fit: dict, disease: str):
    n = len(area_df)
    p = X.shape[1]
    cov = fit['covariance']
    rr_means = []
    rr_lower = []
    rr_upper = []
    exceed = []
    log_rr_sd = []
    design = np.hstack([X, np.eye(n), np.eye(n)])

    for i in range(n):
        d = design[i]
        mean_log_rr = fit['eta'][i]
        var_log_rr = float(d @ cov @ d)
        var_log_rr = max(var_log_rr, 1e-10)
        sd = math.sqrt(var_log_rr)
        log_rr_sd.append(sd)
        rr_means.append(math.exp(mean_log_rr + 0.5 * var_log_rr))
        rr_lower.append(math.exp(mean_log_rr - 1.96 * sd))
        rr_upper.append(math.exp(mean_log_rr + 1.96 * sd))
        exceed.append(float(norm.cdf(mean_log_rr / sd)))

    area_df[f'{disease}_rr_mean'] = rr_means
    area_df[f'{disease}_rr_lower'] = rr_lower
    area_df[f'{disease}_rr_upper'] = rr_upper
    area_df[f'{disease}_exceedance_prob'] = exceed
    area_df[f'{disease}_log_rr_sd'] = log_rr_sd

    beta_summary = []
    for idx, name in enumerate(feature_names):
        mean = float(fit['beta'][idx])
        sd = float(np.sqrt(max(cov[idx, idx], 1e-12)))
        beta_summary.append({
            'name': name,
            'posterior_mean': mean,
            'lower_95': mean - 1.96 * sd,
            'upper_95': mean + 1.96 * sd,
            'rr_multiplier': math.exp(mean),
            'rr_lower_95': math.exp(mean - 1.96 * sd),
            'rr_upper_95': math.exp(mean + 1.96 * sd),
        })
    return area_df, beta_summary


def save_map_figure(area_df: pd.DataFrame, value_columns, titles, output_path: Path, cmap: str, n_rows: int, n_cols: int, cbar_label: str, vmin=None, vmax=None):
    fig, axes = plt.subplots(1, len(value_columns), figsize=(6 * len(value_columns), 5), constrained_layout=True)
    if len(value_columns) == 1:
        axes = [axes]
    for ax, col, title in zip(axes, value_columns, titles):
        grid = area_df.pivot(index='row', columns='col', values=col).sort_index(ascending=True)
        im = ax.imshow(grid.to_numpy(), cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.set_xlabel('Grid column')
        ax.set_ylabel('Grid row')
        cbar = fig.colorbar(im, ax=ax, shrink=0.85)
        cbar.set_label(cbar_label)
    fig.savefig(output_path, dpi=320, bbox_inches='tight')
    plt.close(fig)


def plot_covariate_effects(beta_results: dict, output_path: Path):
    covariates = [x['name'] for x in beta_results['malaria'] if x['name'] != 'Intercept']
    y_positions = np.arange(len(covariates))
    offsets = {'malaria': -0.12, 'dengue': 0.12}
    colors = {'malaria': '#1f77b4', 'dengue': '#ff7f0e'}

    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    for disease in ['malaria', 'dengue']:
        effects = {item['name']: item for item in beta_results[disease] if item['name'] != 'Intercept'}
        means = [effects[name]['posterior_mean'] for name in covariates]
        lowers = [effects[name]['lower_95'] for name in covariates]
        uppers = [effects[name]['upper_95'] for name in covariates]
        xerr = np.vstack([np.array(means) - np.array(lowers), np.array(uppers) - np.array(means)])
        ax.errorbar(means, y_positions + offsets[disease], xerr=xerr, fmt='o', capsize=4, label=disease.title(), color=colors[disease])
    ax.axvline(0, color='black', linestyle='--', linewidth=1)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(covariates)
    ax.set_xlabel('Posterior mean log-relative risk (95% credible interval)')
    ax.set_ylabel('Covariate')
    ax.set_title('Covariate effects from BYM disease mapping models')
    ax.legend(frameon=False)
    fig.savefig(output_path, dpi=320, bbox_inches='tight')
    plt.close(fig)


def plot_temporal_trend(df: pd.DataFrame, output_path: Path):
    monthly = df.groupby('month_label', as_index=False)[['malaria_cases', 'dengue_cases']].sum()
    x = np.arange(len(monthly))
    fig, ax = plt.subplots(figsize=(11, 4.5), constrained_layout=True)
    ax.plot(x, monthly['malaria_cases'], marker='o', linewidth=2, color='#1f77b4', label='Malaria')
    ax.plot(x, monthly['dengue_cases'], marker='s', linewidth=2, color='#ff7f0e', label='Dengue')
    for disease, color, style in [('malaria_cases', '#1f77b4', '--'), ('dengue_cases', '#ff7f0e', ':')]:
        coef = np.polyfit(x, monthly[disease], deg=2)
        ax.plot(x, np.polyval(coef, x), linestyle=style, color=color, alpha=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(monthly['month_label'], rotation=45, ha='right')
    ax.set_xlabel('Month')
    ax.set_ylabel('Cases')
    ax.set_title('Monthly total reported cases')
    ax.legend(frameon=False)
    fig.savefig(output_path, dpi=320, bbox_inches='tight')
    plt.close(fig)


def plot_risk_comparison(area_df: pd.DataFrame, output_path: Path):
    fig, ax = plt.subplots(figsize=(6.5, 6), constrained_layout=True)
    sc = ax.scatter(
        area_df['malaria_rr_mean'],
        area_df['dengue_rr_mean'],
        c=area_df['urbanization_index'],
        cmap='viridis',
        s=55,
        edgecolor='black',
        linewidth=0.25,
        alpha=0.9,
    )
    lims = [
        min(area_df['malaria_rr_mean'].min(), area_df['dengue_rr_mean'].min()) * 0.95,
        max(area_df['malaria_rr_mean'].max(), area_df['dengue_rr_mean'].max()) * 1.05,
    ]
    ax.plot(lims, lims, linestyle='--', color='gray', linewidth=1)
    ax.axvline(1.0, linestyle=':', color='black', linewidth=1)
    ax.axhline(1.0, linestyle=':', color='black', linewidth=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel('Malaria posterior relative risk')
    ax.set_ylabel('Dengue posterior relative risk')
    ax.set_title('Cross-disease comparison of area-level risk')
    cbar = fig.colorbar(sc, ax=ax, shrink=0.85)
    cbar.set_label('Urbanization index')
    fig.savefig(output_path, dpi=320, bbox_inches='tight')
    plt.close(fig)


def create_results_json(area_df: pd.DataFrame, correlation_matrix: pd.DataFrame, model_results: dict, beta_results: dict, feature_names, output_path: Path):
    summary = {
        'metadata': {
            'created_at': datetime.now(timezone.utc).isoformat(),
            'seed': SEED,
            'n_areas': int(area_df.shape[0]),
            'n_months': 24,
            'model_description': 'Poisson BYM-inspired disease mapping with ICAR structured effects, Gaussian unstructured effects, and Laplace approximation over the latent posterior.',
        },
        'covariate_correlation': correlation_matrix.round(4).to_dict(),
        'diseases': {},
    }

    for disease in ['malaria', 'dengue']:
        rr_col = f'{disease}_rr_mean'
        exceed_col = f'{disease}_exceedance_prob'
        smr_col = f'{disease}_smr'
        cases_col = f'{disease}_cases'
        top_hotspots = area_df.nlargest(5, rr_col)[['area_id', rr_col, exceed_col, cases_col]].round(4).to_dict(orient='records')
        summary['diseases'][disease] = {
            'total_cases': int(area_df[cases_col].sum()),
            'overall_rate_per_100k_year': float((area_df[cases_col].sum() / area_df['population'].sum()) * 100000 / 2.0),
            'smr_summary': {
                'mean': float(area_df[smr_col].mean()),
                'median': float(area_df[smr_col].median()),
                'min': float(area_df[smr_col].min()),
                'max': float(area_df[smr_col].max()),
            },
            'rr_summary': {
                'mean': float(area_df[rr_col].mean()),
                'median': float(area_df[rr_col].median()),
                'min': float(area_df[rr_col].min()),
                'max': float(area_df[rr_col].max()),
            },
            'exceedance_summary': {
                'mean': float(area_df[exceed_col].mean()),
                'areas_above_0_8': int((area_df[exceed_col] > 0.8).sum()),
                'areas_above_0_95': int((area_df[exceed_col] > 0.95).sum()),
            },
            'model_parameters': {
                'tau_structured': float(model_results[disease]['tau_structured']),
                'tau_unstructured': float(model_results[disease]['tau_unstructured']),
                'laplace_score': float(model_results[disease]['laplace_score']),
                'objective': float(model_results[disease]['objective']),
                'gradient_norm': float(model_results[disease]['gradient_norm']),
            },
            'covariate_effects': {item['name']: item for item in beta_results[disease]},
            'top_hotspots': top_hotspots,
        }

    summary['comparative_metrics'] = {
        'rr_correlation': float(np.corrcoef(area_df['malaria_rr_mean'], area_df['dengue_rr_mean'])[0, 1]),
        'smr_correlation': float(np.corrcoef(area_df['malaria_smr'], area_df['dengue_smr'])[0, 1]),
        'high_risk_overlap_0_8': int(((area_df['malaria_exceedance_prob'] > 0.8) & (area_df['dengue_exceedance_prob'] > 0.8)).sum()),
    }

    with output_path.open('w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)


def print_summary(area_df: pd.DataFrame, beta_results: dict, model_results: dict):
    print('Disease risk mapping case study completed')
    print('=' * 72)
    print(f"Study area administrative units: {len(area_df)}")
    print('Observation months: 24')
    print(f"Population covered: {int(area_df['population'].sum()):,}")
    print('-' * 72)
    for disease in ['malaria', 'dengue']:
        cases = int(area_df[f'{disease}_cases'].sum())
        print(f"{disease.title()} total cases: {cases:,}")
        print(
            f"  Selected precisions -> structured={model_results[disease]['tau_structured']:.2f}, "
            f"unstructured={model_results[disease]['tau_unstructured']:.2f}"
        )
        print(
            f"  Posterior RR range -> {area_df[f'{disease}_rr_mean'].min():.3f} to "
            f"{area_df[f'{disease}_rr_mean'].max():.3f}"
        )
        print(
            f"  Areas with exceedance probability > 0.80 -> "
            f"{int((area_df[f'{disease}_exceedance_prob'] > 0.80).sum())}"
        )
        print('  Covariate posterior means (log-RR):')
        for item in beta_results[disease]:
            if item['name'] == 'Intercept':
                continue
            print(
                f"    {item['name']}: {item['posterior_mean']:+.3f} "
                f"[{item['lower_95']:+.3f}, {item['upper_95']:+.3f}]"
            )
        print('-' * 72)
    print(
        'Cross-disease RR correlation: '
        f"{np.corrcoef(area_df['malaria_rr_mean'], area_df['dengue_rr_mean'])[0, 1]:.3f}"
    )


def main():
    base_dir = Path(__file__).resolve().parents[1]
    paths = ensure_directories(base_dir)

    coords, adjacency, precision, _ = build_grid(14, 14)
    disease_df = simulate_disease_data(coords, precision, n_months=24)
    data_output = paths['data'] / 'synthetic_disease_data.csv'
    disease_df.to_csv(data_output, index=False)

    area_df = aggregate_area_data(disease_df)
    X, feature_names, means, stds = make_design_matrix(area_df)

    correlation_matrix = area_df[[
        'temperature_mean',
        'precipitation_mean',
        'elevation',
        'population_density',
        'urbanization_index',
        'malaria_rate_per_100k',
        'dengue_rate_per_100k',
    ]].corr()

    model_results = {}
    beta_results = {}
    for disease in ['malaria', 'dengue']:
        fit = fit_bym_laplace(
            y=area_df[f'{disease}_cases'].to_numpy(dtype=float),
            E=np.clip(area_df[f'{disease}_expected'].to_numpy(dtype=float), 1e-6, None),
            X=X,
            Q=precision,
            disease_name=disease,
        )
        model_results[disease] = fit
        area_df, beta_summary = summarize_posterior(area_df, X, feature_names, fit, disease)
        beta_results[disease] = beta_summary

    save_map_figure(
        area_df,
        ['malaria_smr', 'dengue_smr'],
        ['Malaria standardized morbidity ratio', 'Dengue standardized morbidity ratio'],
        paths['figures'] / 'disease_smr_map.png',
        cmap='RdYlGn_r',
        n_rows=14,
        n_cols=14,
        cbar_label='SMR',
    )
    save_map_figure(
        area_df,
        ['malaria_rr_mean', 'dengue_rr_mean'],
        ['Malaria posterior relative risk', 'Dengue posterior relative risk'],
        paths['figures'] / 'disease_rr_map.png',
        cmap='RdYlGn_r',
        n_rows=14,
        n_cols=14,
        cbar_label='Posterior RR',
    )
    save_map_figure(
        area_df,
        ['malaria_exceedance_prob', 'dengue_exceedance_prob'],
        ['Malaria exceedance probability P(RR > 1)', 'Dengue exceedance probability P(RR > 1)'],
        paths['figures'] / 'disease_exceedance_prob.png',
        cmap='viridis',
        n_rows=14,
        n_cols=14,
        cbar_label='Exceedance probability',
        vmin=0.0,
        vmax=1.0,
    )
    plot_covariate_effects(beta_results, paths['figures'] / 'disease_covariate_effects.png')
    plot_temporal_trend(disease_df, paths['figures'] / 'disease_temporal_trend.png')
    plot_risk_comparison(area_df, paths['figures'] / 'disease_risk_comparison.png')

    create_results_json(
        area_df=area_df,
        correlation_matrix=correlation_matrix,
        model_results=model_results,
        beta_results=beta_results,
        feature_names=feature_names,
        output_path=paths['results'] / 'disease_mapping_results.json',
    )

    area_df.to_csv(paths['data'] / 'synthetic_disease_area_summary.csv', index=False)

    print_summary(area_df, beta_results, model_results)


if __name__ == '__main__':
    main()
