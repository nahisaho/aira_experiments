#!/usr/bin/env python3
import json
import math
import random
import warnings
from datetime import datetime, timezone
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp
    from scipy.optimize import minimize
    from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import mean_squared_error, r2_score
    from sklearn.model_selection import KFold, cross_val_predict, cross_val_score, train_test_split
    from sklearn.neural_network import MLPRegressor
    from sklearn.preprocessing import StandardScaler
except ImportError as exc:
    raise SystemExit(f"Missing dependency: {exc}. Please install required packages before running.")

plt.switch_backend('Agg')
plt.style.use('seaborn-v0_8-colorblind')
warnings.filterwarnings('ignore', category=ConvergenceWarning)
warnings.filterwarnings('ignore', message='X has feature names')

SEED = 42
np.random.seed(SEED)
random.seed(SEED)
try:
    import torch
    torch.manual_seed(SEED)
except Exception:
    pass
try:
    import tensorflow as tf
    tf.random.set_seed(SEED)
except Exception:
    pass

BASE_DIR = Path(__file__).resolve().parent
FIG_DIR = BASE_DIR / 'figures'
RESULTS_DIR = BASE_DIR / 'results'
DATA_DIR = BASE_DIR / 'data'
LOG_DIR = BASE_DIR / 'logs'
for directory in [FIG_DIR, RESULTS_DIR, DATA_DIR, LOG_DIR]:
    directory.mkdir(parents=True, exist_ok=True)
LOG_FILE = LOG_DIR / 'process-log.jsonl'
R = 8.314


def now_iso():
    return datetime.now(timezone.utc).isoformat()


def log_event(phase, event_type, skill_or_tool, handoff_in=None, handoff_out=None, files_written=None, status='ok'):
    record = {
        'timestamp': now_iso(),
        'phase': phase,
        'event_type': event_type,
        'actor': 'co-scientist',
        'skill_or_tool': skill_or_tool,
        'handoff_in': handoff_in or {},
        'handoff_out': handoff_out or {},
        'files_written': [str(Path(f).relative_to(BASE_DIR)) if Path(f).is_absolute() else str(f) for f in (files_written or [])],
        'status': status,
    }
    with LOG_FILE.open('a', encoding='utf-8') as fh:
        fh.write(json.dumps(record, ensure_ascii=False) + '\n')


def save_csv(df, filename):
    path = RESULTS_DIR / filename
    df.to_csv(path, index=False)
    log_event('report', 'file_written', 'pandas.to_csv', files_written=[path], handoff_out={'rows': int(len(df)), 'columns': list(df.columns)})
    return path


def save_data_csv(df, filename):
    path = DATA_DIR / filename
    df.to_csv(path, index=False)
    log_event('data', 'file_written', 'pandas.to_csv', files_written=[path], handoff_out={'rows': int(len(df)), 'columns': list(df.columns)})
    return path


def save_fig(filename):
    path = FIG_DIR / filename
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    log_event('visualization', 'file_written', 'matplotlib', files_written=[path])
    return path


def get_polymer_map():
    return {
        'PLA': {'bond_type': 'ester', 'base_crystallinity': 28, 'contact_angle': 78, 'flexibility': 4, 'symmetry': 0.72, 'bde': 82, 'family': 'aliphatic polyester'},
        'PHA': {'bond_type': 'ester', 'base_crystallinity': 52, 'contact_angle': 84, 'flexibility': 5, 'symmetry': 0.68, 'bde': 79, 'family': 'microbial polyester'},
        'PBS': {'bond_type': 'ester', 'base_crystallinity': 38, 'contact_angle': 81, 'flexibility': 6, 'symmetry': 0.64, 'bde': 80, 'family': 'succinate polyester'},
        'PCL': {'bond_type': 'ester', 'base_crystallinity': 48, 'contact_angle': 86, 'flexibility': 8, 'symmetry': 0.58, 'bde': 77, 'family': 'caprolactone polyester'},
        'PGA': {'bond_type': 'ester', 'base_crystallinity': 55, 'contact_angle': 72, 'flexibility': 3, 'symmetry': 0.81, 'bde': 84, 'family': 'glycolide polyester'},
        'PAA': {'bond_type': 'anhydride', 'base_crystallinity': 12, 'contact_angle': 66, 'flexibility': 4, 'symmetry': 0.40, 'bde': 62, 'family': 'polyanhydride'},
        'PU-ester': {'bond_type': 'urethane', 'base_crystallinity': 22, 'contact_angle': 88, 'flexibility': 7, 'symmetry': 0.51, 'bde': 92, 'family': 'polyurethane'},
        'PA-amide': {'bond_type': 'amide', 'base_crystallinity': 33, 'contact_angle': 70, 'flexibility': 3, 'symmetry': 0.76, 'bde': 96, 'family': 'polyamide'},
        'PTMC': {'bond_type': 'carbonate', 'base_crystallinity': 18, 'contact_angle': 74, 'flexibility': 6, 'symmetry': 0.47, 'bde': 74, 'family': 'polycarbonate'},
        'PHU': {'bond_type': 'urethane', 'base_crystallinity': 15, 'contact_angle': 69, 'flexibility': 5, 'symmetry': 0.55, 'bde': 90, 'family': 'hydroxyurethane'},
    }


def encode_bond_type(df):
    return pd.get_dummies(df, columns=['bond_type'], prefix='bond')


def hydrolysis_rate_from_features(row):
    logA_map = {'ester': 7.0, 'amide': 5.2, 'carbonate': 6.1, 'anhydride': 9.0, 'urethane': 5.6}
    ea_map = {'ester': 58_000, 'amide': 72_000, 'carbonate': 63_000, 'anhydride': 46_000, 'urethane': 68_000}
    hydro_map = {'ester': 0.16, 'amide': 0.08, 'carbonate': 0.12, 'anhydride': 0.24, 'urethane': 0.07}
    A = math.exp(logA_map[row['bond_type']])
    temp_k = row['temperature_C'] + 273.15
    arrhenius = A * math.exp(-ea_map[row['bond_type']] / (R * temp_k))
    crystal_factor = math.exp(-0.030 * row['crystallinity_pct'])
    mw_factor = (row['Mw'] / 100000.0) ** -0.24
    pdi_factor = 1.0 + 0.12 * (row['PDI'] - 1.4)
    wetting_factor = 1.0 + hydro_map[row['bond_type']] * ((95.0 - row['contact_angle_deg']) / 35.0)
    descriptor_factor = 1.0 + 0.02 * row['flexibility_proxy'] - 0.09 * row['symmetry_index'] + 0.003 * (85.0 - row['bde_proxy'])
    kh = 1200.0 * arrhenius * crystal_factor * mw_factor * pdi_factor * wetting_factor * descriptor_factor
    return max(kh, 1e-6)


def generate_synthetic_dataset(n_samples=220):
    log_event('plan', 'handoff_started', 'synthetic_dataset', handoff_in={'n_samples': n_samples})
    polymer_map = get_polymer_map()
    polymers = list(polymer_map.keys())
    rows = []
    for idx in range(n_samples):
        polymer = random.choice(polymers)
        info = polymer_map[polymer]
        crystallinity = float(np.clip(np.random.normal(info['base_crystallinity'], 10), 3, 70))
        mn = float(np.random.uniform(35_000, 240_000))
        pdi = float(np.clip(np.random.normal(1.9, 0.35), 1.1, 3.0))
        mw = float(mn * pdi)
        contact_angle = float(np.clip(np.random.normal(info['contact_angle'], 8), 55, 98))
        temperature = float(np.random.uniform(10, 55))
        flexibility = max(1, int(round(np.clip(np.random.normal(info['flexibility'], 1.2), 1, 10))))
        symmetry = float(np.clip(np.random.normal(info['symmetry'], 0.08), 0.25, 0.95))
        bde = float(np.clip(np.random.normal(info['bde'], 5), 58, 102))
        hydrophilicity_idx = float(np.clip((100 - contact_angle) / 100 + np.random.normal(0, 0.03), 0.02, 0.55))
        row = {
            'sample_id': f'S{idx+1:03d}',
            'polymer': polymer,
            'family': info['family'],
            'bond_type': info['bond_type'],
            'crystallinity_pct': crystallinity,
            'Mn': mn,
            'Mw': mw,
            'PDI': pdi,
            'contact_angle_deg': contact_angle,
            'temperature_C': temperature,
            'flexibility_proxy': flexibility,
            'hydrophilicity_index': hydrophilicity_idx,
            'symmetry_index': symmetry,
            'bde_proxy': bde,
        }
        kh = hydrolysis_rate_from_features(row)
        kh *= float(np.exp(np.random.normal(0, 0.18)))
        row['k_h'] = max(kh, 1e-6)
        row['half_life_days'] = math.log(2) / row['k_h']
        rows.append(row)
    df = pd.DataFrame(rows)
    save_data_csv(df, 'synthetic_polymer_dataset.csv')
    log_event('execute', 'handoff_completed', 'synthetic_dataset', handoff_out={'samples': len(df), 'bond_types': df['bond_type'].value_counts().to_dict()}, files_written=[DATA_DIR / 'synthetic_polymer_dataset.csv'])
    return df


def train_hydrolysis_model(df):
    feature_cols = ['crystallinity_pct', 'Mn', 'Mw', 'PDI', 'contact_angle_deg', 'temperature_C', 'flexibility_proxy', 'hydrophilicity_index', 'symmetry_index', 'bde_proxy', 'bond_type']
    model_df = encode_bond_type(df[feature_cols + ['k_h']].copy())
    X = model_df.drop(columns=['k_h'])
    y = model_df['k_h']
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=SEED)
    model = RandomForestRegressor(n_estimators=400, random_state=SEED, min_samples_leaf=2)
    model.fit(X_train, y_train)
    pred = model.predict(X_test)
    metrics = pd.DataFrame([
        {'metric': 'test_rmse', 'value': float(np.sqrt(mean_squared_error(y_test, pred)))},
        {'metric': 'test_r2', 'value': float(r2_score(y_test, pred))},
        {'metric': 'cv_rmse_mean', 'value': float((-cross_val_score(model, X, y, cv=5, scoring='neg_root_mean_squared_error')).mean())},
        {'metric': 'cv_r2_mean', 'value': float(cross_val_score(model, X, y, cv=5, scoring='r2').mean())},
    ])
    save_csv(metrics, 'hydrolysis_model_metrics.csv')
    importance_df = pd.DataFrame({'feature': X.columns, 'importance': model.feature_importances_}).sort_values('importance', ascending=False)
    save_csv(importance_df, 'hydrolysis_feature_importance.csv')

    plt.figure(figsize=(6.4, 4.8))
    plt.scatter(y_test, pred, alpha=0.75, edgecolor='black', linewidth=0.3)
    lims = [min(y_test.min(), pred.min()), max(y_test.max(), pred.max())]
    plt.plot(lims, lims, '--', color='tab:red')
    plt.xlabel('Observed hydrolysis rate k_h (day$^{-1}$)')
    plt.ylabel('Predicted hydrolysis rate k_h (day$^{-1}$)')
    plt.title('Random Forest hydrolysis model')
    save_fig('hydrolysis_model_parity.png')

    plt.figure(figsize=(6.8, 4.8))
    top_imp = importance_df.head(10).iloc[::-1]
    plt.barh(top_imp['feature'], top_imp['importance'], color='tab:blue')
    plt.xlabel('Feature importance')
    plt.ylabel('Descriptor')
    plt.title('Hydrolysis model feature importance')
    save_fig('hydrolysis_feature_importance.png')

    return model, list(X.columns), metrics, importance_df


def predict_kh_for_design(crystallinity, mw, contact_angle=78, temperature=25, bond_type='ester'):
    base = hydrolysis_rate_from_features({
        'bond_type': bond_type,
        'crystallinity_pct': crystallinity,
        'Mw': mw,
        'PDI': 1.8,
        'contact_angle_deg': contact_angle,
        'temperature_C': temperature,
        'flexibility_proxy': 5,
        'symmetry_index': 0.62,
        'bde_proxy': 80,
    })
    return base


def optimize_tradeoff():
    log_event('execute', 'handoff_started', 'tradeoff_optimization')
    alpha, beta = 18.0, 1200.0
    gamma, delta = 1.8, 180.0

    def objectives(x):
        xc, mw = x
        kh = predict_kh_for_design(xc, mw, contact_angle=80)
        sigma_t = alpha * (xc ** 0.5) * ((mw / 1000.0) ** 0.3) - beta * kh
        modulus = gamma * xc * ((mw / 1000.0) ** 0.2) - delta * kh
        return kh, sigma_t, modulus

    def scalar_objective(x, w1, w2, w3):
        kh, sigma_t, modulus = objectives(x)
        kh_n = kh / 0.02
        sigma_n = sigma_t / 250
        modulus_n = modulus / 220
        return -(w1 * kh_n + w2 * sigma_n + w3 * modulus_n)

    sols = []
    weights = np.linspace(0.05, 0.9, 12)
    for w1 in weights:
        for w2 in np.linspace(0.05, 0.9, 12):
            w3 = max(0.05, 1.0 - w1 - w2)
            total = w1 + w2 + w3
            w1n, w2n, w3n = w1 / total, w2 / total, w3 / total
            start = np.array([np.random.uniform(5.0, 65.0), np.random.uniform(50000.0, 300000.0)])
            res = minimize(scalar_objective, x0=start, args=(w1n, w2n, w3n), bounds=[(5.0, 65.0), (50000.0, 300000.0)], method='L-BFGS-B')
            xc, mw = res.x
            kh, sigma_t, modulus = objectives(res.x)
            sols.append({'source': 'scipy_minimize', 'crystallinity_pct': xc, 'Mw': mw, 'k_h': kh, 'tensile_strength_MPa': sigma_t, 'modulus_GPa': modulus})

    for xc in np.linspace(5.0, 65.0, 61):
        for mw in np.linspace(50000.0, 300000.0, 80):
            kh, sigma_t, modulus = objectives((xc, mw))
            sols.append({'source': 'grid_screen', 'crystallinity_pct': xc, 'Mw': mw, 'k_h': kh, 'tensile_strength_MPa': sigma_t, 'modulus_GPa': modulus})

    sol_df = pd.DataFrame(sols).drop_duplicates(subset=['crystallinity_pct', 'Mw', 'k_h', 'tensile_strength_MPa', 'modulus_GPa']).reset_index(drop=True)
    metrics = sol_df[['k_h', 'tensile_strength_MPa', 'modulus_GPa']]
    dominated = []
    for i, row in metrics.iterrows():
        is_dominated = ((metrics['k_h'] >= row['k_h']) & (metrics['tensile_strength_MPa'] >= row['tensile_strength_MPa']) & (metrics['modulus_GPa'] >= row['modulus_GPa']) & ((metrics['k_h'] > row['k_h']) | (metrics['tensile_strength_MPa'] > row['tensile_strength_MPa']) | (metrics['modulus_GPa'] > row['modulus_GPa']))).any()
        dominated.append(is_dominated)
    pareto = sol_df.loc[~pd.Series(dominated)].sort_values(['k_h', 'tensile_strength_MPa']).reset_index(drop=True)
    save_csv(sol_df, 'tradeoff_candidates.csv')
    save_csv(pareto, 'pareto_front.csv')

    plt.figure(figsize=(6.8, 4.8))
    plt.scatter(sol_df['k_h'], sol_df['tensile_strength_MPa'], c=sol_df['modulus_GPa'], cmap='viridis', alpha=0.45)
    plt.plot(pareto['k_h'], pareto['tensile_strength_MPa'], color='tab:red', marker='o', linewidth=1.4, label='Pareto front')
    plt.xlabel('Hydrolysis rate k_h (day$^{-1}$)')
    plt.ylabel('Tensile strength (MPa)')
    cbar = plt.colorbar()
    cbar.set_label('Modulus (GPa)')
    plt.title('Degradability-mechanical tradeoff')
    plt.legend()
    save_fig('pareto_tradeoff.png')
    log_event('execute', 'handoff_completed', 'tradeoff_optimization', handoff_out={'pareto_points': int(len(pareto))}, files_written=[RESULTS_DIR / 'pareto_front.csv'])
    return pareto


def simulate_microbial_degradation():
    log_event('execute', 'handoff_started', 'microbial_degradation')
    Vmax = 0.42
    Km = 1.1
    Kd = 0.8
    tau_biofilm = 3.5
    times = np.linspace(0, 30, 250)
    records = []
    for e_total in [0.5, 1.0, 2.0]:
        for s0 in [1.0, 2.5, 5.0]:
            def rhs(t, y):
                s = max(y[0], 0.0)
                e_active = e_total * (1.0 - math.exp(-t / tau_biofilm))
                rate = Vmax * s / (Km + s) * e_active / (Kd + s)
                return [-rate]
            sol = solve_ivp(rhs, (times.min(), times.max()), [s0], t_eval=times, method='RK45')
            for t, s in zip(sol.t, sol.y[0]):
                e_active = e_total * (1.0 - math.exp(-t / tau_biofilm))
                rate = Vmax * max(s, 0.0) / (Km + max(s, 0.0)) * e_active / (Kd + max(s, 0.0))
                records.append({'time_day': float(t), 'enzyme_total': e_total, 'substrate_initial': s0, 'substrate_conc': float(max(s, 0.0)), 'enzyme_active': e_active, 'instantaneous_rate': rate})
    df = pd.DataFrame(records)
    save_csv(df, 'microbial_degradation_profiles.csv')
    plt.figure(figsize=(7.2, 5.0))
    for (e_total, s0), sub in df.groupby(['enzyme_total', 'substrate_initial']):
        plt.plot(sub['time_day'], sub['substrate_conc'], label=f'E={e_total}, S0={s0}')
    plt.xlabel('Time (day)')
    plt.ylabel('Substrate concentration')
    plt.title('Microbial degradation with biofilm activation')
    plt.legend(ncol=3, fontsize=7)
    save_fig('microbial_degradation_profiles.png')
    log_event('execute', 'handoff_completed', 'microbial_degradation', handoff_out={'profiles': int(df.groupby(['enzyme_total', 'substrate_initial']).ngroups)}, files_written=[RESULTS_DIR / 'microbial_degradation_profiles.csv'])
    return df


def simulate_marine_environment():
    log_event('execute', 'handoff_started', 'marine_simulation')
    polymer_params = {
        'PLA': {'k0': 0.00008, 'alpha_oh': 2.4e5, 'k_enz0': 0.00030, 'diversity': 2.3, 'fraction': 0.14, 'uv': 0.10, 'q10': 1.7},
        'PHA': {'k0': 0.00011, 'alpha_oh': 3.0e5, 'k_enz0': 0.00048, 'diversity': 2.7, 'fraction': 0.19, 'uv': 0.08, 'q10': 1.9},
        'PBS': {'k0': 0.00007, 'alpha_oh': 2.1e5, 'k_enz0': 0.00024, 'diversity': 2.2, 'fraction': 0.11, 'uv': 0.12, 'q10': 1.6},
    }
    days = np.linspace(0, 365, 366)
    rows = []
    for polymer, prm in polymer_params.items():
        def rhs(t, y):
            month_temp = 17 + 7 * math.sin(2 * math.pi * t / 365 - 0.5)
            pH = 8.05 + 0.18 * math.sin(2 * math.pi * t / 365 + 0.3)
            salinity = 32.5 + 1.5 * math.cos(2 * math.pi * t / 365)
            uv_factor = 1.0 + prm['uv'] * (1 + math.sin(2 * math.pi * t / 365 - 1.2))
            oh = 10 ** (pH - 14)
            k_abio = prm['k0'] * (1 + prm['alpha_oh'] * oh) * uv_factor
            f_microbe = (prm['diversity'] / 3.0) * prm['fraction'] * (1 + 0.01 * (salinity - 32.5))
            f_T = prm['q10'] ** ((month_temp - 20) / 10)
            k_enz = prm['k_enz0'] * f_microbe * f_T
            return [-(k_abio + k_enz) * y[0]]
        sol = solve_ivp(rhs, (0, 365), [100.0], t_eval=days, method='RK45')
        for t, w in zip(sol.t, sol.y[0]):
            month_temp = 17 + 7 * math.sin(2 * math.pi * t / 365 - 0.5)
            pH = 8.05 + 0.18 * math.sin(2 * math.pi * t / 365 + 0.3)
            salinity = 32.5 + 1.5 * math.cos(2 * math.pi * t / 365)
            oh = 10 ** (pH - 14)
            uv_factor = 1.0 + prm['uv'] * (1 + math.sin(2 * math.pi * t / 365 - 1.2))
            k_abio = prm['k0'] * (1 + prm['alpha_oh'] * oh) * uv_factor
            f_microbe = (prm['diversity'] / 3.0) * prm['fraction'] * (1 + 0.01 * (salinity - 32.5))
            f_T = prm['q10'] ** ((month_temp - 20) / 10)
            k_enz = prm['k_enz0'] * f_microbe * f_T
            weight_remaining = float(np.clip(w, 0.0, 100.0))
            rows.append({'polymer': polymer, 'day': float(t), 'temperature_C': month_temp, 'pH': pH, 'salinity_ppt': salinity, 'weight_remaining_pct': weight_remaining, 'weight_loss_pct': float(100.0 - weight_remaining), 'k_abio': k_abio, 'k_enz': k_enz, 'k_total': k_abio + k_enz})
    df = pd.DataFrame(rows)
    save_csv(df, 'marine_degradation_profiles.csv')
    plt.figure(figsize=(7.2, 5.0))
    for polymer, sub in df.groupby('polymer'):
        plt.plot(sub['day'], sub['weight_loss_pct'], label=polymer, linewidth=2)
    plt.xlabel('Time (day)')
    plt.ylabel('Weight loss (%)')
    plt.title('Marine degradation under seasonal forcing')
    plt.legend()
    save_fig('marine_degradation_profiles.png')
    log_event('execute', 'handoff_completed', 'marine_simulation', handoff_out={'polymers': list(polymer_params)}, files_written=[RESULTS_DIR / 'marine_degradation_profiles.csv'])
    return df


def ternary_to_cartesian(a, b, c):
    x = 0.5 * (2 * b + c) / (a + b + c)
    y = (math.sqrt(3) / 2) * c / (a + b + c)
    return x, y


def compositional_objective(fracs, monomer_props, chi, triad):
    names = list(triad)
    kh_base = sum(fracs[i] * monomer_props[names[i]]['k_h'] for i in range(len(names)))
    interaction = 0.0
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            pair = tuple(sorted((names[i], names[j])))
            interaction += fracs[i] * fracs[j] * chi.get(pair, 0.0) * 0.006
    kh = kh_base + interaction
    strength = sum(fracs[i] * monomer_props[names[i]]['strength'] for i in range(len(names))) - 8.0 * sum(fracs[i] * fracs[j] for i in range(len(names)) for j in range(i + 1, len(names)))
    modulus = sum(fracs[i] * monomer_props[names[i]]['modulus'] for i in range(len(names))) - 0.4 * sum(fracs[i] * fracs[j] for i in range(len(names)) for j in range(i + 1, len(names)))
    score = 0.50 * (kh / 0.025) + 0.30 * (strength / 60.0) + 0.20 * (modulus / 2.0)
    return kh, strength, modulus, score


def combinatorial_exploration():
    log_event('execute', 'handoff_started', 'copolymer_design')
    monomer_props = {
        'LA': {'k_h': 0.0048, 'strength': 58, 'modulus': 2.2},
        'GA': {'k_h': 0.0105, 'strength': 72, 'modulus': 2.8},
        'ε-CL': {'k_h': 0.0022, 'strength': 32, 'modulus': 0.45},
        'HB': {'k_h': 0.0068, 'strength': 40, 'modulus': 1.8},
        'BS': {'k_h': 0.0035, 'strength': 30, 'modulus': 0.65},
        'SA': {'k_h': 0.0056, 'strength': 36, 'modulus': 0.95},
    }
    chi = {
        ('GA', 'LA'): 0.9, ('LA', 'ε-CL'): -0.3, ('GA', 'ε-CL'): 0.4,
        ('HB', 'LA'): 0.2, ('BS', 'SA'): 0.5, ('HB', 'SA'): 0.1,
        ('BS', 'LA'): -0.1, ('GA', 'HB'): 0.3,
    }
    systems = [('LA', 'GA'), ('LA', 'ε-CL'), ('HB', 'BS', 'SA'), ('LA', 'GA', 'ε-CL')]
    grid_rows = []
    step = 0.05
    for system in systems:
        if len(system) == 2:
            for x1 in np.arange(0, 1 + step / 2, step):
                fracs = [float(round(x1, 4)), float(round(1 - x1, 4))]
                kh, strength, modulus, score = compositional_objective(fracs, monomer_props, chi, system)
                grid_rows.append({'system': '-'.join(system), 'monomer_1': system[0], 'monomer_2': system[1], 'monomer_3': '', 'x1': fracs[0], 'x2': fracs[1], 'x3': 0.0, 'k_h_copol': kh, 'strength_MPa': strength, 'modulus_GPa': modulus, 'score': score})
        else:
            for x1 in np.arange(0, 1 + step / 2, step):
                for x2 in np.arange(0, 1 - x1 + step / 2, step):
                    x3 = round(1 - x1 - x2, 6)
                    if x3 < -1e-9:
                        continue
                    fracs = [float(round(x1, 4)), float(round(x2, 4)), float(round(max(0.0, x3), 4))]
                    kh, strength, modulus, score = compositional_objective(fracs, monomer_props, chi, system)
                    grid_rows.append({'system': '-'.join(system), 'monomer_1': system[0], 'monomer_2': system[1], 'monomer_3': system[2], 'x1': fracs[0], 'x2': fracs[1], 'x3': fracs[2], 'k_h_copol': kh, 'strength_MPa': strength, 'modulus_GPa': modulus, 'score': score})
    grid_df = pd.DataFrame(grid_rows)
    save_csv(grid_df, 'copolymer_grid_search.csv')

    target_system = 'LA-GA-ε-CL'
    ternary = grid_df[grid_df['system'] == target_system].copy().reset_index(drop=True)
    feature_xy = ternary[['x1', 'x2']].values
    y = ternary['score'].values
    kernel = ConstantKernel(1.0, (0.1, 10.0)) * Matern(length_scale=0.2, nu=2.5) + WhiteKernel(noise_level=1e-5)
    gp = GaussianProcessRegressor(kernel=kernel, random_state=SEED, normalize_y=True)
    init_idx = np.linspace(0, len(ternary) - 1, 15, dtype=int)
    observed_idx = set(init_idx.tolist())
    history = []
    best_score = -np.inf
    for step_i in range(12):
        train_idx = sorted(observed_idx)
        gp.fit(feature_xy[train_idx], y[train_idx])
        mean, std = gp.predict(feature_xy, return_std=True)
        acquisition = mean + 1.2 * std
        acquisition[list(observed_idx)] = -np.inf
        next_idx = int(np.argmax(acquisition))
        observed_idx.add(next_idx)
        if y[next_idx] > best_score:
            best_score = float(y[next_idx])
        history.append({'iteration': step_i + 1, 'x1_LA': ternary.loc[next_idx, 'x1'], 'x2_GA': ternary.loc[next_idx, 'x2'], 'x3_eCL': ternary.loc[next_idx, 'x3'], 'predicted_mean': float(mean[next_idx]), 'predicted_std': float(std[next_idx]), 'observed_score': float(y[next_idx]), 'best_score_so_far': best_score})
    history_df = pd.DataFrame(history)
    save_csv(history_df, 'copolymer_gp_optimization.csv')
    top = ternary.nlargest(30, 'score').copy()
    coords = np.array([ternary_to_cartesian(a, b, c) for a, b, c in ternary[['x1', 'x2', 'x3']].values])
    top_coords = np.array([ternary_to_cartesian(a, b, c) for a, b, c in top[['x1', 'x2', 'x3']].values])
    ternary['plot_x'] = coords[:, 0]
    ternary['plot_y'] = coords[:, 1]
    top['plot_x'] = top_coords[:, 0]
    top['plot_y'] = top_coords[:, 1]
    save_csv(top[['system', 'x1', 'x2', 'x3', 'k_h_copol', 'strength_MPa', 'modulus_GPa', 'score']], 'copolymer_top_candidates.csv')

    plt.figure(figsize=(6.2, 5.8))
    plt.scatter(ternary['plot_x'], ternary['plot_y'], c=ternary['score'], cmap='viridis', s=26)
    plt.scatter(top['plot_x'], top['plot_y'], facecolors='none', edgecolors='red', s=48, linewidth=0.8)
    plt.plot([0, 1, 0.5, 0], [0, 0, math.sqrt(3) / 2, 0], color='black', linewidth=1)
    plt.text(-0.03, -0.03, 'LA', fontsize=11)
    plt.text(1.01, -0.03, 'GA', fontsize=11)
    plt.text(0.48, math.sqrt(3) / 2 + 0.03, 'ε-CL', fontsize=11)
    plt.xlabel('Ternary axis X')
    plt.ylabel('Ternary axis Y')
    plt.title('LA-GA-ε-CL compositional landscape')
    cbar = plt.colorbar()
    cbar.set_label('Composite score')
    save_fig('ternary_copolymer_landscape.png')
    log_event('execute', 'handoff_completed', 'copolymer_design', handoff_out={'grid_points': int(len(grid_df)), 'gp_iterations': int(len(history_df))}, files_written=[RESULTS_DIR / 'copolymer_grid_search.csv', RESULTS_DIR / 'copolymer_gp_optimization.csv'])
    return grid_df, history_df, top


def case_study_profiles():
    log_event('execute', 'handoff_started', 'case_studies')
    days = np.linspace(0, 365, 180)
    case_defs = [
        {'polymer': 'PLA', 'variant': 'Standard', 'k': predict_kh_for_design(25, 150000, 78), 'note': 'baseline'},
        {'polymer': 'PLA', 'variant': 'High-D stereodefect', 'k': predict_kh_for_design(18, 140000, 74) * 1.28, 'note': 'higher D-unit content'},
        {'polymer': 'PLA', 'variant': 'Nucleated high crystallinity', 'k': predict_kh_for_design(42, 180000, 80) * 0.76, 'note': 'nucleating agent'},
        {'polymer': 'PHA', 'variant': 'PHB', 'k': predict_kh_for_design(55, 220000, 84) * 0.78, 'note': 'high Tm homopolymer'},
        {'polymer': 'PHA', 'variant': 'PHBV 20% HV', 'k': predict_kh_for_design(35, 180000, 79) * 1.35, 'note': 'HV lowers crystallinity'},
        {'polymer': 'PHA', 'variant': 'PHBV blend', 'k': predict_kh_for_design(28, 160000, 76) * 1.52, 'note': 'blend + compatibilizer'},
        {'polymer': 'PBS', 'variant': 'Standard', 'k': predict_kh_for_design(38, 170000, 81), 'note': 'baseline'},
        {'polymer': 'PBS', 'variant': 'Branched PBS', 'k': predict_kh_for_design(30, 150000, 78) * 1.22, 'note': 'branching increases water uptake'},
        {'polymer': 'PBS', 'variant': 'Chain-extended PBS', 'k': predict_kh_for_design(45, 230000, 84) * 0.71, 'note': 'chain extender'},
    ]
    rows = []
    for case in case_defs:
        for day in days:
            remaining = 100 * math.exp(-case['k'] * day)
            rows.append({'polymer': case['polymer'], 'variant': case['variant'], 'day': float(day), 'weight_remaining_pct': remaining, 'weight_loss_pct': 100 - remaining, 'apparent_k_h': case['k'], 'strategy_note': case['note']})
    df = pd.DataFrame(rows)
    save_csv(df, 'case_study_profiles.csv')
    summary = df.groupby(['polymer', 'variant', 'strategy_note'], as_index=False).agg(final_weight_loss_pct=('weight_loss_pct', 'max'), apparent_k_h=('apparent_k_h', 'first'))
    save_csv(summary, 'case_study_summary.csv')
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8), sharey=True)
    for ax, polymer in zip(axes, ['PLA', 'PHA', 'PBS']):
        subset = df[df['polymer'] == polymer]
        for variant, sub in subset.groupby('variant'):
            ax.plot(sub['day'], sub['weight_loss_pct'], label=variant)
        ax.set_title(f'{polymer} case study')
        ax.set_xlabel('Time (day)')
        ax.grid(alpha=0.2)
    axes[0].set_ylabel('Weight loss (%)')
    axes[-1].legend(fontsize=7, loc='upper left', bbox_to_anchor=(1.02, 1.0))
    save_fig('case_study_profiles.png')
    log_event('execute', 'handoff_completed', 'case_studies', handoff_out={'variants': len(case_defs)}, files_written=[RESULTS_DIR / 'case_study_summary.csv'])
    return df, summary


def df_to_markdown(df):
    cols = list(df.columns)
    header = '| ' + ' | '.join(cols) + ' |'
    divider = '| ' + ' | '.join(['---'] * len(cols)) + ' |'
    rows = []
    for _, row in df.iterrows():
        values = []
        for col in cols:
            val = row[col]
            if isinstance(val, (float, np.floating)):
                values.append(f'{float(val):.5f}')
            else:
                values.append(str(val))
        rows.append('| ' + ' | '.join(values) + ' |')
    return '\n'.join([header, divider] + rows)


def shap_style_importance(model, X, feature_names):
    X_arr = np.asarray(X, dtype=float)
    baseline = model.predict(X_arr)
    values = []
    for feature_idx, feature in enumerate(feature_names):
        X_mod = X_arr.copy()
        X_mod[:, feature_idx] = X_mod[:, feature_idx].mean()
        delta = np.abs(baseline - model.predict(X_mod)).mean()
        values.append({'feature': feature, 'mean_abs_contribution': float(delta)})
    return pd.DataFrame(values).sort_values('mean_abs_contribution', ascending=False)


def ml_structure_relationship(df):
    log_event('execute', 'handoff_started', 'ml_relationship')
    ml_df = df[['flexibility_proxy', 'hydrophilicity_index', 'symmetry_index', 'Mn', 'Mw', 'PDI', 'crystallinity_pct', 'contact_angle_deg', 'temperature_C', 'bde_proxy', 'bond_type', 'k_h']].copy()
    ml_df = encode_bond_type(ml_df)
    X = ml_df.drop(columns=['k_h'])
    y = ml_df['k_h'].values
    kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    models = {
        'LinearRegression': LinearRegression(),
        'RandomForest': RandomForestRegressor(n_estimators=350, random_state=SEED, min_samples_leaf=2),
        'GradientBoosting': GradientBoostingRegressor(random_state=SEED),
        'MLP': MLPRegressor(hidden_layer_sizes=(32, 16), max_iter=4000, random_state=SEED),
    }
    metrics_rows = []
    pred_records = []
    for name, model in models.items():
        X_in = X_scaled if name in {'LinearRegression', 'MLP'} else X.values
        pred = cross_val_predict(model, X_in, y, cv=kf)
        rmse = np.sqrt(mean_squared_error(y, pred))
        r2 = r2_score(y, pred)
        ci_low, ci_high = np.quantile(np.abs(pred - y), [0.025, 0.975])
        metrics_rows.append({'model': name, 'rmse': float(rmse), 'r2': float(r2), 'abs_error_ci_low': float(ci_low), 'abs_error_ci_high': float(ci_high)})
        for obs, prd in zip(y, pred):
            pred_records.append({'model': name, 'observed_k_h': float(obs), 'predicted_k_h': float(prd)})
    metrics_df = pd.DataFrame(metrics_rows).sort_values('rmse')
    save_csv(metrics_df, 'ml_model_comparison.csv')
    save_csv(pd.DataFrame(pred_records), 'ml_crossval_predictions.csv')

    best_name = metrics_df.iloc[0]['model']
    best_model = models[best_name]
    best_X = X_scaled if best_name in {'LinearRegression', 'MLP'} else X.values
    best_model.fit(best_X, y)
    imp_df = shap_style_importance(best_model, best_X, list(X.columns))
    save_csv(imp_df, 'ml_shap_style_importance.csv')

    plt.figure(figsize=(6.8, 4.8))
    plt.bar(metrics_df['model'], metrics_df['rmse'], color=['tab:blue', 'tab:orange', 'tab:green', 'tab:purple'])
    plt.ylabel('RMSE (day$^{-1}$)')
    plt.xlabel('Model')
    plt.title('Cross-validated model comparison')
    save_fig('ml_model_comparison.png')

    plt.figure(figsize=(6.8, 4.8))
    top = imp_df.head(12).iloc[::-1]
    plt.barh(top['feature'], top['mean_abs_contribution'], color='tab:green')
    plt.xlabel('Mean absolute contribution')
    plt.ylabel('Descriptor')
    plt.title(f'SHAP-style importance ({best_name})')
    save_fig('ml_shap_style_importance.png')
    log_event('execute', 'handoff_completed', 'ml_relationship', handoff_out={'best_model': best_name}, files_written=[RESULTS_DIR / 'ml_model_comparison.csv', RESULTS_DIR / 'ml_shap_style_importance.csv'])
    return metrics_df, imp_df


def build_statistical_summary(hydro_metrics, ml_metrics, pareto, marine_df, case_summary, top_candidates):
    marine_summary = marine_df.groupby('polymer', as_index=False).agg(final_weight_loss_pct=('weight_loss_pct', 'max'), mean_k_total=('k_total', 'mean'))
    top_case = case_summary.sort_values('final_weight_loss_pct', ascending=False).head(5)
    top_copoly = top_candidates[['system', 'x1', 'x2', 'x3', 'score', 'k_h_copol', 'strength_MPa', 'modulus_GPa']].head(5)
    lines = [
        '# Statistical Summary',
        '',
        '## Hydrolysis model',
    ]
    for _, row in hydro_metrics.iterrows():
        lines.append(f"- {row['metric']}: {row['value']:.5f}")
    lines += [
        '',
        '## ML comparison',
    ]
    for _, row in ml_metrics.iterrows():
        lines.append(f"- {row['model']}: RMSE={row['rmse']:.5f}, R²={row['r2']:.4f}, 95% abs. error interval=[{row['abs_error_ci_low']:.5f}, {row['abs_error_ci_high']:.5f}]")
    lines += [
        '',
        '## Marine degradation summary',
        df_to_markdown(marine_summary),
        '',
        '## Top case-study variants',
        df_to_markdown(top_case),
        '',
        '## Representative Pareto points',
        df_to_markdown(pareto.head(8)),
        '',
        '## Top copolymer candidates',
        df_to_markdown(top_copoly),
        '',
        'Note: Multi-model comparison involves 4 models; rankings were compared descriptively, and confidence intervals for absolute error are reported for transparency.',
    ]
    path = RESULTS_DIR / 'statistical-summary.md'
    path.write_text('\n'.join(lines), encoding='utf-8')
    log_event('report', 'file_written', 'markdown_summary', files_written=[path])
    return marine_summary


def write_run_summary(summary_dict):
    path = RESULTS_DIR / 'summary_metrics.json'
    path.write_text(json.dumps(summary_dict, indent=2, ensure_ascii=False), encoding='utf-8')
    log_event('report', 'file_written', 'json.dump', files_written=[path], handoff_out=summary_dict)
    return path


def main():
    log_event('plan', 'run_started', 'biodegradable_polymer_framework', handoff_in={'seed': SEED})
    log_event('plan', 'prompt_received', 'biodegradable_polymer_framework', handoff_in={'task': 'biodegradable polymer molecular design framework'})
    log_event('plan', 'skill_selected', 'co-scientist-data-analysis', handoff_out={'components': 7})
    dataset = generate_synthetic_dataset(220)
    hydro_model, feature_names, hydro_metrics, hydro_importance = train_hydrolysis_model(dataset)
    pareto = optimize_tradeoff()
    microbial_df = simulate_microbial_degradation()
    marine_df = simulate_marine_environment()
    copoly_df, gp_history, top_candidates = combinatorial_exploration()
    case_df, case_summary = case_study_profiles()
    ml_metrics, ml_importance = ml_structure_relationship(dataset)
    marine_summary = build_statistical_summary(hydro_metrics, ml_metrics, pareto, marine_df, case_summary, top_candidates)
    summary_payload = {
        'dataset_samples': int(len(dataset)),
        'hydrolysis_test_r2': float(hydro_metrics.loc[hydro_metrics['metric'] == 'test_r2', 'value'].iloc[0]),
        'hydrolysis_test_rmse': float(hydro_metrics.loc[hydro_metrics['metric'] == 'test_rmse', 'value'].iloc[0]),
        'pareto_points': int(len(pareto)),
        'marine_final_weight_loss_pct': {row['polymer']: float(row['final_weight_loss_pct']) for _, row in marine_summary.iterrows()},
        'best_ml_model': str(ml_metrics.iloc[0]['model']),
        'best_ml_r2': float(ml_metrics.iloc[0]['r2']),
        'best_copolymer_score': float(top_candidates.iloc[0]['score']),
        'best_copolymer_system': str(top_candidates.iloc[0]['system']),
        'most_degradable_case_variant': str(case_summary.sort_values('final_weight_loss_pct', ascending=False).iloc[0]['variant']),
    }
    write_run_summary(summary_payload)
    log_event('report', 'report_finalized', 'biodegradable_polymer_framework', handoff_out=summary_payload)
    log_event('report', 'run_completed', 'biodegradable_polymer_framework', handoff_out={'status': 'success'})
    print(json.dumps(summary_payload, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
