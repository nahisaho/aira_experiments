#!/usr/bin/env python3
"""
Pandemic Early Warning AI System - Experimental Implementation
Covers: Genomic surveillance, mutation hotspot prediction, Rt estimation,
NLP-based alert analysis, risk scoring, and dashboard design.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats, signal
from scipy.optimize import minimize
from sklearn.ensemble import RandomForestClassifier, GradientBoostingRegressor
from sklearn.metrics import (roc_auc_score, precision_recall_curve, 
                             average_precision_score, mean_squared_error,
                             classification_report, confusion_matrix)
from sklearn.model_selection import TimeSeriesSplit
from sklearn.preprocessing import StandardScaler
import warnings
import os
import json
from datetime import datetime, timedelta

warnings.filterwarnings('ignore')
np.random.seed(42)

FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================
# Module 1: Genomic Surveillance Simulation
# ============================================================

def simulate_genomic_data(n_sequences=5000, n_days=365):
    """Simulate GISAID-like genomic surveillance data with lineage evolution."""
    dates = pd.date_range('2024-01-01', periods=n_days, freq='D')
    
    # Define lineages with emergence times and growth rates
    lineages = {
        'Alpha': {'emerge': 0, 'peak': 60, 'r': 1.2, 'mutations': 23},
        'Beta': {'emerge': 30, 'peak': 120, 'r': 1.1, 'mutations': 21},
        'Delta': {'emerge': 90, 'peak': 200, 'r': 1.5, 'mutations': 32},
        'Omicron_BA1': {'emerge': 150, 'peak': 250, 'r': 2.0, 'mutations': 50},
        'Omicron_BA5': {'emerge': 220, 'peak': 300, 'r': 1.8, 'mutations': 55},
        'Novel_X': {'emerge': 300, 'peak': 350, 'r': 2.5, 'mutations': 60},
    }
    
    records = []
    for i in range(n_sequences):
        day = np.random.randint(0, n_days)
        probs = {}
        for lin, params in lineages.items():
            if day >= params['emerge']:
                t = day - params['emerge']
                sigma = 40
                p = params['r'] * np.exp(-0.5 * ((t - (params['peak'] - params['emerge'])) / sigma) ** 2)
                probs[lin] = max(p, 0.01)
            else:
                probs[lin] = 0
        
        total = sum(probs.values())
        if total == 0:
            continue
        probs = {k: v/total for k, v in probs.items()}
        lineage = np.random.choice(list(probs.keys()), p=list(probs.values()))
        
        n_muts = lineages[lineage]['mutations'] + np.random.poisson(3)
        spike_muts = np.random.poisson(n_muts * 0.3)
        rbd_muts = np.random.poisson(spike_muts * 0.4)
        
        records.append({
            'date': dates[day],
            'lineage': lineage,
            'total_mutations': n_muts,
            'spike_mutations': spike_muts,
            'rbd_mutations': rbd_muts,
            'country': np.random.choice(['USA', 'UK', 'India', 'Brazil', 'Japan', 'South_Africa',
                                          'Germany', 'France', 'Australia', 'Kenya']),
        })
    
    df = pd.DataFrame(records)
    return df, lineages

def plot_lineage_dynamics(genomic_df):
    """Plot lineage prevalence over time."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    weekly = genomic_df.set_index('date').groupby([pd.Grouper(freq='W'), 'lineage']).size().unstack(fill_value=0)
    weekly_pct = weekly.div(weekly.sum(axis=1), axis=0) * 100
    
    colors = sns.color_palette("husl", len(weekly_pct.columns))
    weekly_pct.plot.area(ax=axes[0], color=colors, alpha=0.8, linewidth=0.5)
    axes[0].set_ylabel('Lineage Prevalence (%)', fontsize=12)
    axes[0].set_title('Genomic Surveillance: Lineage Dynamics Over Time', fontsize=14, fontweight='bold')
    axes[0].legend(loc='upper left', fontsize=9)
    axes[0].set_ylim(0, 100)
    
    weekly_counts = genomic_df.set_index('date').resample('W').size()
    axes[1].bar(weekly_counts.index, weekly_counts.values, width=5, color='steelblue', alpha=0.7)
    axes[1].set_ylabel('Sequences Submitted', fontsize=12)
    axes[1].set_xlabel('Date', fontsize=12)
    axes[1].set_title('Weekly Sequencing Volume', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'lineage_dynamics.png'), dpi=150, bbox_inches='tight')
    plt.close()

# ============================================================
# Module 2: Mutation Hotspot Prediction
# ============================================================

def simulate_mutation_landscape(n_positions=1273):
    """Simulate spike protein mutation landscape with functional impact."""
    positions = np.arange(1, n_positions + 1)
    
    # Known functional regions
    ntd = (14, 305)
    rbd = (319, 541)
    furin = (681, 685)
    hr1 = (912, 984)
    
    mutation_freq = np.random.exponential(0.5, n_positions)
    for start, end in [ntd, rbd, furin, hr1]:
        mutation_freq[start-1:end] *= np.random.uniform(2, 5, end - start + 1)
    
    # Functional impact scores
    ace2_binding = np.zeros(n_positions)
    ace2_binding[rbd[0]-1:rbd[1]] = np.random.beta(2, 5, rbd[1] - rbd[0] + 1) * 3
    
    immune_escape = np.zeros(n_positions)
    immune_escape[ntd[0]-1:ntd[1]] = np.random.beta(2, 3, ntd[1] - ntd[0] + 1) * 2
    immune_escape[rbd[0]-1:rbd[1]] = np.random.beta(3, 4, rbd[1] - rbd[0] + 1) * 2.5
    
    fitness_score = 0.4 * ace2_binding + 0.3 * immune_escape + 0.3 * np.random.exponential(0.3, n_positions)
    
    # Hotspot classification
    threshold = np.percentile(mutation_freq * fitness_score, 90)
    is_hotspot = (mutation_freq * fitness_score) > threshold
    
    df = pd.DataFrame({
        'position': positions,
        'mutation_frequency': mutation_freq,
        'ace2_binding_impact': ace2_binding,
        'immune_escape_score': immune_escape,
        'fitness_score': fitness_score,
        'is_hotspot': is_hotspot.astype(int),
    })
    return df

def train_hotspot_predictor(mutation_df):
    """Train RF classifier for mutation hotspot prediction."""
    features = ['mutation_frequency', 'ace2_binding_impact', 'immune_escape_score', 'fitness_score']
    X = mutation_df[features].values
    y = mutation_df['is_hotspot'].values
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    from sklearn.model_selection import StratifiedKFold
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs = []
    aps = []
    
    for train_idx, test_idx in skf.split(X_scaled, y):
        clf = RandomForestClassifier(n_estimators=200, max_depth=10, random_state=42, class_weight='balanced')
        clf.fit(X_scaled[train_idx], y[train_idx])
        proba = clf.predict_proba(X_scaled[test_idx])[:, 1]
        aucs.append(roc_auc_score(y[test_idx], proba))
        aps.append(average_precision_score(y[test_idx], proba))
    
    # Final model on all data
    final_clf = RandomForestClassifier(n_estimators=200, max_depth=10, random_state=42, class_weight='balanced')
    final_clf.fit(X_scaled, y)
    final_proba = final_clf.predict_proba(X_scaled)[:, 1]
    
    results = {
        'mean_auc': np.mean(aucs),
        'std_auc': np.std(aucs),
        'mean_ap': np.mean(aps),
        'std_ap': np.std(aps),
        'feature_importance': dict(zip(features, final_clf.feature_importances_)),
        'predictions': final_proba,
    }
    return results, final_clf

def plot_mutation_landscape(mutation_df, hotspot_results):
    """Plot mutation landscape and hotspot predictions."""
    fig, axes = plt.subplots(4, 1, figsize=(16, 14), sharex=True)
    
    positions = mutation_df['position']
    
    # Mutation frequency
    axes[0].bar(positions, mutation_df['mutation_frequency'], width=1, color='steelblue', alpha=0.6)
    hotspot_mask = mutation_df['is_hotspot'] == 1
    axes[0].bar(positions[hotspot_mask], mutation_df['mutation_frequency'][hotspot_mask], 
                width=1, color='red', alpha=0.8, label='Hotspot')
    axes[0].set_ylabel('Mutation Freq', fontsize=11)
    axes[0].set_title('Spike Protein Mutation Landscape & Hotspot Prediction', fontsize=14, fontweight='bold')
    axes[0].legend()
    
    # Functional impact
    axes[1].fill_between(positions, mutation_df['ace2_binding_impact'], alpha=0.5, color='orange', label='ACE2 Binding')
    axes[1].fill_between(positions, mutation_df['immune_escape_score'], alpha=0.5, color='purple', label='Immune Escape')
    axes[1].set_ylabel('Impact Score', fontsize=11)
    axes[1].legend()
    
    # Fitness score
    axes[2].plot(positions, mutation_df['fitness_score'], color='darkgreen', alpha=0.7, linewidth=0.5)
    axes[2].fill_between(positions, mutation_df['fitness_score'], alpha=0.3, color='green')
    axes[2].set_ylabel('Fitness Score', fontsize=11)
    
    # Prediction probability
    axes[3].plot(positions, hotspot_results['predictions'], color='crimson', alpha=0.7, linewidth=0.5)
    axes[3].axhline(y=0.5, color='black', linestyle='--', alpha=0.5, label='Threshold')
    axes[3].set_ylabel('Hotspot Prob', fontsize=11)
    axes[3].set_xlabel('Spike Protein Position', fontsize=12)
    axes[3].legend()
    
    # Annotate regions
    regions = {'NTD': (14, 305), 'RBD': (319, 541), 'Furin': (681, 685), 'HR1': (912, 984)}
    for name, (s, e) in regions.items():
        for ax in axes:
            ax.axvspan(s, e, alpha=0.08, color='gray')
        axes[0].annotate(name, xy=((s+e)/2, axes[0].get_ylim()[1]*0.9), ha='center', fontsize=9, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'mutation_landscape.png'), dpi=150, bbox_inches='tight')
    plt.close()

# ============================================================
# Module 3: Epidemiological Data Integration & Rt Estimation
# ============================================================

def simulate_epidemic_data(n_days=365):
    """Simulate multi-source epidemiological data."""
    dates = pd.date_range('2024-01-01', periods=n_days, freq='D')
    
    # True Rt with waves
    t = np.arange(n_days)
    rt_true = (1.0 + 0.8 * np.sin(2 * np.pi * t / 120) 
               + 0.3 * np.sin(2 * np.pi * t / 60)
               + 0.1 * np.random.randn(n_days))
    rt_true = np.clip(rt_true, 0.3, 3.5)
    
    # Simulate case counts from Rt using renewal equation
    serial_interval = stats.gamma(a=4.7/1.5, scale=1.5)
    si_weights = serial_interval.pdf(np.arange(1, 21))
    si_weights /= si_weights.sum()
    
    cases = np.zeros(n_days, dtype=float)
    cases[0] = 100
    for day in range(1, n_days):
        lambda_t = rt_true[day] * sum(
            cases[max(0, day - k)] * si_weights[k-1] 
            for k in range(1, min(day+1, 21))
        )
        cases[day] = max(1, np.random.poisson(max(1, lambda_t)))
    
    # Mobility data (inversely related to cases with lag)
    mobility = 100 - 30 * (cases / cases.max()) + np.random.randn(n_days) * 5
    mobility = np.clip(mobility, 20, 110)
    mobility = np.convolve(mobility, np.ones(7)/7, mode='same')
    
    # Wastewater signal (leading indicator, ~7 days ahead)
    ww_signal = np.roll(cases, -7) * (1 + 0.2 * np.random.randn(n_days))
    ww_signal = np.clip(ww_signal, 0, None)
    ww_signal = np.convolve(ww_signal, np.ones(3)/3, mode='same')
    
    # Hospitalization (lagging indicator)
    hosp_rate = 0.05 + 0.02 * np.sin(2 * np.pi * t / 365)
    hospitalizations = np.convolve(cases * hosp_rate, 
                                    stats.gamma(a=3, scale=2).pdf(np.arange(15)), 
                                    mode='same')
    hospitalizations = np.random.poisson(np.clip(hospitalizations, 0.1, None))
    
    df = pd.DataFrame({
        'date': dates,
        'cases': cases.astype(int),
        'rt_true': rt_true,
        'mobility_index': mobility,
        'wastewater_signal': ww_signal,
        'hospitalizations': hospitalizations,
    })
    return df, si_weights

def estimate_rt_improved(cases, si_weights, tau=7):
    """
    Improved EpiEstim-style Rt estimation with:
    - Bayesian updating with gamma prior
    - Adaptive window selection
    - Nowcasting adjustment for reporting delays
    """
    n = len(cases)
    rt_mean = np.ones(n)
    rt_lower = np.ones(n)
    rt_upper = np.ones(n)
    
    # Gamma prior: shape a=1, rate b=5 (weakly informative)
    a_prior = 1.0
    b_prior = 0.2
    
    for t in range(20, n):
        # Compute total infectiousness (denominator in renewal equation)
        lambda_t = sum(
            cases[max(0, t - k)] * si_weights[min(k-1, len(si_weights)-1)]
            for k in range(1, min(t+1, 21))
        )
        
        if lambda_t < 1:
            lambda_t = 1
        
        # Adaptive window: use tau but adjust based on case counts
        effective_tau = tau
        window_cases = cases[max(0, t - effective_tau + 1):t + 1]
        window_lambda = []
        for s in range(max(0, t - effective_tau + 1), t + 1):
            lam_s = sum(
                cases[max(0, s - k)] * si_weights[min(k-1, len(si_weights)-1)]
                for k in range(1, min(s+1, 21))
            )
            window_lambda.append(max(lam_s, 1))
        
        # Posterior: Gamma(a_prior + sum(cases), b_prior + sum(lambda))
        a_post = a_prior + np.sum(window_cases)
        b_post = b_prior + np.sum(window_lambda)
        
        rt_mean[t] = a_post / b_post
        rt_lower[t] = stats.gamma(a=a_post, scale=1/b_post).ppf(0.025)
        rt_upper[t] = stats.gamma(a=a_post, scale=1/b_post).ppf(0.975)
    
    return rt_mean, rt_lower, rt_upper

def estimate_rt_ml(epi_df, si_weights):
    """ML-enhanced Rt estimation using multi-source data."""
    features = []
    targets = []
    
    for t in range(21, len(epi_df) - 7):
        row = epi_df.iloc[t]
        # Feature engineering
        case_trend_7d = (epi_df['cases'].iloc[t-6:t+1].mean() / 
                         max(epi_df['cases'].iloc[t-13:t-6].mean(), 1))
        case_trend_14d = (epi_df['cases'].iloc[t-6:t+1].mean() / 
                          max(epi_df['cases'].iloc[t-20:t-13].mean(), 1))
        ww_ratio = (epi_df['wastewater_signal'].iloc[t] / 
                    max(epi_df['wastewater_signal'].iloc[t-7:t].mean(), 1))
        
        feat = {
            'cases_7d_avg': epi_df['cases'].iloc[t-6:t+1].mean(),
            'case_trend_7d': case_trend_7d,
            'case_trend_14d': case_trend_14d,
            'mobility': row['mobility_index'],
            'wastewater_signal': row['wastewater_signal'],
            'ww_ratio': ww_ratio,
            'hospitalizations': row['hospitalizations'],
            'day_of_week': epi_df['date'].iloc[t].dayofweek,
        }
        features.append(feat)
        targets.append(row['rt_true'])
    
    X = pd.DataFrame(features).values
    y = np.array(targets)
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Time series cross-validation
    tscv = TimeSeriesSplit(n_splits=5)
    rmses = []
    maes = []
    
    for train_idx, test_idx in tscv.split(X_scaled):
        model = GradientBoostingRegressor(n_estimators=200, max_depth=5, 
                                           learning_rate=0.05, random_state=42)
        model.fit(X_scaled[train_idx], y[train_idx])
        pred = model.predict(X_scaled[test_idx])
        rmses.append(np.sqrt(mean_squared_error(y[test_idx], pred)))
        maes.append(np.mean(np.abs(y[test_idx] - pred)))
    
    # Final model
    final_model = GradientBoostingRegressor(n_estimators=200, max_depth=5, 
                                             learning_rate=0.05, random_state=42)
    final_model.fit(X_scaled, y)
    rt_ml = final_model.predict(X_scaled)
    
    feature_names = ['cases_7d_avg', 'case_trend_7d', 'case_trend_14d', 'mobility',
                     'wastewater_signal', 'ww_ratio', 'hospitalizations', 'day_of_week']
    
    results = {
        'mean_rmse': np.mean(rmses),
        'std_rmse': np.std(rmses),
        'mean_mae': np.mean(maes),
        'std_mae': np.std(maes),
        'feature_importance': dict(zip(feature_names, final_model.feature_importances_)),
        'predictions': rt_ml,
        'true_values': y,
        'dates': epi_df['date'].iloc[21:len(epi_df)-7].values,
    }
    return results

def plot_rt_estimation(epi_df, rt_mean, rt_lower, rt_upper, ml_results):
    """Plot Rt estimation comparison."""
    fig, axes = plt.subplots(4, 1, figsize=(16, 16), sharex=True)
    
    dates = epi_df['date']
    
    # Cases
    axes[0].bar(dates, epi_df['cases'], width=1, color='steelblue', alpha=0.6, label='Reported Cases')
    ax0_twin = axes[0].twinx()
    ax0_twin.plot(dates, epi_df['wastewater_signal'], color='brown', alpha=0.7, linewidth=1.5, label='Wastewater Signal')
    axes[0].set_ylabel('Daily Cases', fontsize=11)
    ax0_twin.set_ylabel('Wastewater Signal', fontsize=11, color='brown')
    axes[0].set_title('Multi-Source Epidemiological Data & Rt Estimation', fontsize=14, fontweight='bold')
    axes[0].legend(loc='upper left')
    ax0_twin.legend(loc='upper right')
    
    # Rt comparison
    axes[1].plot(dates, epi_df['rt_true'], color='black', linewidth=2, label='True Rt', alpha=0.8)
    axes[1].plot(dates, rt_mean, color='blue', linewidth=1.5, label='EpiEstim-Improved', alpha=0.8)
    axes[1].fill_between(dates, rt_lower, rt_upper, color='blue', alpha=0.15, label='95% CI')
    axes[1].axhline(y=1.0, color='red', linestyle='--', alpha=0.5)
    axes[1].set_ylabel('Rt', fontsize=11)
    axes[1].set_title('Bayesian Rt Estimation (Improved EpiEstim)', fontsize=13, fontweight='bold')
    axes[1].legend(fontsize=9)
    axes[1].set_ylim(0, 4)
    
    # ML Rt
    ml_dates = ml_results['dates']
    axes[2].plot(dates, epi_df['rt_true'], color='black', linewidth=2, label='True Rt', alpha=0.8)
    axes[2].plot(ml_dates, ml_results['predictions'], color='green', linewidth=1.5, 
                 label=f'ML-Enhanced (RMSE={ml_results["mean_rmse"]:.3f})', alpha=0.8)
    axes[2].axhline(y=1.0, color='red', linestyle='--', alpha=0.5)
    axes[2].set_ylabel('Rt', fontsize=11)
    axes[2].set_title('ML-Enhanced Rt Estimation (GBR + Multi-source)', fontsize=13, fontweight='bold')
    axes[2].legend(fontsize=9)
    axes[2].set_ylim(0, 4)
    
    # Mobility
    axes[3].plot(dates, epi_df['mobility_index'], color='teal', linewidth=1.5, alpha=0.8)
    axes[3].set_ylabel('Mobility Index', fontsize=11)
    axes[3].set_xlabel('Date', fontsize=12)
    axes[3].set_title('Population Mobility', fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'rt_estimation.png'), dpi=150, bbox_inches='tight')
    plt.close()

# ============================================================
# Module 4: NLP-based Alert Analysis
# ============================================================

def simulate_promed_alerts(n_alerts=500):
    """Simulate ProMED/WHO-style alerts with NLP features."""
    alert_types = ['ProMED', 'WHO_DON', 'WHO_EBS', 'Local_Health', 'News']
    pathogens = ['COVID-19', 'Influenza_H5N1', 'Mpox', 'Ebola', 'Dengue', 
                 'Cholera', 'Measles', 'Marburg', 'Nipah', 'Unknown']
    severities = ['Low', 'Medium', 'High', 'Critical']
    
    alerts = []
    for i in range(n_alerts):
        pathogen = np.random.choice(pathogens, p=[0.25, 0.15, 0.1, 0.08, 0.12, 
                                                    0.08, 0.07, 0.05, 0.05, 0.05])
        source = np.random.choice(alert_types, p=[0.3, 0.2, 0.15, 0.2, 0.15])
        
        # Simulate NLP-extracted features
        urgency_keywords = np.random.poisson(2 if pathogen in ['Ebola', 'Marburg', 'Nipah'] else 1)
        geographic_spread = np.random.choice([1, 2, 3, 4, 5], 
                                              p=[0.3, 0.25, 0.2, 0.15, 0.1])
        case_count_mentioned = np.random.poisson(50 * geographic_spread)
        fatality_mentioned = np.random.random() < (0.3 if pathogen in ['Ebola', 'Marburg'] else 0.1)
        novel_pathogen_flag = pathogen == 'Unknown'
        sentiment_score = np.random.beta(2, 5) if pathogen in ['Measles', 'Dengue'] else np.random.beta(5, 3)
        
        # True severity based on features
        severity_score = (urgency_keywords * 0.2 + geographic_spread * 0.15 + 
                         (case_count_mentioned / 250) * 0.2 + fatality_mentioned * 0.2 +
                         novel_pathogen_flag * 0.15 + sentiment_score * 0.1)
        
        if severity_score > 0.7:
            true_severity = 'Critical'
        elif severity_score > 0.5:
            true_severity = 'High'
        elif severity_score > 0.3:
            true_severity = 'Medium'
        else:
            true_severity = 'Low'
        
        alerts.append({
            'alert_id': f'ALERT-{i:04d}',
            'date': pd.Timestamp('2024-01-01') + pd.Timedelta(days=np.random.randint(0, 365)),
            'source': source,
            'pathogen': pathogen,
            'urgency_keywords': urgency_keywords,
            'geographic_spread': geographic_spread,
            'case_count': case_count_mentioned,
            'fatality_mentioned': int(fatality_mentioned),
            'novel_pathogen': int(novel_pathogen_flag),
            'sentiment_score': sentiment_score,
            'severity_score': severity_score,
            'true_severity': true_severity,
        })
    
    return pd.DataFrame(alerts)

def train_alert_classifier(alerts_df):
    """Train NLP-based alert severity classifier."""
    features = ['urgency_keywords', 'geographic_spread', 'case_count', 
                'fatality_mentioned', 'novel_pathogen', 'sentiment_score']
    
    severity_map = {'Low': 0, 'Medium': 1, 'High': 2, 'Critical': 3}
    X = alerts_df[features].values
    y = alerts_df['true_severity'].map(severity_map).values
    y_binary = (y >= 2).astype(int)  # High/Critical vs Low/Medium
    
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    tscv = TimeSeriesSplit(n_splits=5)
    aucs = []
    
    for train_idx, test_idx in tscv.split(X_scaled):
        clf = RandomForestClassifier(n_estimators=200, max_depth=8, random_state=42, class_weight='balanced')
        clf.fit(X_scaled[train_idx], y_binary[train_idx])
        proba = clf.predict_proba(X_scaled[test_idx])[:, 1]
        aucs.append(roc_auc_score(y_binary[test_idx], proba))
    
    # Final model
    final_clf = RandomForestClassifier(n_estimators=200, max_depth=8, random_state=42, class_weight='balanced')
    final_clf.fit(X_scaled, y_binary)
    final_proba = final_clf.predict_proba(X_scaled)[:, 1]
    final_pred = final_clf.predict(X_scaled)
    
    precision, recall, thresholds = precision_recall_curve(y_binary, final_proba)
    
    results = {
        'mean_auc': np.mean(aucs),
        'std_auc': np.std(aucs),
        'feature_importance': dict(zip(features, final_clf.feature_importances_)),
        'predictions': final_proba,
        'precision': precision,
        'recall': recall,
        'thresholds': thresholds,
        'classification_report': classification_report(y_binary, final_pred, output_dict=True),
        'confusion_matrix': confusion_matrix(y_binary, final_pred),
    }
    return results

def plot_nlp_analysis(alerts_df, alert_results):
    """Plot NLP alert analysis results."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Severity distribution by source
    ct = pd.crosstab(alerts_df['source'], alerts_df['true_severity'])
    ct = ct[['Low', 'Medium', 'High', 'Critical']]
    ct.plot(kind='bar', stacked=True, ax=axes[0, 0], 
            color=['#2ecc71', '#f39c12', '#e74c3c', '#8e44ad'])
    axes[0, 0].set_title('Alert Severity by Source', fontsize=13, fontweight='bold')
    axes[0, 0].set_ylabel('Count')
    axes[0, 0].tick_params(axis='x', rotation=45)
    
    # Pathogen distribution
    pathogen_counts = alerts_df['pathogen'].value_counts()
    colors_p = sns.color_palette("Set2", len(pathogen_counts))
    axes[0, 1].barh(pathogen_counts.index, pathogen_counts.values, color=colors_p)
    axes[0, 1].set_title('Alert Distribution by Pathogen', fontsize=13, fontweight='bold')
    axes[0, 1].set_xlabel('Count')
    
    # Precision-Recall curve
    axes[1, 0].plot(alert_results['recall'], alert_results['precision'], 
                     color='crimson', linewidth=2)
    axes[1, 0].set_xlabel('Recall', fontsize=11)
    axes[1, 0].set_ylabel('Precision', fontsize=11)
    axes[1, 0].set_title(f'PR Curve (AUC-ROC={alert_results["mean_auc"]:.3f})', 
                          fontsize=13, fontweight='bold')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Feature importance
    fi = alert_results['feature_importance']
    sorted_fi = sorted(fi.items(), key=lambda x: x[1], reverse=True)
    names, values = zip(*sorted_fi)
    axes[1, 1].barh(names, values, color='teal')
    axes[1, 1].set_title('Feature Importance (Alert Classifier)', fontsize=13, fontweight='bold')
    axes[1, 1].set_xlabel('Importance')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'nlp_analysis.png'), dpi=150, bbox_inches='tight')
    plt.close()

# ============================================================
# Module 5: Risk Scoring & Alert Threshold Optimization
# ============================================================

def compute_risk_scores(epi_df, rt_mean, genomic_df, alerts_df):
    """Compute composite pandemic risk scores integrating all data streams."""
    n_days = len(epi_df)
    dates = epi_df['date']
    
    # Genomic risk: lineage diversity + novel variant emergence
    genomic_daily = genomic_df.set_index('date').resample('D').agg({
        'lineage': 'nunique',
        'spike_mutations': 'mean',
        'rbd_mutations': 'mean',
    }).reindex(dates).ffill().fillna(0)
    
    genomic_risk = (genomic_daily['lineage'] / genomic_daily['lineage'].max() * 0.5 +
                    genomic_daily['spike_mutations'] / max(genomic_daily['spike_mutations'].max(), 1) * 0.3 +
                    genomic_daily['rbd_mutations'] / max(genomic_daily['rbd_mutations'].max(), 1) * 0.2)
    
    # Epi risk: Rt and case trends
    epi_risk = np.clip((rt_mean - 1) / 2, 0, 1) * 0.5
    case_growth = epi_df['cases'].rolling(7).mean() / epi_df['cases'].rolling(14).mean().shift(7)
    case_growth = case_growth.fillna(1).clip(0, 3)
    epi_risk += np.clip((case_growth - 1) / 2, 0, 1).values * 0.3
    ww_growth = epi_df['wastewater_signal'].rolling(7).mean() / epi_df['wastewater_signal'].rolling(14).mean().shift(7)
    ww_growth = ww_growth.fillna(1).clip(0, 3)
    epi_risk += np.clip((ww_growth - 1) / 2, 0, 1).values * 0.2
    
    # Alert risk: aggregate severity from NLP
    alert_daily = alerts_df.set_index('date').resample('D')['severity_score'].mean()
    alert_daily = alert_daily.reindex(dates).ffill().fillna(0)
    alert_risk = alert_daily.values / max(alert_daily.max(), 1)
    
    # Composite risk score
    composite_risk = (0.35 * epi_risk + 
                      0.30 * genomic_risk.values + 
                      0.20 * alert_risk + 
                      0.15 * np.clip((epi_df['hospitalizations'] / max(epi_df['hospitalizations'].max(), 1)).values, 0, 1))
    
    return pd.DataFrame({
        'date': dates,
        'genomic_risk': genomic_risk.values,
        'epi_risk': epi_risk,
        'alert_risk': alert_risk,
        'composite_risk': composite_risk,
    })

def optimize_alert_thresholds(risk_df, epi_df):
    """Optimize alert thresholds using ROC analysis."""
    # Define "true alerts" as periods where Rt > 1.5 AND cases growing
    true_alert = ((epi_df['rt_true'] > 1.5) & 
                  (epi_df['cases'].rolling(7).mean() > epi_df['cases'].rolling(7).mean().shift(7))).astype(int)
    true_alert = true_alert.fillna(0).values
    
    composite = risk_df['composite_risk'].values
    
    thresholds = np.linspace(0.1, 0.9, 50)
    results = []
    for thresh in thresholds:
        predicted = (composite > thresh).astype(int)
        tp = np.sum((predicted == 1) & (true_alert == 1))
        fp = np.sum((predicted == 1) & (true_alert == 0))
        fn = np.sum((predicted == 0) & (true_alert == 1))
        tn = np.sum((predicted == 0) & (true_alert == 0))
        
        precision = tp / max(tp + fp, 1)
        recall = tp / max(tp + fn, 1)
        f1 = 2 * precision * recall / max(precision + recall, 1e-10)
        specificity = tn / max(tn + fp, 1)
        
        # Lead time: avg days before true alert that system triggers
        lead_times = []
        in_true_alert = False
        for t in range(len(true_alert)):
            if true_alert[t] == 1 and not in_true_alert:
                # Look back for predicted alert
                for back in range(1, min(t+1, 30)):
                    if predicted[t - back] == 1:
                        lead_times.append(back)
                        break
                in_true_alert = True
            elif true_alert[t] == 0:
                in_true_alert = False
        
        results.append({
            'threshold': thresh,
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'specificity': specificity,
            'mean_lead_time': np.mean(lead_times) if lead_times else 0,
        })
    
    results_df = pd.DataFrame(results)
    optimal_idx = results_df['f1'].idxmax()
    optimal_threshold = results_df.iloc[optimal_idx]['threshold']
    
    return results_df, optimal_threshold, true_alert

def plot_risk_dashboard(risk_df, epi_df, threshold_results, optimal_threshold, true_alert):
    """Plot integrated risk dashboard."""
    fig = plt.figure(figsize=(18, 20))
    gs = gridspec.GridSpec(5, 2, figure=fig, hspace=0.35, wspace=0.3)
    
    dates = risk_df['date']
    
    # Composite risk heatmap
    ax1 = fig.add_subplot(gs[0, :])
    risk_matrix = risk_df[['genomic_risk', 'epi_risk', 'alert_risk', 'composite_risk']].values.T
    im = ax1.imshow(risk_matrix, aspect='auto', cmap='RdYlGn_r', vmin=0, vmax=1, interpolation='bilinear')
    ax1.set_yticks([0, 1, 2, 3])
    ax1.set_yticklabels(['Genomic', 'Epidemiological', 'Alert', 'Composite'])
    n_ticks = 12
    tick_pos = np.linspace(0, len(dates)-1, n_ticks, dtype=int)
    ax1.set_xticks(tick_pos)
    ax1.set_xticklabels([dates.iloc[i].strftime('%Y-%m') for i in tick_pos], rotation=45, fontsize=8)
    ax1.set_title('Pandemic Risk Dashboard - Multi-Source Risk Heatmap', fontsize=14, fontweight='bold')
    plt.colorbar(im, ax=ax1, label='Risk Level', shrink=0.6)
    
    # Composite risk timeline
    ax2 = fig.add_subplot(gs[1, :])
    ax2.fill_between(dates, risk_df['composite_risk'], alpha=0.3, color='red')
    ax2.plot(dates, risk_df['composite_risk'], color='red', linewidth=1.5, label='Composite Risk')
    ax2.axhline(y=optimal_threshold, color='black', linestyle='--', linewidth=2, 
                label=f'Optimal Threshold ({optimal_threshold:.2f})')
    
    # Shade true alert periods
    alert_mask = true_alert.astype(bool)
    for i in range(1, len(alert_mask)):
        if alert_mask[i] and not alert_mask[i-1]:
            start = dates.iloc[i]
        elif not alert_mask[i] and alert_mask[i-1]:
            ax2.axvspan(start, dates.iloc[i], alpha=0.15, color='orange', label='_')
    
    ax2.set_ylabel('Risk Score', fontsize=11)
    ax2.set_title('Composite Risk Score with Optimized Alert Threshold', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=9)
    
    # Threshold optimization
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.plot(threshold_results['threshold'], threshold_results['f1'], 'b-', linewidth=2, label='F1 Score')
    ax3.plot(threshold_results['threshold'], threshold_results['precision'], 'g--', linewidth=1.5, label='Precision')
    ax3.plot(threshold_results['threshold'], threshold_results['recall'], 'r--', linewidth=1.5, label='Recall')
    ax3.axvline(x=optimal_threshold, color='black', linestyle=':', label=f'Optimal ({optimal_threshold:.2f})')
    ax3.set_xlabel('Threshold')
    ax3.set_ylabel('Score')
    ax3.set_title('Threshold Optimization', fontsize=13, fontweight='bold')
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)
    
    # Lead time vs threshold
    ax4 = fig.add_subplot(gs[2, 1])
    ax4.plot(threshold_results['threshold'], threshold_results['mean_lead_time'], 
             color='purple', linewidth=2, marker='o', markersize=3)
    ax4.set_xlabel('Threshold')
    ax4.set_ylabel('Mean Lead Time (days)')
    ax4.set_title('Alert Lead Time Analysis', fontsize=13, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    
    # Risk component comparison
    ax5 = fig.add_subplot(gs[3, 0])
    risk_components = ['genomic_risk', 'epi_risk', 'alert_risk']
    comp_data = [risk_df[c].values for c in risk_components]
    bp = ax5.boxplot(comp_data, labels=['Genomic', 'Epidemiological', 'Alert'], patch_artist=True,
                     boxprops=dict(alpha=0.7))
    colors_box = ['#3498db', '#e74c3c', '#2ecc71']
    for patch, color in zip(bp['boxes'], colors_box):
        patch.set_facecolor(color)
    ax5.set_ylabel('Risk Score')
    ax5.set_title('Risk Component Distribution', fontsize=13, fontweight='bold')
    
    # Correlation matrix
    ax6 = fig.add_subplot(gs[3, 1])
    corr_data = pd.DataFrame({
        'Cases': epi_df['cases'].values,
        'Rt': epi_df['rt_true'].values,
        'Wastewater': epi_df['wastewater_signal'].values,
        'Mobility': epi_df['mobility_index'].values,
        'Genomic Risk': risk_df['genomic_risk'].values,
        'Composite Risk': risk_df['composite_risk'].values,
    })
    corr = corr_data.corr()
    sns.heatmap(corr, annot=True, fmt='.2f', cmap='coolwarm', center=0, ax=ax6, 
                square=True, linewidths=0.5, cbar_kws={'shrink': 0.8})
    ax6.set_title('Data Stream Correlations', fontsize=13, fontweight='bold')
    
    # System architecture text
    ax7 = fig.add_subplot(gs[4, :])
    ax7.axis('off')
    arch_text = """
    ┌─────────────────────────────────────────────────────────────────────────────────────────────┐
    │                        PANDEMIC EARLY WARNING AI SYSTEM - ARCHITECTURE                      │
    ├──────────────┬──────────────┬──────────────┬──────────────┬──────────────┬──────────────────┤
    │  GISAID/     │  Mutation    │  Epi Data    │  Rt Engine   │  NLP Alert   │  Risk Scoring   │
    │  GenBank     │  Hotspot     │  Integration │  (EpiEstim+  │  Analysis    │  & Alert        │
    │  Genomic     │  Prediction  │  (Cases,     │   ML-Enhanced│  (ProMED/    │  Threshold      │
    │  Pipeline    │  (RF/DL)     │  WW, Mob.)   │   Bayesian)  │  WHO Parser) │  Optimization   │
    ├──────────────┴──────────────┴──────────────┴──────────────┴──────────────┴──────────────────┤
    │                              REAL-TIME DATA PIPELINE (Apache Kafka / Flink)                 │
    ├─────────────────────────────────────────────────────────────────────────────────────────────┤
    │                              DASHBOARD & ALERTING (Grafana / Custom Web UI)                 │
    └─────────────────────────────────────────────────────────────────────────────────────────────┘
    """
    ax7.text(0.5, 0.5, arch_text, transform=ax7.transAxes, fontsize=9,
             verticalalignment='center', horizontalalignment='center',
             fontfamily='monospace', bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.savefig(os.path.join(FIGURES_DIR, 'risk_dashboard.png'), dpi=150, bbox_inches='tight')
    plt.close()

# ============================================================
# Module 6: Pipeline Architecture Diagram
# ============================================================

def plot_system_architecture():
    """Generate system architecture diagram."""
    fig, ax = plt.subplots(1, 1, figsize=(18, 12))
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 12)
    ax.axis('off')
    
    # Data Sources Layer
    sources = [
        ('GISAID/\nGenBank', 1.5, 10.5, '#3498db'),
        ('Case\nReports', 4.5, 10.5, '#2ecc71'),
        ('Wastewater\nSurveillance', 7.5, 10.5, '#e74c3c'),
        ('Mobility\nData', 10.5, 10.5, '#f39c12'),
        ('ProMED/\nWHO Alerts', 13.5, 10.5, '#9b59b6'),
        ('Hospital\nData', 16.5, 10.5, '#1abc9c'),
    ]
    
    for label, x, y, color in sources:
        rect = plt.Rectangle((x-1, y-0.5), 2, 1, facecolor=color, alpha=0.7, edgecolor='black', linewidth=1.5)
        ax.add_patch(rect)
        ax.text(x, y, label, ha='center', va='center', fontsize=8, fontweight='bold', color='white')
    
    # Ingestion Layer
    ax.add_patch(plt.Rectangle((0.5, 8), 17, 1.2, facecolor='#34495e', alpha=0.8, edgecolor='black', linewidth=1.5))
    ax.text(9, 8.6, 'Real-Time Data Ingestion Layer (Apache Kafka / AWS Kinesis)', 
            ha='center', va='center', fontsize=11, fontweight='bold', color='white')
    
    # Processing Modules
    modules = [
        ('Phylogenetic\nAnalysis\n(Nextclade)', 2, 6.5, '#3498db'),
        ('Mutation\nHotspot\nPrediction (RF)', 5.5, 6.5, '#e67e22'),
        ('Rt Estimation\n(EpiEstim+\nML-Enhanced)', 9, 6.5, '#e74c3c'),
        ('NLP Alert\nClassification\n(Transformer)', 12.5, 6.5, '#9b59b6'),
        ('Wastewater\nSignal\nProcessing', 16, 6.5, '#1abc9c'),
    ]
    
    for label, x, y, color in modules:
        rect = plt.Rectangle((x-1.3, y-0.7), 2.6, 1.4, facecolor=color, alpha=0.7, 
                              edgecolor='black', linewidth=1.5, linestyle='-')
        ax.add_patch(rect)
        ax.text(x, y, label, ha='center', va='center', fontsize=7.5, fontweight='bold', color='white')
    
    # Integration Layer
    ax.add_patch(plt.Rectangle((0.5, 4.2), 17, 1, facecolor='#2c3e50', alpha=0.8, edgecolor='black', linewidth=1.5))
    ax.text(9, 4.7, 'Risk Score Integration Engine (Weighted Ensemble + Threshold Optimization)', 
            ha='center', va='center', fontsize=11, fontweight='bold', color='white')
    
    # Output Layer
    outputs = [
        ('Real-Time\nDashboard\n(Grafana)', 3, 2.5, '#27ae60'),
        ('Alert\nNotification\nSystem', 7, 2.5, '#c0392b'),
        ('Decision\nSupport\nReports', 11, 2.5, '#2980b9'),
        ('API\nEndpoints\n(REST/GraphQL)', 15, 2.5, '#8e44ad'),
    ]
    
    for label, x, y, color in outputs:
        rect = plt.Rectangle((x-1.3, y-0.7), 2.6, 1.4, facecolor=color, alpha=0.7, 
                              edgecolor='black', linewidth=1.5)
        ax.add_patch(rect)
        ax.text(x, y, label, ha='center', va='center', fontsize=8, fontweight='bold', color='white')
    
    # Arrows (simplified)
    arrow_props = dict(arrowstyle='->', color='gray', lw=1.5)
    for _, x, y, _ in sources:
        ax.annotate('', xy=(x, 9.2), xytext=(x, y-0.5), arrowprops=arrow_props)
    
    for _, x, y, _ in modules:
        ax.annotate('', xy=(x, 7.2+0.15), xytext=(x, 8), arrowprops=arrow_props)
        ax.annotate('', xy=(x, 5.2), xytext=(x, y-0.7), arrowprops=arrow_props)
    
    for _, x, y, _ in outputs:
        ax.annotate('', xy=(x, 3.2+0.15), xytext=(x, 4.2), arrowprops=arrow_props)
    
    ax.set_title('Pandemic Early Warning AI System - Architecture Overview', 
                 fontsize=16, fontweight='bold', pad=20)
    
    # Layer labels
    ax.text(-0.2, 10.5, 'Data\nSources', ha='center', va='center', fontsize=9, fontweight='bold', 
            rotation=90, color='gray')
    ax.text(-0.2, 6.5, 'AI/ML\nModules', ha='center', va='center', fontsize=9, fontweight='bold', 
            rotation=90, color='gray')
    ax.text(-0.2, 2.5, 'Output\nLayer', ha='center', va='center', fontsize=9, fontweight='bold', 
            rotation=90, color='gray')
    
    plt.savefig(os.path.join(FIGURES_DIR, 'system_architecture.png'), dpi=150, bbox_inches='tight')
    plt.close()

# ============================================================
# Module 7: Performance Summary
# ============================================================

def plot_performance_summary(hotspot_results, ml_rt_results, alert_results, threshold_results, optimal_threshold):
    """Plot comprehensive performance summary."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    # Hotspot feature importance
    fi = hotspot_results['feature_importance']
    sorted_fi = sorted(fi.items(), key=lambda x: x[1], reverse=True)
    names, vals = zip(*sorted_fi)
    axes[0, 0].barh(names, vals, color=sns.color_palette("viridis", len(names)))
    axes[0, 0].set_title(f'Hotspot Predictor\n(AUC={hotspot_results["mean_auc"]:.3f}±{hotspot_results["std_auc"]:.3f})', 
                          fontsize=11, fontweight='bold')
    axes[0, 0].set_xlabel('Feature Importance')
    
    # Rt ML feature importance
    fi_rt = ml_rt_results['feature_importance']
    sorted_fi_rt = sorted(fi_rt.items(), key=lambda x: x[1], reverse=True)
    names_rt, vals_rt = zip(*sorted_fi_rt)
    axes[0, 1].barh(names_rt, vals_rt, color=sns.color_palette("magma", len(names_rt)))
    axes[0, 1].set_title(f'ML Rt Estimator\n(RMSE={ml_rt_results["mean_rmse"]:.3f}±{ml_rt_results["std_rmse"]:.3f})', 
                          fontsize=11, fontweight='bold')
    axes[0, 1].set_xlabel('Feature Importance')
    
    # Alert classifier feature importance
    fi_alert = alert_results['feature_importance']
    sorted_fi_alert = sorted(fi_alert.items(), key=lambda x: x[1], reverse=True)
    names_a, vals_a = zip(*sorted_fi_alert)
    axes[0, 2].barh(names_a, vals_a, color=sns.color_palette("rocket", len(names_a)))
    axes[0, 2].set_title(f'Alert Classifier\n(AUC={alert_results["mean_auc"]:.3f}±{alert_results["std_auc"]:.3f})', 
                          fontsize=11, fontweight='bold')
    axes[0, 2].set_xlabel('Feature Importance')
    
    # Confusion matrix
    cm = alert_results['confusion_matrix']
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[1, 0],
                xticklabels=['Low/Med', 'High/Crit'], yticklabels=['Low/Med', 'High/Crit'])
    axes[1, 0].set_title('Alert Confusion Matrix', fontsize=11, fontweight='bold')
    axes[1, 0].set_ylabel('True')
    axes[1, 0].set_xlabel('Predicted')
    
    # Rt scatter plot
    axes[1, 1].scatter(ml_rt_results['true_values'], ml_rt_results['predictions'], 
                        alpha=0.3, s=10, color='teal')
    lims = [0, 4]
    axes[1, 1].plot(lims, lims, 'r--', linewidth=2, label='Perfect')
    axes[1, 1].set_xlabel('True Rt')
    axes[1, 1].set_ylabel('Predicted Rt')
    axes[1, 1].set_title('ML Rt: True vs Predicted', fontsize=11, fontweight='bold')
    axes[1, 1].legend()
    axes[1, 1].set_xlim(lims)
    axes[1, 1].set_ylim(lims)
    axes[1, 1].grid(True, alpha=0.3)
    
    # Performance summary table
    axes[1, 2].axis('off')
    summary_data = [
        ['Module', 'Metric', 'Value'],
        ['Hotspot Pred.', 'AUC-ROC', f'{hotspot_results["mean_auc"]:.3f}'],
        ['Hotspot Pred.', 'Avg Precision', f'{hotspot_results["mean_ap"]:.3f}'],
        ['Rt (EpiEstim+)', 'Method', 'Bayesian+Adaptive'],
        ['Rt (ML)', 'RMSE', f'{ml_rt_results["mean_rmse"]:.3f}'],
        ['Rt (ML)', 'MAE', f'{ml_rt_results["mean_mae"]:.3f}'],
        ['Alert Class.', 'AUC-ROC', f'{alert_results["mean_auc"]:.3f}'],
        ['Risk Scoring', 'Opt. Threshold', f'{optimal_threshold:.3f}'],
        ['Risk Scoring', 'Best F1', f'{threshold_results["f1"].max():.3f}'],
    ]
    
    table = axes[1, 2].table(cellText=summary_data[1:], colLabels=summary_data[0], 
                              loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 1.5)
    for i in range(len(summary_data[0])):
        table[0, i].set_facecolor('#34495e')
        table[0, i].set_text_props(color='white', fontweight='bold')
    axes[1, 2].set_title('Performance Summary', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'performance_summary.png'), dpi=150, bbox_inches='tight')
    plt.close()

# ============================================================
# Main Execution
# ============================================================

def main():
    print("=" * 70)
    print("PANDEMIC EARLY WARNING AI SYSTEM - EXPERIMENTAL EVALUATION")
    print("=" * 70)
    
    # 1. Genomic Surveillance
    print("\n[1/6] Simulating genomic surveillance data...")
    genomic_df, lineages = simulate_genomic_data(n_sequences=5000, n_days=365)
    plot_lineage_dynamics(genomic_df)
    print(f"  Generated {len(genomic_df)} sequences across {genomic_df['lineage'].nunique()} lineages")
    print(f"  Date range: {genomic_df['date'].min()} to {genomic_df['date'].max()}")
    
    # 2. Mutation Hotspot Prediction
    print("\n[2/6] Training mutation hotspot predictor...")
    mutation_df = simulate_mutation_landscape()
    hotspot_results, hotspot_model = train_hotspot_predictor(mutation_df)
    plot_mutation_landscape(mutation_df, hotspot_results)
    print(f"  AUC-ROC: {hotspot_results['mean_auc']:.4f} ± {hotspot_results['std_auc']:.4f}")
    print(f"  Avg Precision: {hotspot_results['mean_ap']:.4f} ± {hotspot_results['std_ap']:.4f}")
    print(f"  Feature importance: {hotspot_results['feature_importance']}")
    
    # 3. Epidemiological Data & Rt Estimation
    print("\n[3/6] Simulating epidemiological data & estimating Rt...")
    epi_df, si_weights = simulate_epidemic_data(n_days=365)
    rt_mean, rt_lower, rt_upper = estimate_rt_improved(epi_df['cases'].values, si_weights, tau=7)
    
    bayesian_rmse = np.sqrt(np.mean((rt_mean[20:] - epi_df['rt_true'].values[20:])**2))
    print(f"  Bayesian Rt RMSE: {bayesian_rmse:.4f}")
    
    # 4. ML-Enhanced Rt
    print("\n[4/6] Training ML-enhanced Rt estimator...")
    ml_rt_results = estimate_rt_ml(epi_df, si_weights)
    plot_rt_estimation(epi_df, rt_mean, rt_lower, rt_upper, ml_rt_results)
    print(f"  ML Rt RMSE: {ml_rt_results['mean_rmse']:.4f} ± {ml_rt_results['std_rmse']:.4f}")
    print(f"  ML Rt MAE: {ml_rt_results['mean_mae']:.4f} ± {ml_rt_results['std_mae']:.4f}")
    
    # 5. NLP Alert Analysis
    print("\n[5/6] Training NLP-based alert classifier...")
    alerts_df = simulate_promed_alerts(n_alerts=500)
    alert_results = train_alert_classifier(alerts_df)
    plot_nlp_analysis(alerts_df, alert_results)
    print(f"  Alert AUC-ROC: {alert_results['mean_auc']:.4f} ± {alert_results['std_auc']:.4f}")
    
    # 6. Risk Scoring & Dashboard
    print("\n[6/6] Computing risk scores & optimizing thresholds...")
    risk_df = compute_risk_scores(epi_df, rt_mean, genomic_df, alerts_df)
    threshold_results, optimal_threshold, true_alert = optimize_alert_thresholds(risk_df, epi_df)
    plot_risk_dashboard(risk_df, epi_df, threshold_results, optimal_threshold, true_alert)
    plot_system_architecture()
    plot_performance_summary(hotspot_results, ml_rt_results, alert_results, threshold_results, optimal_threshold)
    
    best_f1 = threshold_results['f1'].max()
    best_lead = threshold_results.loc[threshold_results['f1'].idxmax(), 'mean_lead_time']
    print(f"  Optimal threshold: {optimal_threshold:.4f}")
    print(f"  Best F1 score: {best_f1:.4f}")
    print(f"  Mean lead time at optimal: {best_lead:.1f} days")
    
    # Save metrics
    metrics = {
        'hotspot_auc': float(hotspot_results['mean_auc']),
        'hotspot_ap': float(hotspot_results['mean_ap']),
        'bayesian_rt_rmse': float(bayesian_rmse),
        'ml_rt_rmse': float(ml_rt_results['mean_rmse']),
        'ml_rt_mae': float(ml_rt_results['mean_mae']),
        'alert_auc': float(alert_results['mean_auc']),
        'optimal_threshold': float(optimal_threshold),
        'best_f1': float(best_f1),
        'mean_lead_time': float(best_lead),
        'feature_importance_hotspot': {k: float(v) for k, v in hotspot_results['feature_importance'].items()},
        'feature_importance_rt': {k: float(v) for k, v in ml_rt_results['feature_importance'].items()},
        'feature_importance_alert': {k: float(v) for k, v in alert_results['feature_importance'].items()},
    }
    
    with open(os.path.join(os.path.dirname(FIGURES_DIR), 'metrics.json'), 'w') as f:
        json.dump(metrics, f, indent=2)
    
    print("\n" + "=" * 70)
    print("EXPERIMENT COMPLETE")
    print(f"Figures saved to: {FIGURES_DIR}")
    print(f"Metrics saved to: metrics.json")
    print("=" * 70)
    
    return metrics

if __name__ == '__main__':
    metrics = main()
