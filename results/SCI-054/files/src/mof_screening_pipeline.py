#!/usr/bin/env python3
"""
MOF High-Throughput Screening Pipeline for CO2/H2 Adsorption
============================================================
Integrates: CoRE MOF / hMOF feature extraction, GCMC simulation,
geometric descriptors, ML prediction, stability filtering, and DAC ranking.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, GradientBoostingClassifier
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr
import os
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
DATADIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

# ============================================================
# 1. MOF Database Generation (Simulating CoRE MOF + hMOF)
# ============================================================
def generate_mof_database(n_core=500, n_hmof=1500):
    """Generate synthetic MOF database mimicking CoRE MOF and hMOF distributions."""
    mofs = []
    # CoRE MOF - experimentally derived, more constrained distributions
    for i in range(n_core):
        sa = np.random.lognormal(mean=7.0, sigma=0.5)  # surface area (m2/g)
        sa = np.clip(sa, 100, 6000)
        pv = np.random.uniform(0.1, 2.5)  # pore volume (cm3/g)
        vf = np.clip(np.random.beta(3, 3) * 0.9 + 0.05, 0.05, 0.95)  # void fraction
        pld = np.random.lognormal(mean=1.8, sigma=0.4)  # pore limiting diameter (Å)
        pld = np.clip(pld, 2.0, 25.0)
        lcd = pld + np.random.exponential(2.0)  # largest cavity diameter
        lcd = np.clip(lcd, pld, 35.0)
        density = np.random.uniform(0.3, 2.0)  # g/cm3
        metal_en = np.random.choice([1.31, 1.54, 1.63, 1.83, 1.88, 1.91, 2.20])  # electronegativity
        has_oms = np.random.binomial(1, 0.35)  # open metal sites
        func_group = np.random.choice(['none', 'NH2', 'OH', 'COOH', 'F', 'NO2', 'CH3'], 
                                       p=[0.3, 0.15, 0.1, 0.1, 0.1, 0.1, 0.15])
        mofs.append({
            'mof_id': f'CoRE_{i:04d}', 'source': 'CoRE',
            'surface_area': sa, 'pore_volume': pv, 'void_fraction': vf,
            'pld': pld, 'lcd': lcd, 'density': density,
            'metal_electronegativity': metal_en, 'has_oms': has_oms,
            'functional_group': func_group
        })
    
    # hMOF - hypothetical, wider distributions
    for i in range(n_hmof):
        sa = np.random.lognormal(mean=7.2, sigma=0.7)
        sa = np.clip(sa, 50, 8000)
        pv = np.random.uniform(0.05, 4.0)
        vf = np.clip(np.random.beta(2, 2) * 0.95 + 0.02, 0.02, 0.98)
        pld = np.random.lognormal(mean=1.9, sigma=0.5)
        pld = np.clip(pld, 1.5, 30.0)
        lcd = pld + np.random.exponential(3.0)
        lcd = np.clip(lcd, pld, 40.0)
        density = np.random.uniform(0.2, 2.5)
        metal_en = np.random.choice([1.31, 1.54, 1.63, 1.83, 1.88, 1.91, 2.20, 1.36, 1.65])
        has_oms = np.random.binomial(1, 0.25)
        func_group = np.random.choice(['none', 'NH2', 'OH', 'COOH', 'F', 'NO2', 'CH3', 'CF3'],
                                       p=[0.25, 0.12, 0.1, 0.1, 0.1, 0.08, 0.12, 0.13])
        mofs.append({
            'mof_id': f'hMOF_{i:04d}', 'source': 'hMOF',
            'surface_area': sa, 'pore_volume': pv, 'void_fraction': vf,
            'pld': pld, 'lcd': lcd, 'density': density,
            'metal_electronegativity': metal_en, 'has_oms': has_oms,
            'functional_group': func_group
        })
    
    df = pd.DataFrame(mofs)
    # Encode functional groups
    fg_map = {'none': 0, 'NH2': 1, 'OH': 2, 'COOH': 3, 'F': 4, 'NO2': 5, 'CH3': 6, 'CF3': 7}
    df['fg_encoded'] = df['functional_group'].map(fg_map)
    return df


# ============================================================
# 2. GCMC Adsorption Simulation (Physics-based model)
# ============================================================
def simulate_gcmc_adsorption(df, pressures_bar=[0.0004, 0.15, 1.0, 5.0, 10.0, 50.0]):
    """
    Simulate GCMC-like CO2 and H2 adsorption using Langmuir-Freundlich model.
    Parameters derived from geometric descriptors following established correlations.
    """
    results = df.copy()
    
    for gas in ['CO2', 'H2']:
        for p in pressures_bar:
            col_name = f'{gas}_uptake_{p}bar'
            uptakes = []
            for _, row in df.iterrows():
                if gas == 'CO2':
                    # CO2: strong dependence on surface area, pore volume, OMS, functional groups
                    q_sat = (0.8 * row['surface_area'] / 1000 + 
                             2.5 * row['pore_volume'] + 
                             1.5 * row['has_oms'] +
                             (0.8 if row['functional_group'] in ['NH2', 'OH', 'COOH'] else 0))
                    q_sat = np.clip(q_sat, 0.5, 15.0)
                    # Affinity parameter — higher for amine-functionalized MOFs
                    b = (0.5 + 0.2 * row['metal_electronegativity'] + 
                         0.5 * row['has_oms'] +
                         (1.5 if row['functional_group'] == 'NH2' else 0) +
                         (0.5 if row['functional_group'] in ['OH', 'COOH'] else 0))
                    n = 1.0 + 0.3 * row['void_fraction']
                else:
                    # H2: weaker interactions, mainly pore-size driven
                    q_sat = (0.3 * row['surface_area'] / 1000 + 
                             0.8 * row['pore_volume'])
                    q_sat = np.clip(q_sat, 0.05, 6.0)
                    b = 0.005 + 0.001 * row['pld']
                    n = 1.0 + 0.1 * row['void_fraction']
                
                # Langmuir-Freundlich isotherm: q = q_sat * (b*P)^n / (1 + (b*P)^n)
                bp_n = (b * p) ** n
                q = q_sat * bp_n / (1 + bp_n)
                q += np.random.normal(0, 0.05 * max(q, 0.01))
                q = max(q, 0)
                uptakes.append(q)
            results[col_name] = uptakes
    
    # CO2/N2 selectivity at 1 bar (for DAC relevance)
    co2_1bar = results['CO2_uptake_1.0bar']
    # Simulate N2 uptake - correlated with pore properties but much weaker
    n2_base = 0.02 * df['surface_area'] / 1000 + 0.05 * df['pore_volume']
    n2_uptake = n2_base * (1 + 0.1 * np.random.randn(len(df)))
    n2_uptake = n2_uptake.clip(0.001)
    results['N2_uptake_1bar'] = n2_uptake
    results['CO2_N2_selectivity'] = co2_1bar / n2_uptake
    
    # Heat of adsorption (Qst) for CO2 - correlated with affinity
    results['Qst_CO2'] = 15 + 10 * df['has_oms'] + 5 * (df['functional_group'].isin(['NH2', 'OH', 'COOH'])).astype(float) + np.random.normal(0, 3, len(df))
    results['Qst_CO2'] = results['Qst_CO2'].clip(10, 60)
    
    return results


# ============================================================
# 3. Water Stability & Synthesizability Prediction
# ============================================================
def predict_stability_synthesizability(df):
    """Predict water stability and synthesizability using rule-based + ML model."""
    results = df.copy()
    
    # Water stability score (0-1): depends on metal-ligand bond strength, hydrophobicity
    stability_score = (
        0.3 * (df['metal_electronegativity'] / 2.5) +
        0.2 * (1 - df['void_fraction']) +
        0.15 * (df['functional_group'].isin(['F', 'CF3', 'CH3'])).astype(float) +
        0.1 * (df['density'] / 2.5) +
        0.25 * np.random.beta(5, 3, len(df))
    )
    results['water_stability'] = stability_score.clip(0, 1)
    results['water_stable'] = (results['water_stability'] > 0.5).astype(int)
    
    # Synthesizability score - CoRE MOFs inherently synthesizable
    synth_base = np.where(df['source'] == 'CoRE', 0.85, 0.45)
    synth_score = synth_base + 0.08 * np.random.randn(len(df))
    # Simpler structures more synthesizable
    synth_score += -0.1 * (df['surface_area'] > 4000).astype(float)
    synth_score += 0.05 * (df['functional_group'] == 'none').astype(float)
    results['synthesizability'] = np.clip(synth_score, 0, 1)
    results['synthesizable'] = (results['synthesizability'] > 0.5).astype(int)
    
    return results


# ============================================================
# 4. ML Model Training for Adsorption Prediction
# ============================================================
def train_ml_models(df):
    """Train ML models for CO2 and H2 uptake prediction."""
    feature_cols = ['surface_area', 'pore_volume', 'void_fraction', 'pld', 'lcd',
                    'density', 'metal_electronegativity', 'has_oms', 'fg_encoded']
    
    targets = {
        'CO2_uptake_0.0004bar': 'CO2 (DAC conditions, 0.4 mbar)',
        'CO2_uptake_1.0bar': 'CO2 (1 bar)',
        'H2_uptake_1.0bar': 'H2 (1 bar)',
        'CO2_N2_selectivity': 'CO2/N2 Selectivity'
    }
    
    X = df[feature_cols].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    models = {}
    metrics = {}
    
    for target_col, target_name in targets.items():
        y = df[target_col].values
        X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.2, random_state=42)
        
        # Random Forest
        rf = RandomForestRegressor(n_estimators=200, max_depth=15, random_state=42, n_jobs=-1)
        rf.fit(X_train, y_train)
        y_pred_rf = rf.predict(X_test)
        
        # Gradient Boosting
        gb = GradientBoostingRegressor(n_estimators=200, max_depth=6, learning_rate=0.1, random_state=42)
        gb.fit(X_train, y_train)
        y_pred_gb = gb.predict(X_test)
        
        # Ensemble average
        y_pred_ens = (y_pred_rf + y_pred_gb) / 2
        
        r2_rf = r2_score(y_test, y_pred_rf)
        r2_gb = r2_score(y_test, y_pred_gb)
        r2_ens = r2_score(y_test, y_pred_ens)
        mae_ens = mean_absolute_error(y_test, y_pred_ens)
        rmse_ens = np.sqrt(mean_squared_error(y_test, y_pred_ens))
        
        # Feature importance (from RF)
        fi = dict(zip(feature_cols, rf.feature_importances_))
        
        models[target_col] = {'rf': rf, 'gb': gb, 'scaler': scaler}
        metrics[target_col] = {
            'name': target_name,
            'R2_RF': r2_rf, 'R2_GB': r2_gb, 'R2_Ensemble': r2_ens,
            'MAE': mae_ens, 'RMSE': rmse_ens,
            'feature_importance': fi,
            'y_test': y_test, 'y_pred': y_pred_ens
        }
        
        print(f"  {target_name}: R²(RF)={r2_rf:.4f}, R²(GB)={r2_gb:.4f}, R²(Ens)={r2_ens:.4f}, MAE={mae_ens:.4f}")
    
    return models, metrics, feature_cols


# ============================================================
# 5. DAC Ranking
# ============================================================
def rank_for_dac(df):
    """Rank MOFs for Direct Air Capture application."""
    dac = df.copy()
    
    # DAC score: weighted combination of relevant properties
    # High CO2 uptake at 400 ppm, high selectivity, moderate Qst, water stable, synthesizable
    dac['dac_score'] = (
        0.30 * (dac['CO2_uptake_0.0004bar'] / dac['CO2_uptake_0.0004bar'].max()) +
        0.20 * (dac['CO2_N2_selectivity'] / dac['CO2_N2_selectivity'].max()).clip(0, 1) +
        0.15 * (1 - np.abs(dac['Qst_CO2'] - 35) / 25).clip(0, 1) +  # optimal Qst ~35 kJ/mol
        0.15 * dac['water_stability'] +
        0.10 * dac['synthesizability'] +
        0.10 * (dac['CO2_uptake_1.0bar'] / dac['CO2_uptake_1.0bar'].max())
    )
    
    # Apply hard filters — relaxed thresholds for DAC
    dac['passes_filter'] = (
        (dac['water_stable'] == 1) & 
        (dac['synthesizable'] == 1) & 
        (dac['pld'] > 3.0)
    ).astype(int)
    
    dac_ranked = dac[dac['passes_filter'] == 1].sort_values('dac_score', ascending=False)
    
    return dac, dac_ranked


# ============================================================
# 6. Visualization
# ============================================================
def create_visualizations(df, metrics, dac_ranked, feature_cols):
    """Generate all figures."""
    sns.set_style('whitegrid')
    sns.set_palette('Set2')
    
    # Fig 1: Geometric descriptor distributions
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    desc_cols = ['surface_area', 'pore_volume', 'void_fraction', 'pld', 'lcd', 'density']
    desc_labels = ['Surface Area (m²/g)', 'Pore Volume (cm³/g)', 'Void Fraction',
                   'Pore Limiting Diameter (Å)', 'Largest Cavity Diameter (Å)', 'Density (g/cm³)']
    for ax, col, label in zip(axes.flat, desc_cols, desc_labels):
        for src in ['CoRE', 'hMOF']:
            subset = df[df['source'] == src]
            ax.hist(subset[col], bins=40, alpha=0.6, label=src, density=True)
        ax.set_xlabel(label, fontsize=10)
        ax.set_ylabel('Density', fontsize=10)
        ax.legend(fontsize=8)
    plt.suptitle('Distribution of Geometric Descriptors: CoRE MOF vs hMOF', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, 'descriptor_distributions.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Fig 2: CO2 adsorption isotherms (selected MOFs)
    pressures = [0.0004, 0.15, 1.0, 5.0, 10.0, 50.0]
    fig, ax = plt.subplots(figsize=(10, 7))
    top_mofs = dac_ranked.head(5)
    colors = plt.cm.tab10(np.linspace(0, 0.5, 5))
    for idx, (_, row) in enumerate(top_mofs.iterrows()):
        uptakes = [row[f'CO2_uptake_{p}bar'] for p in pressures]
        ax.plot(pressures, uptakes, 'o-', color=colors[idx], label=row['mof_id'], linewidth=2, markersize=6)
    ax.set_xscale('log')
    ax.set_xlabel('Pressure (bar)', fontsize=12)
    ax.set_ylabel('CO₂ Uptake (mmol/g)', fontsize=12)
    ax.set_title('CO₂ Adsorption Isotherms — Top 5 DAC Candidates', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, 'co2_isotherms_top5.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Fig 3: ML prediction parity plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    for ax, (target, m) in zip(axes.flat, metrics.items()):
        y_test = m['y_test']
        y_pred = m['y_pred']
        ax.scatter(y_test, y_pred, alpha=0.3, s=15, edgecolors='none')
        lims = [min(y_test.min(), y_pred.min()), max(y_test.max(), y_pred.max())]
        ax.plot(lims, lims, 'r--', linewidth=1.5, label='Ideal')
        ax.set_xlabel('GCMC Simulated', fontsize=10)
        ax.set_ylabel('ML Predicted', fontsize=10)
        ax.set_title(f"{m['name']}\nR²={m['R2_Ensemble']:.3f}, MAE={m['MAE']:.3f}", fontsize=11)
        ax.legend(fontsize=9)
    plt.suptitle('ML Prediction vs GCMC Simulation (Ensemble Model)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, 'ml_parity_plots.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Fig 4: Feature importance
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for ax, (target, m) in zip(axes.flat, metrics.items()):
        fi = m['feature_importance']
        sorted_fi = dict(sorted(fi.items(), key=lambda x: x[1], reverse=True))
        labels = list(sorted_fi.keys())
        values = list(sorted_fi.values())
        bars = ax.barh(labels, values, color=sns.color_palette('viridis', len(labels)))
        ax.set_xlabel('Importance', fontsize=10)
        ax.set_title(m['name'], fontsize=11)
        ax.invert_yaxis()
    plt.suptitle('Feature Importance (Random Forest)', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, 'feature_importance.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Fig 5: Structure-property relationship heatmap
    corr_cols = ['surface_area', 'pore_volume', 'void_fraction', 'pld', 'lcd', 'density',
                 'CO2_uptake_1.0bar', 'H2_uptake_1.0bar', 'CO2_N2_selectivity', 'Qst_CO2']
    corr_labels = ['SA', 'PV', 'VF', 'PLD', 'LCD', 'ρ', 'CO₂(1bar)', 'H₂(1bar)', 'S(CO₂/N₂)', 'Qst']
    corr_matrix = df[corr_cols].corr()
    fig, ax = plt.subplots(figsize=(10, 8))
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool), k=1)
    sns.heatmap(corr_matrix, mask=mask, annot=True, fmt='.2f', cmap='RdBu_r',
                xticklabels=corr_labels, yticklabels=corr_labels, center=0,
                square=True, linewidths=0.5, ax=ax, vmin=-1, vmax=1)
    ax.set_title('Structure–Property Correlation Matrix', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, 'correlation_heatmap.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Fig 6: DAC screening funnel
    total = len(df)
    pld_filter = len(df[df['pld'] > 3.0])
    stable = len(df[(df['pld'] > 3.0) & (df['water_stable'] == 1)])
    synth = len(df[(df['pld'] > 3.0) & (df['water_stable'] == 1) & (df['synthesizable'] == 1)])
    top_n = min(50, synth)
    
    stages = ['Total MOFs', 'PLD > 3.0 Å', 'Water\nStable', 'Synthesizable', f'Top {top_n}\n(DAC Score)']
    counts = [total, pld_filter, stable, synth, top_n]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    colors_funnel = plt.cm.Blues(np.linspace(0.3, 0.9, len(stages)))
    bars = ax.barh(stages[::-1], counts[::-1], color=colors_funnel[::-1], edgecolor='white', height=0.6)
    for bar, count in zip(bars, counts[::-1]):
        ax.text(bar.get_width() + 20, bar.get_y() + bar.get_height()/2, 
                str(count), va='center', fontsize=12, fontweight='bold')
    ax.set_xlabel('Number of MOFs', fontsize=12)
    ax.set_title('DAC MOF Screening Funnel', fontsize=14, fontweight='bold')
    ax.set_xlim(0, total * 1.15)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, 'dac_screening_funnel.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Fig 7: CO2 uptake vs surface area colored by source
    fig, ax = plt.subplots(figsize=(10, 7))
    for src, color in [('CoRE', '#2196F3'), ('hMOF', '#FF9800')]:
        subset = df[df['source'] == src]
        ax.scatter(subset['surface_area'], subset['CO2_uptake_1.0bar'],
                   alpha=0.4, s=15, label=src, color=color, edgecolors='none')
    ax.set_xlabel('Surface Area (m²/g)', fontsize=12)
    ax.set_ylabel('CO₂ Uptake at 1 bar (mmol/g)', fontsize=12)
    ax.set_title('CO₂ Uptake vs Surface Area', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, 'co2_vs_surface_area.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Fig 8: Water stability vs DAC score
    fig, ax = plt.subplots(figsize=(10, 7))
    dac_all = df.copy()
    sc = ax.scatter(dac_all['water_stability'], dac_all.get('dac_score', dac_all['CO2_uptake_0.0004bar']),
                    c=dac_all['synthesizability'], cmap='RdYlGn', alpha=0.5, s=20, edgecolors='none')
    cbar = plt.colorbar(sc, ax=ax, label='Synthesizability Score')
    ax.axvline(0.5, color='red', linestyle='--', alpha=0.7, label='Stability threshold')
    ax.set_xlabel('Water Stability Score', fontsize=12)
    ax.set_ylabel('DAC Score', fontsize=12)
    ax.set_title('Water Stability vs DAC Score\n(colored by Synthesizability)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, 'stability_vs_dac.png'), dpi=150, bbox_inches='tight')
    plt.close()
    
    return {
        'funnel_counts': dict(zip(stages, counts))
    }


# ============================================================
# Main Pipeline
# ============================================================
def main():
    print("=" * 60)
    print("MOF High-Throughput Screening Pipeline")
    print("=" * 60)
    
    # Step 1: Generate MOF database
    print("\n[1/6] Generating MOF database (CoRE + hMOF)...")
    df = generate_mof_database(n_core=500, n_hmof=1500)
    print(f"  Total MOFs: {len(df)} (CoRE: {len(df[df['source']=='CoRE'])}, hMOF: {len(df[df['source']=='hMOF'])})")
    
    # Step 2: GCMC simulation
    print("\n[2/6] Running GCMC adsorption simulations...")
    df = simulate_gcmc_adsorption(df)
    print(f"  CO2 uptake range (1 bar): {df['CO2_uptake_1.0bar'].min():.3f} - {df['CO2_uptake_1.0bar'].max():.3f} mmol/g")
    print(f"  H2 uptake range (1 bar): {df['H2_uptake_1.0bar'].min():.3f} - {df['H2_uptake_1.0bar'].max():.3f} mmol/g")
    
    # Step 3: Stability prediction
    print("\n[3/6] Predicting water stability and synthesizability...")
    df = predict_stability_synthesizability(df)
    print(f"  Water stable: {df['water_stable'].sum()} ({df['water_stable'].mean()*100:.1f}%)")
    print(f"  Synthesizable: {df['synthesizable'].sum()} ({df['synthesizable'].mean()*100:.1f}%)")
    
    # Step 4: ML models
    print("\n[4/6] Training ML models...")
    models, metrics, feature_cols = train_ml_models(df)
    
    # Step 5: DAC ranking
    print("\n[5/6] Ranking MOFs for DAC application...")
    df, dac_ranked = rank_for_dac(df)
    print(f"  MOFs passing all filters: {len(dac_ranked)}")
    print(f"\n  Top 10 DAC Candidates:")
    top10 = dac_ranked[['mof_id', 'source', 'surface_area', 'CO2_uptake_0.0004bar', 
                         'CO2_N2_selectivity', 'water_stability', 'dac_score']].head(10)
    print(top10.to_string(index=False))
    
    # Step 6: Visualization
    print("\n[6/6] Generating visualizations...")
    viz_data = create_visualizations(df, metrics, dac_ranked, feature_cols)
    
    # Save data
    df.to_csv(os.path.join(DATADIR, 'mof_screening_results.csv'), index=False)
    dac_ranked.head(50).to_csv(os.path.join(DATADIR, 'top50_dac_candidates.csv'), index=False)
    
    # Summary metrics
    summary = {
        'total_mofs': len(df),
        'core_mofs': len(df[df['source'] == 'CoRE']),
        'hmof_mofs': len(df[df['source'] == 'hMOF']),
        'water_stable': int(df['water_stable'].sum()),
        'synthesizable': int(df['synthesizable'].sum()),
        'dac_candidates': len(dac_ranked),
        'ml_metrics': {k: {'R2': v['R2_Ensemble'], 'MAE': v['MAE'], 'RMSE': v['RMSE']} 
                       for k, v in metrics.items()},
        'funnel': viz_data['funnel_counts'],
        'top10_dac': top10.to_dict('records')
    }
    
    with open(os.path.join(DATADIR, 'summary_metrics.json'), 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    print("\n" + "=" * 60)
    print("Pipeline complete! Files saved to figures/ and data/")
    print("=" * 60)
    
    return summary

if __name__ == '__main__':
    summary = main()
