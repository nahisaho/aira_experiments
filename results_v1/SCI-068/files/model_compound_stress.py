"""
Module 4: Temperature-pH Compound Stress Synergy Model
Models the interactive (synergistic/antagonistic) effects of warming and acidification.
"""

import numpy as np
from scipy.optimize import minimize_scalar
import json


def thermal_performance_curve(T, T_opt=27.0, sigma=3.5, asymmetry=0.7):
    """
    Asymmetric thermal performance curve for coral physiological performance.
    Based on Norberg (2004) thermal niche model.
    
    P(T) = exp(-(T - T_opt)² / (2 * sigma_eff²))
    where sigma_eff is asymmetric (narrower above T_opt).
    """
    T = np.asarray(T, dtype=float)
    sigma_eff = np.where(T <= T_opt, sigma, sigma * asymmetry)
    P = np.exp(-((T - T_opt)**2) / (2 * sigma_eff**2))
    return P


def bleaching_probability(T, T_opt=27.0, DHW_threshold=4.0, duration_weeks=1):
    """
    Bleaching probability based on Degree Heating Weeks (DHW).
    
    Parameters
    ----------
    T : float - Current temperature (°C)
    T_opt : float - MMM (Maximum Monthly Mean) temperature
    DHW_threshold : float - DHW threshold for bleaching
    duration_weeks : int - Duration of thermal stress (weeks)
    """
    T = np.asarray(T, dtype=float)
    hotspot = np.maximum(T - (T_opt + 1.0), 0)
    DHW = hotspot * duration_weeks
    
    # Logistic bleaching probability
    p_bleach = 1.0 / (1.0 + np.exp(-2.5 * (DHW - DHW_threshold)))
    return p_bleach, DHW


def ph_stress_function(pH, pH_opt=8.1, pH_crit=7.6):
    """
    pH-dependent stress function (0 = max stress, 1 = no stress).
    Sigmoid decline below optimal pH.
    """
    pH = np.asarray(pH, dtype=float)
    stress = 1.0 / (1.0 + np.exp(-8.0 * (pH - pH_crit)))
    return stress


def compound_stress_model(T, pH, T_opt=27.0, pH_opt=8.1, 
                           synergy_factor=1.5, model_type='multiplicative'):
    """
    Compound stress model combining thermal and pH stress.
    
    Models:
    - 'additive': S_total = S_T + S_pH
    - 'multiplicative': S_total = S_T × S_pH  
    - 'synergistic': S_total = S_T × S_pH × (1 + γ × S_T × S_pH)
    
    Parameters
    ----------
    T : float - Temperature (°C)
    pH : float - Seawater pH
    synergy_factor : float - γ, synergy coefficient (>0 = synergistic, <0 = antagonistic)
    
    Returns
    -------
    performance : float - Relative performance (0-1)
    stress_components : dict
    """
    T = np.asarray(T, dtype=float)
    pH = np.asarray(pH, dtype=float)
    
    # Individual stress components
    P_thermal = thermal_performance_curve(T, T_opt)
    P_pH = ph_stress_function(pH, pH_opt)
    
    p_bleach, DHW = bleaching_probability(T, T_opt)
    
    if model_type == 'additive':
        S_T = 1 - P_thermal
        S_pH = 1 - P_pH
        S_total = np.minimum(S_T + S_pH, 1.0)
        performance = 1 - S_total
    
    elif model_type == 'multiplicative':
        performance = P_thermal * P_pH
    
    elif model_type == 'synergistic':
        S_T = 1 - P_thermal
        S_pH = 1 - P_pH
        synergy_term = synergy_factor * S_T * S_pH
        S_total = np.minimum(S_T + S_pH + synergy_term, 1.0)
        performance = np.maximum(1 - S_total, 0)
    
    else:
        raise ValueError(f"Unknown model_type: {model_type}")
    
    # Additional bleaching impact
    performance = performance * (1 - p_bleach * 0.8)  # Bleaching reduces performance
    
    return performance, {
        'P_thermal': P_thermal,
        'P_pH': P_pH,
        'p_bleach': p_bleach,
        'DHW': DHW,
        'performance': performance
    }


def species_vulnerability_assessment(T_range, pH_range):
    """
    Assess vulnerability of different coral species to compound stress.
    """
    species_params = {
        'Acropora (branching)': {
            'T_opt': 27.0, 'sigma': 3.0, 'asymmetry': 0.6,
            'pH_crit': 7.7, 'synergy': 1.8,
            'bleaching_threshold': 3.0
        },
        'Porites (massive)': {
            'T_opt': 27.5, 'sigma': 4.0, 'asymmetry': 0.8,
            'pH_crit': 7.5, 'synergy': 1.2,
            'bleaching_threshold': 6.0
        },
        'Montipora (encrusting)': {
            'T_opt': 27.0, 'sigma': 3.5, 'asymmetry': 0.7,
            'pH_crit': 7.6, 'synergy': 1.5,
            'bleaching_threshold': 4.0
        },
        'Stylophora (branching)': {
            'T_opt': 26.5, 'sigma': 3.0, 'asymmetry': 0.65,
            'pH_crit': 7.65, 'synergy': 1.6,
            'bleaching_threshold': 3.5
        },
    }
    
    T_grid, pH_grid = np.meshgrid(T_range, pH_range)
    
    results = {}
    for name, params in species_params.items():
        perf, components = compound_stress_model(
            T_grid, pH_grid,
            T_opt=params['T_opt'],
            pH_opt=8.1,
            synergy_factor=params['synergy'],
            model_type='synergistic'
        )
        
        # Find tipping points
        critical_mask = perf < 0.3
        
        results[name] = {
            'performance_matrix': perf.tolist(),
            'mean_performance': round(float(np.mean(perf)), 4),
            'min_performance': round(float(np.min(perf)), 4),
            'fraction_critical': round(float(np.mean(critical_mask)), 4),
            'params': params
        }
    
    return results


def run_compound_stress_analysis():
    """Run full compound stress analysis."""
    T_range = np.linspace(22, 34, 50)
    pH_range = np.linspace(7.2, 8.4, 50)
    
    # Species vulnerability
    vuln_results = species_vulnerability_assessment(T_range, pH_range)
    
    # Compare stress models
    T_test = np.linspace(24, 33, 30)
    pH_test = np.linspace(7.4, 8.2, 30)
    T_g, pH_g = np.meshgrid(T_test, pH_test)
    
    model_comparison = {}
    for model_type in ['additive', 'multiplicative', 'synergistic']:
        perf, _ = compound_stress_model(T_g, pH_g, model_type=model_type)
        model_comparison[model_type] = {
            'mean_performance': round(float(np.mean(perf)), 4),
            'std_performance': round(float(np.std(perf)), 4),
            'fraction_below_30pct': round(float(np.mean(perf < 0.3)), 4)
        }
    
    # GBR-specific scenarios
    gbr_scenarios = {
        'Present (2020)': {'T': 25.5, 'pH': 8.07},
        'RCP4.5 2050':    {'T': 26.5, 'pH': 7.95},
        'RCP4.5 2100':    {'T': 27.3, 'pH': 7.85},
        'RCP8.5 2050':    {'T': 27.0, 'pH': 7.90},
        'RCP8.5 2100':    {'T': 29.2, 'pH': 7.70},
        'RCP8.5 2100+MHW':{'T': 31.5, 'pH': 7.65},
    }
    
    scenario_results = {}
    species_names = list(vuln_results.keys())
    for scen_name, conditions in gbr_scenarios.items():
        scenario_results[scen_name] = {}
        for sp_name in species_names:
            params = vuln_results[sp_name]['params']
            perf, components = compound_stress_model(
                conditions['T'], conditions['pH'],
                T_opt=params['T_opt'],
                synergy_factor=params['synergy'],
                model_type='synergistic'
            )
            scenario_results[scen_name][sp_name] = {
                'performance': round(float(perf), 4),
                'P_thermal': round(float(components['P_thermal']), 4),
                'P_pH': round(float(components['P_pH']), 4),
                'p_bleach': round(float(components['p_bleach']), 4)
            }
    
    all_results = {
        'model_comparison': model_comparison,
        'species_vulnerability': {k: {kk: vv for kk, vv in v.items() if kk != 'performance_matrix'} 
                                   for k, v in vuln_results.items()},
        'gbr_scenarios': scenario_results,
        'T_range': T_range.tolist(),
        'pH_range': pH_range.tolist()
    }
    
    with open('results/compound_stress.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    
    # Save performance matrices for plotting
    perf_matrices = {}
    for sp_name in vuln_results:
        perf_matrices[sp_name] = vuln_results[sp_name]['performance_matrix']
    
    np.savez('results/stress_matrices.npz',
             T_range=T_range, pH_range=pH_range,
             **{k.replace(' ', '_').replace('(', '').replace(')', ''): np.array(v)
                for k, v in perf_matrices.items()})
    
    return all_results


if __name__ == '__main__':
    results = run_compound_stress_analysis()
    print("=== Compound Stress Analysis ===")
    print("\nModel Comparison:")
    for model, stats in results['model_comparison'].items():
        print(f"  {model}: mean_perf={stats['mean_performance']:.4f}, "
              f"critical_fraction={stats['fraction_below_30pct']:.4f}")
    print("\nGBR Scenario Results:")
    for scen, species_data in results['gbr_scenarios'].items():
        print(f"\n  {scen}:")
        for sp, data in species_data.items():
            print(f"    {sp}: performance={data['performance']:.3f}")
    print("\nResults saved to results/compound_stress.json")
