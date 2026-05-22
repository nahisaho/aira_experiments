"""
Module 6: Great Barrier Reef 2100 Projection Scenarios
Integrated model combining all components for GBR-specific projections.
"""

import numpy as np
import json

from model_carbonate import carbonate_system, project_carbonate_scenarios
from model_calcification import calcification_rate_IpCC, calcification_pH_response
from model_compound_stress import compound_stress_model, thermal_performance_curve
from model_popgen import coral_adaptation_model


# === GBR Regional Parameters ===
GBR_REGIONS = {
    'Northern GBR': {
        'lat_range': (-10, -16),
        'T_baseline': 27.5,
        'pH_baseline': 8.05,
        'coral_cover_2020': 0.28,
        'dominant_species': 'Acropora',
        'bleaching_history': [1998, 2002, 2016, 2017, 2020, 2022],
    },
    'Central GBR': {
        'lat_range': (-16, -20),
        'T_baseline': 26.5,
        'pH_baseline': 8.08,
        'coral_cover_2020': 0.33,
        'dominant_species': 'Acropora',
        'bleaching_history': [2016, 2017, 2020],
    },
    'Southern GBR': {
        'lat_range': (-20, -24),
        'T_baseline': 25.0,
        'pH_baseline': 8.10,
        'coral_cover_2020': 0.36,
        'dominant_species': 'Porites',
        'bleaching_history': [2020],
    },
}


def project_regional_conditions(region_params, scenario, year):
    """Project temperature and pH for a GBR region under a given scenario."""
    delta_year = year - 2020
    
    temp_trajectories = {
        'RCP2.6': lambda dy: 0.01 * dy * (1 - dy/160),
        'RCP4.5': lambda dy: 0.0225 * dy,
        'RCP8.5': lambda dy: 0.0425 * dy + 0.0001 * dy**2,
    }
    
    ph_trajectories = {
        'RCP2.6': lambda dy: -0.001 * dy * (1 - dy/160),
        'RCP4.5': lambda dy: -0.0025 * dy,
        'RCP8.5': lambda dy: -0.004 * dy - 0.00001 * dy**2,
    }
    
    T = region_params['T_baseline'] + temp_trajectories[scenario](delta_year)
    pH = region_params['pH_baseline'] + ph_trajectories[scenario](delta_year)
    
    return T, pH


def estimate_coral_cover(T, pH, current_cover, species, adaptation_factor=1.0):
    """
    Estimate coral cover change based on compound stress and calcification.
    """
    # Compound stress performance
    T_opt = 27.0 if species == 'Acropora' else 27.5
    synergy = 1.8 if species == 'Acropora' else 1.2
    
    perf, components = compound_stress_model(
        T, pH, T_opt=T_opt, synergy_factor=synergy, model_type='synergistic'
    )
    perf = float(perf)
    
    # Calcification capacity
    from model_carbonate import compute_K1, compute_K2, compute_Ksp_aragonite
    T_K = T + 273.15
    S = 35.0
    K1 = compute_K1(T_K, S)
    K2 = compute_K2(T_K, S)
    Ksp = compute_Ksp_aragonite(T_K, S)
    H = 10**(-pH)
    # Approximate Omega from pH
    Ca = 0.01028
    # Rough DIC estimate
    DIC = 2050e-6  # mol/kg
    CO3 = DIC * K1*K2/H**2 / (1 + K1/H + K1*K2/H**2)
    Omega = Ca * CO3 / Ksp
    
    G_rel = float(calcification_rate_IpCC(Omega, species))
    G_ref = float(calcification_rate_IpCC(3.5, species))
    calc_ratio = G_rel / G_ref if G_ref > 0 else 0
    calc_ratio = max(0, min(calc_ratio, 1.5))
    
    # Adaptive capacity
    adaptive_boost = adaptation_factor
    
    # Integrated cover change
    growth_potential = perf * calc_ratio * adaptive_boost
    
    # Logistic cover dynamics
    K_cover = 0.60  # maximum possible cover
    cover = current_cover * growth_potential
    cover = min(cover, K_cover)
    cover = max(cover, 0.01)
    
    return cover, {
        'performance': round(perf, 4),
        'Omega_aragonite': round(float(Omega), 3),
        'calcification_ratio': round(calc_ratio, 4),
        'p_bleach': round(float(components['p_bleach']), 4),
        'growth_potential': round(growth_potential, 4),
    }


def run_gbr_projection():
    """Run complete GBR projection for all regions and scenarios."""
    scenarios = ['RCP2.6', 'RCP4.5', 'RCP8.5']
    years = list(range(2020, 2101, 10))
    
    # Get adaptation factors from popgen model
    # Simplified: RCP2.6 allows full adaptation, RCP8.5 overwhelms
    adaptation_factors = {
        'RCP2.6': {'Acropora': 1.05, 'Porites': 1.08},
        'RCP4.5': {'Acropora': 0.95, 'Porites': 1.02},
        'RCP8.5': {'Acropora': 0.80, 'Porites': 0.90},
    }
    
    projection_results = {}
    
    for region_name, region_params in GBR_REGIONS.items():
        projection_results[region_name] = {}
        
        for scenario in scenarios:
            scenario_timeline = []
            current_cover = region_params['coral_cover_2020']
            species = region_params['dominant_species']
            
            for year in years:
                T, pH = project_regional_conditions(region_params, scenario, year)
                adapt_factor = adaptation_factors[scenario].get(species, 1.0)
                
                # Progressive adaptation (increases slightly over time)
                year_factor = 1 + (adapt_factor - 1) * min((year - 2020) / 80, 1.0)
                
                cover, details = estimate_coral_cover(
                    T, pH, current_cover, species, year_factor
                )
                
                # Marine heatwave events (stochastic)
                np.random.seed(year * 100 + hash(region_name) % 1000)
                mhw_prob = 0.05 + 0.003 * (year - 2020)  # increasing frequency
                if scenario == 'RCP8.5':
                    mhw_prob *= 2
                elif scenario == 'RCP2.6':
                    mhw_prob *= 0.7
                
                if np.random.random() < mhw_prob:
                    mhw_loss = np.random.uniform(0.1, 0.4)
                    cover *= (1 - mhw_loss)
                    details['mhw_event'] = True
                    details['mhw_loss_pct'] = round(mhw_loss * 100, 1)
                else:
                    details['mhw_event'] = False
                
                # Recovery between decades (partial)
                if year > 2020:
                    recovery_rate = 0.02 * details['performance']  # slow recovery
                    cover = min(cover + recovery_rate, region_params['coral_cover_2020'])
                
                current_cover = max(cover, 0.01)
                
                entry = {
                    'year': year,
                    'T_C': round(T, 2),
                    'pH': round(pH, 3),
                    'coral_cover': round(current_cover, 4),
                    **details
                }
                scenario_timeline.append(entry)
            
            projection_results[region_name][scenario] = scenario_timeline
    
    # Summary statistics for 2100
    summary_2100 = {}
    for region in projection_results:
        summary_2100[region] = {}
        for scenario in scenarios:
            timeline = projection_results[region][scenario]
            entry_2100 = timeline[-1]
            entry_2020 = timeline[0]
            summary_2100[region][scenario] = {
                'cover_2020': entry_2020['coral_cover'],
                'cover_2100': entry_2100['coral_cover'],
                'cover_change_pct': round((entry_2100['coral_cover'] / entry_2020['coral_cover'] - 1) * 100, 1),
                'T_2100': entry_2100['T_C'],
                'pH_2100': entry_2100['pH'],
                'Omega_2100': entry_2100['Omega_aragonite'],
                'n_mhw_events': sum(1 for e in timeline if e.get('mhw_event', False)),
            }
    
    # Ecosystem tipping points
    tipping_points = {}
    for region in projection_results:
        tipping_points[region] = {}
        for scenario in scenarios:
            timeline = projection_results[region][scenario]
            tp_year = None
            for entry in timeline:
                if entry['coral_cover'] < 0.10:
                    tp_year = entry['year']
                    break
            tipping_points[region][scenario] = {
                'below_10pct_year': tp_year,
                'functional_collapse': tp_year is not None
            }
    
    all_results = {
        'projection_timelines': projection_results,
        'summary_2100': summary_2100,
        'tipping_points': tipping_points,
        'gbr_regions': {k: {kk: vv for kk, vv in v.items() if kk != 'bleaching_history'}
                        for k, v in GBR_REGIONS.items()}
    }
    
    with open('results/gbr_projection.json', 'w') as f:
        json.dump(all_results, f, indent=2)
    
    return all_results


if __name__ == '__main__':
    results = run_gbr_projection()
    print("=== GBR 2100 Projection Results ===")
    print("\nSummary for 2100:")
    for region, scenarios in results['summary_2100'].items():
        print(f"\n{region}:")
        for scen, data in scenarios.items():
            print(f"  {scen}: cover {data['cover_2020']:.1%} → {data['cover_2100']:.1%} "
                  f"({data['cover_change_pct']:+.1f}%), "
                  f"T={data['T_2100']:.1f}°C, pH={data['pH_2100']:.2f}, "
                  f"Ω={data['Omega_2100']:.2f}")
    
    print("\nTipping Points:")
    for region, scenarios in results['tipping_points'].items():
        for scen, data in scenarios.items():
            if data['functional_collapse']:
                print(f"  {region} [{scen}]: collapse by {data['below_10pct_year']}")
    
    print("\nResults saved to results/gbr_projection.json")
