"""
Module 2: Coral Calcification Rate Model
pH/Ω-dependent calcification modeling based on empirical relationships.
"""

import numpy as np
import json


def calcification_rate_IpCC(Omega_arag, species='Acropora'):
    """
    Calcification rate as a function of aragonite saturation state.
    Based on the IpCC (Internal pH-controlled Calcification) model
    (Comeau et al. 2013; Holcomb et al. 2014).
    
    G = k * (Omega_arag - 1)^n
    
    Parameters
    ----------
    Omega_arag : float or array - Aragonite saturation state
    species : str - Coral species/genus
    
    Returns
    -------
    G : float or array - Calcification rate (µmol CaCO3 cm⁻² h⁻¹)
    """
    # Species-specific parameters (from meta-analysis)
    params = {
        'Acropora':     {'k': 2.8, 'n': 1.7, 'G_max': 12.0, 'threshold': 1.0},
        'Porites':      {'k': 1.5, 'n': 1.3, 'G_max': 6.0,  'threshold': 0.8},
        'Montipora':    {'k': 2.2, 'n': 1.5, 'G_max': 9.0,  'threshold': 0.9},
        'Stylophora':   {'k': 2.0, 'n': 1.4, 'G_max': 8.0,  'threshold': 1.0},
        'Pavona':       {'k': 1.8, 'n': 1.2, 'G_max': 7.0,  'threshold': 0.7},
        'CCA':          {'k': 3.5, 'n': 2.0, 'G_max': 5.0,  'threshold': 1.5},
    }

    p = params.get(species, params['Acropora'])
    Omega_arag = np.asarray(Omega_arag, dtype=float)
    G = np.where(
        Omega_arag > p['threshold'],
        p['k'] * (Omega_arag - p['threshold'])**p['n'],
        -0.5 * (p['threshold'] - Omega_arag)  # dissolution
    )
    G = np.minimum(G, p['G_max'])
    return G


def calcification_pH_response(pH, species='Acropora'):
    """
    Calcification response to pH changes.
    Empirical sigmoid model based on meta-analysis.
    
    G_relative = 1 / (1 + exp(-s * (pH - pH_crit)))
    """
    params = {
        'Acropora':   {'pH_crit': 7.7, 's': 8.0, 'G_ref': 1.0},
        'Porites':    {'pH_crit': 7.5, 's': 6.5, 'G_ref': 1.0},
        'Montipora':  {'pH_crit': 7.6, 's': 7.0, 'G_ref': 1.0},
        'Stylophora': {'pH_crit': 7.65, 's': 7.5, 'G_ref': 1.0},
        'Pavona':     {'pH_crit': 7.55, 's': 6.0, 'G_ref': 1.0},
        'CCA':        {'pH_crit': 7.9, 's': 10.0, 'G_ref': 1.0},
    }
    p = params.get(species, params['Acropora'])
    pH = np.asarray(pH, dtype=float)
    G_rel = 1.0 / (1.0 + np.exp(-p['s'] * (pH - p['pH_crit'])))
    return G_rel


def net_reef_accretion(G_calcification, bioerosion_rate=0.3, dissolution_rate=0.1):
    """
    Net reef accretion = calcification - bioerosion - dissolution.
    
    Parameters
    ----------
    G_calcification : float - Gross calcification (kg CaCO3 m⁻² yr⁻¹)
    bioerosion_rate : float - Bioerosion rate (kg CaCO3 m⁻² yr⁻¹)
    dissolution_rate : float - Chemical dissolution rate
    
    Returns
    -------
    net_accretion : float (kg CaCO3 m⁻² yr⁻¹)
    """
    return G_calcification - bioerosion_rate - dissolution_rate


def run_calcification_analysis():
    """Run comprehensive calcification analysis and save results."""
    species_list = ['Acropora', 'Porites', 'Montipora', 'Stylophora', 'Pavona', 'CCA']
    
    # Omega range
    Omega_range = np.linspace(0.5, 5.0, 100)
    pH_range = np.linspace(7.2, 8.4, 100)
    
    results = {
        'omega_response': {},
        'pH_response': {},
        'critical_thresholds': {},
        'net_accretion_scenarios': {}
    }
    
    for sp in species_list:
        # Omega response
        G_omega = calcification_rate_IpCC(Omega_range, sp)
        results['omega_response'][sp] = {
            'Omega': Omega_range.tolist(),
            'G': G_omega.tolist()
        }
        
        # pH response
        G_pH = calcification_pH_response(pH_range, sp)
        results['pH_response'][sp] = {
            'pH': pH_range.tolist(),
            'G_relative': G_pH.tolist()
        }
        
        # Critical thresholds
        G_present = calcification_rate_IpCC(3.5, sp)
        G_rcp85_2100 = calcification_rate_IpCC(1.8, sp)
        reduction = (1 - G_rcp85_2100 / G_present) * 100 if G_present > 0 else 100
        
        results['critical_thresholds'][sp] = {
            'G_present_day_Omega3.5': round(float(G_present), 3),
            'G_RCP85_2100_Omega1.8': round(float(G_rcp85_2100), 3),
            'percent_reduction': round(float(reduction), 1),
            'Omega_net_dissolution': round(float(
                Omega_range[np.argmin(np.abs(calcification_rate_IpCC(Omega_range, sp)))]
            ), 2)
        }
    
    # Net accretion under different scenarios
    for scenario_name, Omega_val in [('Present', 3.5), ('RCP4.5_2050', 2.8),
                                       ('RCP4.5_2100', 2.3), ('RCP8.5_2050', 2.5),
                                       ('RCP8.5_2100', 1.8)]:
        accretions = {}
        for sp in species_list:
            G = float(calcification_rate_IpCC(Omega_val, sp))
            # Scale to kg/m²/yr (approximate conversion)
            G_annual = G * 8760 * 1e-2 * 100  # rough conversion
            G_annual_kg = G_annual * 1e-4  # to kg
            net = net_reef_accretion(G_annual_kg)
            accretions[sp] = {
                'gross_calcification': round(G_annual_kg, 3),
                'net_accretion': round(net, 3)
            }
        results['net_accretion_scenarios'][scenario_name] = accretions
    
    with open('results/calcification_model.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    return results


if __name__ == '__main__':
    results = run_calcification_analysis()
    print("=== Calcification Model Results ===")
    print("\nCritical Thresholds:")
    for sp, th in results['critical_thresholds'].items():
        print(f"  {sp}: Present G={th['G_present_day_Omega3.5']:.3f}, "
              f"RCP8.5 2100 G={th['G_RCP85_2100_Omega1.8']:.3f}, "
              f"Reduction={th['percent_reduction']:.1f}%")
    print("\nResults saved to results/calcification_model.json")
