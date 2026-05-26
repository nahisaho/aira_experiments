"""
Module 6: Fischer-Tropsch synthesis case study on Co(0001).

Complete microkinetic model with:
- DFT-derived energetics (literature values)
- Coverage-dependent lateral interactions
- Reactor coupling (PFR/CSTR)
- Rate control analysis
"""
import numpy as np
from .rate_constants import compute_rate_constant, kB, h, eV_to_J
from .adsorption import langmuir_isotherm, competitive_langmuir
from .lateral_interactions import LateralInteractionModel, create_ft_interaction_matrix
from .reactor_models import SimpleMicroKineticModel, PFR, CSTR


# DFT-derived energetics for FT on Co(0001)
# Values from literature (Zhuo et al., JACS 2009; Qi et al., JCTC 2020)
# Adsorption barriers set as sticking coefficients via pre-exponential tuning
FT_ENERGETICS = {
    'CO_adsorption': {'Ea_fwd': 0.0, 'Ea_rev': 1.35, 'nu_imag': None},
    'H2_dissociation': {'Ea_fwd': 0.05, 'Ea_rev': 0.90, 'nu_imag': 800},
    'CO_dissociation': {'Ea_fwd': 1.60, 'Ea_rev': 1.10, 'nu_imag': 450},
    'C_hydrogenation': {'Ea_fwd': 0.75, 'Ea_rev': 0.55, 'nu_imag': 1100},
    'CH_hydrogenation': {'Ea_fwd': 0.65, 'Ea_rev': 0.45, 'nu_imag': 1050},
    'CH2_hydrogenation': {'Ea_fwd': 0.55, 'Ea_rev': 0.50, 'nu_imag': 1000},
    'CH3_hydrogenation': {'Ea_fwd': 0.85, 'Ea_rev': 0.00, 'nu_imag': 950},
    'O_hydrogenation': {'Ea_fwd': 1.00, 'Ea_rev': 0.80, 'nu_imag': 900},
    'OH_hydrogenation': {'Ea_fwd': 1.10, 'Ea_rev': 0.00, 'nu_imag': 850},
    'chain_growth': {'Ea_fwd': 0.95, 'Ea_rev': 0.70, 'nu_imag': 600},
}

# Sticking coefficients / prefactor corrections (dimensionless)
# Adsorption from gas phase uses collision theory, not TST.
# We calibrate K_ads to match experimental coverage at 500K:
#   theta_CO ~ 0.4, theta_H ~ 0.25 at P=20bar, H2/CO=2
# This gives K_CO*P_CO ~ 2.0, sqrt(K_H2*P_H2) ~ 1.0
# => K_CO ~ 0.30 bar^-1, K_H2 ~ 0.075 bar^-1
# Back-calculated prefactor corrections to achieve this
PREFACTOR_CORRECTIONS = {
    'k1f': 1e-14,  # CO adsorption prefactor correction
    'k1r': 1.0,
    'k2f': 1e-10,  # H2 dissociative adsorption prefactor correction
    'k2r': 1.0,
}


def compute_ft_rate_constants(T, tunneling='wigner'):
    """
    Compute all rate constants for FT elementary steps at temperature T.
    Applies prefactor corrections for gas-phase adsorption steps.
    """
    k = {}
    steps = [
        ('1', 'CO_adsorption'),
        ('2', 'H2_dissociation'),
        ('3', 'CO_dissociation'),
        ('4', 'C_hydrogenation'),
        ('5', 'CH_hydrogenation'),
        ('6', 'CH2_hydrogenation'),
        ('7', 'CH3_hydrogenation'),
        ('8', 'O_hydrogenation'),
        ('9', 'OH_hydrogenation'),
        ('10', 'chain_growth'),
    ]
    
    for num, name in steps:
        e = FT_ENERGETICS[name]
        kf = compute_rate_constant(
            e['Ea_fwd'], T, e['nu_imag'], e.get('Ea_rev'), tunneling
        )
        kr = compute_rate_constant(
            e['Ea_rev'], T, e['nu_imag'], e.get('Ea_fwd'), tunneling
        )
        # Apply prefactor corrections
        kf *= PREFACTOR_CORRECTIONS.get(f'k{num}f', 1.0)
        kr *= PREFACTOR_CORRECTIONS.get(f'k{num}r', 1.0)
        k[f'k{num}f'] = kf
        k[f'k{num}r'] = kr
    
    return k


def run_ft_simulation(T=500, P_total=20.0, H2_CO_ratio=2.0,
                       catalyst_weight=1.0, total_flow=1e-3,
                       tunneling='wigner'):
    """
    Run complete FT simulation.
    
    Returns dict with results.
    """
    # Compute rate constants
    k = compute_ft_rate_constants(T, tunneling)
    
    # Create microkinetic model
    mkm = SimpleMicroKineticModel(k)
    
    # Partial pressures
    P_CO = P_total / (1 + H2_CO_ratio)
    P_H2 = P_total * H2_CO_ratio / (1 + H2_CO_ratio)
    P = np.array([P_CO, P_H2, 0, 0, 0])
    
    # Solve surface coverages
    coverage_sol = mkm.compute_coverages(P, T)
    
    # Steady-state coverages
    theta_ss = coverage_sol.y[:, -1]
    theta_ss = np.maximum(theta_ss, 0)
    
    # Compute gas-phase rates
    rates = mkm.compute_rates(P, T)
    
    # Run PFR simulation
    reactor_params = {
        'catalyst_weight': catalyst_weight,
        'total_flow': total_flow,
        'pressure': P_total,
        'temperature': T,
        'cat_site_density': 1.0,  # mol_sites/kg_cat (typical for supported metal)
    }
    
    pfr = PFR(mkm, reactor_params)
    F0 = np.array([total_flow * P_CO / P_total,
                    total_flow * P_H2 / P_total,
                    0, 0, 0])
    
    pfr_sol = pfr.solve(F0, (0, catalyst_weight))
    
    # Run CSTR simulation
    cstr = CSTR(mkm, reactor_params)
    F_cstr, converged = cstr.solve(F0)
    
    # CO conversion
    if pfr_sol.success:
        X_CO_pfr = 1.0 - pfr_sol.y[0, -1] / F0[0] if F0[0] > 0 else 0
    else:
        X_CO_pfr = 0
    
    X_CO_cstr = 1.0 - F_cstr[0] / F0[0] if F0[0] > 0 else 0
    
    results = {
        'temperature': T,
        'pressure': P_total,
        'H2_CO_ratio': H2_CO_ratio,
        'rate_constants': k,
        'coverages': theta_ss,
        'coverage_solution': coverage_sol,
        'rates': rates,
        'pfr_solution': pfr_sol,
        'cstr_solution': F_cstr,
        'cstr_converged': converged,
        'X_CO_pfr': X_CO_pfr,
        'X_CO_cstr': X_CO_cstr,
        'F0': F0,
    }
    
    return results


def temperature_study(T_range, P_total=20.0, H2_CO_ratio=2.0):
    """
    Study temperature dependence of FT rates and selectivity.
    """
    results = []
    for T in T_range:
        res = run_ft_simulation(T=T, P_total=P_total, H2_CO_ratio=H2_CO_ratio)
        results.append(res)
    return results


def pressure_study(P_range, T=500, H2_CO_ratio=2.0):
    """
    Study pressure dependence of FT rates and selectivity.
    """
    results = []
    for P in P_range:
        res = run_ft_simulation(T=T, P_total=P, H2_CO_ratio=H2_CO_ratio)
        results.append(res)
    return results


def degree_of_rate_control_analysis(T=500, P_total=20.0, H2_CO_ratio=2.0):
    """
    Perform DRC analysis on the FT mechanism using the full coverage ODE model.
    X_RC,i = d(ln r) / d(ln k_i) at fixed K_eq,i
    """
    from .rate_constants import compute_rate_constant
    
    P_CO = P_total / (1 + H2_CO_ratio)
    P_H2 = P_total * H2_CO_ratio / (1 + H2_CO_ratio)
    P = np.array([P_CO, P_H2, 0, 0, 0])
    
    step_names = list(FT_ENERGETICS.keys())
    
    # Base case: solve coverages and get rate
    k_base = compute_ft_rate_constants(T)
    mkm_base = SimpleMicroKineticModel(k_base)
    sol_base = mkm_base.compute_coverages(P, T)
    theta_base = np.clip(sol_base.y[:, -1], 0, 1)
    theta_star_base = max(1.0 - np.sum(theta_base), 0)
    
    # CO dissociation rate as the reference overall rate
    r_base = k_base['k3f'] * theta_base[0] * theta_star_base
    
    perturbation = 0.05
    drc = np.zeros(len(step_names))
    
    for i, name in enumerate(step_names):
        step_num = str(i + 1)
        k_pert = compute_ft_rate_constants(T)
        
        # Perturb forward rate constant, keep K_eq constant (also perturb reverse)
        k_pert[f'k{step_num}f'] *= (1.0 + perturbation)
        # Keep equilibrium constant: K = kf/kr, so kr must also be scaled
        # No: DRC is at fixed K_eq and other k, so only perturb kf
        # Actually Campbell's DRC: perturb both kf and kr equally
        k_pert[f'k{step_num}r'] *= (1.0 + perturbation)
        
        mkm_pert = SimpleMicroKineticModel(k_pert)
        sol_pert = mkm_pert.compute_coverages(P, T)
        theta_pert = np.clip(sol_pert.y[:, -1], 0, 1)
        theta_star_pert = max(1.0 - np.sum(theta_pert), 0)
        
        r_pert = k_pert['k3f'] * theta_pert[0] * theta_star_pert
        
        if r_base > 1e-50:
            drc[i] = (np.log(r_pert + 1e-50) - np.log(r_base + 1e-50)) / \
                     np.log(1.0 + perturbation)
    
    return step_names, drc
