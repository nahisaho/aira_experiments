"""
Module 2: Adsorption isotherm models
- Langmuir
- Temkin
- Fractal surface (Freundlich-like with fractal dimension)
"""
import numpy as np


def langmuir_isotherm(P, K):
    """
    Langmuir isotherm: theta = K*P / (1 + K*P)
    P: partial pressure [bar]
    K: adsorption equilibrium constant [bar^-1]
    """
    return K * P / (1.0 + K * P)


def competitive_langmuir(P_dict, K_dict):
    """
    Competitive Langmuir adsorption for multiple species.
    P_dict: {species: partial_pressure}
    K_dict: {species: equilibrium_constant}
    Returns: {species: coverage}
    """
    denom = 1.0 + sum(K_dict[s] * P_dict[s] for s in P_dict)
    coverages = {s: K_dict[s] * P_dict[s] / denom for s in P_dict}
    return coverages


def temkin_isotherm(P, K0, alpha, theta_range=None):
    """
    Temkin isotherm: accounts for linear decrease of adsorption energy with coverage.
    dH_ads(theta) = dH_ads_0 * (1 - alpha * theta)
    
    Implicit form solved numerically.
    K0: adsorption equilibrium constant at zero coverage
    alpha: interaction parameter (dimensionless, typically 0-2)
    """
    if theta_range is None:
        theta_range = np.linspace(0.01, 0.99, 200)
    
    # Temkin: K(theta) = K0 * exp(-alpha * theta)
    # theta = K(theta)*P / (1 + K(theta)*P)
    P_calc = theta_range / (K0 * np.exp(-alpha * theta_range) * (1.0 - theta_range))
    return P_calc, theta_range


def fractal_isotherm(P, K, D, d=2):
    """
    Fractal surface isotherm (generalized Freundlich):
      theta = (K*P)^(1/n)  where n relates to fractal dimension D
      n = D/d for a surface with fractal dimension D embedded in d+1 space
    
    D: fractal dimension (2 <= D <= 3 for surfaces)
    d: topological dimension of smooth surface (default 2)
    """
    n = D / d
    theta = (K * P) ** (1.0 / n)
    theta = np.clip(theta, 0, 1)
    return theta


def coverage_dependent_binding_energy(E0, theta, epsilon, model='linear'):
    """
    Coverage-dependent binding energy with lateral interactions.
    
    E0: binding energy at zero coverage [eV]
    theta: surface coverage
    epsilon: lateral interaction parameter [eV]
    model: 'linear', 'quadratic', or 'piecewise'
    """
    if model == 'linear':
        return E0 + epsilon * theta
    elif model == 'quadratic':
        return E0 + epsilon * theta + 0.5 * epsilon * theta**2
    elif model == 'piecewise':
        # Piecewise linear: different slopes in low/high coverage regimes
        theta_c = 0.5  # crossover coverage
        E = np.where(theta < theta_c,
                     E0 + epsilon * theta,
                     E0 + epsilon * theta_c + 2.0 * epsilon * (theta - theta_c))
        return E
    else:
        raise ValueError(f"Unknown model: {model}")
