"""
Module 6: Chemical Reaction Optimization Case Study.

Demonstrates the framework on a realistic chemical reaction optimization problem:
optimizing reaction yield and selectivity by tuning temperature, pressure,
catalyst loading, solvent ratio, residence time, and pH.
"""

import torch
import numpy as np
import json


# Reaction parameter space definition
REACTION_PARAMS = {
    "temperature_C": {"low": 50.0, "high": 200.0, "unit": "°C"},
    "pressure_bar": {"low": 1.0, "high": 50.0, "unit": "bar"},
    "catalyst_loading_mol_pct": {"low": 0.5, "high": 10.0, "unit": "mol%"},
    "solvent_ratio": {"low": 0.1, "high": 1.0, "unit": "v/v"},
    "residence_time_min": {"low": 5.0, "high": 120.0, "unit": "min"},
    "pH": {"low": 4.0, "high": 10.0, "unit": "-"},
}


def get_bounds():
    """Get normalized bounds [0,1]^d and parameter metadata."""
    param_names = list(REACTION_PARAMS.keys())
    dim = len(param_names)
    bounds = torch.stack([
        torch.zeros(dim, dtype=torch.double),
        torch.ones(dim, dtype=torch.double),
    ])
    return bounds, param_names


def denormalize(X_norm, param_names=None):
    """Convert normalized [0,1] values to physical units."""
    if param_names is None:
        param_names = list(REACTION_PARAMS.keys())

    X_phys = torch.zeros_like(X_norm)
    for i, name in enumerate(param_names):
        lo = REACTION_PARAMS[name]["low"]
        hi = REACTION_PARAMS[name]["high"]
        X_phys[:, i] = lo + X_norm[:, i] * (hi - lo)
    return X_phys


def simulate_reaction_yield(X_norm):
    """Simulate reaction yield (%) as a function of normalized parameters.

    Physics-inspired model with:
    - Arrhenius-like temperature dependence
    - Pressure effect on equilibrium
    - Catalyst loading with diminishing returns
    - Solvent effects
    - Residence time approach to equilibrium
    - pH sensitivity window
    """
    X = denormalize(X_norm)
    T = X[:, 0]      # temperature
    P = X[:, 1]      # pressure
    cat = X[:, 2]    # catalyst loading
    solv = X[:, 3]   # solvent ratio
    t_res = X[:, 4]  # residence time
    pH = X[:, 5]     # pH

    # Arrhenius-like (optimal around 140°C)
    T_opt = 140.0
    yield_T = torch.exp(-0.5 * ((T - T_opt) / 30.0) ** 2)

    # Pressure effect (higher is better, diminishing)
    yield_P = 1 - torch.exp(-P / 15.0)

    # Catalyst (diminishing returns, optimal ~5 mol%)
    yield_cat = cat / (cat + 2.0)

    # Solvent ratio effect
    yield_solv = torch.exp(-0.5 * ((solv - 0.6) / 0.2) ** 2)

    # Residence time (approach to equilibrium)
    yield_t = 1 - torch.exp(-t_res / 30.0)

    # pH window (optimal 6.5-7.5)
    yield_pH = torch.exp(-0.5 * ((pH - 7.0) / 1.0) ** 2)

    # Combined yield with interactions
    base_yield = (
        yield_T * yield_P * yield_cat * yield_solv * yield_t * yield_pH
    )

    # Add interaction terms
    interaction = 0.05 * torch.sin(T / 20.0) * torch.cos(P / 10.0)

    # Noise
    noise = 0.02 * torch.randn_like(base_yield)

    yield_pct = 95.0 * base_yield + interaction + noise
    return yield_pct.clamp(0, 100)


def simulate_reaction_selectivity(X_norm):
    """Simulate reaction selectivity (%) — anti-correlated with aggressive conditions."""
    X = denormalize(X_norm)
    T = X[:, 0]
    P = X[:, 1]
    cat = X[:, 2]
    pH = X[:, 5]

    # Selectivity decreases at extreme conditions
    sel_T = torch.exp(-0.5 * ((T - 100.0) / 40.0) ** 2)
    sel_P = torch.exp(-P / 80.0) + 0.5
    sel_cat = 1.0 / (1.0 + 0.1 * cat)
    sel_pH = torch.exp(-0.5 * ((pH - 7.5) / 1.5) ** 2)

    selectivity = 90.0 * sel_T * sel_P * sel_cat * sel_pH
    noise = 0.015 * torch.randn_like(selectivity)
    return (selectivity + noise).clamp(0, 100)


def reaction_objectives(X_norm):
    """Combined objectives: yield and selectivity (both to maximize)."""
    y = simulate_reaction_yield(X_norm)
    s = simulate_reaction_selectivity(X_norm)
    return torch.stack([y, s], dim=-1)


def run_single_objective_optimization(n_init=15, n_iter=40, seed=42):
    """Optimize yield only using standard BO."""
    from bayesopt_framework.acquisition_functions import run_bo_loop

    bounds, param_names = get_bounds()
    result = run_bo_loop(
        objective_fn=simulate_reaction_yield,
        bounds=bounds,
        acq_name="EI",
        n_init=n_init,
        n_iter=n_iter,
        seed=seed,
    )

    best_idx = result["train_Y"].argmax()
    best_X_norm = result["train_X"][best_idx].unsqueeze(0)
    best_X_phys = denormalize(best_X_norm)

    conditions = {}
    for i, name in enumerate(param_names):
        conditions[name] = {
            "value": round(best_X_phys[0, i].item(), 2),
            "unit": REACTION_PARAMS[name]["unit"],
        }

    return {
        "best_yield": round(result["final_best"], 2),
        "best_conditions": conditions,
        "convergence": result["best_values"],
        "n_experiments": n_init + n_iter,
    }


def run_multi_objective_optimization(n_init=15, n_iter=30, seed=42):
    """Optimize yield and selectivity using MOBO (q-EHVI)."""
    from bayesopt_framework.multi_objective import run_mobo_loop

    bounds, param_names = get_bounds()
    ref_point = [0.0, 0.0]  # reference point for HV

    result = run_mobo_loop(
        objective_fn=reaction_objectives,
        bounds=bounds,
        ref_point=ref_point,
        n_objectives=2,
        n_init=n_init,
        n_iter=n_iter,
        batch_size=1,
        seed=seed,
    )

    # Extract Pareto-optimal conditions
    pareto_conditions = []
    for i in range(result["pareto_X"].shape[0]):
        x_norm = result["pareto_X"][i].unsqueeze(0)
        x_phys = denormalize(x_norm)
        cond = {}
        for j, name in enumerate(param_names):
            cond[name] = round(x_phys[0, j].item(), 2)
        cond["yield"] = round(result["pareto_Y"][i, 0].item(), 2)
        cond["selectivity"] = round(result["pareto_Y"][i, 1].item(), 2)
        pareto_conditions.append(cond)

    # Sort by yield
    pareto_conditions.sort(key=lambda c: c["yield"], reverse=True)

    return {
        "final_hv": result["final_hv"],
        "n_pareto": result["n_pareto"],
        "hv_history": result["hv_history"],
        "pareto_front": pareto_conditions[:10],  # top 10
        "n_experiments": n_init + n_iter,
    }


def run_full_case_study(seed=42):
    """Run the complete chemical reaction optimization case study."""
    so_result = run_single_objective_optimization(n_init=15, n_iter=40, seed=seed)
    mo_result = run_multi_objective_optimization(n_init=15, n_iter=30, seed=seed)

    return {
        "single_objective": so_result,
        "multi_objective": mo_result,
        "parameter_space": {
            name: {
                "range": [p["low"], p["high"]],
                "unit": p["unit"],
            }
            for name, p in REACTION_PARAMS.items()
        },
    }
