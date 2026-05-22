"""
Sensitivity Analysis for SOA Box Model
Uses:
  - Morris elementary effect screening (global)
  - Sobol variance decomposition
  - Brute-force one-at-a-time (OAT) for local sensitivity
"""
import numpy as np
from scipy.integrate import solve_ivp
from typing import Dict, List, Tuple
import logging

logger = logging.getLogger(__name__)


def run_boxmodel_scalar(params: Dict[str, float], t_end: float = 3600 * 6) -> float:
    """
    Run simplified box model and return final SOA [μg/m3].
    Used as objective function for sensitivity analysis.
    """
    import sys, os
    sys.path.insert(0, os.path.dirname(__file__))
    from box_model import SimplifiedSOABoxModel, AtmosphericConditions, VOCEmission

    cond = AtmosphericConditions(
        T     = params.get("T",      298.15),
        RH    = params.get("RH",     0.50),
        NOx_ppb = params.get("NOx_ppb", 5.0),
        O3_ppb  = params.get("O3_ppb",  30.0),
        JNO2  = params.get("JNO2",   8e-3),
    )
    voc = VOCEmission("alpha_pinene", params.get("VOC_ppb", 2.0))

    model = SimplifiedSOABoxModel(cond, voc)
    res   = model.run(t_end=t_end)
    return float(res["SOA"][-1])


# ── Parameter ranges for sensitivity analysis ─────────────────────────────────
PARAM_BOUNDS = {
    "T":        (278.0,  318.0),   # K
    "RH":       (0.20,   0.90),
    "NOx_ppb":  (1.0,    20.0),    # ppb
    "O3_ppb":   (10.0,   80.0),    # ppb
    "JNO2":     (1e-3,   2e-2),    # s-1
    "VOC_ppb":  (0.5,    10.0),    # ppb
}

PARAM_NAMES  = list(PARAM_BOUNDS.keys())
PARAM_LABELS = {
    "T":        "Temperature (K)",
    "RH":       "Relative Humidity",
    "NOx_ppb":  "NOx (ppb)",
    "O3_ppb":   "O₃ (ppb)",
    "JNO2":     "J(NO₂) (s⁻¹)",
    "VOC_ppb":  "VOC (ppb)",
}


def _normalize(val: float, lo: float, hi: float) -> float:
    return (val - lo) / (hi - lo)


def _denormalize(x: float, lo: float, hi: float) -> float:
    return lo + x * (hi - lo)


def oat_sensitivity(
    base_params: Dict[str, float],
    perturbation: float = 0.10,
    t_end: float = 3600 * 6,
) -> Dict[str, Dict[str, float]]:
    """One-At-a-Time local sensitivity analysis."""
    base_soa = run_boxmodel_scalar(base_params, t_end)
    results  = {}

    for param in PARAM_NAMES:
        lo, hi  = PARAM_BOUNDS[param]
        delta   = perturbation * (hi - lo)
        p_plus  = {**base_params, param: min(base_params[param] + delta, hi)}
        p_minus = {**base_params, param: max(base_params[param] - delta, lo)}

        soa_plus  = run_boxmodel_scalar(p_plus, t_end)
        soa_minus = run_boxmodel_scalar(p_minus, t_end)

        # Normalized sensitivity index S = (dSOA/SOA) / (dparam/param)
        d_soa    = soa_plus - soa_minus
        d_param  = p_plus[param] - p_minus[param]
        S_abs    = (d_soa / (2 * delta)) if delta > 0 else 0.0
        S_norm   = S_abs * base_params[param] / max(base_soa, 1e-6)

        results[param] = {
            "S_normalized":    S_norm,
            "S_absolute":      S_abs,
            "SOA_plus":        soa_plus,
            "SOA_minus":       soa_minus,
            "SOA_base":        base_soa,
            "delta_SOA":       d_soa,
            "relative_change": d_soa / max(base_soa, 1e-6),
        }

    return results


def morris_screening(
    n_trajectories: int = 20,
    n_levels: int = 4,
    t_end: float = 3600 * 6,
    seed: int = 42,
) -> Dict[str, Dict[str, float]]:
    """Morris elementary effect method for global sensitivity screening."""
    rng   = np.random.default_rng(seed)
    k     = len(PARAM_NAMES)
    delta = n_levels / (2 * (n_levels - 1))

    mu_star = np.zeros(k)     # mean |EE|  (non-monotonic sensitivity)
    sigma   = np.zeros(k)     # std EE      (interaction / non-linearity)
    all_EE  = [[] for _ in range(k)]

    for _ in range(n_trajectories):
        # Random base point
        x0 = rng.uniform(0, 1 - delta, size=k)
        x  = x0.copy()
        perm = rng.permutation(k)

        # Evaluate at base
        params_base = {
            PARAM_NAMES[i]: _denormalize(x[i], *PARAM_BOUNDS[PARAM_NAMES[i]])
            for i in range(k)
        }
        y_prev = run_boxmodel_scalar(params_base, t_end)

        for j in perm:
            x_new = x.copy()
            x_new[j] += delta
            x_new[j]  = min(x_new[j], 1.0)

            params_new = {
                PARAM_NAMES[i]: _denormalize(x_new[i], *PARAM_BOUNDS[PARAM_NAMES[i]])
                for i in range(k)
            }
            y_new = run_boxmodel_scalar(params_new, t_end)
            EE_j  = (y_new - y_prev) / delta

            all_EE[j].append(EE_j)
            y_prev = y_new
            x      = x_new

    for i in range(k):
        ee = np.array(all_EE[i])
        mu_star[i] = float(np.mean(np.abs(ee)))
        sigma[i]   = float(np.std(ee))

    return {
        PARAM_NAMES[i]: {
            "mu_star": mu_star[i],
            "sigma":   sigma[i],
            "ratio":   sigma[i] / max(mu_star[i], 1e-10),
        }
        for i in range(k)
    }


def sobol_first_order(
    n_samples: int = 64,
    t_end: float = 3600 * 6,
    seed: int = 42,
) -> Dict[str, float]:
    """
    Simplified Sobol first-order sensitivity indices via Saltelli sampling.
    S1_i ≈ Var[E(Y|X_i)] / Var[Y]
    """
    rng = np.random.default_rng(seed)
    k   = len(PARAM_NAMES)

    # Sample A and B matrices
    A = rng.uniform(0, 1, (n_samples, k))
    B = rng.uniform(0, 1, (n_samples, k))

    def eval_mat(mat: np.ndarray) -> np.ndarray:
        y = np.zeros(n_samples)
        for j in range(n_samples):
            params = {
                PARAM_NAMES[i]: _denormalize(mat[j, i], *PARAM_BOUNDS[PARAM_NAMES[i]])
                for i in range(k)
            }
            y[j] = run_boxmodel_scalar(params, t_end)
        return y

    yA = eval_mat(A)
    yB = eval_mat(B)

    Var_Y = np.var(np.concatenate([yA, yB]))
    if Var_Y < 1e-12:
        return {name: 0.0 for name in PARAM_NAMES}

    S1 = {}
    for i, name in enumerate(PARAM_NAMES):
        AB_i = A.copy()
        AB_i[:, i] = B[:, i]
        yAB_i = eval_mat(AB_i)
        # Jansen estimator
        S1[name] = float(max(0, (Var_Y - 0.5 * np.mean((yB - yAB_i) ** 2)) / Var_Y))

    return S1
