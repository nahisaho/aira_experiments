"""
Bayesian Optimization of Reaction Conditions
==============================================
Gaussian Process-based Bayesian optimization for continuous flow reaction parameters.
Optimizes: temperature, flow rate, concentration, catalyst loading.
Objective: maximize yield while minimizing impurities.
"""

import numpy as np
import json, os

np.random.seed(42)

# --- Parameter Space ---
PARAM_BOUNDS = {
    "temperature_C":     (50, 150),
    "flow_rate_mL_min":  (0.1, 5.0),
    "concentration_M":   (0.1, 2.0),
    "catalyst_mol_pct":  (1.0, 10.0),
}

class GaussianProcessSurrogate:
    """Simplified GP surrogate using RBF kernel for demonstration."""

    def __init__(self, length_scales=None, noise=0.01):
        self.X_train = None
        self.y_train = None
        self.noise = noise
        self.length_scales = length_scales or np.array([20.0, 1.0, 0.5, 2.0])
        self.K_inv = None

    def rbf_kernel(self, X1, X2):
        sqdist = np.sum(((X1[:, np.newaxis, :] - X2[np.newaxis, :, :]) /
                         self.length_scales) ** 2, axis=2)
        return np.exp(-0.5 * sqdist)

    def fit(self, X, y):
        self.X_train = np.array(X)
        self.y_train = np.array(y)
        K = self.rbf_kernel(self.X_train, self.X_train)
        K += self.noise * np.eye(len(K))
        self.K_inv = np.linalg.inv(K)

    def predict(self, X_new):
        X_new = np.array(X_new)
        if self.X_train is None:
            return np.zeros(len(X_new)), np.ones(len(X_new))
        K_s = self.rbf_kernel(X_new, self.X_train)
        K_ss = self.rbf_kernel(X_new, X_new)
        mu = K_s @ self.K_inv @ self.y_train
        cov = K_ss - K_s @ self.K_inv @ K_s.T
        sigma = np.sqrt(np.maximum(np.diag(cov), 1e-10))
        return mu, sigma

def reaction_model(temp, flow_rate, conc, cat_pct):
    """Synthetic reaction model: Arrhenius kinetics + mixing effects."""
    Ea = 60000  # J/mol
    R = 8.314
    T_K = temp + 273.15
    k = 1e8 * np.exp(-Ea / (R * T_K))

    tau = 30.0 / flow_rate  # residence time proxy
    Da = k * conc * tau     # Damköhler number

    conversion = Da / (1 + Da)           # CSTR-like
    selectivity = 0.95 - 0.005 * (temp - 80)**2 / 1000  # optimum ~80°C
    selectivity = np.clip(selectivity, 0.5, 0.99)

    cat_effect = 1 - np.exp(-0.5 * cat_pct)
    mixing_penalty = 1.0 if flow_rate < 3.0 else 0.95

    yield_pct = conversion * selectivity * cat_effect * mixing_penalty * 100
    yield_pct = np.clip(yield_pct, 0, 99)

    noise = np.random.normal(0, 1.5)
    return yield_pct + noise

def expected_improvement(mu, sigma, y_best, xi=0.01):
    """Expected Improvement acquisition function."""
    from scipy.stats import norm
    imp = mu - y_best - xi
    Z = imp / (sigma + 1e-10)
    ei = imp * norm.cdf(Z) + sigma * norm.pdf(Z)
    ei[sigma < 1e-10] = 0.0
    return ei

def latin_hypercube_sample(bounds, n_samples):
    """Generate Latin Hypercube samples."""
    d = len(bounds)
    result = np.zeros((n_samples, d))
    for i, (lo, hi) in enumerate(bounds.values()):
        cuts = np.linspace(lo, hi, n_samples + 1)
        for j in range(n_samples):
            result[j, i] = np.random.uniform(cuts[j], cuts[j + 1])
    np.random.shuffle(result)
    return result

def run_bayesian_optimization(n_initial=10, n_iterations=30):
    bounds = PARAM_BOUNDS
    param_names = list(bounds.keys())

    # Initial sampling (Latin Hypercube)
    X_init = latin_hypercube_sample(bounds, n_initial)
    y_init = np.array([reaction_model(*x) for x in X_init])

    X_all = X_init.tolist()
    y_all = y_init.tolist()

    gp = GaussianProcessSurrogate()
    optimization_history = []

    for it in range(n_iterations):
        gp.fit(np.array(X_all), np.array(y_all))
        y_best = max(y_all)

        # Generate candidates
        n_candidates = 5000
        candidates = np.zeros((n_candidates, len(bounds)))
        for i, (lo, hi) in enumerate(bounds.values()):
            candidates[:, i] = np.random.uniform(lo, hi, n_candidates)

        mu, sigma = gp.predict(candidates)
        ei = expected_improvement(mu, sigma, y_best)

        best_idx = np.argmax(ei)
        x_next = candidates[best_idx]

        y_next = reaction_model(*x_next)

        X_all.append(x_next.tolist())
        y_all.append(float(y_next))

        record = {
            "iteration": it + 1,
            "parameters": {name: round(float(x_next[i]), 3) for i, name in enumerate(param_names)},
            "predicted_yield": round(float(mu[best_idx]), 2),
            "uncertainty": round(float(sigma[best_idx]), 2),
            "ei_value": round(float(ei[best_idx]), 4),
            "observed_yield": round(float(y_next), 2),
            "best_yield_so_far": round(max(y_all), 2),
        }
        optimization_history.append(record)

    best_idx_all = int(np.argmax(y_all))
    best_params = {name: round(X_all[best_idx_all][i], 3) for i, name in enumerate(param_names)}

    results = {
        "optimization_config": {
            "n_initial_samples": n_initial,
            "n_bo_iterations": n_iterations,
            "total_experiments": n_initial + n_iterations,
            "acquisition_function": "Expected Improvement",
            "surrogate_model": "Gaussian Process (RBF kernel)",
        },
        "optimal_conditions": {
            "parameters": best_params,
            "best_yield_pct": round(max(y_all), 2),
        },
        "convergence": {
            "yield_after_10_experiments": round(max(y_all[:10]), 2),
            "yield_after_20_experiments": round(max(y_all[:20]), 2),
            "yield_after_30_experiments": round(max(y_all[:30]), 2),
            "yield_after_40_experiments": round(max(y_all[:min(40, len(y_all))]), 2),
        },
        "optimization_history": optimization_history,
    }

    return results

if __name__ == "__main__":
    from scipy.stats import norm  # ensure available
    results = run_bayesian_optimization()
    os.makedirs("results", exist_ok=True)

    with open("results/bayesian_optimization_results.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print("=== Bayesian Optimization Results ===")
    print(f"\nOptimal Conditions:")
    for k, v in results["optimal_conditions"]["parameters"].items():
        print(f"  {k}: {v}")
    print(f"\nBest Yield: {results['optimal_conditions']['best_yield_pct']:.2f}%")
    print(f"\nConvergence:")
    for k, v in results["convergence"].items():
        print(f"  {k}: {v}%")
