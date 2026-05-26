"""
Robust Design Module — Parameter uncertainty analysis and
optimization under stochastic variability.
"""

import numpy as np
from typing import Dict, List, Callable, Tuple, Optional
from src.stochastic_sim import StochasticModel, tau_leaping


def latin_hypercube_sample(
    param_ranges: Dict[str, Tuple[float, float]],
    n_samples: int,
    seed: int = 42,
) -> List[Dict[str, float]]:
    """Generate parameter samples using Latin Hypercube Sampling."""
    rng = np.random.RandomState(seed)
    params_list = list(param_ranges.keys())
    n_params = len(params_list)

    result = np.zeros((n_samples, n_params))
    for i in range(n_params):
        lo, hi = param_ranges[params_list[i]]
        intervals = np.linspace(lo, hi, n_samples + 1)
        points = np.array([rng.uniform(intervals[j], intervals[j+1])
                           for j in range(n_samples)])
        rng.shuffle(points)
        result[:, i] = points

    samples = []
    for j in range(n_samples):
        sample = {params_list[i]: result[j, i] for i in range(n_params)}
        samples.append(sample)
    return samples


def robustness_score(
    model_builder: Callable,
    param_ranges: Dict[str, Tuple[float, float]],
    objective_func: Callable,
    n_samples: int = 50,
    t_end: float = 500.0,
    tau: float = 0.5,
    seed: int = 42,
) -> Tuple[float, float, List[float]]:
    """
    Compute robustness score by evaluating the circuit under
    parameter uncertainty (LHS sampling).

    Returns (mean_score, std_score, all_scores).
    """
    samples = latin_hypercube_sample(param_ranges, n_samples, seed)
    scores = []

    for i, param_set in enumerate(samples):
        model = model_builder(param_set)
        times, states = tau_leaping(model, t_end, tau, seed=seed + i)
        score = objective_func(times, states, model.species)
        scores.append(score)

    return float(np.mean(scores)), float(np.std(scores)), scores


def toggle_switch_bistability_score(
    times: np.ndarray,
    states: np.ndarray,
    species: List[str],
) -> float:
    """
    Score toggle switch for bistability.
    High score = clear separation between LacI-high and TetR-high states.
    """
    lacI_idx = species.index("LacI")
    tetR_idx = species.index("TetR")

    # Use final 30% of trajectory
    n = len(times)
    start = int(0.7 * n)
    lacI_ss = states[start:, lacI_idx]
    tetR_ss = states[start:, tetR_idx]

    mean_lacI = np.mean(lacI_ss)
    mean_tetR = np.mean(tetR_ss)

    # Bistability metric: ratio of dominant species to minor
    ratio = abs(mean_lacI - mean_tetR) / (mean_lacI + mean_tetR + 1e-8)
    cv_lacI = np.std(lacI_ss) / (mean_lacI + 1e-8)
    cv_tetR = np.std(tetR_ss) / (mean_tetR + 1e-8)

    # High ratio (clear winner) + low CV (stable) = robust
    score = ratio * (1.0 / (1.0 + cv_lacI + cv_tetR))
    return float(score)


def repressilator_oscillation_score(
    times: np.ndarray,
    states: np.ndarray,
    species: List[str],
) -> float:
    """
    Score repressilator for sustained oscillation quality.
    Uses autocorrelation to detect periodicity.
    """
    lacI_idx = species.index("LacI")
    protein = states[:, lacI_idx]

    if len(protein) < 20:
        return 0.0

    # Detrend
    protein = protein - np.mean(protein)
    if np.std(protein) < 1e-6:
        return 0.0

    # Autocorrelation
    n = len(protein)
    acf = np.correlate(protein, protein, mode='full')
    acf = acf[n-1:]  # positive lags
    acf = acf / (acf[0] + 1e-12)

    # Find first peak after zero crossing
    peaks = []
    for i in range(2, len(acf) - 1):
        if acf[i] > acf[i-1] and acf[i] > acf[i+1] and acf[i] > 0.1:
            peaks.append((i, acf[i]))
            break

    if not peaks:
        return 0.0

    # Score = peak height (strength of oscillation)
    return float(peaks[0][1])


def optimize_circuit_params(
    model_builder: Callable,
    param_ranges: Dict[str, Tuple[float, float]],
    objective_func: Callable,
    n_iterations: int = 30,
    n_samples_per_iter: int = 20,
    t_end: float = 500.0,
    tau: float = 0.5,
    seed: int = 42,
) -> Tuple[Dict[str, float], float, List[float]]:
    """
    Simple evolutionary optimization of circuit parameters
    for robust performance under uncertainty.
    """
    rng = np.random.RandomState(seed)
    params_list = list(param_ranges.keys())

    # Initialize population
    population = latin_hypercube_sample(param_ranges, n_samples_per_iter, seed)
    best_params = None
    best_score = -np.inf
    history = []

    for gen in range(n_iterations):
        scores = []
        for i, params in enumerate(population):
            model = model_builder(params)
            try:
                times, states = tau_leaping(model, t_end, tau, seed=seed + gen * 100 + i)
                score = objective_func(times, states, model.species)
            except Exception:
                score = 0.0
            scores.append(score)

        gen_best_idx = np.argmax(scores)
        gen_best_score = scores[gen_best_idx]
        history.append(gen_best_score)

        if gen_best_score > best_score:
            best_score = gen_best_score
            best_params = population[gen_best_idx].copy()

        # Selection + mutation for next generation
        sorted_idx = np.argsort(scores)[::-1]
        elite = [population[sorted_idx[i]] for i in range(min(5, len(sorted_idx)))]

        new_pop = list(elite)
        while len(new_pop) < n_samples_per_iter:
            parent = elite[rng.randint(len(elite))].copy()
            child = {}
            for p in params_list:
                lo, hi = param_ranges[p]
                mutation = rng.normal(0, 0.1 * (hi - lo))
                child[p] = np.clip(parent[p] + mutation, lo, hi)
            new_pop.append(child)

        population = new_pop

    return best_params, best_score, history
