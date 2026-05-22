"""
Module 2: Acquisition Function Comparison and Problem-Dependent Selection.

Compares EI (Expected Improvement), UCB (Upper Confidence Bound),
and KG (Knowledge Gradient) with problem-dependent selection criteria.
"""

import torch
import numpy as np
import time
from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_mll
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.acquisition import (
    ExpectedImprovement,
    UpperConfidenceBound,
    qKnowledgeGradient,
    qExpectedImprovement,
)
from botorch.optim import optimize_acqf
from botorch.test_functions import Hartmann, Branin


ACQUISITION_REGISTRY = {
    "EI": lambda model, best_f, **kw: ExpectedImprovement(model=model, best_f=best_f),
    "UCB": lambda model, best_f, **kw: UpperConfidenceBound(
        model=model, beta=kw.get("beta", 2.0)
    ),
    "KG": lambda model, best_f, **kw: qKnowledgeGradient(
        model=model, num_fantasies=kw.get("num_fantasies", 64)
    ),
}


SELECTION_GUIDELINES = {
    "EI": {
        "best_for": [
            "Low-noise problems",
            "Late-stage optimization (exploitation)",
            "Well-explored regions",
        ],
        "weaknesses": [
            "Can be too greedy early on",
            "Poor in high-noise settings",
        ],
        "computational_cost": "Low",
    },
    "UCB": {
        "best_for": [
            "Exploration-exploitation trade-off control via β",
            "Early-stage optimization",
            "Theoretical regret guarantees (GP-UCB)",
        ],
        "weaknesses": [
            "Requires β tuning",
            "Can over-explore with high β",
        ],
        "computational_cost": "Low",
    },
    "KG": {
        "best_for": [
            "High-noise problems",
            "Value of information reasoning",
            "Batch settings (look-ahead)",
        ],
        "weaknesses": [
            "High computational cost",
            "Many fantasies needed for accuracy",
        ],
        "computational_cost": "High",
    },
}


def run_bo_loop(
    objective_fn,
    bounds,
    acq_name="EI",
    n_init=10,
    n_iter=40,
    seed=42,
    acq_kwargs=None,
):
    """Run a single BO loop with the specified acquisition function."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    dim = bounds.shape[1]
    acq_kwargs = acq_kwargs or {}

    # Initial random samples
    train_X = bounds[0] + (bounds[1] - bounds[0]) * torch.rand(
        n_init, dim, dtype=torch.double
    )
    train_Y = objective_fn(train_X).unsqueeze(-1)

    best_values = [train_Y.max().item()]
    cumulative_times = [0.0]

    for i in range(n_iter):
        t0 = time.time()

        # Standardize Y
        Y_mean, Y_std = train_Y.mean(), train_Y.std().clamp(min=1e-6)
        Y_norm = (train_Y - Y_mean) / Y_std
        best_f_norm = Y_norm.max().item()

        model = SingleTaskGP(train_X, Y_norm)
        mll = ExactMarginalLogLikelihood(model.likelihood, model)
        fit_gpytorch_mll(mll)

        acq_factory = ACQUISITION_REGISTRY[acq_name]
        acq_fn = acq_factory(model, best_f_norm, **acq_kwargs)

        candidate, _ = optimize_acqf(
            acq_function=acq_fn,
            bounds=bounds,
            q=1,
            num_restarts=5,
            raw_samples=256,
        )

        new_Y = objective_fn(candidate).unsqueeze(-1)
        train_X = torch.cat([train_X, candidate])
        train_Y = torch.cat([train_Y, new_Y])

        best_values.append(train_Y.max().item())
        cumulative_times.append(cumulative_times[-1] + (time.time() - t0))

    return {
        "acq_name": acq_name,
        "best_values": best_values,
        "cumulative_times": cumulative_times,
        "final_best": best_values[-1],
        "train_X": train_X,
        "train_Y": train_Y,
    }


def compare_acquisition_functions(
    objective_fn,
    bounds,
    acq_names=None,
    n_init=10,
    n_iter=30,
    n_trials=3,
    seeds=None,
):
    """Compare multiple acquisition functions across multiple random seeds."""
    if acq_names is None:
        acq_names = ["EI", "UCB"]  # KG excluded by default (slow)
    if seeds is None:
        seeds = list(range(n_trials))

    all_results = {}
    for acq_name in acq_names:
        trial_results = []
        for s in seeds:
            res = run_bo_loop(
                objective_fn, bounds, acq_name=acq_name,
                n_init=n_init, n_iter=n_iter, seed=s,
            )
            trial_results.append(res["best_values"])
        all_results[acq_name] = {
            "trials": trial_results,
            "mean_best": np.mean([t[-1] for t in trial_results]),
            "std_best": np.std([t[-1] for t in trial_results]),
            "mean_trajectory": np.mean(trial_results, axis=0).tolist(),
            "std_trajectory": np.std(trial_results, axis=0).tolist(),
        }

    return all_results


def select_acquisition_function(noise_level, budget, dim, batch_size=1):
    """Problem-dependent acquisition function selection heuristic."""
    score = {"EI": 0, "UCB": 0, "KG": 0}

    # Noise level heuristic
    if noise_level < 0.05:
        score["EI"] += 2
        score["UCB"] += 1
    elif noise_level < 0.2:
        score["UCB"] += 2
        score["EI"] += 1
        score["KG"] += 1
    else:
        score["KG"] += 3
        score["UCB"] += 1

    # Budget heuristic
    if budget < 20:
        score["KG"] += 2  # Value of information critical
        score["EI"] += 1
    elif budget < 50:
        score["UCB"] += 1
        score["EI"] += 1
    else:
        score["EI"] += 2

    # Dimensionality heuristic
    if dim > 15:
        score["UCB"] += 1  # Simpler, faster
        score["KG"] -= 1  # Too expensive

    # Batch setting
    if batch_size > 1:
        score["KG"] += 2
        score["UCB"] += 1

    recommended = max(score, key=score.get)
    return {
        "recommended": recommended,
        "scores": score,
        "reasoning": _build_reasoning(noise_level, budget, dim, batch_size, recommended),
    }


def _build_reasoning(noise_level, budget, dim, batch_size, recommended):
    reasons = []
    if noise_level > 0.2:
        reasons.append(f"High noise (σ={noise_level:.2f}) favors KG for value-of-information reasoning")
    if budget < 20:
        reasons.append(f"Small budget (n={budget}) favors look-ahead methods")
    if dim > 15:
        reasons.append(f"High dimensionality (d={dim}) penalizes expensive acquisitions")
    if batch_size > 1:
        reasons.append(f"Batch setting (q={batch_size}) benefits from KG's fantasy mechanism")
    reasons.append(f"→ Recommended: {recommended}")
    return reasons


def run_acquisition_comparison(seed=42):
    """Run the acquisition function comparison on Hartmann-6."""
    neg_hartmann6 = Hartmann(dim=6, negate=True)
    bounds = torch.stack([
        torch.zeros(6, dtype=torch.double),
        torch.ones(6, dtype=torch.double),
    ])

    results = compare_acquisition_functions(
        objective_fn=neg_hartmann6,
        bounds=bounds,
        acq_names=["EI", "UCB"],
        n_init=10,
        n_iter=30,
        n_trials=3,
        seeds=[seed, seed + 1, seed + 2],
    )

    # Selection examples
    selection_examples = [
        select_acquisition_function(noise_level=0.01, budget=50, dim=6),
        select_acquisition_function(noise_level=0.3, budget=15, dim=6, batch_size=4),
        select_acquisition_function(noise_level=0.1, budget=30, dim=25),
    ]

    return {"comparison": results, "selection_examples": selection_examples}
