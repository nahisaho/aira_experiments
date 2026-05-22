"""
Module 3: Batch Bayesian Optimization (Parallel Experiment Proposals).

Implements q-batch strategies: q-EI, q-UCB, q-KG, and Kriging Believer
for proposing multiple experiments simultaneously.
"""

import torch
import numpy as np
import time
from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_mll
from gpytorch.mlls import ExactMarginalLogLikelihood
from botorch.acquisition import (
    qExpectedImprovement,
    qUpperConfidenceBound,
    qKnowledgeGradient,
)
from botorch.optim import optimize_acqf
from botorch.sampling.normal import SobolQMCNormalSampler
from botorch.test_functions import Hartmann


def propose_batch(model, bounds, batch_size=4, method="qEI", best_f=None, **kwargs):
    """Propose a batch of candidates using the specified batch acquisition method."""
    sampler = SobolQMCNormalSampler(sample_shape=torch.Size([512]))

    if method == "qEI":
        acq_fn = qExpectedImprovement(
            model=model, best_f=best_f, sampler=sampler
        )
    elif method == "qUCB":
        acq_fn = qUpperConfidenceBound(
            model=model, beta=kwargs.get("beta", 2.0), sampler=sampler
        )
    elif method == "qKG":
        acq_fn = qKnowledgeGradient(
            model=model, num_fantasies=kwargs.get("num_fantasies", 32)
        )
    elif method == "kriging_believer":
        return _kriging_believer(model, bounds, batch_size, best_f)
    else:
        raise ValueError(f"Unknown batch method: {method}")

    candidates, acq_value = optimize_acqf(
        acq_function=acq_fn,
        bounds=bounds,
        q=batch_size,
        num_restarts=10,
        raw_samples=512,
    )
    return candidates, acq_value.item()


def _kriging_believer(model, bounds, batch_size, best_f):
    """Kriging Believer: sequential greedy batch construction using GP mean as fantasy.

    Rebuilds the GP at each step to avoid batch-shape issues with get_fantasy_model.
    """
    from botorch.acquisition import ExpectedImprovement

    candidates = []
    # Accumulate training data
    X_all = model.train_inputs[0].clone()
    Y_all = model.train_targets.clone().unsqueeze(-1) if model.train_targets.dim() == 1 else model.train_targets.clone()

    for _ in range(batch_size):
        current_model = SingleTaskGP(X_all, Y_all)
        mll = ExactMarginalLogLikelihood(current_model.likelihood, current_model)
        fit_gpytorch_mll(mll)

        cur_best = Y_all.max().item()
        acq_fn = ExpectedImprovement(model=current_model, best_f=cur_best)
        candidate, _ = optimize_acqf(
            acq_function=acq_fn,
            bounds=bounds,
            q=1,
            num_restarts=5,
            raw_samples=256,
        )
        candidates.append(candidate)

        # Add candidate with GP mean prediction as fantasy observation
        with torch.no_grad():
            fantasy_Y = current_model.posterior(candidate).mean
        X_all = torch.cat([X_all, candidate])
        Y_all = torch.cat([Y_all, fantasy_Y])

    return torch.cat(candidates, dim=0), 0.0


def run_batch_bo_loop(
    objective_fn, bounds, batch_size=4, method="qEI",
    n_init=10, n_batches=10, seed=42,
):
    """Run a batch BO loop."""
    torch.manual_seed(seed)
    dim = bounds.shape[1]

    train_X = bounds[0] + (bounds[1] - bounds[0]) * torch.rand(
        n_init, dim, dtype=torch.double
    )
    train_Y = objective_fn(train_X).unsqueeze(-1)

    best_values = [train_Y.max().item()]
    batch_times = []
    n_evals = [n_init]

    for i in range(n_batches):
        t0 = time.time()

        Y_mean, Y_std = train_Y.mean(), train_Y.std().clamp(min=1e-6)
        Y_norm = (train_Y - Y_mean) / Y_std

        model = SingleTaskGP(train_X, Y_norm)
        mll = ExactMarginalLogLikelihood(model.likelihood, model)
        fit_gpytorch_mll(mll)

        best_f_norm = Y_norm.max().item()
        candidates, _ = propose_batch(
            model, bounds, batch_size=batch_size,
            method=method, best_f=best_f_norm,
        )

        new_Y = objective_fn(candidates).unsqueeze(-1)
        train_X = torch.cat([train_X, candidates])
        train_Y = torch.cat([train_Y, new_Y])

        best_values.append(train_Y.max().item())
        batch_times.append(time.time() - t0)
        n_evals.append(train_X.shape[0])

    return {
        "method": method,
        "batch_size": batch_size,
        "best_values": best_values,
        "batch_times": batch_times,
        "n_evals": n_evals,
        "final_best": best_values[-1],
    }


def compare_batch_methods(seed=42):
    """Compare batch optimization methods on Hartmann-6."""
    neg_hartmann6 = Hartmann(dim=6, negate=True)
    bounds = torch.stack([
        torch.zeros(6, dtype=torch.double),
        torch.ones(6, dtype=torch.double),
    ])

    methods = ["qEI", "qUCB", "kriging_believer"]
    results = {}

    for method in methods:
        res = run_batch_bo_loop(
            neg_hartmann6, bounds, batch_size=4,
            method=method, n_init=10, n_batches=10, seed=seed,
        )
        results[method] = {
            "final_best": res["final_best"],
            "best_values": res["best_values"],
            "n_evals": res["n_evals"],
            "mean_batch_time": np.mean(res["batch_times"]),
            "total_time": np.sum(res["batch_times"]),
        }

    return results
