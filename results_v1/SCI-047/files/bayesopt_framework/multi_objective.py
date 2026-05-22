"""
Module 4: Multi-Objective Bayesian Optimization (Expected Hypervolume Improvement).

Implements EHVI and q-EHVI for simultaneous optimization of multiple objectives
(e.g., yield and selectivity in chemical reactions).
"""

import torch
import numpy as np
import time
from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_mll
from botorch.models.model_list_gp_regression import ModelListGP
from gpytorch.mlls import ExactMarginalLogLikelihood
from gpytorch.mlls.sum_marginal_log_likelihood import SumMarginalLogLikelihood
from botorch.acquisition.multi_objective import (
    qExpectedHypervolumeImprovement,
)
from botorch.utils.multi_objective.box_decompositions.non_dominated import (
    FastNondominatedPartitioning,
)
from botorch.utils.multi_objective.pareto import is_non_dominated
from botorch.optim import optimize_acqf
from botorch.sampling.normal import SobolQMCNormalSampler


def build_multi_objective_model(train_X, train_Y_list):
    """Build independent GP models for each objective."""
    models = []
    for Y in train_Y_list:
        model = SingleTaskGP(train_X, Y)
        models.append(model)
    model_list = ModelListGP(*models)
    return model_list


def fit_multi_objective_model(model_list):
    """Fit all GP models in the model list."""
    mll = SumMarginalLogLikelihood(model_list.likelihood, model_list)
    fit_gpytorch_mll(mll)
    return mll


def compute_hypervolume(Y, ref_point):
    """Compute the hypervolume indicator of Pareto-optimal points."""
    from botorch.utils.multi_objective.hypervolume import Hypervolume

    pareto_mask = is_non_dominated(Y)
    pareto_Y = Y[pareto_mask]

    if pareto_Y.shape[0] == 0:
        return 0.0

    hv = Hypervolume(ref_point=ref_point)
    return hv.compute(pareto_Y)


def run_mobo_loop(
    objective_fn,
    bounds,
    ref_point,
    n_objectives=2,
    n_init=10,
    n_iter=30,
    batch_size=1,
    seed=42,
):
    """Run multi-objective BO loop with q-EHVI."""
    torch.manual_seed(seed)
    dim = bounds.shape[1]

    # Initial samples
    train_X = bounds[0] + (bounds[1] - bounds[0]) * torch.rand(
        n_init, dim, dtype=torch.double
    )
    train_Y = objective_fn(train_X)  # shape: (n, n_objectives)

    ref_point_tensor = torch.tensor(ref_point, dtype=torch.double)
    hv_history = [compute_hypervolume(train_Y, ref_point_tensor)]

    for i in range(n_iter):
        # Build and fit model
        train_Y_list = [train_Y[:, j:j+1] for j in range(n_objectives)]
        model = build_multi_objective_model(train_X, train_Y_list)
        fit_multi_objective_model(model)

        # Partitioning for EHVI
        partitioning = FastNondominatedPartitioning(
            ref_point=ref_point_tensor, Y=train_Y,
        )

        sampler = SobolQMCNormalSampler(sample_shape=torch.Size([128]))
        acq_fn = qExpectedHypervolumeImprovement(
            model=model,
            ref_point=ref_point_tensor.tolist(),
            partitioning=partitioning,
            sampler=sampler,
        )

        candidates, _ = optimize_acqf(
            acq_function=acq_fn,
            bounds=bounds,
            q=batch_size,
            num_restarts=10,
            raw_samples=512,
        )

        new_Y = objective_fn(candidates)
        train_X = torch.cat([train_X, candidates])
        train_Y = torch.cat([train_Y, new_Y])

        hv_history.append(compute_hypervolume(train_Y, ref_point_tensor))

    pareto_mask = is_non_dominated(train_Y)

    return {
        "train_X": train_X,
        "train_Y": train_Y,
        "pareto_X": train_X[pareto_mask],
        "pareto_Y": train_Y[pareto_mask],
        "hv_history": hv_history,
        "final_hv": hv_history[-1],
        "n_pareto": pareto_mask.sum().item(),
    }


def synthetic_mo_objective(X):
    """Synthetic 2-objective function: competing objectives with trade-off."""
    # Objective 1: yield-like (maximize)
    f1 = 1.0 - torch.sum((X - 0.3) ** 2, dim=-1)
    # Objective 2: selectivity-like (maximize, anti-correlated with f1)
    f2 = 1.0 - torch.sum((X - 0.7) ** 2, dim=-1)
    # Add interaction
    f1 = f1 + 0.1 * torch.sin(5 * X[:, 0]) if X.dim() > 1 else f1
    f2 = f2 + 0.1 * torch.cos(5 * X[:, 0]) if X.dim() > 1 else f2
    return torch.stack([f1, f2], dim=-1)


def run_mobo_experiment(seed=42):
    """Run MOBO experiment with synthetic objectives."""
    dim = 4
    bounds = torch.stack([
        torch.zeros(dim, dtype=torch.double),
        torch.ones(dim, dtype=torch.double),
    ])
    ref_point = [-1.0, -1.0]

    results = run_mobo_loop(
        objective_fn=synthetic_mo_objective,
        bounds=bounds,
        ref_point=ref_point,
        n_objectives=2,
        n_init=15,
        n_iter=25,
        batch_size=1,
        seed=seed,
    )

    return {
        "final_hv": results["final_hv"],
        "n_pareto": results["n_pareto"],
        "hv_history": results["hv_history"],
        "pareto_Y": results["pareto_Y"].tolist(),
    }
