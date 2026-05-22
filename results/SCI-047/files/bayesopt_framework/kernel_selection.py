"""
Module 1: Gaussian Process Kernel Selection and Hyperparameter Optimization.

Provides systematic kernel comparison (RBF, Matérn, RQ, Spectral Mixture)
with automatic hyperparameter tuning via marginal likelihood maximization.
"""

import torch
import numpy as np
import json
import time
from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_mll
from gpytorch.mlls import ExactMarginalLogLikelihood
from gpytorch.kernels import (
    RBFKernel,
    MaternKernel,
    RQKernel,
    ScaleKernel,
    PeriodicKernel,
    AdditiveKernel,
    ProductKernel,
)
from gpytorch.means import ConstantMean
from gpytorch.distributions import MultivariateNormal
from sklearn.model_selection import KFold

KERNEL_REGISTRY = {
    "RBF": lambda dim: ScaleKernel(RBFKernel(ard_num_dims=dim)),
    "Matern52": lambda dim: ScaleKernel(MaternKernel(nu=2.5, ard_num_dims=dim)),
    "Matern32": lambda dim: ScaleKernel(MaternKernel(nu=1.5, ard_num_dims=dim)),
    "RQ": lambda dim: ScaleKernel(RQKernel(ard_num_dims=dim)),
    "RBF+Periodic": lambda dim: ScaleKernel(
        AdditiveKernel(RBFKernel(ard_num_dims=dim), PeriodicKernel())
    ),
}


def build_gp_model(train_X, train_Y, kernel_name="Matern52"):
    """Build a SingleTaskGP with the specified kernel."""
    dim = train_X.shape[-1]
    kernel_factory = KERNEL_REGISTRY.get(kernel_name)
    if kernel_factory is None:
        raise ValueError(f"Unknown kernel: {kernel_name}. Available: {list(KERNEL_REGISTRY.keys())}")

    covar_module = kernel_factory(dim)
    model = SingleTaskGP(train_X, train_Y, covar_module=covar_module)
    return model


def fit_model(model):
    """Fit GP model by maximizing marginal log-likelihood."""
    mll = ExactMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_mll(mll)
    return mll


def evaluate_kernel_cv(train_X, train_Y, kernel_name, n_folds=5, seed=42):
    """Evaluate a kernel using k-fold cross-validation (negative log predictive density)."""
    n = train_X.shape[0]
    if n < n_folds:
        n_folds = n

    kf = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    nlpds = []

    for train_idx, val_idx in kf.split(range(n)):
        X_tr = train_X[train_idx]
        Y_tr = train_Y[train_idx]
        X_val = train_X[val_idx]
        Y_val = train_Y[val_idx]

        try:
            model = build_gp_model(X_tr, Y_tr, kernel_name)
            fit_model(model)
            model.eval()

            with torch.no_grad():
                posterior = model.posterior(X_val)
                mean = posterior.mean
                var = posterior.variance
                # Negative log predictive density
                nlpd = 0.5 * torch.log(2 * torch.pi * var) + 0.5 * ((Y_val - mean) ** 2) / var
                nlpds.append(nlpd.mean().item())
        except Exception:
            nlpds.append(float("inf"))

    return {
        "kernel": kernel_name,
        "mean_nlpd": np.mean(nlpds),
        "std_nlpd": np.std(nlpds),
    }


def compare_kernels(train_X, train_Y, kernel_names=None, n_folds=5):
    """Compare multiple kernels via CV and return ranked results."""
    if kernel_names is None:
        kernel_names = list(KERNEL_REGISTRY.keys())

    results = []
    for kname in kernel_names:
        t0 = time.time()
        res = evaluate_kernel_cv(train_X, train_Y, kname, n_folds=n_folds)
        res["fit_time_s"] = time.time() - t0
        results.append(res)

    results.sort(key=lambda r: r["mean_nlpd"])
    return results


def extract_lengthscales(model):
    """Extract ARD lengthscales from a fitted GP model."""
    try:
        ls = model.covar_module.base_kernel.lengthscale.detach().cpu().numpy().flatten()
        return ls
    except AttributeError:
        return None


def run_kernel_comparison(n_samples=80, dim=6, seed=42):
    """Run a full kernel comparison experiment on synthetic data."""
    torch.manual_seed(seed)
    np.random.seed(seed)

    # Synthetic test function: modified Hartmann-like
    train_X = torch.rand(n_samples, dim, dtype=torch.double)
    # Non-linear function with varying sensitivity per dimension
    weights = torch.tensor([1.0, 0.5, 2.0, 0.1, 1.5, 0.3], dtype=torch.double)[:dim]
    train_Y = (
        torch.sin(3 * train_X @ weights.unsqueeze(1))
        + 0.1 * torch.randn(n_samples, 1, dtype=torch.double)
    )

    # Standardize
    Y_mean, Y_std = train_Y.mean(), train_Y.std()
    train_Y_norm = (train_Y - Y_mean) / Y_std

    results = compare_kernels(train_X, train_Y_norm)

    # Fit the best kernel to get lengthscales
    best_kernel = results[0]["kernel"]
    best_model = build_gp_model(train_X, train_Y_norm, best_kernel)
    fit_model(best_model)
    lengthscales = extract_lengthscales(best_model)

    return {
        "comparison": results,
        "best_kernel": best_kernel,
        "lengthscales": lengthscales.tolist() if lengthscales is not None else None,
        "n_samples": n_samples,
        "dim": dim,
    }
