#!/usr/bin/env python3
"""
Bayesian Optimization Framework for High-Dimensional Parameter Spaces
=====================================================================
Experiments covering:
1. GP kernel selection & hyperparameter optimization
2. Acquisition function comparison (EI, UCB, KG)
3. Batch (parallel) optimization
4. Multi-objective BO (EHVI)
5. High-dimensional optimization with REMBO-style dimensionality reduction
6. Chemical reaction optimization case study
"""

import warnings
warnings.filterwarnings("ignore")

import os
import time
import numpy as np
import pandas as pd
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

torch.manual_seed(42)
np.random.seed(42)

FIGURES_DIR = "figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
dtype = torch.double

# ============================================================
# Experiment 1: GP Kernel Selection & Hyperparameter Optimization
# ============================================================
print("=" * 60)
print("Experiment 1: GP Kernel Comparison")
print("=" * 60)

from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_mll
from gpytorch.mlls import ExactMarginalLogLikelihood
from gpytorch.kernels import RBFKernel, MaternKernel, ScaleKernel
from botorch.test_functions import Branin, Hartmann

def evaluate_kernel(train_X, train_Y, test_X, test_Y, kernel_class, kernel_name, **kwargs):
    covar_module = ScaleKernel(kernel_class(ard_num_dims=train_X.shape[-1], **kwargs))
    model = SingleTaskGP(train_X, train_Y, covar_module=covar_module).to(device)
    mll = ExactMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_mll(mll)
    model.eval()
    with torch.no_grad():
        posterior = model.posterior(test_X)
        pred_mean = posterior.mean
        pred_var = posterior.variance
    mse = ((pred_mean - test_Y) ** 2).mean().item()
    nll = -torch.distributions.Normal(pred_mean, pred_var.sqrt()).log_prob(test_Y).mean().item()
    lengthscales = model.covar_module.base_kernel.lengthscale.detach().cpu().numpy().flatten()
    return mse, nll, lengthscales

branin = Branin(negate=True)
bounds_branin = torch.tensor([[-5.0, 0.0], [10.0, 15.0]], device=device, dtype=dtype)

n_train, n_test = 50, 200
train_X = torch.rand(n_train, 2, device=device, dtype=dtype)
train_X = bounds_branin[0] + (bounds_branin[1] - bounds_branin[0]) * train_X
train_Y = branin(train_X).unsqueeze(-1)

test_X = torch.rand(n_test, 2, device=device, dtype=dtype)
test_X = bounds_branin[0] + (bounds_branin[1] - bounds_branin[0]) * test_X
test_Y = branin(test_X).unsqueeze(-1)

kernels = {
    "RBF": (RBFKernel, {}),
    "Matérn-5/2": (MaternKernel, {"nu": 2.5}),
    "Matérn-3/2": (MaternKernel, {"nu": 1.5}),
    "Matérn-1/2": (MaternKernel, {"nu": 0.5}),
}

kernel_results = {}
for name, (kclass, kw) in kernels.items():
    mse, nll, ls = evaluate_kernel(train_X, train_Y, test_X, test_Y, kclass, name, **kw)
    kernel_results[name] = {"MSE": mse, "NLL": nll, "Lengthscales": ls}
    print(f"  {name:12s}: MSE={mse:.4f}, NLL={nll:.4f}, LS={ls}")

# Plot kernel comparison
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
names = list(kernel_results.keys())
mses = [kernel_results[n]["MSE"] for n in names]
nlls = [kernel_results[n]["NLL"] for n in names]

colors = ["#2196F3", "#4CAF50", "#FF9800", "#F44336"]
axes[0].bar(names, mses, color=colors, edgecolor="black", linewidth=0.5)
axes[0].set_ylabel("Mean Squared Error")
axes[0].set_title("(a) Prediction MSE by Kernel")
axes[0].grid(axis="y", alpha=0.3)

axes[1].bar(names, nlls, color=colors, edgecolor="black", linewidth=0.5)
axes[1].set_ylabel("Negative Log-Likelihood")
axes[1].set_title("(b) NLL by Kernel")
axes[1].grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/kernel_comparison.png", dpi=150, bbox_inches="tight")
plt.close()
print("  -> Saved kernel_comparison.png")

# Lengthscale visualization
fig, ax = plt.subplots(figsize=(8, 5))
x_pos = np.arange(len(names))
width = 0.35
for i, name in enumerate(names):
    ls = kernel_results[name]["Lengthscales"]
    ax.bar(x_pos[i] - width/2, ls[0], width, color=colors[i], alpha=0.8, label=f"{name} dim1" if i == 0 else "")
    ax.bar(x_pos[i] + width/2, ls[1], width, color=colors[i], alpha=0.5)
ax.set_xticks(x_pos)
ax.set_xticklabels(names)
ax.set_ylabel("Lengthscale")
ax.set_title("Learned ARD Lengthscales by Kernel (dim1=dark, dim2=light)")
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/ard_lengthscales.png", dpi=150, bbox_inches="tight")
plt.close()
print("  -> Saved ard_lengthscales.png")


# ============================================================
# Experiment 2: Acquisition Function Comparison (EI, UCB, KG)
# ============================================================
print("\n" + "=" * 60)
print("Experiment 2: Acquisition Function Comparison")
print("=" * 60)

from botorch.acquisition import (
    ExpectedImprovement,
    UpperConfidenceBound,
    qKnowledgeGradient,
)
from botorch.optim import optimize_acqf

hartmann6 = Hartmann(dim=6, negate=True)
bounds_h6 = torch.zeros(2, 6, device=device, dtype=dtype)
bounds_h6[1] = 1.0

def run_bo_loop(acq_class, acq_kwargs, n_init=10, n_iter=40, label=""):
    train_X = torch.rand(n_init, 6, device=device, dtype=dtype)
    train_Y = hartmann6(train_X).unsqueeze(-1)
    best_values = [train_Y.max().item()]
    
    for i in range(n_iter):
        model = SingleTaskGP(train_X, train_Y).to(device)
        mll = ExactMarginalLogLikelihood(model.likelihood, model)
        fit_gpytorch_mll(mll)
        
        if acq_class == qKnowledgeGradient:
            acq = acq_class(model, num_fantasies=8, **acq_kwargs)
        elif acq_class == ExpectedImprovement:
            acq = acq_class(model, best_f=train_Y.max(), **acq_kwargs)
        else:
            acq = acq_class(model, **acq_kwargs)
        
        candidate, _ = optimize_acqf(
            acq, bounds=bounds_h6, q=1, num_restarts=5, raw_samples=64
        )
        new_Y = hartmann6(candidate).unsqueeze(-1)
        train_X = torch.cat([train_X, candidate])
        train_Y = torch.cat([train_Y, new_Y])
        best_values.append(train_Y.max().item())
    
    return best_values

acq_configs = {
    "EI": (ExpectedImprovement, {}),
    "UCB (β=2)": (UpperConfidenceBound, {"beta": 2.0}),
    "UCB (β=0.5)": (UpperConfidenceBound, {"beta": 0.5}),
}

n_repeats = 3
acq_results = {}

for name, (cls, kwargs) in acq_configs.items():
    print(f"  Running {name}...")
    all_runs = []
    for r in range(n_repeats):
        torch.manual_seed(42 + r)
        vals = run_bo_loop(cls, kwargs, label=name)
        all_runs.append(vals)
    acq_results[name] = np.array(all_runs)

# Also run KG (fewer iters due to cost)
print("  Running KG (reduced iterations)...")
kg_runs = []
for r in range(n_repeats):
    torch.manual_seed(42 + r)
    vals = run_bo_loop(qKnowledgeGradient, {}, n_iter=20, label="KG")
    kg_runs.append(vals)
acq_results["KG"] = np.array(kg_runs)

# Plot acquisition function comparison
fig, ax = plt.subplots(figsize=(10, 6))
colors_acq = {"EI": "#2196F3", "UCB (β=2)": "#4CAF50", "UCB (β=0.5)": "#FF9800", "KG": "#9C27B0"}

for name, data in acq_results.items():
    mean = data.mean(axis=0)
    std = data.std(axis=0)
    iters = np.arange(len(mean))
    ax.plot(iters, mean, label=name, color=colors_acq[name], linewidth=2)
    ax.fill_between(iters, mean - std, mean + std, alpha=0.15, color=colors_acq[name])

ax.axhline(y=3.32237, color="red", linestyle="--", alpha=0.7, label="Global optimum")
ax.set_xlabel("Iteration")
ax.set_ylabel("Best Value Found")
ax.set_title("Acquisition Function Comparison on Hartmann-6")
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/acquisition_comparison.png", dpi=150, bbox_inches="tight")
plt.close()
print("  -> Saved acquisition_comparison.png")


# ============================================================
# Experiment 3: Batch (Parallel) Optimization
# ============================================================
print("\n" + "=" * 60)
print("Experiment 3: Batch Optimization")
print("=" * 60)

from botorch.acquisition import qExpectedImprovement, qUpperConfidenceBound

def run_batch_bo(batch_size, n_init=10, n_rounds=15):
    train_X = torch.rand(n_init, 6, device=device, dtype=dtype)
    train_Y = hartmann6(train_X).unsqueeze(-1)
    best_values = [train_Y.max().item()]
    wall_times = [0.0]
    
    for i in range(n_rounds):
        model = SingleTaskGP(train_X, train_Y).to(device)
        mll = ExactMarginalLogLikelihood(model.likelihood, model)
        fit_gpytorch_mll(mll)
        
        acq = qExpectedImprovement(model, best_f=train_Y.max())
        t0 = time.time()
        candidates, _ = optimize_acqf(
            acq, bounds=bounds_h6, q=batch_size,
            num_restarts=5, raw_samples=64
        )
        wall_time = time.time() - t0
        
        new_Y = hartmann6(candidates).unsqueeze(-1)
        train_X = torch.cat([train_X, candidates])
        train_Y = torch.cat([train_Y, new_Y])
        best_values.append(train_Y.max().item())
        wall_times.append(wall_time)
    
    total_evals = n_init + n_rounds * batch_size
    return best_values, total_evals, np.sum(wall_times)

batch_sizes = [1, 2, 4, 8]
batch_results = {}
for bs in batch_sizes:
    print(f"  Batch size = {bs}...")
    runs = []
    for r in range(3):
        torch.manual_seed(42 + r)
        bv, te, wt = run_batch_bo(bs)
        runs.append({"best_values": bv, "total_evals": te, "wall_time": wt})
    batch_results[bs] = runs

# Plot batch results
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
batch_colors = {1: "#2196F3", 2: "#4CAF50", 4: "#FF9800", 8: "#F44336"}

for bs in batch_sizes:
    runs = batch_results[bs]
    all_bv = np.array([r["best_values"] for r in runs])
    mean_bv = all_bv.mean(axis=0)
    std_bv = all_bv.std(axis=0)
    rounds = np.arange(len(mean_bv))
    axes[0].plot(rounds, mean_bv, label=f"q={bs}", color=batch_colors[bs], linewidth=2)
    axes[0].fill_between(rounds, mean_bv - std_bv, mean_bv + std_bv, alpha=0.15, color=batch_colors[bs])

axes[0].axhline(y=3.32237, color="red", linestyle="--", alpha=0.7, label="Optimum")
axes[0].set_xlabel("BO Round")
axes[0].set_ylabel("Best Value Found")
axes[0].set_title("(a) Convergence by Batch Size")
axes[0].legend()
axes[0].grid(alpha=0.3)

# Efficiency plot
mean_evals = [np.mean([r["total_evals"] for r in batch_results[bs]]) for bs in batch_sizes]
mean_best = [np.mean([r["best_values"][-1] for r in batch_results[bs]]) for bs in batch_sizes]
mean_wt = [np.mean([r["wall_time"] for r in batch_results[bs]]) for bs in batch_sizes]

axes[1].bar([str(b) for b in batch_sizes], mean_wt, color=[batch_colors[b] for b in batch_sizes],
            edgecolor="black", linewidth=0.5)
axes[1].set_xlabel("Batch Size (q)")
axes[1].set_ylabel("Total Wall Time (s)")
axes[1].set_title("(b) Optimization Wall Time")
axes[1].grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/batch_optimization.png", dpi=150, bbox_inches="tight")
plt.close()
print("  -> Saved batch_optimization.png")


# ============================================================
# Experiment 4: Multi-Objective BO (EHVI)
# ============================================================
print("\n" + "=" * 60)
print("Experiment 4: Multi-Objective Bayesian Optimization (EHVI)")
print("=" * 60)

from botorch.test_functions.multi_objective import BraninCurrin
from botorch.models import SingleTaskGP
from botorch.utils.multi_objective.box_decompositions.non_dominated import (
    FastNondominatedPartitioning,
)
from botorch.acquisition.multi_objective.monte_carlo import (
    qExpectedHypervolumeImprovement,
)
from botorch.utils.multi_objective.pareto import is_non_dominated
from botorch.utils.sampling import draw_sobol_samples
from botorch.models.transforms.outcome import Standardize

branin_currin = BraninCurrin(negate=True)
bc_bounds = torch.zeros(2, 2, device=device, dtype=dtype)
bc_bounds[1] = 1.0
ref_point = torch.tensor([-18.0, -6.0], device=device, dtype=dtype)

n_init_mo = 10
n_iter_mo = 30
train_X_mo = draw_sobol_samples(bounds=bc_bounds, n=1, q=n_init_mo).squeeze(0).to(device, dtype)
train_Y_mo = branin_currin(train_X_mo)

hv_values = []

from botorch.utils.multi_objective.hypervolume import Hypervolume
hv_calc = Hypervolume(ref_point=ref_point)

pareto_mask = is_non_dominated(train_Y_mo)
pareto_Y = train_Y_mo[pareto_mask]
hv_values.append(hv_calc.compute(pareto_Y))

for i in range(n_iter_mo):
    model = SingleTaskGP(train_X_mo, train_Y_mo, outcome_transform=Standardize(m=2)).to(device)
    mll = ExactMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_mll(mll)
    
    partitioning = FastNondominatedPartitioning(ref_point=ref_point, Y=train_Y_mo)
    acq = qExpectedHypervolumeImprovement(
        model=model, ref_point=ref_point,
        partitioning=partitioning,
        sampler=None,
    )
    
    candidate, _ = optimize_acqf(
        acq, bounds=bc_bounds, q=1,
        num_restarts=5, raw_samples=64,
    )
    
    new_Y = branin_currin(candidate)
    train_X_mo = torch.cat([train_X_mo, candidate])
    train_Y_mo = torch.cat([train_Y_mo, new_Y])
    
    pareto_mask = is_non_dominated(train_Y_mo)
    pareto_Y = train_Y_mo[pareto_mask]
    hv_values.append(hv_calc.compute(pareto_Y))

print(f"  Final Hypervolume: {hv_values[-1]:.4f}")
print(f"  # Pareto points: {pareto_mask.sum().item()}")

# Plot MOBO results
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Pareto front
all_Y = train_Y_mo.cpu().numpy()
pareto_mask_np = pareto_mask.cpu().numpy()
axes[0].scatter(all_Y[~pareto_mask_np, 0], all_Y[~pareto_mask_np, 1],
                c="#BBDEFB", alpha=0.5, label="Dominated", edgecolor="gray", s=30)
axes[0].scatter(all_Y[pareto_mask_np, 0], all_Y[pareto_mask_np, 1],
                c="#F44336", s=80, zorder=5, label="Pareto Front", edgecolor="black")
sorted_pareto = all_Y[pareto_mask_np]
sorted_pareto = sorted_pareto[sorted_pareto[:, 0].argsort()]
axes[0].plot(sorted_pareto[:, 0], sorted_pareto[:, 1], "r--", alpha=0.5)
axes[0].set_xlabel("Objective 1 (neg. Branin)")
axes[0].set_ylabel("Objective 2 (neg. Currin)")
axes[0].set_title("(a) Pareto Front Discovery")
axes[0].legend()
axes[0].grid(alpha=0.3)

# Hypervolume convergence
axes[1].plot(range(len(hv_values)), hv_values, "b-o", markersize=3, linewidth=2)
axes[1].set_xlabel("Iteration")
axes[1].set_ylabel("Hypervolume")
axes[1].set_title("(b) Hypervolume Convergence")
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/mobo_ehvi.png", dpi=150, bbox_inches="tight")
plt.close()
print("  -> Saved mobo_ehvi.png")


# ============================================================
# Experiment 5: High-Dimensional BO with Dimensionality Reduction
# ============================================================
print("\n" + "=" * 60)
print("Experiment 5: High-Dimensional BO (REMBO-style)")
print("=" * 60)

def embedded_hartmann(x_high, projection_matrix, target_func, original_bounds):
    """Project high-dim x to low-dim via random matrix, then evaluate."""
    x_low = x_high @ projection_matrix
    x_low = torch.sigmoid(x_low)  # map to [0,1]
    return target_func(x_low)

D_high = 50
d_low = 6
A = torch.randn(D_high, d_low, device=device, dtype=dtype) / np.sqrt(D_high)
bounds_high = torch.zeros(2, D_high, device=device, dtype=dtype)
bounds_high[0] = -1.0
bounds_high[1] = 1.0

def run_rembo(n_init=20, n_iter=50):
    train_X = -1 + 2 * torch.rand(n_init, D_high, device=device, dtype=dtype)
    train_Y = embedded_hartmann(train_X, A, hartmann6, None).unsqueeze(-1)
    best_values = [train_Y.max().item()]
    
    for i in range(n_iter):
        model = SingleTaskGP(train_X, train_Y).to(device)
        mll = ExactMarginalLogLikelihood(model.likelihood, model)
        fit_gpytorch_mll(mll)
        
        acq = ExpectedImprovement(model, best_f=train_Y.max())
        candidate, _ = optimize_acqf(
            acq, bounds=bounds_high, q=1,
            num_restarts=5, raw_samples=128,
        )
        new_Y = embedded_hartmann(candidate, A, hartmann6, None).unsqueeze(-1)
        train_X = torch.cat([train_X, candidate])
        train_Y = torch.cat([train_Y, new_Y])
        best_values.append(train_Y.max().item())
    return best_values

def run_vanilla_high_dim(n_init=20, n_iter=50):
    train_X = torch.rand(n_init, D_high, device=device, dtype=dtype)
    train_Y = torch.tensor(
        [hartmann6(x[:6]).item() for x in train_X],
        device=device, dtype=dtype
    ).unsqueeze(-1)
    best_values = [train_Y.max().item()]
    
    for i in range(n_iter):
        model = SingleTaskGP(train_X, train_Y).to(device)
        mll = ExactMarginalLogLikelihood(model.likelihood, model)
        fit_gpytorch_mll(mll)
        
        acq = ExpectedImprovement(model, best_f=train_Y.max())
        candidate, _ = optimize_acqf(
            acq, bounds=torch.cat([torch.zeros(1, D_high, device=device, dtype=dtype),
                                    torch.ones(1, D_high, device=device, dtype=dtype)]),
            q=1, num_restarts=3, raw_samples=64,
        )
        new_Y = hartmann6(candidate[0, :6]).unsqueeze(-1).unsqueeze(-1)
        train_X = torch.cat([train_X, candidate])
        train_Y = torch.cat([train_Y, new_Y])
        best_values.append(train_Y.max().item())
    return best_values

print("  Running REMBO approach...")
rembo_runs = []
for r in range(3):
    torch.manual_seed(42 + r)
    rembo_runs.append(run_rembo(n_iter=40))

print("  Running Vanilla high-dim BO...")
vanilla_runs = []
for r in range(3):
    torch.manual_seed(42 + r)
    vanilla_runs.append(run_vanilla_high_dim(n_iter=40))

# Random search baseline
print("  Running Random Search baseline...")
random_runs = []
for r in range(3):
    torch.manual_seed(42 + r)
    vals = []
    best = -float("inf")
    for i in range(60):
        x = torch.rand(1, 6, device=device, dtype=dtype)
        y = hartmann6(x).item()
        best = max(best, y)
        vals.append(best)
    random_runs.append(vals)

fig, ax = plt.subplots(figsize=(10, 6))
rembo_arr = np.array(rembo_runs)
vanilla_arr = np.array(vanilla_runs)
random_arr = np.array(random_runs)

for data, label, color in [
    (rembo_arr, "REMBO (D=50→d=6)", "#2196F3"),
    (vanilla_arr, "Vanilla BO (D=50)", "#FF9800"),
    (random_arr, "Random Search", "#9E9E9E"),
]:
    mean = data.mean(axis=0)
    std = data.std(axis=0)
    ax.plot(range(len(mean)), mean, label=label, color=color, linewidth=2)
    ax.fill_between(range(len(mean)), mean - std, mean + std, alpha=0.15, color=color)

ax.axhline(y=3.32237, color="red", linestyle="--", alpha=0.7, label="Global optimum")
ax.set_xlabel("Iteration")
ax.set_ylabel("Best Value Found")
ax.set_title("High-Dimensional BO: REMBO vs Vanilla vs Random (D=50)")
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/high_dim_comparison.png", dpi=150, bbox_inches="tight")
plt.close()
print("  -> Saved high_dim_comparison.png")


# ============================================================
# Experiment 6: Chemical Reaction Optimization Case Study
# ============================================================
print("\n" + "=" * 60)
print("Experiment 6: Chemical Reaction Optimization")
print("=" * 60)

def chemical_reaction_yield(params):
    """Simulated chemical reaction: yield as function of
    temperature (50-200°C), pressure (1-10 atm), catalyst_loading (0.01-0.1 mol%),
    solvent_ratio (0-1), residence_time (1-60 min).
    Returns yield (%) and selectivity (%).
    """
    temp, pressure, cat_load, solv_ratio, res_time = params.unbind(-1)
    
    # Yield model (nonlinear, multi-modal)
    yield_val = (
        60 + 15 * torch.sin(0.05 * (temp - 120))
        + 8 * torch.log1p(pressure)
        - 200 * (cat_load - 0.05) ** 2
        + 10 * solv_ratio * (1 - solv_ratio)
        + 5 * torch.sin(0.1 * res_time)
        - 0.001 * (temp - 150) ** 2
        + 3 * torch.cos(0.3 * pressure * cat_load * 100)
    )
    yield_val = torch.clamp(yield_val, 0, 100)
    
    # Selectivity model
    selectivity = (
        85 - 0.1 * (temp - 130) ** 2 / 100
        + 5 * torch.exp(-((pressure - 5) ** 2) / 8)
        + 50 * cat_load
        - 10 * solv_ratio ** 2
        + 2 * torch.sin(0.15 * res_time)
    )
    selectivity = torch.clamp(selectivity, 0, 100)
    
    return torch.stack([yield_val, selectivity], dim=-1)

# Parameter bounds: temp, pressure, cat_load, solv_ratio, res_time
chem_bounds = torch.tensor(
    [[50.0, 1.0, 0.01, 0.0, 1.0],
     [200.0, 10.0, 0.10, 1.0, 60.0]],
    device=device, dtype=dtype
)
chem_ref_point = torch.tensor([0.0, 0.0], device=device, dtype=dtype)

n_init_chem = 15
n_iter_chem = 40

train_X_chem = torch.rand(n_init_chem, 5, device=device, dtype=dtype)
train_X_chem = chem_bounds[0] + (chem_bounds[1] - chem_bounds[0]) * train_X_chem
train_Y_chem = chemical_reaction_yield(train_X_chem)

hv_chem = []
hv_calc_chem = Hypervolume(ref_point=chem_ref_point)

pareto_mask_c = is_non_dominated(train_Y_chem)
hv_chem.append(hv_calc_chem.compute(train_Y_chem[pareto_mask_c]))

yield_history = [train_Y_chem[:, 0].max().item()]
selectivity_history = [train_Y_chem[:, 1].max().item()]

for i in range(n_iter_chem):
    model_chem = SingleTaskGP(
        train_X_chem, train_Y_chem,
        outcome_transform=Standardize(m=2)
    ).to(device)
    mll_chem = ExactMarginalLogLikelihood(model_chem.likelihood, model_chem)
    fit_gpytorch_mll(mll_chem)
    
    partitioning_chem = FastNondominatedPartitioning(
        ref_point=chem_ref_point, Y=train_Y_chem
    )
    acq_chem = qExpectedHypervolumeImprovement(
        model=model_chem, ref_point=chem_ref_point,
        partitioning=partitioning_chem,
    )
    
    candidate_chem, _ = optimize_acqf(
        acq_chem, bounds=chem_bounds, q=1,
        num_restarts=10, raw_samples=128,
    )
    
    new_Y_chem = chemical_reaction_yield(candidate_chem)
    train_X_chem = torch.cat([train_X_chem, candidate_chem])
    train_Y_chem = torch.cat([train_Y_chem, new_Y_chem])
    
    pareto_mask_c = is_non_dominated(train_Y_chem)
    hv_chem.append(hv_calc_chem.compute(train_Y_chem[pareto_mask_c]))
    yield_history.append(train_Y_chem[:, 0].max().item())
    selectivity_history.append(train_Y_chem[:, 1].max().item())

# Find best conditions
best_yield_idx = train_Y_chem[:, 0].argmax()
best_sel_idx = train_Y_chem[:, 1].argmax()
param_names = ["Temperature (°C)", "Pressure (atm)", "Catalyst (mol%)", "Solvent Ratio", "Residence Time (min)"]

print("\n  Best Yield Conditions:")
for j, pn in enumerate(param_names):
    print(f"    {pn}: {train_X_chem[best_yield_idx, j].item():.2f}")
print(f"    Yield: {train_Y_chem[best_yield_idx, 0].item():.1f}%, Selectivity: {train_Y_chem[best_yield_idx, 1].item():.1f}%")

print("\n  Best Selectivity Conditions:")
for j, pn in enumerate(param_names):
    print(f"    {pn}: {train_X_chem[best_sel_idx, j].item():.2f}")
print(f"    Yield: {train_Y_chem[best_sel_idx, 0].item():.1f}%, Selectivity: {train_Y_chem[best_sel_idx, 1].item():.1f}%")

# Plot chemical reaction results
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Pareto front
all_Y_chem = train_Y_chem.cpu().numpy()
pareto_np = pareto_mask_c.cpu().numpy()
axes[0, 0].scatter(all_Y_chem[~pareto_np, 0], all_Y_chem[~pareto_np, 1],
                    c="#BBDEFB", alpha=0.5, label="Dominated", edgecolor="gray", s=30)
axes[0, 0].scatter(all_Y_chem[pareto_np, 0], all_Y_chem[pareto_np, 1],
                    c="#F44336", s=80, zorder=5, label="Pareto Front", edgecolor="black")
axes[0, 0].set_xlabel("Yield (%)")
axes[0, 0].set_ylabel("Selectivity (%)")
axes[0, 0].set_title("(a) Yield vs Selectivity Pareto Front")
axes[0, 0].legend()
axes[0, 0].grid(alpha=0.3)

# Hypervolume convergence
axes[0, 1].plot(range(len(hv_chem)), hv_chem, "b-o", markersize=3, linewidth=2)
axes[0, 1].set_xlabel("Iteration")
axes[0, 1].set_ylabel("Hypervolume")
axes[0, 1].set_title("(b) Hypervolume Convergence")
axes[0, 1].grid(alpha=0.3)

# Yield & selectivity over time
axes[1, 0].plot(range(len(yield_history)), yield_history, "g-", linewidth=2, label="Best Yield")
axes[1, 0].plot(range(len(selectivity_history)), selectivity_history, "b-", linewidth=2, label="Best Selectivity")
axes[1, 0].set_xlabel("Iteration")
axes[1, 0].set_ylabel("Value (%)")
axes[1, 0].set_title("(c) Best Yield & Selectivity Over Iterations")
axes[1, 0].legend()
axes[1, 0].grid(alpha=0.3)

# Parameter sensitivity (violin-like: pareto vs non-pareto)
pareto_X = train_X_chem[pareto_mask_c].cpu().numpy()
norm_bounds = chem_bounds.cpu().numpy()
pareto_X_norm = (pareto_X - norm_bounds[0]) / (norm_bounds[1] - norm_bounds[0])
bp = axes[1, 1].boxplot(pareto_X_norm, labels=["Temp", "Press", "Cat", "Solv", "Time"],
                         patch_artist=True)
colors_box = ["#2196F3", "#4CAF50", "#FF9800", "#9C27B0", "#F44336"]
for patch, color in zip(bp["boxes"], colors_box):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
axes[1, 1].set_ylabel("Normalized Value")
axes[1, 1].set_title("(d) Pareto-Optimal Parameter Distributions")
axes[1, 1].grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig(f"{FIGURES_DIR}/chemical_optimization.png", dpi=150, bbox_inches="tight")
plt.close()
print("  -> Saved chemical_optimization.png")


# ============================================================
# Summary Table
# ============================================================
print("\n" + "=" * 60)
print("Generating Summary Tables")
print("=" * 60)

# Kernel comparison table
kernel_df = pd.DataFrame({
    "Kernel": list(kernel_results.keys()),
    "MSE": [kernel_results[k]["MSE"] for k in kernel_results],
    "NLL": [kernel_results[k]["NLL"] for k in kernel_results],
})
kernel_df.to_csv("figures/kernel_results.csv", index=False)

# Acquisition function table
acq_df_rows = []
for name, data in acq_results.items():
    final = data[:, -1]
    acq_df_rows.append({
        "Method": name,
        "Final Best (mean)": f"{final.mean():.4f}",
        "Final Best (std)": f"{final.std():.4f}",
        "Iterations": data.shape[1] - 1,
    })
acq_df = pd.DataFrame(acq_df_rows)
acq_df.to_csv("figures/acquisition_results.csv", index=False)

# Batch results table
batch_df_rows = []
for bs in batch_sizes:
    runs = batch_results[bs]
    finals = [r["best_values"][-1] for r in runs]
    wts = [r["wall_time"] for r in runs]
    batch_df_rows.append({
        "Batch Size": bs,
        "Final Best (mean)": f"{np.mean(finals):.4f}",
        "Total Evals": runs[0]["total_evals"],
        "Wall Time (s)": f"{np.mean(wts):.2f}",
    })
batch_df = pd.DataFrame(batch_df_rows)
batch_df.to_csv("figures/batch_results.csv", index=False)

print("  Saved CSV summary tables.")

# Save final numerical results for report
results_summary = {
    "kernel_best": min(kernel_results.items(), key=lambda x: x[1]["MSE"]),
    "acq_best_ei_final": acq_results["EI"][:, -1].mean(),
    "acq_best_ucb2_final": acq_results["UCB (β=2)"][:, -1].mean(),
    "mobo_final_hv": hv_values[-1],
    "mobo_pareto_count": pareto_mask.sum().item(),
    "chem_best_yield": train_Y_chem[best_yield_idx, 0].item(),
    "chem_best_selectivity": train_Y_chem[best_sel_idx, 1].item(),
    "chem_final_hv": hv_chem[-1],
    "rembo_final": rembo_arr[:, -1].mean(),
    "vanilla_final": vanilla_arr[:, -1].mean(),
}

import json
with open("figures/results_summary.json", "w") as f:
    json.dump({k: float(v) if isinstance(v, (int, float, np.floating)) else str(v) for k, v in results_summary.items()}, f, indent=2)

print("\nAll experiments completed successfully!")
print(f"Figures saved to: {FIGURES_DIR}/")
