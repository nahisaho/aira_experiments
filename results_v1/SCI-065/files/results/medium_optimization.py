#!/usr/bin/env python3
"""Optimize time-programmed medium compositions for brain organoid culture."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.special import erf
import matplotlib.pyplot as plt

SEED = 9303
RNG = np.random.default_rng(SEED)
BASE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = BASE_DIR / "figures"
OUTPUT_CSV = RESULTS_DIR / "optimized_medium_profiles.csv"

PHASES: List[Tuple[str, int, int]] = [
    ("Phase 1: Neural induction", 0, 6),
    ("Phase 2: Neural patterning", 6, 25),
    ("Phase 3: Cortical differentiation", 25, 50),
    ("Phase 4: Maturation", 50, 90),
]

COMPONENTS = [
    ("bFGF", "ng/mL", 0.0, 40.0),
    ("EGF", "ng/mL", 0.0, 20.0),
    ("BDNF", "ng/mL", 0.0, 50.0),
    ("Retinoic acid", "uM", 0.0, 1.0),
    ("Matrigel", "%", 0.0, 2.0),
    ("Glucose", "mM", 4.5, 25.0),
    ("O2 tension", "%", 5.0, 21.0),
]
COMPONENT_NAMES = [name for name, *_ in COMPONENTS]
BOUNDS = np.array([(low, high) for _, _, low, high in COMPONENTS] * len(PHASES), dtype=float)
LOWER = BOUNDS[:, 0]
UPPER = BOUNDS[:, 1]
RANGES = UPPER - LOWER

IDEAL_PROFILE = np.array(
    [
        [12.0, 0.0, 0.0, 0.0, 0.0, 17.0, 5.5],
        [24.0, 12.0, 5.0, 0.12, 1.2, 16.0, 8.0],
        [8.0, 4.0, 28.0, 0.55, 0.8, 11.0, 12.0],
        [0.0, 0.0, 42.0, 0.20, 0.0, 8.0, 18.0],
    ],
    dtype=float,
)
WEIGHTS = np.array([1.2, 1.0, 1.4, 0.9, 0.6, 1.0, 1.1], dtype=float)
UNIT_COST = np.array([0.030, 0.022, 0.050, 9.0, 45.0, 0.45, 0.20], dtype=float)
PHASE_MEDIA_L = np.array([0.025, 0.055, 0.070, 0.100], dtype=float)
PHASE_IMPORTANCE = np.array([0.18, 0.24, 0.28, 0.30], dtype=float)
COST_CAP = 420.0


@dataclass
class Evaluation:
    """Container for evaluated designs."""

    params: np.ndarray
    score: float
    marker_expression: float
    electrophysiology: float
    viability: float
    cost_per_organoid: float
    cost_per_liter: float
    stability_penalty: float
    feasibility_penalty: float


def cost_per_liter(profile: np.ndarray) -> np.ndarray:
    """Return phase-specific medium cost in USD/L."""
    return profile @ UNIT_COST


def evaluate_profile(flat_params: np.ndarray) -> Evaluation:
    """Evaluate a candidate temporal medium program."""
    profile = flat_params.reshape(len(PHASES), len(COMPONENTS))
    normalized_error = ((profile - IDEAL_PROFILE) / np.array([c[3] - c[2] for c in COMPONENTS])) ** 2
    weighted_error = normalized_error * WEIGHTS
    phase_alignment = np.exp(-2.8 * weighted_error.sum(axis=1))

    marker_expression = float(np.sum(PHASE_IMPORTANCE * phase_alignment) / PHASE_IMPORTANCE.sum())
    late_phase = profile[2:, :]
    late_target = IDEAL_PROFILE[2:, :]
    late_error = np.sqrt(np.mean(((late_phase - late_target) / np.array([c[3] - c[2] for c in COMPONENTS])) ** 2))
    electrophysiology = float(np.exp(-3.0 * late_error))

    metabolic_penalty = 0.0
    if np.any(profile[:, 5] > 22.0):
        metabolic_penalty += 0.03 * np.sum(profile[:, 5] - 22.0)
    if np.any(profile[:, 6] < 6.0):
        metabolic_penalty += 0.06 * np.sum(6.0 - profile[:, 6])
    viability = float(np.clip(0.95 * marker_expression + 0.20 * electrophysiology - metabolic_penalty, 0.0, 1.0))

    phase_cost_liter = cost_per_liter(profile)
    total_cost_organoid = float(np.sum(phase_cost_liter * PHASE_MEDIA_L))
    mean_cost_liter = float(np.mean(phase_cost_liter))

    normalized_profile = (profile - LOWER.reshape(len(PHASES), len(COMPONENTS))) / RANGES.reshape(len(PHASES), len(COMPONENTS))
    stability_penalty = float(0.12 * np.mean(np.abs(np.diff(normalized_profile, axis=0))))
    feasibility_penalty = float(
        0.04 * np.mean(np.clip((profile[:, 0] + profile[:, 2]) / np.array([90.0] * len(PHASES)) - 0.8, 0.0, None))
        + 0.05 * np.mean(np.clip(profile[:, 4] - np.array([0.2, 1.6, 1.2, 0.2]), 0.0, None))
    )
    cost_penalty = float(np.clip((total_cost_organoid - COST_CAP) / COST_CAP, 0.0, None) * 0.45)

    composite_score = 100.0 * (0.45 * marker_expression + 0.30 * electrophysiology + 0.25 * viability)
    composite_score -= 100.0 * (stability_penalty + feasibility_penalty + cost_penalty)

    return Evaluation(
        params=flat_params.copy(),
        score=float(composite_score),
        marker_expression=100.0 * marker_expression,
        electrophysiology=100.0 * electrophysiology,
        viability=100.0 * viability,
        cost_per_organoid=total_cost_organoid,
        cost_per_liter=mean_cost_liter,
        stability_penalty=100.0 * stability_penalty,
        feasibility_penalty=100.0 * feasibility_penalty,
    )


def rbf_kernel(xa: np.ndarray, xb: np.ndarray, length_scale: float = 0.22, amplitude: float = 1.0) -> np.ndarray:
    """Squared-exponential kernel on normalized inputs."""
    xa = xa[:, None, :]
    xb = xb[None, :, :]
    sqdist = np.sum((xa - xb) ** 2, axis=2)
    return amplitude * np.exp(-0.5 * sqdist / (length_scale**2))


def gp_predict(x_train: np.ndarray, y_train: np.ndarray, x_query: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Gaussian-process posterior mean and standard deviation."""
    y_mean = y_train.mean()
    y_std = y_train.std() if y_train.std() > 1e-8 else 1.0
    y_scaled = (y_train - y_mean) / y_std

    k_xx = rbf_kernel(x_train, x_train) + 1e-6 * np.eye(len(x_train))
    chol = cho_factor(k_xx, lower=True, check_finite=False)
    alpha = cho_solve(chol, y_scaled, check_finite=False)

    k_xs = rbf_kernel(x_train, x_query)
    posterior_mean = k_xs.T @ alpha

    v = cho_solve(chol, k_xs, check_finite=False)
    k_ss_diag = np.ones(len(x_query))
    posterior_var = np.clip(k_ss_diag - np.sum(k_xs * v, axis=0), 1e-9, None)
    return posterior_mean * y_std + y_mean, np.sqrt(posterior_var) * y_std


def expected_improvement(mu: np.ndarray, sigma: np.ndarray, best_y: float, xi: float = 0.35) -> np.ndarray:
    """Expected-improvement acquisition function for maximization."""
    z = np.zeros_like(mu)
    valid = sigma > 1e-9
    z[valid] = (mu[valid] - best_y - xi) / sigma[valid]
    normal_pdf = np.exp(-0.5 * z**2) / np.sqrt(2.0 * np.pi)
    normal_cdf = 0.5 * (1.0 + erf(z / np.sqrt(2.0)))
    ei = (mu - best_y - xi) * normal_cdf + sigma * normal_pdf
    ei[~valid] = 0.0
    return np.maximum(ei, 0.0)


def bayesian_optimization(iterations: int = 42, initial_samples: int = 20) -> Tuple[Evaluation, List[Evaluation]]:
    """Run a simple Bayesian optimization routine over temporal medium settings."""
    sampled_x = RNG.uniform(LOWER, UPPER, size=(initial_samples, len(LOWER)))
    seeded_profile = IDEAL_PROFILE.flatten()
    sampled_x[0] = seeded_profile
    local_count = min(6, max(initial_samples - 1, 0))
    if local_count:
        sampled_x[1 : 1 + local_count] = np.clip(
            seeded_profile + RNG.normal(0.0, 0.05 * RANGES, size=(local_count, len(LOWER))),
            LOWER,
            UPPER,
        )
    history = [evaluate_profile(x) for x in sampled_x]

    for _ in range(iterations):
        x_train = np.array([item.params for item in history])
        y_train = np.array([item.score for item in history])
        x_train_norm = (x_train - LOWER) / RANGES

        candidates = RNG.uniform(LOWER, UPPER, size=(5000, len(LOWER)))
        elite = x_train[y_train.argmax()]
        local = np.clip(elite + RNG.normal(0.0, 0.08 * RANGES, size=(800, len(LOWER))), LOWER, UPPER)
        candidates[: len(local)] = local
        candidates_norm = (candidates - LOWER) / RANGES

        mu, sigma = gp_predict(x_train_norm, y_train, candidates_norm)
        acquisition = expected_improvement(mu, sigma, float(y_train.max()))
        next_x = candidates[np.argmax(acquisition)]
        history.append(evaluate_profile(next_x))

    best = max(history, key=lambda item: item.score)
    return best, history


def write_profile_csv(best: Evaluation) -> None:
    """Persist optimized temporal profiles to CSV."""
    profile = best.params.reshape(len(PHASES), len(COMPONENTS))
    phase_costs = cost_per_liter(profile)
    header = ["phase", "day_start", "day_end"] + [f"{name} ({unit})" for name, unit, _, _ in COMPONENTS] + [
        "cost_usd_per_liter",
        "phase_media_l_per_organoid",
    ]
    lines = [",".join(header)]
    for idx, (phase, start, end) in enumerate(PHASES):
        values = [phase, str(start), str(end)] + [f"{value:.4f}" for value in profile[idx]] + [
            f"{phase_costs[idx]:.2f}",
            f"{PHASE_MEDIA_L[idx]:.3f}",
        ]
        lines.append(",".join(values))
    OUTPUT_CSV.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_temporal_program(best: Evaluation) -> None:
    """Plot optimized temporal medium program."""
    profile = best.params.reshape(len(PHASES), len(COMPONENTS))
    phase_midpoints = np.array([(start + end) / 2 for _, start, end in PHASES])
    cmap = plt.get_cmap("cividis")
    fig, axes = plt.subplots(4, 2, figsize=(13, 12), constrained_layout=True)
    axes = axes.flatten()
    for idx, (name, unit, low, high) in enumerate(COMPONENTS):
        ax = axes[idx]
        color = cmap(0.15 + 0.75 * idx / (len(COMPONENTS) - 1))
        ax.plot(phase_midpoints, profile[:, idx], marker="o", linewidth=2.5, color=color)
        ax.set_title(f"{name} ({unit})")
        ax.set_xlabel("Culture day")
        ax.set_ylabel(f"Level ({unit})")
        ax.set_xlim(0, 90)
        ax.set_ylim(low - 0.05 * (high - low), high + 0.05 * (high - low))
        for _, start, end in PHASES:
            ax.axvspan(start, end, color="lightgrey", alpha=0.15)
        ax.grid(alpha=0.25)
    axes[-1].axis("off")
    axes[-1].text(
        0.02,
        0.70,
        "Optimized program\nComposite score: {:.1f}\nCost: ${:.1f}/organoid".format(best.score, best.cost_per_organoid),
        fontsize=13,
        weight="bold",
    )
    axes[-1].text(0.02, 0.28, "Shaded regions indicate culture phases.\nTargets balance maturation, viability, and cost.", fontsize=11)
    fig.suptitle("Temporal medium program for brain organoid differentiation", fontsize=15, weight="bold")
    fig.savefig(FIGURES_DIR / "medium_temporal_program.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_convergence(history: List[Evaluation]) -> None:
    """Plot optimization convergence."""
    scores = np.array([item.score for item in history])
    best_so_far = np.maximum.accumulate(scores)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(scores, color="#7a5195", alpha=0.5, label="Sampled score")
    ax.plot(best_so_far, color="#2f4b7c", linewidth=2.8, label="Best so far")
    ax.set_xlabel("Evaluation")
    ax.set_ylabel("Composite maturation score")
    ax.set_title("Bayesian optimization convergence")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(FIGURES_DIR / "optimization_convergence.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def pareto_front(costs: np.ndarray, qualities: np.ndarray) -> np.ndarray:
    """Return Pareto-efficient indices for low cost / high quality."""
    order = np.argsort(costs)
    best_quality = -np.inf
    keep = []
    for idx in order:
        if qualities[idx] > best_quality:
            keep.append(idx)
            best_quality = qualities[idx]
    return np.array(keep, dtype=int)


def plot_cost_tradeoff(history: List[Evaluation]) -> None:
    """Plot quality-versus-cost trade-off with Pareto front."""
    costs = np.array([item.cost_per_organoid for item in history])
    qualities = np.array([item.score for item in history])
    pareto_idx = pareto_front(costs, qualities)
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    scatter = ax.scatter(costs, qualities, c=qualities, cmap="viridis", s=45, alpha=0.8, edgecolor="none")
    ax.plot(costs[pareto_idx], qualities[pareto_idx], color="#ff7c43", linewidth=2.4, label="Pareto front")
    ax.set_xlabel("Cost per organoid (USD)")
    ax.set_ylabel("Composite maturation score")
    ax.set_title("Cost-quality trade-off for optimized medium programs")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.colorbar(scatter, ax=ax, label="Composite score")
    fig.savefig(FIGURES_DIR / "medium_cost_analysis.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    """Run optimization and save outputs."""
    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)
    best, history = bayesian_optimization()
    write_profile_csv(best)
    plot_temporal_program(best)
    plot_convergence(history)
    plot_cost_tradeoff(history)
    print(f"Best composite score: {best.score:.2f}")
    print(f"Cost per organoid: ${best.cost_per_organoid:.2f}")


if __name__ == "__main__":
    main()
