"""
Module 3: Physics-Constrained Anomaly Scoring
Incorporates domain-specific physical laws into anomaly scores.
"""
import numpy as np
from dataclasses import dataclass, field
from typing import Callable, List, Optional


@dataclass
class PhysicalConstraint:
    """A single physical constraint with validation logic."""
    name: str
    check_fn: Callable[[np.ndarray], np.ndarray]
    weight: float = 1.0
    description: str = ""
    tolerance: float = 0.0


class PhysicsConstrainedScorer:
    """Anomaly scorer that combines statistical scores with physics-based penalties."""

    def __init__(self, constraints: List[PhysicalConstraint] = None,
                 stat_weight: float = 0.5, phys_weight: float = 0.5):
        self.constraints = constraints or []
        self.stat_weight = stat_weight
        self.phys_weight = phys_weight

    def add_constraint(self, name: str, check_fn: Callable, weight: float = 1.0,
                       description: str = "", tolerance: float = 0.0):
        self.constraints.append(PhysicalConstraint(
            name=name, check_fn=check_fn, weight=weight,
            description=description, tolerance=tolerance
        ))

    def score(self, data: np.ndarray, statistical_scores: np.ndarray) -> dict:
        n = len(data)
        stat_norm = self._normalize(statistical_scores)

        violation_matrix = np.zeros((n, len(self.constraints)))
        constraint_details = []

        for j, c in enumerate(self.constraints):
            violations = c.check_fn(data)
            violations = np.clip(np.abs(violations) - c.tolerance, 0, None)
            violation_matrix[:, j] = violations * c.weight
            n_violated = int(np.sum(violations > 0))
            constraint_details.append({
                "name": c.name,
                "description": c.description,
                "n_violations": n_violated,
                "violation_rate": n_violated / n,
                "max_violation": float(np.max(violations)),
                "mean_violation": float(np.mean(violations[violations > 0])) if n_violated > 0 else 0.0,
            })

        phys_scores = np.sum(violation_matrix, axis=1)
        phys_norm = self._normalize(phys_scores) if np.max(phys_scores) > 0 else phys_scores

        combined = self.stat_weight * stat_norm + self.phys_weight * phys_norm

        return {
            "combined_scores": combined,
            "statistical_scores": stat_norm,
            "physics_scores": phys_norm,
            "violation_matrix": violation_matrix,
            "constraint_details": constraint_details,
            "weights": {"statistical": self.stat_weight, "physics": self.phys_weight},
        }

    @staticmethod
    def _normalize(scores):
        s_min, s_max = np.min(scores), np.max(scores)
        if s_max - s_min < 1e-10:
            return np.zeros_like(scores)
        return (scores - s_min) / (s_max - s_min)


# ── Pre-built constraint library for common physics domains ──

def energy_conservation_constraint(energy_in_col=0, energy_out_col=1, tolerance=0.01):
    """Energy in ≈ Energy out (within tolerance fraction)."""
    def check(data):
        e_in = data[:, energy_in_col]
        e_out = data[:, energy_out_col]
        denom = np.maximum(np.abs(e_in), 1e-10)
        return np.abs(e_in - e_out) / denom - tolerance
    return PhysicalConstraint("energy_conservation", check, weight=2.0,
                              description="E_in ≈ E_out", tolerance=0.0)


def momentum_conservation_constraint(px_col=0, py_col=1, pz_col=2, tolerance=0.05):
    """Total momentum magnitude should be near zero (or conserved)."""
    def check(data):
        p_total = np.sqrt(data[:, px_col]**2 + data[:, py_col]**2 + data[:, pz_col]**2)
        return p_total - tolerance
    return PhysicalConstraint("momentum_conservation", check, weight=2.0,
                              description="|p_total| ≈ 0", tolerance=0.0)


def range_constraint(col, low, high, name=None):
    """Value must be within [low, high]."""
    def check(data):
        vals = data[:, col]
        below = np.clip(low - vals, 0, None)
        above = np.clip(vals - high, 0, None)
        return below + above
    return PhysicalConstraint(
        name or f"range_{col}", check, weight=1.0,
        description=f"col[{col}] ∈ [{low}, {high}]"
    )


def positive_definite_constraint(col, name=None):
    """Value must be strictly positive (e.g., mass, energy)."""
    def check(data):
        return np.clip(-data[:, col], 0, None)
    return PhysicalConstraint(
        name or f"positive_{col}", check, weight=1.5,
        description=f"col[{col}] > 0"
    )


def causality_constraint(time_col, speed_of_light=1.0, x_col=1, y_col=2, z_col=3):
    """Spacelike separation check: Δr ≤ c·Δt."""
    def check(data):
        if len(data) < 2:
            return np.zeros(len(data))
        dt = np.diff(data[:, time_col], prepend=data[0, time_col])
        dx = np.diff(data[:, x_col], prepend=data[0, x_col])
        dy = np.diff(data[:, y_col], prepend=data[0, y_col])
        dz = np.diff(data[:, z_col], prepend=data[0, z_col])
        dr = np.sqrt(dx**2 + dy**2 + dz**2)
        max_r = speed_of_light * np.abs(dt)
        return np.clip(dr - max_r, 0, None)
    return PhysicalConstraint("causality", check, weight=3.0,
                              description="Δr ≤ c·Δt")
