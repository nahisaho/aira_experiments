"""
Differential Privacy (DP) integration for federated learning.

Implements:
  - Gaussian mechanism with (epsilon, delta)-DP
  - Gradient clipping (per-sample and per-client)
  - Privacy budget accounting (RDP -> (epsilon, delta) conversion)
  - Adaptive clipping (Andrew et al., 2021)
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
import math


@dataclass
class PrivacyBudget:
    """Track cumulative privacy expenditure."""
    epsilon: float = 0.0
    delta: float = 0.0
    total_epsilon: float = 10.0
    total_delta: float = 1e-5
    rounds_consumed: int = 0

    @property
    def remaining_epsilon(self) -> float:
        return max(0, self.total_epsilon - self.epsilon)

    @property
    def is_exhausted(self) -> bool:
        return self.epsilon >= self.total_epsilon


class GaussianMechanism:
    """
    Gaussian mechanism for (epsilon, delta)-differential privacy.
    sigma >= Delta_f * sqrt(2 * ln(1.25/delta)) / epsilon
    """

    def __init__(self, epsilon: float, delta: float, sensitivity: float):
        self.epsilon = epsilon
        self.delta = delta
        self.sensitivity = sensitivity
        self.sigma = self._compute_sigma()

    def _compute_sigma(self) -> float:
        return self.sensitivity * math.sqrt(2 * math.log(1.25 / self.delta)) / self.epsilon

    def add_noise(self, value: np.ndarray) -> np.ndarray:
        noise = np.random.normal(0, self.sigma, size=value.shape)
        return value + noise


class GradientClipper:
    """Per-client gradient clipping for DP-FL."""

    def __init__(self, max_norm: float = 1.0):
        self.max_norm = max_norm

    def clip(self, gradients: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        total_norm = 0.0
        for g in gradients.values():
            total_norm += np.sum(g ** 2)
        total_norm = np.sqrt(total_norm)

        clip_factor = min(1.0, self.max_norm / (total_norm + 1e-10))
        clipped = {k: v * clip_factor for k, v in gradients.items()}
        return clipped

    def compute_clip_fraction(self, all_norms: List[float]) -> float:
        clipped_count = sum(1 for n in all_norms if n > self.max_norm)
        return clipped_count / len(all_norms) if all_norms else 0.0


class AdaptiveClipper:
    """Adaptive clipping (Andrew et al., 2021)."""

    def __init__(
        self,
        initial_clip: float = 1.0,
        target_quantile: float = 0.5,
        clip_lr: float = 0.2,
    ):
        self.clip_bound = initial_clip
        self.target_quantile = target_quantile
        self.clip_lr = clip_lr

    def update(self, gradient_norms: List[float]) -> float:
        if not gradient_norms:
            return self.clip_bound
        below = sum(1 for n in gradient_norms if n <= self.clip_bound)
        current_quantile = below / len(gradient_norms)
        self.clip_bound *= math.exp(
            self.clip_lr * (current_quantile - self.target_quantile)
        )
        self.clip_bound = max(self.clip_bound, 0.01)
        return self.clip_bound


class RDPAccountant:
    """
    Renyi Differential Privacy (RDP) accountant.
    RDP of Gaussian mechanism at order alpha: alpha * Delta_f^2 / (2*sigma^2)
    """

    def __init__(self, total_epsilon: float = 10.0, delta: float = 1e-5):
        self.total_epsilon = total_epsilon
        self.delta = delta
        self.rdp_orders = list(range(2, 256))
        self.rdp_eps = np.zeros(len(self.rdp_orders))

    def accumulate(self, sigma: float, sensitivity: float = 1.0, q: float = 1.0):
        for i, alpha in enumerate(self.rdp_orders):
            if q == 1.0:
                rdp = alpha * (sensitivity ** 2) / (2 * sigma ** 2)
            else:
                base_rdp = alpha * (sensitivity ** 2) / (2 * sigma ** 2)
                log_bound = 2 * math.log(q) + base_rdp
                rdp = min(base_rdp, max(0, log_bound))
            self.rdp_eps[i] += rdp

    def get_privacy_spent(self) -> Tuple[float, float]:
        eps_list = []
        for i, alpha in enumerate(self.rdp_orders):
            eps = self.rdp_eps[i] - math.log(self.delta) / (alpha - 1)
            eps_list.append(eps)
        best_eps = min(eps_list) if eps_list else float("inf")
        return max(0, best_eps), self.delta

    def get_remaining_budget(self) -> float:
        spent_eps, _ = self.get_privacy_spent()
        return max(0, self.total_epsilon - spent_eps)

    def can_continue(self) -> bool:
        spent_eps, _ = self.get_privacy_spent()
        return spent_eps < self.total_epsilon


class DPFederatedAggregator:
    """Combines gradient clipping + Gaussian noise + RDP accounting."""

    def __init__(
        self,
        epsilon: float = 10.0,
        delta: float = 1e-5,
        max_grad_norm: float = 1.0,
        noise_multiplier: float = 1.0,
        num_clients: int = 10,
        clients_per_round: int = 5,
    ):
        self.clipper = GradientClipper(max_norm=max_grad_norm)
        self.accountant = RDPAccountant(total_epsilon=epsilon, delta=delta)
        self.noise_multiplier = noise_multiplier
        self.max_grad_norm = max_grad_norm
        self.num_clients = num_clients
        self.clients_per_round = clients_per_round
        self.sampling_rate = clients_per_round / num_clients

    def aggregate_with_dp(
        self,
        global_weights: Dict[str, np.ndarray],
        client_deltas: List[Dict[str, np.ndarray]],
    ) -> Dict[str, np.ndarray]:
        if not self.accountant.can_continue():
            # Return global weights unchanged when budget exhausted
            return {k: v.copy() for k, v in global_weights.items()}

        clipped_deltas = [self.clipper.clip(d) for d in client_deltas]

        K = len(clipped_deltas)
        avg_delta = {}
        for key in global_weights:
            avg_delta[key] = sum(d[key] for d in clipped_deltas) / K

        sigma = self.noise_multiplier * self.max_grad_norm / K
        noised_delta = {}
        for key, val in avg_delta.items():
            noise = np.random.normal(0, sigma, size=val.shape)
            noised_delta[key] = val + noise

        self.accountant.accumulate(
            sigma=self.noise_multiplier,
            sensitivity=self.max_grad_norm,
            q=self.sampling_rate,
        )

        new_weights = {}
        for key in global_weights:
            new_weights[key] = global_weights[key] + noised_delta[key]

        return new_weights

    def get_privacy_report(self) -> Dict[str, float]:
        eps, delta = self.accountant.get_privacy_spent()
        return {
            "epsilon_spent": eps,
            "delta": delta,
            "epsilon_remaining": self.accountant.get_remaining_budget(),
            "budget_utilization": eps / self.accountant.total_epsilon if self.accountant.total_epsilon > 0 else 0,
            "can_continue": self.accountant.can_continue(),
        }
