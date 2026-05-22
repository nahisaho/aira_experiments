"""
Federated aggregation strategies:
  - FedAvg: Weighted averaging of model parameters
  - FedProx: Proximal term for heterogeneous data
  - SCAFFOLD: Variance reduction via control variates
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
import copy


@dataclass
class ClientUpdate:
    """Container for a single client's training result."""
    client_id: str
    weights: Dict[str, np.ndarray]
    num_samples: int
    loss: float
    control_variate: Optional[Dict[str, np.ndarray]] = None
    metrics: Dict[str, float] = field(default_factory=dict)


class FedAvg:
    """
    Federated Averaging (McMahan et al., 2017).

    Convergence guarantee (convex case):
        E[f(w_T) - f(w*)] <= O(1 / (K * T))
    where K = number of clients, T = number of rounds.

    For non-convex with bounded gradient dissimilarity sigma^2:
        (1/T) sum E[||grad f(w_t)||^2] <= O(sigma / sqrt(K*T) + 1/T)
    """

    def __init__(self, learning_rate: float = 1.0):
        self.lr = learning_rate

    def aggregate(
        self,
        global_weights: Dict[str, np.ndarray],
        client_updates: List[ClientUpdate],
    ) -> Dict[str, np.ndarray]:
        total_samples = sum(u.num_samples for u in client_updates)
        new_weights = {}

        for key in global_weights:
            weighted_sum = np.zeros_like(global_weights[key])
            for update in client_updates:
                w = update.num_samples / total_samples
                weighted_sum += w * update.weights[key]
            new_weights[key] = (
                (1 - self.lr) * global_weights[key] + self.lr * weighted_sum
            )

        return new_weights

    def compute_convergence_bound(
        self,
        num_clients: int,
        num_rounds: int,
        gradient_dissimilarity: float,
        local_steps: int,
        learning_rate_local: float,
    ) -> Dict[str, float]:
        """Compute theoretical convergence bounds for FedAvg."""
        convex_bound = 1.0 / (num_clients * num_rounds)
        nonconvex_bound = (
            gradient_dissimilarity / np.sqrt(num_clients * num_rounds)
            + 1.0 / num_rounds
        )
        client_drift = (
            local_steps * learning_rate_local * gradient_dissimilarity
        )
        return {
            "convex_bound": convex_bound,
            "nonconvex_bound": nonconvex_bound,
            "client_drift_bound": client_drift,
            "effective_rate": 1.0 / np.sqrt(num_clients * num_rounds),
        }


class FedProx:
    """
    FedProx (Li et al., 2020): Handles data heterogeneity by adding
    a proximal term mu/2 * ||w - w_global||^2 to the local objective.

    The local objective becomes:
        h_k(w; w_t) = F_k(w) + mu/2 * ||w - w_t||^2

    Convergence: Under (B, sigma)-bounded dissimilarity,
        E[f(w_T)] - f* <= O(B^2/(mu*T))
    """

    def __init__(self, mu: float = 0.01, learning_rate: float = 1.0):
        self.mu = mu
        self.lr = learning_rate

    def compute_local_loss(
        self,
        task_loss: float,
        local_weights: Dict[str, np.ndarray],
        global_weights: Dict[str, np.ndarray],
    ) -> float:
        """Compute FedProx local loss = task_loss + mu/2 * ||w - w_global||^2."""
        proximal_term = 0.0
        for key in local_weights:
            diff = local_weights[key] - global_weights[key]
            proximal_term += np.sum(diff ** 2)
        return task_loss + (self.mu / 2.0) * proximal_term

    def aggregate(
        self,
        global_weights: Dict[str, np.ndarray],
        client_updates: List[ClientUpdate],
    ) -> Dict[str, np.ndarray]:
        total_samples = sum(u.num_samples for u in client_updates)
        new_weights = {}
        for key in global_weights:
            weighted_sum = np.zeros_like(global_weights[key])
            for update in client_updates:
                w = update.num_samples / total_samples
                weighted_sum += w * update.weights[key]
            new_weights[key] = weighted_sum
        return new_weights


class SCAFFOLD:
    """
    SCAFFOLD (Karimireddy et al., 2020): Variance reduction via control variates.

    Each client maintains a control variate c_i, and the server maintains c.
    Local update: w <- w - eta(grad F_i(w) - c_i + c)
    Server update: c <- c + (1/K) sum(c_i_new - c_i_old)

    Convergence: O(1/T) rate regardless of data heterogeneity,
    eliminating the client drift problem of FedAvg.
    """

    def __init__(self, num_clients: int, learning_rate: float = 1.0):
        self.num_clients = num_clients
        self.lr = learning_rate
        self.server_control: Optional[Dict[str, np.ndarray]] = None
        self.client_controls: Dict[str, Dict[str, np.ndarray]] = {}

    def initialize_controls(self, model_weights: Dict[str, np.ndarray]):
        self.server_control = {
            k: np.zeros_like(v) for k, v in model_weights.items()
        }
        for i in range(self.num_clients):
            self.client_controls[str(i)] = {
                k: np.zeros_like(v) for k, v in model_weights.items()
            }

    def compute_corrected_gradient(
        self,
        client_id: str,
        local_gradient: Dict[str, np.ndarray],
    ) -> Dict[str, np.ndarray]:
        corrected = {}
        for key in local_gradient:
            corrected[key] = (
                local_gradient[key]
                - self.client_controls[client_id][key]
                + self.server_control[key]
            )
        return corrected

    def aggregate(
        self,
        global_weights: Dict[str, np.ndarray],
        client_updates: List[ClientUpdate],
        local_steps: int,
        local_lr: float,
    ) -> Dict[str, np.ndarray]:
        K = len(client_updates)
        new_weights = {}

        for key in global_weights:
            delta_sum = np.zeros_like(global_weights[key])
            for update in client_updates:
                delta = update.weights[key] - global_weights[key]
                delta_sum += delta
            new_weights[key] = global_weights[key] + (self.lr / K) * delta_sum

        for update in client_updates:
            cid = update.client_id
            if update.control_variate is not None:
                for key in self.server_control:
                    old_ci = self.client_controls[cid][key].copy()
                    new_ci = update.control_variate[key]
                    self.server_control[key] += (new_ci - old_ci) / self.num_clients
                    self.client_controls[cid][key] = new_ci

        return new_weights
