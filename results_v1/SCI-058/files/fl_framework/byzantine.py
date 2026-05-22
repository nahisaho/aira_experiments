"""
Byzantine-resilient aggregation for federated learning.

Implements:
  - Krum / Multi-Krum (Blanchard et al., 2017)
  - Coordinate-wise median (Yin et al., 2018)
  - Trimmed mean (Yin et al., 2018)
  - FLTrust (Cao et al., 2021)
  - Byzantine anomaly detector
"""

import numpy as np
from typing import Dict, List, Optional, Tuple


def flatten_weights(weights: Dict[str, np.ndarray]) -> np.ndarray:
    return np.concatenate([v.flatten() for v in weights.values()])


def unflatten_weights(
    flat: np.ndarray, template: Dict[str, np.ndarray]
) -> Dict[str, np.ndarray]:
    result = {}
    offset = 0
    for key, val in template.items():
        size = val.size
        result[key] = flat[offset : offset + size].reshape(val.shape)
        offset += size
    return result


class Krum:
    """
    Krum aggregation (Blanchard et al., 2017).
    Tolerates f Byzantine clients where f < (n - 2) / 2.
    """

    def __init__(self, num_byzantine: int = 0, multi_krum_k: int = 1):
        self.num_byzantine = num_byzantine
        self.multi_krum_k = multi_krum_k

    def aggregate(
        self,
        global_weights: Dict[str, np.ndarray],
        client_updates: List[Dict[str, np.ndarray]],
    ) -> Dict[str, np.ndarray]:
        n = len(client_updates)
        f = self.num_byzantine
        n_minus_f_minus_2 = n - f - 2

        if n_minus_f_minus_2 < 1:
            raise ValueError(
                f"Not enough clients ({n}) for {f} Byzantine tolerance."
            )

        flat_updates = [flatten_weights(u) for u in client_updates]

        distances = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                d = np.linalg.norm(flat_updates[i] - flat_updates[j])
                distances[i, j] = d
                distances[j, i] = d

        scores = np.zeros(n)
        for i in range(n):
            sorted_dists = np.sort(distances[i])
            scores[i] = np.sum(sorted_dists[1 : n_minus_f_minus_2 + 1])

        selected_indices = np.argsort(scores)[: self.multi_krum_k]
        selected = [flat_updates[i] for i in selected_indices]
        avg = np.mean(selected, axis=0)

        return unflatten_weights(avg, global_weights)


class CoordinateMedian:
    """
    Coordinate-wise median (Yin et al., 2018).
    Tolerates up to f < n/2 Byzantine clients.
    """

    def aggregate(
        self,
        global_weights: Dict[str, np.ndarray],
        client_updates: List[Dict[str, np.ndarray]],
    ) -> Dict[str, np.ndarray]:
        result = {}
        for key in global_weights:
            stacked = np.stack([u[key] for u in client_updates], axis=0)
            result[key] = np.median(stacked, axis=0)
        return result


class TrimmedMean:
    """
    Trimmed mean (Yin et al., 2018).
    Removes top and bottom beta fraction per coordinate.
    """

    def __init__(self, trim_ratio: float = 0.1):
        self.trim_ratio = trim_ratio

    def aggregate(
        self,
        global_weights: Dict[str, np.ndarray],
        client_updates: List[Dict[str, np.ndarray]],
    ) -> Dict[str, np.ndarray]:
        n = len(client_updates)
        k = max(1, int(n * self.trim_ratio))

        result = {}
        for key in global_weights:
            stacked = np.stack([u[key] for u in client_updates], axis=0)
            sorted_vals = np.sort(stacked, axis=0)
            trimmed = sorted_vals[k : n - k]
            if trimmed.shape[0] == 0:
                result[key] = np.mean(stacked, axis=0)
            else:
                result[key] = np.mean(trimmed, axis=0)
        return result


class FLTrust:
    """
    FLTrust (Cao et al., 2021): Server-guided trust scoring.
    Uses cosine similarity to server's reference gradient for weighting.
    """

    def aggregate(
        self,
        global_weights: Dict[str, np.ndarray],
        client_updates: List[Dict[str, np.ndarray]],
        server_update: Dict[str, np.ndarray],
    ) -> Dict[str, np.ndarray]:
        server_flat = flatten_weights(server_update)
        server_norm = np.linalg.norm(server_flat)

        if server_norm < 1e-10:
            result = {}
            for key in global_weights:
                result[key] = np.mean(
                    [u[key] for u in client_updates], axis=0
                )
            return result

        trust_scores = []
        flat_updates = []
        for u in client_updates:
            flat = flatten_weights(u)
            flat_updates.append(flat)
            cosine_sim = np.dot(flat, server_flat) / (
                np.linalg.norm(flat) * server_norm + 1e-10
            )
            trust_scores.append(max(0, cosine_sim))

        total_trust = sum(trust_scores)
        if total_trust < 1e-10:
            trust_scores = [1.0 / len(client_updates)] * len(client_updates)
        else:
            trust_scores = [ts / total_trust for ts in trust_scores]

        result_flat = np.zeros_like(server_flat)
        for flat, ts in zip(flat_updates, trust_scores):
            client_norm = np.linalg.norm(flat)
            if client_norm > 1e-10:
                normalized = flat * (server_norm / client_norm)
            else:
                normalized = flat
            result_flat += ts * normalized

        return unflatten_weights(result_flat, global_weights)


class ByzantineDetector:
    """Norm-based anomaly detection for Byzantine clients."""

    def __init__(self, z_threshold: float = 3.0):
        self.z_threshold = z_threshold
        self.history: List[List[float]] = []

    def detect(
        self, client_updates: List[Dict[str, np.ndarray]]
    ) -> List[int]:
        norms = [
            np.linalg.norm(flatten_weights(u)) for u in client_updates
        ]
        mean_norm = np.mean(norms)
        std_norm = np.std(norms) + 1e-10

        suspicious = []
        for i, norm in enumerate(norms):
            z_score = abs(norm - mean_norm) / std_norm
            if z_score > self.z_threshold:
                suspicious.append(i)

        self.history.append(norms)
        return suspicious
