"""
Communication efficiency strategies for federated learning.

Implements:
  - Top-K gradient sparsification
  - Stochastic gradient quantization (QSGD)
  - Error feedback (memory) for unbiased compression
  - Knowledge distillation for model-agnostic FL
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass


@dataclass
class CompressionStats:
    """Track compression statistics."""
    original_size: int
    compressed_size: int
    compression_ratio: float
    reconstruction_error: float


class TopKSparsifier:
    """
    Top-K gradient sparsification with error feedback.
    Communication cost: O(K) instead of O(d).
    """

    def __init__(self, compression_ratio: float = 0.01):
        self.compression_ratio = compression_ratio
        self.error_feedback: Dict[str, np.ndarray] = {}

    def compress(
        self, gradients: Dict[str, np.ndarray], use_error_feedback: bool = True
    ) -> Tuple[Dict[str, np.ndarray], CompressionStats]:
        original_size = 0
        compressed_size = 0
        total_error = 0.0

        result = {}
        for key, grad in gradients.items():
            if use_error_feedback and key in self.error_feedback:
                grad = grad + self.error_feedback[key]

            flat = grad.flatten()
            original_size += flat.size
            k = max(1, int(flat.size * self.compression_ratio))

            indices = np.argpartition(np.abs(flat), -k)[-k:]
            mask = np.zeros_like(flat)
            mask[indices] = flat[indices]

            compressed = mask.reshape(grad.shape)
            result[key] = compressed
            compressed_size += k

            if use_error_feedback:
                self.error_feedback[key] = grad - compressed
            total_error += np.sum((grad - compressed) ** 2)

        stats = CompressionStats(
            original_size=original_size,
            compressed_size=compressed_size,
            compression_ratio=compressed_size / max(original_size, 1),
            reconstruction_error=float(np.sqrt(total_error)),
        )
        return result, stats


class StochasticQuantizer:
    """Stochastic gradient quantization (QSGD, Alistarh et al., 2017)."""

    def __init__(self, num_bits: int = 8):
        self.num_bits = num_bits
        self.num_levels = 2 ** num_bits - 1

    def quantize(
        self, gradients: Dict[str, np.ndarray]
    ) -> Tuple[Dict[str, np.ndarray], CompressionStats]:
        original_size = 0
        compressed_size = 0
        total_error = 0.0
        result = {}

        for key, grad in gradients.items():
            original_size += grad.size * 32

            norm = np.linalg.norm(grad)
            if norm < 1e-10:
                result[key] = np.zeros_like(grad)
                compressed_size += grad.size * self.num_bits
                continue

            normalized = np.abs(grad) / norm
            scaled = normalized * self.num_levels
            lower = np.floor(scaled).astype(int)
            prob = scaled - lower

            rng = np.random.random(grad.shape)
            quantized_abs = np.where(rng < prob, lower + 1, lower)

            signs = np.sign(grad)
            reconstructed = signs * (quantized_abs / self.num_levels) * norm
            result[key] = reconstructed

            compressed_size += grad.size * self.num_bits
            total_error += np.sum((grad - reconstructed) ** 2)

        stats = CompressionStats(
            original_size=original_size,
            compressed_size=compressed_size,
            compression_ratio=compressed_size / max(original_size, 1),
            reconstruction_error=float(np.sqrt(total_error)),
        )
        return result, stats


class KnowledgeDistillationFL:
    """
    Knowledge distillation for model-agnostic FL (FedMD/FedDF).
    Shares predictions on proxy dataset instead of model weights.
    """

    def __init__(
        self,
        temperature: float = 3.0,
        alpha: float = 0.5,
        proxy_dataset_size: int = 1000,
    ):
        self.temperature = temperature
        self.alpha = alpha
        self.proxy_dataset_size = proxy_dataset_size

    def softmax_with_temperature(self, logits: np.ndarray) -> np.ndarray:
        scaled = logits / self.temperature
        exp_vals = np.exp(scaled - np.max(scaled, axis=-1, keepdims=True))
        return exp_vals / np.sum(exp_vals, axis=-1, keepdims=True)

    def aggregate_logits(
        self,
        client_logits: List[np.ndarray],
        weights: Optional[List[float]] = None,
    ) -> np.ndarray:
        if weights is None:
            weights = [1.0 / len(client_logits)] * len(client_logits)

        consensus = np.zeros_like(client_logits[0])
        for logit, w in zip(client_logits, weights):
            consensus += w * self.softmax_with_temperature(logit)
        return consensus

    def distillation_loss(
        self,
        student_logits: np.ndarray,
        teacher_soft_labels: np.ndarray,
    ) -> float:
        student_probs = self.softmax_with_temperature(student_logits)
        eps = 1e-10
        kl = np.sum(
            teacher_soft_labels
            * np.log((teacher_soft_labels + eps) / (student_probs + eps))
        )
        return float(kl) * (self.temperature ** 2)

    def compute_communication_savings(
        self, model_params: int, num_classes: int
    ) -> Dict[str, float]:
        weight_comm = model_params * 4
        logit_comm = self.proxy_dataset_size * num_classes * 4
        return {
            "weight_sharing_bytes": weight_comm,
            "logit_sharing_bytes": logit_comm,
            "savings_ratio": 1 - logit_comm / weight_comm,
            "speedup": weight_comm / max(logit_comm, 1),
        }
