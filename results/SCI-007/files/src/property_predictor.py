"""Property prediction models for synthetic therapeutic antibody design."""
from __future__ import annotations

from typing import Dict, Optional, Sequence

import torch
import torch.nn as nn
import torch.nn.functional as F

from data_utils import (
    AA_TO_IDX,
    AA_VOCAB,
    AA_VOCAB_SIZE,
    CHARGE_TABLE,
    HYDROPHOBICITY_TABLE,
    MOLECULAR_WEIGHT_TABLE,
    PAD_IDX,
    encode_sequence,
)


def _feature_tensor(table: Dict[str, float]) -> torch.Tensor:
    return torch.tensor([table.get(token, 0.0) for token in AA_VOCAB], dtype=torch.float32)


class SequenceEncoder(nn.Module):
    """Differentiable encoder for tokenized or soft sequence inputs."""

    def __init__(self, vocab_size: int = AA_VOCAB_SIZE, hidden_dim: int = 128) -> None:
        super().__init__()
        self.embedding = nn.Embedding(vocab_size, hidden_dim)
        self.conv = nn.Conv1d(hidden_dim, hidden_dim, kernel_size=3, padding=1)
        self.norm = nn.LayerNorm(hidden_dim)
        self.proj = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim),
        )

    def _to_probabilities(self, sequence_input: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        if sequence_input.dtype in (torch.int8, torch.int16, torch.int32, torch.int64, torch.long):
            probs = F.one_hot(sequence_input, num_classes=self.embedding.num_embeddings).float()
            mask = (sequence_input != PAD_IDX).float()
        else:
            probs = sequence_input
            if probs.min().item() < 0.0 or torch.max(torch.abs(probs.sum(dim=-1) - 1.0)).item() > 1e-3:
                probs = F.softmax(probs, dim=-1)
            mask = 1.0 - probs[..., PAD_IDX]
        return probs, mask

    def forward(self, sequence_input: torch.Tensor) -> torch.Tensor:
        probs, mask = self._to_probabilities(sequence_input)
        embeddings = probs @ self.embedding.weight
        hidden = self.conv(embeddings.transpose(1, 2)).transpose(1, 2)
        hidden = self.norm(F.gelu(hidden))
        pooled = (hidden * mask.unsqueeze(-1)).sum(dim=1) / mask.sum(dim=1, keepdim=True).clamp_min(1.0)
        return self.proj(pooled)


class PropertyPredictorBase(nn.Module):
    """Base class with differentiable sequence composition features."""

    def __init__(self, hidden_dim: int = 128) -> None:
        super().__init__()
        self.encoder = SequenceEncoder(hidden_dim=hidden_dim)
        self.register_buffer("charge_vector", _feature_tensor(CHARGE_TABLE))
        self.register_buffer("hydro_vector", _feature_tensor(HYDROPHOBICITY_TABLE))
        self.register_buffer("mass_vector", _feature_tensor(MOLECULAR_WEIGHT_TABLE))
        self.register_buffer("aromatic_vector", torch.tensor([1.0 if aa in "FWY" else 0.0 for aa in AA_VOCAB], dtype=torch.float32))
        self.register_buffer("gly_vector", torch.tensor([1.0 if aa == "G" else 0.0 for aa in AA_VOCAB], dtype=torch.float32))
        self.register_buffer("pro_vector", torch.tensor([1.0 if aa == "P" else 0.0 for aa in AA_VOCAB], dtype=torch.float32))
        self.register_buffer("cys_vector", torch.tensor([1.0 if aa == "C" else 0.0 for aa in AA_VOCAB], dtype=torch.float32))
        self.register_buffer("tyr_vector", torch.tensor([1.0 if aa == "Y" else 0.0 for aa in AA_VOCAB], dtype=torch.float32))
        self.register_buffer("human_like_vector", torch.tensor([1.0 if aa in "AGSDTYV" else 0.0 for aa in AA_VOCAB], dtype=torch.float32))
        self.register_buffer("rare_vector", torch.tensor([1.0 if aa in "CMW" else 0.0 for aa in AA_VOCAB], dtype=torch.float32))
        self.register_buffer("hydrophobe_vector", torch.tensor([1.0 if aa in "AILMFWVY" else 0.0 for aa in AA_VOCAB], dtype=torch.float32))

    def sequence_features(self, sequence_input: torch.Tensor) -> Dict[str, torch.Tensor]:
        probs, mask = self.encoder._to_probabilities(sequence_input)
        denom = mask.sum(dim=1).clamp_min(1.0)

        def averaged(vector: torch.Tensor) -> torch.Tensor:
            return ((probs * vector.view(1, 1, -1)).sum(dim=-1) * mask).sum(dim=1) / denom

        charge = averaged(self.charge_vector)
        hydrophobicity = averaged(self.hydro_vector)
        mass = averaged(self.mass_vector)
        aromaticity = averaged(self.aromatic_vector)
        glycine = averaged(self.gly_vector)
        proline = averaged(self.pro_vector)
        cysteine = averaged(self.cys_vector)
        tyrosine = averaged(self.tyr_vector)
        human_like = averaged(self.human_like_vector)
        rare = averaged(self.rare_vector)
        hydrophobes = averaged(self.hydrophobe_vector)
        return {
            "mask": mask,
            "length": denom,
            "charge": charge,
            "abs_charge": charge.abs(),
            "hydrophobicity": hydrophobicity,
            "mass": mass,
            "aromaticity": aromaticity,
            "glycine": glycine,
            "proline": proline,
            "cysteine": cysteine,
            "tyrosine": tyrosine,
            "human_like": human_like,
            "rare": rare,
            "hydrophobes": hydrophobes,
        }


class BindingAffinityPredictor(PropertyPredictorBase):
    """Predict an approximate binding affinity (pKd)."""

    def __init__(self, hidden_dim: int = 128) -> None:
        super().__init__(hidden_dim=hidden_dim)
        self.residual = nn.Sequential(nn.Linear(hidden_dim, hidden_dim), nn.GELU(), nn.Linear(hidden_dim, 1))
        nn.init.zeros_(self.residual[-1].weight)
        nn.init.zeros_(self.residual[-1].bias)

    def forward(self, sequence_input: torch.Tensor) -> torch.Tensor:
        features = self.sequence_features(sequence_input)
        embedding = self.encoder(sequence_input)
        heuristic = (
            6.2
            + 2.6 * features["aromaticity"]
            + 0.7 * features["tyrosine"]
            + 0.22 * features["hydrophobicity"]
            - 0.35 * features["abs_charge"]
        )
        return torch.clamp(heuristic + 0.1 * self.residual(embedding).squeeze(-1), 4.0, 13.0)


class StabilityPredictor(PropertyPredictorBase):
    """Predict an approximate thermal stability (Tm in Celsius)."""

    def __init__(self, hidden_dim: int = 128) -> None:
        super().__init__(hidden_dim=hidden_dim)
        self.residual = nn.Sequential(nn.Linear(hidden_dim, hidden_dim), nn.GELU(), nn.Linear(hidden_dim, 1))
        nn.init.zeros_(self.residual[-1].weight)
        nn.init.zeros_(self.residual[-1].bias)

    def forward(self, sequence_input: torch.Tensor) -> torch.Tensor:
        features = self.sequence_features(sequence_input)
        embedding = self.encoder(sequence_input)
        heuristic = 58.0 + 7.5 * features["hydrophobes"] - 8.0 * features["glycine"] - 5.0 * features["proline"] - 6.0 * features["cysteine"]
        return torch.clamp(heuristic + 0.2 * self.residual(embedding).squeeze(-1), 35.0, 95.0)


class HumanizationScorer(PropertyPredictorBase):
    """Predict a humanness score based on germline-like sequence patterns."""

    def __init__(self, hidden_dim: int = 128) -> None:
        super().__init__(hidden_dim=hidden_dim)
        self.residual = nn.Sequential(nn.Linear(hidden_dim, hidden_dim), nn.GELU(), nn.Linear(hidden_dim, 1))
        nn.init.zeros_(self.residual[-1].weight)
        nn.init.zeros_(self.residual[-1].bias)

    def forward(self, sequence_input: torch.Tensor) -> torch.Tensor:
        features = self.sequence_features(sequence_input)
        embedding = self.encoder(sequence_input)
        logits = 2.8 * features["human_like"] - 1.8 * features["rare"] - 0.8 * features["cysteine"] + 0.1 * self.residual(embedding).squeeze(-1)
        return torch.sigmoid(logits)


class ImmunogenicityPredictor(PropertyPredictorBase):
    """Predict an immunogenicity risk score."""

    def __init__(self, hidden_dim: int = 128) -> None:
        super().__init__(hidden_dim=hidden_dim)
        self.residual = nn.Sequential(nn.Linear(hidden_dim, hidden_dim), nn.GELU(), nn.Linear(hidden_dim, 1))
        nn.init.zeros_(self.residual[-1].weight)
        nn.init.zeros_(self.residual[-1].bias)

    def forward(self, sequence_input: torch.Tensor) -> torch.Tensor:
        features = self.sequence_features(sequence_input)
        humanization = HumanizationScorer.forward(self, sequence_input)
        embedding = self.encoder(sequence_input)
        logits = 1.4 * features["rare"] + 0.7 * features["abs_charge"] + 0.25 * features["hydrophobes"] - 1.6 * humanization + 0.1 * self.residual(embedding).squeeze(-1)
        return torch.sigmoid(logits)


class DevelopabilityPredictor(PropertyPredictorBase):
    """Predict expression level and aggregation propensity."""

    def __init__(self, hidden_dim: int = 128) -> None:
        super().__init__(hidden_dim=hidden_dim)
        self.residual = nn.Sequential(nn.Linear(hidden_dim, hidden_dim), nn.GELU(), nn.Linear(hidden_dim, 2))
        nn.init.zeros_(self.residual[-1].weight)
        nn.init.zeros_(self.residual[-1].bias)

    def forward(self, sequence_input: torch.Tensor) -> Dict[str, torch.Tensor]:
        features = self.sequence_features(sequence_input)
        embedding = self.encoder(sequence_input)
        residual = 0.1 * self.residual(embedding)
        expression = torch.sigmoid(2.1 - 1.2 * features["abs_charge"] - 0.8 * features["rare"] + 0.4 * features["human_like"] + residual[:, 0])
        aggregation = torch.sigmoid(-1.4 + 2.2 * features["hydrophobes"] + 1.1 * features["aromaticity"] - 0.5 * features["glycine"] + residual[:, 1])
        return {"expression_level": expression, "aggregation_propensity": aggregation}


class MultiPropertyOptimizer(nn.Module):
    """Combine individual predictors into a weighted multi-objective score."""

    def __init__(
        self,
        binding_predictor: Optional[BindingAffinityPredictor] = None,
        stability_predictor: Optional[StabilityPredictor] = None,
        humanization_scorer: Optional[HumanizationScorer] = None,
        immunogenicity_predictor: Optional[ImmunogenicityPredictor] = None,
        developability_predictor: Optional[DevelopabilityPredictor] = None,
        weights: Optional[Dict[str, float]] = None,
    ) -> None:
        super().__init__()
        self.binding_predictor = binding_predictor or BindingAffinityPredictor()
        self.stability_predictor = stability_predictor or StabilityPredictor()
        self.humanization_scorer = humanization_scorer or HumanizationScorer()
        self.immunogenicity_predictor = immunogenicity_predictor or ImmunogenicityPredictor()
        self.developability_predictor = developability_predictor or DevelopabilityPredictor()
        self.weights = weights or {
            "binding_affinity": 0.28,
            "stability_tm": 0.18,
            "humanness": 0.17,
            "non_immunogenic": 0.15,
            "expression_level": 0.12,
            "low_aggregation": 0.10,
        }

    def evaluate(self, sequence_input: torch.Tensor | Sequence[str]) -> Dict[str, torch.Tensor]:
        if isinstance(sequence_input, (list, tuple)) and sequence_input and isinstance(sequence_input[0], str):
            sequence_input = torch.stack([encode_sequence(sequence) for sequence in sequence_input])
        binding = self.binding_predictor(sequence_input)
        stability = self.stability_predictor(sequence_input)
        humanness = self.humanization_scorer(sequence_input)
        immunogenicity = self.immunogenicity_predictor(sequence_input)
        developability = self.developability_predictor(sequence_input)
        return {
            "binding_affinity": binding,
            "stability_tm": stability,
            "humanness": humanness,
            "immunogenicity": immunogenicity,
            "expression_level": developability["expression_level"],
            "aggregation_propensity": developability["aggregation_propensity"],
        }

    def objective_vector(self, sequence_input: torch.Tensor | Sequence[str]) -> torch.Tensor:
        metrics = self.evaluate(sequence_input)
        return torch.stack(
            [
                metrics["binding_affinity"] / 12.0,
                metrics["stability_tm"] / 100.0,
                metrics["humanness"],
                1.0 - metrics["immunogenicity"],
                metrics["expression_level"],
                1.0 - metrics["aggregation_propensity"],
            ],
            dim=-1,
        )

    def weighted_score(self, sequence_input: torch.Tensor | Sequence[str]) -> torch.Tensor:
        objectives = self.objective_vector(sequence_input)
        weight_vector = torch.tensor(
            [
                self.weights["binding_affinity"],
                self.weights["stability_tm"],
                self.weights["humanness"],
                self.weights["non_immunogenic"],
                self.weights["expression_level"],
                self.weights["low_aggregation"],
            ],
            dtype=objectives.dtype,
            device=objectives.device,
        )
        return (objectives * weight_vector).sum(dim=-1)


__all__ = [
    "BindingAffinityPredictor",
    "DevelopabilityPredictor",
    "HumanizationScorer",
    "ImmunogenicityPredictor",
    "MultiPropertyOptimizer",
    "SequenceEncoder",
    "StabilityPredictor",
]
