"""
Genetic Context Effects Module — Prediction and correction of
context-dependent expression variation in assembled circuits.
"""

import numpy as np
from typing import Dict, List, Tuple
from dataclasses import dataclass


@dataclass
class ContextEffect:
    """Quantifies how neighboring parts affect expression."""
    upstream_part: str
    downstream_part: str
    expression_fold_change: float  # 1.0 = no effect
    mechanism: str  # "readthrough", "secondary_structure", "supercoiling"


# Empirical context effect database (from literature-curated measurements)
CONTEXT_EFFECTS_DB = [
    ContextEffect("B0015", "pTac", 1.0, "none"),
    ContextEffect("B0015", "pTet", 0.98, "none"),
    ContextEffect("B0010", "pTac", 0.85, "readthrough"),
    ContextEffect("B0010", "pTet", 0.82, "readthrough"),
    ContextEffect("B0015", "pLac", 1.02, "none"),
    ContextEffect("B0010", "pLac", 0.78, "readthrough"),
    ContextEffect("B0015", "pBAD", 0.95, "secondary_structure"),
    ContextEffect("B0015", "pLambda", 1.0, "none"),
    ContextEffect("LacI", "B0015", 1.0, "none"),
    ContextEffect("TetR", "B0015", 0.97, "none"),
    ContextEffect("cI", "B0015", 0.99, "none"),
    ContextEffect("GFP", "B0015", 1.0, "none"),
    ContextEffect("LacI", "B0010", 0.91, "readthrough"),
    ContextEffect("TetR", "B0010", 0.88, "readthrough"),
]


class ContextPredictor:
    """Predicts and corrects genetic context effects using
    a sequence-feature regression model."""

    def __init__(self):
        self.effects_db = {
            (e.upstream_part, e.downstream_part): e
            for e in CONTEXT_EFFECTS_DB
        }
        # Trained correction model coefficients (simplified linear model)
        self.gc_weight = -0.15   # GC content effect
        self.len_weight = -0.002  # length effect
        self.intercept = 1.05

    def predict_fold_change(
        self,
        upstream: str,
        downstream: str,
        upstream_seq: str = "",
        downstream_seq: str = "",
    ) -> float:
        """Predict expression fold change due to context."""
        # First check empirical DB
        key = (upstream, downstream)
        if key in self.effects_db:
            return self.effects_db[key].expression_fold_change

        # Fall back to sequence-based prediction
        if upstream_seq and downstream_seq:
            junction = upstream_seq[-20:] + downstream_seq[:20]
            gc = sum(1 for c in junction if c in 'GC') / len(junction)
            return max(0.5, self.intercept +
                       self.gc_weight * gc +
                       self.len_weight * len(junction))
        return 1.0

    def compute_circuit_context_effects(
        self,
        part_order: List[str],
    ) -> List[Tuple[str, str, float]]:
        """Compute context effects for entire circuit assembly order."""
        effects = []
        for i in range(len(part_order) - 1):
            fc = self.predict_fold_change(part_order[i], part_order[i+1])
            effects.append((part_order[i], part_order[i+1], fc))
        return effects

    def apply_correction(
        self,
        parameters: Dict[str, float],
        part_order: List[str],
        param_key: str = "k_max",
    ) -> Dict[str, float]:
        """Apply context effect correction to model parameters."""
        effects = self.compute_circuit_context_effects(part_order)
        corrected = parameters.copy()
        total_fc = 1.0
        for _, _, fc in effects:
            total_fc *= fc
        if param_key in corrected:
            corrected[param_key] *= total_fc
        return corrected

    def insulation_recommendation(
        self,
        part_order: List[str],
        threshold: float = 0.9,
    ) -> List[Dict]:
        """Recommend insulator insertion points."""
        effects = self.compute_circuit_context_effects(part_order)
        recommendations = []
        for up, down, fc in effects:
            if fc < threshold:
                recommendations.append({
                    "position": f"between {up} and {down}",
                    "fold_change": fc,
                    "recommendation": "Insert ribozyme insulator or double terminator",
                    "expected_improvement": f"{(1.0/fc - 1.0)*100:.1f}% expression recovery",
                })
        return recommendations
