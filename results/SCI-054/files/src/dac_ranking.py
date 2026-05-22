"""
Module 6: DAC (Direct Air Capture) MOF ranking system.

Multi-criteria ranking and optimization for identifying top MOF candidates
for atmospheric CO2 capture applications.
"""
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class DACCandidate:
    """Ranked MOF candidate for DAC application."""
    mof_id: str = ""
    source_db: str = ""
    rank: int = 0

    # Performance metrics
    co2_uptake_dac: float = 0.0        # mmol/g at 420 ppm CO2
    co2_uptake_1bar: float = 0.0       # mmol/g at 1 bar
    working_capacity: float = 0.0       # mmol/g (ads - des)
    co2_n2_selectivity: float = 0.0
    heat_of_adsorption: float = 0.0     # kJ/mol
    regeneration_energy: float = 0.0    # kJ/mol estimated

    # Stability
    water_stability: float = 0.0
    thermal_stability_K: float = 0.0
    synthesizability: float = 0.0

    # Geometric properties
    lcd: float = 0.0
    pld: float = 0.0
    asa: float = 0.0
    porosity: float = 0.0
    metal_type: str = ""
    linker_type: str = ""

    # Composite scores
    performance_score: float = 0.0
    practicality_score: float = 0.0
    overall_score: float = 0.0

    def to_dict(self) -> Dict:
        return {
            "rank": self.rank,
            "mof_id": self.mof_id,
            "source_db": self.source_db,
            "co2_uptake_dac_mmol_g": round(self.co2_uptake_dac, 3),
            "working_capacity_mmol_g": round(self.working_capacity, 3),
            "co2_n2_selectivity": round(self.co2_n2_selectivity, 1),
            "Qst_kJ_mol": round(self.heat_of_adsorption, 1),
            "regeneration_energy_kJ_mol": round(self.regeneration_energy, 1),
            "water_stability": round(self.water_stability, 2),
            "synthesizability": round(self.synthesizability, 2),
            "metal_type": self.metal_type,
            "LCD_A": round(self.lcd, 2),
            "ASA_m2_g": round(self.asa, 1),
            "performance_score": round(self.performance_score, 4),
            "practicality_score": round(self.practicality_score, 4),
            "overall_score": round(self.overall_score, 4),
        }


class DACRanker:
    """
    Multi-criteria ranking system for DAC MOF candidates.

    Ranking criteria (weighted):
    1. CO2 working capacity under DAC conditions (30%)
    2. CO2/N2 selectivity (20%)
    3. Regeneration energy efficiency (15%)
    4. Water stability (15%)
    5. Synthesizability (10%)
    6. Optimal heat of adsorption window (10%)
    """

    DEFAULT_WEIGHTS = {
        "working_capacity": 0.30,
        "selectivity": 0.20,
        "regeneration": 0.15,
        "water_stability": 0.15,
        "synthesizability": 0.10,
        "qst_optimality": 0.10,
    }

    # Optimal Qst range for DAC (kJ/mol)
    # Too low: weak binding, poor selectivity
    # Too high: difficult regeneration
    QST_OPTIMAL_RANGE = (30.0, 50.0)

    # DAC operating conditions
    DAC_CONDITIONS = {
        "T_ads": 298.0,           # K
        "T_des": 373.0,           # K (TSA desorption)
        "P_co2_ads": 0.000420,    # bar (420 ppm)
        "P_co2_des": 1.0,         # bar (pure CO2 product)
        "humidity_rh": 0.50,      # 50% RH typical
    }

    def __init__(self, weights: Optional[Dict[str, float]] = None,
                 min_working_capacity: float = 1.0,
                 min_selectivity: float = 50.0):
        self.weights = weights or self.DEFAULT_WEIGHTS
        self.min_wc = min_working_capacity
        self.min_sel = min_selectivity

    def rank_candidates(self, candidates: List[Dict],
                         top_n: int = 50) -> List[DACCandidate]:
        """Rank MOF candidates for DAC suitability."""
        scored = []
        for cand in candidates:
            dac = DACCandidate(
                mof_id=cand.get("mof_id", ""),
                source_db=cand.get("source_db", ""),
                co2_uptake_dac=cand.get("co2_uptake_dac", 0),
                co2_uptake_1bar=cand.get("co2_uptake_1bar", 0),
                working_capacity=cand.get("working_capacity", 0),
                co2_n2_selectivity=cand.get("selectivity", 0),
                heat_of_adsorption=cand.get("Qst", 0),
                water_stability=cand.get("water_stability", 0),
                thermal_stability_K=cand.get("thermal_stability_K", 0),
                synthesizability=cand.get("synthesizability", 0),
                lcd=cand.get("LCD", 0),
                pld=cand.get("PLD", 0),
                asa=cand.get("ASA", 0),
                porosity=cand.get("porosity", 0),
                metal_type=cand.get("metal_type", ""),
                linker_type=cand.get("linker_type", ""),
            )

            # Estimate regeneration energy
            dac.regeneration_energy = self._estimate_regeneration_energy(
                dac.heat_of_adsorption, dac.working_capacity
            )

            # Compute sub-scores
            dac.performance_score = self._performance_score(dac)
            dac.practicality_score = self._practicality_score(dac)
            dac.overall_score = (
                0.6 * dac.performance_score +
                0.4 * dac.practicality_score
            )
            scored.append(dac)

        # Sort by overall score
        scored.sort(key=lambda x: x.overall_score, reverse=True)

        # Assign ranks
        for i, dac in enumerate(scored):
            dac.rank = i + 1

        return scored[:top_n]

    def _performance_score(self, dac: DACCandidate) -> float:
        """Compute performance sub-score (adsorption metrics)."""
        # Working capacity: normalize to 0-1 (target: 1-5 mmol/g)
        wc_score = min(dac.working_capacity / 5.0, 1.0)

        # Selectivity: normalize (target: 50-1000)
        sel_score = min(np.log10(max(dac.co2_n2_selectivity, 1)) / 3.0, 1.0)

        # Qst optimality: bell curve around optimal range
        qst = dac.heat_of_adsorption
        qst_low, qst_high = self.QST_OPTIMAL_RANGE
        if qst_low <= qst <= qst_high:
            qst_score = 1.0
        elif qst < qst_low:
            qst_score = max(0, 1.0 - (qst_low - qst) / 20.0)
        else:
            qst_score = max(0, 1.0 - (qst - qst_high) / 30.0)

        # Regeneration efficiency
        regen_score = max(0, 1.0 - dac.regeneration_energy / 100.0)

        return (
            self.weights["working_capacity"] * wc_score +
            self.weights["selectivity"] * sel_score +
            self.weights["qst_optimality"] * qst_score +
            self.weights["regeneration"] * regen_score
        ) / (self.weights["working_capacity"] + self.weights["selectivity"] +
             self.weights["qst_optimality"] + self.weights["regeneration"])

    def _practicality_score(self, dac: DACCandidate) -> float:
        """Compute practicality sub-score (stability + synthesizability)."""
        ws = dac.water_stability
        ss = dac.synthesizability
        ts = min(dac.thermal_stability_K / 773.0, 1.0)  # Normalize to UiO-66

        return (
            self.weights["water_stability"] * ws +
            self.weights["synthesizability"] * ss
        ) / (self.weights["water_stability"] + self.weights["synthesizability"])

    def _estimate_regeneration_energy(self, qst: float,
                                       working_capacity: float) -> float:
        """
        Estimate total regeneration energy (kJ/mol CO2).

        Components:
        1. Heat of desorption ≈ Qst
        2. Sensible heat for framework heating
        3. Heat loss (estimated 20%)
        """
        if working_capacity <= 0 or qst <= 0:
            return float("inf")

        desorption_energy = qst
        sensible_heat = 10.0  # kJ/mol estimated for ΔT = 75K
        heat_loss_factor = 1.20

        return (desorption_energy + sensible_heat) * heat_loss_factor

    def apply_hard_filters(self, candidates: List[DACCandidate]) -> List[DACCandidate]:
        """Apply hard DAC constraints."""
        filtered = []
        for dac in candidates:
            reasons = []
            if dac.working_capacity < self.min_wc:
                reasons.append(f"WC={dac.working_capacity:.2f} < {self.min_wc}")
            if dac.co2_n2_selectivity < self.min_sel:
                reasons.append(f"Sel={dac.co2_n2_selectivity:.0f} < {self.min_sel}")
            if dac.water_stability < 0.5:
                reasons.append(f"WS={dac.water_stability:.2f} < 0.5")
            if dac.thermal_stability_K < 373.0:
                reasons.append(f"Td={dac.thermal_stability_K:.0f} K < 373 K")

            if not reasons:
                filtered.append(dac)
            else:
                logger.debug(f"Filtered {dac.mof_id}: {'; '.join(reasons)}")

        logger.info(f"Hard filter: {len(filtered)}/{len(candidates)} passed")
        return filtered

    def generate_report(self, ranked: List[DACCandidate],
                         output_path: Path) -> None:
        """Generate JSON ranking report."""
        report = {
            "ranking_criteria": {
                "weights": self.weights,
                "Qst_optimal_kJ_mol": list(self.QST_OPTIMAL_RANGE),
                "DAC_conditions": self.DAC_CONDITIONS,
            },
            "n_candidates_total": len(ranked),
            "top_candidates": [c.to_dict() for c in ranked],
            "statistics": {
                "mean_overall_score": round(np.mean([c.overall_score for c in ranked]), 4),
                "mean_working_capacity": round(np.mean([c.working_capacity for c in ranked]), 3),
                "mean_selectivity": round(np.mean([c.co2_n2_selectivity for c in ranked]), 1),
                "metal_distribution": self._metal_distribution(ranked),
            },
        }

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(report, f, indent=2)
        logger.info(f"DAC ranking report saved to {output_path}")

    def _metal_distribution(self, candidates: List[DACCandidate]) -> Dict[str, int]:
        dist: Dict[str, int] = {}
        for c in candidates:
            dist[c.metal_type] = dist.get(c.metal_type, 0) + 1
        return dict(sorted(dist.items(), key=lambda x: x[1], reverse=True))


class ParetoFrontAnalysis:
    """Multi-objective Pareto front analysis for MOF selection."""

    @staticmethod
    def compute_pareto_front(objectives: np.ndarray,
                              maximize: Optional[List[bool]] = None) -> np.ndarray:
        """
        Identify Pareto-optimal points.

        Args:
            objectives: (n_points, n_objectives) array
            maximize: list of bool for each objective (True=maximize)

        Returns:
            Boolean mask of Pareto-optimal points
        """
        n = len(objectives)
        if maximize is None:
            maximize = [True] * objectives.shape[1]

        # Convert to maximization
        obj = objectives.copy()
        for j, m in enumerate(maximize):
            if not m:
                obj[:, j] = -obj[:, j]

        is_pareto = np.ones(n, dtype=bool)
        for i in range(n):
            if not is_pareto[i]:
                continue
            for j in range(n):
                if i == j or not is_pareto[j]:
                    continue
                if np.all(obj[j] >= obj[i]) and np.any(obj[j] > obj[i]):
                    is_pareto[i] = False
                    break

        return is_pareto

    @staticmethod
    def pareto_rank(objectives: np.ndarray,
                     maximize: Optional[List[bool]] = None) -> np.ndarray:
        """Assign Pareto ranks (1 = front, 2 = second layer, etc.)."""
        n = len(objectives)
        ranks = np.zeros(n, dtype=int)
        remaining = np.ones(n, dtype=bool)
        rank = 1

        while remaining.any():
            subset = objectives[remaining]
            indices = np.where(remaining)[0]
            pareto = ParetoFrontAnalysis.compute_pareto_front(subset, maximize)
            for k, idx in enumerate(indices):
                if pareto[k]:
                    ranks[idx] = rank
                    remaining[idx] = False
            rank += 1

        return ranks
