"""
ablation_prediction.py
======================
Module 6: Atrial fibrillation ablation effect prediction.
Case study: virtual ablation strategies and outcome prediction.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple
from enum import Enum
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class AblationStrategy(Enum):
    PVI = "pulmonary_vein_isolation"
    PVI_ROOF = "pvi_plus_roof_line"
    PVI_MITRAL = "pvi_plus_mitral_isthmus"
    PVI_POSTERIOR = "pvi_plus_posterior_wall"
    CFAE = "complex_fractionated_electrograms"
    ROTOR = "rotor_guided"
    SUBSTRATE = "substrate_guided"
    HYBRID = "hybrid_pvi_substrate"


@dataclass
class AblationLesion:
    """Represents a single ablation lesion."""
    center: np.ndarray        # (3,) position in mm
    radius: float = 3.5       # Lesion radius (mm)
    depth: float = 4.0        # Lesion depth (mm)
    transmurality: float = 1.0  # 0-1
    conductivity_reduction: float = 0.95  # Fraction of conductivity reduction


@dataclass
class AblationPlan:
    """Complete ablation strategy with lesion set."""
    strategy: AblationStrategy
    lesions: List[AblationLesion] = field(default_factory=list)
    total_ablation_time_s: float = 0.0
    rf_power_w: float = 30.0

    @property
    def n_lesions(self) -> int:
        return len(self.lesions)

    @property
    def total_ablated_area_mm2(self) -> float:
        return sum(np.pi * l.radius**2 for l in self.lesions)


@dataclass
class AblationOutcome:
    """Predicted outcome of ablation procedure."""
    strategy: AblationStrategy
    af_terminated: bool
    time_to_termination_ms: float
    post_ablation_af_inducible: bool
    conduction_gaps: int
    pv_reconnection_risk: float   # 0-1
    recurrence_probability_1yr: float  # 0-1
    metrics: Dict = field(default_factory=dict)


class AtrialGeometryModel:
    """
    Left atrial geometry model with anatomical landmarks.

    Key structures:
    - 4 pulmonary veins (LSPV, LIPV, RSPV, RIPV)
    - Mitral valve annulus
    - Left atrial appendage (LAA)
    - Posterior wall
    - Roof
    """

    def __init__(self, n_elements: int = 2000):
        self.n_elements = n_elements
        self._generate_synthetic_atrium()

    def _generate_synthetic_atrium(self):
        """Generate a synthetic left atrial geometry."""
        rng = np.random.default_rng(42)

        # Ellipsoidal atrial body
        theta = rng.uniform(0, 2 * np.pi, self.n_elements)
        phi = rng.uniform(-np.pi / 2, np.pi / 2, self.n_elements)
        r = np.array([25.0, 20.0, 18.0])  # Semi-axes (mm)

        self.coords = np.column_stack([
            r[0] * np.cos(phi) * np.cos(theta),
            r[1] * np.cos(phi) * np.sin(theta),
            r[2] * np.sin(phi),
        ])

        # Assign regions
        self.regions = np.zeros(self.n_elements, dtype=int)
        self._assign_anatomical_regions()

        # Tissue properties
        self.wall_thickness = 1.5 + rng.exponential(0.5, self.n_elements)  # mm
        self.fibrosis_map = (rng.random(self.n_elements) < 0.15).astype(int)

    def _assign_anatomical_regions(self):
        """Assign anatomical region labels based on coordinates."""
        # Simplified region assignment
        # 0: body, 1: LSPV, 2: LIPV, 3: RSPV, 4: RIPV,
        # 5: roof, 6: posterior, 7: mitral, 8: LAA

        pv_centers = {
            1: np.array([-15, 15, 10]),   # LSPV
            2: np.array([-15, 15, -8]),   # LIPV
            3: np.array([15, 15, 10]),    # RSPV
            4: np.array([15, 15, -8]),    # RIPV
        }

        for region_id, center in pv_centers.items():
            dists = np.linalg.norm(self.coords - center, axis=1)
            self.regions[dists < 8.0] = region_id

        # Roof
        self.regions[self.coords[:, 2] > 14] = 5
        # Posterior wall
        self.regions[(self.coords[:, 1] > 15) & (self.regions == 0)] = 6
        # Mitral annulus
        self.regions[self.coords[:, 2] < -14] = 7
        # LAA
        mask = (self.coords[:, 0] < -18) & (self.coords[:, 2] > 0) & (self.regions == 0)
        self.regions[mask] = 8

    def get_pv_ostia_elements(self) -> Dict[str, List[int]]:
        """Get element indices at PV ostia."""
        return {
            "LSPV": np.where(self.regions == 1)[0].tolist(),
            "LIPV": np.where(self.regions == 2)[0].tolist(),
            "RSPV": np.where(self.regions == 3)[0].tolist(),
            "RIPV": np.where(self.regions == 4)[0].tolist(),
        }


class VirtualAblationSimulator:
    """
    Virtual ablation simulation engine.

    Applies ablation lesions to the computational mesh and
    re-runs electrophysiology simulation to assess effect.
    """

    def __init__(self, atrial_model: AtrialGeometryModel):
        self.atrium = atrial_model
        self.conductivity = np.ones(atrial_model.n_elements) * 0.5  # mS/mm

    def create_pvi_lesion_set(self) -> AblationPlan:
        """Create standard PVI (circumferential) lesion set."""
        lesions = []
        pv_ostia = self.atrium.get_pv_ostia_elements()

        for pv_name, elements in pv_ostia.items():
            if not elements:
                continue
            pv_coords = self.atrium.coords[elements]
            center = pv_coords.mean(axis=0)

            # Circumferential lesions around each PV
            n_points = 12
            radius = 10.0  # Distance from PV center
            for i in range(n_points):
                angle = 2 * np.pi * i / n_points
                lesion_pos = center + radius * np.array([
                    np.cos(angle), np.sin(angle), 0
                ])
                lesions.append(AblationLesion(
                    center=lesion_pos,
                    radius=3.5,
                    depth=4.0,
                    transmurality=0.95,
                ))

        return AblationPlan(
            strategy=AblationStrategy.PVI,
            lesions=lesions,
            total_ablation_time_s=len(lesions) * 20,
        )

    def create_strategy_lesion_set(self, strategy: AblationStrategy
                                     ) -> AblationPlan:
        """Create lesion set for a given strategy."""
        if strategy == AblationStrategy.PVI:
            return self.create_pvi_lesion_set()

        plan = self.create_pvi_lesion_set()
        plan.strategy = strategy

        if strategy in (AblationStrategy.PVI_ROOF, AblationStrategy.HYBRID):
            # Add roof line
            n_roof = 8
            for i in range(n_roof):
                x = -15 + 30 * i / (n_roof - 1)
                plan.lesions.append(AblationLesion(
                    center=np.array([x, 15, 16]),
                    radius=3.5, depth=4.0,
                ))

        if strategy in (AblationStrategy.PVI_MITRAL,):
            # Add mitral isthmus line
            n_mitral = 6
            for i in range(n_mitral):
                y = 0 + 15 * i / (n_mitral - 1)
                plan.lesions.append(AblationLesion(
                    center=np.array([-15, y, -12]),
                    radius=3.5, depth=4.0,
                ))

        if strategy == AblationStrategy.PVI_POSTERIOR:
            # Add posterior wall isolation
            n_post = 10
            for i in range(n_post):
                x = -12 + 24 * i / (n_post - 1)
                plan.lesions.append(AblationLesion(
                    center=np.array([x, 18, 0]),
                    radius=4.0, depth=5.0,
                ))

        plan.total_ablation_time_s = len(plan.lesions) * 20
        return plan

    def apply_ablation(self, plan: AblationPlan) -> np.ndarray:
        """Apply ablation lesions to conductivity field."""
        ablated_conductivity = self.conductivity.copy()

        for lesion in plan.lesions:
            dists = np.linalg.norm(
                self.atrium.coords - lesion.center, axis=1
            )
            affected = dists < lesion.radius
            ablated_conductivity[affected] *= (1 - lesion.conductivity_reduction *
                                                lesion.transmurality)

        n_ablated = np.sum(ablated_conductivity < 0.1 * self.conductivity)
        logger.info(f"Applied {plan.n_lesions} lesions, "
                    f"{n_ablated} elements ablated "
                    f"({n_ablated/self.atrium.n_elements*100:.1f}%)")

        return ablated_conductivity

    def simulate_post_ablation(self, plan: AblationPlan,
                                 duration_ms: float = 5000.0
                                 ) -> AblationOutcome:
        """
        Simulate post-ablation AF behavior.

        Tests:
        1. Does AF terminate?
        2. Can AF be re-induced?
        3. Are there conduction gaps?
        """
        ablated_cond = self.apply_ablation(plan)

        # Conduction gap analysis
        gaps = self._detect_conduction_gaps(ablated_cond, plan)

        # AF termination assessment
        n_ablated_frac = np.sum(ablated_cond < 0.05) / self.atrium.n_elements
        fibrosis_burden = np.mean(self.atrium.fibrosis_map)

        # Simplified outcome prediction model
        # Based on published clinical data correlations
        termination_prob = 0.3 + 0.5 * n_ablated_frac - 0.8 * fibrosis_burden
        termination_prob = np.clip(termination_prob, 0, 1)

        af_terminated = termination_prob > 0.5

        # Time to termination (if terminated)
        t_termination = 3000 * (1 - termination_prob) if af_terminated else duration_ms

        # Re-induction test
        reinduction_prob = 0.3 + 0.5 * gaps - 0.3 * n_ablated_frac
        reinduction_prob = np.clip(reinduction_prob, 0, 1)

        # PV reconnection risk (6-month)
        gap_risk = min(1.0, gaps * 0.15)
        transmurality = np.mean([l.transmurality for l in plan.lesions])
        reconnection_risk = gap_risk + 0.3 * (1 - transmurality)

        # 1-year recurrence
        base_recurrence = 0.30  # ~30% for PVI alone
        strategy_modifier = {
            AblationStrategy.PVI: 1.0,
            AblationStrategy.PVI_ROOF: 0.85,
            AblationStrategy.PVI_MITRAL: 0.80,
            AblationStrategy.PVI_POSTERIOR: 0.75,
            AblationStrategy.HYBRID: 0.70,
            AblationStrategy.SUBSTRATE: 0.80,
            AblationStrategy.ROTOR: 0.85,
            AblationStrategy.CFAE: 0.90,
        }
        modifier = strategy_modifier.get(plan.strategy, 1.0)
        recurrence = base_recurrence * modifier + 0.2 * fibrosis_burden

        return AblationOutcome(
            strategy=plan.strategy,
            af_terminated=af_terminated,
            time_to_termination_ms=t_termination,
            post_ablation_af_inducible=reinduction_prob > 0.5,
            conduction_gaps=gaps,
            pv_reconnection_risk=reconnection_risk,
            recurrence_probability_1yr=min(1.0, recurrence),
            metrics={
                "n_lesions": plan.n_lesions,
                "total_ablated_area_mm2": plan.total_ablated_area_mm2,
                "ablation_time_min": plan.total_ablation_time_s / 60,
                "ablated_fraction": n_ablated_frac,
                "fibrosis_burden": fibrosis_burden,
                "termination_probability": termination_prob,
            },
        )

    def _detect_conduction_gaps(self, ablated_cond: np.ndarray,
                                  plan: AblationPlan) -> int:
        """Detect gaps in ablation lines where conduction persists."""
        gaps = 0
        pv_ostia = self.atrium.get_pv_ostia_elements()

        for pv_name, elements in pv_ostia.items():
            if not elements:
                continue
            pv_cond = ablated_cond[elements]
            # Check if any PV ostium element still conducts
            if np.any(pv_cond > 0.1 * np.max(self.conductivity)):
                gaps += 1

        return gaps


class AblationCaseStudy:
    """
    Complete AF ablation case study pipeline.

    1. Patient-specific atrial model construction
    2. Baseline AF simulation
    3. Virtual ablation (multiple strategies)
    4. Comparative outcome prediction
    5. Optimal strategy recommendation
    """

    def __init__(self, patient_id: str = "AF_001"):
        self.patient_id = patient_id
        self.atrium = AtrialGeometryModel(n_elements=2000)
        self.simulator = VirtualAblationSimulator(self.atrium)

    def run_case_study(self) -> Dict[str, any]:
        """Execute complete case study comparing ablation strategies."""
        logger.info(f"Starting AF ablation case study for patient {self.patient_id}")

        strategies = [
            AblationStrategy.PVI,
            AblationStrategy.PVI_ROOF,
            AblationStrategy.PVI_POSTERIOR,
            AblationStrategy.HYBRID,
        ]

        results = {}

        for strategy in strategies:
            logger.info(f"Testing strategy: {strategy.value}")
            plan = self.simulator.create_strategy_lesion_set(strategy)
            outcome = self.simulator.simulate_post_ablation(plan)

            results[strategy.value] = {
                "outcome": outcome,
                "plan": plan,
            }

            logger.info(
                f"  {strategy.value}: "
                f"terminated={outcome.af_terminated}, "
                f"gaps={outcome.conduction_gaps}, "
                f"recurrence={outcome.recurrence_probability_1yr:.2f}"
            )

        # Find optimal strategy
        best_strategy = min(
            results.items(),
            key=lambda x: x[1]["outcome"].recurrence_probability_1yr
        )

        summary = {
            "patient_id": self.patient_id,
            "n_strategies_tested": len(strategies),
            "results": {
                name: {
                    "af_terminated": r["outcome"].af_terminated,
                    "n_lesions": r["outcome"].metrics["n_lesions"],
                    "ablation_time_min": r["outcome"].metrics["ablation_time_min"],
                    "conduction_gaps": r["outcome"].conduction_gaps,
                    "pv_reconnection_risk": round(r["outcome"].pv_reconnection_risk, 3),
                    "recurrence_1yr": round(r["outcome"].recurrence_probability_1yr, 3),
                }
                for name, r in results.items()
            },
            "optimal_strategy": best_strategy[0],
            "optimal_recurrence": round(
                best_strategy[1]["outcome"].recurrence_probability_1yr, 3
            ),
        }

        return summary
