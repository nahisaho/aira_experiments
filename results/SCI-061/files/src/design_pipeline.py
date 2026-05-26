"""
Cello/SBOL-based Automated Design Pipeline
Integrates all modules into a unified workflow.
"""

import json
from typing import Dict, List, Optional, Tuple
from src.parts_catalog import PartsCatalog, PartType, GeneticPart
from src.circuit_spec import (
    CircuitSpec, LogicGate, GateType, Signal,
    make_toggle_switch, make_repressilator
)
from src.stochastic_sim import (
    StochasticModel, gillespie_ssa, tau_leaping,
    build_toggle_switch_model, build_repressilator_model
)
from src.robust_design import (
    robustness_score, optimize_circuit_params,
    toggle_switch_bistability_score, repressilator_oscillation_score,
    latin_hypercube_sample
)
from src.context_effects import ContextPredictor


class DesignPipeline:
    """
    Automated gene circuit design pipeline (Cello/SBOL-inspired).

    Workflow:
    1. Parse circuit specification
    2. Map to biological parts from catalog
    3. Predict & correct context effects
    4. Stochastic simulation
    5. Robustness analysis under parameter uncertainty
    6. Optimization
    7. Export in SBOL-like format
    """

    def __init__(self):
        self.catalog = PartsCatalog()
        self.context_predictor = ContextPredictor()
        self.results = {}

    def design_circuit(
        self,
        spec: CircuitSpec,
        sim_method: str = "tau_leaping",
        t_end: float = 500.0,
        optimize: bool = True,
        n_opt_iterations: int = 30,
    ) -> Dict:
        """Run complete design pipeline for a circuit specification."""
        result = {
            "circuit_name": spec.name,
            "n_gates": len(spec.gates),
            "n_feedbacks": len(spec.feedbacks),
        }

        # 1. Part mapping
        part_assignments = self._map_parts(spec)
        result["part_assignments"] = part_assignments

        # 2. Assembly order and context effects
        assembly_order = self._get_assembly_order(spec)
        context_effects = self.context_predictor.compute_circuit_context_effects(
            assembly_order
        )
        result["context_effects"] = [
            {"upstream": u, "downstream": d, "fold_change": round(fc, 3)}
            for u, d, fc in context_effects
        ]

        # 3. Insulation recommendations
        insulation = self.context_predictor.insulation_recommendation(
            assembly_order
        )
        result["insulation_recommendations"] = insulation

        # 4. Verilog export
        if not spec.feedbacks:  # Only for combinational
            result["verilog"] = spec.to_verilog()

        # 5. SBOL export
        result["sbol_parts"] = json.loads(self.catalog.to_sbol_json())

        return result

    def _map_parts(self, spec: CircuitSpec) -> List[Dict]:
        """Map gate specifications to catalog parts."""
        assignments = []
        for gate in spec.gates:
            assignment = {"gate_id": gate.gate_id, "gate_type": gate.gate_type.value}
            for attr, ptype in [
                ("promoter", PartType.PROMOTER),
                ("rbs", PartType.RBS),
                ("cds", PartType.CDS),
                ("terminator", PartType.TERMINATOR),
            ]:
                name = getattr(gate, attr, "")
                if name and name in self.catalog.parts:
                    part = self.catalog.get(name)
                    assignment[attr] = {
                        "name": part.name,
                        "sbol_uri": part.sbol_uri,
                    }
                elif name:
                    assignment[attr] = {"name": name, "sbol_uri": ""}
            assignments.append(assignment)
        return assignments

    def _get_assembly_order(self, spec: CircuitSpec) -> List[str]:
        """Determine linear assembly order of parts."""
        order = []
        for gate in spec.gates:
            if gate.promoter:
                order.append(gate.promoter)
            if gate.rbs:
                order.append(gate.rbs)
            if gate.cds:
                order.append(gate.cds)
            if gate.terminator:
                order.append(gate.terminator)
        return order

    def run_simulation(
        self,
        circuit_type: str,
        method: str = "tau_leaping",
        params: Optional[Dict] = None,
        t_end: float = 500.0,
        tau: float = 0.5,
        seed: int = 42,
    ):
        """Run stochastic simulation for a circuit type."""
        if circuit_type == "toggle_switch":
            model = build_toggle_switch_model(params)
        elif circuit_type == "repressilator":
            model = build_repressilator_model(params)
        else:
            raise ValueError(f"Unknown circuit type: {circuit_type}")

        if method == "gillespie":
            return gillespie_ssa(model, t_end, seed)
        else:
            return tau_leaping(model, t_end, tau, seed)

    def run_robustness_analysis(
        self,
        circuit_type: str,
        param_ranges: Dict[str, Tuple[float, float]],
        n_samples: int = 50,
        t_end: float = 500.0,
    ):
        """Run robustness analysis under parameter uncertainty."""
        if circuit_type == "toggle_switch":
            builder = build_toggle_switch_model
            obj_func = toggle_switch_bistability_score
        else:
            builder = build_repressilator_model
            obj_func = repressilator_oscillation_score

        return robustness_score(
            builder, param_ranges, obj_func,
            n_samples=n_samples, t_end=t_end
        )

    def run_optimization(
        self,
        circuit_type: str,
        param_ranges: Dict[str, Tuple[float, float]],
        n_iterations: int = 30,
        t_end: float = 500.0,
    ):
        """Optimize circuit parameters for robust performance."""
        if circuit_type == "toggle_switch":
            builder = build_toggle_switch_model
            obj_func = toggle_switch_bistability_score
        else:
            builder = build_repressilator_model
            obj_func = repressilator_oscillation_score

        return optimize_circuit_params(
            builder, param_ranges, obj_func,
            n_iterations=n_iterations, t_end=t_end
        )
