"""
Module 6: Multi-Objective Optimization for Lead Optimization (Pareto Front)

Implements multi-objective optimization strategies for lead optimization,
balancing potency, selectivity, ADMET properties, and synthetic accessibility.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Callable
from enum import Enum


class ObjectiveType(Enum):
    MINIMIZE = "minimize"
    MAXIMIZE = "maximize"


@dataclass
class Objective:
    """Single optimization objective."""
    name: str
    objective_type: ObjectiveType
    weight: float = 1.0
    target_value: Optional[float] = None
    constraint_min: Optional[float] = None
    constraint_max: Optional[float] = None
    unit: str = ""
    
    def is_satisfied(self, value: float) -> bool:
        if self.constraint_min is not None and value < self.constraint_min:
            return False
        if self.constraint_max is not None and value > self.constraint_max:
            return False
        return True


@dataclass
class CandidateSolution:
    """A candidate in the multi-objective optimization."""
    compound_id: str
    smiles: str
    objective_values: Dict[str, float]
    is_pareto_optimal: bool = False
    crowding_distance: float = 0.0
    domination_count: int = 0
    dominated_solutions: List[str] = field(default_factory=list)
    pareto_rank: int = 0
    
    # Properties
    properties: Dict[str, float] = field(default_factory=dict)


@dataclass
class ParetoFront:
    """Pareto front from multi-objective optimization."""
    solutions: List[CandidateSolution]
    n_objectives: int
    objective_names: List[str]
    
    # Metrics
    hypervolume: float = 0.0
    spacing: float = 0.0
    spread: float = 0.0
    n_pareto_optimal: int = 0
    
    # Optimization history
    generation_history: List[Dict] = field(default_factory=list)


def dominates(a: Dict[str, float], b: Dict[str, float],
              objectives: List[Objective]) -> bool:
    """Check if solution a dominates solution b."""
    at_least_one_better = False
    
    for obj in objectives:
        val_a = a.get(obj.name, float('inf'))
        val_b = b.get(obj.name, float('inf'))
        
        if obj.objective_type == ObjectiveType.MAXIMIZE:
            if val_a < val_b:
                return False
            if val_a > val_b:
                at_least_one_better = True
        else:  # MINIMIZE
            if val_a > val_b:
                return False
            if val_a < val_b:
                at_least_one_better = True
    
    return at_least_one_better


def compute_pareto_front(
    candidates: List[CandidateSolution],
    objectives: List[Objective]
) -> List[CandidateSolution]:
    """Compute Pareto front using non-dominated sorting."""
    n = len(candidates)
    
    for i in range(n):
        candidates[i].domination_count = 0
        candidates[i].dominated_solutions = []
        
        for j in range(n):
            if i == j:
                continue
            if dominates(candidates[i].objective_values,
                        candidates[j].objective_values, objectives):
                candidates[i].dominated_solutions.append(candidates[j].compound_id)
            elif dominates(candidates[j].objective_values,
                          candidates[i].objective_values, objectives):
                candidates[i].domination_count += 1
    
    # First front: non-dominated solutions
    pareto_optimal = []
    for c in candidates:
        if c.domination_count == 0:
            c.is_pareto_optimal = True
            c.pareto_rank = 1
            pareto_optimal.append(c)
    
    # Compute crowding distance
    _compute_crowding_distance(pareto_optimal, objectives)
    
    return pareto_optimal


def _compute_crowding_distance(
    solutions: List[CandidateSolution],
    objectives: List[Objective]
):
    """Compute crowding distance for diversity maintenance."""
    n = len(solutions)
    if n <= 2:
        for s in solutions:
            s.crowding_distance = float('inf')
        return
    
    for s in solutions:
        s.crowding_distance = 0.0
    
    for obj in objectives:
        sorted_solutions = sorted(
            solutions,
            key=lambda s: s.objective_values.get(obj.name, 0)
        )
        
        sorted_solutions[0].crowding_distance = float('inf')
        sorted_solutions[-1].crowding_distance = float('inf')
        
        obj_range = (
            sorted_solutions[-1].objective_values.get(obj.name, 0) -
            sorted_solutions[0].objective_values.get(obj.name, 0)
        )
        
        if obj_range == 0:
            continue
        
        for i in range(1, n - 1):
            val_next = sorted_solutions[i + 1].objective_values.get(obj.name, 0)
            val_prev = sorted_solutions[i - 1].objective_values.get(obj.name, 0)
            sorted_solutions[i].crowding_distance += (val_next - val_prev) / obj_range


def compute_hypervolume(
    pareto_front: List[CandidateSolution],
    objectives: List[Objective],
    reference_point: Optional[Dict[str, float]] = None
) -> float:
    """Compute hypervolume indicator for Pareto front quality."""
    if not pareto_front:
        return 0.0
    
    if reference_point is None:
        reference_point = {}
        for obj in objectives:
            values = [s.objective_values.get(obj.name, 0) for s in pareto_front]
            if obj.objective_type == ObjectiveType.MAXIMIZE:
                reference_point[obj.name] = min(values) - 1.0
            else:
                reference_point[obj.name] = max(values) + 1.0
    
    # 2D hypervolume calculation (for visualization)
    if len(objectives) == 2:
        obj1, obj2 = objectives[0].name, objectives[1].name
        points = sorted(
            [(s.objective_values[obj1], s.objective_values[obj2])
             for s in pareto_front]
        )
        
        hv = 0.0
        ref_x = reference_point[obj1]
        ref_y = reference_point[obj2]
        
        for i, (x, y) in enumerate(points):
            if i < len(points) - 1:
                width = points[i + 1][0] - x
            else:
                width = ref_x - x
            height = ref_y - y
            hv += abs(width * height)
        
        return hv
    
    # Approximate for higher dimensions
    n_samples = 10000
    rng = np.random.RandomState(42)
    
    ranges = {}
    for obj in objectives:
        values = [s.objective_values.get(obj.name, 0) for s in pareto_front]
        ranges[obj.name] = (min(values), reference_point[obj.name])
    
    total_volume = 1.0
    for r in ranges.values():
        total_volume *= abs(r[1] - r[0])
    
    dominated_count = 0
    for _ in range(n_samples):
        point = {}
        for obj in objectives:
            r = ranges[obj.name]
            point[obj.name] = rng.uniform(min(r), max(r))
        
        for sol in pareto_front:
            if dominates(sol.objective_values, point, objectives):
                dominated_count += 1
                break
    
    return total_volume * dominated_count / n_samples


def nsga2_optimize(
    initial_population: List[CandidateSolution],
    objectives: List[Objective],
    n_generations: int = 100,
    population_size: int = 100,
    crossover_rate: float = 0.9,
    mutation_rate: float = 0.1,
    seed: int = 42
) -> ParetoFront:
    """
    NSGA-II multi-objective optimization (simplified simulation).
    """
    rng = np.random.RandomState(seed)
    population = initial_population[:population_size]
    
    generation_history = []
    
    for gen in range(n_generations):
        # Non-dominated sorting
        pareto = compute_pareto_front(population, objectives)
        
        # Track history
        hv = compute_hypervolume(pareto, objectives)
        generation_history.append({
            "generation": gen,
            "n_pareto": len(pareto),
            "hypervolume": hv,
            "best_values": {
                obj.name: (
                    max(s.objective_values.get(obj.name, 0) for s in pareto)
                    if obj.objective_type == ObjectiveType.MAXIMIZE
                    else min(s.objective_values.get(obj.name, float('inf')) for s in pareto)
                )
                for obj in objectives
            }
        })
        
        # Simulated evolution (for demonstration)
        new_population = []
        for sol in population:
            new_sol = CandidateSolution(
                compound_id=f"{sol.compound_id}_g{gen}",
                smiles=sol.smiles,
                objective_values={},
                properties=sol.properties.copy(),
            )
            
            for obj in objectives:
                old_val = sol.objective_values.get(obj.name, 0)
                # Slight improvement with noise
                if obj.objective_type == ObjectiveType.MAXIMIZE:
                    new_val = old_val + rng.normal(0.02, 0.1)
                else:
                    new_val = old_val - rng.normal(0.02, 0.1)
                
                new_sol.objective_values[obj.name] = float(new_val)
            
            new_population.append(new_sol)
        
        # Selection: keep best from combined population
        combined = population + new_population
        pareto_all = compute_pareto_front(combined, objectives)
        
        if len(pareto_all) >= population_size:
            pareto_all.sort(key=lambda s: s.crowding_distance, reverse=True)
            population = pareto_all[:population_size]
        else:
            remaining = [s for s in combined if not s.is_pareto_optimal]
            remaining.sort(key=lambda s: s.domination_count)
            population = pareto_all + remaining[:population_size - len(pareto_all)]
    
    # Final Pareto front
    final_pareto = compute_pareto_front(population, objectives)
    final_hv = compute_hypervolume(final_pareto, objectives)
    
    return ParetoFront(
        solutions=final_pareto,
        n_objectives=len(objectives),
        objective_names=[obj.name for obj in objectives],
        hypervolume=final_hv,
        n_pareto_optimal=len(final_pareto),
        generation_history=generation_history,
    )


def generate_lead_optimization_candidates(
    n_candidates: int = 200,
    seed: int = 42
) -> Tuple[List[CandidateSolution], List[Objective]]:
    """Generate synthetic lead optimization dataset."""
    rng = np.random.RandomState(seed)
    
    objectives = [
        Objective("pKi", ObjectiveType.MAXIMIZE, weight=1.0,
                  constraint_min=6.0, unit=""),
        Objective("selectivity", ObjectiveType.MAXIMIZE, weight=0.8,
                  constraint_min=1.0, unit="log"),
        Objective("clearance", ObjectiveType.MINIMIZE, weight=0.6,
                  constraint_max=50.0, unit="mL/min/kg"),
        Objective("hERG_pIC50", ObjectiveType.MINIMIZE, weight=0.7,
                  constraint_max=5.0, unit=""),
        Objective("SA_score", ObjectiveType.MINIMIZE, weight=0.4,
                  constraint_max=6.0, unit=""),
    ]
    
    candidates = []
    for i in range(n_candidates):
        pki = rng.uniform(4, 10)
        # Anti-correlations to make optimization challenging
        selectivity = max(0, pki * 0.3 + rng.normal(0, 0.5))
        clearance = max(1, 80 - pki * 5 + rng.normal(0, 10))
        herg = max(3, pki * 0.5 + rng.normal(0, 0.5))
        sa_score = max(1, rng.uniform(1, 8))
        
        candidates.append(CandidateSolution(
            compound_id=f"OPT_{i+1:04d}",
            smiles=f"C{i}",  # Placeholder
            objective_values={
                "pKi": float(pki),
                "selectivity": float(selectivity),
                "clearance": float(clearance),
                "hERG_pIC50": float(herg),
                "SA_score": float(sa_score),
            },
            properties={
                "MW": float(rng.uniform(300, 550)),
                "LogP": float(rng.uniform(1, 4)),
                "TPSA": float(rng.uniform(50, 130)),
            }
        ))
    
    return candidates, objectives
