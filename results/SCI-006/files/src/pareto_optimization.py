"""
Module 6: Multi-Objective Pareto Optimization for Lead Optimization

Implements NSGA-II-based multi-objective optimization for lead compound
optimization considering binding affinity, ADMET properties, and synthetic accessibility.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import List, Tuple, Dict
from copy import deepcopy


@dataclass
class Candidate:
    mol_id: str
    objectives: Dict[str, float] = field(default_factory=dict)
    rank: int = 0
    crowding_distance: float = 0.0
    fingerprint: np.ndarray = field(default_factory=lambda: np.array([]))


class NSGA2Optimizer:
    """NSGA-II based multi-objective optimizer for lead optimization."""

    def __init__(self, population_size: int = 100, n_generations: int = 50,
                 mutation_rate: float = 0.1, crossover_rate: float = 0.8):
        self.pop_size = population_size
        self.n_gen = n_generations
        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate
        self.objective_names = ['binding_affinity', 'logP', 'synthetic_accessibility',
                                'selectivity', 'metabolic_stability']

    def evaluate_objectives(self, candidate: Candidate, rng: np.random.RandomState) -> Dict[str, float]:
        """Evaluate multiple drug-like objectives (simulated)."""
        fp = candidate.fingerprint
        objs = {}
        objs['binding_affinity'] = -(5 + np.sum(fp[:20]) * 0.25 + rng.normal(0, 0.3))
        objs['logP'] = 2.0 + np.sum(fp[20:40]) * 0.1 + rng.normal(0, 0.2)
        objs['synthetic_accessibility'] = 2.0 + np.sum(fp[40:60]) * 0.08 + rng.normal(0, 0.15)
        objs['selectivity'] = -(np.sum(fp[60:80]) * 0.15 + rng.normal(0, 0.2))
        objs['metabolic_stability'] = np.sum(fp[80:100]) * 0.12 + rng.normal(0, 0.1)
        return objs

    def dominates(self, obj_a: Dict[str, float], obj_b: Dict[str, float]) -> bool:
        """Check if solution a dominates solution b (all objectives minimized)."""
        at_least_one_better = False
        for key in self.objective_names:
            if obj_a[key] > obj_b[key]:
                return False
            if obj_a[key] < obj_b[key]:
                at_least_one_better = True
        return at_least_one_better

    def fast_non_dominated_sort(self, population: List[Candidate]) -> List[List[int]]:
        """NSGA-II fast non-dominated sorting."""
        n = len(population)
        domination_count = [0] * n
        dominated_set = [[] for _ in range(n)]
        fronts = [[]]

        for i in range(n):
            for j in range(i + 1, n):
                if self.dominates(population[i].objectives, population[j].objectives):
                    dominated_set[i].append(j)
                    domination_count[j] += 1
                elif self.dominates(population[j].objectives, population[i].objectives):
                    dominated_set[j].append(i)
                    domination_count[i] += 1

        for i in range(n):
            if domination_count[i] == 0:
                population[i].rank = 0
                fronts[0].append(i)

        k = 0
        while fronts[k]:
            next_front = []
            for i in fronts[k]:
                for j in dominated_set[i]:
                    domination_count[j] -= 1
                    if domination_count[j] == 0:
                        population[j].rank = k + 1
                        next_front.append(j)
            k += 1
            fronts.append(next_front)

        return [f for f in fronts if f]

    def crowding_distance(self, population: List[Candidate], front: List[int]):
        """Calculate crowding distance for a front."""
        n = len(front)
        if n <= 2:
            for idx in front:
                population[idx].crowding_distance = float('inf')
            return

        for idx in front:
            population[idx].crowding_distance = 0.0

        for obj_name in self.objective_names:
            sorted_front = sorted(front, key=lambda i: population[i].objectives[obj_name])
            obj_range = (population[sorted_front[-1]].objectives[obj_name] -
                         population[sorted_front[0]].objectives[obj_name])
            if obj_range == 0:
                continue

            population[sorted_front[0]].crowding_distance = float('inf')
            population[sorted_front[-1]].crowding_distance = float('inf')

            for k in range(1, n - 1):
                population[sorted_front[k]].crowding_distance += (
                    (population[sorted_front[k + 1]].objectives[obj_name] -
                     population[sorted_front[k - 1]].objectives[obj_name]) / obj_range
                )

    def mutate(self, fp: np.ndarray, rng: np.random.RandomState) -> np.ndarray:
        new_fp = fp.copy()
        mask = rng.random(len(fp)) < self.mutation_rate
        new_fp[mask] = 1 - new_fp[mask]
        return new_fp

    def crossover(self, fp1: np.ndarray, fp2: np.ndarray,
                  rng: np.random.RandomState) -> Tuple[np.ndarray, np.ndarray]:
        point = rng.randint(1, len(fp1))
        child1 = np.concatenate([fp1[:point], fp2[point:]])
        child2 = np.concatenate([fp2[:point], fp1[point:]])
        return child1, child2

    def optimize(self, seed: int = 42) -> Tuple[List[Candidate], List[dict]]:
        """Run NSGA-II optimization."""
        rng = np.random.RandomState(seed)
        fp_size = 128

        # Initialize population
        population = []
        for i in range(self.pop_size):
            fp = (rng.random(fp_size) > 0.5).astype(float)
            cand = Candidate(mol_id=f"GEN0-{i:04d}", fingerprint=fp)
            cand.objectives = self.evaluate_objectives(cand, rng)
            population.append(cand)

        history = []

        for gen in range(self.n_gen):
            # Create offspring
            offspring = []
            while len(offspring) < self.pop_size:
                i, j = rng.randint(0, len(population), 2)
                if rng.random() < self.crossover_rate:
                    fp1, fp2 = self.crossover(population[i].fingerprint,
                                               population[j].fingerprint, rng)
                else:
                    fp1, fp2 = population[i].fingerprint.copy(), population[j].fingerprint.copy()

                fp1 = self.mutate(fp1, rng)
                fp2 = self.mutate(fp2, rng)

                for fp in [fp1, fp2]:
                    cand = Candidate(mol_id=f"GEN{gen + 1}-{len(offspring):04d}", fingerprint=fp)
                    cand.objectives = self.evaluate_objectives(cand, rng)
                    offspring.append(cand)

            # Combine and select
            combined = population + offspring[:self.pop_size]
            fronts = self.fast_non_dominated_sort(combined)

            new_pop = []
            for front in fronts:
                if len(new_pop) + len(front) <= self.pop_size:
                    new_pop.extend([combined[i] for i in front])
                else:
                    self.crowding_distance(combined, front)
                    remaining = sorted(front, key=lambda i: combined[i].crowding_distance, reverse=True)
                    needed = self.pop_size - len(new_pop)
                    new_pop.extend([combined[i] for i in remaining[:needed]])
                    break

            population = new_pop

            # Record history
            pareto_front = [p for p in population if p.rank == 0]
            best_affinity = min(p.objectives['binding_affinity'] for p in population)
            mean_affinity = np.mean([p.objectives['binding_affinity'] for p in population])

            history.append({
                'generation': gen,
                'pareto_size': len(pareto_front),
                'best_affinity': best_affinity,
                'mean_affinity': mean_affinity,
            })

            if (gen + 1) % 10 == 0:
                print(f"  Gen {gen + 1}/{self.n_gen}: Pareto size={len(pareto_front)}, "
                      f"Best affinity={-best_affinity:.2f}")

        return population, history


def run_pareto_optimization(output_dir: str = "figures"):
    """Run multi-objective optimization and generate Pareto front figures."""
    print("=" * 60)
    print("Module 6: Multi-Objective Pareto Optimization")
    print("=" * 60)

    optimizer = NSGA2Optimizer(population_size=100, n_generations=50,
                                mutation_rate=0.1, crossover_rate=0.8)

    population, history = optimizer.optimize(seed=42)

    # Extract Pareto front
    pareto = [p for p in population if p.rank == 0]
    print(f"\nFinal Pareto front size: {len(pareto)}")

    # Top candidates
    best_by_affinity = min(population, key=lambda p: p.objectives['binding_affinity'])
    print(f"Best binding affinity: {-best_by_affinity.objectives['binding_affinity']:.2f}")

    # Figure 7: Pareto optimization results
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 2D Pareto front (affinity vs logP)
    ax = axes[0, 0]
    non_pareto = [p for p in population if p.rank > 0]
    if non_pareto:
        ax.scatter([-p.objectives['binding_affinity'] for p in non_pareto],
                   [p.objectives['logP'] for p in non_pareto],
                   c='lightgray', s=20, alpha=0.5, label='Dominated')
    ax.scatter([-p.objectives['binding_affinity'] for p in pareto],
               [p.objectives['logP'] for p in pareto],
               c='red', s=40, zorder=5, label='Pareto front', edgecolors='darkred', linewidths=0.5)
    ax.set_xlabel('Binding Affinity (pKd)')
    ax.set_ylabel('logP')
    ax.set_title('Pareto Front: Affinity vs lipophilicity')
    ax.legend()

    # Affinity vs Synthetic Accessibility
    ax = axes[0, 1]
    if non_pareto:
        ax.scatter([-p.objectives['binding_affinity'] for p in non_pareto],
                   [p.objectives['synthetic_accessibility'] for p in non_pareto],
                   c='lightgray', s=20, alpha=0.5, label='Dominated')
    ax.scatter([-p.objectives['binding_affinity'] for p in pareto],
               [p.objectives['synthetic_accessibility'] for p in pareto],
               c='blue', s=40, zorder=5, label='Pareto front', edgecolors='darkblue', linewidths=0.5)
    ax.set_xlabel('Binding Affinity (pKd)')
    ax.set_ylabel('Synthetic Accessibility Score')
    ax.set_title('Pareto Front: Affinity vs SA')
    ax.legend()

    # Optimization history
    ax = axes[1, 0]
    gens = [h['generation'] for h in history]
    ax.plot(gens, [-h['best_affinity'] for h in history], 'b-', label='Best affinity', linewidth=2)
    ax.plot(gens, [-h['mean_affinity'] for h in history], 'r--', label='Mean affinity', alpha=0.7)
    ax.set_xlabel('Generation')
    ax.set_ylabel('Binding Affinity (pKd)')
    ax.set_title('Optimization Progress')
    ax.legend()

    # Pareto front size evolution
    ax = axes[1, 1]
    ax.plot(gens, [h['pareto_size'] for h in history], 'g-', linewidth=2)
    ax.set_xlabel('Generation')
    ax.set_ylabel('Pareto Front Size')
    ax.set_title('Pareto Front Size Over Generations')
    ax.fill_between(gens, [h['pareto_size'] for h in history], alpha=0.2, color='green')

    plt.suptitle('NSGA-II Multi-Objective Lead Optimization', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pareto_optimization.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Radar chart for top candidates
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    categories = ['Affinity', 'LogP\n(lower=better)', 'SA Score\n(lower=better)',
                   'Selectivity', 'Met. Stability']
    n_cats = len(categories)
    angles = np.linspace(0, 2 * np.pi, n_cats, endpoint=False).tolist()
    angles += angles[:1]

    top_3 = sorted(pareto, key=lambda p: p.objectives['binding_affinity'])[:3]
    colors = ['red', 'blue', 'green']

    for idx, cand in enumerate(top_3):
        values = [
            -cand.objectives['binding_affinity'] / 12,
            1 - cand.objectives['logP'] / 5,
            1 - cand.objectives['synthetic_accessibility'] / 5,
            -cand.objectives['selectivity'] / 5,
            cand.objectives['metabolic_stability'] / 5,
        ]
        values += values[:1]
        ax.plot(angles, values, 'o-', color=colors[idx], linewidth=2,
                label=f'Candidate {idx + 1}', markersize=4)
        ax.fill(angles, values, alpha=0.1, color=colors[idx])

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=9)
    ax.set_ylim(0, 1)
    ax.set_title('Top 3 Pareto-Optimal Candidates\n(normalized properties)', pad=20)
    ax.legend(loc='lower right', bbox_to_anchor=(1.2, -0.05))
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pareto_radar.png', dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\nFigures saved to {output_dir}/")

    return {
        'pareto_size': len(pareto),
        'best_affinity': -best_by_affinity.objectives['binding_affinity'],
        'history': history,
    }


if __name__ == '__main__':
    run_pareto_optimization()
