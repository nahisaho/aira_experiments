from __future__ import annotations

import itertools
import math
from collections import Counter
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from . import FIGURES_DIR, RESULTS_DIR, ROOT, log_event, save_json, seed_everything
from .fitness_landscape import NKFitnessLandscape


class EvolutionaryPathFinder:
    def __init__(self, landscape: NKFitnessLandscape, antibiotic: str = "ampicillin", seed: int = 42) -> None:
        self.landscape = landscape
        self.antibiotic = antibiotic
        self.rng = np.random.default_rng(seed)
        seed_everything(seed)

    def _fitness(self, genotype: str) -> float:
        return self.landscape.compute_fitness(genotype, self.antibiotic)

    def _neighbors_towards(self, current: str, end: str) -> list[str]:
        neighbors = []
        for idx, (bit_cur, bit_end) in enumerate(zip(current, end)):
            if bit_cur != bit_end:
                mutated = list(current)
                mutated[idx] = bit_end
                neighbors.append("".join(mutated))
        return neighbors

    def enumerate_accessible_paths(self, start: str, end: str, landscape: NKFitnessLandscape | None = None) -> list[list[str]]:
        _ = landscape or self.landscape
        paths: list[list[str]] = []

        def dfs(current: str, visited: list[str]) -> None:
            if current == end:
                paths.append(visited.copy())
                return
            current_fitness = self._fitness(current)
            next_nodes = []
            for neighbor in self._neighbors_towards(current, end):
                if self._fitness(neighbor) > current_fitness + 1e-9:
                    next_nodes.append(neighbor)
            for neighbor in sorted(next_nodes, key=self._fitness):
                dfs(neighbor, [*visited, neighbor])

        dfs(start, [start])
        return paths

    def greedy_adaptive_walk(self, start: str, landscape: NKFitnessLandscape | None = None) -> dict[str, Any]:
        _ = landscape or self.landscape
        current = start
        path = [current]
        while True:
            current_fitness = self._fitness(current)
            better_neighbors = [nbr for nbr in self.landscape.get_neighbors(current) if self._fitness(nbr) > current_fitness + 1e-9]
            if not better_neighbors:
                break
            current = max(better_neighbors, key=self._fitness)
            path.append(current)
        return {
            "path": path,
            "endpoint": current,
            "fitness": [round(self._fitness(node), 4) for node in path],
        }

    def monte_carlo_walks(self, n_walks: int = 200) -> dict[str, Any]:
        endpoints = []
        lengths = []
        for _ in range(n_walks):
            start = format(int(self.rng.integers(0, 2 ** self.landscape.N)), f"0{self.landscape.N}b")
            current = start
            length = 0
            while True:
                current_fitness = self._fitness(current)
                better_neighbors = [nbr for nbr in self.landscape.get_neighbors(current) if self._fitness(nbr) > current_fitness + 1e-9]
                if not better_neighbors:
                    break
                gains = np.array([self._fitness(nbr) - current_fitness for nbr in better_neighbors], dtype=float)
                probabilities = gains / gains.sum()
                current = str(self.rng.choice(better_neighbors, p=probabilities))
                length += 1
            endpoints.append(current)
            lengths.append(length)
        endpoint_counts = Counter(endpoints)
        return {
            "n_walks": n_walks,
            "mean_length": round(float(np.mean(lengths)), 4),
            "std_length": round(float(np.std(lengths)), 4),
            "endpoint_counts": dict(endpoint_counts.most_common(8)),
            "convergence": round(max(endpoint_counts.values()) / max(1, n_walks), 4),
        }

    def compute_path_accessibility(self) -> dict[str, Any]:
        start = "0" * self.landscape.N
        end = "1" * self.landscape.N
        accessible_paths = self.enumerate_accessible_paths(start, end)
        total_direct_paths = math.factorial(sum(s != e for s, e in zip(start, end)))
        return {
            "start": start,
            "end": end,
            "accessible_path_count": len(accessible_paths),
            "total_direct_paths": total_direct_paths,
            "accessibility_ratio": round(len(accessible_paths) / max(1, total_direct_paths), 6),
        }

    def predict_evolutionary_trajectory(self, population_size: int = 400, n_generations: int = 60) -> dict[str, Any]:
        genotypes = [format(idx, f"0{self.landscape.N}b") for idx in range(2 ** self.landscape.N)]
        frequencies = np.zeros(len(genotypes), dtype=float)
        frequencies[0] = 1.0
        fitness = np.array([self._fitness(genotype) for genotype in genotypes], dtype=float)
        mu = 0.01
        mutation_matrix = np.zeros((len(genotypes), len(genotypes)), dtype=float)
        for i, genotype_i in enumerate(genotypes):
            for j, genotype_j in enumerate(genotypes):
                distance = sum(a != b for a, b in zip(genotype_i, genotype_j))
                mutation_matrix[i, j] = (mu ** distance) * ((1 - mu) ** (self.landscape.N - distance))
            mutation_matrix[i] /= mutation_matrix[i].sum()

        history = []
        for generation in range(n_generations + 1):
            top_indices = np.argsort(frequencies)[-5:][::-1]
            history.append(
                {
                    "generation": generation,
                    "dominant_genotypes": {genotypes[idx]: round(float(frequencies[idx]), 4) for idx in top_indices if frequencies[idx] > 0.001},
                }
            )
            selected = frequencies * fitness
            selected /= selected.sum()
            post_mutation = selected @ mutation_matrix
            counts = self.rng.multinomial(population_size, post_mutation)
            frequencies = counts / population_size
        positive_mask = frequencies > 0
        diversity = float(-np.sum(frequencies[positive_mask] * np.log2(frequencies[positive_mask])))
        dominant = genotypes[int(np.argmax(frequencies))]
        return {
            "population_size": population_size,
            "n_generations": n_generations,
            "final_diversity": round(diversity, 4),
            "dominant_genotype": dominant,
            "dominant_frequency": round(float(frequencies.max()), 4),
            "history": history,
        }

    def visualize_paths(self, greedy_result: dict[str, Any], trajectory: dict[str, Any]) -> str:
        fig, axes = plt.subplots(2, 1, figsize=(10, 8), constrained_layout=True)
        axes[0].plot(range(len(greedy_result["fitness"])), greedy_result["fitness"], marker="o", color="tab:blue")
        axes[0].set_title("Greedy adaptive walk")
        axes[0].set_xlabel("Step")
        axes[0].set_ylabel("Fitness")
        axes[0].grid(alpha=0.25)

        generations = [entry["generation"] for entry in trajectory["history"]]
        tracked = Counter()
        for entry in trajectory["history"]:
            tracked.update(entry["dominant_genotypes"].keys())
        top_genotypes = [item for item, _count in tracked.most_common(5)]
        for genotype in top_genotypes:
            values = [entry["dominant_genotypes"].get(genotype, 0.0) for entry in trajectory["history"]]
            axes[1].plot(generations, values, marker=".", label=genotype)
        axes[1].set_title("Wright-Fisher genotype frequencies")
        axes[1].set_xlabel("Generation")
        axes[1].set_ylabel("Frequency")
        axes[1].legend(fontsize=8, ncol=2)
        axes[1].grid(alpha=0.25)

        output_path = FIGURES_DIR / "evolutionary_paths.png"
        fig.savefig(output_path, dpi=300)
        plt.close(fig)
        log_event(
            phase="component_3",
            event_type="file_written",
            skill_or_tool="EvolutionaryPathFinder",
            files_written=[str(output_path.relative_to(ROOT))],
        )
        return str(output_path)


def run_component(seed: int = 42) -> dict[str, Any]:
    landscape = NKFitnessLandscape(seed=seed)
    finder = EvolutionaryPathFinder(landscape, antibiotic="ampicillin", seed=seed)
    accessibility = finder.compute_path_accessibility()
    greedy = finder.greedy_adaptive_walk("00000000")
    monte_carlo = finder.monte_carlo_walks(n_walks=300)
    trajectory = finder.predict_evolutionary_trajectory(population_size=500, n_generations=80)
    finder.visualize_paths(greedy, trajectory)
    summary = {
        "accessibility": accessibility,
        "greedy_walk": {
            "endpoint": greedy["endpoint"],
            "path_length": len(greedy["path"]) - 1,
            "final_fitness": greedy["fitness"][-1],
        },
        "monte_carlo": monte_carlo,
        "trajectory": {
            "dominant_genotype": trajectory["dominant_genotype"],
            "dominant_frequency": trajectory["dominant_frequency"],
            "final_diversity": trajectory["final_diversity"],
        },
    }
    save_json(RESULTS_DIR / "evolutionary_paths.json", summary)
    log_event(
        phase="component_3",
        event_type="handoff_completed",
        skill_or_tool="evolutionary_paths",
        handoff_out=summary,
    )
    return summary
