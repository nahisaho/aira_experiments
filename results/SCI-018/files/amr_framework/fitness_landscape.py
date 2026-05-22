from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Any, Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from . import FIGURES_DIR, RESULTS_DIR, ROOT, log_event, save_json, seed_everything

try:
    import seaborn as sns
except Exception:  # pragma: no cover
    sns = None


@dataclass
class LandscapeStatistics:
    roughness: float
    neutrality: float
    peaks_count: int


class NKFitnessLandscape:
    def __init__(self, N: int = 8, K: int = 3, seed: int = 42) -> None:
        if K >= N:
            raise ValueError("K must be smaller than N.")
        seed_everything(seed)
        self.N = N
        self.K = K
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.antibiotics = ["ampicillin", "ciprofloxacin", "gentamicin"]
        self.antibiotic_loci = {
            "ampicillin": [0, 1, 4],
            "ciprofloxacin": [2, 3, 6],
            "gentamicin": [1, 5, 7],
        }
        self.interactions = self._build_interactions()
        self.tables = self._build_tables()
        self.genotypes = [format(idx, f"0{self.N}b") for idx in range(2 ** self.N)]

    def _build_interactions(self) -> dict[int, list[int]]:
        interactions: dict[int, list[int]] = {}
        for locus in range(self.N):
            neighbors = [((locus + step + 1) % self.N) for step in range(self.K)]
            interactions[locus] = [locus, *neighbors]
        return interactions

    def _build_tables(self) -> dict[str, dict[int, np.ndarray]]:
        tables: dict[str, dict[int, np.ndarray]] = {}
        for antibiotic in self.antibiotics:
            tables[antibiotic] = {}
            for locus in range(self.N):
                local_rng = np.random.default_rng(self.seed + 100 * (locus + 1) + len(antibiotic))
                tables[antibiotic][locus] = local_rng.random(2 ** (self.K + 1))
        return tables

    def genotype_to_tuple(self, genotype: str | Iterable[int]) -> tuple[int, ...]:
        if isinstance(genotype, str):
            return tuple(int(bit) for bit in genotype)
        return tuple(int(bit) for bit in genotype)

    def _local_state_index(self, genotype_bits: tuple[int, ...], locus: int) -> int:
        state = 0
        for bit in [genotype_bits[idx] for idx in self.interactions[locus]]:
            state = (state << 1) | bit
        return state

    def compute_fitness(
        self,
        genotype: str | Iterable[int],
        antibiotic: str,
        antibiotic_concentration: float = 0.7,
    ) -> float:
        bits = self.genotype_to_tuple(genotype)
        contributions = []
        for locus in range(self.N):
            state_index = self._local_state_index(bits, locus)
            contributions.append(self.tables[antibiotic][locus][state_index])
        base_fitness = float(np.mean(contributions))
        resistance_loci = self.antibiotic_loci[antibiotic]
        resistance_score = float(np.mean([bits[idx] for idx in resistance_loci]))
        mutational_burden = sum(bits) / self.N
        antibiotic_pressure = antibiotic_concentration * (0.65 * resistance_score - 0.40 * (1 - resistance_score))
        resistance_cost = 0.12 * mutational_burden * (1 - 0.5 * antibiotic_concentration)
        fitness = 0.60 * base_fitness + 0.40 * (0.55 + antibiotic_pressure - resistance_cost)
        return float(np.clip(fitness, 0.0, 1.0))

    def all_fitness(self, antibiotic: str, antibiotic_concentration: float = 0.7) -> dict[str, float]:
        return {
            genotype: self.compute_fitness(genotype, antibiotic, antibiotic_concentration)
            for genotype in self.genotypes
        }

    def get_neighbors(self, genotype: str) -> list[str]:
        neighbors = []
        genotype_list = list(genotype)
        for idx in range(self.N):
            mutated = genotype_list.copy()
            mutated[idx] = "1" if genotype[idx] == "0" else "0"
            neighbors.append("".join(mutated))
        return neighbors

    def get_landscape_statistics(self, antibiotic: str = "ampicillin") -> dict[str, float]:
        analyzer = FitnessLandscapeAnalyzer(self, antibiotic)
        return {
            "roughness": round(analyzer.compute_roughness(), 4),
            "neutrality": round(analyzer.compute_neutrality(), 4),
            "peaks_count": int(len(analyzer.find_local_peaks())),
            "mean_fitness": round(float(np.mean(list(self.all_fitness(antibiotic).values()))), 4),
        }


class FitnessLandscapeAnalyzer:
    def __init__(self, landscape: NKFitnessLandscape, antibiotic: str = "ampicillin") -> None:
        self.landscape = landscape
        self.antibiotic = antibiotic
        self.fitness_cache = landscape.all_fitness(antibiotic)

    def find_local_peaks(self) -> list[str]:
        peaks = []
        for genotype, fitness in self.fitness_cache.items():
            neighbor_fitness = [self.fitness_cache[neighbor] for neighbor in self.landscape.get_neighbors(genotype)]
            if all(fitness >= val for val in neighbor_fitness):
                peaks.append(genotype)
        return peaks

    def compute_roughness(self) -> float:
        roughness_terms = []
        for genotype, fitness in self.fitness_cache.items():
            neighbor_mean = float(np.mean([self.fitness_cache[nbr] for nbr in self.landscape.get_neighbors(genotype)]))
            roughness_terms.append(abs(fitness - neighbor_mean))
        return float(np.mean(roughness_terms))

    def compute_neutrality(self, tolerance: float = 0.01) -> float:
        neutral_pairs = 0
        total_pairs = 0
        for genotype, fitness in self.fitness_cache.items():
            for neighbor in self.landscape.get_neighbors(genotype):
                total_pairs += 1
                if abs(fitness - self.fitness_cache[neighbor]) <= tolerance:
                    neutral_pairs += 1
        return neutral_pairs / max(1, total_pairs)

    def visualize_landscape(self) -> str:
        fig, axes = plt.subplots(1, len(self.landscape.antibiotics), figsize=(16, 4.8), constrained_layout=True)
        for axis, antibiotic in zip(axes, self.landscape.antibiotics):
            values = np.array(list(self.landscape.all_fitness(antibiotic).values())).reshape(16, 16)
            if sns is not None:
                sns.heatmap(values, cmap="viridis", ax=axis, cbar=True)
            else:
                image = axis.imshow(values, cmap="viridis", aspect="auto")
                fig.colorbar(image, ax=axis)
            axis.set_title(f"{antibiotic.title()} landscape")
            axis.set_xlabel("Genotype block")
            axis.set_ylabel("Genotype block")
        output_path = FIGURES_DIR / "fitness_landscape.png"
        fig.savefig(output_path, dpi=300)
        plt.close(fig)
        log_event(
            phase="component_2",
            event_type="file_written",
            skill_or_tool="FitnessLandscapeAnalyzer",
            files_written=[str(output_path.relative_to(ROOT))],
        )
        return str(output_path)


def run_component(seed: int = 42) -> dict[str, Any]:
    landscape = NKFitnessLandscape(seed=seed)
    per_antibiotic = {}
    for antibiotic in landscape.antibiotics:
        analyzer = FitnessLandscapeAnalyzer(landscape, antibiotic)
        per_antibiotic[antibiotic] = {
            "roughness": round(analyzer.compute_roughness(), 4),
            "neutrality": round(analyzer.compute_neutrality(), 4),
            "peaks_count": int(len(analyzer.find_local_peaks())),
            "max_fitness": round(float(max(analyzer.fitness_cache.values())), 4),
            "mean_fitness": round(float(np.mean(list(analyzer.fitness_cache.values()))), 4),
        }
    FitnessLandscapeAnalyzer(landscape).visualize_landscape()
    summary = {
        "N": landscape.N,
        "K": landscape.K,
        "antibiotics": per_antibiotic,
    }
    save_json(RESULTS_DIR / "fitness_landscape_stats.json", summary)
    log_event(
        phase="component_2",
        event_type="handoff_completed",
        skill_or_tool="fitness_landscape",
        handoff_out=summary,
    )
    return summary
