from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import FIGURES_DIR, RESULTS_DIR, ROOT, log_event, save_json, seed_everything


@dataclass
class TreatmentStrategy:
    strategy_type: str
    drugs: list[str]
    doses: list[float]
    schedule: list[float] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class CombinationTherapy(TreatmentStrategy):
    def loewe_synergy(self, drug_profiles: dict[str, dict[str, float]]) -> float:
        potency = np.mean([drug_profiles[drug]["potency"] * dose for drug, dose in zip(self.drugs, self.doses)])
        diversity_bonus = 0.08 * len(set(self.drugs))
        overlap_penalty = 0.04 * max(0, len(self.drugs) - len(set(self.drugs)))
        return float(np.clip(0.85 + 0.25 * potency + diversity_bonus - overlap_penalty, 0.7, 1.5))


@dataclass
class AntibioticCycling(TreatmentStrategy):
    cycle_period: int = 7


class TreatmentOptimizer:
    def __init__(self, seed: int = 42) -> None:
        seed_everything(seed)
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.drug_profiles = {
            "ampicillin": {"potency": 0.72, "liability": 0.62, "toxicity": 0.12},
            "ciprofloxacin": {"potency": 0.81, "liability": 0.54, "toxicity": 0.19},
            "gentamicin": {"potency": 0.76, "liability": 0.48, "toxicity": 0.16},
        }
        self.archive: list[dict[str, Any]] = []

    def _build_strategy(self, candidate: dict[str, Any]) -> TreatmentStrategy:
        if candidate["strategy_type"] == "combination":
            return CombinationTherapy(**candidate)
        if candidate["strategy_type"] == "cycling":
            return AntibioticCycling(**candidate)
        return TreatmentStrategy(**candidate)

    def random_candidate(self) -> dict[str, Any]:
        strategy_type = str(self.rng.choice(["monotherapy", "combination", "cycling", "mixing"], p=[0.25, 0.30, 0.25, 0.20]))
        all_drugs = list(self.drug_profiles)
        if strategy_type == "monotherapy":
            drug = str(self.rng.choice(all_drugs))
            return {"strategy_type": strategy_type, "drugs": [drug], "doses": [round(float(self.rng.uniform(0.6, 1.0)), 3)], "schedule": [1.0]}
        if strategy_type == "combination":
            drugs = list(self.rng.choice(all_drugs, size=2, replace=False))
            doses = [round(float(x), 3) for x in self.rng.uniform(0.4, 0.9, size=2)]
            return {"strategy_type": strategy_type, "drugs": drugs, "doses": doses, "schedule": [1.0, 1.0]}
        if strategy_type == "cycling":
            drugs = list(self.rng.choice(all_drugs, size=3, replace=False))
            doses = [round(float(x), 3) for x in self.rng.uniform(0.45, 0.8, size=3)]
            return {"strategy_type": strategy_type, "drugs": drugs, "doses": doses, "schedule": [1.0, 1.0, 1.0], "cycle_period": int(self.rng.integers(3, 15))}
        drugs = list(self.rng.choice(all_drugs, size=2, replace=False))
        doses = [round(float(x), 3) for x in self.rng.uniform(0.35, 0.75, size=2)]
        schedule = [round(float(x), 3) for x in self.rng.uniform(0.4, 0.8, size=2)]
        return {"strategy_type": strategy_type, "drugs": drugs, "doses": doses, "schedule": schedule}

    def evaluate_strategy(self, strategy: TreatmentStrategy) -> dict[str, float]:
        potencies = np.array([self.drug_profiles[drug]["potency"] for drug in strategy.drugs], dtype=float)
        liabilities = np.array([self.drug_profiles[drug]["liability"] for drug in strategy.drugs], dtype=float)
        toxicities = np.array([self.drug_profiles[drug]["toxicity"] for drug in strategy.drugs], dtype=float)
        doses = np.array(strategy.doses, dtype=float)
        schedule = np.array(strategy.schedule if strategy.schedule else [1.0] * len(strategy.drugs), dtype=float)
        exposure = float(np.mean(doses * schedule))
        efficacy = 0.45 + 0.45 * float(np.mean(potencies * doses)) - 0.10 * float(np.mean(toxicities * doses))
        resistance_rate = 0.18 + 0.42 * float(np.mean(liabilities * doses))

        if isinstance(strategy, CombinationTherapy):
            synergy = strategy.loewe_synergy(self.drug_profiles)
            efficacy += 0.10 * (synergy - 1.0)
            resistance_rate *= 0.86
        elif isinstance(strategy, AntibioticCycling):
            period = max(2, strategy.cycle_period)
            diversity = len(set(strategy.drugs)) / len(self.drug_profiles)
            efficacy += 0.03 * diversity - 0.02 * abs(period - 7) / 7
            resistance_rate *= 0.78 + 0.03 * abs(period - 7) / 7
        elif strategy.strategy_type == "mixing":
            heterogeneity = float(np.std(schedule))
            efficacy += 0.02 * len(set(strategy.drugs))
            resistance_rate *= 0.82 + 0.08 * heterogeneity

        efficacy = float(np.clip(efficacy, 0.0, 0.99))
        resistance_rate = float(np.clip(resistance_rate, 0.01, 0.99))
        return {"efficacy": efficacy, "resistance_rate": resistance_rate, "objective": resistance_rate + (1 - efficacy)}

    def objective_function(self, strategy: TreatmentStrategy) -> float:
        return self.evaluate_strategy(strategy)["objective"]

    def _candidate_with_metrics(self, candidate: dict[str, Any]) -> dict[str, Any]:
        strategy = self._build_strategy(candidate)
        metrics = self.evaluate_strategy(strategy)
        combined = {**strategy.to_dict(), **metrics}
        self.archive.append(combined)
        return combined

    def _crossover(self, parent_a: dict[str, Any], parent_b: dict[str, Any]) -> dict[str, Any]:
        if parent_a["strategy_type"] != parent_b["strategy_type"]:
            return self.random_candidate()
        child = {key: parent_a[key] for key in parent_a if key not in {"efficacy", "resistance_rate", "objective"}}
        if "doses" in child:
            child["doses"] = [round(float((a + b) / 2), 3) for a, b in zip(parent_a["doses"], parent_b["doses"])]
        if "schedule" in child:
            child["schedule"] = [round(float((a + b) / 2), 3) for a, b in zip(parent_a.get("schedule", []), parent_b.get("schedule", []))]
        if child["strategy_type"] == "cycling":
            child["cycle_period"] = int(round((parent_a["cycle_period"] + parent_b["cycle_period"]) / 2))
        return child

    def _mutate(self, candidate: dict[str, Any], rate: float = 0.25) -> dict[str, Any]:
        mutated = {key: value if not isinstance(value, list) else value.copy() for key, value in candidate.items()}
        if self.rng.random() < rate and "doses" in mutated:
            mutated["doses"] = [round(float(np.clip(dose + self.rng.normal(0, 0.08), 0.2, 1.0)), 3) for dose in mutated["doses"]]
        if self.rng.random() < rate and "schedule" in mutated:
            mutated["schedule"] = [round(float(np.clip(val + self.rng.normal(0, 0.08), 0.2, 1.0)), 3) for val in mutated["schedule"]]
        if mutated["strategy_type"] == "cycling" and self.rng.random() < rate:
            mutated["cycle_period"] = int(np.clip(mutated["cycle_period"] + self.rng.integers(-2, 3), 2, 16))
        if self.rng.random() < rate * 0.4:
            replacement = self.random_candidate()
            if replacement["strategy_type"] == mutated["strategy_type"]:
                mutated["drugs"] = replacement["drugs"]
        return mutated

    def genetic_algorithm_optimize(self, population_size: int = 40, n_generations: int = 30) -> dict[str, Any]:
        population = [self.random_candidate() for _ in range(population_size)]
        history = []
        best_solution = None
        for generation in range(n_generations):
            scored = [self._candidate_with_metrics(candidate) for candidate in population]
            scored.sort(key=lambda item: item["objective"])
            elite = scored[: max(2, population_size // 4)]
            best_solution = elite[0]
            history.append({"generation": generation, "best_objective": round(float(elite[0]["objective"]), 4), "mean_objective": round(float(np.mean([item["objective"] for item in scored])), 4)})
            next_population = [
                {key: value for key, value in elite_member.items() if key not in {"efficacy", "resistance_rate", "objective"}}
                for elite_member in elite
            ]
            while len(next_population) < population_size:
                parents = list(self.rng.choice(elite, size=2, replace=True))
                child = self._crossover(parents[0], parents[1])
                child = self._mutate(child)
                next_population.append(child)
            population = next_population[:population_size]
        assert best_solution is not None
        return {"best_strategy": best_solution, "history": history}

    def pareto_front_analysis(self) -> list[dict[str, Any]]:
        pareto_front = []
        for candidate in self.archive:
            dominated = False
            for other in self.archive:
                if other is candidate:
                    continue
                if other["resistance_rate"] <= candidate["resistance_rate"] and other["efficacy"] >= candidate["efficacy"] and (
                    other["resistance_rate"] < candidate["resistance_rate"] or other["efficacy"] > candidate["efficacy"]
                ):
                    dominated = True
                    break
            if not dominated:
                pareto_front.append(candidate)
        unique = {(tuple(item["drugs"]), tuple(item["doses"]), item["strategy_type"], item.get("cycle_period", 0)): item for item in pareto_front}
        return sorted(unique.values(), key=lambda item: (item["resistance_rate"], -item["efficacy"]))

    def compare_strategies(self) -> pd.DataFrame:
        strategies = [
            TreatmentStrategy("monotherapy", ["ampicillin"], [0.9], [1.0]),
            TreatmentStrategy("monotherapy", ["ciprofloxacin"], [0.85], [1.0]),
            CombinationTherapy("combination", ["ampicillin", "gentamicin"], [0.65, 0.70], [1.0, 1.0]),
            AntibioticCycling("cycling", ["ampicillin", "ciprofloxacin", "gentamicin"], [0.60, 0.65, 0.60], [1.0, 1.0, 1.0], cycle_period=7),
            TreatmentStrategy("mixing", ["ciprofloxacin", "gentamicin"], [0.55, 0.55], [0.55, 0.70]),
        ]
        rows = []
        for strategy in strategies:
            metrics = self.evaluate_strategy(strategy)
            rows.append({"label": f"{strategy.strategy_type}:{'+'.join(strategy.drugs)}", **metrics})
        return pd.DataFrame(rows)

    def visualize(self, comparison_df: pd.DataFrame, pareto_front: list[dict[str, Any]]) -> tuple[str, str]:
        archive_df = pd.DataFrame(self.archive)
        plt.figure(figsize=(7, 5))
        plt.scatter(archive_df["resistance_rate"], archive_df["efficacy"], alpha=0.35, c="tab:blue", label="Candidates")
        pareto_df = pd.DataFrame(pareto_front)
        if not pareto_df.empty:
            plt.scatter(pareto_df["resistance_rate"], pareto_df["efficacy"], c="tab:red", label="Pareto front")
            pareto_df = pareto_df.sort_values("resistance_rate")
            plt.plot(pareto_df["resistance_rate"], pareto_df["efficacy"], c="tab:red", alpha=0.7)
        plt.xlabel("Resistance rate")
        plt.ylabel("Clinical efficacy")
        plt.title("Treatment optimization Pareto front")
        plt.legend()
        plt.grid(alpha=0.25)
        pareto_path = FIGURES_DIR / "pareto_front.png"
        plt.savefig(pareto_path, dpi=300, bbox_inches="tight")
        plt.close()

        plt.figure(figsize=(10, 5))
        x = np.arange(len(comparison_df))
        width = 0.36
        plt.bar(x - width / 2, comparison_df["efficacy"], width=width, label="Efficacy", color="tab:green")
        plt.bar(x + width / 2, comparison_df["resistance_rate"], width=width, label="Resistance rate", color="tab:orange")
        plt.xticks(x, comparison_df["label"], rotation=30, ha="right")
        plt.ylabel("Score")
        plt.title("Strategy comparison")
        plt.legend()
        plt.grid(axis="y", alpha=0.25)
        comparison_path = FIGURES_DIR / "strategy_comparison.png"
        plt.savefig(comparison_path, dpi=300, bbox_inches="tight")
        plt.close()
        log_event(
            phase="component_6",
            event_type="file_written",
            skill_or_tool="TreatmentOptimizer",
            files_written=[str(pareto_path.relative_to(ROOT)), str(comparison_path.relative_to(ROOT))],
        )
        return str(pareto_path), str(comparison_path)


def run_component(seed: int = 42) -> dict[str, Any]:
    optimizer = TreatmentOptimizer(seed=seed)
    ga_result = optimizer.genetic_algorithm_optimize(population_size=48, n_generations=32)
    pareto_front = optimizer.pareto_front_analysis()
    comparison_df = optimizer.compare_strategies()
    comparison_path = RESULTS_DIR / "treatment_strategy_comparison.csv"
    comparison_df.to_csv(comparison_path, index=False)
    optimizer.visualize(comparison_df, pareto_front)
    summary = {
        "best_strategy": ga_result["best_strategy"],
        "pareto_front_size": len(pareto_front),
        "ga_history": ga_result["history"],
        "comparison": comparison_df.to_dict(orient="records"),
    }
    save_json(RESULTS_DIR / "treatment_optimization.json", summary)
    log_event(
        phase="component_6",
        event_type="handoff_completed",
        skill_or_tool="treatment_optimizer",
        handoff_out={"best_strategy_type": ga_result["best_strategy"]["strategy_type"], "pareto_front_size": len(pareto_front)},
        files_written=[str(comparison_path.relative_to(ROOT))],
    )
    return summary
