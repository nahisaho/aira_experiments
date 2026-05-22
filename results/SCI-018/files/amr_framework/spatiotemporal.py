from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import FIGURES_DIR, RESULTS_DIR, ROOT, log_event, save_json, seed_everything

try:
    import seaborn as sns
except Exception:  # pragma: no cover
    sns = None


@dataclass
class SpatialPatch:
    name: str
    population: float
    beta: float
    antibiotic_baseline: float
    selection_pressure: float
    mu: float = 0.001
    sigma: float = 0.18
    gamma: float = 0.12
    alpha: float = 0.10
    omega: float = 0.03
    mutation: float = 0.01
    initial_state: tuple[float, float, float, float, float] = (0.96, 0.01, 0.02, 0.005, 0.005)


class SpatiotemporalAMRModel:
    def __init__(self, seed: int = 42) -> None:
        seed_everything(seed)
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.patches = [
            SpatialPatch("hospital", 1500.0, beta=0.48, antibiotic_baseline=0.80, selection_pressure=0.28),
            SpatialPatch("community", 5000.0, beta=0.24, antibiotic_baseline=0.38, selection_pressure=0.12),
            SpatialPatch("livestock", 3000.0, beta=0.31, antibiotic_baseline=0.58, selection_pressure=0.20),
        ]
        self.n_patches = len(self.patches)
        self.migration = np.array(
            [
                [0.00, 0.010, 0.004],
                [0.008, 0.00, 0.006],
                [0.003, 0.005, 0.00],
            ]
        )
        self.last_simulation: pd.DataFrame | None = None

    def _default_schedule(self, T: int) -> dict[str, np.ndarray]:
        time = np.arange(T + 1)
        return {
            "hospital": np.clip(0.80 + 0.08 * np.sin(2 * np.pi * time / 28), 0.55, 0.95),
            "community": np.clip(0.38 + 0.06 * np.sin(2 * np.pi * (time + 6) / 40), 0.20, 0.60),
            "livestock": np.clip(0.58 + 0.10 * np.cos(2 * np.pi * time / 35), 0.30, 0.80),
        }

    def simulate(self, T: int = 120, antibiotic_use_schedule: dict[str, np.ndarray] | None = None) -> pd.DataFrame:
        schedule = antibiotic_use_schedule or self._default_schedule(T)
        states = {}
        for patch in self.patches:
            proportions = np.array(patch.initial_state, dtype=float)
            proportions = proportions / proportions.sum()
            states[patch.name] = proportions * patch.population

        records = []
        for t in range(T + 1):
            for patch in self.patches:
                S, E, I_S, I_R, R = states[patch.name]
                total = max(1.0, S + E + I_S + I_R + R)
                resistance_prevalence = I_R / max(1e-9, I_S + I_R)
                records.append(
                    {
                        "time": t,
                        "patch": patch.name,
                        "S": S,
                        "E": E,
                        "I_S": I_S,
                        "I_R": I_R,
                        "R": R,
                        "N": total,
                        "antibiotic_use": float(schedule[patch.name][t]),
                        "resistance_prevalence": resistance_prevalence,
                    }
                )
            if t == T:
                break
            next_states: dict[str, np.ndarray] = {}
            for patch in self.patches:
                S, E, I_S, I_R, R = states[patch.name]
                N = max(1.0, S + E + I_S + I_R + R)
                antibiotic_use = float(schedule[patch.name][t])
                p_resist = np.clip(0.03 + 0.6 * antibiotic_use + patch.selection_pressure, 0.01, 0.95)
                force = patch.beta * (I_S + I_R) / N
                dS = patch.mu * N - force * S - patch.mu * S + patch.gamma * R
                dE = force * S - patch.sigma * E - patch.mu * E
                dI_S = patch.sigma * E * (1 - p_resist) - patch.gamma * I_S - patch.mu * I_S - patch.alpha * antibiotic_use * I_S
                dI_R = patch.sigma * E * p_resist + patch.mutation * I_S - patch.gamma * I_R - patch.mu * I_R
                dR = patch.gamma * (I_S + I_R) - patch.mu * R - patch.omega * R
                updated = np.array([S + dS, E + dE, I_S + dI_S, I_R + dI_R, R + dR], dtype=float)
                next_states[patch.name] = np.clip(updated, 0.0, None)

            infected_matrix = np.array([
                [next_states[patch.name][1], next_states[patch.name][2], next_states[patch.name][3]]
                for patch in self.patches
            ], dtype=float)
            for i, patch in enumerate(self.patches):
                outgoing = self.migration[i].sum()
                moved = infected_matrix[i] * outgoing
                next_states[patch.name][1:4] -= moved
                for j, target_patch in enumerate(self.patches):
                    if i == j:
                        continue
                    next_states[target_patch.name][1:4] += infected_matrix[i] * self.migration[i, j]
            states = next_states
        self.last_simulation = pd.DataFrame(records)
        return self.last_simulation.copy()

    def compute_resistance_prevalence(self) -> pd.DataFrame:
        if self.last_simulation is None:
            self.simulate()
        assert self.last_simulation is not None
        return self.last_simulation[["time", "patch", "resistance_prevalence", "antibiotic_use"]].copy()

    def sensitivity_analysis(self, parameter: str, values: list[float]) -> pd.DataFrame:
        base_values = {patch.name: getattr(patch, parameter) for patch in self.patches}
        results = []
        for value in values:
            for patch in self.patches:
                setattr(patch, parameter, value)
            simulation = self.simulate(T=80)
            final_resistance = simulation.groupby("patch")["resistance_prevalence"].last().mean()
            results.append({"parameter": parameter, "value": value, "final_resistance_mean": final_resistance})
        for patch in self.patches:
            setattr(patch, parameter, base_values[patch.name])
        return pd.DataFrame(results)

    def visualize(self, simulation: pd.DataFrame) -> tuple[str, str]:
        plt.figure(figsize=(10, 5))
        for patch, subset in simulation.groupby("patch"):
            plt.plot(subset["time"], subset["resistance_prevalence"], label=patch)
        plt.xlabel("Time")
        plt.ylabel("Resistance prevalence")
        plt.title("Resistance dynamics across spatial patches")
        plt.legend()
        plt.grid(alpha=0.3)
        line_path = FIGURES_DIR / "resistance_dynamics.png"
        plt.savefig(line_path, dpi=300, bbox_inches="tight")
        plt.close()

        pivot = simulation.pivot_table(index="patch", columns="time", values="resistance_prevalence")
        plt.figure(figsize=(12, 3.8))
        if sns is not None:
            sns.heatmap(pivot, cmap="viridis")
        else:
            plt.imshow(pivot.values, aspect="auto", cmap="viridis")
            plt.yticks(range(len(pivot.index)), pivot.index)
            plt.colorbar(label="Resistance prevalence")
        plt.title("Spatial heatmap of resistance prevalence")
        plt.xlabel("Time")
        plt.ylabel("Patch")
        heatmap_path = FIGURES_DIR / "spatial_heatmap.png"
        plt.savefig(heatmap_path, dpi=300, bbox_inches="tight")
        plt.close()
        log_event(
            phase="component_5",
            event_type="file_written",
            skill_or_tool="SpatiotemporalAMRModel",
            files_written=[str(line_path.relative_to(ROOT)), str(heatmap_path.relative_to(ROOT))],
        )
        return str(line_path), str(heatmap_path)


def run_component(seed: int = 42) -> dict[str, Any]:
    model = SpatiotemporalAMRModel(seed=seed)
    simulation = model.simulate(T=120)
    sim_path = RESULTS_DIR / "spatiotemporal_simulation.csv"
    simulation.to_csv(sim_path, index=False)
    sensitivity = model.sensitivity_analysis("selection_pressure", [0.08, 0.14, 0.20, 0.26, 0.32])
    sensitivity_path = RESULTS_DIR / "spatiotemporal_sensitivity.csv"
    sensitivity.to_csv(sensitivity_path, index=False)
    model.visualize(simulation)
    summary = {
        "final_resistance_by_patch": simulation.groupby("patch")["resistance_prevalence"].last().round(4).to_dict(),
        "mean_resistance": round(float(simulation.groupby("time")["resistance_prevalence"].mean().mean()), 4),
        "peak_resistance": round(float(simulation["resistance_prevalence"].max()), 4),
        "sensitivity_max_final": round(float(sensitivity["final_resistance_mean"].max()), 4),
    }
    save_json(RESULTS_DIR / "spatiotemporal_summary.json", summary)
    log_event(
        phase="component_5",
        event_type="handoff_completed",
        skill_or_tool="spatiotemporal",
        handoff_out=summary,
        files_written=[str(sim_path.relative_to(ROOT)), str(sensitivity_path.relative_to(ROOT))],
    )
    return summary
