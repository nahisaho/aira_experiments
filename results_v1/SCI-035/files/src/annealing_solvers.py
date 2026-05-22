"""
Annealing Solvers Module
- OpenJij Simulated Annealing (SA)
- OpenJij Simulated Quantum Annealing (SQA)
- Reverse Annealing protocol
- D-Wave mock / Ocean SDK interface
"""
from __future__ import annotations
import time
import numpy as np
import openjij as oj
import dimod
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field


@dataclass
class SolverResult:
    solver_name: str
    best_energy: float
    best_sample: dict
    all_energies: List[float]
    elapsed_sec: float
    num_reads: int
    feasible: bool = False
    metadata: dict = field(default_factory=dict)

    @property
    def mean_energy(self) -> float:
        return float(np.mean(self.all_energies))

    @property
    def std_energy(self) -> float:
        return float(np.std(self.all_energies))

    @property
    def success_rate(self) -> float:
        """Fraction of reads that achieved within 1% of best energy."""
        threshold = self.best_energy * 1.01 if self.best_energy < 0 else self.best_energy * 0.99
        return sum(1 for e in self.all_energies if e <= threshold) / len(self.all_energies)

    def to_dict(self) -> dict:
        return {
            "solver": self.solver_name,
            "best_energy": self.best_energy,
            "mean_energy": self.mean_energy,
            "std_energy": self.std_energy,
            "elapsed_sec": self.elapsed_sec,
            "num_reads": self.num_reads,
            "success_rate": self.success_rate,
            "feasible": self.feasible,
            "metadata": self.metadata,
        }


# ------------------------------------------------------------------ #
#  Simulated Annealing via OpenJij                                    #
# ------------------------------------------------------------------ #
class SARunner:
    """OpenJij Simulated Annealing with configurable schedules."""

    def __init__(
        self,
        num_reads: int = 100,
        num_sweeps: int = 1000,
        beta_range: Tuple[float, float] = (0.1, 10.0),
        schedule: Optional[List[List[float]]] = None,
    ):
        self.num_reads = num_reads
        self.num_sweeps = num_sweeps
        self.beta_range = beta_range
        self.schedule = schedule  # custom [(beta, num_sweeps), ...]

    def solve(self, Q: dict, label: str = "SA") -> SolverResult:
        sampler = oj.SASampler()
        t0 = time.perf_counter()
        if self.schedule:
            response = sampler.sample_qubo(
                Q,
                schedule=self.schedule,
                num_reads=self.num_reads,
            )
        else:
            response = sampler.sample_qubo(
                Q,
                num_reads=self.num_reads,
                num_sweeps=self.num_sweeps,
                beta_min=self.beta_range[0],
                beta_max=self.beta_range[1],
            )
        elapsed = time.perf_counter() - t0

        energies = [s.energy for s in response.record]
        best_idx = int(np.argmin(energies))
        best_sample = dict(zip(response.variables, response.record[best_idx].sample))

        return SolverResult(
            solver_name=label,
            best_energy=energies[best_idx],
            best_sample=best_sample,
            all_energies=energies,
            elapsed_sec=elapsed,
            num_reads=self.num_reads,
            metadata={"num_sweeps": self.num_sweeps, "beta_range": list(self.beta_range)},
        )


# ------------------------------------------------------------------ #
#  Simulated Quantum Annealing via OpenJij                            #
# ------------------------------------------------------------------ #
class SQARunner:
    """OpenJij Simulated Quantum Annealing."""

    def __init__(
        self,
        num_reads: int = 100,
        num_sweeps: int = 1000,
        trotter: int = 4,
        beta: float = 5.0,
        schedule: Optional[List[List[float]]] = None,
    ):
        self.num_reads = num_reads
        self.num_sweeps = num_sweeps
        self.trotter = trotter
        self.beta = beta
        self.schedule = schedule

    def solve(self, Q: dict, label: str = "SQA") -> SolverResult:
        sampler = oj.SQASampler()
        t0 = time.perf_counter()
        if self.schedule:
            response = sampler.sample_qubo(
                Q,
                schedule=self.schedule,
                num_reads=self.num_reads,
                trotter=self.trotter,
            )
        else:
            response = sampler.sample_qubo(
                Q,
                num_reads=self.num_reads,
                num_sweeps=self.num_sweeps,
                trotter=self.trotter,
                beta=self.beta,
            )
        elapsed = time.perf_counter() - t0

        energies = [s.energy for s in response.record]
        best_idx = int(np.argmin(energies))
        best_sample = dict(zip(response.variables, response.record[best_idx].sample))

        return SolverResult(
            solver_name=label,
            best_energy=energies[best_idx],
            best_sample=best_sample,
            all_energies=energies,
            elapsed_sec=elapsed,
            num_reads=self.num_reads,
            metadata={"num_sweeps": self.num_sweeps, "trotter": self.trotter, "beta": self.beta},
        )


# ------------------------------------------------------------------ #
#  Reverse Annealing Protocol                                         #
# ------------------------------------------------------------------ #
class ReverseAnnealingRunner:
    """
    Reverse annealing: start from a known good solution, anneal back toward
    s=s_target (increase quantum fluctuations), then forward anneal to s=1.
    Implemented via custom schedule in OpenJij SQA.
    """

    def __init__(
        self,
        initial_solution: Optional[dict] = None,
        s_target: float = 0.3,
        hold_time: int = 100,
        num_reads: int = 50,
        trotter: int = 8,
    ):
        self.initial_solution = initial_solution
        self.s_target = s_target
        self.hold_time = hold_time
        self.num_reads = num_reads
        self.trotter = trotter

    def _build_reverse_schedule(self) -> List[List[float]]:
        """
        Schedule: s goes 1.0 → s_target (hold) → 1.0
        OpenJij SQA uses s ∈ [0,1]: s=1 classical, s=0 full quantum.
        """
        schedule = []
        steps_down = 20
        for i in range(steps_down):
            s = 1.0 - (1.0 - self.s_target) * (i + 1) / steps_down
            schedule.append([round(max(0.001, min(0.999, s)), 4), 5])
        schedule.append([round(max(0.001, min(0.999, self.s_target)), 4), self.hold_time])
        steps_up = 30
        for i in range(steps_up):
            s = self.s_target + (1.0 - self.s_target) * (i + 1) / steps_up
            schedule.append([round(max(0.001, min(0.999, s)), 4), max(1, int(5 * (1 - i / steps_up)))])
        return schedule

    def solve(self, Q: dict, initial_solution: Optional[dict] = None) -> SolverResult:
        sampler = oj.SQASampler()
        schedule = self._build_reverse_schedule()
        init_sol = initial_solution or self.initial_solution

        t0 = time.perf_counter()
        try:
            response = sampler.sample_qubo(
                Q,
                schedule=schedule,
                num_reads=self.num_reads,
                trotter=self.trotter,
                initial_state=init_sol,
            )
        except Exception:
            # Fallback without initial state
            response = sampler.sample_qubo(
                Q,
                schedule=schedule,
                num_reads=self.num_reads,
                trotter=self.trotter,
            )
        elapsed = time.perf_counter() - t0

        energies = [s.energy for s in response.record]
        best_idx = int(np.argmin(energies))
        best_sample = dict(zip(response.variables, response.record[best_idx].sample))

        return SolverResult(
            solver_name="ReverseAnnealing",
            best_energy=energies[best_idx],
            best_sample=best_sample,
            all_energies=energies,
            elapsed_sec=elapsed,
            num_reads=self.num_reads,
            metadata={
                "s_target": self.s_target,
                "hold_time": self.hold_time,
                "schedule_steps": len(schedule),
            },
        )


# ------------------------------------------------------------------ #
#  Schedule Tuning Utility                                            #
# ------------------------------------------------------------------ #
def geometric_beta_schedule(
    beta_min: float, beta_max: float, num_steps: int, sweeps_per_step: int = 10
) -> List[List[float]]:
    """Geometric (exponential) cooling schedule for SA."""
    betas = np.geomspace(beta_min, beta_max, num_steps)
    return [[float(b), sweeps_per_step] for b in betas]


def linear_beta_schedule(
    beta_min: float, beta_max: float, num_steps: int, sweeps_per_step: int = 10
) -> List[List[float]]:
    """Linear cooling schedule for SA."""
    betas = np.linspace(beta_min, beta_max, num_steps)
    return [[float(b), sweeps_per_step] for b in betas]


def parabolic_s_schedule(
    s_min: float = 0.0, s_max: float = 1.0, num_steps: int = 50, sweeps_per_step: int = 10
) -> List[List[float]]:
    """Parabolic annealing parameter schedule for SQA (s ∈ [0,1])."""
    t = np.linspace(0, 1, num_steps)
    # Parabolic: slow start, fast middle, slow end
    s_vals = s_min + (s_max - s_min) * (3 * t**2 - 2 * t**3)
    s_vals = np.clip(s_vals, 0.0, 1.0)
    return [[float(s), sweeps_per_step] for s in s_vals]


def parabolic_gamma_schedule(
    gamma_max: float, gamma_min: float, num_steps: int, sweeps_per_step: int = 10
) -> List[List[float]]:
    """Legacy alias — now maps to s-schedule for SQA (s ∈ [0,1])."""
    return parabolic_s_schedule(0.0, 1.0, num_steps, sweeps_per_step)
