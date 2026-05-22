"""
Benchmark Suite: Problem Scaling & Quantum Advantage Analysis
Runs all solvers across problem sizes and records metrics.
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import time
import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Dict, Tuple
from tqdm import tqdm

from src.qubo_formulation import QUBOFormulator, VRPQUBOFormulator
from src.annealing_solvers import SARunner, SQARunner, ReverseAnnealingRunner
from src.annealing_solvers import geometric_beta_schedule, linear_beta_schedule, parabolic_gamma_schedule
from src.classical_solvers import GreedyLocalSearch, QAOASimulator
from src.minor_embedding import MinorEmbeddingAnalyzer

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)


# ------------------------------------------------------------------ #
#  Random QUBO Generator                                              #
# ------------------------------------------------------------------ #
def generate_random_qubo(n: int, density: float = 0.5, seed: int = 42) -> dict:
    rng = np.random.default_rng(seed)
    Q = {}
    variables = [f"x{i}" for i in range(n)]
    for i, vi in enumerate(variables):
        Q[(vi, vi)] = float(rng.uniform(-1, 1))
        for j, vj in enumerate(variables):
            if j > i and rng.random() < density:
                Q[(vi, vj)] = float(rng.uniform(-1, 1))
    return Q


# ------------------------------------------------------------------ #
#  Schedule Comparison Benchmark                                      #
# ------------------------------------------------------------------ #
def benchmark_schedules(Q: dict, num_reads: int = 50) -> pd.DataFrame:
    """Compare geometric vs linear vs parabolic annealing schedules."""
    rows = []

    schedules = {
        "geometric_fast":  geometric_beta_schedule(0.1, 10.0, 50,  sweeps_per_step=10),
        "geometric_slow":  geometric_beta_schedule(0.1, 20.0, 100, sweeps_per_step=20),
        "linear_fast":     linear_beta_schedule(0.1, 10.0, 50,     sweeps_per_step=10),
        "linear_slow":     linear_beta_schedule(0.1, 20.0, 100,    sweeps_per_step=20),
    }

    for sched_name, sched in schedules.items():
        runner = SARunner(num_reads=num_reads, schedule=sched)
        result = runner.solve(Q, label=f"SA_{sched_name}")
        rows.append({
            "schedule": sched_name,
            "best_energy": result.best_energy,
            "mean_energy": result.mean_energy,
            "std_energy": result.std_energy,
            "elapsed_sec": result.elapsed_sec,
            "success_rate": result.success_rate,
        })

    # SQA schedules
    sqa_schedules = {
        "parabolic_standard": parabolic_gamma_schedule(3.0, 0.01, 50, sweeps_per_step=10),
        "parabolic_aggressive": parabolic_gamma_schedule(5.0, 0.001, 100, sweeps_per_step=20),
    }
    for sched_name, sched in sqa_schedules.items():
        runner = SQARunner(num_reads=num_reads, schedule=sched, beta=5.0)
        result = runner.solve(Q, label=f"SQA_{sched_name}")
        rows.append({
            "schedule": f"SQA_{sched_name}",
            "best_energy": result.best_energy,
            "mean_energy": result.mean_energy,
            "std_energy": result.std_energy,
            "elapsed_sec": result.elapsed_sec,
            "success_rate": result.success_rate,
        })

    return pd.DataFrame(rows)


# ------------------------------------------------------------------ #
#  Solver Comparison Benchmark                                        #
# ------------------------------------------------------------------ #
def benchmark_solvers(Q: dict, num_reads: int = 100) -> pd.DataFrame:
    """Compare all solvers on the same QUBO instance."""
    rows = []

    solvers = [
        SARunner(num_reads=num_reads, num_sweeps=1000),
        SQARunner(num_reads=num_reads, num_sweeps=1000),
        GreedyLocalSearch(num_restarts=num_reads, max_iter=2000),
        QAOASimulator(p_layers=2, num_reads=num_reads, max_exact_n=10),
    ]

    # ReverseAnnealing (needs initial solution from SA)
    sa_result = SARunner(num_reads=20).solve(Q)

    for solver in solvers:
        try:
            result = solver.solve(Q)
            rows.append(result.to_dict())
        except Exception as e:
            rows.append({"solver": str(solver.__class__.__name__), "error": str(e)})

    # Reverse annealing
    try:
        ra = ReverseAnnealingRunner(
            initial_solution=sa_result.best_sample,
            s_target=0.3,
            hold_time=50,
            num_reads=50,
        )
        result = ra.solve(Q, sa_result.best_sample)
        rows.append(result.to_dict())
    except Exception as e:
        rows.append({"solver": "ReverseAnnealing", "error": str(e)})

    return pd.DataFrame(rows)


# ------------------------------------------------------------------ #
#  Scaling Benchmark                                                  #
# ------------------------------------------------------------------ #
def benchmark_scaling(
    sizes: List[int] = [5, 8, 10, 15, 20, 30],
    density: float = 0.5,
    num_reads: int = 50,
) -> pd.DataFrame:
    """Measure time-to-solution vs problem size for each solver."""
    rows = []

    for n in tqdm(sizes, desc="Problem sizes"):
        Q = generate_random_qubo(n, density=density)

        for SolverCls, kwargs, name in [
            (SARunner, {"num_reads": num_reads, "num_sweeps": 500}, "SA"),
            (SQARunner, {"num_reads": min(num_reads, 20), "num_sweeps": 200, "trotter": 4}, "SQA"),
            (GreedyLocalSearch, {"num_restarts": num_reads}, "Greedy"),
        ]:
            try:
                solver = SolverCls(**kwargs)
                result = solver.solve(Q)
                rows.append({
                    "n": n,
                    "solver": name,
                    "best_energy": result.best_energy,
                    "mean_energy": result.mean_energy,
                    "elapsed_sec": result.elapsed_sec,
                    "success_rate": result.success_rate,
                    "density": density,
                })
            except Exception as e:
                rows.append({"n": n, "solver": name, "error": str(e)})

        # QAOA only for small instances
        if n <= 10:
            try:
                qaoa = QAOASimulator(p_layers=2, num_reads=num_reads, max_exact_n=10)
                result = qaoa.solve(Q)
                rows.append({
                    "n": n,
                    "solver": "QAOA(p=2)",
                    "best_energy": result.best_energy,
                    "mean_energy": result.mean_energy,
                    "elapsed_sec": result.elapsed_sec,
                    "success_rate": result.success_rate,
                    "density": density,
                })
            except Exception as e:
                rows.append({"n": n, "solver": "QAOA(p=2)", "error": str(e)})

    return pd.DataFrame(rows)


# ------------------------------------------------------------------ #
#  VRP Case Study                                                     #
# ------------------------------------------------------------------ #
def vrp_case_study(
    city_counts: List[int] = [4, 5, 6],
    num_vehicles: int = 2,
    num_reads: int = 50,
) -> Dict:
    """Full VRP QUBO benchmark with solver comparison."""
    results = {}

    for N in city_counts:
        np.random.seed(N * 7)
        coords = np.random.rand(N, 2) * 100
        dist = np.sqrt(((coords[:, None] - coords[None, :]) ** 2).sum(axis=2))

        formulator = VRPQUBOFormulator(
            num_cities=N,
            num_vehicles=num_vehicles,
            distance_matrix=dist,
        )
        Q, meta = formulator.build_qubo()

        print(f"  VRP N={N}: {meta['num_variables']} vars, {len(Q)} QUBO terms")

        solver_rows = []
        for SolverCls, kwargs, name in [
            (SARunner, {"num_reads": num_reads, "num_sweeps": 1000}, "SA"),
            (SQARunner, {"num_reads": min(num_reads, 20), "num_sweeps": 300, "trotter": 4}, "SQA"),
            (GreedyLocalSearch, {"num_restarts": num_reads}, "Greedy"),
        ]:
            try:
                solver = SolverCls(**kwargs)
                result = solver.solve(Q)
                d = result.to_dict()
                d.update({"N": N, "num_vehicles": num_vehicles})
                solver_rows.append(d)
            except Exception as e:
                solver_rows.append({"solver": name, "N": N, "error": str(e)})

        results[N] = {
            "meta": meta,
            "coords": coords.tolist(),
            "distance_matrix": dist.tolist(),
            "solver_results": solver_rows,
        }

    return results


# ------------------------------------------------------------------ #
#  Embedding Analysis                                                 #
# ------------------------------------------------------------------ #
def embedding_analysis(sizes: List[int] = [10, 20, 30, 50]) -> pd.DataFrame:
    """Compare embedding strategies across problem sizes."""
    analyzer = MinorEmbeddingAnalyzer()
    rows = []
    for n in sizes:
        Q = generate_random_qubo(n, density=0.4)
        results = analyzer.compare_strategies(Q, hardware_size=300)
        for r in results:
            d = r.to_dict()
            d["problem_size"] = n
            d["quality_score"] = analyzer.embedding_quality_score(r)
            rows.append(d)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ #
#  Main                                                               #
# ------------------------------------------------------------------ #
if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.INFO)

    print("=" * 60)
    print("Quantum Annealing Evaluation Framework — Benchmark Suite")
    print("=" * 60)

    # 1. Schedule comparison
    print("\n[1] Annealing schedule comparison...")
    Q_bench = generate_random_qubo(20, density=0.5)
    df_schedules = benchmark_schedules(Q_bench, num_reads=50)
    df_schedules.to_csv(RESULTS_DIR / "schedule_comparison.csv", index=False)
    print(df_schedules.to_string(index=False))

    # 2. Solver comparison
    print("\n[2] Solver comparison (n=15)...")
    Q_15 = generate_random_qubo(15, density=0.5)
    df_solvers = benchmark_solvers(Q_15, num_reads=100)
    df_solvers.to_csv(RESULTS_DIR / "solver_comparison.csv", index=False)
    print(df_solvers[["solver", "best_energy", "mean_energy", "elapsed_sec", "success_rate"]].to_string(index=False))

    # 3. Scaling analysis
    print("\n[3] Scaling analysis...")
    df_scaling = benchmark_scaling(sizes=[5, 8, 10, 15, 20, 30], num_reads=50)
    df_scaling.to_csv(RESULTS_DIR / "scaling_analysis.csv", index=False)
    print(df_scaling.groupby(["n", "solver"])["elapsed_sec"].mean().unstack().to_string())

    # 4. VRP case study
    print("\n[4] VRP case study...")
    vrp_results = vrp_case_study(city_counts=[4, 5, 6], num_vehicles=2, num_reads=50)
    with open(RESULTS_DIR / "vrp_results.json", "w") as f:
        json.dump(vrp_results, f, indent=2)
    for N, res in vrp_results.items():
        print(f"  VRP(N={N}): {len(res['solver_results'])} solver runs")

    # 5. Embedding analysis
    print("\n[5] Minor embedding analysis...")
    df_embedding = embedding_analysis(sizes=[10, 20, 30, 50])
    df_embedding.to_csv(RESULTS_DIR / "embedding_analysis.csv", index=False)
    print(df_embedding[["problem_size", "strategy", "avg_chain_length", "overhead", "quality_score"]].to_string(index=False))

    print("\n✓ All benchmarks complete. Results saved to results/")
