#!/usr/bin/env python3
"""
Main evaluation runner for the Quantum Annealing Performance Framework.
Executes all benchmarks and generates figures + report.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import json
import time
import datetime
import numpy as np
import pandas as pd
from pathlib import Path

from src.benchmark_suite import (
    generate_random_qubo,
    benchmark_schedules,
    benchmark_solvers,
    benchmark_scaling,
    vrp_case_study,
    embedding_analysis,
)
from src.visualization import (
    plot_solver_comparison,
    plot_scaling_analysis,
    plot_schedule_comparison,
    plot_embedding_analysis,
    plot_vrp_routes,
    plot_reverse_annealing_schedule,
    plot_qubo_distribution,
)
from src.qubo_formulation import VRPQUBOFormulator, QUBOFormulator

RESULTS_DIR = Path("results")
FIGURES_DIR = Path("figures")
LOGS_DIR = Path("logs")
for d in [RESULTS_DIR, FIGURES_DIR, LOGS_DIR]:
    d.mkdir(exist_ok=True)

LOG_FILE = LOGS_DIR / "process-log.jsonl"


def log_event(phase: str, event_type: str, skill: str, files: list = None, extra: dict = None):
    entry = {
        "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill,
        "files_written": files or [],
        "status": "ok",
    }
    if extra:
        entry.update(extra)
    with open(LOG_FILE, "a") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")


def main():
    start_time = time.perf_counter()
    log_event("PLAN", "run_started", "benchmark_suite")

    summary = {}

    # ============================================================
    # 1. QUBO Formulation Analysis
    # ============================================================
    print("\n" + "=" * 60)
    print("Phase 1: QUBO Formulation Analysis")
    print("=" * 60)
    log_event("EXECUTE", "skill_selected", "qubo_formulation")

    N_vrp = 5
    np.random.seed(42)
    coords = np.random.rand(N_vrp, 2) * 100
    dist = np.sqrt(((coords[:, None] - coords[None, :]) ** 2).sum(axis=2))
    vrp_formulator = VRPQUBOFormulator(N_vrp, num_vehicles=2, distance_matrix=dist)
    Q_vrp, meta_vrp = vrp_formulator.build_qubo()

    print(f"  VRP QUBO (N=5): {meta_vrp['num_variables']} variables, {len(Q_vrp)} terms")
    qubo_stats = QUBOFormulator.qubo_stats(Q_vrp)
    print(f"  Density: {qubo_stats['density']:.3f}")
    print(f"  Coeff ratio: {qubo_stats['coeff_ratio']:.2f}")

    fig7_path = plot_qubo_distribution(Q_vrp, title="VRP (N=5) QUBO")
    log_event("EXECUTE", "file_written", "visualization", [fig7_path])
    summary["qubo_meta"] = meta_vrp
    summary["qubo_stats"] = qubo_stats

    # ============================================================
    # 2. Annealing Schedule Comparison
    # ============================================================
    print("\n" + "=" * 60)
    print("Phase 2: Annealing Schedule Comparison")
    print("=" * 60)
    log_event("EXECUTE", "skill_selected", "annealing_solvers")

    Q_bench = generate_random_qubo(20, density=0.5)
    df_schedules = benchmark_schedules(Q_bench, num_reads=50)
    df_schedules.to_csv(RESULTS_DIR / "schedule_comparison.csv", index=False)

    fig3_path = plot_schedule_comparison(df_schedules)
    fig6_path = plot_reverse_annealing_schedule()

    print(df_schedules[["schedule", "best_energy", "mean_energy", "elapsed_sec"]].to_string(index=False))
    log_event("EXECUTE", "file_written", "visualization", [fig3_path, fig6_path])
    summary["best_schedule"] = df_schedules.loc[df_schedules["best_energy"].idxmin(), "schedule"]
    summary["schedule_df"] = df_schedules.to_dict(orient="records")

    # ============================================================
    # 3. Solver Comparison
    # ============================================================
    print("\n" + "=" * 60)
    print("Phase 3: Solver Comparison (n=15)")
    print("=" * 60)
    log_event("EXECUTE", "skill_selected", "classical_solvers")

    Q_15 = generate_random_qubo(15, density=0.5)
    df_solvers = benchmark_solvers(Q_15, num_reads=100)
    df_solvers.to_csv(RESULTS_DIR / "solver_comparison.csv", index=False)

    fig1_path = plot_solver_comparison(df_solvers)

    print(df_solvers[["solver", "best_energy", "mean_energy", "elapsed_sec"]].dropna(
        subset=["best_energy"]).to_string(index=False))
    log_event("EXECUTE", "file_written", "visualization", [fig1_path])
    summary["solver_df"] = df_solvers.dropna(subset=["best_energy"]).to_dict(orient="records")

    # ============================================================
    # 4. Scaling Analysis
    # ============================================================
    print("\n" + "=" * 60)
    print("Phase 4: Problem Scaling Analysis")
    print("=" * 60)
    log_event("EXECUTE", "skill_selected", "benchmark_suite")

    df_scaling = benchmark_scaling(sizes=[5, 8, 10, 15, 20, 30], num_reads=50)
    df_scaling.to_csv(RESULTS_DIR / "scaling_analysis.csv", index=False)

    fig2_path = plot_scaling_analysis(df_scaling)
    log_event("EXECUTE", "file_written", "visualization", [fig2_path])
    summary["scaling_df"] = df_scaling.to_dict(orient="records")

    # ============================================================
    # 5. VRP Case Study
    # ============================================================
    print("\n" + "=" * 60)
    print("Phase 5: VRP Case Study")
    print("=" * 60)
    log_event("EXECUTE", "skill_selected", "vrp_qubo_formulation")

    vrp_results = vrp_case_study(city_counts=[4, 5, 6], num_vehicles=2, num_reads=50)
    with open(RESULTS_DIR / "vrp_results.json", "w") as f:
        json.dump(vrp_results, f, indent=2, ensure_ascii=False)

    fig5_path = plot_vrp_routes(vrp_results, N_key=5)
    log_event("EXECUTE", "file_written", "visualization", [fig5_path])

    for N_key, res in vrp_results.items():
        rows = res.get("solver_results", [])
        valid = [r for r in rows if "best_energy" in r]
        if valid:
            best = min(valid, key=lambda r: r["best_energy"])
            print(f"  VRP(N={N_key}): best={best['best_energy']:.2f} by {best['solver']} in {best['elapsed_sec']:.2f}s")
    summary["vrp_results"] = {
        str(k): {
            "meta": v["meta"],
            "best_solvers": [r for r in v["solver_results"] if "best_energy" in r],
        }
        for k, v in vrp_results.items()
    }

    # ============================================================
    # 6. Minor Embedding Analysis
    # ============================================================
    print("\n" + "=" * 60)
    print("Phase 6: Minor Embedding Analysis")
    print("=" * 60)
    log_event("EXECUTE", "skill_selected", "minor_embedding")

    df_embedding = embedding_analysis(sizes=[10, 20, 30, 50])
    df_embedding.to_csv(RESULTS_DIR / "embedding_analysis.csv", index=False)

    fig4_path = plot_embedding_analysis(df_embedding)
    log_event("EXECUTE", "file_written", "visualization", [fig4_path])
    print(df_embedding[["problem_size", "strategy", "avg_chain_length", "overhead"]].to_string(index=False))
    summary["embedding_df"] = df_embedding.to_dict(orient="records")

    # ============================================================
    # Save summary JSON
    # ============================================================
    with open(RESULTS_DIR / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    elapsed = time.perf_counter() - start_time
    log_event("REPORT", "run_completed", "evaluation_framework", extra={"elapsed_sec": elapsed})

    print(f"\n✓ All phases complete in {elapsed:.1f}s")
    print(f"  Results: {RESULTS_DIR}/")
    print(f"  Figures: {FIGURES_DIR}/")
    print(f"  Logs:    {LOG_FILE}")

    return summary


if __name__ == "__main__":
    summary = main()
