#!/usr/bin/env python3
"""
Main runner: Execute all experiments and generate figures/results/report.
"""

import json
import os
import sys
import time
import traceback
import numpy as np

# Ensure workspace is on path
WORKSPACE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, WORKSPACE)

RESULTS_DIR = os.path.join(WORKSPACE, "results")
LOGS_DIR = os.path.join(WORKSPACE, "logs")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(LOGS_DIR, exist_ok=True)

LOG_FILE = os.path.join(LOGS_DIR, "process-log.jsonl")


def log_event(phase, event_type, skill="adaptive-experiments", **kwargs):
    entry = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill,
        **kwargs,
    }
    with open(LOG_FILE, "a") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")


def run_all():
    log_event("init", "run_started")
    all_results = {}

    # ========== 1. Kernel Comparison ==========
    print("=" * 60)
    print("[1/6] Kernel Selection & Hyperparameter Optimization")
    print("=" * 60)
    try:
        log_event("kernel", "experiment_started")
        from bayesopt_framework.kernel_selection import run_kernel_comparison
        kernel_res = run_kernel_comparison(n_samples=80, dim=6, seed=42)
        all_results["kernel_comparison"] = kernel_res
        print(f"  Best kernel: {kernel_res['best_kernel']}")
        for r in kernel_res["comparison"]:
            print(f"    {r['kernel']:15s} NLPD={r['mean_nlpd']:.4f} ± {r['std_nlpd']:.4f}  ({r['fit_time_s']:.1f}s)")

        # Visualization
        from bayesopt_framework.visualization import plot_kernel_comparison, plot_lengthscales
        plot_kernel_comparison(kernel_res["comparison"])
        if kernel_res["lengthscales"]:
            plot_lengthscales(kernel_res["lengthscales"])
        log_event("kernel", "experiment_completed", files_written=["figures/kernel_comparison.png", "figures/lengthscales.png"])
    except Exception as e:
        print(f"  ERROR: {e}")
        traceback.print_exc()
        log_event("kernel", "experiment_failed", error=str(e))

    # ========== 2. Acquisition Function Comparison ==========
    print("\n" + "=" * 60)
    print("[2/6] Acquisition Function Comparison")
    print("=" * 60)
    try:
        log_event("acquisition", "experiment_started")
        from bayesopt_framework.acquisition_functions import run_acquisition_comparison
        acq_res = run_acquisition_comparison(seed=42)
        all_results["acquisition_comparison"] = {
            k: {kk: vv for kk, vv in v.items() if kk != "trials"}
            for k, v in acq_res["comparison"].items()
        }
        all_results["acquisition_selection"] = acq_res["selection_examples"]

        for name, data in acq_res["comparison"].items():
            print(f"  {name:5s}: final_best={data['mean_best']:.4f} ± {data['std_best']:.4f}")

        from bayesopt_framework.visualization import plot_acquisition_comparison
        plot_acquisition_comparison(acq_res["comparison"])
        log_event("acquisition", "experiment_completed", files_written=["figures/acquisition_comparison.png"])
    except Exception as e:
        print(f"  ERROR: {e}")
        traceback.print_exc()
        log_event("acquisition", "experiment_failed", error=str(e))

    # ========== 3. Batch Optimization ==========
    print("\n" + "=" * 60)
    print("[3/6] Batch Optimization (Parallel Proposals)")
    print("=" * 60)
    try:
        log_event("batch", "experiment_started")
        from bayesopt_framework.batch_optimization import compare_batch_methods
        batch_res = compare_batch_methods(seed=42)
        all_results["batch_comparison"] = {
            m: {k: v for k, v in d.items() if k not in ("best_values", "n_evals")}
            for m, d in batch_res.items()
        }

        for method, data in batch_res.items():
            print(f"  {method:20s}: best={data['final_best']:.4f}, time={data['total_time']:.1f}s")

        from bayesopt_framework.visualization import plot_batch_comparison
        plot_batch_comparison(batch_res)
        log_event("batch", "experiment_completed", files_written=["figures/batch_comparison.png"])
    except Exception as e:
        print(f"  ERROR: {e}")
        traceback.print_exc()
        log_event("batch", "experiment_failed", error=str(e))

    # ========== 4. Multi-Objective BO ==========
    print("\n" + "=" * 60)
    print("[4/6] Multi-Objective Bayesian Optimization (EHVI)")
    print("=" * 60)
    try:
        log_event("mobo", "experiment_started")
        from bayesopt_framework.multi_objective import run_mobo_experiment
        mobo_res = run_mobo_experiment(seed=42)
        all_results["multi_objective"] = mobo_res
        print(f"  Final HV: {mobo_res['final_hv']:.4f}")
        print(f"  Pareto points: {mobo_res['n_pareto']}")

        from bayesopt_framework.visualization import plot_pareto_front
        plot_pareto_front(mobo_res["pareto_Y"], mobo_res["hv_history"])
        log_event("mobo", "experiment_completed", files_written=["figures/pareto_front.png"])
    except Exception as e:
        print(f"  ERROR: {e}")
        traceback.print_exc()
        log_event("mobo", "experiment_failed", error=str(e))

    # ========== 5. High-Dimensional BO ==========
    print("\n" + "=" * 60)
    print("[5/6] High-Dimensional BO (REMBO, HeSBO)")
    print("=" * 60)
    try:
        log_event("highdim", "experiment_started")
        from bayesopt_framework.high_dimensional import run_highdim_comparison
        hd_res = run_highdim_comparison(dim_high=25, dim_low=6, n_init=15, n_iter=40, seed=42)
        all_results["high_dimensional"] = {
            m: {"final_best": d["final_best"]}
            for m, d in hd_res.items()
        }
        for method, data in hd_res.items():
            print(f"  {method:15s}: final_best={data['final_best']:.4f}")

        from bayesopt_framework.visualization import plot_highdim_comparison
        plot_highdim_comparison(hd_res)
        log_event("highdim", "experiment_completed", files_written=["figures/highdim_comparison.png"])
    except Exception as e:
        print(f"  ERROR: {e}")
        traceback.print_exc()
        log_event("highdim", "experiment_failed", error=str(e))

    # ========== 6. Chemical Reaction Case Study ==========
    print("\n" + "=" * 60)
    print("[6/6] Chemical Reaction Optimization Case Study")
    print("=" * 60)
    try:
        log_event("chemical", "experiment_started")
        from bayesopt_framework.chemical_optimization import run_full_case_study
        chem_res = run_full_case_study(seed=42)
        all_results["chemical_optimization"] = chem_res

        so = chem_res["single_objective"]
        print(f"  Single-obj best yield: {so['best_yield']}%")
        print(f"  Best conditions:")
        for name, val in so["best_conditions"].items():
            print(f"    {name}: {val['value']} {val['unit']}")

        mo = chem_res["multi_objective"]
        print(f"\n  Multi-obj HV: {mo['final_hv']:.2f}")
        print(f"  Pareto solutions: {mo['n_pareto']}")

        from bayesopt_framework.visualization import plot_chemical_convergence, plot_chemical_pareto
        plot_chemical_convergence(so["convergence"])
        if mo["pareto_front"]:
            plot_chemical_pareto(mo["pareto_front"])
        log_event("chemical", "experiment_completed",
                  files_written=["figures/chemical_convergence.png", "figures/chemical_pareto.png"])
    except Exception as e:
        print(f"  ERROR: {e}")
        traceback.print_exc()
        log_event("chemical", "experiment_failed", error=str(e))

    # ========== Save All Results ==========
    results_path = os.path.join(RESULTS_DIR, "all_results.json")

    def make_serializable(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, dict):
            return {k: make_serializable(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [make_serializable(i) for i in obj]
        return obj

    with open(results_path, "w") as f:
        json.dump(make_serializable(all_results), f, indent=2, ensure_ascii=False)

    log_event("finalize", "run_completed", files_written=[results_path])
    print("\n" + "=" * 60)
    print("All experiments completed. Results saved to results/all_results.json")
    print("=" * 60)

    return all_results


if __name__ == "__main__":
    run_all()
