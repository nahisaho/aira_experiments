"""
Main runner: Execute all PINN benchmarks and generate results.
"""

import json
import sys
import os
import time
import numpy as np
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def save_results(results, filename):
    """Save results as JSON (convert numpy arrays to lists)."""
    def convert(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [convert(v) for v in obj]
        return obj

    with open(filename, "w") as f:
        json.dump(convert(results), f, indent=2)


def main():
    all_results = {}
    start_time = time.time()

    # Module 1: Fourier Features
    print("=" * 60)
    print("MODULE 1: Multi-scale Fourier Feature Embedding")
    print("=" * 60)
    try:
        from src.fourier_features import run_multiscale_benchmark
        results = run_multiscale_benchmark()
        all_results["fourier_features"] = results
        save_results(results, "results/fourier_features.json")
        print("[OK] Fourier features benchmark complete")
    except Exception as e:
        print(f"[ERROR] {e}")
        traceback.print_exc()
        all_results["fourier_features"] = {"error": str(e)}

    # Module 2: Inverse Problem
    print("\n" + "=" * 60)
    print("MODULE 2: Inverse Problem & Uncertainty Quantification")
    print("=" * 60)
    try:
        from src.inverse_problem import run_inverse_problem_benchmark
        results = run_inverse_problem_benchmark()
        all_results["inverse_problem"] = results
        save_results(results, "results/inverse_problem.json")
        print("[OK] Inverse problem benchmark complete")
    except Exception as e:
        print(f"[ERROR] {e}")
        traceback.print_exc()
        all_results["inverse_problem"] = {"error": str(e)}

    # Module 3: Causal Training
    print("\n" + "=" * 60)
    print("MODULE 3: Causal Training")
    print("=" * 60)
    try:
        from src.causal_training import run_causal_training_benchmark
        results = run_causal_training_benchmark()
        all_results["causal_training"] = results
        save_results(results, "results/causal_training.json")
        print("[OK] Causal training benchmark complete")
    except Exception as e:
        print(f"[ERROR] {e}")
        traceback.print_exc()
        all_results["causal_training"] = {"error": str(e)}

    # Module 4: Adaptive Collocation
    print("\n" + "=" * 60)
    print("MODULE 4: Adaptive Collocation Points")
    print("=" * 60)
    try:
        from src.adaptive_collocation import run_adaptive_collocation_benchmark
        results = run_adaptive_collocation_benchmark()
        all_results["adaptive_collocation"] = results
        save_results(results, "results/adaptive_collocation.json")
        print("[OK] Adaptive collocation benchmark complete")
    except Exception as e:
        print(f"[ERROR] {e}")
        traceback.print_exc()
        all_results["adaptive_collocation"] = {"error": str(e)}

    # Module 5: Operator Learning
    print("\n" + "=" * 60)
    print("MODULE 5: DeepONet vs FNO Comparison")
    print("=" * 60)
    try:
        from src.operator_learning import run_operator_learning_benchmark
        results = run_operator_learning_benchmark()
        all_results["operator_learning"] = results
        save_results(results, "results/operator_learning.json")
        print("[OK] Operator learning benchmark complete")
    except Exception as e:
        print(f"[ERROR] {e}")
        traceback.print_exc()
        all_results["operator_learning"] = {"error": str(e)}

    # Module 6: Navier-Stokes
    print("\n" + "=" * 60)
    print("MODULE 6: Navier-Stokes Case Study")
    print("=" * 60)
    try:
        from src.navier_stokes import run_navier_stokes_benchmark
        results = run_navier_stokes_benchmark()
        all_results["navier_stokes"] = results
        save_results(results, "results/navier_stokes.json")
        print("[OK] Navier-Stokes benchmark complete")
    except Exception as e:
        print(f"[ERROR] {e}")
        traceback.print_exc()
        all_results["navier_stokes"] = {"error": str(e)}

    elapsed = time.time() - start_time
    print(f"\n{'=' * 60}")
    print(f"All benchmarks complete in {elapsed:.1f} seconds")
    print(f"{'=' * 60}")

    save_results(all_results, "results/all_results.json")
    return all_results


if __name__ == "__main__":
    main()
