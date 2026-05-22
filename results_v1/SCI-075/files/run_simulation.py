"""
Main entry point for the semi-autonomous suturing simulation.
Runs all subsystems and generates results.
"""

import sys
import os
import json
import time
import numpy as np

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.simulation.dvrk_sim import SuturingSimulation, SimulationConfig
from src.lfd.gmm_gmr import SuturingLfDPipeline


def run_simulation():
    """Execute the complete suturing simulation pipeline."""
    print("=" * 70)
    print("Semi-Autonomous Suturing System — dVRK Simulation")
    print("=" * 70)

    config = SimulationConfig(
        dt=0.001,
        sim_duration=60.0,
        tissue_model_type="mass_spring",
        lfd_method="gmr",
        visual_servo_mode="pbvs",
        enable_safety=True,
    )

    sim = SuturingSimulation(config)

    # Phase 1: Generate synthetic demonstrations
    print("\n[1/4] Generating synthetic expert demonstrations...")
    sim.generate_synthetic_demonstrations(n_demos=5)
    print(f"  Generated {5} demonstrations for each of {len(SuturingLfDPipeline.PHASES)} phases")

    # Phase 2: Learn from demonstrations
    print("\n[2/4] Learning from demonstrations (GMM/GMR)...")
    stats = sim.learn_from_demonstrations()
    for phase, s in stats.items():
        print(f"  {phase}: {s['n_demonstrations']} demos, "
              f"avg duration={s['avg_duration']:.2f}s, "
              f"avg max force={s['avg_max_force']:.2f}N")

    # Phase 3: Run full suturing simulation
    print("\n[3/4] Running full suturing simulation...")
    results = sim.run_full_suturing()

    # Print results
    print("\n" + "-" * 50)
    print("Phase Results:")
    print("-" * 50)
    for phase, metrics in results['phase_metrics'].items():
        if 'error' in metrics:
            print(f"  {phase}: {metrics['error']}")
            continue
        print(f"  {phase}:")
        print(f"    Duration: {metrics['duration']*1000:.1f} ms")
        print(f"    Max Force: {metrics['max_force']:.3f} N")
        print(f"    Mean Tracking Error: {metrics['mean_tracking_error']:.3f} mm")
        print(f"    Max Tracking Error: {metrics['max_tracking_error']:.3f} mm")
        print(f"    Max Tissue Strain: {metrics['max_strain']:.4f}")
        print(f"    Safety Violations: {metrics['safety_violations']}")

    overall = results['overall']
    print(f"\n{'=' * 50}")
    print(f"Overall Results:")
    print(f"  Total Time: {overall['total_time']:.3f} s")
    print(f"  Success: {overall['success']}")
    print(f"  Total Safety Violations: {overall['safety_violations']}")
    print(f"  Max Tissue Strain: {overall['tissue_max_strain']:.4f}")

    # Phase 4: Save results
    print("\n[4/4] Saving results...")
    workspace_dir = os.path.dirname(os.path.abspath(__file__))
    metrics_path, log_path = sim.save_results(workspace_dir)
    print(f"  Metrics: {metrics_path}")
    print(f"  Logs: {log_path}")

    # Save detailed results for report
    detailed_results = {
        'config': {
            'dt': config.dt,
            'tissue_model': config.tissue_model_type,
            'lfd_method': config.lfd_method,
            'visual_servo': config.visual_servo_mode,
            'safety_enabled': config.enable_safety,
        },
        'lfd_stats': {k: {sk: float(sv) for sk, sv in v.items()} for k, v in stats.items()},
        'phase_results': {},
        'overall': {
            'total_time': overall['total_time'],
            'success': overall['success'],
            'safety_violations': overall['safety_violations'],
            'tissue_max_strain': overall['tissue_max_strain'],
        }
    }

    for phase, metrics in results['phase_metrics'].items():
        detailed_results['phase_results'][phase] = {
            k: float(v) if isinstance(v, (int, float, np.floating, np.integer)) else v
            for k, v in metrics.items()
        }

    results_path = os.path.join(workspace_dir, 'results', 'detailed_results.json')
    with open(results_path, 'w') as f:
        json.dump(detailed_results, f, indent=2, default=str)
    print(f"  Detailed results: {results_path}")

    print("\n✓ Simulation complete.")
    return detailed_results


if __name__ == '__main__':
    run_simulation()
