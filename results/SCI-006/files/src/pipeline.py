"""
Main Pipeline: AlphaFold2-based Protein-Ligand Binding Affinity Prediction System

Orchestrates all six modules:
1. pLDDT-based docking suitability assessment
2. MD binding pose refinement
3. FEP vs metadynamics comparison
4. GNN binding affinity prediction
5. Activity cliff detection
6. Multi-objective Pareto optimization
"""

import sys
import os
import json
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plddt_assessment import run_plddt_analysis
from md_refinement import run_md_refinement
from fep_metadynamics import run_fep_metadynamics_comparison
from gnn_affinity import run_gnn_training
from activity_cliff import run_activity_cliff_analysis
from pareto_optimization import run_pareto_optimization


def main():
    output_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")
    os.makedirs(output_dir, exist_ok=True)

    results = {}
    start_time = time.time()

    print("\n" + "=" * 70)
    print("  AlphaFold2-based Protein-Ligand Binding Affinity Prediction Pipeline")
    print("=" * 70 + "\n")

    # Module 1
    print("\n[1/6] Running pLDDT Assessment...")
    results['plddt'] = run_plddt_analysis(output_dir)

    # Module 2
    print("\n[2/6] Running MD Refinement...")
    results['md'] = run_md_refinement(output_dir)

    # Module 3
    print("\n[3/6] Running FEP vs Metadynamics Comparison...")
    results['fep_meta'] = run_fep_metadynamics_comparison(output_dir)

    # Module 4
    print("\n[4/6] Running GNN Training...")
    results['gnn'] = run_gnn_training(output_dir)

    # Module 5
    print("\n[5/6] Running Activity Cliff Analysis...")
    results['activity_cliff'] = run_activity_cliff_analysis(output_dir)

    # Module 6
    print("\n[6/6] Running Pareto Optimization...")
    results['pareto'] = run_pareto_optimization(output_dir)

    elapsed = time.time() - start_time
    print(f"\n{'=' * 70}")
    print(f"Pipeline completed in {elapsed:.1f} seconds")
    print(f"{'=' * 70}")

    # Save results summary
    summary = {
        'plddt': {k: v for k, v in list(results['plddt'].items())[:3]},
        'md': {k: {kk: vv for kk, vv in v.items() if not isinstance(vv, list)}
               for k, v in results['md'].items()},
        'fep_meta': results['fep_meta'],
        'gnn': {k: v for k, v in results['gnn'].items() if k not in ('train_losses', 'val_losses')},
        'activity_cliff': results['activity_cliff'],
        'pareto': {k: v for k, v in results['pareto'].items() if k != 'history'},
    }

    with open(os.path.join(os.path.dirname(output_dir), 'results_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    print("\nResults saved to results_summary.json")

    return results


if __name__ == '__main__':
    main()
