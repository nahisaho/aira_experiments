"""
Main Pipeline: AlphaFold2-Enhanced Protein-Ligand Binding Affinity Prediction

Orchestrates all modules into a unified computational pipeline.
"""

import sys
import os
import json
import numpy as np
from datetime import datetime

# Add src to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.plddt_assessment import (
    generate_synthetic_plddt_profile, assess_binding_site,
    generate_docking_recommendation, compute_local_confidence_map
)
from src.md_refinement import (
    MDRefinementPipeline, MDParameters, MDProtocol
)
from src.free_energy import (
    simulate_fep_results, simulate_metadynamics_results, compare_methods
)
from src.gnn_predictor import (
    simulate_gnn_training, simulate_gnn_predictions, evaluate_predictions
)
from src.activity_cliff import (
    generate_synthetic_compounds, detect_activity_cliffs,
    analyze_chemical_space, generate_exploration_strategy
)
from src.multi_objective import (
    generate_lead_optimization_candidates, nsga2_optimize
)


def run_full_pipeline():
    """Execute the complete binding affinity prediction pipeline."""
    results = {}
    timestamp = datetime.now().isoformat()
    
    print("=" * 70)
    print("AlphaFold2-Enhanced Binding Affinity Prediction Pipeline")
    print(f"Started: {timestamp}")
    print("=" * 70)
    
    # === Phase 1: pLDDT Assessment ===
    print("\n[Phase 1] pLDDT-based Docking Suitability Assessment...")
    residues = generate_synthetic_plddt_profile(n_residues=300, seed=42)
    binding_site_ids = list(range(130, 170))  # Residues 130-169
    assessment = assess_binding_site(residues, binding_site_ids)
    recommendation = generate_docking_recommendation(assessment)
    local_confidence = compute_local_confidence_map(residues)
    
    results["plddt_assessment"] = {
        "n_residues": len(residues),
        "binding_site_residues": len(binding_site_ids),
        "recommendation": recommendation,
        "plddt_profile": [(r.residue_id, r.plddt) for r in residues],
        "local_confidence": local_confidence,
    }
    print(f"  Overall suitability: {recommendation['overall_suitability']}")
    print(f"  Recommended strategy: {recommendation['recommended_strategy']}")
    print(f"  Mean binding site pLDDT: {recommendation['binding_site_metrics']['mean_plddt']:.1f}")
    
    # === Phase 2: MD Refinement ===
    print("\n[Phase 2] Molecular Dynamics Refinement...")
    md_pipeline = MDRefinementPipeline(MDParameters.from_protocol(MDProtocol.STANDARD))
    
    plddt_dict = {r.residue_id: r.plddt for r in residues}
    system_info = md_pipeline.prepare_system("protein.pdb", "ligand.sdf", plddt_dict)
    metrics = md_pipeline.simulate_trajectory_metrics(seed=42)
    traj_data = md_pipeline.get_trajectory_data()
    poses = md_pipeline.cluster_poses(n_clusters=5, seed=42)
    
    results["md_refinement"] = {
        "system_info": system_info,
        "metrics": {
            "simulation_time_ns": metrics.simulation_time_ns,
            "protein_rmsd": f"{metrics.protein_rmsd_mean:.3f} ± {metrics.protein_rmsd_std:.3f} nm",
            "ligand_rmsd": f"{metrics.ligand_rmsd_mean:.3f} ± {metrics.ligand_rmsd_std:.3f} nm",
            "mm_pbsa": f"{metrics.mm_pbsa_mean:.1f} ± {metrics.mm_pbsa_std:.1f} kcal/mol",
            "converged": metrics.is_converged,
            "convergence_time_ns": metrics.convergence_time_ns,
        },
        "n_clusters": len(poses),
        "trajectory_data": traj_data,
        "poses": [
            {
                "cluster_id": p.cluster_id,
                "population": f"{p.population_fraction:.1%}",
                "binding_energy": f"{p.binding_energy:.1f} kcal/mol",
            }
            for p in poses
        ],
    }
    print(f"  Simulation time: {metrics.simulation_time_ns:.0f} ns")
    print(f"  Protein RMSD: {metrics.protein_rmsd_mean:.3f} ± {metrics.protein_rmsd_std:.3f} nm")
    print(f"  Found {len(poses)} pose clusters")
    
    # === Phase 3: Free Energy Calculations ===
    print("\n[Phase 3] Free Energy Calculations (FEP vs Metadynamics)...")
    fep_results = simulate_fep_results(n_perturbations=10, seed=42)
    metad_results = simulate_metadynamics_results(n_ligands=10, seed=42)
    comparison = compare_methods(fep_results, metad_results)
    
    results["free_energy"] = {
        "fep": {
            "n_perturbations": len(fep_results),
            "rmse": f"{comparison.fep_rmse:.2f} kcal/mol",
            "mae": f"{comparison.fep_mae:.2f} kcal/mol",
            "r_squared": f"{comparison.fep_r_squared:.3f}",
            "kendall_tau": f"{comparison.fep_kendall_tau:.3f}",
            "gpu_hours": f"{comparison.fep_total_gpu_hours:.0f}",
        },
        "metadynamics": {
            "n_ligands": len(metad_results),
            "rmse": f"{comparison.metad_rmse:.2f} kcal/mol",
            "mae": f"{comparison.metad_mae:.2f} kcal/mol",
            "r_squared": f"{comparison.metad_r_squared:.3f}",
            "kendall_tau": f"{comparison.metad_kendall_tau:.3f}",
            "gpu_hours": f"{comparison.metad_total_gpu_hours:.0f}",
        },
        "recommendation": comparison.recommended_method,
        "reasoning": comparison.reasoning,
        "fep_data": [(r.ddg_fep, r.ddg_experimental) for r in fep_results],
        "metad_data": [(r.dg_binding, r.dg_experimental) for r in metad_results],
    }
    print(f"  FEP RMSE: {comparison.fep_rmse:.2f} kcal/mol | R²: {comparison.fep_r_squared:.3f}")
    print(f"  Metadynamics RMSE: {comparison.metad_rmse:.2f} kcal/mol | R²: {comparison.metad_r_squared:.3f}")
    print(f"  Recommended: {comparison.recommended_method}")
    
    # === Phase 4: GNN Predictions ===
    print("\n[Phase 4] GNN Binding Affinity Prediction...")
    training_curves = simulate_gnn_training(n_epochs=200, seed=42)
    predictions = simulate_gnn_predictions(n_compounds=200, seed=42)
    evaluation = evaluate_predictions(predictions)
    
    results["gnn"] = {
        "training": {
            "epochs": len(training_curves["epoch"]),
            "final_train_loss": f"{training_curves['train_loss'][-1]:.4f}",
            "final_val_loss": f"{training_curves['val_loss'][-1]:.4f}",
            "best_val_rmse": f"{min(training_curves['val_rmse']):.4f}",
        },
        "evaluation": {
            "rmse": f"{evaluation.rmse:.3f}",
            "mae": f"{evaluation.mae:.3f}",
            "r_squared": f"{evaluation.r_squared:.3f}",
            "pearson_r": f"{evaluation.pearson_r:.3f}",
            "spearman_rho": f"{evaluation.spearman_rho:.3f}",
        },
        "training_curves": training_curves,
        "predictions": [(p.predicted_pki, p.experimental_pki, p.uncertainty) for p in predictions],
    }
    print(f"  Test RMSE: {evaluation.rmse:.3f} pKi units")
    print(f"  R²: {evaluation.r_squared:.3f} | Pearson r: {evaluation.pearson_r:.3f}")
    
    # === Phase 5: Activity Cliff Detection ===
    print("\n[Phase 5] Activity Cliff Detection...")
    compounds = generate_synthetic_compounds(n_compounds=100, seed=42)
    cliffs = detect_activity_cliffs(compounds, similarity_threshold=0.65,
                                     activity_threshold=1.0, sali_threshold=3.0)
    space_analysis = analyze_chemical_space(compounds)
    exploration = generate_exploration_strategy(cliffs, space_analysis, compounds)
    
    results["activity_cliffs"] = {
        "n_compounds": len(compounds),
        "n_cliffs_detected": len(cliffs),
        "top_cliffs": [
            {
                "pair": (c.compound_a, c.compound_b),
                "similarity": f"{c.similarity:.3f}",
                "activity_diff": f"{c.activity_diff:.2f} pKi",
                "sali": f"{c.sali:.1f}",
            }
            for c in cliffs[:10]
        ],
        "chemical_space": {
            "n_clusters": space_analysis.n_clusters,
            "coverage": f"{space_analysis.coverage_score:.3f}",
            "diversity": f"{space_analysis.diversity_score:.3f}",
        },
        "exploration_strategy": exploration,
        "pca_coords": space_analysis.pca_coordinates.tolist() if space_analysis.pca_coordinates is not None else [],
        "compounds_data": [(c.compound_id, c.pki, c.scaffold) for c in compounds],
    }
    print(f"  Detected {len(cliffs)} activity cliffs")
    print(f"  Chemical space: {space_analysis.n_clusters} clusters, "
          f"diversity={space_analysis.diversity_score:.3f}")
    
    # === Phase 6: Multi-Objective Optimization ===
    print("\n[Phase 6] Multi-Objective Lead Optimization...")
    candidates, objectives = generate_lead_optimization_candidates(n_candidates=200, seed=42)
    pareto_result = nsga2_optimize(candidates, objectives, n_generations=50,
                                    population_size=100, seed=42)
    
    results["optimization"] = {
        "n_candidates": len(candidates),
        "n_objectives": len(objectives),
        "objectives": [o.name for o in objectives],
        "n_pareto_optimal": pareto_result.n_pareto_optimal,
        "hypervolume": f"{pareto_result.hypervolume:.2f}",
        "generation_history": pareto_result.generation_history,
        "pareto_solutions": [
            {
                "id": s.compound_id,
                "values": {k: round(v, 3) for k, v in s.objective_values.items()},
                "crowding_distance": round(s.crowding_distance, 3),
            }
            for s in pareto_result.solutions[:20]
        ],
    }
    print(f"  Pareto-optimal solutions: {pareto_result.n_pareto_optimal}")
    print(f"  Hypervolume: {pareto_result.hypervolume:.2f}")
    
    # Save results
    print("\n" + "=" * 70)
    print("Saving results...")
    
    # Save full results as JSON
    results_serializable = _make_serializable(results)
    with open("results/pipeline_results.json", "w") as f:
        json.dump(results_serializable, f, indent=2, default=str)
    
    print(f"Results saved to results/pipeline_results.json")
    print(f"Pipeline completed: {datetime.now().isoformat()}")
    print("=" * 70)
    
    return results


def _make_serializable(obj):
    """Convert numpy types to Python native types for JSON serialization."""
    if isinstance(obj, dict):
        return {k: _make_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [_make_serializable(v) for v in obj]
    elif isinstance(obj, tuple):
        return [_make_serializable(v) for v in obj]
    elif isinstance(obj, (np.integer,)):
        return int(obj)
    elif isinstance(obj, (np.floating,)):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.bool_):
        return bool(obj)
    else:
        return obj


if __name__ == "__main__":
    results = run_full_pipeline()
