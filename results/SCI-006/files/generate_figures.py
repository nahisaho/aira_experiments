"""
Figure Generation Script

Generates all publication-quality figures for the binding affinity prediction pipeline.
"""

import sys
import os
import json
import numpy as np

# Ensure matplotlib works headless
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap

# Style
plt.rcParams.update({
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.family': 'DejaVu Sans',
})

# Colorblind-friendly palette (viridis-inspired)
COLORS = ['#440154', '#31688e', '#35b779', '#fde725', '#e76f51', '#264653']
CMAP = plt.cm.viridis


def load_or_generate_data():
    """Load pipeline results or generate fresh."""
    try:
        with open("results/pipeline_results.json", "r") as f:
            return json.load(f)
    except FileNotFoundError:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from pipeline import run_full_pipeline
        return run_full_pipeline()


def fig1_plddt_profile(results, output_dir="figures"):
    """Figure 1: pLDDT profile and binding site assessment."""
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), gridspec_kw={'height_ratios': [2, 1]})
    
    plddt_data = results["plddt_assessment"]["plddt_profile"]
    residue_ids = [d[0] for d in plddt_data]
    plddts = [d[1] for d in plddt_data]
    
    # Top: Full pLDDT profile
    ax = axes[0]
    colors_by_score = []
    for p in plddts:
        if p >= 90:
            colors_by_score.append('#35b779')
        elif p >= 70:
            colors_by_score.append('#31688e')
        elif p >= 50:
            colors_by_score.append('#fde725')
        else:
            colors_by_score.append('#e76f51')
    
    ax.bar(residue_ids, plddts, color=colors_by_score, width=1.0, alpha=0.8)
    
    # Binding site region
    ax.axvspan(130, 170, alpha=0.15, color='red', label='Binding site')
    
    # Threshold lines
    for threshold, label, ls in [(90, 'Very high (≥90)', '--'),
                                   (70, 'Confident (≥70)', '-.'),
                                   (50, 'Low (<50)', ':')]:
        ax.axhline(y=threshold, color='gray', linestyle=ls, alpha=0.5, linewidth=0.8)
    
    ax.set_xlabel('Residue Number')
    ax.set_ylabel('pLDDT Score')
    ax.set_title('AlphaFold2 Per-Residue Confidence (pLDDT) Profile')
    ax.set_ylim(0, 105)
    ax.legend(loc='lower right')
    
    # Bottom: Binding site detail
    ax2 = axes[1]
    bs_data = [(d[0], d[1]) for d in plddt_data if 120 <= d[0] <= 180]
    bs_ids = [d[0] for d in bs_data]
    bs_plddts = [d[1] for d in bs_data]
    
    bs_colors = []
    for p in bs_plddts:
        if p >= 90:
            bs_colors.append('#35b779')
        elif p >= 70:
            bs_colors.append('#31688e')
        elif p >= 50:
            bs_colors.append('#fde725')
        else:
            bs_colors.append('#e76f51')
    
    ax2.bar(bs_ids, bs_plddts, color=bs_colors, width=1.0, alpha=0.9)
    ax2.axvspan(130, 170, alpha=0.15, color='red')
    ax2.axhline(y=70, color='gray', linestyle='--', alpha=0.5)
    ax2.axhline(y=90, color='gray', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Residue Number')
    ax2.set_ylabel('pLDDT Score')
    ax2.set_title('Binding Site Region (Residues 120-180)')
    ax2.set_ylim(0, 105)
    
    # Add assessment text
    metrics = results["plddt_assessment"]["recommendation"]["binding_site_metrics"]
    info_text = (f"Mean pLDDT: {metrics['mean_plddt']:.1f}\n"
                 f"Min pLDDT: {metrics['min_plddt']:.1f}\n"
                 f"Strategy: {results['plddt_assessment']['recommendation']['recommended_strategy']}")
    ax2.text(0.98, 0.95, info_text, transform=ax2.transAxes,
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
             fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/fig1_plddt_profile.png")
    plt.close()
    print("  Saved fig1_plddt_profile.png")


def fig2_md_trajectory(results, output_dir="figures"):
    """Figure 2: MD trajectory analysis."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    traj = results["md_refinement"]["trajectory_data"]
    time_ns = np.array(traj["time_ns"])
    prot_rmsd = np.array(traj["protein_rmsd"])
    lig_rmsd = np.array(traj["ligand_rmsd"])
    rmsf = np.array(traj["rmsf"])
    
    # Top-left: Protein RMSD
    ax = axes[0, 0]
    ax.plot(time_ns, prot_rmsd, color=COLORS[1], linewidth=0.5, alpha=0.6)
    # Running average
    window = max(1, len(prot_rmsd) // 100)
    prot_smooth = np.convolve(prot_rmsd, np.ones(window)/window, mode='valid')
    ax.plot(time_ns[:len(prot_smooth)], prot_smooth, color=COLORS[0], linewidth=2)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (nm)')
    ax.set_title('Protein Backbone RMSD')
    ax.axhline(y=np.mean(prot_rmsd), color='gray', linestyle='--', alpha=0.5)
    
    # Top-right: Ligand RMSD
    ax = axes[0, 1]
    ax.plot(time_ns, lig_rmsd, color=COLORS[2], linewidth=0.5, alpha=0.6)
    lig_smooth = np.convolve(lig_rmsd, np.ones(window)/window, mode='valid')
    ax.plot(time_ns[:len(lig_smooth)], lig_smooth, color=COLORS[4], linewidth=2)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (nm)')
    ax.set_title('Ligand RMSD')
    
    # Bottom-left: RMSF
    ax = axes[1, 0]
    residue_ids = np.arange(1, len(rmsf) + 1)
    ax.fill_between(residue_ids, 0, rmsf, alpha=0.4, color=COLORS[1])
    ax.plot(residue_ids, rmsf, color=COLORS[0], linewidth=1)
    ax.axvspan(130, 170, alpha=0.15, color='red', label='Binding site')
    ax.set_xlabel('Residue Number')
    ax.set_ylabel('RMSF (nm)')
    ax.set_title('Per-Residue RMSF')
    ax.legend()
    
    # Bottom-right: RMSD distribution
    ax = axes[1, 1]
    ax.hist(prot_rmsd, bins=50, alpha=0.6, color=COLORS[1], label='Protein', density=True)
    ax.hist(lig_rmsd, bins=50, alpha=0.6, color=COLORS[2], label='Ligand', density=True)
    ax.set_xlabel('RMSD (nm)')
    ax.set_ylabel('Density')
    ax.set_title('RMSD Distribution')
    ax.legend()
    
    plt.suptitle('Molecular Dynamics Trajectory Analysis (100 ns)', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/fig2_md_trajectory.png")
    plt.close()
    print("  Saved fig2_md_trajectory.png")


def fig3_free_energy_comparison(results, output_dir="figures"):
    """Figure 3: FEP vs Metadynamics comparison."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # Left: FEP correlation
    ax = axes[0]
    fep_data = results["free_energy"]["fep_data"]
    fep_calc = [d[0] for d in fep_data]
    fep_exp = [d[1] for d in fep_data]
    
    ax.scatter(fep_exp, fep_calc, c=COLORS[1], s=80, alpha=0.7, edgecolors='white', zorder=3)
    lims = [min(min(fep_exp), min(fep_calc)) - 0.5,
            max(max(fep_exp), max(fep_calc)) + 0.5]
    ax.plot(lims, lims, 'k--', alpha=0.3, label='y = x')
    ax.fill_between(lims, [l - 1 for l in lims], [l + 1 for l in lims],
                     alpha=0.1, color='gray', label='±1 kcal/mol')
    ax.set_xlabel('Experimental ΔΔG (kcal/mol)')
    ax.set_ylabel('FEP ΔΔG (kcal/mol)')
    ax.set_title(f'FEP (RMSE={results["free_energy"]["fep"]["rmse"]})')
    ax.legend(loc='lower right', fontsize=9)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    
    # Middle: Metadynamics correlation
    ax = axes[1]
    metad_data = results["free_energy"]["metad_data"]
    metad_calc = [d[0] for d in metad_data]
    metad_exp = [d[1] for d in metad_data]
    
    ax.scatter(metad_exp, metad_calc, c=COLORS[2], s=80, alpha=0.7, edgecolors='white', zorder=3)
    lims2 = [min(min(metad_exp), min(metad_calc)) - 0.5,
             max(max(metad_exp), max(metad_calc)) + 0.5]
    ax.plot(lims2, lims2, 'k--', alpha=0.3, label='y = x')
    ax.fill_between(lims2, [l - 1 for l in lims2], [l + 1 for l in lims2],
                     alpha=0.1, color='gray', label='±1 kcal/mol')
    ax.set_xlabel('Experimental ΔG (kcal/mol)')
    ax.set_ylabel('Metadynamics ΔG (kcal/mol)')
    ax.set_title(f'Metadynamics (RMSE={results["free_energy"]["metadynamics"]["rmse"]})')
    ax.legend(loc='lower right', fontsize=9)
    ax.set_aspect('equal')
    ax.set_xlim(lims2)
    ax.set_ylim(lims2)
    
    # Right: Method comparison bar chart
    ax = axes[2]
    methods = ['FEP', 'Metadynamics']
    fep_info = results["free_energy"]["fep"]
    metad_info = results["free_energy"]["metadynamics"]
    
    rmse_vals = [float(fep_info["rmse"].split()[0]), float(metad_info["rmse"].split()[0])]
    mae_vals = [float(fep_info["mae"].split()[0]), float(metad_info["mae"].split()[0])]
    
    x = np.arange(len(methods))
    width = 0.3
    ax.bar(x - width/2, rmse_vals, width, label='RMSE', color=COLORS[1], alpha=0.8)
    ax.bar(x + width/2, mae_vals, width, label='MAE', color=COLORS[2], alpha=0.8)
    ax.set_xlabel('Method')
    ax.set_ylabel('Error (kcal/mol)')
    ax.set_title('Method Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(methods)
    ax.legend()
    
    # Add R² annotations
    for i, (r2, tau) in enumerate([(fep_info["r_squared"], fep_info["kendall_tau"]),
                                     (metad_info["r_squared"], metad_info["kendall_tau"])]):
        ax.annotate(f'R²={r2}\nτ={tau}',
                    xy=(i, max(rmse_vals[i], mae_vals[i]) + 0.05),
                    ha='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/fig3_free_energy_comparison.png")
    plt.close()
    print("  Saved fig3_free_energy_comparison.png")


def fig4_gnn_performance(results, output_dir="figures"):
    """Figure 4: GNN training and prediction performance."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    curves = results["gnn"]["training_curves"]
    preds = results["gnn"]["predictions"]
    
    # Top-left: Training curves
    ax = axes[0, 0]
    epochs = curves["epoch"]
    ax.plot(epochs, curves["train_loss"], color=COLORS[1], label='Train Loss', alpha=0.8)
    ax.plot(epochs, curves["val_loss"], color=COLORS[4], label='Val Loss', alpha=0.8)
    ax.set_xlabel('Epoch')
    ax.set_ylabel('Loss')
    ax.set_title('Training and Validation Loss')
    ax.legend()
    ax.set_yscale('log')
    
    # Top-right: Validation RMSE
    ax = axes[0, 1]
    ax.plot(epochs, curves["val_rmse"], color=COLORS[2], linewidth=2)
    best_epoch = np.argmin(curves["val_rmse"])
    ax.axvline(x=best_epoch, color='gray', linestyle='--', alpha=0.5,
               label=f'Best: {curves["val_rmse"][best_epoch]:.3f} (epoch {best_epoch})')
    ax.set_xlabel('Epoch')
    ax.set_ylabel('RMSE (pKi)')
    ax.set_title('Validation RMSE')
    ax.legend()
    
    # Bottom-left: Predicted vs Experimental
    ax = axes[1, 0]
    pred_vals = [p[0] for p in preds]
    exp_vals = [p[1] for p in preds]
    unc_vals = [p[2] for p in preds]
    
    scatter = ax.scatter(exp_vals, pred_vals, c=unc_vals, cmap='viridis',
                         s=30, alpha=0.6, edgecolors='none')
    plt.colorbar(scatter, ax=ax, label='Uncertainty')
    
    lims = [3.5, 11]
    ax.plot(lims, lims, 'k--', alpha=0.3, label='y = x')
    ax.fill_between(lims, [l - 0.5 for l in lims], [l + 0.5 for l in lims],
                     alpha=0.08, color='gray')
    ax.set_xlabel('Experimental pKi')
    ax.set_ylabel('Predicted pKi')
    ax.set_title(f'Prediction Performance (R²={results["gnn"]["evaluation"]["r_squared"]})')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect('equal')
    ax.legend(loc='lower right')
    
    # Bottom-right: Error distribution
    ax = axes[1, 1]
    errors = [p[0] - p[1] for p in preds]
    ax.hist(errors, bins=40, color=COLORS[1], alpha=0.7, edgecolor='white')
    ax.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    ax.axvline(x=np.mean(errors), color='red', linestyle='--',
               label=f'Mean: {np.mean(errors):.3f}')
    ax.set_xlabel('Prediction Error (pKi)')
    ax.set_ylabel('Count')
    ax.set_title('Error Distribution')
    ax.legend()
    
    plt.suptitle('GNN Binding Affinity Prediction', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/fig4_gnn_performance.png")
    plt.close()
    print("  Saved fig4_gnn_performance.png")


def fig5_activity_cliffs(results, output_dir="figures"):
    """Figure 5: Activity cliff analysis and chemical space."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    cliff_data = results["activity_cliffs"]
    
    # Left: Chemical space PCA
    ax = axes[0]
    if cliff_data.get("pca_coords") and len(cliff_data["pca_coords"]) > 0:
        coords = np.array(cliff_data["pca_coords"])
        compounds = cliff_data["compounds_data"]
        pkis = [c[1] for c in compounds[:len(coords)]]
        scaffolds = [c[2] for c in compounds[:len(coords)]]
        
        scatter = ax.scatter(coords[:, 0], coords[:, 1], c=pkis,
                            cmap='viridis', s=40, alpha=0.7, edgecolors='white',
                            linewidth=0.5)
        plt.colorbar(scatter, ax=ax, label='pKi')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title('Chemical Space (PCA)')
    
    # Middle: SALI values for top cliffs
    ax = axes[1]
    top_cliffs = cliff_data["top_cliffs"][:10]
    if top_cliffs:
        sali_vals = [float(c["sali"]) for c in top_cliffs]
        sim_vals = [float(c["similarity"]) for c in top_cliffs]
        act_diffs = [float(c["activity_diff"].split()[0]) for c in top_cliffs]
        
        labels = [f'{c["pair"][0]}\nvs\n{c["pair"][1]}' for c in top_cliffs[:8]]
        
        ax.barh(range(len(labels)), sali_vals[:len(labels)],
                color=[plt.cm.RdYlGn_r(s / max(sali_vals)) for s in sali_vals[:len(labels)]],
                alpha=0.8)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel('SALI Score')
        ax.set_title('Top Activity Cliffs')
        ax.invert_yaxis()
    
    # Right: Similarity vs Activity difference
    ax = axes[2]
    if top_cliffs:
        all_salis = sali_vals
        ax.scatter(sim_vals, act_diffs, c=all_salis[:len(sim_vals)],
                  cmap='hot_r', s=100, alpha=0.7, edgecolors='black', linewidth=0.5)
        cbar = plt.colorbar(ax.collections[0], ax=ax, label='SALI')
        
        ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.4, label='Activity threshold')
        ax.axvline(x=0.7, color='gray', linestyle=':', alpha=0.4, label='Similarity threshold')
    ax.set_xlabel('Tanimoto Similarity')
    ax.set_ylabel('|ΔpKi|')
    ax.set_title('Activity Cliff Landscape')
    ax.legend(fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/fig5_activity_cliffs.png")
    plt.close()
    print("  Saved fig5_activity_cliffs.png")


def fig6_pareto_optimization(results, output_dir="figures"):
    """Figure 6: Multi-objective optimization and Pareto front."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    opt_data = results["optimization"]
    
    # Left: Pareto front (pKi vs Clearance)
    ax = axes[0]
    pareto_sols = opt_data["pareto_solutions"]
    if pareto_sols:
        pkis = [s["values"].get("pKi", 0) for s in pareto_sols]
        clearances = [s["values"].get("clearance", 0) for s in pareto_sols]
        selectivities = [s["values"].get("selectivity", 0) for s in pareto_sols]
        
        scatter = ax.scatter(pkis, clearances, c=selectivities, cmap='viridis',
                           s=80, alpha=0.8, edgecolors='black', linewidth=0.5)
        plt.colorbar(scatter, ax=ax, label='Selectivity')
        
        # Sort and draw Pareto front line
        pki_cl = sorted(zip(pkis, clearances), key=lambda x: x[0])
        ax.plot([x[0] for x in pki_cl], [x[1] for x in pki_cl],
                'r--', alpha=0.5, linewidth=1.5, label='Pareto front')
    
    ax.set_xlabel('pKi (maximize →)')
    ax.set_ylabel('Clearance (minimize →)')
    ax.set_title('Pareto Front: Potency vs Clearance')
    ax.legend()
    
    # Middle: Hypervolume convergence
    ax = axes[1]
    gen_history = opt_data["generation_history"]
    if gen_history:
        gens = [g["generation"] for g in gen_history]
        hvs = [g["hypervolume"] for g in gen_history]
        n_paretos = [g["n_pareto"] for g in gen_history]
        
        ax.plot(gens, hvs, color=COLORS[0], linewidth=2)
        ax.set_xlabel('Generation')
        ax.set_ylabel('Hypervolume')
        ax.set_title('Optimization Convergence')
        
        ax2 = ax.twinx()
        ax2.plot(gens, n_paretos, color=COLORS[2], linewidth=1.5, linestyle='--')
        ax2.set_ylabel('# Pareto-optimal', color=COLORS[2])
    
    # Right: Radar/spider chart of best solutions
    axes[2].remove()
    ax = fig.add_subplot(133, polar=True)
    if pareto_sols:
        objectives = list(pareto_sols[0]["values"].keys())
        n_obj = len(objectives)
        
        top3 = pareto_sols[:3]
        
        angles = np.linspace(0, 2 * np.pi, n_obj, endpoint=False).tolist()
        angles += angles[:1]
        
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        
        for idx, sol in enumerate(top3):
            values = []
            for obj_name in objectives:
                val = sol["values"].get(obj_name, 0)
                # Normalize to 0-1 range (approximate)
                all_vals = [s["values"].get(obj_name, 0) for s in pareto_sols]
                if max(all_vals) != min(all_vals):
                    norm_val = (val - min(all_vals)) / (max(all_vals) - min(all_vals))
                else:
                    norm_val = 0.5
                values.append(norm_val)
            values += values[:1]
            
            ax.plot(angles, values, 'o-', linewidth=1.5, alpha=0.7,
                   color=COLORS[idx], label=f'Solution {idx+1}')
            ax.fill(angles, values, alpha=0.1, color=COLORS[idx])
        
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(objectives, fontsize=8)
        ax.set_title('Top Pareto Solutions', pad=20)
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/fig6_pareto_optimization.png")
    plt.close()
    print("  Saved fig6_pareto_optimization.png")


def fig7_pipeline_overview(output_dir="figures"):
    """Figure 7: Pipeline architecture overview."""
    fig, ax = plt.subplots(1, 1, figsize=(16, 8))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 8)
    ax.axis('off')
    
    # Pipeline boxes
    boxes = [
        (1, 6, 'AlphaFold2\npLDDT Assessment', '#440154'),
        (4, 6, 'MD Simulation\nPose Refinement', '#31688e'),
        (7, 6, 'Free Energy\nFEP/Metadynamics', '#35b779'),
        (10, 6, 'GNN Predictor\nBinding Affinity', '#fde725'),
        (4, 3, 'Activity Cliff\nDetection', '#e76f51'),
        (7, 3, 'Chemical Space\nExploration', '#264653'),
        (10, 3, 'Multi-Objective\nOptimization', '#e76f51'),
    ]
    
    for x, y, label, color in boxes:
        rect = plt.Rectangle((x - 1.2, y - 0.7), 2.4, 1.4,
                              fill=True, facecolor=color, alpha=0.3,
                              edgecolor=color, linewidth=2, zorder=2)
        ax.add_patch(rect)
        ax.text(x, y, label, ha='center', va='center',
                fontsize=10, fontweight='bold', zorder=3)
    
    # Arrows
    arrow_params = dict(arrowstyle='->', color='gray', lw=2, mutation_scale=15)
    connections = [
        ((2.2, 6), (2.8, 6)),
        ((5.2, 6), (5.8, 6)),
        ((8.2, 6), (8.8, 6)),
        ((4, 5.3), (4, 3.7)),
        ((5.2, 3), (5.8, 3)),
        ((8.2, 3), (8.8, 3)),
        ((10, 5.3), (10, 3.7)),
    ]
    
    for (x1, y1), (x2, y2) in connections:
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=arrow_params)
    
    # Title
    ax.text(8, 7.5, 'AlphaFold2-Enhanced Binding Affinity Prediction Pipeline',
            ha='center', va='center', fontsize=14, fontweight='bold')
    
    # Data flow labels
    ax.text(3.5, 6.5, 'Structure', ha='center', fontsize=8, style='italic', color='gray')
    ax.text(6.5, 6.5, 'Refined\nPoses', ha='center', fontsize=8, style='italic', color='gray')
    ax.text(9.5, 6.5, 'ΔG values', ha='center', fontsize=8, style='italic', color='gray')
    
    plt.savefig(f"{output_dir}/fig7_pipeline_overview.png")
    plt.close()
    print("  Saved fig7_pipeline_overview.png")


def main():
    """Generate all figures."""
    os.makedirs("figures", exist_ok=True)
    
    print("Generating figures...")
    results = load_or_generate_data()
    
    fig1_plddt_profile(results)
    fig2_md_trajectory(results)
    fig3_free_energy_comparison(results)
    fig4_gnn_performance(results)
    fig5_activity_cliffs(results)
    fig6_pareto_optimization(results)
    fig7_pipeline_overview()
    
    print("\nAll figures generated successfully!")


if __name__ == "__main__":
    main()
