"""
Visualization module for AutoLCA pipeline results.
Generates all figures for report and paper.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from lca_pipeline import AutoLCAPipeline

def set_style():
    plt.rcParams.update({
        'figure.figsize': (10, 6),
        'font.size': 11,
        'axes.titlesize': 13,
        'axes.labelsize': 11,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.dpi': 150,
        'savefig.dpi': 150,
        'savefig.bbox': 'tight',
    })

def plot_process_tree(pipeline, figdir):
    """Fig 1: Process tree visualization."""
    import networkx as nx
    fig, ax = plt.subplots(figsize=(14, 8))
    G = pipeline.tree_builder.graph
    procs = pipeline.tree_builder.processes

    pos = {
        'lithium_mining': (0, 4), 'nickel_mining': (2, 4),
        'cobalt_mining': (4, 4), 'manganese_mining': (6, 4),
        'graphite_production': (8, 4),
        'cathode_production': (3, 3), 'anode_production': (7, 3),
        'electrolyte_production': (1, 2), 'separator_production': (5, 2),
        'cell_assembly': (4, 1),
        'module_assembly': (4, 0),
        'pack_assembly': (4, -1),
    }

    cat_colors = {
        'raw_material': '#e74c3c',
        'component': '#3498db',
        'manufacturing': '#2ecc71',
    }

    node_colors = [cat_colors.get(procs[n].category, '#95a5a6') for n in G.nodes()]
    node_sizes = [max(800, procs[n].emissions.get('CO2', 1) * 120) for n in G.nodes()]

    labels = {n: procs[n].name.replace(' ', '\n')[:30] for n in G.nodes()}

    nx.draw_networkx_edges(G, pos, ax=ax, edge_color='#7f8c8d',
                           arrows=True, arrowsize=20, width=2, alpha=0.7)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.85, edgecolors='black', linewidths=1.5)
    nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=7, font_weight='bold')

    patches = [mpatches.Patch(color=c, label=l.replace('_', ' ').title())
               for l, c in cat_colors.items()]
    ax.legend(handles=patches, loc='upper right', framealpha=0.9)
    ax.set_title('EV Battery (NMC 811) Manufacturing Process Tree\n(Node size ∝ CO₂ emissions)', fontsize=14)
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(f'{figdir}/process_tree.png')
    plt.close(fig)
    print("  -> process_tree.png")

def plot_hotspot_analysis(results, figdir):
    """Fig 2: Hotspot contribution analysis."""
    df = pd.DataFrame(results['hotspot'])
    df = df.sort_values('co2_kg', ascending=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    cat_colors = {'raw_material': '#e74c3c', 'component': '#3498db', 'manufacturing': '#2ecc71'}
    colors = [cat_colors.get(r['category'], '#95a5a6') for _, r in df.iterrows()]

    ax1.barh(range(len(df)), df['co2_kg'], color=colors, edgecolor='black', linewidth=0.5)
    ax1.set_yticks(range(len(df)))
    ax1.set_yticklabels([n[:25] for n in df['process_name']], fontsize=9)
    ax1.set_xlabel('GHG Emissions (kg CO₂-eq)')
    ax1.set_title('Process-Level GHG Contributions')

    df_sorted = df.sort_values('contribution_pct', ascending=False)
    cumulative = df_sorted['contribution_pct'].cumsum().values
    ax2.bar(range(len(df_sorted)), df_sorted['contribution_pct'], color='#3498db',
            edgecolor='black', linewidth=0.5, label='Individual')
    ax2.plot(range(len(df_sorted)), cumulative, 'r-o', markersize=5, label='Cumulative')
    ax2.axhline(y=80, color='orange', linestyle='--', label='80% threshold')
    ax2.set_xticks(range(len(df_sorted)))
    ax2.set_xticklabels([n[:12] for n in df_sorted['process_name']], rotation=45, ha='right', fontsize=8)
    ax2.set_ylabel('Contribution (%)')
    ax2.set_title('Pareto Analysis of GHG Hotspots')
    ax2.legend()

    fig.tight_layout()
    fig.savefig(f'{figdir}/hotspot_analysis.png')
    plt.close(fig)
    print("  -> hotspot_analysis.png")

def plot_uncertainty(results, figdir):
    """Fig 3: Monte Carlo uncertainty distribution."""
    samples = results['mc_samples']
    stats = results['uncertainty']['monte_carlo']
    taylor = results['uncertainty']['taylor']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    ax1.hist(samples, bins=80, density=True, color='#3498db', alpha=0.7, edgecolor='black', linewidth=0.3)
    ax1.axvline(stats['mean'], color='red', linestyle='-', linewidth=2, label=f"Mean: {stats['mean']:.1f}")
    ax1.axvline(stats['p5'], color='orange', linestyle='--', linewidth=1.5, label=f"P5: {stats['p5']:.1f}")
    ax1.axvline(stats['p95'], color='orange', linestyle='--', linewidth=1.5, label=f"P95: {stats['p95']:.1f}")
    ax1.set_xlabel('Total GWP (kg CO₂-eq)')
    ax1.set_ylabel('Probability Density')
    ax1.set_title(f'Monte Carlo GWP Distribution (n={len(samples):,})')
    ax1.legend()

    labels = ['Mean', 'Std Dev', 'CV']
    mc_vals = [stats['mean'], stats['std'], stats['cv']]
    taylor_vals = [taylor['mean'], taylor['std'], taylor['cv']]
    x = np.arange(len(labels))
    w = 0.35
    ax2.bar(x - w/2, mc_vals, w, label='Monte Carlo', color='#3498db', edgecolor='black')
    ax2.bar(x + w/2, taylor_vals, w, label='Taylor Expansion', color='#e74c3c', edgecolor='black')
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels)
    ax2.set_ylabel('Value')
    ax2.set_title('MC vs Taylor Expansion Comparison')
    ax2.legend()

    fig.tight_layout()
    fig.savefig(f'{figdir}/uncertainty_analysis.png')
    plt.close(fig)
    print("  -> uncertainty_analysis.png")

def plot_scenarios(results, figdir):
    """Fig 4: Scenario comparison."""
    scenarios = results['scenarios']
    fig, ax = plt.subplots(figsize=(10, 6))

    names = list(scenarios.keys())
    values = list(scenarios.values())
    colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6', '#1abc9c']

    bars = ax.bar(range(len(names)), values, color=colors[:len(names)],
                  edgecolor='black', linewidth=0.8)

    baseline = values[0]
    for i, (bar, val) in enumerate(zip(bars, values)):
        reduction = (1 - val / baseline) * 100
        label = f'{val:.0f}\n({reduction:+.0f}%)' if i > 0 else f'{val:.0f}\n(baseline)'
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                label, ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([n.replace('_', '\n') for n in names], fontsize=9)
    ax.set_ylabel('Total GWP (kg CO₂-eq per battery pack)')
    ax.set_title('Scenario Comparison: EV Battery Manufacturing GWP')
    ax.set_ylim(0, max(values) * 1.25)

    fig.tight_layout()
    fig.savefig(f'{figdir}/scenario_comparison.png')
    plt.close(fig)
    print("  -> scenario_comparison.png")

def plot_scope3(results, figdir):
    """Fig 5: Scope 3 estimation results."""
    model_perf = results['scope3']['model_performance']
    estimates = results['scope3']['ev_battery_estimates']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    targets = list(model_perf.keys())
    rf_scores = [model_perf[t]['rf_r2_mean'] for t in targets]
    gb_scores = [model_perf[t]['gb_r2_mean'] for t in targets]
    x = np.arange(len(targets))
    w = 0.35
    ax1.bar(x - w/2, rf_scores, w, label='Random Forest', color='#3498db', edgecolor='black')
    ax1.bar(x + w/2, gb_scores, w, label='Gradient Boosting', color='#e74c3c', edgecolor='black')
    ax1.set_xticks(x)
    ax1.set_xticklabels([t.replace('_', '\n')[:15] for t in targets], fontsize=8)
    ax1.set_ylabel('R² Score (5-fold CV)')
    ax1.set_title('ML Model Performance for Scope 3 Estimation')
    ax1.legend()
    ax1.set_ylim(0, 1.0)
    ax1.axhline(y=0.8, color='green', linestyle='--', alpha=0.5, label='Good threshold')

    est_labels = [k.replace('_', ' ').title() for k in estimates.keys()]
    est_values = list(estimates.values())
    bars = ax2.bar(range(len(est_labels)), est_values, color='#2ecc71', edgecolor='black')
    ax2.set_xticks(range(len(est_labels)))
    ax2.set_xticklabels(est_labels, rotation=30, ha='right', fontsize=9)
    ax2.set_ylabel('Estimated Emissions (tCO₂-eq)')
    ax2.set_title('Scope 3 Estimates for EV Battery Manufacturer')
    for bar, val in zip(bars, est_values):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                 f'{val:,.0f}', ha='center', va='bottom', fontsize=9)

    fig.tight_layout()
    fig.savefig(f'{figdir}/scope3_analysis.png')
    plt.close(fig)
    print("  -> scope3_analysis.png")

def plot_matching_results(results, figdir):
    """Fig 6: Ecoinvent matching confidence."""
    details = results['matching']['details']
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    names = [d['process'][:20] for d in details.values()]
    sims = [d['similarity'] for d in details.values()]
    confs = [d['confidence'] for d in details.values()]
    conf_colors = {'high': '#2ecc71', 'medium': '#f39c12', 'low': '#e74c3c'}
    colors = [conf_colors.get(c, '#95a5a6') for c in confs]

    ax1.barh(range(len(names)), sims, color=colors, edgecolor='black', linewidth=0.5)
    ax1.set_yticks(range(len(names)))
    ax1.set_yticklabels(names, fontsize=8)
    ax1.set_xlabel('Cosine Similarity Score')
    ax1.set_title('Ecoinvent Auto-Matching Scores')
    ax1.axvline(x=0.5, color='green', linestyle='--', alpha=0.5, label='High threshold')
    ax1.axvline(x=0.3, color='orange', linestyle='--', alpha=0.5, label='Medium threshold')
    ax1.legend()

    stats = results['matching']['statistics']
    conf_labels = list(stats.keys())
    conf_vals = list(stats.values())
    ax2.pie(conf_vals, labels=[f'{l}\n({v})' for l, v in zip(conf_labels, conf_vals)],
            colors=[conf_colors.get(l.split('_')[0], '#95a5a6') for l in conf_labels],
            autopct='%1.0f%%', startangle=90)
    ax2.set_title('Matching Confidence Distribution')

    fig.tight_layout()
    fig.savefig(f'{figdir}/matching_results.png')
    plt.close(fig)
    print("  -> matching_results.png")

def plot_process_uncertainty_contribution(results, figdir):
    """Fig 7: Per-process uncertainty contribution."""
    mc_contribs = results['mc_process_contributions']

    proc_stats = {}
    for pid, samples in mc_contribs.items():
        proc_stats[pid] = {
            'mean': np.mean(samples),
            'std': np.std(samples),
            'cv': np.std(samples) / np.mean(samples) if np.mean(samples) > 0 else 0,
        }

    df = pd.DataFrame(proc_stats).T.sort_values('std', ascending=False)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    ax1.barh(range(len(df)), df['std'], color='#e74c3c', edgecolor='black', linewidth=0.5)
    ax1.set_yticks(range(len(df)))
    ax1.set_yticklabels([n.replace('_', ' ')[:20] for n in df.index], fontsize=9)
    ax1.set_xlabel('Standard Deviation (kg CO₂-eq)')
    ax1.set_title('Uncertainty Contribution by Process')

    ax2.scatter(df['mean'], df['cv'], s=df['std']*50+50, c='#3498db',
                alpha=0.7, edgecolors='black', linewidths=1)
    for i, (idx, row) in enumerate(df.iterrows()):
        ax2.annotate(idx.replace('_', ' ')[:15], (row['mean'], row['cv']),
                     fontsize=7, ha='center', va='bottom')
    ax2.set_xlabel('Mean Emission (kg CO₂-eq)')
    ax2.set_ylabel('Coefficient of Variation')
    ax2.set_title('Mean vs. Uncertainty (bubble size = std)')

    fig.tight_layout()
    fig.savefig(f'{figdir}/uncertainty_contribution.png')
    plt.close(fig)
    print("  -> uncertainty_contribution.png")


def main():
    set_style()
    figdir = os.path.join(os.path.dirname(__file__), '..', 'figures')
    os.makedirs(figdir, exist_ok=True)

    print("Running AutoLCA pipeline...")
    pipeline = AutoLCAPipeline()
    results = pipeline.run()

    print("\nGenerating figures...")
    plot_process_tree(pipeline, figdir)
    plot_hotspot_analysis(results, figdir)
    plot_uncertainty(results, figdir)
    plot_scenarios(results, figdir)
    plot_scope3(results, figdir)
    plot_matching_results(results, figdir)
    plot_process_uncertainty_contribution(results, figdir)
    print("\nAll figures generated!")

    return results

if __name__ == '__main__':
    results = main()
