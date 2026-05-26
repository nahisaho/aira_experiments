#!/usr/bin/env python3
"""
Visualization module for perovskite screening results.
Generates all figures for report.md and paper.md.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from perovskite_screening import ScreeningPipeline, IonMigrationCalculator

# Style
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 150,
})

COLORS = {
    'Sn': '#2196F3', 'Ge': '#4CAF50', 'Bi': '#FF9800',
    'Sb': '#9C27B0', 'Ag': '#607D8B',
    'Single': '#3F51B5', 'Double': '#E91E63',
}

def get_b_color(formula):
    for el in ['Sn', 'Ge', 'Bi', 'Sb']:
        if el in formula:
            return COLORS[el]
    return '#999999'


def fig1_tolerance_stability(df, savepath):
    """Figure 1: Tolerance factor vs stability score."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # (a) Classical tolerance factor vs stability
    ax = axes[0]
    for _, row in df.iterrows():
        c = get_b_color(row['Formula'])
        marker = 's' if row['Type'] == 'Double' else 'o'
        ax.scatter(row['Tolerance Factor (t)'], row['Stability Score'],
                   c=c, marker=marker, s=60, alpha=0.8, edgecolors='k', linewidth=0.5)
    ax.axvspan(0.8, 1.0, alpha=0.1, color='green', label='Stable range')
    ax.set_xlabel('Goldschmidt Tolerance Factor (t)')
    ax.set_ylabel('Stability Score')
    ax.set_title('(a) Tolerance Factor vs Stability')
    ax.legend(loc='upper left', fontsize=8)

    # (b) Octahedral factor vs tolerance factor
    ax = axes[1]
    for _, row in df.iterrows():
        c = get_b_color(row['Formula'])
        marker = 's' if row['Type'] == 'Double' else 'o'
        ax.scatter(row['Tolerance Factor (t)'], row['Octahedral Factor (μ)'],
                   c=c, marker=marker, s=60, alpha=0.8, edgecolors='k', linewidth=0.5)
    ax.axvspan(0.8, 1.0, alpha=0.1, color='green')
    ax.axhspan(0.25, 0.70, alpha=0.1, color='blue', label='Stable μ range')
    ax.set_xlabel('Tolerance Factor (t)')
    ax.set_ylabel('Octahedral Factor (μ)')
    ax.set_title('(b) Structural Stability Map')
    ax.legend(loc='upper left', fontsize=8)

    # (c) New tolerance factor tau distribution
    ax = axes[2]
    sn_tau = df[df['Formula'].str.contains('Sn')]['New τ']
    ge_tau = df[df['Formula'].str.contains('Ge')]['New τ']
    bi_tau = df[df['Formula'].str.contains('Bi')]['New τ']

    data_to_plot = []
    labels_to_plot = []
    colors_to_plot = []
    for d, l, c in [(sn_tau, 'Sn-based', COLORS['Sn']),
                     (ge_tau, 'Ge-based', COLORS['Ge']),
                     (bi_tau, 'Bi-based', COLORS['Bi'])]:
        if len(d) > 0:
            data_to_plot.append(d.values)
            labels_to_plot.append(l)
            colors_to_plot.append(c)

    if data_to_plot:
        bp = ax.boxplot(data_to_plot, labels=labels_to_plot, patch_artist=True)
        for patch, color in zip(bp['boxes'], colors_to_plot):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)

    ax.axhline(y=4.18, color='r', linestyle='--', label='τ = 4.18 (stability limit)')
    ax.set_ylabel('New Tolerance Factor (τ)')
    ax.set_title('(c) τ Distribution by B-site')
    ax.legend(fontsize=8)

    # Custom legend for B-site elements
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['Sn'],
                   markersize=8, label='Sn-based'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['Ge'],
                   markersize=8, label='Ge-based'),
        plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=COLORS['Bi'],
                   markersize=8, label='Bi-based (double)'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3,
               bbox_to_anchor=(0.5, -0.02), fontsize=9)
    plt.tight_layout()
    plt.savefig(savepath, dpi=150)
    plt.close()
    print(f"  Saved: {savepath}")


def fig2_bandgap_analysis(df, savepath):
    """Figure 2: Bandgap prediction comparison and absorption."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # (a) DFT vs ML bandgap comparison
    ax = axes[0]
    for _, row in df.iterrows():
        c = get_b_color(row['Formula'])
        ax.scatter(row['Bandgap DFT (eV)'], row['Bandgap ML (eV)'],
                   c=c, s=60, alpha=0.8, edgecolors='k', linewidth=0.5)
    lims = [0.5, 3.5]
    ax.plot(lims, lims, 'k--', alpha=0.5, label='y = x')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel('DFT Bandgap (eV)')
    ax.set_ylabel('ML Bandgap (eV)')
    ax.set_title('(a) DFT vs ML Bandgap')
    ax.legend()

    # (b) Hybrid bandgap distribution
    ax = axes[1]
    sn_eg = df[df['Formula'].str.contains('Sn')]['Bandgap Hybrid (eV)']
    ge_eg = df[df['Formula'].str.contains('Ge')]['Bandgap Hybrid (eV)']
    bi_eg = df[df['Formula'].str.contains('Bi|Sb')]['Bandgap Hybrid (eV)']

    bins = np.linspace(0.5, 3.5, 20)
    if len(sn_eg) > 0:
        ax.hist(sn_eg, bins=bins, alpha=0.6, color=COLORS['Sn'], label='Sn-based')
    if len(ge_eg) > 0:
        ax.hist(ge_eg, bins=bins, alpha=0.6, color=COLORS['Ge'], label='Ge-based')
    if len(bi_eg) > 0:
        ax.hist(bi_eg, bins=bins, alpha=0.6, color=COLORS['Bi'], label='Bi/Sb-based')
    ax.axvspan(1.1, 1.5, alpha=0.15, color='green', label='Optimal (SQ)')
    ax.set_xlabel('Hybrid Bandgap (eV)')
    ax.set_ylabel('Count')
    ax.set_title('(b) Bandgap Distribution')
    ax.legend(fontsize=8)

    # (c) Bandgap vs PCE
    ax = axes[2]
    for _, row in df.iterrows():
        c = get_b_color(row['Formula'])
        marker = 's' if row['Type'] == 'Double' else 'o'
        ax.scatter(row['Bandgap Hybrid (eV)'], row['PCE (%)'],
                   c=c, marker=marker, s=60, alpha=0.8, edgecolors='k', linewidth=0.5)

    # SQ limit curve
    eg_range = np.linspace(0.5, 3.0, 100)
    sq_pce = []
    from perovskite_screening import DeviceSimulator
    ds = DeviceSimulator()
    for eg in eg_range:
        sq = ds.shockley_queisser_limit(eg)
        sq_pce.append(sq['PCE'])
    ax.plot(eg_range, sq_pce, 'r--', alpha=0.5, label='SQ Limit')
    ax.set_xlabel('Hybrid Bandgap (eV)')
    ax.set_ylabel('Simulated PCE (%)')
    ax.set_title('(c) Bandgap vs PCE')
    ax.legend()

    plt.tight_layout()
    plt.savefig(savepath, dpi=150)
    plt.close()
    print(f"  Saved: {savepath}")


def fig3_defect_migration(df, savepath):
    """Figure 3: Defect formation energy and ion migration."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # (a) Defect formation energy vs non-radiative loss
    ax = axes[0]
    for _, row in df.iterrows():
        c = get_b_color(row['Formula'])
        ax.scatter(row['Defect E_f (eV)'], row['Non-rad Loss (eV)'],
                   c=c, s=60, alpha=0.8, edgecolors='k', linewidth=0.5)
    ax.set_xlabel('Defect Formation Energy (eV)')
    ax.set_ylabel('Non-radiative V_OC Loss (eV)')
    ax.set_title('(a) Defects vs Recombination Loss')

    # (b) Ion migration barrier comparison
    ax = axes[1]
    categories = ['Sn-based', 'Ge-based', 'Bi/Sb-based']
    barriers_data = [
        df[df['Formula'].str.contains('Sn')]['Migration Barrier (eV)'].values,
        df[df['Formula'].str.contains('Ge')]['Migration Barrier (eV)'].values,
        df[df['Formula'].str.contains('Bi|Sb')]['Migration Barrier (eV)'].values,
    ]
    valid_data = [(d, l, c) for d, l, c in zip(
        barriers_data, categories,
        [COLORS['Sn'], COLORS['Ge'], COLORS['Bi']]
    ) if len(d) > 0]

    if valid_data:
        bp = ax.boxplot([d[0] for d in valid_data],
                        labels=[d[1] for d in valid_data], patch_artist=True)
        for patch, item in zip(bp['boxes'], valid_data):
            patch.set_facecolor(item[2])
            patch.set_alpha(0.6)

    ax.set_ylabel('Migration Barrier (eV)')
    ax.set_title('(b) Ion Migration Barriers')

    # (c) NEB-like energy profiles
    ax = axes[2]
    calc = IonMigrationCalculator()
    systems = [
        ('I⁻ in CsSnI₃', 0.28, COLORS['Sn']),
        ('I⁻ in CsGeI₃', 0.32, COLORS['Ge']),
        ('I⁻ in Cs₂AgBiI₆', 0.45, COLORS['Bi']),
    ]
    for label, barrier, color in systems:
        path = calc.neb_interpolated_path(barrier, n_images=20)
        ax.plot(path[:, 0], path[:, 1], '-o', color=color, markersize=3,
                label=f'{label} ({barrier:.2f} eV)')
    ax.set_xlabel('Reaction Coordinate')
    ax.set_ylabel('Energy (eV)')
    ax.set_title('(c) NEB Migration Energy Profiles')
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(savepath, dpi=150)
    plt.close()
    print(f"  Saved: {savepath}")


def fig4_device_performance(df, savepath):
    """Figure 4: Device simulation results."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    top10 = df.head(10)

    # (a) PCE bar chart for top 10
    ax = axes[0]
    colors = [get_b_color(f) for f in top10['Formula']]
    bars = ax.barh(range(len(top10)), top10['PCE (%)'].values, color=colors, alpha=0.8,
                   edgecolor='k', linewidth=0.5)
    ax.set_yticks(range(len(top10)))
    ax.set_yticklabels(top10['Formula'].values, fontsize=9)
    ax.set_xlabel('Simulated PCE (%)')
    ax.set_title('(a) Top 10 PCE')
    ax.invert_yaxis()

    # (b) J-V characteristics for top 3
    ax = axes[1]
    top3 = df.head(3)
    for idx, (_, row) in enumerate(top3.iterrows()):
        voc = row['V_OC (V)']
        jsc = row['J_SC (mA/cm²)']
        ff = row['FF']

        # Generate J-V curve (single diode model)
        v = np.linspace(0, voc * 1.05, 200)
        kbT = 0.02585
        n_ideal = 1.5
        j0 = jsc / (np.exp(voc / (n_ideal * kbT)) - 1)
        j = jsc - j0 * (np.exp(v / (n_ideal * kbT)) - 1)
        j = np.clip(j, 0, jsc * 1.1)

        c = get_b_color(row['Formula'])
        ax.plot(v, j, color=c, linewidth=2,
                label=f"{row['Formula']} (PCE={row['PCE (%)']:.1f}%)")

    ax.set_xlabel('Voltage (V)')
    ax.set_ylabel('Current Density (mA/cm²)')
    ax.set_title('(b) J-V Curves (Top 3)')
    ax.legend(fontsize=8)
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)

    # (c) Radar/spider chart - overall score breakdown for top 3
    ax = axes[2]
    categories_radar = ['Stability', 'Bandgap\nOptimality', 'Defect\nTolerance',
                        'Ion\nStability', 'PCE']
    n_cats = len(categories_radar)
    angles = np.linspace(0, 2 * np.pi, n_cats, endpoint=False).tolist()
    angles += angles[:1]

    for idx, (_, row) in enumerate(top3.iterrows()):
        # Calculate individual scores
        from perovskite_screening import CandidateRanker
        ranker = CandidateRanker()
        values = [
            row['Stability Score'],
            ranker.bandgap_optimality(row['Bandgap Hybrid (eV)']),
            ranker.defect_tolerance_score(row['Defect E_f (eV)']),
            ranker.ion_stability_score(row['Migration Barrier (eV)']),
            row['PCE (%)'] / 25.0,
        ]
        values += values[:1]
        c = get_b_color(row['Formula'])
        ax.plot(angles, values, 'o-', color=c, linewidth=2, markersize=4,
                label=row['Formula'])
        ax.fill(angles, values, alpha=0.1, color=c)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories_radar, fontsize=8)
    ax.set_ylim(0, 1.0)
    ax.set_title('(c) Multi-objective Scores')
    ax.legend(loc='upper right', fontsize=7)

    plt.tight_layout()
    plt.savefig(savepath, dpi=150)
    plt.close()
    print(f"  Saved: {savepath}")


def fig5_ranking_heatmap(df, savepath):
    """Figure 5: Comprehensive ranking heatmap."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    top15 = df.head(15)

    # (a) Property heatmap
    ax = axes[0]
    from perovskite_screening import CandidateRanker
    ranker = CandidateRanker()

    props = np.zeros((len(top15), 5))
    for i, (_, row) in enumerate(top15.iterrows()):
        props[i, 0] = row['Stability Score']
        props[i, 1] = ranker.bandgap_optimality(row['Bandgap Hybrid (eV)'])
        props[i, 2] = ranker.defect_tolerance_score(row['Defect E_f (eV)'])
        props[i, 3] = ranker.ion_stability_score(row['Migration Barrier (eV)'])
        props[i, 4] = row['PCE (%)'] / 25.0

    im = ax.imshow(props, cmap='YlOrRd', aspect='auto', vmin=0, vmax=1)
    ax.set_yticks(range(len(top15)))
    ax.set_yticklabels(top15['Formula'].values, fontsize=9)
    ax.set_xticks(range(5))
    ax.set_xticklabels(['Stability', 'Bandgap\nOpt.', 'Defect\nTol.',
                        'Ion\nStab.', 'PCE'], fontsize=9)
    ax.set_title('(a) Property Scores Heatmap')
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    # Add text annotations
    for i in range(props.shape[0]):
        for j in range(props.shape[1]):
            ax.text(j, i, f'{props[i,j]:.2f}', ha='center', va='center',
                    fontsize=7, color='black' if props[i,j] < 0.6 else 'white')

    # (b) Overall score bar chart
    ax = axes[1]
    colors = [get_b_color(f) for f in top15['Formula']]
    bars = ax.barh(range(len(top15)), top15['Overall Score'].values,
                   color=colors, alpha=0.8, edgecolor='k', linewidth=0.5)
    ax.set_yticks(range(len(top15)))
    ax.set_yticklabels(top15['Formula'].values, fontsize=9)
    ax.set_xlabel('Overall Score')
    ax.set_title('(b) Overall Ranking')
    ax.invert_yaxis()

    # Add score labels
    for i, v in enumerate(top15['Overall Score'].values):
        ax.text(v + 0.005, i, f'{v:.3f}', va='center', fontsize=8)

    plt.tight_layout()
    plt.savefig(savepath, dpi=150)
    plt.close()
    print(f"  Saved: {savepath}")


def fig6_workflow_diagram(savepath):
    """Figure 6: AiiDA workflow pipeline diagram."""
    fig, ax = plt.subplots(1, 1, figsize=(14, 6))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 6)
    ax.set_aspect('equal')
    ax.axis('off')

    steps = [
        (1, 4.5, 'Structure\nGeneration', '#E3F2FD'),
        (3.5, 4.5, 'DFT Geometry\nOptimization', '#E8F5E9'),
        (6, 4.5, 'Electronic\nStructure', '#FFF3E0'),
        (8.5, 4.5, 'ML Property\nPrediction', '#F3E5F5'),
        (1, 1.5, 'Defect\nCalculation', '#E8F5E9'),
        (3.5, 1.5, 'NEB Ion\nMigration', '#FFF3E0'),
        (6, 1.5, 'SCAPS-1D\nDevice Sim.', '#F3E5F5'),
        (8.5, 1.5, 'Multi-objective\nRanking', '#FFEBEE'),
        (11.5, 3, 'Candidate\nDatabase', '#E0F7FA'),
    ]

    for x, y, text, color in steps:
        box = FancyBboxPatch((x - 0.95, y - 0.55), 1.9, 1.1,
                              boxstyle="round,pad=0.1",
                              facecolor=color, edgecolor='black', linewidth=1.5)
        ax.add_patch(box)
        ax.text(x, y, text, ha='center', va='center', fontsize=8, fontweight='bold')

    # Arrows - top row
    for i in range(3):
        x1 = steps[i][0] + 1.0
        x2 = steps[i+1][0] - 1.0
        y = 4.5
        ax.annotate('', xy=(x2, y), xytext=(x1, y),
                     arrowprops=dict(arrowstyle='->', lw=1.5, color='#333'))

    # Arrows - down from row 1
    for i, j in [(2, 4), (3, 5)]:
        ax.annotate('', xy=(steps[j][0], steps[j][1] + 0.6),
                     xytext=(steps[i][0], steps[i][1] - 0.6),
                     arrowprops=dict(arrowstyle='->', lw=1.5, color='#333'))

    # Arrows - bottom row
    for i in range(4, 7):
        x1 = steps[i][0] + 1.0
        x2 = steps[i+1][0] - 1.0
        y = 1.5
        ax.annotate('', xy=(x2, y), xytext=(x1, y),
                     arrowprops=dict(arrowstyle='->', lw=1.5, color='#333'))

    # Arrows to final
    ax.annotate('', xy=(steps[8][0] - 1.0, steps[8][1] + 0.3),
                 xytext=(steps[3][0] + 1.0, steps[3][1]),
                 arrowprops=dict(arrowstyle='->', lw=1.5, color='#333', ls='--'))
    ax.annotate('', xy=(steps[8][0] - 1.0, steps[8][1] - 0.3),
                 xytext=(steps[7][0] + 1.0, steps[7][1]),
                 arrowprops=dict(arrowstyle='->', lw=1.5, color='#333'))

    ax.set_title('AiiDA/Fireworks Automated Screening Pipeline', fontsize=14, fontweight='bold', pad=10)
    plt.tight_layout()
    plt.savefig(savepath, dpi=150)
    plt.close()
    print(f"  Saved: {savepath}")


def main():
    """Generate all figures."""
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    fig_dir = os.path.join(base_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)

    print("Running screening pipeline...")
    pipeline = ScreeningPipeline()
    candidates = pipeline.generate_candidates()
    df, results = pipeline.screen(candidates)

    print(f"\nGenerating figures for {len(df)} candidates...")

    fig1_tolerance_stability(df, os.path.join(fig_dir, 'fig1_tolerance_stability.png'))
    fig2_bandgap_analysis(df, os.path.join(fig_dir, 'fig2_bandgap_analysis.png'))
    fig3_defect_migration(df, os.path.join(fig_dir, 'fig3_defect_migration.png'))
    fig4_device_performance(df, os.path.join(fig_dir, 'fig4_device_performance.png'))
    fig5_ranking_heatmap(df, os.path.join(fig_dir, 'fig5_ranking_heatmap.png'))
    fig6_workflow_diagram(os.path.join(fig_dir, 'fig6_workflow_pipeline.png'))

    print("\nAll figures generated successfully!")
    return df


if __name__ == '__main__':
    df = main()
