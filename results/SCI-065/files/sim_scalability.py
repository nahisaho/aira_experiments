#!/usr/bin/env python3
"""
Scalability Analysis: Batch → Perfusion → Continuous Bioreactor Design.
Compares throughput, quality, and cost across scales.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- Bioreactor configurations ---
configs = {
    'Static Well Plate': {
        'volume_mL': 2, 'organoids': 6, 'medium_change_h': 48,
        'O2_delivery': 'diffusion', 'shear_Pa': 0.0,
        'cost_per_organoid': 5.0, 'quality_score': 0.4,
        'scalability_score': 0.2, 'days_viable': 30,
        'type': 'batch'
    },
    'Spinner Flask': {
        'volume_mL': 125, 'organoids': 50, 'medium_change_h': 48,
        'O2_delivery': 'convection', 'shear_Pa': 0.05,
        'cost_per_organoid': 3.0, 'quality_score': 0.6,
        'scalability_score': 0.4, 'days_viable': 60,
        'type': 'batch'
    },
    'Orbital Shaker': {
        'volume_mL': 50, 'organoids': 24, 'medium_change_h': 48,
        'O2_delivery': 'convection', 'shear_Pa': 0.03,
        'cost_per_organoid': 4.0, 'quality_score': 0.55,
        'scalability_score': 0.3, 'days_viable': 60,
        'type': 'batch'
    },
    'Perfusion Bioreactor': {
        'volume_mL': 500, 'organoids': 200, 'medium_change_h': 0,  # continuous
        'O2_delivery': 'perfusion', 'shear_Pa': 0.02,
        'cost_per_organoid': 2.0, 'quality_score': 0.85,
        'scalability_score': 0.7, 'days_viable': 120,
        'type': 'perfusion'
    },
    'Microfluidic Perfusion': {
        'volume_mL': 0.5, 'organoids': 8, 'medium_change_h': 0,
        'O2_delivery': 'perfusion', 'shear_Pa': 0.01,
        'cost_per_organoid': 8.0, 'quality_score': 0.9,
        'scalability_score': 0.5, 'days_viable': 90,
        'type': 'perfusion'
    },
    'Stirred-Tank (CSTR)': {
        'volume_mL': 2000, 'organoids': 1000, 'medium_change_h': 0,
        'O2_delivery': 'forced', 'shear_Pa': 0.08,
        'cost_per_organoid': 1.5, 'quality_score': 0.75,
        'scalability_score': 0.9, 'days_viable': 90,
        'type': 'continuous'
    },
    'Automated Continuous': {
        'volume_mL': 5000, 'organoids': 5000, 'medium_change_h': 0,
        'O2_delivery': 'forced+perfusion', 'shear_Pa': 0.04,
        'cost_per_organoid': 1.0, 'quality_score': 0.80,
        'scalability_score': 1.0, 'days_viable': 120,
        'type': 'continuous'
    }
}

names = list(configs.keys())
n = len(names)

# --- Figure 1: Scalability comparison ---
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Throughput comparison
ax = axes[0, 0]
organoid_counts = [configs[n]['organoids'] for n in names]
colors_type = {'batch': '#4ECDC4', 'perfusion': '#FF6B6B', 'continuous': '#45B7D1'}
bar_colors = [colors_type[configs[n]['type']] for n in names]
bars = ax.barh(names, organoid_counts, color=bar_colors, alpha=0.8, edgecolor='black')
ax.set_xlabel('Organoids per Batch')
ax.set_title('(A) Throughput by Bioreactor Type')
ax.set_xscale('log')
# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=c, label=t.capitalize()) for t, c in colors_type.items()]
ax.legend(handles=legend_elements, loc='lower right')

# Cost vs Quality
ax = axes[0, 1]
costs = [configs[n]['cost_per_organoid'] for n in names]
qualities = [configs[n]['quality_score'] for n in names]
scatter_colors = [colors_type[configs[n]['type']] for n in names]
for i, name in enumerate(names):
    ax.scatter(costs[i], qualities[i], c=scatter_colors[i], s=200,
              edgecolors='black', zorder=5)
    ax.annotate(name, (costs[i], qualities[i]), textcoords="offset points",
               xytext=(5, 5), fontsize=7)
ax.set_xlabel('Cost per Organoid [USD]')
ax.set_ylabel('Quality Score')
ax.set_title('(B) Cost-Quality Trade-off')
ax.grid(True, alpha=0.3)

# Radar chart data
ax = axes[1, 0]
categories = ['Throughput', 'Quality', 'Scalability', 'Viability', 'Cost-eff.']
selected = ['Static Well Plate', 'Perfusion Bioreactor', 'Automated Continuous']
for name in selected:
    c = configs[name]
    values = [
        c['organoids'] / 5000,
        c['quality_score'],
        c['scalability_score'],
        c['days_viable'] / 120,
        1 - c['cost_per_organoid'] / 10
    ]
    ax.plot(categories, values, 'o-', linewidth=2, markersize=8, label=name)
ax.set_ylim(0, 1.1)
ax.set_ylabel('Normalized Score')
ax.set_title('(C) Multi-criteria Comparison')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Scale-up trajectory
ax = axes[1, 1]
scales = ['Lab\n(1-10)', 'Pilot\n(10-100)', 'Production\n(100-1000)', 'Industrial\n(1000+)']
batch_eff = [0.9, 0.6, 0.3, 0.1]
perfusion_eff = [0.7, 0.8, 0.85, 0.7]
continuous_eff = [0.5, 0.7, 0.9, 0.95]

x = np.arange(len(scales))
ax.plot(x, batch_eff, 'o-', linewidth=2.5, markersize=10, color='#4ECDC4', label='Batch')
ax.plot(x, perfusion_eff, 's-', linewidth=2.5, markersize=10, color='#FF6B6B', label='Perfusion')
ax.plot(x, continuous_eff, '^-', linewidth=2.5, markersize=10, color='#45B7D1', label='Continuous')
ax.set_xticks(x)
ax.set_xticklabels(scales)
ax.set_ylabel('Efficiency Score')
ax.set_title('(D) Scale-up Efficiency Trajectory')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/scalability_analysis.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/scalability_analysis.png")

# --- Figure 2: Bioreactor design schematic comparison ---
fig, ax = plt.subplots(figsize=(12, 6))
# Timeline/roadmap
phases_scale = {
    'Phase 1: R&D (Batch)': (0, 6, '#4ECDC4'),
    'Phase 2: Optimization (Perfusion)': (4, 12, '#FF6B6B'),
    'Phase 3: Scale-up (Continuous)': (10, 18, '#45B7D1'),
    'Phase 4: Production': (16, 24, '#FFD93D')
}
for i, (name, (start, end, color)) in enumerate(phases_scale.items()):
    ax.barh(i, end - start, left=start, height=0.6, color=color, alpha=0.8,
            edgecolor='black', linewidth=1.5)
    ax.text(start + (end - start) / 2, i, name, ha='center', va='center',
            fontsize=9, fontweight='bold')

ax.set_xlabel('Months')
ax.set_title('Bioreactor Development Roadmap')
ax.set_yticks([])
ax.set_xlim(-1, 25)
ax.grid(True, alpha=0.3, axis='x')
plt.tight_layout()
plt.savefig('figures/scalability_roadmap.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/scalability_roadmap.png")

# --- Print summary ---
print("\n=== Scalability Analysis Summary ===")
print(f"{'Configuration':<25} {'Type':<12} {'Organoids':>10} {'Cost/org':>10} {'Quality':>10}")
print("-" * 70)
for name in names:
    c = configs[name]
    print(f"{name:<25} {c['type']:<12} {c['organoids']:>10} ${c['cost_per_organoid']:>8.1f} {c['quality_score']:>10.2f}")
