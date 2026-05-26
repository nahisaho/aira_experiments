from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / 'figures'
FIG_DIR.mkdir(exist_ok=True)

steps = [
    ('Structure\nPreparation', 'analysis', (0.06, 0.62)),
    ('DFT\nOptimization\n(VASP)', 'dft', (0.24, 0.62)),
    ('Interface\nConstruction', 'analysis', (0.42, 0.62)),
    ('NEB\nCalculation', 'dft', (0.60, 0.62)),
    ('AIMD\nSimulation', 'md', (0.78, 0.62)),
    ('Space Charge\nAnalysis', 'analysis', (0.17, 0.20)),
    ('Coating\nScreening', 'analysis', (0.44, 0.20)),
    ('Performance\nPrediction', 'analysis', (0.71, 0.20)),
]
colors = {'dft': '#91bfdb', 'md': '#99d594', 'analysis': '#fdae61'}

fig, ax = plt.subplots(figsize=(12, 4.8))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

box_w, box_h = 0.14, 0.18
for label, kind, (x, y) in steps:
    patch = FancyBboxPatch((x, y), box_w, box_h, boxstyle='round,pad=0.02,rounding_size=0.03',
                           facecolor=colors[kind], edgecolor='black', linewidth=1.3)
    ax.add_patch(patch)
    ax.text(x + box_w / 2, y + box_h / 2, label, ha='center', va='center', fontsize=10, weight='bold')

connections = [
    ((0.20, 0.71), (0.24, 0.71)),
    ((0.38, 0.71), (0.42, 0.71)),
    ((0.56, 0.71), (0.60, 0.71)),
    ((0.74, 0.71), (0.78, 0.71)),
    ((0.85, 0.62), (0.78, 0.38)),
    ((0.78, 0.29), (0.58, 0.29)),
    ((0.51, 0.29), (0.31, 0.29)),
]
for start, end in connections:
    arrow = FancyArrowPatch(start, end, arrowstyle='-|>', mutation_scale=16, lw=1.5, color='dimgray')
    ax.add_patch(arrow)

ax.text(0.50, 0.93, 'VASP/LAMMPS Workflow for Interface Resistance Prediction', ha='center', fontsize=15, weight='bold')
ax.text(0.10, 0.06, 'Blue: DFT calculations', color=colors['dft'], fontsize=10, weight='bold')
ax.text(0.38, 0.06, 'Green: Molecular dynamics', color=colors['md'], fontsize=10, weight='bold')
ax.text(0.67, 0.06, 'Orange: Data analysis / screening', color=colors['analysis'], fontsize=10, weight='bold')

outfile = FIG_DIR / 'simulation_workflow.png'
fig.tight_layout()
fig.savefig(outfile, dpi=300, bbox_inches='tight')
print(f'Saved {outfile}')
