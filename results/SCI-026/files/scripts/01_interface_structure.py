from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / 'figures'
FIG_DIR.mkdir(exist_ok=True)


def draw_stack(ax, x0, width, layers, title, plane_label, colors):
    y = 0.0
    for i, (name, h) in enumerate(layers):
        rect = Rectangle((x0, y), width, h, facecolor=colors[i % len(colors)], edgecolor='black', lw=1.2)
        ax.add_patch(rect)
        ax.text(x0 + width / 2, y + h / 2, name, ha='center', va='center', fontsize=9, weight='bold')
        y += h
    ax.text(x0 + width / 2, y + 0.28, title, ha='center', va='bottom', fontsize=12, weight='bold')
    ax.text(x0 + width / 2, -0.42, plane_label, ha='center', va='top', fontsize=10, color='dimgray')
    ax.plot([x0, x0 + width], [y, y], color='black', ls='--', lw=1.1)
    return y


fig = plt.figure(figsize=(12, 6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.28)

ax0 = fig.add_subplot(gs[0, 0])
ax0.set_xlim(0, 10)
ax0.set_ylim(-0.8, 5.8)
ax0.axis('off')

lco_layers = [('O', 0.45), ('Li', 0.35), ('CoO$_2$', 0.75), ('Li', 0.35), ('CoO$_2$', 0.75), ('Li', 0.35), ('O', 0.45)]
lpscl_layers = [('Cl', 0.50), ('Li cage', 0.55), ('PS$_4$', 0.85), ('Li cage', 0.55), ('S', 0.50)]

h1 = draw_stack(ax0, 1.0, 2.5, lco_layers, r'Layered LiCoO$_2$ (R$\bar{3}$m)', '(104) surface', ['#d73027', '#fee090', '#4575b4'])
h2 = draw_stack(ax0, 6.0, 2.5, lpscl_layers, r'Argyrodite Li$_6$PS$_5$Cl (F$\bar{4}$3m)', '(100) surface', ['#74add1', '#fdae61', '#abd9e9'])

ax0.text(5.0, 5.2, 'Interface model', ha='center', va='center', fontsize=14, weight='bold')
ax0.text(5.0, 4.8, '(104) LiCoO$_2$  ||  (100) Li$_6$PS$_5$Cl', ha='center', va='center', fontsize=11)
ax0.add_patch(FancyArrowPatch((3.7, 2.3), (5.8, 2.3), arrowstyle='<->', mutation_scale=16, lw=1.5, color='black'))
ax0.text(4.75, 2.55, 'interfacial contact', ha='center', fontsize=10)
ax0.annotate('Li diffusion bottleneck\nat pristine interface', xy=(4.8, 1.7), xytext=(4.8, 0.6),
             ha='center', arrowprops=dict(arrowstyle='-|>', lw=1.2, color='dimgray'), fontsize=9, color='dimgray')

ax0.text(1.0, 5.45, 'c-axis slabs', fontsize=10, color='#555555')
ax0.text(6.0, 5.45, 'Li-rich argyrodite framework', fontsize=10, color='#555555')

ax1 = fig.add_subplot(gs[0, 1])
orientations = ['(104)||(100)', '(003)||(110)', '(012)||(111)']
mismatch = np.array([3.2, 5.8, 7.1])
colors = ['#1b9e77', '#d95f02', '#7570b3']
bars = ax1.bar(orientations, mismatch, color=colors, edgecolor='black', linewidth=1.0)
ax1.axhline(5.0, color='firebrick', ls='--', lw=1.2, label='5% heuristic threshold')
for bar, val in zip(bars, mismatch):
    ax1.text(bar.get_x() + bar.get_width() / 2, val + 0.18, f'{val:.1f}%', ha='center', va='bottom', fontsize=10, weight='bold')
ax1.set_ylabel('Lattice mismatch (%)')
ax1.set_title('Interface orientation matching')
ax1.set_ylim(0, 8.4)
ax1.legend(frameon=False, loc='upper left')
ax1.grid(axis='y', alpha=0.25)
ax1.text(0, 7.8, 'Best registry:', color='#1b9e77', fontsize=10, weight='bold')
ax1.text(0, 7.45, '(104)LiCoO$_2$ || (100)Li$_6$PS$_5$Cl', fontsize=9)

fig.suptitle('Li$_6$PS$_5$Cl / LiCoO$_2$ Interface Structure Schematic', fontsize=15, weight='bold', y=0.98)
fig.subplots_adjust(left=0.05, right=0.98, bottom=0.08, top=0.90, wspace=0.28)
outfile = FIG_DIR / 'interface_structure.png'
fig.savefig(outfile, dpi=300, bbox_inches='tight')
print(f'Saved {outfile}')
