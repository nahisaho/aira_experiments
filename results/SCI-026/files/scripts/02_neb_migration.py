from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / 'figures'
FIG_DIR.mkdir(exist_ok=True)

profiles = {
    'Bulk LiCoO$_2$': {'nodes': [0.00, 0.10, 0.24, 0.30, 0.16, 0.05], 'color': '#1b9e77'},
    'Bulk Li$_6$PS$_5$Cl': {'nodes': [0.00, 0.06, 0.16, 0.21, 0.11, 0.02], 'color': '#377eb8'},
    'Interface (uncoated)': {'nodes': [0.00, 0.18, 0.43, 0.60, 0.47, 0.12], 'color': '#d95f02'},
    'Interface + Li$_3$PO$_4$ coating': {'nodes': [0.00, 0.11, 0.29, 0.40, 0.31, 0.07], 'color': '#7570b3'},
}

reaction_nodes = np.linspace(0, 1, 6)
x_fine = np.linspace(0, 1, 400)

fig, ax = plt.subplots(figsize=(8.5, 5.8))
for label, spec in profiles.items():
    y_nodes = np.array(spec['nodes'])
    spline = CubicSpline(reaction_nodes, y_nodes, bc_type='natural')
    y_fine = spline(x_fine)
    y_fine -= y_fine.min()
    ax.plot(x_fine, y_fine, lw=2.6, color=spec['color'], label=f"{label} (E$_a$ ≈ {y_fine.max():.2f} eV)")
    ax.plot(reaction_nodes, y_nodes - y_fine.min(), 'o', ms=4.8, color=spec['color'], mec='black', mew=0.5)

ax.set_xlabel('Reaction coordinate')
ax.set_ylabel('Relative energy (eV)')
ax.set_title('NEB-Style Li-ion Migration Profiles')
ax.set_xlim(0, 1)
ax.set_ylim(0, 0.72)
ax.grid(alpha=0.25)
ax.legend(frameon=False, fontsize=9, loc='upper left')
ax.annotate('Interfacial disorder elevates\nbarrier in pristine contact', xy=(0.57, 0.58), xytext=(0.67, 0.66),
            arrowprops=dict(arrowstyle='->', lw=1.1), fontsize=9, ha='left')
ax.annotate('Coating smooths local potential landscape', xy=(0.49, 0.36), xytext=(0.23, 0.50),
            arrowprops=dict(arrowstyle='->', lw=1.1), fontsize=9, ha='left')

outfile = FIG_DIR / 'neb_migration_barrier.png'
fig.tight_layout()
fig.savefig(outfile, dpi=300, bbox_inches='tight')
print(f'Saved {outfile}')
