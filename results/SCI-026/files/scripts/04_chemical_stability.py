from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / 'figures'
FIG_DIR.mkdir(exist_ok=True)

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))

ax = axes[0]
reactions = [
    'Li$_6$PS$_5$Cl →\nLi$_2$S + Li$_3$P + LiCl',
    'LCO + LPS\nmutual reaction',
    'Surface oxidation\nproducts',
    'Interfacial\noxysulfides',
    'Coating-mediated\npassivation',
]
energies = np.array([-0.18, -0.26, -0.11, -0.21, -0.04])
colors = ['#66c2a5' if e > -0.1 else '#fc8d62' for e in energies]
bars = ax.bar(range(len(reactions)), energies, color=colors, edgecolor='black', linewidth=1.0)
ax.axhline(0, color='black', lw=1.0)
for bar, val in zip(bars, energies):
    ax.text(bar.get_x() + bar.get_width() / 2, val - 0.015, f'{val:.2f}', ha='center', va='top', fontsize=9, weight='bold')
ax.set_xticks(range(len(reactions)), reactions)
ax.set_ylabel('Reaction energy (eV/atom)')
ax.set_title('(a) Interfacial decomposition thermodynamics')
ax.set_ylim(-0.32, 0.06)
ax.grid(axis='y', alpha=0.25)

ax = axes[1]
mu = np.array([-4.3, -3.8, -3.3, -2.8, -2.3, -2.0, -1.7, -1.2])
hull = np.array([0.00, -0.06, -0.10, -0.08, -0.05, -0.03, -0.01, 0.02])
ax.plot(mu, hull, '-o', color='#2c7fb8', lw=2.2, ms=4.5)
ax.axvspan(-2.1, -1.7, color='#a6d96a', alpha=0.45, label='Li$_6$PS$_5$Cl stability window')
ax.axvspan(-4.2, -3.0, color='#74add1', alpha=0.25, label='LiCoO$_2$ operating range')
ax.axvspan(-3.0, -2.1, color='#f46d43', alpha=0.22, label='Instability gap')
ax.text(-1.9, 0.012, '1.7–2.1 V', ha='center', va='bottom', fontsize=9, weight='bold')
ax.text(-3.6, -0.115, '3.0–4.2 V', ha='center', va='bottom', fontsize=9, weight='bold', color='#2c7fb8')
ax.set_xlabel(r'Li chemical potential, $\mu_{Li}$ (eV)')
ax.set_ylabel('Relative grand potential / hull distance (eV/atom)')
ax.set_title('(b) Voltage-chemical potential stability map')
ax.grid(alpha=0.25)
secax = ax.secondary_xaxis('top', functions=(lambda m: -m, lambda v: -v))
secax.set_xlabel('Voltage vs Li/Li$^+$ (V)')
ax.legend(frameon=False, fontsize=8, loc='lower left')

fig.suptitle('Chemical Stability Analysis of the Li$_6$PS$_5$Cl / LiCoO$_2$ Interface', fontsize=14, weight='bold')
fig.tight_layout(rect=[0, 0, 1, 0.95])
outfile = FIG_DIR / 'chemical_stability.png'
fig.savefig(outfile, dpi=300, bbox_inches='tight')
print(f'Saved {outfile}')
