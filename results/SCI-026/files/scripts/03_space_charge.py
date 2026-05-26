from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / 'figures'
FIG_DIR.mkdir(exist_ok=True)

x = np.linspace(-30, 30, 800)
potential = 0.42 / (1 + np.exp((x - 2.0) / 3.0)) + 0.03 * np.exp(-((x + 8.0) / 6.0) ** 2)
potential += 0.02 * np.sin(0.22 * x) * np.exp(-(x / 22.0) ** 2)
potential -= potential.min()

fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharex=True)

ax = axes[0]
ax.plot(x, potential, color='#2c7fb8', lw=2.5)
ax.axvline(0, color='black', ls='--', lw=1.1)
ax.axvspan(0, 12, color='#fee08b', alpha=0.35, label='Electrolyte depletion region')
ax.text(-22, potential.max() * 0.90, 'LiCoO$_2$', fontsize=11, weight='bold')
ax.text(14, potential.max() * 0.90, 'Li$_6$PS$_5$Cl$', fontsize=11, weight='bold')
ax.annotate('Potential drop ≈ 0.4 V', xy=(5, 0.18), xytext=(-16, 0.33),
            arrowprops=dict(arrowstyle='->', lw=1.1), fontsize=9)
ax.set_xlabel('Distance from interface (Å)')
ax.set_ylabel('Electrostatic potential (V)')
ax.set_title('(a) Potential profile across space-charge layer')
ax.grid(alpha=0.25)
ax.legend(frameon=False, loc='upper right')

ax = axes[1]
widths = {
    '(104)||(100), w ≈ 8 Å': (8.0, '#1b9e77'),
    '(003)||(110), w ≈ 12 Å': (12.0, '#d95f02'),
    '(012)||(111), w ≈ 16 Å': (16.0, '#7570b3'),
}
for label, (w, color) in widths.items():
    electrolyte_depletion = 1 - 0.68 * np.exp(-((x - 5.0) / w) ** 2) * (x > -2)
    cathode_accumulation = 1 + 0.12 * np.exp(-((x + 4.0) / 6.0) ** 2) * (x < 8)
    conc = electrolyte_depletion * cathode_accumulation
    ax.plot(x, conc, lw=2.3, color=color, label=label)
ax.axvline(0, color='black', ls='--', lw=1.1)
ax.set_xlabel('Distance from interface (Å)')
ax.set_ylabel('Normalized Li concentration')
ax.set_title('(b) Li carrier redistribution')
ax.set_ylim(0.2, 1.22)
ax.grid(alpha=0.25)
ax.legend(frameon=False, fontsize=8, loc='lower right')

fig.suptitle('Space-Charge Layer Analysis at the Li$_6$PS$_5$Cl / LiCoO$_2$ Interface', fontsize=14, weight='bold')
fig.tight_layout(rect=[0, 0, 1, 0.95])
outfile = FIG_DIR / 'space_charge_layer.png'
fig.savefig(outfile, dpi=300, bbox_inches='tight')
print(f'Saved {outfile}')
