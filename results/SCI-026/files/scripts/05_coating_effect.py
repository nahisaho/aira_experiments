from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / 'figures'
FIG_DIR.mkdir(exist_ok=True)

fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8))

ax = axes[0]
labels = ['Uncoated', 'Li$_3$PO$_4$\n5 nm', 'Li$_3$PO$_4$\n10 nm', 'LiNbO$_3$', 'Li$_2$ZrO$_3$']
resistance = np.array([1200, 350, 250, 280, 320])
errors = np.array([120, 45, 30, 35, 40])
colors = ['#d73027', '#4575b4', '#74add1', '#66bd63', '#fdae61']
ax.bar(labels, resistance, yerr=errors, capsize=4, color=colors, edgecolor='black', linewidth=1.0)
ax.set_ylabel('Interface resistance (Ω·cm$^2$)')
ax.set_title('(a) Coating-induced resistance reduction')
ax.set_ylim(0, 1400)
ax.grid(axis='y', alpha=0.25)

ax = axes[1]
kb = 8.617333262145e-5
T = np.array([300, 325, 350, 375, 400, 425, 450])
x = 1000 / T
systems = {
    'Uncoated': {'Ea': 0.55, 'sigma0': 2.2e1, 'color': '#d73027'},
    'Li$_3$PO$_4$ coated': {'Ea': 0.35, 'sigma0': 5.5e1, 'color': '#4575b4'},
}
for name, spec in systems.items():
    sigma = spec['sigma0'] * np.exp(-spec['Ea'] / (kb * T))
    y = np.log10(sigma)
    ax.plot(x, y, '-o', color=spec['color'], lw=2.4, ms=4.8, label=name)
    coeff = np.polyfit(x, y, 1)
    x_mid = x[len(x) // 2]
    y_mid = np.polyval(coeff, x_mid)
    ax.text(x_mid + (0.03 if 'Uncoated' in name else -0.24), y_mid + (0.10 if 'Uncoated' in name else -0.18),
            f"E$_a$ ≈ {spec['Ea']:.2f} eV", color=spec['color'], fontsize=10, weight='bold')
ax.set_xlabel('1000/T (K$^{-1}$)')
ax.set_ylabel('log$_{10}$(interfacial conductivity / S cm$^{-1}$)')
ax.set_title('(b) Arrhenius response of interfacial transport')
ax.grid(alpha=0.25)
ax.legend(frameon=False, loc='lower left')
ax.invert_xaxis()

fig.suptitle('Effectiveness of Interfacial Coating Layers', fontsize=14, weight='bold')
fig.tight_layout(rect=[0, 0, 1, 0.95])
outfile = FIG_DIR / 'coating_effectiveness.png'
fig.savefig(outfile, dpi=300, bbox_inches='tight')
print(f'Saved {outfile}')
