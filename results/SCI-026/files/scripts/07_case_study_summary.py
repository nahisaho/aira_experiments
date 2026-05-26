from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / 'figures'
FIG_DIR.mkdir(exist_ok=True)


def gauss(x, mu, sigma, amp):
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.4))

ax = axes[0, 0]
E = np.linspace(-8, 4, 1200)
co = gauss(E, -1.7, 0.55, 3.3) + gauss(E, 1.1, 0.45, 2.4)
o = gauss(E, -4.6, 0.75, 2.7) + gauss(E, -2.8, 0.65, 2.1)
s = gauss(E, -3.2, 0.60, 1.8) + gauss(E, -0.8, 0.55, 1.1)
p = gauss(E, -5.3, 0.70, 1.3)
interface = gauss(E, -0.2, 0.18, 0.75) + gauss(E, 0.35, 0.15, 0.55)
for y, label, color in [
    (co, 'Co-3d', '#d73027'),
    (o, 'O-2p', '#4575b4'),
    (s, 'S-3p', '#66bd63'),
    (p, 'P-3s', '#fdae61'),
]:
    ax.plot(E, y, lw=2.0, label=label, color=color)
ax.fill_between(E, interface, color='#984ea3', alpha=0.35, label='Interface states')
ax.axvline(0, color='black', ls='--', lw=1.1)
ax.set_xlim(-8, 4)
ax.set_xlabel('Energy - E$_F$ (eV)')
ax.set_ylabel('PDOS (arb. units)')
ax.set_title('(a) Interface PDOS and band alignment')
ax.legend(frameon=False, ncol=2, fontsize=8)
ax.grid(alpha=0.2)

ax = axes[0, 1]
x = np.linspace(-20, 20, 600)
rho = 0.016 * np.sin(0.95 * x) * np.exp(-(x / 9.5) ** 2) + 0.035 * np.exp(-((x + 3) / 2.4) ** 2) - 0.042 * np.exp(-((x - 2) / 2.8) ** 2)
ax.plot(x, rho, color='#5e3c99', lw=2.3)
ax.fill_between(x, 0, rho, where=rho >= 0, color='#1b9e77', alpha=0.35, label='Electron accumulation')
ax.fill_between(x, 0, rho, where=rho < 0, color='#d95f02', alpha=0.35, label='Electron depletion')
ax.axhline(0, color='black', lw=1.0)
ax.axvline(0, color='black', ls='--', lw=1.0)
ax.set_xlabel('Distance from interface (Å)')
ax.set_ylabel('Δρ (e Å$^{-3}$)')
ax.set_title('(b) Charge density difference profile')
ax.legend(frameon=False, fontsize=8)
ax.grid(alpha=0.2)

ax = axes[1, 0]
kb = 8.617333262145e-5
T = np.array([300, 325, 350, 375, 400, 425, 450])
invT = 1000 / T
bulk_sigma = 8.0e1 * np.exp(-0.20 / (kb * T))
int_sigma = 4.0e1 * np.exp(-0.42 / (kb * T))
ax.plot(invT, np.log10(bulk_sigma), '-o', color='#1b9e77', lw=2.2, label='Bulk Li$_6$PS$_5$Cl')
ax.plot(invT, np.log10(int_sigma), '-s', color='#d95f02', lw=2.2, label='Interface region')
ax.set_xlabel('1000/T (K$^{-1}$)')
ax.set_ylabel('log$_{10}$(ionic conductivity / S cm$^{-1}$)')
ax.set_title('(c) Temperature-dependent ionic transport')
ax.legend(frameon=False, fontsize=8, loc='lower left')
ax.grid(alpha=0.2)
ax.invert_xaxis()

ax = axes[1, 1]
cycles = np.arange(0, 201)
uncoated = 160 - 0.25 * cycles - 6 * np.log1p(cycles / 15)
coated = 162 - 0.10 * cycles - 2.5 * np.log1p(cycles / 25)
ax.plot(cycles, uncoated, color='#d73027', lw=2.4, label='Uncoated interface')
ax.plot(cycles, coated, color='#4575b4', lw=2.4, label='Coated interface')
ax.fill_between(cycles, coated - 2.5, coated + 2.5, color='#4575b4', alpha=0.18)
ax.fill_between(cycles, uncoated - 3.5, uncoated + 3.5, color='#d73027', alpha=0.14)
ax.set_xlabel('Cycle number')
ax.set_ylabel('Discharge capacity (mAh g$^{-1}$)')
ax.set_title('(d) Cycling performance prediction')
ax.legend(frameon=False, fontsize=8)
ax.grid(alpha=0.2)
ax.set_xlim(0, 200)
ax.set_ylim(95, 166)

fig.suptitle('Case Study Summary: Li$_6$PS$_5$Cl / LiCoO$_2$ Interface', fontsize=15, weight='bold')
fig.tight_layout(rect=[0, 0, 1, 0.96])
outfile = FIG_DIR / 'case_study_summary.png'
fig.savefig(outfile, dpi=300, bbox_inches='tight')
print(f'Saved {outfile}')
