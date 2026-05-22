"""
Module 1: Receptor-level ligand binding and signal initiation model
PTI: flg22-FLS2/BAK1 complex formation
ETI: AvrPita-Pita R protein recognition (rice blast)
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

# --- PTI: flg22-FLS2-BAK1 binding model ---
def pti_receptor_odes(y, t, p):
    flg22, FLS2, BAK1, C1, C2, Signal = y
    kon1, koff1, kon2, koff2, kcat = p
    dflg22 = -kon1*flg22*FLS2 + koff1*C1
    dFLS2  = -kon1*flg22*FLS2 + koff1*C1
    dBAK1  = -kon2*C1*BAK1 + koff2*C2
    dC1    = kon1*flg22*FLS2 - koff1*C1 - kon2*C1*BAK1 + koff2*C2
    dC2    = kon2*C1*BAK1 - koff2*C2 - kcat*C2
    dSignal = kcat*C2
    return [dflg22, dFLS2, dBAK1, dC1, dC2, dSignal]

# --- ETI: Guard model ---
def eti_receptor_odes(y, t, p):
    Eff, Target, ModTarget, NLR, NLR_act, ETI_sig = y
    kmod, kact, ksig, kdeg = p
    dEff = -kmod*Eff*Target
    dTarget = -kmod*Eff*Target
    dModTarget = kmod*Eff*Target - kact*ModTarget*NLR
    dNLR = -kact*ModTarget*NLR
    dNLR_act = kact*ModTarget*NLR - kdeg*NLR_act
    dETI_sig = ksig*NLR_act
    return [dEff, dTarget, dModTarget, dNLR, dNLR_act, dETI_sig]

t = np.linspace(0, 120, 1000)

pti_params = [0.01, 0.001, 0.005, 0.0005, 0.02]
pti_y0 = [10.0, 5.0, 3.0, 0, 0, 0]
eti_params = [0.008, 0.015, 0.03, 0.001]
eti_y0 = [5.0, 8.0, 0, 4.0, 0, 0]

pti_sol = odeint(pti_receptor_odes, pti_y0, t, args=(pti_params,))
eti_sol = odeint(eti_receptor_odes, eti_y0, t, args=(eti_params,))

# Dose-response
flg22_doses = np.logspace(-2, 2, 50)
max_signals = []
for dose in flg22_doses:
    y0 = [dose, 5.0, 3.0, 0, 0, 0]
    sol = odeint(pti_receptor_odes, y0, t, args=(pti_params,))
    max_signals.append(sol[-1, 5])

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

ax = axes[0, 0]
ax.plot(t, pti_sol[:, 0], 'b-', label='flg22 (free)', lw=2)
ax.plot(t, pti_sol[:, 3], 'g--', label='flg22:FLS2', lw=2)
ax.plot(t, pti_sol[:, 4], 'r-.', label='flg22:FLS2:BAK1', lw=2)
ax.plot(t, pti_sol[:, 5], 'k-', label='PTI Signal', lw=2.5)
ax.set_xlabel('Time (min)'); ax.set_ylabel('Concentration (a.u.)')
ax.set_title('PTI: flg22-FLS2/BAK1 Complex Formation')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

ax = axes[0, 1]
ax.plot(t, eti_sol[:, 0], 'b-', label='Effector (AvrPita)', lw=2)
ax.plot(t, eti_sol[:, 2], 'g--', label='Modified Target', lw=2)
ax.plot(t, eti_sol[:, 4], 'r-.', label='NLR Activated (Pita)', lw=2)
ax.plot(t, eti_sol[:, 5], 'k-', label='ETI Signal', lw=2.5)
ax.set_xlabel('Time (min)'); ax.set_ylabel('Concentration (a.u.)')
ax.set_title('ETI: Guard Model (AvrPita-Pita)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

ax = axes[1, 0]
ax.semilogx(flg22_doses, max_signals, 'bo-', lw=2, markersize=4)
ax.set_xlabel('flg22 Concentration (nM)'); ax.set_ylabel('Maximum PTI Signal')
ax.set_title('PTI Dose-Response Curve'); ax.grid(True, alpha=0.3)

ax = axes[1, 1]
pti_norm = pti_sol[:, 5] / max(pti_sol[:, 5]) if max(pti_sol[:, 5]) > 0 else pti_sol[:, 5]
eti_norm = eti_sol[:, 5] / max(eti_sol[:, 5]) if max(eti_sol[:, 5]) > 0 else eti_sol[:, 5]
ax.plot(t, pti_norm, 'b-', label='PTI Signal', lw=2)
ax.plot(t, eti_norm, 'r-', label='ETI Signal', lw=2)
ax.fill_between(t, pti_norm, alpha=0.15, color='blue')
ax.fill_between(t, eti_norm, alpha=0.15, color='red')
ax.set_xlabel('Time (min)'); ax.set_ylabel('Normalized Signal')
ax.set_title('PTI vs ETI Signal Dynamics')
ax.legend(); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/01_receptor_binding.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/01_receptor_binding.svg', bbox_inches='tight')
plt.close()

results = {
    'pti_half_max_time': float(t[np.argmin(np.abs(pti_norm - 0.5))]),
    'eti_half_max_time': float(t[np.argmin(np.abs(eti_norm - 0.5))]),
    'pti_max_signal': float(max(pti_sol[:, 5])),
    'eti_max_signal': float(max(eti_sol[:, 5])),
    'pti_EC50_approx': float(flg22_doses[np.argmin(np.abs(np.array(max_signals) - max(max_signals)/2))]),
    'pti_params': {'kon1': 0.01, 'koff1': 0.001, 'kon2': 0.005, 'koff2': 0.0005, 'kcat': 0.02},
    'eti_params': {'kmod': 0.008, 'kact': 0.015, 'ksig': 0.03, 'kdeg': 0.001}
}
with open('results/01_receptor_binding.json', 'w') as f:
    json.dump(results, f, indent=2)

print("Module 1 completed.")
print(f"PTI half-max time: {results['pti_half_max_time']:.1f} min")
print(f"ETI half-max time: {results['eti_half_max_time']:.1f} min")
print(f"PTI EC50 (approx): {results['pti_EC50_approx']:.2f} nM")
