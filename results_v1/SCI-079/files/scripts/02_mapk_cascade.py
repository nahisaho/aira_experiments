"""
Module 2: MAPK cascade dynamics simulation
Three-tier cascade: MAPKKK -> MAPKK -> MAPK
PTI: MEKK1 -> MKK4/5 -> MPK3/6
ETI: stronger/sustained activation
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

def mapk_cascade(y, t, p):
    M3, M3p, M2, M2pp, M1, M1pp = y
    (Vm1, Km1, Vm2, Km2, Vm3, Km3,
     Vp1, Kp1, Vp2, Kp2, Vp3, Kp3, Signal) = p
    dM3  = -Vm1*Signal*M3/(Km1+M3) + Vp1*M3p/(Kp1+M3p)
    dM3p =  Vm1*Signal*M3/(Km1+M3) - Vp1*M3p/(Kp1+M3p)
    dM2   = -Vm2*M3p*M2/(Km2+M2) + Vp2*M2pp/(Kp2+M2pp)
    dM2pp =  Vm2*M3p*M2/(Km2+M2) - Vp2*M2pp/(Kp2+M2pp)
    dM1   = -Vm3*M2pp*M1/(Km3+M1) + Vp3*M1pp/(Kp3+M1pp)
    dM1pp =  Vm3*M2pp*M1/(Km3+M1) - Vp3*M1pp/(Kp3+M1pp)
    return [dM3, dM3p, dM2, dM2pp, dM1, dM1pp]

t = np.linspace(0, 60, 2000)
pti_params = [1.0, 0.5, 0.8, 0.5, 0.6, 0.5, 0.3, 0.5, 0.3, 0.5, 0.3, 0.5, 1.0]
eti_params = [1.5, 0.3, 1.2, 0.3, 1.0, 0.3, 0.15, 0.5, 0.15, 0.5, 0.15, 0.5, 1.5]
y0 = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0]

pti_sol = odeint(mapk_cascade, y0, t, args=(pti_params,))
eti_sol = odeint(mapk_cascade, y0, t, args=(eti_params,))

# Ultrasensitivity
signals = np.linspace(0, 2.0, 100)
steady_mapk = []
for sig in signals:
    p = list(pti_params); p[-1] = sig
    sol = odeint(mapk_cascade, y0, np.linspace(0, 200, 5000), args=(p,))
    steady_mapk.append(sol[-1, 5])

steady_arr = np.array(steady_mapk)
half_max = max(steady_arr) / 2
idx_half = np.argmin(np.abs(steady_arr - half_max))
if 0 < idx_half < len(signals)-1:
    dx = signals[idx_half+1] - signals[idx_half-1]
    dy = steady_arr[idx_half+1] - steady_arr[idx_half-1]
    nH = 4 * half_max * dy / (max(steady_arr) * dx) if dx > 0 else 1.0
else:
    nH = 1.0

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

ax = axes[0, 0]
ax.plot(t, pti_sol[:, 1], 'g-', label='MAPKKK-P (MEKK1)', lw=2)
ax.plot(t, pti_sol[:, 3], 'b--', label='MAPKK-PP (MKK4/5)', lw=2)
ax.plot(t, pti_sol[:, 5], 'r-', label='MAPK-PP (MPK3/6)', lw=2.5)
ax.set_xlabel('Time (min)'); ax.set_ylabel('Active fraction')
ax.set_title('PTI: MAPK Cascade Dynamics'); ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

ax = axes[0, 1]
ax.plot(t, eti_sol[:, 1], 'g-', label='MAPKKK-P', lw=2)
ax.plot(t, eti_sol[:, 3], 'b--', label='MAPKK-PP', lw=2)
ax.plot(t, eti_sol[:, 5], 'r-', label='MAPK-PP', lw=2.5)
ax.set_xlabel('Time (min)'); ax.set_ylabel('Active fraction')
ax.set_title('ETI: MAPK Cascade (Sustained)'); ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

ax = axes[1, 0]
ax.plot(signals, steady_arr, 'ko-', markersize=3, lw=2)
ax.axhline(y=half_max, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Input Signal Strength'); ax.set_ylabel('Steady-state MAPK-PP')
ax.set_title(f'Ultrasensitive Response (Hill coeff ≈ {nH:.2f})'); ax.grid(True, alpha=0.3)

ax = axes[1, 1]
ax.plot(t, pti_sol[:, 5], 'b-', label='PTI: MPK3/6', lw=2)
ax.plot(t, eti_sol[:, 5], 'r-', label='ETI: MAPK', lw=2)
ax.fill_between(t, pti_sol[:, 5], alpha=0.15, color='blue')
ax.fill_between(t, eti_sol[:, 5], alpha=0.15, color='red')
ax.set_xlabel('Time (min)'); ax.set_ylabel('Active MAPK fraction')
ax.set_title('PTI vs ETI: MAPK Activation'); ax.legend(); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/02_mapk_cascade.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/02_mapk_cascade.svg', bbox_inches='tight')
plt.close()

results = {
    'pti_mapk_peak': float(max(pti_sol[:, 5])),
    'pti_mapk_peak_time': float(t[np.argmax(pti_sol[:, 5])]),
    'eti_mapk_peak': float(max(eti_sol[:, 5])),
    'eti_mapk_peak_time': float(t[np.argmax(eti_sol[:, 5])]),
    'hill_coefficient': float(nH),
    'ec50_signal': float(signals[idx_half]),
    'amplification_ratio': float(max(eti_sol[:, 5]) / max(pti_sol[:, 5])) if max(pti_sol[:, 5]) > 0 else 0
}
with open('results/02_mapk_cascade.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Module 2 completed.")
for k, v in results.items(): print(f"  {k}: {v:.4f}")
