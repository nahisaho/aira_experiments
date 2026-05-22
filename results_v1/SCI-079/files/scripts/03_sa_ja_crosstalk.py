"""
Module 3: SA/JA pathway crosstalk model
Mutual antagonism mediated by NPR1, WRKY70, GRX480
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

def sa_ja_crosstalk(y, t, p):
    SA, JA, NPR1, JAZ, PR1, PDF12, W70 = y
    (k_sa_prod, k_ja_prod, k_sa_deg, k_ja_deg,
     k_npr1_act, k_npr1_deg, k_jaz_deg, k_jaz_prod,
     k_pr1, k_pdf12, k_sa_inhibit_ja, k_ja_inhibit_sa,
     k_w70_act, k_w70_deg, SA_input, JA_input) = p
    dSA = k_sa_prod * SA_input / (1 + k_ja_inhibit_sa * JA) - k_sa_deg * SA
    dJA = k_ja_prod * JA_input / (1 + k_sa_inhibit_ja * SA * NPR1) - k_ja_deg * JA
    dNPR1 = k_npr1_act * SA**2 / (1 + SA**2) - k_npr1_deg * NPR1
    dJAZ = k_jaz_prod - k_jaz_deg * JA * JAZ
    dPR1 = k_pr1 * NPR1 * W70 / (1 + JAZ) - 0.1 * PR1
    dPDF12 = k_pdf12 * JA / (1 + NPR1 * W70 + JAZ*0.1) - 0.1 * PDF12
    dW70 = k_w70_act * NPR1 * SA / (1 + SA) - k_w70_deg * W70
    return [dSA, dJA, dNPR1, dJAZ, dPR1, dPDF12, dW70]

t = np.linspace(0, 48, 2000)
y0 = [0.1, 0.1, 0.1, 1.0, 0.0, 0.0, 0.1]

base_p = [1.0, 1.0, 0.05, 0.05, 0.5, 0.1, 0.3, 0.2,
          0.5, 0.5, 2.0, 2.0, 0.3, 0.1]

params_sa = base_p + [1.0, 0.0]
params_ja = base_p + [0.0, 1.0]
params_both = base_p + [1.0, 1.0]

sol_sa = odeint(sa_ja_crosstalk, y0, t, args=(params_sa,))
sol_ja = odeint(sa_ja_crosstalk, y0, t, args=(params_ja,))
sol_both = odeint(sa_ja_crosstalk, y0, t, args=(params_both,))

# Sequential SA->JA
sol_seq = np.zeros((len(t), 7)); sol_seq[0] = y0
for i in range(1, len(t)):
    dt = t[i] - t[i-1]
    p_seq = base_p + ([1.0, 0.0] if t[i] < 12 else [0.5, 1.0])
    sol_step = odeint(sa_ja_crosstalk, sol_seq[i-1], [0, dt], args=(p_seq,))
    sol_seq[i] = sol_step[-1]

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

for ax, sol, title in [(axes[0,0], sol_sa, 'SA-only (Biotrophic)'),
                        (axes[0,1], sol_ja, 'JA-only (Necrotrophic)'),
                        (axes[0,2], sol_both, 'SA+JA (Mixed)')]:
    ax.plot(t, sol[:, 0], 'r-', label='SA', lw=2)
    ax.plot(t, sol[:, 1], 'b-', label='JA', lw=2)
    ax.plot(t, sol[:, 4], 'r--', label='PR1', lw=2)
    ax.plot(t, sol[:, 5], 'b--', label='PDF1.2', lw=2)
    ax.set_title(title); ax.set_xlabel('Time (h)'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
axes[0,0].set_ylabel('Concentration')

ax = axes[1, 0]
ax.plot(t, sol_seq[:, 0], 'r-', label='SA', lw=2)
ax.plot(t, sol_seq[:, 1], 'b-', label='JA', lw=2)
ax.plot(t, sol_seq[:, 4], 'r--', label='PR1', lw=2)
ax.plot(t, sol_seq[:, 5], 'b--', label='PDF1.2', lw=2)
ax.axvline(x=12, color='gray', linestyle=':', alpha=0.7)
ax.set_title('Sequential SA→JA'); ax.set_xlabel('Time (h)'); ax.set_ylabel('Concentration')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

ax = axes[1, 1]
inhibit_strengths = np.linspace(0, 5, 30)
pr1_f, pdf12_f = [], []
for k_inh in inhibit_strengths:
    p_test = base_p + [1.0, 1.0]; p_test[10] = k_inh; p_test[11] = k_inh
    sol_test = odeint(sa_ja_crosstalk, y0, t, args=(p_test,))
    pr1_f.append(sol_test[-1, 4]); pdf12_f.append(sol_test[-1, 5])
ax.plot(inhibit_strengths, pr1_f, 'r-o', label='PR1 (SA)', markersize=4, lw=2)
ax.plot(inhibit_strengths, pdf12_f, 'b-s', label='PDF1.2 (JA)', markersize=4, lw=2)
ax.set_xlabel('Crosstalk Inhibition Strength'); ax.set_ylabel('Steady-state Expression')
ax.set_title('Crosstalk Strength Analysis'); ax.legend(); ax.grid(True, alpha=0.3)

ax = axes[1, 2]
ax.plot(t, sol_both[:, 2], 'g-', label='NPR1', lw=2)
ax.plot(t, sol_both[:, 3], 'm--', label='JAZ', lw=2)
ax.plot(t, sol_both[:, 6], color='orange', label='WRKY70', lw=2)
ax.set_xlabel('Time (h)'); ax.set_ylabel('Concentration')
ax.set_title('Key Regulators (SA+JA)'); ax.legend(); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/03_sa_ja_crosstalk.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/03_sa_ja_crosstalk.svg', bbox_inches='tight')
plt.close()

results = {
    'sa_only': {'PR1': float(sol_sa[-1, 4]), 'PDF12': float(sol_sa[-1, 5])},
    'ja_only': {'PR1': float(sol_ja[-1, 4]), 'PDF12': float(sol_ja[-1, 5])},
    'both': {'PR1': float(sol_both[-1, 4]), 'PDF12': float(sol_both[-1, 5])},
    'antagonism_index': float(1 - sol_both[-1, 5]/sol_ja[-1, 5]) if sol_ja[-1, 5] > 0 else 0,
    'sa_dominance': float(sol_both[-1, 4]/sol_sa[-1, 4]) if sol_sa[-1, 4] > 0 else 0
}
with open('results/03_sa_ja_crosstalk.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Module 3 completed.")
print(f"SA antagonism of JA: {results['antagonism_index']:.3f}")
print(f"SA dominance ratio: {results['sa_dominance']:.3f}")
