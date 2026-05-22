"""
Module 4: Transcription regulatory network (WRKY/TGA TFs)
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

def tf_network_ode(y, t, p):
    W33, W70, W29, TGA, MYC2, PR1, PR2, PDF12, FRK1, W53 = y
    SA, JA, MAPK_pti, MAPK_eti, NPR1, k_deg, n = p
    dW33 = 0.5*MAPK_pti**n/(0.5**n+MAPK_pti**n) + 0.3*JA/(1+JA) - k_deg*W33
    dW70 = 0.8*SA*NPR1/(1+SA*NPR1) - 0.2*JA/(1+JA) - k_deg*W70
    dW29 = 0.6*MAPK_pti**n/(0.3**n+MAPK_pti**n) - k_deg*W29
    dTGA = 0.7*NPR1*SA/(1+NPR1*SA) - k_deg*TGA
    dMYC2 = 0.6*JA/(1+JA)*1/(1+W70) - k_deg*MYC2
    dPR1 = 0.8*TGA*NPR1/(1+TGA*NPR1) + 0.2*W70 - 0.05*PR1
    dPR2 = 0.5*TGA/(1+TGA) + 0.3*SA/(1+SA) - 0.05*PR2
    dPDF12 = 0.7*MYC2/(1+MYC2)*W33/(1+W33)*1/(1+W70) - 0.05*PDF12
    dFRK1 = 0.9*W29**n/(0.5**n+W29**n)*MAPK_pti/(1+MAPK_pti) - 0.05*FRK1
    dW53 = 0.4*MAPK_eti/(1+MAPK_eti) + 0.2*SA/(1+SA) - k_deg*W53
    return [dW33, dW70, dW29, dTGA, dMYC2, dPR1, dPR2, dPDF12, dFRK1, dW53]

t = np.linspace(0, 24, 1500)
y0 = [0.01]*10
tf_names = ['WRKY33', 'WRKY70', 'WRKY29', 'TGA2/5/6', 'MYC2']
gene_names = ['PR1', 'PR2', 'PDF1.2', 'FRK1', 'WRKY53']

scenarios = {
    'PTI': [0.5, 0.1, 0.8, 0.1, 0.3, 0.1, 2],
    'ETI': [1.0, 0.2, 0.3, 0.9, 0.8, 0.1, 2],
    'JA': [0.1, 1.0, 0.2, 0.1, 0.1, 0.1, 2]
}
solutions = {}
for name, params in scenarios.items():
    solutions[name] = odeint(tf_network_ode, y0, t, args=(params,))

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
plot_configs = [
    (axes[0,0], 'PTI', tf_names, range(5), 'PTI: TF Dynamics'),
    (axes[0,1], 'PTI', gene_names, range(5,10), 'PTI: Defense Genes'),
    (axes[0,2], 'ETI', tf_names, range(5), 'ETI: TF Dynamics'),
    (axes[1,0], 'ETI', gene_names, range(5,10), 'ETI: Defense Genes'),
    (axes[1,1], 'JA', gene_names, range(5,10), 'JA: Defense Genes'),
]
for ax, sc, names, idxs, title in plot_configs:
    for i, (idx, n) in enumerate(zip(idxs, names)):
        ax.plot(t, solutions[sc][:, idx], lw=2, label=n)
    ax.set_title(title); ax.set_xlabel('Time (h)'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

ax = axes[1, 2]
all_names = tf_names + gene_names
final_data = np.array([solutions[s][-1, :] for s in ['PTI', 'ETI', 'JA']])
im = ax.imshow(final_data, aspect='auto', cmap='viridis')
ax.set_xticks(range(len(all_names)))
ax.set_xticklabels(all_names, rotation=45, ha='right', fontsize=8)
ax.set_yticks([0, 1, 2]); ax.set_yticklabels(['PTI', 'ETI', 'JA'])
ax.set_title('Final Activity Heatmap')
plt.colorbar(im, ax=ax, label='Activity Level')

plt.tight_layout()
plt.savefig('figures/04_transcription_network.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/04_transcription_network.svg', bbox_inches='tight')
plt.close()

results = {sc: {n: float(solutions[sc][-1, i]) for i, n in enumerate(all_names)} for sc in scenarios}
results['network_edges'] = [
    ['WRKY33', 'PDF1.2', 'activation'], ['WRKY70', 'PR1', 'activation'],
    ['WRKY70', 'PDF1.2', 'repression'], ['WRKY70', 'MYC2', 'repression'],
    ['WRKY29', 'FRK1', 'activation'], ['TGA2/5/6', 'PR1', 'activation'],
    ['TGA2/5/6', 'PR2', 'activation'], ['MYC2', 'PDF1.2', 'activation'],
    ['NPR1', 'TGA2/5/6', 'activation'], ['MAPK', 'WRKY33', 'activation'],
    ['MAPK', 'WRKY29', 'activation'],
]
with open('results/04_transcription_network.json', 'w') as f:
    json.dump(results, f, indent=2)
print("Module 4 completed.")
print(f"PTI FRK1: {solutions['PTI'][-1, 8]:.4f}, ETI PR1: {solutions['ETI'][-1, 5]:.4f}")
