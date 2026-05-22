"""
Module 6: Rice blast resistance case study
Magnaporthe oryzae - Oryza sativa interaction
Integrated model: Pita recognition -> MAPK -> SA/JA -> defense genes -> HR
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

def rice_blast_model(y, t, p):
    """
    Integrated model for rice-M. oryzae interaction
    y = [AvrPita, Pita_act, ROS, MAPK, SA, JA, PR1a, PBZ1, POX, HR, Fungal_load]
    """
    Avr, Pita, ROS, MAPK, SA, JA, PR1a, PBZ1, POX, HR, Fungal = y
    (k_recog, k_ros, k_mapk_act, k_mapk_deg,
     k_sa_prod, k_sa_deg, k_ja_prod, k_ja_deg,
     k_pr1a, k_pbz1, k_pox,
     k_hr, k_hr_kill,
     k_fungal_grow, k_fungal_max,
     has_pita) = p

    # Recognition
    dAvr = -0.01 * Avr  # slow degradation
    recognition = k_recog * Avr * has_pita
    dPita = recognition / (1 + recognition) - 0.05 * Pita

    # ROS burst
    dROS = k_ros * Pita + 0.1 * MAPK - 0.2 * ROS

    # MAPK cascade
    dMAPK = k_mapk_act * (Pita + 0.3*ROS) / (1 + Pita + 0.3*ROS) - k_mapk_deg * MAPK

    # Hormone signaling
    dSA = k_sa_prod * MAPK / (1 + JA*0.5) - k_sa_deg * SA
    dJA = k_ja_prod * 0.3 / (1 + SA*2) - k_ja_deg * JA

    # Defense genes (rice-specific)
    dPR1a = k_pr1a * SA * MAPK / (1 + SA*MAPK) - 0.05 * PR1a
    dPBZ1 = k_pbz1 * SA / (1 + SA) - 0.05 * PBZ1  # Probenazole-inducible
    dPOX = k_pox * ROS * MAPK / (1 + ROS*MAPK) - 0.05 * POX  # Peroxidase

    # Hypersensitive response
    dHR = k_hr * ROS * Pita * MAPK - 0.01 * HR
    hr_effect = HR / (1 + HR)

    # Fungal growth
    dFungal = k_fungal_grow * Fungal * (1 - Fungal/k_fungal_max) * (1 - hr_effect) \
              - 0.1 * (PR1a + PBZ1 + POX) * Fungal / (1 + PR1a + PBZ1 + POX)

    return [dAvr, dPita, dROS, dMAPK, dSA, dJA, dPR1a, dPBZ1, dPOX, dHR, dFungal]

t = np.linspace(0, 96, 3000)  # 96 hours

# Resistant cultivar (has Pita)
params_R = [0.5, 0.8, 0.6, 0.1, 0.4, 0.05, 0.2, 0.05,
            0.6, 0.4, 0.5, 0.3, 0.5, 0.15, 100, 1.0]
y0_R = [5.0, 0, 0, 0, 0.1, 0.1, 0, 0, 0, 0, 1.0]

# Susceptible cultivar (no Pita)
params_S = list(params_R); params_S[-1] = 0.0
y0_S = list(y0_R)

# Partial resistance (reduced Pita expression)
params_P = list(params_R); params_P[-1] = 0.3

sol_R = odeint(rice_blast_model, y0_R, t, args=(params_R,))
sol_S = odeint(rice_blast_model, y0_S, t, args=(params_S,))
sol_P = odeint(rice_blast_model, y0_P := y0_R, t, args=(params_P,))

# R gene pyramiding: Pita + Pi9 (additive effect)
def rice_pyramid(y, t, p):
    vals = rice_blast_model(y, t, p)
    # Pi9 adds extra defense boost
    vals_list = list(vals)
    vals_list[3] *= 1.3  # enhanced MAPK
    vals_list[9] *= 1.2  # enhanced HR
    return vals_list

sol_pyramid = odeint(rice_pyramid, y0_R, t, args=(params_R,))

# Field simulation: disease progress curve
# Logistic model with resistance effect
def disease_progress(y, t, r, K, resistance_factor):
    severity = y[0]
    dsev = r * severity * (1 - severity/K) * (1 - resistance_factor)
    return [dsev]

t_field = np.linspace(0, 120, 500)  # days
sev_none = odeint(disease_progress, [0.01], t_field, args=(0.08, 100, 0))
sev_pita = odeint(disease_progress, [0.01], t_field, args=(0.08, 100, 0.7))
sev_pyramid = odeint(disease_progress, [0.01], t_field, args=(0.08, 100, 0.9))
sev_partial = odeint(disease_progress, [0.01], t_field, args=(0.08, 100, 0.4))

# R gene durability simulation
np.random.seed(42)
n_seasons = 20
resistance_breakdown = np.zeros((n_seasons, 3))
virulence_freq = [0.01, 0.01, 0.01]  # initial avr freq for Pita, Pi9, Pyramid
for season in range(n_seasons):
    for i, (freq, s) in enumerate(zip(virulence_freq,
                                       [0.15, 0.12, 0.03])):
        virulence_freq[i] = min(0.99, freq + s*freq*(1-freq) + np.random.normal(0, 0.01))
        resistance_breakdown[season, i] = 1 - virulence_freq[i]

# --- Figures ---
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Cellular response comparison
ax = axes[0, 0]
ax.plot(t, sol_R[:, 10], 'g-', label='Resistant (Pita+)', lw=2.5)
ax.plot(t, sol_S[:, 10], 'r-', label='Susceptible (pita)', lw=2.5)
ax.plot(t, sol_P[:, 10], 'orange', label='Partial (weak Pita)', lw=2)
ax.plot(t, sol_pyramid[:, 10], 'b--', label='Pyramid (Pita+Pi9)', lw=2)
ax.set_xlabel('Time (h)'); ax.set_ylabel('Fungal Load')
ax.set_title('M. oryzae Growth: R vs S Cultivars')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Defense response timeline (resistant)
ax = axes[0, 1]
ax.plot(t, sol_R[:, 2], 'r-', label='ROS burst', lw=2)
ax.plot(t, sol_R[:, 3], 'b-', label='MAPK', lw=2)
ax.plot(t, sol_R[:, 9], 'k-', label='HR', lw=2.5)
ax.plot(t, sol_R[:, 4], 'g--', label='SA', lw=1.5)
ax.set_xlabel('Time (h)'); ax.set_ylabel('Activity')
ax.set_title('Resistant: Defense Timeline')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Defense genes
ax = axes[0, 2]
ax.plot(t, sol_R[:, 6], 'r-', label='PR1a (R)', lw=2)
ax.plot(t, sol_R[:, 7], 'b-', label='PBZ1 (R)', lw=2)
ax.plot(t, sol_R[:, 8], 'g-', label='POX (R)', lw=2)
ax.plot(t, sol_S[:, 6], 'r--', label='PR1a (S)', lw=1.5, alpha=0.7)
ax.plot(t, sol_S[:, 7], 'b--', label='PBZ1 (S)', lw=1.5, alpha=0.7)
ax.set_xlabel('Time (h)'); ax.set_ylabel('Expression')
ax.set_title('Defense Gene Expression')
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

# Field disease progress
ax = axes[1, 0]
ax.plot(t_field, sev_none[:,0], 'r-', label='No resistance', lw=2)
ax.plot(t_field, sev_partial[:,0], 'orange', label='Partial (QTL)', lw=2)
ax.plot(t_field, sev_pita[:,0], 'b-', label='Pita', lw=2)
ax.plot(t_field, sev_pyramid[:,0], 'g-', label='Pyramid (Pita+Pi9)', lw=2)
ax.set_xlabel('Days after planting'); ax.set_ylabel('Disease Severity (%)')
ax.set_title('Field Disease Progress Curves')
ax.legend(); ax.grid(True, alpha=0.3)

# R gene durability
ax = axes[1, 1]
seasons = range(1, n_seasons+1)
ax.plot(seasons, resistance_breakdown[:,0], 'b-o', label='Pita alone', markersize=4, lw=2)
ax.plot(seasons, resistance_breakdown[:,1], 'g-s', label='Pi9 alone', markersize=4, lw=2)
ax.plot(seasons, resistance_breakdown[:,2], 'r-^', label='Pita+Pi9 pyramid', markersize=4, lw=2)
ax.set_xlabel('Growing Season'); ax.set_ylabel('Effective Resistance')
ax.set_title('R Gene Durability Over Seasons')
ax.legend(); ax.grid(True, alpha=0.3)
ax.set_ylim(0, 1.05)

# Summary heatmap
ax = axes[1, 2]
labels_t = ['0h', '6h', '12h', '24h', '48h', '72h', '96h']
time_points = [0, 6, 12, 24, 48, 72, 96]
comp_names = ['Pita', 'ROS', 'MAPK', 'SA', 'PR1a', 'PBZ1', 'HR', 'Fungal']
comp_idx = [1, 2, 3, 4, 6, 7, 9, 10]
heatmap_data = np.zeros((len(comp_names), len(time_points)))
for j, tp in enumerate(time_points):
    idx = np.argmin(np.abs(t - tp))
    for i, ci in enumerate(comp_idx):
        heatmap_data[i, j] = sol_R[idx, ci]

im = ax.imshow(heatmap_data, aspect='auto', cmap='YlOrRd')
ax.set_xticks(range(len(labels_t))); ax.set_xticklabels(labels_t)
ax.set_yticks(range(len(comp_names))); ax.set_yticklabels(comp_names)
ax.set_title('Resistant: Response Timeline')
plt.colorbar(im, ax=ax, label='Level')

plt.tight_layout()
plt.savefig('figures/06_rice_blast.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/06_rice_blast.svg', bbox_inches='tight')
plt.close()

results = {
    'resistant_fungal_96h': float(sol_R[-1, 10]),
    'susceptible_fungal_96h': float(sol_S[-1, 10]),
    'pyramid_fungal_96h': float(sol_pyramid[-1, 10]),
    'resistance_ratio': float(sol_S[-1, 10]/sol_R[-1, 10]) if sol_R[-1, 10] > 0.01 else float('inf'),
    'hr_peak_time': float(t[np.argmax(sol_R[:, 9])]),
    'ros_peak_time': float(t[np.argmax(sol_R[:, 2])]),
    'field_severity_120d': {
        'none': float(sev_none[-1, 0]),
        'partial': float(sev_partial[-1, 0]),
        'pita': float(sev_pita[-1, 0]),
        'pyramid': float(sev_pyramid[-1, 0])
    },
    'durability_season20': {
        'pita': float(resistance_breakdown[-1, 0]),
        'pi9': float(resistance_breakdown[-1, 1]),
        'pyramid': float(resistance_breakdown[-1, 2])
    }
}
with open('results/06_rice_blast.json', 'w') as f:
    json.dump(results, f, indent=2)

print("Module 6 completed.")
print(f"Fungal load at 96h - R: {results['resistant_fungal_96h']:.2f}, S: {results['susceptible_fungal_96h']:.2f}")
print(f"Resistance ratio: {results['resistance_ratio']:.1f}x")
print(f"Pyramid durability (season 20): {results['durability_season20']['pyramid']:.3f}")
