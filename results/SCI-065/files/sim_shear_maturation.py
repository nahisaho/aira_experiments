#!/usr/bin/env python3
"""
Shear Stress - Tissue Maturation Relationship Model.
Models the effect of mechanical stimulation on brain organoid development.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --- Shear stress response model ---
# Based on literature: organoid viability & maturation follow biphasic response
# Low shear: beneficial (nutrient mixing), High shear: detrimental (damage)

def maturation_index(tau, tau_opt, sigma, alpha, M_max):
    """
    Biphasic maturation model:
    M(tau) = M_max * (tau/tau_opt) * exp(1 - tau/tau_opt) * exp(-alpha*(tau-tau_opt)^2/sigma^2)
    """
    x = tau / tau_opt
    return M_max * x * np.exp(1 - x) * np.exp(-alpha * ((tau - tau_opt)**2) / sigma**2)

def viability(tau, tau_crit, n):
    """Cell viability decreasing with shear: V(tau) = 1 / (1 + (tau/tau_crit)^n)"""
    return 1.0 / (1.0 + (tau / tau_crit)**n)

def neural_marker_expression(tau, tau_opt, k1, k2):
    """Neural maturation markers (MAP2, NeuN) vs shear stress."""
    return k1 * np.exp(-((np.log(tau/tau_opt))**2) / (2*k2**2))

# --- Parameters ---
tau_range = np.linspace(0.001, 1.0, 500)  # Pa

# Maturation index
tau_opt = 0.05      # Optimal shear stress [Pa]
sigma = 0.15
alpha = 2.0
M_max = 1.0

# Viability
tau_crit = 0.3      # Critical shear for damage [Pa]
n_hill = 4          # Hill coefficient

# Neural markers
k1_MAP2, k2_MAP2 = 0.85, 0.8
k1_NeuN, k2_NeuN = 0.75, 0.7
k1_SYN, k2_SYN = 0.70, 0.9

M = maturation_index(tau_range, tau_opt, sigma, alpha, M_max)
V = viability(tau_range, tau_crit, n_hill)
MAP2 = neural_marker_expression(tau_range, tau_opt, k1_MAP2, k2_MAP2)
NeuN = neural_marker_expression(tau_range, tau_opt, k1_NeuN, k2_NeuN)
SYN = neural_marker_expression(tau_range, tau_opt, k1_SYN, k2_SYN)

# --- Figure: Shear stress response ---
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Maturation index
ax = axes[0, 0]
ax.plot(tau_range*1000, M, 'b-', linewidth=2.5, label='Maturation Index')
ax.axvline(x=tau_opt*1000, color='g', linestyle='--', alpha=0.7, label=f'τ_opt = {tau_opt*1000:.0f} mPa')
ax.fill_between(tau_range*1000, 0, M, alpha=0.1, color='blue')
ax.set_xlabel('Shear Stress [mPa]')
ax.set_ylabel('Maturation Index (normalized)')
ax.set_title('(A) Organoid Maturation vs Shear Stress')
ax.legend()
ax.grid(True, alpha=0.3)

# Viability
ax = axes[0, 1]
ax.plot(tau_range*1000, V*100, 'r-', linewidth=2.5, label='Cell Viability')
ax.axvline(x=tau_crit*1000, color='k', linestyle='--', alpha=0.7, label=f'τ_crit = {tau_crit*1000:.0f} mPa')
ax.axhspan(0, 50, alpha=0.1, color='red')
ax.set_xlabel('Shear Stress [mPa]')
ax.set_ylabel('Viability [%]')
ax.set_title('(B) Cell Viability vs Shear Stress')
ax.legend()
ax.grid(True, alpha=0.3)

# Neural markers
ax = axes[1, 0]
ax.plot(tau_range*1000, MAP2, 'g-', linewidth=2, label='MAP2')
ax.plot(tau_range*1000, NeuN, 'm-', linewidth=2, label='NeuN')
ax.plot(tau_range*1000, SYN, 'c-', linewidth=2, label='Synaptophysin')
ax.axvline(x=tau_opt*1000, color='k', linestyle='--', alpha=0.5, label=f'τ_opt = {tau_opt*1000:.0f} mPa')
ax.set_xlabel('Shear Stress [mPa]')
ax.set_ylabel('Marker Expression (normalized)')
ax.set_title('(C) Neural Maturation Markers vs Shear Stress')
ax.legend()
ax.grid(True, alpha=0.3)

# Combined quality metric
ax = axes[1, 1]
Q_total = M * V * (MAP2 + NeuN + SYN) / 3
Q_total_norm = Q_total / np.max(Q_total) if np.max(Q_total) > 0 else Q_total
ax.plot(tau_range*1000, Q_total_norm, 'k-', linewidth=2.5, label='Overall Quality')
optimal_idx = np.argmax(Q_total_norm)
ax.axvline(x=tau_range[optimal_idx]*1000, color='g', linestyle='--',
           label=f'Optimal τ = {tau_range[optimal_idx]*1000:.1f} mPa')
ax.fill_between(tau_range*1000, 0, Q_total_norm, alpha=0.1, color='green')
ax.set_xlabel('Shear Stress [mPa]')
ax.set_ylabel('Quality Score (normalized)')
ax.set_title('(D) Combined Organoid Quality Score')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/shear_maturation.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/shear_maturation.png")

# --- Time-course maturation under different shear conditions ---
days = np.arange(0, 91)
shear_conditions = [0.01, 0.03, 0.05, 0.10, 0.30]  # Pa

fig, ax = plt.subplots(figsize=(10, 6))
for tau in shear_conditions:
    M_val = maturation_index(tau, tau_opt, sigma, alpha, M_max)
    V_val = viability(tau, tau_crit, n_hill)
    # Logistic maturation curve modulated by shear
    maturation_curve = M_val * V_val / (1 + np.exp(-0.08 * (days - 30)))
    ax.plot(days, maturation_curve, linewidth=2, label=f'τ = {tau*1000:.0f} mPa')

ax.set_xlabel('Culture Day')
ax.set_ylabel('Maturation Score')
ax.set_title('Temporal Maturation Under Different Shear Conditions')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('figures/temporal_maturation.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/temporal_maturation.png")

# --- Print summary ---
print("\n=== Shear-Maturation Model Summary ===")
print(f"Optimal shear stress: {tau_range[optimal_idx]*1000:.1f} mPa")
print(f"Critical shear for damage: {tau_crit*1000:.0f} mPa")
print(f"Safe operating range: {tau_opt*0.5*1000:.0f} - {tau_opt*3*1000:.0f} mPa")
print(f"Peak maturation index at τ_opt: {M_max:.2f}")
print(f"Viability at τ_crit: {viability(tau_crit, tau_crit, n_hill)*100:.1f}%")
