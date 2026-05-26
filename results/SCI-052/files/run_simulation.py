#!/usr/bin/env python3
"""
Main simulation script for the microkinetic modeling framework.
Generates all figures for the report and paper.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from src.rate_constants import (compute_rate_constant, eyring_rate_constant,
                                 wigner_tunneling_correction, arrhenius_parameters,
                                 kB, h, eV_to_J)
from src.adsorption import (langmuir_isotherm, competitive_langmuir,
                             temkin_isotherm, fractal_isotherm,
                             coverage_dependent_binding_energy)
from src.lateral_interactions import LateralInteractionModel, create_ft_interaction_matrix
from src.fischer_tropsch import (compute_ft_rate_constants, run_ft_simulation,
                                  temperature_study, degree_of_rate_control_analysis,
                                  FT_ENERGETICS)
from src.reactor_models import SimpleMicroKineticModel, PFR, CSTR

os.makedirs('figures', exist_ok=True)

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'figure.dpi': 150,
    'savefig.bbox': 'tight',
    'savefig.dpi': 150,
})

print("=" * 60)
print("Microkinetic Modeling Framework - Fischer-Tropsch Synthesis")
print("=" * 60)

# ============================================================
# Figure 1: Arrhenius plots with tunneling corrections
# ============================================================
print("\n[1] Generating Arrhenius plots...")

T_range = np.linspace(400, 800, 200)
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Rate constants for key steps
step_data = [
    ('CO dissociation', 1.60, 450),
    ('C hydrogenation', 0.75, 1100),
    ('CH₃ hydrogenation', 0.85, 950),
    ('O hydrogenation', 1.00, 900),
]

for name, Ea, nu in step_data:
    k_no_tunnel = [eyring_rate_constant(Ea, T) for T in T_range]
    k_wigner = [compute_rate_constant(Ea, T, nu, tunneling='wigner') for T in T_range]
    axes[0].semilogy(1000/T_range, k_no_tunnel, '--', alpha=0.5)
    axes[0].semilogy(1000/T_range, k_wigner, '-', label=name)

axes[0].set_xlabel('1000/T [K⁻¹]')
axes[0].set_ylabel('Rate constant k [s⁻¹]')
axes[0].set_title('(a) Arrhenius plots for FT elementary steps')
axes[0].legend(fontsize=9)
axes[0].grid(True, alpha=0.3)

# Panel B: Tunneling correction factors
for name, Ea, nu in step_data:
    kappa = [wigner_tunneling_correction(nu, T) for T in T_range]
    axes[1].plot(T_range, kappa, label=name)

axes[1].set_xlabel('Temperature [K]')
axes[1].set_ylabel('Wigner correction factor κ')
axes[1].set_title('(b) Tunneling correction factors')
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/arrhenius_plots.png')
plt.close()
print("   -> figures/arrhenius_plots.png")

# ============================================================
# Figure 2: Adsorption isotherms comparison
# ============================================================
print("\n[2] Generating adsorption isotherm plots...")

P_range = np.linspace(0.01, 50, 500)
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: Langmuir isotherms at different K
for K, label in [(0.1, 'K=0.1'), (0.5, 'K=0.5'), (1.0, 'K=1.0'), (5.0, 'K=5.0')]:
    theta = langmuir_isotherm(P_range, K)
    axes[0].plot(P_range, theta, label=label)

axes[0].set_xlabel('Pressure [bar]')
axes[0].set_ylabel('Coverage θ')
axes[0].set_title('(a) Langmuir isotherm')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Panel B: Temkin isotherms
for alpha, label in [(0.5, 'α=0.5'), (1.0, 'α=1.0'), (2.0, 'α=2.0')]:
    P_temkin, theta_temkin = temkin_isotherm(P_range, K0=1.0, alpha=alpha)
    axes[1].plot(P_temkin, theta_temkin, label=label)

axes[1].set_xlabel('Pressure [bar]')
axes[1].set_ylabel('Coverage θ')
axes[1].set_title('(b) Temkin isotherm')
axes[1].set_xlim(0, 50)
axes[1].legend()
axes[1].grid(True, alpha=0.3)

# Panel C: Fractal isotherms
for D, label in [(2.0, 'D=2.0 (smooth)'), (2.3, 'D=2.3'), (2.6, 'D=2.6'), (2.9, 'D=2.9')]:
    theta = fractal_isotherm(P_range/50, 1.0, D)
    axes[2].plot(P_range, theta, label=label)

axes[2].set_xlabel('Pressure [bar]')
axes[2].set_ylabel('Coverage θ')
axes[2].set_title('(c) Fractal surface isotherm')
axes[2].legend()
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/adsorption_isotherms.png')
plt.close()
print("   -> figures/adsorption_isotherms.png")

# ============================================================
# Figure 3: Lateral interaction effects
# ============================================================
print("\n[3] Generating lateral interaction plots...")

species, epsilon = create_ft_interaction_matrix()
lim = LateralInteractionModel(species, epsilon, coordination_number=4)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Coverage-dependent binding energy
theta_range = np.linspace(0, 0.8, 100)
E0_CO = -1.35  # eV
for model_type, label in [('linear', 'Linear'), ('quadratic', 'Quadratic'), ('piecewise', 'Piecewise')]:
    E = coverage_dependent_binding_energy(E0_CO, theta_range, 0.4, model=model_type)
    axes[0].plot(theta_range, E, label=label)

axes[0].set_xlabel('Coverage θ_CO')
axes[0].set_ylabel('Binding energy [eV]')
axes[0].set_title('(a) Coverage-dependent CO binding energy')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Panel B: Interaction matrix heatmap
im = axes[1].imshow(epsilon, cmap='RdBu_r', vmin=-0.05, vmax=0.20)
axes[1].set_xticks(range(len(species)))
axes[1].set_yticks(range(len(species)))
axes[1].set_xticklabels(species, rotation=45, ha='right', fontsize=9)
axes[1].set_yticklabels(species, fontsize=9)
axes[1].set_title('(b) Lateral interaction matrix [eV]')
plt.colorbar(im, ax=axes[1], label='ε [eV]')
for i in range(len(species)):
    for j in range(len(species)):
        axes[1].text(j, i, f'{epsilon[i,j]:.2f}', ha='center', va='center', fontsize=7)

plt.tight_layout()
plt.savefig('figures/lateral_interactions.png')
plt.close()
print("   -> figures/lateral_interactions.png")

# ============================================================
# Figure 4: Surface coverage evolution
# ============================================================
print("\n[4] Running FT simulation at T=500K...")

results_500 = run_ft_simulation(T=500, P_total=20.0, H2_CO_ratio=2.0)
cov_sol = results_500['coverage_solution']

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

species_labels = ['CO*', 'H*', 'C*', 'O*', 'CH*', 'CH₂*', 'CH₃*', 'OH*']
colors = plt.cm.tab10(np.linspace(0, 1, 8))

for i, (label, color) in enumerate(zip(species_labels, colors)):
    axes[0].plot(cov_sol.t, np.maximum(cov_sol.y[i], 0), label=label, color=color)

axes[0].set_xlabel('Time [s]')
axes[0].set_ylabel('Surface coverage θ')
axes[0].set_title('(a) Coverage evolution at T=500K, P=20bar')
axes[0].legend(ncol=2, fontsize=8)
axes[0].grid(True, alpha=0.3)
axes[0].set_xscale('log')

# Panel B: Steady-state coverage bar chart
theta_ss = np.maximum(results_500['coverages'], 0)
theta_star = max(1.0 - np.sum(theta_ss), 0)
all_labels = species_labels + ['*']
all_theta = list(theta_ss) + [theta_star]

bars = axes[1].bar(range(len(all_labels)), all_theta, color=list(colors) + ['gray'])
axes[1].set_xticks(range(len(all_labels)))
axes[1].set_xticklabels(all_labels, rotation=45, ha='right')
axes[1].set_ylabel('Coverage θ')
axes[1].set_title('(b) Steady-state surface coverages')
axes[1].grid(True, alpha=0.3, axis='y')

for bar, val in zip(bars, all_theta):
    if val > 0.01:
        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                    f'{val:.3f}', ha='center', va='bottom', fontsize=8)

plt.tight_layout()
plt.savefig('figures/surface_coverages.png')
plt.close()
print("   -> figures/surface_coverages.png")

# ============================================================
# Figure 5: Temperature dependence study
# ============================================================
print("\n[5] Running temperature study...")

T_study = np.linspace(450, 650, 30)
temp_results = temperature_study(T_study)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: TOF vs temperature
TOFs = [abs(r['rates'][0]) for r in temp_results]
axes[0, 0].plot(T_study, TOFs, 'b-o', markersize=3)
axes[0, 0].set_xlabel('Temperature [K]')
axes[0, 0].set_ylabel('CO consumption rate [mol/s]')
axes[0, 0].set_title('(a) Turnover frequency vs temperature')
axes[0, 0].set_yscale('log')
axes[0, 0].grid(True, alpha=0.3)

# Panel B: CO conversion in PFR
X_CO_pfr = [r['X_CO_pfr'] for r in temp_results]
X_CO_cstr = [r['X_CO_cstr'] for r in temp_results]
axes[0, 1].plot(T_study, X_CO_pfr, 'r-o', markersize=3, label='PFR')
axes[0, 1].plot(T_study, X_CO_cstr, 'b-s', markersize=3, label='CSTR')
axes[0, 1].set_xlabel('Temperature [K]')
axes[0, 1].set_ylabel('CO conversion X_CO')
axes[0, 1].set_title('(b) CO conversion: PFR vs CSTR')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

# Panel C: Coverage vs temperature
for i, label in enumerate(species_labels[:4]):
    coverages_T = [np.maximum(r['coverages'][i], 0) for r in temp_results]
    axes[1, 0].plot(T_study, coverages_T, '-o', markersize=3, label=label)

axes[1, 0].set_xlabel('Temperature [K]')
axes[1, 0].set_ylabel('Surface coverage θ')
axes[1, 0].set_title('(c) Key species coverages vs temperature')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)

# Panel D: Selectivity (CH4 vs C2H4)
CH4_rates = [abs(r['rates'][2]) for r in temp_results]
C2H4_rates = [abs(r['rates'][4]) for r in temp_results]
total_C = [ch4 + 2*c2 for ch4, c2 in zip(CH4_rates, C2H4_rates)]
S_CH4 = [ch4/tot if tot > 1e-30 else 0 for ch4, tot in zip(CH4_rates, total_C)]
S_C2 = [2*c2/tot if tot > 1e-30 else 0 for c2, tot in zip(C2H4_rates, total_C)]

axes[1, 1].plot(T_study, S_CH4, 'r-o', markersize=3, label='CH₄')
axes[1, 1].plot(T_study, S_C2, 'g-s', markersize=3, label='C₂H₄')
axes[1, 1].set_xlabel('Temperature [K]')
axes[1, 1].set_ylabel('Carbon selectivity')
axes[1, 1].set_title('(d) Product selectivity vs temperature')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/temperature_study.png')
plt.close()
print("   -> figures/temperature_study.png")

# ============================================================
# Figure 6: Degree of Rate Control analysis
# ============================================================
print("\n[6] Running DRC analysis...")

step_names, drc = degree_of_rate_control_analysis(T=500)

fig, ax = plt.subplots(figsize=(10, 5))

clean_names = ['CO ads.', 'H₂ diss.', 'CO diss.', 'C+H', 'CH+H',
               'CH₂+H', 'CH₃+H', 'O+H', 'OH+H', 'Chain growth']
colors_drc = ['red' if d > 0.3 else 'steelblue' for d in np.abs(drc)]
bars = ax.barh(range(len(drc)), drc, color=colors_drc, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(drc)))
ax.set_yticklabels(clean_names)
ax.set_xlabel('Degree of Rate Control (X_RC)')
ax.set_title('Degree of Rate Control at T=500K, P=20bar, H₂/CO=2')
ax.axvline(x=0, color='black', linewidth=0.5)
ax.grid(True, alpha=0.3, axis='x')

for i, (bar, val) in enumerate(zip(bars, drc)):
    ax.text(val + 0.02 if val >= 0 else val - 0.02, i,
            f'{val:.3f}', ha='left' if val >= 0 else 'right',
            va='center', fontsize=9)

plt.tight_layout()
plt.savefig('figures/degree_of_rate_control.png')
plt.close()
print("   -> figures/degree_of_rate_control.png")

# ============================================================
# Figure 7: PFR concentration profiles
# ============================================================
print("\n[7] Generating PFR profiles...")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

pfr_sol = results_500['pfr_solution']
F0 = results_500['F0']
gas_labels = ['CO', 'H₂', 'CH₄', 'H₂O', 'C₂H₄']
gas_colors = ['blue', 'red', 'green', 'cyan', 'orange']

if pfr_sol.success:
    for i, (label, color) in enumerate(zip(gas_labels, gas_colors)):
        axes[0].plot(pfr_sol.t, pfr_sol.y[i] * 1000, label=label, color=color)
    
    axes[0].set_xlabel('Catalyst weight [kg]')
    axes[0].set_ylabel('Molar flow rate [mmol/s]')
    axes[0].set_title('(a) PFR molar flow profiles')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Panel B: Conversion along PFR
    X_CO_profile = 1.0 - pfr_sol.y[0] / F0[0] if F0[0] > 0 else np.zeros_like(pfr_sol.t)
    axes[1].plot(pfr_sol.t, X_CO_profile * 100, 'b-', linewidth=2)
    axes[1].set_xlabel('Catalyst weight [kg]')
    axes[1].set_ylabel('CO conversion [%]')
    axes[1].set_title('(b) CO conversion along PFR')
    axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/pfr_profiles.png')
plt.close()
print("   -> figures/pfr_profiles.png")

# ============================================================
# Figure 8: Comparison with/without lateral interactions
# ============================================================
print("\n[8] Generating lateral interaction comparison...")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

T_comp = np.linspace(450, 650, 25)

# Without lateral interactions
rates_no_lat = []
for T in T_comp:
    res = run_ft_simulation(T=T)
    rates_no_lat.append(abs(res['rates'][0]))

# With lateral interactions (modify binding energies)
rates_with_lat = []
species_li, eps_li = create_ft_interaction_matrix()
li_model = LateralInteractionModel(species_li, eps_li, coordination_number=4)

for T in T_comp:
    # Approximate effect: increase Ea by mean-field correction
    theta_approx = np.array([0.3, 0.2, 0.05, 0.05, 0.02, 0.01, 0.005, 0.03])
    dE = li_model.mean_field_correction(theta_approx)
    
    k = compute_ft_rate_constants(T)
    # Increase CO dissociation barrier by lateral interaction correction
    Ea_corr = 1.60 + 0.5 * abs(dE[0])  # BEP correction
    k['k3f'] = compute_rate_constant(Ea_corr, T, 450, tunneling='wigner')
    
    mkm = SimpleMicroKineticModel(k)
    P_CO = 20.0 / 3.0
    P_H2 = 40.0 / 3.0
    P = np.array([P_CO, P_H2, 0, 0, 0])
    rates = mkm.compute_rates(P, T)
    rates_with_lat.append(abs(rates[0]))

axes[0].semilogy(T_comp, rates_no_lat, 'b-o', markersize=4, label='Without lateral int.')
axes[0].semilogy(T_comp, rates_with_lat, 'r-s', markersize=4, label='With lateral int.')
axes[0].set_xlabel('Temperature [K]')
axes[0].set_ylabel('CO consumption rate [mol/s]')
axes[0].set_title('(a) Effect of lateral interactions on rate')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Panel B: Apparent activation energy comparison
from numpy.polynomial.polynomial import polyfit

inv_T = 1000.0 / T_comp
ln_r_no = np.log(np.array(rates_no_lat) + 1e-50)
ln_r_with = np.log(np.array(rates_with_lat) + 1e-50)

mask_no = np.isfinite(ln_r_no)
mask_with = np.isfinite(ln_r_with)

if np.sum(mask_no) > 2:
    coeffs_no = np.polyfit(inv_T[mask_no], ln_r_no[mask_no], 1)
    Ea_app_no = -coeffs_no[0] * 8.314 / 1000  # kJ/mol
    axes[1].plot(inv_T[mask_no], ln_r_no[mask_no], 'bo', markersize=4)
    axes[1].plot(inv_T, np.polyval(coeffs_no, inv_T), 'b--',
                label=f'No lat. int. (Ea={Ea_app_no:.1f} kJ/mol)')

if np.sum(mask_with) > 2:
    coeffs_with = np.polyfit(inv_T[mask_with], ln_r_with[mask_with], 1)
    Ea_app_with = -coeffs_with[0] * 8.314 / 1000
    axes[1].plot(inv_T[mask_with], ln_r_with[mask_with], 'rs', markersize=4)
    axes[1].plot(inv_T, np.polyval(coeffs_with, inv_T), 'r--',
                label=f'With lat. int. (Ea={Ea_app_with:.1f} kJ/mol)')

axes[1].set_xlabel('1000/T [K⁻¹]')
axes[1].set_ylabel('ln(rate)')
axes[1].set_title('(b) Apparent activation energy')
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/lateral_interaction_effect.png')
plt.close()
print("   -> figures/lateral_interaction_effect.png")

# ============================================================
# Figure 9: Energy diagram
# ============================================================
print("\n[9] Generating energy diagram...")

fig, ax = plt.subplots(figsize=(14, 6))

# Simplified FT energy profile on Co(0001)
states = [
    'CO(g)+H₂(g)\n+*', 'CO*+H₂(g)', 'TS₁\n(CO diss.)', 'C*+O*+H₂',
    'C*+O*\n+2H*', 'TS₂\n(C+H)', 'CH*+O*+H*',
    'TS₃\n(CH+H)', 'CH₂*+O*', 'TS₄\n(chain gr.)',
    'C₂H₄(g)+O*', 'TS₅\n(O+H)', 'OH*', 'H₂O(g)+*'
]

energies = [0, -1.35, 0.25, -0.25, -1.15, -0.40, -1.15,
            -0.50, -1.10, -0.15, -0.80, 0.20, -0.90, -1.20]

x = np.arange(len(energies)) * 1.5
width = 0.8

for i in range(len(energies)):
    color = 'red' if 'TS' in states[i] else 'steelblue'
    ax.plot([x[i] - width/2, x[i] + width/2], [energies[i], energies[i]],
            color=color, linewidth=2.5)
    if i > 0:
        ax.plot([x[i-1] + width/2, x[i] - width/2],
                [energies[i-1], energies[i]],
                'k--', alpha=0.4, linewidth=0.8)
    
    va = 'bottom' if energies[i] >= (energies[i-1] if i > 0 else 0) else 'top'
    offset = 0.08 if va == 'bottom' else -0.08
    ax.text(x[i], energies[i] + offset, f'{energies[i]:.2f}',
            ha='center', va=va, fontsize=8, color='black')

ax.set_xticks(x)
ax.set_xticklabels(states, fontsize=7, rotation=30, ha='right')
ax.set_ylabel('Relative energy [eV]')
ax.set_title('Potential energy surface for Fischer-Tropsch on Co(0001)')
ax.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
ax.grid(True, alpha=0.2, axis='y')

plt.tight_layout()
plt.savefig('figures/energy_diagram.png')
plt.close()
print("   -> figures/energy_diagram.png")

# ============================================================
# Print summary
# ============================================================
print("\n" + "=" * 60)
print("SIMULATION SUMMARY")
print("=" * 60)

print(f"\nFischer-Tropsch on Co(0001) at T=500K, P=20bar, H2/CO=2")
print(f"\nSteady-state coverages:")
for i, label in enumerate(species_labels):
    print(f"  {label:>6s}: {theta_ss[i]:.4f}")
print(f"  {'*':>6s}: {max(1-np.sum(theta_ss), 0):.4f}")

print(f"\nGas-phase production rates [mol/s]:")
for i, label in enumerate(gas_labels):
    print(f"  {label:>6s}: {results_500['rates'][i]:.4e}")

print(f"\nCO conversion:")
print(f"  PFR:  {results_500['X_CO_pfr']*100:.2f}%")
print(f"  CSTR: {results_500['X_CO_cstr']*100:.2f}%")

print(f"\nDegree of Rate Control:")
for name, val in zip(clean_names, drc):
    marker = " <-- RDS" if abs(val) > 0.3 else ""
    print(f"  {name:>15s}: {val:.4f}{marker}")

print(f"\nGenerated figures:")
for f in sorted(os.listdir('figures')):
    print(f"  figures/{f}")

print("\nSimulation complete!")
