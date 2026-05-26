"""
Main experiment script: Model Selection Framework for Epidemic Models.
Generates all figures and numerical results.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import odeint
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

from epidemic_models import (
    sir_ode, seir_ode, seir_age_structured_ode, seir_spatial_ode,
    seir_vaccination_ode, ABMEpidemic,
    generate_synthetic_data, generate_covid_wave_data
)
from bayesian_inference import (
    estimate_sir_mle, estimate_seir_mle, model_selection_criteria,
    abc_rejection_seir, particle_filter_seir
)

plt.rcParams.update({'font.size': 11, 'figure.dpi': 150, 'savefig.bbox': 'tight'})
FIGDIR = 'figures'
os.makedirs(FIGDIR, exist_ok=True)
results = {}


# ============================================================
# Experiment 1: SIR vs SEIR Model Comparison
# ============================================================
print("=" * 60)
print("Experiment 1: SIR vs SEIR Model Comparison")
print("=" * 60)

N = 1e6
n_days = 150
t = np.arange(n_days)

# Generate data from SEIR (ground truth)
t_data, I_true, observed, new_cases_true = generate_synthetic_data('seir', N, n_days)

# Fit both models
sir_est = estimate_sir_mle(observed, t_data, N)
seir_est = estimate_seir_mle(observed, t_data, N)

print(f"SIR  estimates: beta={sir_est['beta']:.4f}, gamma={sir_est['gamma']:.4f}, R0={sir_est['R0']:.2f}")
print(f"SEIR estimates: beta={seir_est['beta']:.4f}, sigma={seir_est['sigma']:.4f}, gamma={seir_est['gamma']:.4f}, R0={seir_est['R0']:.2f}")

# Model selection
selection = model_selection_criteria(sir_est, seir_est, len(observed))
results['model_selection'] = selection
print(f"AIC preferred: {selection['preferred_AIC']} (ΔAIC={selection['delta_AIC']:.1f})")
print(f"BIC preferred: {selection['preferred_BIC']} (ΔBIC={selection['delta_BIC']:.1f})")

# Figure 1: SIR vs SEIR fit comparison
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Plot SIR fit
sol_sir = odeint(sir_ode, [N-100, 100, 0], t_data,
                 args=(sir_est['beta'], sir_est['gamma'], N))
new_sir = np.diff(np.concatenate([[0], np.cumsum(sol_sir[:, 1] * sir_est['gamma'])]))
new_sir = np.maximum(new_sir, 0)

# Plot SEIR fit
sol_seir = odeint(seir_ode, [N-100, 50, 50, 0], t_data,
                  args=(seir_est['beta'], seir_est['sigma'], seir_est['gamma'], N))
new_seir = np.diff(np.concatenate([[0], np.cumsum(sol_seir[:, 2] * seir_est['gamma'])]))
new_seir = np.maximum(new_seir, 0)

axes[0].plot(t_data, observed, 'k.', alpha=0.3, markersize=2, label='Observed')
axes[0].plot(t_data[:len(new_sir)], new_sir[:len(t_data)], 'b-', linewidth=2, label=f'SIR (R₀={sir_est["R0"]:.2f})')
axes[0].plot(t_data[:len(new_seir)], new_seir[:len(t_data)], 'r-', linewidth=2, label=f'SEIR (R₀={seir_est["R0"]:.2f})')
axes[0].set_xlabel('Days')
axes[0].set_ylabel('Daily New Cases')
axes[0].set_title('Model Fit Comparison')
axes[0].legend()
axes[0].set_xlim(0, n_days)

# Compartment dynamics (SEIR)
axes[1].plot(t_data, sol_seir[:, 0]/N, label='S', color='blue')
axes[1].plot(t_data, sol_seir[:, 1]/N, label='E', color='orange')
axes[1].plot(t_data, sol_seir[:, 2]/N, label='I', color='red')
axes[1].plot(t_data, sol_seir[:, 3]/N, label='R', color='green')
axes[1].set_xlabel('Days')
axes[1].set_ylabel('Proportion')
axes[1].set_title('SEIR Compartment Dynamics')
axes[1].legend()

# Model selection criteria
models = ['SIR', 'SEIR']
aic_vals = [selection['SIR']['AIC'], selection['SEIR']['AIC']]
bic_vals = [selection['SIR']['BIC'], selection['SEIR']['BIC']]
x = np.arange(len(models))
width = 0.35
bars1 = axes[2].bar(x - width/2, aic_vals, width, label='AIC', color='steelblue')
bars2 = axes[2].bar(x + width/2, bic_vals, width, label='BIC', color='coral')
axes[2].set_ylabel('Information Criterion')
axes[2].set_title('Model Selection (lower = better)')
axes[2].set_xticks(x)
axes[2].set_xticklabels(models)
axes[2].legend()

plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig1_model_comparison.png', dpi=150)
plt.close()
print("Saved: fig1_model_comparison.png")


# ============================================================
# Experiment 2: Age-Structured SEIR
# ============================================================
print("\n" + "=" * 60)
print("Experiment 2: Age-Structured SEIR Model")
print("=" * 60)

n_age_groups = 3
age_labels = ['0-19', '20-59', '60+']
N_groups = np.array([2e5, 5e5, 3e5])

# Contact matrix (from POLYMOD-like data)
contact_matrix = np.array([
    [8.0, 3.0, 1.5],
    [3.0, 6.0, 2.0],
    [1.5, 2.0, 4.0]
])
beta_matrix = 0.05 * contact_matrix

sigma, gamma = 0.2, 0.1
I0 = np.array([10, 30, 10])
y0 = np.concatenate([N_groups - I0*2, I0, I0, np.zeros(3)])

t_age = np.arange(200)
sol_age = odeint(seir_age_structured_ode, y0, t_age,
                 args=(beta_matrix, sigma, gamma, N_groups))

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
colors = ['#2196F3', '#FF9800', '#F44336']
for i, (label, color) in enumerate(zip(age_labels, colors)):
    axes[0].plot(t_age, sol_age[:, 2*3+i], color=color, linewidth=2, label=label)
axes[0].set_xlabel('Days')
axes[0].set_ylabel('Infected')
axes[0].set_title('Infected by Age Group')
axes[0].legend(title='Age Group')

# Attack rates
total_infected = sol_age[-1, 3*3:4*3]
attack_rates = total_infected / N_groups * 100
axes[1].bar(age_labels, attack_rates, color=colors)
axes[1].set_ylabel('Attack Rate (%)')
axes[1].set_title('Final Attack Rate by Age')
for i, v in enumerate(attack_rates):
    axes[1].text(i, v + 0.5, f'{v:.1f}%', ha='center', fontweight='bold')

# Peak timing
peak_days = np.argmax(sol_age[:, 2*3:3*3], axis=0)
axes[2].barh(age_labels, peak_days, color=colors)
axes[2].set_xlabel('Peak Day')
axes[2].set_title('Peak Infection Timing')
for i, v in enumerate(peak_days):
    axes[2].text(v + 1, i, f'Day {v}', va='center', fontweight='bold')

plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig2_age_structured.png', dpi=150)
plt.close()
print(f"Attack rates: {dict(zip(age_labels, [f'{ar:.1f}%' for ar in attack_rates]))}")
print(f"Peak days: {dict(zip(age_labels, peak_days))}")
results['age_structured'] = {
    'attack_rates': dict(zip(age_labels, attack_rates.tolist())),
    'peak_days': dict(zip(age_labels, peak_days.tolist()))
}
print("Saved: fig2_age_structured.png")


# ============================================================
# Experiment 3: Spatial SEIR Model
# ============================================================
print("\n" + "=" * 60)
print("Experiment 3: Spatial SEIR with Mobility")
print("=" * 60)

n_regions = 4
region_labels = ['Tokyo', 'Osaka', 'Hokkaido', 'Fukuoka']
N_regions = np.array([14e6, 8.8e6, 5.2e6, 5.1e6])

# Mobility matrix (normalized)
mobility = np.array([
    [0.0, 0.02, 0.005, 0.01],
    [0.02, 0.0, 0.003, 0.008],
    [0.005, 0.003, 0.0, 0.002],
    [0.01, 0.008, 0.002, 0.0]
])

beta_local = 0.4
beta_travel = 0.1
sigma_s, gamma_s = 0.2, 0.1

# Initial: infection starts in Tokyo
I0_spatial = np.array([1000, 0, 0, 0])
y0_spatial = np.concatenate([
    N_regions - I0_spatial*2,
    I0_spatial,
    I0_spatial,
    np.zeros(n_regions)
])

t_spatial = np.arange(180)
sol_spatial = odeint(seir_spatial_ode, y0_spatial, t_spatial,
                     args=(beta_local, beta_travel, sigma_s, gamma_s, N_regions, mobility))

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
region_colors = ['#E53935', '#1E88E5', '#43A047', '#FB8C00']
for i, (label, color) in enumerate(zip(region_labels, region_colors)):
    I_region = sol_spatial[:, 2*n_regions+i]
    axes[0].plot(t_spatial, I_region / N_regions[i] * 100, color=color, linewidth=2, label=label)
axes[0].set_xlabel('Days')
axes[0].set_ylabel('Infection Rate (%)')
axes[0].set_title('Spatial Spread Across Regions')
axes[0].legend()

# Arrival time (when I > threshold)
threshold = 100
arrival_times = []
for i in range(n_regions):
    I_region = sol_spatial[:, 2*n_regions+i]
    idx = np.where(I_region > threshold)[0]
    arrival_times.append(idx[0] if len(idx) > 0 else 180)

axes[1].barh(region_labels, arrival_times, color=region_colors)
axes[1].set_xlabel('Arrival Day (I > 100)')
axes[1].set_title('Epidemic Arrival Times')
for i, v in enumerate(arrival_times):
    axes[1].text(v + 1, i, f'Day {v}', va='center', fontweight='bold')

plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig3_spatial_spread.png', dpi=150)
plt.close()
print(f"Arrival times: {dict(zip(region_labels, arrival_times))}")
results['spatial'] = {'arrival_times': dict(zip(region_labels, arrival_times))}
print("Saved: fig3_spatial_spread.png")


# ============================================================
# Experiment 4: ABM vs ODE Comparison
# ============================================================
print("\n" + "=" * 60)
print("Experiment 4: Agent-Based Model vs ODE")
print("=" * 60)

n_agents = 2000
beta_abm = 0.15
sigma_abm = 0.2
gamma_abm = 0.1
n_steps = 150

# Run ABM multiple times
n_runs = 5
abm_histories = []
for run in range(n_runs):
    abm = ABMEpidemic(n_agents, beta_abm, sigma_abm, gamma_abm,
                      contact_radius=0.05, seed=42+run)
    abm.seed_infection(10)
    history = abm.run(n_steps)
    abm_histories.append(history)
    print(f"  ABM run {run+1}: peak I = {max(history['I'])}")

# Corresponding ODE
t_abm = np.arange(n_steps)
y0_ode = [n_agents - 10, 0, 10, 0]
R0_abm_approx = beta_abm * n_agents * np.pi * 0.05**2 / gamma_abm
beta_ode_equiv = R0_abm_approx * gamma_abm
sol_ode_abm = odeint(seir_ode, y0_ode, t_abm,
                     args=(beta_ode_equiv, sigma_abm, gamma_abm, n_agents))

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# ABM trajectories
for i, h in enumerate(abm_histories):
    alpha = 0.5 if i > 0 else 1.0
    axes[0].plot(h['I'], color='red', alpha=alpha, linewidth=1,
                label='ABM runs' if i == 0 else None)
axes[0].plot(sol_ode_abm[:, 2], 'b--', linewidth=2, label='ODE (mean-field)')
axes[0].set_xlabel('Days')
axes[0].set_ylabel('Infected')
axes[0].set_title('ABM vs ODE: Infected')
axes[0].legend()

# ABM mean vs ODE
abm_mean_I = np.mean([h['I'] for h in abm_histories], axis=0)
abm_std_I = np.std([h['I'] for h in abm_histories], axis=0)
axes[1].fill_between(range(n_steps), abm_mean_I - abm_std_I,
                     abm_mean_I + abm_std_I, alpha=0.3, color='red')
axes[1].plot(abm_mean_I, 'r-', linewidth=2, label='ABM mean ± SD')
axes[1].plot(sol_ode_abm[:, 2], 'b--', linewidth=2, label='ODE')
axes[1].set_xlabel('Days')
axes[1].set_ylabel('Infected')
axes[1].set_title('ABM Mean vs ODE')
axes[1].legend()

# Stochastic variability
peak_times = [np.argmax(h['I']) for h in abm_histories]
peak_sizes = [max(h['I']) for h in abm_histories]
axes[2].scatter(peak_times, peak_sizes, c='red', s=100, zorder=5, label='ABM peaks')
axes[2].axvline(np.argmax(sol_ode_abm[:, 2]), color='blue', linestyle='--', label='ODE peak time')
axes[2].axhline(np.max(sol_ode_abm[:, 2]), color='blue', linestyle=':', label='ODE peak size')
axes[2].set_xlabel('Peak Day')
axes[2].set_ylabel('Peak Infected')
axes[2].set_title('Stochastic Variability (ABM)')
axes[2].legend()

plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig4_abm_vs_ode.png', dpi=150)
plt.close()
results['abm_vs_ode'] = {
    'abm_peak_mean': float(np.mean(peak_sizes)),
    'abm_peak_std': float(np.std(peak_sizes)),
    'ode_peak': float(np.max(sol_ode_abm[:, 2])),
    'abm_peak_day_mean': float(np.mean(peak_times)),
    'ode_peak_day': int(np.argmax(sol_ode_abm[:, 2]))
}
print(f"ABM peak: {np.mean(peak_sizes):.0f} ± {np.std(peak_sizes):.0f} (day {np.mean(peak_times):.0f})")
print(f"ODE peak: {np.max(sol_ode_abm[:, 2]):.0f} (day {np.argmax(sol_ode_abm[:, 2])})")
print("Saved: fig4_abm_vs_ode.png")


# ============================================================
# Experiment 5: Parameter Estimation (ABC + Particle Filter)
# ============================================================
print("\n" + "=" * 60)
print("Experiment 5: Parameter Estimation")
print("=" * 60)

# ABC
print("Running ABC rejection sampler...")
abc_particles = abc_rejection_seir(observed, t_data, N, n_particles=200, epsilon=0.15, seed=42)
print(f"  Accepted {len(abc_particles)} particles")

if len(abc_particles) > 0:
    abc_betas = [p['beta'] for p in abc_particles]
    abc_sigmas = [p['sigma'] for p in abc_particles]
    abc_gammas = [p['gamma'] for p in abc_particles]
    abc_R0s = [p['R0'] for p in abc_particles]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    axes[0, 0].hist(abc_betas, bins=30, color='steelblue', alpha=0.7, density=True)
    axes[0, 0].axvline(0.35, color='red', linestyle='--', label='True β=0.35')
    axes[0, 0].set_xlabel('β')
    axes[0, 0].set_title('ABC Posterior: β')
    axes[0, 0].legend()

    axes[0, 1].hist(abc_sigmas, bins=30, color='orange', alpha=0.7, density=True)
    axes[0, 1].axvline(0.2, color='red', linestyle='--', label='True σ=0.2')
    axes[0, 1].set_xlabel('σ')
    axes[0, 1].set_title('ABC Posterior: σ')
    axes[0, 1].legend()

    axes[1, 0].hist(abc_gammas, bins=30, color='green', alpha=0.7, density=True)
    axes[1, 0].axvline(0.1, color='red', linestyle='--', label='True γ=0.1')
    axes[1, 0].set_xlabel('γ')
    axes[1, 0].set_title('ABC Posterior: γ')
    axes[1, 0].legend()

    axes[1, 1].hist(abc_R0s, bins=30, color='purple', alpha=0.7, density=True)
    axes[1, 1].axvline(3.5, color='red', linestyle='--', label='True R₀=3.5')
    axes[1, 1].set_xlabel('R₀')
    axes[1, 1].set_title('ABC Posterior: R₀')
    axes[1, 1].legend()

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig5_abc_posteriors.png', dpi=150)
    plt.close()
    print(f"  ABC R₀ estimate: {np.mean(abc_R0s):.2f} ± {np.std(abc_R0s):.2f}")
    results['abc'] = {
        'beta': f'{np.mean(abc_betas):.4f} ± {np.std(abc_betas):.4f}',
        'sigma': f'{np.mean(abc_sigmas):.4f} ± {np.std(abc_sigmas):.4f}',
        'gamma': f'{np.mean(abc_gammas):.4f} ± {np.std(abc_gammas):.4f}',
        'R0': f'{np.mean(abc_R0s):.2f} ± {np.std(abc_R0s):.2f}',
        'n_accepted': len(abc_particles)
    }
    print("Saved: fig5_abc_posteriors.png")

# Particle Filter
print("Running Particle Filter...")
filtered_states, ess = particle_filter_seir(observed, N, n_particles=300, seed=42)

fig, axes = plt.subplots(1, 3, figsize=(16, 5))
axes[0].plot(t_data, I_true, 'r-', linewidth=2, label='True I(t)')
axes[0].plot(t_data, filtered_states[:, 2], 'b--', linewidth=2, label='Filtered I(t)')
axes[0].fill_between(t_data, filtered_states[:, 2]*0.8, filtered_states[:, 2]*1.2,
                     alpha=0.2, color='blue')
axes[0].set_xlabel('Days')
axes[0].set_ylabel('Infected')
axes[0].set_title('Particle Filter: State Estimation')
axes[0].legend()

axes[1].plot(t_data, ess, 'g-', linewidth=1.5)
axes[1].axhline(300 * 0.5, color='red', linestyle='--', label='ESS threshold (50%)')
axes[1].set_xlabel('Days')
axes[1].set_ylabel('ESS')
axes[1].set_title('Effective Sample Size')
axes[1].legend()

# Residual
residual = I_true - filtered_states[:, 2]
axes[2].plot(t_data, residual, 'k-', linewidth=1)
axes[2].axhline(0, color='gray', linestyle='--')
axes[2].set_xlabel('Days')
axes[2].set_ylabel('Residual')
axes[2].set_title('Filter Residual (True - Filtered)')

plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig6_particle_filter.png', dpi=150)
plt.close()
print(f"  Mean ESS: {np.mean(ess):.1f}")
results['particle_filter'] = {
    'mean_ess': float(np.mean(ess)),
    'rmse': float(np.sqrt(np.mean(residual**2)))
}
print("Saved: fig6_particle_filter.png")


# ============================================================
# Experiment 6: Intervention Scenarios
# ============================================================
print("\n" + "=" * 60)
print("Experiment 6: Intervention Scenario Analysis")
print("=" * 60)

N_int = 1e6
t_int = np.arange(300)

scenarios = {
    'No Intervention': {'beta': 0.35, 'vax_rate': 0, 'vax_eff': 0},
    'Behavioral (β×0.6)': {'beta': 0.21, 'vax_rate': 0, 'vax_eff': 0},
    'Vaccine (70% eff)': {'beta': 0.35, 'vax_rate': 0.005, 'vax_eff': 0.7},
    'Combined': {'beta': 0.21, 'vax_rate': 0.005, 'vax_eff': 0.7},
}

scenario_colors = ['#D32F2F', '#FF9800', '#4CAF50', '#1976D2']
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

scenario_results = {}
for idx, (name, params) in enumerate(scenarios.items()):
    if params['vax_rate'] > 0:
        y0_v = [N_int - 200, 100, 100, 0, 0]
        sol = odeint(seir_vaccination_ode, y0_v, t_int,
                     args=(params['beta'], 0.2, 0.1, N_int,
                           params['vax_rate'], params['vax_eff']))
        I_scenario = sol[:, 2]
        R_scenario = sol[:, 3]
    else:
        y0_s = [N_int - 200, 100, 100, 0]
        sol = odeint(seir_ode, y0_s, t_int,
                     args=(params['beta'], 0.2, 0.1, N_int))
        I_scenario = sol[:, 2]
        R_scenario = sol[:, 3]

    peak_I = np.max(I_scenario)
    peak_day = np.argmax(I_scenario)
    total_infected = R_scenario[-1]

    scenario_results[name] = {
        'peak_I': float(peak_I),
        'peak_day': int(peak_day),
        'total_infected': float(total_infected),
        'attack_rate': float(total_infected / N_int * 100)
    }

    ax = axes[idx // 2, idx % 2]
    color = scenario_colors[idx]
    if params['vax_rate'] > 0:
        ax.plot(t_int, sol[:, 0]/N_int, label='S', color='blue')
        ax.plot(t_int, sol[:, 1]/N_int, label='E', color='orange')
        ax.plot(t_int, sol[:, 2]/N_int, label='I', color='red')
        ax.plot(t_int, sol[:, 3]/N_int, label='R', color='green')
        ax.plot(t_int, sol[:, 4]/N_int, label='V', color='purple')
    else:
        ax.plot(t_int, sol[:, 0]/N_int, label='S', color='blue')
        ax.plot(t_int, sol[:, 1]/N_int, label='E', color='orange')
        ax.plot(t_int, sol[:, 2]/N_int, label='I', color='red')
        ax.plot(t_int, sol[:, 3]/N_int, label='R', color='green')
    ax.set_title(f'{name}\nPeak: {peak_I/1000:.0f}K (day {peak_day}), AR: {total_infected/N_int*100:.1f}%')
    ax.set_xlabel('Days')
    ax.set_ylabel('Proportion')
    ax.legend(loc='right')

plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig7_interventions.png', dpi=150)
plt.close()

for name, res in scenario_results.items():
    print(f"  {name}: peak={res['peak_I']/1000:.0f}K (day {res['peak_day']}), AR={res['attack_rate']:.1f}%")
results['interventions'] = scenario_results
print("Saved: fig7_interventions.png")

# Intervention comparison summary figure
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
names = list(scenario_results.keys())
peaks = [scenario_results[n]['peak_I']/1000 for n in names]
attack_rates = [scenario_results[n]['attack_rate'] for n in names]
peak_days = [scenario_results[n]['peak_day'] for n in names]

axes[0].barh(names, peaks, color=scenario_colors)
axes[0].set_xlabel('Peak Infected (thousands)')
axes[0].set_title('Peak Infected')

axes[1].barh(names, attack_rates, color=scenario_colors)
axes[1].set_xlabel('Attack Rate (%)')
axes[1].set_title('Final Attack Rate')

axes[2].barh(names, peak_days, color=scenario_colors)
axes[2].set_xlabel('Peak Day')
axes[2].set_title('Peak Timing')

plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig8_intervention_comparison.png', dpi=150)
plt.close()
print("Saved: fig8_intervention_comparison.png")


# ============================================================
# Experiment 7: COVID-19 Wave Case Study (6th/7th)
# ============================================================
print("\n" + "=" * 60)
print("Experiment 7: COVID-19 Wave Case Study")
print("=" * 60)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

wave_results = {}
for w_idx, wave in enumerate(['6th', '7th']):
    t_w, I_true_w, obs_w, new_cases_w = generate_covid_wave_data(wave)
    N_japan = 1.26e8

    # Fit SEIR
    seir_w = estimate_seir_mle(obs_w, t_w, N_japan)
    sol_w = odeint(seir_ode, [N_japan - 10000, 5000, 5000, 0], t_w,
                   args=(seir_w['beta'], seir_w['sigma'], seir_w['gamma'], N_japan))
    new_w_fit = np.diff(np.concatenate([[0], np.cumsum(sol_w[:, 2] * seir_w['gamma'])]))
    new_w_fit = np.maximum(new_w_fit, 0) * 0.3  # reporting rate

    wave_results[wave] = {
        'beta': seir_w['beta'],
        'sigma': seir_w['sigma'],
        'gamma': seir_w['gamma'],
        'R0': seir_w['R0']
    }

    row = w_idx
    axes[row, 0].bar(t_w, obs_w, color='gray', alpha=0.5, label='Reported cases')
    axes[row, 0].plot(t_w[:len(new_w_fit)], new_w_fit[:len(t_w)], 'r-', linewidth=2,
                     label=f'SEIR fit (R₀={seir_w["R0"]:.2f})')
    axes[row, 0].set_xlabel('Days')
    axes[row, 0].set_ylabel('Daily Cases')
    axes[row, 0].set_title(f'COVID-19 {wave} Wave (Japan)')
    axes[row, 0].legend()

    axes[row, 1].plot(t_w, sol_w[:, 0]/N_japan, label='S', color='blue')
    axes[row, 1].plot(t_w, sol_w[:, 1]/N_japan, label='E', color='orange')
    axes[row, 1].plot(t_w, sol_w[:, 2]/N_japan, label='I', color='red')
    axes[row, 1].plot(t_w, sol_w[:, 3]/N_japan, label='R', color='green')
    axes[row, 1].set_xlabel('Days')
    axes[row, 1].set_ylabel('Proportion')
    axes[row, 1].set_title(f'{wave} Wave: Compartments')
    axes[row, 1].legend()

    print(f"  {wave} wave: β={seir_w['beta']:.4f}, σ={seir_w['sigma']:.4f}, γ={seir_w['gamma']:.4f}, R₀={seir_w['R0']:.2f}")

results['covid_waves'] = wave_results
plt.tight_layout()
plt.savefig(f'{FIGDIR}/fig9_covid_waves.png', dpi=150)
plt.close()
print("Saved: fig9_covid_waves.png")


# ============================================================
# Experiment 8: Model Selection Decision Framework
# ============================================================
print("\n" + "=" * 60)
print("Experiment 8: Model Selection Decision Framework")
print("=" * 60)

fig, ax = plt.subplots(1, 1, figsize=(14, 8))
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.set_aspect('equal')
ax.axis('off')

# Decision tree
boxes = [
    (5, 9.2, 'Research Question', '#E3F2FD', 'bold'),
    (2.5, 7.5, 'Individual\nheterogeneity\nneeded?', '#FFF3E0', 'normal'),
    (7.5, 7.5, 'Population-level\ndynamics\nsufficient?', '#E8F5E9', 'normal'),
    (1, 5.5, 'Agent-Based\nModel (ABM)', '#FFCDD2', 'bold'),
    (4, 5.5, 'Network/\nSpatial?', '#FFF9C4', 'normal'),
    (6, 5.5, 'Age\nstructure?', '#FFF9C4', 'normal'),
    (8.5, 5.5, 'Basic SIR/\nSEIR', '#C8E6C9', 'bold'),
    (3, 3.5, 'Spatial SEIR\n+ Mobility', '#B3E5FC', 'bold'),
    (5.5, 3.5, 'Age-structured\nSEIR', '#B3E5FC', 'bold'),
    (5, 1.5, 'Parameter Estimation\n(MCMC / ABC / PF)', '#E1BEE7', 'bold'),
    (5, 0.3, 'Model Selection\n(AIC/BIC/WAIC/LOO-CV)', '#F3E5F5', 'bold'),
]

for x, y, text, color, weight in boxes:
    bbox = dict(boxstyle='round,pad=0.3', facecolor=color, edgecolor='gray')
    ax.text(x, y, text, ha='center', va='center', fontsize=9,
            fontweight=weight, bbox=bbox)

# Arrows
arrows = [
    (5, 8.8, 2.5, 8.1), (5, 8.8, 7.5, 8.1),
    (2.5, 6.9, 1, 6.1), (2.5, 6.9, 4, 6.1),
    (7.5, 6.9, 6, 6.1), (7.5, 6.9, 8.5, 6.1),
    (4, 5.0, 3, 4.1), (6, 5.0, 5.5, 4.1),
    (1, 5.0, 5, 2.0), (3, 3.0, 5, 2.0),
    (5.5, 3.0, 5, 2.0), (8.5, 5.0, 5, 2.0),
    (5, 1.0, 5, 0.7),
]

for x1, y1, x2, y2 in arrows:
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
               arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))

ax.set_title('Model Structure Selection Framework', fontsize=14, fontweight='bold', pad=20)
plt.savefig(f'{FIGDIR}/fig10_decision_framework.png', dpi=150)
plt.close()
print("Saved: fig10_decision_framework.png")


# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
print("ALL EXPERIMENTS COMPLETE")
print("=" * 60)
print(f"\nGenerated figures in {FIGDIR}/:")
for f in sorted(os.listdir(FIGDIR)):
    if f.endswith('.png'):
        size = os.path.getsize(os.path.join(FIGDIR, f)) / 1024
        print(f"  {f} ({size:.0f} KB)")

import json
with open('experiment_results.json', 'w') as f:
    json.dump(results, f, indent=2, default=str)
print("\nResults saved to experiment_results.json")
