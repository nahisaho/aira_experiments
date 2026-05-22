"""
Digital Twin Main Pipeline - Run all modules and generate figures
"""

import sys
import os
import json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
import seaborn as sns

sns.set_theme(style='whitegrid', font_scale=1.1)
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'

from src.flow_simulation import (ResinProperties, MoldGeometry, ProcessConditions,
                                  HeleShawSolver, run_3d_flow_analysis)
from src.cooling_crystallization import (CrystallizationParams, CoolingGeometry,
                                          CoolingSimulator)
from src.residual_stress import run_stress_warpage_analysis
from src.process_quality_model import ProcessQualityModel
from src.data_assimilation import EnsembleKalmanFilter
from src.case_study import AutomotiveCaseStudy


def save_json(data, filepath):
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2, default=str)


# ============================================================
# Module 1: Flow Simulation
# ============================================================
print("=" * 60)
print("Module 1: Resin Flow Simulation")
print("=" * 60)

resin = ResinProperties()
geom = MoldGeometry()
cond = ProcessConditions()

solver = HeleShawSolver(geom, resin, cond)
fill_history = solver.simulate_filling(n_steps=20)
flow_results = solver.get_results()

results_3d = run_3d_flow_analysis(geom, resin, cond)

print(f"  Hele-Shaw: Max Pressure = {flow_results['max_pressure_MPa']:.1f} MPa")
print(f"  Hele-Shaw: Avg Melt Temp = {flow_results['avg_melt_temp_C']:.1f} °C")
print(f"  3D Flow: Max Velocity = {results_3d['max_velocity_m_s']:.4f} m/s")
print(f"  3D Flow: Temp Range = {results_3d['temperature_range_C'][0]:.1f} - {results_3d['temperature_range_C'][1]:.1f} °C")

# Figure 1: Flow simulation results
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Module 1: Resin Flow Simulation Results', fontsize=14, fontweight='bold')

ax = axes[0, 0]
ax.plot(fill_history['time'], fill_history['fill_percent'], 'b-o', markersize=4)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Fill Percentage (%)')
ax.set_title('Cavity Filling Progress')
ax.grid(True, alpha=0.3)

ax = axes[0, 1]
ax.plot(fill_history['time'], fill_history['max_pressure'], 'r-s', markersize=4)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Max Pressure (MPa)')
ax.set_title('Injection Pressure History')
ax.grid(True, alpha=0.3)

ax = axes[1, 0]
ax.plot(fill_history['time'], fill_history['avg_temperature'], 'g-^', markersize=4)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Average Temperature (°C)')
ax.set_title('Melt Temperature During Filling')
ax.grid(True, alpha=0.3)

ax = axes[1, 1]
nz_3d = len(results_3d['velocity_profile'])
z_norm = np.linspace(-1, 1, nz_3d)
ax.plot(results_3d['velocity_profile'], z_norm, 'b-', linewidth=2, label='Velocity')
ax.set_xlabel('Velocity (m/s)')
ax.set_ylabel('Normalized Thickness (z/h)')
ax.set_title('3D: Through-Thickness Velocity Profile')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/fig01_flow_simulation.png')
plt.savefig('figures/fig01_flow_simulation.svg')
plt.close()
print("  -> Saved figures/fig01_flow_simulation.png")

# Figure 2: Pressure and temperature fields
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle('Module 1: Pressure & Temperature Fields', fontsize=14, fontweight='bold')

P_field = np.array(flow_results['pressure_field']) / 1e6
im1 = axes[0].imshow(P_field.T, aspect='auto', cmap='hot', origin='lower',
                      extent=[0, 200, 0, 100])
axes[0].set_xlabel('Length (mm)')
axes[0].set_ylabel('Width (mm)')
axes[0].set_title('Pressure Field (MPa)')
plt.colorbar(im1, ax=axes[0], label='Pressure (MPa)')

T_field = np.array(flow_results['temperature_field']) - 273.15
im2 = axes[1].imshow(T_field.T, aspect='auto', cmap='coolwarm', origin='lower',
                      extent=[0, 200, 0, 100])
axes[1].set_xlabel('Length (mm)')
axes[1].set_ylabel('Width (mm)')
axes[1].set_title('Temperature Field (°C)')
plt.colorbar(im2, ax=axes[1], label='Temperature (°C)')

plt.tight_layout()
plt.savefig('figures/fig02_pressure_temperature_fields.png')
plt.savefig('figures/fig02_pressure_temperature_fields.svg')
plt.close()
print("  -> Saved figures/fig02_pressure_temperature_fields.png")

save_json({'hele_shaw': flow_results, '3d_flow': results_3d, 'fill_history': fill_history},
          'results/flow_simulation_results.json')

# ============================================================
# Module 2: Cooling & Crystallization
# ============================================================
print("\n" + "=" * 60)
print("Module 2: Cooling & Crystallization")
print("=" * 60)

cool_geom = CoolingGeometry()
cryst = CrystallizationParams()
cooling_sim = CoolingSimulator(cool_geom, cryst)
cooling_history = cooling_sim.simulate(total_time=25.0, dt=0.005)
cooling_summary = cooling_sim.get_summary()

for k, v in cooling_summary.items():
    print(f"  {k}: {v:.2f}")
if cooling_history['ejection_time']:
    print(f"  Ejection time: {cooling_history['ejection_time']:.1f} s")

# Figure 3: Cooling curves
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Module 2: Cooling & Crystallization Analysis', fontsize=14, fontweight='bold')

ax = axes[0, 0]
ax.plot(cooling_history['time'], cooling_history['center_temp_C'], 'r-', label='Center', linewidth=2)
ax.plot(cooling_history['time'], cooling_history['surface_temp_C'], 'b-', label='Surface', linewidth=2)
ax.axhline(y=80, color='g', linestyle='--', alpha=0.5, label='Ejection Temp')
if cooling_history['ejection_time']:
    ax.axvline(x=cooling_history['ejection_time'], color='orange', linestyle=':', alpha=0.7)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Temperature (°C)')
ax.set_title('Cooling Curves')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[0, 1]
ax.plot(cooling_history['time'], cooling_history['avg_crystallinity'], 'purple', label='Average', linewidth=2)
ax.plot(cooling_history['time'], cooling_history['center_crystallinity'], 'orange', label='Center', linewidth=2)
ax.set_xlabel('Time (s)')
ax.set_ylabel('Crystallinity')
ax.set_title('Crystallization Kinetics')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1, 0]
ax.plot(cooling_history['final_z_mm'], cooling_history['final_temperature_C'], 'r-', linewidth=2)
ax.set_xlabel('Through-Thickness Position (mm)')
ax.set_ylabel('Temperature (°C)')
ax.set_title('Final Temperature Profile')
ax.grid(True, alpha=0.3)

ax = axes[1, 1]
ax.plot(cooling_history['final_z_mm'], cooling_history['final_crystallinity'], 'purple', linewidth=2)
ax.set_xlabel('Through-Thickness Position (mm)')
ax.set_ylabel('Crystallinity')
ax.set_title('Final Crystallinity Profile')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/fig03_cooling_crystallization.png')
plt.savefig('figures/fig03_cooling_crystallization.svg')
plt.close()
print("  -> Saved figures/fig03_cooling_crystallization.png")

save_json({'summary': cooling_summary, 'ejection_time': cooling_history['ejection_time']},
          'results/cooling_results.json')

# ============================================================
# Module 3: Residual Stress & Warpage
# ============================================================
print("\n" + "=" * 60)
print("Module 3: Residual Stress & Warpage")
print("=" * 60)

T_profile = np.array(cooling_history['final_temperature_C'])
X_profile = np.array(cooling_history['final_crystallinity'])
stress_warpage = run_stress_warpage_analysis(T_profile, X_profile)

print(f"  Max Tensile Stress: {stress_warpage['stress']['max_tensile_MPa']:.1f} MPa")
print(f"  Max Compressive Stress: {stress_warpage['stress']['max_compressive_MPa']:.1f} MPa")
print(f"  Curvature: {stress_warpage['stress']['curvature_1_m']:.4f} 1/m")
print(f"  Warpage X: {stress_warpage['warpage']['warpage_x_mm']:.3f} mm")
print(f"  Warpage Y: {stress_warpage['warpage']['warpage_y_mm']:.3f} mm")
print(f"  Total Warpage: {stress_warpage['warpage']['warpage_total_mm']:.3f} mm")

# Figure 4: Stress and warpage
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
fig.suptitle('Module 3: Residual Stress & Warpage Prediction', fontsize=14, fontweight='bold')

ax = axes[0]
ax.plot(stress_warpage['stress']['z_mm'], stress_warpage['stress']['stress_MPa'], 'r-', linewidth=2)
ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
ax.fill_between(stress_warpage['stress']['z_mm'], stress_warpage['stress']['stress_MPa'],
                alpha=0.3, color='red')
ax.set_xlabel('Through-Thickness Position (mm)')
ax.set_ylabel('Residual Stress (MPa)')
ax.set_title('Through-Thickness Stress Profile')
ax.grid(True, alpha=0.3)

ax = axes[1]
ax.plot(stress_warpage['stress']['z_mm'], stress_warpage['stress']['modulus_GPa'], 'b-', linewidth=2)
ax.set_xlabel('Through-Thickness Position (mm)')
ax.set_ylabel('Elastic Modulus (GPa)')
ax.set_title('Modulus Distribution')
ax.grid(True, alpha=0.3)

ax = axes[2]
W = np.array(stress_warpage['warpage']['deformation_field_mm'])
x_coords = stress_warpage['warpage']['x_coords_mm']
y_coords = stress_warpage['warpage']['y_coords_mm']
im = ax.imshow(W, aspect='auto', cmap='RdBu_r', origin='lower',
               extent=[min(x_coords), max(x_coords), min(y_coords), max(y_coords)])
ax.set_xlabel('Length (mm)')
ax.set_ylabel('Width (mm)')
ax.set_title('Warpage Deformation (mm)')
plt.colorbar(im, ax=ax, label='Deformation (mm)')

plt.tight_layout()
plt.savefig('figures/fig04_stress_warpage.png')
plt.savefig('figures/fig04_stress_warpage.svg')
plt.close()
print("  -> Saved figures/fig04_stress_warpage.png")

save_json(stress_warpage, 'results/stress_warpage_results.json')

# ============================================================
# Module 4: Process-Quality Surrogate Model
# ============================================================
print("\n" + "=" * 60)
print("Module 4: Process-Quality Surrogate Model")
print("=" * 60)

pq_model = ProcessQualityModel()
training_result = pq_model.train(n_samples=80)

for name, scores in training_result['cv_scores'].items():
    print(f"  {name}: R² = {scores['r2_mean']:.3f} ± {scores['r2_std']:.3f}")

# Sensitivity analysis
print("\n  Computing Sobol sensitivity indices...")
sensitivities = pq_model.sobol_sensitivity(n_samples=100)

# Optimization
print("  Running process optimization...")
opt_result = pq_model.optimize()
print(f"  Optimization success: {opt_result['optimization_success']}")
print("  Optimal parameters:")
for k, v in opt_result['optimal_parameters'].items():
    print(f"    {k}: {v:.1f}")

# Figure 5: Sensitivity heatmap
fig, ax = plt.subplots(figsize=(12, 6))
fig.suptitle('Module 4: Process Parameter Sensitivity Analysis (Sobol Indices)',
             fontsize=14, fontweight='bold')

quality_labels = list(sensitivities.keys())
param_labels = list(list(sensitivities.values())[0].keys())
S_matrix = np.zeros((len(quality_labels), len(param_labels)))
for i, q in enumerate(quality_labels):
    for j, p in enumerate(param_labels):
        S_matrix[i, j] = sensitivities[q][p]

short_params = ['P_inj', 'P_pack', 't_cool', 'T_melt', 'T_mold', 'v_inj']
short_quality = ['Warpage', 'Sink', 'Weight', 'Shrinkage', 'Stress']

sns.heatmap(S_matrix, annot=True, fmt='.2f', cmap='YlOrRd',
            xticklabels=short_params, yticklabels=short_quality,
            ax=ax, vmin=0, vmax=0.5, cbar_kws={'label': 'Sobol Index (S1)'})
ax.set_xlabel('Process Parameters')
ax.set_ylabel('Quality Metrics')

plt.tight_layout()
plt.savefig('figures/fig05_sensitivity_analysis.png')
plt.savefig('figures/fig05_sensitivity_analysis.svg')
plt.close()
print("  -> Saved figures/fig05_sensitivity_analysis.png")

# Figure 6: Surrogate model predictions
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle('Module 4: Surrogate Model – Predicted vs Actual',
             fontsize=14, fontweight='bold')

for idx, name in enumerate(pq_model.quality_names):
    ax = axes[idx // 3, idx % 3]
    y_actual = pq_model.y_train[name]
    y_pred = []
    for i in range(len(pq_model.X_train)):
        pred = pq_model.predict(pq_model.X_train[i])
        y_pred.append(pred[name]['mean'])
    y_pred = np.array(y_pred)

    ax.scatter(y_actual, y_pred, alpha=0.5, s=20)
    lims = [min(y_actual.min(), y_pred.min()), max(y_actual.max(), y_pred.max())]
    ax.plot(lims, lims, 'r--', linewidth=1)
    ax.set_xlabel('Actual')
    ax.set_ylabel('Predicted')
    ax.set_title(name.replace('_', ' ').title())
    ax.grid(True, alpha=0.3)

axes[1, 2].set_visible(False)
plt.tight_layout()
plt.savefig('figures/fig06_surrogate_model.png')
plt.savefig('figures/fig06_surrogate_model.svg')
plt.close()
print("  -> Saved figures/fig06_surrogate_model.png")

save_json({'training': training_result, 'sensitivities': sensitivities,
           'optimization': opt_result}, 'results/process_quality_results.json')

# ============================================================
# Module 5: Data Assimilation (EnKF)
# ============================================================
print("\n" + "=" * 60)
print("Module 5: Data Assimilation (EnKF)")
print("=" * 60)

np.random.seed(42)
enkf = EnsembleKalmanFilter(n_ensemble=50)
da_results = enkf.run_assimilation(n_cycles=25)

print(f"  Convergence: {da_results['convergence_achieved']}")
print(f"  Final RMSE: {da_results['final_rmse']:.3f}")
print("  Parameter estimation errors:")
for name in enkf.state_names[:4]:
    err = da_results['estimation_error'][name]
    unc = da_results['uncertainty'][name]
    print(f"    {name}: error={err:.3f}, uncertainty={unc:.3f}")

# Figure 7: Data assimilation convergence
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Module 5: Ensemble Kalman Filter – Data Assimilation',
             fontsize=14, fontweight='bold')

ax = axes[0, 0]
ax.plot(da_results['rmse_history'], 'b-o', markersize=4)
ax.set_xlabel('Assimilation Cycle')
ax.set_ylabel('RMSE')
ax.set_title('Observation RMSE Convergence')
ax.grid(True, alpha=0.3)

state_history = np.array(da_results['state_mean_history'])
std_history = np.array(da_results['state_std_history'])
true_state = np.array([da_results['true_state'][n] for n in enkf.state_names])
cycles = np.arange(len(state_history))

for plot_idx, param_idx in enumerate([0, 2, 4]):
    ax = axes[(plot_idx + 1) // 2, (plot_idx + 1) % 2]
    name = enkf.state_names[param_idx]
    mean_vals = state_history[:, param_idx]
    std_vals = std_history[:, param_idx]

    ax.plot(cycles, mean_vals, 'b-', linewidth=2, label='Estimated')
    ax.fill_between(cycles, mean_vals - 2 * std_vals, mean_vals + 2 * std_vals,
                    alpha=0.2, color='blue', label='95% CI')
    ax.axhline(y=true_state[param_idx], color='r', linestyle='--', linewidth=1.5,
               label=f'True ({true_state[param_idx]:.2f})')
    ax.set_xlabel('Assimilation Cycle')
    ax.set_ylabel('Value')
    ax.set_title(name.replace('_', ' ').title())
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/fig07_data_assimilation.png')
plt.savefig('figures/fig07_data_assimilation.svg')
plt.close()
print("  -> Saved figures/fig07_data_assimilation.png")

save_json(da_results, 'results/data_assimilation_results.json')

# ============================================================
# Module 6: Automotive Case Study
# ============================================================
print("\n" + "=" * 60)
print("Module 6: Automotive Case Study")
print("=" * 60)

case = AutomotiveCaseStudy()
scenario_results = case.run_all_scenarios()

for name, result in scenario_results.items():
    status = "✓ PASS" if result['quality_pass'] else "✗ FAIL"
    print(f"  [{status}] {name}:")
    print(f"    Warpage: {result['quality_metrics']['warpage_mm']:.4f} mm")
    print(f"    Sink: {result['quality_metrics']['sink_depth_mm']:.4f} mm")
    print(f"    Weight: {result['quality_metrics']['weight_g']:.2f} g")
    print(f"    Cycle Time: {result['process_metrics']['cycle_time_s']:.1f} s")

# Monte Carlo simulation
print("\n  Running Monte Carlo analysis (500 samples)...")
mc_results = case.run_monte_carlo('nominal', n_samples=500)
print(f"  Pass Rate: {mc_results['pass_rate_pct']:.1f}%")

# Figure 8: Scenario comparison
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Module 6: Automotive Case Study – Scenario Comparison',
             fontsize=14, fontweight='bold')

scenarios = list(scenario_results.keys())
colors_map = {'nominal': '#2196F3', 'high_speed': '#F44336', 'low_stress': '#4CAF50', 'optimized': '#FF9800'}

metrics_to_plot = ['warpage_mm', 'sink_depth_mm', 'weight_g', 'shrinkage_pct']
limits = [0.5, 0.02, None, 1.0]
titles = ['Warpage (mm)', 'Sink Depth (mm)', 'Part Weight (g)', 'Shrinkage (%)']

for idx, (metric, limit, title) in enumerate(zip(metrics_to_plot, limits, titles)):
    ax = axes[idx // 2, idx % 2]
    values = [scenario_results[s]['quality_metrics'][metric] for s in scenarios]
    colors = [colors_map[s] for s in scenarios]
    bars = ax.bar(scenarios, values, color=colors, alpha=0.8, edgecolor='black')
    if limit is not None:
        ax.axhline(y=limit, color='red', linestyle='--', linewidth=1.5, label='Spec Limit')
        ax.legend()
    ax.set_ylabel(title)
    ax.set_title(title)
    ax.grid(True, alpha=0.3, axis='y')
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f'{val:.4f}', ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.savefig('figures/fig08_scenario_comparison.png')
plt.savefig('figures/fig08_scenario_comparison.svg')
plt.close()
print("  -> Saved figures/fig08_scenario_comparison.png")

# Figure 9: Monte Carlo distributions
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle('Module 6: Monte Carlo Quality Prediction (n=500)',
             fontsize=14, fontweight='bold')

mc_metrics = ['warpage', 'sink', 'weight', 'shrinkage', 'stress']
mc_titles = ['Warpage (mm)', 'Sink Depth (mm)', 'Weight (g)', 'Shrinkage (%)', 'Residual Stress (MPa)']
mc_limits_vals = [0.5, 0.02, None, 1.0, None]

for idx, (key, title, lim) in enumerate(zip(mc_metrics, mc_titles, mc_limits_vals)):
    ax = axes[idx // 3, idx % 3]
    data = mc_results['raw_data'][key]
    ax.hist(data, bins=30, alpha=0.7, color='steelblue', edgecolor='white')
    mean_val = mc_results['statistics'][key]['mean']
    ax.axvline(x=mean_val, color='red', linestyle='--', linewidth=1.5, label=f'Mean: {mean_val:.4f}')
    if lim is not None:
        ax.axvline(x=lim, color='orange', linestyle='-', linewidth=2, label=f'Spec: {lim}')
    ax.set_xlabel(title)
    ax.set_ylabel('Frequency')
    ax.set_title(title)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

axes[1, 2].set_visible(False)
plt.tight_layout()
plt.savefig('figures/fig09_monte_carlo.png')
plt.savefig('figures/fig09_monte_carlo.svg')
plt.close()
print("  -> Saved figures/fig09_monte_carlo.png")

save_json({'scenarios': {k: v for k, v in scenario_results.items()},
           'monte_carlo': {k: v for k, v in mc_results.items() if k != 'raw_data'},
           'mc_statistics': mc_results['statistics']},
          'results/case_study_results.json')

# ============================================================
# Figure 10: Digital Twin Architecture Diagram
# ============================================================
print("\n  Generating architecture diagram...")

fig, ax = plt.subplots(figsize=(16, 10))
ax.set_xlim(0, 16)
ax.set_ylim(0, 10)
ax.axis('off')
ax.set_title('Injection Molding Digital Twin Architecture\n(Moldflow / OpenFOAM Integration)',
             fontsize=16, fontweight='bold', pad=20)

# Physical layer
rect_phys = plt.Rectangle((0.5, 7.5), 15, 2, facecolor='#E3F2FD', edgecolor='#1565C0', linewidth=2)
ax.add_patch(rect_phys)
ax.text(8, 9.2, 'PHYSICAL LAYER', fontsize=12, fontweight='bold', ha='center', color='#1565C0')
boxes_phys = [('Injection\nMachine', 1.5), ('Mold\n(Sensors)', 4.5), ('Cooling\nSystem', 7.5),
              ('Quality\nInspection', 10.5), ('MES/SCADA\nData', 13.5)]
for label, x in boxes_phys:
    r = plt.Rectangle((x - 0.8, 7.7), 1.6, 1.0, facecolor='white', edgecolor='#1565C0',
                       linewidth=1.5, zorder=5)
    ax.add_patch(r)
    ax.text(x, 8.2, label, fontsize=8, ha='center', va='center', zorder=6)

# Data layer
rect_data = plt.Rectangle((0.5, 5.5), 15, 1.5, facecolor='#FFF3E0', edgecolor='#E65100', linewidth=2)
ax.add_patch(rect_data)
ax.text(8, 6.7, 'DATA INTEGRATION LAYER', fontsize=11, fontweight='bold', ha='center', color='#E65100')
data_items = ['Pressure\nSensors', 'Temp\nSensors', 'Displacement\nSensors', 'Process\nLog', 'Historian\nDB']
for i, label in enumerate(data_items):
    x = 1.5 + i * 3
    ax.text(x, 5.9, label, fontsize=8, ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#E65100'))

# Simulation layer
rect_sim = plt.Rectangle((0.5, 3.0), 15, 2.0, facecolor='#E8F5E9', edgecolor='#2E7D32', linewidth=2)
ax.add_patch(rect_sim)
ax.text(8, 4.7, 'SIMULATION ENGINE LAYER', fontsize=11, fontweight='bold', ha='center', color='#2E7D32')
sim_items = [('Moldflow\nHele-Shaw', 2), ('OpenFOAM\n3D CFD', 5), ('Crystallization\nKinetics', 8),
             ('FEA: Stress\n& Warpage', 11), ('GP Surrogate\nModel', 14)]
for label, x in sim_items:
    r = plt.Rectangle((x - 1, 3.2), 2, 1.2, facecolor='white', edgecolor='#2E7D32',
                       linewidth=1.5, zorder=5)
    ax.add_patch(r)
    ax.text(x, 3.8, label, fontsize=8, ha='center', va='center', zorder=6)

# Intelligence layer
rect_int = plt.Rectangle((0.5, 0.5), 15, 2.0, facecolor='#F3E5F5', edgecolor='#6A1B9A', linewidth=2)
ax.add_patch(rect_int)
ax.text(8, 2.2, 'INTELLIGENCE & DECISION LAYER', fontsize=11, fontweight='bold', ha='center', color='#6A1B9A')
int_items = [('EnKF Data\nAssimilation', 2.5), ('Quality\nPrediction', 5.5),
             ('Process\nOptimization', 8.5), ('Anomaly\nDetection', 11.5), ('Dashboard\n& Alerts', 14)]
for label, x in int_items:
    ax.text(x, 1.2, label, fontsize=8, ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#6A1B9A', linewidth=1.5))

# Arrows between layers
for x in [3, 6, 9, 12]:
    ax.annotate('', xy=(x, 7.5), xytext=(x, 7.0),
                arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
    ax.annotate('', xy=(x, 5.5), xytext=(x, 5.0),
                arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
    ax.annotate('', xy=(x, 3.0), xytext=(x, 2.5),
                arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))

plt.savefig('figures/fig10_architecture.png')
plt.savefig('figures/fig10_architecture.svg')
plt.close()
print("  -> Saved figures/fig10_architecture.png")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
print("PIPELINE COMPLETE")
print("=" * 60)

all_results = {
    'flow': {'max_pressure_MPa': flow_results['max_pressure_MPa'],
             'avg_temp_C': flow_results['avg_melt_temp_C'],
             '3d_max_velocity': results_3d['max_velocity_m_s']},
    'cooling': cooling_summary,
    'stress_warpage': {
        'max_tensile_MPa': stress_warpage['stress']['max_tensile_MPa'],
        'max_compressive_MPa': stress_warpage['stress']['max_compressive_MPa'],
        'total_warpage_mm': stress_warpage['warpage']['warpage_total_mm']},
    'surrogate_model': training_result['cv_scores'],
    'data_assimilation': {
        'convergence': da_results['convergence_achieved'],
        'final_rmse': da_results['final_rmse']},
    'case_study': {
        'pass_rate_pct': mc_results['pass_rate_pct'],
        'scenarios': {k: v['quality_pass'] for k, v in scenario_results.items()}},
}
save_json(all_results, 'results/summary.json')
print("  All results saved to results/")
print("  All figures saved to figures/")
