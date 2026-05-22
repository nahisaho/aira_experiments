#!/usr/bin/env python3
"""
食品テクスチャ予測フレームワーク — 全モジュール実行スクリプト
"""
import sys, os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.viscoelastic_gel import (
    GeneralizedMaxwellModel, MaxwellElement, GeneralizedKelvinVoigtModel,
    KelvinVoigtElement, FractionalSpringPot, FractionalMaxwell,
    FEMesh2D, fem_uniaxial_compression, build_gel_model, fit_generalized_maxwell,
    POLYSACCHARIDE_PARAMS
)
from src.emulsion_rheology import (
    microstructure_to_rheology, krieger_dougherty, CGMDSimulation,
    droplet_size_distribution
)
from src.tpa_prediction import (
    TPAPredictionModel, simulate_tpa_curve, compute_tpa_from_curve
)
from src.oral_processing import simulate_oral_processing
from src.food_printing import (
    FoodInkRheology, PrintabilityPredictor, ExtrusionModel,
    optimize_printing_parameters, rheology_suitability_map
)
from src.plant_meat_design import (
    PlantProteinFormulation, HMECProcessConditions,
    hmec_texture_prediction, optimize_formulation, run_case_study,
    texture_similarity_score, REFERENCE_MEATS
)

os.makedirs('figures', exist_ok=True)
os.makedirs('results', exist_ok=True)
os.makedirs('data', exist_ok=True)
os.makedirs('logs', exist_ok=True)

all_results = {}

# ============================================================
# Module 1: 多糖類ゲルの粘弾性モデリング
# ============================================================
print("=" * 60)
print("Module 1: 多糖類ゲルの粘弾性モデリング")
print("=" * 60)

# 1a: 各多糖類のMaxwellモデル構築と比較
gel_types = ['κ-carrageenan', 'agar', 'gellan', 'alginate', 'pectin_LM']
t = np.logspace(-3, 3, 200)
omega = np.logspace(-2, 3, 200)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for gel_name in gel_types:
    model = build_gel_model(gel_name, concentration=1.5, temperature=25.0)
    G_t = model.relaxation_modulus(t)
    axes[0, 0].loglog(t, G_t, label=gel_name, linewidth=2)
    Gp = model.storage_modulus(omega)
    Gpp = model.loss_modulus(omega)
    axes[0, 1].loglog(omega, Gp, '-', label=f"{gel_name} G'", linewidth=1.5)
    axes[0, 1].loglog(omega, Gpp, '--', label=f"{gel_name} G''", linewidth=1.0, alpha=0.7)

axes[0, 0].set_xlabel('Time [s]'); axes[0, 0].set_ylabel('G(t) [Pa]')
axes[0, 0].set_title('Relaxation Modulus'); axes[0, 0].legend(fontsize=7)
axes[0, 1].set_xlabel('ω [rad/s]'); axes[0, 1].set_ylabel('Modulus [Pa]')
axes[0, 1].set_title("Storage (G') and Loss (G'') Moduli"); axes[0, 1].legend(fontsize=6)

# 1b: 濃度依存性
concentrations = np.linspace(0.5, 3.0, 6)
colors_conc = plt.cm.viridis(np.linspace(0, 1, len(concentrations)))
for c, col in zip(concentrations, colors_conc):
    m = build_gel_model('κ-carrageenan', concentration=c)
    axes[1, 0].loglog(omega, m.storage_modulus(omega), color=col, label=f'{c:.1f}%')
axes[1, 0].set_xlabel('ω [rad/s]'); axes[1, 0].set_ylabel("G' [Pa]")
axes[1, 0].set_title('κ-Carrageenan: Concentration Dependence'); axes[1, 0].legend(fontsize=8)

# 1c: FEM圧縮シミュレーション
mesh = FEMesh2D.generate_rectangle(10, 10, 8, 8)
model_fem = build_gel_model('agar', concentration=2.0)
fem_result = fem_uniaxial_compression(mesh, model_fem, strain_max=0.4, n_steps=30, dt=0.05)
axes[1, 1].plot(fem_result['strain']*100, fem_result['stress']/1000, 'b-', linewidth=2)
axes[1, 1].set_xlabel('Strain [%]'); axes[1, 1].set_ylabel('Stress [kPa]')
axes[1, 1].set_title(f'FEM Uniaxial Compression (Agar 2%)\nMesh: {fem_result["mesh_nodes"]} nodes, {fem_result["mesh_elements"]} elements')

plt.tight_layout()
plt.savefig('figures/fig1_viscoelastic_gel.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig1_viscoelastic_gel.svg', bbox_inches='tight')
plt.close()
print("  → figures/fig1_viscoelastic_gel.png saved")

# 1d: 分数階微分モデル
fig2, ax2 = plt.subplots(1, 2, figsize=(12, 5))
sp_a = FractionalSpringPot(V=1000, alpha=0.15)
sp_b = FractionalSpringPot(V=500, alpha=0.85)
fm = FractionalMaxwell(sp_a, sp_b)
ax2[0].loglog(omega, fm.storage_modulus(omega), 'b-', label="G' (Fractional Maxwell)", linewidth=2)
ax2[0].loglog(omega, fm.loss_modulus(omega), 'r--', label="G'' (Fractional Maxwell)", linewidth=2)
m_std = build_gel_model('κ-carrageenan', concentration=1.5)
ax2[0].loglog(omega, m_std.storage_modulus(omega), 'b:', label="G' (Prony series)", alpha=0.6)
ax2[0].loglog(omega, m_std.loss_modulus(omega), 'r:', label="G'' (Prony series)", alpha=0.6)
ax2[0].set_xlabel('ω [rad/s]'); ax2[0].set_ylabel('Modulus [Pa]')
ax2[0].set_title('Fractional Maxwell vs Prony Series'); ax2[0].legend()

# Fitting demo
t_fit = np.logspace(-2, 2, 50)
model_true = build_gel_model('gellan', concentration=1.0)
G_true = model_true.relaxation_modulus(t_fit)
noise = np.random.default_rng(42).normal(1, 0.03, len(G_true))
G_noisy = G_true * noise
fitted = fit_generalized_maxwell(t_fit, G_noisy, n_elements=3)
G_fitted = fitted.relaxation_modulus(t_fit)
ax2[1].loglog(t_fit, G_noisy, 'ko', markersize=4, label='Noisy data', alpha=0.6)
ax2[1].loglog(t_fit, G_fitted, 'r-', linewidth=2, label='Fitted (3 elements)')
ax2[1].loglog(t_fit, G_true, 'b--', linewidth=1, label='True model', alpha=0.7)
ax2[1].set_xlabel('Time [s]'); ax2[1].set_ylabel('G(t) [Pa]')
ax2[1].set_title('Model Fitting: Gellan Gum'); ax2[1].legend()
plt.tight_layout()
plt.savefig('figures/fig1b_fractional_fitting.png', dpi=300, bbox_inches='tight')
plt.close()
print("  → figures/fig1b_fractional_fitting.png saved")

# Save Module 1 results
mod1_results = {
    'gel_types_modeled': gel_types,
    'fem_mesh': {'nodes': int(fem_result['mesh_nodes']), 'elements': int(fem_result['mesh_elements'])},
    'fem_max_stress_kPa': float(np.max(fem_result['stress'])/1000),
    'fitting_R2': float(1 - np.sum((G_noisy - G_fitted)**2) / np.sum((G_noisy - G_noisy.mean())**2)),
}
all_results['module1_viscoelastic_gel'] = mod1_results
print(f"  Fitting R² = {mod1_results['fitting_R2']:.4f}")

# ============================================================
# Module 2: 乳化系の微視的構造とレオロジー
# ============================================================
print("\n" + "=" * 60)
print("Module 2: 乳化系のレオロジー")
print("=" * 60)

fig3, axes3 = plt.subplots(2, 2, figsize=(14, 10))

# 2a: 体積分率依存性
phis = [0.1, 0.2, 0.3, 0.4, 0.5]
colors_phi = plt.cm.cividis(np.linspace(0, 1, len(phis)))
for phi, col in zip(phis, colors_phi):
    rr = microstructure_to_rheology(phi, d_mean=5e-6, d_std=1e-6)
    axes3[0, 0].loglog(rr['omega'], rr['G_storage'], '-', color=col, label=f'φ={phi}')
    axes3[0, 0].loglog(rr['omega'], rr['G_loss'], '--', color=col, alpha=0.5)
axes3[0, 0].set_xlabel('ω [rad/s]'); axes3[0, 0].set_ylabel('Modulus [Pa]')
axes3[0, 0].set_title('Palierne Model: Volume Fraction Effect'); axes3[0, 0].legend(fontsize=8)

# 2b: Krieger-Dougherty
phi_arr = np.linspace(0.01, 0.6, 100)
for phi_max, ls in [(0.58, '--'), (0.64, '-'), (0.74, ':')]:
    eta_kd = krieger_dougherty(phi_arr, eta_0=1e-3, phi_max=phi_max)
    axes3[0, 1].semilogy(phi_arr, eta_kd, ls, linewidth=2, label=f'φ_max={phi_max}')
axes3[0, 1].set_xlabel('Volume Fraction φ'); axes3[0, 1].set_ylabel('η [Pa·s]')
axes3[0, 1].set_title('Krieger-Dougherty Viscosity Model'); axes3[0, 1].legend()

# 2c: 液滴サイズ分布の影響
d_means = [1e-6, 5e-6, 10e-6, 20e-6]
colors_d = plt.cm.magma(np.linspace(0.2, 0.8, len(d_means)))
for dm, col in zip(d_means, colors_d):
    rr = microstructure_to_rheology(0.3, d_mean=dm, d_std=dm*0.3)
    axes3[1, 0].loglog(rr['omega'], rr['G_storage'], '-', color=col,
                       label=f'd={dm*1e6:.0f}μm', linewidth=2)
axes3[1, 0].set_xlabel('ω [rad/s]'); axes3[1, 0].set_ylabel("G' [Pa]")
axes3[1, 0].set_title('Effect of Droplet Size on Rheology'); axes3[1, 0].legend()

# 2d: CG-MD simulation
print("  Running CG-MD simulation...")
cgmd = CGMDSimulation()
cgmd.initialize_emulsion(n_oil=100, n_water=400, n_surfactant=50, droplet_radius=4.0)
traj = cgmd.run(n_steps=500, save_interval=50)
axes3[1, 1].plot(traj['step'], traj['droplet_radius'], 'b-o', linewidth=2, markersize=4)
axes3[1, 1].set_xlabel('Step'); axes3[1, 1].set_ylabel('Mean Droplet Radius')
axes3[1, 1].set_title('CG-MD: Droplet Radius Evolution')

plt.tight_layout()
plt.savefig('figures/fig2_emulsion_rheology.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig2_emulsion_rheology.svg', bbox_inches='tight')
plt.close()
print("  → figures/fig2_emulsion_rheology.png saved")

mod2_results = {
    'volume_fractions_tested': phis,
    'droplet_sizes_um': [d*1e6 for d in d_means],
    'cgmd_final_radius': float(traj['droplet_radius'][-1]),
    'cgmd_n_beads': len(cgmd.beads),
}
all_results['module2_emulsion'] = mod2_results

# ============================================================
# Module 3: TPA予測モデル
# ============================================================
print("\n" + "=" * 60)
print("Module 3: TPA予測モデル")
print("=" * 60)

# 3a: TPA曲線シミュレーション
fig4, axes4 = plt.subplots(2, 2, figsize=(14, 10))

gels_tpa = {
    'Agar 2%': {'G_inf': 1000, 'G_el': [3000, 1000, 300], 'tau': [0.01, 0.5, 10]},
    'κ-Carrageenan 1.5%': {'G_inf': 200, 'G_el': [800, 400, 100], 'tau': [0.1, 1.0, 50]},
    'Pectin LM 2%': {'G_inf': 150, 'G_el': [500, 250, 80], 'tau': [0.05, 0.8, 25]},
}

tpa_results_summary = {}
for name, params in gels_tpa.items():
    disp, force = simulate_tpa_curve(params['G_inf'], params['G_el'], params['tau'],
                                     n_points=200, compression_ratio=0.5)
    axes4[0, 0].plot(disp, force, linewidth=2, label=name)
    tpa = compute_tpa_from_curve(disp, force)
    tpa_results_summary[name] = {
        'hardness_N': float(tpa.hardness),
        'cohesiveness': float(tpa.cohesiveness),
        'springiness': float(tpa.springiness),
        'gumminess_N': float(tpa.gumminess),
        'chewiness_N': float(tpa.chewiness),
        'resilience': float(tpa.resilience),
    }

axes4[0, 0].set_xlabel('Displacement [mm]'); axes4[0, 0].set_ylabel('Force [N]')
axes4[0, 0].set_title('Simulated TPA Curves'); axes4[0, 0].legend()

# 3b: ML model training
print("  Training TPA prediction ML model...")
tpa_model = TPAPredictionModel()
X_train, y_train = tpa_model.generate_training_data(n_samples=500)
ml_scores = tpa_model.train(X_train, y_train)

targets = list(ml_scores.keys())
r2_means = [ml_scores[t]['cv_r2_mean'] for t in targets]
r2_stds = [ml_scores[t]['cv_r2_std'] for t in targets]
bars = axes4[0, 1].bar(range(len(targets)), r2_means, yerr=r2_stds,
                        color=plt.cm.viridis(np.linspace(0.2, 0.8, len(targets))),
                        capsize=4)
axes4[0, 1].set_xticks(range(len(targets)))
axes4[0, 1].set_xticklabels(targets, rotation=45, ha='right', fontsize=8)
axes4[0, 1].set_ylabel('CV R² Score'); axes4[0, 1].set_title('ML Model Performance (5-fold CV)')
axes4[0, 1].set_ylim(0, 1.1)
for bar, v in zip(bars, r2_means):
    axes4[0, 1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                     f'{v:.3f}', ha='center', fontsize=8)

# 3c: Feature importance (hardness)
fi = ml_scores['hardness']['feature_importance']
features = list(fi.keys())
importances = list(fi.values())
sorted_idx = np.argsort(importances)
axes4[1, 0].barh([features[i] for i in sorted_idx],
                  [importances[i] for i in sorted_idx],
                  color='steelblue')
axes4[1, 0].set_xlabel('Feature Importance')
axes4[1, 0].set_title('Feature Importance for Hardness Prediction')

# 3d: Prediction example
test_compositions = [
    {'protein': 15, 'fat': 10, 'carbohydrate': 20, 'moisture': 50,
     'salt': 1.5, 'pH': 5.5, 'temperature': 120, 'heating_time': 30, 'cooling_rate': 5},
    {'protein': 25, 'fat': 5, 'carbohydrate': 10, 'moisture': 55,
     'salt': 1.0, 'pH': 6.0, 'temperature': 160, 'heating_time': 60, 'cooling_rate': 10},
]
pred_labels = ['Sample A', 'Sample B']
pred_hardness = []
pred_cohes = []
for comp in test_compositions:
    pred = tpa_model.predict(comp)
    pred_hardness.append(pred.hardness)
    pred_cohes.append(pred.cohesiveness)

x_pos = np.arange(len(pred_labels))
w = 0.35
axes4[1, 1].bar(x_pos - w/2, pred_hardness, w, label='Hardness [N]', color='coral')
ax_twin = axes4[1, 1].twinx()
ax_twin.bar(x_pos + w/2, pred_cohes, w, label='Cohesiveness', color='steelblue')
axes4[1, 1].set_xticks(x_pos)
axes4[1, 1].set_xticklabels(pred_labels)
axes4[1, 1].set_ylabel('Hardness [N]', color='coral')
ax_twin.set_ylabel('Cohesiveness', color='steelblue')
axes4[1, 1].set_title('TPA Prediction for Test Samples')
axes4[1, 1].legend(loc='upper left'); ax_twin.legend(loc='upper right')

plt.tight_layout()
plt.savefig('figures/fig3_tpa_prediction.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig3_tpa_prediction.svg', bbox_inches='tight')
plt.close()
print("  → figures/fig3_tpa_prediction.png saved")

mod3_results = {
    'tpa_simulation': tpa_results_summary,
    'ml_cv_r2': {t: float(ml_scores[t]['cv_r2_mean']) for t in targets},
    'prediction_examples': {
        'Sample_A': {'hardness': pred_hardness[0], 'cohesiveness': pred_cohes[0]},
        'Sample_B': {'hardness': pred_hardness[1], 'cohesiveness': pred_cohes[1]},
    }
}
all_results['module3_tpa'] = mod3_results

# ============================================================
# Module 4: 口腔内プロセシング
# ============================================================
print("\n" + "=" * 60)
print("Module 4: 口腔内プロセシング")
print("=" * 60)

fig5, axes5 = plt.subplots(2, 2, figsize=(14, 10))

# 4a: 異なる食品硬さでの咀嚼シミュレーション
food_types = {
    'Soft (H=50N)': {'hardness': 50, 'cohesiveness': 0.3, 'viscosity': 0.05},
    'Medium (H=200N)': {'hardness': 200, 'cohesiveness': 0.5, 'viscosity': 0.1},
    'Hard (H=500N)': {'hardness': 500, 'cohesiveness': 0.7, 'viscosity': 0.5},
}

oral_results = {}
for name, params in food_types.items():
    result = simulate_oral_processing(
        food_hardness=params['hardness'],
        food_cohesiveness=params['cohesiveness'],
        bolus_viscosity=params['viscosity'],
        n_chews=30
    )
    mast = result['mastication']
    axes5[0, 0].plot(mast['chew_number'], mast['mean_size'], '-o', label=name,
                     markersize=3, linewidth=2)
    oral_results[name] = {
        'swallow_trigger': result['swallow_trigger_chew'],
        'final_mean_size': float(mast['mean_size'][-1]),
        'final_d90': float(mast['d90'][-1]),
        'final_moisture': float(mast['moisture'][-1]),
        'transit_time': float(result['swallowing']['transit_time']),
        'residue_fraction': float(result['swallowing']['residue_fraction']),
    }

axes5[0, 0].set_xlabel('Chew Number'); axes5[0, 0].set_ylabel('Mean Particle Size [mm]')
axes5[0, 0].set_title('Particle Size Reduction During Mastication'); axes5[0, 0].legend()

# 4b: ボーラス特性
result_med = simulate_oral_processing(food_hardness=200, food_cohesiveness=0.5,
                                       bolus_viscosity=0.1, n_chews=30)
mast_med = result_med['mastication']
axes5[0, 1].plot(mast_med['chew_number'], mast_med['moisture'], 'b-', label='Moisture [%]', linewidth=2)
ax5_twin = axes5[0, 1].twinx()
ax5_twin.plot(mast_med['chew_number'], mast_med['bolus_cohesion'], 'r--', label='Bolus Cohesion', linewidth=2)
axes5[0, 1].set_xlabel('Chew Number'); axes5[0, 1].set_ylabel('Moisture [%]', color='b')
ax5_twin.set_ylabel('Bolus Cohesion', color='r')
axes5[0, 1].set_title('Bolus Formation During Mastication')
axes5[0, 1].legend(loc='upper left'); ax5_twin.legend(loc='center right')

# 4c: 嚥下シミュレーション
swal = result_med['swallowing']
axes5[1, 0].plot(swal['time']*1000, swal['position'], 'g-', linewidth=2)
axes5[1, 0].axhline(y=120, color='r', linestyle='--', alpha=0.5, label='Pharynx length')
axes5[1, 0].set_xlabel('Time [ms]'); axes5[1, 0].set_ylabel('Bolus Position [mm]')
axes5[1, 0].set_title(f"Swallowing: Transit Time = {swal['transit_time']*1000:.1f} ms")
axes5[1, 0].legend()

# 4d: 嚥下トリガーの比較
names_oral = list(oral_results.keys())
triggers = [oral_results[n]['swallow_trigger'] for n in names_oral]
axes5[1, 1].bar(range(len(names_oral)), triggers,
                color=['#2ecc71', '#3498db', '#e74c3c'])
axes5[1, 1].set_xticks(range(len(names_oral)))
axes5[1, 1].set_xticklabels(names_oral, fontsize=8)
axes5[1, 1].set_ylabel('Chews to Swallow Threshold')
axes5[1, 1].set_title('Swallow Trigger Point by Food Type')
for i, v in enumerate(triggers):
    axes5[1, 1].text(i, v + 0.3, str(v), ha='center', fontweight='bold')

plt.tight_layout()
plt.savefig('figures/fig4_oral_processing.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig4_oral_processing.svg', bbox_inches='tight')
plt.close()
print("  → figures/fig4_oral_processing.png saved")

all_results['module4_oral_processing'] = oral_results

# ============================================================
# Module 5: 3Dフードプリンティング
# ============================================================
print("\n" + "=" * 60)
print("Module 5: 3Dフードプリンティング")
print("=" * 60)

fig6, axes6 = plt.subplots(2, 2, figsize=(14, 10))

# 5a: インクレオロジー適性マップ
print("  Computing rheology suitability map...")
smap = rheology_suitability_map(n_points=25)
im = axes6[0, 0].contourf(np.log10(smap['yield_stress']),
                            np.log10(smap['G_storage']),
                            smap['printability_score'],
                            levels=20, cmap='viridis')
plt.colorbar(im, ax=axes6[0, 0], label='Printability Score')
axes6[0, 0].set_xlabel('log₁₀(Yield Stress [Pa])')
axes6[0, 0].set_ylabel("log₁₀(G' [Pa])")
axes6[0, 0].set_title('Printability Suitability Map')

# Mark optimal region
axes6[0, 0].axvline(x=np.log10(200), color='r', linestyle='--', alpha=0.5)
axes6[0, 0].axvline(x=np.log10(500), color='r', linestyle='--', alpha=0.5)

# 5b: 印刷パラメータ最適化
inks = {
    'Starch paste': FoodInkRheology(300, 5.0, 0.35, 2000, 400, recovery_time=8),
    'Cheese': FoodInkRheology(500, 10.0, 0.3, 5000, 1000, recovery_time=15),
    'Chocolate': FoodInkRheology(200, 3.0, 0.4, 1500, 300, recovery_time=5),
    'Puree': FoodInkRheology(50, 1.0, 0.5, 500, 150, recovery_time=3),
}

print_results = {}
for name, ink in inks.items():
    opt = optimize_printing_parameters(ink, target_layers=10)
    print_results[name] = opt
    print(f"  {name}: Score={opt['printability_scores']['overall']:.3f}, "
          f"Nozzle={opt['optimal_nozzle_diameter']:.2f}mm")

ink_names = list(print_results.keys())
overall_scores = [print_results[n]['printability_scores']['overall'] for n in ink_names]
axes6[0, 1].barh(ink_names, overall_scores,
                  color=plt.cm.viridis(np.linspace(0.2, 0.8, len(ink_names))))
axes6[0, 1].set_xlabel('Overall Printability Score')
axes6[0, 1].set_title('Ink Printability Comparison')
for i, v in enumerate(overall_scores):
    axes6[0, 1].text(v + 0.01, i, f'{v:.3f}', va='center')

# 5c: 形状保持性 vs 層数
ink_test = inks['Starch paste']
layers_range = np.arange(1, 30)
pred = PrintabilityPredictor()
sri_values = [pred.shape_retention_index(ink_test, n) for n in layers_range]
axes6[1, 0].plot(layers_range, sri_values, 'b-o', linewidth=2, markersize=3)
axes6[1, 0].set_xlabel('Number of Layers'); axes6[1, 0].set_ylabel('Shape Retention Index')
axes6[1, 0].set_title('Shape Retention vs Layer Count (Starch Paste)')
axes6[1, 0].axhline(y=0.9, color='r', linestyle='--', label='Threshold (SRI=0.9)')
axes6[1, 0].legend()

# 5d: 押出し力 vs 流量
flow_rates = np.linspace(1, 50, 50)
ext = ExtrusionModel(nozzle_diameter=1.0)
for name, ink in list(inks.items())[:3]:
    forces = [ext.extrusion_force(ink, fr) for fr in flow_rates]
    axes6[1, 1].plot(flow_rates, forces, linewidth=2, label=name)
axes6[1, 1].set_xlabel('Flow Rate [mm³/s]'); axes6[1, 1].set_ylabel('Extrusion Force [N]')
axes6[1, 1].set_title('Extrusion Force vs Flow Rate'); axes6[1, 1].legend()

plt.tight_layout()
plt.savefig('figures/fig5_food_printing.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig5_food_printing.svg', bbox_inches='tight')
plt.close()
print("  → figures/fig5_food_printing.png saved")

# Serialize print_results
mod5_serializable = {}
for name, pr in print_results.items():
    mod5_serializable[name] = {
        'optimal_nozzle_mm': float(pr['optimal_nozzle_diameter']),
        'optimal_speed_mm_s': float(pr['optimal_print_speed']),
        'optimal_flow_mm3_s': float(pr['optimal_flow_rate']),
        'extrusion_force_N': float(pr['extrusion_force']),
        'overall_score': float(pr['printability_scores']['overall']),
    }
all_results['module5_printing'] = mod5_serializable

# ============================================================
# Module 6: 植物性代替肉ケーススタディ
# ============================================================
print("\n" + "=" * 60)
print("Module 6: 植物性代替肉ケーススタディ")
print("=" * 60)

print("  Running optimization for 3 meat targets...")
case_results = run_case_study()

fig7, axes7 = plt.subplots(2, 2, figsize=(14, 10))

# 6a: 最適組成の比較
meat_names = ['beef_patty', 'chicken_breast', 'pork_sausage']
ingredients = ['soy_protein', 'pea_protein', 'wheat_gluten', 'starch',
               'fat_content', 'fiber', 'methylcellulose', 'moisture']
colors_ing = plt.cm.Set3(np.linspace(0, 1, len(ingredients)))

x_pos = np.arange(len(meat_names))
bottom = np.zeros(len(meat_names))
for ing, col in zip(ingredients, colors_ing):
    values = [case_results[m]['optimal_formulation'][ing] for m in meat_names]
    axes7[0, 0].bar(x_pos, values, 0.6, bottom=bottom, label=ing, color=col)
    bottom += values
axes7[0, 0].set_xticks(x_pos)
axes7[0, 0].set_xticklabels([m.replace('_', ' ').title() for m in meat_names], fontsize=9)
axes7[0, 0].set_ylabel('Content [%]')
axes7[0, 0].set_title('Optimized Plant-Based Formulations')
axes7[0, 0].legend(fontsize=7, loc='upper right')

# 6b: 類似度スコア (レーダーチャート)
categories = ['hardness', 'cohesiveness', 'springiness', 'juiciness', 'fiber_alignment']
n_cats = len(categories)
angles = np.linspace(0, 2*np.pi, n_cats, endpoint=False).tolist()
angles += angles[:1]

ax_radar = fig7.add_subplot(2, 2, 2, projection='polar')
for meat, col in zip(meat_names, ['#e74c3c', '#3498db', '#2ecc71']):
    sim = case_results[meat]['similarity_scores']
    values = [sim[c] for c in categories]
    values += values[:1]
    ax_radar.plot(angles, values, 'o-', linewidth=2, color=col,
                 label=meat.replace('_', ' ').title())
    ax_radar.fill(angles, values, alpha=0.1, color=col)
ax_radar.set_xticks(angles[:-1])
ax_radar.set_xticklabels(categories, fontsize=8)
ax_radar.set_ylim(0, 1)
ax_radar.set_title('Texture Similarity to Target Meat', pad=20)
ax_radar.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=7)

# 6c: 感度分析
sens = case_results['sensitivity_analysis']
axes7[1, 0].plot(sens['soy_protein'], sens['hardness'], 'b-o', linewidth=2, markersize=4)
axes7[1, 0].set_xlabel('Soy Protein Content [%]')
axes7[1, 0].set_ylabel('Predicted Hardness [N]', color='b')
ax7_twin = axes7[1, 0].twinx()
ax7_twin.plot(sens['soy_protein'], sens['similarity'], 'r--s', linewidth=2, markersize=4)
ax7_twin.set_ylabel('Similarity to Beef Patty', color='r')
axes7[1, 0].set_title('Sensitivity Analysis: Soy Protein Content')

# 6d: テクスチャ比較 (予測 vs 参照)
for i, meat in enumerate(meat_names):
    pred = case_results[meat]['predicted_texture']
    ref = REFERENCE_MEATS[meat]
    x = i
    axes7[1, 1].bar(x - 0.15, pred['hardness'], 0.3, color='steelblue',
                    label='Predicted' if i == 0 else None)
    axes7[1, 1].bar(x + 0.15, ref['hardness'], 0.3, color='coral',
                    label='Reference' if i == 0 else None)
axes7[1, 1].set_xticks(range(len(meat_names)))
axes7[1, 1].set_xticklabels([m.replace('_', ' ').title() for m in meat_names], fontsize=9)
axes7[1, 1].set_ylabel('Hardness [N]')
axes7[1, 1].set_title('Predicted vs Reference Meat Hardness')
axes7[1, 1].legend()

plt.tight_layout()
plt.savefig('figures/fig6_plant_meat.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/fig6_plant_meat.svg', bbox_inches='tight')
plt.close()
print("  → figures/fig6_plant_meat.png saved")

# Serialize case study results
mod6_serializable = {}
for meat in meat_names:
    cr = case_results[meat]
    mod6_serializable[meat] = {
        'optimal_formulation': {k: float(v) for k, v in cr['optimal_formulation'].items()},
        'similarity_overall': float(cr['similarity_scores']['overall']),
        'similarity_detail': {k: float(v) for k, v in cr['similarity_scores'].items()},
        'predicted_hardness': float(cr['predicted_texture']['hardness']),
        'predicted_DT': float(cr['predicted_texture']['degree_of_texturization']),
        'predicted_FAI': float(cr['predicted_texture']['fiber_alignment_index']),
    }
all_results['module6_plant_meat'] = mod6_serializable

# Print summary
for meat in meat_names:
    sim = case_results[meat]['similarity_scores']['overall']
    print(f"  {meat}: Similarity = {sim:.3f}")

# ============================================================
# Save all results
# ============================================================
print("\n" + "=" * 60)
print("Saving results...")
print("=" * 60)

with open('results/all_results.json', 'w') as f:
    json.dump(all_results, f, indent=2, default=str)
print("  → results/all_results.json saved")

# Save numerical data
np.savez('data/gel_relaxation.npz',
         time=t, omega=omega,
         gel_types=gel_types)
np.savez('data/fem_compression.npz', **fem_result)
print("  → data/gel_relaxation.npz saved")
print("  → data/fem_compression.npz saved")

# Process log
log_entries = [
    {"timestamp": datetime.now().isoformat(), "phase": "execute",
     "event_type": "run_completed", "actor": "co-scientist",
     "skill_or_tool": "food-texture-framework",
     "modules_executed": ["viscoelastic_gel", "emulsion_rheology",
                          "tpa_prediction", "oral_processing",
                          "food_printing", "plant_meat_design"],
     "files_written": [
         "figures/fig1_viscoelastic_gel.png", "figures/fig1b_fractional_fitting.png",
         "figures/fig2_emulsion_rheology.png", "figures/fig3_tpa_prediction.png",
         "figures/fig4_oral_processing.png", "figures/fig5_food_printing.png",
         "figures/fig6_plant_meat.png", "results/all_results.json",
         "data/gel_relaxation.npz", "data/fem_compression.npz"
     ],
     "status": "ok"}
]
with open('logs/process-log.jsonl', 'w') as f:
    for entry in log_entries:
        f.write(json.dumps(entry) + '\n')
print("  → logs/process-log.jsonl saved")

print("\n" + "=" * 60)
print("ALL MODULES COMPLETED SUCCESSFULLY")
print("=" * 60)
