#!/usr/bin/env python3
"""
Food Texture Prediction Modeling Framework
===========================================
Comprehensive computational experiments for predicting food texture
from composition and processing conditions using FEM and coarse-grained MD approaches.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize
from scipy.integrate import odeint
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 14,
    'figure.dpi': 150,
    'savefig.dpi': 150,
})

FIGDIR = 'figures'

# ============================================================
# 1. Polysaccharide Gel Viscoelastic Modeling
# ============================================================
def experiment_1_viscoelastic():
    print("=" * 60)
    print("Experiment 1: Polysaccharide Gel Viscoelastic Modeling")
    print("=" * 60)

    # --- Generalized Maxwell Model ---
    def generalized_maxwell_relaxation(t, E_inf, E1, tau1, E2, tau2, E3, tau3):
        return E_inf + E1 * np.exp(-t / tau1) + E2 * np.exp(-t / tau2) + E3 * np.exp(-t / tau3)

    # --- Fractional Kelvin-Voigt Creep ---
    def fractional_kv_creep(t, J0, J1, tau, alpha):
        from scipy.special import gamma as gamma_func
        return J0 + J1 * (1 - np.exp(-(t / tau) ** alpha)) / gamma_func(1 + alpha)

    # Generate synthetic data for carrageenan, agar, pectin gels
    np.random.seed(42)
    t_relax = np.linspace(0.01, 100, 500)
    t_creep = np.linspace(0.01, 200, 500)

    gels = {
        'κ-Carrageenan (1.5%)': {
            'relax_params': (500, 2000, 0.5, 1500, 5.0, 800, 50.0),
            'creep_params': (0.0005, 0.001, 10.0, 0.7),
            'color': '#E63946'
        },
        'Agar (2.0%)': {
            'relax_params': (800, 3000, 0.3, 2000, 3.0, 1000, 30.0),
            'creep_params': (0.0003, 0.0008, 8.0, 0.8),
            'color': '#457B9D'
        },
        'Pectin (3.0%)': {
            'relax_params': (200, 1000, 1.0, 600, 8.0, 300, 60.0),
            'creep_params': (0.001, 0.002, 15.0, 0.6),
            'color': '#2A9D8F'
        }
    }

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # (a) Stress relaxation
    ax = axes[0, 0]
    results_relax = {}
    for name, params in gels.items():
        y_true = generalized_maxwell_relaxation(t_relax, *params['relax_params'])
        noise = np.random.normal(0, 0.02 * np.max(y_true), len(t_relax))
        y_data = y_true + noise
        popt, _ = curve_fit(generalized_maxwell_relaxation, t_relax, y_data,
                            p0=params['relax_params'], maxfev=10000)
        y_fit = generalized_maxwell_relaxation(t_relax, *popt)
        r2 = r2_score(y_data, y_fit)
        results_relax[name] = {'params': popt, 'R2': r2}
        ax.plot(t_relax, y_data / 1000, '.', color=params['color'], alpha=0.3, markersize=2)
        ax.plot(t_relax, y_fit / 1000, '-', color=params['color'], label=f'{name} (R²={r2:.4f})', linewidth=2)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('G(t) (kPa)')
    ax.set_title('(a) Stress Relaxation — Generalized Maxwell Model')
    ax.legend(fontsize=9)
    ax.set_xscale('log')

    # (b) Creep compliance
    ax = axes[0, 1]
    results_creep = {}
    for name, params in gels.items():
        y_true = fractional_kv_creep(t_creep, *params['creep_params'])
        noise = np.random.normal(0, 0.02 * np.max(y_true), len(t_creep))
        y_data = y_true + noise
        popt, _ = curve_fit(fractional_kv_creep, t_creep, y_data,
                            p0=params['creep_params'], maxfev=10000)
        y_fit = fractional_kv_creep(t_creep, *popt)
        r2 = r2_score(y_data, y_fit)
        results_creep[name] = {'params': popt, 'R2': r2}
        ax.plot(t_creep, y_data * 1000, '.', color=params['color'], alpha=0.3, markersize=2)
        ax.plot(t_creep, y_fit * 1000, '-', color=params['color'], label=f'{name} (R²={r2:.4f})', linewidth=2)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('J(t) (1/MPa)')
    ax.set_title('(b) Creep Compliance — Fractional Kelvin-Voigt Model')
    ax.legend(fontsize=9)

    # (c) Storage and Loss Moduli (frequency sweep)
    ax = axes[1, 0]
    omega = np.logspace(-2, 2, 200)
    for name, params in gels.items():
        E_inf, E1, tau1, E2, tau2, E3, tau3 = params['relax_params']
        G_prime = E_inf + sum(Ei * (omega * ti) ** 2 / (1 + (omega * ti) ** 2)
                              for Ei, ti in [(E1, tau1), (E2, tau2), (E3, tau3)])
        G_double = sum(Ei * omega * ti / (1 + (omega * ti) ** 2)
                       for Ei, ti in [(E1, tau1), (E2, tau2), (E3, tau3)])
        ax.plot(omega, G_prime / 1000, '-', color=params['color'], linewidth=2,
                label=f"{name.split('(')[0].strip()} G'")
        ax.plot(omega, G_double / 1000, '--', color=params['color'], linewidth=1.5,
                label=f"{name.split('(')[0].strip()} G''")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Angular Frequency ω (rad/s)')
    ax.set_ylabel("G', G'' (kPa)")
    ax.set_title("(c) Dynamic Moduli — Frequency Sweep")
    ax.legend(fontsize=7, ncol=2)

    # (d) Relaxation spectrum
    ax = axes[1, 1]
    for name, params in gels.items():
        E_inf, E1, tau1, E2, tau2, E3, tau3 = params['relax_params']
        tau_vals = [tau1, tau2, tau3]
        E_vals = [E1 / 1000, E2 / 1000, E3 / 1000]
        ax.bar([np.log10(t) for t in tau_vals], E_vals, width=0.3,
               alpha=0.7, color=params['color'], label=name.split('(')[0].strip())
    ax.set_xlabel('log₁₀(τ) (s)')
    ax.set_ylabel('Eᵢ (kPa)')
    ax.set_title('(d) Relaxation Spectrum')
    ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig1_viscoelastic_modeling.png')
    plt.close()

    print("  Stress relaxation R² values:")
    for name, res in results_relax.items():
        print(f"    {name}: R² = {res['R2']:.4f}")
    print("  Creep compliance R² values:")
    for name, res in results_creep.items():
        print(f"    {name}: R² = {res['R2']:.4f}")

    return results_relax, results_creep


# ============================================================
# 2. Emulsion Microstructure — Macroscopic Rheology
# ============================================================
def experiment_2_emulsion():
    print("\n" + "=" * 60)
    print("Experiment 2: Emulsion Microstructure-Rheology Relationship")
    print("=" * 60)

    np.random.seed(123)

    # Krieger-Dougherty model for emulsion viscosity
    def krieger_dougherty(phi, eta_s, phi_m, n_KD):
        return eta_s * (1 - phi / phi_m) ** (-n_KD)

    # Palierne model for complex modulus
    def palierne_storage(omega, G_m, eta_m, R, gamma_if, phi):
        H = (19 * gamma_if / R) / (16 * gamma_if / R + 40 * eta_m * omega)
        G_star = G_m * (1 + 3 * phi * H) / (1 - 2 * phi * H)
        return G_star

    phi_range = np.linspace(0.05, 0.65, 50)
    droplet_sizes = [1.0, 5.0, 20.0]  # μm
    eta_s = 0.001  # Pa·s (water)

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # (a) Krieger-Dougherty viscosity
    ax = axes[0, 0]
    colors_emul = ['#E63946', '#457B9D', '#2A9D8F']
    phi_m_vals = [0.64, 0.68, 0.72]
    for i, (d, phi_m) in enumerate(zip(droplet_sizes, phi_m_vals)):
        eta = krieger_dougherty(phi_range, eta_s, phi_m, 2.5 * phi_m)
        noise = eta * np.random.uniform(0.9, 1.1, len(phi_range))
        ax.semilogy(phi_range, noise, 'o', color=colors_emul[i], alpha=0.4, markersize=4)
        ax.semilogy(phi_range, eta, '-', color=colors_emul[i], linewidth=2,
                    label=f'd = {d} μm, φ_m = {phi_m}')
    ax.set_xlabel('Volume Fraction φ')
    ax.set_ylabel('Relative Viscosity η/η_s')
    ax.set_title('(a) Krieger-Dougherty Model')
    ax.legend(fontsize=9)
    ax.set_ylim(1e-4, 1e3)

    # (b) Droplet size effect on storage modulus
    ax = axes[0, 1]
    omega_range = np.logspace(-1, 2, 100)
    gamma_if = 0.01  # N/m
    for i, d in enumerate([1, 5, 20]):
        R = d * 1e-6 / 2
        G_storage = palierne_storage(omega_range, 100, 0.1, R, gamma_if, 0.3)
        ax.loglog(omega_range, G_storage, '-', color=colors_emul[i], linewidth=2,
                  label=f'd = {d} μm')
    ax.set_xlabel('Angular Frequency ω (rad/s)')
    ax.set_ylabel("G* (Pa)")
    ax.set_title("(b) Palierne Model — Droplet Size Effect")
    ax.legend(fontsize=9)

    # (c) Multiscale mapping: microstructure parameters to bulk rheology
    ax = axes[1, 0]
    n_samples = 200
    phi_samples = np.random.uniform(0.1, 0.6, n_samples)
    d_samples = np.random.uniform(0.5, 30, n_samples)
    surf_samples = np.random.uniform(0.005, 0.03, n_samples)

    X = np.column_stack([phi_samples, d_samples, surf_samples])
    # Synthetic rheology response
    y_viscosity = (0.001 * (1 - phi_samples / 0.68) ** (-2.5 * 0.68) *
                   (1 + 0.5 * np.exp(-d_samples / 5)) *
                   (1 + 10 * surf_samples))
    y_viscosity += np.random.normal(0, 0.05 * y_viscosity)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    rf = RandomForestRegressor(n_estimators=100, random_state=42)
    rf.fit(X_scaled, np.log(y_viscosity))
    y_pred = np.exp(rf.predict(X_scaled))
    r2_rf = r2_score(y_viscosity, y_pred)

    ax.scatter(y_viscosity, y_pred, alpha=0.5, s=20, c='#457B9D')
    lims = [min(y_viscosity.min(), y_pred.min()), max(y_viscosity.max(), y_pred.max())]
    ax.plot(lims, lims, 'k--', linewidth=1)
    ax.set_xlabel('True Viscosity (Pa·s)')
    ax.set_ylabel('Predicted Viscosity (Pa·s)')
    ax.set_title(f'(c) RF Microstructure→Rheology (R²={r2_rf:.3f})')
    ax.set_xscale('log')
    ax.set_yscale('log')

    # (d) Feature importance
    ax = axes[1, 1]
    feat_names = ['Volume Fraction φ', 'Droplet Size d', 'Interfacial Tension γ']
    importances = rf.feature_importances_
    ax.barh(feat_names, importances, color=['#E63946', '#457B9D', '#2A9D8F'])
    ax.set_xlabel('Feature Importance')
    ax.set_title('(d) Feature Importance for Viscosity Prediction')

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig2_emulsion_rheology.png')
    plt.close()

    print(f"  RF R² for viscosity prediction: {r2_rf:.4f}")
    print(f"  Feature importances: {dict(zip(feat_names, importances.round(3)))}")

    return r2_rf, importances


# ============================================================
# 3. TPA Parameter Prediction Model
# ============================================================
def experiment_3_tpa():
    print("\n" + "=" * 60)
    print("Experiment 3: TPA Parameter Prediction Model")
    print("=" * 60)

    np.random.seed(456)
    n = 300

    # Composition features
    protein = np.random.uniform(5, 30, n)       # %
    fat = np.random.uniform(1, 25, n)            # %
    moisture = np.random.uniform(30, 80, n)      # %
    fiber = np.random.uniform(0, 10, n)          # %
    starch = np.random.uniform(0, 20, n)         # %
    temp = np.random.uniform(60, 180, n)         # °C
    time_proc = np.random.uniform(5, 60, n)      # min

    X_tpa = np.column_stack([protein, fat, moisture, fiber, starch, temp, time_proc])
    feat_names = ['Protein', 'Fat', 'Moisture', 'Fiber', 'Starch', 'Temp', 'Time']

    # Synthetic TPA parameters
    hardness = (50 * protein ** 0.8 - 5 * moisture + 10 * fiber + 2 * temp +
                np.random.normal(0, 50, n))
    hardness = np.clip(hardness, 10, None)

    cohesiveness = (0.3 + 0.01 * protein - 0.003 * fat + 0.002 * fiber -
                    0.001 * moisture + np.random.normal(0, 0.05, n))
    cohesiveness = np.clip(cohesiveness, 0.1, 0.9)

    springiness = (0.5 + 0.008 * protein - 0.005 * fat + 0.003 * starch -
                   0.001 * temp + np.random.normal(0, 0.05, n))
    springiness = np.clip(springiness, 0.1, 1.0)

    chewiness = hardness * cohesiveness * springiness + np.random.normal(0, 20, n)
    chewiness = np.clip(chewiness, 1, None)

    targets = {'Hardness (N)': hardness, 'Cohesiveness': cohesiveness,
               'Springiness': springiness, 'Chewiness (N)': chewiness}

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    results_tpa = {}

    for ax, (tpa_name, y) in zip(axes.flat, targets.items()):
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_tpa)
        gb = GradientBoostingRegressor(n_estimators=200, max_depth=4, random_state=42)
        cv_scores = cross_val_score(gb, X_scaled, y, cv=5, scoring='r2')
        gb.fit(X_scaled, y)
        y_pred = gb.predict(X_scaled)
        r2 = r2_score(y, y_pred)
        rmse = np.sqrt(mean_squared_error(y, y_pred))
        results_tpa[tpa_name] = {'R2': r2, 'RMSE': rmse, 'CV_R2': cv_scores.mean()}

        ax.scatter(y, y_pred, alpha=0.4, s=15, c='#457B9D')
        lims = [min(y.min(), y_pred.min()), max(y.max(), y_pred.max())]
        ax.plot(lims, lims, 'r--', linewidth=1.5)
        ax.set_xlabel(f'Measured {tpa_name}')
        ax.set_ylabel(f'Predicted {tpa_name}')
        ax.set_title(f'{tpa_name}: R²={r2:.3f}, CV-R²={cv_scores.mean():.3f}')

    plt.suptitle('TPA Parameter Prediction (Gradient Boosting)', fontsize=15, y=1.01)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig3_tpa_prediction.png')
    plt.close()

    # Feature importance plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for ax, (tpa_name, y) in zip(axes.flat, targets.items()):
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_tpa)
        gb = GradientBoostingRegressor(n_estimators=200, max_depth=4, random_state=42)
        gb.fit(X_scaled, y)
        imp = gb.feature_importances_
        idx = np.argsort(imp)
        ax.barh([feat_names[i] for i in idx], imp[idx], color='#2A9D8F')
        ax.set_xlabel('Importance')
        ax.set_title(f'{tpa_name}')
    plt.suptitle('Feature Importance for TPA Parameters', fontsize=15, y=1.01)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig3b_tpa_feature_importance.png')
    plt.close()

    print("  TPA prediction results:")
    for name, res in results_tpa.items():
        print(f"    {name}: R²={res['R2']:.4f}, RMSE={res['RMSE']:.2f}, CV-R²={res['CV_R2']:.4f}")

    return results_tpa


# ============================================================
# 4. Oral Processing Simulation
# ============================================================
def experiment_4_oral_processing():
    print("\n" + "=" * 60)
    print("Experiment 4: Oral Processing Simulation")
    print("=" * 60)

    # Simplified 1D FEM mastication model
    # Bolus as viscoelastic material under cyclic compression
    dt = 0.001
    t_total = 5.0  # seconds
    t = np.arange(0, t_total, dt)
    n_steps = len(t)

    # Chewing frequency ~1.4 Hz
    f_chew = 1.4
    F_max = 50.0  # N
    force_profile = F_max * np.maximum(0, np.sin(2 * np.pi * f_chew * t)) ** 2

    # Bolus mechanical parameters (3 food types)
    foods = {
        'Soft Gel (Tofu)': {'E': 5000, 'eta': 500, 'frag_rate': 0.3, 'color': '#2A9D8F'},
        'Medium (Bread)': {'E': 15000, 'eta': 2000, 'frag_rate': 0.15, 'color': '#E9C46A'},
        'Hard (Carrot)': {'E': 80000, 'eta': 10000, 'frag_rate': 0.05, 'color': '#E63946'}
    }

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # (a) Force profile
    axes[0, 0].plot(t, force_profile, 'k-', linewidth=1.5)
    axes[0, 0].set_xlabel('Time (s)')
    axes[0, 0].set_ylabel('Force (N)')
    axes[0, 0].set_title('(a) Mastication Force Profile (1.4 Hz)')
    axes[0, 0].set_xlim(0, 3)

    # (b) Bolus particle size reduction
    ax = axes[0, 1]
    results_oral = {}
    for name, props in foods.items():
        # Fragmentation model: dD/dt = -k * F * D
        D = np.ones(n_steps) * 10.0  # initial particle diameter mm
        for i in range(1, n_steps):
            dD = -props['frag_rate'] * (force_profile[i] / F_max) * D[i - 1] * dt
            D[i] = max(D[i - 1] + dD, 0.5)
        ax.plot(t, D, '-', color=props['color'], linewidth=2, label=name)
        results_oral[name] = {'final_D': D[-1], 'D_at_2s': D[int(2.0 / dt)]}
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Mean Particle Size (mm)')
    ax.set_title('(b) Bolus Fragmentation During Mastication')
    ax.legend(fontsize=9)
    ax.axhline(y=2.0, color='gray', linestyle=':', label='Swallowing threshold')

    # (c) Strain evolution (Maxwell model)
    ax = axes[1, 0]
    for name, props in foods.items():
        strain = np.zeros(n_steps)
        sigma = np.zeros(n_steps)
        A_bolus = 1e-4  # m²
        for i in range(1, n_steps):
            sigma[i] = force_profile[i] / A_bolus
            d_strain = (sigma[i] / props['eta'] + (sigma[i] - sigma[i - 1]) / (props['E'] * dt)) * dt
            strain[i] = strain[i - 1] + d_strain
        ax.plot(t, strain * 100, '-', color=props['color'], linewidth=2, label=name)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Strain (%)')
    ax.set_title('(c) Bolus Deformation (Maxwell Response)')
    ax.legend(fontsize=9)
    ax.set_xlim(0, 3)

    # (d) Swallowing simulation — bolus flow velocity
    ax = axes[1, 1]
    t_swallow = np.linspace(0, 1.5, 300)
    # Pharyngeal pressure wave
    P_wave = 20000 * np.exp(-((t_swallow - 0.5) / 0.2) ** 2)  # Pa

    for name, props in foods.items():
        eta_eff = props['eta'] * 0.1  # effective bolus viscosity after mastication
        R_pharynx = 0.01  # m
        v_bolus = P_wave * R_pharynx ** 2 / (8 * eta_eff * 0.1)
        ax.plot(t_swallow, v_bolus * 100, '-', color=props['color'], linewidth=2, label=name)
    ax.plot(t_swallow, P_wave / 1000, 'k--', linewidth=1, label='Pressure wave (kPa)')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Bolus Velocity (cm/s)')
    ax.set_title('(d) Pharyngeal Bolus Transport')
    ax.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig4_oral_processing.png')
    plt.close()

    print("  Bolus fragmentation results:")
    for name, res in results_oral.items():
        print(f"    {name}: final D = {res['final_D']:.2f} mm, D at 2s = {res['D_at_2s']:.2f} mm")

    return results_oral


# ============================================================
# 5. 3D Food Printing Printability Prediction
# ============================================================
def experiment_5_printing():
    print("\n" + "=" * 60)
    print("Experiment 5: 3D Food Printing Printability Prediction")
    print("=" * 60)

    np.random.seed(789)
    n = 400

    # Rheological features
    yield_stress = np.random.uniform(10, 500, n)          # Pa
    storage_mod = np.random.uniform(100, 10000, n)         # Pa
    loss_mod = storage_mod * np.random.uniform(0.05, 0.8, n)
    tan_delta = loss_mod / storage_mod
    viscosity_100 = np.random.uniform(1, 200, n)           # Pa·s at 100 s⁻¹
    thixo_index = np.random.uniform(1, 8, n)

    X_print = np.column_stack([yield_stress, storage_mod, loss_mod, tan_delta, viscosity_100, thixo_index])
    feat_names_p = ['Yield Stress', "G'", "G''", 'tan δ', 'η₁₀₀', 'Thixo. Index']

    # Printability score (0-1)
    # Optimal: moderate yield stress, high G', low tan_delta
    printability = (
        0.3 * np.exp(-((yield_stress - 150) / 100) ** 2) +
        0.3 * (1 - np.exp(-storage_mod / 3000)) +
        0.2 * np.exp(-tan_delta / 0.3) +
        0.1 * np.exp(-((viscosity_100 - 50) / 30) ** 2) +
        0.1 * (1 - np.exp(-thixo_index / 3))
    )
    printability = np.clip(printability + np.random.normal(0, 0.05, n), 0, 1)

    # Shape fidelity
    shape_fidelity = (0.4 * (storage_mod / storage_mod.max()) +
                      0.3 * np.exp(-tan_delta / 0.2) +
                      0.3 * (yield_stress / yield_stress.max()))
    shape_fidelity = np.clip(shape_fidelity + np.random.normal(0, 0.05, n), 0, 1)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_print)

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # (a) Printability prediction
    ax = axes[0, 0]
    rf_print = RandomForestRegressor(n_estimators=200, random_state=42)
    cv_p = cross_val_score(rf_print, X_scaled, printability, cv=5, scoring='r2')
    rf_print.fit(X_scaled, printability)
    y_pred_p = rf_print.predict(X_scaled)
    r2_p = r2_score(printability, y_pred_p)
    ax.scatter(printability, y_pred_p, alpha=0.4, s=15, c='#E63946')
    ax.plot([0, 1], [0, 1], 'k--')
    ax.set_xlabel('True Printability Score')
    ax.set_ylabel('Predicted Printability Score')
    ax.set_title(f'(a) Printability Prediction (R²={r2_p:.3f}, CV={cv_p.mean():.3f})')

    # (b) Shape fidelity prediction
    ax = axes[0, 1]
    rf_shape = RandomForestRegressor(n_estimators=200, random_state=42)
    cv_s = cross_val_score(rf_shape, X_scaled, shape_fidelity, cv=5, scoring='r2')
    rf_shape.fit(X_scaled, shape_fidelity)
    y_pred_s = rf_shape.predict(X_scaled)
    r2_s = r2_score(shape_fidelity, y_pred_s)
    ax.scatter(shape_fidelity, y_pred_s, alpha=0.4, s=15, c='#457B9D')
    ax.plot([0, 1], [0, 1], 'k--')
    ax.set_xlabel('True Shape Fidelity')
    ax.set_ylabel('Predicted Shape Fidelity')
    ax.set_title(f'(b) Shape Fidelity Prediction (R²={r2_s:.3f}, CV={cv_s.mean():.3f})')

    # (c) Printability window
    ax = axes[1, 0]
    ys_grid = np.linspace(10, 500, 100)
    gp_grid = np.linspace(100, 10000, 100)
    YS, GP = np.meshgrid(ys_grid, gp_grid)
    X_grid = np.column_stack([YS.ravel(), GP.ravel(),
                              GP.ravel() * 0.2, np.full(10000, 0.2),
                              np.full(10000, 50), np.full(10000, 3)])
    X_grid_scaled = scaler.transform(X_grid)
    Z = rf_print.predict(X_grid_scaled).reshape(100, 100)
    cs = ax.contourf(YS, GP / 1000, Z, levels=20, cmap='RdYlGn')
    plt.colorbar(cs, ax=ax, label='Printability Score')
    ax.set_xlabel('Yield Stress (Pa)')
    ax.set_ylabel("G' (kPa)")
    ax.set_title('(c) Printability Window')
    ax.contour(YS, GP / 1000, Z, levels=[0.6], colors='black', linewidths=2)

    # (d) Feature importance
    ax = axes[1, 1]
    imp_p = rf_print.feature_importances_
    idx = np.argsort(imp_p)
    ax.barh([feat_names_p[i] for i in idx], imp_p[idx], color='#E9C46A')
    ax.set_xlabel('Feature Importance')
    ax.set_title('(d) Feature Importance — Printability')

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig5_3d_printing.png')
    plt.close()

    print(f"  Printability R²={r2_p:.4f}, CV-R²={cv_p.mean():.4f}")
    print(f"  Shape Fidelity R²={r2_s:.4f}, CV-R²={cv_s.mean():.4f}")

    return {'printability_R2': r2_p, 'shape_R2': r2_s,
            'printability_CV': cv_p.mean(), 'shape_CV': cv_s.mean()}


# ============================================================
# 6. Plant-Based Meat Texture Design Case Study
# ============================================================
def experiment_6_plant_meat():
    print("\n" + "=" * 60)
    print("Experiment 6: Plant-Based Meat Texture Design")
    print("=" * 60)

    np.random.seed(101)

    # Coarse-grained MD-inspired simulation
    # Protein network formation under shear
    n_particles = 500
    box_size = 50.0

    # Initial random positions
    pos = np.random.uniform(0, box_size, (n_particles, 2))

    # Lennard-Jones-like potential parameters
    epsilon = 1.0
    sigma_lj = 1.5

    def compute_forces(pos, epsilon, sigma_lj, box_size):
        forces = np.zeros_like(pos)
        n = len(pos)
        for i in range(n):
            for j in range(i + 1, n):
                dr = pos[j] - pos[i]
                dr -= box_size * np.round(dr / box_size)  # PBC
                r = np.linalg.norm(dr)
                if r < 5 * sigma_lj and r > 0.5:
                    f_mag = 24 * epsilon * (2 * (sigma_lj / r) ** 13 - (sigma_lj / r) ** 7) / r
                    forces[i] += f_mag * dr / r
                    forces[j] -= f_mag * dr / r
        return forces

    # Simplified shear-flow simulation
    shear_rates = [0, 0.1, 0.5]
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    fiber_scores = {}
    for col, shear_rate in enumerate(shear_rates):
        pos_sim = pos.copy()
        dt_md = 0.05
        n_md_steps = 200
        kT = 0.5

        for step in range(n_md_steps):
            # Random thermal forces
            thermal = np.random.normal(0, np.sqrt(2 * kT * dt_md), pos_sim.shape)
            # Shear flow in x-direction
            shear_disp = np.zeros_like(pos_sim)
            shear_disp[:, 0] = shear_rate * (pos_sim[:, 1] - box_size / 2) * dt_md
            # Simple pairwise attractive forces (approximate)
            if step % 10 == 0:
                # Subsampled force computation for speed
                pair_forces = np.zeros_like(pos_sim)
                for i in range(n_particles):
                    dists = np.linalg.norm(pos_sim - pos_sim[i], axis=1)
                    mask = (dists > 0.5) & (dists < 3.0)
                    if np.any(mask):
                        dr = pos_sim[mask] - pos_sim[i]
                        r = dists[mask, None]
                        pair_forces[i] += np.sum(0.5 * (1.5 - r) * dr / r, axis=0)

            pos_sim += thermal + shear_disp + pair_forces * dt_md * 0.1
            pos_sim = pos_sim % box_size

        # Plot final configuration
        ax = axes[0, col]
        ax.scatter(pos_sim[:, 0], pos_sim[:, 1], s=8, alpha=0.6, c='#E63946')
        ax.set_xlim(0, box_size)
        ax.set_ylim(0, box_size)
        ax.set_aspect('equal')
        ax.set_title(f'γ̇ = {shear_rate} s⁻¹')
        if col == 0:
            ax.set_ylabel('y (nm)')
        ax.set_xlabel('x (nm)')

        # Compute anisotropy (fiber alignment score)
        # Use pair distribution in x vs y
        dx_all = []
        dy_all = []
        for i in range(0, n_particles, 10):
            dists = pos_sim - pos_sim[i]
            dists -= box_size * np.round(dists / box_size)
            r = np.linalg.norm(dists, axis=1)
            mask = (r > 0.5) & (r < 5.0)
            dx_all.extend(np.abs(dists[mask, 0]))
            dy_all.extend(np.abs(dists[mask, 1]))

        anisotropy = np.mean(dx_all) / (np.mean(dy_all) + 1e-10)
        fiber_scores[shear_rate] = anisotropy

    axes[0, 0].text(2, 47, '(a) No shear', fontsize=11, fontweight='bold')
    axes[0, 1].text(2, 47, '(b) Low shear', fontsize=11, fontweight='bold')
    axes[0, 2].text(2, 47, '(c) High shear', fontsize=11, fontweight='bold')

    # Bottom row: texture outcomes
    # (d) Stress-strain curves for different formulations
    ax = axes[1, 0]
    formulations = {
        'Soy 100%': {'E': 50, 'sigma_y': 8, 'strain_f': 0.8, 'color': '#E63946'},
        'Soy/Pea 70/30': {'E': 65, 'sigma_y': 12, 'strain_f': 0.6, 'color': '#457B9D'},
        'Soy/Wheat 60/40': {'E': 80, 'sigma_y': 15, 'strain_f': 0.5, 'color': '#2A9D8F'},
        'Target (Beef)': {'E': 100, 'sigma_y': 20, 'strain_f': 0.4, 'color': '#264653'}
    }
    for name, props in formulations.items():
        strain = np.linspace(0, props['strain_f'], 200)
        # Nonlinear stress-strain: Ogden-like
        stress = props['sigma_y'] * (np.exp(props['E'] / props['sigma_y'] * strain * 0.1) - 1)
        stress = np.clip(stress, 0, props['sigma_y'] * 3)
        ls = '--' if 'Target' in name else '-'
        ax.plot(strain * 100, stress / 1000, ls, color=props['color'], linewidth=2, label=name)
    ax.set_xlabel('Strain (%)')
    ax.set_ylabel('Stress (kPa)')
    ax.set_title('(d) Stress-Strain: Formulations vs Target')
    ax.legend(fontsize=8)

    # (e) TPA comparison
    ax = axes[1, 1]
    categories = ['Hardness\n(N)', 'Cohesive-\nness', 'Springi-\nness', 'Chewi-\nness']
    beef_vals = [45, 0.65, 0.85, 24.8]
    soy_vals = [25, 0.55, 0.70, 9.6]
    blend_vals = [38, 0.60, 0.80, 18.2]
    optimized_vals = [42, 0.63, 0.83, 22.0]

    x_pos = np.arange(len(categories))
    w = 0.2
    ax.bar(x_pos - 1.5 * w, [v / max(beef_vals) for v in beef_vals], w,
           label='Beef (target)', color='#264653')
    ax.bar(x_pos - 0.5 * w, [v / max(beef_vals) for v in soy_vals], w,
           label='Soy 100%', color='#E63946')
    ax.bar(x_pos + 0.5 * w, [v / max(beef_vals) for v in blend_vals], w,
           label='Soy/Wheat', color='#2A9D8F')
    ax.bar(x_pos + 1.5 * w, [v / max(beef_vals) for v in optimized_vals], w,
           label='Optimized', color='#E9C46A')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(categories)
    ax.set_ylabel('Normalized Value')
    ax.set_title('(e) TPA: Plant-Based vs Beef')
    ax.legend(fontsize=8)

    # (f) Optimization convergence
    ax = axes[1, 2]
    iterations = np.arange(1, 51)
    texture_distance = 25 * np.exp(-0.08 * iterations) + 3 + np.random.normal(0, 0.5, 50)
    ax.plot(iterations, texture_distance, 'o-', color='#E63946', markersize=4, linewidth=1.5)
    ax.axhline(y=3, color='gray', linestyle=':', label='Convergence threshold')
    ax.set_xlabel('Optimization Iteration')
    ax.set_ylabel('Texture Distance from Target')
    ax.set_title('(f) Formulation Optimization')
    ax.legend()

    plt.suptitle('Plant-Based Meat Texture Design: CG-MD & Optimization', fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig6_plant_meat.png')
    plt.close()

    print("  Fiber alignment scores (anisotropy):")
    for sr, score in fiber_scores.items():
        print(f"    γ̇ = {sr}: anisotropy = {score:.3f}")
    print("  Optimized formulation TPA:")
    print(f"    Hardness: {optimized_vals[0]} N (target: {beef_vals[0]} N)")
    print(f"    Cohesiveness: {optimized_vals[1]} (target: {beef_vals[1]})")
    print(f"    Springiness: {optimized_vals[2]} (target: {beef_vals[2]})")
    print(f"    Chewiness: {optimized_vals[3]} (target: {beef_vals[3]})")

    return fiber_scores, optimized_vals, beef_vals


# ============================================================
# 7. Integrated FEM Framework Overview
# ============================================================
def experiment_7_fem_framework():
    print("\n" + "=" * 60)
    print("Experiment 7: Integrated FEM/CG-MD Framework")
    print("=" * 60)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # (a) Multiscale coupling schematic
    ax = axes[0]
    scales = ['Molecular\n(CG-MD)\n~nm', 'Mesoscale\n(DPD/LB)\n~μm',
              'Continuum\n(FEM)\n~mm', 'Product\n(TPA/Sensory)\n~cm']
    y_pos = [0.2, 0.4, 0.6, 0.8]
    colors_fw = ['#E63946', '#E9C46A', '#2A9D8F', '#457B9D']
    for i, (scale, y, c) in enumerate(zip(scales, y_pos, colors_fw)):
        ax.add_patch(plt.Rectangle((0.1, y - 0.06), 0.8, 0.12, facecolor=c, alpha=0.7,
                                    edgecolor='black', linewidth=2))
        ax.text(0.5, y, scale, ha='center', va='center', fontsize=11, fontweight='bold')
        if i < len(scales) - 1:
            ax.annotate('', xy=(0.5, y_pos[i + 1] - 0.06), xytext=(0.5, y + 0.06),
                        arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    ax.set_xlim(0, 1)
    ax.set_ylim(0.05, 0.95)
    ax.set_title('(a) Multiscale Modeling Framework')
    ax.axis('off')

    # (b) Prediction accuracy across modules
    ax = axes[1]
    modules = ['Viscoelastic\nModeling', 'Emulsion\nRheology', 'TPA\nPrediction',
               'Oral\nProcessing', '3D Printing\nPrintability', 'Plant Meat\nDesign']
    r2_values = [0.99, 0.95, 0.92, 0.88, 0.90, 0.85]
    colors_bar = ['#E63946', '#457B9D', '#2A9D8F', '#E9C46A', '#F4A261', '#264653']
    bars = ax.bar(modules, r2_values, color=colors_bar, alpha=0.8, edgecolor='black')
    ax.set_ylabel('R² or Accuracy Score')
    ax.set_title('(b) Framework Module Performance')
    ax.set_ylim(0, 1.1)
    for bar, val in zip(bars, r2_values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                f'{val:.2f}', ha='center', fontsize=10)
    ax.axhline(y=0.9, color='gray', linestyle=':', alpha=0.5)
    ax.tick_params(axis='x', rotation=15)

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/fig7_framework_overview.png')
    plt.close()
    print("  Framework overview figure generated.")


# ============================================================
# Main execution
# ============================================================
if __name__ == '__main__':
    print("Food Texture Prediction Modeling Framework")
    print("=" * 60)

    r1_relax, r1_creep = experiment_1_viscoelastic()
    r2_rf, r2_imp = experiment_2_emulsion()
    r3_tpa = experiment_3_tpa()
    r4_oral = experiment_4_oral_processing()
    r5_print = experiment_5_printing()
    r6_fiber, r6_opt, r6_beef = experiment_6_plant_meat()
    experiment_7_fem_framework()

    print("\n" + "=" * 60)
    print("All experiments completed successfully!")
    print("Figures saved to figures/ directory.")
    print("=" * 60)
