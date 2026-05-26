#!/usr/bin/env python3
"""
Automated Optimization System for Continuous Flow Synthesis Reactions
=====================================================================
Covers: CFD simulation, RTD analysis, Bayesian optimization, online analytics
feedback control, scale-up design, and pharmaceutical case study.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import integrate, optimize, stats, signal
from scipy.special import erfc
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel, WhiteKernel
import os, json, warnings
warnings.filterwarnings('ignore')

FIGDIR = 'figures'
os.makedirs(FIGDIR, exist_ok=True)

np.random.seed(42)

# ============================================================
# 1. Microreactor CFD Flow Field Simulation (2D simplified)
# ============================================================
def simulate_cfd():
    """Simulate 2D velocity field in a serpentine microreactor channel."""
    # Channel dimensions
    L = 0.05   # channel length (m)
    W = 0.001  # channel width (m)
    nx, ny = 200, 50
    x = np.linspace(0, L, nx)
    y = np.linspace(-W/2, W/2, ny)
    X, Y = np.meshgrid(x, y)

    # Parabolic (Poiseuille) flow profile
    u_max = 0.05  # m/s
    U = u_max * (1 - (2*Y/W)**2)
    V = np.zeros_like(U)

    # Add secondary flow (Dean vortices) in curved sections
    n_bends = 4
    for i in range(n_bends):
        x_bend = L * (i + 0.5) / n_bends
        sigma_x = L / (n_bends * 4)
        dean_strength = 0.008 * (-1)**i
        V += dean_strength * np.sin(2 * np.pi * Y / W) * np.exp(-((X - x_bend)**2) / (2 * sigma_x**2))

    speed = np.sqrt(U**2 + V**2)

    # Pressure field (simplified)
    dp_dx = -12 * 1e-3 * u_max / (W**2)  # mu=1e-3 Pa.s (water)
    P = dp_dx * X + 101325

    # --- Figure 1: Velocity field ---
    fig, axes = plt.subplots(2, 1, figsize=(12, 6))

    ax = axes[0]
    skip = 5
    c = ax.contourf(X*1000, Y*1000, speed*1000, levels=30, cmap='viridis')
    ax.quiver(X[::3, ::skip]*1000, Y[::3, ::skip]*1000,
              U[::3, ::skip], V[::3, ::skip], color='white', alpha=0.7, scale=0.5)
    plt.colorbar(c, ax=ax, label='Velocity (mm/s)')
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_title('(a) Velocity Field in Serpentine Microreactor')

    ax = axes[1]
    c2 = ax.contourf(X*1000, Y*1000, (P - 101325)/1000, levels=30, cmap='coolwarm')
    plt.colorbar(c2, ax=ax, label='Gauge Pressure (kPa)')
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_title('(b) Pressure Distribution')

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/cfd_velocity_field.png', dpi=150, bbox_inches='tight')
    plt.close()

    # --- Figure 2: Velocity profiles at different positions ---
    fig, ax = plt.subplots(figsize=(8, 5))
    positions = [0.1, 0.3, 0.5, 0.7, 0.9]
    for pos in positions:
        idx = int(pos * (nx - 1))
        ax.plot(y*1000, U[:, idx]*1000, label=f'x = {x[idx]*1000:.1f} mm')
    ax.set_xlabel('y (mm)')
    ax.set_ylabel('Velocity (mm/s)')
    ax.set_title('Velocity Profiles at Different Axial Positions')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/velocity_profiles.png', dpi=150, bbox_inches='tight')
    plt.close()

    print("[CFD] Velocity field simulation complete.")
    print(f"  Max velocity: {speed.max()*1000:.2f} mm/s")
    print(f"  Pressure drop: {abs(dp_dx)*L/1000:.2f} kPa")
    return {'u_max': speed.max(), 'dp': abs(dp_dx)*L}


# ============================================================
# 2. Residence Time Distribution (RTD)
# ============================================================
def simulate_rtd():
    """Compute RTD for different reactor models and compare with experimental data."""
    t = np.linspace(0, 60, 1000)  # time in seconds
    tau = 20.0  # mean residence time (s)

    # Plug Flow Reactor (PFR)
    E_pfr = np.zeros_like(t)
    idx_tau = np.argmin(np.abs(t - tau))
    E_pfr[idx_tau] = 1.0 / (t[1] - t[0])

    # CSTR
    E_cstr = (1/tau) * np.exp(-t/tau)

    # Laminar flow reactor
    E_laminar = np.where(t >= tau/2, tau**2 / (2 * t**3), 0)

    # Tanks-in-series (N=5, 10, 20)
    def tanks_in_series(t, N, tau):
        theta = t / tau
        from math import factorial
        return (N/tau) * (N*theta)**(N-1) / factorial(N-1) * np.exp(-N*theta)

    E_tis5 = np.array([tanks_in_series(ti, 5, tau) for ti in t])
    E_tis10 = np.array([tanks_in_series(ti, 10, tau) for ti in t])
    E_tis20 = np.array([tanks_in_series(ti, 20, tau) for ti in t])

    # Axial dispersion model
    Pe = 50  # Peclet number
    E_disp = np.where(t > 0,
        1/(2*np.sqrt(np.pi*t/tau/Pe)) * np.exp(-Pe*(1-t/tau)**2/(4*t/tau)),
        0)

    # Simulated experimental data (axial dispersion + noise)
    Pe_exp = 40
    E_exp_true = np.where(t > 0,
        1/(2*np.sqrt(np.pi*t/tau/Pe_exp)) * np.exp(-Pe_exp*(1-t/tau)**2/(4*t/tau)),
        0)
    E_exp_true /= np.trapz(E_exp_true, t)
    noise = np.random.normal(0, 0.001, len(t))
    E_exp = np.maximum(E_exp_true + noise, 0)
    E_exp /= np.trapz(E_exp, t)

    # --- Figure 3: RTD comparison ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.plot(t, E_cstr, '--', label='CSTR (N=1)', alpha=0.8)
    ax.plot(t, E_tis5, '-', label='Tanks-in-Series (N=5)', alpha=0.8)
    ax.plot(t, E_tis10, '-', label='Tanks-in-Series (N=10)', alpha=0.8)
    ax.plot(t, E_tis20, '-', label='Tanks-in-Series (N=20)', alpha=0.8)
    ax.plot(t, E_laminar, '-.', label='Laminar Flow', alpha=0.8)
    ax.axvline(tau, color='k', linestyle=':', alpha=0.5, label=f'τ = {tau} s')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('E(t) (1/s)')
    ax.set_title('(a) RTD for Different Reactor Models')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 60)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(t, E_disp, 'b-', label=f'Axial Dispersion (Pe={Pe})', lw=2)
    ax.scatter(t[::20], E_exp[::20], c='red', s=15, label='Experimental Data', zorder=5)
    ax.plot(t, E_exp_true, 'r--', alpha=0.5, label=f'Fit (Pe={Pe_exp})')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('E(t) (1/s)')
    ax.set_title('(b) Experimental vs Model RTD')
    ax.legend(fontsize=8)
    ax.set_xlim(0, 60)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/rtd_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Compute moments
    mean_exp = np.trapz(t * E_exp, t)
    var_exp = np.trapz((t - mean_exp)**2 * E_exp, t)

    # --- Figure 4: Cumulative RTD ---
    fig, ax = plt.subplots(figsize=(8, 5))
    F_cstr = 1 - np.exp(-t/tau)
    F_tis10 = np.array([np.trapz(E_tis10[:i+1], t[:i+1]) for i in range(len(t))])
    F_disp = np.array([np.trapz(E_disp[:i+1], t[:i+1]) for i in range(len(t))])
    F_exp = np.array([np.trapz(E_exp[:i+1], t[:i+1]) for i in range(len(t))])

    ax.plot(t, F_cstr, '--', label='CSTR')
    ax.plot(t, F_tis10, '-', label='TIS (N=10)')
    ax.plot(t, F_disp, '-', label=f'Axial Dispersion (Pe={Pe})')
    ax.scatter(t[::30], F_exp[::30], c='red', s=15, label='Experimental', zorder=5)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('F(t)')
    ax.set_title('Cumulative Residence Time Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/rtd_cumulative.png', dpi=150, bbox_inches='tight')
    plt.close()

    print("[RTD] Analysis complete.")
    print(f"  Mean residence time: {mean_exp:.2f} s")
    print(f"  Variance: {var_exp:.2f} s²")
    print(f"  Estimated Pe: {Pe_exp}")
    return {'tau': mean_exp, 'var': var_exp, 'Pe': Pe_exp}


# ============================================================
# 3. Bayesian Optimization of Reaction Conditions
# ============================================================
def reaction_yield(params):
    """Simulated reaction yield as function of (T, flow_rate, conc, catalyst_loading).
    Models an amide coupling reaction in continuous flow."""
    T, flow_rate, conc, cat = params
    # Normalized inputs
    T_n = (T - 80) / 40       # 60-120 °C
    fr_n = (flow_rate - 0.5) / 0.5  # 0.1-1.0 mL/min
    c_n = (conc - 0.3) / 0.3  # 0.05-0.6 M
    cat_n = (cat - 0.05) / 0.05  # 0.01-0.10 equiv

    # Complex response surface
    yield_val = (
        75.0
        + 12 * np.exp(-0.5 * (T_n - 0.3)**2)
        - 5 * fr_n**2
        + 8 * np.exp(-2 * (c_n - 0.2)**2)
        + 6 * cat_n * (1 - cat_n)
        - 3 * T_n * fr_n
        + 4 * c_n * cat_n
        - 2 * T_n**2 * c_n
    )
    # Add noise
    yield_val += np.random.normal(0, 1.5)
    return min(max(yield_val, 0), 100)


def bayesian_optimization():
    """Run Bayesian optimization for reaction conditions."""
    # Parameter bounds
    bounds = {
        'Temperature (°C)': (60, 120),
        'Flow rate (mL/min)': (0.1, 1.0),
        'Concentration (M)': (0.05, 0.6),
        'Catalyst (equiv)': (0.01, 0.10)
    }
    param_names = list(bounds.keys())
    lb = np.array([b[0] for b in bounds.values()])
    ub = np.array([b[1] for b in bounds.values()])

    # Initial random samples
    n_init = 10
    n_iter = 40
    X_init = np.random.uniform(lb, ub, size=(n_init, 4))
    Y_init = np.array([reaction_yield(x) for x in X_init])

    X_all = list(X_init)
    Y_all = list(Y_init)

    # GP model
    kernel = ConstantKernel(1.0) * Matern(length_scale=np.ones(4), nu=2.5) + WhiteKernel(noise_level=1.0)

    best_yields = [max(Y_all)]
    acq_values = []

    for i in range(n_iter):
        X_train = np.array(X_all)
        Y_train = np.array(Y_all)

        # Normalize
        X_norm = (X_train - lb) / (ub - lb)

        gpr = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5, normalize_y=True)
        gpr.fit(X_norm, Y_train)

        # Expected Improvement acquisition function
        y_best = max(Y_train)

        def neg_ei(x_norm):
            x_norm = x_norm.reshape(1, -1)
            mu, sigma = gpr.predict(x_norm, return_std=True)
            if sigma < 1e-6:
                return 0.0
            z = (mu - y_best) / sigma
            ei = sigma * (z * stats.norm.cdf(z) + stats.norm.pdf(z))
            return -ei.item()

        # Multi-start optimization of acquisition
        best_ei = 0
        best_x = None
        for _ in range(50):
            x0 = np.random.uniform(0, 1, 4)
            res = optimize.minimize(neg_ei, x0, bounds=[(0,1)]*4, method='L-BFGS-B')
            if -res.fun > best_ei:
                best_ei = -res.fun
                best_x = res.x

        acq_values.append(best_ei)
        x_new = best_x * (ub - lb) + lb
        y_new = reaction_yield(x_new)

        X_all.append(x_new)
        Y_all.append(y_new)
        best_yields.append(max(Y_all))

    X_all = np.array(X_all)
    Y_all = np.array(Y_all)
    best_idx = np.argmax(Y_all)

    # --- Figure 5: Optimization convergence ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    ax.plot(range(len(Y_all)), Y_all, 'o', alpha=0.5, markersize=4, label='Observed')
    ax.plot(range(n_init, n_init + n_iter + 1), best_yields, 'r-', lw=2, label='Best so far')
    ax.axvline(n_init, color='gray', linestyle='--', alpha=0.5, label='BO starts')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Yield (%)')
    ax.set_title('(a) Bayesian Optimization Convergence')
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.semilogy(range(len(acq_values)), acq_values, 'g-', lw=1.5)
    ax.set_xlabel('BO Iteration')
    ax.set_ylabel('Expected Improvement')
    ax.set_title('(b) Acquisition Function Value')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/bayesian_optimization.png', dpi=150, bbox_inches='tight')
    plt.close()

    # --- Figure 6: Parameter exploration ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for idx, (ax, name) in enumerate(zip(axes.flat, param_names)):
        sc = ax.scatter(range(len(X_all)), X_all[:, idx], c=Y_all, cmap='RdYlGn',
                       vmin=60, vmax=95, s=20, edgecolors='k', linewidth=0.3)
        ax.axvline(n_init, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel('Iteration')
        ax.set_ylabel(name)
        ax.set_title(f'Parameter: {name}')
        ax.grid(True, alpha=0.3)
    plt.colorbar(sc, ax=axes, label='Yield (%)', shrink=0.6)
    plt.suptitle('Parameter Exploration During Bayesian Optimization', fontsize=13)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/parameter_exploration.png', dpi=150, bbox_inches='tight')
    plt.close()

    # --- Figure 7: Response surface (T vs flow_rate) ---
    T_grid = np.linspace(60, 120, 50)
    fr_grid = np.linspace(0.1, 1.0, 50)
    TT, FR = np.meshgrid(T_grid, fr_grid)
    # Fix conc and cat at optimum
    opt_conc, opt_cat = X_all[best_idx, 2], X_all[best_idx, 3]

    Z = np.zeros_like(TT)
    for i in range(50):
        for j in range(50):
            x_test = np.array([[TT[i,j], FR[i,j], opt_conc, opt_cat]])
            x_norm = (x_test - lb) / (ub - lb)
            Z[i,j] = gpr.predict(x_norm)[0]

    fig, ax = plt.subplots(figsize=(8, 6))
    c = ax.contourf(TT, FR, Z, levels=20, cmap='RdYlGn')
    ax.scatter(X_all[:, 0], X_all[:, 1], c='black', s=15, alpha=0.6, zorder=5)
    ax.scatter(X_all[best_idx, 0], X_all[best_idx, 1], c='red', s=100, marker='*',
              edgecolors='black', zorder=10, label='Optimum')
    plt.colorbar(c, label='Predicted Yield (%)')
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Flow Rate (mL/min)')
    ax.set_title(f'Response Surface (Conc={opt_conc:.3f} M, Cat={opt_cat:.3f} equiv)')
    ax.legend()
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/response_surface.png', dpi=150, bbox_inches='tight')
    plt.close()

    opt_params = X_all[best_idx]
    opt_yield = Y_all[best_idx]
    print("[BO] Bayesian Optimization complete.")
    print(f"  Best yield: {opt_yield:.1f}%")
    print(f"  Optimal T: {opt_params[0]:.1f}°C")
    print(f"  Optimal flow rate: {opt_params[1]:.3f} mL/min")
    print(f"  Optimal concentration: {opt_params[2]:.3f} M")
    print(f"  Optimal catalyst: {opt_params[3]:.3f} equiv")
    return {
        'best_yield': opt_yield,
        'opt_T': opt_params[0], 'opt_fr': opt_params[1],
        'opt_conc': opt_params[2], 'opt_cat': opt_params[3],
        'n_experiments': len(Y_all),
        'convergence_iter': np.argmax(np.array(best_yields) > opt_yield - 1.0)
    }


# ============================================================
# 4. Online Analytics & Feedback Control
# ============================================================
def simulate_feedback_control():
    """Simulate PID feedback control with online HPLC/IR monitoring."""
    dt = 1.0  # seconds
    t_total = 600  # seconds
    t = np.arange(0, t_total, dt)
    n = len(t)

    # Setpoints
    T_sp = 92.0  # °C
    conc_sp = 0.35  # M (product concentration target)

    # PID parameters for temperature
    Kp_T, Ki_T, Kd_T = 2.0, 0.1, 0.5
    # PID for concentration (via flow rate adjustment)
    Kp_c, Ki_c, Kd_c = 0.5, 0.02, 0.1

    T = np.zeros(n); T[0] = 85.0
    conc = np.zeros(n); conc[0] = 0.28
    T_heater = np.zeros(n); T_heater[0] = 90.0
    flow_adj = np.zeros(n); flow_adj[0] = 0.5

    # HPLC measurements (every 30s)
    hplc_interval = 30
    hplc_times = []
    hplc_readings = []

    # IR measurements (continuous, every 5s)
    ir_interval = 5
    ir_times = []
    ir_readings = []

    integral_T = 0
    integral_c = 0
    prev_err_T = 0
    prev_err_c = 0

    # Disturbances
    dist_T = np.zeros(n)
    dist_T[100:150] = -5  # cooling disturbance
    dist_T[300:350] = 3   # heating disturbance
    dist_c = np.zeros(n)
    dist_c[200:250] = -0.05  # feed concentration drop

    for i in range(1, n):
        # Temperature dynamics (first order + delay)
        tau_T = 15.0  # thermal time constant
        T[i] = T[i-1] + (T_heater[i-1] - T[i-1] + dist_T[i]) * dt / tau_T
        T[i] += np.random.normal(0, 0.2)

        # Concentration dynamics
        tau_c = 25.0
        k_rxn = 0.1 * np.exp(-3000 * (1/(T[i]+273.15) - 1/(T_sp+273.15)))
        conc[i] = conc[i-1] + (conc_sp * flow_adj[i-1] - conc[i-1] * k_rxn + dist_c[i]) * dt / tau_c
        conc[i] += np.random.normal(0, 0.005)

        # PID for temperature
        err_T = T_sp - T[i]
        integral_T += err_T * dt
        integral_T = np.clip(integral_T, -50, 50)
        deriv_T = (err_T - prev_err_T) / dt
        T_heater[i] = T_sp + Kp_T * err_T + Ki_T * integral_T + Kd_T * deriv_T
        T_heater[i] = np.clip(T_heater[i], 60, 130)
        prev_err_T = err_T

        # PID for concentration (adjusted via flow rate)
        err_c = conc_sp - conc[i]
        integral_c += err_c * dt
        integral_c = np.clip(integral_c, -5, 5)
        deriv_c = (err_c - prev_err_c) / dt
        flow_adj[i] = 0.5 + Kp_c * err_c + Ki_c * integral_c + Kd_c * deriv_c
        flow_adj[i] = np.clip(flow_adj[i], 0.1, 1.0)
        prev_err_c = err_c

        # HPLC readings
        if i % hplc_interval == 0:
            hplc_times.append(t[i])
            hplc_readings.append(conc[i] + np.random.normal(0, 0.01))

        # IR readings
        if i % ir_interval == 0:
            ir_times.append(t[i])
            ir_readings.append(conc[i] * 1.0 + np.random.normal(0, 0.015))

    # --- Figure 8: Feedback control ---
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

    ax = axes[0]
    ax.plot(t, T, 'b-', alpha=0.7, label='Reactor Temperature')
    ax.axhline(T_sp, color='r', linestyle='--', label=f'Setpoint ({T_sp}°C)')
    ax.fill_between(t, T_sp-2, T_sp+2, alpha=0.1, color='green', label='±2°C band')
    ax.plot(t, T_heater, 'orange', alpha=0.5, label='Heater Output')
    ax.set_ylabel('Temperature (°C)')
    ax.set_title('(a) Temperature Control with PID')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(t, conc, 'b-', alpha=0.7, label='Product Concentration')
    ax.scatter(hplc_times, hplc_readings, c='red', s=30, zorder=5, label='HPLC Readings')
    ax.scatter(ir_times, ir_readings, c='green', s=10, marker='x', zorder=4, label='IR Readings')
    ax.axhline(conc_sp, color='r', linestyle='--', label=f'Target ({conc_sp} M)')
    ax.set_ylabel('Concentration (M)')
    ax.set_title('(b) Product Concentration Monitoring')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    ax.plot(t, flow_adj, 'g-', alpha=0.7, label='Flow Rate Adjustment')
    ax.set_ylabel('Flow Rate (mL/min)')
    ax.set_xlabel('Time (s)')
    ax.set_title('(c) Flow Rate Control Action')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/feedback_control.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Compute control performance
    ss_start = 200  # after settling
    T_ss = T[ss_start:]
    conc_ss = conc[ss_start:]
    T_rmse = np.sqrt(np.mean((T_ss - T_sp)**2))
    conc_rmse = np.sqrt(np.mean((conc_ss - conc_sp)**2))

    print("[Control] Feedback control simulation complete.")
    print(f"  Temperature RMSE: {T_rmse:.2f}°C")
    print(f"  Concentration RMSE: {conc_rmse:.4f} M")
    return {'T_rmse': T_rmse, 'conc_rmse': conc_rmse}


# ============================================================
# 5. Scale-up Design: Numbering Up vs Scaling Up
# ============================================================
def scaleup_analysis():
    """Compare numbering-up and scaling-up strategies."""
    # Single reactor parameters
    V_single = 0.5e-6  # 0.5 mL
    Q_single = 0.5e-6 / 60  # 0.5 mL/min -> m³/s
    yield_single = 92.0
    d_h = 0.001  # 1 mm hydraulic diameter

    # Target throughputs (g/h)
    targets = np.array([1, 5, 10, 50, 100, 500])

    # Numbering up
    n_reactors = np.ceil(targets / 1.0)  # 1 g/h per reactor
    nu_yield = yield_single * np.ones_like(targets, dtype=float)
    nu_yield -= np.random.uniform(0, 1.5, len(targets))  # slight distribution loss
    nu_cost_relative = n_reactors * 1.0 + n_reactors * 0.1  # linear + distribution cost

    # Scaling up (increasing channel diameter)
    scale_factor = targets / 1.0
    su_diameter = d_h * scale_factor**(1/3) * 1000  # mm
    # Reynolds number effect on mixing
    Re_base = 50
    Re_scaled = Re_base * scale_factor**(1/3)
    # Yield drops due to poorer mixing at larger scales
    su_yield = yield_single - 2 * np.log(scale_factor) - 0.05 * (Re_scaled - Re_base)
    su_yield = np.maximum(su_yield, 50)
    su_cost_relative = scale_factor**0.7  # economies of scale

    # Hybrid approach
    hybrid_n = np.ceil(np.sqrt(n_reactors))
    hybrid_scale = targets / hybrid_n
    hybrid_yield = yield_single - 0.5 * np.log(hybrid_scale) - 0.5
    hybrid_cost = hybrid_n * 1.5 + hybrid_scale**0.5

    # --- Figure 9: Scale-up comparison ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    ax = axes[0]
    ax.semilogx(targets, nu_yield, 'bo-', label='Numbering Up', lw=2)
    ax.semilogx(targets, su_yield, 'rs-', label='Scaling Up', lw=2)
    ax.semilogx(targets, hybrid_yield, 'g^-', label='Hybrid', lw=2)
    ax.set_xlabel('Target Throughput (g/h)')
    ax.set_ylabel('Yield (%)')
    ax.set_title('(a) Yield vs Throughput')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(40, 100)

    ax = axes[1]
    ax.loglog(targets, nu_cost_relative, 'bo-', label='Numbering Up', lw=2)
    ax.loglog(targets, su_cost_relative, 'rs-', label='Scaling Up', lw=2)
    ax.loglog(targets, hybrid_cost, 'g^-', label='Hybrid', lw=2)
    ax.set_xlabel('Target Throughput (g/h)')
    ax.set_ylabel('Relative Cost')
    ax.set_title('(b) Cost vs Throughput')
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    efficiency_nu = nu_yield / (nu_cost_relative / nu_cost_relative[0])
    efficiency_su = su_yield / (su_cost_relative / su_cost_relative[0])
    efficiency_hy = hybrid_yield / (hybrid_cost / hybrid_cost[0])
    ax.semilogx(targets, efficiency_nu, 'bo-', label='Numbering Up', lw=2)
    ax.semilogx(targets, efficiency_su, 'rs-', label='Scaling Up', lw=2)
    ax.semilogx(targets, efficiency_hy, 'g^-', label='Hybrid', lw=2)
    ax.set_xlabel('Target Throughput (g/h)')
    ax.set_ylabel('Cost-Efficiency (Yield/Relative Cost)')
    ax.set_title('(c) Cost-Efficiency Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/scaleup_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()

    print("[Scale-up] Analysis complete.")
    return {
        'nu_yield_500': nu_yield[-1], 'su_yield_500': su_yield[-1],
        'hybrid_yield_500': hybrid_yield[-1]
    }


# ============================================================
# 6. Pharmaceutical Case Study: Ibuprofen Intermediate
# ============================================================
def pharma_case_study():
    """Continuous flow synthesis of an ibuprofen intermediate via Friedel-Crafts acylation."""
    # Simulate optimization trajectory
    n_exp = 35

    # Grid of conditions explored
    temps = np.random.uniform(40, 120, n_exp)
    res_times = np.random.uniform(5, 120, n_exp)
    cat_loadings = np.random.uniform(1.0, 3.0, n_exp)

    # Yield model for Friedel-Crafts acylation
    def fc_yield(T, rt, cat):
        y = (60 + 25 * np.exp(-0.5 * ((T-85)/15)**2)
             * (1 - np.exp(-rt/30))
             * np.minimum(cat/1.5, 1.0))
        return y + np.random.normal(0, 2)

    yields = np.array([fc_yield(T, rt, cat) for T, rt, cat in zip(temps, res_times, cat_loadings)])

    # Selectivity model
    def fc_selectivity(T, rt):
        return 95 - 0.15 * (T - 80)**2 / 100 - 0.1 * rt + np.random.normal(0, 1.5)

    selectivities = np.array([fc_selectivity(T, rt) for T, rt in zip(temps, res_times)])

    # Purity via simulated HPLC
    purities = 0.6 * selectivities / 100 + 0.35 + np.random.normal(0, 0.02, n_exp)
    purities = np.clip(purities, 0.5, 0.999)

    # --- Figure 10: Case study results ---
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    ax = axes[0, 0]
    sc = ax.scatter(temps, res_times, c=yields, cmap='RdYlGn', s=50, edgecolors='k', linewidth=0.5)
    plt.colorbar(sc, ax=ax, label='Yield (%)')
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Residence Time (s)')
    ax.set_title('(a) Yield Map')

    ax = axes[0, 1]
    sc = ax.scatter(temps, res_times, c=selectivities, cmap='RdYlBu', s=50, edgecolors='k', linewidth=0.5)
    plt.colorbar(sc, ax=ax, label='Selectivity (%)')
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Residence Time (s)')
    ax.set_title('(b) Selectivity Map')

    ax = axes[1, 0]
    # Pareto front: yield vs selectivity
    ax.scatter(yields, selectivities, c=temps, cmap='coolwarm', s=50, edgecolors='k', linewidth=0.5)
    # Find Pareto front
    pareto = []
    for i in range(len(yields)):
        dominated = False
        for j in range(len(yields)):
            if yields[j] > yields[i] and selectivities[j] > selectivities[i]:
                dominated = True
                break
        if not dominated:
            pareto.append(i)
    pareto_y = yields[pareto]
    pareto_s = selectivities[pareto]
    order = np.argsort(pareto_y)
    ax.plot(pareto_y[order], pareto_s[order], 'r--', lw=2, label='Pareto Front')
    c2 = ax.scatter(yields, selectivities, c=temps, cmap='coolwarm', s=50, edgecolors='k', linewidth=0.5)
    plt.colorbar(c2, ax=ax, label='Temperature (°C)')
    ax.set_xlabel('Yield (%)')
    ax.set_ylabel('Selectivity (%)')
    ax.set_title('(c) Yield-Selectivity Trade-off')
    ax.legend()

    ax = axes[1, 1]
    # Production rate over optimization
    prod_rates = yields * 0.15 * 60 / 1000  # g/h estimate
    sorted_idx = np.argsort(yields)[::-1]
    cummax_prod = np.maximum.accumulate(prod_rates)
    ax.plot(range(n_exp), cummax_prod, 'b-', lw=2, label='Best Production Rate')
    ax.bar(range(n_exp), prod_rates, alpha=0.4, color='steelblue', label='Individual Runs')
    ax.set_xlabel('Experiment Number')
    ax.set_ylabel('Production Rate (g/h)')
    ax.set_title('(d) Production Rate Optimization')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.suptitle('Pharmaceutical Case Study: Ibuprofen Intermediate (Friedel-Crafts Acylation)', fontsize=13)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/pharma_case_study.png', dpi=150, bbox_inches='tight')
    plt.close()

    best_idx = np.argmax(yields)
    print("[Pharma] Case study complete.")
    print(f"  Best yield: {yields[best_idx]:.1f}%")
    print(f"  Best selectivity: {selectivities[best_idx]:.1f}%")
    print(f"  Optimal T: {temps[best_idx]:.1f}°C, τ: {res_times[best_idx]:.1f}s, cat: {cat_loadings[best_idx]:.2f} equiv")
    return {
        'best_yield': yields[best_idx],
        'best_selectivity': selectivities[best_idx],
        'best_T': temps[best_idx],
        'best_rt': res_times[best_idx],
        'best_cat': cat_loadings[best_idx],
        'purity': purities[best_idx]
    }


# ============================================================
# 7. System Architecture Diagram
# ============================================================
def create_system_diagram():
    """Create system architecture diagram for process control integration."""
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8)
    ax.set_aspect('equal')
    ax.axis('off')

    # Helper functions
    def draw_box(x, y, w, h, text, color='lightblue', fontsize=8):
        rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor='black', lw=1.5, zorder=2)
        ax.add_patch(rect)
        ax.text(x+w/2, y+h/2, text, ha='center', va='center', fontsize=fontsize,
               fontweight='bold', wrap=True, zorder=3)

    def draw_arrow(x1, y1, x2, y2, text='', color='black'):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color=color, lw=1.5), zorder=1)
        if text:
            mx, my = (x1+x2)/2, (y1+y2)/2
            ax.text(mx, my+0.15, text, ha='center', va='bottom', fontsize=6, color=color)

    # Process layer
    draw_box(0.5, 4.5, 2.5, 1.2, 'Feed Pumps\n& Mixers', '#FFE0B2')
    draw_box(3.5, 4.5, 2.5, 1.2, 'Microreactor\n(Serpentine)', '#C8E6C9')
    draw_box(6.5, 4.5, 2.5, 1.2, 'Heat\nExchanger', '#FFCDD2')
    draw_box(9.5, 4.5, 2.5, 1.2, 'Product\nCollection', '#E1BEE7')

    # Analytics layer
    draw_box(3.5, 2.5, 2, 1.0, 'Online HPLC', '#B3E5FC')
    draw_box(6.5, 2.5, 2, 1.0, 'Inline FTIR', '#B3E5FC')
    draw_box(9.5, 2.5, 2, 1.0, 'Temperature\nSensors', '#B3E5FC')

    # Control layer
    draw_box(1, 0.5, 3, 1.2, 'Bayesian Optimization\nEngine (Python/GPyOpt)', '#FFF9C4', fontsize=7)
    draw_box(5, 0.5, 3, 1.2, 'PID Controller\n(LabVIEW/OPC-UA)', '#FFF9C4', fontsize=7)
    draw_box(9, 0.5, 3.5, 1.2, 'Data Acquisition\n& SCADA System', '#FFF9C4', fontsize=7)

    # Top layer - supervisory
    draw_box(3, 6.5, 4, 1, 'Process Control\nSoftware (MES/ERP)', '#D1C4E9', fontsize=8)
    draw_box(8, 6.5, 4, 1, 'Digital Twin\n& Simulation', '#D1C4E9', fontsize=8)

    # Arrows - process flow
    draw_arrow(3.0, 5.1, 3.5, 5.1, 'Flow', 'darkgreen')
    draw_arrow(6.0, 5.1, 6.5, 5.1, 'Flow', 'darkgreen')
    draw_arrow(9.0, 5.1, 9.5, 5.1, 'Flow', 'darkgreen')

    # Arrows - analytics
    draw_arrow(4.5, 4.5, 4.5, 3.5, 'Sample', 'blue')
    draw_arrow(7.5, 4.5, 7.5, 3.5, 'Signal', 'blue')
    draw_arrow(10.5, 4.5, 10.5, 3.5, 'T data', 'blue')

    # Arrows - control
    draw_arrow(4.5, 2.5, 2.5, 1.7, 'Yield data', 'red')
    draw_arrow(7.5, 2.5, 6.5, 1.7, 'Spectral data', 'red')
    draw_arrow(10.5, 2.5, 10.5, 1.7, 'Sensor data', 'red')

    # Arrows - actuation
    draw_arrow(2.5, 1.7, 1.75, 4.5, 'Set flow', 'purple')
    draw_arrow(6.5, 1.7, 7.75, 4.5, 'Set T', 'purple')

    # Arrows - supervisory
    draw_arrow(5.0, 6.5, 4.0, 5.7, '', 'gray')
    draw_arrow(10.0, 6.5, 10.0, 5.7, '', 'gray')
    draw_arrow(7.0, 7.0, 8.0, 7.0, 'Model update', 'gray')

    ax.set_title('Integrated Continuous Flow Synthesis: Process Control Architecture', fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/system_architecture.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("[Diagram] System architecture created.")


# ============================================================
# Main execution
# ============================================================
if __name__ == '__main__':
    print("=" * 60)
    print("Continuous Flow Synthesis Automated Optimization System")
    print("=" * 60)

    results = {}

    print("\n--- Step 1: CFD Simulation ---")
    results['cfd'] = simulate_cfd()

    print("\n--- Step 2: RTD Analysis ---")
    results['rtd'] = simulate_rtd()

    print("\n--- Step 3: Bayesian Optimization ---")
    results['bo'] = bayesian_optimization()

    print("\n--- Step 4: Feedback Control ---")
    results['control'] = simulate_feedback_control()

    print("\n--- Step 5: Scale-up Analysis ---")
    results['scaleup'] = scaleup_analysis()

    print("\n--- Step 6: Pharmaceutical Case Study ---")
    results['pharma'] = pharma_case_study()

    print("\n--- Step 7: System Architecture ---")
    create_system_diagram()

    # Save results
    # Convert numpy types to Python types for JSON serialization
    def convert(obj):
        if isinstance(obj, (np.integer,)): return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        return obj

    serializable = {}
    for k, v in results.items():
        serializable[k] = {kk: convert(vv) for kk, vv in v.items()}

    with open('results.json', 'w') as f:
        json.dump(serializable, f, indent=2)

    print("\n" + "=" * 60)
    print("All simulations complete. Results saved.")
    print(f"Figures: {os.listdir(FIGDIR)}")
    print("=" * 60)
