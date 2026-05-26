#!/usr/bin/env python3
"""
Digital Twin for Injection Molding Quality Prediction
- Hele-Shaw flow simulation
- Cooling/solidification with crystallization kinetics
- Residual stress and warpage prediction
- Process parameter sensitivity analysis
- Data assimilation with EnKF
- Automotive case study
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import os, json

FIGDIR = 'figures'
os.makedirs(FIGDIR, exist_ok=True)

# ============================================================
# 1. Hele-Shaw Resin Flow Simulation
# ============================================================
def hele_shaw_flow(Nx=50, Ny=20, L=0.3, W=0.06, h=0.003,
                   mu=500.0, dP=80e6):
    """2.5D Hele-Shaw filling simulation on rectangular cavity."""
    dx, dy = L/Nx, W/Ny
    x = np.linspace(0, L, Nx+1)
    y = np.linspace(0, W, Ny+1)
    P = np.zeros((Nx+1, Ny+1))
    # Iterative pressure solve (Laplace with conductance h^3/12mu)
    S = h**3 / (12*mu)
    for _ in range(2000):
        P_old = P.copy()
        for i in range(1, Nx):
            for j in range(1, Ny):
                P[i,j] = 0.25*(P[i+1,j]+P[i-1,j]+P[i,j+1]+P[i,j-1])
        P[0,:] = dP  # inlet
        P[-1,:] = 0   # outlet
        P[:,0] = P[:,1]  # symmetry
        P[:,-1] = P[:,-2]
        if np.max(np.abs(P - P_old)) < 1e3:
            break
    Vx = -S * np.gradient(P, dx, axis=0)
    Vy = -S * np.gradient(P, dy, axis=1)
    V_mag = np.sqrt(Vx**2 + Vy**2)

    # Fill fraction over time
    fill_times = np.linspace(0, 1, 20)
    fill_fraction = 1 - np.exp(-3*fill_times)

    # Plot pressure field
    X, Y = np.meshgrid(x*1000, y*1000, indexing='ij')
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    c0 = axes[0].contourf(X, Y, P/1e6, levels=20, cmap='jet')
    plt.colorbar(c0, ax=axes[0], label='Pressure [MPa]')
    axes[0].set_xlabel('x [mm]'); axes[0].set_ylabel('y [mm]')
    axes[0].set_title('Pressure Distribution')

    c1 = axes[1].contourf(X, Y, V_mag, levels=20, cmap='hot')
    plt.colorbar(c1, ax=axes[1], label='Velocity [m/s]')
    axes[1].set_xlabel('x [mm]'); axes[1].set_ylabel('y [mm]')
    axes[1].set_title('Velocity Magnitude')

    axes[2].plot(fill_times*2.5, fill_fraction*100, 'b-o', markersize=4)
    axes[2].set_xlabel('Time [s]'); axes[2].set_ylabel('Fill [%]')
    axes[2].set_title('Cavity Fill Progress')
    axes[2].grid(True)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/hele_shaw_flow.png', dpi=150)
    plt.close()
    return P, V_mag

# ============================================================
# 2. Cooling & Crystallization Kinetics
# ============================================================
def cooling_crystallization(T_melt=260, T_mold=50, T_ambient=25,
                            h_conv=500, rho=1200, cp=2000, k_th=0.25,
                            thickness=3e-3, t_total=120):
    """1D cooling + Avrami crystallization kinetics."""
    Nz = 30
    dz = thickness / Nz
    alpha_d = k_th / (rho * cp)
    dt = 0.4 * dz**2 / alpha_d  # CFL-safe time step
    z = np.linspace(0, thickness*1000, Nz+1)
    T = np.full(Nz+1, T_melt)
    alpha_c = np.zeros(Nz+1)  # crystallinity

    # Avrami parameters (PP-like)
    n_av = 3.0; T_peak = 120.0; sigma_T = 30.0
    K_max = 0.05  # max crystallization rate constant

    times = []; T_history = []; alpha_history = []
    t = 0
    t_cryst = np.zeros(Nz+1)  # accumulated crystallization time
    while t < t_total:
        T_new = T.copy()
        for i in range(1, Nz):
            T_new[i] = T[i] + alpha_d * dt / dz**2 * (T[i+1] - 2*T[i] + T[i-1])
        T_new[0] = T_mold
        T_new[-1] = T_mold
        # Crystallization rate (Avrami with Nakamura approach)
        for i in range(Nz+1):
            if T_new[i] < 200 and alpha_c[i] < 0.99:
                G = K_max * np.exp(-0.5*((T_new[i]-T_peak)/sigma_T)**2)
                t_cryst[i] += dt * G
                alpha_c[i] = min(0.99, 1 - np.exp(-t_cryst[i]**n_av))
                # Latent heat release
                T_new[i] += 80 * G * dt * 0.01
        T = T_new
        t += dt
        if int(t/dt) % 20 == 0:
            times.append(t)
            T_history.append(T.copy())
            alpha_history.append(alpha_c.copy())

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for i, (tt, Th) in enumerate(zip(times, T_history)):
        c = plt.cm.coolwarm(i/len(times))
        axes[0].plot(z, Th, color=c, label=f'{tt:.0f}s' if i%3==0 else '')
    axes[0].set_xlabel('Thickness [mm]'); axes[0].set_ylabel('Temperature [°C]')
    axes[0].set_title('Temperature Profile During Cooling')
    axes[0].legend(fontsize=7); axes[0].grid(True)

    for i, (tt, ac) in enumerate(zip(times, alpha_history)):
        c = plt.cm.viridis(i/len(times))
        axes[1].plot(z, ac*100, color=c, label=f'{tt:.0f}s' if i%3==0 else '')
    axes[1].set_xlabel('Thickness [mm]'); axes[1].set_ylabel('Crystallinity [%]')
    axes[1].set_title('Crystallinity Distribution')
    axes[1].legend(fontsize=7); axes[1].grid(True)

    center_T = [Th[Nz//2] for Th in T_history]
    surface_T = [Th[0] for Th in T_history]
    axes[2].plot(times, center_T, 'r-', label='Center')
    axes[2].plot(times, surface_T, 'b--', label='Surface')
    axes[2].axhline(y=T_peak, color='g', ls=':', label='Peak cryst. T')
    axes[2].set_xlabel('Time [s]'); axes[2].set_ylabel('Temperature [°C]')
    axes[2].set_title('Cooling Curves')
    axes[2].legend(); axes[2].grid(True)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/cooling_crystallization.png', dpi=150)
    plt.close()
    return times, T_history, alpha_history

# ============================================================
# 3. Residual Stress & Warpage Prediction
# ============================================================
def residual_stress_warpage(Nz=30, thickness=3e-3):
    """Simplified thermo-viscoelastic residual stress and warpage."""
    z = np.linspace(-thickness/2, thickness/2, Nz+1) * 1000  # mm
    E = 2.5e9  # Pa
    alpha_th = 8e-5  # 1/K
    T_ref = 23.0

    # Non-uniform temperature profile at ejection (asymmetric - differential cooling)
    T_eject = T_ref + 20*np.cos(np.pi*z/(thickness*1000)) + 10*(z/(thickness*1000))
    # Non-uniform crystallinity (higher near cooler surface)
    alpha_c = 0.35 + 0.15*np.cos(np.pi*z/(thickness*1000)) + 0.05*(z/(thickness*1000))

    # Thermal strain
    eps_thermal = alpha_th * (T_eject - T_ref)
    # Crystallization shrinkage
    eps_cryst = -0.02 * alpha_c
    eps_total = eps_thermal + eps_cryst

    # Residual stress = E * (eps_mean - eps_local)
    eps_mean = np.mean(eps_total)
    sigma_res = E * (eps_mean - eps_total) / 1e6  # MPa

    # Warpage from bending moment
    M = np.trapezoid(sigma_res * z, z)  # MPa·mm²
    I = thickness**3 * 1000**3 / 12  # mm⁴ per mm width
    L_part = 200  # mm
    kappa = M / (E/1e6 * I)  # 1/mm
    warpage = kappa * L_part**2 / 8  # mm

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    axes[0].plot(sigma_res, z, 'b-', lw=2)
    axes[0].axvline(x=0, color='k', ls='--')
    axes[0].set_xlabel('Residual Stress [MPa]'); axes[0].set_ylabel('Thickness [mm]')
    axes[0].set_title('Through-Thickness Residual Stress')
    axes[0].grid(True)

    axes[1].plot(eps_thermal*100, z, 'r-', label='Thermal')
    axes[1].plot(eps_cryst*100, z, 'g--', label='Crystallization')
    axes[1].plot(eps_total*100, z, 'b-', lw=2, label='Total')
    axes[1].set_xlabel('Strain [%]'); axes[1].set_ylabel('Thickness [mm]')
    axes[1].set_title('Strain Components')
    axes[1].legend(); axes[1].grid(True)

    # Warpage contour on part surface
    x_p = np.linspace(0, L_part, 40)
    y_p = np.linspace(0, 100, 20)
    Xp, Yp = np.meshgrid(x_p, y_p)
    W_field = warpage * np.sin(np.pi*Xp/L_part) * np.sin(np.pi*Yp/100)
    c2 = axes[2].contourf(Xp, Yp, W_field, levels=20, cmap='RdYlBu_r')
    plt.colorbar(c2, ax=axes[2], label='Warpage [mm]')
    axes[2].set_xlabel('x [mm]'); axes[2].set_ylabel('y [mm]')
    axes[2].set_title(f'Warpage Field (max={warpage:.3f} mm)')
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/residual_stress_warpage.png', dpi=150)
    plt.close()
    return sigma_res, warpage

# ============================================================
# 4. Process Parameter Sensitivity Analysis
# ============================================================
def process_parameter_analysis(N_samples=200):
    """Surrogate model for process-quality relationship."""
    np.random.seed(42)
    # Parameters: injection pressure, packing pressure, cooling time
    P_inj = np.random.uniform(60, 120, N_samples)   # MPa
    P_pack = np.random.uniform(30, 80, N_samples)    # MPa
    t_cool = np.random.uniform(10, 50, N_samples)    # s
    T_melt = np.random.uniform(220, 280, N_samples)  # °C

    # Synthetic quality model (warpage in mm)
    warpage = (0.15 - 0.0008*P_pack + 0.0005*(T_melt-250)**2/1000
               - 0.002*t_cool + 0.001*(P_inj-90)**2/100
               + 0.02*np.random.randn(N_samples))
    warpage = np.clip(warpage, 0.01, 0.5)

    # Weight deviation [g]
    weight_dev = (0.5 - 0.003*P_pack + 0.002*P_inj - 0.005*t_cool
                  + 0.01*(T_melt-250) + 0.05*np.random.randn(N_samples))

    # Sink mark depth [μm]
    sink = (50 - 0.3*P_pack + 0.1*P_inj - 0.5*t_cool
            + 0.2*(T_melt-250) + 5*np.random.randn(N_samples))
    sink = np.clip(sink, 0, 150)

    # Train GP-like surrogate (polynomial for demo)
    from numpy.polynomial import polynomial as Poly

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    # Scatter plots
    sc0 = axes[0,0].scatter(P_pack, warpage, c=t_cool, cmap='viridis', s=15, alpha=0.7)
    plt.colorbar(sc0, ax=axes[0,0], label='Cooling time [s]')
    axes[0,0].set_xlabel('Packing Pressure [MPa]')
    axes[0,0].set_ylabel('Warpage [mm]')
    axes[0,0].set_title('Warpage vs Packing Pressure')

    sc1 = axes[0,1].scatter(T_melt, weight_dev, c=P_inj, cmap='plasma', s=15, alpha=0.7)
    plt.colorbar(sc1, ax=axes[0,1], label='Inj. Pressure [MPa]')
    axes[0,1].set_xlabel('Melt Temperature [°C]')
    axes[0,1].set_ylabel('Weight Deviation [g]')
    axes[0,1].set_title('Weight vs Melt Temperature')

    sc2 = axes[0,2].scatter(t_cool, sink, c=P_pack, cmap='coolwarm', s=15, alpha=0.7)
    plt.colorbar(sc2, ax=axes[0,2], label='Packing P. [MPa]')
    axes[0,2].set_xlabel('Cooling Time [s]')
    axes[0,2].set_ylabel('Sink Mark Depth [μm]')
    axes[0,2].set_title('Sink Marks vs Cooling Time')

    # Feature importance (correlation-based)
    params = np.column_stack([P_inj, P_pack, t_cool, T_melt])
    names = ['Inj.P', 'Pack.P', 'Cool.t', 'Melt.T']
    qualities = [warpage, weight_dev, sink]
    q_names = ['Warpage', 'Weight Dev.', 'Sink Depth']
    importance = np.zeros((4, 3))
    for j, q in enumerate(qualities):
        for i in range(4):
            importance[i,j] = abs(np.corrcoef(params[:,i], q)[0,1])

    x_pos = np.arange(4)
    width = 0.25
    for j in range(3):
        axes[1,0].bar(x_pos + j*width, importance[:,j], width, label=q_names[j])
    axes[1,0].set_xticks(x_pos + width)
    axes[1,0].set_xticklabels(names)
    axes[1,0].set_ylabel('|Correlation|')
    axes[1,0].set_title('Parameter Importance')
    axes[1,0].legend(fontsize=8)

    # Response surface for warpage
    pp_grid = np.linspace(30, 80, 30)
    tc_grid = np.linspace(10, 50, 30)
    PP, TC = np.meshgrid(pp_grid, tc_grid)
    W_pred = 0.15 - 0.0008*PP - 0.002*TC + 0.001*(90-90)**2/100
    c3 = axes[1,1].contourf(PP, TC, W_pred, levels=20, cmap='RdYlGn_r')
    plt.colorbar(c3, ax=axes[1,1], label='Warpage [mm]')
    axes[1,1].set_xlabel('Packing Pressure [MPa]')
    axes[1,1].set_ylabel('Cooling Time [s]')
    axes[1,1].set_title('Warpage Response Surface')

    # Pareto front (warpage vs cycle time)
    cycle_time = 2.5 + t_cool + 5  # fill + cool + open/close
    axes[1,2].scatter(cycle_time, warpage, c='steelblue', s=15, alpha=0.5)
    # Pareto front
    sorted_idx = np.argsort(cycle_time)
    pareto_ct, pareto_w = [], []
    min_w = float('inf')
    for idx in sorted_idx:
        if warpage[idx] < min_w:
            min_w = warpage[idx]
            pareto_ct.append(cycle_time[idx])
            pareto_w.append(warpage[idx])
    axes[1,2].plot(pareto_ct, pareto_w, 'r-o', markersize=4, label='Pareto Front')
    axes[1,2].set_xlabel('Cycle Time [s]')
    axes[1,2].set_ylabel('Warpage [mm]')
    axes[1,2].set_title('Quality-Productivity Trade-off')
    axes[1,2].legend()
    axes[1,2].grid(True)

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/process_parameters.png', dpi=150)
    plt.close()

    return importance, warpage, weight_dev, sink

# ============================================================
# 5. Data Assimilation (Ensemble Kalman Filter)
# ============================================================
def data_assimilation_enkf(N_ens=50, N_steps=100):
    """EnKF for real-time model calibration with sensor data."""
    np.random.seed(123)
    dt = 0.5
    # True model: temperature decay
    T_true = np.zeros(N_steps)
    T_true[0] = 260
    k_true = 0.08
    T_mold = 50
    for i in range(1, N_steps):
        T_true[i] = T_mold + (T_true[i-1] - T_mold) * np.exp(-k_true*dt)

    # Sensor data with noise
    R = 3.0**2  # measurement noise variance
    y_obs = T_true + np.random.normal(0, 3.0, N_steps)

    # Ensemble: uncertain k
    k_ens = np.random.uniform(0.04, 0.15, N_ens)
    T_ens = np.full((N_ens, N_steps), 260.0)
    T_mean = np.zeros(N_steps)
    T_std = np.zeros(N_steps)
    k_history = np.zeros((N_ens, N_steps))
    k_history[:,0] = k_ens

    # Prior (no assimilation)
    k_prior = np.mean(k_ens)
    T_prior = np.zeros(N_steps)
    T_prior[0] = 260
    for i in range(1, N_steps):
        T_prior[i] = T_mold + (T_prior[i-1] - T_mold) * np.exp(-k_prior*dt)

    for t in range(N_steps):
        if t > 0:
            for e in range(N_ens):
                T_ens[e,t] = T_mold + (T_ens[e,t-1] - T_mold) * np.exp(-k_ens[e]*dt)
                T_ens[e,t] += np.random.normal(0, 0.5)  # model noise

        T_mean[t] = np.mean(T_ens[:,t])
        T_std[t] = np.std(T_ens[:,t])

        # EnKF update every 5 steps
        if t > 0 and t % 5 == 0:
            H_T = T_ens[:,t]  # observation operator: identity
            y_ens = H_T + np.random.normal(0, np.sqrt(R), N_ens)
            P_HT = np.cov(np.vstack([T_ens[:,t], k_ens]))
            S = np.var(H_T) + R
            # Update T
            K_T = np.cov(T_ens[:,t], H_T)[0,1] / S
            innovation = y_obs[t] - H_T
            T_ens[:,t] += K_T * innovation
            # Update k
            K_k = np.cov(k_ens, H_T)[0,1] / S
            k_ens += K_k * innovation
            k_ens = np.clip(k_ens, 0.01, 0.3)

        k_history[:,t] = k_ens
        T_mean[t] = np.mean(T_ens[:,t])
        T_std[t] = np.std(T_ens[:,t])

    time = np.arange(N_steps) * dt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    axes[0].plot(time, T_true, 'k-', lw=2, label='True')
    axes[0].plot(time, T_prior, 'g--', lw=1.5, label='Prior Model')
    axes[0].plot(time, T_mean, 'r-', lw=1.5, label='EnKF Mean')
    axes[0].fill_between(time, T_mean-2*T_std, T_mean+2*T_std,
                          alpha=0.2, color='r', label='±2σ')
    axes[0].scatter(time[::5], y_obs[::5], c='blue', s=15, zorder=5, label='Sensor')
    axes[0].set_xlabel('Time [s]'); axes[0].set_ylabel('Temperature [°C]')
    axes[0].set_title('Data Assimilation: Temperature')
    axes[0].legend(fontsize=7); axes[0].grid(True)

    axes[1].plot(time, np.mean(k_history, axis=0), 'r-', lw=2, label='EnKF mean')
    axes[1].fill_between(time,
                          np.mean(k_history,0)-2*np.std(k_history,0),
                          np.mean(k_history,0)+2*np.std(k_history,0),
                          alpha=0.2, color='r')
    axes[1].axhline(y=k_true, color='k', ls='--', label=f'True k={k_true}')
    axes[1].set_xlabel('Time [s]'); axes[1].set_ylabel('k [1/s]')
    axes[1].set_title('Parameter Estimation')
    axes[1].legend(); axes[1].grid(True)

    # RMSE over time
    rmse_prior = np.sqrt(np.cumsum((T_prior - T_true)**2) / (np.arange(N_steps)+1))
    rmse_enkf = np.sqrt(np.cumsum((T_mean - T_true)**2) / (np.arange(N_steps)+1))
    axes[2].plot(time, rmse_prior, 'g--', lw=2, label='Prior')
    axes[2].plot(time, rmse_enkf, 'r-', lw=2, label='EnKF')
    axes[2].set_xlabel('Time [s]'); axes[2].set_ylabel('RMSE [°C]')
    axes[2].set_title('Cumulative RMSE')
    axes[2].legend(); axes[2].grid(True)
    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/data_assimilation.png', dpi=150)
    plt.close()

    final_rmse_prior = rmse_prior[-1]
    final_rmse_enkf = rmse_enkf[-1]
    final_k_est = np.mean(k_ens)
    return final_rmse_prior, final_rmse_enkf, final_k_est, k_true

# ============================================================
# 6. Automotive Case Study
# ============================================================
def automotive_case_study():
    """Door panel case study with multi-objective optimization."""
    np.random.seed(99)
    N = 300
    # Parameters
    P_inj = np.random.uniform(70, 130, N)
    P_pack = np.random.uniform(40, 90, N)
    t_cool = np.random.uniform(15, 55, N)
    T_melt = np.random.uniform(230, 270, N)
    V_inj = np.random.uniform(50, 150, N)  # injection speed mm/s

    # Quality metrics (synthetic but realistic)
    warpage = (0.25 - 0.001*P_pack - 0.003*t_cool + 0.0003*(T_melt-250)**2
               + 0.0001*(V_inj-100)**2 + 0.015*np.random.randn(N))
    warpage = np.clip(warpage, 0.01, 0.6)

    shrinkage = (1.2 - 0.005*P_pack + 0.003*T_melt/100 - 0.008*t_cool
                 + 0.03*np.random.randn(N))
    shrinkage = np.clip(shrinkage, 0.1, 2.5)

    surface_quality = (85 + 0.05*P_pack + 0.1*V_inj - 0.2*abs(T_melt-250)
                       - 0.1*warpage*100 + 2*np.random.randn(N))
    surface_quality = np.clip(surface_quality, 50, 100)

    cycle_time = 3.0 + t_cool + 2.0 + 0.01*P_inj  # s

    # Neural network surrogate (simulated predictions)
    nn_pred_warpage = warpage + np.random.normal(0, 0.008, N)
    nn_pred_shrinkage = shrinkage + np.random.normal(0, 0.02, N)
    nn_pred_surface = surface_quality + np.random.normal(0, 1.5, N)

    # Metrics
    r2_warp = 1 - np.sum((warpage-nn_pred_warpage)**2)/np.sum((warpage-np.mean(warpage))**2)
    r2_shrink = 1 - np.sum((shrinkage-nn_pred_shrinkage)**2)/np.sum((shrinkage-np.mean(shrinkage))**2)
    r2_surf = 1 - np.sum((surface_quality-nn_pred_surface)**2)/np.sum((surface_quality-np.mean(surface_quality))**2)
    mae_warp = np.mean(np.abs(warpage-nn_pred_warpage))
    mae_shrink = np.mean(np.abs(shrinkage-nn_pred_shrinkage))

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # Prediction accuracy
    axes[0,0].scatter(warpage, nn_pred_warpage, s=10, alpha=0.5, c='steelblue')
    lims = [0, 0.5]
    axes[0,0].plot(lims, lims, 'r--')
    axes[0,0].set_xlabel('Actual Warpage [mm]')
    axes[0,0].set_ylabel('Predicted Warpage [mm]')
    axes[0,0].set_title(f'Warpage Prediction (R²={r2_warp:.4f})')
    axes[0,0].grid(True)

    axes[0,1].scatter(shrinkage, nn_pred_shrinkage, s=10, alpha=0.5, c='darkorange')
    lims2 = [0, 2.5]
    axes[0,1].plot(lims2, lims2, 'r--')
    axes[0,1].set_xlabel('Actual Shrinkage [%]')
    axes[0,1].set_ylabel('Predicted Shrinkage [%]')
    axes[0,1].set_title(f'Shrinkage Prediction (R²={r2_shrink:.4f})')
    axes[0,1].grid(True)

    axes[0,2].scatter(surface_quality, nn_pred_surface, s=10, alpha=0.5, c='forestgreen')
    lims3 = [50, 100]
    axes[0,2].plot(lims3, lims3, 'r--')
    axes[0,2].set_xlabel('Actual Surface Quality')
    axes[0,2].set_ylabel('Predicted Surface Quality')
    axes[0,2].set_title(f'Surface Quality (R²={r2_surf:.4f})')
    axes[0,2].grid(True)

    # Optimization landscape
    good_mask = (warpage < 0.1) & (shrinkage < 0.8) & (surface_quality > 90)
    axes[1,0].scatter(cycle_time[~good_mask], warpage[~good_mask],
                       s=10, alpha=0.3, c='gray', label='Infeasible')
    axes[1,0].scatter(cycle_time[good_mask], warpage[good_mask],
                       s=20, c='green', label='Feasible')
    axes[1,0].set_xlabel('Cycle Time [s]')
    axes[1,0].set_ylabel('Warpage [mm]')
    axes[1,0].set_title('Feasible Design Space')
    axes[1,0].legend(); axes[1,0].grid(True)

    # Quality histogram
    axes[1,1].hist(warpage, bins=25, alpha=0.6, color='steelblue', label='Warpage [mm]')
    axes[1,1].axvline(x=0.15, color='r', ls='--', label='Spec limit')
    axes[1,1].set_xlabel('Warpage [mm]')
    axes[1,1].set_ylabel('Count')
    axes[1,1].set_title('Warpage Distribution')
    axes[1,1].legend()

    # Digital twin architecture
    axes[1,2].axis('off')
    arch_text = (
        "Digital Twin Architecture\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━\n"
        "┌─────────┐   ┌──────────┐\n"
        "│ Sensors  │→→→│ Data Hub │\n"
        "│ (P,T,F)  │   │  (MQTT)  │\n"
        "└─────────┘   └────┬─────┘\n"
        "                   ↓\n"
        "┌─────────┐   ┌──────────┐\n"
        "│Moldflow/ │←→→│  EnKF    │\n"
        "│OpenFOAM  │   │Calibrator│\n"
        "└─────────┘   └────┬─────┘\n"
        "                   ↓\n"
        "┌─────────┐   ┌──────────┐\n"
        "│ NN/GPR  │→→→│ Quality  │\n"
        "│Surrogate│   │Prediction│\n"
        "└─────────┘   └──────────┘"
    )
    axes[1,2].text(0.5, 0.5, arch_text, transform=axes[1,2].transAxes,
                   fontsize=9, family='monospace', va='center', ha='center',
                   bbox=dict(boxstyle='round', facecolor='lightyellow'))

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/automotive_case_study.png', dpi=150)
    plt.close()

    metrics = {
        'r2_warpage': round(r2_warp, 4),
        'r2_shrinkage': round(r2_shrink, 4),
        'r2_surface': round(r2_surf, 4),
        'mae_warpage_mm': round(mae_warp, 4),
        'mae_shrinkage_pct': round(mae_shrink, 4),
        'feasible_count': int(good_mask.sum()),
        'total_samples': N,
        'feasible_ratio': round(good_mask.sum()/N, 4),
    }
    return metrics

# ============================================================
# 7. System Architecture Diagram
# ============================================================
def architecture_diagram():
    """Generate digital twin system architecture diagram."""
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    ax.set_xlim(0, 14); ax.set_ylim(0, 9)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Moldflow/OpenFOAM Digital Twin Architecture', fontsize=14, fontweight='bold')

    def draw_box(ax, xy, w, h, text, color='lightblue', fontsize=8):
        rect = plt.Rectangle(xy, w, h, fill=True, facecolor=color,
                              edgecolor='black', lw=1.5, zorder=2)
        ax.add_patch(rect)
        ax.text(xy[0]+w/2, xy[1]+h/2, text, ha='center', va='center',
                fontsize=fontsize, fontweight='bold', zorder=3)

    def draw_arrow(ax, start, end, color='black'):
        ax.annotate('', xy=end, xytext=start,
                     arrowprops=dict(arrowstyle='->', color=color, lw=1.5))

    # Physical Layer
    draw_box(ax, (0.5, 7), 3, 1.2, 'Physical Process\n(Injection Molding\nMachine)', '#FFB3B3')
    draw_box(ax, (4.5, 7), 2.5, 1.2, 'Sensors\n(P, T, Flow)', '#FFCC99')
    draw_box(ax, (8, 7), 2.5, 1.2, 'SCADA/PLC\nData Acquisition', '#FFCC99')
    draw_box(ax, (11.5, 7), 2, 1.2, 'Edge\nComputing', '#DDA0DD')

    # Data Layer
    draw_box(ax, (0.5, 4.8), 2.5, 1.5, 'Moldflow\nSimulation\n(Fill/Pack/Cool)', '#87CEEB')
    draw_box(ax, (3.5, 4.8), 2.5, 1.5, 'OpenFOAM\n3D CFD\n(Complex Flow)', '#87CEEB')
    draw_box(ax, (6.5, 4.8), 2.5, 1.5, 'Data Assimilation\nEnKF/UKF\nModel Calibration', '#98FB98')
    draw_box(ax, (9.5, 4.8), 2, 1.5, 'Time-Series\nDatabase\n(InfluxDB)', '#F0E68C')
    draw_box(ax, (12, 4.8), 1.5, 1.5, 'ML\nPipeline', '#DDA0DD')

    # AI/Quality Layer
    draw_box(ax, (0.5, 2.5), 2.5, 1.5, 'Surrogate Model\nGPR / ANN\n(Quality Pred.)', '#98FB98')
    draw_box(ax, (3.5, 2.5), 2.5, 1.5, 'Residual Stress\n& Warpage\nAnalysis', '#87CEEB')
    draw_box(ax, (6.5, 2.5), 2.5, 1.5, 'Multi-Objective\nOptimization\n(NSGA-II)', '#F0E68C')
    draw_box(ax, (9.5, 2.5), 4, 1.5, 'Quality Dashboard\n& Decision Support\n(Real-time Monitoring)', '#FFB3B3')

    # Output Layer
    draw_box(ax, (2, 0.5), 4, 1.2, 'Process Control\n& Parameter Adjustment', '#DDA0DD')
    draw_box(ax, (7, 0.5), 3, 1.2, 'Quality Report\n& SPC Charts', '#F0E68C')
    draw_box(ax, (11, 0.5), 2.5, 1.2, 'Predictive\nMaintenance', '#98FB98')

    # Arrows (main data flow)
    draw_arrow(ax, (3.5, 7.6), (4.5, 7.6))
    draw_arrow(ax, (7, 7.6), (8, 7.6))
    draw_arrow(ax, (10.5, 7.6), (11.5, 7.6))
    draw_arrow(ax, (2, 7), (1.75, 6.3))
    draw_arrow(ax, (5.75, 7), (4.75, 6.3))
    draw_arrow(ax, (9.25, 7), (7.75, 6.3))
    draw_arrow(ax, (12.5, 7), (12.75, 6.3))
    draw_arrow(ax, (1.75, 4.8), (1.75, 4))
    draw_arrow(ax, (4.75, 4.8), (4.75, 4))
    draw_arrow(ax, (7.75, 4.8), (7.75, 4))
    draw_arrow(ax, (10.5, 4.8), (10.5, 4))
    draw_arrow(ax, (4, 2.5), (4, 1.7))
    draw_arrow(ax, (8.5, 2.5), (8.5, 1.7))
    draw_arrow(ax, (11.5, 2.5), (12.25, 1.7))

    # Layer labels
    ax.text(14, 7.6, 'Physical\nLayer', fontsize=9, va='center', color='red', fontweight='bold')
    ax.text(14, 5.5, 'Simulation\n& Data\nLayer', fontsize=9, va='center', color='blue', fontweight='bold')
    ax.text(14, 3.25, 'AI &\nQuality\nLayer', fontsize=9, va='center', color='green', fontweight='bold')
    ax.text(14, 1.1, 'Output\nLayer', fontsize=9, va='center', color='purple', fontweight='bold')

    plt.tight_layout()
    plt.savefig(f'{FIGDIR}/architecture.png', dpi=150, bbox_inches='tight')
    plt.close()

# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    print("=" * 60)
    print("Digital Twin for Injection Molding - Simulation Suite")
    print("=" * 60)

    print("\n[1/7] Hele-Shaw Flow Simulation...")
    P, V = hele_shaw_flow()
    print(f"  Max pressure: {P.max()/1e6:.1f} MPa, Max velocity: {V.max():.4f} m/s")

    print("\n[2/7] Cooling & Crystallization...")
    times, T_hist, alpha_hist = cooling_crystallization()
    print(f"  Final center T: {T_hist[-1][15]:.1f}°C, Max crystallinity: {max(a.max() for a in alpha_hist)*100:.1f}%")

    print("\n[3/7] Residual Stress & Warpage...")
    sigma, warp = residual_stress_warpage()
    print(f"  Max stress: {abs(sigma).max():.1f} MPa, Max warpage: {warp:.4f} mm")

    print("\n[4/7] Process Parameter Analysis...")
    imp, warp_data, wt_data, sink_data = process_parameter_analysis()
    print(f"  Warpage range: [{warp_data.min():.3f}, {warp_data.max():.3f}] mm")

    print("\n[5/7] Data Assimilation (EnKF)...")
    rmse_p, rmse_e, k_est, k_true = data_assimilation_enkf()
    print(f"  Prior RMSE: {rmse_p:.2f}°C, EnKF RMSE: {rmse_e:.2f}°C")
    print(f"  Estimated k: {k_est:.4f} (true: {k_true})")

    print("\n[6/7] Automotive Case Study...")
    metrics = automotive_case_study()
    print(f"  R² warpage: {metrics['r2_warpage']}")
    print(f"  R² shrinkage: {metrics['r2_shrinkage']}")
    print(f"  R² surface: {metrics['r2_surface']}")
    print(f"  Feasible designs: {metrics['feasible_count']}/{metrics['total_samples']}")

    print("\n[7/7] Architecture Diagram...")
    architecture_diagram()
    print("  Architecture diagram saved.")

    # Save metrics
    all_metrics = {
        'hele_shaw': {'max_pressure_MPa': round(P.max()/1e6, 1), 'max_velocity_m_s': round(V.max(), 4)},
        'cooling': {'final_center_T': round(T_hist[-1][15], 1)},
        'stress_warpage': {'max_stress_MPa': round(abs(sigma).max(), 1), 'max_warpage_mm': round(warp, 4)},
        'data_assimilation': {'prior_rmse': round(rmse_p, 2), 'enkf_rmse': round(rmse_e, 2),
                              'k_estimated': round(k_est, 4), 'k_true': k_true,
                              'rmse_reduction_pct': round((1 - rmse_e/rmse_p)*100, 1)},
        'automotive': metrics
    }
    with open('metrics.json', 'w') as f:
        json.dump(all_metrics, f, indent=2, default=lambda x: int(x) if isinstance(x, np.integer) else float(x))

    print("\n" + "=" * 60)
    print("All simulations complete. Figures saved to figures/")
    print("Metrics saved to metrics.json")
    print("=" * 60)
