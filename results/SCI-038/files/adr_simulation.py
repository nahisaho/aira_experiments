#!/usr/bin/env python3
"""
Active Debris Removal (ADR) Mission Optimal Trajectory Design System

Modules:
1. Debris catalog scoring (collision risk × removal effectiveness)
2. Multi-target low-thrust trajectory optimization
3. Rendezvous & proximity operations (Hill/CW equations)
4. Tumbling debris attitude estimation
5. Capture mechanism dynamics (robotic arm / net / harpoon)
6. Mission sequence optimization (cost minimization)
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, differential_evolution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import permutations
import os
import json

np.random.seed(42)
FIGURES_DIR = "figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

# Constants
MU_EARTH = 3.986004418e14   # m^3/s^2
R_EARTH = 6371e3             # m
J2 = 1.08263e-3

# ============================================================
# Module 1: Debris Catalog & Target Scoring
# ============================================================

def generate_debris_catalog(n=30):
    """Generate synthetic debris catalog with orbital elements and properties."""
    catalog = []
    for i in range(n):
        alt = np.random.uniform(600e3, 1200e3)
        a = R_EARTH + alt
        e = np.random.uniform(0.0001, 0.02)
        inc = np.random.uniform(60, 100)  # deg - SSO-like
        raan = np.random.uniform(0, 360)
        aop = np.random.uniform(0, 360)
        mass = np.random.uniform(100, 3000)  # kg
        area = np.random.uniform(1, 30)  # m^2 cross-section
        collision_prob = np.random.uniform(1e-6, 1e-3)
        catalog.append({
            'id': f'DEB-{i+1:04d}',
            'a': a, 'e': e, 'inc': inc, 'raan': raan, 'aop': aop,
            'mass': mass, 'area': area,
            'collision_prob': collision_prob,
            'altitude_km': alt / 1e3,
        })
    return catalog


def compute_debris_scores(catalog):
    """Score debris by collision risk × removal effectiveness.
    
    Environmental impact score = collision_prob × mass × (area / mass) × lifetime_factor
    """
    scores = []
    for d in catalog:
        n = np.sqrt(MU_EARTH / d['a']**3)
        # Approximate atmospheric lifetime factor (higher alt → longer lived → more dangerous)
        lifetime_factor = np.exp((d['altitude_km'] - 600) / 200)
        # Collision risk component
        collision_risk = d['collision_prob'] * lifetime_factor
        # Removal effectiveness (larger, more massive → more fragments if hit)
        removal_effect = d['mass'] * d['area'] / 1000.0
        score = collision_risk * removal_effect
        d['score'] = score
        d['collision_risk'] = collision_risk
        d['removal_effect'] = removal_effect
        d['mean_motion'] = n
        scores.append(score)
    return catalog


def select_targets(catalog, n_targets=5):
    """Select top-N targets based on combined score."""
    sorted_cat = sorted(catalog, key=lambda x: x['score'], reverse=True)
    return sorted_cat[:n_targets]


def plot_scoring(catalog, selected, filename="debris_scoring.png"):
    """Plot debris scoring results."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # All debris
    alts = [d['altitude_km'] for d in catalog]
    scores = [d['score'] for d in catalog]
    masses = [d['mass'] for d in catalog]
    
    sel_alts = [d['altitude_km'] for d in selected]
    sel_scores = [d['score'] for d in selected]
    sel_masses = [d['mass'] for d in selected]
    
    # Score vs altitude
    axes[0].scatter(alts, scores, c='gray', alpha=0.5, label='All debris')
    axes[0].scatter(sel_alts, sel_scores, c='red', s=100, marker='*', label='Selected targets')
    axes[0].set_xlabel('Altitude (km)')
    axes[0].set_ylabel('Combined Score')
    axes[0].set_title('Debris Scoring: Altitude vs Score')
    axes[0].legend()
    axes[0].set_yscale('log')
    axes[0].grid(True, alpha=0.3)
    
    # Score vs mass
    axes[1].scatter(masses, scores, c='gray', alpha=0.5, label='All debris')
    axes[1].scatter(sel_masses, sel_scores, c='red', s=100, marker='*', label='Selected targets')
    axes[1].set_xlabel('Mass (kg)')
    axes[1].set_ylabel('Combined Score')
    axes[1].set_title('Debris Scoring: Mass vs Score')
    axes[1].legend()
    axes[1].set_yscale('log')
    axes[1].grid(True, alpha=0.3)
    
    # Bar chart of selected
    ids = [d['id'] for d in selected]
    axes[2].barh(ids, sel_scores, color='crimson')
    axes[2].set_xlabel('Combined Score')
    axes[2].set_title('Top-5 Selected Targets')
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved {filename}")


# ============================================================
# Module 2: Low-Thrust Trajectory Optimization
# ============================================================

def hohmann_delta_v(a1, a2):
    """Compute Hohmann transfer ΔV between two circular orbits."""
    v1 = np.sqrt(MU_EARTH / a1)
    v2 = np.sqrt(MU_EARTH / a2)
    at = (a1 + a2) / 2
    vt1 = np.sqrt(MU_EARTH * (2/a1 - 1/at))
    vt2 = np.sqrt(MU_EARTH * (2/a2 - 1/at))
    dv1 = abs(vt1 - v1)
    dv2 = abs(v2 - vt2)
    return dv1 + dv2


def plane_change_delta_v(a, inc1_deg, inc2_deg):
    """ΔV for simple plane change."""
    v = np.sqrt(MU_EARTH / a)
    di = np.radians(abs(inc2_deg - inc1_deg))
    return 2 * v * np.sin(di / 2)


def transfer_delta_v(d1, d2):
    """Total ΔV to transfer between two debris orbits (in-plane + out-of-plane)."""
    dv_ip = hohmann_delta_v(d1['a'], d2['a'])
    dv_oop = plane_change_delta_v((d1['a'] + d2['a'])/2, d1['inc'], d2['inc'])
    # RAAN change cost approximation
    raan_diff = abs(d1['raan'] - d2['raan']) % 360
    if raan_diff > 180:
        raan_diff = 360 - raan_diff
    dv_raan = 0.01 * raan_diff  # simplified
    return np.sqrt(dv_ip**2 + dv_oop**2) + dv_raan


def transfer_time(d1, d2):
    """Approximate transfer time (Hohmann + drift)."""
    at = (d1['a'] + d2['a']) / 2
    t_hohmann = np.pi * np.sqrt(at**3 / MU_EARTH)
    raan_diff = abs(d1['raan'] - d2['raan']) % 360
    if raan_diff > 180:
        raan_diff = 360 - raan_diff
    # J2 drift rate
    n1 = np.sqrt(MU_EARTH / d1['a']**3)
    raan_dot = -1.5 * n1 * J2 * (R_EARTH / d1['a'])**2 * np.cos(np.radians(d1['inc']))
    if abs(raan_dot) > 1e-12:
        t_drift = np.radians(raan_diff) / abs(raan_dot)
    else:
        t_drift = 0
    return t_hohmann + min(t_drift, 365*86400)


def compute_cost_matrix(targets):
    """Compute ΔV cost matrix between all target pairs."""
    n = len(targets)
    dv_matrix = np.zeros((n, n))
    time_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                dv_matrix[i, j] = transfer_delta_v(targets[i], targets[j])
                time_matrix[i, j] = transfer_time(targets[i], targets[j])
    return dv_matrix, time_matrix


def low_thrust_trajectory_sim(a_start, a_end, thrust_acc=1e-4, dt=60):
    """Simulate low-thrust spiral transfer."""
    a = a_start
    trajectory = [a]
    times = [0]
    t = 0
    direction = 1.0 if a_end > a_start else -1.0
    while (direction > 0 and a < a_end) or (direction < 0 and a > a_end):
        v = np.sqrt(MU_EARTH / a)
        da = 2 * a**2 * thrust_acc * direction / (MU_EARTH / a)**0.5 * dt
        a += da
        t += dt
        trajectory.append(a)
        times.append(t)
        if t > 1e7:
            break
    return np.array(times), np.array(trajectory)


def optimize_sequence_ga(targets, pop_size=100, generations=200):
    """Genetic algorithm for debris visit sequence optimization."""
    n = len(targets)
    dv_matrix, time_matrix = compute_cost_matrix(targets)
    
    # Initialize population
    population = [np.random.permutation(n) for _ in range(pop_size)]
    
    def fitness(seq):
        total_dv = sum(dv_matrix[seq[i], seq[i+1]] for i in range(len(seq)-1))
        total_time = sum(time_matrix[seq[i], seq[i+1]] for i in range(len(seq)-1))
        return total_dv + 0.001 * total_time / 86400  # weighted
    
    best_fitness_history = []
    
    for gen in range(generations):
        fitnesses = [fitness(ind) for ind in population]
        sorted_idx = np.argsort(fitnesses)
        population = [population[i] for i in sorted_idx]
        fitnesses = [fitnesses[i] for i in sorted_idx]
        best_fitness_history.append(fitnesses[0])
        
        # Elitism + crossover + mutation
        new_pop = population[:pop_size // 5]  # elites
        while len(new_pop) < pop_size:
            p1 = population[np.random.randint(pop_size // 3)]
            p2 = population[np.random.randint(pop_size // 3)]
            # Order crossover
            c1, c2 = np.random.randint(0, n, 2)
            if c1 > c2: c1, c2 = c2, c1
            child = np.full(n, -1)
            child[c1:c2] = p1[c1:c2]
            fill = [g for g in p2 if g not in child[c1:c2]]
            idx = 0
            for i in range(n):
                if child[i] == -1:
                    child[i] = fill[idx]
                    idx += 1
            # Mutation
            if np.random.random() < 0.3:
                i1, i2 = np.random.randint(0, n, 2)
                child[i1], child[i2] = child[i2], child[i1]
            new_pop.append(child)
        population = new_pop[:pop_size]
    
    best_seq = population[0]
    return best_seq, best_fitness_history, dv_matrix, time_matrix


def plot_trajectory_optimization(targets, best_seq, fitness_history, dv_matrix,
                                  filename="trajectory_optimization.png"):
    """Plot trajectory optimization results."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Fitness convergence
    axes[0, 0].plot(fitness_history, 'b-', linewidth=1.5)
    axes[0, 0].set_xlabel('Generation')
    axes[0, 0].set_ylabel('Total Cost (ΔV + time penalty)')
    axes[0, 0].set_title('GA Convergence')
    axes[0, 0].grid(True, alpha=0.3)
    
    # ΔV matrix heatmap
    im = axes[0, 1].imshow(dv_matrix, cmap='YlOrRd')
    axes[0, 1].set_xlabel('Target Index')
    axes[0, 1].set_ylabel('Target Index')
    axes[0, 1].set_title('ΔV Transfer Cost Matrix (m/s)')
    plt.colorbar(im, ax=axes[0, 1])
    
    # Optimal sequence in orbital element space
    alts = [targets[i]['altitude_km'] for i in best_seq]
    incs = [targets[i]['inc'] for i in best_seq]
    ids = [targets[i]['id'] for i in best_seq]
    axes[1, 0].plot(alts, incs, 'ro-', markersize=10, linewidth=2)
    for k, (alt, inc, tid) in enumerate(zip(alts, incs, ids)):
        axes[1, 0].annotate(f'{k+1}:{tid}', (alt, inc), fontsize=8,
                           textcoords="offset points", xytext=(5, 5))
    axes[1, 0].set_xlabel('Altitude (km)')
    axes[1, 0].set_ylabel('Inclination (deg)')
    axes[1, 0].set_title('Optimal Visit Sequence')
    axes[1, 0].grid(True, alpha=0.3)
    
    # ΔV per leg
    dvs = [dv_matrix[best_seq[i], best_seq[i+1]] for i in range(len(best_seq)-1)]
    leg_labels = [f'{ids[i]}→{ids[i+1]}' for i in range(len(ids)-1)]
    axes[1, 1].bar(range(len(dvs)), dvs, color='steelblue')
    axes[1, 1].set_xticks(range(len(dvs)))
    axes[1, 1].set_xticklabels(leg_labels, rotation=45, ha='right', fontsize=8)
    axes[1, 1].set_ylabel('ΔV (m/s)')
    axes[1, 1].set_title('ΔV per Transfer Leg')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved {filename}")


def plot_low_thrust_transfer(targets, best_seq, filename="low_thrust_transfer.png"):
    """Simulate and plot low-thrust spiral transfers."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(best_seq)-1))
    all_times = []
    all_alts = []
    cumulative_t = 0
    
    for k in range(len(best_seq) - 1):
        d1 = targets[best_seq[k]]
        d2 = targets[best_seq[k+1]]
        t, a = low_thrust_trajectory_sim(d1['a'], d2['a'], thrust_acc=5e-5)
        axes[0].plot(t / 86400, (a - R_EARTH) / 1e3, color=colors[k],
                    label=f"Leg {k+1}: {d1['id']}→{d2['id']}")
        all_times.append(t + cumulative_t)
        all_alts.append((a - R_EARTH) / 1e3)
        cumulative_t += t[-1]
    
    axes[0].set_xlabel('Time (days)')
    axes[0].set_ylabel('Altitude (km)')
    axes[0].set_title('Low-Thrust Spiral Transfers (per leg)')
    axes[0].legend(fontsize=7)
    axes[0].grid(True, alpha=0.3)
    
    # Combined timeline
    for k in range(len(all_times)):
        axes[1].plot(all_times[k] / 86400, all_alts[k], color=colors[k], linewidth=2)
    axes[1].set_xlabel('Mission Time (days)')
    axes[1].set_ylabel('Altitude (km)')
    axes[1].set_title('Complete Mission Altitude Profile')
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved {filename}")


# ============================================================
# Module 3: Hill Equation Rendezvous & Proximity Operations
# ============================================================

def hill_equations(t, state, n):
    """Clohessy-Wiltshire (Hill) equations of relative motion."""
    x, y, z, vx, vy, vz = state
    ax = 3*n**2*x + 2*n*vy
    ay = -2*n*vx
    az = -n**2*z
    return [vx, vy, vz, ax, ay, az]


def hill_equations_controlled(t, state, n, thrust_profile):
    """Hill equations with thrust input for RPO maneuvers."""
    x, y, z, vx, vy, vz = state
    fx, fy, fz = thrust_profile(t, state)
    ax = 3*n**2*x + 2*n*vy + fx
    ay = -2*n*vx + fy
    az = -n**2*z + fz
    return [vx, vy, vz, ax, ay, az]


def plan_rpo_approach(n, x0, v0, xf, vf, t_maneuver):
    """Plan two-impulse rendezvous using CW state transition matrix."""
    nt = n * t_maneuver
    sn = np.sin(nt)
    cn = np.cos(nt)
    
    # CW state transition matrix
    phi_rr = np.array([
        [4 - 3*cn, 0, 0],
        [6*(sn - nt), 1, 0],
        [0, 0, cn]
    ])
    phi_rv = np.array([
        [sn/n, 2*(1-cn)/n, 0],
        [-2*(1-cn)/n, (4*sn - 3*nt)/n, 0],
        [0, 0, sn/n]
    ])
    
    x0_vec = np.array(x0)
    xf_vec = np.array(xf)
    
    # Required initial velocity for rendezvous
    v0_req = np.linalg.solve(phi_rv, xf_vec - phi_rr @ x0_vec)
    
    dv1 = v0_req - np.array(v0)
    
    # Compute final velocity
    phi_vr = np.array([
        [3*n*sn, 0, 0],
        [-6*n*(1-cn), 0, 0],
        [0, 0, -n*sn]
    ])
    phi_vv = np.array([
        [cn, 2*sn, 0],
        [-2*sn, 4*cn - 3, 0],
        [0, 0, cn]
    ])
    vf_actual = phi_vr @ x0_vec + phi_vv @ v0_req
    dv2 = np.array(vf) - vf_actual
    
    return dv1, dv2, v0_req


def simulate_rpo(n, x0, v0, duration, thrust_profile=None, dt_out=1.0):
    """Simulate relative motion for RPO."""
    state0 = list(x0) + list(v0)
    if thrust_profile is None:
        sol = solve_ivp(hill_equations, [0, duration], state0, args=(n,),
                       max_step=dt_out, dense_output=True)
    else:
        sol = solve_ivp(lambda t, s: hill_equations_controlled(t, s, n, thrust_profile),
                       [0, duration], state0, max_step=dt_out, dense_output=True)
    t_eval = np.linspace(0, duration, int(duration / dt_out))
    states = sol.sol(t_eval).T
    return t_eval, states


def plot_rpo_results(t_eval, states, dv1, dv2, n, filename="rpo_simulation.png"):
    """Plot RPO simulation results."""
    fig = plt.figure(figsize=(16, 12))
    
    # 3D trajectory
    ax1 = fig.add_subplot(2, 2, 1, projection='3d')
    ax1.plot(states[:, 0], states[:, 1], states[:, 2], 'b-', linewidth=1.5)
    ax1.scatter(*states[0, :3], c='green', s=100, marker='o', label='Start')
    ax1.scatter(*states[-1, :3], c='red', s=100, marker='*', label='End')
    ax1.scatter(0, 0, 0, c='orange', s=200, marker='s', label='Target')
    ax1.set_xlabel('X - Radial (m)')
    ax1.set_ylabel('Y - Along-track (m)')
    ax1.set_zlabel('Z - Cross-track (m)')
    ax1.set_title('RPO 3D Trajectory (LVLH Frame)')
    ax1.legend(fontsize=8)
    
    # XY plane (V-bar approach)
    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(states[:, 1], states[:, 0], 'b-', linewidth=1.5)
    ax2.scatter(states[0, 1], states[0, 0], c='green', s=100, marker='o', label='Start')
    ax2.scatter(states[-1, 1], states[-1, 0], c='red', s=100, marker='*', label='End')
    ax2.scatter(0, 0, c='orange', s=200, marker='s', label='Target')
    ax2.set_xlabel('Along-track Y (m)')
    ax2.set_ylabel('Radial X (m)')
    ax2.set_title('RPO - XY Plane (V-bar Approach)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_aspect('equal')
    
    # Distance over time
    ax3 = fig.add_subplot(2, 2, 3)
    dist = np.sqrt(states[:, 0]**2 + states[:, 1]**2 + states[:, 2]**2)
    ax3.plot(t_eval / 60, dist, 'b-', linewidth=1.5)
    ax3.set_xlabel('Time (min)')
    ax3.set_ylabel('Range (m)')
    ax3.set_title('Range to Target vs Time')
    ax3.grid(True, alpha=0.3)
    
    # Velocity components
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.plot(t_eval / 60, states[:, 3], label='Vx (radial)')
    ax4.plot(t_eval / 60, states[:, 4], label='Vy (along-track)')
    ax4.plot(t_eval / 60, states[:, 5], label='Vz (cross-track)')
    ax4.set_xlabel('Time (min)')
    ax4.set_ylabel('Velocity (m/s)')
    ax4.set_title('Relative Velocity Components')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.suptitle(f'ΔV₁ = {np.linalg.norm(dv1):.3f} m/s, ΔV₂ = {np.linalg.norm(dv2):.3f} m/s', 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved {filename}")


# ============================================================
# Module 4: Tumbling Debris Attitude Estimation
# ============================================================

def euler_equations(t, state, I_body, torque_func=None):
    """Euler's equations for rigid body rotational dynamics.
    state = [q0, q1, q2, q3, wx, wy, wz]
    """
    q = state[:4]
    w = state[4:]
    Ix, Iy, Iz = I_body
    
    # External torque
    if torque_func:
        tx, ty, tz = torque_func(t)
    else:
        tx, ty, tz = 0, 0, 0
    
    # Quaternion kinematics
    q0, q1, q2, q3 = q
    wx, wy, wz = w
    dq0 = 0.5 * (-q1*wx - q2*wy - q3*wz)
    dq1 = 0.5 * (q0*wx + q2*wz - q3*wy)
    dq2 = 0.5 * (q0*wy - q1*wz + q3*wx)
    dq3 = 0.5 * (q0*wz + q1*wy - q2*wx)
    
    # Euler equations
    dwx = ((Iy - Iz) * wy * wz + tx) / Ix
    dwy = ((Iz - Ix) * wz * wx + ty) / Iy
    dwz = ((Ix - Iy) * wx * wy + tz) / Iz
    
    return [dq0, dq1, dq2, dq3, dwx, dwy, dwz]


def simulate_tumbling(I_body, w0, duration, dt=0.1, torque_func=None):
    """Simulate tumbling debris rotation."""
    q0 = [1, 0, 0, 0]  # Initial quaternion (identity)
    state0 = q0 + list(w0)
    
    sol = solve_ivp(lambda t, s: euler_equations(t, s, I_body, torque_func),
                   [0, duration], state0, max_step=dt, dense_output=True)
    
    t_eval = np.linspace(0, duration, int(duration / dt))
    states = sol.sol(t_eval).T
    
    # Normalize quaternions
    for i in range(len(states)):
        qnorm = np.linalg.norm(states[i, :4])
        states[i, :4] /= qnorm
    
    return t_eval, states


def estimate_tumbling_rate(t, states, noise_std=0.02):
    """Estimate tumbling rate from noisy measurements (Extended Kalman Filter-like)."""
    w_true = states[:, 4:]
    w_measured = w_true + np.random.randn(*w_true.shape) * noise_std
    
    # Simple moving average filter as baseline estimator
    window = 20
    w_estimated = np.zeros_like(w_measured)
    for i in range(len(w_measured)):
        start = max(0, i - window)
        w_estimated[i] = np.mean(w_measured[start:i+1], axis=0)
    
    # Compute estimation error
    errors = np.linalg.norm(w_true - w_estimated, axis=1)
    
    return w_measured, w_estimated, errors


def plot_tumbling(t, states, w_measured, w_estimated, errors, 
                  filename="tumbling_analysis.png"):
    """Plot tumbling debris analysis."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Angular velocities - true
    axes[0, 0].plot(t, np.degrees(states[:, 4]), label='ωx')
    axes[0, 0].plot(t, np.degrees(states[:, 5]), label='ωy')
    axes[0, 0].plot(t, np.degrees(states[:, 6]), label='ωz')
    axes[0, 0].set_xlabel('Time (s)')
    axes[0, 0].set_ylabel('Angular Rate (deg/s)')
    axes[0, 0].set_title('True Angular Velocity')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Measured vs estimated (x-component)
    axes[0, 1].plot(t, np.degrees(w_measured[:, 0]), 'gray', alpha=0.3, label='Measured')
    axes[0, 1].plot(t, np.degrees(w_estimated[:, 0]), 'r-', linewidth=2, label='Estimated')
    axes[0, 1].plot(t, np.degrees(states[:, 4]), 'b--', linewidth=1, label='True')
    axes[0, 1].set_xlabel('Time (s)')
    axes[0, 1].set_ylabel('ωx (deg/s)')
    axes[0, 1].set_title('Estimation: ωx Component')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Estimation error
    axes[1, 0].plot(t, np.degrees(errors), 'r-')
    axes[1, 0].set_xlabel('Time (s)')
    axes[1, 0].set_ylabel('Error (deg/s)')
    axes[1, 0].set_title('Estimation Error (L2 norm)')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Quaternion evolution
    axes[1, 1].plot(t, states[:, 0], label='q0')
    axes[1, 1].plot(t, states[:, 1], label='q1')
    axes[1, 1].plot(t, states[:, 2], label='q2')
    axes[1, 1].plot(t, states[:, 3], label='q3')
    axes[1, 1].set_xlabel('Time (s)')
    axes[1, 1].set_ylabel('Quaternion Component')
    axes[1, 1].set_title('Attitude Quaternion Evolution')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved {filename}")


# ============================================================
# Module 5: Capture Mechanism Dynamics
# ============================================================

def simulate_robotic_arm(joint_angles_init, target_pos, 
                          link_lengths=[2.0, 1.5, 1.0],
                          duration=30, dt=0.1):
    """Simulate 3-link planar robotic arm capture maneuver."""
    n_steps = int(duration / dt)
    n_links = len(link_lengths)
    angles = np.array(joint_angles_init, dtype=float)
    angles_history = [angles.copy()]
    ee_history = []
    torque_history = []
    
    def forward_kinematics(angles):
        x, y = 0, 0
        positions = [(x, y)]
        cumulative_angle = 0
        for i, (a, l) in enumerate(zip(angles, link_lengths)):
            cumulative_angle += a
            x += l * np.cos(cumulative_angle)
            y += l * np.sin(cumulative_angle)
            positions.append((x, y))
        return np.array(positions), np.array([x, y])
    
    for step in range(n_steps):
        positions, ee_pos = forward_kinematics(angles)
        ee_history.append(ee_pos.copy())
        
        error = target_pos - ee_pos
        dist = np.linalg.norm(error)
        
        if dist < 0.05:
            for _ in range(n_steps - step - 1):
                ee_history.append(ee_pos.copy())
                angles_history.append(angles.copy())
                torque_history.append(np.zeros(n_links))
            break
        
        # Jacobian (numerical)
        J = np.zeros((2, n_links))
        eps = 1e-6
        for i in range(n_links):
            angles_plus = angles.copy()
            angles_plus[i] += eps
            _, ee_plus = forward_kinematics(angles_plus)
            J[:, i] = (ee_plus - ee_pos) / eps
        
        # Damped least squares IK
        lam = 0.1
        JT = J.T
        dq = JT @ np.linalg.solve(J @ JT + lam**2 * np.eye(2), error)
        
        max_rate = 0.1
        dq = np.clip(dq, -max_rate, max_rate)
        
        torques = dq * 10  # simplified torque
        torque_history.append(torques)
        
        angles += dq
        angles_history.append(angles.copy())
    
    return (np.array(angles_history), np.array(ee_history), 
            np.array(torque_history), link_lengths)


def simulate_net_capture(net_size=5.0, target_pos=np.array([10.0, 0.0, 0.0]),
                         target_vel=np.array([0.5, 0.1, 0.0]),
                         launch_vel=2.0, duration=15, dt=0.01):
    """Simulate tethered net deployment and capture dynamics."""
    n_corners = 4
    # Net corner initial positions (folded)
    offsets = np.array([[0.1, 0.1, 0], [-0.1, 0.1, 0],
                        [-0.1, -0.1, 0], [0.1, -0.1, 0]])
    
    launch_dir = target_pos / np.linalg.norm(target_pos)
    spread_rate = 0.3
    
    corner_pos = [offsets[i].copy() for i in range(n_corners)]
    corner_vel = [launch_dir * launch_vel + offsets[i] * spread_rate * 10 
                  for i in range(n_corners)]
    
    history = {'corners': [], 'center': [], 'target': [], 'time': []}
    t = 0
    captured = False
    capture_time = None
    
    while t < duration:
        center = np.mean(corner_pos, axis=0)
        history['corners'].append([c.copy() for c in corner_pos])
        history['center'].append(center.copy())
        history['target'].append(target_pos.copy())
        history['time'].append(t)
        
        # Check capture
        target_in_net = np.linalg.norm(target_pos - center) < net_size / 2
        if target_in_net and not captured:
            captured = True
            capture_time = t
        
        if captured:
            # Post-capture: net wraps, deceleration
            for i in range(n_corners):
                corner_vel[i] *= 0.95
                corner_pos[i] += corner_vel[i] * dt
            target_vel *= 0.98
            target_pos += target_vel * dt
        else:
            # Free flight with drag
            for i in range(n_corners):
                corner_pos[i] += corner_vel[i] * dt
            target_pos += target_vel * dt
        
        t += dt
    
    return history, captured, capture_time


def simulate_harpoon(target_distance=15, target_vel=0.3, 
                      harpoon_vel=5.0, duration=10, dt=0.001):
    """Simulate harpoon capture dynamics."""
    # State: [harpoon_pos, harpoon_vel, target_pos, target_vel, tether_tension]
    h_pos = 0
    h_vel = harpoon_vel
    t_pos = target_distance
    t_vel = target_vel
    tether_length = 0
    tether_stiffness = 500  # N/m
    tether_damping = 50     # Ns/m
    harpoon_mass = 2.0      # kg
    target_mass = 500       # kg
    
    history = {'time': [], 'h_pos': [], 't_pos': [], 'h_vel': [], 't_vel': [],
               'tension': [], 'separation': []}
    
    t = 0
    penetrated = False
    penetration_time = None
    
    while t < duration:
        separation = t_pos - h_pos
        history['time'].append(t)
        history['h_pos'].append(h_pos)
        history['t_pos'].append(t_pos)
        history['h_vel'].append(h_vel)
        history['t_vel'].append(t_vel)
        history['separation'].append(separation)
        
        if not penetrated and separation <= 0:
            penetrated = True
            penetration_time = t
            # Momentum transfer at impact
            v_combined = (harpoon_mass * h_vel + target_mass * t_vel) / (harpoon_mass + target_mass)
            h_vel = v_combined
            t_vel = v_combined
        
        if penetrated:
            # Tether dynamics
            tether_stretch = max(0, (h_pos - 0) - tether_length)  # simplified
            tension = 0
            if tether_stretch > 0:
                tension = tether_stiffness * tether_stretch + tether_damping * (h_vel - 0)
            history['tension'].append(tension)
            h_pos += h_vel * dt
            t_pos += t_vel * dt
        else:
            history['tension'].append(0)
            h_pos += h_vel * dt
            t_pos += t_vel * dt
        
        t += dt
    
    return history, penetrated, penetration_time


def plot_capture_mechanisms(arm_data, net_data, harpoon_data, 
                            filename="capture_mechanisms.png"):
    """Plot capture mechanism simulation results."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    # ---- Robotic Arm ----
    angles_hist, ee_hist, torque_hist, link_lengths = arm_data
    ax = axes[0, 0]
    ax.plot(ee_hist[:, 0], ee_hist[:, 1], 'b-', linewidth=1.5, label='End effector path')
    ax.scatter(ee_hist[0, 0], ee_hist[0, 1], c='green', s=100, marker='o', label='Start')
    ax.scatter(ee_hist[-1, 0], ee_hist[-1, 1], c='red', s=100, marker='*', label='Target')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_title('Robotic Arm: End Effector Path')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    ax = axes[1, 0]
    if len(torque_hist) > 0:
        for i in range(torque_hist.shape[1]):
            ax.plot(np.arange(len(torque_hist)) * 0.1, torque_hist[:, i], 
                   label=f'Joint {i+1}')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Torque (Nm)')
    ax.set_title('Robotic Arm: Joint Torques')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # ---- Net Capture ----
    net_history, net_captured, net_cap_time = net_data
    ax = axes[0, 1]
    centers = np.array(net_history['center'])
    targets = np.array(net_history['target'])
    ax.plot(centers[:, 0], centers[:, 1], 'b-', label='Net center')
    ax.plot(targets[:, 0], targets[:, 1], 'r--', label='Target')
    if net_captured:
        cap_idx = np.argmin(np.abs(np.array(net_history['time']) - net_cap_time))
        ax.axvline(x=centers[cap_idx, 0], color='g', linestyle=':', alpha=0.5)
        ax.scatter(centers[cap_idx, 0], centers[cap_idx, 1], c='green', s=200, 
                  marker='D', zorder=5, label=f'Capture t={net_cap_time:.1f}s')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_title(f'Net Capture {"✓" if net_captured else "✗"}')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    
    ax = axes[1, 1]
    times = net_history['time']
    separations = [np.linalg.norm(np.array(c) - np.array(t)) 
                   for c, t in zip(net_history['center'], net_history['target'])]
    ax.plot(times, separations, 'b-')
    if net_captured:
        ax.axvline(x=net_cap_time, color='r', linestyle='--', label=f'Capture @ {net_cap_time:.1f}s')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Separation (m)')
    ax.set_title('Net-Target Separation')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # ---- Harpoon ----
    h_hist, h_penetrated, h_pen_time = harpoon_data
    ax = axes[0, 2]
    ax.plot(h_hist['time'], h_hist['h_pos'], 'b-', label='Harpoon')
    ax.plot(h_hist['time'], h_hist['t_pos'], 'r-', label='Target')
    if h_penetrated:
        ax.axvline(x=h_pen_time, color='g', linestyle='--', 
                  label=f'Impact @ {h_pen_time:.2f}s')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Position (m)')
    ax.set_title(f'Harpoon Capture {"✓" if h_penetrated else "✗"}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    ax = axes[1, 2]
    ax.plot(h_hist['time'], h_hist['tension'], 'r-')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Tension (N)')
    ax.set_title('Tether Tension')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved {filename}")


# ============================================================
# Module 6: Mission Sequence Optimization
# ============================================================

def mission_cost_function(sequence, targets, dv_matrix, time_matrix,
                          capture_time=3600*4, deorbit_dv=100):
    """Compute total mission cost for a given sequence."""
    n = len(sequence)
    total_dv = 0
    total_time = 0
    total_fuel = 0
    
    Isp = 3000  # s (electric propulsion)
    g0 = 9.81
    m_spacecraft = 2000  # kg initial wet mass
    m_current = m_spacecraft
    
    for i in range(n - 1):
        leg_dv = dv_matrix[sequence[i], sequence[i+1]]
        leg_time = time_matrix[sequence[i], sequence[i+1]]
        total_dv += leg_dv
        total_time += leg_time + capture_time
        # Tsiolkovsky fuel mass
        m_fuel = m_current * (1 - np.exp(-leg_dv / (Isp * g0)))
        total_fuel += m_fuel
        m_current -= m_fuel
    
    # Add deorbit ΔV for each target
    total_dv += deorbit_dv * n
    
    # Cost model: propellant_cost + time_cost + operations_cost
    fuel_cost = total_fuel * 50000  # $/kg to orbit
    time_cost = total_time / 86400 * 10000  # $/day operations
    launch_cost = 50e6  # base launch cost
    
    total_cost = launch_cost + fuel_cost + time_cost
    
    return {
        'total_dv': total_dv,
        'total_time_days': total_time / 86400,
        'total_fuel_kg': total_fuel,
        'fuel_cost': fuel_cost,
        'time_cost': time_cost,
        'launch_cost': launch_cost,
        'total_cost': total_cost,
        'final_mass': m_current,
    }


def optimize_mission_cost(targets, dv_matrix, time_matrix):
    """Exhaustive + heuristic mission sequence cost optimization."""
    n = len(targets)
    
    if n <= 7:
        # Exact for small n
        best_cost = float('inf')
        best_seq = None
        all_results = []
        for perm in permutations(range(n)):
            result = mission_cost_function(list(perm), targets, dv_matrix, time_matrix)
            all_results.append((list(perm), result))
            if result['total_cost'] < best_cost:
                best_cost = result['total_cost']
                best_seq = list(perm)
        return best_seq, all_results
    else:
        # GA for larger sets
        best_seq, _, _, _ = optimize_sequence_ga(targets)
        result = mission_cost_function(list(best_seq), targets, dv_matrix, time_matrix)
        return list(best_seq), [(list(best_seq), result)]


def plot_mission_optimization(targets, best_seq, all_results, dv_matrix,
                              filename="mission_optimization.png"):
    """Plot mission sequence optimization results."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Cost distribution
    costs = [r[1]['total_cost'] / 1e6 for r in all_results]
    axes[0, 0].hist(costs, bins=30, color='steelblue', alpha=0.7, edgecolor='black')
    best_result = mission_cost_function(best_seq, targets, dv_matrix, 
                                        compute_cost_matrix(targets)[1])
    axes[0, 0].axvline(x=best_result['total_cost'] / 1e6, color='red', 
                       linestyle='--', linewidth=2, label=f'Best: ${best_result["total_cost"]/1e6:.1f}M')
    axes[0, 0].set_xlabel('Total Mission Cost ($M)')
    axes[0, 0].set_ylabel('Count')
    axes[0, 0].set_title('Mission Cost Distribution (All Sequences)')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Cost breakdown
    labels = ['Launch', 'Fuel', 'Operations']
    values = [best_result['launch_cost']/1e6, best_result['fuel_cost']/1e6, 
              best_result['time_cost']/1e6]
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']
    axes[0, 1].pie(values, labels=labels, colors=colors, autopct='%1.1f%%',
                   startangle=90, textprops={'fontsize': 10})
    axes[0, 1].set_title(f'Cost Breakdown (Total: ${best_result["total_cost"]/1e6:.1f}M)')
    
    # ΔV budget
    dvs = [dv_matrix[best_seq[i], best_seq[i+1]] for i in range(len(best_seq)-1)]
    labels = [f'{targets[best_seq[i]]["id"]}→\n{targets[best_seq[i+1]]["id"]}' 
              for i in range(len(best_seq)-1)]
    axes[1, 0].bar(range(len(dvs)), dvs, color='coral')
    axes[1, 0].set_xticks(range(len(dvs)))
    axes[1, 0].set_xticklabels(labels, fontsize=8)
    axes[1, 0].set_ylabel('ΔV (m/s)')
    axes[1, 0].set_title(f'ΔV Budget (Total: {sum(dvs):.0f} m/s)')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Mission timeline
    times = []
    cumulative = 0
    for i in range(len(best_seq) - 1):
        t_transfer = compute_cost_matrix(targets)[1][best_seq[i], best_seq[i+1]] / 86400
        times.append(('Transfer', cumulative, t_transfer))
        cumulative += t_transfer
        times.append(('Capture', cumulative, 4/24))  # 4 hours capture
        cumulative += 4/24
        times.append(('Deorbit', cumulative, 1))  # 1 day deorbit
        cumulative += 1
    
    colors_tl = {'Transfer': '#4ECDC4', 'Capture': '#FF6B6B', 'Deorbit': '#45B7D1'}
    for phase, start, dur in times:
        axes[1, 1].barh(0, dur, left=start, height=0.5, 
                        color=colors_tl[phase], alpha=0.8)
    axes[1, 1].set_xlabel('Mission Time (days)')
    axes[1, 1].set_title(f'Mission Timeline (Total: {cumulative:.0f} days)')
    axes[1, 1].set_yticks([])
    # Legend
    for phase, color in colors_tl.items():
        axes[1, 1].barh([], [], color=color, label=phase)
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved {filename}")


def plot_orbital_overview(targets, best_seq, filename="orbital_overview.png"):
    """Plot 3D orbital overview of the mission."""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Draw Earth
    u = np.linspace(0, 2*np.pi, 50)
    v = np.linspace(0, np.pi, 30)
    xe = R_EARTH/1e3 * np.outer(np.cos(u), np.sin(v))
    ye = R_EARTH/1e3 * np.outer(np.sin(u), np.sin(v))
    ze = R_EARTH/1e3 * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(xe, ye, ze, alpha=0.3, color='blue')
    
    # Draw orbits for selected targets
    colors = plt.cm.Set1(np.linspace(0, 1, len(best_seq)))
    theta = np.linspace(0, 2*np.pi, 200)
    
    for idx, ti in enumerate(best_seq):
        t = targets[ti]
        a_km = t['a'] / 1e3
        inc_r = np.radians(t['inc'])
        raan_r = np.radians(t['raan'])
        
        # Orbital positions
        x = a_km * np.cos(theta)
        y = a_km * np.sin(theta)
        z = np.zeros_like(theta)
        
        # Rotation matrices
        Rz1 = np.array([[np.cos(raan_r), -np.sin(raan_r), 0],
                        [np.sin(raan_r), np.cos(raan_r), 0],
                        [0, 0, 1]])
        Rx = np.array([[1, 0, 0],
                       [0, np.cos(inc_r), -np.sin(inc_r)],
                       [0, np.sin(inc_r), np.cos(inc_r)]])
        
        for j in range(len(theta)):
            pos = Rz1 @ Rx @ np.array([x[j], y[j], z[j]])
            x[j], y[j], z[j] = pos
        
        ax.plot(x, y, z, color=colors[idx], alpha=0.7, linewidth=1.5,
               label=f'{idx+1}: {t["id"]} ({t["altitude_km"]:.0f} km)')
    
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_title('ADR Mission - Orbital Overview')
    ax.legend(fontsize=8, loc='upper left')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, filename), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved {filename}")


# ============================================================
# Main Execution
# ============================================================

def main():
    print("=" * 70)
    print("  Active Debris Removal (ADR) Mission Design System")
    print("=" * 70)
    
    results = {}
    
    # ---- Module 1: Debris Catalog & Scoring ----
    print("\n[1/6] Debris Catalog Generation & Target Scoring...")
    catalog = generate_debris_catalog(n=30)
    catalog = compute_debris_scores(catalog)
    selected = select_targets(catalog, n_targets=5)
    
    print(f"  Generated {len(catalog)} debris objects")
    print(f"  Selected top-{len(selected)} targets:")
    for d in selected:
        print(f"    {d['id']}: alt={d['altitude_km']:.0f}km, mass={d['mass']:.0f}kg, "
              f"score={d['score']:.4f}")
    
    plot_scoring(catalog, selected)
    results['scoring'] = {
        'n_catalog': len(catalog),
        'n_selected': len(selected),
        'selected_ids': [d['id'] for d in selected],
        'selected_scores': [d['score'] for d in selected],
    }
    
    # ---- Module 2: Trajectory Optimization ----
    print("\n[2/6] Multi-Target Trajectory Optimization (GA)...")
    best_seq, fitness_history, dv_matrix, time_matrix = optimize_sequence_ga(selected)
    
    total_dv = sum(dv_matrix[best_seq[i], best_seq[i+1]] for i in range(len(best_seq)-1))
    print(f"  Optimal sequence: {[selected[i]['id'] for i in best_seq]}")
    print(f"  Total ΔV: {total_dv:.1f} m/s")
    print(f"  GA converged to cost: {fitness_history[-1]:.2f}")
    
    plot_trajectory_optimization(selected, best_seq, fitness_history, dv_matrix)
    plot_low_thrust_transfer(selected, best_seq)
    plot_orbital_overview(selected, best_seq)
    
    results['trajectory'] = {
        'best_sequence': [selected[i]['id'] for i in best_seq],
        'total_dv_ms': total_dv,
        'ga_final_cost': fitness_history[-1],
        'ga_generations': len(fitness_history),
        'dv_per_leg': [dv_matrix[best_seq[i], best_seq[i+1]] for i in range(len(best_seq)-1)],
    }
    
    # ---- Module 3: RPO Simulation ----
    print("\n[3/6] Rendezvous & Proximity Operations (Hill Equations)...")
    target_debris = selected[0]
    n_orbit = target_debris['mean_motion']
    
    x0 = [100, -500, 20]        # Initial relative position (m)
    v0 = [0, 0.5, 0]            # Initial relative velocity (m/s)
    xf = [0, 0, 0]              # Target: origin
    vf = [0, 0, 0]              # Target: zero relative velocity
    t_maneuver = 2 * np.pi / n_orbit * 0.75  # 3/4 orbit
    
    dv1, dv2, v0_req = plan_rpo_approach(n_orbit, x0, v0, xf, vf, t_maneuver)
    t_eval, states = simulate_rpo(n_orbit, x0, list(v0_req), t_maneuver, dt_out=5.0)
    
    print(f"  Initial separation: {np.linalg.norm(x0):.0f} m")
    print(f"  ΔV₁ (departure): {np.linalg.norm(dv1):.3f} m/s")
    print(f"  ΔV₂ (arrival): {np.linalg.norm(dv2):.3f} m/s")
    print(f"  Total RPO ΔV: {np.linalg.norm(dv1)+np.linalg.norm(dv2):.3f} m/s")
    print(f"  Maneuver time: {t_maneuver/60:.1f} min")
    
    plot_rpo_results(t_eval, states, dv1, dv2, n_orbit)
    
    results['rpo'] = {
        'initial_range_m': np.linalg.norm(x0),
        'dv1_ms': np.linalg.norm(dv1),
        'dv2_ms': np.linalg.norm(dv2),
        'total_rpo_dv_ms': np.linalg.norm(dv1) + np.linalg.norm(dv2),
        'maneuver_time_min': t_maneuver / 60,
        'final_range_m': np.linalg.norm(states[-1, :3]),
    }
    
    # ---- Module 4: Tumbling Debris ----
    print("\n[4/6] Tumbling Debris Attitude Estimation...")
    I_body = [500, 800, 300]    # kg⋅m² (asymmetric body)
    w0 = [0.05, 0.1, 0.03]     # rad/s initial tumble rates
    tumble_duration = 120       # seconds
    
    t_tumble, tumble_states = simulate_tumbling(I_body, w0, tumble_duration)
    w_measured, w_estimated, est_errors = estimate_tumbling_rate(t_tumble, tumble_states)
    
    mean_rate = np.mean(np.linalg.norm(tumble_states[:, 4:], axis=1))
    mean_error = np.mean(est_errors)
    print(f"  Inertia tensor: diag({I_body}) kg⋅m²")
    print(f"  Mean tumble rate: {np.degrees(mean_rate):.2f} deg/s")
    print(f"  Mean estimation error: {np.degrees(mean_error):.4f} deg/s")
    
    plot_tumbling(t_tumble, tumble_states, w_measured, w_estimated, est_errors)
    
    results['tumbling'] = {
        'inertia_kgm2': I_body,
        'initial_rates_rads': w0,
        'mean_tumble_rate_degs': np.degrees(mean_rate),
        'mean_estimation_error_degs': np.degrees(mean_error),
    }
    
    # ---- Module 5: Capture Mechanisms ----
    print("\n[5/6] Capture Mechanism Dynamics...")
    
    # Robotic arm
    print("  Simulating robotic arm...")
    arm_data = simulate_robotic_arm(
        joint_angles_init=[0.5, -0.3, 0.2],
        target_pos=np.array([3.5, 2.0])
    )
    
    # Net capture
    print("  Simulating net capture...")
    net_data = simulate_net_capture()
    net_history, net_captured, net_cap_time = net_data
    
    # Harpoon capture
    print("  Simulating harpoon capture...")
    harpoon_data = simulate_harpoon()
    h_hist, h_penetrated, h_pen_time = harpoon_data
    
    print(f"  Robotic arm: reached target in {len(arm_data[1])*0.1:.1f}s")
    print(f"  Net capture: {'Success' if net_captured else 'Failed'} "
          f"(t={net_cap_time:.2f}s)" if net_captured else "")
    print(f"  Harpoon: {'Impact' if h_penetrated else 'Missed'} "
          f"(t={h_pen_time:.3f}s)" if h_penetrated else "")
    
    plot_capture_mechanisms(arm_data, net_data, harpoon_data)
    
    results['capture'] = {
        'robotic_arm_time_s': len(arm_data[1]) * 0.1,
        'net_captured': net_captured,
        'net_capture_time_s': net_cap_time,
        'harpoon_impact': h_penetrated,
        'harpoon_impact_time_s': h_pen_time,
    }
    
    # ---- Module 6: Mission Cost Optimization ----
    print("\n[6/6] Mission Sequence Cost Optimization...")
    best_cost_seq, all_results_cost = optimize_mission_cost(
        selected, dv_matrix, time_matrix)
    best_cost_result = mission_cost_function(best_cost_seq, selected, dv_matrix, time_matrix)
    
    print(f"  Optimal sequence: {[selected[i]['id'] for i in best_cost_seq]}")
    print(f"  Total mission cost: ${best_cost_result['total_cost']/1e6:.1f}M")
    print(f"  Total ΔV: {best_cost_result['total_dv']:.0f} m/s")
    print(f"  Total fuel: {best_cost_result['total_fuel_kg']:.1f} kg")
    print(f"  Mission duration: {best_cost_result['total_time_days']:.0f} days")
    
    plot_mission_optimization(selected, best_cost_seq, all_results_cost, dv_matrix)
    
    results['mission_cost'] = {
        'optimal_sequence': [selected[i]['id'] for i in best_cost_seq],
        'total_cost_M': best_cost_result['total_cost'] / 1e6,
        'total_dv_ms': best_cost_result['total_dv'],
        'total_fuel_kg': best_cost_result['total_fuel_kg'],
        'mission_duration_days': best_cost_result['total_time_days'],
        'cost_breakdown': {
            'launch_M': best_cost_result['launch_cost'] / 1e6,
            'fuel_M': best_cost_result['fuel_cost'] / 1e6,
            'operations_M': best_cost_result['time_cost'] / 1e6,
        },
        'n_sequences_evaluated': len(all_results_cost),
    }
    
    # Save results
    with open('simulation_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to simulation_results.json")
    
    print("\n" + "=" * 70)
    print("  Simulation Complete!")
    print("=" * 70)
    
    return results


if __name__ == '__main__':
    main()
