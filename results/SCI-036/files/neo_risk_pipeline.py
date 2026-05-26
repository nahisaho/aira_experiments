#!/usr/bin/env python3
"""
NEO Impact Risk Assessment Pipeline
====================================
Bayesian impact probability evaluation system using REBOUND N-body integrator.

Modules:
1. Monte Carlo orbital uncertainty propagation
2. Gravitational perturbations + Yarkovsky effect
3. Keyhole systematic search
4. Bayesian probability update
5. Impact energy / damage estimation
6. DART/Hera-type deflection simulation
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.special import gamma as gamma_func
import rebound
import os
import json
import time

np.random.seed(42)
FIGURES_DIR = "figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================
# Constants
# ============================================================
AU = 1.496e11        # m
YEAR = 365.25        # days
G_SI = 6.674e-11     # m^3 kg^-1 s^-2
M_SUN = 1.989e30     # kg
M_EARTH = 5.972e24   # kg
R_EARTH = 6.371e6    # m
RHO_ROCK = 2500.0    # kg/m^3 (typical asteroid density)
C_LIGHT = 3e8        # m/s


# ============================================================
# 1. Monte Carlo Orbital Uncertainty Propagation
# ============================================================
def create_neo_simulation(a=1.1, e=0.3, inc=0.1, omega=0.5, Omega=1.0, f=0.0):
    """Create a REBOUND simulation with Sun, Earth, and a NEO."""
    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')
    sim.integrator = "ias15"
    sim.dt = 0.001

    # Sun
    sim.add(m=1.0)
    # Earth (simplified circular orbit)
    sim.add(m=3.003e-6, a=1.0, e=0.0167, inc=0.0, omega=1.993, Omega=0.0, f=1.75)
    # Jupiter
    sim.add(m=9.543e-4, a=5.2026, e=0.0489, inc=0.0227, omega=4.779, Omega=1.755, f=0.6)
    # NEO
    sim.add(m=0.0, a=a, e=e, inc=inc, omega=omega, Omega=Omega, f=f)

    sim.move_to_com()
    return sim


def monte_carlo_propagation(nominal_elements, cov_matrix, n_clones=500, t_end=10.0):
    """
    Propagate orbital uncertainty via Monte Carlo sampling.
    Returns minimum distances to Earth for each clone.
    """
    a0, e0, inc0, omega0, Omega0, f0 = nominal_elements

    # Sample clones from multivariate normal
    samples = np.random.multivariate_normal(nominal_elements, cov_matrix, size=n_clones)
    # Ensure physical validity
    samples[:, 0] = np.clip(samples[:, 0], 0.5, 3.0)   # a in AU
    samples[:, 1] = np.clip(samples[:, 1], 0.01, 0.99)  # e
    samples[:, 2] = np.clip(samples[:, 2], 0.0, np.pi)  # inc

    min_distances = []
    final_positions = []
    n_steps = 200

    for i, s in enumerate(samples):
        try:
            sim = create_neo_simulation(*s)
            times = np.linspace(0, t_end, n_steps)
            min_d = np.inf
            for t in times:
                sim.integrate(t)
                ps = sim.particles
                dx = ps[3].x - ps[1].x
                dy = ps[3].y - ps[1].y
                dz = ps[3].z - ps[1].z
                d = np.sqrt(dx**2 + dy**2 + dz**2)
                if d < min_d:
                    min_d = d
            min_distances.append(min_d)
            final_positions.append([ps[3].x, ps[3].y])
        except Exception:
            min_distances.append(np.inf)
            final_positions.append([0, 0])

    return np.array(min_distances), np.array(final_positions), samples


# ============================================================
# 2. Yarkovsky Effect Model
# ============================================================
def yarkovsky_acceleration(r_helio_AU, diameter_m, albedo=0.1, thermal_inertia=200,
                            spin_obliquity=0.0, rotation_period=6*3600):
    """
    Compute Yarkovsky acceleration magnitude (diurnal component).
    Uses the linearized model from Vokrouhlicky et al.
    """
    r_m = r_helio_AU * AU
    L_sun = 3.828e26  # W
    flux = L_sun / (4 * np.pi * r_m**2)

    # Thermal parameter
    omega_rot = 2 * np.pi / rotation_period
    Theta = thermal_inertia * np.sqrt(omega_rot) / (4 * 5.67e-8 *
            ((1 - albedo) * flux / (4 * 5.67e-8))**0.75)

    # Yarkovsky acceleration (diurnal)
    mass = RHO_ROCK * (4/3) * np.pi * (diameter_m/2)**3
    area = np.pi * (diameter_m/2)**2
    F_rad = flux * area * (1 - albedo) / C_LIGHT

    # Simplified seasonal + diurnal
    a_yarko = (F_rad / mass) * np.cos(spin_obliquity) * Theta / (1 + Theta + 0.5 * Theta**2)
    return a_yarko


def propagate_with_yarkovsky(nominal_elements, diameter_m=300, n_clones=200, t_end=10.0):
    """
    Propagate with Yarkovsky effect as stochastic perturbation.
    """
    cov = np.diag([1e-6, 1e-5, 1e-5, 1e-4, 1e-4, 1e-3])
    samples = np.random.multivariate_normal(nominal_elements, cov, size=n_clones)
    samples[:, 0] = np.clip(samples[:, 0], 0.5, 3.0)
    samples[:, 1] = np.clip(samples[:, 1], 0.01, 0.99)
    samples[:, 2] = np.clip(samples[:, 2], 0.0, np.pi)

    # Yarkovsky da/dt samples (uncertainty in thermal properties)
    diameters = np.random.normal(diameter_m, diameter_m * 0.1, n_clones)
    obliquities = np.random.uniform(0, np.pi, n_clones)

    semi_major_drift = []
    for i in range(n_clones):
        a_yarko = yarkovsky_acceleration(samples[i, 0], max(diameters[i], 10),
                                          spin_obliquity=obliquities[i])
        # Convert to da/dt in AU/yr (approximate)
        da_dt = a_yarko * (YEAR * 86400)**2 / AU * 2 * samples[i, 0]
        semi_major_drift.append(da_dt)

    return np.array(semi_major_drift), samples


# ============================================================
# 3. Keyhole Search Algorithm
# ============================================================
def compute_b_plane(sim, neo_idx=3, earth_idx=1):
    """Compute b-plane coordinates (xi, zeta) for close approach."""
    ps = sim.particles
    # Relative position and velocity
    dx = ps[neo_idx].x - ps[earth_idx].x
    dy = ps[neo_idx].y - ps[earth_idx].y
    dz = ps[neo_idx].z - ps[earth_idx].z
    dvx = ps[neo_idx].vx - ps[earth_idx].vx
    dvy = ps[neo_idx].vy - ps[earth_idx].vy
    dvz = ps[neo_idx].vz - ps[earth_idx].vz

    r = np.array([dx, dy, dz])
    v = np.array([dvx, dvy, dvz])

    v_hat = v / np.linalg.norm(v)
    # B-vector (impact parameter vector)
    b_vec = r - np.dot(r, v_hat) * v_hat
    xi = np.dot(b_vec, np.array([1, 0, 0]))
    zeta = np.dot(b_vec, np.array([0, 0, 1]))

    return xi, zeta, np.linalg.norm(b_vec)


def systematic_keyhole_search(nominal_elements, n_scan=300, t_end=15.0):
    """
    Systematically search for keyholes in the b-plane.
    Scan over initial true anomaly and semi-major axis offsets.
    """
    a0, e0, inc0, omega0, Omega0, f0 = nominal_elements

    da_range = np.linspace(-5e-4, 5e-4, n_scan)
    keyholes = []
    b_plane_data = []

    for da in da_range:
        try:
            sim = create_neo_simulation(a0 + da, e0, inc0, omega0, Omega0, f0)
            # Integrate to find close approach
            min_d = np.inf
            min_t = 0
            n_steps = 150
            for t in np.linspace(0, t_end, n_steps):
                sim.integrate(t)
                ps = sim.particles
                dx = ps[3].x - ps[1].x
                dy = ps[3].y - ps[1].y
                dz = ps[3].z - ps[1].z
                d = np.sqrt(dx**2 + dy**2 + dz**2)
                if d < min_d:
                    min_d = d
                    min_t = t

            # Compute b-plane at closest approach
            sim2 = create_neo_simulation(a0 + da, e0, inc0, omega0, Omega0, f0)
            sim2.integrate(min_t)
            xi, zeta, b = compute_b_plane(sim2)
            b_plane_data.append((xi, zeta, b, da))

            # Check for keyhole: close approach within Earth Hill sphere
            r_hill_earth = 1.0 * (M_EARTH / (3 * M_SUN))**(1.0/3.0)  # AU
            if min_d < 3 * r_hill_earth:
                keyholes.append({
                    'da': da,
                    'min_distance_AU': min_d,
                    'b_plane_xi': xi,
                    'b_plane_zeta': zeta,
                    'b_norm': b,
                    'approach_time_yr': min_t
                })
        except Exception:
            continue

    return keyholes, np.array(b_plane_data)


# ============================================================
# 4. Bayesian Impact Probability Update
# ============================================================
def compute_impact_probability(min_distances, impact_threshold_AU=4.26e-5):
    """Compute impact probability from Monte Carlo samples.
    impact_threshold_AU ~ 1 Earth radius in AU.
    """
    n_impacts = np.sum(min_distances < impact_threshold_AU)
    return n_impacts / len(min_distances)


def bayesian_update(prior_prob, likelihood_ratio):
    """
    Bayesian update: P(impact|data) = P(data|impact)*P(impact) / P(data)
    Using odds form for numerical stability.
    """
    prior_odds = prior_prob / (1 - prior_prob + 1e-30)
    posterior_odds = prior_odds * likelihood_ratio
    posterior_prob = posterior_odds / (1 + posterior_odds)
    return posterior_prob


def sequential_bayesian_update(initial_prob, observations_lr):
    """
    Sequential Bayesian update as new observations arrive.
    observations_lr: list of likelihood ratios from new data.
    """
    probs = [initial_prob]
    current = initial_prob
    for lr in observations_lr:
        current = bayesian_update(current, lr)
        probs.append(current)
    return np.array(probs)


# ============================================================
# 5. Impact Energy & Damage Estimation
# ============================================================
def impact_energy_megatons(diameter_m, velocity_kms=20.0, density=RHO_ROCK):
    """Compute impact energy in megatons TNT."""
    mass = density * (4/3) * np.pi * (diameter_m/2)**3
    v = velocity_kms * 1e3
    KE = 0.5 * mass * v**2
    MT_TNT = 4.184e15  # Joules per megaton TNT
    return KE / MT_TNT


def damage_radius_km(energy_MT):
    """
    Estimate damage radius using scaling law (Glasstone & Dolan).
    Returns: (fireball_km, overpressure_4psi_km, thermal_km)
    """
    # Scaling relations for airbursts/ground impacts
    fireball = 0.35 * energy_MT**0.4        # km
    overpressure = 4.7 * energy_MT**(1/3)   # 4 psi overpressure, km
    thermal = 8.0 * energy_MT**(1/3)        # thermal radiation, km
    return fireball, overpressure, thermal


def tsunami_wave_height(energy_MT, distance_km=100):
    """Estimate tsunami wave height for ocean impact."""
    # Ward & Asphaug (2000) scaling
    h = 10 * (energy_MT / 1000)**0.54 * (100 / max(distance_km, 1))**1.0
    return h  # meters


# ============================================================
# 6. DART/Hera Deflection Simulation
# ============================================================
def simulate_deflection(nominal_elements, spacecraft_mass_kg=610,
                         impact_velocity_kms=6.6, beta=3.61,
                         asteroid_diameter_m=160, deflection_time_yr=5.0,
                         t_total=15.0):
    """
    Simulate kinetic impactor deflection (DART-type).
    beta: momentum enhancement factor
    """
    a0, e0, inc0, omega0, Omega0, f0 = nominal_elements

    # Asteroid mass
    asteroid_mass = RHO_ROCK * (4/3) * np.pi * (asteroid_diameter_m/2)**3
    # Delta-v from impactor
    v_impact = impact_velocity_kms * 1e3
    delta_v = beta * spacecraft_mass_kg * v_impact / asteroid_mass
    delta_v_AU_yr = delta_v * YEAR * 86400 / AU

    # Undeflected trajectory
    sim_no = create_neo_simulation(*nominal_elements)
    min_d_no = []
    times = np.linspace(0, t_total, 300)
    for t in times:
        sim_no.integrate(t)
        ps = sim_no.particles
        dx = ps[3].x - ps[1].x
        dy = ps[3].y - ps[1].y
        dz = ps[3].z - ps[1].z
        d = np.sqrt(dx**2 + dy**2 + dz**2)
        min_d_no.append(d)

    # Deflected trajectory - apply delta-v at deflection_time
    sim_def = create_neo_simulation(*nominal_elements)
    min_d_def = []
    deflection_applied = False
    for t in times:
        sim_def.integrate(t)
        if t >= deflection_time_yr and not deflection_applied:
            ps = sim_def.particles
            # Apply delta-v along velocity direction
            v_mag = np.sqrt(ps[3].vx**2 + ps[3].vy**2 + ps[3].vz**2)
            ps[3].vx += delta_v_AU_yr * ps[3].vx / v_mag
            ps[3].vy += delta_v_AU_yr * ps[3].vy / v_mag
            ps[3].vz += delta_v_AU_yr * ps[3].vz / v_mag
            deflection_applied = True
        ps = sim_def.particles
        dx = ps[3].x - ps[1].x
        dy = ps[3].y - ps[1].y
        dz = ps[3].z - ps[1].z
        d = np.sqrt(dx**2 + dy**2 + dz**2)
        min_d_def.append(d)

    return times, np.array(min_d_no), np.array(min_d_def), delta_v, delta_v_AU_yr


# ============================================================
# MAIN: Run all experiments
# ============================================================
def main():
    results = {}
    print("=" * 60)
    print("NEO Impact Risk Assessment Pipeline")
    print("=" * 60)

    # Define nominal NEO (Apophis-like)
    nominal = [1.0924, 0.1912, np.radians(3.33), np.radians(126.4),
               np.radians(204.4), np.radians(45.0)]

    # --------------------------------------------------------
    # Experiment 1: Monte Carlo Orbital Propagation
    # --------------------------------------------------------
    print("\n[1/6] Monte Carlo Orbital Uncertainty Propagation...")
    cov = np.diag([2e-6, 5e-6, 1e-5, 5e-4, 5e-4, 1e-3])
    t0 = time.time()
    min_dists, final_pos, clones = monte_carlo_propagation(nominal, cov, n_clones=300, t_end=10.0)
    t1 = time.time()
    print(f"  Completed {len(min_dists)} clones in {t1-t0:.1f}s")
    print(f"  Min distance range: [{np.min(min_dists):.6f}, {np.max(min_dists):.6f}] AU")
    print(f"  Mean min distance: {np.mean(min_dists):.6f} AU")
    print(f"  Std min distance: {np.std(min_dists):.6f} AU")

    results['mc_propagation'] = {
        'n_clones': 300,
        'min_dist_mean': float(np.mean(min_dists)),
        'min_dist_std': float(np.std(min_dists)),
        'min_dist_min': float(np.min(min_dists)),
        'min_dist_max': float(np.max(min_dists)),
        'compute_time_s': t1-t0
    }

    # Plot 1: Minimum distance distribution
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].hist(min_dists, bins=40, color='steelblue', edgecolor='black', alpha=0.8, density=True)
    axes[0].axvline(4.26e-5, color='red', linestyle='--', linewidth=2, label='Earth radius (4.26e-5 AU)')
    axes[0].axvline(np.mean(min_dists), color='orange', linestyle='-', linewidth=2, label=f'Mean = {np.mean(min_dists):.5f} AU')
    axes[0].set_xlabel('Minimum Distance to Earth (AU)', fontsize=12)
    axes[0].set_ylabel('Probability Density', fontsize=12)
    axes[0].set_title('Monte Carlo: Minimum Approach Distance Distribution', fontsize=13)
    axes[0].legend(fontsize=10)
    axes[0].set_yscale('log')

    # Plot 2: Clone final positions
    axes[1].scatter(final_pos[:, 0], final_pos[:, 1], s=3, alpha=0.5, c=min_dists,
                   cmap='hot_r', label='NEO clones')
    axes[1].scatter([0], [0], s=100, c='yellow', edgecolors='orange', marker='*', zorder=5, label='Sun')
    theta_e = np.linspace(0, 2*np.pi, 100)
    axes[1].plot(np.cos(theta_e), np.sin(theta_e), 'b--', alpha=0.3, label='Earth orbit')
    axes[1].set_xlabel('x (AU)', fontsize=12)
    axes[1].set_ylabel('y (AU)', fontsize=12)
    axes[1].set_title('Clone Final Positions (color = min dist)', fontsize=13)
    axes[1].legend(fontsize=9)
    axes[1].set_aspect('equal')
    cb = plt.colorbar(axes[1].collections[0], ax=axes[1])
    cb.set_label('Min Distance (AU)')
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/mc_propagation.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> Saved figures/mc_propagation.png")

    # --------------------------------------------------------
    # Experiment 2: Yarkovsky Effect Analysis
    # --------------------------------------------------------
    print("\n[2/6] Yarkovsky Effect Modeling...")
    diameters = [50, 100, 200, 300, 500, 1000]
    yarko_accel = []
    yarko_da_dt = []
    for d in diameters:
        a_y = yarkovsky_acceleration(1.1, d)
        da = a_y * (YEAR * 86400)**2 / AU * 2 * 1.1
        yarko_accel.append(a_y)
        yarko_da_dt.append(da)
        print(f"  D={d}m: a_yarko={a_y:.3e} m/s^2, da/dt={da:.3e} AU/yr")

    drift_rates, yarko_samples = propagate_with_yarkovsky(nominal, diameter_m=300, n_clones=200)
    print(f"  Drift rate distribution: mean={np.mean(drift_rates):.3e}, std={np.std(drift_rates):.3e} AU/yr")

    results['yarkovsky'] = {
        'diameters': diameters,
        'accelerations': [float(a) for a in yarko_accel],
        'da_dt_AU_yr': [float(d) for d in yarko_da_dt],
        'drift_mean': float(np.mean(drift_rates)),
        'drift_std': float(np.std(drift_rates))
    }

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].loglog(diameters, yarko_accel, 'bo-', markersize=8, linewidth=2)
    axes[0].set_xlabel('Asteroid Diameter (m)', fontsize=12)
    axes[0].set_ylabel('Yarkovsky Acceleration (m/s²)', fontsize=12)
    axes[0].set_title('Yarkovsky Effect vs Asteroid Size', fontsize=13)
    axes[0].grid(True, alpha=0.3)

    axes[1].hist(drift_rates * 1e4, bins=30, color='coral', edgecolor='black', alpha=0.8, density=True)
    axes[1].set_xlabel('Semi-major Axis Drift Rate (×10⁻⁴ AU/yr)', fontsize=12)
    axes[1].set_ylabel('Probability Density', fontsize=12)
    axes[1].set_title('Yarkovsky da/dt Distribution (D=300m)', fontsize=13)
    axes[1].axvline(np.mean(drift_rates)*1e4, color='red', linestyle='--', label=f'Mean={np.mean(drift_rates)*1e4:.2f}')
    axes[1].legend()
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/yarkovsky_effect.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> Saved figures/yarkovsky_effect.png")

    # --------------------------------------------------------
    # Experiment 3: Keyhole Search
    # --------------------------------------------------------
    print("\n[3/6] Systematic Keyhole Search...")
    t0 = time.time()
    keyholes, b_plane = systematic_keyhole_search(nominal, n_scan=250, t_end=12.0)
    t1 = time.time()
    print(f"  Scanned 250 trajectory variants in {t1-t0:.1f}s")
    print(f"  Found {len(keyholes)} potential keyhole regions")
    for kh in keyholes[:5]:
        print(f"    da={kh['da']:.6f} AU, min_d={kh['min_distance_AU']:.6f} AU, "
              f"t={kh['approach_time_yr']:.2f} yr")

    results['keyholes'] = {
        'n_scanned': 250,
        'n_found': len(keyholes),
        'keyhole_details': keyholes[:10],
        'compute_time_s': t1-t0
    }

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    if len(b_plane) > 0:
        axes[0].scatter(b_plane[:, 0], b_plane[:, 1], s=5, c=b_plane[:, 2],
                       cmap='viridis', alpha=0.7)
        if len(keyholes) > 0:
            kh_xi = [k['b_plane_xi'] for k in keyholes]
            kh_zeta = [k['b_plane_zeta'] for k in keyholes]
            axes[0].scatter(kh_xi, kh_zeta, s=50, c='red', marker='x', linewidths=2,
                          label=f'Keyholes ({len(keyholes)})', zorder=5)
        theta = np.linspace(0, 2*np.pi, 100)
        r_hill = (M_EARTH / (3 * M_SUN))**(1.0/3.0)
        axes[0].plot(r_hill*np.cos(theta), r_hill*np.sin(theta), 'g--', alpha=0.5,
                    label='Earth Hill sphere')
        axes[0].set_xlabel('ξ (AU)', fontsize=12)
        axes[0].set_ylabel('ζ (AU)', fontsize=12)
        axes[0].set_title('B-plane Map & Keyhole Locations', fontsize=13)
        axes[0].legend(fontsize=9)

    if len(b_plane) > 0:
        axes[1].plot(b_plane[:, 3], b_plane[:, 2], 'b-', linewidth=1.5)
        axes[1].set_xlabel('Δa (AU)', fontsize=12)
        axes[1].set_ylabel('b-plane distance (AU)', fontsize=12)
        axes[1].set_title('Close Approach Distance vs Orbital Offset', fontsize=13)
        if len(keyholes) > 0:
            kh_da = [k['da'] for k in keyholes]
            kh_b = [k['b_norm'] for k in keyholes]
            axes[1].scatter(kh_da, kh_b, c='red', s=30, zorder=5, label='Keyholes')
            axes[1].legend()
        axes[1].grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/keyhole_search.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> Saved figures/keyhole_search.png")

    # --------------------------------------------------------
    # Experiment 4: Bayesian Probability Update
    # --------------------------------------------------------
    print("\n[4/6] Bayesian Impact Probability Update...")
    impact_threshold = 4.26e-5  # ~1 Earth radius in AU
    extended_threshold = 1e-3   # ~0.001 AU for "close approach"

    ip_strict = compute_impact_probability(min_dists, impact_threshold)
    ip_extended = compute_impact_probability(min_dists, extended_threshold)
    print(f"  Impact probability (strict, 1 R_Earth): {ip_strict:.6f}")
    print(f"  Close approach probability (<0.001 AU): {ip_extended:.6f}")

    # Simulate sequential Bayesian updates with new observations
    initial_prob = 1e-4  # Prior from initial orbit determination
    # Simulate likelihood ratios from successive observations
    np.random.seed(123)
    n_obs = 20
    # Gradually improving observations reduce uncertainty
    likelihood_ratios = np.exp(np.random.normal(-0.3, 0.5, n_obs))
    # Add a "recovery" observation that increases probability
    likelihood_ratios[7] = 3.5
    likelihood_ratios[12] = 0.1  # Strong constraint reduces probability

    posterior_probs = sequential_bayesian_update(initial_prob, likelihood_ratios)
    print(f"  Initial probability: {initial_prob:.2e}")
    print(f"  After {n_obs} observations: {posterior_probs[-1]:.2e}")

    results['bayesian'] = {
        'impact_prob_strict': float(ip_strict),
        'impact_prob_extended': float(ip_extended),
        'initial_prior': initial_prob,
        'final_posterior': float(posterior_probs[-1]),
        'n_observations': n_obs
    }

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    obs_indices = np.arange(len(posterior_probs))
    axes[0].semilogy(obs_indices, posterior_probs, 'b-o', markersize=5, linewidth=2)
    axes[0].axhline(initial_prob, color='gray', linestyle='--', alpha=0.5, label='Initial prior')
    axes[0].set_xlabel('Observation Number', fontsize=12)
    axes[0].set_ylabel('Impact Probability', fontsize=12)
    axes[0].set_title('Bayesian Sequential Update of Impact Probability', fontsize=13)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Also show Palermo/Torino scale
    # Torino scale: 0 (no hazard) to 10 (certain collision)
    # Simplified mapping
    torino = np.zeros_like(posterior_probs)
    for i, p in enumerate(posterior_probs):
        E = impact_energy_megatons(300) * p
        if E < 1e-6:
            torino[i] = 0
        elif E < 1e-3:
            torino[i] = 1
        elif E < 1:
            torino[i] = 2
        elif E < 100:
            torino[i] = 4
        else:
            torino[i] = 7

    # Palermo scale
    f_b = 0.03 * 300**(-2.6)  # background impact frequency (yr^-1)
    palermo = np.log10(posterior_probs / (f_b * 100 + 1e-30))

    axes[1].plot(obs_indices, palermo, 'r-s', markersize=5, linewidth=2)
    axes[1].axhline(0, color='black', linestyle='-', alpha=0.3, label='Palermo = 0 (background)')
    axes[1].axhline(-2, color='green', linestyle='--', alpha=0.5, label='Palermo = -2 (monitor)')
    axes[1].set_xlabel('Observation Number', fontsize=12)
    axes[1].set_ylabel('Palermo Scale', fontsize=12)
    axes[1].set_title('Palermo Technical Impact Hazard Scale', fontsize=13)
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/bayesian_update.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> Saved figures/bayesian_update.png")

    # --------------------------------------------------------
    # Experiment 5: Impact Energy & Damage Model
    # --------------------------------------------------------
    print("\n[5/6] Impact Energy & Damage Estimation...")
    diameters_damage = np.logspace(1, 4, 50)  # 10m to 10km
    velocities = [15, 20, 25, 30]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Energy vs diameter
    for v in velocities:
        energies = [impact_energy_megatons(d, v) for d in diameters_damage]
        axes[0, 0].loglog(diameters_damage, energies, linewidth=2, label=f'v={v} km/s')
    axes[0, 0].axhline(15, color='gray', linestyle=':', label='Hiroshima (15 kT)')
    axes[0, 0].axhline(50, color='orange', linestyle=':', label='Tunguska (~50 MT)')
    axes[0, 0].axhline(1e5, color='red', linestyle=':', label='Chicxulub (~10⁵ MT)')
    axes[0, 0].set_xlabel('Diameter (m)', fontsize=11)
    axes[0, 0].set_ylabel('Energy (MT TNT)', fontsize=11)
    axes[0, 0].set_title('Impact Energy vs Asteroid Diameter', fontsize=12)
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(True, alpha=0.3)

    # Damage radius vs diameter
    for v in velocities:
        fb_list, op_list, th_list = [], [], []
        for d in diameters_damage:
            E = impact_energy_megatons(d, v)
            fb, op, th = damage_radius_km(E)
            fb_list.append(fb)
            op_list.append(op)
            th_list.append(th)
        if v == 20:
            axes[0, 1].loglog(diameters_damage, op_list, 'r-', linewidth=2, label='Overpressure (4 psi)')
            axes[0, 1].loglog(diameters_damage, th_list, 'orange', linewidth=2, label='Thermal radiation')
            axes[0, 1].loglog(diameters_damage, fb_list, 'b-', linewidth=2, label='Fireball')
    axes[0, 1].set_xlabel('Diameter (m)', fontsize=11)
    axes[0, 1].set_ylabel('Damage Radius (km)', fontsize=11)
    axes[0, 1].set_title('Damage Radius vs Asteroid Size (v=20 km/s)', fontsize=12)
    axes[0, 1].legend(fontsize=9)
    axes[0, 1].grid(True, alpha=0.3)

    # Tsunami wave height
    distances_coast = np.logspace(1, 3, 50)
    for E_MT in [10, 100, 1000, 10000]:
        heights = [tsunami_wave_height(E_MT, d) for d in distances_coast]
        axes[1, 0].loglog(distances_coast, heights, linewidth=2, label=f'E={E_MT} MT')
    axes[1, 0].axhline(10, color='red', linestyle='--', alpha=0.5, label='Devastating (>10m)')
    axes[1, 0].set_xlabel('Distance from Impact (km)', fontsize=11)
    axes[1, 0].set_ylabel('Wave Height (m)', fontsize=11)
    axes[1, 0].set_title('Tsunami Wave Height (Ocean Impact)', fontsize=12)
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].grid(True, alpha=0.3)

    # Summary table as text
    table_data = []
    for d in [30, 50, 100, 200, 500, 1000, 5000, 10000]:
        E = impact_energy_megatons(d, 20)
        fb, op, th = damage_radius_km(E)
        freq = 1 / (0.03 * d**(-2.6) + 1e-30) if d > 0 else 0
        table_data.append([d, E, op, th])
        print(f"  D={d:>5d}m: E={E:>12.1f} MT, Overpressure={op:>8.1f} km, Thermal={th:>8.1f} km")

    axes[1, 1].axis('off')
    col_labels = ['Diameter\n(m)', 'Energy\n(MT)', 'Blast\n(km)', 'Thermal\n(km)']
    table = axes[1, 1].table(cellText=[[f'{r[0]}', f'{r[1]:.1f}', f'{r[2]:.1f}', f'{r[3]:.1f}']
                                        for r in table_data],
                             colLabels=col_labels, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    axes[1, 1].set_title('Impact Damage Summary (v=20 km/s)', fontsize=12, pad=20)

    results['damage'] = {
        'table': [{'diameter_m': r[0], 'energy_MT': r[1], 'blast_km': r[2], 'thermal_km': r[3]}
                  for r in table_data]
    }

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/impact_damage.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> Saved figures/impact_damage.png")

    # --------------------------------------------------------
    # Experiment 6: DART/Hera Deflection Simulation
    # --------------------------------------------------------
    print("\n[6/6] DART/Hera Deflection Simulation...")
    times, dist_no, dist_def, dv, dv_au = simulate_deflection(
        nominal, spacecraft_mass_kg=610, impact_velocity_kms=6.6,
        beta=3.61, asteroid_diameter_m=160, deflection_time_yr=5.0, t_total=15.0
    )
    print(f"  Delta-v applied: {dv:.6f} m/s ({dv_au:.3e} AU/yr)")
    print(f"  Undeflected min distance: {np.min(dist_no):.6f} AU")
    print(f"  Deflected min distance: {np.min(dist_def):.6f} AU")
    print(f"  Distance change at closest: {np.min(dist_def) - np.min(dist_no):.6f} AU")

    # Multi-beta comparison
    betas = [1.0, 2.0, 3.61, 5.0, 7.0]
    deflection_results = []
    for b in betas:
        _, _, d_def, dv_b, _ = simulate_deflection(
            nominal, beta=b, deflection_time_yr=5.0, t_total=15.0
        )
        deflection_results.append({
            'beta': b,
            'min_dist_AU': float(np.min(d_def)),
            'delta_v_ms': float(dv_b * AU / (YEAR * 86400))
        })
        print(f"  β={b:.2f}: min_d={np.min(d_def):.6f} AU, Δv={dv_b * AU / (YEAR * 86400):.4f} m/s")

    # Lead time comparison
    lead_times = [1, 2, 5, 10, 15, 20]
    lead_time_results = []
    for lt in lead_times:
        try:
            _, d_no_lt, d_def_lt, _, _ = simulate_deflection(
                nominal, beta=3.61, deflection_time_yr=float(lt),
                t_total=max(lt + 10, 20)
            )
            miss_diff = np.min(d_def_lt) - np.min(d_no_lt)
            lead_time_results.append({'lead_time_yr': lt, 'miss_distance_change_AU': float(miss_diff)})
        except Exception:
            lead_time_results.append({'lead_time_yr': lt, 'miss_distance_change_AU': 0.0})

    results['deflection'] = {
        'delta_v_ms': float(dv),
        'min_dist_undeflected': float(np.min(dist_no)),
        'min_dist_deflected': float(np.min(dist_def)),
        'beta_comparison': deflection_results,
        'lead_time_comparison': lead_time_results
    }

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Trajectory comparison
    axes[0, 0].plot(times, dist_no, 'r-', linewidth=2, label='Undeflected')
    axes[0, 0].plot(times, dist_def, 'b-', linewidth=2, label='Deflected (β=3.61)')
    axes[0, 0].axvline(5.0, color='green', linestyle='--', alpha=0.5, label='Deflection time')
    axes[0, 0].axhline(4.26e-5, color='gray', linestyle=':', alpha=0.5, label='1 R_Earth')
    axes[0, 0].set_xlabel('Time (years)', fontsize=11)
    axes[0, 0].set_ylabel('Distance to Earth (AU)', fontsize=11)
    axes[0, 0].set_title('Trajectory: Deflected vs Undeflected', fontsize=12)
    axes[0, 0].legend(fontsize=9)
    axes[0, 0].set_yscale('log')
    axes[0, 0].grid(True, alpha=0.3)

    # Beta comparison
    beta_vals = [r['beta'] for r in deflection_results]
    min_dists_beta = [r['min_dist_AU'] for r in deflection_results]
    axes[0, 1].bar(range(len(betas)), min_dists_beta, color=['#ff6b6b', '#ffa07a', '#98fb98', '#87ceeb', '#6495ed'],
                   edgecolor='black')
    axes[0, 1].set_xticks(range(len(betas)))
    axes[0, 1].set_xticklabels([f'β={b}' for b in betas])
    axes[0, 1].set_ylabel('Minimum Distance (AU)', fontsize=11)
    axes[0, 1].set_title('Effect of Momentum Enhancement (β)', fontsize=12)
    axes[0, 1].grid(True, alpha=0.3, axis='y')

    # Lead time effect
    lt_vals = [r['lead_time_yr'] for r in lead_time_results]
    md_change = [abs(r['miss_distance_change_AU']) for r in lead_time_results]
    axes[1, 0].plot(lt_vals, md_change, 'go-', markersize=8, linewidth=2)
    axes[1, 0].set_xlabel('Lead Time (years)', fontsize=11)
    axes[1, 0].set_ylabel('|ΔMiss Distance| (AU)', fontsize=11)
    axes[1, 0].set_title('Deflection Effectiveness vs Lead Time', fontsize=12)
    axes[1, 0].grid(True, alpha=0.3)

    # Delta-v vs asteroid size
    ast_diameters = np.logspace(1.5, 3.5, 30)
    dv_list = []
    for d in ast_diameters:
        m = RHO_ROCK * (4/3) * np.pi * (d/2)**3
        dv_val = 3.61 * 610 * 6600 / m
        dv_list.append(dv_val)
    axes[1, 1].loglog(ast_diameters, dv_list, 'purple', linewidth=2)
    axes[1, 1].axhline(0.01, color='red', linestyle='--', label='Typical deflection threshold')
    axes[1, 1].set_xlabel('Asteroid Diameter (m)', fontsize=11)
    axes[1, 1].set_ylabel('Delta-v (m/s)', fontsize=11)
    axes[1, 1].set_title('Achievable Δv vs Asteroid Size (DART-type)', fontsize=12)
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/deflection_simulation.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> Saved figures/deflection_simulation.png")

    # --------------------------------------------------------
    # Summary Figure: Pipeline Overview
    # --------------------------------------------------------
    print("\n[Summary] Generating pipeline overview figure...")
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.axis('off')

    # Pipeline flow diagram as text
    steps = [
        ("1. Orbit\nDetermination", "Astrometric\nobservations\n→ Orbital elements\n& covariance"),
        ("2. MC\nPropagation", f"N={results['mc_propagation']['n_clones']} clones\nMean min dist:\n{results['mc_propagation']['min_dist_mean']:.5f} AU"),
        ("3. Yarkovsky\nEffect", f"da/dt mean:\n{results['yarkovsky']['drift_mean']:.2e}\nAU/yr"),
        ("4. Keyhole\nSearch", f"Found {results['keyholes']['n_found']}\nkeyhole\nregions"),
        ("5. Bayesian\nUpdate", f"Prior: {results['bayesian']['initial_prior']:.1e}\nPosterior:\n{results['bayesian']['final_posterior']:.2e}"),
        ("6. Deflection\nAssessment", f"Δv={results['deflection']['delta_v_ms']:.4f} m/s\nβ=3.61")
    ]

    colors = ['#3498db', '#2ecc71', '#e74c3c', '#f39c12', '#9b59b6', '#1abc9c']
    for i, (title, desc) in enumerate(steps):
        x = 0.08 + i * 0.155
        rect = plt.Rectangle((x, 0.3), 0.13, 0.4, facecolor=colors[i], alpha=0.3,
                             edgecolor=colors[i], linewidth=2, transform=ax.transAxes)
        ax.add_patch(rect)
        ax.text(x + 0.065, 0.62, title, fontsize=11, fontweight='bold',
               ha='center', va='center', transform=ax.transAxes)
        ax.text(x + 0.065, 0.42, desc, fontsize=8,
               ha='center', va='center', transform=ax.transAxes)
        if i < 5:
            ax.annotate('', xy=(x + 0.145, 0.5), xytext=(x + 0.13, 0.5),
                       arrowprops=dict(arrowstyle='->', color='black', lw=2),
                       transform=ax.transAxes)

    ax.text(0.5, 0.85, 'NEO Impact Risk Assessment Pipeline', fontsize=16,
           fontweight='bold', ha='center', transform=ax.transAxes)
    ax.text(0.5, 0.15, 'Bayesian framework with REBOUND N-body integration',
           fontsize=12, ha='center', style='italic', transform=ax.transAxes)

    plt.savefig(f'{FIGURES_DIR}/pipeline_overview.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> Saved figures/pipeline_overview.png")

    # Save results to JSON
    with open('experiment_results.json', 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print("\n  -> Saved experiment_results.json")

    print("\n" + "=" * 60)
    print("All experiments completed successfully!")
    print("=" * 60)

    return results


if __name__ == '__main__':
    results = main()
