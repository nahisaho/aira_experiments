"""
DART/Hera-type Deflection Mission Simulation Module
Models kinetic impactor + gravity tractor deflection effectiveness.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Optional
import rebound


@dataclass
class NEOPhysicalParams:
    """NEO physical parameters for deflection modeling."""
    diameter_km: float
    density_kg_m3: float
    mass_kg: Optional[float] = None
    beta_min: float = 1.0    # Momentum enhancement factor (β) lower bound
    beta_max: float = 5.0    # Momentum enhancement factor upper bound (ejecta)
    beta_mean: float = 2.5   # DART measured β for Dimorphos ≈ 2.5-4.5

    def __post_init__(self):
        if self.mass_kg is None:
            R = self.diameter_km * 500
            self.mass_kg = self.density_kg_m3 * (4/3) * np.pi * R**3


@dataclass
class KineticImpactorMission:
    """Parameters defining a kinetic impactor mission."""
    spacecraft_mass_kg: float = 570.0   # DART: 570 kg
    impact_velocity_km_s: float = 6.14  # DART impact velocity
    lead_time_years: float = 10.0       # Years before Earth encounter
    mission_name: str = "DART-like"


@dataclass
class GravityTractorMission:
    """Parameters for gravity tractor deflection."""
    spacecraft_mass_kg: float = 1000.0
    hover_distance_m: float = 300.0   # Distance above surface
    thrust_N: float = 0.04            # 40 mN ion thruster
    operating_time_years: float = 5.0


@dataclass
class DeflectionResult:
    """Results from deflection simulation."""
    mission_type: str
    delta_v_imparted_m_s: float          # Δv imparted to NEO [m/s]
    b_plane_displacement_km: float       # B-plane miss distance change [km]
    miss_distance_change_km: float       # Earth miss distance change [km]
    deflection_success_probability: float
    beta_uncertainty_factor: float       # Spread due to β uncertainty
    notes: str = ""


class KineticImpactorSimulator:
    """
    Simulates kinetic impactor deflection using momentum transfer mechanics.
    Models momentum enhancement factor β and N-body trajectory modification.
    """

    def __init__(self, neo: NEOPhysicalParams, mission: KineticImpactorMission):
        self.neo = neo
        self.mission = mission

    def compute_delta_v(self, beta: float) -> float:
        """
        Δv imparted to NEO from kinetic impactor.
        Δv = β × m_sc × v_impact / m_neo

        β = 1: pure momentum transfer (no ejecta)
        β > 1: ejecta amplifies momentum transfer
        """
        p_impactor = self.mission.spacecraft_mass_kg * self.mission.impact_velocity_km_s * 1e3
        delta_v = beta * p_impactor / self.neo.mass_kg
        return delta_v   # m/s

    def compute_miss_distance(self, delta_v_m_s: float) -> float:
        """
        Approximate miss distance change from Δv applied at lead time T.
        d_miss ≈ Δv × T × (some orbital geometry factor)
        Full calculation requires orbit integration; here use linear approximation.
        """
        T_s = self.mission.lead_time_years * 3.156e7  # seconds
        # For a ~1 AU orbit, 1 m/s Δv × 10 yr ≈ ~3×10^8 m ≈ 300 km miss distance
        # More precise: use Valsecchi corridor width sensitivity
        d_miss_m = delta_v_m_s * T_s * 0.5    # factor 0.5 for geometry
        return d_miss_m / 1e3  # km

    def b_plane_displacement(self, delta_v_m_s: float,
                               encounter_geometry_factor: float = 1.0) -> float:
        """B-plane displacement [km] from applied Δv."""
        T_s = self.mission.lead_time_years * 3.156e7
        # Öpik linearized sensitivity
        disp = delta_v_m_s * T_s / 1e3 * encounter_geometry_factor
        return disp  # km

    def simulate_with_uncertainty(self, n_mc: int = 1000, seed: int = 42) -> dict:
        """
        Monte Carlo deflection simulation accounting for β uncertainty.
        """
        rng = np.random.default_rng(seed)
        # Sample β from truncated normal
        beta_mean = self.neo.beta_mean
        beta_std = (self.neo.beta_max - self.neo.beta_min) / 4.0
        betas = rng.normal(beta_mean, beta_std, n_mc)
        betas = np.clip(betas, self.neo.beta_min, self.neo.beta_max)

        delta_vs = np.array([self.compute_delta_v(b) for b in betas])
        miss_dists = np.array([self.compute_miss_distance(dv) for dv in delta_vs])
        b_disps = np.array([self.b_plane_displacement(dv) for dv in delta_vs])

        # Success: deflection moves NEO outside Earth's collision radius (6371 km)
        earth_radius_km = 6371.0
        success_threshold_km = earth_radius_km * 1.5  # need >1.5× Earth radius displacement
        success_prob = np.mean(b_disps > success_threshold_km)

        return {
            'betas': betas,
            'delta_vs_m_s': delta_vs,
            'miss_distances_km': miss_dists,
            'b_plane_displacements_km': b_disps,
            'mean_delta_v_m_s': np.mean(delta_vs),
            'std_delta_v_m_s': np.std(delta_vs),
            'mean_miss_km': np.mean(miss_dists),
            'success_probability': success_prob,
            'p10_miss_km': np.percentile(miss_dists, 10),
            'p90_miss_km': np.percentile(miss_dists, 90),
        }


class GravityTractorSimulator:
    """
    Gravity tractor: spacecraft hovering above NEO uses mutual gravity
    to slowly drag the asteroid off course.
    """
    G = 6.674e-11

    def __init__(self, neo: NEOPhysicalParams, mission: GravityTractorMission):
        self.neo = neo
        self.mission = mission

    def gravitational_force_N(self) -> float:
        """Mutual gravitational force between spacecraft and NEO."""
        r = self.mission.hover_distance_m + self.neo.diameter_km * 500
        F = self.G * self.mission.spacecraft_mass_kg * self.neo.mass_kg / r**2
        return F

    def net_thrust_efficiency(self) -> float:
        """Effective thrust considering off-axis gravity component."""
        theta = np.arctan2(self.mission.hover_distance_m,
                            self.neo.diameter_km * 500)
        return self.mission.thrust_N * np.sin(theta) / self.gravitational_force_N()

    def compute_delta_v(self) -> float:
        """Total Δv accumulated over mission duration."""
        t_s = self.mission.operating_time_years * 3.156e7
        F_net = self.gravitational_force_N() * self.net_thrust_efficiency()
        a_neo = F_net / self.neo.mass_kg
        return a_neo * t_s  # m/s


class HeraFollowOnSimulator:
    """
    Models Hera-type characterization + kinetic impactor precision improvement.
    Hera measures: mass, β, internal structure → reduces deflection uncertainty.
    """

    def __init__(self, neo: NEOPhysicalParams):
        self.neo = neo

    def pre_hera_uncertainty(self) -> dict:
        """Uncertainty in deflection outcome before Hera characterization."""
        return {
            'beta_uncertainty_1sigma': (self.neo.beta_max - self.neo.beta_min) / 4,
            'mass_uncertainty_fraction': 0.30,      # 30% mass uncertainty
            'delta_v_uncertainty_fraction': 0.35,
        }

    def post_hera_uncertainty(self) -> dict:
        """Uncertainty after Hera characterization mission."""
        return {
            'beta_uncertainty_1sigma': 0.2,          # Hera reduces β uncertainty to ~±0.2
            'mass_uncertainty_fraction': 0.05,       # 5% mass uncertainty
            'delta_v_uncertainty_fraction': 0.08,    # 8% Δv uncertainty
        }

    def compare(self) -> dict:
        pre = self.pre_hera_uncertainty()
        post = self.post_hera_uncertainty()
        return {
            'pre_hera': pre,
            'post_hera': post,
            'beta_improvement_factor': pre['beta_uncertainty_1sigma'] / post['beta_uncertainty_1sigma'],
            'mass_improvement_factor': pre['mass_uncertainty_fraction'] / post['mass_uncertainty_fraction'],
            'delta_v_improvement_factor': pre['delta_v_uncertainty_fraction'] / post['delta_v_uncertainty_fraction'],
        }


def simulate_full_deflection_campaign(neo: NEOPhysicalParams,
                                       kinetic: KineticImpactorMission,
                                       gt: GravityTractorMission,
                                       save_path: str) -> dict:
    """
    Full deflection scenario comparison: DART-like, Gravity Tractor, Combined.
    """
    dart_sim = KineticImpactorSimulator(neo, kinetic)
    dart_results = dart_sim.simulate_with_uncertainty()

    gt_sim = GravityTractorSimulator(neo, gt)
    gt_dv = gt_sim.compute_delta_v()
    gt_miss = dart_sim.compute_miss_distance(gt_dv)

    hera_sim = HeraFollowOnSimulator(neo)
    hera_comparison = hera_sim.compare()

    # Combined: kinetic + gravity tractor
    combined_dv_mean = dart_results['mean_delta_v_m_s'] + gt_dv

    print(f"[Deflection] DART-like: ΔV = {dart_results['mean_delta_v_m_s']:.4f} m/s, "
          f"miss = {dart_results['mean_miss_km']:.1f} km, "
          f"P_success = {dart_results['success_probability']:.3f}")
    print(f"[Deflection] Gravity Tractor: ΔV = {gt_dv:.4f} m/s, miss = {gt_miss:.1f} km")
    print(f"[Deflection] Combined ΔV = {combined_dv_mean:.4f} m/s")
    print(f"[Deflection] Hera β improvement: {hera_comparison['beta_improvement_factor']:.1f}×")

    # --- Plots ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. β distribution and Δv scatter
    ax = axes[0, 0]
    ax.hist(dart_results['betas'], bins=30, color='steelblue', edgecolor='white',
             density=True, alpha=0.7, label='β samples')
    ax.axvline(neo.beta_mean, color='crimson', linestyle='--', linewidth=2,
                label=f'β mean = {neo.beta_mean}')
    ax.set_xlabel('Momentum Enhancement Factor β')
    ax.set_ylabel('Density')
    ax.set_title('β Distribution (DART-like Impact)')
    ax.legend()

    # 2. Miss distance distribution
    ax = axes[0, 1]
    ax.hist(dart_results['b_plane_displacements_km'], bins=30, color='seagreen',
             edgecolor='white', density=True, alpha=0.7, label='B-plane displacement')
    ax.axvline(6371, color='crimson', linestyle='--', linewidth=2,
                label='Earth radius (6371 km)')
    ax.axvline(dart_results['mean_miss_km'], color='navy', linestyle='-', linewidth=1.5,
                label=f'Mean: {dart_results["mean_miss_km"]:.0f} km')
    ax.set_xlabel('B-plane Displacement [km]')
    ax.set_ylabel('Density')
    ax.set_title(f'Deflection Miss Distance Distribution\n'
                  f'P(success) = {dart_results["success_probability"]:.1%}')
    ax.legend()

    # 3. Mission comparison bar chart
    ax = axes[1, 0]
    mission_names = ['Kinetic\nImpactor', 'Gravity\nTractor', 'Combined\n(KI+GT)']
    delta_vs = [dart_results['mean_delta_v_m_s'], gt_dv, combined_dv_mean]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    bars = ax.bar(mission_names, delta_vs, color=colors, edgecolor='white', linewidth=1.5)
    for bar, dv in zip(bars, delta_vs):
        ax.text(bar.get_x() + bar.get_width() / 2, dv + max(delta_vs) * 0.01,
                 f'{dv:.4f}\nm/s', ha='center', va='bottom', fontsize=9)
    ax.set_ylabel('Δv Imparted to NEO [m/s]')
    ax.set_title('Deflection Mission Comparison: Δv')

    # 4. Pre vs post Hera uncertainty
    ax = axes[1, 1]
    metrics = ['β uncertainty\n(1σ)', 'Mass uncertainty\n(%)', 'Δv uncertainty\n(%)']
    pre_vals = [hera_comparison['pre_hera']['beta_uncertainty_1sigma'],
                 hera_comparison['pre_hera']['mass_uncertainty_fraction'] * 100,
                 hera_comparison['pre_hera']['delta_v_uncertainty_fraction'] * 100]
    post_vals = [hera_comparison['post_hera']['beta_uncertainty_1sigma'],
                  hera_comparison['post_hera']['mass_uncertainty_fraction'] * 100,
                  hera_comparison['post_hera']['delta_v_uncertainty_fraction'] * 100]
    x = np.arange(len(metrics))
    width = 0.35
    ax.bar(x - width/2, pre_vals, width, label='Pre-Hera', color='#d62728', alpha=0.8)
    ax.bar(x + width/2, post_vals, width, label='Post-Hera', color='#2ca02c', alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, fontsize=9)
    ax.set_ylabel('Uncertainty Value')
    ax.set_title('Hera Characterization: Uncertainty Reduction')
    ax.legend()

    plt.suptitle('DART/Hera-type Deflection Mission Simulation', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[Plot] Saved: {save_path}")

    return {
        'dart': dart_results,
        'gravity_tractor': {'delta_v_m_s': gt_dv, 'miss_km': gt_miss},
        'combined_dv_m_s': combined_dv_mean,
        'hera_improvement': hera_comparison,
    }
