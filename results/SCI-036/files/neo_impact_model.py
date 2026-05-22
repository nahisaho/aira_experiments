"""
NEO Impact Energy & Damage Estimation Module
Coupled atmospheric entry / ground damage model for NEO collision scenarios.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass
from scipy.integrate import odeint
from typing import Optional


@dataclass
class ImpactScenario:
    """Parameters defining a potential NEO impact."""
    diameter_km: float          # Object diameter [km]
    density_kg_m3: float        # Bulk density [kg/m^3]
    velocity_km_s: float        # Impact velocity [km/s]
    entry_angle_deg: float      # Entry angle from horizontal [deg]
    target_type: str            # 'land' or 'ocean'
    latitude_deg: float = 45.0  # Impact latitude (for ejecta model)


@dataclass
class DamageEstimate:
    """Estimated damage zones from an impact."""
    kinetic_energy_Mt: float         # Total kinetic energy [Mt TNT]
    crater_diameter_km: float        # Final crater diameter [km]
    blast_radius_km: dict            # Radii for different overpressure levels
    thermal_radius_km: float         # Thermal radiation 3rd-degree burn radius [km]
    ejecta_radius_km: float          # Significant ejecta blanket radius [km]
    tsunami_wave_height_m: Optional[float]   # For ocean impacts [m]
    affected_area_km2: float         # Total affected area [km^2]
    casualties_estimate: dict        # Low/median/high fatality estimates
    torino_scale: int
    palermo_scale: float


class AtmosphericEntryModel:
    """
    Meteoroid/asteroid atmospheric entry model.
    Integrates the equations of motion through the atmosphere.
    Based on Chyba et al. (1993) and Collins et al. (2005).
    """
    ATMO_SCALE_HEIGHT = 8.0e3     # m
    SURFACE_DENSITY = 1.2         # kg/m^3
    DRAG_COEFF = 0.9
    LIFT_COEFF = 0.0
    ABLATION_COEFF = 1e-9         # kg/J (low ablation for rocky body)
    G = 9.81                      # m/s^2

    def __init__(self, scenario: ImpactScenario):
        self.s = scenario
        self.radius_m = scenario.diameter_km * 500.0
        self.mass_kg = scenario.density_kg_m3 * (4/3) * np.pi * self.radius_m**3
        self.area_m2 = np.pi * self.radius_m**2
        self.v0 = scenario.velocity_km_s * 1e3
        self.theta = np.radians(scenario.entry_angle_deg)

    def _air_density(self, altitude_m: float) -> float:
        return self.SURFACE_DENSITY * np.exp(-altitude_m / self.ATMO_SCALE_HEIGHT)

    def _equations_of_motion(self, state: list, t: float) -> list:
        """ODE: [v, h, m] → [dv/dt, dh/dt, dm/dt]."""
        v, h, m = state
        if h < 0:
            return [0, 0, 0]
        rho = self._air_density(max(h, 0))
        area = self.area_m2 * (self.mass_kg / m if m > 0 else 1)**(2/3)
        # Drag deceleration
        dvdt = -0.5 * rho * v**2 * self.DRAG_COEFF * area / max(m, 1) + self.G * np.sin(self.theta)
        # Vertical velocity component
        dhdt = -v * np.sin(self.theta)
        # Mass ablation
        dmdt = -self.ABLATION_COEFF * 0.5 * rho * v**3 * area
        return [dvdt, dhdt, dmdt]

    def integrate(self) -> dict:
        """Integrate entry trajectory from 100 km altitude to surface."""
        h0 = 100e3  # 100 km
        t = np.linspace(0, 200.0, 20000)
        state0 = [self.v0, h0, self.mass_kg]

        with np.errstate(over='ignore', invalid='ignore'):
            sol = odeint(self._equations_of_motion, state0, t,
                         rtol=1e-6, atol=1e-8, full_output=False)

        v_traj = sol[:, 0]
        h_traj = sol[:, 1]
        m_traj = sol[:, 2]

        # Find surface impact point
        impact_idx = np.argmax(h_traj < 0)
        if impact_idx == 0:
            impact_idx = len(t) - 1

        # Check for airburst (Mach < 1 or fragmentation)
        airburst_altitude = None
        ram_pressure = 0.5 * self._air_density(max(h_traj[impact_idx], 0)) * v_traj[impact_idx]**2
        # Crude fragmentation: when dynamic pressure > compressive strength (~10 MPa)
        for i in range(len(t)):
            if h_traj[i] > 0:
                q = 0.5 * self._air_density(h_traj[i]) * v_traj[i]**2
                if q > 1e7:  # 10 MPa
                    airburst_altitude = h_traj[i]
                    break

        v_impact = v_traj[impact_idx]
        m_impact = max(m_traj[impact_idx], self.mass_kg * 0.01)   # floor at 1% remaining

        return {
            'v_impact_km_s': v_impact / 1e3,
            'm_impact_kg': m_impact,
            'airburst_altitude_m': airburst_altitude,
            'v_trajectory_km_s': v_traj / 1e3,
            'h_trajectory_km': h_traj / 1e3,
            't': t,
        }


class ImpactDamageModel:
    """
    Empirical scaling laws for crater size and damage radii.
    Combines: Collins et al. (2005) ImpactEarth model;
              Holsapple (1993) crater scaling; Glasstone & Dolan blast model.
    """
    EARTH_RADIUS_KM = 6371.0

    def __init__(self, scenario: ImpactScenario):
        self.s = scenario

    def _kinetic_energy_Mt(self, v_km_s: float, mass_kg: float) -> float:
        KE_J = 0.5 * mass_kg * (v_km_s * 1e3)**2
        return KE_J / 4.184e15

    def _crater_diameter_km(self, KE_Mt: float, target_type: str) -> float:
        """Holsapple (1993) / Collins scaling for crater diameter."""
        # Simple power-law: D [km] ~ 0.074 * (E [Mt])^0.294 (land)
        if target_type == 'land':
            return 0.074 * KE_Mt**0.294
        else:
            return 0.09 * KE_Mt**0.294   # ocean: slightly larger transient

    def _blast_radii(self, KE_Mt: float) -> dict:
        """
        Blast wave radii for different overpressures.
        From Glasstone & Dolan nuclear blast scaling (adapted for asteroid).
        r [km] ~ r0 * E^(1/3)
        """
        E_kT = KE_Mt * 1000  # kilotons
        return {
            'total_destruction_1_4_bar':  0.28 * E_kT**(1/3),
            'severe_damage_0_6_bar':      0.56 * E_kT**(1/3),
            'moderate_damage_0_3_bar':    0.90 * E_kT**(1/3),
            'window_breakage_0_007_bar':  5.50 * E_kT**(1/3),
        }

    def _thermal_radius_km(self, KE_Mt: float) -> float:
        """Radius for 3rd-degree burns (~125 kJ/m^2)."""
        # Thermal energy ~ 35% of total; E_th = 0.35 * KE; 4π r^2 * 125e3 = E_th
        E_J = KE_Mt * 4.184e15
        E_thermal = 0.35 * E_J
        r_m = np.sqrt(E_thermal / (4 * np.pi * 125e3))
        return r_m / 1e3

    def _tsunami_height_km_km(self, KE_Mt: float, depth_m: float = 4000.0) -> float:
        """
        Coastal tsunami wave height at 1000 km for ocean impact.
        Simple energy scaling.
        """
        # Wave height ∝ (E / depth / R^2)^0.5 — very crude
        H = 0.3 * (KE_Mt / 1000)**0.5  # meters at 500 km
        return H

    def _casualties(self, blast_radii: dict, thermal_radius: float) -> dict:
        """Very rough casualty estimates based on affected areas."""
        # World average population density ~60 / km^2 (land)
        rho_pop_land = 60.0
        r_sev = blast_radii['severe_damage_0_6_bar']
        r_mod = blast_radii['moderate_damage_0_3_bar']
        r_thm = thermal_radius

        area_severe_km2 = np.pi * r_sev**2
        area_moderate_km2 = np.pi * r_mod**2 - area_severe_km2
        area_thermal_km2 = np.pi * r_thm**2

        fatalities_low    = int(area_severe_km2 * rho_pop_land * 0.5)
        fatalities_median = int(area_severe_km2 * rho_pop_land * 1.0 +
                                 area_moderate_km2 * rho_pop_land * 0.1)
        fatalities_high   = int((area_severe_km2 + area_moderate_km2) * rho_pop_land * 1.5)

        return {
            'low': fatalities_low,
            'median': fatalities_median,
            'high': fatalities_high,
            'total_affected_km2': area_moderate_km2 + area_thermal_km2,
        }

    def estimate_damage(self, v_km_s: Optional[float] = None) -> DamageEstimate:
        """Compute full damage estimate for the scenario."""
        v = v_km_s or self.s.velocity_km_s
        R_m = self.s.diameter_km * 500
        mass = self.s.density_kg_m3 * (4/3) * np.pi * R_m**3

        KE_Mt = self._kinetic_energy_Mt(v, mass)
        crater_diam = self._crater_diameter_km(KE_Mt, self.s.target_type)
        blast = self._blast_radii(KE_Mt)
        thermal = self._thermal_radius_km(KE_Mt)
        ejecta_r = crater_diam * 3.0   # ejecta blanket extends ~3× crater radius
        tsunami = self._tsunami_height_km_km(KE_Mt) if self.s.target_type == 'ocean' else None
        casualties = self._casualties(blast, thermal)

        # Palermo / Torino (use 1e-4 default prob if standalone)
        f_bg = 0.03 * KE_Mt**(-0.8) if KE_Mt > 0 else 1e-4
        palermo = np.log10(1e-4 / (f_bg * 50))

        if KE_Mt < 1:
            torino = 1
        elif KE_Mt < 100:
            torino = 4
        elif KE_Mt < 1e5:
            torino = 7
        else:
            torino = 10

        return DamageEstimate(
            kinetic_energy_Mt=KE_Mt,
            crater_diameter_km=crater_diam,
            blast_radius_km=blast,
            thermal_radius_km=thermal,
            ejecta_radius_km=ejecta_r,
            tsunami_wave_height_m=tsunami,
            affected_area_km2=casualties['total_affected_km2'],
            casualties_estimate=casualties,
            torino_scale=torino,
            palermo_scale=palermo,
        )


def plot_damage_map(damage: DamageEstimate, scenario: ImpactScenario, save_path: str) -> None:
    """Visualize damage zones as concentric circles and energy breakdown."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    # Damage zone map
    ax = axes[0]
    radii_dict = {
        'Total destruction (1.4 bar)': damage.blast_radius_km['total_destruction_1_4_bar'],
        'Severe damage (0.6 bar)':     damage.blast_radius_km['severe_damage_0_6_bar'],
        'Moderate damage (0.3 bar)':   damage.blast_radius_km['moderate_damage_0_3_bar'],
        'Thermal 3rd-degree burns':    damage.thermal_radius_km,
        'Ejecta blanket':              damage.ejecta_radius_km,
        'Window breakage (0.007 bar)': damage.blast_radius_km['window_breakage_0_007_bar'],
    }
    colors = ['#8B0000', '#DC143C', '#FF6347', '#FF8C00', '#FFD700', '#ADFF2F']

    sorted_radii = sorted(radii_dict.items(), key=lambda x: x[1], reverse=True)
    for (label, r), color in zip(sorted_radii, colors[::-1]):
        circle = plt.Circle((0, 0), r, color=color, alpha=0.3, label=f'{label}: {r:.1f} km')
        ax.add_patch(circle)

    if damage.crater_diameter_km > 0.001:
        crater = plt.Circle((0, 0), damage.crater_diameter_km / 2,
                              color='black', alpha=0.8, label=f'Crater: Ø{damage.crater_diameter_km:.2f} km')
        ax.add_patch(crater)

    max_r = max(r for _, r in radii_dict.items())
    ax.set_xlim(-max_r * 1.1, max_r * 1.1)
    ax.set_ylim(-max_r * 1.1, max_r * 1.1)
    ax.set_aspect('equal')
    ax.set_xlabel('Distance [km]')
    ax.set_ylabel('Distance [km]')
    ax.set_title(f'Impact Damage Zones\n(D={scenario.diameter_km} km, '
                  f'v={scenario.velocity_km_s} km/s, E={damage.kinetic_energy_Mt:.1f} Mt)')
    ax.legend(loc='upper right', fontsize=8)

    # Energy and casualties bar chart
    ax2 = axes[1]
    categories = ['Low', 'Median', 'High']
    values = [damage.casualties_estimate['low'],
               damage.casualties_estimate['median'],
               damage.casualties_estimate['high']]
    bar_colors = ['#90EE90', '#FFA500', '#DC143C']
    bars = ax2.bar(categories, values, color=bar_colors, edgecolor='white')
    ax2.set_ylabel('Estimated Fatalities')
    ax2.set_title(f'Casualty Estimates\n(Torino Scale: {damage.torino_scale}, '
                   f'Palermo: {damage.palermo_scale:.1f})')
    ax2.set_yscale('log' if max(values) > 0 else 'linear')
    for bar, v in zip(bars, values):
        ax2.text(bar.get_x() + bar.get_width() / 2, max(v * 1.1, 1),
                  f'{v:,}', ha='center', va='bottom', fontsize=10)

    # Add energy annotation
    ax2.text(0.05, 0.95, f'KE = {damage.kinetic_energy_Mt:.1f} Mt\n'
                           f'Crater Ø = {damage.crater_diameter_km:.2f} km\n'
                           f'Thermal r = {damage.thermal_radius_km:.1f} km',
              transform=ax2.transAxes, fontsize=9, va='top',
              bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[Plot] Saved: {save_path}")


def compare_scenarios(scenarios: list[ImpactScenario], save_path: str) -> list[DamageEstimate]:
    """Compare multiple impact scenarios (diameter classes)."""
    damages = []
    for s in scenarios:
        m = ImpactDamageModel(s)
        damages.append(m.estimate_damage())

    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    diameters = [s.diameter_km for s in scenarios]
    KEs = [d.kinetic_energy_Mt for d in damages]
    blasts = [d.blast_radius_km['severe_damage_0_6_bar'] for d in damages]
    fatalities = [d.casualties_estimate['median'] for d in damages]

    axes[0].loglog(diameters, KEs, 'o-', color='crimson', markersize=8)
    axes[0].set_xlabel('Diameter [km]')
    axes[0].set_ylabel('Kinetic Energy [Mt TNT]')
    axes[0].set_title('Impact Energy vs. Size')
    axes[0].grid(True, alpha=0.4)

    axes[1].loglog(diameters, blasts, 's-', color='steelblue', markersize=8)
    axes[1].set_xlabel('Diameter [km]')
    axes[1].set_ylabel('Severe Damage Radius [km]')
    axes[1].set_title('Blast Damage Radius vs. Size')
    axes[1].grid(True, alpha=0.4)

    axes[2].loglog(diameters, [max(f, 1) for f in fatalities], '^-', color='seagreen', markersize=8)
    axes[2].set_xlabel('Diameter [km]')
    axes[2].set_ylabel('Estimated Fatalities (median)')
    axes[2].set_title('Casualty Estimate vs. Size')
    axes[2].grid(True, alpha=0.4)

    plt.suptitle('NEO Impact Scenario Comparison', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"[Plot] Saved: {save_path}")
    return damages
