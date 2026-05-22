"""
NEO Gravitational Perturbation & Yarkovsky Effect Module
Models planetary perturbations and non-gravitational forces acting on NEOs.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass
import rebound
from typing import Optional


@dataclass
class YarkovskyParameters:
    """Parameters for the Yarkovsky effect model."""
    diameter: float         # Object diameter [m]
    density: float          # Bulk density [kg/m^3]
    albedo: float           # Bond albedo
    emissivity: float = 0.9
    thermal_conductivity: float = 0.01   # W/m/K  (regolith-like)
    heat_capacity: float = 680.0          # J/kg/K
    rotation_period: float = 6.0 * 3600  # s (6 hours)
    obliquity: float = 0.0                # rad (prograde spin → max drift)

    @property
    def thermal_inertia(self) -> float:
        """Thermal inertia Γ = sqrt(κ ρ C) [J m^-2 K^-1 s^-0.5]."""
        return np.sqrt(self.thermal_conductivity * self.density * self.heat_capacity)

    @property
    def radius(self) -> float:
        return self.diameter / 2.0


class YarkovskyModel:
    """
    Diurnal Yarkovsky effect: transverse acceleration due to
    asymmetric thermal re-radiation from a rotating body.

    Reference: Vokrouhlický et al. (2000), Farnocchia et al. (2013).
    """
    AU = 1.496e11       # m
    AU_PER_YR = 0.0     # will be computed
    L_SUN = 3.828e26    # W
    STEFAN = 5.670e-8   # W/m^2/K^4
    c = 3e8             # m/s

    def __init__(self, params: YarkovskyParameters):
        self.p = params

    def _thermal_parameter(self, heliocentric_dist_au: float) -> float:
        """Thermal parameter Θ at given distance."""
        F = self.L_SUN / (4 * np.pi * (heliocentric_dist_au * self.AU)**2)
        T_sub = ((1 - self.p.albedo) * F / (self.p.emissivity * self.STEFAN))**0.25
        omega_rot = 2 * np.pi / self.p.rotation_period
        Gamma = self.p.thermal_inertia
        return Gamma * np.sqrt(omega_rot) / (self.p.emissivity * self.STEFAN * T_sub**3)

    def da_dt_au_per_yr(self, a_au: float) -> float:
        """
        Mean semimajor axis drift da/dt [AU/yr] from diurnal Yarkovsky.
        Sign depends on obliquity: prograde→positive (outward), retrograde→negative.
        """
        r = self.p.radius       # m
        rho = self.p.density    # kg/m^3
        F = self.L_SUN / (4 * np.pi * (a_au * self.AU)**2)  # W/m^2
        Theta = self._thermal_parameter(a_au)

        # Diurnal Yarkovsky acceleration (transverse component)
        kappa = 4 * self.p.emissivity * self.STEFAN
        T_eq = ((1 - self.p.albedo) * F / kappa)**0.25
        a_yark = -(8 * self.p.albedo * F) / (9 * self.c * rho * r) * \
                  Theta / (1 + Theta + 0.5 * Theta**2) * np.cos(self.p.obliquity)

        # Convert transverse acceleration [m/s^2] → da/dt [AU/yr]
        v_orb = np.sqrt(1.327e20 / (a_au * self.AU))  # m/s
        da_dt_m_s = 2 * a_yark / v_orb * (a_au * self.AU)   # rough adiabatic approx
        yr_to_s = 3.156e7
        da_dt_au_yr = da_dt_m_s * yr_to_s / self.AU
        return da_dt_au_yr

    def cumulative_drift(self, a0_au: float, t_years: float) -> float:
        """Total semimajor axis drift after t_years [AU]."""
        return self.da_dt_au_per_yr(a0_au) * t_years


class GravitationalPerturbationAnalyzer:
    """
    Analyzes secular and resonant gravitational perturbations using
    REBOUND orbit integration with planet ephemerides.
    """

    def __init__(self, neo, yarkovsky_params: Optional[YarkovskyParameters] = None):
        self.neo = neo
        self.yp = yarkovsky_params
        self.yark_model = YarkovskyModel(yarkovsky_params) if yarkovsky_params else None

    def _build_sim_with_yarkovsky(self) -> rebound.Simulation:
        """Build REBOUND simulation with Yarkovsky as extra force."""
        sim = rebound.Simulation()
        sim.integrator = "ias15"
        sim.units = ('yr', 'AU', 'Msun')
        sim.add(m=1.0)  # Sun

        # Major planets
        planet_data = [
            (9.55e-4, 5.203, 0.0489, np.radians(1.30), np.radians(273.9), np.radians(100.5), np.radians(20.0)),
            (2.86e-4, 9.537, 0.0539, np.radians(2.49), np.radians(339.4), np.radians(113.7), np.radians(317.0)),
            (4.37e-5, 19.19, 0.0473, np.radians(0.77), np.radians(98.0),  np.radians(74.0),  np.radians(142.2)),
            (5.15e-5, 30.07, 0.0086, np.radians(1.77), np.radians(276.3), np.radians(131.8), np.radians(256.2)),
        ]
        for pm, a, e, inc, omega, Omega, M in planet_data:
            sim.add(m=pm, a=a, e=e, inc=inc, omega=omega, Omega=Omega, M=M, primary=sim.particles[0])

        # Earth (separate for perturbation tracking)
        sim.add(m=3.003e-6, a=1.0, e=0.0167, inc=0.0,
                omega=np.radians(114.2), Omega=0.0, M=np.radians(357.5),
                primary=sim.particles[0])

        # NEO
        sim.add(m=0.0, a=self.neo.a, e=self.neo.e, inc=self.neo.inc,
                omega=self.neo.omega, Omega=self.neo.Omega, M=self.neo.M,
                primary=sim.particles[0])
        sim.move_to_com()
        return sim

    def analyze_secular_evolution(self, t_years: float = 200.0,
                                   n_outputs: int = 500) -> dict:
        """
        Track orbital element evolution and detect resonances/perturbations.
        """
        print(f"[Perturbation] Secular evolution over {t_years} yr")
        sim = self._build_sim_with_yarkovsky()
        neo_idx = len(sim.particles) - 1
        earth_idx = neo_idx - 1

        times = np.linspace(0, t_years, n_outputs)
        a_arr = np.zeros(n_outputs)
        e_arr = np.zeros(n_outputs)
        inc_arr = np.zeros(n_outputs)
        dist_earth = np.zeros(n_outputs)

        # Yarkovsky cumulative drift
        yark_da = np.zeros(n_outputs)

        for i, t in enumerate(times):
            sim.integrate(t)
            neo_p = sim.particles[neo_idx]
            orb = neo_p.orbit(primary=sim.particles[0])
            a_arr[i] = orb.a
            e_arr[i] = orb.e
            inc_arr[i] = orb.inc

            # Earth distance
            ep = sim.particles[earth_idx]
            dist_earth[i] = np.sqrt((neo_p.x - ep.x)**2 +
                                     (neo_p.y - ep.y)**2 +
                                     (neo_p.z - ep.z)**2)

            # Add Yarkovsky drift analytically
            if self.yark_model:
                yark_da[i] = self.yark_model.cumulative_drift(self.neo.a, t)

        # Compute resonance proximity: check n:m resonances with Jupiter
        # Jupiter period ≈ 11.86 yr, so check if a_neo corresponds to integer ratio
        T_jup = 11.86
        resonances_found = []
        for n in range(1, 7):
            for m in range(1, 7):
                a_res = (n / m * T_jup)**(2/3)  # Kepler's 3rd law
                if abs(self.neo.a - a_res) < 0.05:
                    resonances_found.append(f"{n}:{m} at {a_res:.3f} AU")

        results = {
            'times': times,
            'a': a_arr,
            'e': e_arr,
            'inc': inc_arr,
            'dist_earth_au': dist_earth,
            'yarkovsky_da': yark_da,
            'resonances_found': resonances_found,
        }
        if self.yark_model:
            da_total = self.yark_model.cumulative_drift(self.neo.a, t_years)
            results['yarkovsky_total_drift_au'] = da_total
            results['yarkovsky_rate_au_per_yr'] = self.yark_model.da_dt_au_per_yr(self.neo.a)
        print(f"[Perturbation] Done. Resonances: {resonances_found or 'none detected'}")
        return results

    def plot_perturbations(self, results: dict, save_path: str) -> None:
        """Visualize orbital evolution with perturbation contributions."""
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        t = results['times']

        axes[0, 0].plot(t, results['a'], color='steelblue', linewidth=0.8)
        if 'yarkovsky_da' in results:
            a_yark = results['a'][0] + results['yarkovsky_da']
            axes[0, 0].plot(t, a_yark, color='orange', linewidth=1.2,
                           linestyle='--', label='Yarkovsky drift only')
            axes[0, 0].legend(fontsize=9)
        axes[0, 0].set_xlabel('Time [yr]')
        axes[0, 0].set_ylabel('Semi-major axis [AU]')
        axes[0, 0].set_title('Orbital Semi-major Axis Evolution')

        axes[0, 1].plot(t, results['e'], color='crimson', linewidth=0.8)
        axes[0, 1].set_xlabel('Time [yr]')
        axes[0, 1].set_ylabel('Eccentricity')
        axes[0, 1].set_title('Eccentricity Evolution')

        axes[1, 0].plot(t, np.degrees(results['inc']), color='seagreen', linewidth=0.8)
        axes[1, 0].set_xlabel('Time [yr]')
        axes[1, 0].set_ylabel('Inclination [deg]')
        axes[1, 0].set_title('Inclination Evolution')

        axes[1, 1].semilogy(t, results['dist_earth_au'] / 0.00257,
                            color='purple', linewidth=0.8)
        axes[1, 1].axhline(1.0, color='gold', linestyle='--', label='1 lunar distance')
        axes[1, 1].set_xlabel('Time [yr]')
        axes[1, 1].set_ylabel('Earth Distance [LD]')
        axes[1, 1].set_title('NEO-Earth Distance')
        axes[1, 1].legend()

        plt.suptitle('Gravitational Perturbation Analysis', fontsize=13, fontweight='bold')
        plt.tight_layout()
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"[Plot] Saved: {save_path}")


def compute_yarkovsky_uncertainty_band(yp: YarkovskyParameters,
                                       a0_au: float, t_years: float,
                                       n_samples: int = 1000,
                                       seed: int = 42) -> dict:
    """
    Monte Carlo uncertainty on Yarkovsky drift due to
    thermal parameter uncertainties.
    """
    rng = np.random.default_rng(seed)
    # Assume ±30% uncertainty on thermal inertia and ±20% on density
    tc_samples = rng.normal(yp.thermal_conductivity, 0.3 * yp.thermal_conductivity, n_samples)
    rho_samples = rng.normal(yp.density, 0.2 * yp.density, n_samples)
    obliquity_samples = rng.uniform(0, np.pi, n_samples)

    drifts = []
    for tc, rho, obl in zip(tc_samples, rho_samples, obliquity_samples):
        yp_s = YarkovskyParameters(
            diameter=yp.diameter, density=max(100, rho),
            albedo=yp.albedo, thermal_conductivity=max(0.001, tc),
            obliquity=obl
        )
        m = YarkovskyModel(yp_s)
        drifts.append(m.cumulative_drift(a0_au, t_years))

    drifts = np.array(drifts)
    return {
        'mean': np.mean(drifts),
        'std': np.std(drifts),
        'p5': np.percentile(drifts, 5),
        'p95': np.percentile(drifts, 95),
        'samples': drifts,
    }
