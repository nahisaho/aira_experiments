"""
NEO Orbital Propagation Module
Monte Carlo uncertainty propagation using REBOUND n-body integrator.
"""

import numpy as np
import rebound
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import Optional
import warnings
warnings.filterwarnings('ignore')


@dataclass
class OrbitalElements:
    """Keplerian orbital elements with uncertainties (1-sigma)."""
    a: float        # Semi-major axis [AU]
    e: float        # Eccentricity
    inc: float      # Inclination [rad]
    omega: float    # Argument of perihelion [rad]
    Omega: float    # Longitude of ascending node [rad]
    M: float        # Mean anomaly [rad]
    # Uncertainties (1-sigma)
    da: float = 1e-8
    de: float = 1e-8
    dinc: float = 1e-6
    domega: float = 1e-6
    dOmega: float = 1e-6
    dM: float = 1e-6
    # Physical properties
    diameter: float = 0.14    # km
    density: float = 1500.0   # kg/m^3
    H: float = 22.0           # Absolute magnitude


class MonteCarloOrbitalPropagator:
    """
    Monte Carlo orbital propagation with REBOUND n-body integrator.
    Propagates an ensemble of virtual asteroids to assess close approach statistics.
    """

    PLANETS = {
        'mercury': (3.302e23, 0.387098, 0.205630, np.radians(7.005), np.radians(29.124), np.radians(48.331), np.radians(174.796)),
        'venus':   (4.869e24, 0.723332, 0.006772, np.radians(3.395), np.radians(54.884), np.radians(76.680), np.radians(50.416)),
        'earth':   (5.972e24, 1.000000, 0.016708, np.radians(0.000), np.radians(114.208), np.radians(0.0), np.radians(357.517)),
        'mars':    (6.417e23, 1.523679, 0.093400, np.radians(1.850), np.radians(286.502), np.radians(49.558), np.radians(19.373)),
        'jupiter': (1.899e27, 5.202887, 0.048900, np.radians(1.303), np.radians(273.867), np.radians(100.464), np.radians(20.020)),
        'saturn':  (5.685e26, 9.536676, 0.053862, np.radians(2.489), np.radians(339.391), np.radians(113.665), np.radians(317.020)),
    }
    AU_TO_KM = 1.496e8
    EARTH_RADIUS_AU = 6.371e3 / 1.496e8
    HILL_SPHERE_EARTH_AU = 0.01  # ~1% of 1 AU

    def __init__(self, neo: OrbitalElements, n_clones: int = 500, seed: int = 42):
        self.neo = neo
        self.n_clones = n_clones
        self.rng = np.random.default_rng(seed)
        self.close_approach_data: list[dict] = []

    def _sample_orbital_elements(self) -> np.ndarray:
        """Sample n_clones orbital element sets from Gaussian distributions."""
        n = self.n_clones
        elements = np.zeros((n, 6))
        elements[:, 0] = self.rng.normal(self.neo.a,     self.neo.da,     n)
        elements[:, 1] = np.clip(self.rng.normal(self.neo.e, self.neo.de, n), 0.0, 0.999)
        elements[:, 2] = self.rng.normal(self.neo.inc,   self.neo.dinc,   n)
        elements[:, 3] = self.rng.normal(self.neo.omega, self.neo.domega, n)
        elements[:, 4] = self.rng.normal(self.neo.Omega, self.neo.dOmega, n)
        elements[:, 5] = self.rng.normal(self.neo.M,     self.neo.dM,     n)
        return elements

    def _build_simulation(self, orbital_elements: np.ndarray) -> rebound.Simulation:
        """Build a REBOUND simulation with Sun, planets, and NEO clones."""
        sim = rebound.Simulation()
        sim.integrator = "ias15"
        sim.units = ('yr', 'AU', 'Msun')
        sim.dt = 0.01  # years

        # Add Sun
        sim.add(m=1.0)

        # Add perturbing planets
        for name, (mass, a, e, inc, omega, Omega, M) in self.PLANETS.items():
            m_solar = mass / 1.989e30
            sim.add(m=m_solar, a=a, e=e, inc=inc, omega=omega, Omega=Omega, M=M, primary=sim.particles[0])

        # Add NEO clones
        for i in range(len(orbital_elements)):
            a, e, inc, omega, Omega, M = orbital_elements[i]
            sim.add(m=0.0, a=a, e=e, inc=inc, omega=omega, Omega=Omega, M=M, primary=sim.particles[0])

        sim.move_to_com()
        return sim

    def _detect_close_approaches(self, sim: rebound.Simulation,
                                  n_clones: int, t_yr: float) -> list[dict]:
        """Check minimum distances of all clones to Earth."""
        approaches = []
        earth = sim.particles[4]  # Earth is index 4 (after Sun + 3 inner planets)
        neo_start = len(self.PLANETS) + 1  # offset after Sun + planets

        for i in range(n_clones):
            neo_p = sim.particles[neo_start + i]
            dx = neo_p.x - earth.x
            dy = neo_p.y - earth.y
            dz = neo_p.z - earth.z
            dist_au = np.sqrt(dx**2 + dy**2 + dz**2)
            if dist_au < self.HILL_SPHERE_EARTH_AU:
                approaches.append({
                    'clone_id': i,
                    't_yr': t_yr,
                    'dist_au': dist_au,
                    'dist_ld': dist_au / 0.00257,  # lunar distances
                    'impact': dist_au < self.EARTH_RADIUS_AU * 1.1
                })
        return approaches

    def propagate(self, t_years: float = 100.0, n_outputs: int = 200) -> dict:
        """
        Propagate all clones for t_years and collect close approach statistics.
        Returns dict with impact probabilities and orbital dispersion.
        """
        print(f"[MC Propagation] {self.n_clones} clones × {t_years} yr integration")
        orbital_elements = self._sample_orbital_elements()
        sim = self._build_simulation(orbital_elements)

        impacts = []
        min_distances = np.full(self.n_clones, np.inf)
        times = np.linspace(0, t_years, n_outputs)

        for i, t in enumerate(times):
            sim.integrate(t)
            approaches = self._detect_close_approaches(sim, self.n_clones, t)
            for app in approaches:
                ci = app['clone_id']
                if app['dist_au'] < min_distances[ci]:
                    min_distances[ci] = app['dist_au']
                if app['impact']:
                    impacts.append(app)

            if (i + 1) % 50 == 0:
                print(f"  t = {t:.1f} yr  |  impacts so far: {len(impacts)}")

        n_impact = len(set(a['clone_id'] for a in impacts))
        impact_prob = n_impact / self.n_clones

        results = {
            'n_clones': self.n_clones,
            't_years': t_years,
            'n_impacts': n_impact,
            'impact_probability': impact_prob,
            'min_distances_au': min_distances,
            'close_approaches': impacts,
            'orbital_elements_samples': orbital_elements,
        }
        self.close_approach_data = impacts
        print(f"[MC Propagation] Done. P_impact = {impact_prob:.4e} ({n_impact}/{self.n_clones})")
        return results

    def plot_minimum_distances(self, results: dict, save_path: str) -> None:
        """Plot histogram of minimum Earth-approach distances."""
        dists = results['min_distances_au']
        dists_ld = dists / 0.00257

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Histogram
        finite_dists = dists_ld[dists_ld < 1000]
        axes[0].hist(finite_dists, bins=50, color='steelblue', edgecolor='white', alpha=0.8)
        axes[0].axvline(1.0, color='gold', linestyle='--', linewidth=2, label='Moon distance')
        axes[0].axvline(0.00257 / 0.00257, color='crimson', linestyle='--', linewidth=2, label='Earth radius limit')
        axes[0].set_xlabel('Minimum Distance [Lunar Distances]')
        axes[0].set_ylabel('Count')
        axes[0].set_title('Distribution of Minimum Earth Approach Distances')
        axes[0].legend()
        axes[0].set_yscale('log')

        # Scatter: a vs e of sampled clones
        elems = results['orbital_elements_samples']
        axes[1].scatter(elems[:, 0], elems[:, 1], alpha=0.3, s=8, c='steelblue', label='Virtual asteroids')
        impact_ids = list(set(a['clone_id'] for a in results['close_approaches']))
        if impact_ids:
            axes[1].scatter(elems[impact_ids, 0], elems[impact_ids, 1],
                            alpha=0.9, s=30, c='crimson', zorder=5, label=f'Impactors ({len(impact_ids)})')
        axes[1].set_xlabel('Semi-major axis [AU]')
        axes[1].set_ylabel('Eccentricity')
        axes[1].set_title('Orbital Element Dispersion of Virtual Asteroids')
        axes[1].legend()

        plt.tight_layout()
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"[Plot] Saved: {save_path}")


def demo_neo() -> OrbitalElements:
    """Return a representative Apophis-like NEO for demonstration."""
    return OrbitalElements(
        a=0.9226, e=0.1911, inc=np.radians(3.331),
        omega=np.radians(126.4), Omega=np.radians(204.4), M=np.radians(222.2),
        da=2e-7, de=2e-7, dinc=np.radians(0.001),
        domega=np.radians(0.01), dOmega=np.radians(0.01), dM=np.radians(0.01),
        diameter=0.37, density=2600.0, H=19.7
    )
