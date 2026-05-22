"""
NEO Keyhole Analysis Module
Systematic search for keyholes (resonant return collision corridors) in b-plane.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle
from dataclasses import dataclass
from typing import Optional


@dataclass
class KeyholeResult:
    """Represents a found keyhole region."""
    resonance_n: int           # Orbit resonance numerator (returns in n orbits)
    resonance_m: int           # Orbit resonance denominator
    center_xi: float           # b-plane ξ coordinate [km]
    center_zeta: float         # b-plane ζ coordinate [km]
    half_width_km: float       # Keyhole half-width [km]
    impact_prob: float         # Probability of falling in this keyhole
    encounter_year: float      # Year of resonant return
    delta_v_deflect: float     # Required Δv for deflection [m/s]


class BPlaneAnalyzer:
    """
    Öpik-Valsecchi b-plane (close approach plane) analysis for keyhole search.
    The b-plane is perpendicular to the incoming asymptote; (ξ, ζ) coords.

    Reference: Valsecchi et al. (2003) A&A; Milani et al. (2005).
    """
    AU_TO_KM = 1.496e8
    EARTH_RADIUS_KM = 6371.0
    EARTH_MASS_KG = 5.972e24
    EARTH_MU = 3.986e14  # m^3/s^2

    def __init__(self, neo, uncertainty_sigma_km: float = 1000.0):
        """
        neo: OrbitalElements
        uncertainty_sigma_km: 1-sigma positional uncertainty in b-plane [km]
        """
        self.neo = neo
        self.sigma_km = uncertainty_sigma_km

    def _b_plane_ellipse(self, sigma_xi_km: float, sigma_zeta_km: float,
                          rho: float = 0.0) -> tuple:
        """Return (a, b, angle) of 3-sigma uncertainty ellipse in b-plane."""
        cov = np.array([[sigma_xi_km**2, rho * sigma_xi_km * sigma_zeta_km],
                         [rho * sigma_xi_km * sigma_zeta_km, sigma_zeta_km**2]])
        vals, vecs = np.linalg.eigh(cov)
        order = np.argsort(vals)[::-1]
        vals, vecs = vals[order], vecs[:, order]
        a = 3 * np.sqrt(vals[0])   # 3-sigma semi-major axis
        b = 3 * np.sqrt(vals[1])
        angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
        return a, b, angle

    def _keyhole_width(self, a_au: float, encounter_vel_kms: float,
                        resonance_n: int, resonance_m: int) -> float:
        """
        Estimate keyhole width for n:m resonant return [km].
        Analytic approximation from Valsecchi et al. (2003).
        """
        G = 6.674e-11
        M_earth = self.EARTH_MASS_KG
        v_inf = encounter_vel_kms * 1e3  # m/s
        mu_earth = G * M_earth
        b_earth = self.EARTH_RADIUS_KM * 1e3 * np.sqrt(1 + 2 * mu_earth / (v_inf**2 * self.EARTH_RADIUS_KM * 1e3))
        # Resonance condition: period ratio close to n/m
        T_orb = a_au**1.5  # years
        delta_a_resonance = a_au * (2/3) * (resonance_m / resonance_n - 1) / resonance_n
        # Width ∝ gravitational focusing × orbital resonance width
        width_km = b_earth * abs(delta_a_resonance) / a_au * 1e-3
        width_km = max(width_km, 1.0)   # physical floor
        return width_km

    def _probability_in_keyhole(self, center_xi: float, center_zeta: float,
                                  half_width_km: float,
                                  sigma_xi_km: float, sigma_zeta_km: float) -> float:
        """
        Fraction of Gaussian probability within a 1D stripe (keyhole) in b-plane.
        Approximates keyhole as a rectangle centered at (center_xi, center_zeta).
        """
        from scipy import stats
        if sigma_xi_km <= 0 or sigma_zeta_km <= 0 or half_width_km <= 0:
            return 0.0
        z_xi = abs(center_xi) / sigma_xi_km
        z_zeta_lo = (center_zeta - half_width_km) / sigma_zeta_km
        z_zeta_hi = (center_zeta + half_width_km) / sigma_zeta_km
        # Only compute if within ±6σ
        if z_xi > 6 or (abs(z_zeta_lo) > 6 and abs(z_zeta_hi) > 6):
            return 0.0
        prob_xi = stats.norm.pdf(center_xi, 0, sigma_xi_km)
        prob_zeta = stats.norm.cdf(z_zeta_hi) - stats.norm.cdf(z_zeta_lo)
        prob = prob_xi * prob_zeta * 2 * half_width_km
        if not np.isfinite(prob):
            return 0.0
        return min(prob, 1.0)

    def search_keyholes(self, encounter_vel_kms: float = 5.87,
                         n_resonances: int = 8,
                         b_plane_center: tuple = (0.0, 0.0),
                         sigma_ratio: float = 2.0) -> list[KeyholeResult]:
        """
        Systematic search for all keyholes within ±3σ of b-plane nominal.
        Returns list of KeyholeResult sorted by impact probability.
        """
        print(f"[Keyhole] Searching resonant return keyholes (v_∞={encounter_vel_kms} km/s)")
        sigma_xi = self.sigma_km
        sigma_zeta = self.sigma_km * sigma_ratio

        keyholes = []
        # Scan resonances n:m for n=1..n_resonances, m=1..n_resonances
        for n in range(1, n_resonances + 1):
            for m in range(1, n_resonances + 1):
                # Skip unphysical ratios
                if abs(n / m - self.neo.a**1.5) > 0.5:  # rough Kepler filter
                    continue
                # b-plane ζ position of resonant keyhole strip
                # (based on resonant semi-major axis, simplified geometry)
                a_res = (n / m)**(2/3)
                if abs(a_res - 1.0) > 2.0:  # filter unreasonable resonances
                    continue

                zeta_k = (a_res - 1.0) * 1e5 + b_plane_center[1]   # scale to km
                xi_k = b_plane_center[0] + self.sigma_km * np.sin(np.pi * n / n_resonances)

                hw = self._keyhole_width(self.neo.a, encounter_vel_kms, n, m)
                prob = self._probability_in_keyhole(xi_k, zeta_k, hw, sigma_xi, sigma_zeta)
                if prob < 1e-12:
                    continue

                # Estimate deflection Δv needed to miss keyhole
                # Δv ~ (keyhole_width / σ_b) × (v_inf) / sqrt(t_remaining)
                t_remaining = 10.0  # years until encounter (placeholder)
                dv_deflect = (hw / sigma_xi) * encounter_vel_kms * 1e3 / (3.156e7 * t_remaining) * 1e3  # m/s
                encounter_year = 2029.0 + n * self.neo.a**1.5 * m  # simplified

                kh = KeyholeResult(
                    resonance_n=n, resonance_m=m,
                    center_xi=xi_k, center_zeta=zeta_k,
                    half_width_km=hw,
                    impact_prob=prob,
                    encounter_year=encounter_year,
                    delta_v_deflect=dv_deflect,
                )
                keyholes.append(kh)

        keyholes.sort(key=lambda k: k.impact_prob, reverse=True)
        print(f"[Keyhole] Found {len(keyholes)} candidate keyholes")
        for kh in keyholes[:5]:
            print(f"  {kh.resonance_n}:{kh.resonance_m} | "
                  f"ζ={kh.center_zeta:.1f} km | w={kh.half_width_km:.2f} km | "
                  f"P={kh.impact_prob:.2e}")
        return keyholes

    def plot_b_plane(self, keyholes: list[KeyholeResult],
                      encounter_vel_kms: float,
                      save_path: str,
                      sigma_ratio: float = 2.0) -> None:
        """Visualize b-plane with uncertainty ellipse and keyhole locations."""
        fig, axes = plt.subplots(1, 2, figsize=(14, 7))

        sigma_xi = max(self.sigma_km, 1.0)
        sigma_zeta = max(self.sigma_km * sigma_ratio, 1.0)

        # Left: b-plane overview
        ax = axes[0]
        # 1σ, 2σ, 3σ uncertainty ellipses
        for s, alpha in [(1, 0.5), (2, 0.35), (3, 0.2)]:
            ell = Ellipse((0, 0), 2 * s * sigma_xi, 2 * s * sigma_zeta,
                          facecolor='steelblue', alpha=alpha, edgecolor='steelblue',
                          label=f'{s}σ uncertainty' if s == 3 else '')
            ax.add_patch(ell)

        # Earth cross-section (collision disk)
        G, M = 6.674e-11, 5.972e24
        v_inf = encounter_vel_kms * 1e3
        b_earth_km = self.EARTH_RADIUS_KM * np.sqrt(1 + 2 * G * M / (v_inf**2 * self.EARTH_RADIUS_KM * 1e3))
        earth_circ = Circle((0, 0), b_earth_km, facecolor='none',
                             edgecolor='crimson', linewidth=2, label=f'Earth disk ({b_earth_km:.0f} km)')
        ax.add_patch(earth_circ)

        # Plot keyholes
        top_kh = keyholes[:10]
        for kh in top_kh:
            ax.barh(kh.center_zeta, 2 * sigma_xi * 3, left=-3 * sigma_xi,
                    height=2 * kh.half_width_km, color='gold', alpha=0.3,
                    edgecolor='darkorange', linewidth=1.2)
            ax.text(sigma_xi * 3.1, kh.center_zeta,
                    f'{kh.resonance_n}:{kh.resonance_m}', fontsize=7, va='center', color='darkorange')

        ax.set_xlim(-4 * sigma_xi, 5 * sigma_xi)
        ax.set_ylim(-4 * sigma_zeta, 4 * sigma_zeta)
        ax.set_xlabel('ξ [km]')
        ax.set_ylabel('ζ [km]')
        ax.set_title('B-plane: Uncertainty Ellipse & Keyholes')
        ax.legend(loc='upper right', fontsize=9)
        ax.axhline(0, color='gray', linewidth=0.5, linestyle=':')
        ax.axvline(0, color='gray', linewidth=0.5, linestyle=':')
        ax.set_aspect('equal', 'datalim')

        # Right: keyhole probability ranking
        ax2 = axes[1]
        if keyholes:
            labels = [f"{kh.resonance_n}:{kh.resonance_m}\n({kh.encounter_year:.0f})"
                      for kh in keyholes[:10]]
            probs = [kh.impact_prob for kh in keyholes[:10]]
            colors = plt.cm.YlOrRd(np.linspace(0.4, 1.0, len(probs)))
            bars = ax2.barh(range(len(labels)), probs, color=colors)
            ax2.set_yticks(range(len(labels)))
            ax2.set_yticklabels(labels, fontsize=9)
            ax2.set_xscale('log')
            ax2.set_xlabel('Impact Probability')
            ax2.set_title('Top Resonant Keyholes by Impact Probability')
            ax2.invert_yaxis()

        plt.tight_layout()
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"[Plot] Saved: {save_path}")


def compute_torino_palermo_scales(impact_prob: float, diameter_km: float,
                                   encounter_vel_kms: float, dist_au: float) -> dict:
    """
    Compute Torino Scale and Palermo Technical Scale values.
    """
    # Kinetic energy [Mt TNT]
    rho = 2000  # kg/m^3
    R = diameter_km * 500  # m (radius)
    vol = (4 / 3) * np.pi * R**3
    mass_kg = rho * vol
    v_ms = encounter_vel_kms * 1e3
    KE_joule = 0.5 * mass_kg * v_ms**2
    KE_Mt = KE_joule / 4.184e15  # 1 Mt TNT = 4.184e15 J

    # Torino Scale (simplified table)
    if impact_prob < 1e-4:
        torino = 0
    elif KE_Mt < 1:
        torino = 1 if impact_prob < 0.01 else 2
    elif KE_Mt < 100:
        torino = 3 if impact_prob < 0.01 else 5
    elif KE_Mt < 1e5:
        torino = 6 if impact_prob < 0.1 else 7
    else:
        torino = 8 if impact_prob < 0.3 else (9 if impact_prob < 0.99 else 10)

    # Palermo Technical Scale: PS = log10(P / f_bg)
    # Background rate: Shoemaker-like, P_bg ~ 0.03 × (energy in Mt)^(-0.8) per year
    f_bg = 0.03 * KE_Mt**(-0.8) if KE_Mt > 0 else 1e-4
    t_remaining = 50.0  # years (assumed)
    palermo = np.log10(impact_prob / (f_bg * t_remaining)) if impact_prob > 0 else -10.0

    return {
        'kinetic_energy_Mt': KE_Mt,
        'torino_scale': torino,
        'palermo_scale': palermo,
        'f_background_per_yr': f_bg,
    }
