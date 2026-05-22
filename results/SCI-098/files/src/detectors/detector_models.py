"""
Detector models for dark matter experiments.
Implements directional detectors (CYGNUS/MIMAC), liquid noble gas,
solid-state, and crystal detectors.
"""
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple, List
from ..core.constants import NuclearTarget, TARGETS, V_0, V_EARTH, V_ESC


@dataclass
class DetectorConfig:
    """Generic detector configuration."""
    name: str
    target_name: str
    mass_kg: float
    threshold_kev: float
    max_energy_kev: float
    energy_resolution: float      # fractional σ/E at 1 keV
    efficiency: float              # detection efficiency (0-1)
    exposure_days: float
    background_rate: float         # events/keV/kg/day (flat background)

    @property
    def exposure_kg_day(self) -> float:
        return self.mass_kg * self.exposure_days

    @property
    def target(self) -> NuclearTarget:
        return TARGETS[self.target_name]

    def energy_resolution_sigma(self, E_kev: float) -> float:
        """Energy resolution σ(E) in keV."""
        return self.energy_resolution * np.sqrt(E_kev)


class DirectionalDetector:
    """CYGNUS/MIMAC-type directional detector with angular sensitivity.

    Models gas TPC with 3D track reconstruction capability.
    """

    def __init__(self, config: DetectorConfig,
                 angular_resolution_deg: float = 30.0,
                 head_tail_recognition: bool = True,
                 gas_pressure_torr: float = 50.0,
                 track_length_threshold_mm: float = 1.0):
        self.config = config
        self.angular_res = np.radians(angular_resolution_deg)
        self.head_tail = head_tail_recognition
        self.pressure = gas_pressure_torr
        self.track_threshold = track_length_threshold_mm

    def recoil_direction_distribution(self, cos_theta: np.ndarray,
                                       m_dm_gev: float) -> np.ndarray:
        """Angular distribution of nuclear recoils from DM wind.

        dR/d(cosθ) ∝ exp(-(v_min² + v_e² sin²θ)/(2σ²))

        Args:
            cos_theta: cosine of angle from DM wind direction
            m_dm_gev: dark matter mass

        Returns:
            relative angular distribution
        """
        sin2_theta = 1 - cos_theta**2
        v_e = V_EARTH  # km/s

        # Anisotropy from Earth motion through DM halo
        sigma_v = V_0 / np.sqrt(2)
        exponent = -v_e**2 * sin2_theta / (2 * sigma_v**2)

        dist = np.exp(exponent)

        # Smear with angular resolution
        if self.angular_res > 0:
            from scipy.ndimage import gaussian_filter1d
            sigma_bins = self.angular_res / (np.pi / len(cos_theta))
            dist = gaussian_filter1d(dist, sigma=max(1, sigma_bins))

        # Head-tail asymmetry
        if self.head_tail:
            dist *= (1 + 0.5 * cos_theta)  # forward-backward asymmetry

        return dist / np.max(dist)

    def track_length_mm(self, Er_kev: float, A: int) -> float:
        """Estimated recoil track length in gas TPC (mm).

        Simplified Lindhard/SRIM scaling.
        """
        # Approximate range in CF4 gas at given pressure
        # R ∝ E^0.7 / (ρ × A^0.5)
        rho_rel = self.pressure / 760.0  # relative to STP
        return 0.5 * (Er_kev / 50.0)**0.7 / (rho_rel * np.sqrt(A / 19.0))

    def angular_sensitivity(self, Er_kev: float) -> float:
        """Angular resolution achievable at given recoil energy."""
        track_len = self.track_length_mm(Er_kev, self.config.target.A)
        if track_len < self.track_threshold:
            return np.pi  # no directionality below threshold

        # Resolution improves with track length
        sigma = self.angular_res * (self.track_threshold / track_len)**0.5
        return min(sigma, np.pi)

    def discovery_reach_directional(self, m_dm_gev: float,
                                      n_signal: int = 30) -> Dict:
        """Calculate directional discovery reach.

        With head-tail: uses dipole (Rayleigh) test.
        Without head-tail: uses quadrupole test (cos²θ).
        """
        cos_theta = np.linspace(-1, 1, 200)
        angular_dist = self.recoil_direction_distribution(cos_theta, m_dm_gev)
        angular_dist /= np.trapz(angular_dist, cos_theta)

        if self.head_tail:
            # Dipole: mean cos(θ)
            mean_cos = np.trapz(cos_theta * angular_dist, cos_theta)
            stat = mean_cos**2
        else:
            # Quadrupole: <cos²θ> - 1/3
            mean_cos2 = np.trapz(cos_theta**2 * angular_dist, cos_theta)
            stat = (mean_cos2 - 1./3.)**2
            mean_cos = np.sqrt(stat)

        n_3sigma = int(np.ceil(9.0 / max(stat, 1e-6)))
        n_5sigma = int(np.ceil(25.0 / max(stat, 1e-6)))

        return {
            'mean_cosine': float(mean_cos),
            'n_events_3sigma': min(n_3sigma, int(1e8)),
            'n_events_5sigma': min(n_5sigma, int(1e8)),
            'angular_resolution_deg': np.degrees(self.angular_res),
            'head_tail': self.head_tail
        }


class LiquidNobleDetector:
    """Liquid noble gas detector (LXe/LAr dual-phase TPC)."""

    def __init__(self, config: DetectorConfig,
                 s1_threshold_pe: float = 3.0,
                 s2_threshold_pe: float = 100.0,
                 light_yield_pe_per_kev: float = 8.0,
                 electron_lifetime_us: float = 1000.0):
        self.config = config
        self.s1_threshold = s1_threshold_pe
        self.s2_threshold = s2_threshold_pe
        self.light_yield = light_yield_pe_per_kev
        self.e_lifetime = electron_lifetime_us

    def nuclear_recoil_efficiency(self, Er_kev: np.ndarray) -> np.ndarray:
        """Detection efficiency for nuclear recoils.

        Includes Lindhard quenching and S1/S2 thresholds.
        """
        A = self.config.target.A
        Z = self.config.target.Z

        # Lindhard quenching factor
        epsilon = 11.5 * Er_kev * Z**(-7./3.)
        k = 0.133 * Z**(2./3.) * A**(-1./2.)
        g = 3 * epsilon**0.15 + 0.7 * epsilon**0.6 + epsilon
        L = k * g / (1 + k * g)

        # Effective energy (keVee)
        E_ee = Er_kev * L

        # S1 detection efficiency (Poisson threshold)
        mean_pe = E_ee * self.light_yield
        from scipy.stats import poisson
        eff_s1 = 1 - poisson.cdf(self.s1_threshold - 1, np.clip(mean_pe, 1e-10, None))

        # Combined with energy threshold
        eff_threshold = np.where(Er_kev > self.config.threshold_kev, 1.0, 0.0)

        return eff_s1 * eff_threshold * self.config.efficiency

    def discrimination_power(self, Er_kev: float) -> float:
        """ER/NR discrimination power at given energy.

        Returns log10 of ER rejection factor.
        """
        if Er_kev < self.config.threshold_kev:
            return 0.0
        # Typical dual-phase TPC discrimination
        # Better at higher energies
        base_rejection = 3.0  # log10, i.e., 99.9% ER rejection
        energy_factor = min(1.0, (Er_kev / 5.0)**0.5)
        return base_rejection * energy_factor


class SolidStateDetector:
    """Solid-state detector (Ge/Si bolometer, CCD)."""

    def __init__(self, config: DetectorConfig,
                 phonon_resolution_ev: float = 10.0,
                 ionization_yield: float = 0.3):
        self.config = config
        self.phonon_res = phonon_resolution_ev * 1e-3  # keV
        self.ion_yield = ionization_yield

    def nuclear_recoil_efficiency(self, Er_kev: np.ndarray) -> np.ndarray:
        """Detection efficiency including ionization quenching."""
        # Lindhard quenching for Ge/Si
        A = self.config.target.A
        Z = self.config.target.Z
        epsilon = 11.5 * Er_kev * Z**(-7./3.)
        k = 0.133 * Z**(2./3.) * A**(-1./2.)
        g = 3 * epsilon**0.15 + 0.7 * epsilon**0.6 + epsilon
        L = k * g / (1 + k * g)

        E_ion = Er_kev * L * self.ion_yield
        eff = np.where(E_ion > 3 * self.phonon_res, 1.0,
                      np.where(E_ion > self.phonon_res,
                              (E_ion - self.phonon_res) / (2 * self.phonon_res),
                              0.0))
        return eff * self.config.efficiency


class CrystalDetector:
    """NaI(Tl) scintillation detector (DAMA/LIBRA type)."""

    def __init__(self, config: DetectorConfig,
                 quenching_na: float = 0.3,
                 quenching_i: float = 0.09):
        self.config = config
        self.q_na = quenching_na
        self.q_i = quenching_i

    def quenched_energy(self, Er_kev: np.ndarray,
                        recoil_species: str = 'Na') -> np.ndarray:
        """Convert nuclear recoil to electron-equivalent energy."""
        q = self.q_na if recoil_species == 'Na' else self.q_i
        return Er_kev * q


# === Pre-configured detector specifications ===

def get_xenon_nt() -> DetectorConfig:
    """XENON-nT like detector."""
    return DetectorConfig(
        name='XENON-nT', target_name='Xe131', mass_kg=5900,
        threshold_kev=1.0, max_energy_kev=70.0,
        energy_resolution=0.4, efficiency=0.85,
        exposure_days=365, background_rate=1.5e-5
    )

def get_darwin() -> DetectorConfig:
    """DARWIN/XLZD next-generation xenon detector."""
    return DetectorConfig(
        name='DARWIN', target_name='Xe131', mass_kg=40000,
        threshold_kev=0.5, max_energy_kev=70.0,
        energy_resolution=0.35, efficiency=0.90,
        exposure_days=3650, background_rate=5e-6
    )

def get_darkside20k() -> DetectorConfig:
    """DarkSide-20k liquid argon detector."""
    return DetectorConfig(
        name='DarkSide-20k', target_name='Ar40', mass_kg=20000,
        threshold_kev=0.6, max_energy_kev=200.0,
        energy_resolution=0.15, efficiency=0.80,
        exposure_days=3650, background_rate=1e-6
    )

def get_supercdms() -> DetectorConfig:
    """SuperCDMS germanium detector."""
    return DetectorConfig(
        name='SuperCDMS', target_name='Ge76', mass_kg=30,
        threshold_kev=0.04, max_energy_kev=50.0,
        energy_resolution=0.05, efficiency=0.70,
        exposure_days=1825, background_rate=1e-4
    )

def get_cygnus() -> DetectorConfig:
    """CYGNUS directional detector (CF4 gas)."""
    return DetectorConfig(
        name='CYGNUS', target_name='CF4', mass_kg=1000,
        threshold_kev=3.0, max_energy_kev=200.0,
        energy_resolution=0.3, efficiency=0.50,
        exposure_days=3650, background_rate=1e-3
    )

def get_cosine100() -> DetectorConfig:
    """COSINE-100 NaI detector."""
    return DetectorConfig(
        name='COSINE-100', target_name='Na23', mass_kg=106,
        threshold_kev=1.0, max_energy_kev=20.0,
        energy_resolution=0.5, efficiency=0.65,
        exposure_days=1095, background_rate=2.7
    )
