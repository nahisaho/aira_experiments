"""
Background models for dark matter direct detection experiments.
Includes neutrino floor (CEνNS), radiogenic, cosmogenic, and
surface backgrounds.
"""
import numpy as np
from typing import Dict, List, Tuple, Optional
from ..core.constants import (
    NuclearTarget, TARGETS, M_PROTON, N_AVOGADRO,
    helm_form_factor
)


class NeutrinoFloor:
    """Coherent elastic neutrino-nucleus scattering (CEνNS) background.

    Calculates the irreducible neutrino background (neutrino fog)
    from solar, atmospheric, and DSNB neutrinos.
    """

    # Neutrino source fluxes at Earth (cm^-2 s^-1) and energies
    NEUTRINO_SOURCES = {
        'pp':       {'flux': 5.98e10, 'E_max_mev': 0.423, 'type': 'solar'},
        'pep':      {'flux': 1.44e8,  'E_max_mev': 1.442, 'type': 'solar'},
        'hep':      {'flux': 7.98e3,  'E_max_mev': 18.77, 'type': 'solar'},
        '7Be_384':  {'flux': 4.56e8,  'E_max_mev': 0.384, 'type': 'solar_line'},
        '7Be_862':  {'flux': 4.56e9,  'E_max_mev': 0.862, 'type': 'solar_line'},
        '8B':       {'flux': 5.46e6,  'E_max_mev': 16.36, 'type': 'solar'},
        '13N':      {'flux': 2.17e8,  'E_max_mev': 1.199, 'type': 'solar'},
        '15O':      {'flux': 1.56e8,  'E_max_mev': 1.732, 'type': 'solar'},
        '17F':      {'flux': 3.40e6,  'E_max_mev': 1.740, 'type': 'solar'},
        'atm':      {'flux': 10.5,    'E_max_mev': 1e4,   'type': 'atmospheric'},
        'dsnb':     {'flux': 85.5,    'E_max_mev': 50.0,  'type': 'dsnb'},
    }

    def __init__(self, target: NuclearTarget):
        self.target = target
        self.sin2_theta_w = 0.2312

    def weak_charge(self) -> float:
        """Weak nuclear charge Q_W."""
        N = self.target.A - self.target.Z
        Z = self.target.Z
        return N - Z * (1 - 4 * self.sin2_theta_w)

    def cevns_cross_section(self, E_nu_mev: float, Er_kev: float) -> float:
        """CEνNS differential cross section dσ/dE_r (cm²/keV).

        dσ/dE_r = G_F² Q_W² M_N / (4π) × (1 - M_N E_r / (2 E_ν²)) × F²(q)
        """
        G_F = 1.1664e-5  # GeV^-2
        hbarc = 197.3269804  # MeV·fm
        Q_w = self.weak_charge()
        m_N = self.target.mass_gev
        A = self.target.A

        E_nu_gev = E_nu_mev * 1e-3
        Er_gev = Er_kev * 1e-6

        # Kinematic check
        Er_max = 2 * E_nu_gev**2 / (m_N + 2 * E_nu_gev)
        if Er_gev > Er_max:
            return 0.0

        # Form factor: q² = 2 m_N E_R, with m_N in GeV and E_R in keV → q in MeV
        q_mev = np.sqrt(2 * m_N * Er_kev)  # GeV·keV = MeV²
        F2 = helm_form_factor(np.array([q_mev]), A)[0]

        # Cross section in natural units (GeV^-2)
        dsigma = G_F**2 * Q_w**2 * m_N / (4 * np.pi) * \
                 (1 - m_N * Er_gev / (2 * E_nu_gev**2)) * F2

        # Convert to cm²/keV
        gev2_to_cm2 = 0.389379e-27  # GeV^-2 to cm²  (×10^{-27})
        dsigma_cm2 = dsigma * gev2_to_cm2 * 1e-6  # per keV

        return dsigma_cm2

    def neutrino_spectrum(self, source: str, E_nu_mev: np.ndarray) -> np.ndarray:
        """Neutrino energy spectrum for given source.

        Returns dΦ/dE (cm^-2 s^-1 MeV^-1).
        """
        info = self.NEUTRINO_SOURCES[source]
        E_max = info['E_max_mev']
        flux = info['flux']

        if info['type'] == 'solar_line':
            # Monoenergetic line (delta function approximation)
            sigma = 0.01 * E_max  # small width
            return flux * np.exp(-0.5 * ((E_nu_mev - E_max) / sigma)**2) / \
                   (sigma * np.sqrt(2 * np.pi))

        elif info['type'] == 'solar':
            # Approximate continuous spectrum (beta-like)
            mask = E_nu_mev < E_max
            spec = np.zeros_like(E_nu_mev)
            spec[mask] = flux * E_nu_mev[mask]**2 * \
                        (1 - E_nu_mev[mask] / E_max)**2 / (E_max**3 / 6)
            return spec

        elif info['type'] == 'atmospheric':
            # Power law with exponential cutoff
            E0 = 500.0  # MeV reference
            return flux * (E_nu_mev / E0)**(-2.7) * np.exp(-E_nu_mev / E_max)

        elif info['type'] == 'dsnb':
            # Fermi-Dirac with T ~ 5 MeV
            T = 5.0  # MeV
            return flux * E_nu_mev**2 / (np.exp(E_nu_mev / T) + 1) / (2 * T**3 * 0.9)

        return np.zeros_like(E_nu_mev)

    def recoil_rate_source(self, source: str, Er_kev: np.ndarray,
                           n_Enu: int = 500) -> np.ndarray:
        """Recoil spectrum from a single neutrino source (events/keV/kg/day)."""
        info = self.NEUTRINO_SOURCES[source]
        E_max_nu = info['E_max_mev']

        n_atoms = N_AVOGADRO / self.target.A * 1e3  # atoms/kg

        rates = np.zeros_like(Er_kev)
        E_nu = np.linspace(0.01, E_max_nu, n_Enu)
        dEnu = E_nu[1] - E_nu[0]
        spectrum = self.neutrino_spectrum(source, E_nu)

        for i, Er in enumerate(Er_kev):
            integrand = np.zeros(n_Enu)
            for j, Enu in enumerate(E_nu):
                integrand[j] = spectrum[j] * self.cevns_cross_section(Enu, Er)

            rates[i] = np.trapz(integrand, E_nu) * n_atoms * 86400

        return rates

    def total_neutrino_rate(self, Er_kev: np.ndarray,
                            sources: Optional[List[str]] = None) -> np.ndarray:
        """Total CEνNS recoil rate from all neutrino sources."""
        if sources is None:
            sources = list(self.NEUTRINO_SOURCES.keys())

        total = np.zeros_like(Er_kev)
        for src in sources:
            total += self.recoil_rate_source(src, Er_kev)

        return total

    def neutrino_floor_cross_section(self, m_dm_gev: float,
                                      exposure_kg_day: float,
                                      cl: float = 0.90) -> float:
        """Estimate neutrino floor cross section for given DM mass.

        Returns σ_SI (cm²) at which neutrino background limits discovery.
        Uses the criterion S/√B = z_cl for the cross section where
        signal equals statistical fluctuation of neutrino background.
        """
        from scipy.stats import norm

        Er = np.linspace(max(0.5, 1.0), 50.0, 100)

        # Compute neutrino background events
        nu_rate = self.total_neutrino_rate(Er, ['8B', 'hep', 'atm'])
        n_nu = np.trapz(nu_rate, Er) * exposure_kg_day

        z = norm.ppf(cl)
        if n_nu < 0.01:
            return 1e-50

        # Compute signal rate at reference σ = 1e-45 cm²
        from ..signals.dm_signals import WIMPSignal
        ref_sigma = 1e-45
        sig = WIMPSignal(m_dm_gev, ref_sigma, self.target, 1.0)
        dR = sig.differential_rate(Er)
        signal_per_exposure = np.trapz(dR, Er)

        if signal_per_exposure < 1e-30:
            return 1e-40

        # σ_floor: S = z × √B
        # S = signal_per_exposure × exposure × (σ/ref_σ)
        # z√B = z × √n_nu
        sigma_floor = z * np.sqrt(n_nu) / (signal_per_exposure * exposure_kg_day) * ref_sigma

        return max(sigma_floor, 1e-50)


class RadiogenicBackground:
    """Radiogenic background model from detector materials."""

    ISOTOPES = {
        'U238': {'activity_bq_kg': 1e-3, 'E_range_kev': (0.1, 3000)},
        'Th232': {'activity_bq_kg': 5e-4, 'E_range_kev': (0.1, 2800)},
        'K40': {'activity_bq_kg': 1e-2, 'E_range_kev': (0.1, 1460)},
        'Rn222': {'activity_bq_kg': 1e-4, 'E_range_kev': (0.1, 5500)},
        'Kr85': {'activity_bq_kg': 1e-5, 'E_range_kev': (0.1, 687)},
        'Pb210': {'activity_bq_kg': 5e-4, 'E_range_kev': (0.1, 46.5)},
        'Ar39': {'activity_bq_kg': 1.0, 'E_range_kev': (0.1, 565)},
    }

    def __init__(self, target_name: str, shielding_factor: float = 1e-6):
        self.target_name = target_name
        self.shielding = shielding_factor

    def flat_background_rate(self, Er_kev: np.ndarray,
                             isotope: str = 'all') -> np.ndarray:
        """Approximate flat background rate (events/keV/kg/day).

        Simplified model with exponential + flat components.
        """
        rate = np.zeros_like(Er_kev)

        if isotope == 'all':
            isotopes = self.ISOTOPES.keys()
        else:
            isotopes = [isotope]

        for iso in isotopes:
            info = self.ISOTOPES[iso]
            activity = info['activity_bq_kg'] * self.shielding
            E_min, E_max = info['E_range_kev']

            mask = (Er_kev >= E_min) & (Er_kev <= E_max)
            # Exponential + flat spectrum
            rate[mask] += activity * 86400 / (E_max - E_min) * \
                         (0.7 * np.exp(-Er_kev[mask] / (0.1 * E_max)) + 0.3)

        return rate


class CosmogenicBackground:
    """Cosmogenic activation backgrounds."""

    ACTIVATION_PRODUCTS = {
        'Xe': {'H3': 2e-4, 'Xe127': 1e-3, 'I125': 5e-4},
        'Ar': {'Ar39': 1.0, 'Ar37': 0.05},
        'Ge': {'H3': 1e-3, 'Ge68': 2e-3, 'Co57': 1e-4},
        'Na': {'H3': 5e-4, 'Na22': 1e-3},
    }

    def __init__(self, target_element: str, cooling_days: float = 180.0):
        self.element = target_element
        self.cooling = cooling_days

    def activation_rate(self, Er_kev: np.ndarray) -> np.ndarray:
        """Cosmogenic activation background rate (events/keV/kg/day)."""
        products = self.ACTIVATION_PRODUCTS.get(self.element, {})
        rate = np.zeros_like(Er_kev)

        for product, activity_bq_kg in products.items():
            # Decay during cooling
            half_lives = {'H3': 4500, 'Xe127': 36.4, 'I125': 59.4,
                         'Ar39': 9.8e9, 'Ar37': 35.0, 'Ge68': 271,
                         'Co57': 272, 'Na22': 950}
            t_half = half_lives.get(product, 365)
            decay_factor = np.exp(-0.693 * self.cooling / t_half)

            rate += activity_bq_kg * decay_factor * 86400 / 1000  # Flat approx

        return rate


class BackgroundBudget:
    """Complete background budget combining all sources."""

    def __init__(self, target: NuclearTarget,
                 shielding_factor: float = 1e-6,
                 cooling_days: float = 180.0):
        self.target = target
        self.nu_floor = NeutrinoFloor(target)
        self.radiogenic = RadiogenicBackground(target.name[:2], shielding_factor)
        self.cosmogenic = CosmogenicBackground(target.name[:2], cooling_days)

    def total_background(self, Er_kev: np.ndarray,
                         include_neutrinos: bool = True) -> Dict[str, np.ndarray]:
        """Calculate all background components."""
        result = {}

        result['radiogenic'] = self.radiogenic.flat_background_rate(Er_kev)
        result['cosmogenic'] = self.cosmogenic.activation_rate(Er_kev)

        if include_neutrinos:
            result['neutrino_cevns'] = self.nu_floor.total_neutrino_rate(
                Er_kev, ['pp', '7Be_862', '8B']
            )

        result['total'] = sum(result.values())
        return result

    def evaluate_reduction_strategies(self, Er_kev: np.ndarray) -> Dict:
        """Evaluate background reduction strategies systematically."""
        strategies = {
            'baseline': {'shielding': 1e-6, 'cooling': 180, 'fiducial': 1.0},
            'enhanced_shielding': {'shielding': 1e-8, 'cooling': 180, 'fiducial': 1.0},
            'underground_cooling': {'shielding': 1e-6, 'cooling': 730, 'fiducial': 1.0},
            'fiducial_cut': {'shielding': 1e-6, 'cooling': 180, 'fiducial': 0.5},
            'optimal': {'shielding': 1e-8, 'cooling': 730, 'fiducial': 0.5},
        }

        results = {}
        for name, params in strategies.items():
            radio = RadiogenicBackground(self.target.name[:2], params['shielding'])
            cosmo = CosmogenicBackground(self.target.name[:2], params['cooling'])

            bg = radio.flat_background_rate(Er_kev) + cosmo.activation_rate(Er_kev)
            bg *= params['fiducial']

            results[name] = {
                'params': params,
                'rate': bg,
                'total_rate': float(np.trapz(bg, Er_kev)),
                'reduction_factor': None
            }

        baseline_total = results['baseline']['total_rate']
        for name in results:
            if baseline_total > 0:
                results[name]['reduction_factor'] = \
                    results[name]['total_rate'] / baseline_total

        return results
