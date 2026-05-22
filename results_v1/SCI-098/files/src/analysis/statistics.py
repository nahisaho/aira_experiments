"""
Statistical analysis tools for dark matter detection.
Includes sensitivity calculations, annual modulation analysis,
multi-target complementarity, and discovery reach estimation.
"""
import numpy as np
from typing import Dict, List, Tuple, Optional
from ..core.constants import (
    NuclearTarget, TARGETS, RHO_DM_LOCAL, V_0, V_ESC, V_EARTH, V_SUN,
    helm_form_factor, v_min, eta_integral, M_PROTON, N_AVOGADRO,
    C_LIGHT
)


from ..signals.dm_signals import WIMPSignal
from ..backgrounds.background_models import NeutrinoFloor, BackgroundBudget


class SensitivityCalculator:
    """Calculate exclusion limits and discovery reach."""

    def __init__(self, target: NuclearTarget, exposure_kg_day: float,
                 threshold_kev: float = 1.0, max_energy_kev: float = 50.0,
                 background_rate: float = 0.0, efficiency: float = 1.0):
        self.target = target
        self.exposure = exposure_kg_day
        self.threshold = threshold_kev
        self.max_energy = max_energy_kev
        self.bg_rate = background_rate
        self.efficiency = efficiency

    def exclusion_limit_90cl(self, m_dm_values: np.ndarray) -> np.ndarray:
        """90% CL exclusion limit on σ_SI using Poisson statistics.

        Optimum interval method approximation.
        """
        sigma_limits = np.zeros_like(m_dm_values)
        Er = np.linspace(self.threshold, self.max_energy, 200)

        for i, m_dm in enumerate(m_dm_values):
            # Reference signal at σ = 1e-45 cm²
            sig = WIMPSignal(m_dm, 1e-45, self.target, 1.0)
            dR = sig.differential_rate(Er) * self.efficiency

            total_signal_per_sigma = np.trapz(dR, Er) * self.exposure
            total_bg = self.bg_rate * (self.max_energy - self.threshold) * self.exposure

            if total_signal_per_sigma < 1e-30:
                sigma_limits[i] = np.inf
                continue

            # 90% CL Poisson upper limit
            if total_bg < 0.5:
                n_ul = 2.3  # Feldman-Cousins, 0 observed
            else:
                n_ul = 1.64 * np.sqrt(total_bg) + 1.0  # Gaussian approx

            sigma_limits[i] = n_ul / total_signal_per_sigma * 1e-45

        return sigma_limits

    def discovery_reach_3sigma(self, m_dm_values: np.ndarray) -> np.ndarray:
        """3σ discovery reach on σ_SI."""
        sigma_disc = np.zeros_like(m_dm_values)
        Er = np.linspace(self.threshold, self.max_energy, 200)

        for i, m_dm in enumerate(m_dm_values):
            sig = WIMPSignal(m_dm, 1e-45, self.target, 1.0)
            dR = sig.differential_rate(Er) * self.efficiency
            total_signal_per_sigma = np.trapz(dR, Er) * self.exposure
            total_bg = self.bg_rate * (self.max_energy - self.threshold) * self.exposure

            if total_signal_per_sigma < 1e-30:
                sigma_disc[i] = np.inf
                continue

            n_disc = 3 * np.sqrt(max(total_bg, 1.0)) + 3
            sigma_disc[i] = n_disc / total_signal_per_sigma * 1e-45

        return sigma_disc


class AnnualModulation:
    """Annual modulation signal analysis.

    The Earth's velocity relative to the DM halo varies annually,
    creating a ~2-6% modulation in the event rate.
    """

    def __init__(self, target: NuclearTarget, m_dm_gev: float,
                 sigma_si_cm2: float, exposure_kg_day: float):
        self.target = target
        self.m_dm = m_dm_gev
        self.sigma_si = sigma_si_cm2
        self.exposure = exposure_kg_day

    def earth_velocity(self, t_days: float) -> float:
        """Earth velocity in galactic frame as function of time.

        v_E(t) ≈ v_sun + v_orb × cos(2π(t - t_peak)/T)

        Args:
            t_days: days from Jan 1

        Returns:
            velocity in km/s
        """
        v_sun = V_SUN
        v_orb = 29.8  # km/s (Earth orbital speed)
        t_peak = 152.0  # June 2 (day of year, peak)
        omega = 2 * np.pi / 365.25

        return v_sun + v_orb * np.cos(omega * (t_days - t_peak))

    def modulated_rate(self, Er_kev: np.ndarray,
                       t_days: float) -> np.ndarray:
        """Event rate at specific time of year."""
        v_e = self.earth_velocity(t_days)

        m_N = self.target.mass_gev
        mu_p = (self.m_dm * M_PROTON * 1e-3) / (self.m_dm + M_PROTON * 1e-3)
        mu_N = self.target.reduced_mass_dm(self.m_dm)
        A = self.target.A
        sigma_N = self.sigma_si * (mu_N / mu_p)**2 * A**2

        N_T = N_AVOGADRO * 1e3 / A
        n_chi = RHO_DM_LOCAL / self.m_dm
        c_cm = C_LIGHT * 100

        rates = np.zeros_like(Er_kev)
        for i, Er in enumerate(Er_kev):
            q_mev = np.sqrt(2 * m_N * Er)  # GeV·keV = MeV²
            F2 = helm_form_factor(np.array([q_mev]), A)[0]
            vm = v_min(Er, self.m_dm, m_N)
            eta = eta_integral(vm, v_e=v_e)
            eta_cgs = eta * 1e-5

            dR_atom = n_chi * sigma_N * m_N * c_cm**2 * F2 * eta_cgs / (2 * mu_N**2)
            rates[i] = dR_atom * 1e-6 * N_T * 86400

        return rates

    def modulation_amplitude(self, Er_kev: np.ndarray) -> np.ndarray:
        """Modulation amplitude S_m(E) = (R_max - R_min) / 2."""
        R_june = self.modulated_rate(Er_kev, 152.0)  # June peak
        R_dec = self.modulated_rate(Er_kev, 335.0)   # December minimum
        return (R_june - R_dec) / 2

    def modulation_fraction(self, Er_range: Tuple[float, float] = (2, 6)) -> float:
        """Fractional modulation amplitude integrated over energy range."""
        Er = np.linspace(Er_range[0], Er_range[1], 100)
        R_avg = self.modulated_rate(Er, 0)  # average rate
        S_m = self.modulation_amplitude(Er)

        avg_rate = np.trapz(R_avg, Er)
        mod_amp = np.trapz(S_m, Er)

        if avg_rate < 1e-30:
            return 0.0
        return mod_amp / avg_rate

    def detection_significance(self, n_years: float = 3.0,
                                Er_range: Tuple[float, float] = (2, 6),
                                bg_rate: float = 0.0) -> Dict:
        """Statistical significance of annual modulation detection.

        Uses likelihood ratio test for modulation vs. constant rate.
        Scales exposure with observation period.
        """
        Er = np.linspace(Er_range[0], Er_range[1], 50)

        # Scale exposure with observation period
        daily_exposure_kg = self.exposure / 365.25  # kg (assuming exposure given as kg×1year)
        total_exposure = daily_exposure_kg * n_years * 365.25

        R_june = self.modulated_rate(Er, 152.0)
        R_dec = self.modulated_rate(Er, 335.0)
        R_avg = 0.5 * (R_june + R_dec)

        S_m = np.trapz(self.modulation_amplitude(Er), Er)
        R_mean = np.trapz(R_avg, Er)

        total_signal = R_mean * total_exposure
        total_bg = bg_rate * (Er_range[1] - Er_range[0]) * total_exposure
        total_events = total_signal + total_bg

        if total_events > 0 and R_mean > 0:
            A_m = S_m / R_mean
            significance = A_m * np.sqrt(total_events / 2)
        else:
            significance = 0

        return {
            'total_events': total_events,
            'modulation_fraction': S_m / (R_mean + 1e-30),
            'significance_sigma': significance,
            'n_years': n_years,
            'exposure_kg_day': total_exposure,
            'energy_range_kev': Er_range,
        }


class MultiTargetComplementarity:
    """Analyze complementarity of multiple detector targets."""

    def __init__(self, targets: Dict[str, Dict]):
        """
        Args:
            targets: dict of {name: {'target': NuclearTarget,
                                      'exposure': float,
                                      'threshold': float,
                                      'max_energy': float,
                                      'background': float,
                                      'efficiency': float}}
        """
        self.targets = targets

    def combined_sensitivity(self, m_dm_values: np.ndarray) -> Dict:
        """Calculate individual and combined sensitivities."""
        results = {}

        for name, config in self.targets.items():
            calc = SensitivityCalculator(
                config['target'], config['exposure'],
                config['threshold'], config['max_energy'],
                config['background'], config['efficiency']
            )
            results[name] = calc.exclusion_limit_90cl(m_dm_values)

        # Combined: take best limit at each mass
        combined = np.full_like(m_dm_values, 1e-30)
        for name in results:
            combined = np.minimum(combined, results[name])

        results['combined'] = combined
        return results

    def target_response_matrix(self, m_dm: float,
                                sigma_si: float = 1e-46) -> Dict:
        """Response comparison across targets at fixed DM parameters."""
        response = {}
        for name, config in self.targets.items():
            Er = np.linspace(config['threshold'], config['max_energy'], 200)
            sig = WIMPSignal(m_dm, sigma_si, config['target'], config['exposure'])
            dR = sig.differential_rate(Er)
            total = np.trapz(dR, Er) * config['exposure'] * config['efficiency']
            E_peak = Er[np.argmax(dR)] if np.max(dR) > 0 else 0

            response[name] = {
                'total_events': total,
                'peak_energy_kev': E_peak,
                'max_rate': float(np.max(dR)),
                'A': config['target'].A,
                'Z': config['target'].Z,
            }
        return response

    def si_sd_discrimination(self, m_dm_values: np.ndarray) -> Dict:
        """Assess SI vs SD discrimination capability across targets.

        Targets with different Z/N ratios help distinguish
        spin-independent from spin-dependent interactions.
        """
        results = {}
        for name, config in self.targets.items():
            t = config['target']
            N = t.A - t.Z
            results[name] = {
                'Z': t.Z, 'N': N, 'A': t.A,
                'Z_over_A': t.Z / t.A,
                'N_over_A': N / t.A,
                'spin': t.spin,
                'Sp': t.sp, 'Sn': t.sn,
                'SI_enhancement': t.A**2,
                'SD_sensitivity': 'yes' if t.spin > 0 else 'no',
            }

        return results


class NeutrinoFloorCalculator:
    """Calculate neutrino floor for different targets and exposures."""

    def __init__(self, target: NuclearTarget):
        self.target = target
        self.nu_floor = NeutrinoFloor(target)

    def compute_floor(self, m_dm_values: np.ndarray,
                      exposure_kg_day: float) -> np.ndarray:
        """Compute neutrino floor cross section vs DM mass."""
        floor = np.zeros_like(m_dm_values)

        for i, m_dm in enumerate(m_dm_values):
            floor[i] = self.nu_floor.neutrino_floor_cross_section(
                m_dm, exposure_kg_day
            )

        return floor

    def exposure_to_floor(self, m_dm_gev: float,
                           target_exposures: np.ndarray) -> np.ndarray:
        """Sensitivity vs exposure, showing approach to neutrino floor."""
        sensitivities = np.zeros_like(target_exposures)

        for i, exp in enumerate(target_exposures):
            calc = SensitivityCalculator(
                self.target, exp,
                threshold_kev=1.0, max_energy_kev=50.0,
                background_rate=1e-5, efficiency=0.85
            )
            sensitivities[i] = calc.exclusion_limit_90cl(
                np.array([m_dm_gev])
            )[0]

        return sensitivities
