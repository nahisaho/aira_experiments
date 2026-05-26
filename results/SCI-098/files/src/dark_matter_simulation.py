#!/usr/bin/env python3
"""
Next-Generation Dark Matter Direct Detection Simulation Framework
=================================================================
Monte Carlo simulation framework for evaluating detection strategies
for various dark matter candidates, directional sensitivity,
neutrino floor projections, background reduction, multi-target
complementarity, and annual modulation statistical power.

Designed as a lightweight Python implementation inspired by
GEANT4/ROOT analysis workflows.
"""

import numpy as np
from scipy import integrate, special, stats, interpolate
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import json

# Physical constants
C_LIGHT = 3.0e10       # cm/s
HBAR_C = 197.3e-13     # MeV·cm
GF = 1.166e-11         # MeV^-2 (Fermi constant)
M_PROTON = 938.272      # MeV/c^2
M_NEUTRON = 939.565     # MeV/c^2
M_ELECTRON = 0.511      # MeV/c^2
RHO_DM = 0.3            # GeV/cm^3 (local DM density)
V0 = 220.0              # km/s (local circular velocity)
V_ESC = 544.0            # km/s (escape velocity)
V_EARTH = 232.0          # km/s (Earth velocity through halo)
YEAR_SECONDS = 3.156e7   # seconds in a year

# Output directory
FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)


# ============================================================
# Module 1: Dark Matter Velocity Distribution & Kinematics
# ============================================================
class DMHalo:
    """Standard Halo Model for dark matter velocity distribution."""

    def __init__(self, v0=V0, v_esc=V_ESC, v_earth=V_EARTH, rho=RHO_DM):
        self.v0 = v0
        self.v_esc = v_esc
        self.v_earth = v_earth
        self.rho = rho
        self.k0 = self._normalization()

    def _normalization(self):
        x_esc = self.v_esc / self.v0
        return (np.pi * self.v0**2)**(3/2) * (
            special.erf(x_esc) - 2*x_esc/np.sqrt(np.pi) * np.exp(-x_esc**2)
        )

    def velocity_distribution(self, v):
        """Maxwell-Boltzmann truncated at escape velocity."""
        if v > self.v_esc + self.v_earth:
            return 0.0
        return (4*np.pi*v**2 / self.k0) * np.exp(-v**2 / self.v0**2)

    def mean_inverse_speed(self, v_min):
        """Compute <1/v> integral for rate calculation (km/s)^-1."""
        x = v_min / self.v0
        y = self.v_earth / self.v0
        z = self.v_esc / self.v0

        if x > y + z:
            return 0.0

        k1 = special.erf(z) - 2*z/np.sqrt(np.pi)*np.exp(-z**2)

        if x < abs(y - z):
            eta = special.erf(x + y) - special.erf(x - y) - \
                  4*y/(np.sqrt(np.pi)) * np.exp(-z**2)
            return eta / (2 * self.v0 * y * k1)
        else:
            eta = special.erf(z) - special.erf(x - y) - \
                  2*(y + z - x)/(np.sqrt(np.pi)) * np.exp(-z**2)
            return eta / (2 * self.v0 * y * k1)


# ============================================================
# Module 2: WIMP Differential Rate
# ============================================================
class WIMPDetector:
    """WIMP-nucleon scattering rate calculator."""

    # Target properties: (Z, A, atomic_mass_GeV)
    TARGETS = {
        'Xe': (54, 131, 121.76),
        'Ar': (18, 40, 37.21),
        'Ge': (32, 73, 67.93),
        'NaI': (11, 23, 21.41),  # Na component
        'F': (9, 19, 17.69),     # for directional (CF4)
        'He': (2, 4, 3.73),      # for CYGNUS
    }

    def __init__(self, target='Xe', exposure_kg_yr=1000.0):
        self.target = target
        Z, A, mA = self.TARGETS[target]
        self.Z = Z
        self.A = A
        self.mA_GeV = A * 0.9315  # GeV
        self.exposure = exposure_kg_yr
        self.halo = DMHalo()

    def reduced_mass(self, m_chi):
        """Reduced mass of DM-nucleus system in GeV."""
        return (m_chi * self.mA_GeV) / (m_chi + self.mA_GeV)

    def reduced_mass_nucleon(self, m_chi):
        """Reduced mass of DM-nucleon system in GeV."""
        m_n = 0.9396
        return (m_chi * m_n) / (m_chi + m_n)

    def v_min(self, E_r, m_chi):
        """Minimum velocity for recoil energy E_r (keV), in km/s."""
        E_r_GeV = E_r * 1e-6
        mu = self.reduced_mass(m_chi)
        return np.sqrt(self.mA_GeV * E_r_GeV / (2 * mu**2)) * 3e5  # km/s

    def helm_form_factor(self, E_r):
        """Helm nuclear form factor squared."""
        q = np.sqrt(2 * self.mA_GeV * E_r * 1e-6)  # GeV
        s = 0.9 / 197.3e-3  # fm -> GeV^-1 (0.9 fm skin thickness)
        r_n = np.sqrt((1.23 * self.A**(1/3) - 0.6)**2 + 7/3 * (0.52*np.pi)**2 - 5 * 0.9**2)
        r_n_gev = r_n / 197.3e-3  # fm -> GeV^-1
        qr = q * r_n_gev
        if qr < 1e-10:
            return 1.0
        return (3 * (np.sin(qr) - qr*np.cos(qr)) / qr**3)**2 * np.exp(-(q*s)**2)

    def differential_rate(self, E_r, m_chi, sigma_SI):
        """
        dR/dE_r in counts/keV/kg/day.
        m_chi: DM mass in GeV
        sigma_SI: SI cross-section in cm^2
        E_r: recoil energy in keV
        """
        mu_n = self.reduced_mass_nucleon(m_chi)
        mu_A = self.reduced_mass(m_chi)

        # Coherent enhancement: WIMP-nucleus cross-section
        sigma_A = sigma_SI * (mu_A / mu_n)**2 * self.A**2

        v_min_val = self.v_min(E_r, m_chi)
        eta_km = self.halo.mean_inverse_speed(v_min_val)  # (km/s)^-1

        # Number of target nuclei per gram
        m_nucleus_g = self.A * 1.66054e-24  # grams per nucleus
        N_T = 1.0 / m_nucleus_g  # nuclei per gram

        # DM number density
        n_chi = self.halo.rho / m_chi  # cm^-3

        # eta in natural units (c=1): eta_nat = eta_km * c_km
        eta_nat = eta_km * 3.0e5

        # dR/dE_R [events/(GeV g s)] = N_T * n_chi * c_cgs * sigma_A * F^2
        #                                * m_A/(2*mu_A^2) * eta_nat
        c_cgs = 3.0e10  # cm/s
        F2 = self.helm_form_factor(E_r)
        diff_factor = self.mA_GeV / (2.0 * mu_A**2)  # 1/GeV

        rate_GeV_g_s = N_T * n_chi * c_cgs * sigma_A * F2 * diff_factor * eta_nat

        # Convert to events/(keV kg day): * 86.4
        rate = rate_GeV_g_s * 86.4

        return max(rate, 0.0)

    def total_rate(self, m_chi, sigma_SI, E_thr=1.0, E_max=100.0, n_bins=200):
        """Total rate in counts/kg/day above threshold."""
        energies = np.linspace(E_thr, E_max, n_bins)
        rates = np.array([self.differential_rate(E, m_chi, sigma_SI) for E in energies])
        return np.trapz(rates, energies)

    def exclusion_limit(self, m_chi_values, n_observed=0, E_thr=1.0, CL=0.9):
        """
        Calculate 90% CL exclusion limit on sigma_SI.
        """
        upper_limit_counts = 2.3 if n_observed == 0 else stats.poisson.isf(1-CL, n_observed)
        limits = []
        for m_chi in m_chi_values:
            rate_per_sigma = self.total_rate(m_chi, 1e-45, E_thr)
            if rate_per_sigma > 0:
                sigma_limit = upper_limit_counts / (rate_per_sigma * self.exposure * 365.25) * 1e-45
            else:
                sigma_limit = 1e-30
            limits.append(sigma_limit)
        return np.array(limits)


# ============================================================
# Module 3: Non-WIMP DM Candidates
# ============================================================
class AxionDetector:
    """Axion/ALP detection sensitivity via axio-electric effect."""

    def __init__(self, target='Xe', exposure_kg_yr=1000.0):
        self.target = target
        self.exposure = exposure_kg_yr
        Z, A, mA = WIMPDetector.TARGETS[target]
        self.Z = Z
        self.A = A

    def axioelectric_rate(self, m_axion_eV, g_ae, sigma_pe_barns=1e5):
        """
        Axio-electric absorption rate (counts/kg/year).
        m_axion_eV: axion mass in eV
        g_ae: axion-electron coupling
        sigma_pe_barns: photoelectric cross-section in barns at E=m_a
        """
        rho_a = 0.3  # GeV/cm^3 local DM density
        m_a_GeV = m_axion_eV * 1e-9
        v_a = 220.0  # km/s, ~7.3e-4 c

        # Rate ~ (rho_a / m_a) * g_ae^2 * sigma_pe * v_a
        n_a = rho_a / m_a_GeV  # number density in 1/cm^3 * GeV units
        # Simplified rate
        rate = 1.2e19 * g_ae**2 * (sigma_pe_barns / 1e5) * (1.0 / (m_axion_eV / 1e3))
        return rate * self.exposure

    def sensitivity_curve(self, mass_range_eV, target_counts=2.3):
        """Compute 90% CL sensitivity to g_ae."""
        g_ae_limits = []
        for m_a in mass_range_eV:
            rate_per_g2 = self.axioelectric_rate(m_a, 1.0) / self.exposure
            if rate_per_g2 > 0:
                g_limit = np.sqrt(target_counts / (rate_per_g2 * self.exposure))
            else:
                g_limit = 1.0
            g_ae_limits.append(g_limit)
        return np.array(g_ae_limits)


class DarkPhotonDetector:
    """Dark photon (A') absorption sensitivity."""

    def __init__(self, target='Xe', exposure_kg_yr=1000.0):
        self.target = target
        self.exposure = exposure_kg_yr

    def absorption_rate(self, m_dp_eV, kinetic_mixing):
        """Dark photon absorption rate."""
        rate = 3.6e20 * kinetic_mixing**2 * (m_dp_eV / 1e3)**(-1)
        return rate * self.exposure

    def sensitivity(self, mass_range_eV, target_counts=2.3):
        kappa_limits = []
        for m in mass_range_eV:
            rate_per_k2 = self.absorption_rate(m, 1.0) / self.exposure
            if rate_per_k2 > 0:
                kappa = np.sqrt(target_counts / (rate_per_k2 * self.exposure))
            else:
                kappa = 1.0
            kappa_limits.append(kappa)
        return np.array(kappa_limits)


class PBHDetector:
    """Primordial Black Hole detection via microlensing / gravitational effects."""

    @staticmethod
    def pbh_constraint(m_pbh_solar, f_dm=1.0):
        """
        Constraints on PBH fraction of dark matter.
        Returns allowed fraction based on mass range.
        """
        log_m = np.log10(m_pbh_solar)
        # Simplified constraint landscape
        if log_m < -16:
            return min(1.0, 10**(log_m + 17))
        elif log_m < -12:
            return 0.1
        elif log_m < -10:
            return min(1.0, 10**(-(log_m + 10)))
        elif log_m < 1:
            return 0.01 + 0.05 * np.abs(np.sin(log_m * np.pi))
        else:
            return min(1.0, 0.3 * 10**(log_m - 1))


# ============================================================
# Module 4: Directional Sensitivity (CYGNUS/MIMAC)
# ============================================================
class DirectionalDetector:
    """Directional dark matter detector sensitivity calculator."""

    def __init__(self, gas='He:SF6', pressure_torr=40, volume_m3=10.0):
        self.gas = gas
        self.pressure = pressure_torr
        self.volume = volume_m3
        # Effective target mass estimation
        if 'He' in gas:
            self.target_mass_kg = volume_m3 * pressure_torr/760 * 4.0/22.4 * 1e-3
        else:
            self.target_mass_kg = volume_m3 * pressure_torr/760 * 19.0/22.4 * 1e-3
        self.angular_resolution_deg = 20.0  # typical for gas TPC

    def directional_rate(self, E_r, m_chi, sigma_SI, cos_theta):
        """
        Directional differential rate: dR/dE_r/dOmega.
        cos_theta: angle relative to WIMP wind direction.
        """
        det = WIMPDetector('He', self.target_mass_kg)
        isotropic_rate = det.differential_rate(E_r, m_chi, sigma_SI)
        # Dipole anisotropy approximation
        directional_factor = (1 + cos_theta) / (4 * np.pi)
        return isotropic_rate * directional_factor * 2

    def angular_distribution(self, m_chi, sigma_SI, E_thr=5.0, n_angles=50):
        """Compute angular distribution of recoils."""
        cos_theta = np.linspace(-1, 1, n_angles)
        rates = []
        for ct in cos_theta:
            r = sum(self.directional_rate(E, m_chi, sigma_SI, ct)
                    for E in np.linspace(E_thr, 50, 20))
            rates.append(r)
        return cos_theta, np.array(rates)

    def discovery_reach(self, m_chi_values, exposure_yr=3.0, n_sigma=3.0):
        """
        Minimum cross-section for directional discovery.
        Uses head-tail asymmetry as discriminator.
        """
        limits = []
        for m_chi in m_chi_values:
            det = WIMPDetector('He', self.target_mass_kg)
            rate_ref = det.total_rate(m_chi, 1e-45, E_thr=5.0)
            if rate_ref > 0:
                # Need ~O(10) events for directional discrimination
                n_needed = max(10, n_sigma**2 * 4)
                sigma = n_needed / (rate_ref * self.target_mass_kg * exposure_yr * 365.25) * 1e-45
            else:
                sigma = 1e-30
            limits.append(sigma)
        return np.array(limits)


# ============================================================
# Module 5: Neutrino Floor Calculation
# ============================================================
class NeutrinoFloor:
    """
    Neutrino floor/fog calculation for various targets.
    Based on coherent elastic neutrino-nucleus scattering (CEνNS).
    """

    # Neutrino flux sources: (name, max_energy_MeV, flux_cm2_s)
    NEUTRINO_SOURCES = {
        'pp': (0.423, 5.98e10),
        '7Be_384': (0.384, 4.56e8),
        '7Be_862': (0.862, 4.56e9),
        'pep': (1.442, 1.44e8),
        '8B': (16.36, 5.46e6),
        'hep': (18.77, 7.98e3),
        'atm': (1000.0, 10.5),
        'DSNB': (50.0, 85.7),
    }

    def __init__(self, target='Xe'):
        Z, A, mA = WIMPDetector.TARGETS[target]
        self.Z = Z
        self.A = A
        self.N = A - Z
        self.mA_GeV = A * 0.9315
        self.target = target

    def cevns_cross_section(self, E_nu_MeV):
        """CEνNS cross-section in cm^2."""
        Q_w = self.N - (1 - 4*0.2312) * self.Z  # weak charge
        E_nu_GeV = E_nu_MeV * 1e-3
        sigma = (GF**2 / (4*np.pi)) * Q_w**2 * E_nu_GeV**2 * 1e-26  # simplified
        return sigma

    def neutrino_recoil_rate(self, E_r_keV, source='8B'):
        """Nuclear recoil rate from neutrino source in events/keV/kg/year."""
        E_max, flux = self.NEUTRINO_SOURCES[source]
        E_r_MeV = E_r_keV * 1e-3

        # Maximum recoil energy from neutrino of energy E_nu
        E_nu_min = np.sqrt(self.mA_GeV * E_r_MeV / 2) * 1e3  # MeV (approximate)

        if E_nu_min > E_max:
            return 0.0

        sigma = self.cevns_cross_section(E_max / 2)
        rate = flux * sigma * YEAR_SECONDS * 1e3 / self.mA_GeV
        rate *= np.exp(-E_r_keV / (2 * E_max**2 / (self.mA_GeV * 1e3)))

        return max(rate, 0.0)

    def compute_floor(self, m_chi_values, exposure_range=(1, 1e6)):
        """
        Compute neutrino floor for given target.
        Returns sigma_SI values at which neutrino BG becomes limiting.
        """
        floor = []
        for m_chi in m_chi_values:
            det = WIMPDetector(self.target)
            # Estimate neutrino background
            nu_bg = sum(self.neutrino_recoil_rate(E, src)
                        for src in ['8B', 'atm', 'hep']
                        for E in np.linspace(1, 50, 20))
            nu_bg_total = nu_bg * 50.0 / 20.0  # integrated

            # Floor: where DM signal rate equals sqrt(nu_bg) for large exposure
            if nu_bg_total > 0:
                rate_ref = det.total_rate(m_chi, 1e-45, E_thr=1.0)
                if rate_ref > 0:
                    # Systematic floor limited by neutrino flux uncertainty (~5%)
                    sigma_floor = (0.05 * nu_bg_total) / (rate_ref * 365.25) * 1e-45
                else:
                    sigma_floor = 1e-40
            else:
                sigma_floor = 1e-40
            floor.append(max(sigma_floor, 1e-52))

        return np.array(floor)


# ============================================================
# Module 6: Background Model
# ============================================================
class BackgroundModel:
    """Systematic background evaluation for various strategies."""

    BACKGROUND_SOURCES = {
        'radon': {'rate_per_kg_yr': 1e-2, 'reducible': True, 'reduction_factor': 0.01},
        'krypton': {'rate_per_kg_yr': 5e-3, 'reducible': True, 'reduction_factor': 0.001},
        'neutron': {'rate_per_kg_yr': 1e-3, 'reducible': True, 'reduction_factor': 0.1},
        'surface': {'rate_per_kg_yr': 2e-3, 'reducible': True, 'reduction_factor': 0.05},
        'neutrino_cevns': {'rate_per_kg_yr': 5e-4, 'reducible': False, 'reduction_factor': 1.0},
        'detector_noise': {'rate_per_kg_yr': 1e-4, 'reducible': True, 'reduction_factor': 0.01},
    }

    STRATEGIES = {
        'baseline': {'radon': 1.0, 'krypton': 1.0, 'neutron': 1.0,
                     'surface': 1.0, 'neutrino_cevns': 1.0, 'detector_noise': 1.0},
        'distillation': {'radon': 0.01, 'krypton': 0.001, 'neutron': 1.0,
                         'surface': 1.0, 'neutrino_cevns': 1.0, 'detector_noise': 1.0},
        'active_veto': {'radon': 0.01, 'krypton': 0.001, 'neutron': 0.1,
                        'surface': 0.05, 'neutrino_cevns': 1.0, 'detector_noise': 0.1},
        'fiducialization': {'radon': 0.1, 'krypton': 0.1, 'neutron': 0.5,
                            'surface': 0.01, 'neutrino_cevns': 1.0, 'detector_noise': 0.5},
        'combined': {'radon': 0.001, 'krypton': 0.0001, 'neutron': 0.05,
                     'surface': 0.005, 'neutrino_cevns': 1.0, 'detector_noise': 0.01},
    }

    def total_background(self, strategy='baseline', exposure_kg_yr=1000.0):
        """Total background counts for given strategy and exposure."""
        factors = self.STRATEGIES[strategy]
        total = 0
        breakdown = {}
        for src, props in self.BACKGROUND_SOURCES.items():
            rate = props['rate_per_kg_yr'] * factors[src] * exposure_kg_yr
            total += rate
            breakdown[src] = rate
        return total, breakdown

    def evaluate_all_strategies(self, exposure_kg_yr=1000.0):
        """Evaluate all background reduction strategies."""
        results = {}
        for strategy in self.STRATEGIES:
            total, breakdown = self.total_background(strategy, exposure_kg_yr)
            results[strategy] = {'total': total, 'breakdown': breakdown}
        return results


# ============================================================
# Module 7: Multi-Target Complementarity
# ============================================================
class MultiTargetAnalysis:
    """Evaluate complementarity of multiple target materials."""

    def __init__(self, targets=None, exposures=None):
        if targets is None:
            targets = ['Xe', 'Ar', 'Ge', 'NaI']
        if exposures is None:
            exposures = {'Xe': 1000, 'Ar': 3000, 'Ge': 100, 'NaI': 250}
        self.targets = targets
        self.exposures = exposures
        self.detectors = {t: WIMPDetector(t, exposures.get(t, 100))
                          for t in targets}

    def combined_exclusion(self, m_chi_values, E_thresholds=None):
        """Combined exclusion from multiple targets."""
        if E_thresholds is None:
            E_thresholds = {'Xe': 1.0, 'Ar': 10.0, 'Ge': 0.5, 'NaI': 2.0}

        individual = {}
        for t in self.targets:
            individual[t] = self.detectors[t].exclusion_limit(
                m_chi_values, E_thr=E_thresholds.get(t, 1.0))

        # Combined: take minimum (best) limit at each mass
        combined = np.full_like(m_chi_values, 1e-30, dtype=float)
        for t in self.targets:
            combined = np.minimum(combined, individual[t])

        return individual, combined

    def mass_reconstruction(self, m_true, sigma_true, n_experiments=1000):
        """
        Test mass reconstruction capability with multiple targets.
        Returns reconstructed mass distribution for each target.
        """
        results = {}
        for t in self.targets:
            det = self.detectors[t]
            # Simulate observed rates
            expected = det.total_rate(m_true, sigma_true) * det.exposure * 365.25
            observed = np.random.poisson(max(0.1, expected), n_experiments)
            # Reconstruct mass from rate (simplified)
            m_recon = m_true * (1 + 0.3 * (observed - expected) / max(expected, 1))
            results[t] = m_recon
        return results


# ============================================================
# Module 8: Annual Modulation Analysis
# ============================================================
class AnnualModulation:
    """Statistical analysis of annual modulation signal."""

    def __init__(self, target='NaI', exposure_kg_yr=250.0):
        self.target = target
        self.exposure = exposure_kg_yr
        self.det = WIMPDetector(target, exposure_kg_yr)

    def modulated_rate(self, t_days, m_chi, sigma_SI, E_thr=2.0):
        """
        Rate including annual modulation.
        t: time in days from Jan 1
        Modulation amplitude ~7% of unmodulated rate.
        """
        R0 = self.det.total_rate(m_chi, sigma_SI, E_thr)
        # Phase: maximum around June 2 (day ~152)
        omega = 2 * np.pi / 365.25
        phase = 152.0  # days
        modulation_fraction = 0.07  # ~7% for typical WIMP
        R = R0 * (1 + modulation_fraction * np.cos(omega * (t_days - phase)))
        return R

    def generate_data(self, m_chi, sigma_SI, n_years=5, bins_per_year=12):
        """Generate binned annual modulation data."""
        n_bins = n_years * bins_per_year
        bin_width = 365.25 / bins_per_year
        times = np.array([(i + 0.5) * bin_width for i in range(n_bins)])
        rates = np.array([self.modulated_rate(t % 365.25, m_chi, sigma_SI)
                          for t in times])
        expected = rates * self.exposure * bin_width / 365.25
        observed = np.random.poisson(np.maximum(expected, 0.1))
        return times, observed, expected

    def modulation_significance(self, m_chi, sigma_SI, n_years=5,
                                n_mc=50, bins_per_year=12):
        """
        Compute statistical significance of modulation detection.
        Returns: median significance in sigma units.
        """
        significances = []
        for _ in range(n_mc):
            times, observed, expected = self.generate_data(
                m_chi, sigma_SI, n_years, bins_per_year)

            # Fit: R = A + B*cos(omega*(t - phase))
            omega = 2 * np.pi / 365.25
            phase = 152.0
            cos_vals = np.cos(omega * ((times % 365.25) - phase))

            mean_rate = np.mean(observed)
            if mean_rate < 0.1:
                significances.append(0.0)
                continue

            # Test statistic: amplitude of cosine component
            amplitude = 2 * np.mean((observed - mean_rate) * cos_vals)
            sigma_amp = np.sqrt(2 * np.mean(observed)) / np.sqrt(len(observed))

            if sigma_amp > 0:
                significances.append(abs(amplitude) / sigma_amp)
            else:
                significances.append(0.0)

        return np.median(significances), np.percentile(significances, [16, 84])

    def power_vs_exposure(self, m_chi, sigma_SI, year_range=None):
        """Detection power as function of exposure time."""
        if year_range is None:
            year_range = [1, 2, 3, 5, 7]
        powers = []
        for n_yr in year_range:
            sig, _ = self.modulation_significance(m_chi, sigma_SI, n_yr, n_mc=30)
            powers.append(sig)
        return np.array(year_range), np.array(powers)


# ============================================================
# Plotting Functions
# ============================================================
def plot_exclusion_curves():
    """Fig 1: Exclusion limits for multiple targets."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    m_chi = np.logspace(0, 3, 80)

    targets = {'Xe': ('b', '-', 1000), 'Ar': ('r', '--', 3000),
               'Ge': ('g', '-.', 100), 'NaI': ('orange', ':', 250)}

    for t, (color, ls, exp) in targets.items():
        det = WIMPDetector(t, exp)
        limits = det.exclusion_limit(m_chi, E_thr=1.0 if t != 'Ar' else 10.0)
        ax.plot(m_chi, limits, color=color, ls=ls, lw=2,
                label=f'{t} ({exp} kg·yr)')

    # Neutrino floor
    nf = NeutrinoFloor('Xe')
    floor = nf.compute_floor(m_chi)
    ax.fill_between(m_chi, floor, 1e-30, alpha=0.2, color='yellow',
                    label='Neutrino floor (Xe)')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Dark Matter Mass $m_\chi$ [GeV/$c^2$]', fontsize=14)
    ax.set_ylabel(r'SI Cross-section $\sigma_{SI}$ [cm$^2$]', fontsize=14)
    ax.set_title('WIMP-Nucleon Exclusion Limits: Multi-Target Comparison', fontsize=14)
    ax.set_xlim(1, 1000)
    ax.set_ylim(1e-50, 1e-40)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'exclusion_limits.png'), dpi=150)
    plt.close()
    print("  -> exclusion_limits.png")


def plot_directional_sensitivity():
    """Fig 2: Directional detector sensitivity and angular distribution."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Left: Angular distribution
    dd = DirectionalDetector(volume_m3=10)
    for m_chi, color in [(10, 'b'), (50, 'r'), (200, 'g')]:
        cos_th, rates = dd.angular_distribution(m_chi, 1e-45, E_thr=5.0)
        rates_norm = rates / np.max(rates) if np.max(rates) > 0 else rates
        ax1.plot(np.degrees(np.arccos(cos_th)), rates_norm, color=color, lw=2,
                 label=f'$m_\\chi$ = {m_chi} GeV')

    ax1.set_xlabel('Recoil angle [degrees]', fontsize=13)
    ax1.set_ylabel('Normalized rate', fontsize=13)
    ax1.set_title('Directional Recoil Distribution', fontsize=13)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)

    # Right: Discovery reach
    m_chi = np.logspace(0.3, 3, 60)
    for vol, ls in [(10, '-'), (100, '--'), (1000, ':')]:
        dd = DirectionalDetector(volume_m3=vol)
        reach = dd.discovery_reach(m_chi, exposure_yr=3.0)
        ax2.plot(m_chi, reach, ls=ls, lw=2,
                 label=f'CYGNUS {vol} m³, 3 yr')

    # Neutrino floor
    nf = NeutrinoFloor('He')
    floor = nf.compute_floor(m_chi)
    ax2.fill_between(m_chi, floor, 1e-30, alpha=0.15, color='yellow',
                     label='ν floor (He)')

    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel(r'$m_\chi$ [GeV/$c^2$]', fontsize=13)
    ax2.set_ylabel(r'$\sigma_{SI}$ [cm$^2$]', fontsize=13)
    ax2.set_title('Directional Discovery Reach (3σ)', fontsize=13)
    ax2.set_xlim(2, 1000)
    ax2.set_ylim(1e-50, 1e-38)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'directional_sensitivity.png'), dpi=150)
    plt.close()
    print("  -> directional_sensitivity.png")


def plot_neutrino_floor():
    """Fig 3: Neutrino floor comparison across targets."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    m_chi = np.logspace(0, 3, 80)

    colors = {'Xe': 'blue', 'Ar': 'red', 'Ge': 'green', 'He': 'purple'}
    for t, c in colors.items():
        nf = NeutrinoFloor(t)
        floor = nf.compute_floor(m_chi)
        ax.plot(m_chi, floor, color=c, lw=2, label=f'{t} neutrino floor')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$m_\chi$ [GeV/$c^2$]', fontsize=14)
    ax.set_ylabel(r'$\sigma_{SI}$ [cm$^2$]', fontsize=14)
    ax.set_title('Neutrino Floor for Different Target Materials', fontsize=14)
    ax.set_xlim(1, 1000)
    ax.set_ylim(1e-52, 1e-42)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'neutrino_floor.png'), dpi=150)
    plt.close()
    print("  -> neutrino_floor.png")


def plot_background_strategies():
    """Fig 4: Background reduction strategy comparison."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    bg = BackgroundModel()
    results = bg.evaluate_all_strategies(1000.0)

    # Left: Total background bar chart
    strategies = list(results.keys())
    totals = [results[s]['total'] for s in strategies]
    colors = ['gray', 'skyblue', 'salmon', 'lightgreen', 'gold']
    bars = ax1.bar(strategies, totals, color=colors, edgecolor='black')
    ax1.set_ylabel('Total Background [counts / 1000 kg·yr]', fontsize=12)
    ax1.set_title('Background Reduction Strategies', fontsize=13)
    ax1.set_yscale('log')
    ax1.set_ylim(0.01, 100)
    for bar, val in zip(bars, totals):
        ax1.text(bar.get_x() + bar.get_width()/2, val*1.2,
                 f'{val:.2f}', ha='center', fontsize=9)
    ax1.tick_params(axis='x', rotation=30)

    # Right: Breakdown by source for combined strategy
    sources = list(bg.BACKGROUND_SOURCES.keys())
    x = np.arange(len(strategies))
    width = 0.13
    for i, src in enumerate(sources):
        vals = [results[s]['breakdown'][src] for s in strategies]
        ax2.bar(x + i*width, vals, width, label=src, alpha=0.8)

    ax2.set_ylabel('Background [counts / 1000 kg·yr]', fontsize=12)
    ax2.set_title('Background Breakdown by Source', fontsize=13)
    ax2.set_yscale('log')
    ax2.set_ylim(1e-4, 100)
    ax2.set_xticks(x + width*2.5)
    ax2.set_xticklabels(strategies, rotation=30)
    ax2.legend(fontsize=9, ncol=2)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'background_strategies.png'), dpi=150)
    plt.close()
    print("  -> background_strategies.png")


def plot_multi_target():
    """Fig 5: Multi-target complementarity."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    m_chi = np.logspace(0, 3, 80)
    mta = MultiTargetAnalysis()
    individual, combined = mta.combined_exclusion(m_chi)

    colors = {'Xe': 'blue', 'Ar': 'red', 'Ge': 'green', 'NaI': 'orange'}
    for t, c in colors.items():
        ax1.plot(m_chi, individual[t], color=c, lw=1.5, ls='--',
                 label=f'{t} only', alpha=0.7)
    ax1.plot(m_chi, combined, 'k-', lw=3, label='Combined')

    nf = NeutrinoFloor('Xe')
    floor = nf.compute_floor(m_chi)
    ax1.fill_between(m_chi, floor, 1e-30, alpha=0.15, color='yellow',
                     label='ν floor')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$m_\chi$ [GeV/$c^2$]', fontsize=13)
    ax1.set_ylabel(r'$\sigma_{SI}$ [cm$^2$]', fontsize=13)
    ax1.set_title('Multi-Target Complementarity', fontsize=13)
    ax1.set_xlim(1, 1000)
    ax1.set_ylim(1e-50, 1e-40)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Right: Mass reconstruction test
    np.random.seed(42)
    m_true = 50.0
    sigma_true = 1e-46
    recon = mta.mass_reconstruction(m_true, sigma_true)
    for t, c in colors.items():
        ax2.hist(recon[t], bins=30, color=c, alpha=0.5, label=t, density=True)
    ax2.axvline(m_true, color='k', ls='--', lw=2, label='True mass')
    ax2.set_xlabel(r'Reconstructed $m_\chi$ [GeV/$c^2$]', fontsize=13)
    ax2.set_ylabel('Probability density', fontsize=13)
    ax2.set_title(f'Mass Reconstruction ($m_\\chi$ = {m_true} GeV)', fontsize=13)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'multi_target.png'), dpi=150)
    plt.close()
    print("  -> multi_target.png")


def plot_annual_modulation():
    """Fig 6: Annual modulation analysis."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    np.random.seed(42)
    am = AnnualModulation('NaI', 250.0)
    m_chi, sigma_SI = 50.0, 5e-44

    # Top-left: Modulated rate vs time
    ax = axes[0, 0]
    t = np.linspace(0, 365.25 * 5, 500)
    rates = [am.modulated_rate(ti % 365.25, m_chi, sigma_SI) for ti in t]
    ax.plot(t / 365.25, rates, 'b-', lw=1.5)
    ax.set_xlabel('Time [years]', fontsize=12)
    ax.set_ylabel('Rate [counts/kg/day]', fontsize=12)
    ax.set_title('Expected Modulated Rate (NaI, 250 kg·yr)', fontsize=12)
    ax.grid(True, alpha=0.3)

    # Top-right: Simulated data with fit
    ax = axes[0, 1]
    times, observed, expected = am.generate_data(m_chi, sigma_SI, n_years=5)
    bin_width = 365.25 / 12
    ax.errorbar(times / 365.25, observed, yerr=np.sqrt(np.maximum(observed, 1)),
                fmt='ko', markersize=4, label='Observed')
    ax.plot(times / 365.25, expected, 'r-', lw=1.5, label='Expected')
    ax.set_xlabel('Time [years]', fontsize=12)
    ax.set_ylabel('Counts per bin', fontsize=12)
    ax.set_title('Simulated Annual Modulation Data', fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Bottom-left: Detection power vs exposure
    ax = axes[1, 0]
    for sigma, color, label in [(5e-44, 'b', '5×10⁻⁴⁴'), (1e-43, 'r', '10⁻⁴³'),
                                  (5e-43, 'g', '5×10⁻⁴³')]:
        am_test = AnnualModulation('NaI', 250.0)
        years, powers = am_test.power_vs_exposure(m_chi, sigma,
                                                   year_range=[1, 2, 3, 5, 7])
        ax.plot(years, powers, 'o-', color=color, lw=2, label=f'σ = {label} cm²')

    ax.axhline(3.0, color='gray', ls='--', label='3σ threshold')
    ax.axhline(5.0, color='gray', ls=':', label='5σ threshold')
    ax.set_xlabel('Observation time [years]', fontsize=12)
    ax.set_ylabel('Significance [σ]', fontsize=12)
    ax.set_title('Modulation Detection Power', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Bottom-right: Significance heatmap
    ax = axes[1, 1]
    masses = [10, 30, 50, 100, 200]
    sigmas = [1e-44, 5e-44, 1e-43, 5e-43, 1e-42]
    sig_matrix = np.zeros((len(masses), len(sigmas)))
    for i, m in enumerate(masses):
        for j, s in enumerate(sigmas):
            am_test = AnnualModulation('NaI', 250.0)
            sig_val, _ = am_test.modulation_significance(m, s, n_years=5, n_mc=30)
            sig_matrix[i, j] = sig_val

    im = ax.imshow(sig_matrix, aspect='auto', cmap='YlOrRd',
                   extent=[-0.5, len(sigmas)-0.5, -0.5, len(masses)-0.5])
    ax.set_xticks(range(len(sigmas)))
    ax.set_xticklabels([f'{s:.0e}' for s in sigmas], rotation=45, fontsize=8)
    ax.set_yticks(range(len(masses)))
    ax.set_yticklabels([str(m) for m in masses])
    ax.set_xlabel(r'$\sigma_{SI}$ [cm²]', fontsize=12)
    ax.set_ylabel(r'$m_\chi$ [GeV]', fontsize=12)
    ax.set_title('Modulation Significance [σ] (5 yr, NaI)', fontsize=12)
    plt.colorbar(im, ax=ax, label='Significance [σ]')

    # Add text annotations
    for i in range(len(masses)):
        for j in range(len(sigmas)):
            ax.text(j, i, f'{sig_matrix[i,j]:.1f}', ha='center', va='center',
                    fontsize=9, color='white' if sig_matrix[i,j] > 3 else 'black')

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'annual_modulation.png'), dpi=150)
    plt.close()
    print("  -> annual_modulation.png")


def plot_non_wimp_candidates():
    """Fig 7: Non-WIMP dark matter candidate sensitivity."""
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

    # Axion sensitivity
    m_axion = np.logspace(-1, 2, 60)  # eV
    for target, color in [('Xe', 'blue'), ('Ge', 'green')]:
        ad = AxionDetector(target, 1000)
        g_limits = ad.sensitivity_curve(m_axion)
        ax1.plot(m_axion, g_limits, color=color, lw=2, label=f'{target} (1 t·yr)')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Axion mass [eV]', fontsize=13)
    ax1.set_ylabel(r'$g_{ae}$ coupling', fontsize=13)
    ax1.set_title('Axion-Electron Coupling Sensitivity', fontsize=13)
    ax1.set_ylim(1e-14, 1e-10)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)

    # Dark photon sensitivity
    m_dp = np.logspace(-1, 2, 60)
    for target, color in [('Xe', 'blue'), ('Ge', 'green')]:
        dpd = DarkPhotonDetector(target, 1000)
        kappa = dpd.sensitivity(m_dp)
        ax2.plot(m_dp, kappa, color=color, lw=2, label=f'{target} (1 t·yr)')

    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Dark photon mass [eV]', fontsize=13)
    ax2.set_ylabel(r'Kinetic mixing $\kappa$', fontsize=13)
    ax2.set_title('Dark Photon Sensitivity', fontsize=13)
    ax2.set_ylim(1e-17, 1e-12)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)

    # PBH constraints
    m_pbh = np.logspace(-18, 5, 200)
    f_constraints = [PBHDetector.pbh_constraint(m) for m in m_pbh]
    ax3.plot(m_pbh, f_constraints, 'k-', lw=2)
    ax3.fill_between(m_pbh, f_constraints, 1.0, alpha=0.2, color='red',
                     label='Excluded')
    ax3.axhline(1.0, color='gray', ls='--', alpha=0.5)
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlabel(r'PBH mass [$M_\odot$]', fontsize=13)
    ax3.set_ylabel(r'$f_{PBH}$ (DM fraction)', fontsize=13)
    ax3.set_title('Primordial Black Hole Constraints', fontsize=13)
    ax3.set_ylim(1e-3, 2)
    ax3.legend(fontsize=11)
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'non_wimp_candidates.png'), dpi=150)
    plt.close()
    print("  -> non_wimp_candidates.png")


def plot_recoil_spectra():
    """Fig 8: Recoil energy spectra for different targets and DM masses."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    E_r = np.linspace(0.5, 80, 200)

    # Left: Different masses, Xe target
    det = WIMPDetector('Xe', 1000)
    for m_chi, color in [(10, 'blue'), (50, 'red'), (100, 'green'), (500, 'purple')]:
        rates = [det.differential_rate(E, m_chi, 1e-45) for E in E_r]
        ax1.plot(E_r, rates, color=color, lw=2, label=f'$m_\\chi$ = {m_chi} GeV')

    ax1.set_xlabel('Recoil energy [keV]', fontsize=13)
    ax1.set_ylabel('dR/dE [counts/keV/kg/day]', fontsize=13)
    ax1.set_title('Xe Target: Recoil Spectra', fontsize=13)
    ax1.set_yscale('log')
    ax1.set_ylim(1e-12, 1e-4)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)

    # Right: Different targets, 50 GeV WIMP
    m_chi = 50.0
    for t, color in [('Xe', 'blue'), ('Ar', 'red'), ('Ge', 'green'), ('NaI', 'orange')]:
        det = WIMPDetector(t, 1000)
        rates = [det.differential_rate(E, m_chi, 1e-45) for E in E_r]
        ax2.plot(E_r, rates, color=color, lw=2, label=f'{t}')

    ax2.set_xlabel('Recoil energy [keV]', fontsize=13)
    ax2.set_ylabel('dR/dE [counts/keV/kg/day]', fontsize=13)
    ax2.set_title(f'50 GeV WIMP: Target Comparison', fontsize=13)
    ax2.set_yscale('log')
    ax2.set_ylim(1e-12, 1e-4)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'recoil_spectra.png'), dpi=150)
    plt.close()
    print("  -> recoil_spectra.png")


# ============================================================
# Main Execution
# ============================================================
def run_all_simulations():
    """Execute complete simulation suite."""
    print("=" * 60)
    print("Dark Matter Direct Detection Simulation Framework")
    print("=" * 60)

    print("\n[1/8] Generating exclusion limit curves...")
    plot_exclusion_curves()

    print("[2/8] Computing directional sensitivity...")
    plot_directional_sensitivity()

    print("[3/8] Computing neutrino floor...")
    plot_neutrino_floor()

    print("[4/8] Evaluating background strategies...")
    plot_background_strategies()

    print("[5/8] Analyzing multi-target complementarity...")
    plot_multi_target()

    print("[6/8] Running annual modulation analysis...")
    plot_annual_modulation()

    print("[7/8] Computing non-WIMP candidate sensitivity...")
    plot_non_wimp_candidates()

    print("[8/8] Generating recoil spectra...")
    plot_recoil_spectra()

    # Print summary statistics
    print("\n" + "=" * 60)
    print("SIMULATION RESULTS SUMMARY")
    print("=" * 60)

    # Background analysis
    bg = BackgroundModel()
    results = bg.evaluate_all_strategies(1000.0)
    print("\nBackground Reduction Analysis (1000 kg·yr):")
    for strategy, data in results.items():
        print(f"  {strategy:15s}: {data['total']:.4f} counts")

    # Multi-target analysis
    m_test = np.array([10, 50, 100, 500])
    mta = MultiTargetAnalysis()
    individual, combined = mta.combined_exclusion(m_test)
    print("\nMulti-Target Exclusion Limits (σ_SI [cm²]):")
    print(f"  {'Mass [GeV]':>12s} {'Xe':>12s} {'Ar':>12s} {'Ge':>12s} {'Combined':>12s}")
    for i, m in enumerate(m_test):
        print(f"  {m:12.0f} {individual['Xe'][i]:12.2e} {individual['Ar'][i]:12.2e} "
              f"{individual['Ge'][i]:12.2e} {combined[i]:12.2e}")

    # Annual modulation
    print("\nAnnual Modulation Significance (5 yr, NaI 250 kg):")
    am = AnnualModulation('NaI', 250.0)
    for sigma in [1e-44, 1e-43, 1e-42]:
        sig, ci = am.modulation_significance(50.0, sigma, n_years=5, n_mc=50)
        print(f"  σ = {sigma:.0e} cm²: {sig:.1f}σ ({ci[0]:.1f}–{ci[1]:.1f}σ)")

    print("\n" + "=" * 60)
    print("All simulations complete. Figures saved to figures/")
    print("=" * 60)

    return results


if __name__ == '__main__':
    run_all_simulations()
