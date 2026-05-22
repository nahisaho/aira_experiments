"""
Dark matter signal models for various DM candidates.
Implements differential event rates for WIMPs, axions, dark photons,
and primordial black holes.
"""
import numpy as np
from typing import Dict, Optional, Tuple
from ..core.constants import (
    NuclearTarget, TARGETS, RHO_DM_LOCAL, V_0, V_ESC, V_EARTH,
    helm_form_factor, v_min, eta_integral, M_PROTON, ALPHA_EM,
    M_ELECTRON, C_LIGHT, HBAR_EV, K_BOLTZMANN_EV, N_AVOGADRO,
    GEV_TO_KG, PB_TO_CM2, V_SUN
)


class WIMPSignal:
    """Standard WIMP-nucleus elastic scattering signal model."""

    def __init__(self, m_dm_gev: float, sigma_si_cm2: float,
                 target: NuclearTarget, exposure_kg_day: float):
        self.m_dm = m_dm_gev
        self.sigma_si = sigma_si_cm2
        self.target = target
        self.exposure = exposure_kg_day

    def differential_rate(self, Er_kev: np.ndarray) -> np.ndarray:
        """dR/dE_r in events/(keV·kg·day).

        Uses Lewin & Smith formulation with explicit unit tracking:
        dR/dE = N_T × n_χ × m_N σ_N c² F² η / (2μ_N²) × unit_conversions

        Args:
            Er_kev: recoil energies in keV

        Returns:
            differential rate array
        """
        A = self.target.A
        m_N = self.target.mass_gev
        mu_p = (self.m_dm * M_PROTON * 1e-3) / (self.m_dm + M_PROTON * 1e-3)
        mu_N = self.target.reduced_mass_dm(self.m_dm)

        sigma_N = self.sigma_si * (mu_N / mu_p)**2 * A**2  # cm²
        N_T = N_AVOGADRO * 1e3 / A     # atoms/kg
        n_chi = RHO_DM_LOCAL / self.m_dm  # cm⁻³
        c_cm = C_LIGHT * 100            # cm/s

        rates = np.zeros_like(Er_kev, dtype=float)
        for i, Er in enumerate(Er_kev):
            q_mev = np.sqrt(2 * m_N * Er)  # GeV·keV = MeV²
            F2 = helm_form_factor(np.array([q_mev]), A)[0]

            vm = v_min(Er, self.m_dm, m_N)  # km/s
            eta_kms = eta_integral(vm)       # (km/s)⁻¹
            eta_cgs = eta_kms * 1e-5         # s/cm

            # Rate per atom [1/(GeV·s)]
            dR_atom = n_chi * sigma_N * m_N * c_cm**2 * F2 * eta_cgs / (2 * mu_N**2)
            # → per keV (×1e-6) per kg (×N_T) per day (×86400)
            rates[i] = dR_atom * 1e-6 * N_T * 86400

        return rates

    def total_events(self, Er_range: Tuple[float, float],
                     n_bins: int = 100) -> float:
        """Total expected events in energy range."""
        Er = np.linspace(Er_range[0], Er_range[1], n_bins)
        dR = self.differential_rate(Er)
        return np.trapz(dR, Er) * self.exposure


class AxionSignal:
    """Axion-electron coupling signal model (axioelectric effect).

    Based on solar axion flux and axioelectric absorption in detector.
    """

    def __init__(self, m_axion_ev: float, g_ae: float,
                 target: NuclearTarget, exposure_kg_day: float):
        self.m_a = m_axion_ev
        self.g_ae = g_ae
        self.target = target
        self.exposure = exposure_kg_day

    def solar_axion_flux(self, E_kev: np.ndarray) -> np.ndarray:
        """Primakoff solar axion flux at Earth (axions/cm²/s/keV).

        Simplified parametric form from solar model.
        """
        E_ev = E_kev * 1e3
        T_core = 1.3e7 * K_BOLTZMANN_EV  # Solar core temperature ~1.1 keV
        T_kev = T_core * 1e-3

        flux_norm = 6.02e10  # axions/cm²/s/keV (for g_a_gamma = 1e-10 GeV^-1)
        flux = flux_norm * (E_kev**2.481 / np.exp(E_kev / (0.8 * T_kev))) * \
               (self.g_ae / 1e-13)**2

        return flux

    def axioelectric_cross_section(self, E_kev: np.ndarray) -> np.ndarray:
        """Axioelectric absorption cross section (cm²).

        σ_ae ≈ σ_pe × (g_ae²/β) × (3E²/(16πα m_e²))
        Simplified using photoelectric data scaling.
        """
        Z = self.target.Z
        # Approximate photoelectric cross section (barn) scaling
        sigma_pe = 1e-24 * Z**5 * (13.6e-3 / E_kev)**3.5  # cm², rough scaling

        beta = np.sqrt(np.clip(1.0 - (self.m_a * 1e-6)**2 / (E_kev * 1e-3)**2, 0, 1))
        beta = np.where(beta < 1e-10, 1e-10, beta)

        sigma_ae = sigma_pe * (self.g_ae**2) / beta * \
                   3 * (E_kev * 1e-3)**2 / (16 * np.pi * ALPHA_EM * (M_ELECTRON * 1e-3)**2)

        return sigma_ae

    def differential_rate(self, E_kev: np.ndarray) -> np.ndarray:
        """Detection rate via axioelectric effect (events/keV/kg/day)."""
        flux = self.solar_axion_flux(E_kev)
        sigma = self.axioelectric_cross_section(E_kev)
        n_atoms = N_AVOGADRO / self.target.A * 1e3  # atoms/kg

        return flux * sigma * n_atoms * 86400  # per day

    def total_events(self, E_range: Tuple[float, float],
                     n_bins: int = 100) -> float:
        E = np.linspace(E_range[0], E_range[1], n_bins)
        dR = self.differential_rate(E)
        return np.trapz(dR, E) * self.exposure


class DarkPhotonSignal:
    """Dark photon (A') absorption signal model.

    Dark photon absorbed by target electrons, depositing full mass energy.
    """

    def __init__(self, m_dp_kev: float, kappa: float,
                 target: NuclearTarget, exposure_kg_day: float):
        self.m_dp = m_dp_kev          # dark photon mass in keV
        self.kappa = kappa            # kinetic mixing parameter
        self.target = target
        self.exposure = exposure_kg_day

    def dm_absorption_rate(self) -> float:
        """Total absorption rate (events/kg/day).

        R = (ρ_DM / m_A') × κ² × σ_pe(E=m_A') × n_atoms
        """
        rho_dm = RHO_DM_LOCAL * 1e6  # keV/cm³
        n_atoms = N_AVOGADRO / self.target.A * 1e3

        Z = self.target.Z
        sigma_pe = 1e-24 * Z**5 * (13.6e-3 / self.m_dp)**3.5

        rate = (rho_dm / self.m_dp) * self.kappa**2 * sigma_pe * n_atoms * 86400
        # Apply velocity-averaged factor
        rate *= V_0 / C_LIGHT * 1e3

        return rate

    def spectral_shape(self, E_kev: np.ndarray,
                       sigma_E: float = 0.1) -> np.ndarray:
        """Energy spectrum (Gaussian peak at m_dp with detector resolution)."""
        return np.exp(-0.5 * ((E_kev - self.m_dp) / sigma_E)**2) / \
               (sigma_E * np.sqrt(2 * np.pi))

    def differential_rate(self, E_kev: np.ndarray,
                          sigma_E: float = 0.1) -> np.ndarray:
        """dR/dE (events/keV/kg/day)."""
        R_total = self.dm_absorption_rate()
        return R_total * self.spectral_shape(E_kev, sigma_E)


class PrimordialBHSignal:
    """Primordial black hole (PBH) signal via Hawking radiation.

    PBH evaporation produces particles detectable as nuclear recoils.
    Relevant for asteroid-mass PBHs (10^{16}-10^{17} g).
    """

    def __init__(self, m_pbh_g: float, f_dm: float,
                 target: NuclearTarget, exposure_kg_day: float):
        self.m_pbh = m_pbh_g          # PBH mass in grams
        self.f_dm = f_dm              # fraction of DM as PBHs
        self.target = target
        self.exposure = exposure_kg_day

    def hawking_temperature_gev(self) -> float:
        """Hawking temperature T_H = ℏc³/(8πGM) in GeV."""
        M_kg = self.m_pbh * 1e-3
        T_H = (HBAR_EV * C_LIGHT**3) / (8 * np.pi * 6.674e-11 * M_kg)
        return T_H * 1e-9  # eV to GeV

    def hawking_flux(self, E_gev: np.ndarray) -> np.ndarray:
        """Hawking radiation particle flux (particles/GeV/s/sr).

        Simplified blackbody spectrum with graybody factor.
        """
        T_H = self.hawking_temperature_gev()

        # Graybody factor (simplified, spin-0)
        r_s = 2 * 6.674e-11 * self.m_pbh * 1e-3 / C_LIGHT**2  # Schwarzschild radius
        sigma_abs = 27 * np.pi * r_s**2 / 4  # geometric cross section

        # Blackbody spectrum
        with np.errstate(over='ignore'):
            exp_factor = np.clip(E_gev / T_H, 0, 500)
            flux = sigma_abs * E_gev**2 / (2 * np.pi**2) / (np.exp(exp_factor) - 1)

        return flux

    def neutrino_induced_recoils(self, Er_kev: np.ndarray) -> np.ndarray:
        """Nuclear recoil spectrum from PBH Hawking neutrinos.

        Neutrinos from PBH → coherent elastic neutrino-nucleus scattering.
        """
        rho_pbh = RHO_DM_LOCAL * self.f_dm  # GeV/cm³
        n_pbh = rho_pbh / (self.m_pbh * 1e-3 * GEV_TO_KG * 1e3)  # cm^-3 approx

        A = self.target.A
        Z = self.target.Z
        m_N = self.target.mass_gev

        # Weak charge
        Q_w = A - Z - Z * (1 - 4 * 0.2312)  # sin²θ_W ≈ 0.2312
        G_F_natural = 1.1664e-5  # GeV^-2

        T_H = self.hawking_temperature_gev()

        rates = np.zeros_like(Er_kev)
        for i, Er in enumerate(Er_kev):
            E_nu_min = np.sqrt(m_N * Er * 1e-6 / 2)  # GeV
            if E_nu_min > 10 * T_H:
                continue
            # Integrate neutrino spectrum from Hawking radiation
            E_nu = np.linspace(E_nu_min, min(50 * T_H, 100), 200)
            with np.errstate(over='ignore'):
                exp_f = np.clip(E_nu / T_H, 0, 500)
                spectrum = E_nu**2 / (np.exp(exp_f) + 1)  # Fermi-Dirac

            # CEνNS cross section (approximate)
            sigma_cevns = G_F_natural**2 * Q_w**2 * m_N * (1e-6 * Er) / (4 * np.pi)

            rates[i] = np.trapz(spectrum * sigma_cevns, E_nu) * n_pbh * 86400

        n_atoms = N_AVOGADRO / self.target.A * 1e3
        return rates * n_atoms


class SignalFactory:
    """Factory for creating signal models."""

    @staticmethod
    def create_wimp(m_dm: float, sigma_si: float,
                    target_name: str, exposure: float) -> WIMPSignal:
        return WIMPSignal(m_dm, sigma_si, TARGETS[target_name], exposure)

    @staticmethod
    def create_axion(m_a: float, g_ae: float,
                     target_name: str, exposure: float) -> AxionSignal:
        return AxionSignal(m_a, g_ae, TARGETS[target_name], exposure)

    @staticmethod
    def create_dark_photon(m_dp: float, kappa: float,
                           target_name: str, exposure: float) -> DarkPhotonSignal:
        return DarkPhotonSignal(m_dp, kappa, TARGETS[target_name], exposure)

    @staticmethod
    def create_pbh(m_pbh: float, f_dm: float,
                   target_name: str, exposure: float) -> PrimordialBHSignal:
        return PrimordialBHSignal(m_pbh, f_dm, TARGETS[target_name], exposure)
