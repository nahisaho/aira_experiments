"""
Core physics constants and utility functions for dark matter simulations.
Uses natural units where appropriate with conversion factors.
"""
import numpy as np
from math import erf as _erf
from dataclasses import dataclass, field
from typing import Optional, Dict, Tuple

# === Physical Constants (SI + natural unit conversions) ===
C_LIGHT = 2.99792458e8        # m/s
HBAR = 1.054571817e-34        # J·s
HBAR_EV = 6.582119569e-16     # eV·s
K_BOLTZMANN = 1.380649e-23    # J/K
K_BOLTZMANN_EV = 8.617333262e-5  # eV/K
G_NEWTON = 6.67430e-11        # m³/(kg·s²)
M_PROTON = 938.272046         # MeV/c²
M_NEUTRON = 939.565379        # MeV/c²
M_ELECTRON = 0.510999         # MeV/c²
ALPHA_EM = 1.0 / 137.036      # fine structure constant
G_FERMI = 1.1663787e-5        # GeV^-2
N_AVOGADRO = 6.02214076e23    # mol^-1

# Conversion factors
GEV_TO_KG = 1.78266192e-27
KG_TO_GEV = 1.0 / GEV_TO_KG
CM2_TO_PB = 1e36              # cm² to pb
PB_TO_CM2 = 1e-36
KEV_TO_JOULE = 1.602176634e-16
GEV_TO_MEV = 1e3
TEV_TO_GEV = 1e3

# Astrophysical parameters
RHO_DM_LOCAL = 0.3            # GeV/cm³ (local DM density)
V_0 = 220.0                   # km/s (circular velocity)
V_ESC = 544.0                 # km/s (escape velocity)
V_EARTH = 232.0               # km/s (Earth velocity, annual average)
V_SUN = 232.0                 # km/s (Sun velocity in galactic frame)

# Standard Halo Model parameters
SIGMA_V = V_0 / np.sqrt(2)    # velocity dispersion


@dataclass
class NuclearTarget:
    """Properties of a nuclear target for DM scattering."""
    name: str
    Z: int                     # atomic number
    A: int                     # mass number
    mass_gev: float            # nuclear mass in GeV
    abundance: float           # natural abundance (fraction)
    spin: float = 0.0          # nuclear spin
    sp: float = 0.0            # proton spin expectation value <Sp>
    sn: float = 0.0            # neutron spin expectation value <Sn>

    @property
    def reduced_mass(self):
        """Reduced mass with proton (GeV)."""
        mp = M_PROTON * 1e-3  # GeV
        return (self.mass_gev * mp) / (self.mass_gev + mp)

    def reduced_mass_dm(self, m_dm_gev: float) -> float:
        """Reduced mass with dark matter particle (GeV)."""
        return (self.mass_gev * m_dm_gev) / (self.mass_gev + m_dm_gev)


# Standard nuclear targets
TARGETS = {
    'Xe131': NuclearTarget('Xe-131', 54, 131, 121.99, 0.212, 1.5, -0.009, -0.272),
    'Xe129': NuclearTarget('Xe-129', 54, 129, 120.12, 0.264, 0.5, 0.010, 0.329),
    'Ar40':  NuclearTarget('Ar-40',  18, 40,  37.21,  0.996, 0.0, 0.0, 0.0),
    'Ge76':  NuclearTarget('Ge-76',  32, 76,  70.73,  0.078, 0.0, 0.0, 0.0),
    'Ge73':  NuclearTarget('Ge-73',  32, 73,  67.93,  0.076, 4.5, 0.030, 0.378),
    'Na23':  NuclearTarget('Na-23',  11, 23,  21.41,  1.000, 1.5, 0.248, 0.020),
    'I127':  NuclearTarget('I-127',  53, 127, 118.21, 1.000, 2.5, 0.309, 0.075),
    'F19':   NuclearTarget('F-19',   9,  19,  17.70,  1.000, 0.5, -0.109, 0.469),
    'CF4':   NuclearTarget('CF4-F',  9,  19,  17.70,  1.000, 0.5, -0.109, 0.469),
}


def helm_form_factor(q_mev: np.ndarray, A: int) -> np.ndarray:
    """Helm nuclear form factor F²(q).

    Args:
        q_mev: momentum transfer in MeV
        A: mass number

    Returns:
        F²(q) array
    """
    hbarc = 197.3269804  # MeV·fm
    s = 0.9  # fm (skin thickness)
    r_n = np.sqrt((1.23 * A**(1./3.) - 0.6)**2 + 7./3. * (np.pi * 0.52)**2 - 5 * s**2)

    qr = q_mev * r_n / hbarc
    qs = q_mev * s / hbarc

    # Avoid division by zero
    with np.errstate(divide='ignore', invalid='ignore'):
        f = np.where(
            qr < 1e-6,
            1.0,
            3.0 * (np.sin(qr) - qr * np.cos(qr)) / qr**3
        )
    return (f * np.exp(-0.5 * qs**2))**2


def maxwell_boltzmann_velocity(v: np.ndarray, v0: float = V_0,
                                v_esc: float = V_ESC) -> np.ndarray:
    """Truncated Maxwell-Boltzmann velocity distribution in galactic frame.

    Args:
        v: velocities in km/s
        v0: characteristic velocity
        v_esc: escape velocity

    Returns:
        f(v) (unnormalized)
    """
    N_esc = (_erf(v_esc / v0) -
             2 * v_esc / (np.sqrt(np.pi) * v0) * np.exp(-(v_esc / v0)**2))
    fv = np.where(
        v < v_esc,
        (1.0 / N_esc) * (np.pi * v0**2)**(-1.5) * 4 * np.pi * v**2 *
        np.exp(-v**2 / v0**2),
        0.0
    )
    return fv


def v_min(Er_kev: float, m_dm_gev: float, m_nucleus_gev: float) -> float:
    """Minimum velocity for a given recoil energy (km/s).

    Args:
        Er_kev: recoil energy in keV
        m_dm_gev: DM mass in GeV
        m_nucleus_gev: nucleus mass in GeV

    Returns:
        v_min in km/s
    """
    Er_gev = Er_kev * 1e-6
    mu = (m_dm_gev * m_nucleus_gev) / (m_dm_gev + m_nucleus_gev)
    v = np.sqrt(m_nucleus_gev * Er_gev / (2 * mu**2))  # natural units
    return v * C_LIGHT * 1e-3  # km/s


def eta_integral(vmin_kms: float, v_e: float = V_EARTH,
                 v0: float = V_0, v_esc: float = V_ESC) -> float:
    """Mean inverse speed integral η(v_min) for SHM.

    Args:
        vmin_kms: minimum velocity in km/s
        v_e: Earth velocity
        v0: characteristic velocity
        v_esc: escape velocity

    Returns:
        η in s/km
    """
    x = vmin_kms / v0
    y = v_e / v0
    z = v_esc / v0

    N = _erf(z) - 2 * z * np.exp(-z**2) / np.sqrt(np.pi)

    if x > y + z:
        return 0.0
    elif x < z - y:
        eta = (_erf(x + y) - _erf(x - y) -
               4 * y * np.exp(-z**2) / np.sqrt(np.pi))
        return eta / (2 * N * v0 * y)
    else:
        eta = (_erf(z) - _erf(x - y) -
               2 * (y + z - x) * np.exp(-z**2) / np.sqrt(np.pi))
        return eta / (2 * N * v0 * y)
