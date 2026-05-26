"""
Module 1: DFT-based rate constant calculation using Transition State Theory
with Wigner tunneling correction.
"""
import numpy as np

# Physical constants
kB = 1.380649e-23      # Boltzmann constant [J/K]
h = 6.62607015e-34     # Planck constant [J·s]
R = 8.314462           # Gas constant [J/(mol·K)]
eV_to_J = 1.602176634e-19  # eV to Joule conversion


def eyring_rate_constant(Ea, T, A_prefactor=None):
    """
    Eyring equation (TST):
      k = (kB*T/h) * exp(-Ea/(kB*T))
    Ea in eV, T in K.
    """
    Ea_J = Ea * eV_to_J
    k_tst = (kB * T / h) * np.exp(-Ea_J / (kB * T))
    return k_tst


def wigner_tunneling_correction(nu_imag, T):
    """
    Wigner tunneling correction factor:
      kappa = 1 + (1/24) * (h*nu / (kB*T))^2
    nu_imag: imaginary frequency of TS in cm^-1 (positive value).
    """
    nu_Hz = nu_imag * 2.998e10  # cm^-1 to Hz
    u = h * nu_Hz / (kB * T)
    kappa = 1.0 + (1.0 / 24.0) * u**2
    return kappa


def eckart_tunneling_correction(Ea_fwd, Ea_rev, nu_imag, T):
    """
    Simplified Eckart tunneling correction (asymmetric barrier).
    Uses an analytical approximation for the 1D Eckart barrier.
    """
    Ef = Ea_fwd * eV_to_J
    Er = Ea_rev * eV_to_J
    nu_Hz = nu_imag * 2.998e10
    u = h * nu_Hz / (kB * T)
    alpha = 2 * np.pi * Ef / (h * nu_Hz)
    beta = 2 * np.pi * Er / (h * nu_Hz)
    if alpha < 0.01 or beta < 0.01:
        return wigner_tunneling_correction(nu_imag, T)
    # Approximate Eckart transmission
    kappa = 1.0 + (1.0 / 24.0) * u**2 + (7.0 / 5760.0) * u**4
    return max(kappa, 1.0)


def compute_rate_constant(Ea, T, nu_imag=None, Ea_rev=None, tunneling='wigner'):
    """
    Compute rate constant with optional tunneling correction.
    
    Parameters:
        Ea: Forward activation energy [eV]
        T: Temperature [K]
        nu_imag: Imaginary frequency of transition state [cm^-1]
        Ea_rev: Reverse activation energy [eV] (for Eckart)
        tunneling: 'wigner', 'eckart', or None
    
    Returns:
        k: Rate constant [s^-1]
    """
    k_tst = eyring_rate_constant(Ea, T)
    
    if tunneling == 'wigner' and nu_imag is not None:
        kappa = wigner_tunneling_correction(nu_imag, T)
    elif tunneling == 'eckart' and nu_imag is not None and Ea_rev is not None:
        kappa = eckart_tunneling_correction(Ea, Ea_rev, nu_imag, T)
    else:
        kappa = 1.0
    
    return kappa * k_tst


def compute_equilibrium_constant(dG, T):
    """
    Compute equilibrium constant from Gibbs free energy change.
    dG in eV, T in K.
    """
    dG_J = dG * eV_to_J
    return np.exp(-dG_J / (kB * T))


def arrhenius_parameters(Ea, T_range):
    """
    Compute Arrhenius plot data for a range of temperatures.
    Returns 1/T and ln(k).
    """
    inv_T = 1.0 / T_range
    ln_k = np.log(kB * T_range / h) - Ea * eV_to_J / (kB * T_range)
    return inv_T, ln_k
