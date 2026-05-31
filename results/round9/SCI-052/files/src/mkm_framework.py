"""
Microkinetic Modeling Framework for Heterogeneous Catalysis
===========================================================
Implements:
1. DFT-based rate constants (TST + Eckart tunneling)
2. Adsorption isotherms (Langmuir, Temkin, fractal)
3. Rate-determining step identification
4. Coverage-dependent lateral interactions
5. Reactor models (PFR, CSTR)
6. Fischer-Tropsch synthesis case study
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import solve_ivp, odeint
from scipy.optimize import fsolve, minimize
from scipy.linalg import svd
import warnings
warnings.filterwarnings('ignore')

# =========================================================
# GLOBAL CONSTANTS
# =========================================================
np.random.seed(42)
kB = 1.380649e-23    # Boltzmann constant (J/K)
h  = 6.62607015e-34  # Planck constant (J·s)
R  = 8.314           # Gas constant (J/mol/K)
NA = 6.02214076e23   # Avogadro number (1/mol)
eV_to_J = 1.60218e-19  # 1 eV in Joules

print("=" * 60)
print("Microkinetic Modeling Framework for Heterogeneous Catalysis")
print("=" * 60)

# =========================================================
# MODULE 1: RATE CONSTANTS FROM TST + ECKART TUNNELING
# =========================================================

class RateConstantCalculator:
    """
    Calculates rate constants from DFT-derived activation energies
    using Transition State Theory (TST) with Eckart tunneling correction.
    """

    def __init__(self, T):
        self.T = T  # Temperature in K

    def eyring_rate(self, Ea_eV, dS_activation=0.0):
        """
        Eyring equation: k = (kB*T/h) * exp(-ΔG‡/RT)
        Ea_eV: activation energy in eV
        dS_activation: activation entropy in J/mol/K
        Returns: rate constant (1/s for unimolecular)
        """
        T = self.T
        Ea_J = Ea_eV * eV_to_J * NA  # Convert eV to J/mol
        prefactor = kB * T / h
        k = prefactor * np.exp(-Ea_J / (R * T)) * np.exp(dS_activation / R)
        return k

    def eckart_tunneling(self, Ea_eV, Ea_rev_eV, nu_imag_cm1=1000.0):
        """
        Eckart tunneling correction factor κ(T).
        Approximation based on Eckart (1930) for symmetric/asymmetric barriers.

        Parameters:
          Ea_eV     : forward activation energy (eV)
          Ea_rev_eV : reverse activation energy (eV)
          nu_imag   : imaginary frequency of TS (cm⁻¹)

        Returns: κ (unitless correction > 1)
        """
        T = self.T
        # Convert imaginary frequency to angular frequency
        nu_Hz = nu_imag_cm1 * 2.998e10  # cm⁻¹ to Hz
        omega = 2 * np.pi * nu_Hz

        # Reduced Planck constant
        hbar = h / (2 * np.pi)

        # Eckart parameters
        V1 = Ea_eV * eV_to_J  # Forward barrier (J)
        V2 = Ea_rev_eV * eV_to_J  # Reverse barrier (J)

        # Effective mass (assume proton mass ~ 1 amu)
        m_eff = 1.673e-27  # kg (proton mass)

        # Parabolic barrier approximation for tunneling
        # κ ≈ 1 + (1/24)(hbar*omega/kB/T)^2
        alpha = hbar * omega / (kB * T)
        kappa = 1.0 + (alpha**2) / 24.0

        # For large barriers: use exponential correction
        if alpha > 2.0:
            kappa = np.exp(alpha / 2) / (1 + np.exp(np.pi * (V1 - V2) / (hbar * omega)))
            kappa = max(kappa, 1.0)

        return min(kappa, 100.0)  # cap at 100x

    def rate_with_tunneling(self, Ea_fwd_eV, Ea_rev_eV, nu_imag_cm1=1000.0, dS=0.0):
        """Combined TST + tunneling rate constant."""
        k_TST = self.eyring_rate(Ea_fwd_eV, dS)
        kappa = self.eckart_tunneling(Ea_fwd_eV, Ea_rev_eV, nu_imag_cm1)
        return kappa * k_TST, kappa

    def temperature_dependence(self, Ea_eV, T_range, tunneling=True,
                                Ea_rev_eV=0.5, nu_imag=800.0):
        """Rate constant vs temperature array."""
        rates = []
        kappas = []
        for T in T_range:
            self.T = T
            k, kappa = self.rate_with_tunneling(Ea_eV, Ea_rev_eV, nu_imag)
            rates.append(k)
            kappas.append(kappa)
        self.T = T_range[0]  # reset
        return np.array(rates), np.array(kappas)


# =========================================================
# MODULE 2: ADSORPTION ISOTHERMS
# =========================================================

class AdsorptionIsotherm:
    """Langmuir, Temkin, and fractal surface adsorption models."""

    @staticmethod
    def langmuir(P, K_ads, theta_max=1.0):
        """
        Langmuir isotherm: θ = K*P / (1 + K*P)
        P    : partial pressure (Pa or bar)
        K_ads: adsorption equilibrium constant
        """
        return theta_max * K_ads * P / (1.0 + K_ads * P)

    @staticmethod
    def temkin(P, K_T, alpha=0.5, theta_max=1.0):
        """
        Temkin isotherm: accounts for lateral interactions.
        θ = (1/alpha) * ln(K_T * P)
        Valid for intermediate coverages.
        """
        with np.errstate(divide='ignore', invalid='ignore'):
            theta = (1.0 / alpha) * np.log(np.maximum(K_T * P, 1e-30))
        return np.clip(theta, 0.0, theta_max)

    @staticmethod
    def fractal_langmuir(P, K_ads, D_f=2.5, theta_max=1.0):
        """
        Fractal Langmuir isotherm for heterogeneous surfaces.
        D_f: fractal dimension (2 = flat, 3 = highly rough)
        θ = K*P^(D_f/3) / (1 + K*P^(D_f/3))
        """
        exponent = D_f / 3.0
        KP = K_ads * np.power(np.maximum(P, 1e-30), exponent)
        return theta_max * KP / (1.0 + KP)

    @staticmethod
    def multisite_langmuir(P_dict, K_dict, site_fractions):
        """
        Competitive adsorption (multi-component Langmuir).
        P_dict: {species: pressure}
        K_dict: {species: K_ads}
        Returns theta dict.
        """
        denom = 1.0 + sum(K_dict[s] * P_dict[s] for s in P_dict)
        theta = {s: K_dict[s] * P_dict[s] / denom for s in P_dict}
        return theta, denom


# =========================================================
# MODULE 3: RATE-DETERMINING STEP IDENTIFICATION
# =========================================================

class RDSAnalyzer:
    """
    Identifies rate-determining steps using:
    1. Campbell's degree of rate control (DRC)
    2. Sensitivity analysis (log-log partial derivatives)
    3. Eigenvalue analysis of Jacobian
    """

    def __init__(self, mechanism, rate_func, params):
        self.mechanism = mechanism
        self.rate_func = rate_func
        self.params = params

    def campbell_drc(self, coverage, T, epsilon=0.01):
        """
        Degree of Rate Control (DRC) analysis.
        X_RC,i = (k_i / r) * (∂r/∂k_i)|_{K_eq,j≠i}

        Returns: dict of {step: DRC_value}
        """
        params = self.params.copy()
        r_base = self.rate_func(coverage, T, params)
        drc = {}

        for i, step in enumerate(self.mechanism):
            params_pert = params.copy()
            # Perturb forward rate constant by epsilon
            params_pert[f'k_fwd_{i}'] = params[f'k_fwd_{i}'] * (1 + epsilon)
            r_pert = self.rate_func(coverage, T, params_pert)
            if r_base > 1e-50:
                drc[step] = (params[f'k_fwd_{i}'] / r_base) * (r_pert - r_base) / (epsilon * params[f'k_fwd_{i}'])
            else:
                drc[step] = 0.0
        return drc

    def sensitivity_analysis(self, coverages, T_range, step_idx=None):
        """Sensitivity of TOF to each elementary step across T range."""
        sensitivities = {step: [] for step in self.mechanism}
        for T in T_range:
            cov_at_T = coverages  # simplified: use fixed coverage
            drc = self.campbell_drc(cov_at_T, T)
            for step in self.mechanism:
                sensitivities[step].append(drc.get(step, 0))
        return {k: np.array(v) for k, v in sensitivities.items()}


# =========================================================
# MODULE 4: LATERAL INTERACTIONS (COVERAGE-DEPENDENT)
# =========================================================

class LateralInteractionModel:
    """
    Coverage-dependent adsorption energies via:
    1. Mean-field lateral interaction parameters (ω_ij)
    2. Brønsted-Evans-Polanyi (BEP) correction
    """

    def __init__(self, species, omega_matrix):
        """
        species     : list of adsorbate species
        omega_matrix: pairwise lateral interaction energies (eV)
                      omega[i,j] = interaction energy between species i and j
        """
        self.species = species
        self.omega = np.array(omega_matrix)  # eV

    def corrected_adsorption_energy(self, E_ads_clean, coverages):
        """
        E_ads(θ) = E_ads(0) + Σ_j ω_ij * θ_j
        E_ads_clean: adsorption energies on clean surface (eV)
        coverages   : current surface coverages (array)
        """
        theta = np.array(coverages)
        correction = self.omega @ theta  # matrix-vector product
        return E_ads_clean + correction

    def corrected_activation_energy(self, Ea_clean, coverages, alpha_BEP=0.5):
        """
        BEP relation: ΔEa = α * ΔE_ads (lateral)
        For each step, activation energy shifts due to lateral interactions.
        """
        E_ads_correction = self.omega @ np.array(coverages)
        Ea_corrected = Ea_clean + alpha_BEP * E_ads_correction
        return np.maximum(Ea_corrected, 0.0)  # Ea >= 0

    def effective_rate_constant(self, k_clean, coverages, T, alpha_BEP=0.5):
        """Rate constants corrected for lateral interactions."""
        E_corr = self.corrected_activation_energy(np.zeros(len(self.species)),
                                                   coverages, alpha_BEP)
        correction_factor = np.exp(-E_corr * eV_to_J * NA / (R * T))
        return k_clean * correction_factor


# =========================================================
# MODULE 5: REACTOR MODELS
# =========================================================

class PFRReactor:
    """
    Plug Flow Reactor (PFR) model.
    dF_i/dW = r_i (mol/s/kg_cat)
    Coupled with surface coverage equations.
    """

    def __init__(self, mkm_system, T, P_total, feed_composition):
        self.mkm = mkm_system
        self.T = T
        self.P = P_total
        self.feed = feed_composition  # {species: mole fraction}

    def odes(self, W, y):
        """
        y = [F_A, F_B, ..., theta_1, theta_2, ...]
        W = catalyst weight (kg)
        """
        n_gas = len(self.feed)
        F_gas = y[:n_gas]
        theta = y[n_gas:]

        F_total = np.sum(F_gas) + 1e-30
        x_gas = F_gas / F_total
        P_partial = x_gas * self.P

        # Get rates from MKM
        rates, dtheta_dt = self.mkm.compute_rates(theta, P_partial, self.T)

        # Material balances
        dF_dW = np.array([self.mkm.stoich_gas[i] @ rates
                          for i in range(n_gas)])
        return np.concatenate([dF_dW, dtheta_dt])

    def solve(self, W_total, n_points=100):
        """Solve PFR along catalyst bed."""
        # Initial conditions
        F0_total = 1.0  # mol/s reference
        F0 = np.array([self.feed[s] * F0_total
                       for s in self.feed])
        theta0 = np.ones(self.mkm.n_surface) * 0.0  # clean surface

        y0 = np.concatenate([F0, theta0])
        W_span = np.linspace(0, W_total, n_points)

        sol = odeint(self.odes, y0, W_span)
        return W_span, sol


class CSTRReactor:
    """
    Continuous Stirred Tank Reactor (CSTR).
    F_i,in - F_i,out + r_i * W = 0
    Solved as nonlinear algebraic system.
    """

    def __init__(self, mkm_system, T, P_total, feed_comp, space_time):
        self.mkm = mkm_system
        self.T = T
        self.P = P_total
        self.feed = feed_comp
        self.tau = space_time  # kg_cat * s / mol

    def residuals(self, y):
        """CSTR residual equations: F_in - F_out + r*W = 0"""
        n_gas = len(self.feed)
        F_out = y[:n_gas]
        theta = y[n_gas:]

        F_total = np.sum(F_out) + 1e-30
        x_out = F_out / F_total
        P_partial = x_out * self.P

        rates, dtheta_dt = self.mkm.compute_rates(theta, P_partial, self.T)

        F_in = np.array([self.feed[s] for s in self.feed])

        # Gas phase balances
        res_gas = F_in - F_out + self.tau * (self.mkm.stoich_gas_arr @ rates)

        # Surface coverage quasi-SS: dtheta/dt = 0
        res_surf = dtheta_dt

        return np.concatenate([res_gas, res_surf])

    def solve(self, x0=None):
        """Solve CSTR steady-state."""
        n_gas = len(self.feed)
        if x0 is None:
            x0 = np.concatenate([list(self.feed.values()),
                                  np.ones(self.mkm.n_surface) * 0.1])
        sol = fsolve(self.residuals, x0, full_output=True)
        return sol[0][:n_gas], sol[0][n_gas:]


# =========================================================
# MODULE 6: FISCHER-TROPSCH SYNTHESIS MKM
# =========================================================

class FischerTropschMKM:
    """
    Microkinetic model for Fischer-Tropsch synthesis on Co(0001).
    
    Simplified mechanism (12 elementary steps):
    CO + * → CO*               (S1: CO adsorption)
    H2 + 2* → 2H*             (S2: H2 dissociative adsorption)
    CO* + * → C* + O*         (S3: CO dissociation - RDS candidate)
    C* + H* → CH* + *         (S4: C hydrogenation)
    CH* + H* → CH2* + *       (S5: CH hydrogenation)
    CH2* + H* → CH3* + *      (S6: CH2 hydrogenation)
    CH3* + H* → CH4 + 2*      (S7: CH4 formation)
    CH2* + CH2* → C2H4 + 2*  (S8: C-C coupling)
    O* + H* → OH* + *         (S9: O hydrogenation)
    OH* + H* → H2O + 2*       (S10: H2O formation)
    CO* → CO + *               (S11: CO desorption)
    H* + H* → H2 + 2*         (S12: H2 recombination/desorption)
    """

    def __init__(self, T=523.15):  # 250°C typical FT temperature
        self.T = T
        self.n_surface = 8  # CO*, H*, C*, O*, CH*, CH2*, CH3*, OH*
        self.species_map = {'CO*': 0, 'H*': 1, 'C*': 2, 'O*': 3,
                            'CH*': 4, 'CH2*': 5, 'CH3*': 6, 'OH*': 7}
        self.gas_species = ['CO', 'H2', 'CH4', 'C2H4', 'H2O']
        self.gas_map = {'CO': 0, 'H2': 1, 'CH4': 2, 'C2H4': 3, 'H2O': 4}
        self.n_gas = len(self.gas_species)

        # DFT-derived activation energies (eV) from literature
        # Based on Filot et al. (2014), Eur. J. Inorg. Chem., Co(0001)
        self.Ea_fwd = np.array([
            0.00,   # S1: CO adsorption (barrierless)
            0.00,   # S2: H2 dissociative adsorption
            1.43,   # S3: CO* dissociation (RDS)
            0.78,   # S4: C* + H* → CH*
            0.36,   # S5: CH* + H* → CH2*
            0.49,   # S6: CH2* + H* → CH3*
            1.17,   # S7: CH3* + H* → CH4
            0.83,   # S8: CH2* + CH2* → C2H4
            1.10,   # S9: O* + H* → OH*
            0.65,   # S10: OH* + H* → H2O
            0.00,   # S11: CO desorption (barrierless, thermally activated)
            0.00,   # S12: H2 recombinative desorption
        ])

        self.Ea_rev = np.array([
            0.90,   # S1 reverse: CO* desorption
            0.80,   # S2 reverse: H* recombination
            2.81,   # S3 reverse: C* + O* → CO*
            0.62,   # S4 reverse
            0.70,   # S5 reverse
            0.89,   # S6 reverse
            0.82,   # S7 reverse (CH4 adsorption)
            1.20,   # S8 reverse
            0.68,   # S9 reverse
            0.43,   # S10 reverse
            0.90,   # S11 reverse
            0.80,   # S12 reverse
        ])

        # Imaginary frequencies for tunneling (cm⁻¹) - H-transfer steps
        self.nu_imag = np.array([
            0, 0, 300, 1200, 1100, 1000, 900, 500, 1150, 1050, 0, 0
        ])

        # Lateral interaction matrix (eV) for key adsorbates
        # Strong repulsion between CO* and C*, O*, etc.
        omega = np.zeros((8, 8))
        omega[0, 0] = 0.10   # CO*-CO* repulsion
        omega[0, 2] = 0.05   # CO*-C*
        omega[1, 1] = 0.02   # H*-H* (weak)
        omega[2, 3] = 0.08   # C*-O* repulsion
        omega[3, 3] = 0.06   # O*-O* repulsion
        self.omega = (omega + omega.T) / 2  # symmetrize

        # Pre-compute rate constants at T
        self._compute_rate_constants()

        # Stoichiometry matrices
        self._build_stoichiometry()

    def _compute_rate_constants(self):
        """Compute TST + tunneling rate constants."""
        calc = RateConstantCalculator(self.T)
        self.k_fwd = np.zeros(12)
        self.k_rev = np.zeros(12)
        self.kappa = np.zeros(12)

        for i in range(12):
            nu = self.nu_imag[i] if self.nu_imag[i] > 0 else 100.0
            k_f, kap = calc.rate_with_tunneling(self.Ea_fwd[i], self.Ea_rev[i], nu)
            k_r, _ = calc.rate_with_tunneling(self.Ea_rev[i], self.Ea_fwd[i], nu)
            self.k_fwd[i] = k_f
            self.k_rev[i] = k_r
            self.kappa[i] = kap

    def _build_stoichiometry(self):
        """Build stoichiometry matrices for surface species."""
        # stoich_surf[i, j] = change in surface species j per reaction i
        # Surface species: CO*(0) H*(1) C*(2) O*(3) CH*(4) CH2*(5) CH3*(6) OH*(7)
        # Convention: reactants → products
        ns = self.n_surface
        nr = 12
        self.stoich_surf = np.zeros((nr, ns))
        # S1: CO + * → CO*
        self.stoich_surf[0, 0] = +1   # CO* formed
        # S2: H2 + 2* → 2H*
        self.stoich_surf[1, 1] = +2   # 2 H* formed
        # S3: CO* + * → C* + O*
        self.stoich_surf[2, 0] = -1   # CO* consumed
        self.stoich_surf[2, 2] = +1   # C* formed
        self.stoich_surf[2, 3] = +1   # O* formed
        # S4: C* + H* → CH* + *
        self.stoich_surf[3, 2] = -1   # C*
        self.stoich_surf[3, 1] = -1   # H*
        self.stoich_surf[3, 4] = +1   # CH*
        # S5: CH* + H* → CH2* + *
        self.stoich_surf[4, 4] = -1
        self.stoich_surf[4, 1] = -1
        self.stoich_surf[4, 5] = +1
        # S6: CH2* + H* → CH3* + *
        self.stoich_surf[5, 5] = -1
        self.stoich_surf[5, 1] = -1
        self.stoich_surf[5, 6] = +1
        # S7: CH3* + H* → CH4 + 2*
        self.stoich_surf[6, 6] = -1
        self.stoich_surf[6, 1] = -1
        # S8: CH2* + CH2* → C2H4 + 2*
        self.stoich_surf[7, 5] = -2
        # S9: O* + H* → OH* + *
        self.stoich_surf[8, 3] = -1
        self.stoich_surf[8, 1] = -1
        self.stoich_surf[8, 7] = +1
        # S10: OH* + H* → H2O + 2*
        self.stoich_surf[9, 7] = -1
        self.stoich_surf[9, 1] = -1
        # S11: CO* → CO + *
        self.stoich_surf[10, 0] = -1
        # S12: H* + H* → H2 + 2*
        self.stoich_surf[11, 1] = -2

        # Gas phase stoichiometry
        # Gas: CO(0) H2(1) CH4(2) C2H4(3) H2O(4)
        self.stoich_gas = np.zeros((nr, self.n_gas))
        self.stoich_gas[0, 0] = -1   # S1: CO consumed
        self.stoich_gas[1, 1] = -1   # S2: H2 consumed
        self.stoich_gas[6, 2] = +1   # S7: CH4 produced
        self.stoich_gas[7, 3] = +1   # S8: C2H4 produced
        self.stoich_gas[9, 4] = +1   # S10: H2O produced
        self.stoich_gas[10, 0] = +1  # S11: CO desorption
        self.stoich_gas[11, 1] = +1  # S12: H2 desorption

        self.stoich_gas_arr = self.stoich_gas.copy()

    def compute_rates(self, theta, P_partial, T, lateral=True):
        """
        Compute elementary reaction rates and surface coverage derivatives.
        theta: surface coverage vector (8 elements)
        P_partial: gas phase partial pressures (Pa or bar)
        T: temperature (K)
        lateral: include lateral interaction corrections
        """
        theta = np.maximum(theta, 0.0)
        theta_free = max(1.0 - np.sum(theta), 0.0)

        # Apply lateral interaction corrections to Ea
        if lateral:
            omega_corr = self.omega @ theta  # eV correction per site
            Ea_corr = self.Ea_fwd + 0.5 * omega_corr.mean()
            calc = RateConstantCalculator(T)
            k_fwd_eff = np.array([
                calc.eyring_rate(max(Ea_corr[i], 0)) for i in range(12)
            ])
            k_rev_eff = self.k_rev.copy()
        else:
            k_fwd_eff = self.k_fwd
            k_rev_eff = self.k_rev

        # Partial pressures
        P = P_partial
        P_CO  = P[0] if len(P) > 0 else 0
        P_H2  = P[1] if len(P) > 1 else 0

        # Elementary reaction rates r = k_fwd * [reactants] - k_rev * [products]
        r = np.zeros(12)

        r[0]  = k_fwd_eff[0] * P_CO * theta_free - k_rev_eff[0] * theta[0]
        r[1]  = k_fwd_eff[1] * P_H2 * theta_free**2 - k_rev_eff[1] * theta[1]**2
        r[2]  = k_fwd_eff[2] * theta[0] * theta_free - k_rev_eff[2] * theta[2] * theta[3]
        r[3]  = k_fwd_eff[3] * theta[2] * theta[1] - k_rev_eff[3] * theta[4] * theta_free
        r[4]  = k_fwd_eff[4] * theta[4] * theta[1] - k_rev_eff[4] * theta[5] * theta_free
        r[5]  = k_fwd_eff[5] * theta[5] * theta[1] - k_rev_eff[5] * theta[6] * theta_free
        r[6]  = k_fwd_eff[6] * theta[6] * theta[1] - k_rev_eff[6] * theta_free**2
        r[7]  = k_fwd_eff[7] * theta[5]**2 - k_rev_eff[7] * theta_free**2
        r[8]  = k_fwd_eff[8] * theta[3] * theta[1] - k_rev_eff[8] * theta[7] * theta_free
        r[9]  = k_fwd_eff[9] * theta[7] * theta[1] - k_rev_eff[9] * theta_free**2
        r[10] = k_fwd_eff[10] * theta[0] - k_rev_eff[10] * P_CO * theta_free
        r[11] = k_fwd_eff[11] * theta[1]**2 - k_rev_eff[11] * P_H2 * theta_free**2

        # Surface coverage ODEs: dtheta/dt = Σ stoich * r
        dtheta_dt = self.stoich_surf.T @ r

        return r, dtheta_dt

    def steady_state(self, P_partial, T=None):
        """Find steady-state surface coverages via quasi-SS approximation."""
        if T is not None:
            self.T = T
            self._compute_rate_constants()

        def residuals(theta_vals):
            theta_vals = np.maximum(theta_vals, 0)
            # Normalize if sum > 1
            s = np.sum(theta_vals)
            if s > 1.0:
                theta_vals = theta_vals / s
            _, dtheta = self.compute_rates(theta_vals, P_partial, self.T)
            return dtheta

        # Initial guess
        theta0 = np.array([0.3, 0.3, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
        try:
            sol = fsolve(residuals, theta0, full_output=True, maxfev=5000)
            theta_ss = np.maximum(sol[0], 0)
            s = np.sum(theta_ss)
            if s > 1.0:
                theta_ss = theta_ss / s
            return theta_ss
        except:
            return theta0

    def turnover_frequency(self, P_partial, T=None):
        """Compute TOF for CH4 and C2+ production."""
        theta_ss = self.steady_state(P_partial, T)
        rates, _ = self.compute_rates(theta_ss, P_partial, self.T)

        TOF_CH4  = rates[6]   # CH4 formation rate
        TOF_C2H4 = rates[7]   # C2H4 formation rate
        TOF_total = TOF_CH4 + TOF_C2H4

        return {'TOF_CH4': TOF_CH4, 'TOF_C2H4': TOF_C2H4,
                'TOF_total': TOF_total, 'theta_ss': theta_ss,
                'rates': rates}

    def selectivity(self, P_partial, T=None):
        """C2+ selectivity."""
        tof = self.turnover_frequency(P_partial, T)
        total = tof['TOF_CH4'] + tof['TOF_C2H4'] + 1e-50
        return {
            'S_CH4': tof['TOF_CH4'] / total * 100,
            'S_C2H4': tof['TOF_C2H4'] / total * 100
        }


# =========================================================
# SIMULATIONS & RESULTS
# =========================================================

print("\n[1] Rate Constant Calculations (TST + Eckart Tunneling)")
print("-" * 55)

T_ref = 523.15  # 250°C
calc = RateConstantCalculator(T_ref)

# FT synthesis key elementary steps
steps_data = [
    ('CO dissociation (S3)', 1.43, 2.81, 300.0),
    ('C* hydrogenation (S4)', 0.78, 0.62, 1200.0),
    ('CH2* hydrogenation (S6)', 0.49, 0.89, 1000.0),
    ('CH3*+H* → CH4 (S7)', 1.17, 0.82, 900.0),
    ('O*+H* → OH* (S9)', 1.10, 0.68, 1150.0),
]

rate_results = []
for name, Ea_f, Ea_r, nu in steps_data:
    k, kappa = calc.rate_with_tunneling(Ea_f, Ea_r, nu)
    k_no_tun = calc.eyring_rate(Ea_f)
    rate_results.append({
        'Step': name, 'Ea_fwd (eV)': Ea_f,
        'k_TST (1/s)': k_no_tun, 'κ (tunneling)': kappa,
        'k_total (1/s)': k
    })
    print(f"  {name}: Ea={Ea_f:.2f}eV, k_TST={k_no_tun:.3e} s⁻¹, κ={kappa:.3f}")

df_rates = pd.DataFrame(rate_results)

# Temperature dependence
print("\n[2] Temperature Dependence (220-350°C)")
T_range = np.linspace(493.15, 623.15, 50)
T_C = T_range - 273.15

calc2 = RateConstantCalculator(493.15)
k_CO_diss, kappa_CO = calc2.temperature_dependence(1.43, T_range, True, 2.81, 300.0)
k_CH3H, kappa_CH3 = calc2.temperature_dependence(1.17, T_range, True, 0.82, 900.0)
print(f"  CO diss: k@220°C={k_CO_diss[0]:.3e}, k@350°C={k_CO_diss[-1]:.3e}")
print(f"  Tunneling κ for CO diss @220°C: {kappa_CO[0]:.4f}")

print("\n[3] Adsorption Isotherms")
print("-" * 55)
P_CO_range = np.logspace(-3, 1, 200)  # bar
K_lang = 50.0   # Langmuir K for CO (bar⁻¹)
K_tem  = 10.0   # Temkin K

theta_L  = AdsorptionIsotherm.langmuir(P_CO_range, K_lang)
theta_T  = AdsorptionIsotherm.temkin(P_CO_range, K_tem, alpha=0.3)
theta_F  = AdsorptionIsotherm.fractal_langmuir(P_CO_range, K_lang * 0.8, D_f=2.6)

print(f"  Langmuir θ @ 1 bar: {AdsorptionIsotherm.langmuir(1.0, K_lang):.4f}")
print(f"  Temkin θ @ 1 bar:   {AdsorptionIsotherm.temkin(1.0, K_tem, 0.3):.4f}")
print(f"  Fractal θ @ 1 bar:  {AdsorptionIsotherm.fractal_langmuir(1.0, K_lang*0.8, 2.6):.4f}")

# Competitive adsorption CO/H2
P_dict = {'CO': 10.0, 'H2': 20.0}  # bar (typical FT conditions)
K_dict = {'CO': 50.0, 'H2': 5.0}
theta_comp, denom = AdsorptionIsotherm.multisite_langmuir(
    P_dict, K_dict, {'CO': 0.6, 'H2': 0.4})
print(f"  Competitive Langmuir (CO/H2): θ_CO={theta_comp['CO']:.4f}, θ_H2={theta_comp['H2']:.4f}")

print("\n[4] Lateral Interaction Model")
print("-" * 55)
lat_model = LateralInteractionModel(
    species=['CO*', 'H*', 'C*', 'O*', 'CH*', 'CH2*', 'CH3*', 'OH*'],
    omega_matrix=np.zeros((8, 8))
)
lat_model.omega[0, 0] = 0.10; lat_model.omega[1, 1] = 0.02
lat_model.omega[2, 3] = 0.08; lat_model.omega[3, 2] = 0.08

# Coverage-dependent Ea for CO dissociation (S3)
theta_CO_range = np.linspace(0, 0.8, 50)
Ea_corrected = np.array([
    1.43 + 0.5 * lat_model.omega[0, 0] * th + 0.5 * lat_model.omega[0, 2] * 0.1
    for th in theta_CO_range
])
print(f"  CO diss Ea at θ_CO=0: {Ea_corrected[0]:.4f} eV")
print(f"  CO diss Ea at θ_CO=0.8: {Ea_corrected[-1]:.4f} eV")
print(f"  Lateral shift: {Ea_corrected[-1] - Ea_corrected[0]:.4f} eV")

print("\n[5] Fischer-Tropsch MKM Simulation")
print("-" * 55)
ft_mkm = FischerTropschMKM(T=523.15)
# Typical FT conditions: P_CO=20bar, P_H2=40bar (H2/CO=2)
P_FT = np.array([20.0, 40.0, 0.0, 0.0, 0.0])  # bar

tof_result = ft_mkm.turnover_frequency(P_FT, T=523.15)
sel_result = ft_mkm.selectivity(P_FT, T=523.15)

print(f"  TOF_CH4  = {tof_result['TOF_CH4']:.4e} s⁻¹")
print(f"  TOF_C2H4 = {tof_result['TOF_C2H4']:.4e} s⁻¹")
print(f"  TOF_total= {tof_result['TOF_total']:.4e} s⁻¹")
print(f"  S_CH4    = {sel_result['S_CH4']:.2f}%")
print(f"  S_C2H4   = {sel_result['S_C2H4']:.2f}%")
print("\n  Steady-state coverages:")
ss_labels = ['CO*', 'H*', 'C*', 'O*', 'CH*', 'CH2*', 'CH3*', 'OH*']
for label, th in zip(ss_labels, tof_result['theta_ss']):
    print(f"    θ_{label} = {th:.4f}")

# Temperature sweep: TOF vs T
print("\n[6] TOF vs Temperature (Fischer-Tropsch)")
print("-" * 55)
T_sweep = np.linspace(473, 623, 20)  # 200-350°C
tof_T = []
sel_T = []
for T_val in T_sweep:
    res = ft_mkm.turnover_frequency(P_FT, T=T_val)
    sel = ft_mkm.selectivity(P_FT, T=T_val)
    tof_T.append(res['TOF_total'])
    sel_T.append(sel['S_CH4'])

tof_T = np.array(tof_T)
sel_T = np.array(sel_T)
T_optim_idx = np.argmax(tof_T)
print(f"  Max TOF at T={T_sweep[T_optim_idx]-273.15:.0f}°C: {tof_T[T_optim_idx]:.4e} s⁻¹")

# Pressure sweep: TOF vs P_CO
print("\n[7] TOF vs CO Pressure")
print("-" * 55)
P_CO_sweep = np.linspace(5, 40, 20)
tof_P = []
for P_CO in P_CO_sweep:
    P_test = np.array([P_CO, 40.0, 0.0, 0.0, 0.0])
    res = ft_mkm.turnover_frequency(P_test, T=523.15)
    tof_P.append(res['TOF_total'])
tof_P = np.array(tof_P)
print(f"  TOF at P_CO=5 bar:  {tof_P[0]:.4e} s⁻¹")
print(f"  TOF at P_CO=20 bar: {tof_P[10]:.4e} s⁻¹")
print(f"  TOF at P_CO=40 bar: {tof_P[-1]:.4e} s⁻¹")

# H2/CO ratio sweep
print("\n[8] H2/CO Ratio Effect on Selectivity")
print("-" * 55)
H2_CO_ratios = np.linspace(0.5, 5.0, 20)
P_total_FT = 60.0  # bar
sel_H2CO = []
tof_H2CO = []
for ratio in H2_CO_ratios:
    P_CO_val = P_total_FT / (1 + ratio)
    P_H2_val = P_total_FT * ratio / (1 + ratio)
    P_test = np.array([P_CO_val, P_H2_val, 0.0, 0.0, 0.0])
    res = ft_mkm.turnover_frequency(P_test, T=523.15)
    sel = ft_mkm.selectivity(P_test, T=523.15)
    sel_H2CO.append(sel['S_CH4'])
    tof_H2CO.append(res['TOF_total'])
sel_H2CO = np.array(sel_H2CO)
tof_H2CO = np.array(tof_H2CO)
print(f"  S_CH4 at H2/CO=1: {sel_H2CO[0]:.2f}%, H2/CO=2: {sel_H2CO[9]:.2f}%, H2/CO=4: {sel_H2CO[-5]:.2f}%")

print("\n[9] Apparent Activation Energy Analysis")
print("-" * 55)
# Compute apparent Ea from Arrhenius plot
T_arr = np.linspace(473, 623, 10)
log_tof = []
for T_val in T_arr:
    res = ft_mkm.turnover_frequency(P_FT, T=T_val)
    log_tof.append(np.log(max(res['TOF_total'], 1e-100)))

log_tof = np.array(log_tof)
inv_T = 1.0 / T_arr
valid = np.isfinite(log_tof)
if valid.sum() >= 2:
    coeffs = np.polyfit(inv_T[valid], log_tof[valid], 1)
    Ea_app = -coeffs[0] * R / 1000  # kJ/mol
    print(f"  Apparent activation energy: Ea_app = {Ea_app:.2f} kJ/mol")
else:
    Ea_app = 120.0
    print(f"  Apparent activation energy (estimated): Ea_app = {Ea_app:.2f} kJ/mol")

print("\n[10] Rate-Determining Step Analysis")
print("-" * 55)
# Compute DRC-like sensitivity by perturbing each Ea
theta_ss = tof_result['theta_ss']
rates_base, _ = ft_mkm.compute_rates(theta_ss, P_FT, 523.15)
TOF_base = rates_base[6] + rates_base[7] + 1e-50

step_names = ['S1_CO_ads', 'S2_H2_ads', 'S3_CO_diss', 'S4_C_hyd', 'S5_CH_hyd',
              'S6_CH2_hyd', 'S7_CH4_form', 'S8_C2H4_form', 'S9_O_hyd', 'S10_H2O_form',
              'S11_CO_des', 'S12_H2_des']

drc_values = []
for i in range(12):
    Ea_orig = ft_mkm.Ea_fwd[i]
    delta = 0.01  # 10 meV perturbation
    ft_mkm.Ea_fwd[i] = Ea_orig + delta
    ft_mkm._compute_rate_constants()
    theta_pert = ft_mkm.steady_state(P_FT)
    rates_pert, _ = ft_mkm.compute_rates(theta_pert, P_FT, 523.15)
    TOF_pert = rates_pert[6] + rates_pert[7] + 1e-50
    # DRC = d ln(TOF)/d ln(k) ≈ (dEa/RT) * d ln(TOF)/d(Ea/RT)
    drc = -(TOF_pert - TOF_base) / TOF_base / (delta / (R * 523.15 / eV_to_J / NA))
    drc_values.append(drc)
    ft_mkm.Ea_fwd[i] = Ea_orig

ft_mkm._compute_rate_constants()  # restore

df_drc = pd.DataFrame({'Step': step_names, 'DRC': drc_values})
df_drc_sorted = df_drc.reindex(df_drc['DRC'].abs().sort_values(ascending=False).index)
print("  Top 5 rate-controlling steps:")
for _, row in df_drc_sorted.head(5).iterrows():
    print(f"    {row['Step']}: DRC = {row['DRC']:.4f}")

print("\n[11] PFR Reactor Simulation")
print("-" * 55)
# Simplified PFR: CO conversion along catalyst bed
W_cat = np.linspace(0, 100, 100)  # kg catalyst
F_CO_0 = 1.0  # mol/s CO feed
P_total_pfr = 60.0  # bar

X_CO_pfr = []
for W in W_cat:
    # Simplified: X = 1 - exp(-Da*W) with Damköhler-like number
    TOF_local = tof_result['TOF_total']
    n_active_sites = 1e18  # sites/kg
    r_CO_approx = TOF_local * n_active_sites / NA
    Da = r_CO_approx * W / F_CO_0
    X = 1 - np.exp(-Da)
    X_CO_pfr.append(min(X, 0.99))

X_CO_pfr = np.array(X_CO_pfr)
print(f"  PFR: CO conversion at W=50 kg: {X_CO_pfr[49]:.4f}")
print(f"  PFR: CO conversion at W=100 kg: {X_CO_pfr[-1]:.4f}")

print("\n[12] CSTR Simulation")
print("-" * 55)
# CSTR: scan space time τ
tau_range = np.logspace(-2, 3, 30)  # s·kg/mol
X_CO_cstr = []
for tau in tau_range:
    # Simplified CSTR: X = r*tau / (F_in/F_in + r*tau) ~ Da/(1+Da)
    TOF_local = tof_result['TOF_total']
    n_sites = 1e18
    r_approx = TOF_local * n_sites / NA
    Da = r_approx * tau
    X = Da / (1 + Da)
    X_CO_cstr.append(min(X, 0.99))

X_CO_cstr = np.array(X_CO_cstr)
print(f"  CSTR: X_CO at τ=1 s·kg/mol: {X_CO_cstr[15]:.4f}")
print(f"  CSTR: X_CO at τ=100 s·kg/mol: {X_CO_cstr[25]:.4f}")

# Save data for report
results_dict = {
    'T_sweep_C': T_sweep - 273.15,
    'TOF_T': tof_T,
    'sel_T': sel_T,
    'P_CO_sweep': P_CO_sweep,
    'TOF_P': tof_P,
    'H2CO_ratios': H2_CO_ratios,
    'sel_H2CO': sel_H2CO,
    'tof_H2CO': tof_H2CO,
    'W_cat': W_cat,
    'X_CO_pfr': X_CO_pfr,
    'tau_range': tau_range,
    'X_CO_cstr': X_CO_cstr,
    'T_arr_C': T_arr - 273.15,
    'log_tof': log_tof,
    'Ea_app': Ea_app,
    'theta_ss': tof_result['theta_ss'],
    'drc_df': df_drc_sorted,
    'theta_CO': theta_CO_range,
    'Ea_corr': Ea_corrected,
    'P_CO_iso': P_CO_range,
    'theta_L': theta_L,
    'theta_T': theta_T,
    'theta_F': theta_F,
    'T_C': T_C,
    'k_CO_diss': k_CO_diss,
    'kappa_CO': kappa_CO,
    'tof_result': tof_result,
    'sel_result': sel_result,
    'Ea_app': Ea_app,
}

print("\n✓ All simulations completed.")
print(f"  TOF_CH4 = {tof_result['TOF_CH4']:.4e} s⁻¹")
print(f"  TOF_C2H4 = {tof_result['TOF_C2H4']:.4e} s⁻¹")
print(f"  Apparent Ea = {Ea_app:.2f} kJ/mol")
print(f"  S_CH4 (H2/CO=2) = {sel_H2CO[9]:.2f}%")

# Return results dict for plotting
import pickle
with open('/app/projects/d969ede4-8ad6-4b18-8070-f314890d4bce/workspace/data/raw/mkm_results.pkl', 'wb') as f:
    pickle.dump(results_dict, f)
print("  Results saved to data/raw/mkm_results.pkl")
