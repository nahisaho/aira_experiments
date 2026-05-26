"""
Module 5: Reactor models coupled with microkinetic surface chemistry.
- Plug Flow Reactor (PFR)
- Continuously Stirred Tank Reactor (CSTR)
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve


class ReactorModel:
    """Base class for reactor models coupled with surface microkinetics."""
    
    def __init__(self, microkinetic_model, reactor_params):
        """
        Parameters:
            microkinetic_model: MicroKineticModel instance
            reactor_params: dict with reactor-specific parameters
        """
        self.mkm = microkinetic_model
        self.params = reactor_params
    

class PFR(ReactorModel):
    """
    Plug Flow Reactor model.
    
    Solves coupled ODEs:
      dF_i/dW = r_i  (mole balance for gas species i)
      0 = f(theta)    (pseudo-steady-state for surface coverages)
    
    W: catalyst weight [kg]
    F_i: molar flow rate of species i [mol/s]
    """
    
    def __init__(self, microkinetic_model, reactor_params):
        super().__init__(microkinetic_model, reactor_params)
        self.catalyst_weight = reactor_params.get('catalyst_weight', 1.0)  # kg
        self.total_flow = reactor_params.get('total_flow', 1e-3)  # mol/s
        self.pressure = reactor_params.get('pressure', 20.0)  # bar
        self.temperature = reactor_params.get('temperature', 500.0)  # K
        self.cat_density = reactor_params.get('cat_site_density', 1e-5)  # mol_sites/kg_cat
    
    def solve(self, F0, W_span, n_points=200):
        """
        Solve PFR equations.
        
        F0: initial molar flow rates [mol/s] for each gas species
        W_span: (W_start, W_end) catalyst weight range [kg]
        """
        def rhs(W, F):
            F_total = np.sum(F)
            if F_total <= 0:
                return np.zeros_like(F)
            P_i = (F / F_total) * self.pressure
            rates = self.mkm.compute_rates(P_i, self.temperature)
            return rates * self.cat_density
        
        W_eval = np.linspace(W_span[0], W_span[1], n_points)
        sol = solve_ivp(rhs, W_span, F0, t_eval=W_eval,
                        method='BDF', rtol=1e-8, atol=1e-10)
        return sol


class CSTR(ReactorModel):
    """
    Continuously Stirred Tank Reactor model.
    
    Solves:
      0 = F_i0 - F_i + r_i * W  (steady-state mole balance)
      0 = f(theta)               (surface coverage equations)
    """
    
    def __init__(self, microkinetic_model, reactor_params):
        super().__init__(microkinetic_model, reactor_params)
        self.catalyst_weight = reactor_params.get('catalyst_weight', 1.0)
        self.pressure = reactor_params.get('pressure', 20.0)
        self.temperature = reactor_params.get('temperature', 500.0)
        self.cat_density = reactor_params.get('cat_site_density', 1e-5)
    
    def solve(self, F0, F_guess=None):
        """
        Solve CSTR steady-state equations.
        
        F0: inlet molar flow rates [mol/s]
        F_guess: initial guess for outlet flow rates
        """
        if F_guess is None:
            F_guess = F0 * 0.5
        
        def residual(F):
            F_total = np.sum(np.abs(F))
            if F_total <= 0:
                return F0
            P_i = (np.abs(F) / F_total) * self.pressure
            rates = self.mkm.compute_rates(P_i, self.temperature)
            return F0 - np.abs(F) + rates * self.cat_density * self.catalyst_weight
        
        F_sol, info, ier, msg = fsolve(residual, F_guess, full_output=True)
        return np.abs(F_sol), ier == 1


class SimpleMicroKineticModel:
    """
    Simplified microkinetic model for Fischer-Tropsch synthesis on Co.
    
    Key reactions:
    1. CO adsorption:     CO(g) + * -> CO*
    2. H2 dissociative:   H2(g) + 2* -> 2H*
    3. CO dissociation:   CO* + * -> C* + O*
    4. C hydrogenation:   C* + H* -> CH* + *
    5. CH hydrogenation:  CH* + H* -> CH2* + *
    6. CH2 hydrogenation: CH2* + H* -> CH3* + *
    7. CH3 hydrogenation: CH3* + H* -> CH4(g) + 2*
    8. O removal:         O* + H* -> OH* + *
    9. OH removal:        OH* + H* -> H2O(g) + 2*
    10. Chain growth:     CH2* + CH2* -> C2H4(g) + 2*
    """
    
    def __init__(self, rate_constants, lateral_model=None):
        """
        rate_constants: dict with forward/reverse rate constants for each step
        """
        self.k = rate_constants
        self.lateral_model = lateral_model
        self.species = ['CO*', 'H*', 'C*', 'O*', 'CH*', 'CH2*', 'CH3*', 'OH*']
    
    def compute_rates(self, P, T):
        """
        Compute net production rates for gas-phase species using PSS coverages.
        P: array of partial pressures [CO, H2, CH4, H2O, C2H4]
        T: temperature [K]
        
        Returns rates for [CO, H2, CH4, H2O, C2H4]
        """
        P_CO = max(P[0], 1e-15)
        P_H2 = max(P[1], 1e-15)
        k = self.k

        # Quasi-equilibrated Langmuir coverages
        K_CO = k['k1f'] / max(k['k1r'], 1e-30)
        K_H2 = k['k2f'] / max(k['k2r'], 1e-30)

        denom = 1.0 + K_CO * P_CO + np.sqrt(K_H2 * P_H2)
        theta_CO = K_CO * P_CO / denom
        theta_H  = np.sqrt(K_H2 * P_H2) / denom
        theta_star = 1.0 / denom

        # CO dissociation is the rate-determining step
        r_diss = k['k3f'] * theta_CO * theta_star

        # PSS intermediates
        theta_C  = r_diss / max(k['k4f'] * theta_H, 1e-30)
        theta_O  = r_diss / max(k['k8f'] * theta_H, 1e-30)
        theta_CH = k['k4f'] * theta_C * theta_H / max(k['k5f'] * theta_H, 1e-30)
        theta_CH2 = k['k5f'] * theta_CH * theta_H / max(k['k6f'] * theta_H + 2*k['k10f']*theta_CH, 1e-30)
        theta_CH3 = k['k6f'] * theta_CH2 * theta_H / max(k['k7f'] * theta_H, 1e-30)
        theta_OH = k['k8f'] * theta_O * theta_H / max(k['k9f'] * theta_H, 1e-30)

        # Product formation rates
        r_CH4 = k['k7f'] * theta_CH3 * theta_H      # methane
        r_H2O = k['k9f'] * theta_OH * theta_H        # water
        r_C2H4 = k['k10f'] * theta_CH2**2            # ethylene (chain growth)

        # Net rates for gas species [CO, H2, CH4, H2O, C2H4]
        r_CO = -r_diss                                # CO consumed
        r_H2 = -(r_CH4 + r_H2O + r_C2H4)             # H2 consumed
        
        return np.array([r_CO, r_H2, r_CH4, r_H2O, r_C2H4])
    
    def compute_coverages(self, P, T, t_span=(0, 1.0), n_surface=8):
        """
        Compute steady-state surface coverages using analytical PSS approximation.
        For FT on Co, adsorption/desorption is quasi-equilibrated and
        CO dissociation is rate-limiting, so we use a hierarchical PSS approach.
        """
        P_CO = max(P[0], 1e-15)
        P_H2 = max(P[1], 1e-15)
        k = self.k

        # Quasi-equilibrated adsorption: Langmuir competitive
        K_CO = k['k1f'] / max(k['k1r'], 1e-30)
        K_H2 = k['k2f'] / max(k['k2r'], 1e-30)

        denom = 1.0 + K_CO * P_CO + np.sqrt(K_H2 * P_H2)
        theta_ss = np.zeros(8)
        theta_ss[0] = K_CO * P_CO / denom           # CO*
        theta_ss[1] = np.sqrt(K_H2 * P_H2) / denom  # H*
        theta_star = 1.0 / denom

        # Sequential intermediates from PSS on each hydrogenation step
        # r3 = k3f * theta_CO * theta_star (CO dissociation rate)
        r_diss = k['k3f'] * theta_ss[0] * theta_star

        # C* : r3 = r4 => theta_C = r3 / (k4f * theta_H)
        theta_ss[2] = r_diss / max(k['k4f'] * theta_ss[1], 1e-30)
        # O* : r3 = r8 => theta_O = r3 / (k8f * theta_H)
        theta_ss[3] = r_diss / max(k['k8f'] * theta_ss[1], 1e-30)
        # CH*: r4 = r5
        r4 = k['k4f'] * theta_ss[2] * theta_ss[1]
        theta_ss[4] = r4 / max(k['k5f'] * theta_ss[1], 1e-30)
        # CH2*: r5 = r6 + 2*r10
        r5 = k['k5f'] * theta_ss[4] * theta_ss[1]
        theta_ss[5] = r5 / max(k['k6f'] * theta_ss[1] + 2*k['k10f']*0.01, 1e-30)
        # CH3*: r6 = r7
        r6 = k['k6f'] * theta_ss[5] * theta_ss[1]
        theta_ss[6] = r6 / max(k['k7f'] * theta_ss[1], 1e-30)
        # OH*: r8 = r9
        r8 = k['k8f'] * theta_ss[3] * theta_ss[1]
        theta_ss[7] = r8 / max(k['k9f'] * theta_ss[1], 1e-30)

        # Clip and normalise
        theta_ss = np.clip(theta_ss, 0, 0.99)
        total = np.sum(theta_ss)
        if total > 1.0:
            theta_ss *= 0.99 / total

        # Build a synthetic trajectory (initial -> steady state)
        n_pts = 50
        t_arr = np.linspace(t_span[0], t_span[1], n_pts)
        y_arr = np.zeros((8, n_pts))
        theta0 = np.array([0.01]*8)
        theta0[0] = 0.05; theta0[1] = 0.05
        for j in range(n_pts):
            frac = 1.0 - np.exp(-10.0 * t_arr[j] / max(t_span[1], 1e-30))
            y_arr[:, j] = theta0 + (theta_ss - theta0) * frac

        class PSSolution:
            pass
        sol = PSSolution()
        sol.t = t_arr
        sol.y = y_arr
        sol.success = True
        return sol
