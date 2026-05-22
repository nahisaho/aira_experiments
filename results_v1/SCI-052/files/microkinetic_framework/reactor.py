"""
Module 5: Reactor Models (PFR/CSTR) Coupling
=============================================
Couples microkinetic surface models with ideal reactor models.
"""

import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Tuple, Optional


@dataclass
class ReactorConditions:
    """Operating conditions for the reactor."""
    T: float                    # Temperature [K]
    P_total: float              # Total pressure [bar]
    feed_composition: Dict[str, float]   # {species: mole fraction}
    F_total: float = 1.0        # Total molar flow rate [mol/s]
    catalyst_mass: float = 1.0  # Catalyst mass [kg]
    site_density: float = 1.5e-5  # Surface site density [mol/kg_cat]
    void_fraction: float = 0.4  # Bed void fraction
    reactor_volume: float = 1.0 # Reactor volume [m^3]


@dataclass
class ReactorResult:
    """Results from reactor simulation."""
    reactor_type: str
    positions: np.ndarray        # Axial positions (PFR) or time steps
    gas_compositions: Dict[str, np.ndarray]  # Mole fractions along reactor
    surface_coverages: Dict[str, np.ndarray] # Surface coverages
    conversion: float            # Reactant conversion
    selectivities: Dict[str, float]  # Product selectivities
    turnover_frequency: float    # TOF [1/s]
    space_time_yield: float      # STY [mol/(kg_cat·s)]
    conditions: ReactorConditions


def solve_steady_state_coverages(rate_constants_list: list,
                                 gas_pressures: Dict[str, float],
                                 T: float,
                                 surface_species: List[str],
                                 stoich_matrix: np.ndarray,
                                 rate_expressions: Callable,
                                 max_iter: int = 5000,
                                 tol: float = 1e-12) -> np.ndarray:
    """
    Solve for steady-state surface coverages (dθ/dt = 0).

    Uses pseudo-transient continuation method.

    Parameters
    ----------
    rate_constants_list : list
        Rate constants for each elementary step.
    gas_pressures : dict
        {species: partial pressure [bar]}
    T : float
        Temperature [K].
    surface_species : list of str
        Surface species names.
    stoich_matrix : np.ndarray
        Stoichiometric matrix [n_species x n_reactions].
    rate_expressions : callable
        Function(rate_constants, coverages, pressures, T) -> rates array.
    max_iter : int
        Maximum iterations.
    tol : float
        Convergence tolerance.

    Returns
    -------
    np.ndarray
        Steady-state coverages.
    """
    n_species = len(surface_species)
    theta = np.ones(n_species) * 0.01  # Initial guess
    dt = 1e-3  # pseudo time step

    for iteration in range(max_iter):
        rates = rate_expressions(rate_constants_list, theta, gas_pressures, T)
        d_theta = stoich_matrix @ rates

        # Ensure coverage constraint: sum(theta) + theta_vacant = 1
        theta_new = theta + dt * d_theta

        # Project onto feasible region
        theta_new = np.clip(theta_new, 1e-15, 1.0)
        total = np.sum(theta_new)
        if total > 1.0:
            theta_new /= total

        residual = np.max(np.abs(d_theta))
        if residual < tol:
            break

        # Adaptive time stepping
        if residual < 1e-6:
            dt = min(dt * 1.5, 1e3)
        elif residual > 1e-2:
            dt = max(dt * 0.5, 1e-8)

        theta = theta_new

    return theta


class PFRReactor:
    """Plug Flow Reactor model coupled with microkinetic model."""

    def __init__(self, conditions: ReactorConditions,
                 surface_species: List[str],
                 gas_species: List[str],
                 stoich_gas: np.ndarray,
                 stoich_surface: np.ndarray,
                 rate_expressions: Callable):
        """
        Parameters
        ----------
        conditions : ReactorConditions
        surface_species : list of str
        gas_species : list of str
        stoich_gas : np.ndarray
            Gas-phase stoichiometric matrix [n_gas x n_reactions].
        stoich_surface : np.ndarray
            Surface stoichiometric matrix [n_surface x n_reactions].
        rate_expressions : callable
            Function(rate_constants, coverages, pressures, T) -> rates.
        """
        self.conditions = conditions
        self.surface_species = surface_species
        self.gas_species = gas_species
        self.stoich_gas = stoich_gas
        self.stoich_surface = stoich_surface
        self.rate_expressions = rate_expressions

    def _ode_system(self, W, y, rate_constants_list):
        """ODE system for PFR: dF/dW = r(θ, P)"""
        n_gas = len(self.gas_species)
        n_surf = len(self.surface_species)

        F = y[:n_gas]  # Molar flow rates
        theta = y[n_gas:n_gas + n_surf]  # Surface coverages

        # Partial pressures from molar flows
        F_total = np.sum(np.maximum(F, 1e-30))
        pressures = {sp: max(F[i] / F_total * self.conditions.P_total, 1e-30)
                     for i, sp in enumerate(self.gas_species)}

        # Calculate rates
        theta_clipped = np.clip(theta, 1e-15, 1.0)
        rates = self.rate_expressions(rate_constants_list, theta_clipped,
                                      pressures, self.conditions.T)

        # Gas phase: dF_i/dW = Σ ν_ij * r_j * site_density
        dF_dW = self.stoich_gas @ rates * self.conditions.site_density

        # Surface: pseudo-steady-state (dθ/dt ≈ 0 at each position)
        d_theta = self.stoich_surface @ rates
        tau_surface = 1e-6  # Fast surface relaxation
        d_theta_dW = d_theta * tau_surface

        return np.concatenate([dF_dW, d_theta_dW])

    def solve(self, rate_constants_list: list,
              W_span: Tuple[float, float] = None,
              n_points: int = 200) -> ReactorResult:
        """
        Solve PFR equations.

        Parameters
        ----------
        rate_constants_list : list
            Rate constants for all elementary steps.
        W_span : tuple
            (W_start, W_end) catalyst weight range [kg].
        n_points : int
            Number of output points.

        Returns
        -------
        ReactorResult
        """
        if W_span is None:
            W_span = (0, self.conditions.catalyst_mass)

        n_gas = len(self.gas_species)
        n_surf = len(self.surface_species)

        # Initial conditions
        F0 = np.array([self.conditions.feed_composition.get(sp, 0.0) *
                        self.conditions.F_total
                       for sp in self.gas_species])
        theta0 = np.ones(n_surf) * 0.01
        y0 = np.concatenate([F0, theta0])

        W_eval = np.linspace(W_span[0], W_span[1], n_points)

        sol = solve_ivp(
            lambda W, y: self._ode_system(W, y, rate_constants_list),
            W_span, y0, t_eval=W_eval,
            method='BDF', rtol=1e-8, atol=1e-12,
            max_step=(W_span[1] - W_span[0]) / 20
        )

        # Extract results
        gas_comp = {}
        for i, sp in enumerate(self.gas_species):
            F_total = np.sum(np.maximum(sol.y[:n_gas], 1e-30), axis=0)
            gas_comp[sp] = np.maximum(sol.y[i], 0) / F_total

        surf_cov = {}
        for i, sp in enumerate(self.surface_species):
            surf_cov[sp] = np.clip(sol.y[n_gas + i], 0, 1)

        # Conversion of first reactant
        F_in = F0[0]
        F_out = max(sol.y[0, -1], 0)
        conversion = (F_in - F_out) / F_in if F_in > 0 else 0

        # Selectivities (product moles produced / reactant moles consumed)
        selectivities = {}
        reactant_consumed = F_in - F_out
        for i, sp in enumerate(self.gas_species):
            if i == 0:
                continue
            F_prod = max(sol.y[i, -1] - F0[i], 0)
            selectivities[sp] = F_prod / reactant_consumed if reactant_consumed > 0 else 0

        # TOF
        tof = reactant_consumed / (self.conditions.catalyst_mass *
                                    self.conditions.site_density) if self.conditions.site_density > 0 else 0

        sty = reactant_consumed / self.conditions.catalyst_mass if self.conditions.catalyst_mass > 0 else 0

        return ReactorResult(
            reactor_type="PFR",
            positions=sol.t,
            gas_compositions=gas_comp,
            surface_coverages=surf_cov,
            conversion=conversion,
            selectivities=selectivities,
            turnover_frequency=tof,
            space_time_yield=sty,
            conditions=self.conditions
        )


class CSTRReactor:
    """Continuous Stirred Tank Reactor model."""

    def __init__(self, conditions: ReactorConditions,
                 surface_species: List[str],
                 gas_species: List[str],
                 stoich_gas: np.ndarray,
                 stoich_surface: np.ndarray,
                 rate_expressions: Callable):
        self.conditions = conditions
        self.surface_species = surface_species
        self.gas_species = gas_species
        self.stoich_gas = stoich_gas
        self.stoich_surface = stoich_surface
        self.rate_expressions = rate_expressions

    def solve(self, rate_constants_list: list,
              max_iter: int = 5000, tol: float = 1e-10) -> ReactorResult:
        """
        Solve CSTR equations at steady state.

        F_i,in - F_i,out + r_i * W_cat * n_sites = 0
        """
        n_gas = len(self.gas_species)
        n_surf = len(self.surface_species)

        F_in = np.array([self.conditions.feed_composition.get(sp, 0.0) *
                          self.conditions.F_total
                         for sp in self.gas_species])

        # Initial guess
        F_out = F_in.copy() * 0.9
        theta = np.ones(n_surf) * 0.01

        W = self.conditions.catalyst_mass
        n_s = self.conditions.site_density

        for iteration in range(max_iter):
            F_total = np.sum(np.maximum(F_out, 1e-30))
            pressures = {sp: max(F_out[i] / F_total * self.conditions.P_total, 1e-30)
                         for i, sp in enumerate(self.gas_species)}

            theta_clipped = np.clip(theta, 1e-15, 1.0)
            rates = self.rate_expressions(rate_constants_list, theta_clipped,
                                          pressures, self.conditions.T)

            # Gas phase residual
            R_gas = self.stoich_gas @ rates * n_s
            F_out_new = F_in + R_gas * W
            F_out_new = np.maximum(F_out_new, 1e-30)

            # Surface residual
            d_theta = self.stoich_surface @ rates
            theta_new = theta + 0.01 * d_theta
            theta_new = np.clip(theta_new, 1e-15, 1.0)
            total = np.sum(theta_new)
            if total > 1.0:
                theta_new /= total

            # Check convergence
            res_gas = np.max(np.abs(F_out_new - F_out) / np.maximum(F_in, 1e-30))
            res_surf = np.max(np.abs(d_theta))

            if res_gas < tol and res_surf < tol:
                break

            # Damped update
            F_out = 0.3 * F_out_new + 0.7 * F_out
            theta = 0.3 * theta_new + 0.7 * theta

        # Results
        gas_comp = {}
        F_total_out = np.sum(np.maximum(F_out, 1e-30))
        for i, sp in enumerate(self.gas_species):
            gas_comp[sp] = np.array([F_in[i] / np.sum(F_in),
                                     max(F_out[i], 0) / F_total_out])

        surf_cov = {sp: np.array([0.01, theta[i]])
                    for i, sp in enumerate(self.surface_species)}

        conversion = (F_in[0] - F_out[0]) / F_in[0] if F_in[0] > 0 else 0
        reactant_consumed = F_in[0] - F_out[0]

        selectivities = {}
        for i, sp in enumerate(self.gas_species):
            if i == 0:
                continue
            F_prod = max(F_out[i] - F_in.get(i, 0) if isinstance(F_in, dict) else F_out[i] - F_in[i], 0)
            selectivities[sp] = F_prod / reactant_consumed if reactant_consumed > 0 else 0

        tof = reactant_consumed / (W * n_s) if n_s > 0 else 0
        sty = reactant_consumed / W if W > 0 else 0

        return ReactorResult(
            reactor_type="CSTR",
            positions=np.array([0, 1]),
            gas_compositions=gas_comp,
            surface_coverages=surf_cov,
            conversion=conversion,
            selectivities=selectivities,
            turnover_frequency=tof,
            space_time_yield=sty,
            conditions=self.conditions
        )
