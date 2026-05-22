"""
Compartmental Epidemiological Models Framework
================================================
SIR, SEIR, Age-Structured SEIR, Metapopulation SEIR, Vaccination SEIR
with intervention modeling and reproduction number calculation.

Scientific reference:
- Keeling & Rohani (2008) "Modeling Infectious Diseases in Humans and Animals"
- Diekmann, Heesterbeek & Britton (2013) "Mathematical Tools for Understanding
  Infectious Disease Dynamics"
"""

import numpy as np
from scipy.integrate import solve_ivp
from typing import Dict, List, Optional, Tuple, Callable
from dataclasses import dataclass, field


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class SimulationResult:
    """Container for simulation output."""
    t: np.ndarray
    y: Dict[str, np.ndarray]
    R0: float
    Reff_series: Optional[np.ndarray] = None
    metadata: Dict = field(default_factory=dict)


@dataclass
class InterventionSchedule:
    """Time-varying intervention specification."""
    name: str
    start_day: float
    end_day: float
    contact_reduction: float = 0.0        # fraction reduction in beta
    vaccination_rate: float = 0.0          # daily per-capita vaccination rate
    vaccine_efficacy: float = 0.0          # VE against infection


# ---------------------------------------------------------------------------
# 1. Basic SIR Model
# ---------------------------------------------------------------------------

class SIRModel:
    """
    Standard SIR (Susceptible-Infected-Recovered) model.

    dS/dt = -beta * S * I / N
    dI/dt =  beta * S * I / N - gamma * I
    dR/dt =  gamma * I

    Parameters
    ----------
    beta  : transmission rate (contacts × probability per contact per day)
    gamma : recovery rate (1/gamma = infectious period in days)
    N     : total population
    """

    def __init__(self, beta: float, gamma: float, N: float):
        self.beta = beta
        self.gamma = gamma
        self.N = N

    def compute_R0(self) -> float:
        return self.beta / self.gamma

    def compute_Reff(self, S: float) -> float:
        return self.compute_R0() * (S / self.N)

    def _ode(self, t: float, y: np.ndarray,
             beta_t: Callable[[float], float]) -> np.ndarray:
        S, I, R = y
        b = beta_t(t)
        dS = -b * S * I / self.N
        dI = b * S * I / self.N - self.gamma * I
        dR = self.gamma * I
        return [dS, dI, dR]

    def simulate(self, S0: float, I0: float, R0: float,
                 t_span: Tuple[float, float], t_eval: Optional[np.ndarray] = None,
                 interventions: Optional[List[InterventionSchedule]] = None
                 ) -> SimulationResult:
        if t_eval is None:
            t_eval = np.arange(t_span[0], t_span[1] + 1)

        beta_t = self._build_beta_t(interventions)
        sol = solve_ivp(self._ode, t_span, [S0, I0, R0],
                        t_eval=t_eval, args=(beta_t,), method="RK45",
                        max_step=1.0)

        Reff = np.array([self.compute_Reff(s) for s in sol.y[0]])
        return SimulationResult(
            t=sol.t,
            y={"S": sol.y[0], "I": sol.y[1], "R": sol.y[2]},
            R0=self.compute_R0(),
            Reff_series=Reff,
            metadata={"model": "SIR", "beta": self.beta, "gamma": self.gamma},
        )

    def _build_beta_t(self, interventions: Optional[List[InterventionSchedule]]
                      ) -> Callable[[float], float]:
        if not interventions:
            return lambda t: self.beta

        def beta_t(t: float) -> float:
            reduction = 0.0
            for iv in interventions:
                if iv.start_day <= t <= iv.end_day:
                    reduction = max(reduction, iv.contact_reduction)
            return self.beta * (1.0 - reduction)
        return beta_t


# ---------------------------------------------------------------------------
# 2. SEIR Model
# ---------------------------------------------------------------------------

class SEIRModel:
    """
    SEIR model adding an Exposed (latent) compartment.

    dS/dt = -beta * S * I / N
    dE/dt =  beta * S * I / N - sigma * E
    dI/dt =  sigma * E - gamma * I
    dR/dt =  gamma * I

    Parameters
    ----------
    beta  : transmission rate
    sigma : incubation rate (1/sigma = mean latent period)
    gamma : recovery rate (1/gamma = mean infectious period)
    N     : total population
    """

    def __init__(self, beta: float, sigma: float, gamma: float, N: float):
        self.beta = beta
        self.sigma = sigma
        self.gamma = gamma
        self.N = N

    def compute_R0(self) -> float:
        return self.beta / self.gamma

    def compute_Reff(self, S: float) -> float:
        return self.compute_R0() * (S / self.N)

    def _ode(self, t: float, y: np.ndarray,
             beta_t: Callable[[float], float]) -> np.ndarray:
        S, E, I, R = y
        b = beta_t(t)
        force = b * S * I / self.N
        dS = -force
        dE = force - self.sigma * E
        dI = self.sigma * E - self.gamma * I
        dR = self.gamma * I
        return [dS, dE, dI, dR]

    def simulate(self, S0, E0, I0, R0, t_span, t_eval=None,
                 interventions=None) -> SimulationResult:
        if t_eval is None:
            t_eval = np.arange(t_span[0], t_span[1] + 1)

        beta_t = self._build_beta_t(interventions)
        sol = solve_ivp(self._ode, t_span, [S0, E0, I0, R0],
                        t_eval=t_eval, args=(beta_t,), method="RK45",
                        max_step=1.0)

        Reff = np.array([self.compute_Reff(s) for s in sol.y[0]])
        return SimulationResult(
            t=sol.t,
            y={"S": sol.y[0], "E": sol.y[1], "I": sol.y[2], "R": sol.y[3]},
            R0=self.compute_R0(),
            Reff_series=Reff,
            metadata={"model": "SEIR", "beta": self.beta,
                       "sigma": self.sigma, "gamma": self.gamma},
        )

    def _build_beta_t(self, interventions):
        if not interventions:
            return lambda t: self.beta

        def beta_t(t):
            reduction = 0.0
            for iv in interventions:
                if iv.start_day <= t <= iv.end_day:
                    reduction = max(reduction, iv.contact_reduction)
            return self.beta * (1.0 - reduction)
        return beta_t


# ---------------------------------------------------------------------------
# 3. Age-Structured SEIR
# ---------------------------------------------------------------------------

class AgeStructuredSEIR:
    """
    Age-structured SEIR with contact matrix.

    Age groups: 0-19, 20-39, 40-64, 65+  (configurable)

    The force of infection for age group i:
        lambda_i = sum_j  beta * C_ij * I_j / N_j

    R0 is computed from the next-generation matrix (NGM):
        K_ij = beta * C_ij * (S_i / N_i) / gamma
        R0   = spectral_radius(K)

    Parameters
    ----------
    beta          : baseline transmission probability per contact
    sigma         : incubation rate
    gamma         : recovery rate
    N_age         : array of population sizes per age group
    contact_matrix: square matrix C_ij (avg daily contacts from i to j)
    susceptibility: relative susceptibility by age group (default all 1)
    """

    AGE_LABELS = ["0-19", "20-39", "40-64", "65+"]

    # Typical contact matrix (Mossong et al. 2008 POLYMOD, Japan-adjusted)
    DEFAULT_CONTACT_MATRIX = np.array([
        [8.0, 3.0, 2.5, 1.0],
        [3.0, 7.0, 4.0, 1.5],
        [2.5, 4.0, 6.0, 2.0],
        [1.0, 1.5, 2.0, 3.0],
    ])

    def __init__(self, beta: float, sigma: float, gamma: float,
                 N_age: np.ndarray,
                 contact_matrix: Optional[np.ndarray] = None,
                 susceptibility: Optional[np.ndarray] = None):
        self.beta = beta
        self.sigma = sigma
        self.gamma = gamma
        self.n_groups = len(N_age)
        self.N_age = np.asarray(N_age, dtype=float)
        self.N_total = self.N_age.sum()
        self.C = contact_matrix if contact_matrix is not None else self.DEFAULT_CONTACT_MATRIX
        self.susceptibility = (susceptibility if susceptibility is not None
                               else np.ones(self.n_groups))

    def compute_R0(self) -> float:
        """R0 from spectral radius of the next-generation matrix."""
        K = np.zeros((self.n_groups, self.n_groups))
        for i in range(self.n_groups):
            for j in range(self.n_groups):
                K[i, j] = (self.beta * self.susceptibility[i]
                           * self.C[i, j] / self.gamma)
        return float(np.max(np.abs(np.linalg.eigvals(K))))

    def compute_Reff(self, S_age: np.ndarray) -> float:
        K = np.zeros((self.n_groups, self.n_groups))
        for i in range(self.n_groups):
            for j in range(self.n_groups):
                K[i, j] = (self.beta * self.susceptibility[i]
                           * self.C[i, j] * (S_age[i] / self.N_age[i])
                           / self.gamma)
        return float(np.max(np.abs(np.linalg.eigvals(K))))

    def _ode(self, t, y, beta_t):
        n = self.n_groups
        S = y[0:n]
        E = y[n:2*n]
        I = y[2*n:3*n]
        R = y[3*n:4*n]

        b = beta_t(t)
        force = np.zeros(n)
        for i in range(n):
            for j in range(n):
                force[i] += b * self.susceptibility[i] * self.C[i, j] * I[j] / self.N_age[j]

        dS = -force * S
        dE = force * S - self.sigma * E
        dI = self.sigma * E - self.gamma * I
        dR = self.gamma * I
        return np.concatenate([dS, dE, dI, dR])

    def simulate(self, S0, E0, I0, R0, t_span, t_eval=None,
                 interventions=None) -> SimulationResult:
        S0, E0, I0, R0 = [np.asarray(x, dtype=float) for x in [S0, E0, I0, R0]]
        if t_eval is None:
            t_eval = np.arange(t_span[0], t_span[1] + 1)

        beta_t = self._build_beta_t(interventions)
        y0 = np.concatenate([S0, E0, I0, R0])
        sol = solve_ivp(self._ode, t_span, y0, t_eval=t_eval,
                        args=(beta_t,), method="RK45", max_step=1.0)

        n = self.n_groups
        result_y = {}
        for k, label in enumerate(self.AGE_LABELS[:n]):
            result_y[f"S_{label}"] = sol.y[k]
            result_y[f"E_{label}"] = sol.y[n + k]
            result_y[f"I_{label}"] = sol.y[2*n + k]
            result_y[f"R_{label}"] = sol.y[3*n + k]
        result_y["I_total"] = sum(sol.y[2*n + k] for k in range(n))
        result_y["S_total"] = sum(sol.y[k] for k in range(n))

        Reff = np.array([
            self.compute_Reff(sol.y[0:n, idx]) for idx in range(len(sol.t))
        ])

        return SimulationResult(
            t=sol.t, y=result_y,
            R0=self.compute_R0(), Reff_series=Reff,
            metadata={"model": "AgeStructuredSEIR",
                       "age_groups": self.AGE_LABELS[:n]},
        )

    def _build_beta_t(self, interventions):
        if not interventions:
            return lambda t: self.beta
        def beta_t(t):
            reduction = 0.0
            for iv in interventions:
                if iv.start_day <= t <= iv.end_day:
                    reduction = max(reduction, iv.contact_reduction)
            return self.beta * (1.0 - reduction)
        return beta_t


# ---------------------------------------------------------------------------
# 4. Metapopulation SEIR (spatial heterogeneity)
# ---------------------------------------------------------------------------

class MetapopulationSEIR:
    """
    Multi-patch SEIR with mobility coupling.

    Each patch i has its own SEIR dynamics. Patches are coupled via
    a mobility matrix M where M_ij = fraction of i's population
    commuting to j daily.

    Effective force of infection in patch i mixes local and commuter
    contributions.

    Parameters
    ----------
    n_patches       : number of spatial patches
    beta            : array of per-patch transmission rates
    sigma, gamma    : scalar (uniform across patches for simplicity)
    N_patch         : array of population per patch
    mobility_matrix : M_ij, row-stochastic (rows sum to 1)
    """

    def __init__(self, beta: np.ndarray, sigma: float, gamma: float,
                 N_patch: np.ndarray, mobility_matrix: np.ndarray):
        self.n = len(N_patch)
        self.beta = np.asarray(beta, dtype=float)
        self.sigma = sigma
        self.gamma = gamma
        self.N_patch = np.asarray(N_patch, dtype=float)
        self.M = np.asarray(mobility_matrix, dtype=float)

    def compute_R0(self) -> float:
        """Spatial R0 from the next-generation matrix across patches."""
        K = np.zeros((self.n, self.n))
        for i in range(self.n):
            for j in range(self.n):
                K[i, j] = self.beta[j] * self.M[i, j] / self.gamma
        return float(np.max(np.abs(np.linalg.eigvals(K))))

    def _ode(self, t, y, beta_t_func):
        n = self.n
        S = y[0:n]
        E = y[n:2*n]
        I = y[2*n:3*n]
        R = y[3*n:4*n]

        # Effective I* in each patch accounting for commuters
        I_eff = np.zeros(n)
        N_eff = np.zeros(n)
        for j in range(n):
            for i in range(n):
                I_eff[j] += self.M[i, j] * I[i]
                N_eff[j] += self.M[i, j] * self.N_patch[i]

        dS = np.zeros(n)
        dE = np.zeros(n)
        dI = np.zeros(n)
        dR = np.zeros(n)

        for i in range(n):
            force = 0.0
            for j in range(n):
                if N_eff[j] > 0:
                    force += self.M[i, j] * beta_t_func(t, j) * I_eff[j] / N_eff[j]
            dS[i] = -force * S[i]
            dE[i] = force * S[i] - self.sigma * E[i]
            dI[i] = self.sigma * E[i] - self.gamma * I[i]
            dR[i] = self.gamma * I[i]

        return np.concatenate([dS, dE, dI, dR])

    def simulate(self, S0, E0, I0, R0, t_span, t_eval=None,
                 interventions=None) -> SimulationResult:
        if t_eval is None:
            t_eval = np.arange(t_span[0], t_span[1] + 1)

        def beta_t_func(t, patch_idx):
            reduction = 0.0
            if interventions:
                for iv in interventions:
                    if iv.start_day <= t <= iv.end_day:
                        reduction = max(reduction, iv.contact_reduction)
            return self.beta[patch_idx] * (1.0 - reduction)

        y0 = np.concatenate([S0, E0, I0, R0])
        sol = solve_ivp(self._ode, t_span, y0, t_eval=t_eval,
                        args=(beta_t_func,), method="RK45", max_step=1.0)

        n = self.n
        result_y = {}
        for k in range(n):
            result_y[f"S_patch{k}"] = sol.y[k]
            result_y[f"I_patch{k}"] = sol.y[2*n + k]
        result_y["I_total"] = sum(sol.y[2*n + k] for k in range(n))

        return SimulationResult(
            t=sol.t, y=result_y, R0=self.compute_R0(),
            metadata={"model": "MetapopulationSEIR", "n_patches": n},
        )


# ---------------------------------------------------------------------------
# 5. Vaccination SEIR
# ---------------------------------------------------------------------------

class VaccinationSEIR:
    """
    SEIR with vaccination, waning immunity, and two-dose schedule.

    Compartments: S, E, I, R, V1 (one dose), V2 (two doses)

    Transitions:
        S  → E   : infection (beta, reduced by VE for V)
        S  → V1  : first dose vaccination
        V1 → V2  : second dose (after delay)
        V2 → S   : waning immunity
        E  → I   : incubation (sigma)
        I  → R   : recovery (gamma)

    Parameters
    ----------
    beta         : transmission rate
    sigma        : incubation rate
    gamma        : recovery rate
    N            : total population
    ve_dose1     : vaccine efficacy after dose 1
    ve_dose2     : vaccine efficacy after dose 2
    waning_rate  : rate of immunity loss (1/waning_rate = duration)
    dose2_delay  : days between dose 1 and dose 2
    """

    def __init__(self, beta, sigma, gamma, N,
                 ve_dose1=0.5, ve_dose2=0.9,
                 waning_rate=1/365, dose2_delay=28):
        self.beta = beta
        self.sigma = sigma
        self.gamma = gamma
        self.N = N
        self.ve1 = ve_dose1
        self.ve2 = ve_dose2
        self.waning = waning_rate
        self.dose2_rate = 1.0 / dose2_delay

    def compute_R0(self) -> float:
        return self.beta / self.gamma

    def compute_Reff(self, S, V1, V2) -> float:
        effective_susceptible = S + V1 * (1 - self.ve1) + V2 * (1 - self.ve2)
        return self.compute_R0() * (effective_susceptible / self.N)

    def _ode(self, t, y, beta_t, vax_rate_t):
        S, E, I, R, V1, V2 = y
        b = beta_t(t)
        v = vax_rate_t(t)

        force = b * I / self.N
        dS = -force * S - v * S + self.waning * V2
        dE = force * (S + V1 * (1 - self.ve1) + V2 * (1 - self.ve2)) - self.sigma * E
        dI = self.sigma * E - self.gamma * I
        dR = self.gamma * I
        dV1 = v * S - self.dose2_rate * V1 - force * V1 * (1 - self.ve1)
        dV2 = self.dose2_rate * V1 - self.waning * V2 - force * V2 * (1 - self.ve2)
        return [dS, dE, dI, dR, dV1, dV2]

    def simulate(self, S0, E0, I0, R0, V1_0, V2_0, t_span, t_eval=None,
                 interventions=None) -> SimulationResult:
        if t_eval is None:
            t_eval = np.arange(t_span[0], t_span[1] + 1)

        beta_t = self._build_beta_t(interventions)
        vax_rate_t = self._build_vax_t(interventions)

        sol = solve_ivp(self._ode, t_span, [S0, E0, I0, R0, V1_0, V2_0],
                        t_eval=t_eval, args=(beta_t, vax_rate_t),
                        method="RK45", max_step=1.0)

        Reff = np.array([
            self.compute_Reff(sol.y[0, i], sol.y[4, i], sol.y[5, i])
            for i in range(len(sol.t))
        ])

        return SimulationResult(
            t=sol.t,
            y={"S": sol.y[0], "E": sol.y[1], "I": sol.y[2],
               "R": sol.y[3], "V1": sol.y[4], "V2": sol.y[5]},
            R0=self.compute_R0(), Reff_series=Reff,
            metadata={"model": "VaccinationSEIR",
                       "ve_dose1": self.ve1, "ve_dose2": self.ve2},
        )

    def _build_beta_t(self, interventions):
        if not interventions:
            return lambda t: self.beta
        def beta_t(t):
            reduction = 0.0
            for iv in interventions:
                if iv.start_day <= t <= iv.end_day:
                    reduction = max(reduction, iv.contact_reduction)
            return self.beta * (1.0 - reduction)
        return beta_t

    def _build_vax_t(self, interventions):
        if not interventions:
            return lambda t: 0.0
        def vax_rate(t):
            rate = 0.0
            for iv in interventions:
                if iv.start_day <= t <= iv.end_day:
                    rate = max(rate, iv.vaccination_rate)
            return rate
        return vax_rate


# ---------------------------------------------------------------------------
# Helper / dispatcher functions
# ---------------------------------------------------------------------------

def simulate_scenario(model_type: str, params: dict,
                      interventions: Optional[List[InterventionSchedule]] = None,
                      t_span: Tuple[float, float] = (0, 180)
                      ) -> SimulationResult:
    """Dispatch to the appropriate model class."""
    N = params.get("N", 1_000_000)

    if model_type == "SIR":
        m = SIRModel(params["beta"], params["gamma"], N)
        return m.simulate(N - params.get("I0", 10), params.get("I0", 10), 0,
                          t_span, interventions=interventions)

    elif model_type == "SEIR":
        m = SEIRModel(params["beta"], params["sigma"], params["gamma"], N)
        return m.simulate(N - params.get("I0", 10), 0, params.get("I0", 10), 0,
                          t_span, interventions=interventions)

    elif model_type == "AgeStructuredSEIR":
        fracs = params.get("age_fractions", [0.15, 0.25, 0.35, 0.25])
        N_age = np.array(fracs) * N
        m = AgeStructuredSEIR(params["beta"], params["sigma"],
                              params["gamma"], N_age)
        I0_age = np.array(fracs) * params.get("I0", 10)
        S0_age = N_age - I0_age
        return m.simulate(S0_age, np.zeros(4), I0_age, np.zeros(4),
                          t_span, interventions=interventions)

    elif model_type == "VaccinationSEIR":
        m = VaccinationSEIR(
            params["beta"], params["sigma"], params["gamma"], N,
            ve_dose1=params.get("ve_dose1", 0.5),
            ve_dose2=params.get("ve_dose2", 0.9),
        )
        return m.simulate(N - params.get("I0", 10), 0, params.get("I0", 10),
                          0, 0, 0, t_span, interventions=interventions)

    else:
        raise ValueError(f"Unknown model_type: {model_type}")


def compute_R0(model_type: str, params: dict) -> float:
    """Quick R0 computation for any model type."""
    N = params.get("N", 1_000_000)
    if model_type in ("SIR",):
        return params["beta"] / params["gamma"]
    elif model_type in ("SEIR", "VaccinationSEIR"):
        return params["beta"] / params["gamma"]
    elif model_type == "AgeStructuredSEIR":
        fracs = params.get("age_fractions", [0.15, 0.25, 0.35, 0.25])
        N_age = np.array(fracs) * N
        m = AgeStructuredSEIR(params["beta"], params["sigma"],
                              params["gamma"], N_age)
        return m.compute_R0()
    else:
        return params["beta"] / params["gamma"]


# ---------------------------------------------------------------------------
# Demonstration
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("Compartmental Models Framework — Demonstration")
    print("=" * 60)

    # SEIR baseline
    params = {"beta": 0.5, "sigma": 1/3, "gamma": 1/7, "N": 1_000_000, "I0": 100}
    model = SEIRModel(params["beta"], params["sigma"], params["gamma"], params["N"])
    print(f"\nSEIR R0 = {model.compute_R0():.2f}")

    res = model.simulate(params["N"] - 100, 0, 100, 0, (0, 300))
    peak_I = np.max(res.y["I"])
    peak_day = res.t[np.argmax(res.y["I"])]
    print(f"Peak infected: {peak_I:,.0f} on day {peak_day:.0f}")

    # With 50% contact reduction (lockdown)
    lockdown = InterventionSchedule("lockdown", start_day=30, end_day=90,
                                     contact_reduction=0.5)
    res_iv = model.simulate(params["N"] - 100, 0, 100, 0, (0, 300),
                            interventions=[lockdown])
    peak_I_iv = np.max(res_iv.y["I"])
    peak_day_iv = res_iv.t[np.argmax(res_iv.y["I"])]
    print(f"With lockdown — Peak: {peak_I_iv:,.0f} on day {peak_day_iv:.0f}")
    print(f"Peak reduction: {(1 - peak_I_iv/peak_I)*100:.1f}%")

    # Age-structured
    N_age = np.array([15e6, 25e6, 35e6, 25e6])  # ~100M
    age_model = AgeStructuredSEIR(0.03, 1/3, 1/7, N_age)
    print(f"\nAge-Structured R0 = {age_model.compute_R0():.2f}")

    print("\nDone.")
