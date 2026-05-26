"""
Epidemic Modeling Framework: SIR, SEIR, Age-Structured SEIR, and Agent-Based Models
with Bayesian parameter estimation and model selection.
"""

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')


# ============================================================
# 1. Compartmental Models (ODE-based)
# ============================================================

def sir_ode(y, t, beta, gamma, N):
    S, I, R = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I
    return [dSdt, dIdt, dRdt]


def seir_ode(y, t, beta, sigma, gamma, N):
    S, E, I, R = y
    dSdt = -beta * S * I / N
    dEdt = beta * S * I / N - sigma * E
    dIdt = sigma * E - gamma * I
    dRdt = gamma * I
    return [dSdt, dEdt, dIdt, dRdt]


def seir_age_structured_ode(y, t, beta_matrix, sigma, gamma, N_groups):
    """Age-structured SEIR with contact matrix."""
    n = len(N_groups)
    S = y[0:n]
    E = y[n:2*n]
    I = y[2*n:3*n]
    R = y[3*n:4*n]

    # Force of infection for each age group
    lambda_i = np.zeros(n)
    for i in range(n):
        for j in range(n):
            lambda_i[i] += beta_matrix[i, j] * I[j] / N_groups[j]

    dSdt = -lambda_i * S
    dEdt = lambda_i * S - sigma * E
    dIdt = sigma * E - gamma * I
    dRdt = gamma * I
    return np.concatenate([dSdt, dEdt, dIdt, dRdt])


def seir_spatial_ode(y, t, beta_local, beta_travel, sigma, gamma, N_regions, mobility_matrix):
    """Spatially heterogeneous SEIR with mobility."""
    n = len(N_regions)
    S = y[0:n]
    E = y[n:2*n]
    I = y[2*n:3*n]
    R = y[3*n:4*n]

    dSdt = np.zeros(n)
    dEdt = np.zeros(n)
    dIdt = np.zeros(n)
    dRdt = np.zeros(n)

    for i in range(n):
        # Local transmission
        local_foi = beta_local * I[i] / N_regions[i]
        # Travel-related transmission
        travel_foi = 0
        for j in range(n):
            if i != j:
                travel_foi += beta_travel * mobility_matrix[j, i] * I[j] / N_regions[j]

        total_foi = local_foi + travel_foi
        dSdt[i] = -total_foi * S[i]
        dEdt[i] = total_foi * S[i] - sigma * E[i]
        dIdt[i] = sigma * E[i] - gamma * I[i]
        dRdt[i] = gamma * I[i]

    return np.concatenate([dSdt, dEdt, dIdt, dRdt])


def seir_vaccination_ode(y, t, beta, sigma, gamma, N, vaccination_rate, vaccine_efficacy):
    """SEIR with vaccination."""
    S, E, I, R, V = y
    effective_vax = vaccination_rate * vaccine_efficacy
    dSdt = -beta * S * I / N - effective_vax * S
    dEdt = beta * S * I / N + beta * (1 - vaccine_efficacy) * V * I / N - sigma * E
    dIdt = sigma * E - gamma * I
    dRdt = gamma * I
    dVdt = effective_vax * S - beta * (1 - vaccine_efficacy) * V * I / N
    return [dSdt, dEdt, dIdt, dRdt, dVdt]


# ============================================================
# 2. Agent-Based Model
# ============================================================

class Agent:
    def __init__(self, agent_id, age_group=0, x=0.0, y=0.0):
        self.id = agent_id
        self.state = 'S'
        self.age_group = age_group
        self.x = x
        self.y = y
        self.days_infected = 0
        self.days_exposed = 0
        self.vaccinated = False


class ABMEpidemic:
    def __init__(self, n_agents, beta, sigma, gamma, contact_radius=0.05,
                 n_age_groups=3, seed=42):
        self.rng = np.random.RandomState(seed)
        self.beta = beta
        self.sigma = sigma
        self.gamma = gamma
        self.contact_radius = contact_radius

        self.agents = []
        for i in range(n_agents):
            age = self.rng.choice(n_age_groups)
            x = self.rng.uniform(0, 1)
            y_pos = self.rng.uniform(0, 1)
            self.agents.append(Agent(i, age, x, y_pos))

    def seed_infection(self, n_initial=5):
        indices = self.rng.choice(len(self.agents), n_initial, replace=False)
        for idx in indices:
            self.agents[idx].state = 'I'

    def step(self):
        # Find infected agents
        infected = [a for a in self.agents if a.state == 'I']
        susceptible = [a for a in self.agents if a.state == 'S']

        # Transmission
        for s_agent in susceptible:
            for i_agent in infected:
                dist = np.sqrt((s_agent.x - i_agent.x)**2 + (s_agent.y - i_agent.y)**2)
                if dist < self.contact_radius:
                    if self.rng.random() < self.beta:
                        s_agent.state = 'E'
                        break

        # State transitions
        for agent in self.agents:
            if agent.state == 'E':
                agent.days_exposed += 1
                if self.rng.random() < self.sigma:
                    agent.state = 'I'
                    agent.days_exposed = 0
            elif agent.state == 'I':
                agent.days_infected += 1
                if self.rng.random() < self.gamma:
                    agent.state = 'R'
                    agent.days_infected = 0

        # Random movement
        for agent in self.agents:
            agent.x = np.clip(agent.x + self.rng.normal(0, 0.01), 0, 1)
            agent.y = np.clip(agent.y + self.rng.normal(0, 0.01), 0, 1)

    def get_counts(self):
        counts = {'S': 0, 'E': 0, 'I': 0, 'R': 0}
        for a in self.agents:
            counts[a.state] += 1
        return counts

    def run(self, n_steps):
        history = {'S': [], 'E': [], 'I': [], 'R': []}
        for _ in range(n_steps):
            self.step()
            c = self.get_counts()
            for k in history:
                history[k].append(c[k])
        return history


# ============================================================
# 3. Simulate observed data (synthetic COVID-19-like)
# ============================================================

def generate_synthetic_data(model_type='seir', N=1e6, n_days=120, noise_scale=0.05, seed=42):
    """Generate synthetic epidemic data with observation noise."""
    rng = np.random.RandomState(seed)
    t = np.arange(n_days)

    if model_type == 'sir':
        beta, gamma = 0.3, 0.1
        y0 = [N - 100, 100, 0]
        sol = odeint(sir_ode, y0, t, args=(beta, gamma, N))
        I_true = sol[:, 1]
    elif model_type == 'seir':
        beta, sigma, gamma = 0.35, 0.2, 0.1
        y0 = [N - 100, 50, 50, 0]
        sol = odeint(seir_ode, y0, t, args=(beta, sigma, gamma, N))
        I_true = sol[:, 2]
    else:
        raise ValueError(f"Unknown model type: {model_type}")

    # Daily new cases (approximate)
    new_cases = np.diff(np.concatenate([[0], np.cumsum(I_true * gamma)]))
    new_cases = np.maximum(new_cases, 0)
    observed = new_cases + rng.normal(0, noise_scale * np.maximum(new_cases, 1), len(new_cases))
    observed = np.maximum(observed, 0).astype(int)

    return t, I_true, observed, new_cases


def generate_covid_wave_data(wave='6th', seed=42):
    """Generate synthetic data mimicking COVID-19 6th/7th wave in Japan."""
    rng = np.random.RandomState(seed)

    if wave == '6th':
        # 6th wave (Omicron BA.1, Jan-Mar 2022): R0~8, shorter incubation
        N = 1.26e8
        beta, sigma, gamma = 0.8, 0.33, 0.14
        n_days = 90
        I0 = 5000
    else:
        # 7th wave (Omicron BA.5, Jul-Sep 2022): R0~10, high transmissibility
        N = 1.26e8
        beta, sigma, gamma = 1.0, 0.33, 0.14
        n_days = 100
        I0 = 3000

    t = np.arange(n_days)
    y0 = [N - I0 * 2, I0, I0, 0]
    sol = odeint(seir_ode, y0, t, args=(beta, sigma, gamma, N))

    I_true = sol[:, 2]
    new_cases_true = np.diff(np.concatenate([[0], np.cumsum(sol[:, 1] * sigma)]))
    new_cases_true = np.maximum(new_cases_true, 0)

    reporting_rate = 0.3
    observed = rng.poisson(new_cases_true * reporting_rate)

    return t, I_true, observed, new_cases_true
