"""
Bayesian Parameter Estimation and Model Selection for Epidemic Models.
Uses PyMC for MCMC inference, WAIC and LOO-CV for model comparison.
"""

import numpy as np
import pymc as pm
import arviz as az
from scipy.integrate import odeint
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')


def seir_deterministic(params, t, N):
    beta, sigma, gamma = params
    y0 = [N - 100, 50, 50, 0]
    def ode(y, t_):
        S, E, I, R = y
        return [-beta*S*I/N, beta*S*I/N - sigma*E, sigma*E - gamma*I, gamma*I]
    sol = odeint(ode, y0, t)
    return sol


def sir_deterministic(params, t, N):
    beta, gamma = params
    y0 = [N - 100, 100, 0]
    def ode(y, t_):
        S, I, R = y
        return [-beta*S*I/N, beta*S*I/N - gamma*I, gamma*I]
    sol = odeint(ode, y0, t)
    return sol


def fit_sir_pymc(observed, t, N, n_samples=1000, n_tune=500, seed=42):
    """Bayesian SIR model fitting with PyMC."""
    with pm.Model() as sir_model:
        beta = pm.Lognormal('beta', mu=np.log(0.3), sigma=0.5)
        gamma = pm.Lognormal('gamma', mu=np.log(0.1), sigma=0.5)

        # Solve ODE deterministically
        sol = sir_deterministic([beta.eval(), gamma.eval()], t, N)
        I_pred_approx = np.maximum(sol[:, 1], 1)

        # Use approximate likelihood
        sigma_obs = pm.HalfNormal('sigma_obs', sigma=1000)

        # Simplified: use log-normal likelihood on I curve
        I_pred = pm.Deterministic('I_pred',
            pm.math.abs(beta * (N - 100) / N * 100 * np.exp((beta - gamma) * t[:len(observed)])))

        likelihood = pm.Normal('obs', mu=I_pred, sigma=sigma_obs,
                              observed=observed[:len(t)])

        trace = pm.sample(n_samples, tune=n_tune, random_seed=seed,
                         cores=1, progressbar=False, return_inferencedata=True)
    return trace, sir_model


def fit_seir_pymc(observed, t, N, n_samples=1000, n_tune=500, seed=42):
    """Bayesian SEIR model fitting with PyMC."""
    with pm.Model() as seir_model:
        beta = pm.Lognormal('beta', mu=np.log(0.35), sigma=0.5)
        sigma_e = pm.Lognormal('sigma_e', mu=np.log(0.2), sigma=0.5)
        gamma = pm.Lognormal('gamma', mu=np.log(0.1), sigma=0.5)

        sigma_obs = pm.HalfNormal('sigma_obs', sigma=1000)

        I_pred = pm.Deterministic('I_pred',
            pm.math.abs(beta * sigma_e / (sigma_e + gamma) *
                       (N - 100) / N * 50 * np.exp((beta * sigma_e / (sigma_e + gamma) - gamma) * t[:len(observed)])))

        likelihood = pm.Normal('obs', mu=I_pred, sigma=sigma_obs,
                              observed=observed[:len(t)])

        trace = pm.sample(n_samples, tune=n_tune, random_seed=seed,
                         cores=1, progressbar=False, return_inferencedata=True)
    return trace, seir_model


def compute_model_comparison(traces, model_names):
    """Compute WAIC and LOO-CV for model comparison."""
    results = {}
    for name, trace in zip(model_names, traces):
        try:
            waic = az.waic(trace)
            results[name] = {
                'waic': waic.elpd_waic,
                'waic_se': waic.se,
                'p_waic': waic.p_waic
            }
        except Exception:
            results[name] = {'waic': np.nan, 'waic_se': np.nan, 'p_waic': np.nan}
        try:
            loo = az.loo(trace)
            results[name]['loo'] = loo.elpd_loo
            results[name]['loo_se'] = loo.se
            results[name]['p_loo'] = loo.p_loo
        except Exception:
            results[name]['loo'] = np.nan
            results[name]['loo_se'] = np.nan
            results[name]['p_loo'] = np.nan
    return results


# ============================================================
# Simplified Parameter Estimation (MLE + Bootstrap for speed)
# ============================================================

def estimate_sir_mle(observed, t, N):
    """Maximum likelihood estimation for SIR."""
    def neg_log_lik(params):
        beta, gamma = np.exp(params)
        y0 = [N - 100, 100, 0]
        try:
            sol = odeint(lambda y, t_: [-beta*y[0]*y[1]/N,
                                         beta*y[0]*y[1]/N - gamma*y[1],
                                         gamma*y[1]], y0, t)
            I_pred = sol[:, 1]
            new_cases = np.diff(np.concatenate([[0], np.cumsum(I_pred * gamma)]))
            new_cases = np.maximum(new_cases, 1)
            obs_trimmed = observed[:len(new_cases)]
            return np.sum((obs_trimmed - new_cases[:len(obs_trimmed)])**2 / new_cases[:len(obs_trimmed)])
        except Exception:
            return 1e10

    result = minimize(neg_log_lik, [np.log(0.3), np.log(0.1)], method='Nelder-Mead')
    beta, gamma = np.exp(result.x)
    return {'beta': beta, 'gamma': gamma, 'R0': beta/gamma, 'neg_log_lik': result.fun}


def estimate_seir_mle(observed, t, N):
    """Maximum likelihood estimation for SEIR."""
    def neg_log_lik(params):
        beta, sigma, gamma = np.exp(params)
        y0 = [N - 100, 50, 50, 0]
        try:
            sol = odeint(lambda y, t_: [-beta*y[0]*y[2]/N,
                                         beta*y[0]*y[2]/N - sigma*y[1],
                                         sigma*y[1] - gamma*y[2],
                                         gamma*y[2]], y0, t)
            I_pred = sol[:, 2]
            new_cases = np.diff(np.concatenate([[0], np.cumsum(I_pred * gamma)]))
            new_cases = np.maximum(new_cases, 1)
            obs_trimmed = observed[:len(new_cases)]
            return np.sum((obs_trimmed - new_cases[:len(obs_trimmed)])**2 / new_cases[:len(obs_trimmed)])
        except Exception:
            return 1e10

    result = minimize(neg_log_lik, [np.log(0.35), np.log(0.2), np.log(0.1)], method='Nelder-Mead')
    beta, sigma, gamma = np.exp(result.x)
    return {'beta': beta, 'sigma': sigma, 'gamma': gamma, 'R0': beta/gamma, 'neg_log_lik': result.fun}


def compute_bic(neg_log_lik, n_params, n_data):
    """Bayesian Information Criterion."""
    return n_params * np.log(n_data) + 2 * neg_log_lik


def compute_aic(neg_log_lik, n_params):
    """Akaike Information Criterion."""
    return 2 * n_params + 2 * neg_log_lik


def model_selection_criteria(sir_result, seir_result, n_data):
    """Compare SIR vs SEIR using AIC/BIC."""
    sir_aic = compute_aic(sir_result['neg_log_lik'], 2)
    sir_bic = compute_bic(sir_result['neg_log_lik'], 2, n_data)
    seir_aic = compute_aic(seir_result['neg_log_lik'], 3)
    seir_bic = compute_bic(seir_result['neg_log_lik'], 3, n_data)

    return {
        'SIR': {'AIC': sir_aic, 'BIC': sir_bic, 'R0': sir_result['R0']},
        'SEIR': {'AIC': seir_aic, 'BIC': seir_bic, 'R0': seir_result['R0']},
        'preferred_AIC': 'SIR' if sir_aic < seir_aic else 'SEIR',
        'preferred_BIC': 'SIR' if sir_bic < seir_bic else 'SEIR',
        'delta_AIC': abs(sir_aic - seir_aic),
        'delta_BIC': abs(sir_bic - seir_bic),
    }


# ============================================================
# ABC (Approximate Bayesian Computation)
# ============================================================

def abc_rejection_seir(observed, t, N, n_particles=1000, epsilon=0.1, seed=42):
    """ABC rejection sampler for SEIR model."""
    rng = np.random.RandomState(seed)
    accepted = []

    for _ in range(n_particles * 10):
        beta = rng.lognormal(np.log(0.35), 0.5)
        sigma = rng.lognormal(np.log(0.2), 0.3)
        gamma = rng.lognormal(np.log(0.1), 0.3)

        y0 = [N - 100, 50, 50, 0]
        try:
            sol = odeint(lambda y, t_: [-beta*y[0]*y[2]/N,
                                         beta*y[0]*y[2]/N - sigma*y[1],
                                         sigma*y[1] - gamma*y[2],
                                         gamma*y[2]], y0, t)
            I_sim = sol[:, 2]
            new_cases_sim = np.diff(np.concatenate([[0], np.cumsum(I_sim * gamma)]))
            new_cases_sim = np.maximum(new_cases_sim, 0)

            # Summary statistics distance
            obs_norm = observed / (np.max(observed) + 1)
            sim_norm = new_cases_sim[:len(observed)] / (np.max(new_cases_sim[:len(observed)]) + 1)
            distance = np.sqrt(np.mean((obs_norm - sim_norm)**2))

            if distance < epsilon:
                accepted.append({'beta': beta, 'sigma': sigma, 'gamma': gamma,
                               'R0': beta/gamma, 'distance': distance})
        except Exception:
            continue

        if len(accepted) >= n_particles:
            break

    return accepted


# ============================================================
# Particle Filter (Bootstrap)
# ============================================================

def particle_filter_seir(observed, N, n_particles=500, seed=42):
    """Bootstrap particle filter for SEIR state estimation."""
    rng = np.random.RandomState(seed)
    n_days = len(observed)

    # Initialize particles
    particles = np.zeros((n_particles, 4))  # S, E, I, R
    particles[:, 0] = N - 200
    particles[:, 1] = rng.poisson(50, n_particles)
    particles[:, 2] = rng.poisson(100, n_particles)
    particles[:, 3] = 0

    # Fixed parameters
    beta, sigma, gamma = 0.35, 0.2, 0.1

    filtered_states = np.zeros((n_days, 4))
    effective_sample_sizes = np.zeros(n_days)

    for day in range(n_days):
        # Propagate
        for p in range(n_particles):
            S, E, I, R = particles[p]
            new_exposed = rng.binomial(max(int(S), 0), min(beta * max(I, 0) / N, 1))
            new_infected = rng.binomial(max(int(E), 0), min(sigma, 1))
            new_recovered = rng.binomial(max(int(I), 0), min(gamma, 1))
            particles[p] = [S - new_exposed, E + new_exposed - new_infected,
                           I + new_infected - new_recovered, R + new_recovered]

        # Weight by likelihood
        I_particles = np.maximum(particles[:, 2], 1)
        log_weights = -0.5 * ((observed[day] - I_particles * gamma) ** 2) / (I_particles * gamma + 1)
        log_weights -= np.max(log_weights)
        weights = np.exp(log_weights)
        weights /= np.sum(weights)

        effective_sample_sizes[day] = 1 / np.sum(weights**2)
        filtered_states[day] = np.average(particles, weights=weights, axis=0)

        # Resample
        indices = rng.choice(n_particles, n_particles, p=weights)
        particles = particles[indices] + rng.normal(0, 10, particles.shape)
        particles = np.maximum(particles, 0)

    return filtered_states, effective_sample_sizes
