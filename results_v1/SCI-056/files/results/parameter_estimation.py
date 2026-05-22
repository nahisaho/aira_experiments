from __future__ import annotations

"""Parameter estimation utilities for epidemiological models.

This module collects three complementary Bayesian-style estimation strategies plus
basic model selection utilities for compartmental epidemic models such as SIR and
SEIR.

Method selection guide
----------------------
- MCMC (``MCMCEstimator``): the gold standard when you can evaluate the likelihood.
  It targets the exact posterior distribution asymptotically and is the right
  default for ODE-based epidemic models with a tractable observation model.
- Particle filtering (``ParticleFilter``): best for state-space models with latent
  states, noisy observations, and sequential/real-time updates. It estimates hidden
  states and a particle marginal likelihood for PMCMC workflows.
- ABC (``ABCEstimator``): use when the likelihood is unavailable or too expensive to
  derive, for example with agent-based models or complex stochastic simulators.
  Inference is simulation-based and depends on chosen summary statistics.
- Model selection (``ModelSelector``): WAIC and PSIS-LOO are useful for predictive
  within-dataset comparison, whereas Bayes factors target between-model evidence and
  are more sensitive to priors.

The implementations are dependency-light and rely only on NumPy plus the Python
standard library so they remain portable in restricted environments.
"""

from dataclasses import dataclass
from math import lgamma
from typing import Any, Callable, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple
import json
import statistics

import numpy as np

Array = np.ndarray
PriorSpec = Mapping[str, Any]


def _to_numpy(data: Any) -> Array:
    return np.asarray(data, dtype=float)


def _as_2d(data: Any) -> Array:
    arr = _to_numpy(data)
    if arr.ndim == 1:
        return arr[:, None]
    return arr


def _logsumexp(values: Array, axis: Optional[int] = None) -> Array:
    values = np.asarray(values, dtype=float)
    max_v = np.max(values, axis=axis, keepdims=True)
    stabilized = np.exp(values - max_v)
    summed = np.sum(stabilized, axis=axis, keepdims=True)
    result = max_v + np.log(summed + 1e-300)
    if axis is None:
        return np.asarray(result).reshape(())
    return np.squeeze(result, axis=axis)


def _normal_logpdf(x: Array, mean: float, sd: float) -> Array:
    sd = max(float(sd), 1e-12)
    return -0.5 * np.log(2.0 * np.pi * sd * sd) - 0.5 * ((x - mean) / sd) ** 2


def _gamma_logpdf(x: float, shape: float, rate: float) -> float:
    if x <= 0 or shape <= 0 or rate <= 0:
        return -np.inf
    return shape * np.log(rate) - lgamma(shape) + (shape - 1.0) * np.log(x) - rate * x


def _prior_sample(spec: PriorSpec, rng: np.random.Generator) -> float:
    if hasattr(spec, "sample"):
        return float(spec.sample(rng))
    if callable(spec.get("sample")):
        return float(spec["sample"](rng))
    dist = spec.get("dist", "normal")
    if dist == "normal":
        return float(rng.normal(spec.get("mu", 0.0), spec.get("sigma", 1.0)))
    if dist == "uniform":
        return float(rng.uniform(spec.get("low", 0.0), spec.get("high", 1.0)))
    if dist == "lognormal":
        return float(rng.lognormal(spec.get("mu", 0.0), spec.get("sigma", 1.0)))
    if dist == "gamma":
        return float(rng.gamma(spec.get("shape", 2.0), 1.0 / spec.get("rate", 1.0)))
    raise ValueError(f"Unsupported prior distribution: {dist}")


def _prior_logpdf(spec: PriorSpec, x: float) -> float:
    if hasattr(spec, "logpdf"):
        return float(spec.logpdf(x))
    if callable(spec.get("logpdf")):
        return float(spec["logpdf"](x))
    dist = spec.get("dist", "normal")
    if dist == "normal":
        return float(_normal_logpdf(np.asarray(x), spec.get("mu", 0.0), spec.get("sigma", 1.0)))
    if dist == "uniform":
        low, high = spec.get("low", 0.0), spec.get("high", 1.0)
        if x < low or x > high:
            return -np.inf
        return -np.log(high - low + 1e-300)
    if dist == "lognormal":
        if x <= 0:
            return -np.inf
        mu, sigma = spec.get("mu", 0.0), max(spec.get("sigma", 1.0), 1e-12)
        return -np.log(x) + float(_normal_logpdf(np.log(x), mu, sigma))
    if dist == "gamma":
        return _gamma_logpdf(x, spec.get("shape", 2.0), spec.get("rate", 1.0))
    raise ValueError(f"Unsupported prior distribution: {dist}")


def _hdi(samples: Array, cred_mass: float = 0.95) -> Tuple[float, float]:
    values = np.sort(np.asarray(samples, dtype=float).ravel())
    if values.size == 0:
        return (np.nan, np.nan)
    interval = max(1, int(np.floor(cred_mass * values.size)))
    n_intervals = values.size - interval
    if n_intervals <= 0:
        return float(values[0]), float(values[-1])
    widths = values[interval:] - values[:n_intervals]
    idx = int(np.argmin(widths))
    return float(values[idx]), float(values[idx + interval])


def _autocorrelation_1d(x: Array, lag: int) -> float:
    x = np.asarray(x, dtype=float)
    if lag <= 0:
        return 1.0
    if lag >= x.size:
        return 0.0
    x_centered = x - np.mean(x)
    denom = np.dot(x_centered, x_centered)
    if denom <= 0:
        return 0.0
    return float(np.dot(x_centered[:-lag], x_centered[lag:]) / denom)


def _effective_sample_size(chains: Array) -> Array:
    chains = np.asarray(chains, dtype=float)
    if chains.ndim == 2:
        chains = chains[:, :, None]
    m, n, p = chains.shape
    ess = np.empty(p, dtype=float)
    for j in range(p):
        chain_vals = chains[:, :, j]
        rho_sum = 0.0
        for lag in range(1, max(2, n // 2)):
            rho_lag = np.mean([_autocorrelation_1d(chain_vals[c], lag) for c in range(m)])
            if rho_lag <= 0:
                break
            rho_sum += rho_lag
        ess[j] = m * n / max(1.0 + 2.0 * rho_sum, 1e-9)
    return ess


def _gelman_rubin(chains: Array) -> Array:
    chains = np.asarray(chains, dtype=float)
    if chains.ndim == 2:
        chains = chains[:, :, None]
    m, n, p = chains.shape
    if m < 2:
        return np.full(p, np.nan)
    chain_means = np.mean(chains, axis=1)
    chain_vars = np.var(chains, axis=1, ddof=1)
    w = np.mean(chain_vars, axis=0)
    b = n * np.var(chain_means, axis=0, ddof=1)
    var_hat = (n - 1.0) / n * w + b / n
    return np.sqrt(np.maximum(var_hat / np.maximum(w, 1e-12), 1e-12))


def _weighted_quantile(values: Array, weights: Array, quantile: float) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cumulative = np.cumsum(weights) / np.sum(weights)
    idx = np.searchsorted(cumulative, quantile)
    idx = min(max(idx, 0), len(values) - 1)
    return float(values[idx])


def _dict_samples_to_matrix(samples: Sequence[Mapping[str, float]], param_names: Sequence[str]) -> Array:
    return np.asarray([[sample[name] for name in param_names] for sample in samples], dtype=float)


@dataclass
class PosteriorSummary:
    mean: float
    median: float
    hdi_2_5: float
    hdi_97_5: float


class MCMCEstimator:
    """Metropolis-Hastings parameter estimation for deterministic epidemic ODE models.

    Parameters
    ----------
    ode_model:
        Callable with signature ``ode_model(params, times, initial_state)`` returning
        a simulated trajectory as a NumPy array.
    observed_data:
        Observed epidemic time series. For multi-output models, use
        ``observed_indices`` to select the observed compartments from the simulated
        trajectory.
    priors:
        Mapping from parameter name to prior specification. Built-in prior specs can
        use ``dist`` values ``normal``, ``uniform``, ``lognormal``, or ``gamma``.
    likelihood:
        ``"poisson"`` for count data or ``"gaussian"`` for approximately continuous
        observations.

    Notes
    -----
    MCMC is the gold standard when the likelihood is available. It targets the exact
    posterior in the large-sample limit, supports direct uncertainty quantification,
    and is especially appropriate for ODE-based SIR/SEIR models with a tractable
    observation model.

    PyMC-style setup template
    -------------------------
    ``pymc_wrapper_template()`` returns commented code showing how the same model
    would be specified in PyMC. The actual implementation here is written from
    scratch for portability and does not depend on PyMC.
    """

    def __init__(
        self,
        ode_model: Callable[[Mapping[str, float], Array, Array], Array],
        observed_data: Any,
        priors: Mapping[str, PriorSpec],
        times: Optional[Sequence[float]] = None,
        initial_state: Optional[Sequence[float]] = None,
        likelihood: str = "poisson",
        observed_indices: Optional[Sequence[int] | int] = None,
        sigma: float = 1.0,
    ) -> None:
        self.ode_model = ode_model
        self.observed_data = _as_2d(observed_data)
        self.times = np.asarray(times if times is not None else np.arange(len(self.observed_data)), dtype=float)
        self.initial_state = _to_numpy(initial_state) if initial_state is not None else np.array([])
        self.priors = dict(priors)
        self.param_names = list(priors.keys())
        self.likelihood = likelihood
        self.observed_indices = observed_indices
        self.sigma = sigma
        self.last_fit: Dict[str, Any] = {}

    @staticmethod
    def pymc_wrapper_template() -> str:
        """Return commented PyMC code showing an equivalent model specification."""

        return """\
# Example PyMC formulation for the same estimator concept:
# import pymc as pm
# import pytensor.tensor as pt
#
# with pm.Model() as model:
#     beta = pm.LogNormal('beta', mu=np.log(0.4), sigma=0.4)
#     gamma = pm.LogNormal('gamma', mu=np.log(0.2), sigma=0.4)
#     sigma = pm.LogNormal('sigma', mu=np.log(0.3), sigma=0.4)
#
#     theta = pt.stack([beta, sigma, gamma])
#     mu = seir_solver_op(theta)  # custom PyTensor Op wrapping the ODE solver
#     y = pm.Poisson('y', mu=pt.clip(mu[:, 2], 1e-6, np.inf), observed=observed_cases)
#
#     trace = pm.sample(tune=1000, draws=2000, chains=4, target_accept=0.9)
"""

    def _vector_to_params(self, vector: Sequence[float]) -> Dict[str, float]:
        return {name: float(value) for name, value in zip(self.param_names, vector)}

    def _initial_vector(self, rng: np.random.Generator) -> Array:
        vector = []
        for name in self.param_names:
            spec = self.priors[name]
            if "init" in spec:
                value = float(spec["init"])
            else:
                for _ in range(100):
                    value = _prior_sample(spec, rng)
                    if np.isfinite(_prior_logpdf(spec, value)):
                        break
                else:
                    raise RuntimeError(f"Could not initialize parameter {name} from prior")
            vector.append(value)
        return np.asarray(vector, dtype=float)

    def _simulate(self, params: Mapping[str, float]) -> Array:
        simulated = _as_2d(self.ode_model(params, self.times, self.initial_state))
        if self.observed_indices is None:
            if simulated.shape[1] == self.observed_data.shape[1]:
                return simulated
            if self.observed_data.shape[1] == 1:
                idx = min(simulated.shape[1] - 1, 2)
                return simulated[:, [idx]]
            raise ValueError("Observed data shape does not match simulation output")
        indices = [self.observed_indices] if isinstance(self.observed_indices, int) else list(self.observed_indices)
        return simulated[:, indices]

    def pointwise_log_likelihood(self, params: Mapping[str, float]) -> Array:
        simulated = np.clip(self._simulate(params), 1e-9, None)
        observed = self.observed_data
        if self.likelihood == "poisson":
            log_fact = np.vectorize(lgamma)(observed + 1.0)
            return observed * np.log(simulated) - simulated - log_fact
        if self.likelihood == "gaussian":
            return _normal_logpdf(observed, simulated, self.sigma)
        raise ValueError(f"Unsupported likelihood: {self.likelihood}")

    def log_likelihood(self, params: Mapping[str, float]) -> float:
        ll = self.pointwise_log_likelihood(params)
        if not np.all(np.isfinite(ll)):
            return -np.inf
        return float(np.sum(ll))

    def log_prior(self, params: Mapping[str, float]) -> float:
        value = 0.0
        for name, x in params.items():
            lp = _prior_logpdf(self.priors[name], float(x))
            if not np.isfinite(lp):
                return -np.inf
            value += lp
        return float(value)

    def log_posterior(self, params: Mapping[str, float]) -> float:
        lp = self.log_prior(params)
        if not np.isfinite(lp):
            return -np.inf
        ll = self.log_likelihood(params)
        return lp + ll

    def _run_single_chain(
        self,
        rng: np.random.Generator,
        n_samples: int,
        warmup: int,
        initial: Optional[Array],
        proposal_scale: Optional[Array],
        adapt_interval: int,
        target_accept: float,
    ) -> Dict[str, Any]:
        total_steps = n_samples + warmup
        position = np.array(initial if initial is not None else self._initial_vector(rng), dtype=float)
        log_post = self.log_posterior(self._vector_to_params(position))
        if not np.isfinite(log_post):
            position = self._initial_vector(rng)
            log_post = self.log_posterior(self._vector_to_params(position))
        scales = np.array(proposal_scale if proposal_scale is not None else np.maximum(np.abs(position) * 0.1, 0.05), dtype=float)
        trace = np.zeros((n_samples, len(position)), dtype=float)
        log_posts = np.zeros(n_samples, dtype=float)
        accepted = 0
        accept_window = 0
        window_history: List[Array] = []
        for step in range(total_steps):
            proposal = position + rng.normal(0.0, scales, size=position.size)
            proposal_params = self._vector_to_params(proposal)
            proposal_log_post = self.log_posterior(proposal_params)
            if np.isfinite(proposal_log_post):
                log_alpha = proposal_log_post - log_post
                if np.log(rng.uniform()) < min(0.0, log_alpha):
                    position = proposal
                    log_post = proposal_log_post
                    accepted += 1
                    accept_window += 1
            if step < warmup:
                window_history.append(position.copy())
                if (step + 1) % max(adapt_interval, 1) == 0:
                    acc_rate = accept_window / max(adapt_interval, 1)
                    scales *= np.exp(acc_rate - target_accept)
                    empirical = np.std(np.asarray(window_history), axis=0, ddof=1)
                    empirical = np.nan_to_num(empirical, nan=0.0, posinf=0.0, neginf=0.0)
                    scales = 0.5 * scales + 0.5 * np.maximum(2.38 * empirical / np.sqrt(position.size + 1e-12), 1e-3)
                    accept_window = 0
                    window_history = []
            else:
                idx = step - warmup
                trace[idx] = position
                log_posts[idx] = log_post
        return {
            "samples": trace,
            "log_posterior": log_posts,
            "acceptance_rate": accepted / max(total_steps, 1),
            "proposal_scale": scales,
        }

    def sample(
        self,
        n_samples: int = 1000,
        warmup: int = 500,
        n_chains: int = 4,
        seed: int = 42,
        proposal_scale: Optional[Sequence[float]] = None,
        adapt_interval: int = 50,
        target_accept: float = 0.234,
    ) -> Dict[str, Any]:
        """Run multiple Metropolis-Hastings chains with adaptive warmup tuning."""

        chains: List[Array] = []
        log_posteriors: List[Array] = []
        acceptance_rates: List[float] = []
        proposal = None if proposal_scale is None else np.asarray(proposal_scale, dtype=float)
        for chain_id in range(n_chains):
            rng = np.random.default_rng(seed + chain_id)
            result = self._run_single_chain(
                rng=rng,
                n_samples=n_samples,
                warmup=warmup,
                initial=None,
                proposal_scale=proposal,
                adapt_interval=adapt_interval,
                target_accept=target_accept,
            )
            chains.append(result["samples"])
            log_posteriors.append(result["log_posterior"])
            acceptance_rates.append(float(result["acceptance_rate"]))
        chains_arr = np.asarray(chains, dtype=float)
        diagnostics = convergence_report(chains_arr, self.param_names)
        summary = self.posterior_summary(chains_arr)
        trace_data = self.trace_plot_data(chains_arr)
        posterior_samples = [self._vector_to_params(row) for row in chains_arr.reshape(-1, chains_arr.shape[-1])]
        self.last_fit = {
            "chains": chains_arr,
            "log_posterior": np.asarray(log_posteriors),
            "acceptance_rates": acceptance_rates,
            "diagnostics": diagnostics,
            "summary": summary,
            "trace_data": trace_data,
            "posterior_samples": posterior_samples,
        }
        return self.last_fit

    def posterior_summary(self, chains: Array) -> Dict[str, Dict[str, float]]:
        flat = np.asarray(chains, dtype=float).reshape(-1, len(self.param_names))
        output: Dict[str, Dict[str, float]] = {}
        for idx, name in enumerate(self.param_names):
            values = flat[:, idx]
            hdi_low, hdi_high = _hdi(values)
            output[name] = {
                "mean": float(np.mean(values)),
                "median": float(np.median(values)),
                "hdi_2.5%": hdi_low,
                "hdi_97.5%": hdi_high,
            }
        return output

    def trace_plot_data(self, chains: Array) -> Dict[str, Dict[str, Any]]:
        chains = np.asarray(chains, dtype=float)
        x = np.arange(chains.shape[1]).tolist()
        data: Dict[str, Dict[str, Any]] = {}
        for idx, name in enumerate(self.param_names):
            data[name] = {
                "x": x,
                "chains": [chains[c, :, idx].tolist() for c in range(chains.shape[0])],
            }
        return data


class ParticleFilter:
    """Bootstrap particle filter for latent epidemic state estimation.

    Particle filtering is preferable when the epidemic process is naturally treated as
    a state-space model with latent states, streaming observations, or online
    forecasting needs. It supports particle marginal likelihood estimation and a
    simple PMCMC wrapper for joint state-parameter inference.
    """

    def __init__(
        self,
        transition_function: Callable[[Array, int, Mapping[str, float], np.random.Generator], Array],
        observation_loglikelihood: Callable[[float | Array, Array, Mapping[str, float], int], Array],
        n_particles: int = 200,
        resampling_threshold: float = 0.5,
    ) -> None:
        self.transition_function = transition_function
        self.observation_loglikelihood = observation_loglikelihood
        self.n_particles = int(n_particles)
        self.resampling_threshold = float(resampling_threshold)

    @staticmethod
    def systematic_resampling(weights: Array, rng: np.random.Generator) -> Array:
        weights = np.asarray(weights, dtype=float)
        n = weights.size
        positions = (rng.uniform() + np.arange(n)) / n
        cumulative = np.cumsum(weights)
        indexes = np.zeros(n, dtype=int)
        i = j = 0
        while i < n:
            if positions[i] < cumulative[j]:
                indexes[i] = j
                i += 1
            else:
                j += 1
        return indexes

    @staticmethod
    def effective_sample_size(weights: Array) -> float:
        weights = np.asarray(weights, dtype=float)
        return float(1.0 / np.sum(np.square(weights) + 1e-300))

    def _propagate(self, particles: Array, t: int, params: Mapping[str, float], rng: np.random.Generator) -> Array:
        try:
            proposed = self.transition_function(particles, t, params, rng)
        except Exception:
            proposed = np.asarray([self.transition_function(p, t, params, rng) for p in particles], dtype=float)
        return np.asarray(proposed, dtype=float)

    def run(
        self,
        observations: Sequence[float] | Array,
        initial_particles: Array | Callable[[np.random.Generator, int], Array],
        params: Mapping[str, float],
        seed: int = 42,
    ) -> Dict[str, Any]:
        rng = np.random.default_rng(seed)
        observations_arr = np.asarray(observations, dtype=float)
        particles = (
            np.asarray(initial_particles(rng, self.n_particles), dtype=float)
            if callable(initial_particles)
            else np.asarray(initial_particles, dtype=float)
        )
        if particles.shape[0] != self.n_particles:
            raise ValueError("initial_particles must contain n_particles rows")
        weights = np.full(self.n_particles, 1.0 / self.n_particles, dtype=float)
        filtered_states: List[Array] = []
        ess_history: List[float] = []
        log_marginal_likelihood = 0.0
        particle_history: List[Array] = []
        for t, observation in enumerate(observations_arr):
            particles = self._propagate(particles, t, params, rng)
            logw = np.asarray(self.observation_loglikelihood(observation, particles, params, t), dtype=float).reshape(-1)
            max_logw = np.max(logw)
            w_unnorm = np.exp(logw - max_logw)
            normalizer = np.sum(w_unnorm)
            log_marginal_likelihood += float(max_logw + np.log(normalizer + 1e-300) - np.log(self.n_particles))
            weights = w_unnorm / (normalizer + 1e-300)
            ess = self.effective_sample_size(weights)
            ess_history.append(ess)
            filtered_states.append(np.average(particles, axis=0, weights=weights))
            particle_history.append(particles.copy())
            if ess < self.resampling_threshold * self.n_particles:
                idx = self.systematic_resampling(weights, rng)
                particles = particles[idx]
                weights = np.full(self.n_particles, 1.0 / self.n_particles)
        return {
            "filtered_states": np.asarray(filtered_states, dtype=float),
            "marginal_likelihood": float(log_marginal_likelihood),
            "ess_history": ess_history,
            "particle_history": particle_history,
            "weights": weights,
        }

    def pmcmc(
        self,
        param_priors: Mapping[str, PriorSpec],
        observations: Sequence[float] | Array,
        initial_particles: Array | Callable[[np.random.Generator, int], Array],
        n_iterations: int = 300,
        warmup: int = 100,
        proposal_scale: Optional[Mapping[str, float]] = None,
        seed: int = 42,
    ) -> Dict[str, Any]:
        """Particle marginal Metropolis-Hastings wrapper for joint state-parameter inference."""

        rng = np.random.default_rng(seed)
        names = list(param_priors.keys())
        scales = {name: (proposal_scale or {}).get(name, 0.05) for name in names}
        current = {name: _prior_sample(param_priors[name], rng) for name in names}
        current_run = self.run(observations, initial_particles, current, seed=seed)
        current_log_like = current_run["marginal_likelihood"]
        current_log_prior = sum(_prior_logpdf(param_priors[name], current[name]) for name in names)
        posterior_samples: List[Dict[str, float]] = []
        filtered_states: List[Array] = []
        loglik_trace: List[float] = []
        accepted = 0
        for iteration in range(n_iterations):
            proposal = {name: current[name] + rng.normal(0.0, scales[name]) for name in names}
            proposal_log_prior = sum(_prior_logpdf(param_priors[name], proposal[name]) for name in names)
            if np.isfinite(proposal_log_prior):
                proposed_run = self.run(observations, initial_particles, proposal, seed=seed + iteration + 1)
                proposal_log_like = proposed_run["marginal_likelihood"]
                log_alpha = (proposal_log_like + proposal_log_prior) - (current_log_like + current_log_prior)
                if np.log(rng.uniform()) < min(0.0, log_alpha):
                    current = proposal
                    current_log_like = proposal_log_like
                    current_log_prior = proposal_log_prior
                    current_run = proposed_run
                    accepted += 1
            if iteration >= warmup:
                posterior_samples.append(dict(current))
                filtered_states.append(current_run["filtered_states"])
                loglik_trace.append(current_log_like)
        matrix = _dict_samples_to_matrix(posterior_samples, names) if posterior_samples else np.empty((0, len(names)))
        return {
            "posterior_samples": posterior_samples,
            "parameter_posteriors": {name: matrix[:, i].tolist() for i, name in enumerate(names)} if posterior_samples else {},
            "filtered_states": filtered_states[-1] if filtered_states else np.empty((0, 0)),
            "marginal_likelihood_trace": loglik_trace,
            "acceptance_rate": accepted / max(n_iterations, 1),
        }


class ABCEstimator:
    """Approximate Bayesian Computation with a sequential Monte Carlo schedule.

    ABC-SMC is appropriate when simulation is easy but the likelihood is intractable,
    such as in agent-based epidemic models, stochastic contact networks, or bespoke
    simulators with complex observation processes.
    """

    def __init__(
        self,
        simulator: Callable[[Mapping[str, float]], Array],
        priors: Mapping[str, PriorSpec],
        population_size: int = 200,
        n_generations: int = 4,
        summary_weights: Optional[Mapping[str, float]] = None,
        seed: int = 42,
    ) -> None:
        self.simulator = simulator
        self.priors = dict(priors)
        self.param_names = list(priors.keys())
        self.population_size = int(population_size)
        self.n_generations = int(n_generations)
        self.summary_weights = dict(summary_weights or {
            "peak_timing": 1.0,
            "peak_height": 1.0,
            "final_size": 1.0,
            "r0_estimate": 1.0,
        })
        self.seed = seed

    @staticmethod
    def _extract_curve(data: Array) -> Array:
        arr = _as_2d(data)
        if arr.shape[1] == 1:
            return arr[:, 0]
        idx = min(arr.shape[1] - 1, 2)
        return arr[:, idx]

    def summary_statistics(self, series: Array, params: Optional[Mapping[str, float]] = None) -> Dict[str, float]:
        arr = _as_2d(series)
        curve = np.maximum(self._extract_curve(arr), 0.0)
        peak_idx = int(np.argmax(curve))
        peak_height = float(np.max(curve))
        if arr.shape[1] >= 2 and np.all(arr[0] >= 0):
            total_pop = float(np.sum(arr[0]))
            final_size = float(max(total_pop - arr[-1, 0], 0.0))
        else:
            final_size = float(np.sum(curve))
        if params and "beta" in params and "gamma" in params and params["gamma"] != 0:
            r0_estimate = float(params["beta"] / params["gamma"])
        else:
            start = max(curve[0], 1e-9)
            stop = max(curve[min(peak_idx, len(curve) - 1)], 1e-9)
            r0_estimate = float(max(1.0, 1.0 + np.log(stop / start + 1e-9) / max(peak_idx + 1, 1)))
        return {
            "peak_timing": float(peak_idx),
            "peak_height": peak_height,
            "final_size": final_size,
            "r0_estimate": r0_estimate,
        }

    def distance(self, simulated_stats: Mapping[str, float], observed_stats: Mapping[str, float]) -> float:
        total = 0.0
        for key, weight in self.summary_weights.items():
            diff = simulated_stats[key] - observed_stats[key]
            total += float(weight) * diff * diff
        return float(np.sqrt(total))

    def _kernel_pdf(self, point: Mapping[str, float], center: Mapping[str, float], scales: Mapping[str, float]) -> float:
        density = 1.0
        for name in self.param_names:
            sd = max(scales[name], 1e-6)
            density *= float(np.exp(_normal_logpdf(point[name], center[name], sd)))
        return density

    def fit(
        self,
        observed_data: Array,
        epsilon_quantile: float = 0.5,
        max_attempts: int = 100000,
    ) -> Dict[str, Any]:
        rng = np.random.default_rng(self.seed)
        observed_stats = self.summary_statistics(observed_data)
        populations: List[List[Dict[str, float]]] = []
        weights_history: List[List[float]] = []
        epsilon_history: List[float] = []
        acceptance_rates: List[float] = []
        previous_population: List[Dict[str, float]] = []
        previous_weights: Array = np.array([])
        kernel_scales = {name: 0.1 for name in self.param_names}
        for generation in range(self.n_generations):
            accepted: List[Dict[str, float]] = []
            accepted_weights: List[float] = []
            distances: List[float] = []
            attempts = 0
            epsilon = np.inf if generation == 0 else epsilon_history[-1]
            while len(accepted) < self.population_size and attempts < max_attempts:
                attempts += 1
                if generation == 0:
                    candidate = {name: _prior_sample(self.priors[name], rng) for name in self.param_names}
                else:
                    idx = rng.choice(len(previous_population), p=previous_weights)
                    anchor = previous_population[int(idx)]
                    candidate = {name: anchor[name] + rng.normal(0.0, kernel_scales[name]) for name in self.param_names}
                    if not np.isfinite(sum(_prior_logpdf(self.priors[name], candidate[name]) for name in self.param_names)):
                        continue
                simulated = self.simulator(candidate)
                sim_stats = self.summary_statistics(simulated, candidate)
                dist = self.distance(sim_stats, observed_stats)
                distances.append(dist)
                if generation == 0 and len(distances) >= self.population_size:
                    epsilon = float(np.quantile(distances, epsilon_quantile))
                if dist <= epsilon:
                    accepted.append(candidate)
                    if generation == 0:
                        accepted_weights.append(1.0)
                    else:
                        numerator = np.exp(sum(_prior_logpdf(self.priors[name], candidate[name]) for name in self.param_names))
                        denominator = 0.0
                        for w, previous in zip(previous_weights, previous_population):
                            denominator += float(w) * self._kernel_pdf(candidate, previous, kernel_scales)
                        accepted_weights.append(float(numerator / max(denominator, 1e-300)))
            if not accepted:
                raise RuntimeError("ABC-SMC failed to accept any particles")
            distances_arr = np.asarray(
                [self.distance(self.summary_statistics(self.simulator(theta), theta), observed_stats) for theta in accepted],
                dtype=float,
            )
            order = np.argsort(distances_arr)[: self.population_size]
            accepted = [accepted[i] for i in order]
            accepted_weights = [accepted_weights[i] for i in order]
            distances_arr = distances_arr[order]
            weights = np.asarray(accepted_weights, dtype=float)
            weights /= np.sum(weights)
            previous_population = accepted
            previous_weights = weights
            populations.append(accepted)
            weights_history.append(weights.tolist())
            epsilon = float(np.quantile(distances_arr, epsilon_quantile))
            epsilon_history.append(epsilon)
            acceptance_rates.append(len(accepted) / max(attempts, 1))
            matrix = _dict_samples_to_matrix(accepted, self.param_names)
            spread = np.std(matrix, axis=0, ddof=1) if len(accepted) > 1 else np.ones(len(self.param_names)) * 0.1
            kernel_scales = {name: float(max(2.0 * spread[i], 1e-3)) for i, name in enumerate(self.param_names)}
        final_population = populations[-1]
        final_weights = np.asarray(weights_history[-1], dtype=float)
        return {
            "posterior_samples": final_population,
            "weights": final_weights.tolist(),
            "acceptance_rates": acceptance_rates,
            "epsilon_history": epsilon_history,
            "populations": populations,
            "observed_summary_statistics": observed_stats,
        }


class ModelSelector:
    """Posterior predictive and evidence-based model comparison helpers.

    Interpretation guidelines
    -------------------------
    - WAIC: lower is better. Differences <2 are usually weak, 4-7 indicate moderate
      support, and >10 often suggest a practically important predictive difference.
    - PSIS-LOO: use the Pareto-k diagnostics to assess reliability. Values below 0.5
      are generally safe, 0.5-0.7 deserve caution, and >0.7 indicate unstable
      importance sampling.
    - Bayes factors: values above 3 are often described as moderate evidence and above
      10 as strong evidence, but results are sensitive to the prior and the marginal
      likelihood estimator.
    """

    @staticmethod
    def waic(log_likelihood_matrix: Array) -> Dict[str, float]:
        ll = np.asarray(log_likelihood_matrix, dtype=float)
        lppd = float(np.sum(_logsumexp(ll, axis=0) - np.log(ll.shape[0])))
        p_waic = float(np.sum(np.var(ll, axis=0, ddof=1)))
        waic = -2.0 * (lppd - p_waic)
        return {"waic": waic, "lppd": lppd, "p_waic": p_waic}

    @staticmethod
    def _psis_smooth(weights: Array) -> Tuple[Array, float]:
        weights = np.asarray(weights, dtype=float)
        n = weights.size
        if n < 5:
            return weights / np.sum(weights), np.nan
        order = np.argsort(weights)
        tail_size = max(3, int(0.2 * n))
        threshold = weights[order[-tail_size]]
        tail = weights[order[-tail_size:]]
        if np.any(tail <= 0):
            smoothed = np.minimum(weights, np.quantile(weights, 0.9))
            return smoothed / np.sum(smoothed), np.nan
        hill = np.mean(np.log(tail) - np.log(threshold + 1e-300))
        k_hat = float(max(hill, 0.0))
        u = (np.arange(1, tail_size + 1) - 0.5) / tail_size
        if k_hat < 1e-8:
            smooth_tail = threshold * np.exp(u)
        else:
            smooth_tail = threshold * np.power(1.0 - u, -k_hat)
        smoothed = weights.copy()
        smoothed[order[-tail_size:]] = smooth_tail
        smoothed /= np.sum(smoothed)
        return smoothed, k_hat

    @classmethod
    def psis_loo(cls, log_likelihood_matrix: Array) -> Dict[str, Any]:
        ll = np.asarray(log_likelihood_matrix, dtype=float)
        n_draws, n_obs = ll.shape
        elpd_terms: List[float] = []
        pareto_k: List[float] = []
        for i in range(n_obs):
            raw_weights = np.exp(-ll[:, i] - np.max(-ll[:, i]))
            smoothed, k_hat = cls._psis_smooth(raw_weights)
            pareto_k.append(k_hat)
            elpd_terms.append(float(_logsumexp(ll[:, i] + np.log(smoothed + 1e-300))))
        elpd_loo = float(np.sum(elpd_terms))
        looic = -2.0 * elpd_loo
        return {"elpd_loo": elpd_loo, "looic": looic, "pareto_k": pareto_k}

    @staticmethod
    def harmonic_mean_log_marginal(log_likelihood_matrix: Array) -> float:
        total_log_lik = np.sum(np.asarray(log_likelihood_matrix, dtype=float), axis=1)
        return float(-_logsumexp(-total_log_lik) + np.log(len(total_log_lik)))

    @staticmethod
    def bridge_sampling_concept(log_likelihood_matrix: Array) -> Dict[str, Any]:
        total_log_lik = np.sum(np.asarray(log_likelihood_matrix, dtype=float), axis=1)
        center = float(np.mean(total_log_lik))
        stabilized = total_log_lik - center
        bridge_proxy = float(center + np.log(np.mean(1.0 / (1.0 + np.exp(-stabilized)))))
        return {
            "bridge_proxy_log_evidence": bridge_proxy,
            "note": "Proxy based on a logistic bridge identity; use a full bridge sampler for publication-grade evidence estimates.",
        }

    @staticmethod
    def bayes_factor(log_marginal_a: float, log_marginal_b: float) -> float:
        return float(np.exp(log_marginal_a - log_marginal_b))

    @staticmethod
    def _resolve_log_likelihood_matrix(model_spec: Mapping[str, Any], data: Any) -> Array:
        if "log_likelihood_matrix" in model_spec:
            return np.asarray(model_spec["log_likelihood_matrix"], dtype=float)
        if "samples" not in model_spec or "log_likelihood_fn" not in model_spec:
            raise ValueError("Model spec must provide either log_likelihood_matrix or samples + log_likelihood_fn")
        rows = []
        for sample in model_spec["samples"]:
            pointwise = np.asarray(model_spec["log_likelihood_fn"](sample, data), dtype=float)
            rows.append(pointwise.ravel())
        return np.asarray(rows, dtype=float)

    @classmethod
    def compare_models(cls, models_dict: Mapping[str, Mapping[str, Any]], data: Any) -> List[Dict[str, Any]]:
        comparison: List[Dict[str, Any]] = []
        for name, spec in models_dict.items():
            ll_matrix = cls._resolve_log_likelihood_matrix(spec, data)
            waic = cls.waic(ll_matrix)
            loo = cls.psis_loo(ll_matrix)
            log_marginal = cls.harmonic_mean_log_marginal(ll_matrix)
            bridge = cls.bridge_sampling_concept(ll_matrix)
            comparison.append({
                "model": name,
                "waic": waic["waic"],
                "p_waic": waic["p_waic"],
                "looic": loo["looic"],
                "elpd_loo": loo["elpd_loo"],
                "pareto_k_max": float(np.nanmax(loo["pareto_k"])),
                "log_marginal_hmean": log_marginal,
                "bridge_proxy_log_evidence": bridge["bridge_proxy_log_evidence"],
            })
        comparison.sort(key=lambda row: (row["waic"], row["looic"]))
        if comparison:
            best = comparison[0]
            for row in comparison:
                row["delta_waic"] = float(row["waic"] - best["waic"])
                row["bayes_factor_vs_best"] = 1.0 if row is best else cls.bayes_factor(best["log_marginal_hmean"], row["log_marginal_hmean"])
        return comparison


def generate_synthetic_data(
    model: Callable[..., Array],
    true_params: Mapping[str, float],
    noise_type: str = "poisson",
    seed: int = 42,
    **model_kwargs: Any,
) -> Array:
    """Generate noisy synthetic observations from a deterministic epidemic model.

    Additional keyword arguments are passed through to ``model``. Typical examples are
    ``times``, ``initial_state``, or ``observed_index``.
    """

    rng = np.random.default_rng(seed)
    trajectory = _as_2d(model(true_params, model_kwargs.get("times"), model_kwargs.get("initial_state")))
    observed_index = model_kwargs.get("observed_index")
    if observed_index is not None:
        signal = trajectory[:, [observed_index]]
    else:
        signal = trajectory
    signal = np.clip(signal, 1e-9, None)
    if noise_type == "poisson":
        return rng.poisson(signal)
    if noise_type == "gaussian":
        sd = float(model_kwargs.get("sigma", 1.0))
        return signal + rng.normal(0.0, sd, size=signal.shape)
    raise ValueError(f"Unsupported noise_type: {noise_type}")


def plot_posterior_data(samples: Array | Mapping[str, Sequence[float]] | Sequence[Mapping[str, float]], param_names: Sequence[str]) -> Dict[str, Any]:
    """Return plotting-ready posterior summaries instead of drawing figures.

    The returned dictionary contains trace arrays, histogram counts, bin edges, and
    key quantiles so a plotting layer can remain separate from the inference code.
    """

    if isinstance(samples, Mapping):
        matrix = np.asarray([samples[name] for name in param_names], dtype=float).T
    elif isinstance(samples, Sequence) and samples and isinstance(samples[0], Mapping):
        matrix = _dict_samples_to_matrix(samples, param_names)
    else:
        matrix = np.asarray(samples, dtype=float)
        if matrix.ndim == 3:
            matrix = matrix.reshape(-1, matrix.shape[-1])
    output = {}
    for idx, name in enumerate(param_names):
        values = matrix[:, idx]
        counts, bins = np.histogram(values, bins="auto", density=True)
        output[name] = {
            "trace": values.tolist(),
            "histogram_density": counts.tolist(),
            "histogram_bins": bins.tolist(),
            "quantiles": {
                "2.5%": float(np.quantile(values, 0.025)),
                "50%": float(np.quantile(values, 0.5)),
                "97.5%": float(np.quantile(values, 0.975)),
            },
        }
    return output


def convergence_report(chains: Array | Mapping[str, Any], param_names: Optional[Sequence[str]] = None) -> Dict[str, Any]:
    """Compute convergence diagnostics for MCMC output.

    Returns R-hat, effective sample size, chain means, chain standard deviations, and
    trace-plot-ready data. Provide either a raw ``(chains, draws, parameters)`` array
    or the dictionary returned by ``MCMCEstimator.sample``.
    """

    if isinstance(chains, Mapping):
        param_names = list(chains.get("summary", {}).keys()) or list(param_names or [])
        chains_arr = np.asarray(chains["chains"], dtype=float)
    else:
        chains_arr = np.asarray(chains, dtype=float)
    if chains_arr.ndim == 2:
        chains_arr = chains_arr[:, :, None]
    n_params = chains_arr.shape[-1]
    names = list(param_names or [f"param_{i}" for i in range(n_params)])
    rhat = _gelman_rubin(chains_arr)
    ess = _effective_sample_size(chains_arr)
    report = {
        "r_hat": {name: float(rhat[i]) for i, name in enumerate(names)},
        "ess": {name: float(ess[i]) for i, name in enumerate(names)},
        "chain_means": {name: [float(np.mean(chains_arr[c, :, i])) for c in range(chains_arr.shape[0])] for i, name in enumerate(names)},
        "chain_sds": {name: [float(np.std(chains_arr[c, :, i], ddof=1)) for c in range(chains_arr.shape[0])] for i, name in enumerate(names)},
        "trace_data": {name: [chains_arr[c, :, i].tolist() for c in range(chains_arr.shape[0])] for i, name in enumerate(names)},
    }
    return report


def _euler_integrate(rhs: Callable[[Array, Mapping[str, float]], Array], initial_state: Sequence[float], times: Sequence[float], params: Mapping[str, float]) -> Array:
    times_arr = np.asarray(times, dtype=float)
    state = np.asarray(initial_state, dtype=float)
    states = [state.copy()]
    for t0, t1 in zip(times_arr[:-1], times_arr[1:]):
        dt = float(t1 - t0)
        state = np.maximum(state + dt * rhs(state, params), 0.0)
        states.append(state.copy())
    return np.asarray(states, dtype=float)


def seir_model(params: Mapping[str, float], times: Sequence[float], initial_state: Sequence[float]) -> Array:
    beta = params["beta"]
    sigma = params["sigma"]
    gamma = params["gamma"]

    def rhs(state: Array, theta: Mapping[str, float]) -> Array:
        s, e, i, r = state
        n = max(np.sum(state), 1e-9)
        force = beta * s * i / n
        return np.array([
            -force,
            force - sigma * e,
            sigma * e - gamma * i,
            gamma * i,
        ])

    return _euler_integrate(rhs, initial_state, times, params)


def sir_model(params: Mapping[str, float], times: Sequence[float], initial_state: Sequence[float]) -> Array:
    beta = params["beta"]
    gamma = params["gamma"]

    def rhs(state: Array, theta: Mapping[str, float]) -> Array:
        s, i, r = state
        n = max(np.sum(state), 1e-9)
        force = beta * s * i / n
        return np.array([
            -force,
            force - gamma * i,
            gamma * i,
        ])

    return _euler_integrate(rhs, initial_state, times, params)


def _demo_particle_transition(particles: Array, t: int, params: Mapping[str, float], rng: np.random.Generator) -> Array:
    beta, gamma = params["beta"], params["gamma"]
    if particles.ndim == 1:
        particles = particles[None, :]
    next_particles = []
    for state in particles:
        s, i, r = state
        n = max(s + i + r, 1e-9)
        force = beta * s * i / n
        drift = np.array([-force, force - gamma * i, gamma * i])
        noise = rng.normal(0.0, 0.5, size=3)
        next_particles.append(np.maximum(state + drift + noise, 0.0))
    return np.asarray(next_particles, dtype=float)


def _demo_particle_loglik(observation: float, particles: Array, params: Mapping[str, float], t: int) -> Array:
    mu = np.clip(particles[:, 1], 1e-6, None)
    obs = float(observation)
    return obs * np.log(mu) - mu - lgamma(obs + 1.0)


if __name__ == "__main__":
    times = np.arange(0, 61, 1.0)
    seir_initial = np.array([990.0, 5.0, 3.0, 2.0])
    true_params = {"beta": 0.42, "sigma": 0.22, "gamma": 0.14}
    observed = generate_synthetic_data(seir_model, true_params, noise_type="poisson", seed=42, times=times, initial_state=seir_initial, observed_index=2)

    seir_priors = {
        "beta": {"dist": "uniform", "low": 0.1, "high": 0.8, "init": 0.35},
        "sigma": {"dist": "uniform", "low": 0.05, "high": 0.5, "init": 0.2},
        "gamma": {"dist": "uniform", "low": 0.05, "high": 0.4, "init": 0.15},
    }
    sir_priors = {
        "beta": {"dist": "uniform", "low": 0.1, "high": 0.8, "init": 0.35},
        "gamma": {"dist": "uniform", "low": 0.05, "high": 0.4, "init": 0.15},
    }

    seir_estimator = MCMCEstimator(
        ode_model=seir_model,
        observed_data=observed,
        priors=seir_priors,
        times=times,
        initial_state=seir_initial,
        observed_indices=2,
        likelihood="poisson",
    )
    seir_fit = seir_estimator.sample(n_samples=120, warmup=120, n_chains=2, seed=11)
    diagnostics = convergence_report(seir_fit, seir_estimator.param_names)

    sir_initial = np.array([995.0, 3.0, 2.0])
    sir_estimator = MCMCEstimator(
        ode_model=sir_model,
        observed_data=observed,
        priors=sir_priors,
        times=times,
        initial_state=sir_initial,
        observed_indices=1,
        likelihood="poisson",
    )
    sir_fit = sir_estimator.sample(n_samples=120, warmup=120, n_chains=2, seed=17)

    selector = ModelSelector()
    comparison = selector.compare_models(
        {
            "SEIR": {
                "samples": seir_fit["posterior_samples"],
                "log_likelihood_fn": lambda sample, data: seir_estimator.pointwise_log_likelihood(sample).ravel(),
            },
            "SIR": {
                "samples": sir_fit["posterior_samples"],
                "log_likelihood_fn": lambda sample, data: sir_estimator.pointwise_log_likelihood(sample).ravel(),
            },
        },
        observed,
    )

    posterior_plot_data = plot_posterior_data(seir_fit["chains"], seir_estimator.param_names)

    print("Synthetic SEIR parameter estimation demo")
    print("Posterior summary (SEIR):")
    for name, stats in seir_fit["summary"].items():
        print(f"  {name}: mean={stats['mean']:.3f}, median={stats['median']:.3f}, 95% HDI=({stats['hdi_2.5%']:.3f}, {stats['hdi_97.5%']:.3f})")
    print("Convergence diagnostics:")
    for name in seir_estimator.param_names:
        print(f"  {name}: R-hat={diagnostics['r_hat'][name]:.3f}, ESS={diagnostics['ess'][name]:.1f}")
    print("Model comparison (lower WAIC is better):")
    for row in comparison:
        print(f"  {row['model']}: WAIC={row['waic']:.2f}, delta_WAIC={row['delta_waic']:.2f}, LOOIC={row['looic']:.2f}, BF_vs_best={row['bayes_factor_vs_best']:.2f}")
    print("Trace data keys:", list(posterior_plot_data.keys()))
