"""System identification tools for sim-to-real parameter transfer."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Callable, Mapping, Optional, Protocol, Sequence

import torch

from .domain_randomization import NumericRange


@dataclass(frozen=True)
class ParameterSpec:
    """Named simulator parameter with search bounds."""

    name: str
    bounds: NumericRange
    initial_value: Optional[float] = None

    @property
    def default_value(self) -> float:
        return self.bounds.center if self.initial_value is None else self.initial_value


@dataclass
class Trajectory:
    """Observed or simulated trajectory under optional partial observability."""

    positions: torch.Tensor
    velocities: Optional[torch.Tensor] = None
    observation_mask: Optional[torch.Tensor] = None

    def __post_init__(self) -> None:
        self.positions = _ensure_2d(self.positions)
        if self.velocities is not None:
            self.velocities = _ensure_2d(self.velocities)
            if self.velocities.shape != self.positions.shape:
                raise ValueError("positions and velocities must have identical shapes.")
        if self.observation_mask is not None:
            self.observation_mask = self.observation_mask.to(dtype=torch.bool)

    def inferred_velocities(self) -> torch.Tensor:
        if self.velocities is not None:
            return self.velocities
        if self.positions.shape[0] == 1:
            return torch.zeros_like(self.positions)
        velocity = torch.zeros_like(self.positions)
        velocity[1:] = self.positions[1:] - self.positions[:-1]
        velocity[0] = velocity[1]
        return velocity


class SimulatorProtocol(Protocol):
    """Simulator interface used by the identification pipeline."""

    def rollout(
        self,
        parameters: Mapping[str, float],
        actions: torch.Tensor,
        observation_mask: Optional[torch.Tensor] = None,
    ) -> Trajectory:
        """Run the simulator with candidate parameters."""


@dataclass(frozen=True)
class BayesianOptimizationConfig:
    """Configuration for Gaussian-process Bayesian optimization."""

    iterations: int = 40
    init_points: int = 8
    candidate_count: int = 256
    acquisition_xi: float = 0.01
    kernel_length_scale: float = 0.25
    kernel_variance: float = 1.0
    noise: float = 1e-6


@dataclass(frozen=True)
class CMAESConfig:
    """Configuration for CMA-ES parameter search."""

    generations: int = 50
    population_size: Optional[int] = None
    sigma: float = 0.25
    min_sigma: float = 1e-4


@dataclass(frozen=True)
class SystemIdentificationConfig:
    """Configuration bundle for the identification pipeline."""

    optimizer: str = "bayesian"
    position_weight: float = 1.0
    velocity_weight: float = 0.5
    bayesian: BayesianOptimizationConfig = field(default_factory=BayesianOptimizationConfig)
    cmaes: CMAESConfig = field(default_factory=CMAESConfig)


@dataclass(frozen=True)
class IdentificationResult:
    """Result of simulator parameter identification."""

    parameters: Mapping[str, float]
    loss: float
    optimizer: str
    history: Sequence[tuple[Mapping[str, float], float]]
    simulated_trajectory: Trajectory


class TrajectoryMatchingLoss:
    """Weighted trajectory matching loss for position and velocity alignment."""

    def __init__(self, position_weight: float = 1.0, velocity_weight: float = 0.5) -> None:
        self.position_weight = position_weight
        self.velocity_weight = velocity_weight

    def __call__(self, real: Trajectory, simulated: Trajectory) -> torch.Tensor:
        real_positions, sim_positions = _align_length(real.positions, simulated.positions)
        real_velocities, sim_velocities = _align_length(real.inferred_velocities(), simulated.inferred_velocities())
        mask = _broadcast_mask(real.observation_mask, real_positions.shape)
        position_loss = _masked_mse(real_positions, sim_positions, mask)
        velocity_loss = _masked_mse(real_velocities, sim_velocities, mask)
        return self.position_weight * position_loss + self.velocity_weight * velocity_loss


class GaussianProcessSurrogate:
    """Minimal Gaussian-process surrogate used by Bayesian optimization."""

    def __init__(
        self,
        length_scale: float = 0.25,
        variance: float = 1.0,
        noise: float = 1e-6,
    ) -> None:
        self.length_scale = length_scale
        self.variance = variance
        self.noise = noise
        self._x: Optional[torch.Tensor] = None
        self._y: Optional[torch.Tensor] = None
        self._chol: Optional[torch.Tensor] = None
        self._alpha: Optional[torch.Tensor] = None

    def fit(self, x: torch.Tensor, y: torch.Tensor) -> None:
        if x.ndim != 2 or y.ndim != 1:
            raise ValueError("Expected x to be 2-D and y to be 1-D.")
        self._x = x
        self._y = y
        kernel = self._kernel(x, x)
        kernel = kernel + torch.eye(x.shape[0], device=x.device, dtype=x.dtype) * self.noise
        self._chol = torch.linalg.cholesky(kernel)
        self._alpha = torch.cholesky_solve(y[:, None], self._chol).squeeze(-1)

    def predict(self, x_star: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        if self._x is None or self._y is None or self._chol is None or self._alpha is None:
            raise RuntimeError("The surrogate must be fitted before prediction.")
        cross_kernel = self._kernel(x_star, self._x)
        mean = cross_kernel @ self._alpha
        solve = torch.cholesky_solve(cross_kernel.transpose(0, 1), self._chol)
        variance = self._kernel(x_star, x_star).diagonal() - (cross_kernel * solve.transpose(0, 1)).sum(dim=1)
        return mean, variance.clamp_min(1e-9)

    def _kernel(self, x1: torch.Tensor, x2: torch.Tensor) -> torch.Tensor:
        distances = torch.cdist(x1 / self.length_scale, x2 / self.length_scale, p=2) ** 2
        return self.variance * torch.exp(-0.5 * distances)


class BayesianParameterOptimizer:
    """Bayesian optimizer over bounded simulator parameters."""

    def __init__(
        self,
        parameter_specs: Sequence[ParameterSpec],
        config: BayesianOptimizationConfig,
        *,
        seed: Optional[int] = None,
        device: Optional[torch.device] = None,
    ) -> None:
        self.parameter_specs = list(parameter_specs)
        self.config = config
        self.device = device or torch.device("cpu")
        self.generator = torch.Generator(device=self.device)
        if seed is not None:
            self.generator.manual_seed(seed)
        self.surrogate = GaussianProcessSurrogate(
            length_scale=config.kernel_length_scale,
            variance=config.kernel_variance,
            noise=config.noise,
        )

    def optimize(self, objective: Callable[[torch.Tensor], float]) -> tuple[torch.Tensor, float, list[tuple[torch.Tensor, float]]]:
        dim = len(self.parameter_specs)
        x_samples = self._sample_uniform(self.config.init_points)
        history: list[tuple[torch.Tensor, float]] = []
        y_samples = []
        for candidate in x_samples:
            value = float(objective(candidate))
            history.append((candidate.clone(), value))
            y_samples.append(value)
        x_data = x_samples.clone()
        y_data = torch.tensor(y_samples, dtype=torch.float32, device=self.device)

        for _ in range(self.config.iterations):
            self.surrogate.fit(x_data, y_data)
            candidates = self._sample_uniform(self.config.candidate_count)
            scores = self._expected_improvement(candidates, y_data.min())
            next_candidate = candidates[torch.argmax(scores)]
            next_value = float(objective(next_candidate))
            history.append((next_candidate.clone(), next_value))
            x_data = torch.cat([x_data, next_candidate.unsqueeze(0)], dim=0)
            y_data = torch.cat([y_data, torch.tensor([next_value], device=self.device)], dim=0)

        best_index = int(torch.argmin(y_data).item())
        return x_data[best_index], float(y_data[best_index].item()), history

    def _sample_uniform(self, count: int) -> torch.Tensor:
        lows = torch.tensor([spec.bounds.low for spec in self.parameter_specs], device=self.device)
        highs = torch.tensor([spec.bounds.high for spec in self.parameter_specs], device=self.device)
        samples = torch.rand((count, len(self.parameter_specs)), generator=self.generator, device=self.device)
        return lows + samples * (highs - lows)

    def _expected_improvement(self, candidates: torch.Tensor, best_value: torch.Tensor) -> torch.Tensor:
        mean, variance = self.surrogate.predict(candidates)
        std = variance.sqrt().clamp_min(1e-6)
        improvement = best_value - mean - self.config.acquisition_xi
        z = improvement / std
        standard_normal = torch.distributions.Normal(torch.tensor(0.0, device=std.device), torch.tensor(1.0, device=std.device))
        return improvement * standard_normal.cdf(z) + std * torch.exp(standard_normal.log_prob(z))


class CMAESParameterOptimizer:
    """Covariance Matrix Adaptation Evolution Strategy optimizer."""

    def __init__(
        self,
        parameter_specs: Sequence[ParameterSpec],
        config: CMAESConfig,
        *,
        seed: Optional[int] = None,
        device: Optional[torch.device] = None,
    ) -> None:
        self.parameter_specs = list(parameter_specs)
        self.config = config
        self.device = device or torch.device("cpu")
        self.dimension = len(self.parameter_specs)
        self.population_size = config.population_size or (4 + int(3 * math.log(self.dimension)))
        self.mu = self.population_size // 2
        weights = torch.log(torch.tensor(self.mu + 0.5, device=self.device)) - torch.log(
            torch.arange(1, self.mu + 1, device=self.device, dtype=torch.float32)
        )
        self.weights = weights / weights.sum()
        self.mueff = self.weights.sum() ** 2 / (self.weights.square().sum())
        self.cc = (4.0 + self.mueff / self.dimension) / (self.dimension + 4.0 + 2.0 * self.mueff / self.dimension)
        self.cs = (self.mueff + 2.0) / (self.dimension + self.mueff + 5.0)
        self.c1 = 2.0 / (((self.dimension + 1.3) ** 2) + self.mueff)
        self.cmu = min(
            1.0 - self.c1,
            2.0 * (self.mueff - 2.0 + 1.0 / self.mueff) / (((self.dimension + 2.0) ** 2) + self.mueff),
        )
        self.damps = 1.0 + 2.0 * max(0.0, math.sqrt((self.mueff - 1.0) / (self.dimension + 1.0)) - 1.0) + self.cs
        self.chi_n = math.sqrt(self.dimension) * (1.0 - 1.0 / (4.0 * self.dimension) + 1.0 / (21.0 * self.dimension**2))
        self.mean = torch.tensor([spec.default_value for spec in self.parameter_specs], device=self.device, dtype=torch.float32)
        self.sigma = config.sigma
        self.covariance = torch.eye(self.dimension, device=self.device)
        self.path_c = torch.zeros(self.dimension, device=self.device)
        self.path_sigma = torch.zeros(self.dimension, device=self.device)
        self.generator = torch.Generator(device=self.device)
        if seed is not None:
            self.generator.manual_seed(seed)

    def optimize(self, objective: Callable[[torch.Tensor], float]) -> tuple[torch.Tensor, float, list[tuple[torch.Tensor, float]]]:
        history: list[tuple[torch.Tensor, float]] = []
        best_candidate = self.mean.clone()
        best_value = float("inf")

        for generation in range(self.config.generations):
            candidates, directions = self._ask()
            fitness = torch.tensor([float(objective(candidate)) for candidate in candidates], device=self.device)
            history.extend((candidate.clone(), float(value)) for candidate, value in zip(candidates, fitness.tolist()))
            order = torch.argsort(fitness)
            candidates = candidates[order]
            directions = directions[order]
            fitness = fitness[order]
            if float(fitness[0].item()) < best_value:
                best_candidate = candidates[0].clone()
                best_value = float(fitness[0].item())
            self._tell(candidates, directions, generation)
            if self.sigma <= self.config.min_sigma:
                break
        return best_candidate, best_value, history

    def _ask(self) -> tuple[torch.Tensor, torch.Tensor]:
        eigenvalues, eigenvectors = torch.linalg.eigh(self.covariance)
        eigenvalues = eigenvalues.clamp_min(1e-9)
        sqrt_cov = eigenvectors @ torch.diag(eigenvalues.sqrt())
        gaussian = torch.randn((self.population_size, self.dimension), generator=self.generator, device=self.device)
        directions = gaussian @ sqrt_cov.transpose(0, 1)
        candidates = self.mean.unsqueeze(0) + self.sigma * directions
        candidates = self._clip_to_bounds(candidates)
        bounded_directions = (candidates - self.mean.unsqueeze(0)) / max(self.sigma, 1e-9)
        return candidates, bounded_directions

    def _tell(self, candidates: torch.Tensor, directions: torch.Tensor, generation: int) -> None:
        selected_candidates = candidates[: self.mu]
        selected_directions = directions[: self.mu]
        old_mean = self.mean.clone()
        weighted_direction = (self.weights[:, None] * selected_directions).sum(dim=0)
        self.mean = (self.weights[:, None] * selected_candidates).sum(dim=0)

        cov_inv_sqrt = _inverse_sqrt_matrix(self.covariance)
        self.path_sigma = (1.0 - self.cs) * self.path_sigma + math.sqrt(
            self.cs * (2.0 - self.cs) * self.mueff
        ) * (cov_inv_sqrt @ weighted_direction)
        norm_path_sigma = torch.linalg.norm(self.path_sigma)
        h_sigma_condition = norm_path_sigma / math.sqrt(
            1.0 - (1.0 - self.cs) ** (2.0 * (generation + 1))
        ) < (1.4 + 2.0 / (self.dimension + 1.0)) * self.chi_n
        h_sigma = 1.0 if h_sigma_condition else 0.0
        self.path_c = (1.0 - self.cc) * self.path_c + h_sigma * math.sqrt(
            self.cc * (2.0 - self.cc) * self.mueff
        ) * weighted_direction

        rank_one = torch.outer(self.path_c, self.path_c)
        rank_mu = sum(
            weight * torch.outer(direction, direction)
            for weight, direction in zip(self.weights, selected_directions)
        )
        correction = (1.0 - h_sigma) * self.cc * (2.0 - self.cc)
        self.covariance = (
            (1.0 - self.c1 - self.cmu + correction) * self.covariance
            + self.c1 * rank_one
            + self.cmu * rank_mu
        )
        self.covariance = 0.5 * (self.covariance + self.covariance.transpose(0, 1))
        self.sigma *= math.exp((norm_path_sigma / self.chi_n - 1.0) * self.cs / self.damps)
        if torch.isnan(self.mean).any() or torch.isnan(self.covariance).any():
            self.mean = old_mean
            self.covariance = torch.eye(self.dimension, device=self.device)
            self.sigma = self.config.sigma

    def _clip_to_bounds(self, candidates: torch.Tensor) -> torch.Tensor:
        lows = torch.tensor([spec.bounds.low for spec in self.parameter_specs], device=self.device)
        highs = torch.tensor([spec.bounds.high for spec in self.parameter_specs], device=self.device)
        return torch.max(torch.min(candidates, highs), lows)


class SystemIdentifier:
    """End-to-end real-to-sim parameter alignment pipeline."""

    def __init__(
        self,
        simulator: SimulatorProtocol,
        parameter_specs: Sequence[ParameterSpec],
        config: Optional[SystemIdentificationConfig] = None,
        *,
        seed: Optional[int] = None,
        device: Optional[torch.device] = None,
    ) -> None:
        self.simulator = simulator
        self.parameter_specs = list(parameter_specs)
        self.config = config or SystemIdentificationConfig()
        self.seed = seed
        self.device = device or torch.device("cpu")
        self.loss = TrajectoryMatchingLoss(
            position_weight=self.config.position_weight,
            velocity_weight=self.config.velocity_weight,
        )

    def identify(self, real_trajectory: Trajectory, actions: torch.Tensor) -> IdentificationResult:
        def objective(candidate: torch.Tensor) -> float:
            parameters = self._tensor_to_parameters(candidate)
            simulated = self.simulator.rollout(
                parameters,
                actions,
                observation_mask=real_trajectory.observation_mask,
            )
            return float(self.loss(real_trajectory, simulated).item())

        if self.config.optimizer.lower() == "bayesian":
            optimizer = BayesianParameterOptimizer(
                self.parameter_specs,
                self.config.bayesian,
                seed=self.seed,
                device=self.device,
            )
        elif self.config.optimizer.lower() == "cmaes":
            optimizer = CMAESParameterOptimizer(
                self.parameter_specs,
                self.config.cmaes,
                seed=self.seed,
                device=self.device,
            )
        else:
            raise ValueError(f"Unsupported optimizer: {self.config.optimizer}")

        best_tensor, best_loss, history = optimizer.optimize(objective)
        best_parameters = self._tensor_to_parameters(best_tensor)
        best_trajectory = self.simulator.rollout(
            best_parameters,
            actions,
            observation_mask=real_trajectory.observation_mask,
        )
        typed_history = [(self._tensor_to_parameters(candidate), value) for candidate, value in history]
        return IdentificationResult(
            parameters=best_parameters,
            loss=best_loss,
            optimizer=self.config.optimizer,
            history=typed_history,
            simulated_trajectory=best_trajectory,
        )

    def _tensor_to_parameters(self, values: torch.Tensor) -> dict[str, float]:
        return {
            spec.name: float(value.item())
            for spec, value in zip(self.parameter_specs, values.to(dtype=torch.float32), strict=True)
        }


def _ensure_2d(tensor: torch.Tensor) -> torch.Tensor:
    tensor = tensor.to(dtype=torch.float32)
    if tensor.ndim == 1:
        return tensor.unsqueeze(-1)
    if tensor.ndim != 2:
        raise ValueError("Trajectory tensors must be 1-D or 2-D.")
    return tensor


def _align_length(first: torch.Tensor, second: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    steps = min(first.shape[0], second.shape[0])
    return first[:steps], second[:steps]


def _broadcast_mask(mask: Optional[torch.Tensor], shape: torch.Size) -> Optional[torch.Tensor]:
    if mask is None:
        return None
    if mask.shape == shape:
        return mask
    if mask.ndim == 1 and mask.shape[0] == shape[-1]:
        return mask.unsqueeze(0).expand(shape[0], -1)
    if mask.ndim == 1 and mask.shape[0] == shape[0]:
        return mask.unsqueeze(-1).expand(-1, shape[-1])
    raise ValueError("Observation mask shape is incompatible with trajectory shape.")


def _masked_mse(first: torch.Tensor, second: torch.Tensor, mask: Optional[torch.Tensor]) -> torch.Tensor:
    error = (first - second) ** 2
    if mask is None:
        return error.mean()
    weights = mask.to(dtype=error.dtype)
    normalizer = weights.sum().clamp_min(1.0)
    return (error * weights).sum() / normalizer


def _inverse_sqrt_matrix(matrix: torch.Tensor) -> torch.Tensor:
    eigenvalues, eigenvectors = torch.linalg.eigh(matrix)
    inv_sqrt = eigenvectors @ torch.diag(eigenvalues.clamp_min(1e-9).rsqrt()) @ eigenvectors.transpose(0, 1)
    return 0.5 * (inv_sqrt + inv_sqrt.transpose(0, 1))
