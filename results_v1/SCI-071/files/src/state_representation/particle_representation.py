"""Particle-based state representation for deformable objects."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float32]
IntArray = NDArray[np.int64]


def _as_float_array(value: Any, *, name: str, ndim: int) -> FloatArray:
    array = np.asarray(value, dtype=np.float32)
    if array.ndim != ndim:
        raise ValueError(f"{name} must have {ndim} dimensions, got {array.ndim}.")
    return np.ascontiguousarray(array)


@dataclass(slots=True)
class ParticleState:
    """Structured particle state for a deformable object."""

    positions: FloatArray
    velocities: FloatArray
    masses: FloatArray
    connectivity_radius: float
    neighbors: list[IntArray] = field(default_factory=list)
    neighbor_distances: list[FloatArray] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.positions.ndim != 2 or self.positions.shape[1] != 3:
            raise ValueError(
                f"positions must have shape (N, 3), got {self.positions.shape}."
            )
        if self.velocities.shape != self.positions.shape:
            raise ValueError("velocities must match positions shape.")
        if self.masses.ndim != 1 or self.masses.shape[0] != self.positions.shape[0]:
            raise ValueError("masses must have shape (N,).")
        if np.any(self.masses <= 0.0):
            raise ValueError("masses must be strictly positive.")
        if self.connectivity_radius <= 0.0:
            raise ValueError("connectivity_radius must be positive.")
        if len(self.neighbors) != self.positions.shape[0]:
            raise ValueError("neighbors must contain one entry per particle.")
        if len(self.neighbor_distances) != self.positions.shape[0]:
            raise ValueError("neighbor_distances must contain one entry per particle.")
        for idx, (neighbor_ids, distances) in enumerate(
            zip(self.neighbors, self.neighbor_distances)
        ):
            if neighbor_ids.ndim != 1:
                raise ValueError(f"neighbors[{idx}] must be one-dimensional.")
            if distances.ndim != 1:
                raise ValueError(
                    f"neighbor_distances[{idx}] must be one-dimensional."
                )
            if neighbor_ids.shape[0] != distances.shape[0]:
                raise ValueError(
                    f"neighbors[{idx}] and neighbor_distances[{idx}] must align."
                )

    def __repr__(self) -> str:
        avg_degree = (
            float(np.mean([neighbor.shape[0] for neighbor in self.neighbors]))
            if self.neighbors
            else 0.0
        )
        return (
            "ParticleState("
            f"positions={self.positions.shape}, "
            f"velocities={self.velocities.shape}, "
            f"masses={self.masses.shape}, "
            f"connectivity_radius={self.connectivity_radius:.4f}, "
            f"avg_degree={avg_degree:.2f})"
        )


class ParticleRepresentation:
    """Encodes deformable object observations into a particle state."""

    def __init__(
        self,
        *,
        connectivity_radius: float = 0.05,
        max_neighbors: int | None = None,
        dtype: type[np.float32] = np.float32,
    ) -> None:
        if connectivity_radius <= 0.0:
            raise ValueError("connectivity_radius must be positive.")
        self.connectivity_radius = connectivity_radius
        self.max_neighbors = max_neighbors
        self.dtype = dtype

    def find_neighbors(
        self,
        positions: FloatArray,
        connectivity_radius: float | None = None,
    ) -> tuple[list[IntArray], list[FloatArray]]:
        """Perform radius-based neighbor search."""
        radius = connectivity_radius or self.connectivity_radius
        if radius <= 0.0:
            raise ValueError("connectivity_radius must be positive.")

        deltas = positions[:, None, :] - positions[None, :, :]
        distances = np.linalg.norm(deltas, axis=-1)
        within_radius = (distances <= radius) & (distances > 0.0)

        neighbors: list[IntArray] = []
        neighbor_distances: list[FloatArray] = []
        for idx in range(positions.shape[0]):
            indices = np.nonzero(within_radius[idx])[0]
            dists = distances[idx, indices]
            order = np.argsort(dists)
            if self.max_neighbors is not None:
                order = order[: self.max_neighbors]
            neighbors.append(indices[order].astype(np.int64, copy=False))
            neighbor_distances.append(dists[order].astype(np.float32, copy=False))
        return neighbors, neighbor_distances

    def encode(self, observation: ParticleState | Mapping[str, Any]) -> ParticleState:
        """Encode an observation into :class:`ParticleState`."""
        if isinstance(observation, ParticleState):
            return observation
        if not isinstance(observation, Mapping):
            raise TypeError("Particle observation must be a ParticleState or mapping.")
        if "positions" not in observation:
            raise KeyError("Particle observation must include 'positions'.")

        positions = _as_float_array(observation["positions"], name="positions", ndim=2)
        if positions.shape[1] != 3:
            raise ValueError(f"positions must have shape (N, 3), got {positions.shape}.")
        positions = positions.astype(self.dtype, copy=False)

        velocities_input = observation.get("velocities")
        if velocities_input is None:
            velocities = np.zeros_like(positions, dtype=self.dtype)
        else:
            velocities = _as_float_array(velocities_input, name="velocities", ndim=2)
            if velocities.shape != positions.shape:
                raise ValueError("velocities must match positions shape.")
            velocities = velocities.astype(self.dtype, copy=False)

        masses_input = observation.get("masses")
        if masses_input is None:
            masses = np.ones((positions.shape[0],), dtype=self.dtype)
        else:
            masses = np.asarray(masses_input, dtype=self.dtype)
            if masses.ndim != 1 or masses.shape[0] != positions.shape[0]:
                raise ValueError("masses must have shape (N,).")

        radius = float(observation.get("connectivity_radius", self.connectivity_radius))
        if radius <= 0.0:
            raise ValueError("connectivity_radius must be positive.")

        raw_neighbors = observation.get("neighbors")
        raw_neighbor_distances = observation.get("neighbor_distances")
        if raw_neighbors is None or raw_neighbor_distances is None:
            neighbors, neighbor_distances = self.find_neighbors(positions, radius)
        else:
            neighbors = [np.asarray(ids, dtype=np.int64) for ids in raw_neighbors]
            neighbor_distances = [
                np.asarray(dists, dtype=np.float32) for dists in raw_neighbor_distances
            ]

        metadata = dict(observation.get("metadata", {}))
        return ParticleState(
            positions=positions,
            velocities=velocities,
            masses=masses.astype(self.dtype, copy=False),
            connectivity_radius=radius,
            neighbors=neighbors,
            neighbor_distances=neighbor_distances,
            metadata=metadata,
        )

    def decode(self, state: ParticleState | Mapping[str, Any]) -> dict[str, Any]:
        """Decode a particle state into a dictionary observation."""
        particle_state = self.encode(state)
        return {
            "positions": particle_state.positions.copy(),
            "velocities": particle_state.velocities.copy(),
            "masses": particle_state.masses.copy(),
            "connectivity_radius": particle_state.connectivity_radius,
            "neighbors": [neighbor.copy() for neighbor in particle_state.neighbors],
            "neighbor_distances": [dist.copy() for dist in particle_state.neighbor_distances],
            "metadata": dict(particle_state.metadata),
        }

    def get_feature_dim(self) -> int:
        """Return the per-particle feature dimensionality."""
        return 7

    def __repr__(self) -> str:
        return (
            "ParticleRepresentation("
            f"connectivity_radius={self.connectivity_radius:.4f}, "
            f"max_neighbors={self.max_neighbors}, "
            f"feature_dim={self.get_feature_dim()})"
        )
