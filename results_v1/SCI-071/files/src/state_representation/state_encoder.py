"""Unified state encoder for deformable object representations."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Literal, cast

import numpy as np
import torch

from .latent_representation import LatentRepresentation, LatentState
from .mesh_representation import MeshRepresentation, MeshState
from .particle_representation import ParticleRepresentation, ParticleState


RepresentationName = Literal["mesh", "particle", "latent"]


class StateEncoder:
    """Convert deformable object states between mesh, particle, and latent spaces."""

    def __init__(
        self,
        *,
        mesh_representation: MeshRepresentation | None = None,
        particle_representation: ParticleRepresentation | None = None,
        latent_representation: LatentRepresentation | None = None,
        default_representation: RepresentationName = "latent",
    ) -> None:
        if default_representation not in {"mesh", "particle", "latent"}:
            raise ValueError("default_representation must be one of 'mesh', 'particle', or 'latent'.")
        self.mesh_representation = mesh_representation or MeshRepresentation()
        self.particle_representation = particle_representation or ParticleRepresentation()
        self.latent_representation = latent_representation or LatentRepresentation()
        self.default_representation = default_representation

    def encode(self, observation: Any) -> MeshState | ParticleState | LatentState | list[Any]:
        """Encode an observation into the configured default representation."""
        if isinstance(observation, (MeshState, ParticleState, LatentState)):
            source = self._infer_representation(observation)
            return self.convert(observation, source=source, target=self.default_representation)
        return self._encode_to(self.default_representation, observation)

    def decode(self, state: Any) -> Any:
        """Decode a state back into an observation."""
        if isinstance(state, list):
            return [self.decode(item) for item in state]
        representation = self._infer_representation(state)
        return self._representation_for(representation).decode(state)

    def convert(
        self,
        data: Any,
        *,
        source: RepresentationName,
        target: RepresentationName,
    ) -> MeshState | ParticleState | LatentState | list[Any]:
        """Convert data between any pair of supported representations."""
        if source not in {"mesh", "particle", "latent"}:
            raise ValueError(f"Unsupported source representation: {source}.")
        if target not in {"mesh", "particle", "latent"}:
            raise ValueError(f"Unsupported target representation: {target}.")
        if source == target:
            return self._encode_to(source, data)

        batched = self._split_batch(data, source)
        if batched is not None:
            return [self.convert(item, source=source, target=target) for item in batched]

        source_state = self._encode_to(source, data)
        return self._convert_state(source_state, source=source, target=target)

    def get_feature_dim(self) -> int:
        """Return the feature dimension of the default representation."""
        return self._representation_for(self.default_representation).get_feature_dim()

    def _convert_state(
        self,
        state: MeshState | ParticleState | LatentState,
        *,
        source: RepresentationName,
        target: RepresentationName,
    ) -> MeshState | ParticleState | LatentState:
        if source == "mesh" and target == "particle":
            return self._mesh_to_particle(cast(MeshState, state))
        if source == "particle" and target == "mesh":
            return self._particle_to_mesh(cast(ParticleState, state))
        if source == "mesh" and target == "latent":
            return self.latent_representation.encode(cast(MeshState, state).vertices)
        if source == "particle" and target == "latent":
            return self.latent_representation.encode(cast(ParticleState, state).positions)
        if source == "latent" and target == "particle":
            return self._latent_to_particle(cast(LatentState, state))
        if source == "latent" and target == "mesh":
            particle_state = self._latent_to_particle(cast(LatentState, state))
            return self._particle_to_mesh(particle_state)
        raise RuntimeError(f"Unsupported conversion: {source} -> {target}.")

    def _encode_to(
        self,
        representation: RepresentationName,
        data: Any,
    ) -> MeshState | ParticleState | LatentState:
        encoder = self._representation_for(representation)
        return encoder.encode(data)

    def _representation_for(self, name: RepresentationName) -> Any:
        if name == "mesh":
            return self.mesh_representation
        if name == "particle":
            return self.particle_representation
        if name == "latent":
            return self.latent_representation
        raise ValueError(f"Unknown representation: {name}.")

    def _infer_representation(self, state: Any) -> RepresentationName:
        if isinstance(state, MeshState):
            return "mesh"
        if isinstance(state, ParticleState):
            return "particle"
        if isinstance(state, LatentState):
            return "latent"
        raise TypeError(
            "Unable to infer representation type. Expected MeshState, ParticleState, or LatentState."
        )

    def _split_batch(self, data: Any, source: RepresentationName) -> list[Any] | None:
        if source == "latent":
            if isinstance(data, LatentState) and data.latent.shape[0] > 1:
                return [self._slice_latent_state(data, index) for index in range(data.latent.shape[0])]
            return None

        if isinstance(data, Sequence) and not isinstance(data, (str, bytes, bytearray, np.ndarray, torch.Tensor)):
            return list(data)

        if not isinstance(data, Mapping):
            return None

        if source == "mesh" and "vertices" in data:
            vertices = np.asarray(data["vertices"])
            if vertices.ndim == 3:
                batch_size = vertices.shape[0]
                return [self._slice_mapping_batch(data, batch_size, index) for index in range(batch_size)]
        if source == "particle" and "positions" in data:
            positions = np.asarray(data["positions"])
            if positions.ndim == 3:
                batch_size = positions.shape[0]
                return [self._slice_mapping_batch(data, batch_size, index) for index in range(batch_size)]
        return None

    def _slice_mapping_batch(
        self,
        data: Mapping[str, Any],
        batch_size: int,
        index: int,
    ) -> dict[str, Any]:
        sliced: dict[str, Any] = {}
        for key, value in data.items():
            if isinstance(value, np.ndarray) and value.ndim > 0 and value.shape[0] == batch_size:
                sliced[key] = value[index]
            elif isinstance(value, torch.Tensor) and value.ndim > 0 and value.shape[0] == batch_size:
                sliced[key] = value[index]
            elif isinstance(value, list) and len(value) == batch_size:
                sliced[key] = value[index]
            else:
                sliced[key] = value
        return sliced

    def _slice_latent_state(self, state: LatentState, index: int) -> LatentState:
        return LatentState(
            latent=state.latent[index : index + 1],
            mean=state.mean[index : index + 1],
            logvar=state.logvar[index : index + 1],
            num_points=(state.num_points[index],),
            input_dim=state.input_dim,
            batched=False,
            mask=state.mask[index : index + 1] if state.mask is not None else None,
        )

    def _mesh_to_particle(self, mesh_state: MeshState) -> ParticleState:
        velocities = self._extract_mesh_velocities(mesh_state)
        masses = self._estimate_vertex_masses(mesh_state)
        radius = self._infer_mesh_radius(mesh_state)
        metadata = {
            "faces": mesh_state.faces.copy(),
            "vertex_features": mesh_state.vertex_features.copy(),
            "edge_index": mesh_state.edge_index.copy(),
            "edge_features": mesh_state.edge_features.copy(),
        }
        return self.particle_representation.encode(
            {
                "positions": mesh_state.vertices,
                "velocities": velocities,
                "masses": masses,
                "connectivity_radius": radius,
                "metadata": metadata,
            }
        )

    def _particle_to_mesh(self, particle_state: ParticleState) -> MeshState:
        metadata = particle_state.metadata
        faces = np.asarray(
            metadata.get("faces", self._generate_faces_from_particles(particle_state)),
            dtype=np.int64,
        )
        extra_vertex_features = self._extract_extra_vertex_features(metadata)
        mesh_observation: dict[str, Any] = {
            "vertices": particle_state.positions,
            "faces": faces,
            "velocities": particle_state.velocities,
        }
        if extra_vertex_features is not None:
            mesh_observation["vertex_features"] = extra_vertex_features
        if "edge_index" in metadata:
            mesh_observation["edge_index"] = np.asarray(metadata["edge_index"], dtype=np.int64)
        if "edge_features" in metadata:
            mesh_observation["edge_features"] = np.asarray(
                metadata["edge_features"],
                dtype=np.float32,
            )
        return self.mesh_representation.encode(mesh_observation)

    def _latent_to_particle(self, latent_state: LatentState) -> ParticleState:
        decoded = self.latent_representation.decode(latent_state)["points"]
        if isinstance(decoded, list):
            if len(decoded) != 1:
                raise ValueError("Batched latent states must be split before latent->particle conversion.")
            points_tensor = decoded[0]
        else:
            points_tensor = decoded if decoded.ndim == 2 else decoded.squeeze(0)
        positions = points_tensor.detach().cpu().numpy().astype(np.float32, copy=False)
        connectivity_radius = self._infer_particle_radius(positions)
        return self.particle_representation.encode(
            {
                "positions": positions,
                "velocities": np.zeros_like(positions, dtype=np.float32),
                "masses": np.ones((positions.shape[0],), dtype=np.float32),
                "connectivity_radius": connectivity_radius,
            }
        )

    def _extract_mesh_velocities(self, mesh_state: MeshState) -> np.ndarray:
        vertex_features = mesh_state.vertex_features
        if vertex_features.shape[1] >= 6:
            return vertex_features[:, 3:6].copy()
        return np.zeros_like(mesh_state.vertices, dtype=np.float32)

    def _extract_extra_vertex_features(self, metadata: Mapping[str, Any]) -> np.ndarray | None:
        vertex_features = metadata.get("vertex_features")
        if vertex_features is None:
            return None
        features = np.asarray(vertex_features, dtype=np.float32)
        if features.ndim != 2:
            raise ValueError("Stored vertex_features metadata must be two-dimensional.")
        if features.shape[1] <= 6:
            return None
        return features[:, 6:]

    def _estimate_vertex_masses(self, mesh_state: MeshState) -> np.ndarray:
        masses = np.ones((mesh_state.vertices.shape[0],), dtype=np.float32)
        if mesh_state.faces.size == 0:
            return masses

        triangle_vertices = mesh_state.vertices[mesh_state.faces]
        edge_ab = triangle_vertices[:, 1] - triangle_vertices[:, 0]
        edge_ac = triangle_vertices[:, 2] - triangle_vertices[:, 0]
        areas = 0.5 * np.linalg.norm(np.cross(edge_ab, edge_ac), axis=1)
        masses = np.zeros((mesh_state.vertices.shape[0],), dtype=np.float32)
        for corner in range(3):
            np.add.at(masses, mesh_state.faces[:, corner], areas / 3.0)
        masses[masses <= 1e-8] = 1.0
        return masses

    def _infer_mesh_radius(self, mesh_state: MeshState) -> float:
        if mesh_state.edge_index.size == 0:
            return self.particle_representation.connectivity_radius
        edge_vectors = (
            mesh_state.vertices[mesh_state.edge_index[:, 1]]
            - mesh_state.vertices[mesh_state.edge_index[:, 0]]
        )
        edge_lengths = np.linalg.norm(edge_vectors, axis=1)
        median_length = float(np.median(edge_lengths)) if edge_lengths.size else 0.0
        if median_length <= 1e-8:
            return self.particle_representation.connectivity_radius
        return max(median_length * 1.5, 1e-6)

    def _infer_particle_radius(self, positions: np.ndarray) -> float:
        if positions.shape[0] < 2:
            return self.particle_representation.connectivity_radius
        deltas = positions[:, None, :] - positions[None, :, :]
        distances = np.linalg.norm(deltas, axis=-1)
        np.fill_diagonal(distances, np.inf)
        nearest = np.min(distances, axis=1)
        finite = nearest[np.isfinite(nearest)]
        if finite.size == 0:
            return self.particle_representation.connectivity_radius
        return max(float(np.median(finite)) * 1.5, 1e-6)

    def _generate_faces_from_particles(self, particle_state: ParticleState) -> np.ndarray:
        triangles: set[tuple[int, int, int]] = set()
        for i, neighbors_i in enumerate(particle_state.neighbors):
            if neighbors_i.shape[0] < 2:
                continue
            limited_neighbors = neighbors_i[: min(neighbors_i.shape[0], 12)]
            neighbor_set_i = set(int(value) for value in limited_neighbors.tolist())
            for j in sorted(neighbor_set_i):
                if j <= i:
                    continue
                neighbor_set_j = set(int(value) for value in particle_state.neighbors[j].tolist())
                common = sorted(value for value in neighbor_set_i.intersection(neighbor_set_j) if value > j)
                for k in common[:12]:
                    triangle = (i, j, k)
                    area = self._triangle_area(particle_state.positions[list(triangle)])
                    if area > 1e-10:
                        triangles.add(triangle)
        if not triangles:
            return np.empty((0, 3), dtype=np.int64)
        return np.asarray(sorted(triangles), dtype=np.int64)

    def _triangle_area(self, triangle_vertices: np.ndarray) -> float:
        edge_ab = triangle_vertices[1] - triangle_vertices[0]
        edge_ac = triangle_vertices[2] - triangle_vertices[0]
        return 0.5 * float(np.linalg.norm(np.cross(edge_ab, edge_ac)))

    def __repr__(self) -> str:
        return (
            "StateEncoder("
            f"default_representation={self.default_representation!r}, "
            f"mesh={self.mesh_representation!r}, "
            f"particle={self.particle_representation!r}, "
            f"latent={self.latent_representation!r})"
        )
