"""Mesh-based state representation for deformable objects."""

from __future__ import annotations

from dataclasses import dataclass
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


def _as_int_array(value: Any, *, name: str, ndim: int) -> IntArray:
    array = np.asarray(value, dtype=np.int64)
    if array.ndim != ndim:
        raise ValueError(f"{name} must have {ndim} dimensions, got {array.ndim}.")
    return np.ascontiguousarray(array)


def _validate_vertices(vertices: FloatArray) -> None:
    if vertices.shape[1] != 3:
        raise ValueError(f"vertices must have shape (N, 3), got {vertices.shape}.")


def _validate_faces(faces: IntArray, num_vertices: int) -> None:
    if faces.size == 0:
        return
    if faces.shape[1] != 3:
        raise ValueError(f"faces must have shape (M, 3), got {faces.shape}.")
    if np.any(faces < 0) or np.any(faces >= num_vertices):
        raise ValueError("faces contain invalid vertex indices.")


def _validate_vertex_features(vertex_features: FloatArray, num_vertices: int) -> None:
    if vertex_features.shape[0] != num_vertices:
        raise ValueError(
            "vertex_features must have the same number of rows as vertices."
        )


def _validate_edge_index(edge_index: IntArray, num_vertices: int) -> None:
    if edge_index.size == 0:
        return
    if edge_index.shape[1] != 2:
        raise ValueError(f"edge_index must have shape (E, 2), got {edge_index.shape}.")
    if np.any(edge_index < 0) or np.any(edge_index >= num_vertices):
        raise ValueError("edge_index contains invalid vertex indices.")


def _validate_edge_features(edge_features: FloatArray, num_edges: int) -> None:
    if edge_features.shape[0] != num_edges:
        raise ValueError("edge_features must align with edge_index.")


def _compute_unique_edges(faces: IntArray) -> IntArray:
    if faces.size == 0:
        return np.empty((0, 2), dtype=np.int64)

    undirected_edges = np.concatenate(
        [
            faces[:, [0, 1]],
            faces[:, [1, 2]],
            faces[:, [2, 0]],
        ],
        axis=0,
    )
    undirected_edges = np.sort(undirected_edges, axis=1)
    return np.unique(undirected_edges, axis=0)


def _compute_vertex_normals(vertices: FloatArray, faces: IntArray) -> FloatArray:
    normals = np.zeros_like(vertices, dtype=np.float32)
    if faces.size == 0:
        return normals

    tri_vertices = vertices[faces]
    edge_ab = tri_vertices[:, 1] - tri_vertices[:, 0]
    edge_ac = tri_vertices[:, 2] - tri_vertices[:, 0]
    face_normals = np.cross(edge_ab, edge_ac)
    face_magnitudes = np.linalg.norm(face_normals, axis=1, keepdims=True)
    valid_faces = face_magnitudes.squeeze(-1) > 1e-12
    if np.any(valid_faces):
        face_normals[valid_faces] /= face_magnitudes[valid_faces]
    face_normals[~valid_faces] = 0.0

    for corner in range(3):
        np.add.at(normals, faces[:, corner], face_normals)

    magnitudes = np.linalg.norm(normals, axis=1, keepdims=True)
    valid_vertices = magnitudes.squeeze(-1) > 1e-12
    normals[valid_vertices] /= magnitudes[valid_vertices]
    normals[~valid_vertices] = 0.0
    return normals.astype(np.float32, copy=False)


def _compute_edge_features(vertices: FloatArray, edge_index: IntArray) -> FloatArray:
    if edge_index.size == 0:
        return np.empty((0, 4), dtype=np.float32)

    edge_vectors = vertices[edge_index[:, 1]] - vertices[edge_index[:, 0]]
    lengths = np.linalg.norm(edge_vectors, axis=1, keepdims=True)
    directions = np.divide(
        edge_vectors,
        np.clip(lengths, a_min=1e-12, a_max=None),
        out=np.zeros_like(edge_vectors),
    )
    return np.concatenate([lengths.astype(np.float32), directions.astype(np.float32)], axis=1)


@dataclass(slots=True)
class MeshState:
    """Structured mesh state for a deformable object."""

    vertices: FloatArray
    faces: IntArray
    vertex_features: FloatArray
    edge_index: IntArray
    edge_features: FloatArray

    def __post_init__(self) -> None:
        _validate_vertices(self.vertices)
        _validate_faces(self.faces, self.vertices.shape[0])
        _validate_vertex_features(self.vertex_features, self.vertices.shape[0])
        _validate_edge_index(self.edge_index, self.vertices.shape[0])
        _validate_edge_features(self.edge_features, self.edge_index.shape[0])

    def __repr__(self) -> str:
        return (
            "MeshState("
            f"vertices={self.vertices.shape}, "
            f"faces={self.faces.shape}, "
            f"vertex_features={self.vertex_features.shape}, "
            f"edge_index={self.edge_index.shape}, "
            f"edge_features={self.edge_features.shape})"
        )


class MeshRepresentation:
    """Encodes deformable object observations into a mesh state."""

    def __init__(
        self,
        *,
        include_normals: bool = True,
        include_velocities: bool = True,
        dtype: type[np.float32] = np.float32,
    ) -> None:
        self.include_normals = include_normals
        self.include_velocities = include_velocities
        self.dtype = dtype
        self._feature_dim = (3 if include_normals else 0) + (3 if include_velocities else 0)

    def encode(self, observation: MeshState | Mapping[str, Any]) -> MeshState:
        """Encode a mesh observation into :class:`MeshState`."""
        if isinstance(observation, MeshState):
            self._feature_dim = int(observation.vertex_features.shape[1])
            return observation
        if not isinstance(observation, Mapping):
            raise TypeError("Mesh observation must be a MeshState or mapping.")

        if "vertices" not in observation:
            raise KeyError("Mesh observation must include 'vertices'.")
        vertices = _as_float_array(observation["vertices"], name="vertices", ndim=2)
        _validate_vertices(vertices)
        vertices = vertices.astype(self.dtype, copy=False)

        raw_faces = observation.get("faces", np.empty((0, 3), dtype=np.int64))
        faces = _as_int_array(raw_faces, name="faces", ndim=2)
        _validate_faces(faces, vertices.shape[0])

        feature_blocks: list[FloatArray] = []
        if self.include_normals:
            normals = observation.get("normals")
            if normals is None:
                normals_array = _compute_vertex_normals(vertices, faces)
            else:
                normals_array = _as_float_array(normals, name="normals", ndim=2)
                _validate_vertex_features(normals_array, vertices.shape[0])
                if normals_array.shape[1] != 3:
                    raise ValueError("normals must have shape (N, 3).")
            feature_blocks.append(normals_array.astype(self.dtype, copy=False))

        if self.include_velocities:
            velocities = observation.get("velocities")
            if velocities is None:
                velocities_array = np.zeros_like(vertices, dtype=self.dtype)
            else:
                velocities_array = _as_float_array(velocities, name="velocities", ndim=2)
                _validate_vertex_features(velocities_array, vertices.shape[0])
                if velocities_array.shape[1] != 3:
                    raise ValueError("velocities must have shape (N, 3).")
                velocities_array = velocities_array.astype(self.dtype, copy=False)
            feature_blocks.append(velocities_array)

        extra_vertex_features = observation.get("vertex_features")
        if extra_vertex_features is not None:
            extra_features_array = _as_float_array(
                extra_vertex_features,
                name="vertex_features",
                ndim=2,
            )
            _validate_vertex_features(extra_features_array, vertices.shape[0])
            feature_blocks.append(extra_features_array.astype(self.dtype, copy=False))

        vertex_features = (
            np.concatenate(feature_blocks, axis=1)
            if feature_blocks
            else np.empty((vertices.shape[0], 0), dtype=self.dtype)
        )

        edge_index_input = observation.get("edge_index")
        if edge_index_input is None:
            edge_index = _compute_unique_edges(faces)
        else:
            edge_index = _as_int_array(edge_index_input, name="edge_index", ndim=2)
        _validate_edge_index(edge_index, vertices.shape[0])

        edge_features_input = observation.get("edge_features")
        if edge_features_input is None:
            edge_features = _compute_edge_features(vertices, edge_index)
        else:
            edge_features = _as_float_array(
                edge_features_input,
                name="edge_features",
                ndim=2,
            ).astype(self.dtype, copy=False)
        _validate_edge_features(edge_features, edge_index.shape[0])

        state = MeshState(
            vertices=vertices,
            faces=faces,
            vertex_features=vertex_features,
            edge_index=edge_index,
            edge_features=edge_features,
        )
        self._feature_dim = int(state.vertex_features.shape[1])
        return state

    def decode(self, state: MeshState | Mapping[str, Any]) -> dict[str, FloatArray | IntArray]:
        """Decode a mesh state into a dictionary observation."""
        mesh_state = self.encode(state)
        observation: dict[str, FloatArray | IntArray] = {
            "vertices": mesh_state.vertices.copy(),
            "faces": mesh_state.faces.copy(),
            "vertex_features": mesh_state.vertex_features.copy(),
            "edge_index": mesh_state.edge_index.copy(),
            "edge_features": mesh_state.edge_features.copy(),
        }
        offset = 0
        if self.include_normals and mesh_state.vertex_features.shape[1] >= offset + 3:
            observation["normals"] = mesh_state.vertex_features[:, offset : offset + 3].copy()
            offset += 3
        if self.include_velocities and mesh_state.vertex_features.shape[1] >= offset + 3:
            observation["velocities"] = mesh_state.vertex_features[:, offset : offset + 3].copy()
        return observation

    def get_feature_dim(self) -> int:
        """Return the mesh vertex feature dimensionality."""
        return self._feature_dim

    def __repr__(self) -> str:
        return (
            "MeshRepresentation("
            f"include_normals={self.include_normals}, "
            f"include_velocities={self.include_velocities}, "
            f"feature_dim={self._feature_dim})"
        )
