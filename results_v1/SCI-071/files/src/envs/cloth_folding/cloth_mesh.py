"""Cloth mesh generation and manipulation utilities."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Tuple

import numpy as np


@dataclass
class MaterialProperties:
    """Default material parameters for cloth simulation."""

    density: float = 0.2
    stretching_stiffness: float = 1.0
    bending_stiffness: float = 0.05
    damping: float = 0.01
    friction: float = 0.6
    thickness: float = 0.002


@dataclass
class ClothMesh:
    """Structured representation of a rectangular cloth mesh."""

    vertices: np.ndarray
    faces: np.ndarray
    edges: np.ndarray
    uvs: np.ndarray
    resolution: Tuple[int, int]
    size: Tuple[float, float]
    material: MaterialProperties = field(default_factory=MaterialProperties)

    def copy(self) -> "ClothMesh":
        """Return a deep copy of the mesh."""

        return ClothMesh(
            vertices=self.vertices.copy(),
            faces=self.faces.copy(),
            edges=self.edges.copy(),
            uvs=self.uvs.copy(),
            resolution=self.resolution,
            size=self.size,
            material=MaterialProperties(**self.material.__dict__),
        )


def generate_uv_coordinates(resolution_x: int, resolution_y: int) -> np.ndarray:
    """Generate UV coordinates for a rectangular grid."""

    u = np.linspace(0.0, 1.0, resolution_x, dtype=np.float32)
    v = np.linspace(0.0, 1.0, resolution_y, dtype=np.float32)
    uu, vv = np.meshgrid(u, v, indexing="xy")
    return np.stack([uu, vv], axis=-1).reshape(-1, 2)


def _grid_faces(resolution_x: int, resolution_y: int) -> np.ndarray:
    faces = []
    for row in range(resolution_y - 1):
        for col in range(resolution_x - 1):
            top_left = row * resolution_x + col
            top_right = top_left + 1
            bottom_left = top_left + resolution_x
            bottom_right = bottom_left + 1
            faces.append([top_left, bottom_left, top_right])
            faces.append([top_right, bottom_left, bottom_right])
    return np.asarray(faces, dtype=np.int32)


def _compute_edges(faces: np.ndarray) -> np.ndarray:
    edges = set()
    for tri in faces:
        pairs = ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0]))
        for start, end in pairs:
            edge = tuple(sorted((int(start), int(end))))
            edges.add(edge)
    return np.asarray(sorted(edges), dtype=np.int32)


def generate_rectangular_cloth_mesh(
    width: float = 0.6,
    height: float = 0.6,
    resolution_x: int = 21,
    resolution_y: int = 21,
    origin: np.ndarray | None = None,
    material: MaterialProperties | None = None,
) -> ClothMesh:
    """Generate a rectangular cloth mesh centered around the origin."""

    if resolution_x < 2 or resolution_y < 2:
        raise ValueError("resolution must be at least 2 in each dimension")
    origin = np.zeros(3, dtype=np.float32) if origin is None else np.asarray(origin, dtype=np.float32)
    x_coords = np.linspace(-width / 2.0, width / 2.0, resolution_x, dtype=np.float32)
    y_coords = np.linspace(-height / 2.0, height / 2.0, resolution_y, dtype=np.float32)
    xx, yy = np.meshgrid(x_coords, y_coords, indexing="xy")
    zz = np.zeros_like(xx, dtype=np.float32)
    vertices = np.stack([xx, yy, zz], axis=-1).reshape(-1, 3) + origin[None, :]
    faces = _grid_faces(resolution_x, resolution_y)
    edges = _compute_edges(faces)
    uvs = generate_uv_coordinates(resolution_x, resolution_y)
    return ClothMesh(
        vertices=vertices,
        faces=faces,
        edges=edges,
        uvs=uvs,
        resolution=(resolution_x, resolution_y),
        size=(width, height),
        material=material or MaterialProperties(),
    )


def subdivide_mesh(mesh: ClothMesh, levels: int = 1) -> ClothMesh:
    """Subdivide a rectangular cloth mesh by regenerating a denser grid."""

    result = mesh.copy()
    for _ in range(levels):
        rx, ry = result.resolution
        result = generate_rectangular_cloth_mesh(
            width=result.size[0],
            height=result.size[1],
            resolution_x=2 * rx - 1,
            resolution_y=2 * ry - 1,
            origin=result.vertices.mean(axis=0),
            material=result.material,
        )
    return result


def simplify_mesh(mesh: ClothMesh, target_face_count: int) -> ClothMesh:
    """Simplify a rectangular cloth mesh while preserving aspect ratio."""

    if target_face_count <= 0:
        raise ValueError("target_face_count must be positive")
    current_face_count = mesh.faces.shape[0]
    if current_face_count <= target_face_count:
        return mesh.copy()
    aspect_ratio = mesh.size[0] / max(mesh.size[1], 1e-6)
    quad_cells = max(target_face_count // 2, 1)
    resolution_y = max(int(np.sqrt(quad_cells / max(aspect_ratio, 1e-6))) + 1, 2)
    resolution_x = max(int(aspect_ratio * (resolution_y - 1)) + 1, 2)
    return generate_rectangular_cloth_mesh(
        width=mesh.size[0],
        height=mesh.size[1],
        resolution_x=resolution_x,
        resolution_y=resolution_y,
        origin=mesh.vertices.mean(axis=0),
        material=mesh.material,
    )


def compute_face_features(mesh: ClothMesh) -> Dict[str, np.ndarray]:
    """Compute face normals, areas, and centers."""

    triangles = mesh.vertices[mesh.faces]
    edge_a = triangles[:, 1] - triangles[:, 0]
    edge_b = triangles[:, 2] - triangles[:, 0]
    normals = np.cross(edge_a, edge_b)
    double_area = np.linalg.norm(normals, axis=1, keepdims=True)
    normals = normals / np.maximum(double_area, 1e-8)
    areas = 0.5 * double_area[:, 0]
    centers = triangles.mean(axis=1)
    return {"normals": normals.astype(np.float32), "areas": areas.astype(np.float32), "centers": centers.astype(np.float32)}


def compute_edge_features(mesh: ClothMesh) -> Dict[str, np.ndarray]:
    """Compute edge lengths and boundary indicators."""

    edge_vertices = mesh.vertices[mesh.edges]
    lengths = np.linalg.norm(edge_vertices[:, 1] - edge_vertices[:, 0], axis=1)
    edge_occurrences = {tuple(edge): 0 for edge in map(tuple, mesh.edges.tolist())}
    for tri in mesh.faces:
        pairs = [tuple(sorted((int(tri[0]), int(tri[1])))), tuple(sorted((int(tri[1]), int(tri[2])))), tuple(sorted((int(tri[2]), int(tri[0]))))]
        for pair in pairs:
            edge_occurrences[pair] += 1
    boundary_mask = np.asarray([edge_occurrences[tuple(edge.tolist())] == 1 for edge in mesh.edges], dtype=bool)
    return {"lengths": lengths.astype(np.float32), "boundary_mask": boundary_mask}


def compute_mesh_features(mesh: ClothMesh) -> Dict[str, np.ndarray]:
    """Compute edge and face features for a cloth mesh."""

    features: Dict[str, np.ndarray] = {}
    features.update({f"face_{key}": value for key, value in compute_face_features(mesh).items()})
    features.update({f"edge_{key}": value for key, value in compute_edge_features(mesh).items()})
    return features


__all__ = [
    "ClothMesh",
    "MaterialProperties",
    "compute_edge_features",
    "compute_face_features",
    "compute_mesh_features",
    "generate_rectangular_cloth_mesh",
    "generate_uv_coordinates",
    "simplify_mesh",
    "subdivide_mesh",
]
