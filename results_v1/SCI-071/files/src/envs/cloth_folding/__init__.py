"""Cloth folding environment package."""

from .cloth_folding_env import ClothFoldingEnv
from .cloth_mesh import (
    ClothMesh,
    MaterialProperties,
    compute_edge_features,
    compute_face_features,
    compute_mesh_features,
    generate_rectangular_cloth_mesh,
    generate_uv_coordinates,
    simplify_mesh,
    subdivide_mesh,
)
from .fold_primitives import (
    flatten,
    fold_edge,
    interpolate_trajectory,
    pick_and_place,
    plan_collision_aware_trajectory,
    shake,
)
from .reward_functions import (
    chamfer_distance_reward,
    composite_reward,
    coverage_reward,
    fold_line_reward,
    smoothness_penalty,
)

__all__ = [
    "ClothFoldingEnv",
    "ClothMesh",
    "MaterialProperties",
    "compute_edge_features",
    "compute_face_features",
    "compute_mesh_features",
    "generate_rectangular_cloth_mesh",
    "generate_uv_coordinates",
    "simplify_mesh",
    "subdivide_mesh",
    "flatten",
    "fold_edge",
    "interpolate_trajectory",
    "pick_and_place",
    "plan_collision_aware_trajectory",
    "shake",
    "chamfer_distance_reward",
    "composite_reward",
    "coverage_reward",
    "fold_line_reward",
    "smoothness_penalty",
]
