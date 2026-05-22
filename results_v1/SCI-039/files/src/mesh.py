"""
mesh.py — Multi-scale icosahedral mesh construction for GraphCast-style models.

Builds hierarchical meshes at 3 resolutions (0.25°, 1°, 2.5°) and creates
grid-to-mesh / mesh-to-grid bipartite edges plus intra-mesh edges.
"""

import numpy as np
import torch
from torch_geometric.data import Data
from typing import Tuple, Dict, List


def lat_lon_to_xyz(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Convert lat/lon (degrees) to 3D Cartesian coordinates on unit sphere."""
    lat_r = np.deg2rad(lat)
    lon_r = np.deg2rad(lon)
    x = np.cos(lat_r) * np.cos(lon_r)
    y = np.cos(lat_r) * np.sin(lon_r)
    z = np.sin(lat_r)
    return np.stack([x, y, z], axis=-1)


def create_regular_grid(resolution: float) -> Tuple[np.ndarray, np.ndarray]:
    """Create a regular lat-lon grid at given resolution (degrees)."""
    lats = np.arange(-90, 90 + resolution, resolution)
    lons = np.arange(0, 360, resolution)
    lat_grid, lon_grid = np.meshgrid(lats, lons, indexing='ij')
    return lat_grid.flatten(), lon_grid.flatten()


def build_knn_edges(src_xyz: np.ndarray, dst_xyz: np.ndarray, k: int = 6) -> np.ndarray:
    """Build k-nearest-neighbor edges from src to dst nodes."""
    from scipy.spatial import cKDTree
    tree = cKDTree(dst_xyz)
    _, indices = tree.query(src_xyz, k=k)
    src_indices = np.repeat(np.arange(len(src_xyz)), k)
    dst_indices = indices.flatten()
    return np.stack([src_indices, dst_indices], axis=0)


def build_radius_edges(xyz: np.ndarray, radius: float) -> np.ndarray:
    """Build edges within a radius on the unit sphere."""
    from scipy.spatial import cKDTree
    tree = cKDTree(xyz)
    pairs = tree.query_pairs(r=radius, output_type='ndarray')
    if len(pairs) == 0:
        return np.zeros((2, 0), dtype=np.int64)
    edges = np.concatenate([pairs, pairs[:, ::-1]], axis=0)
    return edges.T


class MultiScaleMesh:
    """
    Multi-scale mesh for weather prediction.

    Resolutions:
      - Fine:   0.25° (721 x 1440 = ~1M nodes, subsampled for demo)
      - Medium: 1.0°  (181 x 360 = ~65K nodes)
      - Coarse: 2.5°  (73 x 144 = ~10K nodes)
    """

    RESOLUTIONS = {
        'fine': 0.25,
        'medium': 1.0,
        'coarse': 2.5,
    }

    def __init__(self, use_subset: bool = True, subset_factor: int = 4):
        self.use_subset = use_subset
        self.subset_factor = subset_factor
        self.grids = {}
        self.xyz = {}
        self.n_nodes = {}

        for name, res in self.RESOLUTIONS.items():
            effective_res = res * subset_factor if use_subset else res
            lat, lon = create_regular_grid(effective_res)
            self.grids[name] = (lat, lon)
            self.xyz[name] = lat_lon_to_xyz(lat, lon)
            self.n_nodes[name] = len(lat)

    def build_mesh_graph(self) -> Dict:
        """Build full multi-scale graph with inter-scale and intra-scale edges."""
        # Intra-scale edges (within each resolution)
        intra_edges = {}
        radii = {'fine': 0.05, 'medium': 0.06, 'coarse': 0.15}
        for name in self.RESOLUTIONS:
            intra_edges[name] = build_radius_edges(self.xyz[name], radii[name])

        # Inter-scale edges (between resolutions)
        inter_edges = {}
        inter_edges['fine_to_medium'] = build_knn_edges(
            self.xyz['fine'], self.xyz['medium'], k=4
        )
        inter_edges['medium_to_coarse'] = build_knn_edges(
            self.xyz['medium'], self.xyz['coarse'], k=4
        )
        inter_edges['coarse_to_medium'] = build_knn_edges(
            self.xyz['coarse'], self.xyz['medium'], k=4
        )
        inter_edges['medium_to_fine'] = build_knn_edges(
            self.xyz['medium'], self.xyz['fine'], k=4
        )

        return {
            'n_nodes': self.n_nodes,
            'grids': self.grids,
            'xyz': self.xyz,
            'intra_edges': intra_edges,
            'inter_edges': inter_edges,
        }

    def get_grid_info(self) -> str:
        """Return summary string of mesh configuration."""
        lines = ["Multi-Scale Mesh Configuration:"]
        for name, res in self.RESOLUTIONS.items():
            eff = res * self.subset_factor if self.use_subset else res
            lines.append(f"  {name:8s}: {eff:5.2f}° resolution, {self.n_nodes[name]:>6d} nodes")
        return "\n".join(lines)


def create_pyg_data(mesh_graph: Dict, features: Dict[str, torch.Tensor]) -> Data:
    """
    Convert multi-scale mesh graph to PyG Data object.

    Features dict maps scale name -> tensor of shape (n_nodes, n_features).
    Nodes are concatenated: [fine_nodes, medium_nodes, coarse_nodes].
    """
    offsets = {}
    offset = 0
    all_features = []
    for name in ['fine', 'medium', 'coarse']:
        offsets[name] = offset
        offset += mesh_graph['n_nodes'][name]
        all_features.append(features[name])

    x = torch.cat(all_features, dim=0)

    # Collect all edges with proper offsets
    all_edges = []

    # Intra-scale edges
    for name in ['fine', 'medium', 'coarse']:
        edges = mesh_graph['intra_edges'][name]
        if edges.shape[1] > 0:
            edges_t = torch.tensor(edges, dtype=torch.long) + offsets[name]
            all_edges.append(edges_t)

    # Inter-scale edges
    scale_map = {
        'fine_to_medium': ('fine', 'medium'),
        'medium_to_coarse': ('medium', 'coarse'),
        'coarse_to_medium': ('coarse', 'medium'),
        'medium_to_fine': ('medium', 'fine'),
    }
    for key, (src_scale, dst_scale) in scale_map.items():
        edges = mesh_graph['inter_edges'][key]
        src_edges = torch.tensor(edges[0], dtype=torch.long) + offsets[src_scale]
        dst_edges = torch.tensor(edges[1], dtype=torch.long) + offsets[dst_scale]
        all_edges.append(torch.stack([src_edges, dst_edges], dim=0))

    edge_index = torch.cat(all_edges, dim=1) if all_edges else torch.zeros(2, 0, dtype=torch.long)

    return Data(x=x, edge_index=edge_index)
