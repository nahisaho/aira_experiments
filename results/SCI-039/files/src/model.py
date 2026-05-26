"""
GraphWeatherNet: A Graph Neural Network for Data-Driven Weather Prediction
Inspired by GraphCast (Lam et al., 2023) and Pangu-Weather (Bi et al., 2023)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import MessagePassing, global_mean_pool
from torch_geometric.data import Data, Batch
import numpy as np
import math


class SphericalGraphBuilder:
    """Build multi-resolution icosahedral mesh graphs on the sphere."""

    def __init__(self, resolutions=[2.5, 1.0, 0.25]):
        self.resolutions = resolutions

    def lat_lon_to_xyz(self, lat, lon):
        lat_r, lon_r = np.radians(lat), np.radians(lon)
        x = np.cos(lat_r) * np.cos(lon_r)
        y = np.cos(lat_r) * np.sin(lon_r)
        z = np.sin(lat_r)
        return np.stack([x, y, z], axis=-1)

    def build_grid(self, resolution):
        lats = np.arange(-90, 90 + resolution, resolution)
        lons = np.arange(0, 360, resolution)
        lat_grid, lon_grid = np.meshgrid(lats, lons, indexing='ij')
        return lat_grid.flatten(), lon_grid.flatten()

    def build_edges(self, lats, lons, resolution, k_neighbors=6):
        n = len(lats)
        xyz = self.lat_lon_to_xyz(lats, lons)
        # Use angular distance for nearest neighbors on sphere
        edges_src, edges_dst = [], []
        # For efficiency, use grid structure
        n_lat = len(np.unique(lats))
        n_lon = len(np.unique(lons))

        for i in range(n):
            lat_idx = i // n_lon
            lon_idx = i % n_lon
            neighbors = []
            for dlat in [-1, 0, 1]:
                for dlon in [-1, 0, 1]:
                    if dlat == 0 and dlon == 0:
                        continue
                    nlat = lat_idx + dlat
                    nlon = (lon_idx + dlon) % n_lon
                    if 0 <= nlat < n_lat:
                        j = nlat * n_lon + nlon
                        neighbors.append(j)
            for j in neighbors:
                edges_src.append(i)
                edges_dst.append(j)

        return np.array(edges_src), np.array(edges_dst)

    def build_graph(self, resolution):
        lats, lons = self.build_grid(resolution)
        src, dst = self.build_edges(lats, lons, resolution)
        pos = self.lat_lon_to_xyz(lats, lons)
        return {
            'lats': lats, 'lons': lons,
            'pos': pos,
            'edge_index': np.stack([src, dst]),
            'n_nodes': len(lats)
        }


class PressureLevelEncoder(nn.Module):
    """Encode atmospheric variables across pressure levels.
    Variables: temperature (T), u-wind, v-wind, specific humidity (q)
    Pressure levels: 50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000 hPa
    """

    def __init__(self, n_pressure_levels=13, n_variables=4, embed_dim=128):
        super().__init__()
        self.n_levels = n_pressure_levels
        self.n_vars = n_variables
        self.embed_dim = embed_dim

        # Learnable pressure level embeddings
        self.level_embedding = nn.Embedding(n_pressure_levels, 32)

        # Variable-specific encoders
        self.var_encoders = nn.ModuleList([
            nn.Sequential(
                nn.Linear(n_pressure_levels, 64),
                nn.GELU(),
                nn.Linear(64, 64)
            ) for _ in range(n_variables)
        ])

        # Surface variable encoder (2m temp, 10m u/v wind, mslp)
        self.surface_encoder = nn.Sequential(
            nn.Linear(4, 64),
            nn.GELU(),
            nn.Linear(64, 64)
        )

        # Fusion layer
        input_dim = 64 * n_variables + 64 + 32 * n_pressure_levels
        self.fusion = nn.Sequential(
            nn.Linear(input_dim, embed_dim * 2),
            nn.GELU(),
            nn.LayerNorm(embed_dim * 2),
            nn.Linear(embed_dim * 2, embed_dim)
        )

    def forward(self, pressure_vars, surface_vars):
        """
        pressure_vars: (B, N, n_vars, n_levels)
        surface_vars: (B, N, 4)
        """
        B, N = pressure_vars.shape[:2]

        # Encode each variable across levels
        var_features = []
        for i, enc in enumerate(self.var_encoders):
            var_features.append(enc(pressure_vars[:, :, i, :]))  # (B, N, 64)

        # Encode surface
        surf_feat = self.surface_encoder(surface_vars)  # (B, N, 64)

        # Pressure level embeddings
        level_idx = torch.arange(self.n_levels, device=pressure_vars.device)
        level_emb = self.level_embedding(level_idx).flatten()  # (n_levels * 32,)
        level_emb = level_emb.unsqueeze(0).unsqueeze(0).expand(B, N, -1)

        # Concatenate all features
        features = torch.cat(var_features + [surf_feat, level_emb], dim=-1)
        return self.fusion(features)  # (B, N, embed_dim)


class GraphProcessorBlock(MessagePassing):
    """Message-passing block with residual connections and layer normalization."""

    def __init__(self, embed_dim=128, edge_dim=3):
        super().__init__(aggr='mean')
        self.embed_dim = embed_dim

        self.edge_mlp = nn.Sequential(
            nn.Linear(edge_dim, 64),
            nn.GELU(),
            nn.Linear(64, embed_dim)
        )

        self.message_mlp = nn.Sequential(
            nn.Linear(embed_dim * 3, embed_dim * 2),
            nn.GELU(),
            nn.Linear(embed_dim * 2, embed_dim)
        )

        self.update_mlp = nn.Sequential(
            nn.Linear(embed_dim * 2, embed_dim * 2),
            nn.GELU(),
            nn.Linear(embed_dim * 2, embed_dim)
        )

        self.norm1 = nn.LayerNorm(embed_dim)
        self.norm2 = nn.LayerNorm(embed_dim)

    def forward(self, x, edge_index, edge_attr):
        edge_feat = self.edge_mlp(edge_attr)
        out = self.propagate(edge_index, x=x, edge_feat=edge_feat)
        x = self.norm1(x + out)
        x = self.norm2(x + self.update_mlp(torch.cat([x, out], dim=-1)))
        return x

    def message(self, x_i, x_j, edge_feat):
        return self.message_mlp(torch.cat([x_i, x_j, edge_feat], dim=-1))


class MultiScaleProcessor(nn.Module):
    """Multi-scale graph processor with cross-resolution connections."""

    def __init__(self, embed_dim=128, n_blocks=4):
        super().__init__()
        self.blocks = nn.ModuleList([
            GraphProcessorBlock(embed_dim) for _ in range(n_blocks)
        ])

    def forward(self, x, edge_index, edge_attr):
        for block in self.blocks:
            x = block(x, edge_index, edge_attr)
        return x


class PressureLevelDecoder(nn.Module):
    """Decode embeddings back to atmospheric variables."""

    def __init__(self, embed_dim=128, n_pressure_levels=13, n_variables=4):
        super().__init__()
        self.n_levels = n_pressure_levels
        self.n_vars = n_variables

        self.var_decoders = nn.ModuleList([
            nn.Sequential(
                nn.Linear(embed_dim, 128),
                nn.GELU(),
                nn.Linear(128, n_pressure_levels)
            ) for _ in range(n_variables)
        ])

        self.surface_decoder = nn.Sequential(
            nn.Linear(embed_dim, 64),
            nn.GELU(),
            nn.Linear(64, 4)
        )

    def forward(self, x):
        """x: (B*N, embed_dim) -> pressure_vars: (B*N, n_vars, n_levels), surface: (B*N, 4)"""
        pressure_vars = torch.stack([dec(x) for dec in self.var_decoders], dim=1)
        surface_vars = self.surface_decoder(x)
        return pressure_vars, surface_vars


class PhysicsConstraintLayer(nn.Module):
    """Soft physics constraints for mass and energy conservation."""

    def __init__(self, n_pressure_levels=13):
        super().__init__()
        self.n_levels = n_pressure_levels
        # Pressure level weights for vertical integration (dp in Pa)
        pressure_levels = torch.tensor(
            [50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000],
            dtype=torch.float32
        ) * 100  # Convert to Pa
        dp = torch.zeros(n_pressure_levels)
        for i in range(n_pressure_levels):
            if i == 0:
                dp[i] = (pressure_levels[1] - pressure_levels[0]) / 2
            elif i == n_pressure_levels - 1:
                dp[i] = (pressure_levels[-1] - pressure_levels[-2]) / 2
            else:
                dp[i] = (pressure_levels[i + 1] - pressure_levels[i - 1]) / 2
        self.register_buffer('dp', dp)
        self.g = 9.80665

    def mass_conservation_loss(self, q_pred, q_input):
        """Ensure column-integrated specific humidity is approximately conserved."""
        col_q_pred = (q_pred * self.dp).sum(dim=-1) / self.g
        col_q_input = (q_input * self.dp).sum(dim=-1) / self.g
        return F.mse_loss(col_q_pred, col_q_input)

    def energy_conservation_loss(self, T_pred, T_input, u_pred, u_input, v_pred, v_input):
        """Approximate total energy conservation."""
        cp = 1004.0  # specific heat at constant pressure
        KE_pred = 0.5 * (u_pred ** 2 + v_pred ** 2)
        KE_input = 0.5 * (u_input ** 2 + v_input ** 2)
        TE_pred = (cp * T_pred + KE_pred) * self.dp
        TE_input = (cp * T_input + KE_input) * self.dp
        col_TE_pred = TE_pred.sum(dim=-1) / self.g
        col_TE_input = TE_input.sum(dim=-1) / self.g
        return F.mse_loss(col_TE_pred, col_TE_input)


class GraphWeatherNet(nn.Module):
    """
    Complete Graph Neural Network for Weather Prediction.
    Architecture: Encoder -> Processor (multi-scale GNN) -> Decoder
    with physics-informed constraints.
    """

    def __init__(self, embed_dim=128, n_pressure_levels=13, n_variables=4, n_blocks=4):
        super().__init__()
        self.encoder = PressureLevelEncoder(n_pressure_levels, n_variables, embed_dim)
        self.processor = MultiScaleProcessor(embed_dim, n_blocks)
        self.decoder = PressureLevelDecoder(embed_dim, n_pressure_levels, n_variables)
        self.physics = PhysicsConstraintLayer(n_pressure_levels)
        self.embed_dim = embed_dim

    def forward(self, pressure_vars, surface_vars, edge_index, edge_attr):
        """
        pressure_vars: (B, N, 4, 13) - T, u, v, q at 13 levels
        surface_vars: (B, N, 4) - 2m T, 10m u, 10m v, MSLP
        """
        B, N = pressure_vars.shape[:2]

        # Normalize inputs for stable training
        p_mean = pressure_vars.mean(dim=(0, 1), keepdim=True)
        p_std = pressure_vars.std(dim=(0, 1), keepdim=True).clamp(min=1e-6)
        pressure_normed = (pressure_vars - p_mean) / p_std

        s_mean = surface_vars.mean(dim=(0, 1), keepdim=True)
        s_std = surface_vars.std(dim=(0, 1), keepdim=True).clamp(min=1e-6)
        surface_normed = (surface_vars - s_mean) / s_std

        # Encode
        x = self.encoder(pressure_normed, surface_normed)  # (B, N, embed_dim)
        x = x.view(B * N, -1)

        # Adjust edge_index for batched graph
        batch_edge_index = []
        batch_edge_attr = []
        for b in range(B):
            batch_edge_index.append(edge_index + b * N)
            batch_edge_attr.append(edge_attr)
        edge_index_batched = torch.cat(batch_edge_index, dim=1)
        edge_attr_batched = torch.cat(batch_edge_attr, dim=0)

        # Process
        x = self.processor(x, edge_index_batched, edge_attr_batched)

        # Decode
        pred_pressure, pred_surface = self.decoder(x)
        pred_pressure = pred_pressure.view(B, N, 4, -1)
        pred_surface = pred_surface.view(B, N, 4)

        # Denormalize: output as residual added to input
        pred_pressure = pressure_vars + pred_pressure * p_std.squeeze(0)
        pred_surface = surface_vars + pred_surface * s_std.squeeze(0)

        return pred_pressure, pred_surface

    def compute_loss(self, pred_p, pred_s, target_p, target_s, input_p,
                     physics_weight=0.1):
        """Combined MSE + physics constraint loss."""
        mse_p = F.mse_loss(pred_p, target_p)
        mse_s = F.mse_loss(pred_s, target_s)

        # Physics constraints
        mass_loss = self.physics.mass_conservation_loss(
            pred_p[:, :, 3, :], input_p[:, :, 3, :])
        energy_loss = self.physics.energy_conservation_loss(
            pred_p[:, :, 0, :], input_p[:, :, 0, :],
            pred_p[:, :, 1, :], input_p[:, :, 1, :],
            pred_p[:, :, 2, :], input_p[:, :, 2, :])

        total = mse_p + mse_s + physics_weight * (mass_loss + energy_loss)
        return total, {
            'mse_pressure': mse_p.item(),
            'mse_surface': mse_s.item(),
            'mass_conservation': mass_loss.item(),
            'energy_conservation': energy_loss.item(),
            'total': total.item()
        }
