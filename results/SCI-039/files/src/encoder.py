"""
encoder.py — Atmospheric variable encoding for pressure-level fields.

Encodes temperature (T), u/v wind components, and specific humidity (q)
across multiple pressure levels with positional (lat/lon/level) embeddings.
"""

import torch
import torch.nn as nn
import numpy as np
from typing import Dict, List, Optional


# Standard pressure levels (hPa) used in ERA5 / GraphCast
PRESSURE_LEVELS = [50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000]
SURFACE_VARIABLES = ['mslp', 't2m', 'u10', 'v10', 'tp']
PRESSURE_VARIABLES = ['t', 'u', 'v', 'q', 'z']  # per level

N_PRESSURE_LEVELS = len(PRESSURE_LEVELS)
N_SURFACE_VARS = len(SURFACE_VARIABLES)
N_PRESSURE_VARS = len(PRESSURE_VARIABLES)
N_INPUT_FEATURES = N_SURFACE_VARS + N_PRESSURE_VARS * N_PRESSURE_LEVELS  # 5 + 5*13 = 70


class SinusoidalPositionEncoding(nn.Module):
    """Sinusoidal encoding for lat/lon/pressure-level coordinates."""

    def __init__(self, d_model: int = 64, max_freq: int = 32):
        super().__init__()
        self.d_model = d_model
        self.max_freq = max_freq
        freqs = 2.0 ** torch.arange(0, max_freq).float()
        self.register_buffer('freqs', freqs)

    def forward(self, coords: torch.Tensor) -> torch.Tensor:
        """
        coords: (N, C) where C = number of coordinate dimensions.
        Returns: (N, C * 2 * max_freq)
        """
        encoded = []
        for i in range(coords.shape[-1]):
            c = coords[:, i:i+1]  # (N, 1)
            c_scaled = c * self.freqs.unsqueeze(0)  # (N, max_freq)
            encoded.append(torch.sin(c_scaled))
            encoded.append(torch.cos(c_scaled))
        return torch.cat(encoded, dim=-1)


class AtmosphericEncoder(nn.Module):
    """
    Encoder for atmospheric state vectors.

    Takes raw atmospheric variables at each grid point and produces
    a latent representation suitable for GNN message passing.

    Input per node: [surface_vars (5), pressure_level_vars (5 × 13)] = 70
    Also incorporates positional encoding from lat/lon/time.
    """

    def __init__(
        self,
        n_input_features: int = N_INPUT_FEATURES,
        d_model: int = 256,
        d_position: int = 64,
        n_layers: int = 2,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.n_input_features = n_input_features
        self.d_model = d_model

        # Positional encoding for (lat, lon, day_of_year, hour)
        self.pos_encoder = SinusoidalPositionEncoding(d_model=d_position)
        pos_dim = 4 * 2 * 32  # 4 coords × 2 (sin+cos) × 32 freqs = 256

        # Variable normalization (learnable per-variable scale/shift)
        self.var_norm = nn.LayerNorm(n_input_features)

        # Main encoder MLP
        layers = []
        in_dim = n_input_features + pos_dim
        for i in range(n_layers):
            out_dim = d_model
            layers.extend([
                nn.Linear(in_dim, out_dim),
                nn.LayerNorm(out_dim),
                nn.SiLU(),
                nn.Dropout(dropout),
            ])
            in_dim = out_dim
        self.mlp = nn.Sequential(*layers)

        # Pressure-level attention (learn importance weighting across levels)
        self.level_attention = nn.Sequential(
            nn.Linear(N_PRESSURE_LEVELS, N_PRESSURE_LEVELS),
            nn.Softmax(dim=-1),
        )

    def forward(
        self,
        x: torch.Tensor,
        lat: torch.Tensor,
        lon: torch.Tensor,
        time_features: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Args:
            x: (N, n_input_features) atmospheric variables
            lat: (N,) latitude in degrees
            lon: (N,) longitude in degrees
            time_features: (N, 2) optional [day_of_year_frac, hour_frac]

        Returns:
            (N, d_model) encoded representation
        """
        # Normalize variables
        x_norm = self.var_norm(x)

        # Positional encoding
        lat_norm = lat / 90.0
        lon_norm = lon / 180.0 - 1.0
        if time_features is None:
            time_features = torch.zeros(x.shape[0], 2, device=x.device)
        coords = torch.stack([lat_norm, lon_norm, time_features[:, 0], time_features[:, 1]], dim=-1)
        pos_enc = self.pos_encoder(coords)

        # Concatenate and encode
        h = torch.cat([x_norm, pos_enc], dim=-1)
        h = self.mlp(h)

        return h


class AtmosphericDecoder(nn.Module):
    """
    Decoder that maps latent GNN features back to atmospheric variables.
    Predicts residuals (delta) from the input state for stable training.
    """

    def __init__(
        self,
        d_model: int = 256,
        n_output_features: int = N_INPUT_FEATURES,
        n_layers: int = 2,
        dropout: float = 0.1,
    ):
        super().__init__()
        layers = []
        in_dim = d_model
        for i in range(n_layers - 1):
            layers.extend([
                nn.Linear(in_dim, d_model),
                nn.LayerNorm(d_model),
                nn.SiLU(),
                nn.Dropout(dropout),
            ])
            in_dim = d_model
        layers.append(nn.Linear(in_dim, n_output_features))
        self.mlp = nn.Sequential(*layers)

    def forward(self, h: torch.Tensor) -> torch.Tensor:
        """
        Args:
            h: (N, d_model) latent features from GNN processor

        Returns:
            (N, n_output_features) predicted residuals
        """
        return self.mlp(h)
