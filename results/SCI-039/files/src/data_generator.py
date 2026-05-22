"""
data_generator.py — Synthetic ERA5-like data generator for model development.

Generates realistic atmospheric fields with proper spatial correlations,
vertical structure, and temporal evolution for training and evaluation.
"""

import numpy as np
import torch
from typing import Tuple, Dict, Optional
import json
import os

PRESSURE_LEVELS = [50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000]


class ERA5SyntheticGenerator:
    """
    Generate synthetic ERA5-like atmospheric data.

    Creates physically plausible atmospheric fields with:
      - Realistic mean profiles and variability
      - Spatial correlations (large-scale patterns)
      - Temporal evolution (smooth transitions)
      - Proper vertical structure
    """

    # Climatological mean profiles (approximate)
    TEMP_PROFILE = np.array([
        205, 210, 215, 218, 222, 228, 240, 252, 262, 272, 282, 288, 293
    ])  # K, from 50 hPa to 1000 hPa

    WIND_SCALE = np.array([
        15, 20, 25, 30, 35, 30, 20, 15, 12, 10, 8, 7, 6
    ])  # m/s typical magnitudes

    HUMIDITY_PROFILE = np.array([
        0.00001, 0.00005, 0.0001, 0.0003, 0.0005, 0.001,
        0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.014
    ])  # kg/kg

    GEO_PROFILE = np.array([
        20000, 16000, 14000, 12000, 10500, 9000,
        7000, 5500, 4000, 3000, 1500, 800, 100
    ])  # m²/s² (geopotential)

    def __init__(self, n_lat: int = 46, n_lon: int = 90, seed: int = 42):
        self.n_lat = n_lat
        self.n_lon = n_lon
        self.n_nodes = n_lat * n_lon
        self.rng = np.random.RandomState(seed)

        self.lats = np.linspace(-90, 90, n_lat)
        self.lons = np.linspace(0, 360, n_lon, endpoint=False)
        self.lat_grid, self.lon_grid = np.meshgrid(self.lats, self.lons, indexing='ij')
        self.lat_flat = self.lat_grid.flatten()
        self.lon_flat = self.lon_grid.flatten()

    def _add_spatial_pattern(self, shape: Tuple, scale: float = 1.0) -> np.ndarray:
        """Add large-scale spatial patterns using low-frequency harmonics."""
        field = np.zeros(shape)
        n_modes = 5
        for _ in range(n_modes):
            kx = self.rng.randint(1, 6)
            ky = self.rng.randint(1, 6)
            phase = self.rng.uniform(0, 2 * np.pi)
            amp = scale * self.rng.uniform(0.3, 1.0)
            lat_r = np.deg2rad(self.lat_grid)
            lon_r = np.deg2rad(self.lon_grid)
            pattern = amp * np.cos(kx * lat_r + ky * lon_r + phase)
            field += pattern.flatten().reshape(-1, 1) if len(shape) > 1 else pattern.flatten()
        return field

    def generate_state(self, t: float = 0.0) -> Dict[str, np.ndarray]:
        """
        Generate one atmospheric state.

        Args:
            t: time parameter (0-1 range for temporal variation)

        Returns:
            Dict with 'surface', 'pressure_level', 'lat', 'lon' arrays
        """
        n = self.n_nodes
        n_lev = len(PRESSURE_LEVELS)

        # Temperature: climatological profile + spatial + noise
        temp = np.tile(self.TEMP_PROFILE.astype(np.float64), (n, 1))
        temp += self._add_spatial_pattern((n, n_lev), scale=10.0)
        lat_effect = -30 * np.abs(self.lat_flat / 90.0).reshape(-1, 1)
        temp += lat_effect
        temp += self.rng.normal(0, 2, (n, n_lev))
        # Temporal variation
        temp += 3 * np.sin(2 * np.pi * t) * np.cos(np.deg2rad(self.lat_flat)).reshape(-1, 1)

        # Wind (u, v)
        u = np.tile(self.WIND_SCALE.astype(np.float64), (n, 1)) * self._add_spatial_pattern((n, n_lev), scale=0.3)
        jet_lat = 30 + 10 * np.sin(2 * np.pi * t)
        jet_effect = np.exp(-((self.lat_flat - jet_lat) / 15) ** 2).reshape(-1, 1)
        u += 30 * jet_effect * np.array([0.5, 1, 1.5, 2, 2.5, 2, 1.5, 1, 0.5, 0.3, 0.2, 0.1, 0.05])
        u += self.rng.normal(0, 3, (n, n_lev))

        v = np.tile(self.WIND_SCALE.astype(np.float64) * 0.3, (n, 1)) * self._add_spatial_pattern((n, n_lev), scale=0.3)
        v += self.rng.normal(0, 2, (n, n_lev))

        # Specific humidity (non-negative, decrease with height)
        q = np.tile(self.HUMIDITY_PROFILE, (n, 1))
        q *= (1 + 0.3 * self._add_spatial_pattern((n, n_lev), scale=1.0))
        q *= (1 + 0.5 * np.cos(np.deg2rad(self.lat_flat)).reshape(-1, 1))
        q = np.maximum(q + self.rng.normal(0, 0.001, (n, n_lev)), 1e-7)

        # Geopotential
        z = np.tile(self.GEO_PROFILE.astype(np.float64), (n, 1))
        z += self._add_spatial_pattern((n, n_lev), scale=200.0)
        z += self.rng.normal(0, 50, (n, n_lev))

        # Surface variables: mslp, t2m, u10, v10, tp
        mslp = 101325 + self._add_spatial_pattern((n,), scale=2000) + self.rng.normal(0, 500, n)
        t2m = temp[:, -1] + self.rng.normal(0, 2, n)
        u10 = u[:, -1] * 0.7 + self.rng.normal(0, 1, n)
        v10 = v[:, -1] * 0.7 + self.rng.normal(0, 1, n)
        tp = np.maximum(q[:, -1] * 10 * np.abs(self.rng.normal(0, 1, n)), 0)

        surface = np.stack([mslp, t2m, u10, v10, tp], axis=-1)  # (n, 5)

        # Concatenate all pressure-level variables: t, u, v, q, z (each n_lev)
        pressure = np.concatenate([temp, u, v, q, z], axis=-1)  # (n, 5*13=65)

        # Full state vector: surface (5) + pressure (65) = 70
        state = np.concatenate([surface, pressure], axis=-1)

        return {
            'state': state.astype(np.float32),
            'surface': surface.astype(np.float32),
            'temperature': temp.astype(np.float32),
            'u_wind': u.astype(np.float32),
            'v_wind': v.astype(np.float32),
            'humidity': q.astype(np.float32),
            'geopotential': z.astype(np.float32),
            'lat': self.lat_flat.astype(np.float32),
            'lon': self.lon_flat.astype(np.float32),
        }

    def generate_sequence(self, n_steps: int = 20, dt: float = 0.05) -> list:
        """Generate a temporal sequence of atmospheric states."""
        states = []
        for i in range(n_steps):
            t = i * dt
            state = self.generate_state(t)
            state['time_index'] = i
            states.append(state)
        return states

    def generate_training_data(
        self,
        n_samples: int = 50,
        lead_times: list = [1, 4, 20],
    ) -> Dict:
        """
        Generate training pairs (input_state, target_state) for multiple lead times.

        lead_times: in units of 6-hour steps (1=6h, 4=24h, 20=120h)
        """
        max_lt = max(lead_times)
        seq_len = n_samples + max_lt

        sequence = self.generate_sequence(seq_len, dt=0.02)

        data = {lt: {'inputs': [], 'targets': []} for lt in lead_times}

        for i in range(n_samples):
            input_state = sequence[i]['state']
            for lt in lead_times:
                target_state = sequence[i + lt]['state']
                data[lt]['inputs'].append(input_state)
                data[lt]['targets'].append(target_state)

        # Convert to tensors
        for lt in lead_times:
            data[lt]['inputs'] = torch.tensor(np.array(data[lt]['inputs']))
            data[lt]['targets'] = torch.tensor(np.array(data[lt]['targets']))

        data['lat'] = torch.tensor(self.lat_flat)
        data['lon'] = torch.tensor(self.lon_flat)

        return data
