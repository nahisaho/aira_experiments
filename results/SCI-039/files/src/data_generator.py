"""
Synthetic ERA5-like data generator for training and evaluation.
Generates realistic atmospheric fields with proper spatial correlations
and physical relationships between variables.
"""

import numpy as np
import torch


class ERA5SyntheticGenerator:
    """Generate synthetic atmospheric data mimicking ERA5 reanalysis."""

    PRESSURE_LEVELS = [50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000]

    # Climatological mean profiles (approximate)
    T_MEAN_PROFILE = np.array([
        210, 210, 215, 218, 222, 228, 240, 252, 262, 270, 282, 287, 290
    ], dtype=np.float32)

    T_STD_PROFILE = np.array([
        8, 10, 12, 12, 10, 8, 7, 6, 5, 5, 6, 5, 5
    ], dtype=np.float32)

    def __init__(self, resolution=2.5, seed=42):
        self.resolution = resolution
        self.rng = np.random.RandomState(seed)
        self.lats = np.arange(-90, 90 + resolution, resolution)
        self.lons = np.arange(0, 360, resolution)
        self.n_lat = len(self.lats)
        self.n_lon = len(self.lons)
        self.n_nodes = self.n_lat * self.n_lon
        self.n_levels = len(self.PRESSURE_LEVELS)

    def _generate_correlated_field(self, mean, std, correlation_length=10):
        """Generate spatially correlated random field using spectral method."""
        field = self.rng.randn(self.n_lat, self.n_lon)
        # Simple spatial smoothing
        from scipy.ndimage import gaussian_filter
        sigma = correlation_length / self.resolution
        field = gaussian_filter(field, sigma=sigma, mode='wrap')
        # Normalize and scale
        field = (field - field.mean()) / (field.std() + 1e-8)
        field = mean + std * field
        return field

    def generate_temperature_field(self):
        """Generate 3D temperature field (lat, lon, level)."""
        T = np.zeros((self.n_lat, self.n_lon, self.n_levels), dtype=np.float32)
        for k in range(self.n_levels):
            # Add latitude dependence
            lat_factor = np.cos(np.radians(self.lats))[:, None] * np.ones((1, self.n_lon))
            base = self.T_MEAN_PROFILE[k] + 20 * lat_factor * (1 if k > 6 else 0.3)
            anomaly = self._generate_correlated_field(0, self.T_STD_PROFILE[k])
            T[:, :, k] = base + anomaly
        return T.reshape(self.n_nodes, self.n_levels)

    def generate_wind_field(self):
        """Generate u and v wind components with jet stream structure."""
        u = np.zeros((self.n_lat, self.n_lon, self.n_levels), dtype=np.float32)
        v = np.zeros((self.n_lat, self.n_lon, self.n_levels), dtype=np.float32)

        for k in range(self.n_levels):
            p = self.PRESSURE_LEVELS[k]
            # Jet stream around 200-300 hPa
            jet_strength = 30 * np.exp(-((p - 250) / 100) ** 2)
            lat_jet = np.exp(-((self.lats - 30) / 15) ** 2) + \
                      0.5 * np.exp(-((self.lats + 30) / 15) ** 2)

            u_base = jet_strength * lat_jet[:, None] * np.ones((1, self.n_lon))
            u[:, :, k] = u_base + self._generate_correlated_field(0, 5)
            v[:, :, k] = self._generate_correlated_field(0, 3)

        return (u.reshape(self.n_nodes, self.n_levels),
                v.reshape(self.n_nodes, self.n_levels))

    def generate_humidity_field(self, T):
        """Generate specific humidity field consistent with temperature."""
        T_2d = T.reshape(self.n_lat, self.n_lon, self.n_levels)
        q = np.zeros_like(T_2d)

        for k in range(self.n_levels):
            p = self.PRESSURE_LEVELS[k]
            # Clausius-Clapeyron approximation
            es = 611.2 * np.exp(17.67 * (T_2d[:, :, k] - 273.15) / (T_2d[:, :, k] - 29.65))
            qs = 0.622 * es / (p * 100)
            # Relative humidity varies
            rh = np.clip(self._generate_correlated_field(0.5, 0.2), 0.05, 0.95)
            q[:, :, k] = rh * qs

        return q.reshape(self.n_nodes, self.n_levels)

    def generate_surface_vars(self, T, u, v):
        """Generate surface variables: 2m T, 10m u, 10m v, MSLP."""
        T_2d = T.reshape(self.n_lat, self.n_lon, self.n_levels)
        u_2d = u.reshape(self.n_lat, self.n_lon, self.n_levels)
        v_2d = v.reshape(self.n_lat, self.n_lon, self.n_levels)

        t2m = T_2d[:, :, -1] + self._generate_correlated_field(0, 2)
        u10 = u_2d[:, :, -1] * 0.7 + self._generate_correlated_field(0, 1)
        v10 = v_2d[:, :, -1] * 0.7 + self._generate_correlated_field(0, 1)
        mslp = self._generate_correlated_field(101325, 1500)

        surface = np.stack([t2m, u10, v10, mslp], axis=-1)
        return surface.reshape(self.n_nodes, 4).astype(np.float32)

    def generate_sample(self):
        """Generate one complete atmospheric state."""
        T = self.generate_temperature_field()
        u, v = self.generate_wind_field()
        q = self.generate_humidity_field(T)
        surface = self.generate_surface_vars(T, u, v)

        pressure_vars = np.stack([T, u, v, q], axis=1)  # (N, 4, 13)
        return pressure_vars, surface

    def generate_forecast_pair(self, lead_time_hours=6):
        """Generate input/target pair with realistic temporal evolution."""
        # Input state
        p_in, s_in = self.generate_sample()

        # Evolve state (simplified advection + noise)
        dt_factor = lead_time_hours / 6.0
        noise_scale = 0.02 * dt_factor

        p_out = p_in.copy()
        # Temperature tendency
        p_out[:, 0, :] += self.rng.randn(*p_in[:, 0, :].shape) * self.T_STD_PROFILE * noise_scale
        # Wind tendency
        p_out[:, 1, :] += self.rng.randn(*p_in[:, 1, :].shape) * 2.0 * noise_scale
        p_out[:, 2, :] += self.rng.randn(*p_in[:, 2, :].shape) * 1.5 * noise_scale
        # Humidity tendency
        p_out[:, 3, :] = np.clip(
            p_out[:, 3, :] + self.rng.randn(*p_in[:, 3, :].shape) * 0.001 * noise_scale,
            0, None
        )

        s_out = s_in.copy()
        s_out[:, 0] += self.rng.randn(self.n_nodes) * 1.0 * noise_scale
        s_out[:, 1] += self.rng.randn(self.n_nodes) * 0.5 * noise_scale
        s_out[:, 2] += self.rng.randn(self.n_nodes) * 0.5 * noise_scale
        s_out[:, 3] += self.rng.randn(self.n_nodes) * 100 * noise_scale

        return (p_in.astype(np.float32), s_in.astype(np.float32),
                p_out.astype(np.float32), s_out.astype(np.float32))

    def generate_dataset(self, n_samples=50, lead_times=[6, 24, 120]):
        """Generate full dataset for multiple lead times."""
        datasets = {}
        for lt in lead_times:
            data = []
            for i in range(n_samples):
                self.rng = np.random.RandomState(42 + i * 100 + lt)
                p_in, s_in, p_out, s_out = self.generate_forecast_pair(lt)
                data.append({
                    'input_pressure': torch.from_numpy(p_in),
                    'input_surface': torch.from_numpy(s_in),
                    'target_pressure': torch.from_numpy(p_out),
                    'target_surface': torch.from_numpy(s_out),
                })
            datasets[lt] = data
        return datasets
