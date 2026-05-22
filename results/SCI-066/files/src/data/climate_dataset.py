"""
Data pipeline for CMIP6/ClimateBench climate data.

Handles loading, preprocessing, and batching of climate model outputs
for training the ESM emulator. Supports xarray-based data access.
"""

import numpy as np
import torch
from torch.utils.data import Dataset
from typing import Dict, Tuple, Optional


# SSP scenario mapping
SSP_SCENARIOS = {
    "ssp126": 0,
    "ssp245": 1,
    "ssp370": 2,
    "ssp585": 3,
}

# Climate variable names in CMIP6
CLIMATE_VARS = {
    "tas": "Near-Surface Air Temperature",
    "pr": "Precipitation",
    "zos": "Sea Surface Height",
}


def generate_synthetic_cmip6_data(
    n_years: int = 100,
    spatial_size: Tuple[int, int] = (64, 128),
    n_scenarios: int = 4,
    seed: int = 42,
) -> Dict:
    """
    Generate synthetic CMIP6-like data for development and testing.

    Produces physically plausible climate fields with:
    - Latitudinal temperature gradient
    - ENSO-like variability patterns
    - Scenario-dependent warming trends
    - Realistic spatial correlation structure
    """
    rng = np.random.RandomState(seed)
    H, W = spatial_size

    lat = np.linspace(-90, 90, H)
    lon = np.linspace(0, 360, W)
    lat_grid, lon_grid = np.meshgrid(lat, lon, indexing="ij")

    # Base temperature: latitudinal gradient + seasonal pattern
    base_temp = 288.0 - 40.0 * np.abs(np.sin(np.radians(lat_grid)))

    # Warming trends per scenario (K/century)
    warming_rates = {0: 1.5, 1: 2.7, 2: 3.6, 3: 4.8}

    data = {}
    for scenario_id in range(n_scenarios):
        years = np.arange(n_years)
        fields = np.zeros((n_years, 3, H, W), dtype=np.float32)

        for t in range(n_years):
            trend = warming_rates[scenario_id] * t / 100.0

            # Temperature: base + trend + noise
            temp_noise = rng.randn(H, W).astype(np.float32) * 2.0
            # Spatial smoothing
            for _ in range(5):
                temp_noise[1:-1, 1:-1] = 0.25 * (
                    temp_noise[:-2, 1:-1] + temp_noise[2:, 1:-1]
                    + temp_noise[1:-1, :-2] + temp_noise[1:-1, 2:]
                )
            fields[t, 0] = base_temp + trend + temp_noise

            # Precipitation: exponential relationship with temperature
            precip_base = 3.0 * np.exp(-0.5 * (lat_grid / 15) ** 2)  # ITCZ
            precip_noise = np.abs(rng.randn(H, W).astype(np.float32)) * 0.5
            fields[t, 1] = precip_base * (1 + 0.07 * trend) + precip_noise

            # Sea level: global mean rise + regional patterns
            sl_trend = 0.003 * warming_rates[scenario_id] * t  # m/decade
            sl_pattern = 0.02 * np.sin(np.radians(lon_grid * 2))
            sl_noise = rng.randn(H, W).astype(np.float32) * 0.01
            fields[t, 2] = sl_trend + sl_pattern + sl_noise

        # Forcing vector: CO2, CH4, N2O, aerosol, solar, etc.
        forcing = np.zeros((n_years, 8), dtype=np.float32)
        forcing[:, 0] = 5.35 * np.log(
            (400 + warming_rates[scenario_id] * 10 * years / 100) / 280
        )  # CO2 forcing
        forcing[:, 1] = 0.036 * (1800 + 5 * years) ** 0.5  # CH4
        forcing[:, 2] = 0.12 * (320 + 0.5 * years) ** 0.5  # N2O
        forcing[:, 3] = -0.5 + 0.002 * years  # Aerosol
        forcing[:, 4] = 0.1 * np.sin(2 * np.pi * years / 11)  # Solar cycle
        forcing[:, 5:] = rng.randn(n_years, 3).astype(np.float32) * 0.05

        data[scenario_id] = {
            "fields": fields,
            "forcing": forcing,
            "years": years,
        }

    return data


class ClimateDataset(Dataset):
    """
    PyTorch Dataset for climate field sequences.

    Produces (sequence, target, scenario_id, forcing) tuples
    for training the ESM emulator.
    """

    def __init__(self, data: Dict, seq_length: int = 10,
                 normalize: bool = True):
        self.seq_length = seq_length
        self.samples = []
        self.stats = None

        all_fields = []
        for scenario_id, scenario_data in data.items():
            fields = scenario_data["fields"]
            forcing = scenario_data["forcing"]
            n_years = len(fields)

            for t in range(n_years - seq_length):
                self.samples.append({
                    "sequence": fields[t:t + seq_length],
                    "target": fields[t + seq_length],
                    "scenario_id": scenario_id,
                    "forcing": forcing[t + seq_length],
                })
                all_fields.append(fields[t:t + seq_length + 1])

        if normalize:
            all_fields = np.concatenate(all_fields, axis=0)
            self.stats = {
                "mean": all_fields.mean(axis=(0, 2, 3), keepdims=True),
                "std": all_fields.std(axis=(0, 2, 3), keepdims=True) + 1e-8,
            }

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, idx: int) -> Tuple:
        sample = self.samples[idx]
        seq = torch.from_numpy(sample["sequence"].copy())
        target = torch.from_numpy(sample["target"].copy())

        if self.stats is not None:
            mean = torch.from_numpy(self.stats["mean"])
            std = torch.from_numpy(self.stats["std"])
            seq = (seq - mean) / std
            target = (target - mean.squeeze(0)) / std.squeeze(0)

        return (
            seq,
            target,
            torch.tensor(sample["scenario_id"], dtype=torch.long),
            torch.from_numpy(sample["forcing"].copy()),
        )


class ClimateBenchDataLoader:
    """
    Interface for loading ClimateBench benchmark data.

    ClimateBench provides standardized CMIP6 outputs for emulator
    benchmarking across NorESM2, CanESM5, MIROC6, etc.
    """

    SUPPORTED_MODELS = [
        "NorESM2-LM", "CanESM5", "MIROC6", "MPI-ESM1-2-LR",
        "UKESM1-0-LL", "IPSL-CM6A-LR", "ACCESS-ESM1-5",
    ]

    SUPPORTED_VARIABLES = ["tas", "pr", "zos", "psl"]

    def __init__(self, data_dir: str = "data/climatebench"):
        self.data_dir = data_dir

    def load_xarray(self, model: str, variable: str,
                    scenario: str) -> Optional[object]:
        """
        Load ClimateBench data as xarray Dataset.

        In production: loads from NetCDF files via xarray.
        Returns None if xarray/data not available.
        """
        try:
            import xarray as xr
            path = f"{self.data_dir}/{model}/{variable}_{scenario}.nc"
            return xr.open_dataset(path)
        except (ImportError, FileNotFoundError):
            return None

    def get_synthetic_benchmark(self) -> Dict:
        """Fallback to synthetic data for testing."""
        return generate_synthetic_cmip6_data()
