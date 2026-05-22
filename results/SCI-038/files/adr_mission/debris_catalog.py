from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import DATA_DIR, GM_EARTH, R_EARTH, SEED

CATALOG_PATH = DATA_DIR / "debris_catalog.csv"


def _atmospheric_density(altitude_m: np.ndarray) -> np.ndarray:
    rho0 = 3.5e-12
    scale_height = 85e3
    return rho0 * np.exp(-(altitude_m - 400e3) / scale_height)


def generate_debris_catalog(num_objects: int = 50, seed: int = SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)

    semi_major_axis_km = rng.uniform(6571.0, 7371.0, num_objects)
    eccentricity = np.clip(rng.beta(1.5, 18.0, num_objects), 0.0, 0.08)
    inclination_deg = np.concatenate(
        [
            rng.normal(53.0, 4.0, num_objects // 3),
            rng.normal(74.0, 6.0, num_objects // 3),
            rng.normal(98.0, 3.0, num_objects - 2 * (num_objects // 3)),
        ]
    )
    rng.shuffle(inclination_deg)
    inclination_deg = np.clip(inclination_deg, 0.0, 120.0)
    raan_deg = rng.uniform(0.0, 360.0, num_objects)
    aop_deg = rng.uniform(0.0, 360.0, num_objects)
    true_anomaly_deg = rng.uniform(0.0, 360.0, num_objects)

    mass_kg = np.exp(rng.uniform(np.log(10.0), np.log(2000.0), num_objects))
    size_m = np.exp(rng.uniform(np.log(0.1), np.log(5.0), num_objects))
    area_m2 = np.pi * (size_m / 2.0) ** 2
    radar_cross_section_m2 = area_m2 * rng.uniform(0.7, 1.8, num_objects)

    altitude_m = semi_major_axis_km * 1e3 - R_EARTH
    orbital_velocity = np.sqrt(GM_EARTH / (semi_major_axis_km * 1e3))
    cd = 2.2
    ballistic_coefficient = mass_kg / (cd * area_m2)
    rho = _atmospheric_density(altitude_m)

    decay_rate_m_per_day = 0.5 * cd * area_m2 / mass_kg * rho * orbital_velocity**2 * 86400.0
    decay_rate_km_per_day = np.clip(decay_rate_m_per_day / 1000.0, 1e-6, None)
    decay_lifetime_days = np.clip(altitude_m / np.clip(decay_rate_m_per_day, 1e-9, None) / 86400.0, 30.0, 25000.0)

    catalog = pd.DataFrame(
        {
            "debris_id": [f"DEBRIS-{i:03d}" for i in range(1, num_objects + 1)],
            "semi_major_axis_km": semi_major_axis_km,
            "eccentricity": eccentricity,
            "inclination_deg": inclination_deg,
            "raan_deg": raan_deg,
            "aop_deg": aop_deg,
            "true_anomaly_deg": true_anomaly_deg,
            "mass_kg": mass_kg,
            "size_m": size_m,
            "area_m2": area_m2,
            "radar_cross_section_m2": radar_cross_section_m2,
            "ballistic_coefficient": ballistic_coefficient,
            "orbital_velocity_m_s": orbital_velocity,
            "decay_rate_km_day": decay_rate_km_per_day,
            "decay_lifetime_days": decay_lifetime_days,
            "altitude_km": altitude_m / 1000.0,
        }
    )

    catalog.to_csv(CATALOG_PATH, index=False)
    return catalog


if __name__ == "__main__":
    df = generate_debris_catalog()
    print(json.dumps({"catalog_path": str(CATALOG_PATH), "objects": int(len(df))}, indent=2))
