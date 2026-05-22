"""
Case study data and validation for Sakurajima and Aso volcanoes.

Generates synthetic but realistic test datasets based on published
geodetic observations from:
  - Sakurajima: Iguchi et al. (2013), Hotta et al. (2016)
  - Aso: Ohkura et al. (2009), Sudo et al. (2006)

All coordinates are in local ENU frame centered on the volcano summit.
"""

import numpy as np
from typing import Dict, Tuple
from dataclasses import dataclass

from .source_models import (
    MogiSource, SpheroidSource,
    mogi_displacement, mogi_gravity,
    spheroid_displacement
)
from .bayesian_inversion import InversionData


# ==============================================================================
# Sakurajima Volcano
# ==============================================================================

def sakurajima_source_params() -> Dict:
    """
    Published source parameters for Sakurajima volcano.

    Based on Hotta et al. (2016) and Iguchi et al. (2013):
    - Shallow source (Showa crater): ~1.5 km depth, beneath summit
    - Deep source (Aira caldera): ~10 km depth, NE of volcano

    Returns dict with 'shallow' and 'deep' Mogi sources.
    """
    shallow = MogiSource(
        x=500.0,       # 500 m east of summit
        y=-200.0,      # 200 m south
        d=1500.0,      # 1.5 km depth
        dV=0.5e6,      # 0.5 × 10^6 m^3
        nu=0.25
    )

    deep = MogiSource(
        x=3000.0,      # 3 km NE (Aira caldera center)
        y=5000.0,      # 5 km north
        d=10000.0,     # 10 km depth
        dV=8.0e6,      # 8 × 10^6 m^3
        nu=0.25
    )

    spheroid_deep = SpheroidSource(
        x=3000.0, y=5000.0, d=10000.0,
        a=3000.0,      # 3 km semi-major
        b=1500.0,      # 1.5 km semi-minor (prolate)
        dP=15e6,       # 15 MPa
        strike=30.0,   # NE-SW
        dip=80.0,      # near-vertical
        nu=0.25,
        mu=3.0e10
    )

    return {
        'shallow': shallow,
        'deep': deep,
        'spheroid_deep': spheroid_deep,
        'volcano_name': 'Sakurajima',
        'summit_elevation': 1117,  # m
        'latitude': 31.585,
        'longitude': 130.657
    }


def sakurajima_gnss_network() -> Tuple[np.ndarray, np.ndarray, list]:
    """
    Synthetic GNSS station coordinates mimicking the actual network.

    Returns (obs_x, obs_y, station_names) in local ENU [m].
    """
    stations = {
        'SAKR': (0, 0),
        'HARU': (-3000, 2000),
        'KITA': (1000, 4000),
        'FURU': (-2500, -1500),
        'ARIM': (5000, 3000),
        'KAGO': (-8000, 8000),
        'TARU': (2000, -3000),
        'YOSH': (-5000, -4000),
        'KOKU': (4000, -2000),
        'SBTK': (3000, 7000),
        'AIRA': (4000, 6000),
        'KURA': (-1500, 5000),
    }
    names = list(stations.keys())
    coords = np.array(list(stations.values()))
    return coords[:, 0], coords[:, 1], names


def generate_sakurajima_data(
    noise_level: float = 0.003,
    include_insar: bool = True,
    include_gravity: bool = True,
    seed: int = 42
) -> Tuple[InversionData, Dict]:
    """
    Generate synthetic observation data for Sakurajima.

    Parameters
    ----------
    noise_level : GNSS noise std [m] (default 3 mm)
    include_insar : generate InSAR data
    include_gravity : generate gravity data
    seed : random seed

    Returns
    -------
    data : InversionData
    truth : dict with true source parameters
    """
    rng = np.random.default_rng(seed)
    params = sakurajima_source_params()

    # GNSS network
    gnss_x, gnss_y, stn_names = sakurajima_gnss_network()
    n_gnss = len(gnss_x)

    # True displacement (two-source model)
    disp_shallow = mogi_displacement(gnss_x, gnss_y, params['shallow'])
    disp_deep = mogi_displacement(gnss_x, gnss_y, params['deep'])
    true_disp = disp_shallow + disp_deep

    gnss_disp = true_disp + rng.normal(0, noise_level, true_disp.shape)
    gnss_sigma = np.full_like(gnss_disp, noise_level)
    # Vertical is noisier
    gnss_sigma[:, 2] *= 3.0

    # All observation coordinates
    all_x = gnss_x.copy()
    all_y = gnss_y.copy()
    gnss_idx = np.arange(n_gnss)

    # InSAR data
    insar_los = None
    insar_sigma = None
    insar_idx = None
    insar_look = None

    if include_insar:
        # Regular grid
        ix = np.linspace(-10000, 10000, 40)
        iy = np.linspace(-10000, 10000, 40)
        ixx, iyy = np.meshgrid(ix, iy)
        insar_x = ixx.flatten()
        insar_y = iyy.flatten()
        n_insar = len(insar_x)

        # Look vector (ascending, ~34° incidence, heading=-12°)
        inc = np.radians(34)
        head = np.radians(-12)
        look = np.array([
            -np.sin(inc) * np.sin(head),
            np.sin(inc) * np.cos(head),
            np.cos(inc)
        ])
        insar_look_arr = np.tile(look, (n_insar, 1))

        disp_s = mogi_displacement(insar_x, insar_y, params['shallow'])
        disp_d = mogi_displacement(insar_x, insar_y, params['deep'])
        true_insar_3d = disp_s + disp_d
        true_los = np.sum(true_insar_3d * insar_look_arr, axis=1)

        insar_noise = 0.005  # 5 mm
        insar_los = true_los + rng.normal(0, insar_noise, n_insar)
        insar_sigma = np.full(n_insar, insar_noise)
        insar_look = insar_look_arr

        # Add to coordinate arrays
        offset = len(all_x)
        all_x = np.concatenate([all_x, insar_x])
        all_y = np.concatenate([all_y, insar_y])
        insar_idx = np.arange(offset, offset + n_insar)

    # Gravity data
    gravity = None
    gravity_sigma = None
    gravity_idx = None

    if include_gravity:
        grav_x = gnss_x[:6]
        grav_y = gnss_y[:6]
        n_grav = len(grav_x)

        grav_s = mogi_gravity(grav_x, grav_y, params['shallow'])
        grav_d = mogi_gravity(grav_x, grav_y, params['deep'])
        true_grav = grav_s + grav_d

        grav_noise = 5.0  # 5 µGal
        gravity = true_grav + rng.normal(0, grav_noise, n_grav)
        gravity_sigma = np.full(n_grav, grav_noise)
        gravity_idx = np.arange(n_grav)  # same coords as first 6 GNSS

    data = InversionData(
        obs_x=all_x,
        obs_y=all_y,
        gnss_disp=gnss_disp,
        gnss_sigma=gnss_sigma,
        gnss_idx=gnss_idx,
        insar_los=insar_los,
        insar_sigma=insar_sigma,
        insar_idx=insar_idx,
        insar_look=insar_look,
        gravity=gravity,
        gravity_sigma=gravity_sigma,
        gravity_idx=gravity_idx,
    )

    truth = {
        'shallow': params['shallow'],
        'deep': params['deep'],
        'station_names': stn_names,
    }

    return data, truth


# ==============================================================================
# Aso Volcano
# ==============================================================================

def aso_source_params() -> Dict:
    """
    Published source parameters for Aso volcano.

    Based on Ohkura et al. (2009) and Sudo et al. (2006):
    - Shallow source: ~1-2 km depth beneath Nakadake crater
    - Intermediate source: ~4-6 km depth (Kusasenri area)

    Returns dict with source models.
    """
    shallow = MogiSource(
        x=0.0,
        y=0.0,
        d=1500.0,       # 1.5 km
        dV=0.3e6,       # 0.3 × 10^6 m^3
        nu=0.25
    )

    intermediate = MogiSource(
        x=-2000.0,      # 2 km west (Kusasenri)
        y=1000.0,       # 1 km north
        d=5000.0,       # 5 km depth
        dV=3.0e6,       # 3 × 10^6 m^3
        nu=0.25
    )

    spheroid = SpheroidSource(
        x=-2000.0, y=1000.0, d=5000.0,
        a=2000.0,       # 2 km semi-major
        b=1000.0,       # 1 km semi-minor
        dP=10e6,        # 10 MPa
        strike=45.0,    # NE-SW
        dip=70.0,
        nu=0.25,
        mu=2.5e10
    )

    return {
        'shallow': shallow,
        'intermediate': intermediate,
        'spheroid': spheroid,
        'volcano_name': 'Aso',
        'summit_elevation': 1592,
        'latitude': 32.884,
        'longitude': 131.104
    }


def aso_gnss_network() -> Tuple[np.ndarray, np.ndarray, list]:
    """GNSS station coordinates for Aso volcano."""
    stations = {
        'ASON': (0, 3000),
        'ASOE': (3000, 0),
        'ASOS': (0, -3000),
        'ASOW': (-3000, 0),
        'NAKA': (500, 200),
        'KUSA': (-2000, 1500),
        'TAKA': (5000, 5000),
        'SUGO': (-5000, -3000),
        'CHOY': (2000, -5000),
        'OTTA': (-6000, 4000),
    }
    names = list(stations.keys())
    coords = np.array(list(stations.values()))
    return coords[:, 0], coords[:, 1], names


def generate_aso_data(
    noise_level: float = 0.003,
    include_insar: bool = True,
    include_gravity: bool = True,
    seed: int = 123
) -> Tuple[InversionData, Dict]:
    """Generate synthetic observation data for Aso volcano."""
    rng = np.random.default_rng(seed)
    params = aso_source_params()

    gnss_x, gnss_y, stn_names = aso_gnss_network()
    n_gnss = len(gnss_x)

    disp_s = mogi_displacement(gnss_x, gnss_y, params['shallow'])
    disp_i = mogi_displacement(gnss_x, gnss_y, params['intermediate'])
    true_disp = disp_s + disp_i

    gnss_disp = true_disp + rng.normal(0, noise_level, true_disp.shape)
    gnss_sigma = np.full_like(gnss_disp, noise_level)
    gnss_sigma[:, 2] *= 3.0

    all_x = gnss_x.copy()
    all_y = gnss_y.copy()
    gnss_idx = np.arange(n_gnss)

    insar_los = insar_sigma = insar_idx = insar_look = None
    if include_insar:
        ix = np.linspace(-8000, 8000, 32)
        iy = np.linspace(-8000, 8000, 32)
        ixx, iyy = np.meshgrid(ix, iy)
        insar_x = ixx.flatten()
        insar_y = iyy.flatten()
        n_insar = len(insar_x)

        inc = np.radians(38)
        head = np.radians(-10)
        look = np.array([
            -np.sin(inc) * np.sin(head),
            np.sin(inc) * np.cos(head),
            np.cos(inc)
        ])
        insar_look_arr = np.tile(look, (n_insar, 1))

        disp_s2 = mogi_displacement(insar_x, insar_y, params['shallow'])
        disp_i2 = mogi_displacement(insar_x, insar_y, params['intermediate'])
        true_los = np.sum((disp_s2 + disp_i2) * insar_look_arr, axis=1)

        insar_noise = 0.004
        insar_los = true_los + rng.normal(0, insar_noise, n_insar)
        insar_sigma = np.full(n_insar, insar_noise)
        insar_look = insar_look_arr

        offset = len(all_x)
        all_x = np.concatenate([all_x, insar_x])
        all_y = np.concatenate([all_y, insar_y])
        insar_idx = np.arange(offset, offset + n_insar)

    gravity = gravity_sigma = gravity_idx = None
    if include_gravity:
        grav_x = gnss_x[:5]
        grav_y = gnss_y[:5]
        n_grav = len(grav_x)

        grav_s = mogi_gravity(grav_x, grav_y, params['shallow'])
        grav_i = mogi_gravity(grav_x, grav_y, params['intermediate'])
        true_grav = grav_s + grav_i

        grav_noise = 5.0
        gravity = true_grav + rng.normal(0, grav_noise, n_grav)
        gravity_sigma = np.full(n_grav, grav_noise)
        gravity_idx = np.arange(n_grav)

    data = InversionData(
        obs_x=all_x, obs_y=all_y,
        gnss_disp=gnss_disp, gnss_sigma=gnss_sigma, gnss_idx=gnss_idx,
        insar_los=insar_los, insar_sigma=insar_sigma,
        insar_idx=insar_idx, insar_look=insar_look,
        gravity=gravity, gravity_sigma=gravity_sigma, gravity_idx=gravity_idx,
    )

    truth = {
        'shallow': params['shallow'],
        'intermediate': params['intermediate'],
        'station_names': stn_names,
    }

    return data, truth


# ==============================================================================
# Time series generation for Kalman filter testing
# ==============================================================================

def generate_timeseries_data(
    volcano: str = "sakurajima",
    n_epochs: int = 365,
    dt_days: float = 1.0,
    inflation_rate: float = 1e4,   # m^3/day
    noise_level: float = 0.002,
    eruption_day: int = 200,
    eruption_volume: float = -2e6,
    seed: int = 42
) -> Dict:
    """
    Generate synthetic time series data with inflation/deflation events.

    Parameters
    ----------
    volcano : "sakurajima" or "aso"
    n_epochs : number of time steps
    dt_days : time step in days
    inflation_rate : steady inflation rate [m^3/day]
    noise_level : displacement noise [m]
    eruption_day : day of eruption (deflation event)
    eruption_volume : volume loss during eruption [m^3]

    Returns
    -------
    dict with 'times', 'observations', 'true_states', etc.
    """
    rng = np.random.default_rng(seed)

    if volcano == "sakurajima":
        params = sakurajima_source_params()
        source = params['deep']
        obs_x, obs_y, names = sakurajima_gnss_network()
    else:
        params = aso_source_params()
        source = params['intermediate']
        obs_x, obs_y, names = aso_gnss_network()

    times = np.arange(n_epochs) * dt_days
    n_obs = len(obs_x)

    true_dV = np.zeros(n_epochs)
    true_rate = np.zeros(n_epochs)

    # Build volume change history
    for i in range(n_epochs):
        if i < eruption_day:
            true_rate[i] = inflation_rate
            true_dV[i] = inflation_rate * times[i]
        elif i == eruption_day:
            true_rate[i] = eruption_volume / dt_days
            true_dV[i] = true_dV[i-1] + eruption_volume
        else:
            # Post-eruption re-inflation (slower)
            true_rate[i] = inflation_rate * 0.5
            true_dV[i] = true_dV[i-1] + true_rate[i] * dt_days

    # Generate observations
    observations = []
    true_states = []

    for i in range(n_epochs):
        src = MogiSource(
            x=source.x, y=source.y, d=source.d,
            dV=true_dV[i], nu=source.nu
        )
        true_disp = mogi_displacement(obs_x, obs_y, src)
        obs_disp = true_disp + rng.normal(0, noise_level, true_disp.shape)

        observations.append(obs_disp.flatten())
        true_states.append([source.x, source.y, source.d,
                           true_dV[i], true_rate[i]])

    return {
        'times': times,
        'observations': observations,
        'true_states': np.array(true_states),
        'true_dV': true_dV,
        'true_rate': true_rate,
        'obs_x': obs_x,
        'obs_y': obs_y,
        'station_names': names,
        'noise_level': noise_level,
        'n_obs_per_epoch': n_obs * 3,
    }
