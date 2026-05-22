"""
Module 1: Resin Flow Simulation
- Hele-Shaw approximation for thin-wall injection molding
- 3D Navier-Stokes based flow analysis (simplified FDM)
- Cross-WLF viscosity model for polymer melts
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict


@dataclass
class ResinProperties:
    """Polymer melt properties for flow simulation."""
    name: str = "PA66-GF30"
    density: float = 1350.0
    specific_heat: float = 1700.0
    thermal_conductivity: float = 0.30
    n_power: float = 0.3
    tau_star: float = 1.5e5
    D1: float = 1.5e12
    D2: float = 263.15
    D3: float = 0.0
    A1: float = 25.0
    A2: float = 51.6
    melt_temperature: float = 553.15
    no_flow_temperature: float = 493.15


@dataclass
class MoldGeometry:
    """Mold cavity geometry (simplified rectangular plate)."""
    length: float = 0.200
    width: float = 0.100
    thickness: float = 0.003
    gate_width: float = 0.005
    gate_thickness: float = 0.002
    nx: int = 50
    ny: int = 25
    nz: int = 10


@dataclass
class ProcessConditions:
    """Injection molding process parameters."""
    injection_pressure: float = 80e6
    packing_pressure: float = 50e6
    injection_speed: float = 0.05
    melt_temperature: float = 553.15
    mold_temperature: float = 353.15
    injection_time: float = 1.5
    packing_time: float = 8.0
    cooling_time: float = 20.0


def cross_wlf_viscosity(shear_rate: np.ndarray, temperature: float,
                         resin: ResinProperties) -> np.ndarray:
    """Cross-WLF viscosity model."""
    T_ref = resin.D2
    if temperature <= T_ref:
        return np.full_like(shear_rate, 1e6)
    eta_0 = resin.D1 * np.exp(-resin.A1 * (temperature - T_ref) /
                               (resin.A2 + (temperature - T_ref)))
    eta = eta_0 / (1.0 + (eta_0 * shear_rate / resin.tau_star) ** (1.0 - resin.n_power))
    return eta


class HeleShawSolver:
    """Hele-Shaw approximation solver for thin-wall cavity filling."""

    def __init__(self, geometry: MoldGeometry, resin: ResinProperties,
                 conditions: ProcessConditions):
        self.geom = geometry
        self.resin = resin
        self.cond = conditions
        self.dx = geometry.length / geometry.nx
        self.dy = geometry.width / geometry.ny
        self.dz = geometry.thickness / geometry.nz
        self.pressure = np.zeros((geometry.nx, geometry.ny))
        self.temperature = np.full((geometry.nx, geometry.ny, geometry.nz),
                                    conditions.melt_temperature)
        self.fill_fraction = np.zeros((geometry.nx, geometry.ny))
        self.viscosity = np.zeros((geometry.nx, geometry.ny))
        self.shear_rate = np.zeros((geometry.nx, geometry.ny))

    def compute_fluidity(self, i: int, j: int) -> float:
        h = self.geom.thickness / 2.0
        z_points = np.linspace(-h, h, self.geom.nz)
        dz = z_points[1] - z_points[0]
        T_avg = np.mean(self.temperature[i, j, :])
        gamma_dot = self.shear_rate[i, j] if self.shear_rate[i, j] > 0 else 1.0
        eta = cross_wlf_viscosity(np.array([gamma_dot]), T_avg, self.resin)[0]
        self.viscosity[i, j] = eta
        S = np.sum(z_points ** 2 / eta) * dz
        return S

    def solve_pressure_field(self) -> np.ndarray:
        P = self.pressure.copy()
        S = np.zeros((self.geom.nx, self.geom.ny))
        for i in range(self.geom.nx):
            for j in range(self.geom.ny):
                if self.fill_fraction[i, j] > 0.01:
                    S[i, j] = self.compute_fluidity(i, j)
        P[0, self.geom.ny // 2 - 1:self.geom.ny // 2 + 1] = self.cond.injection_pressure
        for iteration in range(200):
            P_old = P.copy()
            for i in range(1, self.geom.nx - 1):
                for j in range(1, self.geom.ny - 1):
                    if self.fill_fraction[i, j] < 0.01:
                        continue
                    sx = S[i, j] / self.dx ** 2
                    sy = S[i, j] / self.dy ** 2
                    denom = 2 * (sx + sy)
                    if denom < 1e-20:
                        continue
                    P[i, j] = (sx * (P[i + 1, j] + P[i - 1, j]) +
                               sy * (P[i, j + 1] + P[i, j - 1])) / denom
            residual = np.max(np.abs(P - P_old)) / (np.max(np.abs(P)) + 1e-10)
            if residual < 1e-4:
                break
        self.pressure = P
        return P

    def simulate_filling(self, n_steps: int = 30) -> Dict:
        dt = self.cond.injection_time / n_steps
        flow_front_x = 0.0
        v_inj = self.cond.injection_speed
        history = {
            'time': [], 'fill_percent': [], 'max_pressure': [],
            'avg_temperature': [], 'flow_front': [], 'max_shear_rate': []
        }
        for step in range(n_steps):
            t = step * dt
            flow_front_x = min(flow_front_x + v_inj * dt, self.geom.length)
            front_idx = int(flow_front_x / self.dx)
            for i in range(min(front_idx + 1, self.geom.nx)):
                for j in range(self.geom.ny):
                    self.fill_fraction[i, j] = 1.0
            Q = v_inj * self.geom.width * self.geom.thickness
            for i in range(min(front_idx + 1, self.geom.nx)):
                for j in range(self.geom.ny):
                    self.shear_rate[i, j] = 6 * Q / (self.geom.width * self.geom.thickness ** 2)
            alpha = self.resin.thermal_conductivity / (self.resin.density * self.resin.specific_heat)
            for i in range(min(front_idx + 1, self.geom.nx)):
                for j in range(self.geom.ny):
                    for k in range(self.geom.nz):
                        z_rel = k / (self.geom.nz - 1)
                        wall_factor = 4 * z_rel * (1 - z_rel)
                        cooling = alpha * dt / (self.geom.thickness / 2) ** 2
                        self.temperature[i, j, k] -= cooling * (
                            self.temperature[i, j, k] - self.cond.mold_temperature) * (1 - wall_factor * 0.3)
                        gamma = self.shear_rate[i, j]
                        eta = cross_wlf_viscosity(np.array([gamma]),
                                                   self.temperature[i, j, k], self.resin)[0]
                        q_visc = eta * gamma ** 2
                        dT_visc = q_visc * dt / (self.resin.density * self.resin.specific_heat)
                        self.temperature[i, j, k] += min(dT_visc, 2.0)
            self.solve_pressure_field()
            fill_pct = np.sum(self.fill_fraction) / (self.geom.nx * self.geom.ny) * 100
            history['time'].append(t)
            history['fill_percent'].append(fill_pct)
            history['max_pressure'].append(np.max(self.pressure) / 1e6)
            history['avg_temperature'].append(np.mean(self.temperature[
                self.fill_fraction > 0.5]) - 273.15 if np.any(self.fill_fraction > 0.5) else 280)
            history['flow_front'].append(flow_front_x * 1000)
            history['max_shear_rate'].append(np.max(self.shear_rate))
        return history

    def get_results(self) -> Dict:
        return {
            'pressure_field': self.pressure.tolist(),
            'temperature_field': np.mean(self.temperature, axis=2).tolist(),
            'fill_fraction': self.fill_fraction.tolist(),
            'viscosity_field': self.viscosity.tolist(),
            'max_pressure_MPa': float(np.max(self.pressure) / 1e6),
            'avg_melt_temp_C': float(np.mean(self.temperature[self.fill_fraction > 0.5]) - 273.15),
            'max_shear_rate': float(np.max(self.shear_rate)),
        }


def run_3d_flow_analysis(geometry: MoldGeometry, resin: ResinProperties,
                          conditions: ProcessConditions) -> Dict:
    """Simplified 3D flow analysis using finite difference method."""
    nx, ny, nz = 20, 10, 5
    dx = geometry.length / nx
    dy = geometry.width / ny
    dz = geometry.thickness / nz
    u = np.zeros((nx, ny, nz))
    T = np.full((nx, ny, nz), conditions.melt_temperature)
    P = np.zeros((nx, ny))
    for k in range(nz):
        z_norm = (k - nz / 2) / (nz / 2)
        u[:, :, k] = conditions.injection_speed * (1 - z_norm ** 2)
    for i in range(nx):
        P[i, :] = conditions.injection_pressure * (1 - i / nx)
    n_time_steps = 50
    dt = conditions.injection_time / n_time_steps
    alpha = resin.thermal_conductivity / (resin.density * resin.specific_heat)
    for step in range(n_time_steps):
        T_new = T.copy()
        for i in range(1, nx - 1):
            for j in range(1, ny - 1):
                for k in range(1, nz - 1):
                    dT_adv = -u[i, j, k] * (T[i, j, k] - T[i - 1, j, k]) / dx * dt
                    dT_diff = alpha * (T[i, j, k - 1] - 2 * T[i, j, k] + T[i, j, k + 1]) / dz ** 2 * dt
                    T_new[i, j, k] = T[i, j, k] + dT_adv + dT_diff
        T_new[:, :, 0] = conditions.mold_temperature
        T_new[:, :, -1] = conditions.mold_temperature
        T_new[0, :, :] = conditions.melt_temperature
        T = T_new
    return {
        'velocity_field_shape': list(u.shape),
        'max_velocity_m_s': float(np.max(u)),
        'pressure_drop_MPa': float(conditions.injection_pressure / 1e6),
        'temperature_range_C': [float(np.min(T) - 273.15), float(np.max(T) - 273.15)],
        'centerline_temp_C': float(np.mean(T[:, ny // 2, nz // 2]) - 273.15),
        'wall_temp_C': float(np.mean(T[:, :, 0]) - 273.15),
        'velocity_profile': u[nx // 2, ny // 2, :].tolist(),
        'temperature_profile_z': T[nx // 2, ny // 2, :].tolist(),
        'pressure_along_length': P[:, ny // 2].tolist(),
    }
