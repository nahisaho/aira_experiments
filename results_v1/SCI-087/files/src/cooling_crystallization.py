"""
Module 2: Cooling & Solidification Process Modeling
- Heat transfer during cooling phase
- Nakamura crystallization kinetics model
- Crystallinity-dependent material properties
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, Tuple


@dataclass
class CrystallizationParams:
    """Nakamura/Avrami crystallization kinetics parameters for PA66."""
    X_inf: float = 0.35
    n_avrami: float = 2.5
    K0: float = 1.0e6
    Kg: float = 3.5e5
    Tm: float = 533.15
    Tg: float = 323.15
    U_star: float = 6284.0
    T_inf: float = 273.15
    delta_Hf: float = 1.9e5


@dataclass
class CoolingGeometry:
    """Cooling channel layout for mold."""
    part_thickness: float = 0.003
    channel_diameter: float = 0.010
    channel_depth: float = 0.015
    channel_pitch: float = 0.025
    coolant_temperature: float = 323.15
    coolant_flow_rate: float = 10.0
    n_nodes_z: int = 41


def nakamura_rate(T: float, X: float, params: CrystallizationParams) -> float:
    """Nakamura crystallization rate."""
    if X >= params.X_inf * 0.999 or T >= params.Tm or T <= params.T_inf + 10:
        return 0.0
    X_rel = max(X / params.X_inf, 0.001)
    if X_rel >= 0.999:
        return 0.0
    R = 8.314
    dT = params.Tm - T
    if dT <= 0:
        return 0.0
    K_T = params.K0 * np.exp(-params.U_star / (R * (T - params.T_inf))) * \
          np.exp(-params.Kg / (T * dT))
    n = params.n_avrami
    term = max(-np.log(1.0 - X_rel), 1e-10)
    rate = n * K_T * (1 - X_rel) * term ** ((n - 1) / n)
    return rate * params.X_inf


class CoolingSimulator:
    """1D through-thickness cooling simulation with crystallization."""

    def __init__(self, geom: CoolingGeometry, cryst: CrystallizationParams,
                 initial_temp: float = 553.15):
        self.geom = geom
        self.cryst = cryst
        self.nz = geom.n_nodes_z
        self.dz = geom.part_thickness / (self.nz - 1)
        self.rho = 1350.0
        self.cp_amorphous = 1700.0
        self.cp_crystalline = 1500.0
        self.k_amorphous = 0.24
        self.k_crystalline = 0.38
        self.T = np.full(self.nz, initial_temp)
        self.X = np.zeros(self.nz)

    def effective_properties(self, T: float, X: float) -> Tuple[float, float]:
        X_rel = X / self.cryst.X_inf if self.cryst.X_inf > 0 else 0
        cp = self.cp_amorphous * (1 - X_rel) + self.cp_crystalline * X_rel
        k = self.k_amorphous * (1 - X_rel) + self.k_crystalline * X_rel
        return cp, k

    def compute_htc(self) -> float:
        k_water = 0.65
        Re = 5000
        Pr = 5.0
        Nu = 0.023 * Re ** 0.8 * Pr ** 0.4
        h_coolant = Nu * k_water / self.geom.channel_diameter
        k_steel = 40.0
        R_mold = self.geom.channel_depth / k_steel
        R_coolant = 1.0 / h_coolant
        return 1.0 / (R_mold + R_coolant)

    def simulate(self, total_time: float = 30.0, dt: float = 0.01) -> Dict:
        n_steps = int(total_time / dt)
        h_wall = self.compute_htc()
        record_interval = max(1, n_steps // 200)
        history = {
            'time': [], 'center_temp_C': [], 'surface_temp_C': [],
            'avg_crystallinity': [], 'center_crystallinity': [],
            'ejection_ready': False, 'ejection_time': None,
            'solidification_front': []
        }
        for step in range(n_steps):
            t = step * dt
            T_new = self.T.copy()
            X_new = self.X.copy()
            for k in range(self.nz):
                cp, kk = self.effective_properties(self.T[k], self.X[k])
                alpha = kk / (self.rho * cp)
                if 0 < k < self.nz - 1:
                    d2T = (self.T[k + 1] - 2 * self.T[k] + self.T[k - 1]) / self.dz ** 2
                elif k == 0:
                    d2T = 2 * (self.T[1] - self.T[0] -
                               h_wall * self.dz / kk * (self.T[0] - self.geom.coolant_temperature)) / self.dz ** 2
                else:
                    d2T = 2 * (self.T[-2] - self.T[-1] -
                               h_wall * self.dz / kk * (self.T[-1] - self.geom.coolant_temperature)) / self.dz ** 2
                dXdt = nakamura_rate(self.T[k], self.X[k], self.cryst)
                X_new[k] = min(self.X[k] + dXdt * dt, self.cryst.X_inf)
                dT = alpha * d2T * dt + self.cryst.delta_Hf / cp * dXdt * dt
                T_new[k] = self.T[k] + dT
            self.T = T_new
            self.X = X_new
            solid_fraction = np.sum(self.T < self.cryst.Tg + 50) / self.nz
            if step % record_interval == 0:
                history['time'].append(t)
                history['center_temp_C'].append(float(self.T[self.nz // 2] - 273.15))
                history['surface_temp_C'].append(float(self.T[0] - 273.15))
                history['avg_crystallinity'].append(float(np.mean(self.X)))
                history['center_crystallinity'].append(float(self.X[self.nz // 2]))
                history['solidification_front'].append(float(solid_fraction))
            ejection_temp = self.cryst.Tg + 30
            if self.T[self.nz // 2] < ejection_temp and not history['ejection_ready']:
                history['ejection_ready'] = True
                history['ejection_time'] = t
        z_positions = np.linspace(0, self.geom.part_thickness * 1000, self.nz)
        history['final_z_mm'] = z_positions.tolist()
        history['final_temperature_C'] = (self.T - 273.15).tolist()
        history['final_crystallinity'] = self.X.tolist()
        return history

    def get_summary(self) -> Dict:
        return {
            'final_center_temp_C': float(self.T[self.nz // 2] - 273.15),
            'final_surface_temp_C': float(self.T[0] - 273.15),
            'avg_crystallinity_pct': float(np.mean(self.X) * 100),
            'max_crystallinity_pct': float(np.max(self.X) * 100),
            'surface_crystallinity_pct': float(self.X[0] * 100),
            'center_crystallinity_pct': float(self.X[self.nz // 2] * 100),
            'crystallinity_gradient': float((self.X[0] - self.X[self.nz // 2]) * 100),
        }
