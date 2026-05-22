"""
Module 5: Real-time Sensor Data Assimilation (Ensemble Kalman Filter)
"""

import numpy as np
from typing import Dict
from dataclasses import dataclass


@dataclass
class SensorConfig:
    n_pressure_sensors: int = 3
    pressure_noise_std: float = 0.5
    pressure_locations: list = None
    n_temp_sensors: int = 4
    temp_noise_std: float = 1.0
    temp_locations: list = None
    displacement_noise_std: float = 0.01

    def __post_init__(self):
        if self.pressure_locations is None:
            self.pressure_locations = [0.2, 0.5, 0.8]
        if self.temp_locations is None:
            self.temp_locations = [0.1, 0.3, 0.6, 0.9]


class EnsembleKalmanFilter:
    def __init__(self, n_ensemble: int = 50, n_state: int = 8, n_obs: int = 7):
        self.Ne = n_ensemble
        self.Ns = n_state
        self.No = n_obs
        self.state_names = [
            'injection_pressure_correction', 'packing_pressure_correction',
            'melt_temp_correction', 'mold_temp_correction',
            'viscosity_correction_factor', 'htc_correction_factor',
            'shrinkage_correction_factor', 'stress_correction_factor',
        ]
        self.ensemble = np.zeros((self.Ne, self.Ns))
        self._initialize_ensemble()
        self.history = {'time': [], 'state_mean': [], 'state_std': [], 'innovation': [], 'rmse': []}

    def _initialize_ensemble(self):
        prior_mean = np.array([0, 0, 0, 0, 1.0, 1.0, 1.0, 1.0])
        prior_std = np.array([5, 3, 5, 3, 0.2, 0.2, 0.15, 0.15])
        for i in range(self.Ne):
            self.ensemble[i] = prior_mean + prior_std * np.random.randn(self.Ns)
            self.ensemble[i, 4:] = np.maximum(self.ensemble[i, 4:], 0.3)

    def forward_model(self, state: np.ndarray, nominal_params: np.ndarray) -> np.ndarray:
        P_inj = nominal_params[0] + state[0]
        T_melt = nominal_params[2] + state[2]
        T_mold = nominal_params[3] + state[3]
        visc_corr = state[4]
        htc_corr = state[5]
        P_sensors = np.zeros(3)
        for k, x_loc in enumerate([0.2, 0.5, 0.8]):
            P_sensors[k] = max(0, P_inj - P_inj * (1 - x_loc) * visc_corr)
        T_sensors = np.zeros(4)
        for k, x_loc in enumerate([0.1, 0.3, 0.6, 0.9]):
            T_sensors[k] = T_melt - (T_melt - T_mold) * x_loc * 0.3 * htc_corr
        return np.concatenate([P_sensors, T_sensors])

    def update(self, observations: np.ndarray, obs_noise: np.ndarray,
               nominal_params: np.ndarray, time: float) -> Dict:
        R = np.diag(obs_noise ** 2)
        Y_pred = np.zeros((self.Ne, self.No))
        for i in range(self.Ne):
            Y_pred[i] = self.forward_model(self.ensemble[i], nominal_params)
        x_mean = np.mean(self.ensemble, axis=0)
        y_mean = np.mean(Y_pred, axis=0)
        X_anom = self.ensemble - x_mean
        Y_anom = Y_pred - y_mean
        Pxy = X_anom.T @ Y_anom / (self.Ne - 1)
        Pyy = Y_anom.T @ Y_anom / (self.Ne - 1) + R
        K = Pxy @ np.linalg.inv(Pyy)
        innovation = observations - y_mean
        for i in range(self.Ne):
            obs_perturbed = observations + obs_noise * np.random.randn(self.No)
            self.ensemble[i] += K @ (obs_perturbed - Y_pred[i])
            self.ensemble[i, 4:] = np.clip(self.ensemble[i, 4:], 0.3, 2.0)
        x_updated = np.mean(self.ensemble, axis=0)
        x_std = np.std(self.ensemble, axis=0)
        rmse = np.sqrt(np.mean(innovation ** 2))
        self.history['time'].append(time)
        self.history['state_mean'].append(x_updated.tolist())
        self.history['state_std'].append(x_std.tolist())
        self.history['innovation'].append(innovation.tolist())
        self.history['rmse'].append(float(rmse))
        return {
            'state_mean': {self.state_names[k]: float(x_updated[k]) for k in range(self.Ns)},
            'state_std': {self.state_names[k]: float(x_std[k]) for k in range(self.Ns)},
            'rmse': float(rmse),
        }

    def run_assimilation(self, n_cycles: int = 20, nominal_params: np.ndarray = None) -> Dict:
        if nominal_params is None:
            nominal_params = np.array([80.0, 50.0, 280.0, 80.0, 50.0])
        true_state = np.array([3.0, -2.0, 5.0, -3.0, 1.15, 0.90, 1.05, 0.95])
        sensor_noise = np.array([0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0])
        for cycle in range(n_cycles):
            true_obs = self.forward_model(true_state, nominal_params)
            noisy_obs = true_obs + sensor_noise * np.random.randn(self.No)
            self.update(noisy_obs, sensor_noise, nominal_params, cycle * 1.0)
        final_mean = np.array(self.history['state_mean'][-1])
        final_std = np.array(self.history['state_std'][-1])
        error = np.abs(final_mean - true_state)
        return {
            'n_cycles': n_cycles,
            'true_state': {self.state_names[k]: float(true_state[k]) for k in range(self.Ns)},
            'estimated_state': {self.state_names[k]: float(final_mean[k]) for k in range(self.Ns)},
            'estimation_error': {self.state_names[k]: float(error[k]) for k in range(self.Ns)},
            'uncertainty': {self.state_names[k]: float(final_std[k]) for k in range(self.Ns)},
            'rmse_history': self.history['rmse'],
            'convergence_achieved': bool(np.all(error < 2 * final_std + 1.0)),
            'final_rmse': float(self.history['rmse'][-1]),
            'state_mean_history': self.history['state_mean'],
            'state_std_history': self.history['state_std'],
        }
