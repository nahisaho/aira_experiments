"""
Module 4: Process Parameter – Quality Relationship Model
"""

import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
from scipy.optimize import differential_evolution
from typing import Dict


class ProcessQualityModel:
    def __init__(self):
        self.param_names = [
            'injection_pressure_MPa', 'packing_pressure_MPa',
            'cooling_time_s', 'melt_temp_C', 'mold_temp_C',
            'injection_speed_mm_s'
        ]
        self.quality_names = [
            'warpage_mm', 'sink_depth_mm', 'weight_g',
            'shrinkage_pct', 'residual_stress_MPa'
        ]
        self.param_bounds = np.array([
            [50, 120], [30, 80], [10, 40],
            [260, 300], [60, 100], [30, 80],
        ])
        self.scaler_X = StandardScaler()
        self.scaler_y = {}
        self.models = {}
        self.X_train = None
        self.y_train = {}

    def generate_doe(self, n_samples: int = 80) -> np.ndarray:
        n_params = len(self.param_names)
        samples = np.zeros((n_samples, n_params))
        for j in range(n_params):
            points = np.linspace(0, 1, n_samples + 1)
            lower, upper = points[:-1], points[1:]
            samples[:, j] = np.random.uniform(lower, upper)
            np.random.shuffle(samples[:, j])
        for j in range(n_params):
            samples[:, j] = self.param_bounds[j, 0] + \
                            samples[:, j] * (self.param_bounds[j, 1] - self.param_bounds[j, 0])
        return samples

    def physics_model(self, params: np.ndarray) -> Dict[str, float]:
        P_inj, P_pack, t_cool, T_melt, T_mold, v_inj = params
        dT = T_melt - T_mold
        warpage = 0.05 * (dT / 200) ** 1.5 * (1 - P_pack / 100) * \
                  np.exp(-t_cool / 30) * (1 + 0.3 * np.random.randn())
        warpage = max(0.01, abs(warpage))
        sink = 0.02 * (1 - P_pack / P_inj) * np.exp(-t_cool / 25) * \
               (1 + T_mold / 300) * (1 + 0.2 * np.random.randn())
        sink = max(0.001, abs(sink))
        pvT_shrinkage = 0.005 * (T_melt / 280) * (80 / max(P_pack, 30))
        weight = 45.0 * (1 - pvT_shrinkage) * (1 + 0.001 * np.random.randn())
        shrinkage = 0.8 * (1 - P_pack / 100) * (dT / 200) * np.exp(-t_cool / 35) + 0.3
        shrinkage = max(0.1, shrinkage + 0.05 * np.random.randn())
        stress = 5.0 * (P_inj / 80) * (dT / 200) * (1 - t_cool / 50) + 3.0 * (v_inj / 50)
        stress = max(1.0, stress + 0.5 * np.random.randn())
        return {
            'warpage_mm': warpage, 'sink_depth_mm': sink, 'weight_g': weight,
            'shrinkage_pct': shrinkage, 'residual_stress_MPa': stress,
        }

    def train(self, n_samples: int = 80) -> Dict:
        np.random.seed(42)
        X = self.generate_doe(n_samples)
        self.X_train = X
        y_data = {name: [] for name in self.quality_names}
        for i in range(n_samples):
            result = self.physics_model(X[i])
            for name in self.quality_names:
                y_data[name].append(result[name])
        X_scaled = self.scaler_X.fit_transform(X)
        cv_scores = {}
        for name in self.quality_names:
            y = np.array(y_data[name])
            self.y_train[name] = y
            scaler_y = StandardScaler()
            y_scaled = scaler_y.fit_transform(y.reshape(-1, 1)).ravel()
            self.scaler_y[name] = scaler_y
            kernel = ConstantKernel(1.0) * Matern(
                length_scale=np.ones(len(self.param_names)), nu=2.5)
            gp = GaussianProcessRegressor(
                kernel=kernel, n_restarts_optimizer=5, alpha=0.01, normalize_y=False)
            gp.fit(X_scaled, y_scaled)
            self.models[name] = gp
            scores = cross_val_score(gp, X_scaled, y_scaled, cv=5, scoring='r2')
            cv_scores[name] = {'r2_mean': float(np.mean(scores)), 'r2_std': float(np.std(scores))}
        return {'n_training_samples': n_samples, 'cv_scores': cv_scores}

    def predict(self, params: np.ndarray) -> Dict:
        X_scaled = self.scaler_X.transform(params.reshape(1, -1))
        predictions = {}
        for name in self.quality_names:
            y_pred_scaled, y_std_scaled = self.models[name].predict(X_scaled, return_std=True)
            y_pred = self.scaler_y[name].inverse_transform(y_pred_scaled.reshape(-1, 1)).ravel()[0]
            y_std = y_std_scaled[0] * self.scaler_y[name].scale_[0]
            predictions[name] = {
                'mean': float(y_pred), 'std': float(y_std),
                'ci_95': [float(y_pred - 1.96 * y_std), float(y_pred + 1.96 * y_std)]
            }
        return predictions

    def sobol_sensitivity(self, n_samples: int = 200) -> Dict:
        np.random.seed(123)
        n_params = len(self.param_names)
        A = np.random.rand(n_samples, n_params)
        B = np.random.rand(n_samples, n_params)
        for j in range(n_params):
            A[:, j] = self.param_bounds[j, 0] + A[:, j] * (self.param_bounds[j, 1] - self.param_bounds[j, 0])
            B[:, j] = self.param_bounds[j, 0] + B[:, j] * (self.param_bounds[j, 1] - self.param_bounds[j, 0])
        sensitivities = {}
        for qname in self.quality_names:
            yA = np.array([self.predict(A[i])[qname]['mean'] for i in range(n_samples)])
            yB = np.array([self.predict(B[i])[qname]['mean'] for i in range(n_samples)])
            var_total = np.var(np.concatenate([yA, yB]))
            S1 = np.zeros(n_params)
            for j in range(n_params):
                AB_j = A.copy()
                AB_j[:, j] = B[:, j]
                yAB_j = np.array([self.predict(AB_j[i])[qname]['mean'] for i in range(n_samples)])
                S1[j] = np.mean(yB * (yAB_j - yA)) / (var_total + 1e-10)
                S1[j] = max(0, min(1, S1[j]))
            s = np.sum(S1)
            if s > 0:
                S1 = S1 / s
            sensitivities[qname] = {self.param_names[j]: float(S1[j]) for j in range(n_params)}
        return sensitivities

    def optimize(self) -> Dict:
        def objective(params):
            pred = self.predict(np.array(params))
            return (pred['warpage_mm']['mean'] * 10) ** 2 + \
                   (pred['sink_depth_mm']['mean'] * 20) ** 2 + \
                   (pred['shrinkage_pct']['mean'] * 2) ** 2
        result = differential_evolution(objective, bounds=list(self.param_bounds), seed=42, maxiter=100)
        optimal_quality = self.predict(result.x)
        return {
            'optimal_parameters': {self.param_names[i]: float(result.x[i]) for i in range(len(self.param_names))},
            'predicted_quality': optimal_quality,
            'optimization_success': result.success,
        }
