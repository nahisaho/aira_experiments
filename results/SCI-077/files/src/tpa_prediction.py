"""
Module 3: テクスチャプロファイル分析（TPA）パラメータの予測モデル
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score


@dataclass
class TPAResult:
    hardness: float
    cohesiveness: float
    springiness: float
    gumminess: float
    chewiness: float
    resilience: float
    adhesiveness: float
    fracturability: Optional[float] = None


def compute_tpa_from_curve(displacement, force):
    n = len(force)
    mid = n // 2
    force_1, disp_1 = force[:mid], displacement[:mid]
    hardness = np.max(force_1)
    A1 = np.trapz(np.maximum(force_1, 0), disp_1)
    force_2, disp_2 = force[mid:], displacement[mid:]
    A2 = np.trapz(np.maximum(force_2, 0), disp_2 - disp_2[0])
    cohesiveness = A2 / A1 if A1 > 0 else 0
    peak1_idx = np.argmax(force_1)
    peak2_idx = np.argmax(force_2)
    d1 = disp_1[peak1_idx]
    d2 = (disp_2[peak2_idx] - disp_2[0]) if len(disp_2) > 0 else 0
    springiness = d2 / d1 if d1 > 0 else 0
    gumminess = hardness * cohesiveness
    chewiness = gumminess * springiness
    neg_force = np.minimum(force, 0)
    adhesiveness = np.trapz(neg_force, displacement)
    A1_up = np.trapz(np.maximum(force_1[:peak1_idx+1], 0), disp_1[:peak1_idx+1])
    A1_down = np.trapz(np.maximum(force_1[peak1_idx:], 0), disp_1[peak1_idx:])
    resilience = A1_down / A1_up if A1_up > 0 else 0
    return TPAResult(hardness=hardness, cohesiveness=cohesiveness,
                     springiness=springiness, gumminess=gumminess,
                     chewiness=chewiness, resilience=resilience,
                     adhesiveness=adhesiveness)


def simulate_tpa_curve(G_inf, G_elements, tau_elements, sample_height=10.0,
                       compression_ratio=0.5, speed=1.0, rest_time=2.0,
                       n_points=200, contact_area=100.0):
    max_disp = sample_height * compression_ratio
    t_compress = max_disp / speed
    t_total = 4 * t_compress + rest_time
    t = np.linspace(0, t_total, n_points)

    def get_strain(ti):
        if ti < t_compress:
            return (ti / t_compress) * compression_ratio
        elif ti < 2 * t_compress:
            return compression_ratio * (2 - ti / t_compress)
        elif ti < 2 * t_compress + rest_time:
            return 0.0
        elif ti < 3 * t_compress + rest_time:
            return ((ti - 2*t_compress - rest_time) / t_compress) * compression_ratio
        elif ti < 4 * t_compress + rest_time:
            return compression_ratio * (4 + rest_time/t_compress - ti/t_compress)
        return 0.0

    strain_history = np.array([get_strain(ti) for ti in t])
    dt_val = t[1] - t[0]
    strain_rate = np.gradient(strain_history, dt_val)
    displacement = strain_history * sample_height
    force = np.zeros(n_points)

    for i in range(n_points):
        sigma = 0.0
        for j in range(i + 1):
            t_diff = t[i] - t[j]
            G_relax = G_inf
            for Gk, tauk in zip(G_elements, tau_elements):
                G_relax += Gk * np.exp(-t_diff / tauk)
            sigma += G_relax * strain_rate[j] * dt_val
        force[i] = sigma * contact_area * 1e-6
    return displacement, force


@dataclass
class TPAPredictionModel:
    models: Dict[str, GradientBoostingRegressor] = field(default_factory=dict)
    scalers: Dict[str, StandardScaler] = field(default_factory=dict)
    feature_names: List[str] = field(default_factory=lambda: [
        'protein', 'fat', 'carbohydrate', 'moisture', 'salt',
        'pH', 'temperature', 'heating_time', 'cooling_rate'
    ])
    targets: List[str] = field(default_factory=lambda: [
        'hardness', 'cohesiveness', 'springiness',
        'gumminess', 'chewiness', 'resilience'
    ])

    def generate_training_data(self, n_samples=500, seed=42):
        rng = np.random.default_rng(seed)
        X = np.zeros((n_samples, len(self.feature_names)))
        X[:, 0] = rng.uniform(5, 30, n_samples)
        X[:, 1] = rng.uniform(1, 40, n_samples)
        X[:, 2] = rng.uniform(5, 50, n_samples)
        X[:, 3] = rng.uniform(20, 80, n_samples)
        X[:, 4] = rng.uniform(0, 3, n_samples)
        X[:, 5] = rng.uniform(3, 8, n_samples)
        X[:, 6] = rng.uniform(60, 200, n_samples)
        X[:, 7] = rng.uniform(1, 120, n_samples)
        X[:, 8] = rng.uniform(0.5, 20, n_samples)
        y = {}
        protein, fat, carb, moisture = X[:,0], X[:,1], X[:,2], X[:,3]
        t_heat, cool_rate, pH = X[:,7], X[:,8], X[:,5]
        y['hardness'] = 50*protein**0.8*(1+0.1*carb)*np.exp(-0.02*moisture)*np.exp(-0.01*fat)*(1+0.3*np.log1p(t_heat))*(1+0.01*cool_rate)*np.abs(pH-5.0)**0.3*(1+rng.normal(0,0.05,n_samples))
        y['cohesiveness'] = np.clip(0.3+0.01*fat+0.005*protein-0.001*(moisture-50)**2/50+rng.normal(0,0.02,n_samples),0.1,0.95)
        y['springiness'] = np.clip(0.4+0.01*protein+0.005*carb-0.003*moisture+0.01*cool_rate+rng.normal(0,0.03,n_samples),0.2,0.99)
        y['gumminess'] = y['hardness'] * y['cohesiveness']
        y['chewiness'] = y['gumminess'] * y['springiness']
        y['resilience'] = np.clip(0.2+0.005*protein+0.003*carb-0.002*moisture+rng.normal(0,0.02,n_samples),0.05,0.8)
        return X, y

    def train(self, X, y):
        scores = {}
        for target in self.targets:
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            model = GradientBoostingRegressor(n_estimators=200, max_depth=5,
                                              learning_rate=0.05, min_samples_leaf=5, random_state=42)
            cv_scores = cross_val_score(model, X_scaled, y[target], cv=5, scoring='r2')
            model.fit(X_scaled, y[target])
            self.models[target] = model
            self.scalers[target] = scaler
            scores[target] = {'cv_r2_mean': cv_scores.mean(), 'cv_r2_std': cv_scores.std(),
                              'feature_importance': dict(zip(self.feature_names, model.feature_importances_))}
        return scores

    def predict(self, composition):
        X = np.array([[composition.get(f, 0) for f in self.feature_names]])
        preds = {}
        for target in self.targets:
            X_scaled = self.scalers[target].transform(X)
            preds[target] = float(self.models[target].predict(X_scaled)[0])
        return TPAResult(hardness=max(preds['hardness'],0),
                         cohesiveness=np.clip(preds['cohesiveness'],0,1),
                         springiness=np.clip(preds['springiness'],0,1),
                         gumminess=max(preds['gumminess'],0),
                         chewiness=max(preds['chewiness'],0),
                         resilience=np.clip(preds['resilience'],0,1),
                         adhesiveness=0.0)
