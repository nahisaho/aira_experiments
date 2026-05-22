"""
Surrogate Models for HEA Property Prediction
Gaussian Process Regression with compositional kernels.
Targets: yield strength (MPa), elongation (%), pitting corrosion potential (V_SCE)
"""

import numpy as np
import torch
import gpytorch
from typing import Dict, List, Tuple, Optional
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, KFold
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, RBF, WhiteKernel, ConstantKernel
import warnings
warnings.filterwarnings("ignore")


# -----------------------------------------------------------------------
# Physics-informed property models (semi-empirical for data generation)
# -----------------------------------------------------------------------
class HEAPropertySimulator:
    """
    Semi-empirical property simulator based on descriptor-property relationships
    established in literature. Used to generate synthetic 'DFT+experiment' data
    for training surrogate models.

    References:
    - Li et al. (2019) NPJ Comput. Mater. — strength-descriptor correlations
    - Senkov et al. (2018) Nature Rev. Mater. — refractory HEA properties
    - Peng et al. (2020) Acta Mater. — ML prediction of HEA hardness
    """

    def __init__(self, noise_level: float = 0.05, random_state: int = 42):
        self.noise_level = noise_level
        self.rng = np.random.default_rng(random_state)

    def yield_strength(self, desc: Dict) -> float:
        """
        σ_y estimate [MPa]:
        Combines lattice distortion hardening (Labusch), solid solution
        strengthening, and Hall-Petch grain-boundary contribution (fixed grain size).

        σ_y ≈ A * G * δ^(4/3) * c^(2/3) + B * VEC_factor + C
        """
        G = desc["G_Voigt"]
        delta = desc["delta_r"]
        vec = desc["VEC"]
        dH = abs(desc["dH_mix"])
        Tm = desc["Tm_mean"]

        # Solid solution strengthening (Labusch-type)
        ss_hard = 0.04 * G * (delta / 100) ** (4/3)

        # VEC-based structure factor (BCC ≈ harder than FCC at low T)
        vec_factor = 1.4 if vec < 6.87 else (1.0 if vec > 8.0 else 1.2)

        # Chemical ordering contribution (proportional to |ΔH_mix|)
        order_hard = 2.5 * dH

        # Melting-point scaled strength (refractory contribution)
        Tm_factor = (Tm / 1728) ** 0.5  # normalized to pure Ni

        base = 350 * vec_factor * Tm_factor + ss_hard + order_hard
        noise = self.rng.normal(0, base * self.noise_level)
        return max(50, base + noise)

    def elongation(self, desc: Dict) -> float:
        """
        Elongation to fracture [%]:
        FCC alloys (VEC > 8) → high ductility; BCC (VEC < 6.87) → lower ductility.
        Large δ reduces ductility (lattice distortion → dislocation pinning).
        """
        vec = desc["VEC"]
        delta = desc["delta_r"]
        dH = abs(desc["dH_mix"])

        if vec > 8.0:
            base_el = 55.0
        elif vec > 6.87:
            base_el = 35.0
        else:
            base_el = 18.0

        # Lattice distortion penalty
        delta_pen = max(0, delta - 3.0) * 2.5
        # Ordering penalty
        order_pen = max(0, dH - 5.0) * 0.8

        el = max(0.5, base_el - delta_pen - order_pen)
        noise = self.rng.normal(0, el * self.noise_level)
        return max(0.5, el + noise)

    def pitting_potential(self, desc: Dict) -> float:
        """
        Pitting corrosion potential E_pit [V_SCE] in 3.5 wt% NaCl:
        Cr and Mo strongly enhance passivity; Mn is detrimental.
        Based on empirical fit from: Shi et al. (2017) Corrosion Sci.
        """
        def _get(key):
            v = desc.get(key, 0.0)
            return 0.0 if (v is None or (isinstance(v, float) and np.isnan(v))) else float(v)

        x_Cr = _get("x_Cr")
        x_Mo = _get("x_Mo")
        x_Mn = _get("x_Mn")
        x_Al = _get("x_Al")
        x_Ni = _get("x_Ni")

        E_pit = (-0.30                          # base (austenitic base)
                 + 1.80 * x_Cr                 # Cr passivation
                 + 3.00 * x_Mo                 # Mo pitting resistance
                 - 0.90 * x_Mn                 # Mn sulfide inclusion
                 + 0.40 * x_Ni                 # Ni stabilization
                 - 0.30 * x_Al)                # Al oxide formation (ambiguous)

        noise = self.rng.normal(0, 0.05)
        return E_pit + noise

    def simulate_dataset(self, descriptors_df) -> "pd.DataFrame":
        """Augment descriptor DataFrame with simulated property values."""
        import pandas as pd
        df = descriptors_df.copy()
        props = [self.yield_strength(row.to_dict()) for _, row in df.iterrows()]
        elos = [self.elongation(row.to_dict()) for _, row in df.iterrows()]
        epits = [self.pitting_potential(row.to_dict()) for _, row in df.iterrows()]

        df["yield_strength_MPa"] = props
        df["elongation_pct"] = elos
        df["pitting_potential_V"] = epits
        return df


# -----------------------------------------------------------------------
# Gaussian Process Surrogate using scikit-learn
# -----------------------------------------------------------------------
DESCRIPTOR_COLS = [
    "delta_r", "VEC", "dS_mix", "dH_mix", "Omega",
    "delta_chi", "Tm_mean", "density", "B_Voigt", "G_Voigt",
]

PROPERTY_COLS = ["yield_strength_MPa", "elongation_pct", "pitting_potential_V"]


class HEASurrogateModel:
    """
    Multi-output surrogate models (one GP per property).
    Kernel: Matern(ν=5/2) * ConstantKernel + WhiteKernel (noise).
    """

    def __init__(self, features: List[str] = None, random_state: int = 42):
        self.features = features or DESCRIPTOR_COLS
        self.random_state = random_state
        self.scalers: Dict[str, StandardScaler] = {}
        self.models: Dict[str, GaussianProcessRegressor] = {}
        self._feature_scaler = StandardScaler()
        self._fitted = False

    def _make_kernel(self):
        return (ConstantKernel(1.0, (0.01, 100))
                * Matern(length_scale=np.ones(len(self.features)),
                         length_scale_bounds=(0.01, 100),
                         nu=2.5)
                + WhiteKernel(noise_level=0.1, noise_level_bounds=(1e-4, 1.0)))

    def fit(self, X_df, y_df, n_restarts: int = 5):
        """Fit one GP per property."""
        X = X_df[self.features].values
        X_scaled = self._feature_scaler.fit_transform(X)

        for prop in PROPERTY_COLS:
            y = y_df[prop].values
            scaler = StandardScaler()
            y_scaled = scaler.fit_transform(y.reshape(-1, 1)).ravel()
            self.scalers[prop] = scaler

            gp = GaussianProcessRegressor(
                kernel=self._make_kernel(),
                n_restarts_optimizer=n_restarts,
                normalize_y=False,
                random_state=self.random_state
            )
            gp.fit(X_scaled, y_scaled)
            self.models[prop] = gp

        self._fitted = True
        return self

    def predict(self, X_df, return_std: bool = True):
        """Predict all properties with uncertainty."""
        X = X_df[self.features].values if hasattr(X_df, "columns") else X_df
        X_scaled = self._feature_scaler.transform(X)

        predictions = {}
        for prop, gp in self.models.items():
            if return_std:
                mu_s, std_s = gp.predict(X_scaled, return_std=True)
                mu = self.scalers[prop].inverse_transform(mu_s.reshape(-1, 1)).ravel()
                std = std_s * self.scalers[prop].scale_[0]
                predictions[prop] = (mu, std)
            else:
                mu_s = gp.predict(X_scaled, return_std=False)
                mu = self.scalers[prop].inverse_transform(mu_s.reshape(-1, 1)).ravel()
                predictions[prop] = mu

        return predictions

    def cross_validate(self, X_df, y_df, cv: int = 5) -> Dict[str, float]:
        """Return mean R² per property via k-fold CV."""
        X = X_df[self.features].values
        X_scaled = self._feature_scaler.fit_transform(X)
        cv_scores = {}

        for prop in PROPERTY_COLS:
            y = y_df[prop].values
            gp = GaussianProcessRegressor(
                kernel=self._make_kernel(),
                n_restarts_optimizer=3,
                normalize_y=True,
                random_state=self.random_state
            )
            scores = cross_val_score(gp, X_scaled, y,
                                     cv=KFold(cv, shuffle=True, random_state=42),
                                     scoring="r2")
            cv_scores[prop] = float(np.mean(scores))

        return cv_scores
