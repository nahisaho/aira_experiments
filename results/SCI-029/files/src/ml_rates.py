"""
ML-Based Photochemical Rate Constant Predictor
Extends Evans-Polanyi (Bell-Evans-Polanyi) relationship with
molecular descriptors and Gaussian Process Regression.

Evans-Polanyi: log(k) = log(A) - alpha * Ea / (RT)
Extended: log(k) = f(structural features, BDE, ionization energy, ...)
"""
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error, r2_score
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings("ignore")

R_GAS   = 8.314e-3   # kJ/(mol·K)
T_ATM   = 298.15     # K


@dataclass
class MoleculeDescriptors:
    """Molecular features for rate constant prediction."""
    name: str
    n_carbons: int
    n_oxygens: int
    n_nitrogens: int
    n_double_bonds: int
    n_aromatic_rings: int
    BDE: float          # C-H bond dissociation energy [kJ/mol]
    IP: float           # ionization potential [eV]
    EA: float           # electron affinity [eV]
    dipole_moment: float  # Debye
    log_Psat: float     # log10(Psat/Pa)
    delta_H_rxn: float  # reaction enthalpy [kJ/mol]
    alpha_EP: float     # Evans-Polanyi alpha coefficient [0-1]


# ── Training dataset from literature ─────────────────────────────────────────
# OH reactions: k in cm3 molecule-1 s-1  (Atkinson 2003, NIST database)
TRAINING_DATA: List[Tuple[MoleculeDescriptors, float]] = [
    # (descriptor, log10(k_OH))
    (MoleculeDescriptors("alpha-pinene",  10, 0, 0, 1, 0, 413.0, 8.82, -0.7, 1.25, 2.80,  -88.0, 0.47), np.log10(5.33e-11)),
    (MoleculeDescriptors("beta-pinene",   10, 0, 0, 1, 0, 410.0, 8.79, -0.8, 1.20, 2.93,  -92.0, 0.47), np.log10(7.89e-11)),
    (MoleculeDescriptors("limonene",      10, 0, 0, 2, 0, 405.0, 8.61, -0.9, 1.30, 2.30,  -98.0, 0.48), np.log10(1.71e-10)),
    (MoleculeDescriptors("isoprene",       5, 0, 0, 2, 0, 403.0, 8.85, -1.0, 1.40, 4.87, -100.0, 0.49), np.log10(1.00e-10)),
    (MoleculeDescriptors("toluene",        7, 0, 0, 4, 1, 432.0, 8.83, -0.3, 0.36, 3.57,  -54.0, 0.52), np.log10(5.63e-12)),
    (MoleculeDescriptors("benzene",        6, 0, 0, 3, 1, 473.0, 9.25, -1.1, 0.00, 3.85,  -32.0, 0.55), np.log10(1.22e-12)),
    (MoleculeDescriptors("propene",        3, 0, 0, 1, 0, 410.0, 9.73, -2.0, 0.35, 5.30,  -85.0, 0.47), np.log10(4.85e-12)),
    (MoleculeDescriptors("ethene",         2, 0, 0, 1, 0, 452.0,10.51, -2.2, 0.00, 5.59,  -77.0, 0.50), np.log10(8.52e-12)),
    (MoleculeDescriptors("n-hexane",       6, 0, 0, 0, 0, 420.0,10.18, -0.1, 0.00, 4.03,  -65.0, 0.45), np.log10(5.45e-12)),
    (MoleculeDescriptors("methanol",       1, 1, 0, 0, 0, 401.0,10.85, -1.8, 1.70, 5.97,  -71.0, 0.46), np.log10(9.44e-13)),
    (MoleculeDescriptors("acetaldehyde",   2, 1, 0, 1, 0, 368.0, 9.76, -0.6, 2.74, 5.04,  -60.0, 0.43), np.log10(1.50e-11)),
    (MoleculeDescriptors("formaldehyde",   1, 1, 0, 1, 0, 363.0,10.88, -0.6, 2.33, 5.26,  -58.0, 0.42), np.log10(8.50e-12)),
    (MoleculeDescriptors("acrolein",       3, 1, 0, 2, 0, 385.0, 9.95, -0.8, 3.12, 4.85,  -72.0, 0.44), np.log10(1.99e-11)),
    (MoleculeDescriptors("pinaldehyde",   10, 1, 0, 1, 0, 370.0, 9.65, -0.5, 2.80, 0.95,  -63.0, 0.43), np.log10(3.20e-11)),
    (MoleculeDescriptors("pinonaldehyde", 10, 1, 0, 1, 0, 372.0, 9.63, -0.5, 2.75, 0.94,  -65.0, 0.43), np.log10(4.00e-11)),
    (MoleculeDescriptors("cresol",         7, 1, 0, 3, 1, 375.0, 8.17, -0.2, 1.45, 2.22,  -80.0, 0.44), np.log10(4.71e-11)),
    (MoleculeDescriptors("MVK",            4, 1, 0, 2, 0, 393.0, 9.65, -0.7, 3.05, 4.09,  -70.0, 0.45), np.log10(1.99e-11)),
    (MoleculeDescriptors("methacrolein",   4, 1, 0, 2, 0, 388.0, 9.91, -0.9, 3.47, 3.96,  -75.0, 0.45), np.log10(2.89e-11)),
    (MoleculeDescriptors("m-xylene",       8, 0, 0, 4, 1, 431.0, 8.56, -0.3, 0.37, 3.21,  -60.0, 0.52), np.log10(2.31e-11)),
    (MoleculeDescriptors("naphthalene",   10, 0, 0, 5, 2, 468.0, 8.14, -0.3, 0.00, 2.35,  -46.0, 0.54), np.log10(2.44e-11)),
]


def build_feature_matrix(data: List[Tuple[MoleculeDescriptors, float]]) -> Tuple[np.ndarray, np.ndarray]:
    """Convert descriptor list to numpy arrays."""
    X, y = [], []
    for desc, log_k in data:
        features = [
            desc.n_carbons,
            desc.n_oxygens,
            desc.n_double_bonds,
            desc.n_aromatic_rings,
            desc.BDE,
            desc.IP,
            desc.EA,
            desc.dipole_moment,
            desc.log_Psat,
            desc.delta_H_rxn,
            desc.alpha_EP,
        ]
        X.append(features)
        y.append(log_k)
    return np.array(X, dtype=float), np.array(y, dtype=float)


class EvansPolanyiGPR:
    """
    Extended Evans-Polanyi model using Gaussian Process Regression.
    Prior mean function: log(k) = log(A_EP) - alpha * delta_H_rxn / (RT)
    Kernel: RBF on molecular descriptors
    """

    def __init__(self):
        kernel = (
            ConstantKernel(1.0) *
            RBF(length_scale=np.ones(11), length_scale_bounds=(0.1, 100)) +
            WhiteKernel(noise_level=0.1)
        )
        self.gpr = GaussianProcessRegressor(
            kernel=kernel,
            n_restarts_optimizer=10,
            normalize_y=True,
            alpha=1e-6,
        )
        self.scaler = StandardScaler()
        self.ridge = Ridge(alpha=1.0)
        self.trained = False
        self.cv_scores: Dict[str, float] = {}
        self.feature_names = [
            "n_carbons", "n_oxygens", "n_double_bonds", "n_aromatic_rings",
            "BDE", "IP", "EA", "dipole_moment", "log_Psat",
            "delta_H_rxn", "alpha_EP",
        ]

    def _evans_polanyi_prior(self, X_raw: np.ndarray) -> np.ndarray:
        """Evans-Polanyi mean function: log(k) = log(A) - alpha * dH / (RT)"""
        log_A = -10.0   # pre-exponential (cm3 molecule-1 s-1)
        alpha = X_raw[:, 10]       # alpha_EP column
        dH    = X_raw[:, 9]        # delta_H_rxn column
        return log_A - alpha * dH / (R_GAS * T_ATM)

    def fit(self, X: np.ndarray, y: np.ndarray):
        # Remove Evans-Polanyi prior from target
        y_residual = y - self._evans_polanyi_prior(X)

        # Scale features
        X_scaled = self.scaler.fit_transform(X)

        # Fit GPR on residuals
        self.gpr.fit(X_scaled, y_residual)

        # Also fit linear ridge as fallback
        self.ridge.fit(X_scaled, y)

        # Cross-validation
        cv = cross_val_score(self.ridge, X_scaled, y, cv=5, scoring="r2")
        self.cv_scores = {"R2_cv_mean": float(cv.mean()), "R2_cv_std": float(cv.std())}
        self.trained = True

    def predict(self, X: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Returns (log_k_pred, log_k_std)."""
        prior = self._evans_polanyi_prior(X)
        X_scaled = self.scaler.transform(X)
        residual_mean, residual_std = self.gpr.predict(X_scaled, return_std=True)
        return prior + residual_mean, residual_std

    def predict_new_species(self, desc: MoleculeDescriptors) -> Tuple[float, float]:
        features = np.array([[
            desc.n_carbons, desc.n_oxygens, desc.n_double_bonds,
            desc.n_aromatic_rings, desc.BDE, desc.IP, desc.EA,
            desc.dipole_moment, desc.log_Psat, desc.delta_H_rxn, desc.alpha_EP,
        ]])
        log_k_mean, log_k_std = self.predict(features)
        return float(log_k_mean[0]), float(log_k_std[0])

    def get_metrics(self, X: np.ndarray, y: np.ndarray) -> Dict[str, float]:
        log_k_pred, _ = self.predict(X)
        rmse  = float(np.sqrt(mean_squared_error(y, log_k_pred)))
        r2    = float(r2_score(y, log_k_pred))
        mae   = float(np.mean(np.abs(y - log_k_pred)))
        return {"RMSE": rmse, "R2": r2, "MAE": mae, **self.cv_scores}

    def feature_importance(self, X: np.ndarray, y: np.ndarray) -> Dict[str, float]:
        """Permutation-based feature importance."""
        X_s = self.scaler.transform(X)
        base_r2 = r2_score(y, self.ridge.predict(X_s))
        importances = {}
        rng = np.random.default_rng(42)
        for i, name in enumerate(self.feature_names):
            X_perm = X_s.copy()
            X_perm[:, i] = rng.permutation(X_perm[:, i])
            perm_r2 = r2_score(y, self.ridge.predict(X_perm))
            importances[name] = float(base_r2 - perm_r2)
        return importances


def train_rate_predictor() -> Tuple[EvansPolanyiGPR, np.ndarray, np.ndarray, Dict]:
    """Train the ML model and return model + data + metrics."""
    X, y = build_feature_matrix(TRAINING_DATA)
    model = EvansPolanyiGPR()
    model.fit(X, y)
    metrics = model.get_metrics(X, y)
    importances = model.feature_importance(X, y)
    metrics["feature_importances"] = importances
    return model, X, y, metrics


# ── Novel species predictions ─────────────────────────────────────────────────
NEW_SPECIES_DESCRIPTORS = [
    MoleculeDescriptors("camphene",         10, 0, 0, 1, 0, 412.0, 8.76, -0.7, 1.15, 2.88, -89.0, 0.47),
    MoleculeDescriptors("delta-3-carene",   10, 0, 0, 1, 0, 409.0, 8.80, -0.8, 1.28, 2.75, -91.0, 0.47),
    MoleculeDescriptors("myrcene",          10, 0, 0, 3, 0, 400.0, 8.70, -1.0, 1.38, 3.12, -97.0, 0.48),
    MoleculeDescriptors("linalool",         10, 1, 0, 2, 0, 395.0, 8.55, -0.8, 1.60, 0.83, -95.0, 0.47),
    MoleculeDescriptors("p-cymene",         10, 0, 0, 3, 1, 430.0, 8.45, -0.3, 0.31, 2.35, -58.0, 0.51),
    MoleculeDescriptors("sabinene",         10, 0, 0, 1, 0, 411.0, 8.78, -0.7, 1.22, 2.81, -87.0, 0.47),
    MoleculeDescriptors("o-cresol",          7, 1, 0, 3, 1, 372.0, 8.14, -0.2, 1.45, 2.22, -82.0, 0.44),
    MoleculeDescriptors("glyoxal",           2, 2, 0, 2, 0, 355.0, 9.94, -0.5, 4.00, 5.45, -52.0, 0.42),
]
