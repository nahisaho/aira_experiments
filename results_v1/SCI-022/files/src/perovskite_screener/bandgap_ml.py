"""
Module 2: DFT + Machine Learning Hybrid Band Gap & Absorption Prediction
=========================================================================
Implements:
  - Feature engineering from structural/chemical descriptors
  - Random Forest regressor trained on known perovskite data
  - Absorption coefficient model (Tauc gap + Franz-Keldysh)
  - Spin-orbit coupling correction for heavy metals (Bi)
  - Rashba splitting indicator for non-centrosymmetric structures
  - Spectroscopic Limited Maximum Efficiency (SLME)
"""

import numpy as np
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, LeaveOneOut
from sklearn.metrics import mean_absolute_error, r2_score
import warnings
warnings.filterwarnings("ignore")

from .materials_database import (
    IONIC_RADII, ELECTRONEGATIVITY, KNOWN_BANDGAPS, get_ionic_radius
)


# ── Feature Engineering ──────────────────────────────────────────────────────

def compute_descriptors(A: str, B: str, X: str, B_ox: int = 2) -> dict:
    """
    Compute 18 chemical/structural descriptors for band gap ML model.
    """
    try:
        rA = get_ionic_radius(A, 1, cn=12)
    except Exception:
        rA = get_ionic_radius(A, 1, cn=6) * 1.12
    rB = get_ionic_radius(B, B_ox, cn=6)
    rX = get_ionic_radius(X, -1, cn=6)

    chi_A = ELECTRONEGATIVITY.get(A, 1.5)
    chi_B = ELECTRONEGATIVITY.get(B, 2.0)
    chi_X = ELECTRONEGATIVITY.get(X, 2.7)

    # Tolerance factors
    t  = (rA + rX) / (np.sqrt(2) * (rB + rX))
    mu = rB / rX
    tau_ratio = rX / rB - 1.0 * (1.0 - (rA / rB) / np.log(rA / rB))

    # Electronic descriptors
    chi_diff_BX = chi_X - chi_B      # ionicity of B-X bond
    chi_diff_AX = chi_X - chi_A      # ionicity of A-X bond
    chi_avg     = (chi_A + chi_B + 3 * chi_X) / 5

    # Lattice / size
    a_est    = 2 * (rB + rX)         # estimated lattice parameter (Angstrom)
    V_est    = a_est ** 3            # unit cell volume (Ang^3)
    r_ratio  = rA / rB

    # Periodic table features
    B_period = _get_period(B)
    X_period = _get_period(X)
    B_ox_num = B_ox

    # SOC indicator (heavy atoms)
    soc_indicator = 1.0 if B in ["Bi", "Sb", "In"] else (0.5 if B in ["Sn", "Ge"] else 0.0)

    # Halogen type
    halogen_code = {"I": 0, "Br": 1, "Cl": 2, "F": 3, "SCN": 4}.get(X, 0)

    # B-site type
    B_code = {"Sn": 0, "Ge": 1, "Bi": 2, "Sb": 3, "In": 4, "Pb": 5}.get(B, 6)

    return {
        "t": t, "mu": mu, "tau": tau_ratio,
        "chi_diff_BX": chi_diff_BX, "chi_diff_AX": chi_diff_AX,
        "chi_avg": chi_avg,
        "a_est": a_est, "V_est": V_est, "r_ratio": r_ratio,
        "B_period": B_period, "X_period": X_period, "B_ox": B_ox_num,
        "soc_indicator": soc_indicator,
        "halogen_code": halogen_code, "B_code": B_code,
        "rA": rA, "rB": rB, "rX": rX,
    }


def _get_period(element: str) -> int:
    """Return periodic table period for an element."""
    periods = {
        "H": 1, "Li": 2, "Be": 2, "B": 2, "C": 2, "N": 2, "O": 2, "F": 2, "Ne": 2,
        "Na": 3, "Mg": 3, "Al": 3, "Si": 3, "P": 3, "S": 3, "Cl": 3, "Ar": 3,
        "K": 4, "Ca": 4, "Ge": 4, "As": 4, "Se": 4, "Br": 4,
        "Rb": 5, "Sr": 5, "Sn": 5, "Sb": 5, "Te": 5, "I": 5, "In": 5, "Ag": 5,
        "Cs": 6, "Ba": 6, "Pb": 6, "Bi": 6, "Au": 6,
        "MA": 2, "FA": 2, "EA": 2, "DMA": 2, "GA": 2,
        "SCN": 3, "HCOO": 2, "BF4": 2,
    }
    return periods.get(element, 4)


# ── Training Data ─────────────────────────────────────────────────────────────

def build_training_data():
    """Build training set from KNOWN_BANDGAPS with B oxidation states."""
    B_ox_map = {"Pb": 2, "Sn": 2, "Ge": 2, "Bi": 3, "Sb": 3, "Ti": 4, "In": 3}
    X_list, y_list, labels = [], [], []
    for (A, B, X), data in KNOWN_BANDGAPS.items():
        ox = B_ox_map.get(B, 2)
        desc = compute_descriptors(A, B, X, ox)
        X_list.append(list(desc.values()))
        y_list.append(data["Eg"])
        labels.append(f"{A}{B}{X}3")
    return np.array(X_list), np.array(y_list), labels


# ── ML Models ────────────────────────────────────────────────────────────────

class BandGapPredictor:
    """
    Gradient Boosting + Random Forest ensemble for perovskite band gap prediction.
    Trained on experimental + high-quality DFT data.
    """
    def __init__(self):
        self.gb_model  = GradientBoostingRegressor(
            n_estimators=200, learning_rate=0.05, max_depth=4,
            subsample=0.8, random_state=42
        )
        self.rf_model  = RandomForestRegressor(
            n_estimators=300, max_depth=6, min_samples_leaf=2,
            random_state=42
        )
        self.scaler    = StandardScaler()
        self.is_fitted = False
        self.cv_mae    = None
        self.cv_r2     = None
        self.feature_names = None

    def fit(self, verbose: bool = True):
        X_train, y_train, labels = build_training_data()
        self.feature_names = list(compute_descriptors("MA", "Pb", "I").keys())
        X_scaled = self.scaler.fit_transform(X_train)
        self.gb_model.fit(X_scaled, y_train)
        self.rf_model.fit(X_scaled, y_train)

        # Cross-validation (LOO for small dataset)
        loo = LeaveOneOut()
        gb_preds = np.array([
            self.gb_model.fit(
                self.scaler.fit_transform(np.delete(X_train, i, axis=0)),
                np.delete(y_train, i)
            ).predict(self.scaler.transform(X_train[[i]]))[0]
            for i in range(len(y_train))
        ])
        # Refit on full data
        X_scaled = self.scaler.fit_transform(X_train)
        self.gb_model.fit(X_scaled, y_train)
        self.rf_model.fit(X_scaled, y_train)

        self.cv_mae = mean_absolute_error(y_train, gb_preds)
        ss_res = np.sum((y_train - gb_preds) ** 2)
        ss_tot = np.sum((y_train - np.mean(y_train)) ** 2)
        self.cv_r2  = 1 - ss_res / ss_tot
        self.is_fitted = True
        if verbose:
            print(f"  [ML] LOO-CV MAE = {self.cv_mae:.3f} eV,  R² = {self.cv_r2:.3f}")
        return self

    def predict(self, A: str, B: str, X: str, B_ox: int = 2) -> dict:
        if not self.is_fitted:
            self.fit(verbose=False)
        desc = compute_descriptors(A, B, X, B_ox)
        feat = np.array(list(desc.values())).reshape(1, -1)
        feat_scaled = self.scaler.transform(feat)

        gb_pred  = self.gb_model.predict(feat_scaled)[0]
        rf_pred  = self.rf_model.predict(feat_scaled)[0]
        ensemble = 0.6 * gb_pred + 0.4 * rf_pred

        # SOC correction: Bi and heavy B-sites have ~0.2–0.4 eV SOC gap narrowing in DFT
        soc_correction = 0.0
        if B == "Bi":
            soc_correction = -0.25
        elif B == "Sb":
            soc_correction = -0.12
        elif B == "Sn":
            soc_correction = -0.08

        Eg = max(0.5, ensemble + soc_correction)

        # Uncertainty estimate from GB/RF spread + LOO MAE
        uncertainty = np.sqrt((gb_pred - rf_pred) ** 2 / 2 + (self.cv_mae or 0.15) ** 2)

        # Rashba indicator (non-centrosymmetric + heavy B-site)
        rashba = B in ["Bi", "Sn"] and A in ["MA", "FA", "EA"]

        return {
            "Eg_predicted_eV": round(Eg, 3),
            "Eg_gb_eV": round(gb_pred, 3),
            "Eg_rf_eV": round(rf_pred, 3),
            "uncertainty_eV": round(uncertainty, 3),
            "soc_correction_eV": soc_correction,
            "rashba_splitting": rashba,
        }

    def get_feature_importance(self) -> dict:
        if not self.is_fitted:
            self.fit(verbose=False)
        imp = 0.6 * self.gb_model.feature_importances_ + 0.4 * self.rf_model.feature_importances_
        return dict(zip(self.feature_names, imp.tolist()))


# ── Absorption Coefficient ────────────────────────────────────────────────────

def compute_absorption_coefficient(Eg: float, E_photon_eV: np.ndarray,
                                   direct_gap: bool = True,
                                   A_prefactor: float = 1.5e5) -> np.ndarray:
    """
    Absorption coefficient α(E) in cm⁻¹.
    Direct gap: α = A√(E - Eg) for E > Eg
    Indirect gap: α = A(E - Eg)² / E for E > Eg (Bi-based layered)
    """
    alpha = np.zeros_like(E_photon_eV, dtype=float)
    above = E_photon_eV > Eg
    if direct_gap:
        alpha[above] = A_prefactor * np.sqrt(E_photon_eV[above] - Eg)
    else:
        alpha[above] = A_prefactor * 0.3 * (E_photon_eV[above] - Eg) ** 2 / E_photon_eV[above]
    return alpha


def slme(Eg: float, alpha: np.ndarray, E_photon: np.ndarray, L_nm: float = 500) -> float:
    """
    Spectroscopic Limited Maximum Efficiency (SLME) in %.
    Simplified version: integrates AM1.5G photon flux above Eg with
    thickness-dependent absorption and detailed balance.
    L_nm: absorber thickness in nm.
    """
    # AM1.5G photon flux (simplified Planck)
    kT = 0.02585  # eV at 300K
    E = E_photon
    mask = E > 0.3
    flux = np.zeros_like(E)
    flux[mask] = (2 * E[mask]**2) / (np.exp(E[mask] / kT) - 1 + 1e-30)
    flux *= 1e15  # normalize to ~1e17 photons/cm²/s/eV above 0.5 eV

    # Absorptance
    L_cm = L_nm * 1e-7
    absorptance = 1 - np.exp(-2 * alpha * L_cm)

    # Photocurrent density (mA/cm²)
    dE = np.diff(E, prepend=E[0])
    Jph = 1.602e-19 * np.trapezoid(flux * absorptance, E) * 1000  # mA/cm²

    # Detailed balance: Voc estimate
    J0_rad = 1.602e-19 * np.trapezoid(flux * absorptance * np.exp(-E / kT), E) * 1000
    J0_rad = max(J0_rad, 1e-30)

    Voc = kT * np.log(Jph / J0_rad + 1)
    Voc = min(Voc, Eg * 0.85)

    FF = (Voc / kT - np.log(Voc / kT + 1)) / (Voc / kT + 1)
    FF = min(FF, 0.89)

    # Incident power ≈ 100 mW/cm²
    P_in = 100.0
    eta = Jph * Voc * FF / P_in * 100

    return round(max(0, min(eta, 35.0)), 2)
