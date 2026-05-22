"""
Module 4: Machine learning prediction of adsorption isotherms.

Uses molecular/structural descriptors to predict CO2/H2 uptake
via gradient boosting, random forest, and neural network models.
"""
import json
import logging
import pickle
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


class IsothermModel:
    """Base class for isotherm fitting models (Langmuir, Dual-site Langmuir, BET)."""

    @staticmethod
    def langmuir(P: np.ndarray, q_sat: float, K: float) -> np.ndarray:
        """Single-site Langmuir: q = q_sat * K * P / (1 + K * P)"""
        return q_sat * K * P / (1.0 + K * P)

    @staticmethod
    def dual_site_langmuir(P: np.ndarray, q1: float, K1: float,
                            q2: float, K2: float) -> np.ndarray:
        """Dual-site Langmuir model."""
        return (q1 * K1 * P / (1.0 + K1 * P) +
                q2 * K2 * P / (1.0 + K2 * P))

    @staticmethod
    def freundlich(P: np.ndarray, K: float, n: float) -> np.ndarray:
        """Freundlich model: q = K * P^(1/n)"""
        return K * np.power(P, 1.0 / n)

    @staticmethod
    def fit_langmuir(P: np.ndarray, q: np.ndarray) -> Tuple[float, float]:
        """Fit Langmuir parameters via linearization: P/q = 1/(q_sat*K) + P/q_sat."""
        mask = (P > 0) & (q > 0)
        P_f, q_f = P[mask], q[mask]
        if len(P_f) < 2:
            return 0.0, 0.0
        y = P_f / q_f
        coeffs = np.polyfit(P_f, y, 1)
        q_sat = 1.0 / coeffs[0] if coeffs[0] != 0 else 0.0
        K = coeffs[0] / coeffs[1] if coeffs[1] != 0 else 0.0
        return q_sat, K


class AdsorptionPredictor:
    """
    ML-based predictor for gas adsorption in MOFs.

    Supports multiple targets:
    - CO2 uptake at 0.15 bar (flue gas)
    - CO2 uptake at 0.0004 bar (DAC)
    - H2 uptake at 100 bar
    - CO2/N2 selectivity
    - Heat of adsorption
    """

    def __init__(self, model_type: str = "gradient_boosting",
                 n_estimators: int = 500, max_depth: int = 8,
                 learning_rate: float = 0.05, random_state: int = 42):
        self.model_type = model_type
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.learning_rate = learning_rate
        self.random_state = random_state
        self.models: Dict[str, Any] = {}
        self.scalers: Dict[str, Any] = {}
        self.feature_names: List[str] = []
        self.metrics: Dict[str, Dict] = {}

    def _create_model(self):
        """Create ML model based on configuration."""
        if self.model_type == "gradient_boosting":
            from sklearn.ensemble import GradientBoostingRegressor
            return GradientBoostingRegressor(
                n_estimators=self.n_estimators,
                max_depth=self.max_depth,
                learning_rate=self.learning_rate,
                random_state=self.random_state,
                subsample=0.8,
                min_samples_split=5,
                min_samples_leaf=3,
            )
        elif self.model_type == "random_forest":
            from sklearn.ensemble import RandomForestRegressor
            return RandomForestRegressor(
                n_estimators=self.n_estimators,
                max_depth=self.max_depth,
                random_state=self.random_state,
                n_jobs=-1,
            )
        elif self.model_type == "xgboost":
            from xgboost import XGBRegressor
            return XGBRegressor(
                n_estimators=self.n_estimators,
                max_depth=self.max_depth,
                learning_rate=self.learning_rate,
                random_state=self.random_state,
                tree_method="hist",
            )
        else:
            raise ValueError(f"Unknown model type: {self.model_type}")

    def train(self, X: np.ndarray, y: np.ndarray,
              target_name: str, feature_names: List[str],
              test_size: float = 0.2, cv_folds: int = 5) -> Dict:
        """Train model with cross-validation and evaluation."""
        from sklearn.model_selection import train_test_split, cross_val_score
        from sklearn.preprocessing import StandardScaler
        from sklearn.metrics import mean_absolute_error, r2_score

        self.feature_names = feature_names

        # Remove NaN/Inf
        mask = np.all(np.isfinite(X), axis=1) & np.isfinite(y)
        X_clean, y_clean = X[mask], y[mask]
        logger.info(f"Training {target_name}: {len(y_clean)} samples "
                     f"({len(y) - len(y_clean)} removed)")

        # Scale features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_clean)
        self.scalers[target_name] = scaler

        # Train/test split
        X_train, X_test, y_train, y_test = train_test_split(
            X_scaled, y_clean, test_size=test_size,
            random_state=self.random_state
        )

        # Train model
        model = self._create_model()
        model.fit(X_train, y_train)
        self.models[target_name] = model

        # Evaluate
        y_pred_train = model.predict(X_train)
        y_pred_test = model.predict(X_test)

        # Cross-validation
        cv_scores = cross_val_score(
            self._create_model(), X_scaled, y_clean,
            cv=cv_folds, scoring="r2"
        )

        metrics = {
            "target": target_name,
            "n_samples": int(len(y_clean)),
            "n_features": int(X_clean.shape[1]),
            "train_r2": round(float(r2_score(y_train, y_pred_train)), 4),
            "test_r2": round(float(r2_score(y_test, y_pred_test)), 4),
            "train_mae": round(float(mean_absolute_error(y_train, y_pred_train)), 4),
            "test_mae": round(float(mean_absolute_error(y_test, y_pred_test)), 4),
            "cv_r2_mean": round(float(cv_scores.mean()), 4),
            "cv_r2_std": round(float(cv_scores.std()), 4),
        }

        # Feature importance
        if hasattr(model, "feature_importances_"):
            importances = model.feature_importances_
            feat_imp = sorted(
                zip(feature_names, importances),
                key=lambda x: x[1], reverse=True
            )
            metrics["top_features"] = [
                {"name": n, "importance": round(float(v), 4)}
                for n, v in feat_imp[:10]
            ]

        self.metrics[target_name] = metrics
        logger.info(f"{target_name}: R²_test={metrics['test_r2']}, "
                     f"MAE_test={metrics['test_mae']}")
        return metrics

    def predict(self, X: np.ndarray, target_name: str) -> np.ndarray:
        """Predict using trained model."""
        if target_name not in self.models:
            raise ValueError(f"No trained model for {target_name}")
        X_scaled = self.scalers[target_name].transform(X)
        return self.models[target_name].predict(X_scaled)

    def predict_isotherm(self, X_single: np.ndarray,
                          pressure_points: List[float]) -> np.ndarray:
        """Predict full isotherm by predicting at multiple pressures."""
        # Augment features with log(P) for pressure-dependent prediction
        predictions = []
        for p in pressure_points:
            X_aug = np.append(X_single, [np.log10(max(p, 1e-10))])
            target = f"loading_P{p:.4f}"
            if target in self.models:
                pred = self.predict(X_aug.reshape(1, -1), target)
                predictions.append(pred[0])
            else:
                predictions.append(0.0)
        return np.array(predictions)

    def save_models(self, output_dir: Path):
        """Save trained models to disk."""
        output_dir.mkdir(parents=True, exist_ok=True)
        for target, model in self.models.items():
            model_path = output_dir / f"model_{target}.pkl"
            with open(model_path, "wb") as f:
                pickle.dump(model, f)
            scaler_path = output_dir / f"scaler_{target}.pkl"
            with open(scaler_path, "wb") as f:
                pickle.dump(self.scalers[target], f)

        meta = {
            "model_type": self.model_type,
            "feature_names": self.feature_names,
            "targets": list(self.models.keys()),
            "metrics": self.metrics,
        }
        with open(output_dir / "model_metadata.json", "w") as f:
            json.dump(meta, f, indent=2)

    def load_models(self, model_dir: Path):
        """Load pre-trained models."""
        with open(model_dir / "model_metadata.json") as f:
            meta = json.load(f)
        self.feature_names = meta["feature_names"]
        self.metrics = meta.get("metrics", {})
        for target in meta["targets"]:
            with open(model_dir / f"model_{target}.pkl", "rb") as f:
                self.models[target] = pickle.load(f)
            with open(model_dir / f"scaler_{target}.pkl", "rb") as f:
                self.scalers[target] = pickle.load(f)


class MultiTargetPredictor:
    """
    Multi-target adsorption predictor for screening.

    Predicts multiple properties simultaneously:
    1. CO2 uptake at DAC conditions (0.0004 bar, 298 K)
    2. CO2 uptake at flue gas conditions (0.15 bar, 298 K)
    3. CO2 uptake at 1 bar, 298 K
    4. H2 uptake at 100 bar, 77 K
    5. CO2/N2 selectivity
    6. Isosteric heat of adsorption (Qst)
    """

    TARGETS = [
        "CO2_0.0004bar_298K",
        "CO2_0.15bar_298K",
        "CO2_1bar_298K",
        "H2_100bar_77K",
        "CO2_N2_selectivity",
        "Qst_CO2_kJ_mol",
    ]

    def __init__(self, **kwargs):
        self.predictors = {
            target: AdsorptionPredictor(**kwargs)
            for target in self.TARGETS
        }

    def train_all(self, X: np.ndarray, targets: Dict[str, np.ndarray],
                   feature_names: List[str], **kwargs) -> Dict:
        """Train all target models."""
        all_metrics = {}
        for target_name in self.TARGETS:
            if target_name in targets:
                metrics = self.predictors[target_name].train(
                    X, targets[target_name], target_name,
                    feature_names, **kwargs
                )
                all_metrics[target_name] = metrics
        return all_metrics

    def predict_all(self, X: np.ndarray) -> Dict[str, np.ndarray]:
        """Predict all targets."""
        results = {}
        for target_name, predictor in self.predictors.items():
            if target_name in predictor.models:
                results[target_name] = predictor.predict(X, target_name)
        return results

    def save_all(self, output_dir: Path):
        for target, pred in self.predictors.items():
            if pred.models:
                pred.save_models(output_dir / target)
