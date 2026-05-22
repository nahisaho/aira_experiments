"""
Module 2: Multivariate Outlier Detection
- Isolation Forest (ensemble-based)
- Deep SVDD (deep one-class classification)
"""
import numpy as np
from sklearn.ensemble import IsolationForest
from sklearn.preprocessing import StandardScaler


class IsolationForestDetector:
    """Enhanced Isolation Forest with feature importance tracking."""

    def __init__(self, contamination=0.05, n_estimators=200, random_state=42):
        self.model = IsolationForest(
            contamination=contamination,
            n_estimators=n_estimators,
            random_state=random_state,
            n_jobs=-1,
        )
        self.scaler = StandardScaler()
        self.contamination = contamination

    def fit_predict(self, X: np.ndarray, feature_names=None) -> dict:
        X_scaled = self.scaler.fit_transform(X)
        labels = self.model.fit_predict(X_scaled)
        scores = self.model.decision_function(X_scaled)

        anomaly_mask = labels == -1
        n_anomalies = int(np.sum(anomaly_mask))

        feature_importance = self._compute_feature_importance(X_scaled, feature_names)

        return {
            "labels": labels,
            "scores": scores,
            "n_anomalies": n_anomalies,
            "anomaly_rate": n_anomalies / len(X),
            "anomaly_indices": np.where(anomaly_mask)[0].tolist(),
            "feature_importance": feature_importance,
            "score_stats": {
                "mean": float(np.mean(scores)),
                "std": float(np.std(scores)),
                "min": float(np.min(scores)),
                "max": float(np.max(scores)),
                "threshold": float(np.percentile(scores, self.contamination * 100)),
            },
        }

    def _compute_feature_importance(self, X, feature_names):
        importances = {}
        base_scores = self.model.decision_function(X)
        names = feature_names or [f"feature_{i}" for i in range(X.shape[1])]
        for i, name in enumerate(names):
            X_perm = X.copy()
            np.random.shuffle(X_perm[:, i])
            perm_scores = self.model.decision_function(X_perm)
            importances[name] = float(np.mean(np.abs(base_scores - perm_scores)))
        total = sum(importances.values()) or 1.0
        return {k: v / total for k, v in importances.items()}


class DeepSVDDDetector:
    """Deep Support Vector Data Description (numpy-only lightweight version).
    Uses a simple autoencoder-based approach as Deep SVDD proxy without PyTorch dependency.
    """

    def __init__(self, encoding_dim=8, nu=0.05, random_state=42):
        self.encoding_dim = encoding_dim
        self.nu = nu
        self.rng = np.random.RandomState(random_state)
        self.center_ = None
        self.radius_ = None
        self.W1_ = None
        self.W2_ = None

    def fit(self, X: np.ndarray, n_epochs=100, lr=0.01):
        n_features = X.shape[1]
        self.W1_ = self.rng.randn(n_features, self.encoding_dim) * 0.1
        self.W2_ = self.rng.randn(self.encoding_dim, n_features) * 0.1

        scaler = StandardScaler()
        X_s = scaler.fit_transform(X)
        self._scaler = scaler

        for epoch in range(n_epochs):
            h = np.maximum(0, X_s @ self.W1_)  # ReLU
            recon = h @ self.W2_
            error = recon - X_s
            # Backprop (simplified gradient descent)
            grad_W2 = h.T @ error / len(X_s)
            grad_h = error @ self.W2_.T
            grad_h[h <= 0] = 0
            grad_W1 = X_s.T @ grad_h / len(X_s)
            self.W1_ -= lr * grad_W1
            self.W2_ -= lr * grad_W2

        embeddings = np.maximum(0, X_s @ self.W1_)
        self.center_ = np.mean(embeddings, axis=0)
        dists = np.linalg.norm(embeddings - self.center_, axis=1)
        self.radius_ = np.percentile(dists, (1 - self.nu) * 100)
        return self

    def predict(self, X: np.ndarray) -> dict:
        X_s = self._scaler.transform(X)
        embeddings = np.maximum(0, X_s @ self.W1_)
        dists = np.linalg.norm(embeddings - self.center_, axis=1)
        labels = np.where(dists > self.radius_, -1, 1)
        scores = -(dists - self.radius_)  # negative = anomaly

        return {
            "labels": labels,
            "distances": dists,
            "scores": scores,
            "n_anomalies": int(np.sum(labels == -1)),
            "radius": float(self.radius_),
            "center_norm": float(np.linalg.norm(self.center_)),
        }
