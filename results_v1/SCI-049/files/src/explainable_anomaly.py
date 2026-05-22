"""
Module 5: Explainable Anomaly Detection
- SHAP-based feature attribution for anomalies
- Rule extraction from anomaly decisions
- Automated root cause analysis
"""
import numpy as np
from sklearn.ensemble import IsolationForest
from sklearn.tree import DecisionTreeClassifier
from sklearn.preprocessing import StandardScaler


class ExplainableAnomalyDetector:
    """Anomaly detector with built-in explainability."""

    def __init__(self, contamination=0.05, n_estimators=200, random_state=42):
        self.model = IsolationForest(
            contamination=contamination, n_estimators=n_estimators,
            random_state=random_state, n_jobs=-1
        )
        self.scaler = StandardScaler()
        self.feature_names = None

    def fit_predict_explain(self, X: np.ndarray, feature_names=None) -> dict:
        self.feature_names = feature_names or [f"feat_{i}" for i in range(X.shape[1])]
        X_scaled = self.scaler.fit_transform(X)

        labels = self.model.fit_predict(X_scaled)
        scores = self.model.decision_function(X_scaled)
        anomaly_mask = labels == -1

        # Permutation-based feature importance
        global_importance = self._permutation_importance(X_scaled, scores)

        # Per-anomaly local explanations
        local_explanations = self._local_explanations(X_scaled, anomaly_mask)

        # Rule extraction
        rules = self._extract_rules(X_scaled, labels)

        # Root cause clustering
        root_causes = self._root_cause_analysis(X, anomaly_mask)

        return {
            "labels": labels,
            "scores": scores,
            "n_anomalies": int(np.sum(anomaly_mask)),
            "global_feature_importance": global_importance,
            "local_explanations": local_explanations,
            "decision_rules": rules,
            "root_causes": root_causes,
        }

    def _permutation_importance(self, X, base_scores):
        importance = {}
        for i, name in enumerate(self.feature_names):
            X_perm = X.copy()
            np.random.RandomState(i).shuffle(X_perm[:, i])
            perm_scores = self.model.decision_function(X_perm)
            importance[name] = float(np.mean(np.abs(base_scores - perm_scores)))
        total = sum(importance.values()) or 1.0
        return {k: round(v / total, 4) for k, v in
                sorted(importance.items(), key=lambda x: -x[1])}

    def _local_explanations(self, X, anomaly_mask, max_samples=50):
        anomaly_indices = np.where(anomaly_mask)[0]
        if len(anomaly_indices) > max_samples:
            anomaly_indices = anomaly_indices[:max_samples]

        normal_mean = np.mean(X[~anomaly_mask], axis=0)
        normal_std = np.std(X[~anomaly_mask], axis=0) + 1e-10

        explanations = []
        for idx in anomaly_indices:
            z_scores = (X[idx] - normal_mean) / normal_std
            top_features = np.argsort(np.abs(z_scores))[::-1][:5]
            explanation = {
                "index": int(idx),
                "anomaly_score": float(self.model.decision_function(X[idx:idx+1])[0]),
                "top_contributing_features": [
                    {
                        "feature": self.feature_names[f],
                        "z_score": float(z_scores[f]),
                        "value": float(X[idx, f]),
                        "normal_mean": float(normal_mean[f]),
                        "direction": "high" if z_scores[f] > 0 else "low",
                    }
                    for f in top_features
                ],
            }
            explanations.append(explanation)
        return explanations

    def _extract_rules(self, X, labels, max_depth=4):
        """Train a surrogate decision tree to extract interpretable rules."""
        binary_labels = (labels == -1).astype(int)
        tree = DecisionTreeClassifier(max_depth=max_depth, random_state=42)
        tree.fit(X, binary_labels)

        rules = []
        self._traverse_tree(tree.tree_, 0, [], rules, tree)
        anomaly_rules = [r for r in rules if r["prediction"] == "anomaly"]
        return anomaly_rules[:10]

    def _traverse_tree(self, tree, node_id, path, rules, clf):
        if tree.children_left[node_id] == -1:  # leaf
            counts = tree.value[node_id][0]
            pred = "anomaly" if np.argmax(counts) == 1 else "normal"
            if pred == "anomaly" and counts[1] > 0:
                rules.append({
                    "conditions": list(path),
                    "prediction": pred,
                    "support": int(counts[1]),
                    "confidence": float(counts[1] / sum(counts)) if sum(counts) > 0 else 0,
                })
            return

        feat = self.feature_names[tree.feature[node_id]]
        thresh = round(float(tree.threshold[node_id]), 4)

        left_path = path + [f"{feat} <= {thresh}"]
        self._traverse_tree(tree, tree.children_left[node_id], left_path, rules, clf)

        right_path = path + [f"{feat} > {thresh}"]
        self._traverse_tree(tree, tree.children_right[node_id], right_path, rules, clf)

    def _root_cause_analysis(self, X, anomaly_mask):
        """Cluster anomalies to identify common root causes."""
        if np.sum(anomaly_mask) < 3:
            return [{"type": "insufficient_anomalies", "count": int(np.sum(anomaly_mask))}]

        anomalies = X[anomaly_mask]
        normals = X[~anomaly_mask]

        normal_mean = np.mean(normals, axis=0)
        normal_std = np.std(normals, axis=0) + 1e-10

        z_matrix = (anomalies - normal_mean) / normal_std
        dominant_features = np.argmax(np.abs(z_matrix), axis=1)

        from collections import Counter
        feature_counts = Counter(dominant_features)

        root_causes = []
        for feat_idx, count in feature_counts.most_common(5):
            feat_mask = dominant_features == feat_idx
            root_causes.append({
                "feature": self.feature_names[feat_idx],
                "count": count,
                "fraction": round(count / len(anomalies), 3),
                "mean_z_score": float(np.mean(z_matrix[feat_mask, feat_idx])),
                "direction": "high" if np.mean(z_matrix[feat_mask, feat_idx]) > 0 else "low",
            })
        return root_causes
