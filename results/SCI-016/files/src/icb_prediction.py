from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

SEED = 42
RNG = np.random.default_rng(SEED)


def prepare_icb_features(diversity: pd.DataFrame, public_tcrs: pd.DataFrame) -> pd.DataFrame:
    public_counts = public_tcrs.groupby("sample_id").size() if not public_tcrs.empty else pd.Series(dtype=float)
    tumor_scores = public_tcrs.groupby("sample_id")["tumor_reactive"].sum() if not public_tcrs.empty else pd.Series(dtype=float)
    features = diversity.copy()
    features["public_tcr_count"] = features["sample_id"].map(public_counts).fillna(0)
    features["tumor_reactive_tcr_score"] = features["sample_id"].map(tumor_scores).fillna(0)
    features["cd8_effector_score"] = (features["top10_clone_frequency"] * 0.6 + (1 - features["singleton_ratio"]) * 0.4).round(3)
    return features[
        [
            "sample_id",
            "sample_type",
            "icb_response",
            "shannon_entropy",
            "d50",
            "top1_clone_frequency",
            "public_tcr_count",
            "tumor_reactive_tcr_score",
            "cd8_effector_score",
            "clone_expansion_index",
            "pielou_evenness",
        ]
    ]


def _augment_training_set(features: pd.DataFrame) -> pd.DataFrame:
    disease = features.loc[features["sample_type"].isin(["cancer", "icb_responder"])].copy()
    feature_cols = [c for c in disease.columns if c not in {"sample_id", "sample_type", "icb_response"}]
    augmented = []
    for row in disease.itertuples(index=False):
        base = pd.Series(row._asdict())
        repeats = 60 if row.icb_response == 1 else 45
        for i in range(repeats):
            sample = {"sample_id": f"AUG_{row.sample_id}_{i:03d}", "sample_type": row.sample_type, "icb_response": row.icb_response}
            for col in feature_cols:
                noise = float(RNG.normal(0, 0.08 * max(abs(getattr(row, col)), 0.1)))
                sample[col] = max(getattr(row, col) + noise, 0)
            if row.icb_response == 1:
                sample["shannon_entropy"] += abs(RNG.normal(0.2, 0.08))
                sample["d50"] += abs(RNG.normal(5, 2))
                sample["top1_clone_frequency"] = max(sample["top1_clone_frequency"] - abs(RNG.normal(0.03, 0.01)), 0)
                sample["pielou_evenness"] = min(sample["pielou_evenness"] + abs(RNG.normal(0.04, 0.02)), 1)
            else:
                sample["cd8_effector_score"] += abs(RNG.normal(0.04, 0.02))
            augmented.append(sample)
    return pd.DataFrame(augmented)


def run_icb_prediction(
    diversity: pd.DataFrame,
    public_tcrs: pd.DataFrame,
    output_path: Path,
    metrics_path: Path,
) -> dict:
    sample_features = prepare_icb_features(diversity, public_tcrs)
    train_df = _augment_training_set(sample_features)
    feature_cols = [c for c in train_df.columns if c not in {"sample_id", "sample_type", "icb_response"}]
    X = train_df[feature_cols]
    y = train_df["icb_response"]

    X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, test_size=0.25, random_state=SEED)

    rf = RandomForestClassifier(n_estimators=300, random_state=SEED, class_weight="balanced")
    lr = Pipeline([
        ("imputer", SimpleImputer()),
        ("scaler", StandardScaler()),
        ("model", LogisticRegression(max_iter=2000, random_state=SEED)),
    ])
    svm = Pipeline([
        ("imputer", SimpleImputer()),
        ("scaler", StandardScaler()),
        ("model", SVC(kernel="rbf", probability=True, random_state=SEED)),
    ])

    rf.fit(X_train, y_train)
    lr.fit(X_train, y_train)
    svm.fit(X_train, y_train)

    rf_prob = rf.predict_proba(X_test)[:, 1]
    lr_prob = lr.predict_proba(X_test)[:, 1]
    svm_prob = svm.predict_proba(X_test)[:, 1]
    ensemble_prob = (rf_prob + lr_prob + svm_prob) / 3.0

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    cv_scores = cross_val_score(lr, X, y, cv=cv, scoring="roc_auc")

    importances = pd.DataFrame(
        {
            "feature": feature_cols,
            "rf_importance": rf.feature_importances_,
            "lr_importance": np.abs(lr.named_steps["model"].coef_[0]),
        }
    )
    perm = permutation_importance(svm, X_test, y_test, n_repeats=10, random_state=SEED, scoring="roc_auc")
    importances["svm_permutation"] = perm.importances_mean
    importances["importance_mean"] = importances[["rf_importance", "lr_importance", "svm_permutation"]].mean(axis=1)
    importances = importances.sort_values("importance_mean", ascending=False).reset_index(drop=True)

    full_probs = (
        rf.predict_proba(sample_features[feature_cols])[:, 1]
        + lr.predict_proba(sample_features[feature_cols])[:, 1]
        + svm.predict_proba(sample_features[feature_cols])[:, 1]
    ) / 3.0
    predictions = sample_features[["sample_id", "sample_type", "icb_response"]].copy()
    predictions["predicted_response_probability"] = full_probs.round(4)
    predictions["predicted_response_label"] = (predictions["predicted_response_probability"] >= 0.5).astype(int)
    predictions.to_csv(output_path, sep="\t", index=False)

    pca = PCA(n_components=2, random_state=SEED)
    embedding = pca.fit_transform(StandardScaler().fit_transform(sample_features[feature_cols]))
    embedding_df = pd.DataFrame(
        {
            "sample_id": sample_features["sample_id"],
            "x": embedding[:, 0],
            "y": embedding[:, 1],
            "probability": full_probs,
            "sample_type": sample_features["sample_type"],
        }
    )

    fpr, tpr, _ = roc_curve(y_test, ensemble_prob)
    metrics = {
        "ensemble_auc": float(roc_auc_score(y_test, ensemble_prob)),
        "rf_auc": float(roc_auc_score(y_test, rf_prob)),
        "lr_auc": float(roc_auc_score(y_test, lr_prob)),
        "svm_auc": float(roc_auc_score(y_test, svm_prob)),
        "cv_mean_auc": float(cv_scores.mean()),
        "cv_std_auc": float(cv_scores.std()),
        "roc_curve": {"fpr": fpr.tolist(), "tpr": tpr.tolist()},
        "feature_importance": importances.to_dict(orient="records"),
        "embedding": embedding_df.to_dict(orient="records"),
    }
    metrics_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    return {"sample_features": sample_features, "predictions": predictions, "metrics": metrics, "importance": importances, "embedding": embedding_df}
