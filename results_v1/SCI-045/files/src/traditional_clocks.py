"""
Traditional Epigenetic Clocks: Horvath-like and GrimAge-like implementations.
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNet, ElasticNetCV
from sklearn.model_selection import KFold
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
from scipy import stats


class HorvathStyleClock:
    def __init__(self, alpha=0.1, l1_ratio=0.5, max_iter=10000):
        self.scaler = StandardScaler()
        self.model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, max_iter=max_iter, random_state=42)

    def fit(self, X, age):
        X_scaled = self.scaler.fit_transform(X)
        self.model.fit(X_scaled, np.log(age + 1))
        self.n_nonzero_ = np.sum(self.model.coef_ != 0)
        return self

    def predict(self, X):
        return np.exp(self.model.predict(self.scaler.transform(X))) - 1

    def get_clock_cpgs(self, cpg_names):
        return list(np.array(cpg_names)[self.model.coef_ != 0])


class GrimAgeStyleClock:
    def __init__(self, n_surrogate_models=7, alpha=0.05, l1_ratio=0.7):
        self.n_surrogates = n_surrogate_models
        self.surrogate_models = []
        self.meta_model = ElasticNet(alpha=alpha, l1_ratio=l1_ratio, max_iter=10000, random_state=42)
        self.scaler = StandardScaler()
        self.meta_scaler = StandardScaler()

    def fit(self, X, age):
        n = X.shape[0]
        X_scaled = self.scaler.fit_transform(X)
        surrogate_preds = np.zeros((n, self.n_surrogates))
        for i in range(self.n_surrogates):
            rng = np.random.RandomState(i + 100)
            target = age * (0.5 + 0.1 * i) + rng.normal(0, 5, n)
            model = ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=10000, random_state=42 + i)
            model.fit(X_scaled, target)
            surrogate_preds[:, i] = model.predict(X_scaled)
            self.surrogate_models.append(model)
        meta_features = np.column_stack([surrogate_preds, age])
        meta_scaled = self.meta_scaler.fit_transform(meta_features)
        self.meta_model.fit(meta_scaled, age)
        return self

    def predict(self, X, age_hint=None):
        X_scaled = self.scaler.transform(X)
        surrogate_preds = np.column_stack([m.predict(X_scaled) for m in self.surrogate_models])
        if age_hint is None:
            age_hint = np.mean(surrogate_preds, axis=1)
        meta_features = np.column_stack([surrogate_preds, age_hint])
        return self.meta_model.predict(self.meta_scaler.transform(meta_features))


class ImprovedElasticNetClock:
    def __init__(self):
        self.scaler = StandardScaler()
        self.model = None

    def _engineer_features(self, X):
        row_mean = X.mean(axis=1).reshape(-1, 1)
        row_std = X.std(axis=1).reshape(-1, 1)
        top_var_idx = np.argsort(X.var(axis=0))[-50:]
        return np.hstack([X, row_mean, row_std, X[:, top_var_idx] ** 2])

    def fit(self, X, age):
        X_eng = self._engineer_features(X)
        X_scaled = self.scaler.fit_transform(X_eng)
        self.model = ElasticNetCV(
            l1_ratio=[0.1, 0.3, 0.5, 0.7, 0.9],
            alphas=np.logspace(-4, 0, 20),
            cv=5, max_iter=20000, random_state=42, n_jobs=-1
        )
        self.model.fit(X_scaled, np.log(age + 1))
        self.n_nonzero_ = np.sum(self.model.coef_ != 0)
        return self

    def predict(self, X):
        X_eng = self._engineer_features(X)
        return np.exp(self.model.predict(self.scaler.transform(X_eng))) - 1


def evaluate_clock(y_true, y_pred, label=""):
    mae = mean_absolute_error(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    r2 = r2_score(y_true, y_pred)
    r, p = stats.pearsonr(y_true, y_pred)
    return {
        "model": label, "MAE": round(mae, 3), "RMSE": round(rmse, 3),
        "R2": round(r2, 4), "Pearson_r": round(r, 4),
        "Pearson_p": float(f"{p:.2e}"), "Median_AE": round(float(np.median(np.abs(y_true - y_pred))), 3),
        "n_samples": len(y_true),
    }


def cross_validate_clock(clock_class, X, y, n_folds=5, **kwargs):
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)
    all_preds = np.zeros(len(y))
    fold_metrics = []
    for fold, (tr, te) in enumerate(kf.split(X)):
        clock = clock_class(**kwargs)
        clock.fit(X[tr], y[tr])
        all_preds[te] = clock.predict(X[te])
        fold_metrics.append(evaluate_clock(y[te], all_preds[te], f"fold_{fold}"))
    return all_preds, evaluate_clock(y, all_preds, "CV_overall"), fold_metrics
