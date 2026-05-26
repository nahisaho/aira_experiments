"""
Baseline models for comparison.
1. Persistence forecast (no change)
2. Climatology forecast (mean state)
3. Linear regression model
"""

import torch
import torch.nn as nn
import numpy as np
from sklearn.linear_model import Ridge


class PersistenceModel:
    """Predict current state = previous state (no change)."""

    def predict(self, input_pressure, input_surface):
        return input_pressure.clone(), input_surface.clone()


class ClimatologyModel:
    """Predict climatological mean state."""

    def __init__(self):
        self.mean_pressure = None
        self.mean_surface = None

    def fit(self, dataset):
        all_p = torch.stack([d['target_pressure'] for d in dataset])
        all_s = torch.stack([d['target_surface'] for d in dataset])
        self.mean_pressure = all_p.mean(dim=0)
        self.mean_surface = all_s.mean(dim=0)

    def predict(self, input_pressure, input_surface):
        B = input_pressure.shape[0]
        return (self.mean_pressure.unsqueeze(0).expand(B, -1, -1, -1),
                self.mean_surface.unsqueeze(0).expand(B, -1, -1))


class LinearRegressionModel:
    """Linear regression baseline."""

    def __init__(self):
        self.model_p = None
        self.model_s = None

    def fit(self, dataset):
        X_p = torch.stack([d['input_pressure'] for d in dataset]).numpy()
        Y_p = torch.stack([d['target_pressure'] for d in dataset]).numpy()
        X_s = torch.stack([d['input_surface'] for d in dataset]).numpy()
        Y_s = torch.stack([d['target_surface'] for d in dataset]).numpy()

        n = X_p.shape[0]
        X_p_flat = X_p.reshape(n, -1)
        Y_p_flat = Y_p.reshape(n, -1)
        X_s_flat = X_s.reshape(n, -1)
        Y_s_flat = Y_s.reshape(n, -1)

        self.model_p = Ridge(alpha=1.0)
        self.model_p.fit(X_p_flat, Y_p_flat)
        self.model_s = Ridge(alpha=1.0)
        self.model_s.fit(X_s_flat, Y_s_flat)

        self.p_shape = X_p.shape[1:]
        self.s_shape = X_s.shape[1:]

    def predict(self, input_pressure, input_surface):
        B = input_pressure.shape[0]
        X_p = input_pressure.numpy().reshape(B, -1)
        X_s = input_surface.numpy().reshape(B, -1)

        Y_p = self.model_p.predict(X_p).reshape(B, *self.p_shape)
        Y_s = self.model_s.predict(X_s).reshape(B, *self.s_shape)

        return torch.from_numpy(Y_p.astype(np.float32)), torch.from_numpy(Y_s.astype(np.float32))
