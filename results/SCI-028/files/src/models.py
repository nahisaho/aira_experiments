"""
Disruption Prediction Models for Tokamak Plasmas
Implements LSTM, CNN-LSTM, and Physics-Informed Neural Network architectures.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np


class LSTMPredictor(nn.Module):
    """LSTM-based disruption predictor (baseline)."""
    def __init__(self, input_dim=11, hidden_dim=64, num_layers=2, dropout=0.3):
        super().__init__()
        self.lstm = nn.LSTM(input_dim, hidden_dim, num_layers,
                           batch_first=True, dropout=dropout, bidirectional=False)
        self.fc = nn.Sequential(
            nn.Linear(hidden_dim, 32),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(32, 1)
        )

    def forward(self, x):
        out, _ = self.lstm(x)
        out = self.fc(out[:, -1, :])
        return out.squeeze(-1)


class CNNLSTMPredictor(nn.Module):
    """CNN-LSTM hybrid for multi-scale temporal feature extraction."""
    def __init__(self, input_dim=11, cnn_channels=32, lstm_hidden=64, num_layers=2, dropout=0.3):
        super().__init__()
        # Multi-scale CNN
        self.conv1 = nn.Conv1d(input_dim, cnn_channels, kernel_size=3, padding=1)
        self.conv2 = nn.Conv1d(input_dim, cnn_channels, kernel_size=7, padding=3)
        self.conv3 = nn.Conv1d(input_dim, cnn_channels, kernel_size=15, padding=7)
        self.bn = nn.BatchNorm1d(cnn_channels * 3)
        self.lstm = nn.LSTM(cnn_channels * 3, lstm_hidden, num_layers,
                           batch_first=True, dropout=dropout)
        self.fc = nn.Sequential(
            nn.Linear(lstm_hidden, 32),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(32, 1)
        )

    def forward(self, x):
        # x: (batch, seq_len, features) -> CNN expects (batch, features, seq_len)
        x_t = x.transpose(1, 2)
        c1 = F.relu(self.conv1(x_t))
        c2 = F.relu(self.conv2(x_t))
        c3 = F.relu(self.conv3(x_t))
        c = torch.cat([c1, c2, c3], dim=1)
        c = self.bn(c)
        c = c.transpose(1, 2)  # back to (batch, seq_len, channels)
        out, _ = self.lstm(c)
        out = self.fc(out[:, -1, :])
        return out.squeeze(-1)


class PhysicsInformedPredictor(nn.Module):
    """
    Physics-Informed ML model for disruption prediction.
    Incorporates MHD stability constraints as auxiliary loss terms.
    """
    def __init__(self, input_dim=11, hidden_dim=64, num_layers=2, dropout=0.3):
        super().__init__()
        # Physics feature extractor
        self.physics_net = nn.Sequential(
            nn.Linear(input_dim, 32),
            nn.ReLU(),
            nn.Linear(32, 4)  # stability indicators: Troyon, Greenwald, q95, locked mode
        )
        # Main predictor with physics features
        self.lstm = nn.LSTM(input_dim + 4, hidden_dim, num_layers,
                           batch_first=True, dropout=dropout)
        self.fc = nn.Sequential(
            nn.Linear(hidden_dim, 32),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(32, 1)
        )

    def compute_physics_features(self, x):
        """Compute physics-based stability indicators per timestep."""
        batch, seq, feat = x.shape
        x_flat = x.reshape(-1, feat)
        physics = self.physics_net(x_flat)
        return physics.reshape(batch, seq, -1)

    def physics_loss(self, x, physics_features):
        """
        Physics-informed regularization based on MHD stability criteria.
        Penalizes predictions inconsistent with known stability boundaries.
        """
        # Extract raw features
        beta_N = x[:, :, 3]  # normalized beta
        q95 = x[:, :, 6]    # safety factor
        n_G = x[:, :, 7]    # Greenwald fraction
        Ip = x[:, :, 0]     # plasma current

        # Troyon limit: beta_N < C_T * Ip / (a * B) ~ approximate
        troyon_indicator = physics_features[:, :, 0]
        troyon_target = torch.clamp(beta_N / 3.5, 0, 1)  # normalized distance to Troyon limit
        loss_troyon = F.mse_loss(troyon_indicator, troyon_target)

        # Greenwald limit: n_G should stay below 1
        gw_indicator = physics_features[:, :, 1]
        gw_target = torch.clamp(n_G, 0, 2)
        loss_gw = F.mse_loss(gw_indicator, gw_target)

        # q95 stability: q95 > 2 for stability
        q_indicator = physics_features[:, :, 2]
        q_target = torch.clamp(1.0 / (q95 + 1e-6), 0, 1)
        loss_q = F.mse_loss(q_indicator, q_target)

        return 0.1 * (loss_troyon + loss_gw + loss_q)

    def forward(self, x):
        physics_feats = self.compute_physics_features(x)
        combined = torch.cat([x, physics_feats], dim=-1)
        out, _ = self.lstm(combined)
        out = self.fc(out[:, -1, :])
        return out.squeeze(-1), physics_feats


class TearingModeDetector(nn.Module):
    """Multi-task model for TM/NTM detection alongside disruption prediction."""
    def __init__(self, input_dim=11, hidden_dim=48, dropout=0.3):
        super().__init__()
        self.shared_lstm = nn.LSTM(input_dim, hidden_dim, 2,
                                  batch_first=True, dropout=dropout)
        # Disruption head
        self.disrupt_head = nn.Sequential(
            nn.Linear(hidden_dim, 16), nn.ReLU(), nn.Linear(16, 1)
        )
        # TM/NTM detection head
        self.tm_head = nn.Sequential(
            nn.Linear(hidden_dim, 16), nn.ReLU(), nn.Linear(16, 1)
        )

    def forward(self, x):
        out, _ = self.shared_lstm(x)
        last = out[:, -1, :]
        disrupt_out = self.disrupt_head(last).squeeze(-1)
        tm_out = self.tm_head(last).squeeze(-1)
        return disrupt_out, tm_out


class TransferablePredictor(nn.Module):
    """
    Model with domain-adaptive architecture for cross-device transfer learning.
    Uses a shared feature extractor + device-specific adaptation layers.
    """
    def __init__(self, input_dim=11, shared_hidden=64, adapt_hidden=32, dropout=0.3):
        super().__init__()
        # Shared feature extractor (frozen after pretraining)
        self.shared_encoder = nn.LSTM(input_dim, shared_hidden, 2,
                                     batch_first=True, dropout=dropout)
        # Domain adaptation layer
        self.adapt = nn.Sequential(
            nn.Linear(shared_hidden, adapt_hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        # Prediction head
        self.head = nn.Linear(adapt_hidden, 1)

    def forward(self, x):
        out, _ = self.shared_encoder(x)
        adapted = self.adapt(out[:, -1, :])
        pred = self.head(adapted).squeeze(-1)
        return pred

    def freeze_shared(self):
        for param in self.shared_encoder.parameters():
            param.requires_grad = False

    def unfreeze_shared(self):
        for param in self.shared_encoder.parameters():
            param.requires_grad = True
