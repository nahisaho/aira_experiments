"""
model.py — GraphCast-style weather prediction model.

Full encode-process-decode architecture:
  1. Encoder: atmospheric variables → latent node features
  2. Processor: multi-layer GNN message passing on multi-scale mesh
  3. Decoder: latent features → predicted atmospheric residuals
"""

import torch
import torch.nn as nn
from typing import Dict, Optional, Tuple

from .encoder import AtmosphericEncoder, AtmosphericDecoder, N_INPUT_FEATURES
from .gnn_processor import GNNProcessor
from .physics import PhysicsLoss


class WeatherGNN(nn.Module):
    """
    GraphCast-style weather prediction model.

    Architecture:
      - Encoder: MLP with positional encoding → d_model latent space
      - Processor: N layers of GNN message passing
      - Decoder: MLP from latent space → atmospheric variable residuals

    Training:
      - Predicts residuals: X(t+dt) = X(t) + Decoder(Processor(Encoder(X(t))))
      - Loss: MSE on target + physics constraint losses
      - Autoregressive rollout for multi-step forecasts
    """

    def __init__(
        self,
        d_model: int = 256,
        n_encoder_layers: int = 2,
        n_processor_layers: int = 8,
        n_decoder_layers: int = 2,
        n_heads: int = 4,
        dropout: float = 0.1,
        n_input_features: int = N_INPUT_FEATURES,
        physics_weights: Optional[Dict[str, float]] = None,
    ):
        super().__init__()
        self.d_model = d_model

        # Encode-Process-Decode
        self.encoder = AtmosphericEncoder(
            n_input_features=n_input_features,
            d_model=d_model,
            n_layers=n_encoder_layers,
            dropout=dropout,
        )
        self.processor = GNNProcessor(
            d_model=d_model,
            n_layers=n_processor_layers,
            n_heads=n_heads,
            dropout=dropout,
        )
        self.decoder = AtmosphericDecoder(
            d_model=d_model,
            n_output_features=n_input_features,
            n_layers=n_decoder_layers,
            dropout=dropout,
        )

        # Physics loss
        pw = physics_weights or {}
        self.physics_loss = PhysicsLoss(**pw)

        # Learnable per-variable loss weights
        self.log_var_weights = nn.Parameter(torch.zeros(n_input_features))

    def forward(
        self,
        x: torch.Tensor,
        edge_index: torch.Tensor,
        lat: torch.Tensor,
        lon: torch.Tensor,
        time_features: Optional[torch.Tensor] = None,
        edge_attr: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Single-step forward pass.

        Args:
            x: (N, F) input atmospheric state
            edge_index: (2, E) graph edges
            lat, lon: (N,) coordinates in degrees
            time_features: (N, 2) optional temporal features
            edge_attr: (E, 3) optional edge features

        Returns:
            (N, F) predicted next atmospheric state (input + residual)
        """
        h = self.encoder(x, lat, lon, time_features)
        h = self.processor(h, edge_index, edge_attr)
        residual = self.decoder(h)
        return x + residual

    def rollout(
        self,
        x: torch.Tensor,
        edge_index: torch.Tensor,
        lat: torch.Tensor,
        lon: torch.Tensor,
        n_steps: int,
        time_features: Optional[torch.Tensor] = None,
        edge_attr: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Autoregressive multi-step rollout.

        Args:
            n_steps: number of 6-hour steps to predict

        Returns:
            (n_steps, N, F) predicted states at each step
        """
        predictions = []
        current = x
        for step in range(n_steps):
            current = self.forward(
                current, edge_index, lat, lon, time_features, edge_attr
            )
            predictions.append(current)
        return torch.stack(predictions, dim=0)

    def compute_loss(
        self,
        pred: torch.Tensor,
        target: torch.Tensor,
        lat: torch.Tensor,
        meta: Optional[Dict] = None,
    ) -> Dict[str, torch.Tensor]:
        """
        Compute weighted MSE loss with latitude weighting and physics constraints.

        Latitude weighting: cos(lat) to account for grid cell area variation.
        """
        # Latitude-weighted MSE
        cos_lat = torch.cos(torch.deg2rad(lat)).unsqueeze(-1)
        weights = cos_lat / cos_lat.mean()

        # Per-variable precision weighting (learned)
        precision = torch.exp(-self.log_var_weights)
        mse = ((pred - target) ** 2 * weights * precision.unsqueeze(0)).mean()
        reg = self.log_var_weights.mean()  # regularization on uncertainty
        data_loss = mse + reg

        losses = {'data_loss': data_loss, 'mse': mse}

        # Physics losses (simplified — use subset of variables)
        if meta is not None:
            physics = self.physics_loss(
                pred={'q': pred[:, -13:]},  # last 13 = humidity across levels
                target={'q': target[:, -13:]},
                meta=meta,
            )
            losses.update(physics)
            losses['total'] = data_loss + physics.get('total_physics', 0)
        else:
            losses['total'] = data_loss

        return losses

    def get_model_size(self) -> Dict[str, int]:
        """Return model parameter counts."""
        total = sum(p.numel() for p in self.parameters())
        trainable = sum(p.numel() for p in self.parameters() if p.requires_grad)
        return {
            'total_parameters': total,
            'trainable_parameters': trainable,
            'encoder_params': sum(p.numel() for p in self.encoder.parameters()),
            'processor_params': sum(p.numel() for p in self.processor.parameters()),
            'decoder_params': sum(p.numel() for p in self.decoder.parameters()),
        }
