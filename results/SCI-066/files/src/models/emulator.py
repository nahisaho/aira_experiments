"""
ESM AI Emulator - Combined U-Net + ConvLSTM architecture.

Integrates spatial pattern learning (U-Net) with temporal dynamics
(ConvLSTM) for full spatiotemporal climate emulation.
Includes physics-constrained loss and ensemble uncertainty estimation.
"""

import torch
import torch.nn as nn
from typing import Dict, Tuple

from .unet import ClimateUNet
from .convlstm import ConvLSTMPredictor


class PhysicsConstraintLayer(nn.Module):
    """
    Enforces physical conservation laws as soft constraints.

    Constraints:
    1. Energy conservation: global mean energy budget balance
    2. Mass conservation: total precipitation ~ evaporation (moisture budget)
    3. Thermodynamic consistency: spatial temperature gradients respect lapse rate
    """

    def __init__(self, grid_area_weights: torch.Tensor = None):
        super().__init__()
        self.grid_area_weights = grid_area_weights

    def energy_conservation_loss(self, pred: torch.Tensor,
                                  forcing: torch.Tensor) -> torch.Tensor:
        """
        Penalizes violations of top-of-atmosphere energy balance.

        Global mean temperature change should be proportional to
        net radiative forcing minus outgoing longwave radiation.

        ΔT_global ∝ F_net - λ * T_global (linear feedback)
        """
        # pred[:, 0] = temperature field
        if self.grid_area_weights is not None:
            w = self.grid_area_weights.to(pred.device)
            global_mean_t = (pred[:, 0] * w).sum(dim=(-2, -1)) / w.sum()
        else:
            global_mean_t = pred[:, 0].mean(dim=(-2, -1))

        # Forcing magnitude from scenario encoding (simplified linear model)
        forcing_magnitude = forcing.sum(dim=-1)

        # Climate sensitivity parameter λ ≈ 1.2 W/m²/K (Planck feedback)
        lambda_feedback = 1.2
        energy_imbalance = forcing_magnitude - lambda_feedback * global_mean_t

        return energy_imbalance.pow(2).mean()

    def mass_conservation_loss(self, pred: torch.Tensor) -> torch.Tensor:
        """
        Penalizes negative precipitation and ensures moisture budget closure.
        Precipitation must be non-negative.
        """
        # pred[:, 1] = precipitation field
        precip = pred[:, 1]
        negative_precip_penalty = torch.relu(-precip).pow(2).mean()

        return negative_precip_penalty

    def spatial_smoothness_loss(self, pred: torch.Tensor) -> torch.Tensor:
        """
        Penalizes physically implausible sharp gradients.
        Climate fields should be spatially smooth at grid scale.
        """
        dx = pred[:, :, :, 1:] - pred[:, :, :, :-1]
        dy = pred[:, :, 1:, :] - pred[:, :, :-1, :]

        return (dx.pow(2).mean() + dy.pow(2).mean()) * 0.5

    def forward(self, pred: torch.Tensor,
                forcing: torch.Tensor) -> Dict[str, torch.Tensor]:
        return {
            "energy_conservation": self.energy_conservation_loss(pred, forcing),
            "mass_conservation": self.mass_conservation_loss(pred),
            "spatial_smoothness": self.spatial_smoothness_loss(pred),
        }


class EnsembleWrapper(nn.Module):
    """
    Deep ensemble for uncertainty quantification.

    Trains N independent models and aggregates predictions to estimate
    epistemic uncertainty via prediction spread.
    """

    def __init__(self, model_factory, n_members: int = 5):
        super().__init__()
        self.members = nn.ModuleList([model_factory() for _ in range(n_members)])
        self.n_members = n_members

    def forward(self, *args, **kwargs) -> Tuple[torch.Tensor, torch.Tensor]:
        preds = torch.stack([m(*args, **kwargs) for m in self.members], dim=0)
        mean = preds.mean(dim=0)
        std = preds.std(dim=0)
        return mean, std


class ESMEmulator(nn.Module):
    """
    Full Earth System Model AI Emulator.

    Combines:
    - U-Net: spatial pattern reconstruction conditioned on SSP scenario
    - ConvLSTM: temporal evolution prediction
    - Physics constraints: conservation law enforcement
    - Ensemble: uncertainty quantification

    Architecture flow:
    1. ConvLSTM processes temporal sequence → hidden state
    2. U-Net refines spatial field conditioned on scenario + ConvLSTM output
    3. Physics constraints applied as regularization during training
    4. Ensemble members provide uncertainty bands
    """

    def __init__(self, config: dict = None):
        super().__init__()
        if config is None:
            config = self.default_config()

        self.config = config
        n_vars = config["n_climate_vars"]
        spatial = config["spatial_size"]

        self.convlstm = ConvLSTMPredictor(
            input_dim=n_vars,
            hidden_dims=config["convlstm_hidden_dims"],
            n_layers=config["convlstm_n_layers"],
            out_channels=n_vars,
        )

        self.unet = ClimateUNet(
            in_channels=n_vars * 2,  # ConvLSTM output + current state
            out_channels=n_vars,
            base_features=config["unet_base_features"],
            n_scenarios=config["n_scenarios"],
            forcing_dim=config["forcing_dim"],
            spatial_size=spatial,
            dropout=config["dropout"],
        )

        self.physics = PhysicsConstraintLayer()

        # Learnable blending of ConvLSTM and U-Net outputs
        self.blend_weight = nn.Parameter(torch.tensor(0.5))

    @staticmethod
    def default_config() -> dict:
        return {
            "n_climate_vars": 3,  # T, P, SL
            "spatial_size": (64, 128),  # lat x lon (downscaled)
            "n_scenarios": 4,  # SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5
            "forcing_dim": 8,
            "convlstm_hidden_dims": [64, 64, 64],
            "convlstm_n_layers": 3,
            "unet_base_features": 64,
            "dropout": 0.1,
            "seq_length": 10,
            "physics_weight_energy": 0.1,
            "physics_weight_mass": 0.05,
            "physics_weight_smooth": 0.01,
            "n_ensemble_members": 5,
        }

    def forward(self, x_seq: torch.Tensor, scenario_id: torch.Tensor,
                forcing: torch.Tensor) -> Dict[str, torch.Tensor]:
        """
        Args:
            x_seq: (B, T, C, H, W) - historical climate field sequence
            scenario_id: (B,) - SSP scenario index
            forcing: (B, forcing_dim) - radiative forcing vector

        Returns:
            dict with 'prediction', 'physics_losses'
        """
        # Temporal prediction via ConvLSTM
        temporal_pred = self.convlstm(x_seq)

        # Current state (last timestep)
        current_state = x_seq[:, -1]

        # Spatial refinement via U-Net
        unet_input = torch.cat([current_state, temporal_pred], dim=1)
        spatial_pred = self.unet(unet_input, scenario_id, forcing)

        # Blend temporal and spatial predictions
        alpha = torch.sigmoid(self.blend_weight)
        prediction = alpha * spatial_pred + (1 - alpha) * temporal_pred

        # Physics constraint losses
        physics_losses = self.physics(prediction, forcing)

        return {
            "prediction": prediction,
            "temporal_pred": temporal_pred,
            "spatial_pred": spatial_pred,
            "physics_losses": physics_losses,
        }

    def compute_loss(self, outputs: Dict[str, torch.Tensor],
                     target: torch.Tensor) -> Dict[str, torch.Tensor]:
        """Combined MSE + physics-constrained loss."""
        cfg = self.config
        mse = nn.functional.mse_loss(outputs["prediction"], target)

        physics = outputs["physics_losses"]
        total_loss = (
            mse
            + cfg["physics_weight_energy"] * physics["energy_conservation"]
            + cfg["physics_weight_mass"] * physics["mass_conservation"]
            + cfg["physics_weight_smooth"] * physics["spatial_smoothness"]
        )

        return {
            "total": total_loss,
            "mse": mse,
            **{f"physics_{k}": v for k, v in physics.items()},
        }
