"""
physics.py — Physical consistency constraints for weather prediction.

Implements differentiable physics constraints for:
  1. Mass conservation (continuity equation)
  2. Energy conservation (total atmospheric energy budget)
  3. Moisture conservation (non-negative humidity)
  4. Hydrostatic balance approximation
"""

import torch
import torch.nn as nn
import numpy as np
from typing import Dict, Tuple

# Physical constants
R_d = 287.05    # Dry air gas constant (J/kg/K)
C_p = 1004.0    # Specific heat at constant pressure (J/kg/K)
g = 9.81        # Gravitational acceleration (m/s²)
L_v = 2.5e6     # Latent heat of vaporization (J/kg)
R_earth = 6.371e6  # Earth radius (m)

PRESSURE_LEVELS_PA = torch.tensor([
    5000, 10000, 15000, 20000, 25000, 30000,
    40000, 50000, 60000, 70000, 85000, 92500, 100000
], dtype=torch.float32)


class MassConservationLoss(nn.Module):
    """
    Enforces approximate mass conservation via the continuity equation.

    In pressure coordinates, mass conservation requires:
        ∂u/∂x + ∂v/∂y + ∂ω/∂p ≈ 0

    We approximate this as: the column-integrated divergence should be near zero.
    """

    def __init__(self, n_levels: int = 13):
        super().__init__()
        self.n_levels = n_levels

    def forward(
        self,
        u_pred: torch.Tensor,
        v_pred: torch.Tensor,
        lat: torch.Tensor,
        lon: torch.Tensor,
        dp: torch.Tensor,
    ) -> torch.Tensor:
        """
        Args:
            u_pred: (N, n_levels) predicted u-wind
            v_pred: (N, n_levels) predicted v-wind
            lat: (N,) latitude in radians
            lon: (N,) longitude in radians
            dp: (n_levels,) pressure thickness per level

        Returns:
            Scalar mass conservation loss
        """
        cos_lat = torch.cos(lat).unsqueeze(-1).clamp(min=1e-6)

        # Column-integrated mass flux (simplified)
        mass_flux_u = (u_pred * dp.unsqueeze(0) / (R_earth * cos_lat)).sum(dim=-1)
        mass_flux_v = (v_pred * dp.unsqueeze(0) / R_earth).sum(dim=-1)

        # Divergence should be near zero globally
        divergence = mass_flux_u + mass_flux_v
        return torch.mean(divergence ** 2)


class EnergyConservationLoss(nn.Module):
    """
    Enforces approximate energy conservation.

    Total atmospheric energy per unit mass:
        E = C_p * T + g * z + 0.5 * (u² + v²) + L_v * q

    The global integral of E should be approximately conserved.
    """

    def __init__(self):
        super().__init__()

    def forward(
        self,
        t_pred: torch.Tensor,
        t_input: torch.Tensor,
        u_pred: torch.Tensor,
        u_input: torch.Tensor,
        v_pred: torch.Tensor,
        v_input: torch.Tensor,
        q_pred: torch.Tensor,
        q_input: torch.Tensor,
        z_pred: torch.Tensor,
        z_input: torch.Tensor,
        weights: torch.Tensor,
    ) -> torch.Tensor:
        """
        Compute energy conservation loss as normalized difference
        in global mean total energy.
        """
        def total_energy(t, u, v, q, z):
            thermal = C_p * t
            kinetic = 0.5 * (u**2 + v**2)
            potential = z  # geopotential already in m²/s²
            latent = L_v * q
            return (thermal + kinetic + potential + latent).mean(dim=-1)

        e_pred = total_energy(t_pred, u_pred, v_pred, q_pred, z_pred)
        e_input = total_energy(t_input, u_input, v_input, q_input, z_input)

        # Weighted global mean energy
        e_pred_global = (e_pred * weights).sum() / weights.sum()
        e_input_global = (e_input * weights).sum() / weights.sum()

        # Relative energy change should be small
        relative_change = (e_pred_global - e_input_global) / (e_input_global.abs() + 1e-8)
        return relative_change ** 2


class MoistureConstraint(nn.Module):
    """Enforce non-negative specific humidity with soft penalty."""

    def forward(self, q_pred: torch.Tensor) -> torch.Tensor:
        return torch.mean(torch.relu(-q_pred) ** 2)


class HydrostaticBalance(nn.Module):
    """
    Approximate hydrostatic balance: ∂Φ/∂p ≈ -R_d * T_v / p
    where T_v = T * (1 + 0.608 * q) is virtual temperature.
    """

    def __init__(self):
        super().__init__()

    def forward(
        self,
        z_pred: torch.Tensor,
        t_pred: torch.Tensor,
        q_pred: torch.Tensor,
        pressure_levels: torch.Tensor,
    ) -> torch.Tensor:
        """
        Args:
            z_pred: (N, n_levels) geopotential
            t_pred: (N, n_levels) temperature
            q_pred: (N, n_levels) specific humidity
            pressure_levels: (n_levels,) pressure in Pa
        """
        t_v = t_pred * (1 + 0.608 * q_pred)

        # Finite difference approximation
        dz = z_pred[:, 1:] - z_pred[:, :-1]
        dp = pressure_levels[1:] - pressure_levels[:-1]
        p_mid = 0.5 * (pressure_levels[1:] + pressure_levels[:-1])
        t_v_mid = 0.5 * (t_v[:, 1:] + t_v[:, :-1])

        # Hydrostatic: dz/dp ≈ -R_d * T_v / p
        lhs = dz / dp.unsqueeze(0)
        rhs = -R_d * t_v_mid / p_mid.unsqueeze(0)

        return torch.mean((lhs - rhs) ** 2)


class PhysicsLoss(nn.Module):
    """Combined physics-informed loss with configurable weights."""

    def __init__(
        self,
        w_mass: float = 0.1,
        w_energy: float = 0.1,
        w_moisture: float = 0.05,
        w_hydrostatic: float = 0.05,
    ):
        super().__init__()
        self.w_mass = w_mass
        self.w_energy = w_energy
        self.w_moisture = w_moisture
        self.w_hydrostatic = w_hydrostatic

        self.mass_loss = MassConservationLoss()
        self.energy_loss = EnergyConservationLoss()
        self.moisture_loss = MoistureConstraint()
        self.hydrostatic_loss = HydrostaticBalance()

    def forward(self, pred: Dict, target: Dict, meta: Dict) -> Dict[str, torch.Tensor]:
        """
        Compute all physics losses.

        Returns dict with individual losses and weighted total.
        """
        losses = {}

        # Moisture non-negativity
        if 'q' in pred:
            losses['moisture'] = self.moisture_loss(pred['q']) * self.w_moisture

        # Hydrostatic balance
        if all(k in pred for k in ['z', 't', 'q']):
            losses['hydrostatic'] = self.hydrostatic_loss(
                pred['z'], pred['t'], pred['q'],
                meta.get('pressure_levels', PRESSURE_LEVELS_PA)
            ) * self.w_hydrostatic

        losses['total_physics'] = sum(losses.values())
        return losses
