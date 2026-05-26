"""
Evaluation metrics for weather prediction.
RMSE, ACC (Anomaly Correlation Coefficient), and physics constraint metrics.
"""

import torch
import numpy as np


VARIABLE_NAMES = ['Temperature', 'U-wind', 'V-wind', 'Specific Humidity']
SURFACE_NAMES = ['2m Temperature', '10m U-wind', '10m V-wind', 'MSLP']
PRESSURE_LEVELS = [50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000]


def compute_rmse(pred, target):
    """Root Mean Square Error."""
    return torch.sqrt(torch.mean((pred - target) ** 2)).item()


def compute_rmse_per_variable(pred_p, target_p, pred_s, target_s):
    """RMSE per variable and level."""
    results = {}
    for i, name in enumerate(VARIABLE_NAMES):
        level_rmse = []
        for k in range(len(PRESSURE_LEVELS)):
            rmse = compute_rmse(pred_p[:, :, i, k], target_p[:, :, i, k])
            level_rmse.append(rmse)
        results[name] = {
            'levels': dict(zip(PRESSURE_LEVELS, level_rmse)),
            'mean': np.mean(level_rmse)
        }
    for i, name in enumerate(SURFACE_NAMES):
        results[name] = compute_rmse(pred_s[:, :, i], target_s[:, :, i])
    return results


def compute_acc(pred, target, climatology):
    """Anomaly Correlation Coefficient."""
    pred_anom = pred - climatology
    target_anom = target - climatology

    numerator = torch.sum(pred_anom * target_anom)
    denominator = torch.sqrt(
        torch.sum(pred_anom ** 2) * torch.sum(target_anom ** 2) + 1e-10
    )
    return (numerator / denominator).item()


def compute_acc_per_variable(pred_p, target_p, clim_p):
    """ACC per variable."""
    results = {}
    for i, name in enumerate(VARIABLE_NAMES):
        acc = compute_acc(pred_p[:, :, i, :], target_p[:, :, i, :], clim_p[:, :, i, :])
        results[name] = acc
    return results


def compute_skill_score(model_rmse, baseline_rmse):
    """Skill score relative to baseline (1 = perfect, 0 = same as baseline)."""
    return 1.0 - model_rmse / (baseline_rmse + 1e-10)


def compute_physics_metrics(pred_p, input_p):
    """Compute mass and energy conservation metrics."""
    dp = torch.tensor(
        [50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000],
        dtype=torch.float32
    ) * 100
    # Layer thicknesses
    dp_weights = torch.zeros(13)
    for i in range(13):
        if i == 0:
            dp_weights[i] = (dp[1] - dp[0]) / 2
        elif i == 12:
            dp_weights[i] = (dp[-1] - dp[-2]) / 2
        else:
            dp_weights[i] = (dp[i + 1] - dp[i - 1]) / 2

    g = 9.80665

    # Column moisture (precipitable water)
    q_pred = pred_p[:, :, 3, :]
    q_input = input_p[:, :, 3, :]
    pw_pred = (q_pred * dp_weights).sum(dim=-1) / g
    pw_input = (q_input * dp_weights).sum(dim=-1) / g
    mass_error = torch.abs(pw_pred - pw_input).mean().item()

    # Column energy
    T_pred, T_input = pred_p[:, :, 0, :], input_p[:, :, 0, :]
    u_pred, u_input = pred_p[:, :, 1, :], input_p[:, :, 1, :]
    v_pred, v_input = pred_p[:, :, 2, :], input_p[:, :, 2, :]

    cp = 1004.0
    KE_pred = 0.5 * (u_pred ** 2 + v_pred ** 2)
    KE_input = 0.5 * (u_input ** 2 + v_input ** 2)
    E_pred = ((cp * T_pred + KE_pred) * dp_weights).sum(dim=-1) / g
    E_input = ((cp * T_input + KE_input) * dp_weights).sum(dim=-1) / g
    energy_error = torch.abs(E_pred - E_input).mean().item()
    energy_rel_error = (torch.abs(E_pred - E_input) / (torch.abs(E_input) + 1e-10)).mean().item()

    return {
        'mass_error_kg_m2': mass_error,
        'energy_error_J_m2': energy_error,
        'energy_relative_error': energy_rel_error,
    }
