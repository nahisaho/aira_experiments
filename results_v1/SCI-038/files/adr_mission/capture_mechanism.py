from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import FIGURES_DIR, RESULTS_DIR

CAPTURE_CSV = RESULTS_DIR / "capture_analysis.csv"
CAPTURE_FIG = FIGURES_DIR / "capture_mechanisms.png"


plt.style.use("seaborn-v0_8-whitegrid")


def forward_kinematics(theta1: np.ndarray, theta2: np.ndarray, l1: float = 2.5, l2: float = 2.0) -> tuple[np.ndarray, np.ndarray]:
    x = l1 * np.cos(theta1) + l2 * np.cos(theta1 + theta2)
    y = l1 * np.sin(theta1) + l2 * np.sin(theta1 + theta2)
    return x, y


def robotic_arm_success(rotation_deg_s: np.ndarray) -> np.ndarray:
    end_effector_speed = np.deg2rad(rotation_deg_s) * 1.5
    contact_force_proxy = np.exp(-0.18 * end_effector_speed**2)
    workspace_margin = 0.92
    return np.clip(workspace_margin * contact_force_proxy, 0.0, 1.0)


def net_spread_radius(time_s: np.ndarray, deployment_speed: float = 4.5) -> np.ndarray:
    return deployment_speed * time_s * np.exp(-0.04 * time_s)


def net_success(rotation_deg_s: np.ndarray) -> np.ndarray:
    spread_time = 1.8
    radius = net_spread_radius(np.array([spread_time]))[0]
    geometric_factor = 1.0 - np.exp(-radius / 3.5)
    rotation_penalty = np.exp(-0.06 * rotation_deg_s)
    return np.clip(0.88 * geometric_factor * rotation_penalty + 0.08, 0.0, 1.0)


def harpoon_success(rotation_deg_s: np.ndarray) -> np.ndarray:
    relative_angle_error = np.deg2rad(rotation_deg_s) * 0.12
    penetration = 1.0 / (1.0 + np.exp(5.0 * (relative_angle_error - 0.18)))
    ballistic_factor = np.exp(-0.015 * rotation_deg_s)
    return np.clip(0.82 * penetration * ballistic_factor, 0.0, 1.0)


def analyze_capture_mechanisms() -> tuple[pd.DataFrame, dict[str, float]]:
    rotation_rates = np.linspace(0.0, 30.0, 121)
    arm_prob = robotic_arm_success(rotation_rates)
    net_prob = net_success(rotation_rates)
    harpoon_prob = harpoon_success(rotation_rates)

    theta1 = np.linspace(-np.pi / 2, np.pi / 2, 60)
    theta2 = np.linspace(-np.pi, np.pi, 60)
    th1, th2 = np.meshgrid(theta1, theta2)
    workspace_x, workspace_y = forward_kinematics(th1, th2)
    workspace_radius = np.sqrt(workspace_x**2 + workspace_y**2)

    df = pd.DataFrame(
        {
            "rotation_rate_deg_s": rotation_rates,
            "robotic_arm_success": arm_prob,
            "net_success": net_prob,
            "harpoon_success": harpoon_prob,
        }
    )
    df.to_csv(CAPTURE_CSV, index=False)

    fig, ax = plt.subplots(figsize=(10, 7), constrained_layout=True)
    ax.plot(rotation_rates, arm_prob, label="Robotic arm", color="#1f77b4", linewidth=2.2)
    ax.plot(rotation_rates, net_prob, label="Net", color="#ff7f0e", linewidth=2.2)
    ax.plot(rotation_rates, harpoon_prob, label="Harpoon", color="#2ca02c", linewidth=2.2)
    ax.set_xlabel("Debris rotation rate (deg/s)")
    ax.set_ylabel("Capture success probability")
    ax.set_title("Capture mechanism performance versus debris rotation")
    ax.legend()
    fig.savefig(CAPTURE_FIG, dpi=300)
    plt.close(fig)

    summary = {
        "workspace_radius_max_m": float(workspace_radius.max()),
        "robotic_arm_threshold_deg_s": float(rotation_rates[np.argmin(np.abs(arm_prob - 0.5))]),
        "net_threshold_deg_s": float(rotation_rates[np.argmin(np.abs(net_prob - 0.5))]),
        "harpoon_threshold_deg_s": float(rotation_rates[np.argmin(np.abs(harpoon_prob - 0.5))]),
    }
    (RESULTS_DIR / "capture_summary.json").write_text(json.dumps(summary, indent=2))
    return df, summary


if __name__ == "__main__":
    _, summary = analyze_capture_mechanisms()
    print(json.dumps(summary, indent=2))
