from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from . import FIGURES_DIR, RESULTS_DIR, SEED

ROTATION_CSV = RESULTS_DIR / "rotation_analysis.csv"
ROTATION_FIG = FIGURES_DIR / "debris_rotation.png"


plt.style.use("seaborn-v0_8-whitegrid")


def cylinder_inertia(mass_kg: float, radius_m: float, length_m: float) -> tuple[float, float, float]:
    i_axial = 0.5 * mass_kg * radius_m**2
    i_transverse = (1.0 / 12.0) * mass_kg * (3 * radius_m**2 + length_m**2)
    return i_transverse, i_transverse, i_axial


def euler_torque_free(_: float, omega: np.ndarray, inertia: tuple[float, float, float]) -> np.ndarray:
    i1, i2, i3 = inertia
    w1, w2, w3 = omega
    return np.array(
        [
            ((i2 - i3) / i1) * w2 * w3,
            ((i3 - i1) / i2) * w3 * w1,
            ((i1 - i2) / i3) * w1 * w2,
        ]
    )


def analyze_rotation() -> tuple[pd.DataFrame, dict[str, float]]:
    rng = np.random.default_rng(SEED)
    inertia = cylinder_inertia(mass_kg=750.0, radius_m=1.2, length_m=6.5)
    t_eval = np.linspace(0.0, 7200.0, 2400)
    initial_omega_deg_s = np.array([2.5, 3.8, 12.0])
    initial_omega = np.deg2rad(initial_omega_deg_s)

    sol = solve_ivp(
        lambda t, y: euler_torque_free(t, y, inertia),
        (t_eval[0], t_eval[-1]),
        initial_omega,
        t_eval=t_eval,
        rtol=1e-9,
        atol=1e-9,
    )

    omega_deg_s = np.rad2deg(sol.y.T)
    spin_mag = np.linalg.norm(omega_deg_s, axis=1)
    true_period_s = 360.0 / np.median(np.abs(omega_deg_s[:, 2]))
    light_curve = 1.0 + 0.25 * np.sin(2 * np.pi * t_eval / true_period_s) + 0.08 * np.sin(4 * np.pi * t_eval / true_period_s + 0.6)
    light_curve += rng.normal(0.0, 0.05, len(t_eval))

    dt = t_eval[1] - t_eval[0]
    centered = light_curve - light_curve.mean()
    fft_vals = np.fft.rfft(centered)
    freqs = np.fft.rfftfreq(len(centered), d=dt)
    positive = freqs > 0
    peak_idx = np.argmax(np.abs(fft_vals[positive]))
    dominant_freq = freqs[positive][peak_idx]
    estimated_period_s = 1.0 / dominant_freq

    window_mask = spin_mag < 5.0
    capture_windows = []
    if np.any(window_mask):
        indices = np.where(window_mask)[0]
        start = indices[0]
        prev = indices[0]
        for idx in indices[1:]:
            if idx != prev + 1:
                capture_windows.append((t_eval[start], t_eval[prev]))
                start = idx
            prev = idx
        capture_windows.append((t_eval[start], t_eval[prev]))

    df = pd.DataFrame(
        {
            "time_s": t_eval,
            "omega1_deg_s": omega_deg_s[:, 0],
            "omega2_deg_s": omega_deg_s[:, 1],
            "omega3_deg_s": omega_deg_s[:, 2],
            "spin_rate_deg_s": spin_mag,
            "light_curve": light_curve,
        }
    )
    df.to_csv(ROTATION_CSV, index=False)

    fig, axes = plt.subplots(3, 1, figsize=(11, 12), constrained_layout=True)
    axes[0].plot(t_eval / 60.0, omega_deg_s[:, 0], label="ω1", color="#1f77b4")
    axes[0].plot(t_eval / 60.0, omega_deg_s[:, 1], label="ω2", color="#ff7f0e")
    axes[0].plot(t_eval / 60.0, omega_deg_s[:, 2], label="ω3", color="#2ca02c")
    axes[0].axhline(5.0, color="crimson", linestyle="--", linewidth=1.0, label="Capture threshold")
    axes[0].set_xlabel("Time (min)")
    axes[0].set_ylabel("Angular velocity (deg/s)")
    axes[0].set_title("Torque-free tumbling dynamics")
    axes[0].legend()

    axes[1].plot(t_eval / 60.0, light_curve, color="#4c72b0")
    axes[1].set_xlabel("Time (min)")
    axes[1].set_ylabel("Normalized brightness")
    axes[1].set_title("Simulated debris light curve")

    axes[2].plot(freqs[positive] * 1000.0, np.abs(fft_vals[positive]), color="#55a868")
    axes[2].axvline(dominant_freq * 1000.0, color="crimson", linestyle="--", linewidth=1.0)
    axes[2].set_xlabel("Frequency (mHz)")
    axes[2].set_ylabel("FFT magnitude")
    axes[2].set_title("Rotation period estimation via FFT")
    fig.savefig(ROTATION_FIG, dpi=300)
    plt.close(fig)

    summary = {
        "estimated_period_s": float(estimated_period_s),
        "true_period_s": float(true_period_s),
        "mean_spin_rate_deg_s": float(spin_mag.mean()),
        "capture_windows_count": int(len(capture_windows)),
        "longest_capture_window_s": float(max((end - start) for start, end in capture_windows) if capture_windows else 0.0),
    }
    (RESULTS_DIR / "rotation_summary.json").write_text(json.dumps({"summary": summary, "capture_windows_s": capture_windows}, indent=2))
    return df, summary


if __name__ == "__main__":
    _, summary = analyze_rotation()
    print(json.dumps(summary, indent=2))
