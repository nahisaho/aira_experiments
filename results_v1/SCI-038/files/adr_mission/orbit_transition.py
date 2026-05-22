from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import FIGURES_DIR, GM_EARTH, R_EARTH, RESULTS_DIR

DELTA_V_MATRIX_PATH = RESULTS_DIR / "delta_v_matrix.csv"
HEATMAP_PATH = FIGURES_DIR / "delta_v_heatmap.png"


plt.style.use("seaborn-v0_8-whitegrid")


def circular_velocity(radius_m: float) -> float:
    return np.sqrt(GM_EARTH / radius_m)


def edelbaum_delta_v(a1_km: float, inc1_deg: float, a2_km: float, inc2_deg: float) -> tuple[float, float]:
    r1 = a1_km * 1e3
    r2 = a2_km * 1e3
    v1 = circular_velocity(r1)
    v2 = circular_velocity(r2)
    delta_i = np.deg2rad(abs(inc2_deg - inc1_deg))
    delta_v = np.sqrt(v1**2 + v2**2 - 2.0 * v1 * v2 * np.cos(0.5 * np.pi * delta_i))
    return float(delta_v), float(delta_i)


def compute_transfer_time(delta_v_m_s: float, thrust_acceleration: float = 1e-4, isp_s: float = 3000.0) -> float:
    _ = isp_s
    return float(delta_v_m_s / thrust_acceleration / 86400.0)


def build_delta_v_matrix(targets: pd.DataFrame, initial_orbit: dict[str, float] | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    if initial_orbit is None:
        initial_orbit = {"name": "CHASER-INITIAL", "semi_major_axis_km": 6971.0, "inclination_deg": 97.4}

    nodes = pd.concat(
        [
            pd.DataFrame([initial_orbit]),
            targets[["debris_id", "semi_major_axis_km", "inclination_deg"]].rename(columns={"debris_id": "name"}),
        ],
        ignore_index=True,
    )

    n = len(nodes)
    delta_v = np.zeros((n, n))
    transfer_days = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            dv, _ = edelbaum_delta_v(
                nodes.iloc[i]["semi_major_axis_km"],
                nodes.iloc[i]["inclination_deg"],
                nodes.iloc[j]["semi_major_axis_km"],
                nodes.iloc[j]["inclination_deg"],
            )
            delta_v[i, j] = dv
            transfer_days[i, j] = compute_transfer_time(dv)

    dv_df = pd.DataFrame(delta_v, index=nodes["name"], columns=nodes["name"])
    dv_df.to_csv(DELTA_V_MATRIX_PATH)
    time_df = pd.DataFrame(transfer_days, index=nodes["name"], columns=nodes["name"])
    time_df.to_csv(RESULTS_DIR / "transfer_time_matrix_days.csv")

    fig, ax = plt.subplots(figsize=(9, 8), constrained_layout=True)
    im = ax.imshow(dv_df.values, cmap="cividis")
    ax.set_xticks(range(n), dv_df.columns, rotation=90)
    ax.set_yticks(range(n), dv_df.index)
    ax.set_title("Pairwise low-thrust transfer ΔV matrix")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("ΔV (m/s)")
    fig.savefig(HEATMAP_PATH, dpi=300)
    plt.close(fig)

    summary = {
        "max_delta_v_m_s": float(dv_df.values.max()),
        "mean_delta_v_m_s": float(dv_df.values[dv_df.values > 0].mean()),
        "max_transfer_days": float(time_df.values.max()),
    }
    (RESULTS_DIR / "orbit_transition_summary.json").write_text(json.dumps(summary, indent=2))
    return dv_df, time_df


if __name__ == "__main__":
    from .debris_catalog import generate_debris_catalog
    from .target_selection import score_targets

    catalog = generate_debris_catalog()
    _, top_targets = score_targets(catalog)
    dv_df, tt_df = build_delta_v_matrix(top_targets)
    print(dv_df.head().to_string())
