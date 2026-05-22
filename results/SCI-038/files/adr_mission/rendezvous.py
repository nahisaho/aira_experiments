from __future__ import annotations

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import FIGURES_DIR, GM_EARTH, R_EARTH, RESULTS_DIR

RENDEZVOUS_CSV = RESULTS_DIR / "rendezvous_trajectories.csv"
RENDEZVOUS_FIG = FIGURES_DIR / "rendezvous_trajectory.png"


plt.style.use("seaborn-v0_8-whitegrid")


def mean_motion(altitude_km: float = 700.0) -> float:
    return np.sqrt(GM_EARTH / (R_EARTH + altitude_km * 1e3) ** 3)


def cw_rr(n: float, t: np.ndarray) -> np.ndarray:
    ct = np.cos(n * t)
    st = np.sin(n * t)
    mats = np.zeros((len(t), 3, 3))
    mats[:, 0, 0] = 4 - 3 * ct
    mats[:, 1, 0] = 6 * (st - n * t)
    mats[:, 1, 1] = 1
    mats[:, 2, 2] = ct
    return mats


def cw_rv(n: float, t: np.ndarray) -> np.ndarray:
    ct = np.cos(n * t)
    st = np.sin(n * t)
    mats = np.zeros((len(t), 3, 3))
    mats[:, 0, 0] = st / n
    mats[:, 0, 1] = 2 * (1 - ct) / n
    mats[:, 1, 0] = 2 * (ct - 1) / n
    mats[:, 1, 1] = (4 * st - 3 * n * t) / n
    mats[:, 2, 2] = st / n
    return mats


def cw_vr(n: float, t: np.ndarray) -> np.ndarray:
    ct = np.cos(n * t)
    st = np.sin(n * t)
    mats = np.zeros((len(t), 3, 3))
    mats[:, 0, 0] = 3 * n * st
    mats[:, 1, 0] = 6 * n * (ct - 1)
    mats[:, 2, 2] = -n * st
    return mats


def cw_vv(n: float, t: np.ndarray) -> np.ndarray:
    ct = np.cos(n * t)
    st = np.sin(n * t)
    mats = np.zeros((len(t), 3, 3))
    mats[:, 0, 0] = ct
    mats[:, 0, 1] = 2 * st
    mats[:, 1, 0] = -2 * st
    mats[:, 1, 1] = 4 * ct - 3
    mats[:, 2, 2] = ct
    return mats


def propagate_cw(r0: np.ndarray, v0: np.ndarray, n: float, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rr = cw_rr(n, t)
    rv = cw_rv(n, t)
    vr = cw_vr(n, t)
    vv = cw_vv(n, t)
    r = np.einsum("tij,j->ti", rr, r0) + np.einsum("tij,j->ti", rv, v0)
    v = np.einsum("tij,j->ti", vr, r0) + np.einsum("tij,j->ti", vv, v0)
    return r, v


def two_impulse_rendezvous(r0: np.ndarray, tf: float, n: float, drift_rate_m_s: float = 0.0) -> dict[str, np.ndarray | float]:
    rr_tf = cw_rr(n, np.array([tf]))[0]
    rv_tf = cw_rv(n, np.array([tf]))[0]
    v0_req = -np.linalg.solve(rv_tf, rr_tf @ r0)
    v0_req[1] += drift_rate_m_s

    t = np.linspace(0.0, tf, 400)
    r, v = propagate_cw(r0, v0_req, n, t)
    vf = v[-1]
    delta_v1 = np.linalg.norm(v0_req)
    delta_v2 = np.linalg.norm(vf)
    return {
        "t": t,
        "r": r,
        "v": v,
        "delta_v1": float(delta_v1),
        "delta_v2": float(delta_v2),
        "delta_v_total": float(delta_v1 + delta_v2),
    }


def simulate_rendezvous_scenarios() -> pd.DataFrame:
    n = mean_motion(700.0)
    scenarios = [
        {"name": "Scenario-A", "r0": np.array([0.0, -5000.0, 150.0]), "tf": 3.5 * 3600.0, "drift": -0.02},
        {"name": "Scenario-B", "r0": np.array([200.0, -5000.0, -120.0]), "tf": 4.0 * 3600.0, "drift": 0.00},
        {"name": "Scenario-C", "r0": np.array([-180.0, -5000.0, 90.0]), "tf": 4.5 * 3600.0, "drift": 0.03},
    ]

    frames = []
    summaries = []
    fig = plt.figure(figsize=(11, 8), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    cmap = plt.get_cmap("tab10")

    for idx, scenario in enumerate(scenarios):
        result = two_impulse_rendezvous(scenario["r0"], scenario["tf"], n, scenario["drift"])
        r = result["r"]
        t = result["t"]
        df = pd.DataFrame(
            {
                "scenario": scenario["name"],
                "time_s": t,
                "x_m": r[:, 0],
                "y_m": r[:, 1],
                "z_m": r[:, 2],
                "range_m": np.linalg.norm(r, axis=1),
                "delta_v_total_m_s": result["delta_v_total"],
            }
        )
        frames.append(df)
        summaries.append(
            {
                "scenario": scenario["name"],
                "delta_v1_m_s": result["delta_v1"],
                "delta_v2_m_s": result["delta_v2"],
                "delta_v_total_m_s": result["delta_v_total"],
                "final_range_m": float(np.linalg.norm(r[-1])),
            }
        )
        ax.plot(r[:, 1], r[:, 0], r[:, 2], color=cmap(idx), label=scenario["name"])
        ax.scatter(r[0, 1], r[0, 0], r[0, 2], color=cmap(idx), marker="o")
        ax.scatter(r[-1, 1], r[-1, 0], r[-1, 2], color=cmap(idx), marker="^")

    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(RENDEZVOUS_CSV, index=False)
    (RESULTS_DIR / "rendezvous_summary.json").write_text(json.dumps(summaries, indent=2))

    ax.set_xlabel("Along-track y (m)")
    ax.set_ylabel("Radial x (m)")
    ax.set_zlabel("Cross-track z (m)")
    ax.set_title("CW rendezvous trajectories for V-bar approach scenarios")
    ax.legend()
    fig.savefig(RENDEZVOUS_FIG, dpi=300)
    plt.close(fig)
    return combined


if __name__ == "__main__":
    df = simulate_rendezvous_scenarios()
    print(df.groupby("scenario")["range_m"].agg(["min", "max"]).to_string())
