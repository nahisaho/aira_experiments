#!/usr/bin/env python3
"""Simplified steady axisymmetric CFD simulation for a perfusion bioreactor.

The solver assembles 2D r-z fields by solving a radial finite-difference
Brinkman momentum balance at each axial station under steady perfusion.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

RHO = 1007.0
MU = 0.001
POROSITY = 0.4
PERMEABILITY = 1.0e-10
FORCHHEIMER = 1.0e5
FLOW_RATES_ML_MIN = [0.5, 1.0, 2.0, 5.0]
REPRESENTATIVE_FLOW = 5.0
NR = 90
NZ = 180


def logistic_window(values: np.ndarray, lower: float, upper: float, width: float) -> np.ndarray:
    rise = 1.0 / (1.0 + np.exp(-(values - lower) / width))
    fall = 1.0 / (1.0 + np.exp((values - upper) / width))
    return rise * fall


def load_geometry(root: Path) -> Dict:
    geometry_path = root / "data" / "geometry_params.json"
    if not geometry_path.exists():
        raise FileNotFoundError(
            f"Geometry file not found: {geometry_path}. Run data/bioreactor_geometry.py first."
        )
    with geometry_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def solve_radial_profile(radius: np.ndarray, drag: np.ndarray) -> np.ndarray:
    dr = radius[1] - radius[0]
    n = radius.size
    matrix = sparse.lil_matrix((n, n), dtype=float)
    rhs = np.ones(n, dtype=float)

    matrix[0, 0] = -2.0 * MU / dr**2 - drag[0]
    matrix[0, 1] = 2.0 * MU / dr**2
    rhs[0] = -1.0

    for i in range(1, n - 1):
        r_i = radius[i]
        coeff_minus = MU * (1.0 / dr**2 - 1.0 / (2.0 * r_i * dr))
        coeff_center = -2.0 * MU / dr**2 - drag[i]
        coeff_plus = MU * (1.0 / dr**2 + 1.0 / (2.0 * r_i * dr))
        matrix[i, i - 1] = coeff_minus
        matrix[i, i] = coeff_center
        matrix[i, i + 1] = coeff_plus
        rhs[i] = -1.0

    matrix[-1, -1] = 1.0
    rhs[-1] = 0.0

    profile = spsolve(matrix.tocsr(), rhs)
    if np.any(~np.isfinite(profile)):
        raise RuntimeError("Radial velocity solver produced non-finite values.")
    return profile


def compute_fields(geometry: Dict, flow_rate_ml_min: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    vessel_radius = geometry["vessel"]["radius_m"]
    vessel_height = geometry["vessel"]["height_m"]
    basket_radius = geometry["basket"]["radius_m"]
    basket_z_min = geometry["basket"]["z_min_m"]
    basket_z_max = geometry["basket"]["z_max_m"]

    radius = np.linspace(0.0, vessel_radius, NR)
    axial = np.linspace(0.0, vessel_height, NZ)
    radial_mask = 1.0 / (1.0 + np.exp((radius - basket_radius) / (2.5 * (radius[1] - radius[0]))))
    axial_mask = logistic_window(axial, basket_z_min, basket_z_max, width=1.5e-3)

    flow_rate_m3_s = flow_rate_ml_min * 1.0e-6 / 60.0
    superficial_velocity = flow_rate_m3_s / (np.pi * basket_radius**2)
    inertial_drag = RHO * FORCHHEIMER * superficial_velocity / max(POROSITY, 1e-9)
    porous_drag = (MU / PERMEABILITY + inertial_drag) * POROSITY

    velocity_field = np.zeros((NZ, NR), dtype=float)
    dp_dz = np.zeros(NZ, dtype=float)

    for j, mask in enumerate(axial_mask):
        drag_profile = porous_drag * mask * radial_mask
        unit_profile = solve_radial_profile(radius, drag_profile)
        unit_flow = 2.0 * np.pi * np.trapezoid(unit_profile * radius, radius)
        if unit_flow <= 0.0:
            raise RuntimeError(f"Non-positive unit flow at axial index {j}.")
        gradient = flow_rate_m3_s / unit_flow
        velocity_field[j, :] = gradient * unit_profile
        dp_dz[j] = -gradient

    pressure = np.zeros(NZ, dtype=float)
    for j in range(NZ - 2, -1, -1):
        dz = axial[j + 1] - axial[j]
        pressure[j] = pressure[j + 1] - 0.5 * (dp_dz[j] + dp_dz[j + 1]) * dz

    dr = radius[1] - radius[0]
    wall_shear = MU * np.abs((velocity_field[:, -1] - velocity_field[:, -2]) / dr)
    basket_idx = min(max(np.searchsorted(radius, basket_radius), 1), NR - 2)
    basket_shear = MU * np.abs((velocity_field[:, basket_idx + 1] - velocity_field[:, basket_idx - 1]) / (2.0 * dr))

    return radius, axial, velocity_field, pressure, np.column_stack([wall_shear, basket_shear])


def write_csv(path: Path, fieldnames: List[str], rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def save_results(
    root: Path,
    geometry: Dict,
    fields: Dict[float, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
) -> List[Dict[str, float]]:
    velocity_rows: List[Dict[str, object]] = []
    shear_rows: List[Dict[str, float]] = []

    for flow_rate, (radius, axial, velocity_field, pressure, shear) in fields.items():
        basket_radius = geometry["basket"]["radius_m"]
        basket_z_min = geometry["basket"]["z_min_m"]
        basket_z_max = geometry["basket"]["z_max_m"]

        for j, z_val in enumerate(axial):
            for i, r_val in enumerate(radius):
                velocity_rows.append(
                    {
                        "flow_rate_ml_min": flow_rate,
                        "z_m": z_val,
                        "r_m": r_val,
                        "u_z_m_s": velocity_field[j, i],
                        "velocity_magnitude_m_s": abs(velocity_field[j, i]),
                        "pressure_Pa": pressure[j],
                        "region": "basket" if (r_val <= basket_radius and basket_z_min <= z_val <= basket_z_max) else "bulk",
                    }
                )
            shear_rows.append(
                {
                    "flow_rate_ml_min": flow_rate,
                    "z_m": z_val,
                    "wall_shear_Pa": float(shear[j, 0]),
                    "basket_surface_shear_Pa": float(shear[j, 1]),
                    "pressure_Pa": float(pressure[j]),
                }
            )

    write_csv(
        root / "results" / "velocity_field.csv",
        ["flow_rate_ml_min", "z_m", "r_m", "u_z_m_s", "velocity_magnitude_m_s", "pressure_Pa", "region"],
        velocity_rows,
    )
    write_csv(
        root / "results" / "shear_stress.csv",
        ["flow_rate_ml_min", "z_m", "wall_shear_Pa", "basket_surface_shear_Pa", "pressure_Pa"],
        shear_rows,
    )
    return shear_rows


def add_basket_outline(ax, geometry: Dict) -> None:
    basket_radius = geometry["basket"]["radius_m"] * 1e3
    z_min = geometry["basket"]["z_min_m"] * 1e3
    z_max = geometry["basket"]["z_max_m"] * 1e3
    ax.plot(
        [z_min, z_max, z_max, z_min, z_min],
        [basket_radius, basket_radius, 0.0, 0.0, basket_radius],
        color="white",
        linestyle="--",
        linewidth=1.0,
        alpha=0.9,
    )


def create_figures(
    root: Path,
    geometry: Dict,
    results: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    shear_rows: List[Dict[str, float]],
) -> None:
    radius, axial, velocity_field, pressure, _ = results
    z_mm, r_mm = np.meshgrid(axial * 1e3, radius * 1e3, indexing="ij")
    pressure_field = np.repeat(pressure[:, None], radius.size, axis=1)

    plt.style.use("seaborn-v0_8-whitegrid")

    fig1, ax1 = plt.subplots(figsize=(7.2, 4.8))
    c1 = ax1.contourf(z_mm, r_mm, velocity_field, levels=30, cmap="viridis")
    add_basket_outline(ax1, geometry)
    ax1.set_xlabel("Axial position (mm)")
    ax1.set_ylabel("Radius (mm)")
    ax1.set_title("Velocity magnitude field")
    cb1 = fig1.colorbar(c1, ax=ax1)
    cb1.set_label("Velocity (m/s)")
    fig1.tight_layout()
    fig1.savefig(root / "figures" / "velocity_field.png", dpi=300)
    fig1.savefig(root / "figures" / "velocity_field.svg")
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(7.2, 4.8))
    c2 = ax2.contourf(z_mm, r_mm, pressure_field, levels=30, cmap="cividis")
    add_basket_outline(ax2, geometry)
    ax2.set_xlabel("Axial position (mm)")
    ax2.set_ylabel("Radius (mm)")
    ax2.set_title("Pressure distribution")
    cb2 = fig2.colorbar(c2, ax=ax2)
    cb2.set_label("Pressure (Pa)")
    fig2.tight_layout()
    fig2.savefig(root / "figures" / "pressure_field.png", dpi=300)
    fig2.savefig(root / "figures" / "pressure_field.svg")
    plt.close(fig2)

    fig3, ax3 = plt.subplots(figsize=(7.2, 4.8))
    for flow_rate in FLOW_RATES_ML_MIN:
        subset = [row for row in shear_rows if abs(row["flow_rate_ml_min"] - flow_rate) < 1e-12]
        z_vals = np.array([row["z_m"] * 1e3 for row in subset])
        shear_vals = np.array([row["basket_surface_shear_Pa"] for row in subset])
        ax3.plot(z_vals, shear_vals, linewidth=2, label=f"{flow_rate:.1f} mL/min")
    ax3.axvspan(geometry["basket"]["z_min_m"] * 1e3, geometry["basket"]["z_max_m"] * 1e3, color="0.9", zorder=0)
    ax3.set_xlabel("Axial position (mm)")
    ax3.set_ylabel("Shear stress (Pa)")
    ax3.set_title("Shear stress along basket surface")
    ax3.legend(frameon=True)
    fig3.tight_layout()
    fig3.savefig(root / "figures" / "shear_stress_distribution.png", dpi=300)
    fig3.savefig(root / "figures" / "shear_stress_distribution.svg")
    plt.close(fig3)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    try:
        geometry = load_geometry(root)
        all_fields: Dict[float, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
        for flow_rate in FLOW_RATES_ML_MIN:
            all_fields[flow_rate] = compute_fields(geometry, flow_rate)
        shear_rows = save_results(root, geometry, all_fields)
        representative = all_fields[REPRESENTATIVE_FLOW]
        create_figures(root, geometry, representative, shear_rows)
        print("Simplified CFD results written to results/ and figures/.")
        return 0
    except Exception as exc:  # pragma: no cover - runtime guard
        print(f"Simulation failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
