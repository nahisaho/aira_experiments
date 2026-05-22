#!/usr/bin/env python3
"""Brain organoid bioreactor oxygen and nutrient transport modeling.

This script solves steady-state and time-dependent reaction-diffusion equations
for oxygen and glucose inside a spherical organoid. Michaelis-Menten kinetics
are used for oxygen and glucose consumption, and lactate is modeled as a
byproduct of glucose utilization.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

np.random.seed(42)

BASE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = BASE_DIR / "figures"

D_O2 = 2.5e-9
VMAX_O2 = 5.0e-3
KM_O2 = 4.6e-3
C_EXT_O2 = 0.20
O2_CRITICAL = 0.01

D_GLUCOSE = 6.7e-10
VMAX_GLUCOSE = 1.2e-2
KM_GLUCOSE = 0.5
C_EXT_GLUCOSE = 5.0

D_LACTATE = 1.0e-9
LACTATE_YIELD = 1.8
RADIUS_SET = np.array([0.5e-3, 1.0e-3, 1.5e-3, 2.0e-3])

plt.style.use("seaborn-v0_8-whitegrid")


@dataclass(frozen=True)
class SpeciesParams:
    diffusion: float
    vmax: float
    km: float
    boundary: float


def michaelis_menten(concentration: np.ndarray, vmax: float, km: float) -> np.ndarray:
    """Return Michaelis-Menten consumption rate."""
    conc = np.clip(concentration, 0.0, None)
    return vmax * conc / (km + conc + 1e-12)


def build_linear_system(radius_m: float, n_points: int, diffusion: float, sink_coeff: np.ndarray, boundary: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a linearized steady-state system for spherical diffusion-reaction."""
    radial_grid = np.linspace(0.0, radius_m, n_points)
    dr = radial_grid[1] - radial_grid[0]
    n_unknowns = n_points - 1
    matrix = np.zeros((n_unknowns, n_unknowns), dtype=float)
    rhs = np.zeros(n_unknowns, dtype=float)

    matrix[0, 0] = -6.0 * diffusion / dr**2 - sink_coeff[0]
    matrix[0, 1] = 6.0 * diffusion / dr**2

    for i in range(1, n_unknowns):
        radius = radial_grid[i]
        lower = diffusion * (1.0 / dr**2 - 1.0 / (radius * dr))
        diag = -2.0 * diffusion / dr**2 - sink_coeff[i]
        upper = diffusion * (1.0 / dr**2 + 1.0 / (radius * dr))
        matrix[i, i - 1] = lower
        matrix[i, i] = diag
        if i < n_unknowns - 1:
            matrix[i, i + 1] = upper
        else:
            rhs[i] -= upper * boundary

    return radial_grid, matrix, rhs


def solve_consumption_profile(
    radius_m: float,
    params: SpeciesParams,
    n_points: int = 120,
    max_iter: int = 80,
    tol: float = 1e-8,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve steady-state profile by Picard linearization."""
    profile = np.full(n_points, params.boundary, dtype=float)
    radial_grid = np.linspace(0.0, radius_m, n_points)

    for _ in range(max_iter):
        sink_coeff = params.vmax / (params.km + np.clip(profile[:-1], 0.0, None) + 1e-12)
        radial_grid, matrix, rhs = build_linear_system(radius_m, n_points, params.diffusion, sink_coeff, params.boundary)
        updated = np.linalg.solve(matrix, rhs)
        new_profile = np.clip(np.concatenate([updated, [params.boundary]]), 0.0, None)
        if np.max(np.abs(new_profile - profile)) < tol:
            return radial_grid, new_profile
        profile = 0.65 * profile + 0.35 * new_profile

    return radial_grid, profile


def solve_lactate_profile(radial_grid: np.ndarray, glucose_profile: np.ndarray, boundary_value: float = 0.0) -> np.ndarray:
    """Solve steady-state lactate diffusion with glucose-derived production."""
    n_points = radial_grid.size
    radius_m = radial_grid[-1]
    source = LACTATE_YIELD * michaelis_menten(glucose_profile[:-1], VMAX_GLUCOSE, KM_GLUCOSE)
    _, matrix, rhs = build_linear_system(radius_m, n_points, D_LACTATE, np.zeros(n_points - 1), boundary_value)
    rhs -= source
    interior = np.linalg.solve(matrix, rhs)
    return np.concatenate([interior, [boundary_value]])


def radial_laplacian(profile: np.ndarray, dr: float) -> np.ndarray:
    """Radial Laplacian for transient integration."""
    lap = np.zeros_like(profile)
    lap[0] = 6.0 * (profile[1] - profile[0]) / dr**2
    idx = np.arange(1, profile.size - 1)
    radius = idx * dr
    lap[idx] = (
        (profile[idx + 1] - 2.0 * profile[idx] + profile[idx - 1]) / dr**2
        + (profile[idx + 1] - profile[idx - 1]) / (radius * dr)
    )
    return lap


def time_dependent_profile(radius_m: float, params: SpeciesParams, t_end_s: float = 4.0 * 3600.0, n_points: int = 90) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Integrate the transient diffusion-consumption problem."""
    radial_grid = np.linspace(0.0, radius_m, n_points)
    dr = radial_grid[1] - radial_grid[0]
    initial_state = np.full(n_points - 1, params.boundary, dtype=float)
    sample_times = np.linspace(0.0, t_end_s, 6)

    def rhs(_: float, interior: np.ndarray) -> np.ndarray:
        full = np.concatenate([np.clip(interior, 0.0, None), [params.boundary]])
        diffusion_term = params.diffusion * radial_laplacian(full, dr)[:-1]
        reaction_term = michaelis_menten(full[:-1], params.vmax, params.km)
        return diffusion_term - reaction_term

    solution = solve_ivp(
        rhs,
        (0.0, t_end_s),
        initial_state,
        method="BDF",
        t_eval=sample_times,
        atol=1e-8,
        rtol=1e-6,
    )
    if not solution.success:
        raise RuntimeError("Transient diffusion solver failed")
    profiles = np.vstack([np.concatenate([solution.y[:, i], [params.boundary]]) for i in range(solution.y.shape[1])])
    return radial_grid, solution.t, profiles


def center_concentration(radius_m: float, params: SpeciesParams) -> float:
    """Return center concentration at steady state."""
    _, profile = solve_consumption_profile(radius_m, params, n_points=90)
    return float(profile[0])


def find_critical_radius(params: SpeciesParams, threshold: float = O2_CRITICAL, lower: float = 0.3e-3, upper: float = 3.0e-3, iterations: int = 12) -> float:
    """Estimate the radius where center oxygen crosses the critical threshold."""
    if center_concentration(upper, params) >= threshold:
        return upper
    lo, hi = lower, upper
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        if center_concentration(mid, params) >= threshold:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def save_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    """Write rows to CSV."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def make_profile_figure(profiles: dict[float, tuple[np.ndarray, np.ndarray]], ylabel: str, filename: str, cmap_name: str) -> None:
    """Save radial concentration profile figure."""
    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.get_cmap(cmap_name)
    for idx, (radius_m, (radial_grid, profile)) in enumerate(profiles.items()):
        ax.plot(radial_grid * 1e3, profile, color=cmap(idx / max(1, len(profiles) - 1)), linewidth=2.4, label=f"R = {radius_m * 1e3:.1f} mm")
    ax.set_xlabel("Radius from center (mm)")
    ax.set_ylabel(ylabel)
    ax.set_title(ylabel.split(" (")[0] + " radial profiles")
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / filename, dpi=300)
    plt.close(fig)


def make_critical_radius_figure(critical_rows: list[dict]) -> None:
    """Plot critical radius sensitivity to metabolic demand and external oxygen."""
    fig, ax = plt.subplots(figsize=(7, 5))
    ext_levels = sorted({row["external_o2_mol_m3"] for row in critical_rows})
    cmap = plt.get_cmap("viridis")
    for idx, ext in enumerate(ext_levels):
        subset = [row for row in critical_rows if row["external_o2_mol_m3"] == ext]
        ax.plot(
            [row["vmax_scale"] for row in subset],
            [row["critical_radius_mm"] for row in subset],
            marker="o",
            linewidth=2.2,
            color=cmap(idx / max(1, len(ext_levels) - 1)),
            label=f"External O2 = {ext:.2f} mol m$^{{-3}}$",
        )
    ax.set_xlabel("Oxygen Vmax scaling factor")
    ax.set_ylabel("Critical radius (mm)")
    ax.set_title("Critical radius sensitivity analysis")
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "critical_radius_analysis.png", dpi=300)
    plt.close(fig)


def make_time_evolution_figure(radial_grid: np.ndarray, times_s: np.ndarray, profiles: np.ndarray) -> None:
    """Plot transient oxygen concentration profiles."""
    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = plt.get_cmap("cividis")
    for idx, (time_s, profile) in enumerate(zip(times_s, profiles)):
        ax.plot(radial_grid * 1e3, profile, color=cmap(idx / max(1, len(times_s) - 1)), linewidth=2.2, label=f"t = {time_s / 3600.0:.1f} h")
    ax.axhline(O2_CRITICAL, linestyle="--", color="black", linewidth=1.4, label="Critical O2 threshold")
    ax.set_xlabel("Radius from center (mm)")
    ax.set_ylabel("Oxygen concentration (mol m$^{-3}$)")
    ax.set_title("Transient oxygen diffusion-consumption")
    ax.legend(frameon=True, ncol=2)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "oxygen_time_evolution.png", dpi=300)
    plt.close(fig)


def main() -> None:
    """Generate transport-model CSV files and figures."""
    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)

    oxygen_params = SpeciesParams(D_O2, VMAX_O2, KM_O2, C_EXT_O2)
    glucose_params = SpeciesParams(D_GLUCOSE, VMAX_GLUCOSE, KM_GLUCOSE, C_EXT_GLUCOSE)

    oxygen_profiles: dict[float, tuple[np.ndarray, np.ndarray]] = {}
    glucose_profiles: dict[float, tuple[np.ndarray, np.ndarray]] = {}
    oxygen_rows: list[dict] = []
    nutrient_rows: list[dict] = []

    for radius_m in RADIUS_SET:
        radial_o2, o2_profile = solve_consumption_profile(radius_m, oxygen_params)
        radial_glc, glucose_profile = solve_consumption_profile(radius_m, glucose_params)
        lactate_profile = solve_lactate_profile(radial_glc, glucose_profile)

        oxygen_profiles[radius_m] = (radial_o2, o2_profile)
        glucose_profiles[radius_m] = (radial_glc, glucose_profile)

        for radial_position, concentration in zip(radial_o2, o2_profile):
            oxygen_rows.append(
                {
                    "species": "oxygen",
                    "organoid_radius_mm": f"{radius_m * 1e3:.3f}",
                    "radial_position_mm": f"{radial_position * 1e3:.6f}",
                    "radial_fraction": f"{radial_position / radius_m:.6f}",
                    "concentration_mol_m3": f"{concentration:.8f}",
                }
            )

        for radial_position, glucose_c, lactate_c in zip(radial_glc, glucose_profile, lactate_profile):
            nutrient_rows.extend(
                [
                    {
                        "species": "glucose",
                        "organoid_radius_mm": f"{radius_m * 1e3:.3f}",
                        "radial_position_mm": f"{radial_position * 1e3:.6f}",
                        "radial_fraction": f"{radial_position / radius_m:.6f}",
                        "concentration_mol_m3": f"{glucose_c:.8f}",
                    },
                    {
                        "species": "lactate",
                        "organoid_radius_mm": f"{radius_m * 1e3:.3f}",
                        "radial_position_mm": f"{radial_position * 1e3:.6f}",
                        "radial_fraction": f"{radial_position / radius_m:.6f}",
                        "concentration_mol_m3": f"{lactate_c:.8f}",
                    },
                ]
            )

    critical_radius_mm = find_critical_radius(oxygen_params) * 1e3
    vmax_scales = np.linspace(0.6, 1.4, 5)
    external_o2_levels = [0.10, 0.15, 0.20, 0.25, 0.30]
    critical_rows = []
    for external_o2 in external_o2_levels:
        for scale in vmax_scales:
            trial_params = SpeciesParams(D_O2, VMAX_O2 * scale, KM_O2, external_o2)
            critical_rows.append(
                {
                    "external_o2_mol_m3": external_o2,
                    "vmax_scale": scale,
                    "critical_radius_mm": find_critical_radius(trial_params) * 1e3,
                }
            )

    radial_grid_t, times_s, profiles_t = time_dependent_profile(1.5e-3, oxygen_params)

    save_csv(
        RESULTS_DIR / "oxygen_profiles.csv",
        oxygen_rows,
        ["species", "organoid_radius_mm", "radial_position_mm", "radial_fraction", "concentration_mol_m3"],
    )
    save_csv(
        RESULTS_DIR / "nutrient_profiles.csv",
        nutrient_rows,
        ["species", "organoid_radius_mm", "radial_position_mm", "radial_fraction", "concentration_mol_m3"],
    )
    save_csv(
        RESULTS_DIR / "critical_radius_scan.csv",
        [
            {
                "baseline_critical_radius_mm": f"{critical_radius_mm:.6f}",
                "external_o2_mol_m3": f"{row['external_o2_mol_m3']:.3f}",
                "vmax_scale": f"{row['vmax_scale']:.3f}",
                "critical_radius_mm": f"{row['critical_radius_mm']:.6f}",
            }
            for row in critical_rows
        ],
        ["baseline_critical_radius_mm", "external_o2_mol_m3", "vmax_scale", "critical_radius_mm"],
    )

    make_profile_figure(oxygen_profiles, "Oxygen concentration (mol m$^{-3}$)", "oxygen_radial_profile.png", "viridis")
    make_profile_figure(glucose_profiles, "Glucose concentration (mol m$^{-3}$)", "glucose_radial_profile.png", "cividis")
    make_critical_radius_figure(critical_rows)
    make_time_evolution_figure(radial_grid_t, times_s, profiles_t)

    print(f"Baseline critical radius: {critical_radius_mm:.3f} mm")
    print("Saved oxygen and nutrient transport outputs.")


if __name__ == "__main__":
    main()
