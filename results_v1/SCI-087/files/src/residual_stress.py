"""
Module 3: Residual Stress & Warpage Prediction
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict


@dataclass
class MechanicalProperties:
    """Temperature-dependent mechanical properties for PA66-GF30."""
    E_fiber_dir: float = 11000e6
    E_transverse: float = 5500e6
    nu: float = 0.38
    CTE_fiber: float = 2.0e-5
    CTE_transverse: float = 7.0e-5
    crystallization_shrinkage: float = 0.012


def compute_thermal_stress_profile(
    temperature_profile: np.ndarray,
    crystallinity_profile: np.ndarray,
    z_positions: np.ndarray,
    mold_temp: float,
    mech: MechanicalProperties,
    reference_temp: float = 220.0
) -> Dict:
    nz = len(z_positions)
    E_profile = np.zeros(nz)
    for i in range(nz):
        T = temperature_profile[i]
        E_ref = mech.E_fiber_dir
        if T > 200:
            E_profile[i] = E_ref * 0.3
        elif T > 150:
            E_profile[i] = E_ref * (0.3 + 0.7 * (200 - T) / 50)
        else:
            E_profile[i] = E_ref * (1.0 - 0.002 * T)
    alpha_avg = (mech.CTE_fiber + mech.CTE_transverse) / 2
    thermal_strain = alpha_avg * (reference_temp - temperature_profile)
    cryst_strain = crystallinity_profile * mech.crystallization_shrinkage / 0.35
    free_strain = thermal_strain + cryst_strain
    eps_avg = np.trapz(E_profile * free_strain, z_positions) / np.trapz(E_profile, z_positions)
    z_mid = (z_positions[-1] + z_positions[0]) / 2
    z_centered = z_positions - z_mid
    M = np.trapz(E_profile * (free_strain - eps_avg) * z_centered, z_positions)
    I = np.trapz(E_profile * z_centered ** 2, z_positions)
    kappa = M / I if abs(I) > 1e-20 else 0.0
    sigma = E_profile / (1 - mech.nu ** 2) * (free_strain - eps_avg - kappa * z_centered)
    return {
        'z_mm': z_positions.tolist(),
        'stress_MPa': (sigma / 1e6).tolist(),
        'thermal_strain': thermal_strain.tolist(),
        'crystallization_strain': cryst_strain.tolist(),
        'modulus_GPa': (E_profile / 1e9).tolist(),
        'curvature_1_m': float(kappa),
        'max_tensile_MPa': float(np.max(sigma) / 1e6),
        'max_compressive_MPa': float(np.min(sigma) / 1e6),
        'avg_stress_MPa': float(np.mean(sigma) / 1e6),
    }


def predict_warpage(part_length: float, part_width: float,
                     curvature_x: float, curvature_y: float) -> Dict:
    warp_x = curvature_x * part_length ** 2 / 8 * 1000
    warp_y = curvature_y * part_width ** 2 / 8 * 1000
    warp_total = np.sqrt(warp_x ** 2 + warp_y ** 2)
    nx, ny = 20, 10
    x = np.linspace(0, part_length, nx)
    y = np.linspace(0, part_width, ny)
    X, Y = np.meshgrid(x, y)
    x_mid, y_mid = part_length / 2, part_width / 2
    W = curvature_x / 2 * (X - x_mid) ** 2 + curvature_y / 2 * (Y - y_mid) ** 2
    return {
        'warpage_x_mm': float(warp_x),
        'warpage_y_mm': float(warp_y),
        'warpage_total_mm': float(warp_total),
        'deformation_field_mm': (W * 1000).tolist(),
        'x_coords_mm': (x * 1000).tolist(),
        'y_coords_mm': (y * 1000).tolist(),
        'max_deformation_mm': float(np.max(np.abs(W)) * 1000),
    }


def run_stress_warpage_analysis(
    temperature_profile_C: np.ndarray,
    crystallinity_profile: np.ndarray,
    part_length: float = 0.200,
    part_width: float = 0.100,
    part_thickness: float = 0.003,
) -> Dict:
    nz = len(temperature_profile_C)
    z_positions = np.linspace(0, part_thickness * 1000, nz)
    mech = MechanicalProperties()
    stress_result = compute_thermal_stress_profile(
        temperature_profile_C, crystallinity_profile, z_positions,
        mold_temp=80.0, mech=mech
    )
    kappa_x = stress_result['curvature_1_m']
    kappa_y = kappa_x * (mech.CTE_transverse / mech.CTE_fiber)
    warpage_result = predict_warpage(part_length, part_width, kappa_x, kappa_y)
    return {'stress': stress_result, 'warpage': warpage_result}
