"""
Topological Insulator Design Framework
Module 4: Surface State Dirac Dispersion — Slab Calculation

Computes surface spectral function A(k,ω) for a finite slab,
revealing the topological surface states (TSS) as Dirac cones.

Method: Transfer matrix / finite-slab diagonalization
"""

import numpy as np
from numpy import linalg as LA
import json
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
from _model_utils import build_tb_slab_2d


def compute_slab_bands(n_layers: int = 30, t: float = 0.45,
                       soc: float = 0.30, delta: float = 0.28,
                       n_kpts: int = 100) -> dict:
    """
    Compute band structure of a finite slab along Γ-K direction.
    Returns energies and identifies surface vs bulk bands.
    """
    H_slab = build_tb_slab_2d(n_layers=n_layers, t=t, lam=soc, delta=delta)
    norb = 4 * n_layers

    k_arr = np.linspace(-0.4, 0.4, n_kpts)  # Å⁻¹
    energies = np.zeros((n_kpts, norb))
    # Surface weight: project onto surface layers (first + last 2 layers)
    surf_weight = np.zeros((n_kpts, norb))
    surf_orbs = list(range(8)) + list(range(norb - 8, norb))

    for i, kx in enumerate(k_arr):
        H = H_slab(kx, 0.0)
        evals, evecs = LA.eigh(H)
        energies[i] = evals
        for n in range(norb):
            surf_weight[i, n] = np.sum(np.abs(evecs[surf_orbs, n])**2)

    # Identify Dirac cone: surface bands near zero energy around Gamma
    mid = n_kpts // 2
    mid_bands = energies[mid]
    # Surface states are near E=0 with high surface weight
    surf_mask = (np.abs(mid_bands) < 0.15) & (surf_weight[mid] > 0.3)
    dirac_band_indices = np.where(surf_mask)[0]

    # Estimate Dirac velocity: linear fit near Gamma
    dirac_velocity = None
    if len(dirac_band_indices) >= 2:
        # Use uppermost surface state below zero
        above_zero = [b for b in dirac_band_indices if mid_bands[b] > 0]
        if above_zero:
            bn = above_zero[0]
            # Linear fit over small k range
            k_fit = k_arr[mid - 10:mid + 10]
            e_fit = energies[mid - 10:mid + 10, bn]
            coeffs = np.polyfit(k_fit, e_fit, 1)
            # Dirac velocity in eV·Å
            dirac_velocity = float(abs(coeffs[0]))

    return {
        "k_arr": k_arr.tolist(),
        "energies_eV": energies.tolist(),
        "surf_weight": surf_weight.tolist(),
        "n_layers": n_layers,
        "n_bands": norb,
        "dirac_band_indices": dirac_band_indices.tolist(),
        "dirac_velocity_eVA": dirac_velocity,
        "params": {"t": t, "soc": soc, "delta": delta},
    }


def compute_spectral_function(H_slab_func, k_arr, omega_arr,
                              eta: float = 0.02) -> np.ndarray:
    """
    Compute surface spectral function A_surf(k, ω) via Green's function.
    A(k,ω) = -1/π · Im Tr[G_surf(k, ω+iη)]
    G(k,ω) = [(ω+iη)I - H(k)]⁻¹
    Only surface projected (first 2 layers).
    """
    nk = len(k_arr)
    nw = len(omega_arr)
    norb = None
    A = np.zeros((nk, nw))

    n_surf_orbs = 8  # first 2 layers × 4 orbitals

    for ik, kx in enumerate(k_arr):
        H = H_slab_func(kx, 0.0)
        if norb is None:
            norb = H.shape[0]
        I = np.eye(norb, dtype=complex)

        for iw, omega in enumerate(omega_arr):
            G = LA.inv((omega + 1j * eta) * I - H)
            # Surface projected spectral weight
            A[ik, iw] = -np.imag(np.trace(G[:n_surf_orbs, :n_surf_orbs])) / np.pi

    return A


def run_slab_calculations():
    """Run slab band structure and spectral function calculations."""
    os.makedirs("results", exist_ok=True)

    print("=" * 60)
    print("SLAB SURFACE STATE CALCULATION")
    print("=" * 60)

    # 1. Slab band structure for different materials/parameters
    materials = {
        "Bi2Se3": {"t": 0.45, "soc": 0.30, "delta": 0.28},
        "Bi2Te3": {"t": 0.42, "soc": 0.38, "delta": 0.15},
        "Sb2Te3": {"t": 0.40, "soc": 0.22, "delta": 0.21},
    }

    slab_results = {}
    for mat, params in materials.items():
        print(f"  Computing slab bands for {mat}...")
        res = compute_slab_bands(
            n_layers=24, n_kpts=80, **params
        )
        slab_results[mat] = res

        vD = res["dirac_velocity_eVA"]
        vD_str = f"{vD:.3f} eV·Å" if vD else "N/A"
        n_ss = len(res["dirac_band_indices"])
        print(f"    Surface bands near EF: {n_ss}  |  vDirac ≈ {vD_str}")

    # 2. Spectral function for Bi2Se3 (higher resolution)
    print("\n  Computing spectral function A(k,ω) for Bi2Se3...")
    H_surf = build_tb_slab_2d(n_layers=20, t=0.45, lam=0.30, delta=0.28)
    k_sf = np.linspace(-0.35, 0.35, 60)
    omega_sf = np.linspace(-0.8, 0.8, 80)
    A_surf = compute_spectral_function(H_surf, k_sf, omega_sf, eta=0.025)

    spectral_data = {
        "k_arr": k_sf.tolist(),
        "omega_arr": omega_sf.tolist(),
        "A_surface": A_surf.tolist(),
        "material": "Bi2Se3",
    }

    with open("results/slab_bands.json", "w") as f:
        json.dump(slab_results, f, indent=2)
    with open("results/spectral_function.json", "w") as f:
        json.dump(spectral_data, f, indent=2)

    print("Saved: results/slab_bands.json")
    print("Saved: results/spectral_function.json")
    return slab_results, spectral_data


if __name__ == "__main__":
    run_slab_calculations()
