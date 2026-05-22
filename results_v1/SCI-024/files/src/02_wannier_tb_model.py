"""
Topological Insulator Design Framework
Module 2: Wannier Function / Tight-Binding Model (Bi2Se3 4-band k·p + TB)

Based on:
  Liu et al., PRB 82, 045122 (2010)  — effective model for Bi2Se3
  Zhang et al., Nature Physics 5, 438 (2009) — original band inversion
"""

import numpy as np
from numpy import linalg as LA
import json
import os

# ---------------------------------------------------------------------------
# Physical constants and model parameters for Bi2Se3
# ---------------------------------------------------------------------------

# k·p model parameters (Liu et al. PRB 2010)
# Units: eV for energies, Å for lengths, eV·Å for A parameters
BI2SE3_PARAMS = {
    "M0":  0.28,    # eV  — Dirac mass at Gamma
    "M1": -10.0,    # eV·Å² — kz² coefficient in M(k)
    "M2": -56.6,    # eV·Å² — k_perp² coefficient in M(k)
    "A1":  2.2,     # eV·Å — velocity in kz direction
    "A2":  4.1,     # eV·Å — velocity in k_perp direction
    "C":   -0.0068, # eV  — overall shift
    "D1":  1.3,     # eV·Å² — particle-hole asymmetry kz
    "D2":  19.6,    # eV·Å² — particle-hole asymmetry k_perp
    "soc_scale": 1.0,  # dimensionless SOC scaling factor
}

# ---------------------------------------------------------------------------
# Gamma (Dirac) matrices for the 4-band model
# Basis: |p1+_z↑⟩, |p2-_z↑⟩, |p1+_z↓⟩, |p2-_z↓⟩
# ---------------------------------------------------------------------------

def gamma_matrices():
    """
    4x4 Gamma matrices for the Bi2Se3 effective Hamiltonian.
    Γ1 = σx⊗σx,  Γ2 = σx⊗σy,  Γ3 = σx⊗σz,  Γ4 = σy⊗I,  Γ5 = σz⊗I
    """
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    I2 = np.eye(2, dtype=complex)

    G = {}
    G[1] = np.kron(sx, sx)
    G[2] = np.kron(sx, sy)
    G[3] = np.kron(sx, sz)
    G[4] = np.kron(sy, I2)
    G[5] = np.kron(sz, I2)
    return G


# ---------------------------------------------------------------------------
# 4-band k·p Hamiltonian
# ---------------------------------------------------------------------------

class Bi2Se3KPModel:
    """
    Low-energy 4-band k·p model for Bi2Se3 (and family).
    H(k) = ε(k)·I₄ + M(k)·Γ₅ + A₁kz·Γ₄ + A₂(kx·Γ₁ + ky·Γ₂)
    """

    def __init__(self, params: dict = None):
        self.p = {**BI2SE3_PARAMS, **(params or {})}
        self.G = gamma_matrices()

    def hamiltonian(self, kx: float, ky: float, kz: float) -> np.ndarray:
        p = self.p
        lam = p["soc_scale"]

        k2_perp = kx**2 + ky**2
        k2_z = kz**2

        eps = p["C"] + p["D1"] * k2_z + p["D2"] * k2_perp
        M = p["M0"] - lam * (p["M1"] * k2_z + p["M2"] * k2_perp)

        H = eps * np.eye(4, dtype=complex)
        H += M * self.G[5]
        H += lam * p["A1"] * kz * self.G[4]
        H += lam * p["A2"] * (kx * self.G[1] + ky * self.G[2])

        return H

    def band_energies(self, k_points: np.ndarray) -> np.ndarray:
        """
        Compute eigenvalues along a k-path.
        k_points: (N, 3) array of [kx, ky, kz] in Å⁻¹
        Returns (N, 4) array of energies in eV.
        """
        nk = len(k_points)
        energies = np.zeros((nk, 4))

        for i, (kx, ky, kz) in enumerate(k_points):
            H = self.hamiltonian(kx, ky, kz)
            evals = np.sort(LA.eigvalsh(H))
            energies[i] = evals

        return energies

    def band_gap(self, n_kpts: int = 50) -> float:
        """Estimate band gap over a dense k-grid."""
        kx = np.linspace(-0.3, 0.3, n_kpts)
        ky = np.linspace(-0.3, 0.3, n_kpts)
        kz = np.linspace(-0.3, 0.3, n_kpts)

        min_gap = np.inf
        for kxi in kx:
            for kyi in ky:
                for kzi in kz:
                    H = self.hamiltonian(kxi, kyi, kzi)
                    evals = np.sort(LA.eigvalsh(H))
                    gap = evals[2] - evals[1]  # gap between bands 2 and 3
                    if gap < min_gap:
                        min_gap = gap
        return float(min_gap)


# ---------------------------------------------------------------------------
# Wannier-based tight-binding on a 3D lattice (simplified)
# ---------------------------------------------------------------------------

class WannierTBModel:
    """
    Wannier-function-based tight-binding model.
    Mimics the output of Wannier90 for Bi2Se3.

    The 3D model uses:
      - 5 Wannier functions per unit cell (pz-like on Bi and Se sites)
      - Nearest/next-nearest neighbor hopping
      - On-site SOC terms

    Lattice: rhombohedral, approximated as hexagonal supercell
    """

    def __init__(self, material: str = "Bi2Se3"):
        self.material = material
        self._set_hopping_parameters(material)

    def _set_hopping_parameters(self, material: str):
        """Set effective hopping parameters for Bi2Se3 and analogs."""
        hop_db = {
            "Bi2Se3": {
                "t1": 0.45,   # nearest-neighbor QL coupling (eV)
                "t2": 0.15,   # next-QL coupling
                "t3": 0.05,   # third-neighbor
                "soc": 0.30,  # on-site SOC lambda (eV)
                "delta": 0.28, # Dirac mass (band inversion gap) (eV)
                "mu": 0.0,    # chemical potential
            },
            "Bi2Te3": {
                "t1": 0.42, "t2": 0.13, "t3": 0.04,
                "soc": 0.38, "delta": 0.15, "mu": 0.0,
            },
            "Sb2Te3": {
                "t1": 0.40, "t2": 0.12, "t3": 0.04,
                "soc": 0.22, "delta": 0.21, "mu": 0.0,
            },
            "TlBiSe2": {
                "t1": 0.48, "t2": 0.14, "t3": 0.05,
                "soc": 0.35, "delta": 0.35, "mu": 0.0,
            },
            "MnBi2Te4": {
                "t1": 0.44, "t2": 0.14, "t3": 0.05,
                "soc": 0.36, "delta": 0.20, "mu": 0.0,
                "mag": 0.05,  # magnetic exchange
            },
        }
        params = hop_db.get(material, hop_db["Bi2Se3"])
        self.p = params

    def hamiltonian_2d(self, kx: float, ky: float) -> np.ndarray:
        """
        Effective 2D Hamiltonian from integrating out kz.
        4-orbital model: |s,↑⟩, |px+ipy,↑⟩, |s,↓⟩, |px-ipy,↓⟩
        """
        p = self.p
        t1, t2, soc, delta = p["t1"], p["t2"], p["soc"], p["delta"]
        mu = p.get("mu", 0.0)

        # k-dependent hopping on triangular lattice (a=1)
        g = (np.exp(1j * kx)
             + np.exp(1j * (-0.5 * kx + np.sqrt(3) / 2 * ky))
             + np.exp(1j * (-0.5 * kx - np.sqrt(3) / 2 * ky)))

        gstar = np.conj(g)

        # Next-neighbor contribution
        g2 = (np.exp(1j * (kx + ky * np.sqrt(3) / 2))
              + np.exp(-1j * ky * np.sqrt(3))
              + np.exp(1j * (-kx + ky * np.sqrt(3) / 2)))

        H = np.zeros((4, 4), dtype=complex)
        # On-site: mass + SOC
        H[0, 0] = delta / 2 - mu
        H[1, 1] = -delta / 2 + soc - mu
        H[2, 2] = delta / 2 - mu
        H[3, 3] = -delta / 2 + soc - mu

        # Off-diagonal hoppings
        H[0, 1] = t1 * g + t2 * g2
        H[1, 0] = np.conj(H[0, 1])
        H[2, 3] = -(t1 * g + t2 * g2)
        H[3, 2] = np.conj(H[2, 3])

        # Spin-flip SOC
        H[0, 3] = 1j * soc * 0.1 * (kx - 1j * ky)
        H[3, 0] = np.conj(H[0, 3])
        H[1, 2] = 1j * soc * 0.1 * (kx + 1j * ky)
        H[2, 1] = np.conj(H[1, 2])

        return H

    def band_structure_path(self, path_name: str = "hex") -> tuple:
        """
        Return k-path and band energies for standard high-symmetry path.
        """
        # Hexagonal BZ high-symmetry points (units: 2pi/a, a=1)
        if path_name == "hex":
            sym_pts = {
                "Γ": np.array([0.0, 0.0]),
                "M": np.array([np.pi, np.pi / np.sqrt(3)]),
                "K": np.array([4 * np.pi / 3, 0.0]),
            }
            segments = [("Γ", "M"), ("M", "K"), ("K", "Γ")]
        else:
            sym_pts = {
                "Γ": np.array([0.0, 0.0]),
                "X": np.array([np.pi, 0.0]),
                "M": np.array([np.pi, np.pi]),
            }
            segments = [("Γ", "X"), ("X", "M"), ("M", "Γ")]

        npts_seg = 60
        k_path = []
        k_dist = []
        ticks = []
        tick_labels = []
        dist = 0.0

        for start_label, end_label in segments:
            start = sym_pts[start_label]
            end = sym_pts[end_label]
            ticks.append(dist)
            tick_labels.append(start_label)

            for i in range(npts_seg):
                t = i / npts_seg
                k = start + t * (end - start)
                k_path.append(k)
                k_dist.append(dist + LA.norm(end - start) * t)

            dist += LA.norm(end - start)

        ticks.append(dist)
        tick_labels.append(segments[-1][1])
        k_path = np.array(k_path)

        energies = np.zeros((len(k_path), 4))
        for i, k in enumerate(k_path):
            H = self.hamiltonian_2d(k[0], k[1])
            energies[i] = np.sort(LA.eigvalsh(H))

        return np.array(k_dist), energies, ticks, tick_labels


def run_wannier_tb():
    """Run tight-binding model calculations."""
    results = {}
    os.makedirs("results", exist_ok=True)

    print("=" * 60)
    print("WANNIER TIGHT-BINDING MODEL — BAND STRUCTURE")
    print("=" * 60)

    for mat in ["Bi2Se3", "Bi2Te3", "Sb2Te3", "TlBiSe2"]:
        model = WannierTBModel(mat)
        k_dist, energies, ticks, tick_labels = model.band_structure_path("hex")
        gap = energies[:, 2].min() - energies[:, 1].max()

        print(f"  {mat:12s}: Band gap = {gap:.3f} eV")
        results[mat] = {
            "k_dist": k_dist.tolist(),
            "energies_eV": energies.tolist(),
            "ticks": ticks,
            "tick_labels": tick_labels,
            "band_gap_eV": float(gap),
            "params": model.p,
        }

    with open("results/wannier_tb_bands.json", "w") as f:
        json.dump(results, f, indent=2)

    print("\nSaved: results/wannier_tb_bands.json")
    return results


if __name__ == "__main__":
    run_wannier_tb()
