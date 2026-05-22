"""
Shared model utilities used across all pipeline modules.
Provides Bi2Se3 k·p Hamiltonian factory functions.
"""
import numpy as np
from numpy import linalg as LA


def gamma_matrices():
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    I2 = np.eye(2, dtype=complex)
    return {
        1: np.kron(sx, sx),
        2: np.kron(sx, sy),
        3: np.kron(sx, sz),
        4: np.kron(sy, I2),
        5: np.kron(sz, I2),
    }


BI2SE3_PARAMS = {
    "M0": 0.28, "M1": -10.0, "M2": -56.6,
    "A1": 2.2, "A2": 4.1,
    "C": -0.0068, "D1": 1.3, "D2": 19.6,
}


def build_bi2se3_kp(soc_scale: float = 1.0):
    """
    Return a callable H(kx, ky) for the Bi2Se3 4-band k·p model
    (kz=0 plane, 2D effective model).
    """
    G = gamma_matrices()
    p = BI2SE3_PARAMS

    def H_func(kx: float, ky: float) -> np.ndarray:
        lam = soc_scale
        k2 = kx**2 + ky**2

        eps = p["C"] + p["D2"] * k2
        M = p["M0"] - lam * p["M2"] * k2

        H = eps * np.eye(4, dtype=complex)
        H += M * G[5]
        H += lam * p["A2"] * (kx * G[1] + ky * G[2])
        return H

    return H_func


def build_tb_slab_2d(n_layers: int = 20, t: float = 0.45, lam: float = 0.30,
                     delta: float = 0.28):
    """
    Build 2D slab Hamiltonian for surface state calculation.
    Returns callable H(kx, ky) for (4*n_layers x 4*n_layers) slab.
    """
    norb = 4 * n_layers  # 4 orbitals per layer

    def H_slab(kx: float, ky: float) -> np.ndarray:
        H = np.zeros((norb, norb), dtype=complex)

        k2 = kx**2 + ky**2
        vk = lam * np.sqrt(k2) if k2 > 0 else 0.0
        phi = np.arctan2(ky, kx) if k2 > 0 else 0.0

        for il in range(n_layers):
            base = il * 4
            # On-site block: 4x4 k·p at this layer
            eps = 0.0
            M = delta * (1 - k2 * 0.5)  # simplified dispersion
            soc = lam

            # 4x4 block
            H[base, base] = eps + M
            H[base+1, base+1] = eps - M + soc
            H[base+2, base+2] = eps + M
            H[base+3, base+3] = eps - M + soc

            # In-plane velocity terms (spin-momentum locking)
            H[base, base+1] = vk * np.exp(-1j * phi)
            H[base+1, base] = np.conj(H[base, base+1])
            H[base+2, base+3] = -vk * np.exp(1j * phi)
            H[base+3, base+2] = np.conj(H[base+2, base+3])

            # Interlayer coupling
            if il < n_layers - 1:
                base2 = (il + 1) * 4
                for a in range(4):
                    H[base + a, base2 + a] = -t * 0.5
                    H[base2 + a, base + a] = -t * 0.5

        return H

    return H_slab
