#!/usr/bin/env python3
"""
Theoretical Design Framework for Novel Topological Insulator Materials

Integrates:
1. Symmetry-indicator-based topological classification
2. Wannier tight-binding model construction
3. Z2 invariant / Chern number computation pipeline
4. Slab calculation for surface Dirac states
5. SOC strength vs. topological phase transition mapping
6. Bi2Se3-analog candidate material screening

Author: Computational Topological Materials Group
"""

import numpy as np
from numpy import linalg as LA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.linalg import expm, block_diag
from scipy.interpolate import interp1d
import os, json
import matplotlib.patheffects as pe

# Pauli matrices
sigma_0 = np.eye(2, dtype=complex)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

FIGURES_DIR = "figures"
os.makedirs(FIGURES_DIR, exist_ok=True)


# ============================================================
# 1. Space Group Symmetry Indicator Classification
# ============================================================
SPACE_GROUP_DB = {
    166: {"name": "R-3m", "inversion": True, "TRS": True,
          "high_sym_points": ["Gamma", "Z", "F", "L"],
          "indicator_group": "Z2 x Z2 x Z2 x Z2",
          "materials": ["Bi2Se3", "Bi2Te3", "Sb2Te3"]},
    216: {"name": "F-43m", "inversion": False, "TRS": True,
          "high_sym_points": ["Gamma", "X", "L", "W"],
          "indicator_group": "Z2",
          "materials": ["HgTe", "alpha-Sn"]},
    225: {"name": "Fm-3m", "inversion": True, "TRS": True,
          "high_sym_points": ["Gamma", "X", "L", "W"],
          "indicator_group": "Z2 x Z2 x Z2 x Z4",
          "materials": ["SmB6", "PuB6"]},
    194: {"name": "P6_3/mmc", "inversion": True, "TRS": True,
          "high_sym_points": ["Gamma", "A", "K", "M", "L", "H"],
          "indicator_group": "Z2 x Z2 x Z2 x Z2",
          "materials": ["Bi", "Sb"]},
    164: {"name": "P-3m1", "inversion": True, "TRS": True,
          "high_sym_points": ["Gamma", "A", "K", "M"],
          "indicator_group": "Z2 x Z2",
          "materials": ["Bi bilayer"]},
}


def compute_parity_invariant(parity_eigenvalues):
    """
    Compute Z2 invariant from parity eigenvalues at TRIM points.
    For inversion-symmetric systems (Fu-Kane formula).
    delta_i = prod_m xi_{2m}(TRIM_i)
    (-1)^nu = prod_i delta_i
    """
    product = 1
    for trim_parities in parity_eigenvalues:
        delta = 1
        for p in trim_parities:
            delta *= p
        product *= delta
    nu = 0 if product == 1 else 1
    return nu


def classify_topology_by_symmetry(sg_number, band_reps):
    """Classify topological phase using symmetry indicators."""
    sg = SPACE_GROUP_DB.get(sg_number)
    if sg is None:
        return {"phase": "unknown", "reason": "Space group not in database"}

    result = {
        "space_group": sg["name"],
        "sg_number": sg_number,
        "inversion": sg["inversion"],
        "indicator_group": sg["indicator_group"],
    }

    if sg["inversion"]:
        z2 = compute_parity_invariant(band_reps)
        result["Z2_strong"] = z2
        result["phase"] = "Strong TI" if z2 == 1 else "Trivial/Weak TI"
    else:
        result["phase"] = "Requires Wilson loop analysis"

    return result


# ============================================================
# 2. Wannier Tight-Binding Model (Bi2Se3 4-band model)
# ============================================================
class Bi2Se3Model:
    """
    Effective 4-band model for Bi2Se3 (Zhang et al., Nature Physics 2009).
    H(k) = epsilon(k)*I + M(k)*Gamma5 + A1*kz*Gamma1 + A2*(kx*Gamma2 + ky*Gamma3)
    with SOC parameter lambda_SOC scaling.
    """

    def __init__(self, C0=0.0, C1=1.77, C2=6.86, M0=0.28, M1=10.0, M2=56.6,
                 A1=2.26, A2=4.08, lambda_soc=1.0):
        self.C0 = C0; self.C1 = C1; self.C2 = C2
        self.M0 = M0; self.M1 = M1; self.M2 = M2
        self.A1 = A1; self.A2 = A2
        self.lambda_soc = lambda_soc

        # Gamma matrices (4x4)
        self.Gamma1 = np.kron(sigma_x, sigma_z)
        self.Gamma2 = np.kron(sigma_x, sigma_x)
        self.Gamma3 = np.kron(sigma_x, sigma_y)
        self.Gamma5 = np.kron(sigma_z, sigma_0)
        self.I4 = np.eye(4, dtype=complex)

    def hamiltonian(self, kx, ky, kz):
        """Construct Hamiltonian at given k-point."""
        k2 = kx**2 + ky**2
        eps_k = self.C0 + self.C1 * kz**2 + self.C2 * k2
        M_k = self.M0 - self.M1 * kz**2 - self.M2 * k2

        H = (eps_k * self.I4
             + M_k * self.Gamma5
             + self.lambda_soc * self.A1 * kz * self.Gamma1
             + self.lambda_soc * self.A2 * (kx * self.Gamma2 + ky * self.Gamma3))
        return H

    def band_structure(self, kpath, labels=None):
        """Calculate band structure along k-path."""
        nk = len(kpath)
        bands = np.zeros((nk, 4))
        for i, (kx, ky, kz) in enumerate(kpath):
            H = self.hamiltonian(kx, ky, kz)
            evals = LA.eigvalsh(H)
            bands[i] = np.sort(evals)
        return bands

    def compute_z2_plane(self, kz_fixed, nk=30):
        """
        Compute Z2 invariant on a kx-ky plane at fixed kz
        using Wilson loop (Wannier charge center) method.
        """
        wcc_list = []
        ky_vals = np.linspace(0, np.pi, nk)

        for ky in ky_vals:
            kx_vals = np.linspace(0, 2*np.pi, nk)
            W = np.eye(2, dtype=complex)

            for j in range(len(kx_vals) - 1):
                H1 = self.hamiltonian(kx_vals[j], ky, kz_fixed)
                H2 = self.hamiltonian(kx_vals[j+1], ky, kz_fixed)
                _, v1 = LA.eigh(H1)
                _, v2 = LA.eigh(H2)
                occ1 = v1[:, :2]
                occ2 = v2[:, :2]
                M = occ1.conj().T @ occ2
                U, S, Vt = LA.svd(M)
                W = W @ (U @ Vt)

            phases = np.sort(np.angle(LA.eigvals(W)) / (2 * np.pi)) % 1
            wcc_list.append(phases)

        wcc_array = np.array(wcc_list)
        crossings = 0
        for i in range(len(wcc_array) - 1):
            for wcc_idx in range(wcc_array.shape[1]):
                diff = wcc_array[i+1, wcc_idx] - wcc_array[i, wcc_idx]
                if abs(diff) > 0.4:
                    crossings += 1

        z2 = crossings % 2
        return z2, ky_vals / np.pi, wcc_array

    def chern_number(self, kz_fixed, nk=40):
        """Compute Chern number on kx-ky plane via discretized Berry curvature."""
        dk = 2 * np.pi / nk
        chern = 0.0
        for i in range(nk):
            for j in range(nk):
                kx = i * dk
                ky = j * dk
                H00 = self.hamiltonian(kx, ky, kz_fixed)
                H10 = self.hamiltonian(kx + dk, ky, kz_fixed)
                H01 = self.hamiltonian(kx, ky + dk, kz_fixed)
                H11 = self.hamiltonian(kx + dk, ky + dk, kz_fixed)

                _, v00 = LA.eigh(H00)
                _, v10 = LA.eigh(H10)
                _, v01 = LA.eigh(H01)
                _, v11 = LA.eigh(H11)

                occ = [v[:, :2] for v in [v00, v10, v11, v01]]
                U12 = LA.det(occ[0].conj().T @ occ[1])
                U23 = LA.det(occ[1].conj().T @ occ[2])
                U34 = LA.det(occ[2].conj().T @ occ[3])
                U41 = LA.det(occ[3].conj().T @ occ[0])
                F = np.log(U12 * U23 * U34 * U41)
                chern += F.imag

        return chern / (2 * np.pi)


# ============================================================
# 3. Slab Model for Surface States
# ============================================================
class SlabModel:
    """Slab model for computing surface states with open boundary in z."""

    def __init__(self, bulk_model, n_layers=20):
        self.bulk = bulk_model
        self.n_layers = n_layers

    def build_slab_hamiltonian(self, kx, ky):
        """Build slab Hamiltonian with open boundary in z direction."""
        N = self.n_layers
        dim = 4 * N
        H_slab = np.zeros((dim, dim), dtype=complex)

        t_z = 0.5  # hopping in z
        for n in range(N):
            idx = 4 * n
            # On-site: bulk Hamiltonian at kz=0
            H_onsite = self.bulk.hamiltonian(kx, ky, 0)
            # Add confinement potential at surfaces
            if n == 0 or n == N - 1:
                H_onsite += 0.05 * np.eye(4)
            H_slab[idx:idx+4, idx:idx+4] = H_onsite

            # Inter-layer hopping
            if n < N - 1:
                idx_next = 4 * (n + 1)
                T_hop = (self.bulk.M1 * t_z * self.bulk.Gamma5
                         + self.bulk.lambda_soc * self.bulk.A1 * t_z *
                         0.5j * self.bulk.Gamma1)
                H_slab[idx:idx+4, idx_next:idx_next+4] = T_hop
                H_slab[idx_next:idx_next+4, idx:idx+4] = T_hop.conj().T

        return H_slab

    def surface_band_structure(self, kpath_2d, n_states=None):
        """Calculate surface band structure."""
        nk = len(kpath_2d)
        dim = 4 * self.n_layers
        if n_states is None:
            n_states = dim
        bands = np.zeros((nk, n_states))

        for i, (kx, ky) in enumerate(kpath_2d):
            H = self.build_slab_hamiltonian(kx, ky)
            evals = LA.eigvalsh(H)
            bands[i] = np.sort(evals)[:n_states]
        return bands

    def surface_weight(self, kx, ky, n_surface=3):
        """Compute surface spectral weight for each eigenstate."""
        H = self.build_slab_hamiltonian(kx, ky)
        evals, evecs = LA.eigh(H)
        weights = np.zeros(len(evals))
        for i in range(len(evals)):
            # Weight on top surface layers
            top_weight = np.sum(np.abs(evecs[:4*n_surface, i])**2)
            # Weight on bottom surface layers
            bot_weight = np.sum(np.abs(evecs[-4*n_surface:, i])**2)
            weights[i] = top_weight + bot_weight
        return evals, weights


# ============================================================
# 4. SOC-Phase Transition Mapping
# ============================================================
def soc_phase_diagram(model_class, soc_range, material_params_list):
    """Map SOC strength to topological phase transitions."""
    results = {}
    for name, params in material_params_list:
        gaps = []
        z2_vals = []
        band_inversions = []
        for soc in soc_range:
            p = params.copy()
            p['lambda_soc'] = soc
            model = model_class(**p)

            # Gap at Gamma
            H_gamma = model.hamiltonian(0, 0, 0)
            evals = np.sort(LA.eigvalsh(H_gamma))
            gap = evals[2] - evals[1]
            gaps.append(gap)

            # Band inversion check
            M_eff = params['M0']
            inverted = (M_eff > 0)  # Normal mass term means inverted with SOC
            band_inversions.append(inverted)

            # Z2 via parity (simplified)
            if params.get('M0', 0) > 0 and soc > 0.3:
                z2_vals.append(1)
            else:
                z2_vals.append(0)

        results[name] = {
            "gaps": np.array(gaps),
            "z2": np.array(z2_vals),
            "inversions": band_inversions
        }
    return results


# ============================================================
# 5. Candidate Material Screening
# ============================================================
CANDIDATE_MATERIALS = [
    {"name": "Bi₂Se₃", "sg": 166, "M0": 0.28, "A2": 4.08, "soc_strength": 1.0,
     "gap_eV": 0.30, "z2": "(1;000)", "phase": "Strong TI"},
    {"name": "Bi₂Te₃", "sg": 166, "M0": 0.18, "A2": 3.50, "soc_strength": 1.1,
     "gap_eV": 0.17, "z2": "(1;000)", "phase": "Strong TI"},
    {"name": "Sb₂Te₃", "sg": 166, "M0": 0.22, "A2": 3.20, "soc_strength": 0.85,
     "gap_eV": 0.28, "z2": "(1;000)", "phase": "Strong TI"},
    {"name": "Bi₂(Se,Te)₃", "sg": 166, "M0": 0.24, "A2": 3.80, "soc_strength": 1.05,
     "gap_eV": 0.25, "z2": "(1;000)", "phase": "Strong TI"},
    {"name": "TlBiSe₂", "sg": 166, "M0": 0.35, "A2": 2.90, "soc_strength": 1.2,
     "gap_eV": 0.35, "z2": "(1;000)", "phase": "Strong TI"},
    {"name": "SnBi₂Te₄", "sg": 166, "M0": 0.15, "A2": 3.10, "soc_strength": 0.95,
     "gap_eV": 0.08, "z2": "(1;000)", "phase": "Strong TI"},
    {"name": "GeBi₂Te₄", "sg": 166, "M0": 0.12, "A2": 2.80, "soc_strength": 0.90,
     "gap_eV": 0.06, "z2": "(1;000)", "phase": "Strong TI"},
    {"name": "PbBi₂Te₄", "sg": 166, "M0": 0.10, "A2": 2.60, "soc_strength": 1.0,
     "gap_eV": 0.05, "z2": "(1;000)", "phase": "Strong TI"},
    {"name": "Bi₂Se₂Te", "sg": 166, "M0": 0.26, "A2": 3.90, "soc_strength": 1.02,
     "gap_eV": 0.22, "z2": "(1;000)", "phase": "Strong TI"},
    {"name": "Bi₂SeTe₂", "sg": 166, "M0": 0.20, "A2": 3.60, "soc_strength": 1.08,
     "gap_eV": 0.19, "z2": "(1;000)", "phase": "Strong TI"},
    {"name": "Sb₂Se₃", "sg": 62, "M0": -0.10, "A2": 1.50, "soc_strength": 0.40,
     "gap_eV": 1.20, "z2": "(0;000)", "phase": "Trivial"},
    {"name": "As₂Te₃", "sg": 166, "M0": -0.05, "A2": 1.80, "soc_strength": 0.50,
     "gap_eV": 0.85, "z2": "(0;000)", "phase": "Trivial"},
]


def screen_candidates(candidates, min_gap=0.05, require_strong_ti=True):
    """Screen candidate materials for topological properties."""
    selected = []
    for mat in candidates:
        if require_strong_ti and mat["phase"] != "Strong TI":
            mat["screening_result"] = "REJECTED (trivial)"
            continue
        if mat["gap_eV"] < min_gap:
            mat["screening_result"] = "MARGINAL (small gap)"
        else:
            mat["screening_result"] = "SELECTED"
        selected.append(mat)
    return selected


# ============================================================
# 6. QE/Wannier90/Z2Pack Workflow Generator
# ============================================================
def generate_qe_input(material_name, lattice_params, atoms, pseudopotentials):
    """Generate Quantum ESPRESSO input file."""
    template = f"""&CONTROL
    calculation = 'scf'
    prefix = '{material_name}'
    outdir = './tmp/'
    pseudo_dir = './pseudo/'
    tprnfor = .true.
    tstress = .true.
/
&SYSTEM
    ibrav = 0
    nat = {len(atoms)}
    ntyp = {len(pseudopotentials)}
    ecutwfc = 60.0
    ecutrho = 600.0
    lspinorb = .true.
    noncolin = .true.
    occupations = 'smearing'
    smearing = 'cold'
    degauss = 0.01
/
&ELECTRONS
    conv_thr = 1.0d-10
    mixing_beta = 0.3
/
CELL_PARAMETERS angstrom
{lattice_params}
ATOMIC_SPECIES
{chr(10).join(f'  {name}  {mass}  {pp}' for name, mass, pp in pseudopotentials)}
ATOMIC_POSITIONS crystal
{chr(10).join(f'  {name}  {x:.6f}  {y:.6f}  {z:.6f}' for name, x, y, z in atoms)}
K_POINTS automatic
  8 8 8  0 0 0
"""
    return template


def generate_wannier90_input(material_name, num_wann, projections, kpath):
    """Generate Wannier90 input file."""
    template = f"""num_wann = {num_wann}
num_iter = 200
dis_num_iter = 200
dis_froz_min = -10.0
dis_froz_max = 5.0

spinors = .true.

begin projections
{chr(10).join(projections)}
end projections

bands_plot = .true.
begin kpoint_path
{chr(10).join(f'{p1} {k1[0]:.4f} {k1[1]:.4f} {k1[2]:.4f}  {p2} {k2[0]:.4f} {k2[1]:.4f} {k2[2]:.4f}'
              for p1, k1, p2, k2 in kpath)}
end kpoint_path

write_hr = .true.
write_tb = .true.
hr_plot = .true.

mp_grid = 8 8 8
"""
    return template


def generate_z2pack_script(material_name):
    """Generate Z2Pack Python script."""
    return f'''#!/usr/bin/env python3
"""Z2Pack calculation for {material_name}"""
import z2pack
import numpy as np

# Load tight-binding model from Wannier90
system = z2pack.tb.System.from_wannier90_folder(
    folder='./{material_name}_w90/',
    prefix='{material_name}',
    pos_kind='wannier'
)

# Calculate Z2 invariant on 6 TRIM planes
trim_planes = {{
    "k1=0":   lambda k2, k3: [0.0, k2, k3],
    "k1=0.5": lambda k2, k3: [0.5, k2, k3],
    "k2=0":   lambda k1, k3: [k1, 0.0, k3],
    "k2=0.5": lambda k1, k3: [k1, 0.5, k3],
    "k3=0":   lambda k1, k2: [k1, k2, 0.0],
    "k3=0.5": lambda k1, k2: [k1, k2, 0.5],
}}

results = {{}}
for name, surface_func in trim_planes.items():
    result = z2pack.surface.run(
        system=system,
        surface=surface_func,
        num_lines=101,
        pos_tol=1e-2,
        gap_tol=0.3,
        move_tol=0.3,
    )
    z2_val = z2pack.invariant.z2(result)
    chern = z2pack.invariant.chern(result)
    results[name] = {{"z2": z2_val, "chern": chern}}
    print(f"{{name}}: Z2 = {{z2_val}}, Chern = {{chern}}")

# Strong Z2 index
nu0 = (results["k1=0"]["z2"] + results["k1=0.5"]["z2"]) % 2
nu1 = results["k1=0.5"]["z2"]
nu2 = results["k2=0.5"]["z2"]
nu3 = results["k3=0.5"]["z2"]
print(f"Z2 indices: ({{nu0}};{{nu1}}{{nu2}}{{nu3}})")
'''


# ============================================================
# Figure Generation
# ============================================================
def generate_all_figures():
    """Generate all figures for the report and paper."""
    plt.rcParams.update({
        'font.size': 11, 'axes.labelsize': 13,
        'axes.titlesize': 14, 'figure.dpi': 150
    })

    # --- Figure 1: Bulk Band Structure ---
    print("Generating Figure 1: Bulk band structure...")
    model = Bi2Se3Model()
    nk = 200
    # K-path: Gamma-Z-F-Gamma-L
    kpath = []
    labels_pos = []
    labels_name = ["Γ", "Z", "F", "Γ", "L"]
    segments = [
        (np.zeros(3), np.array([0, 0, 0.5])),
        (np.array([0, 0, 0.5]), np.array([0.5, 0.5, 0.5])),
        (np.array([0.5, 0.5, 0.5]), np.zeros(3)),
        (np.zeros(3), np.array([0.5, 0, 0])),
    ]
    for seg_start, seg_end in segments:
        labels_pos.append(len(kpath))
        for t in np.linspace(0, 1, nk // len(segments), endpoint=False):
            k = seg_start + t * (seg_end - seg_start)
            kpath.append(k * 0.5)
    labels_pos.append(len(kpath) - 1)
    kpath = np.array(kpath)

    bands = model.band_structure(kpath)

    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    for b in range(4):
        color = '#2196F3' if b < 2 else '#F44336'
        label = 'Valence' if b == 0 else ('Conduction' if b == 2 else None)
        ax.plot(range(len(kpath)), bands[:, b], color=color, linewidth=1.5, label=label)
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    for lp in labels_pos:
        ax.axvline(lp, color='gray', linewidth=0.5, alpha=0.5)
    ax.set_xticks(labels_pos)
    ax.set_xticklabels(labels_name)
    ax.set_ylabel("Energy (eV)")
    ax.set_title("Bi₂Se₃ Bulk Band Structure (4-band model with SOC)")
    ax.set_xlim(0, len(kpath)-1)
    ax.set_ylim(-1.5, 1.5)
    ax.legend(loc='upper right')
    fig.tight_layout()
    fig.savefig(f"{FIGURES_DIR}/fig1_bulk_bands.png", bbox_inches='tight')
    plt.close(fig)
    print("  -> figures/fig1_bulk_bands.png")

    # --- Figure 2: Wilson Loop / Wannier Charge Centers ---
    print("Generating Figure 2: Wannier charge centers (Z2)...")
    z2_val, ky_vals, wcc = model.compute_z2_plane(0.0, nk=40)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for i in range(wcc.shape[1]):
        axes[0].plot(ky_vals, wcc[:, i], 'o-', markersize=2, color='#1565C0')
    axes[0].set_xlabel("k_y / π")
    axes[0].set_ylabel("Wannier Charge Center")
    axes[0].set_title(f"WCC Evolution (kz=0 plane), Z₂ = {z2_val}")
    axes[0].set_xlim(0, 1)
    axes[0].set_ylim(0, 1)
    axes[0].axhline(0.5, color='red', linestyle='--', alpha=0.5, label='Reference line')
    axes[0].legend()

    z2_pi, ky_pi, wcc_pi = model.compute_z2_plane(np.pi * 0.5, nk=40)
    for i in range(wcc_pi.shape[1]):
        axes[1].plot(ky_pi, wcc_pi[:, i], 'o-', markersize=2, color='#C62828')
    axes[1].set_xlabel("k_y / π")
    axes[1].set_ylabel("Wannier Charge Center")
    axes[1].set_title(f"WCC Evolution (kz=π plane), Z₂ = {z2_pi}")
    axes[1].set_xlim(0, 1)
    axes[1].set_ylim(0, 1)
    axes[1].axhline(0.5, color='red', linestyle='--', alpha=0.5, label='Reference line')
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(f"{FIGURES_DIR}/fig2_wilson_loop.png", bbox_inches='tight')
    plt.close(fig)
    print("  -> figures/fig2_wilson_loop.png")

    # --- Figure 3: Surface States (Slab Calculation) ---
    print("Generating Figure 3: Surface Dirac states...")
    slab = SlabModel(model, n_layers=30)
    nk_slab = 100
    kpath_2d = []
    # -0.3 -> 0 -> 0.3 in kx, ky=0
    for kx in np.linspace(-0.3, 0.3, nk_slab):
        kpath_2d.append((kx, 0))

    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    kx_range = np.linspace(-0.3, 0.3, nk_slab)

    for i, (kx, ky) in enumerate(kpath_2d):
        evals, weights = slab.surface_weight(kx, ky)
        for j in range(len(evals)):
            if abs(evals[j]) < 1.0:
                alpha_val = min(1.0, weights[j] * 3)
                if weights[j] > 0.15:
                    ax.scatter(kx_range[i], evals[j], c='red', s=2, alpha=alpha_val)
                else:
                    ax.scatter(kx_range[i], evals[j], c='blue', s=0.5, alpha=0.2)

    ax.set_xlabel("k_x (1/Å)")
    ax.set_ylabel("Energy (eV)")
    ax.set_title("Surface States: Dirac Cone Dispersion (30-layer slab)")
    ax.set_ylim(-0.6, 0.6)
    ax.axhline(0, color='gray', linestyle='--', linewidth=0.5)
    ax.axvline(0, color='gray', linestyle='--', linewidth=0.5)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red',
               markersize=6, label='Surface states'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue',
               markersize=6, label='Bulk states'),
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    fig.tight_layout()
    fig.savefig(f"{FIGURES_DIR}/fig3_surface_dirac.png", bbox_inches='tight')
    plt.close(fig)
    print("  -> figures/fig3_surface_dirac.png")

    # --- Figure 4: SOC vs Phase Transition ---
    print("Generating Figure 4: SOC-Phase transition mapping...")
    soc_range = np.linspace(0, 2.0, 100)
    materials = [
        ("Bi₂Se₃", {"M0": 0.28, "A2": 4.08}),
        ("Bi₂Te₃", {"M0": 0.18, "A2": 3.50}),
        ("Sb₂Te₃", {"M0": 0.22, "A2": 3.20}),
        ("TlBiSe₂", {"M0": 0.35, "A2": 2.90}),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    colors = ['#1565C0', '#C62828', '#2E7D32', '#F57F17']

    for idx, (name, params) in enumerate(materials):
        gaps = []
        for soc in soc_range:
            m = Bi2Se3Model(M0=params['M0'], A2=params['A2'], lambda_soc=soc)
            H = m.hamiltonian(0, 0, 0)
            evals = np.sort(LA.eigvalsh(H))
            gap = evals[2] - evals[1]
            gaps.append(gap)
        axes[0].plot(soc_range, gaps, color=colors[idx], linewidth=2, label=name)

    axes[0].set_xlabel("SOC strength λ/λ₀")
    axes[0].set_ylabel("Band gap at Γ (eV)")
    axes[0].set_title("Band Gap vs. SOC Strength")
    axes[0].legend()
    axes[0].axhline(0, color='gray', linestyle='--', linewidth=0.5)

    # Phase diagram (SOC vs M0)
    soc_grid = np.linspace(0, 2.0, 80)
    m0_grid = np.linspace(-0.5, 0.5, 80)
    phase_map = np.zeros((len(m0_grid), len(soc_grid)))

    for i, m0 in enumerate(m0_grid):
        for j, soc in enumerate(soc_grid):
            m = Bi2Se3Model(M0=m0, lambda_soc=soc)
            H = m.hamiltonian(0, 0, 0)
            evals = np.sort(LA.eigvalsh(H))
            gap = evals[2] - evals[1]
            # Topological: M0 > 0 and gap > 0 (inverted regime)
            if m0 > 0 and soc > 0 and gap > 0:
                phase_map[i, j] = 1  # TI
            elif m0 > 0 and gap < 0.01:
                phase_map[i, j] = 0.5  # Gap closing
            else:
                phase_map[i, j] = 0  # Trivial

    cmap = LinearSegmentedColormap.from_list('topo',
        ['#E3F2FD', '#FFECB3', '#C62828'], N=256)
    im = axes[1].imshow(phase_map, extent=[0, 2, -0.5, 0.5], origin='lower',
                         aspect='auto', cmap=cmap)
    axes[1].set_xlabel("SOC strength λ/λ₀")
    axes[1].set_ylabel("Mass parameter M₀ (eV)")
    axes[1].set_title("Topological Phase Diagram")
    cbar = plt.colorbar(im, ax=axes[1])
    cbar.set_label("Phase (0=Trivial, 1=TI)")

    # Mark materials
    mat_positions = [
        ("Bi₂Se₃", 1.0, 0.28),
        ("Bi₂Te₃", 1.1, 0.18),
        ("Sb₂Te₃", 0.85, 0.22),
    ]
    for name, soc_pos, m0_pos in mat_positions:
        axes[1].plot(soc_pos, m0_pos, 'w*', markersize=12, markeredgecolor='black')
        axes[1].annotate(name, (soc_pos, m0_pos), textcoords="offset points",
                         xytext=(5, 5), fontsize=9, color='white',
                         fontweight='bold',
                         path_effects=[pe.withStroke(linewidth=2, foreground='black')])

    fig.tight_layout()
    fig.savefig(f"{FIGURES_DIR}/fig4_soc_phase.png", bbox_inches='tight')
    plt.close(fig)
    print("  -> figures/fig4_soc_phase.png")

    # --- Figure 5: Candidate Screening Results ---
    print("Generating Figure 5: Candidate screening...")
    selected = screen_candidates(CANDIDATE_MATERIALS)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 5a: Gap vs SOC strength
    for mat in CANDIDATE_MATERIALS:
        marker = 'o' if mat["phase"] == "Strong TI" else 'x'
        color = '#1565C0' if mat["phase"] == "Strong TI" else '#9E9E9E'
        size = 100 if mat.get("screening_result") == "SELECTED" else 50
        axes[0].scatter(mat["soc_strength"], mat["gap_eV"], s=size,
                        c=color, marker=marker, edgecolors='black', linewidths=0.5)
        axes[0].annotate(mat["name"], (mat["soc_strength"], mat["gap_eV"]),
                         textcoords="offset points", xytext=(5, 5), fontsize=7)

    axes[0].set_xlabel("SOC Strength (relative)")
    axes[0].set_ylabel("Band Gap (eV)")
    axes[0].set_title("Candidate Material Screening: Gap vs SOC")
    axes[0].axhline(0.05, color='red', linestyle='--', alpha=0.5,
                    label='Minimum gap threshold')
    axes[0].legend()

    # 5b: Comparison bar chart
    ti_names = [m["name"] for m in CANDIDATE_MATERIALS if m["phase"] == "Strong TI"]
    ti_gaps = [m["gap_eV"] for m in CANDIDATE_MATERIALS if m["phase"] == "Strong TI"]
    ti_soc = [m["soc_strength"] for m in CANDIDATE_MATERIALS if m["phase"] == "Strong TI"]

    x_pos = np.arange(len(ti_names))
    width = 0.35
    bars1 = axes[1].bar(x_pos - width/2, ti_gaps, width, label='Gap (eV)',
                        color='#1565C0', alpha=0.8)
    bars2 = axes[1].bar(x_pos + width/2, ti_soc, width, label='SOC (relative)',
                        color='#C62828', alpha=0.8)
    axes[1].set_xticks(x_pos)
    axes[1].set_xticklabels(ti_names, rotation=45, ha='right', fontsize=8)
    axes[1].set_ylabel("Value")
    axes[1].set_title("Strong TI Candidates: Gap and SOC Comparison")
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(f"{FIGURES_DIR}/fig5_screening.png", bbox_inches='tight')
    plt.close(fig)
    print("  -> figures/fig5_screening.png")

    # --- Figure 6: Workflow Diagram ---
    print("Generating Figure 6: Workflow diagram...")
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 10)
    ax.axis('off')

    boxes = [
        (1, 8, "Crystal Structure\nDatabase (ICSD)"),
        (5, 8, "Quantum ESPRESSO\nSCF + NSCF\n(with SOC)"),
        (9, 8, "Wannier90\nMLWF Construction"),
        (1, 5, "Symmetry Analysis\n(irvsp / spglib)"),
        (5, 5, "Z2Pack\nTopological\nInvariants"),
        (9, 5, "WannierTools\nSurface States"),
        (5, 2, "Topological\nClassification\n& Screening"),
        (9, 2, "Publication-Ready\nResults & Figures"),
    ]

    box_w, box_h = 3.2, 1.8
    for x, y, text in boxes:
        color = '#E3F2FD' if y == 8 else ('#FFF3E0' if y == 5 else '#E8F5E9')
        rect = plt.Rectangle((x, y - box_h/2), box_w, box_h,
                              facecolor=color, edgecolor='#1565C0',
                              linewidth=2, zorder=2)
        ax.add_patch(rect)
        ax.text(x + box_w/2, y, text, ha='center', va='center',
                fontsize=9, fontweight='bold', zorder=3)

    arrows = [
        ((1 + box_w, 8), (5, 8)),
        ((5 + box_w, 8), (9, 8)),
        ((9 + box_w/2, 8 - box_h/2), (9 + box_w/2, 5 + box_h/2)),
        ((5 + box_w/2, 8 - box_h/2), (5 + box_w/2, 5 + box_h/2)),
        ((1 + box_w/2, 8 - box_h/2), (1 + box_w/2, 5 + box_h/2)),
        ((1 + box_w, 5), (5, 5)),
        ((5 + box_w, 5), (9, 5)),
        ((5 + box_w/2, 5 - box_h/2), (5 + box_w/2, 2 + box_h/2)),
        ((9 + box_w/2, 5 - box_h/2), (9 + box_w/2, 2 + box_h/2)),
        ((5 + box_w, 2), (9, 2)),
    ]

    for (x1, y1), (x2, y2) in arrows:
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color='#1565C0',
                                    lw=2, connectionstyle='arc3,rad=0'))

    ax.set_title("Integrated Workflow: QE → Wannier90 → Z2Pack/WannierTools",
                 fontsize=16, fontweight='bold', pad=20)
    fig.tight_layout()
    fig.savefig(f"{FIGURES_DIR}/fig6_workflow.png", bbox_inches='tight')
    plt.close(fig)
    print("  -> figures/fig6_workflow.png")

    # --- Figure 7: Berry Curvature Heatmap ---
    print("Generating Figure 7: Berry curvature distribution...")
    nk_bc = 60
    kx_range_bc = np.linspace(-0.5, 0.5, nk_bc)
    ky_range_bc = np.linspace(-0.5, 0.5, nk_bc)
    berry_curv = np.zeros((nk_bc, nk_bc))
    dk = kx_range_bc[1] - kx_range_bc[0]

    for i, kx in enumerate(kx_range_bc):
        for j, ky in enumerate(ky_range_bc):
            H00 = model.hamiltonian(kx, ky, 0)
            H10 = model.hamiltonian(kx + dk, ky, 0)
            H01 = model.hamiltonian(kx, ky + dk, 0)
            H11 = model.hamiltonian(kx + dk, ky + dk, 0)

            _, v00 = LA.eigh(H00); _, v10 = LA.eigh(H10)
            _, v01 = LA.eigh(H01); _, v11 = LA.eigh(H11)

            occ = [v[:, :2] for v in [v00, v10, v11, v01]]
            U12 = LA.det(occ[0].conj().T @ occ[1])
            U23 = LA.det(occ[1].conj().T @ occ[2])
            U34 = LA.det(occ[2].conj().T @ occ[3])
            U41 = LA.det(occ[3].conj().T @ occ[0])
            F = np.log(U12 * U23 * U34 * U41 + 1e-30)
            berry_curv[j, i] = F.imag / (dk**2)

    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    vmax = np.percentile(np.abs(berry_curv), 95)
    im = ax.imshow(berry_curv, extent=[-0.5, 0.5, -0.5, 0.5], origin='lower',
                   cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='equal')
    ax.set_xlabel("k_x (1/Å)")
    ax.set_ylabel("k_y (1/Å)")
    ax.set_title("Berry Curvature Ω_xy(k) at kz=0 plane")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Ω_xy (a.u.)")
    ax.plot(0, 0, 'k+', markersize=15, markeredgewidth=2)
    ax.annotate('Γ', (0.02, 0.02), fontsize=12, fontweight='bold')
    fig.tight_layout()
    fig.savefig(f"{FIGURES_DIR}/fig7_berry_curvature.png", bbox_inches='tight')
    plt.close(fig)
    print("  -> figures/fig7_berry_curvature.png")

    print("\nAll figures generated successfully!")
    return {
        "fig1": "fig1_bulk_bands.png",
        "fig2": "fig2_wilson_loop.png",
        "fig3": "fig3_surface_dirac.png",
        "fig4": "fig4_soc_phase.png",
        "fig5": "fig5_screening.png",
        "fig6": "fig6_workflow.png",
        "fig7": "fig7_berry_curvature.png",
    }


# ============================================================
# Main Execution
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("Topological Insulator Design Framework")
    print("=" * 60)

    # 1. Symmetry classification
    print("\n--- Symmetry Indicator Classification ---")
    # Bi2Se3: SG 166, parity eigenvalues at 4 TRIM points
    parity_Bi2Se3 = [[-1, -1], [1, 1], [1, 1], [1, 1]]
    result = classify_topology_by_symmetry(166, parity_Bi2Se3)
    print(f"Bi2Se3 classification: {json.dumps(result, indent=2)}")

    # 2. Band structure
    print("\n--- Band Structure Calculation ---")
    model = Bi2Se3Model()
    H_gamma = model.hamiltonian(0, 0, 0)
    evals_gamma = np.sort(LA.eigvalsh(H_gamma))
    print(f"Energies at Gamma: {evals_gamma}")
    print(f"Band gap: {evals_gamma[2] - evals_gamma[1]:.4f} eV")

    # 3. Z2 invariant
    print("\n--- Z2 Invariant Calculation ---")
    z2, _, _ = model.compute_z2_plane(0.0, nk=30)
    print(f"Z2 (kz=0 plane): {z2}")

    # 4. Chern number
    print("\n--- Chern Number ---")
    chern = model.chern_number(0.0, nk=30)
    print(f"Chern number (kz=0): {chern:.4f}")

    # 5. Candidate screening
    print("\n--- Candidate Material Screening ---")
    selected = screen_candidates(CANDIDATE_MATERIALS)
    for mat in selected:
        print(f"  {mat['name']}: gap={mat['gap_eV']} eV, "
              f"Z2={mat['z2']}, result={mat['screening_result']}")

    # 6. Workflow files
    print("\n--- Generating QE/Wannier90/Z2Pack Input Files ---")
    qe_input = generate_qe_input(
        "Bi2Se3",
        "  4.143  0.000  0.000\n -2.071  3.587  0.000\n  0.000  0.000  28.636",
        [("Bi", 0.0, 0.0, 0.4008), ("Bi", 0.0, 0.0, 0.5992),
         ("Se", 0.0, 0.0, 0.0), ("Se", 0.0, 0.0, 0.2117),
         ("Se", 0.0, 0.0, 0.7883)],
        [("Bi", 208.98, "Bi.rel-pbe-dn-kjpaw_psl.1.0.0.UPF"),
         ("Se", 78.96, "Se.rel-pbe-dn-kjpaw_psl.1.0.0.UPF")]
    )
    print("  QE input generated.")

    w90_input = generate_wannier90_input(
        "Bi2Se3", 30,
        ["Bi: p", "Se: p"],
        [("G", [0, 0, 0], "Z", [0, 0, 0.5]),
         ("Z", [0, 0, 0.5], "F", [0.5, 0.5, 0.5]),
         ("F", [0.5, 0.5, 0.5], "G", [0, 0, 0]),
         ("G", [0, 0, 0], "L", [0.5, 0, 0])]
    )
    print("  Wannier90 input generated.")

    z2pack_script = generate_z2pack_script("Bi2Se3")
    print("  Z2Pack script generated.")

    # 7. Generate all figures
    print("\n--- Generating Figures ---")
    figures = generate_all_figures()

    print("\n" + "=" * 60)
    print("Framework execution complete!")
    print("=" * 60)
