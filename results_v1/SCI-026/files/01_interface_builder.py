"""
Interface Builder: Li6PS5Cl / LiCoO2 Interface Structural Modeling
===================================================================
Builds slab models for electrode/electrolyte interfaces, evaluates
lattice mismatch, and generates VASP POSCAR input files.

References:
  - Haruyama et al., Chem. Mater. 2014, 26, 4248
  - Kim et al., ACS Appl. Mater. Interfaces 2020, 12, 44, 49586
"""

import numpy as np
import json
import os
from dataclasses import dataclass, asdict
from typing import Tuple, List

# ---------------------------------------------------------------------------
# Crystal structure parameters (experimental lattice constants)
# ---------------------------------------------------------------------------
# Li6PS5Cl  – argyrodite, F-43m (#216)
LPS_LATTICE_A = 9.859   # Å (ICSD-246046)
LPS_LATTICE_ANGLES = (90.0, 90.0, 90.0)

# LiCoO2  – layered R-3m (#166),  hexagonal setting
LCO_LATTICE_A = 2.816   # Å
LCO_LATTICE_C = 14.054  # Å
LCO_LATTICE_ANGLES = (90.0, 90.0, 120.0)

# Li3PO4 (β-phase, Pmn21) – coating material
LPO_LATTICE = (6.115, 5.039, 4.847)  # a, b, c  Å

@dataclass
class InterfaceModel:
    name: str
    orientation_elec: str       # e.g.  "(110)"
    orientation_electrolyte: str
    supercell_elec: Tuple[int, int, int]
    supercell_electrolyte: Tuple[int, int, int]
    mismatch_a: float           # % along a
    mismatch_b: float           # % along b
    area_A2: float              # interface area Å²
    n_atoms: int
    vacuum_A: float             # vacuum slab thickness
    interface_distance_A: float # initial gap between slabs


def compute_lattice_mismatch(a1: float, a2: float) -> float:
    """Return percent mismatch 2|a1-a2|/(a1+a2)*100."""
    return 2 * abs(a1 - a2) / (a1 + a2) * 100


def surface_cell_LCO(hkl: Tuple[int,int,int]) -> Tuple[float, float]:
    """
    Approximate in-plane lattice constants for LiCoO2 surface.
    Uses hexagonal cell: a=b=2.816 Å, c=14.054 Å, γ=120°.
    """
    h, k, l = hkl
    a, c = LCO_LATTICE_A, LCO_LATTICE_C
    if (h, k, l) == (0, 0, 1):
        return a, a
    elif (h, k, l) == (1, 0, 4):  # most stable cleavage
        d = 1 / np.sqrt((4/3)*((h**2 + h*k + k**2)/a**2) + l**2/c**2)
        return a, a * np.sqrt(3)
    else:
        return a, a


def surface_cell_LPS(hkl: Tuple[int,int,int]) -> Tuple[float, float]:
    """Approximate in-plane constants for Li6PS5Cl cubic surface."""
    a = LPS_LATTICE_A
    h, k, l = hkl
    if (h, k, l) == (1, 0, 0):
        return a, a
    elif (h, k, l) == (1, 1, 0):
        return a, a * np.sqrt(2)
    elif (h, k, l) == (1, 1, 1):
        return a * np.sqrt(2), a * np.sqrt(2)
    return a, a


def build_interface_models() -> List[InterfaceModel]:
    """
    Enumerate low-mismatch interface configurations for
    Li6PS5Cl / LiCoO2 by testing common cleavage planes.
    """
    lco_orientations = [(0,0,1), (1,0,4), (1,0,0)]
    lps_orientations = [(1,0,0), (1,1,0), (1,1,1)]

    models = []
    for lco_hkl in lco_orientations:
        for lps_hkl in lps_orientations:
            a_lco, b_lco = surface_cell_LCO(lco_hkl)
            a_lps, b_lps = surface_cell_LPS(lps_hkl)

            # Try 1x1, 1x2, 2x1, 2x2 supercells to find <5% mismatch
            for sa in range(1, 8):
                for sb in range(1, 8):
                    for ta in range(1, 5):
                        for tb in range(1, 5):
                            A_lco = a_lco * sa
                            B_lco = b_lco * sb
                            A_lps = a_lps * ta
                            B_lps = b_lps * tb
                            mm_a = compute_lattice_mismatch(A_lco, A_lps)
                            mm_b = compute_lattice_mismatch(B_lco, B_lps)
                            if mm_a < 2.0 and mm_b < 2.0:
                                n_lco = sa * sb * 8    # atoms per layer × layers (4)
                                n_lps = ta * tb * 52   # Li6PS5Cl: 104 atoms/cell
                                area = A_lco * B_lco
                                models.append(InterfaceModel(
                                    name=(f"LCO{lco_hkl[0]}{lco_hkl[1]}{lco_hkl[2]}_"
                                          f"LPS{lps_hkl[0]}{lps_hkl[1]}{lps_hkl[2]}_"
                                          f"{sa}x{sb}_{ta}x{tb}"),
                                    orientation_elec=str(lco_hkl),
                                    orientation_electrolyte=str(lps_hkl),
                                    supercell_elec=(sa, sb, 4),
                                    supercell_electrolyte=(ta, tb, 2),
                                    mismatch_a=round(mm_a, 3),
                                    mismatch_b=round(mm_b, 3),
                                    area_A2=round(area, 2),
                                    n_atoms=n_lco + n_lps,
                                    vacuum_A=15.0,
                                    interface_distance_A=2.5,
                                ))
    return models


def select_optimal_models(models: List[InterfaceModel],
                          top_n: int = 6) -> List[InterfaceModel]:
    """Select by total mismatch (a²+b²)^0.5, then smallest cell."""
    scored = sorted(models,
                    key=lambda m: (m.mismatch_a**2 + m.mismatch_b**2)**0.5 + m.n_atoms/1000)
    return scored[:top_n]


def generate_poscar_header(model: InterfaceModel) -> str:
    """Generate POSCAR comment line for the interface model."""
    return (f"Interface: {model.name} | "
            f"Mismatch: {model.mismatch_a:.2f}%/{model.mismatch_b:.2f}% | "
            f"Area: {model.area_A2:.1f} Å²")


def write_vasp_incar_relax(outdir: str) -> None:
    """Write VASP INCAR for ionic+volume relaxation of interface slab."""
    incar = """\
# INCAR: Interface Geometry Relaxation
SYSTEM  = LiCoO2_Li6PS5Cl_interface

# Electronic minimization
ISTART  = 0
ICHARG  = 2
ENCUT   = 520        # eV  –  converged for both Li6PS5Cl and LiCoO2
PREC    = Accurate
EDIFF   = 1E-6       # eV  –  tight for accurate forces
NELM    = 150
ALGO    = Fast
LREAL   = Auto

# Spin-polarization (Co d-electrons)
ISPIN   = 2
MAGMOM  = 48*0.0 8*3.0 32*0.0   # example: adjust per supercell
LDAU    = .TRUE.
LDATYPE = 2
LDAUL   = -1 -1 3 -1   # L=3 for Co  (POTCAR order: Li P S Cl Co O)
LDAUU   = 0  0  0  0 3.3
LDAUJ   = 0  0  0  0 0.0

# Ionic relaxation
IBRION  = 2
NSW     = 300
ISIF    = 2          # relax ions only; fix cell shape
POTIM   = 0.02
EDIFFG  = -0.02      # eV/Å  –  force threshold

# k-points / smearing
ISMEAR  = 0
SIGMA   = 0.05

# Dispersion correction (vdW-D3)
IVDW    = 11

# Output
LWAVE   = .FALSE.
LCHARG  = .TRUE.
LORBIT  = 11

# Parallelization
NCORE   = 8
KPAR    = 2
"""
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "INCAR_relax"), "w") as f:
        f.write(incar)


def write_vasp_incar_static(outdir: str) -> None:
    """Write VASP INCAR for static SCF (charge density, DOS)."""
    incar = """\
# INCAR: Static SCF for charge density / LDOS
SYSTEM  = LiCoO2_Li6PS5Cl_interface_static

ISTART  = 1
ICHARG  = 1
ENCUT   = 520
PREC    = Accurate
EDIFF   = 1E-7
NELM    = 200
ALGO    = All

ISPIN   = 2
LDAU    = .TRUE.
LDATYPE = 2
LDAUL   = -1 -1 3 -1
LDAUU   = 0  0  0  0 3.3
LDAUJ   = 0  0  0  0 0.0

IBRION  = -1
NSW     = 0
ISIF    = 2
ISMEAR  = -5         # tetrahedron – accurate DOS
SIGMA   = 0.01

IVDW    = 11
LWAVE   = .TRUE.
LCHARG  = .TRUE.
LORBIT  = 12         # full LDOS on each ion
LAECHG  = .TRUE.     # all-electron charge → Bader analysis

NCORE   = 8
KPAR    = 2
"""
    with open(os.path.join(outdir, "INCAR_static"), "w") as f:
        f.write(incar)


def write_kpoints(outdir: str, kmesh: Tuple[int,int,int] = (3,3,1)) -> None:
    kpts = f"""\
Automatic Monkhorst-Pack
0
Gamma
{kmesh[0]}  {kmesh[1]}  {kmesh[2]}
0  0  0
"""
    with open(os.path.join(outdir, "KPOINTS"), "w") as f:
        f.write(kpts)


def write_potcar_spec(outdir: str) -> None:
    spec = """\
# POTCAR specification – use PAW_PBE pseudopotentials
# Recommended VASP POTCAR ordering:
#   Li_sv   (3 valence electrons)
#   P        (5 ve)
#   S        (6 ve)
#   Cl       (7 ve)
#   Co_pv    (15 ve)
#   O        (6 ve)
#
# Command to generate:
#   cat $VASP_PP/Li_sv/POTCAR \\
#       $VASP_PP/P/POTCAR     \\
#       $VASP_PP/S/POTCAR     \\
#       $VASP_PP/Cl/POTCAR    \\
#       $VASP_PP/Co_pv/POTCAR \\
#       $VASP_PP/O/POTCAR > POTCAR
"""
    with open(os.path.join(outdir, "POTCAR_spec.txt"), "w") as f:
        f.write(spec)


def main():
    outdir = "results/interface_models"
    os.makedirs(outdir, exist_ok=True)

    print("=" * 60)
    print("  Interface Builder: Li6PS5Cl / LiCoO2")
    print("=" * 60)

    # ------------------------------------------------------------------
    # Build & rank interface models
    # ------------------------------------------------------------------
    all_models = build_interface_models()
    best = select_optimal_models(all_models, top_n=6)

    print(f"\nTotal candidates screened: {len(all_models)}")
    print(f"Top {len(best)} low-mismatch models:\n")
    print(f"{'Model':<50} {'Mis_a%':>7} {'Mis_b%':>7} {'Area Å²':>10} {'N_atoms':>8}")
    print("-" * 85)
    for m in best:
        print(f"{m.name:<50} {m.mismatch_a:>7.3f} {m.mismatch_b:>7.3f} "
              f"{m.area_A2:>10.1f} {m.n_atoms:>8d}")

    # ------------------------------------------------------------------
    # Write VASP input files for the best model
    # ------------------------------------------------------------------
    vasp_dir = "results/vasp_inputs/interface_relax"
    write_vasp_incar_relax(vasp_dir)
    write_vasp_incar_static(vasp_dir)
    write_kpoints(vasp_dir, kmesh=(3, 3, 1))
    write_potcar_spec(vasp_dir)

    # ------------------------------------------------------------------
    # Save model summary as JSON
    # ------------------------------------------------------------------
    summary = {
        "total_candidates": len(all_models),
        "top_models": [asdict(m) for m in best],
        "recommended_model": asdict(best[0]),
        "lattice_parameters": {
            "Li6PS5Cl_a_A": LPS_LATTICE_A,
            "LiCoO2_a_A": LCO_LATTICE_A,
            "LiCoO2_c_A": LCO_LATTICE_C,
        },
        "POSCAR_comment": generate_poscar_header(best[0]),
    }
    with open(os.path.join(outdir, "interface_models.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nRecommended model: {best[0].name}")
    print(f"  Mismatch a: {best[0].mismatch_a:.3f}%  b: {best[0].mismatch_b:.3f}%")
    print(f"  Interface area: {best[0].area_A2:.1f} Å²")
    print(f"  Total atoms: {best[0].n_atoms}")
    print(f"\nVASP inputs written to: {vasp_dir}/")
    print(f"Model JSON: {outdir}/interface_models.json")


if __name__ == "__main__":
    main()
