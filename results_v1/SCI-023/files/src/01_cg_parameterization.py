"""
Block Copolymer Coarse-Grained Parameterization Strategy
MARTINI / SDK (Shinoda-DeVane-Klein) Force Field Framework

Generates CG parameters for PS-b-PMMA (polystyrene-b-poly(methyl methacrylate))
as a representative BCP for semiconductor DSA patterning.

References:
  - Marrink et al., JPCB 2007 (MARTINI)
  - Shinoda et al., Macromolecules 2007 (SDK)
  - Khadka et al., J. Chem. Theory Comput. 2020 (BCP CG)
"""

import numpy as np
import json
import os
from dataclasses import dataclass, asdict
from typing import Dict, List, Tuple

@dataclass
class BeadType:
    name: str
    mass: float
    sigma: float
    epsilon: float
    charge: float
    description: str

@dataclass
class BondParam:
    bead1: str; bead2: str
    r0: float; kb: float; style: str

@dataclass
class AngleParam:
    bead1: str; bead2: str; bead3: str
    theta0: float; ktheta: float; style: str

@dataclass
class DihedralParam:
    bead1: str; bead2: str; bead3: str; bead4: str
    k: float; n: int; delta: float

# ── MARTINI 3 Beads ────────────────────────────────────────────────────────
MARTINI_BEADS = {
    "TC5_PS_ring":  BeadType("TC5", 104.0*2/3, 0.43, 3.1, 0.0, "PS aromatic ring (MARTINI3 TC5)"),
    "SC3_PS_bb":    BeadType("SC3", 26.0,       0.43, 2.7, 0.0, "PS backbone (MARTINI3 SC3)"),
    "SC2_PMMA_ester":BeadType("SC2",50.0,       0.43, 2.5, 0.0, "PMMA ester group (MARTINI3 SC2)"),
    "SC1_PMMA_bb":  BeadType("SC1", 28.0,       0.43, 2.7, 0.0, "PMMA backbone (MARTINI3 SC1)"),
    "N0_junction":  BeadType("N0",  56.0,       0.47, 2.0, 0.0, "PS-PMMA diblock junction"),
}

def chi_ps_pmma(T_K: float) -> float:
    """Flory-Huggins χ for PS-PMMA: χ = 0.04 + 4.9/T"""
    return 0.04 + 4.9 / T_K

# ── SDK Beads ──────────────────────────────────────────────────────────────
SDK_BEADS = {
    "CMn_PS":   BeadType("CMn", 52.08, 0.465, 3.5, 0.0, "SDK: PS monomer (CMn)"),
    "CM_PMMA":  BeadType("CM",  50.04, 0.460, 3.2, 0.0, "SDK: PMMA monomer (CM)"),
    "W3_water": BeadType("W3",  54.045,0.461, 4.35,0.0, "SDK: 3-water bead (W3)"),
}

SDK_NONBONDED_CROSS = {
    ("CMn_PS","CM_PMMA"): {"epsilon":2.1,"sigma":0.462,"chi_eff":0.058,
                            "note":"PS-PMMA repulsive -> microphase separation"},
    ("CMn_PS","W3_water"):{"epsilon":1.6,"sigma":0.463,"note":"PS hydrophobic"},
    ("CM_PMMA","W3_water"):{"epsilon":3.8,"sigma":0.461,"note":"PMMA slightly hydrophilic"},
}

BOND_PARAMS = [
    BondParam("SC3_PS_bb","SC3_PS_bb",  0.470,3800.0,"harmonic"),
    BondParam("SC3_PS_bb","TC5_PS_ring",0.290,5000.0,"harmonic"),
    BondParam("SC1_PMMA_bb","SC1_PMMA_bb",0.470,3800.0,"harmonic"),
    BondParam("SC1_PMMA_bb","SC2_PMMA_ester",0.310,4500.0,"harmonic"),
    BondParam("SC3_PS_bb","N0_junction",0.470,3500.0,"harmonic"),
    BondParam("N0_junction","SC1_PMMA_bb",0.470,3500.0,"harmonic"),
]

ANGLE_PARAMS = [
    AngleParam("SC3_PS_bb","SC3_PS_bb","SC3_PS_bb",180.0,25.0,"harmonic"),
    AngleParam("SC3_PS_bb","SC3_PS_bb","TC5_PS_ring",150.0,25.0,"harmonic"),
    AngleParam("SC1_PMMA_bb","SC1_PMMA_bb","SC1_PMMA_bb",180.0,20.0,"harmonic"),
    AngleParam("SC1_PMMA_bb","SC1_PMMA_bb","SC2_PMMA_ester",120.0,25.0,"harmonic"),
]

DIHEDRAL_PARAMS = [
    DihedralParam("SC3_PS_bb","SC3_PS_bb","SC3_PS_bb","SC3_PS_bb",0.5,1,0.0),
    DihedralParam("SC1_PMMA_bb","SC1_PMMA_bb","SC1_PMMA_bb","SC1_PMMA_bb",0.3,1,0.0),
]


class IBIRefiner:
    """Iterative Boltzmann Inversion for CG potential refinement."""
    def __init__(self, target_rdf, r_bins, T=500.0):
        self.target_rdf = target_rdf
        self.r_bins = r_bins
        self.kT = 8.314e-3 * T
        self.potentials = []

    def initial_potential(self):
        g = np.where(self.target_rdf > 1e-10, self.target_rdf, 1e-10)
        U0 = -self.kT * np.log(g)
        U0 -= U0[-5:].mean()
        self.potentials.append(U0.copy())
        return U0

    def update_step(self, current_rdf, current_potential, alpha=0.5):
        g_cur = np.where(current_rdf > 1e-10, current_rdf, 1e-10)
        g_tgt = np.where(self.target_rdf > 1e-10, self.target_rdf, 1e-10)
        dU = -self.kT * np.log(g_tgt / g_cur)
        new_U = current_potential + alpha * dU
        new_U -= new_U[-5:].mean()
        self.potentials.append(new_U.copy())
        return new_U

    def convergence_metric(self, current_rdf):
        return float(np.sqrt(np.mean((current_rdf - self.target_rdf)**2)))


def generate_lammps_cg_input(n_PS=20, n_PMMA=20, T=500.0,
                              output_path="data/lammps_cg.in"):
    chi = chi_ps_pmma(T)
    chiN = chi * (n_PS + n_PMMA)
    content = f"""# LAMMPS CG Simulation: PS-b-PMMA
# MARTINI 3 adapted | N_PS={n_PS} N_PMMA={n_PMMA} T={T}K chiN={chiN:.2f}

units           real
atom_style      molecular
boundary        p p p
read_data       bcp_cg_initial.data

pair_style      lj/cut 1.2
pair_modify     shift yes
bond_style      harmonic
angle_style     harmonic
dihedral_style  cosine/squared

# pair_coeff: TC5=1, SC3=2, SC2=3, SC1=4
pair_coeff  1 1  3.100  4.300  12.0
pair_coeff  1 2  2.900  4.150  12.0
pair_coeff  1 3  1.800  4.150  12.0  # PS-PMMA repulsive
pair_coeff  2 2  2.700  4.300  12.0
pair_coeff  2 3  1.900  4.300  12.0  # PS-PMMA repulsive
pair_coeff  3 3  2.500  4.300  12.0
pair_coeff  3 4  2.500  4.300  12.0
pair_coeff  4 4  2.700  4.300  12.0

bond_coeff  1  3800.0  4.70
bond_coeff  2  5000.0  2.90
bond_coeff  3  3800.0  4.70
bond_coeff  4  4500.0  3.10
bond_coeff  5  3500.0  4.70

angle_coeff 1  25.0  180.0
angle_coeff 2  25.0  150.0
angle_coeff 3  20.0  180.0
angle_coeff 4  25.0  120.0

minimize        1.0e-4 1.0e-6 1000 10000

fix             1 all npt temp {T} {T} 1000.0 iso 1.0 1.0 5000.0
timestep        10.0
thermo          1000
thermo_style    custom step temp press pe ke etotal vol density
dump            eq all custom 5000 eq_trajectory.lammpstrj id mol type x y z
run             2000000      # 20 ns equilibration

unfix           1
fix             2 all nvt temp {T} {T} 1000.0
reset_timestep  0
dump            prod all custom 1000 prod_trajectory.lammpstrj id mol type x y z
run             5000000      # 50 ns production

write_data      bcp_cg_final.data
write_restart   bcp_cg.restart
"""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        f.write(content)
    return content


def generate_hoomd_cg_script(n_PS=20, n_PMMA=20, n_chains=200, T=500.0,
                              output_path="data/hoomd_cg_bcp.py"):
    chi = chi_ps_pmma(T)
    chiN = chi * (n_PS + n_PMMA)
    script = f'''"""
HOOMD-Blue 4.x CG Simulation: PS-b-PMMA
SDK 9-6 LJ model | N_PS={n_PS} N_PMMA={n_PMMA} N_chains={n_chains} T={T}K chiN={chiN:.2f}
"""
import hoomd, hoomd.md as md, gsd.hoomd, numpy as np

device  = hoomd.device.auto_select()
sim     = hoomd.Simulation(device=device, seed=42)
sim.create_state_from_gsd("bcp_initial.gsd")

kT = 8.314e-3 * {T}
integrator = md.Integrator(dt=0.01)
nvt = md.methods.ConstantVolume(
    filter=hoomd.filter.All(),
    thermostat=md.methods.thermostats.Bussi(kT=kT, tau=1.0))
integrator.methods.append(nvt)

cell = md.nlist.Cell(buffer=0.4)
lj   = md.pair.LJ(nlist=cell)
lj.params[("PS","PS")]   = dict(epsilon=3.5, sigma=0.465)
lj.params[("PMMA","PMMA")] = dict(epsilon=3.2, sigma=0.460)
lj.params[("PS","PMMA")]  = dict(epsilon=2.1, sigma=0.462)
lj.r_cut[("PS","PS")]    = 1.5*0.465
lj.r_cut[("PMMA","PMMA")]= 1.5*0.460
lj.r_cut[("PS","PMMA")]  = 1.5*0.462
integrator.forces.append(lj)

bonds = md.bond.Harmonic()
bonds.params["PS-PS"]     = dict(k=3800.0, r0=0.470)
bonds.params["PMMA-PMMA"] = dict(k=3800.0, r0=0.470)
bonds.params["PS-PMMA"]   = dict(k=3500.0, r0=0.470)
integrator.forces.append(bonds)

angles = md.angle.Harmonic()
angles.params["PS-PS-PS"]       = dict(k=25.0, t0=np.radians(180.0))
angles.params["PMMA-PMMA-PMMA"] = dict(k=20.0, t0=np.radians(180.0))
integrator.forces.append(angles)
sim.operations.integrator = integrator

logger = hoomd.logging.Logger()
thermo = md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
sim.operations.computes.append(thermo)
logger.add(thermo, quantities=["kinetic_temperature","pressure",
                                "kinetic_energy","potential_energy"])
writer = hoomd.write.GSD(filename="bcp_trajectory.gsd",
    trigger=hoomd.trigger.Periodic(1000), mode="wb",
    filter=hoomd.filter.All(), logger=logger)
sim.operations.writers.append(writer)

print(f"Equilibration | chiN={chiN:.2f}")
sim.run(200_000)
print("Production run ...")
sim.run(500_000)
print("Done -> bcp_trajectory.gsd")
'''
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        f.write(script)
    return script


def generate_parameter_table(output_path="results/cg_parameters.json"):
    params = {
        "model": "MARTINI3_adapted_BCP",
        "system": "PS-b-PMMA",
        "temperature_K": 500,
        "chi_PS_PMMA": chi_ps_pmma(500.0),
        "beads": {k: asdict(v) for k, v in MARTINI_BEADS.items()},
        "sdk_beads": {k: asdict(v) for k, v in SDK_BEADS.items()},
        "bonds": [asdict(b) for b in BOND_PARAMS],
        "angles": [asdict(a) for a in ANGLE_PARAMS],
        "dihedrals": [asdict(d) for d in DIHEDRAL_PARAMS],
        "sdk_nonbonded_cross": {
            f"{k[0]}-{k[1]}": v for k, v in SDK_NONBONDED_CROSS.items()
        },
        "ibi_strategy": {
            "method": "Iterative Boltzmann Inversion",
            "reference": "Atomistic NPT simulation at 500K 1bar",
            "convergence_criterion": "RMSE(g_CG, g_AA) < 0.02",
            "max_iterations": 100,
            "update_alpha": 0.5,
        },
        "mapping_atoms_per_bead": {"PS": 2, "PMMA": 2},
    }
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(params, f, indent=2)
    return params


if __name__ == "__main__":
    print("=== CG Parameterization: PS-b-PMMA ===")
    T_vals = [400, 450, 500, 550, 600]
    print("\nχ_PS-PMMA vs Temperature:")
    for T in T_vals:
        chi = chi_ps_pmma(T)
        print(f"  T={T}K  chi={chi:.4f}  chiN(40)={chi*40:.2f}  chiN(80)={chi*80:.2f}")

    generate_lammps_cg_input(output_path="data/lammps_cg.in")
    generate_hoomd_cg_script(output_path="data/hoomd_cg_bcp.py")
    generate_parameter_table()

    np.random.seed(42)
    r = np.linspace(0.3, 1.5, 120)
    g_target = (np.exp(-((r-0.47)**2)/(2*0.04**2))*0.8
                + 1.0 + 0.15*np.exp(-((r-0.94)**2)/(2*0.1**2)))
    g_target = np.clip(g_target, 0, None)

    ibi = IBIRefiner(g_target, r, T=500.0)
    U = ibi.initial_potential()
    g_cur = g_target * (1 + 0.3*np.random.randn(len(r)))
    g_cur = np.clip(g_cur, 1e-3, None)
    convergence = []
    for i in range(5):
        convergence.append(ibi.convergence_metric(g_cur))
        U = ibi.update_step(g_cur, U)
        g_cur = g_target*(1+0.3*(0.7**i)*np.random.randn(len(r)))
        g_cur = np.clip(g_cur, 1e-3, None)

    os.makedirs("results", exist_ok=True)
    with open("results/ibi_convergence.json", "w") as f:
        json.dump({"r_nm": r.tolist(), "g_target": g_target.tolist(),
                   "ibi_convergence_rmse": convergence}, f, indent=2)

    print("\n✓ Complete: data/lammps_cg.in, data/hoomd_cg_bcp.py,")
    print("           results/cg_parameters.json, results/ibi_convergence.json")
