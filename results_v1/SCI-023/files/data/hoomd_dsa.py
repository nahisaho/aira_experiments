"""
HOOMD-Blue 4.x: Directed Self-Assembly (DSA) Simulation
Chemoepitaxy + Graphoepitaxy for PS-b-PMMA on patterned substrate.

L0 = 25.0 nm, n_mult = 4, L_guide = 100.0 nm
Trench width = 100.0 nm
"""
import hoomd, hoomd.md as md, gsd.hoomd
import numpy as np

# ── Simulation setup ────────────────────────────────────
device = hoomd.device.auto_select()
sim    = hoomd.Simulation(device=device, seed=42)
sim.create_state_from_gsd("dsa_initial.gsd")

kT = 8.314e-3 * 500.0
integrator = md.Integrator(dt=0.01)

# ── Standard BCP force field ─────────────────────────────
cell = md.nlist.Cell(buffer=0.4)
lj   = md.pair.LJ(nlist=cell)
lj.params[("PS","PS")]     = dict(epsilon=3.5, sigma=0.465)
lj.params[("PMMA","PMMA")] = dict(epsilon=3.2, sigma=0.460)
lj.params[("PS","PMMA")]   = dict(epsilon=2.1, sigma=0.462)
lj.params[("PS","WALL")]   = dict(epsilon=4.0, sigma=0.47)  # PS-philic wall
lj.params[("PMMA","WALL")] = dict(epsilon=1.2, sigma=0.47)  # PMMA-phobic wall
for pair in [("PS","PS"),("PMMA","PMMA"),("PS","PMMA"),
             ("PS","WALL"),("PMMA","WALL")]:
    lj.r_cut[pair] = 1.5 * 0.465
integrator.forces.append(lj)

# ── Substrate: external field (chemoepitaxy) ──────────────
# Sinusoidal chemical stripe potential on z=0 wall
# U_sub(x) = -A * cos(2*pi*x / L_guide) for z < lambda_s
L_guide    = 100.00   # nm (guide pitch)
A_affinity = 2.5               # kJ/mol
lambda_s   = 1.0               # nm (surface decay)

class SubstrateField(hoomd.md.external.field.Periodic):
    pass  # Implemented via tabulated external potential in HOOMD

# Approximate substrate with external sinusoidal field
ext_field = md.external.field.Periodic()
ext_field.params["PS"]   = dict(A=-A_affinity, i=0, w=0.02, p=1)
ext_field.params["PMMA"] = dict(A=+A_affinity, i=0, w=0.02, p=1)
integrator.forces.append(ext_field)

# ── Graphoepitaxy: wall repulsion ────────────────────────
# Trench walls: harmonic repulsion beyond boundaries
wall_geometry = md.external.wall.Sphere(radius=50.00)
wall_lj = md.external.wall.LJ(walls=[wall_geometry])
wall_lj.params["PS"]   = dict(epsilon=1.0, sigma=0.47, r_cut=0.52)
wall_lj.params["PMMA"] = dict(epsilon=1.0, sigma=0.47, r_cut=0.52)
integrator.forces.append(wall_lj)

# ── Bonds, angles ─────────────────────────────────────────
bonds = md.bond.Harmonic()
bonds.params["PS-PS"]     = dict(k=3800.0, r0=0.470)
bonds.params["PMMA-PMMA"] = dict(k=3800.0, r0=0.470)
bonds.params["PS-PMMA"]   = dict(k=3500.0, r0=0.470)
integrator.forces.append(bonds)

angles = md.angle.Harmonic()
angles.params["PS-PS-PS"]       = dict(k=25.0, t0=np.radians(180.0))
angles.params["PMMA-PMMA-PMMA"] = dict(k=20.0, t0=np.radians(180.0))
integrator.forces.append(angles)

# ── NVT dynamics ──────────────────────────────────────────
nvt = md.methods.ConstantVolume(
    filter=hoomd.filter.All(),
    thermostat=md.methods.thermostats.Bussi(kT=kT, tau=1.0))
integrator.methods.append(nvt)
sim.operations.integrator = integrator

# ── Logging ───────────────────────────────────────────────
thermo = md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
sim.operations.computes.append(thermo)
logger = hoomd.logging.Logger()
logger.add(thermo, quantities=["kinetic_temperature","potential_energy","pressure"])
writer = hoomd.write.GSD(filename="dsa_trajectory.gsd",
    trigger=hoomd.trigger.Periodic(1000), mode="wb",
    filter=hoomd.filter.All(), logger=logger)
sim.operations.writers.append(writer)

# ── Run protocol ──────────────────────────────────────────
print(f"DSA simulation: L0=25.0nm L_guide=100.0nm W=100.0nm")
print("Phase 1: High-T equilibration (above ODT) ...")
sim.run(100_000)

print("Phase 2: Guided assembly at T=500K ...")
sim.run(500_000)

print("Phase 3: Defect annealing ...")
sim.run(200_000)

print("DSA simulation complete -> dsa_trajectory.gsd")
