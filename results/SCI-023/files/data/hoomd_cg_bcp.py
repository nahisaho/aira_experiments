"""
HOOMD-Blue 4.x CG Simulation: PS-b-PMMA
SDK 9-6 LJ model | N_PS=20 N_PMMA=20 N_chains=200 T=500.0K chiN=1.99
"""
import hoomd, hoomd.md as md, gsd.hoomd, numpy as np

device  = hoomd.device.auto_select()
sim     = hoomd.Simulation(device=device, seed=42)
sim.create_state_from_gsd("bcp_initial.gsd")

kT = 8.314e-3 * 500.0
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

print(f"Equilibration | chiN=1.99")
sim.run(200_000)
print("Production run ...")
sim.run(500_000)
print("Done -> bcp_trajectory.gsd")
