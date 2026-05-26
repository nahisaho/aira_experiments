#!/usr/bin/env python3
"""HOOMD-blue: Block Copolymer Self-Assembly (CG-DPD)"""
import hoomd, hoomd.md, gsd.hoomd, numpy as np

N, f_A, n_chains, L, kT = 32, 0.5, 500, 40.0, 1.0
N_A = int(f_A * N)

def create_initial_config(filename="init.gsd"):
    snap = gsd.hoomd.Frame()
    snap.particles.N = n_chains * N
    snap.particles.types = ['A', 'B']
    snap.configuration.box = [L, L, L, 0, 0, 0]
    positions, typeid, bonds_group = [], [], []
    idx = 0
    for chain in range(n_chains):
        pos = np.random.uniform(-L/2, L/2, 3)
        for bead in range(N):
            positions.append(pos + np.random.randn(3)*0.5)
            typeid.append(0 if bead < N_A else 1)
            if bead > 0: bonds_group.append([idx-1, idx])
            idx += 1
    snap.particles.position = np.array(positions)
    snap.particles.typeid = np.array(typeid)
    snap.bonds.N = len(bonds_group)
    snap.bonds.types = ['polymer']
    snap.bonds.typeid = np.zeros(len(bonds_group), dtype=int)
    snap.bonds.group = np.array(bonds_group)
    with gsd.hoomd.open(filename, 'w') as f: f.append(snap)
    return filename

def run_simulation():
    device = hoomd.device.auto_select()
    sim = hoomd.Simulation(device=device, seed=42)
    sim.create_state_from_gsd(create_initial_config())
    dpd = hoomd.md.pair.DPD(nlist=hoomd.md.nlist.Cell(buffer=0.4), kT=kT, default_r_cut=1.5)
    dpd.params[('A','A')] = dict(A=25.0, gamma=4.5)
    dpd.params[('B','B')] = dict(A=25.0, gamma=4.5)
    dpd.params[('A','B')] = dict(A=40.0, gamma=4.5)
    fene = hoomd.md.bond.FENEWCA(nlist=hoomd.md.nlist.Cell(buffer=0.4))
    fene.params['polymer'] = dict(k=30.0, r0=1.5, epsilon=1.0, sigma=1.0, delta=0.0)
    integrator = hoomd.md.Integrator(dt=0.005, methods=[
        hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All(),
            thermostat=hoomd.md.methods.thermostats.Bussi(kT=kT))
    ], forces=[dpd, fene])
    sim.operations.integrator = integrator
    sim.run(100_000)   # Equilibration
    sim.run(2_000_000) # Production
    for T in np.linspace(1.0, 0.5, 5):
        integrator.methods[0].thermostat.kT = T
        sim.run(100_000)
    print("Simulation complete!")

if __name__ == "__main__": run_simulation()
