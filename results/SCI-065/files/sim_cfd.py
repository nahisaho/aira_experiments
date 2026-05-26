#!/usr/bin/env python3
"""
CFD Simulation of Perfusion Bioreactor for Brain Organoid Culture.
Solves 2D axisymmetric Navier-Stokes (simplified) for velocity/pressure fields
and computes wall shear stress on organoid surfaces.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve

# --- Bioreactor geometry parameters ---
R_reactor = 0.025    # reactor radius [m] (25 mm)
L_reactor = 0.10     # reactor length [m] (100 mm)
R_organoid = 0.0005  # organoid radius [m] (500 µm)
mu = 0.001           # dynamic viscosity [Pa·s] (water-like medium)
rho = 1000.0         # density [kg/m³]
Q_inlet = 1e-7       # volumetric flow rate [m³/s] (0.1 mL/s)

# Grid
Nr, Nz = 80, 200
r = np.linspace(0, R_reactor, Nr)
z = np.linspace(0, L_reactor, Nz)
dr = r[1] - r[0]
dz = z[1] - z[0]
R, Z = np.meshgrid(r, z)

# Analytical Poiseuille flow profile for cylindrical pipe
U_avg = Q_inlet / (np.pi * R_reactor**2)
Re = rho * U_avg * 2 * R_reactor / mu
print(f"Reynolds number: {Re:.2f}")
print(f"Average velocity: {U_avg*1000:.3f} mm/s")

# Velocity profile: u_z(r) = 2*U_avg*(1 - (r/R)^2)
Vz = 2 * U_avg * (1 - (R / R_reactor)**2)
Vr = np.zeros_like(Vz)

# Pressure drop (Hagen-Poiseuille)
dP_dz = -8 * mu * Q_inlet / (np.pi * R_reactor**4)
P = np.zeros_like(Z)
for i in range(Nz):
    P[i, :] = -dP_dz * (L_reactor - z[i])

# Compute shear stress field: tau = mu * dVz/dr
tau = np.zeros_like(Vz)
for i in range(Nz):
    tau[i, 1:-1] = mu * (Vz[i, 2:] - Vz[i, :-2]) / (2 * dr)
    tau[i, 0] = mu * (Vz[i, 1] - Vz[i, 0]) / dr
    tau[i, -1] = mu * (Vz[i, -1] - Vz[i, -2]) / dr

# --- Organoid-surface shear stress analysis ---
# Place organoids at different radial positions
r_positions = np.array([0.0, 0.005, 0.010, 0.015, 0.020])
flow_rates = np.array([0.05e-7, 0.1e-7, 0.5e-7, 1e-7, 5e-7])

shear_at_organoid = np.zeros((len(flow_rates), len(r_positions)))
for i, Q in enumerate(flow_rates):
    U = Q / (np.pi * R_reactor**2)
    for j, rp in enumerate(r_positions):
        # Shear rate at radial position rp
        gamma_dot = 4 * U * rp / R_reactor**2
        shear_at_organoid[i, j] = mu * gamma_dot

# --- Figure 1: Velocity field ---
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Velocity magnitude
ax = axes[0, 0]
Vmag = np.sqrt(Vz**2 + Vr**2) * 1000  # mm/s
c = ax.contourf(Z*1000, R*1000, Vmag, levels=30, cmap='viridis')
plt.colorbar(c, ax=ax, label='Velocity [mm/s]')
ax.set_xlabel('Axial position z [mm]')
ax.set_ylabel('Radial position r [mm]')
ax.set_title('(A) Velocity Magnitude in Perfusion Bioreactor')

# Pressure field
ax = axes[0, 1]
c = ax.contourf(Z*1000, R*1000, P, levels=30, cmap='coolwarm')
plt.colorbar(c, ax=ax, label='Pressure [Pa]')
ax.set_xlabel('Axial position z [mm]')
ax.set_ylabel('Radial position r [mm]')
ax.set_title('(B) Pressure Distribution')

# Shear stress field
ax = axes[1, 0]
c = ax.contourf(Z*1000, R*1000, np.abs(tau)*1000, levels=30, cmap='hot')
plt.colorbar(c, ax=ax, label='Shear Stress [mPa]')
ax.set_xlabel('Axial position z [mm]')
ax.set_ylabel('Radial position r [mm]')
ax.set_title('(C) Wall Shear Stress Distribution')

# Velocity profile at mid-section
ax = axes[1, 1]
mid_idx = Nz // 2
ax.plot(r*1000, Vz[mid_idx, :]*1000, 'b-', linewidth=2)
ax.set_xlabel('Radial position r [mm]')
ax.set_ylabel('Axial velocity [mm/s]')
ax.set_title('(D) Velocity Profile at z = L/2')
ax.axhline(y=U_avg*1000, color='r', linestyle='--', label=f'U_avg = {U_avg*1000:.3f} mm/s')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/cfd_velocity_field.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/cfd_velocity_field.png")

# --- Figure 2: Shear stress on organoids ---
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
for i, Q in enumerate(flow_rates):
    label = f'Q = {Q*1e7:.2f} × 10⁻⁷ m³/s'
    ax.plot(r_positions*1000, shear_at_organoid[i, :]*1000, 'o-', label=label, linewidth=2)
ax.axhspan(0, 100, alpha=0.15, color='green', label='Safe zone (<0.1 Pa)')
ax.axhspan(100, 500, alpha=0.15, color='yellow', label='Caution (0.1-0.5 Pa)')
ax.set_xlabel('Radial Position of Organoid [mm]')
ax.set_ylabel('Shear Stress [mPa]')
ax.set_title('(A) Shear Stress on Organoid Surfaces')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Reynolds number vs flow rate
ax = axes[1]
Re_vals = [rho * (Q / (np.pi * R_reactor**2)) * 2 * R_reactor / mu for Q in flow_rates]
ax.semilogy(flow_rates * 1e7, Re_vals, 'rs-', linewidth=2, markersize=8)
ax.set_xlabel('Flow Rate [× 10⁻⁷ m³/s]')
ax.set_ylabel('Reynolds Number')
ax.set_title('(B) Reynolds Number vs Flow Rate')
ax.axhline(y=2300, color='k', linestyle='--', label='Turbulent transition (Re=2300)')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/cfd_shear_analysis.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/cfd_shear_analysis.png")

# --- Print summary ---
print("\n=== CFD Simulation Summary ===")
print(f"Reactor: R={R_reactor*1000:.1f} mm, L={L_reactor*1000:.1f} mm")
print(f"Flow rate: Q={Q_inlet*1e6:.3f} mL/s")
print(f"Re = {Re:.2f} (laminar)")
print(f"Pressure drop: {abs(dP_dz)*L_reactor:.3f} Pa")
print(f"Max shear stress: {np.max(np.abs(tau))*1000:.3f} mPa")
print(f"Max velocity: {np.max(Vmag):.3f} mm/s")
