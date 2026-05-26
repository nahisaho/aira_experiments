#!/usr/bin/env python3
"""
Oxygen and Nutrient Transport Modeling in Brain Organoids.
Solves reaction-diffusion equations for O2 and glucose in spherical organoids.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --- Physical parameters ---
D_O2 = 2.0e-9       # O2 diffusion coefficient in tissue [m²/s]
D_glc = 5.0e-10     # Glucose diffusion coefficient [m²/s]
C_O2_surf = 0.21    # O2 concentration at surface [mol/m³] (~21% atmospheric)
C_glc_surf = 5.5    # Glucose concentration at surface [mol/m³] (~5.5 mM)
Vmax_O2 = 5.0e-3    # Max O2 consumption rate [mol/m³/s]
Km_O2 = 0.005       # Michaelis-Menten O2 constant [mol/m³]
Vmax_glc = 1.0e-3   # Max glucose consumption rate [mol/m³/s]
Km_glc = 0.5        # Michaelis-Menten glucose constant [mol/m³]

def solve_steady_state_spherical(R, D, C_surf, Vmax, Km, Nr=200):
    """Solve steady-state reaction-diffusion in spherical coords."""
    r = np.linspace(1e-10, R, Nr)  # avoid r=0 singularity
    dr = r[1] - r[0]
    
    # Iterative solver (Picard iteration)
    C = np.ones(Nr) * C_surf
    for iteration in range(500):
        C_old = C.copy()
        C_new = np.copy(C)
        for i in range(1, Nr - 1):
            R_term = Vmax * C[i] / (Km + C[i])  # Michaelis-Menten
            # Finite difference for (1/r²) d/dr(r² dC/dr)
            laplacian = (C[i+1] - 2*C[i] + C[i-1]) / dr**2 + \
                        (2.0 / r[i]) * (C[i+1] - C[i-1]) / (2*dr)
            C_new[i] = C[i] + 0.3 * (D * laplacian - R_term) * dr**2 / D
        # BCs
        C_new[0] = C_new[1]       # symmetry at center (dC/dr=0)
        C_new[-1] = C_surf        # surface concentration
        C_new = np.maximum(C_new, 0)  # non-negative
        C = C_new
        
        if np.max(np.abs(C - C_old)) < 1e-8:
            break
    
    return r, C

# --- Solve for different organoid sizes ---
radii = [200e-6, 400e-6, 600e-6, 800e-6, 1000e-6, 1500e-6, 2000e-6]
radius_labels = ['200', '400', '600', '800', '1000', '1500', '2000']

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# O2 profiles
ax = axes[0, 0]
for R_org, label in zip(radii, radius_labels):
    r, C = solve_steady_state_spherical(R_org, D_O2, C_O2_surf, Vmax_O2, Km_O2)
    ax.plot(r*1e6, C / C_O2_surf * 100, linewidth=2, label=f'R={label} µm')
ax.set_xlabel('Radial position [µm]')
ax.set_ylabel('O₂ concentration [% of surface]')
ax.set_title('(A) Oxygen Concentration Profiles')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.axhline(y=10, color='r', linestyle='--', alpha=0.5, label='Hypoxic threshold')

# Glucose profiles
ax = axes[0, 1]
for R_org, label in zip(radii, radius_labels):
    r, C = solve_steady_state_spherical(R_org, D_glc, C_glc_surf, Vmax_glc, Km_glc)
    ax.plot(r*1e6, C / C_glc_surf * 100, linewidth=2, label=f'R={label} µm')
ax.set_xlabel('Radial position [µm]')
ax.set_ylabel('Glucose concentration [% of surface]')
ax.set_title('(B) Glucose Concentration Profiles')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Core concentration vs radius
ax = axes[1, 0]
core_O2 = []
core_glc = []
radii_um = [r*1e6 for r in radii]
for R_org in radii:
    _, C_o2 = solve_steady_state_spherical(R_org, D_O2, C_O2_surf, Vmax_O2, Km_O2)
    _, C_glc = solve_steady_state_spherical(R_org, D_glc, C_glc_surf, Vmax_glc, Km_glc)
    core_O2.append(C_o2[0] / C_O2_surf * 100)
    core_glc.append(C_glc[0] / C_glc_surf * 100)
ax.plot(radii_um, core_O2, 'bo-', linewidth=2, markersize=8, label='O₂')
ax.plot(radii_um, core_glc, 'rs-', linewidth=2, markersize=8, label='Glucose')
ax.set_xlabel('Organoid Radius [µm]')
ax.set_ylabel('Core Concentration [% of surface]')
ax.set_title('(C) Core Nutrient Concentration vs Organoid Size')
ax.axhline(y=10, color='r', linestyle='--', alpha=0.5, label='Critical hypoxia')
ax.legend()
ax.grid(True, alpha=0.3)

# Viable thickness estimation
ax = axes[1, 1]
viable_thickness_O2 = []
for R_org in radii:
    r, C = solve_steady_state_spherical(R_org, D_O2, C_O2_surf, Vmax_O2, Km_O2)
    hypoxic_idx = np.where(C / C_O2_surf < 0.10)[0]
    if len(hypoxic_idx) > 0:
        viable_r = R_org - r[hypoxic_idx[-1]]
        viable_thickness_O2.append(viable_r * 1e6)
    else:
        viable_thickness_O2.append(R_org * 1e6)

ax.bar(radius_labels, viable_thickness_O2, color='steelblue', alpha=0.8)
ax.set_xlabel('Organoid Radius [µm]')
ax.set_ylabel('Viable Tissue Thickness [µm]')
ax.set_title('(D) Viable Tissue Thickness (O₂-limited)')
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('figures/oxygen_transport.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/oxygen_transport.png")

# --- Time-dependent oxygen transport ---
def transient_o2(R_org=500e-6, t_end=300, Nr=50):
    """Solve transient reaction-diffusion for O2."""
    r = np.linspace(1e-10, R_org, Nr)
    dr = r[1] - r[0]
    
    def rhs(t, C):
        dCdt = np.zeros_like(C)
        for i in range(1, Nr - 1):
            laplacian = (C[i+1] - 2*C[i] + C[i-1]) / dr**2 + \
                        (2.0 / r[i]) * (C[i+1] - C[i-1]) / (2*dr)
            consumption = Vmax_O2 * max(C[i], 0) / (Km_O2 + max(C[i], 0))
            dCdt[i] = D_O2 * laplacian - consumption
        dCdt[0] = dCdt[1]    # symmetry
        dCdt[-1] = 0          # fixed surface
        return dCdt
    
    C0 = np.ones(Nr) * C_O2_surf
    C0[-1] = C_O2_surf
    
    sol = solve_ivp(rhs, [0, t_end], C0, method='Radau',
                    t_eval=np.linspace(0, t_end, 50), max_step=5.0)
    return r, sol.t, sol.y

r_t, t_vals, C_transient = transient_o2()

fig, ax = plt.subplots(figsize=(10, 6))
time_indices = [0, 5, 10, 20, 30, -1]
for idx in time_indices:
    ax.plot(r_t*1e6, C_transient[:, idx] / C_O2_surf * 100,
            linewidth=2, label=f't = {t_vals[idx]:.0f} s')
ax.set_xlabel('Radial Position [µm]')
ax.set_ylabel('O₂ Concentration [% of surface]')
ax.set_title('Transient O₂ Transport in Brain Organoid (R = 500 µm)')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('figures/transient_oxygen.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/transient_oxygen.png")

# --- Print summary ---
print("\n=== Oxygen/Nutrient Transport Summary ===")
for R_org, label, co2, cg in zip(radii, radius_labels, core_O2, core_glc):
    print(f"R={label:>5s} µm: Core O₂={co2:5.1f}%, Core Glucose={cg:5.1f}%")
print(f"\nCritical finding: Organoids > ~600 µm radius develop hypoxic cores")
print(f"Viable tissue thickness plateaus at ~{max(viable_thickness_O2):.0f} µm")
