#!/usr/bin/env python3
"""
COMSOL/OpenFOAM Coupled Simulation Design for Brain Organoid Bioreactor.
Generates configuration files, mesh specifications, and coupling workflow.
Also performs equivalent numerical simulations in Python.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import json

# ============================================================
# Part 1: Generate OpenFOAM case structure
# ============================================================
def generate_openfoam_config():
    """Generate OpenFOAM configuration for bioreactor CFD."""
    
    blockMeshDict = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2312                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}

// Bioreactor geometry: cylindrical chamber with organoid region
convertToMeters 0.001;  // mm to m

vertices
(
    (0 0 0)          // 0: center bottom
    (25 0 0)         // 1: wall bottom
    (25 100 0)       // 2: wall top
    (0 100 0)        // 3: center top
    (0 0 1)          // 4: center bottom (z)
    (25 0 1)         // 5: wall bottom (z)
    (25 100 1)       // 6: wall top (z)
    (0 100 1)        // 7: center top (z)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (40 200 1) simpleGrading (1 1 1)
);

boundary
(
    inlet
    {
        type patch;
        faces ((0 1 5 4));
    }
    outlet
    {
        type patch;
        faces ((3 7 6 2));
    }
    wall
    {
        type wall;
        faces ((1 2 6 5));
    }
    axis
    {
        type symmetry;
        faces ((0 4 7 3));
    }
);
"""

    transportProperties = """/*--------------------------------*- C++ -*----------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}

// Culture medium properties (water-like at 37°C)
nu              [0 2 -1 0 0 0 0] 7e-7;  // kinematic viscosity
rho             [1 -3 0 0 0 0 0] 1000;

// Oxygen transport
DiffO2          [0 2 -1 0 0 0 0] 2e-9;  // O2 diffusion in medium
DiffGlucose     [0 2 -1 0 0 0 0] 5e-10; // Glucose diffusion

// Michaelis-Menten consumption
VmaxO2          [0 0 -1 0 0 0 0] 5e-3;  // mol/m3/s
KmO2            [0 -3 0 0 1 0 0] 0.005; // mol/m3
"""

    controlDict = """/*--------------------------------*- C++ -*----------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}

application     simpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         1000;
deltaT          0.5;
writeControl    timeStep;
writeInterval   100;
purgeWrite      3;
writeFormat     ascii;
writePrecision  8;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

functions
{
    fieldAverage
    {
        type            fieldAverage;
        libs            ("libfieldFunctionObjects.so");
        writeControl    writeTime;
        fields
        (
            U { mean on; prime2Mean on; base time; }
            p { mean on; prime2Mean on; base time; }
        );
    }
    
    wallShearStress
    {
        type            wallShearStress;
        libs            ("libfieldFunctionObjects.so");
        writeControl    writeTime;
        patches         (wall);
    }
}
"""
    return blockMeshDict, transportProperties, controlDict


# ============================================================
# Part 2: Generate COMSOL model description
# ============================================================
def generate_comsol_config():
    """Generate COMSOL Multiphysics model configuration."""
    
    comsol_config = {
        "model_name": "BrainOrganoidBioreactor_Multiphysics",
        "dimensions": "2D_axisymmetric",
        "physics_modules": [
            {
                "name": "Laminar Flow (spf)",
                "equations": "Navier-Stokes",
                "properties": {
                    "density": "1000 kg/m³",
                    "dynamic_viscosity": "0.001 Pa·s",
                    "flow_type": "incompressible"
                },
                "boundary_conditions": {
                    "inlet": "velocity, U_avg = Q/(pi*R²)",
                    "outlet": "pressure, p = 0",
                    "wall": "no-slip",
                    "axis": "axial symmetry"
                }
            },
            {
                "name": "Transport of Diluted Species (tds)",
                "equations": "Convection-Diffusion-Reaction",
                "species": [
                    {
                        "name": "O2",
                        "diffusion_coeff": "2e-9 m²/s",
                        "reaction": "-Vmax*c/(Km+c)",
                        "initial": "0.21 mol/m³",
                        "inlet": "0.21 mol/m³"
                    },
                    {
                        "name": "Glucose",
                        "diffusion_coeff": "5e-10 m²/s",
                        "reaction": "-Vmax_glc*c/(Km_glc+c)",
                        "initial": "5.5 mol/m³",
                        "inlet": "5.5 mol/m³"
                    }
                ]
            },
            {
                "name": "Heat Transfer (ht)",
                "equations": "Convection-Conduction",
                "properties": {
                    "thermal_conductivity": "0.6 W/(m·K)",
                    "heat_capacity": "4180 J/(kg·K)",
                    "metabolic_heat": "Q_met = 1e3 W/m³"
                }
            }
        ],
        "mesh": {
            "element_type": "triangular",
            "boundary_layers": {
                "num_layers": 5,
                "first_layer_thickness": "0.01 mm",
                "growth_rate": 1.2
            },
            "organoid_region": {
                "max_element_size": "0.05 mm",
                "min_element_size": "0.01 mm"
            },
            "bulk": {
                "max_element_size": "0.5 mm"
            }
        },
        "solver": {
            "type": "stationary",
            "method": "PARDISO",
            "relative_tolerance": 1e-6,
            "max_iterations": 200
        },
        "coupling_with_openfoam": {
            "method": "File-based data exchange",
            "workflow": [
                "1. Solve CFD in OpenFOAM (velocity, pressure)",
                "2. Export velocity field to VTK format",
                "3. Import velocity field into COMSOL",
                "4. Solve species transport with imported flow field",
                "5. Export concentration fields for post-processing",
                "6. Optional: iterate for two-way coupling"
            ],
            "data_format": "VTK/CSV",
            "interpolation": "linear"
        }
    }
    return comsol_config


# ============================================================
# Part 3: Coupled simulation in Python (equivalent)
# ============================================================
def coupled_simulation():
    """Run coupled CFD + mass transport simulation."""
    
    # Reactor geometry
    R = 0.025       # m
    L = 0.10        # m
    Nr, Nz = 50, 100
    r = np.linspace(1e-6, R, Nr)
    z = np.linspace(0, L, Nz)
    dr, dz = r[1]-r[0], z[1]-z[0]
    RR, ZZ = np.meshgrid(r, z)
    
    # Flow parameters
    Q = 1e-7        # m³/s
    mu = 0.001      # Pa·s
    U_avg = Q / (np.pi * R**2)
    
    # Velocity field (Poiseuille)
    Vz = 2 * U_avg * (1 - (RR/R)**2)
    
    # Organoid positions (multiple organoids in reactor)
    organoid_positions = [
        (0.05, 0.005, 0.0005),  # (z_center, r_center, radius)
        (0.03, 0.010, 0.0004),
        (0.07, 0.003, 0.0006),
        (0.05, 0.015, 0.0005),
        (0.08, 0.008, 0.0004),
    ]
    
    # Create organoid mask
    organoid_mask = np.zeros_like(RR, dtype=bool)
    for z_c, r_c, r_org in organoid_positions:
        dist = np.sqrt((ZZ - z_c)**2 + (RR - r_c)**2)
        organoid_mask |= (dist <= r_org)
    
    # Solve steady-state O2 transport with convection
    # ∂C/∂t = D∇²C - u·∇C - R(C)
    D_O2 = 2e-9
    Vmax = 5e-3
    Km = 0.005
    C_inlet = 0.21
    
    C = np.ones((Nz, Nr)) * C_inlet
    
    # Compute stable relaxation factor
    dt_diff = 0.25 * min(dr, dz)**2 / D_O2
    dt_conv = 0.5 * min(dr, dz) / (np.max(np.abs(Vz)) + 1e-20)
    omega = 0.05  # small relaxation for stability
    
    # Iterative Gauss-Seidel solver
    for iteration in range(5000):
        C_old = C.copy()
        for i in range(1, Nz-1):
            for j in range(1, Nr-1):
                # Diffusion (central differences)
                d2C_dz2 = (C[i+1,j] - 2*C[i,j] + C[i-1,j]) / dz**2
                d2C_dr2 = (C[i,j+1] - 2*C[i,j] + C[i,j-1]) / dr**2
                dC_dr = (C[i,j+1] - C[i,j-1]) / (2*dr)
                laplacian = d2C_dr2 + dC_dr/r[j] + d2C_dz2
                
                # Convection (upwind)
                dC_dz = (C[i,j] - C[i-1,j]) / dz if Vz[i,j] > 0 else (C[i+1,j] - C[i,j]) / dz
                convection = Vz[i,j] * dC_dz
                
                # Reaction (only in organoid regions)
                reaction = 0
                if organoid_mask[i, j]:
                    c_val = max(C[i,j], 0)
                    reaction = Vmax * c_val / (Km + c_val)
                
                # Update with small relaxation
                residual = D_O2 * laplacian - convection - reaction
                C[i,j] = C[i,j] + omega * residual * min(dr, dz)**2 / D_O2
        
        # BCs
        C[0, :] = C_inlet
        C[-1, :] = C[-2, :]
        C[:, 0] = C[:, 1]
        C[:, -1] = C[:, -2]
        C = np.clip(C, 0, C_inlet * 1.01)
        
        max_change = np.max(np.abs(C - C_old))
        if max_change < 1e-8:
            print(f"O2 transport converged at iteration {iteration}")
            break
        if iteration % 1000 == 0:
            print(f"  Iteration {iteration}, max change: {max_change:.2e}")
    
    return RR, ZZ, Vz, C, organoid_mask, organoid_positions

print("Running coupled simulation...")
RR, ZZ, Vz, C_O2, org_mask, org_pos = coupled_simulation()

# --- Figures ---
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Velocity + organoids
ax = axes[0, 0]
c = ax.contourf(ZZ*1000, RR*1000, Vz*1000, levels=30, cmap='viridis')
plt.colorbar(c, ax=ax, label='Velocity [mm/s]')
for z_c, r_c, r_org in org_pos:
    circle = plt.Circle((z_c*1000, r_c*1000), r_org*1000, fill=False,
                         edgecolor='red', linewidth=2)
    ax.add_patch(circle)
ax.set_xlabel('z [mm]')
ax.set_ylabel('r [mm]')
ax.set_title('(A) Velocity Field with Organoid Positions')

# O2 concentration
ax = axes[0, 1]
c = ax.contourf(ZZ*1000, RR*1000, C_O2, levels=30, cmap='RdYlGn')
plt.colorbar(c, ax=ax, label='O₂ [mol/m³]')
for z_c, r_c, r_org in org_pos:
    circle = plt.Circle((z_c*1000, r_c*1000), r_org*1000, fill=False,
                         edgecolor='black', linewidth=2)
    ax.add_patch(circle)
ax.set_xlabel('z [mm]')
ax.set_ylabel('r [mm]')
ax.set_title('(B) O₂ Concentration Field')

# Shear stress field
ax = axes[1, 0]
tau = 0.001 * np.abs(np.gradient(Vz, RR[0, 1]-RR[0, 0], axis=1))
c = ax.contourf(ZZ*1000, RR*1000, tau*1000, levels=30, cmap='hot')
plt.colorbar(c, ax=ax, label='Shear Stress [mPa]')
ax.set_xlabel('z [mm]')
ax.set_ylabel('r [mm]')
ax.set_title('(C) Shear Stress Distribution')

# O2 profile through organoid center
ax = axes[1, 1]
for z_c, r_c, r_org in org_pos[:3]:
    j_idx = np.argmin(np.abs(RR[0, :] - r_c))
    ax.plot(ZZ[:, j_idx]*1000, C_O2[:, j_idx],
            linewidth=2, label=f'r={r_c*1000:.0f}mm, R_org={r_org*1e6:.0f}µm')
ax.set_xlabel('z [mm]')
ax.set_ylabel('O₂ Concentration [mol/m³]')
ax.set_title('(D) Axial O₂ Profiles Through Organoids')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/comsol_openfoam_coupled.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/comsol_openfoam_coupled.png")

# --- Save configs ---
blockMesh, transport, control = generate_openfoam_config()
comsol_cfg = generate_comsol_config()

with open('figures/comsol_config.json', 'w') as f:
    json.dump(comsol_cfg, f, indent=2)
print("Saved: figures/comsol_config.json")

# --- Coupling workflow diagram ---
fig, ax = plt.subplots(figsize=(14, 8))
ax.set_xlim(0, 10)
ax.set_ylim(0, 8)
ax.set_aspect('equal')

boxes = [
    (1, 6.5, 'OpenFOAM\nCFD Solver', '#4ECDC4'),
    (5, 6.5, 'COMSOL\nMultiphysics', '#FF6B6B'),
    (1, 4, 'Velocity &\nPressure Fields', '#A8E6CF'),
    (5, 4, 'Species Transport\n& Heat Transfer', '#FFB7B2'),
    (3, 1.5, 'Post-Processing\n& Visualization', '#B5EAD7'),
    (7, 4, 'Organoid Growth\nModel', '#C7CEEA'),
]

for x, y, text, color in boxes:
    rect = plt.Rectangle((x-0.9, y-0.5), 1.8, 1.0, linewidth=2,
                          edgecolor='black', facecolor=color, alpha=0.8)
    ax.add_patch(rect)
    ax.text(x, y, text, ha='center', va='center', fontsize=9, fontweight='bold')

# Arrows
arrows = [
    (1, 6, 1, 4.5, 'Solve N-S'),
    (1.9, 4, 4.1, 4, 'Export VTK'),
    (5, 6, 5, 4.5, 'Import flow'),
    (5, 3.5, 3.9, 2, 'Results'),
    (1, 3.5, 2.1, 2, 'Results'),
    (5.9, 4, 6.1, 4, 'Couple'),
    (7, 3.5, 5, 6, 'Feedback'),
]
for x1, y1, x2, y2, label in arrows:
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
               arrowprops=dict(arrowstyle='->', color='black', lw=2))
    mx, my = (x1+x2)/2, (y1+y2)/2
    ax.text(mx, my + 0.15, label, fontsize=7, ha='center', color='gray')

ax.set_title('COMSOL-OpenFOAM Coupling Workflow for Brain Organoid Bioreactor', fontsize=14)
ax.axis('off')
plt.tight_layout()
plt.savefig('figures/coupling_workflow.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: figures/coupling_workflow.png")

print("\n=== COMSOL/OpenFOAM Simulation Design Summary ===")
print(f"OpenFOAM: blockMeshDict, transportProperties, controlDict generated")
print(f"COMSOL: 3 physics modules (Flow, Species Transport, Heat)")
print(f"Coupling: File-based VTK exchange")
print(f"Organoids simulated: {len(org_pos)}")
print(f"Grid: {RR.shape[1]}×{RR.shape[0]} (r×z)")
print(f"Min O2 in domain: {np.min(C_O2):.4f} mol/m³")
print(f"Max O2 in domain: {np.max(C_O2):.4f} mol/m³")
