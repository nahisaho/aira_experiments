"""
Directed Self-Assembly (DSA) Simulation
Template-Polymer Interactions for BCP Lithography

Topics:
  1. Chemical epitaxy (chemoepitaxy): stripe-patterned substrates
  2. Graphoepitaxy: topographic guide patterns
  3. Template-polymer interaction potentials
  4. Confinement effects on lamella orientation and period
  5. HOOMD-Blue DSA simulation protocol

References:
  - Kim et al., Science 2010 (chemoepitaxy DSA)
  - Segalman et al., Adv. Mater. 2001 (graphoepitaxy)
  - Liu et al., Macromolecules 2013 (simulation of DSA)
  - Bates et al., Science 2012 (BCP lithography review)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.optimize import minimize
import json, os

# ─────────────────────────────────────────────────────────────────────────────
# Chemical Epitaxy (Chemoepitaxy) Model
# ─────────────────────────────────────────────────────────────────────────────

class ChemoepitaxyModel:
    """
    Substrate with periodic chemical stripes: PS-philic / PMMA-philic alternating.
    Interaction: U_sub(x,y) = -A * cos(2π x / L_guide) * exp(-z / lambda_s)
    
    Parameters:
      L_guide: guide pattern pitch (nm)
      A: surface affinity strength (kJ/mol)
      lambda_s: surface interaction decay length (nm)
      L0: BCP natural period (nm)
      n: integer multiplier (L_guide = n * L0 for n:1 DSA)
    """
    def __init__(self, L0=25.0, n_mult=1, A=2.5, lambda_s=1.0):
        self.L0       = L0
        self.n_mult   = n_mult
        self.L_guide  = L0 * n_mult
        self.A        = A
        self.lambda_s = lambda_s

    def surface_potential(self, x, z=0.0):
        """Substrate interaction energy for PS bead at (x, z)."""
        decay = np.exp(-z / self.lambda_s)
        return -self.A * np.cos(2*np.pi * x / self.L_guide) * decay

    def commensurability_energy(self, L_bcp):
        """
        Free energy penalty for L_bcp ≠ L_guide / n_mult.
        Elastic model: dF ~ B * (L_bcp/L0 - 1)^2
        """
        B = 10.0  # kJ/mol  (elastic modulus of BCP film)
        return B * (L_bcp/self.L0 - 1.0)**2

    def optimal_multiplier(self, L_guide_range=(10, 60)):
        """
        Find best integer multiplier n such that L_guide ≈ n * L0.
        """
        n_values = range(1, 8)
        results = []
        for n in n_values:
            L_g = self.L0 * n
            if L_guide_range[0] <= L_g <= L_guide_range[1]:
                results.append({"n":n, "L_guide_nm":L_g,
                                "L0_nm":self.L0, "ratio":n})
        return results

    def order_parameter(self, phi_profile, period):
        """
        Orientational order parameter Psi2 from density profile.
        """
        x = np.linspace(0, period, len(phi_profile))
        psi2 = 2.0 * np.abs(np.mean(phi_profile * np.exp(2j*np.pi*x/period)))
        return float(np.abs(psi2))


class GraphoepitaxyModel:
    """
    Topographic guide: trench confinement forces BCP lamellar orientation.
    Trench width W, depth D, wall chemistry (neutral or preferential).

    Commensurability: W = n * L0  (integer n for well-ordered DSA)
    """
    def __init__(self, L0=25.0, W=None, D=30.0, wall_chi_S=-0.2):
        self.L0 = L0
        self.W  = W if W else L0  # default: 1:1
        self.D  = D   # trench depth (nm)
        self.wall_chi = wall_chi_S  # surface affinity (negative = PS-philic)

    def confinement_free_energy(self, n_lam):
        """
        Free energy of n lamellar periods in trench of width W.
        F = F_bulk(n) + F_surface + F_elastic
        """
        L_n = self.W / n_lam      # actual lamellar period under confinement
        F_elastic = 0.5 * (L_n - self.L0)**2 / self.L0  # kJ/mol per chain
        F_surface = abs(self.wall_chi) * 2 * self.D / self.L0
        return F_elastic + F_surface

    def optimal_lamellar_count(self, max_n=6):
        """Find integer n_lam that minimizes free energy."""
        n_arr = range(1, max_n+1)
        F_arr = [self.confinement_free_energy(n) for n in n_arr]
        n_opt = n_arr[np.argmin(F_arr)]
        return n_opt, min(F_arr), dict(zip(n_arr, F_arr))

    def critical_width_ratio(self):
        """Ratio W/L0 for optimal commensurability."""
        n_opt, _, F_dict = self.optimal_lamellar_count()
        return self.W / self.L0, n_opt

    def dsa_success_window(self, tolerance=0.10):
        """
        Range of W/L0 ratios that give acceptable order (within tolerance of integer).
        """
        windows = []
        for n in range(1, 7):
            center = n * self.L0
            window = (center*(1-tolerance), center*(1+tolerance))
            windows.append({"n":n, "W_center_nm":center,
                           "W_min_nm":window[0], "W_max_nm":window[1],
                           "tolerance":tolerance})
        return windows


# ─────────────────────────────────────────────────────────────────────────────
# HOOMD-Blue DSA Simulation Script
# ─────────────────────────────────────────────────────────────────────────────

def generate_hoomd_dsa_script(
    L0=25.0, n_mult=4, W_trench=None,
    output_path="data/hoomd_dsa.py"
) -> str:
    if W_trench is None:
        W_trench = L0 * n_mult  # commensurable width

    script = f'''"""
HOOMD-Blue 4.x: Directed Self-Assembly (DSA) Simulation
Chemoepitaxy + Graphoepitaxy for PS-b-PMMA on patterned substrate.

L0 = {L0} nm, n_mult = {n_mult}, L_guide = {L0*n_mult} nm
Trench width = {W_trench} nm
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
L_guide    = {L0*n_mult:.2f}   # nm (guide pitch)
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
wall_geometry = md.external.wall.Sphere(radius={W_trench/2:.2f})
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
print(f"DSA simulation: L0={L0}nm L_guide={L0*n_mult:.1f}nm W={W_trench:.1f}nm")
print("Phase 1: High-T equilibration (above ODT) ...")
sim.run(100_000)

print("Phase 2: Guided assembly at T=500K ...")
sim.run(500_000)

print("Phase 3: Defect annealing ...")
sim.run(200_000)

print("DSA simulation complete -> dsa_trajectory.gsd")
'''
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        f.write(script)
    return script


# ─────────────────────────────────────────────────────────────────────────────
# DSA Pattern Fidelity Analysis
# ─────────────────────────────────────────────────────────────────────────────

def compute_dsa_fidelity(L0=25.0, sigma_LWR=0.5, n_mult_range=range(1,7)):
    """
    Estimate Line Width Roughness (LWR) and placement error for DSA.
    sigma_LWR: intrinsic BCP fluctuation (nm)
    
    Pattern fidelity metrics:
      - LWR (3σ): line width roughness
      - LER: line edge roughness
      - Placement error: deviation from template center
    """
    results = []
    for n in n_mult_range:
        L_guide = L0 * n
        # LWR amplified by n (more degrees of freedom)
        lwr_3sigma = sigma_LWR * np.sqrt(n) * 3
        # Registration error ~ L_guide * 0.02 (2% of pitch)
        placement_err = L_guide * 0.02
        # Defect density (from Cheng et al. empirical)
        rho_def = 1e-3 * np.exp(0.5 * n)  # /100 nm^2
        results.append({
            "n_mult": n,
            "L_guide_nm": L_guide,
            "LWR_3sigma_nm": lwr_3sigma,
            "LER_3sigma_nm": lwr_3sigma * 0.7,
            "placement_error_nm": placement_err,
            "defect_density_per100nm2": rho_def,
            "suitable_for_7nm": lwr_3sigma < 2.0 and placement_err < 1.5,
        })
    return results


def plot_dsa_analysis(L0=25.0, output_path="figures/dsa_analysis.png"):
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # 1. Chemoepitaxy surface potential
    ax = axes[0][0]
    x = np.linspace(0, 100, 500)
    for n, col in zip([1,2,4], ["#2E86C1","#E74C3C","#28B463"]):
        model = ChemoepitaxyModel(L0=L0, n_mult=n, A=2.5)
        U = model.surface_potential(x, z=0)
        ax.plot(x, U, color=col, lw=2, label=f"n={n}, L_guide={L0*n:.0f}nm")
    ax.set_xlabel("Position x (nm)", fontsize=11)
    ax.set_ylabel("Surface Potential U (kJ/mol)", fontsize=11)
    ax.set_title("Chemoepitaxy Surface Potential", fontsize=12)
    ax.legend(fontsize=9); ax.grid(alpha=0.3)

    # 2. Graphoepitaxy commensurability window
    ax = axes[0][1]
    W_arr = np.linspace(5, 150, 1000)
    geo = GraphoepitaxyModel(L0=L0, D=30.0)
    F_arr = []
    for W in W_arr:
        geo.W = W
        F_vals = [geo.confinement_free_energy(n) for n in range(1,8)]
        F_arr.append(min(F_vals))
    ax.plot(W_arr/L0, F_arr, 'b-', lw=2)
    for n in range(1, 7):
        ax.axvline(n, color='red', ls=':', lw=1.0, alpha=0.5)
        ax.text(n, max(F_arr)*0.9, f"n={n}", ha='center', fontsize=9, color='red')
    ax.set_xlabel("Trench Width  W / L₀", fontsize=11)
    ax.set_ylabel("Confinement Free Energy (kJ/mol)", fontsize=11)
    ax.set_title("Graphoepitaxy Commensurability", fontsize=12)
    ax.grid(alpha=0.3)

    # 3. DSA success windows
    ax = axes[1][0]
    geo.W = L0
    windows = geo.dsa_success_window(tolerance=0.10)
    for i, w in enumerate(windows):
        rect = patches.Rectangle(
            (w["W_min_nm"], i*0.8), w["W_max_nm"]-w["W_min_nm"], 0.6,
            linewidth=1, edgecolor='blue', facecolor='lightblue', alpha=0.7)
        ax.add_patch(rect)
        ax.text((w["W_min_nm"]+w["W_max_nm"])/2, i*0.8+0.3,
                f"n={w['n']}", ha='center', va='center', fontsize=9)
    ax.set_xlim(0, 160); ax.set_ylim(-0.5, len(windows)*0.8)
    ax.set_xlabel("Trench Width W (nm)", fontsize=11)
    ax.set_yticks([])
    ax.set_title(f"DSA Success Windows (L₀={L0}nm, ±10%)", fontsize=12)
    ax.axvline(L0, color='red', ls='--', lw=1.5, label=f"L₀={L0}nm")
    ax.legend(); ax.grid(axis='x', alpha=0.3)

    # 4. Pattern fidelity vs n_mult
    ax = axes[1][1]
    fidelity = compute_dsa_fidelity(L0=L0, sigma_LWR=0.5)
    n_arr  = [r["n_mult"] for r in fidelity]
    lwr    = [r["LWR_3sigma_nm"] for r in fidelity]
    place  = [r["placement_error_nm"] for r in fidelity]
    rho    = [r["defect_density_per100nm2"] for r in fidelity]
    ax2    = ax.twinx()
    ax.bar([n-0.2 for n in n_arr], lwr,   0.35, color='#3498DB', alpha=0.8,
           label="LWR 3σ (nm)")
    ax.bar([n+0.2 for n in n_arr], place, 0.35, color='#E74C3C', alpha=0.8,
           label="Placement error (nm)")
    ax2.plot(n_arr, rho, 'ko-', lw=2, ms=6, label="Defect density (/100nm²)")
    ax.axhline(2.0, color='purple', ls='--', lw=1.5, label="LWR limit 2nm")
    ax.set_xlabel("Multiplication Factor  n", fontsize=11)
    ax.set_ylabel("Roughness / Error (nm)", fontsize=11)
    ax2.set_ylabel("Defect Density (/100nm²)", fontsize=11)
    ax.set_title("DSA Pattern Fidelity vs Multiplication", fontsize=12)
    lines1, lab1 = ax.get_legend_handles_labels()
    lines2, lab2 = ax2.get_legend_handles_labels()
    ax.legend(lines1+lines2, lab1+lab2, fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

    fig.suptitle("Directed Self-Assembly (DSA) Analysis\nChemoepitaxy & Graphoepitaxy",
                 fontsize=14, fontweight='bold')
    fig.tight_layout()
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=== DSA Template-Polymer Interaction Analysis ===")
    L0 = 25.0  # nm (14nm node target)

    # Chemoepitaxy analysis
    print("\nChemoepitaxy multipliers:")
    chem = ChemoepitaxyModel(L0=L0, n_mult=4)
    for r in chem.optimal_multiplier():
        print(f"  n={r['n']:d}  L_guide={r['L_guide_nm']:.1f}nm")

    # Graphoepitaxy analysis
    print("\nGraphoepitaxy commensurability:")
    geo = GraphoepitaxyModel(L0=L0, D=30.0)
    windows = geo.dsa_success_window(tolerance=0.08)
    for w in windows[:5]:
        print(f"  n={w['n']}  W={w['W_center_nm']:.1f}nm  "
              f"[{w['W_min_nm']:.1f}, {w['W_max_nm']:.1f}]nm")

    # Pattern fidelity
    print("\nDSA Pattern Fidelity:")
    fidelity = compute_dsa_fidelity(L0=L0, sigma_LWR=0.5)
    for r in fidelity:
        flag = "✓" if r["suitable_for_7nm"] else "✗"
        print(f"  n={r['n_mult']}  LWR={r['LWR_3sigma_nm']:.2f}nm  "
              f"err={r['placement_error_nm']:.2f}nm  {flag}")

    os.makedirs("results", exist_ok=True)
    with open("results/dsa_fidelity.json", "w") as f:
        json.dump(fidelity, f, indent=2, default=lambda x: bool(x) if isinstance(x, (bool, np.bool_)) else str(x))
    with open("results/dsa_commensurability.json", "w") as f:
        json.dump(windows, f, indent=2)

    generate_hoomd_dsa_script(L0=L0, n_mult=4)
    plot_dsa_analysis(L0=L0)

    print("\n✓ DSA analysis complete.")
    print("  data/hoomd_dsa.py, results/dsa_fidelity.json")
    print("  results/dsa_commensurability.json, figures/dsa_analysis.png")
