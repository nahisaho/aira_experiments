"""
Dynamic Self-Assembly Simulation Protocols
Nucleation, Growth, and Defect Annealing

Models:
  1. Classical Nucleation Theory (CNT) for microphase separation onset
  2. Phase-field model for spinodal decomposition
  3. Defect density evolution (dislocations, disclinations)
  4. Annealing protocol design (thermal, solvent vapor)

References:
  - Hohenberg & Halperin, Rev.Mod.Phys 1977 (Model B)
  - Cheng et al., ACS Nano 2010 (DSA defect annealing)
  - Segalman, Matsudaira, Kramer, Macromolecules 2003
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, laplace
from scipy.signal import find_peaks
import json, os

# ─────────────────────────────────────────────────────────────────────────────
# Classical Nucleation Theory (CNT) for BCP Microphase Separation
# ─────────────────────────────────────────────────────────────────────────────

class NucleationTheory:
    """
    Microphase nucleation in the weak-segregation limit.
    Nucleation barrier based on Leibler free energy expansion.
    """
    def __init__(self, N=100, f=0.5, b=0.68e-9, rho=6e26):
        self.N   = N
        self.f   = f
        self.b   = b         # segment length (m)
        self.rho = rho       # bead density (1/m^3)
        self.kB  = 1.38e-23  # J/K
        self.chi_ODT = 10.495 / N   # MF ODT for f=0.5

    def free_energy_barrier(self, chi, T_K):
        """
        Nucleation barrier dF* (in units of kT).
        Using Binder-Frisch approximation for microphase.
        """
        kT = self.kB * T_K
        eps = (chi - self.chi_ODT) / self.chi_ODT   # reduced distance from ODT
        if eps <= 0:
            return np.inf   # below ODT: no barrier
        # Interfacial tension ~ kT * chi^0.5 / b^2
        gamma = kT * np.sqrt(chi) / (self.b**2 * self.rho**(2/3))
        # Critical nucleus size: R* ~ b * N^{1/2} / eps
        R_star = self.b * np.sqrt(self.N) / max(eps, 1e-6)
        # Barrier: dF* ~ gamma * R*^2 / kT
        dF_star = gamma * R_star**2 / kT
        return max(dF_star, 0)

    def nucleation_rate(self, chi, T_K):
        """
        Nucleation rate J ~ exp(-dF*/kT) [relative units].
        """
        dF = self.free_energy_barrier(chi, T_K)
        return np.exp(-min(dF, 700))

    def induction_time(self, chi, T_K, D0=1e-14):
        """
        Induction time tau_ind ~ R*^2 / D   (diffusion-limited)
        D0: diffusion prefactor (m^2/s), Rouse: D ~ kT/(N*xi)
        """
        eps = (chi - self.chi_ODT) / self.chi_ODT
        if eps <= 0:
            return np.inf
        R_star = self.b * np.sqrt(self.N) / max(eps, 1e-6)
        D = D0 / self.N
        return R_star**2 / D


# ─────────────────────────────────────────────────────────────────────────────
# Cahn-Hilliard / Model B Phase-Field for Spinodal Decomposition
# ─────────────────────────────────────────────────────────────────────────────

class CahnHilliardBCP:
    """
    Model B (conserved order parameter) dynamics for BCP.
    Dimensionless units: length in units of Rg, time in units of tau_D.

    dphi/dt = M * laplacian(dF/dphi)
    F = ∫ [a2*phi^2 + a4*phi^4 + kappa*(∇phi)^2] dr
    """
    def __init__(self, Nx=128, Ny=128, dx=0.5, dt=0.01,
                 a2=-0.5, a4=0.25, kappa=1.0, M=1.0,
                 f0=0.0, noise=0.05, seed=42):
        self.Nx, self.Ny = Nx, Ny
        self.dx, self.dt = dx, dt
        self.a2, self.a4 = a2, a4
        self.kappa = kappa
        self.M = M
        np.random.seed(seed)
        self.phi = f0 + noise * np.random.randn(Nx, Ny)
        self.history = []

    def mu(self):
        """Chemical potential: mu = dF/dphi"""
        phi = self.phi
        return (2*self.a2*phi + 4*self.a4*phi**3
                - self.kappa * laplace(phi) / self.dx**2)

    def step(self):
        """One Euler step of Model B."""
        mu = self.mu()
        laplace_mu = laplace(mu) / self.dx**2
        self.phi += self.dt * self.M * laplace_mu
        # No-flux (periodic BC handled by laplace with wrap mode)

    def structure_factor(self):
        """Compute S(q) = |phi_k|^2"""
        phi_k = np.fft.fft2(self.phi)
        Sq = np.abs(phi_k)**2 / (self.Nx * self.Ny)
        return np.fft.fftshift(Sq)

    def peak_wavenumber(self):
        """Return dominant q* from structure factor."""
        Sq = self.structure_factor()
        Nx, Ny = self.Nx, self.Ny
        qx = np.fft.fftshift(np.fft.fftfreq(Nx, d=self.dx)) * 2*np.pi
        qy = np.fft.fftshift(np.fft.fftfreq(Ny, d=self.dx)) * 2*np.pi
        Qx, Qy = np.meshgrid(qx, qy, indexing='ij')
        Q = np.sqrt(Qx**2 + Qy**2)
        # Azimuthal average
        q_bins = np.linspace(0, np.max(Q)/2, 60)
        Sq_avg = np.zeros(len(q_bins)-1)
        for i in range(len(q_bins)-1):
            mask = (Q >= q_bins[i]) & (Q < q_bins[i+1])
            if mask.sum() > 0:
                Sq_avg[i] = Sq[mask].mean()
        q_centers = 0.5*(q_bins[:-1]+q_bins[1:])
        idx = np.argmax(Sq_avg[2:]) + 2
        return q_centers[idx], Sq_avg

    def run(self, n_steps, record_every=100):
        """Run simulation and record snapshots."""
        for step in range(n_steps):
            self.step()
            if step % record_every == 0:
                q_star, _ = self.peak_wavenumber()
                self.history.append({
                    "step": step,
                    "phi_mean": float(self.phi.mean()),
                    "phi_std": float(self.phi.std()),
                    "q_star": float(q_star),
                    "phi_snap": self.phi.copy() if step % (record_every*10) == 0 else None
                })
        return self.history


# ─────────────────────────────────────────────────────────────────────────────
# Defect Analysis and Annealing Model
# ─────────────────────────────────────────────────────────────────────────────

class DefectDynamics:
    """
    Phenomenological model for defect density evolution during annealing.
    Defect types: dislocations (edge, topological charge ±1) and disclinations.

    Model (Gruner & Fredrickson):
    drho_def/dt = -k_ann * rho_def^2 + k_gen * exp(-E_a/kT)
    """
    def __init__(self, N=100, L0_nm=25.0):
        self.N    = N
        self.L0   = L0_nm   # nm
        self.kB   = 8.314e-3  # kJ/mol/K

    def anneal_rate(self, T_K, E_a=50.0):
        """Arrhenius annihilation rate constant."""
        return 1e12 * np.exp(-E_a / (self.kB * T_K))

    def generation_rate(self, T_K, E_gen=100.0):
        """Thermal generation rate."""
        return 1e10 * np.exp(-E_gen / (self.kB * T_K))

    def simulate_anneal(self, T_K, t_max=1e-3, dt=1e-7,
                        rho0=1e13, E_a=50.0):
        """
        Simulate defect density evolution rho(t) during isothermal anneal.
        rho in units of defects/m^2
        """
        k_ann = self.anneal_rate(T_K, E_a)
        k_gen = self.generation_rate(T_K)
        t = 0; rho = rho0
        times, densities = [0], [rho0]
        while t < t_max:
            drho = (-k_ann * rho**2 + k_gen) * dt
            rho = max(rho + drho, 0)
            t += dt
            if len(times) % 100 == 0:
                times.append(t)
                densities.append(rho)
        return np.array(times), np.array(densities)

    def equilibrium_defect_density(self, T_K, E_a=50.0):
        """rho_eq = sqrt(k_gen/k_ann)"""
        k_ann = self.anneal_rate(T_K, E_a)
        k_gen = self.generation_rate(T_K)
        return np.sqrt(k_gen / k_ann)

    def recommend_annealing_protocol(self, target_rho=1e10):
        """
        Find temperature and time for target defect density.
        Returns recommended T(t) profile.
        """
        protocols = []
        for T in [400, 420, 450, 480, 500]:
            t_vals, rho_vals = self.simulate_anneal(T, t_max=1e-2, dt=1e-7,
                                                     rho0=1e13)
            # Find time to reach target
            idx = np.searchsorted(-rho_vals, -target_rho)
            if idx < len(t_vals):
                protocols.append({"T_K":T,"t_anneal_s":t_vals[idx],
                                   "rho_final":rho_vals[-1]})
        return sorted(protocols, key=lambda x: x["t_anneal_s"])


# ─────────────────────────────────────────────────────────────────────────────
# LAMMPS Protocol for Controlled Cooling (Nucleation Study)
# ─────────────────────────────────────────────────────────────────────────────

def generate_lammps_nucleation_protocol(output_path="data/lammps_nucleation.in"):
    content = """# LAMMPS: BCP Nucleation & Growth Protocol
# Controlled cooling from disordered melt (T > T_ODT) to ordered state

units           real
atom_style      molecular
boundary        p p p
read_data       bcp_cg_initial.data

# [Force field same as equilibration run - include here]
include         forcefield.inc

# ── Step 1: Melt at T=600K (above ODT) ──────────────────
fix             melt all nvt temp 600.0 600.0 500.0
timestep        10.0
run             500000    # 5 ns melt
write_restart   melt_600K.restart

# ── Step 2: Quench to T_ODT+5K (near ODT, watch nucleation) ──
unfix           melt
variable        T_ODT equal 530.0     # estimated ODT temperature for this system
fix             quench all nvt temp 600.0 ${T_ODT} 1000.0
run             200000    # 2 ns linear cool
unfix           quench

# ── Step 3: Hold at T_ODT+5K, track S(q) for nucleation ──
fix             hold1 all nvt temp ${T_ODT} ${T_ODT} 500.0
compute         sq all structure/factor 100 1.0 sq_output.dat
fix             sq_avg all ave/time 100 10 1000 c_sq mode vector file sq_time.dat
dump            nuc all custom 2000 nucleation.lammpstrj id mol type x y z
run             2000000   # 20 ns: observe nucleation

# ── Step 4: Cool to T=500K (strong segregation) ──────────
unfix           hold1
fix             cool all nvt temp ${T_ODT} 500.0 2000.0
run             1000000   # 10 ns

# ── Step 5: Anneal at T=500K ─────────────────────────────
unfix           cool
fix             anneal all nvt temp 500.0 500.0 500.0
run             5000000   # 50 ns annealing
write_data      after_anneal.data
write_restart   annealed.restart
"""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        f.write(content)
    return content


def generate_lammps_defect_anneal(output_path="data/lammps_defect_anneal.in"):
    content = """# LAMMPS: Defect Annealing Protocol for DSA BCP
# Starts from a defective lamellar structure, anneals defects

units           real
atom_style      molecular
boundary        p p p
read_restart    defective_lam.restart

include         forcefield.inc

# ── Thermal annealing: stepwise temperature ramp ─────────
variable        T_start equal 440.0
variable        T_ann   equal 500.0
variable        T_end   equal 440.0

fix             ann1 all nvt temp ${T_start} ${T_ann} 500.0
dump            def1 all custom 5000 defect_ann1.lammpstrj id mol type x y z
run             1000000   # 10 ns ramp-up

unfix           ann1
fix             ann2 all nvt temp ${T_ann} ${T_ann} 500.0
run             5000000   # 50 ns hold (defect annihilation)
write_restart   mid_anneal.restart

# ── Solvent vapor anneal (reduced chi): lower T_glass ────
# Simulate SVA by reducing epsilon_cross temporarily
variable        eps_cross equal 1.5   # reduced from 1.9 (solvent swelling)
pair_coeff  2 3  ${eps_cross}  4.300  12.0
fix             sva all nvt temp ${T_ann} ${T_ann} 500.0
run             2000000   # 20 ns SVA

# ── Remove solvent (restore chi) and cool ────────────────
pair_coeff  2 3  1.9  4.300  12.0
fix             cool all nvt temp ${T_ann} ${T_end} 1000.0
run             2000000   # 20 ns cool
write_data      defect_annealed.data
"""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        f.write(content)
    return content


# ─────────────────────────────────────────────────────────────────────────────
# Visualization
# ─────────────────────────────────────────────────────────────────────────────

def plot_dynamics(output_path="figures/dynamics_analysis.png"):
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # 1. Phase-field simulation
    ch = CahnHilliardBCP(Nx=128, Ny=128, dx=0.5, dt=0.005,
                          a2=-0.3, a4=0.25, kappa=0.5, noise=0.05)
    history = ch.run(n_steps=5000, record_every=100)

    # Select 3 snapshots
    snaps_idx = [0, 20, 49]
    ax_labels = ["t=0 (random)", "t=1000 steps\n(nucleation)", "t=5000 steps\n(ordered)"]
    snaps = [h for h in history if h["phi_snap"] is not None]
    for i, (idx, lbl) in enumerate(zip(range(min(3,len(snaps))), ax_labels)):
        ax = axes[0][i]
        snap_idx = min(idx, len(snaps)-1)
        phi_snap = snaps[snap_idx]["phi_snap"]
        if phi_snap is not None:
            im = ax.imshow(phi_snap, cmap="RdBu_r", origin="lower",
                           vmin=-1.0, vmax=1.0,
                           extent=[0,128*0.5,0,128*0.5])
            ax.set_title(lbl, fontsize=11)
            ax.set_xlabel("x (nm)"); ax.set_ylabel("y (nm)")
            plt.colorbar(im, ax=ax, label="φ (order parameter)")

    # 2. q*(t) evolution
    ax = axes[1][0]
    steps = [h["step"] for h in history]
    q_stars = [h["q_star"] for h in history]
    ax.plot(steps, q_stars, 'b-', linewidth=2)
    ax.axhline(np.mean(q_stars[-20:]), color='red', ls='--', lw=1.5,
               label=f"q* = {np.mean(q_stars[-20:]):.3f}")
    ax.set_xlabel("Simulation Step", fontsize=11)
    ax.set_ylabel("Peak Wavenumber q* (1/nm)", fontsize=11)
    ax.set_title("Order Parameter q*(t) Evolution", fontsize=12)
    ax.legend(); ax.grid(alpha=0.3)

    # 3. Defect density annealing
    ax = axes[1][1]
    dd = DefectDynamics(N=100, L0_nm=25.0)
    for T, col in zip([440, 460, 480, 500],
                       plt.cm.YlOrRd(np.linspace(0.3, 0.9, 4))):
        t, rho = dd.simulate_anneal(T, t_max=5e-3, dt=1e-7, rho0=1e13)
        ax.semilogy(t*1e3, rho/1e12, color=col, lw=2, label=f"T={T}K")
    ax.axhline(1.0, color='blue', ls='--', lw=1.5, alpha=0.7,
               label="Target: 1/μm²")
    ax.set_xlabel("Annealing Time (ms)", fontsize=11)
    ax.set_ylabel("Defect Density (×10¹² m⁻²)", fontsize=11)
    ax.set_title("Defect Annealing Kinetics", fontsize=12)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # 4. Nucleation rate vs χ
    ax = axes[1][2]
    nt = NucleationTheory(N=100, f=0.5)
    chi_arr = np.linspace(0.08, 0.20, 100)
    T_K = 500.0
    J_arr = [nt.nucleation_rate(c, T_K) for c in chi_arr]
    tau_arr = [nt.induction_time(c, T_K) for c in chi_arr]
    ax2 = ax.twinx()
    ax.semilogy(chi_arr, J_arr, 'b-', lw=2, label="Nucleation rate J")
    ax2.semilogy(chi_arr, [min(t,1e10) for t in tau_arr],
                 'r--', lw=2, label="Induction time τ (s)")
    ax.set_xlabel("Flory-Huggins χ", fontsize=11)
    ax.set_ylabel("Nucleation Rate J (rel.)", fontsize=11, color='b')
    ax2.set_ylabel("Induction Time τ (s)", fontsize=11, color='r')
    ax.set_title("Nucleation Kinetics vs χ (N=100, T=500K)", fontsize=12)
    ax.axvline(nt.chi_ODT, color='k', ls=':', lw=1.5, label=f"χ_ODT={nt.chi_ODT:.4f}")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1+lines2, labels1+labels2, fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle("BCP Self-Assembly Dynamics: Nucleation, Growth & Defect Annealing",
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
    print("=== Dynamic Self-Assembly Simulation ===")

    # CNT analysis
    nt = NucleationTheory(N=100, f=0.5)
    print("\nNucleation analysis (N=100, T=500K):")
    for chi in [0.10, 0.12, 0.14, 0.16, 0.18, 0.20]:
        dF = nt.free_energy_barrier(chi, 500.0)
        J  = nt.nucleation_rate(chi, 500.0)
        print(f"  chi={chi:.2f}  dF*={dF:.2f}kT  J_rel={J:.3e}")

    # Defect annealing
    dd = DefectDynamics(N=100, L0_nm=25.0)
    protocols = dd.recommend_annealing_protocol(target_rho=1e10)
    print("\nRecommended annealing protocols:")
    for p in protocols[:5]:
        print(f"  T={p['T_K']}K  t≈{p['t_anneal_s']:.3e}s  rho_final={p['rho_final']:.2e}")

    os.makedirs("results", exist_ok=True)
    with open("results/annealing_protocols.json","w") as f:
        json.dump(protocols, f, indent=2, default=float)

    generate_lammps_nucleation_protocol()
    generate_lammps_defect_anneal()
    plot_dynamics()

    print("\n✓ Dynamic simulation protocols complete.")
    print("  data/lammps_nucleation.in, data/lammps_defect_anneal.in")
    print("  results/annealing_protocols.json, figures/dynamics_analysis.png")
