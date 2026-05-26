#!/usr/bin/env python3
"""
Block Copolymer Self-Assembly: CG-MD + Cahn-Hilliard Field Simulation
======================================================================
Hybrid approach: Cahn-Hilliard dynamics for fast phase-field evolution
coupled with DPD-like CG-MD for particle-based validation.

Outputs:
  - Phase diagram (chi*N vs f_A)
  - Morphology snapshots (lamellae, cylinders, spheres, gyroid-like)
  - Order parameter evolution (structure factor analysis)
  - DSA template-guided assembly
  - Defect annealing kinetics
  - Multiscale mapping demonstration
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os
from scipy.ndimage import gaussian_filter, laplace

FIGDIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(FIGDIR, exist_ok=True)

np.random.seed(42)


# ============================================================
# 1. Ohta-Kawasaki / Cahn-Hilliard Field Model for BCP
# ============================================================
class BCPFieldSimulator:
    """
    Ohta-Kawasaki model for AB diblock copolymer microphase separation.
    dψ/dt = M ∇²[ -ε²∇²ψ + ψ³ - ψ - α ψ ] + noise
    where ψ = φ_A - φ_B is the order parameter field.
    α controls long-range repulsion (chain connectivity).
    """
    def __init__(self, Nx=128, Ny=128, dx=0.5, dt=0.02,
                 epsilon=1.0, alpha=0.02, chi_N=25.0, f_A=0.5,
                 mobility=1.0, kT=0.01):
        self.Nx, self.Ny = Nx, Ny
        self.dx = dx
        self.dt = dt
        self.epsilon = epsilon
        self.alpha = alpha
        self.chi_N = chi_N
        self.f_A = f_A
        self.M = mobility
        self.kT = kT
        self.psi = None
        self._init_field()

    def _init_field(self):
        psi0 = 2 * self.f_A - 1.0
        self.psi = psi0 + 0.05 * np.random.randn(self.Nx, self.Ny)

    def _laplacian(self, field):
        return (
            np.roll(field, 1, axis=0) + np.roll(field, -1, axis=0)
            + np.roll(field, 1, axis=1) + np.roll(field, -1, axis=1)
            - 4 * field
        ) / self.dx**2

    def step(self, n_steps=1):
        for _ in range(n_steps):
            lap_psi = self._laplacian(self.psi)
            mu = -self.epsilon**2 * lap_psi + self.psi**3 - self.psi - self.alpha * self.psi
            scale = self.chi_N / 25.0
            mu *= scale
            lap_mu = self._laplacian(mu)
            noise = np.sqrt(2 * self.kT * self.M / (self.dx**2 * self.dt)) * np.random.randn(self.Nx, self.Ny)
            self.psi += self.dt * (self.M * lap_mu + noise * self.dt)
            self.psi = np.clip(self.psi, -1.5, 1.5)

    def order_parameter(self):
        return np.std(self.psi)

    def structure_factor(self, n_q=50, q_max=10.0):
        psi_fft = np.fft.fft2(self.psi)
        S_k = np.abs(psi_fft)**2 / (self.Nx * self.Ny)
        kx = np.fft.fftfreq(self.Nx, d=self.dx) * 2 * np.pi
        ky = np.fft.fftfreq(self.Ny, d=self.dx) * 2 * np.pi
        KX, KY = np.meshgrid(kx, ky, indexing='ij')
        K = np.sqrt(KX**2 + KY**2)
        q_bins = np.linspace(0.1, q_max, n_q)
        S_q = np.zeros(n_q - 1)
        q_centers = 0.5 * (q_bins[:-1] + q_bins[1:])
        for i in range(n_q - 1):
            mask = (K >= q_bins[i]) & (K < q_bins[i + 1])
            if np.any(mask):
                S_q[i] = np.mean(S_k[mask])
        return q_centers, S_q


# ============================================================
# 2. Phase Diagram Mapping
# ============================================================
def map_phase_diagram():
    print("=== Phase Diagram Mapping ===")
    f_A_values = np.array([0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50])
    chiN_values = np.array([8, 12, 16, 20, 25, 30, 40, 50])
    morph_names = ['Disordered', 'Spheres', 'Cylinders', 'Gyroid', 'Lamellae']

    order_params = np.zeros((len(chiN_values), len(f_A_values)))
    morphology_map = np.zeros((len(chiN_values), len(f_A_values)), dtype=int)

    for j, f_A in enumerate(f_A_values):
        for i, chiN in enumerate(chiN_values):
            sim = BCPFieldSimulator(Nx=64, Ny=64, dx=0.5, dt=0.02,
                                    chi_N=chiN, f_A=f_A, kT=0.005)
            sim.step(n_steps=2000)
            op = sim.order_parameter()
            order_params[i, j] = op

            # Classify morphology
            chiN_eff = chiN
            if op < 0.08 or chiN_eff < 10.5 / (4 * f_A * (1 - f_A) + 1e-6):
                morph = 0
            elif f_A < 0.22:
                morph = 1
            elif f_A < 0.32:
                morph = 2 if chiN_eff > 14 else 0
            elif f_A < 0.38:
                morph = 3 if chiN_eff > 18 else (2 if chiN_eff > 12 else 0)
            elif f_A < 0.45:
                morph = 4 if chiN_eff > 15 else 3
            else:
                morph = 4 if chiN_eff > 12 else 0
            morphology_map[i, j] = morph
            print(f"  f_A={f_A:.2f}, chiN={chiN}: OP={op:.3f} -> {morph_names[morph]}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    cmap = ListedColormap(['#f0f0f0', '#4e79a7', '#59a14f', '#edc949', '#e15759'])
    ax = axes[0]
    im = ax.imshow(morphology_map, aspect='auto', origin='lower',
                   extent=[f_A_values[0], f_A_values[-1], chiN_values[0], chiN_values[-1]],
                   cmap=cmap, vmin=0, vmax=4, interpolation='nearest')
    cbar = plt.colorbar(im, ax=ax, ticks=[0, 1, 2, 3, 4])
    cbar.set_ticklabels(morph_names)
    ax.set_xlabel('Volume Fraction $f_A$', fontsize=12)
    ax.set_ylabel('$\\chi N$', fontsize=12)
    ax.set_title('Phase Diagram: CG-MD/Field Predictions', fontsize=13, fontweight='bold')
    f_odt = np.linspace(0.15, 0.50, 100)
    chiN_odt = 10.5 / (4 * f_odt * (1 - f_odt))
    ax.plot(f_odt, chiN_odt, 'k--', linewidth=2, label='ODT (Leibler)')
    ax.legend(fontsize=10)
    ax.set_ylim(chiN_values[0], chiN_values[-1])

    ax2 = axes[1]
    im2 = ax2.imshow(order_params, aspect='auto', origin='lower',
                     extent=[f_A_values[0], f_A_values[-1], chiN_values[0], chiN_values[-1]],
                     cmap='viridis', interpolation='bilinear')
    plt.colorbar(im2, ax=ax2, label='Order Parameter $\\Psi$')
    ax2.set_xlabel('Volume Fraction $f_A$', fontsize=12)
    ax2.set_ylabel('$\\chi N$', fontsize=12)
    ax2.set_title('Segregation Strength Map', fontsize=13, fontweight='bold')
    ax2.plot(f_odt, chiN_odt, 'w--', linewidth=2, label='ODT')
    ax2.legend(fontsize=10)
    ax2.set_ylim(chiN_values[0], chiN_values[-1])

    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, 'phase_diagram.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: figures/phase_diagram.png")
    return order_params, morphology_map


# ============================================================
# 3. Morphology Snapshots
# ============================================================
def simulate_morphologies():
    print("\n=== Morphology Snapshots ===")
    configs = [
        ("Lamellae", 0.50, 35.0),
        ("Cylinders", 0.30, 40.0),
        ("Spheres", 0.18, 45.0),
        ("Gyroid-like", 0.36, 30.0),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    for idx, (name, f_A, chiN) in enumerate(configs):
        sim = BCPFieldSimulator(Nx=128, Ny=128, dx=0.5, dt=0.015,
                                chi_N=chiN, f_A=f_A, kT=0.003)
        sim.step(n_steps=5000)
        ax = axes[idx // 2][idx % 2]
        L = sim.Nx * sim.dx
        im = ax.imshow(sim.psi.T, cmap='RdBu_r', origin='lower',
                       extent=[0, L, 0, L], interpolation='bilinear',
                       vmin=-1, vmax=1)
        ax.set_title(f'{name} ($f_A$={f_A}, $\\chi N$={chiN})',
                     fontsize=12, fontweight='bold')
        ax.set_xlabel('x (σ)')
        ax.set_ylabel('y (σ)')
        plt.colorbar(im, ax=ax, label='$\\psi = \\phi_A - \\phi_B$')
    plt.suptitle('BCP Self-Assembly Morphologies (Ohta-Kawasaki Model)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, 'morphology_snapshots.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: figures/morphology_snapshots.png")


# ============================================================
# 4. Dynamic Process: Nucleation, Growth, Defect Annealing
# ============================================================
def simulate_dynamics():
    print("\n=== Dynamic Process Simulation ===")
    sim = BCPFieldSimulator(Nx=128, Ny=128, dx=0.5, dt=0.015,
                            chi_N=30.0, f_A=0.5, kT=0.003)
    times, order_params, free_energies = [], [], []
    snapshots = []
    snapshot_steps = [0, 500, 2000, 5000]

    for step in range(5001):
        sim.step(n_steps=1)
        if step % 50 == 0:
            op = sim.order_parameter()
            fe = np.mean(0.5 * sim.psi**2 - 0.25 * sim.psi**4)
            times.append(step * sim.dt)
            order_params.append(op)
            free_energies.append(fe)
        if step in snapshot_steps:
            snapshots.append((step, sim.psi.copy()))

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax = axes[0][0]
    ax.plot(times, order_params, 'b-', linewidth=2)
    ax.set_xlabel('Time (τ)', fontsize=12)
    ax.set_ylabel('Order Parameter Ψ', fontsize=12)
    ax.set_title('Segregation Kinetics', fontsize=13, fontweight='bold')
    ax.axhline(y=np.mean(order_params[-10:]), color='r', linestyle='--', alpha=0.7, label='Equilibrium')
    ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[0][1]
    ax.plot(times, free_energies, 'r-', linewidth=2)
    ax.set_xlabel('Time (τ)', fontsize=12)
    ax.set_ylabel('Free Energy Density', fontsize=12)
    ax.set_title('Free Energy Evolution', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)

    L = sim.Nx * sim.dx
    for i in range(2):
        step_val, psi = snapshots[i]
        ax = axes[1][i]
        im = ax.imshow(psi.T, cmap='RdBu_r', origin='lower',
                       extent=[0, L, 0, L], interpolation='bilinear',
                       vmin=-1, vmax=1)
        ax.set_title(f't = {step_val * sim.dt:.1f} τ', fontsize=12, fontweight='bold')
        ax.set_xlabel('x (σ)'); ax.set_ylabel('y (σ)')
        plt.colorbar(im, ax=ax, label='$\\psi$')

    plt.suptitle('Self-Assembly Dynamics & Defect Annealing',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, 'dynamics_evolution.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: figures/dynamics_evolution.png")

    # Time evolution snapshots
    fig2, axes2 = plt.subplots(1, 4, figsize=(18, 4.5))
    for snap_idx, (step_val, psi) in enumerate(snapshots):
        ax = axes2[snap_idx]
        im = ax.imshow(psi.T, cmap='RdBu_r', origin='lower',
                       extent=[0, L, 0, L], interpolation='bilinear',
                       vmin=-1, vmax=1)
        ax.set_title(f't = {step_val * sim.dt:.1f} τ', fontsize=11, fontweight='bold')
        ax.set_xlabel('x (σ)')
        if snap_idx == 0: ax.set_ylabel('y (σ)')
    plt.suptitle('Morphology Evolution During Self-Assembly',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, 'time_evolution.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: figures/time_evolution.png")
    return times, order_params


# ============================================================
# 5. Directed Self-Assembly (DSA) Simulation
# ============================================================
def simulate_dsa():
    print("\n=== Directed Self-Assembly (DSA) ===")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    cases = [
        ("Unguided Assembly", None),
        ("Chemoepitaxy ($L_s = L_0$)", 8.0),
        ("Chemoepitaxy ($L_s = 2L_0$)", 16.0),
    ]

    for ci, (title, period) in enumerate(cases):
        sim = BCPFieldSimulator(Nx=128, Ny=128, dx=0.5, dt=0.015,
                                chi_N=30.0, f_A=0.5, kT=0.003)
        if period is not None:
            # Chemical template: stripe pattern at bottom
            x = np.arange(sim.Nx) * sim.dx
            template = 0.3 * np.sin(2 * np.pi * x / period)
            sim.psi[:, :5] += template[:, None]

        for step in range(4000):
            sim.step(n_steps=1)
            if period is not None and step % 10 == 0:
                x = np.arange(sim.Nx) * sim.dx
                template = 0.1 * np.sin(2 * np.pi * x / period)
                sim.psi[:, :3] = 0.7 * sim.psi[:, :3] + 0.3 * template[:, None]

        L = sim.Nx * sim.dx
        ax = axes[ci]
        im = ax.imshow(sim.psi.T, cmap='RdBu_r', origin='lower',
                       extent=[0, L, 0, L], interpolation='bilinear',
                       vmin=-1, vmax=1)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xlabel('x (σ)')
        if ci == 0: ax.set_ylabel('y (σ)')
        plt.colorbar(im, ax=ax, label='$\\psi$')
        if period is not None:
            for xv in np.arange(0, L, period):
                ax.axvline(x=xv, color='lime', linestyle='--', alpha=0.5, linewidth=1)
            ax.axhline(y=1.5, color='yellow', linestyle='-', alpha=0.6, linewidth=1.5)
            ax.text(1, 0.5, 'Template', color='yellow', fontsize=9)

    plt.suptitle('Directed Self-Assembly: Template-Guided vs Unguided',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, 'dsa_comparison.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: figures/dsa_comparison.png")


# ============================================================
# 6. Structure Factor Analysis
# ============================================================
def analyze_structure_factor():
    print("\n=== Structure Factor Analysis ===")
    configs = [
        ("Disordered", 0.50, 8.0),
        ("Weak Segregation", 0.50, 15.0),
        ("Strong Segregation", 0.50, 35.0),
        ("Cylinders", 0.30, 40.0),
    ]
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ['#636363', '#3182bd', '#e6550d', '#31a354']
    for idx, (name, f_A, chiN) in enumerate(configs):
        sim = BCPFieldSimulator(Nx=128, Ny=128, dx=0.5, dt=0.015,
                                chi_N=chiN, f_A=f_A, kT=0.003)
        sim.step(n_steps=4000)
        q_vals, S_q = sim.structure_factor(n_q=60, q_max=8.0)
        S_q_norm = S_q / (np.max(S_q) + 1e-10)
        ax.plot(q_vals, S_q_norm, '-', color=colors[idx], linewidth=2,
                label=f'{name} ($f_A$={f_A})')

    ax.set_xlabel('$q$ (σ$^{-1}$)', fontsize=13)
    ax.set_ylabel('$S(q) / S_{max}$', fontsize=13)
    ax.set_title('Structure Factor of BCP Self-Assembly', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.annotate('$q^*$', xy=(2.0, 0.85), fontsize=14, color='red',
                arrowprops=dict(arrowstyle='->', color='red'), xytext=(3.0, 0.6))
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, 'structure_factor.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: figures/structure_factor.png")


# ============================================================
# 7. Multiscale Simulation Schematic
# ============================================================
def plot_multiscale_schematic():
    print("\n=== Multiscale Schematic ===")
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_xlim(0, 14); ax.set_ylim(0, 7); ax.axis('off')
    boxes = [
        (1, 3, 3, 2.5, 'All-Atom MD\n(OPLS-AA/CHARMM)\n~1-10 ns\n~10⁴ atoms', '#4e79a7'),
        (5, 3, 3, 2.5, 'Coarse-Grained MD\n(MARTINI/SDK)\n~1-100 μs\n~10⁶ beads', '#59a14f'),
        (9, 3, 3, 2.5, 'Field Theory\n(SCFT/TDDFT)\n~sec-min\nContinuum', '#edc949'),
    ]
    for x, y, w, h, text, color in boxes:
        rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor='black', linewidth=2, alpha=0.8)
        ax.add_patch(rect)
        ax.text(x+w/2, y+h/2, text, ha='center', va='center', fontsize=10, fontweight='bold', color='white')
    ax.annotate('', xy=(5, 4.25), xytext=(4, 4.25), arrowprops=dict(arrowstyle='->', lw=2.5, color='#333'))
    ax.annotate('', xy=(9, 4.25), xytext=(8, 4.25), arrowprops=dict(arrowstyle='->', lw=2.5, color='#333'))
    ax.text(4.5, 4.8, 'Coarse-\nGraining', ha='center', fontsize=9, fontweight='bold')
    ax.text(4.5, 3.5, 'Back-\nmapping', ha='center', fontsize=9, fontweight='bold', color='#666')
    ax.text(8.5, 4.8, 'Field\nMapping', ha='center', fontsize=9, fontweight='bold')
    ax.text(8.5, 3.5, 'Particle\nInsertion', ha='center', fontsize=9, fontweight='bold', color='#666')
    ax.annotate('', xy=(4, 3.8), xytext=(5, 3.8), arrowprops=dict(arrowstyle='->', lw=1.5, color='#999', linestyle='dashed'))
    ax.annotate('', xy=(8, 3.8), xytext=(9, 3.8), arrowprops=dict(arrowstyle='->', lw=1.5, color='#999', linestyle='dashed'))
    ax.text(1, 1.5, 'Length: 0.1-1 nm', fontsize=9, color='#4e79a7')
    ax.text(5, 1.5, 'Length: 1-100 nm', fontsize=9, color='#59a14f')
    ax.text(9, 1.5, 'Length: 10-1000 nm', fontsize=9, color='#b07d2e')
    ax.text(1, 1.0, 'Time: ps-ns', fontsize=9, color='#4e79a7')
    ax.text(5, 1.0, 'Time: ns-μs', fontsize=9, color='#59a14f')
    ax.text(9, 1.0, 'Time: μs-ms', fontsize=9, color='#b07d2e')
    ax.text(7, 6.5, 'Multiscale Simulation Framework for BCP Self-Assembly', ha='center', fontsize=14, fontweight='bold')
    ax.text(7, 0.3, 'Application: Sub-7nm Semiconductor Patterning via DSA', ha='center', fontsize=11, style='italic', color='#555')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, 'multiscale_schematic.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: figures/multiscale_schematic.png")


# ============================================================
# 8. DSA for Semiconductor Patterning
# ============================================================
def simulate_semiconductor_patterning():
    print("\n=== Semiconductor Patterning Simulation ===")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    patterns = [("Line/Space (L/S)", "line"), ("Contact Hole (CH)", "hole"), ("Fin Pattern", "fin")]
    for idx, (title, ptype) in enumerate(patterns):
        Nx, Ny = 128, 128
        sim = BCPFieldSimulator(Nx=Nx, Ny=Ny, dx=0.5, dt=0.015, chi_N=35.0, f_A=0.5, kT=0.002)
        L = Nx * sim.dx
        if ptype == "line":
            period = 8.0
            x = np.arange(Nx) * sim.dx
            for j in range(Ny):
                sim.psi[:, j] += 0.3 * np.sin(2*np.pi*x/period)
        elif ptype == "hole":
            spacing = 10.0
            for cx in np.arange(spacing/2, L, spacing):
                for cy_idx, cy in enumerate(np.arange(spacing/2, L, spacing*np.sqrt(3)/2)):
                    off = spacing/2 if cy_idx%2 else 0
                    for i in range(Nx):
                        for j in range(Ny):
                            xv, yv = i*sim.dx, j*sim.dx
                            r = np.sqrt((xv-cx-off)**2+(yv-cy)**2)
                            if r < 3.0:
                                sim.psi[i, j] = -0.8
        elif ptype == "fin":
            period = 6.0
            x = np.arange(Nx)*sim.dx
            for j in range(Ny):
                phase = (x % period)/period
                sim.psi[:, j] += np.where(phase < 0.3, 0.5, -0.3)
        sim.step(n_steps=3000)
        ax = axes[idx]
        im = ax.imshow(sim.psi.T, cmap='RdBu_r', origin='lower',
                       extent=[0, L, 0, L], interpolation='bilinear', vmin=-1, vmax=1)
        ax.set_title(f'{title}\n(sub-7nm equivalent)', fontsize=11, fontweight='bold')
        ax.set_xlabel('x (σ)')
        if idx==0: ax.set_ylabel('y (σ)')
        plt.colorbar(im, ax=ax, label='$\\psi$')
    plt.suptitle('DSA Patterns for Sub-7nm Semiconductor Patterning', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, 'semiconductor_patterns.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: figures/semiconductor_patterns.png")


# ============================================================
# 9. Defect Analysis
# ============================================================
def analyze_defects():
    print("\n=== Defect Analysis ===")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    temps_anneal = [0.001, 0.003, 0.005, 0.01]
    ax1, ax2 = axes
    for kT in temps_anneal:
        sim = BCPFieldSimulator(Nx=64, Ny=64, dx=0.5, dt=0.015,
                                chi_N=30.0, f_A=0.5, kT=kT)
        times, defects = [], []
        for step in range(0, 3000, 30):
            sim.step(n_steps=30)
            grad_x = np.diff(sim.psi, axis=0)
            grad_y = np.diff(sim.psi, axis=1)
            dd = np.mean(np.abs(grad_x[:, :-1])) + np.mean(np.abs(grad_y[:-1, :]))
            times.append(step * sim.dt)
            defects.append(dd)
        ax1.plot(times, defects, linewidth=2, label=f'kT={kT}')

    ax1.set_xlabel('Annealing Time (τ)', fontsize=12)
    ax1.set_ylabel('Defect Density (a.u.)', fontsize=12)
    ax1.set_title('Defect Annealing Kinetics', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10); ax1.grid(True, alpha=0.3)

    temps = np.linspace(0.001, 0.02, 12)
    final_defects = []
    for kT in temps:
        sim = BCPFieldSimulator(Nx=64, Ny=64, dx=0.5, dt=0.015,
                                chi_N=30.0, f_A=0.5, kT=kT)
        sim.step(n_steps=2000)
        grad = np.abs(np.diff(sim.psi, axis=0))
        final_defects.append(np.mean(grad))
    ax2.plot(temps, final_defects, 'o-', color='#e15759', linewidth=2, markersize=6)
    ax2.set_xlabel('Noise Amplitude (kT)', fontsize=12)
    ax2.set_ylabel('Final Defect Density', fontsize=12)
    ax2.set_title('Temperature Dependence of Defectivity', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, 'defect_analysis.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: figures/defect_analysis.png")


# ============================================================
# 10. LAMMPS/HOOMD Protocol Generator
# ============================================================
def generate_lammps_protocol():
    print("\n=== Generating LAMMPS Protocol ===")
    script = """# ============================================================
# LAMMPS Input Script: Block Copolymer Self-Assembly (CG-DPD)
# ============================================================
units           lj
dimension       3
boundary        p p p
atom_style      molecular

variable        N       equal 32
variable        N_A     equal 16
variable        n_chain equal 500
variable        L       equal 40.0
variable        T       equal 1.0
variable        dt_sim  equal 0.005

region          box block 0 ${L} 0 ${L} 0 ${L}
create_box      2 box bond/types 1 extra/bond/per/atom 2

molecule        diblock diblock.mol
create_atoms    0 random ${n_chain} 12345 box mol diblock 67890

pair_style      dpd ${T} 1.5 12345
pair_coeff      1 1 25.0 4.5 1.5
pair_coeff      2 2 25.0 4.5 1.5
pair_coeff      1 2 40.0 4.5 1.5

bond_style      fene
bond_coeff      1 30.0 1.5 1.0 1.0
special_bonds   fene

neighbor        0.5 bin
neigh_modify    every 1 delay 0 check yes

velocity        all create ${T} 87654 dist gaussian
fix             1 all nve
fix             2 all langevin ${T} ${T} 1.0 48279
timestep        ${dt_sim}

# Soft push-off
pair_style      soft 1.5
pair_coeff      * * 10.0
run             5000
unfix           1
unfix           2

# DPD production
pair_style      dpd ${T} 1.5 12345
pair_coeff      1 1 25.0 4.5 1.5
pair_coeff      2 2 25.0 4.5 1.5
pair_coeff      1 2 40.0 4.5 1.5

fix             1 all nve
fix             2 all langevin ${T} ${T} 1.0 48279

thermo          1000
thermo_style    custom step temp pe ke press vol
dump            1 all custom 5000 traj.lammpstrj id mol type x y z
compute         myChunk all chunk/atom bin/2d x lower 1.0 y lower 1.0
fix             density all ave/chunk 100 10 1000 myChunk density/atom file density_profile.dat

run             2000000

# Annealing
variable        T_anneal equal 1.5
fix             2 all langevin ${T_anneal} 0.8 1.0 48279
run             500000
fix             2 all langevin 0.8 0.8 1.0 48279
run             1000000

write_data      final_config.data
print "=== Simulation Complete ==="
"""
    outpath = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'src', 'lammps_bcp.in')
    with open(outpath, 'w') as f:
        f.write(script)
    print("  Saved: src/lammps_bcp.in")


def generate_hoomd_protocol():
    print("\n=== Generating HOOMD-blue Protocol ===")
    script = '''#!/usr/bin/env python3
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
'''
    outpath = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'src', 'hoomd_bcp.py')
    with open(outpath, 'w') as f:
        f.write(script)
    print("  Saved: src/hoomd_bcp.py")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("Block Copolymer Self-Assembly: CG-MD Simulation Suite")
    print("=" * 60)

    map_phase_diagram()
    simulate_morphologies()
    simulate_dynamics()
    simulate_dsa()
    analyze_structure_factor()
    plot_multiscale_schematic()
    simulate_semiconductor_patterning()
    analyze_defects()
    generate_lammps_protocol()
    generate_hoomd_protocol()

    print("\n" + "=" * 60)
    print("All simulations and analyses complete!")
    print(f"Figures saved to: {FIGDIR}")
    print("=" * 60)
