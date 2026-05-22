from __future__ import annotations

import json
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from pipeline_common import FIGURES_DIR, RESULTS_DIR, append_log, ensure_3d_molecule, save_json, set_global_seed

DEMO_SMILES = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "caffeine": "Cn1cnc2n(C)c(=O)n(C)c(=O)c12",
    "acetaminophen": "CC(=O)NC1=CC=C(O)C=C1",
    "warfarin": "CC(=O)CC(c1ccccc1)c1c(O)oc2ccccc2c1=O",
}


class SimpleMDSimulator:
    def __init__(self, temperature: float = 300.0, dt: float = 0.0015, dielectric: float = 4.0):
        self.temperature = temperature
        self.dt = dt
        self.dielectric = dielectric
        self.k_b = 0.0019872041
        self.epsilon = 0.12
        self.spring_k = 3.0

    def _molecule_arrays(self, mol: Chem.Mol):
        mol3d = ensure_3d_molecule(mol, seed=42)
        AllChem.ComputeGasteigerCharges(mol3d)
        conf = mol3d.GetConformer()
        positions = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z] for i in range(mol3d.GetNumAtoms())], dtype=float)
        charges = np.array([float(atom.GetProp("_GasteigerCharge")) if atom.HasProp("_GasteigerCharge") else 0.0 for atom in mol3d.GetAtoms()], dtype=float)
        masses = np.array([max(atom.GetMass(), 1.0) for atom in mol3d.GetAtoms()], dtype=float)
        sigmas = np.array([Chem.GetPeriodicTable().GetRvdw(atom.GetAtomicNum()) for atom in mol3d.GetAtoms()], dtype=float)
        bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondTypeAsDouble()) for bond in mol3d.GetBonds()]
        eq_lengths = {(min(i, j), max(i, j)): np.linalg.norm(positions[i] - positions[j]) for i, j, _ in bonds}
        return positions, charges, masses, sigmas, bonds, eq_lengths

    def _forces(self, positions, charges, sigmas, bonds, eq_lengths):
        n_atoms = len(positions)
        forces = np.zeros_like(positions)
        potential = 0.0
        bonded_pairs = {(min(i, j), max(i, j)): order for i, j, order in bonds}
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                rij = positions[j] - positions[i]
                dist = np.linalg.norm(rij) + 1e-8
                unit = rij / dist
                sigma = 0.5 * (sigmas[i] + sigmas[j])
                sr6 = (sigma / max(dist, 0.8)) ** 6
                lj_energy = 4 * self.epsilon * (sr6**2 - sr6)
                lj_force_mag = 24 * self.epsilon * (2 * sr6**2 - sr6) / max(dist, 0.8)
                coul_energy = 332.0636 * charges[i] * charges[j] / (self.dielectric * dist)
                coul_force_mag = coul_energy / dist
                pair_force = (lj_force_mag + coul_force_mag) * unit
                if (i, j) not in bonded_pairs:
                    forces[i] -= pair_force
                    forces[j] += pair_force
                    potential += lj_energy + coul_energy
        for i, j, order in bonds:
            key = (min(i, j), max(i, j))
            rij = positions[j] - positions[i]
            dist = np.linalg.norm(rij) + 1e-8
            unit = rij / dist
            delta = dist - eq_lengths[key]
            force = self.spring_k * order * delta * unit
            forces[i] += force
            forces[j] -= force
            potential += 0.5 * self.spring_k * order * delta**2
        return forces, potential

    def _kinetic_energy(self, velocities, masses):
        return 0.5 * np.sum(masses[:, None] * velocities**2)

    def run_simulation(self, mol: Chem.Mol, n_steps: int = 500) -> Dict:
        positions, charges, masses, sigmas, bonds, eq_lengths = self._molecule_arrays(mol)
        rng = np.random.default_rng(42 + len(positions))
        velocities = rng.normal(0.0, np.sqrt(self.k_b * self.temperature / masses)[:, None], size=positions.shape)
        initial_positions = positions.copy()
        forces, potential = self._forces(positions, charges, sigmas, bonds, eq_lengths)
        rmsd_traj, potential_traj, kinetic_traj = [], [], []
        for _ in range(n_steps):
            positions = positions + velocities * self.dt + 0.5 * (forces / masses[:, None]) * self.dt**2
            new_forces, potential = self._forces(positions, charges, sigmas, bonds, eq_lengths)
            velocities = velocities + 0.5 * (forces + new_forces) / masses[:, None] * self.dt
            kinetic = self._kinetic_energy(velocities, masses)
            inst_temp = (2.0 * kinetic) / (3.0 * len(masses) * self.k_b + 1e-8)
            scale = np.sqrt(max(0.2, 1 + self.dt / 0.05 * (self.temperature / max(inst_temp, 1e-6) - 1)))
            velocities *= scale
            forces = new_forces
            centered_initial = initial_positions - initial_positions.mean(axis=0)
            centered_current = positions - positions.mean(axis=0)
            rmsd = float(np.sqrt(np.mean(np.sum((centered_initial - centered_current) ** 2, axis=1))))
            rmsd_traj.append(rmsd)
            potential_traj.append(float(potential))
            kinetic_traj.append(float(kinetic))
        return {
            "initial_coords": initial_positions.tolist(),
            "final_coords": positions.tolist(),
            "rmsd_trajectory": rmsd_traj,
            "potential_energy": potential_traj,
            "kinetic_energy": kinetic_traj,
            "temperature": self.temperature,
            "n_steps": n_steps,
        }

    def compute_binding_pose_rmsd(self, initial_coords, final_coords) -> float:
        initial = np.asarray(initial_coords, dtype=float)
        final = np.asarray(final_coords, dtype=float)
        initial -= initial.mean(axis=0)
        final -= final.mean(axis=0)
        return float(np.sqrt(np.mean(np.sum((initial - final) ** 2, axis=1))))

    def analyze_energy_convergence(self, energy_trajectory) -> Dict:
        arr = np.asarray(energy_trajectory, dtype=float)
        last = arr[-100:] if len(arr) >= 100 else arr
        x = np.arange(len(last))
        slope = float(np.polyfit(x, last, 1)[0]) if len(last) > 1 else 0.0
        return {
            "mean_last_window": float(last.mean()),
            "std_last_window": float(last.std()),
            "slope_last_window": slope,
            "delta_total": float(arr[-1] - arr[0]),
            "converged": bool(abs(slope) < 0.05 and last.std() < max(5.0, abs(last.mean()) * 0.1)),
        }


def save_plots(results: Dict[str, Dict]) -> None:
    cmap = plt.get_cmap("viridis")
    fig1, axes = plt.subplots(3, 1, figsize=(14, 12), dpi=150, sharex=True)
    for idx, (name, data) in enumerate(results.items()):
        color = cmap(0.15 + 0.15 * idx)
        steps = np.arange(len(data["trajectory"]["potential_energy"]))
        axes[0].plot(steps, data["trajectory"]["potential_energy"], label=name, color=color)
        axes[1].plot(steps, data["trajectory"]["kinetic_energy"], label=name, color=color)
        total = np.array(data["trajectory"]["potential_energy"]) + np.array(data["trajectory"]["kinetic_energy"])
        axes[2].plot(steps, total, label=name, color=color)
    axes[0].set_ylabel("Potential energy")
    axes[1].set_ylabel("Kinetic energy")
    axes[2].set_ylabel("Total energy")
    axes[2].set_xlabel("Step")
    axes[0].set_title("MD energy convergence")
    axes[0].legend(ncol=2)
    fig1.tight_layout()
    fig1.savefig(FIGURES_DIR / "md_energy_convergence.png", bbox_inches="tight")
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(12, 8), dpi=150)
    for idx, (name, data) in enumerate(results.items()):
        ax2.plot(data["trajectory"]["rmsd_trajectory"], label=name, color=cmap(0.15 + 0.15 * idx))
    ax2.set_xlabel("Step")
    ax2.set_ylabel("RMSD (Å)")
    ax2.set_title("Binding pose RMSD trajectories")
    ax2.legend(ncol=2)
    fig2.tight_layout()
    fig2.savefig(FIGURES_DIR / "md_rmsd_analysis.png", bbox_inches="tight")
    plt.close(fig2)


def main() -> Dict:
    set_global_seed(42)
    append_log("md_simulation", "run_started", "02_md_simulation.py", {"n_steps": 500})
    simulator = SimpleMDSimulator()
    results = {}
    for name, smiles in DEMO_SMILES.items():
        mol = Chem.MolFromSmiles(smiles)
        trajectory = simulator.run_simulation(mol, n_steps=500)
        results[name] = {
            "trajectory": trajectory,
            "final_rmsd": simulator.compute_binding_pose_rmsd(trajectory["initial_coords"], trajectory["final_coords"]),
            "energy_convergence": simulator.analyze_energy_convergence(trajectory["potential_energy"]),
        }
    save_plots(results)
    payload = {
        "molecules": results,
        "summary": {
            "mean_final_rmsd": float(np.mean([v["final_rmsd"] for v in results.values()])),
            "converged_fraction": float(np.mean([v["energy_convergence"]["converged"] for v in results.values()])),
        },
    }
    json_path = RESULTS_DIR / "md_simulation_results.json"
    save_json(json_path, payload)
    append_log("md_simulation", "run_completed", "02_md_simulation.py", {"molecule_count": len(results)}, payload["summary"], [str(json_path), str(FIGURES_DIR / "md_energy_convergence.png"), str(FIGURES_DIR / "md_rmsd_analysis.png")])
    print(json.dumps({"status": "ok", "result_file": str(json_path)}))
    return payload


if __name__ == "__main__":
    main()
