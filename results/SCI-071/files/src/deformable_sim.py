"""
Deformable Object Manipulation Planning System
- State representations: Mesh, Particle, Latent Space
- Physics simulation: FEM / MPM
- Manipulation sequence planning
- Sim-to-Real domain randomization
- Visual feedback reactive control
- Cloth folding case study
"""

import numpy as np
from scipy.spatial import Delaunay
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
from sklearn.decomposition import PCA
import json, os

np.random.seed(42)

# ============================================================
# 1. State Representations
# ============================================================

class MeshRepresentation:
    """Triangular mesh representation for deformable objects."""
    def __init__(self, n_x=10, n_y=10, size=1.0):
        self.n_x, self.n_y = n_x, n_y
        x = np.linspace(0, size, n_x)
        y = np.linspace(0, size, n_y)
        xx, yy = np.meshgrid(x, y)
        self.vertices = np.column_stack([xx.ravel(), yy.ravel(), np.zeros(n_x * n_y)])
        self.rest_vertices = self.vertices.copy()
        pts_2d = self.vertices[:, :2]
        tri = Delaunay(pts_2d)
        self.triangles = tri.simplices
        self.n_vertices = len(self.vertices)

    def get_state(self):
        return self.vertices.copy()

    def set_state(self, vertices):
        self.vertices = vertices.copy()

    def compute_strain_energy(self):
        energy = 0.0
        for tri in self.triangles:
            rest = self.rest_vertices[tri]
            curr = self.vertices[tri]
            F = self._deformation_gradient(rest[:, :2], curr[:, :2])
            E = 0.5 * (F.T @ F - np.eye(2))
            energy += np.sum(E ** 2)
        return energy

    def _deformation_gradient(self, rest, curr):
        D_rest = np.column_stack([rest[1] - rest[0], rest[2] - rest[0]])
        D_curr = np.column_stack([curr[1] - curr[0], curr[2] - curr[0]])
        try:
            F = D_curr @ np.linalg.inv(D_rest)
        except np.linalg.LinAlgError:
            F = np.eye(2)
        return F


class ParticleRepresentation:
    """Particle-based representation (MPM-style)."""
    def __init__(self, n_particles=200, size=1.0):
        self.n_particles = n_particles
        self.positions = np.random.rand(n_particles, 3) * size
        self.positions[:, 2] = 0
        self.velocities = np.zeros((n_particles, 3))
        self.masses = np.ones(n_particles) * 0.01
        self.rest_positions = self.positions.copy()

    def get_state(self):
        return np.concatenate([self.positions, self.velocities], axis=1)

    def compute_density(self, grid_size=20):
        grid = np.zeros((grid_size, grid_size))
        for p in self.positions:
            ix = int(np.clip(p[0] * grid_size, 0, grid_size - 1))
            iy = int(np.clip(p[1] * grid_size, 0, grid_size - 1))
            grid[ix, iy] += 1
        return grid


class LatentSpaceRepresentation:
    """Learned latent space via PCA (proxy for VAE)."""
    def __init__(self, latent_dim=8):
        self.latent_dim = latent_dim
        self.pca = PCA(n_components=latent_dim)
        self.fitted = False

    def fit(self, states):
        flat = states.reshape(len(states), -1)
        self.pca.fit(flat)
        self.fitted = True
        return self.pca.explained_variance_ratio_

    def encode(self, state):
        flat = state.reshape(1, -1)
        return self.pca.transform(flat)[0]

    def decode(self, latent):
        return self.pca.inverse_transform(latent.reshape(1, -1))[0]

    def reconstruction_error(self, state):
        z = self.encode(state)
        recon = self.decode(z)
        flat = state.ravel()
        return np.mean((flat - recon) ** 2)


# ============================================================
# 2. Physics Simulators
# ============================================================

class FEMSimulator:
    """Simplified 2D FEM for cloth-like deformation (neo-Hookean)."""
    def __init__(self, mesh: MeshRepresentation, youngs=1000.0, poisson=0.3, dt=0.01):
        self.mesh = mesh
        self.E = youngs
        self.nu = poisson
        self.dt = dt
        self.mu = self.E / (2 * (1 + self.nu))
        self.lam = self.E * self.nu / ((1 + self.nu) * (1 - 2 * self.nu))
        self.velocities = np.zeros_like(mesh.vertices)
        self.gravity = np.array([0, 0, -9.81])
        self.damping = 0.98

    def step(self, external_forces=None):
        n = self.mesh.n_vertices
        forces = np.zeros((n, 3))
        forces += self.mesh.masses_array[:, None] * self.gravity if hasattr(self.mesh, 'masses_array') else 0.01 * self.gravity
        for tri_idx in self.mesh.triangles:
            f_elastic = self._elastic_force(tri_idx)
            for i, vi in enumerate(tri_idx):
                forces[vi] += f_elastic[i]
        if external_forces is not None:
            forces += external_forces
        self.velocities += forces * self.dt
        self.velocities *= self.damping
        self.mesh.vertices += self.velocities * self.dt

    def _elastic_force(self, tri_idx):
        rest = self.mesh.rest_vertices[tri_idx, :2]
        curr = self.mesh.vertices[tri_idx, :2]
        D_rest = np.column_stack([rest[1] - rest[0], rest[2] - rest[0]])
        D_curr = np.column_stack([curr[1] - curr[0], curr[2] - curr[0]])
        det = np.linalg.det(D_rest)
        if abs(det) < 1e-10:
            return np.zeros((3, 3))
        F = D_curr @ np.linalg.inv(D_rest)
        P = self.mu * (F - np.linalg.inv(F.T)) + self.lam * np.log(max(np.linalg.det(F), 1e-6)) * np.linalg.inv(F.T)
        H = -abs(det) * 0.5 * P @ np.linalg.inv(D_rest.T)
        forces = np.zeros((3, 3))
        forces[1, :2] = H[:, 0]
        forces[2, :2] = H[:, 1]
        forces[0, :2] = -H[:, 0] - H[:, 1]
        return forces * 0.001


class MPMSimulator:
    """Simplified Material Point Method simulator."""
    def __init__(self, particles: ParticleRepresentation, grid_size=32, dt=0.005):
        self.particles = particles
        self.grid_size = grid_size
        self.dt = dt
        self.gravity = np.array([0, 0, -9.81])
        self.grid_v = np.zeros((grid_size, grid_size, 3))
        self.grid_m = np.zeros((grid_size, grid_size))

    def step(self, external_forces=None):
        gs = self.grid_size
        self.grid_v[:] = 0
        self.grid_m[:] = 0
        # P2G
        for i in range(self.particles.n_particles):
            pos = self.particles.positions[i]
            gx = int(np.clip(pos[0] * gs, 1, gs - 2))
            gy = int(np.clip(pos[1] * gs, 1, gs - 2))
            w = self.particles.masses[i]
            self.grid_m[gx, gy] += w
            self.grid_v[gx, gy] += w * self.particles.velocities[i]
        # Grid update
        for i in range(gs):
            for j in range(gs):
                if self.grid_m[i, j] > 1e-8:
                    self.grid_v[i, j] /= self.grid_m[i, j]
                    self.grid_v[i, j] += self.gravity * self.dt
                    self.grid_v[i, j] *= 0.99
                    if j < 2:
                        self.grid_v[i, j, 2] = max(self.grid_v[i, j, 2], 0)
        # G2P
        for i in range(self.particles.n_particles):
            pos = self.particles.positions[i]
            gx = int(np.clip(pos[0] * gs, 1, gs - 2))
            gy = int(np.clip(pos[1] * gs, 1, gs - 2))
            self.particles.velocities[i] = self.grid_v[gx, gy]
            self.particles.positions[i] += self.particles.velocities[i] * self.dt
        if external_forces is not None:
            self.particles.velocities += external_forces * self.dt


# ============================================================
# 3. Manipulation Sequence Planning
# ============================================================

class ManipulationPlanner:
    """Plans manipulation sequences to reach target states."""
    def __init__(self, state_dim, action_dim=6, horizon=20):
        self.state_dim = state_dim
        self.action_dim = action_dim
        self.horizon = horizon

    def plan_cem(self, current_state, target_state, dynamics_fn, n_iter=50, n_samples=200, elite_frac=0.1):
        """Cross-Entropy Method planning."""
        n_elite = max(int(n_samples * elite_frac), 5)
        mean = np.zeros(self.horizon * self.action_dim)
        std = np.ones(self.horizon * self.action_dim) * 0.5
        best_cost = float('inf')
        best_actions = None
        cost_history = []

        for it in range(n_iter):
            samples = np.random.randn(n_samples, len(mean)) * std + mean
            costs = np.zeros(n_samples)
            for i in range(n_samples):
                actions = samples[i].reshape(self.horizon, self.action_dim)
                state = current_state.copy()
                for a in actions:
                    state = dynamics_fn(state, a)
                costs[i] = np.sum((state - target_state) ** 2)
            elite_idx = np.argsort(costs)[:n_elite]
            elite = samples[elite_idx]
            mean = elite.mean(axis=0)
            std = elite.std(axis=0) + 1e-5
            if costs[elite_idx[0]] < best_cost:
                best_cost = costs[elite_idx[0]]
                best_actions = samples[elite_idx[0]].reshape(self.horizon, self.action_dim)
            cost_history.append(best_cost)

        return best_actions, cost_history

    def plan_rrt(self, current_state, target_state, dynamics_fn, max_iter=1000, step_size=0.1):
        """RRT-based planning in state space."""
        tree = [current_state.copy()]
        parents = [-1]
        actions_list = [None]

        for _ in range(max_iter):
            if np.random.rand() < 0.1:
                sample = target_state
            else:
                sample = current_state + np.random.randn(*current_state.shape) * 2
            dists = [np.linalg.norm(n - sample) for n in tree]
            nearest_idx = np.argmin(dists)
            nearest = tree[nearest_idx]
            direction = sample - nearest
            norm = np.linalg.norm(direction)
            if norm > 1e-6:
                direction = direction / norm * min(step_size, norm)
            action = direction[:self.action_dim] if len(direction) >= self.action_dim else np.zeros(self.action_dim)
            new_state = dynamics_fn(nearest, action)
            tree.append(new_state)
            parents.append(nearest_idx)
            actions_list.append(action)
            if np.linalg.norm(new_state - target_state) < 0.5:
                path = []
                idx = len(tree) - 1
                while idx > 0:
                    path.append(actions_list[idx])
                    idx = parents[idx]
                return list(reversed(path)), True
        return [], False


# ============================================================
# 4. Domain Randomization for Sim-to-Real
# ============================================================

class DomainRandomizer:
    """Applies domain randomization to simulation parameters."""
    def __init__(self, base_params):
        self.base_params = base_params
        self.ranges = {
            'youngs_modulus': (500, 2000),
            'poisson_ratio': (0.2, 0.45),
            'friction': (0.1, 0.8),
            'mass_scale': (0.5, 2.0),
            'gravity_noise': (-0.5, 0.5),
            'action_noise_std': (0.0, 0.05),
            'observation_noise_std': (0.0, 0.02),
        }
        self.history = []

    def sample(self):
        params = {}
        for key, (lo, hi) in self.ranges.items():
            params[key] = np.random.uniform(lo, hi)
        self.history.append(params)
        return params

    def get_randomized_dynamics(self, base_dynamics_fn, params):
        def randomized_fn(state, action):
            noisy_action = action + np.random.randn(*action.shape) * params.get('action_noise_std', 0)
            new_state = base_dynamics_fn(state, noisy_action)
            new_state += np.random.randn(*new_state.shape) * params.get('observation_noise_std', 0)
            return new_state
        return randomized_fn

    def evaluate_transfer_gap(self, policy_fn, real_dynamics, sim_dynamics, initial_state, target_state, n_steps=20):
        """Evaluate sim-to-real gap for a given policy."""
        state_sim = initial_state.copy()
        state_real = initial_state.copy()
        for _ in range(n_steps):
            action = policy_fn(state_sim, target_state)
            state_sim = sim_dynamics(state_sim, action)
            state_real = real_dynamics(state_real, action)
        sim_err = np.linalg.norm(state_sim - target_state)
        real_err = np.linalg.norm(state_real - target_state)
        return sim_err, real_err, abs(sim_err - real_err)


# ============================================================
# 5. Visual Feedback Reactive Controller
# ============================================================

class VisualFeedbackController:
    """Reactive controller using visual (state) feedback."""
    def __init__(self, kp=1.0, kd=0.1, max_force=5.0):
        self.kp = kp
        self.kd = kd
        self.max_force = max_force
        self.prev_error = None

    def compute_action(self, current_state, target_state):
        error = target_state - current_state
        if self.prev_error is None:
            d_error = np.zeros_like(error)
        else:
            d_error = error - self.prev_error
        self.prev_error = error.copy()
        action = self.kp * error + self.kd * d_error
        norm = np.linalg.norm(action)
        if norm > self.max_force:
            action = action / norm * self.max_force
        return action

    def run_episode(self, initial_state, target_state, dynamics_fn, max_steps=100, threshold=0.1):
        state = initial_state.copy()
        trajectory = [state.copy()]
        errors = []
        for step in range(max_steps):
            action = self.compute_action(state, target_state)
            state = dynamics_fn(state, action)
            trajectory.append(state.copy())
            err = np.linalg.norm(state - target_state)
            errors.append(err)
            if err < threshold:
                break
        return np.array(trajectory), errors


# ============================================================
# 6. Cloth Folding Case Study
# ============================================================

class ClothFoldingEnv:
    """Simulated cloth folding environment."""
    def __init__(self, grid_n=8):
        self.grid_n = grid_n
        self.mesh = MeshRepresentation(n_x=grid_n, n_y=grid_n, size=1.0)
        self.target_vertices = self._compute_fold_target()
        self.grasp_points = self._select_grasp_points()

    def _compute_fold_target(self):
        target = self.mesh.rest_vertices.copy()
        n = self.grid_n
        for i in range(n * n):
            x, y, z = target[i]
            if x > 0.5:
                target[i] = [1.0 - x, y, 0.02]
        return target

    def _select_grasp_points(self):
        n = self.grid_n
        top_right = n - 1
        bottom_right = n * n - 1
        return [top_right, bottom_right]

    def get_state(self):
        return self.mesh.get_state()

    def get_target(self):
        return self.target_vertices

    def compute_reward(self):
        diff = self.mesh.vertices - self.target_vertices
        return -np.mean(np.sum(diff ** 2, axis=1))

    def compute_coverage(self):
        diff = np.linalg.norm(self.mesh.vertices - self.target_vertices, axis=1)
        return np.mean(diff < 0.05)

    def step(self, action):
        for gp in self.grasp_points:
            self.mesh.vertices[gp] += action[:3] * 0.1
        for i in range(self.mesh.n_vertices):
            if i not in self.grasp_points:
                for gp in self.grasp_points:
                    dist = np.linalg.norm(self.mesh.rest_vertices[i] - self.mesh.rest_vertices[gp])
                    influence = np.exp(-dist * 3)
                    self.mesh.vertices[i] += action[:3] * 0.1 * influence
        return self.get_state(), self.compute_reward()

    def reset(self):
        self.mesh.vertices = self.mesh.rest_vertices.copy()
        return self.get_state()
