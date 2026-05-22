"""
Learning from Demonstration (LfD) Module
=========================================
GMM/GMR-based trajectory learning from expert surgical demonstrations.
Supports DMP (Dynamic Movement Primitives) as an alternative policy representation.
"""

import numpy as np
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass, field


@dataclass
class Demonstration:
    """Single demonstration trajectory."""
    positions: np.ndarray       # (T, 3) end-effector positions
    orientations: np.ndarray    # (T, 4) quaternions
    gripper_angles: np.ndarray  # (T,) gripper jaw angle
    forces: np.ndarray          # (T, 3) measured forces
    timestamps: np.ndarray      # (T,) time stamps [s]
    metadata: Dict = field(default_factory=dict)

    @property
    def duration(self) -> float:
        return self.timestamps[-1] - self.timestamps[0]

    @property
    def num_steps(self) -> int:
        return len(self.timestamps)


class GaussianMixtureRegression:
    """
    GMM/GMR for encoding and reproducing demonstrated trajectories.
    Input (time) -> Output (pose, force) regression via Gaussian conditioning.
    """

    def __init__(self, n_components: int = 5, reg_covar: float = 1e-4):
        self.n_components = n_components
        self.reg_covar = reg_covar
        self.weights_: Optional[np.ndarray] = None
        self.means_: Optional[np.ndarray] = None
        self.covariances_: Optional[np.ndarray] = None
        self._fitted = False

    def fit(self, demonstrations: List[Demonstration],
            use_dtw_alignment: bool = True) -> 'GaussianMixtureRegression':
        """
        Fit GMM to aligned demonstration data.

        Parameters
        ----------
        demonstrations : list of Demonstration
            Expert demonstrations to learn from.
        use_dtw_alignment : bool
            If True, align demonstrations using DTW before fitting.
        """
        aligned = self._align_demonstrations(demonstrations, use_dtw_alignment)

        # Build data matrix: [time, x, y, z, qw, qx, qy, qz, gripper]
        data_blocks = []
        for demo in aligned:
            t_norm = (demo.timestamps - demo.timestamps[0]) / demo.duration
            block = np.column_stack([
                t_norm, demo.positions, demo.orientations, demo.gripper_angles
            ])
            data_blocks.append(block)
        data = np.vstack(data_blocks)

        # EM algorithm for GMM fitting
        self._fit_em(data)
        self._fitted = True
        return self

    def predict(self, t_query: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Predict pose trajectory via GMR conditioning on time.

        Returns
        -------
        means : (N, 8) predicted [x,y,z, qw,qx,qy,qz, gripper]
        covariances : (N, 8, 8) predicted covariance matrices
        """
        if not self._fitted:
            raise RuntimeError("Model not fitted. Call fit() first.")

        N = len(t_query)
        dim_out = self.means_.shape[1] - 1
        pred_means = np.zeros((N, dim_out))
        pred_covs = np.zeros((N, dim_out, dim_out))

        for i, t in enumerate(t_query):
            mu, sigma = self._gmr_condition(t)
            pred_means[i] = mu
            pred_covs[i] = sigma

        return pred_means, pred_covs

    def _align_demonstrations(self, demos: List[Demonstration],
                               use_dtw: bool) -> List[Demonstration]:
        """Temporally align demonstrations using DTW or linear interpolation."""
        if not use_dtw or len(demos) <= 1:
            return demos

        # Resample all demos to median length
        lengths = [d.num_steps for d in demos]
        target_len = int(np.median(lengths))

        aligned = []
        for demo in demos:
            t_new = np.linspace(demo.timestamps[0], demo.timestamps[-1], target_len)
            positions = np.array([
                np.interp(t_new, demo.timestamps, demo.positions[:, j])
                for j in range(3)
            ]).T
            orientations = np.array([
                np.interp(t_new, demo.timestamps, demo.orientations[:, j])
                for j in range(4)
            ]).T
            # Re-normalize quaternions
            orientations /= np.linalg.norm(orientations, axis=1, keepdims=True)
            gripper = np.interp(t_new, demo.timestamps, demo.gripper_angles)
            forces = np.array([
                np.interp(t_new, demo.timestamps, demo.forces[:, j])
                for j in range(3)
            ]).T

            aligned.append(Demonstration(
                positions=positions, orientations=orientations,
                gripper_angles=gripper, forces=forces,
                timestamps=t_new, metadata=demo.metadata
            ))
        return aligned

    def _fit_em(self, data: np.ndarray, max_iter: int = 100, tol: float = 1e-4):
        """Expectation-Maximization for GMM."""
        N, D = data.shape
        K = self.n_components

        # Initialize with k-means-like heuristic
        indices = np.random.choice(N, K, replace=False)
        self.means_ = data[indices].copy()
        self.covariances_ = np.array([np.eye(D) * 0.1 for _ in range(K)])
        self.weights_ = np.ones(K) / K

        log_likelihood_prev = -np.inf

        for iteration in range(max_iter):
            # E-step
            responsibilities = self._e_step(data)

            # M-step
            self._m_step(data, responsibilities)

            # Check convergence
            log_likelihood = self._log_likelihood(data)
            if abs(log_likelihood - log_likelihood_prev) < tol:
                break
            log_likelihood_prev = log_likelihood

    def _e_step(self, data: np.ndarray) -> np.ndarray:
        """Compute responsibilities."""
        N = data.shape[0]
        K = self.n_components
        resp = np.zeros((N, K))

        for k in range(K):
            resp[:, k] = self.weights_[k] * self._gaussian_pdf(
                data, self.means_[k], self.covariances_[k]
            )

        resp_sum = resp.sum(axis=1, keepdims=True)
        resp_sum = np.maximum(resp_sum, 1e-300)
        return resp / resp_sum

    def _m_step(self, data: np.ndarray, resp: np.ndarray):
        """Update GMM parameters."""
        N = data.shape[0]
        K = self.n_components

        for k in range(K):
            Nk = resp[:, k].sum()
            if Nk < 1e-10:
                continue
            self.weights_[k] = Nk / N
            self.means_[k] = (resp[:, k:k+1].T @ data) / Nk
            diff = data - self.means_[k]
            self.covariances_[k] = (
                (diff * resp[:, k:k+1]).T @ diff
            ) / Nk + np.eye(data.shape[1]) * self.reg_covar

    def _gaussian_pdf(self, x: np.ndarray, mean: np.ndarray,
                      cov: np.ndarray) -> np.ndarray:
        """Multivariate Gaussian probability density."""
        D = len(mean)
        diff = x - mean
        try:
            cov_inv = np.linalg.inv(cov)
            det = np.linalg.det(cov)
        except np.linalg.LinAlgError:
            return np.zeros(len(x))
        det = max(det, 1e-300)
        norm = 1.0 / (np.sqrt((2 * np.pi) ** D * det))
        exponent = -0.5 * np.sum(diff @ cov_inv * diff, axis=1)
        return norm * np.exp(exponent)

    def _gmr_condition(self, t: float) -> Tuple[np.ndarray, np.ndarray]:
        """Condition GMM on time input to predict output."""
        K = self.n_components
        dim_out = self.means_.shape[1] - 1

        # Compute component weights for input t
        beta = np.zeros(K)
        for k in range(K):
            mu_in = self.means_[k, 0]
            sigma_in = self.covariances_[k, 0, 0]
            beta[k] = self.weights_[k] * np.exp(
                -0.5 * (t - mu_in) ** 2 / sigma_in
            ) / np.sqrt(2 * np.pi * sigma_in)
        beta_sum = beta.sum()
        if beta_sum < 1e-300:
            beta = np.ones(K) / K
        else:
            beta /= beta_sum

        # Conditional mean and covariance
        mu_out = np.zeros(dim_out)
        sigma_out = np.zeros((dim_out, dim_out))

        for k in range(K):
            mu_in_k = self.means_[k, 0]
            mu_out_k = self.means_[k, 1:]
            sigma_in_k = self.covariances_[k, 0, 0]
            sigma_out_in_k = self.covariances_[k, 1:, 0]
            sigma_out_k = self.covariances_[k, 1:, 1:]

            # Conditional mean: mu_out + sigma_out_in * sigma_in^-1 * (t - mu_in)
            cond_mu = mu_out_k + sigma_out_in_k * (t - mu_in_k) / sigma_in_k
            cond_sigma = sigma_out_k - np.outer(
                sigma_out_in_k, sigma_out_in_k
            ) / sigma_in_k

            mu_out += beta[k] * cond_mu
            sigma_out += beta[k] * (cond_sigma + np.outer(cond_mu, cond_mu))

        sigma_out -= np.outer(mu_out, mu_out)
        return mu_out, sigma_out

    def _log_likelihood(self, data: np.ndarray) -> float:
        """Compute log-likelihood of data under the model."""
        N = data.shape[0]
        ll = np.zeros(N)
        for k in range(self.n_components):
            ll += self.weights_[k] * self._gaussian_pdf(
                data, self.means_[k], self.covariances_[k]
            )
        return np.sum(np.log(np.maximum(ll, 1e-300)))


class DynamicMovementPrimitive:
    """
    DMP for encoding point-to-point suturing motions.
    Learns nonlinear forcing functions from demonstrations.
    """

    def __init__(self, n_basis: int = 25, alpha_y: float = 25.0,
                 beta_y: float = 6.25, alpha_x: float = 1.0, dt: float = 0.01):
        self.n_basis = n_basis
        self.alpha_y = alpha_y
        self.beta_y = beta_y
        self.alpha_x = alpha_x
        self.dt = dt
        self.weights: Optional[np.ndarray] = None
        self.centers: np.ndarray = np.zeros(0)
        self.widths: np.ndarray = np.zeros(0)

    def learn(self, demo: Demonstration, dim: int = 3):
        """Learn DMP weights from a single demonstration."""
        T = demo.num_steps
        tau = demo.duration

        # Setup basis function centers and widths (log-spaced in phase)
        self.centers = np.exp(-self.alpha_x * np.linspace(0, 1, self.n_basis))
        self.widths = 1.0 / (np.diff(self.centers) ** 2)
        self.widths = np.append(self.widths, self.widths[-1])

        self.weights = np.zeros((dim, self.n_basis))

        for d in range(dim):
            y = demo.positions[:, d]
            dy = np.gradient(y, demo.timestamps) * tau
            ddy = np.gradient(dy, demo.timestamps) * tau

            g = y[-1]
            y0 = y[0]

            # Compute target forcing function
            f_target = ddy - self.alpha_y * (self.beta_y * (g - y) - dy)

            # Compute phase variable
            x = np.exp(-self.alpha_x * np.linspace(0, 1, T))

            # Basis function activations
            Psi = np.exp(-self.widths[None, :] * (x[:, None] - self.centers[None, :]) ** 2)
            Psi_sum = Psi.sum(axis=1, keepdims=True)
            Psi_norm = Psi / np.maximum(Psi_sum, 1e-10)

            # Weighted least-squares
            sx = x * (g - y0)
            for b in range(self.n_basis):
                num = (Psi_norm[:, b] * sx * f_target).sum()
                den = (Psi_norm[:, b] * sx ** 2).sum()
                self.weights[d, b] = num / max(den, 1e-10)

    def generate(self, y0: np.ndarray, goal: np.ndarray,
                 tau: float = 1.0, T: int = 200) -> np.ndarray:
        """
        Generate trajectory from start to goal using learned forcing function.

        Returns
        -------
        trajectory : (T, 3) positions
        """
        if self.weights is None:
            raise RuntimeError("DMP not learned. Call learn() first.")

        dim = len(y0)
        y = y0.copy().astype(float)
        dy = np.zeros(dim)
        x = 1.0

        trajectory = np.zeros((T, dim))

        for t in range(T):
            # Basis function activations
            psi = np.exp(-self.widths * (x - self.centers) ** 2)
            psi_norm = psi / max(psi.sum(), 1e-10)

            for d in range(dim):
                f = (self.weights[d] * psi_norm).sum() * x * (goal[d] - y0[d])

                ddy = self.alpha_y * (self.beta_y * (goal[d] - y[d]) - dy[d]) + f
                dy[d] += ddy * self.dt / tau
                y[d] += dy[d] * self.dt / tau

            x -= self.alpha_x * x * self.dt / tau
            trajectory[t] = y.copy()

        return trajectory


class SuturingLfDPipeline:
    """
    Complete LfD pipeline for suturing task decomposition and learning.

    Suturing is decomposed into phases:
    1. Approach: Move needle to entry point
    2. Insert: Drive needle through tissue
    3. Pull-through: Pull suture thread
    4. Knot-tying: Form and tighten knot
    """

    PHASES = ['approach', 'insert', 'pull_through', 'knot_tying']

    def __init__(self, method: str = 'gmr', n_components: int = 5):
        self.method = method
        self.models: Dict[str, object] = {}
        self.demonstrations: Dict[str, List[Demonstration]] = {
            phase: [] for phase in self.PHASES
        }

        for phase in self.PHASES:
            if method == 'gmr':
                self.models[phase] = GaussianMixtureRegression(n_components)
            elif method == 'dmp':
                self.models[phase] = DynamicMovementPrimitive()
            else:
                raise ValueError(f"Unknown method: {method}")

    def add_demonstration(self, phase: str, demo: Demonstration):
        """Add a demonstration for a specific suturing phase."""
        if phase not in self.PHASES:
            raise ValueError(f"Unknown phase: {phase}. Must be one of {self.PHASES}")
        self.demonstrations[phase].append(demo)

    def learn_all_phases(self):
        """Learn models for all phases with available demonstrations."""
        for phase in self.PHASES:
            demos = self.demonstrations[phase]
            if not demos:
                continue

            if self.method == 'gmr':
                self.models[phase].fit(demos)
            elif self.method == 'dmp':
                self.models[phase].learn(demos[0])

    def generate_trajectory(self, phase: str,
                           t_query: Optional[np.ndarray] = None,
                           goal: Optional[np.ndarray] = None) -> np.ndarray:
        """Generate trajectory for a specific suturing phase."""
        model = self.models[phase]

        if self.method == 'gmr':
            if t_query is None:
                t_query = np.linspace(0, 1, 200)
            means, covs = model.predict(t_query)
            return means[:, :3]  # Return positions only
        elif self.method == 'dmp':
            demos = self.demonstrations[phase]
            if not demos:
                raise RuntimeError(f"No demonstrations for phase: {phase}")
            y0 = demos[0].positions[0]
            if goal is None:
                goal = demos[0].positions[-1]
            return model.generate(y0, goal)

    def get_phase_statistics(self) -> Dict:
        """Get statistics about learned demonstrations."""
        stats = {}
        for phase in self.PHASES:
            demos = self.demonstrations[phase]
            if demos:
                stats[phase] = {
                    'n_demonstrations': len(demos),
                    'avg_duration': np.mean([d.duration for d in demos]),
                    'avg_steps': np.mean([d.num_steps for d in demos]),
                    'avg_max_force': np.mean([
                        np.max(np.linalg.norm(d.forces, axis=1)) for d in demos
                    ])
                }
        return stats
