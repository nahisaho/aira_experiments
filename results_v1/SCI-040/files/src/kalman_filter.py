"""
Kalman filter and Extended Kalman Filter (EKF) for time-varying
volcanic deformation source estimation.

Tracks temporal evolution of source parameters (position, volume change rate)
accounting for process noise (magma supply variability) and measurement noise.

Supports:
  - Standard KF for linearized models
  - EKF for nonlinear Mogi/spheroid models
  - Unscented KF (UKF) for highly nonlinear cases
  - Smoothing (RTS smoother) for retrospective analysis

References:
  Segall (2013), JGR
  Fournier et al. (2009), GJI
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Callable
from dataclasses import dataclass, field


@dataclass
class KalmanState:
    """State vector and covariance for Kalman filter."""
    x: np.ndarray          # State vector
    P: np.ndarray          # State covariance
    param_names: List[str] # Parameter names for state elements
    timestamp: float = 0.0


@dataclass
class KalmanConfig:
    """Configuration for Kalman filter."""
    # State model: x_{k+1} = F * x_k + w_k, w_k ~ N(0, Q)
    # Observation model: z_k = h(x_k) + v_k, v_k ~ N(0, R)
    dt: float = 1.0                     # Time step [days]
    process_noise_dV: float = 1e4       # Volume change rate noise [m^3/day]
    process_noise_pos: float = 10.0     # Position random walk [m/day]
    process_noise_depth: float = 5.0    # Depth random walk [m/day]
    obs_x: Optional[np.ndarray] = None  # Observation coordinates
    obs_y: Optional[np.ndarray] = None
    adaptive_Q: bool = True             # Adaptive process noise
    innovation_window: int = 10         # Window for adaptive Q


def build_state_transition(n_state: int, dt: float, model: str = "constant_rate") -> np.ndarray:
    """
    Build state transition matrix F.

    State vector: [x, y, d, dV, dV_rate] for constant_rate model
    State vector: [x, y, d, dV] for random_walk model

    Parameters
    ----------
    n_state : state dimension
    dt : time step
    model : "constant_rate" or "random_walk"

    Returns
    -------
    F : (n_state, n_state) transition matrix
    """
    F = np.eye(n_state)

    if model == "constant_rate" and n_state >= 5:
        # dV += dV_rate * dt
        F[3, 4] = dt

    return F


def build_process_noise(config: KalmanConfig, n_state: int, model: str = "constant_rate") -> np.ndarray:
    """
    Build process noise covariance Q.

    Returns
    -------
    Q : (n_state, n_state) process noise matrix
    """
    Q = np.zeros((n_state, n_state))
    dt = config.dt

    # Position random walk
    Q[0, 0] = config.process_noise_pos**2 * dt
    Q[1, 1] = config.process_noise_pos**2 * dt
    Q[2, 2] = config.process_noise_depth**2 * dt

    # Volume change noise
    Q[3, 3] = config.process_noise_dV**2 * dt

    if model == "constant_rate" and n_state >= 5:
        # Rate noise (integrated random walk)
        Q[3, 3] = config.process_noise_dV**2 * dt**3 / 3
        Q[3, 4] = config.process_noise_dV**2 * dt**2 / 2
        Q[4, 3] = Q[3, 4]
        Q[4, 4] = config.process_noise_dV**2 * dt

    return Q


def mogi_observation_function(
    state: np.ndarray,
    obs_x: np.ndarray,
    obs_y: np.ndarray
) -> np.ndarray:
    """
    Nonlinear observation function h(x) for Mogi model.

    State: [x_src, y_src, d_src, dV, ...]
    Returns: predicted displacements [ux_1, uy_1, uz_1, ..., ux_N, uy_N, uz_N]
    """
    from .source_models import MogiSource, mogi_displacement

    src = MogiSource(x=state[0], y=state[1], d=state[2], dV=state[3])
    disp = mogi_displacement(obs_x, obs_y, src)
    return disp.flatten()


def compute_jacobian(
    h_func: Callable,
    state: np.ndarray,
    obs_x: np.ndarray,
    obs_y: np.ndarray,
    delta: float = 1e-4
) -> np.ndarray:
    """
    Compute Jacobian of observation function via finite differences.

    Returns
    -------
    H : (n_obs, n_state) Jacobian matrix
    """
    n_state = len(state)
    h0 = h_func(state, obs_x, obs_y)
    n_obs = len(h0)
    H = np.zeros((n_obs, n_state))

    for j in range(n_state):
        dx = np.zeros(n_state)
        dx[j] = max(abs(state[j]) * delta, delta)
        h_plus = h_func(state + dx, obs_x, obs_y)
        H[:, j] = (h_plus - h0) / dx[j]

    return H


class ExtendedKalmanFilter:
    """
    Extended Kalman Filter for time-varying volcanic source estimation.

    Tracks: [x_src, y_src, d_src, dV, dV_rate]
    """

    def __init__(self, config: KalmanConfig, model: str = "constant_rate"):
        self.config = config
        self.model = model
        self.n_state = 5 if model == "constant_rate" else 4

        self.F = build_state_transition(self.n_state, config.dt, model)
        self.Q = build_process_noise(config, self.n_state, model)

        # Storage for smoothing
        self.states_history: List[KalmanState] = []
        self.predictions_history: List[Tuple[np.ndarray, np.ndarray]] = []
        self.innovations: List[np.ndarray] = []

    def initialize(
        self,
        x0: np.ndarray,
        P0: np.ndarray,
        param_names: Optional[List[str]] = None
    ) -> KalmanState:
        """Initialize filter state."""
        if param_names is None:
            param_names = ['x_src', 'y_src', 'd_src', 'dV', 'dV_rate'][:self.n_state]

        state = KalmanState(
            x=x0.copy(),
            P=P0.copy(),
            param_names=param_names,
            timestamp=0.0
        )
        self.states_history = [state]
        return state

    def predict(self, state: KalmanState, dt: Optional[float] = None) -> KalmanState:
        """
        Prediction step: propagate state forward in time.
        """
        if dt is not None and dt != self.config.dt:
            F = build_state_transition(self.n_state, dt, self.model)
            Q = build_process_noise(
                KalmanConfig(dt=dt,
                             process_noise_dV=self.config.process_noise_dV,
                             process_noise_pos=self.config.process_noise_pos,
                             process_noise_depth=self.config.process_noise_depth),
                self.n_state, self.model
            )
        else:
            F = self.F
            Q = self.Q
            dt = self.config.dt

        x_pred = F @ state.x
        P_pred = F @ state.P @ F.T + Q

        pred_state = KalmanState(
            x=x_pred,
            P=P_pred,
            param_names=state.param_names,
            timestamp=state.timestamp + dt
        )

        self.predictions_history.append((x_pred.copy(), P_pred.copy()))
        return pred_state

    def update(
        self,
        pred_state: KalmanState,
        z: np.ndarray,
        R: np.ndarray,
        obs_x: np.ndarray,
        obs_y: np.ndarray,
        h_func: Callable = None
    ) -> KalmanState:
        """
        Update step: incorporate new observations.

        Parameters
        ----------
        pred_state : predicted state
        z : observation vector
        R : observation noise covariance
        obs_x, obs_y : observation coordinates
        h_func : observation function (default: Mogi)
        """
        if h_func is None:
            h_func = mogi_observation_function

        # Predicted observation
        z_pred = h_func(pred_state.x, obs_x, obs_y)

        # Innovation
        innovation = z - z_pred
        self.innovations.append(innovation)

        # Jacobian
        H = compute_jacobian(h_func, pred_state.x, obs_x, obs_y)

        # Innovation covariance
        S = H @ pred_state.P @ H.T + R

        # Kalman gain
        K = pred_state.P @ H.T @ np.linalg.inv(S)

        # State update
        x_upd = pred_state.x + K @ innovation
        P_upd = (np.eye(self.n_state) - K @ H) @ pred_state.P

        # Joseph form for numerical stability
        IKH = np.eye(self.n_state) - K @ H
        P_upd = IKH @ pred_state.P @ IKH.T + K @ R @ K.T

        upd_state = KalmanState(
            x=x_upd,
            P=P_upd,
            param_names=pred_state.param_names,
            timestamp=pred_state.timestamp
        )

        self.states_history.append(upd_state)
        return upd_state

    def filter_sequence(
        self,
        observations: List[np.ndarray],
        R_list: List[np.ndarray],
        obs_x: np.ndarray,
        obs_y: np.ndarray,
        timestamps: Optional[np.ndarray] = None,
        h_func: Callable = None
    ) -> List[KalmanState]:
        """
        Run EKF over a time series of observations.

        Parameters
        ----------
        observations : list of observation vectors
        R_list : list of observation covariances (or single R)
        obs_x, obs_y : station coordinates
        timestamps : observation times
        h_func : observation function

        Returns
        -------
        filtered_states : list of KalmanState
        """
        if len(self.states_history) == 0:
            raise RuntimeError("Call initialize() first")

        current = self.states_history[-1]
        results = []

        for k, z in enumerate(observations):
            R = R_list[k] if isinstance(R_list, list) else R_list

            if timestamps is not None and k > 0:
                dt = timestamps[k] - timestamps[k-1]
            else:
                dt = None

            pred = self.predict(current, dt=dt)
            current = self.update(pred, z, R, obs_x, obs_y, h_func)
            results.append(current)

            # Adaptive Q
            if self.config.adaptive_Q and k >= self.config.innovation_window:
                self._update_process_noise(k)

        return results

    def _update_process_noise(self, k: int):
        """Adaptively update process noise based on innovation sequence."""
        w = self.config.innovation_window
        recent = self.innovations[max(0, k-w):k+1]
        innov_cov = np.mean([np.outer(v, v) for v in recent], axis=0)

        # Simple scaling of Q
        H = compute_jacobian(
            mogi_observation_function,
            self.states_history[-1].x,
            self.config.obs_x, self.config.obs_y
        )
        expected_cov = H @ self.Q @ H.T
        scale = np.trace(innov_cov) / max(np.trace(expected_cov), 1e-20)
        scale = np.clip(scale, 0.1, 10.0)
        self.Q *= scale

    def rts_smoother(self) -> List[KalmanState]:
        """
        Rauch-Tung-Striebel smoother for retrospective analysis.

        Provides improved estimates by using all observations.

        Returns
        -------
        smoothed_states : list of KalmanState (reverse-time smoothed)
        """
        n = len(self.states_history)
        if n < 2:
            return self.states_history

        smoothed = [None] * n
        smoothed[-1] = self.states_history[-1]

        for k in range(n - 2, -1, -1):
            state_k = self.states_history[k]
            pred_k1, P_pred_k1 = self.predictions_history[k]

            # Smoother gain
            G = state_k.P @ self.F.T @ np.linalg.inv(P_pred_k1)

            x_s = state_k.x + G @ (smoothed[k+1].x - pred_k1)
            P_s = state_k.P + G @ (smoothed[k+1].P - P_pred_k1) @ G.T

            smoothed[k] = KalmanState(
                x=x_s,
                P=P_s,
                param_names=state_k.param_names,
                timestamp=state_k.timestamp
            )

        return smoothed


class UnscentedKalmanFilter:
    """
    Unscented Kalman Filter for highly nonlinear source models.
    Uses sigma points to propagate uncertainty through nonlinear functions.
    """

    def __init__(self, config: KalmanConfig, model: str = "constant_rate"):
        self.config = config
        self.model = model
        self.n_state = 5 if model == "constant_rate" else 4

        self.F = build_state_transition(self.n_state, config.dt, model)
        self.Q = build_process_noise(config, self.n_state, model)

        # UKF parameters
        self.alpha = 1e-3
        self.beta = 2.0
        self.kappa = 0.0
        self.lam = self.alpha**2 * (self.n_state + self.kappa) - self.n_state

        self.states_history: List[KalmanState] = []

    def _sigma_points(self, x: np.ndarray, P: np.ndarray) -> np.ndarray:
        """Generate sigma points."""
        n = len(x)
        sigma_pts = np.zeros((2*n + 1, n))
        sigma_pts[0] = x

        sqrt_P = np.linalg.cholesky((n + self.lam) * P)

        for i in range(n):
            sigma_pts[i+1] = x + sqrt_P[i]
            sigma_pts[n+i+1] = x - sqrt_P[i]

        return sigma_pts

    def _weights(self) -> Tuple[np.ndarray, np.ndarray]:
        """Compute UKF weights."""
        n = self.n_state
        Wm = np.full(2*n + 1, 1.0 / (2*(n + self.lam)))
        Wc = np.full(2*n + 1, 1.0 / (2*(n + self.lam)))
        Wm[0] = self.lam / (n + self.lam)
        Wc[0] = self.lam / (n + self.lam) + (1 - self.alpha**2 + self.beta)
        return Wm, Wc

    def initialize(self, x0: np.ndarray, P0: np.ndarray,
                   param_names: Optional[List[str]] = None) -> KalmanState:
        if param_names is None:
            param_names = ['x_src', 'y_src', 'd_src', 'dV', 'dV_rate'][:self.n_state]
        state = KalmanState(x=x0.copy(), P=P0.copy(),
                           param_names=param_names, timestamp=0.0)
        self.states_history = [state]
        return state

    def predict_update(
        self,
        state: KalmanState,
        z: np.ndarray,
        R: np.ndarray,
        obs_x: np.ndarray,
        obs_y: np.ndarray,
        h_func: Callable = None
    ) -> KalmanState:
        """Combined predict-update step using sigma points."""
        if h_func is None:
            h_func = mogi_observation_function

        Wm, Wc = self._weights()

        # Prediction
        x_pred = self.F @ state.x
        P_pred = self.F @ state.P @ self.F.T + self.Q

        # Sigma points
        sigmas = self._sigma_points(x_pred, P_pred)

        # Transform sigma points through observation function
        z_sigmas = np.array([h_func(s, obs_x, obs_y) for s in sigmas])
        z_mean = Wm @ z_sigmas

        # Innovation covariance
        n_z = len(z)
        Pzz = np.zeros((n_z, n_z))
        Pxz = np.zeros((self.n_state, n_z))

        for i in range(2*self.n_state + 1):
            dz = z_sigmas[i] - z_mean
            dx = sigmas[i] - x_pred
            Pzz += Wc[i] * np.outer(dz, dz)
            Pxz += Wc[i] * np.outer(dx, dz)

        Pzz += R

        # Kalman gain
        K = Pxz @ np.linalg.inv(Pzz)

        # Update
        x_upd = x_pred + K @ (z - z_mean)
        P_upd = P_pred - K @ Pzz @ K.T

        upd_state = KalmanState(
            x=x_upd, P=P_upd,
            param_names=state.param_names,
            timestamp=state.timestamp + self.config.dt
        )
        self.states_history.append(upd_state)
        return upd_state
