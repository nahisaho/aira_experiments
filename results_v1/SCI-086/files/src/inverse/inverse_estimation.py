"""
inverse_estimation.py
=====================
Module 4: Patient-specific parameter estimation via inverse problems.
Uses ECG and echocardiography data to personalize model parameters.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple, Callable
from pathlib import Path
import logging
import json

logger = logging.getLogger(__name__)


@dataclass
class ECGData:
    """12-lead ECG measurement data."""
    signals: np.ndarray        # (n_leads, n_samples)
    sampling_rate: float       # Hz
    lead_names: List[str] = field(default_factory=lambda: [
        "I", "II", "III", "aVR", "aVL", "aVF",
        "V1", "V2", "V3", "V4", "V5", "V6"
    ])

    @property
    def duration_ms(self) -> float:
        return self.signals.shape[1] / self.sampling_rate * 1000

    def get_qrs_duration(self) -> float:
        """Estimate QRS duration from lead II."""
        lead_ii = self.signals[1]
        threshold = np.max(np.abs(lead_ii)) * 0.1
        above = np.abs(lead_ii) > threshold
        if not np.any(above):
            return 100.0
        indices = np.where(above)[0]
        return (indices[-1] - indices[0]) / self.sampling_rate * 1000

    def get_qt_interval(self) -> float:
        """Estimate QT interval."""
        return self.get_qrs_duration() + 280.0  # Simplified


@dataclass
class EchoData:
    """Echocardiography measurement data."""
    edv: float              # End-diastolic volume (mL)
    esv: float              # End-systolic volume (mL)
    ef: float               # Ejection fraction (%)
    wall_thickness: Dict[str, float] = field(default_factory=dict)
    gls: float = -20.0      # Global longitudinal strain (%)
    e_prime: float = 10.0   # E' velocity (cm/s)
    e_a_ratio: float = 1.5  # E/A ratio

    @property
    def stroke_volume(self) -> float:
        return self.edv - self.esv


@dataclass
class InverseResult:
    """Result of inverse parameter estimation."""
    estimated_params: Dict[str, float]
    cost_history: List[float]
    n_iterations: int
    converged: bool
    residual: float
    sensitivity: Dict[str, float] = field(default_factory=dict)
    confidence_intervals: Dict[str, Tuple[float, float]] = field(default_factory=dict)


class ForwardECGModel:
    """
    Forward ECG computation from cardiac electrical activity.

    Uses lead field theory and pseudo-ECG computation:
    φ_ECG(r_e) = ∫ σ_i ∇V_m · ∇(1/|r-r_e|) dΩ
    """

    def __init__(self, torso_conductivity: float = 0.2):
        self.sigma_torso = torso_conductivity
        self.electrode_positions = self._standard_12lead_positions()

    def _standard_12lead_positions(self) -> Dict[str, np.ndarray]:
        """Standard 12-lead ECG electrode positions (simplified)."""
        return {
            "RA": np.array([150.0, 100.0, 0.0]),
            "LA": np.array([-150.0, 100.0, 0.0]),
            "LL": np.array([-50.0, -200.0, 0.0]),
            "V1": np.array([30.0, 50.0, 100.0]),
            "V2": np.array([10.0, 50.0, 100.0]),
            "V3": np.array([-10.0, 30.0, 100.0]),
            "V4": np.array([-40.0, 10.0, 80.0]),
            "V5": np.array([-70.0, 10.0, 60.0]),
            "V6": np.array([-100.0, 10.0, 40.0]),
        }

    def compute_pseudo_ecg(self, V_m: np.ndarray,
                            node_positions: np.ndarray
                            ) -> Dict[str, np.ndarray]:
        """
        Compute pseudo-ECG signals from transmembrane potential field.

        V_m: (n_nodes,) transmembrane potentials at current time
        node_positions: (n_nodes, 3) spatial coordinates
        """
        ecg = {}
        for name, pos in self.electrode_positions.items():
            phi = 0.0
            for i in range(len(V_m)):
                r = pos - node_positions[i]
                dist = np.linalg.norm(r)
                if dist > 1e-6:
                    phi += V_m[i] / (4 * np.pi * self.sigma_torso * dist)
            ecg[name] = phi

        # Compute standard leads
        leads = {
            "I": ecg["LA"] - ecg["RA"],
            "II": ecg["LL"] - ecg["RA"],
            "III": ecg["LL"] - ecg["LA"],
        }
        leads["aVR"] = ecg["RA"] - 0.5 * (ecg["LA"] + ecg["LL"])
        leads["aVL"] = ecg["LA"] - 0.5 * (ecg["RA"] + ecg["LL"])
        leads["aVF"] = ecg["LL"] - 0.5 * (ecg["RA"] + ecg["LA"])

        for v in ["V1", "V2", "V3", "V4", "V5", "V6"]:
            wct = (ecg["RA"] + ecg["LA"] + ecg["LL"]) / 3
            leads[v] = ecg[v] - wct

        return leads


class BayesianInverseEstimator:
    """
    Bayesian inverse problem solver for patient-specific parameter estimation.

    Uses Ensemble Kalman Inversion (EKI) or MCMC for parameter estimation.

    Estimated parameters:
    - Tissue conductivities (σ_il, σ_it)
    - Ionic model parameters (G_Na, G_CaL, etc.)
    - Mechanical parameters (a, a_f, T_ref)
    - Conduction system timing
    """

    def __init__(self, forward_model: Callable,
                 param_bounds: Dict[str, Tuple[float, float]],
                 n_ensemble: int = 50):
        self.forward_model = forward_model
        self.param_bounds = param_bounds
        self.n_ensemble = n_ensemble

    def estimate_eki(self, observations: np.ndarray,
                      obs_noise_std: float = 0.1,
                      n_iterations: int = 30
                      ) -> InverseResult:
        """
        Ensemble Kalman Inversion (EKI) for parameter estimation.

        Iglesias et al. (2013): iterative ensemble-based method.
        """
        param_names = list(self.param_bounds.keys())
        n_params = len(param_names)
        n_obs = len(observations)

        # Initialize ensemble
        rng = np.random.default_rng(42)
        ensemble = np.zeros((self.n_ensemble, n_params))
        for j, (name, (lo, hi)) in enumerate(self.param_bounds.items()):
            ensemble[:, j] = rng.uniform(lo, hi, self.n_ensemble)

        cost_history = []

        for iteration in range(n_iterations):
            # Forward model evaluation for each ensemble member
            predictions = np.zeros((self.n_ensemble, n_obs))
            for i in range(self.n_ensemble):
                params = {name: ensemble[i, j]
                         for j, name in enumerate(param_names)}
                predictions[i] = self.forward_model(params)

            # Prediction statistics
            pred_mean = predictions.mean(axis=0)
            pred_anomaly = predictions - pred_mean

            # Parameter anomaly
            param_mean = ensemble.mean(axis=0)
            param_anomaly = ensemble - param_mean

            # Cross-covariance
            C_up = (param_anomaly.T @ pred_anomaly) / (self.n_ensemble - 1)
            C_pp = (pred_anomaly.T @ pred_anomaly) / (self.n_ensemble - 1)
            C_pp += obs_noise_std**2 * np.eye(n_obs)

            # Kalman gain
            K = C_up @ np.linalg.inv(C_pp)

            # Update ensemble
            for i in range(self.n_ensemble):
                noise = rng.normal(0, obs_noise_std, n_obs)
                innovation = observations + noise - predictions[i]
                ensemble[i] += K @ innovation

            # Enforce bounds
            for j, (name, (lo, hi)) in enumerate(self.param_bounds.items()):
                ensemble[:, j] = np.clip(ensemble[:, j], lo, hi)

            # Cost
            cost = np.mean((pred_mean - observations)**2)
            cost_history.append(cost)

            if iteration % 5 == 0:
                logger.info(f"EKI iteration {iteration}: cost={cost:.6f}")

        # Final estimates
        final_params = {name: ensemble[:, j].mean()
                       for j, name in enumerate(param_names)}

        # Confidence intervals (from ensemble spread)
        ci = {}
        for j, name in enumerate(param_names):
            vals = ensemble[:, j]
            ci[name] = (np.percentile(vals, 2.5), np.percentile(vals, 97.5))

        # Sensitivity analysis
        sensitivity = self._compute_sensitivity(final_params, observations)

        converged = len(cost_history) > 1 and \
                    cost_history[-1] < cost_history[0] * 0.01

        return InverseResult(
            estimated_params=final_params,
            cost_history=cost_history,
            n_iterations=n_iterations,
            converged=converged,
            residual=cost_history[-1],
            sensitivity=sensitivity,
            confidence_intervals=ci,
        )

    def _compute_sensitivity(self, params: Dict[str, float],
                               observations: np.ndarray
                               ) -> Dict[str, float]:
        """Local sensitivity analysis via finite differences."""
        sensitivity = {}
        base_pred = self.forward_model(params)
        base_cost = np.mean((base_pred - observations)**2)

        for name, value in params.items():
            perturbed = params.copy()
            delta = abs(value) * 0.01 + 1e-8
            perturbed[name] = value + delta
            pert_pred = self.forward_model(perturbed)
            pert_cost = np.mean((pert_pred - observations)**2)
            sensitivity[name] = abs(pert_cost - base_cost) / delta

        # Normalize
        max_sens = max(sensitivity.values()) if sensitivity else 1.0
        sensitivity = {k: v / max_sens for k, v in sensitivity.items()}

        return sensitivity


class ECGInverseSolver:
    """
    ECG-based inverse problem for conduction parameters.

    Estimates:
    - Conduction velocities (CV_long, CV_trans)
    - Activation sequence timing
    - Repolarization heterogeneity (APD dispersion)
    """

    def __init__(self, mesh_nodes: np.ndarray):
        self.mesh_nodes = mesh_nodes
        self.forward_ecg = ForwardECGModel()

    def estimate_conduction_params(self, ecg_data: ECGData,
                                     initial_guess: Optional[Dict] = None
                                     ) -> InverseResult:
        """Estimate conduction parameters from 12-lead ECG."""
        if initial_guess is None:
            initial_guess = {
                "cv_long": 0.6,      # m/s
                "cv_trans": 0.2,     # m/s
                "apd_base": 280.0,   # ms
                "apd_dispersion": 30.0,  # ms
            }

        param_bounds = {
            "cv_long": (0.3, 1.2),
            "cv_trans": (0.05, 0.5),
            "apd_base": (200.0, 400.0),
            "apd_dispersion": (10.0, 80.0),
        }

        # Target features from ECG
        target_features = np.array([
            ecg_data.get_qrs_duration(),
            ecg_data.get_qt_interval(),
        ])

        def forward_fn(params):
            qrs = 120.0 / params["cv_long"]  # Simplified
            qt = qrs + params["apd_base"]
            return np.array([qrs, qt])

        estimator = BayesianInverseEstimator(
            forward_model=forward_fn,
            param_bounds=param_bounds,
            n_ensemble=30,
        )

        result = estimator.estimate_eki(target_features, obs_noise_std=5.0)
        logger.info(f"ECG inverse: estimated params={result.estimated_params}")
        return result


class MechanicsInverseSolver:
    """
    Echo-based inverse problem for mechanical parameters.

    Estimates:
    - Passive stiffness (a, a_f, b, b_f)
    - Active contractility (T_ref)
    - Windkessel parameters (R_p, C)
    """

    def estimate_from_echo(self, echo_data: EchoData,
                            initial_guess: Optional[Dict] = None
                            ) -> InverseResult:
        """Estimate mechanical parameters from echocardiography."""
        param_bounds = {
            "a": (0.01, 1.0),       # Passive isotropic (kPa)
            "a_f": (1.0, 50.0),     # Fiber stiffness (kPa)
            "T_ref": (50.0, 200.0), # Reference tension (kPa)
            "R_p": (0.3, 2.0),      # Peripheral resistance
        }

        # Target hemodynamic features
        target = np.array([
            echo_data.edv,
            echo_data.esv,
            echo_data.ef,
            echo_data.gls,
        ])

        def forward_fn(params):
            # Simplified forward model
            edv = 130.0 - 10 * params["a"]  # Stiffer → smaller
            T = params["T_ref"]
            ef = 40 + T / 5.0                # More contractile → higher EF
            ef = min(ef, 75)
            esv = edv * (1 - ef / 100)
            gls = -10 - T / 20.0             # More contractile → more strain
            return np.array([edv, esv, ef, gls])

        estimator = BayesianInverseEstimator(
            forward_model=forward_fn,
            param_bounds=param_bounds,
            n_ensemble=40,
        )

        # Normalize target for cost computation
        result = estimator.estimate_eki(target, obs_noise_std=2.0,
                                         n_iterations=25)

        logger.info(f"Mechanics inverse: T_ref={result.estimated_params.get('T_ref', 'N/A'):.1f} kPa, "
                    f"EF error={result.residual:.3f}")
        return result
