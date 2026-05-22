"""
Visual Servoing Module
========================
Stereo vision-based 3D reconstruction, needle/suture tracking,
and image-based/position-based visual servoing for suturing.
"""

import numpy as np
from typing import Optional, Tuple, List, Dict
from dataclasses import dataclass
from enum import Enum


class VisualServoMode(Enum):
    IBVS = "image_based"  # Image-Based Visual Servoing
    PBVS = "position_based"  # Position-Based Visual Servoing
    HYBRID = "hybrid"


@dataclass
class CameraIntrinsics:
    """Stereo endoscope camera parameters."""
    fx: float = 700.0
    fy: float = 700.0
    cx: float = 320.0
    cy: float = 240.0
    baseline: float = 0.005  # 5mm stereo baseline
    width: int = 640
    height: int = 480

    @property
    def K(self) -> np.ndarray:
        """Camera intrinsic matrix."""
        return np.array([
            [self.fx, 0, self.cx],
            [0, self.fy, self.cy],
            [0, 0, 1]
        ])


@dataclass
class TrackedObject:
    """Tracked surgical object state."""
    position_3d: np.ndarray      # (3,) 3D position [m]
    orientation: np.ndarray       # (4,) quaternion
    velocity_3d: np.ndarray       # (3,) velocity [m/s]
    confidence: float             # Tracking confidence [0, 1]
    pixel_coords: np.ndarray      # (2,) image coordinates [u, v]
    bounding_box: np.ndarray      # (4,) [x, y, w, h]
    timestamp: float


class StereoReconstructor:
    """
    3D reconstruction from stereo endoscope images.
    Supports both dense and sparse reconstruction.
    """

    def __init__(self, intrinsics: Optional[CameraIntrinsics] = None):
        self.intrinsics = intrinsics or CameraIntrinsics()
        self.T_cam_to_robot: np.ndarray = np.eye(4)  # Hand-eye calibration

    def set_hand_eye_transform(self, T: np.ndarray):
        """Set camera-to-robot transformation from hand-eye calibration."""
        self.T_cam_to_robot = T.copy()

    def triangulate_point(self, uv_left: np.ndarray,
                          uv_right: np.ndarray) -> np.ndarray:
        """
        Triangulate a 3D point from stereo correspondences.

        Parameters
        ----------
        uv_left : (2,) pixel coordinates in left image
        uv_right : (2,) pixel coordinates in right image

        Returns
        -------
        point_3d : (3,) 3D point in camera frame [m]
        """
        K = self.intrinsics
        disparity = uv_left[0] - uv_right[0]
        if abs(disparity) < 1e-6:
            return np.array([0, 0, 1.0])

        Z = K.fx * K.baseline / disparity
        X = (uv_left[0] - K.cx) * Z / K.fx
        Y = (uv_left[1] - K.cy) * Z / K.fy

        return np.array([X, Y, Z])

    def triangulate_points(self, pts_left: np.ndarray,
                           pts_right: np.ndarray) -> np.ndarray:
        """Triangulate multiple point correspondences."""
        n_points = len(pts_left)
        points_3d = np.zeros((n_points, 3))
        for i in range(n_points):
            points_3d[i] = self.triangulate_point(pts_left[i], pts_right[i])
        return points_3d

    def to_robot_frame(self, points_cam: np.ndarray) -> np.ndarray:
        """Transform points from camera frame to robot base frame."""
        if points_cam.ndim == 1:
            p_homo = np.append(points_cam, 1.0)
            return (self.T_cam_to_robot @ p_homo)[:3]

        ones = np.ones((len(points_cam), 1))
        p_homo = np.hstack([points_cam, ones])
        return (self.T_cam_to_robot @ p_homo.T).T[:, :3]


class NeedleTracker:
    """
    Surgical needle detection and tracking.
    Uses ellipse fitting for circular needle arc detection.
    """

    def __init__(self, needle_radius: float = 0.008):
        """
        Parameters
        ----------
        needle_radius : float
            Radius of the circular suturing needle [m]. Default 8mm.
        """
        self.needle_radius = needle_radius
        self.state: Optional[TrackedObject] = None
        self.kalman_state = np.zeros(9)  # [x,y,z, vx,vy,vz, ax,ay,az]
        self.kalman_P = np.eye(9) * 0.01
        self._initialized = False

    def detect_needle(self, image: np.ndarray,
                      depth_map: Optional[np.ndarray] = None) -> Optional[Dict]:
        """
        Detect needle in endoscope image.

        Uses color segmentation + Hough ellipse detection.
        Returns needle center, orientation, and grasp point.
        """
        # Simulated detection pipeline
        # In production: HSV segmentation → edge detection → ellipse fitting
        h, w = image.shape[:2] if image.ndim >= 2 else (480, 640)

        # Placeholder: simulate detection with known parameters
        detection = {
            'center_2d': np.array([w / 2, h / 2]),
            'major_axis': 50.0,  # pixels
            'minor_axis': 30.0,
            'angle': 45.0,  # degrees
            'confidence': 0.85,
            'tip_2d': np.array([w / 2 + 25, h / 2 - 15]),
            'grasp_point_2d': np.array([w / 2 - 20, h / 2 + 10])
        }
        return detection

    def update_tracking(self, detection: Dict, point_3d: np.ndarray,
                        timestamp: float):
        """Update Kalman filter with new detection."""
        dt = 0.033  # ~30 fps
        if self._initialized:
            dt = timestamp - self.state.timestamp if self.state else dt

        # State transition matrix
        F = np.eye(9)
        F[0:3, 3:6] = np.eye(3) * dt
        F[3:6, 6:9] = np.eye(3) * dt
        F[0:3, 6:9] = np.eye(3) * 0.5 * dt ** 2

        # Process noise
        Q = np.eye(9) * 0.001
        Q[6:9, 6:9] *= 10.0

        # Predict
        self.kalman_state = F @ self.kalman_state
        self.kalman_P = F @ self.kalman_P @ F.T + Q

        # Measurement update
        H = np.zeros((3, 9))
        H[0:3, 0:3] = np.eye(3)
        R = np.eye(3) * 0.0001

        z = point_3d
        y_resid = z - H @ self.kalman_state
        S = H @ self.kalman_P @ H.T + R
        K_gain = self.kalman_P @ H.T @ np.linalg.inv(S)

        self.kalman_state += K_gain @ y_resid
        self.kalman_P = (np.eye(9) - K_gain @ H) @ self.kalman_P

        self.state = TrackedObject(
            position_3d=self.kalman_state[:3].copy(),
            orientation=np.array([1, 0, 0, 0]),
            velocity_3d=self.kalman_state[3:6].copy(),
            confidence=detection.get('confidence', 0.5),
            pixel_coords=detection.get('center_2d', np.zeros(2)),
            bounding_box=np.zeros(4),
            timestamp=timestamp
        )
        self._initialized = True

    def predict_pose(self, dt: float) -> Optional[np.ndarray]:
        """Predict needle position at dt seconds ahead."""
        if not self._initialized:
            return None
        return (
            self.kalman_state[:3] +
            self.kalman_state[3:6] * dt +
            0.5 * self.kalman_state[6:9] * dt ** 2
        )


class SutureTracker:
    """Track suture thread using spline fitting."""

    def __init__(self):
        self.control_points: np.ndarray = np.array([])
        self.thread_length: float = 0.0

    def detect_thread(self, image: np.ndarray) -> Optional[np.ndarray]:
        """
        Detect suture thread and return fitted spline control points.
        Uses color segmentation + B-spline fitting.
        """
        # Simulated detection
        n_points = 10
        t = np.linspace(0, 1, n_points)
        # Simulated catenary-like thread shape
        points_2d = np.column_stack([
            320 + 100 * t,
            240 + 50 * np.sin(np.pi * t)
        ])
        self.control_points = points_2d
        return points_2d

    def get_tension_estimate(self) -> float:
        """Estimate thread tension from thread shape (curvature analysis)."""
        if len(self.control_points) < 3:
            return 0.0

        # Curvature from finite differences
        dx = np.diff(self.control_points[:, 0])
        dy = np.diff(self.control_points[:, 1])
        ddx = np.diff(dx)
        ddy = np.diff(dy)

        curvatures = []
        for i in range(len(ddx)):
            num = abs(dx[i] * ddy[i] - dy[i] * ddx[i])
            den = (dx[i]**2 + dy[i]**2)**1.5
            if den > 1e-10:
                curvatures.append(num / den)

        avg_curvature = np.mean(curvatures) if curvatures else 0.0
        tension = 1.0 / max(avg_curvature, 0.01)
        return min(tension, 10.0)


class VisualServoController:
    """
    Visual servo controller combining IBVS and PBVS.

    IBVS: Directly controls pixel errors (robust to calibration errors)
    PBVS: Controls 3D pose errors (decoupled, intuitive)
    """

    def __init__(self, mode: VisualServoMode = VisualServoMode.PBVS,
                 intrinsics: Optional[CameraIntrinsics] = None,
                 gain: float = 0.5):
        self.mode = mode
        self.intrinsics = intrinsics or CameraIntrinsics()
        self.gain = gain
        self.reconstructor = StereoReconstructor(self.intrinsics)
        self.needle_tracker = NeedleTracker()
        self.suture_tracker = SutureTracker()

    def compute_ibvs(self, features_current: np.ndarray,
                     features_desired: np.ndarray,
                     depth_est: float = 0.1) -> np.ndarray:
        """
        Image-Based Visual Servoing.

        Parameters
        ----------
        features_current : (N, 2) current feature points [u, v]
        features_desired : (N, 2) desired feature points [u, v]
        depth_est : float estimated depth [m]

        Returns
        -------
        velocity : (6,) Cartesian velocity command [vx,vy,vz, wx,wy,wz]
        """
        K = self.intrinsics
        n_features = len(features_current)

        # Build interaction matrix Ls
        Ls = np.zeros((2 * n_features, 6))
        for i in range(n_features):
            u, v = features_current[i]
            x = (u - K.cx) / K.fx
            y = (v - K.cy) / K.fy
            Z = depth_est

            Ls[2*i] = [-1/Z, 0, x/Z, x*y, -(1+x**2), y]
            Ls[2*i+1] = [0, -1/Z, y/Z, 1+y**2, -x*y, -x]

        # Error in image space
        error = (features_current - features_desired).flatten()

        # Velocity command: v = -lambda * Ls^+ * e
        Ls_pinv = np.linalg.pinv(Ls)
        velocity = -self.gain * Ls_pinv @ error

        return velocity

    def compute_pbvs(self, pose_current: np.ndarray,
                     pose_desired: np.ndarray) -> np.ndarray:
        """
        Position-Based Visual Servoing.

        Parameters
        ----------
        pose_current : (6,) current pose [x,y,z, rx,ry,rz]
        pose_desired : (6,) desired pose [x,y,z, rx,ry,rz]

        Returns
        -------
        velocity : (6,) Cartesian velocity command
        """
        error = pose_desired - pose_current
        velocity = self.gain * error
        return velocity

    def compute(self, image_left: np.ndarray,
                image_right: np.ndarray,
                target_features: Optional[np.ndarray] = None,
                target_pose: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Main visual servo computation.

        Returns
        -------
        velocity : (6,) Cartesian velocity command
        """
        if self.mode == VisualServoMode.IBVS:
            if target_features is None:
                raise ValueError("IBVS requires target features")

            detection = self.needle_tracker.detect_needle(image_left)
            if detection is None:
                return np.zeros(6)

            current_features = np.array([
                detection['center_2d'],
                detection['tip_2d'],
                detection['grasp_point_2d']
            ])
            return self.compute_ibvs(current_features, target_features[:3])

        elif self.mode == VisualServoMode.PBVS:
            if target_pose is None:
                raise ValueError("PBVS requires target pose")

            detection = self.needle_tracker.detect_needle(image_left)
            if detection is None:
                return np.zeros(6)

            # 3D reconstruction
            point_3d = self.reconstructor.triangulate_point(
                detection['center_2d'],
                detection['center_2d'] - np.array([5, 0])  # Simulated disparity
            )
            point_robot = self.reconstructor.to_robot_frame(point_3d)

            current_pose = np.zeros(6)
            current_pose[:3] = point_robot

            return self.compute_pbvs(current_pose, target_pose)

        else:
            # Hybrid: PBVS for translation, IBVS for fine orientation
            vel_pbvs = np.zeros(6)
            vel_ibvs = np.zeros(6)

            if target_pose is not None:
                vel_pbvs = self.compute_pbvs(
                    np.zeros(6),  # placeholder
                    target_pose
                )

            velocity = np.zeros(6)
            velocity[:3] = vel_pbvs[:3]
            velocity[3:] = vel_ibvs[3:] if np.any(vel_ibvs[3:]) else vel_pbvs[3:]
            return velocity


class HandEyeCalibrator:
    """
    Automatic hand-eye calibration (AX=XB) for stereo endoscope.
    Uses the Tsai-Lenz method.
    """

    def __init__(self):
        self.A_transforms: List[np.ndarray] = []  # Robot end-effector motions
        self.B_transforms: List[np.ndarray] = []  # Camera motions
        self.X: Optional[np.ndarray] = None  # cam-to-ee transform

    def add_measurement(self, T_robot: np.ndarray, T_camera: np.ndarray):
        """Add a pair of robot and camera pose measurements."""
        self.A_transforms.append(T_robot.copy())
        self.B_transforms.append(T_camera.copy())

    def calibrate(self) -> np.ndarray:
        """
        Solve AX=XB using Tsai-Lenz method.

        Returns
        -------
        X : (4, 4) camera-to-end-effector transformation
        """
        n = len(self.A_transforms)
        if n < 3:
            raise ValueError("Need at least 3 measurement pairs")

        # Compute relative motions
        A_rel = []
        B_rel = []
        for i in range(n - 1):
            A_rel.append(np.linalg.inv(self.A_transforms[i]) @ self.A_transforms[i+1])
            B_rel.append(np.linalg.inv(self.B_transforms[i]) @ self.B_transforms[i+1])

        # Simplified: use SVD-based solution for rotation
        # Full implementation would use the Tsai-Lenz dual quaternion method
        M = np.zeros((3 * len(A_rel), 3))
        b_vec = np.zeros(3 * len(A_rel))

        for i, (A, B) in enumerate(zip(A_rel, B_rel)):
            Ra = A[:3, :3]
            Rb = B[:3, :3]
            # (Ra - I) * Rx = Rb - I (simplified)
            M[3*i:3*i+3] = Ra - np.eye(3)
            b_vec[3*i:3*i+3] = (Rb - np.eye(3))[:, 0]

        # Least squares for rotation
        Rx, _, _, _ = np.linalg.lstsq(M, b_vec, rcond=None)

        self.X = np.eye(4)
        # In practice, use proper rotation estimation
        self.X[:3, 3] = Rx

        return self.X
