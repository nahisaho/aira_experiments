from __future__ import annotations

from collections import deque
from pathlib import Path
from typing import Deque, Optional

import numpy as np
import yaml

import rclpy
from rclpy.node import Node
from sensor_msgs.msg import Image, Imu
from nav_msgs.msg import Odometry
from std_msgs.msg import String

try:
    from px4_msgs.msg import VehicleOdometry
except Exception:  # pragma: no cover - optional PX4 bridge dependency
    VehicleOdometry = None

from .vio_backend import (
    ImuSample,
    OnManifoldPreintegrator,
    adaptive_feature_budget,
    detect_degenerate_motion,
)


class VioNode(Node):
    """ROS2 VIO orchestration node with VINS-Fusion style control surfaces."""

    def __init__(self) -> None:
        super().__init__("vio_module")
        config_path = self.declare_parameter("config", "data/vio_config.yaml").value
        self.config = self._load_config(config_path)

        self.imu_buffer: Deque[ImuSample] = deque(maxlen=4000)
        self.left_image_stamp_ns: Optional[int] = None
        self.right_image_stamp_ns: Optional[int] = None
        self.last_track_ratio: float = 1.0
        self.last_blur_score: float = 100.0
        self.tracking_state: str = "INITIALIZING"
        self.gravity = np.array([0.0, 0.0, -self.config["g_norm"]], dtype=float)
        self.preintegrator = OnManifoldPreintegrator(self.gravity)

        self.create_subscription(Imu, self.config["imu_topic"], self._on_imu, 200)
        self.create_subscription(Image, self.config["image0_topic"], self._on_left_image, 20)
        self.create_subscription(Image, self.config["image1_topic"], self._on_right_image, 20)

        self.odom_pub = self.create_publisher(Odometry, "/vio/odometry", 20)
        self.status_pub = self.create_publisher(String, "/vio/tracking_status", 10)
        self.px4_pub = None
        if self.config.get("px4_bridge", {}).get("enabled", False) and VehicleOdometry is not None:
            self.px4_pub = self.create_publisher(
                VehicleOdometry,
                self.config["px4_bridge"]["publish_px4_topic"],
                20,
            )

        self.process_timer = self.create_timer(1.0 / self.config["camera_rate_hz"], self._process_frame)
        self.get_logger().info("VIO node initialized")

    def _load_config(self, path_str: str) -> dict:
        path = Path(path_str)
        if not path.exists():
            raise FileNotFoundError(f"Configuration file not found: {path}")
        with path.open("r", encoding="utf-8") as handle:
            raw = handle.read().replace("%YAML:1.0", "")
        return yaml.safe_load(raw)

    def _to_seconds(self, stamp) -> float:
        return float(stamp.sec) + float(stamp.nanosec) * 1e-9

    def _on_imu(self, msg: Imu) -> None:
        accel = np.array(
            [msg.linear_acceleration.x, msg.linear_acceleration.y, msg.linear_acceleration.z],
            dtype=float,
        )
        gyro = np.array([msg.angular_velocity.x, msg.angular_velocity.y, msg.angular_velocity.z], dtype=float)
        self.imu_buffer.append(ImuSample(self._to_seconds(msg.header.stamp), accel, gyro))

    def _on_left_image(self, msg: Image) -> None:
        self.left_image_stamp_ns = msg.header.stamp.sec * 1_000_000_000 + msg.header.stamp.nanosec

    def _on_right_image(self, msg: Image) -> None:
        self.right_image_stamp_ns = msg.header.stamp.sec * 1_000_000_000 + msg.header.stamp.nanosec

    def _process_frame(self) -> None:
        if self.left_image_stamp_ns is None or self.right_image_stamp_ns is None or len(self.imu_buffer) < 4:
            return

        tracker_cfg = self.config["feature_tracker"]
        texture_entropy = 0.45 if self.tracking_state != "LOST" else 0.20
        target_features = adaptive_feature_budget(
            texture_entropy=texture_entropy,
            track_ratio=self.last_track_ratio,
            blur_score=self.last_blur_score,
            min_features=tracker_cfg["min_feature_count"],
            max_features=tracker_cfg["max_feature_count"],
        )

        mean_parallax_px = 2.0 if self.tracking_state != "INITIALIZING" else 20.0
        angular_rate = float(np.linalg.norm(self.imu_buffer[-1].gyro))
        linear_speed = 0.5
        degenerate, reasons = detect_degenerate_motion(
            mean_parallax_px=mean_parallax_px,
            angular_rate_rad_s=angular_rate,
            linear_speed_m_s=linear_speed,
            hessian_condition_number=5e4 if self.tracking_state != "DEGRADED" else 2e5,
        )
        self.tracking_state = "TRACKING_DEGRADED" if degenerate else "TRACKING_GOOD"

        preintegrated = self.preintegrator.integrate(
            list(self.imu_buffer)[-20:],
            accel_bias=np.zeros(3),
            gyro_bias=np.zeros(3),
        )
        position = preintegrated.delta_p
        velocity = preintegrated.delta_v
        odom_msg = self._build_ros_odometry(position, velocity, target_features, reasons)
        self.odom_pub.publish(odom_msg)

        status = String()
        status.data = f"state={self.tracking_state}; target_features={target_features}; reasons={','.join(reasons) or 'none'}"
        self.status_pub.publish(status)

        if self.px4_pub is not None:
            self.px4_pub.publish(self._build_px4_odometry(position, velocity))

    def _build_ros_odometry(
        self,
        position: np.ndarray,
        velocity: np.ndarray,
        target_features: int,
        reasons: list[str],
    ) -> Odometry:
        msg = Odometry()
        msg.header.stamp = self.get_clock().now().to_msg()
        msg.header.frame_id = self.config["px4_bridge"]["input_frame"]
        msg.child_frame_id = self.config["px4_bridge"]["body_frame"]
        msg.pose.pose.position.x = float(position[0])
        msg.pose.pose.position.y = float(position[1])
        msg.pose.pose.position.z = float(position[2])
        msg.twist.twist.linear.x = float(velocity[0])
        msg.twist.twist.linear.y = float(velocity[1])
        msg.twist.twist.linear.z = float(velocity[2])
        covariance_scale = 0.02 if self.tracking_state == "TRACKING_GOOD" else 0.15
        msg.pose.covariance[0] = covariance_scale
        msg.pose.covariance[7] = covariance_scale
        msg.pose.covariance[14] = covariance_scale * 1.5
        msg.pose.covariance[21] = 0.05
        msg.pose.covariance[28] = 0.05
        msg.pose.covariance[35] = 0.08
        msg.twist.covariance[0] = 0.05
        msg.twist.covariance[7] = 0.05
        msg.twist.covariance[14] = 0.08
        self.last_track_ratio = min(1.0, target_features / max(self.config["feature_tracker"]["max_feature_count"], 1))
        self.last_blur_score = 80.0 if not reasons else 60.0
        return msg

    def _build_px4_odometry(self, position: np.ndarray, velocity: np.ndarray):
        msg = VehicleOdometry()
        now_us = int(self.get_clock().now().nanoseconds / 1000)
        msg.timestamp = now_us
        msg.timestamp_sample = now_us
        msg.pose_frame = VehicleOdometry.POSE_FRAME_NED
        msg.velocity_frame = VehicleOdometry.VELOCITY_FRAME_NED
        msg.position = [float(position[1]), float(position[0]), float(-position[2])]
        msg.q = [1.0, 0.0, 0.0, 0.0]
        msg.velocity = [float(velocity[1]), float(velocity[0]), float(-velocity[2])]
        msg.angular_velocity = [float(self.imu_buffer[-1].gyro[1]), float(self.imu_buffer[-1].gyro[0]), float(-self.imu_buffer[-1].gyro[2])]
        msg.position_variance = [0.02, 0.02, 0.03]
        msg.orientation_variance = [0.05, 0.05, 0.08]
        msg.velocity_variance = [0.05, 0.05, 0.08]
        msg.reset_counter = 0
        msg.quality = 100 if self.tracking_state == "TRACKING_GOOD" else 60
        return msg


def main() -> None:
    rclpy.init()
    node = VioNode()
    try:
        rclpy.spin(node)
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == "__main__":
    main()
