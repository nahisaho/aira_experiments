from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import yaml

import rclpy
from rclpy.node import Node
from sensor_msgs.msg import PointCloud2
from nav_msgs.msg import Odometry
from std_msgs.msg import String

from .map_backends import OctomapBackend, VdbFusionBackend

try:
    from sensor_msgs_py import point_cloud2
except Exception:  # pragma: no cover - optional helper
    point_cloud2 = None


class MappingNode(Node):
    """ROS2 mapping node exposing VDBFusion primary and Octomap fallback backends."""

    def __init__(self) -> None:
        super().__init__("mapping_module")
        config_path = self.declare_parameter("config", "data/mapping_config.yaml").value
        self.config = self._load_config(config_path)
        self.latest_pose: Optional[np.ndarray] = None
        self.latest_pose_confidence = 1.0
        self.current_backend = "vdbfusion"

        vdb_cfg = self.config["vdbfusion"]
        octo_cfg = self.config["octomap"]
        self.vdb = VdbFusionBackend(
            voxel_size=vdb_cfg["voxel_size_navigation_m"],
            truncation_distance=vdb_cfg["truncation_distance_navigation_m"],
            max_weight=vdb_cfg["max_weight"],
        )
        self.octomap = OctomapBackend(
            resolution=octo_cfg["resolution_realtime_m"],
            prob_hit=octo_cfg["prob_hit"],
            prob_miss=octo_cfg["prob_miss"],
        )

        self.create_subscription(Odometry, self.config["mapping"]["odometry_topic"], self._on_odometry, 20)
        self.create_subscription(PointCloud2, self.config["mapping"]["pointcloud_topic"], self._on_cloud, 10)
        self.status_pub = self.create_publisher(String, "/mapping/status", 10)
        self.create_timer(1.0, self._publish_status)
        self.get_logger().info("Mapping node initialized")

    def _load_config(self, path_str: str) -> dict:
        with Path(path_str).open("r", encoding="utf-8") as handle:
            return yaml.safe_load(handle)

    def _on_odometry(self, msg: Odometry) -> None:
        self.latest_pose = np.array(
            [msg.pose.pose.position.x, msg.pose.pose.position.y, msg.pose.pose.position.z],
            dtype=float,
        )
        covariance_trace = msg.pose.covariance[0] + msg.pose.covariance[7] + msg.pose.covariance[14]
        self.latest_pose_confidence = 1.0 / (1.0 + covariance_trace)

    def _on_cloud(self, msg: PointCloud2) -> None:
        if self.latest_pose is None or point_cloud2 is None:
            return
        points = np.array(
            [[p[0], p[1], p[2]] for p in point_cloud2.read_points(msg, field_names=("x", "y", "z"), skip_nans=True)],
            dtype=float,
        )
        if points.size == 0:
            return

        if self._should_fallback_to_octomap():
            self.current_backend = "octomap"
            self.octomap.insert_point_cloud(points, self.latest_pose)
        else:
            self.current_backend = "vdbfusion"
            self.vdb.integrate_points(points, self.latest_pose, pose_confidence=self.latest_pose_confidence)

    def _should_fallback_to_octomap(self) -> bool:
        switch_cfg = self.config["backend_switching"]
        if not switch_cfg["allow_runtime_switch"]:
            return False
        return self.latest_pose_confidence < 0.55 and self.config["octomap"]["use_as_fallback"]

    def save_maps(self) -> None:
        base = Path(self.config["mapping"]["serialization_directory"])
        self.vdb.serialize(str(base / "vdb_current.yaml"))
        self.octomap.serialize(str(base / "octomap_current.yaml"))

    def _publish_status(self) -> None:
        msg = String()
        occupied_points = (
            self.vdb.occupancy_points().shape[0]
            if self.current_backend == "vdbfusion"
            else self.octomap.occupied_points().shape[0]
        )
        msg.data = (
            f"backend={self.current_backend}; pose_confidence={self.latest_pose_confidence:.3f}; "
            f"occupied_points={occupied_points}"
        )
        self.status_pub.publish(msg)


def main() -> None:
    rclpy.init()
    node = MappingNode()
    try:
        rclpy.spin(node)
    finally:
        node.save_maps()
        node.destroy_node()
        rclpy.shutdown()


if __name__ == "__main__":
    main()
