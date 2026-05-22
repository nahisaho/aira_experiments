# Autonomous Flight System Architecture for GPS-Denied Environments

本設計は warehouse / industrial indoor inspection を主対象とし、GPS-denied 環境での autonomous navigation を ROS2 + PX4 + VSLAM/VIO + dynamic obstacle avoidance で実現する。

## 1. System Overview

### 1.1 Mission Objective
- Robust navigation without GNSS
- Real-time 3D perception and collision avoidance
- PX4-based flight control with ROS2 mission autonomy
- Optional inventory/asset scanning during flight
- Graceful degradation under low texture, dynamic obstacles, or sensor dropouts

### 1.2 High-Level Functional Partition
1. **Perception**: stereo/depth/LiDAR/IMU data ingestion
2. **Localization**: Visual-Inertial Odometry (VIO)
3. **Mapping**: local 3D occupancy / TSDF-style voxel fusion
4. **Scene Understanding**: obstacle detection, tracking, trajectory prediction
5. **Planning**: global mission planning + local kinodynamic replanning
6. **Control Interface**: ROS2 autonomy stack to PX4 offboard control
7. **Safety**: health monitoring, failsafe state machine, degraded mode handling
8. **Mission Payload**: barcode / QR inventory scanning

## 2. Overall System Block Diagram

```text
+----------------------------------------------------------------------------------+
|                                Ground Control Station                            |
|         QGroundControl / Mission Console / Telemetry Logger / Remote Override    |
+-----------------------------------^-------------------------------^--------------+
                                    | WiFi 6E                       | 900 MHz
                                    | high-rate telemetry/video     | low-rate backup
                                    v                               v
+----------------------------------------------------------------------------------+
|                         Companion Computer: Jetson Orin NX                       |
|                                                                                  |
|  +----------------+     +----------------+      +-----------------------------+   |
|  | Sensor Drivers | --> |   vio_node     | ---> | mapping_node                |   |
|  | RealSense/IMU  |     | VINS-Fusion    |      | VDBFusion / OctoMap fallback|   |
|  +----------------+     +----------------+      +-----------------------------+   |
|          |                       |                         |                       |
|          |                       v                         v                       |
|          |            +----------------+        +--------------------------+      |
|          +----------> | detection_node | -----> | tracking_node            |      |
| depth/RGB             | YOLOv8-nano    |        | ByteTrack                |      |
|                       +----------------+        +--------------------------+      |
|                                                           |                      |
|                                                           v                      |
|                                              +-----------------------------+     |
|                                              | prediction_node             |     |
|                                              | Social Force + LSTM         |     |
|                                              +-----------------------------+     |
|                                                           |                      |
|             +----------------------+                      v                      |
| mission --->| global_planner_node  |---------------> +----------------------+   |
|             +----------------------+                 | local_planner_node   |   |
|                                                      | EGO-Planner v2       |   |
|                                                      +----------------------+   |
|                                                                 |               |
|                                                                 v               |
|                                                     +------------------------+   |
|                                                     | flight_controller_node |   |
|                                                     | MAVROS2 / PX4 bridge   |   |
|                                                     +------------------------+   |
|                                                                 | MAVLink/RTPS  |
|                 +-------------------------+                      v               |
|                 | safety_monitor_node     |<--------------+ +-----------------+  |
|                 | health/failsafe manager |               | | Pixhawk 6X      |  |
|                 +-------------------------+               | | PX4 v1.15        |  |
|                             ^                             | +-----------------+  |
|                             |                             |         | ESC/PWM    |
|                 +-------------------------+               |         v            |
|                 | inventory_scanner_node  |---------------+    Motors/Frame      |
|                 | QR/Barcode payload      |                                       |
+----------------------------------------------------------------------------------+
```

## 3. Software Stack Layers

```text
Application Layer
 ├─ Mission logic, warehouse task orchestration, inventory scan workflow
Autonomy Layer
 ├─ Global planning, local planning, prediction, safety state machine
Perception Layer
 ├─ VIO, mapping, detection, tracking, semantic fusion
Middleware Layer
 ├─ ROS2 Humble/Jazzy, DDS (CycloneDDS or Fast DDS), tf2, rosbag2
Vehicle Interface Layer
 ├─ MAVROS2, MAVLink, PX4 uXRCE-DDS / micro-RTPS, parameter sync
Real-Time Control Layer
 ├─ PX4 commander, estimator fallback, rate control, position/attitude loops
Hardware Layer
 ├─ Jetson Orin NX, Pixhawk 6X, RealSense D455, BMI088 IMU, optional LiDAR
```

## 4. ROS2 Node Graph and Interfaces

### 4.1 Node Graph (text)

```text
/realsense_camera_node
  ├─ /camera/infra1/image_rect      ┐
  ├─ /camera/infra2/image_rect      ├─> /vio_node
  ├─ /camera/color/image_raw        ├─> /detection_node
  ├─ /camera/depth/image_rect_raw   ├─> /mapping_node
  └─ /camera/aligned_depth_to_color ┘    /inventory_scanner_node

/imu_driver_node
  └─ /imu/data -------------------------> /vio_node, /safety_monitor_node

/lidar_node (optional)
  └─ /lidar/points ----------------------> /mapping_node, /safety_monitor_node

/vio_node
  ├─ /localization/odometry -------------> /mapping_node, /local_planner_node,
  |                                        /global_planner_node, /flight_controller_node,
  |                                        /safety_monitor_node
  ├─ /localization/pose -----------------> consumers needing low-bandwidth pose
  ├─ /tf, /tf_static --------------------> entire graph
  └─ /localization/status ---------------> /safety_monitor_node

/mapping_node
  ├─ /mapping/local_cloud ---------------> /local_planner_node
  ├─ /mapping/occupancy -----------------> /local_planner_node, /global_planner_node
  ├─ /mapping/octomap -------------------> optional fallback consumers
  └─ /mapping/status --------------------> /safety_monitor_node

/detection_node
  └─ /perception/obstacles_raw ----------> /tracking_node

/tracking_node
  └─ /perception/tracks -----------------> /prediction_node, /local_planner_node,
                                           /inventory_scanner_node, /safety_monitor_node

/prediction_node
  └─ /perception/predicted_tracks -------> /local_planner_node, /safety_monitor_node

/global_planner_node
  ├─ /mission/global_path ---------------> /local_planner_node
  └─ /mission/goal_status ---------------> /safety_monitor_node, GCS

/local_planner_node
  ├─ /planning/local_trajectory ---------> /flight_controller_node
  ├─ /planning/debug_markers ------------> RViz / diagnostics
  └─ /planning/planner_status -----------> /safety_monitor_node

/flight_controller_node
  ├─ /mavros/setpoint_raw/local ---------> MAVROS2/PX4
  ├─ /vehicle/flight_state --------------> all supervisory nodes
  └─ /vehicle/controller_status ---------> /safety_monitor_node

/safety_monitor_node
  ├─ /safety/status ---------------------> all nodes / GCS
  ├─ /mission/failsafe_event ------------> /flight_controller_node, /global_planner_node
  └─ services for arming/mode aborts

/inventory_scanner_node
  ├─ /inventory/items -------------------> warehouse WMS bridge / GCS
  └─ /inventory/scan_debug -------------> operator UI
```

### 4.2 Topic/Service Matrix

| Node | Interface | Direction | Type | Rate | QoS | Purpose |
|---|---|---:|---|---:|---|---|
| `vio_node` | `/localization/odometry` | pub | `nav_msgs/Odometry` | 50 Hz | reliable, volatile, depth 10 | primary state estimate for planning/control |
| `vio_node` | `/localization/pose` | pub | `geometry_msgs/PoseWithCovarianceStamped` | 50 Hz | reliable, volatile, depth 10 | low-overhead pose output |
| `vio_node` | `/localization/path` | pub | `nav_msgs/Path` | 5 Hz | reliable, transient_local, depth 1 | visualization and debugging |
| `vio_node` | `/localization/status` | pub | `diagnostic_msgs/DiagnosticStatus` | 2 Hz | reliable, volatile, depth 5 | estimator health |
| `vio_node` | `/camera/infra1/image_rect` | sub | `sensor_msgs/Image` | 30 Hz | sensor_data (best_effort, depth 5) | left stereo |
| `vio_node` | `/camera/infra2/image_rect` | sub | `sensor_msgs/Image` | 30 Hz | sensor_data | right stereo |
| `vio_node` | `/imu/data` | sub | `sensor_msgs/Imu` | 200 Hz | sensor_data | high-rate inertial fusion |
| `mapping_node` | `/mapping/local_cloud` | pub | `sensor_msgs/PointCloud2` | 10 Hz | best_effort, volatile, depth 5 | local obstacle cloud |
| `mapping_node` | `/mapping/occupancy` | pub | `nav_msgs/OccupancyGrid` or voxel slice | 5 Hz | reliable, transient_local, depth 1 | planner cost representation |
| `mapping_node` | `/mapping/octomap` | pub | `octomap_msgs/Octomap` | 1 Hz | reliable, transient_local, depth 1 | fallback 3D occupancy map |
| `mapping_node` | `/mapping/status` | pub | `autonomy_interfaces/MappingStatus` | 1 Hz | reliable, volatile, depth 5 | map confidence / coverage |
| `mapping_node` | `/camera/depth/image_rect_raw` | sub | `sensor_msgs/Image` | 30 Hz | sensor_data | depth fusion input |
| `mapping_node` | `/localization/odometry` | sub | `nav_msgs/Odometry` | 50 Hz | reliable, volatile, depth 10 | pose for map integration |
| `mapping_node` | `/lidar/points` | sub | `sensor_msgs/PointCloud2` | 10 Hz | sensor_data | optional LiDAR fusion |
| `detection_node` | `/perception/obstacles_raw` | pub | `autonomy_interfaces/ObstacleArray` | 15 Hz | reliable, volatile, depth 10 | dynamic obstacle detections |
| `detection_node` | `/camera/color/image_raw` | sub | `sensor_msgs/Image` | 30 Hz | sensor_data | RGB detection input |
| `detection_node` | `/camera/aligned_depth_to_color` | sub | `sensor_msgs/Image` | 30 Hz | sensor_data | depth association |
| `tracking_node` | `/perception/tracks` | pub | `autonomy_interfaces/TrackedObjectArray` | 15 Hz | reliable, volatile, depth 10 | persistent tracked objects |
| `tracking_node` | `/perception/obstacles_raw` | sub | `autonomy_interfaces/ObstacleArray` | 15 Hz | reliable, volatile, depth 10 | input detections |
| `prediction_node` | `/perception/predicted_tracks` | pub | `autonomy_interfaces/TrackedObjectArray` + `visualization_msgs/MarkerArray` | 10 Hz | reliable, volatile, depth 10 | future motion hypotheses |
| `prediction_node` | `/perception/tracks` | sub | `autonomy_interfaces/TrackedObjectArray` | 15 Hz | reliable, volatile, depth 10 | tracking history |
| `local_planner_node` | `/planning/local_trajectory` | pub | `trajectory_msgs/MultiDOFJointTrajectory` | 20 Hz | reliable, volatile, depth 10 | executable short-horizon trajectory |
| `local_planner_node` | `/planning/planner_status` | pub | `diagnostic_msgs/DiagnosticStatus` | 2 Hz | reliable, volatile, depth 5 | planner state |
| `local_planner_node` | `/mapping/local_cloud` | sub | `sensor_msgs/PointCloud2` | 10 Hz | best_effort, volatile, depth 5 | obstacle map |
| `local_planner_node` | `/perception/predicted_tracks` | sub | `autonomy_interfaces/TrackedObjectArray` | 10 Hz | reliable, volatile, depth 10 | dynamic obstacle forecast |
| `local_planner_node` | `/mission/global_path` | sub | `nav_msgs/Path` | on change / 1 Hz keepalive | reliable, transient_local, depth 1 | route guidance |
| `global_planner_node` | `/mission/global_path` | pub | `nav_msgs/Path` | on change / 1 Hz keepalive | reliable, transient_local, depth 1 | mission route |
| `global_planner_node` | `/mission/waypoints` | sub | `geometry_msgs/PoseArray` | on request | reliable, transient_local, depth 1 | operator mission input |
| `flight_controller_node` | `/mavros/setpoint_raw/local` | pub | `mavros_msgs/PositionTarget` | 30 Hz | reliable, volatile, depth 10 | PX4 offboard setpoint stream |
| `flight_controller_node` | `/vehicle/flight_state` | pub | `autonomy_interfaces/FlightState` | 10 Hz | reliable, volatile, depth 10 | normalized vehicle state |
| `flight_controller_node` | `/planning/local_trajectory` | sub | `trajectory_msgs/MultiDOFJointTrajectory` | 20 Hz | reliable, volatile, depth 10 | tracking input |
| `flight_controller_node` | `/mavros/state` | sub | `mavros_msgs/State` | 10 Hz | reliable, volatile, depth 10 | PX4 mode / arm status |
| `flight_controller_node` | `/mavros/local_position/odom` | sub | `nav_msgs/Odometry` | 30 Hz | reliable, volatile, depth 10 | PX4 estimate/failsafe cross-check |
| `safety_monitor_node` | `/safety/status` | pub | `autonomy_interfaces/SafetyStatus` | 5 Hz | reliable, transient_local, depth 1 | fleet-visible system health |
| `safety_monitor_node` | `/mission/failsafe_event` | pub | `std_msgs/String` | event-driven | reliable, transient_local, depth 5 | abort, land, hold, RTH-like indoor behavior |
| `inventory_scanner_node` | `/inventory/items` | pub | `autonomy_interfaces/InventoryItem` | event-driven / up to 5 Hz | reliable, transient_local, depth 20 | scanned asset output |
| `inventory_scanner_node` | `/camera/color/image_raw` | sub | `sensor_msgs/Image` | 30 Hz | sensor_data | QR/barcode image source |

### 4.3 Recommended QoS Profiles

| Profile Name | Reliability | Durability | Depth | Use Case |
|---|---|---|---:|---|
| `sensor_data` | best_effort | volatile | 5 | images, IMU, depth, LiDAR |
| `state_estimate` | reliable | volatile | 10 | odometry, flight state, planner trajectory |
| `latched_map` | reliable | transient_local | 1 | occupancy map, global path, safety state |
| `event_control` | reliable | transient_local | 5 | mission events, failsafe triggers, inventory scan events |

## 5. PX4-MAVROS Integration Layer

### 5.1 Integration Concept
- **Companion side**: ROS2 autonomy stack on Jetson Orin NX
- **Bridge**: MAVROS2 for service/topic compatibility and PX4 Offboard setpoint streaming
- **Vehicle side**: Pixhawk 6X running PX4 v1.15
- **Transport**: MAVLink over UART/Ethernet; optional uXRCE-DDS bridge for selected PX4 topics

### 5.2 Data Flow
1. `flight_controller_node` converts local planned trajectory into `mavros_msgs/PositionTarget`
2. MAVROS2 publishes to PX4 Offboard interface at ≥30 Hz
3. PX4 executes position/velocity/yaw-rate setpoints
4. PX4 publishes mode, arming state, battery, estimator status, IMU backup state
5. `safety_monitor_node` cross-checks VIO-based odometry with PX4 estimator and link health
6. On autonomy fault, `flight_controller_node` commands HOLD / LAND / MANUAL handover

### 5.3 MAVROS2 Interfaces
- Published:
  - `/mavros/setpoint_raw/local`
  - `/mavros/vision_pose/pose`
  - `/mavros/vision_speed/speed_twist`
- Subscribed:
  - `/mavros/state`
  - `/mavros/local_position/odom`
  - `/mavros/battery`
  - `/mavros/extended_state`
  - `/mavros/imu/data`
- Services:
  - `/mavros/cmd/arming`
  - `/mavros/set_mode`
  - `/mavros/cmd/land`
  - `/mavros/param/set`

### 5.4 PX4 Estimation Strategy
- Primary navigation source for autonomy: `vio_node`
- PX4 internal estimator (EKF2) receives vision pose aiding
- If VIO quality drops below threshold, system transitions:
  1. reduce speed envelope
  2. freeze global progression
  3. attempt hover using last good pose / optical flow fallback if installed
  4. controlled landing if localization not recovered within timeout

## 6. Sensor Suite Specification

### 6.1 Baseline Sensor Set

| Sensor | Model | Placement | Rate | Role |
|---|---|---|---:|---|
| Stereo + RGB-D camera | Intel RealSense D455 | forward-facing, 10-15° down tilt | 30 Hz image / depth | VIO features, depth, obstacle detection, QR scanning |
| IMU | BMI088 | center of gravity, vibration isolated | 200-400 Hz | inertial propagation, vibration monitoring |
| Barometer | Pixhawk onboard | FC internal | 50 Hz | altitude stabilization backup |
| Magnetometer | external if needed | isolated mast | 10-50 Hz | optional yaw aiding; low trust indoors |
| Rangefinder (optional) | LightWare / Benewake | downward | 20-50 Hz | precise indoor landing / floor height |
| LiDAR (optional) | Livox Mid-360 or Ouster OS0-class | top/front mount | 10 Hz | sparse texture fallback, better occupancy coverage |

### 6.2 Sensor Placement Guidance
- Camera baseline aligned with body X-axis
- D455 mounted close to vehicle centerline to minimize rotational parallax error
- BMI088 isolated with soft dampers and rigid frame reference
- Optional LiDAR mounted high enough to minimize propeller occlusion
- QR/barcode FOV can be improved with slight gimbal or fixed 20° tilt depending on aisle geometry

### 6.3 Time Synchronization
- ROS2 timestamps from hardware driver clock
- Jetson monotonic time as system reference
- IMU-camera sync via driver-level alignment or software interpolation in VINS-Fusion
- PX4 timesync through MAVROS2

## 7. Communication Architecture

### 7.1 Middleware
- **ROS2 DDS** on companion computer and edge devices
- Recommended DDS: **CycloneDDS** for predictable discovery and lower overhead indoors
- Use separate DDS domain ID for each vehicle in multi-UAV deployments

### 7.2 micro-ROS / Embedded Extensions
- micro-ROS can be used for auxiliary embedded boards such as:
  - lighting controller
  - gripper / scanning trigger board
  - battery management monitor
  - UWB anchor/beacon interface if later added
- micro-ROS transports: UART or UDP over companion network

### 7.3 Network Partitioning
- High-bandwidth links:
  - onboard MIPI/USB3 camera to Jetson
  - optional CSI encoder to video stream to GCS via WiFi 6E
- Safety/command links:
  - MAVLink telemetry over 900 MHz
  - RC/manual override on dedicated control link
- DDS traffic classes:
  - perception topics local-only when possible
  - compressed telemetry mirrored to GCS bridge

### 7.4 Reliability Strategy
- Critical control topics use `reliable`
- Sensor topics use `best_effort`
- Safety state and mission events use `transient_local`
- Mission planner caches last valid path and safety mode locally
- Telemetry bridge throttles noncritical visualization under bandwidth pressure

## 8. Hardware Platform Specification

### 8.1 Companion Computer
- **NVIDIA Jetson Orin NX 16GB**
- 8-core Arm CPU + Ampere GPU with Tensor Cores
- CUDA / TensorRT acceleration for YOLOv8-nano and LSTM inference
- 512 GB NVMe SSD recommended for rosbag2 and model storage
- Ubuntu 22.04 LTS + ROS2 Humble (or 24.04 + Jazzy if full stack validated)

### 8.2 Flight Stack Hardware
- **Pixhawk 6X** flight controller
- PX4 v1.15
- Triple-redundant IMU on FC for flight stabilization
- Isolated power rail and separate safety switch / buzzer recommended

### 8.3 Mechanical Platform
- 450 mm wheelbase quadrotor
- 10-11 inch propellers
- 4S or 6S Li-ion/LiPo depending payload/endurance tradeoff
- Landing gear clearance sufficient for D455 and optional downward rangefinder

## 9. Detailed Node Responsibilities

### 9.1 `vio_node`
- Based on VINS-Fusion with stereo + IMU
- Outputs drift-bounded local odometry in `map` / `odom` frame
- Provides estimator confidence, feature count, reprojection residual, relocalization status

### 9.2 `mapping_node`
- Primary: VDBFusion for efficient sparse voxel integration
- Fallback: OctoMap for conservative occupancy when compute is constrained or VDB fusion unavailable
- Maintains local rolling map (10-20 m radius) and optional global stitched submaps

### 9.3 `detection_node`
- YOLOv8-nano TensorRT engine for dynamic obstacle classes: person, forklift, pallet jack, drone, cart
- Uses depth alignment for 3D centroid and approximate size

### 9.4 `tracking_node`
- ByteTrack multi-object tracking
- Maintains stable IDs, velocities, occlusion timers, confidence history

### 9.5 `prediction_node`
- Hybrid predictor:
  - Social Force Model for short-horizon reactive avoidance
  - LSTM for learned motion priors in aisle intersections / human traffic
- Publishes multi-hypothesis trajectories with uncertainty

### 9.6 `local_planner_node`
- EGO-Planner v2 for fast kinodynamic replanning
- Inputs static map + dynamic predictions + vehicle limits
- Outputs jerk-bounded local trajectory for 2-3 s horizon

### 9.7 `global_planner_node`
- Mission graph / waypoint sequencing
- Optional aisle graph routing and no-fly zone constraints
- Replans around blocked aisles based on map updates

### 9.8 `flight_controller_node`
- Converts trajectory to PX4-compatible position/velocity/yaw setpoints
- Handles arming, offboard engagement, takeoff, landing, hover, mission state

### 9.9 `safety_monitor_node`
- Watches VIO confidence, map freshness, planner latency, FCU link, battery, CPU/GPU temperature
- Owns escalation ladder: WARN -> DEGRADED -> HOLD -> LAND

### 9.10 `inventory_scanner_node`
- Decodes barcode/QR from RGB frames
- Associates scan with pose and shelf zone
- Publishes inventory events to WMS bridge

## 10. Failsafe and Degraded-Mode Design

### 10.1 Trigger Conditions
- VIO tracking loss > 0.5 s
- Local planner latency > 150 ms sustained
- FCU link loss > 0.5 s
- Battery below reserve threshold
- Propagated obstacle prediction uncertainty above threshold in narrow aisle
- Thermal throttling on Jetson

### 10.2 Response Hierarchy
1. Reduce max velocity and acceleration
2. Hover and hold position
3. Backtrack to last safe waypoint if map remains valid
4. Controlled vertical landing at safe spot / corridor center
5. Manual takeover via RC or GCS

## 11. Recommended Frames and Conventions
- `map`: global local-consistent frame initialized at mission start
- `odom`: continuous VIO frame without discrete jumps
- `base_link`: vehicle body frame
- `camera_link`, `camera_depth_optical_frame`, `imu_link`, `lidar_link`
- ENU used in ROS2; converted as needed for PX4/NED inside `flight_controller_node` or MAVROS2

## 12. Implementation Notes
- Use composition containers for perception nodes to reduce copy overhead
- Enable intra-process communication where practical
- Use TensorRT INT8/FP16 for `detection_node`
- Pin planner and flight interface threads to isolated CPU cores
- Record synchronized rosbag2 profiles for validation: perception-only, autonomy-debug, mission-full

## 13. Design Summary
This architecture provides a practical indoor autonomous flight stack: VIO-centric localization, sparse voxel mapping, dynamic obstacle understanding, fast local replanning, and a strict safety supervisor around PX4 Offboard control. 日本語で言えば、warehouse のような GPS が使えない環境でも、認識・計画・制御を疎結合にしつつ、failsafe を強くした production-oriented 構成である。
