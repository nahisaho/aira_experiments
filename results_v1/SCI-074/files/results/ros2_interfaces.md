# ROS2 Custom Interfaces for GPS-Denied Autonomous Flight

本ドキュメントは custom interface package `autonomy_interfaces` を想定している。実装時には `msg/` と `srv/` に分離し、`rosidl_default_generators` で build する。

## 1. Package Layout

```text
autonomy_interfaces/
├── msg/
│   ├── Obstacle.msg
│   ├── ObstacleArray.msg
│   ├── TrackedObject.msg
│   ├── TrackedObjectArray.msg
│   ├── FlightState.msg
│   ├── InventoryItem.msg
│   ├── MappingStatus.msg
│   └── SafetyStatus.msg
└── srv/
    ├── StartMission.srv
    ├── PauseMission.srv
    ├── ResumeMission.srv
    ├── AbortMission.srv
    ├── SetSafetyMode.srv
    ├── ClearMapRegion.srv
    └── ScanShelf.srv
```

> Note: the user requested `TrackedObject.msg`; however, for topic transport a `TrackedObjectArray.msg` wrapper is strongly recommended. `ObstacleArray.msg` similarly wraps a reusable `Obstacle.msg` element.

## 2. Message Definitions

### 2.1 `Obstacle.msg`

```text
std_msgs/Header header
string id
uint16 class_id
string class_label
float32 confidence
geometry_msgs/Pose pose
geometry_msgs/Vector3 size
geometry_msgs/Vector3 velocity
float32 distance
bool is_dynamic
float32[36] pose_covariance
```

**Field intent**
- `id`: detector-local ID if available
- `class_id` / `class_label`: semantic class from YOLOv8-nano
- `pose`: 3D obstacle centroid in `map` or `odom`
- `size`: approximate bounding box dimensions in meters
- `velocity`: estimated linear velocity in m/s
- `distance`: range from vehicle centroid
- `is_dynamic`: static vs dynamic classification hint

### 2.2 `ObstacleArray.msg`

```text
std_msgs/Header header
string frame_id
uint32 sequence_id
Obstacle[] obstacles
```

**Usage**
- Published by `detection_node`
- Consumed by `tracking_node`, `safety_monitor_node`
- Typical rate: 10-15 Hz

### 2.3 `TrackedObject.msg`

```text
std_msgs/Header header
string track_id
uint16 class_id
string class_label
float32 confidence
geometry_msgs/Pose pose
geometry_msgs/Vector3 size
geometry_msgs/Vector3 velocity
geometry_msgs/Accel acceleration
float32 yaw
float32 yaw_rate
builtin_interfaces/Duration track_age
uint32 hit_streak
uint32 miss_count
bool occluded
float32 existence_probability
geometry_msgs/Pose[] predicted_poses
float32[] prediction_time_offsets
float32[] covariance_diagonal
```

**Usage**
- Core state for `tracking_node` and `prediction_node`
- `predicted_poses` stores short-horizon forecast samples
- `covariance_diagonal` is flattened `[x,y,z,vx,vy,vz,...]` uncertainty summary

### 2.4 `TrackedObjectArray.msg`

```text
std_msgs/Header header
string frame_id
TrackedObject[] objects
```

### 2.5 `FlightState.msg`

```text
std_msgs/Header header
string vehicle_id
string nav_state
string autonomy_mode
bool armed
bool offboard_enabled
bool failsafe_active
geometry_msgs/Pose pose
geometry_msgs/Twist twist
geometry_msgs/Vector3 body_accel
float32 roll
float32 pitch
float32 yaw
float32 battery_voltage
float32 battery_current
float32 battery_remaining
float32 estimated_time_remaining
uint8 localization_quality
uint8 map_quality
uint8 link_quality
string active_controller
string last_warning
```

**Enumerations (recommended constants in code)**
- `nav_state`: `IDLE`, `TAKEOFF`, `MISSION`, `HOLD`, `LAND`, `ABORT`, `MANUAL`
- `autonomy_mode`: `STANDBY`, `ASSISTED`, `FULLY_AUTONOMOUS`, `DEGRADED`
- quality fields: `0=UNKNOWN`, `1=BAD`, `2=FAIR`, `3=GOOD`, `4=EXCELLENT`

### 2.6 `InventoryItem.msg`

```text
std_msgs/Header header
string item_id
string code_type
string raw_code
string decoded_text
float32 confidence
geometry_msgs/Pose item_pose
geometry_msgs/Pose vehicle_pose
string zone_id
string shelf_id
string mission_id
builtin_interfaces/Time scan_time
sensor_msgs/RegionOfInterest image_roi
string image_frame_id
bool verified
string notes
```

**Usage**
- Published by `inventory_scanner_node`
- Can be bridged to WMS / ERP integration layer

### 2.7 `MappingStatus.msg`

```text
std_msgs/Header header
string map_frame
float32 local_map_radius
float32 explored_volume_m3
float32 occupied_ratio
float32 map_resolution
float32 map_latency_ms
float32 integration_rate_hz
bool lidar_fusion_enabled
bool octomap_fallback_active
bool loop_closure_available
uint8 health
string status_text
```

**Health values**
- `0=UNKNOWN`
- `1=DEGRADED`
- `2=NOMINAL`
- `3=RECOVERING`

### 2.8 `SafetyStatus.msg`

```text
std_msgs/Header header
string safety_mode
bool safe_to_arm
bool safe_to_fly
bool localization_ok
bool mapping_ok
bool planning_ok
bool control_link_ok
bool battery_ok
bool thermal_ok
bool obstacle_margin_ok
float32 minimum_obstacle_distance
float32 vio_tracking_score
float32 planner_cycle_ms
float32 cpu_load
float32 gpu_load
float32 soc_temperature_c
string active_failsafe
string recommended_action
string[] active_faults
```

**Recommended `safety_mode` values**
- `NORMAL`
- `WARN`
- `DEGRADED`
- `HOLD`
- `LAND_NOW`
- `MANUAL_REQUIRED`

## 3. Mission Control Services

### 3.1 `StartMission.srv`

```text
# Request
string mission_id
geometry_msgs/PoseArray waypoints
float32 cruise_speed
float32 max_acceleration
bool enable_inventory_scan
bool enable_dynamic_avoidance
---
# Response
bool accepted
string message
builtin_interfaces/Time accepted_at
```

### 3.2 `PauseMission.srv`

```text
# Request
string mission_id
string reason
---
# Response
bool success
string message
```

### 3.3 `ResumeMission.srv`

```text
# Request
string mission_id
---
# Response
bool success
string message
```

### 3.4 `AbortMission.srv`

```text
# Request
string mission_id
string abort_mode    # HOLD, LAND, MANUAL_HANDOVER
string reason
---
# Response
bool success
string message
```

### 3.5 `SetSafetyMode.srv`

```text
# Request
string requested_mode
string rationale
---
# Response
bool success
string active_mode
string message
```

### 3.6 `ClearMapRegion.srv`

```text
# Request
geometry_msgs/Pose center
geometry_msgs/Vector3 size
string reason
---
# Response
bool success
float32 cleared_volume_m3
string message
```

### 3.7 `ScanShelf.srv`

```text
# Request
string zone_id
string shelf_id
bool hover_during_scan
float32 timeout_s
---
# Response
bool success
uint32 items_detected
string message
```

## 4. Interface Mapping to Nodes

| Node | Publishes | Calls / Offers Services |
|---|---|---|
| `detection_node` | `ObstacleArray` | - |
| `tracking_node` | `TrackedObjectArray` | - |
| `prediction_node` | `TrackedObjectArray` with forecasts | - |
| `flight_controller_node` | `FlightState` | consumes mission-control services indirectly |
| `inventory_scanner_node` | `InventoryItem` | offers `ScanShelf` |
| `mapping_node` | `MappingStatus` | offers `ClearMapRegion` |
| `safety_monitor_node` | `SafetyStatus` | offers `SetSafetyMode`; may gate mission services |
| `global_planner_node` | - | offers `StartMission`, `PauseMission`, `ResumeMission`, `AbortMission` |

## 5. Design Rationale
- `FlightState.msg` normalizes PX4 + autonomy stack state into one supervisory message.
- `SafetyStatus.msg` is intentionally explicit so GCS and onboard nodes can make deterministic decisions.
- `ObstacleArray.msg` and `TrackedObject.msg` separate detection from temporal estimation.
- `InventoryItem.msg` binds scanned code to vehicle pose for warehouse traceability.
- Mission services are idempotent where possible, which simplifies GCS retries over unreliable links.

## 6. Versioning Guidance
- Package name: `autonomy_interfaces`
- Semantic versioning recommended: `0.x` during prototyping, `1.x` once mission API stabilizes
- Reserve new fields by appending only; avoid breaking field order changes after deployment

## 7. Example IDL Dependency List

```cmake
find_package(rosidl_default_generators REQUIRED)
find_package(std_msgs REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(builtin_interfaces REQUIRED)
```

日本語メモ: interface は最初から「現場運用」で壊れにくい粒度にしておくのが重要であり、特に `SafetyStatus` と `FlightState` は later-stage integration cost を大きく下げる。
