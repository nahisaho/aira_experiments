# Hardware Specification for ROS2/PX4 Autonomous Quadrotor

本仕様は GPS-denied indoor / warehouse autonomous flight を想定した reference platform である。Priority は `compute margin`, `sensor observability`, `failsafe robustness`, `payload practicality` の4点。

## 1. Target Platform Summary

| Category | Selected Hardware | Notes |
|---|---|---|
| Compute | NVIDIA Jetson Orin NX 16GB | primary AI + ROS2 autonomy computer |
| Vision Sensor | Intel RealSense D455 | stereo + RGB + depth |
| Inertial Sensor | BMI088 IMU | companion-side precision inertial sensing |
| Flight Controller | Pixhawk 6X | PX4 v1.15 |
| Airframe | Custom 450 mm quadrotor | indoor payload-capable platform |
| Telemetry | WiFi 6E + 900 MHz telemetry | high-rate + resilient backup |
| Storage | 512 GB NVMe SSD | rosbag2, models, logs |
| Power System | 6S Li-ion/LiPo | recommended for efficiency and thrust margin |

## 2. Compute Subsystem

### 2.1 Jetson Orin NX 16GB

| Parameter | Spec / Recommendation |
|---|---|
| Module | Jetson Orin NX 16GB |
| CPU | 8-core Arm Cortex-A78AE |
| GPU | NVIDIA Ampere architecture |
| AI acceleration | Tensor Cores; TensorRT FP16/INT8 |
| Memory | 16 GB LPDDR5 |
| Storage | 512 GB NVMe SSD minimum |
| OS | Ubuntu 22.04 LTS |
| ROS2 | Humble recommended |
| Power mode | 15 W to 25 W configurable |
| Thermal | active heatsink + ducted airflow required |

### 2.2 Compute Workload Budget

| Workload | Typical Load | Notes |
|---|---:|---|
| `vio_node` (VINS-Fusion stereo+IMU) | 2-3 CPU cores | depends on feature count and image resolution |
| `mapping_node` (VDBFusion) | 1-2 CPU cores | memory-bandwidth sensitive |
| `detection_node` (YOLOv8-nano TensorRT) | 15-30% GPU | at 640x384 or 640x480, 15-25 FPS |
| `tracking_node` + `prediction_node` | <1 CPU core + small GPU/CPU | lightweight compared with detection |
| `local_planner_node` | 1 CPU core burst | latency-critical |
| ROS2 + logging + visualization | 1-2 CPU cores | variable |

**Recommendation**
- Run Jetson in 20-25 W mode for stable latency.
- Dedicate one isolated CPU core to `flight_controller_node` + safety callback group.
- Keep mean GPU utilization below 70% to preserve headroom for perception spikes.

## 3. Sensor Suite

### 3.1 Intel RealSense D455

| Parameter | Value / Use |
|---|---|
| Sensor type | Stereo RGB-D camera |
| Global shutter | Yes for stereo imagers |
| Baseline | ~95 mm |
| Depth range | practical indoor 0.4 m to 6+ m |
| Frame rate | 30 FPS nominal |
| Mounting | forward-facing, slightly downward |
| Uses | VIO, obstacle detection, depth association, inventory scanning |

**Mounting guidance**
- Mount near CG to reduce rotational motion distortion.
- Use vibration-damped but rigid bracket.
- Prefer unobstructed forward FOV wider than rotor plane.

### 3.2 BMI088 IMU

| Parameter | Value / Use |
|---|---|
| Type | 6-axis industrial IMU |
| Sample rate | 200-400 Hz recommended to companion |
| Placement | close to vehicle CG |
| Uses | VIO inertial propagation, vibration diagnostics |

**Integration note**
- Synchronize BMI088 and camera timestamps as tightly as possible.
- Retain Pixhawk IMUs for flight stabilization; BMI088 is companion-side estimation input, not a replacement for FC safety sensors.

### 3.3 Optional Sensors

| Sensor | Purpose | Decision Trigger |
|---|---|---|
| 3D LiDAR (Livox Mid-360 class) | better mapping in low-texture scenes | use if warehouse texture is insufficient for camera-only mapping |
| Downward rangefinder | precision landing / floor height | use if landing tolerance <10 cm |
| UWB anchors | position fallback | use in large repetitive spaces |
| Event camera | high-speed motion robustness | only for aggressive flight profiles |

## 4. Flight Controller Subsystem

### 4.1 Pixhawk 6X + PX4 v1.15

| Parameter | Spec |
|---|---|
| FCU | Pixhawk 6X |
| Firmware | PX4 v1.15 |
| Estimator | EKF2 with vision pose aiding |
| Interfaces | UART, CAN, I2C, PWM, Ethernet depending carrier |
| Safety | dedicated IO MCU, safety switch, buzzer, RC override |

### 4.2 FCU Role Split
- **Pixhawk 6X**: hard-real-time stabilization, actuator mixing, low-level failsafe
- **Jetson Orin NX**: high-level autonomy, perception, planning, mission logic

This split keeps the aircraft flyable even if companion AI processes fail.

## 5. Airframe Specification

### 5.1 Custom 450 mm Quadrotor

| Parameter | Target |
|---|---|
| Wheelbase | 450 mm |
| Propeller size | 10-11 inch |
| Motor KV | 320-500 KV depending 4S/6S selection |
| ESC rating | 35-45 A |
| Frame material | carbon fiber with modular payload rails |
| Landing gear | raised skid/leg design |

### 5.2 Mission Suitability
- Large enough to carry Jetson + D455 + telemetry + optional LiDAR
- Small enough for warehouse aisle navigation if prop guards or protective ducts are considered
- Better hover efficiency than sub-250 mm platforms for industrial payloads

## 6. Communication Subsystem

### 6.1 WiFi 6E

| Attribute | Role |
|---|---|
| Band | 6 GHz preferred, 5 GHz fallback |
| Use | high-rate telemetry, debug video, map streaming, software deployment |
| Throughput target | 50-200 Mbps practical local link |
| Limitation | line-of-sight and warehouse attenuation |

### 6.2 900 MHz Telemetry

| Attribute | Role |
|---|---|
| Use | low-rate command/telemetry backup |
| Payload | mode, battery, position health, abort commands |
| Advantage | better penetration and range than WiFi |
| Limitation | not suitable for raw video or map streaming |

### 6.3 Recommended Link Policy
- Use WiFi 6E as primary data plane
- Use 900 MHz telemetry as supervisory backup plane
- RC/manual takeover remains independent of both if operationally required

## 7. Power Budget Analysis

Assumed nominal power architecture: **6S battery** with regulated 5 V / 12 V rails.

### 7.1 Power Consumption Table

| Subsystem | Voltage Rail | Avg Power (W) | Peak Power (W) | Notes |
|---|---|---:|---:|---|
| Jetson Orin NX 16GB | 12-19 V input | 20 | 25 | 20-25 W performance mode |
| NVMe SSD + carrier overhead | 5/12 V | 4 | 6 | includes USB hub margin |
| Intel RealSense D455 | USB 5 V | 3.5 | 4.5 | depends on streaming mode |
| BMI088 IMU board | 5 V | 0.5 | 1 | negligible but included |
| Pixhawk 6X | 5 V | 2 | 3 | with peripherals |
| WiFi 6E module | 5 V | 4 | 7 | transmission bursts |
| 900 MHz telemetry | 5 V | 2 | 5 | TX duty cycle dependent |
| Optional LiDAR | 12 V | 8 | 15 | if installed |
| LEDs / buzzer / misc | 5 V | 1.5 | 3 | status accessories |
| Propulsion (hover) | main battery | 280 | 420 | platform dependent |

### 7.2 Total Budget

**Without LiDAR**
- Avionics + compute average: ~37.5 W
- Propulsion hover average: ~280 W
- **Total hover average:** ~317.5 W

**With LiDAR**
- Additional average: +8 W
- **Total hover average:** ~325.5 W

### 7.3 Battery Endurance Estimate

Assume:
- 6S 8000 mAh Li-ion pack
- Nominal energy ≈ 6 x 3.6 V x 8 Ah = **172.8 Wh**
- Usable energy at 80% discharge = **138.2 Wh**

Estimated endurance:
- Without LiDAR: 138.2 / 317.5 ≈ **26.1 min** ideal electrical estimate
- With real maneuvering, reserves, and nonideal efficiency: **14-18 min operational mission time**

**Engineering recommendation**
- Plan missions around 10-12 min usable sortie time with 25-30% reserve for indoor autonomy.

## 8. Weight Budget

### 8.1 Estimated Mass Breakdown

| Component | Mass (g) |
|---|---:|
| 450 mm frame | 650 |
| 4x motors | 320 |
| 4x ESCs | 160 |
| Propellers | 60 |
| Pixhawk 6X + GPS-less accessories | 120 |
| Jetson Orin NX + carrier + heatsink | 280 |
| NVMe SSD | 20 |
| Intel RealSense D455 | 95 |
| BMI088 IMU module | 15 |
| WiFi 6E + antennas | 35 |
| 900 MHz telemetry + antennas | 40 |
| Power distribution / regulators / wiring | 180 |
| Landing gear / mounts / vibration isolation | 180 |
| 6S 8000 mAh battery | 900 |
| Payload margin / protective cage / scanner accessories | 250 |
| **Estimated AUW without LiDAR** | **3305 g** |
| Optional LiDAR | 250 |
| **Estimated AUW with LiDAR** | **3555 g** |

### 8.2 Thrust Margin Requirement
- Minimum recommended thrust-to-weight ratio: **>2.0:1** for safe maneuvering indoors
- Target total thrust:
  - without LiDAR: >6.6 kgf
  - with LiDAR: >7.1 kgf
- Prefer motor/prop combination that can provide 2.0-2.3 kgf per motor at acceptable efficiency

## 9. Thermal and EMC Considerations

### 9.1 Thermal
- Jetson requires forced airflow; passive cooling alone is insufficient in enclosed warehouses.
- Keep intake/exhaust clear of prop wash turbulence zones that cause recirculation.
- Set thermal monitoring alarms in `safety_monitor_node` at 80-85°C SoC temperature.

### 9.2 EMI / Vibration
- Separate telemetry antennas from high-current ESC leads.
- Route camera and USB cables away from power distribution board.
- Use balanced props and soft-mounted IMU/camera assembly.

## 10. Recommended Electrical Architecture

```text
Main 6S Battery
 ├─ Power Distribution Board
 │   ├─ 4x ESC -> Motors
 │   ├─ 5V BEC -> Pixhawk 6X + telemetry
 │   ├─ 12V/19V regulator -> Jetson carrier
 │   └─ 5V regulator -> RealSense / accessories
 └─ Current/voltage sensing -> Pixhawk power module
```

## 11. Integration Checklist
- Confirm D455 extrinsics to `base_link`
- Measure BMI088-to-camera transform and time offset
- Validate Jetson boot behavior under voltage sag
- Tune PX4 for heavier companion-compute mass distribution
- Verify UART/Ethernet bandwidth between Jetson and Pixhawk
- Flight-test with prop guards if operating near shelving or people

## 12. Recommended Procurement / Build Notes
- Prefer industrial USB locking connectors for D455
- Use conformal coating or dust protection if operating in warehouse with particulate matter
- Keep spare battery, props, and telemetry antennas in maintenance kit
- Reserve mounting points for optional LiDAR even if not installed on Rev-A

## 13. Summary
この hardware specification は、indoor autonomous flight に必要な perception-heavy workload を Jetson Orin NX で処理しつつ、Pixhawk 6X に real-time stabilization を分担させる設計である。Power と weight の両面で 450 mm class は妥当で、optional LiDAR を含めても industrial prototype として十分現実的な構成になっている。
