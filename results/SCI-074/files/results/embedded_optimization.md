# Embedded GPU Optimization Design for GPS-Denied Autonomous Drone

## Scope and assumptions
- Platform: NVIDIA Jetson Orin NX 16GB running JetPack 6.x, ROS 2 Humble, PX4 1.14-class autopilot bridge, RealSense D455 RGB-D input, and a forward IMU synchronized to the VIO stack.
- Target operating point: 1280×800 image capture, 640×640 detection crop, 30 Hz perception loop, 20 Hz map update, 20 Hz planning, 100 Hz control.
- Design target: sustain 20–30 fps throughput with <50 ms capture-to-command latency by overlapping CPU, GPU, and DLA work across consecutive frames.
- Important caveat: direct published Orin NX numbers exist for YOLO/TensorRT and Jetson platform thermals; SuperPoint, ByteTrack, VDBFusion, and ESDF values below are **engineering allocations** cross-checked against related Jetson/Isaac/NVBlox/EGO-Planner benchmarks rather than one-to-one vendor measurements.

## Platform summary
| Resource | Value | Design implication |
|---|---:|---|
| CPU | 8× Arm Cortex-A78AE @ 2.0 GHz | Reserve isolated cores for ROS 2 executor, VIO, and planner threads |
| GPU | 1024-core Ampere @ 918 MHz | Primary accelerator for depth, mapping, ESDF, and prediction |
| DLA | 2× NVDLA v2.0 | Offload supported TensorRT subgraphs to preserve GPU headroom |
| Memory | 16 GB LPDDR5, 102.4 GB/s | Enough for onboard dense mapping if peak RSS stays <12 GB |
| Power modes | 10 W / 15 W / 25 W | Switch by mission phase and thermal state |

## Processing budget (target 30 fps pipeline)
The table below shows **per-frame active compute time** for each module. These numbers are **not fully additive** because the design uses pipelining and asynchronous CUDA/DLA execution.

| Module | CPU (ms) | GPU (ms) | DLA (ms) | Memory (MB) | Strategy |
|--------|---------:|---------:|---------:|------------:|----------|
| Image capture | 2.5 | - | - | 96 | CPU with V4L2/Isaac ROS zero-copy ring buffers |
| Feature extraction (SuperPoint) | - | 4.8 | 4.2 | 220 | DLA preferred if exported to supported TensorRT ops; GPU FP16 fallback |
| VIO optimization | 6.5 | - | - | 300 | CPU multi-threaded back-end, one big core reserved for bundle adjustment |
| Depth processing | - | 6.0 | - | 420 | GPU depth alignment, filtering, and point projection |
| YOLOv8-nano detection | - | 1.8 | 11.5 | 430 | Backbone/head on DLA, GPU handles unsupported layers, decode, and NMS |
| ByteTrack tracking | 1.2 | - | - | 80 | CPU; association cost is small once detector output is available |
| Trajectory prediction | - | 2.7 | - | 160 | GPU TensorRT inference on batched obstacle tracks |
| VDBFusion map update | - | 4.5 | - | 1024 | GPU TSDF/VDB integration with voxel hashing |
| ESDF computation | - | 3.6 | - | 768 | GPU incremental ESDF update on local planning window |
| Path planning (EGO-Planner) | 2.0 | - | - | 180 | CPU on isolated planner thread, 20 Hz nominal |
| **Total pipeline reservation** | **12.2** | **23.4** | **15.7** | **9300 peak RSS** | **Pipelined across frames; critical path 41–46 ms** |

### How to interpret the budget
- **Throughput budget**: With DLA handling detection and part of SuperPoint, GPU occupancy stays below saturation, enabling **24–30 fps** sustained operation in 25 W mode.
- **Latency budget**: Critical path is capture → feature/depth → VIO/map update → planner → control publication. Overlapped execution keeps end-to-end latency in the **41–46 ms** range.
- **Peak memory**: The 9.3 GB estimate includes TensorRT engines, ROS 2 buffers, local TSDF/ESDF map, trajectory buffers, and a 20–25% fragmentation margin. This stays below the stated **<12 GB** target.

## Module rationale
### Image capture
- D455 RGB + depth DMA buffers plus double/triple buffering typically consume tens of MB; 96 MB is a safe budget including timestamp metadata and pinned transport buffers.
- Prefer `isaac_ros_nitros_image_type` or equivalent zero-copy transport to avoid repeated host-device copies.

### SuperPoint feature extraction
- ⚠️ Direct Orin NX SuperPoint-on-DLA benchmarks are limited in public literature.
- A **4–5 ms** budget is realistic for a reduced-resolution FP16 TensorRT engine on Orin NX; DLA can be used only if the exported graph avoids unsupported ops. Otherwise keep SuperPoint on the GPU and reserve DLA for detection.

### VIO optimization
- A CPU back-end budget of **6–7 ms** is realistic for a 30 Hz stereo/RGB-D VIO stack on Cortex-A78AE when front-end features are already extracted and IMU preintegration is vectorized.
- Use 3 worker threads: front-end tracking, IMU propagation, and nonlinear optimization.

### Depth, mapping, and ESDF
- The depth + TSDF/ESDF budget is aligned with Isaac ROS/NVBlox-style embedded GPU mapping pipelines, where distance-field updates complete in **single-digit milliseconds to low tens of milliseconds** on Orin-class hardware.
- Keep the local ESDF window bounded (for example 12 m × 12 m × 4 m) to avoid quadratic growth in update time.

### Detection, tracking, and prediction
- Published YOLOv8n TensorRT numbers on Orin-class Jetsons are typically around **15–17 ms raw inference** at 640×640 in FP16/INT8 when run on the main compute path; the proposed split design uses **DLA for the heavy detector path** and only **~1–2 ms GPU time** for remaining decode/NMS work.
- ByteTrack itself is lightweight; the tracker rarely exceeds **1–2 ms** on embedded ARM CPUs when the detector already provides boxes.
- A small trajectory predictor (MLP/LSTM/transformer-lite) fits comfortably inside a **2–3 ms** FP16 TensorRT slot on Ampere.

## Optimization techniques
### 1. TensorRT optimization
- Build all neural networks as TensorRT engines with explicit batch and dynamic shape constraints.
- Default to **FP16** for SuperPoint, depth-related learned models, and trajectory prediction.
- Use **INT8** for YOLOv8-nano only after warehouse-domain calibration; otherwise FP16 may be preferable because INT8 can trade away small-label accuracy.
- Cache prebuilt engines per power mode and per image resolution to avoid runtime rebuilds.

### 2. DLA offloading
- Map `YOLOv8n` to **DLA0** and, if operator coverage allows, `SuperPoint` to **DLA1**.
- Keep fallback GPU engines available because DLA does not support every layer pattern.
- Pin DLA workloads to fixed frequencies during scan mode to avoid jitter.

### 3. CUDA streams
- Use at least three streams:
  1. `stream_depth_map`: depth filtering, point projection, VDBFusion integration
  2. `stream_esdf_predict`: ESDF update and trajectory prediction
  3. `stream_misc`: GPU fallback inference and copy kernels
- Synchronize with CUDA events instead of host blocking.

### 4. CPU-GPU pipelining
- Frame `N`: image capture + IMU propagation on CPU
- Frame `N-1`: detector/SuperPoint on DLA/GPU
- Frame `N-2`: map update + ESDF on GPU
- Frame `N-3`: planner/control publication on CPU
- This four-stage overlap is what makes a 30 fps target feasible despite a >33 ms summed compute budget.

### 5. Memory optimization
- Use **zero-copy buffers first**, unified memory second.
- Recommended policy:
  - Pinned host buffers for camera ingress and PX4 bridge traffic
  - Zero-copy ROS 2 intra-process transport for images, depth maps, and feature tensors
  - Unified memory only for low-touch shared metadata or map indices, with explicit prefetch before kernels to avoid page-fault spikes
- Preallocate CUDA pools and avoid frequent `cudaMalloc` in callbacks.

### 6. Adaptive computation
- If average GPU utilization exceeds 85% for >2 s:
  - Reduce detector rate from 30 Hz to 20 Hz
  - Drop detection crop from 640×640 to 512×512
  - Reduce ESDF update volume or voxel resolution (for example 5 cm → 7.5 cm)
- If temperature or battery pressure persists, hold mapping at 15 Hz while keeping VIO/control unchanged.

### 7. ROS 2 executor optimization
- Use a **multi-threaded executor** with distinct callback groups for control, perception, mapping, and logging.
- Pin control/VIO callbacks to isolated CPU cores and elevate to `SCHED_FIFO`/real-time priority where permitted.
- Prefer CycloneDDS or FastDDS profiles with shared-memory transport enabled for large image topics.
- Run PX4 offboard control in a dedicated high-priority node to avoid perception backpressure affecting command latency.

## Thermal management
### Conservative software thresholds
| SoC temperature | Action |
|---|---|
| < 80 °C | Normal operation |
| 80–85 °C | Fan to 100%, enable reduced logging, watch sustained load |
| 85–90 °C | Reduce detector frequency/resolution; clamp map update to 15–20 Hz |
| 90–95 °C | Switch to 15 W profile, bias toward safe hover/return behavior |
| ~100 °C | Hardware throttling may begin (platform-dependent) |
| 105 °C | Absolute protection/shutdown threshold in NVIDIA thermal guidance |

### Active cooling requirements
- A passive sink is usually insufficient for sustained 25 W autonomy workloads in enclosed drone fuselages.
- Recommended integration target:
  - finned heatsink or vapor-chamber spreader sized for **25 W continuous dissipation**
  - PWM blower/fan providing roughly **5–8 CFM** directed across fins
  - low-impedance thermal interface material and unobstructed exhaust path
- Place the Jetson module away from battery heat soak and ESC exhaust where possible.

### Power mode switching by mission phase
| Mission phase | Recommended mode | Rationale |
|---|---|---|
| Pre-flight checks / idle on pad | 10 W | Sensor health, map load, comms only |
| Takeoff / transit to aisle | 15 W | VIO + obstacle awareness with moderate compute |
| Dense shelf scanning / dynamic avoidance | 25 W | Full detector + mapping + planning stack |
| Return-to-pad with low battery or thermal pressure | 15 W | Preserve control margin while reducing thermal load |
| Final descent / dock alignment | 10 W | Precision landing with limited perception set |

## Expected benchmark results
| Metric | Expected value | Notes |
|---|---:|---|
| End-to-end latency | **41–46 ms** | Capture to control command in nominal 25 W scan mode |
| Sustained throughput | **24–30 fps** | Full stack, scene-dependent |
| Average module power | **15–20 W** | Jetson module only; platform power will be higher |
| Peak memory | **9–11 GB** | Leaves safety margin below 16 GB installed RAM |
| Control rate | **100 Hz** | Offboard control loop separated from perception |

## Implementation notes for ROS 2 / PX4
- Separate perception, mapping, and control into composable nodes so zero-copy transport works inside a single process where beneficial.
- Publish only compact products to PX4 (`trajectory_setpoint`, local obstacle cones, fail-safe flags) rather than raw point clouds.
- Use lock-free queues between barcode/inventory tasks and flight-critical perception to prevent scan bursts from delaying control.

## References
1. NVIDIA Jetson Orin NX Series Data Sheet: https://developer.nvidia.com/downloads/jetson-orin-nx-series-data-sheet
2. NVIDIA Jetson Benchmarks / MLPerf Edge references: https://developer.nvidia.com/embedded/jetson-benchmarks
3. YOLOv8 benchmarking on Jetson-class devices (MDPI): https://www.mdpi.com/2073-431X/15/2/74
4. NVIDIA Jetson Linux Developer Guide, platform power and thermal management: https://docs.nvidia.com/jetson/archives/r35.6.1/DeveloperGuide/SD/PlatformPowerAndPerformance.html
5. Jetson Orin NX and Orin Nano Thermal Design Guide download page: https://developer.nvidia.com/downloads/jetson-orin-nx-orin-nano-series-thermal-design-guide
6. Isaac ROS Visual SLAM documentation: https://nvidia-isaac-ros.github.io/repositories_and_packages/isaac_ros_visual_slam/index.html
7. nvblox paper / project: https://arxiv.org/html/2311.00626v2 and https://github.com/nvidia-isaac/nvblox
8. EGO-Planner project: https://gitee.com/iszhouxin/ego-planner
