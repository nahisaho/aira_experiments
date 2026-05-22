# Visual-Inertial Odometry Design for GPS-Denied Autonomous Drone

## 1. Scope and assumptions

This document specifies a ROS2/PX4-oriented VIO subsystem derived from the VINS-Fusion architecture for an autonomous drone operating without GNSS. The reference sensor suite is:

- **Stereo visual input**: Intel RealSense D455 left/right imagers
- **Inertial input**: BMI088 IMU
- **Target compute**: Jetson Orin NX
- **Primary output**: drift-bounded pose, velocity, and covariance for PX4 external vision fusion
- **Assumptions**:
  1. The D455 stereo pair is hardware-synchronized internally.
  2. The BMI088 is rigidly mounted to the camera frame with offline extrinsic calibration as an initial prior.
  3. Flight envelope includes indoor corridors, warehouses, and inspection spaces with intermittent low-texture regions.
  4. Optional LiDAR constraints are available through a separate lidar-inertial frontend but are not mandatory for the nominal pipeline.

## 2. System architecture

The module is organized into five stages:

1. **Sensor ingestion and time alignment**
2. **Visual frontend** (feature extraction + KLT tracking + outlier rejection)
3. **IMU pre-integration** on the manifold between keyframes
4. **Sliding-window nonlinear optimization** with marginalization
5. **Loop closure / map reuse / optional cross-modal constraints**

A ROS2 deployment should expose the following interfaces:

- Subscriptions:
  - `/camera/infra1/image_rect_raw`
  - `/camera/infra2/image_rect_raw`
  - `/imu/data_raw`
  - optional `/lidar/odometry`
- Publications:
  - `/vio/odometry`
  - `/vio/path`
  - `/vio/tracking_status`
  - `/fmu/in/vehicle_visual_odometry` (PX4 bridge)
- Services/actions:
  - `/vio/reset`
  - `/vio/save_map`
  - `/vio/load_map`

## 3. Core VINS-Fusion pipeline

### 3.1 Front-end

#### Feature extraction: SuperPoint vs ORB

| Criterion | SuperPoint | ORB |
|---|---:|---:|
| Repeatability under viewpoint change | High | Medium |
| Illumination robustness | High | Medium |
| Descriptor type | Learned float/binary projection dependent | Binary BRIEF |
| Embedded compute cost | Higher (GPU/NPU preferred) | Low |
| Memory bandwidth | Moderate-high | Low |
| Cold-start latency on Orin NX | 8-14 ms/frame (TensorRT) | 2-5 ms/frame |
| Recommended use | Primary for difficult indoor scenes | Fallback / low-power mode |

**Design choice**:
- Default to **SuperPoint** when GPU resources are available and perception load is below 70%.
- Fall back to **ORB** when thermal throttling, power-saving mode, or CPU-only execution is required.
- Use **KLT optical flow** after initial detection to reduce per-frame descriptor extraction cost.

#### Tracking

The tracker uses pyramidal Lucas-Kanade with forward-backward consistency checking:

\[
\delta \mathbf{u} = \arg\min_{\delta \mathbf{u}} \sum_{\mathbf{x}\in\Omega} \left(I_t(\mathbf{x}+\delta \mathbf{u}) - I_{t-1}(\mathbf{x})\right)^2
\]

Tracking robustness enhancements:
- image pyramid levels: 4-5
- forward-backward threshold: 0.5-1.0 px
- epipolar gating for stereo matches
- RANSAC essential/fundamental matrix filtering
- track age weighting for keyframe selection

### 3.2 IMU pre-integration on manifold

Between states \(i\) and \(j\), integrate IMU measurements in the local tangent space to avoid repeated reintegration during optimization. Using bias-corrected accelerometer and gyroscope signals:

\[
\hat{\omega}_k = \omega_k - b^g_i - n^g_k, \qquad \hat{a}_k = a_k - b^a_i - n^a_k
\]

The pre-integrated increments are:

\[
\Delta R_{ij} = \prod_{k=i}^{j-1} \exp\left(\hat{\omega}_k \Delta t_k\right)
\]

\[
\Delta v_{ij} = \sum_{k=i}^{j-1} \Delta R_{ik}\hat{a}_k \Delta t_k
\]

\[
\Delta p_{ij} = \sum_{k=i}^{j-1} \left(\Delta v_{ik}\Delta t_k + \frac{1}{2}\Delta R_{ik}\hat{a}_k \Delta t_k^2\right)
\]

The residual terms injected into the optimizer are:

\[
\mathbf{r}_{\Delta R} = \log\left((\Delta \tilde{R}_{ij}\,\exp(J^g_{\Delta R}\delta b_g))^\top R_i^\top R_j\right)
\]

\[
\mathbf{r}_{\Delta v} = R_i^\top(v_j-v_i-g\Delta t_{ij}) - \left(\Delta \tilde{v}_{ij}+J^g_{\Delta v}\delta b_g+J^a_{\Delta v}\delta b_a\right)
\]

\[
\mathbf{r}_{\Delta p} = R_i^\top\left(p_j-p_i-v_i\Delta t_{ij}-\frac{1}{2}g\Delta t_{ij}^2\right) - \left(\Delta \tilde{p}_{ij}+J^g_{\Delta p}\delta b_g+J^a_{\Delta p}\delta b_a\right)
\]

### 3.3 Back-end sliding window optimization

A window of 10-20 keyframes is optimized with Ceres Solver. The state vector is:

\[
\mathcal{X} = \{x_k, \lambda_l, t_d, b^g_k, b^a_k\}
\]

where \(x_k = (R_k, p_k, v_k)\), \(\lambda_l\) are inverse depths, and \(t_d\) is the camera-IMU temporal offset.

The cost function is:

\[
\min_{\mathcal{X}} \; \|r_p-H_p\mathcal{X}\|^2 + \sum_{(i,j)} \|r^{imu}_{ij}\|_{P^{-1}_{ij}}^2 + \sum_{l,j} \rho\left(\|r^{proj}_{l,j}\|_{\Sigma^{-1}_{l,j}}^2\right) + \sum_c \rho\left(\|r^{loop}_c\|_{\Sigma^{-1}_c}^2\right)
\]

Implementation choices:
- robust loss: Huber for reprojection, Cauchy for loop constraints
- marginalization: Schur-complement based prior on oldest keyframe
- solver budget: 8-12 iterations / 20-25 ms maximum per optimization cycle
- asynchronous loop closure thread to avoid front-end stalls

### 3.4 Loop closure with DBoW3

The loop detector uses bag-of-words retrieval followed by geometric verification:

1. Query database with current keyframe descriptor vector.
2. Reject candidates using temporal exclusion and covisibility consistency.
3. Estimate relative pose by PnP or 3D-3D alignment.
4. Accept loop only if inlier ratio, reprojection RMSE, and pose covariance pass thresholds.

Recommended thresholds:
- minimum BoW score margin: 1.15× next-best candidate
- RANSAC inlier count: > 40
- reprojection RMSE: < 1.5 px
- yaw disagreement with IMU prior: < 20 deg unless relocalization mode is active

### 3.5 Multi-camera support

For wider field of view, support additional forward-downward or fisheye cameras. The measurement model for camera \(c\) is:

\[
\mathbf{z}_{l,c,j} = \pi\left(T_{BC_c}^{-1} T_{WB_j}^{-1} P_l\right)
\]

with camera-specific extrinsics \(T_{BC_c}\). Multi-camera operation improves observability during banked turns, wall-following, and texture-sparse forward motion. Recommended policy:
- keep one **primary stereo pair** for scale and robust initialization
- add **downward monocular** or **fisheye side cameras** only if calibration quality remains < 0.5 px reprojection RMSE

## 4. Precision improvement techniques

### 4.1 Adaptive feature management

Use scene texture and track survivability to adapt the feature budget:

\[
N_f = \operatorname{clip}\left(N_{min} + k_H \hat{H}_I + k_T r_{track} - k_B r_{blur},\; N_{min},\; N_{max}\right)
\]

where:
- \(\hat{H}_I\): normalized image gradient entropy
- \(r_{track}\): fraction of successfully propagated tracks
- \(r_{blur}\): blur score from Laplacian variance inversion

Recommended policy:
- low texture: raise target features from 180 to 320 and relax spatial non-max distance
- high texture: cap at 180-220 to limit optimization density
- reserve at least 30% of features in image periphery for yaw/pitch observability

### 4.2 IMU-camera temporal calibration

Estimate an online time offset \(t_d\) jointly in the optimizer. The projection residual becomes:

\[
\mathbf{r}^{proj}_{l,j}(t_d) = \mathbf{z}_{l,j} - \pi\left(T_{CB} T_{BW}(t_j + t_d) P_l\right)
\]

Guidelines:
- initialize \(t_d\) from hardware timestamp characterization
- bound online correction to ±20 ms in nominal mode
- freeze temporal calibration after convergence variance stays below \(0.1\,\text{ms}^2\) for 100 frames

### 4.3 Robust initialization

Use a staged transition:

1. **Vision-only bootstrap** with 5-point geometry / stereo depth
2. Estimate up-to-scale motion and gravity direction hypothesis
3. Solve for gyroscope bias via short batch alignment
4. Recover scale, velocity, and accelerometer bias from inertial consistency
5. Switch to tightly coupled VI optimization once covariance drops below threshold

Transition criteria:
- median parallax > 18 px
- tracked features > 120
- IMU excitation includes both translational and rotational motion for at least 1.5 s

### 4.4 Degenerate motion handling

Detect near-singular motion from parallax, Hessian condition number, and motion type:

\[
\kappa(H_{tri}) = \frac{\sigma_{max}}{\sigma_{min}}
\]

Trigger degeneration handling when any of the following persist for > 0.5 s:
- mean parallax < 1.5 px
- \(\kappa(H_{tri}) > 10^5\)
- pure rotation score \(\|\omega\| / (\|v\| + \epsilon) > \tau\)

Mitigations:
- delay keyframe insertion
- increase IMU weight and reduce inverse-depth updates
- hold scale state from previous covariance-consistent estimate
- request downward / side camera measurements if available
- if LiDAR constraints exist, temporarily increase their factor weight

### 4.5 Map reuse for repeated flights

Persistent map reuse reduces relocalization time and drift in repeated inspection routes.

Design:
- save keyframe poses, descriptors, landmarks, and covisibility graph
- store map UUID, sensor calibration hash, and environment version
- on new mission, perform place recognition before full estimator reset
- if relocalization succeeds, inject a pose prior instead of hard state overwrite

Recommended acceptance checks for reuse:
- calibration hash exact match
- average descriptor similarity above threshold
- geometric verification inlier ratio > 0.35
- relocalized pose covariance below mission-specific bound

### 4.6 Optional LiDAR-inertial constraint integration

LiDAR odometry or scan-to-map alignment can be fused as an external factor:

\[
\mathbf{r}^{lidar}_{ij} = \log\left((T^{L}_{ij})^{-1} (T_{LB}^{-1} T_{B_i}^{-1} T_{B_j} T_{LB})\right)
\]

Use cases:
- glossy or low-light environments where visual tracks collapse
- long corridors with limited visual parallax
- dust/fog where IMU-only drift must be bounded by another geometry source

Fusion policy:
- enable only when LiDAR deskewing and extrinsics are validated
- down-weight LiDAR factors during rapid vibration or scan deformation events
- keep fusion optional to preserve low-SWaP deployments

## 5. PX4 / ROS2 integration design

### 5.1 State publication

Publish VIO in ENU on ROS2 topics and bridge to PX4 NED using a deterministic frame adapter:

- ROS2 estimator frame: `map -> base_link`
- PX4 accepted external vision frame: local NED body-aligned or FRD after conversion
- timestamp source: monotonic synchronized clock with transport delay compensation

### 5.2 Failsafe interaction

Expose estimator health flags:
- tracking state: `INITIALIZING`, `TRACKING_GOOD`, `TRACKING_DEGRADED`, `RELOCALIZING`, `LOST`
- covariance thresholds for PX4 fusion enable/disable
- timeout-driven invalidation if state age > 80 ms

Recommended PX4 behavior:
- fuse external vision yaw/position only when covariance and innovation gates pass
- drop to attitude control / hover contingency if VIO is lost for > 0.5 s indoors

## 6. Performance targets

### 6.1 Primary metrics

- **ATE target**: < 0.5% of trajectory length
- **RPE target**: < 0.01 m per meter traveled
- **Frame compute budget**: < 30 ms per frame on Jetson Orin NX
- **End-to-end latency**: < 50 ms sensor-to-PX4 publication
- **Relocalization latency with map reuse**: < 2.0 s for previously mapped areas

### 6.2 Representative embedded budget

| Pipeline stage | Target latency (ms/frame) | Notes |
|---|---:|---|
| Image rectification / sync | 1-2 | Prefer ISP or CUDA path |
| Feature detection / refresh | 3-10 | ORB at low end, SuperPoint at high end |
| KLT tracking + outlier rejection | 2-4 | CPU SIMD/CUDA optical flow |
| IMU preintegration | 0.5-1.5 | Runs at IMU rate, amortized |
| Sliding-window optimization | 10-14 | 10-frame window, capped iterations |
| Loop closure thread | async | No front-end stall budget |
| ROS2/PX4 serialization | 1-2 | Includes frame conversion |
| **Total** | **17.5-29.5** | Meets 30 ms target |

## 7. Comparative positioning

The following values are representative engineering ranges from public literature, open-source implementations, and embedded deployment experience; exact values depend on sensor quality, dataset, and optimization budget.

| System | Visual frontend | Inertial coupling | Loop closure | Multi-camera | Typical ATE (% length) | Embedded latency on Orin NX | Notes |
|---|---|---|---|---|---:|---:|---|
| **VINS-Fusion** | Point features + KLT | Tight, sliding window | Yes | Stereo / multi-cam capable | 0.3-0.8 | 18-30 ms | Best trade-off for extensibility and map reuse |
| **ORB-SLAM3** | ORB-only | Tight, MAP-based | Yes | Mono/stereo/RGB-D | 0.3-0.7 | 20-35 ms | Excellent place recognition; weaker in severe motion blur |
| **OKVIS2** | Keypoint + patch optimization | Tight, optimization-centric | Limited / external | Stereo-first | 0.4-0.9 | 22-38 ms | Strong estimator rigor; less turnkey loop/map reuse |
| **Basalt** | Direct + optical flow hybrid | Tight, spline-capable | No native global loop stack | Stereo / multi | 0.2-0.6 | 14-24 ms | Very strong local accuracy, but more integration work for persistent maps |

**Selection rationale**:
- choose **VINS-Fusion-derived architecture** when persistent mapping, loop closure, and ROS2/PX4 integration simplicity matter most;
- choose **Basalt** if the sole goal is local odometric accuracy with a tightly controlled sensor stack;
- choose **ORB-SLAM3** if place recognition dominates and ORB-only operation is acceptable;
- choose **OKVIS2** where estimator transparency is preferred over mature map management.

## 8. Verification and acceptance tests

1. **Static bias test**: 5 min stationary IMU run, Allan variance fit, verify configured noise densities.
2. **Excitation test**: handheld figure-eight motion to confirm robust initialization.
3. **Corridor test**: repeated flights with long forward motion to evaluate degenerate handling.
4. **Low-light test**: compare SuperPoint and ORB fallback latency / track survival.
5. **Map reuse test**: load previous map and measure time-to-first-valid-pose.
6. **PX4 flight test**: verify innovation gates and estimator failover logic.

Acceptance criteria:
- ATE and RPE targets met on at least three representative indoor routes
- no estimator divergence under 2 s pure rotation event
- relocalization success rate > 95% in repeated mapped flights
- CPU+GPU utilization leaves > 20% headroom for planning and control

## 9. Implementation guidance

- Keep IMU integration and visual tracking on separate callback groups.
- Pin feature extraction threads to performance cores on Orin NX.
- Use lock-free queues for image/IMU buffering.
- Save estimator priors and pose graph under mission-tagged directories.
- Gate loop closure updates so they never perturb PX4 output with discontinuous jumps; apply global corrections to map frame while maintaining smooth local odometry output.
