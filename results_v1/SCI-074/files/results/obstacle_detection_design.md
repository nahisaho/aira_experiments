# Dynamic Obstacle Detection, Tracking, and Prediction Design

## Scope and assumptions
- Platform: GPS-denied autonomous multirotor using ROS 2, PX4, Jetson Orin NX, stereo camera or Intel RealSense D455, and a local map/planning stack.
- Operational environment: indoor warehouse with aisles, shelves, forklifts, pallets, carts, and pedestrians.
- Primary body frame: `base_link`; world frame: `map`; PX4 control frame bridge via `px4_ros_com` / `microRTPS`.
- Target outputs: 3D obstacle tracks with uncertainty, short/medium/long horizon trajectory predictions, and zone-based safety alerts.

## 1. System architecture

### 1.1 ROS 2 node graph
1. `camera_driver_node`
   - Publishes `/camera/left/image_rect`, `/camera/right/image_rect`, `/camera/depth/image`, `/camera/*/camera_info`.
2. `image_preproc_node`
   - Rectification, timestamp alignment, CLAHE/gamma correction when illumination degrades.
3. `yolov8_detector_node`
   - TensorRT INT8 inference of YOLOv8-nano at 30 Hz.
   - Publishes `/perception/detections_2d`.
4. `depth_localization_node`
   - Fuses detection ROIs with stereo/D455 depth to produce 3D detections.
   - Publishes `/perception/detections_3d`.
5. `bytetrack_3d_node`
   - Multi-object tracking with score-stratified association and 3D constant-acceleration Kalman filtering.
   - Publishes `/perception/tracks`.
6. `trajectory_prediction_node`
   - Generates 0-10 s predictions and covariance envelopes.
   - Publishes `/perception/predicted_trajectories`.
7. `safety_zone_monitor_node`
   - Computes zone occupancy, time-to-collision, and collision probability.
   - Publishes `/safety/dynamic_obstacle_alerts` and stop/replan triggers.

### 1.2 Key interfaces
| Topic | Type | Rate | Purpose |
|---|---|---:|---|
| `/perception/detections_2d` | custom `Detection2DArray` | 30 Hz | Class, bbox, confidence |
| `/perception/detections_3d` | custom `Detection3DArray` | 30 Hz | 3D centroid, extent, covariance |
| `/perception/tracks` | custom `TrackArray` | 30 Hz | Tracked state, ID, lifecycle |
| `/perception/predicted_trajectories` | custom `PredictedTrajectoryArray` | 10-20 Hz | Mean/covariance over horizon |
| `/safety/dynamic_obstacle_alerts` | custom `ObstacleAlertArray` | 20 Hz | Zone, TTC, collision probability |
| `/planning/replan_trigger` | `std_msgs/Bool` | event | Active replanning request |
| `/px4/collision_halt` | bridge message/service | event | Emergency halt / hover command |

Recommended QoS: `SensorDataQoS` for raw images/depth, `KeepLast(5)` reliable for tracks and alerts.

## 2. Detection pipeline

### 2.1 Detector selection and deployment
- Model: YOLOv8-nano trained/fine-tuned on warehouse classes.
- Classes: `person`, `forklift`, `pallet`, `cart`, `shelf`, `unknown_dynamic`.
- Runtime target: 30 Hz on Jetson Orin NX using TensorRT INT8.
- Detection range: 0.3-15 m.
- Confidence thresholds:
  - High-score set: `s >= 0.5`
  - Low-score set: `0.1 <= s < 0.5`
  - Suppress below `0.1`
- Input resolution recommendation: `640x384` or `640x480` depending on camera FOV and latency budget.

Expected pipeline latency budget:
- Image acquisition and rectification: 6-8 ms
- TensorRT inference: 10-14 ms
- NMS and decoding: 1-2 ms
- Depth projection and 3D fusion: 4-6 ms
- Total perception latency: 21-30 ms

### 2.2 Training and calibration
- Start from COCO-pretrained YOLOv8-nano, then fine-tune with warehouse-specific data.
- Add synthetic augmentation for pallets, forklifts, shelving occlusion, motion blur, low light, reflective floors.
- Perform INT8 calibration with representative indoor frames covering 0.3-15 m and all classes.
- Maintain per-camera extrinsic calibration `T_base_camera` and stereo rectification parameters.

### 2.3 2D-to-3D localization
For each 2D detection with bounding box center `(u, v)` and depth `z`, compute the 3D point in the camera frame:

```math
\mathbf{p}_c = z K^{-1} \begin{bmatrix}u \\ v \\ 1\end{bmatrix}
```

where `K` is the camera intrinsic matrix. Transform to the map frame:

```math
\mathbf{p}_m = T_{map}^{base} T_{base}^{camera} \mathbf{p}_c
```

Depth estimation policy:
1. Extract robust depth from the lower 40% of the ROI for ground-contact objects (`person`, `forklift`, `pallet`, `cart`).
2. Use median or trimmed mean over valid depth pixels to reject holes/outliers.
3. Use ROI depth variance to estimate measurement covariance.
4. For shelves, estimate front-face plane by RANSAC on depth points and output oriented extent.
5. If valid depth ratio `< 0.3`, mark detection as range-uncertain and keep only 2D evidence for one cycle.

Measurement vector:

```math
\mathbf{z}_k = [x, y, z, w, h, d]^T
```

Approximate covariance:

```math
R_k = \operatorname{diag}(\sigma_x^2, \sigma_y^2, \sigma_z^2, \sigma_w^2, \sigma_h^2, \sigma_d^2)
```

with `\sigma_z` increasing quadratically with range for passive stereo and piecewise for D455 beyond ~6 m.

### 2.4 Detection filtering
Apply the following gates before tracking:
- Range gate: `0.3 <= ||p|| <= 15 m`
- Height gate using warehouse priors
- Ground consistency check for floor-supported objects
- Temporal persistence check for `unknown_dynamic`
- Non-maximum suppression in image space plus 3D duplicate suppression

## 3. Multi-object tracking with ByteTrack

## 3.1 Tracker state
Use a 3D constant-acceleration Kalman filter per track:

```math
\mathbf{x}_k = [x, y, z, v_x, v_y, v_z, a_x, a_y, a_z]^T
```

Discrete-time transition for sample period `\Delta t`:

```math
\mathbf{x}_{k+1} = F(\Delta t)\mathbf{x}_k + \mathbf{w}_k
```

with block form:

```math
F(\Delta t) = \begin{bmatrix}
I_3 & \Delta t I_3 & \tfrac{1}{2}\Delta t^2 I_3 \\
0 & I_3 & \Delta t I_3 \\
0 & 0 & I_3
\end{bmatrix}
```

and process noise:

```math
Q = q \begin{bmatrix}
\tfrac{\Delta t^5}{20}I_3 & \tfrac{\Delta t^4}{8}I_3 & \tfrac{\Delta t^3}{6}I_3 \\
\tfrac{\Delta t^4}{8}I_3 & \tfrac{\Delta t^3}{3}I_3 & \tfrac{\Delta t^2}{2}I_3 \\
\tfrac{\Delta t^3}{6}I_3 & \tfrac{\Delta t^2}{2}I_3 & \Delta t I_3
\end{bmatrix}
```

Measurement model for localized detections:

```math
\mathbf{z}_k = H\mathbf{x}_k + \mathbf{v}_k, \quad H = [I_3\; 0\; 0]
```

### 3.2 ByteTrack association flow
At frame `k`:
1. Predict all active tracks with the Kalman filter.
2. Split detections into high-score set `D_h` and low-score set `D_l`.
3. First association: match active tracks to `D_h` using Hungarian assignment.
4. Second association: unmatched confirmed tracks matched to `D_l` to recover weak/partially occluded observations.
5. Tentative-track association: newly initialized tracks matched after confirmed tracks.
6. Unmatched high-score detections spawn tentative tracks.

Association cost should combine 2D and 3D consistency:

```math
C_{ij} = \lambda_{iou}(1 - IoU_{ij}) + \lambda_m d_M(\hat{\mathbf{p}}_i, \mathbf{z}_j) + \lambda_c \mathbb{1}[class_i \neq class_j]
```

where the Mahalanobis distance is

```math
d_M^2 = (\mathbf{z} - H\hat{\mathbf{x}})^T S^{-1}(\mathbf{z} - H\hat{\mathbf{x}}), \quad S = HPH^T + R
```

Recommended initial weights: `\lambda_{iou}=0.35`, `\lambda_m=0.6`, `\lambda_c=0.05`.

### 3.3 Track lifecycle management
- **Initialization**: spawn tentative track from unmatched high-score detection.
- **Confirmation**: require `N_init = 3` hits within `M = 5` frames.
- **Deletion**:
  - tentative tracks deleted after `2` consecutive misses;
  - confirmed tracks deleted after `T_lost = 15` frames or if covariance grows beyond limit.
- **State labels**: `tentative`, `confirmed`, `occluded`, `deleted`.
- **Covariance inflation** during occlusion to reflect rising uncertainty.

### 3.4 Re-identification during occlusion
For occluded objects, add lightweight appearance embeddings from a small ReID head or feature extractor attached to the detector backbone.
- Maintain an exponential moving average feature vector per track.
- When an unmatched detection appears after occlusion, perform gated association using:
  - motion gate in 3D,
  - class compatibility,
  - cosine similarity of appearance embeddings.

Cosine similarity:

```math
s_{reid}(\mathbf{f}_i, \mathbf{f}_j) = \frac{\mathbf{f}_i^T \mathbf{f}_j}{\|\mathbf{f}_i\|\|\mathbf{f}_j\|}
```

Only allow re-identification if `d_M^2 < \chi^2_{3,0.99}` and `s_reid > 0.75`.

### 3.5 Performance targets and evaluation
- MOTA `> 70%`
- ID switches `< 5%` of total trajectories or `< 5` per standard benchmark sequence
- MOTP and HOTA should also be tracked for engineering validation
- Runtime target: tracking update `< 5 ms/frame`

Validation datasets should include:
- warehouse forklift-pedestrian interaction,
- narrow aisle crossings,
- heavy shelf occlusion,
- reflective floor artifacts,
- low-light loading zones.

## 4. Trajectory prediction

### 4.1 Horizon partitioning
- **Short-term (0-2 s):** Kalman extrapolation
- **Medium-term (2-5 s):** Social Force Model (pedestrians) or constant-turn-rate/acceleration for vehicles
- **Long-term (5-10 s):** LSTM-based sequence model, TensorRT optimized

The predictor should output at `10-20 Hz` with a prediction sample step of `0.2 s`.

### 4.2 Short-term prediction: Kalman extrapolation
For tracks with state mean `\hat{x}_k` and covariance `P_k`:

```math
\hat{x}_{k+n|k} = F^n \hat{x}_k
```

```math
P_{k+n|k} = F^n P_k (F^n)^T + \sum_{i=0}^{n-1} F^i Q (F^i)^T
```

This mode is robust for immediate collision checking and emergency response.

### 4.3 Medium-term prediction: Social Force Model
For pedestrians, model acceleration as:

```math
m_i \frac{d\mathbf{v}_i}{dt} = m_i \frac{\mathbf{v}_i^0 - \mathbf{v}_i}{\tau_i} + \sum_j \mathbf{f}_{ij} + \sum_o \mathbf{f}_{io}
```

where:
- `\mathbf{v}_i^0` is desired velocity,
- `\tau_i` is relaxation time,
- `\mathbf{f}_{ij}` captures repulsion from other people/robots,
- `\mathbf{f}_{io}` captures repulsion from shelves, walls, pallets.

Use the local ESDF and aisle topology to infer candidate walking direction and constrain desired velocity along traversable corridors.

For forklifts/carts, a kinematic model is preferable:

```math
\dot{x}=v\cos\psi, \quad \dot{y}=v\sin\psi, \quad \dot{\psi}=\omega, \quad \dot{v}=a
```

### 4.4 Long-term prediction: LSTM model
Input sequence per track:
- past `T_h = 2-3 s` states,
- class embedding,
- local map context (aisle direction, obstacles),
- nearby agent encodings.

Output:
- future means `\mu_t` for `5-10 s`,
- covariance or mixture parameters.

Recommended deployment:
- Pre-train offline on warehouse trajectory logs.
- Export ONNX, optimize with TensorRT FP16/INT8.
- Run only on confirmed dynamic tracks within 10 m to preserve compute budget.

### 4.5 Uncertainty propagation and collision probability
Represent each predicted obstacle state at time `t` as Gaussian:

```math
\mathbf{o}_t \sim \mathcal{N}(\mu_{o,t}, \Sigma_{o,t})
```

and the drone planned state as:

```math
\mathbf{r}_t \sim \mathcal{N}(\mu_{r,t}, \Sigma_{r,t})
```

Relative state:

```math
\mathbf{d}_t = \mathbf{r}_t - \mathbf{o}_t \sim \mathcal{N}(\mu_{d,t}, \Sigma_{d,t}), \quad \Sigma_{d,t}=\Sigma_{r,t}+\Sigma_{o,t}
```

Collision likelihood can be approximated through Mahalanobis distance to an inflated safety ellipsoid:

```math
D_t^2 = \mu_{d,t}^T \Sigma_{d,t}^{-1} \mu_{d,t}
```

If `D_t^2 <= \chi^2_{3,\alpha}`, the overlap probability exceeds confidence level `\alpha`.
A practical collision risk score is:

```math
P_{coll}(t) \approx 1 - F_{\chi^2_3}(D_t^2)
```

where `F_{\chi^2_3}` is the CDF of the `\chi^2` distribution with 3 DoF.

### 4.6 Prediction arbitration
Use a mode switch by horizon and class:
- all objects: Kalman 0-2 s,
- pedestrians: Social Force 2-5 s, LSTM 5-10 s,
- forklifts/carts: kinematic extrapolation 2-5 s, LSTM only if dataset quality supports it,
- shelves/pallets: static unless track velocity exceeds threshold.

## 5. Safety zones and alerting logic

### 5.1 Zone definitions
- **Emergency stop zone:** `0-1.5 m` → immediate halt/hover/brake
- **Avoidance zone:** `1.5-5 m` → active replanning and speed reduction
- **Monitoring zone:** `5-15 m` → track and predict only

Distance should be evaluated against the predicted minimum separation over the planning horizon, not only current range.

### 5.2 Decision logic
For each predicted obstacle track, compute:
- current range `d_0`,
- predicted minimum distance `d_min`,
- minimum time-to-collision `TTC_min`,
- peak collision probability `P_coll^max`.

Trigger hierarchy:
1. **Emergency halt** if `d_min < 1.5 m` or `TTC_min < 0.7 s` or `P_coll^max > 0.6`.
2. **Forced replan** if `1.5 <= d_min < 5 m` or `P_coll^max > 0.2`.
3. **Monitor** if `5 <= d_min < 15 m`.

Use hysteresis to avoid chattering:
- exit emergency only when `d_min > 2.0 m` for `0.5 s`,
- exit avoidance only when `d_min > 5.5 m` for `1.0 s`.

### 5.3 Unknown dynamic objects
`unknown_dynamic` is a safety-first class.
- Inflate covariance and object dimensions.
- Force conservative velocity prior.
- Treat as dynamic for planning until `N_static = 20` consecutive frames indicate near-zero motion.

## 6. Implementation plan

### 6.1 ROS 2 package breakdown
- `perception_msgs`: `Detection2D`, `Detection3D`, `Track`, `PredictedTrajectory`, `ObstacleAlert`
- `warehouse_detector`: TensorRT YOLOv8 wrapper
- `depth_localizer`: ROI depth fusion and 3D covariance estimation
- `dynamic_tracker`: ByteTrack + Kalman + ReID
- `trajectory_predictor`: Kalman/Social Force/LSTM ensemble
- `dynamic_safety_monitor`: zone logic and planner/PX4 interfacing

### 6.2 Execution model
- Detection path pinned to dedicated CUDA stream.
- Tracker and predictor on separate callback groups under a multithreaded executor.
- Use intra-process communications or loaned messages to reduce copies.
- Synchronize via timestamp tolerances `< 20 ms`.

### 6.3 Parameterization
Critical tunables:
- detector thresholds,
- depth validity ratio,
- Kalman process noise by class,
- ByteTrack score thresholds and lost age,
- Social Force parameters `\tau`, repulsion gains,
- collision probability thresholds,
- zone hysteresis times.

### 6.4 Verification tests
1. **Unit tests**: projection math, covariance propagation, assignment logic.
2. **Playback tests**: rosbag warehouse sequences with labeled obstacles.
3. **HIL/SITL**: Gazebo/Ignition or Isaac Sim with PX4 SITL.
4. **Flight envelope tests**: staged indoor runs with safety tether and soft geofence.

Acceptance criteria:
- 30 Hz perception sustained on Orin NX,
- detection recall sufficient at 0.3-15 m,
- MOTA > 70%, ID switches < 5%,
- emergency stop latency from sensor timestamp to PX4 command < 120 ms.

## 7. Failure modes and mitigations
| Failure mode | Mitigation |
|---|---|
| Depth dropout on reflective surfaces | ROI robust statistics, stereo fallback, covariance inflation |
| Detector misses due to motion blur | shorter exposure, blur augmentation, low-score recovery in ByteTrack |
| ID switch after occlusion | 3D gating + ReID + class prior |
| Overconfident long-horizon predictions | covariance growth floor, mode arbitration, planner chance constraints |
| Perception overload | class-conditional prediction budget, skip long-horizon model for distant/static tracks |

## 8. Recommended engineering defaults
- Detector input: `640x384`, batch size `1`, INT8
- Tracking update: 30 Hz
- Prediction update: 15 Hz
- Max simultaneously predicted tracks: `20`
- Emergency object inflation radius: `0.4-0.6 m` beyond physical extent
- Warehouse evaluation KPI set: recall, precision, MOTA, HOTA, IDF1, mean TTC warning lead time
