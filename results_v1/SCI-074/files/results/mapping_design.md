# 3D Mapping System Design for GPS-Denied Autonomous Drone

## 1. Scope and operating assumptions

This design specifies a ROS2/PX4-compatible 3D mapping subsystem for indoor and near-indoor autonomous drone flight. The primary sensor is the **Intel RealSense D455 depth stream**, spatially registered to the VIO state estimate. The mapper maintains both a **dense TSDF map** for planning / inspection and a **fallback occupancy map** for degraded compute modes.

Target roles:
- local obstacle avoidance
- traversability / free-space extraction
- inspection-grade surface reconstruction in selected regions
- map persistence and reuse across repeated flights

## 2. Primary mapper: VDBFusion

### 2.1 Why OpenVDB TSDF

VDBFusion stores a sparse truncated signed distance field (TSDF) on top of OpenVDB's hierarchical tree. For each voxel center \(\mathbf{v}\), the signed distance update is:

\[
D_t(\mathbf{v}) = \frac{W_{t-1}(\mathbf{v}) D_{t-1}(\mathbf{v}) + w_t(\mathbf{v}) d_t(\mathbf{v})}{W_{t-1}(\mathbf{v}) + w_t(\mathbf{v})}
\]

\[
W_t(\mathbf{v}) = \min\left(W_{max}, W_{t-1}(\mathbf{v}) + w_t(\mathbf{v})\right)
\]

where:
- \(d_t(\mathbf{v})\) is the projective signed distance from the current depth measurement,
- \(w_t(\mathbf{v})\) is a confidence weight derived from range, incidence angle, and pose covariance.

Compared with Octomap, TSDF provides:
- continuous zero-crossing surfaces for smoother planning and meshing
- better surface normal estimation for inspection and contact-aware control
- faster nearest-surface queries and ray marching
- natural GPU acceleration path for projective integration and raycasting

### 2.2 Integration with D455 point clouds

Recommended pipeline:
1. Acquire rectified depth or point cloud from D455 at 640×480 or 848×480.
2. Transform points into the `map` frame using the current VIO pose.
3. Reject invalid ranges and apply temporal deskewing if motion exceeds threshold.
4. Fuse points into the active VDB tiles.
5. Publish local ESDF/occupancy slices for planning.

Depth preprocessing:
- valid range: 0.3-8.0 m (navigation), 0.2-4.0 m (inspection)
- voxel downsample before fusion for bounded compute
- bilateral or temporal filter only if latency budget permits
- use confidence weight reduction for grazing incidence and depth edge pixels

### 2.3 Resolution tiers

- **5 cm voxels**: default navigation map
  - suitable for corridor flight, obstacle inflation, and real-time replanning
  - lower memory growth and faster updates
- **1 cm voxels**: inspection submaps
  - enabled only in operator-selected ROI or autonomous close-inspection mode
  - should be bounded to local submaps to avoid memory blowup

### 2.4 Confidence-weighted incremental fusion

A practical measurement weight is:

\[
w_t = w_r \cdot w_\theta \cdot w_p
\]

with:

\[
w_r = \exp\left(-\alpha_r (z-z_0)^2\right), \qquad w_\theta = \max(0, \mathbf{n}_{ray}^\top \mathbf{n}_{surf})^\gamma, \qquad w_p = \frac{1}{1 + \operatorname{tr}(\Sigma_{pose})}
\]

This down-weights distant points, grazing rays, and points fused while localization covariance is high.

Recommended truncation distances:
- navigation map: \(\mu = 3v\) to \(5v\), where \(v\) is voxel size
- inspection map: \(\mu = 6v\) for smoother meshing

### 2.5 Memory management with hierarchical VDB tiles

OpenVDB's tree structure naturally partitions active and inactive tiles:
- root node: sparse global region index
- internal nodes: mid-level tile grouping
- leaf nodes: active voxel blocks near observed surfaces

Memory policy:
- keep **active tiles** around the robot and recently observed surfaces in RAM
- mark distant stable tiles as **inactive / compressed**
- support tile-level serialization for submap streaming and multi-flight reuse
- maintain a separate LRU cache for 1 cm inspection tiles

## 3. Fallback mapper: Octomap

Octomap is maintained as a compute-light fallback and a compatibility layer for planners expecting probabilistic occupancy.

### 3.1 Occupancy update model

Octomap stores log-odds occupancy:

\[
L_t(n) = L_{t-1}(n) + \log\frac{p(n|z_t)}{1-p(n|z_t)} - L_0
\]

with clamping:

\[
L_{min} \le L_t(n) \le L_{max}
\]

Advantages:
- explicit free/occupied state through ray-casting
- mature ROS ecosystem support
- compact octree compression for large but coarse maps

Recommended modes:
- **10 cm** voxels for low-latency obstacle mapping
- **5 cm** voxels for detail mode if VDBFusion is unavailable

### 3.2 Free-space clearing

For each ray from sensor origin to measured hit:
- nodes before the hit receive free-space updates
- the terminal node receives an occupied update
- max ray length must match depth reliability envelope

This is useful for navigation when dynamic obstacles move out of the sensor frustum and free space must be explicitly cleared.

## 4. VDBFusion vs Octomap comparison

Representative values below assume indoor mapping on Jetson Orin NX with D455 depth at 640×480 / 30 Hz and a 5-10 m sensing horizon. Actual numbers depend on scene sparsity, filtering, and implementation quality.

| Feature | VDBFusion | Octomap |
|---------|-----------|---------|
| Update speed | 4-10 M points/s CPU, 15-30 M points/s with CUDA/OpenVDB-accelerated fusion | 0.5-1.5 M points/s CPU |
| Memory usage | 0.4-1.2 GB for active 100×100×10 m map at 5 cm; 1 cm limited to local ROI | 0.2-0.6 GB at 10 cm, 0.6-2.0 GB at 5 cm for similar environment |
| Raycasting | 20-60 Hz local ray marching / surface queries; continuous SDF gradients available | 5-20 Hz depending on octree depth; binary free/occupied semantics |
| GPU support | Good fit for CUDA kernels and GPU ray marching | No native GPU path in mainstream implementation |

Interpretation:
- **VDBFusion** is the better primary map for trajectory optimization, dense inspection, and smooth collision cost generation.
- **Octomap** remains valuable as a robust, low-complexity backup and for consumers that expect occupancy probabilities.

## 5. Map management strategy

### 5.1 Submapping for large environments

Use overlapping submaps rather than a monolithic global structure.

Recommended policy:
- create a new submap every 15-25 m translation or 90 deg accumulated heading change
- maintain one active local submap, one neighboring submap, and a compressed global catalog
- use loop closures from VIO to refine submap origins without rewriting every voxel immediately
- run background merge when the vehicle is landed or compute margin is available

Submap metadata should include:
- submap UUID
- origin pose in global frame
- calibration hash
- timestamp interval
- voxel size and truncation distance
- map quality metrics (coverage, covariance, density)

### 5.2 Serialization / deserialization

Serialization requirements:
- VDBFusion: save OpenVDB grids plus metadata JSON/YAML sidecar
- Octomap: save `.bt` or `.ot` files with frame metadata and covariance summary
- support partial tile/submap serialization for incremental mission checkpoints

Operational recommendations:
- save local navigation submaps every 5-10 s or on key mission events
- checkpoint inspection submaps when ROI acquisition finishes
- verify calibration hash before reloading previously saved maps

### 5.3 Map sharing between flights

For repeated operations in the same facility:
- preload the global submap index during startup
- lazy-load only nearby VDB tiles when relocalization confidence is high
- keep a shared occupancy layer for conservative obstacle memory
- tag stale areas using observation age and relocalization consistency checks

Conflict handling:
- if a previously free corridor becomes occupied repeatedly, downgrade confidence of old free-space cells
- if large map disagreement persists, fork a new environment version rather than overwriting blindly

## 6. ROS2 / PX4 integration

Subscriptions:
- `/camera/depth/color/points` or rectified depth image
- `/vio/odometry`
- optional `/semantic/labels` for class-aware fusion

Publications:
- `/mapping/vdb/local_grid`
- `/mapping/octomap`
- `/planning/esdf_slice`
- `/planning/occupied_cloud`
- `/mapping/status`

Service interfaces:
- `/mapping/save_submap`
- `/mapping/load_submap`
- `/mapping/clear_local`
- `/mapping/switch_backend`

PX4 itself typically consumes obstacle abstractions indirectly through planners rather than full 3D maps, so the ROS2 planner node should transform map products into trajectory constraints or local obstacle messages.

## 7. Compute budget and deployment modes

| Mode | Backend | Voxel size | Target latency | Use case |
|---|---|---:|---:|---|
| Navigation | VDBFusion | 5 cm | 20-35 ms / depth frame | Real-time collision avoidance |
| Inspection | VDBFusion | 1 cm local ROI | 35-70 ms / ROI update | Close-range asset reconstruction |
| Fallback | Octomap | 10 cm | 10-20 ms / depth frame | Low-power degraded mode |
| Detailed fallback | Octomap | 5 cm | 20-40 ms / depth frame | Compatibility mode |

## 8. Validation plan

1. **Static wall test**: verify TSDF zero-crossing stability and surface noise.
2. **Repeat loop flight**: check submap reuse and map alignment across flights.
3. **Narrow corridor test**: measure free-space clearing correctness.
4. **Inspection hover test**: validate 1 cm local ROI fidelity and memory capping.
5. **Backend switch test**: hot-swap VDBFusion to Octomap without planner failure.

Acceptance criteria:
- local obstacle map update rate stays above 20 Hz in navigation mode
- submap load/save success rate > 99%
- no memory exhaustion during 30 min warehouse mission
- map drift between repeated flights remains within VIO global consistency bounds

## 9. Design recommendation

Use **VDBFusion as the default mapper** and **Octomap as a fallback / interoperability layer**. This combination gives:
- dense geometry where inspection and smooth planning need it,
- probabilistic occupancy when compute is constrained or legacy planners are used,
- persistent submaps that can be shared across repeated flights in GPS-denied environments.
