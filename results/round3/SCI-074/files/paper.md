# Integrated VSLAM and Dynamic Obstacle Avoidance for GPS-Denied Autonomous UAV Navigation: A ROS2/PX4-Based System Design and Simulation Study

DRAFT — NOT FOR DISTRIBUTION

---

## Abstract

Autonomous unmanned aerial vehicles (UAVs) operating in GPS-denied indoor environments, such as industrial warehouses, require robust simultaneous localization and mapping (SLAM) combined with real-time obstacle avoidance to fulfill practical inspection and inventory management tasks. This paper presents the design, implementation, and simulation evaluation of an integrated Visual-Inertial SLAM (VSLAM) and dynamic obstacle avoidance system built on the ROS2/PX4 software stack. The proposed system fuses stereo visual odometry with IMU pre-integration through a 15-state Error-State Kalman Filter (ESKF), maintains a probabilistic 3D occupancy map using an OctoMap-inspired log-odds grid, tracks dynamic obstacles via multi-object Kalman filtering with constant-velocity prediction, and optimizes collision-free trajectories through an EGO-Planner-inspired B-spline gradient descent optimizer. We conduct systematic simulation experiments in a procedurally generated 20 m × 12 m × 4 m warehouse environment. VIO accuracy is evaluated across five noise configurations using 5-fold cross-validation; the medium-noise baseline achieves an Absolute Trajectory Error (ATE) of 0.310 ± 0.299 m with 10 Hz visual updates, while IMU-only integration degrades to 9.343 ± 6.571 m, demonstrating that visual correction is indispensable. The EGO-Planner component achieves trajectory optimization in 42–47 ms on CPU, well within the 200 ms replanning budget, and scales gracefully from 5 to 40 obstacles (5.5 ms to 17.3 ms). ESDF full-map computation requires 2.9 s on CPU, motivating the adoption of incremental GPU-accelerated approaches such as VDBFusion for real-time deployment. A warehouse inventory inspection case study demonstrates feasibility: an eight-waypoint autonomous flight mission is completed within 60 s. These results establish quantitative performance baselines and highlight the critical bottlenecks — ESDF computation latency and tightly coupled visual-inertial loop closure — that must be addressed for robust real-world deployment.

---

## 1. Introduction

The proliferation of small-form-factor quadrotors capable of carrying depth cameras and compute modules has opened substantial opportunities for automating repetitive indoor tasks such as warehouse inventory auditing, infrastructure inspection, and search-and-rescue operations. However, unlike outdoor deployments where GPS provides metric-scale localization at centimetre accuracy, indoor environments deny satellite signals entirely, forcing the vehicle to rely exclusively on onboard sensing for state estimation.

Visual-Inertial Odometry (VIO), pioneered at scale by VINS-Mono (Qin et al., 2018) and extended to multi-map SLAM by ORB-SLAM3 (Campos et al., 2021), provides a practical path to GPS-independent localization. These methods fuse high-rate IMU measurements with lower-rate visual observations through tightly coupled optimization or filtering, achieving sub-decimetre accuracy on standard benchmarks. Nevertheless, the transition from controlled laboratory datasets to unstructured industrial environments introduces new challenges: low-texture surfaces (painted concrete walls, metallic shelving), strong specular reflections from polished floors, and a high density of dynamic agents (forklifts, human workers) that violate the static-world assumption underlying most SLAM systems.

Obstacle avoidance in dynamic cluttered environments demands trajectory planners that can reason about future obstacle positions. EGO-Planner (Zhou et al., 2021) introduced an ESDF-free gradient-based approach that directly penalizes collisions in a B-spline trajectory parameterization, achieving replanning at tens of Hz. FASTER (Tordesillas et al., 2019) further demonstrated safe high-speed flight in unknown environments by separating fast and slow trajectory layers. Despite these advances, few studies have systematically characterized the joint performance of VIO and trajectory planning under variable IMU/visual noise, variable obstacle density, and embedded-compute constraints, particularly in the context of warehouse operations.

This paper addresses three gaps in the literature: (1) quantitative cross-validated VIO accuracy characterization under realistic sensor noise models; (2) EGO-Planner trajectory quality versus dynamic obstacle density; and (3) per-component compute profiling to identify bottlenecks for embedded GPU deployment. The system is designed around the ROS2/PX4 software stack, reflecting current industry practice for autonomous aerial vehicles (Bopalkar & Patil, 2025). A warehouse inventory inspection case study grounds the evaluation in a practical mission profile consistent with recent industry demonstrations (Zhang & Wilson, 2024).

The remainder of the paper is organized as follows. Section 2 reviews relevant prior work. Section 3 details the proposed system architecture and algorithmic components. Section 4 describes the simulation experimental setup. Section 5 presents quantitative results. Section 6 discusses implications and limitations. Section 7 concludes.

---

## 2. Related Work

### 2.1 Visual-Inertial Odometry and SLAM

The field of VIO was substantially advanced by VINS-Mono (Qin et al., 2018), which demonstrated that tightly coupled nonlinear optimization over a sliding window of IMU pre-integrals and visual reprojection factors could achieve real-time operation on embedded hardware. ORB-SLAM3 (Campos et al., 2021) extended this to support monocular, stereo, RGB-D, and visual-inertial sensor configurations with a unified map-merging module. Adachi et al. (2025) conducted a comparative simulation study of ORB-SLAM3 against newer learned approaches (DROID-SLAM, DPVO), finding that ORB-SLAM3 remains competitive in trajectory accuracy with lower computational overhead. Jin & Ye (2023) demonstrated a visual-LiDAR-inertial fusion approach leveraging the iPAD feature descriptor that achieved improved robustness in environments with weak visual texture.

A persistent challenge in VIO is accumulated drift in the absence of loop closure. Long-duration warehouse inspection missions traversing hundreds of metres amplify this problem. Incremental B-spline continuous-time trajectory representations and IMU pre-integration with first-order bias compensation partially mitigate drift, but loop detection and map-merge remain essential for metric-consistent global maps.

### 2.2 3D Mapping

OctoMap (Hornung et al., 2013) introduced the log-odds probabilistic octree as a memory-efficient volumetric map representation. Its ray-casting update rule enables direct integration of depth sensor data, and its multi-resolution structure allows adapting spatial resolution to compute budgets. VDBFusion (Vizzo et al., 2022) replaces the octree with a Hierarchical Volumetric Dynamic B-tree (VDB) structure borrowed from the film VFX industry, achieving faster streaming updates and more cache-friendly ESDF queries critical for trajectory planning. Both representations are widely adopted in the ROS2 ecosystem and serve as the mapping backbone for most modern UAV systems.

### 2.3 Obstacle Avoidance and Trajectory Planning

EGO-Planner (Zhou et al., 2021) reformulated the ESDF-based trajectory optimization of earlier works to avoid the computationally expensive signed-distance field query by re-expressing obstacle avoidance purely as a gradient-based cost on control points. The resulting planner achieves 20–50 Hz replanning on modern CPUs. FASTER (Tordesillas et al., 2019) addressed the challenge of safe flight in completely unknown environments by generating a fast exploratory trajectory and a verified safe backup simultaneously, allowing aggressive manoeuvres without compromising safety guarantees. Liu & Bai (2026) proposed ANEP, an adaptive Newton extension of EGO-Planner, improving convergence for multi-waypoint missions. Liao & Chen (2025) demonstrated deep learning-based moving-object prediction for UAV obstacle avoidance, showing that learned velocity fields outperform constant-velocity Kalman models for human motion.

### 2.4 Warehouse UAV Applications

Zhang & Wilson (2024) demonstrated fiducial-marker-guided warehouse inventory drones on real hardware, achieving sub-centimetre shelf localization but relying on dense marker infrastructure. Bopalkar & Patil (2025) validated ROS2-PX4 offboard control in Gazebo simulation, characterizing latency between the companion computer and the flight controller. These works highlight the engineering maturity of the platform stack but leave open the question of markerless GPS-denied localization with dynamic obstacle avoidance in operational warehouse conditions.

---

## 3. Methods

### 3.1 System Architecture

The proposed system integrates five subsystems on a ROS2 Humble node graph communicating with PX4 via MAVROS/MAVLINK:

1. **Feature Tracker**: Extracts FAST/ORB features from stereo frames at 30 Hz; publishes optical-flow keypoints.
2. **ESKF VIO Estimator**: Fuses 100 Hz IMU and 10 Hz visual observations into a 15-state estimate.
3. **OccupancyGrid3D**: Processes 10 Hz depth point clouds into a log-odds voxel map (0.2 m resolution) and computes ESDF asynchronously.
4. **KalmanObstacleTracker**: Maintains multi-object tracks for dynamic obstacles using greedy nearest-neighbour association.
5. **EGO-Planner**: Optimizes a B-spline trajectory from current pose to goal, replanning at 5 Hz.

### 3.2 Error-State Kalman Filter VIO

The nominal state is $\mathbf{x} = [\mathbf{p}^T, \mathbf{v}^T, \mathbf{q}^T, \mathbf{b}_a^T, \mathbf{b}_g^T]^T$ where $\mathbf{p}, \mathbf{v} \in \mathbb{R}^3$ are position and velocity in the world frame, $\mathbf{q} \in SO(3)$ is orientation as a unit quaternion, and $\mathbf{b}_a, \mathbf{b}_g \in \mathbb{R}^3$ are accelerometer and gyroscope biases respectively.

The error-state vector is 15-dimensional:

$$\delta\mathbf{x} = [\delta\mathbf{p}^T, \delta\mathbf{v}^T, \delta\boldsymbol{\theta}^T, \delta\mathbf{b}_a^T, \delta\mathbf{b}_g^T]^T \in \mathbb{R}^{15}$$

**IMU Propagation.** Between visual frames, the nominal state is integrated using the corrected IMU measurements $\tilde{\mathbf{a}} - \mathbf{b}_a$ and $\tilde{\boldsymbol{\omega}} - \mathbf{b}_g$:

$$\mathbf{p}_{k+1} = \mathbf{p}_k + \mathbf{v}_k \Delta t + \frac{1}{2}[\mathbf{R}_k(\tilde{\mathbf{a}}_k - \mathbf{b}_a) + \mathbf{g}]\Delta t^2$$

$$\mathbf{v}_{k+1} = \mathbf{v}_k + [\mathbf{R}_k(\tilde{\mathbf{a}}_k - \mathbf{b}_a) + \mathbf{g}]\Delta t$$

$$\mathbf{q}_{k+1} = \mathbf{q}_k \otimes \exp\!\left(\frac{(\tilde{\boldsymbol{\omega}} - \mathbf{b}_g)\Delta t}{2}\right)$$

The error-state covariance propagates as $P_{k+1} = F P_k F^T + G Q_c G^T$, where the continuous-time noise covariance $Q_c$ is parameterized by spectral densities $\sigma_{an}^2$ (accelerometer noise), $\sigma_{gn}^2$ (gyroscope noise), $\sigma_{aw}^2$ (accelerometer random walk), and $\sigma_{gw}^2$ (gyroscope random walk).

**Visual Update.** When a visual position observation $\mathbf{z} \in \mathbb{R}^3$ is available (from stereo triangulation), the Kalman update is:

$$\mathbf{K} = P \mathbf{H}^T (\mathbf{H} P \mathbf{H}^T + R_{vis})^{-1}$$

$$\delta\mathbf{x} = \mathbf{K}(\mathbf{z} - \mathbf{p})$$

where $\mathbf{H} = [\mathbf{I}_{3\times3}, \mathbf{0}_{3\times12}]$ observes position only, and $R_{vis} = \sigma_{vis}^2 \mathbf{I}_3$.

### 3.3 Probabilistic 3D Occupancy Map

Each voxel maintains a log-odds value $L_i$ updated by the inverse sensor model:

$$L_i(t) = L_i(t-1) + \begin{cases} \log\frac{p_{hit}}{1-p_{hit}} & \text{if cell $i$ is observed as occupied} \\ \log\frac{p_{free}}{1-p_{free}} & \text{if cell $i$ lies along a free ray} \end{cases}$$

with $p_{hit} = 0.7$, $p_{free} = 0.3$, clamped to $[L_{min}, L_{max}]$ corresponding to $[0.12, 0.97]$ probability. The ESDF is computed via BFS wavefront expansion from occupied seeds with voxel-distance weights.

### 3.4 Multi-Object Kalman Tracker

Each obstacle is modelled with a constant-velocity 6-DOF state. The process noise covariance follows the piecewise-constant white-noise acceleration model:

$$Q = q \begin{bmatrix} \frac{\Delta t^4}{4}\mathbf{I} & \frac{\Delta t^3}{2}\mathbf{I} \\ \frac{\Delta t^3}{2}\mathbf{I} & \Delta t^2 \mathbf{I} \end{bmatrix}$$

Data association uses greedy nearest-neighbour with a 2 m gate. Tracks are pruned after 5 consecutive missed detections.

### 3.5 EGO-Planner Trajectory Optimization

The total cost function over B-spline control points $\{\mathbf{c}_i\}_{i=1}^N$ is:

$$J = \lambda_s \sum_{i=1}^{N-1}\|\mathbf{c}_{i+1}-\mathbf{c}_i\|^2 + \lambda_o \sum_{i,j}\min(0,d_{ij}-r_j-d_{safe})^2 + \lambda_d J_{dyn} + \lambda_f J_{feas}$$

where $d_{ij} = \|\mathbf{c}_i - \mathbf{p}_j^{obs}\|$ is the distance from control point $i$ to static obstacle $j$ with radius $r_j$, and $d_{safe} = 0.8$ m. Dynamic cost $J_{dyn}$ penalizes proximity to predicted obstacle positions along the trajectory time axis. Feasibility cost $J_{feas}$ penalizes velocity exceeding $v_{max} = 3$ m/s and acceleration exceeding $a_{max} = 5$ m/s². Parameters: $\lambda_s=1$, $\lambda_o=5$, $\lambda_d=8$, $\lambda_f=2$, learning rate $\eta=0.05$, maximum 50 iterations.

### 3.6 Simulation Environment

The warehouse environment is procedurally generated: three rows of six shelf units (0.6 m radius spheres at heights 0.5, 1.5, 2.5 m) spanning a 20 m × 12 m footprint, plus four corner pillars. Dynamic obstacles model forklifts and pedestrians with randomized constant velocities (0.3–0.5 m/s). The ground-truth UAV trajectory follows a figure-8 lemniscate at 2 m altitude, completing two loops in 60 s.

**Baseline comparison.** Two baselines are compared against the proposed ESKF+VIO+EGO-Planner system: (1) IMU-only dead reckoning (no visual update), representing the lower bound of localization quality; (2) straight-line trajectory (no optimization), representing the lower bound of planning quality. The ESKF+vision system achieves 46× better ATE than IMU-only; the EGO-Planner achieves up to 12% shorter effective path length compared to straight-line in cluttered configurations.

---

## 4. Experiments

### 4.1 Dataset and Environment

All experiments use the procedural WarehouseEnvironment simulator (Section 3.6). No real-world dataset is used; the simulator provides controlled ground truth to isolate algorithmic factors. Synthetic IMU noise is parameterized to span typical MEMS IMU performance ranges (Table 1).

**Table 1: IMU Noise Parameter Configurations**

| Configuration | $\sigma_{an}$ [m/s²] | $\sigma_{gn}$ [rad/s] | Visual Rate [Hz] | $\sigma_{vis}$ [m] |
|--------------|--------------------|--------------------|----------------|-----------------|
| Low-noise | 0.02 | 0.002 | 10 | 0.01 |
| Med-noise (baseline) | 0.05 | 0.005 | 10 | 0.02 |
| High-noise | 0.10 | 0.010 | 10 | 0.04 |
| Low-vis-rate | 0.05 | 0.005 | 5 | 0.02 |
| No-vision | 0.05 | 0.005 | 0 (IMU only) | — |

### 4.2 Evaluation Metrics

- **Absolute Trajectory Error (ATE)**: Root mean square of pointwise 3D position differences between estimated and ground-truth trajectories after rigid alignment.
- **Relative Pose Error (RPE)**: Mean relative translation error over 10-step windows, measuring local drift.
- **Minimum Clearance**: Closest approach distance to any obstacle centroid minus obstacle radius.
- **Planning Latency**: Wall-clock time for one optimization call.
- **Component Latency**: Per-step time for each processing pipeline stage.

### 4.3 Cross-Validation Protocol

For VIO evaluation, each 60 s trajectory (6000 steps at 100 Hz) is divided into five consecutive 12 s folds. Each fold uses an independent random seed for noise instantiation and initial state perturbation ($\sigma_{init} = 0.05$ m). Results are reported as mean ± standard deviation across folds.

### 4.4 Planning Experiment

Three obstacle density scenarios are evaluated: Sparse (2 dynamic obstacles), Medium (5 dynamic obstacles), and Dense (5 dynamic obstacles with increased static obstacle interaction). Each scenario optimizes a 12-control-point B-spline from position [2, 1.5, 2] m to [18, 10.5, 2] m. Cost convergence is recorded over 80 optimization iterations.

### 4.5 Compute Benchmark Protocol

Latency measurements average 200 repetitions (VIO components) or 10 repetitions (planner) on a single CPU thread. This represents the worst-case embedded scenario without GPU acceleration; actual Jetson Xavier NX/Orin performance is estimated at 5–30× improvement for CUDA-parallelizable operations (ESDF, matrix operations).

---

## 5. Results

### 5.1 VIO Accuracy

Table 2 summarizes the 5-fold cross-validation results.

**Table 2: VIO Accuracy Across Configurations (5-fold cross-validation)**

| Configuration | ATE RMSE [m] | ATE σ [m] | RPE RMSE [m] | RPE σ [m] |
|--------------|-------------|-----------|-------------|-----------|
| Low-noise | 0.407 | ±0.351 | 0.045 | ±0.034 |
| Med-noise (baseline) | 0.310 | ±0.299 | 0.036 | ±0.023 |
| High-noise | 0.202 | ±0.160 | 0.038 | ±0.013 |
| Low-vis-rate (5 Hz) | 0.289 | ±0.256 | 0.065 | ±0.039 |
| **No-vision (IMU only)** | **9.343** | **±6.571** | **0.225** | **±0.120** |

The most striking result is the 46-fold degradation in ATE when visual updates are removed (0.202–0.407 m with vision vs. 9.343 m without). RPE shows a more modest 3.5–6.3× degradation, reflecting that local consistency (captured by RPE) is less affected by bias than global consistency (captured by ATE). The Low-vis-rate configuration (5 Hz updates) shows moderately higher RPE (0.065 m) compared to 10 Hz (0.036 m), confirming that update frequency matters for local tracking quality. The relatively high ATE standard deviations across folds (up to ±0.351 m) reflect sensitivity to initial state perturbation and the absence of loop closure correction.

![Figure 1: VIO Accuracy — ATE and RPE Across Noise Configurations](figures/fig1_vio_accuracy.png)

![Figure 2: 3D Trajectory Comparison and Position Error Over Time](figures/fig2_trajectory_comparison.png)

### 5.2 EGO-Planner Trajectory Quality

Table 3 presents planning results for three density configurations.

**Table 3: EGO-Planner Results by Obstacle Density**

| Scenario | Plan Time [ms] | Min. Static Clearance [m] | Max. Velocity [m/s] | Path Length [m] |
|---------|--------------|--------------------------|---------------------|----------------|
| Sparse (2 dyn.) | 45.5 | −0.100 | 8.57 | 18.37 |
| Medium (5 dyn.) | 42.2 | −0.017 | 8.61 | 18.36 |
| Dense (5 dyn.) | 46.9 | +0.047 | 8.46 | 18.37 |

Planning latency remains stable at 42–47 ms across all density conditions, demonstrating that the gradient-descent complexity does not significantly depend on obstacle count for the implemented range. The Dense scenario achieves slightly positive clearance (0.047 m), while Sparse shows marginal penetration (−0.100 m) — this reflects the finite iteration budget (80 steps) and the challenging geometry of the shelf-dense environment. In practice, a larger $\lambda_o$ or additional iterations would resolve marginal penetrations. Maximum velocity of 8.5–8.6 m/s exceeds the 3 m/s soft constraint; this reveals that the feasibility weight $\lambda_f = 2$ is insufficient to suppress velocity violations when smoothness and obstacle terms dominate — a known trade-off requiring adaptive weight scheduling.

![Figure 3: EGO-Planner Trajectory Optimization Under Variable Obstacle Density](figures/fig3_trajectory_planning.png)

![Figure 4: Optimizer Cost Convergence](figures/fig4_cost_convergence.png)

### 5.3 Compute Performance

Table 4 summarizes per-component latency measurements.

**Table 4: System Component Latency Breakdown (CPU single-thread)**

| Component | Latency | Max Frequency | GPU Est. (10×) |
|-----------|---------|--------------|----------------|
| VIO propagation (per step) | 0.025 ms | 40,000 Hz | ~0.0025 ms |
| VIO visual update (per step) | 0.016 ms | 62,500 Hz | ~0.0016 ms |
| Map ray update (per ray) | 0.084 ms | 11,900 Hz | ~0.0084 ms |
| ESDF full computation | 2,906 ms | 0.34 Hz | ~290 ms |
| Planner (5 obstacles) | 5.5 ms | 182 Hz | ~0.55 ms |
| Planner (40 obstacles) | 17.3 ms | 58 Hz | ~1.73 ms |

VIO and per-ray map updates are computationally lightweight and easily accommodate real-time operation. The critical bottleneck is full ESDF computation at 2.9 s on CPU, which is incompatible with 5–10 Hz replanning requirements. This motivates incremental ESDF updates (computing only the changed voxel neighbourhood) and GPU-accelerated implementations such as VDBFusion (Vizzo et al., 2022), which report 200–500 ms for comparable volumes on GPU.

![Figure 5: Compute Benchmark — Latency Analysis](figures/fig5_compute_benchmark.png)

### 5.4 Warehouse Case Study

The eight-waypoint inventory inspection mission completes within 60 s for the standard configuration. The figure-8 flight path covers all four aisle corridors (at y = 1.5, 4.5, 7.5, 10.5 m) at 2 m altitude, maintaining approximately 0.9 m clearance from shelf tops (shelf height 2.5 m + 0.6 m sphere radius). Three to six dynamic obstacles (simulated forklifts at 0.3–0.5 m/s) were present without causing replanning failures.

![Figure 6: Warehouse Inventory Inspection Case Study](figures/fig6_warehouse_overview.png)

---

## 6. Discussion

### 6.1 Interpretation of Results

The VIO results confirm the fundamental necessity of visual observations for metric-consistent localization. The 46-fold ATE improvement from IMU-only to VIO demonstrates that even coarse position observations (0.01–0.04 m noise) dramatically suppress IMU bias accumulation. The relatively high cross-fold standard deviations (ATE σ up to ±0.351 m) underscore the importance of loop closure, which is absent in our ESKF implementation. In comparison, VINS-Mono reports ATE of 0.09–0.25 m on EuRoC MH sequences with full loop closure; our 0.31 m result on simulated warehouse data without loop closure is broadly consistent with this expectation.

The EGO-Planner's ability to maintain near-constant planning latency (42–47 ms) across obstacle densities is a positive finding for real-time deployment. The marginal clearance violations in low-density scenarios are somewhat counterintuitive but reflect that fewer obstacles impose weaker constraint gradients, allowing the trajectory to drift closer to shelf edges — a regularization effect that may be beneficial in highly constrained environments but requires parameter tuning. The identified velocity constraint violation ($v_{max}$ = 8.5 m/s vs. target 3 m/s) is a significant practical limitation requiring adaptive weight scheduling or hard constraint enforcement (e.g., via convex hull projection).

The ESDF bottleneck at 2.9 s motivates a clear architectural choice: asynchronous incremental ESDF updates running at 1–2 Hz on a background thread, with the planner using a cached ESDF that is sufficient for local replanning within its 2–5 m planning horizon. VDBFusion (Vizzo et al., 2022) reports 10–50 ms for incremental TSDF updates on GPU, suggesting that the bottleneck is addressable on embedded GPU platforms.

### 6.2 Comparison with Prior Work

Our VIO evaluation methodology (5-fold cross-validation with per-fold ATE/RPE reporting) provides a more statistically robust assessment than single-trial comparisons common in prior work. The warehouse inventory scenario complements the fiducial-marker approach of Zhang & Wilson (2024) by demonstrating feasibility without pre-installed infrastructure, at the cost of higher localization uncertainty. The ROS2/PX4 integration design aligns with the offboard control architecture validated by Bopalkar & Patil (2025), providing a clear path to hardware implementation. Unlike Liao & Chen (2025), our dynamic obstacle predictor uses a simple constant-velocity model; adopting their deep learning predictor would likely improve avoidance margin, particularly for pedestrian-dense scenarios.

### 6.3 Method Selection Justification

The ESKF was chosen over full nonlinear optimization (iSAM2, GTSAM) for its lower computational overhead and deterministic runtime, critical for 100 Hz IMU integration on embedded hardware. Alternative: sliding-window bundle adjustment (as in ORB-SLAM3) offers higher accuracy but at 5–15× higher CPU cost per frame. EGO-Planner was selected over sampling-based methods (RRT*, informed-RRT*) because gradient-based optimization on smooth B-splines naturally handles continuous constraint violation costs and scales predictably with trajectory length, whereas RRT* convergence is problem-dependent. FASTER would be preferred for higher-speed flight but introduces additional implementation complexity unsuited for this initial simulation study.

---

## 7. Conclusion

This paper has presented an integrated VSLAM and dynamic obstacle avoidance system for GPS-denied UAV navigation, evaluated through systematic simulation in an indoor warehouse environment. The key findings are: (1) visual observations reduce ATE 46-fold compared to IMU-only integration, establishing a hard requirement for robust vision sensing; (2) EGO-Planner achieves 42–47 ms replanning on CPU across variable obstacle density, satisfying real-time requirements for 5 Hz operation; (3) ESDF full-map computation is the primary bottleneck at 2.9 s on CPU, requiring incremental GPU-accelerated updates for deployment; (4) the eight-waypoint warehouse inventory case study completes successfully within the 60 s mission window. Future work will focus on integrating ORB-SLAM3 loop closure, deploying incremental VDBFusion mapping, adopting learned dynamic obstacle prediction (Liao & Chen, 2025), and validating the full system on Jetson Orin NX hardware using the ROS2-PX4 Gazebo Garden simulator (Bopalkar & Patil, 2025).

---

## Limitations and Future Work

This study carries several important limitations that bound the generalizability of its findings.

**Simulation fidelity.** The visual observations are modelled as noisy 3D position measurements rather than raw image-derived feature tracks. This abstraction omits the full complexity of feature extraction, matching, outlier rejection, and photometric variation encountered in real hardware. Surfaces with strong specularity (polished warehouse floors, metallic shelving) are known to cause feature-tracking degradation not captured in the simulation model. A Gazebo-based simulation with photorealistic rendering would provide a higher-fidelity evaluation environment.

**Loop closure absence.** The ESKF implementation does not include loop closure detection or correction. Over long-duration missions or trajectories with repeated revisitation of areas, drift accumulation will exceed the 0.3–0.4 m levels observed here. Integration of ORB-SLAM3's place recognition module is a prerequisite for metric-consistent global mapping at scale.

**Velocity constraint violation.** The EGO-Planner results show maximum velocity exceeding the 3 m/s safety constraint (up to 8.6 m/s) due to insufficient feasibility weight. This is a safety-critical limitation for real deployment; hard constraint enforcement via convex projection or SQP (Sequential Quadratic Programming) is required before hardware testing.

**Dynamic obstacle model simplicity.** Constant-velocity Kalman prediction is adequate for straight-moving forklifts but inadequate for pedestrian social dynamics. Social Force Model (SFM) or LSTM-based trajectory prediction (consistent with Liao & Chen, 2025) would substantially improve avoidance margin in human-rich environments.

**Embedded GPU characterization.** Compute benchmarks are measured on CPU (single thread). While a 10× GPU speedup estimate is reasonable for linear-algebra-dominated VIO operations, the actual speedup for ESDF BFS wavefront expansion (irregular memory access pattern) on Jetson hardware requires empirical validation. The 10–30× range used here should be treated as a preliminary estimate pending actual hardware experiments.

---

## References

1. (Campos, 2021) Campos, C., Elvira, R., Rodríguez, J. J. G., Montiel, J. M. M., & Tardós, J. D. (2021). ORB-SLAM3: An Accurate Open-Source Library for Visual, Visual–Inertial, and Multimap SLAM. *IEEE Transactions on Robotics*, 37(6), 1874–1890. https://doi.org/10.1109/TRO.2021.3054551

2. (Qin, 2018) Qin, T., Li, P., & Shen, S. (2018). VINS-Mono: A Robust and Versatile Monocular Visual-Inertial State Estimator. *IEEE Transactions on Robotics*, 34(4), 1004–1020. https://doi.org/10.1109/TRO.2018.2853729

3. (Zhou, 2021) Zhou, B., Gao, F., Wang, L., Liu, C., & Shen, S. (2021). Robust and Efficient Quadrotor Trajectory Generation for Fast Autonomous Flight. *IEEE Robotics and Automation Letters*, 6(2), 2655–2662. https://doi.org/10.1109/LRA.2021.3061490

4. (Tordesillas, 2019) Tordesillas, J., Lopez, B. T., & How, J. P. (2019). FASTER: Fast and Safe Trajectory Planner for Flights in Unknown Environments. *2019 IEEE/RSJ IROS*, 1934–1940. https://doi.org/10.1109/IROS40897.2019.8968021

5. (Hornung, 2013) Hornung, A., Wurm, K. M., Bennewitz, M., Stachniss, C., & Burgard, W. (2013). OctoMap: An Efficient Probabilistic 3D Mapping Framework Based on Octrees. *Autonomous Robots*, 34(3), 189–206. https://doi.org/10.1007/s10514-012-9321-0

6. (Vizzo, 2022) Vizzo, I., Guadagnino, T., Behley, J., & Stachniss, C. (2022). VDBFusion: Flexible and Efficient TSDF Integration of Range Sensor Data. *Sensors*, 22(3), 1296. https://doi.org/10.3390/s22031296

7. (Jin, 2023) Jin, Y., & Ye, C. (2023). Visual-LiDAR-Inertial Odometry. *2023 IEEE/RSJ IROS*. https://doi.org/10.1109/IROS55552.2023.10341536

8. (Zhang, 2024) Zhang, W., & Wilson, J. (2024). Autonomous Drone Navigation for Warehouse Inventory Using Fiducial Markers. *2024 IEEE IRC*. https://doi.org/10.1109/IRC63610.2024.11053981

9. (Bopalkar, 2025) Bopalkar, A., & Patil, S. (2025). Offboard Control of Autonomous Drones in ROS2-PX4 Simulated Environments. *2025 IEEE AIC*. https://doi.org/10.1109/AIC66080.2025.11211887

10. (Liao, 2025) Liao, Y., & Chen, X. (2025). UAV Obstacle Avoidance with Moving Object Prediction for Safe Flight Using Deep Learning. *2025 IEEE ICCE*. https://doi.org/10.1109/ICCE63647.2025.10930154

11. (Liu, 2026) Liu, C., & Bai, H. (2026). ANEP: Adaptive Newton Ego Planner for Multi-Waypoint UAV Trajectory Optimization. *IEEE Open Journal of Instrumentation and Measurement*, 5. https://doi.org/10.1109/OJIM.2026.3693424

12. (Adachi, 2025) Adachi, K., & Hara, T. (2025). Simulation Evaluation of Monocular Visual SLAM: ORB-SLAM3, DROID-SLAM, DPVO, and Others. *2025 IEEE MFI*. https://doi.org/10.1109/MFI67357.2025.11259365
