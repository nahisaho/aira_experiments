# Local Path Planning and Mission Planning Design

## Scope and assumptions
- Vehicle: quadrotor-class autonomous drone using ROS 2 + PX4 with onboard perception in a GPS-denied warehouse.
- Local planner consumes ESDF, free-space topology, robot state estimate, and dynamic obstacle predictions.
- Control stack: local planner generates dynamically feasible trajectories; PX4 executes via offboard position/trajectory interfaces.

## 1. Planning stack overview

### 1.1 Hierarchical planning architecture
1. **Global mission layer**
   - Warehouse shelf coverage route generation
   - TSP-based shelf ordering
   - battery-aware segmentation and return-to-home (RTH)
2. **Local motion layer**
   - EGO-Planner v2 B-spline trajectory optimization against ESDF and dynamic obstacles
   - replanning at 10 Hz over a 5-8 s horizon
3. **Control layer**
   - SE(3) geometric tracking controller with feed-forward and feedback

### 1.2 Data interfaces
| Topic | Source | Consumer | Purpose |
|---|---|---|---|
| `/mapping/esdf` | `vdbfusion_node` | `ego_planner_node` | Signed distance and gradient queries |
| `/perception/predicted_trajectories` | predictor | `ego_planner_node` | Dynamic obstacle future occupancy |
| `/mission/global_waypoints` | mission planner | local planner | Current sub-goal sequence |
| `/planner/local_bspline` | local planner | tracker/controller bridge | Executable local trajectory |
| `/fmu/in/trajectory_setpoint` | ROS 2 bridge | PX4 | Offboard setpoints |

## 2. Local path planning with EGO-Planner v2

### 2.1 Why EGO-Planner v2
EGO-Planner v2 is well suited to onboard quadrotor planning because it:
- directly optimizes B-spline control points,
- uses ESDF gradients for fast collision costs,
- supports frequent replanning without graph search at every cycle,
- is lighter than sampling-heavy MPPI and more natural for aerial systems than TEB.

### 2.2 Trajectory representation
Use a uniform B-spline of degree `p=3` or `p=4` with control points `\mathbf{Q}_i` and knot interval `\Delta t`:

```math
\mathbf{r}(t) = \sum_{i=0}^{n} N_{i,p}(t) \mathbf{Q}_i
```

Velocity and acceleration follow directly from spline derivatives:

```math
\dot{\mathbf{r}}(t) = \sum_i \dot{N}_{i,p}(t) \mathbf{Q}_i, \quad
\ddot{\mathbf{r}}(t) = \sum_i \ddot{N}_{i,p}(t) \mathbf{Q}_i
```

Advantages:
- local support enables fast incremental updates,
- smoothness constraints are naturally encoded,
- control points provide compact optimization variables.

### 2.3 ESDF from VDBFusion
VDBFusion incrementally fuses depth/point cloud measurements into a sparse voxel structure. The planner queries:
- signed distance `d(\mathbf{x})`,
- spatial gradient `\nabla d(\mathbf{x})`,
- free/unknown/occupied status.

Recommended map parameters:
- voxel size: `0.1-0.2 m`,
- truncation distance: `0.4-0.6 m`,
- local map rolling window: `20 x 20 x 6 m` around the drone.

Use TSDF-to-ESDF conversion or direct ESDF maintenance in the local planning volume. Unknown space policy should be configurable:
- conservative: unknown treated as occupied during mission-critical operation,
- exploratory: unknown penalized but not fully blocked in controlled validation.

### 2.4 Optimization objective
Optimize control points by minimizing

```math
J = w_s J_{smooth} + w_c J_{coll} + w_d J_{dyn} + w_o J_{dynobs} + w_g J_{goal}
```

where:
1. **Smoothness**

```math
J_{smooth} = \int_0^T \left( \|\dddot{\mathbf{r}}(t)\|^2 + \lambda_{snap}\|\mathbf{r}^{(4)}(t)\|^2 \right) dt
```

2. **Static collision avoidance** using ESDF

```math
J_{coll} = \int_0^T \phi(d(\mathbf{r}(t))) dt
```

with hinge penalty

```math
\phi(d) = \begin{cases}
(d_{safe} - d)^2, & d < d_{safe} \\
0, & d \ge d_{safe}
\end{cases}
```

and gradient

```math
\frac{\partial J_{coll}}{\partial \mathbf{r}} \propto -2(d_{safe}-d)\nabla d(\mathbf{r})
```

3. **Dynamic feasibility**

```math
J_{dyn} = \int_0^T \psi_v(\|\dot{\mathbf{r}}(t)\| - v_{max}) + \psi_a(\|\ddot{\mathbf{r}}(t)\| - a_{max}) dt
```

with soft penalties on velocity and acceleration limit violations.

4. **Dynamic obstacle avoidance**

For obstacle `j` with predicted mean `\mu_j(t)` and covariance `\Sigma_j(t)`, define inflated relative covariance `\Sigma_{rel,j}(t)` and penalize trajectories violating a chance-constrained margin:

```math
J_{dynobs} = \sum_j \int_0^T \phi_o\left( \rho - \sqrt{(\mathbf{r}(t)-\mu_j(t))^T \Sigma_{rel,j}^{-1}(t) (\mathbf{r}(t)-\mu_j(t))} \right) dt
```

where `\rho` is the required Mahalanobis clearance threshold.

5. **Goal progress**

```math
J_{goal} = \|\mathbf{r}(T) - \mathbf{r}_{goal}\|^2 + \lambda_\psi \|\psi(T)-\psi_{goal}\|^2
```

### 2.5 Dynamic obstacle incorporation
At each replan cycle:
1. Sample predicted obstacle trajectories over the 5-8 s horizon.
2. Convert each sample to an inflated occupancy tube using obstacle dimensions, drone radius, and uncertainty margin.
3. Add analytic penalty terms to `J_dynobs`.
4. For imminent threats, also modify the terminal goal to bias lateral or vertical escape corridors.

Recommended inflation:

```math
r_{infl}(t) = r_{geom} + k_\sigma \sqrt{\lambda_{max}(\Sigma_j(t))}
```

with `k_\sigma = 2` for approximately 95% confidence.

### 2.6 Replanning triggers
Trigger replanning if any of the following occurs:
- fixed timer at 10 Hz,
- local trajectory collision with updated ESDF,
- new dynamic obstacle enters avoidance zone,
- deviation from nominal trajectory exceeds threshold,
- goal update from mission planner,
- velocity/state estimator reset,
- map change introduces new occupied space near the active path.

Recommended thresholds:
- state deviation: position error `> 0.35 m` or yaw error `> 10 deg`,
- predicted collision probability along active path `> 0.2`,
- ESDF safety margin violation `d < 0.5 m`.

### 2.7 Replanning frequency and horizon
- Replanning frequency: **10 Hz**
- Planning horizon: **5-8 s**
- Execution commitment window: `0.5-1.0 s`

Use a receding-horizon policy: optimize over the full horizon while only committing the initial trajectory segment. This balances foresight and responsiveness.

### 2.8 Initialization and fallback behavior
- Warm-start with shifted previous B-spline control points.
- If optimization fails within cycle budget, execute a safety fallback:
  1. use last valid safe prefix,
  2. reduce speed,
  3. if no safe prefix exists, transition to hover/brake.

Maximum local optimization budget should be `<= 40-50 ms` to remain compatible with 10 Hz replanning.

## 3. Planner comparison

| Planner | Strengths | Weaknesses | Fit for this system |
|---|---|---|---|
| **EGO-Planner v2** | Fast gradient-based B-spline optimization, strong fit for quadrotors, good onboard performance | Needs reliable ESDF and good warm starts | **Recommended default** |
| **FASTER** | Explicit free-space/unknown-space reasoning, robust in exploration | More complex corridor generation and integration overhead | Good backup if unknown-space exploration dominates |
| **TEB** | Mature local planner, strong for differential-drive robots | Less natural for full 3D aerial dynamics; parameter tuning can be difficult | Not preferred for warehouse drones |
| **MPPI** | Handles complex dynamics and nonconvex objectives | High compute demand, many rollouts, latency-sensitive on embedded systems | Good research option, not first deployment choice |

Decision summary:
- Use **EGO-Planner v2** as the production local planner.
- Consider **FASTER** if the mission expands into unknown-space exploration.
- Reserve **MPPI** for future work when higher compute or learned cost shaping becomes available.

## 4. Global mission planning

### 4.1 Coverage path planning for warehouse inventory
Model the warehouse as aisles and shelf faces to be inspected. Each shelf face is a coverage target with a required observation pose set.

Workflow:
1. Build a graph of aisles, intersections, and shelf observation viewpoints.
2. Generate candidate viewpoints satisfying camera FOV, stand-off distance, and occlusion constraints.
3. Solve visit ordering with TSP or clustered TSP.
4. For each aisle segment, produce nominal corridor waypoints for the local planner.

A shelf viewpoint `v_i` should satisfy:
- distance to shelf `d_i` in sensor sweet spot,
- normal alignment within camera incidence constraint,
- visibility score above threshold.

### 4.2 TSP-based optimal shelf visit ordering
Define graph `G=(V,E)` where `V` are shelf targets plus depot/home node and edge weights are expected flight costs:

```math
c_{ij} = \alpha L_{ij} + \beta E_{ij} + \gamma R_{ij}
```

where:
- `L_{ij}` = travel distance/time,
- `E_{ij}` = estimated energy cost,
- `R_{ij}` = risk penalty from congestion/no-fly proximity.

Solve a TSP or VRP variant:
- exact solver for small missions,
- Lin-Kernighan or OR-Tools heuristic for larger shelf sets.

### 4.3 Battery-aware mission segmentation
Let battery state of charge be `SOC`. Partition the mission into segments so that for each segment `s`:

```math
E_{seg,s} + E_{reserve} + E_{RTH,s} \le E_{avail}(SOC)
```

Use a conservative energy model:

```math
E = \int_0^T P(t) dt \approx \int_0^T (P_{hover} + k_v \|v\| + k_a \|a\| + k_w \|w\|) dt
```

where `w` reflects estimated wind disturbance or ventilation load.

Recommended policy:
- mission launch allowed only if predicted full segment plus reserve is feasible,
- reserve margin: `20-25%` SOC,
- preemptive return if remaining energy falls below `E_RTH + E_reserve`.

### 4.4 Return-to-home planning with safety margin
RTH should select the safest feasible route to a landing or docking station.
- maintain one or more home nodes,
- choose route minimizing risk-weighted travel cost,
- enforce extra clearance around shelves and dynamic congestion areas.

RTH trigger examples:
- low battery,
- perception degradation,
- mission timeout,
- repeated local planning failure.

### 4.5 Geofencing for warehouse no-fly zones
Represent forbidden regions as 3D polyhedra or voxel masks.
- permanent no-fly zones: charging infrastructure, humans-only corridors, sprinkler clearances,
- temporary no-fly zones: forklift work cells, blocked aisles, maintenance areas.

Enforce geofence in both global and local layers:
- global mission planner excludes these regions from route generation,
- local planner adds hard collision penalties or hard constraints,
- PX4 geofence/flight termination layer acts as final backstop.

## 5. Trajectory tracking controller

### 5.1 Geometric controller on SE(3)
For desired trajectory `(\mathbf{x}_d, \dot{\mathbf{x}}_d, \ddot{\mathbf{x}}_d, R_d, \Omega_d)`, define position and velocity errors:

```math
\mathbf{e}_x = \mathbf{x} - \mathbf{x}_d, \quad \mathbf{e}_v = \mathbf{v} - \dot{\mathbf{x}}_d
```

Attitude error:

```math
e_R = \frac{1}{2}(R_d^T R - R^T R_d)^\vee
```

Angular velocity error:

```math
e_\Omega = \Omega - R^T R_d \Omega_d
```

Thrust command:

```math
f = \left( -k_x e_x - k_v e_v - mg e_3 + m\ddot{x}_d + \Delta_{ff} \right) \cdot Re_3
```

Moment command:

```math
M = -k_R e_R - k_\Omega e_\Omega + \Omega \times J\Omega - J(\hat{\Omega}R^TR_d\Omega_d - R^TR_d\dot{\Omega}_d)
```

where `\Delta_ff` is a feed-forward term from planned acceleration/jerk and optional disturbance estimate.

### 5.2 Feed-forward + feedback architecture
- **Feed-forward:** desired acceleration, jerk, yaw rate from the B-spline trajectory.
- **Feedback:** position/velocity/attitude stabilization from onboard state estimates.
- **Interface to PX4:**
  - either use PX4 native position/trajectory setpoints and let PX4 low-level loops close the remaining gap,
  - or embed the SE(3) controller in ROS 2 and command body rates/thrust if certification and timing allow.

For production robustness, keep inner-rate stabilization in PX4 and use the ROS 2 controller as a high-rate outer loop or trajectory setpoint generator.

### 5.3 Wind disturbance rejection
Indoor warehouses may have ventilation jets and downwash recirculation. Add disturbance rejection via:
- velocity-integral or disturbance observer term,
- drag model compensation,
- adaptive bias estimate updated from tracking residuals.

Example disturbance observer:

```math
\dot{\hat{d}} = -k_d \hat{d} + k_d m (\ddot{x}_{meas} - \ddot{x}_{model})
```

then include `-\hat{d}` in the thrust command.

### 5.4 Tracking error bounds
Under bounded disturbance `\|d(t)\| \le \bar{d}` and feasible reference trajectories respecting `v_max`, `a_max`, and jerk bounds, the closed loop should satisfy practical boundedness:

```math
\|e_x(t)\| \le \epsilon_x, \quad \|e_v(t)\| \le \epsilon_v
```

after transient convergence, where `\epsilon_x` and `\epsilon_v` depend on controller gains, estimator latency, and disturbance bound.

Engineering targets:
- RMS position tracking error `< 0.15 m` in nominal warehouse flight,
- peak error `< 0.35 m` during aggressive avoidance,
- yaw tracking error `< 8 deg`.

## 6. Implementation blueprint

### 6.1 ROS 2 packages
- `mission_planner`: coverage graph, TSP ordering, segmentation, RTH
- `vdbfusion_mapping`: depth fusion and ESDF publication
- `ego_planner_ros2`: B-spline optimizer and dynamic obstacle cost plugin
- `trajectory_bridge_px4`: converts spline samples to PX4 trajectory setpoints
- `tracking_controller`: optional outer-loop geometric controller and monitoring

### 6.2 Execution timing
- ESDF update: `10-20 Hz`
- local replanning: `10 Hz`
- trajectory sample publication: `30-50 Hz`
- controller outer loop: `50-100 Hz`
- PX4 inner loop remains on FMU at native high rate

### 6.3 Verification plan
1. Static-map collision tests against shelf corridors.
2. Dynamic obstacle playback with crossing pedestrians/forklifts.
3. Planning latency benchmark on Jetson Orin NX.
4. SITL/HITL warehouse missions including battery-triggered RTH.
5. Real indoor validation with safety pilot and progressive speed limits.

Acceptance criteria:
- 10 Hz replanning sustained,
- zero hard-collision violations in validation scenarios,
- dynamic obstacle near-miss probability below threshold,
- successful mission segmentation and RTH under low battery conditions.
