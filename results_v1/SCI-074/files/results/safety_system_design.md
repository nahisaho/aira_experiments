# Safety System and Failsafe Architecture Design

## Scope and assumptions
- System: indoor GPS-denied autonomous drone using ROS 2 + PX4 with onboard perception, mapping, planning, and offboard mission execution.
- Safety objective: prevent collision, loss of control, battery depletion incidents, and geofence violations while maintaining graceful degradation.
- Safety philosophy: layered defense with independent monitors and deterministic escalation paths.

## 1. Multi-level failsafe architecture

### 1.1 Layered safety stack
#### A. Software layer (ROS 2 application level)
Responsibilities:
- dynamic obstacle risk assessment,
- map/planner safety validation,
- watchdog supervision of perception/planning/control nodes,
- mission management and degraded-mode transitions.

Core nodes:
- `safety_supervisor_node`
- `sensor_health_monitor_node`
- `planner_guard_node`
- `battery_safety_node`
- `comm_watchdog_node`
- `flight_state_machine_node`

#### B. Firmware layer (PX4 flight stack)
Responsibilities:
- arming checks,
- estimator validity checks,
- RC/manual override handling,
- offboard loss failsafe,
- geofence backup actions,
- low-level attitude/rate stabilization.

#### C. Hardware layer
Responsibilities:
- motor kill / hardware interlock,
- independent power supervision,
- watchdog reset line,
- battery protection circuitry,
- optional prop guards and collision cage.

Design principle:

```math
Safety = Software\ mitigation \prec Firmware\ containment \prec Hardware\ last\ resort
```

Meaning: higher layers attempt graceful mitigation; lower layers guarantee containment if software fails.

### 1.2 Independence and escalation
Each layer should be able to force a safer state without waiting for upstream recovery:
- ROS 2 can request hover, brake, return, or land.
- PX4 can exit offboard and enter Hold/Land on heartbeat loss.
- Hardware can cut propulsion if commanded by certified external interlock or kill switch.

## 2. Watchdog timers and node supervision

### 2.1 Heartbeat design
Each critical node publishes heartbeat:
- topic: `/health/<node_name>`
- message: timestamp, sequence number, state, optional diagnostics

Critical nodes and nominal periods:
- perception pipeline: `33 ms`
- tracker/predictor: `50-100 ms`
- local planner: `100 ms`
- state estimator bridge: `20-50 ms`
- mission manager: `200 ms`

Supervisor timeout rule:

```math
T_{timeout,i} = \max(3T_i, T_{min})
```

where `T_i` is nominal period and `T_min` is a floor such as `100 ms`.

A node is unhealthy if:

```math
t_{now} - t_{last,i} > T_{timeout,i}
```

### 2.2 Recovery policy
- single missed heartbeat: warning only,
- repeated timeout for noncritical node: isolate subsystem and degrade functionality,
- timeout for critical node (state estimate, planner, obstacle monitor): switch to hover or controlled land depending on altitude and map confidence,
- repeated multi-node failure: trigger PX4 failsafe mode.

## 3. Sensor health monitoring

### 3.1 Monitored metrics
For each sensor stream measure:
- frequency,
- end-to-end latency,
- timestamp jitter,
- dropout ratio,
- field-of-view coverage or valid pixel ratio,
- estimator innovation consistency where applicable.

Examples:
- stereo/depth camera: valid depth fraction in central ROI, exposure status, frame age,
- IMU: sample rate, clipping, vibration level,
- barometer/rangefinder: residual spikes,
- VIO/LiDAR odometry: covariance growth, innovation residuals.

### 3.2 Health score
Define normalized health score:

```math
H_s = w_f h_f + w_l h_l + w_c h_c + w_d h_d
```

where:
- `h_f`: frequency compliance,
- `h_l`: latency compliance,
- `h_c`: coverage/quality,
- `h_d`: dropout score.

Example decision thresholds:
- `H_s >= 0.8`: healthy,
- `0.5 <= H_s < 0.8`: degraded,
- `H_s < 0.5`: failed.

### 3.3 Multi-sensor voting and redundancy
- If D455 depth degrades but VIO and obstacle tracks remain stable, reduce speed and continue.
- If both depth and obstacle localization fail, suspend forward flight and hover/land.
- If state estimation covariance exceeds threshold, block aggressive maneuvers and initiate return or land.

## 4. Battery management and emergency landing

### 4.1 Battery state logic
Monitor:
- state of charge `SOC`,
- terminal voltage under load,
- current and power,
- estimated remaining flight time.

Energy-aware guard condition:

```math
E_{remain} - E_{RTH} - E_{land} - E_{reserve} \ge 0
```

If violated, mission continuation is unsafe.

### 4.2 Battery thresholds
Recommended thresholds:
- **Advisory**: `SOC <= 35%` → no new mission segment
- **Return trigger**: `SOC <= 25%` or reserve violation predicted
- **Emergency land**: `SOC <= 15%`, rapid voltage sag, or battery fault

Use both SOC and under-load voltage to avoid false confidence from flat discharge curves.

### 4.3 Emergency landing selection
Choose landing site by risk score:

```math
J_{land}(s) = \alpha d(s) + \beta o(s) + \gamma g(s) + \delta c(s)
```

where:
- `d(s)`: distance/time to site,
- `o(s)`: obstacle density,
- `g(s)`: ground suitability / flatness,
- `c(s)`: crowd or traffic penalty.

Preferred order:
1. designated home pad,
2. pre-mapped emergency landing zone,
3. local hover and controlled descent area.

## 5. Communication loss handling

### 5.1 Communication channels
- ROS 2 intra-system communications,
- companion computer to PX4 offboard heartbeat,
- optional operator uplink / RC override.

### 5.2 Loss scenarios and actions
1. **ROS 2 graph degradation but PX4 alive**
   - if planner heartbeats lost: PX4 Hold or Brake.
2. **Offboard link loss to PX4**
   - PX4 exits Offboard after timeout and enters Hold / Return / Land according to configuration.
3. **Operator uplink loss**
   - continue autonomously only if mission is approved for disconnected operation and all health scores remain green; otherwise return/land.

Recommended PX4 configuration concept:
- short offboard loss timeout (`0.3-0.5 s`),
- Hold first if map and localization are valid,
- Land if loss persists or localization degrades.

## 6. Collision-imminent response hierarchy

### 6.1 Risk metrics
Compute at each cycle:
- minimum predicted distance `d_min`,
- time-to-collision `TTC_min`,
- collision probability `P_coll^max`,
- stopping distance estimate `d_stop`.

Stopping distance estimate for current speed `v` and max deceleration `a_{brake}`:

```math
d_{stop} = \frac{v^2}{2a_{brake}} + v t_{latency}
```

### 6.2 Escalation policy
1. **Monitor**: track only, no maneuver.
2. **Avoid**: local replan, reduce speed cap.
3. **Brake/Hover**: command zero-velocity or hover setpoint.
4. **Emergency stop / climb / descend**: execute prevalidated evasive primitive if braking alone is insufficient.
5. **Emergency land**: if no safe flight corridor remains.
6. **Hardware kill**: only if uncontrolled thrust or imminent catastrophic contact cannot be mitigated otherwise.

Decision example:
- if `d_min > 5 m`: monitor,
- if `1.5 < d_min <= 5 m` or `P_coll^max > 0.2`: avoid,
- if `d_min <= max(1.5 m, d_stop)` or `TTC_min < 0.7 s`: brake/hover,
- if attitude/control diverges simultaneously: land or kill depending on altitude and environment.

## 7. Geofence enforcement

### 7.1 Geofence representation
Use layered geofences:
- mission soft geofence in ROS 2 planner,
- firmware geofence in PX4 as backup,
- optional hard exclusion via safety MCU/interlock for special zones.

Represent zones as 3D polygon prisms or voxel masks. Add inflation margin:

```math
G_{inflated} = G \oplus B(r_{margin})
```

where `\oplus` denotes Minkowski sum.

### 7.2 Enforcement actions
- approaching boundary: reduce speed and bias planner inward,
- predicted crossing within horizon: immediate replan,
- actual crossing due to disturbance or estimator jump: PX4 geofence action (Hold/Land) and supervisor event.

## 8. Flight state machine

### 8.1 Primary nominal states
```text
IDLE → ARMED → TAKEOFF → MISSION → RETURN → LAND → DISARMED
```

### 8.2 Error and degraded states
Additional states:
- `HOVER_HOLD`
- `EMERGENCY_BRAKE`
- `EMERGENCY_LAND`
- `FAULT`
- `MANUAL_OVERRIDE`

### 8.3 State semantics
| State | Entry conditions | Exit conditions / transitions |
|---|---|---|
| `IDLE` | powered, disarmed, self-check running | `ARMED` if operator command and all pre-arm checks pass |
| `ARMED` | motors armed, ready for takeoff | `TAKEOFF` on takeoff command; `FAULT` on failed checks |
| `TAKEOFF` | autonomous ascent to safe hover altitude | `MISSION` on altitude/pose stabilization; `HOVER_HOLD` on transient issue; `EMERGENCY_LAND` on severe fault |
| `MISSION` | executing mission segments | `RETURN` on mission complete, low battery, or operator recall; `EMERGENCY_BRAKE` on collision threat; `HOVER_HOLD` on planner/perception degradation |
| `RETURN` | navigating to home or recovery point | `LAND` on arrival; `EMERGENCY_LAND` if energy/perception insufficient |
| `LAND` | controlled descent and touchdown | `DISARMED` after touchdown and spin-down; `EMERGENCY_BRAKE` if obstacle enters landing path |
| `DISARMED` | landed and safe | `IDLE` after reset/ready cycle |
| `HOVER_HOLD` | temporary safe stop | `MISSION` or `RETURN` if health recovers; `EMERGENCY_LAND` if timeout or worsening fault |
| `EMERGENCY_BRAKE` | imminent collision | `HOVER_HOLD`, `RETURN`, or `EMERGENCY_LAND` depending on post-brake health |
| `EMERGENCY_LAND` | controlled immediate landing | `DISARMED` after touchdown; `FAULT` if landing impossible |
| `FAULT` | unrecoverable safety failure | `MANUAL_OVERRIDE` or `DISARMED` |
| `MANUAL_OVERRIDE` | pilot/RC takes control | `IDLE`/`DISARMED` after manual recovery |

### 8.4 Transition guards
Examples:
- `IDLE -> ARMED` only if:

```math
H_{state\ estimator} > 0.8 \land H_{depth} > 0.8 \land SOC > SOC_{min\_launch} \land geofence\_loaded = true
```

- `MISSION -> RETURN` if:

```math
mission\_complete \lor SOC \le 25\% \lor comm\_policy\_requires\_return
```

- `MISSION -> EMERGENCY_BRAKE` if:

```math
d_{min} \le \max(1.5, d_{stop}) \lor TTC_{min} < 0.7
```

- `HOVER_HOLD -> EMERGENCY_LAND` if hold timeout exceeded or localization remains degraded.

### 8.5 State machine implementation notes
- Implement state machine as a deterministic ROS 2 lifecycle or explicit hierarchical state machine.
- Every transition should emit an event log and reason code.
- Safety transitions must preempt mission transitions.
- Manual override must be accepted from any autonomous state.

## 9. Verification and validation

### 9.1 Required tests
1. watchdog timeout injection,
2. sensor dropout and latency injection,
3. low battery replay and RTH trigger test,
4. offboard link loss in PX4 SITL/HITL,
5. dynamic obstacle emergency brake scenario,
6. geofence breach attempt in simulation,
7. landing-site selection with blocked home pad.

### 9.2 Acceptance targets
- detection of critical node timeout within configured threshold,
- safe hover or land after planner/perception failure,
- no mission continuation after reserve-energy violation,
- deterministic state transition traces with reason codes,
- PX4 backup failsafe successfully triggered when companion stack is unavailable.

## 10. Recommended implementation artifacts
- `safety_msgs`: heartbeat, health status, transition event, failsafe command
- `safety_supervisor`: central arbitration and event logger
- `flight_state_machine`: nominal/degraded/error state coordination
- `px4_safety_bridge`: explicit mapping from supervisor commands to PX4 modes/actions
- `safety_test_suite`: fault injection launch files and scenario tests
