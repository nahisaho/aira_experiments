# Warehouse Inventory Management Case Study

## Scenario description
This case study assumes a **GPS-denied autonomous drone** using **ROS 2 + PX4** on **Jetson Orin NX 16GB** to perform inventory counting inside a medium-sized warehouse.

### Warehouse profile
- Footprint: **50 m × 30 m × 8 m**
- Storage: **5 aisles** with **4-tier metal shelving** to **6 m** height
- Inventory scale: **~5000 SKUs** with barcode/QR labels on shelf faces, pallets, or tote bins
- Illumination: **300–500 lux** industrial LED with localized shadowing under upper shelves
- Dynamics: **2–3 forklifts** and **5–10 workers** during regular daytime operation
- Flight envelope: **1.5–6.0 m AGL**, typically **0.35–0.80 m** stand-off from labels during scans

## Design assumptions
- Warehouse labels are reasonably maintained and mostly in the **60–120 mm** class for primary shelf labels; smaller labels are handled by closer hover-and-scan passes.
- The drone has a front RGB-D camera for navigation and a scan-focused RGB path using the RealSense D455 RGB stream at **1280×800**.
- Global localization is map-relative (VIO + loop closure + geometric shelf constraints), not GPS-based.
- During business hours, the system favors conservative obstacle separation and yields to humans rather than tightly squeezing past them.

## Mission planning concept
### 1. Pre-flight
- Load the warehouse 3D map, shelf graph, no-fly polygons, and aisle/tier waypoint library from `data/warehouse_layout.yaml`.
- Run health checks for camera, IMU, battery, prop guards, obstacle sensors, PX4 offboard link, and Jetson thermal margin.
- Precompute aisle coverage order to minimize deadheading and ensure the return path always remains within the remaining battery budget.

### 2. Takeoff
- Start from a marked charging/landing pad with a protected takeoff cylinder.
- Climb vertically to **1.8 m AGL**, verify stable VIO quality, then translate to a staging waypoint before entering the first aisle.

### 3. Coverage flight
Use a **boustrophedon (lawn-mower) pattern** with altitude-layered shelf passes:
1. Enter aisle centerline at transit speed.
2. Align laterally to the target shelf face.
3. Slow to scan speed and move along the face.
4. At each bay or waypoint cluster, **hover 1.0–1.5 s** for sharp barcode/QR capture.
5. Change altitude to the next tier and sweep the aisle in the reverse direction.
6. Exit into the cross-aisle, then transition to the next aisle.

### 4. Dynamic replanning
- **Forklift detected**: retreat to nearest pull-out waypoint or hold in cross-aisle until the aisle is clear.
- **Worker detected**: stop or back off to maintain a human-drone separation bubble, then continue after clearance.
- **Temporary obstruction** (ladder, pallet jack, open dock door airflow): mark the shelf segment incomplete and schedule a revisit at the end of the mission or during a quieter interval.

### 5. Battery management
- Split the mission into **15-minute flight segments**.
- Trigger auto-return at **30% battery** or earlier if the planner predicts insufficient energy to complete the next aisle and return.
- Opportunistically recharge between segments and resume from the last completed aisle-tier checkpoint.

### 6. Landing
- Return to the nearest charging pad, switch to precision landing mode, descend using vision-based landing alignment, and dock for charging or battery swap.

## Inventory scanning system
### Sensor and decoding stack
- **Camera**: RealSense D455 RGB at **1280×800** for label reading, depth stream used mainly for stand-off control and obstacle awareness
- **Barcode**: `ZBar` for fast first-pass decode plus a fine-tuned lightweight detector for locating partially occluded/angled labels
- **QR**: OpenCV QR decoder for structured codes and shelf-location markers
- **OCR fallback**: PaddleOCR for damaged or partially torn labels when barcode decoding fails
- **Execution model**: dedicated scanning thread pool and non-blocking queue; scan results never block flight-critical control callbacks

### Why stand-off distance matters
At 1280×800, small 25–50 mm codes are difficult to decode much beyond about **0.2–0.5 m**. For robust warehouse performance, the preferred operating model is:
- **60–80 mm labels**: scan at **0.35–0.60 m** stand-off
- **100 mm+ shelf markers / QR location tags**: scan at **0.6–1.2 m**
- **Damaged or reflective labels**: reduce speed and increase dwell time, or trigger revisit with OCR fallback

## End-to-end autonomy stack for the warehouse
| Layer | Main function | Notes |
|---|---|---|
| PX4 offboard control | Stabilization, failsafes, landing, battery monitoring | Runs independently from high-level scanning logic |
| ROS 2 mission manager | Waypoint sequencing, aisle state machine, WMS transaction batching | Owns mission progress and retries |
| VIO + mapping | GPS-denied localization inside repetitive aisles | Use loop closures, shelf geometry, and fiducial aids if available |
| Dynamic avoidance | Worker/forklift tracking and local replan | Prefer yielding over aggressive bypass |
| Scan pipeline | Barcode / QR / OCR extraction and shelf association | Dedicated threads and bounded queues |
| WMS integration | Inventory reconciliation and exception handling | Upload incrementally, not only at end-of-flight |

## Performance estimates
### Throughput model
A realistic single-drone cycle looks like this:
- **Shelf-face scan speed**: 0.6–0.8 m/s
- **Dwell per scan cluster**: 1.0–1.5 s
- **Average successful reads per cluster**: 3–6 labels depending on bay density
- **Revisit overhead**: 5–8% of mission time for occlusions, glare, and human traffic

### KPI table
| Metric | Expected value | Basis |
|---|---:|---|
| Scan rate | **180–220 items/min** | Derived from hover-and-scan cadence and 3–6 labels per cluster |
| Full warehouse scan | **23–28 min** | ~5000 SKUs / 180–220 items/min plus transit overhead |
| Position accuracy | **±10 cm** | Adequate for shelf-face association and aisle/tier assignment |
| Scanning success per flyover | **98.0–98.8%** | With revisit logic and OCR fallback |
| Missed item rate | **1.2–1.8%** | Mostly due to occlusion, glare, or damaged labels |
| False reading rate | **0.03–0.08%** | Controlled by multi-frame confirmation and WMS plausibility checks |
| Jetson module power | **15–20 W average** | Scan mode in 15–25 W power profiles |

## Shelf association logic
To avoid reading the right label but assigning it to the wrong location, the system should fuse:
1. Drone pose in the warehouse frame
2. Shelf-row polygon and tier height band
3. Camera optical axis and label depth
4. Temporal consistency over 2–3 adjacent frames
5. WMS plausibility constraints (expected SKU family or location ID format)

A decode is only committed when both **visual confidence** and **spatial consistency** pass threshold.

## Operational considerations
### Night operation mode
- Best operational window for a first deployment because it removes most worker/forklift conflicts.
- Increase transit speed slightly and reduce required avoidance detours.
- Use static aisle closure rules for areas under maintenance.

### Multi-drone extension
- Future architecture should partition by aisle groups or height bands rather than free-for-all coverage.
- Share occupancy maps and claimed aisle segments over ROS 2 DDS or a supervisor service.
- Reserve emergency hover points at both ends of each aisle to avoid deadlock.

### WMS integration
- Preferred interface: REST or message-bus transaction containing timestamp, location ID, SKU/label payload, confidence, image snippet reference, and reconciliation state.
- Upload reads incrementally so failed flights do not lose completed inventory data.
- Generate exception queues for unreadable labels, duplicate reads, and location mismatches.

### Regulatory and compliance
- Indoor warehouse operation typically avoids FAA airspace constraints, but site safety rules, insurance requirements, and corporate EHS policies still apply.
- Define operator training, emergency-stop procedures, maintenance logs, prop-guard requirements, and human-notification signage.

### Human-drone safety protocol
- Announce mission start on facility HMI or beacon light.
- Maintain conservative stand-off from workers and stop when a person enters the protected aisle segment.
- Cap lateral speed during worker-adjacent operation and prefer backing out into cross-aisles rather than overtaking.
- Keep prop guards and remote emergency kill available at all times.

## Failure handling strategy
| Failure mode | Immediate action | Recovery |
|---|---|---|
| VIO quality drops in repetitive aisle | Slow down, hold position, switch to marker-assisted relocalization if available | Backtrack to last high-confidence keyframe |
| Barcode decode failures rise | Increase dwell, reduce yaw rate, enable OCR fallback | Schedule localized rescan pass |
| Forklift blocks aisle | Retreat to holding waypoint | Resume after occupancy clears |
| Thermal pressure on Jetson | Lower detector frequency and switch to 15 W mode | Pause mission or return to pad if sustained |
| Battery forecast insufficient | Abort remaining aisle | Return and resume from checkpoint after recharge |

## Why this design is realistic
- The throughput target is deliberately lower than highly optimized commercial pallet-scanning systems because this case study assumes **shelf-face SKU scanning**, not only coarse pallet occupancy.
- The ±10 cm navigation tolerance is achievable for indoor drone inventory work and is sufficient when paired with shelf geometry and multi-frame label confirmation.
- The mission plan is built around **safe, repeatable aisle behavior**, which matters more in warehouses than maximizing raw flight speed.

## References
1. RealSense D455 product page/specifications: https://www.intelrealsense.com/depth-camera-d455/
2. NVIDIA Isaac ROS Visual SLAM documentation: https://nvidia-isaac-ros.github.io/repositories_and_packages/isaac_ros_visual_slam/index.html
3. Corvus Robotics case studies: https://blog.corvus-robotics.com/corvus-robotics-case-studies
4. Gather AI / warehouse drone scanning article: https://www.insidelogistics.ca/products/drone-powered-scanning-for-inferred-case-counting-and-location-occupancy/
5. SupplyChainBrain warehouse drone inventory case article: https://www.supplychainbrain.com/articles/40614-how-a-3pl-eliminated-shrink-with-drone-powered-inventory-counting
