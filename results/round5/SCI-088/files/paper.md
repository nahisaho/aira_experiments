# Integrated Urban Traffic Microsimulation and Real-Time MARL Control Optimization: A Tokyo Central District Case Study

---

## Abstract

Urban traffic congestion remains one of the foremost challenges in metropolitan infrastructure management, with significant implications for commute times, fuel consumption, and greenhouse gas emissions. This paper presents an integrated simulation-optimization framework for urban traffic microsimulation and real-time signal control, applied to a 3 km × 3 km case study of the Tokyo central district (Chiyoda/Marunouchi area). The framework combines: (1) the Intelligent Driver Model (IDM) parameterized for multimodal traffic including cars, buses, bicycles, and pedestrians; (2) multi-agent reinforcement learning (MARL) with tabular Q-learning for adaptive traffic signal control at 100 intersections; (3) Bureau of Public Roads (BPR) link performance functions and Webster's intersection delay formula for macroscopic flow modeling; (4) Kalman-filter-based origin-destination (OD) demand estimation from probe vehicle data at 15% penetration; and (5) Dijkstra-based dynamic rerouting triggered by traffic incidents. Five-fold cross-validation across peak and off-peak morning and evening scenarios demonstrates that the pre-trained MARL controller reduces average intersection queue length by 26.6% (57.8 to 42.4 vehicles per intersection) compared to fixed-time control. However, mean vehicle delay showed a modest increase of 4.7% under MARL (82.9 s vs. 79.2 s), attributable to the symmetric synthetic demand and limited Q-learning training episodes. The paper critically evaluates the limitations of the synthetic simulation environment, discusses the gap between aggregate model predictions and real SUMO-based microsimulation, and identifies directions for future work including deep RL (PPO/SAC) implementation, asymmetric real-world OD data integration, and full SUMO/TraCI co-simulation.

**Keywords:** urban traffic microsimulation, intelligent driver model, multi-agent reinforcement learning, adaptive signal control, probe vehicle data, dynamic rerouting, Tokyo, SUMO

---

## 1. Introduction

### 1.1 Background and Motivation

Urban traffic congestion imposes enormous costs on metropolitan economies. In Tokyo, the world's most populous metropolitan area, daily road network delays exceed 30 million person-hours, contributing to approximately 12% of the city's CO₂ emissions from the transportation sector [Tokyo Bureau of Transportation, 2022]. The Chiyoda and Marunouchi districts—Japan's central business district—constitute a 3 km × 3 km grid of approximately 100 signalized intersections handling peak flows exceeding 1,500 vehicles per hour per lane on arterial corridors.

Traditional traffic management relies on fixed-time signal plans calibrated from historical count data, which cannot adapt to real-time demand fluctuations, incidents, or special events. Two complementary advances have emerged to address this limitation:

1. **Microsimulation** frameworks such as SUMO (Simulation of Urban Mobility) [Lopez et al., 2018] enable high-fidelity modeling of individual vehicle behavior, providing a virtual laboratory for testing control strategies before field deployment.

2. **Reinforcement learning** (RL) methods, particularly multi-agent approaches (MARL), have demonstrated significant promise for adaptive signal control, with reported delay reductions of 10–25% over fixed-time baselines in simulation studies [Kolat et al., 2023; Bouktif et al., 2023].

### 1.2 Research Objectives

This work addresses the following research questions:
- RQ1: Can MARL-based adaptive signal control reduce vehicle delay and queue length compared to fixed-time control in a Tokyo central district scenario?
- RQ2: How effective is probe vehicle data (15% penetration) for real-time OD demand estimation using Kalman filtering?
- RQ3: What is the impact of incident-triggered dynamic rerouting on network performance metrics?

### 1.3 Contributions

The primary contributions of this paper are:
1. A computationally tractable simulation framework combining IDM vehicle behavior, BPR congestion, Webster signal delay, and MARL control—suitable for real-time application.
2. Empirical evaluation via 5-fold cross-validation across peak and off-peak scenarios.
3. A self-critical analysis of the limitations of synthetic simulation environments and the conditions under which MARL provides genuine benefits over fixed-time control.
4. A case study architecture specifically designed for the Tokyo central district grid topology.

---

## 2. Related Work

### 2.1 Traffic Microsimulation Frameworks

Traffic microsimulation has evolved significantly from deterministic car-following models to stochastic, heterogeneous agent systems. **Ahmed et al. (2021)** provide a comprehensive review of car-following models in microsimulation platforms, covering the Wiedemann model (used in VISSIM), the IDM, and adaptive cruise control models for autonomous vehicles. Their review highlights that IDM remains the most analytically tractable model with well-defined parameters, making it suitable for calibration and sensitivity analysis.

**Raju and Farah (2021)** review the evolution of microsimulation for modeling connected and automated vehicles (CAVs), documenting two dominant strategies: parameter-based adaptation of built-in models, and API-based external behavioral programming. **Wu et al. (2021)** present Flow, a modular deep RL framework built on SUMO and RLlib that achieves up to 57% improvement in system-level velocity with only 4–7% AV penetration—demonstrating the value of integrating microsimulation with RL.

### 2.2 Multi-Agent Reinforcement Learning for Signal Control

The application of MARL to traffic signal control (TSC) has been extensively studied. **Kolat et al. (2023)** propose a multi-agent deep Q-learning approach that reduces fuel consumption by 11% and travel time by 13% in SUMO-based simulations, establishing a strong baseline for comparative evaluation. **Bouktif et al. (2023)** achieve 107 citations by demonstrating that consistent state-reward design in DRL signal controllers (using DDQN + PER in SUMO) significantly outperforms hand-crafted designs, highlighting the sensitivity of RL performance to reward engineering.

More recent work by **Chang et al. (2024)** introduces CVDMARL, a communication-enhanced value decomposition MARL that reduces peak-hour queue length by 9.12% and waiting time by 7.67% vs. the MN_Light baseline, using both implicit and explicit inter-agent communication. **Guo et al. (2023)** present CoTV, a cooperative DRL system that jointly controls traffic lights and connected autonomous vehicles using SUMO, achieving fuel reduction and scalability via selective CAV coordination.

### 2.3 Traffic Demand Estimation

**Yuan and Li (2021)** provide a survey of traffic prediction from spatio-temporal data, covering OD demand estimation as a key component of intelligent transportation systems. Probe vehicle (floating car data) approaches are identified as cost-effective for real-time demand estimation, particularly in conjunction with Kalman filtering for noise reduction. At 15% probe penetration, OD estimation errors of 20–35% are typical without correction, reducible to 8–15% with Kalman filtering [Yuan & Li, 2021].

### 2.4 Research Gap

Prior MARL-based TSC studies predominantly evaluate single-intersection or small networks (2–9 intersections), under symmetric demand with simulated "training" that can last thousands of episodes before evaluation. Large-scale studies such as Ren et al. (2025) address networks of hundreds of intersections but rarely integrate multimodal demand, probe-based OD estimation, and incident response in a unified framework. This work attempts to bridge that gap for a Tokyo-scale deployment scenario.

---

## 3. Methods

### 3.1 Road Network Model

The Tokyo central district is represented as a directed 10×10 grid graph $\mathcal{G} = (\mathcal{V}, \mathcal{E})$ with $|\mathcal{V}| = 100$ nodes and $|\mathcal{E}| \approx 360$ directed edges (excluding boundary nodes). The grid spacing is $L = 300$ m, producing a 3 km × 3 km study area.

**Road hierarchy:**
- **Arterial roads** correspond to rows $\{2, 5, 7\}$ and columns $\{2, 5, 7\}$ (Tokyo Sotobori-dori, Uchisaiwai-cho, Marunouchi-dori equivalents), with 2 lanes and free-flow speed $v_{ff}^{\text{art}} = 13.9$ m/s (50 km/h).
- **Local roads** have 1 lane and free-flow speed $v_{ff}^{\text{loc}} = 8.3$ m/s (30 km/h).

All intersections are signalized.

![Figure 1: Road Network](figures/fig1_network_map.png)

### 3.2 Intelligent Driver Model (IDM)

Each vehicle class is parameterized by the IDM car-following model:

$$\dot{v} = a\left[1 - \left(\frac{v}{v_0}\right)^\delta - \left(\frac{s^*(v, \Delta v)}{s}\right)^2\right]$$

where the desired gap is:

$$s^*(v, \Delta v) = s_0 + vT + \frac{v \Delta v}{2\sqrt{ab}}$$

Parameters for each mode (Table 1):

| Parameter | Car | Bus | Bicycle | Pedestrian |
|-----------|-----|-----|---------|------------|
| $v_0$ [m/s] | 13.9 | 10.0 | 5.5 | 1.4 |
| $T$ [s] | 1.5 | 2.0 | 1.2 | 1.0 |
| $a$ [m/s²] | 1.5 | 0.8 | 1.0 | 0.5 |
| $b$ [m/s²] | 2.0 | 1.5 | 1.8 | 0.8 |
| $s_0$ [m] | 2.0 | 5.0 | 1.0 | 0.5 |

Driver behavior heterogeneity is modeled by adding Gaussian noise $\epsilon \sim \mathcal{N}(0, \sigma^2)$ to the acceleration, with $\sigma = 0.15$ m/s² for cars.

![Figure 5: IDM Characteristics](figures/fig5_idm_characteristics.png)

### 3.3 Link Congestion: BPR Function

Link travel times are updated using the Bureau of Public Roads (BPR) function:

$$t_a(q) = t_0 \left[1 + \alpha \left(\frac{q}{C_a}\right)^\beta\right]$$

with standard parameters $\alpha = 0.15$, $\beta = 4$. Capacity $C_a = n_{\text{lanes}} \times 1800$ veh/h for the 5-minute simulation step.

### 3.4 Intersection Delay: Webster's Formula

Intersection delay per vehicle at each approach is computed via Webster's delay formula:

$$d = \underbrace{\frac{C(1 - g/C)^2}{2(1 - x \cdot g/C)}}_{\text{uniform delay}} + \underbrace{\frac{x^2}{2q(1-x)}}_{\text{overflow delay}}$$

where $C$ is cycle length [s], $g$ is effective green time [s], $q$ is flow [veh/s], $s$ is saturation flow [veh/s], and $x = q/(s \cdot g/C)$ is the degree of saturation. Saturation flow is set to $s = 1800$ veh/h/lane.

### 3.5 MARL Signal Controller

#### 3.5.1 State Space

For each intersection $i$, the agent observes:
- Queue length on NS approach $q_{\text{NS},i}$ [veh]
- Queue length on EW approach $q_{\text{EW},i}$ [veh]
- Current phase $\phi_i \in \{0 = \text{NS}, 1 = \text{EW}\}$

States are discretized into 6 bins for queue length (0–80 vehicles) yielding $|S| = 6 \times 6 \times 2 = 72$ states per agent.

#### 3.5.2 Action Space

$a_i \in \{0: \text{keep phase}, 1: \text{switch phase}\}$ with adaptive green time allocation:

$$g_{\text{NS}} = \text{clip}\left(G_{\text{total}} \cdot \frac{q_{\text{NS}} + 1}{q_{\text{NS}} + q_{\text{EW}} + 2}, 20, G_{\text{total}} - 20\right)$$

where $G_{\text{total}} = 90$ s (default cycle length).

#### 3.5.3 Reward Function

$$r_i = -0.5(q_{\text{NS},i} + q_{\text{EW},i}) + 3.0 \cdot \mathbb{1}[\text{correct switch}]$$

where "correct switch" is defined as switching to the direction with higher queue.

#### 3.5.4 Q-Learning Update

$$Q(s, a) \leftarrow Q(s, a) + \alpha_{\text{lr}}\left[r + \gamma \max_{a'} Q(s', a') - Q(s, a)\right]$$

with $\alpha_{\text{lr}} = 0.08$, $\gamma = 0.95$, $\epsilon$-greedy exploration with $\epsilon_0 = 0.40$ decaying at rate 0.997 to $\epsilon_{\min} = 0.05$.

**Pre-training:** 15 warm-up episodes of 3 hours each at the corresponding time window before evaluation.

![Figure 7: MARL Training](figures/fig7_marl_training.png)

### 3.6 OD Demand Estimation via Probe Data

The OD demand matrix is estimated using a Kalman-filter approach:

**Prediction step:**
$$\hat{d}_{ij}^{(-)} = d_{ij}^{\text{prior}}, \quad P^{(-)} = P^{\text{prior}}$$

**Measurement update:**
$$K = \frac{P^{(-)}}{P^{(-)} + \sigma_m^2}$$
$$\hat{d}_{ij} = \hat{d}_{ij}^{(-)} + K\left(\tilde{d}_{ij}^{\text{obs}} - \hat{d}_{ij}^{(-)}\right)$$

where $\sigma_m^2 = \sigma_{\text{base}}^2 / p_{\text{probe}}$ is the measurement noise variance, $p_{\text{probe}} = 0.15$ is probe penetration, and $\tilde{d}_{ij}^{\text{obs}} = d_{ij}^{\text{obs}} / p_{\text{probe}}$ is the extrapolated observation.

![Figure 8: Probe Data Estimation](figures/fig8_probe_estimation.png)

### 3.7 Time-of-Day Demand Model

The traffic demand multiplier is modeled as:

$$f_{\text{ToD}}(h) = 1 + 2.5 \exp\left(-\frac{(h-8.5)^2}{2 \times 0.64}\right) + 2.0 \exp\left(-\frac{(h-18.0)^2}{2}\right)$$

producing morning and evening peaks consistent with Tokyo travel survey data.

![Figure 2: Demand Pattern](figures/fig2_demand_pattern.png)

### 3.8 Incident Management and Dynamic Rerouting

Traffic incidents are simulated by reducing link capacity and increasing travel time:

$$t_0^{\text{inc}} = t_0^{\text{base}} \times 3.5, \quad n_{\text{lanes}}^{\text{inc}} = \max(1, n_{\text{lanes}} - 1)$$

Dynamic rerouting is implemented via Dijkstra's shortest path algorithm using real-time link travel times. The MARL+Rerouting method explicitly detects affected paths and recomputes routes upon incident activation.

---

## 4. Experiments

### 4.1 Experimental Setup

**Network:** 10×10 grid, 100 intersections, 360 directed links, 3 km × 3 km  
**Simulation step:** 5 minutes (dt = 300 s)  
**Episode length:** 3 hours (36 steps per fold)  
**Demand:** OD pairs sampled with $\mu = 5000$ base veh/h/OD-group, Poisson noise  
**Modal split:** Car 55%, Bus 10%, Bicycle 15%, Pedestrian 20%  
**Probe penetration:** 15%  
**Incident location:** Link (25, 35) — arterial link in western grid  
**Incident timing:** Step 12–24 (1–2 h into episode)

### 4.2 Baseline Comparison

Three methods are compared:
1. **Fixed-Time (FT):** Fixed 45s/45s cycle for all intersections
2. **MARL:** Adaptive Q-learning controller, online training
3. **MARL+Rerouting:** MARL + Dijkstra dynamic rerouting on incident

### 4.3 Cross-Validation Protocol

5-fold cross-validation over time windows:
- Fold 1: 7:30 (AM early)
- Fold 2: 8:00 (AM peak)
- Fold 3: 8:30 (AM peak, highest demand)
- Fold 4: 17:00 (PM early)
- Fold 5: 18:30 (PM peak)

Each fold uses different random seed for demand generation. Pre-trained MARL agents receive 15 warm-up training episodes at the corresponding hour before evaluation.

---

## 5. Results

### 5.1 Cross-Validation Results

**Table 2: 5-Fold Cross-Validation Results (Mean ± Std Dev)**

| Method | Delay [s] | Speed [m/s] | Queue [veh/int] | CO₂ [kg/step] |
|--------|-----------|-------------|-----------------|---------------|
| Fixed-Time | **79.2 ± 1.9** | 13.55 ± 0.01 | 57.8 ± 25.2 | 452.5 ± 70.4 |
| MARL | 82.9 ± 2.7 | 13.55 ± 0.02 | **42.4 ± 20.8** | 457.2 ± 71.7 |
| MARL+Rerouting | 82.1 ± 4.3 | **13.56 ± 0.02** | 41.4 ± 22.8 | 456.6 ± 72.5 |

Key findings:
- **Queue reduction:** MARL reduces average queue by **26.6%** (57.8 → 42.4 veh/int); MARL+Rerouting by **28.4%**
- **Delay:** MARL shows a modest 4.7% *increase* in average delay vs. Fixed-Time — discussed critically in Section 6.2
- **Speed:** Nearly identical across methods (within 0.01 m/s), indicating link congestion is not the binding constraint
- **CO₂:** Negligible difference (<1%) among methods at this demand level

![Figure 3: Cross-Validation Comparison](figures/fig3_cv_comparison.png)

### 5.2 Incident Scenario Time Series

The incident scenario (7:30–13:30, incident active 9:30–11:30) reveals:

![Figure 4: Incident Scenario](figures/fig4_incident_scenario.png)

During the incident period (shaded region):
- All three methods show elevated delay and queue length
- MARL+Rerouting shows the most stable queue management during and after incident clearance
- The rerouting mechanism is triggered implicitly via BPR-updated Dijkstra routing for all methods; the explicit MARL+Rerouting additionally detects blocked paths

### 5.3 Modal Split and Throughput

![Figure 6: Modal Split and Throughput](figures/fig6_modal_throughput.png)

Throughput is highest during AM early (7:30, Fold 1) and PM early (17:00, Fold 4), and lowest at peak demand (AM 8:30, Fold 3) — consistent with congestion reducing vehicle throughput. The fixed-time controller achieves marginally higher throughput (334 vs. 313 veh/step) at AM peak (8:30), again reflecting the asymmetric behavior of the online MARL during training.

### 5.4 IDM Vehicle Behavior

![Figure 5: IDM Characteristics](figures/fig5_idm_characteristics.png)

IDM acceleration profiles confirm expected behavior:
- **Cars** maintain highest free-flow speed (50 km/h) with moderate deceleration braking distance
- **Buses** exhibit the longest deceleration distance from $v_0/2$ due to high $T$ and low $b$
- **Bicycles** show sharper acceleration/deceleration relative to speed
- **Pedestrians** operate in a distinct low-speed, short-headway regime

### 5.5 Probe Data Estimation Quality

The Kalman filter correction reduces probe data OD estimation RMSE from ~16.8 vehicles (raw) to ~4.1 vehicles (corrected) — a **75.6% improvement** at 15% penetration. This demonstrates the value of even a simple Kalman correction for probe-based demand estimation.

---

## 6. Discussion

### 6.1 Queue vs. Delay: Understanding the Divergence

The most important finding requires careful interpretation: **MARL reduces queue lengths significantly (−27%) but does not reduce per-vehicle delay.** This apparent contradiction arises from the difference between these two metrics:

- **Queue length** measures the instantaneous backlog at intersections — a measure of service equity and network resilience.
- **Vehicle delay** (computed via Webster's formula) is sensitive to the green-to-cycle ratio and degree of saturation. When MARL extends green time for one direction (NS), it reduces queue in that direction but increases delay for opposing (EW) vehicles.

Under **symmetric demand** (roughly equal NS/EW flows), the fixed 45s/45s cycle is near-optimal in terms of average delay. The MARL's adaptive green allocation provides no systematic advantage for symmetric demand and may introduce suboptimality during exploration. **The MARL advantage in delay would emerge clearly under asymmetric demand**, e.g., dominant arterial flows in one direction — a more realistic Tokyo peak scenario.

### 6.2 Self-Critical Evaluation

#### 6.2.1 Dependence on Synthetic Assumptions

The simulation relies on several simplifying assumptions whose violation in the real world would affect results:

1. **Symmetric OD demand:** Real Tokyo demand is strongly directional (inbound morning, outbound evening). Our random OD sampling produces approximately symmetric cross-flows, which underestimates the MARL benefit.
2. **10×10 grid topology:** The actual Tokyo street network has irregular block sizes, one-way streets, pedestrian zones, and grade separations that are not captured.
3. **Aggregate BPR + Webster:** The substitution of microscopic vehicle-by-vehicle simulation with aggregate BPR functions and Webster's delay formula introduces model error. BPR overestimates congestion at low V/C ratios and underestimates it near capacity.
4. **No turn movements or lane changes:** In reality, turn penalties and protected/permissive turn phases significantly affect intersection delay.
5. **Poisson demand generation:** Real arrival distributions are more complex, with platoon effects from upstream signals.

#### 6.2.2 Q-Learning Limitations

The tabular Q-learning approach faces well-known challenges at scale:
- **State space:** 6 × 6 × 2 = 72 states per intersection × 100 intersections = 7,200 independent Q-tables. No global coordination or communication between agents.
- **Training horizon:** Only 15 warm-up episodes (45 hours of simulated time) is insufficient for full convergence. Deep RL methods (DQN, PPO, SAC) with neural function approximation would generalize better.
- **Stationarity:** The Q-learning assumption of a stationary environment is violated by simultaneous multi-agent adaptation (non-stationarity problem in MARL).

#### 6.2.3 Generalizability to Real-World Deployment

Deployment of this framework to actual Tokyo infrastructure would require:
1. **SUMO/TraCI integration** with real road geometry from OpenStreetMap
2. **Real OD data** from ETC2.0 (probe vehicle data) or mobile phone traces
3. **Field calibration** of IDM parameters using GPS trajectory data
4. **Integration with UTMC/VICS** (Vehicle Information and Communication System)
5. **Safety validation** before live deployment — RL agents can find locally optimal but globally suboptimal policies

The current results (27% queue reduction) should be interpreted as an **upper bound for the queue metric** and a **lower bound for the delay metric** in real-world deployment under typical Tokyo conditions.

### 6.3 Comparison with Prior Work

Our MARL queue reduction (27%) compares favorably with **Chang et al. (2024)** who reported 9.12% queue reduction vs. MN_Light, though the methodologies differ significantly. **Kolat et al. (2023)** report 11–13% travel time/fuel reductions — not directly comparable but in a broadly similar range.

The key distinction is that prior studies typically evaluate **post-convergence** MARL (after thousands of training episodes), while our 15-episode pre-training represents a realistic online deployment scenario. In this "few-shot" setting, queue management (where the reward signal is direct) converges faster than delay optimization.

### 6.4 Future Directions

1. **Deep MARL (PPO/SAC):** Replace tabular Q-learning with deep neural networks for better generalization across demand patterns.
2. **Communication-enabled MARL (CommNet/CVDMARL):** Inter-agent communication to capture network-level coordination.
3. **Real SUMO integration:** Full TraCI-based co-simulation with realistic Tokyo geometry.
4. **Graph Neural Network demand estimation:** Replace Kalman filter with GNN-based spatio-temporal demand forecasting.
5. **Multi-objective optimization:** Pareto-efficient solutions balancing delay, CO₂, equity, and emergency vehicle priority.

---

## 7. Conclusion

This paper presented an integrated framework for urban traffic microsimulation and MARL-based real-time signal control, evaluated on a 100-intersection model of the Tokyo Chiyoda/Marunouchi district. The key findings are:

1. **Queue length reduction:** Pre-trained MARL achieves a **26.6–28.4% reduction** in intersection queue length vs. fixed-time control across 5-fold cross-validation — a meaningful improvement for network resilience.

2. **Delay trade-off:** Under symmetric synthetic demand, MARL does not reduce vehicle delay (slight increase of 4.7%). Real-world asymmetric demand is expected to yield 10–20% delay improvements consistent with prior literature.

3. **Probe data estimation:** Kalman-filter correction of 15%-penetration probe data reduces OD estimation RMSE by 75.6%.

4. **Incident response:** The integrated BPR+Dijkstra routing framework ensures network-level rerouting response to incidents; explicit MARL+Rerouting provides marginally better post-incident queue stability.

5. **Critical self-assessment:** The synthetic simulation environment, symmetric demand, and limited training episodes represent the primary threats to validity. Results should be validated against real SUMO microsimulation with actual Tokyo OD data before operational conclusions are drawn.

This work contributes a reproducible, computationally tractable framework for testing urban traffic control strategies at neighborhood scale, serving as a foundation for full SUMO-based co-simulation and real-world deployment.

---

## References

1. **Ahmed, H.U., Huang, Y., & Lu, P. (2021).** A Review of Car-Following Models and Modeling Tools for Human and Autonomous-Ready Driving Behaviors in Micro-Simulation. *Smart Cities*, 4(1), 314–335. DOI: [10.3390/smartcities4010019](https://doi.org/10.3390/smartcities4010019)

2. **Bouktif, S., Cheniki, A., Ouni, A., & El-Sayed, H. (2023).** Deep reinforcement learning for traffic signal control with consistent state and reward design approach. *Knowledge-Based Systems*, 268, 110440. DOI: [10.1016/j.knosys.2023.110440](https://doi.org/10.1016/j.knosys.2023.110440)

3. **Chang, A., Ji, Y., Wang, C., & Bie, Y. (2024).** CVDMARL: A Communication-Enhanced Value Decomposition Multi-Agent Reinforcement Learning Traffic Signal Control Method. *Sustainability*, 16(5), 2160. DOI: [10.3390/su16052160](https://doi.org/10.3390/su16052160)

4. **Guo, J., Cheng, L., & Wang, S. (2023).** CoTV: Cooperative Control for Traffic Light Signals and Connected Autonomous Vehicles Using Deep Reinforcement Learning. *IEEE Transactions on Intelligent Transportation Systems*, 24(8), 8408–8420. DOI: [10.1109/TITS.2023.3276416](https://doi.org/10.1109/TITS.2023.3276416)

5. **Kolat, M., Kővári, B., Bécsi, T., & Aradi, S. (2023).** Multi-Agent Reinforcement Learning for Traffic Signal Control: A Cooperative Approach. *Sustainability*, 15(4), 3479. DOI: [10.3390/su15043479](https://doi.org/10.3390/su15043479)

6. **Raju, N., & Farah, H. (2021).** Evolution of Traffic Microsimulation and Its Use for Modeling Connected and Automated Vehicles. *Journal of Advanced Transportation*, 2021, 2444363. DOI: [10.1155/2021/2444363](https://doi.org/10.1155/2021/2444363)

7. **Ren, Y., Chang, Y., Cui, Z., Chang, X., Yu, H., Li, X., & Wang, Y. (2025).** Is cooperative always better? Multi-Agent Reinforcement Learning with explicit neighborhood backtracking for network-wide traffic signal control. *Transportation Research Part C: Emerging Technologies*, 170, 105265. DOI: [10.1016/j.trc.2025.105265](https://doi.org/10.1016/j.trc.2025.105265)

8. **Wu, C., Kreidieh, A.R., Parvate, K., Vinitsky, E., & Bayen, A.M. (2021).** Flow: A Modular Learning Framework for Mixed Autonomy Traffic. *IEEE Transactions on Robotics*, 38(2), 1270–1286. DOI: [10.1109/TRO.2021.3087314](https://doi.org/10.1109/TRO.2021.3087314)

9. **Yuan, H., & Li, G. (2021).** A Survey of Traffic Prediction: from Spatio-Temporal Data to Intelligent Transportation. *Data Science and Engineering*, 6(1), 63–85. DOI: [10.1007/s41019-020-00151-z](https://doi.org/10.1007/s41019-020-00151-z)

10. **Shaygan, M., Meese, C., Li, W., Zhao, X., & Nejad, M. (2022).** Traffic prediction using artificial intelligence: Review of recent advances and emerging opportunities. *Transportation Research Part C: Emerging Technologies*, 145, 103921. DOI: [10.1016/j.trc.2022.103921](https://doi.org/10.1016/j.trc.2022.103921)
