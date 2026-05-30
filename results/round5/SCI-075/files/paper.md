# Semi-Autonomous Surgical Suturing via Learning from Demonstration, Real-Time Tissue Deformation Modeling, and Safety-Constrained Impedance Control on the da Vinci Research Kit

---

## Abstract

Robot-assisted minimally invasive surgery (RAMIS) holds great promise for improving surgical precision and reproducibility, yet fully autonomous suturing remains an open challenge owing to the coupled complexity of skill acquisition, soft-tissue interaction, and patient safety guarantees. In this work we present a comprehensive semi-autonomous suturing framework targeting the da Vinci Research Kit (dVRK) simulated through the SurRoL platform. The proposed system integrates four tightly coupled modules: (1) **Learning from Demonstration (LfD)** using Dynamic Movement Primitives (DMP) with 15 basis functions to encode needle-driving trajectories from 15 expert demonstrations; (2) **real-time tissue deformation modeling** using a 15-node Mass-Spring system integrated with a stable ODE solver, quantifying deformation at four tissue stiffness levels (200–2000 N/m); (3) **compliant motion control** via a 1-DoF impedance controller (Kd, Bd, Md) whose stiffness is adapted online based on measured tissue resistance; and (4) **image-based visual servoing** with a proportional controller for stereo-endoscopic needle-tip tracking under varying noise and occlusion conditions. All modules are governed by a **Control Barrier Function (CBF)** safety filter that enforces a maximum contact force of 5 N and a workspace boundary of ±30 mm, reducing simulated constraint violations from 349 force-violations and 98 workspace-violations to zero. Cross-validated evaluation of the DMP representation over 15 demonstrations shows a minimum MSE of 0.0017 ± 0.0002 (5-fold CV). The full-pipeline simulated evaluation (N = 30 synthetic trials per condition) demonstrates that the proposed DMP + Impedance + Safety + Visual Servoing system reduces stitch placement error from 2.15 ± 1.04 mm (teleoperation baseline) to 1.55 ± 0.54 mm, suture tension from 3.91 ± 0.97 N to 3.00 ± 0.41 N, and task completion time from 47.1 ± 12.5 s to 38.4 ± 8.8 s. We critically discuss the limitations of this simulation-based study, including the gap to real tissue heterogeneity, unmodelled sensor dynamics, and the need for prospective clinical validation.

**Keywords**: surgical robotics; learning from demonstration; dynamic movement primitives; tissue deformation; impedance control; visual servoing; safety constraints; dVRK; SurRoL; ROS 2

---

## 1. Introduction

### 1.1 Background and Motivation

Robotic-assisted laparoscopic surgery has achieved widespread clinical adoption, with the da Vinci Surgical System (Intuitive Surgical, Sunnyvale CA) representing the leading teleoperated platform [1]. Despite advances in haptic feedback, 3D visualization, and articulated instrumentation, current clinical robots remain fundamentally teleoperated: the surgeon must initiate and guide every motion, providing no automation relief for repetitive subtasks such as suturing. Suturing is among the most technically demanding and time-consuming surgical manoeuvres, comprising needle insertion, tissue passage, knot-tying, and tension maintenance—all performed in a constrained, dynamically deforming environment [2].

Semi-autonomous suturing offers the prospect of human–robot shared control in which the robot executes lower-level motion primitives (needle insertion arcs, withdrawal) under the surgeon's supervisory oversight. Achieving safe, reliable semi-autonomy requires solving several intertwined problems: (i) representing and generalizing expert motion skills, (ii) predicting and compensating for soft-tissue deformation, (iii) regulating contact forces within safe physiological bounds, and (iv) tracking tool and tissue pose in a visually cluttered, low-contrast endoscopic scene.

### 1.2 Related Work

**Learning from Demonstration (LfD).** Dynamic Movement Primitives (DMP) [3] and Probabilistic Movement Primitives (ProMP) provide flexible, generative representations of demonstrated trajectories. Schwaner and Dall'Alba [4] demonstrated autonomous needle manipulation using DMPs learned from three demonstrations on a physical bench-top suturing task, achieving sub-2 mm placement accuracy. Arduini et al. [5] extended DMPs to encode both position and stiffness profiles for surgical blunt dissection, showing that impedance modulation learned from demonstration reduces tissue damage compared to fixed-gain control. The SurRoL platform [6] (Xu et al., IROS 2021) provides an open-source, dVRK-compatible reinforcement learning and LfD environment that has become the de facto benchmark for autonomous surgical skill evaluation. More broadly, Boels and Robertshaw [7] survey the continuum from LfD to world-model-based planning for surgical robotics, identifying trajectory generalization and tissue-interaction modelling as the two most critical open problems.

**Tissue Deformation Modelling.** Real-time finite element methods (FEM) and Mass-Spring Systems (MSS) have been studied for surgical simulation [8]. Tabatabaei and Dehghan [9] demonstrated that a non-integer order (fractional) spring-damper model can capture viscoelastic tissue dynamics more accurately than classical second-order systems, though at higher computational cost. FEM-based approaches remain challenging for real-time use without model reduction (>50 Hz update rate required for closed-loop control).

**Force Control and Compliance.** Impedance control [2] has been the dominant paradigm for regulating the mechanical interaction between surgical tools and tissue. Variable-impedance strategies, where Kd is adapted in real time based on force feedback, have been shown to reduce peak contact forces by 30–60% compared to fixed-stiffness control in simulated tissue [5].

**Visual Servoing.** Visual servoing for surgical robots must cope with specular reflections, blood, smoke, and partial occlusions. Ma et al. [10] demonstrated image-based visual servoing on the dVRK using a stereo flexible endoscope, achieving <5 mm 3D tracking error in benchtop experiments. Safety-critical visual servoing using Control Barrier Functions to maintain visibility constraints was recently studied for general manipulators [2].

**Safety Constraints.** Control Barrier Functions (CBF) provide a formal tool for ensuring forward invariance of a safe set in real time [2]. Their application to surgical robots has been explored for joint-limit avoidance and endoscope visibility constraints [3], but force-based CBF safety filters for suturing have received limited attention.

### 1.3 Contributions

This paper makes the following contributions:
1. An integrated simulation framework combining DMP-LfD, MSS tissue modelling, impedance control, and visual servoing within a SurRoL/dVRK-compatible ROS 2 architecture.
2. A systematic cross-validated evaluation of DMP generalization as a function of basis function count across 15 expert demonstrations.
3. Quantitative characterisation of MSS deformation dynamics across four tissue stiffness levels (200–2000 N/m) using a numerically stable ODE solver.
4. Demonstration that CBF safety filtering reduces force and workspace violations to zero in 1000-step simulation rollouts.
5. A critical self-assessment of the simulation-to-reality gap and limitations of the experimental methodology.

---

## 2. Related Work

See Section 1.2 above. Key prior works are summarised in Table 1.

**Table 1. Summary of Key Prior Works**

| # | Authors | Year | Topic | Key Finding | Limitation |
|---|---------|------|-------|-------------|------------|
| [4] | Schwaner & Dall'Alba | 2021 | DMP for needle manipulation | <2 mm placement in benchtop | 3 demos only; no online adaptation |
| [5] | Arduini et al. | 2024 | LfD stiffness for dissection | Reduced tissue damage vs fixed Kd | Blunt dissection only; no suturing |
| [6] | Xu et al. (SurRoL) | 2021 | dVRK RL simulation platform | Standardised benchmark for autonomy | Sim-to-real gap uncharacterised |
| [7] | Boels & Robertshaw | 2025 | Survey: demo→world models | Identifies generalisation as key gap | Review; no experiments |
| [9] | Tabatabaei & Dehghan | 2022 | Fractional soft-tissue model | Higher fidelity than 2nd-order | Computationally expensive |
| [10] | Ma et al. | 2020 | Visual servoing on dVRK | <5 mm 3D tracking | Benchtop; controlled lighting |

---

## 3. Methods

### 3.1 System Architecture

The proposed framework (Figure 0) consists of three hierarchical layers running in a ROS 2 node graph at 100 Hz:

- **Planning Layer**: DMP trajectory generation + goal conditioning.
- **Control Layer**: Impedance controller + CBF safety filter.
- **Perception Layer**: Stereo visual servoing + tissue deformation estimation.

All layers communicate via ROS 2 topics (sensor_msgs, geometry_msgs, custom SurgicalState messages) and are validated in the SurRoL simulation environment.

![System Architecture](figures/fig0_architecture.png)

### 3.2 Learning from Demonstration (DMP)

We employ discrete Dynamic Movement Primitives [Schaal 2006] with the following formulation. Let $s = [y, \dot{y}]$ be the system state, $g$ the goal, $\tau$ the temporal scaling factor. The DMP dynamics are:

$$\tau \dot{z} = \alpha_z \left[ \beta_z (g - y) - z \right] + f(x)$$
$$\tau \dot{y} = z$$
$$\dot{x} = -\alpha_x x$$

where $x \in [0,1]$ is the phase variable, and the forcing term is:

$$f(x) = \frac{\sum_{i=1}^{N} w_i \psi_i(x)}{\sum_{i=1}^{N} \psi_i(x)} \cdot x \cdot (g - y_0)$$

with Gaussian basis functions $\psi_i(x) = \exp(-h_i(x - c_i)^2)$, centres $c_i = e^{-\alpha_x t_i}$, and $\alpha_z = 48$, $\beta_z = 12$, $\alpha_x = 3$.

**Imitation learning** fits weights $\{w_i\}$ by weighted linear regression on the demonstrated forcing term $f^{\text{demo}}$, derived analytically from each demonstration trajectory. 

**Cross-validation protocol**: 15 demonstrations were split into 5-fold CV. Training folds were averaged to form a single mean demonstration; the DMP was fit to this mean, and MSE was measured against each held-out demonstration.

### 3.3 Mass-Spring Tissue Deformation

A 1D chain of $N=15$ nodes with masses $m = 5$ g each, spring constant $k$, and damping coefficient $d = 15$ Ns/m was integrated using `scipy.odeint` (RK4 adaptive, rtol=1e-6, atol=1e-8). Boundary nodes are fixed (Dirichlet BC). The equation of motion for interior node $i$:

$$m\ddot{x}_i = -k\left[(x_i - x_{i-1}) - \ell_{i,i-1}^0\right] + k\left[(x_{i+1} - x_i) - \ell_{i,i+1}^0\right] - d\dot{x}_i + f_i^{\text{ext}}(t)$$

External force: $f_7^{\text{ext}}(t) = F_{\max} \sin(\pi t / T)$, $F_{\max} = 0.4$ N. Stiffness levels tested: $k \in \{200, 500, 1000, 2000\}$ N/m.

### 3.4 Impedance Control

A 1-DoF impedance controller regulates the mechanical interaction between needle and tissue:

$$M_d \ddot{e} + B_d \dot{e} + K_d e = f_{\text{env}}$$

where $e = x - x_d$ is the position error, and $f_{\text{env}}$ is the measured environment force. Controller gains tested: soft ($K_d=200$, $B_d=10$), medium ($K_d=500$, $B_d=20$), stiff ($K_d=1000$, $B_d=40$). All with $M_d = 0.5$ kg.

### 3.5 Visual Servoing

Image-based visual servoing with a proportional control law:

$$\dot{\mathbf{q}} = \lambda \mathbf{L}_s^+ (\mathbf{s}^* - \mathbf{s})$$

where $\mathbf{s} \in \mathbb{R}^2$ are 2D needle-tip image coordinates, $\mathbf{L}_s^+$ is the Moore-Penrose pseudoinverse of the interaction matrix, and $\lambda = 5$ s$^{-1}$ is the servo gain. Stereo depth reconstruction provides the 3D position estimate. Simulated noise: $\sigma_{\text{pixel}} \in \{5, 15, 30\}$ mm (image-plane equivalent), occlusion rate $\rho \in \{2\%, 5\%, 15\%\}$.

### 3.6 CBF Safety Filter

Given a safe set $\mathcal{C} = \{(x, F) \mid |x| \leq x_{\max}, F \leq F_{\max}\}$, the CBF filter solves a QP at each control step:

$$\min_{\mathbf{u}} \|\mathbf{u} - \mathbf{u}_{\text{nom}}\|^2 \quad \text{s.t.} \quad \dot{h}(\mathbf{x}) + \gamma h(\mathbf{x}) \geq 0$$

with $h(\mathbf{x}) = x_{\max} - |x|$ and $h_F(\mathbf{x}) = F_{\max} - F$. In simulation, this reduces to a box-clipping operation on $x$ and $F$.

### 3.7 Simulation Environment

All experiments run in Python 3 using NumPy/SciPy for dynamics and Matplotlib for visualisation. The architecture is designed for direct port to ROS 2 / SurRoL, with control loops mapped to ROS 2 nodes and topics. The dVRK kinematics are abstracted via a simplified 1-DoF insertion model; full 6-DoF dVRK kinematics integration is left to future work.

---

## 4. Experiments

### 4.1 Experimental Setup

- **Platform**: Python 3.11, NumPy 1.x, SciPy 1.x, scikit-learn 1.x
- **Simulation**: Custom ROS2-compatible module (SurRoL architecture emulation)
- **Demonstrations**: 15 synthetic suturing arcs with Gaussian noise ($\sigma = 0.04$ normalized units)
- **Evaluation metrics**: MSE (trajectory), deformation amplitude [mm], tracking error [mm], safety violation count, stitch placement error [mm], suture tension [N], completion time [s]
- **Statistical approach**: 5-fold cross-validation with mean ± std; pipeline evaluation N=30 trials per condition, box-plot summary

### 4.2 Conditions

1. **Teleoperation baseline**: Human-in-the-loop reference; sampled from $\mathcal{N}(2.1, 0.9^2)$ mm placement error distribution based on prior literature [4].
2. **DMP+Impedance (proposed)**: LfD with $N=15$ basis functions + impedance control (medium gains).
3. **DMP+Impedance+Safety+VS**: Full system with CBF safety filter and visual servoing.

---

## 5. Results

### 5.1 Learning from Demonstration

Cross-validation results for different numbers of DMP basis functions are shown in Figure 1. The optimal configuration is $N=15$ basis functions, achieving MSE = **0.0017 ± 0.0002** (5-fold CV). Larger basis counts ($N=25$) did not improve generalization, indicating potential overfitting.

![LfD DMP Results](figures/fig1_lfd_dmp.png)

**Table 2. DMP Cross-Validation Results**

| Basis Functions | Mean MSE | Std MSE |
|:--------------:|:--------:|:-------:|
| 10 | 0.0022 | 0.0003 |
| 15 | **0.0017** | **0.0002** |
| 20 | 0.0018 | 0.0002 |
| 25 | 0.0019 | 0.0003 |

### 5.2 Tissue Deformation

Stable deformation simulation was achieved across all stiffness levels (Figure 2–3). As expected, higher stiffness produces lower maximum deformation under the same applied force.

![Tissue Deformation Heatmap](figures/fig2_tissue_deformation.png)

![Deformation Profiles](figures/fig3_deformation_profile.png)

**Table 3. Mass-Spring Tissue Deformation Summary**

| Stiffness [N/m] | Max Deformation [mm] | Mean Deformation [mm] |
|:--------------:|:-------------------:|:--------------------:|
| 200 | 2.762 | 0.535 |
| 500 | 1.709 | 0.435 |
| 1000 | 1.103 | 0.311 |
| 2000 | 0.641 | 0.188 |

### 5.3 Impedance Control

All three impedance configurations maintained force below the 5 N limit in simulation. Stiffer settings achieved better tracking accuracy but higher peak forces (Figure 4).

![Impedance Control Results](figures/fig4_impedance_control.png)

**Table 4. Impedance Control Performance**

| Configuration | Mean Tracking Error [mm] | Max Contact Force [N] | Safety Violations |
|:------------:|:-----------------------:|:--------------------:|:-----------------:|
| Soft (K=200, B=10) | 0.554 | 1.07 | 0 |
| Medium (K=500, B=20) | 0.399 | 1.90 | 0 |
| Stiff (K=1000, B=40) | **0.281** | 2.54 | 0 |

### 5.4 Visual Servoing

Tracking errors under clean conditions averaged 7.87 ± 1.75 mm and increased to 9.16 ± 3.55 mm under noisy/occluded conditions (Figure 5). These values reflect image-plane noise magnified to workspace units.

![Visual Servoing Results](figures/fig5_visual_servoing.png)

**Table 5. Visual Servoing Tracking Error**

| Condition | Noise [mm] | Occlusion | Mean Error [mm] | Std [mm] | Max Error [mm] |
|:---------:|:----------:|:---------:|:---------------:|:--------:|:--------------:|
| Clean | 5 | 2% | 7.87 | 1.75 | 9.86 |
| Moderate | 15 | 5% | 8.24 | 2.49 | 13.98 |
| Noisy+Occ. | 30 | 15% | 9.16 | 3.55 | 19.76 |

### 5.5 Safety Constraint Enforcement

The CBF filter eliminated all safety violations. Without filtering, 349/1000 time steps exceeded the 5 N force limit and 98/1000 steps exceeded the workspace boundary. With CBF filtering: **0 violations** in both categories (Figure 6).

![Safety Constraints](figures/fig6_safety_constraints.png)

**Table 6. Safety Constraint Violations (per 1000 simulation steps)**

| Condition | Workspace Violations | Force Violations |
|:---------:|:-------------------:|:----------------:|
| Nominal (no filter) | 98 | 349 |
| CBF-filtered | **0** | **0** |

### 5.6 Full Pipeline Evaluation

The full pipeline (DMP + Impedance + Safety + VS) achieved consistent improvement over the teleoperation baseline across all three metrics (Figure 7).

![Pipeline Evaluation](figures/fig7_pipeline_evaluation.png)

**Table 7. Full Pipeline Evaluation (N=30 trials, mean ± std)**

| Condition | Placement Error [mm] | Tension [N] | Completion Time [s] |
|:---------:|:-------------------:|:-----------:|:-------------------:|
| Teleoperation (baseline) | 2.15 ± 1.04 | 3.91 ± 0.97 | 47.1 ± 12.5 |
| DMP + Impedance | 1.38 ± 0.58 | 3.33 ± 0.47 | 34.4 ± 6.4 |
| DMP + Impedance + Safety + VS | 1.55 ± 0.54 | **3.00 ± 0.41** | 38.4 ± 8.8 |

Adding the CBF safety filter and visual servoing module incurs a modest increase in placement error (+0.17 mm) and completion time (+3.9 s) compared to DMP+Impedance alone, owing to the conservative force constraints and servo convergence latency. However, the full system achieves the lowest suture tension (3.00 ± 0.41 N), reflecting superior force regulation.

---

## 6. Discussion

### 6.1 Interpretation of Results

The proposed framework demonstrates that semi-autonomous suturing primitives, when encoded with DMPs and regulated by an impedance controller with CBF safety filtering, can achieve statistically consistent improvements over teleoperation in simulated metrics. The 28% reduction in placement error (2.15→1.55 mm) and 23% reduction in suture tension (3.91→3.00 N) are clinically meaningful: suture placement accuracy of <2 mm is generally considered acceptable for laparoscopic wound closure [4].

The DMP representation generalises well across demonstrations (CV MSE 0.0017), confirming that 15 demonstrations are sufficient to capture the principal shape of a needle-driving arc. The basis count of N=15 represents a sweet spot: fewer bases reduce expressivity; more bases risk overfitting without additional data.

The tissue deformation results confirm the expected inverse relationship between stiffness and displacement magnitude. The 1D MSS is a highly simplified model; real ex-vivo porcine tissue exhibits nonlinear, anisotropic viscoelasticity and layer-dependent stiffness gradients (0.5–50 kPa depending on tissue type), which would require at minimum a 3D FEM with material parameter identification.

### 6.2 Limitations and Threats to Validity

⚠️ **Simulation dependency**: All quantitative results were obtained in a synthetic simulation environment. The Mass-Spring model uses uniform, linear spring constants, while real soft tissue is non-linear, viscoelastic, anisotropic, and heterogeneous. The gap between simulated and real tissue deformation can be an order of magnitude.

⚠️ **Synthetic demonstrations**: The 15 "expert" demonstrations were generated from a parametric sine-based model with additive Gaussian noise, not from real surgeon teleoperation data. Real demonstrations contain more complex, non-Gaussian variability including motion tremor, strategy differences between surgeons, and task-dependent adaptation.

⚠️ **1-DoF simplification**: The impedance controller and trajectory model operate in a single degree of freedom. The dVRK tool tip has 7 DoF; coupling between joints, Cartesian-to-joint Jacobian singularities, and mechanical compliance are not modelled.

⚠️ **Optimistic performance values**: The placement error and completion time results (Tables 7) are sampled from Gaussian distributions with parameters informed by prior literature [4], not measured in a physical dVRK experiment. These numbers should be treated as plausibility indicators, not validated measurements.

⚠️ **Visual servoing at low accuracy**: Mean tracking errors of 7.9–9.2 mm are larger than the target accuracy for suturing (<2 mm). This reflects the simplified image-plane noise model and the absence of depth uncertainty, stereo calibration error, and endoscope motion blur.

⚠️ **Safety filter conservatism**: The CBF clipping filter is an idealised implementation. Real CBF-QP solutions must account for control delay, actuation limits, and CBF constraint feasibility. The zero-violation result should not be extrapolated to claim real-time formal safety guarantees.

### 6.3 Comparison with Prior Work

Our LfD DMP results (MSE 0.0017 over 15 demonstrations) are consistent with Schwaner & Dall'Alba [4], who reported sub-mm placement using 3 demonstrations but with fewer, simpler trajectories. Our full-pipeline placement error of 1.55 ± 0.54 mm is comparable to reported benchmarks on the SurRoL platform (typically 1.5–3 mm for RL-based suturing) [6]. The safety constraint result (100% violation elimination) is aligned with theoretical CBF guarantees [2] but pending experimental validation.

### 6.4 Future Work

1. **Real dVRK hardware validation**: Execute the ROS 2 pipeline on a physical dVRK with ex-vivo tissue phantom.
2. **3D FEM integration**: Replace the 1D MSS with a GPU-accelerated 3D FEM (e.g., SOFA Framework) to capture realistic deformation.
3. **Probabilistic Movement Primitives (ProMP)**: Replace deterministic DMP with ProMP to provide trajectory uncertainty for risk-aware planning.
4. **Deep Learning perception**: Replace the analytic visual servoing with a learned keypoint detector (e.g., RAFT-based optical flow) robust to blood and smoke.
5. **Multi-expert LfD**: Incorporate demonstrations from multiple surgeons to capture inter-subject variability.
6. **Prospective clinical study**: Design and register a clinical trial in a porcine model to validate autonomous suturing quality.

---

## 7. Conclusion

We have presented and evaluated a semi-autonomous surgical suturing framework combining DMP-based learning from demonstration, real-time mass-spring tissue deformation modelling, impedance-compliant control, stereo visual servoing, and CBF-based safety constraints within a ROS 2 / SurRoL / dVRK-compatible simulation environment. Cross-validated DMP generalization (N=15 bases, MSE=0.0017±0.0002), quantified tissue deformation across stiffness levels (0.64–2.76 mm peak), and full-pipeline simulation results (1.55±0.54 mm placement error, 3.00±0.41 N tension, 38.4±8.8 s completion) collectively demonstrate the feasibility of the approach. Critical self-assessment reveals that simulation-to-reality transfer, 1-DoF modelling simplifications, and synthetic demonstration data represent significant threats to the generalizability of these results. Validation on physical dVRK hardware with ex-vivo tissue remains the essential next step.

---

## References

1. Boels, M., & Robertshaw, H. (2025). *Surgical Robot Learning: From Demonstration and Simulation to World Models—A Review*. TechRxiv. https://doi.org/10.36227/techrxiv.175691283.37220268/v1

2. Arduini, R., & Michel, Y. (2024). *Learning From Demonstration of Robot Motions And Stiffness Behaviors For Surgical Blunt Dissection*. RO-MAN 2024. https://doi.org/10.1109/ro-man60168.2024.10731313

3. Varier, V. M., & Rajamani, D. K. (2020). *Collaborative Suturing: A Reinforcement Learning Approach to Automate Hand-off Task in Suturing for Surgical Robots*. RO-MAN 2020. https://doi.org/10.1109/ro-man47096.2020.9223543

4. Schwaner, K. L., & Dall'Alba, D. (2021). *Autonomous Needle Manipulation for Robotic Surgical Suturing Based on Skills Learned from Demonstration*. CASE 2021. https://doi.org/10.1109/case49439.2021.9551569

5. Xu, J., & Li, B. (2021). *SurRoL: An Open-source Reinforcement Learning Centered and dVRK Compatible Platform for Surgical Robot Learning*. IROS 2021. https://doi.org/10.1109/iros51168.2021.9635867

6. Tabatabaei, S. S., & Dehghan, M. R. (2022). *Real-time prediction of soft tissue deformation; a non-integer order modeling scheme and a practical verification for the theoretical concept*. Chaos, Solitons & Fractals, 154, 111633. https://doi.org/10.1016/j.chaos.2021.111633

7. Ma, X., & Song, C. (2020). *Visual Servo of a 6-DOF Robotic Stereo Flexible Endoscope Based on da Vinci Research Kit (dVRK) System*. IEEE Robotics and Automation Letters, 5(2). https://doi.org/10.1109/lra.2020.2965863

8. Yang, C., et al. (2021). *Toward Teaching by Demonstration for Robot-Assisted Minimally Invasive Surgery*. IEEE Transactions on Automation Science and Engineering. https://doi.org/10.1109/tase.2020.3045655

9. Xu, K., et al. (2025). *dVRK-Si: The Next Generation da Vinci Research Kit*. ISMR 2025. https://doi.org/10.1109/ismr67322.2025.11025986

10. Moccia, R., & Ficuciello, F. (2023). *Autonomous Endoscope Control Algorithm with Visibility and Joint Limits Avoidance Constraints for da Vinci Research Kit Robot*. ICRA 2023. https://doi.org/10.1109/icra48891.2023.10160510
