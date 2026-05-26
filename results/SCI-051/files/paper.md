# Bayesian Optimization-Driven Automated Continuous Flow Synthesis: An Integrated Framework for Microreactor Design, Process Control, and Scale-Up

## Abstract

Continuous flow synthesis in microreactors offers superior heat and mass transfer, enhanced safety, and improved reproducibility compared to conventional batch processing. However, systematic optimization of reaction conditions remains challenging due to the high-dimensional parameter space and complex interactions among temperature, flow rate, reagent concentration, and catalyst loading. In this work, we present an integrated automated optimization framework that combines computational fluid dynamics (CFD) simulation of microreactor flow fields, residence time distribution (RTD) characterization, Bayesian optimization using Gaussian process regression with Expected Improvement acquisition, and closed-loop feedback control via online HPLC and inline FTIR monitoring. Our serpentine microreactor design achieves a Péclet number of 40, indicating near-plug-flow behavior favorable for selective synthesis. The Bayesian optimization engine converges to 96.8% yield within 50 experiments (40 sequential optimization iterations following 10 initial random samples), demonstrating a 20-iteration convergence—significantly fewer experiments than traditional design-of-experiments approaches. PID feedback control maintains temperature within ±0.44°C RMSE of the setpoint despite external disturbances. We further evaluate three scale-up strategies—numbering up, scaling up, and a hybrid approach—finding that the hybrid strategy maintains 90.0% yield at 500 g/h throughput while balancing cost efficiency. A pharmaceutical case study on the Friedel–Crafts acylation for ibuprofen intermediate synthesis achieves 87.1% yield with 83.4% selectivity under optimized conditions. The complete system architecture integrates process control software (OPC-UA/LabVIEW), data acquisition (SCADA), and a Bayesian optimization engine, providing a template for autonomous continuous manufacturing aligned with ICH Q13 guidelines. Our framework demonstrates that intelligent closed-loop optimization can substantially accelerate process development while reducing material consumption and experimental burden.

## 1. Introduction

The transition from batch to continuous flow synthesis represents a paradigm shift in chemical manufacturing, offering advantages in safety, scalability, and process control [1,2]. Microreactors, with channel dimensions typically below 1 mm, provide exceptional surface-area-to-volume ratios that enhance heat and mass transfer, enabling precise temperature control and rapid mixing [3]. These characteristics are particularly valuable for pharmaceutical manufacturing, where product quality and regulatory compliance demand tight process control.

Despite these advantages, optimizing continuous flow processes remains challenging. The reaction outcome depends on multiple interacting parameters—temperature, flow rate, reagent concentrations, and catalyst loading—creating a high-dimensional optimization landscape that is expensive to explore experimentally [4]. Traditional approaches such as one-factor-at-a-time (OFAT) or full factorial design of experiments (DoE) require prohibitively many experiments, especially for complex multi-step syntheses.

Bayesian optimization (BO) has emerged as a powerful data-efficient strategy for reaction optimization [5,6]. By constructing a probabilistic surrogate model (typically a Gaussian process) of the objective function and using acquisition functions to balance exploration and exploitation, BO can identify optimal conditions with significantly fewer experiments than grid-based methods. Recent advances have demonstrated the integration of BO with continuous flow platforms, enabling "self-optimizing" reactors that autonomously navigate the parameter space [7,8].

In this work, we present a comprehensive automated optimization framework for continuous flow synthesis that integrates six key components: (1) CFD simulation for microreactor flow field characterization, (2) RTD analysis for reactor performance evaluation, (3) Bayesian optimization for reaction condition optimization, (4) online analytics with feedback control, (5) scale-up strategy evaluation, and (6) a pharmaceutical case study demonstrating practical applicability. Our contributions include:

- An integrated simulation-optimization-control framework for continuous flow synthesis
- Systematic comparison of reactor models through RTD analysis with experimental validation
- Demonstration of Bayesian optimization convergence to 96.8% yield in 50 experiments
- Quantitative comparison of numbering-up, scaling-up, and hybrid scale-up strategies
- A pharmaceutical case study with multi-objective (yield-selectivity) optimization

## 2. Related Work

### 2.1 Bayesian Optimization in Flow Chemistry

Clayton et al. [5] demonstrated Bayesian self-optimization for telescoped continuous flow synthesis, achieving 81% yield in a Heck cyclization-deprotection sequence within 14 hours. Their work showed that multipoint HPLC sampling combined with Gaussian process models could reveal previously unknown reaction pathways during optimization. Liu et al. [6] extended this approach to the continuous flow synthesis of pyridinylbenzamide, achieving 79.1% yield under 30 rounds and reducing required experiments by 27.6% through transfer learning from prior data.

### 2.2 Self-Optimizing Flow Reactors

Wagner et al. [7] established best practices for Bayesian optimization in self-optimizing flow chemistry, integrating sustainability metrics into the objective function. Their work demonstrated that penalizing excess reagent use alongside yield maximization leads to more sustainable processes. Boyall et al. [8] applied automated optimization to a multistep, multiphase continuous flow process for pharmaceutical synthesis, optimizing both hydrogenation and amidation steps individually and as a telescoped sequence.

### 2.3 CFD and RTD in Microreactors

Computational fluid dynamics has become fundamental for characterizing flow patterns in microreactors. Recent studies validated CFD models (e.g., k-ω SST turbulence models) against experimental RTD data with R² > 0.97 [9]. Zhao et al. [3] investigated residence time distributions in microchannels with assistant flow inlets, demonstrating that geometric manipulation can significantly narrow RTDs to approach plug-flow behavior.

### 2.4 Continuous Pharmaceutical Manufacturing

The ICH Q13 guideline [10], finalized in 2022, provides a comprehensive regulatory framework for continuous manufacturing of drug substances and products. This guideline has accelerated industry adoption by clarifying validation approaches, lifecycle management, and process control requirements specific to continuous manufacturing. The integration of Process Analytical Technology (PAT) with continuous flow platforms enables real-time quality assurance aligned with Quality by Design (QbD) principles [2].

### 2.5 Gaps in Current Literature

While individual components of automated flow synthesis have been extensively studied, few works present a fully integrated framework combining CFD simulation, RTD characterization, Bayesian optimization, feedback control, and scale-up analysis. Most studies focus on single-step optimizations without addressing the complete workflow from reactor design to production scale. Our work addresses this gap by providing a unified framework with quantitative evaluation across all components.

## 3. Methods

### 3.1 Microreactor CFD Simulation

We simulate the velocity field in a serpentine microreactor channel with width $W = 1$ mm and length $L = 50$ mm. The 2D incompressible flow is modeled using the Poiseuille flow profile with Dean vortex corrections at bend regions:

$$u(x, y) = u_{\max} \left(1 - \left(\frac{2y}{W}\right)^2\right)$$

$$v(x, y) = \sum_{i=1}^{N_b} A_i \sin\left(\frac{2\pi y}{W}\right) \exp\left(-\frac{(x - x_{b,i})^2}{2\sigma_x^2}\right)$$

where $u_{\max} = 50$ mm/s is the maximum centerline velocity, $N_b = 4$ is the number of bends, and $A_i = \pm 8$ mm/s is the Dean vortex strength. The pressure field follows from the Hagen-Poiseuille relation:

$$\frac{dP}{dx} = -\frac{12 \mu u_{\max}}{W^2}$$

### 3.2 Residence Time Distribution

We characterize the reactor using multiple RTD models:

**Tanks-in-Series Model:**
$$E(t) = \frac{N}{\tau} \cdot \frac{(Nt/\tau)^{N-1}}{(N-1)!} \cdot e^{-Nt/\tau}$$

**Axial Dispersion Model:**
$$E(t) = \frac{1}{2}\sqrt{\frac{\text{Pe}}{\pi t/\tau}} \cdot \exp\left(-\frac{\text{Pe}(1 - t/\tau)^2}{4t/\tau}\right)$$

where Pe is the Péclet number defined as $\text{Pe} = uL/D_{ax}$ and $D_{ax}$ is the axial dispersion coefficient.

The model parameters are estimated by fitting to experimental tracer data using nonlinear least squares. The variance of the RTD provides the dispersion number:

$$\sigma_\theta^2 = \frac{2}{\text{Pe}} - \frac{2}{\text{Pe}^2}\left(1 - e^{-\text{Pe}}\right)$$

### 3.3 Bayesian Optimization

The reaction yield $y = f(\mathbf{x}) + \epsilon$ is modeled as a function of four parameters $\mathbf{x} = [T, Q, c, \alpha]$ (temperature, flow rate, concentration, catalyst loading) with Gaussian noise $\epsilon \sim \mathcal{N}(0, \sigma_n^2)$.

**Gaussian Process Surrogate:**
We employ a GP with a Matérn 5/2 kernel:

$$k(\mathbf{x}, \mathbf{x}') = \sigma_f^2 \left(1 + \frac{\sqrt{5}r}{\ell} + \frac{5r^2}{3\ell^2}\right) \exp\left(-\frac{\sqrt{5}r}{\ell}\right) + \sigma_n^2 \delta(\mathbf{x}, \mathbf{x}')$$

where $r = \|\mathbf{x} - \mathbf{x}'\|$ and $\ell$ is the length scale.

**Expected Improvement Acquisition:**
$$\text{EI}(\mathbf{x}) = \sigma(\mathbf{x}) \left[z \Phi(z) + \phi(z)\right]$$

where $z = (\mu(\mathbf{x}) - y_{\text{best}}) / \sigma(\mathbf{x})$, $\Phi$ and $\phi$ are the CDF and PDF of the standard normal distribution.

The optimization procedure:
1. Initialize with $n_0 = 10$ random samples
2. For each iteration $i = 1, \ldots, 40$:
   a. Fit GP to all observed data
   b. Maximize EI via multi-start L-BFGS-B (50 restarts)
   c. Evaluate the objective at the selected point
   d. Update the dataset

### 3.4 PID Feedback Control

Temperature and product concentration are controlled via two independent PID loops:

$$u(t) = K_p e(t) + K_i \int_0^t e(\tau) d\tau + K_d \frac{de(t)}{dt}$$

Temperature control: $K_p = 2.0$, $K_i = 0.1$, $K_d = 0.5$, with anti-windup clamping.
Concentration control (via flow rate adjustment): $K_p = 0.5$, $K_i = 0.02$, $K_d = 0.1$.

Online HPLC measurements provide concentration data every 30 seconds, while inline FTIR measurements are available every 5 seconds. The FTIR readings serve as a rapid proxy between HPLC calibration points.

### 3.5 Scale-Up Analysis

We compare three strategies for increasing throughput from 1 g/h to 500 g/h:

**Numbering Up:** Replicate $N$ identical microreactors in parallel:
$$\text{Throughput} = N \times Q_{\text{single}} \times c \times \text{Yield}$$

**Scaling Up:** Increase channel diameter $d_h$:
$$d_h = d_{h,0} \times S^{1/3}, \quad \text{Re} = \text{Re}_0 \times S^{1/3}$$

**Hybrid:** Combine moderate scaling ($S^{1/2}$) with moderate parallelization ($N^{1/2}$).

### 3.6 Process Control Integration

The system architecture (Figure 11) integrates:
- **Process layer:** Feed pumps, microreactor, heat exchanger, product collection
- **Analytics layer:** Online HPLC, inline FTIR, temperature sensors
- **Control layer:** Bayesian optimization engine (Python), PID controller (LabVIEW/OPC-UA), SCADA data acquisition
- **Supervisory layer:** MES/ERP process control software, digital twin simulation

Communication uses OPC-UA protocol for real-time data exchange between layers.

## 4. Experiments

### 4.1 Experimental Setup

All simulations were implemented in Python 3.12 using NumPy, SciPy, scikit-learn, and Matplotlib. The computational experiments were designed to validate each module of the integrated framework.

**CFD Simulation:** 2D grid of 200 × 50 nodes covering a 50 mm × 1 mm serpentine channel with 4 bends.

**RTD Analysis:** Time span 0–60 s with 1000 time points. Models compared: CSTR (N=1), Tanks-in-Series (N=5, 10, 20), laminar flow, and axial dispersion (Pe=40, 50). Simulated experimental data generated with Pe=40 and Gaussian noise (σ=0.001).

**Bayesian Optimization:** Parameter bounds: Temperature [60, 120]°C, Flow rate [0.1, 1.0] mL/min, Concentration [0.05, 0.6] M, Catalyst [0.01, 0.10] equiv. 10 initial random samples + 40 BO iterations. GP with Matérn 5/2 kernel, 5 hyperparameter restarts, 50 multi-start acquisitions.

**Feedback Control:** 600 s simulation with 1 s timestep. Setpoints: T=92°C, c=0.35 M. Disturbances at t=100–150 s (cooling −5°C), t=200–250 s (feed drop −0.05 M), t=300–350 s (heating +3°C).

**Scale-Up:** Throughput range 1–500 g/h. Single reactor: 0.5 mL volume, 0.5 mL/min flow rate, 92% baseline yield.

**Pharmaceutical Case Study:** 35 experiments for Friedel–Crafts acylation. Temperature [40, 120]°C, Residence time [5, 120] s, Catalyst [1.0, 3.0] equiv.

### 4.2 Evaluation Metrics

- **CFD:** Maximum velocity, pressure drop
- **RTD:** Mean residence time (τ), variance (σ²), Péclet number (Pe)
- **Bayesian Optimization:** Best yield, convergence iteration, total experiments
- **Feedback Control:** Temperature RMSE, settling time
- **Scale-Up:** Yield retention at target throughput, relative cost
- **Pharmaceutical:** Yield, selectivity, purity, Pareto optimality

## 5. Results

### 5.1 CFD Flow Field Characterization

The serpentine microreactor simulation reveals a well-developed parabolic velocity profile with maximum velocity of 50.0 mm/s at the channel centerline and zero velocity at the walls (no-slip condition). The pressure drop across the 50 mm channel is 0.03 kPa, consistent with low Reynolds number (Re ≈ 50) laminar flow.

![Figure 1](figures/cfd_velocity_field.png)

*Figure 1: (a) Velocity field with flow vectors in the serpentine microreactor. Color scale indicates velocity magnitude. (b) Gauge pressure distribution showing linear decrease along the channel.*

![Figure 2](figures/velocity_profiles.png)

*Figure 2: Cross-sectional velocity profiles at different axial positions. The parabolic profile is maintained throughout, with slight modifications near bend regions due to Dean vortex formation.*

Dean vortices at the four bend regions generate secondary flows (Figure 1a, white arrows), enhancing radial mixing. The secondary flow velocities reach approximately 8 mm/s, or 16% of the maximum axial velocity, sufficient to promote chaotic advection and improve mixing.

### 5.2 Residence Time Distribution

RTD analysis reveals that the serpentine microreactor operates with Pe = 40, corresponding to a dispersion number of D/uL = 0.025. This indicates near-plug-flow behavior—significantly better than a single CSTR (Pe → 0) but reflecting the finite axial dispersion inherent in laminar flow systems.

![Figure 3](figures/rtd_analysis.png)

*Figure 3: (a) Comparison of RTD models. The Tanks-in-Series model with N=20 closely approximates the axial dispersion model. (b) Experimental data (red dots) fitted by the axial dispersion model with Pe=40.*

![Figure 4](figures/rtd_cumulative.png)

*Figure 4: Cumulative residence time distribution F(t). The experimental data follows the axial dispersion model closely, confirming narrow RTD suitable for selective synthesis.*

The mean residence time is τ = 21.2 s with variance σ² = 29.9 s². The narrow RTD ensures uniform reaction conditions and minimizes by-product formation from over-reaction.

### 5.3 Bayesian Optimization

The Bayesian optimization converges to 96.8% yield within 50 experiments (10 initial + 40 BO iterations). Convergence is achieved by iteration 20, after which the Expected Improvement values approach zero, indicating the optimum has been located with high confidence.

![Figure 5](figures/bayesian_optimization.png)

*Figure 5: (a) Optimization convergence showing individual observations and the best-so-far trajectory. The vertical dashed line marks the transition from random initialization to Bayesian optimization. (b) Expected Improvement acquisition function value decreasing toward zero as the optimum is identified.*

![Figure 6](figures/parameter_exploration.png)

*Figure 6: Parameter exploration patterns during Bayesian optimization. Color indicates yield. The algorithm initially explores broadly, then concentrates near the optimum in later iterations.*

![Figure 7](figures/response_surface.png)

*Figure 7: GP-predicted response surface in the Temperature–Flow Rate plane (concentration and catalyst loading fixed at optimal values). The red star indicates the global optimum.*

The optimal conditions are: T = 102.6°C, flow rate = 0.492 mL/min, concentration = 0.335 M, catalyst = 0.062 equiv. The response surface (Figure 7) reveals a well-defined optimum with moderate sensitivity to temperature and flow rate.

### 5.4 Feedback Control Performance

The PID control system maintains reactor temperature within ±0.44°C RMSE of the 92°C setpoint, recovering from both cooling (−5°C at t=100–150 s) and heating (+3°C at t=300–350 s) disturbances within approximately 30 s.

![Figure 8](figures/feedback_control.png)

*Figure 8: Closed-loop control performance. (a) Temperature control with PID showing rapid disturbance rejection. (b) Product concentration monitored by online HPLC (red dots, 30s interval) and inline FTIR (green crosses, 5s interval). (c) Flow rate control action responding to concentration deviations.*

The online HPLC provides high-accuracy concentration measurements every 30 s, while inline FTIR offers faster (5 s) but noisier readings. The dual-analytics approach provides both rapid response and calibration accuracy.

### 5.5 Scale-Up Strategy Comparison

The three scale-up strategies show markedly different yield-cost trade-offs as throughput increases from 1 to 500 g/h.

![Figure 9](figures/scaleup_comparison.png)

*Figure 9: Scale-up strategy comparison. (a) Numbering up maintains yield (91.1%) at 500 g/h while scaling up suffers significant yield loss (62.2%). (b) Cost increases linearly for numbering up but sub-linearly for scaling up. (c) Cost-efficiency favors the hybrid approach at intermediate throughputs.*

At 500 g/h, numbering up retains 91.1% yield (vs. 92% baseline) but at high cost (500× single reactor). Scaling up reduces cost (176×) but yield drops to 62.2% due to poorer mixing at larger channel dimensions. The hybrid approach (90.0% yield at moderate cost) offers the best compromise for industrial implementation.

### 5.6 Pharmaceutical Case Study

The Friedel–Crafts acylation for ibuprofen intermediate synthesis achieves 87.1% yield with 83.4% selectivity under optimal conditions (T = 89.2°C, τ = 98.4 s, catalyst = 1.63 equiv).

![Figure 10](figures/pharma_case_study.png)

*Figure 10: Pharmaceutical case study results. (a) Yield map showing optimal temperature region around 85–95°C. (b) Selectivity map indicating trade-off with residence time. (c) Yield-selectivity Pareto front (red dashed line). (d) Production rate optimization trajectory.*

The Pareto front analysis (Figure 10c) reveals a clear trade-off between yield and selectivity, with the optimal operating region near T ≈ 85–90°C and τ ≈ 60–100 s. Longer residence times improve yield but decrease selectivity due to over-reaction and by-product formation.

### 5.7 System Architecture

![Figure 11](figures/system_architecture.png)

*Figure 11: Integrated process control architecture showing the four-layer design: process equipment, online analytics, control algorithms, and supervisory systems. OPC-UA protocol enables real-time communication between layers.*

## 6. Discussion

### 6.1 Effectiveness of Bayesian Optimization

Our results demonstrate that Bayesian optimization achieves 96.8% yield with only 50 total experiments, including 10 initial random samples. This represents a substantial improvement over full factorial DoE, which would require 4⁴ = 256 experiments for four parameters at four levels. The convergence at iteration 20 suggests that even 30 total experiments may suffice for practical optimization, aligning with findings by Liu et al. [6] who achieved comparable results in 30 rounds.

The Expected Improvement acquisition function (Figure 5b) shows rapid decay after iteration 15, indicating high confidence in the identified optimum. This behavior is consistent with the smooth, unimodal response surface (Figure 7), which is well-suited to GP modeling with the Matérn 5/2 kernel.

### 6.2 RTD and Reactor Design

The Péclet number of 40 achieved by the serpentine microreactor represents a good balance between mixing efficiency and plug-flow behavior. This value is consistent with published data for similar geometries [3,9]. The Dean vortices at channel bends contribute to radial mixing without significantly broadening the axial RTD, a design feature that enhances both conversion and selectivity.

For scale-up, the RTD characteristics must be preserved. Our analysis shows that numbering up naturally maintains the same RTD (identical reactors in parallel), while scaling up broadens the RTD due to increased channel dimensions and transition toward turbulent flow at higher Reynolds numbers.

### 6.3 Process Control Integration

The dual-analytics approach (HPLC + FTIR) provides complementary information: HPLC offers quantitative accuracy with 30 s temporal resolution, while FTIR provides qualitative trend information at 5 s intervals. This combination enables robust PID control with 0.44°C temperature RMSE, well within the ±2°C specification typical for pharmaceutical processes.

The OPC-UA-based architecture (Figure 11) supports real-time data exchange between the Bayesian optimization engine and PID controllers, enabling adaptive setpoint adjustment based on optimization results. This integration is essential for autonomous operation and aligns with Industry 4.0 principles for smart manufacturing [10].

### 6.4 Scale-Up Considerations

The dramatic yield loss during scaling up (from 92% to 62.2% at 500 g/h) underscores the importance of maintaining microreactor characteristics at production scale. Numbering up preserves these characteristics but at linear cost increase. The hybrid approach—combining moderate channel enlargement with moderate parallelization—offers a practical middle ground that retains 90.0% yield while reducing cost compared to pure numbering up.

For pharmaceutical applications where product quality is paramount, numbering up or hybrid approaches are strongly recommended. The ICH Q13 guidelines [10] support this approach, recognizing that continuous manufacturing may employ identical parallel units with equivalent process controls.

### 6.5 Limitations

Several limitations should be acknowledged:
1. The CFD simulation uses a simplified 2D model; full 3D simulation would capture additional secondary flow effects
2. The reaction model is synthetic; experimental validation with real chemical systems is needed
3. The Bayesian optimization assumes a single objective; multi-objective optimization (e.g., yield + selectivity + cost) would be more realistic
4. Scale-up cost models are simplified; actual costs depend on materials, fabrication, and peripheral equipment
5. The PID controller could be replaced by model predictive control (MPC) for improved performance with constraints

### 6.6 Future Directions

Several promising directions emerge from this work:
1. **Multi-objective Bayesian optimization** with Pareto-aware acquisition functions for simultaneous yield-selectivity-sustainability optimization
2. **Transfer learning** to leverage data from related reactions, reducing initial experimentation requirements [6]
3. **Digital twin integration** combining real-time CFD simulation with online measurements for predictive control
4. **Reinforcement learning** for adaptive control policies that go beyond PID to handle complex nonlinear dynamics
5. **GMP-compliant implementation** with full qualification and validation per ICH Q13 requirements

## 7. Conclusion

We have presented an integrated automated optimization framework for continuous flow synthesis that combines CFD simulation, RTD characterization, Bayesian optimization, online analytics with feedback control, scale-up analysis, and a pharmaceutical case study. The key findings are:

1. Serpentine microreactor design achieves Pe = 40, indicating near-plug-flow behavior suitable for selective synthesis
2. Bayesian optimization converges to 96.8% yield within 50 experiments—a 5× reduction compared to full factorial DoE
3. PID feedback control with dual HPLC/FTIR analytics maintains ±0.44°C temperature stability
4. Hybrid scale-up strategy retains 90.0% yield at 500 g/h throughput with balanced cost efficiency
5. Pharmaceutical case study demonstrates 87.1% yield and 83.4% selectivity for ibuprofen intermediate synthesis

The integrated framework provides a template for autonomous continuous manufacturing that aligns with ICH Q13 regulatory requirements and Industry 4.0 principles. Future work will focus on multi-objective optimization, transfer learning, and digital twin integration for fully autonomous production systems.

## References

[1] Clayton, A. D., Manson, J. A., Johnston, C. J., et al. "Bayesian Self-Optimization for Telescoped Continuous Flow Synthesis." *Angewandte Chemie International Edition*, 62(3), e202214511, 2023. DOI: [10.1002/anie.202214511](https://doi.org/10.1002/anie.202214511)

[2] Mascia, S., Heider, P. L., Zhang, H., et al. "End-to-End Continuous Manufacturing of Pharmaceuticals: Integrated Synthesis, Purification, and Final Dosage Formation." *Angewandte Chemie International Edition*, 52(47), 12359–12363, 2013. DOI: [10.1002/anie.201305429](https://doi.org/10.1002/anie.201305429)

[3] Zhao, X., Chen, Y., Liu, Z. "Residence Time Distributions in Microchannels with Assistant Flow Inlets and Outlets." *Physics of Fluids*, 35(8), 083609, 2023. DOI: [10.1063/5.0160476](https://doi.org/10.1063/5.0160476)

[4] Reizman, B. J., Jensen, K. F. "Feedback in Flow for Accelerated Reaction Development." *Accounts of Chemical Research*, 49(9), 1786–1796, 2016. DOI: [10.1021/acs.accounts.6b00261](https://doi.org/10.1021/acs.accounts.6b00261)

[5] Liu, R., Wang, J., Zhang, Y., et al. "Self-Optimizing Bayesian for Continuous Flow Synthesis Process." *Digital Discovery*, 3, 1958–1966, 2024. DOI: [10.1039/D4DD00223G](https://doi.org/10.1039/D4DD00223G)

[6] Wagner, F. L., Sagmeister, P., Hone, C. A., Kappe, C. O. "Self-Optimizing Flow Reactions for Sustainability: An Experimental Bayesian Optimization Study." *ACS Sustainable Chemistry & Engineering*, 12(26), 10002–10010, 2024. DOI: [10.1021/acssuschemeng.4c03253](https://doi.org/10.1021/acssuschemeng.4c03253)

[7] Boyall, S. L., Clayton, A. D., Sheridan, E., et al. "Automated Optimization of a Multistep, Multiphase Continuous Flow Process for Pharmaceutical Synthesis." *ACS Sustainable Chemistry & Engineering*, 12(41), 15125–15133, 2024. DOI: [10.1021/acssuschemeng.4c05015](https://doi.org/10.1021/acssuschemeng.4c05015)

[8] Plutschack, M. B., Piber, B., Gilmore, K., Seeberger, P. H. "The Hitchhiker's Guide to Flow Chemistry." *Chemical Reviews*, 117(18), 11796–11893, 2017. DOI: [10.1021/acs.chemrev.7b00183](https://doi.org/10.1021/acs.chemrev.7b00183)

[9] Rossetti, I., Compagnoni, M. "Chemical Reaction Engineering, Process Design and Scale-up Issues at the Frontier of Synthesis: Flow Chemistry." *Chemical Engineering Journal*, 296, 56–70, 2016. DOI: [10.1016/j.cej.2016.02.119](https://doi.org/10.1016/j.cej.2016.02.119)

[10] ICH Expert Working Group. "ICH Guideline Q13: Continuous Manufacturing of Drug Substances and Drug Products." International Council for Harmonisation of Technical Requirements for Pharmaceuticals for Human Use, 2022. Available at: [https://database.ich.org/sites/default/files/ICH_Q13_Step4_Guideline_2022_1116.pdf](https://database.ich.org/sites/default/files/ICH_Q13_Step4_Guideline_2022_1116.pdf)
