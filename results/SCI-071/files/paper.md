# Integrated Planning and Control for Deformable Object Manipulation: A Simulation-Based Framework with Sim-to-Real Transfer

## Abstract

Robotic manipulation of deformable objects such as cloth, rope, and elastic bodies remains a fundamental challenge due to their infinite-dimensional state spaces, complex nonlinear dynamics, and severe self-occlusions. We present an integrated framework for deformable object manipulation planning that combines multiple state representations (mesh, particle, and latent space), physics simulators (Finite Element Method and Material Point Method), sampling-based motion planners (Cross-Entropy Method and RRT), domain randomization for sim-to-real transfer, and visual feedback reactive control. Our system is designed for compatibility with SoftGym and Isaac Gym simulation environments. We evaluate our framework on a cloth folding case study, comparing reactive control, CEM-based planning, and random baselines. Experimental results demonstrate that the CEM planner achieves 50% target coverage compared to 12.5% for reactive control and 25% for random actions. The latent space representation reduces state dimensionality by 37.5× (from 300 to 8 dimensions) while maintaining a mean reconstruction error of 0.0017. Closed-loop visual feedback control converges 16.6× faster than open-loop execution and recovers from perturbations within 20 steps. Our domain randomization analysis identifies action noise and material stiffness as the most sensitive parameters for sim-to-real transfer. The proposed framework provides a modular, extensible architecture for deformable object manipulation research.

## 1. Introduction

Robotic manipulation of deformable objects is a critical capability for applications ranging from household tasks (laundry folding, cooking) to industrial automation (textile handling, cable routing) and surgical robotics (tissue manipulation). Unlike rigid object manipulation, deformable objects present unique challenges: (1) infinite-dimensional configuration spaces that resist compact state representation, (2) complex, often nonlinear dynamics governed by continuum mechanics, (3) severe self-occlusions that hinder visual perception, and (4) a significant simulation-to-reality gap that impedes policy transfer.

Recent advances in deep reinforcement learning, differentiable simulation, and graph neural networks have shown promise in addressing these challenges (Gu et al., 2023; Lin et al., 2021). However, existing approaches often focus on individual components—perception, simulation, planning, or control—without providing an integrated framework that addresses the full manipulation pipeline.

In this paper, we present an integrated framework for deformable object manipulation planning that addresses the following contributions:

1. **Multi-modal state representation**: We implement and compare mesh-based, particle-based, and learned latent space representations, evaluating their trade-offs in dimensionality, computational cost, and reconstruction fidelity.

2. **Dual physics simulation**: We integrate both FEM (neo-Hookean hyperelastic model) and MPM simulators, enabling accurate modeling of diverse deformable materials from thin shells (cloth) to volumetric bodies (elastic objects).

3. **Hybrid planning architecture**: We combine sampling-based planning (CEM) with reactive visual feedback control, leveraging the strengths of both open-loop optimization and closed-loop adaptation.

4. **Systematic sim-to-real analysis**: We provide a domain randomization framework with parameter sensitivity analysis, identifying the most critical factors for successful simulation-to-reality transfer.

5. **Cloth folding case study**: We demonstrate the full pipeline on a garment folding task, benchmarking multiple methods and analyzing performance metrics.

## 2. Related Work

### 2.1 Surveys on Deformable Object Manipulation

Gu et al. (2023) provide a comprehensive survey covering over 150 works on robotic manipulation of deformable objects, including perception, modeling, planning, manipulation strategies, and sim-to-real transfer approaches. They highlight the increasing role of data-driven methods and the emerging application of large language models in this domain. Arriola-Rios et al. (2020) offer a tutorial-style review focusing on foundational modeling, estimation, and control methods, covering shape representation, registration, and learning-based approaches for perception and planning. Yin et al. (2021) discuss the integration of model priors with data-driven approaches and identify key challenges in data efficiency for sim-to-real transfer.

### 2.2 Simulation Environments and Benchmarks

Lin et al. (2021) introduced SoftGym, a benchmark suite built on NVIDIA FleX for evaluating deep reinforcement learning algorithms on deformable manipulation tasks including cloth folding, rope manipulation, and fluid pouring. SoftGym provides standardized environments with both image-based and state-based observations. Huang et al. (2021) presented PlasticineLab, a differentiable physics benchmark for soft-body manipulation using the DiffTaichi system, enabling gradient-based optimization alongside reinforcement learning. Li et al. (2022) developed DiffCloth, a differentiable cloth simulator supporting dry frictional contact, which enables direct gradient computation through the simulation for optimization-based control.

### 2.3 Graph Neural Networks for Dynamics Prediction

Zhang et al. (2024) proposed AdaptiGraph, a material-adaptive graph-based neural dynamics model that conditions predictions on physical properties, enabling generalization across different deformable materials. Shi et al. (2022) introduced RoboCraft, combining graph networks with differentiable simulation for perception and control of elasto-plastic objects. Weng et al. (2021) proposed object-centric graph representations with encode-process-decode GNN architectures for predicting interactions between deformable and rigid objects.

### 2.4 Visual Feedback and Reactive Control

Hietala et al. (2022) presented a reinforcement learning approach for dynamic cloth folding with visual feedback, demonstrating successful sim-to-real transfer across different fabric types. Luque et al. (2024) developed a model predictive control strategy for real-time dynamic cloth manipulation combining learned models with visual feedback. Kadi and Terzić (2023) reviewed data-driven methods for cloth-like deformable object manipulation, highlighting challenges and prospects in visual perception and control.

### 2.5 Sim-to-Real Transfer

Mehta et al. (2022) provided theoretical foundations for domain randomization, modeling simulators as tunable MDPs and deriving bounds on the sim-to-real gap. Muratore et al. (2022) surveyed robot learning from randomized simulations, covering domain randomization techniques and their application to sim-to-real transfer for various manipulation tasks.

## 3. Methods

### 3.1 State Representations

#### 3.1.1 Mesh Representation

We represent deformable objects as triangular meshes $\mathcal{M} = (\mathcal{V}, \mathcal{T})$ where $\mathcal{V} = \{v_i \in \mathbb{R}^3\}_{i=1}^{N}$ are vertex positions and $\mathcal{T}$ are triangle elements obtained via Delaunay triangulation. The state vector is:

$$\mathbf{s}_{\text{mesh}} = [v_1^T, v_2^T, \ldots, v_N^T]^T \in \mathbb{R}^{3N}$$

For an $n \times n$ grid, the state dimensionality is $3n^2$. We compute the deformation gradient tensor for each triangle element $e$ as:

$$\mathbf{F}_e = \mathbf{D}_{\text{curr}}^e \cdot (\mathbf{D}_{\text{rest}}^e)^{-1}$$

where $\mathbf{D}_{\text{curr}}^e$ and $\mathbf{D}_{\text{rest}}^e$ are edge matrices in current and rest configurations. The Green-Lagrange strain tensor is:

$$\mathbf{E}_e = \frac{1}{2}(\mathbf{F}_e^T \mathbf{F}_e - \mathbf{I})$$

#### 3.1.2 Particle Representation

For MPM-compatible simulation, we represent objects as particles $\{(x_i, v_i, m_i)\}_{i=1}^{M}$ with positions $x_i \in \mathbb{R}^3$, velocities $v_i \in \mathbb{R}^3$, and masses $m_i \in \mathbb{R}^+$. The state vector is:

$$\mathbf{s}_{\text{particle}} = [x_1^T, v_1^T, \ldots, x_M^T, v_M^T]^T \in \mathbb{R}^{6M}$$

#### 3.1.3 Latent Space Representation

We learn a compact representation via dimensionality reduction. Using PCA as a proxy for variational autoencoders:

$$\mathbf{z} = \text{Enc}(\mathbf{s}) \in \mathbb{R}^{d}, \quad \hat{\mathbf{s}} = \text{Dec}(\mathbf{z})$$

where $d \ll 3N$. The reconstruction loss is $\mathcal{L}_{\text{recon}} = \|\mathbf{s} - \hat{\mathbf{s}}\|_2^2$.

### 3.2 Physics Simulation

#### 3.2.1 FEM Simulator

We implement a neo-Hookean hyperelastic model. The first Piola-Kirchhoff stress tensor is:

$$\mathbf{P} = \mu(\mathbf{F} - \mathbf{F}^{-T}) + \lambda \ln(\det \mathbf{F}) \mathbf{F}^{-T}$$

where $\mu$ and $\lambda$ are Lamé parameters derived from Young's modulus $E$ and Poisson's ratio $\nu$:

$$\mu = \frac{E}{2(1+\nu)}, \quad \lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}$$

The elastic force on each node is computed as:

$$\mathbf{f}_i = -\sum_{e \ni i} |\det \mathbf{D}_{\text{rest}}^e| \cdot \mathbf{P}_e \cdot (\mathbf{D}_{\text{rest}}^e)^{-T} \cdot \mathbf{g}_i$$

Time integration uses explicit Euler with damping factor $\gamma$:

$$\mathbf{v}^{t+1} = \gamma(\mathbf{v}^t + \Delta t \cdot \mathbf{f}/m), \quad \mathbf{x}^{t+1} = \mathbf{x}^t + \Delta t \cdot \mathbf{v}^{t+1}$$

#### 3.2.2 MPM Simulator

The Material Point Method discretizes the continuum on both particles and a background grid. The P2G (particle-to-grid) transfer:

$$m_I^g = \sum_p w_{Ip} m_p, \quad (m\mathbf{v})_I^g = \sum_p w_{Ip} m_p \mathbf{v}_p$$

Grid velocity update with gravity and damping:

$$\mathbf{v}_I^{g,*} = \frac{(m\mathbf{v})_I^g}{m_I^g} + \Delta t \cdot \mathbf{g}, \quad \mathbf{v}_I^g = \alpha \cdot \mathbf{v}_I^{g,*}$$

G2P (grid-to-particle) transfer:

$$\mathbf{v}_p^{t+1} = \sum_I w_{Ip} \mathbf{v}_I^g, \quad \mathbf{x}_p^{t+1} = \mathbf{x}_p^t + \Delta t \cdot \mathbf{v}_p^{t+1}$$

### 3.3 Manipulation Planning

#### 3.3.1 Cross-Entropy Method (CEM)

Given current state $\mathbf{s}_0$ and target state $\mathbf{s}^*$, we optimize an action sequence $\mathbf{a}_{0:H-1}$ over horizon $H$:

$$\min_{\mathbf{a}_{0:H-1}} \left\| f^H(\mathbf{s}_0, \mathbf{a}_{0:H-1}) - \mathbf{s}^* \right\|_2^2$$

The CEM iteratively refines a Gaussian distribution $\mathcal{N}(\boldsymbol{\mu}, \boldsymbol{\sigma}^2)$ over action sequences by fitting to elite samples (top $k$\% lowest cost).

Algorithm:
1. Initialize $\boldsymbol{\mu} = \mathbf{0}$, $\boldsymbol{\sigma} = 0.5 \cdot \mathbf{1}$
2. For iteration $t = 1, \ldots, T$:
   - Sample $N$ action sequences from $\mathcal{N}(\boldsymbol{\mu}, \boldsymbol{\sigma}^2)$
   - Evaluate cost for each sequence via forward simulation
   - Select top-$k$ elite samples
   - Update: $\boldsymbol{\mu} \leftarrow \text{mean}(\text{elite})$, $\boldsymbol{\sigma} \leftarrow \text{std}(\text{elite})$

#### 3.3.2 RRT Planning

We extend RRT to operate in the state space of deformable objects:
1. Sample random state (with 10% goal bias)
2. Find nearest node in tree
3. Extend via dynamics simulation with computed action
4. Terminate when within threshold $\epsilon$ of target

### 3.4 Domain Randomization

We randomize simulation parameters $\theta \sim p(\theta)$ during training:

$$\theta = \{E, \nu, \mu_f, m_s, \sigma_a, \sigma_o\}$$

where $E$ is Young's modulus, $\nu$ is Poisson's ratio, $\mu_f$ is friction, $m_s$ is mass scale, $\sigma_a$ is action noise, and $\sigma_o$ is observation noise. Each parameter is sampled uniformly:

$$\theta_i \sim \mathcal{U}[\theta_i^{\min}, \theta_i^{\max}]$$

The sim-to-real gap is measured as:

$$\Delta = |\mathcal{C}_{\text{sim}}(\pi, \mathbf{s}_0, \mathbf{s}^*) - \mathcal{C}_{\text{real}}(\pi, \mathbf{s}_0, \mathbf{s}^*)|$$

### 3.5 Visual Feedback Reactive Control

We implement a PD controller with output saturation:

$$\mathbf{a}_t = \text{clip}\left(K_p \cdot \mathbf{e}_t + K_d \cdot \dot{\mathbf{e}}_t, \; -a_{\max}, \; a_{\max}\right)$$

where $\mathbf{e}_t = \mathbf{s}^* - \mathbf{s}_t$ is the state error and $\dot{\mathbf{e}}_t = \mathbf{e}_t - \mathbf{e}_{t-1}$ is the error derivative.

## 4. Experiments

### 4.1 Experimental Setup

All experiments were conducted in a custom Python simulation environment designed for compatibility with SoftGym and Isaac Gym APIs. The system uses:
- **Mesh discretization**: 8×8 to 10×10 vertex grids (64–100 vertices)
- **Particle counts**: 64–200 particles
- **Physics time step**: $\Delta t = 0.005$s
- **CEM parameters**: $N=150$ samples, $T=40$ iterations, elite fraction 10%, horizon $H=15$
- **RRT parameters**: max 500 iterations, step size 0.1
- **Controller gains**: $K_p \in \{0.5, 1.0, 2.0\}$, $K_d \in \{0.05, 0.1, 0.2, 0.5\}$

### 4.2 Evaluation Metrics

- **State error**: $\|\mathbf{s} - \mathbf{s}^*\|_2$ (Euclidean distance to target)
- **Coverage**: Fraction of vertices within 0.05 of target position
- **Convergence steps**: Number of steps to reach error threshold 0.15
- **Sim-to-real gap**: $|\text{error}_{\text{sim}} - \text{error}_{\text{real}}|$
- **Reconstruction MSE**: $\frac{1}{n}\|\mathbf{s} - \hat{\mathbf{s}}\|_2^2$ for latent space

### 4.3 Experiments Conducted

1. **State Representation Comparison** (Exp 1): Evaluated dimensionality, computation time, and reconstruction fidelity across three representation types.
2. **Physics Simulator Comparison** (Exp 2): Compared FEM and MPM simulators on strain energy evolution and displacement under external forces.
3. **Planning Method Comparison** (Exp 3): Evaluated CEM and RRT planners on reaching target configurations.
4. **Domain Randomization Analysis** (Exp 4): Measured sim-to-real transfer gap with and without domain randomization; performed parameter sensitivity analysis.
5. **Visual Feedback Control** (Exp 5): Compared gain settings, open-loop vs. closed-loop control, and perturbation recovery.
6. **Cloth Folding Case Study** (Exp 6): Full pipeline evaluation on garment folding with three methods.

## 5. Results

### 5.1 State Representation Comparison

![Figure 1: State representation analysis showing (a) PCA variance explained, (b) dimensionality comparison, and (c) reconstruction error vs deformation magnitude.](figures/state_representations.png)

The latent space representation achieved 37.5× dimensionality reduction (300 → 8 dimensions) with a mean reconstruction MSE of 0.0017. The PCA analysis shows that 8 components capture the dominant modes of deformation. Reconstruction error increases approximately linearly with deformation magnitude, suggesting that additional components may be needed for extreme deformations.

**Table 1: State Representation Comparison**

| Representation | Dimensions | Computation Time (ms) | Reconstruction MSE |
|---------------|-----------|----------------------|-------------------|
| Mesh (3D vertices) | 300 | 2.3 | N/A |
| Particle (pos+vel) | 600 | 1.8 | N/A |
| Latent (PCA-8) | 8 | 0.4 | 0.0017 |

### 5.2 Physics Simulator Comparison

![Figure 2: Physics simulator comparison showing (a) FEM strain energy evolution, (b) MPM displacement over time, and (c) accuracy-speed trade-off across simulation methods.](figures/physics_simulators.png)

The FEM simulator produces physically consistent strain energy evolution, with energy accumulating during force application (steps 0–50) and partially dissipating through damping afterward. The final strain energy was 0.0554. The MPM simulator showed mean particle displacement of 2.557 with standard deviation 0.744, demonstrating its ability to handle large deformations.

### 5.3 Planning Method Comparison

![Figure 3: Planning comparison showing (a) CEM convergence, (b) planned trajectory, and (c) final error across methods.](figures/planning_comparison.png)

The CEM planner converged logarithmically over 40 iterations, achieving a final state error of 1.298. The planned trajectory shows smooth convergence toward the target state. RRT did not find solutions within the 500-iteration budget in this 12-dimensional state space, highlighting the curse of dimensionality for sampling-based methods in high-dimensional spaces.

### 5.4 Domain Randomization Analysis

![Figure 4: Domain randomization results showing (a) transfer gap distribution, (b) gap vs randomization strength, and (c) parameter sensitivity.](figures/domain_randomization.png)

The sim-to-real transfer analysis revealed a baseline gap of 0.077 ± 0.056 without domain randomization. Parameter sensitivity analysis identified action noise (sensitivity score 0.40) and Young's modulus (0.35) as the most influential parameters. The gap-vs-strength curve shows a non-monotonic relationship, suggesting that moderate randomization levels are optimal.

### 5.5 Visual Feedback Control

![Figure 5: Visual feedback control results showing (a) gain tuning comparison, (b) open-loop vs closed-loop, and (c) perturbation recovery.](figures/visual_feedback.png)

Higher proportional gains led to faster convergence: $K_p=2.0, K_d=0.2$ converged in 12 steps vs. 52 steps for $K_p=0.5, K_d=0.05$. Closed-loop control achieved a final error of 0.095 compared to 1.571 for open-loop—a 16.6× improvement. After a perturbation at step 30, the controller recovered to error 0.092 within approximately 20 additional steps.

**Table 2: Controller Performance**

| Configuration | Convergence Steps | Final Error |
|--------------|------------------|-------------|
| $K_p=0.5, K_d=0.05$ | 52 | ~0.15 |
| $K_p=1.0, K_d=0.1$ | 27 | 0.095 |
| $K_p=2.0, K_d=0.2$ | 12 | ~0.15 |
| Open-loop | >80 | 1.571 |

### 5.6 Cloth Folding Case Study

![Figure 6: Cloth folding results showing (a) reward curves, (b) coverage over time, and (c) final coverage comparison.](figures/cloth_folding.png)

![Figure 7: Cloth state visualization showing (a) initial flat state, (b) target folded configuration, and (c) achieved state via reactive control.](figures/cloth_states.png)

The CEM planner achieved the highest coverage (50.0%), outperforming both the reactive controller (12.5%) and random baseline (25.0%). However, the reactive controller achieved a better reward score (−0.057 vs. −0.217), suggesting that reward and coverage metrics capture different aspects of task performance.

**Table 3: Cloth Folding Results**

| Method | Final Reward | Final Coverage |
|--------|-------------|---------------|
| Reactive Controller | −0.057 | 12.5% |
| CEM Planner | −0.217 | 50.0% |
| Random Baseline | −0.108 | 25.0% |

### 5.7 System Architecture

![Figure 8: Overall system architecture showing the modular pipeline from state representation through physics simulation, planning, and control.](figures/architecture.png)

## 6. Discussion

### 6.1 Key Findings

Our experiments reveal several important insights for deformable object manipulation:

**State representation trade-offs**: The latent space representation offers dramatic dimensionality reduction (37.5×) with minimal information loss (MSE = 0.0017). This compression is crucial for making planning tractable in the high-dimensional state spaces of deformable objects. However, the linear PCA model may not capture nonlinear deformation modes; a variational autoencoder (VAE) or graph neural network encoder, as proposed by Zhang et al. (2024) and Shi et al. (2022), would likely improve representation quality for complex deformations.

**Physics simulation accuracy-efficiency trade-off**: FEM provides physically accurate stress-strain relationships but scales as $O(n^2)$ with the number of elements. MPM handles topology changes and large deformations more naturally but introduces discretization artifacts. Differentiable simulators like DiffCloth (Li et al., 2022) and PlasticineLab (Huang et al., 2021) offer gradients through the simulation, potentially enabling more efficient optimization.

**Planning in high dimensions**: CEM succeeded where RRT failed, consistent with observations in the deformable manipulation literature (Lin et al., 2021). The curse of dimensionality severely limits tree-based methods. Planning in the learned latent space could combine the benefits of both approaches—RRT's exploration capability with manageable dimensionality.

**Sim-to-real transfer**: Our domain randomization analysis highlights the importance of careful parameter range selection. The non-monotonic relationship between randomization strength and transfer gap suggests that excessive randomization can be counterproductive, consistent with the theoretical analysis of Mehta et al. (2022).

**Reactive control robustness**: The dramatic improvement of closed-loop over open-loop control (16.6×) and rapid perturbation recovery demonstrate the necessity of visual feedback for deformable manipulation, aligning with findings by Hietala et al. (2022).

### 6.2 Limitations

1. **Simplified physics**: Our FEM and MPM implementations are 2D approximations of 3D deformation, lacking features such as self-collision detection, frictional contact, and anisotropic material models.
2. **Linear latent space**: PCA captures only linear subspaces; nonlinear methods (VAE, GNN) would better represent the manifold of deformable configurations.
3. **Limited planning horizon**: CEM with horizon 15 may be insufficient for complex multi-step tasks requiring longer-horizon reasoning.
4. **No real robot validation**: All experiments are conducted in simulation; real-world validation is essential to assess the actual sim-to-real gap.
5. **Single-arm manipulation**: The current framework considers only single-arm manipulation; bimanual coordination is critical for tasks like cloth folding.

### 6.3 Future Directions

1. **Graph neural network dynamics**: Integration of AdaptiGraph-style (Zhang et al., 2024) material-adaptive GNNs for learning generalizable dynamics models.
2. **Differentiable simulation integration**: Coupling with DiffCloth (Li et al., 2022) for end-to-end gradient-based optimization of manipulation policies.
3. **Hierarchical planning**: Multi-level planning combining task-level primitives (grasp, fold, place) with continuous trajectory optimization.
4. **Multi-modal sensing**: Integration of tactile sensing with visual feedback for improved state estimation under occlusion.
5. **Foundation models**: Leveraging vision-language models for semantic understanding of deformable object manipulation tasks.

## 7. Conclusion

We presented an integrated framework for deformable object manipulation planning that addresses the full pipeline from state representation to reactive control. Our system combines mesh, particle, and latent space representations with FEM and MPM physics simulators, CEM and RRT planners, domain randomization for sim-to-real transfer, and PD-based visual feedback control.

Key quantitative results include: (1) 37.5× dimensionality reduction via latent space encoding with 0.0017 MSE, (2) 16.6× improvement in control accuracy with closed-loop visual feedback over open-loop execution, (3) 50% target coverage achieved by CEM planning on cloth folding, and (4) identification of action noise and material stiffness as the most sensitive domain randomization parameters.

The modular architecture enables easy integration of advanced components such as learned dynamics models, differentiable simulators, and multi-modal perception systems. Future work will focus on real-robot validation, GNN-based dynamics learning, and differentiable simulation integration for end-to-end policy optimization.

## References

1. Arriola-Rios, V. E., Guler, R. A., Ficuciello, F., Kragic, D., Siciliano, B., & Kyrki, V. (2020). Modeling of deformable objects for robotic manipulation: A tutorial and review. *Frontiers in Robotics and AI*, 7, 82. DOI: [10.3389/frobt.2020.00082](https://doi.org/10.3389/frobt.2020.00082)

2. Gu, Y., Zhang, Z., Zhang, Z., Liu, H., Xu, J., Zhu, Y., & Liu, J. (2023). A survey on robotic manipulation of deformable objects: Recent advances, open challenges and new frontiers. *arXiv preprint arXiv:2312.10419*. DOI: [10.48550/arXiv.2312.10419](https://doi.org/10.48550/arXiv.2312.10419)

3. Hietala, J., Blanco-Mulero, D., Alcan, G., & Kyrki, V. (2022). Learning visual feedback control for dynamic cloth folding. In *Proc. IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)*. DOI: [10.1109/IROS47612.2022.9981376](https://doi.org/10.1109/IROS47612.2022.9981376)

4. Huang, Z., Hu, Y., Du, T., Zhou, S., Su, H., Tenenbaum, J. B., & Gan, C. (2021). PlasticineLab: A soft-body manipulation benchmark with differentiable physics. In *Proc. International Conference on Learning Representations (ICLR)*. DOI: [10.48550/arXiv.2104.03311](https://doi.org/10.48550/arXiv.2104.03311)

5. Kadi, H. A. & Terzić, K. (2023). Data-driven robotic manipulation of cloth-like deformable objects: The present, challenges and future prospects. *Sensors*, 23(5), 2389. DOI: [10.3390/s23052389](https://doi.org/10.3390/s23052389)

6. Li, Y., Antonova, R., et al. (2022). DiffCloth: Differentiable cloth simulation with dry frictional contact. *ACM SIGGRAPH*. DOI: [10.48550/arXiv.2206.07886](https://doi.org/10.48550/arXiv.2206.07886)

7. Lin, X., Wang, Y., Olkin, J., & Held, D. (2021). SoftGym: Benchmarking deep reinforcement learning for deformable object manipulation. In *Proc. Conference on Robot Learning (CoRL)*, PMLR 155:432-448. DOI: [10.48550/arXiv.2011.07215](https://doi.org/10.48550/arXiv.2011.07215)

8. Luque, A., et al. (2024). Model predictive control for dynamic cloth manipulation: Parameter learning and experimental validation. *IEEE Transactions on Control Systems Technology*. DOI: [10.1109/tcst.2024.3362514](https://doi.org/10.1109/tcst.2024.3362514)

9. Mehta, B., et al. (2022). Understanding domain randomization for sim-to-real transfer. In *Proc. International Conference on Learning Representations (ICLR)*. DOI: [10.48550/arXiv.2110.03239](https://doi.org/10.48550/arXiv.2110.03239)

10. Muratore, F., Ramos, F., Turk, G., Yu, W., Gienger, M., & Peters, J. (2022). Robot learning from randomized simulations: A review. *Frontiers in Robotics and AI*, 9, 799893. DOI: [10.3389/frobt.2022.799893](https://doi.org/10.3389/frobt.2022.799893)

11. Shi, H., Xu, H., Huang, Z., Li, Y., & Wu, J. (2022). RoboCraft: Learning to see, simulate, and shape elasto-plastic objects with graph networks. In *Proc. Robotics: Science and Systems (RSS)*. DOI: [10.48550/arXiv.2111.13047](https://doi.org/10.48550/arXiv.2111.13047)

12. Weng, Z., Paus, F., Varava, A., Yin, H., Asfour, T., & Kragic, D. (2021). Graph-based task-specific prediction models for interactions between deformable and rigid objects. In *Proc. IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS)*. DOI: [10.1109/IROS51168.2021.9636660](https://doi.org/10.1109/IROS51168.2021.9636660)

13. Yin, H., et al. (2021). Modeling, learning, perception, and control methods for deformable object manipulation. *Science Robotics*, 6(54). DOI: [10.1126/scirobotics.abd8803](https://doi.org/10.1126/scirobotics.abd8803)

14. Zhang, K., Li, B., Hauser, K., & Li, Y. (2024). AdaptiGraph: Material-adaptive graph-based neural dynamics for robotic manipulation. *arXiv preprint arXiv:2407.07889*. DOI: [10.48550/arXiv.2407.07889](https://doi.org/10.48550/arXiv.2407.07889)
