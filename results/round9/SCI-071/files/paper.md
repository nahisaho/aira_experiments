# Deformable Object Manipulation Planning via Mass-Spring Simulation, Latent State Representation, and Domain-Randomized Sim-to-Real Transfer

---

## Abstract

Robotic manipulation of deformable objects—including cloth, rope, and elastic materials—remains one of the most challenging open problems in robotics, owing to their infinite-dimensional configuration spaces, complex contact dynamics, and large sim-to-real gaps. This paper presents a comprehensive manipulation planning system for deformable objects designed around four interconnected modules: (1) a mass-spring mesh state representation with PCA-based latent compression, (2) physics simulation via an explicit semi-implicit Euler integrator, (3) trajectory optimization through the Cross-Entropy Method (CEM), and (4) domain-randomized sim-to-real transfer with supervised performance prediction.

We implement and evaluate each component on a cloth-folding benchmark task using an 8×8 mass-spring cloth mesh (64 nodes, 306 springs) and a 20-link deformable linear object (DLO/rope). PCA latent compression reduces the 128-dimensional cloth state to 2 principal components capturing 95.3% of variance [cell:4], demonstrating the intrinsic low-dimensionality of cloth configurations. CEM planning achieves a 19.0% reduction in mean node position error (from 0.2064 m to 0.1672 m) over 30 iterations [cell:6]. Domain randomization over stiffness (±30%), friction (±20%), and mass (±20%) yields a statistically significant 21.2% improvement in robustness (t=15.52, p<0.001) [cell:7]. A Gradient Boosting Tree model trained on 600 synthetic sim-to-real experiments predicts real-world task performance with R²=0.796±0.021 and RMSE=0.069±0.003 (5-fold cross-validation) [cell:10b], identifying simulation stiffness (importance=0.570) as the dominant factor governing sim-to-real transfer quality.

These results are obtained from reproducible simulations (seed=42) and are intentionally contextualized with their limitations: the mass-spring model is a linear approximation of real cloth physics, the "real-world performance" targets are synthetically generated, and the 19% CEM improvement is modest—reflecting the difficulty of cloth folding even in simulation. We discuss how these findings align with and diverge from prior benchmarks (SoftGym, GarmentLab) and provide a roadmap for future work integrating FEM/MPM simulators, neural state estimators, and end-to-end sim-to-real pipelines.

**Keywords**: deformable object manipulation, sim-to-real transfer, domain randomization, Cross-Entropy Method, mass-spring simulation, cloth folding, latent state representation

---

## 1. Introduction

Deformable object manipulation (DOM) is a fundamental challenge in robotics with applications in textile manufacturing, surgical assistance, household automation, and logistics. Unlike rigid objects, deformable objects exhibit infinite-dimensional configuration spaces governed by complex material properties (stiffness, friction, damping), making it difficult to apply classical planning algorithms directly.

Recent advances have leveraged deep reinforcement learning (RL) with physics simulators to learn manipulation policies. SoftGym [Lin et al., 2021] introduced a benchmark suite using particle-based simulation for cloth, rope, and fluid manipulation, demonstrating that RL agents can acquire manipulation skills in simulation. However, transferring these policies to real robots—the **sim-to-real gap**—remains a major barrier. Physical simulators cannot perfectly capture real-world material properties, sensor noise, or contact dynamics, leading to policy degradation when deployed.

This work addresses three interconnected research questions:
1. **State Representation**: What is the minimal latent dimensionality required to faithfully represent deformable object states for planning?
2. **Planning Efficiency**: How effectively can sampling-based trajectory optimization (CEM) solve cloth-folding under deformable dynamics?
3. **Sim-to-Real Gap**: Which simulation parameters most strongly predict real-world transfer performance, and how much does domain randomization (DR) help?

### 1.1 Contributions
- A complete, reproducible deformable manipulation pipeline implemented in Python (mass-spring cloth + rope simulation, PCA state compression, CEM planning, DR evaluation)
- Quantitative analysis of PCA latent space structure across stiffness configurations, showing 95.3% variance explained in 2 PCs
- Domain randomization experiments demonstrating 21.2% robustness improvement (p<0.001)
- Gradient Boosting model for predicting sim-to-real performance, identifying stiffness calibration as the dominant factor (importance=0.570)
- Self-critical analysis of limitations of mass-spring models vs. FEM/MPM approaches

---

## 2. Related Work

### 2.1 Deformable Object Simulation

The choice of simulation model critically affects both fidelity and computational cost. Three main paradigms exist:

**Mass-Spring Systems (MSS)**: Approximate deformable objects as networks of point masses connected by Hookean springs [Nealen et al., 2006]. Computationally efficient and easy to implement, but lack accuracy for large deformations and do not satisfy constitutive equations of continuum mechanics. Used as the primary simulation model in this work.

**Finite Element Method (FEM)**: Discretizes the object into finite elements and solves continuum mechanics equations. Provides physically accurate simulation for elastic solids and thin-shell cloth [Li et al., 2024; DOI: 10.1145/3641519.3657449]. Higher computational cost; used in GarmentLab [Lu et al., 2023].

**Material Point Method (MPM)**: Hybrid Lagrangian-Eulerian method well-suited for highly deformable or fracturing materials [Li et al., 2024]. Used in FlexibleNet and similar systems.

### 2.2 Prior Work on Deformable Manipulation

**SoftGym** [Lin et al., 2021]: Established benchmark using Position-Based Dynamics (PBD) via PyFlex for cloth, rope, and fluid tasks. Key finding: RL agents can solve cloth-spreading, rope-routing, and fluid-pouring in simulation, but sim-to-real transfer requires explicit adaptation.

**Visual Feedback for Cloth Folding** [Hietala et al., 2022; DOI: 10.1109/IROS47612.2022.9981376]: Visual feedback control with MuJoCo simulation, demonstrating dynamic cloth folding with successful real-robot transfer. Key contribution: closed-loop visual feedback significantly outperforms open-loop planned trajectories.

**Cloth Manipulation Planning via Mesh Representations** [Frontiers Neurorobotics, 2023; DOI: 10.3389/fnbot.2022.1045747]: Neural network estimating mesh representations from voxel input with uncertainty-aware planning. Key finding: mesh representation outperforms keypoint representations for planning quality.

**Bayesian Sim-to-Real** [Antonova et al., 2022]: Bayesian approach to real-to-sim adaptation, using posterior updates over simulation parameters to minimize reality gap. Demonstrates that stiffness and friction are the most critical parameters for real-robot transfer—consistent with our findings.

**DextAIRity** [Xu et al., 2022; arXiv:2206.04091]: Non-contact manipulation using airflow, achieving state-of-the-art performance on cloth/bag manipulation. Highlights that contact-free approaches can circumvent sim-to-real contact modeling challenges.

**Diffusion Dynamics Models** [CoRL 2023; arXiv:2308.05444]: Transformer-based diffusion model for cloth state estimation and dynamics, enabling model-predictive control with uncertainty quantification.

### 2.3 Domain Randomization for Sim-to-Real

Domain randomization (DR) [Tobin et al., 2017] has become the dominant technique for sim-to-real transfer, randomly varying simulation parameters during training to expose the policy to a wide range of conditions. For deformable objects, key parameters include material stiffness, friction coefficients, gravity, and simulation step size. Isaac Gym and OmniIsaacGymEnvs provide native DR support through YAML-configurable parameter randomization [NVIDIA, 2023].

### 2.4 Research Gaps Addressed

Existing work lacks: (1) systematic analysis of which simulation parameters dominate sim-to-real transfer for cloth specifically, (2) quantitative comparison of latent compression methods, and (3) open-source tools for predicting sim-to-real performance without running expensive real-robot experiments. This work provides preliminary evidence on all three.

---

## 3. Methods

### 3.1 Mass-Spring Cloth Model

We implement a 2D mass-spring cloth mesh with $N = r \times c$ nodes (default: $8 \times 8 = 64$ nodes). The cloth is discretized into structural, shear, and bending springs:

- **Structural springs**: Connect adjacent nodes horizontally/vertically, rest length $\ell_0$, stiffness $k_s$
- **Shear springs**: Connect diagonal neighbors, rest length $\ell_0\sqrt{2}$, stiffness $k_{sh}$
- **Bending springs**: Connect nodes two steps apart, rest length $2\ell_0$, stiffness $k_b$

The spring force between nodes $i$ and $j$ is:
$$\mathbf{f}_{ij} = k \left( \|\mathbf{x}_j - \mathbf{x}_i\| - \ell_{ij}^0 \right) \frac{\mathbf{x}_j - \mathbf{x}_i}{\|\mathbf{x}_j - \mathbf{x}_i\|}$$

Integration uses semi-implicit Euler with velocity damping:
$$\mathbf{v}_{t+1} = \alpha \left(\mathbf{v}_t + \Delta t \cdot \mathbf{a}_t \right), \quad \mathbf{x}_{t+1} = \mathbf{x}_t + \Delta t \cdot \mathbf{v}_{t+1}$$

where $\alpha = 0.98$ is the damping coefficient and $\Delta t = 0.005$ s.

**Default parameters**: $k_s = 100$ N/m, $k_{sh} = 50$ N/m, $k_b = 10$ N/m, node mass $m = 0.1$ kg.

### 3.2 PCA Latent State Representation

Given a dataset of $M$ cloth states $\{\mathbf{s}^{(i)}\}_{i=1}^M$ where $\mathbf{s}^{(i)} \in \mathbb{R}^{2N}$ (XY positions of all nodes), we apply PCA to obtain a latent embedding $\mathbf{z}^{(i)} = V_d^T (\mathbf{s}^{(i)} - \bar{\mathbf{s}})$ where $V_d$ contains the top-$d$ principal components.

States are collected from 5 cloth configurations with different stiffness settings (Soft/Medium/Stiff/Med-Soft/Med-Stiff), generating 125 state snapshots. StandardScaler normalization is applied before PCA.

### 3.3 Cross-Entropy Method (CEM) Planning

CEM treats trajectory planning as a stochastic optimization problem [Rubinstein, 1997]. An action sequence $\mathbf{a} = (a_1, \ldots, a_K)$ with $a_k = (\hat{n}_k, \Delta x_k, \Delta y_k)$ (normalized node index, displacement components) is optimized to minimize mean node position error to target state.

**Algorithm**:
1. Initialize $\mu = \mathbf{0}$, $\sigma = 0.3 \cdot \mathbf{1}$ (action distribution)
2. Sample $N_{pop}$ candidates from $\mathcal{N}(\mu, \text{diag}(\sigma^2))$
3. Evaluate cost: $C(\mathbf{a}) = \frac{1}{N} \sum_{n=1}^N \|\hat{\mathbf{x}}_n - \mathbf{x}_n^*\|_2$ where $\hat{\mathbf{x}}$ is final simulated state
4. Select top $N_e = 0.2 N_{pop}$ elite samples; update $\mu, \sigma$
5. Apply noise decay: $\sigma \leftarrow 0.97 \sigma$
6. Repeat for $T=30$ iterations

**Parameters**: $N_{pop} = 60$, $N_e = 12$, $K = 5$ actions, $T = 30$ iterations.

The action effect model uses Gaussian influence: displaced node $n$ propagates displacement $(\Delta x, \Delta y)$ to all nodes with influence $\exp(-d_{in}/0.3)$ where $d_{in}$ is the Euclidean distance from node $i$ to the grasped node $n$.

### 3.4 Domain Randomization Protocol

We evaluate policy robustness under domain randomization matching the IsaacGym DR convention:

| Parameter | Baseline | DR Range |
|-----------|----------|----------|
| Stiffness | 1.0 (normalized) | Uniform[0.7, 1.3] |
| Mass | 1.0 (normalized) | Uniform[0.8, 1.2] |
| Friction | 1.0 (normalized) | Uniform[0.8, 1.2] |
| Perception noise | 0 | σ=0.01 m |

For each of 100 DR trials, simulation parameters are drawn independently, the planned action sequence is executed, and the final state error is recorded.

### 3.5 Sim-to-Real Performance Prediction Model

We generate a synthetic dataset of $N=600$ simulated experiments with features:
- Simulation stiffness (30–350 N/m), friction (0.1–0.9), gravity (8.5–11.0 m/s²)
- DR range magnitude (0.05–0.60)
- Number of parallel environments (50–3000)
- Training steps (10k–800k)

The target variable is "real-world task completion rate" synthesized via a physically-motivated function incorporating exponential stiffness matching, logarithmic training saturation, and additive Gaussian noise (σ=0.06). Three models are compared with 5-fold cross-validation: Ridge regression, Random Forest (100 trees, max_depth=8), and Gradient Boosting Trees (GBT, 100 estimators, lr=0.1, max_depth=4).

### 3.6 NatureLM and GALACTICA MCP Tools

**Attempted tools and outcomes**:
- **NatureLM MCP** (`ask_naturelm`): Tool not available in the current ToolUniverse environment. Search query `ask_naturelm, NatureLM, naturelm` returned 0 matches. No quantitative predictions could be obtained.
- **GALACTICA MCP** (`scientific_qa`, `predict_citations`): Tool not available. Search for `galactica, scientific_qa, predict_citations` returned 0 matches.
- **Semantic Scholar API** (`SemanticScholar_search_papers`): Attempted 4 queries, all returned HTTP 429 (rate limit exceeded). Literature survey was conducted via web search as alternative.

These tool failures are documented for scientific transparency. The experimental design and parameter choices are grounded in the literature survey findings rather than model-generated predictions.

### 3.7 Python Implementation

All experiments use Python 3.11.2 with numpy (random_state=42 / np.random.seed(42)) for full reproducibility.

```python
# Core simulation (cell 2)
class MassSpringCloth:
    def __init__(self, rows=8, cols=8, rest_length=0.1,
                 k_structural=100.0, k_shear=50.0, k_bend=10.0,
                 damping=0.98, mass=0.1, dt=0.005):
        ...

# CEM Planner (cell 6)
class FastClothPlanner:
    def plan_cem(self, n_iter=30, population=60, elite_frac=0.2):
        mu = np.zeros(action_dim); sigma = np.ones(action_dim) * 0.3
        for iter in range(n_iter):
            samples = mu + sigma * np.random.randn(population, action_dim)
            costs = [self.cost(s) for s in samples]
            elite = samples[np.argsort(costs)[:n_elite]]
            mu, sigma = elite.mean(0), elite.std(0) + 1e-6
            sigma *= 0.97

# ML Model (cell 10b)
gb_model = GradientBoostingRegressor(n_estimators=100, max_depth=4,
                                      learning_rate=0.1, random_state=42)
cv_scores = cross_val_score(gb_model, X, y, cv=KFold(5, shuffle=True, random_state=42))
```

---

## 4. Experiments

### 4.1 Simulation Setup

All experiments run in a Jupyter kernel (Python 3.11.2, numpy 2.3.5, scikit-learn 1.6.1, scipy 1.17.1). No GPU is required; computation runs on CPU. Total runtime is under 2 minutes.

**Cloth benchmark task**: Fold an 8×8 cloth mesh in half along the horizontal axis (bottom 4 rows folded onto top 4 rows). The task evaluates both the planning module (CEM) and the DR robustness of the resulting policy.

**Rope benchmark**: Simulate a 20-link rope hanging under gravity from a fixed endpoint, serving as a deformable linear object (DLO) test case.

### 4.2 Evaluation Metrics

- **Mean node position error (m)**: $\frac{1}{N}\sum_{n=1}^N \|\hat{\mathbf{x}}_n - \mathbf{x}_n^*\|_2$ — primary planning metric
- **PCA explained variance ratio**: $\sum_{k=1}^d \lambda_k / \sum_{k=1}^N \lambda_k$
- **Model R²** (cross-validation): coefficient of determination for sim-to-real prediction
- **Robustness improvement**: relative reduction in mean error under DR
- **Pearson correlation**: between stiffness ratio and final error

### 4.3 Domain Randomization Ranges

Matching Isaac Gym / SoftGym conventions:
- Stiffness: ±30% (represents typical manufacturing tolerance for textile materials)
- Friction: ±20% (surface finish variation)
- Mass: ±20% (density variation in cloth)

---

## 5. Results

### 5.1 Cloth Simulation [cell:2]

The 8×8 mass-spring cloth mesh (64 nodes, 306 springs) correctly exhibits gravity-induced draping behavior. Under free-hang boundary conditions (top two corners fixed), the maximum node displacement is **1.3799 m** after 300 simulation steps (dt=5ms, total 1.5 s). The cloth reaches a stable catenary-like equilibrium, validating the mass-spring implementation.

![Figure 1: Cloth simulation snapshots](figures/fig1_cloth_simulation.png)
**Figure 1**: Mass-spring cloth simulation snapshots at t=0, 50, 150, 300 steps. Red triangles indicate fixed nodes. Structural springs (k_s=100 N/m) maintain cloth connectivity while shear springs (k_sh=50 N/m) prevent unrealistic in-plane deformation.

### 5.2 Latent State Representation [cell:4, cell:5]

PCA applied to 125 cloth state snapshots from 5 stiffness configurations reveals strong low-dimensional structure:

| Principal Components | Explained Variance |
|---------------------|--------------------|
| PC1 | **83.2%** |
| PC1 + PC2 | **95.3%** |
| PC1 + PC2 + PC3 | **98.3%** |
| PC1 through PC10 | 100.0% |

The 2D latent space cleanly separates stiffness conditions (Fig. 2), confirming that cloth deformation modes are dominated by 1–2 primary eigenmodes (rigid-body-like translation and first bending mode). This finding supports the use of 2–4 dimensional latent representations for planning, consistent with [Frontiers Neurorobotics, 2023].

![Figure 2: Latent space visualization](figures/fig2_latent_space.png)
**Figure 2**: Left: PCA latent space (PC1 vs PC2) colored by stiffness category. Right: explained variance by component. PC1 and PC2 together capture 95.3% of variance [cell:4], justifying low-dimensional state representations for cloth manipulation planning.

### 5.3 CEM Planning Results [cell:6]

CEM planning with 30 iterations (population=60, elite=12, 5 actions):

| Metric | Value |
|--------|-------|
| Initial cost (before planning) | 0.2064 m |
| Final best cost (after 30 iter) | **0.1672 m** |
| Cost reduction | **19.0%** |
| Convergence pattern | Monotonically decreasing |

The CEM optimizer achieves a 19.0% reduction in mean node position error for the cloth folding task [cell:6]. The convergence curve (Fig. 3, left) shows rapid improvement in the first 10 iterations followed by plateau, suggesting the action influence model (Gaussian propagation) saturates the available search space. Larger populations or more sophisticated action parameterizations (e.g., Bézier curves) would likely yield further improvement.

### 5.4 Domain Randomization Robustness [cell:7]

| Policy | Mean Error (m) | Std Error (m) |
|--------|---------------|---------------|
| Baseline (no DR) | 0.1694 | 0.0094 |
| DR-augmented | **0.1334** | 0.0212 |
| Improvement | **21.2%** | — |

Two-sample t-test confirms the improvement is statistically significant: t=15.52, p<0.001 [cell:7]. The DR policy achieves lower mean error at the cost of higher variance (std 0.0212 vs 0.0094), consistent with the known bias-variance tradeoff in domain randomization.

Stiffness ratio is strongly anti-correlated with performance error (Pearson r = −0.839, p<0.001), indicating that **higher stiffness ratios (stiffer cloth) lead to lower position errors** under the action influence model used. This is a model artifact reflecting how the Gaussian influence function interacts with stiffness scaling—not a universal physical law.

![Figure 3: Planning and domain randomization results](figures/fig3_planning_dr.png)
**Figure 3**: Left: CEM convergence curve. Center: Error distribution under DR (baseline vs DR-augmented). Right: Stiffness sensitivity scatter plot with Pearson r=−0.839.

### 5.5 Sim-to-Real Transfer Prediction [cell:10b, cell:11]

| Model | R² (5-fold CV) | RMSE (5-fold CV) |
|-------|---------------|-----------------|
| Ridge Regression | 0.367 ± 0.047 | 0.1215 ± 0.0071 |
| Random Forest | 0.756 ± 0.024 | 0.0753 ± 0.0041 |
| **Gradient Boosting** | **0.796 ± 0.021** | **0.0688 ± 0.0027** |

Test set performance (20% holdout): R²=0.769, RMSE=0.0704 [cell:11].

Feature importance analysis from the Gradient Boosting model:

| Feature | Importance |
|---------|-----------|
| **Stiffness (N/m)** | **0.570** |
| Friction | 0.204 |
| DR_Range | 0.086 |
| Gravity (m/s²) | 0.067 |
| N_Envs | 0.038 |
| Train_Steps | 0.035 |

Stiffness calibration dominates (importance=0.570), followed by friction (0.204) [cell:11]. Training steps and environment count have surprisingly small effects in this model, likely because the synthetic data was generated with stiffness and friction as primary determinants.

![Figure 4: ML analysis for sim-to-real prediction](figures/fig4_ml_analysis.png)
**Figure 4**: Left: Gradient Boosting predicted vs. actual performance (test R²=0.769). Center: Feature importance. Right: Model comparison (5-fold CV).

### 5.6 System Architecture and Method Comparison [cell:13]

![Figure 5: Comprehensive results](figures/fig5_comprehensive.png)
**Figure 5**: Top-left: Reactive control error (10 trials ± 1σ). Top-right: System architecture diagram. Bottom-left: Method comparison. Bottom-right: Performance heatmap (stiffness × DR range).

### 5.7 Rope/DLO Simulation [cell:14]

The 20-link rope simulation correctly models DLO hanging dynamics:
- Maximum tip Y displacement: −1.306 m (catenary-like profile)
- Final tip position: [0.759, −1.306] m from fixed endpoint

![Figure 6: Rope simulation](figures/fig6_rope_simulation.png)
**Figure 6**: Rope/DLO simulation snapshots (20 links, k=80 N/m). The rope hangs naturally under gravity from its fixed left endpoint.

---

## 6. Discussion

### 6.1 Interpretation of Results

**Latent representation**: The finding that 95.3% of cloth state variance is captured in 2 PCs [cell:4] is consistent with prior work showing that cloth configurations lie on low-dimensional manifolds. However, this result is obtained from quasi-static configurations of a single cloth mesh under gravity—dynamic manipulation trajectories would require more components. Real cloth manipulation systems (e.g., GarmentLab) typically use 8–32 dimensional representations.

**CEM planning**: The 19.0% cost reduction [cell:6] is modest. This reflects several limitations: (1) the Gaussian influence action model is a crude approximation of cloth physics, (2) CEM with 30 iterations and population 60 has limited search budget, (3) the cloth folding task requires coordinated multi-point manipulation that 5 sequential grasps cannot fully achieve. State-of-the-art methods (e.g., diffusion-based planning [arXiv:2308.05444]) achieve much higher folding completion rates (~60–80%) using learned dynamics models.

**Domain randomization**: The DR improvement (21.2%, p<0.001) is statistically robust but quantitatively modest. The higher variance of the DR policy (std 0.0212 vs 0.0094 for baseline) reflects the known cost of DR: by training for many conditions, the policy becomes less specialized. In real applications, adaptive domain randomization (ADR) [OpenAI, 2019] adjusts DR ranges during training to avoid policy degradation.

**Stiffness dominance in sim-to-real**: The Gradient Boosting model assigns 57% of importance to stiffness, consistent with Antonova et al. [2022]'s finding that stiffness calibration is critical for real-robot transfer. However, the synthetic data generation function was designed with stiffness as a primary factor—this creates circular reasoning. A properly designed study would require actual real-robot experiments across varied stiffness conditions.

### 6.2 NatureLM and GALACTICA Tool Assessment

Both NatureLM MCP (`ask_naturelm`) and GALACTICA MCP (`scientific_qa`, `predict_citations`) were unavailable in the current ToolUniverse environment, returning 0 search matches. This prevented:
- NatureLM quantitative predictions (e.g., expected cloth deformation parameters, stiffness ranges for cotton/polyester)
- GALACTICA citation predictions for validating experimental design
- Cross-validation of our results against model-generated benchmarks

As a result, experimental parameters were set based on the literature survey (SoftGym, GarmentLab, Antonova et al.) rather than model predictions. The planned cross-validation between NatureLM and GALACTICA outputs could not be performed.

### 6.3 Critical Limitations

1. **Synthetic data**: All "real-world performance" values are generated by a mathematical function, not measured on physical robots. Conclusions about which parameters matter for real transfer cannot be validated without actual experiments.

2. **Mass-spring vs. FEM/MPM**: The mass-spring model lacks proper constitutive modeling, fails to conserve energy correctly, and cannot model large deformations accurately. FEM [Li et al., 2024] or PBD (SoftGym) would be more physically accurate.

3. **Action model simplification**: The Gaussian influence propagation used in CEM planning is a first-order approximation of cloth dynamics. A learned neural dynamics model would provide better predictions at the cost of training data.

4. **2D cloth simulation**: We simulate cloth in 2D (XY positions only). Real cloth manipulation requires 3D state representation with surface normals and contact forces.

5. **No real robot validation**: All results are from simulation. Sim-to-real transfer quality cannot be assessed without physical experiments.

### 6.4 Future Directions

1. **FEM/MPM integration**: Replace mass-spring with differentiable FEM (e.g., DiffTaichi, Warp) for physics-accurate gradients
2. **Neural state estimation**: Learn $\phi: I \to \mathbf{z}$ from RGB-D images to latent state
3. **Model-based RL**: Train neural dynamics model $f: (\mathbf{z}, \mathbf{a}) \to \mathbf{z}'$ for model-predictive control
4. **Isaac Gym/Lab integration**: Implement DR-based policy training using GPU-parallel simulation
5. **Real robot experiments**: Validate sim-to-real predictions on physical cloth folding platform

---

## 7. Conclusion

We present a deformable object manipulation planning system comprising mass-spring simulation, PCA latent state representation, CEM trajectory optimization, domain randomization robustness evaluation, and ML-based sim-to-real performance prediction. Key quantitative findings are:

- 2 PCA components capture 95.3% of cloth state variance [cell:4]
- CEM achieves 19.0% cost reduction over 30 iterations [cell:6]
- Domain randomization improves robustness by 21.2% (p<0.001) [cell:7]
- Gradient Boosting predicts sim-to-real transfer with R²=0.796 ± 0.021 [cell:10b]
- Stiffness calibration is the dominant sim-to-real factor (importance=0.570) [cell:11]

Critically, these results are obtained from fully synthetic experiments without real-robot validation. The practical significance of the 19% CEM improvement and 21% DR improvement is uncertain without comparison to state-of-the-art methods on standardized benchmarks (SoftGym, GarmentLab). Future work should integrate differentiable physics, neural state estimators, and GPU-parallel RL training to approach the performance of current state-of-the-art cloth manipulation systems.

---

## References

1. **Lin, X., Wang, Y., Olkin, J., & Held, D.** (2021). SoftGym: Benchmarking Deep Reinforcement Learning for Deformable Object Manipulation. *Proceedings of the 2020 Conference on Robot Learning (CoRL)*. PMLR 155:432-448. URL: https://proceedings.mlr.press/v155/lin21a.html

2. **Hietala, J., Blanco-Mulero, D., Alcan, G., & Kyrki, V.** (2022). Learning Visual Feedback Control for Dynamic Cloth Folding. *IEEE/RSJ IROS 2022*. DOI: 10.1109/IROS47612.2022.9981376

3. **Antonova, R., Yang, J., Sundaresan, P., Fox, D., Ramos, F., & Bohg, J.** (2022). A Bayesian Treatment of Real-to-Sim for Deformable Object Manipulation. *IEEE Robotics and Automation Letters (RA-L)*. URL: https://research.nvidia.com/labs/srl/publication/antonova-2022-a/

4. **De Gusseme, V.-L., & Wyffels, F.** (2022). Effective cloth folding trajectories in simulation with only two parameters. *Frontiers in Neurorobotics*, 16. DOI: 10.3389/fnbot.2022.989702

5. **Blanco-Mulero, D., Barbany, O., Alenyà, G., Kheddar, A., Moreno-Noguer, F., & Kyrki, V.** (2023). Cloth manipulation planning on basis of mesh representations with incomplete domain knowledge and voxel-to-mesh estimation. *Frontiers in Neurorobotics*. DOI: 10.3389/fnbot.2022.1045747

6. **Li, X., et al.** (2024). A Dynamic Duo of Finite Elements and Material Points. *ACM SIGGRAPH 2024*. DOI: 10.1145/3641519.3657449

7. **Xu, Z., Chi, C., Burchfiel, B., Cousineau, E., Feng, S., & Song, S.** (2022). DextAIRity: Deformable Manipulation Can be a Breeze. *Robotics: Science and Systems (RSS) 2022*. arXiv:2206.04091

8. **Lu, H., et al.** (2023/2024). GarmentLab: A Unified Simulation and Benchmark for Garment Manipulation. *NeurIPS 2024*. URL: https://garmentlab.github.io/

9. **Scheikl, P. M., et al.** (2023). Sim-To-Real Transfer for Visual Reinforcement Learning of Deformable Object Manipulation for Robot-Assisted Surgery. *IEEE Robotics and Automation Letters*. arXiv:2406.06092

10. **Mitrano, P., McConachie, D., & Berenson, D.** (2021). Learning where to trust unreliable models in an unstructured world for deformable object manipulation. *Science Robotics*, 6(54).

---

## Reproducibility

| Item | Value |
|------|-------|
| Random seed (numpy) | `np.random.seed(42)` |
| Random seed (Python) | `random.seed(42)` |
| Python version | 3.11.2 |
| numpy | 2.3.5 |
| scipy | 1.17.1 |
| scikit-learn | 1.6.1 |
| pandas | 2.3.3 |
| matplotlib | 3.10.9 |
| seaborn | 0.13.2 |
| Jupyter kernel | Python 3 (ipykernel) |
| Experiment data | `data/raw/experiment_results.json` |
| Notebook | `deformable_manipulation.ipynb` |

All simulation parameters are set as constants in the respective class constructors. The synthetic data generation functions use deterministic mathematical operations seeded at seed=42. Total computation time: < 2 minutes on CPU.
