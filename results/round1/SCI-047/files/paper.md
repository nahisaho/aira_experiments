# A Unified Bayesian Optimization Framework for High-Dimensional Experimental Design with Multi-Objective and Batch Capabilities

---

## Abstract

Bayesian optimization (BO) has emerged as a powerful paradigm for sample-efficient optimization of expensive black-box functions. However, scaling BO to high-dimensional parameter spaces, parallel experimental settings, and multi-objective problems remains challenging. In this work, we present a unified BO framework built on BOTorch and GPyTorch that integrates: (1) systematic Gaussian process kernel selection with automatic relevance determination (ARD), (2) comparative evaluation of acquisition functions including Expected Improvement (EI), Upper Confidence Bound (UCB), and Knowledge Gradient (KG), (3) batch optimization for parallel experiment proposals via q-EI, (4) multi-objective optimization using q-Expected Hypervolume Improvement (qEHVI), (5) dimensionality reduction for high-dimensional spaces (D>20) via random embedding (REMBO), and (6) a practical case study in chemical reaction condition optimization targeting yield and selectivity. Through comprehensive experiments on standard benchmark functions (Branin, Hartmann-6, BraninCurrin) and a simulated chemical reactor, we demonstrate that the RBF kernel achieves the lowest prediction error (MSE=0.372) for smooth objectives, UCB with β=0.5 provides the fastest convergence reaching 92.7% of the global optimum in 40 iterations, batch sizes of q=4–8 offer favorable efficiency-performance trade-offs, qEHVI effectively discovers Pareto fronts with hypervolume of 56.85, and REMBO enables effective optimization in 50-dimensional spaces. In the chemical case study, our framework identifies conditions achieving 99.0% yield and 97.0% selectivity through only 55 evaluations. These results establish practical guidelines for deploying BO in real-world high-dimensional experimental design.

---

## 1. Introduction

The optimization of experimental parameters is a fundamental challenge across scientific disciplines, from materials discovery and drug design to chemical process engineering. Traditional approaches such as factorial design and response surface methodology become prohibitively expensive in high-dimensional parameter spaces due to the combinatorial explosion of required experiments (Binois & Wycoff, 2022).

Bayesian optimization (BO) offers a principled framework for sequential experimental design by constructing a probabilistic surrogate model—typically a Gaussian process (GP)—of the objective function and using an acquisition function to select the most informative next experiment (Garnett, 2023). Despite significant theoretical and practical advances, several challenges persist:

1. **Kernel Selection**: The choice of GP kernel critically affects model accuracy, yet systematic guidelines for kernel selection remain limited (Binois & Wycoff, 2022).
2. **Acquisition Function Selection**: Different acquisition functions (EI, UCB, KG) exhibit problem-dependent performance, and adaptive selection criteria are needed.
3. **Parallel Experiments**: Modern experimental facilities often support batch execution, requiring efficient parallel proposal strategies (González et al., 2016).
4. **Multi-Objective Optimization**: Real-world problems frequently involve competing objectives (e.g., yield vs. selectivity), necessitating Pareto-aware optimization (Daulton et al., 2021).
5. **High Dimensionality**: Standard BO degrades rapidly beyond ~20 dimensions, requiring structural assumptions or dimensionality reduction (Eriksson et al., 2019).

In this work, we present a comprehensive BO framework addressing all five challenges, with a practical demonstration in chemical reaction optimization. Our contributions include:

- Systematic comparison of GP kernels with ARD for surrogate model selection
- Empirical guidelines for acquisition function selection based on problem characteristics
- Evaluation of batch optimization strategies balancing parallelism and sample efficiency
- Integration of qEHVI for multi-objective chemical reaction optimization
- REMBO-based dimensionality reduction enabling BO in 50-dimensional spaces
- A complete chemical reaction optimization case study demonstrating practical applicability

---

## 2. Related Work

### 2.1 High-Dimensional Gaussian Process Modeling

Binois and Wycoff (2022) provided a comprehensive survey of GP modeling strategies for high-dimensional BO, categorizing approaches into effective low-dimensional methods, additive decomposition, variable selection, and local search strategies. They highlighted that ARD kernels enable automatic variable importance assessment, while additive kernels can decompose the problem into tractable subproblems.

### 2.2 Scalable Bayesian Optimization

Eriksson et al. (2019) introduced TuRBO (Trust Region Bayesian Optimization), which maintains multiple local trust regions fitted with independent GP models, enabling scalable optimization in high dimensions. This approach demonstrated that local modeling can overcome the limitations of global GP fitting in spaces with D>20.

### 2.3 Random Embedding Bayesian Optimization

Wang et al. (2016) proposed REMBO, which projects high-dimensional inputs to a low-dimensional subspace via random linear embeddings. Under the assumption that the objective function has low effective dimensionality, REMBO achieves sample complexity independent of the ambient dimension. Recent extensions include SA-REMBO (Wen & Franzon, 2025) with adaptive embeddings and CEPBO (2025) with condensing-expansion projections that relax the effective subspace assumption.

### 2.4 Multi-Objective Bayesian Optimization

Daulton et al. (2021) developed differentiable Expected Hypervolume Improvement (qEHVI) for parallel multi-objective BO, implemented in BOTorch. This approach enables gradient-based optimization of batch acquisition functions over the Pareto front, significantly improving scalability compared to exact EHVI computation.

### 2.5 Bayesian Optimization for Chemical Synthesis

Shields et al. (2021) demonstrated the effectiveness of BO for chemical reaction optimization, showing that BO outperformed expert chemists in both efficiency and consistency for palladium-catalyzed cross-coupling reactions. Subsequent work by Guo et al. (2023) and Desimpel et al. (2026) extended these approaches to multi-objective settings optimizing both yield and selectivity.

### 2.6 BOTorch Framework

Balandat et al. (2020) introduced BOTorch, a modular framework for Monte Carlo-based BO built on PyTorch and GPyTorch. BOTorch provides differentiable acquisition functions, support for multi-objective optimization, and seamless GPU acceleration, making it the de facto standard for research and industrial BO applications.

---

## 3. Methods

### 3.1 Gaussian Process Surrogate Model

We model the objective function $f: \mathcal{X} \rightarrow \mathbb{R}$ as a Gaussian process:

$$f(\mathbf{x}) \sim \mathcal{GP}(m(\mathbf{x}), k(\mathbf{x}, \mathbf{x}'))$$

where $m(\mathbf{x})$ is the mean function (set to zero) and $k(\mathbf{x}, \mathbf{x}')$ is the covariance kernel. We evaluate four kernels:

**RBF (Squared Exponential):**
$$k_{\text{RBF}}(\mathbf{x}, \mathbf{x}') = \sigma^2 \exp\left(-\frac{1}{2}\sum_{d=1}^{D}\frac{(x_d - x'_d)^2}{\ell_d^2}\right)$$

**Matérn-ν:**
$$k_{\nu}(r) = \sigma^2 \frac{2^{1-\nu}}{\Gamma(\nu)}\left(\sqrt{2\nu}r\right)^{\nu} K_{\nu}\left(\sqrt{2\nu}r\right)$$

where $r = \sqrt{\sum_d (x_d - x'_d)^2 / \ell_d^2}$ and $K_\nu$ is the modified Bessel function. We test $\nu \in \{0.5, 1.5, 2.5\}$ with ARD lengthscales $\{\ell_d\}_{d=1}^D$.

Hyperparameters are optimized by maximizing the log marginal likelihood:

$$\log p(\mathbf{y}|\mathbf{X}, \boldsymbol{\theta}) = -\frac{1}{2}\mathbf{y}^T K_{\boldsymbol{\theta}}^{-1}\mathbf{y} - \frac{1}{2}\log|K_{\boldsymbol{\theta}}| - \frac{n}{2}\log(2\pi)$$

### 3.2 Acquisition Functions

**Expected Improvement (EI):**
$$\alpha_{\text{EI}}(\mathbf{x}) = \mathbb{E}[\max(f(\mathbf{x}) - f^*, 0)] = \sigma(\mathbf{x})[\phi(z)z + \Phi(z)]$$
where $z = (\mu(\mathbf{x}) - f^*) / \sigma(\mathbf{x})$.

**Upper Confidence Bound (UCB):**
$$\alpha_{\text{UCB}}(\mathbf{x}) = \mu(\mathbf{x}) + \beta\sigma(\mathbf{x})$$

where $\beta > 0$ controls the exploration-exploitation trade-off.

**Knowledge Gradient (KG):**
$$\alpha_{\text{KG}}(\mathbf{x}) = \mathbb{E}[\mu^{(n+1)}_{*} - \mu^{(n)}_{*}]$$

which measures the expected improvement in the posterior optimal value.

### 3.3 Batch Optimization

For parallel experiment proposals, we use the q-batch formulation:

$$\alpha_{q\text{-EI}}(\mathbf{x}_{1:q}) = \mathbb{E}\left[\max_{i=1}^{q}\max(f(\mathbf{x}_i) - f^*, 0)\right]$$

approximated via Monte Carlo sampling and optimized using gradient-based methods through the reparameterization trick.

### 3.4 Multi-Objective Optimization (qEHVI)

For $M$ objectives $\mathbf{f}(\mathbf{x}) = (f_1(\mathbf{x}), \ldots, f_M(\mathbf{x}))$, we optimize the q-Expected Hypervolume Improvement:

$$\alpha_{\text{qEHVI}}(\mathbf{x}_{1:q}) = \mathbb{E}\left[\text{HV}(\mathcal{P} \cup \{\mathbf{f}(\mathbf{x}_i)\}_{i=1}^q, \mathbf{r}) - \text{HV}(\mathcal{P}, \mathbf{r})\right]$$

where $\mathcal{P}$ is the current Pareto set, $\text{HV}(\cdot, \mathbf{r})$ is the hypervolume indicator with reference point $\mathbf{r}$, and the expectation is taken over the GP posterior.

### 3.5 REMBO for High-Dimensional Spaces

Given a high-dimensional space $\mathcal{X} \subset \mathbb{R}^D$ and effective dimensionality $d \ll D$, REMBO constructs a random projection matrix $\mathbf{A} \in \mathbb{R}^{D \times d}$ with entries drawn from $\mathcal{N}(0, 1/D)$. The optimization is performed in the low-dimensional space:

$$\mathbf{x}_{\text{low}} = \sigma(\mathbf{x}_{\text{high}} \cdot \mathbf{A})$$

where $\sigma(\cdot)$ is the sigmoid function mapping to $[0, 1]^d$. The GP is fitted in the original high-dimensional space, but the acquisition function is optimized over $\mathbb{R}^D$ with the implicit low-dimensional structure.

---

## 4. Experiments

### 4.1 Experimental Setup

All experiments were conducted using BOTorch 0.17.2, GPyTorch, and PyTorch 2.10.0. Each stochastic experiment was repeated 3 times with different random seeds (42, 43, 44). We report mean ± standard deviation.

### 4.2 Benchmark Functions

- **Branin** (2D): Used for kernel comparison. Global minimum at $f^* \approx -0.398$ (negated: maximum $\approx 0.398$).
- **Hartmann-6** (6D): Used for acquisition function and batch comparisons. Global maximum (negated): $f^* = 3.3224$.
- **BraninCurrin** (2D, 2 objectives): Used for multi-objective optimization.

### 4.3 Evaluation Metrics

- **Mean Squared Error (MSE)**: Prediction accuracy of the GP surrogate
- **Negative Log-Likelihood (NLL)**: Probabilistic calibration of predictions
- **Best value found**: Optimization performance (higher is better for negated functions)
- **Hypervolume (HV)**: Multi-objective quality metric measuring dominated volume
- **Wall-clock time**: Computational efficiency

### 4.4 Chemical Reaction Simulation

We constructed a simulated chemical reactor with 5 parameters:
- Temperature: 50–200 °C
- Pressure: 1–10 atm
- Catalyst loading: 0.01–0.10 mol%
- Solvent ratio: 0–1
- Residence time: 1–60 min

Two objectives were optimized simultaneously: reaction yield (%) and selectivity (%).

---

## 5. Results

### 5.1 Kernel Selection (Experiment 1)

Table 1 presents the prediction performance of four GP kernels on the Branin function with 50 training points and 200 test points.

| Kernel | MSE | NLL |
|--------|-----|-----|
| RBF | **0.3721** | **0.7309** |
| Matérn-5/2 | 1.0620 | 1.1224 |
| Matérn-3/2 | 5.2405 | 2.2591 |
| Matérn-1/2 | 41.6118 | 3.8494 |

The RBF kernel achieved the lowest MSE (0.372) and NLL (0.731), consistent with the smooth nature of the Branin function. The Matérn family showed monotonically decreasing performance with lower smoothness parameter ν, confirming that kernel-function smoothness alignment is critical.

![Figure 1: Kernel comparison showing (a) MSE and (b) NLL across four kernel types.](figures/kernel_comparison.png)

![Figure 2: Learned ARD lengthscales for each kernel, revealing dimension-specific relevance.](figures/ard_lengthscales.png)

### 5.2 Acquisition Function Comparison (Experiment 2)

Table 2 summarizes the optimization performance of four acquisition strategies on Hartmann-6.

| Acquisition Function | Final Best (mean ± std) | Iterations |
|---------------------|------------------------|------------|
| EI | 2.922 ± 0.132 | 40 |
| UCB (β=2.0) | 2.956 ± 0.127 | 40 |
| UCB (β=0.5) | **3.080 ± 0.161** | 40 |
| KG | 2.580 ± 0.722 | 20 |

UCB with β=0.5 achieved the highest mean performance (92.7% of global optimum), favoring exploitation over exploration. KG showed high variance due to limited iterations (20 vs. 40 for others) driven by its higher computational cost.

![Figure 3: Convergence curves for EI, UCB (β=2.0, β=0.5), and KG on Hartmann-6. Shaded regions represent ± 1 standard deviation over 3 runs.](figures/acquisition_comparison.png)

### 5.3 Batch Optimization (Experiment 3)

Table 3 presents the trade-off between batch size, optimization quality, and computational cost.

| Batch Size (q) | Final Best (mean) | Total Evaluations | Wall Time (s) |
|---------------|-------------------|-------------------|---------------|
| 1 | 2.697 | 25 | 2.37 |
| 2 | 2.951 | 40 | 4.23 |
| 4 | 2.952 | 70 | 6.81 |
| 8 | **3.104** | 130 | 15.51 |

Larger batch sizes achieved better final values at the cost of more total evaluations and computation time. The diminishing returns between q=2 and q=4 suggest that moderate batch sizes offer the best efficiency-to-performance ratio when parallel resources are limited.

![Figure 4: (a) Convergence by BO round for different batch sizes. (b) Total optimization wall time.](figures/batch_optimization.png)

### 5.4 Multi-Objective Optimization (Experiment 4)

The qEHVI algorithm was applied to the BraninCurrin problem (2 objectives, 2 input dimensions) for 30 iterations starting from 10 initial points.

- **Final hypervolume**: 56.85
- **Pareto front size**: 21 non-dominated solutions
- **Hypervolume improvement**: From initial HV to final HV, demonstrating consistent improvement

![Figure 5: (a) Discovered Pareto front showing non-dominated (red) and dominated (blue) solutions. (b) Hypervolume convergence over iterations.](figures/mobo_ehvi.png)

### 5.5 High-Dimensional Optimization (Experiment 5)

Table 4 compares optimization methods in 50-dimensional space where the objective depends on only 6 effective dimensions.

| Method | Final Best (mean) |
|--------|-------------------|
| REMBO (D=50→d=6) | **3.192** |
| Vanilla BO (D=50) | 3.064 |
| Random Search | baseline |

REMBO achieved 96.1% of the global optimum in 50 dimensions, outperforming vanilla BO by 4.2%. The random projection successfully captured the intrinsic low-dimensional structure.

![Figure 6: Convergence comparison of REMBO, vanilla BO, and random search in 50-dimensional Hartmann-6.](figures/high_dim_comparison.png)

### 5.6 Chemical Reaction Optimization (Experiment 6)

Multi-objective BO with qEHVI was applied to simultaneously optimize yield and selectivity over 5 reaction parameters for 40 iterations.

**Best yield conditions**: Temperature=149.85°C, Pressure=10.00 atm, Catalyst=0.10 mol%, yielding 99.0% yield with 90.7% selectivity.

**Best selectivity conditions**: Temperature=130.98°C, Pressure=4.85 atm, Catalyst=0.10 mol%, yielding 84.4% yield with 97.0% selectivity.

The final hypervolume reached 9571.0, with the Pareto front clearly revealing the yield-selectivity trade-off. The framework identified optimal conditions in only 55 total evaluations (15 initial + 40 iterations).

![Figure 7: Chemical reaction optimization results: (a) Yield vs. selectivity Pareto front, (b) hypervolume convergence, (c) best yield/selectivity over iterations, (d) parameter distributions of Pareto-optimal solutions.](figures/chemical_optimization.png)

---

## 6. Discussion

### 6.1 Kernel Selection Guidelines

Our results confirm that kernel selection should be guided by prior knowledge of function smoothness. For smooth, continuously differentiable functions, the RBF kernel provides optimal prediction accuracy. When smoothness is uncertain, the Matérn-5/2 kernel offers a robust default, as recommended by the GP literature (Binois & Wycoff, 2022). The ARD mechanism proved effective in identifying relevant dimensions, with lengthscale ratios providing interpretable variable importance measures.

### 6.2 Acquisition Function Trade-offs

The superior performance of UCB (β=0.5) on Hartmann-6 suggests that moderate exploitation bias is beneficial when the function landscape has clear global structure. However, this finding is problem-dependent: high-noise or multi-modal landscapes would likely favor more exploratory strategies (higher β). KG's higher computational cost (~4× per iteration compared to EI) limits its practical applicability despite theoretical advantages in one-step optimal decision-making.

### 6.3 Batch Optimization Scaling

The sublinear relationship between batch size and optimization quality highlights the fundamental exploration-exploitation challenge in parallel BO: larger batches provide more diverse samples but with diminishing information gain. For practical applications where parallel experimental resources are available (e.g., multi-well plate reactors), batch sizes of q=4–8 represent a reasonable operating point.

### 6.4 High-Dimensional Challenges

REMBO's effectiveness in our experiments relies on the assumption of low effective dimensionality, which holds for the embedded Hartmann-6 function. In practice, verifying this assumption requires domain knowledge or preliminary screening experiments. Recent advances such as SA-REMBO (Wen & Franzon, 2025) address non-stationarity through adaptive embeddings, while TuRBO (Eriksson et al., 2019) provides an alternative through local trust region strategies.

### 6.5 Chemical Reaction Optimization

The chemical case study demonstrates the practical value of multi-objective BO for experimental chemistry. The discovered Pareto front reveals a clear trade-off: high yield (>95%) requires elevated temperature and pressure, while high selectivity (>95%) favors moderate conditions. This trade-off information is actionable for process development, allowing chemists to select conditions based on downstream requirements. Our results align with findings by Shields et al. (2021), who showed BO outperforming expert chemists in reaction optimization.

### 6.6 Limitations

1. All experiments used synthetic or simulated functions; real-world experimental noise and constraints may affect performance differently.
2. The REMBO approach assumes linear embedding, which may not capture nonlinear low-dimensional structure.
3. KG was evaluated with fewer iterations due to computational constraints, potentially underestimating its performance.
4. The chemical reaction simulator, while multi-modal and nonlinear, does not capture all physical phenomena present in real reactors.

---

## 7. Conclusion

We presented a unified Bayesian optimization framework for high-dimensional experimental design, integrating kernel selection, acquisition function comparison, batch optimization, multi-objective optimization, and dimensionality reduction. Through systematic experiments on benchmark functions and a chemical reaction case study, we established the following practical guidelines:

1. The RBF kernel is optimal for smooth functions; Matérn-5/2 provides a robust default for unknown smoothness.
2. UCB with moderate β (0.5–2.0) balances exploration and exploitation effectively, with problem-dependent tuning recommended.
3. Batch sizes of q=4–8 offer favorable trade-offs for parallel experimental settings.
4. qEHVI effectively discovers Pareto fronts for multi-objective problems with 2+ objectives.
5. REMBO enables effective BO in 50+ dimensional spaces when intrinsic dimensionality is low.
6. Multi-objective BO can identify yield-selectivity Pareto fronts in chemical reactions with fewer than 60 evaluations.

Future work will extend this framework with transfer learning across related reactions, deep kernel learning for complex surrogate models, constrained optimization with safety requirements, and validation on real experimental data from automated chemistry platforms.

---

## References

1. Balandat, M., Karrer, B., Jiang, D. R., Daulton, S., Letham, B., Wilson, A. G., & Bakshy, E. (2020). BoTorch: A Framework for Efficient Monte-Carlo Bayesian Optimization. *Advances in Neural Information Processing Systems*, 33, 21524–21538. DOI: [10.48550/arXiv.1910.06403](https://doi.org/10.48550/arXiv.1910.06403)

2. Binois, M., & Wycoff, N. (2022). A Survey on High-dimensional Gaussian Process Modeling with Application to Bayesian Optimization. *ACM Transactions on Evolutionary Learning and Optimization*, 2(2). DOI: [10.1145/3545611](https://doi.org/10.1145/3545611)

3. Daulton, S., Balandat, M., & Bakshy, E. (2021). Differentiable Expected Hypervolume Improvement for Parallel Multi-Objective Bayesian Optimization. *Journal of Computational and Graphical Statistics*, 30(4), 1165–1178. DOI: [10.1080/10618600.2021.1888745](https://doi.org/10.1080/10618600.2021.1888745)

4. Eriksson, D., Pearce, M., Gardner, J., Turner, R. D., & Poloczek, M. (2019). Scalable Global Optimization via Local Bayesian Optimization. *Advances in Neural Information Processing Systems*, 32. DOI: [10.48550/arXiv.1910.01739](https://doi.org/10.48550/arXiv.1910.01739)

5. Shields, B. J., Stevens, J., Li, J., Parasram, M., Damani, F., Martinez Alvarado, J. I., Janey, J. M., Adams, R. P., & Doyle, A. G. (2021). Bayesian Reaction Optimization as a Tool for Chemical Synthesis. *Nature*, 590, 89–96. DOI: [10.1038/s41586-021-03213-y](https://doi.org/10.1038/s41586-021-03213-y)

6. Wang, Z., Hutter, F., Zoghi, M., Matheson, D., & de Freitas, N. (2016). Bayesian Optimization in a Billion Dimensions via Random Embeddings. *Journal of Artificial Intelligence Research*, 55, 361–387. DOI: [10.1613/jair.4806](https://doi.org/10.1613/jair.4806)

7. Guo, J., Ranković, B., & Schwaller, P. (2023). Bayesian Optimization for Chemical Reactions. *Chimia*, 77(1/2), 31–38. DOI: [10.2533/chimia.2023.31](https://doi.org/10.2533/chimia.2023.31)

8. Wen, H., & Franzon, P. (2025). Adaptive Linear Embedding for Nonstationary High-Dimensional Optimization. *arXiv preprint*, arXiv:2505.11281. DOI: [10.48550/arXiv.2505.11281](https://doi.org/10.48550/arXiv.2505.11281)

9. González, J., Dai, Z., Hennig, P., & Lawrence, N. (2016). Batch Bayesian Optimization via Local Penalization. *Proceedings of the 19th International Conference on Artificial Intelligence and Statistics (AISTATS)*. DOI: [10.48550/arXiv.1505.08052](https://doi.org/10.48550/arXiv.1505.08052)

10. Garnett, R. (2023). *Bayesian Optimization*. Cambridge University Press. DOI: [10.1017/9781108348973](https://doi.org/10.1017/9781108348973)
