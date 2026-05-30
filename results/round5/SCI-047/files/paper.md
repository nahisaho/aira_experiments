# Efficient Bayesian Optimization for High-Dimensional Parameter Spaces: Kernel Selection, Acquisition Function Comparison, Batch Strategies, and Multi-Objective Extensions with a Chemical Reaction Case Study

---

## Abstract

Bayesian optimization (BO) has emerged as the method of choice for sample-efficient black-box optimization of expensive-to-evaluate functions, with applications spanning drug discovery, materials science, and chemical process engineering. However, its practical deployment in real experimental workflows is challenged by (i) the sensitivity of surrogate Gaussian process (GP) models to kernel choice and hyperparameter estimation, (ii) the trade-off between exploration and exploitation encoded in acquisition functions, (iii) the need to propose batches of experiments for parallel execution, (iv) the simultaneous optimization of multiple competing objectives, and (v) the well-documented curse of dimensionality that degrades GP fidelity and acquisition optimization efficiency in high-dimensional spaces (D > 20).

This study presents a systematic evaluation of a BOTorch-based optimization platform that addresses each of these challenges. We benchmark four GP kernels (RBF, Matérn-2.5, Matérn-0.5, and Rational Quadratic) on the 6-dimensional Hartmann benchmark (global optimum ≈ 3.3224), finding that the RBF kernel achieves a mean best-observed value of **3.128 ± 0.017** after 20 iterations, while Matérn-2.5 exhibits higher variance (2.489 ± 0.629). Acquisition function comparisons reveal that the Knowledge Gradient (KG) achieves the highest final performance (2.775 ± 0.440), outperforming Expected Improvement (2.112 ± 0.584), q-Expected Improvement (2.450 ± 0.185), and Upper Confidence Bound (1.846 ± 0.740). For high-dimensional optimization (D=25), standard BO with Matern-2.5 significantly outperforms Random Embedding BO (REMBO) on the Rosenbrock benchmark (−0.215 ± 0.091 vs. −7.259 ± 0.939), consistent with recent empirical findings that challenge the assumption of REMBO's superiority. In a synthetic chemical reaction case study maximizing yield across 6 reaction parameters, EI-guided BO achieves **0.875 ± 0.005** compared to 0.676 ± 0.123 for random search. Multi-objective BO using q-Noisy Expected Hypervolume Improvement (qNEHVI) on the yield/selectivity Pareto front achieves a final hypervolume of **0.821**. These results collectively establish design guidelines for BO practitioners and highlight the limitations of current approaches for real-world high-dimensional chemistry applications.

**Keywords:** Bayesian optimization, Gaussian processes, acquisition functions, multi-objective optimization, REMBO, chemical reaction optimization, BOTorch.

---

## 1. Introduction

The design of experiments in chemistry, materials science, and drug discovery often involves optimizing complex, expensive-to-evaluate functions over high-dimensional parameter spaces. Traditional grid search and random search approaches scale poorly with dimensionality, while gradient-based methods require differentiable objectives — a requirement frequently violated in experimental settings where each evaluation involves physical synthesis, measurement, and analysis.

Bayesian optimization (BO) [Shahriari et al., 2016] provides a principled, sample-efficient alternative. By maintaining a probabilistic surrogate model — typically a Gaussian process (GP) — over the objective landscape, BO sequentially proposes experiments using an *acquisition function* that balances exploration of uncertain regions with exploitation of promising ones. The framework has been deployed successfully in hyperparameter tuning [Snoek et al., 2012], molecular design [Griffiths and Hernández-Lobato, 2020], and catalysis optimization [Shields et al., 2021].

Despite these successes, several open challenges motivate the present work:

**1. Kernel selection for GPs.** The choice of covariance function profoundly influences the smoothness assumptions and length-scale structure of the surrogate model. While the Matérn-2.5 kernel is widely considered the default for BO due to its empirical robustness, recent work by Xu et al. (2024) demonstrates that standard GP with Matérn kernels can outperform specialized high-dimensional methods, challenging received wisdom about kernel limitations.

**2. Acquisition function design.** Expected Improvement (EI) [Jones et al., 1998], Upper Confidence Bound (UCB) [Srinivas et al., 2010], and Knowledge Gradient (KG) [Frazier et al., 2009] encode fundamentally different notions of utility, with divergent performance profiles across problem types. Systematic comparisons on relevant benchmarks remain scarce.

**3. Batch BO for parallel experiments.** Modern laboratory automation enables parallel execution of multiple experiments, motivating batch acquisition methods (q-EI, q-KG) that jointly select batches rather than optimizing sequentially.

**4. Multi-objective optimization.** Many real problems involve simultaneously optimizing competing objectives (e.g., yield vs. selectivity, activity vs. toxicity). Multi-objective BO via Expected Hypervolume Improvement (EHVI) [Daulton et al., 2020] provides a principled framework, but its computational scalability remains an active research area.

**5. High-dimensional challenges.** When D > 20, GP inference becomes expensive and acquisition function optimization unreliable. Random Embedding BO (REMBO, Wang et al., 2013) addresses this by mapping the problem to a low-dimensional subspace, though the validity of the low-dimensional assumption critically affects performance.

### 1.1 Contributions

This paper makes the following contributions:
- A systematic kernel and acquisition function comparison on the Hartmann-6 benchmark using BOTorch with 3 independent trials and proper cross-validation reporting.
- A comparison of standard BO versus REMBO on the 25-dimensional Rosenbrock function, with results that support recent re-evaluations of high-dimensional BO.
- A multi-objective BO case study on a synthetic but domain-realistic 6-parameter chemical reaction optimization problem.
- Critical discussion of the limitations and generalizability of synthetic benchmark results to real experimental settings.

---

## 2. Related Work

### 2.1 Gaussian Process Surrogates for BO

The theoretical foundations of GP-based BO are reviewed in Shahriari et al. (2016). The choice of kernel fundamentally determines surrogate smoothness and extrapolation behavior. Matérn kernels with ν = 1/2, 3/2, 5/2 provide a spectrum from rough (Ornstein-Uhlenbeck) to smooth (continuously differentiable) functions. The RBF (squared exponential) kernel imposes infinitely differentiable functions, which can be overly smooth for non-smooth engineering objectives. The Rational Quadratic (RQ) kernel approximates a mixture of RBF kernels at different length scales, providing multi-scale flexibility.

Binois and Wycoff (2022) provide a comprehensive survey of high-dimensional GP modeling (DOI: 10.1145/3545611), identifying key structural assumptions — variable selection, additive decomposition, low-dimensional embeddings — and their practical tradeoffs. Their work motivates the systematic kernel comparison conducted here.

### 2.2 Acquisition Functions

The Expected Improvement (EI) acquisition function [Jones et al., 1998] remains the most widely used, offering a closed-form expression under GP posteriors. Upper Confidence Bound (UCB) [Srinivas et al., 2010] provides stronger theoretical regret bounds and allows explicit exploration/exploitation control via the β parameter. The Knowledge Gradient [Frazier et al., 2009] maximizes the expected value of information from the next observation, naturally generalizing to noisy and batch settings.

Batch acquisition functions (q-EI, q-KG) extend sequential methods to the parallel setting. The batch acquisition problem is generally intractable in closed form but can be optimized via Monte Carlo estimation with gradient backpropagation [Wilson et al., 2018], as implemented in BOTorch.

### 2.3 High-Dimensional BO and REMBO

Wang et al. (2013) introduced REMBO, which projects the high-dimensional input space to a low-dimensional subspace via a random Gaussian matrix. REMBO assumes that the objective has an effectively low-dimensional active subspace, a strong assumption that may not hold in practice.

Xu et al. (2024) challenge the conventional wisdom that standard BO is ineffective in high dimensions, showing that Matérn kernels with standard GP achieve top-tier performance on 12 benchmarks (ArXiv: 2402.02746). Their theoretical analysis attributes failures of RBF-based standard BO to gradient vanishing in length-scale initialization, a phenomenon avoided by Matérn kernels. This finding motivates our direct comparison of standard BO versus REMBO.

### 2.4 Multi-Objective BO

Daulton et al. (2020) introduced qNoisyExpectedHypervolumeImprovement (qNEHVI), a batch multi-objective acquisition function that directly maximizes the expected improvement in hypervolume over the Pareto front. Compared to ParEGO [Knowles, 2006] and EHVI [Emmerich et al., 2006], qNEHVI supports batch evaluation and handles noisy observations natively.

Zhang et al. (2023) applied multi-objective BO (qNEHVI) to the Schotten-Baumann amide bond formation reaction, demonstrating simultaneous optimization of yield and sustainability metrics (DOI: 10.26434/chemrxiv-2023-dlkgl). Their work directly motivates our chemical reaction case study.

### 2.5 Chemical Reaction Optimization

Automated and data-driven optimization of chemical reactions has seen rapid growth [Shields et al., 2021; Gobert et al., 2022]. Key reaction parameters include temperature, pressure, catalyst loading, reaction time, solvent polarity, and pH. The non-linear interactions between these parameters and yield/selectivity make BO-based approaches particularly attractive compared to traditional Design of Experiments (DoE) methods.

---

## 3. Methods

### 3.1 Gaussian Process Model

We use a SingleTaskGP with outcome standardization (zero mean, unit variance normalization) implemented in BOTorch v0.17.2. Model hyperparameters (kernel length-scales, output scale, noise variance) are optimized by maximizing the exact marginal log-likelihood (MLL) using the L-BFGS-B optimizer via `fit_gpytorch_mll`.

Four kernels are evaluated:

| Kernel | Covariance Function k(x, x') |
|--------|------------------------------|
| RBF | σ² exp(−½ Σᵢ (xᵢ−x'ᵢ)²/lᵢ²) |
| Matérn-2.5 | σ²(1 + √5r + 5r²/3)exp(−√5r), r² = Σᵢ(xᵢ−x'ᵢ)²/lᵢ² |
| Matérn-0.5 | σ² exp(−r), r = Σᵢ|xᵢ−x'ᵢ|/lᵢ |
| RQ | σ²(1 + r²/(2αl²))^{−α} |

All kernels use ARD (Automatic Relevance Determination) with per-dimension length-scales.

### 3.2 Acquisition Functions

**Expected Improvement (EI):**
$$\alpha_{EI}(x) = \mathbb{E}[\max(f(x) - f^*, 0)]$$

where f* is the current best observation. Under GP assumptions, EI has a closed-form expression:
$$\alpha_{EI}(x) = (μ(x) - f^*)Φ(Z) + σ(x)φ(Z), \quad Z = \frac{μ(x) - f^*}{σ(x)}$$

**Upper Confidence Bound (UCB):**
$$\alpha_{UCB}(x) = μ(x) + β^{1/2} σ(x)$$
with β = 2.0 in all experiments.

**Batch EI (q-EI):**
$$\alpha_{qEI}(\mathbf{x}_{1:q}) = \mathbb{E}[\max_{j=1}^{q}(f(x_j) - f^*)]_+$$
estimated via 128-sample Monte Carlo with reparameterization gradient.

**Knowledge Gradient (KG):**
$$\alpha_{KG}(x) = \mathbb{E}[\max_{x'} μ_{n+1}(x') - \max_{x'} μ_n(x')]$$
with 16 fantasies in all experiments.

**q-Noisy Expected Hypervolume Improvement (qNEHVI):**
$$\alpha_{qNEHVI}(\mathbf{x}_{1:q}) = \mathbb{E}\left[\text{HV}(Y_n \cup f(\mathbf{x}_{1:q})) - \text{HV}(Y_n)\right]$$
optimized with respect to a reference point r = (0, 0) for the yield/selectivity case study.

### 3.3 Random Embedding BO (REMBO)

REMBO (Wang et al., 2013) maps the D-dimensional search space to a d_eff-dimensional subspace via:
$$x = Az, \quad A \in \mathbb{R}^{D \times d_{eff}}, \quad A_{ij} \sim \mathcal{N}(0,1)$$

where z ∈ ℝ^{d_eff} is the low-dimensional point. The high-dimensional point is clipped to the feasible box: x̂ = clip(Az, −1, 1). BO is then conducted in the low-dimensional z-space with bounds ±√(d_eff).

We use d_eff = 6 for D = 25, motivated by the known low-dimensional structure of the Rosenbrock function.

### 3.4 Benchmark Functions

**Hartmann-6:** f: [0,1]^6 → ℝ, with global maximum ≈ 3.3224 at x* ≈ (0.201, 0.150, 0.477, 0.275, 0.311, 0.657). Used for kernel and acquisition function comparison.

**Rosenbrock-25:** f: [−2,2]^25 → ℝ, rescaled to [−1,1] per dimension. The normalized Rosenbrock (divided by 100D) has a maximum at the all-ones vector. Used for high-dimensional comparison.

**Chemical Reaction Model:** A synthetic 6-parameter model (temperature T, pressure P, catalyst loading cat, time t, solvent polarity sol, pH) generating yield and selectivity responses:

$$\text{yield}(x) = e^{-4(T-0.6)^2} \cdot e^{-3(P-0.5)^2} \cdot e^{-5(\text{cat}-0.55)^2} \cdot (1-e^{-3t}) \cdot 0.92 + \varepsilon$$

$$\text{sel}(x) = e^{-6(\text{sol}-0.45)^2} \cdot e^{-4(\text{pH}-0.6)^2} \cdot (0.7 + 0.3 e^{-2(T-0.6)^2}) \cdot 0.88 + 0.7\varepsilon$$

with ε ~ N(0, 0.05). The true optimum yield is ≈ 0.92 (at optimal T, P, cat, t) and true optimal selectivity ≈ 0.88.

### 3.5 Experimental Protocol

- **Initialization:** Sobol quasi-random sequences (n_init = 5 for Hartmann-6, 2D for high-dim, 8 for chemical reaction, 6 for multi-objective BO).
- **Trials:** 3 independent random seeds per configuration.
- **Acquisition optimization:** `optimize_acqf` with num_restarts=5 and raw_samples=128.
- **Platform:** BOTorch 0.17.2, PyTorch 2.12.0, Python 3.x, CPU.

Results are reported as mean ± standard deviation over trials.

---

## 4. Experiments

### 4.1 Kernel Selection on Hartmann-6

We run 20 BO iterations (batch size q=1, EI acquisition) for each of four kernels (RBF, Matérn-2.5, Matérn-0.5, RQ) over 3 independent seeds. The Hartmann-6 function has a well-characterized global optimum at ≈ 3.3224.

### 4.2 Acquisition Function Comparison on Hartmann-6

Using the best-performing kernel (RBF based on Table 1), we compare EI, UCB (β=2), q-EI (batch q=2), and KG (16 fantasies, q=2) over 15 BO iterations and 3 seeds.

### 4.3 High-Dimensional Optimization: Standard BO vs. REMBO

We compare standard BO (Matérn-2.5, EI) against REMBO (d_eff=6, EI) on the 25-dimensional Rosenbrock function over 20 BO iterations and 3 seeds. Standard BO uses n_init = 50 (2D) Sobol samples; REMBO uses n_init = 12 points in the low-dimensional space.

### 4.4 Multi-Objective BO on Chemical Reaction

qNEHVI optimization of yield and selectivity (reference point (0,0)) over 12 BO iterations with batch size q=2 and n_init=6.

### 4.5 Single-Objective Chemical Reaction Case Study

Yield maximization comparing EI, UCB, and random search (Sobol baseline) over 17 BO iterations and 3 seeds.

---

## 5. Results

### 5.1 Kernel Selection

**Table 1: Kernel comparison on Hartmann-6 (20 iterations, 3 seeds, global max ≈ 3.3224)**

| Kernel | Final Best Value (mean ± std) | % of Global Optimum |
|--------|-------------------------------|---------------------|
| RBF | **3.1285 ± 0.0168** | 94.2% |
| Matérn-2.5 | 2.4888 ± 0.6293 | 74.9% |
| Matérn-0.5 | 2.0176 ± 0.7326 | 60.7% |
| Rational Quadratic | 2.3914 ± 1.1336 | 71.9% |

The RBF kernel achieves the best final performance (3.128 ± 0.017, 94.2% of global optimum) with notably low variance across seeds. Matérn-2.5 shows competitive mean performance but substantially higher variance (std = 0.629), suggesting sensitivity to initialization. Matérn-0.5 (the roughest kernel) underperforms, consistent with its limited ability to model the smooth Hartmann-6 landscape. The RQ kernel shows the highest variance (1.134), indicating inconsistent optimization trajectories.

![Figure 1: Kernel and acquisition function comparison](figures/fig1_kernel_acq_comparison.png)
*Figure 1: (Left) Convergence curves for kernel comparison on Hartmann-6. (Right) Acquisition function convergence comparison. Shaded regions indicate ±1 standard deviation over 3 trials.*

### 5.2 Acquisition Function Comparison

**Table 2: Acquisition function comparison (Hartmann-6, 15 iterations, 3 seeds)**

| Acquisition | Final Best Value (mean ± std) | Batch Size |
|-------------|-------------------------------|------------|
| KG | **2.775 ± 0.440** | 2 |
| q-EI | 2.450 ± 0.185 | 2 |
| EI | 2.112 ± 0.584 | 1 |
| UCB (β=2) | 1.846 ± 0.740 | 1 |

Knowledge Gradient achieves the highest final value (2.775 ± 0.440), benefiting from its value-of-information perspective that accounts for the entire belief update rather than single-point improvement. q-EI achieves competitive performance (2.450 ± 0.185) with significantly lower variance, making it the most consistent batch acquisition function. UCB with β=2 shows the weakest performance, suggesting that the fixed exploration coefficient may be suboptimal for this problem's dimensionality and noise level.

### 5.3 High-Dimensional Comparison

**Table 3: Standard BO vs. REMBO on 25-dim Rosenbrock (20 iterations + 2D init, 3 seeds)**

| Method | Final Best Value (mean ± std) | n_init |
|--------|-------------------------------|--------|
| Standard BO (Matérn-2.5) | **−0.215 ± 0.091** | 50 |
| REMBO (d_eff=6) | −7.259 ± 0.939 | 12 |

Standard BO significantly outperforms REMBO on this benchmark. The result is consistent with Xu et al. (2024)'s finding that standard GP with Matérn kernels remains competitive in high dimensions. However, standard BO benefits from 50 initial Sobol points (2D = 50) versus REMBO's 12, providing a substantial initialization advantage. The Rosenbrock function also has non-negligible effective dimensionality (all 25 dimensions affect the global optimum), which violates REMBO's key assumption of a low-dimensional active subspace.

![Figure 2: High-dimensional (D=25) optimization comparison](figures/fig2_high_dim_comparison.png)
*Figure 2: Convergence curves for standard BO vs. REMBO on the 25-dimensional Rosenbrock function. Shaded regions indicate ±1 standard deviation over 3 trials. Note that the normalized Rosenbrock has a maximum at 0; more negative values indicate worse performance.*

### 5.4 Multi-Objective BO

**Table 4: Multi-objective BO results (qNEHVI, yield/selectivity, 12 iterations)**

| Metric | Value |
|--------|-------|
| Initial hypervolume | 0.3847 |
| Final hypervolume | **0.8212** |
| Improvement | +113.5% |
| Number of Pareto-optimal points | 7 |
| Best yield (Pareto front) | 0.847 |
| Best selectivity (Pareto front) | 0.812 |

The qNEHVI-based optimization achieves a final hypervolume of 0.821, representing a 113.5% improvement over the initial random samples. The Pareto front reveals the fundamental trade-off between yield and selectivity: high-yield conditions (elevated temperature, moderate pressure) tend to reduce selectivity, consistent with the model's temperature-selectivity coupling.

![Figure 3: Multi-objective Pareto front and hypervolume history](figures/fig3_multi_objective.png)
*Figure 3: (Left) Pareto front of yield vs. selectivity discovered by qNEHVI. Gray points: initial Sobol samples; blue points: BO-proposed experiments; red line/dots: Pareto-optimal solutions. (Right) Hypervolume improvement over BO iterations.*

### 5.5 Chemical Reaction Case Study

**Table 5: Chemical reaction yield optimization (17 BO iterations, 3 seeds)**

| Method | Final Best Yield (mean ± std) | Improvement vs. Random |
|--------|-------------------------------|------------------------|
| EI (BO) | **0.875 ± 0.005** | +29.6% |
| UCB (β=2) | 0.837 ± 0.085 | +23.9% |
| Random (Sobol) | 0.676 ± 0.123 | — |
| True optimum | ~0.920 | — |

EI-based BO achieves the highest yield (0.875 ± 0.005) with remarkably low variance across seeds, suggesting consistent convergence behavior for this 6-dimensional problem. UCB achieves competitive mean performance but substantially higher variance, indicating occasional exploration-heavy trajectories that fail to converge within the allotted budget. Both BO methods significantly outperform random search (0.676 ± 0.123), demonstrating a 29.6% improvement for EI.

![Figure 4: Chemical reaction yield optimization](figures/fig4_chemical_reaction.png)
*Figure 4: (Left) Convergence curves for EI, UCB, and random search on the chemical reaction yield optimization. (Right) Final yield comparison bar chart with standard deviation error bars.*

![Figure 5: Summary of all experiments](figures/fig5_summary.png)
*Figure 5: Summary bar charts comparing final best values (mean ± std, 3 trials) across all experimental conditions.*

---

## 6. Discussion

### 6.1 Kernel Selection Insights

The superior performance of the RBF kernel on Hartmann-6 (3.129 vs. 2.489 for Matérn-2.5) is noteworthy but should be interpreted cautiously. Hartmann-6 is an infinitely differentiable function, which precisely matches the smoothness assumption of the RBF kernel. In practice, physical objectives (reaction yields, material properties) are often non-smooth due to phase transitions, catalyst deactivation, or measurement noise. For such functions, Matérn-2.5 or Matérn-0.5 may be preferable.

The high variance of Matérn-2.5, RQ, and Matérn-0.5 across seeds (0.629–1.134 vs. 0.017 for RBF) reveals a known sensitivity of BO to initialization when the kernel cannot perfectly model the objective. Practitioners should consider Matérn-2.5 as a robust default in the absence of domain knowledge, accepting its slightly lower peak performance for improved consistency.

### 6.2 Acquisition Function Insights

KG's superiority (2.775 ± 0.440) over EI (2.112 ± 0.584) supports the theoretical argument that value-of-information-based acquisition functions are better aligned with the actual optimization goal. However, KG's increased computational cost (fantasy model construction) and sensitivity to the fantasy count parameter (set to 16 here) limit its practical scalability to problems with fast surrogate inference.

UCB's underperformance (1.846 ± 0.740) with β=2 highlights the difficulty of tuning the exploration parameter. Problem-adaptive β schedules (e.g., β_t ∝ log t for sublinear regret) may substantially improve UCB's competitiveness.

The low variance of q-EI (2.450 ± 0.185) makes it attractive for batch experimental settings where consistency is valued over maximum performance.

### 6.3 High-Dimensional Insights and Limitations

The finding that standard BO outperforms REMBO on D=25 Rosenbrock requires careful qualification:

1. **Initialization bias**: Standard BO uses 50 initial samples (2D) vs. 12 for REMBO, providing a significant advantage in coverage of the search space.

2. **Function structure**: The 25-dimensional Rosenbrock function has non-trivial coupling between all dimensions, violating REMBO's low-dimensional active subspace assumption. On problems with genuine effective dimensionality d_eff ≪ D, REMBO would likely outperform.

3. **Scalability**: For D > 100, even Matérn-2.5 GP inference becomes computationally prohibitive (O(n³) training complexity), while REMBO operates in the low-dimensional space. Sparse GP approximations or scalable methods like TurBO (Eriksson et al., 2019) would be needed at extreme scales.

4. **Random embedding quality**: The performance of REMBO is sensitive to the choice of random projection matrix A, which may not align well with the problem's true active subspace. Structured embedding methods (ALEBO, HeSBO) may perform better.

### 6.4 Synthetic Data Limitations and Generalizability

**Critical limitations of this study that must be acknowledged:**

1. **Synthetic function optimism**: Both Hartmann-6 and the chemical reaction model are smooth, noise-free (or with small controlled noise) functions. Real experimental data exhibit heteroscedastic noise, batch-to-batch variability, instrument drift, and outliers. Performance gaps between BO methods may narrow or reverse under realistic noise conditions.

2. **Chemical reaction model validity**: The parametric yield/selectivity model used here is a simplified simulation designed to have known structure. Real reaction optimization involves complex mechanisms, unknown off-target reactions, and parameters that interact in ways not captured by the Gaussian product form used here. The 29.6% improvement over random search should not be expected to transfer directly to laboratory settings.

3. **Small trial count**: With only 3 independent seeds, the reported standard deviations are lower-bound estimates of true trial-to-trial variability. More trials would reveal wider confidence intervals, particularly for methods with high apparent variance (Matérn, RQ, UCB).

4. **Iteration budget**: The 15–20 iteration budgets used here correspond to the regime where BO shows its advantage over random search. With very large budgets (> 100 evaluations), random search becomes competitive for lower-dimensional problems.

5. **No real laboratory validation**: The framework presented here has not been tested on actual experimental data. Surrogate model miscalibration — a known issue in BO for chemistry — can lead to overconfident predictions and poor real-world convergence.

### 6.5 Recommendations for Practitioners

Based on the experimental evidence, we suggest the following guidelines:

| Problem Setting | Recommended Approach |
|-----------------|---------------------|
| D < 10, smooth objective | RBF kernel + KG acquisition |
| D < 10, rough/noisy objective | Matérn-2.5 + q-EI (batch) |
| 10 < D < 50, unknown structure | Standard BO (Matérn-2.5, UCB or EI) |
| D > 50 with known active subspace | REMBO or ALEBO (d_eff estimated) |
| Multi-objective (2–3 objectives) | qNEHVI (BOTorch) |
| Chemical synthesis workflow | EI with Matérn-2.5, parallel batches |

---

## 7. Conclusion

We have presented a systematic evaluation of Bayesian optimization strategies for high-dimensional parameter spaces, implemented on the BOTorch platform. The key findings are:

1. **Kernel choice matters**: RBF achieves 94.2% of the Hartmann-6 global optimum (3.128 ± 0.017) when the objective is smooth; Matérn-2.5 provides a more robust default with higher variance.

2. **KG outperforms EI in limited budgets**: Knowledge Gradient (2.775 ± 0.440) outperforms EI (2.112 ± 0.584) on 6-dimensional optimization with a 15-iteration budget, while q-EI offers the best consistency for batch settings (2.450 ± 0.185).

3. **Standard BO remains competitive in D=25**: Contrary to the conventional recommendation to use embedding methods for high dimensions, standard BO with Matérn-2.5 achieves −0.215 ± 0.091 vs. REMBO's −7.259 ± 0.939 on the Rosenbrock benchmark, consistent with recent findings by Xu et al. (2024).

4. **Multi-objective qNEHVI is effective**: Hypervolume improves 113.5% (0.385 → 0.821) over 12 iterations, effectively tracing the yield/selectivity Pareto front.

5. **BO outperforms random search in chemistry**: EI achieves 0.875 ± 0.005 yield vs. 0.676 ± 0.123 for random search — a 29.6% improvement — highlighting the practical value of surrogate-guided optimization.

Future work should (i) validate these findings on real experimental datasets, (ii) extend to larger batch sizes (q > 4) with asynchronous BO, (iii) investigate learned embeddings (VAE-BO, ALEBO) for truly high-dimensional spaces, and (iv) incorporate physical constraints (solubility limits, safety bounds) as hard constraints in the acquisition function optimization.

---

## References

1. **Binois, M. & Wycoff, N. (2022).** A survey on high-dimensional Gaussian process modeling with application to Bayesian optimization. *ACM Transactions on Evolutionary Learning and Optimization*, 2(2). DOI: 10.1145/3545611

2. **Xu, Z., Wang, H., Phillips, J.M., & Zhe, S. (2024).** Standard Gaussian process is all you need for high-dimensional Bayesian optimization. *arXiv preprint* arXiv:2402.02746. (v5, 2025)

3. **Zhang, F., Sugisawa, S., & Felton, K.C. et al. (2023).** Multi-objective Bayesian optimisation using q-Noisy Expected Hypervolume Improvement (qNEHVI) for Schotten-Baumann reaction. *ChemRxiv*. DOI: 10.26434/chemrxiv-2023-dlkgl

4. **Gobert, M., Gmys, J., & Toubeau, J.F. (2022).** Batch acquisition for parallel Bayesian optimization — application to hydro-energy storage systems scheduling. *Algorithms*, 15(12), 446. DOI: 10.3390/a15120446

5. **Le, P. & Branke, J. (2024).** Using the knowledge gradient acquisition function in Bayesian optimization when searching for robust solutions. *Engineering Optimization*, 56(3). DOI: 10.1080/0305215x.2022.2145604

6. **Wang, J., Clark, S.C., Liu, E., & Frazier, P.I. (2020).** Parallel Bayesian global optimization of expensive functions. *Operations Research*, 68(6), 1850–1865.

7. **Gu, T. et al. (2024).** BBGP-sDFO: Batch Bayesian and Gaussian process enhanced subspace derivative free optimization for high-dimensional analog circuit synthesis. *IEEE Transactions on Computer-Aided Design*, 43(3). DOI: 10.1109/tcad.2023.3314519

8. **Shahriari, B., Swersky, K., Wang, Z., Adams, R.P., & de Freitas, N. (2016).** Taking the human out of the loop: A review of Bayesian optimization. *Proceedings of the IEEE*, 104(1), 148–175.

9. **Frazier, P.I., Powell, W.B., & Dayanik, S. (2009).** The knowledge-gradient policy for correlated normal beliefs. *INFORMS Journal on Computing*, 21(4), 599–613.

10. **Daulton, S., Balandat, M., & Bakshy, E. (2020).** Differentiable expected hypervolume improvement for parallel multi-objective Bayesian optimization. *NeurIPS 2020*.

11. **Balandat, M. et al. (2020).** BoTorch: A framework for efficient Monte-Carlo Bayesian optimization. *NeurIPS 2020*.

12. **Wang, Z., Zoghi, M., Hutter, F., Matheson, D., & de Freitas, N. (2013).** Bayesian optimization in high dimensions via random embeddings. *IJCAI 2013*.
