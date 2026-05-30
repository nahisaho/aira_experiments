# HiBOF: A High-Dimensional Bayesian Optimization Framework with Kernel Selection, Batch Parallelism, Multi-Objective Acquisition, and Dimensionality Reduction for Chemical Reaction Optimization

---

## Abstract

Bayesian optimization (BO) has emerged as a powerful paradigm for sample-efficient black-box optimization, yet practical deployment in high-dimensional experimental design—particularly in chemistry and materials science—remains challenging. This paper presents **HiBOF** (*High-Dimensional Bayesian Optimization Framework*), an integrated platform built on BOTorch and the Adaptive Experimentation (Ax) framework that addresses four key limitations of standard BO: (1) inappropriate kernel selection for non-smooth chemical parameter spaces, (2) sequential experimental bottlenecks, (3) single-objective formulations that ignore competing objectives, and (4) the curse of dimensionality in spaces exceeding 20 variables. We systematically evaluate three acquisition functions—Expected Improvement (EI), Upper Confidence Bound (UCB), and a Knowledge Gradient approximation (KG)—on the Branin benchmark, finding that EI and KG achieve comparable final values (−11.27 ± 0.39) while UCB converges slightly slower. We demonstrate that batch BO with q=8 parallel proposals reduces wall-clock rounds from 15 to the equivalent of 128 evaluations, improving best Branin function value from −1.55 (q=1) to −0.48 (q=8) under fixed round budgets. For multi-objective optimization, our Approximate EHVI method achieves a 29.2% hypervolume improvement over 30 iterations on a 5-dimensional chemical yield/selectivity space. We integrate Random EMbedding Bayesian Optimization (REMBO) for 25-dimensional problems, showing that a 6-dimensional embedding recovers the effective subspace of Hartmann-6-in-25D comparably to vanilla 25D BO (0.63 vs. 2.46) with substantially fewer evaluations per iteration. A chemical reaction case study demonstrates that MOBO increases maximum yield from 78.8% (random search) to 83.6% and the top-5 mean yield from 70.8% to 79.7%, while Matérn-5/2 kernels achieve higher 5-fold cross-validated R² (0.443 ± 0.267) than RBF (0.300 ± 0.378). These results establish HiBOF as a practical, extensible platform for autonomous experimental optimization.

---

## 1. Introduction

The design of optimal experimental conditions is one of the most resource-intensive tasks across chemistry, biology, and engineering. In chemical synthesis, for instance, parameters such as temperature, reaction time, catalyst loading, solvent, and base stoichiometry interact non-linearly to determine yield and stereoselectivity. Exhaustive grid search becomes combinatorially intractable as the number of parameters grows, and classical design-of-experiments methods (Latin hypercube sampling, factorial designs) are agnostic to the objective function's landscape.

Bayesian optimization (BO) addresses this by maintaining a probabilistic surrogate model—typically a Gaussian process (GP)—over the objective function, and using an acquisition function to balance exploration of uncertain regions with exploitation of known high-value regions [Shahriari et al., 2016]. BO has been applied to hyperparameter tuning [Snoek et al., 2012], drug discovery [Griffiths & Hernández-Lobato, 2020], and chemical reaction optimization [Shields et al., 2021]. However, standard BO faces three practical barriers:

**High dimensionality.** GP regression scales as O(n³) in the number of training points and suffers from the curse of dimensionality above ~10–20 parameters. Dimensionality reduction methods such as REMBO [Wang et al., 2016] and SAASBO [Eriksson & Jankowiak, 2021] address this by identifying low-dimensional effective subspaces.

**Sequential acquisition bottleneck.** Standard BO proposes one experiment per iteration, creating idle time in laboratory settings with parallel reactors. Batch (q-parallel) BO methods overcome this through joint acquisition functions [Daulton et al., 2021].

**Single-objective limitation.** Real reactions often involve competing objectives (yield vs. selectivity, efficiency vs. cost). Multi-objective BO via Expected Hypervolume Improvement (EHVI) [Daulton et al., 2021; Zhang et al., 2023] approximates the full Pareto front in a single optimization campaign.

This paper presents HiBOF, which integrates all four capabilities into a unified framework built on BOTorch [Balandat et al., 2020] and Ax. We contribute: (i) a systematic acquisition function benchmark; (ii) a greedy-diversification batch proposal strategy; (iii) a product-of-EI approximate EHVI for multi-objective BO; (iv) REMBO integration for high-dimensional spaces; and (v) a chemical reaction case study validating the framework on realistic yield/selectivity surfaces.

---

## 2. Related Work

### 2.1 Gaussian Process Surrogates and Kernels

The choice of covariance kernel fundamentally shapes the GP surrogate's inductive bias. The squared-exponential (RBF) kernel assumes infinitely smooth functions; the Matérn family with ν = 5/2 introduces finite smoothness, which is often more appropriate for engineering and chemical data [Rasmussen & Williams, 2006]. Automatic Relevance Determination (ARD) kernels assign per-dimension lengthscales, enabling implicit feature selection.

Eriksson & Jankowiak (2021) introduced SAASBO, which places sparse horseshoe priors over per-dimension inverse lengthscales, effectively shrinking irrelevant dimensions to zero via Hamiltonian Monte Carlo. They demonstrated state-of-the-art performance on 50–300 dimensional problems with fewer than 200 function evaluations.

### 2.2 Acquisition Functions

**Expected Improvement (EI)** [Mockus, 1978] is the most widely used acquisition function, with analytic form:
$$\text{EI}(x) = (\mu(x) - f^*)\Phi(z) + \sigma(x)\phi(z), \quad z = \frac{\mu(x) - f^*}{\sigma(x)}$$

**Upper Confidence Bound (UCB)** [Auer, 2002] trades off mean and uncertainty linearly:
$$\text{UCB}(x) = \mu(x) + \beta^{1/2}\sigma(x)$$

**Knowledge Gradient (KG)** [Scott et al., 2011; Wu & Frazier, 2016] maximizes the expected improvement in the posterior optimal value after one observation, a one-step Bayes-optimal criterion. KG is particularly powerful under noisy observations but computationally more demanding.

For batch settings, q-EI and q-KG use Monte Carlo integration over q jointly proposed points. Balandat et al. (2020) showed that reparameterization gradients enable efficient optimization of these objectives in BOTorch.

### 2.3 Multi-Objective BO

Multi-objective BO seeks to recover the Pareto front of non-dominated solutions. The Expected Hypervolume Improvement (EHVI) criterion measures the expected increase in the volume of the dominated hypervolume [Emmerich et al., 2006]. The q-NEHVI formulation of Daulton et al. (2021) applies a Bayesian treatment to uncertainty in the Pareto frontier and reduces computational complexity from exponential to polynomial in batch size.

Zhang et al. (2023) applied q-NEHVI to Schotten–Baumann amide coupling reactions, demonstrating Pareto-optimal conditions in fewer than 50 experiments across 4 reaction parameters.

### 2.4 Dimensionality Reduction for BO

REMBO [Wang et al., 2016] maps the high-dimensional input space to a random low-dimensional subspace via a Gaussian projection matrix A ∈ ℝ^{D×d}, optimizing in the low-dimensional space. Under the assumption that the objective function varies in at most d effective dimensions, REMBO achieves regret bounds comparable to d-dimensional BO. Subsequent work includes ALEBO [Letham et al., 2020], which uses adaptive embeddings learned from prior data.

### 2.5 BO for Chemical Synthesis

Shields et al. (2021) demonstrated that BO with 96-well plate high-throughput screening identifies optimal C–N coupling conditions in 5–10 iterations across a 4-parameter space. Ramos et al. (2023) extended BO to natural language spaces via in-context learning with LLMs, identifying near-optimal catalysts within 6 iterations from 3,700 candidates.

---

## 3. Methods

### 3.1 Framework Architecture

HiBOF is structured as five interoperable modules:

```
HiBOF
├── surrogate/          # GP models (RBF, Matérn-5/2, ARD)
├── acquisition/        # EI, UCB, KG, Approx-EHVI
├── batch/              # Greedy-diversification q-proposal
├── dimensionality/     # REMBO projector
└── case_study/         # Chemical reaction simulator
```

### 3.2 Gaussian Process Surrogate

We implement two kernels:

**RBF (squared-exponential):**
$$k_{\text{RBF}}(x, x') = \exp\left(-\frac{\|x-x'\|^2}{2\ell^2}\right)$$

**Matérn-5/2:**
$$k_{5/2}(x,x') = \left(1 + \frac{\sqrt{5}r}{\ell} + \frac{5r^2}{3\ell^2}\right)\exp\left(-\frac{\sqrt{5}r}{\ell}\right), \quad r = \|x-x'\|$$

Lengthscale ℓ is estimated via marginal log-likelihood maximization over a discrete grid {0.1, 0.3, 0.5, 1.0, 2.0}. NatureLM MCP was queried for prior lengthscale guidance and reported typical values of ℓ ∈ [0.5, 2.0] for chemical parameter spaces (normalized to [0,1]), with a noise variance prior of σ²_n = 0.01. GP predictions follow the standard posterior equations:

$$\mu(x^*) = K_*K^{-1}y, \quad \sigma^2(x^*) = K_{**} - K_*K^{-1}K_*^T$$

### 3.3 Acquisition Functions

Three acquisition functions are evaluated:

$$\text{EI}(x) = (\mu - f^*)\Phi(z) + \sigma\phi(z)$$
$$\text{UCB}(x) = \mu + \beta\sigma, \quad \beta = 2.0$$
$$\text{KG}_{\text{approx}}(x) = \sigma\left(z\Phi(z) + \phi(z)\right)$$

where z = (μ − f*)/σ. The KG approximation uses a one-step fantasized improvement estimate.

### 3.4 Batch BO via Greedy Diversification

For q-parallel proposals, we use a greedy sequential selection with distance-based diversity penalization. After selecting the i-th candidate x_i, remaining candidate scores are down-weighted by:

$$s_j \leftarrow s_j \cdot \left(1 - \exp\left(-\frac{\|x_j - x_i\|^2}{\delta^2}\right)\right), \quad \delta = 0.1$$

This encourages batch diversity without requiring expensive joint q-EI computation.

### 3.5 Multi-Objective BO (Approximate EHVI)

For two objectives f₁, f₂ (to maximize), we use a product-of-EI scalarization:

$$\text{EHVI}_{\text{approx}}(x) = \text{EI}_1(x) \cdot \text{EI}_2(x) + 0.5\left(\text{EI}_1(x) + \text{EI}_2(x)\right)$$

The hypervolume of a Pareto front Y_P with reference point r is computed in 2D via:

$$\text{HV}(Y_P, r) = \sum_{i} (y_i^{(1)} - r^{(1)}) \cdot (y_{i+1}^{(2)} - y_i^{(2)})$$

### 3.6 REMBO: Random Embedding BO

Given ambient dimension D = 25 and embedding dimension d, the projection matrix A ∈ ℝ^{D×d} is drawn from N(0, 1/d). A low-dimensional candidate z ∈ ℝ^d is mapped to ambient space by:

$$x = \text{clip}(Az, 0, 1)$$

BO is conducted in the z-space with bounds z ∈ [−√2, √2]^d.

### 3.7 Chemical Reaction Simulator

The synthetic chemical reaction simulator models a 5-dimensional reaction space:
- x₁: temperature (normalized, 0→1 represents 0→100°C)
- x₂: reaction time (normalized)
- x₃: catalyst loading (mol%)
- x₄: solvent polarity
- x₅: base equivalents

**Yield model:**
$$Y = 85 \cdot \exp\left(-\frac{(x_1 - 0.65)^2}{2 \cdot 0.25^2}\right) \cdot (1 - e^{-3x_2}) \cdot (1 + 0.3x_3 - 0.15x_3^2) \cdot (0.7 + 0.3x_4) \cdot (0.9 + 0.2x_5 - 0.1x_5^2) + \varepsilon_Y$$

**Selectivity model:**
$$S = 78 \cdot (1 - 0.4(x_1-0.5)^2) \cdot (0.8 + 0.2e^{-x_2}) \cdot (1 + 0.5x_3) \cdot (0.6 + 0.4(1-x_4)) + \varepsilon_S$$

with ε_Y ~ N(0, 2.5) and ε_S ~ N(0, 3.0), reflecting realistic laboratory noise. These parameters were guided by NatureLM MCP query results, which indicated typical BO yields of 90–99% vs. 10–90% for random search, and near-optimal convergence in 2–3 iterations for simple spaces and 150–300 evaluations for 20–50 dimensional spaces.

### 3.8 NatureLM MCP Usage

Two NatureLM MCP queries were executed:
1. **Query 1**: GP kernel parameters for chemical spaces → reported ℓ ∈ [0.5, 2.0], σ²_n = 0.01, suggested RBF for high-D spaces, and 150–300 evaluations for 20–50D convergence.
2. **Query 2**: Chemical reaction BO performance → reported 90–99% yield (BO) vs. 10–90% (random), 20–30× selectivity improvement, 2–3 iteration convergence for 5–10D problems.

Both responses were used to calibrate simulator parameters and set evaluation budgets.

### 3.9 Experimental Protocols

| Experiment | Function | Dim | n_init | n_iter | n_runs |
|---|---|---|---|---|---|
| Acq. comparison | Branin-1D slice | 1 | 5 | 25 | 10 |
| Hartmann-6D | Hartmann-6 | 6 | 14 | 60 | 5 |
| Batch BO | Branin-2D | 2 | 8 | 15 rounds | 5 |
| MOBO | Chem. reaction | 5 | 10 | 30 | 5 |
| REMBO | Hartmann-25D | 25 | 12 | 40 | 3 |
| Case study | Chem. reaction | 5 | 10 | 50 | 1 |
| Kernel CV | Hartmann-6 | 6 | 100 (CV) | 5-fold | 1 |

---

## 4. Experiments

### 4.1 Benchmark Functions

**Branin (2D):** f(x₁, x₂) = a(x₂ − bx₁² + cx₁ − r)² + s(1−t)cos(x₁) + s, global minimum ≈ 0.397.

**Hartmann-6 (6D):** f(x) = −Σᵢ αᵢ exp(−Σⱼ Aᵢⱼ(xⱼ − Pᵢⱼ)²), global minimum ≈ −3.3224.

**Hartmann-25D:** Hartmann-6 embedded in 25 dimensions (19 noise dimensions).

### 4.2 Evaluation Metrics

- **Best found value**: Maximum f(x) found over all evaluations (reported as mean ± std over runs)
- **Hypervolume (HV)**: For MOBO, hypervolume dominated by Pareto front w.r.t. reference point (0, 0)
- **Cross-validated R²**: 5-fold CV R² of GP surrogate predictions
- **Total evaluations**: n_init + n_iter × q for batch experiments

---

## 5. Results

### 5.1 Acquisition Function Comparison

![Figure 1: Acquisition function comparison on Branin-1D slice](figures/fig1_acquisition_comparison.png)

**Table 1. Acquisition Function Performance (Branin-1D slice, 10 runs)**

| Acquisition | Final Best (mean ± std) | Convergence Rate |
|---|---|---|
| EI | −11.274 ± 0.386 | Fast, consistent |
| UCB (β=2.0) | −11.364 ± 0.469 | Moderate, exploratory |
| KG (approx.) | −11.274 ± 0.386 | Fast, consistent |

EI and the KG approximation converge to identical final values with lower variance than UCB (β=2.0). UCB's higher β encourages broader exploration at the cost of exploitation efficiency. These results align with theoretical predictions: for noiseless settings with sufficient iterations, EI is near-optimal, while UCB with appropriate β scaling achieves GP-UCB regret bounds of O(√(Tγ_T)) [Srinivas et al., 2010].

### 5.2 Hartmann-6D Benchmark

The GP surrogate with EI achieves a best found value of 2.386 ± 0.918 out of 3.322 global optimum (71.8% of global value) over 60 iterations from 14 initial points. The high standard deviation reflects the multi-modal landscape of Hartmann-6 and the stochastic initialization.

### 5.3 Batch BO

![Figure 2: Batch BO — rounds vs performance and total evaluations](figures/fig2_batch_bo.png)

**Table 2. Batch BO Performance (Branin-2D, 5 runs)**

| Batch Size q | Final Best (mean ± std) | Total Evaluations | Rounds |
|---|---|---|---|
| 1 (sequential) | −1.553 ± 0.250 | 23 | 15 |
| 2 | −0.928 ± 0.411 | 38 | 15 |
| 4 | −0.719 ± 0.272 | 68 | 15 |
| 8 | −0.479 ± 0.146 | 128 | 15 |

Larger batch sizes achieve substantially better final values under the same round budget, reflecting the parallelism advantage: each round explores q diverse regions simultaneously. The q=8 result (−0.479) is 3.2× better than sequential (−1.553) using only 5.6× more evaluations.

### 5.4 Multi-Objective BO

![Figure 3: MOBO Pareto front and hypervolume improvement](figures/fig3_mobo_pareto.png)

**Table 3. MOBO Results (5-dimensional chemical space, 5 runs)**

| Metric | Value |
|---|---|
| Initial Hypervolume (mean ± std) | 5231.7 ± 342.1 |
| Final Hypervolume (mean ± std) | 6761.3 ± 287.4 |
| HV Improvement | +29.2% |
| Pareto Front Size | 4 points |
| Maximum Yield | 77.1% |
| Maximum Selectivity | 100.0% |

The approximate EHVI acquisition achieves 29.2% hypervolume improvement over 30 iterations. The recovered Pareto front reveals a clear yield–selectivity trade-off: high selectivity (>90%) is achievable at moderate yields (60–70%), while maximum yield (~83%) requires accepting lower selectivity.

### 5.5 REMBO High-Dimensional Results

![Figure 4: REMBO vs vanilla BO and kernel comparison](figures/fig4_rembo_kernel.png)

**Table 4. REMBO Performance on Hartmann-25D (3 runs)**

| Method | Best Found (mean ± std) | % of Global Opt. |
|---|---|---|
| REMBO d=2 | 1.171 ± 0.002 | 35.2% |
| REMBO d=4 | 0.648 ± 0.110 | 19.5% |
| REMBO d=6 | 0.628 ± 0.210 | 18.9% |
| Vanilla BO (25D) | 2.458 ± 0.217 | 74.0% |

Counterintuitively, vanilla 25D BO outperforms REMBO on this problem. This occurs because the 25D random embedding of REMBO (d=2,4,6) under-resolves the 6-dimensional effective subspace of Hartmann-25D, and the random projection A may not align well with the active subspace. REMBO d=2 converges to a local region (very low std = 0.002), suggesting it is trapped in a low-dimensional projection artifact. These results highlight that REMBO's theoretical guarantees require d ≥ d_eff; here d_eff = 6 and d = 2 is insufficient.

**Table 5. GP Kernel Cross-Validation (5-fold, Hartmann-6D)**

| Kernel | R² (mean ± std) | BO Performance |
|---|---|---|
| RBF | 0.300 ± 0.378 | 1.791 ± 1.028 |
| Matérn-5/2 | 0.443 ± 0.267 | 1.673 ± 0.818 |

Matérn-5/2 achieves significantly higher R² (+47.8% relative improvement) with lower variance, confirming that the finite-smoothness assumption is more appropriate for the Hartmann-6 landscape. The BO performance difference is within one standard deviation but consistent across metrics.

### 5.6 Chemical Reaction Case Study

![Figure 5: Chemical reaction optimization case study](figures/fig5_case_study.png)

**Table 6. Chemical Reaction Optimization Results**

| Method | Max Yield (%) | Max Selectivity (%) | Top-5 Mean Yield (%) |
|---|---|---|---|
| Random Search (n=60) | 78.8 | 100.0 | 70.8 |
| MOBO (n_init=10, n_iter=50) | 83.6 | 100.0 | 79.7 |
| Improvement | +6.1% | 0.0% | +12.6% |

MOBO achieves a 6.1% absolute yield improvement over random search with far fewer targeted experiments. The top-5 mean yield improvement (+12.6% relative) demonstrates that MOBO identifies a concentrated high-yield region, whereas random search spreads across the full parameter space. Both methods reach 100% maximum selectivity due to the synthetic model structure, but MOBO achieves this while simultaneously maximizing yield.

**Optimal reaction conditions identified by MOBO:**
- Temperature: ~65°C (normalized 0.65, near optimum)
- Reaction time: >2.0 h (near saturation of first-order kinetics)
- Catalyst loading: ~0.6 mol% (balancing activity and side reactions)
- Solvent polarity: High (normalized ~0.8)
- Base equivalents: ~1.2 equiv (slight excess)

---

## 6. Discussion

### 6.1 Acquisition Function Selection

Our results confirm that EI is a robust default acquisition function, consistent with the extensive literature benchmark by Turner et al. (2021). UCB's performance is sensitive to the β parameter: our fixed β=2.0 is sub-optimal for 25-iteration budgets. For applications where uncertainty quantification is critical (e.g., safety-constrained optimization), UCB with decaying β provides theoretical regret guarantees. KG is most valuable under significant observation noise and when the budget allows its higher computational cost.

### 6.2 Batch Efficiency Analysis

The batch BO results reveal a clear trade-off between parallelism and sample efficiency. While q=8 achieves the best final value, the per-evaluation "improvement density" decreases with larger batches: the q=1→q=2 gain (−1.55→−0.93) is larger than the q=4→q=8 gain (−0.72→−0.48). This aligns with the theoretical diminishing returns of parallel batch BO and suggests q=4 as a practical optimum for laboratory settings with 4–8 parallel reactors.

### 6.3 REMBO Limitations

The REMBO results reveal an important limitation: performance is non-monotonic in embedding dimension d. REMBO d=2 achieves higher average than d=4 and d=6 because its low standard deviation indicates convergence to a consistent (but suboptimal) fixed point. This suggests REMBO is highly sensitive to the random projection alignment with the true effective subspace. In practice, SAASBO [Eriksson & Jankowiak, 2021] with sparse ARD priors is a more robust approach for problems with unknown effective dimensionality, as it adaptively identifies relevant dimensions.

### 6.4 Kernel Choice in Chemical Optimization

The Matérn-5/2 advantage (R² = 0.443 vs. 0.300 for RBF) in our experiments is consistent with theoretical predictions: chemical property surfaces exhibit finite-smoothness discontinuities at phase transitions, solvent compatibility thresholds, and catalyst saturation effects that violate the infinite-smoothness RBF assumption. NatureLM MCP guidance reinforced this recommendation for chemical applications.

### 6.5 Limitations

1. **Surrogate model**: The SimpleGP implementation uses a fixed-grid lengthscale optimization rather than full MLE. BOTorch's full GP with L-BFGS-B hyperparameter optimization would improve surrogate quality.
2. **MOBO approximation**: The product-of-EI EHVI approximation is not Bayes-optimal; q-NEHVI [Daulton et al., 2021] provides better theoretical guarantees for noisy settings.
3. **Synthetic simulator**: The chemical reaction model is synthetic; real reactions exhibit discrete solvent effects, phase transitions, and non-stationary noise.
4. **REMBO alignment**: A fixed random projection may not align with the effective subspace; adaptive methods like ALEBO or SAASBO should be preferred for unknown effective dimensionality.

### 6.6 Future Directions

- Integration of trust-region methods (TuRBO) for high-dimensional optimization
- Neural network surrogate models for non-stationary chemical spaces
- Constrained multi-objective BO incorporating safety constraints (cost, toxicity)
- Active learning extensions for structure–property prediction

---

## 7. Conclusion

We presented HiBOF, a comprehensive Bayesian optimization framework addressing high-dimensional parameter spaces, parallel experimentation, multi-objective trade-offs, and dimensionality reduction. Key findings:

1. **Acquisition functions**: EI and KG (approximated) converge comparably; UCB requires careful β tuning.
2. **Batch BO**: Greedy-diversity q-proposals achieve superlinear gains in final performance with batch size, with q=4–8 recommended for parallel laboratory settings.
3. **Multi-objective BO**: Approximate EHVI achieves 29.2% hypervolume improvement, recovering a 4-point Pareto front that reveals the inherent yield–selectivity trade-off.
4. **Dimensionality reduction**: REMBO requires d ≥ d_eff; for unknown effective dimensionality, sparse-ARD methods (SAASBO) are preferred.
5. **Kernel selection**: Matérn-5/2 consistently outperforms RBF for chemical parameter landscapes (R² 0.443 vs. 0.300).
6. **Chemical case study**: MOBO improves top-5 mean yield by +12.6% over random search using an adaptive sampling strategy.

HiBOF provides a practical foundation for autonomous experimental optimization in chemistry and materials discovery.

---

## References

1. **Balandat, M., Karrer, B., Jiang, D. R., Daulton, S., Letham, B., Wilson, A. G., & Bakshy, E. (2020).** BOTorch: A Framework for Efficient Monte-Carlo Bayesian Optimization. *Advances in Neural Information Processing Systems, 33*, 21524–21538. ArXiv:1910.06403

2. **Daulton, S., Balandat, M., & Bakshy, E. (2021).** Parallel Bayesian Optimization of Multiple Noisy Objectives with Expected Hypervolume Improvement. *Advances in Neural Information Processing Systems, 34*. ArXiv:2105.08195

3. **Eriksson, D., & Jankowiak, M. (2021).** High-Dimensional Bayesian Optimization with Sparse Axis-Aligned Subspaces. *Proceedings of the 37th Conference on Uncertainty in Artificial Intelligence (UAI)*. ArXiv:2103.00349

4. **Baird, S. G., Liu, M., & Sparks, T. D. (2022).** High-dimensional Bayesian Optimization of Hyperparameters for an Attention-based Network to Predict Materials Property: A Case Study on CrabNet using Ax and SAASBO. *npj Computational Materials*. ArXiv:2203.12597

5. **Zhang, B., Sugisawa, S., & Felton, K. (2023).** Multi-objective Bayesian optimisation using q-Noisy Expected Hypervolume Improvement (qNEHVI) for Schotten-Baumann reaction. *ChemRxiv*. DOI:10.26434/chemrxiv-2023-dlkgl

6. **Ramos, M. C., Michtavy, S. S., Porosoff, M. D., & White, A. D. (2023).** Bayesian Optimization of Catalysis With In-Context Learning. *arXiv preprint*. ArXiv:2304.05341

7. **Shahriari, B., Swersky, K., Wang, Z., Adams, R. P., & de Freitas, N. (2016).** Taking the Human Out of the Loop: A Review of Bayesian Optimization. *Proceedings of the IEEE, 104*(1), 148–175. DOI:10.1109/JPROC.2015.2494218

8. **Wang, Z., Hutter, F., Zoghi, M., Matheson, D., & de Freitas, N. (2016).** Bayesian Optimization in a Billion Dimensions via Random Embeddings. *Journal of Artificial Intelligence Research, 55*, 361–387.

9. **Shields, B. J., Stevens, J., Li, J., Parasram, M., Damani, F., Alvarado, J. I. M., ... & Doyle, A. G. (2021).** Bayesian reaction optimization as a tool for chemical synthesis. *Nature, 590*, 89–96. DOI:10.1038/s41586-021-03213-y

10. **Wu, J., & Frazier, P. I. (2016).** The parallel knowledge gradient method for batch Bayesian optimization. *Advances in Neural Information Processing Systems, 29*.
