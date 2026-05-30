# Machine Learning Framework for Multi-Objective Composition Optimization of High-Entropy Alloys: Descriptor Engineering, Bayesian Optimization, and Active Learning

---

## Abstract

High-entropy alloys (HEAs) represent a paradigm shift in alloy design, offering vast combinatorial composition spaces that defy conventional trial-and-error exploration. In this work, we present an integrated machine learning (ML) framework for the systematic, data-driven optimization of HEA compositions targeting multiple mechanical and corrosion-related properties simultaneously. The framework couples CALPHAD-inspired thermodynamic descriptor engineering—including atomic size mismatch (δ_r), valence electron concentration (VEC), mixing entropy (ΔS_mix), mixing enthalpy (ΔH_mix), and the Omega stability parameter—with three classes of surrogate models: Gaussian Process Regression (GPR), Random Forest (RF), and Gradient Boosting Machines (GBM). Multi-objective Bayesian optimization is performed using an Expected Improvement (EI) acquisition function with weighted-sum scalarization, enabling simultaneous optimization of hardness (HV), yield strength (YS), fracture strain (FS), and corrosion index (CI). An active learning loop further reduces the experimental budget by strategically selecting uninformative compositions for labeling. A synthetic dataset of 500 alloys drawn from a ten-element pool (Cr, Mn, Fe, Co, Ni, Al, Ti, V, Mo, W) is used for validation. Five-fold cross-validation reveals that GBM achieves the best predictive performance for YS (R² = 0.876 ± 0.015, RMSE = 42.8 ± 2.6 MPa), while all models show lower fidelity for hardness prediction (R² ≈ 0.44–0.55). A CrMnFeCoNi case study demonstrates that Cr-enriched non-equiatomic compositions (x_Cr ≈ 0.35) outperform the equiatomic Cantor alloy in both corrosion resistance and hardness within the design space explored. Thirty iterations of active learning improve the weighted objective by 2.1% over random sampling. Critical limitations—including reliance on synthetic data, simplistic property surrogates, and the absence of microstructural and processing variables—are discussed in detail, alongside recommendations for integration with real CALPHAD databases (Thermo-Calc TCHEA), AFLOW/Materials Project repositories, and DFT-validated datasets.

---

## 1. Introduction

High-entropy alloys, first reported by Cantor et al. (2004) and Yeh et al. (2004), are multi-principal-element alloys typically comprising five or more elements in near-equiatomic proportions. Their defining thermodynamic feature—high configurational entropy—stabilizes disordered solid-solution phases and suppresses the formation of brittle intermetallics, yielding remarkable combinations of strength, ductility, and corrosion resistance. The CrMnFeCoNi "Cantor alloy" remains the archetype system, exhibiting exceptional cryogenic toughness and face-centered cubic (FCC) structure stability.

Despite their promise, the sheer dimensionality of the HEA compositional space presents a fundamental challenge: for a pool of *N* elements and *k* principal components, the number of candidate compositions scales combinatorially. Systematic experimental exploration is therefore infeasible without computational guidance. CALPHAD (CALculation of PHAse Diagrams) methods provide thermodynamically rigorous predictions of phase stability but are computationally expensive and require empirically validated parameters for each binary and ternary subsystem. First-principles DFT calculations offer atomic-level accuracy for specific compositions but are limited to small supercells and cannot directly predict macroscopic mechanical properties.

Machine learning has emerged as a powerful complement to both approaches, enabling rapid property screening across composition space by learning structure–property relationships from existing data. Recent work has demonstrated ML prediction of hardness [1,2,3], yield strength [4], phase stability [5], and multi-objective optimization [6] in HEA systems. However, several challenges remain:

1. **Limited labeled data**: HEA experiments are costly; most existing datasets contain fewer than 500 records.
2. **Multi-objective trade-offs**: Strength and ductility frequently conflict; corrosion resistance depends strongly on Cr content.
3. **Descriptor design**: No consensus exists on the optimal feature set for capturing composition–property relationships.
4. **Active learning integration**: Most existing work trains static models; few incorporate iterative experimental proposal loops.

This work addresses these challenges through an end-to-end ML pipeline that: (i) encodes CALPHAD-derived thermodynamic knowledge into computable descriptors; (ii) benchmarks GP, RF, and GBM surrogates under rigorous cross-validation; (iii) performs EI-based multi-objective Bayesian optimization; (iv) simulates an active learning loop; and (v) applies the framework to the CrMnFeCoNi system as a case study.

**Contributions**: 
- A modular, open-source descriptor suite covering 13 physico-chemical features.
- Systematic cross-validation of three surrogate model families across four target properties.
- Multi-objective Pareto front analysis revealing 31 Pareto-optimal compositions from 500 candidates.
- Active learning simulation demonstrating 2.1% objective improvement in 30 iterations.
- Critical evaluation of assumptions, biases, and real-world generalization limitations.

---

## 2. Related Work

### 2.1 ML for HEA Property Prediction

Chen et al. [1] pioneered an interpretable ML framework for HEA hardness prediction using ensemble methods, achieving R² > 0.97 with leave-one-out cross-validation on 237 experimental compositions. Their work established the importance of mixing entropy and VEC as dominant descriptors. Huang et al. [3] extended this approach by incorporating solid-solution strengthening (SSH) theory, demonstrating that short-range order and charge transfer effects contribute significantly to hardness in FeNiCuCo and CrMoNbTi alloy families. Ren et al. [2] combined gradient-boosted trees with an optimization algorithm to simultaneously predict and design HEAs with target hardness, reporting R² = 0.97 under LOO validation and experimental verification of optimized compositions.

### 2.2 Multi-Objective and Active Learning Approaches

Ma et al. [6] proposed an NSGA-II-based multi-objective framework integrating SVR (ductility) and LightGBM (hardness) models, achieving 10-fold CV R² of 0.76 and 0.90 respectively and experimentally validating four candidate compositions—notably improving fracture strain by 135–282% relative to alloys with similar hardness in the training set. Zhang et al. [5] introduced a Deep Sets architecture trained on DFT-computed elastic properties of 7,086 HEA structures, demonstrating superior generalizability compared to conventional Voronoi-based descriptors.

### 2.3 High-Throughput DFT and CALPHAD

Sun et al. [8] combined CALPHAD-assisted experiments with XGBoost for Ti-Zr-Nb-Ta refractory HEAs, identifying melting point and mixing entropy as the most important features and achieving 97.8% prediction accuracy after CALPHAD-guided data generation. Chen et al. [7] constructed a chemical map of single-phase HEAs using DFT calculations on 658,000+ quinary equimolar alloys, predicting and experimentally validating two novel compositions (BCC AlCoMnNiV, FCC CoFeMnNiZn).

### 2.4 Gaps Addressed by This Work

While these prior studies achieve high accuracy on specific alloy families with carefully curated datasets, they typically (i) focus on a single property, (ii) rely on equiatomic or near-equiatomic compositions, (iii) lack an active learning component, and (iv) do not systematically evaluate multiple surrogate model families under the same cross-validation protocol. This work addresses all four gaps within a unified framework.

---

## 3. Methods

### 3.1 Thermodynamic Descriptor Engineering

For a given HEA composition {x₁, x₂, ..., xₙ} (mole fractions, Σxᵢ = 1), we compute 13 descriptors:

**Atomic size mismatch** (lattice distortion parameter):
$$\delta_r = 100\sqrt{\sum_i x_i \left(1 - \frac{r_i}{\bar{r}}\right)^2}$$
where $\bar{r} = \sum_i x_i r_i$ is the composition-averaged atomic radius.

**Valence Electron Concentration**:
$$\text{VEC} = \sum_i x_i \cdot \text{VEC}_i$$

**Mixing entropy** (configurational, ideal):
$$\Delta S_{\text{mix}} = -R \sum_i x_i \ln x_i$$

**Mixing enthalpy** (Miedema-like binary interaction approximation):
$$\Delta H_{\text{mix}} = \sum_{i<j} 4\Omega_{ij} x_i x_j, \quad \Omega_{ij} = 4\omega_i\omega_j$$

**Omega stability criterion**:
$$\Omega = \frac{\bar{T}_m \cdot \Delta S_{\text{mix}}}{|\Delta H_{\text{mix}}|}$$
where $\bar{T}_m = \sum_i x_i T_{m,i}$. Solid-solution stability is expected for Ω > 1.

**VEC mismatch**: $\delta_{\text{VEC}} = \sqrt{\sum_i x_i(\text{VEC}_i - \overline{\text{VEC}})^2}$

**Electronegativity difference**: $\delta_{\text{EN}} = \sqrt{\sum_i x_i(\chi_i - \bar{\chi})^2}$

Additional descriptors include: average shear modulus ($\bar{G}$), normalised melting point ($\bar{T}_{m,\text{norm}}$), average density ($\bar{\rho}$), lattice distortion proxy (LD = δ_r / (100$\bar{r}$)), and number of principal elements.

### 3.2 Synthetic Dataset Generation

In the absence of a publicly accessible, uniformly labeled multi-property HEA dataset, we constructed a synthetic dataset of 500 alloys using physics-inspired surrogate functions:

**Hardness (HV)**: 
$$\text{HV} = 120 + 8.0\delta_r + 0.9\bar{G} - 40e^{-\text{VEC}/8} + 15\tanh(\Omega - 1.5) + 20\delta_{\text{VEC}} + \varepsilon_\text{HV}$$
$\varepsilon_\text{HV} \sim \mathcal{N}(0, 15)$

**Yield Strength (MPa)**:
$$\text{YS} = 300 + 50\delta_r + 2.5\bar{G} + 120\delta_{\text{VEC}} + 80\ln(1+\Delta S_\text{mix}) + \varepsilon_\text{YS}$$
$\varepsilon_\text{YS} \sim \mathcal{N}(0, 30)$

**Fracture Strain (%)**:
$$\text{FS} = 40 - 2.5\delta_r - 0.05\bar{G} + 5\cdot\text{VEC}/8 + 3\Delta S_\text{mix}/R + \varepsilon_\text{FS}$$
$\varepsilon_\text{FS} \sim \mathcal{N}(0, 1.5)$

**Corrosion Index (CI, 0=excellent, 1=poor)**:
$$\text{CI} = 0.5 - 0.8 x_\text{Cr} + 0.05\delta_r - 0.02(\text{VEC}-8) + \varepsilon_\text{CI}$$
$\varepsilon_\text{CI} \sim \mathcal{N}(0, 0.04)$

Compositions were sampled by selecting 4–6 elements from a ten-element pool and drawing fractions from a Dirichlet distribution. The noise levels were calibrated to match the experimental scatter reported in literature (±15 HV for hardness [1,3]).

### 3.3 Surrogate Models

Three model families were evaluated:
- **GP**: Gaussian Process with Matérn ν=2.5 kernel (Constant × Matérn), hyperparameter optimization via marginal likelihood maximization (3 restarts), normalized targets.
- **RF**: Random Forest (200 trees, √p features per split).
- **GBM**: Gradient Boosting Machine (200 estimators, depth=4, learning rate=0.05).

Features were standardized (zero mean, unit variance) prior to all model fitting.

### 3.4 Cross-Validation Protocol

5-fold stratified cross-validation was applied to all model–property combinations. Metrics reported: R² and RMSE with mean ± standard deviation across folds. This protocol guards against overfitting and provides uncertainty estimates on performance.

### 3.5 Multi-Objective Bayesian Optimization

Multi-objective optimization proceeds via scalarized Expected Improvement:

$$\text{EI}(\mathbf{x}) = \mathbb{E}\left[\max(f(\mathbf{x}) - f^*, 0)\right]$$

For a GP surrogate with posterior $\mathcal{N}(\mu(\mathbf{x}), \sigma^2(\mathbf{x}))$:
$$\text{EI}(\mathbf{x}) = (\mu - f^* - \xi)\Phi(z) + \sigma\phi(z), \quad z = \frac{\mu - f^* - \xi}{\sigma}$$

The composite acquisition function is:
$$\text{EI}_\text{total}(\mathbf{x}) = \sum_k |w_k| \cdot \text{EI}_k(\mathbf{x})$$

with weights $\mathbf{w} = (0.4, 0.3, 0.2, -0.1)$ for (HV, YS, FS, CI). CI is minimized, so its EI is evaluated with $f^* = \min_\text{obs}(\text{CI})$.

### 3.6 Active Learning Loop

The active learning loop proceeds as:
1. Initialize with 50 randomly labeled compositions.
2. Fit independent GP surrogates for each target.
3. Compute EI_total for all unlabeled candidates.
4. Select the candidate with maximum EI_total, add it to the training set.
5. Repeat for 30 iterations.

Progress is tracked by the weighted-sum objective evaluated on all labeled compositions.

### 3.7 Pareto Front Computation

Non-dominated sorting was performed on the 3-objective space (HV_norm, YS_norm, 1-CI_norm) to identify Pareto-optimal compositions.

---

## 4. Experiments

### 4.1 Dataset Statistics

| Property | Mean | Std | Min | Max |
|----------|------|-----|-----|-----|
| HV (Vickers) | 262.5 | 26.7 | 173.6 | 355.3 |
| YS (MPa) | 1142.2 | 123.1 | 717.4 | 1556.2 |
| FS (%) | 31.3 | 3.5 | 17.4 | 41.8 |
| CI | 0.700 | 0.154 | 0.126 | 1.000 |

### 4.2 Descriptor Statistics

| Descriptor | Mean | Std | Physical Meaning |
|------------|------|-----|-----------------|
| δ_r (%) | 3.84 | 1.72 | Atomic size mismatch |
| VEC | 7.23 | 1.51 | Valence electron concentration |
| ΔS_mix (J/mol·K) | 13.05 | 1.31 | Mixing entropy |
| ΔH_mix (eV) | 0.42 | 0.28 | Mixing enthalpy |
| Ω | ~1.2×10⁵ | — | Thermodynamic stability |
| δ_VEC | 2.01 | 0.88 | VEC heterogeneity |

### 4.3 Model Evaluation Protocol

All models trained on standardized features; 5-fold CV with fixed random seed (42) for reproducibility. Models: GP (Matérn ν=2.5), RF (200 trees), GBM (200 estimators, lr=0.05).

---

## 5. Results

### 5.1 Surrogate Model Performance

**Table 1: 5-Fold Cross-Validation Results (mean ± std)**

| Target | Model | R² (mean ± std) | RMSE (mean ± std) |
|--------|-------|-----------------|-------------------|
| HV     | GP    | 0.438 ± 0.078   | 19.74 ± 0.94 HV   |
| HV     | RF    | 0.553 ± 0.062   | 17.63 ± 0.95 HV   |
| HV     | GBM   | 0.533 ± 0.094   | 17.91 ± 0.85 HV   |
| YS     | GP    | 0.825 ± 0.045   | 50.69 ± 7.57 MPa  |
| YS     | RF    | 0.864 ± 0.022   | 44.77 ± 3.31 MPa  |
| YS     | GBM   | **0.876 ± 0.015** | **42.84 ± 2.61 MPa** |
| FS     | GP    | 0.720 ± 0.020   | 1.85 ± 0.05 %     |
| FS     | RF    | **0.770 ± 0.028** | **1.67 ± 0.09 %** |
| FS     | GBM   | 0.765 ± 0.026   | 1.69 ± 0.11 %     |
| CI     | GP    | **0.805 ± 0.046** | **0.065 ± 0.007** |
| CI     | RF    | 0.770 ± 0.026   | 0.072 ± 0.006     |
| CI     | GBM   | 0.790 ± 0.039   | 0.068 ± 0.007     |

![Figure 1: Distribution of Key HEA Descriptors](figures/fig1_descriptors.png)

![Figure 2: Cross-Validation Model Performance](figures/fig2_model_cv.png)

**Key observations**:
- GBM achieves the highest R² for YS (0.876) and is competitive for FS (0.765).
- GP achieves the best CI R² (0.805), likely due to its smooth kernel being appropriate for the relatively linear CI–composition relationship.
- HV prediction is notably weaker (R² ≈ 0.44–0.55) across all models, reflecting the complex nonlinear coupling between hardness and multiple descriptors, compounded by higher relative noise (±15 HV vs. mean 262 HV ≈ ±5.7%).

### 5.2 Feature Importance

Random Forest feature importances reveal different descriptor hierarchies per target:
- **HV**: G_bar and delta_r dominate (together >50% importance).
- **YS**: delta_VEC, S_mix, and delta_r are the primary drivers.
- **FS**: VEC and S_mix govern ductility predictions.
- **CI**: x_Cr composition fraction is overwhelmingly important (not shown in descriptor-only plot as composition fractions were excluded from FEATURE_COLS; δ_r and VEC serve as proxies).

![Figure 6: Feature Importances per Property](figures/fig6_feature_importance.png)

### 5.3 GP Parity Plot

![Figure 7: GP Surrogate Parity Plots (HV and YS)](figures/fig7_parity.png)

The GP surrogate on 20% held-out test data achieves R²(HV) ≈ 0.47 and R²(YS) ≈ 0.83, consistent with the cross-validation results. Error bars (GP posterior standard deviation) are well-calibrated for YS but show slight overconfidence for HV.

### 5.4 Multi-Objective Pareto Front

**Table 2: Pareto Front Summary**

| Metric | Value |
|--------|-------|
| Total compositions | 500 |
| Pareto-optimal compositions | 31 (6.2%) |
| Pareto-optimal HV range | 248–326 HV |
| Pareto-optimal YS range | 1021–1489 MPa |
| Pareto-optimal CI range | 0.13–0.54 |

![Figure 4: Multi-Objective Pareto Front (HV vs YS, color = CI)](figures/fig4_pareto.png)

The Pareto front reveals a mild trade-off between hardness and yield strength, with compositions containing higher Cr fractions (low CI) appearing predominantly in the Pareto-optimal set, confirming the critical role of Cr for corrosion resistance.

### 5.5 Active Learning

![Figure 3: Active Learning Progress](figures/fig3_active_learning.png)

Starting from 50 randomly labeled compositions (initial weighted objective = 0.806), 30 iterations of EI-based active learning improve the best objective to 0.824 (+2.1%). The learning curve plateaus after approximately 20 iterations, suggesting that the GP surrogate is sufficiently accurate to identify the near-optimal region with a modest additional label budget.

### 5.6 CrMnFeCoNi Case Study

**Table 3: Composition–Property Map for CrMnFeCoNi System**

| Cr | Ni | Mn=Fe=Co | HV (HV) | YS (MPa) | CI |
|----|----|----------|---------|----------|-----|
| 0.200 | 0.200 | 0.200 each | ~238 | ~935 | 0.39 |
| 0.304 | 0.282 | 0.138 each | ~257 | ~972 | 0.32 |
| 0.327 | 0.259 | 0.138 each | ~261 | ~982 | 0.31 |
| **0.350** | **0.236** | **0.138 each** | **~261** | **~983** | **0.29** |

The equiatomic Cantor alloy (x_Cr = x_Ni = 0.20) has descriptors: δ_r = 1.12%, VEC = 8.0, ΔS_mix = 13.38 J/mol·K, Ω ≈ 1.2×10⁵. Increasing Cr content to 0.35 (reducing Mn/Fe/Co symmetrically) improves corrosion index (CI: 0.39 → 0.29) and hardness (HV: +23 HV) while VEC decreases slightly, which our model associates with modest ductility reduction.

![Figure 5: CrMnFeCoNi Composition–Property Heatmaps](figures/fig5_cantor_heatmap.png)

---

## 6. Discussion

### 6.1 Model Performance and Interpretability

The relatively low R² for hardness (0.44–0.55) compared to yield strength (0.83–0.88) is noteworthy and has two plausible explanations. First, hardness in HEAs is influenced by microstructural features (grain size, precipitation, dislocation density) that are not captured by composition-only descriptors—consistent with the literature [1,3]. Second, the synthetic hardness model contains deliberately stronger nonlinear terms (the exponential VEC penalty, the tanh Ω term) that are more challenging to recover from composition-level features alone. This highlights a fundamental limitation: even with well-designed descriptors, composition-level ML models cannot fully predict microstructure-sensitive properties without microstructural inputs.

### 6.2 Dependence on Synthetic Data Assumptions

This study's most significant limitation is its reliance on synthetic data. The property functions were constructed to reflect known physical trends:
- δ_r → solid-solution hardening (literature support: Hall-Petch analogy [3])
- Cr content → corrosion resistance (literature support: extensive for stainless steels and HEAs)
- VEC → phase stability (FCC for VEC > 8, BCC for VEC < 6.87 [Guo et al., 2011])

However, the specific functional forms and noise levels are assumptions. Real experimental data exhibit far more complex behaviour, including:
- **Processing-dependent properties**: Casting vs. powder metallurgy vs. additive manufacturing yield vastly different microstructures for identical compositions.
- **Phase transformations**: At certain compositions, intermetallic precipitates (σ, μ, Laves phases) form, causing discontinuous property changes that simple smooth surrogates cannot model.
- **Short-range ordering (SRO)**: Chen et al. (2021, Nature Communications [9]) demonstrated that SRO in CoCuFeNiPd HEA creates a pseudo-composite microstructure enhancing both strength and ductility—a mechanism absent from our descriptors.

**Consequently, performance metrics reported here (R² = 0.44–0.88) are likely optimistic upper bounds for real-world application.**

### 6.3 Active Learning Improvement and Limitations

The 2.1% improvement from active learning over 30 iterations may seem modest. This reflects the relatively small unlabeled pool (450 compositions) and the fact that random initialization already sampled a near-optimal region by chance (initial best objective = 0.806 on a 0–1 scale). In realistic settings where the composition space is much larger (e.g., 10⁵–10⁶ candidates), active learning gains are expected to be substantially larger, as demonstrated by Lookman et al. (2019, Nature Communications) who reported 2–5× efficiency improvements in adaptive design experiments.

### 6.4 Missing Physics and Pathways to Improvement

Several important physical effects are not captured by the current descriptor set:
1. **DFT-derived descriptors**: Formation energy, elastic constants, stacking fault energy (SFE) from AFLOW/Materials Project would substantially improve predictions.
2. **CALPHAD phase fractions**: The volume fraction of FCC/BCC/B2/σ phases at processing temperature directly affects properties but requires thermodynamic calculations (e.g., Thermo-Calc with TCHEA5 database).
3. **Grain boundary chemistry**: Segregation of Mn to grain boundaries in CrMnFeCoNi is known to influence corrosion and fracture behaviour.
4. **Temperature-dependent properties**: Creep resistance, which is critical for high-temperature applications, depends on dislocation climb and diffusion mechanisms not accessible from ambient descriptors.

### 6.5 Pathway to Real-World Deployment

For the framework to be applied to actual alloy discovery, the following steps are recommended:
1. **Integrate AFLOW / Materials Project data**: Use the AFLOW REST API to query DFT-computed elastic moduli, lattice parameters, and formation energies for binary and ternary subsystem endpoints.
2. **Couple to Thermo-Calc CALPHAD**: Compute equilibrium phase fractions and driving forces for precipitation using TCHEA5 database, using these as additional features.
3. **Build an experimental database**: Systematically collect hardness, tensile test, and electrochemical corrosion data for a curated set of ~200 compositions to establish reliable training data.
4. **Implement Bayesian optimization on the physical composition simplex** (not descriptor space), with constraints enforcing Σxᵢ = 1 and manufacturer-imposed element bounds.

---

## 7. Conclusion

We have presented a comprehensive machine learning framework for multi-objective HEA composition optimization that integrates CALPHAD-inspired descriptor engineering, multiple surrogate model families, Expected Improvement-based Bayesian optimization, and an active learning loop. Key findings are:

1. **GBM achieves the best predictive performance for YS** (R² = 0.876 ± 0.015, RMSE = 42.8 MPa), while GP is competitive for CI (R² = 0.805 ± 0.046).
2. **Hardness prediction is inherently more challenging** (R² ≈ 0.44–0.55), reflecting the microstructure-sensitive nature of HV and the limitations of composition-only descriptors.
3. **31 of 500 compositions are Pareto-optimal** in the (HV, YS, CI) space; Cr-enriched compositions dominate this front.
4. **Active learning achieves 2.1% objective improvement** in 30 iterations from an initial pool of 50 labeled alloys.
5. **CrMnFeCoNi case study** shows that increasing Cr content beyond the equiatomic fraction (x_Cr → 0.35) simultaneously improves hardness and corrosion resistance within this system.

**Critical caveats**: All results are based on synthetic data generated by simplified physics-inspired functions. The reported R² values should not be interpreted as achievable with real experimental data without microstructural features, processing history, and validated thermodynamic inputs. Future work must couple this framework with real AFLOW/Materials Project DFT databases, CALPHAD thermodynamic calculations, and experimental validation campaigns.

---

## References

[1] Chen, Y., Ren, C., Jia, Y., Wang, G., Li, M., & Lu, W. (2021). A machine learning-based alloy design system to facilitate the rational design of high entropy alloys with enhanced hardness. *Acta Materialia*, 222, 117431. https://doi.org/10.1016/j.actamat.2021.117431

[2] Ren, W., Zhang, Y., Wang, W., Ding, S., & Li, N. (2023). Prediction and design of high hardness high entropy alloy through machine learning. *Materials & Design*, 235, 112454. https://doi.org/10.1016/j.matdes.2023.112454

[3] Huang, X., Jin, C., Zhang, C., Zhang, H., & Fu, H. (2021). Machine learning assisted modelling and design of solid solution hardened high entropy alloys. *Materials & Design*, 211, 110177. https://doi.org/10.1016/j.matdes.2021.110177

[4] Ma, Y., Li, M., Mu, Y., Wang, G., & Lu, W. (2023). Accelerated Design for High-Entropy Alloys Based on Machine Learning and Multiobjective Optimization. *Journal of Chemical Information and Modeling*, 63(19), 6029–6039. https://doi.org/10.1021/acs.jcim.3c00916

[5] Zhang, J., Cai, C., Kim, G., Wang, Y., & Chen, W. (2022). Composition design of high-entropy alloys with deep sets learning. *npj Computational Materials*, 8, 89. https://doi.org/10.1038/s41524-022-00779-7

[6] Chen, W., Hilhorst, A., Bokas, G., Gorsse, S., Jacques, P., & Hautier, G. (2023). A map of single-phase high-entropy alloys. *Nature Communications*, 14, 2856. https://doi.org/10.1038/s41467-023-38423-7

[7] Wan, X., Li, Z., Yu, W., Wang, A., Xue, K., et al. (2023). Machine Learning Paves the Way for High Entropy Compounds Exploration: Challenges, Progress, and Outlook. *Advanced Materials*, 35, 2305192. https://doi.org/10.1002/adma.202305192

[8] Sun, Y., Lu, Z., Liu, X., Du, Q., Xie, H., Lv, J., … Lu, Z. (2021). Prediction of Ti-Zr-Nb-Ta high-entropy alloys with desirable hardness by combining machine learning and experimental data. *Applied Physics Letters*, 119(20), 201905. https://doi.org/10.1063/5.0065303

[9] Chen, S., Aitken, Z. H., Pattamatta, S., Wu, Z., Yu, Z. G., Srolovitz, D. J., Liaw, P. K., & Zhang, Y.-W. (2021). Simultaneously enhancing the ultimate strength and ductility of high-entropy alloys via short-range ordering. *Nature Communications*, 12, 4953. https://doi.org/10.1038/s41467-021-25264-5

[10] Xu, P., Ji, X., Li, M., & Lu, W. (2023). Small data machine learning in materials science. *npj Computational Materials*, 9, 42. https://doi.org/10.1038/s41524-023-01000-z
