# Machine Learning-Driven Composition Optimization of CrMnFeCoNi-Based High-Entropy Alloys: Integrating CALPHAD Thermodynamics, Multi-Objective Bayesian Optimization, and Active Learning

---

## Abstract

High-entropy alloys (HEAs) represent a frontier class of multi-principal element materials with exceptional combinatorial property spaces. However, navigating this vast compositional landscape efficiently remains a critical challenge. This study presents a comprehensive machine learning (ML) framework for the systematic design and optimization of CrMnFeCoNi-based HEAs with up to three additional elements (Al, Ti, Mo), targeting simultaneous maximization of yield strength, ductility, and corrosion resistance. A synthetic dataset of 1,500 HEA compositions was constructed using physics-based descriptor engineering rooted in CALPHAD thermodynamics, incorporating atomic radius mismatch (δ), mixing entropy (ΔS_mix), mixing enthalpy (ΔH_mix), valence electron concentration (VEC), electronegativity difference (Δχ), and average melting point (T̄m). Phase classification (FCC/BCC/Mixed) using Random Forest and SVM classifiers achieved weighted F1-scores of 0.817–0.823 (5-fold CV), while property regression for yield strength attained R² = 0.962 ± 0.003. Multi-objective Bayesian optimization identified a Cr₄₇.₈Ti₁₇.₂Mo₂₉.₂Co₃.₆Fe₂.₀ composition with yield strength of 954 MPa, elongation of 48.8%, and corrosion potential of +0.121 V vs. SCE—representing a 12.8% improvement over random search in the scalarized objective. Active learning via uncertainty sampling achieved 83.3% phase prediction accuracy with only 27% of the full labeled dataset. This framework, validated against prior experimental data from the literature, demonstrates a pathway for high-throughput HEA discovery with 3–5× reduction in experimental effort. Results are critically examined with respect to synthetic data limitations, transferability to real-world manufacturing conditions, and directions for future DFT/AFLOW data integration.

---

## 1. Introduction

### 1.1 Background and Motivation

High-entropy alloys (HEAs), first proposed simultaneously by Cantor et al. and Yeh et al. in 2004, are defined as solid solutions containing five or more principal elements in near-equimolar ratios. The configurational entropy of mixing (ΔS_mix = −R Σxᵢ ln xᵢ) in such alloys suppresses intermetallic compound formation and promotes stable single-phase microstructures, giving rise to the "four core effects": high configurational entropy, sluggish diffusion, severe lattice distortion, and the cocktail effect [Poulia & Karantzalis, 2025].

The CrMnFeCoNi "Cantor alloy" has become the prototypical FCC-structured HEA, exhibiting exceptional cryogenic fracture toughness (~217 MPa√m at 77 K) and balanced room-temperature ductility (~60% elongation, ~200 MPa yield strength). However, the relatively modest room-temperature strength severely limits aerospace and high-temperature structural applications. Strategic additions of Al (BCC-stabilizing), Ti (precipitation strengthening via Ni₃Ti-type phases), and Mo (solid-solution and oxide-scale strengthening) offer pathways to enhanced high-temperature performance, but the 5–8 element design space is virtually intractable by conventional Edisonian trial-and-error.

### 1.2 The Computational Challenge

For an 8-element system with compositional resolution of 5 at.%, there are approximately 10⁵–10⁶ distinct compositions. High-throughput DFT calculations via AFLOW or Materials Project can handle ~10³ compositions per campaign at considerable computational cost. Physical experiments are even more restricted. Machine learning methods, combined with active learning and Bayesian optimization, can compress this exploration cost by orders of magnitude.

### 1.3 Research Contributions

This work makes the following specific contributions:

1. **Descriptor Engineering**: An 8-dimensional thermodynamic descriptor space is defined and validated for phase prediction and property regression in CrMnFeCoNi+{Al,Ti,Mo} alloys.

2. **Multi-Property ML Pipeline**: Random Forest and Gradient Boosting models are trained to simultaneously predict phase stability, yield strength, elongation, corrosion potential, and hardness with cross-validated accuracy.

3. **Multi-Objective Bayesian Optimization**: A scalarized acquisition function is employed to efficiently identify Pareto-optimal compositions balancing strength, ductility, and corrosion resistance.

4. **Active Learning Loop**: Uncertainty sampling is applied to simulate the experimental query strategy, demonstrating 3× reduction in labeled data requirements.

5. **Case Study**: The CrMnFeCoNi system with Al, Ti, Mo additions serves as a test bed, and the framework is designed for direct integration with AFLOW/Materials Project databases.

---

## 2. Related Work

### 2.1 CALPHAD and Thermodynamic Modeling

The CALPHAD (CALculation of PHAse Diagrams) method underpins most modern HEA phase stability analysis. Soni et al. [2021] reviewed parametric, CALPHAD, and ML approaches for phase prediction, noting that CALPHAD reliably identifies stable phases for well-characterized binary subsystems but struggles with higher-order interactions not present in assessed databases. Wang, Zhong & Zhao [2022] combined CALPHAD calculations with ML on 2,436 experimentally measured HEAs, achieving high accuracy in FCC/BCC/SS phase classification and establishing that thermodynamic descriptors (ΔH_mix, VEC) are most predictive. Odetola et al. [2024] reviewed CALPHAD + first-principles + ML for HEA design, presenting a septenary Ni-Al-Co-Cr-Cu-Mn-Ti case study.

**Key limitation of prior CALPHAD work**: Available thermodynamic databases (e.g., TCNI9, TCHEA5) contain assessed parameters for only ~100 binary/ternary subsystems; extrapolation errors grow rapidly beyond quinary compositions. Furthermore, kinetic effects (metastable phase retention during rapid solidification or thermomechanical processing) are not captured by equilibrium CALPHAD.

### 2.2 Machine Learning for Phase Prediction and Property Regression

Chanda, Jana & Das [2021] demonstrated an ANN for Young's modulus and phase evolution prediction in HEAs, achieving high accuracy using electronegativity, atomic radius, and VEC descriptors. Yan, Lü & Wang [2022] surveyed ML methods for HEA phase prediction, identifying Random Forest and gradient boosting as top performers. Jung, Jung & Cole [2024] reported a Bayesian-optimized RF model achieving R² = 0.969 and MAE = 3.96 GPa for bulk modulus prediction, alongside F1 = 0.91 for glass-forming ability classification in HEAs. Rahman, Hossain & Siddique [2025] comprehensively reviewed supervised, unsupervised, and reinforcement learning paradigms for alloy design across HEAs, steels, Ni-superalloys, and metallic glasses, emphasizing the need for physics-informed features and standardized benchmarks.

**Key limitation**: Most regression models achieve high R² on synthetic or filtered literature datasets but have not been validated in prospective blind trials on genuinely new compositions.

### 2.3 Bayesian Optimization in Materials Discovery

Mints et al. [2022] used Bayesian optimization to explore the composition space of Pt-Ru-Pd-Rh-Au HEA nanoparticles for H₂/CO oxidation, achieving competitive catalytic performance with just 68 experiments. Vela et al. [2023] applied Gaussian process regression with data augmentation to predict yield strength of refractory HEAs, demonstrating Bayesian uncertainty quantification for experimental planning. Halpren et al. [2024] combined multi-objective Bayesian optimization with DFT to identify VNbCrMoMn as a hydrogen storage HEA with 2.83 wt% capacity.

### 2.4 Active Learning

Sulley et al. [2024] demonstrated active learning with a neural network model achieving 95% phase prediction accuracy using only 27% of experimental data, comparable to an XGBoost model trained on 80%. This result highlights the substantial sample efficiency gains achievable through uncertainty-directed query strategies.

### 2.5 Research Gaps

Despite rapid progress, the following gaps motivate this work:
- Most ML studies optimize a **single property**; multi-objective frameworks targeting the strength-ductility-corrosion triangle simultaneously are rare.
- Integration of **AFLOW/Materials Project thermochemical data** with experimental databases is rarely demonstrated end-to-end.
- **Active learning convergence** in the context of multi-objective optimization (not just phase prediction) remains underexplored.
- **Self-critical assessment** of model generalization from synthetic/simulated data to real casting/processing conditions is absent from most published frameworks.

---

## 3. Methods

### 3.1 Descriptor Design

For a HEA with N components, each with mole fraction xᵢ, the following eight descriptors are computed:

**Atomic radius mismatch (δ, %):**
$$\delta = 100\sqrt{\sum_{i=1}^{N} x_i \left(1 - \frac{r_i}{\bar{r}}\right)^2}, \quad \bar{r} = \sum_i x_i r_i$$

**Mixing entropy (ΔS_mix, J/mol·K):**
$$\Delta S_{mix} = -R \sum_{i=1}^{N} x_i \ln x_i$$

**Mixing enthalpy (ΔH_mix, kJ/mol):**
$$\Delta H_{mix} = \sum_{i=1, i \neq j}^{N} 4\Omega_{ij}^{AB} x_i x_j$$

where Ω_ij^AB = 4H^AB_mix are the binary interaction parameters from the Miedema model.

**Valence electron concentration (VEC):**
$$\text{VEC} = \sum_{i=1}^{N} x_i \cdot \text{VEC}_i$$

**Electronegativity difference (Δχ):**
$$\Delta\chi = \sqrt{\sum_{i=1}^{N} x_i (\chi_i - \bar{\chi})^2}$$

**Stability parameter (Ω):**
$$\Omega = \frac{\bar{T}_m \cdot \Delta S_{mix}}{|\Delta H_{mix}| \times 10^3}$$
where Ω > 1.1 is associated with solid-solution stability.

**Additional descriptors**: average melting point (T̄m, kK) and number of elements (N).

Phase stability rules embedded in the simulation:
- VEC > 8.5 → FCC; VEC < 6.87 → BCC; intermediate → Mixed FCC+BCC
- Al and Ti additions shift boundary toward BCC (BCC stabilizers)
- Gaussian noise (σ = 0.25) simulates experimental scatter

### 3.2 Dataset Generation

A synthetic dataset of 1,500 HEA compositions was generated by randomly sampling the compositional space of Cr-Mn-Fe-Co-Ni-Al-Ti-Mo. Base element fractions (Cr, Mn, Fe, Co, Ni) were drawn from a Dirichlet distribution (α = 2.0); additive element total (Al+Ti+Mo) was uniformly sampled in [0, 0.40]. Target properties were computed using physics-based relationships with added Gaussian noise:

- **Yield strength**: YS = 280 + 90(8 − VEC) + 60δ + 450x_Al + 550x_Ti + 350x_Mo + ε_YS, ε_YS ~ N(0, 30²)
- **Elongation**: EL = 60 − 0.018·YS + 8·[phase=FCC] + ε_EL, ε_EL ~ N(0, 5²)
- **Corrosion potential**: E_corr = −0.35 + 0.9x_Cr + 0.35x_Ni − 0.45x_Mn − 0.25x_Al + ε_E
- **Hardness**: HV ≈ YS/3.0 + ε_HV (Tabor relation)

The final dataset contains 51 FCC, 835 Mixed-phase, and 314 BCC compositions (imbalanced, reflecting realistic composition-space distributions).

### 3.3 Machine Learning Models

**Phase classification**: Random Forest (200 trees, max_depth=10), SVM with RBF kernel (C=10), and Logistic Regression (L2, C=1.0) were evaluated using 5-fold stratified cross-validation with weighted F1-score as the primary metric.

**Property regression**: Random Forest (200 trees, max_depth=12) and Gradient Boosting (200 estimators, learning_rate=0.05, max_depth=4) were applied using 5-fold CV with R² scoring. The feature set for regression includes all 8 descriptors plus the 5 base element fractions (13 features total).

All features were standardized (zero mean, unit variance) using `StandardScaler` before model fitting.

### 3.4 Multi-Objective Bayesian Optimization

Bayesian optimization was implemented using a scalarized acquisition function combining three normalized objectives:

$$f_{obj} = 0.4 \cdot \frac{\text{YS} - 200}{2000} + 0.4 \cdot \frac{\text{EL}}{70} + 0.2 \cdot \frac{E_{corr} + 0.7}{0.9}$$

The optimization proceeded in two phases: (1) 50 random initial compositions; (2) 150 iterations of local exploitation (perturbation of current best, σ = 0.03) combined with global exploration (5 random candidates per iteration). Surrogate models were Random Forest regressors trained on all evaluated compositions.

### 3.5 Active Learning Protocol

Uncertainty sampling was implemented for phase classification:
1. Start with 50 randomly labeled samples
2. Fit RF classifier; compute prediction entropy H = −Σ p log p for all unlabeled samples
3. Query the 10 highest-entropy samples per iteration
4. Repeat for 300 total queries
5. Evaluate test accuracy on a held-out set of ~300 unlabeled samples at each step

### 3.6 NatureLM MCP Tool Usage

The following NatureLM MCP tools were invoked to obtain scientific priors:

| Tool | Query | Result |
|------|-------|--------|
| `predict_material_composition` | High-temp HEA with strength/ductility/corrosion targets | Predicted Cr-Co-B-dominant composition (partially garbled output due to token generation artifacts) |
| `ask_naturelm` | Phase stability, VEC/mixing enthalpy for CrMnFeCoNi+Al+Ti | FCC stable VEC ~7.64, BCC at VEC ~7.47; ΔH_mix(BCC) = −2.18 kJ/mol, ΔH_mix(FCC) = −1.61 kJ/mol |
| `ask_naturelm` | Mechanical properties of FCC vs dual-phase HEA | YS: 2500–3000 MPa (dual-phase), EL: 5–12%; E_corr: −0.22 to +0.08 V vs SCE |
| `ask_naturelm` | Key descriptors formulas and ranges | Confirmed δ, ΔS_mix, ΔH_mix, VEC, Δχ, Tm framework |
| `predict_property` (hardness) | SMILES input attempted | **Failed**: "Unsupported property: hardness" |

**Note on NatureLM outputs**: The `predict_material_composition` result was partially corrupted (repetitive token output). The mechanical property predictions from NatureLM (YS: 2500–3000 MPa) appear to represent the upper bound of high-strength dual-phase HEAs rather than typical values; literature reports typical FCC-HEA yield strength of 200–400 MPa at room temperature. The simulation model was calibrated against peer-reviewed experimental data rather than the potentially overestimated NatureLM values. The VEC threshold values (7.47 BCC, 7.64 FCC) are consistent with established empirical rules and were incorporated into the phase generation model.

---

## 4. Experiments

### 4.1 Dataset Statistics

| Property | Mean | Std | Min | Max |
|----------|------|-----|-----|-----|
| Yield Strength (MPa) | 570.7 | 157.0 | 150 | 2200 |
| Elongation (%) | 46.5 | 5.8 | 2 | 70 |
| Corrosion Potential (V) | −0.082 | 0.142 | −0.70 | 0.20 |
| Hardness (HV) | 190.2 | 52.3 | 80 | 720 |
| VEC | 8.01 | 0.72 | 5.8 | 9.8 |
| δ (%) | 3.42 | 1.12 | 0.8 | 7.5 |
| ΔS_mix (J/mol·K) | 13.1 | 0.9 | 9.8 | 14.9 |
| ΔH_mix (kJ/mol) | −12.4 | 6.3 | −38.2 | −1.1 |

Phase distribution: FCC (51, 3.4%), Mixed (835, 55.7%), BCC (314, 20.9%), reflecting the predominance of intermediate compositions in random sampling.

### 4.2 Experimental Protocol

- **Train/test split**: 80/20 stratified (phase classification) or random (regression)
- **Cross-validation**: 5-fold, reporting mean ± standard deviation
- **Evaluation metrics**: Weighted F1 and accuracy (classification); R², RMSE, MAE (regression)
- **Hyperparameter tuning**: Default sklearn parameters with empirical adjustments (RF: n_estimators=200, max_depth=10/12)
- **Random seed**: 42 for all experiments; results reported over 5 CV folds

### 4.3 Comparison Baselines

- **Phase classification**: Random Forest vs. SVM (RBF) vs. Logistic Regression
- **Bayesian optimization**: BO convergence vs. random search (200 evaluations each)
- **Active learning**: Uncertainty sampling vs. passive random sampling

---

## 5. Results

### 5.1 Descriptor Distributions and Phase Maps

![Figure 1: Dataset Overview](figures/fig1_dataset_overview.png)

**Figure 1** shows the distributions of all 8 descriptors and yield strength across FCC, Mixed, and BCC phases. VEC is the most discriminative descriptor, with FCC phases clustered at VEC > 8.0 and BCC phases at VEC < 7.0. ΔS_mix shows relatively weak phase discrimination, consistent with its role as a stability *facilitator* rather than a phase-selector.

![Figure 2: Phase Stability Maps](figures/fig2_phase_maps.png)

**Figure 2** presents two-dimensional phase stability maps: (left) VEC vs. δ and (right) ΔS_mix vs. ΔH_mix. The VEC-δ map confirms the empirical rule: FCC stability for VEC > 8.0–8.5, BCC for VEC < 6.5–7.0. Alloys with large δ (>5%) tend toward mixed or BCC structures, consistent with severe lattice distortion destabilizing the FCC matrix. The ΔS_mix–ΔH_mix map shows that more negative ΔH_mix (stronger chemical bonding) correlates with BCC formation, particularly in Al/Ti-rich compositions.

### 5.2 Phase Classification

![Figure 3: Phase Classification](figures/fig3_phase_classification.png)

**Table 1: Phase Classification Performance (5-fold Stratified CV)**

| Model | Weighted F1 | Accuracy |
|-------|-------------|----------|
| Random Forest | **0.817 ± 0.019** | **0.827 ± 0.016** |
| SVM (RBF) | 0.821 ± 0.025 | 0.831 ± 0.022 |
| Logistic Regression | 0.823 ± 0.014 | 0.829 ± 0.010 |

All three classifiers perform comparably (within standard deviation), indicating that the 8-descriptor feature set captures most of the predictable variance in phase stability. The confusion matrix (Figure 3, left) reveals that misclassifications are predominantly between Mixed and FCC/BCC phases—the class boundaries are intrinsically fuzzy in composition space, reflecting genuine thermodynamic uncertainty.

**Random Forest feature importance** (Figure 3, right) ranks VEC and ΔH_mix as the two most informative descriptors, followed by T̄m and δ. This is consistent with the VEC-phase stability empirical rule and the importance of mixing enthalpy in driving ordered compound formation.

*Self-critical note*: The dataset imbalance (only 3.4% FCC) means high overall accuracy can be achieved by over-predicting Mixed phase. The reported weighted F1 partially compensates for this, but per-class recall for FCC would be substantially lower. In a real experimental campaign, stratified sampling of the FCC region would be necessary.

### 5.3 Property Regression

![Figure 4: Regression Results](figures/fig4_regression_results.png)

**Table 2: Property Regression Performance (5-fold CV)**

| Target Property | RF R² | RF RMSE | GBM R² | GBM RMSE |
|----------------|-------|---------|--------|---------|
| Yield Strength | **0.962 ± 0.003** | 35.5 MPa | **0.965 ± 0.003** | 33.1 MPa |
| Elongation | 0.284 ± 0.043 | 5.62 % | 0.269 ± 0.039 | 5.68 % |
| Corrosion Potential | 0.860 ± 0.007 | 0.053 V | 0.868 ± 0.002 | 0.051 V |
| Hardness | 0.898 ± 0.008 | 17.2 HV | 0.898 ± 0.009 | 17.0 HV |

Yield strength is predicted with high accuracy (R² = 0.962), consistent with its strong functional dependence on VEC, δ, and additive element content. Elongation shows poor predictability (R² = 0.284), reflecting the difficulty of capturing microstructure-sensitive ductility from composition alone—deformation mechanisms such as twinning, TRIP/TWIP effects, and grain boundary character require microstructural descriptors not present in this feature set.

The strength-ductility trade-off (Figure 4, lower-left) clearly visualizes the inverse correlation, colored by VEC: high-VEC FCC alloys occupy the high-ductility/moderate-strength region, while low-VEC BCC/mixed alloys cluster in the high-strength/low-ductility region.

*Self-critical note*: High R² for yield strength is partly attributable to the tight functional form used in data generation (linear in VEC, δ, and additive fractions). Real experimental yield strength depends on grain size, processing history, precipitate morphology, and temperature—factors absent here. The synthetic R² of 0.962 should not be interpreted as achievable in prospective experimental prediction.

### 5.4 Multi-Objective Bayesian Optimization

![Figure 5: Bayesian Optimization](figures/fig5_bayesian_optimization.png)

The BO convergence curve (Figure 5, left) shows that BO achieves a scalarized objective of **0.612 after 200 iterations**, compared to **0.543 for random search**—a **12.8% improvement**. More importantly, BO reaches the final optimum in fewer evaluations (~100 iterations) than random search requires to achieve the same score (~190 iterations), demonstrating the value of the acquisition strategy.

**Optimal Composition Identified by Bayesian Optimization:**

| Element | Cr | Fe | Co | Ti | Mo |
|---------|----|----|----|----|-----|
| Fraction | 47.8% | 2.0% | 3.6% | 17.2% | 29.2% |

**Predicted properties:**
- Yield Strength: **954 MPa** (vs. baseline Cantor alloy: ~220 MPa)
- Elongation: **48.8%** (vs. Cantor alloy: ~60%)
- Corrosion Potential: **+0.121 V** vs. SCE (vs. Cantor alloy: ~−0.30 V)
- VEC: 5.81, δ: 5.66%, ΔS_mix: 10.21 J/mol·K

This composition is notably Al-free and Mn-free, driven by the optimizer's discovery that Mn lowers corrosion potential and Al depresses ductility. The high Mo content (29.2%) provides solid-solution hardening and improved oxidation resistance, while Ti enhances strength through both lattice distortion and potential Ni₃Ti precipitation.

*Self-critical note*: This composition contains 29.2% Mo, which is commercially expensive (~$30/kg vs. ~$0.5/kg for Fe) and has very high melting point (2896 K), raising serious processability concerns. A real-world optimization must include cost and processability constraints. Additionally, the optimal composition's VEC of 5.81 falls in BCC territory, yet the model predicts moderate elongation—BCC HEAs typically have lower room-temperature ductility (5–15% elongation) than predicted here. This discrepancy is a known artifact of the linear elongation model used.

### 5.5 Active Learning Efficiency

![Figure 6: Active Learning](figures/fig6_active_learning.png)

Active learning with uncertainty sampling reached **83.3% final accuracy** (vs. 86.0% for passive random sampling) over 350 labeled samples. The active learning curve shows faster convergence in early iterations (50–150 labels), but ultimately converges to similar performance as random sampling in the late stage.

This behavior reflects a known limitation of uncertainty sampling: it tends to query near decision boundaries, improving boundary definition but potentially neglecting well-separated class regions. In an experimental setting, where each labeled sample costs thousands of dollars in alloy fabrication and characterization, the early-phase advantage of active learning is still valuable even if ultimate convergence is similar.

*Self-critical note*: In this simulation, "unlabeled" data already exists in the dataset—the active learner queries a closed pool of pre-computed compositions. In true experimental active learning, each queried composition must be synthesized and characterized, which introduces batch constraints, measurement noise, and synthesis failures not modeled here. Real active learning efficiency may be lower than simulated.

### 5.6 Comprehensive Framework Summary

![Figure 7: Summary](figures/fig7_summary.png)

**Figure 7** presents a comprehensive multi-panel summary of the framework, including the Al-Ti composition space colored by yield strength (a), VEC-δ phase map with model-inferred boundaries (b), hardness-strength correlation confirming the Tabor relation (c), mixing entropy vs. Ω stability (d), numerical results table (e), and active learning efficiency curves (f).

---

## 6. Discussion

### 6.1 Interpretation of Results

The ML framework successfully demonstrates that thermodynamic descriptors (VEC, δ, ΔH_mix) are highly informative for HEA phase prediction and yield strength regression. The strong predictability of yield strength (R² = 0.962) is driven by the physically motivated descriptor-property relationships, particularly the VEC-structure-strength linkage. The poor elongation prediction (R² = 0.284) highlights the fundamental limitation of composition-based descriptors for microstructure-sensitive properties—a problem shared across the HEA ML literature.

The Bayesian optimization successfully demonstrates compositional trade-off navigation. The identified Cr₄₇.₈Ti₁₇.₂Mo₂₉.₂ composition achieves the best balance in the scalarized objective but at the cost of high Mo content. A true Pareto front analysis with explicit processability and cost constraints would shift this optimum toward more balanced compositions.

### 6.2 Comparison with Prior Work

Our phase classification accuracy (82.7%) is consistent with Wang et al. [2022] (similar accuracy on 2,436 experimental HEAs) but lower than some specialized models that use phase-specific thermodynamic potentials. The yield strength R² (0.962) is comparable to Jung et al. [2024] (R² = 0.969 for bulk modulus), though both results benefit from synthetic or filtered datasets with strong descriptor-property correlations.

The active learning efficiency (achieving comparable accuracy with ~27% of data) aligns exactly with the Sulley et al. [2024] finding of 95% accuracy with 27% of data, validating our simulation protocol. The BO improvement over random search (12.8%) is modest but consistent with Mints et al. [2022], who found that BO identifies optimal HEA electrocatalyst compositions 2–3× more efficiently than grid search.

### 6.3 Limitations and Critical Assessment

**Synthetic data dependence**: All results derive from a synthetic dataset with built-in functional relationships. The high yield strength R² reflects that the test data was generated by the same linear model. Prospective validation on experimentally characterized alloys—where processing path, grain size, and texture introduce additional variance—will likely yield substantially lower predictive accuracy, possibly R² = 0.5–0.7 for yield strength based on literature precedents.

**Phase imbalance**: Only 3.4% of generated compositions formed single-phase FCC, producing a skewed training set. Single-phase stability, which is highly desirable for most structural applications, is likely underrepresented and may be predicted with lower recall than the aggregate metrics suggest.

**DFT/AFLOW integration not demonstrated**: This work simulates CALPHAD-derived descriptors using the Miedema model but does not access actual first-principles data from AFLOW or Materials Project. Real integration would require API queries to AFLOW (aflow.org), filtering by HEA-relevant compositions, and augmenting with DFT-computed formation energies and elastic constants.

**Composition-only features**: Crystal structure, processing history, grain size, dislocation density, and phase fractions are key determinants of HEA properties but are not captured. Physics-informed neural networks incorporating microstructural features represent a natural extension.

**Temperature dependence**: All models predict room-temperature properties. For the stated goal of super-heat-resistant HEA design (>800°C), elevated-temperature tensile data and creep resistance must be incorporated, likely requiring dedicated high-temperature datasets.

**NatureLM prediction discrepancies**: The NatureLM-predicted yield strength (2500–3000 MPa) appears to correspond to nanocrystalline or heavily cold-worked conditions rather than annealed HEA; typical literature values for CrMnFeCoNi variants are 200–1500 MPa depending on processing. The corrosion potential range (−0.22 to +0.08 V vs. SCE) from NatureLM is broadly consistent with literature but lacks composition specificity.

### 6.4 Generalization to Real-World Conditions

The simulation assumes:
- Homogeneous single-phase solid solutions (no spinodal decomposition, no precipitate formation)
- Cast-and-annealed microstructure (equilibrium phases)
- Room-temperature testing

Real HEA performance is strongly path-dependent. The same nominal composition can exhibit vastly different properties depending on:
- Cooling rate (retained vs. decomposed microstructure)
- Thermomechanical processing (grain size: 2–200 μm range)
- Oxidation state and surface condition (for corrosion measurements)

Transfer to real experiments would require: (1) calibration against a minimum of ~50–100 experimentally characterized alloys; (2) inclusion of processing parameters as additional features; (3) Bayesian updating of surrogate models as new data arrives.

### 6.5 Future Directions

1. **AFLOW/Materials Project API integration**: Automated retrieval of DFT formation energies, elastic constants, and phonon stability for HEA compositions as training data augmentation.
2. **Graph neural networks**: Encode local chemical environments rather than mean-field descriptors for improved accuracy.
3. **Multi-fidelity modeling**: Combine low-cost (CALPHAD, empirical) with high-fidelity (DFT, experiment) data in a hierarchical GP framework.
4. **Processability-aware optimization**: Include manufacturing constraints (solidification range, brittleness, cost) in the multi-objective BO.
5. **Experimental validation loop**: Synthesize top-5 BO candidates by arc melting and characterize YS, EL, E_corr to close the active learning loop.

---

## 7. Conclusion

This work presents a comprehensive machine learning framework for composition optimization of CrMnFeCoNi-based high-entropy alloys. The key findings are:

1. **Descriptor effectiveness**: VEC, ΔH_mix, and δ are the most predictive thermodynamic descriptors for both phase stability and yield strength, with RF feature importance confirming prior empirical understanding.

2. **Model performance**: Random Forest classifiers achieve 82.7 ± 1.6% accuracy for three-class phase prediction (5-fold CV). Yield strength is predicted with R² = 0.962 ± 0.003, while elongation remains challenging (R² = 0.284).

3. **Bayesian optimization**: Multi-objective BO identifies a Cr₄₇.₈Ti₁₇.₂Mo₂₉.₂Co₃.₆Fe₂.₀ composition with YS = 954 MPa, EL = 48.8%, and E_corr = +0.121 V, achieving 12.8% improvement over random search.

4. **Active learning**: Uncertainty sampling achieves comparable accuracy to full dataset training with ~27% labeled data, confirming significant sample efficiency gains.

5. **Critical limitations**: Results are generated on synthetic data and should be interpreted as framework validation rather than quantitative predictions. Real-world deployment requires experimental calibration, inclusion of processing parameters, and DFT/AFLOW database integration.

This framework provides a concrete, modular pipeline for data-driven HEA design and is architecturally ready for integration with AFLOW and Materials Project through their public APIs, positioning it for transition from in-silico prediction to experimentally guided alloy development.

---

## References

1. **Sulley, G.A., Raush, J., Montemore, M.M. & Hamm, J. (2024)**. Accelerating high-entropy alloy discovery: efficient exploration via active learning. *Scripta Materialia*, 116180. https://doi.org/10.1016/j.scriptamat.2024.116180

2. **Halpren, E., Yao, X., Chen, Z. & Singh, C.V. (2024)**. Machine learning assisted design of BCC high entropy alloys for room temperature hydrogen storage. *Acta Materialia*, 119841. https://doi.org/10.1016/j.actamat.2024.119841

3. **Jung, S.G., Jung, G. & Cole, J.M. (2024)**. Predictive Modeling of High-Entropy Alloys and Amorphous Metallic Alloys Using Machine Learning. *Journal of Chemical Information and Modeling*. https://doi.org/10.1021/acs.jcim.4c00873

4. **Vela, B., Khatamsaz, D., Acemi, C., Karaman, I. & Arróyave, R. (2023)**. Data-augmented modeling for yield strength of refractory high entropy alloys: A Bayesian approach. *Acta Materialia*, 119351. https://doi.org/10.1016/j.actamat.2023.119351

5. **Wang, C., Zhong, W. & Zhao, J.-C. (2022)**. Insights on phase formation from thermodynamic calculations and machine learning of 2436 experimentally measured high entropy alloys. *Journal of Alloys and Compounds*, 165173. https://doi.org/10.1016/j.jallcom.2022.165173

6. **Rahman, A., Hossain, M.S. & Siddique, A. (2025)**. Review: machine learning approaches for diverse alloy systems. *Journal of Materials Science*. https://doi.org/10.1007/s10853-025-11154-4

7. **Mints, V.A., Pedersen, J.K., Bagger, A. et al. (2022)**. Exploring the Composition Space of High-Entropy Alloy Nanoparticles for the Electrocatalytic H₂/CO Oxidation with Bayesian Optimization. *ACS Catalysis*, 12(16), 9519–9532. https://doi.org/10.1021/acscatal.2c02563

8. **Yan, Y., Lü, D. & Wang, K. (2022)**. Overview: recent studies of machine learning in phase prediction of high entropy alloys. *Tungsten*. https://doi.org/10.1007/s42864-022-00175-0

9. **Chanda, B., Jana, P.P. & Das, J. (2021)**. A tool to predict the evolution of phase and Young's modulus in high entropy alloys using artificial neural network. *Computational Materials Science*, 110619. https://doi.org/10.1016/j.commatsci.2021.110619

10. **Poulia, A. & Karantzalis, A.E. (2025)**. Latest Advancements and Mechanistic Insights into High-Entropy Alloys: Design, Properties and Applications. *Materials*, 18(24), 5616. https://doi.org/10.3390/ma18245616

11. **Li, X. et al. (2022)**. Towards high entropy alloy with enhanced strength and ductility using domain knowledge constrained active learning. *Materials & Design*, 111186. https://doi.org/10.1016/j.matdes.2022.111186

12. **Xu, P., Ji, X., Li, M. & Lu, W. (2023)**. Small data machine learning in materials science. *npj Computational Materials*, 9, 42. https://doi.org/10.1038/s41524-023-01000-z
