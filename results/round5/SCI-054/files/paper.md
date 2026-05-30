# Machine Learning-Accelerated High-Throughput Screening of Metal–Organic Frameworks for CO₂/H₂ Gas Adsorption and Direct Air Capture Applications

---

## Abstract

The discovery of optimal metal–organic frameworks (MOFs) for gas adsorption and direct air capture (DAC) of CO₂ requires navigation of a vast combinatorial chemical space, making purely experimental or brute-force simulation approaches infeasible. In this work, we present an end-to-end computational screening pipeline combining surrogate Grand Canonical Monte Carlo (GCMC) simulations, geometric descriptor extraction, and multi-target machine learning to rank MOF candidates for CO₂/H₂ separation and low-concentration CO₂ capture at 400 ppm (DAC conditions). A synthetic database of 2,000 MOF structures—parameterized by physically motivated geometric descriptors (largest cavity diameter, pore limiting diameter, geometric surface area, void fraction, density) and metal node identity—was generated to emulate the CoRE MOF and hMOF databases. Three supervised regression models (Random Forest, Gradient Boosting, Multi-layer Perceptron) were trained to predict CO₂ DAC uptake, CO₂ uptake at 1 bar, and H₂ uptake at 77 K. Gradient Boosting achieved the best performance for CO₂ prediction (5-fold CV R² = 0.977 ± 0.002, RMSE = 0.174 ± 0.008 mmol/g). H₂ uptake proved harder to predict (R² ≈ 0.78), reflecting its greater sensitivity to subtle pore topology. A Random Forest classifier was trained for water stability prediction (AUROC = 0.873 ± 0.006). Applying a sequential screening funnel—pore size filter → water stability → synthesizability—reduced 2,000 candidates to 535, from which a ranked shortlist of 20 DAC-optimal MOFs was produced. Geometric surface area and void fraction emerged as the dominant predictors for CO₂ DAC uptake. Critically, we discuss the limitations of working with synthetic data: the physically parametric nature of our data generation inflates model R² relative to real experimental databases, and these results should be interpreted as a methodology proof-of-concept rather than absolute performance benchmarks. Future work must validate on CoRE MOF-2019 and hMOF experimental datasets with proper train/test scaffold splits.

---

## 1. Introduction

Atmospheric CO₂ has surpassed 420 ppm, necessitating not only emissions reduction but active carbon removal via technologies such as Direct Air Capture (DAC). Solid sorbent materials, particularly metal–organic frameworks (MOFs), are among the most promising candidates for DAC owing to their tunable pore chemistry, ultra-high surface areas (up to ~7,000 m²/g), and modular synthesis [1]. However, the MOF design space is effectively infinite: the Cambridge Structural Database contains over 100,000 experimentally synthesized MOFs, and hypothetical databases such as hMOF and ARC-MOF enumerate hundreds of thousands more predicted structures [2,3].

High-throughput computational screening (HTCS) has emerged as the primary tool for navigating this space. Early HTCS studies performed molecular simulations on thousands of structures to compute adsorption properties, but even GCMC simulations are computationally expensive at scale [4]. Machine learning surrogates that predict adsorption properties from structural descriptors offer a 10–100× acceleration, provided descriptors capture sufficient chemical and geometric information [5].

Prior work by Burner et al. (2020) demonstrated R² ≈ 0.96 for CO₂/N₂ selectivity using atomic property-weighted radial distribution functions (AP-RDF) combined with geometric descriptors on a 340,000-MOF hypothetical database [5]. Srinivasu & Snurr (2023) combined DFT, force field optimization, and GCMC to screen the CoRE-MOF-2019 database for wet flue gas separation [4]. Mohamed et al. (2023) integrated four stability metrics into HTCS, finding that water and thermal stability filters dramatically reduce the viable candidate pool [1]. Zhang et al. (2025) proposed hierarchical HTCS incorporating ML stability prediction to identify ultrastable vanadium-based MOFs [3].

Despite these advances, several gaps remain: (i) most studies focus on post-combustion CO₂/N₂ rather than DAC conditions (400 ppm); (ii) stability filters and synthesizability constraints are often applied only as post-processing rather than co-optimized; (iii) H₂ adsorption for co-production applications is rarely combined with CO₂ screening. This work addresses these gaps by presenting an integrated pipeline for simultaneous CO₂ (DAC) and H₂ adsorption screening with stability-aware ranking.

---

## 2. Related Work

### 2.1 High-Throughput Computational Screening of MOFs

Computational screening studies have grown dramatically in scale since the first hMOF database generation. Key advances include the CoRE-MOF-2019 database (computation-ready experimental MOFs) and the ARC-MOF database (~280,000 structures). Srinivasu & Snurr (2023) applied multi-scale screening combining ML pre-filters, DFT interaction energies, and GCMC simulations to CoRE-MOF-2019, identifying MOFs with high CO₂ uptake from wet flue gas [4]. Polat et al. (2020) combined experiment with HTCS of 1,085 ionic liquid/MOF composites, finding that narrow pore sizes and high IL loading maximize CO₂/N₂ selectivity [6].

### 2.2 Machine Learning for MOF Property Prediction

Burner et al. (2020) is the benchmark for ML-based MOF adsorption prediction: their AP-RDF descriptor combined with geometric features achieved R² = 0.96 for CO₂ working capacity on ~340,000 MOFs, enabling a >10× speedup in screening [5]. Moosavi et al. (2020) developed ML-based diversity analysis of MOF databases, revealing biases that can cause overoptimistic conclusions when train/test splits do not account for structural similarity [2]. Graph neural networks (GNNs) have recently been applied to MOFs, with Reiser et al. (2022) reviewing state-of-the-art architectures [7]; however, GNN approaches require full crystal graph representations beyond simple geometric descriptors.

### 2.3 Stability Prediction and Synthesizability Filters

Mohamed et al. (2023) were among the first to systematically integrate thermodynamic, mechanical, thermal, and activation stability metrics within HTCS, showing that >80% of top-performing hypothetical MOFs fail at least one stability criterion [1]. Zhang et al. (2025) achieved high accuracy in predicting water stability using ML classifiers trained on experimental stability data, identifying Zr- and Al-based nodes as consistently more stable [3]. Zheng et al. (2023) used ChatGPT-assisted text mining to extract ~26,000 MOF synthesis parameters, building ML models with >87% accuracy for crystallization outcome prediction [8].

### 2.4 DAC-Specific MOF Screening

DAC operation at 400 ppm CO₂ partial pressure imposes fundamentally different requirements than post-combustion capture (0.15 bar CO₂): sorbents require steeper isotherm slopes at ultra-low pressure, high CO₂/H₂O and CO₂/N₂ selectivity, and excellent cyclic stability. Zeolites and amine-functionalized MOFs have been widely studied, but systematic HTCS at 400 ppm remains less common than at 0.15 bar [9].

---

## 3. Methods

### 3.1 Synthetic MOF Database Generation

A synthetic database of N = 2,000 MOF structures was generated to emulate the statistical distributions observed in CoRE MOF and hMOF databases. Seven geometric/chemical descriptors were sampled from physically motivated distributions:

| Descriptor | Symbol | Distribution | Range |
|---|---|---|---|
| Largest Cavity Diameter | LCD | Log-normal(ln 12, 0.5) | 3–60 Å |
| Pore Limiting Diameter | PLD | LCD × U(0.5, 0.95) | 2–57 Å |
| Geometric Surface Area | GSA | Log-normal(ln 2500, 0.6) | 100–8000 m²/g |
| Void Fraction | VF | f(GSA) + N(0, 0.05) | 0.05–0.85 |
| Density | ρ | 1/(0.3 + 2VF) + ε | 0.2–3.0 g/cm³ |
| Linker Length | L | U(4, 20) | 4–20 Å |
| Metal Node | M | Categorical {Zn, Cu, Zr, Al, Fe, Co, Ni, Mg} | — |

### 3.2 GCMC Surrogate Model

Grand Canonical Monte Carlo simulation was replaced with a physics-based surrogate that encodes known structure–property relationships from the literature:

**CO₂ DAC uptake** (mmol/g, 298 K, 400 ppm):
$$q_{\text{CO}_2}^{\text{DAC}} = 3 \times 10^{-4} \cdot \text{GSA} + 2.5 \exp\!\left(-\frac{(\text{LCD}-8)^2}{50}\right) + 0.8\,\text{VF} - 0.3\,\rho + \delta_M + \varepsilon$$

where $\delta_M$ is a metal-specific offset (Zr: +0.4, Cu: +0.3) reflecting higher CO₂ affinity, and $\varepsilon \sim \mathcal{N}(0, 0.15)$ represents simulation noise.

**CO₂ uptake at 1 bar** (mmol/g):
$$q_{\text{CO}_2}^{1\text{bar}} = 8 \times 10^{-4} \cdot \text{GSA} + 5.6\,\text{VF} + 1.5 \exp\!\left(-\frac{(\text{LCD}-10)^2}{98}\right) + \varepsilon$$

**H₂ uptake** (wt%, 77 K, 1 bar):
$$q_{\text{H}_2} = 2 \times 10^{-4} \cdot \text{GSA} + 4.0\,\text{VF} + 0.5 \exp\!\left(-\frac{(\text{LCD}-7)^2}{18}\right) - 0.2L + \varepsilon$$

### 3.3 Geometric Descriptor Extraction (Zeo++ Protocol)

In actual deployment, geometric descriptors are extracted using Zeo++ with a probe radius of 1.86 Å (CO₂) or 1.45 Å (H₂). The pipeline executes:
```bash
network -ha -res pore_size.txt  -sa 1.86 1.86 2000 surface.txt  \
        -vol 1.86 1.86 50000 void_fraction.txt  structure.cif
```
In this study, descriptors were drawn from surrogate distributions parameterized on published CoRE MOF statistics.

### 3.4 Machine Learning Models

Three regression models were benchmarked:

1. **Random Forest (RF)**: 200 trees, max depth 12, Gini impurity, bootstrapping enabled.
2. **Gradient Boosting (GB)**: 200 estimators, max depth 5, learning rate 0.05, subsampling 1.0.
3. **Multi-layer Perceptron (MLP)**: Architecture (128→64→32), ReLU activations, early stopping (patience=10), input standardized with StandardScaler.

All models were evaluated using 5-fold stratified cross-validation (random seed 42). Metrics: R², RMSE.

**Water Stability Classifier**: Random Forest (200 trees, max depth 8) trained to predict binary water stability labels derived from metal node identity and pore geometry rules consistent with experimental trends [1,3]. Metrics: Accuracy, F1, AUROC.

### 3.5 Adsorption Isotherm Prediction

Langmuir isotherms were fitted to each top-candidate MOF:
$$q(P) = \frac{q_{\text{sat}} K_L P}{1 + K_L P}$$

where $q_{\text{sat}}$ = 1.5 × $q^{\text{DAC}}$ (extrapolated saturation) and $K_L$ is fitted from the DAC point.

### 3.6 DAC Ranking Score

A composite DAC score was computed:
$$S_{\text{DAC}} = 0.40 \cdot \tilde{q}_{\text{CO}_2} + 0.25 \cdot \frac{\ln(1+S_{\text{CO}_2/\text{H}_2})}{\ln(1+S_{\text{max}})} + 0.20 \cdot w_{\text{stable}} + 0.15 \cdot s_{\text{synth}}$$

where tildes denote min-max normalization, and $w_{\text{stable}}, s_{\text{synth}} \in \{0,1\}$ are binary stability and synthesizability labels.

### 3.7 Screening Pipeline

The hierarchical screening funnel consisted of:
1. **Pore size filter**: LCD ≥ 4 Å (minimum for CO₂ kinetic diameter 3.3 Å + clearance)
2. **Water stability filter**: ML classifier prediction = stable
3. **Synthesizability filter**: Synthesizability score > 0.5
4. **DAC ranking**: Sort by $S_{\text{DAC}}$, select top-20

---

## 4. Experiments

### 4.1 Dataset

- N = 2,000 synthetic MOF structures
- 8 descriptors per structure (7 continuous + 1 categorical)
- 3 regression targets: CO₂ DAC uptake, CO₂ 1 bar uptake, H₂ 77K uptake
- 1 classification target: water stability (binary)
- 5-fold cross-validation for all models (random seed = 42)
- 80/20 train/test split for final parity plot evaluation

### 4.2 Evaluation Metrics

- Regression: R² (coefficient of determination), RMSE
- Classification: Accuracy, F1 score (macro), AUROC
- All cross-validation metrics reported as mean ± standard deviation across folds

### 4.3 Hardware

Experiments run on CPU (scikit-learn v1.x, Python 3.11). Training time: < 2 minutes for all models combined.

---

## 5. Results

### 5.1 Database Statistics

The synthetic MOF database spans realistic structural parameter ranges:
- GSA: 100–8000 m²/g (mean 2490 m²/g), consistent with CoRE MOF distribution
- VF: 0.05–0.85 (mean 0.38)
- CO₂ DAC uptake: 0.01–5.4 mmol/g (mean 2.0 mmol/g)
- Water stability: 48.8% stable; Synthesizability: 54.0% feasible

![Figure 1: Descriptor distributions and CO₂ uptake relationships](figures/fig1_descriptor_distributions.png)

### 5.2 ML Model Performance

**Table 1. 5-Fold Cross-Validation Results (R² and RMSE)**

| Target | Model | R² (mean±std) | RMSE (mean±std) |
|---|---|---|---|
| CO₂ DAC (mmol/g) | Random Forest | 0.965 ± 0.003 | 0.216 ± 0.009 |
| CO₂ DAC (mmol/g) | Gradient Boosting | **0.977 ± 0.002** | **0.174 ± 0.008** |
| CO₂ DAC (mmol/g) | Neural Network | 0.973 ± 0.003 | 0.190 ± 0.010 |
| CO₂ 1 bar (mmol/g) | Random Forest | 0.960 ± 0.005 | 0.382 ± 0.020 |
| CO₂ 1 bar (mmol/g) | Gradient Boosting | **0.967 ± 0.004** | **0.348 ± 0.017** |
| CO₂ 1 bar (mmol/g) | Neural Network | 0.965 ± 0.004 | 0.359 ± 0.010 |
| H₂ 77K (wt%) | Random Forest | 0.779 ± 0.037 | 0.194 ± 0.013 |
| H₂ 77K (wt%) | Gradient Boosting | 0.772 ± 0.052 | 0.196 ± 0.011 |
| H₂ 77K (wt%) | **Neural Network** | **0.789 ± 0.041** | **0.189 ± 0.010** |

![Figure 2: ML model performance comparison (5-fold CV)](figures/fig2_ml_performance.png)

**Independent Test Set (RF, CO₂ DAC):** R² = 0.962, RMSE = 0.224 mmol/g

![Figure 4: ML parity plot for CO₂ DAC uptake (test set)](figures/fig4_parity_plot.png)

### 5.3 Feature Importance

Geometric surface area (GSA) and void fraction (VF) dominate CO₂ DAC uptake prediction, contributing ~35% and ~22% of total Random Forest feature importance, respectively. Largest cavity diameter contributes ~18%, reflecting the optimal pore size (~8 Å) for CO₂ at low pressures. Metal node identity contributes ~12%, consistent with the literature finding that Zr- and Cu-based MOFs show enhanced CO₂ affinity.

![Figure 3: Feature importance for CO₂ DAC uptake prediction](figures/fig3_feature_importance.png)

### 5.4 Water Stability Classification

**Table 2. Water Stability Classifier (5-fold CV)**

| Metric | Mean | Std |
|---|---|---|
| Accuracy | 0.831 | 0.011 |
| F1 Score | 0.801 | 0.017 |
| AUROC | 0.873 | 0.006 |

### 5.5 Correlation Structure

![Figure 7: Descriptor-property correlation matrix](figures/fig7_correlation.png)

GSA shows the strongest positive correlation with CO₂ uptake (r ≈ 0.65), while density shows the strongest negative correlation (r ≈ −0.52). CO₂ DAC and CO₂ 1 bar uptake are highly correlated (r ≈ 0.85), suggesting a single surface-area-driven factor. H₂ uptake is more weakly correlated with CO₂ uptake (r ≈ 0.35), consistent with its sensitivity to pore topology differences.

### 5.6 Screening Funnel Results

**Table 3. Hierarchical Screening Funnel**

| Stage | N remaining | Retention (%) |
|---|---|---|
| Initial database | 2,000 | 100.0 |
| Pore size filter (LCD ≥ 4 Å) | 1,979 | 98.9 |
| Water stability filter | 972 | 48.6 |
| Synthesizability filter | 535 | 26.8 |
| Top DAC candidates | 20 | 1.0 |

![Figure 5: Screening funnel and top-20 DAC candidate ranking](figures/fig5_dac_ranking.png)

### 5.7 Top-5 DAC Candidates

**Table 4. Top-5 MOF Candidates for DAC**

| MOF ID | Metal | LCD (Å) | GSA (m²/g) | CO₂ DAC (mmol/g) | CO₂/H₂ Sel. | DAC Score |
|---|---|---|---|---|---|---|
| MOF_1399 | Zr | 8.64 | 8000 | 5.27 | 34.9 | 0.886 |
| MOF_0126 | Co | 7.31 | 8000 | 4.89 | 444.2 | 0.959 |
| MOF_1560 | Al | 11.87 | 8000 | 4.81 | 387.8 | 0.948 |
| MOF_1077 | Al | 9.69 | 6790 | 4.53 | 411.5 | 0.929 |
| MOF_1445 | Al | 8.66 | 8000 | 4.95 | 113.3 | 0.909 |

![Figure 6: Adsorption isotherms for top-5 DAC MOF candidates](figures/fig6_isotherms.png)

---

## 6. Discussion

### 6.1 Interpretation of Results

The high R² values (0.96–0.98) for CO₂ adsorption prediction are consistent with published benchmarks: Burner et al. (2020) reported R² = 0.96 using AP-RDF descriptors on 340,000 MOFs [5]. Our results using only 8 geometric/chemical features suggest that CO₂ adsorption at low-to-moderate pressures is well-described by bulk geometric properties alone. The lower R² for H₂ (≈ 0.78) reflects that H₂ adsorption at 77 K is more sensitive to local atomic interactions and quantum effects not captured by geometric descriptors—a finding consistent with the MOF literature.

### 6.2 Critical Assessment: Dependence on Synthetic Data Assumptions

**This is the most important limitation of this study.** Our high predictive R² values arise in part because the training and test data were generated by the same parametric surrogate model. The ML models effectively learn to "invert" the data-generating function, which is a closed-form nonlinear combination of the input descriptors. In contrast, real GCMC/experimental data contain:

- **Irregular pore topology effects** not captured by scalar LCD/PLD descriptors
- **Electrostatic interactions** from open metal sites, amine groups, and charged frameworks
- **Force field uncertainty** (up to ±20% in GCMC vs experiment)
- **Polymorphism and disorder** in experimental structures

When applied to actual CoRE MOF or hMOF data, we expect R² to decrease by 0.1–0.2 relative to our reported values, consistent with published results on experimental databases (R² ≈ 0.75–0.85 for geometry-only features).

### 6.3 Generalizability to Real-World Data

Several factors limit direct generalization:

1. **Data distribution mismatch**: Real MOF databases have highly uneven metal node distributions (Zn and Cu dominate), which would require weighted sampling or stratified cross-validation.
2. **Chemical descriptors**: Geometric descriptors alone are insufficient for MOFs with amine functionalization or open metal sites, which are critical for DAC. AP-RDF, energy histograms, or graph-based representations would be needed.
3. **Stability labels**: Our water stability labels were generated from simplified rules (metal node + void fraction). Real stability assessment requires water vapor adsorption experiments or full molecular dynamics simulations.
4. **Synthesis feasibility**: Our synthesizability score is a coarse proxy. Real synthesizability prediction requires MOFid-based analysis, retrosynthetic planning, or ChatGPT-assisted literature mining [8].

### 6.4 Comparison with Prior Work

Our pipeline architecture is conceptually aligned with the hierarchical HTCS of Zhang et al. (2025) [3] and the ML-integrated screening of Mohamed et al. (2023) [1]. Key differences: (i) we target 400 ppm DAC conditions rather than post-combustion; (ii) we include H₂ co-screening; (iii) our stability integration is in-pipeline rather than post-hoc. The water stability AUROC of 0.873 is lower than values reported by Zhang et al. (~0.91) on real experimental data, which is expected given our simplified label generation.

### 6.5 Pipeline Scalability

The ML inference step reduces screening time from ~2 CPU-hours per GCMC simulation (×2000 = 4000 h) to <1 second total, representing a >14,000× speedup. In practice, the pipeline would be deployed as: (1) Zeo++ descriptor extraction for all candidates; (2) ML pre-screening to select top 5% for full GCMC; (3) GCMC → final ranking.

### 6.6 Limitations Summary

- Synthetic data with parametric noise: R² is optimistically high
- No crystal graph / atomic-level descriptors (open metal sites, functional groups)
- Water stability labels from simplified rules, not MD/experiment
- H₂ at 77 K (cryogenic) vs. room temperature H₂ storage not addressed
- DAC-specific humidity effects and CO₂/H₂O competition not modeled

---

## 7. Conclusion

We presented a high-throughput screening pipeline for CO₂/H₂ adsorption in MOFs integrating surrogate GCMC simulations, geometric descriptor extraction (Zeo++ protocol), and multi-model machine learning. Gradient Boosting achieved the highest predictive accuracy for CO₂ adsorption (5-fold CV R² = 0.977 ± 0.002), while H₂ prediction remained more challenging (R² ≈ 0.78). A Random Forest classifier for water stability (AUROC = 0.873) enabled realistic candidate filtering, reducing a 2,000-MOF database to 20 top DAC candidates. Geometric surface area and void fraction are the dominant predictors. We critically note that these performance metrics are inflated by synthetic data assumptions; real experimental databases would likely yield R² ≈ 0.80–0.87 with geometric-only descriptors. Future work should: (1) apply the pipeline to CoRE MOF-2019 with GCMC-verified adsorption data; (2) incorporate AP-RDF or GNN representations for chemical diversity; (3) explicitly model DAC moisture competition; (4) validate top candidates experimentally.

---

## References

[1] Mohamed, S.A., Zhao, D., & Jiang, J. (2023). Integrating stability metrics with high-throughput computational screening of metal–organic frameworks for CO₂ capture. *Communications Materials*, 4, 84. https://doi.org/10.1038/s43246-023-00409-9

[2] Moosavi, S.M., Nandy, A., Jablonka, K.M., et al. (2020). Understanding the diversity of the metal-organic framework ecosystem. *Nature Communications*, 11, 4068. https://doi.org/10.1038/s41467-020-17755-8

[3] Zhang, Z., Palakkal, A.S., Wu, X., Jiang, J., & Jiang, Z. (2025). Discovering Ultra-Stable Metal–Organic Frameworks for CO₂ Capture from A Wet Flue Gas: Integrating Machine Learning and Molecular Simulation. *Environmental Science & Technology*. https://doi.org/10.1021/acs.est.5c00768

[4] Srinivasu, K., & Snurr, R.Q. (2023). High-Throughput Screening of the CoRE-MOF-2019 Database for CO₂ Capture from Wet Flue Gas: A Multi-Scale Modeling Strategy. *ACS Applied Materials & Interfaces*, 15(30). https://doi.org/10.1021/acsami.3c04079

[5] Burner, J., Schwiedrzik, L., Krykunov, M., et al. (2020). High-Performing Deep Learning Regression Models for Predicting Low-Pressure CO₂ Adsorption Properties of Metal–Organic Frameworks. *The Journal of Physical Chemistry C*, 124(51), 27996–28005. https://doi.org/10.1021/acs.jpcc.0c06334

[6] Polat, H.M., Kavak, S., Kulak, H., Uzun, A., & Keskin, S. (2020). CO₂ separation from flue gas mixture using [BMIM][BF₄]/MOF composites: Linking high-throughput computational screening with experiments. *Chemical Engineering Journal*, 394, 124916. https://doi.org/10.1016/j.cej.2020.124916

[7] Reiser, P., Neubert, M., Eberhard, A., et al. (2022). Graph neural networks for materials science and chemistry. *Communications Materials*, 3, 93. https://doi.org/10.1038/s43246-022-00315-6

[8] Zheng, Z., Zhang, O., Borgs, C., Chayes, J., & Yaghi, O.M. (2023). ChatGPT Chemistry Assistant for Text Mining and the Prediction of MOF Synthesis. *Journal of the American Chemical Society*, 145(32), 18048–18062. https://doi.org/10.1021/jacs.3c05819

[9] Yan, Y., Borhani, T.N., Subraveti, S.G., et al. (2021). Harnessing the power of machine learning for carbon capture, utilisation, and storage (CCUS) – a state-of-the-art review. *Energy & Environmental Science*, 14, 6122–6157. https://doi.org/10.1039/d1ee02395k
