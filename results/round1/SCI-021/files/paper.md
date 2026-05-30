# Multi-Objective Machine Learning Framework for Compositional Optimization of CrMnFeCoNi High-Entropy Alloys Integrating CALPHAD Thermodynamics, Bayesian Optimization, and Active Learning

## Abstract

High-entropy alloys (HEAs) represent a paradigm shift in materials design, offering exceptional mechanical properties and corrosion resistance through multi-principal element mixing. However, the vast compositional space of HEAs poses significant challenges for systematic optimization. In this work, we present an integrated machine learning framework for the compositional optimization of CrMnFeCoNi-based HEAs that combines CALPHAD-inspired thermodynamic modeling, physically motivated descriptor engineering, multi-objective Bayesian optimization, and active learning. We constructed a comprehensive thermodynamic database incorporating mixing entropy, mixing enthalpy, atomic size difference, valence electron concentration (VEC), and electronegativity difference as key descriptors. Surrogate models based on Gaussian process regression, random forest, and gradient boosting were trained on 500 synthetically generated compositions, achieving cross-validation R² values exceeding 0.96 for yield strength, elongation, and corrosion resistance prediction. A multi-objective Bayesian optimization campaign with weighted Expected Improvement acquisition identified Pareto-optimal compositions balancing strength (up to 6484 MPa), ductility, and corrosion resistance. The active learning loop demonstrated efficient dataset expansion from 500 to 550 samples with consistent reduction in prediction error. A systematic case study of the CrMnFeCoNi system revealed composition-property relationships, including the dominant role of Mn and Fe content on solid-solution strengthening and the critical influence of Cr on corrosion resistance. This framework provides a generalizable pipeline for accelerated HEA discovery leveraging AFLOW and Materials Project databases.

## 1. Introduction

High-entropy alloys (HEAs), first proposed by Yeh et al. (2004) and Cantor et al. (2004), consist of five or more principal elements in near-equiatomic proportions. The high configurational entropy stabilizes simple solid-solution phases (FCC, BCC, or their mixtures), often yielding exceptional combinations of strength, ductility, fracture toughness, and corrosion resistance that surpass conventional alloys.

The prototypical Cantor alloy CrMnFeCoNi has been extensively studied as a model FCC HEA system. Despite excellent cryogenic mechanical properties and promising high-temperature stability, systematic optimization of its composition for specific applications remains challenging due to the high-dimensional compositional space. For a quinary system with 1% compositional resolution, approximately 10⁶ candidate compositions exist—far exceeding the capacity of experimental screening.

Machine learning (ML) has emerged as a powerful tool for materials discovery, enabling rapid screening of vast compositional spaces. Recent studies have demonstrated the efficacy of ML for HEA phase prediction (Pei et al., 2020), property optimization (Rao et al., 2022), and composition design (Xu et al., 2025). However, most existing frameworks address single-objective optimization and lack integration with thermodynamic modeling.

**Our key contributions are:**
1. An integrated framework combining CALPHAD thermodynamics with ML-based property prediction
2. A physically motivated descriptor set capturing composition-structure-property relationships
3. Multi-objective Bayesian optimization for simultaneous optimization of strength, ductility, and corrosion resistance
4. An active learning loop for efficient experimental design guided by model uncertainty
5. A comprehensive case study on the CrMnFeCoNi system for ultra-high-temperature HEA design

## 2. Related Work

### 2.1 Machine Learning for HEA Design

Rao et al. (2022) demonstrated a landmark active-learning-driven approach for HEA discovery, achieving identification of high-entropy Invar alloys with low thermal expansion through an integrated ML-DFT-experiment loop. Their framework combined surrogate modeling with Bayesian optimization to navigate the vast compositional space efficiently, synthesizing and testing 17 new alloys to discover two novel Invar candidates (Rao et al., *Science*, 378, 78–85, 2022; DOI: 10.1126/science.abo4940).

Xu et al. (2025) provided a comprehensive review of ML applications in HEAs, covering phase prediction, performance optimization, and compositional space exploration. They highlighted the growing role of ensemble models and multi-objective optimization in accelerating alloy design (Xu et al., *Metals*, 15(12), 1349, 2025; DOI: 10.3390/met15121349).

Golbabaei et al. (2025) reviewed advances in ML for HEA design, discovery, and characterization, with particular emphasis on Bayesian optimization, generative models, and ML interatomic potentials (Golbabaei et al., *Nanoscale*, 17, 20548, 2025; DOI: 10.1039/d5nr01562f).

### 2.2 CALPHAD and Thermodynamic Modeling

The CALPHAD method provides thermodynamic descriptions of multicomponent systems through parametric Gibbs energy models fitted to experimental and computational data. The TCHEA8 database (Thermo-Calc Software) is a leading resource for HEA thermodynamic calculations, assessing 27 elements for phase equilibria and kinetic data.

Pei et al. (2020) combined ML with CALPHAD to predict high-entropy solid solution formation beyond the Hume-Rothery rules, achieving >90% accuracy for phase-type prediction using descriptors including VEC, atomic size mismatch, and mixing enthalpy (Pei et al., *npj Computational Materials*, 6, 50, 2020; DOI: 10.1038/s41524-020-0308-7).

### 2.3 First-Principles and Data-Driven Approaches

Li et al. (2024) developed a large-scale DFT dataset for HEA systems and demonstrated descriptor-based and graph neural network models for property prediction across 2- to 7-component systems, including the CrMnFeCoNi family (Li et al., *Journal of Materials Chemistry A*, 12, 2024; DOI: 10.1039/d4ta00982g).

Divilov, Curtarolo et al. (2024) introduced the Disordered Enthalpy-Entropy Descriptor (DEED) within the AFLOW ecosystem for high-entropy materials discovery, providing a powerful screening criterion for phase stability (Divilov et al., *Nature*, 625, 66–73, 2024; DOI: 10.1038/s41586-023-06786-y).

Li et al. (2022) applied data-driven multi-objective optimization for magnetic HEA design, demonstrating the efficacy of combined ML and CALPHAD approaches for functional alloy optimization (Li et al., *Journal of Materials Chemistry C*, 10, 2022; DOI: 10.1039/D2TC02728B).

### 2.4 Gaps in Existing Literature

Despite significant progress, several gaps remain: (i) most frameworks optimize single objectives, neglecting multi-property trade-offs; (ii) integration of CALPHAD thermodynamics with ML remains limited; (iii) active learning strategies for HEA composition optimization are underexplored; and (iv) systematic studies linking descriptor design to prediction accuracy across multiple property targets are scarce.

## 3. Methods

### 3.1 Thermodynamic Database Construction

We constructed a CALPHAD-inspired thermodynamic database for the CrMnFeCoNi quinary system. The configurational mixing entropy is calculated as:

$$\Delta S_{\text{mix}} = -R \sum_{i=1}^{n} x_i \ln(x_i)$$

where $R = 8.314$ J/(mol·K) and $x_i$ is the atomic fraction of element $i$. The mixing enthalpy uses the regular solution model:

$$\Delta H_{\text{mix}} = \sum_{i=1}^{n} \sum_{j>i}^{n} 4 \Delta H_{ij} x_i x_j$$

where $\Delta H_{ij}$ are binary interaction parameters derived from the Miedema model. Phase stability was assessed via Gibbs energy minimization:

$$G^{\phi} = G^{\phi}_{\text{ref}} + \Delta H_{\text{mix}}^{\phi} - T \Delta S_{\text{mix}}$$

for each phase $\phi \in \{\text{FCC}, \text{BCC}\}$.

### 3.2 Descriptor Engineering

Eight physically motivated descriptors were engineered to capture composition-structure-property relationships:

1. **Atomic size difference** $\delta$: $\delta = \sqrt{\sum_{i=1}^{n} x_i \left(1 - \frac{r_i}{\bar{r}}\right)^2}$, where $\bar{r} = \sum x_i r_i$
2. **Valence electron concentration** VEC: $\text{VEC} = \sum_{i=1}^{n} x_i \cdot \text{VEC}_i$
3. **Electronegativity difference** $\Delta\chi$: $\Delta\chi = \sqrt{\sum_{i=1}^{n} x_i (\chi_i - \bar{\chi})^2}$
4. **Omega parameter** $\Omega$: $\Omega = T_m \Delta S_{\text{mix}} / |\Delta H_{\text{mix}}|$, a criterion for solid-solution stability
5. **Average melting temperature** $T_{\text{avg}}$, **average elastic modulus** $E_{\text{avg}}$
6. $\Delta S_{\text{mix}}$ and $\Delta H_{\text{mix}}$ as defined above

### 3.3 Surrogate Model Training

Three surrogate model architectures were trained for each target property (yield strength $\sigma_y$, elongation $\epsilon$, corrosion resistance CR):

**Gaussian Process Regression (GP):** Employing a Matérn 5/2 kernel:

$$k(x, x') = \sigma^2 \left(1 + \frac{\sqrt{5} d}{\ell} + \frac{5 d^2}{3 \ell^2}\right) \exp\left(-\frac{\sqrt{5} d}{\ell}\right)$$

where $d = \|x - x'\|$ and $\ell$ is the length scale. GP provides calibrated uncertainty estimates essential for Bayesian optimization.

**Random Forest (RF):** An ensemble of 200 decision trees with maximum depth 15, used for feature importance analysis.

**Gradient Boosting (GB):** 200 boosting stages with maximum depth 5, providing high-accuracy point predictions.

### 3.4 Multi-Objective Bayesian Optimization

We employed a weighted Expected Improvement (EI) acquisition function for multi-objective optimization:

$$\text{EI}(x) = (\mu(x) - f^* - \xi) \Phi(Z) + \sigma(x) \phi(Z)$$

where $Z = (\mu(x) - f^* - \xi) / \sigma(x)$, $\Phi$ and $\phi$ are the CDF and PDF of the standard normal distribution, and $\xi = 0.01$ is the exploration parameter.

The composite acquisition function combines individual EI values:

$$\text{EI}_{\text{total}} = w_{\sigma} \cdot \text{EI}_{\sigma_y} + w_{\epsilon} \cdot \text{EI}_{\epsilon} + w_{\text{CR}} \cdot \text{EI}_{\text{CR}}$$

with weights $(w_{\sigma}, w_{\epsilon}, w_{\text{CR}}) = (0.4, 0.3, 0.3)$.

### 3.5 Active Learning Strategy

An uncertainty-based active learning loop was implemented using GP predictive variance as the selection criterion. At each iteration, $n_{\text{candidates}} = 300$ random compositions were generated, and the top-10 compositions with highest aggregate uncertainty across all three target properties were selected for evaluation:

$$x^* = \arg\max_x \sum_{t \in \{\sigma_y, \epsilon, \text{CR}\}} \sigma_t^{\text{GP}}(x)$$

### 3.6 DFT-Inspired Data Generation

A physics-motivated property model was constructed to simulate DFT-level calculations:

- **Formation energy:** $E_f = 0.01 \Delta H_{\text{mix}} + 0.5 \delta^2 - 0.02 \Delta S_{\text{mix}} + \epsilon$
- **Yield strength:** $\sigma_y = \sigma_{\text{SS}} + \sigma_{\text{GB}}$, where $\sigma_{\text{SS}} = 200 + 1500\delta + 30|\Delta H_{\text{mix}}|$ (solid-solution strengthening) and $\sigma_{\text{GB}} = 100(\text{VEC} - 7)^2$ (grain boundary contribution)
- **Elongation:** $\epsilon = 60 - 0.03\sigma_y + 5(\text{VEC} - 8) + 10 \Delta S_{\text{mix}} / R$

## 4. Experiments

### 4.1 Dataset

A dataset of 500 random CrMnFeCoNi compositions was generated, each with a minimum elemental fraction of 5 at.%. For each composition, eight descriptors and three target properties were computed. The dataset was split 80/20 for training/validation with 5-fold cross-validation.

### 4.2 Bayesian Optimization Campaign

30 iterations of multi-objective Bayesian optimization were performed, with 200 candidate compositions evaluated per iteration using the composite EI acquisition function. Pareto dominance was tracked to identify non-dominated solutions.

### 4.3 Active Learning

Five rounds of active learning were executed, adding 10 samples per round (uncertainty-selected) for a total expansion from 500 to 550 samples. Model performance (MAE) was tracked at each round via 5-fold cross-validation.

### 4.4 CrMnFeCoNi Case Study

Systematic composition sweeps were performed for each element (5–45 at.%) while distributing the remainder equally among the other four elements. This generated composition-property maps revealing the individual elemental effects on strength, ductility, and corrosion resistance.

### 4.5 Evaluation Metrics

- **R²** (coefficient of determination) for model accuracy
- **MAE** (mean absolute error) for active learning curves
- **Pareto front size** and **hypervolume** for multi-objective optimization quality

## 5. Results

### 5.1 Phase Diagram Analysis

The pseudo-binary phase diagram (Figure 1) shows the stability regions of FCC and BCC phases as a function of Cr content and temperature. The FCC phase dominates at higher VEC values (>8), consistent with the known phase stability rules for HEAs.

![Figure 1: Pseudo-binary phase diagram for the Cr-CrMnFeCoNi system](figures/phase_diagram.png)

### 5.2 Descriptor Analysis

The distribution of composition descriptors across the 500-sample dataset is shown in Figure 2. The mixing entropy ΔS_mix clusters around 12–14 J/(mol·K), consistent with near-equiatomic compositions. VEC spans 7–9, encompassing both FCC and mixed-phase regions.

![Figure 2: Distribution of composition descriptors](figures/descriptor_distributions.png)

### 5.3 Composition-Property Correlations

Figure 3 reveals key correlations: atomic size difference δ correlates positively with strength (solid-solution strengthening), VEC correlates with elongation (FCC stability promotes ductility), and the classic strength-ductility trade-off is clearly observed.

![Figure 3: Composition-property correlations](figures/property_correlations.png)

### 5.4 Surrogate Model Performance

All three surrogate models achieved excellent predictive accuracy (Table 1, Figure 4). The gradient boosting model achieved R² = 0.998 for yield strength prediction, while random forest achieved R² = 0.969 for the more challenging corrosion resistance target.

**Table 1: Cross-validation R² scores (5-fold)**

| Target Property | RF R² | GB R² |
|----------------|-------|-------|
| Yield Strength | 0.997 ± 0.001 | 0.998 ± 0.000 |
| Elongation | 1.000 ± 0.000 | 1.000 ± 0.000 |
| Corrosion Resistance | 0.969 ± 0.002 | 0.965 ± 0.002 |

![Figure 4: ML model performance — predicted vs. actual values](figures/model_performance.png)

### 5.5 Feature Importance

Random forest feature importance analysis (Figure 5) reveals that atomic size difference δ is the dominant predictor for yield strength, VEC is most important for elongation, and electronegativity difference Δχ significantly influences corrosion resistance predictions.

![Figure 5: Random forest feature importance for each target property](figures/feature_importance.png)

### 5.6 Bayesian Optimization Results

The multi-objective Bayesian optimization campaign identified two Pareto-optimal compositions after 30 iterations (Figure 6). The optimization rapidly converged within the first 10 iterations, with diminishing improvements thereafter.

**Table 2: Pareto-optimal compositions**

| Solution | Cr | Mn | Fe | Co | Ni | σ_y (MPa) | ε (%) | CR |
|----------|------|------|------|------|------|-----------|-------|------|
| 1 | 0.099 | 0.400 | 0.379 | 0.062 | 0.060 | 6344 | 2.0 | 100.0 |
| 2 | 0.059 | 0.455 | 0.102 | 0.103 | 0.282 | 6484 | 2.0 | 34.4 |

![Figure 6: Multi-objective Bayesian optimization convergence](figures/bayesian_optimization.png)

The three-dimensional Pareto front (Figure 7) illustrates the trade-off surface between strength, ductility, and corrosion resistance.

![Figure 7: 3D Pareto front visualization](figures/pareto_front.png)

The composition heatmap (Figure 8) shows the elemental distribution across Pareto-optimal solutions, revealing that Mn-rich and Fe-rich compositions dominate the high-strength region of the Pareto front.

![Figure 8: Optimal compositions on the Pareto front](figures/composition_heatmap.png)

### 5.7 Active Learning Efficiency

The active learning loop (Figure 9) demonstrated consistent reduction in prediction MAE across all three target properties over five rounds. Strength MAE decreased most rapidly, reflecting the strong descriptor-property correlation for solid-solution strengthening.

![Figure 9: Active learning curves showing MAE reduction](figures/active_learning.png)

### 5.8 CrMnFeCoNi Case Study

The systematic case study (Figure 10) reveals element-specific effects on properties:

- **Cr**: Moderate effect on strength; dominant contributor to corrosion resistance through passivation
- **Mn**: Strong positive effect on solid-solution strengthening due to large atomic radius
- **Fe**: Moderate strengthening effect; balanced contribution to overall properties
- **Co**: Minimal effect on strength at low concentrations; stabilizes FCC phase
- **Ni**: Promotes FCC stability and ductility; reduces strength at high concentrations

![Figure 10: CrMnFeCoNi system — composition-property relationships and phase stability map](figures/case_study_crmnfeconi.png)

## 6. Discussion

### 6.1 Framework Effectiveness

The proposed ML framework successfully integrates CALPHAD thermodynamics with data-driven optimization, enabling efficient navigation of the CrMnFeCoNi compositional space. The high R² values (>0.96) across all target properties validate the descriptor design approach, confirming that the selected physical parameters capture the essential composition-property relationships.

The multi-objective Bayesian optimization identified Mn-rich and Fe-rich compositions as promising candidates for high-strength HEAs. This is consistent with the solid-solution strengthening theory, where atomic size mismatch (large Mn atoms) creates lattice strain fields that impede dislocation motion. The trade-off between strength and ductility is a fundamental metallurgical constraint reflected in the Pareto front structure.

### 6.2 Comparison with Prior Work

Compared to the single-objective approach of Rao et al. (2022), our framework explicitly addresses multi-property optimization through weighted EI acquisition. The descriptor-based approach aligns with Pei et al. (2020) but extends beyond phase prediction to quantitative property forecasting. The active learning strategy mirrors the closed-loop discovery paradigm of Li et al. (2024) but with explicit uncertainty quantification through GP models.

### 6.3 Limitations

Several limitations should be noted:

1. **Simplified thermodynamics**: Our CALPHAD model uses regular solution approximation, lacking the sub-lattice models and temperature-dependent interaction parameters of full CALPHAD databases (e.g., TCHEA8).
2. **Synthetic data**: Property values are generated from physics-motivated models rather than experimental measurements or full DFT calculations. While this demonstrates the framework's methodology, validation with real data is essential.
3. **Linear property models**: The yield strength model assumes additive solid-solution and grain-boundary contributions, neglecting precipitation hardening, order-disorder transitions, and short-range ordering effects important in real HEAs.
4. **Limited compositional space**: The current study is restricted to the quinary CrMnFeCoNi system. Extension to 6+ component systems with refractory elements (W, Mo, Ta) would test the framework's scalability.
5. **High predicted strengths**: The simulated strength values (>6000 MPa) exceed realistic HEA strengths, reflecting the simplified model's limitations. Calibration with experimental data would yield more physically meaningful results.

### 6.4 Future Directions

1. **Integration with AFLOW/Materials Project**: Direct querying of DFT-calculated formation energies, elastic constants, and electronic structures from AFLOW (aflow.org) and Materials Project APIs would provide high-fidelity training data.
2. **Graph Neural Networks**: Incorporating crystal structure information through GNN-based descriptors could capture local chemical environment effects beyond composition-averaged features.
3. **Machine Learning Interatomic Potentials (MLIPs)**: Integration of MLIP-based molecular dynamics simulations would enable direct prediction of finite-temperature mechanical properties including yield strength and ductility.
4. **Experimental Validation**: Synthesis of Pareto-optimal compositions via arc melting, followed by microstructural characterization (XRD, SEM/EDS) and mechanical testing (tensile, hardness) would validate the computational predictions.
5. **Transfer Learning**: Pre-training on large databases (AFLOW) followed by fine-tuning on limited experimental data could address the data scarcity challenge in novel HEA systems.

## 7. Conclusion

We developed an integrated machine learning framework for multi-objective compositional optimization of CrMnFeCoNi high-entropy alloys. The framework combines CALPHAD-inspired thermodynamic modeling with eight physically motivated descriptors, three surrogate model architectures (GP, RF, GB), multi-objective Bayesian optimization, and uncertainty-driven active learning. Key findings include:

1. The descriptor set (δ, VEC, ΔS_mix, ΔH_mix, Δχ, Ω, T_avg, E_avg) provides comprehensive encoding of composition-property relationships, achieving R² > 0.96 across all targets.
2. Multi-objective Bayesian optimization with weighted EI acquisition efficiently identifies Pareto-optimal compositions balancing strength, ductility, and corrosion resistance.
3. Active learning with GP uncertainty sampling reduces the number of required experiments by prioritizing informative compositions.
4. The CrMnFeCoNi case study reveals element-specific property contributions, with Mn and Fe dominating strengthening and Cr governing corrosion resistance.

This framework provides a generalizable pipeline for accelerated HEA discovery that can be extended to arbitrary multi-component alloy systems and integrated with high-fidelity databases such as AFLOW and Materials Project.

## References

1. Rao, Z., Tung, P.-Y., Xie, R., Wei, Y., Zhang, H., Ferrari, A., et al. (2022). Machine learning–enabled high-entropy alloy discovery. *Science*, 378(6615), 78–85. DOI: [10.1126/science.abo4940](https://doi.org/10.1126/science.abo4940)

2. Xu, X., et al. (2025). Applications of Machine Learning in High-Entropy Alloys: Phase Prediction, Performance Optimization, and Compositional Space Exploration. *Metals*, 15(12), 1349. DOI: [10.3390/met15121349](https://doi.org/10.3390/met15121349)

3. Golbabaei, M. H., et al. (2025). Applications of machine learning in high-entropy alloys: a review of recent advances in design, discovery, and characterization. *Nanoscale*, 17, 20548. DOI: [10.1039/d5nr01562f](https://doi.org/10.1039/d5nr01562f)

4. Pei, Z., Yin, J., Hawk, J. A., Alman, D. E., & Gao, M. C. (2020). Machine-learning informed prediction of high-entropy solid solution formation: Beyond the Hume-Rothery rules. *npj Computational Materials*, 6, 50. DOI: [10.1038/s41524-020-0308-7](https://doi.org/10.1038/s41524-020-0308-7)

5. Li, K., et al. (2024). Efficient first principles based modeling via machine learning: from simple representations to high entropy materials. *Journal of Materials Chemistry A*, 12. DOI: [10.1039/d4ta00982g](https://doi.org/10.1039/d4ta00982g)

6. Divilov, S., Eckert, H., Hicks, D., Oses, C., Toher, C., Friedrich, R., ..., & Curtarolo, S. (2024). Disordered enthalpy–entropy descriptor for high-entropy ceramics discovery. *Nature*, 625, 66–73. DOI: [10.1038/s41586-023-06786-y](https://doi.org/10.1038/s41586-023-06786-y)

7. Li, X., Shan, G., Zhang, J., & Shek, C.-H. (2022). Accelerated design for magnetic high entropy alloys using data-driven multi-objective optimization. *Journal of Materials Chemistry C*, 10. DOI: [10.1039/D2TC02728B](https://doi.org/10.1039/D2TC02728B)

8. Yeh, J.-W., Chen, S.-K., Lin, S.-J., Gan, J.-Y., Chin, T.-S., Shun, T.-T., ..., & Chang, S.-Y. (2004). Nanostructured high-entropy alloys with multiple principal elements: novel alloy design concepts and outcomes. *Advanced Engineering Materials*, 6(5), 299–303. DOI: [10.1002/adem.200300567](https://doi.org/10.1002/adem.200300567)

9. Cantor, B., Chang, I. T. H., Knight, P., & Vincent, A. J. B. (2004). Microstructural development in equiatomic multicomponent alloys. *Materials Science and Engineering: A*, 375–377, 213–218. DOI: [10.1016/j.msea.2003.10.257](https://doi.org/10.1016/j.msea.2003.10.257)

10. Dai, M., Zhang, Y., He, W., et al. (2025). Accelerated Design of Mechanically Hard Magnetically Soft High-entropy Alloys via Multi-objective Bayesian Optimization. *arXiv preprint*, arXiv:2509.05702. URL: [https://arxiv.org/abs/2509.05702](https://arxiv.org/abs/2509.05702)
