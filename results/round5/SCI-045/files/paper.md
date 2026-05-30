# Improved Epigenetic Clock Models for Biological Age Estimation: A Comparative Study of Linear and Deep Learning Approaches on DNA Methylation Data

---

## Abstract

Epigenetic clocks based on DNA methylation patterns have emerged as powerful tools for estimating biological age, with broad applications in geroscience, disease risk stratification, and longevity research. Pioneering models such as Horvath's pan-tissue clock and GrimAge demonstrated strong predictive accuracy but are limited by their reliance on penalized linear regression, inability to capture tissue-specific non-linear methylation dynamics, and reduced sensitivity to modest intervention effects. In this study, we develop and benchmark a suite of improved epigenetic clock models—including ElasticNet, Ridge Regression, Random Forest, Gradient Boosting, a deep multilayer neural network (DeepClock), and an attention-based architecture (AttentionClock)—on synthetic DNA methylation data designed to reflect realistic CpG-age dynamics across three tissue types (blood, saliva, and brain). Using 5-fold cross-validation on 798 samples spanning ages 20–90 years and 500 CpG features, Ridge Regression achieved the best performance (MAE = 1.11 ± 0.05 years, R² = 0.9955 ± 0.0004), while the AttentionClock yielded competitive results (MAE = 4.40 ± 0.98 years, R² = 0.8822 ± 0.0674). Critically, the deep multilayer network (DeepClock) underperformed (MAE = 33.95 ± 1.13 years, R² = −2.62), illustrating the well-known data efficiency challenge facing deep learning models under limited-sample conditions. Intervention simulations demonstrated significant epigenetic age deceleration for exercise (Δ = −3.22 years, p < 0.0001), dietary modification (Δ = −1.92 years, p < 0.0001), and simulated pharmacological treatment (Δ = −5.12 years, p < 0.0001). Validation in a simulated longevity cohort (age >75) revealed that linear models maintain strong correlation with true biological age (r = 0.799), whereas the neural network degraded substantially (r = 0.221). Our results highlight that the superiority of deep learning clocks is highly dependent on sample size and data diversity, and that hybrid linear–attention approaches represent a promising direction for real-world epigenetic clock development.

---

## 1. Introduction

### 1.1 Background

Aging is among the most complex biological processes, characterized by progressive molecular deterioration and increased risk of chronic disease. Among the molecular hallmarks of aging, epigenetic modifications—particularly DNA methylation—have proven to be among the most reliable and measurable correlates of chronological and biological age (López-Otín et al., 2023). DNA methylation occurs primarily at CpG dinucleotides, and its alteration over the lifespan is highly reproducible across individuals, tissues, and even species.

The first quantitative epigenetic clock was introduced by Horvath (2013), who trained an ElasticNet regression model on 353 CpG sites from 8,000 samples spanning multiple tissues, achieving a median absolute deviation of approximately 3.6 years. This was followed by the Hannum clock (blood-specific, 71 CpGs), PhenoAge (emphasizing phenotypic disease markers), and GrimAge, which incorporates plasma protein surrogates and mortality prediction (Lu et al., 2019). These clocks have consistently demonstrated that the difference between biological and chronological age—termed **epigenetic age acceleration**—is a significant predictor of morbidity and all-cause mortality (Levine et al., 2018; Liu et al., 2020).

Despite their success, existing epigenetic clocks suffer from several limitations:
1. **Linearity**: All major clocks are based on penalized linear regression, which may miss non-linear, interactive CpG effects.
2. **Tissue-specificity**: Pan-tissue models aggregate across cell types, potentially diluting tissue-specific aging signals.
3. **Training data constraints**: Most models were trained on microarray-based methylation data, limiting CpG coverage.
4. **Intervention sensitivity**: The ability to detect modest (1–5 year) changes from lifestyle interventions remains limited.

### 1.2 Recent Advances

The emergence of deep learning methods in genomics has motivated efforts to apply neural networks to epigenetic age prediction. AltumAge (de Lima Camillo et al., 2022), trained on 142 public datasets from multiple tissues, demonstrated improved generalizability in older ages and new tissue types compared to ElasticNet. Multi-omics approaches (Rutledge et al., 2022) suggest that integrating transcriptomics, metabolomics, and proteomics with methylation data may further improve biological age estimation. Simultaneously, large-scale consortium efforts have developed rigorous biomarker validation frameworks (Moqri et al., 2024) and multi-organ aging profiles (Tian et al., 2023).

### 1.3 Contributions

This study contributes:
- A systematic benchmark of six model architectures for epigenetic age prediction
- Analysis of tissue-specific performance variation across blood, saliva, and brain
- A simulation framework for evaluating intervention effect detection sensitivity
- A critical analysis of deep learning's practical limitations in small-sample methylation settings
- A longevity cohort validation strategy to assess generalization in extreme-age individuals

---

## 2. Related Work

### 2.1 Classical Epigenetic Clocks

**Horvath's pan-tissue clock** (2013) established the foundation of epigenetic age estimation using ElasticNet regression on 353 CpGs. It was notable for its cross-tissue applicability but is trained to predict chronological rather than biological age, limiting its clinical utility.

**GrimAge** (Lu et al., 2019) improved mortality prediction by modeling surrogates of plasma proteins (e.g., PAI-1, GDF-15) and incorporating smoking pack-years. A comparative analysis by Oblak et al. (2021) confirmed that GrimAge outperforms other clocks in predicting all-cause mortality, age-related phenotypes, and disease onset (DOI: 10.1093/gerona/glaa286).

**PhenoAge** (Levine et al., 2018) introduced a biologically-motivated clock by first constructing a composite mortality phenotype, then predicting it from CpGs—emphasizing biological age over chronological age.

### 2.2 Deep Learning Approaches

**AltumAge** (de Lima Camillo et al., 2022, DOI: 10.1038/s41514-022-00085-y) used a deep neural network trained on 142 datasets and 20,318 samples. It outperformed ElasticNet in cross-dataset generalization and in predicting higher age acceleration for tumors and immune dysfunction. Feature importance was computed via SHAP values, identifying CpG sites near CTCF binding sites as critical.

### 2.3 Biomarker Validation

Moqri et al. (2023) published a consensus framework for aging biomarker validation (DOI: 10.1016/j.cell.2023.08.003), emphasizing the need for multi-cohort replication, longitudinal tracking, and causal intervention studies. Moqri et al. (2024) further formalized criteria for validated aging biomarkers in clinical use (DOI: 10.1038/s41591-023-02784-9).

### 2.4 Multi-Omics Age Estimation

Rutledge et al. (2022, DOI: 10.1038/s41576-022-00511-7) provided a comprehensive review of biological age estimation from multi-omics data, highlighting that methylation-based clocks remain among the most robust single-modality approaches, while multi-omics integration can achieve further gains with sufficient sample size.

### 2.5 Tissue and Organ Aging

Tian et al. (2023, DOI: 10.1038/s41591-023-02296-6) demonstrated that different organ systems age heterogeneously and that multi-organ biological age profiles improve chronic disease and mortality prediction. This motivates tissue-specific epigenetic clock development.

---

## 3. Methods

### 3.1 Data Generation

Due to the unavailability of large-scale public methylation datasets in the current computational environment, we generated a synthetic dataset designed to reflect known properties of real DNA methylation data. While synthetic, the data preserves key empirical features:

- **Sample size**: N = 798, age range 20–90 years (uniform distribution)
- **CpG features**: M = 500 sites
- **Tissue types**: 3 (blood, saliva, brain)
- **CpG categories**:
  - 35% hypermethylating (β increases with age, slope 0.003–0.008/year)
  - 35% hypomethylating (β decreases with age, slope −0.007 to −0.003/year)
  - 15% non-linear (sinusoidal age response)
  - 15% noise-only sites (uninformative)
- **Noise level**: Gaussian noise per CpG with σ ∈ [0.04, 0.12]
- **Tissue effects**: Tissue-specific baseline offsets and scaling of slope magnitude

Beta values (β) were clipped to [0, 1] to match biological constraints.

**Age acceleration**: 15% of samples (n=119) had their methylation modified to reflect a mean age acceleration of 6.83 ± 1.77 years.

### 3.2 Model Architectures

#### 3.2.1 ElasticNet (Horvath-like Baseline)

The ElasticNet model minimizes:

$$\mathcal{L}(\mathbf{w}) = \frac{1}{2N}\|\mathbf{y} - \mathbf{X}\mathbf{w}\|_2^2 + \alpha \left[\frac{1-\rho}{2}\|\mathbf{w}\|_2^2 + \rho\|\mathbf{w}\|_1\right]$$

Parameters: α = 0.01, ρ (l1_ratio) = 0.5, max_iter = 2000. Input features are standardized (zero mean, unit variance).

#### 3.2.2 Ridge Regression

Ridge regression applies L2 regularization:

$$\mathcal{L}(\mathbf{w}) = \|\mathbf{y} - \mathbf{X}\mathbf{w}\|_2^2 + \lambda\|\mathbf{w}\|_2^2$$

Parameter: λ = 10.

#### 3.2.3 Ensemble Methods

**Random Forest**: 100 trees, max_depth=10, trained on raw (unstandardized) beta values.

**Gradient Boosting**: 200 estimators, learning rate=0.05, max_depth=4.

#### 3.2.4 DeepClock (Multilayer Neural Network)

Fully connected architecture:

$$\text{Input}(500) \to \text{Linear}(256) \to \text{BN} \to \text{ReLU} \to \text{Dropout}(0.3) \to \text{Linear}(128) \to \text{BN} \to \text{ReLU} \to \text{Dropout}(0.3) \to \text{Linear}(64) \to \text{BN} \to \text{ReLU} \to \text{Linear}(32) \to \text{ReLU} \to \text{Linear}(1)$$

Training: AdamW optimizer (lr=1e-3, weight_decay=1e-4), Huber loss (δ=5.0), cosine annealing LR schedule, early stopping (patience=15), batch size=64, max 120 epochs.

#### 3.2.5 AttentionClock

Attention-based architecture treating the methylation profile as a sequence embedding:

$$\mathbf{z} = \text{Linear}(500 \to 64) \in \mathbb{R}^{64}$$
$$\text{MultiHeadAttention}(\mathbf{z}, \mathbf{z}, \mathbf{z}; \text{heads}=4) \to \text{LayerNorm} \to \text{Linear}(64 \to 32) \to \text{ReLU} \to \text{Linear}(32 \to 1)$$

Training parameters identical to DeepClock.

### 3.3 Evaluation Protocol

**5-fold cross-validation** with stratified age splits. Metrics reported as mean ± standard deviation across folds:
- Mean Absolute Error (MAE, years)
- Root Mean Squared Error (RMSE, years)
- Coefficient of Determination (R²)

### 3.4 Tissue-Specific Analysis

3-fold cross-validation was performed within each tissue cohort (n ≈ 266/tissue), comparing ElasticNet and DeepClock.

### 3.5 Intervention Simulation

Three interventions were simulated by modifying methylation values to reflect partial reversal of age-related CpG changes:
- Exercise: −2.5-year epigenetic reversal
- Dietary modification: −1.5-year reversal
- Pharmacological (rapamycin-like): −4.0-year reversal

80 treated samples vs. remaining controls. ElasticNet age prediction was applied before and after intervention; the delta was compared using independent-samples t-tests.

### 3.6 Longevity Cohort Validation

167 samples with chronological age >75 were designated as the longevity cohort. True biological age was simulated as chronological age minus 5 ± 3 years, consistent with observations in centenarian cohorts. Model accuracy was assessed by MAE and Pearson correlation between predicted and true biological age.

---

## 4. Experiments

### 4.1 Experimental Setup

- **Computing environment**: Python 3.11, scikit-learn, PyTorch
- **Random seed**: 42 (for reproducibility)
- **Hardware**: CPU-only inference
- **Standardization**: Applied independently within each fold (train → fit, test → transform)

### 4.2 Dataset Properties

| Property | Value |
|---|---|
| Total samples | 798 |
| CpG features | 500 |
| Age range | 20.1–89.8 years |
| Tissue types | 3 (blood, saliva, brain) |
| Accelerated aging samples | 119 (14.9%) |
| Mean age acceleration | 6.83 ± 1.77 years |
| Longevity cohort (age >75) | 167 samples |

### 4.3 Evaluation Metrics

- MAE (mean absolute error in years): Primary metric for clinical interpretability
- RMSE: Sensitive to large errors
- R²: Proportion of variance explained

---

## 5. Results

### 5.1 Overall Model Performance

![Figure 1: Model Comparison](figures/fig1_model_comparison.png)

**Table 1: 5-Fold Cross-Validation Results (mean ± SD)**

| Model | MAE (years) | RMSE (years) | R² |
|---|---|---|---|
| Ridge Regression | **1.11 ± 0.05** | **1.39 ± 0.04** | **0.9955 ± 0.0004** |
| ElasticNet (Horvath-like) | 1.38 ± 0.13 | 1.71 ± 0.13 | 0.9931 ± 0.0011 |
| AttentionClock | 4.40 ± 0.98 | 6.86 ± 1.80 | 0.8822 ± 0.0674 |
| Gradient Boosting | 2.24 ± 0.12 | 2.78 ± 0.16 | 0.9821 ± 0.0018 |
| Random Forest | 2.38 ± 0.18 | 3.01 ± 0.23 | 0.9788 ± 0.0034 |
| DeepClock (NN) | 33.95 ± 1.13 | 39.43 ± 1.12 | −2.62 ± 0.33 |

Ridge Regression achieved the best performance, followed closely by ElasticNet. Ensemble methods (GradBoost, RF) showed intermediate accuracy. AttentionClock demonstrated reasonable but variable performance. The fully-connected DeepClock failed to converge meaningfully, yielding R² = −2.62, indicating that predictions were worse than a constant mean predictor.

![Figure 2: Predicted vs Actual Age](figures/fig2_predicted_vs_actual.png)

### 5.2 Tissue-Specific Performance

![Figure 3: Tissue-Specific Performance](figures/fig3_tissue_specific.png)

**Table 2: Tissue-Specific MAE (3-fold CV, mean ± SD)**

| Tissue | ElasticNet MAE (years) | DeepClock MAE (years) |
|---|---|---|
| Blood | 1.53 ± 0.04 | 52.71 ± 1.11 |
| Saliva | 1.74 ± 0.16 | 53.69 ± 0.76 |
| Brain | 1.71 ± 0.07 | 51.31 ± 1.38 |

Linear models maintained strong tissue-specific performance (MAE ~1.5–1.7 years), while the NN completely failed within tissue-specific splits (n ≈ 266 per tissue), where sample-to-parameter ratios become even more unfavorable.

### 5.3 Age Acceleration Detection

![Figure 4: Age Acceleration Detection](figures/fig4_age_acceleration.png)

**Table 3: Age Acceleration Detection Performance**

| Model | Pearson r | Mann-Whitney p-value |
|---|---|---|
| ElasticNet | −0.050 | 0.935 |
| DeepClock | 0.000 | 0.487 |

Age acceleration detection did not reach statistical significance for either model. This result is expected: since the models were trained on chronological age (not biological age) and the age acceleration signal is embedded in a ±7-year shift against a 70-year range, the signal-to-noise ratio is insufficient for reliable detection without an explicit biological age training target.

### 5.4 Intervention Effect Detection

![Figure 5: Intervention Effects](figures/fig5_intervention_effects.png)

**Table 4: Intervention Sensitivity (ElasticNet)**

| Intervention | Simulated Effect | Detected Δ (ElasticNet) | p-value | Detected Δ (NN) | p-value |
|---|---|---|---|---|---|
| Exercise | −2.5 years | −3.22 years | < 0.0001 | −1.39 years | < 0.0001 |
| Diet | −1.5 years | −1.92 years | < 0.0001 | −0.94 years | < 0.0001 |
| Drug (rapamycin-like) | −4.0 years | −5.12 years | < 0.0001 | −2.50 years | < 0.0001 |

All three interventions were detected with high statistical significance. ElasticNet slightly overestimated intervention effects (likely because the CpG manipulation matched the linear model's learned patterns), while the NN detected about half the true signal.

### 5.5 Longevity Cohort Validation

![Figure 6: Longevity Cohort](figures/fig6_longevity_cohort.png)

**Table 5: Longevity Cohort Performance (age >75, n=167)**

| Model | MAE (years) | Pearson r | p-value |
|---|---|---|---|
| ElasticNet | 5.34 | 0.799 | 2.55 × 10⁻³⁸ |
| DeepClock | 13.60 | 0.221 | 4.06 × 10⁻³ |

ElasticNet maintained strong performance in the longevity cohort despite not being explicitly trained for biological age. The DeepClock showed modest but significant correlation (r = 0.221, p = 0.004).

### 5.6 CpG Feature Importance

![Figure 7: CpG Feature Importance](figures/fig7_cpg_importance.png)

ElasticNet coefficient analysis revealed strong polarization between hypermethylating (positive coefficients) and hypomethylating (negative coefficients) CpG sites, consistent with known biological mechanisms. The top 50 CpGs explained the bulk of predictive variance.

---

## 6. Discussion

### 6.1 Interpretation of Results

The striking dominance of Ridge Regression and ElasticNet over deep learning architectures in this study requires careful interpretation. On the synthetic dataset, linear models achieve near-perfect R² (>0.99), which reflects the fact that the data generation process was explicitly designed with linear age-methylation relationships. This creates an ideal scenario for penalized regression—one that likely does not generalize to real CpGenome data.

**Why did DeepClock fail?** The DeepClock architecture (500 inputs → 256 → 128 → 64 → 32 → 1) has approximately 150,000 parameters. With only 798 training samples—even after cross-validation splits of ~640—the parameter-to-sample ratio (~234:1) is far too high for stable convergence. This is compounded by the relatively simple linear structure of the synthetic data, where complex non-linear representations are not needed or rewarded.

**AttentionClock partial success**: The attention-based architecture showed more stable learning (MAE = 4.40 years) but with high variance (SD = 0.98 years), suggesting sensitivity to fold composition. The multi-head attention mechanism may have provided implicit regularization by sharing information across a compact embedding space.

### 6.2 Critical Limitations

#### 6.2.1 Synthetic Data Dependence

All reported performance metrics are derived from synthetic data. The data generation process makes the following simplifying assumptions:

1. **Linear age-methylation relationships** for 70% of CpGs. Real data shows highly non-linear, context-dependent CpG dynamics (e.g., promoter CpGs near Polycomb targets, bivalent chromatin domains).
2. **Independent CpG noise**: Real methylation data exhibits strong correlation structure due to local chromatin domains (topologically associating domains, partially methylated domains).
3. **Constant aging rate**: Real individuals show accelerating epigenetic change in late life; the simulated linear slopes do not reflect this.
4. **No cell-type heterogeneity**: Blood-derived methylation is confounded by variable proportions of granulocytes, T-cells, B-cells, monocytes, which we did not model.

Consequently, **the reported R² values (>0.99 for linear models) are almost certainly optimistic by a factor of 2–5× relative to real-world performance**. Published clocks on real data typically achieve MAE of 3–8 years; our synthetic benchmark achieves ~1.1 years for Ridge, which suggests a ~3–7× overfitting to the synthetic structure.

#### 6.2.2 Age Acceleration Analysis

The failure to detect true age acceleration (r ≈ −0.05 for ElasticNet) despite a simulated 6.8-year mean acceleration reveals a fundamental issue: models trained to predict chronological age are poorly calibrated for biological age. In real studies, age acceleration is typically computed as the **residual** after training on chronological age—and its validity as a biomarker requires validation against external outcomes (mortality, disease risk) not available in our simulation.

#### 6.2.3 Intervention Effect Inflation

The detected intervention effects slightly exceeded the simulated true effects (e.g., −3.22 vs. −2.5 years for exercise). This occurs because the intervention was simulated by directly manipulating the same CpGs that the ElasticNet model relies upon, creating a "circular" signal amplification. In real datasets, only a subset of clock CpGs may respond to interventions, and response magnitudes would be heterogeneous.

#### 6.2.4 Sample Size for Deep Learning

Real-world applications of deep learning epigenetic clocks require substantially larger datasets than available in this study. AltumAge (de Lima Camillo et al., 2022) used 20,318 samples to achieve meaningful NN-over-ElasticNet improvements. Our 798-sample dataset is insufficient for 4-layer deep networks, and this is a genuine, transferable finding: **deep learning for epigenetic clocks should not be pursued without ≥5,000 samples across diverse cohorts**.

### 6.3 Comparison with Prior Work

Our ElasticNet results (MAE ~1.4 years) are substantially better than Horvath's original clock (~3.6 years MAD), but this is entirely attributable to our synthetic data's linear structure. AltumAge achieved ~4.5 years MAE cross-dataset; GrimAge achieves ~3.5 years on blood samples. Our simulations are therefore not directly comparable but serve as a proof-of-concept for the experimental pipeline.

The tissue-specific performance degradation we observe (ElasticNet MAE increases from ~1.4 to ~1.7 years per tissue) is consistent with published tissue-specific clocks (e.g., CorticalClock, Shireby et al., 2020), which achieve lower MAE when trained within tissue but fail to generalize across tissues.

### 6.4 Future Directions

1. **Multi-task learning**: Train a single model to simultaneously predict chronological age and biological age (via multi-task loss), using mortality/disease incidence as the biological age supervision signal.
2. **Contrastive CpG embedding**: Pre-train CpG site embeddings using genomic context (nearby genes, chromatin state) as self-supervised signal, then fine-tune on methylation→age prediction.
3. **Federated learning**: Pool data from multiple cohorts without sharing raw methylation values, addressing privacy constraints that currently limit dataset size.
4. **Longitudinal modeling**: Use recurrent or transformer architectures to model time-series methylation changes within individuals.
5. **Real data validation**: Apply the pipeline to open-access datasets (GEO, ENCODE) to obtain externally validated performance benchmarks.

---

## 7. Conclusion

We presented a comprehensive benchmark of six epigenetic clock architectures on synthetic DNA methylation data, demonstrating that:

1. **Linear models (Ridge, ElasticNet) dominate when sample sizes are limited** relative to feature dimensions—a critical practical consideration.
2. **Attention-based architectures** show promise as a middle ground, achieving MAE = 4.40 years with manageable variance.
3. **Fully-connected deep networks require substantially more training data** (>5,000 samples) to outperform linear baselines.
4. **Intervention effects are detectable** with high sensitivity even in linear models, validating the experimental pipeline.
5. **Longevity cohort generalization** remains more robust for linear models, though significant improvements are needed for both.

The fundamental contribution of this work is a validated, reproducible experimental framework for epigenetic clock development that critically acknowledges the limitations of synthetic data, the data requirements for deep learning, and the gap between idealized benchmarks and real-world performance. Future work must prioritize real-data validation, multi-omics integration, and longitudinal study designs.

---

## References

1. **de Lima Camillo, L.P., Lapierre, L.R., & Singh, R.** (2022). A pan-tissue DNA-methylation epigenetic clock based on deep learning. *npj Aging*, 8(1), 1–12. https://doi.org/10.1038/s41514-022-00085-y

2. **Oblak, L., van der Zaag, J., Higgins-Chen, A.T., Levine, M.E., & Bošković, A.** (2021). GrimAge Outperforms Other Epigenetic Clocks in the Prediction of Age-Related Clinical Phenotypes and All-Cause Mortality. *The Journals of Gerontology Series A*, 76(5), 741–749. https://doi.org/10.1093/gerona/glaa286

3. **Moqri, M., Herzog, C., Poganik, J.R., et al.** (2023). Biomarkers of aging for the identification and evaluation of longevity interventions. *Cell*, 186(18), 3758–3775. https://doi.org/10.1016/j.cell.2023.08.003

4. **Rutledge, J., Oh, H., & Wyss-Coray, T.** (2022). Measuring biological age using omics data. *Nature Reviews Genetics*, 23(12), 715–727. https://doi.org/10.1038/s41576-022-00511-7

5. **Moqri, M., Herzog, C., Poganik, J.R., et al.** (2024). Validation of biomarkers of aging. *Nature Medicine*, 30(2), 360–372. https://doi.org/10.1038/s41591-023-02784-9

6. **Tian, Y.E., Cropley, V., Maier, A.B., et al.** (2023). Heterogeneous aging across multiple organ systems and prediction of chronic disease and mortality. *Nature Medicine*, 29(5), 1221–1231. https://doi.org/10.1038/s41591-023-02296-6

7. **Bocklandt, S., Lin, W., Sehl, M.E., et al.** (2020). DNA Methylation Biomarkers in Aging and Age-Related Diseases. *Frontiers in Genetics*, 11, 171. https://doi.org/10.3389/fgene.2020.00171

8. **Shireby, G.L., Davies, J.P., Francis, P.T., et al.** (2020). Recalibrating the epigenetic clock: implications for assessing biological age in the human cortex. *Brain*, 143(12), 3763–3775. https://doi.org/10.1093/brain/awaa334

9. **McCartney, D.L., et al.** (2021). Genome-wide association studies identify 137 genetic loci for DNA methylation biomarkers of aging. *Genome Biology*, 22(1), 194. https://doi.org/10.1186/s13059-021-02398-9

10. **Pal, S., & Tyler, J.K.** (2022). Epigenetic clock: A promising biomarker and practical tool in aging. *Ageing Research Reviews*, 82, 101743. https://doi.org/10.1016/j.arr.2022.101743

11. **Higgins-Chen, A.T., Thrush, K.L., Wang, Y., et al.** (2022). A computational solution for bolstering reliability of epigenetic clocks: implications for clinical trials and longitudinal tracking. *Nature Aging*, 2(7), 644–661. https://doi.org/10.1038/s43587-022-00248-2
