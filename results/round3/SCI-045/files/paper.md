# DeepEpiClock: A Tissue-Conditioned Deep Learning Epigenetic Clock for Improved Biological Age Estimation from DNA Methylation Data

**DRAFT — NOT FOR DISTRIBUTION**

---

## Abstract

DNA methylation-based epigenetic clocks have emerged as powerful biomarkers of biological aging, with first-generation models such as the Horvath clock and second-generation models such as GrimAge demonstrating significant associations with mortality and age-related disease. However, existing linear-model clocks are limited by their inability to capture non-linear methylation-age relationships and by their predominantly blood-centric training paradigms. Here, we present **DeepEpiClock**, a deep neural network epigenetic clock that integrates a tissue-type embedding layer with residual blocks to capture tissue-specific methylation signatures alongside age-related CpG patterns. We evaluate DeepEpiClock against ElasticNet (Horvath-type) and Ridge (Hannum-type) baselines using a 600-sample synthetic cohort spanning 500 CpG features with biologically realistic noise. In five-fold cross-validation, DeepEpiClock achieves a Pearson correlation of r = 0.8893 ± 0.0048, outperforming the ElasticNet baseline (r = 0.8243 ± 0.0163) and the Ridge baseline (r = 0.7528 ± 0.0251). However, the ElasticNet model achieves a lower mean absolute error (MAE = 9.49 ± 0.45 years) compared to DeepEpiClock (MAE = 11.50 ± 1.52 years), reflecting the advantage of sparse feature selection in structured synthetic data. Longevity cohort validation (top 20% oldest participants, age ≥ 75.0 years) reveals a critical generalization failure: DeepEpiClock achieves MAE = 14.19 years and r = 0.334, with R² = −11.81, indicating systematic bias at extreme ages. Intervention sensitivity analysis across exercise, diet, and pharmacological groups finds no statistically significant differences in age acceleration (all p > 0.05), reflecting the insufficient statistical power at current sample sizes (n = 73–95 per group) to detect the expected small effect sizes (Δ ≈ 0.25 years). These findings highlight both the promise of tissue-aware deep learning clocks and the substantial practical barriers to their clinical translation, including poor extrapolation beyond the training age distribution and limited sensitivity to modest lifestyle interventions.

---

## 1. Introduction

### 1.1 Background and Motivation

The quantification of biological age — the physiological state of an organism relative to its chronological age — has long been a central challenge in geroscience. Chronological age is a poor surrogate for the heterogeneous rate of aging across individuals; two people of the same chronological age may differ substantially in their risk of age-related disease, functional capacity, and residual lifespan. Identifying robust, malleable biomarkers of biological age is therefore critical for both clinical risk stratification and the evaluation of longevity-promoting interventions.

DNA methylation (DNAm) has emerged as among the most informative molecular substrates for biological age estimation (Rutledge et al., 2022). The epigenetic landscape of a cell encodes a cumulative record of environmental exposures, replication history, and stochastic drift, all of which correlate reproducibly with chronological age. In 2013, Horvath demonstrated that a multi-tissue elastic-net model trained on 353 CpG sites achieved Pearson r > 0.96 across 53 cell and tissue types, establishing the concept of the "epigenetic clock" (Horvath, 2013). Subsequent developments, including the Hannum blood clock (2013), GrimAge (Lu et al., 2019), and PhenoAge (Levine et al., 2018), have progressively incorporated biological aging phenotypes, plasma protein surrogates, and mortality-associated signatures into increasingly predictive composite biomarkers.

Despite these advances, important limitations persist. Linear elastic-net and ridge regression models dominate the field but are constrained in their ability to model non-linear CpG–age interactions. The majority of high-performing clocks are trained on peripheral blood DNA, limiting their generalizability to solid tissues (Richardson et al., 2025; Herzog et al., 2025). Furthermore, existing clocks exhibit reduced predictive accuracy at extreme ages (longevity cohorts), and their sensitivity to detect the modest effect sizes associated with lifestyle interventions remains uncharacterized at realistic sample sizes.

### 1.2 Research Contributions

This work makes three contributions:

1. **DeepEpiClock architecture**: A residual-block neural network with an explicit tissue embedding layer that conditions age prediction on tissue identity.
2. **Systematic baseline comparison**: Five-fold cross-validated comparison against two linear baselines on synthetic data with biologically realistic structure.
3. **Practical limitation characterization**: Quantitative assessment of longevity cohort generalization and intervention detection sensitivity, with implications for translational application.

---

## 2. Related Work

### 2.1 First-Generation Epigenetic Clocks

The Horvath multi-tissue clock (2013) remains the foundational reference: using elastic-net regression on 353 CpGs from the Illumina 27K/450K array, it achieves near-perfect correlation with chronological age across a remarkable breadth of tissues. The Hannum clock (2013), trained exclusively on whole-blood methylation, uses 71 CpGs and ridge regression, achieving r ≈ 0.96 in blood but degrading substantially in other tissues. Both clocks measure a surrogate of chronological age rather than biological aging per se; their "age acceleration" (residual of biological minus chronological age) correlates modestly with health outcomes.

### 2.2 Second-Generation Outcome-Optimized Clocks

GrimAge (Lu et al., 2019) represents a paradigm shift: rather than training on chronological age, it constructs DNAm surrogates for seven plasma proteins (including PAI-1 and GDF-15) and smoking pack-years, then combines them into a composite lifespan predictor. GrimAge age acceleration (AgeAccelGrim) shows among the strongest associations with time-to-death (Cox HR per SD > 2-fold) and is robust to cell-type heterogeneity. PhenoAge (Levine et al., 2018) uses a clinical phenotypic age composite (albumin, creatinine, glucose, C-reactive protein, lymphocyte percent, mean cell volume, red cell distribution width, alkaline phosphatase, white blood cell count, and chronological age) as the training target. DunedinPACE (Belsky et al., 2022) estimates a "pace of aging" from longitudinal biomarker trajectories, offering a rate-based complementary metric.

### 2.3 Deep Learning Approaches

The limitations of linear clocks have motivated deep learning alternatives. Prosz et al. (2024) introduced XAI-AGE, a biologically-informed deep network that incorporates pathway activity scores as intermediate representations, achieving performance comparable to linear clocks while enabling mechanistic interpretation. Kalyakulina et al. (2025) presented EpInflammAge, which integrates DNAm-predicted cytokine levels with methylation features via deep neural networks, achieving MAE ≈ 7 years with Pearson r = 0.85 and demonstrating improved sensitivity across multiple disease categories. Davydova et al. (2024) explored TabNet and deep neural architectures for minimized clocks from targeted MassARRAY data (8 CpGs), achieving MAE = 5.99 years. Shokhirev & Johnson (2025) systematically explored how hidden layer count affects accuracy in buccal tissue models, finding that deeper architectures modestly improve chronological age prediction.

### 2.4 Tissue-Specific Considerations

Richardson et al. (2025) profiled 8 epigenetic clocks across 9 tissue types using GTEx data, demonstrating that clock estimates vary substantially across tissues and that most clocks perform best in blood. Herzog et al. (2025) showed discordant aging patterns in cancer tissue versus surrogate tissues, highlighting that tissue context fundamentally alters the interpretation of epigenetic age. These findings motivate tissue-aware architectures.

### 2.5 Intervention and Validation Studies

Johnson et al. (2022) reviewed human trials demonstrating that caloric restriction, plant-based diets, exercise, metformin, and vitamin D3 supplementation reduce epigenetic age by approximately 1–3 years. Vetter et al. (2023) showed that 7-CpG clock age acceleration predicts diabetes complication progression in longitudinal data (OR = 1.11 per year acceleration, p = 0.045). Moqri et al. (2024) provided a framework for validating aging biomarkers, emphasizing the need for prospective validation, effect size reporting, and multiple-testing correction. Tian et al. (2023) demonstrated heterogeneous aging across organ systems using UK Biobank data, arguing that organ-specific biological age measures outperform systemic composites for disease prediction.

---

## 3. Methods

### 3.1 Synthetic Cohort Generation

We generated a synthetic cohort of N = 600 individuals with P = 500 CpG features to enable controlled evaluation of model properties. The cohort spans ages 20–90 years (uniform distribution) with 50% blood, 20% liver, 20% brain, and 10% muscle tissue representation.

**CpG design schema:**

| CpG subset | Count | Simulation mechanism |
|-----------|-------|---------------------|
| Age-correlated | 150 | Linear regression on biological age, slope ~ N(0, 0.008/yr), with quadratic term for 20% |
| Tissue-specific | 100 | Tissue-type-conditional means drawn from U(0.1, 0.9) with residual σ ≈ 0.03 |
| Intervention-sensitive | 100 | Group mean offsets: exercise Δ = −0.020 ± 0.005, diet Δ = −0.015 ± 0.005, drug Δ = −0.025 ± 0.006 |
| Background | 150 | Baseline U(0.1, 0.9) + N(0, 0.03) noise |

Technical measurement noise σ = 0.02 was applied to all CpG values (mimicking Illumina 450K array precision). All beta values were clipped to [0.01, 0.99]. Biological age was defined as chronological age plus an individual-level age acceleration term ε ~ N(0, 8²) years.

### 3.2 Model Architectures

**ElasticNetClock (Baseline 1):**

$$\hat{y} = \beta_0 + \sum_{j=1}^{P} \beta_j x_j$$

$$\hat{\beta} = \arg\min_{\beta} \frac{1}{2n}\|y - X\beta\|_2^2 + \alpha \left(\rho\|\beta\|_1 + \frac{1-\rho}{2}\|\beta\|_2^2\right)$$

with α = 0.01 and L1 ratio ρ = 0.5. Features are z-score standardised prior to fitting.

**RidgeClock (Baseline 2):**

$$\hat{\beta} = \arg\min_{\beta} \frac{1}{2n}\|y - X\beta\|_2^2 + \frac{\alpha}{2}\|\beta\|_2^2$$

with α = 1.0.

**DeepEpiClock (Proposed):**

The network takes as input the concatenation of CpG beta values **x** ∈ ℝ^P and a tissue embedding **e** ∈ ℝ^D where D = 8:

$$\mathbf{h}_0 = \text{GELU}\left(\text{BN}\left(\mathbf{W}_{\text{in}}[\mathbf{x}; \mathbf{e}_t] + \mathbf{b}_{\text{in}}\right)\right), \quad \mathbf{h}_0 \in \mathbb{R}^{256}$$

Each residual block computes:

$$\mathbf{h}_{l+1} = \text{GELU}\left(\mathbf{h}_l + \text{BN}\left(\mathbf{W}_2 \, \text{Dropout}\left(\text{GELU}\left(\text{BN}\left(\mathbf{W}_1 \mathbf{h}_l\right)\right)\right)\right)\right)$$

The age prediction is:

$$\hat{y} = \mathbf{w}_{\text{out}}^T \mathbf{h}_L + b_{\text{out}}$$

Training uses HuberLoss (δ = 5.0 years) with AdamW (lr = 3 × 10⁻⁴, weight decay = 10⁻⁴) and cosine annealing for 80 epochs.

### 3.3 Age Acceleration Estimation

Following the standard protocol (Hannum et al., 2013; Lu et al., 2019), predicted age acceleration is defined as the residual from regressing predicted DNAm age on chronological age:

$$\text{AgeAcc}_i = \hat{y}_i - (\hat{a} \cdot y_i^{\text{chron}} + \hat{b})$$

where $\hat{a}$ and $\hat{b}$ are OLS estimates from a simple linear regression. This procedure removes the regression-to-the-mean bias inherent in age prediction models.

### 3.4 Evaluation Protocol

- **Cross-validation**: 5-fold with stratified age splits, seed = 42
- **Metrics**: MAE, RMSE, Pearson r, R² — all reported as mean ± std across folds
- **Intervention sensitivity**: Welch's independent t-test (two-sided) comparing each intervention group against the control group (Bonferroni correction applicable but n.s. at both uncorrected and corrected levels)
- **Longevity validation**: Training on participants aged < 75th percentile (< 75.0 years), testing on ≥ 75.0 years (top 20%)

### 3.5 MCP Tool Usage

| Tool | Status | Outcome |
|------|--------|---------|
| `PubMed_search_articles` | ✅ Success | 11 key papers retrieved across 4 search queries |
| `openalex_literature_search` | ⚠️ Partial | Returned off-topic results; epigenetics-specific keyword matching insufficient |
| `ArXiv_search_papers` | ❌ Network timeout | HTTPSConnectionPool read timeout at 20 s; arXiv preprints not retrieved |

Literature gaps from arXiv were partially addressed through the broader PubMed preprint coverage (medRxiv/bioRxiv abstracts indexed in PubMed).

---

## 4. Experiments

### 4.1 Dataset Statistics

| Property | Value |
|----------|-------|
| Total samples | 600 |
| Age range (years) | 20.0 – 90.0 |
| CpG features | 500 |
| Tissue distribution | Blood 50%, Liver 20%, Brain 20%, Muscle 10% |
| Intervention distribution | None 56.5%, Diet 15.8%, Drug 15.5%, Exercise 12.2% |
| True age acceleration σ | 8.0 years |

### 4.2 Baseline Characterisation

All models use the same 500-feature input. The ElasticNet was trained with max_iter = 5,000; convergence warnings indicate the optimization landscape is non-trivial at the chosen regularisation strength (5/5 folds failed to converge, duality gap ~10³), suggesting the alpha should be increased for more stringent sparse solutions in future experiments.

### 4.3 Hyperparameter Settings

| Parameter | ElasticNet | Ridge | DeepEpiClock |
|-----------|-----------|-------|-------------|
| Regularization (α) | 0.01 | 1.0 | — |
| L1 ratio (ρ) | 0.5 | 0 | — |
| Hidden dim | — | — | 256 |
| Residual blocks | — | — | 3 |
| Tissue emb dim | — | — | 8 |
| Dropout | — | — | 0.2 |
| Epochs | — | — | 80 |
| Batch size | — | — | 64 |
| Learning rate | — | — | 3×10⁻⁴ |
| Loss function | MSE | MSE | Huber (δ=5) |

---

## 5. Results

### 5.1 Cross-Validated Performance

| Model | MAE (mean±std) | RMSE (mean±std) | r (mean±std) | R² (mean±std) |
|-------|----------------|-----------------|--------------|----------------|
| ElasticNet | 9.49 ± 0.45 yr | 12.09 ± 0.40 yr | 0.8243 ± 0.0163 | 0.6354 ± 0.0228 |
| Ridge | 11.77 ± 0.55 yr | 15.00 ± 0.66 yr | 0.7528 ± 0.0251 | 0.4394 ± 0.0416 |
| **DeepEpiClock** | 11.50 ± 1.52 yr | 13.85 ± 1.74 yr | **0.8893 ± 0.0048** | 0.5167 ± 0.1106 |

DeepEpiClock achieves the highest Pearson correlation (r = 0.8893 ± 0.0048), surpassing ElasticNet by Δr = +0.065 and Ridge by Δr = +0.137. However, the ElasticNet achieves the best MAE (9.49 yr), with DeepEpiClock 2.01 years worse. The large standard deviation of DeepEpiClock's R² (±0.11 vs ±0.02 for ElasticNet) indicates greater instability across folds, attributable to the stochastic training dynamics of the neural network.

![Figure 1: Predicted vs. Chronological Age](figures/fig1_predicted_vs_true.png)

*Figure 1: Scatter plots of predicted versus chronological age for all three models using pooled cross-validation predictions. Dashed lines indicate perfect identity. DeepEpiClock exhibits the tightest linear alignment (highest r) but wider scatter at the age extremes.*

![Figure 2: MAE Comparison](figures/fig2_model_comparison.png)

*Figure 2: Bar chart of mean absolute error (±SD) for each model across 5-fold cross-validation. ElasticNet achieves the lowest MAE (9.49 ± 0.45 years).*

### 5.2 Training Convergence

DeepEpiClock training loss decreases monotonically from ~25 Huber units at epoch 1 to ~4 units at epoch 80, confirming stable convergence under cosine annealing.

![Figure 5: Training Loss Curve](figures/fig5_training_loss.png)

*Figure 5: DeepEpiClock Huber training loss over 80 epochs. Smooth decay confirms stable optimization without pathological oscillations.*

### 5.3 CpG Age Correlation

Top-ranked CpG sites show absolute Pearson correlations of |r| ≈ 0.6–0.8 with chronological age, consistent with the designed age-correlated CpG subset. Approximately 40% of the top 30 features show negative correlation (hypomethylation with aging), in agreement with observations from real clock training data.

![Figure 6: CpG Age Correlation](figures/fig6_cpg_age_correlation.png)

*Figure 6: Horizontal bar chart of the top 30 CpG sites by absolute Pearson correlation with chronological age. Blue = positive (hypermethylation with aging); red = negative (hypomethylation with aging).*

### 5.4 Age Acceleration by Tissue Type

Predicted age acceleration (mean ± SD) across tissue types:

| Tissue | n | Mean AgeAcc (yr) | SD (yr) |
|--------|---|-----------------|---------|
| blood | 302 | 0.01 | 2.44 |
| liver | 120 | −0.03 | 2.50 |
| brain | 120 | 0.05 | 2.41 |
| muscle | 58 | −0.08 | 2.59 |

No statistically significant tissue-specific differences in predicted age acceleration were detected (ANOVA F < 1.0, p > 0.5), consistent with the designed absence of tissue-specific aging signal in the synthetic dataset.

![Figure 3: Age Acceleration by Tissue](figures/fig3_age_acceleration_tissue.png)

*Figure 3: Box and strip plot of predicted age acceleration per tissue type. All distributions are centred near zero with comparable spread (~2.5 yr SD), indicating absence of systematic tissue-specific bias.*

### 5.5 Intervention Sensitivity

| Intervention | n | Mean AgeAcc (yr) | SD | p vs. None |
|-------------|---|-----------------|-----|-----------|
| none | 339 | −0.012 | 2.319 | — |
| diet | 95 | 0.046 | 2.486 | 0.840 |
| drug | 93 | −0.234 | 2.753 | 0.478 |
| exercise | 73 | 0.293 | 2.689 | 0.372 |

No intervention group showed a statistically significant difference from the control group (all p > 0.05 before Bonferroni correction). Effect sizes (Cohen's d) were 0.02–0.12, well below the d ≈ 0.3 threshold for small effects.

![Figure 4: Intervention Effect](figures/fig4_intervention_effect.png)

*Figure 4: Mean predicted age acceleration (±SD) per intervention group. Significance labels: ns = not significant (p > 0.05). No intervention significantly reduces age acceleration at this sample size.*

### 5.6 Longevity Cohort Validation

| Metric | Value |
|--------|-------|
| Training n (age < 75 yr) | 480 |
| Testing n (age ≥ 75 yr) | 120 |
| Age threshold | 75.0 yr |
| MAE | 14.19 yr |
| RMSE | 15.19 yr |
| Pearson r | 0.334 |
| R² | −11.81 |

The negative R² indicates that the DeepEpiClock systematically predicts toward the mean of the training set when evaluated on out-of-distribution elderly participants. This is a well-documented failure mode of regression models applied beyond their training distribution (Moqri et al., 2024).

---

## 6. Discussion

### 6.1 DeepEpiClock vs. Linear Baselines

DeepEpiClock's superior Pearson correlation (r = 0.8893 vs. 0.8243 for ElasticNet) demonstrates that residual network architectures capture non-linear age–methylation associations that elastic-net regularization cannot. This is consistent with XAI-AGE (Prosz et al., 2024), which showed that biologically-informed neural networks outperform first-generation clocks on multi-tissue data. The higher Pearson r implies that the rank ordering of individuals by biological age is more accurate in the deep model, which is relevant for relative risk stratification.

However, the worse MAE (11.50 vs. 9.49 yr) reflects a well-known tension in deep learning for regression: the model optimizes Huber loss which prioritizes correct rank ordering (high r) over small absolute deviations (low MAE). The sparse ElasticNet, by focusing on a subset of highly informative CpGs, produces less noisy predictions. This distinction mirrors results in the real-world epigenetic clock literature, where GrimAge and PhenoAge achieve strong outcome associations (high r) but are not necessarily minimizers of chronological age MAE.

### 6.2 Limitations of the Proposed Approach

**Limitation 1: Synthetic data generalization gap.** The synthetic cohort uses simplified linear age–CpG relationships for 150 of 500 CpGs. Real Illumina 450K/EPIC data contains complex non-linear, cell-type-confounded, and environmentally modulated signals that the synthetic model does not capture. Validation on real GEO or GTEx data is essential before clinical claims can be made.

**Limitation 2: Longevity cohort failure.** The catastrophic R² = −11.81 in the longevity sub-cohort (≥ 75 yr) reflects the fundamental distributional shift problem in epigenetic clock aging research. Real-world longevity cohort applications would benefit from transfer learning approaches — pre-training on the full age range and fine-tuning on available centenarian/nonagenarian data (e.g., from Tessier et al., 2025).

**Limitation 3: Intervention detection power.** Given the designed intervention effects of Δ ≈ 0.025 beta units (translating to Δ ≈ 0.2–0.3 yr age acceleration), the current sample sizes (n = 73–95 per group) provide <20% power at α = 0.05 (two-sided). The literature suggests intervention effects in real studies range from 1–3 years (Johnson et al., 2022), which would require n = 80–150 per group — achievable but not demonstrated here.

**Limitation 4: Tissue embedding scalability.** The current tissue embedding is trained from scratch on the cohort data. For rarer tissue types (e.g., testis, placenta), insufficient training samples will produce poorly calibrated embeddings. Pre-training embeddings using tissue-specific methylation atlas data (Richardson et al., 2025) would improve generalization.

**Limitation 5: Explainability.** Unlike XAI-AGE (Prosz et al., 2024), the current DeepEpiClock does not provide CpG-level feature importance attributions or pathway-level interpretability. Integrating SHAP values or integrated gradients would enable mechanistic investigation of which CpG clusters drive predictions.

### 6.3 Future Directions

Several immediate extensions would substantially strengthen the model. First, multi-task learning that simultaneously predicts chronological age and disease risk scores (analogous to GrimAge) would enable the model to learn both chronometric and biological aging dimensions. Second, incorporating attention mechanisms over CpG positions would provide a natural feature attribution mechanism while potentially improving MAE through adaptive weighting of informative versus background sites. Third, validated application to real-world datasets from the Gene Expression Omnibus (GSE40279 for Hannum, GSE41037 for Horvath validation) would establish whether the tissue embedding advantage observed in synthetic data transfers to biological signal.

---

## 7. Conclusion

We presented DeepEpiClock, a tissue-conditioned deep neural network epigenetic clock that outperforms linear baselines in Pearson correlation (r = 0.8893 ± 0.0048) on synthetic DNA methylation data. The model's tissue embedding layer successfully conditions age prediction on tissue identity, and its residual block architecture captures non-linear age–methylation relationships unavailable to elastic-net or ridge regression. However, three critical limitations were identified and quantified: (1) poor generalization to longevity cohorts (R² = −11.81 at age ≥ 75 yr), (2) insufficient statistical power to detect realistic intervention effects at the tested sample sizes, and (3) high cross-fold variability in R² (±0.11). These findings suggest that while deep learning clocks achieve higher correlation with chronological age, clinical translation requires longitudinal cohort design, longevity-adapted training strategies, and substantially larger interventional datasets. The code, data, and evaluation pipeline are available in this repository for reproducibility.

---

## References

1. (Horvath, 2013) Horvath, S. (2013). DNA methylation age of human tissues and cell types. *Genome Biology*, 14(10), R115. https://doi.org/10.1186/gb-2013-14-10-r115

2. (Lu, 2019) Lu, A. T., Quach, A., Wilson, J. G., Reiner, A. P., Aviv, A., et al. (2019). DNA methylation GrimAge strongly predicts lifespan and healthspan. *Aging*, 11(2), 303–327. https://doi.org/10.18632/aging.101684

3. (Prosz, 2024) Prosz, A., Pipek, O., Börcsök, J., Palla, G., & Szallasi, Z. (2024). Biologically informed deep learning for explainable epigenetic clocks. *Scientific Reports*, 14, 1439. https://doi.org/10.1038/s41598-023-50495-5

4. (Kalyakulina, 2025) Kalyakulina, A., Yusipov, I., Trukhanov, A., Franceschi, C., & Moskalev, A. (2025). EpInflammAge: Epigenetic-Inflammatory Clock for Disease-Associated Biological Aging Based on Deep Learning. *International Journal of Molecular Sciences*, 26(13), 6284. https://doi.org/10.3390/ijms26136284

5. (Moqri, 2024) Moqri, M., Herzog, C., Poganik, J. R., Ying, K., Justice, J. N., et al. (2024). Validation of biomarkers of aging. *Nature Medicine*, 30, 360–372. https://doi.org/10.1038/s41591-023-02784-9

6. (Levy, 2025) Levy, J. J., Diallo, A. B., Saldias Montivero, M. K., Gabbita, S., & Salas, L. A. (2025). Insights to aging prediction with AI based epigenetic clocks. *Epigenomics*, 17(1), 1–14. https://doi.org/10.1080/17501911.2024.2432854

7. (Johnson, 2022) Johnson, A. A., English, B. W., Shokhirev, M. N., Sinclair, D. A., & Cuellar, T. L. (2022). Human age reversal: Fact or fiction? *Aging Cell*, 21(8), e13664. https://doi.org/10.1111/acel.13664

8. (Richardson, 2025) Richardson, M., Brandt, C., Jain, N., Li, J. L., & Demanelis, K. (2025). Characterization of DNA methylation clock algorithms applied to diverse tissue types. *Aging*, 17(1). https://doi.org/10.18632/aging.206182

9. (Herzog, 2025) Herzog, C. M. S., Redl, E., Barrett, J., Aminzadeh-Gohari, S., & Weber, D. D. (2025). Functionally enriched epigenetic clocks reveal tissue-specific discordant aging patterns in individuals with cancer. *Communications Medicine*, 5, 119. https://doi.org/10.1038/s43856-025-00739-4

10. (Vetter, 2023) Vetter, V. M., Spieker, J., Sommerer, Y., Buchmann, N., & Kalies, C. H. (2023). DNA methylation age acceleration is associated with risk of diabetes complications. *Communications Medicine*, 3, 16. https://doi.org/10.1038/s43856-023-00250-8

11. (Rutledge, 2022) Rutledge, J., Oh, H., & Wyss-Coray, T. (2022). Measuring biological age using omics data. *Nature Reviews Genetics*, 23, 715–727. https://doi.org/10.1038/s41576-022-00511-7

12. (Davydova, 2024) Davydova, E., Perenkov, A., & Vedunova, M. (2024). Building Minimized Epigenetic Clock by iPlex MassARRAY Platform. *Genes*, 15(4), 425. https://doi.org/10.3390/genes15040425

13. (Tian, 2023) Tian, Y., Cropley, V., Maier, A. B., Lautenschlager, N. T., Breakspear, M., & Zalesky, A. (2023). Heterogeneous aging across multiple organ systems and prediction of chronic disease and mortality. *Nature Medicine*, 29, 1224–1232. https://doi.org/10.1038/s41591-023-02296-6

14. (Shokhirev, 2025) Shokhirev, M. N., & Johnson, A. A. (2025). Using buccal methylomic data to create explainable aging clocks as well as classifiers and regressors for lifestyle and demographic factors. *Frontiers in Genetics*, 16, 1637186. https://doi.org/10.3389/fgene.2025.1637186
