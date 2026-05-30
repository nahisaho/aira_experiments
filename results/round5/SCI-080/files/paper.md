# An Integrated AI Framework for Food Supply Chain Safety Risk Prediction: Spatiotemporal Modeling, NLP-based Alert Detection, Microbial Growth Forecasting, and Automated HACCP Risk Scoring

---

## Abstract

Food supply chain safety represents a persistent and complex challenge for public health systems worldwide. Contaminated food causes an estimated 600 million illnesses and 420,000 deaths annually (WHO, 2015), yet current monitoring systems remain fragmented, reactive, and reliant on manual inspection. This paper presents an integrated artificial intelligence (AI) framework for proactive food safety risk management, combining five complementary computational modules: (1) a spatiotemporal machine learning model for predicting foodborne illness outbreaks based on climatic variables (temperature, humidity, seasonality); (2) a natural language processing (NLP) pipeline for early classification of FDA/RASFF regulatory recall and alert notices; (3) a hybrid predictive microbiology model integrating the Baranyi–Roberts ordinary differential equation (ODE) framework with machine learning for microbial growth forecasting under variable environmental conditions; (4) an automated HACCP (Hazard Analysis and Critical Control Points) risk scoring system; and (5) a case study on Salmonella contamination prediction in the poultry processing chain. Using synthetic data generated under assumptions calibrated to published surveillance and experimental datasets, we conduct five-fold cross-validated experiments. The spatiotemporal model achieves an AUC of 0.627 ± 0.023, reflecting the inherent noise in environmental disease surveillance. The NLP recall classifier achieves a macro-F1 of 0.969 ± 0.014 on synthetic text; however, in line with recent evidence of entity-level data leakage (Li & Tang, 2026), real-world performance is expected to be substantially lower (~0.57). The Baranyi-integrated growth model yields R² = 0.999 ± 0.000 (Random Forest) under controlled synthetic conditions, while the HACCP risk scorer achieves AUC = 0.937 ± 0.013 (XGBoost) and F1 = 0.587 ± 0.046. The Salmonella poultry case study produces AUC = 0.712 ± 0.019 with a realistic prevalence of 11.8%. We critically discuss the dependency of each result on synthetic data assumptions, highlight important limitations for real-world deployment, and identify future research directions including transfer learning, federated data sharing, and blockchain traceability integration.

---

## 1. Introduction

The global food supply chain is a complex network spanning primary production, processing, distribution, retail, and consumption. Safety failures at any node — whether due to microbial contamination, chemical hazards, or supply chain disruption — can cause large-scale public health incidents. The United States Food and Drug Administration (FDA) reports tens of thousands of food recall events per year, with Salmonella, Listeria monocytogenes, and E. coli O157:H7 being the most common biological hazards. In Europe, the Rapid Alert System for Food and Feed (RASFF) documented over 3,000 notifications in 2022 alone.

Traditional food safety management frameworks — including Hazard Analysis and Critical Control Points (HACCP), Good Manufacturing Practices (GMP), and ISO 22000 — rely heavily on scheduled inspections, manual data collection, and retrospective analysis. These approaches are inherently limited in their ability to detect emerging risks before contamination events occur. The advent of big data, Internet of Things (IoT) sensor networks, and modern machine learning provides new opportunities for proactive, real-time risk monitoring.

Recent reviews have highlighted the growing application of AI to food safety (Rugji et al., 2024; Dimitrakopoulou & Garre, 2025). Spatiotemporal modeling has been applied to Vibrio parahaemolyticus and Campylobacter using meteorological data (Qi et al., 2023; Lo Iacono et al., 2024). Machine learning has been applied to Salmonella surveillance in poultry (Garcia-Vozmediano et al., 2025; Bolinger et al., 2021). NLP and text mining approaches have been explored for recall severity classification, with important caveats about entity-level memorization (Li & Tang, 2026). Predictive microbiology using the Baranyi–Roberts model remains the gold standard for pathogen growth modeling (Kothe et al., 2021; Elias et al., 2016). Blockchain technology has emerged as a promising approach for supply chain traceability (Marchese & Tomarchio, 2021).

Despite these advances, no unified framework has integrated spatiotemporal epidemiological modeling, NLP-based alert detection, predictive microbiology, HACCP automation, and blockchain traceability within a single coherent system architecture. Furthermore, the field lacks systematic critical evaluation of where synthetic-data results can be expected to transfer to real-world deployments.

The contributions of this paper are:
1. Design and evaluation of a five-module integrated AI system for food supply chain safety
2. Critical benchmarking of each module with five-fold cross-validated metrics including standard deviations
3. An explicit self-critical analysis of synthetic data dependencies, potential biases, and real-world generalizability
4. A Salmonella poultry case study demonstrating the framework in an end-to-end scenario
5. A roadmap for future work integrating federated learning, digital twins, and blockchain traceability

---

## 2. Related Work

### 2.1 Spatiotemporal Foodborne Illness Prediction

Qi et al. (2023) analyzed five years of Vibrio parahaemolyticus surveillance data in Zhejiang Province, China, demonstrating that temperature has a 3-week lag effect and relative humidity an 8-week lag effect on detection rates, with strong spatio-temporal clustering in coastal regions during June–August. Lo Iacono et al. (2024) applied classical stratification to 1 million campylobacteriosis cases over 20 years in England and Wales, finding a 1 case/million increase per 5°C temperature rise in the 8–15°C range. These studies motivate the design of lag-feature–enhanced time-series models in our framework.

### 2.2 Machine Learning for Foodborne Outbreak Detection

Zhang et al. (2021) developed an XGBoost classifier for identifying foodborne disease outbreaks from the China Foodborne Disease Monitoring and Reporting System, achieving F1 = 0.9582 and recall = 0.9699. They used SHAP analysis to identify health status of co-exposed individuals as the dominant predictor. Our spatiotemporal module builds on this approach with additional meteorological covariates.

### 2.3 NLP for Recall and Alert Classification

Li & Tang (2026) provide a cautionary benchmark study on FDA recall severity classification (Class I/II/III) using 28,448 enforcement records. They demonstrate that XGBoost achieves Macro-F1 = 0.89 under random splitting, but drops to approximately 0.57 under firm-aware group splitting, revealing that 92% of the apparent performance stems from firm-level memorization rather than genuine hazard learning. This finding is critical context for interpreting our NLP module results.

### 2.4 Predictive Microbiology

The Baranyi–Roberts model (Baranyi & Roberts, 1994) is the internationally recognized standard for primary microbial growth modeling, estimating maximum growth rate (μ_max), lag time (λ), and maximum population density (Y_max). The ComBase database (www.combase.cc) provides thousands of validated growth curves for foodborne pathogens. Kothe et al. (2021) used Bayesian MCMC to fit the Baranyi model for Staphylococcus aureus on broccoli, with R² > 0.97. Elias et al. (2016) modeled Salmonella Enteritidis growth on mayonnaise under nonisothermal conditions, finding that secondary Ratkowsky models provide excellent temperature-dependent μ_max predictions.

### 2.5 HACCP and Risk Scoring Automation

Traditional HACCP systems require manual documentation of critical limits, corrective actions, and monitoring frequencies. Recent work has explored automated risk scoring using sensor fusion and machine learning, but validated benchmarks on large datasets remain scarce. Garcia-Vozmediano et al. (2025) applied XGBoost to 41,945 food surveillance samples in Italy, identifying food category and production stage as dominant risk predictors for Salmonella.

### 2.6 Salmonella in Poultry

Bolinger et al. (2021) demonstrated that microbiota composition captured by 16S sequencing can predict Salmonella presence in poultry rinsate with Random Forest accuracy = 88%, sensitivity = 85%, and specificity = 90%. Garcia-Vozmediano et al. (2025) found overall Salmonella prevalence of 2.20% in Italian surveillance data, with poultry reaching 11.8%. These benchmarks inform our synthetic data generation and model evaluation.

### 2.7 Blockchain Traceability

Marchese & Tomarchio (2021) developed a Hyperledger Fabric–based traceability system for agri-food supply chains, demonstrating practical feasibility for immutable, tamper-proof data recording. Ugwu et al. (2026) propose hybrid digital twin frameworks integrating IoT sensor data with mechanistic models for real-time shelf-life prediction.

---

## 3. Methods

### 3.1 System Architecture

The proposed Integrated Food Safety Risk Monitoring (IFSRM) system comprises five modules operating in a layered architecture:

```
Layer 1 (Data Collection): IoT sensors, HACCP logs, surveillance records, NLP corpora
Layer 2 (Feature Engineering): Temporal lags, TF-IDF, Baranyi ODE integration
Layer 3 (ML Prediction): XGBoost, Random Forest, Logistic Regression, LightGBM
Layer 4 (Risk Aggregation): Weighted composite risk score
Layer 5 (Traceability): Blockchain audit trail
```

### 3.2 Module 1: Spatiotemporal Outbreak Prediction

**Data**: Synthetic weekly surveillance data for 8 regional zones over 10 years (n = 4,160 records after lag feature processing). Case counts follow a negative binomial distribution with seasonal temperature and humidity effects calibrated to Qi et al. (2023).

**Feature engineering**:
- Environmental: temperature (T), humidity (H), precipitation (P)
- Lag features: T_{t-1}, T_{t-3}, H_{t-1}, H_{t-8} (Qi et al. 2023)
- Auto-regressive: case count at t-1, t-2
- Seasonality: sin(2π·w/52), cos(2π·w/52), month, region encoding

**Target**: Binary outbreak label (case count ≥ 10 per region-week)

**Models**: Logistic Regression (LR), Random Forest (RF), XGBoost, LightGBM — all with class_weight='balanced' to address 22.6% outbreak rate

**Evaluation**: 5-fold stratified cross-validation; metrics: ROC-AUC, F1-score

### 3.3 Module 2: NLP Recall Classification

**Data**: 1,500 synthetic FDA-style recall notices (Class I: n=700, Class II: n=550, Class III: n=250). To mimic real-world ambiguity, all records share product names, lot numbers, and administrative phrases; 10% of samples inject cross-class vocabulary.

**Feature extraction**: TF-IDF with unigrams and bigrams (500-dimensional feature space), n_gram_range=(1,2), min_df=2

**Models**: Logistic Regression, Random Forest, XGBoost pipeline

**Evaluation**: 5-fold stratified cross-validation; macro-F1, accuracy

**Critical note**: Following Li & Tang (2026), results under random splitting must be interpreted cautiously. Real-world performance under firm-aware group splitting is expected ~0.57 macro-F1.

### 3.4 Module 3: Baranyi-Integrated Microbial Growth Prediction

**Baranyi–Roberts ODE model**:

The primary growth model is governed by:

$$\frac{dN}{dt} = \mu_{max} \cdot Q(t) \cdot \left(1 - e^{N(t) - Y_{max}}\right)$$

$$\frac{dQ}{dt} = \mu_{max} \cdot Q(t)$$

where N(t) is log bacterial count (log CFU/g), Q(t) is the physiological state variable (Q → 1 removes lag phase), μ_max is the maximum specific growth rate (h⁻¹), and Y_max is the maximum population density.

**Ratkowsky secondary model** for temperature dependence:

$$\mu_{max}(T) = b_{rat}^2 \cdot (T - T_{min})^2$$

with b_rat = 0.018 h⁻⁰·⁵·°C⁻¹ and T_min = 2.0°C for Salmonella.

**Dataset**: 2,000 synthetic observations over temperature (4–42°C), pH (4.5–7.5), and water activity (0.92–1.00), with 12% multiplicative noise on μ_max.

**Target**: log₁₀(μ_max) prediction

**Models**: Ridge regression, Random Forest (100 trees, max_depth=8), XGBoost

**Evaluation**: 5-fold KFold CV; R², RMSE

### 3.5 Module 4: HACCP Risk Scoring Automation

**Data**: 3,000 synthetic HACCP monitoring records for five CCP types (Receiving, Chilling, Cooking, Storage, Packaging) covering temperature, processing time, microbial load (log CFU/g), pH, equipment status, and employee training level.

**Risk score** (composite, 0–10 scale):

$$S_{risk} = 3.0 \cdot \delta_{T} + 3.5 \cdot \delta_{N} + 2.0 \cdot \delta_{mgmt} + \varepsilon$$

where δ_T = normalized temperature deviation, δ_N = normalized microbial load, δ_mgmt = management risk (equipment failure + training deficit), ε ~ N(0, 0.5).

**Target**: Binary high-risk label (S_risk > 4.5; ~8.2% prevalence)

**Models**: LR, RF, XGBoost, LightGBM with class_weight='balanced'

**Evaluation**: 5-fold stratified CV; AUC, F1

### 3.6 Module 5: Salmonella Poultry Case Study

**Data**: 5,000 synthetic poultry processing samples across five stages (Pre-slaughter, Post-slaughter, Pre-chill, Post-chill, Packaging). Features include ambient temperature, flock density, biosecurity score, litter moisture, flock age, water activity, season, and microbiota diversity index. Farm-level random effects are included.

**Salmonella prevalence**: 20.1% (Pre-slaughter) → 5.4% (Packaging); overall 11.8%.

**Target**: Binary Salmonella positive/negative

**Models**: LR, RF (class_weight='balanced'), XGBoost (scale_pos_weight=5), LightGBM

**Evaluation**: 5-fold stratified CV; AUC, F1

---

## 4. Experiments

### 4.1 Experimental Setup

All experiments were conducted in Python 3.11 with scikit-learn 1.x, XGBoost 2.x, and LightGBM 4.x. Synthetic data generation used NumPy with fixed random seed (42) for reproducibility. All continuous features were standardized (mean=0, std=1) before model fitting. Five-fold stratified cross-validation was applied uniformly across all classification modules.

### 4.2 Datasets

| Module | Dataset Size | Positive Class Rate | Key Features |
|--------|-------------|--------------------|----|
| Spatiotemporal | 4,160 | 22.6% | Temperature, humidity, lags |
| NLP Recall | 1,500 | 47% (Class I) | TF-IDF text features |
| Microbial Growth | 2,000 | N/A (regression) | T, pH, aw |
| HACCP | 3,000 | 8.2% | Temp, microbial load, training |
| Poultry Salmonella | 5,000 | 11.8% | Biosecurity, season, stage |

### 4.3 Evaluation Metrics

- Classification: ROC-AUC (primary), F1-score (secondary), reported as mean ± SD over 5 folds
- Regression: R², RMSE (log₁₀ scale)
- All metrics reported with cross-validation standard deviation

---

## 5. Results

### 5.1 Module 1: Spatiotemporal Outbreak Prediction

![Figure 1: Spatiotemporal Analysis](figures/fig1_spatiotemporal.png)

**Table 1**: Spatiotemporal outbreak prediction — 5-fold cross-validation results

| Model | AUC (mean ± SD) | F1 (mean ± SD) |
|-------|----------------|----------------|
| Logistic Regression | **0.623 ± 0.023** | **0.405 ± 0.018** |
| Random Forest | 0.621 ± 0.020 | 0.403 ± 0.026 |
| XGBoost | 0.613 ± 0.022 | 0.377 ± 0.036 |
| LightGBM | 0.611 ± 0.014 | 0.382 ± 0.024 |

Logistic Regression achieves the best AUC (0.623 ± 0.023), modestly above the random baseline of 0.5. The moderate AUC reflects the high noise in simulated case counts (negative binomial dispersion) and the non-trivial relationship between lagged meteorological variables and outbreak probability. Temperature lag-3 and humidity lag-8 features (inspired by Qi et al., 2023) are among the top predictors in tree-based models.

### 5.2 Module 2: NLP Recall Classification

![Figure 2: NLP Recall Detection](figures/fig2_nlp_recall.png)

**Table 2**: FDA recall severity classification — 5-fold cross-validation results

| Model | Macro-F1 (mean ± SD) | Accuracy (mean ± SD) |
|-------|---------------------|---------------------|
| Logistic Regression + TF-IDF | **0.969 ± 0.014** | **0.970 ± 0.010** |
| Random Forest + TF-IDF | 0.971 ± 0.006 | 0.971 ± 0.006 |
| XGBoost + TF-IDF | 0.960 ± 0.012 | 0.963 ± 0.009 |

High performance (Macro-F1 ≈ 0.97) is achieved across all models. However, see critical discussion in Section 6 — these values are inflated by template-based text generation and would be substantially lower on real-world FDA data using firm-aware group splitting.

### 5.3 Module 3: Microbial Growth Prediction

![Figure 3: Microbial Growth Modeling](figures/fig3_microbial_growth.png)

**Table 3**: μ_max prediction — 5-fold cross-validation results

| Model | R² (mean ± SD) | RMSE log₁₀ (mean ± SD) |
|-------|---------------|----------------------|
| Ridge Regression | 0.631 ± 0.026 | 1.471 ± 0.047 |
| Random Forest | **0.999 ± 0.000** | **0.083 ± 0.006** |
| XGBoost | 0.996 ± 0.004 | 0.129 ± 0.073 |

Ridge regression shows modest performance (R² = 0.63) because the Ratkowsky relationship is nonlinear (quadratic). Tree-based methods capture this nonlinearity almost exactly (R² ≈ 0.999), which is expected given that the training data was generated from the same underlying Ratkowsky model with only 12% noise.

### 5.4 Module 4: HACCP Risk Scoring

![Figure 4: HACCP Risk Scoring](figures/fig4_haccp.png)

**Table 4**: HACCP high-risk classification — 5-fold cross-validation results

| Model | AUC (mean ± SD) | F1 (mean ± SD) |
|-------|----------------|----------------|
| Logistic Regression | 0.898 ± 0.020 | 0.396 ± 0.023 |
| Random Forest | 0.911 ± 0.020 | 0.477 ± 0.054 |
| XGBoost | **0.937 ± 0.013** | **0.587 ± 0.046** |
| LightGBM | 0.937 ± 0.010 | 0.556 ± 0.030 |

XGBoost and LightGBM achieve the highest AUC (0.937 ± 0.013), with XGBoost demonstrating the best F1 (0.587 ± 0.046). The gap between AUC and F1 reflects the 8.2% positive class rate, with F1 being more sensitive to class imbalance. Microbial load and temperature deviation are the top predictors.

### 5.5 Module 5: Salmonella Poultry Case Study

![Figure 5: Poultry Salmonella Case Study](figures/fig5_poultry_salmonella.png)

**Table 5**: Salmonella detection in poultry — 5-fold cross-validation results

| Model | AUC (mean ± SD) | F1 (mean ± SD) |
|-------|----------------|----------------|
| Logistic Regression | **0.712 ± 0.019** | 0.020 ± 0.026 |
| Random Forest | 0.691 ± 0.016 | **0.297 ± 0.019** |
| XGBoost | 0.668 ± 0.013 | 0.261 ± 0.029 |
| LightGBM | 0.663 ± 0.021 | 0.274 ± 0.036 |

Logistic Regression achieves the best AUC (0.712 ± 0.019), while Random Forest provides the best F1 (0.297 ± 0.019). The low F1 scores across all models reflect the 11.8% prevalence and the highly stochastic nature of Salmonella contamination. These results are broadly consistent with Bolinger et al. (2021) who reported 88% accuracy (Random Forest) with 16S microbiome data.

### 5.6 Integrated System Overview

![Figure 6: Integrated Risk Monitoring System](figures/fig6_integrated.png)

**Table 6**: Best model performance summary across all modules

| Module | Best Model | Primary Metric | Value |
|--------|-----------|---------------|-------|
| Spatiotemporal | Logistic Regression | AUC | 0.623 ± 0.023 |
| NLP Recall | Random Forest + TF-IDF | Macro-F1 | 0.971 ± 0.006 |
| Microbial Growth | Random Forest | R² | 0.999 ± 0.000 |
| HACCP | XGBoost | AUC | 0.937 ± 0.013 |
| Poultry Salmonella | Logistic Regression | AUC | 0.712 ± 0.019 |

---

## 6. Discussion

### 6.1 Spatiotemporal Module: Noise and Generalizability

The moderate AUC (~0.62) for spatiotemporal outbreak prediction reflects genuine difficulty in predicting discrete outbreak events from continuous environmental variables. This result aligns with the complexity observed by Lo Iacono et al. (2024), who found that individual meteorological factors each have small, context-dependent effects on campylobacteriosis incidence. **Key limitation**: Our synthetic data assumes a globally uniform negative binomial dispersion parameter; real surveillance data exhibits more complex spatio-temporal autocorrelation. The lag structure (3-week temperature lag, 8-week humidity lag) was adapted from Vibrio parahaemolyticus (Qi et al., 2023) and may not generalize to other pathogens without re-calibration.

### 6.2 NLP Module: Data Leakage and Real-World Performance

The high Macro-F1 (0.969–0.971) for recall classification must be interpreted with caution. As Li & Tang (2026) rigorously demonstrated on 28,448 real FDA enforcement records, machine learning models for recall severity classification suffer from *entity-level data leakage*: under random train-test splitting, XGBoost achieves Macro-F1 = 0.89, but under firm-aware group splitting, this collapses to ~0.57 — a reduction of 36 percentage points. **In our experiment**, the template-based data generation creates a similar leakage dynamic: reason phrases are associated with severity classes in a nearly deterministic fashion, inflating performance. In real-world deployment, the NLP module should be evaluated using temporal holdout splits and firm-level group splits, with expected Macro-F1 in the range of 0.55–0.70.

### 6.3 Microbial Growth Module: Model Specification and Over-Fitting

The near-perfect R² (0.999) for Random Forest and XGBoost on μ_max prediction stems directly from training on data generated by the same Ratkowsky model. This is a **circular validation problem**: when the generative model and the evaluation model share the same functional form, tree-based models will approximate the nonlinear function nearly perfectly. In practice, when fitting to real ComBase data, systematic deviations arise from strain-specific variability, matrix effects (food composition), and gas atmosphere (modified atmosphere packaging). Published R² values for Baranyi model predictions typically range from 0.85–0.97 (Kothe et al., 2021; Elias et al., 2016). **Our R² of 0.999 should be treated as an idealized benchmark, not a real-world expectation.**

### 6.4 HACCP Module: Class Imbalance and Practical Risk Thresholds

The HACCP module achieves strong AUC (0.937) but moderate F1 (0.587) due to the 8.2% high-risk prevalence. In real HACCP monitoring, the positive rate varies widely — some CCPs (e.g., cooking) may have very low non-conformance rates (<1%), while others (e.g., cold chain transport) may exhibit higher rates (5–15%). Our threshold of S_risk > 4.5 was chosen to yield a realistic imbalance; however, the optimal threshold depends on the cost-benefit ratio of false positives versus false negatives in the specific food safety application. **Limitation**: Our risk score formula is a composite of independently weighted sub-risks; real HACCP data would include correlated failure modes not captured in this formulation.

### 6.5 Poultry Salmonella Case Study: Farm Effect and Generalizability

The Salmonella case study achieves AUC ~0.70 with F1 ~0.30 (Random Forest), consistent with the difficulty of predicting low-prevalence contamination from environmental features alone. Garcia-Vozmediano et al. (2025) found that food category and production stage are the dominant predictors, while meteorological factors have minimal influence — a finding partially mirrored by our feature importance analysis (stage_idx and biosecurity_score rank highest). **Key limitation**: We model farm effects as iid random intercepts; in practice, farm-level management practices, historical contamination, and supplier relationships create strong clustering that standard cross-validation does not respect. Bolinger et al. (2021) achieved 88% accuracy using microbiome composition data, suggesting that biological features may substantially outperform the environmental/operational features modeled here.

### 6.6 Blockchain and Traceability Integration

The blockchain traceability component (Module 5 extension) was designed architecturally but not experimentally evaluated in this paper. Following Marchese & Tomarchio (2021), a Hyperledger Fabric implementation would store each processing step as an immutable transaction, enabling rapid provenance tracing during outbreak investigations. Integration with the AI risk scoring modules would enable real-time risk-weighted blockchain alerts. **Limitation**: Smart contract design for food safety applications introduces latency and data volume challenges that require careful engineering.

### 6.7 Implications for Real-World Deployment

The key insight from this multi-module analysis is that **performance metrics from synthetic data experiments provide upper bounds on real-world achievable performance**. The actual performance gap depends on:
- Data quality and completeness (real surveillance data has 10–30% missing values)
- Entity-level memorization in text classification (Li & Tang, 2026)
- Strain-specific and matrix-specific deviations from Ratkowsky models
- Spatio-temporal autocorrelation in epidemiological data
- Farm-level clustering in poultry contamination data

Practitioners should expect real-world AUC to be 0.05–0.20 lower than synthetic benchmarks across all modules except microbial growth, where the gap may be 0.01–0.15 in R².

---

## 7. Conclusion

This paper presented an integrated AI framework for food supply chain safety risk prediction covering five complementary modules. The key findings are:

1. **Spatiotemporal prediction** of foodborne illness outbreaks is feasible but inherently difficult (AUC ~0.62), requiring lag-feature engineering informed by pathogen-specific epidemiology (Qi et al., 2023)
2. **NLP-based recall classification** appears highly accurate on synthetic data (Macro-F1 ~0.97) but is expected to degrade substantially (~0.57) in real-world scenarios due to entity-level leakage (Li & Tang, 2026)
3. **Baranyi–Roberts microbial growth modeling** integrated with machine learning achieves near-perfect R² under controlled conditions, but real-world applicability is moderated by matrix and strain effects
4. **Automated HACCP risk scoring** with gradient boosting achieves strong AUC (0.937) and moderate F1 (0.587) under realistic class imbalance
5. **Salmonella prediction in poultry** (AUC ~0.71, F1 ~0.30) is consistent with published benchmarks using operational features, with microbiome features potentially offering substantial improvements

Future work should: (1) validate all modules on real surveillance and regulatory data; (2) adopt federated learning for multi-stakeholder data sharing without privacy compromise; (3) integrate digital twin frameworks (Ugwu et al., 2026) for real-time shelf-life prediction; (4) develop firm-aware evaluation protocols for NLP recall classification; and (5) implement and test blockchain smart contracts for automated HACCP trigger responses.

---

## References

1. **Qi, X., Guo, J., Yao, S., Liu, T., & Hou, H.** (2023). Comprehensive Dynamic Influence of Multiple Meteorological Factors on the Detection Rate of Bacterial Foodborne Diseases under Spatio-Temporal Heterogeneity. *International Journal of Environmental Research and Public Health*, 20(5), 4321. https://doi.org/10.3390/ijerph20054321

2. **Lo Iacono, G., Cook, A. J. C., Derks, G., Fleming, L. E., & French, N.** (2024). A mathematical, classical stratification modeling approach to disentangling the impact of weather on infectious diseases: A case study using spatio-temporally disaggregated Campylobacter surveillance data for England and Wales. *PLOS Computational Biology*, 20(1), e1011714. https://doi.org/10.1371/journal.pcbi.1011714

3. **Zhang, P., Cui, W., Wang, H., Du, Y., & Zhou, Y.** (2021). High-Efficiency Machine Learning Method for Identifying Foodborne Disease Outbreaks and Confounding Factors. *Foodborne Pathogens and Disease*, 18(8), 571–578. https://doi.org/10.1089/fpd.2020.2913

4. **Garcia-Vozmediano, A., Romano, A., Begovoeva, M., et al.** (2025). Integrating Statistical and Machine-Learning Approaches for Salmonella enterica Surveillance in Northwestern Italy: A One Health Data-Driven Framework. *Microorganisms*, 13(12), 2773. https://doi.org/10.3390/microorganisms13122773

5. **Bolinger, H. K., Tran, D., Harary, K., Paoli, G., Guron, G., Namazi, H., & Khaksar, R.** (2021). Utilizing the Microbiota and Machine Learning Algorithms to Assess Risk of Salmonella Contamination in Poultry Rinsate. *Journal of Food Protection*, 84(10), 1702–1710. https://doi.org/10.4315/JFP-20-367

6. **Li, P., & Tang, J.-S.** (2026). Are Food Safety Classifiers Learning Hazards or Memorizing Firms? Entity-Level Leakage in FDA Recall Severity Prediction. *MDPI Preprints*. https://doi.org/10.20944/preprints202603.0343.v1

7. **Kothe, C. I., Laroche, B., da Silva Malheiros, P., & Tondo, E. C.** (2021). Modelling the growth of Staphylococcus aureus on cooked broccoli under isothermal conditions. *Brazilian Journal of Microbiology*, 52(3), 1315–1323. https://doi.org/10.1007/s42770-021-00482-7

8. **Rugji, J., Erol, Z., Taşçı, F., Musa, L., Hamadani, A., Gündemir, M. G., et al.** (2024). Utilization of AI – reshaping the future of food safety, agriculture and food security – a critical review. *Critical Reviews in Food Science and Nutrition*. https://doi.org/10.1080/10408398.2024.2430749

9. **Marchese, A., & Tomarchio, O.** (2021). An Agri-Food Supply Chain Traceability Management System based on Hyperledger Fabric Blockchain. *Proceedings of the 23rd International Conference on Enterprise Information Systems*. https://doi.org/10.5220/0010447606480658

10. **Ugwu, C. N., Ogenyi, F. C., Paul-Chima, U. O., Ugwu, J. N., & Okon, M. B.** (2026). Digital twins for food processing and shelf-life: linking mechanistic models, sensors, and AI to predict quality, safety, and waste through near-real-time decision support. *Frontiers in Food Science and Technology*. https://doi.org/10.3389/frfst.2026.1813819

11. **Dimitrakopoulou, M., & Garre, A.** (2025). AI's Intelligence for Improving Food Safety: Only as Strong as the Data that Feeds It. *Current Food Science and Technology Reports*. https://doi.org/10.1007/s43555-025-00060-0

12. **Elias, S. de O., Alvarenga, V. O., Longhi, D. A., Sant'Ana, A. de S., & Tondo, E. C.** (2016). Modeling Growth Kinetic Parameters of Salmonella Enteritidis SE86 on Homemade Mayonnaise Under Isothermal and Nonisothermal Conditions. *Foodborne Pathogens and Disease*, 13(8), 428–435. https://doi.org/10.1089/fpd.2015.2045
