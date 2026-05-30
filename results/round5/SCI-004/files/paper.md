# Integrative Pharmacogenomics Modeling: From CYP Enzyme Polymorphisms to Clinical Decision Support via Machine Learning and Mendelian Randomization

---

## Abstract

Pharmacogenomics—the study of how genetic variation shapes drug response—holds transformative potential for precision medicine, yet its systematic computational modeling across multiple biological layers remains underexplored. This work presents an integrated computational framework addressing six interconnected pharmacogenomics challenges: (1) CYP2D6/CYP2C19 metabolizer phenotype prediction and plasma concentration modeling, (2) HLA-B\*1502–mediated carbamazepine adverse drug reaction (ADR) prediction, (3) Mendelian randomization (MR)–based drug target validation using simulated GWAS summary statistics, (4) anticancer drug sensitivity prediction using GDSC/CCLE-like genomic features, (5) deep learning drug–gene interaction network learning, and (6) prototype evaluation of a pharmacogenomics-guided clinical decision support system (CDSS). Using synthetic datasets with realistic noise calibrated to published epidemiological parameters, we demonstrate that CYP phenotype classification achieves high but potentially overfit accuracy (RF: 1.000±0.000, SVM: 0.980±0.011) when activity-score features closely approximate the label-generation function—a critical methodological caveat. HLA-mediated ADR prediction achieves moderate AUROC (0.681±0.047) under realistic class imbalance (10.1% prevalence), consistent with real-world observations. MR analysis recovers causal effects close to the ground truth (IVW: β=−0.366, SE=0.001 vs. true β=−0.35), though overly precise SE reflects clean synthetic data. Drug sensitivity prediction attains R²=0.038–0.347 across eight drugs—comparable to published GDSC benchmarks. Drug–gene interaction learning from feature-pair representations yields near-random AUROC (≈0.49–0.50), reflecting insufficient structural information in the representation. CDSS simulation shows a 40.7% relative efficacy improvement (0.536→0.754, p<0.001) and 35.4% ADR reduction (0.158→0.102, p=0.011). These findings highlight both the promise and the substantial methodological pitfalls of pharmacogenomics machine learning, particularly overfitting, data leakage, and limited generalizability to unseen real-world populations. We provide a reproducible evaluation pipeline and discuss translational requirements for clinical implementation.

**Keywords**: pharmacogenomics, CYP2D6, HLA-B\*1502, Mendelian randomization, drug sensitivity, deep learning, CDSS, precision medicine

---

## 1. Introduction

Interindividual variability in drug efficacy and toxicity is a fundamental challenge in clinical pharmacology. A substantial fraction of this variability—estimated at 34–98% for some drugs—has a heritable genetic basis [Zhou et al., 2022]. Pharmacogenomics exploits this genetic basis to predict drug response at the individual level, enabling genotype-guided prescribing. Key pharmacogenes include the cytochrome P450 enzymes (notably CYP2D6 and CYP2C19, collectively involved in the metabolism of ~25% of clinically used drugs), human leukocyte antigen (HLA) genes mediating immune-mediated hypersensitivity reactions, and germline variants modulating tumor drug sensitivity.

Despite decades of research, clinical implementation of pharmacogenomics remains limited. A key bottleneck is the computational framework for integrating heterogeneous genetic data into actionable predictions [Taylor et al., 2020; Sperber et al., 2021]. Machine learning methods offer a promising path but require rigorous validation to avoid common pitfalls: data leakage from phenotype-derived features, inflated performance from artificial datasets, and poor generalizability across ancestries [Zhou et al., 2022].

This paper makes the following contributions:

1. **A multi-task computational framework** spanning CYP metabolism, HLA immunogenetics, causal inference (MR), cancer pharmacogenomics, and CDSS design.
2. **Critical self-evaluation**: for each module, we explicitly analyze whether performance is genuine or an artifact of simulation assumptions.
3. **Calibrated synthetic datasets** with realistic noise levels and class prevalence based on published epidemiological data.
4. **A reproducible evaluation pipeline** using stratified k-fold cross-validation with standard deviation reporting.

Our work is directly motivated by clinical need: CYP2D6 alone accounts for >95% of pharmacogenomic drug prescribing in primary care [Magavern et al., 2021], and HLA-B\*1502–carbamazepine is mandated for pre-prescription genotyping in many Asian countries [Kloypan et al., 2021].

---

## 2. Related Work

### 2.1 CYP Pharmacogenomics

Taylor et al. (2020) provided a comprehensive review of CYP2D6's role in metabolizing ~20% of commonly used drugs and catalogued its >100 star alleles with differing functional consequences. Machine learning approaches to variant effect prediction have been proposed but remain limited by training data size and structural variant complexity [Zhou et al., 2022]. Lauschke et al. (2022) discussed how rare variants account for a substantial proportion of unexplained heritability in CYP-mediated drug response.

### 2.2 HLA and Adverse Drug Reactions

Kloypan et al. (2021) reviewed the strong genetic association between HLA-B\*1502 and carbamazepine-induced Stevens-Johnson syndrome/toxic epidermal necrolysis (SJS/TEN) in Asian populations, with odds ratios exceeding 1000 in some studies. Predictive models must account for population stratification, as allele frequency and risk magnitude differ substantially across ancestries.

### 2.3 Mendelian Randomization for Drug Target Validation

MR uses genetic variants as instrumental variables (IVs) to establish causal relationships, providing a natural experiment framework analogous to randomized trials [Liu et al., 2025; Ke et al., 2026]. Multi-omics MR integrating GWAS, eQTL, and pQTL data has emerged as a powerful approach for drug target identification, as demonstrated for cardiac arrest [Liu et al., 2025] and IBS [Ke et al., 2026].

### 2.4 Cancer Drug Sensitivity Prediction

Wang et al. (2021) introduced DeepDRK, achieving AUROC=0.84 on patient-derived cell lines by integrating multi-omics features through kernel-based similarity matrices. Park et al. (2023) benchmarked deep learning vs. traditional ML for 24 GDSC drugs and found comparable performance (R²=−8 to 0.47), with ridge regression for panobinostat achieving the best R²=0.470. Taherdoost & Ghofrani (2024) reviewed AI's role across pharmacogenomics tasks, emphasizing the need for explainability.

### 2.5 Clinical Implementation

Sperber et al. (2021) analyzed implementation strategies across IGNITE Network projects, finding that pharmacogenomics programs require more strategies than disease-focused projects. Gill et al. (2021) reported a clinical PGx testing system at Arkansas Children's Hospital integrated with EPIC-EHR, covering 66 pediatric medications via 174 SNPs in 23 genes.

---

## 3. Methods

### 3.1 Module 1: CYP2D6/CYP2C19 Phenotype Prediction

#### 3.1.1 Data Simulation

We simulated N=1,200 patients with CYP2D6 and CYP2C19 star allele pairs drawn from population-frequency distributions consistent with European ancestry. CYP2D6 star alleles included \*1 (normal, f=0.35), \*2 (normal, 0.15), \*4 (non-functional, 0.20), \*5 (deletion, 0.05), \*10 (reduced, 0.10), \*41 (reduced, 0.10), and xN (duplication/gain, 0.05). An activity score was computed as:

$$AS_{CYP2D6} = \text{score}(A_1) + \text{score}(A_2)$$

where scores map alleles to activity values: non-functional→0, reduced→0.5, normal→1.0, duplicated→2.0. Metabolizer phenotype was defined by thresholding: AS≤0.5 = Poor Metabolizer (PM), 0.5<AS≤1.25 = Intermediate (IM), 1.25<AS≤2.0 = Normal (NM), AS>2.0 = Ultra-rapid (UM).

#### 3.1.2 Drug Plasma Concentration Modeling

Warfarin/risperidone-like Cmax was simulated as:

$$C_{max} = \text{PhenoFactor} \times \frac{D}{100} \times (1 + 0.5 \cdot I_{inh}) \times (1 - 0.003(Age - 50)) \times e^{\epsilon}$$

where $\epsilon \sim \mathcal{N}(0, 0.35)$ captures lognormal variability, and $I_{inh}$ is a co-medication inhibitor indicator.

#### 3.1.3 Feature Matrix

Features: encoded star allele identifiers (label-encoded), activity scores, age, weight, dose, inhibitor status, specific null-allele indicators. StandardScaler normalization was applied.

#### 3.1.4 Models and Evaluation

Classifiers: Random Forest (200 trees, depth=8), Logistic Regression (C=1), MLP(64,32), SVM-RBF. Regression: Random Forest, Ridge (α=1), Gradient Boosting (lr=0.05, 200 trees). 5-fold stratified/k-fold CV with accuracy/R²/RMSE metrics.

### 3.2 Module 2: HLA-B\*1502 ADR Prediction

#### 3.2.1 Data Simulation

N=1,500 patients. HLA-B\*1502 carrier frequency: Asian=8%, European=1%, African=0.5%, Other=3%. ADR (SJS/TEN) logistic model:

$$\text{logit}(P_{ADR}) = -4.0 + 4.5 \cdot HLA_{1502} + 0.8 \cdot Eth_{Asian} \cdot HLA_{1502} + 0.5 \cdot HLA_{A3101} + 0.3 \cdot HLA_{B5801} + 0.002 \cdot D_{CBZ} + \epsilon$$

#### 3.2.2 Evaluation

Stratified 5-fold CV. Metrics: AUROC, F1 (class-balanced). Class imbalance handled with class_weight='balanced'.

### 3.3 Module 3: Mendelian Randomization

#### 3.3.1 Instrument Variable Generation

N=30 cis-eQTL SNPs for CYP2C9 (exposure: CYP2C9 expression); outcome: warfarin dose requirement. SNP-exposure effects γ̂ ~ N(0.25, 0.05²); outcome effects Γ̂ = β_true · γ̂ + δ_pleiotropic + ε, where β_true = −0.35, δ ~ N(0, 0.03²), ε ~ N(0, 0.04²).

#### 3.3.2 MR Methods

- **IVW (Inverse-Variance Weighted)**:
  $$\hat{\beta}_{IVW} = \frac{\sum_j w_j r_j}{\sum_j w_j}; \quad r_j = \frac{\hat{\Gamma}_j}{\hat{\gamma}_j}; \quad w_j = \frac{\hat{\gamma}_j^2}{\hat{\sigma}_{\Gamma_j}^2}$$
  
- **Weighted Median**: median of ratio estimates weighted by instrument strength.

- **MR-Egger**: weighted linear regression allowing non-zero intercept to test directional pleiotropy.

### 3.4 Module 4: Cancer Drug Sensitivity

#### 3.4.1 Data

N=800 simulated cancer cell lines × 8 drugs. Feature matrix: 200-dimensional gene expression + 50 binary mutation indicators + 30 copy number profiles (total 280 features). IC50 (log-transformed) generated as sparse linear combination of 8–20 gene expression features plus lognormal noise (σ=1.5).

#### 3.4.2 Models

Random Forest Regressor (100 trees, depth=6) with 5-fold CV. Features: top 100 of 280 (standardized).

### 3.5 Module 5: Drug–Gene Interaction Network

#### 3.5.1 Data

N=2,000 drug-gene pairs (20 drugs × 100 genes, randomly sampled). Drug representation: 50-dim ECFP-like binary fingerprint + 10 physicochemical properties. Gene representation: 30-dim feature vector (expression + SNP profiles). Pair features: concatenated 90-dim vector. Interaction label: Bernoulli(sigmoid(true_interaction_score)).

#### 3.5.2 Models

MLP(64,32), MLP(128,64,32), Random Forest, Logistic Regression. 5-fold stratified CV with AUROC and F1.

### 3.6 Module 6: CDSS Simulation

Binary outcome simulation (N=500/group): efficacy (Bernoulli 0.60 vs. 0.75) and ADR rate (0.15 vs. 0.09) for standard vs. PGx-guided care. Statistical testing: χ² test.

---

## 4. Experiments

### 4.1 Experimental Setup

All experiments used Python 3.11 with scikit-learn 1.x. Random seed 42 for reproducibility. 5-fold cross-validation for all supervised learning tasks.

### 4.2 Evaluation Metrics

- **Classification**: Accuracy, AUROC, F1-score (macro/binary)
- **Regression**: R² (coefficient of determination), RMSE
- **Causal inference**: Effect estimate (β), SE, p-value
- **Clinical simulation**: Efficacy rate, ADR rate, χ² p-value

### 4.3 Baseline Comparisons

For each task, results compared against: random baseline (AUROC=0.5, R²=0), single-feature (HLA-B\*1502 only for ADR), and published benchmark values from GDSC literature.

---

## 5. Results

### 5.1 CYP2D6/CYP2C19 Phenotype Prediction

**Table 1: CYP2D6 Metabolizer Phenotype Classification (5-fold CV, N=1,200)**

| Model | Accuracy ± SD | Notes |
|---|---|---|
| Random Forest | **1.000 ± 0.000** | ⚠️ Likely data leakage (see Discussion) |
| Logistic Regression | **1.000 ± 0.000** | ⚠️ Same concern |
| MLP(64,32) | 0.996 ± 0.008 | Near-perfect |
| SVM-RBF | 0.980 ± 0.011 | — |

**Table 2: Drug Plasma Concentration Prediction (log-Cmax, 5-fold CV)**

| Model | R² ± SD | RMSE ± SD |
|---|---|---|
| Gradient Boosting | **0.780 ± 0.020** | 0.235 ± 0.014 |
| Random Forest | 0.768 ± 0.023 | 0.241 ± 0.015 |
| Ridge Regression | 0.738 ± 0.023 | 0.257 ± 0.014 |

![Figure 1: CYP2D6 Phenotype Prediction](figures/exp1_cyp_phenotype.png)

*Figure 1: (Left) Confusion matrix for Random Forest CYP2D6 metabolizer phenotype prediction on training data. (Right) Top feature importances. Note: near-perfect training accuracy reflects feature construction, not generalization.*

### 5.2 HLA-B\*1502 / Carbamazepine ADR Prediction

Dataset characteristics: N=1,500, ADR prevalence=10.1% (HLA-B\*1502 positive: 82.1%, negative: 8.1%).

**Table 3: HLA-ADR Prediction Performance (5-fold CV)**

| Model | AUROC ± SD | F1 ± SD |
|---|---|---|
| Logistic Regression | **0.683 ± 0.061** | 0.289 ± 0.030 |
| Random Forest | 0.681 ± 0.047 | **0.325 ± 0.073** |
| MLP(64,32) | 0.620 ± 0.068 | 0.293 ± 0.076 |
| Random Baseline | 0.500 | ~0.18 |

![Figure 2: HLA-B*1502 ADR Prediction](figures/exp2_hla_adr.png)

*Figure 2: (Left) ROC curve for HLA-ADR prediction. (Right) ADR rate stratified by HLA-B\*1502 status and ethnicity, showing the strong interaction between HLA carrier status and Asian ancestry.*

### 5.3 Mendelian Randomization Drug Target Validation

**Table 4: MR Results – CYP2C9 Expression → Warfarin Dose (True β = −0.35)**

| Method | β Estimate | SE | p-value | Notes |
|---|---|---|---|---|
| IVW | −0.366 | 0.001 | <0.0001 | Consistent with truth |
| Weighted Median | −0.366 | — | — | Robust to pleiotropy |
| MR-Egger | −0.198 | — | — | Intercept = −0.042 |
| True Causal Effect | −0.350 | — | — | Simulation ground truth |

![Figure 3: Mendelian Randomization Analysis](figures/exp3_mr_analysis.png)

*Figure 3: (Left) Forest plot of MR causal effect estimates across three methods. (Right) Scatter plot of SNP–exposure vs. SNP–outcome effects for 30 instrumental variables.*

### 5.4 Cancer Drug Sensitivity Prediction

**Table 5: Drug Sensitivity Prediction (5-fold CV, GDSC-like Data, N=800 cell lines)**

| Drug | R² ± SD | RMSE ± SD |
|---|---|---|
| Imatinib | **0.347 ± 0.056** | 2.321 ± 0.134 |
| Erlotinib | 0.312 ± 0.069 | 2.453 ± 0.233 |
| Olaparib | 0.263 ± 0.034 | 2.192 ± 0.104 |
| Vemurafenib | 0.234 ± 0.080 | 2.339 ± 0.160 |
| Tamoxifen | 0.231 ± 0.057 | 3.572 ± 0.389 |
| Cisplatin | 0.204 ± 0.034 | 2.458 ± 0.074 |
| Docetaxel | 0.200 ± 0.052 | 2.636 ± 0.128 |
| Gemcitabine | 0.038 ± 0.061 | 2.795 ± 0.120 |

![Figure 4: Cancer Drug Sensitivity Prediction](figures/exp4_gdsc_sensitivity.png)

*Figure 4: (Left) Drug sensitivity prediction R² by drug. (Right) Predicted vs. actual log(IC50) for the best-performing drug (Imatinib).*

### 5.5 Drug–Gene Interaction Network Learning

**Table 6: Drug–Gene Interaction Prediction (5-fold CV, N=2,000 pairs)**

| Model | AUROC ± SD | F1 ± SD |
|---|---|---|
| Random Forest | 0.502 ± 0.011 | 0.176 ± 0.029 |
| LR | 0.500 ± 0.022 | 0.245 ± 0.039 |
| MLP(128,64,32) | 0.490 ± 0.023 | 0.387 ± 0.022 |
| MLP(64,32) | 0.490 ± 0.018 | 0.386 ± 0.027 |
| Random Baseline | 0.500 | — |

![Figure 5: Drug-Gene Network Learning](figures/exp5_network_learning.png)

*Figure 5: (Left) AUROC and F1 comparison across models for drug–gene interaction prediction. (Right) Learning curve showing that both training and validation AUROC remain near 0.5 regardless of training set size, indicating insufficient signal in the feature representation.*

### 5.6 CDSS Simulation

**Table 7: CDSS – Standard vs. PGx-Guided Care (N=500)**

| Outcome | Standard Care | PGx-Guided | Δ | p-value |
|---|---|---|---|---|
| Drug Efficacy | 0.536 | 0.754 | +40.7% | <0.001 |
| ADR Rate | 0.158 | 0.102 | −35.4% | 0.011 |

![Figure 6: CDSS Summary](figures/exp6_cdss_summary.png)

*Figure 6: (Left) Clinical outcome comparison for standard vs. PGx-guided care. (Right) Summary of cross-validation performance across all experimental modules.*

---

## 6. Discussion

### 6.1 CYP Phenotype Classification: Data Leakage Analysis

The near-perfect accuracy (1.000±0.000) for Random Forest and Logistic Regression in CYP2D6 phenotype classification is **not evidence of a robust predictive model**. The metabolizer phenotype labels were defined as a deterministic threshold function of the activity score, which was itself computed from the same star allele features used for training. This constitutes a structural data leakage: the features essentially encode the answer. In real-world scenarios:

- Novel or rare star alleles with unknown functional effects would be misclassified
- Post-translational regulation, protein–protein interactions, and environmental factors (e.g., dietary inhibitors) are not captured
- The phenotype-genotype concordance in clinical data is substantially lower (kappa ≈ 0.7–0.85 in validation studies)

**Expected real-world performance**: AUROC ~0.75–0.85 for metabolizer classification, with PM/UM being easier to classify than IM. The Cmax regression results (R²=0.74–0.78) are more credible, as they use realistic noise from lognormal pharmacokinetic variability.

### 6.2 HLA-ADR Prediction: Realistic but Limited

The AUROC of 0.68±0.05 for HLA-mediated ADR prediction is consistent with real-world challenges:

- **Class imbalance** (10.1% prevalence) severely limits F1-score despite reasonable discrimination
- **Ancestry confounding**: HLA-B\*1502 frequency varies ~16-fold between Asian and European populations; ancestry must be included as a covariate
- **Sensitivity vs. specificity trade-off**: In clinical practice, pre-prescription screening requires high sensitivity (>95%) at the cost of many false positives, a regime where AUROC alone is insufficient

Real-world observations from mandated screening in Taiwan showed positive predictive values of only 2–5%, meaning ~95% of carriers screened out would not have developed SJS/TEN. This limits the utility of expanded biomarker panels beyond HLA-B\*1502 alone.

### 6.3 Mendelian Randomization: Overly Precise Simulation

The IVW p-value of effectively 0 and SE=0.001 reflect an unrealistically clean simulation. Real MR studies with 30 IVs typically achieve p~10⁻⁶ to 10⁻¹⁵, not machine zero. The divergence between IVW (β=−0.366) and MR-Egger (β=−0.198) correctly reflects simulated directional pleiotropy (intercept=−0.042), but the absolute magnitude of separation would be larger in real data with stronger pleiotropic pathways. The qualitative conclusion—CYP2C9 expression causally affects warfarin dose requirement—is consistent with well-established CYP2C9\*2/\*3 pharmacogenomics literature.

### 6.4 Drug Sensitivity Prediction: Moderate, Realistic Performance

R²=0.038–0.347 across eight simulated drugs is well-aligned with published GDSC benchmarks. Park et al. (2023) reported R²=−8 to 0.47 across 24 GDSC drugs, and Wang et al.'s DeepDRK (2021) achieved AUROC=0.84 only with multi-omics feature integration. Our single-omics Random Forest approach—using only gene expression, mutation, and copy number—mirrors typical unimodal performance in real data. Notably, Gemcitabine shows near-zero R²=0.038±0.061, consistent with published findings that some chemotherapy agents have weak genomic predictors.

**Key assumption dependency**: Performance depends on how many informative genomic features drive the true signal (here: 8–20 genes). Real GDSC data shows that many drugs have diffuse, polygenic response architecture, limiting predictive accuracy even with large training sets.

### 6.5 Drug–Gene Interaction: Negative Result and Its Informative Value

All models achieved AUROC≈0.49–0.50, effectively random. This is a **meaningful negative result** arising from an important experimental design flaw: the interaction labels were generated using drug-specific weight vectors that are inaccessible to the model. In other words, the concatenated drug fingerprint + gene feature representation does not preserve the interaction structure encoded in the true labeling function. This mirrors a common failure mode in drug–target interaction (DTI) prediction: generic molecular descriptors cannot capture the specific binding geometry or allosteric mechanisms that determine activity.

**Implications**: DTI prediction requires either: (a) co-evolutionary sequence features (e.g., protein language models), (b) 3D structural docking scores, or (c) pre-trained embeddings (e.g., ESM-2 for proteins, ChemBERTa for drugs). Simple concatenation of independent embeddings fails even when the signal is deterministic.

### 6.6 CDSS Simulation: Optimistic Assumptions

The CDSS simulation assumes an idealized setting where: (1) genotype testing is always available and accurate, (2) clinicians perfectly adhere to PGx recommendations, (3) no implementation delays or cost constraints. Real-world implementation studies [Sperber et al., 2021; Gill et al., 2021] report substantially smaller benefits due to:

- Prescriber non-adherence to genomic recommendations (50–70% in early implementation)
- Limited evidence for many drug-gene pairs
- Ethnic diversity in pharmacogene frequencies requiring ancestry-stratified guidelines

### 6.7 Limitations

1. **Synthetic data**: All results are derived from simulated datasets. Validation in cohort studies (e.g., UK Biobank, BioMe Biobank) is required before clinical translation.
2. **Single-ancestry focus**: Population allele frequencies were modeled on European-ancestry distributions; multi-ancestry modeling is essential for equitable implementation.
3. **No longitudinal outcomes**: Real-world treatment outcomes are time-dependent and confounded by sequential therapeutic decisions.
4. **Feature engineering**: Several modules use hand-crafted features (e.g., activity scores) that embed domain knowledge unavailable in de novo discovery settings.
5. **Sample sizes**: N=800–2,000 are small relative to typical GDSC (>900 cell lines × 400 drugs) and biobank datasets (>100,000 participants).

---

## 7. Conclusion

We presented an integrated pharmacogenomics modeling framework addressing CYP enzyme phenotyping, HLA-mediated ADR prediction, MR-based causal inference, cancer drug sensitivity modeling, drug–gene interaction learning, and CDSS prototype evaluation. Our results reveal a spectrum of outcomes: from artificially inflated accuracy due to feature-label coupling (CYP phenotype), through realistic moderate performance (HLA-ADR, drug sensitivity), to a meaningful negative result (drug-gene interaction with insufficient representations). 

Critical methodological lessons include: (1) synthetic accuracy gains from deterministic label construction must be distinguished from genuine predictive modeling; (2) AUROC alone is insufficient for rare-event prediction tasks requiring high sensitivity; (3) multi-omics integration and structural molecular representations are essential for drug–gene interaction learning; (4) MR-based causal inference provides a valuable complement to association-based approaches for target validation.

Future work should focus on: multi-ancestry cohort validation, integration of protein language model embeddings for drug–target representation, prospective clinical trials of PGx-guided CDSS, and federated learning approaches to enable privacy-preserving multi-center training.

---

## References

1. **Taylor CE, Crosby IT, Yip V, Maguire P, Pirmohamed M, Turner RM** (2020). A Review of the Important Role of CYP2D6 in Pharmacogenomics. *Genes*, 11(11):1295. DOI: [10.3390/genes11111295](https://doi.org/10.3390/genes11111295)

2. **Kloypan C, Koomdee N, Satapornpong P, et al.** (2021). A Comprehensive Review of HLA and Severe Cutaneous Adverse Drug Reactions: Implication for Clinical Pharmacogenomics and Precision Medicine. *Pharmaceuticals*, 14(11):1077. DOI: [10.3390/ph14111077](https://doi.org/10.3390/ph14111077)

3. **Wang Y, Yang Y, Chen S, Wang J** (2021). DeepDRK: a deep learning framework for drug repurposing through kernel-based multi-omics integration. *Briefings in Bioinformatics*, 22(5):bbab048. DOI: [10.1093/bib/bbab048](https://doi.org/10.1093/bib/bbab048)

4. **Sperber NR, Dong OM, Roberts MC, et al.** (2021). Strategies to Integrate Genomic Medicine into Clinical Care: Evidence from the IGNITE Network. *Journal of Personalized Medicine*, 11(7):647. DOI: [10.3390/jpm11070647](https://doi.org/10.3390/jpm11070647)

5. **Zhou Y, Tremmel R, Schaeffeler E, Schwab M, Lauschke VM** (2022). Challenges and opportunities associated with rare-variant pharmacogenomics. *Trends in Pharmacological Sciences*, 43(10):852–865. DOI: [10.1016/j.tips.2022.07.002](https://doi.org/10.1016/j.tips.2022.07.002)

6. **Park A, Lee Y, Nam S** (2023). A performance evaluation of drug response prediction models for individual drugs. *Scientific Reports*, 13:12316. DOI: [10.1038/s41598-023-39179-2](https://doi.org/10.1038/s41598-023-39179-2)

7. **Taherdoost H, Ghofrani A** (2024). AI's role in revolutionizing personalized medicine by reshaping pharmacogenomics and drug therapy. *Intelligent Pharmacy*, 2(5):660–669. DOI: [10.1016/j.ipha.2024.08.005](https://doi.org/10.1016/j.ipha.2024.08.005)

8. **Gill PS, Yu F, Porter-Gill P, et al.** (2021). Implementing Pharmacogenomics Testing: Single Center Experience at Arkansas Children's Hospital. *Journal of Personalized Medicine*, 11(5):394. DOI: [10.3390/jpm11050394](https://doi.org/10.3390/jpm11050394)

9. **Liu W, Liu C, Liu D, et al.** (2025). GATM Identified as a Potential Drug Target for Cardiac Arrest through Multi-Omics Mendelian Randomization Integrating GWAS, eQTL, and pQTL Data. *Resuscitation*. DOI: [10.1016/s0300-9572(25)00513-1](https://doi.org/10.1016/s0300-9572(25)00513-1)

10. **Magavern EF, Daly AK, Gilchrist A, Hughes D** (2021). Pharmacogenomics spotlight commentary: From the United Kingdom to global populations. *British Journal of Clinical Pharmacology*, 87(9):3430–3434. DOI: [10.1111/bcp.14917](https://doi.org/10.1111/bcp.14917)
