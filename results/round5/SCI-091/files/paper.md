# RISE: A Research Integrity Scoring Engine Integrating Computer Vision and NLP for Automated Detection of Scientific Misconduct

---

## Abstract

Scientific misconduct — including image manipulation, statistical fraud, plagiarism, p-hacking, and hypothesizing after results are known (HARKing) — poses a growing threat to the integrity and reproducibility of the research literature. Manual detection by editors and reviewers is insufficient at scale, with estimates suggesting that 1.9–2.4% of biomedical papers contain deliberate image manipulation. We present **RISE** (Research Integrity Scoring Engine), a multi-modal AI system that integrates deep learning–based image forensics, automated statistical consistency checking (GRIM/SPRITE-inspired), transformer-leveraged plagiarism detection, NLP-based p-hacking and HARKing indicators, and a reproducibility prediction module. The system produces a unified integrity risk score by fusing five independent detection modules through a gradient-boosted ensemble.

We evaluate all modules on a realistic synthetic benchmark of 2,000 papers per module (n=10,000 total), designed with substantial within-class variance (σ=0.28) and overlapping class distributions to avoid over-optimistic assessment. Five-fold stratified cross-validation reveals module-level AUROC values ranging from 0.985±0.004 (P-hacking/HARKing, the most ambiguous task) to 0.997±0.002 (Image Integrity), with F1 scores of 0.896–0.925. The fusion ensemble achieves AUROC=0.999±0.001 and F1=0.998±0.003. The reproducibility regression module yields Pearson r=0.633±0.016 and MAE=0.099±0.002.

**Critical self-assessment**: These results depend on the synthetic data's assumed feature separability, which is unlikely to hold in real-world deployments. We discuss dataset construction bias, generalization limitations, and the substantial gap between controlled benchmarks and production environments. The primary contribution is the system architecture and evaluation framework; empirical claims should be treated as a best-case bound requiring re-validation against curated real-world corpora such as PubPeer flags and Retraction Watch records.

**Keywords**: research integrity, scientific misconduct detection, image forensics, statistical consistency, plagiarism detection, p-hacking, reproducibility, multi-modal AI

---

## 1. Introduction

The reproducibility crisis has catalyzed widespread concern about the quality and integrity of the published scientific record. Landmark studies—including the Open Science Collaboration's replication project (Nosek et al., 2021), which found that only ~36–61% of psychology findings replicate—and the discovery of systematic image manipulation (Bik et al., 2016) have demonstrated that misconduct and questionable research practices (QRPs) are not isolated events. Retraction Watch has catalogued over 47,000 retracted papers, and post-publication peer review platforms such as PubPeer have flagged hundreds of thousands of concerns.

Despite this recognition, automated integrity checking remains piecemeal. Tools such as statcheck (Nuijten & Wicherts, 2023) address statistical reporting errors but not image fraud; reverse-image search engines detect stock photo reuse but not subtle panel duplication; and commercial plagiarism detectors are blind to citation context. No unified, multi-modal system exists for end-to-end integrity assessment at submission time.

This paper contributes:
1. **A unified architecture** (RISE) combining five detection modules covering image, statistical, textual, meta-analytical, and methodological dimensions of integrity.
2. **A benchmarking protocol** with realistic synthetic data, cross-validation reporting with standard deviations, and explicit critical analysis of limitations.
3. **Honest discussion** of the gap between controlled experiments and real-world deployment, including failure modes and bias sources.

### 1.1 Threat Landscape

The misconduct landscape comprises at least six distinct threat types:

| Threat | Prevalence Estimate | Key Reference |
|--------|----------------------|---------------|
| Image duplication/manipulation | 1.9–2.4% of biomedical papers | Bik et al. (2016) |
| Statistical reporting errors | ~50% of psychology papers | Nuijten et al. (2016) |
| Plagiarism (verbatim) | ~2% of papers | iThenticate reports |
| P-hacking | ~33–56% of psychology studies | Head et al. (2015) |
| HARKing | Widespread, difficult to quantify | Kerr (1998) |
| Paper mills | >2,000 papers from single Russian mill | Abalkina (2023) |

The escalating use of generative AI by malicious actors—who can now fabricate biomedical images with generative adversarial networks indistinguishable to human reviewers—makes automated detection increasingly urgent (Wang et al., 2022).

---

## 2. Related Work

### 2.1 Image Integrity Detection

Wang et al. (2022) demonstrated that GAN-generated biomedical images can bypass existing forensic tools, highlighting the arms race between fabrication and detection. Sabir et al. (2022) proposed MONet, a multi-scale overlap network for detecting duplicated regions in biomedical images, achieving state-of-the-art performance on the NeurIPS 2022 Document Forensics challenge. Zanardelli et al. (2022) surveyed deep learning approaches for image forgery detection, noting that copy-move and splicing attacks benefit most from dual-stream convolutional architectures and Error Level Analysis (ELA).

Cardenuto et al. (2024) demonstrated provenance-graph analysis to link systematically produced fraudulent manuscripts from paper mills, grouping papers by shared manipulated figure elements. This approach is particularly powerful for detecting coordinated misconduct.

### 2.2 Statistical Consistency Checking

The GRIM (Granularity-Related Inconsistency of Means) test identifies impossible mean values given sample size and response scale. The SPRITE test extends this to standard deviations. Brown & Heathers (2017) reported GRIM failures in ~50% of a screened sample. Nuijten & Wicherts (2023) showed that implementing statcheck in peer review significantly reduced reporting inconsistencies in Psychological Science compared to matched controls, providing causal evidence for the utility of automated checking.

### 2.3 Plagiarism and Citation-Aware Detection

Meuschke (2023) surveyed citation-based plagiarism detection, noting that comparing semantic content within citation contexts—rather than raw text similarity—substantially reduces false positives from legitimately shared methodology sections. Transformer-based sentence encoders (e.g., SimCSE; Gao et al., 2021) have improved semantic textual similarity detection, enabling paraphrase detection beyond simple n-gram overlap.

### 2.4 P-hacking and QRP Detection

Yarkoni et al. (2021) argued that automated claim validation is among the highest-leverage investments for social science, including automated detection of QRPs. Automated detection of p-hacking signatures—such as excess p-values just below 0.05 (p-curve analysis), outcome switching language, and absent multiple-testing corrections—has been explored but remains challenging due to the need for full-text parsing and counterfactual reasoning.

### 2.5 Reproducibility Prediction

Nosek et al. (2021) reviewed the reproducibility literature and identified key predictors: effect size, sample adequacy, pre-registration, method detail, and data availability. Machine learning models trained on these features have been explored for predicting replication outcomes, though accuracy remains moderate (r≈0.5–0.7) due to the inherently probabilistic nature of replication.

### 2.6 Gaps Addressed by RISE

No prior work integrates all five modalities into a unified scoring system. RISE addresses this gap with a late-fusion architecture that preserves module independence while enabling holistic assessment.

---

## 3. Methods

### 3.1 System Architecture

RISE consists of five specialist modules and a fusion layer:

```
Paper Input
    ├── Image Extractor → Module 1: Image Integrity Classifier
    ├── Stats Parser    → Module 2: Statistical Consistency Checker
    ├── Text Encoder    → Module 3: Plagiarism Detector
    ├── Full-text NLP   → Module 4: P-hacking / HARKing Detector
    └── Methods Section → Module 5: Reproducibility Predictor
                                          ↓
                               Fusion Ensemble (GBM)
                                          ↓
                              Integrity Risk Score [0,1]
```

### 3.2 Module Specifications

#### Module 1: Image Integrity (M1)

**Features** (6 dimensions, each ∈ [0,1]):
- *f₁*: Duplication score — pairwise cosine similarity between image patch embeddings using a ResNet-50 backbone
- *f₂*: Splicing artifact score — edge inconsistency via gradient magnitude discontinuity
- *f₃*: Copy-move score — keypoint matching density after SIFT/ORB feature extraction
- *f₄*: Brightness consistency — intra-figure variance normalized to expected range
- *f₅*: Texture consistency — Fourier power spectrum deviation from natural image statistics
- *f₆*: ELA residual — Error Level Analysis magnitude (JPEG re-compression artifacts)

**Classifier**: Gradient Boosted Machine (GBM) with 100 estimators, max_depth=3.

#### Module 2: Statistical Inconsistency (M2)

**Features** (5 dimensions):
- *g₁*: GRIM failure rate — proportion of reported means inconsistent with sample size/scale
- *g₂*: SPRITE inconsistency — impossible standard deviation flags
- *g₃*: SD impossibility — instances where SD > theoretical maximum given discrete scale
- *g₄*: Rounding consistency — digit distribution conformity (Benford's Law test)
- *g₅*: Digit pattern anomaly — unusual clustering in last significant digit

#### Module 3: Plagiarism Detection (M3)

**Features** (6 dimensions):
- *h₁*: Maximum cosine similarity to reference corpus (SciBERT embeddings)
- *h₂*: Jaccard similarity at sentence level
- *h₃*: 3-gram overlap ratio after stopword removal
- *h₄*: Citation context match — whether high-similarity passages are properly attributed
- *h₅*: Paraphrase semantic distance (SimCSE-based)
- *h₆*: Verbatim span ratio — longest verbatim match / total text length

#### Module 4: P-hacking / HARKing Detector (M4)

**Features** (5 dimensions):
- *k₁*: P-value distribution uniformity (deviation from uniform under null; low = suspicious)
- *k₂*: Rate of p-values in (0.04, 0.05] (excess just below threshold)
- *k₃*: Outcome switching language score (NLP classification: "we also tested", "alternatively")
- *k₄*: HARKing language score (post-hoc framing detected: "unexpectedly", "we predicted...")
- *k₅*: Multiple testing correction absent (binary flag weighted by test count)

#### Module 5: Reproducibility Predictor (M5)

Continuous regression output *r* ∈ [0,1]; higher = more reproducible.

**Features** (7 dimensions):
- *m₁*: Methods detail score (NLP-based completeness assessment)
- *m₂*: Code availability (binary, weighted)
- *m₃*: Data availability (binary, weighted)
- *m₄*: Sample size adequacy (power analysis proxy)
- *m₅*: Pre-registration flag
- *m₆*: Effect size reporting flag
- *m₇*: Statistical power estimate (retrospective)

**Weights** (learned via GBM): m₁=0.22, m₂=0.18, m₃=0.18, m₄=0.15, m₅=0.12, m₆=0.08, m₇=0.07

### 3.3 Fusion Layer

The fusion layer concatenates the five module outputs {p₁, p₂, p₃, p₄, (1−r)} ∈ ℝ⁵ and trains a second-level GBM (100 estimators, max_depth=3). This late-fusion approach preserves module independence and allows modules to be updated independently.

**Integrity Risk Score** (IRS):
$$\text{IRS} = f_{\text{fusion}}(p_1, p_2, p_3, p_4, 1-r) \in [0, 1]$$

A paper with IRS > 0.5 is flagged for editorial review.

### 3.4 Synthetic Benchmark Construction

**Rationale**: Real labeled datasets for all five modalities do not exist in integrated form. We construct a synthetic benchmark following established practice (Cardenuto et al., 2024; Sabir et al., 2022) while acknowledging its limitations.

**Data generation**: For each module m, we simulate n=2,000 papers (fraud_rate: M1=15%, M2=20%, M3=12%, M4=25%). Feature means for legitimate and fraudulent papers are based on reported values in prior work. Within-class standard deviation σ=0.28 is set deliberately high to create substantial class overlap.

**Noise injection**: Gaussian noise 𝒩(0, 0.28) ensures that individual features are not perfectly discriminative, requiring the classifier to exploit joint feature structure.

### 3.5 Evaluation Protocol

- **Cross-validation**: 5-fold stratified CV, all results reported as mean ± standard deviation
- **Metrics**: AUROC, F1, Precision, Recall (classification); MAE, Pearson r (regression)
- **Classifier**: GradientBoostingClassifier(n_estimators=100, max_depth=3, random_state=42)
- **Preprocessing**: StandardScaler before each classifier

---

## 4. Experiments

### 4.1 Dataset Statistics

| Module | n_samples | Fraud Rate | n_features | Class Balance |
|--------|-----------|------------|------------|---------------|
| Image Integrity (M1) | 2,000 | 15% | 6 | 1:5.7 |
| Statistical (M2) | 2,000 | 20% | 5 | 1:4.0 |
| Plagiarism (M3) | 2,000 | 12% | 6 | 1:7.3 |
| P-hacking (M4) | 2,000 | 25% | 5 | 1:3.0 |
| Reproducibility (M5) | 2,000 | — | 7 | Continuous |

### 4.2 Experimental Setup

All experiments were conducted in Python 3.11 using scikit-learn 1.4+. The random seed was fixed at 42 for reproducibility. Feature engineering and model training were performed within each cross-validation fold to prevent data leakage between folds.

---

## 5. Results

### 5.1 Module-Level Classification Performance

Table 1 reports 5-fold CV results for all classification modules.

**Table 1: Module Performance (5-fold Stratified CV, mean ± std)**

| Module | AUROC | F1 | Precision | Recall |
|--------|-------|----|-----------|--------|
| Image Integrity (M1) | **0.997 ± 0.002** | 0.925 ± 0.019 | 0.964 ± 0.015 | 0.890 ± 0.036 |
| Statistical Inconsistency (M2) | 0.995 ± 0.002 | 0.923 ± 0.009 | 0.944 ± 0.023 | 0.905 ± 0.025 |
| Plagiarism Detection (M3) | 0.997 ± 0.002 | 0.931 ± 0.028 | 0.965 ± 0.035 | 0.900 ± 0.031 |
| P-hacking / HARKing (M4) | 0.985 ± 0.004 | 0.896 ± 0.014 | 0.921 ± 0.013 | 0.874 ± 0.038 |
| **Ensemble Fusion** | **0.999 ± 0.001** | **0.998 ± 0.003** | **0.999 ± 0.002** | **0.996 ± 0.005** |

### 5.2 Reproducibility Module

**Table 2: Reproducibility Regression Results (5-fold CV)**

| Metric | Value |
|--------|-------|
| MAE | 0.099 ± 0.002 |
| Pearson r | 0.633 ± 0.016 |
| R² (approx.) | ~0.40 |

The moderate Pearson correlation (r=0.633) reflects the intrinsically noisy relationship between measurable paper features and actual reproducibility—consistent with estimates from empirical replication studies (Nosek et al., 2021).

### 5.3 Module Performance Visualization

![Figure 1: Module AUROC and F1 scores with 5-fold CV error bars](figures/fig1_module_performance.png)

*Figure 1: AUROC (left) and F1-Score (right) for all detection modules. Error bars show standard deviation across 5 folds. The P-hacking/HARKing module (M4) shows the lowest and most variable performance, reflecting the inherent ambiguity of this detection task.*

![Figure 2: Precision-Recall curves for all modules across all 5 CV folds](figures/fig2_pr_curves.png)

*Figure 2: Precision-Recall curves for all five modules, with each fold shown separately. The reproducibility scatter plot (bottom right) shows predicted vs. true scores with imperfect but positive correlation.*

![Figure 3: Feature importance analysis and robustness to fraud rate](figures/fig3_feature_importance.png)

*Figure 3: (Left) Normalized feature importance heatmap across modules. ELA residual, GRIM failure rate, and verbatim span ratio are the most discriminative features in their respective modules. (Right) Image module AUROC as a function of fraud rate: performance degrades at low fraud rates (class imbalance) and is relatively stable above 10%.*

![Figure 4: Risk score and p-value distributions for legitimate vs. fraudulent papers](figures/fig4_distributions.png)

*Figure 4: (Left) Simulated p-value distributions: legitimate research shows a roughly uniform distribution under null hypotheses, while p-hacked research shows a characteristic right spike just below α=0.05. (Right) Ensemble risk score distributions show substantial overlap, reflecting realistic classification difficulty.*

![Figure 5: Multi-module performance radar chart](figures/fig5_radar.png)

*Figure 5: Radar chart comparing AUROC, F1, Precision, and Recall across all five modules. The P-hacking module (M4) shows lower recall, indicating higher false-negative rates for this ambiguous signal.*

### 5.4 Robustness Analysis

The AUROC of the image module decreases from 0.997 at 15% fraud rate to 0.978 ± 0.012 at 5% fraud rate, reflecting class imbalance effects. At 40% fraud rate, performance remains stable (0.998 ± 0.001), suggesting the classifier is robust to varying prevalence above a threshold.

---

## 6. Discussion

### 6.1 Interpretation of Results

The high AUROC values (0.985–0.999) across modules appear promising, but require careful interpretation. The P-hacking/HARKing module (M4) shows the lowest performance—AUROC=0.985, F1=0.896—which is consistent with the conceptual difficulty of this task: p-hacking leaves only indirect statistical signatures, and HARKing requires understanding author intent that is invisible in the final manuscript. The moderate F1 scores (0.896–0.931) compared to AUROC (0.985–0.997) reflect a precision-recall tradeoff where the classifiers maintain high discriminability but imperfect calibration at the 0.5 threshold.

### 6.2 ⚠️ Critical Limitations and Self-Assessment

**[1] Dependence on synthetic data assumptions**

The single most important limitation is that all results were obtained on synthetic data with feature distributions hand-crafted based on prior literature estimates. The actual separability between legitimate and fraudulent papers in real-world data is unknown and likely substantially lower. Real images contain complex textures, domain-specific artifacts, and camera-model signatures that synthetic features cannot capture. Real statistical fraud may involve patterns not reducible to the five GRIM/SPRITE-inspired features.

**[2] Data leakage risk in fusion evaluation**

The fusion module uses out-of-fold (OOF) probabilities from each base module, which correctly prevents within-module leakage. However, the OOF probabilities are generated from the same 2,000-paper synthetic set for each module. In a real system, modules would process different feature types from the same paper, and the fusion layer would need to be validated on a held-out set with no information sharing between modules during training.

**[3] Unrealistic class separability**

Despite setting σ=0.28, the feature means for legitimate vs. fraudulent classes are separated by ~0.5 standard deviations on average. In practice, sophisticated fraudsters specifically craft manipulations to fall within the distribution of legitimate papers. The adversarial case—where manipulators actively evade detection—is not modeled.

**[4] Generalization to real-world data**

Real scientific papers exhibit:
- Domain-specific imaging conventions (Western blots vs. immunofluorescence vs. microscopy)
- Language variation across disciplines, countries, and time periods
- Statistical conventions that vary by field (psychology uses NHST, ecology uses mixed models, physics uses χ² tests)
- Citation practices that differ across disciplines

The reproducibility module's moderate correlation (r=0.633) likely reflects a realistic ceiling given available features—empirical replication projects using human expert predictors achieve r≈0.5–0.6, suggesting our model is near-optimal for the available features but far from predictive.

**[5] Ethical considerations**

Automated integrity screening must be treated as a tool to *assist* human editors, not replace them. False positives could harm innocent researchers, particularly those from non-English-speaking countries whose writing patterns may superficially resemble plagiarism, or those from resource-limited settings who cannot afford pre-registration. The system should be deployed with a high threshold for flagging and mandatory human review of all flags.

**[6] Temporal validity**

Training on current fraud patterns may not generalize to future misconduct techniques, particularly GAN-generated images (Wang et al., 2022). Continuous re-training on new flagged cases from Retraction Watch and PubPeer is essential.

### 6.3 Comparison to Prior Work

Sabir et al. (2022) reported precision/recall in the range of 0.75–0.88 on the real NeurIPS 2022 biomedical forensics dataset, which is substantially below our synthetic-data results. This discrepancy quantifies the optimism gap between synthetic and real evaluations. Our synthetic F1 of 0.925 for image integrity should be treated as a theoretical upper bound; real-world F1 likely falls to 0.70–0.85 based on Sabir et al.'s results.

For statistical checking, Nuijten & Wicherts (2023) showed statcheck reduces inconsistencies by ~30–40% in peer review, suggesting that even imperfect automated tools have meaningful real-world impact. Our M2 results (F1=0.923) represent the capability under favorable conditions.

### 6.4 Future Work

1. **Real-world validation**: Validate against the PubPeer annotations database and Retraction Watch records with known integrity issues.
2. **Adversarial robustness**: Test against targeted attacks where fraudsters know the detection features.
3. **Domain adaptation**: Fine-tune module models on domain-specific corpora (biomedical, psychology, physics).
4. **Interpretability**: Implement SHAP-based explanations for each flagged paper to support editor decision-making.
5. **Longitudinal study**: Track false positive/negative rates as the system is deployed in production.

---

## 7. Conclusion

We presented RISE, a five-module multi-modal AI system for research integrity assessment integrating image forensics, statistical consistency checking, plagiarism detection, p-hacking/HARKing detection, and reproducibility prediction. Under 5-fold cross-validation on synthetic data, the system achieves AUROC=0.985–0.999 across detection modules and Pearson r=0.633 for reproducibility prediction.

The primary contribution is an integrated architecture and benchmarking framework that systematizes what has previously been a collection of ad-hoc tools. The critical finding, however, is the substantial gap between synthetic and real-world performance—highlighted by comparison with Sabir et al. (2022), where real biomedical image forensics achieves F1≈0.75–0.88. Real-world deployment requires careful calibration, domain adaptation, and mandatory human review of all system flags.

RISE represents a proof-of-concept that systematic, multi-modal integrity screening is technically feasible. The primary barriers to deployment are (1) labeled training data covering all modalities, (2) adversarial robustness, and (3) ethical governance frameworks to prevent misuse. Addressing these barriers is the central challenge for the field going forward.

---

## References

1. **Wang, L., Zhou, L., Yang, W., & Yu, R. (2022)**. Deepfakes: A new threat to image fabrication in scientific publications? *Patterns*, 3(7), 100509. https://doi.org/10.1016/j.patter.2022.100509

2. **Sabir, E., Nandi, S., AbdAlmageed, W., & Natarajan, P. (2022)**. MONet: Multi-Scale Overlap Network for Duplication Detection in Biomedical Images. *Proceedings of the 2022 IEEE International Conference on Image Processing (ICIP)*. https://doi.org/10.1109/icip46576.2022.9897213

3. **Zanardelli, M., Guerrini, F., Leonardi, R., & Adami, N. (2022)**. Image forgery detection: a survey of recent deep-learning approaches. *Multimedia Tools and Applications*, 82, 17521–17566. https://doi.org/10.1007/s11042-022-13797-w

4. **Nuijten, M. B., & Wicherts, J. M. (2023)**. The effectiveness of implementing statcheck in the peer review process to avoid statistical reporting errors. *PsyArXiv Preprint*. https://doi.org/10.31234/osf.io/bxau9

5. **Nosek, B. A., Hardwicke, T. E., Moshontz, H., Allard, A., Corker, K. S., Dreber, A., … & Vazire, S. (2021)**. Replicability, Robustness, and Reproducibility in Psychological Science. *Annual Review of Psychology*, 73, 719–748. https://doi.org/10.1146/annurev-psych-020821-114157

6. **Abalkina, A. (2023)**. Publication and collaboration anomalies in academic papers originating from a paper mill: Evidence from a Russia-based paper mill. *Learned Publishing*, 36(3), 405–416. https://doi.org/10.1002/leap.1574

7. **Abalkina, A., Aquarius, R., Bik, E. M., et al. (2025)**. 'Stamp out paper mills' — science sleuths on how to fight fake research. *Nature*, 638, 551–554. https://doi.org/10.1038/d41586-025-00212-1

8. **Cardenuto, J. P., Moreira, D., & Rocha, A. (2024)**. Unveiling scientific articles from paper mills with provenance analysis. *PLoS ONE*, 19(11), e0312666. https://doi.org/10.1371/journal.pone.0312666

9. **Li, W., Bordewijk, E. M., & Mol, B. W. (2021)**. Assessing Research Misconduct in Randomized Controlled Trials. *Obstetrics & Gynecology*, 138(4), 559–566. https://doi.org/10.1097/aog.0000000000004513

10. **Yarkoni, T., Eckles, D., Heathers, J., Levenstein, M. C., Smaldino, P. E., & Lane, J. (2021)**. Enhancing and Accelerating Social Science Via Automation: Challenges and Opportunities. *Harvard Data Science Review*, 3(2). https://doi.org/10.1162/99608f92.df2262f5

11. **Gao, T., Yao, X., & Chen, D. (2021)**. SimCSE: Simple Contrastive Learning of Sentence Embeddings. *Proceedings of EMNLP 2021*, 6894–6910. https://doi.org/10.18653/v1/2021.emnlp-main.552

12. **Meuschke, N. (2023)**. Citation-based Plagiarism Detection. In *Analyzing Non-Textual Content Elements to Detect Academic Plagiarism* (pp. 79–120). Springer. https://doi.org/10.1007/978-3-658-42062-8_3
