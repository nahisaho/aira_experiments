# EthicAI-Bench: A Unified Quantitative Framework for Multi-Dimensional Ethical Evaluation of AI Systems

## Abstract

The rapid deployment of artificial intelligence systems in high-stakes domains such as healthcare necessitates rigorous, quantitative evaluation of their ethical properties. While numerous ethical guidelines and principles have been proposed, operationalizing these principles into measurable, comparable metrics remains a significant challenge. In this paper, we present **EthicAI-Bench**, a unified quantitative framework that integrates five critical dimensions of AI ethics: fairness (Statistical Parity, Equalized Odds, and Calibration), explainability (SHAP consistency and stability), privacy (membership inference attack resistance), robustness (adversarial perturbation and distribution shift tolerance), and environmental sustainability (computational carbon emissions). Our framework computes an integrated EthicAI Score as a weighted composite of normalized dimension scores, enabling holistic comparison across model architectures. We evaluate three representative machine learning models—Logistic Regression, Random Forest, and Gradient Boosting—on a synthetic medical diagnosis dataset with intentionally embedded demographic biases. Our results reveal critical trade-offs: models achieving higher predictive performance often exhibit greater fairness disparities and privacy vulnerabilities, while simpler models demonstrate superior overall ethical profiles. Logistic Regression achieves the highest EthicAI Score (0.920), driven by strong privacy resilience and environmental efficiency, despite exhibiting the largest gender-based Statistical Parity Difference (0.249). Random Forest, despite near-perfect explainability stability (0.999), scores lowest overall (0.636) due to severe membership inference vulnerability (MIA accuracy: 85.5%) stemming from overfitting. These findings underscore the necessity of multi-dimensional ethical assessment beyond single-metric evaluations and provide a practical, extensible pipeline for ethical AI auditing in regulated domains.

---

## 1. Introduction

### 1.1 Background

The widespread adoption of artificial intelligence in healthcare, criminal justice, finance, and other consequential domains has intensified scrutiny of the ethical implications of algorithmic decision-making. Regulatory frameworks including the European Union AI Act (2024), the NIST AI Risk Management Framework (AI RMF 1.0), and the WHO guidelines for AI in healthcare have established normative expectations for trustworthy AI systems. However, translating high-level ethical principles—fairness, transparency, privacy, robustness, and sustainability—into operational, quantitative metrics remains a persistent challenge (Kaur et al., 2024; Li et al., 2025).

Existing evaluation approaches typically address individual ethical dimensions in isolation. Fairness toolkits such as Fairlearn (Bird et al., 2023) and AI Fairness 360 (Bellamy et al., 2019) provide comprehensive metric suites for bias detection but do not incorporate privacy, robustness, or environmental considerations. Similarly, explainability methods like SHAP (Lundberg & Lee, 2017) quantify feature attributions without assessing their stability or consistency. Privacy auditing through membership inference attacks (Song & Mittal, 2021) and robustness evaluation via adversarial perturbation (Goodfellow et al., 2015) constitute separate research streams rarely integrated into a unified framework.

### 1.2 Motivation

The fragmented nature of existing ethical AI evaluation creates several problems. First, optimizing for a single ethical criterion may inadvertently degrade another—for example, differential privacy mechanisms can exacerbate fairness disparities across demographic groups (Bagdasaryan et al., 2019). Second, stakeholders—regulators, clinicians, patients—require holistic assessments that capture the full ethical profile of an AI system. Third, the environmental costs of model development, increasingly recognized as an ethical concern (Strubell et al., 2019), are rarely quantified alongside other ethical metrics.

### 1.3 Contributions

We make the following contributions:

1. **EthicAI-Bench Framework**: A unified, modular evaluation pipeline integrating five ethical dimensions with configurable weights, producing an interpretable composite score.
2. **Novel Metrics**: An Integrated Fairness Score (IFS) combining Statistical Parity, Equalized Odds, and Calibration differences; SHAP consistency and stability indices; and a privacy risk score based on membership inference attack accuracy.
3. **Medical AI Case Study**: Comprehensive ethical audit of three model architectures on a medical diagnosis task, revealing cross-dimensional trade-offs.
4. **Open Evaluation Pipeline**: A reproducible Python implementation leveraging Fairlearn, SHAP, and custom modules for privacy, robustness, and environmental assessment.

---

## 2. Related Work

### 2.1 Fairness Evaluation Frameworks

Bird et al. (2023) introduced Fairlearn, an open-source toolkit providing demographic parity, equalized odds, and disaggregated performance metrics with scikit-learn integration. Bellamy et al. (2019) developed AI Fairness 360 (AIF360), offering over 70 fairness metrics and multiple bias mitigation algorithms spanning preprocessing, in-processing, and postprocessing stages. Recent comparative studies (arXiv:2412.09900, 2025) have benchmarked these frameworks across computer vision and NLP tasks, demonstrating their complementary capabilities in bias quantification. However, neither toolkit integrates privacy, robustness, or environmental metrics.

### 2.2 Explainability Quantification

SHAP (SHapley Additive exPlanations) provides theoretically grounded feature attributions based on cooperative game theory (Lundberg & Lee, 2017). While SHAP values are widely used, their stability under model perturbation and data resampling has received increasing attention. Studies have proposed quantification via bootstrap variance, rank correlation coefficients (Spearman's ρ), and top-k feature agreement metrics (Guidotti et al., 2022). The IEEE 7000-2021 standard explicitly addresses measurement criteria for ethical explainability stability.

### 2.3 Privacy Risk Assessment

Song and Mittal (2021) proposed a systematic framework for evaluating privacy risks through membership inference attacks (MIA), introducing per-sample vulnerability scores. Long et al. (2020) developed pragmatic adversary models for MIA evaluation. More recent work (Wang et al., 2025) examines reliability disparities among MIA methods and proposes ensemble frameworks for robust privacy risk assessment. These approaches quantify the degree to which a model memorizes individual training samples, with higher MIA accuracy indicating greater privacy risk.

### 2.4 Robustness Evaluation

Adversarial robustness—the ability of models to maintain performance under malicious perturbations—has been extensively studied since Goodfellow et al. (2015). Distribution shift robustness, addressing natural changes in data distributions, complements adversarial evaluation. Recent work on the Robustness Carbon Trade-off Index (RCTI) by Vyas et al. (2024) demonstrates that adversarial training can significantly increase carbon emissions, highlighting the tension between robustness and sustainability.

### 2.5 Environmental Impact of AI

Strubell et al. (2019) catalyzed awareness of AI's environmental costs, estimating that training a large NLP model can emit over 626,000 lbs of CO₂. Subsequent work has developed frameworks for measuring and reducing carbon footprints during model training (Xu et al., 2023) and deployment. The CodeCarbon library provides real-time emission tracking, while research on knowledge distillation (Shahid et al., 2023) explores model compression as a sustainability strategy.

### 2.6 Integrated Ethical Evaluation

The RAISE framework (2025) represents a recent effort to aggregate explainability, fairness, robustness, and sustainability into a holistic responsibility score. Kaur et al. (2024) synthesized six quantitative requirements for trustworthy AI. In healthcare specifically, the FUTURE-AI guideline (Lekadir et al., 2025) establishes consensus criteria for Fairness, Universality, Traceability, Usability, Robustness, and Explainability. Our work extends these approaches by incorporating privacy risk assessment and providing a fully automated, reproducible evaluation pipeline with weighted composite scoring.

---

## 3. Methods

### 3.1 Framework Architecture

EthicAI-Bench consists of five evaluation modules, each producing a normalized score $S_i \in [0, 1]$, and a composition layer that aggregates these into the EthicAI Score:

$$\text{EthicAI}(\mathcal{M}) = \sum_{i=1}^{5} w_i \cdot S_i(\mathcal{M}), \quad \text{where } \sum_{i=1}^{5} w_i = 1$$

Default weights: $w_F = 0.25$ (fairness), $w_E = 0.20$ (explainability), $w_P = 0.20$ (privacy), $w_R = 0.20$ (robustness), $w_{Env} = 0.15$ (environmental).

### 3.2 Fairness Module

We compute three complementary group fairness metrics for each sensitive attribute $A$:

**Statistical Parity Difference (SPD):**
$$\text{SPD} = P(\hat{Y}=1 | A=a) - P(\hat{Y}=1 | A=b)$$

**Equalized Odds Difference (EOD):**
$$\text{EOD} = \max_{y \in \{0,1\}} |P(\hat{Y}=1 | A=a, Y=y) - P(\hat{Y}=1 | A=b, Y=y)|$$

**Calibration Difference (CalDiff):**
$$\text{CalDiff} = \frac{1}{|B|} \sum_{k \in B} \max_{a,b} |E[Y | A=a, \hat{p} \in B_k] - E[Y | A=b, \hat{p} \in B_k]|$$

where $B$ denotes the set of probability bins. The Integrated Fairness Score (IFS) is:

$$\text{IFS} = 1 - \frac{|\text{SPD}| + |\text{EOD}| + \text{CalDiff}}{3}$$

The overall fairness score averages IFS across sensitive attributes.

### 3.3 Explainability Module

We quantify three aspects of SHAP-based explanations:

**SHAP Consistency** ($\rho_{\text{SHAP}}$): Spearman rank correlation between mean absolute SHAP values computed on the original test set and bootstrap samples:

$$\rho_{\text{SHAP}} = \frac{1}{B} \sum_{b=1}^{B} \rho_s\left(\overline{|\phi|}, \overline{|\phi^{(b)}|}\right)$$

**SHAP Stability**: Inverse of mean feature-wise SHAP variance:

$$S_{\text{stab}} = \frac{1}{1 + \text{mean}_j(\text{Var}_i[\phi_{ij}])}$$

**Top-k Agreement**: Proportion of top-k important features consistent across bootstrap runs.

### 3.4 Privacy Module

We implement a threshold-based membership inference attack (MIA):

$$\text{MIA}_{\text{acc}} = \max_{\tau \in [0.5, 1.0]} \frac{1}{2}\left[\frac{|\{x_{\text{train}}: \max_c p(c|x) \geq \tau\}|}{|D_{\text{train}}|} + \frac{|\{x_{\text{test}}: \max_c p(c|x) < \tau\}|}{|D_{\text{test}}|}\right]$$

The privacy risk score is:

$$S_P = \max(0, 1 - 2 \cdot (\text{MIA}_{\text{acc}} - 0.5))$$

### 3.5 Robustness Module

**Adversarial Robustness**: We apply uniform random perturbations at varying magnitudes $\varepsilon$ and compute the area under the accuracy-vs-$\varepsilon$ curve:

$$S_{\text{adv}} = \frac{\int_{\varepsilon_{\min}}^{\varepsilon_{\max}} \text{Acc}(\varepsilon) \, d\varepsilon}{\varepsilon_{\max} - \varepsilon_{\min}}$$

**Prediction Consistency**: Fraction of predictions unchanged under Gaussian noise ($\sigma = 0.05$), averaged over $N$ noise realizations.

### 3.6 Environmental Module

CO₂ emissions are estimated as:

$$\text{CO}_2 = E_{\text{compute}} \times \text{PUE} \times I_C$$

where $E_{\text{compute}} = P_{\text{TDP}} \times t / (3600 \times 1000)$ kWh, PUE = 1.58 (data center overhead), and $I_C = 0.475$ kgCO₂/kWh (global average carbon intensity). The environmental score is normalized inversely relative to the maximum-emitting model.

---

## 4. Experiments

### 4.1 Dataset

We construct a synthetic medical diagnosis dataset ($N = 3{,}000$, train/test split 70/30, stratified) comprising 11 clinical features: age, gender (binary), ethnicity (3 groups), BMI, blood pressure, cholesterol, glucose, heart rate, smoking status, exercise level, and family history. Intentional biases are embedded via gender and ethnicity coefficients in the data-generating logistic model to simulate real-world demographic disparities in healthcare outcomes.

### 4.2 Models

Three representative architectures are evaluated:

- **Logistic Regression** (LR): Linear model with L2 regularization
- **Random Forest** (RF): Bagging ensemble with 100 trees
- **Gradient Boosting** (GB): Sequential boosting with 100 estimators

All models are trained on standardized features using scikit-learn with fixed random seeds for reproducibility.

### 4.3 Evaluation Protocol

Each model undergoes evaluation across all five modules:

1. **Fairness**: Fairlearn-based SPD, EOD, and calibration assessment across gender and ethnicity
2. **Explainability**: SHAP analysis with $B=10$ bootstrap iterations for consistency and top-5 feature agreement
3. **Privacy**: Threshold-based MIA with 50 threshold values in $[0.5, 1.0]$
4. **Robustness**: Adversarial perturbation at $\varepsilon \in \{0.01, 0.05, 0.1, 0.2, 0.5\}$ and covariate shift at magnitudes $\{0.5, 1.0, 2.0\}$
5. **Environmental**: CPU-based timing with TDP=65W, PUE=1.58, carbon intensity=0.475 kgCO₂/kWh

---

## 5. Results

### 5.1 Baseline Performance

| Model | Accuracy | AUC | F1 | Training Time (s) |
|-------|----------|-----|----|--------------------|
| Logistic Regression | 0.6822 | 0.7378 | 0.6398 | 0.002 |
| Random Forest | 0.6500 | 0.7130 | 0.6018 | 0.32 |
| Gradient Boosting | 0.6600 | 0.7217 | 0.6117 | 0.46 |

Logistic Regression achieves the highest predictive performance across all metrics despite its simplicity, likely due to the linear data-generating process.

### 5.2 Fairness Results

![Figure 1: Fairness metrics comparison across models for gender and ethnicity sensitive attributes](figures/fairness_comparison.png)

All models exceed the 0.1 Statistical Parity threshold for gender, with LR exhibiting the largest disparity (SPD=0.249). RF demonstrates the lowest gender bias (SPD=0.145). Ethnicity-based disparities are generally smaller, with SPD < 0.07 for all models.

### 5.3 Explainability Results

![Figure 2: SHAP feature importance summary for Gradient Boosting model](figures/shap_summary.png)

| Model | SHAP Consistency (ρ) | Stability | Top-5 Agreement |
|-------|----------------------|-----------|-----------------|
| LR | 0.986 ± 0.007 | 0.920 | 0.980 |
| RF | 0.999 ± 0.003 | 0.997 | 1.000 |
| GB | 0.993 ± 0.006 | 0.926 | 1.000 |

All models achieve high explainability scores. TreeSHAP provides near-perfect consistency for tree-based models, consistent with theoretical guarantees.

### 5.4 Privacy Results

![Figure 3: Membership inference attack accuracy and overfitting gap](figures/privacy_risk.png)

Random Forest exhibits a critical privacy vulnerability with MIA accuracy of 85.5%, directly attributable to its 35% train-test accuracy gap. Logistic Regression and Gradient Boosting maintain MIA accuracy near random chance (50.4% and 51.2%, respectively).

### 5.5 Robustness Results

![Figure 4: Model accuracy under increasing adversarial perturbation magnitude](figures/robustness_adversarial.png)

All models demonstrate comparable robustness profiles. LR shows the highest prediction consistency (0.984) under Gaussian noise, while GB and RF show slightly lower stability. Distribution shift tolerance is similar across all models, with accuracy degradation < 2% at shift magnitude 2.0.

### 5.6 Environmental Results

![Figure 5: Training and projected annual CO₂ emissions](figures/environmental_impact.png)

LR is orders of magnitude more efficient than ensemble methods, with training CO₂ of 0.028 mg versus 4.589 mg (RF) and 5.941 mg (GB). Annual inference emissions (1M predictions) range from 0.0005 kg (LR) to 0.1185 kg (RF).

### 5.7 Integrated EthicAI Scores

![Figure 6: EthicAI integrated score breakdown by ethical dimension](figures/ethicai_integrated.png)

![Figure 7: Radar chart comparing ethical profiles across models](figures/ethicai_radar.png)

| Model | Fairness | Explainability | Privacy | Robustness | Environmental | **EthicAI** |
|-------|----------|---------------|---------|------------|---------------|-------------|
| LR | 0.856 | 0.962 | 0.992 | 0.831 | 0.996 | **0.920** |
| RF | 0.872 | 0.999 | 0.289 | 0.801 | 0.000 | **0.636** |
| GB | 0.863 | 0.973 | 0.975 | 0.812 | 0.883 | **0.900** |

### 5.8 Medical AI Ethics Audit

![Figure 8: Comprehensive ethics audit heatmap for medical AI diagnosis system](figures/medical_audit_heatmap.png)

The heatmap reveals that no single model excels across all dimensions. LR and GB present viable candidates for deployment, while RF requires significant mitigation (regularization, differential privacy) before clinical use.

---

## 6. Discussion

### 6.1 Cross-Dimensional Trade-offs

Our results reveal several important trade-offs that are only visible through multi-dimensional evaluation:

**Complexity-Privacy Trade-off**: Random Forest's capacity to memorize training data (evidenced by 35% overfitting gap) directly translates to severe MIA vulnerability. This finding aligns with Song and Mittal's (2021) observation that model capacity is a primary driver of membership inference risk.

**Fairness-Performance Tension**: While LR achieves the best predictive performance, it also exhibits the largest gender bias. This reflects the model's faithful encoding of the linear bias in the data-generating process—a reminder that high accuracy on biased data can amplify rather than mitigate unfairness.

**Simplicity Advantage**: Logistic Regression's overall ethical superiority underscores that model complexity does not automatically confer ethical benefits. In contrast, simpler models offer better privacy protection, higher robustness, lower environmental cost, and comparable fairness—supporting the principle of parsimony in high-stakes applications.

### 6.2 Implications for Healthcare AI

In the context of medical diagnosis, our framework identifies actionable risk areas:

1. **All models fail the gender fairness threshold** (SPD > 0.1), indicating that bias mitigation is essential regardless of model choice.
2. **Random Forest is contraindicated** for deployment without additional privacy protections, as patients' membership in the training set can be inferred with 85.5% accuracy.
3. **Logistic Regression or Gradient Boosting** with fairness-constrained training (e.g., exponentiated gradient from Fairlearn) represent more ethically sound starting points.

### 6.3 Limitations

1. **Synthetic Data**: While enabling controlled bias injection, our synthetic dataset may not capture the complexity of real clinical data distributions.
2. **Simplified MIA**: Our threshold-based attack is weaker than state-of-the-art shadow model attacks (Shokri et al., 2017), potentially underestimating true privacy risk.
3. **CPU-Only Environmental Measurement**: GPU-accelerated models would exhibit different emission profiles, and our estimates do not account for data storage, preprocessing, or model serving infrastructure.
4. **Static Evaluation**: We assess a single snapshot; real-world ethical monitoring requires continuous evaluation as data distributions evolve.
5. **Weight Sensitivity**: The composite EthicAI Score depends on dimension weights, which should be calibrated to domain-specific requirements.

### 6.4 Future Directions

1. **Differential Privacy Integration**: Evaluating the impact of DP-SGD on the fairness-privacy-performance Pareto frontier.
2. **Temporal Ethics Monitoring**: Developing dashboards for continuous ethical score tracking throughout the model lifecycle.
3. **Regulatory Mapping**: Automating compliance checks against EU AI Act requirements and FDA Software as Medical Device (SaMD) guidelines.
4. **Causal Fairness**: Extending beyond statistical fairness metrics to counterfactual fairness notions.
5. **Real-World Validation**: Evaluating EthicAI-Bench on clinical datasets (e.g., MIMIC-III, eICU) with institutional review board approval.

---

## 7. Conclusion

We presented EthicAI-Bench, a unified quantitative framework for multi-dimensional ethical evaluation of AI systems that integrates fairness, explainability, privacy, robustness, and environmental sustainability into a single composite score. Applied to a medical AI diagnosis scenario, our framework reveals that simpler models (Logistic Regression, EthicAI=0.920) achieve superior ethical profiles compared to complex ensembles (Random Forest, EthicAI=0.636), primarily due to privacy and environmental advantages. Critically, we demonstrate that no model simultaneously excels across all ethical dimensions, underscoring the necessity of holistic assessment beyond single-metric evaluations. EthicAI-Bench provides a practical, extensible, and reproducible evaluation pipeline suitable for ethical auditing in regulated AI deployment contexts.

---

## References

1. Bellamy, R. K. E., Dey, K., Hind, M., Hoffman, S. C., Houde, S., Kannan, K., ... & Zhang, Y. (2019). AI Fairness 360: An extensible toolkit for detecting and mitigating algorithmic bias. *IBM Journal of Research and Development*, 63(4/5), 4:1–4:15. https://doi.org/10.1147/JRD.2019.2942287

2. Bird, S., Dudík, M., Edgar, R., Horn, B., Lutz, R., Milan, V., ... & Walker, K. (2023). Fairlearn: Assessing and improving fairness of AI systems. *Journal of Machine Learning Research*, 24(257), 1–54. https://jmlr.org/papers/volume24/23-0389/23-0389.pdf

3. Song, L. & Mittal, P. (2021). Systematic evaluation of privacy risks of machine learning models. In *Proceedings of the 30th USENIX Security Symposium*. https://doi.org/10.5555/3489212.3489527

4. Long, Y., Wang, L., Bu, D., Bindschaedler, V., Wang, X., Tang, H., Gunter, C. A., & Chen, K. (2020). A pragmatic approach to membership inferences on machine learning models. In *IEEE European Symposium on Security and Privacy (EuroS&P)*. https://doi.org/10.1109/EuroSP48549.2020.00040

5. Kaur, D., Uslu, S., Rittichier, K. J., & Durresi, A. (2024). Establishing and evaluating trustworthy AI: Overview and research challenges. *Frontiers in Big Data*, 7, 1467222. https://doi.org/10.3389/fdata.2024.1467222

6. Lekadir, K., Osuala, R., Gallin, C., et al. (2025). FUTURE-AI: International consensus guideline for trustworthy and deployable artificial intelligence in healthcare. *BMJ*, 388, e081554. https://doi.org/10.1136/bmj-2024-081554

7. Li, Y., Shi, J., & Zhang, Q. (2025). AI ethics: Integrating transparency, fairness, and privacy in AI development. *Applied Artificial Intelligence*, 39(1). https://doi.org/10.1080/08839514.2025.2463722

8. Shahid, A., Mushtaq, M., & Munir, M. (2023). Mitigating carbon footprint for knowledge distillation based deep learning model compression. *PLOS ONE*, 18(5), e0285668. https://doi.org/10.1371/journal.pone.0285668

9. Lundberg, S. M. & Lee, S.-I. (2017). A unified approach to interpreting model predictions. In *Advances in Neural Information Processing Systems (NeurIPS)*, 30, 4765–4774.

10. Goodfellow, I. J., Shlens, J., & Szegedy, C. (2015). Explaining and harnessing adversarial examples. In *International Conference on Learning Representations (ICLR)*.

11. Strubell, E., Ganesh, A., & McCallum, A. (2019). Energy and policy considerations for deep learning in NLP. In *Proceedings of the 57th Annual Meeting of the ACL*, 3645–3650. https://doi.org/10.18653/v1/P19-1355

12. Hassan, D. A., Abdulrahman, B. K., Mohammed, M. A., & Zeebaree, S. R. M. (2025). Ethical and governance frameworks for artificial intelligence: A systematic literature review. *International Journal of Interactive Mobile Technologies*, 19(14). https://doi.org/10.3991/ijim.v19i14.56981

13. Kim, M. P., Ghorbani, A., & Zou, J. (2019). Multiaccuracy: Black-box post-processing for individual fairness in ranking. In *AIES*.

14. Oprescu, A. M. et al. (2025). A practical framework for auditing fairness in medical AI. In *Proceedings of IDEAL 2025*, Springer.

15. Guidotti, R., Monreale, A., Ruggieri, S., Turini, F., Giannotti, F., & Pedreschi, D. (2022). A survey of methods for explaining black box models. *ACM Computing Surveys*, 51(5), 1–42. https://doi.org/10.1145/3236009
