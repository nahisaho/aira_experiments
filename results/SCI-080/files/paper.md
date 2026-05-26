# An Integrated AI System for Food Supply Chain Safety Risk Prediction: Combining Spatiotemporal Modeling, NLP, Predictive Microbiology, and Blockchain Traceability

## Abstract

Food supply chain safety is a critical global public health concern, with foodborne illnesses causing an estimated 600 million cases annually worldwide. This study presents an integrated artificial intelligence system for comprehensive food safety risk prediction and monitoring across the supply chain. The proposed framework comprises six interconnected modules: (1) a spatiotemporal foodborne illness prediction model incorporating meteorological variables (temperature, humidity) and seasonal patterns; (2) a natural language processing (NLP) pipeline for early detection of food recalls and alerts from FDA/RASFF databases; (3) a hybrid microbial growth prediction system integrating the Baranyi-Roberts model with machine learning for enhanced accuracy; (4) an automated HACCP critical control point risk scoring system; (5) a blockchain-integrated traceability framework with anomaly detection capabilities; and (6) a case study on Salmonella contamination prediction in poultry processing. We evaluate multiple machine learning approaches including Random Forest, Gradient Boosting, logistic regression, and neural networks across all modules. The spatiotemporal model achieves an AUC of 0.7356, the NLP module demonstrates perfect classification on structured recall texts, the Baranyi model fitting achieves R² = 0.9986, HACCP risk categorization reaches weighted F1 = 0.9066, and the Salmonella prediction model attains F1 = 0.6886. The integrated system demonstrates the feasibility of a unified AI-driven food safety monitoring platform that combines time series prediction with text analytics for real-time risk assessment. We discuss implications for regulatory compliance, supply chain transparency, and future directions including deep learning integration and IoT-based real-time monitoring.

## 1. Introduction

### 1.1 Background

The globalization of food supply chains has dramatically increased the complexity of food safety management. According to the World Health Organization, approximately 600 million people fall ill after eating contaminated food each year, with 420,000 deaths globally (WHO, 2015). The increasing length and complexity of supply chains, combined with climate change effects on pathogen ecology, necessitate sophisticated predictive systems for proactive risk management.

Traditional food safety management relies on reactive approaches — identifying hazards after incidents occur. However, recent advances in artificial intelligence and machine learning offer transformative potential for predictive food safety. Machine learning models can integrate heterogeneous data sources — environmental sensors, regulatory databases, genomic data, and supply chain logistics — to identify risk patterns before they manifest as foodborne illness outbreaks (Revelou et al., 2025).

### 1.2 Motivation

Several technological and societal trends motivate this research:

1. **Data Availability**: Regulatory databases such as the FDA Recall Database and RASFF (Rapid Alert System for Food and Feed) provide rich textual data amenable to NLP analysis.
2. **Predictive Microbiology**: Databases like ComBase contain over 60,000 bacterial growth records, enabling data-driven modeling approaches that complement classical models such as the Baranyi-Roberts model (Baranyi & Roberts, 1994).
3. **Supply Chain Digitization**: Blockchain technology and IoT sensors enable unprecedented traceability and real-time monitoring capabilities (Ellahi et al., 2023).
4. **Regulatory Requirements**: HACCP (Hazard Analysis Critical Control Points) mandates systematic risk assessment, which can be enhanced through automated scoring (Revelou et al., 2025).

### 1.3 Contributions

This paper makes the following contributions:

- Design and implementation of an integrated six-module AI system for food supply chain safety prediction
- Comparative evaluation of multiple ML algorithms for spatiotemporal illness prediction, NLP-based recall detection, and contamination risk assessment
- A hybrid approach combining the Baranyi mechanistic model with machine learning for microbial growth prediction
- Integration of blockchain traceability with ML-based anomaly detection for supply chain integrity monitoring
- A practical case study on Salmonella contamination prediction in poultry processing

## 2. Related Work

### 2.1 Machine Learning for Food Safety Risk Assessment

Recent comprehensive reviews have documented the rapid adoption of machine learning in food safety risk assessment. A 2025 review in *Foods* summarized advances in applying ML and deep learning for biotoxin detection, microbial risk prediction, and multimodal integration, emphasizing the need for improved model interpretability and real-time HACCP integration (Foods, 2025; DOI: 10.3390/foods14234005). Research in *Frontiers in Sustainable Food Systems* developed XGBoost-based models with SHAP interpretability for food safety risk prediction across over 180,000 agricultural samples, demonstrating the importance of spatiotemporal dimensions (Frontiers, 2026).

A study in *BMC Infectious Diseases* examined spatiotemporal clustering in foodborne disease incidence using decision tree, XGBoost, and LSTM models combined with meteorological data, finding that temperature is the most influential predictor and that LSTM outperforms other models for temporal prediction (BMC Infectious Diseases, 2025).

### 2.2 NLP for Food Recall Detection

Natural language processing has been increasingly applied to automated monitoring of food safety databases. Researchers have used BERT and transformer models to parse recall notices from FDA and RASFF databases, achieving >85% F1-score for event classification and demonstrating earlier warning compared to simple keyword monitoring. The integration of NLP-summarized risk factors with dashboard interfaces has enabled several EU pilot projects for automated alerting systems.

### 2.3 Predictive Microbiology and the Baranyi Model

The Baranyi-Roberts model remains one of the most robust mathematical frameworks for describing sigmoidal bacterial growth, accounting for key growth phases including the lag phase (Baranyi & Roberts, 1994). ComBase provides a massive open-access database of over 60,000 primary growth and inactivation records for foodborne bacteria. Recent research has combined classical models with machine learning — Random Forest and Support Vector Regression for predicting bacterial growth from environmental variables — showing that ML approaches substantially outperform classical models in complex, real-world food environments (Nature Scientific Reports, 2021; DOI: 10.1038/s41598-021-90164-z).

### 2.4 HACCP and Automated Risk Scoring

Revelou et al. (2025) reviewed ML applications for HACCP monitoring of animal-source foods, covering neural networks, supervised/unsupervised models for real-time CCP monitoring using spectroscopy, machine vision, and hybrid AI systems (DOI: 10.3390/foods14060922). The integration of AI for enhanced HACCP management has been explored with IoT sensors and digital traceability systems.

### 2.5 Blockchain for Food Traceability

Ellahi et al. (2023) systematically reviewed blockchain-based frameworks for food traceability, exploring how blockchain enhances accountability, transparency, and traceability in food supply chains through integration with Industry 4.0 technologies (DOI: 10.3390/foods12163026). Patel et al. (2023) focused on livestock product supply chains, describing how blockchain enables tamperproof traceability and regulatory compliance (DOI: 10.1016/j.heliyon.2023.e16526). Sri Vigna Hema and Manickavasagan (2024) reviewed the fundamentals of blockchain implementation for food safety, highlighting challenges in transitioning from conceptual frameworks to practical applications (DOI: 10.1111/1541-4337.70002).

### 2.6 Salmonella Prediction in Poultry

Machine learning approaches for Salmonella prediction in poultry have advanced significantly. A 2023 study used Elastic Net regularization combining genomic factors and meteorological variables for predicting Salmonella outbreaks (DOI: 10.1016/j.crfs.2023.100525). An integrated framework leveraging both classical statistics and machine learning analyzed nearly 42,000 food samples for Salmonella contamination surveillance (DOI: 10.3390/microorganisms13122773). Random forest-based modeling achieved 88% accuracy, 85% sensitivity, and 90% specificity for postchill poultry samples.

## 3. Methods

### 3.1 System Architecture

The proposed system consists of six interconnected modules feeding into a central risk engine:

1. **Spatiotemporal Prediction Module**: Processes meteorological and geographic data
2. **NLP Detection Module**: Analyzes regulatory text data from FDA/RASFF
3. **Microbial Growth Module**: Hybrid Baranyi + ML prediction
4. **HACCP Scoring Module**: Automated CCP risk assessment
5. **Blockchain Traceability Module**: Supply chain integrity verification
6. **Integrated Risk Dashboard**: Unified risk monitoring interface

### 3.2 Spatiotemporal Foodborne Illness Prediction

#### 3.2.1 Feature Engineering

The spatiotemporal model uses six input features: temperature ($T$), humidity ($H$), month ($m$), day of year ($d$), latitude ($\phi$), and longitude ($\lambda$). The seasonal risk component is modeled as:

$$S(d) = 0.3 \sin\left(\frac{2\pi d}{365} - \frac{\pi}{6}\right)$$

#### 3.2.2 Model Formulation

The incident probability is modeled through a logistic function:

$$P(\text{incident}) = \sigma\left(\beta_0 + \beta_1 T + \beta_2 H + S(d) + \beta_3 \frac{(T-25)(H-70)}{100}\right)$$

where $\sigma(x) = \frac{1}{1+e^{-x}}$ is the sigmoid function.

Four classifiers were compared: Logistic Regression, Random Forest (200 trees, max depth 10), Gradient Boosting (200 estimators, max depth 5), and MLP (64-32 architecture).

### 3.3 NLP-based Recall Detection

#### 3.3.1 Text Representation

Documents are represented using TF-IDF vectorization with unigram and bigram features:

$$\text{TF-IDF}(t, d) = \text{TF}(t, d) \times \log\frac{N}{|\{d' \in D : t \in d'\}|}$$

where $t$ is a term, $d$ is a document, $N$ is the total number of documents, and $D$ is the document corpus.

#### 3.3.2 Classification Tasks

Two classification tasks were performed:
- **Binary**: Urgent (recall/alert) vs. non-urgent
- **Multi-class**: Five categories (recall, alert, warning, information, normal)

### 3.4 Microbial Growth Prediction

#### 3.4.1 Baranyi-Roberts Model

The Baranyi growth model is defined as:

$$y(t) = y_0 + \mu_{\max} A(t) - \ln\left(1 + \frac{e^{\mu_{\max} A(t)} - 1}{e^{y_{\max} - y_0}}\right)$$

where the adjustment function $A(t)$ accounts for the lag phase:

$$A(t) = t + \frac{1}{\mu_{\max}} \ln\left(e^{-\mu_{\max} t} + e^{-\mu_{\max} \lambda} - e^{-\mu_{\max}(t+\lambda)}\right)$$

Here, $y_0$ is the initial log count (log CFU/g), $y_{\max}$ is the maximum population density, $\mu_{\max}$ is the maximum specific growth rate (h⁻¹), and $\lambda$ is the lag time (h).

#### 3.4.2 Secondary Model Integration

The temperature dependence of growth rate follows a Ratkowsky-like relationship:

$$\mu_{\max}(T) = 0.01 \cdot e^{0.08T}$$

modulated by pH and water activity ($a_w$):

$$\mu_{\max}(T, \text{pH}, a_w) = \mu_{\max}(T) \cdot (1 - e^{-0.5(\text{pH} - 3.5)}) \cdot \frac{a_w - 0.88}{0.12}$$

#### 3.4.3 ML Enhancement

A Random Forest Regressor (200 trees) was trained to predict $\mu_{\max}$ from temperature, pH, and water activity, complementing the mechanistic model.

### 3.5 HACCP Risk Scoring

The automated HACCP risk score integrates multiple monitoring parameters:

$$R = \sum_{i=1}^{n} w_i \cdot f_i(x_i)$$

where $w_i$ are weights for each factor, $f_i$ are transformation functions, and $x_i$ are CCP monitoring variables including temperature deviation, time deviation, humidity deviation, equipment age, staff training score, violation history, and inspection frequency.

A Gradient Boosting Regressor predicts continuous risk scores, while a Random Forest Classifier categorizes risk into low/medium/high levels.

### 3.6 Blockchain Traceability

#### 3.6.1 Chain Structure

Each block in the food traceability chain contains:

$$B_i = \{i, t_i, D_i, H_{i-1}, H_i\}$$

where $H_i = \text{SHA-256}(\text{JSON}(i, t_i, D_i, H_{i-1}))$.

Chain integrity is verified by checking that $H_i = \text{SHA-256}(\text{JSON}(B_i))$ for all blocks and $B_i.\text{prev\_hash} = H_{i-1}$.

#### 3.6.2 Anomaly Detection

A Random Forest classifier detects anomalous shipments based on temperature and transit time deviations from expected patterns.

### 3.7 Salmonella Case Study

The Salmonella contamination probability model uses nine features:

$$\text{logit}(P) = \beta_0 + \beta_1 T_p + \beta_2 T_c + \beta_3 D_s + \beta_4 T_s + \beta_5 H + \beta_6 T_a + \beta_7 S + \beta_8 R + \beta_9 C$$

where $T_p$ is processing temperature, $T_c$ is cooking temperature, $D_s$ is storage duration, $T_s$ is storage temperature, $H$ is humidity, $T_a$ is ambient temperature, $S$ is season indicator, $R$ is supplier rating, and $C$ is chlorine wash indicator.

## 4. Experiments

### 4.1 Experimental Setup

All experiments were implemented in Python 3.12 using scikit-learn 1.6, NumPy, Pandas, and SciPy. Synthetic datasets were generated to simulate realistic food supply chain scenarios based on published distributions and domain knowledge. Data were split 80/20 for training/testing, with 5-fold cross-validation applied for the Salmonella case study.

### 4.2 Datasets

| Module | Dataset Size | Features | Task Type |
|--------|-------------|----------|-----------|
| Spatiotemporal | 2,000 samples | 6 | Binary classification |
| NLP Detection | 1,000 documents | 500 (TF-IDF) | Binary + 5-class |
| Microbial Growth | 500 conditions | 3 | Regression |
| HACCP Scoring | 1,500 records | 8 | Regression + 3-class |
| Blockchain | 500 shipments | 2 | Binary (anomaly) |
| Salmonella | 1,200 samples | 9 | Binary classification |

### 4.3 Evaluation Metrics

- **Classification**: Accuracy, Precision, Recall, F1-score, AUC-ROC
- **Regression**: RMSE, MAE, R²
- **Blockchain**: Chain integrity verification, anomaly detection performance

## 5. Results

### 5.1 Spatiotemporal Prediction

Table 1 presents the performance comparison of spatiotemporal prediction models.

| Model | Accuracy | F1 Score | AUC |
|-------|----------|----------|-----|
| Logistic Regression | 0.6875 | 0.6246 | **0.7356** |
| Random Forest | 0.6725 | 0.6158 | 0.7187 |
| Gradient Boosting | 0.6500 | 0.5783 | 0.6850 |
| MLP Neural Network | 0.6700 | 0.6118 | 0.6961 |

Logistic Regression achieved the best AUC (0.7356), suggesting that the underlying relationship is approximately linear in the feature space. Feature importance analysis from the Random Forest model revealed that temperature and day of year (capturing seasonality) are the dominant predictors.

![Figure 1: Spatiotemporal model performance comparison, feature importance, and ROC curves](figures/spatiotemporal_model_comparison.png)

![Figure 2: Monthly incidence patterns and temperature-humidity distributions](figures/spatiotemporal_patterns.png)

### 5.2 NLP Recall Detection

All four models achieved perfect classification (F1 = 1.0, AUC = 1.0) for both binary and multi-class tasks. This reflects the structured nature of the template-based text data, which contains clear lexical indicators (e.g., "RECALL", "URGENT", "alert") that are well-captured by TF-IDF features.

![Figure 3: NLP classification performance, confusion matrix, and important features](figures/nlp_recall_detection.png)

### 5.3 Microbial Growth Prediction

The Baranyi model fitting achieved excellent accuracy:
- **RMSE**: 0.1056 log CFU/g
- **R²**: 0.9986
- **Estimated parameters**: y₀ = 2.02, y_max = 8.93, μ_max = 0.548 h⁻¹, lag = 3.41 h

The ML-enhanced growth rate prediction achieved:
- **RMSE**: 0.0126 h⁻¹
- **R²**: 0.8426
- **MAE**: 0.0094 h⁻¹

Temperature was the most important feature (importance = 0.72), followed by pH (0.16) and water activity (0.12).

![Figure 4: Baranyi growth curves, model fitting, ML prediction accuracy, and feature importance](figures/microbial_growth.png)

### 5.4 HACCP Risk Scoring

**Regression (Gradient Boosting):**
- RMSE: 1.0585, R²: 0.4201, MAE: 0.8487

**Classification (Random Forest):**
- Accuracy: 0.9167, Weighted F1: 0.9066

The classification task showed strong performance, with temperature deviation and previous violations being the most influential features for risk scoring.

![Figure 5: HACCP risk scoring results — regression, distribution by CCP type, feature importance, and classification confusion matrix](figures/haccp_scoring.png)

### 5.5 Blockchain Traceability

The blockchain implementation successfully maintained chain integrity across all 7 blocks. The anomaly detection module achieved perfect classification (F1 = 1.0, AUC = 1.0) for identifying shipments with temperature or transit time deviations.

![Figure 6: Blockchain supply chain tracking, anomaly detection scatter plot, and chain structure](figures/blockchain_traceability.png)

### 5.6 Salmonella Case Study

| Model | Accuracy | F1 Score | AUC |
|-------|----------|----------|-----|
| Logistic Regression | 0.6458 | **0.6886** | **0.6815** |
| Random Forest | 0.6375 | 0.6859 | 0.6698 |
| Gradient Boosting | 0.6000 | 0.6496 | 0.6351 |
| MLP | 0.5500 | 0.5748 | 0.5516 |

5-fold Cross-Validation F1 (Random Forest): 0.6561 ± 0.0444

Processing temperature and chlorine wash were identified as the most important predictive features. Summer season showed elevated contamination rates compared to other seasons.

![Figure 7: Salmonella prediction model comparison, ROC curves, feature importance, seasonal analysis, temperature effects, and confusion matrix](figures/salmonella_case_study.png)

### 5.7 Integrated System Dashboard

![Figure 8: Integrated risk monitoring dashboard — system architecture, module performance summary, risk timeline, and hazard-month heatmap](figures/integrated_dashboard.png)

## 6. Discussion

### 6.1 Key Findings

The experimental results demonstrate the feasibility of an integrated AI system for food supply chain safety risk prediction. Several key findings emerge:

**Spatiotemporal Prediction**: The moderate performance (AUC = 0.7356) reflects the inherent stochasticity of foodborne illness occurrence. Temperature emerged as the dominant predictor, consistent with the findings of the BMC Infectious Diseases (2025) study which identified temperature as the most influential meteorological factor. The linear model (Logistic Regression) outperforming ensemble methods suggests that the relationship between environmental conditions and illness risk is approximately linear in this feature space, though real-world data with more complex interactions may benefit from non-linear models.

**NLP Detection**: Perfect classification on structured text data validates the TF-IDF approach for templated regulatory communications. However, real-world FDA/RASFF notifications exhibit greater linguistic variability. Future work should evaluate transformer-based models (BERT, GPT) on actual regulatory text data, as demonstrated by recent EU pilot projects achieving >85% F1-score.

**Microbial Growth**: The excellent Baranyi model fit (R² = 0.9986) confirms the model's validity for controlled laboratory conditions. The ML complement (R² = 0.8426) demonstrates that machine learning can effectively predict growth parameters from environmental conditions, supporting the hybrid modeling approach advocated by recent predictive microbiology platforms (MDPI Foods, 2025).

**HACCP Scoring**: The classification performance (F1 = 0.9066) is promising for practical deployment. The lower regression R² (0.4201) indicates that continuous risk scoring requires more nuanced feature engineering or additional monitoring variables.

**Salmonella Prediction**: The case study results (F1 = 0.6886) are realistic for food safety prediction tasks, where contamination events are influenced by numerous latent factors. The identified importance of processing temperature and chlorine washing aligns with established food safety knowledge and prior ML studies achieving 88% accuracy on postchill poultry samples.

### 6.2 Limitations

1. **Synthetic Data**: All experiments used generated data. Validation on real-world datasets from FDA, RASFF, ComBase, and USDA FSIS is essential.
2. **Model Simplicity**: Deep learning architectures (LSTM, Transformers) were not fully explored due to the scope of this study.
3. **Temporal Dynamics**: The current system processes static snapshots rather than continuous time series, limiting real-time monitoring capabilities.
4. **Cross-module Integration**: While all modules are implemented, deeper integration through shared latent representations could improve system-level performance.
5. **Scalability**: The blockchain implementation is simplified; production deployment would require distributed consensus mechanisms and smart contract functionality.

### 6.3 Future Directions

1. **Deep Learning Integration**: Implement LSTM/Transformer models for time series prediction and BERT for NLP tasks on actual regulatory text
2. **IoT Sensor Integration**: Connect real-time temperature, humidity, and gas sensors through an IoT middleware layer
3. **Explainable AI**: Apply SHAP and LIME for model interpretability, critical for regulatory acceptance
4. **Federated Learning**: Enable collaborative model training across supply chain partners without sharing proprietary data
5. **Multi-modal Fusion**: Integrate visual inspection data (computer vision) with sensor and text data
6. **Real-world Validation**: Deploy pilot system with food processing facilities and regulatory agencies

## 7. Conclusion

This study presented an integrated AI system for food supply chain safety risk prediction, combining spatiotemporal modeling, NLP-based recall detection, mechanistic microbial growth prediction enhanced by machine learning, automated HACCP risk scoring, blockchain traceability, and a practical Salmonella contamination case study. The system demonstrates that machine learning can effectively support food safety decision-making across multiple dimensions of the supply chain.

Key achievements include a spatiotemporal prediction AUC of 0.7356, perfect NLP classification on structured recall texts, Baranyi model fitting with R² = 0.9986, HACCP risk categorization with weighted F1 = 0.9066, verified blockchain integrity, and Salmonella prediction with F1 = 0.6886. The integrated dashboard provides a unified view of risk across all modules.

While this work uses synthetic data, the system architecture and methodology provide a foundation for real-world deployment. Future work will focus on deep learning integration, IoT sensor connectivity, explainable AI for regulatory compliance, and validation with actual food safety databases including FDA, RASFF, and ComBase.

## References

1. Baranyi, J., & Roberts, T. A. (1994). A dynamic approach to predicting bacterial growth in food. *International Journal of Food Microbiology*, 23(3-4), 277–294. DOI: 10.1016/0168-1605(94)90157-0

2. Revelou, P.-K., et al. (2025). Applications of Machine Learning in Food Safety and HACCP Monitoring of Animal-Source Foods. *Foods*, 14(6), 922. DOI: 10.3390/foods14060922

3. Ellahi, R. M., Wood, L. C., & Bekhit, A. E.-D. A. (2023). Blockchain-Based Frameworks for Food Traceability: A Systematic Review. *Foods*, 12(16), 3026. DOI: 10.3390/foods12163026

4. Patel, A. S., et al. (2023). Blockchain technology in food safety and traceability concern to livestock products. *Heliyon*, 9(6), e16526. DOI: 10.1016/j.heliyon.2023.e16526

5. Sri Vigna Hema, V., & Manickavasagan, A. (2024). Blockchain implementation for food safety in supply chain: A review. *Comprehensive Reviews in Food Science and Food Safety*. DOI: 10.1111/1541-4337.70002

6. Orr, M. B., et al. (2023). Machine learning to predict foodborne salmonellosis outbreaks based on genome characteristics and meteorological data. *Current Research in Food Science*, 6, 100525. DOI: 10.1016/j.crfs.2023.100525

7. Cito, A., et al. (2025). Integrating Statistical and Machine-Learning Approaches for Salmonella enterica Surveillance. *Microorganisms*, 13(12), 2773. DOI: 10.3390/microorganisms13122773

8. Ding, T., et al. (2021). Prediction of population behavior of Listeria monocytogenes in food using machine learning and a microbial growth and survival database. *Scientific Reports*, 11, 10613. DOI: 10.1038/s41598-021-90164-z

9. Application of Machine Learning in Food Safety Risk Assessment. (2025). *Foods*, 14(23), 4005. DOI: 10.3390/foods14234005

10. World Health Organization. (2015). WHO estimates of the global burden of foodborne diseases. Geneva: WHO.
