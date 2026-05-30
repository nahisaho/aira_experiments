# Privacy-Preserving Federated Learning for Multi-Site Clinical Survival Analysis: A Framework Integrating Differential Privacy, Byzantine Robustness, and Communication Efficiency

---

## Abstract

Federated learning (FL) has emerged as a promising paradigm for collaborative training of machine learning models across multiple healthcare institutions without centralizing sensitive patient data. However, deploying FL in real-world clinical settings faces several fundamental challenges: data heterogeneity (non-IID distributions across sites), privacy vulnerabilities against inference attacks, susceptibility to Byzantine (malicious) client behavior, and high communication overhead. This paper presents a comprehensive FL framework tailored for privacy-preserving medical data analysis, demonstrated through a multi-site survival analysis case study using Cox proportional hazards (Cox PH) models.

We design and empirically evaluate five components: (1) a convergence-guaranteed FedAvg baseline with formal analysis under non-IID conditions; (2) FedProx and SCAFFOLD for heterogeneous data stabilization; (3) differentially private federated averaging (DP-FedAvg) with Gaussian noise calibrated via Rényi DP accounting; (4) Top-K gradient sparsification for communication compression; and (5) coordinate-wise median and trimmed-mean aggregation rules for Byzantine fault tolerance. Experiments were conducted on a synthetic multi-site dataset of 1,488 patients across six clinical sites with non-IID age distributions and varying event rates (0.64–0.83).

Key findings include: (a) FedAvg achieves C-index 0.6920 ± 0.0080, nearly matching the centralized upper bound (0.6928 ± 0.0084) and substantially outperforming local-only training (0.6764 ± 0.0184); (b) DP-FedAvg with ε = 5 retains 99.3% of non-private performance, while tight budgets (ε = 0.5) cause severe degradation (C-index 0.5480); (c) Byzantine-robust aggregation (coordinate-wise median) maintains C-index 0.6894 even under sign-flip attacks from 1-in-6 malicious clients, whereas vanilla FedAvg collapses to 0.5515; (d) Top-K (70%) gradient compression reduces communication volume by 30% with negligible performance loss (−0.0001). We critically discuss the limitations arising from synthetic data assumptions and the gap toward real-world deployment. This framework provides a principled foundation for privacy-compliant federated clinical research.

---

## 1. Introduction

The proliferation of electronic health records (EHRs) and multi-institutional clinical registries has created unprecedented opportunities for building large-scale predictive models in medicine. Centralized pooling of patient data from multiple hospitals, however, is routinely prohibited by legal frameworks such as HIPAA, GDPR, and institutional ethics boards, creating a tension between data utility and privacy protection [1].

Federated learning, originally proposed by McMahan et al. [2], addresses this tension by enabling distributed model training where only model parameters—not raw patient data—are transmitted to a central aggregation server. Despite its appeal, naive FL deployment in healthcare settings faces multiple barriers that have been identified and partially addressed in the recent literature.

**Non-IID data heterogeneity** is particularly acute in clinical settings: hospitals differ in patient demographics, treatment protocols, imaging equipment, and disease prevalence [3]. Standard FedAvg assumes independent and identically distributed (IID) data, and its convergence degrades as data divergence grows. FedProx [4] and SCAFFOLD [5] were proposed to mitigate this through proximal regularization and variance reduction, respectively. FedBN [3] introduced local batch normalization to handle feature-level shifts, though it was designed primarily for imaging tasks.

**Privacy amplification** beyond encryption is essential. Differential privacy (DP) provides rigorous formal guarantees through bounded information leakage, quantified by the privacy budget (ε, δ) [6]. Integrating DP with FL—typically by adding calibrated Gaussian noise to client updates—degrades model accuracy in proportion to privacy stringency, creating a fundamental privacy–utility trade-off [7].

**Byzantine resilience** is critical because FL's distributed nature makes it susceptible to malicious clients who inject crafted gradients to corrupt the global model. Robust aggregation rules such as coordinate-wise median [8] and FLTrust [9] have been proposed, but their effectiveness under non-IID conditions remains an active area of investigation.

**Communication efficiency** becomes a bottleneck as the number of clients grows. Gradient compression techniques, including Top-K sparsification and knowledge distillation, reduce the per-round communication cost while preserving convergence [1].

Despite these individual advances, few works have comprehensively evaluated all four challenges together in a clinical survival analysis context. Survival analysis—specifically the Cox proportional hazards model—is a cornerstone of clinical outcomes research (e.g., cancer prognosis, post-transplant monitoring), yet federated implementations remain scarce. This paper makes the following contributions:

1. A unified Flower/PySyft-inspired framework integrating FedAvg, FedProx, SCAFFOLD, DP, Byzantine robustness, and gradient compression for Cox PH models.
2. Empirical evaluation with realistic non-IID multi-site clinical simulation (6 sites, 1,488 patients).
3. Systematic privacy–utility trade-off analysis across five privacy budget levels (ε ∈ {0.5, 1, 2, 5, 10}).
4. Critical self-assessment of the framework's limitations and generalizability to real-world deployment.

---

## 2. Related Work

### 2.1 Federated Learning Fundamentals

Kairouz et al. [1] provide the definitive survey of FL challenges, identifying non-IID data, privacy, communication efficiency, and security as the four primary open problems. Their analysis motivates the design choices in this work. Xu et al. [10] specifically reviewed FL applications in healthcare informatics, covering EHR analysis, medical imaging, and genomics, noting that privacy regulations are the primary driver of FL adoption in medicine.

### 2.2 Non-IID Federated Optimization

Li et al. (FedProx) [4] proved convergence guarantees for FL with non-IID data by adding a proximal term μ||w − w*||² to the local objective, bounding the objective dissimilarity. The proximal term prevents local models from drifting too far from the global model. Karimireddy et al. (SCAFFOLD) [5] identified "client drift" as the root cause of FedAvg's degradation on non-IID data and proposed control variates to correct for it, achieving linear speedup without the smoothness assumptions required by FedProx. Li et al. (FedBN) [3] addressed feature-shift non-IID—common in multi-scanner medical imaging—through local batch normalization, achieving faster convergence than FedAvg and FedProx.

### 2.3 Differential Privacy in Federated Learning

Wei et al. [6] proposed a DP-PFL framework with Rényi DP composition for personalized FL, deriving convergence bounds that reveal an optimal model size for a given privacy level. Zhang et al. [7] analyzed the clipping operation in DP-FedAvg—necessary to bound gradient sensitivity—providing the first rigorous convergence analysis of clipped FedAvg and demonstrating that clipping bias depends on the distribution of client updates. Adnan et al. [11] applied DP-FL to histopathology image analysis on the TCGA dataset, empirically demonstrating that distributed DP training achieves similar performance to centralized training under strong privacy guarantees (ε ≤ 10).

### 2.4 Byzantine Robustness

Cao et al. (FLTrust) [9] proposed bootstrapping trust via a clean root dataset maintained by the server, assigning trust scores to client updates based on directional agreement with the server model update. This approach was shown to be secure against existing attacks and adaptive attacks. Xu et al. [8] proposed SignGuard, exploiting element-wise sign information of gradient vectors to detect model poisoning attacks that evade distance-based defenses. Colosimo and De Rango [12] proposed Median-Krum, combining statistical median with distance-based Krum for improved robustness.

### 2.5 Research Gap

No prior work comprehensively benchmarks FedAvg, FedProx, SCAFFOLD, DP integration, Byzantine robustness, and gradient compression together in a clinical survival analysis setting with Cox PH models. This work addresses that gap.

---

## 3. Methods

### 3.1 Problem Formulation

Consider K clinical sites, each holding a local dataset D_k = {(x_i, t_i, δ_i)}_{i∈D_k}, where x_i ∈ ℝ^d are feature vectors, t_i > 0 are observed event/censoring times, and δ_i ∈ {0,1} indicates event occurrence. The global objective is:

$$\min_{w \in \mathbb{R}^d} F(w) = \sum_{k=1}^{K} \frac{n_k}{n} F_k(w), \quad F_k(w) = \ell_{\text{Cox}}(w; D_k) + \frac{\lambda}{2}\|w\|^2$$

where ℓ_Cox is the negative partial log-likelihood of the Cox PH model, and n = Σ_k n_k. The non-IID assumption means the local data distributions P_k differ across sites.

### 3.2 FedAvg Baseline

In each communication round r, the server broadcasts the global model w^r to all clients. Each client k performs local optimization:

$$w_k^{r+1} = w_k^r - \eta \nabla F_k(w^r)$$

The server aggregates via weighted averaging:

$$w^{r+1} = \sum_{k=1}^{K} \frac{n_k}{n} w_k^{r+1}$$

**Convergence:** Under L-smoothness and bounded gradient dissimilarity (ζ²), FedAvg converges at rate O(1/√(RT)) for non-convex objectives, where R is rounds and T is local steps.

### 3.3 FedProx

To address non-IID divergence, each client solves:

$$\min_w \left[ F_k(w) + \frac{\mu}{2} \|w - w^r\|^2 \right]$$

The proximal term μ/2 ||w - w^r||² prevents excessive local deviation. Convergence is guaranteed for μ > 0 even when clients perform inexact local optimization.

### 3.4 SCAFFOLD

SCAFFOLD introduces client control variates c_k and a server control variate c to correct for client drift:

$$w_k^{r+1} = w_k^r - \eta (\nabla F_k(w_k^r) - c_k + c)$$

$$c_k^{r+1} = c_k - c + \frac{1}{\eta T}(w^r - w_k^{r+1})$$

This achieves linear speedup and eliminates the bounded heterogeneity assumption required by FedProx.

### 3.5 Differential Privacy Integration

We apply the Gaussian mechanism to client model updates. Each client clips their update:

$$\tilde{w}_k = w_k / \max\left(1, \frac{\|w_k\|_2}{C}\right)$$

and adds calibrated noise:

$$\hat{w}_k = \tilde{w}_k + \mathcal{N}(0, \sigma^2 C^2 \mathbf{I})$$

with σ = √(2 ln(1.25/δ)) / ε. This guarantees (ε, δ)-DP per round. Privacy budget composition across R rounds uses Rényi DP accounting [6]:

$$\varepsilon_{\text{total}} \approx \sigma \sqrt{2R \ln(1/\delta)} \cdot q$$

where q = n_k/n is the sampling rate.

**Privacy budget management:** We track cumulative ε per client. When ε_total exceeds a preset threshold, the client opts out of further participation.

### 3.6 Gradient Compression (Top-K Sparsification)

The gradient update Δw_k = w_k - w^r is compressed by retaining only the top-K elements by absolute magnitude:

$$\text{TopK}(\Delta w, k) = \Delta w \odot \mathbf{1}[|\Delta w_i| \geq \text{threshold}_k]$$

This reduces communication from O(d) to O(Kd) floats per round, with compression ratio ρ = K/d = 0.7 in our experiments.

### 3.7 Byzantine-Robust Aggregation

**Coordinate-wise median:**
$$w^{r+1}_j = \text{median}\{w_{k,j}^{r+1} : k=1,\ldots,K\} \quad \forall j \in [d]$$

Robust to up to ⌊(K-1)/2⌋ Byzantine clients under IID conditions.

**Trimmed mean:** Remove the top and bottom ⌊αK⌋ entries per coordinate before averaging (α = 0.2 in our experiments).

### 3.8 Simulation Design

We generated a synthetic multi-site dataset with K=6 sites and n=1,488 total patients using the following data generating process:

- **Covariates:** age ~ N(μ_k, σ_k²) with site-specific means (48–65 years), cancer stage ∈ {1,2,3,4}, WBC count, creatinine, binary treatment, binary comorbidity
- **True log-hazard:** h(x) = 0.04(age-55) + 0.5(stage-2.5) + 0.15·WBC - 1.2·treatment + 0.5·comorbidity + ε, ε ~ N(0, 0.15²)
- **Survival times:** Weibull(shape=1.5, scale=exp(-h(x)·0.5)) × 24 months
- **Censoring:** Exponential(45+5k) months, yielding site-specific event rates 0.64–0.83

Evaluation: 5-fold cross-validation stratified by site, reporting C-index (concordance index) ± standard deviation. All experiments repeated with fixed random seed (42) for reproducibility.

---

## 4. Experiments

### 4.1 Experimental Setup

| Parameter | Value |
|-----------|-------|
| Number of sites (K) | 6 |
| Total patients | 1,488 |
| Features | 6 (age, stage, WBC, creatinine, treatment, comorbidity) |
| Communication rounds | 20–40 |
| Penalization (λ) | 0.05 |
| FedProx μ | 0.1 |
| SCAFFOLD learning rate | 0.005 |
| DP sensitivity (C) | 0.3 |
| DP δ | 10⁻⁵ |
| Compression ratio (ρ) | 0.70 |
| Trimmed mean α | 0.20 |
| CV folds | 5 |
| Random seed | 42 |

### 4.2 Baselines

- **Centralized Cox PH:** All data pooled at a single server (privacy-violating upper bound)
- **Local-only Cox PH:** Single-site training on site-0 data only (lower bound)

### 4.3 Evaluation Metric

The primary metric is Harrell's C-index (concordance index), ranging from 0.5 (random) to 1.0 (perfect). In clinical survival models, C-index > 0.65 is generally considered clinically useful.

---

## 5. Results

### 5.1 Convergence Behavior

![Figure 1: Convergence and cross-validated performance comparison](figures/fig1_convergence_comparison.png)

FedAvg, FedProx, and SCAFFOLD all converge to C-index ≈ 0.693 within 10–15 rounds, nearly identical to the centralized benchmark (0.694 in-sample). DP-FedAvg (ε=5) converges slightly slower but reaches 0.678, while DP-FedAvg (ε=2) shows noticeable oscillation, stabilizing at 0.651. The Byzantine+Median variant converges smoothly to 0.690 despite 1 malicious client.

### 5.2 Cross-Validated Performance

**Table 1: Cross-validated C-index (5-fold, mean ± SD with 95% CI)**

| Method | Mean C-index | SD | 95% CI |
|--------|-------------|-----|--------|
| Centralized (upper bound) | **0.6928** | 0.0084 | [0.6855, 0.7002] |
| FedAvg | 0.6920 | 0.0080 | [0.6850, 0.6991] |
| FedProx (μ=0.1) | 0.6920 | 0.0080 | [0.6850, 0.6991] |
| SCAFFOLD | 0.6920 | 0.0080 | [0.6850, 0.6991] |
| DP-FedAvg (ε=5) | 0.6882 | 0.0068 | [0.6822, 0.6942] |
| DP-FedAvg (ε=2) | 0.5895 | 0.0825 | [0.5172, 0.6618] |
| FedAvg + Byz. Median | 0.6875 | 0.0077 | [0.6807, 0.6943] |
| Local-only (lower bound) | 0.6764 | 0.0184 | [0.6603, 0.6925] |

FedAvg achieves 99.9% of centralized performance (0.6920 vs. 0.6928), a gap of only 0.0008, demonstrating that federation introduces negligible accuracy loss when data are correctly modeled. All FL methods substantially outperform local-only training (Δ = +0.016).

### 5.3 Privacy–Utility Trade-off

![Figure 2: Privacy–Utility trade-off in DP-FedAvg](figures/fig2_privacy_utility_tradeoff.png)

**Table 2: DP-FedAvg performance across privacy budgets**

| Privacy Budget ε | Mean C-index | SD | Δ vs. FedAvg |
|-----------------|-------------|-----|-------------|
| 0.5 (strongest) | 0.5480 | 0.0807 | −0.1440 |
| 1.0 | 0.5919 | 0.0362 | −0.1001 |
| 2.0 | 0.6259 | 0.0449 | −0.0661 |
| 5.0 | 0.6839 | 0.0073 | −0.0081 |
| 10.0 | 0.6807 | 0.0146 | −0.0113 |
| ∞ (no DP) | 0.6920 | 0.0080 | — |

A sharp phase transition occurs between ε=2 and ε=5. At ε≥5, performance is within 1.2% of non-private FL. Below ε=2, noise dominates and variance increases dramatically (SD=0.0825 at ε=2, SD=0.0807 at ε=0.5).

### 5.4 Byzantine Robustness

![Figure 3: Byzantine robustness of aggregation methods](figures/fig3_byzantine_robustness.png)

**Table 3: C-index under Byzantine attacks (5-fold CV)**

| Aggregation | 0 Byzantine | 1 Byzantine | 2 Byzantine |
|------------|-------------|-------------|-------------|
| FedAvg | 0.6920 ± 0.0080 | 0.5515 ± 0.0334 | 0.5100 ± 0.0493 |
| Coord. Median | 0.6915 ± 0.0089 | **0.6894 ± 0.0099** | **0.6862 ± 0.0108** |
| Trimmed Mean (α=0.2) | 0.6916 ± 0.0083 | 0.6882 ± 0.0077 | 0.5858 ± 0.0499 |

FedAvg is catastrophically vulnerable: a single sign-flip Byzantine client (1/6 = 16.7% of clients) reduces C-index from 0.6920 to 0.5515 (−0.1405). Coordinate-wise median retains 99.6% of clean performance with 1 Byzantine client and 99.2% with 2. Trimmed mean (α=0.2) is robust to 1 Byzantine but becomes unstable with 2 (losing 1 client from each tail of 6 effectively removes non-Byzantine clients).

### 5.5 Communication Efficiency

![Figure 5: Communication efficiency and compression](figures/fig5_communication_efficiency.png)

Top-K compression (ρ=0.70) reduces per-round communication by 30% with negligible performance loss: C-index 0.6931 vs. 0.6932 without compression (Δ = −0.0001). This near-lossless compression is expected given the low dimensionality (d=6 features) of the Cox PH model; compression benefits scale with model dimensionality.

### 5.6 Non-IID Data Heterogeneity

![Figure 4: Non-IID data heterogeneity across sites](figures/fig4_noniid_heterogeneity.png)

Site-specific age distributions vary from mean 48 (Site 3) to 65 (Site 4), with event rates ranging from 0.64 to 0.83. Despite this heterogeneity, FL algorithms converge to nearly identical C-indices, suggesting the Cox PH model's structural linearity makes it relatively robust to covariate shifts compared to non-linear models.

---

## 6. Discussion

### 6.1 Interpretation of Results

The near-perfect match between federated and centralized Cox PH performance (ΔC-index = 0.0008) is a notable finding. It implies that for linear models with moderate non-IID shifts, FedAvg's convergence to the global optimum is essentially exact over sufficient rounds. This contrasts with findings on deep neural networks, where client drift causes more significant degradation [3, 5].

FedProx and SCAFFOLD produce identical results to FedAvg in our setting. This is expected: with only 6 features and 6 rounds of local optimization per communication step, client drift is minimal. The advantage of FedProx and SCAFFOLD would be more pronounced with larger models, more local epochs, or more extreme data heterogeneity.

The privacy–utility cliff at ε≈2–5 is consistent with the theoretical analysis of Wei et al. [6], who showed that optimal DP-FL model size and privacy budget are tightly coupled. Below ε=2 in our setting (d=6, n≈250 per site), the noise magnitude σ = 0.3·√(2·ln(1.25/10⁻⁵))/ε ≈ 1.4 at ε=0.5 substantially exceeds the true coefficient magnitudes (~0.3–1.2), making accurate gradient estimation impossible.

### 6.2 Limitations and Critical Self-Assessment

**Synthetic data dependence:** All results derive from a controlled simulation with a known true log-hazard function. The Cox PH model is correctly specified by design—a luxury not available in real clinical data. In practice, proportional hazards assumption violations (time-varying effects, unmeasured confounders) would reduce C-index values by 5–15% and potentially alter the relative performance of FL algorithms.

**Non-IID severity:** The heterogeneity we imposed (age shift of ~17 years across sites) is moderate. More severe distribution shifts—e.g., different cancer types across sites—could cause FedAvg divergence not observed here. Our results may be optimistic for highly heterogeneous federated settings.

**Convergence equivalence of FedAvg/FedProx/SCAFFOLD:** The identical results across three algorithms reflect the low model dimensionality (d=6) and relatively few local steps. For deep learning models (d > 10⁶), the distinction becomes critical.

**Privacy accounting:** We used a simplified Gaussian mechanism analysis. Tighter privacy accounting via DP-SGD with Moments Accountant or PRV accountant would yield lower ε estimates for the same noise level, potentially improving the privacy–utility frontier.

**Byzantine assumption:** We assume a known number of Byzantine clients (1 or 2 out of 6), which is unrealistic. Real Byzantine defenses must be robust without knowledge of the attacker count. FLTrust [9] addresses this through trust bootstrapping, but requires a clean validation dataset at the server—which may itself raise privacy concerns in clinical settings.

**Scalability:** Our simulation with 6 sites and 1,488 patients is small by FL standards. Real federated healthcare networks (e.g., TriNetX [13]) span 100+ institutions with millions of records. Scaling gradient compression, privacy budget management, and Byzantine detection to this scale remains an open problem.

**Real-world generalizability:** We caution against directly interpreting our C-index values as achievable in real clinical deployment. Factors such as data quality variability, missing data patterns, label noise in EHR-derived outcomes, and heterogeneous clinical coding practices would all reduce predictive performance and potentially alter the relative rankings of FL algorithms.

### 6.3 Comparison with Prior Work

Adnan et al. [11] reported that DP-FL on TCGA histopathology data achieved performance comparable to centralized training at ε=10, consistent with our finding that ε≥5 incurs minimal degradation. However, their task (image classification) and model (CNN) differ substantially from survival analysis, limiting direct comparison.

Zhang et al. [7] showed that clipped FedAvg performs surprisingly well even with substantial data heterogeneity for deep architectures, which aligns with our finding that FedAvg nearly matches centralized performance for Cox PH. Their theoretical framework for clipping bias could be extended to the survival analysis setting.

### 6.4 Future Directions

1. **Flower/PySyft deployment:** Migrate the simulation to Flower's actual federated simulation framework for realistic asynchronous communication and system heterogeneity effects.
2. **Personalized FL for survival analysis:** Implement per-FedAvg or pFedMe to learn site-specific Cox models while leveraging shared global structure.
3. **Secure aggregation:** Combine Byzantine robustness with cryptographic secure aggregation to protect against an honest-but-curious server.
4. **Time-varying privacy budgets:** Adaptive ε allocation across rounds using privacy odometers.
5. **Real-world validation:** Apply the framework to publicly available federated datasets such as TCGA or UK Biobank.

---

## 7. Conclusion

This paper presented a comprehensive federated learning framework for privacy-preserving multi-site clinical survival analysis. Through systematic experimentation on a realistic synthetic dataset of 1,488 patients across 6 non-IID sites, we demonstrated that:

1. **FedAvg achieves near-centralized performance** for Cox PH survival analysis (C-index gap: 0.0008), validating FL as a viable privacy-preserving alternative.
2. **Differential privacy with ε≥5 preserves utility** (within 1.2% of non-private performance), while ε<2 causes severe degradation—defining a practical operating range for clinical FL.
3. **Byzantine-robust aggregation (coordinate-wise median) is essential**: a single malicious client causes FedAvg to collapse, while median aggregation retains 99.6% of clean performance.
4. **Top-K gradient compression (70%) achieves 30% communication reduction** with negligible performance loss for low-dimensional models.
5. **FedProx and SCAFFOLD offer no advantage over FedAvg for low-dimensional linear models** but remain important for high-dimensional deep learning settings.

Critically, we emphasize that these results are conditioned on a well-specified synthetic model. Translation to real EHR data requires addressing model misspecification, unmeasured confounding, data quality heterogeneity, and regulatory compliance beyond DP. Nevertheless, this framework provides a rigorous empirical and theoretical foundation for privacy-compliant federated clinical research.

---

## References

[1] Kairouz, P., McMahan, H. B., et al. (2020). Advances and Open Problems in Federated Learning. *Foundations and Trends® in Machine Learning*, 14(1–2), 1–210. https://doi.org/10.1561/2200000083

[2] McMahan, H. B., Moore, E., Ramage, D., Hampson, S., & Agüera y Arcas, B. (2017). Communication-Efficient Learning of Deep Networks from Decentralized Data. *AISTATS 2017*.

[3] Li, X., Jiang, M., Zhang, X., Kamp, M., & Dou, Q. (2021). FedBN: Federated Learning on Non-IID Features via Local Batch Normalization. *ICLR 2021*. https://arxiv.org/abs/2102.07623

[4] Li, T., Sahu, A. K., Zaheer, M., Sanjabi, M., Smola, A., & Smith, V. (2020). Federated Optimization in Heterogeneous Networks (FedProx). *MLSys 2020*. https://arxiv.org/abs/1812.06127

[5] Karimireddy, S. P., Kale, S., Mohri, M., Reddi, S., Stich, S., & Suresh, A. T. (2020). SCAFFOLD: Stochastic Controlled Averaging for Federated Learning. *ICML 2020*. https://arxiv.org/abs/1910.06378

[6] Wei, K., Li, J., Ma, C., Ding, M., Chen, W., Wu, J., Tao, M., & Poor, H. V. (2023). Personalized Federated Learning With Differential Privacy and Convergence Guarantee. *IEEE Transactions on Information Forensics and Security*. https://doi.org/10.1109/TIFS.2023.3293417

[7] Zhang, X., Chen, X., Hong, M., Wu, Z. S., & Yi, J. (2021). Understanding Clipping for Federated Learning: Convergence and Client-Level Differential Privacy. *ICML 2021*. https://www.semanticscholar.org/paper/06757c9d1bf42bad757f8883605f35ad8f02c557

[8] Xu, J., Huang, S.-L., Song, L., & Lan, T. (2022). Byzantine-robust Federated Learning through Collaborative Malicious Gradient Filtering. *IEEE ICDCS 2022*. https://doi.org/10.1109/icdcs54860.2022.00120

[9] Cao, X., Fang, M., Liu, J., & Gong, N. Z. (2021). FLTrust: Byzantine-robust Federated Learning via Trust Bootstrapping. *NDSS 2021*. https://doi.org/10.14722/ndss.2021.24434

[10] Xu, J., Glicksberg, B. S., Su, C., Walker, P., Bian, J., & Wang, F. (2020). Federated Learning for Healthcare Informatics. *Journal of Healthcare Informatics Research*. https://doi.org/10.1007/s41666-020-00082-4

[11] Adnan, M., Kalra, S., Cresswell, J. C., Taylor, G. W., & Tizhoosh, H. R. (2022). Federated learning and differential privacy for medical image analysis. *Scientific Reports*, 12, 1–10. https://doi.org/10.1038/s41598-022-05539-7

[12] Colosimo, F., & De Rango, F. (2023). Median-Krum: A Joint Distance-Statistical Based Byzantine-Robust Algorithm in Federated Learning. *ACM MobiCom Workshop 2023*. https://doi.org/10.1145/3616390.3618283

[13] Ludwig, R. J., et al. (2025). A comprehensive review of methodologies and application to use the real-world data and analytics platform TriNetX. *Frontiers in Pharmacology*. https://doi.org/10.3389/fphar.2025.1516126
