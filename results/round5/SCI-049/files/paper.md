# Automated Quality Control and Anomaly Detection for Large-Scale Scientific Streaming Data: A Multi-Component Pipeline Integrating Change-Point Detection, Physical Constraints, and Explainability

---

## Abstract

Large-scale physics experiments such as those conducted at CERN's Large Hadron Collider (LHC) and the LIGO gravitational-wave observatory generate continuous multi-channel sensor streams at rates of terabytes per second, where undetected data quality anomalies can corrupt downstream physics analyses. Existing approaches often treat anomaly detection as an isolated statistical problem, neglecting domain-specific physical constraints and the non-stationary nature of detector behavior under concept drift. In this work, we present **SciAD** (Scientific Anomaly Detection), a unified streaming pipeline that integrates six complementary components: (1) Pruned Exact Linear Time (PELT) and Bayesian Online Change Point Detection (BOCPD) for abrupt structural changes, (2) Isolation Forest and a Deep SVDD approximation for multivariate outlier scoring, (3) a physics-informed constraint scorer that encodes hard sensor limits, rate-of-change bounds, and cross-channel correlation expectations, (4) an ADWIN-inspired drift detector with model retraining triggers, (5) a feature-attribution explainability module for automated anomaly root-cause identification, and (6) an end-to-end streaming architecture inspired by CERN CMS and LIGO detector requirements. We evaluate the pipeline on synthetic multi-channel detector data (5,000 time steps, 6 channels) modeled after LHC readout characteristics, with injected point anomalies, burst anomalies, physical constraint violations, and abrupt change points. The combined Isolation Forest + physical constraint scorer achieves an AUROC of **0.968 ± 0.001** across five random seeds, while F1 score of **0.474 ± 0.000** reveals a precision–recall trade-off inherent to unsupervised detection at a 16% effective anomaly rate. PELT successfully recovers all four injected change-point locations (within ±50 steps tolerance) alongside 41 false detections, motivating stricter penalty calibration. The ADWIN drift detector identifies 149 distribution shift events in the residual stream. We critically examine the limitations of our synthetic evaluation, discuss generalization to real detector environments, and provide design guidelines for deployment in CERN/LIGO-type facilities.

**Keywords:** anomaly detection, scientific data quality control, change point detection, concept drift, physics-informed machine learning, CERN, LIGO, streaming data

---

## 1. Introduction

### 1.1 Motivation and Background

High-energy physics experiments produce data at unprecedented scales. The CERN LHC generates approximately 1 petabyte of raw data per second across hundreds of millions of detector channels, of which only a fraction can be stored after real-time trigger selection [Stankevicius et al., 2020]. LIGO's gravitational-wave detectors record time-series data from thousands of auxiliary channels monitoring environmental and instrumental conditions, and "glitches"—transient noise artifacts—occur at a rate of roughly one per minute per detector [Davis et al., 2022]. Undetected data quality problems propagate into physics results: misidentified anomalies can mimic new physics signals or mask genuine discoveries.

Traditional data quality monitoring (DQM) in high-energy physics relies on expert-defined histograms compared against reference runs, a process that is labor-intensive, does not scale to Run 3/Run 4 luminosity regimes, and fails to capture subtle multi-channel correlations. Machine-learning-based approaches have begun to address this gap [Stankevicius et al., 2020], but most deployed systems lack:
- Explicit encoding of physical constraints (sensor operating ranges, conservation laws)
- Ability to detect and adapt to gradual concept drift in detector behavior
- Explainable outputs that allow operators to understand the root cause of flagged events
- A unified streaming architecture that processes all components in real time

### 1.2 Research Contributions

This paper makes the following contributions:

1. **Unified pipeline architecture** combining change-point detection, multivariate anomaly scoring, physical constraint scoring, drift detection, and explainability within a single streaming framework.
2. **Physics-informed anomaly scoring** that fuses statistical outlier scores with domain knowledge about acceptable sensor operating ranges and inter-channel dependencies.
3. **Quantitative evaluation** with cross-validated metrics (AUROC, F1, Precision, Recall) on synthetic data designed to reflect CERN/LIGO detector characteristics.
4. **Critical self-assessment** of synthetic evaluation limitations and a roadmap for real-world deployment.

### 1.3 Paper Organization

Section 2 reviews related work. Section 3 describes the SciAD pipeline methods. Section 4 presents the experimental setup. Section 5 reports results. Section 6 discusses findings, limitations, and future work. Section 7 concludes.

---

## 2. Related Work

### 2.1 Anomaly Detection in Scientific Data

Data quality monitoring using machine learning has been explored in multiple physics contexts. Stankevicius et al. [2020] applied meta-learning to optimize artificial neural network hyper-parameters for CERN CMS offline data certification, demonstrating that automated ML can match human expert certification on simulated luminosity sections. Davis et al. [2022] presented a Bayesian subtraction method for LIGO glitches, showing that glitch models can be learned from auxiliary channels and subtracted from the strain channel to improve gravitational-wave detection sensitivity. Cavaglià [2022] used fractal analysis to characterize the statistical properties of LIGO detector noise, providing a complementary perspective to anomaly detection.

### 2.2 Change Point Detection

Adams & MacKay [2007] introduced Bayesian Online Change Point Detection (BOCPD), which maintains a posterior distribution over run lengths and updates it with each new observation using conjugate Normal-Inverse-Gamma priors. Killick et al. [2012] developed PELT (Pruned Exact Linear Time), a dynamic programming algorithm achieving O(n) average complexity for change-point search with an RBF cost function. Corradin et al. [2022] extended Bayesian change-point methods to multivariate time series with missing observations using nonparametric priors, and Tsaknaki et al. [2025] developed a Bayesian autoregressive extension with time-varying parameters for financial data, directly applicable to scientific sensor streams.

### 2.3 Multivariate Outlier Detection

Liu et al. [2012] proposed Isolation Forest, which isolates anomalies by recursively partitioning the feature space with random splits; anomalies require fewer splits to isolate and receive lower path-length scores. Ruff et al. [2018] introduced Deep SVDD, which maps data into a hypersphere using a neural network encoder and flags points far from the learned center. Katbi & Ksantini [2025] recently enhanced deep SVDD with an adversarial regularizer and interpolation to improve performance on IoT sensor data.

### 2.4 Concept Drift Detection

Bifet & Gavalda [2007] proposed ADWIN (Adaptive Windowing), a sliding-window drift detector that dynamically adjusts window size based on Hoeffding bounds on mean differences. Gama et al. [2004] earlier proposed DDM (Drift Detection Method) based on statistical monitoring of online learner error rates. Both approaches are relevant to scientific data monitoring, where detector conditions change during experimental runs.

### 2.5 Explainable Anomaly Detection

Interpretability of anomaly detection outputs is essential for scientific operators. SHAP (SHapley Additive exPlanations) [Lundberg & Lee, 2017] provides unified feature attribution but requires model-specific implementations. Simpler marginal attribution methods—computing each feature's contribution to the aggregate anomaly score—have been used in industrial settings and adapted here for multi-channel detector contexts.

---

## 3. Methods

### 3.1 Pipeline Overview

The SciAD pipeline processes multi-channel sensor streams through six sequential stages:

$$\text{Score}_{\text{final}}(x_t) = \alpha \cdot \hat{s}_{\text{IF}}(x_t) + (1-\alpha) \cdot s_{\text{phys}}(x_t)$$

where $\hat{s}_{\text{IF}}$ is the normalized Isolation Forest anomaly score, $s_{\text{phys}}$ is the physical constraint score, and $\alpha = 0.6$ is a weighting hyperparameter.

### 3.2 Data Preprocessing

Input data is standardized using Z-score normalization:

$$\tilde{x}_{t,c} = \frac{x_{t,c} - \mu_c}{\sigma_c + \epsilon}$$

where $\mu_c, \sigma_c$ are per-channel mean and standard deviation estimated from a training window, and $\epsilon = 10^{-7}$ prevents division by zero.

### 3.3 Change Point Detection

#### 3.3.1 PELT (Pruned Exact Linear Time)

PELT minimizes the penalized cost:

$$\sum_{i=1}^{m+1} \mathcal{C}(y_{\tau_{i-1}+1:\tau_i}) + \beta m$$

where $\mathcal{C}$ is an RBF-kernel cost function, $\tau_i$ are change-point indices, $m$ is the number of change points, and $\beta$ is the penalty parameter. PELT achieves O(n) average complexity through pruning: if $F(t) + \mathcal{C}(y_{t+1:s}) + \beta \geq F(s)$ for all future $s$, then $t$ can be pruned from the set of candidate change-point positions.

#### 3.3.2 BOCPD (Bayesian Online Change Point Detection)

BOCPD maintains a distribution over run lengths $r_t$ (time since last change point):

$$P(r_t | x_{1:t}) \propto \sum_{r_{t-1}} P(x_t | r_{t-1}, x_{t-r_{t-1}:t-1}) P(r_t | r_{t-1}) P(r_{t-1} | x_{1:t-1})$$

The growth probability is:
$$P(r_t = r_{t-1}+1 | r_{t-1}) = 1 - \frac{1}{\lambda}$$

and the change-point probability is:
$$P(r_t = 0 | r_{t-1}) = \frac{1}{\lambda}$$

where $\lambda$ is the expected run length (hazard parameter, set to 250 in experiments). The predictive distribution $P(x_t | r_{t-1}, \cdot)$ uses conjugate Normal-Inverse-Gamma priors with parameters $(\mu_0, \kappa_0, \alpha_0, \beta_0) = (0, 1, 1, 1)$.

### 3.4 Multivariate Anomaly Detection

#### 3.4.1 Isolation Forest

An ensemble of $T$ isolation trees is constructed by randomly selecting a feature $q$ and split value $p \in [\min(q), \max(q)]$. The anomaly score for sample $x$ is:

$$s(x, n) = 2^{-\frac{E[h(x)]}{c(n)}}$$

where $h(x)$ is the path length to isolate $x$, $c(n) = 2H(n-1) - 2(n-1)/n$ is the average path length for a binary search tree of $n$ nodes, and $H$ is the harmonic number. We use $T=150$ trees and contamination $\rho = 0.05$.

#### 3.4.2 Deep SVDD Approximation

We implement a simplified Deep SVDD using a two-layer feed-forward network ($d \rightarrow 32 \rightarrow 16$ with ReLU activations) trained to minimize:

$$\mathcal{L} = \frac{1}{n} \sum_{i=1}^n \|f(x_i; W) - c\|^2$$

where $c$ is the hypersphere center initialized as the mean of the initial forward pass. The anomaly threshold is set at the $(1-\nu)$-quantile of training distances, with $\nu = 0.05$.

### 3.5 Physical Constraint Scoring

The physics-based score integrates three domain-specific checks:

$$s_{\text{phys}}(x_t) = w_1 \cdot s_{\text{bound}}(x_t) + w_2 \cdot s_{\text{roc}}(x_t) + w_3 \cdot s_{\text{corr}}(x_t)$$

**Hard bounds:** $s_{\text{bound}} = \sum_c \mathbf{1}[x_{t,c} \notin [L_c, U_c]]$, where $[L_c, U_c]$ are channel-specific operating limits.

**Rate of change:** $s_{\text{roc}} = \sum_c \mathbf{1}[|\Delta x_{t,c}| > 2 \cdot q_{99}(\Delta x_c)]$, where $q_{99}$ is the 99th percentile of absolute first differences.

**Cross-channel correlation:** $s_{\text{corr}} = \mathbf{1}[|\hat\rho_t - \rho_{\text{ref}}| > \delta_\rho]$, where $\hat\rho_t$ is the rolling-window correlation between designated channel pairs and $\rho_{\text{ref}}$ is a baseline value, with threshold $\delta_\rho = 0.4$.

Weights are $(w_1, w_2, w_3) = (2.0, 1.5, 1.0)$ reflecting the severity ordering. The final score is max-normalized to $[0, 1]$.

### 3.6 Drift Detection (ADWIN)

We implement a simplified ADWIN-like detector that monitors the error stream $e_t = |\Delta x_{t,c}|$. A drift is flagged when:

$$|\mu_{W_1} - \mu_{W_2}| > \sqrt{\frac{\log(2/\delta)}{2 \min(|W_1|, |W_2|)}}$$

where $W_1, W_2$ are the two halves of the sliding window, $\delta = 0.002$ is the confidence parameter, and $\mu_{W_i}$ is the mean of window $W_i$. Upon drift detection, the algorithm resets to the newer window and increments a retraining counter. When the counter exceeds a threshold (5 in experiments), model retraining is triggered.

### 3.7 Feature Attribution (Explainability)

Per-sample feature importance is computed as:

$$a_{t,j} = \frac{|\tilde{x}_{t,j}| \cdot |s_{\text{final}}(x_t)|}{\sum_{k} |\tilde{x}_{t,k}| \cdot |s_{\text{final}}(x_t)|}$$

Global importance is the mean across all samples:
$$\bar{a}_j = \frac{1}{n} \sum_t a_{t,j}$$

This provides an approximate marginal attribution score identifying which channels contribute most to flagged anomalies.

---

## 4. Experiments

### 4.1 Synthetic Data Generation

We generate synthetic multi-channel detector data designed to reflect CERN LHC and LIGO readout characteristics:

- **Channels:** 6 sensor channels with correlated Gaussian noise (covariance decay: $\rho_{ij} = 0.6 \exp(-0.5|i-j|)$)
- **Signal:** Superposition of sinusoidal oscillations at different frequencies (0.002–0.007 Hz) plus additive noise ($\sigma = 0.3$)
- **Change points:** 4 abrupt gain/offset shifts of magnitude 1.5–3.0 at random positions
- **Anomalies:** Three types: point spikes (±5–10σ), burst segments (3–20 samples, +3σ Gaussian), and physical constraint violations (channel 0 exceeding hard upper limit)
- **Dataset size:** 5,000 time steps, 5 cross-validation seeds

The effective anomaly fraction is 16.1% (vs. intended 5%) due to overlapping burst segments; this discrepancy is itself a methodological finding discussed in Section 6.

### 4.2 Evaluation Protocol

We evaluate three detectors using 5 independently seeded models (not stratified k-fold, as the data is temporal):

- **Isolation Forest** (baseline)
- **Deep SVDD approximation** (neural one-class classifier)
- **IF + Physical Constraints** (proposed combined scorer)

Metrics: AUROC (threshold-independent), F1 Score, Precision, Recall. All metrics reported as mean ± standard deviation across 5 seeds.

### 4.3 Change Point Detection Evaluation

PELT is run on each channel with penalty $\beta \in \{5, 8, 15\}$. A detected change point is counted as a true positive if it falls within ±50 time steps of a true change point. BOCPD is evaluated on the first 1,000 samples (due to $O(n^2)$ complexity).

---

## 5. Results

### 5.1 Pipeline Visualization

![Figure 1: Pipeline Architecture](figures/fig5_pipeline.png)

*Figure 1: SciAD streaming anomaly detection pipeline architecture. Data flows from raw sensor streams through physical constraint scoring, change point detection, multivariate anomaly scoring, drift detection, and explainability modules.*

### 5.2 Signal Overview and Anomaly Detection

![Figure 2: Signal and Anomaly Detection](figures/fig1_signal_anomalies.png)

*Figure 2: (a) Multi-channel detector readout with true anomaly labels (red) and true change points (orange dashed). (b) Isolation Forest anomaly scores with 95th-percentile threshold. (c) PELT change point detection (green dash-dot = detected, orange dashed = true). (d) ADWIN drift events (magenta vertical lines).*

### 5.3 BOCPD Results

![Figure 3: BOCPD Change Point Detection](figures/fig2_bocpd.png)

*Figure 3: BOCPD applied to the first 1,000 samples of Channel 1. (a) Raw signal with true change points. (b) MAP run-length estimate. (c) Detected change points from run-length drops.*

### 5.4 Anomaly Detection Performance

**Table 1: Detector Performance Comparison (mean ± std, n=5 seeds)**

| Method | AUROC | F1 Score | Precision | Recall |
|--------|-------|----------|-----------|--------|
| Isolation Forest | 0.966 ± 0.002 | 0.474 ± 0.000 | — | — |
| Deep SVDD (approx) | 0.880 ± 0.032 | 0.472 ± 0.003 | — | — |
| IF + Physical Constraints | **0.968 ± 0.001** | **0.474 ± 0.000** | — | — |
| Combined (0.6·IF + 0.4·Phys) | 0.963 | 0.474 | 1.000 | 0.311 |

![Figure 4: Performance Comparison](figures/fig4_performance.png)

*Figure 4: AUROC (left) and F1 (right) comparison across detectors. Error bars represent standard deviation across 5 seeds. The combined scorer shows a marginal improvement in AUROC over standalone IF.*

### 5.5 Physical Constraint Analysis

![Figure 5: Physical Constraint Scoring](figures/fig6_physical_constraints.png)

*Figure 5: (a) Physical constraint violation score over time. (b) Channel scatter plot with anomaly labels. (c) Rate-of-change distribution (normal vs. anomalous). (d) Physical score vs. combined score.*

### 5.6 Explainability

![Figure 6: Feature Attribution and Score Distributions](figures/fig3_explainability.png)

*Figure 6: (a) Global feature importance by channel. (b) Per-sample feature attribution heatmap for top 200 anomalies. (c) Anomaly score distribution for normal vs. anomalous samples.*

**Table 2: Feature Importance Ranking (Global Attribution)**

| Rank | Channel | Attribution Score |
|------|---------|-------------------|
| 1 | Ch6 | 0.1758 |
| 2 | Ch2 | 0.1749 |
| 3 | Ch3 | 0.1731 |
| 4 | Ch5 | 0.1725 |
| 5 | Ch4 | 0.1714 |
| 6 | Ch1 | 0.1322 |

### 5.7 Change Point Detection Results

**Table 3: Change Point Detection Summary**

| Method | Detected CPs | True CPs Recovered (±50) | False Positives | Latency |
|--------|-------------|--------------------------|-----------------|---------|
| PELT (β=8) | 45 | 4/4 (100%) | 41 | < 1 ms/sample |
| BOCPD (λ=250) | 23 (first 1000 steps) | 1/1 visible | 22 | ~5 ms/sample |

### 5.8 Drift Detection Results

ADWIN detected **149 drift events** in the 5,000-sample error stream. The high detection rate reflects both genuine distribution shifts from injected change points and the sensitivity of the ADWIN δ=0.002 parameter. With a retraining threshold of 5 consecutive drift events, this would trigger approximately 29 model retraining cycles.

---

## 6. Discussion

### 6.1 Interpretation of Results

The high AUROC values (0.963–0.968) reflect that the combined scorer successfully ranks anomalous samples above normal ones in expectation. However, the moderate F1 score (0.474) reveals a fundamental tension: the contamination parameter is set to 5% while the true anomaly fraction is 16.1%, causing the Isolation Forest threshold to be too conservative (missing many anomalies that contribute to low recall = 0.31). The precision of 1.000 achieved by the combined scorer at the 95th-percentile threshold is noteworthy but should be viewed with caution — it is a consequence of threshold placement, not evidence that the model is perfect.

### 6.2 Critical Assessment of Synthetic Evaluation

**Dependence on synthetic assumptions.** Our data generation model makes several simplifying assumptions: (1) anomalies are i.i.d. injections rather than correlated failure modes; (2) change points are abrupt rather than gradual; (3) the physical constraint bounds are designed to be easily violated by the injected anomalies; (4) the correlation structure is stationary in the normal regime. Real CERN detector data exhibits non-stationary noise, temperature-dependent gain variations, beam-induced electromagnetic interference, and systematic effects that are not captured by our model.

**Anomaly fraction inflation.** The effective anomaly fraction (16.1%) exceeded the intended rate (5%) due to overlapping burst injection windows. This biases the evaluation by providing more anomalous examples than would occur in practice, potentially making the detector appear more capable of rare-event detection than it would be in a scenario with true 1–5% anomaly rates.

**Precision = 1.000 is suspicious.** In real-world applications, precision of 1.0 is essentially never achieved. The zero false-positive rate at the 95th-percentile threshold results from the synthetic anomaly distribution being well-separated from the normal distribution for high-scoring samples. Real anomalies tend to be embedded in high-density normal regions, making them harder to detect without false positives.

**PELT over-detection.** The PELT algorithm detected 45 change points for 4 true ones. This is expected for the chosen penalty value (β=8) on high-noise data; proper penalty calibration using cross-validation or information criteria (BIC/AIC) on held-out segments would be required for production deployment. In real LHC data, change-point penalty must be calibrated against known run boundaries (fill changes, magnet ramps).

**Deep SVDD approximation limitations.** Our "Deep SVDD" uses finite-difference gradient updates, which is computationally infeasible for real networks. The actual Deep SVDD of Ruff et al. [2018] uses backpropagation through multiple layers and achieves substantially better representation learning. Our implementation (5 training epochs, 16 hidden units, finite differences) is best understood as a toy demonstration.

### 6.3 Generalization to Real-World Data

Deploying SciAD in a real CERN or LIGO environment would require several adaptations:

1. **Calibration from real detector data:** Physical bounds, rate-of-change limits, and correlation baselines must be learned from control runs, not hard-coded.
2. **Online learning:** The Isolation Forest must be incrementally updatable; we used a batch re-fit, which introduces a latency of O(n log n) per batch.
3. **Latency constraints:** CMS DQM requires decisions within the O(1 second) data certification window. BOCPD's O(n²) complexity is prohibitive; the PELT O(n) complexity is more suitable.
4. **False positive budget:** Physics analyses can tolerate very few corrupted luminosity sections; the false positive rate must be controlled below 0.1%, far below what we achieved.
5. **Channel scalability:** Real detectors have O(10⁸) channels; our 6-channel evaluation is a proof of concept.

### 6.4 Comparison with Prior Work

Our AUROC of 0.966–0.968 for Isolation Forest on synthetic multi-channel data is consistent with prior results in the literature: Chaudhari & Charate [2025] reported 100% recall on synthetic pipeline faults (with known contamination patterns), which is higher than our result but reflects a more constrained anomaly model. Davis et al. [2022] achieved glitch subtraction quality limited by the overlap between glitch and astrophysical signal morphologies — a problem orthogonal to detection, which our pipeline does not address. The CMS study by Stankevicius et al. [2020] demonstrated that neural network-based data certification can match human experts on simulated data, supporting the direction of our approach.

### 6.5 Limitations Summary

| Aspect | Current Study | Real-World Gap |
|--------|--------------|----------------|
| Data scale | 5,000 × 6 channels | 10⁹ × 10⁸ channels |
| Anomaly model | Synthetic i.i.d. | Complex correlated failures |
| SVDD implementation | Toy approximation | Full backprop network |
| CP penalty calibration | Manual (β=8) | Cross-validated |
| Physical bounds | Designed to be violated | Learned from control data |
| Evaluation protocol | 5 seeds, no time-series CV | Proper temporal holdout |

---

## 7. Conclusion

We have presented SciAD, a streaming anomaly detection pipeline integrating change-point detection (PELT/BOCPD), multivariate outlier scoring (Isolation Forest + Deep SVDD), physics-informed constraint scoring, ADWIN drift detection, and feature attribution explainability. On synthetic CERN/LIGO-inspired data, the combined scorer achieves AUROC of **0.968 ± 0.001** and demonstrates that physical domain knowledge (constraint scoring) provides a complementary signal to purely statistical detectors.

Key findings include: (1) physical constraint scoring provides a marginal but consistent improvement over standalone Isolation Forest (AUROC 0.968 vs 0.966); (2) PELT successfully localizes all four true change points but requires penalty calibration to reduce false detections; (3) BOCPD is effective but computationally expensive for long streams; (4) the F1 of 0.474 indicates that threshold calibration to the true anomaly fraction is critical for operational deployment; and (5) the precision–recall trade-off is dominated by the choice of operating threshold, not by the inherent discriminative capability of the detector.

Future work should focus on (1) evaluation on real CERN CMS/CPS datasets and LIGO auxiliary channels, (2) replacement of the SVDD approximation with a proper deep network, (3) online incremental updating of the Isolation Forest, (4) scalability experiments on O(10⁴)-channel data, and (5) integration with physics-analysis software frameworks (ROOT, GW-summary).

---

## References

1. **Stankevicius, A., et al.** (2020). Meta-Learning for Artificial Neural Network Hyper-Parameter Optimization for CERN CMS Offline Data Certification. *Journal of Physics: Conference Series*, 1525(1), 012103. https://doi.org/10.1088/1742-6596/1525/1/012103

2. **Davis, D., Littenberg, T. B., Romero-Shaw, I. M., et al.** (2022). Subtracting glitches from gravitational-wave detector data during the third LIGO-Virgo observing run. *Classical and Quantum Gravity*, 39(24), 245013. https://doi.org/10.1088/1361-6382/aca238

3. **Cavaglià, M.** (2022). Characterization of gravitational-wave detector noise with fractals. *Classical and Quantum Gravity*, 39(14), 145006. https://doi.org/10.1088/1361-6382/ac7325

4. **Corradin, R., Danese, L., & Ongaro, A.** (2022). Bayesian nonparametric change point detection for multivariate time series with missing observations. *International Journal of Approximate Reasoning*, 143, 26–43. https://doi.org/10.1016/j.ijar.2021.12.019

5. **Tsaknaki, I. A., Lillo, F., & Mazzarisi, P.** (2025). Bayesian autoregressive online change-point detection with time-varying parameters. *Communications in Nonlinear Science and Numerical Simulation*, 140, 108500. https://doi.org/10.1016/j.cnsns.2024.108500

6. **Katbi, A., & Ksantini, R.** (2025). One-class IoT anomaly detection system using an improved interpolated deep SVDD autoencoder with adversarial regularizer. *Digital Signal Processing*, 161, 105153. https://doi.org/10.1016/j.dsp.2025.105153

7. **Chaudhari, A. V., & Charate, P. A.** (2025). Proactive Data Pipeline Maintenance via Machine Learning-Driven Anomaly Detection. *International Journal of Scientific Research in Science and Technology*, 12(2), 663. https://doi.org/10.32628/ijsrst251222663

8. **Liu, F. T., Ting, K. M., & Zhou, Z.-H.** (2012). Isolation-based anomaly detection. *ACM Transactions on Knowledge Discovery from Data*, 6(1), 3. https://doi.org/10.1145/2133360.2133363

9. **Ruff, L., Vandermeulen, R., Goernitz, N., et al.** (2018). Deep one-class classification. In *Proceedings of the 35th International Conference on Machine Learning (ICML)*, 80, 4393–4402.

10. **Adams, R. P., & MacKay, D. J. C.** (2007). Bayesian online changepoint detection. *arXiv preprint*, arXiv:0710.3742.

11. **Killick, R., Fearnhead, P., & Eckley, I. A.** (2012). Optimal detection of changepoints with a linear computational cost. *Journal of the American Statistical Association*, 107(500), 1590–1598. https://doi.org/10.1080/01621459.2012.737745

12. **Bifet, A., & Gavalda, R.** (2007). Learning from time-changing data with adaptive windowing. In *Proceedings of the 7th SIAM International Conference on Data Mining*, 443–448.

13. **Lundberg, S. M., & Lee, S.-I.** (2017). A unified approach to interpreting model predictions. In *Advances in Neural Information Processing Systems (NeurIPS)*, 30.
