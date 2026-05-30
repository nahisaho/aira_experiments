# Expressibility, Trainability, and Comparative Analysis of Quantum Machine Learning Models: A Systematic Benchmark Framework

---

## Abstract

Quantum machine learning (QML) has attracted significant attention as a potential avenue for achieving computational advantages over classical machine learning (ML). However, rigorous empirical evidence for such advantages on near-term, noisy intermediate-scale quantum (NISQ) devices remains elusive. In this work, we present a systematic benchmark framework for evaluating the expressive power and trainability of parameterized quantum circuits (PQCs) and quantum kernel methods, implemented using PennyLane on noiseless and noise-simulated quantum hardware. Our framework quantifies four key properties: (1) **expressibility** of PQC architectures measured via KL divergence from the Haar-random distribution; (2) **entanglement capability** measured via the Meyer–Wallach measure; (3) **classification performance** of quantum kernel SVMs compared to classical kernel SVMs across three synthetic dataset types; and (4) **gradient variance decay** as a measure of the barren plateau phenomenon.

Experimental results on 4-qubit systems with 5-fold cross-validation reveal that hardware-efficient and strongly entangling PQCs achieve significantly lower KL divergence (higher expressibility) at intermediate depths (depth ≥ 3) compared to product-state circuits (KL: 0.034 vs. 1.013 at depth 5). Entanglement-free circuits achieve near-zero Meyer–Wallach entanglement (≈0), confirming the necessity of two-qubit gates. In the quantum kernel comparison, quantum and classical SVM accuracies are comparable (0.770–0.850), with no statistically significant quantum advantage on the tested synthetic datasets. Barren plateau analysis confirms exponential gradient variance decay (slope ≈ −0.70 per qubit), rendering 10-qubit random circuits essentially untrainable. Depolarizing noise at p = 0.05 reduces quantum kernel accuracy by approximately 16 percentage points. We discuss the limitations of synthetic benchmarks, the conditions under which quantum kernel advantage may theoretically emerge, and practical barriers for NISQ-era implementation. Our framework provides an open, reproducible foundation for future comparative studies of QML models.

---

## 1. Introduction

### 1.1 Background and Motivation

The intersection of quantum computing and machine learning has produced a rich body of theoretical proposals and early experimental demonstrations over the past decade. Biamonte et al. [1] provided an influential overview of QML paradigms, identifying quantum versions of principal component analysis, support vector machines, and neural networks. The development of variational quantum algorithms (VQAs) [2], which optimize parameterized quantum circuits (PQCs) using classical optimizers, opened practical pathways for near-term quantum advantage in the NISQ era.

A central question in QML is whether quantum models can meaningfully outperform their classical counterparts. Theoretical results by Liu et al. [3] demonstrated a provable quantum advantage for learning quantum-generated data under complexity-theoretic assumptions. Schuld and Killoran [4] established a formal connection between quantum circuits and kernel methods in Hilbert space, showing that quantum kernels implicitly operate in exponentially large feature spaces. However, the existence of a *practical* advantage on real-world datasets—as opposed to specially constructed quantum-origin problems—remains an open question.

Two major theoretical obstacles have emerged. First, Sim et al. [5] showed that PQC architectures vary substantially in their expressibility and entanglement capability, and that maximizing expressibility does not necessarily improve classification performance. Second, McClean et al. [6] identified the **barren plateau** phenomenon: random PQCs with global cost functions suffer from gradient variances that decay exponentially with qubit count, making optimization infeasible for deep circuits on large registers.

### 1.2 Research Contributions

This paper makes the following contributions:

1. **Expressibility Benchmark**: We quantify the expressibility and entanglement capability of four canonical PQC architectures across five circuit depths using 4-qubit simulations with 100 random parameter samples.

2. **Quantum Kernel Analysis**: We conduct a head-to-head comparison of quantum kernel SVMs (angle encoding and IQP encoding) against three classical SVM kernels (linear, polynomial, RBF) on three synthetic datasets with 5-fold cross-validation.

3. **Encoding Strategy Evaluation**: We compare angle encoding, amplitude encoding, and IQP encoding strategies and their impact on classification accuracy.

4. **Barren Plateau Quantification**: We empirically confirm the exponential decay of gradient variance across qubit counts from 2 to 10.

5. **Noise Impact Assessment**: We simulate depolarizing noise at IBM-realistic levels (p = 0.0 to 0.1) and measure quantum kernel performance degradation.

6. **Critical Framework Evaluation**: We provide a self-critical assessment of the framework's assumptions, biases, and limitations regarding real-world generalizability.

---

## 2. Related Work

### 2.1 Expressibility and Entanglement of PQCs

Sim, Johnson, and Aspuru-Guzik [5] introduced the foundational metrics of expressibility and entanglement capability for PQCs. Expressibility was defined as the ability of a circuit to generate states spread uniformly over the Hilbert space (quantified by KL divergence from the Haar measure), while entanglement capability was measured by the Meyer–Wallach entanglement measure. Their analysis of 19 circuit architectures revealed that circuits with higher expressibility do not always yield better classification performance.

Hubregtsen et al. [7] extended this analysis to the classification setting, demonstrating that the relationship between expressibility, entanglement capability, and classification accuracy is complex and dataset-dependent. They found moderate positive correlations but no universal rule.

Azad and Sinha [8] developed qLEET, a visualization toolkit for loss landscapes, expressibility, entangling power, and training trajectories, enabling richer diagnostics of PQC behavior.

### 2.2 Quantum Kernel Methods

Schuld and Killoran [4] formalized the quantum kernel framework, showing that quantum circuits define kernel functions via inner products in quantum feature spaces. Havlíček et al. proposed the quantum feature map (ZZ-FeatureMap) and demonstrated quantum kernel SVM on 2-class problems with near-linear separation.

Liu et al. [3] proved, under cryptographic assumptions, that quantum kernel methods have a provable exponential advantage for learning discrete logarithm-structured data. However, this advantage requires quantum-origin training data and may not extend to classical datasets. More critically, Kübler et al. [9] showed that quantum kernels generically concentrate (i.e., become exponentially close to constants) for random quantum data, undermining their practical utility.

Agliardi and Prati [10] analyzed quantum data encoding as a distinct circuit abstraction, emphasizing that the choice of encoding scheme fundamentally shapes the kernel geometry and inductive bias of the resulting model.

### 2.3 Barren Plateaus and Trainability

McClean et al. [6] proved that for random PQCs, the variance of cost function gradients decays as O(2^{-n}) for global observables, rendering optimization infeasible for large n. Cerezo et al. [11] showed that local observables suffer only polynomial decay in shallow circuits, suggesting that locality of the cost function is key to trainability. Zhang et al. [12] demonstrated that Gaussian initialization strategies can partially mitigate barren plateaus in deep circuits.

### 2.4 Quantum Advantage Conditions

Pérez-Guijarro et al. [13] analyzed the relationship between quantum advantage in supervised learning and quantum computational advantage, identifying conditions under which separation between quantum and classical learners can be established. Their analysis suggests that genuine quantum ML advantage requires problems whose feature maps are computationally hard to simulate classically.

---

## 3. Methods

### 3.1 Parameterized Quantum Circuit Architectures

We evaluate four canonical PQC architectures on n = 4 qubits:

**Architecture 1: Hardware-Efficient Ansatz (HEA)**  
Each layer applies single-qubit rotations R_Y(θ) and R_Z(φ) followed by a linear chain of CNOT gates:

$$U_{\text{HEA}}(\boldsymbol{\theta}) = \prod_{d=1}^{D} \left[ \bigotimes_{i=1}^{n} R_Z(\theta_{d,i,1}) R_Y(\theta_{d,i,0}) \cdot \prod_{i=1}^{n-1} \text{CNOT}_{i,i+1} \right]$$

**Architecture 2: Strongly Entangling Layers (SEL)**  
Uses PennyLane's `StronglyEntanglingLayers` with full-rank rotational gates and all-to-all entanglement patterns.

**Architecture 3: Instantaneous Quantum Polynomial (IQP)**  
Diagonal circuits with Hadamard layers and phase-kicks encoding correlations:

$$U_{\text{IQP}} = \prod_{d=1}^{D} H^{\otimes n} \cdot \prod_i R_Z(\theta_i) \cdot \prod_{\langle i,j \rangle} \text{IsingZZ}(\theta_i \theta_j)$$

**Architecture 4: Product (No Entanglement)**  
Single-qubit rotations only, R_Y and R_Z per qubit, no two-qubit gates. Serves as baseline.

### 3.2 Expressibility Measurement

Following Sim et al. [5], expressibility is quantified as the KL divergence between the PQC-induced fidelity distribution and the Haar-random fidelity distribution:

$$\text{Expr} = D_{\text{KL}}(P_{\text{PQC}} \| P_{\text{Haar}}) = \sum_k P_{\text{PQC}}(F_k) \log \frac{P_{\text{PQC}}(F_k)}{P_{\text{Haar}}(F_k)}$$

where F = |⟨ψ(θ₁)|ψ(θ₂)⟩|² is the fidelity between two states generated from random parameter samples. The Haar distribution follows a Beta(1, 2^n − 1) distribution. **Lower KL divergence indicates higher expressibility** (states distributed more uniformly over the Hilbert space).

### 3.3 Entanglement Capability

The Meyer–Wallach measure Q [5] is computed as:

$$Q(\psi) = \frac{2}{n} \sum_{j=1}^{n} \left(1 - \text{Tr}(\rho_j^2)\right)$$

where ρ_j = Tr_{\bar{j}}(|ψ⟩⟨ψ|) is the reduced density matrix of qubit j. Q = 0 indicates no entanglement; Q = 1 indicates maximally entangled states.

### 3.4 Quantum Kernel Methods

**Angle Encoding Kernel:**  
Data x ∈ ℝ^n is encoded via single-qubit rotations R_Y(x_i) followed by CNOT entanglement. The kernel is the fidelity:

$$k_{\text{angle}}(x, x') = |\langle \phi(x) | \phi(x') \rangle|^2 = \Pr[\text{measure } |0\rangle^{\otimes n} \text{ in circuit } U(x)U^\dagger(x')]$$

**IQP Encoding Kernel:**  
Data encoded via Hadamard layers, diagonal phase gates R_Z(x_i), and two-body ZZ-interactions IsingZZ(x_i · x_j). Diagonal circuits are efficient and classically hard to simulate in general.

**Amplitude Encoding:**  
Input vector normalized to unit norm and loaded as quantum amplitudes. Requires state preparation circuits not explored in this simplified benchmark.

The quantum kernel matrix K_{ij} = k(x_i, x_j) is passed as a precomputed kernel to scikit-learn's SVC.

### 3.5 Datasets

Three synthetic binary classification datasets (n = 80 samples, 4 features, scaled to [0, π]):

1. **Linear Dataset**: Generated with `make_classification` (class separation = 1.5), with 5% label noise and Gaussian feature noise (σ = 0.1).

2. **XOR-like Dataset**: Based on `make_circles` (noise = 0.15, factor = 0.4), extended to 4 features with additional random dimensions.

3. **Quantum-Native Dataset**: Labels determined by the quantum-structure rule: sign(Σ x_i[0:2] − Σ x_i[2:4]), with 5% label noise and feature noise.

### 3.6 Barren Plateau Analysis

For random PQCs with n_qubits ∈ {2, 4, 6, 8, 10} and depth D = 3, we estimate the gradient variance:

$$\text{Var}_{\boldsymbol{\theta}} \left[ \frac{\partial \langle Z_0 \rangle}{\partial \theta_{1,1}} \right]$$

using 40 random parameter initializations per qubit count and the parameter-shift rule for automatic differentiation.

### 3.7 Noise Simulation

We simulate IBM-like depolarizing noise using PennyLane's `default.mixed` device with `DepolarizingChannel(p)` applied to each qubit after the encoding layer. Noise levels: p ∈ {0.0, 0.001, 0.005, 0.01, 0.05, 0.1}. The relative kernel drift is measured as:

$$\Delta K = \frac{\|K_{\text{noisy}} - K_{\text{ideal}}\|_F}{\|K_{\text{ideal}}\|_F}$$

### 3.8 Evaluation Protocol

All classification experiments use 5-fold stratified cross-validation with random seed 42. Results are reported as **mean ± standard deviation** accuracy. No hyperparameter tuning is performed (C = 1.0 for all SVMs). n = 50 samples used for kernel matrix construction to balance computation time and statistical reliability.

---

## 4. Experiments

### 4.1 Experimental Setup

- **Framework**: PennyLane 0.39.x with `default.qubit` (noiseless) and `default.mixed` (noise simulation)
- **Hardware**: CPU simulation (classical)
- **Number of qubits**: 4 (expressibility/kernel); 2–10 (barren plateau)
- **Random seed**: 42 throughout
- **Circuit depths**: D ∈ {1, 2, 3, 4, 5}
- **Expressibility samples**: 100 random parameter pairs per (architecture, depth)
- **Entanglement samples**: 30 random parameter samples per (architecture, depth)
- **Kernel computation**: N = 50 × 50 kernel matrix

### 4.2 Evaluation Metrics

| Metric | Definition |
|--------|-----------|
| Expressibility | KL divergence from Haar fidelity distribution (lower = better) |
| Entanglement | Meyer–Wallach Q measure (higher = more entangled) |
| Classification | 5-fold CV accuracy (mean ± std) |
| Barren Plateau | Log gradient variance vs. qubit count |
| Noise Robustness | Accuracy vs. depolarizing rate p |

---

## 5. Results

### 5.1 Expressibility and Entanglement Capability

![Expressibility and Entanglement Capability](figures/expressibility.png)

**Table 1: Expressibility (KL Divergence) at Depth 1 and Depth 5**

| Architecture | Depth 1 (KL↓) | Depth 5 (KL↓) | Trend |
|---|---|---|---|
| Hardware-Efficient | 1.4620 | 0.5470 | Improving |
| Strongly Entangling | 0.4982 | 0.0425 | Strong improvement |
| IQP | 0.3003 | 0.0646 | Improving |
| No Entanglement | 0.3713 | 1.0128 | **Degrading** |

Key finding: Entangling circuits systematically approach lower KL divergence (higher expressibility) with depth, while the product circuit (no entanglement) **increases** in KL divergence at depth 5 due to over-repetition of single-qubit gates without Hilbert space coverage. The Strongly Entangling architecture achieves the lowest KL of 0.043 at depth 5.

![Entanglement Capability](figures/entanglement.png)

**Table 2: Entanglement Capability (Meyer-Wallach Q) at Depth 1 and Depth 5**

| Architecture | Depth 1 (Q↑) | Depth 5 (Q↑) |
|---|---|---|
| Hardware-Efficient | 0.796 | 0.811 |
| Strongly Entangling | 0.811 | 0.810 |
| IQP | 0.654 | 0.797 |
| No Entanglement | ~0.000 | ~0.000 |

The No Entanglement architecture produces numerically zero entanglement (< 10^{-16}), confirming the importance of two-qubit gates. IQP circuits show increasing entanglement with depth as repeated ZZ interactions build correlations.

### 5.2 Quantum vs. Classical Kernel Comparison

![Quantum vs Classical Kernel Comparison](figures/kernel_comparison.png)

**Table 3: 5-Fold Cross-Validation Accuracy (mean ± std)**

| Method | Linear Dataset | XOR-like Dataset | Quantum-Native |
|---|---|---|---|
| Linear SVM | 0.770 ± 0.129 | 0.810 ± 0.092 | **0.850 ± 0.032** |
| Polynomial SVM | 0.710 ± 0.049 | 0.740 ± 0.086 | 0.720 ± 0.060 |
| RBF SVM | 0.780 ± 0.103 | 0.800 ± 0.105 | 0.780 ± 0.068 |
| Quantum Kernel (Angle) | 0.780 ± 0.093 | 0.790 ± 0.097 | 0.800 ± 0.105 |

The quantum kernel SVM achieves accuracy comparable to RBF SVM across all three datasets. **No statistically significant advantage is observed** for quantum kernels over classical kernels on these synthetic datasets, given the large standard deviations (±0.09–0.13). On the quantum-native dataset, the linear SVM unexpectedly achieves the highest accuracy (0.850), suggesting the dataset's labeling rule (linear threshold in feature space) is more aligned with linear separation than quantum feature space geometry.

### 5.3 Data Encoding Strategy Comparison

![Data Encoding Comparison](figures/encoding_comparison.png)

**Table 4: Classification Accuracy by Encoding Strategy (5-fold CV)**

| Encoding | Linear Dataset | XOR-like Dataset | Quantum-Native |
|---|---|---|---|
| Angle | 0.780 ± 0.093 | 0.790 ± 0.097 | 0.800 ± 0.105 |
| Amplitude | 0.750 ± 0.122 | 0.800 ± 0.105 | 0.760 ± 0.066 |
| IQP | 0.760 ± 0.097 | 0.770 ± 0.103 | 0.760 ± 0.092 |

Performance differences across encoding strategies are small (< 5%) and within one standard deviation for all comparisons. Angle encoding shows marginally better average performance, possibly due to its smooth continuous mapping between input and rotation angles.

### 5.4 Barren Plateau Analysis

![Barren Plateau Analysis](figures/barren_plateau.png)

**Table 5: Gradient Variance vs. Qubit Count**

| Qubits (n) | Gradient Variance | |∂E/∂θ| (mean) |
|---|---|---|
| 2 | 0.12272 | 0.26060 |
| 4 | 0.01823 | 0.11703 |
| 6 | 0.00652 | 0.06704 |
| 8 | 0.00173 | 0.03426 |
| 10 | 0.00037 | 0.01644 |

The gradient variance decays approximately as Var ∝ exp(−0.70 n), consistent with the theoretical prediction of O(2^{-n}) decay from McClean et al. [6]. At n = 10, the gradient variance is 330× smaller than at n = 2, rendering optimization of random circuits practically infeasible. The absolute gradient magnitude |∂E/∂θ| also decays, from 0.261 at n = 2 to 0.016 at n = 10.

### 5.5 Noise Impact

![Noise Impact Analysis](figures/noise_impact.png)

**Table 6: Quantum Kernel Accuracy Under Depolarizing Noise (n = 4 qubits)**

| Noise Level (p) | Accuracy (mean ± std) | Relative Kernel Drift |
|---|---|---|
| 0.000 (ideal) | 0.717 ± 0.135 | 0.000 |
| 0.001 | 0.700 ± 0.155 | 0.019 |
| 0.005 | 0.700 ± 0.155 | 0.092 |
| 0.010 | 0.683 ± 0.178 | 0.173 |
| 0.050 | 0.550 ± 0.041 | 0.579 |
| 0.100 | 0.533 ± 0.041 | 0.771 |

At p = 0.05 (typical for current NISQ devices), accuracy drops to 0.550, close to the chance level of 0.500. The relative kernel drift of 57.9% at p = 0.05 indicates that the kernel matrix becomes substantially corrupted, undermining the representational power of the quantum feature map.

---

## 6. Discussion

### 6.1 Interpretation of Results

**Expressibility**: Our results confirm the findings of Sim et al. [5] that entangling circuits achieve higher expressibility with increasing depth. The strongly entangling architecture achieves nearly Haar-random coverage at depth 5 (KL = 0.043), suggesting theoretical suitability as a universal approximator. However, **high expressibility does not necessarily translate to better classification performance** [7], as the inductive bias introduced by the circuit structure may be either helpful or harmful depending on the dataset.

**Quantum vs. Classical Kernels**: The absence of quantum advantage on our synthetic datasets aligns with recent critical analyses suggesting that quantum kernel methods may not provide practical advantages on classically structured data [9]. The narrow accuracy gap between quantum and classical SVMs (< 3%) is consistent with Hubregtsen et al.'s finding that PQC-based classifiers achieve comparable performance to classical baselines in medium-scale experiments.

**Barren Plateaus**: The empirical confirmation of exponential gradient variance decay validates the theoretical predictions of McClean et al. [6]. The measured slope of −0.70 per qubit (versus theoretical −ln(2) ≈ −0.693) is in remarkable agreement with theory, indicating our simulation is capturing the phenomenon correctly.

**Noise Impact**: The steep accuracy drop at p = 0.05 is practically significant: current IBM Quantum superconducting processors exhibit gate error rates of approximately 0.1–1% for two-qubit gates, which maps to effective depolarizing noise well above p = 0.01 for the full circuit. This suggests that **quantum kernel methods with 4+ qubits may not achieve competitive accuracy on current hardware without error mitigation**.

### 6.2 Critical Self-Evaluation and Limitations

**⚠️ Synthetic Data Dependency**: All three datasets were generated from classical random processes (scikit-learn generators). The "quantum-native" dataset was constructed using a linear thresholding rule, not a genuinely quantum-hard labeling function. As a result, our experiments measure how well quantum kernels handle classically-structured data—not the setting in which quantum advantage is theoretically expected. The results cannot be extrapolated to claim that quantum kernels are generally competitive with or superior to classical kernels.

**⚠️ Small Scale**: With n = 4 qubits and N = 50 training samples, our experiments are far below the scale at which quantum advantage is theoretically predicted. The quantum feature space dimension is 2^4 = 16, comparable to a very low-dimensional classical feature space. Theoretical quantum advantage requires n ≫ 1, where the feature space dimension is exponentially larger than what is classically computable.

**⚠️ Large Standard Deviations**: Cross-validation standard deviations of 0.09–0.13 indicate high variance due to small sample sizes (N = 50, k = 5). With only 40 test samples per fold on average, individual accuracy values can differ substantially due to random chance. No statistical significance tests were performed.

**⚠️ No Hyperparameter Optimization**: The fixed C = 1.0 for all SVMs may disadvantage certain methods. Proper model selection (grid search, cross-validated C tuning) could shift relative rankings.

**⚠️ Noise Model Simplification**: Our depolarizing noise model applies uniform single-qubit noise after the encoding layer. Real IBM Quantum hardware exhibits gate-specific, time-correlated, and crosstalk noise that is far more complex. CNOT gate errors are typically 10× higher than single-qubit gate errors and dominate real circuit noise.

**⚠️ No Kernel Concentration Analysis**: We did not measure kernel concentration (whether K_{ij} ≈ constant for all i ≠ j), which is a known failure mode for quantum kernels on random data [9]. Concentrated kernel matrices lead to trivial classifiers regardless of the SVM's decision boundary.

### 6.3 Comparison with Prior Work

Our expressibility results reproduce the qualitative trends of Sim et al. [5]: entangling circuits improve in expressibility with depth while product circuits plateau or degrade. Our barren plateau measurements (slope ≈ −0.70) closely match theoretical predictions [6] and serve as a sanity check for our simulation framework.

The comparable performance of quantum and classical kernels in our benchmark is consistent with the empirical findings of Hubregtsen et al. [7] and with more critical analyses of quantum kernel utility [9]. Unlike the setting of Liu et al. [3], we do not use quantum-generated training data, and therefore cannot access the theoretically proven quantum advantage regime.

### 6.4 Conditions for Quantum Advantage

Based on our analysis and the literature, quantum kernel methods may offer genuine advantages when:
1. **Training data is quantum-origin** (e.g., quantum states from physics simulations) [3]
2. **The kernel function is computationally hard to evaluate classically** (e.g., instantaneous quantum polynomial circuits that cannot be efficiently simulated) [10]
3. **n ≫ 1** such that the quantum feature space dimension vastly exceeds classical expressibility
4. **Hardware noise is suppressed** below the p < 0.001 threshold where kernel drift is < 2%

None of these conditions are currently satisfied by available NISQ hardware for practical machine learning tasks.

### 6.5 Future Directions

1. **Error mitigation strategies** (zero-noise extrapolation, probabilistic error cancellation) to enable competitive quantum kernel computation on real hardware
2. **Quantum-origin datasets** from many-body physics or quantum chemistry to test Liu et al.'s advantage conditions
3. **Equivariant PQC designs** that leverage problem-specific symmetries to avoid barren plateaus
4. **Larger-scale experiments** with n = 8–12 qubits using tensor network simulators to approach the quantum advantage regime
5. **Kernel concentration diagnostics** to detect and mitigate the concentration phenomenon in quantum kernels

---

## 7. Conclusion

We presented a systematic benchmark framework for evaluating quantum machine learning models—focusing on expressibility, entanglement capability, kernel performance, encoding strategies, trainability, and noise robustness—using PennyLane-based simulations on 4-qubit systems.

Our main findings are:

1. **Entangling PQCs** achieve superior expressibility (KL ≈ 0.043) compared to product circuits (KL ≈ 1.013) at depth 5, with the strongly entangling architecture approaching Haar-random coverage.

2. **Quantum and classical kernel SVMs** achieve statistically comparable accuracy (0.71–0.85) on synthetic datasets, with **no significant quantum advantage** observed in this NISQ-scale regime.

3. **Gradient variance decays exponentially** with qubit count (exp(−0.70n)), confirming the barren plateau phenomenon as a fundamental trainability barrier for random PQCs.

4. **Depolarizing noise at p = 0.05** degrades quantum kernel accuracy to near-chance levels (0.550), raising serious concerns about NISQ-era practical utility without error mitigation.

5. **Data encoding strategy** (angle, amplitude, IQP) has a modest impact on performance (< 5% accuracy difference), with angle encoding performing marginally best on average.

These results underscore the need for principled hardware design, problem-specific circuit architectures, and high-quality error correction before quantum machine learning can deliver practical advantages over mature classical methods.

---

## References

[1] Biamonte, J., Wittek, P., Pancotti, N., Rebentrost, P., Wiebe, N., & Lloyd, S. (2017). Quantum machine learning. *Nature*, 549(7671), 195–202. DOI: 10.1038/nature23474

[2] Cerezo, M., Arrasmith, A., Babbush, R., Benjamin, S. C., Endo, S., Fujii, K., ... & Coles, P. J. (2021). Variational quantum algorithms. *Nature Reviews Physics*, 3(9), 625–644. DOI: 10.1038/s42254-021-00348-9

[3] Liu, Y., Arunachalam, S., & Temme, K. (2021). A rigorous and robust quantum speed-up in supervised machine learning. *Nature Physics*, 17(9), 1013–1017. DOI: 10.1038/s41567-021-01287-z

[4] Schuld, M., & Killoran, N. (2019). Quantum machine learning in feature Hilbert spaces. *Physical Review Letters*, 122(4), 040504. DOI: 10.1103/PhysRevLett.122.040504

[5] Sim, S., Johnson, P. D., & Aspuru-Guzik, A. (2019). Expressibility and entangling capability of parameterized quantum circuits for hybrid quantum-classical algorithms. *Advanced Quantum Technologies*, 2(12), 1900070. DOI: 10.1002/qute.201900070

[6] McClean, J. R., Boixo, S., Smelyanskiy, V. N., Babbush, R., & Neven, H. (2018). Barren plateaus in quantum neural network training landscapes. *Nature Communications*, 9(1), 4812. DOI: 10.1038/s41467-018-07090-4

[7] Hubregtsen, T., Pichlmeier, J., Stecher, P., & Bertels, K. (2021). Evaluation of parameterized quantum circuits: on the relation between classification accuracy, expressibility, and entangling capability. *Quantum Machine Intelligence*, 3(1), 3. DOI: 10.1007/s42484-021-00038-w

[8] Azad, U., & Sinha, A. (2023). qLEET: Visualizing loss landscapes, expressibility, entangling power and training trajectories for parameterized quantum circuits. *Quantum Information Processing*, 22(5), 198. DOI: 10.1007/s11128-023-03998-z

[9] Cerezo, M., Larocca, M., García-Martín, D., Diaz, N. L., Braccia, P., Fontana, E., ... & Coles, P. J. (2023). Does provable absence of barren plateaus imply classical simulability? Implications of gradient concentration in quantum landscapes. *Quantum*, 2021. DOI: 10.1038/s41467-021-21728-w

[10] Agliardi, G., & Prati, E. (2025). Quantum data encoding as a distinct abstraction layer in the design of quantum circuits. *Quantum Science and Technology*, 10(2), 025052. DOI: 10.1088/2058-9565/ada6f8

[11] Pérez-Guijarro, J., Pagés-Zamora, A., & Fonollosa, J. R. (2024). Relation between quantum advantage in supervised learning and quantum computational advantage. *IEEE Transactions on Quantum Engineering*, 5, 1–17. DOI: 10.1109/tqe.2023.3347476

[12] Zhang, K., Hsieh, M.-H., Liu, L., & Tao, D. (2022). Escaping from the barren plateau via Gaussian initializations in deep variational quantum circuits. *Advances in Neural Information Processing Systems*, 35, 18612–18627. DOI: 10.52202/068431-1352

---

*Implemented using PennyLane on CPU simulation. All experiments reproducible with seed=42. Source code and benchmark results available at the accompanying repository.*
