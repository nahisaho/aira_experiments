# Comparative Analysis Framework for Quantum Machine Learning Expressibility and Classical Baselines

## Abstract
We present a comparative research framework for quantum machine learning (QML) focused on six axes: (1) parameterized quantum circuit (PQC) expressibility and entanglement capability, (2) theoretical conditions for quantum kernel advantage, (3) the impact of angle, amplitude, and IQP data encodings, (4) dataset characteristics associated with potential quantum advantage, (5) barren plateaus and trainability, and (6) practical robustness under IBM-style noise proxies. We combined ToolUniverse-based literature search, NatureLM responses, and executable PennyLane/Qiskit benchmarks. Experimentally, deep entangled circuits achieved the highest expressibility (KL divergence to Haar: 0.0040) and entanglement capability (Meyer-Wallach: 0.9424 ± 0.0419), but literature and gradient analysis jointly indicate that such expressibility can correlate with degraded trainability. In kernel classification, the quantum kernel was competitive and slightly superior on the Circles dataset (0.960 ± 0.037 vs 0.950 ± 0.055 for RBF-SVM), but underperformed RBF-SVM on Linear and Moons. For encoding strategies, amplitude encoding outperformed IQP and angle encoding in the tested setting, yet all quantum encodings lagged behind the classical RBF baseline. Gradient variance decayed strongly with circuit depth, consistent with barren plateau theory. The noise experiment, implemented with a depolarizing-channel proxy, showed expectation-value degradation but also exposed a methodological limitation in the supplied classification setup. Overall, our results support a nuanced view: QML advantage is conditional on feature-map structure, data geometry, trainability, and noise, rather than guaranteed by high expressibility alone.

**Keywords:** quantum machine learning, parameterized quantum circuits, expressibility, quantum kernels, barren plateaus, data encoding, noise robustness

## 1. Introduction
Quantum machine learning promises access to feature representations and hypothesis classes that may be difficult for classical models to emulate efficiently. However, practical QML performance depends on multiple interacting factors: the expressibility of the ansatz, the induced kernel or feature map, the data encoding strategy, the geometry of the dataset, the trainability of variational parameters, and the effects of noise on NISQ hardware. These considerations motivate a comparative framework rather than single-metric evaluation.

This paper develops such a framework and applies it to executable benchmarks. Our goal is not to claim a universal quantum advantage, but to clarify when QML is theoretically and practically promising, and when classical baselines remain stronger.

## 2. Literature Search via ToolUniverse

### 2.1 Search Procedure
Following the user-specified workflow, we first used `tooluniverse-find_tools` to identify literature-search tools corresponding to Semantic Scholar, PubMed, and Crossref. We then queried these sources with the following topics:

- quantum machine learning expressibility parameterized quantum circuits
- quantum kernel methods advantage classification
- barren plateaus quantum neural networks trainability
- quantum data encoding angle amplitude IQP
- quantum advantage machine learning benchmark

### 2.2 Search Outcomes
- **Semantic Scholar:** repeated API failures on the prescribed queries (HTTP 400), followed by a rate-limit event (HTTP 429) on retry.
- **PubMed:** successful retrieval of multiple relevant papers.
- **Crossref:** successful retrieval of multiple relevant papers and DOI verification.

### 2.3 Selected Papers and Key Findings

| Topic | Paper | Authors | Year | DOI | Key finding |
|---|---|---|---:|---|---|
| PQC expressibility | *Expressibility and Entangling Capability of Parameterized Quantum Circuits for Hybrid Quantum-Classical Algorithms* | Sim, Johnson, Aspuru-Guzik | 2019 | 10.1002/qute.201900070 | Expressibility and entangling capability can be estimated statistically; connectivity and gate choices materially affect ansatz quality; expressibility saturates with depth. |
| Expressibility-trainability trade-off | *Connecting Ansatz Expressibility to Gradient Magnitudes and Barren Plateaus* | Holmes, Sharma, Cerezo, Coles | 2022 | 10.1103/PRXQuantum.3.010313 | Higher ansatz expressibility is linked to smaller gradients and can worsen trainability. |
| Encoding-dependent expressivity | *Effect of data encoding on the expressive power of variational quantum-machine-learning models* | Schuld, Sweke, Meyer | 2021 | 10.1103/PhysRevA.103.032430 | Data encoding, not just the variational block, determines the accessible function class of a QML model. |
| Quantum kernels | *Quantum Machine Learning in Feature Hilbert Spaces* | Schuld, Killoran | 2019 | 10.1103/PhysRevLett.122.040504 | Quantum models can be interpreted as kernel methods in feature Hilbert spaces. |
| Quantum-enhanced feature maps | *Supervised learning with quantum-enhanced feature spaces* | Havlíček et al. | 2019 | 10.1038/s41586-019-0980-2 | Quantum feature maps may provide classically hard-to-evaluate kernels and a route to near-term kernel-based QML. |
| Data-driven limits of advantage | *Power of data in quantum machine learning* | Huang et al. | 2021 | 10.1038/s41467-021-22539-9 | With enough data, classical learners can become competitive even on tasks motivated by quantum structure; advantage is data-dependent. |
| Information-theoretic limits | *Information-Theoretic Bounds on Quantum Advantage in Machine Learning* | Huang, Kueng, Preskill | 2021 | 10.1103/PhysRevLett.126.190505 | Quantum advantage requires stringent conditions on data distribution, learnability, and classical approximability. |
| Barren plateaus | *Barren plateaus in quantum neural network training landscapes* | McClean et al. | 2018 | 10.1038/s41467-018-07090-4 | Gradients can vanish exponentially with system size for broad classes of random or highly expressive circuits. |
| Higher-order trainability | *Higher order derivatives of quantum neural networks with barren plateaus* | Cerezo, Coles | 2021 | 10.1088/2058-9565/abf51a | Higher-order derivatives are also exponentially suppressed; second-order optimization does not automatically resolve barren plateaus. |
| Specific QNN trainability | *Trainability of Dissipative Perceptron-Based Quantum Neural Networks* | Sharma, Cerezo, Cincio, Coles | 2022 | 10.1103/PhysRevLett.128.180505 | Even specialized dissipative QNN architectures can exhibit barren plateaus, depending on depth and cost structure. |
| Robust encoding | *Robust data encodings for quantum classifiers* | LaRose, Coyle | 2020 | 10.1103/PhysRevA.102.032420 | Robustness of encoding strongly affects downstream classifier stability and utility. |
| Quantum SVM theory | *Quantum Support Vector Machine for Big Data Classification* | Rebentrost, Mohseni, Lloyd | 2014 | 10.1103/PhysRevLett.113.130503 | Presents a theoretical quantum speedup for SVM-like classification under strong algorithmic assumptions. |
| Practical benchmark | *On the Applicability of Quantum Machine Learning* | Raubitzek, Mallinger | 2023 | 10.3390/e25070992 | Quantum models can beat weak classical baselines, but strong classical learners often remain superior in practice. |

## 3. NatureLM Scientific Validation
We queried NatureLM with four fixed scientific prompts. No tool execution errors occurred, but the quality of the answers was mixed and, in two cases, materially inaccurate.

### 3.1 Raw NatureLM Outputs

**Q1.** “What is the expressibility of parameterized quantum circuits and how is it quantified using frame potential?”  
**Output:** “The expressibility of parameterized quantum circuits is the ability to describe a quantum circuit using a gate that depends on a parameter. The parameter can be any real number, complex number, or even a function of time. The expressibility of parameterized quantum circuits is quantified using frame potential.”

**Q2.** “What are the theoretical conditions for quantum kernel methods to show advantage over classical kernels?”  
**Output:** A generic description of quantum kernels and feature-space mapping, without explicit hardness or statistical conditions.

**Q3.** “What causes barren plateaus in quantum neural networks and how does it scale with system size?”  
**Output:** Included incorrect references to dropout and dropout regularization as mitigation mechanisms, suggesting conflation with classical deep learning.

**Q4.** “How does quantum data encoding strategy (angle encoding vs amplitude encoding vs IQP encoding) affect model expressibility?”  
**Output:** Claimed a ranking of encoding expressibility but without a precise theoretical basis or citation.

### 3.2 Assessment
NatureLM was useful only as a weak auxiliary summarizer. It did not reliably reproduce the accepted literature-level statements. Accordingly, all conclusions in this paper prioritize literature evidence and executable experiments over NatureLM prose.

## 4. Experimental Setup

### 4.1 Implementation
The complete benchmark script was saved to:

`/app/projects/cc17ec7d-6222-4328-b918-5ee4a1869ce0/workspace/experiments/qml_benchmark.py`

The script was executed in a local virtual environment with PennyLane 0.45.0 and the requested scientific Python stack.

### 4.2 Experiments
We ran six experiments:
1. Expressibility measurement via KL divergence to a Haar-inspired fidelity distribution.
2. Quantum kernel SVM vs classical linear/RBF SVM.
3. Encoding comparison: angle, IQP, amplitude.
4. Barren plateau analysis via gradient variance vs depth/width.
5. Noise robustness under depolarizing-channel simulation.
6. Entanglement capability via the Meyer-Wallach measure.

### 4.3 Datasets
Synthetic datasets included Linear, Moons, and Circles, all small enough for explicit kernel evaluation. This design supports controlled comparison but limits claims about large-scale generalization.

## 5. Results

### 5.1 PQC Expressibility and Entanglement Capability
| Ansatz | KL divergence | Meyer-Wallach mean ± SD |
|---|---:|---:|
| Shallow (L=1) | 0.5588 | 0.4926 ± 0.1665 |
| HEA (L=2) | 0.0264 | 0.8210 ± 0.1314 |
| Deep (L=3) | 0.0040 | 0.9424 ± 0.0419 |
| IQP | 0.0735 | 0.8037 ± 0.1434 |

Interpretation: the deeper entangled circuit was the most expressive and most entangling. This agrees with the literature: more expressive ansätze approach Haar-random statistics more closely. However, literature also warns that this same property often correlates with vanishing gradients.

![Expressibility analysis](figures/expressibility_analysis.png)

![Expressibility vs entanglement](figures/expressibility_vs_entanglement.png)

### 5.2 Quantum Kernels vs Classical Kernels
| Dataset | RBF-SVM | Linear-SVM | Quantum kernel SVM |
|---|---:|---:|---:|
| Linear | 0.850 ± 0.084 | 0.650 ± 0.130 | 0.790 ± 0.097 |
| Moons | 0.920 ± 0.068 | 0.860 ± 0.086 | 0.810 ± 0.058 |
| Circles | 0.950 ± 0.055 | 0.430 ± 0.081 | **0.960 ± 0.037** |

Interpretation: the quantum kernel was not uniformly superior. It slightly outperformed the classical RBF baseline on Circles, but underperformed on Linear and Moons. This supports the view that quantum kernel advantage is conditional on dataset geometry rather than generic.

![Kernel comparison](figures/kernel_comparison.png)

### 5.3 Data Encoding Comparison
| Encoding | Accuracy |
|---|---:|
| Angle | 0.354 ± 0.118 |
| IQP | 0.542 ± 0.078 |
| Amplitude | 0.583 ± 0.029 |
| Classical RBF | 0.900 ± 0.019 |

Interpretation: amplitude encoding produced the best quantum result in this setup, followed by IQP and angle encoding. Nonetheless, the classical RBF baseline remained substantially stronger. This result is consistent with the theoretical observation that encoding choice defines the accessible hypothesis class, but also highlights that better expressivity does not automatically translate into better trainability or accuracy.

![Encoding comparison](figures/encoding_comparison.png)

### 5.4 Barren Plateau Analysis
**Gradient variance vs depth (n = 4):**

| Depth | Variance |
|---:|---:|
| 1 | 5.30e-01 |
| 2 | 2.99e-01 |
| 4 | 8.85e-02 |
| 6 | 5.62e-02 |
| 8 | 4.02e-02 |
| 10 | 3.21e-02 |

**Gradient variance vs width (L = 3):**

| Width | Variance |
|---:|---:|
| 2 | 2.16e-01 |
| 3 | 1.67e-01 |
| 4 | 1.55e-01 |
| 5 | 1.60e-01 |
| 6 | 1.70e-01 |
| 7 | 1.52e-01 |
| 8 | 1.28e-01 |

Interpretation: depth scaling showed a clear decay in gradient variance, consistent with barren plateau onset. Width scaling was less clean, likely because the experiment remained in a small-qubit, local-cost, finite-sample regime. Thus, the numerical trend supports the theory qualitatively but does not fully reproduce the asymptotic exponential width dependence predicted in the literature.

![Barren plateau](figures/barren_plateau.png)

### 5.5 Noise Robustness Under IBM-Style Noise Proxy
| Depolarizing rate p | Accuracy |
|---:|---:|
| 0.000 | 0.400 |
| 0.005 | 0.400 |
| 0.010 | 0.400 |
| 0.020 | 0.400 |
| 0.050 | 0.400 |
| 0.100 | 0.400 |

Interpretation: the supplied benchmark uses a depolarizing-channel simulation rather than a calibrated IBM backend model. Moreover, in the classification loop, the noisy circuit does not actually encode the input sample `x`, so the accuracy stays flat at a poor baseline. Consequently, the most meaningful output of this section is the expectation-value degradation plot, not the classification accuracy. This is still useful as a simplified IBM-style noise proxy, but not as a faithful hardware benchmark.

![Noise analysis](figures/noise_analysis.png)

## 6. Discussion

### 6.1 Expressibility Is Necessary but Not Sufficient
Our experiments and the literature agree that deeper, more entangling circuits can better approximate Haar-like state distributions. However, the same circuits become harder to optimize. Thus, expressibility should be treated as one axis of model design, not the sole objective.

### 6.2 Conditions for Quantum Kernel Advantage
Theoretical work suggests that a quantum kernel can show advantage only when:
1. the induced feature map is classically hard to simulate or approximate,
2. the dataset aligns with that feature map,
3. the sample complexity remains favorable, and
4. classical competitors do not already capture the relevant geometry.

Our benchmark reflects exactly this conditionality: the quantum kernel helped on Circles, but not elsewhere.

### 6.3 Data Encoding as a First-Class Design Choice
Encoding changes the function class available to the model. Amplitude encoding was strongest among the tested quantum choices, but practical amplitude state preparation remains expensive. Angle encoding is hardware-friendly, while IQP-style encoding can induce richer nonlinear interactions. There is no universally best encoding; the right choice depends on data dimension, hardware constraints, and optimization behavior.

### 6.4 Dataset Characteristics for Quantum Advantage
Based on the literature and the benchmark, datasets most favorable to quantum models tend to exhibit:
- nonlinear decision structure,
- higher-order correlations,
- feature interactions naturally represented by quantum phase relationships,
- weak alignment with standard low-complexity classical kernels.

By contrast, when an RBF kernel already captures the geometry well, the margin for quantum advantage narrows considerably.

### 6.5 Practical NISQ Implications
In near-term hardware, noise and trainability jointly constrain QML. Even before hardware noise dominates, barren plateaus can make optimization impractical. Therefore, practical QML design should emphasize structured ansätze, shallow circuits, problem-informed feature maps, and backend-aware error mitigation.

## 7. Limitations
1. Semantic Scholar queries failed, so the literature synthesis relied mainly on Crossref and PubMed.
2. NatureLM outputs were partially inaccurate and cannot be treated as authoritative.
3. The noise experiment used a depolarizing proxy rather than calibrated IBM backend noise.
4. The width-scaling barren plateau experiment remained too small to show a strong asymptotic law.
5. The encoding benchmark used simple numerical training on small subsets and is best interpreted qualitatively.

## 8. Conclusion
This study supports a conditional, design-centric view of QML. Highly expressive circuits can generate powerful feature spaces, but excessive expressibility risks barren plateaus. Quantum kernels can be competitive, but only on suitably aligned datasets. Data encoding is not a preprocessing detail; it is a core determinant of model capacity. Under realistic NISQ constraints, the most promising direction is not “maximum expressibility,” but rather a careful co-design of ansatz structure, encoding, dataset geometry, and noise robustness.

## References
1. Sim, S., Johnson, P. D., & Aspuru-Guzik, A. (2019). *Expressibility and Entangling Capability of Parameterized Quantum Circuits for Hybrid Quantum-Classical Algorithms*. DOI: 10.1002/qute.201900070.
2. Holmes, Z., Sharma, K., Cerezo, M., & Coles, P. J. (2022). *Connecting Ansatz Expressibility to Gradient Magnitudes and Barren Plateaus*. DOI: 10.1103/PRXQuantum.3.010313.
3. Schuld, M., Sweke, R., & Meyer, J. J. (2021). *Effect of data encoding on the expressive power of variational quantum-machine-learning models*. DOI: 10.1103/PhysRevA.103.032430.
4. Schuld, M., & Killoran, N. (2019). *Quantum Machine Learning in Feature Hilbert Spaces*. DOI: 10.1103/PhysRevLett.122.040504.
5. Havlíček, V., et al. (2019). *Supervised learning with quantum-enhanced feature spaces*. DOI: 10.1038/s41586-019-0980-2.
6. Huang, H.-Y., Broughton, M., Mohseni, M., Babbush, R., Boixo, S., Neven, H., & McClean, J. R. (2021). *Power of data in quantum machine learning*. DOI: 10.1038/s41467-021-22539-9.
7. Huang, H.-Y., Kueng, R., & Preskill, J. (2021). *Information-Theoretic Bounds on Quantum Advantage in Machine Learning*. DOI: 10.1103/PhysRevLett.126.190505.
8. McClean, J. R., Boixo, S., Smelyanskiy, V. N., Babbush, R., & Neven, H. (2018). *Barren plateaus in quantum neural network training landscapes*. DOI: 10.1038/s41467-018-07090-4.
9. Cerezo, M., & Coles, P. J. (2021). *Higher order derivatives of quantum neural networks with barren plateaus*. DOI: 10.1088/2058-9565/abf51a.
10. Sharma, K., Cerezo, M., Cincio, L., & Coles, P. J. (2022). *Trainability of Dissipative Perceptron-Based Quantum Neural Networks*. DOI: 10.1103/PhysRevLett.128.180505.
11. LaRose, R., & Coyle, B. (2020). *Robust data encodings for quantum classifiers*. DOI: 10.1103/PhysRevA.102.032420.
12. Rebentrost, P., Mohseni, M., & Lloyd, S. (2014). *Quantum Support Vector Machine for Big Data Classification*. DOI: 10.1103/PhysRevLett.113.130503.
13. Raubitzek, S., & Mallinger, K. (2023). *On the Applicability of Quantum Machine Learning*. DOI: 10.3390/e25070992.
