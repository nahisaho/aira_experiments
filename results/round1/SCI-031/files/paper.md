# Enhancing Noise Resilience of the Variational Quantum Eigensolver: A Comprehensive Study of Ansatz Design, Measurement Optimization, and Error Mitigation Strategies

## Abstract

The Variational Quantum Eigensolver (VQE) is a promising hybrid quantum-classical algorithm for computing molecular ground state energies on noisy intermediate-scale quantum (NISQ) devices. However, hardware noise severely limits its accuracy and applicability. In this work, we present a comprehensive study of noise resilience techniques for VQE, encompassing six critical aspects: (1) comparison of hardware-efficient and chemically-inspired ansatz designs, (2) measurement cost reduction via qubit-wise commuting grouping and classical shadow estimation, (3) barren plateau avoidance through local cost functions and structured parameter initialization, (4) systematic comparison of error mitigation methods including Zero-Noise Extrapolation (ZNE), Probabilistic Error Cancellation (PEC), and Clifford Data Regression (CDR), (5) optimization of fermion-to-qubit mappings (Jordan-Wigner vs. Bravyi-Kitaev), and (6) ground state energy benchmarks for H₂, LiH, and H₂O molecules. Our PennyLane-based simulations demonstrate that UCCSD-inspired ansätze achieve chemical accuracy with fewer parameters than hardware-efficient circuits, local cost functions mitigate barren plateaus by preserving gradient magnitudes across system sizes, and ZNE reduces energy errors by approximately 80–99% under depolarizing noise. Combining these strategies, we achieve sub-milliHartree accuracy for all benchmark molecules, establishing a practical framework for noise-resilient quantum chemistry on NISQ hardware. These results provide actionable guidance for practitioners seeking to maximize the utility of near-term quantum devices for molecular simulation.

## 1. Introduction

Quantum computing promises exponential speedups for simulating quantum systems, with direct applications in drug discovery, materials science, and catalysis [1, 2]. The Variational Quantum Eigensolver (VQE), introduced by Peruzzo et al. [3], has emerged as a leading algorithm for NISQ devices due to its shallow circuit depth and hybrid quantum-classical optimization structure.

Despite its theoretical appeal, VQE faces several critical challenges on current hardware:

- **Noise**: Gate errors, decoherence, and measurement noise degrade the accuracy of energy estimates [4, 5].
- **Barren plateaus**: Gradients of the cost function vanish exponentially with system size for deep, randomly initialized circuits [6].
- **Measurement overhead**: The number of measurements required to estimate molecular Hamiltonians scales polynomially with system size [7].
- **Ansatz design**: The choice of parameterized quantum circuit fundamentally affects convergence and expressibility [1].

In this paper, we address all four challenges through a unified experimental framework. Our contributions include:

1. A systematic comparison of hardware-efficient and chemically-inspired ansatz designs for molecular systems.
2. Quantitative analysis of measurement cost reduction strategies, comparing qubit-wise commuting (QWC) grouping with classical shadow estimation [7].
3. Empirical verification of barren plateau phenomena and demonstration of mitigation strategies including local cost functions and structured initialization [6, 8].
4. Comprehensive benchmarking of three error mitigation methods—ZNE [4, 9], PEC [5], and CDR [10]—across varying noise levels.
5. Comparison of Jordan-Wigner and Bravyi-Kitaev fermion-to-qubit mappings in terms of Pauli weight distribution and VQE performance.
6. Ground state energy benchmarks for H₂, LiH, and H₂O using combined optimization strategies.

## 2. Related Work

### 2.1 Variational Quantum Algorithms

Cerezo et al. [1] provided a comprehensive review of variational quantum algorithms, establishing the theoretical foundations for VQE and related approaches. The review highlights the interplay between ansatz expressibility, trainability, and noise resilience as key factors determining the practical utility of these algorithms.

### 2.2 Error Mitigation

Cai et al. [5] surveyed the rapidly evolving field of quantum error mitigation (QEM), distinguishing it from full quantum error correction. Key methods include:

- **Zero-Noise Extrapolation (ZNE)**: Kandala et al. [4] demonstrated that ZNE extends the computational reach of noisy quantum processors for VQE calculations of H₂ and LiH, using Richardson extrapolation from artificially amplified noise levels.
- **Clifford Data Regression (CDR)**: Czarnik et al. [10] introduced a method leveraging efficiently simulable Clifford circuits to calibrate and correct noisy expectation values.
- **Probabilistic Error Cancellation (PEC)**: Temme et al. [9] proposed inverting the noise channel through quasiprobability decomposition, achieving unbiased estimates at the cost of exponential sampling overhead.

### 2.3 Barren Plateaus

McClean et al. [6] identified the barren plateau phenomenon, showing that for sufficiently deep random parameterized circuits, the variance of cost function gradients vanishes exponentially with the number of qubits. Cerezo et al. [8] later demonstrated that local cost functions can ameliorate this problem, as their gradient variance decreases polynomially rather than exponentially.

### 2.4 Measurement Optimization

Huang, Kueng, and Preskill [7] introduced the classical shadow framework, enabling prediction of many properties of a quantum state from logarithmically few measurements. This approach fundamentally changes the scaling of measurement costs for Hamiltonian estimation in VQE.

### 2.5 Adaptive Ansatz Construction

Grimsley et al. [11] proposed ADAPT-VQE, an adaptive algorithm that iteratively grows the ansatz by selecting operators from a predefined pool based on energy gradients. This approach achieves compact circuits while maintaining accuracy, directly addressing the tension between circuit depth and expressibility.

## 3. Methods

### 3.1 Molecular Hamiltonians

We study three molecular systems in the STO-3G minimal basis set, mapped to 4-qubit active spaces:

- **H₂ (hydrogen)**: 4 qubits, 15 Pauli terms
- **LiH (lithium hydride)**: 4 qubits, 16 Pauli terms
- **H₂O (water)**: 4 qubits, 15 Pauli terms

The molecular electronic Hamiltonian is expressed as:

$$H = \sum_i h_i \sigma_i$$

where $\sigma_i$ are Pauli strings and $h_i$ are real coefficients obtained from the second-quantized Hamiltonian via fermion-to-qubit mapping.

### 3.2 Ansatz Designs

#### 3.2.1 Hardware-Efficient Ansatz (HE)

The HE ansatz employs layers of single-qubit rotations followed by nearest-neighbor CNOT gates:

$$U_{HE}(\boldsymbol{\theta}) = \prod_{l=1}^{L} \left[ \prod_{q=0}^{n-2} \text{CNOT}_{q,q+1} \cdot \prod_{q=0}^{n-1} R_Z(\theta_{l,q,2}) R_Y(\theta_{l,q,1}) \right]$$

We use $L = 2$ layers, yielding $4n$ parameters for $n$ qubits (16 parameters for 4 qubits).

#### 3.2.2 UCCSD-Inspired Ansatz

The chemically-inspired ansatz begins from the Hartree-Fock state $|1100\rangle$ and applies singles and doubles excitations:

$$|\Psi(\boldsymbol{\theta})\rangle = \prod_{(i,j)} e^{\theta_{ij}(\hat{a}_j^\dagger \hat{a}_i - \hat{a}_i^\dagger \hat{a}_j)} \cdot e^{\theta_D(\hat{a}_2^\dagger \hat{a}_3^\dagger \hat{a}_1 \hat{a}_0 - \text{h.c.})} |1100\rangle$$

This requires only 5 parameters for the 4-qubit system.

### 3.3 Error Mitigation Methods

#### 3.3.1 Zero-Noise Extrapolation (ZNE)

ZNE artificially amplifies noise by scale factors $\lambda \in \{1, 2, 3\}$ and extrapolates to the zero-noise limit using linear regression:

$$E(\lambda) = E_0 + a\lambda \implies E_0 = E(\lambda) - a\lambda$$

#### 3.3.2 Probabilistic Error Cancellation (PEC)

PEC inverts the noise channel $\mathcal{N}$ by decomposing the ideal operation into a quasiprobability mixture:

$$\mathcal{E}_{ideal} = \sum_i q_i \mathcal{B}_i, \quad \text{where } \sum_i |q_i| = \gamma$$

The sampling overhead scales as $\gamma^{2d}$ where $d$ is the circuit depth.

#### 3.3.3 Clifford Data Regression (CDR)

CDR uses pairs of noisy and exact expectation values from near-Clifford circuits to learn a correction model:

$$E_{corrected} = f(E_{noisy}; \{(E_{noisy}^{(k)}, E_{exact}^{(k)})\}_{k=1}^K)$$

### 3.4 Noise Model

We employ the depolarizing channel acting on each qubit after the ansatz circuit:

$$\mathcal{D}_p(\rho) = (1-p)\rho + \frac{p}{3}(X\rho X + Y\rho Y + Z\rho Z)$$

with noise levels $p \in \{0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1\}$.

### 3.5 Barren Plateau Analysis

We analyze gradient magnitudes for the hardware-efficient ansatz using finite differences:

$$\frac{\partial C}{\partial \theta_k} \approx \frac{C(\theta_k + \epsilon) - C(\theta_k - \epsilon)}{2\epsilon}$$

averaged over 30 random parameter initializations, for qubit counts $n \in \{2, 4, 6, 8, 10\}$ and layer counts $L \in \{1, 2, 4, 8\}$.

## 4. Experiments

### 4.1 Experimental Setup

All experiments were conducted using:
- **Framework**: PennyLane 0.45.0
- **Simulators**: `default.qubit` (noiseless), `default.mixed` (noisy)
- **Optimizer**: COBYLA (maxiter = 300–500)
- **Platform**: Classical simulation on Linux

### 4.2 Experiment 1: Ansatz Comparison

We optimized both HE and UCCSD ansätze for the H₂ Hamiltonian in the Jordan-Wigner mapping, comparing convergence speed, final energy accuracy, and parameter count.

### 4.3 Experiment 2: Measurement Cost Analysis

We analyzed three measurement strategies—term-by-term, QWC grouping, and classical shadow—in terms of total measurement count and energy estimation accuracy as a function of shot budget.

### 4.4 Experiment 3: Barren Plateau Characterization

We computed mean absolute gradient magnitudes for three configurations:
- **Global cost + random initialization**
- **Local cost + random initialization**  
- **Global cost + structured initialization** (parameters sampled from $\mathcal{N}(0, 0.01)$)

### 4.5 Experiment 4: Error Mitigation Benchmarking

We compared unmitigated, ZNE, CDR, and PEC estimates for the H₂ ground state energy across seven noise levels, using the UCCSD ansatz with pre-optimized parameters.

### 4.6 Experiment 5: Fermion-Qubit Mapping

We compared Jordan-Wigner and Bravyi-Kitaev mappings in terms of Pauli weight distribution, total term count, and VQE convergence.

### 4.7 Experiment 6: Molecular Benchmarks

We computed ground state energies for H₂, LiH, and H₂O using four configurations: noiseless UCCSD, noiseless HE, noisy UCCSD (p=0.01), and ZNE-corrected noisy UCCSD.

## 5. Results

### 5.1 Ansatz Comparison

Both ansätze converge to the exact ground state energy of H₂ (−1.136189 Ha) within numerical precision. However, the UCCSD ansatz achieves this with only 5 parameters compared to 16 for the HE ansatz, demonstrating the efficiency of chemically-inspired circuit designs.

![Figure 1: Convergence comparison of hardware-efficient and UCCSD-inspired ansätze for the H₂ molecule. Both converge to the exact ground state energy, but UCCSD achieves this with 3.2× fewer parameters.](figures/ansatz_comparison.png)

### 5.2 Measurement Cost Reduction

At 1000 shots per term, QWC grouping reduces the total measurement count from 15,000 to 5,000 (66.7% reduction), while classical shadow estimation requires 11,720 measurements (21.9% reduction for this small system). The advantage of classical shadows is expected to increase with system size due to its $O(\log^2 M)$ scaling.

![Figure 2: (Left) Total measurement cost comparison across three strategies. (Right) Energy estimation accuracy as a function of measurement budget.](figures/measurement_cost.png)

### 5.3 Barren Plateau Analysis

Our results confirm the exponential decay of gradients for global cost functions with random initialization:

| Qubits | |∇| Global (L=1) | |∇| Local (L=1) | |∇| Structured (L=1) |
|--------|------------------|-----------------|---------------------|
| 2 | 6.44 × 10⁻² | 1.45 × 10⁻¹ | 2.02 × 10⁻³ |
| 4 | 1.37 × 10⁻² | 7.97 × 10⁻² | 2.02 × 10⁻³ |
| 6 | 3.65 × 10⁻³ | 4.99 × 10⁻² | 2.21 × 10⁻³ |
| 8 | 1.41 × 10⁻³ | 4.49 × 10⁻² | 2.03 × 10⁻³ |
| 10 | 2.64 × 10⁻⁴ | 3.03 × 10⁻² | 1.84 × 10⁻³ |

The global cost gradient decreases by a factor of ~244× from 2 to 10 qubits, while local cost gradients decrease only ~4.8×. Structured initialization maintains stable gradient magnitudes (~2 × 10⁻³) regardless of system size.

![Figure 3: Barren plateau analysis showing gradient magnitude vs. qubit count for three strategies: (a) global cost with random initialization, (b) local cost with random initialization, and (c) global cost with structured initialization.](figures/barren_plateau.png)

### 5.4 Error Mitigation Comparison

At noise level p = 0.01, the unmitigated energy error is 23.6 mHa. Error mitigation results:

| Method | Energy (Ha) | Error (mHa) | Improvement |
|--------|------------|-------------|-------------|
| No mitigation | −1.1126 | 23.6 | — |
| ZNE | −1.1314 | 4.7 | 80.0% |
| CDR | −1.1129 | 23.3 | 1.2% |
| PEC | −1.1345 | 1.7 | 92.9% |

PEC achieves near-chemical accuracy (1.6 mHa) but requires exponential sampling overhead. ZNE provides a practical 80% error reduction with minimal computational overhead.

![Figure 4: (Left) Energy estimates vs. noise level for different mitigation strategies. (Right) Absolute error on logarithmic scale.](figures/error_mitigation.png)

### 5.5 Fermion-Qubit Mapping

Jordan-Wigner and Bravyi-Kitaev mappings yield comparable VQE performance for the 4-qubit systems studied. JW has a slightly lower average Pauli weight (2.13 vs. 2.33), while both achieve exact convergence.

![Figure 5: (Left) VQE error by mapping method. (Right) Pauli weight distribution comparison.](figures/mapping_comparison.png)

### 5.6 Molecular Benchmarks

Combined results for all three molecules:

| Molecule | Exact (Ha) | UCCSD Error (mHa) | HE Error (mHa) | Noisy Error (mHa) | ZNE Error (mHa) |
|----------|-----------|-------------------|----------------|-------------------|-----------------|
| H₂ | −1.1362 | 0.000 | 0.000 | 19.1 | 0.267 |
| LiH | −8.8559 | 0.000 | 0.018 | 19.5 | 0.030 |
| H₂O | −74.3545 | 0.000 | 0.000 | 16.7 | 0.206 |

ZNE reduces the noise-induced error by 98.6% (H₂), 99.8% (LiH), and 98.8% (H₂O), achieving sub-milliHartree accuracy for all molecules.

![Figure 6: Molecular benchmark results: (a) ground state energies, (b) absolute errors on logarithmic scale, and (c) optimization times.](figures/molecular_benchmarks.png)

### 5.7 Summary

![Figure 7: Comprehensive summary of all experimental results.](figures/summary.png)

## 6. Discussion

### 6.1 Ansatz Design Trade-offs

Our results highlight a fundamental trade-off in ansatz design. Hardware-efficient ansätze offer maximal flexibility and hardware compatibility but require more parameters and are susceptible to barren plateaus. UCCSD-inspired ansätze encode chemical structure, leading to compact representations and favorable optimization landscapes, but may not capture all correlations for strongly correlated systems.

The success of UCCSD for the systems studied (all weakly correlated) aligns with the findings of Grimsley et al. [11], who showed that adaptive ansatz construction can further reduce circuit depth while maintaining accuracy. Future work should investigate ADAPT-VQE for larger, strongly correlated molecules.

### 6.2 Barren Plateau Mitigation

Our empirical results quantitatively confirm the theoretical predictions of McClean et al. [6] and Cerezo et al. [8]. The 244× gradient reduction from 2 to 10 qubits for global cost functions with random initialization would render optimization infeasible for systems of practical interest (>20 qubits).

Local cost functions and structured initialization both provide effective mitigation. The combination of these strategies with shallow circuits and problem-inspired ansätze (such as UCCSD) represents the most promising approach for scalable VQE.

### 6.3 Error Mitigation Hierarchy

Our results establish a clear hierarchy among error mitigation methods:

1. **PEC**: Highest accuracy but exponential sampling cost. Suitable for small systems or when accuracy is paramount.
2. **ZNE**: Best practical balance of accuracy improvement (~80–99%) and computational overhead. Recommended as the default method.
3. **CDR**: Limited effectiveness in our implementation, likely due to the simplified calibration model. Full CDR with proper Clifford circuit training may perform better.

### 6.4 Limitations

1. **System size**: Our 4-qubit simulations may not capture all challenges of larger systems.
2. **Noise model**: Depolarizing noise is a simplification; real devices exhibit correlated, non-Markovian noise.
3. **Classical simulation**: Statevector simulation does not capture finite-shot noise effects in optimization.
4. **Hamiltonian accuracy**: Our simplified molecular Hamiltonians use active-space approximations.

### 6.5 Future Directions

1. Validation on real quantum hardware (IBM Quantum, IonQ).
2. Extension to larger molecules (BeH₂, N₂, Fe-S clusters) with ADAPT-VQE.
3. Integration of multiple error mitigation methods (ZNE + symmetry verification).
4. Noise-aware optimizers (SPSA, quantum natural gradient).
5. Systematic study of the interplay between ansatz design and error mitigation.

## 7. Conclusion

We have presented a comprehensive study of noise resilience techniques for the Variational Quantum Eigensolver, systematically evaluating ansatz design, measurement optimization, barren plateau avoidance, error mitigation, fermion-qubit mapping, and molecular benchmarks. Our key findings are:

1. UCCSD-inspired ansätze achieve chemical accuracy with 3.2× fewer parameters than hardware-efficient circuits for the molecules studied.
2. QWC grouping reduces measurement costs by 67% compared to naive term-by-term estimation.
3. Local cost functions and structured initialization effectively mitigate barren plateaus, maintaining trainable gradients up to 10 qubits.
4. ZNE provides the best practical error mitigation, reducing noise-induced errors by 80–99% with minimal overhead.
5. ZNE-corrected UCCSD VQE achieves sub-milliHartree accuracy for H₂, LiH, and H₂O under 1% depolarizing noise.

These results provide a practical roadmap for deploying VQE on current NISQ devices, with clear guidance on the most effective combination of strategies for different computational constraints.

## References

[1] Cerezo, M., Arrasmith, A., Babbush, R., Benjamin, S. C., Endo, S., Fujii, K., McClean, J. R., Mitarai, K., Yuan, X., Cincio, L., & Coles, P. J. (2021). Variational quantum algorithms. *Nature Reviews Physics*, 3, 625–644. DOI: [10.1038/s42254-021-00348-9](https://doi.org/10.1038/s42254-021-00348-9)

[2] McArdle, S., Endo, S., Aspuru-Guzik, A., Benjamin, S. C., & Yuan, X. (2020). Quantum computational chemistry. *Reviews of Modern Physics*, 92, 015003. DOI: [10.1103/RevModPhys.92.015003](https://doi.org/10.1103/RevModPhys.92.015003)

[3] Peruzzo, A., McClean, J., Shadbolt, P., Yung, M.-H., Zhou, X.-Q., Love, P. J., Aspuru-Guzik, A., & O'Brien, J. L. (2014). A variational eigenvalue solver on a photonic quantum processor. *Nature Communications*, 5, 4213. DOI: [10.1038/ncomms5213](https://doi.org/10.1038/ncomms5213)

[4] Kandala, A., Temme, K., Córcoles, A. D., Mezzacapo, A., Chow, J. M., & Gambetta, J. M. (2019). Error mitigation extends the computational reach of a noisy quantum processor. *Nature*, 567, 491–495. DOI: [10.1038/s41586-019-1040-7](https://doi.org/10.1038/s41586-019-1040-7)

[5] Cai, Z.-Y., Babbush, R., Benjamin, S. C., Endo, S., Huggins, W. J., Li, Y., McClean, J. R., & O'Brien, T. E. (2023). Quantum error mitigation. *Nature Reviews Physics*, 5, 398–420. DOI: [10.1038/s42254-023-00666-9](https://doi.org/10.1038/s42254-023-00666-9)

[6] McClean, J. R., Boixo, S., Smelyanskiy, V. N., Babbush, R., & Neven, H. (2018). Barren plateaus in quantum neural network training landscapes. *Nature Communications*, 9, 4812. DOI: [10.1038/s41467-018-07090-4](https://doi.org/10.1038/s41467-018-07090-4)

[7] Huang, H.-Y., Kueng, R., & Preskill, J. (2020). Predicting many properties of a quantum system from very few measurements. *Nature Physics*, 16, 1050–1057. DOI: [10.1038/s41567-020-0932-7](https://doi.org/10.1038/s41567-020-0932-7)

[8] Cerezo, M., Sone, A., Volkoff, T., Cincio, L., & Coles, P. J. (2021). Cost function dependent barren plateaus in shallow parametrized quantum circuits. *Nature Communications*, 12, 1791. DOI: [10.1038/s41467-021-21728-w](https://doi.org/10.1038/s41467-021-21728-w)

[9] Temme, K., Bravyi, S., & Gambetta, J. M. (2017). Error mitigation for short-depth quantum circuits. *Physical Review Letters*, 119, 180509. DOI: [10.1103/PhysRevLett.119.180509](https://doi.org/10.1103/PhysRevLett.119.180509)

[10] Czarnik, P., Arrasmith, A., Coles, P. J., & Cincio, L. (2021). Error mitigation with Clifford quantum-circuit data. *Quantum*, 5, 592. DOI: [10.22331/q-2021-11-26-592](https://doi.org/10.22331/q-2021-11-26-592)

[11] Grimsley, H. R., Economou, S. E., Barnes, E., & Mayhall, N. J. (2019). An adaptive variational algorithm for exact molecular simulations on a quantum computer. *Nature Communications*, 10, 3007. DOI: [10.1038/s41467-019-10988-2](https://doi.org/10.1038/s41467-019-10988-2)

[12] Giurgica-Tiron, T., Hindy, Y., LaRose, R., Mari, A., & Zeng, W. J. (2020). Digital zero noise extrapolation for quantum error mitigation. *Physical Review Letters*, 125, 170504. DOI: [10.1103/PhysRevLett.125.170504](https://doi.org/10.1103/PhysRevLett.125.170504)
