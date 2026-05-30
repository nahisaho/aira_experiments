# Real-Time EEG Signal Processing and Decoding for Non-Invasive Brain–Computer Interfaces: An Integrated Pipeline Combining Artifact Removal, CSP, Deep Learning, and Online Adaptation

---

## Abstract

Non-invasive brain–computer interfaces (BCIs) based on electroencephalography (EEG) provide a critical communication pathway for individuals with severe motor impairments, including those with locked-in syndrome. Despite decades of progress, real-time BCI deployment remains challenged by poor signal quality, high inter-subject variability, non-stationarity of neural signals, and the computational cost of online processing. This paper presents an integrated real-time EEG processing pipeline combining: (1) artifact removal via Artifact Subspace Reconstruction (ASR) and Independent Component Analysis (ICA)-inspired denoising; (2) Common Spatial Pattern (CSP) feature extraction paired with Support Vector Machine (SVM) classification for motor imagery; (3) compact convolutional neural networks (EEGNet) and a Convolutional Transformer architecture (EEG Conformer) for end-to-end decoding; (4) an Euclidean Alignment transfer learning approach for cross-subject P300 speller adaptation; and (5) an incremental online learner with sliding-window adaptation for concept drift mitigation. We evaluate the pipeline on synthetic EEG datasets designed to emulate BCI Competition IV Dataset 2a conditions (22 channels, 250 Hz, 4-class motor imagery) and P300 speller paradigms. Five-fold cross-validation yields: CSP+SVM: 0.813±0.102 accuracy (κ=0.752±0.135); EEGNet: 0.570±0.164 accuracy (κ=0.427±0.218); EEG Conformer: 0.531±0.036 accuracy (κ=0.374±0.048). Transfer learning for P300 classification achieves 0.845±0.017 accuracy (vs. 0.970±0.015 without alignment, suggesting alignment sensitivity to domain gap). Adaptive online learning recovers post-drift accuracy to 0.769 versus 0.303 for a static classifier. We critically discuss the limitations of synthetic evaluation, the risk of optimistic bias, and the challenges of real-world deployment. This work provides a modular reference architecture for EEG-based assistive communication systems with particular relevance to locked-in syndrome applications.

**Keywords:** brain–computer interface, EEG, motor imagery, P300 speller, CSP, EEGNet, EEG Conformer, transfer learning, online learning, concept drift, locked-in syndrome

---

## 1. Introduction

Brain–computer interfaces that bypass the peripheral motor system and provide direct communication from neural activity to external devices hold transformative potential for individuals with severe neurological impairments [1]. For patients with amyotrophic lateral sclerosis (ALS), brainstem stroke, or late-stage locked-in syndrome (LIS), EEG-based BCIs may represent the sole viable communication channel [2]. The non-invasive nature of EEG—with millisecond temporal resolution, moderate spatial resolution, and no surgical risk—makes it the most clinically accessible BCI modality [3].

Despite significant advances over the past two decades, translating laboratory EEG-BCI systems into reliable real-world assistive devices remains an open challenge [4]. Key obstacles include:

- **Signal quality degradation**: Scalp EEG is contaminated by electromyographic (EMG) artifacts, ocular artifacts (EOG), and environmental electrical noise, often exceeding the neural signal amplitude by an order of magnitude.
- **Non-stationarity**: EEG statistics shift across sessions due to electrode impedance changes, fatigue, medication, and evolving neural plasticity—a phenomenon termed *concept drift*.
- **Inter-subject variability**: Spatial patterns of motor imagery ERD/ERS differ substantially across individuals, limiting the generalizability of trained models.
- **Real-time latency constraints**: Assistive communication applications require end-to-end processing latency below 500 ms to be responsive.
- **Limited labeled data**: Calibration sessions are cognitively and physically demanding for patients with motor impairments.

The emergence of deep learning—particularly compact CNN architectures such as EEGNet [5] and transformer-based models such as EEG Conformer [6]—has opened new avenues for automatic feature learning from raw EEG. Transfer learning approaches, including Euclidean Alignment (EA) [7] and domain adaptation, address the subject-variability problem by leveraging knowledge from source subjects to accelerate target-subject calibration.

This paper makes the following contributions:

1. **Integrated real-time pipeline**: We present a modular, end-to-end pipeline from raw EEG acquisition through online classification, covering preprocessing, feature extraction, classification, and adaptation.
2. **Comparative evaluation**: We systematically compare classical (CSP+SVM), deep learning (EEGNet), and transformer (EEG Conformer) approaches under identical experimental conditions with proper cross-validation.
3. **Transfer learning for P300**: We implement and evaluate Euclidean Alignment for cross-subject P300 speller adaptation with limited target-domain labels.
4. **Online adaptation**: We demonstrate that incremental adaptive classifiers significantly outperform static models under concept drift conditions relevant to long-duration BCI operation.
5. **Critical analysis**: We explicitly address the limitations of synthetic evaluation and the challenges of real-world deployment.

---

## 2. Related Work

### 2.1 EEG-Based BCI Paradigms

Rashid et al. [1] provide a comprehensive review of EEG-based BCI systems, identifying motor imagery, P300 visual evoked potentials, and SSVEP as the three dominant paradigms. Motor imagery BCIs exploit event-related desynchronization (ERD) in the mu (8–12 Hz) and beta (18–26 Hz) bands, which exhibit characteristic contra-lateral lateralization during imagined limb movements. P300 BCIs exploit the late positive deflection (~300 ms post-stimulus) elicited by rare target stimuli in an oddball paradigm, making them particularly suitable for spelling applications.

### 2.2 Signal Processing and Feature Extraction

Common Spatial Pattern (CSP) filtering remains the most widely used spatial filtering method for motor imagery BCI [3]. CSP finds a linear projection maximizing variance ratio between two classes, producing discriminative log-variance features. Extensions to multi-class settings use one-versus-rest (OVR) or one-versus-one (OVO) strategies. Saha and Baumert [8] highlight that intra-session and inter-session variability cause covariate shift in CSP feature distributions, motivating adaptive and transfer learning approaches.

### 2.3 Deep Learning for EEG Decoding

Lawhern et al. [5] introduced EEGNet, a compact depthwise separable CNN achieving cross-paradigm generalization with as few as 2.7K parameters. Gu et al. [3] survey deep learning methods for EEG-BCI, including CNNs, RNNs, and hybrid architectures. Nallani and Ramachandran [9] propose RLEEGNet, integrating deep Q-networks with 1D-CNN-LSTM for adaptive motor imagery classification, though the reported 100% accuracy on BCI Competition IV datasets raises concerns about data leakage.

### 2.4 Transformer Architectures for EEG

Song et al. [6] proposed EEG Conformer, combining temporal-spatial convolution with a self-attention module to capture both local and global temporal dependencies. The architecture achieves state-of-the-art performance on BCI Competition IV 2a (accuracy: 0.7758) and is now a widely used benchmark. Pfeffer et al. [10] review transformer-based models for EEG processing, noting their particular strength in capturing long-range temporal dependencies but highlighting the challenge of training with limited EEG data. Keutayeva et al. [11] demonstrate that compact convolutional transformers (EEGCCT) achieve 70.12% accuracy on leave-one-subject-out evaluation of BCI IV 2a, outperforming EEGNet under strict subject-independence conditions.

### 2.5 Transfer Learning and Adaptation

Cross-subject transfer is essential for practical deployment. Euclidean Alignment whitens EEG epochs by the subject-specific covariance matrix, projecting data toward a common distribution [8]. Saha and Baumert [8] discuss how both psychological and neurophysiological predictors of inter-subject variability may augment transfer learning. Lee et al. [12] propose continual learning of transformer models, initializing with action observation EEG and progressively adapting through motor imagery feedback, achieving improvement from 0.63 to 0.77 over three training sessions.

### 2.6 Gaps in Prior Work

Existing studies share several limitations this work seeks to address: (1) most deep learning results come from offline post-hoc evaluation, not real-time streaming scenarios; (2) concept drift during extended BCI operation is rarely studied in online settings; (3) few studies integrate artifact removal, classification, and adaptation into a unified testable pipeline; (4) evaluation on synthetic data with realistic noise characteristics is underutilized as a controlled benchmark.

---

## 3. Methods

### 3.1 Synthetic EEG Data Generation

We generated synthetic EEG emulating BCI Competition IV Dataset 2a conditions: 22 channels, 250 Hz sampling rate, 4.0-second trials, 4 motor imagery classes (left hand, right hand, feet, tongue), 288 trials per recording (72 per class), balanced design.

**Signal model**: For each trial of class $c$, the EEG signal $\mathbf{X} \in \mathbb{R}^{C \times T}$ was generated as:

$$\mathbf{X} = \underbrace{\mathbf{a}_c \cdot s(t)^{\top}}_{\text{neural}} + \underbrace{\mathbf{f} \cdot b(t)^{\top}}_{\text{blink}} + \underbrace{\mathbf{e}_{k} \cdot m(t)^{\top}}_{\text{EMG}} + \underbrace{\boldsymbol{\eta}(t)}_{\text{pink noise}}$$

where $\mathbf{a}_c \in \mathbb{R}^C$ is the class-specific topographic pattern, $s(t) = A_\mu \sin(2\pi f_\mu t + \phi) + A_\beta \sin(2\pi f_\beta t + \psi)$ is the combined mu+beta source signal with per-trial amplitude variability ($\sigma = 0.3$), $b(t)$ is a Gaussian blink transient ($\sim$ 100 μV) at a random time in [0.5, 1.5] s, $m(t)$ is an exponentially decaying EMG burst on a random channel, and $\boldsymbol{\eta}(t)$ is $1/f$ pink noise (noise_level = 3.0, SNR ≈ 3 dB). Class-specific topographic patterns were handcrafted to reflect known ERD lateralization: left-hand MI produces C4-dominant activation; right-hand produces C3-dominant; feet produces central midline; tongue produces frontal.

### 3.2 Preprocessing Pipeline

**Bandpass filtering**: Zero-phase 4th-order Butterworth filter (4–40 Hz), implemented via `scipy.signal.filtfilt`.

**Artifact Subspace Reconstruction (ASR)**: We implement a simplified ASR-like procedure that clips channel amplitudes exceeding $\theta \cdot \sigma_{\text{robust}}$ where $\theta = 4.5$ and $\sigma_{\text{robust}}$ is the median channel standard deviation per trial. A Hann window is applied to reduce edge artifacts.

**ICA-inspired denoising**: PCA whitening is applied to the multi-trial data matrix; the top 2 principal components (which capture frontal ocular variance) are attenuated by a factor of 20 before back-projection, simulating ICA-based artifact removal.

**Normalization**: Each trial is z-score normalized per channel.

### 3.3 CSP + SVM for Motor Imagery

**Multi-class CSP (OVR)**: For 4-class classification, we train one binary CSP per class using one-versus-rest decomposition. Each binary CSP solves:

$$\mathbf{W}^* = \arg\max_{\mathbf{W}} \frac{\mathbf{W}^{\top} \mathbf{C}_c \mathbf{W}}{\mathbf{W}^{\top} (\mathbf{C}_c + \mathbf{C}_{\text{rest}}) \mathbf{W}}$$

via generalized eigendecomposition of $(\mathbf{C}_c, \mathbf{C}_c + \mathbf{C}_{\text{rest}})$ with Tikhonov regularization ($\lambda = 10^{-4}$). The 4 most discriminative spatial filters (2 per eigenvalue extreme) are retained per class, yielding 32 spatial filters total.

**Feature extraction**: Log-variance of spatially filtered signals:

$$f_k = \log \frac{\text{var}(\mathbf{w}_k^{\top} \mathbf{X})}{\sum_j \text{var}(\mathbf{w}_j^{\top} \mathbf{X}) + \epsilon}$$

**Classifier**: RBF-kernel SVM (C=1.0, γ='scale').

### 3.4 EEGNet

We implement EEGNet [5] with the following architecture:
- **Temporal convolution**: Conv2D(1→F1=8, kernel=(1,64), padding=(0,32)) + BN
- **Spatial (depthwise) convolution**: DepthwiseConv2D(F1, D=2, kernel=(C,1)) + BN + ELU + AvgPool(1,4) + Dropout(0.4)
- **Separable convolution**: SeparableConv2D(F1·D→F2=16, kernel=(1,16)) + BN + ELU + AvgPool(1,8) + Dropout(0.4)
- **Classifier**: Dense(F2_flat → 4)

Total parameters: ~5.8K. Trained with Adam (lr=5×10⁻⁴, weight decay=10⁻⁴), cosine annealing LR schedule, batch size 32, 60 epochs.

### 3.5 EEG Conformer

We implement a simplified EEG Conformer [6] with:
- **Patch embedding**: Conv2D(1→emb=32, kernel=(1,40), stride=(1,40)) + BN + ELU, followed by spatial convolution Conv2D(emb, emb, (C,1))
- **Transformer encoder**: 2 layers, 4 attention heads, FFN dim=128, Dropout(0.4)
- **CLS token** aggregation
- **Classifier**: LayerNorm + Dense(32→4)

Total parameters: ~45K. Same training protocol as EEGNet.

### 3.6 P300 Speller with Transfer Learning

**P300 data simulation**: Binary classification (target vs. non-target), 600 source trials + 400 target trials. Target ERP: N200 component (−2 μV, 200 ms) + P300 component (+6 μV, 300 ms, Gaussian width 80 ms) on central-parietal channels, with per-trial amplitude variability (σ=0.3) and background pink noise. Target:non-target ratio = 1:5 (matching real P300 speller).

**Euclidean Alignment (EA)**: The source covariance matrix $\mathbf{R}_{\text{src}}$ is estimated from source-domain data. Target trials are aligned via:

$$\tilde{\mathbf{X}}_i = \mathbf{R}_{\text{src}}^{-1/2} \mathbf{X}_i$$

using eigendecomposition: $\mathbf{R}_{\text{src}}^{-1/2} = \mathbf{V} \boldsymbol{\Lambda}^{-1/2} \mathbf{V}^{\top}$.

**Classifier**: Linear Discriminant Analysis (LDA) on downsampled (4×) flattened EEG features, preceded by StandardScaler. Only 50 labeled target trials are used for training (low-data regime simulation).

### 3.7 Online Learning with Concept Drift

**Drift simulation**: After 50% of the stream, a gradual linear interpolation between the identity and a random orthogonal rotation is applied to CSP features:

$$\mathbf{x}'_i = [(1-\alpha_i)\mathbf{I} + \alpha_i \mathbf{Q}]\,\mathbf{x}_i, \quad \alpha_i = 0.4 \cdot \frac{i - i_{\text{drift}}}{N - i_{\text{drift}}}$$

where $\mathbf{Q}$ is a random orthogonal matrix and $\alpha_i$ ramps from 0 to 0.4.

**Adaptive classifier**: Sliding window LDA (window=50, update frequency=10 trials). Drift detection uses a simple heuristic: if rolling accuracy drops by >15% over 20 consecutive trials, a drift event is flagged and the buffer is reset.

**Comparison**: Static SVM trained on calibration data (first 10% of stream) is compared to the adaptive learner.

### 3.8 Evaluation

All classification experiments use 5-fold stratified cross-validation. We report:
- **Accuracy** (Acc)
- **Macro-averaged F1** (F1)
- **Cohen's κ** (chance-corrected agreement)
- **Standard deviation** across folds

Baseline chance accuracy is 0.25 for 4-class MI and 0.167 for P300 (1:5 target ratio).

---

## 4. Experiments

### 4.1 Dataset

Synthetic EEG: N=288 trials, 22 channels, 250 Hz, 4.0 s duration (1,000 time points per trial), 4 balanced classes. Pink noise SNR ≈ 3 dB (deliberately challenging). P300 synthetic data: 600 source + 400 target binary trials (0.8 s epochs). All data generated with fixed random seeds for reproducibility.

### 4.2 Experimental Conditions

| Experiment | Paradigm | N_trials | N_classes | Evaluation |
|------------|----------|----------|-----------|------------|
| A1: CSP+SVM | Motor imagery | 288 | 4 | 5-fold CV |
| A2: EEGNet | Motor imagery | 288 | 4 | 5-fold CV |
| A3: EEG Conformer | Motor imagery | 288 | 4 | 5-fold CV |
| B1: P300 w/o TL | P300 speller | 400 | 2 | 5-fold CV |
| B2: P300 w/ EA-TL | P300 speller | 400+600 | 2 | 5-fold CV |
| C: Online learning | Motor imagery | 259 (stream) | 4 | Rolling acc. |

### 4.3 Figures

**Figure 1** shows the effect of the preprocessing pipeline on a representative EEG trial (channel 13, class: left-hand MI). Raw EEG contains high-amplitude blink artifacts and EMG bursts; bandpass filtering suppresses out-of-band noise; ASR+ICA+normalization yields a signal with substantially reduced artifact amplitude while preserving neural oscillatory structure.

![Figure 1: Preprocessing Pipeline](figures/fig1_preprocessing.png)

**Figure 2** shows the first four CSP spatial filters learned from the 4-class OVR decomposition. Components 1 and 2 show contralateral activation patterns consistent with known motor cortex lateralization; components 3 and 4 reflect foot-related central midline and tongue-related frontal patterns.

![Figure 2: CSP Spatial Filters](figures/fig2_csp_patterns.png)

---

## 5. Results

### 5.1 Motor Imagery Classification

Table 1 summarizes the 5-fold cross-validation results for all three classifiers.

**Table 1: Motor Imagery Classification Results (5-Fold CV, 4-class, Synthetic EEG)**

| Model | Accuracy (mean±SD) | Macro F1 (mean±SD) | Cohen's κ (mean±SD) | Params |
|-------|-------------------|-------------------|---------------------|--------|
| Chance (1/4) | 0.250 | 0.250 | 0.000 | — |
| CSP + SVM | **0.813 ± 0.102** | **0.763 ± 0.136** | **0.752 ± 0.135** | — |
| EEGNet | 0.570 ± 0.164 | 0.553 ± 0.163 | 0.427 ± 0.218 | ~5.8K |
| EEG Conformer | 0.531 ± 0.036 | 0.528 ± 0.037 | 0.374 ± 0.048 | ~45K |

![Figure 3: Model Comparison (Accuracy, F1, Cohen's κ)](figures/fig3_model_comparison.png)

**Key observations**:

1. CSP+SVM achieves the highest mean accuracy (0.813) but with substantial fold-to-fold variance (SD=0.102), including one fold achieving 1.000 accuracy (see Discussion for analysis).
2. EEGNet achieves moderate accuracy (0.570) with high variance (SD=0.164), reflecting the challenge of training a neural network with only ~230 training trials per fold.
3. EEG Conformer achieves the lowest mean accuracy (0.531) but lowest variance (SD=0.036), suggesting more stable but less optimized learning.
4. All models substantially exceed chance (0.250), confirming the discriminability of the synthetic neural patterns.

**Figure 4** shows the average confusion matrix for EEGNet across 5 folds. Off-diagonal confusions are most prevalent between left-hand and right-hand classes (spatially antagonistic patterns), and between feet and tongue (both involving distal or non-limb imagery).

![Figure 4: EEGNet Confusion Matrix](figures/fig4_confusion_matrix.png)

### 5.2 P300 Speller Classification

**Table 2: P300 Speller Classification (5-Fold CV, Binary)**

| Method | Accuracy (mean±SD) | Notes |
|--------|-------------------|-------|
| Chance | 0.167 | Target rate 1:5 |
| LDA w/o TL | **0.970 ± 0.015** | 50 target labels only |
| LDA w/ EA-TL | 0.845 ± 0.017 | EA from source domain |

![Figure 5: P300 ERP and Transfer Learning](figures/fig5_p300_tl.png)

The left panel of Figure 5 shows the grand-average ERP for target vs. non-target stimuli at the pseudo-Pz channel, clearly revealing the N200–P300 complex with a peak at ~300 ms. Notably, Euclidean Alignment (EA-TL) *decreases* accuracy relative to no-TL (0.845 vs. 0.970). This counter-intuitive result is analyzed in the Discussion section.

### 5.3 Online Learning and Concept Drift

**Table 3: Online Learning Under Concept Drift**

| Classifier | Pre-Drift Accuracy | Post-Drift Accuracy | Relative Drop |
|-----------|-------------------|--------------------|-|
| Static SVM | 0.769 | 0.303 | −60.6% |
| Adaptive LDA | 0.756 | 0.769 | +1.7% |

![Figure 6: Concept Drift and Online Adaptation](figures/fig6_online_learning.png)

The static classifier suffers catastrophic accuracy degradation post-drift (0.769 → 0.303), nearly collapsing to chance level. The adaptive classifier maintains stability (0.756 → 0.769), with the sliding-window retraining successfully tracking the distributional shift. Five drift detection events were flagged, though some of these may be false positives triggered by the inherent variance of the 4-class problem.

**Figure 7** shows training dynamics for EEGNet and EEG Conformer on the final cross-validation fold. Both models converge within 40 epochs; EEGNet shows faster initial convergence while EEG Conformer stabilizes at a slightly lower accuracy.

![Figure 7: Training Curves](figures/fig7_training_curves.png)

---

## 6. Discussion

### 6.1 Interpretation of Results

**CSP+SVM dominance**: The superior performance of CSP+SVM (0.813) over deep learning methods reflects a fundamental constraint: the number of training trials (230 per fold) is insufficient to fully leverage the representational capacity of neural networks. Deep learning methods for EEG typically require hundreds to thousands of trials per subject to match CSP-based pipelines [5]. The one perfect fold (CSP+SVM = 1.000) warrants investigation—with only 57 test trials, a fold accuracy of 1.0 has probability (0.813)^57 ≈ 0.0001 under independent sampling, suggesting that this fold's train-test split may have produced particularly favorable class separability by chance. This highlights the importance of reporting standard deviations alongside mean accuracy.

**Deep learning underperformance**: EEGNet's high fold-to-fold variance (SD=0.164) and EEG Conformer's low accuracy (0.531) are consistent with findings in the literature for small datasets. Keutayeva et al. [11] achieve 70.12% with CCT on BCI IV 2a under LOSO evaluation, which involves more training data per test than 5-fold CV. The EEG Conformer's lower variance (SD=0.036) compared to EEGNet may reflect its larger parameter count providing more stable but underfitted solutions.

**P300 TL result**: The finding that EA-TL reduces P300 accuracy (0.845 vs. 0.970) is counter-intuitive but has a plausible explanation: Euclidean Alignment uses the source covariance as a reference, and when the source and target domains differ substantially (different random seeds, noise levels), the alignment can distort target-domain geometry rather than normalize it. This is consistent with the observation of Saha and Baumert [8] that transfer learning can be detrimental when domain gap is large or the alignment method does not match the type of distributional shift. In practice, iterative or alignment-aware approaches (e.g., CORAL, DANN) may be more robust.

**Online adaptation**: The dramatic recovery of post-drift accuracy from 0.303 to 0.769 with adaptive LDA demonstrates the importance of continual adaptation for long-duration BCI operation. In real clinical settings, concept drift arises from electrode impedance changes, fatigue, and neural plasticity—all continuous and unpredictable processes [4].

### 6.2 Critical Analysis: Limitations and Validity Threats

This study has several important limitations that constrain its generalizability:

**1. Synthetic data dependency**: All results are obtained on synthetically generated EEG, not real physiological recordings. The ground-truth class-specific topographic patterns are simplified (superposition of two sinusoids with additive noise), lacking the full complexity of real motor imagery EEG which involves: distributed cortical networks, subcortical modulations, trial-to-trial phase variability, non-Gaussian amplitude distributions, and complex artifact morphologies. Real EEG typically yields substantially lower CSP+SVM accuracy (0.55–0.75 on BCI IV 2a [5,6]) compared to our 0.813.

**2. Optimistic bias risk**: With only 288 trials and 5-fold CV (230 training / 57 test), confidence intervals are wide and results are sensitive to the specific train-test split. The one fold achieving 1.000 accuracy on CSP+SVM is anomalous and would not survive in a larger trial budget. We report standard deviations to quantify this uncertainty.

**3. No real-time latency evaluation**: While we describe an online processing architecture, the experiment does not measure actual processing latency. Real-time constraints (<500 ms for BCI, <200 ms for smooth motor control) require profiling on target hardware.

**4. Single-subject simulation**: All synthetic data is generated from a single statistical model with fixed SNR. Real BCIs face extreme inter-subject variability (BCI illiteracy rates of 15–30% [1]) that cannot be captured by a single synthetic subject.

**5. Artifact removal validation**: Our ASR/ICA simulation is highly simplified. Real ICA requires convergence over multiple trials, component selection criteria (e.g., dipole fitting, scalp topography), and is sensitive to the number of channels and data length. The artifact removal efficacy in our pipeline is not independently validated.

**6. P300 evaluation**: The high no-TL accuracy (0.970) likely reflects the idealized P300 generation with clean Gaussian ERP template—real P300 classification with LDA on 50 training trials typically yields 60–75% accuracy for untrained subjects, rising to 85–90% after practice [1].

### 6.3 Implications for Locked-In Syndrome

For patients with complete LIS, voluntary muscle control is absent, rendering ERP-based paradigms (P300, SSVEP) dependent on residual visual attention. For ALS patients in the late stage, even controlled eye movements may be compromised. Several design implications follow from this work:

- **Minimal-fatigue paradigms**: P300 spellers require extended recording sessions; the low-data transfer learning approach (50 training trials) is critical for clinical feasibility.
- **Robust drift handling**: LIS patients cannot easily re-calibrate; the online adaptation component is essential for multi-hour operation.
- **Multi-modal fusion**: Combining EEG with fNIRS or EMG (when residual) can improve classification robustness.
- **Personalized architectures**: The high inter-subject variability suggests that patient-specific fine-tuning is preferable to zero-shot cross-subject models.

### 6.4 Comparison with Prior Work

Our CSP+SVM accuracy (0.813 on synthetic) is higher than typical real-data benchmarks (0.55–0.75 on BCI IV 2a), confirming that synthetic data is easier. EEGNet's 0.570 is below published results on real data (0.67–0.72 in Lawhern et al. [5]), which may seem paradoxical but reflects our small trial count. EEG Conformer's published accuracy of 0.775 on BCI IV 2a [6] far exceeds our 0.531, primarily because Song et al. train with cross-subject data (N=9 subjects combined) whereas we use only a single-subject trial budget.

### 6.5 Future Directions

1. **Validation on real EEG**: Apply pipeline to BCI Competition IV 2a (MNE-Python compatible format), PhysioNet MI dataset, or BCI2000 P300 dataset.
2. **Riemannian geometry**: Covariance-based Riemannian methods (MDM, FgMDM) have shown state-of-the-art performance for small datasets and are natural extensions of CSP [Pan et al., 2025].
3. **Federated learning**: For patient privacy in clinical settings, federated or privacy-preserving transfer learning across multiple BCI users.
4. **Hardware-in-the-loop**: Integration with MNE-Realtime or BrainFlow for genuine real-time streaming, profiling end-to-end latency.
5. **Attention visualization**: EEG Conformer's class activation mapping could identify clinically relevant temporal-spatial patterns for LIS patients.

---

## 7. Conclusion

We presented an integrated real-time EEG processing pipeline for non-invasive BCI, combining artifact removal (ASR/ICA), CSP+SVM, EEGNet, EEG Conformer, Euclidean Alignment transfer learning for P300, and incremental online adaptation. On challenging synthetic 4-class motor imagery data (SNR ≈ 3 dB, 288 trials), CSP+SVM achieves the highest accuracy (0.813±0.102, κ=0.752), while deep learning methods underperform due to limited training data. P300 classification achieves 0.970±0.015 accuracy; Euclidean Alignment reduces accuracy in this scenario, warranting caution in applying alignment under large domain gaps. The adaptive online learner recovers post-drift accuracy to 0.769 versus 0.303 for the static classifier, demonstrating the critical importance of concept drift adaptation for clinical BCI deployment.

These results, however, must be interpreted with caution: synthetic data evaluation overestimates real-world performance, and the deep learning models are severely data-limited. The path from laboratory demonstration to clinical assistive device for locked-in syndrome patients requires validation on real physiological data, real-time latency profiling, and longitudinal evaluation across sessions. The modular pipeline presented here provides a reproducible foundation for such future work.

---

## References

[1] Rashid, M., Sulaiman, N., Abdul Majeed, A. P. P., et al. (2020). Current Status, Challenges, and Possible Solutions of EEG-Based Brain-Computer Interface: A Comprehensive Review. *Frontiers in Neurorobotics*, 14, 25. https://doi.org/10.3389/fnbot.2020.00025

[2] Saha, S., Mamun, K. A., Ahmed, K., et al. (2021). Progress in Brain Computer Interface: Challenges and Opportunities. *Frontiers in Systems Neuroscience*, 15, 578875. https://doi.org/10.3389/fnsys.2021.578875

[3] Gu, X., Cao, Z., Jolfaei, A., et al. (2021). EEG-based Brain-Computer Interfaces (BCIs): A Survey of Recent Studies on Signal Sensing Technologies and Computational Intelligence Approaches and their Applications. *IEEE Transactions on Computational Biology and Bioinformatics*, 18(5), 2232–2249. https://doi.org/10.1109/tcbb.2021.3052811

[4] Värbu, K., Muhammad, N., & Muhammad, Y. (2022). Past, Present, and Future of EEG-Based BCI Applications. *Sensors*, 22(9), 3331. https://doi.org/10.3390/s22093331

[5] Lawhern, V. J., Solon, A. J., Waytowich, N. R., et al. (2018). EEGNet: a compact convolutional neural network for EEG-based brain–computer interfaces. *Journal of Neural Engineering*, 15(5), 056013. https://doi.org/10.1088/1741-2552/aace8c

[6] Song, Y., Zheng, Q., Liu, B., & Gao, X. (2022). EEG Conformer: Convolutional Transformer for EEG Decoding and Visualization. *IEEE Transactions on Neural Systems and Rehabilitation Engineering*, 31, 710–719. https://doi.org/10.1109/tnsre.2022.3230250

[7] Saha, S., & Baumert, M. (2020). Intra- and Inter-subject Variability in EEG-Based Sensorimotor Brain Computer Interface: A Review. *Frontiers in Computational Neuroscience*, 13, 87. https://doi.org/10.3389/fncom.2019.00087

[8] Pfeffer, M. A., Ling, S. H., & Wong, J. K. W. (2024). Exploring the frontier: Transformer-based models in EEG signal analysis for brain-computer interfaces. *Computers in Biology and Medicine*, 176, 108705. https://doi.org/10.1016/j.compbiomed.2024.108705

[9] Keutayeva, A., Fakhrutdinov, N., & Abibullaev, B. (2024). Compact convolutional transformer for subject-independent motor imagery EEG-based BCIs. *Scientific Reports*, 14, 22982. https://doi.org/10.1038/s41598-024-73755-4

[10] Lee, P.-L., Chen, S.-H., Chang, T.-C., et al. (2023). Continual Learning of a Transformer-Based Deep Learning Classifier Using an Initial Model from Action Observation EEG Data to Online Motor Imagery Classification. *Bioengineering*, 10(2), 186. https://doi.org/10.3390/bioengineering10020186

[11] Nallani, S. V. C., & Ramachandran, G. (2024). RLEEGNet: Integrating Brain-Computer Interfaces with Adaptive AI for Intuitive Responsiveness and High-Accuracy Motor Imagery Classification. *TechRxiv*. https://doi.org/10.36227/techrxiv.172651885.52358678/v1

[12] Pan, L., Wang, K., Huang, Y., et al. (2025). Enhancing motor imagery EEG classification with a Riemannian geometry-based spatial filtering (RSF) method. *Neural Networks*, 107511. https://doi.org/10.1016/j.neunet.2025.107511
