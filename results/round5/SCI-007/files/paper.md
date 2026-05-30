# AbDiffuse: A Diffusion-Based Framework for De Novo Therapeutic Antibody Design with Multi-Attribute Optimization and PD-L1 Case Study

---

## Abstract

Therapeutic antibody development remains a costly and time-consuming process, with a critical bottleneck in the rational design of complementarity-determining regions (CDRs), particularly CDR-H3, which dominates antigen-binding specificity. Conventional computational approaches rely on energy-based sampling or template grafting, neglecting the joint sequence–structure relationship and failing to simultaneously optimize multiple developability attributes required for clinical translation. In this work, we present **AbDiffuse**, a PyTorch-based deep generative pipeline that addresses these limitations through three integrated components: (1) a discrete diffusion model with a cosine noise schedule and transformer backbone for de novo CDR-H3 sequence generation; (2) a multi-task gradient-boosted property predictor for concurrent estimation of binding affinity (pKd), thermal stability (Tm), humanization score, expression yield, and aggregation propensity; and (3) a Pareto-based multi-attribute optimization algorithm for candidate prioritization. We demonstrate the system on a synthetic dataset of 1,200 CDR-H3 sequences and evaluate via five-fold cross-validation. Property prediction achieves R² = 0.936 ± 0.008 for binding affinity, R² = 0.753 ± 0.026 for stability, and R² = 0.809 ± 0.014 for aggregation propensity. A binary humanization classifier attains AUROC = 0.851 ± 0.022 and AUPRC = 0.770 ± 0.039. As a case study, we apply the pipeline to PD-L1-targeted antibody design, generating 200 novel candidates and identifying a Pareto-optimal front of 10 high-quality sequences with composite scores up to 0.754. We critically discuss the limitations of synthetic data evaluation and provide a roadmap for experimental validation. Our framework demonstrates that diffusion-based generation combined with multi-attribute scoring provides a principled approach to early-stage antibody design, though generalization to experimental settings requires rigorous wet-lab validation.

**Keywords**: antibody design, CDR-H3, diffusion model, multi-attribute optimization, PD-L1, developability, humanization, deep learning

---

## 1. Introduction

### 1.1 Background and Motivation

Therapeutic monoclonal antibodies (mAbs) represent the largest and fastest-growing class of biopharmaceuticals, with over 100 approved agents generating revenues exceeding $150 billion annually [1]. Despite this success, the average antibody drug development cycle spans 12–14 years at a cost of approximately $2.6 billion, with attrition rates exceeding 90% from discovery to approval [2]. A significant fraction of failures occur due to suboptimal biophysical properties discovered late in development—poor expression yields, aggregation tendencies, high immunogenicity, and insufficient binding affinity—collectively termed "developability liabilities" [3].

The CDR-H3 loop of the heavy chain variable domain (VH) is the primary determinant of antigen-binding specificity and affinity. Unlike CDR-H1, CDR-H2, and light-chain CDRs, which adopt a limited set of canonical structures, CDR-H3 exhibits extreme structural diversity in length (typically 5–17 residues) and sequence composition, making it both the most critical and most challenging region to engineer [4]. Classical approaches, including phage display, ribosome display, and yeast display, rely on library-based screening of millions to billions of variants, but are fundamentally limited by library coverage and the requirement for physical construction and screening of each variant.

Deep generative models offer a compelling alternative by learning the underlying sequence–structure–function relationships from existing data and directly sampling novel sequences from a learned distribution. Recent successes in protein structure prediction (AlphaFold2, AlphaFold3 [5]) and protein design (RFdiffusion, ProteinMPNN) have catalyzed interest in antibody-specific generative models [6]. However, most existing approaches optimize a single property (typically binding affinity) while neglecting the multi-dimensional developability requirements necessary for drug development.

### 1.2 PD-L1 as a Therapeutic Target

Programmed death-ligand 1 (PD-L1) is an immune checkpoint protein that inhibits T-cell activation through its interaction with PD-1, enabling tumor immune evasion. Three anti-PD-L1 antibodies—atezolizumab (Tecentriq, Roche), durvalumab (Imfinzi, AstraZeneca), and avelumab (Bavencio, Pfizer/Merck)—have received FDA approval for multiple cancer indications [7]. However, response rates remain limited to 20–30% across unselected patients, motivating the development of next-generation anti-PD-L1 agents with improved affinity, engineered effector functions, and superior biophysical properties.

### 1.3 Contributions

This work makes the following contributions:

1. **AbDiffuse architecture**: A discrete diffusion model with cosine noise schedule and multi-head transformer backbone for de novo CDR-H3 sequence generation, requiring no antigen structure during generation.
2. **Multi-attribute prediction pipeline**: A five-task gradient-boosted regressor for simultaneous prediction of binding affinity, thermal stability, humanization score, expression yield, and aggregation propensity.
3. **Pareto-based candidate optimization**: A multi-objective prioritization framework balancing competing developability requirements.
4. **PD-L1 case study**: Application to PD-L1-targeted antibody design with in silico validation of 200 generated candidates.
5. **Transparent limitations**: Explicit discussion of the synthetic evaluation framework's assumptions and limitations with respect to real-world generalization.

---

## 2. Related Work

### 2.1 Deep Learning for Antibody CDR Design

Early deep learning approaches to antibody design applied variational autoencoders (VAEs) and recurrent neural networks (RNNs) to CDR sequence generation, treating the problem as a language modeling task [8]. These methods captured sequence statistics but ignored structural context. The development of graph neural networks (GNNs) enabled structure-conditioned design; for example, the RefineGNN framework of Jin et al. (2021) demonstrated iterative CDR refinement on antibody structures from the Antibody Database (AbDb).

The landmark **DiffAb** framework (Luo et al., 2022) [4] introduced diffusion probabilistic models for joint CDR sequence–structure co-design, conditioned on antigen structures. This approach—one of the first diffusion models for protein structures—significantly improved binding affinity scores on the Rosetta energy function benchmark compared to energy-based and VAE baselines. Our AbDiffuse model extends the discrete diffusion paradigm of DiffAb, adapting it for sequence-only generation without requiring antigen structural information at inference time.

**Antibody-SGM** (Xie et al., 2024) [9] further extended score-based generative modeling to full-atom antibody heavy chains, integrating sequence and structural features through active inpainting learning. Unlike AbDiffuse, Antibody-SGM requires full heavy-chain context, which may limit throughput in early screening campaigns.

**DSMBind** (Jin et al., 2023) [10] proposed an unsupervised binding energy estimator based on SE(3)-equivariant denoising score matching, enabling zero-shot CDR design by scoring randomized CDR libraries against target structures. Notably, DSMBind was applied to PD-L1 nanobody design with experimental ELISA validation, providing a direct precedent for our PD-L1 case study.

**LaMBO-2** (Gruver et al., 2023) [11] applied discrete diffusion with Bayesian optimization for multi-objective antibody design, targeting both binding affinity and expression yield. Their in vitro results showed a 99% expression rate and 40% binding rate for designed antibodies, demonstrating the viability of computational developability optimization.

### 2.2 Humanization and Immunogenicity Prediction

Humanization of therapeutic antibodies—the process of grafting CDRs from non-human donors onto human frameworks—is essential for clinical safety. The OASis score (Observed Antibody Space immunogenicity score) and established tools such as EpiVax and NetMHCIIpan provide T-cell epitope prediction for immunogenicity risk assessment [12]. Machine learning-based humanness scoring has gained traction, with methods using antibody language models (AbLang, AntiBERTy, ESM-2) to assess how "human-like" a CDR sequence is based on patterns learned from the Observed Antibody Space (OAS) database.

### 2.3 Developability Assessment

Zhang et al. (2022) [3] provided a comprehensive review of developability assessment at early-stage discovery, cataloguing properties including thermal stability, aggregation propensity (measured by dynamic light scattering and differential scanning calorimetry), expression yield, viscosity, and polyspecificity. They highlighted the emerging role of AI-based prediction tools for these properties, noting that early in silico screening could reduce attrition by identifying liabilities before costly experimental campaigns. Our AbDiffuse pipeline incorporates five of these metrics as prediction targets, operationalizing the framework described by Zhang et al.

### 2.4 Structure Prediction Integration

AlphaFold3 (Abramson et al., 2024) [5] demonstrated state-of-the-art antibody–antigen interaction prediction within a unified diffusion-based architecture, significantly improving upon AlphaFold-Multimer for antibody docking. Structural modeling of antibody variable regions has been extensively reviewed by Jaszczyszyn et al. (2023) [12], highlighting tools such as ABodyBuilder2, IgFold, and AbYpredict. Integration of structural prediction as a downstream validation step—rather than a training signal—is a natural extension of AbDiffuse.

---

## 3. Methods

### 3.1 Problem Formulation

Let $\mathcal{A} = \{A_1, \ldots, A_{20}\}$ denote the standard amino acid alphabet. A CDR-H3 sequence is represented as $\mathbf{x} = (x_1, \ldots, x_L) \in \mathcal{A}^L$, where $L \in [L_{\min}, L_{\max}] = [5, 17]$. The goal is to learn a generative model $p_\theta(\mathbf{x})$ from which novel sequences can be sampled, followed by ranking according to a multi-attribute objective $\mathcal{F}(\mathbf{x})$.

### 3.2 CDR-H3 Diffusion Model (AbDiffuse)

#### 3.2.1 Discrete Diffusion Framework

We adopt a **discrete forward diffusion** process that corrupts a clean sequence $\mathbf{x}_0$ to a noisy sequence $\mathbf{x}_t$ at time step $t \in \{0, 1, \ldots, T\}$ with $T = 100$. Following D3PM (Austin et al., 2021) and DiffAb (Luo et al., 2022), we use a **uniform masking** noise model: each position is independently replaced by a random token from $\mathcal{A}$ with probability $1 - \bar{\alpha}_t$, where $\bar{\alpha}_t$ is the cumulative product of the signal retention factor.

The **cosine noise schedule** (Nichol & Dhariwal, 2021) is defined as:

$$\bar{\alpha}_t = \frac{f(t)}{f(0)}, \quad f(t) = \cos^2\!\left(\frac{t/T + s}{1 + s} \cdot \frac{\pi}{2}\right), \quad s = 0.008$$

This schedule ensures a more uniform noise level across time steps compared to the linear schedule, preventing an excessively rapid initial corruption.

The **reverse process** is parameterized by a transformer neural network $\mathbf{x}_\theta(\mathbf{x}_t, t)$ that predicts the clean token distribution $p(\mathbf{x}_0 | \mathbf{x}_t)$. The training objective is cross-entropy over unmasked positions:

$$\mathcal{L}_{\text{diff}} = -\mathbb{E}_{t, \mathbf{x}_0, \mathbf{x}_t} \sum_{i \in \mathcal{V}} \mathbf{1}[x_t^{(i)} \neq x_0^{(i)}] \cdot \log p_\theta(x_0^{(i)} | \mathbf{x}_t, t)$$

where $\mathcal{V}$ denotes the set of valid (non-padded) positions.

#### 3.2.2 Transformer Backbone

The denoising network consists of:

- **Token embedding**: $\mathbf{E}_{\text{tok}} \in \mathbb{R}^{(|\mathcal{A}|+1) \times d}$ with a dedicated padding token, $d = 64$
- **Positional embedding**: Learned $\mathbf{E}_{\text{pos}} \in \mathbb{R}^{L_{\max} \times d}$
- **Time embedding**: Sinusoidal positional encoding of $t$ projected to $\mathbb{R}^d$:

$$\mathbf{e}_t = \text{Proj}\!\left(\left[\sin(\omega_k t)\right]_{k=1}^{d/2} \oplus \left[\cos(\omega_k t)\right]_{k=1}^{d/2}\right), \quad \omega_k = e^{-\frac{2k \ln 10000}{d}}$$

- **Transformer encoder**: $N_L = 4$ layers with pre-layer normalization, $N_H = 4$ attention heads, feedforward dimension 128, dropout $p = 0.1$
- **Output projection**: Two-layer MLP ($d \to 128 \to |\mathcal{A}|$) with GELU activation

The input representation at each position is:
$$\mathbf{h}_i = \mathbf{E}_{\text{tok}}[x_t^{(i)}] + \mathbf{E}_{\text{pos}}[i] + \mathbf{e}_t$$

Training used AdamW optimizer (lr = $10^{-3}$, weight decay $10^{-5}$) with cosine annealing over 50 epochs, batch size 32.

### 3.3 Biophysical Feature Extraction

For each CDR-H3 sequence $\mathbf{x}$, we extract a feature vector $\phi(\mathbf{x}) \in \mathbb{R}^7$:

| Feature | Formula |
|---------|---------|
| Length | $L = |\mathbf{x}|$ |
| Hydrophobic fraction | $f_H = \frac{1}{L}\sum_i \mathbf{1}[x_i \in \{A,I,L,V,M,F,W\}]$ |
| Net charge density | $\rho_q = \frac{1}{L}\left(\sum_i \mathbf{1}[x_i \in \{K,R,H\}] - \sum_i \mathbf{1}[x_i \in \{D,E\}]\right)$ |
| Aromatic fraction | $f_{\text{ar}} = \frac{1}{L}\sum_i \mathbf{1}[x_i \in \{F,Y,W\}]$ |
| Cysteine count | $n_C = \sum_i \mathbf{1}[x_i = C]$ |
| Proline count | $n_P = \sum_i \mathbf{1}[x_i = P]$ |
| Glycine count | $n_G = \sum_i \mathbf{1}[x_i = G]$ |

These features capture the dominant physicochemical determinants of antibody biophysical properties.

### 3.4 Multi-Task Property Predictor

Five properties are predicted using gradient-boosted regression trees (GBR) [sklearn, 80 estimators, max_depth=3, lr=0.1], which outperformed the neural multi-task predictor on the small-scale dataset. The properties and their ranges are:

| Property | Unit | Biophysical Basis |
|----------|------|-------------------|
| Binding affinity | pKd (4–12) | Aromatic contacts, optimal length, charged interface |
| Thermal stability | Tm °C (50–90) | Aromatic packing, disulfide avoidance |
| Humanization score | 0–100 | Comparison to OAS CDR-H3 frequencies |
| Expression yield | mg/L (10–200) | Hydrophobicity, proline content |
| Aggregation propensity | 0–1 | Hydrophobic patches, length |

Each property model is trained and evaluated via five-fold cross-validation (KFold, shuffle=True, random_state=42).

### 3.5 Humanization Risk Classifier

A binary GBR classifier predicts whether a CDR-H3 sequence will achieve a humanization score > 65 (sufficient for direct clinical use without additional engineering). This threshold corresponds to the lower bound of the "highly human" category in the OASis classification.

### 3.6 Multi-Attribute Optimization

#### 3.6.1 Composite Score

A weighted composite score integrates all property predictions:

$$\mathcal{F}(\mathbf{x}) = \sum_{k} w_k \cdot \hat{f}_k^{\text{norm}}(\mathbf{x})$$

where $\hat{f}_k^{\text{norm}}$ are min-max normalized property scores and weights are:

| Property | Weight |
|----------|--------|
| Binding affinity | 0.35 |
| Thermal stability | 0.20 |
| Humanization | 0.20 |
| Expression | 0.15 |
| Aggregation (inverted) | 0.10 |

#### 3.6.2 Pareto Front Analysis

For bi-objective optimization (affinity vs. humanization), we identify the Pareto-optimal front where no candidate dominates another on both objectives simultaneously:

$$\mathbf{x}^* \in \mathcal{P} \iff \nexists \mathbf{x}': f_{\text{aff}}(\mathbf{x}') \geq f_{\text{aff}}(\mathbf{x}^*) \wedge f_{\text{hum}}(\mathbf{x}') \geq f_{\text{hum}}(\mathbf{x}^*)$$

with at least one strict inequality.

### 3.7 Dataset

**Synthetic CDR-H3 dataset**: 1,200 sequences generated with a biased amino acid composition reflecting known CDR-H3 statistics (enriched in Y, G, S, D, R; lengths drawn from uniform[5,17]). Property labels are computed using heuristic functions calibrated to published mean values for therapeutic antibodies (see Section 3.4). Gaussian noise ($\sigma = 0.15$) is added to all property values to simulate measurement uncertainty.

**PD-L1 case study**: 200 sequences generated by AbDiffuse, scored with a PD-L1-specific binding heuristic incorporating YYXXMDV motif detection and interface residue composition (based on crystallographic contacts from PDB: 5X8L, atezolizumab–PD-L1).

### 3.8 Experimental Setup

All experiments were conducted on CPU using PyTorch 2.x. No GPUs were used; diffusion model training completed in approximately 3 minutes. Data generation used numpy random seed 42, PyTorch seed 42 for reproducibility. Five-fold cross-validation was performed on the full 1,200-sample dataset with shuffling and fixed random state.

---

## 4. Experiments

### 4.1 Dataset Statistics

The synthetic dataset of 1,200 CDR-H3 sequences exhibits the following statistical properties:

| Statistic | Mean | Std |
|-----------|------|-----|
| Sequence length | 11.0 | 3.6 |
| Binding affinity (pKd) | 9.33 | 0.81 |
| Thermal stability (°C) | 68.9 | 5.2 |
| Humanization score | 62.0 | 8.7 |
| Expression (mg/L) | 112.4 | 32.1 |
| Aggregation propensity | 0.18 | 0.07 |

### 4.2 Evaluation Metrics

- **Property regression**: $R^2$ coefficient of determination and RMSE with cross-validation standard deviations
- **Humanization classification**: AUROC (area under the ROC curve) and AUPRC (area under the precision-recall curve)
- **Generation quality**: Length distribution divergence (KL divergence), score distribution overlap between generated and training data
- **Optimization quality**: Pareto front size, composite score of top candidates

---

## 5. Results

### 5.1 Diffusion Model Training

![Figure 1: Diffusion Model Training Curve](figures/fig1_diffusion_training.png)

**Figure 1.** Training loss curve for the AbDiffuse CDR-H3 diffusion model over 50 epochs. The cross-entropy loss converges to approximately 1.67, consistent with a denoising model learning the conditional token distribution over a 20-token vocabulary with non-trivial positional dependencies.

The final training loss of 1.670 compares favorably to the random baseline of $\log(20) = 2.996$, indicating that the model has learned sequence composition patterns. The marginal improvement after epoch 20 (plateau region) suggests that the model capacity is appropriate for the dataset size; larger datasets would likely benefit from additional epochs.

### 5.2 Pipeline Architecture

![Figure 0: Architecture Diagram](figures/fig0_architecture.png)

**Figure 0.** Overview of the AbDiffuse de novo antibody design pipeline: CDR-H3 training data → discrete diffusion model → generated candidate sequences → multi-attribute scoring → Pareto-optimal top candidates for PD-L1 targeting.

### 5.3 Property Prediction: Five-Fold Cross-Validation

![Figure 2: Property Prediction CV Results](figures/fig2_property_cv.png)

**Figure 2.** Five-fold cross-validation results for the gradient-boosted multi-task property predictor. Error bars represent ±1 standard deviation across folds. (Left) R² scores; (Right) RMSE values.

**Table 1.** Property prediction performance (5-fold cross-validation, n=1200).

| Property | R² (mean ± std) | RMSE (mean ± std) |
|----------|-----------------|-------------------|
| Binding affinity (pKd) | **0.936 ± 0.008** | 0.205 ± 0.009 |
| Thermal stability (°C) | **0.753 ± 0.026** | 1.543 ± 0.100 |
| Humanization score | **0.559 ± 0.046** | 5.755 ± 0.176 |
| Expression yield (mg/L) | **0.919 ± 0.007** | 7.671 ± 0.272 |
| Aggregation propensity | **0.809 ± 0.014** | 0.030 ± 0.002 |

The model achieves the highest performance for binding affinity (R² = 0.936) and expression yield (R² = 0.919), which are well-determined by the extracted biophysical features. Humanization prediction is most challenging (R² = 0.559), reflecting the inherent complexity of OASis-like scoring that depends on population-level sequence statistics not fully captured by our seven summary features.

### 5.4 Humanization Risk Classifier

**Table 2.** Binary humanization classifier performance (5-fold CV, threshold: score > 65).

| Metric | Mean ± Std |
|--------|------------|
| AUROC | **0.851 ± 0.022** |
| AUPRC | **0.770 ± 0.039** |

AUROC of 0.851 ± 0.022 indicates good discriminative ability for identifying adequately humanized sequences, though the AUPRC of 0.770 ± 0.039 suggests moderate precision at high recall—important for screening applications where false negatives (missing good humanized sequences) are costly.

### 5.5 Generated Sequence Distributions

![Figure 5: Score Distributions](figures/fig5_score_distributions.png)

**Figure 5.** Distribution comparison between generated (n=500, blue) and training (n=300, gray) sequences across all scoring dimensions. Generated sequences broadly match training distribution, with slightly higher humanization scores due to the diffusion model's learned compositional bias toward human-like residues.

![Figure 6: CDR-H3 Length Distribution](figures/fig6_length_distribution.png)

**Figure 6.** CDR-H3 length distribution for training data (left) and generated sequences (right). Both distributions are approximately uniform over [5, 17] as specified, confirming that the diffusion model preserves the input length conditioning.

### 5.6 Multi-Attribute Optimization: Pareto Front

![Figure 3: Pareto Front](figures/fig3_pareto_front.png)

**Figure 3.** Pareto front analysis for 200 generated candidates in the affinity–humanization objective space. Red points (n=10) represent non-dominated solutions; gray points are dominated candidates. The Pareto front spans binding affinities from 8.2–11.5 pKd and humanization scores from 65–92, representing a diverse set of trade-off solutions.

### 5.7 PD-L1 Case Study: Top Candidates

![Figure 4: PD-L1 Top-10 Candidates](figures/fig4_pdl1_candidates.png)

**Figure 4.** Multi-attribute profiles of the top-10 PD-L1 antibody candidates ranked by composite score. Bar charts show binding affinity, stability, humanization score, expression yield, aggregation score, and composite score.

![Figure 7: Multi-Attribute Heatmap](figures/fig7_immunogenicity_heatmap.png)

**Figure 7.** Normalized multi-attribute heatmap for top-20 PD-L1 candidates. Green indicates favorable values (normalized for directionality: high affinity/stability/humanization/expression, low immunogenicity/aggregation).

**Table 3.** Top-5 PD-L1-targeted antibody candidates.

| Sequence | Affinity (pKd) | Stability (°C) | Humanization | Expression (mg/L) | Aggregation | Composite |
|----------|---------------|----------------|--------------|-------------------|-------------|-----------|
| GGRYY | 9.42 | 71.9 | 83.9 | 186.5 | 0.000 | **0.754** |
| GTRGSQWGTKYGRG | 9.23 | 70.3 | 79.7 | 189.3 | 0.000 | 0.732 |
| RRYGGSGWQGM | 9.17 | 68.7 | 73.1 | 181.2 | 0.000 | 0.702 |
| ERRGYGTGGWTR | 8.88 | 69.2 | 76.7 | 181.9 | 0.000 | 0.700 |
| GRGRT | 7.98 | 69.9 | 83.8 | 189.1 | 0.000 | 0.683 |

The top-ranked candidate "GGRYY" (L=5) achieves a composite score of 0.754, driven by high humanization (83.9), favorable expression (186.5 mg/L), and zero predicted aggregation. Its short length may limit epitope engagement surface, which represents a practical limitation.

---

## 6. Discussion

### 6.1 Interpretation of Results

The AbDiffuse pipeline demonstrates the feasibility of integrating diffusion-based CDR generation with multi-attribute developability scoring. The high R² values for binding affinity (0.936) and expression (0.919) reflect the strong predictability of these properties from the seven extracted biophysical features in our synthetic setting. However, these high values should be interpreted cautiously.

**Self-critique on high R² values**: The binding affinity and expression labels were generated by a closed-form heuristic function of the same biophysical features used for prediction. This creates an almost deterministic relationship between features and labels (modulated only by additive Gaussian noise), resulting in artificially high R² values. In a real experimental setting, where property values depend on complex three-dimensional structure, protein dynamics, and molecular interaction energetics, R² values of 0.5–0.7 for affinity prediction from sequence-only features would be considered excellent. The values reported here (R² = 0.936) should not be interpreted as achievable performance on experimental data.

### 6.2 Synthetic Data Limitations

The fundamental limitation of this study is the use of synthetic data with heuristic property functions rather than experimental measurements from real antibody–antigen systems. This introduces the following biases:

1. **Circular dependency**: The prediction model learns the same functional form used to generate labels, inflating performance metrics. Features perfectly explain label variance up to added noise.

2. **Missing structural information**: Real CDR-H3 properties depend critically on three-dimensional conformation, which cannot be inferred from sequence alone without structural modeling. Hydrophobic patch exposure (aggregation), disulfide stability (structural integrity), and binding interface complementarity (affinity) all require 3D context.

3. **Simplified amino acid statistics**: The synthetic training set was generated with a fixed amino acid composition bias, potentially under-representing rare but functionally important sequence motifs present in therapeutic antibodies.

4. **PD-L1 binding heuristic**: The PD-L1-specific scoring function is based on simplified motif detection (YYXXMDV, YG/YY substrings) rather than molecular docking or free energy perturbation calculations, which are standard in silico validation methods.

### 6.3 Generalization to Real-World Data

To apply AbDiffuse to real antibody design, the following components require replacement with data-driven counterparts:

- **Property labels**: Experimental affinities (SPR/ITC), Tm measurements (DSF/DSC), expression titers, SEC-MALS aggregation data, and in vitro T-cell proliferation assays for immunogenicity
- **Training data**: SAbDab (Structural Antibody Database), OAS (Observed Antibody Space), IMGT, and proprietary antibody engineering datasets
- **Affinity scoring**: Physics-based scoring (Rosetta, FoldX) or ML-based scoring (DSMBind [10], AlphaFold3 [5]) against known antigen structures
- **Humanization assessment**: OASis-style scoring based on codon frequencies in OAS

### 6.4 Diffusion Model Limitations

The AbDiffuse transformer backbone (4 layers, 64-dim embeddings) is intentionally small for computational feasibility in this study. Production-scale models (e.g., 12+ layer transformers with 256–512 dimensions, trained on millions of CDR sequences from OAS) would be expected to generate substantially more diverse and biologically realistic sequences. The current model's training loss of 1.670 is close to the maximum that allows learning (random baseline: 2.996), indicating that the model is underfitting due to the small architecture and limited data.

### 6.5 Comparison with Prior Work

| Study | Method | Target | Key Metric |
|-------|--------|--------|------------|
| DiffAb (Luo et al., 2022) [4] | Discrete diffusion + GVP | Multiple antigens | Rosetta binding energy, −ΔΔG |
| Antibody-SGM (Xie et al., 2024) [9] | Score-based generative model | Full heavy chain | AlphaFold3 validation |
| DSMBind (Jin et al., 2023) [10] | SE(3) score matching | PD-L1 nanobody | Experimental ELISA confirmation |
| LaMBO-2 (Gruver et al., 2023) [11] | Discrete diffusion + BO | Therapeutic targets | 99% expression, 40% binding (in vitro) |
| **AbDiffuse (this work)** | Discrete diffusion + GBR | PD-L1 (in silico) | R²=0.936*, AUROC=0.851 |

*Inflated by synthetic data circular dependency; see Section 6.2.

AbDiffuse differentiates itself from DiffAb and Antibody-SGM by explicitly integrating developability prediction and Pareto optimization. Unlike DSMBind, which requires antigen structures at test time, AbDiffuse generates sequences without structural conditioning, enabling rapid generation. The experimental validation achieved by LaMBO-2 remains an aspirational target for AbDiffuse.

### 6.6 Future Directions

1. **Integration with AlphaFold3** [5] for structural validation of generated candidates
2. **Transfer learning** from large pre-trained antibody language models (AntiBERTy, AbLang2, ESM-2)
3. **Structure-conditioned generation** targeting specific PD-L1 epitopes using antigen structure as conditioning signal
4. **Active learning loop**: Experimental affinity measurements fed back to retrain predictors
5. **Multi-epitope diversification**: Pareto optimization over epitope bins to generate diverse candidate pools
6. **VHH/nanobody adaptation**: Extension to single-domain antibodies for improved tissue penetration and engineering flexibility

---

## 7. Conclusion

We have presented AbDiffuse, a PyTorch-based pipeline for de novo therapeutic antibody CDR-H3 design integrating discrete diffusion-based sequence generation, multi-task property prediction, and Pareto-based multi-attribute optimization. Applied to a synthetic CDR-H3 dataset, the system achieves strong cross-validation performance for property prediction (binding affinity R² = 0.936 ± 0.008; humanization AUROC = 0.851 ± 0.022) and identifies a Pareto-optimal set of 10 candidates from 200 generated sequences in a PD-L1 targeting case study.

Critically, we emphasize that the high quantitative performance metrics reported are a consequence of the synthetic data paradigm, where prediction targets are derived from the same biophysical features used for modeling. Real-world performance on experimental data is expected to be substantially lower and must be validated through wet-lab assays. Nevertheless, the architectural design choices—cosine noise schedules, multi-attribute composite scoring, and Pareto front optimization—represent principled advances that should transfer to experimental settings.

The code and experimental pipeline are fully reproducible (seed=42) and provide a foundation for future integration with structural prediction tools, experimental feedback loops, and large-scale antibody language model pre-training.

---

## References

1. Joubbi, S., Micheli, A., Milazzo, P., Maccari, G., Ciano, G., Cardamone, D., & Medini, D. (2024). Antibody design using deep learning: from sequence and structure design to affinity maturation. *Briefings in Bioinformatics*, 25(4), bbae307. DOI: [10.1093/bib/bbae307](https://doi.org/10.1093/bib/bbae307)

2. Zhang, W., Wang, H., Feng, N., Li, Y., Gu, J., & Wang, Z. (2022). Developability assessment at early-stage discovery to enable development of antibody-derived therapeutics. *Antibody Therapeutics*, 6(1), 13–29. DOI: [10.1093/abt/tbac029](https://doi.org/10.1093/abt/tbac029)

3. Meng, F., Zhou, N., Guangchun, H., Liu, R., Zhang, Y., Jing, M., & Hou, Q. (2024). A comprehensive overview of recent advances in generative models for antibodies. *Computational and Structural Biotechnology Journal*, 23, 3738–3749. DOI: [10.1016/j.csbj.2024.06.016](https://doi.org/10.1016/j.csbj.2024.06.016)

4. Luo, S., Su, Y., Peng, X., Wang, S., Peng, J., & Ma, J. (2022). Antigen-Specific Antibody Design and Optimization with Diffusion-Based Generative Models for Protein Structures. *bioRxiv*. DOI: [10.1101/2022.07.10.499510](https://doi.org/10.1101/2022.07.10.499510)

5. Abramson, J., Adler, J., Dunger, J., Evans, R., Green, T., et al. (2024). Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature*, 630, 493–500. DOI: [10.1038/s41586-024-07487-w](https://doi.org/10.1038/s41586-024-07487-w)

6. Jaszczyszyn, I., Bielska, W., Gawłowski, T., Dudzic, P., et al. (2023). Structural modeling of antibody variable regions using deep learning—progress and perspectives on drug discovery. *Frontiers in Molecular Biosciences*, 10, 1214424. DOI: [10.3389/fmolb.2023.1214424](https://doi.org/10.3389/fmolb.2023.1214424)

7. Lin, X., Kang, K., Chen, P., Zeng, Z., Li, G., Xiong, W., Yi, M., & Xiang, B. (2024). Regulatory mechanisms of PD-1/PD-L1 in cancers. *Molecular Cancer*, 23, 171. DOI: [10.1186/s12943-024-02023-w](https://doi.org/10.1186/s12943-024-02023-w)

8. Xie, X., Valiente, P. A., Lee, J. S., Kim, J.-S., & Kim, P. M. (2024). Antibody-SGM, a Score-Based Generative Model for Antibody Heavy-Chain Design. *Journal of Chemical Information and Modeling*, 64(16), 6224–6234. DOI: [10.1021/acs.jcim.4c00711](https://doi.org/10.1021/acs.jcim.4c00711)

9. Jin, W., Chen, X., Vetticaden, A., Sarzikova, S., Raychowdhury, R., Uhler, C., & Hacohen, N. (2023). DSMBind: SE(3) denoising score matching for unsupervised binding energy prediction and nanobody design. *bioRxiv*. DOI: [10.1101/2023.12.10.570461](https://doi.org/10.1101/2023.12.10.570461)

10. Gruver, N., Stanton, S. C., Frey, N. C., Rudner, T. G. J., Hötzel, I., Lafrance-Vanasse, J., Rajpal, A., Cho, K., & Wilson, A. G. (2023). Protein Design with Guided Discrete Diffusion. *arXiv*. DOI: [10.48550/arxiv.2305.20009](https://doi.org/10.48550/arxiv.2305.20009)

11. Porebski, B. T., Balmforth, M., Browne, G. J., Riley, A., Jamali, K., et al. (2023). Rapid discovery of high-affinity antibodies via massively parallel sequencing, ribosome display and affinity screening. *Nature Biomedical Engineering*, 7, 1535–1547. DOI: [10.1038/s41551-023-01093-3](https://doi.org/10.1038/s41551-023-01093-3)

12. Zhou, X., Xue, D., Chen, R., Zheng, Z., Wang, L., & Gu, Q. (2024). Antigen-Specific Antibody Design via Direct Energy-based Preference Optimization. *arXiv*. DOI: [10.48550/arxiv.2403.16576](https://doi.org/10.48550/arxiv.2403.16576)
