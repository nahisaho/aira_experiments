# Improving Cross-Ancestry Polygenic Risk Score Transferability via Bayesian LD Correction, Multi-Ancestry Meta-Analysis, and Local Ancestry Integration

## Abstract

Polygenic risk scores (PRS) derived from genome-wide association studies (GWAS) in European-ancestry populations exhibit substantially reduced predictive accuracy when applied to non-European populations, a critical barrier to equitable precision medicine. This study develops and evaluates an integrated statistical framework for improving PRS transferability from European (EUR) to East Asian (EAS) populations, motivated by the UK Biobank-to-BioBank Japan transfer problem. We propose three complementary approaches: (1) a Bayesian linkage disequilibrium (LD) correction method that projects EUR GWAS effect estimates into EAS LD space using continuous shrinkage priors, (2) a multi-ancestry fixed-effect meta-analysis for re-estimating SNP effects across populations, and (3) a local ancestry-informed PRS weighting scheme. Through comprehensive simulation experiments modeling type 2 diabetes (T2D) genetic architecture with varying population divergence (Fst = 0.01–0.20), sample sizes (N = 500–10,000), and heritabilities (h² = 0.1–0.7), we demonstrate that the combined approach improves prediction accuracy by 2.8% in AUC and 19.9% in liability-scale R² compared to direct EUR-to-EAS PRS transfer. Multi-ancestry meta-analysis provides the largest single-method improvement, while Bayesian LD correction offers complementary gains. These results highlight the importance of accounting for population-specific LD structure and leveraging multi-ancestry GWAS data for equitable genetic risk prediction. Our open-source simulation framework enables systematic evaluation of cross-ancestry PRS methods under controlled conditions.

## 1. Introduction

Polygenic risk scores (PRS) aggregate the effects of many genetic variants to predict individual susceptibility to complex diseases and traits (Khera et al., 2018). As GWAS sample sizes have grown, PRS have shown increasing clinical promise for risk stratification in conditions such as cardiovascular disease, type 2 diabetes, and cancer. However, the overwhelming majority of GWAS participants are of European ancestry (Sirugo et al., 2019), creating a fundamental equity problem: PRS derived from European GWAS perform significantly worse in non-European populations (Martin et al., 2019).

This "portability gap" arises from multiple sources. First, linkage disequilibrium (LD) patterns differ across populations due to distinct demographic histories, causing EUR-optimized SNP weights to be miscalibrated in other populations (Privé et al., 2022). Second, allele frequencies diverge between populations (quantified by Fst), altering the contribution of each variant to phenotypic variance. Third, causal effect sizes may vary due to gene-environment interactions or epistatic effects that differ across ancestries.

Recent methodological advances have begun to address this challenge. PRS-CSx (Ruan et al., 2022) employs a Bayesian framework that jointly models GWAS summary statistics across ancestries with continuous shrinkage priors coupled across populations. BridgePRS (Hoggart et al., 2024) leverages shared genetic effects across ancestries using a Bayesian ridge regression framework. TL-PRS (Zhao et al., 2022) applies transfer learning to fine-tune EUR-based models for target populations. CT-SLEB (Zhao et al., 2023) combines clumping-thresholding with empirical Bayes and superlearning for multi-ancestry prediction.

Despite these advances, several gaps remain. First, few methods explicitly model the transformation between population-specific LD structures when projecting effect sizes. Second, the integration of local ancestry information with cross-ancestry PRS methods has been underexplored. Third, systematic evaluation across realistic parameter spaces (Fst, sample size, heritability) is limited.

In this study, we make the following contributions:

1. We formalize the EUR-to-EAS PRS transfer problem as an LD-space projection and develop a Bayesian correction method with continuous shrinkage priors.
2. We implement a multi-ancestry meta-analysis pipeline for SNP effect re-estimation across populations.
3. We propose a local ancestry-informed PRS correction model that applies population-specific weights based on genomic segment ancestry.
4. We integrate these three approaches into a combined framework and evaluate it through comprehensive simulation experiments.
5. We demonstrate the framework using type 2 diabetes as a case study, with realistic genetic architecture parameters.

## 2. Related Work

### 2.1 PRS Portability Problem

Martin et al. (2019) demonstrated that PRS accuracy decreases with increasing genetic distance from the discovery population, with European-derived PRS explaining only ~25% of the variance in East Asian populations compared to European populations. Privé et al. (2022) systematically evaluated the portability of 245 polygenic scores across 9 ancestry groups within the UK Biobank, confirming pervasive accuracy loss and identifying LD mismatch as a primary driver.

### 2.2 Cross-Ancestry PRS Methods

**PRS-CSx** (Ruan et al., 2022) extends the PRS-CS framework to multiple ancestries by placing a shared continuous shrinkage prior on SNP effects across populations, coupled through population-specific LD reference panels. Applied to 17 quantitative traits and 5 diseases across four ancestry groups, PRS-CSx demonstrated 50–70% relative improvement in prediction accuracy for non-European populations compared to single-ancestry PRS.

**BridgePRS** (Hoggart et al., 2024) uses a two-stage Bayesian approach: first estimating ancestry-specific effects via ridge regression, then combining them using shared effect priors. BridgePRS showed 61% larger average R² in African ancestry samples compared to PRS-CSx when evaluated in the BioMe Biobank.

**TL-PRS** (Zhao et al., 2022) applies transfer learning by fine-tuning European-derived PRS models for target populations using gradient descent optimization. This approach achieved 25–29% relative improvement in South Asian and African ancestry UK Biobank participants.

**CT-SLEB** (Zhao et al., 2023) combines clumping and thresholding, empirical Bayes shrinkage, and super learning to construct multi-ancestry PRS. Evaluated on 5 million individuals across 13 traits, CT-SLEB showed improved accuracy across all ancestry groups, with particular gains in underrepresented populations.

### 2.3 Local Ancestry Approaches

Local ancestry inference identifies the ancestral origin of specific genomic segments in admixed individuals. Several studies have shown that incorporating local ancestry can improve PRS accuracy by applying ancestry-appropriate effect sizes at each locus (Marnetto et al., 2020). However, integration of local ancestry with cross-ancestry Bayesian methods remains limited.

## 3. Methods

### 3.1 Problem Formulation

Consider two populations: a source population (EUR, indexed by subscript $s$) and a target population (EAS, indexed by subscript $t$). For $p$ SNPs, let:

- $\hat{\boldsymbol{\beta}}_s$: GWAS marginal effect estimates from the source population
- $\mathbf{R}_s, \mathbf{R}_t$: LD correlation matrices for source and target populations
- $n_s, n_t$: sample sizes

The standard PRS for an individual with genotype vector $\mathbf{g}$ is:

$$\text{PRS} = \mathbf{g}^\top \hat{\boldsymbol{\beta}}_s$$

The transferability problem arises because $\hat{\boldsymbol{\beta}}_s$ is estimated under LD structure $\mathbf{R}_s$, but applied to genotypes generated under $\mathbf{R}_t$.

### 3.2 Bayesian LD Correction

We develop a Bayesian estimator that projects source effect estimates into the target LD space. Under a Gaussian prior with global shrinkage parameter $\phi$:

$$\boldsymbol{\beta} \sim \mathcal{N}(\mathbf{0}, \phi \mathbf{I})$$

The posterior mean in the source LD space is:

$$\hat{\boldsymbol{\beta}}_{\text{post}} = \left(n_s \tilde{\mathbf{R}}_s + \phi^{-1} \mathbf{I}\right)^{-1} n_s \tilde{\mathbf{R}}_s \hat{\boldsymbol{\beta}}_s$$

where $\tilde{\mathbf{R}}_s = \mathbf{R}_s + \phi \mathbf{I}$ is the regularized LD matrix.

The LD-adjusted estimate for the target population is obtained via the transformation:

$$\hat{\boldsymbol{\beta}}_{\text{adj}} = \tilde{\mathbf{R}}_t \tilde{\mathbf{R}}_s^{-1} \hat{\boldsymbol{\beta}}_{\text{post}}$$

This transformation effectively "deconvolves" the source LD structure and "reconvolves" with the target LD structure, accounting for differences in correlation patterns between populations.

### 3.3 Multi-Ancestry Meta-Analysis

For SNP $j$, we combine effect estimates across populations using inverse-variance weighted fixed-effect meta-analysis:

$$\hat{\beta}_{j,\text{meta}} = \frac{w_{j,s} \hat{\beta}_{j,s} + w_{j,t} \hat{\beta}_{j,t}}{w_{j,s} + w_{j,t}}$$

where $w_{j,k} = 1/\text{SE}_{j,k}^2$ are precision weights. The meta-analysis standard error is:

$$\text{SE}_{j,\text{meta}} = \sqrt{\frac{1}{w_{j,s} + w_{j,t}}}$$

We also implement a random-effects variant using the DerSimonian-Laird estimator to account for effect heterogeneity across populations.

### 3.4 Local Ancestry-Informed PRS

For individuals with known local ancestry at each locus, we compute a hybrid PRS:

$$\text{PRS}_{\text{LA}} = \sum_{j=1}^{p} g_j \cdot \left[ \ell_j \hat{\beta}_{j,s} + (1 - \ell_j) \hat{\beta}_{j,\text{adj}} \right]$$

where $\ell_j \in \{0, 1\}$ indicates local EUR ancestry at locus $j$, $\hat{\beta}_{j,s}$ is the EUR effect estimate, and $\hat{\beta}_{j,\text{adj}}$ is the LD-adjusted estimate.

### 3.5 Combined Framework

The proposed combined method integrates all three approaches:

1. Compute Bayesian LD-corrected effects $\hat{\boldsymbol{\beta}}_{\text{adj}}$
2. Compute meta-analysis effects $\hat{\boldsymbol{\beta}}_{\text{meta}}$
3. Form combined weights: $\hat{\boldsymbol{\beta}}_{\text{comb}} = \alpha \hat{\boldsymbol{\beta}}_{\text{meta}} + (1-\alpha) \hat{\boldsymbol{\beta}}_{\text{adj}}$ with $\alpha = 0.5$
4. Apply local ancestry correction using combined weights

### 3.6 Simulation Design

**Population genetics model:**
- $p = 500$ SNPs with $p_c = 50$ causal variants
- Allele frequencies simulated via Balding-Nichols model with varying Fst
- Block-diagonal LD matrices with exponential decay (EUR: $\rho = 0.85$; EAS: $\rho = 0.75$)
- Genotypes simulated via Cholesky decomposition preserving LD structure

**Phenotype model:**
- Continuous liability: $y = \mathbf{g}^\top \boldsymbol{\beta} + \epsilon$, $\epsilon \sim \mathcal{N}(0, \sigma_e^2)$
- Binary trait via liability threshold model with prevalence $K = 0.10$ (T2D)
- Cross-population effect correlation $\sim 0.8$ (partial sharing)

**Parameter grid:**
- Fst $\in \{0.01, 0.05, 0.10, 0.15, 0.20\}$
- $n_{\text{EAS}} \in \{500, 1000, 2000, 5000, 10000\}$
- $h^2 \in \{0.1, 0.2, 0.3, 0.5, 0.7\}$

## 4. Experiments

### 4.1 Experimental Setup

All simulations were implemented in Python using NumPy, SciPy, and scikit-learn. Random seed was fixed at 42 for reproducibility. The primary simulation parameters were: $p = 500$ SNPs, $p_c = 50$ causal variants, $n_{\text{EUR}} = 10{,}000$, $n_{\text{EAS}} = 5{,}000$, Fst = 0.10, $h^2 = 0.50$, and disease prevalence $K = 0.10$.

### 4.2 Evaluation Metrics

- **AUC (Area Under the ROC Curve):** Discrimination accuracy for binary trait classification
- **R² (liability scale):** Squared correlation between PRS and underlying genetic liability
- **R² (observed scale):** Squared correlation between PRS and observed binary phenotype

### 4.3 Baseline Methods

Six methods were compared:
1. **Direct Transfer:** EUR C+T PRS applied directly to EAS genotypes (baseline)
2. **Target Population:** EAS-specific C+T PRS (oracle upper bound)
3. **Bayesian LD Correction:** Our LD-space projection method
4. **Multi-Ancestry Meta-Analysis:** Fixed-effect meta-analysis of EUR and EAS GWAS
5. **Local Ancestry PRS:** Ancestry-specific weights based on local ancestry
6. **Combined (Proposed):** Integration of methods 3–5

### 4.4 Parameter Sweep Experiments

Three systematic parameter sweeps were conducted:
- **Fst sweep:** Population divergence from 0.01 to 0.20 (5 levels)
- **Sample size sweep:** Target population N from 500 to 10,000 (5 levels)
- **Heritability sweep:** Trait heritability from 0.1 to 0.7 (5 levels)

## 5. Results

### 5.1 Main Simulation Results

Under default parameters (Fst = 0.10, $n_{\text{EUR}} = 10{,}000$, $n_{\text{EAS}} = 5{,}000$, $h^2 = 0.50$), the proposed combined method achieved AUC = 0.8135 and R²(liability) = 0.3582, representing improvements of +2.8% and +19.9% over direct transfer (AUC = 0.7914, R² = 0.2987), respectively.

| Method | AUC | R²(liability) | R²(observed) |
|--------|-----|---------------|--------------|
| Direct Transfer (EUR→EAS) | 0.7914 | 0.2987 | 0.1009 |
| Target Pop (EAS GWAS) | 0.8335 | 0.3911 | 0.1319 |
| Bayesian LD Correction | 0.7980 | 0.3288 | 0.1045 |
| Multi-Ancestry Meta | 0.8131 | 0.3463 | 0.1159 |
| Local Ancestry PRS | 0.7927 | 0.3145 | 0.1003 |
| **Combined (Proposed)** | **0.8135** | **0.3582** | **0.1167** |

The multi-ancestry meta-analysis provided the largest single-method improvement in AUC (+0.0217), while the Bayesian LD correction offered the second-largest gain in R²(liability) (+0.0301). The combined method captured complementary information from both approaches.

![Figure 1: Method comparison showing AUC and R² for all six PRS methods](figures/method_comparison.png)

### 5.2 PRS Distribution Analysis

Figure 2 shows PRS distributions stratified by case/control status for each method. The combined method exhibits the greatest separation between cases and controls, approaching (but not reaching) the target population oracle.

![Figure 2: PRS distributions by T2D case/control status for each method](figures/prs_distributions.png)

### 5.3 LD Structure Differences

The simulated EUR and EAS LD matrices exhibit distinct block structures reflecting different decay rates. The difference matrix reveals systematic deviations that motivate LD-aware PRS correction methods.

![Figure 3: LD matrix comparison between EUR and EAS populations](figures/ld_comparison.png)

### 5.4 Effect Size Concordance

Cross-population effect size correlation was r = 0.458 for GWAS marginal estimates, reflecting both true effect sharing and noise. Causal variants showed stronger concordance than non-causal variants, consistent with shared genetic architecture with population-specific LD confounding.

![Figure 4: Effect size comparison between EUR and EAS GWAS estimates](figures/effect_sizes.png)

### 5.5 Allele Frequency Divergence

At Fst = 0.10, allele frequencies between populations showed moderate divergence (r ≈ 0.87), with individual SNPs exhibiting frequency differences up to 0.3.

![Figure 5: Allele frequency divergence between EUR and EAS populations](figures/allele_freq_divergence.png)

### 5.6 Impact of Population Divergence (Fst)

Increasing Fst from 0.01 to 0.20 reduced direct transfer AUC from 0.7954 to 0.7732, while the combined method maintained higher accuracy across all divergence levels (AUC 0.7983 at Fst = 0.20, ΔAUC = +0.0251).

![Figure 6: PRS method performance across population divergence levels (Fst)](figures/fst_sweep.png)

### 5.7 Impact of Target Sample Size

The combined method showed consistent improvement over direct transfer across all target sample sizes, with the largest absolute gains at moderate sample sizes (N_EAS = 2,000: ΔAUC = +0.0304).

![Figure 7: PRS method performance across target population sample sizes](figures/sample_size_sweep.png)

### 5.8 Impact of Heritability

All methods showed monotonically increasing performance with heritability. The combined method's advantage was most pronounced at moderate heritabilities (h² = 0.2–0.5), consistent with the regime where LD correction and meta-analysis provide the most value.

![Figure 8: PRS method performance across heritability values](figures/heritability_sweep.png)

## 6. Discussion

### 6.1 Key Findings

Our simulation study demonstrates that integrating Bayesian LD correction, multi-ancestry meta-analysis, and local ancestry information can substantially improve PRS transferability from European to East Asian populations. The combined approach achieved a 19.9% improvement in liability-scale R² compared to naive direct transfer, closing approximately 63% of the gap between direct transfer and the target-population oracle.

The multi-ancestry meta-analysis was the single most impactful component, consistent with findings from PRS-CSx (Ruan et al., 2022) and CT-SLEB (Zhao et al., 2023) that leveraging multi-ancestry GWAS data is critical. The Bayesian LD correction provided complementary gains by explicitly modeling the transformation between population-specific LD structures, an approach conceptually aligned with BridgePRS (Hoggart et al., 2024) but using a direct LD-space projection rather than ridge regression.

### 6.2 Comparison with Existing Methods

Our Bayesian LD correction shares conceptual foundations with PRS-CSx but uses an explicit LD transformation matrix rather than coupled shrinkage priors. The meta-analysis component is simpler than CT-SLEB's super-learning framework but provides a strong baseline that outperforms direct transfer. The local ancestry integration is analogous to approaches used in admixed populations but applied here to correct for residual LD mismatch even in non-admixed individuals.

### 6.3 Implications for Type 2 Diabetes

T2D has a heritability of approximately 0.3–0.5 and prevalence differences across populations. Our simulations at h² = 0.3–0.5 and prevalence K = 0.10 suggest that the combined approach could provide clinically meaningful improvements in risk stratification for Japanese and East Asian populations using existing European GWAS data augmented with smaller target-population studies.

### 6.4 Limitations

Several limitations should be noted. First, our simulation uses 500 SNPs, far fewer than genome-wide analyses involving millions of variants. The computational tractability of the LD matrix inversion approach would require block-wise or low-rank approximations at genomic scale. Second, we model only two discrete populations without continuous admixture structure. Third, environmental factors and gene-environment interactions are not modeled. Fourth, the fixed mixing weight ($\alpha = 0.5$) in the combined method could be optimized through cross-validation. Fifth, our simulations assume known LD matrices, whereas in practice these must be estimated from reference panels.

### 6.5 Future Directions

Future work should extend this framework to: (1) genome-wide scale using sparse LD representations, (2) multiple target populations simultaneously, (3) adaptive mixing weights learned from validation data, (4) integration with functional genomic annotations for improved variant prioritization, and (5) application to real GWAS summary statistics from UK Biobank and BioBank Japan.

## 7. Conclusion

We developed an integrated statistical framework for improving cross-ancestry PRS transferability, combining Bayesian LD correction, multi-ancestry meta-analysis, and local ancestry-informed weighting. Through comprehensive simulation experiments modeling type 2 diabetes genetic architecture, we demonstrated consistent improvements over direct EUR-to-EAS PRS transfer across varying levels of population divergence, sample size, and heritability. The proposed framework provides a principled and extensible approach to addressing the PRS portability problem, contributing to the goal of equitable precision medicine across diverse populations.

## References

1. Martin, A.R., Kanai, M., Kamatani, Y., Okada, Y., Neale, B.M., & Daly, M.J. (2019). Clinical use of current polygenic risk scores may exacerbate health disparities. *Nature Genetics*, 51, 584–591. https://doi.org/10.1038/s41588-019-0379-x

2. Ruan, Y., Lin, Y.F., Feng, Y.C.A., et al. (2022). Improving polygenic prediction in ancestrally diverse populations. *Nature Genetics*, 54(5), 573–580. https://doi.org/10.1038/s41588-021-01054-8

3. Privé, F., Aschard, H., Carmi, S., Folkersen, L., Hoggart, C., O'Reilly, P.F., & Vilhjálmsson, B.J. (2022). Portability of 245 polygenic scores when derived from the UK Biobank and applied to 9 ancestry groups from the same cohort. *American Journal of Human Genetics*, 109(1), 12–23. https://doi.org/10.1016/j.ajhg.2021.11.008

4. Zhao, Z., Fritsche, L.G., Smith, J.A., et al. (2022). The construction of cross-population polygenic risk scores using transfer learning. *American Journal of Human Genetics*, 109(4), 716–729. https://doi.org/10.1016/j.ajhg.2022.02.011

5. Zhao, Z., Zheng, W., Kraft, P., et al. (2023). A new method for multiancestry polygenic prediction improves performance across diverse populations. *Nature Genetics*, 55, 1757–1768. https://doi.org/10.1038/s41588-023-01501-z

6. Hoggart, C.J., Choi, S.W., García-González, J., Souaiaia, T., Preuss, M., & O'Reilly, P. (2024). BridgePRS leverages shared genetic effects across ancestries to increase polygenic risk score portability. *Nature Genetics*, 56(1), 180–186. https://doi.org/10.1038/s41588-023-01583-9

7. Sirugo, G., Williams, S.M., & Tishkoff, S.A. (2019). The missing diversity in human genetic studies. *Cell*, 177(1), 26–31. https://doi.org/10.1016/j.cell.2019.02.048

8. Khera, A.V., Chaffin, M., Aragam, K.G., et al. (2018). Genome-wide polygenic scores for common diseases identify individuals with risk equivalent to monogenic mutations. *Nature Genetics*, 50(9), 1219–1224. https://doi.org/10.1038/s41588-018-0183-z
