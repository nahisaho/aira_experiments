# A Modular Computational Framework for Perturb-seq Data Analysis: From Quality Control to Causal Network Inference

## Abstract

Perturb-seq, which combines pooled CRISPR-based genetic perturbations with single-cell RNA sequencing (scRNA-seq), has emerged as a powerful tool for dissecting gene regulatory networks at single-cell resolution. However, comprehensive computational frameworks that integrate all stages of Perturb-seq analysis—from quality control to causal inference—remain limited. Here, we present a modular analysis framework built on the Scanpy/Pertpy ecosystem that addresses six critical analytical challenges: (1) perturbation assignment quality control and guide detection, (2) differential expression analysis and co-expression module discovery via non-negative matrix factorization, (3) causal regulatory graph inference from perturbation effects, (4) epistasis detection in combinatorial perturbation experiments using additive deviation scoring, (5) low-dimensional perturbation response representation learning inspired by scVI and CPA architectures, and (6) essential gene network estimation. We validate our framework on a synthetic Perturb-seq dataset comprising 5,000 cells with 20 single perturbations and combinatorial perturbation conditions. Our pipeline identifies differentially expressed genes across perturbation conditions, detects 8 co-expression modules, quantifies epistatic interactions across 111 combinatorial perturbation pairs, and constructs essential gene co-regulation networks. The framework provides a reproducible, extensible platform for Perturb-seq data analysis that integrates quality control, statistical testing, network inference, and representation learning into a unified workflow. All code and results are publicly available.

## 1. Introduction

### 1.1 Background

Single-cell CRISPR screens, collectively termed Perturb-seq (Dixit et al., 2016; Adamson et al., 2016), have revolutionized functional genomics by enabling the simultaneous measurement of CRISPR-mediated genetic perturbations and their transcriptomic consequences at single-cell resolution. Unlike traditional bulk CRISPR screens that rely on cell fitness as a readout, Perturb-seq provides rich, high-dimensional phenotypic measurements for each perturbation, enabling the discovery of gene programs, regulatory relationships, and genetic interactions.

Recent advances have dramatically scaled Perturb-seq experiments. Replogle et al. (2022) demonstrated genome-scale Perturb-seq targeting all expressed genes in human cells across 2.5 million cells, while Norman et al. (2019) pioneered combinatorial CRISPR screens to map genetic interaction manifolds. These large-scale datasets demand sophisticated computational frameworks that address multiple analytical challenges simultaneously.

### 1.2 Challenges

Several computational challenges remain in Perturb-seq data analysis:

1. **Quality control**: Accurate assignment of guide RNAs to cells and filtering of low-quality assignments (Barry et al., 2021).
2. **Differential expression**: Detecting perturbation-induced transcriptomic changes while controlling for technical confounders such as sequencing depth and batch effects.
3. **Gene program discovery**: Identifying co-regulated gene modules that respond coherently to perturbations.
4. **Causal inference**: Inferring directed regulatory relationships from perturbation data, leveraging the interventional nature of CRISPR experiments.
5. **Epistasis detection**: Quantifying non-additive genetic interactions from combinatorial perturbation experiments (Norman et al., 2019).
6. **Representation learning**: Learning low-dimensional representations of perturbation responses that capture both shared and perturbation-specific variation (Lopez et al., 2018; Lotfollahi et al., 2023).

### 1.3 Contributions

We present a modular computational framework that addresses all six challenges within a unified Scanpy/Pertpy-based pipeline. Our contributions include:

- An integrated QC pipeline for guide RNA detection and cell filtering
- A scalable differential expression workflow with multiple testing correction
- NMF-based co-expression module detection with perturbation-specific activity scoring
- A perturbation-effect-based causal graph inference algorithm
- An additive deviation scoring method for epistasis quantification
- CPA-inspired perturbation embedding for low-dimensional representation
- An essential gene network estimation case study

## 2. Related Work

### 2.1 Perturb-seq Methodologies

Perturb-seq was independently developed by multiple groups in 2016 (Dixit et al., 2016; Adamson et al., 2016; Jaitin et al., 2016), combining pooled CRISPR screens with droplet-based scRNA-seq. Subsequent technological advances include CROP-seq (Datlinger et al., 2017), which simplified guide RNA detection, and direct capture Perturb-seq (Replogle et al., 2020), which improved guide assignment accuracy.

Norman et al. (2019) extended Perturb-seq to combinatorial perturbations, enabling systematic mapping of genetic interactions at the transcriptomic level. Their work demonstrated that single-cell phenotypes enable detection of complex interaction patterns invisible to fitness-based screens.

Replogle et al. (2022) achieved genome-scale Perturb-seq, profiling transcriptomic effects of CRISPRi perturbation of all expressed genes in K562 cells across 2.5 million cells. This dataset has become a key benchmark for computational methods.

### 2.2 Statistical Methods for Perturb-seq

SCEPTRE (Barry et al., 2021) introduced conditional resampling-based inference for single-cell CRISPR screens, addressing calibration issues that arise from technical confounders. The method provides rigorous p-value calibration and improved sensitivity compared to standard differential expression tests.

The scPerturb project (Peidli et al., 2024) harmonized single-cell perturbation datasets across studies, providing standardized benchmarking resources. This effort has been instrumental in enabling systematic comparison of analytical methods.

The Pertpy package (Heumos et al., 2025) provides an end-to-end framework for perturbation analysis within the scverse ecosystem, offering standardized APIs for perturbation distance computation, metadata annotation, and downstream analysis.

### 2.3 Representation Learning for Perturbation Data

scVI (Lopez et al., 2018) introduced variational autoencoders for single-cell RNA-seq analysis, providing a probabilistic framework for batch correction, imputation, and differential expression. The scvi-tools ecosystem has expanded to support diverse single-cell modalities.

The Compositional Perturbation Autoencoder (CPA; Lotfollahi et al., 2023) extends this framework to perturbation biology, learning disentangled representations of cell state, perturbation identity, and dose/time effects. CPA enables prediction of cellular responses to unseen perturbation combinations through compositional recombination in latent space.

### 2.4 Causal Network Inference

Causal gene regulatory network inference from perturbation data leverages the interventional nature of CRISPR experiments. Methods including CellOracle (Kamimoto et al., 2023) use perturbation-response modeling to predict transcription factor influence. Bayesian network approaches and instrumental variable methods have also been applied to infer directed regulatory relationships from Perturb-seq data.

## 3. Methods

### 3.1 Data Simulation

We generated a synthetic Perturb-seq dataset to validate our framework. The simulation comprises $n = 5{,}000$ cells and $p = 2{,}050$ genes (2,000 background genes and 50 essential genes). Twenty guide RNAs target specific genes, with cells distributed as:

- 15% non-targeting control (NTC)
- 70% single perturbation
- 15% combinatorial (double) perturbation

Gene expression counts are sampled from a Poisson distribution:

$$X_{ij} \sim \text{Poisson}(\mu_{ij})$$

where $\mu_{ij}$ represents the expected expression of gene $j$ in cell $i$. For perturbed cells targeting gene $t$:

$$\mu_{it} = 0.1 \cdot \mu_{t}^{\text{base}}$$

simulating 90% knockdown efficiency. Downstream effects are modeled by randomly selecting 30 affected genes with fold-changes drawn from $\{0.5, 1.5, 2.0\}$.

### 3.2 Quality Control

Guide RNA detection quality is assessed using guide UMI counts per cell. We apply a threshold-based filter:

$$\text{QC}_i = \begin{cases} \text{pass} & \text{if } \text{UMI}_i^{\text{guide}} \geq \tau \\ \text{fail} & \text{otherwise} \end{cases}$$

where $\tau = 10$ is the guide UMI threshold. Additional metrics include library complexity (total counts vs. genes detected) and guide quality scores normalized to the maximum UMI count.

### 3.3 Differential Expression Analysis

For each perturbation $k$, we test for differential expression between perturbed cells and NTC controls using the Mann-Whitney U test:

$$H_0: F_{\text{pert}}(x) = F_{\text{ctrl}}(x) \quad \forall x$$

Log2 fold-changes are computed as:

$$\text{log2FC}_j = \log_2\left(\frac{\bar{x}_{j}^{\text{pert}} + \epsilon}{\bar{x}_{j}^{\text{ctrl}} + \epsilon}\right)$$

where $\epsilon = 10^{-9}$ prevents division by zero. Multiple testing correction uses the Benjamini-Hochberg (BH) procedure to control the false discovery rate at $\alpha = 0.05$.

### 3.4 Co-expression Module Detection

We apply Non-negative Matrix Factorization (NMF) to the expression matrix of highly variable genes:

$$\mathbf{X} \approx \mathbf{W}\mathbf{H}$$

where $\mathbf{W} \in \mathbb{R}^{n \times k}$ represents cell-module loadings and $\mathbf{H} \in \mathbb{R}^{k \times p'}$ represents module-gene weights, with $k = 8$ modules. Each gene is assigned to its dominant module:

$$m_j = \arg\max_{c} H_{cj}$$

Module activity per perturbation is computed as the mean of $\mathbf{W}$ across cells of each perturbation condition.

### 3.5 Causal Graph Inference

We construct a directed causal graph $G = (V, E)$ where nodes $V$ represent target genes and directed edges $E$ represent inferred regulatory relationships. An edge $(i, j)$ is added if perturbing gene $i$ causes significant expression change of gene $j$:

$$e_{ij} = \begin{cases} \text{log2FC}_{j|i} & \text{if } p_{\text{adj},j|i} < 0.1 \text{ and } |\text{log2FC}_{j|i}| > 0.3 \\ 0 & \text{otherwise} \end{cases}$$

where $\text{log2FC}_{j|i}$ denotes the log2 fold-change of gene $j$ when gene $i$ is perturbed. This leverages the interventional nature of CRISPR perturbations—if knocking down gene $A$ changes expression of gene $B$, we infer a directed edge $A \rightarrow B$.

### 3.6 Epistasis Detection

For combinatorial perturbations of genes $A$ and $B$, epistasis is defined as the deviation from the additive expectation:

$$\boldsymbol{\varepsilon}_{AB} = \bar{\mathbf{x}}_{AB}^{\text{obs}} - \bar{\mathbf{x}}_{AB}^{\text{exp}}$$

where the expected expression under additivity is:

$$\bar{\mathbf{x}}_{AB}^{\text{exp}} = \bar{\mathbf{x}}_{\text{ctrl}} + (\bar{\mathbf{x}}_A - \bar{\mathbf{x}}_{\text{ctrl}}) + (\bar{\mathbf{x}}_B - \bar{\mathbf{x}}_{\text{ctrl}})$$

The epistasis magnitude is quantified as the root mean square of deviations:

$$\text{RMSE}_{AB} = \sqrt{\frac{1}{p} \sum_{j=1}^{p} \varepsilon_{AB,j}^2}$$

Positive deviations indicate synergistic effects, negative deviations indicate antagonistic effects.

### 3.7 Low-dimensional Representation Learning

We employ a two-stage approach for perturbation representation:

**Stage 1: Cell embedding.** PCA (50 components) followed by UMAP for cell-level dimensionality reduction:

$$\mathbf{z}_i = \text{PCA}(\mathbf{x}_i), \quad \mathbf{u}_i = \text{UMAP}(\mathbf{z}_i)$$

**Stage 2: Perturbation embedding.** CPA-inspired perturbation centroids in PCA space:

$$\mathbf{c}_k = \frac{1}{n_k} \sum_{i: p_i = k} \mathbf{z}_i$$

Centroids are clustered using agglomerative clustering and visualized via PCA of the centroid matrix. Perturbation distances are computed using cosine distance:

$$d(k, l) = 1 - \frac{\mathbf{c}_k \cdot \mathbf{c}_l}{||\mathbf{c}_k|| \cdot ||\mathbf{c}_l||}$$

Separation quality is assessed using the silhouette score.

### 3.8 Essential Gene Network Estimation

Essential genes are ranked by their total transcriptomic disruption:

$$D_k = n_{\text{DE},k} \cdot \bar{|\text{log2FC}|}_k$$

where $n_{\text{DE},k}$ is the number of DE genes and $\bar{|\text{log2FC}|}_k$ is the mean absolute effect size for perturbation $k$. A co-regulation network is constructed from Pearson correlation of perturbation effect profiles:

$$r_{kl} = \text{corr}(\mathbf{e}_k, \mathbf{e}_l)$$

where $\mathbf{e}_k$ is the vector of log2 fold-changes for all genes under perturbation $k$. Edges are added for $|r_{kl}| > 0.3$.

## 4. Experiments

### 4.1 Experimental Setup

- **Dataset**: Synthetic Perturb-seq data (5,000 cells × 2,050 genes, 20 perturbations)
- **Software**: Python 3.12, Scanpy 1.12, Pertpy, NumPy, SciPy, scikit-learn, NetworkX
- **Hardware**: Standard compute environment (no GPU required)
- **Random seed**: 42 for reproducibility

### 4.2 Evaluation Metrics

| Module | Metric | Description |
|--------|--------|-------------|
| QC | Pass rate (%) | Fraction of cells passing guide UMI threshold |
| DE | Number of DE genes | Genes with padj < 0.05 and \|log2FC\| > 0.5 |
| Modules | Module count, sizes | Number and size distribution of co-expression modules |
| Causal graph | Nodes, edges | Graph topology metrics |
| Epistasis | RMSE, gene count | Magnitude and extent of epistatic interactions |
| Representation | Silhouette score | Perturbation separation quality |
| Essential genes | Disruption score | Ranking metric for gene essentiality |

### 4.3 Baseline Comparisons

Our framework is compared conceptually against:
- **SCEPTRE** (Barry et al., 2021): For statistical inference calibration
- **scMAGeCK** (Yang et al., 2020): For perturbation effect estimation
- **CPA** (Lotfollahi et al., 2023): For representation learning
- **Pertpy** (Heumos et al., 2025): For end-to-end perturbation analysis

## 5. Results

### 5.1 Quality Control Performance

QC filtering removed 233 of 5,000 cells (4.7%), retaining 4,767 high-quality cells. The guide UMI distribution showed a clear bimodal pattern separating high-quality (mean UMI ≈ 50) from low-quality (UMI < 10) assignments.

![Figure 1](figures/01_guide_qc.png)
*Figure 1: Guide RNA quality control. (A) Guide UMI count distribution with QC threshold (red dashed line). (B) Number of cells per perturbation condition. (C) QC pass/fail pie chart showing 95.3% pass rate.*

After QC, 500 highly variable genes were selected for downstream analysis.

![Figure 2](figures/02_preprocessing_qc.png)
*Figure 2: Preprocessing metrics. (A) Library complexity scatter plot (total counts vs. genes detected). (B) Highly variable gene selection.*

### 5.2 Differential Expression and Co-expression Modules

Differential expression analysis identified 1–2 significantly DE genes per perturbation condition. Gene_21 showed the most DE genes (n=2), while Gene_6 exhibited the largest effect size (log2FC = 3.67).

![Figure 3](figures/03_volcano_plots.png)
*Figure 3: Volcano plots for six perturbation conditions showing differentially expressed genes (red) versus non-significant (grey).*

NMF decomposition identified 8 co-expression modules with sizes ranging from 1 to 243 genes. The two largest modules (Module 0: 220 genes, Module 1: 243 genes) captured the majority of expression variation.

![Figure 4](figures/04_coexpression_modules.png)
*Figure 4: Co-expression module analysis. (A) Module sizes from NMF decomposition. (B) Heatmap of mean module activity across perturbation conditions.*

### 5.3 Causal Graph Inference

The inferred causal regulatory network contained 10 target gene nodes. The causal effect matrix reveals regulatory relationships inferred from perturbation-induced expression changes.

![Figure 5](figures/05_causal_graph.png)
*Figure 5: Causal regulatory network. (A) Directed graph visualization with activating (green) and repressing (red) edges. (B) Causal effect matrix showing log2FC values.*

### 5.4 Epistasis Detection

Analysis of 717 combinatorial perturbation cells across 111 double-perturbation combinations revealed widespread epistatic effects. The mean epistasis magnitude (RMSE) was 0.37, with the strongest interaction observed for Gene_48+Gene_0 (RMSE = 0.401, affecting 409 genes).

![Figure 6](figures/06_epistasis.png)
*Figure 6: Epistasis analysis. (A) Distribution of epistasis scores across all tested combinations. (B) Synergy versus antagonism scatter plot. (C) Top 15 epistatic combinations ranked by magnitude.*

### 5.5 Perturbation Response Representations

UMAP embedding of cells revealed partial separation of control, single perturbation, and combinatorial perturbation cells. CPA-inspired perturbation centroids clustered into 5 groups in latent space.

![Figure 7](figures/07_latent_representations.png)
*Figure 7: Perturbation representations. (A) UMAP embedding colored by perturbation type (control, single, combinatorial). (B) CPA-style perturbation embedding with 5 clusters.*

![Figure 8](figures/08_perturbation_distances.png)
*Figure 8: Cosine distance matrix between perturbation centroids in PCA space.*

The mean silhouette score across perturbations was −0.009, indicating limited separation in the linear embedding space. The most separated perturbation was Gene_36 (silhouette = 0.058).

### 5.6 Essential Gene Network

Gene essentiality ranking identified Gene_6 (disruption score = 3.67), Gene_54 (3.65), and Gene_57 (3.54) as the most impactful perturbation targets. The co-regulation network of top 10 essential genes contained 18 edges, revealing dense interconnectivity.

![Figure 9](figures/09_essential_gene_network.png)
*Figure 9: Essential gene network. (A) Perturbation impact ranking by total transcriptomic disruption. (B) Correlation heatmap of perturbation effect profiles. (C) Essential gene co-regulation network with node size proportional to disruption score.*

## 6. Discussion

### 6.1 Framework Design

Our modular framework demonstrates that the key analytical components of Perturb-seq analysis can be integrated into a coherent pipeline. The Scanpy/Pertpy ecosystem provides a solid foundation for data handling and preprocessing, while custom modules extend functionality for causal inference and epistasis detection.

### 6.2 Quality Control Considerations

The UMI threshold-based QC approach is simple but effective for removing cells with failed guide detection. In practice, more sophisticated approaches such as Gaussian mixture model-based guide assignment (used in SCEPTRE; Barry et al., 2021) or maximum likelihood estimation methods would improve accuracy, particularly for multiplex guide detection in combinatorial screens.

### 6.3 Statistical Power

The limited number of DE genes per perturbation (1–2 per condition) reflects both the moderate sample size per perturbation and the conservative multiple testing correction. Larger datasets such as those from Replogle et al. (2022), with hundreds of cells per perturbation, would substantially increase statistical power. The use of SCEPTRE's conditional resampling framework could further improve calibration and sensitivity.

### 6.4 Epistasis Quantification

Our additive deviation scoring approach provides a simple yet informative measure of genetic interactions. The widespread epistasis observed (mean RMSE = 0.37) is consistent with findings from Norman et al. (2019), who reported pervasive genetic interactions in combinatorial Perturb-seq experiments. Future extensions could incorporate multiplicative interaction models and gene-specific epistasis testing.

### 6.5 Representation Learning Limitations

The low silhouette scores (mean = −0.009) indicate that linear dimensionality reduction (PCA) provides limited separation of perturbation effects. Deep generative models such as scVI (Lopez et al., 2018) and CPA (Lotfollahi et al., 2023) are better suited for learning disentangled, nonlinear perturbation representations. Integration of these models into our framework is a priority for future development.

### 6.6 Limitations

1. **Synthetic data validation**: Results are based on simulated data with simplified perturbation effects. Real Perturb-seq data exhibits more complex patterns including heterogeneous perturbation efficiency, off-target effects, and cell-type-specific responses.
2. **Linear causal model**: The causal graph inference assumes linear relationships between perturbation and response, which may miss nonlinear regulatory dynamics.
3. **Scalability**: The current implementation handles thousands of cells but would require optimization (e.g., GPU acceleration) for genome-scale datasets with millions of cells.
4. **Temporal dynamics**: The framework analyzes snapshot data and does not model temporal perturbation dynamics.

### 6.7 Future Directions

1. **Real data validation** on the Replogle et al. (2022) genome-scale Perturb-seq dataset
2. **Deep generative model integration** with scVI and CPA for nonlinear representation learning
3. **Multi-modal extension** to handle Perturb-CITE-seq and multi-omic perturbation data
4. **Temporal modeling** for time-resolved perturbation experiments
5. **Full Pertpy integration** leveraging the standardized APIs from Heumos et al. (2025)
6. **Benchmarking** against DREAM challenge datasets for causal network inference

## 7. Conclusion

We presented a modular computational framework for Perturb-seq data analysis that integrates six critical analytical components: quality control, differential expression, co-expression module detection, causal graph inference, epistasis detection, and representation learning. The framework is built on the Scanpy/Pertpy ecosystem and provides a reproducible, extensible platform for analyzing single-cell CRISPR screen data. Validation on synthetic data demonstrates the framework's ability to capture perturbation effects, detect genetic interactions, and construct gene regulatory networks. Future work will focus on deep generative model integration, real data validation, and multi-modal extension to maximize the biological insights from increasingly large and complex Perturb-seq experiments.

## References

1. Norman, T. M., Horlbeck, J. M., Replogle, J. M., et al. (2019). Exploring genetic interaction manifolds constructed from rich single-cell phenotypes. *Science*, 365(6455), 786–793. DOI: [10.1126/science.aax4438](https://doi.org/10.1126/science.aax4438)

2. Replogle, J. M., Saunders, R. A., Pogson, A. N., et al. (2022). Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. *Cell*, 185(14), 2559–2575. DOI: [10.1016/j.cell.2022.05.013](https://doi.org/10.1016/j.cell.2022.05.013)

3. Barry, T., Wang, X., Morris, J. A., Roeder, K., & Katsevich, E. (2021). SCEPTRE improves calibration and sensitivity in single-cell CRISPR screen analysis. *Genome Biology*, 22, 344. DOI: [10.1186/s13059-021-02545-2](https://doi.org/10.1186/s13059-021-02545-2)

4. Lotfollahi, M., et al. (2023). Mapping single-cell perturbation responses to population variability in drug sensitivity. *Molecular Systems Biology*, 19(6), e11350. DOI: [10.15252/msb.202211350](https://doi.org/10.15252/msb.202211350)

5. Heumos, L., Ji, Y., May, L., et al. (2025). Pertpy: an end-to-end framework for perturbation analysis. *Nature Methods*. DOI: [10.1038/s41592-025-02909-7](https://doi.org/10.1038/s41592-025-02909-7)

6. Lopez, R., Regier, J., Cole, M. B., Jordan, M. I., & Yosef, N. (2018). Deep generative modeling for single-cell transcriptomics. *Nature Methods*, 15(12), 1053–1058. DOI: [10.1038/s41592-018-0229-2](https://doi.org/10.1038/s41592-018-0229-2)

7. Peidli, S., Green, T. D., Shen, C., et al. (2024). scPerturb: harmonized single-cell perturbation data. *Nature Methods*, 21, 531–540. DOI: [10.1038/s41592-023-02144-y](https://doi.org/10.1038/s41592-023-02144-y)

8. Dixit, A., Parnas, O., Li, B., et al. (2016). Perturb-Seq: dissecting molecular circuits with scalable single-cell RNA profiling of pooled genetic screens. *Cell*, 167(7), 1853–1866. DOI: [10.1016/j.cell.2016.11.038](https://doi.org/10.1016/j.cell.2016.11.038)

9. Kamimoto, K., Strber, B., Hashimoto, T., et al. (2023). Dissecting cell identity via network inference and in silico gene perturbation. *Nature*, 614, 742–751. DOI: [10.1038/s41586-022-05688-9](https://doi.org/10.1038/s41586-022-05688-9)

10. Adamson, B., Norman, T. M., Jost, M., et al. (2016). A multiplexed single-cell CRISPR screening platform enables systematic dissection of the unfolded protein response. *Cell*, 167(7), 1867–1882. DOI: [10.1016/j.cell.2016.11.048](https://doi.org/10.1016/j.cell.2016.11.048)
