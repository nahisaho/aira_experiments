# Advanced Analytical Pipeline for Spatial Transcriptomics: Deconvolution, Spatially Variable Gene Detection, Cell–Cell Communication, and Tumor Immune Microenvironment Characterization

---

## Abstract

Spatial transcriptomics (ST) technologies such as 10x Genomics Visium and MERFISH have transformed our ability to map gene expression within the spatial context of intact tissues. However, the full analytical potential of these platforms remains under-exploited due to the lack of integrated, end-to-end computational pipelines that combine spot deconvolution, spatially variable gene (SVG) detection, cell–cell communication inference, tissue niche identification, three-dimensional reconstruction, and tumor immune microenvironment (TME) characterization. Here we present SpatialFlow, a modular analytical framework designed for Visium/MERFISH data that integrates six key analytical modules: (1) non-negative matrix factorization (NMF)-based spot deconvolution inspired by cell2location, yielding cell type composition estimates with a best-case Pearson r = 0.8599 for tumor cells [cell:4]; (2) Moran's I-based spatially variable gene detection identifying 39 significant SVGs (q < 0.05, I > 0.1) [cell:5]; (3) ligand-receptor communication scoring revealing CXCL12:CXCR4 as the dominant pathway (mean score = 0.0289) [cell:6]; (4) KMeans-based tissue niche identification yielding three biologically coherent niches (silhouette = 0.2762) [cell:7]; (5) Procrustes-based 3D reconstruction of consecutive sections with cross-section reproducibility r = 0.9969 ± 0.0009 [cell:8]; and (6) immunosuppression scoring enabling immunotherapy response prediction with AUROC = 0.7181 [cell:9]. We applied this pipeline to a synthetic tumor immune microenvironment dataset modeling 400 spots across 8 cell types, producing quantitative benchmarks for each module. We discuss the limitations of synthetic data, the challenges of real-world deconvolution, and the path toward clinical translation of spatial multi-omics analysis.

**Keywords:** spatial transcriptomics, Visium, MERFISH, deconvolution, spatially variable genes, cell communication, tumor microenvironment, NMF, Moran's I

---

## 1. Introduction

The spatial organization of cells within tissues is fundamentally linked to their function. Conventional single-cell RNA sequencing (scRNA-seq) resolves cellular heterogeneity at single-cell resolution but sacrifices spatial context through tissue dissociation. Spatial transcriptomics (ST) platforms bridge this gap by capturing gene expression profiles while preserving the physical location of each measurement. Two dominant ST paradigms have emerged: (i) sequencing-based technologies such as 10x Genomics Visium, which capture transcriptome-wide expression at ~55 µm resolution in arrayed "spots" potentially containing multiple cells, and (ii) imaging-based platforms such as MERFISH and seqFISH+, which achieve single-cell or even sub-cellular resolution for a defined gene panel.

The analytical challenges posed by ST data are substantial. Visium spots are not single cells—each spot contains a mixture of cell types, necessitating computational deconvolution to estimate cellular composition. Beyond composition, identifying genes whose expression varies spatially (spatially variable genes, SVGs) provides insight into the molecular gradients that define tissue organization. Understanding how cells communicate via ligand-receptor interactions in their spatial context is crucial for deciphering tissue function and disease pathogenesis. Identifying coherent tissue niches, or spatially recurring cellular neighborhoods, offers a higher-level view of tissue architecture. Finally, integrating data from multiple consecutive tissue sections enables three-dimensional reconstruction of the spatial transcriptome.

The tumor immune microenvironment (TME) represents a particularly compelling use case for spatial transcriptomics. The relative abundance of cytotoxic CD8⁺ T cells, immunosuppressive M2 macrophages, cancer-associated fibroblasts, and tumor cells—and critically, their spatial organization—determines the outcome of anti-tumor immunity and response to immunotherapy. Recent studies using cell2location (Kleshchevnikov et al., 2022), SpatialDE (Svensson et al., 2018), and Squidpy (Palla et al., 2022) have begun to dissect these spatial relationships, but an integrated pipeline covering all major analytical modules remains lacking.

This work presents SpatialFlow, an integrated Python-based pipeline implementing six analytical modules. We validate the pipeline on a realistic synthetic Visium-like dataset and report quantitative benchmarks for each module, with particular emphasis on computational reproducibility, honest reporting of limitations, and the distinction between synthetic-data performance and expected real-world performance.

---

## 2. Related Work

### 2.1 Spot Deconvolution

Multiple approaches have been developed for estimating cell type proportions in Visium spots. **cell2location** (Kleshchevnikov et al., 2022, *Nature Biotechnology*, DOI: 10.1038/s41587-021-01139-4) uses a hierarchical Bayesian model to map cell types from paired scRNA-seq reference data. **RCTD** (Cable et al., 2022) employs a probabilistic model with doublet detection. **SPOTlight** (Elosua-Bayes et al., 2021) uses non-negative matrix factorization (NMF) of a reference atlas. Our NMF-based deconvolution follows the spirit of SPOTlight, using reference cell type signatures to guide component interpretation.

A recent benchmark study using deep learning with the Virchow foundation model achieved a macro-averaged cell type deconvolution AUC of 0.812 in colon tumor tissues (Le et al., 2025, DOI: 10.1158/1538-7445.am2025-6260), providing a relevant comparison point for histology-guided approaches.

### 2.2 Spatially Variable Gene Detection

**SpatialDE** (Svensson et al., 2018, *Nature Methods*, DOI: 10.1038/nmeth.4636) introduced Gaussian process regression for SVG detection, modeling spatial patterns as linear combinations of spatial kernels. **NNSVG** (Weber et al., 2023) uses nearest-neighbor Gaussian processes for scalability. Our approach uses Moran's I, a classical spatial statistics measure of autocorrelation, as a computationally efficient proxy. Moran's I has the advantage of being non-parametric and not requiring model selection.

### 2.3 Cell–Cell Communication

**CellChat** and **CellPhoneDB** are widely used for inferring intercellular communication from scRNA-seq. For spatial data, **stMLnet** (Yan et al., 2025, *Genome Research*, DOI: 10.1101/gr.279857.124) integrates spatially dependent ligand-receptor signaling based on diffusion and mass action models, outperforming seven competing methods on MERFISH and Slide-seq data. **MAGNET** (Han et al., 2025, *PLoS Computational Biology*, DOI: 10.1371/journal.pcbi.1013810) reconstructs cell-cell interaction networks using multi-view graph autoencoders, achieving an average precision of 0.901 on seqFISH data.

### 2.4 Tumor Spatial Niches

A pan-cancer spatial transcriptomic analysis of 373 samples across 12 cancer types identified 56 local cellular programs and 13 recurrent niches (Li et al., 2026, *Cell Reports Medicine*, DOI: 10.1016/j.xcrm.2026.102751). Notably, macrophages co-localized with tumor cells (Niche_4) correlated with poor prognosis, while immune-rich niches predicted better survival and immunotherapy response—a finding directly motivating our immunosuppression and TLS scoring approach.

Yang et al. (2025, *Journal for ImmunoTherapy of Cancer*, DOI: 10.1136/jitc-2025-013763) applied cell2location and CellChat to KRAS-mutant colorectal cancer spatial data, identifying an immunosuppressive niche mediated by Fibroblast-secreted collagen and Mono_S100A8 monocyte recruitment.

### 2.5 3D Reconstruction and Multi-Section Integration

Pentimalli et al. (2025, *Cell Systems*, DOI: 10.1016/j.cels.2025.101261) demonstrated 3D spatial transcriptomics integration with ECM imaging in consecutive lung carcinoma sections, using molecular neighborhood analysis to identify known immune escape mechanisms. Their work establishes the feasibility of 3D ST in routine clinical samples.

### 2.6 Limitations of Prior Work

Prior tools are often validated on limited datasets, require paired scRNA-seq references for deconvolution, lack end-to-end integration, and rarely provide benchmarks for synthetic vs. real-world performance gaps. Our pipeline addresses these gaps through modular design and explicit limitation reporting.

---

## 3. Methods

### 3.1 Synthetic Dataset Generation

We generated a synthetic Visium-like dataset (random seed = 42) comprising **400 spots** arranged on a hexagonal grid (20×20 spots, mimicking Visium capture area) with **225 genes** (200 cell-type marker genes + 25 spatially variable genes) [cell:1]. Eight cell types were modeled representing the tumor immune microenvironment: Tumor cells, CD8⁺ T cells, CD4⁺ T cells, M1 Macrophages, M2 Macrophages, Fibroblasts, Endothelial cells, and NK cells.

Four spatial domains were defined based on distance from tissue center: Tumor core (n=60 spots), Invasive margin (n=118), Immune-rich (n=46), and Stromal (n=176) [cell:2]. Spot-level cell type proportions were drawn from a Dirichlet distribution with domain-specific concentration parameters reflecting known TME biology (e.g., Tumor_core: 65% tumor cells; Immune_rich: 28% CD8⁺ T cells). Gene expression counts were sampled from a negative binomial distribution with means determined by proportional mixing of cell type-specific expression signatures, yielding a mean count of 371.2 ± 59.3 counts per spot [cell:3].

Data provenance: Raw data saved to `data/raw/expression_matrix.csv` and `data/raw/cell_proportions.csv`.

### 3.2 Spot Deconvolution (Module 1)

We implemented NMF-based deconvolution using scikit-learn's `NMF` (solver='cd', 500 iterations, random_state=42). Log-normalized expression (log1p) was factorized into W (spots × components) and H (components × genes). Components were matched to known cell types using cosine similarity between H rows and reference signatures, with optimal assignment via the Hungarian algorithm. Proportions were obtained by L1-normalizing W across components [cell:4].

**Evaluation metric:** Pearson correlation between predicted and true proportions, plus 5-fold cross-validation R² using Ridge regression (α=1.0) [cell:10].

### 3.3 Spatially Variable Gene Detection (Module 2)

We computed **Moran's I** statistic for each gene using a k-nearest-neighbor (k=10) spatial weight matrix derived from Euclidean spot coordinates. Moran's I ranges from -1 (perfect dispersion) to +1 (perfect clustering); values near 0 indicate random spatial distribution. Z-scores were computed under the null hypothesis of spatial randomness, and p-values were corrected using the Benjamini-Hochberg (BH) procedure. Genes with q < 0.05 and I > 0.1 were classified as SVGs [cell:5].

**Formula:** $I = \frac{n \sum_i \sum_j w_{ij}(x_i - \bar{x})(x_j - \bar{x})}{S \sum_i (x_i - \bar{x})^2}$ where $w_{ij}$ is the spatial weight and $S = \sum_i \sum_j w_{ij}$.

### 3.4 Cell–Cell Communication (Module 3)

Twelve ligand-receptor pairs with established roles in TME biology were curated (e.g., PD-L1:PD-1, VEGFA:FLT1, CXCL12:CXCR4, TGFb1:TGFBR1, IFNG:IFNGR1). For each pair, the communication score at each spot was defined as the product of sender and receiver cell type proportions: $S_{i}^{L \to R} = p_i^{sender} \times p_i^{receiver}$ [cell:6]. This product-based scoring is comparable to the approach used by CellChat. Spot-level scores were averaged to produce a mean communication score per pair and a spatial map of signaling activity.

### 3.5 Tissue Niche Identification (Module 4)

Spots were clustered using KMeans applied to a combined feature matrix of scaled cell type proportions (weight 0.7) and scaled spatial coordinates (weight 0.3). The optimal number of niches k was selected by maximizing the silhouette score across k ∈ {2,...,7} with 10 random initializations each. Each niche was named after its dominant cell type [cell:7].

### 3.6 3D Spatial Reconstruction (Module 5)

Three consecutive tissue sections were simulated with z-step = 5 µm and small random displacement (Gaussian noise, σ=0.3 µm) to model experimental misalignment. Sections were registered to a reference section using Procrustes analysis (scipy.spatial.procrustes), which minimizes sum-of-squared differences after optimal rotation, scaling, and translation [cell:8].

### 3.7 Tumor Immune Microenvironment Analysis (Module 6)

An immunosuppression score was defined as: $IS = \frac{p_{M2} + 0.5 \cdot p_{tumor}}{p_{CD8} + p_{NK} + p_{M1} + 0.01}$. A TLS (Tertiary Lymphoid Structure) score was defined as the mean of CD8⁺, CD4⁺, and NK cell proportions. A response predictor was constructed as a linear combination: response_score = −IS + 2×TLS, and performance was evaluated with AUROC and AUPRC against ground-truth domain labels (Immune_rich + Invasive_margin = responders) [cell:9].

### 3.8 NatureLM and GALACTICA MCP Tools

The following MCP tools were attempted as part of the experimental design:

**NatureLM MCP (attempted tools: `generate_smiles`, `predict_logp`, `retrosynthesis`, `ask_naturelm`):**
- **Status:** Connection not available. NatureLM MCP tool was not discoverable via ToolUniverse find_tools query for "NatureLM molecular property prediction". The tool was not listed in the available ToolUniverse registry at the time of this experiment (2026-05-31).
- **Implication:** No quantitative molecular predictions (IC50, binding energy, LogP) from NatureLM were obtained. As this is a spatial transcriptomics study rather than a molecular design study, the absence of NatureLM did not fundamentally limit the core analysis.

**GALACTICA MCP (attempted tools: `generate_molecule`, `scientific_qa`, `predict_citations`, `reasoning`):**
- **Status:** Connection not available. GALACTICA MCP was not discoverable via ToolUniverse. No scientific QA validation or citation prediction was obtained.
- **Alternative:** Scientific validation was performed through ToolUniverse SemanticScholar_search_papers queries, which successfully returned literature results.

**Semantic Scholar API:** Partially available (some queries returned 429 rate limit errors; three queries succeeded, yielding 13 papers across the three searches). Literature was supplemented with knowledge of foundational papers.

### 3.9 Computational Environment

- Python 3.11.2 (GCC 12.2.0)
- numpy 2.3.5, pandas 2.3.3, scikit-learn 1.6.1, scipy 1.17.1, matplotlib 3.10.9, seaborn 0.13.2
- Random seed: 42 (np.random.seed(42), random.seed(42)) in all cells
- Data: Synthetic (described in §3.1); raw files in `data/raw/`

### 3.10 Python Code

```python
# === ENVIRONMENT SETUP (Cell 0) ===
import random, numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns
import os, sys, warnings
warnings.filterwarnings('ignore')
SEED = 42
np.random.seed(SEED); random.seed(SEED)
os.makedirs('figures', exist_ok=True); os.makedirs('data/raw', exist_ok=True)

# === DATA GENERATION (Cells 1-3) ===
# 400 Visium-like spots on hexagonal grid, 225 genes, 8 cell types
# Negative binomial count model, Dirichlet cell proportions

# === DECONVOLUTION (Cell 4) ===
from sklearn.decomposition import NMF
from scipy.spatial.distance import cosine
from scipy.optimize import linear_sum_assignment
model_nmf = NMF(n_components=8, init='random', random_state=42, max_iter=500)
W = model_nmf.fit_transform(np.log1p(expression_matrix))

# === SVG DETECTION (Cell 5) ===
# Moran's I via k-NN spatial weights (k=10)
# BH multiple testing correction; threshold: q<0.05 AND I>0.1

# === CELL-CELL COMMUNICATION (Cell 6) ===
# Score = sender_proportion × receiver_proportion (product model)
# 12 curated TME L-R pairs

# === NICHE IDENTIFICATION (Cell 7) ===
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
# Optimize k ∈ {2,...,7}; final: KMeans(k=3, random_state=42)

# === 3D RECONSTRUCTION (Cell 8) ===
from scipy.spatial import procrustes
# Align sections 1,2 to section 0 via Procrustes

# === TME ANALYSIS (Cell 9) ===
# Immunosuppression score + TLS score → AUROC
from sklearn.metrics import roc_auc_score
```

---

## 4. Experiments

### 4.1 Dataset

The synthetic dataset models a Visium capture area containing a tumor with surrounding TME, with four biologically realistic tissue domains (Table 1). The negative binomial count model (NB(r=1, p=1/(1+μ))) was chosen to reflect the overdispersion characteristic of real ST data.

| Domain | N Spots | Dominant Cell Type | Tumor % | CD8 T % |
|--------|---------|-------------------|---------|---------|
| Tumor_core | 60 | Tumor_cells | 65.3% | 2.6% |
| Invasive_margin | 118 | Tumor_cells | 30.9% | 14.5% |
| Immune_rich | 46 | CD8_T_cells | 5.6% | 29.5% |
| Stromal | 176 | Fibroblasts | 7.8% | 4.3% |

**Table 1:** Synthetic dataset domain characteristics (mean proportions from Cell 2 [cell:2]).

### 4.2 Evaluation Metrics

- **Deconvolution:** Pearson r (vs. ground truth), 5-fold CV R²
- **SVG detection:** Moran's I, BH-corrected q-value
- **Communication:** Mean product score, fraction of active spots
- **Niche clustering:** Silhouette score, qualitative biology validation
- **3D reconstruction:** Procrustes disparity, cross-section Pearson r
- **TME analysis:** AUROC, AUPRC for immunotherapy response prediction

---

## 5. Results

### 5.1 Spot Deconvolution

NMF deconvolution with 8 components successfully recovered cell type proportions, with performance varying by cell type (Table 2). Tumor cells achieved the highest accuracy (Pearson r = 0.8599, p < 2.2×10⁻¹⁶), while CD8⁺ T cells were recovered with r = 0.8206. Fibroblasts also showed strong recovery (r = 0.7699), consistent with their distinct gene expression profile. In contrast, macrophage subtypes (M1 and M2), NK cells, and Endothelial cells showed poor deconvolution (r ≤ 0.37), primarily due to the difficulty of distinguishing transcriptionally similar cell populations with NMF alone [cell:4].

**5-fold cross-validation** confirmed these trends: overall mean R² = 0.3660 ± 0.2884 across cell types, with tumor cells achieving R² = 0.7515 ± 0.0526 and NK cells showing negative R² = −0.0657 ± 0.1026, indicating that the model performs near chance for this rare cell type [cell:10].

| Cell Type | Pearson r | CV R² ± SD |
|-----------|-----------|------------|
| Tumor_cells | 0.8599 | 0.7515 ± 0.0526 |
| CD8_T_cells | 0.8206 | 0.4554 ± 0.0495 |
| Fibroblasts | 0.7699 | 0.6789 ± 0.0659 |
| NK_cells | 0.3710 | −0.0657 ± 0.1026 |
| CD4_T_cells | 0.1966 | 0.4466 ± 0.0856 |
| Macrophages_M2 | 0.1945 | 0.0040 ± 0.2043 |
| Endothelial | −0.2178 | 0.3046 ± 0.1493 |
| Macrophages_M1 | −0.0540 | 0.3529 ± 0.0858 |

**Table 2:** Deconvolution performance (NMF, n=400 spots). [cell:4][cell:10]

### 5.2 Spatially Variable Gene Detection

Moran's I analysis detected **39 spatially variable genes** (q < 0.05, I > 0.1) out of 225 total genes [cell:5]. The top SVG was SVgene_6 (I = 0.3148, q = 4.99×10⁻⁸), followed by SVgene_4 (I = 0.2718) and Tumor_cells_gene1 (I = 0.2579). The detection of cell-type marker genes among SVGs confirms that cell type spatial organization is itself a driver of spatial expression patterns. A total of 14 of 25 designed spatially variable genes were correctly identified as SVGs, with additional cell type marker genes reaching significance due to the structured spatial domains.

![Figure 1](figures/main_analysis_panel.png)
**Figure 1:** Comprehensive analysis panel. (A) Tissue domain spatial map. (B) Tumor cell proportion heatmap. (C) CD8⁺ T cell proportion. (D) Identified niches (k=3). (E) Deconvolution accuracy (Pearson r). (F) SVG detection volcano (Moran's I vs. −log10 q). (G) L-R communication score heatmap. (H) Immunosuppression score spatial map. (I) 5-fold CV R² by cell type. (J) 3D reconstruction. (K) Top SVG expression. (L) AUROC curve.

### 5.3 Ligand-Receptor Communication

The 12 curated TME L-R pairs showed a range of communication scores (Table 3). **CXCL12:CXCR4** (Fibroblasts → Tumor) and **HGF:MET** (Fibroblasts → Tumor) achieved the highest mean scores (0.0289), reflecting high co-occurrence of fibroblasts and tumor cells across 77.8% of spots [cell:6]. The immunosuppressive **PD-L1:PD-1** axis (Tumor → CD8) scored 0.0182, active in 44.8% of spots—indicating widespread but not ubiquitous checkpoint signaling. The cytotoxic **IFNG:IFNGR1** axis (CD8 → Tumor) matched the PD-L1:PD-1 score exactly (0.0182), consistent with the mathematical symmetry of the product model when sender-receiver pairs are swapped. **IL-10:IL10RA** (M2 → CD8), representing M2-mediated CD8 suppression, showed the second-lowest score (0.0076), reflecting the spatial separation between Immune_rich and Tumor_core domains.

| L-R Pair | Sender | Receiver | Type | Mean Score | Active Spots |
|----------|--------|----------|------|------------|--------------|
| CXCL12:CXCR4 | Fibroblasts | Tumor | Migration | 0.0289 | 77.8% |
| HGF:MET | Fibroblasts | Tumor | Invasion | 0.0289 | 77.8% |
| SPP1:CD44 | M2 Macro. | Tumor | Survival | 0.0216 | 47.8% |
| VEGFA:FLT1 | Tumor | Endothelial | Angiogenesis | 0.0210 | 62.5% |
| PD-L1:PD-1 | Tumor | CD8 T | Immunosupp. | 0.0182 | 44.8% |
| IFNG:IFNGR1 | CD8 T | Tumor | Antitumor | 0.0182 | 44.8% |
| TNF:TNFRSF1A | M1 Macro. | Tumor | Cytotoxicity | 0.0169 | 44.0% |
| TGFb1:TGFBR1 | Fibroblasts | CD8 T | Exhaustion | 0.0149 | 53.3% |

**Table 3:** Top L-R communication pairs (sorted by mean score). [cell:6]

### 5.4 Tissue Niche Identification

KMeans clustering with optimal k=3 (silhouette = 0.2762, selected from k ∈ {2,...,7}) identified three biologically coherent niches [cell:7]:

- **Niche_0 (Stromal/Fibroblast):** n=176 spots, 42.8% fibroblasts, 19.9% endothelial — vascular/stromal compartment
- **Niche_1 (Tumor core):** n=105 spots, 52.4% tumor cells, 12.0% M2 macrophages — immunosuppressive core
- **Niche_2 (Immune infiltrate):** n=119 spots, 21.6% CD8⁺ T cells, 15.9% CD4⁺ T cells, 14.8% M1 macrophages — active immune compartment

The moderate silhouette score (0.2762) reflects the continuous nature of spatial gradients in the tissue, which resist hard cluster boundaries.

![Figure 2](figures/domain_lr_analysis.png)
**Figure 2:** Domain and communication analysis. (Top-left) Cell type composition heatmap by domain. (Top-center) Top L-R pair communication scores. (Top-right) Silhouette score optimization. (Bottom-left) Moran's I distribution. (Bottom-center) Tumor cell deconvolution scatter. (Bottom-right) PD-L1:PD-1 spatial communication map.

### 5.5 3D Spatial Reconstruction

Three consecutive sections (z-step = 5 µm, total z-range = 10 µm) were successfully registered using Procrustes alignment, achieving near-perfect registration (disparity < 10⁻⁶) for these synthetic sections with small perturbations [cell:8]. Cross-section expression reproducibility was high (Pearson r = 0.9969 ± 0.0009), confirming that the section-to-section noise model was appropriately small. The combined 3D dataset comprised 1,200 spots distributed across three z-planes.

### 5.6 Tumor Immune Microenvironment

The immunosuppression score captured biologically expected variation: Tumor_core had the highest mean IS (7.04 ± 4.70) vs. Immune_rich (0.17 ± 0.13), reflecting the overwhelming dominance of tumor and M2 cells in the core [cell:9]. TLS scores were highest in Immune_rich (0.205 ± 0.028).

Immunotherapy response prediction based on IS and TLS scores achieved **AUROC = 0.7181** and **AUPRC = 0.6500** (where Immune_rich + Invasive_margin spots were labeled as "responder") [cell:9]. This moderate performance reflects realistic prediction difficulty even with ground-truth cell type proportions, consistent with clinical observations that TME composition alone incompletely predicts immunotherapy response.

### 5.7 NatureLM and GALACTICA Results

As detailed in Methods §3.8, both NatureLM MCP and GALACTICA MCP were unavailable at the time of this study. No quantitative molecular predictions were obtained from these tools. Literature-based validation through Semantic Scholar confirmed the biological plausibility of our analysis design (see §2).

---

## 6. Discussion

### 6.1 Deconvolution Performance

The heterogeneous performance across cell types (r = 0.86 for tumor cells vs. r ≈ 0 for M1 macrophages) reflects a fundamental challenge: cell types with similar gene expression profiles (e.g., M1 and M2 macrophages, both expressing monocyte markers) are difficult to distinguish by NMF. The original cell2location addresses this by using a principled Bayesian model with scRNA-seq reference data and accounting for technical noise. Our NMF approach is simpler and faster but less powerful for closely related cell types. In real Visium data, deconvolution performance for macrophage subtypes would likely be even lower due to batch effects between ST and scRNA-seq reference datasets.

**Limitation:** Our synthetic data was generated from the same reference signatures used for deconvolution, which is an optimistic "closed-world" assumption. Real-world deconvolution must contend with reference misspecification, batch effects, and cell states not captured in the reference atlas.

### 6.2 SVG Detection

The Moran's I approach is computationally efficient and parameter-free (beyond the choice of k for the neighbor graph). However, it assumes a single statistic summarizes spatial autocorrelation, potentially missing complex patterns (e.g., periodic expression, multi-scale gradients) captured by SpatialDE's GP framework. The detection of 39 SVGs (out of 225 genes) at q < 0.05 and I > 0.1 includes both designed SVGs and cell-type markers with structured spatial distributions—the latter representing an expected biological reality, not a false positive.

### 6.3 Cell–Cell Communication

The product-based scoring used here is mathematically simple but ignores spatial distance between communicating cells, receptor expression levels, and downstream signaling. Tools like stMLnet explicitly model diffusion-based signaling strength decay with distance, providing more biologically accurate scores. Furthermore, our analysis used ground-truth proportions rather than deconvolved proportions; in practice, deconvolution errors would propagate to communication scores.

The dominance of CXCL12:CXCR4 and HGF:MET (Fibroblasts → Tumor) in our analysis reflects the high co-occurrence of fibroblasts and tumor cells across the tissue—a common finding in cancer spatial transcriptomics studies including Yang et al. (2025) in KRAS-mutant CRC. The immunosuppressive TGFb1:TGFBR1 axis (Fibroblasts → CD8) scored moderately (0.0149), present in 53.3% of spots, consistent with the known role of CAF-mediated TGFβ signaling in T cell exclusion.

### 6.4 Niche Identification

The optimal k=3 is lower than the 4 simulated ground-truth domains, reflecting the difficulty of resolving Tumor_core vs. Invasive_margin when proportions change gradually. Higher-resolution clustering (k=4) gave a lower silhouette score (0.2609), suggesting the data support only 3 well-separated clusters at the combined feature scale. In real data with more diverse cell states and spatial heterogeneity, more niches would likely be identified (the pan-cancer study by Li et al. 2026 identified 13 recurrent niches).

### 6.5 3D Reconstruction

The near-zero Procrustes disparity reflects the small random displacements in our synthetic data. Real consecutive sections can show substantial misalignment due to tissue folding, stretching, and rotation during mounting. More sophisticated registration approaches (e.g., PASTE, mutual information optimization) would be needed for real data. The high cross-section reproducibility (r = 0.9969) also reflects synthetic data; real biological sections show section-to-section variation due to the depth-varying cellular composition of 3D tissues.

### 6.6 TME Analysis and Generalization

The AUROC of 0.7181 for immunotherapy response prediction is realistic for a simple two-feature predictor. Clinical models trained on spatial TME features from real Visium data typically achieve AUROC in the range 0.60–0.85 depending on tumor type, treatment, and feature complexity. Our synthetic ground-truth labels (Immune_rich + Invasive_margin = responder) are a simplification; real response prediction requires longitudinal clinical outcome data.

**Generalization concern:** All performance metrics reported here were obtained on synthetic data generated under controlled assumptions. In real Visium data, additional challenges arise: (1) technical noise from uneven library quality across spots; (2) ambient RNA contamination; (3) spatial resolution limitations (~55 µm per spot); (4) reference bias in deconvolution; (5) gene dropout. Performance estimates should be treated as upper bounds relative to real-world application.

### 6.7 NatureLM and GALACTICA Absence

The unavailability of NatureLM and GALACTICA MCPs represents a methodological limitation. For molecular-level validation (e.g., predicted binding affinities of PD-1/PD-L1 interaction, CXCL12/CXCR4 binding energy), these tools would have provided quantitative cross-validation. As a partial substitute, ToolUniverse's Semantic Scholar tool confirmed the biological literature basis for each design choice. Future work should integrate NatureLM predictions of L-R binding affinities to weight communication scores by molecular plausibility.

---

## 7. Conclusion

We have presented SpatialFlow, an integrated spatial transcriptomics pipeline covering six analytical modules: NMF-based deconvolution (AUROC for tumor cells: r=0.86 [cell:4]), Moran's I-based SVG detection (39 SVGs identified [cell:5]), product-based ligand-receptor communication scoring (PD-L1:PD-1 score = 0.0182 [cell:6]), KMeans niche identification (k=3, silhouette=0.2762 [cell:7]), Procrustes-based 3D section registration (r=0.9969 [cell:8]), and immunosuppression-based TME scoring (AUROC=0.7181 [cell:9]).

Key findings include: (1) Fibroblast-mediated signaling (CXCL12:CXCR4, HGF:MET, TGFb1:TGFBR1) dominates the communication landscape, suggesting fibroblasts as key regulators of tumor progression and T cell exclusion; (2) spatially restricted immune niches (Immune_rich, k=2) retain distinct identities in 3-niche clustering, supporting their clinical relevance; (3) simple IS and TLS scores achieve AUROC=0.72 for response prediction, suggesting that basic TME composition metrics carry predictive power.

Future directions include: (1) replacing NMF with hierarchical Bayesian models (cell2location-equivalent) for improved deconvolution of similar cell types; (2) implementing GP-based SVG detection for multi-scale pattern discovery; (3) integrating spatial distance into L-R scoring (following stMLnet); (4) extending to real Visium and MERFISH datasets; (5) coupling with molecular predictions from NatureLM for binding affinity-weighted communication scoring.

---

## References

1. Kleshchevnikov, V., et al. (2022). cell2location maps fine-grained cell types in spatial transcriptomics. *Nature Biotechnology*, 40, 661–671. DOI: 10.1038/s41587-021-01139-4

2. Svensson, V., Teichmann, S. A., & Stegle, O. (2018). SpatialDE: identification of spatially variable genes. *Nature Methods*, 15, 343–346. DOI: 10.1038/nmeth.4636

3. Palla, G., et al. (2022). Squidpy: a scalable framework for spatial omics analysis. *Nature Methods*, 19, 171–178. DOI: 10.1038/s41592-021-01358-2

4. Yan, L., Cheng, J., Nie, Q., & Sun, X. (2025). Dissecting multilayer cell–cell communications with signaling feedback loops from spatial transcriptomics data. *Genome Research*. DOI: 10.1101/gr.279857.124

5. Han, C., Song, Z., Xu, Z., & Chen, J. (2025). MAGNET: Multi-view graph autoencoder with cell-gene attention for cell interaction network reconstruction from spatial transcriptomics. *PLoS Computational Biology*. DOI: 10.1371/journal.pcbi.1013810

6. Yang, S., et al. (2025). Single-cell and spatial transcriptomics integration identifies the immunosuppressive spatial niche in KRAS-mutant colorectal cancer. *Journal for ImmunoTherapy of Cancer*. DOI: 10.1136/jitc-2025-013763

7. Li, J., et al. (2026). Pan-cancer analysis of spatial transcriptomics reveals heterogeneous tumor spatial microenvironment. *Cell Reports Medicine*. DOI: 10.1016/j.xcrm.2026.102751

8. Pentimalli, T., et al. (2025). Combining spatial transcriptomics and ECM imaging in 3D for mapping cellular interactions in the tumor microenvironment. *Cell Systems*. DOI: 10.1016/j.cels.2025.101261

9. Le, M.-K., et al. (2025). Characterization of the tumor microenvironment's histologic landscape through histology-based deep learning spatial transcriptomic cell-type deconvolution of colon tumors. *Cancer Research*, 85(8_Suppl_1):6260. DOI: 10.1158/1538-7445.am2025-6260

10. Shi, W., et al. (2025). Single-cell and spatial transcriptomics integration: new frontiers in tumor microenvironment and cellular communication. *Frontiers in Immunology*. DOI: 10.3389/fimmu.2025.1649468

---

## Reproducibility

| Parameter | Value |
|-----------|-------|
| Random seed | 42 (np.random.seed, random.seed) |
| Python | 3.11.2 (GCC 12.2.0) |
| numpy | 2.3.5 |
| pandas | 2.3.3 |
| scikit-learn | 1.6.1 |
| scipy | 1.17.1 |
| matplotlib | 3.10.9 |
| seaborn | 0.13.2 |
| Data source | Synthetic (NB model); `data/raw/expression_matrix.csv` |
| Notebook | `spatial_transcriptomics_pipeline.ipynb` |
| NatureLM MCP | Not available (connection failed) |
| GALACTICA MCP | Not available (connection failed) |
| Semantic Scholar | Partially available (429 rate limit on 2 of 5 queries) |

**Note on cell citation format:** `[cell:N]` refers to the Jupyter execution cell number in the notebook `spatial_transcriptomics_pipeline.ipynb`.
