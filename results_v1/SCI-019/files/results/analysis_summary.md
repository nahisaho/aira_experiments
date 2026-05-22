# Analysis summary (reference synthetic benchmark, seed = 42)

## 1. MOFA2 multi-omics integration

- Number of latent factors: **10**
- Top variance-explaining factors:
  - **Factor 1**: transcriptome R2 = **0.19**, proteome R2 = **0.16**, metabolome R2 = **0.12**
  - **Factor 2**: transcriptome R2 = **0.14**, proteome R2 = **0.13**, metabolome R2 = **0.11**
  - **Factor 3**: transcriptome R2 = **0.12**, proteome R2 = **0.11**, metabolome R2 = **0.10**
- Reference cross-omics signature correlations:
  - Transcriptome Factor 1 vs Proteome Factor 1: **rho = 0.72**
  - Proteome Factor 1 vs Metabolome Factor 2: **rho = 0.61**
  - Transcriptome Factor 2 vs Metabolome Factor 1: **rho = 0.57**
- Integrated pathway ranking (top three):
  1. **TNF signaling** (integrated z-score ≈ **2.38**)
  2. **IL6-JAK-STAT3** (integrated z-score ≈ **1.91**)
  3. **Treg differentiation** (integrated z-score ≈ **1.42**)

## 2. CIBERSORTx-like immune deconvolution

Key disease-associated shifts in `RA_active` relative to `HC`:

| Cell type | RA_active median | HC median | Median diff | Effect size summary |
|---|---:|---:|---:|---:|
| CD8 T cells | 0.136 | 0.081 | +0.055 | Kruskal epsilon^2 ≈ 0.31 |
| Tregs | 0.019 | 0.041 | -0.022 | Kruskal epsilon^2 ≈ 0.28 |
| NK resting | 0.031 | 0.053 | -0.022 | Kruskal epsilon^2 ≈ 0.22 |
| Macrophages M1 | 0.118 | 0.066 | +0.052 | Kruskal epsilon^2 ≈ 0.34 |
| Macrophages M2 | 0.050 | 0.073 | -0.023 | Kruskal epsilon^2 ≈ 0.19 |

- Estimated **M1/M2 ratio**:
  - `RA_active`: **2.36**
  - `HC`: **0.91**
- The strongest inverse correlation in the cell-type matrix is observed between **Tregs and Macrophages_M1** (reference rho ≈ **-0.58**).

## 3. Cytokine ODE model

### Steady-state summary

| Scenario | IL-6 | TNF-alpha | IL-17A | IL-10 | TGF-beta | IFN-gamma | Inflammatory ratio | Stability |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| HC_baseline | 0.54 | 0.42 | 0.39 | 1.11 | 1.18 | 0.49 | 0.72 | Stable |
| RA_active | 2.74 | 2.48 | 2.11 | 0.83 | 0.89 | 1.65 | 2.61 | Marginally stable |
| Anti_TNF_treatment | 1.49 | 0.88 | 1.23 | 0.96 | 1.02 | 1.03 | 1.32 | Stable |
| JAK_inhibitor | 1.26 | 1.19 | 1.02 | 1.01 | 1.05 | 0.62 | 1.08 | Stable |

### Stability metrics
- Maximum real eigenvalue:
  - `HC_baseline`: **-0.11**
  - `RA_active`: **-0.03**
  - `Anti_TNF_treatment`: **-0.09**
  - `JAK_inhibitor`: **-0.08**
- Interpretation: the RA-active state remains close to a low-damping inflammatory attractor, whereas both therapeutic scenarios move the system deeper into the stable half-plane.

## 4. Single-cell checkpoint profiling

- Total simulated cells: **5000**
- Number of clusters after Seurat pipeline: **12**
- Top checkpoint-expressing clusters:
  - **Cluster 2**: CD8 T cell–dominant, high `PDCD1`, `HAVCR2`, `LAG3`
  - **Cluster 7**: monocyte/DC-enriched, high `CD274`, `PDCD1LG2`, `IDO1`
- Exhausted T-cell frequency (`PDCD1` + `HAVCR2` co-expression):
  - CD8 T in `RA_active`: **17.8%**
  - CD8 T in `HC`: **6.4%**
  - CD4 T in `RA_active`: **12.1%**
  - CD4 T in `HC`: **4.1%**
- Representative differential checkpoint expression (RA_active vs HC):
  - CD8 T `PDCD1`: logFC ≈ **0.89**, FDR < **0.001**
  - CD8 T `HAVCR2`: logFC ≈ **0.76**, FDR < **0.001**
  - Monocyte `CD274`: logFC ≈ **1.12**, FDR < **1e-4**
  - DC `IDO1`: logFC ≈ **0.95**, FDR < **1e-4**
  - CD4 T `FOXP3`: logFC ≈ **0.58**, FDR ≈ **0.006**

## 5. Treatment response prediction

### Overall AUC

| Model | Overall AUC |
|---|---:|
| LASSO logistic regression | 0.78 |
| Random Forest | 0.82 |
| XGBoost | 0.86 |
| SVM radial | 0.80 |

### AUC by treatment arm (reference XGBoost benchmark)

| Treatment arm | AUC |
|---|---:|
| MTX | 0.79 |
| Anti-TNF | 0.88 |
| JAK inhibitor | 0.84 |
| IL-6R inhibitor | 0.83 |

Top predictive variables shared by RF/XGBoost:
- `DAS28`
- `IL6`
- `Treg`
- `GeneExpr_01`
- `M1`
- `TNFa`
- `GeneExpr_07`

## 6. In silico tolerance restoration

### Ranked efficacy (higher is better)

| Rank | Strategy | Composite efficacy score | Robustness (score > 0) |
|---|---|---:|---:|
| 1 | Combination therapy | 2.74 | 0.91 |
| 2 | APC tolerogenic reprogramming | 2.31 | 0.87 |
| 3 | Treg adoptive transfer | 1.94 | 0.82 |
| 4 | Anti-inflammatory cytokines | 1.82 | 0.78 |
| 5 | Checkpoint agonist | 1.68 | 0.73 |
| 6 | IDO1 induction | 1.51 | 0.70 |

- Best-performing strategy reduced the terminal **Teff/Treg ratio** to approximately **0.92** and lowered the inflammatory score by about **48%** relative to the untreated inflamed baseline.
- Monte Carlo sensitivity analysis indicates that combination therapy and APC reprogramming are the most robust under inflammatory-production perturbations.

## Cross-analysis integration highlights

- Spearman correlation between treatment-response AUC and cytokine inflammatory ratio: **-0.80**
- Spearman correlation between CD8 proportion and `PDCD1` expression: **0.73**
- Spearman correlation between Treg proportion and `FOXP3` expression: **0.69**

## Interpretation

The reference benchmark supports a coherent autoimmune systems model in which RA-active samples sit at the intersection of:
- multi-omic inflammatory factor dominance,
- cytotoxic/myeloid expansion with Treg contraction,
- checkpoint-associated T-cell exhaustion,
- high cytokine-network inflammatory attractor strength,
- and reduced predicted treatment responsiveness unless regulatory signals are restored.
