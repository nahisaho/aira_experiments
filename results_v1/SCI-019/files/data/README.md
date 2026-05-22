# Synthetic data documentation for the autoimmune systems immunology framework

## Overview

This workspace contains a fully synthetic, seed-controlled benchmark dataset designed for end-to-end development of a systems immunology analysis pipeline focused on autoimmune disease, with emphasis on rheumatoid arthritis (RA). All synthetic datasets are generated with `set.seed(42)` inside the corresponding R scripts to maximize reproducibility and to provide deterministic reference outputs for method development, figure generation, and pipeline integration.

The design mimics a realistic translational immunology workflow spanning bulk multi-omics, immune deconvolution, cytokine dynamical systems modeling, single-cell checkpoint profiling, treatment response prediction, and in silico tolerance restoration modeling.

## Synthetic data generation methodology

### 1. Multi-omics dataset (`scripts/01_multiomics_integration.R`)

Three omics layers are simulated for 80 samples divided into four equally sized groups:

- `RA_active` (n = 20)
- `RA_remission` (n = 20)
- `SLE` (n = 20)
- `HC` (healthy controls, n = 20)

#### Transcriptome
- Dimensions: **500 genes × 80 samples**
- Generation strategy:
  - Ten latent immunobiology factors are constructed from group-specific templates plus Gaussian noise.
  - Gene-level weights project latent structure into transcript abundance space.
  - Count-like data are generated using a negative binomial distribution.
  - Counts are normalized to log-CPM using `edgeR` for downstream `limma`/`voom` analysis.

#### Proteome
- Dimensions: **200 proteins × 80 samples**
- Generation strategy:
  - Same latent factor framework as transcriptome.
  - Continuous abundance values are generated with lower technical noise than transcriptome.
  - Values mimic log-scale protein abundance profiles measured by LC-MS or multiplex proteomics.

#### Metabolome
- Dimensions: **150 metabolites × 80 samples**
- Generation strategy:
  - Latent factors are propagated into metabolite abundance space.
  - Positive abundance constraints are enforced.
  - Signal amplitudes are smaller than transcriptomic effects to reflect typical metabolomics compression.

### 2. Immune deconvolution dataset (`scripts/02_cibersortx_deconvolution.R`)

- Dimensions: **22 LM22-like immune cell types × 80 samples**
- Generation strategy:
  - Dirichlet-like sampling is used to generate cell fraction vectors summing to 1 per sample.
  - Group-specific concentration parameters are modified to reproduce biologically plausible shifts.
  - Simulated clinical metadata include `DAS28` and `CRP`.

### 3. Cytokine dynamical system (`scripts/03_cytokine_ode_model.R`)

- State variables:
  - `IL6`
  - `TNFa`
  - `IL17A`
  - `IL10`
  - `TGFb`
  - `IFNg`
- Time window: **0 to 72 hours**
- Scenarios:
  - `HC_baseline`
  - `RA_active`
  - `Anti_TNF_treatment`
  - `JAK_inhibitor`

The ODE model uses Hill-type activation and inhibition terms to emulate nonlinear inflammatory and regulatory feedback.

### 4. Single-cell checkpoint dataset (`scripts/04_singlecell_checkpoint.R`)

- Dimensions: **5000 cells × 2000 genes**
- Conditions:
  - `RA_active` = 2500 cells
  - `HC` = 2500 cells
- Cell type composition:
  - CD4 T: 30%
  - CD8 T: 25%
  - B cells: 15%
  - NK cells: 10%
  - Monocytes: 12%
  - DCs: 8%

Checkpoint genes intentionally embedded in the matrix:
- `PDCD1` (PD-1)
- `CD274` (PD-L1)
- `CTLA4`
- `LAG3`
- `HAVCR2` (TIM-3)
- `TIGIT`
- `PDCD1LG2` (PD-L2)
- `IDO1`
- `FOXP3`

Cell-type-specific effects are added so that inhibitory receptors are stronger in T cells and ligand expression is enriched in antigen-presenting cells (APCs).

### 5. Treatment response prediction dataset (`scripts/05_treatment_response_prediction.R`)

- Cohort size: **200 RA patients**
- Response endpoint: **binary ACR50 at 6 months**
- Treatment arms:
  - `MTX`
  - `Anti-TNF`
  - `JAK inhibitor`
  - `IL-6R inhibitor`

Feature classes:
- Clinical: 8 variables
- Transcriptomic signatures: 50 variables
- Deconvolution-derived immune variables: 10 variables
- Baseline cytokines: 4 variables
- Treatment assignment: 1 categorical variable

The response probability is generated from a structured logistic signal combining disease activity, inflammatory tone, regulatory tone, transcriptomic biomarkers, and treatment-specific interaction effects.

### 6. Tolerance restoration simulations (`scripts/06_tolerance_insilico.R`)

- State variables:
  - `Teff`
  - `Treg`
  - `APC_active`
  - `APC_tolerogenic`
  - `IL10`
  - `TGFb`
  - `IL6`
  - `TNFa`
  - `IDO1`
  - `PDL1`
- Simulation window: **30 days**
- Strategies:
  1. Treg adoptive transfer
  2. Anti-inflammatory cytokines
  3. Checkpoint agonist
  4. IDO1 induction
  5. APC tolerogenic reprogramming
  6. Combination therapy

Monte Carlo perturbation analysis (1000 runs) is included to estimate robustness to parameter uncertainty.

## Data formats and column descriptions

### Common metadata columns
- `sample_id`: synthetic bulk sample identifier
- `patient_id`: synthetic patient identifier
- `cell_id`: synthetic single-cell identifier
- `group`: disease class for bulk multi-omics / deconvolution
- `condition`: RA-active vs HC label for single-cell data
- `cell_type`: manually simulated cell identity

### Key result table conventions
- `comparison`: group contrast (for example `RA_active_vs_HC`)
- `adj.P.Val` / `padj`: multiple-testing adjusted p-value (BH/FDR)
- `effect_size`: standardized or model-derived magnitude of change
- `AUC`: area under the ROC curve
- `tolerance_score`: composite regulatory minus inflammatory score
- `inflammatory_score`: summary inflammatory burden metric
- `regulatory_score`: summary suppressive cytokine metric

### Matrix orientation
- Bulk omics matrices: **rows = features, columns = samples**
- Immune deconvolution sample-level export: **rows = samples, columns = cell types**
- scRNA-seq count matrix stored in RDS object: **rows = genes, columns = cells**

## Biological basis for simulated group differences

The simulation explicitly encodes plausible immunopathology:

### RA active disease
- Increased CD8 T-cell abundance and inflammatory myeloid skewing
- Reduced Treg representation and lower compensatory tolerance tone
- Strong TNF/IL-6/IL-17 axis activity
- Enhanced checkpoint receptor expression in T cells consistent with chronic activation/exhaustion
- Worse baseline disease activity and lower treatment response when inflammatory burden is high

### RA remission
- Intermediate inflammatory burden
- Partial restoration of NK cells and Treg fractions
- Reduced but still detectable immune dysregulation

### SLE
- Strong B-cell/plasma-cell and activated dendritic-cell signatures
- Persistent Treg suppression
- Distinct but partially overlapping inflammatory architecture compared with RA

### Healthy control
- Balanced innate/adaptive composition
- Lower inflammatory cytokines
- Lower checkpoint induction outside homeostatic expression ranges

## How to replace synthetic data with real data

To adapt the framework for real experiments:

1. **Bulk multi-omics**
   - Replace synthetic matrices with normalized transcriptome/proteome/metabolome matrices.
   - Keep the same orientation: features in rows, samples in columns.
   - Update `sample_annot` with real phenotype and batch metadata.

2. **Immune deconvolution**
   - Replace simulated proportion matrix with actual CIBERSORTx or `immunedeconv` output.
   - Ensure sample IDs match downstream bulk omics metadata.

3. **Cytokine ODE model**
   - Fit production/degradation parameters using serum cytokine time-course data, perturbation studies, or literature priors.
   - Replace arbitrary units with calibrated concentrations if quantitative assays are available.

4. **scRNA-seq**
   - Replace synthetic counts with a raw count matrix (genes × cells).
   - Update cell metadata and, if needed, cell annotations after clustering.
   - Preserve checkpoint gene symbols in HGNC format.

5. **Treatment response model**
   - Replace synthetic patient table with baseline clinical/omics features and a binary responder label.
   - Refit all models and recalibrate probabilities.
   - Recompute subgroup AUCs and feature importance values.

6. **Tolerance restoration model**
   - Re-estimate parameters from mechanistic experiments, co-culture assays, or clinical biomarker dynamics.
   - Adjust efficacy metrics to reflect the therapeutic decision context.

## Reproducibility notes

- All scripts set `set.seed(42)`.
- All figures are written to `figures/`.
- All tabular results are written to `results/`.
- Synthetic processed data are exported to `data/` when the script is executed.
- Benchmark summaries in documentation assume the reference seed-controlled run described in `results/analysis_summary.md` and `report.md`.
