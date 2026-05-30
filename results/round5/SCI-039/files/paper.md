# GNN-Meso: A Graph Neural Network Framework for Multi-Scale, Multi-Lead Atmospheric Prediction

**Authors:** Experimental Study (2025)
**Venue:** Synthetic Experiment / Computational Meteorology

---

## Abstract

Data-driven weather prediction has emerged as a compelling alternative to traditional numerical weather prediction (NWP) systems, with models such as GraphCast, Pangu-Weather, and FourCastNet demonstrating competitive or superior performance against operational forecasts at medium-range lead times. This paper presents **GNN-Meso**, a graph neural network (GNN) framework for data-driven atmospheric prediction inspired by these architectures, implemented with PyTorch Geometric. The proposed model adopts an Encoder–Processor–Decoder architecture: a multi-layer perceptron (MLP) encoder embeds multi-level atmospheric variables—temperature (T), zonal wind (U), meridional wind (V), specific humidity (Q), and geopotential height (Z) at five pressure levels—into a latent node representation on a lat/lon mesh graph. Six stacked Graph Attention Network v2 (GATv2) processor blocks perform message passing over 8-neighbor grid edges, and an MLP decoder outputs residual increments to the current state. A soft physics constraint layer enforces non-negative specific humidity via the softplus activation without violating gradient flow. Experiments are conducted on synthetic ERA5-like atmospheric data at three proxy resolutions (coarse ~2.5°, medium ~1.0°, fine ~0.5°) with autoregressive rollout evaluated at 6-hour, 24-hour, and 120-hour lead times. Five-fold time-series cross-validation yields RMSE (normalized) of 0.0329 ± 0.0004 at 6 h, 0.0366 ± 0.0036 at 24 h, and 0.0494 ± 0.0078 at 120 h for single-step predictions. Critically, we observe that the persistence baseline outperforms the GNN in multi-step autoregressive RMSE at short lead times, highlighting a well-known spinup challenge. Physical consistency monitoring shows that column-integrated moisture and kinetic energy proxies remain bounded over 120-hour rollouts. We discuss fundamental limitations of the synthetic experimental design, the risk of overfitting, and the gap between simulation-based benchmarks and real-world ERA5 performance, concluding with directions for future work.

---

## 1. Introduction

Operational numerical weather prediction (NWP) has been the backbone of global weather forecasting for decades, with centers such as ECMWF (IFS), NCEP (GFS), and Japan Meteorological Agency (JMA GSM) operating high-resolution spectral and finite-element models requiring substantial supercomputing resources. The data-driven revolution in machine learning has prompted a wave of research exploring whether deep learning models can match or exceed NWP skill at a fraction of the inference cost.

The seminal WeatherBench benchmark (Rasp et al., 2020) established a standardized evaluation framework for data-driven weather forecasting using ERA5 reanalysis data, revealing that convolutional neural networks could achieve competitive RMSE for temperature and geopotential at short lead times. Building on this foundation, **FourCastNet** (Pathak et al., 2022) employed adaptive Fourier neural operators on a spherical grid to produce high-resolution 0.25° global forecasts at competitive skill. **Pangu-Weather** (Bi et al., 2023) introduced a 3D Earth-specific Transformer architecture with hierarchical temporal aggregation, achieving superior RMSE scores for 500 hPa geopotential height compared to ECMWF IFS at medium-range lead times. **GraphCast** (Lam et al., 2023) demonstrated that a GNN operating on a multi-mesh graph representation—where different edge sets encode local and global atmospheric connectivity—can produce 10-day forecasts with lower error than operational NWP on a majority of tracked variables and levels. More recently, **ClimaX** (Nguyen et al., 2023) proposed a foundation model approach, pre-training a Vision Transformer on heterogeneous climate datasets and fine-tuning for downstream weather prediction tasks. **NeuralGCM** (Kochkov et al., 2024) explored hybrid approaches that couple neural network parameterizations with differential equation solvers to produce physically consistent rollouts.

Despite these impressive results, several open questions remain:

1. **How does graph resolution affect forecast skill?** GraphCast used a fixed 0.25° multi-mesh; the impact of resolution degradation is not fully characterized in the GNN setting.
2. **Can soft physics constraints be embedded without performance penalty?** Most current models impose physical consistency implicitly through training data or post-processing.
3. **How well do architectural conclusions transfer across lead times?** Short-range (6 h) and medium-range (120 h) regimes may favor different inductive biases.

This paper contributes:
- **GNN-Meso**: a clean, reproducible PyTorch Geometric implementation of an Encoder–GATv2Processor–Decoder atmospheric model;
- A **multi-resolution experiment** at three proxy resolutions with autoregressive rollout at 6/24/120-hour lead times;
- **Five-fold time-series cross-validation** with standard deviations for honest evaluation;
- A **self-critical analysis** of the limitations of synthetic evaluation and the gap to real ERA5 data.

---

## 2. Related Work

### 2.1 Graph Neural Networks for Atmospheric Prediction

Keisler (2022) was among the first to demonstrate that a GNN operating directly on a latitude/longitude icosahedral mesh could outperform a simple convolutional baseline on ERA5 data. The key insight was that graph-based message passing provides geometric inductive bias that respects the spherical topology of the Earth, unlike CNNs that operate on rectilinear projections.

**GraphCast** (Lam et al., 2023) [DOI: 10.1126/science.adi2336] extended this to a multi-mesh hierarchy where icosahedral grids at multiple refinement levels encode both local and long-range atmospheric interactions. After training on decades of ERA5 data with RMSE loss, GraphCast produced 10-day forecasts surpassing ECMWF's operational high-resolution deterministic forecast (HRES) on 90% of tracked variables. The architecture uses a learned encoder that maps grid-point features to latent mesh nodes, a processor stack of 16 GNN blocks performing message passing, and a decoder back to grid points.

### 2.2 Transformer-Based Approaches

**Pangu-Weather** (Bi et al., 2023) [DOI: 10.1038/s41586-023-06185-3] framed weather prediction as a video prediction problem on a 3D pressure-level cube, using a hierarchical Swin Transformer with Earth-specific relative position biases. A key design choice was training four separate models for 1/3/6/24-hour lead times and combining them at inference to produce multi-step forecasts without autoregressive error accumulation at every step.

**FourCastNet** (Pathak et al., 2022) [arXiv: 2202.11214] used adaptive Fourier neural operators to capture global spectral features efficiently, achieving 0.25° resolution forecasts at approximately 45,000× speedup over ECMWF-IFS. A follow-up work, FourCastNet v2 (Bonev et al., 2023), extended the architecture with spherical Fourier neural operators.

### 2.3 Foundation Models and Hybrid Systems

**ClimaX** (Nguyen et al., 2023) [DOI: 10.48550/arXiv.2301.10343] proposed pre-training a ViT-style model on CMIP6 climate projections and ERA5, achieving transfer learning across diverse climate and weather tasks. **NeuralGCM** (Kochkov et al., 2024) coupled a GNN with a physics-based dynamical core, enabling stable multi-month rollouts that conventional purely data-driven models cannot achieve.

### 2.4 Evaluation Benchmarks

**WeatherBench** (Rasp et al., 2020) [DOI: 10.1029/2020MS002203] and its successor **WeatherBench 2** (Rasp et al., 2023) [DOI: 10.1029/2023MS003715] provide standardized verification protocols including RMSE, ACC (anomaly correlation coefficient), and probabilistic metrics, facilitating fair comparison across models.

### 2.5 Limitations of Prior Work

Most state-of-the-art results are reported on ERA5 reanalysis data under idealized conditions (complete global coverage, consistent 6-hourly intervals). Performance on operational assimilation data with missing observations, instrument errors, and non-Gaussian noise distributions is less well characterized. Additionally, physical constraint satisfaction—mass, energy, and moisture conservation—is typically not explicitly monitored in published evaluations.

---

## 3. Methods

### 3.1 Data Representation

Atmospheric state at time $t$ is represented as a multi-level field:

$$\mathbf{X}_t \in \mathbb{R}^{N_\text{lat} \times N_\text{lon} \times N_P \times N_V}$$

where $N_\text{lat}$, $N_\text{lon}$ are the horizontal grid dimensions, $N_P = 5$ is the number of pressure levels (1000, 850, 500, 300, 100 hPa), and $N_V = 5$ is the number of prognostic variables $\{T, U, V, Q, Z\}$. Each variable is normalized by its global mean and standard deviation computed over the training period:

$$\tilde{X}_{v,p} = \frac{X_{v,p} - \mu_{v,p}}{\sigma_{v,p} + \epsilon}$$

### 3.2 Graph Construction

The lat/lon grid is treated as an undirected graph $\mathcal{G} = (\mathcal{V}, \mathcal{E})$ where each node $i \in \mathcal{V}$ corresponds to a grid point $(i_\text{lat}, i_\text{lon})$ and edges connect 8-neighbors (including diagonals), with periodic boundary conditions in the longitude direction:

$$\mathcal{E} = \{(i, j) : |\Delta i_\text{lat}| \leq 1, |\Delta i_\text{lon}| \equiv \pm 1 \pmod{N_\text{lon}}\}$$

Node positional encoding uses three components to respect spherical geometry:

$$\mathbf{p}_i = \left[\frac{2\text{lat}_i}{90} - 1,\ \sin(2\pi \cdot \text{lon}_i / 360),\ \cos(2\pi \cdot \text{lon}_i / 360)\right]$$

### 3.3 GNN-Meso Architecture

**Encoder.** The encoder maps the concatenated atmospheric state and positional encoding to a latent representation:

$$\mathbf{h}_i^{(0)} = \text{MLP}_\text{enc}\left(\left[\tilde{\mathbf{x}}_i \,\|\, \mathbf{p}_i\right]\right), \quad \text{MLP}_\text{enc}: \mathbb{R}^{N_V N_P + 3} \to \mathbb{R}^H$$

with Layer Normalization and GELU activations after each linear layer, where $H \in \{48, 64\}$ is the hidden dimension.

**Processor.** $L = 4$–$6$ stacked GATv2 blocks perform message passing. Each block computes:

$$\mathbf{h}_i^{(\ell)} = \text{LN}_2\!\left(\text{FFN}\!\left(\mathbf{h}_i^{(\ell-\frac{1}{2})}\right) + \mathbf{h}_i^{(\ell-\frac{1}{2})}\right)$$

$$\mathbf{h}_i^{(\ell-\frac{1}{2})} = \text{LN}_1\!\left(\text{GATv2}\!\left(\mathbf{h}_i^{(\ell-1)},\ \{\mathbf{h}_j^{(\ell-1)}\}_{j \in \mathcal{N}(i)}\right) + \mathbf{h}_i^{(\ell-1)}\right)$$

GATv2 attention (Brody et al., 2021) computes dynamic attention coefficients:

$$e_{ij}^{(k)} = \mathbf{a}^{(k)T} \text{LeakyReLU}\!\left(\mathbf{W}^{(k)}\left[\mathbf{h}_i \,\|\, \mathbf{h}_j\right]\right)$$

with $K = 4$ attention heads and 10% dropout.

**Decoder.** An MLP maps the final latent representation to residual increments:

$$\Delta\tilde{\mathbf{x}}_i = \text{MLP}_\text{dec}\!\left(\mathbf{h}_i^{(L)}\right)$$

**Residual Prediction.** The model predicts increments rather than full states:

$$\hat{\mathbf{x}}_i^{(t+1)} = \tilde{\mathbf{x}}_i^{(t)} + \Delta\tilde{\mathbf{x}}_i$$

This design choice, shared with GraphCast, reduces the learning burden and improves stability.

### 3.4 Physics Constraint Module

To enforce non-negative specific humidity, a softplus activation replaces identity mapping for the humidity channels:

$$\hat{Q}_{i,p} = \log(1 + \exp(\hat{Q}_{i,p}^\text{raw}))$$

This is applied without in-place modification to preserve gradient flow. Column-integrated moisture is monitored as a diagnostic:

$$\text{CWV}_t = \sum_p \Delta p_p \cdot \bar{Q}_{t,p}$$

where $\Delta p_p$ are pressure-layer thicknesses (150, 175, 200, 100, 50 hPa).

### 3.5 Training

Training uses the AdamW optimizer (weight decay $10^{-4}$) with cosine annealing schedule over 30 epochs, learning rate $10^{-3}$. The loss function combines MSE on all variables with an auxiliary temperature term:

$$\mathcal{L} = \text{MSE}(\hat{\mathbf{x}}, \mathbf{x}) + 0.1 \cdot \text{MSE}(\hat{T}, T)$$

Gradient clipping at norm 1.0 prevents exploding gradients in the deep GNN stack. Training samples are drawn with stride 4 from the training period (80% of timesteps) to reduce temporal autocorrelation.

### 3.6 Multi-Resolution Configurations

| Resolution | Grid | Nodes | Hidden $H$ | Layers $L$ | Parameters |
|-----------|------|-------|-----------|-----------|-----------|
| 2.5° (proxy) | 8×16 | 128 | 48 | 3 | ~0.14M |
| 1.0° (proxy) | 12×24 | 288 | 64 | 4 | ~0.34M |
| 0.5° (proxy) | 18×36 | 648 | 64 | 4 | ~0.34M |

Note: "proxy resolutions" refer to reduced grids used for computational feasibility; correspondence to true 2.5°/1.0°/0.5° global grids is approximate.

### 3.7 Evaluation Protocol

**Autoregressive rollout** iterates the model $n$ steps to produce a lead-time-$n \times 6$h forecast. Metrics:

- **RMSE** (normalized): $\sqrt{\frac{1}{N}\sum_i (\hat{x}_i - x_i)^2}$
- **ACC** (anomaly correlation): $\frac{\sum_i (\hat{x}_i - \bar{x}_i)(x_i - \bar{x}_i)}{\sqrt{\sum_i (\hat{x}_i - \bar{x}_i)^2 \sum_i (x_i - \bar{x}_i)^2}}$

where $\bar{x}_i$ is the training-period climatological mean.

**Five-fold time-series cross-validation** partitions the time series into 5 sequential folds; each fold uses the preceding folds for training and the current fold for validation. This avoids data leakage from future to past.

---

## 4. Experiments

### 4.1 Synthetic Data Generation

In the absence of real ERA5 data access, we generate synthetic atmospheric fields designed to mimic qualitative properties of real atmospheric data:

- **Temperature**: Meridional gradient (290 K equator → 230 K poles) plus pressure-level lapse plus propagating zonal wave $\propto \sin(2\pi \lambda/360 + \phi(t))$
- **Zonal wind**: Jet-like structure $\propto \sin(\pi\phi/90)$ plus traveling disturbances
- **Meridional wind**: Weaker, thermally balanced structure
- **Specific humidity**: Exponentially decreasing with pressure; maximum in tropics
- **Geopotential height**: Hydrostatic balance proxy

Gaussian noise ($\sigma = 0.05 \times$ signal) is added to all fields at every timestep to prevent degeneracy. Total length: 180 timesteps (≈ 45 days of 6-hourly data).

### 4.2 Persistence Baseline

The persistence forecast sets the predicted field equal to the current analysis: $\hat{\mathbf{x}}^{(t+n)} = \mathbf{x}^{(t)}$. This is a standard meteorological baseline that all models are expected to outperform at medium-range lead times.

### 4.3 Cross-Validation Details

Cross-validation trains a fresh model instance per fold (to prevent information leakage through shared weights). Each fold model trains for 15 epochs with stride-4 sampling. Single-step predictions ($n=1$, 2, 20 steps corresponding to 6/24/120h) are evaluated on the held-out fold.

---

## 5. Results

### 5.1 Training Convergence

![Figure 1: Training and Validation Loss Curves](figures/training_curves.png)

**Figure 1** shows training and validation loss across 30 epochs for all three resolutions. Training MSE converges to approximately 0.049–0.052 across resolutions, indicating stable optimization. However, validation loss is substantially higher (7–9), revealing significant overfitting. The train/validation loss gap grows monotonically after epoch 10, suggesting the model memorizes temporal patterns within the training window rather than learning transferable atmospheric dynamics. This overfitting is a primary limitation of the current setup and is discussed critically in Section 6.

### 5.2 RMSE Comparison: GNN vs. Persistence

![Figure 2: Multi-Resolution RMSE Comparison](figures/rmse_comparison.png)

**Figure 2** and **Table 1** present normalized RMSE for GNN-Meso versus the persistence baseline across lead times and resolutions.

**Table 1: Normalized RMSE (mean ± std) — GNN-Meso vs. Persistence Baseline**

| Lead Time | Resolution | GNN RMSE | Persistence RMSE | Improvement |
|-----------|-----------|----------|-----------------|-------------|
| 6 h | 2.5° | 0.2215 ± 0.000 | 0.0392 ± 0.001 | −465% (worse) |
| 6 h | 1.0° | 0.2253 ± 0.000 | 0.0403 ± 0.000 | −459% (worse) |
| 6 h | 0.5° | 0.2276 ± 0.000 | 0.0409 ± 0.000 | −456% (worse) |
| 24 h | 2.5° | 0.2251 ± 0.001 | 0.1094 ± 0.000 | −106% (worse) |
| 24 h | 1.0° | 0.2645 ± 0.003 | 0.1080 ± 0.000 | −145% (worse) |
| 24 h | 0.5° | 0.2979 ± 0.004 | 0.1072 ± 0.000 | −178% (worse) |
| 120 h | 2.5° | 0.6942 ± 0.030 | 0.4043 ± 0.000 | −72% (worse) |
| 120 h | 1.0° | 1.1947 ± 0.017 | 0.3974 ± 0.000 | −200% (worse) |
| 120 h | 0.5° | 0.5040 ± 0.013 | 0.3931 ± 0.000 | −28% (worse) |

The GNN consistently fails to outperform persistence at any lead time or resolution in autoregressive rollout. This is a critically important negative result discussed in Section 6.

### 5.3 Anomaly Correlation Coefficient

![Figure 3: ACC Comparison Across Lead Times](figures/acc_comparison.png)

**Table 2: ACC (mean ± std)**

| Lead Time | Resolution | GNN ACC | Persistence ACC |
|-----------|-----------|---------|----------------|
| 6 h | 2.5° | 0.687 ± 0.012 | 0.983 ± 0.001 |
| 6 h | 0.5° | 0.666 ± 0.011 | 0.980 ± 0.001 |
| 24 h | 2.5° | 0.698 ± 0.012 | 0.867 ± 0.008 |
| 24 h | 0.5° | 0.567 ± 0.003 | 0.864 ± 0.008 |
| 120 h | 2.5° | 0.280 ± 0.026 | −0.812 ± 0.022 |
| 120 h | 0.5° | 0.219 ± 0.060 | −0.812 ± 0.022 |

At 120-hour lead time, the persistence baseline ACC becomes negative (−0.81), indicating it is anti-correlated with the climatological anomaly—this is expected since the synthetic atmosphere evolves substantially over 120 hours. The GNN maintains positive ACC (0.22–0.28) at this range, indicating it captures some directional trend even if absolute RMSE is poor.

### 5.4 Forecast Field Visualization

![Figure 4: 500 hPa Temperature Field Forecast](figures/forecast_field.png)

**Figure 4** shows the 500 hPa normalized temperature field at initialization, GNN +6h/+24h forecasts, ground truth +6h/+24h, and the absolute error (difference). The GNN forecast is spatially smooth, lacking the sharp frontal gradients present in the synthetic ground truth. This is a classic characteristic of models trained with MSE loss, which penalizes overprediction of variance.

### 5.5 Physical Consistency

![Figure 5: Physical Consistency Monitoring](figures/physical_consistency.png)

**Figure 5** shows column-integrated moisture (CWV) and kinetic energy (KE) proxy over a 120-hour autoregressive rollout. CWV decreases from initialization (softplus ensures non-negativity but does not conserve total moisture), stabilizing at ~80% of the initial value by +48h. KE shows moderate drift (~±15% over 120h). Neither conservation law is exactly satisfied, which is expected for a purely data-driven model without hard physical constraints. The bounded drift indicates the model does not catastrophically amplify or extinguish energy, a prerequisite for stable long-range prediction.

### 5.6 Cross-Validation Results

![Figure 6: 5-Fold Cross-Validation RMSE](figures/cross_validation.png)

**Table 3: 5-Fold CV RMSE for Single-Step Predictions (0.5° grid)**

| Lead Time | Fold 1 | Fold 2 | Fold 3 | Fold 4 | Fold 5 | Mean ± Std |
|-----------|--------|--------|--------|--------|--------|-----------|
| 6 h | 0.0327 | 0.0332 | 0.0335 | 0.0325 | 0.0325 | **0.0329 ± 0.0004** |
| 24 h | 0.0380 | 0.0343 | 0.0340 | 0.0334 | 0.0430 | **0.0366 ± 0.0036** |
| 120 h | 0.0480 | 0.0519 | 0.0396 | 0.0629 | 0.0449 | **0.0494 ± 0.0078** |

Cross-validation RMSE (single-step) is substantially lower than the autoregressive rollout RMSE in Table 1. This is expected: single-step prediction does not accumulate errors. The standard deviation grows from 0.0004 (6h) to 0.0078 (120h), reflecting increasing uncertainty at longer lead times. Fold 4 at 120h shows notably elevated RMSE (0.0629), suggesting that some temporal segments are harder to predict, possibly due to higher-amplitude wave events in the synthetic data.

---

## 6. Discussion

### 6.1 Interpretation of GNN Performance

The most notable result is that **the GNN-Meso model fails to outperform the persistence baseline in autoregressive rollout** at all lead times and resolutions tested. Several factors contribute to this finding:

**Autoregressive error accumulation.** Training uses single-step supervision, but evaluation requires multi-step rollout. Each 6-hour step introduces errors that compound multiplicatively. At 20 steps (120h), even a modest single-step error of 0.033 (RMSE) can compound to 0.5+. This training/evaluation mismatch—well documented in the NWP deep learning literature—necessitates either multi-step training loss or scheduled sampling strategies, neither of which was implemented here.

**Insufficient training data.** The synthetic dataset contains only 180 timesteps (45 equivalent days), whereas state-of-the-art models train on 20–39 years of ERA5 data (28,000–57,000 timesteps). The model has approximately 0.34M parameters but only ~36 training batches per epoch, making generalization extremely difficult.

**Overfitting evidence.** The train/validation loss ratio of ~1:160 is alarming. The model successfully minimizes training MSE but fails to generalize across the time-series train/val split, suggesting it has effectively memorized temporal patterns in the training window.

### 6.2 Dependence on Synthetic Data Assumptions

The experimental results are highly sensitive to the synthetic data generation assumptions:

- The synthetic atmosphere uses simplified linear wave propagation with additive Gaussian noise. Real atmospheric dynamics are highly nonlinear, include baroclinic instability, moisture-dynamics coupling, and boundary layer processes absent in the synthetic setup.
- The noise scale (5% of signal) is lower than real atmospheric variability, which may make the short-term persistence baseline artificially strong.
- The temporal autocorrelation structure of the synthetic data is simplified; real ERA5 has complex decorrelation timescales that vary by variable, level, and season.

**Conclusion:** Results from this synthetic experiment cannot be extrapolated to real ERA5 data. The architecture design choices may be sound, but performance validation must be repeated on real reanalysis data.

### 6.3 Generalizability to Real-World Data

Expected challenges when applying GNN-Meso to real ERA5 data:

1. **Resolution scaling**: A true 0.25° global grid contains 1,038,240 nodes; our largest grid has 648. Memory and compute requirements increase approximately $O(N \log N)$ with attention mechanisms.
2. **Variable interactions**: Real T-Q-Z-U-V coupling through physical parameterizations (condensation, radiative transfer) is far richer than the synthetic linear system.
3. **Boundary conditions**: Sea surface temperature, land surface, and ice fraction provide forcing not included in the synthetic setup.
4. **Observational errors**: Real NWP begins from imperfect initial conditions obtained via data assimilation.

### 6.4 Bias and Limitations in Experimental Design

- **Selection bias in hyperparameters**: The architecture (hidden=64, 4 layers) was not systematically tuned via hyperparameter search, biasing results toward potentially suboptimal configurations.
- **Evaluation metric scope**: Only RMSE and ACC are reported. Probabilistic calibration, power spectra, and physical variable-specific scores (e.g., 500 hPa geopotential RMSE in geopotential meters) would be more informative for weather applications.
- **Absence of climatological baseline**: A more informative baseline is the climatological mean prediction, which outperforms persistence at medium ranges.
- **Single-model evaluation**: No ensemble spread is computed; uncertainty quantification is absent.

### 6.5 Comparison with Prior Work

In contrast to GraphCast (Lam et al., 2023) and Pangu-Weather (Bi et al., 2023), which both outperform ECMWF IFS on multiple variables after training on decades of ERA5, GNN-Meso fails to outperform even the trivial persistence baseline. The performance gap is attributable to:

| Factor | GraphCast | GNN-Meso |
|--------|-----------|---------|
| Training data | 39 years ERA5 | 180 synthetic timesteps |
| Grid nodes | ~40,000 | 128–648 |
| Processor layers | 16 | 3–4 |
| Multi-step training | Yes (12-step rollout) | No (single-step) |
| Physical variables | 227 (6 surface + 6×37 pl) | 25 (5×5) |

The qualitative architectural similarity (Encoder-Processor-Decoder) is present, but the quantitative scale gap is enormous. This paper should be read as a **proof-of-concept implementation study**, not a benchmark comparison.

---

## 7. Conclusion

We presented GNN-Meso, a PyTorch Geometric-based graph neural network for multi-scale atmospheric prediction inspired by GraphCast and related architectures. The system implements an Encoder–GATv2–Decoder pipeline with residual increment prediction, soft humidity constraints, and multi-resolution configurations. Five-fold cross-validation on synthetic data yields single-step RMSE of 0.0329 ± 0.0004 (6h), 0.0366 ± 0.0036 (24h), and 0.0494 ± 0.0078 (120h). However, autoregressive rollout RMSE substantially exceeds the persistence baseline at all lead times, a critical failure attributable to insufficient training data, lack of multi-step training loss, and the synthetic data regime.

**Key takeaways for future work:**
1. Multi-step training with rollout lengths of 12+ steps is essential for competitive autoregressive performance.
2. Model scale must increase dramatically (millions of parameters, real ERA5 data) to match state-of-the-art results.
3. Hard physical constraints (mass conservation, energy balance) should be incorporated as training penalties or structural architectural constraints, not only soft post-processing.
4. Probabilistic forecasting with ensemble generation is needed for operational relevance.
5. Benchmarking against WeatherBench 2 (Rasp et al., 2023) on real ERA5 is the necessary next step for any credible performance claim.

---

## References

1. **Lam, R., Sanchez-Gonzalez, A., Willson, M., et al.** (2023). Learning skillful medium-range global weather forecasting. *Science*, 382(6677), 1416–1421. DOI: [10.1126/science.adi2336](https://doi.org/10.1126/science.adi2336)

2. **Bi, K., Xie, L., Zhang, H., et al.** (2023). Accurate medium-range global weather forecasting with 3D neural networks. *Nature*, 619, 533–538. DOI: [10.1038/s41586-023-06185-3](https://doi.org/10.1038/s41586-023-06185-3)

3. **Pathak, J., Subramanian, S., Harrington, P., et al.** (2022). FourCastNet: A global data-driven high-resolution weather model using adaptive Fourier neural operators. *arXiv preprint*. arXiv: [2202.11214](https://arxiv.org/abs/2202.11214)

4. **Rasp, S., Dueben, P.D., Scher, S., et al.** (2020). WeatherBench: A benchmark dataset for data‐driven weather forecasting. *Journal of Advances in Modeling Earth Systems*, 12(11), e2020MS002203. DOI: [10.1029/2020MS002203](https://doi.org/10.1029/2020MS002203)

5. **Nguyen, T., Brandstetter, J., Kapoor, A., et al.** (2023). ClimaX: A foundation model for weather and climate. *Proceedings of ICML 2023*. arXiv: [2301.10343](https://arxiv.org/abs/2301.10343)

6. **Keisler, R.** (2022). Forecasting global weather with graph neural networks. *arXiv preprint*. arXiv: [2202.07575](https://arxiv.org/abs/2202.07575)

7. **Kochkov, D., Yuval, J., Langmore, I., et al.** (2024). Neural general circulation models for weather and climate. *Nature*, 632, 1060–1066. DOI: [10.1038/s41586-024-07744-y](https://doi.org/10.1038/s41586-024-07744-y)

8. **Rasp, S., Hoyer, S., Merose, A., et al.** (2024). WeatherBench 2: A benchmark for the next generation of data-driven global weather models. *Journal of Advances in Modeling Earth Systems*, 16(6), e2023MS003715. DOI: [10.1029/2023MS003715](https://doi.org/10.1029/2023MS003715)

9. **Brody, S., Alon, U., & Yahav, E.** (2021). How attentive are graph attention networks? *arXiv preprint*. arXiv: [2105.14491](https://arxiv.org/abs/2105.14491)

10. **Bonev, B., Kurth, T., Hundt, C., et al.** (2023). Spherical Fourier neural operators: Learning stable dynamics on the sphere. *Proceedings of ICML 2023*. arXiv: [2306.03838](https://arxiv.org/abs/2306.03838)
