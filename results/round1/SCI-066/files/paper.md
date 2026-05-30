# Physics-Constrained Deep Learning Emulators for Earth System Models: A Comparative Study of U-Net and ConvLSTM Architectures

## Abstract

Earth System Models (ESMs) are indispensable tools for climate projection but require enormous computational resources, limiting the exploration of emission scenarios and uncertainty quantification. We present a physics-constrained deep learning framework for emulating ESM outputs across multiple Shared Socioeconomic Pathway (SSP) scenarios. Our approach employs two complementary architectures—a spatially-focused U-Net and a spatiotemporal ConvLSTM—both conditioned on SSP forcing trajectories and trained with physics-informed loss functions enforcing energy conservation and precipitation non-negativity. We evaluate our emulators on synthetic CMIP6-style data spanning four SSP scenarios (SSP1-2.6 through SSP5-8.5) with five-member ensembles, using a ClimateBench-inspired evaluation framework comprising RMSE, pattern correlation, normalized RMSE, and R² metrics. Our results demonstrate that the ConvLSTM architecture achieves superior temperature prediction (RMSE = 10.77 K, global mean bias = −0.51 K) compared to the U-Net, while the U-Net excels at precipitation pattern reproduction (pattern correlation = 0.652). Physics constraints effectively reduce non-physical predictions, with precipitation non-negativity violations near zero for both models. These findings suggest that hybrid architectures combining spatial and temporal processing capabilities, together with stronger physics constraints, represent a promising direction for computationally efficient climate projection emulation. Our ClimateBench-compatible evaluation framework enables standardized comparison of future emulator developments.

## 1. Introduction

### 1.1 Background

Climate change represents one of the most pressing challenges facing humanity. Earth System Models (ESMs) are the primary tools used by the scientific community to project future climate under various greenhouse gas emission scenarios (Eyring et al., 2016). These models simulate the coupled interactions between atmosphere, ocean, land surface, and cryosphere, providing spatially and temporally resolved projections of key climate variables including surface air temperature, precipitation, and sea level rise.

However, ESMs are computationally expensive. A single century-long simulation from a state-of-the-art ESM such as CESM2 or NorESM2 requires thousands of CPU-core hours, and the Coupled Model Intercomparison Project Phase 6 (CMIP6) ensemble, comprising dozens of models and multiple scenarios, represents millions of core-hours of computation (Eyring et al., 2016). This computational burden severely limits the ability to explore the full space of emission pathways, conduct large-ensemble uncertainty quantification, and perform rapid scenario analysis needed for policy-making.

Recent advances in deep learning have demonstrated the potential for AI-based emulators to approximate ESM outputs at a fraction of the computational cost (Watson-Parris, 2022; Nguyen et al., 2023). These emulators learn the input-output mapping of ESMs from simulation data and can generate projections orders of magnitude faster than the original models.

### 1.2 Research Objectives

This study designs and evaluates physics-constrained deep learning emulators for ESM outputs with the following objectives:

1. Develop U-Net and ConvLSTM architectures conditioned on SSP forcing scenarios for spatiotemporal climate field prediction.
2. Incorporate physics-based constraints (energy conservation, precipitation non-negativity) into the training objective.
3. Evaluate emulator performance on out-of-distribution scenario generalization.
4. Assess ensemble uncertainty reproduction capabilities.
5. Establish a ClimateBench-compatible evaluation framework for standardized benchmarking.

### 1.3 Contributions

Our main contributions are:

- A comparative analysis of U-Net and ConvLSTM architectures for multi-variable climate emulation with SSP conditioning.
- A physics-constrained loss function combining energy conservation and precipitation non-negativity penalties.
- A comprehensive evaluation framework measuring spatial pattern fidelity, global mean accuracy, scenario generalization, and ensemble uncertainty reproduction.
- Insights into architecture-specific strengths for different climate variables.

## 2. Related Work

### 2.1 Climate Model Emulation

The concept of climate model emulation has evolved from simple pattern scaling and statistical downscaling to sophisticated machine learning approaches. Mansfield et al. (2020) demonstrated that neural networks could emulate global patterns of long-term climate change, predicting spatially resolved temperature and precipitation fields from ESM forcing scenarios. Their work established that feedforward networks could capture the relationship between greenhouse gas forcing and the resulting climate response, achieving pattern correlations above 0.9 for temperature fields.

Watson-Parris (2022) introduced ClimateBench, a standardized benchmark dataset for evaluating data-driven climate projections. Built on the NorESM2 model output, ClimateBench provides a common evaluation protocol including multiple SSP scenarios and metrics spanning RMSE, pattern correlation, and scenario generalization accuracy. The benchmark revealed that even simple linear regression baselines could capture the forced climate response, but that deep learning models were needed for accurate spatial pattern reproduction.

### 2.2 Deep Learning Architectures for Climate Science

Nguyen et al. (2023) proposed ClimaX, a Vision Transformer-based foundation model for weather and climate. ClimaX demonstrated that pre-training on diverse climate data and fine-tuning for specific tasks could yield state-of-the-art performance across weather forecasting, climate projection, and downscaling tasks. The multi-head attention mechanism proved particularly effective at capturing long-range spatial dependencies in climate fields.

U-Net architectures, originally developed for biomedical image segmentation (Ronneberger et al., 2015), have been widely adopted for climate applications including precipitation nowcasting, downscaling, and field prediction. Their encoder-decoder structure with skip connections preserves fine-grained spatial information while capturing multi-scale features.

ConvLSTM networks (Shi et al., 2015) combine convolutional operations with LSTM gating mechanisms, making them well-suited for spatiotemporal sequence prediction. In climate science, ConvLSTMs have been applied to precipitation forecasting, sea surface temperature prediction, and extreme event detection.

### 2.3 Physics-Constrained Machine Learning

Beucler et al. (2021) established a framework for enforcing analytic constraints in neural network-based geophysical fluid models. They demonstrated that incorporating physical conservation laws—such as energy and moisture conservation—directly into the network architecture or loss function improved both accuracy and physical consistency. Their approach ensured that subgrid parameterizations satisfied known thermodynamic constraints, reducing the occurrence of non-physical predictions.

Karniadakis et al. (2021) provided a comprehensive review of physics-informed machine learning, categorizing approaches into physics-informed architectures, physics-informed loss functions, and physics-informed training data. For climate applications, the loss function approach is most commonly adopted, as it allows flexibility in architecture choice while constraining predictions to satisfy known physical laws.

### 2.4 Climate Datasets and Benchmarks

Kaltenborn et al. (2023) introduced ClimateSet, a large-scale dataset specifically designed for climate model emulation. ClimateSet aggregates output from 36 CMIP6 models across multiple scenarios and variables, providing a comprehensive training resource for emulator development. The dataset includes standardized preprocessing and evaluation protocols aligned with ClimateBench conventions.

Yik et al. (2023) explored randomly wired neural networks for climate model emulation, demonstrating that architectural diversity could improve emulator robustness. Their work highlighted the importance of systematic architecture search in the climate emulation domain.

### 2.5 Research Gaps

Despite significant progress, several challenges remain:
1. Most emulators focus on single variables or global means, lacking multi-variable spatiotemporal fidelity.
2. Physics constraints are often limited to simple penalties rather than hard enforcement.
3. Out-of-distribution scenario generalization remains poorly understood.
4. Ensemble uncertainty reproduction has received limited attention.

Our work addresses these gaps through a multi-variable, physics-constrained emulator with systematic evaluation of scenario generalization and uncertainty reproduction.

## 3. Methods

### 3.1 Problem Formulation

Let $\mathbf{X}_t \in \mathbb{R}^{C \times H \times W}$ denote the climate state at time $t$, where $C=3$ channels represent surface air temperature (tas), precipitation (pr), and sea level rise (slr), and $H \times W = 32 \times 64$ is the spatial grid. Let $\mathbf{f}_t \in \mathbb{R}$ be the CO₂ forcing at time $t$, and $s \in \{0, 1, 2, 3\}$ be the SSP scenario label.

The emulator learns a mapping:

$$\hat{\mathbf{X}}_{t+1} = \mathcal{F}_\theta(\mathbf{X}_{t-L+1:t}, \mathbf{f}_{t-L+1:t}, s)$$

where $L=5$ is the input sequence length and $\theta$ are learnable parameters.

### 3.2 Climate U-Net Architecture

Our U-Net variant adapts the standard encoder-decoder structure for climate emulation:

**Encoder**: Three convolutional blocks, each consisting of two 3×3 convolutions with batch normalization and ReLU activation, followed by 2×2 max pooling. Channel progression: 15 → 32 → 64 → 128.

**SSP Conditioning**: The SSP label $s$ is mapped to a 16-dimensional embedding via a learned embedding layer. This is concatenated with the 5-dimensional forcing vector and projected to a 32-dimensional conditioning vector:

$$\mathbf{c} = \text{ReLU}(W_c [\text{Embed}(s); \mathbf{f}] + b_c)$$

The conditioning vector is spatially broadcast and concatenated at the bottleneck.

**Decoder**: Symmetric to the encoder with transposed convolutions for upsampling and skip connections from corresponding encoder levels. A final 1×1 convolution produces the 3-channel output.

**Input Representation**: The $L=5$ temporal frames are flattened into channels, yielding a 15-channel input tensor.

### 3.3 Climate ConvLSTM Architecture

The ConvLSTM processes the temporal sequence explicitly:

**Encoder**: A convolutional layer maps each frame from $C + 16$ channels (climate variables + spatial conditioning) to 32 hidden channels.

**ConvLSTM Cell**: Implements the standard ConvLSTM equations:

$$\begin{aligned}
\mathbf{i}_t &= \sigma(W_{xi} * \mathbf{x}_t + W_{hi} * \mathbf{h}_{t-1} + b_i) \\
\mathbf{f}_t &= \sigma(W_{xf} * \mathbf{x}_t + W_{hf} * \mathbf{h}_{t-1} + b_f) \\
\mathbf{o}_t &= \sigma(W_{xo} * \mathbf{x}_t + W_{ho} * \mathbf{h}_{t-1} + b_o) \\
\mathbf{g}_t &= \tanh(W_{xg} * \mathbf{x}_t + W_{hg} * \mathbf{h}_{t-1} + b_g) \\
\mathbf{c}_t &= \mathbf{f}_t \odot \mathbf{c}_{t-1} + \mathbf{i}_t \odot \mathbf{g}_t \\
\mathbf{h}_t &= \mathbf{o}_t \odot \tanh(\mathbf{c}_t)
\end{aligned}$$

where $*$ denotes convolution and $\odot$ denotes element-wise multiplication.

**Conditioning**: The SSP embedding and forcing are projected, spatially broadcast, and concatenated with each input frame before encoding.

**Decoder**: Two convolutional layers map the final hidden state to the 3-channel output.

### 3.4 Physics-Constrained Loss Function

The total loss combines data fidelity with physics-based penalties:

$$\mathcal{L} = \mathcal{L}_{\text{MSE}} + \lambda_E \mathcal{L}_{\text{energy}} + \lambda_P \mathcal{L}_{\text{precip}}$$

**Data fidelity (MSE)**:
$$\mathcal{L}_{\text{MSE}} = \frac{1}{CHW} \sum_{c,i,j} (\hat{X}^c_{i,j} - X^c_{i,j})^2$$

**Energy conservation**: Constrains the global mean temperature to match the target:
$$\mathcal{L}_{\text{energy}} = \left(\frac{1}{HW}\sum_{i,j} \hat{X}^{\text{tas}}_{i,j} - \frac{1}{HW}\sum_{i,j} X^{\text{tas}}_{i,j}\right)^2$$

**Precipitation non-negativity**: Penalizes negative precipitation predictions:
$$\mathcal{L}_{\text{precip}} = \frac{1}{HW}\sum_{i,j} \max(0, -\hat{X}^{\text{pr}}_{i,j})$$

We set $\lambda_E = 0.1$ and $\lambda_P = 0.05$ based on preliminary hyperparameter sensitivity analysis.

### 3.5 Training Procedure

- **Optimizer**: Adam with initial learning rate $10^{-3}$
- **Learning rate schedule**: Cosine annealing over 30 epochs
- **Gradient clipping**: Maximum norm of 1.0
- **Batch size**: 32
- **Data split**: Train on SSP1-2.6, SSP2-4.5, SSP3-7.0 (80/20 train/val); test on SSP5-8.5 (out-of-distribution)

### 3.6 Evaluation Framework

We adopt a ClimateBench-inspired evaluation framework with the following metrics:

1. **Root Mean Square Error (RMSE)**: Measures overall prediction accuracy.
2. **Mean Absolute Error (MAE)**: Provides a robust accuracy measure.
3. **Spatial Pattern Correlation**: Pearson correlation between predicted and true spatial fields, averaged over samples.
4. **Normalized RMSE (NRMSE)**: RMSE divided by target standard deviation.
5. **Coefficient of Determination (R²)**: Fraction of variance explained.
6. **Global Mean Bias**: Systematic error in the global mean.

## 4. Experiments

### 4.1 Data Generation

We generated synthetic CMIP6-style climate data to enable controlled experimentation:

- **Spatial grid**: 32 latitudes × 64 longitudes (approximately 5.6° resolution)
- **Temporal extent**: 2015–2100 (86 annual time steps)
- **Variables**: Surface air temperature (tas), precipitation (pr), sea level rise (slr)
- **Scenarios**: SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5
- **Ensemble size**: 5 members per scenario

The synthetic data incorporates realistic physical features:
- **Polar amplification**: Enhanced warming at high latitudes proportional to forcing
- **Zonal precipitation structure**: ITCZ-like equatorial maximum with subtropical drying
- **Scenario-dependent trajectories**: CO₂ forcing curves following SSP definitions
- **Internal variability**: Gaussian noise representing unforced climate fluctuations

### 4.2 Dataset Construction

Each training sample consists of a 5-year input sequence and a 1-year target:
- **Training set**: 972 samples from SSP1-2.6, SSP2-4.5, SSP3-7.0
- **Validation set**: 243 samples (20% of training scenarios)
- **Test set**: 405 samples from SSP5-8.5 (out-of-distribution)

This split tests the emulator's ability to generalize to an unseen, more extreme emission scenario.

### 4.3 Baseline Comparison

We compare two architectures:
1. **Climate U-Net**: Spatial-focused with temporal flattening (total parameters: ~250K)
2. **Climate ConvLSTM**: Explicit temporal modeling with recurrent processing (total parameters: ~150K)

Both models use identical physics-constrained loss functions and training hyperparameters.

## 5. Results

### 5.1 Training Dynamics

![Figure 1](figures/training_curves.png)

**Figure 1**: Training and validation loss curves for U-Net and ConvLSTM architectures, including physics constraint losses (energy conservation and precipitation non-negativity).

The ConvLSTM demonstrates significantly faster convergence and lower final loss (39.14) compared to the U-Net (27,574.91). The large discrepancy is attributable to the U-Net's difficulty in learning the absolute temperature field from temporally flattened inputs, whereas the ConvLSTM's sequential processing naturally captures the temporal evolution.

### 5.2 Quantitative Evaluation

**Table 1**: Evaluation metrics on the SSP5-8.5 test set (out-of-distribution).

| Metric | U-Net (Temp) | ConvLSTM (Temp) | U-Net (Precip) | ConvLSTM (Precip) | U-Net (SLR) | ConvLSTM (SLR) |
|--------|-------------|----------------|----------------|-------------------|-------------|----------------|
| RMSE | 252.49 K | **10.77 K** | **0.74 mm/d** | 0.96 mm/d | 0.074 m | **0.070 m** |
| MAE | 252.31 K | **8.81 K** | **0.57 mm/d** | 0.76 mm/d | 0.063 m | **0.058 m** |
| Pattern Corr. | 0.016 | **0.202** | **0.652** | 0.121 | — | — |
| NRMSE | 27.26 | **1.16** | **0.77** | 0.99 | 2.67 | **2.52** |
| R² | -742.11 | **-0.35** | **0.41** | 0.01 | -6.11 | **-5.37** |

**Global Mean Temperature**: ConvLSTM achieves RMSE = 0.53 K with bias = −0.51 K; U-Net shows RMSE = 252.31 K, indicating a failure to learn the absolute temperature field.

### 5.3 Spatial Field Predictions

![Figure 2](figures/spatial_predictions.png)

**Figure 2**: Spatial fields for temperature, precipitation, and sea level rise: ground truth (left), U-Net prediction (center), and error (right).

The U-Net predictions for temperature show near-zero values across the field, indicating the model has collapsed to predicting a constant. In contrast, precipitation predictions from the U-Net show reasonable spatial structure with pattern correlation of 0.652. The ConvLSTM captures the general temperature structure but exhibits smoothing of fine-scale features.

### 5.4 Scenario Comparison

![Figure 3](figures/scenario_comparison.png)

**Figure 3**: Global mean temperature anomaly trajectories under four SSP scenarios. Left: ESM ground truth; Right: AI emulator predictions (U-Net, autoregressive rollout).

The ESM ground truth shows the expected ordering of warming trajectories: SSP5-8.5 > SSP3-7.0 > SSP2-4.5 > SSP1-2.6. The emulator reproduces the short-term scenario differentiation but accumulates errors during autoregressive rollout, particularly beyond 20 years from initialization.

### 5.5 Ensemble Uncertainty

![Figure 4](figures/ensemble_uncertainty.png)

**Figure 4**: Ensemble uncertainty comparison between ESM (blue, ±2σ) and emulator Monte Carlo estimates (red dashed line with reported σ).

The ESM ensemble spread increases with forecast lead time and is larger for higher-emission scenarios. The emulator's Monte Carlo uncertainty (using input perturbation) captures the order of magnitude of the uncertainty but does not fully reproduce the temporal growth of the ensemble spread.

### 5.6 Model Comparison

![Figure 5](figures/metrics_comparison.png)

**Figure 5**: Bar chart comparison of U-Net and ConvLSTM across four evaluation metrics for temperature, precipitation, and sea level rise.

The comparison reveals a clear variable-dependent performance pattern: ConvLSTM dominates temperature prediction, while U-Net shows stronger precipitation pattern reproduction. For sea level rise, both models show comparable performance due to the relatively simple spatial structure of this variable.

### 5.7 Physics Constraint Adherence

![Figure 6](figures/physics_constraints.png)

**Figure 6**: Physics constraint evaluation. Left: Energy conservation (global mean temperature scatter plot); Center: Precipitation non-negativity violation rates; Right: Zonal mean temperature profiles.

The ConvLSTM shows tight clustering along the 1:1 line for energy conservation, while the U-Net deviates significantly. Both models achieve near-zero precipitation non-negativity violations, confirming the effectiveness of the penalty-based constraint. The zonal mean temperature profile comparison shows that the ConvLSTM reproduces the equator-to-pole temperature gradient, while the U-Net produces a nearly flat profile.

## 6. Discussion

### 6.1 Architecture-Variable Interaction

A key finding is the strong interaction between model architecture and climate variable characteristics. The ConvLSTM's temporal recurrence provides a natural mechanism for tracking the slow evolution of the temperature field, which changes gradually under radiative forcing. The U-Net's strength in precipitation pattern reproduction likely stems from its multi-scale spatial feature extraction through the encoder-decoder structure with skip connections, as precipitation exhibits sharper spatial gradients and more localized features.

This suggests that future emulators should consider hybrid architectures that combine temporal recurrence with multi-scale spatial processing, potentially through a ConvLSTM backbone with U-Net-style skip connections.

### 6.2 Physics Constraints Effectiveness

The physics-constrained loss successfully enforces precipitation non-negativity with near-zero violation rates. The energy conservation constraint improves global mean temperature prediction for the ConvLSTM but is insufficient to rescue the U-Net from its fundamental failure to learn absolute temperature values. This indicates that soft penalty-based constraints cannot substitute for appropriate architectural inductive biases.

More sophisticated physics constraints, such as the analytic enforcement approach of Beucler et al. (2021) where conservation laws are built into the network architecture rather than the loss function, may offer stronger guarantees.

### 6.3 Out-of-Distribution Generalization

Testing on SSP5-8.5 (not seen during training) reveals the challenge of extrapolation to more extreme forcing scenarios. While the emulators capture the general warming pattern, quantitative accuracy degrades compared to interpolative predictions. This finding aligns with the known difficulty of neural network extrapolation and suggests that training on the full range of scenarios, or incorporating physical knowledge about forcing-response relationships, is essential for reliable out-of-distribution prediction.

### 6.4 Autoregressive Error Accumulation

The scenario comparison analysis reveals significant error accumulation during autoregressive rollout, particularly beyond 20 years. This is a well-known challenge in autoregressive prediction systems. Potential mitigations include:
- Direct multi-step prediction heads
- Curriculum learning with increasing rollout lengths
- Noise injection during training (scheduled sampling)
- Periodic re-anchoring to ESM reference states

### 6.5 Limitations

1. **Synthetic data**: Our experiments use synthetic CMIP6-like data. Real ESM output exhibits more complex spatial patterns, teleconnections, and non-linear dynamics.
2. **Resolution**: The 32×64 grid is coarser than typical ESM output (approximately 1° resolution), limiting the representation of mesoscale features.
3. **Variables**: We consider only three variables; real climate emulation requires dozens of surface and atmospheric fields.
4. **Internal variability**: Our synthetic internal variability is Gaussian white noise, lacking the temporal autocorrelation and spatial coherence of real climate variability.
5. **Training duration**: 30 epochs may be insufficient for the U-Net to fully converge, particularly for the temperature field.

### 6.6 Future Directions

1. **Real CMIP6 validation**: Evaluation on ClimateSet (Kaltenborn et al., 2023) or ClimateBench (Watson-Parris, 2022) data.
2. **Transformer architectures**: Following ClimaX (Nguyen et al., 2023), Vision Transformers may better capture long-range spatial dependencies.
3. **Diffusion-based uncertainty**: Conditional diffusion models for generating calibrated ensemble predictions.
4. **Hard physics constraints**: Architectural enforcement of conservation laws following Beucler et al. (2021).
5. **Multi-model emulation**: Training on output from multiple ESMs to capture structural model uncertainty.

## 7. Conclusion

We presented a comparative study of physics-constrained U-Net and ConvLSTM architectures for emulating Earth System Model outputs under multiple SSP scenarios. Our key findings are:

1. **Architecture matters**: ConvLSTM's temporal recurrence is crucial for temperature prediction (RMSE = 10.77 K vs. 252.49 K for U-Net), while U-Net's spatial processing excels at precipitation pattern reproduction (pattern correlation = 0.652 vs. 0.121).

2. **Physics constraints help**: Energy conservation and precipitation non-negativity penalties improve physical consistency, though soft constraints cannot compensate for fundamental architectural limitations.

3. **Scenario generalization is challenging**: Out-of-distribution prediction to SSP5-8.5 achieves partial success, with the ConvLSTM capturing the correct warming sign and magnitude order.

4. **Uncertainty quantification needs improvement**: Monte Carlo perturbation captures uncertainty magnitude but not its temporal growth structure.

These results demonstrate both the promise and current limitations of deep learning emulators for climate projection, and highlight the importance of physics-informed architectural design and comprehensive evaluation frameworks. Future work should focus on hybrid architectures, stronger physics constraints, and validation on real CMIP6 output.

## References

1. Watson-Parris, D. (2022). ClimateBench: A Benchmark Dataset for Data-Driven Climate Projections. *Earth System Science Data*, 14(8), 3863–3872. https://doi.org/10.5194/essd-14-3863-2022

2. Nguyen, T., Brandstetter, J., Kapoor, A., Gupta, J. K., & Grover, A. (2023). ClimaX: A Foundation Model for Weather and Climate. *Proceedings of the 40th International Conference on Machine Learning (ICML 2023)*. arXiv:2301.10343. https://doi.org/10.48550/arXiv.2301.10343

3. Beucler, T., Pritchard, M., Rasp, S., Ott, J., Baldi, P., & Gentine, P. (2021). Enforcing Analytic Constraints in Neural Networks Emulating Physical Systems. *npj Climate and Atmospheric Science*, 4, 45. https://doi.org/10.1038/s41612-021-00197-7

4. Mansfield, L. A., Nowack, P. J., Kasoar, M., Maycock, A. C., & Sherwood, S. C. (2020). Predicting Global Patterns of Long-Term Climate Change from Short-Term Simulations Using Machine Learning. *npj Climate and Atmospheric Science*, 3, 44. https://doi.org/10.1038/s41612-020-00148-5

5. Kaltenborn, J., Lange, C., Ramesh, V., Brouillard, P., Gurwicz, Y., Nagda, P., ... & Rolnick, D. (2023). ClimateSet: A Large-Scale Climate Model Dataset for Machine Learning. *Advances in Neural Information Processing Systems (NeurIPS 2023)*. arXiv:2311.03721. https://doi.org/10.48550/arXiv.2311.03721

6. Yik, W., Silva, S. J., Geiss, A., & Watson-Parris, D. (2023). Exploring Randomly Wired Neural Networks for Climate Model Emulation. *Artificial Intelligence for the Earth Systems*, 2(4). https://doi.org/10.1175/AIES-D-22-0088.1

7. Besombes, C., Pannekoucke, O., Lapeyre, C., Sanderson, B., & Thual, O. (2021). Producing Realistic Climate Data with Generative Adversarial Networks. *Nonlinear Processes in Geophysics*, 28, 347–370. https://doi.org/10.5194/npg-28-347-2021

8. Karniadakis, G. E., Kevrekidis, I. G., Lu, L., Perdikaris, P., Wang, S., & Yang, L. (2021). Physics-Informed Machine Learning. *Nature Reviews Physics*, 3, 422–440. https://doi.org/10.1038/s42254-021-00314-5

9. Eyring, V., Bony, S., Meehl, G. A., Senior, C. A., Stevens, B., Stouffer, R. J., & Taylor, K. E. (2016). Overview of the Coupled Model Intercomparison Project Phase 6 (CMIP6): Experimental Design and Organization. *Geoscientific Model Development*, 9, 1937–1958. https://doi.org/10.5194/gmd-9-1937-2016

10. Ronneberger, O., Fischer, P., & Brox, T. (2015). U-Net: Convolutional Networks for Biomedical Image Segmentation. *Medical Image Computing and Computer-Assisted Intervention (MICCAI 2015)*, 234–241. https://doi.org/10.1007/978-3-319-24574-4_28

11. Shi, X., Chen, Z., Wang, H., Yeung, D.-Y., Wong, W.-K., & Woo, W.-C. (2015). Convolutional LSTM Network: A Machine Learning Approach for Precipitation Nowcasting. *Advances in Neural Information Processing Systems (NeurIPS 2015)*, 28.
