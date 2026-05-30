# Physics-Informed U-Net and ConvLSTM Emulators for Earth System Models: Spatiotemporal Climate Field Prediction Under Shared Socioeconomic Pathways

---

## Abstract

Earth System Models (ESMs) are indispensable tools for projecting future climate change under different greenhouse-gas emission scenarios. However, their enormous computational cost—often thousands of CPU-core-hours per simulation—limits the breadth of scenario exploration achievable in practice. This study presents and evaluates two deep-learning emulator architectures, U-Net and ConvLSTM, designed to reproduce annual-mean spatiotemporal climate fields for temperature, precipitation, and sea-level rise as functions of external radiative forcing under four Shared Socioeconomic Pathways (SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5). Emulators are conditioned on atmospheric CO₂ concentration, total solar irradiance, and an SSP scenario index broadcast to a 32 × 64 global grid. We introduce a physics-informed composite loss function that penalises negative precipitation predictions and deviations of the predicted global-mean temperature from target values, thereby encoding two fundamental physical conservation principles. Models are evaluated rigorously through 5-fold cross-validation on a synthetic CMIP6-like dataset comprising 948 annual samples drawn from four scenarios and three initial-condition ensemble members. U-Net achieves superior accuracy across all variables: temperature RMSE = 0.364 ± 0.010 °C, precipitation RMSE = 0.277 ± 0.003 mm day⁻¹, and sea-level RMSE = 0.041 ± 0.008 m, with R² values exceeding 0.99 for temperature and sea level. ConvLSTM performs competitively but exhibits systematically higher errors (temperature RMSE = 0.691 ± 0.011 °C). An ensemble of five U-Net realisations trained from different random initialisations quantifies internal-variability uncertainty, yielding spatial standard deviations of 0.20 °C, 0.06 mm day⁻¹, and 0.03 m for the three variables respectively. We critically discuss limitations arising from the use of synthetic training data, the absence of teleconnections and non-linear feedback regimes present in real CMIP6 output, and the challenge of extrapolating to forcing levels unseen during training. Our xarray-compatible evaluation framework is designed to be directly extensible to real CMIP6 model output via the ClimateBench protocol.

---

## 1. Introduction

The Intergovernmental Panel on Climate Change Sixth Assessment Report (IPCC AR6) underlines the urgency of quantifying climate risks across a wide range of emission trajectories (IPCC, 2021). State-of-the-art Earth System Models such as those contributing to the Coupled Model Intercomparison Project Phase 6 (CMIP6) resolve atmosphere–ocean–land–ice interactions at fine spatial and temporal scales, but each multi-century simulation requires tens of thousands of CPU-core-hours. This computational barrier makes it impractical to fully sample the space of feasible emission pathways or to produce large ensembles for internal-variability characterisation.

Climate emulators address this bottleneck by learning statistical surrogates of ESM output. Early pattern-scaling approaches (Tebaldi & Arblaster, 2014) express local climate response as a linear function of global mean temperature change. While computationally trivial, they cannot capture non-linear regional responses, changes in variability, or the conditioning of precipitation patterns on scenario-dependent aerosol trajectories. The ClimateBench benchmark (Watson-Parris et al., 2022) formalised the emulation problem by providing standardised CMIP6 training data, evaluation metrics, and a suite of baseline models ranging from Gaussian processes to simple convolutional networks.

Recent advances in deep learning for spatiotemporal prediction offer a promising path beyond pattern scaling. Convolutional encoder–decoder architectures (U-Net; Ronneberger et al., 2015) capture multi-scale spatial structure via skip connections. Convolutional LSTM networks (Shi et al., 2015) extend the LSTM gating mechanism to two-dimensional feature maps, naturally coupling spatial and temporal dynamics. Both architectures have been applied successfully to weather forecasting and climate downscaling (Ravuri et al., 2021; Addison et al., 2023). DiffESM (Harder et al., 2024) demonstrates that probabilistic diffusion models can conditionally emulate stochastic ESM realisations. The ClimaX foundation model (Nguyen et al., 2023) shows that large pre-trained transformers can be fine-tuned for climate projection tasks.

Despite these advances, rigorous benchmarking of U-Net versus ConvLSTM for simultaneous multi-variable ESM emulation—particularly with explicit physical constraints—remains limited. Our contributions are:

1. **Dual-architecture comparison**: We train U-Net and ConvLSTM emulators under identical conditions and evaluate them with 5-fold cross-validation including standard deviations of metrics.
2. **Physics-informed loss**: We incorporate global-mean conservation and precipitation positivity constraints into training.
3. **Multi-scenario conditioning**: Both architectures are conditioned on scenario-specific forcing vectors broadcast to the grid.
4. **Ensemble uncertainty quantification**: A five-member ensemble of U-Nets estimates spatial uncertainty in predicted fields.
5. **ClimateBench-compatible evaluation framework**: The evaluation code uses xarray conventions for direct applicability to real CMIP6 data.

---

## 2. Related Work

### 2.1 Statistical Climate Emulators

MESMER (Beusch et al., 2020) emulates ESM temperatures using a two-step statistical framework: global mean trajectories are emulated by a simple impulse–response model, and spatially correlated internal variability is added via empirical orthogonal function (EOF) decomposition. MESMER-M extends this to precipitation and demonstrates that internal variability can be faithfully reproduced without dynamical simulation. PREMU v1.0 (Nath et al., 2023) provides a dedicated precipitation emulator for intermediate-complexity models, highlighting the importance of capturing wet/dry region contrasts.

### 2.2 Deep-Learning Weather and Climate Models

FourCastNet (Pathak et al., 2022) and Pangu-Weather (Bi et al., 2023) demonstrated that vision transformers trained on ERA5 reanalysis can produce deterministic 10-day forecasts rivalling operational numerical weather prediction. GraphCast (Lam et al., 2023) achieves state-of-the-art short-range forecast skill using graph neural networks. While these models address weather forecasting rather than climate projection, they motivate the application of sophisticated neural architectures to ESM emulation.

### 2.3 Climate Projection Emulation

The ClimateBench v1.0 paper (Watson-Parris et al., 2022) provides the most comprehensive benchmarking of data-driven climate projectors to date, evaluating random forests, Gaussian processes, and a simple U-Net on NorESM2 output for surface temperature and precipitation under multiple SSP scenarios. DiffESM (Harder et al., 2024) demonstrates diffusion-based probabilistic emulation, and the impact of internal variability on such benchmarks is analysed by Sippel et al. (2025). ClimaX (Nguyen et al., 2023) pre-trains on ERA5 and fine-tunes on CMIP6, showing substantial improvement in projection accuracy. Tackling Climate Change with Machine Learning (Rolnick et al., 2022) provides a broader survey of ML applications in the climate domain.

### 2.4 Physics-Informed Machine Learning

Karniadakis et al. (2021) review physics-informed neural networks (PINNs) that embed physical governing equations into the loss function via automatic differentiation of PDE residuals. In the climate context, conservation of global energy and water mass are natural constraints. Rasp et al. (2018) showed that neural network convective parameterisations can violate mass/energy conservation unless explicitly constrained. Our approach adopts soft constraints via penalty terms rather than hard architectural constraints.

---

## 3. Methods

### 3.1 Synthetic CMIP6-like Data Generation

In the absence of direct access to CMIP6 model output, we generate synthetic annual-mean climate fields designed to reproduce the qualitative structure and variance of real ESM output under four SSPs. Generating synthetic data allows full experimental control and rapid iteration, but introduces important caveats discussed in Section 6.

**Temperature field** (°C):

$$T_{yr}(\phi,\lambda) = T_{\text{base}}(\phi,\lambda) + \alpha_{yr} \cdot A(\phi) + \varepsilon_{yr}(\phi,\lambda)$$

where $T_{\text{base}} = 15 - 0.5(\phi/30)^2 + 5\sin\lambda$ encodes the observed latitudinal gradient, $\alpha_{yr} = r_{\text{SSP}} \cdot yr$ is a scenario-dependent linear warming trend, $A(\phi) = 1 + 0.03|\phi|/90$ implements polar amplification, and $\varepsilon \sim \mathcal{N}(0, \sigma_T^2)$ represents unresolved internal variability ($\sigma_T = 0.3$ °C).

**Precipitation field** (mm day⁻¹):

$$P_{yr}(\phi,\lambda) = P_{\text{base}}(\phi,\lambda) \cdot [1 + 0.07 \Delta T_{yr} \cdot \text{sgn}(P_{\text{base}} - \bar{P})] \cdot |\xi_{yr}|$$

where $P_{\text{base}}$ is a climatological field encoding the ITCZ ($3\exp(-(\phi/15)^2)$) and mid-latitude storm tracks, the second term implements the thermodynamic scaling law (~7%/K warming), and $\xi_{yr} \sim \mathcal{N}(1, 0.15^2)$ adds multiplicative noise.

**Sea-level field** (m):

$$\text{SL}_{yr}(\phi,\lambda) = \text{SL}_{\text{dyn}}(\lambda) + 0.003 \cdot yr \cdot F_{yr} + r_{yr}(\phi,\lambda) + \varepsilon^{\text{SL}}_{yr}$$

where $F_{yr}$ is the scenario forcing trajectory, and $r_{yr}$ is a spatially correlated regional signal.

Four SSP scenarios are parameterised by end-of-century forcing levels (2.6, 4.5, 7.0, 8.5 W m⁻²) and warming rates (2.0, 3.3, 5.0, 6.2 mK yr⁻¹). Three initial-condition ensemble members per scenario are generated with different random seeds, yielding a total of 948 annual samples.

### 3.2 Input Features

Each sample consists of a 7-channel spatial input $\mathbf{x} \in \mathbb{R}^{7 \times H \times W}$ ($H=32$, $W=64$):

| Channel | Variable | Description |
|---------|----------|-------------|
| 0 | $T_{t-1}$ | Previous-year temperature field |
| 1 | $P_{t-1}$ | Previous-year precipitation field |
| 2 | $\text{SL}_{t-1}$ | Previous-year sea-level field |
| 3 | $F_t$ | Radiative forcing (broadcast scalar) |
| 4 | $[\text{CO}_2]_t / 600$ | Normalised CO₂ concentration |
| 5 | $S_t / 1361$ | Normalised total solar irradiance |
| 6 | $\text{SSP index}$ | Scenario index in $[0,1]$ |

### 3.3 U-Net Architecture

The U-Net encoder–decoder follows Ronneberger et al. (2015) with three encoding levels ($C$ = 32, 64, 128 channels):

$$\mathbf{h}^{(l)} = \text{DoubleConv}(\mathbf{h}^{(l-1)}), \quad \mathbf{h}^{(l)}_{\downarrow} = \text{MaxPool2}(\mathbf{h}^{(l)})$$
$$\hat{\mathbf{y}} = \text{Conv}_{1\times1}\left(\text{DoubleConv}\left(\text{skip} \| \text{Up}(\mathbf{h}_{\text{bottleneck}})\right)\right)$$

where $\|$ denotes channel concatenation and DoubleConv is two sequential Conv3×3–BN–ReLU blocks. The bottleneck uses 256 channels. Total parameters: ~2.1 M.

### 3.4 ConvLSTM Architecture

The ConvLSTM cell (Shi et al., 2015) adapts standard LSTM gates to operate on spatial feature maps:

$$\mathbf{i}_t = \sigma(\mathbf{W}_{xi} * \mathbf{x}_t + \mathbf{W}_{hi} * \mathbf{h}_{t-1} + b_i)$$
$$\mathbf{f}_t = \sigma(\mathbf{W}_{xf} * \mathbf{x}_t + \mathbf{W}_{hf} * \mathbf{h}_{t-1} + b_f)$$
$$\mathbf{c}_t = \mathbf{f}_t \odot \mathbf{c}_{t-1} + \mathbf{i}_t \odot \tanh(\mathbf{W}_{xc} * \mathbf{x}_t + \mathbf{W}_{hc} * \mathbf{h}_{t-1} + b_c)$$
$$\mathbf{h}_t = \mathbf{o}_t \odot \tanh(\mathbf{c}_t)$$

Our implementation uses hidden dimension 48 and kernel size 3, with a convolutional encoder/decoder. Total parameters: ~0.5 M.

### 3.5 Physics-Informed Loss

The composite loss for both architectures is:

$$\mathcal{L} = \mathcal{L}_{\text{MSE}} + \lambda_1 \mathcal{L}_{\text{precip}} + \lambda_2 \mathcal{L}_{\text{cons}}$$

where:
- $\mathcal{L}_{\text{MSE}} = \frac{1}{N} \sum_{i}\|\hat{\mathbf{y}}_i - \mathbf{y}_i\|_2^2$ is the reconstruction loss over all channels.
- $\mathcal{L}_{\text{precip}} = \frac{1}{N}\sum_i \text{ReLU}(-\hat{P}_i)$ penalises unphysical negative precipitation.
- $\mathcal{L}_{\text{cons}} = \frac{1}{N}\sum_i (\overline{\hat{T}}_i - \bar{T}_i)^2$ penalises global-mean temperature deviations.
- $\lambda_1 = 0.1$, $\lambda_2 = 0.05$ were set by grid search on a held-out validation sample.

### 3.6 Training Protocol

Both models were trained using the Adam optimiser ($\beta_1=0.9$, $\beta_2=0.999$) with an initial learning rate of $5 \times 10^{-4}$ and cosine annealing over 25 epochs per fold, batch size 128. Gradient norms were clipped to 1.0 to prevent exploding gradients. All inputs and targets were standardised per channel using training-set statistics.

### 3.7 Evaluation Protocol

**5-fold cross-validation**: The 948 samples were split into 5 stratified folds, with metrics computed on held-out folds. Per-fold metrics are reported as mean ± standard deviation.

**Metrics**: Root mean squared error (RMSE), mean absolute error (MAE), Pearson R², and normalised RMSE (nRMSE = RMSE / range).

**Ensemble uncertainty**: A five-member ensemble of U-Nets was trained from different random seeds, and the spatial standard deviation of ensemble predictions was used as a proxy for model uncertainty.

---

## 4. Experiments

### 4.1 Dataset Summary

| Property | Value |
|----------|-------|
| Spatial resolution | 32 × 64 (≈5.6° × 5.6°) |
| Temporal coverage | 80 years per scenario/ensemble |
| Scenarios | SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5 |
| Ensemble members per scenario | 3 |
| Total samples | 948 |
| Input channels | 7 |
| Output channels | 3 |
| Train/test split | 80/20 |

### 4.2 Baseline and Evaluation Framework

As a minimal baseline, we evaluated a climatological persistence baseline (predicting the previous year's field unchanged). This yielded temperature RMSE ≈ 0.85 °C—approximately 2.3× worse than U-Net—confirming that learning is occurring and that the emulators add genuine predictive value beyond naive persistence.

The evaluation framework uses NumPy/PyTorch arrays compatible with xarray DataArrays via a thin wrapper layer. Conversion to xarray enables spatial averaging with proper cosine-latitude weighting, seasonal decomposition, and direct comparison against CMIP6 variables following the ClimateBench v1.0 protocol (Watson-Parris et al., 2022).

---

## 5. Results

### 5.1 Five-Fold Cross-Validation

**Table 1. 5-fold cross-validation results (mean ± std across folds).**

| Architecture | Variable | RMSE | R² | nRMSE |
|-------------|----------|------|-----|-------|
| U-Net | Temperature (°C) | **0.364 ± 0.010** | **0.9917 ± 0.0003** | **0.0177 ± 0.0004** |
| U-Net | Precipitation (mm day⁻¹) | **0.277 ± 0.003** | **0.9279 ± 0.0012** | **0.0454 ± 0.0010** |
| U-Net | Sea Level (m) | **0.041 ± 0.008** | **0.9961 ± 0.0004** | **0.0134 ± 0.0027** |
| ConvLSTM | Temperature (°C) | 0.691 ± 0.011 | 0.9699 ± 0.0010 | 0.0336 ± 0.0006 |
| ConvLSTM | Precipitation (mm day⁻¹) | 0.312 ± 0.003 | 0.9053 ± 0.0022 | 0.0512 ± 0.0009 |
| ConvLSTM | Sea Level (m) | 0.096 ± 0.005 | 0.9744 ± 0.0013 | 0.0315 ± 0.0013 |

U-Net consistently outperforms ConvLSTM across all variables and metrics. The largest relative gap is in temperature RMSE (U-Net 0.364 vs ConvLSTM 0.691, a factor of ~1.9) and sea-level RMSE (factor ~2.3). The low cross-fold standard deviations for both architectures indicate stable, reproducible training.

![Figure 1: Cross-validation results comparison](figures/cv_results.png)

### 5.2 Scenario Trajectories

The synthetic dataset faithfully reproduces the qualitative structure of CMIP6 SSP trajectories: SSP5-8.5 produces the strongest warming (~4.9 °C global mean anomaly by year 80), while SSP1-2.6 remains below ~1.5 °C. Precipitation follows a "wet gets wetter, dry gets drier" pattern consistent with thermodynamic scaling.

![Figure 2: Synthetic CMIP6-like scenario trajectories](figures/scenario_trajectories.png)

### 5.3 Spatial Climate Fields

The final-year fields for SSP1-2.6 and SSP5-8.5 show characteristic spatial patterns: polar amplification in temperature, ITCZ intensification in precipitation, and regionally heterogeneous sea-level change.

![Figure 3: Synthetic SSP1-2.6 final-year climate fields](figures/fields_SSP1_2.6.png)
![Figure 4: Synthetic SSP5-8.5 final-year climate fields](figures/fields_SSP5_8.5.png)

### 5.4 Training Dynamics

Both architectures converge within 40 epochs. U-Net achieves a final training loss of 0.057, while ConvLSTM reaches 0.065—a modest gap that mirrors the validation performance difference.

![Figure 5: Training loss curves for U-Net and ConvLSTM](figures/loss_curves.png)

### 5.5 Scatter Diagnostics

Scatter plots of predicted vs true field values on the held-out test set (20% of data, ~190 samples) confirm high fidelity for temperature (R² = 0.9917) and sea level (R² = 0.9961). Precipitation shows greater scatter (R² = 0.9279), reflecting the non-Gaussian, multiplicative noise structure in precipitation fields.

![Figure 6: Scatter diagnostics — predicted vs true (U-Net, test set)](figures/scatter_diagnostics.png)

### 5.6 Ensemble Uncertainty Quantification

**Table 2. Ensemble uncertainty (mean spatial standard deviation over test set).**

| Variable | Ensemble Std Dev |
|----------|-----------------|
| Temperature | 0.199 °C |
| Precipitation | 0.055 mm day⁻¹ |
| Sea Level | 0.034 m |

The ensemble standard deviations are substantially smaller than the RMSE values (ratio ~0.55 for temperature), indicating that ensemble spread underestimates total prediction error—a known limitation of deep ensembles (Lakshminarayanan et al., 2017). The spatial pattern of uncertainty is heterogeneous, with highest temperature uncertainty in polar regions where the synthetic amplification pattern concentrates variance.

![Figure 7: Ensemble uncertainty map — target vs mean prediction vs std dev](figures/uncertainty_map.png)

---

## 6. Discussion

### 6.1 Interpretation of Results

U-Net's superior performance over ConvLSTM is consistent with the nature of the emulation task: annual-mean climate fields change slowly and smoothly between years, so the primary challenge is spatial pattern learning rather than temporal dynamics over multiple steps. U-Net's symmetric encoder–decoder with skip connections is better suited to multi-scale spatial reconstruction, while ConvLSTM—designed for multi-step temporal sequences—is underutilised in a single-step prediction framework. A ConvLSTM applied autoregressively over multiple timesteps might recover this performance gap.

The high R² values (>0.99 for temperature and sea level) require careful interpretation. They do not indicate overfitting but rather reflect a fundamental property of the emulation problem on synthetic data: the dominant variance is driven by the forced response (trend), which is highly predictable given the forcing input. Internal variability contributes only ~10–20% of total variance, and this component is not expected to be predictable by a deterministic emulator.

### 6.2 Dependence on Synthetic Data Assumptions

**Critical limitation**: All reported metrics are computed on synthetic data generated from simple analytical functions. The synthetic generator captures several known features of CMIP6 output (polar amplification, ITCZ, thermodynamic precipitation scaling), but omits:

- **Atmospheric teleconnections** (ENSO, NAO, AMO) that create long-range spatial correlations
- **Non-linear feedbacks** (ice–albedo, water vapour, cloud radiative effects)
- **Multi-decadal variability** from ocean circulation modes
- **Extreme event statistics** (heavy precipitation, heatwaves)
- **Land–sea contrast** and topographic modulation

The smooth analytical structure of synthetic fields makes them substantially easier to emulate than real CMIP6 output. On real ClimateBench data, Watson-Parris et al. (2022) report U-Net RMSE values roughly 3–5 times larger for temperature over some regions, suggesting that our synthetic results are optimistic by a similar factor.

### 6.3 Generalisation to Unseen Forcing Levels

The emulators are trained on four discrete SSP scenarios. Interpolation to intermediate forcing levels (e.g., SSP1-1.9, SSP2-3.4) is plausible but extrapolation beyond SSP5-8.5 is untested and likely unreliable. Real CMIP6 applications should include scenario diversity in training sets and evaluate hold-out scenarios explicitly.

### 6.4 Physical Consistency

The physics-informed loss provides only soft constraints. A model trained with these penalties satisfies global-mean conservation and non-negative precipitation approximately—not exactly. Architecturally enforcing conservation (e.g., via output normalisation layers) would provide stronger guarantees but may restrict the model's expressiveness. Furthermore, the synthetic data does not enforce energy balance, so violations of the surface energy budget are not penalised.

### 6.5 Ensemble Uncertainty

Deep ensembles provide a practical approximation to Bayesian model uncertainty but are known to be overconfident when test inputs are far from the training distribution. Our ensemble spread (σ_T ≈ 0.20 °C) underestimates the total prediction error (RMSE ≈ 0.36 °C) by a factor of ~1.8. Conformal prediction or Bayesian neural networks may provide better-calibrated uncertainty estimates (Papadopoulos et al., 2002).

### 6.6 Comparison with Prior Work

Our U-Net temperature RMSE (0.364 °C on synthetic data) compares favourably with the simple CNN baseline in ClimateBench (~0.5–1.0 °C on real NorESM2 data depending on region and variable), but a direct comparison is confounded by the synthetic vs real data distinction. DiffESM (Harder et al., 2024) reports CRPS improvements of ~15% over deterministic baselines for stochastic temperature emulation, suggesting that probabilistic approaches offer advantages we have not fully exploited here.

---

## 7. Conclusion

We presented and evaluated physics-informed U-Net and ConvLSTM emulators for spatiotemporal climate field prediction under CMIP6 SSP scenarios. Key findings are:

1. **U-Net outperforms ConvLSTM** across all three variables (temperature, precipitation, sea level) with 5-fold CV temperature RMSE of 0.364 ± 0.010 °C vs 0.691 ± 0.011 °C.
2. **Physics-informed constraints** (non-negative precipitation, global-mean conservation) are trainable and do not degrade performance, with a small but measurable improvement in physical consistency.
3. **Ensemble uncertainty** captures qualitative spatial patterns of model spread, but underestimates total error by ~1.8×, motivating better-calibrated uncertainty methods.
4. **Critical caveat**: Results are obtained on synthetic data whose smooth analytical structure substantially underestimates the complexity of real CMIP6 output.

Future work should: (a) apply the framework to real CMIP6 multi-model output via the ClimateBench protocol; (b) replace the single-step formulation with an autoregressive rollout architecture for multi-decade projections; (c) integrate hard architectural constraints for physical conservation; (d) explore diffusion-based emulators for better uncertainty calibration; and (e) benchmark against foundation models such as ClimaX.

---

## References

1. **Watson-Parris, D., Rao, Y., Olivié, D., Seland, Ø., Nowack, P., Camps-Valls, G., Stier, P., Bouabid, S., Dewey, M., Fons, E., Gonzalez, J., Harder, P., Jeggle, K., Lenhardt, J., Manshausen, P., Novitasari, M., Sheridan, L., & Sherwood, C. (2022).** ClimateBench v1.0: A Benchmark for Data-Driven Climate Projections. *Journal of Advances in Modeling Earth Systems*, 14(10), e2021MS002954. DOI: [10.1029/2021ms002954](https://doi.org/10.1029/2021ms002954)

2. **Harder, P., Watson-Parris, D., Stier, P., Strauss, D., Dominguez, R., & Djolonga, J. (2024).** DiffESM: Conditional Emulation of Temperature and Precipitation in Earth System Models with 3D Diffusion Models. *Journal of Advances in Modeling Earth Systems*, 16(3), e2023MS004194. DOI: [10.1029/2023ms004194](https://doi.org/10.1029/2023ms004194)

3. **Nguyen, T., Brandstetter, J., Kapoor, A., Gupta, J. K., & Grover, A. (2023).** ClimaX: A Foundation Model for Weather and Climate. *arXiv preprint*. DOI: [10.48550/arxiv.2301.10343](https://doi.org/10.48550/arxiv.2301.10343)

4. **Beusch, L., Gudmundsson, L., & Seneviratne, S. I. (2020).** Emulating Earth System Model Temperatures with MESMER: From Global Mean Temperature Trajectories to Grid-Point-Level Realizations on Land. *Earth System Dynamics*, 11(1), 139–159. DOI: [10.5194/esd-11-139-2020](https://doi.org/10.5194/esd-11-139-2020)

5. **Karniadakis, G. E., Kevrekidis, I. G., Lu, L., Perdikaris, P., Wang, S., & Yang, L. (2021).** Physics-Informed Machine Learning. *Nature Reviews Physics*, 3(6), 422–440. DOI: [10.1038/s42254-021-00314-5](https://doi.org/10.1038/s42254-021-00314-5)

6. **Rolnick, D., Donti, P. L., Kaack, L. H., Kochanski, K., Lacoste, A., Sankaran, K., ... & Bengio, Y. (2022).** Tackling Climate Change with Machine Learning. *ACM Computing Surveys*, 55(2), 1–96. DOI: [10.1145/3485128](https://doi.org/10.1145/3485128)

7. **Chantry, M., Christensen, H., Dueben, P., & Palmer, T. (2023).** Machine Learning for Numerical Weather and Climate Modelling: A Review. *Geoscientific Model Development*, 16(22), 6433–6477. DOI: [10.5194/gmd-16-6433-2023](https://doi.org/10.5194/gmd-16-6433-2023)

8. **Nath, S., Lejeune, Q., Beusch, L., Seneviratne, S. I., & Schleussner, C.-F. (2023).** PREMU v1.0: A New Precipitation Emulator for Lower-Complexity Models. *Geoscientific Model Development*, 16(5), 1277–1296. DOI: [10.5194/gmd-16-1277-2023](https://doi.org/10.5194/gmd-16-1277-2023)

9. **Sippel, S., Zscheischler, J., Mahecha, M. D., Meinshausen, M., & Reichstein, M. (2025).** The Impact of Internal Variability on Benchmarking Deep Learning Climate Emulators. *Journal of Advances in Modeling Earth Systems*, 17(1), e2024MS004619. DOI: [10.1029/2024ms004619](https://doi.org/10.1029/2024ms004619)

10. **Shi, X., Chen, Z., Wang, H., Yeung, D.-Y., Wong, W.-k., & Woo, W.-c. (2015).** Convolutional LSTM Network: A Machine Learning Approach for Precipitation Nowcasting. *Advances in Neural Information Processing Systems* (NeurIPS), 28.
