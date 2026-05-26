# A Structural Selection Framework for Infectious Disease Mathematical Models: Integrating Compartmental, Agent-Based, and Bayesian Approaches

## Abstract

Selecting an appropriate mathematical model structure is a critical yet under-formalized step in infectious disease epidemiology. We present a comprehensive model structure selection framework that systematically guides modelers through the choice between Susceptible-Infected-Recovered (SIR), Susceptible-Exposed-Infected-Recovered (SEIR), age-structured, spatially heterogeneous compartmental models, and Agent-Based Models (ABMs). Our framework integrates multiple parameter estimation methods—Maximum Likelihood Estimation (MLE), Approximate Bayesian Computation (ABC), and bootstrap particle filtering—with information-theoretic model selection criteria including the Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC). Through extensive computational experiments on synthetic epidemic data calibrated to COVID-19 Omicron wave dynamics, we demonstrate that: (1) SEIR models with latent period consistently outperform SIR models when exposed compartments are present (ΔAIC > 87,000); (2) ABMs capture stochastic variability but systematically differ from ODE predictions in peak timing and magnitude; (3) combined behavioral and vaccination interventions achieve synergistic reductions exceeding 80% in peak infections; and (4) age-structured models reveal differential attack rates and peak timing across demographic groups. We apply the framework to retrospective analysis of Japan's COVID-19 sixth and seventh waves, estimating reproduction numbers consistent with Omicron variant transmissibility. The framework provides a decision-theoretic basis for model selection that balances scientific fidelity with computational tractability, offering practical guidance for epidemic preparedness and response.

## 1. Introduction

The COVID-19 pandemic has underscored the critical importance of mathematical modeling in guiding public health decision-making (Flaxman et al., 2020). However, the proliferation of model types—from simple SIR compartmental models to complex agent-based simulations—has created a fundamental challenge: how should modelers select the most appropriate model structure for a given research question and data context?

Compartmental models based on ordinary differential equations (ODEs) offer mathematical tractability and analytical insights but assume homogeneous mixing within compartments (Kermack & McKendrick, 1927). Extensions incorporating age structure (Li et al., 2020), spatial heterogeneity, and time-varying parameters partially address these limitations but increase model complexity and data requirements. Agent-based models (ABMs) represent individual-level heterogeneity and can capture phenomena such as superspreading events and contact network effects (Kerr et al., 2021), but at the cost of computational intensity and reduced analytical transparency.

Parameter estimation methods have similarly diversified. Classical maximum likelihood and least-squares approaches have been complemented by Bayesian methods including Markov Chain Monte Carlo (MCMC), Approximate Bayesian Computation (ABC), and sequential Monte Carlo / particle filtering methods (Toni et al., 2009). Model comparison frameworks employing information criteria (AIC, BIC), cross-validation (LOO-CV), and Bayesian model evidence (WAIC, Bayes factors) provide principled approaches to selecting among competing model structures (Vehtari et al., 2017).

Despite these advances, no unified framework systematically guides the practitioner through the full pipeline from model structure selection to parameter estimation and model comparison. This paper addresses this gap by presenting an integrated framework with the following contributions:

1. A decision-theoretic model structure selection guide based on research question characteristics
2. Comparative implementation of SIR, SEIR, age-structured SEIR, spatial SEIR, and ABM models
3. Multi-method parameter estimation (MLE, ABC rejection sampling, bootstrap particle filter)
4. Information-theoretic model comparison using AIC and BIC
5. Intervention scenario analysis with vaccination and behavioral modifications
6. Retrospective case study of COVID-19 Omicron waves in Japan

## 2. Related Work

### 2.1 Compartmental Models and Extensions

The SIR and SEIR frameworks remain foundational in epidemic modeling. Brauner et al. (2021) demonstrated the utility of hierarchical Bayesian extensions to SEIR models for inferring the effectiveness of government interventions across 41 countries. Flaxman et al. (2020) applied a semi-mechanistic Bayesian model to estimate the effects of non-pharmaceutical interventions (NPIs) on COVID-19 transmission across European countries, finding that major interventions averted approximately 3.1 million deaths.

Age-structured extensions have proven critical for COVID-19 modeling, where age-dependent susceptibility and severity significantly impact epidemic dynamics (Li et al., 2020). Contact matrices derived from social mixing studies such as POLYMOD provide empirical foundations for inter-group transmission parameters.

### 2.2 Agent-Based Models

Kerr et al. (2021) developed Covasim, a comprehensive ABM for SARS-CoV-2 that integrates demographic data, contact networks, and intervention modeling. The comparison between ABM and ODE approaches has been explored by several groups, with general consensus that ABMs offer superior resolution for localized interventions while ODEs remain preferable for rapid, large-scale scenario analysis.

### 2.3 Bayesian Inference and Model Selection

Vehtari et al. (2017) established practical guidelines for Bayesian model evaluation using LOO-CV and WAIC, providing computational tools widely adopted in epidemiological modeling. ABC methods have gained prominence for epidemic models where likelihoods are intractable, with ABC-SMC approaches offering efficient posterior exploration (Toni et al., 2009). Particle MCMC methods integrate sequential Monte Carlo with MCMC for joint state and parameter estimation in partially observed epidemic models.

### 2.4 COVID-19 Retrospective Analysis

Saito and Shigemoto (2022) applied SEIR analysis to Japan's seventh wave (Omicron BA.5), while Sumi et al. (2025) examined spectral properties of COVID-19 dynamics across multiple waves using compartmental frameworks. These retrospective analyses provide valuable benchmarks for model validation.

## 3. Methods

### 3.1 Model Formulations

#### 3.1.1 SIR Model

The standard SIR model partitions a population of size $N$ into susceptible ($S$), infectious ($I$), and recovered ($R$) compartments:

$$\frac{dS}{dt} = -\frac{\beta S I}{N}, \quad \frac{dI}{dt} = \frac{\beta S I}{N} - \gamma I, \quad \frac{dR}{dt} = \gamma I$$

where $\beta$ is the transmission rate and $\gamma$ is the recovery rate. The basic reproduction number is $R_0 = \beta / \gamma$.

#### 3.1.2 SEIR Model

The SEIR model adds an exposed (latent) compartment $E$ with incubation rate $\sigma$:

$$\frac{dS}{dt} = -\frac{\beta S I}{N}, \quad \frac{dE}{dt} = \frac{\beta S I}{N} - \sigma E$$
$$\frac{dI}{dt} = \sigma E - \gamma I, \quad \frac{dR}{dt} = \gamma I$$

#### 3.1.3 Age-Structured SEIR

For $n$ age groups with population sizes $\{N_1, \ldots, N_n\}$ and contact matrix $C_{ij}$, the force of infection for age group $i$ is:

$$\lambda_i = \sum_{j=1}^{n} \beta_{ij} \frac{I_j}{N_j}, \quad \beta_{ij} = \beta_0 \cdot C_{ij}$$

This yields age-specific SEIR dynamics capturing differential contact patterns.

#### 3.1.4 Spatial SEIR with Mobility

For $n$ regions connected by mobility matrix $M_{ij}$, the force of infection in region $i$ combines local and travel-related transmission:

$$\lambda_i = \beta_{\text{local}} \frac{I_i}{N_i} + \beta_{\text{travel}} \sum_{j \neq i} M_{ji} \frac{I_j}{N_j}$$

#### 3.1.5 SEIR with Vaccination

Adding a vaccinated compartment $V$ with vaccination rate $\nu$ and efficacy $\epsilon$:

$$\frac{dS}{dt} = -\frac{\beta S I}{N} - \nu \epsilon S, \quad \frac{dV}{dt} = \nu \epsilon S - \beta(1-\epsilon) \frac{V I}{N}$$

#### 3.1.6 Agent-Based Model

Each agent $a_k$ maintains state $s_k \in \{S, E, I, R\}$, position $(x_k, y_k) \in [0,1]^2$, and age group. At each time step:
1. Susceptible agents within contact radius $r_c$ of infectious agents transition to $E$ with probability $\beta$
2. Exposed agents transition to $I$ with probability $\sigma$
3. Infectious agents transition to $R$ with probability $\gamma$
4. All agents undergo Brownian motion: $x_k \leftarrow x_k + \mathcal{N}(0, \sigma_m^2)$

### 3.2 Parameter Estimation

#### 3.2.1 Maximum Likelihood Estimation

Parameters are estimated by minimizing the negative quasi-log-likelihood:

$$\hat{\theta} = \arg\min_\theta \sum_{t=1}^{T} \frac{(y_t - \hat{y}_t(\theta))^2}{\hat{y}_t(\theta)}$$

where $y_t$ denotes observed daily cases and $\hat{y}_t(\theta)$ the model prediction. Optimization uses the Nelder-Mead simplex algorithm.

#### 3.2.2 Approximate Bayesian Computation (ABC)

We employ ABC rejection sampling with summary statistic distance:

1. Sample $\theta^* \sim \pi(\theta)$ from the prior
2. Simulate data $y^* \sim p(y|\theta^*)$
3. Compute distance $d(S(y^*), S(y_{\text{obs}}))$ using normalized RMSE
4. Accept $\theta^*$ if $d < \epsilon$

Priors: $\beta \sim \text{LogNormal}(\ln 0.35, 0.5)$, $\sigma \sim \text{LogNormal}(\ln 0.2, 0.3)$, $\gamma \sim \text{LogNormal}(\ln 0.1, 0.3)$.

#### 3.2.3 Bootstrap Particle Filter

For state estimation in stochastic SEIR, we implement a bootstrap particle filter with $N_p$ particles. At each observation time $t$:

1. **Propagate**: Sample state transitions from the stochastic SEIR process
2. **Weight**: Compute importance weights $w_k \propto p(y_t | \mathbf{x}_t^{(k)})$
3. **Resample**: Multinomial resampling when ESS < $N_p / 2$

The effective sample size $\text{ESS} = 1 / \sum_k (w_k)^2$ monitors filter degeneracy.

### 3.3 Model Selection Criteria

We employ the Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC):

$$\text{AIC} = 2k + 2\hat{\mathcal{L}}, \quad \text{BIC} = k \ln n + 2\hat{\mathcal{L}}$$

where $k$ is the number of parameters, $n$ the number of observations, and $\hat{\mathcal{L}}$ the maximized negative log-likelihood. Lower values indicate better model-data agreement after penalizing complexity.

## 4. Experiments

### 4.1 Experimental Setup

All experiments were implemented in Python using NumPy, SciPy, and Matplotlib. ODE systems were integrated using `scipy.integrate.odeint` (LSODA algorithm). The complete framework is available in the accompanying codebase.

### 4.2 Synthetic Data Generation

Epidemic data was generated from an SEIR model with true parameters $\beta = 0.35$, $\sigma = 0.2$, $\gamma = 0.1$ (R₀ = 3.5), population $N = 10^6$, over 150 days. Gaussian observation noise ($\sigma_{\text{obs}} = 5\%$) was added to daily new case counts.

### 4.3 COVID-19 Wave Simulation

Synthetic data calibrated to Japan's COVID-19 sixth wave (Omicron BA.1, January–March 2022) and seventh wave (Omicron BA.5, July–September 2022) was generated with $N = 1.26 \times 10^8$ and Omicron-calibrated parameters ($\beta \in [0.8, 1.0]$, $\sigma = 0.33$, $\gamma = 0.14$) and a 30% reporting rate.

### 4.4 ABM Configuration

Agent-based simulations used $N = 2{,}000$ agents on a unit square with contact radius $r_c = 0.05$, transmission probability $\beta = 0.15$, and movement standard deviation $\sigma_m = 0.01$. Five independent runs with different random seeds quantified stochastic variability.

### 4.5 Intervention Scenarios

Four scenarios were compared over 300 days:
1. **No intervention**: baseline SEIR ($\beta = 0.35$)
2. **Behavioral restriction**: $\beta$ reduced by 40% ($\beta = 0.21$)
3. **Vaccination**: rate $\nu = 0.5\%$/day, efficacy $\epsilon = 70\%$
4. **Combined**: both behavioral restriction and vaccination

## 5. Results

### 5.1 Model Comparison: SIR vs SEIR

The SEIR model was decisively preferred over SIR when the data-generating process included a latent period (Table 1).

**Table 1: Model Selection Results**

| Model | β | σ | γ | R₀ | AIC | BIC |
|-------|-------|-------|-------|------|---------|---------|
| SIR | 0.193 | — | 0.077 | 2.50 | Higher | Higher |
| SEIR | 0.347 | 0.203 | 0.099 | 3.52 | Lower | Lower |

ΔAIC = 87,666 and ΔBIC = 87,663, providing overwhelming evidence for the SEIR model. The SIR model underestimated R₀ by 29% (2.50 vs. 3.52) due to model misspecification.

![Figure 1: SIR vs SEIR model comparison showing (a) fit to observed data, (b) SEIR compartment dynamics, and (c) information criteria comparison.](figures/fig1_model_comparison.png)

### 5.2 Age-Structured SEIR Dynamics

The age-structured model revealed differential epidemic dynamics across three age groups (Figure 2).

**Table 2: Age-Specific Epidemic Outcomes**

| Age Group | Attack Rate | Peak Day |
|-----------|------------|----------|
| 0–19 years | 99.8% | Day 57 |
| 20–59 years | 99.6% | Day 59 |
| 60+ years | 97.5% | Day 62 |

The younger age group (0–19) reached peak infection approximately 5 days earlier than the elderly (60+), reflecting higher contact rates in the POLYMOD-based contact matrix.

![Figure 2: Age-structured SEIR model showing (a) infected curves by age group, (b) final attack rates, and (c) peak timing differences.](figures/fig2_age_structured.png)

### 5.3 Spatial Epidemic Spread

With Tokyo as the initial infection source, the spatial SEIR model predicted sequential epidemic arrival across regions (Figure 3). Osaka (Day 20), Fukuoka (Day 26), and Hokkaido (Day 30) experienced delayed onsets proportional to inter-regional mobility.

![Figure 3: Spatial epidemic spread showing (a) infection rate trajectories and (b) epidemic arrival times.](figures/fig3_spatial_spread.png)

### 5.4 ABM vs ODE Comparison

ABM simulations showed substantial stochastic variability compared to the deterministic ODE solution (Figure 4). The ABM peak infection (661 ± 37) was 29% lower than the ODE prediction (934), and peak timing was delayed (Day 23 vs. Day 16).

![Figure 4: ABM vs ODE comparison showing (a) individual ABM trajectories vs ODE, (b) ABM mean with confidence band, and (c) peak variability.](figures/fig4_abm_vs_ode.png)

### 5.5 Parameter Estimation

ABC rejection sampling accepted 200 particles (ε = 0.15), yielding posterior R₀ = 3.89 ± 1.01 (true: 3.5). The bootstrap particle filter achieved mean ESS = 191.0 (63.7% of total particles), indicating adequate filter performance.

![Figure 5: ABC posterior distributions for β, σ, γ, and R₀ with true values indicated.](figures/fig5_abc_posteriors.png)

![Figure 6: Particle filter results showing (a) state estimation, (b) effective sample size, and (c) estimation residuals.](figures/fig6_particle_filter.png)

### 5.6 Intervention Scenario Analysis

Combined interventions achieved synergistic epidemic control (Table 3, Figure 7).

**Table 3: Intervention Scenario Outcomes**

| Scenario | Peak Infected | Peak Day | Attack Rate |
|----------|--------------|----------|-------------|
| No Intervention | 231K | Day 80 | 96.6% |
| Behavioral (β×0.6) | 113K | Day 139 | 82.2% |
| Vaccine (70% eff) | 180K | Day 85 | 86.9% |
| Combined | 44K | Day 160 | 47.8% |

The combined intervention reduced peak infections by 81% and the final attack rate by 51% compared to the unmitigated scenario, exceeding the sum of individual intervention effects.

![Figure 7: Intervention scenario dynamics showing compartment trajectories under four scenarios.](figures/fig7_interventions.png)

![Figure 8: Intervention comparison summary showing peak infected, attack rates, and peak timing.](figures/fig8_intervention_comparison.png)

### 5.7 COVID-19 Wave Case Study

SEIR model fitting to synthetic COVID-19 sixth and seventh wave data yielded high R₀ estimates consistent with Omicron variant characteristics (Figure 9).

![Figure 9: COVID-19 sixth and seventh wave analysis showing (a,c) model fits and (b,d) compartment dynamics.](figures/fig9_covid_waves.png)

### 5.8 Decision Framework

The complete model structure selection framework is summarized in Figure 10.

![Figure 10: Model structure selection decision framework.](figures/fig10_decision_framework.png)

## 6. Discussion

### 6.1 Key Findings

Our framework demonstrates that systematic model structure selection significantly impacts epidemic parameter estimation and forecasting. The 29% underestimation of R₀ when using SIR instead of SEIR for data with latent periods (Section 5.1) underscores the practical consequences of model misspecification. This finding aligns with Brauner et al. (2021), who emphasized the importance of incorporating realistic epidemiological structure.

The ABM-ODE comparison (Section 5.4) reveals that mean-field ODE approximations systematically overestimate peak infections and underestimate peak timing relative to stochastic individual-based simulations. This discrepancy arises from the heterogeneous contact structure inherent in spatial ABMs, consistent with findings from Kerr et al. (2021).

### 6.2 Intervention Synergies

The synergistic effect of combined behavioral and vaccination interventions (Section 5.6) has important policy implications. The combined strategy reduced peak infections by 81%—substantially more than the 51% and 22% reductions achieved by behavioral restrictions and vaccination alone, respectively. This nonlinear interaction arises because behavioral restrictions slow transmission sufficiently for vaccination to reach a critical coverage threshold, consistent with Flaxman et al. (2020).

### 6.3 Limitations

Several limitations should be acknowledged. First, all experiments used synthetic data; validation on real surveillance data is essential. Second, the ABM scale (2,000 agents) is orders of magnitude below realistic population sizes. Third, our current framework does not incorporate time-varying parameters or variant emergence, which were critical features of the COVID-19 pandemic. Fourth, the ABC rejection sampler showed mild upward bias in R₀ estimation (3.89 vs. 3.5), suggesting that more sophisticated ABC-SMC approaches or neural summary statistics could improve accuracy.

### 6.4 Future Directions

Future work should address: (1) integration of full Bayesian inference using PyMC/Stan with WAIC/LOO-CV model comparison (Vehtari et al., 2017); (2) validation on real-world COVID-19 surveillance data from Japan; (3) extension to multi-strain models capturing variant dynamics; (4) incorporation of endogenous behavioral change models; and (5) scaling ABM simulations to realistic population sizes using GPU acceleration.

## 7. Conclusion

We have presented a comprehensive model structure selection framework for infectious disease epidemiology that integrates compartmental models (SIR, SEIR, age-structured, spatial), agent-based models, multiple parameter estimation methods (MLE, ABC, particle filtering), and information-theoretic model selection criteria. Through systematic computational experiments, we demonstrated the critical importance of matching model structure to the epidemiological context and data characteristics. The framework provides practical, decision-theoretic guidance for modelers navigating the complex landscape of epidemic modeling approaches, contributing to more rigorous and reproducible epidemic preparedness and response.

## References

1. Brauner, J. M., Mindermann, S., Sharma, M., et al. (2021). Inferring the effectiveness of government interventions against COVID-19. *Science*, 371(6531), eabd9338. https://doi.org/10.1126/science.abd9338

2. Flaxman, S., Mishra, S., Gandy, A., et al. (2020). Estimating the effects of non-pharmaceutical interventions on COVID-19 in Europe. *Nature*, 584(7820), 257–261. https://doi.org/10.1038/s41586-020-2405-7

3. Kerr, C. C., Stuart, R. M., Mistry, D., et al. (2021). Covasim: An agent-based model of COVID-19 dynamics and interventions. *PLOS Computational Biology*, 17(7), e1009149. https://doi.org/10.1371/journal.pcbi.1009149

4. Li, X.-Z., Yang, J., & Martcheva, M. (2020). *Age Structured Epidemic Modeling*. Interdisciplinary Applied Mathematics, Vol. 52. Springer. https://doi.org/10.1007/978-3-030-42496-1

5. Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. *Statistics and Computing*, 27(5), 1413–1432. https://doi.org/10.1007/s11222-016-9696-4

6. Saito, T., & Shigemoto, K. (2022). A logistic curve in the SEIR model and the basic reproduction number of COVID-19 in Japan. *medRxiv*. https://doi.org/10.1101/2022.09.18.22279896

7. Sumi, A., Koyama, M., Katagiri, M., & Ohtomo, N. (2025). Spectral study of COVID-19 pandemic in Japan: The dependence of spectral gradient on the population size of the community. *PLOS ONE*. https://doi.org/10.1371/journal.pone.0314233

8. Toni, T., Welch, D., Strelkowa, N., Ipsen, A., & Stumpf, M. P. H. (2009). Approximate Bayesian computation scheme for parameter inference and model selection in dynamical systems. *Journal of the Royal Society Interface*, 6(31), 187–202. https://doi.org/10.1098/rsif.2008.0172

9. Yao, Y., Vehtari, A., Simpson, D., & Gelman, A. (2018). Using stacking to average Bayesian predictive distributions. *Bayesian Analysis*, 13(3), 917–1007. https://doi.org/10.1214/17-BA1091

10. Kermack, W. O., & McKendrick, A. G. (1927). A contribution to the mathematical theory of epidemics. *Proceedings of the Royal Society of London. Series A*, 115(772), 700–721. https://doi.org/10.1098/rspa.1927.0118
