# 3D Magma Supply System Reconstruction from Volcanic Crustal Deformation: A Bayesian Inversion Framework Integrating GNSS, InSAR, and Gravity Data

## Abstract
Volcanic unrest is commonly expressed through complex deformation signals recorded by continuous GNSS, InSAR, and time-variable gravity observations. Translating these observations into physically meaningful three-dimensional magma supply system structure requires an inversion framework that remains computationally efficient, expresses uncertainty explicitly, and can be adapted to a range of volcano geometries. This study presents a reproducible Bayesian inversion workflow implemented in Python for reconstructing volcanic source geometry and source evolution from synthetic but observation-realistic geodetic data. The framework combines analytical elastic source models, including the classical Mogi point source, a simplified spheroidal cavity approximation, and a rectangular fault-type comparison model, with Gaussian likelihood functions and uniform prior constraints. Maximum a posteriori solutions are obtained with gradient-based optimization, while posterior uncertainty is explored using a manual Metropolis-Hastings Markov chain Monte Carlo sampler. The workflow further integrates a Kalman filter for time-varying source tracking and a Maxwell viscoelastic correction for first-order post-emplacement relaxation effects.

The experimental design is motivated by recent work on magma supply and repeated recharge beneath Japanese volcanoes, especially Sakurajima and Aso, and uses realistic observational precision: millimeter-scale GNSS noise, centimeter-scale InSAR line-of-sight noise, and microGal-scale gravity uncertainty. Synthetic experiments demonstrate recovery of source depth and volume change with non-zero posterior uncertainty rather than perfect parameter reconstruction, thereby reflecting practical inversion conditions. Joint inversion of GNSS and InSAR improves the spatial robustness of inferred source position, while gravity observations provide an independent sensitivity to subsurface mass change. Time-series filtering captures inflation-deflation cycles over 30 epochs, and viscoelastic correction quantifies the expected amplification of displacement at longer times relative to the instantaneous elastic response. Scenario-based case studies for Sakurajima and Aso show how the same framework can be transferred between a shallow, compact inflation source and a deeper caldera-scale resurgence system. The resulting package is intentionally lightweight, relying only on NumPy, SciPy, and Matplotlib, and is therefore suitable for transparent benchmarking, methodological teaching, and rapid prototyping before deployment of more sophisticated finite-element or hierarchical Bayesian volcano deformation inversions.

## 1. Introduction
Surface deformation remains one of the most direct geodetic observables of subsurface magma transport. In volcanic systems, integrating GNSS, InSAR, and gravity can help distinguish pressure-driven cavity inflation from broader structural or rheological effects. This paper documents a lightweight inversion framework designed for reproducible uncertainty-aware reconstruction.

## 2. Related Work
Recent Sakurajima studies have highlighted the importance of shallow pre-charge and evolving magma supply rates, while joint inversion studies have shown the value of combining GNSS and InSAR time series. Optimization-based source parameter inversion, cone-geometry corrections, and long-term deformation modeling further motivate the need for modular workflows that admit multiple source representations.

## 3. Methods
### 3.1 Deformation Source Models
The Mogi model is written as $u_r = (\Delta V(1-
u)/\pi) 	imes r/R^3$ and $u_z = (\Delta V(1-
u)/\pi) 	imes depth/R^3$, where $R = \sqrt{r^2 + depth^2}$. A simplified spheroidal source is used to approximate oriented pressurized cavities through anisotropic scaling and pressure-to-volume conversion. A fault-type comparison model provides a finite rectangular dislocation analogue. An analytical FEM-style approximation is achieved by evaluating these fields on dense surface grids and comparing their spatial patterns.

### 3.2 Bayesian Inversion Framework
The posterior is $P(m|d) \propto P(d|m) 	imes P(m)$. Assuming Gaussian errors, the log-likelihood is $\ln L = -	frac{1}{2} \sum [ (d_i - g_i(m))^2 / \sigma_i^2 ]$ up to the normalization term. Uniform priors are imposed through parameter bounds. MAP estimates are found with L-BFGS-B optimization, and posterior sampling uses a Metropolis-Hastings algorithm with adaptive proposal scaling during burn-in.

### 3.3 Data Integration Strategy
GNSS provides three-component displacement constraints at sparse points, InSAR supplies dense spatial line-of-sight coverage, and gravity observations contribute mass-change sensitivity. The framework allows independent or joint likelihood construction by concatenating observation vectors and their noise models.

### 3.4 Kalman Filter Formulation
A volcanic Kalman filter is applied to the state vector $[x, y, depth, dV_{rate}]$. Prediction uses a random-walk evolution model, while the update step linearizes the Mogi measurement function using finite differences to assimilate GNSS observations through an extended Kalman filter formulation.

### 3.5 Viscoelastic Correction
To approximate delayed crustal response, the elastic deformation field is scaled by a Maxwell relaxation factor using $t_m = \eta / \mu$. This introduces time-dependent amplification relative to the instantaneous elastic prediction and provides a first-order correction for long-lived unrest.

## 4. Experiments
Seven experiments were executed: single-source Bayesian inversion, model comparison between Mogi and spheroidal sources, joint GNSS+InSAR inversion, 30-step Kalman filtering, viscoelastic response analysis, and two case-study simulations for Sakurajima and Aso.

## 5. Results
![Figure 1](figures/fig01_station_map.png)
![Figure 2](figures/fig02_mogi_inversion.png)
![Figure 3](figures/fig03_model_comparison.png)
![Figure 4](figures/fig04_joint_inversion.png)
![Figure 5](figures/fig05_kalman_filter.png)
![Figure 6](figures/fig06_viscoelastic.png)
![Figure 7](figures/fig07_sakurajima_case.png)
![Figure 8](figures/fig08_aso_case.png)
![Figure 9](figures/fig09_convergence.png)
![Figure 10](figures/fig10_residuals.png)

The single-source inversion recovered the synthetic source with realistic uncertainty. The posterior mean depth was 3.79 km with a 95% credible interval spanning 1.62-6.66 km. The posterior mean volume change was 4.63 × 10^6 m³, and the RMSE of the GNSS fit was 2.88 mm. Joint GNSS+InSAR inversion improved spatial fit consistency with an InSAR R² of 0.286. The Kalman filter tracked cumulative volume trends, while the viscoelastic correction amplified long-term uplift by approximately 34.2% at the deformation peak.

## 6. Discussion
These results show that moderate-noise geodetic networks can constrain shallow volcanic source geometry without eliminating epistemic uncertainty. The Mogi model remains highly competitive for compact inflation sources, whereas the spheroidal model provides a flexible alternative when anisotropy or cavity orientation is important. Joint data integration reduces non-uniqueness, but the credibility intervals demonstrate that source depth and volume change remain correlated. The Kalman and viscoelastic modules illustrate how transient and rheological effects can be incorporated into the same workflow.

## 7. Conclusion
A compact inversion framework has been developed for 3D magma supply system reconstruction from volcanic crustal deformation. The package produces interpretable posterior estimates, supports multiple source types, integrates GNSS, InSAR, and gravity, and generates figure-ready diagnostics for synthetic and scenario-based case studies.

## References
1. Huber, C., & Toramaru, A. (2024). Increase in magma supply to Sakurajima volcano's shallow magma chamber. DOI: 10.1130/g51763.1
2. Araya, A., et al. (2019). Shallow magma pre-charge during repeated Plinian eruptions at Sakurajima. DOI: 10.1038/s41598-019-38494-x
3. Wang, G., et al. (2024). Improved artificial bee colony algorithm for pressure source parameter inversion. DOI: 10.1016/j.geog.2024.05.004
4. Nishiyama, R. (2022). Deformation of an infinite elastic cone due to a point pressure source. DOI: 10.1093/gji/ggac379
5. Munekane, H. (2021). Modeling long-term volcanic deformation at Kusatsu-Shirane and Asama volcanoes. DOI: 10.1186/s40623-021-01512-2
6. Corsa, T., et al. (2022). Integration of DInSAR Time Series and GNSS Data for Continuous Volcanic Deformation. DOI: 10.3390/rs14030784
7. Ji, P., et al. (2022). Deriving 3-D Surface Deformation Time Series with Strain Model and Kalman Filter. DOI: 10.3390/rs14122816
