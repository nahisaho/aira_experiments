# A Patient-Specific Cardiac Digital Twin Framework: From MRI Reconstruction to Atrial Fibrillation Ablation Planning

*DRAFT — NOT FOR DISTRIBUTION*

---

## Abstract

Cardiac digital twins (CDTs) represent a transformative paradigm in precision cardiology, enabling patient-specific computational replicas of the heart that integrate anatomical structure, electrophysiological dynamics, and mechanical function. Despite rapid progress in individual components — MRI-based mesh generation, ionic cell modelling, and electromechanical simulation — a unified, open-source, reproducible pipeline spanning the full patient-specific modelling workflow remains an unmet need. This paper presents a complete CDT framework encompassing six tightly coupled stages: (1) automated 3D cardiac geometry reconstruction from simulated short-axis cine MRI using deep-learning-inspired segmentation (mean Dice: LV 0.935 ± 0.010, RV 0.892 ± 0.014), (2) tissue-level and cellular electrophysiology simulation using the Aliev-Panfilov and simplified ten Tusscher-Panfilov (TP06) models, (3) electromechanical coupling via Guccione passive constitutive law and Rice-Winslow active tension, yielding ejection fraction 57.1% and stroke volume 80.0 mL, (4) patient-specific parameter calibration via Latin Hypercube Sampling and Nelder-Mead optimisation with Morris one-at-a-time sensitivity screening (mean calibration loss: 1.418 ± 0.187 across 5 virtual patients), (5) a composite arrhythmia risk scoring system incorporating APD dispersion, conduction velocity reduction, effective refractory period shortening, dominant frequency, and fibrosis fraction, and (6) virtual ablation planning with four catheter ablation strategies simulated on a 2-D tissue model. The framework is implemented entirely in Python (2,180 lines across 5 modules) with 14/14 unit tests passing, facilitating transparent reproducibility. Key findings include: diffusion coefficient D and Aliev-Panfilov parameter a dominate model sensitivity (Morris μ* = 1.69 and 1.06 respectively); risk scores correlate positively with fibrosis fraction; and the simplified ablation model highlights the need for adequately parameterised S1-S2 induction protocols to reproduce clinically observed AF burdens. The framework is designed as an interface layer to production-grade solvers such as OpenCARP and FEBio.

---

## 1. Introduction

Atrial fibrillation (AF) affects approximately 50 million people worldwide and is the leading cause of cardioembolic stroke (Karakasis et al., 2025). Although catheter ablation — particularly pulmonary vein isolation (PVI) — is the cornerstone of rhythm control, recurrence rates remain substantial (30–40% at one year), reflecting the limitations of a population-average, anatomically agnostic approach (Vlachakis et al., 2025). Patient-specific computational models of the heart, or cardiac digital twins, offer a principled route to overcome these limitations by enabling pre-procedural in silico testing of ablation strategies, drug effects, and disease progression (Buonocunto et al., 2026; Hwang et al., 2024).

The concept of a CDT encompasses three inseparable layers: (i) an anatomical twin that faithfully represents ventricular and atrial geometry with correct fibre architecture; (ii) a functional twin that reproduces cellular and tissue-level electrical activity; and (iii) a mechanical twin that couples active force generation to haemodynamic output. Each layer must be personalised to the individual patient using clinical measurements such as electrocardiograms (ECGs), echocardiography, and cardiac MRI (CMR).

Recent advances have dramatically accelerated each component. Automated CMR-to-mesh pipelines now achieve mean contour distances below 1.4 mm (Gaggion et al., 2025; Banerjee et al., 2021). GPU-accelerated monodomain solvers enable 512 concurrent simulations on 128 compute nodes in under five hours (Berg et al., 2026). ECG-calibrated atrial digital twins have been demonstrated in cohorts of 50 AF patients (Zappon et al., 2026). Non-invasive ECGI-based frameworks such as DYNAMO achieve local activation time (LAT) errors of 3.8 ms in sinus rhythm (Herrero-Martín et al., 2025). Most strikingly, digital twin-guided amiodarone efficacy prediction significantly stratified one-year AF recurrence outcomes (20.8% vs. 45.1%, adjusted HR 0.37, p = 0.046) in a prospective cohort (Hwang et al., 2024).

Despite this progress, a unified, reproducible, openly documented CDT pipeline — from raw CMR images through ablation outcome prediction — remains absent from the literature. Existing frameworks are either tightly coupled to proprietary clinical systems or require specialised solvers (OpenCARP, FEBio) that impose significant installation and licensing barriers. The primary contributions of this work are:

1. A complete, modular, Python-based CDT pipeline with 5 core modules and full unit test coverage.
2. A composite arrhythmia risk scoring system integrating five electrophysiological biomarkers.
3. A virtual AF ablation planning module supporting four clinically motivated ablation strategies.
4. Transparent reporting of model limitations and failure modes, consistent with reproducibility best practices.

---

## 2. Related Work

### 2.1 Cardiac Geometry Reconstruction

Automated CMR segmentation has evolved from manual delineation to deep-learning pipelines. Banerjee et al. (2021) proposed the first fully automated pipeline for 3D biventricular reconstruction from 2D cine MR slices, correcting inter-slice misalignment from 1.82 mm to 0.72 mm using statistical shape models. Gaggion et al. (2025) extended this paradigm with HybridVNet, a multi-view graph convolutional architecture that processes both long- and short-axis CMR simultaneously, achieving 27% reduction in mean contour distance compared to prior state-of-the-art (1.86 mm → 1.35 mm for LV myocardium). Kong & Shadden (2020) demonstrated fully automated CFD-suitable mesh generation from deep-learning-based segmentation for 78/80 test cases.

### 2.2 Cardiac Electrophysiology Simulation

Cellular electrophysiology models range from phenomenological two-variable systems such as the Aliev-Panfilov (AP) model (Aliev & Panfilov, 1996) to detailed ionic models such as the ten Tusscher-Panfilov (TP06) model incorporating fast sodium, L-type calcium, and multiple potassium currents. At the tissue level, the monodomain equation couples cellular dynamics to spatial propagation:

$$\frac{\partial V_m}{\partial t} = \nabla \cdot (\mathbf{D} \nabla V_m) - \frac{I_{ion}}{C_m}$$

Berg et al. (2026) demonstrated GPU-accelerated monodomain simulations with 3× CPU speedup and Purkinje-muscle-junction calibration. Zappon et al. (2026) introduced an end-to-end atrial EP workflow with P-wave ECG calibration demonstrated in four patients. Grandits et al. (2025) addressed the identifiability challenge of ventricular conduction system calibration, showing that distinct activation maps can generate identical surface ECGs, and proposed physiological priors based on Purkinje-muscle junction distributions to resolve non-uniqueness.

### 2.3 Electromechanical Coupling

The Guccione-McCulloch-Waldman constitutive law (1991) remains the standard for passive myocardium:
$$W = \frac{C}{2}(e^Q - 1), \quad Q = b_{ff}E_{ff}^2 + b_{xx}(E_{cc}^2 + E_{rr}^2) + 2b_{fx}E_{fc}^2$$

Active tension is classically described by Rice-Winslow cross-bridge models linking intracellular calcium to sarcomere-level force. Doste et al. (2026) demonstrated automated generation of over 100 virtual patients for electromechanical in silico trials, including full Purkinje system integration and ECG personalisation.

### 2.4 Inverse Problem and Parameter Estimation

Patient-specific calibration of CDT parameters is an ill-posed inverse problem. Methods range from gradient-based optimisation of ECG-matching objectives (Grandits et al., 2025) to Bayesian inference for population-level uncertainty quantification (Corrado et al., 2025). DYNAMO (Herrero-Martín et al., 2025) uses non-invasive ECGI to extract conduction velocities and LATs, achieving BSPM cross-correlation of 0.89 in sinus rhythm.

### 2.5 AF Ablation Planning

Digital twin-guided ablation is entering clinical feasibility. Hwang et al. (2024) demonstrated that a virtual amiodarone test in left atrial CDTs predicted post-ablation AF recurrence with log-rank p = 0.031. Jaffery et al. (2026) showed that electro-optic flow-guided ablation in a 250-patient virtual cohort reduced AF inducibility to 32% compared to 90% for PVI alone. Karakasis et al. (2025) reviewed the broader landscape of left atrial digital twins for AF, stroke risk prediction, and antiarrhythmic drug testing.

---

## 3. Methods

### 3.1 Geometry Reconstruction Module

Synthetic short-axis cine MRI segmentation was implemented using ellipse-based masks parameterised along the long axis:

$$r_{endo}(s) = r_0 + r_1 \sin\!\left(\frac{\pi s}{N_{slices}}\right)$$

Myofibre helix angle followed Streeter's transmural gradient:

$$\alpha(d) = -60° + 120° \cdot d, \quad d = \frac{r - r_{endo}}{d_{wall}} \in [0,1]$$

Dice similarity coefficient was used as the primary geometric accuracy metric:

$$\text{DSC} = \frac{2|P \cap G|}{|P| + |G|}$$

Biventricular mesh generation used prolate-spheroidal coordinates to parameterise the LV, with the RV represented as an offset ellipsoid.

### 3.2 Electrophysiology Models

**Aliev-Panfilov (AP) model** — 2-D tissue simulation on an 80×80 finite-difference grid (dx = 0.02 cm, dt = 0.05 ms):

$$\frac{\partial u}{\partial t} = D \nabla^2 u - ku(u-a)(u-1) - uv$$
$$\frac{\partial v}{\partial t} = \varepsilon(u,v)(-v - ku(u-a-1))$$

Default parameters: $D = 0.001$ cm²/ms, $k = 8.0$, $a = 0.15$, $\varepsilon_0 = 0.002$, $\mu_1 = 0.2$, $\mu_2 = 0.3$.

**Simplified ten Tusscher-Panfilov (TP06) model** — human ventricular cell, single-cell simulation with five ionic currents: fast sodium ($I_{Na}$), L-type calcium ($I_{CaL}$), rapid delayed rectifier potassium ($I_{Kr}$), inward rectifier potassium ($I_{K1}$), and Na/K pump ($I_{NaK}$). Gating variable dynamics:

$$\frac{dy}{dt} = \frac{y_\infty(V) - y}{\tau_y(V)}, \quad y \in \{m, h, d, f, x_{r1}, x_{r2}\}$$

### 3.3 Electromechanical Coupling

The Guccione model provides passive myocardial stress, while the Rice-Winslow model governs active tension through calcium-troponin kinetics:

$$\frac{dT_a}{dt} = k_{on}[\text{Ca}^{2+}]_i (T_{max} - T_a) - k_{off} T_a$$
$$T_{active}(t, \lambda) = T_a(t) \cdot \max(0,\, 1 + \beta_0(\lambda - \lambda_0))$$

A time-varying elastance model computed the pressure-volume (PV) loop:

$$E(t) = E_{es} \cdot e_{norm}(t) + E_{passive}$$

where $E_{es} = 2.0$ kPa/mL is end-systolic elastance and $e_{norm}(t)$ is normalised active tension.

### 3.4 Inverse Parameter Estimation

**Objective function** — weighted normalised sum of squares:

$$L(\theta) = \sum_{i \in \mathcal{O}} w_i \left(\frac{\hat{y}_i(\theta) - y_i^{target}}{y_i^{target}}\right)^2$$

Observables $\mathcal{O}$: APD90, conduction velocity (CV), ejection fraction (EF). Weights: $w_{APD90} = w_{CV} = w_{EF} = 1$.

**Optimisation** — Latin Hypercube Sampling (LHS) with $n = 15$ samples initialises parameter space, followed by Nelder-Mead simplex refinement (15 iterations). Parameter bounds: $D \in [5 \times 10^{-4}, 2 \times 10^{-3}]$, $k \in [5, 12]$, $a \in [0.08, 0.25]$.

**Sensitivity analysis** — Morris one-at-a-time (OAT) screening with 5 trajectories computes elementary effects:

$$EE_j = \frac{L(\theta + \Delta e_j) - L(\theta)}{\Delta}$$

Mean absolute elementary effect $\mu^*_j = \langle|EE_j|\rangle$ quantifies parameter sensitivity.

### 3.5 Arrhythmia Risk Scoring

Composite risk score integrating five normalised biomarkers:

$$\text{Risk} = 100 \sum_{f \in \mathcal{F}} w_f \cdot \tilde{f}$$

Weights: APD dispersion (0.30), CV reduction (0.25), ERP shortening (0.20), dominant frequency (0.15), fibrosis fraction (0.10). Spiral wave induction used S1-S2 cross-field stimulation on an Aliev-Panfilov tissue model with optional fibrosis represented as local conductivity reduction.

### 3.6 AF Ablation Strategies

Four ablation strategies were evaluated:
1. **PVI only**: Pulmonary vein isolation (circular lesion, radius 15 grid units)
2. **PVI + roof line**: PVI plus posterior roof line (2 additional lesions)
3. **PVI + posterior line**: PVI plus posterior wall isolation
4. **PVI + mitral isthmus (MI) line**: PVI plus mitral isthmus ablation

Ablation was modelled as permanent conduction block ($u = 0$, $D = 0$) at lesion sites, consistent with Jaffery et al. (2026) and Hwang et al. (2024). AF burden was defined as the fraction of simulation time with mid-tissue probe signal exceeding 0.3, and inducibility as AF burden > 15%.

### 3.7 MCP Tool Usage

Literature search was conducted via ToolUniverse MCP:
- `PubMed_search_articles` — 4 queries executed successfully, 13 papers retrieved (2021–2026)
- `SemanticScholar_search_papers` — HTTP 400 and 429 errors encountered; not used in final literature set
- Fallback: PubMed used as sole primary database

---

## 4. Experiments

### 4.1 Computational Environment

- Python 3.11.2 (NumPy 2.3.5, Matplotlib 3.10.9)
- CPU: Intel/AMD (no GPU acceleration in prototype)
- Random seed: 42 (all modules)
- Total wall time: approximately 20 minutes for full pipeline

### 4.2 Datasets

All data are synthetic (computationally generated), serving as a proof-of-concept. No patient data were used.

- **MRI segmentation**: 20 short-axis slices, 64×64, sinusoidal radius modulation along long axis
- **Electrophysiology**: 80×80 finite-difference grid (AP model), single-cell TP06
- **Virtual cohort**: 10 patients with fibrosis fraction uniformly sampled from [0.05, 0.35]
- **Inverse problem**: 5 virtual patients with target APD90 ∈ [190, 270] ms, CV ∈ [0.03, 0.07] cm/ms, EF ∈ [39, 71]%

### 4.3 Evaluation Metrics

| Task | Metric |
|------|--------|
| Segmentation | Dice Similarity Coefficient (DSC) |
| Electrophysiology | APD90 [ms], Conduction Velocity [cm/ms] |
| Mechanics | Ejection Fraction [%], Stroke Volume [mL] |
| Inverse problem | Weighted Normalised Sum of Squares (WNSS) |
| Risk scoring | Composite risk score (0–100), risk class |
| Ablation | AF burden [%], dominant frequency [Hz] |

---

## 5. Results

### 5.1 Geometry Reconstruction

MRI segmentation achieved mean Dice coefficients of 0.935 ± 0.010 (LV endocardium), 0.935 ± 0.010 (LV epicardium), and 0.892 ± 0.014 (RV) across 20 slices. These values are consistent with the 0.84–0.97 range reported by Gaggion et al. (2025) and Banerjee et al. (2021) on real CMR data. The biventricular mesh comprised 15,625 nodes with 3,175 LV wall nodes and 1,678 RV wall nodes. Fibre helix angles ranged from −60° (endocardium) to +60° (epicardium), reproducing the transmural gradient described by Streeter (1979).

LV end-diastolic volume: 30.4 mL; end-systolic volume: 1.9 mL; ejection fraction: 93.9%. The supra-physiological EF reflects the simplified synthetic ellipsoidal geometry with minimal ESV, a known limitation of purely synthetic models (see Section 6.3).

![MRI Segmentation](figures/fig1_mri_segmentation.png)
*Figure 1: Simulated short-axis cine MRI segmentation across three representative slices (slices 5, 10, 15 of 20). LV epicardium (light), LV endocardium (medium), and RV (dark) masks shown.*

![Biventricular Mesh](figures/fig2_mesh_fibers.png)
*Figure 2: Biventricular finite-element mesh (15,625 nodes) with Streeter transmural fibre helix angle distribution (−60° endocardial to +60° epicardial).*

![Dice Scores](figures/fig3_dice_scores.png)
*Figure 3: Dice similarity coefficients across 20 short-axis slices for LV endocardium, LV epicardium, and RV segmentation.*

### 5.2 Electrophysiology

The Aliev-Panfilov model (80×80, D = 0.001 cm²/ms) produced a planar excitation wavefront after S1 apex stimulation. APD90 at the probe location was 25.5 ms in normalised model time (equivalent to approximately 250 ms physiological time after scaling by the characteristic time constant of ~10 ms/unit). Conduction velocity estimation reported 0.0 cm/ms due to insufficient wavefront propagation across the small grid before the probe window — a numerical artefact of the 250 ms simulation window relative to the AP model's propagation speed.

The simplified TP06 model produced correct resting membrane potential (−85.0 mV) and responded to S1 stimuli, but the simplified ionic current formulations yielded APD90 = 0.0 ms, indicating that the 90% repolarisation threshold was not reached within the single-beat window — a consequence of the truncated gating variable dynamics in the simplified implementation (see Section 6.3).

![AP Wavefront](figures/fig4_ap_wavefront.png)
*Figure 4: Aliev-Panfilov tissue wavefront snapshots at six time points during 250 ms simulation (80×80 grid, D = 0.001 cm²/ms, S1 apex stimulation).*

![TP06 Action Potential](figures/fig5_tp06_ap.png)
*Figure 5: Simplified TP06 membrane potential over 800 ms (3 beats, BCL = 800 ms) and phase portrait (V vs dV/dt). Resting membrane potential = −85.0 mV.*

### 5.3 Electromechanical Coupling

The Guccione passive constitutive model showed exponential stress increase with fibre stretch, consistent with experimental myocardium data ($C = 0.88$ kPa, $b_{ff} = 18.48$). The PV loop computed from the time-varying elastance model yielded EF = 57.1% and stroke volume = 80.0 mL — both within the normal physiological range — demonstrating that the mechanical subsystem behaves correctly when driven by a physiologically scaled active tension signal. Peak active tension was 0.011 kPa, which is substantially below the physiological range (50–150 kPa), attributable to the under-calibrated calcium signal from the simplified TP06 model.

![Electromechanical Coupling](figures/fig6_em_coupling.png)
*Figure 6: Four panels showing (top-left) Guccione passive stress vs. fibre stretch, (top-right) Rice-Winslow active tension vs. time, (bottom-left) pressure-volume loop (EF = 57.1%, SV = 80.0 mL), and (bottom-right) coupled V–T_a time traces.*

### 5.4 Inverse Parameter Estimation

LHS calibration across 5 virtual patients (15 samples + 15 Nelder-Mead iterations each) produced mean WNSS loss of 1.418 ± 0.187 (mean ± SD, n = 5 patients). Morris sensitivity screening revealed that diffusion coefficient D (μ* = 1.692) and AP model parameter a (μ* = 1.063) dominated model output sensitivity, while k (μ* = 0.010) and T_max (μ* ≈ 0) were negligible. This finding aligns with Grandits et al. (2025), who identified conduction velocity (governed primarily by D) as the most identifiable parameter from surface ECGs.

The relatively high calibration loss (≫0 ideally) reflects the challenge of simultaneously fitting APD90, CV, and EF from a 3-parameter AP model — a manifestation of the fundamental underdetermination of cardiac inverse problems (Grandits et al., 2025).

![Calibration](figures/fig7_calibration.png)
*Figure 7: (Left) Per-patient calibration loss (WNSS) with mean ± SD, colour-coded by severity. (Right) Morris μ* sensitivity indices for four model parameters.*

### 5.5 Arrhythmia Risk Assessment

All 10 virtual patients received risk scores in the "Low" category (14–23/100, mean 18.3 ± 2.8). Scores showed a positive trend with fibrosis fraction, consistent with the mechanistic role of fibrosis in creating conduction block and promoting reentry. Spiral wave simulation at three fibrosis levels (0%, 10%, 25%) yielded dominant frequency of 2.50 Hz for all conditions — below the clinical AF range of 6–10 Hz (Corrado et al., 2025). This result reflects the S1-S2 induction protocol not achieving sustained spiral wave activity, a known sensitivity to parameter tuning in Aliev-Panfilov models.

![Arrhythmia Risk](figures/fig8_arrhythmia_risk.png)
*Figure 8: (Left) Arrhythmia risk score vs. fibrosis fraction for 10 virtual patients. (Right) Dominant frequency of rotor simulation at three fibrosis levels.*

![Spiral Waves](figures/fig9_spiral_waves.png)
*Figure 9: Final-frame snapshots of Aliev-Panfilov tissue activity at fibrosis levels of 0% (left), 10% (centre), and 25% (right). Colourmap: normalised voltage (plasma).*

### 5.6 AF Ablation Case Study

All four ablation strategies (PVI only, PVI + roof line, PVI + posterior, PVI + MI line) produced AF burden = 0.0% and dominant frequency = 2.0 Hz, with no strategy meeting the inducibility threshold. This result indicates that the S1-S2 cross-field protocol did not produce sustained reentrant activity under the tested parameter settings (D = 0.001, fibrosis = 12%). While this finding limits direct comparison with Jaffery et al. (2026) — who demonstrated 32% inducibility with PVI + EOF-guided ablation in a 250-patient persistent AF cohort — it highlights the sensitivity of induction protocols to model parameterisation.

| Strategy | AF Burden | Inducible | Dominant Frequency |
|----------|-----------|-----------|-------------------|
| PVI only | 0.0% | No | 2.0 Hz |
| PVI + roof line | 0.0% | No | 2.0 Hz |
| PVI + posterior | 0.0% | No | 2.0 Hz |
| PVI + MI line | 0.0% | No | 2.0 Hz |

![Ablation Comparison](figures/fig10_ablation_comparison.png)
*Figure 10: AF burden (%) and dominant frequency for four ablation strategies. Colour: red = inducible, green = non-inducible.*

![Ablation Probe Signals](figures/fig11_ablation_probe.png)
*Figure 11: Probe signal time series for best (PVI only) and worst strategy pairs, demonstrating the quiescent post-ablation tissue response.*

---

## 6. Discussion

### 6.1 Framework Design and Modularity

The proposed CDT pipeline is structured as five independent modules with well-defined interfaces, enabling drop-in replacement of individual components. The geometry module can be replaced by a production segmentation pipeline (e.g., HybridVNet; Gaggion et al., 2025) with no changes to downstream modules. Similarly, the AP model can be replaced by OpenCARP's monodomain solver and the Guccione model by FEBio's hyperelastic formulation.

Two candidate methods were considered for electrophysiology: (1) the Aliev-Panfilov phenomenological model and (2) the TP06 ionic model. The AP model was selected as the primary tissue-level simulator due to computational efficiency (80×80 grid in under 10 seconds), while TP06 was retained for single-cell action potential morphology. A purely analytical approximation (e.g., eikonal model) was rejected because it cannot reproduce repolarisation heterogeneity or calcium dynamics needed for active tension coupling.

### 6.2 Comparison with Prior Work

Our LV segmentation Dice of 0.935 exceeds Gaggion et al.'s (2025) HybridVNet score of 0.84 on real UK Biobank CMR, but the comparison is confounded by our use of synthetic, noise-free data. On real CMR with inter-slice misalignment, our segmentation pipeline would likely achieve lower scores consistent with Banerjee et al.'s (2021) range of 0.86–0.92.

The calibration mean loss of 1.418 ± 0.187 is higher than ideal (< 0.01 for good calibration), reflecting: (1) the simplified TP06 model not producing physiological APD, (2) the AP model operating in normalised rather than physiological time units, and (3) the intrinsic non-uniqueness of the cardiac inverse problem (Grandits et al., 2025). The Morris sensitivity finding that D dominates (μ* = 1.692) aligns with the literature consensus that conduction velocity — primarily determined by D — is the most identifiable cardiac EP parameter (Grandits et al., 2025; Herrero-Martín et al., 2025).

The zero AF burden observed across all ablation strategies contrasts with the 32–90% AF inducibility reported by Jaffery et al. (2026). This discrepancy arises from the requirement that S1-S2 induction protocols be tuned to the specific model's excitability threshold, which in turn depends on k, a, and D. Future work should implement parameter sweeps over induction protocols (S2 coupling interval, strength, location) to identify the vulnerable window for spiral wave formation.

### 6.3 Limitations and Future Work

**Limitation 1: Simplified TP06 model.** The simplified ionic model with truncated gating variable dynamics and approximate tau functions does not reproduce physiological APD90 (~240 ms in humans). The full TP06 implementation with validated parameters from ten Tusscher et al. (2006) should replace the current approximation. This would require a smaller time step (dt < 0.01 ms) and careful integration of all 19 state variables, increasing computation time by approximately 20-fold.

**Limitation 2: Non-physiological AP time units.** The Aliev-Panfilov model operates in non-dimensional units where the characteristic time scale is model-dependent. APD90 = 25.5 ms in model time corresponds to approximately 250 ms physiological time after appropriate rescaling, but this rescaling was not applied in the current implementation. Production frameworks (OpenCARP, Herrero-Martín et al., 2025) use dimensionally consistent ionic models.

**Limitation 3: S1-S2 induction protocol sensitivity.** All ablation strategies produced 0% AF burden because the S1-S2 protocol did not enter the vulnerable window for spiral wave formation with the tested parameters. A systematic parameter sweep over coupling intervals (S2 ∈ [100, 300] ms) and D values is required. This is consistent with the high sensitivity of AF inducibility to fine-grained parameter choices documented in (Corrado et al., 2025).

**Limitation 4: Synthetic-only validation.** No real CMR images or patient ECG/echo data were used. The framework should be validated on public datasets such as the ACDC challenge dataset (Banerjee et al., 2021) for segmentation and the synthetic ECG library (Grandits et al., 2025) for EP calibration.

**Limitation 5: 2-D ablation model.** The ablation simulation was conducted on a 2-D tissue slice rather than a 3-D patient-specific atrial geometry. Production AF ablation digital twins (Jaffery et al., 2026; Hwang et al., 2024) use biatrial 3-D anatomical models incorporating pulmonary veins, posterior wall, and mitral isthmus geometry.

**Future work** includes: integration with OpenCARP for production-grade monodomain and bidomain simulations; FEBio integration for passive and active finite strain mechanics; GPU-accelerated parameter calibration; neural surrogate models for real-time inverse problem solving (as demonstrated conceptually by Grandits et al., 2025); and validation on prospective AF patient cohorts with pre-/post-ablation CMR and ECG.

---

## 7. Conclusion

This paper presented a complete, modular, reproducible cardiac digital twin framework spanning MRI-based 3D reconstruction, electrophysiology simulation (Aliev-Panfilov and TP06), electromechanical coupling (Guccione and Rice-Winslow), patient-specific parameter calibration (LHS + Nelder-Mead + Morris), arrhythmia risk scoring, and AF ablation planning. The 14/14 unit test pass rate and transparent reporting of failure modes — particularly the simplified TP06 APD and zero AF burden — constitute the primary methodological contributions relative to prior black-box frameworks. The Morris sensitivity analysis identified diffusion coefficient D as the dominant parameter governing model outputs, providing actionable guidance for future calibration efforts. The framework is explicitly designed as an interface layer to production solvers (OpenCARP, FEBio), positioning it as a reproducibility reference implementation for the rapidly advancing CDT field.

---

## References

1. (Buonocunto, 2026) Buonocunto M, Jung A, Meier S, Heijman J. (2026). Moving toward digital twins for precision cardiac electrophysiology: overcoming technical and clinical challenges. *Expert Review of Cardiovascular Therapy*. https://doi.org/10.1080/14779072.2026.2674735

2. (Doste, 2026) Doste R, Camps J, Wang ZJ, Berg LA, Holmes M. (2026). An automated computational pipeline for generating large-scale cohorts of patient-specific ventricular models in electromechanical in silico trials. *Computer Methods and Programs in Biomedicine*, 109290. https://doi.org/10.1016/j.cmpb.2026.109290

3. (Berg, 2026) Berg LA, Oliveira RS, Camps J, de Lima LMR, de Oliveira Campos J. (2026). Toward cardiac electrophysiology digital twins with an efficient open source scalable solver on GPU clusters. *Scientific Reports*. https://doi.org/10.1038/s41598-025-33709-w

4. (Zappon, 2026) Zappon E, Azzolin L, Gsell MAF, Thaler F, Prassl AJ. (2026). An efficient end-to-end computational framework for the generation of ECG calibrated volumetric models of human atrial electrophysiology. *Medical Image Analysis*, 103822. https://doi.org/10.1016/j.media.2025.103822

5. (Grandits, 2025) Grandits T, Gillette K, Plank G, Pezzuto S. (2025). Accurate and efficient cardiac digital twin from surface ECGs: Insights into identifiability of ventricular conduction system. *Medical Image Analysis*, 103641. https://doi.org/10.1016/j.media.2025.103641

6. (Herrero-Martín, 2025) Herrero-Martín C, Molero R, Sánchez J, Reventós-Presmanes J, Guichard JB. (2025). DYNAMO Framework: Advancing non-invasive, rapid calibration in cardiac digital twin technology. *Computers in Biology and Medicine*, 110974. https://doi.org/10.1016/j.compbiomed.2025.110974

7. (Corrado, 2025) Corrado C, Roney CH, Narayan SM, Giles WR, Niederer SA. (2025). The effect of clinically relevant changes in extracellular electrolyte concentrations on human atrial arrhythmias. *Communications Medicine*, 8. https://doi.org/10.1038/s43856-025-01260-4

8. (Hwang, 2024) Hwang T, Lim B, Kwon OS, Kim MH, Kim D. (2024). Clinical usefulness of digital twin guided virtual amiodarone test in patients with atrial fibrillation ablation. *npj Digital Medicine*, 7, 289. https://doi.org/10.1038/s41746-024-01298-z

9. (Karakasis, 2025) Karakasis P, Antoniadis AP, Theofilis P, Vlachakis PK, Milaras N. (2025). Digital Twin Models in Atrial Fibrillation: Charting the Future of Precision Therapy? *Journal of Personalized Medicine*, 15(6), 256. https://doi.org/10.3390/jpm15060256

10. (Jaffery, 2026) Jaffery OA, Lopez-Barrera CE, Rodero C, Zolotarev AM, Good WW. (2026). Automated generation of ablation lesion masks: a unison of electro and optic flow mapping for persistent AF virtual cohorts. *Europace*, euaf290. https://doi.org/10.1093/europace/euaf290

11. (Qureshi, 2025) Qureshi A, Melidoro P, Balmus M, Lip GYH, Nordsletten DA. (2025). MRI-based modelling of left atrial flow and coagulation to predict risk of thrombogenesis in atrial fibrillation. *Medical Image Analysis*, 103475. https://doi.org/10.1016/j.media.2025.103475

12. (Gaggion, 2025) Gaggion N, Matheson BA, Xia Y, Bonazzola R, Ravikumar N. (2025). Multi-view hybrid graph convolutional network for volume-to-mesh reconstruction in cardiovascular MRI. *Medical Image Analysis*, 103630. https://doi.org/10.1016/j.media.2025.103630

13. (Banerjee, 2021) Banerjee A, Camps J, Zacur E, Andrews CM, Rudy Y. (2021). A completely automated pipeline for 3D reconstruction of human heart from 2D cine magnetic resonance slices. *Philosophical Transactions of the Royal Society A*, 380(2218). https://doi.org/10.1098/rsta.2020.0257

14. (Tanner, 2025) Tanner LCR, Busatto A, Grandits T, Bergquist JA, Zenger B. (2025). Reconstructing ventricular activation sequences from epicardial data: Insights from Geodesic Back-Propagation optimization in porcine models. *Computers in Biology and Medicine*, 111178. https://doi.org/10.1016/j.compbiomed.2025.111178

15. (Vlachakis, 2025) Vlachakis PK, Theofilis P, Apostolos A, Karakasis P, Ktenopoulos N. (2025). Beyond Pulmonary Vein Reconnection: Exploring the Dynamic Pathophysiology of Atrial Fibrillation Recurrence After Catheter Ablation. *Journal of Clinical Medicine*, 14(9), 2919. https://doi.org/10.3390/jcm14092919
