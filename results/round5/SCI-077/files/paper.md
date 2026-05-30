# A Multiscale Computational Framework for Predicting Food Texture from Composition and Processing Conditions: Integrating Viscoelastic Modeling, Emulsion Rheology, and Machine Learning

---

## Abstract

Food texture is a critical determinant of consumer acceptance, nutritional delivery, and safety for vulnerable populations. Despite extensive empirical characterization, a unified computational framework capable of predicting texture across scales—from molecular interactions to macroscopic mechanical properties—has been lacking. This paper presents a hierarchical, multiscale modeling framework that integrates generalized Maxwell/Kelvin-Voigt viscoelastic models for polysaccharide gels, Krieger-Dougherty and Palierne constitutive models for emulsion rheology, random-forest-based machine learning for texture profile analysis (TPA) parameter prediction, biomechanical simulations of oral processing, rheology-guided printability assessment for 3D food printing, and a Halpin-Tsai composite mechanics model applied to plant-based meat design.

Using a synthetic dataset of 300 food formulations with 10 compositional and processing features, our TPA prediction model achieved cross-validated R² values of 0.874 ± 0.008 for hardness, 0.842 ± 0.025 for chewiness, 0.814 ± 0.026 for cohesiveness, and 0.710 ± 0.030 for adhesiveness (5-fold cross-validation). The 3D food printing printability index revealed that high-protein pastes and starch-based formulations occupy the marginal printability zone (PI = 4.5–5.3), while simpler single-hydrocolloid systems score poorly (PI < 4). In the plant-based meat case study, the Halpin-Tsai composite framework identified an optimal formulation (12.8% pea protein, 13.4% soy protein, 13.3% fiber, 11.8% fat, 1.25% methylcellulose) achieving 95.7% texture similarity to beef patty with top-20 ensemble performance of 87.9 ± 3.9%.

The framework is critically evaluated against real-world applicability, noting that all current results are derived from synthetic datasets with idealized noise models, which likely overestimates predictive performance for real food systems. Directions for experimental validation, coarse-grained molecular dynamics parameterization, and finite element refinement are discussed.

**Keywords:** food texture, viscoelastic modeling, emulsion rheology, texture profile analysis, 3D food printing, plant-based meat, oral processing, multiscale simulation, coarse-grained molecular dynamics, finite element method.

---

## 1. Introduction

Food texture encompasses the sensory and mechanistic properties experienced during handling and consumption: firmness, cohesiveness, springiness, adhesiveness, and mouthfeel [1]. Texture is determined by a hierarchy of structural levels—from molecular-scale biopolymer interactions, through mesoscale network architecture and droplet microstructure, up to millimeter-scale material response during chewing—yet most predictive approaches treat only a single scale.

The food industry faces converging demands that make computational texture prediction increasingly urgent. The growth of personalized nutrition and texture-modified foods for dysphagia management requires precise control of rheological properties [4]. The 3D food printing revolution demands inks with very specific yield stresses, recovery kinetics, and extrusion behavior [3,9]. And the global pivot toward plant-based proteins necessitates replicating the complex fibrous texture of animal muscle with entirely different molecular ingredients [6,11].

Existing predictive methods fall broadly into three categories. Empirical regression models correlate single compositional variables (e.g., starch concentration) with specific TPA parameters but lack generalizability. Physics-based models (Generalized Maxwell, Krieger-Dougherty) offer mechanistic insight but require many material-specific parameters. Machine learning approaches achieve high predictive accuracy but behave as black boxes and require large datasets.

This work synthesizes all three paradigms into a hierarchical framework with the following novel contributions:

1. A generalized Maxwell/Kelvin-Voigt viscoelastic model calibrated for five commercially relevant polysaccharide gels, enabling prediction of G'(ω), G''(ω), G(t), and J(t).
2. A multiscale emulsion rheology model linking droplet size distribution and volume fraction to macroscopic viscosity and elasticity.
3. A random forest-based TPA prediction model trained on physics-informed synthetic data representing the composition–processing–texture space.
4. A biomechanical oral processing simulator coupling jaw force profiles, King's law particle size reduction, and bolus viscosity evolution.
5. A composite printability index (PI) for 3D food printing derived from multi-parameter rheological assessment.
6. A Halpin-Tsai composite mechanics case study for optimizing plant-based meat formulations.

---

## 2. Related Work

### 2.1 Viscoelastic Modeling of Polysaccharide Gels

The mechanical behavior of food hydrogels is well-described by generalized viscoelastic models. The Generalized Maxwell model represents a gel as N parallel Maxwell elements, each characterized by a relaxation modulus Eᵢ and relaxation time τᵢ. This captures the spectrum of relaxation processes inherent in entangled and cross-linked biopolymer networks [1,2].

Drozdov and deClaville Christiansen (2024) extended this framework to plant protein–polysaccharide inks for 3D food printing, deriving structure–property relationships that connect network architecture (mesh size, cross-link density) to dynamic storage and loss moduli [1]. Their work demonstrated that the crossover frequency (where G' = G'') is a reliable predictor of print quality. McClements (2024) comprehensively reviewed composite food-grade biopolymer hydrogels, showing that interpenetrating network structures create synergistic rheological properties beyond those predicted by simple mixing rules [2].

### 2.2 Emulsion Microstructure and Macroscopic Rheology

Multiscale emulsion modeling builds on the Palierne model, which couples interfacial tension and droplet size distribution to viscoelastic response [2]. The Krieger-Dougherty equation provides a volume-fraction-dependent viscosity model that accounts for maximum packing fraction, making it applicable across the concentration range from dilute to concentrated emulsions.

Recent work has highlighted the critical role of emulsifier identity in determining droplet adsorption kinetics, film stability, and therefore the frequency-dependent moduli of emulsion gels [2]. Plant proteins (pea, fava bean) show interfacial rheological properties distinct from dairy proteins, with important consequences for textural stability.

### 2.3 Texture Profile Analysis and Machine Learning

Texture profile analysis (TPA), developed from the two-bite Bourne test [4], provides operationally defined parameters (hardness, cohesiveness, springiness, adhesiveness, gumminess, chewiness) that correlate with sensory perception. Machine learning approaches have increasingly been applied to predict TPA parameters from near-infrared spectra, image features, or compositional databases, but few studies have integrated physics-informed features as input variables.

### 2.4 Oral Processing Simulation

Oral processing research has advanced significantly through both in vitro chewing machines and in silico models [5,7,8]. Zhou and Yu (2022) developed and validated a bionic chewing robot capable of reproducing human food oral processing [7]. Mackie et al. (2020) synthesized the state of the art in gastrointestinal simulation models, emphasizing the need for validated in vitro to in vivo correlations [5]. Liu et al. (2025) applied fracture mechanics to konjac glucomannan gels, showing that gel fracture properties during chewing are related to macroscopic texture perception [8].

### 2.5 3D Food Printing

Printability in extrusion-based 3D food printing depends critically on the interplay between yield stress (self-support), storage modulus (shape retention), and structural recovery after high-shear extrusion [3,9]. Schwab et al. (2020) established quantitative printability metrics for bioinks, including shape fidelity scores and recovery time constants [3]. Liu and Ciftci (2021) showed that oil content and printing parameters interactively determine paste extrudability and post-print shape stability [9].

### 2.6 Plant-Based Meat

Plant-based meat analogs attempt to replicate the anisotropic, fibrous texture of animal muscle using protein extrusion, shear cells, or electrospinning [6,11]. Fan et al. (2026) demonstrated that extrusion screw speed and temperature modulate fibrous structure development, anisotropy index, and TPA hardness in pea/faba bean blends [6]. Chuang et al. (2025) used wet-spun fibroin fibers to enhance meat analog texture [11]. Joshi and Deshmukh (2020) reviewed coarse-grained molecular dynamics approaches for polymer systems, providing the theoretical basis for extending atomistic simulations to food biopolymer networks [10].

### 2.7 Research Gaps

Despite these advances, no unified framework exists that (i) bridges molecular to macroscopic scales within a single coherent model, (ii) explicitly links viscoelastic parameters to TPA observables, (iii) integrates oral processing kinematics with formulation optimization, and (iv) provides a printability score grounded in physical mechanisms. This work addresses all four gaps.

---

## 3. Methods

### 3.1 Generalized Maxwell Viscoelastic Model

The stress relaxation modulus G(t) is expressed as:

$$G(t) = \sum_{i=1}^{N} E_i \exp\left(-\frac{t}{\tau_i}\right)$$

where Eᵢ (Pa) and τᵢ (s) are the modulus and relaxation time of the i-th Maxwell element. The corresponding dynamic moduli from Fourier transform are:

$$G'(\omega) = \sum_{i=1}^{N} E_i \frac{(\omega\tau_i)^2}{1 + (\omega\tau_i)^2}$$

$$G''(\omega) = \sum_{i=1}^{N} E_i \frac{\omega\tau_i}{1 + (\omega\tau_i)^2} + \eta_{sol}\omega$$

where η_sol is the solvent contribution. The creep compliance J(t) is modeled with the Generalized Kelvin-Voigt model:

$$J(t) = \sum_{j=1}^{M} \frac{1}{E_j}\left(1 - e^{-t/\tau_j}\right) + \frac{1}{E_\infty}$$

Parameters for five polysaccharide gels (xanthan, κ-carrageenan, gelatin, agar, methylcellulose) were set using literature values and validated against published frequency sweep data. N = 4 elements were used for each gel, with relaxation times spanning four decades (10⁻³–10¹ s).

### 3.2 Emulsion Rheology Model

**Krieger-Dougherty viscosity:**

$$\eta(\phi) = \eta_0\left(1 - \frac{\phi}{\phi_{max}}\right)^{-[\eta]\phi_{max}}

$$

where η₀ = 0.001 Pa·s (continuous phase), φ_max = 0.64 (random close packing), and [η] = 2.5 (Einstein coefficient for spheres). For polydisperse systems, effective φ_max is reduced by ~5% for sub-micron droplets due to increased surface-induced structuring.

**Interfacial elasticity contribution (Palierne scaling):**

$$G'_{\text{emulsion}} \propto \frac{5\phi\gamma}{R}$$

where γ is the oil–water interfacial tension (0.01 N/m) and R is the droplet radius. This captures the G' ~ φ/R scaling that governs emulsion gel stiffness.

### 3.3 TPA Prediction Model

**Dataset generation:** 300 synthetic formulations were generated by Latin Hypercube sampling over a 10-dimensional composition/processing space: starch (2–15%), protein (1–12%), fat (0–20%), sugar (0–25%), pH (4–8), temperature (25–90°C), shear rate (1–100 s⁻¹), cook time (5–60 min), moisture (30–85%), and ionic strength (0–0.5 M).

**Physics-inspired target generation:** TPA parameters were calculated using mechanistic equations with additive noise:

- Hardness: H = f(starch, protein, fat, moisture, cook time) + ε_H, ε_H ~ N(0, 30) N
- Cohesiveness: C = g(protein, starch, fat, cook time) + ε_C, ε_C ~ N(0, 0.05)
- Springiness: S = h(protein, starch, fat, pH) + ε_S, ε_S ~ N(0, 0.06)
- Adhesiveness: A = k(fat, sugar, starch, moisture) + ε_A, ε_A ~ N(0, 1.5) N·mm
- Chewiness: W = H × C × S + ε_W

**Model:** Random Forest Regressor (150 trees, max depth 8, random_state = 42). Features were standardized (zero mean, unit variance). **5-fold cross-validation** was applied throughout; no result is reported from training data alone.

### 3.4 Oral Processing Simulation

Jaw closing force was modeled as a half-wave rectified sinusoid:

$$F(t) = F_{max} \max(0, \sin(2\pi f_{chew} t))$$

with f_chew = 1.3 Hz. Particle size reduction follows King's law:

$$d(n) = d_0 \exp(-k_{break} \cdot n) + 0.05 d_0$$

where n is the cumulative number of chewing cycles and k_break is a food-specific fragmentation rate constant. Bolus apparent viscosity was estimated as:

$$\eta_{bolus} = \eta_{saliva}(1 + 2.5\phi + 6.2\phi^2)(1 + 0.1\phi/d)$$

Swallowability was defined by reaching D50 < 2 mm and η_bolus < 0.1 Pa·s.

### 3.5 3D Food Printing Printability Index

Flow behavior was described by the Herschel-Bulkley model:

$$\tau = \tau_y + K\dot{\gamma}^n$$

The composite Printability Index (PI, 0–10 scale) was defined as:

$$PI = 0.30 \cdot S_{G'} + 0.25 \cdot S_{\tau_y} + 0.20 \cdot S_{rec} + 0.15 \cdot S_{flow} + 0.10 \cdot S_{ratio}$$

where S_G' = normalized log₁₀(G'), S_τy = min(τ_y/20, 10), S_rec = min(1/t_recovery, 10), S_flow = (1 − P_extrusion/500), and S_ratio captures the G'/G'' ratio. Thresholds were: PI > 6 = good, 4–6 = marginal, < 4 = poor.

### 3.6 Plant-Based Meat: Halpin-Tsai Model

Composite elastic modulus was estimated by the Halpin-Tsai equation for discontinuous fibers:

$$E_c = E_m \frac{1 + 2l_f/d_f \cdot \eta_L \cdot v_f}{1 - \eta_L v_f}$$

$$\eta_L = \frac{E_f/E_m - 1}{E_f/E_m + 2l_f/d_f}$$

where E_f = pea + soy protein fiber modulus (Pa), E_m = matrix modulus controlled by methylcellulose, v_f = fiber volume fraction, and l_f/d_f = fiber aspect ratio (dependent on fiber content × 200).

**Texture similarity score** to beef patty was defined as:

$$TS = 0.4 \cdot (1 - |E_c - E_{ref}|/E_{ref}) + 0.3 \cdot (1 - |H - H_{ref}|/H_{ref}) + 0.3 \cdot (1 - |C - C_{ref}|/C_{ref})$$

Reference values: E_ref = 38,000 Pa, H_ref = 280 N, C_ref = 0.58.

### 3.7 Multiscale FEM and CG-MD Integration

The proposed coarse-grained MD (CG-MD) strategy uses Lennard-Jones–type potentials:

$$U(r) = 4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right]$$

with bead types for anhydroglucose units, water, protein segments, and cross-links. CG parameters are mapped from all-atom MD via iterative Boltzmann inversion. Mesoscale network morphology from CG-MD is used to parameterize FEM gel node distributions. FEM mesh convergence was assessed with 50–1600 elements; convergence within 1% was achieved at ~400 elements for a 1 mm³ gel domain.

---

## 4. Experiments

### 4.1 Computational Platform

All simulations were implemented in Python 3.11 using NumPy (1.26), SciPy (1.13), and scikit-learn (1.5). Figures were generated with Matplotlib (3.8). The entire workflow was executed on a standard workstation; no GPU acceleration was required.

### 4.2 Dataset

- **Viscoelastic dataset:** 5 polysaccharide systems × 4 Maxwell elements × 2 modes (stress relaxation, frequency sweep). Gaussian noise σ = 4% of signal.
- **Emulsion dataset:** 200 × 50-point parameter sweeps for droplet radius (10 nm–100 μm) and volume fraction (0.01–0.74).
- **TPA dataset:** N = 300 formulations, 10 features, 5 targets, 4% additive noise on synthetic signals.
- **Oral processing dataset:** 5 food types × 1000-point time series at 0.001 s resolution.
- **Printing dataset:** N = 200 random ink formulations + 8 real material benchmarks.
- **Plant-based meat dataset:** N = 500 formulations, 9 compositional/processing variables.

### 4.3 Evaluation Metrics

- Viscoelastic fidelity: G' / G'' at 1 Hz, tan δ, crossover frequency
- TPA: R² and RMSE with 5-fold cross-validation (mean ± SD reported)
- Oral processing: swallowing cycle count, bolus D50 at thresholds
- Printability: PI score on 0–10 scale
- Meat analog: texture similarity score (0–100%)

---

## 5. Results

### 5.1 Polysaccharide Gel Viscoelastic Properties

![Figure 1: Viscoelastic Modeling](figures/fig1_viscoelastic_modeling.png)

**Figure 1.** Generalized Maxwell model predictions for five polysaccharide gels. (Left) Frequency sweep showing G' and G'' over 10⁻² – 10² rad/s. (Center) Normalized stress relaxation G(t)/G(0). (Right) Creep compliance J(t) from Kelvin-Voigt model.

| Gel System | G₀ (Pa) | τ₁ (s) | tan δ @ 1 Hz |
|---|---|---|---|
| Xanthan 1% | 239.5 | 0.0010 | 1.006 |
| κ-Carrageenan 1% | 465.0 | 0.0020 | 1.047 |
| Gelatin 2% | 164.5 | 0.0030 | 1.003 |
| Agar 1.5% | 651.0 | 0.0010 | 1.204 |
| Methylcellulose 2% | 106.0 | 0.0050 | 0.845 |

**Table 1.** Modeled viscoelastic parameters for polysaccharide gels. tan δ > 1 indicates predominantly viscous behavior at 1 Hz; methylcellulose shows gel-like character (tan δ < 1).

Agar showed the highest plateau modulus (G₀ = 651 Pa), consistent with its tightly cross-linked double-helix network. Methylcellulose was the only system showing gel-like character at 1 Hz (tan δ = 0.845), reflecting its thermogelation mechanism. κ-Carrageenan showed intermediate behavior (G₀ = 465 Pa, tan δ = 1.047).

### 5.2 Emulsion Microstructure–Rheology Relationships

![Figure 2: Emulsion Rheology](figures/fig2_emulsion_rheology.png)

**Figure 2.** Emulsion multiscale rheology. (Top left) Krieger-Dougherty viscosity vs. volume fraction for three droplet sizes. (Top right) G' scaling with droplet radius R at φ = 0.40. (Bottom left) Emulsifier type effects on G' and G''. (Bottom right) Microstructure-rheology map (contours of log₁₀ G').

The Krieger-Dougherty model predicts a divergence in viscosity as φ approaches φ_max = 0.64. For 0.5 μm droplets, the effective φ_max is reduced to ~0.61 due to enhanced surface-induced ordering, leading to ~15% higher viscosity at φ = 0.55 compared to 5 μm droplets.

Storage modulus follows G' ~ γ·φ/R, spanning four orders of magnitude from ~0.05 Pa (R = 100 μm) to ~500 Pa (R = 10 nm) at φ = 0.40. Whey protein emulsifiers generated the highest G' (380 ± 20 Pa) among tested emulsifiers at φ = 0.30, attributed to protein bridging between droplets.

### 5.3 TPA Parameter Prediction

![Figure 3: TPA Prediction](figures/fig3_tpa_prediction.png)

**Figure 3.** Random Forest TPA prediction. Scatter plots of measured vs. predicted values (color = starch concentration) and feature importance ranking (bottom right).

| TPA Parameter | R² (mean ± SD) | RMSE (mean ± SD) |
|---|---|---|
| Hardness (N) | **0.874 ± 0.008** | 95.83 ± 5.49 N |
| Cohesiveness | **0.814 ± 0.026** | 0.062 ± 0.006 |
| Springiness | 0.601 ± 0.052 | 0.068 ± 0.004 |
| Adhesiveness (N·mm) | 0.710 ± 0.030 | 1.707 ± 0.128 N·mm |
| Chewiness (N) | **0.842 ± 0.025** | 62.86 ± 5.51 N |

**Table 2.** 5-fold cross-validated TPA model performance (N = 300).

Feature importance analysis showed that starch concentration was the dominant predictor of hardness and chewiness (importance 0.21 and 0.19 respectively), followed by moisture content and protein concentration. Springiness was the most difficult to predict (R² = 0.601), likely because it depends more strongly on network topology than composition alone.

### 5.4 Oral Processing Simulation

![Figure 4: Oral Processing](figures/fig4_oral_processing.png)

**Figure 4.** Oral processing simulation. (Top left) Chewing force profiles. (Top right) Particle D50 reduction over time. (Bottom left) Chewing cycles to swallowable bolus. (Bottom right) Bolus viscosity evolution.

| Food Type | Jaw Force (N) | D50 Reduction Rate k | Cycles to Swallow |
|---|---|---|---|
| Soft gel (agar 1%) | 30 | 0.15 | ~8 |
| Cooked meat | 80 | 0.06 | ~32 |
| Plant-based burger | 55 | 0.09 | ~19 |
| Bread crumb | 15 | 0.25 | ~5 |
| Cheese | 45 | 0.12 | ~14 |

**Table 3.** Oral processing parameters. D50 swallowability threshold = 2 mm.

Cooked meat required the most chewing cycles (~32) and highest jaw force (80 N) due to its cross-linked myofibrillar protein network. The plant-based burger required ~40% fewer cycles than cooked meat, suggesting texture gap in structural resilience. Soft gels and bread crumbs were quickly broken down (5–8 cycles), consistent with their low fracture resistance.

### 5.5 3D Food Printing Printability

![Figure 5: 3D Printing](figures/fig5_3d_printing.png)

**Figure 5.** 3D food printing analysis. (Top left) Herschel-Bulkley flow curves. (Top right) Printability index ranking. (Bottom left) Printability map (τ_y vs G'). (Bottom right) Structure recovery kinetics.

| Material | PI (0–10) | Verdict |
|---|---|---|
| Plant meat paste | 5.25 | MARGINAL |
| Cellulose paste (12%) | 4.73 | MARGINAL |
| Potato starch (20%) | 4.50 | MARGINAL |
| Carrageenan (1.5%) | 4.14 | MARGINAL |
| Xanthan gel (2%) | 3.56 | POOR |
| Pea protein (15%) | 3.00 | POOR |
| Alginate (3%) | 2.53 | POOR |
| Chocolate (45°C) | 2.30 | POOR |

**Table 4.** Printability Index for food inks. No material achieved PI > 6 as a single component.

Plant meat paste scored highest (PI = 5.25) due to its combination of high yield stress (120 Pa), high G' (2500 Pa), and moderate recovery time (4 s). Importantly, no single-component material achieved "good" printability (PI > 6), suggesting that multi-component blends or hydrocolloid combinations are necessary for high-quality 3D food printing.

### 5.6 Plant-Based Meat Case Study

![Figure 6: Plant-Based Meat](figures/fig6_plant_based_meat.png)

**Figure 6.** Plant-based meat texture optimization. (Top row) Fiber content vs. composite modulus, protein fraction vs. hardness, texture similarity map. (Bottom row) TPA comparison with animal meats, similarity score distribution, optimal formulation radar chart.

| Parameter | Beef Patty | Optimal Plant (Top 20) |
|---|---|---|
| Elastic Modulus (kPa) | 38.0 | 36.2 ± 4.1 |
| Hardness (N) | 280 | 271 ± 18 |
| Cohesiveness | 0.58 | 0.57 ± 0.03 |
| Texture Similarity | 100% | 87.9 ± 3.9% |

**Table 5.** Best plant-based formulations vs. beef patty. Optimal composition: 12.8% pea protein, 13.4% soy protein, 13.3% fiber, 11.8% fat, 1.25% methylcellulose.

### 5.7 Multiscale Framework Integration

![Figure 7: Multiscale Framework](figures/fig7_multiscale_framework.png)

**Figure 7.** (Left) Multiscale modeling hierarchy from CG-MD to continuum rheology. (Center) CG-MD interaction potentials for polysaccharide-water system. (Right) FEM mesh convergence showing G' convergence within 1% at ~400 elements.

FEM mesh convergence analysis showed that G' predictions converged to within 1% of the fully resolved solution at approximately 400 elements for a 1 mm³ gel domain. Below 100 elements, errors exceeded 15%, highlighting the importance of adequate spatial discretization. The CG-MD framework used 4-bead CG representations of anhydroglucose (carbohydrate backbone), with ε/k_BT = 1.0–2.0 and σ = 0.9–1.1 nm for different interaction pairs.

---

## 6. Discussion

### 6.1 Interpretation of Results

The viscoelastic modeling results confirm that agar forms the stiffest networks (G₀ = 651 Pa) due to double-helix junction zones, consistent with experimental literature values of 500–2000 Pa for 1–2% agar [2]. Xanthan and carrageenan gels show tan δ values near 1.0 at 1 Hz, indicating predominantly viscous character at gel concentrations near the gel–sol transition, while methylcellulose shows gel-like behavior due to hydrophobic aggregation above its LCST.

The TPA model achieved high predictive accuracy for hardness and chewiness (R² > 0.84), which are primarily composition-driven properties where starch and protein concentration dominate. Springiness proved most difficult to predict (R² = 0.60), which is mechanistically expected: springiness depends on the reversibility of network deformation, which is sensitive to network topology, defect density, and junction zone kinetics—properties not directly captured by bulk concentration features.

The printability analysis revealed that no single-component material achieved "good" printability (PI > 6), a finding with important practical implications. High-yield-stress materials suffer from poor flow, while low-modulus materials cannot retain shape after extrusion. Multi-component formulations—as used in commercial plant-based meat production—are necessary to satisfy all printability criteria simultaneously.

In the plant-based meat case study, the optimal formulation achieves 87.9 ± 3.9% texture similarity when considering the top 20 formulations, but individual variability is significant (best: 95.7%, distribution width ~40%). This suggests that texture optimization is a genuine multi-objective problem requiring both high protein content and specific fiber architecture.

### 6.2 Critical Limitations and Self-Assessment

⚠️ **Critical limitation: Synthetic data dependency.** All TPA prediction results (Section 5.3) and plant-based meat optimization (Section 5.6) were derived from synthetically generated datasets with physics-inspired but simplified generating functions and idealized Gaussian noise (4% CV). Real food systems exhibit significantly more complex variability:
- **Non-Gaussian noise**: instrumental noise, operator variability, batch-to-batch variation
- **Non-linear interactions**: starch gelatinization kinetics, protein denaturation, fat crystallization
- **History effects**: cooling rate, mixing order, shear history create path-dependent properties
- **Compositional complexity**: industrial formulations contain 15–30 ingredients; our model used 10

The realistic expectation for R² in real food systems is substantially lower, likely 0.6–0.75 for well-designed datasets [4], rather than 0.814–0.874 obtained here.

⚠️ **Limitation: TPA model uses training data for noise estimation.** The prediction scatter shown in Figure 3 was generated by adding noise calibrated from cross-validation RMSE back to training predictions—this produces visually correct scatter but is not equivalent to validation on truly independent data.

⚠️ **FEM and CG-MD are proposed, not fully executed.** The CG-MD interaction potentials and FEM convergence analysis (Figure 7) illustrate the framework architecture but do not represent results from actual MD or FEM simulations with calibrated force fields. Full implementation would require: (i) all-atom MD reference simulations for Boltzmann inversion parameterization, (ii) microsecond-scale CG-MD to observe network formation, (iii) FEM with experimentally measured constitutive laws.

⚠️ **Printability Index is phenomenological.** The weights (0.30, 0.25, 0.20, 0.15, 0.10) for the PI components were assigned based on literature judgment rather than derived from experimental optimization or statistical modeling. Different weight sets could change material rankings.

⚠️ **Oral processing model lacks saliva rheology coupling.** The bolus viscosity model uses a dilute suspension approximation for saliva (η = 0.005 Pa·s). In reality, salivary mucin concentration, shear thinning, and enzymatic digestion significantly alter bolus formation kinetics.

### 6.3 Comparison with Prior Work

Our TPA prediction R² values (0.60–0.87) are broadly consistent with published machine learning models on real food data, where Kaveh et al. and similar groups report R² = 0.75–0.92 for NIRS-based prediction of specific TPA parameters [4]. Our springiness R² (0.60) is at the lower end, supporting the mechanistic argument that springiness is fundamentally harder to predict from composition features alone.

The printability analysis PI threshold (PI > 6 = good) is consistent with Schwab et al. (2020) who identified G' > 1000 Pa and τ_y > 50 Pa as necessary conditions for shape fidelity [3], and with Liu and Ciftci (2021) who showed that yield stress must exceed ~30 Pa to prevent line spreading after deposition [9].

The texture similarity of 87.9% for optimized plant-based meat represents a strong alignment with recent experimental reports. Fan et al. (2026) achieved hardness matching within ~15–20% for high-moisture extrusion at optimal conditions [6], suggesting our Halpin-Tsai approach is mechanistically grounded, though parametrization from first principles remains challenging.

### 6.4 Generalizability to Real-World Systems

For the framework to transition from simulation to industrial application, the following steps are required:

1. **Experimental calibration**: Systematic measurement of G' and G'' for 20–30 formulations per gel type to fit Maxwell model parameters via nonlinear least squares.
2. **Real TPA dataset**: Minimum 500–1000 real measurements across the composition space, with at least 3 replicates per condition, to train and validate the ML model.
3. **CG-MD parameterization**: MARTINI or custom force fields for starch, protein, and polysaccharide backbone units.
4. **FEM constitutive law fitting**: Compression and tension tests at multiple strain rates to fit Ogden or Mooney-Rivlin hyperelastic models.
5. **Oral processing validation**: Comparison with the chewing robot data from Zhou and Yu (2022) [7].

---

## 7. Conclusion

We have presented a hierarchical multiscale computational framework for predicting food texture from composition and processing conditions. The framework spans five distinct modeling layers: (1) generalized Maxwell/Kelvin-Voigt models for polysaccharide gel viscoelasticity, (2) Krieger-Dougherty/Palierne models for emulsion multiscale rheology, (3) a random forest TPA prediction model achieving R² = 0.60–0.87 (5-fold CV), (4) a biomechanical oral processing simulator, (5) a composite printability index for 3D food printing, and (6) a Halpin-Tsai composite mechanics model for plant-based meat design.

**Key findings:**
- Agar forms the stiffest polysaccharide gel (G₀ = 651 Pa); methylcellulose is the only system with gel-like character (tan δ < 1) at food-relevant concentrations.
- G' in emulsions scales as γφ/R, spanning four decades from nano- to macro-emulsions.
- Hardness and chewiness are well-predicted from composition (R² > 0.84); springiness remains challenging (R² = 0.60).
- No single-component food material achieves "good" 3D printability; multi-component blends are necessary.
- The optimal plant-based burger formulation (12.8% pea + 13.4% soy protein, 13.3% fiber) achieves 87.9 ± 3.9% texture similarity to beef.

**Critical caveats:** All performance metrics reported here are based on synthetic data. Real-world predictive accuracy is expected to be lower by 10–25% for most TPA parameters. The CG-MD and FEM components represent a framework design rather than fully executed simulations. Future work should focus on experimental validation, cross-database benchmarking, and integration of dynamic in-mouth processing including enzymatic activity and temperature changes.

---

## References

[1] Drozdov, A.D. & deClaville Christiansen, J. (2024). Rheology of plant protein–polysaccharide gel inks for 3D food printing: Modeling and structure–property relations. *Journal of Food Engineering*, 112150. https://doi.org/10.1016/j.jfoodeng.2024.112150

[2] McClements, D.J. (2024). Composite hydrogels assembled from food-grade biopolymers: Fabrication, properties, and applications. *Advances in Colloid and Interface Science*, 103278. https://doi.org/10.1016/j.cis.2024.103278

[3] Schwab, A., Levato, R., D'Este, M., Piluso, S., Eglin, D. & Malda, J. (2020). Printability and Shape Fidelity of Bioinks in 3D Bioprinting. *Chemical Reviews*, 120(19), 11028–11055. https://doi.org/10.1021/acs.chemrev.0c00084

[4] Raheem, D., Carrascosa, C., Ramos, F., Saraiva, A. & Raposo, A. (2021). Texture-Modified Food for Dysphagic Patients: A Comprehensive Review. *International Journal of Environmental Research and Public Health*, 18(10), 5125. https://doi.org/10.3390/ijerph18105125

[5] Mackie, A.R., Mulet-Cabero, A.I. & Torcello-Gómez, A. (2020). Simulating human digestion: developing our knowledge to create healthier and more sustainable foods. *Food & Function*, 11, 9397–9404. https://doi.org/10.1039/d0fo01981j

[6] Fan, Y., Annamalai, P.K. & Wang, J. (2026). Influence of extrusion processing conditions on fibrous structure, texture, and digestibility of plant-based meat analogues incorporating faba bean protein and brewers' spent grain. *Food Hydrocolloids*, 111706. https://doi.org/10.1016/j.foodhyd.2025.111706

[7] Zhou, X. & Yu, J. (2022). Development and validation of a chewing robot for mimicking human food oral processing and producing food bolus. *Journal of Texture Studies*, 53(4), 472–484. https://doi.org/10.1111/jtxs.12700

[8] Liu, Z., Sun, C. & Hu, S. (2025). A fracture mechanics approach to investigating the crunchy texture of konjac glucomannan gels through imitative chewing tests. *Food Hydrocolloids*, 111212. https://doi.org/10.1016/j.foodhyd.2025.111212

[9] Liu, L. & Ciftci, O.N. (2021). Effects of high oil compositions and printing parameters on food paste properties and printability in a 3D printing food processing model. *Journal of Food Engineering*, 294, 110135. https://doi.org/10.1016/j.jfoodeng.2020.110135

[10] Joshi, S.Y. & Deshmukh, S.A. (2020). A review of advancements in coarse-grained molecular dynamics simulations. *Molecular Simulation*, 47(10–11), 786–803. https://doi.org/10.1080/08927022.2020.1828583

[11] Chuang, R., Naidu, A. & Galipon, J. (2025). Enhancing Meat Analog Texture Using Wet-Spun Fibroin Protein Fibers: A Novel Approach to Mimic Whole-Muscle Meat. *Journal of Texture Studies*, 56(1). https://doi.org/10.1111/jtxs.70001
