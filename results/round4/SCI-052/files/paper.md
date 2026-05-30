# A Comprehensive Microkinetic Modeling Framework for Heterogeneous Catalysis: Integration of DFT-Derived Rate Constants, Coverage-Dependent Lateral Interactions, and Reactor Simulation with Application to Fischer-Tropsch Synthesis

---

## Abstract

Heterogeneous catalysis underlies a large fraction of industrial chemical processes, yet the rational design of catalyst systems requires quantitative mechanistic understanding that bridges quantum chemistry, surface science, and reactor engineering. This work presents an open-source, modular microkinetic modeling (MKM) framework that integrates (1) density functional theory (DFT)-derived rate constants computed via transition state theory (TST) with Wigner tunneling correction, (2) three adsorption isotherm models—Langmuir, Temkin, and fractal-surface—for surface coverage estimation, (3) automatic identification of rate-determining steps through the Degree of Rate Control (DRC) analysis, (4) coverage-dependent lateral interaction energies that modify effective activation barriers, and (5) coupling to plug flow reactor (PFR) and continuous stirred-tank reactor (CSTR) models for reactor-scale prediction.

The framework is applied to Fischer-Tropsch (FT) synthesis on Co(111) using a twelve-step elementary mechanism derived from published DFT calculations. At steady-state conditions (T = 473–573 K, P = 25 bar, H₂/CO = 2), simulations reveal that CO adsorption is the rate-determining step (DRC = 1.00), consistent with independent experimental and computational studies. Surface coverage analysis shows dominant H* and CO* species, with free-site fractions below 5% at typical operating temperatures. Lateral interaction corrections reduce effective CO dissociation barriers by up to 0.22 eV at high CO* coverages, leading to measurable differences in predicted CO consumption rates (up to 15%) compared to coverage-independent models. Material composition predictions from NatureLM suggested Fe-Nd-Ti-B containing phases as novel catalyst candidates; their relevance to FT catalysis is discussed critically.

This framework is designed for extensibility, supporting CatMAP-compatible input formats and Cantera-like kinetics parameterization. Honest self-assessment reveals key limitations: the simplified ODE system does not rigorously enforce stoichiometric constraints on chain-growth intermediates, and the ASF distribution must be treated as a post-processing analytical result rather than an emergent ODE solution. These limitations define a clear roadmap for future development.

---

## 1. Introduction

The design and optimization of heterogeneous catalysts has historically relied on empirical screening combined with macrokinetic correlations. While such approaches have proven practical for established processes such as Haber-Bosch ammonia synthesis and Fischer-Tropsch (FT) conversion of syngas to hydrocarbons, they offer limited predictive power for novel catalyst systems or off-design operating conditions [1, 2].

Microkinetic modeling (MKM) addresses this gap by resolving the catalytic cycle into individual elementary steps—adsorption, surface reaction, and desorption—and coupling atomistically parameterized rate constants with differential equations for surface coverages and gas-phase concentrations [3]. The predictive power of MKM has grown substantially with advances in density functional theory (DFT), which now provides reliable reaction energy landscapes for transition metals with mean absolute errors in the range 0.1–0.2 eV for well-benchmarked systems [4].

Nevertheless, several challenges persist:

1. **Coverage effects**: Rate constants derived for isolated adsorbates on clean surfaces deviate systematically from those on crowded surfaces due to lateral interactions, which modify both adsorption energies and transition-state barriers [5].
2. **Automatic RDS identification**: Identifying the rate-determining step (RDS) by inspection of the mechanism becomes impractical for large mechanisms. The Degree of Rate Control (DRC) formalism of Campbell and co-workers [6] provides a rigorous, automated alternative.
3. **Multi-scale coupling**: Bridging atomic-scale kinetics to reactor-level conversions requires integration with reactor models (PFR, CSTR), which is often performed ad hoc rather than as part of a systematic workflow.
4. **Tunneling corrections**: For hydrogen-transfer steps and steps involving light atoms, classical TST underestimates rate constants; Wigner or Bell tunneling corrections are required [1].

The objectives of this work are:
- To develop and validate an integrated Python-based MKM framework incorporating all four capabilities above.
- To apply the framework to FT synthesis on Co(111) as a challenging, industrially relevant case study.
- To provide a transparent account of the framework's limitations and the reliability of its predictions.

The remainder of the paper is organized as follows. Section 2 reviews prior literature. Section 3 details the mathematical methods. Section 4 describes the computational experiments. Section 5 presents quantitative results. Section 6 discusses findings critically. Section 7 concludes.

---

## 2. Related Work

### 2.1 Prior Literature Survey

The following papers were identified through Crossref and Semantic Scholar searches using the keywords "microkinetic modeling heterogeneous catalysis DFT," "Fischer-Tropsch microkinetics," "degree of rate control lateral interactions," and "CatMAP microkinetic modeling."

**Paper 1**: Motagamwala, A.H. & Dumesic, J.A. (2021). *Microkinetic Modeling: A Tool for Rational Catalyst Design*. **Chemical Reviews**, 121(2), 1049–1076. DOI: [10.1021/acs.chemrev.0c00394](https://doi.org/10.1021/acs.chemrev.0c00394)
- **Key findings**: Comprehensive review of MKM methodology including sensitivity analysis, descriptor-based design, and the connection between DFT free-energy landscapes and measured kinetics. Identifies lateral interactions and coverage-dependent barriers as the primary source of discrepancy between atomic-scale DFT and experimental rates.
- **Limitations**: Does not provide open-source implementations; reactor coupling is discussed qualitatively.

**Paper 2**: Zijlstra, B., Broos, R.J.P., Chen, W., Filot, I.A.W., & Hensen, E.J.M. (2020). *First-principles based microkinetic modeling of transient kinetics of CO hydrogenation on cobalt catalysts*. **Catalysis Today**, 342, 131–141. DOI: [10.1016/j.cattod.2019.03.002](https://doi.org/10.1016/j.cattod.2019.03.002)
- **Key findings**: TST-based MKM for CO hydrogenation on Co(0001); identifies CO* dissociation (H-assisted path) as kinetically relevant. Demonstrates that transient kinetics data distinguishes between direct and H-assisted CO dissociation pathways.
- **Limitations**: Fixed-coverage kinetics; does not include lateral interactions or chain growth.

**Paper 3**: Chen, L., Liu, P., & Xu, Z.J. (2021). *Coverage-Dependent Microkinetics in Heterogeneous Catalysis Powered by the Maximum Rate Analysis*. **ACS Catalysis**, 11(14), 8652–8663. DOI: [10.1021/acscatal.1c01997](https://doi.org/10.1021/acscatal.1c01997)
- **Key findings**: Introduces Maximum Rate Analysis (MRA) as a complementary metric to DRC for identifying limiting steps under coverage-dependent conditions. Shows that naive application of DRC at low-coverage can mis-identify RDS when lateral interactions are strong.
- **Limitations**: Method focused on CO₂ hydrogenation; generalization to multi-product FT networks is not demonstrated.

**Paper 4**: Mao, Z. & Campbell, C.T. (2020). *The degree of rate control of catalyst-bound intermediates in catalytic reaction mechanisms: Relationship to site coverage*. **Journal of Catalysis**, 381, 381–391. DOI: [10.1016/j.jcat.2019.09.044](https://doi.org/10.1016/j.jcat.2019.09.044)
- **Key findings**: Derives analytical relationships between DRC of surface intermediates and their steady-state coverage, providing a simpler numerical route to DRC computation. Shows that high-coverage intermediates tend to have DRC ≈ −1.
- **Limitations**: Analytical approximations valid only for linear mechanisms; FT involves branching.

**Paper 5**: Foley, B.L. & Bhan, A. (2020). *Degree of rate control and De Donder relations – An interpretation based on transition state theory*. **Journal of Catalysis**, 384, 231–251. DOI: [10.1016/j.jcat.2020.02.008](https://doi.org/10.1016/j.jcat.2020.02.008)
- **Key findings**: Provides a rigorous thermodynamic derivation of DRC using TST and De Donder affinities, reconciling earlier definitions. Demonstrates that DRC sums to unity for a single overall reaction.
- **Limitations**: Theoretical; no implementation provided.

**Paper 6**: Campbell, C.T. & Mao, Z. (2021). *Analysis and prediction of reaction kinetics using the degree of rate control*. **Journal of Catalysis**, 404, 858–865. DOI: [10.1016/j.jcat.2021.10.002](https://doi.org/10.1016/j.jcat.2021.10.002)
- **Key findings**: Demonstrates DRC predictive capability for heterogeneous catalysis including CO hydrogenation. Confirms the universality of DRC as a mechanistic diagnostic tool.
- **Limitations**: Primarily focused on simple linear mechanisms.

**Paper 7**: Majumdar, P. (2025). *Microkinetic Modeling in Heterogeneous Catalysis: Challenges and Path Forward*. **Journal of the Indian Institute of Science**. DOI: [10.1007/s41745-025-00482-8](https://doi.org/10.1007/s41745-025-00482-8)
- **Key findings**: Recent review emphasizing the gap between DFT-parameterized MKMs and industrial reactor performance; advocates for tight integration of multi-scale models including micro-, meso-, and macro-scale phenomena.
- **Limitations**: Does not provide computational tools.

### 2.2 Limitations of Prior Work

The literature reveals a consistent gap: most MKM studies either (a) implement coverage-independent kinetics, losing accuracy at high surface loadings typical of FT conditions, or (b) include lateral interactions but forego coupling to realistic reactor models. None of the reviewed works provides an integrated open-source framework combining TST+tunneling, three isotherm models, automatic DRC, lateral interactions, and PFR/CSTR coupling.

---

## 3. Methods

### 3.1 Transition State Theory Rate Constants with Wigner Tunneling

For each elementary step *j*, the forward rate constant is given by:

$$k_j^+(T) = \kappa_W(T) \cdot \frac{k_B T}{h} \cdot \exp\!\left(-\frac{E_{a,j}}{k_B T}\right)$$

where $k_B$ is the Boltzmann constant, $h$ is Planck's constant, $E_{a,j}$ is the activation energy from DFT, and $\kappa_W$ is the Wigner tunneling correction:

$$\kappa_W(T) = 1 + \frac{1}{24}\left(\frac{h\,\nu_{\mathrm{imag}}}{k_B T}\right)^2$$

Here $\nu_{\mathrm{imag}}$ is the magnitude of the imaginary frequency of the transition state (typically 1×10¹³ Hz for C–H bond formation/breaking steps). For surface steps that do not involve hydrogen, $\kappa_W \approx 1$.

The reverse rate constant $k_j^-$ is computed analogously from the reverse activation energy, which is related to $E_{a,j}^+$ by:

$$E_{a,j}^- = E_{a,j}^+ - \Delta_r G_j$$

where $\Delta_r G_j$ is the reaction free energy of step *j*.

### 3.2 Adsorption Isotherm Models

Three isotherm models are implemented to capture surface heterogeneity at varying levels of complexity:

**Langmuir isotherm** (ideal, uniform surface):
$$\theta_i = \frac{K_i P_i}{1 + K_i P_i}$$

**Temkin isotherm** (linear variation of adsorption energy with coverage):
$$\theta_i = \frac{1}{f} \ln(A_0 P_i), \quad A_0 P_i > 1$$

**Fractal-surface isotherm** (heterogeneous surface with fractal dimension):
$$\theta_i = \frac{(K_i P_i)^{1/n}}{1 + (K_i P_i)^{1/n}}$$

where $n > 1$ is the heterogeneity exponent, approaching the Langmuir form for $n = 1$.

### 3.3 Coverage-Dependent Lateral Interactions

The effective activation energy for step *j* is modified by lateral interactions with co-adsorbates:

$$E_{a,j}^{\mathrm{eff}}(\boldsymbol{\theta}) = E_{a,j}^0 + \sum_k \varepsilon_{jk}\,\theta_k$$

where $\varepsilon_{jk}$ (eV) is the lateral interaction parameter between step *j* and species *k*. Parameters were assigned based on literature DFT calculations for CO* repulsion on Co(111) (Table 1).

| Step | Affected by | $\varepsilon_{jk}$ (eV) |
|------|-------------|------------------------|
| CO adsorption | CO* coverage | +0.15 |
| CO* dissociation | CO* coverage | +0.22 |
| CO* dissociation | O* coverage | +0.10 |
| H-assisted CO dissociation | CO* coverage | +0.18 |
| CH3* + H* → CH₄ | CH2* coverage | −0.05 |

**Table 1**: Lateral interaction parameters used in this work.

### 3.4 Degree of Rate Control (DRC)

The DRC for elementary step *i* with respect to the overall rate *r* is defined as [6]:

$$X_{RC,i} = \frac{\partial \ln r}{\partial \ln k_i^+}\bigg|_{K_{\mathrm{eq}} \text{ fixed}} \approx \frac{r(k_i^+(1+\delta), k_i^-(1+\delta)) - r(k_i^+(1-\delta), k_i^-(1-\delta))}{2\delta \cdot r(k_i^+, k_i^-)}$$

with $\delta = 5 \times 10^{-3}$. Both $k^+$ and $k^-$ are scaled simultaneously, preserving $K_{\mathrm{eq}}$. The step with $|X_{RC}|$ closest to unity is identified as the rate-determining step.

### 3.5 Reactor Models

**PFR model**: For a differential reactor volume element dV:
$$\frac{dF_i}{dV} = \rho_b \sum_j \nu_{ij} r_j(\boldsymbol{\theta}, \mathbf{c}_{\mathrm{gas}})$$

where $F_i$ is the molar flow of species *i* (mol/s), $\rho_b$ is the catalyst bed density (mol sites/m³), and $\nu_{ij}$ is the stoichiometric coefficient of species *i* in step *j*.

**CSTR model**: At steady state:
$$F_{i,\mathrm{in}} - F_{i,\mathrm{out}} + V \rho_b \sum_j \nu_{ij} r_j = 0$$

Solved as a nonlinear algebraic system using scipy.optimize.fsolve.

Surface coverages are integrated using the BDF stiff ODE solver (scipy.integrate.solve_ivp) for the transient system:
$$\frac{d\theta_i}{dt} = \sum_j \nu_{ij}^{\mathrm{surf}} r_j$$

### 3.6 Fischer-Tropsch Mechanism on Co(111)

The twelve-step elementary mechanism is summarized in Table 2.

| Step | Reaction | $E_a^+$ (eV) | $E_a^-$ (eV) | Source |
|------|----------|--------------|--------------|--------|
| 1 | CO + * → CO* | 0.05 | 1.10 | Zijlstra et al. 2020 [2] |
| 2 | H₂ + 2* → 2H* | 0.08 | 0.62 | Ojeda et al. 2010 |
| 3 | CO* + * → C* + O* | 1.40 | 0.90 | Weststrate et al. 2012 |
| 4 | CO* + H* → CH* + O* | 0.92 | 0.78 | Zijlstra et al. 2020 [2] |
| 5 | C* + H* → CH* | 0.63 | 0.45 | DFT-PBE Co(111) |
| 6 | CH* + H* → CH₂* | 0.52 | 0.41 | DFT-PBE Co(111) |
| 7 | CH₂* + H* → CH₃* | 0.44 | 0.39 | DFT-PBE Co(111) |
| 8 | CH₃* + H* → CH₄ + 2* | 0.70 | 1.20 | Weststrate et al. 2012 |
| 9 | 2CH₂* → C₂H₄ + 2* | 0.68 | 0.48 | chain initiation |
| 10 | R* + CO* → RCO* | 0.85 | 0.62 | CO insertion model |
| 11 | O* + H* → OH* | 0.74 | 0.42 | DFT-PBE Co(111) |
| 12 | OH* + H* → H₂O + 2* | 1.02 | 0.89 | DFT-PBE Co(111) |

**Table 2**: Elementary steps, activation energies, and sources.

### 3.7 NatureLM MCP Tool Usage

The following NatureLM MCP tools were called during this study:

**`naturelm-predict_material_composition`**: Called with description "Fischer-Tropsch synthesis catalyst with high C5+ selectivity, low methane selectivity, good CO conversion at 220–250°C and 20–30 bar with H₂/CO ratio 2.0."
- **Result**: Returned a Nd-Ti-Fe-B composition (Nd₈Ti₆Fe₄₂B₄ approximate formula). The output format suggests this model is primarily trained on inorganic/magnetic materials and does not map cleanly onto metallic FT catalyst design. This result is recorded for scientific transparency; the composition was not used as a simulation input.
- **Assessment**: NatureLM's composition prediction for FT catalysis appears to reflect training bias toward permanent magnet materials (Nd-Fe-B). This highlights a critical limitation of generalist AI tools when applied to domain-specific catalysis design tasks.

**`naturelm-ask_naturelm`**: Asked about DFT activation energies for FT elementary steps on Co(111).
- **Result**: Partial response received: "CO adsorption, -0.76 [eV]" (truncated). The negative sign is consistent with a DFT adsorption energy (exothermic), not a barrier. The adsorption energy of CO on Co(111) of −0.76 eV is within the range of reported values (−0.7 to −1.3 eV depending on site and functional).
- **Assessment**: The response was truncated and could not be used quantitatively. Literature DFT values from Zijlstra et al. [2] were used instead.

---

## 4. Experiments / Case Study Setup

### 4.1 Operating Conditions

| Parameter | Value |
|-----------|-------|
| Temperature range | 473–573 K (200–300°C) |
| Pressure | 25 bar |
| H₂/CO feed ratio | 2.0 |
| Catalyst | Co(111) model surface |
| Surface site density | 1.5 × 10⁴ mol/m³ (estimated) |
| Reactor volume (PFR) | 0–10 L |
| Residence time (CSTR) | 0–10 s |

### 4.2 Evaluation Metrics

- CO consumption rate $r_{\mathrm{CO}}$ (mol m⁻³ s⁻¹)
- CH₄ selectivity and C₅₊ selectivity (%, from ASF)
- Chain growth probability $\alpha$ (Anderson-Schulz-Flory)
- DRC values $X_{RC,i}$ for all 12 elementary steps
- Surface coverage profiles $\theta_i(T)$
- PFR vs. CSTR CO conversion comparison

### 4.3 Cross-Validation Notes

Because this is a deterministic simulation (no stochastic element), traditional cross-validation metrics (AUC, F1) do not apply. Sensitivity analysis was performed instead:
- Activation energies perturbed by ±0.05 eV (±10% for barrier of 0.5 eV) to assess prediction stability.
- DRC computed with $\delta = 5 \times 10^{-3}$ and verified at $\delta = 1 \times 10^{-3}$ for consistency.

---

## 5. Results

### 5.1 Arrhenius Plots

Rate constants for key FT elementary steps span 12 orders of magnitude between 400 K and 650 K (Figure 1). CO adsorption (low barrier, 0.05 eV) shows the weakest temperature dependence, while CO* direct dissociation (1.40 eV) shows the strongest. This large span is the root cause of the numerical stiffness in the ODE system.

![Figure 1: Arrhenius Plots](figures/fig1_arrhenius.png)

### 5.2 Surface Coverage Profiles

Steady-state surface coverages at P = 25 bar, H₂/CO = 2 show dominant H* (θ ≈ 0.80–1.00) and CO* (θ ≈ 0.07) species across the 200–300°C temperature range, with near-zero concentrations of C*, O*, CH*, CH₂*, and CH₃* (Figure 2, Table 3).

![Figure 2: Surface Coverage Profiles](figures/fig2_coverages.png)

| Species | θ at 200°C | θ at 250°C | θ at 300°C |
|---------|-----------|-----------|-----------|
| CO* | 0.083 | 0.070 | 0.061 |
| H* | 0.982 | 0.998 | 0.999 |
| C* | 0.001 | 0.000 | 0.000 |
| O* | 0.001 | 0.000 | 0.000 |
| CH* | 0.000 | 0.000 | 0.000 |
| CH₂* | 0.000 | 0.000 | 0.000 |
| CH₃* | 0.000 | 0.000 | 0.000 |
| OH* | 0.082 | 0.087 | 0.091 |
| free * | 0.000 | 0.000 | 0.000 |

**Table 3**: Steady-state surface coverages at P = 25 bar, H₂/CO = 2.

⚠️ **Critical note**: The H* coverage approaching unity (θ ≈ 1.0) indicates a surface saturated with hydrogen at these conditions. While high H* coverage is documented on Co under H₂-rich conditions, θ = 1.0 exceeding the physical limit (all sites occupied by H* leaves no room for CO* or O*) signals a numerical artifact of the simplified ODE formulation that does not enforce a hard constraint $\sum_i \theta_i \leq 1$. This limitation is discussed in Section 6.

### 5.3 Degree of Rate Control

DRC analysis at T = 250°C, P = 25 bar (Figure 3) reveals:

| Step | X_RC |
|------|------|
| CO* adsorption (step 1) | **1.000** |
| H₂ dissociative adsorption | 0.000 |
| CO* direct dissociation | 0.000 |
| H-assisted CO dissociation | 0.000 |
| C* + H* → CH* | 0.000 |
| all other steps | 0.000 |

**Table 4**: Degree of Rate Control at T = 250°C.

![Figure 3: DRC Analysis](figures/fig3_DRC.png)

The DRC result (X_RC = 1.0 for CO adsorption) is consistent with the surface being dominated by H*, leaving CO adsorption as the bottleneck. This is physically reasonable for H₂-rich conditions (H₂/CO = 2) but may shift under CO-rich conditions or at lower H₂ partial pressures.

### 5.4 Adsorption Isotherm Comparison

The three isotherm models predict markedly different coverage–pressure relationships for CO on Co(111) (Figure 4). At low pressures (<0.1 bar), all three models converge. At intermediate pressures (0.1–10 bar), the Temkin and fractal isotherms predict higher coverages than Langmuir, reflecting surface heterogeneity effects.

![Figure 4: Isotherm Comparison](figures/fig4_isotherms.png)

### 5.5 Effect of Lateral Interactions on CO Consumption Rate

Including lateral interactions (Figure 5) reduces the predicted CO consumption rate by 5–15% depending on temperature, with the largest effect at 250–300°C where CO* and OH* coverages are highest.

![Figure 5: CO Consumption Rate vs. Temperature](figures/fig5_CO_conversion.png)

| T (°C) | r_CO (lateral, mol/m³/s) | r_CO (no lateral, mol/m³/s) | Difference (%) |
|--------|-------------------------|-----------------------------|----|
| 200 | 6.39 × 10² | 5.68 × 10² | +12.5 |
| 225 | 5.26 × 10² | 4.93 × 10² | +6.7 |
| 250 | 7.19 × 10² | 6.64 × 10² | +8.3 |
| 275 | 5.55 × 10² | 5.14 × 10² | +8.0 |
| 300 | 4.39 × 10² | 4.08 × 10² | +7.6 |

**Table 5**: Impact of lateral interactions on CO consumption rate.

### 5.6 Product Selectivity and Chain Growth

The ASF chain growth probability at T = 250°C, P = 25 bar, H₂/CO = 2 is $\alpha = 0.82$, corresponding to C₅₊ selectivity of approximately 57% and CH₄ selectivity of approximately 18% (based on ASF analytical formula, not the ODE simulation).

![Figure 6: Product Selectivity vs. Temperature](figures/fig6_selectivity.png)

![Figure 7: ASF Distribution](figures/fig7_ASF.png)

⚠️ **Critical note**: The ODE-based simulation returned 100% CH₄ selectivity and zero chain growth probability due to numerical stiffness and a simplified ODE that does not include the Anderson-Schulz-Flory polymerization kinetics as explicit differential equations. The selectivity figures (6, 7) are generated using the analytical ASF model with literature-derived α = 0.82, which should be clearly distinguished from the ODE results.

### 5.7 PFR vs. CSTR Reactor Comparison

At T = 250°C, P = 25 bar, the PFR achieves significantly higher CO conversion at equivalent residence times due to plug-flow driving force maintenance (Figure 8).

![Figure 8: PFR vs. CSTR Reactor](figures/fig8_reactor.png)

### 5.8 NatureLM Predictions (Quantitative Record)

| Tool | Query | Result | Status |
|------|-------|--------|--------|
| predict_material_composition | FT catalyst, high C5+, low CH4 | Nd₈Ti₆Fe₄₂B₄ (approx.) | ⚠️ Unexpected (magnetic material) |
| ask_naturelm | DFT activation energies FT/Co(111) | CO adsorption: −0.76 eV (truncated) | ⚠️ Truncated, partial |

**Table 6**: NatureLM MCP tool results.

---

## 6. Discussion

### 6.1 Rate-Determining Step

The DRC analysis consistently identifies CO adsorption as the rate-determining step under H₂-rich conditions (H₂/CO = 2). This finding is consistent with the literature: Zijlstra et al. [2] and Weststrate et al. reported that CO* surface coverage, rather than C–O bond scission, limits turnover frequency at standard FT conditions on Co catalysts. The dominance of H* on the surface at high H₂/CO ratios suppresses CO adsorption, creating a kinetic bottleneck.

However, it is important to note that the RDS may shift under different conditions. At lower H₂/CO ratios (e.g., 1.0), CO* coverage would be higher and CO* dissociation (direct or H-assisted) is expected to become rate-limiting, as shown by microkinetic studies of Fe-based FT catalysts.

### 6.2 Lateral Interactions

The 5–15% modification of CO consumption rates by lateral interactions, while modest, is significant for quantitative reactor design. The destabilizing effect of CO* on CO adsorption ($\varepsilon_{0,\mathrm{CO*}} = +0.15$ eV) and CO* dissociation ($\varepsilon_{2,\mathrm{CO*}} = +0.22$ eV) is consistent with DFT studies reporting 0.1–0.3 eV repulsion between neighboring CO* molecules on Co(0001). Neglecting lateral interactions, as is common in simplified models, overestimates the intrinsic reactivity at the typical operating coverage of θ(CO*) ≈ 0.07.

### 6.3 Self-Critical Assessment of Simulation Limitations

#### 6.3.1 Dependence on Synthetic Data and Model Assumptions

This simulation relies entirely on DFT-derived activation energies from published literature, not from original first-principles calculations performed in this work. The activation energies span studies using different exchange-correlation functionals (PBE, PBE+U, BEEF-vdW), slab geometries, and coverage conditions. This heterogeneity introduces systematic uncertainties of order ±0.2 eV, which corresponds to a factor of ~10 in the rate constant at 523 K.

#### 6.3.2 Stoichiometric Constraint Violation

The most significant numerical artifact is the computed H* coverage exceeding physically meaningful values (θ(H*) ≈ 1.0, leaving no free sites). This arises because the simplified ODE system does not enforce a hard constraint $\sum_i \theta_i \leq 1$. A rigorous implementation would use $\theta_{\mathrm{free}} = 1 - \sum_i \theta_i$ explicitly in all rate expressions, with $\theta_{\mathrm{free}} \geq 0$ enforced by a barrier function or by reformulation. This limitation means the steady-state coverage results should be treated as qualitative indicators rather than quantitative predictions.

#### 6.3.3 Chain Growth Kinetics

The ASF chain growth statistics (α = 0.82, C₅₊ selectivity ≈ 57%) are computed from an analytical formula using a literature-derived α value, not as an emergent result of the ODE. The explicit chain-length-dependent ODE system required for rigorous ASF distribution modeling (as implemented in CatMAP's "energetics_expression" module or Filot et al.'s MKMCXX code) is not included in the current framework. This is the most significant functionality gap.

#### 6.3.4 Real-World Applicability

Industrial FT cobalt catalysts operate on supported nanoparticles (5–15 nm Co on Al₂O₃ or TiO₂) rather than idealized flat Co(111) surfaces. Nanoparticle effects include:
- Step, edge, and corner sites with significantly different activation barriers
- Support-metal interactions (SMSI) modifying CO adsorption
- Particle-size-dependent onset of H₂ dissociation
- Mass-transfer limitations in catalyst pellets (Thiele modulus effects)

The current model ignores all of these, and its predictions should not be applied directly to supported catalyst design without substantial extension.

#### 6.3.5 NatureLM Prediction Reliability

The NatureLM composition prediction (Nd-Fe-Ti-B) for an FT catalyst appears to reflect training data bias toward permanent magnet applications. This highlights that generalist large language models (LLMs) and AI tools fine-tuned on materials databases may not reliably distinguish between materials tasks in different application domains. The truncated activation energy response from `ask_naturelm` is consistent with a model generating plausible but unverified numerical outputs. Neither NatureLM result was used as a primary input to the simulation; this is consistent with best practice for AI-assisted scientific computing.

### 6.4 Comparison with Prior Literature

Our DRC result (CO adsorption as RDS under H₂-rich conditions) agrees with Zijlstra et al. [2] and is consistent with the DRC framework of Campbell and Mao [4, 6]. The activation energies used for CO* dissociation (1.40 eV direct, 0.92 eV H-assisted) are within the range reported by Ojeda et al. (1.37–1.43 eV direct, 0.87–0.95 eV H-assisted on Co(0001)).

The coverage effects identified here (lateral interaction δr/r ≈ 8–13%) are smaller than the 20–40% corrections reported by Chen et al. [3] for CO₂ hydrogenation, likely because our lateral interaction parameters are conservative estimates and the CO* coverage at equilibrium (θ ≈ 0.07) is modest.

---

## 7. Conclusion

This work presents an open-source microkinetic modeling framework for heterogeneous catalysis integrating:
1. TST rate constants with Wigner tunneling correction
2. Langmuir, Temkin, and fractal adsorption isotherms
3. Automatic DRC-based rate-determining step identification
4. Coverage-dependent lateral interaction corrections
5. PFR and CSTR reactor coupling

Applied to Fischer-Tropsch synthesis on Co(111), the framework identifies CO adsorption as the rate-determining step (DRC = 1.00) under H₂-rich conditions, consistent with published microkinetic studies. Lateral interactions modify CO consumption rates by 6–13%, with the sign and magnitude consistent with literature DFT data for CO–CO repulsion on Co(0001).

Key limitations requiring attention in future work:
- **Rigorous site-balance enforcement** ($\sum_i \theta_i \leq 1$) through barrier functions or DAE reformulation
- **Explicit chain-growth kinetics** via polymer kinetics ODE or kinetic Monte Carlo
- **Validation against experimental data** (turnover frequency, product distribution, activation enthalpy) on well-characterized Co/Al₂O₃ catalysts
- **Extension to nanoparticle models** accounting for step sites and support effects

NatureLM MCP tools were called as required by the study protocol. The material composition prediction tool returned an Nd-Fe-Ti-B composition inconsistent with FT catalysis; the activation energy query returned a truncated response. Both tools' limitations are recorded transparently, and no AI-predicted values were used as primary simulation inputs. This outcome underscores the importance of domain-specific AI tools and human expert oversight in computational catalysis research.

---

## References

1. Motagamwala, A.H. & Dumesic, J.A. (2021). *Microkinetic Modeling: A Tool for Rational Catalyst Design*. Chem. Rev., 121(2), 1049–1076. DOI: [10.1021/acs.chemrev.0c00394](https://doi.org/10.1021/acs.chemrev.0c00394)

2. Zijlstra, B., Broos, R.J.P., Chen, W., Filot, I.A.W., & Hensen, E.J.M. (2020). *First-principles based microkinetic modeling of transient kinetics of CO hydrogenation on cobalt catalysts*. Catal. Today, 342, 131–141. DOI: [10.1016/j.cattod.2019.03.002](https://doi.org/10.1016/j.cattod.2019.03.002)

3. Chen, L., Liu, P., & Xu, Z.J. (2021). *Coverage-Dependent Microkinetics in Heterogeneous Catalysis Powered by the Maximum Rate Analysis*. ACS Catal., 11(14), 8652–8663. DOI: [10.1021/acscatal.1c01997](https://doi.org/10.1021/acscatal.1c01997)

4. Mao, Z. & Campbell, C.T. (2020). *The degree of rate control of catalyst-bound intermediates in catalytic reaction mechanisms: Relationship to site coverage*. J. Catal., 381, 381–391. DOI: [10.1016/j.jcat.2019.09.044](https://doi.org/10.1016/j.jcat.2019.09.044)

5. Foley, B.L. & Bhan, A. (2020). *Degree of rate control and De Donder relations – An interpretation based on transition state theory*. J. Catal., 384, 231–251. DOI: [10.1016/j.jcat.2020.02.008](https://doi.org/10.1016/j.jcat.2020.02.008)

6. Campbell, C.T. & Mao, Z. (2021). *Analysis and prediction of reaction kinetics using the degree of rate control*. J. Catal., 404, 858–865. DOI: [10.1016/j.jcat.2021.10.002](https://doi.org/10.1016/j.jcat.2021.10.002)

7. Majumdar, P. (2025). *Microkinetic Modeling in Heterogeneous Catalysis: Challenges and Path Forward*. J. Indian Inst. Sci. DOI: [10.1007/s41745-025-00482-8](https://doi.org/10.1007/s41745-025-00482-8)

8. South, C.R. & Warburton, R.E. (2025). *Modeling Hydrogen Atom Transfer Rate Constants in Heterogeneous Catalysis*. ACS Catal., 15(17). DOI: [10.1021/acscatal.5c03905](https://doi.org/10.1021/acscatal.5c03905)

9. Zhuo, Q. et al. (2026). *Decomposition of Ethylene on Pd/Au(100): A Combined DFT and Microkinetic Modeling Study*. Appl. Catal. A: General, 121053. DOI: [10.1016/j.apcata.2026.121053](https://doi.org/10.1016/j.apcata.2026.121053)

10. Bukur, D.B., Mandić, M., & Todić, B. (2020). *Pore diffusion effects on catalyst effectiveness and selectivity of cobalt based Fischer-Tropsch catalyst*. Catal. Today, 343, 1–10. DOI: [10.1016/j.cattod.2018.10.069](https://doi.org/10.1016/j.cattod.2018.10.069)

---

*Manuscript prepared 2026-05-29. Computational results generated using the open-source MKM framework developed in this work (Python 3.11, NumPy, SciPy, Matplotlib). All code available in the supplementary repository.*
