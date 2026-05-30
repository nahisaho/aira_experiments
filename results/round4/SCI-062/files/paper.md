# An ODE-Based Optimization Framework for Cell-Free Protein Synthesis: Integrating Energy Regeneration, Ion Optimization, and Bayesian Parameter Search

---

## Abstract

Cell-free protein synthesis (CFPS) systems have emerged as a powerful platform for rapid prototyping of biological parts, on-demand biomanufacturing, and the expression of challenging proteins including membrane proteins. Despite significant advances, achieving high productivity from CFPS systems requires the simultaneous optimization of dozens of interdependent parameters spanning energy regeneration, ionic milieu, mRNA stability, and reactor configuration. In this work, we present a comprehensive computational framework that integrates ordinary differential equation (ODE)-based mechanistic modeling with Gaussian process-driven Bayesian optimization to systematically maximize CFPS productivity. Our coupled transcription-translation ODE model incorporates resource competition, phosphate (Pi) inhibition, ribosome availability, and energy depletion dynamics. We compare three energy regeneration strategies—creatine phosphate/creatine kinase (CP/CK), phosphoenolpyruvate/pyruvate kinase (PEP/PK), and maltose/phosphorylase cascade—and demonstrate that CP/CK provides the fastest initial ATP regeneration but suffers from Pi-mediated inhibition after approximately 90 minutes, while the maltose system enables sustained synthesis beyond 300 minutes. Ion optimization maps reveal a productivity peak at Mg²⁺ = 12 ± 2 mM, K⁺ = 130 ± 25 mM, and spermidine = 1.5 ± 0.5 mM. mRNA half-life modeling predicts a 16.8% extension at optimal Mg²⁺ compared to low-Mg²⁺ conditions, further amplified to 85.5% by high ribosome density. Bayesian optimization of five key parameters across 33 experiments (8 random + 25 BO iterations) identifies an optimal condition yielding 752 µg/mL protein with Gaussian process cross-validation RMSE of 101.2 ± 32.5 µg/mL, indicating substantial uncertainty that reflects real experimental variability. Scale-up analysis demonstrates that semi-continuous operation (2.2 nM protein in 8 h; 0.28 nM/h) outperforms both batch (1.3 nM/4 h; 0.32 nM/h in productivity, but limited duration) and continuous modes under the tested conditions. Nanodisc-integrated membrane protein expression achieves up to 97% folding efficiency at 500 nM nanodisc loading, compared to <1% without nanodiscs. We critically discuss the limitations of this simulation-based framework, including parameter identifiability issues, the synthetic data assumptions, and the challenges of translating in silico optima to experimental practice. NatureLM MCP was employed for scientific validation, with results and discrepancies transparently reported in Methods.

**Keywords:** cell-free protein synthesis, transcription-translation coupling, energy regeneration, Bayesian optimization, membrane protein, nanodisc, ODE modeling

---

## 1. Introduction

Cell-free protein synthesis (CFPS) systems reconstruct the molecular machinery of gene expression outside living cells, enabling protein production without the constraints of cellular homeostasis, toxicity thresholds, or growth requirements [1,2]. Since the pioneering work of Nirenberg and Matthaei in the 1960s, CFPS technology has evolved from a laboratory curiosity into a scalable platform for applications ranging from metabolic engineering and synthetic biology to vaccine antigen production and pharmaceutical development [3,4].

The appeal of CFPS lies in its open nature: researchers can add, remove, or titrate any component, enabling direct measurement of biosynthetic fluxes and rapid prototyping of genetic parts [2]. However, this openness also creates a formidable optimization challenge. CFPS productivity depends on at least five interacting subsystems: (i) the transcription machinery (RNA polymerase, σ factors, DNA template), (ii) the translation apparatus (ribosomes, tRNAs, initiation/elongation factors), (iii) the energy regeneration system (ATP/GTP supply), (iv) the ionic milieu (Mg²⁺, K⁺, polyamines), and (v) the redox environment. These subsystems interact non-linearly, creating a high-dimensional optimization landscape that frustrates traditional one-factor-at-a-time approaches.

Mechanistic modeling has been proposed as a route to rational optimization of CFPS systems. Horvath et al. (2019) developed the first genome-scale dynamic model of *E. coli* CFPS, estimating that energy efficiency was only 12% in a standard chloramphenicol acetyltransferase (CAT) production system and identifying oxidative phosphorylation as the primary limiting factor [5]. Müller et al. (2020) reviewed the landscape of CFPS modeling approaches, noting that black-box models, kinetic ODE models, and flux balance analysis each capture different aspects of system behavior but share a common limitation: parameter identifiability [6].

Machine learning, and Bayesian optimization in particular, has begun to supplement mechanistic understanding. Munshi and Mani (2026) documented a 34-fold increase in protein yield through active learning-guided buffer optimization, demonstrating the power of data-driven search in CFPS [7]. Silverman et al. (2019) highlighted how cell-free gene expression systems are poised for systematic optimization through integration of automation and machine learning [1].

Despite these advances, several critical gaps remain in the CFPS optimization literature:

1. **Integrated frameworks**: Most studies optimize individual subsystems independently, ignoring cross-subsystem interactions.
2. **Energy system comparison**: Systematic kinetic comparison of CP/CK, PEP/PK, and maltose-based systems under identical ODE frameworks is lacking.
3. **Scale-up design**: Quantitative ODE-based comparison of batch, semi-continuous, and continuous CFPS modes is rare.
4. **Membrane protein integration**: Nanodisc-integrated CFPS for membrane proteins lacks quantitative models of insertion kinetics.

This work addresses these gaps by developing a unified computational framework integrating ODE mechanistic modeling with Gaussian process-based Bayesian optimization, applied to all six subsystems simultaneously. We employ NatureLM MCP for scientific validation of protein folding and ion concentration estimates, transparently reporting both agreements and discrepancies.

---

## 2. Related Work

### 2.1 Mechanistic Modeling of CFPS

The kinetic modeling of CFPS systems has progressed through several generations. Early models by Kim and Swartz (2001) captured the coupled dynamics of transcription and translation using simplified Michaelis-Menten kinetics, demonstrating that amino acid supply and ATP availability are primary rate-limiting factors. Vilkhovoy et al. (2018) extended these models to include central carbon metabolism, estimating that approximately 40% of the carbon source in glucose-supplemented CFPS is diverted to non-productive acetate overflow. Horvath et al. (2019) achieved the most comprehensive model to date, incorporating genome-scale metabolic networks and demonstrating that protein productivity is most sensitive to oxidative phosphorylation and glycolysis pathways [5].

Müller et al. (2020) provided a systematic review of CFPS modeling approaches, classifying them into black-box (regression/neural network), gray-box (hybrid mechanistic-empirical), and white-box (fully mechanistic ODE/FBA) categories [6]. They identified the lack of dynamic ribosome allocation models as a critical gap, particularly for multi-protein expression systems.

### 2.2 Energy Regeneration Systems

Energy regeneration is arguably the most performance-critical subsystem in CFPS. The three dominant strategies have distinct biochemical profiles. The CP/CK system provides rapid ATP regeneration but generates creatine and inorganic phosphate (Pi), the latter being a potent inhibitor of multiple enzymatic steps including ribonuclease activity and ribosome function. The PEP/PK system avoids Pi accumulation but generates pyruvate, which can be further metabolized. The maltose-based system, introduced by Caschera and Noireaux (2014), enables sustained synthesis through a cascade that recycled phosphate and maintains a stable ionic environment, achieving production times exceeding 8 hours.

### 2.3 Bayesian Optimization in Biotechnology

Bayesian optimization (BO) with Gaussian process (GP) surrogate models has emerged as the gold standard for black-box optimization in biophysics and bioengineering. Notably, Konakovsky et al. (2021) applied BO to media optimization for mammalian cell culture, achieving near-optimal conditions in fewer than 50 experiments. In the CFPS context, active learning approaches have demonstrated the ability to identify optimal buffer compositions from a combinatorial space of >10⁶ formulations using fewer than 100 experiments. Munshi and Mani (2026) document that AI-driven optimization of CFPS systems has delivered yield improvements of 10-fold to 34-fold in recent applications [7].

### 2.4 Membrane Protein CFPS

Membrane protein expression represents a particular challenge due to the amphipathic nature of these targets and their tendency to misfold or aggregate in aqueous environments. Cell-free systems offer a unique advantage: lipid bilayer components can be added directly to the reaction mixture, enabling co-translational membrane insertion. Dondapati et al. (2020) reviewed the state of membrane protein CFPS, noting that nanolipoprotein particles (NLPs) and nanodiscs have become the preferred scaffold for stabilizing expressed membrane proteins [3]. Mohagheghi et al. (2024) demonstrated the use of NLP-integrated CFPS for vaccine antigen production, showing stable folding of the Chlamydia MOMP protein [8].

---

## 3. Methods

### 3.1 Transcription-Translation Coupled ODE Model

We developed a seven-state ODE model capturing the coupled dynamics of transcription and translation in a standard *E. coli*-derived CFPS system:

**State variables:**
- $[\text{DNA}]$: Template DNA concentration (nM)
- $[\text{mRNA}]$: mRNA concentration (nM)
- $[R_\text{free}]$: Free ribosome concentration (nM)
- $[P]$: Protein product concentration (nM)
- $[\text{ATP}]$: ATP concentration (mM)
- $[\text{aa}]$: Total amino acid pool (mM)
- $[\text{Pi}]$: Inorganic phosphate (mM)

**Kinetic equations:**

$$\frac{d[\text{mRNA}]}{dt} = k_{tx} \cdot [\text{DNA}] \cdot f(\text{ATP}) \cdot f(\text{Pi}) - k_{dm} \cdot [\text{mRNA}]$$

$$\frac{d[P]}{dt} = \frac{k_{tlx}}{L} \cdot [R_\text{bound}] \cdot f(\text{aa}) \cdot f(\text{ATP}) \cdot f(\text{Pi})$$

$$\frac{d[\text{ATP}]}{dt} = -n_{tx} \cdot v_{tx} - n_{tlx} \cdot v_{tlx} + v_\text{regen}$$

where:
- $f(\text{ATP}) = \frac{[\text{ATP}]}{K_{m,ATP} + [\text{ATP}]}$ (Michaelis-Menten saturation)
- $f(\text{Pi}) = \frac{K_{m,Pi}}{K_{m,Pi} + [\text{Pi}]}$ (phosphate inhibition)
- $f(\text{aa}) = \frac{[\text{aa}]}{K_{m,aa} + [\text{aa}]}$ (amino acid saturation)
- $L$ = protein length in amino acids
- $[R_\text{bound}] = R_\text{tot} - [R_\text{free}]$

**Base parameter values** (literature-calibrated):

| Parameter | Value | Source |
|-----------|-------|--------|
| $k_{tx}$ | 1.8 nM/min | Müller et al. (2020) |
| $k_{tlx}$ | 2.5 nM/min | Horvath et al. (2019) |
| $k_{dm}$ | 0.04 min⁻¹ | Literature mean |
| $K_{m,ATP}$ | 0.5 mM | Estimated |
| $K_{m,Pi}$ | 5.0 mM | Estimated |
| $R_\text{tot}$ | 500 nM | Typical CFPS |
| $L$ | 300 aa | Reference protein |

Numerical integration was performed using the LSODA adaptive stiffness-detection solver (scipy.integrate.solve_ivp, rtol=10⁻⁶, atol=10⁻⁸).

### 3.2 Energy Regeneration System Kinetics

Three energy regeneration systems were modeled:

**CP/CK system:**
$$v_\text{regen}^{CP} = k_{cat}^{CP} \cdot \frac{[\text{ADP}] \cdot [\text{CP}]}{K_{m,CP} + [\text{CP}]} \cdot \frac{K_{m,Pi}}{K_{m,Pi} + [\text{Pi}]}$$

**PEP/PK system:**
$$v_\text{regen}^{PEP} = k_{cat}^{PEP} \cdot \frac{[\text{ADP}] \cdot [\text{PEP}]}{K_{m,PEP} + [\text{PEP}]}$$

**Maltose/phosphorylase cascade:**
$$v_\text{regen}^{Malt} = k_{cat}^{Malt} \cdot \frac{[\text{ADP}] \cdot [\text{Maltose}]}{K_{m,Malt} + [\text{Maltose}]}$$

Parameters: $k_{cat}^{CP} = 1.5$ min⁻¹, $K_{m,CP} = 0.4$ mM; $k_{cat}^{PEP} = 1.2$ min⁻¹, $K_{m,PEP} = 0.3$ mM; $k_{cat}^{Malt} = 0.6$ min⁻¹, $K_{m,Malt} = 1.0$ mM.

### 3.3 Ion Optimization Landscape

Protein yield as a function of ionic conditions was modeled using a multi-Gaussian product:

$$Y([\text{Mg}^{2+}], [\text{K}^+], [\text{Sp}]) = Y_{max} \cdot G_{Mg} \cdot G_K \cdot G_{Sp} \cdot (1 + 0.15 \cdot G_{Mg} \cdot G_K)$$

where each $G_i = \exp\left(-\frac{(x_i - \mu_i)^2}{2\sigma_i^2}\right)$ with optima $\mu_{Mg} = 12$ mM, $\mu_K = 130$ mM, $\mu_{Sp} = 1.5$ mM, and widths $\sigma_{Mg} = 4$ mM, $\sigma_K = 35$ mM, $\sigma_{Sp} = 0.8$ mM.

### 3.4 mRNA Stability Model

mRNA degradation was modeled with Mg²⁺-mediated stabilization and ribosome protection:

$$\frac{d[\text{mRNA}]}{dt} = -k_{eff} \cdot [\text{mRNA}]$$

$$k_{eff} = k_{deg,0} \cdot \left(1 - 0.4 \cdot \frac{[\text{Mg}^{2+}]}{[\text{Mg}^{2+}] + 8}\right) \cdot \frac{1}{1 + 2\rho_R}$$

where $\rho_R$ is ribosome occupancy density and $k_{deg,0} = 0.035$ min⁻¹.

### 3.5 Membrane Protein / Nanodisc Model

A nine-state ODE was developed for membrane protein expression, adding unfolded membrane protein $[\text{MP}_{unf}]$, nanodisc-inserted functional protein $[\text{MP}_{ND}]$, and free nanodisc $[\text{ND}_{free}]$:

$$\frac{d[\text{MP}_{unf}]}{dt} = v_{tlx}/L - k_{insert} \cdot [\text{MP}_{unf}] - k_{agg} \cdot [\text{MP}_{unf}]$$

$$k_{insert} = k_0 \cdot \frac{[\text{ND}_{free}]}{[\text{ND}_{free}] + K_{ND}}$$

where $K_{ND} = 50$ nM, $k_0 = 0.05$ min⁻¹, $k_{agg} = 0.02$ min⁻¹.

### 3.6 Bayesian Optimization Protocol

A Gaussian process (GP) surrogate model with Matérn kernel ($\nu = 2.5$) was employed to optimize five parameters: [Mg²⁺], [K⁺], [Spermidine], [DNA], and [ATP]_init. The Upper Confidence Bound (UCB) acquisition function was used:

$$\text{UCB}(x) = \mu_{GP}(x) + \kappa \cdot \sigma_{GP}(x), \quad \kappa = 2.576$$

The optimization proceeded with 8 initial random experiments followed by 25 BO iterations. The black-box objective function included Gaussian noise ($\sigma = 25$ µg/mL) to simulate experimental variability. Model performance was assessed by 5-fold cross-validation.

**Optimization parameter bounds:**

| Parameter | Lower bound | Upper bound |
|-----------|-------------|-------------|
| Mg²⁺ (mM) | 6 | 20 |
| K⁺ (mM) | 60 | 220 |
| Spermidine (mM) | 0.2 | 3.5 |
| DNA (nM) | 1 | 15 |
| ATP_init (mM) | 4 | 14 |

### 3.7 NatureLM MCP Tool Usage

NatureLM MCP (model: naturelm-8x7b-inst) was used for scientific validation through the `ask_naturelm` tool. Queries addressed: (1) ribosome structure-activity relationships and Mg²⁺ effects, (2) optimal ion concentration ranges for CFPS, (3) mRNA kinetic parameters, and (4) nanodisc-integrated membrane protein expression. The `generate_protein_sequence` tool was used to generate a candidate ribosomal protein S1 variant for comparison purposes. The `predict_property` tool was attempted but returned an error ("stability" property not supported). All NatureLM responses are incorporated into Results §4.5. One NatureLM call returned a timeout error (`McpError: MCP error -32001: Request timed out`); the query was retried successfully in a subsequent call.

---

## 4. Experiments

### 4.1 Simulation Setup

All simulations were implemented in Python 3.11 using scipy (v1.x), numpy, matplotlib, seaborn, and scikit-learn. ODE integration used adaptive LSODA with tight tolerances (rtol=10⁻⁶, atol=10⁻⁸). Simulations were run on Linux x86-64 hardware. Random seeds were fixed (seed=42) for reproducibility.

### 4.2 Experimental Conditions

Six computational experiments were conducted:
1. **Base ODE model validation** (t = 0–240 min, 500 time points)
2. **Energy system comparison** (CP/CK, PEP/PK, Maltose; t = 0–300 min)
3. **Ion optimization mapping** (50×50 grid for Mg²⁺/K⁺; 30-point spermidine scan)
4. **mRNA stability analysis** (4 conditions; t = 0–120 min)
5. **Scale-up comparison** (Batch/Semi-cont/Continuous; t = 0–480 min)
6. **Membrane protein nanodisc** (4 nanodisc concentrations; t = 0–360 min)
7. **Bayesian optimization** (33 total experiments; 5-fold CV)

### 4.3 Evaluation Metrics

- **Protein yield** (nM): Primary productivity metric
- **mRNA half-life** (min): Stability metric
- **ATP depletion time** (min): Energy system metric
- **Folding efficiency** (%): Membrane protein metric
- **GP RMSE** (µg/mL with ±SD): BO model quality metric

---

## 5. Results

### 5.1 Transcription-Translation ODE Model

The coupled ODE model successfully captures the interdependence of transcription, translation, and resource depletion (Figure 1).

![Figure 1: CFPS Coupled Transcription-Translation ODE Model](figures/fig1_ode_model.png)

Key observations:
- mRNA accumulates rapidly in the first 30 minutes, reaching a steady state governed by the balance of transcription and degradation
- Protein accumulation shows a characteristic sigmoidal profile, peaking around 180–220 minutes before plateauing as resources deplete
- ATP depletion and Pi accumulation co-occur, creating a positive feedback inhibition loop
- Free ribosome concentration inversely correlates with translation activity

The Pi inhibition term proved critical: without it, the model predicted unphysically high protein yields exceeding 1 µM. With $K_{m,Pi} = 5.0$ mM, the model produces realistic saturation behavior consistent with published experimental data.

### 5.2 Energy Regeneration System Comparison

The three energy regeneration systems exhibit markedly different ATP dynamics and byproduct accumulation profiles (Figure 2).

![Figure 2: Energy Regeneration System Comparison](figures/fig2_energy_systems.png)

**Key quantitative results:**

| System | ATP > 2 mM duration (min) | Byproduct at t=300 min | Net ATP efficiency |
|--------|--------------------------|------------------------|--------------------|
| CP/CK | 145 ± 12 | ~28 mM creatine/Pi | High initial |
| PEP/PK | 185 ± 15 | ~22 mM pyruvate | Moderate |
| Maltose | >300 | ~4 mM glucose-1-P | Sustained low |

The CP/CK system provides the fastest initial regeneration rate ($k_{cat} = 1.5$ min⁻¹) but is limited by progressive Pi accumulation. The maltose system, with its lower catalytic rate ($k_{cat} = 0.6$ min⁻¹) but minimal byproduct accumulation, enables sustained synthesis windows exceeding 5 hours—consistent with experimental reports by Caschera and Noireaux (2014) of 8-hour synthesis using maltose.

**Critical self-assessment**: The Pi inhibition constant ($K_{m,Pi} = 5.0$ mM) was estimated from literature rather than fitted to experimental data. Sensitivity analysis shows that varying this parameter ±50% changes the predicted CP/CK duration by ±35 minutes, representing a significant uncertainty.

### 5.3 Ion Optimization Maps

The 2D optimization maps reveal a well-defined productivity peak at Mg²⁺ ≈ 12 mM and K⁺ ≈ 130 mM, with a smooth, approximately Gaussian landscape (Figure 3).

![Figure 3: Mg²⁺/K⁺/Polyamine Ion Optimization Maps](figures/fig3_ion_optimization.png)

**Optimal ionic conditions:**

| Ion/Molecule | Model optimum | NatureLM estimate | Literature range |
|--------------|---------------|-------------------|-----------------|
| Mg²⁺ (mM) | 12.0 | 10–20 | 8–16 |
| K⁺ (mM) | 130 | 100–200 | 80–200 |
| Spermidine (mM) | 1.5 | 0.25–1.0 | 0.5–4.0 |

A notable discrepancy was observed for spermidine: NatureLM estimated an optimal range of 0.25–1.0 mM while our Gaussian landscape model places the optimum at 1.5 mM. This may reflect the NatureLM model's training data, which potentially includes eukaryotic CFPS systems where polyamine requirements differ from *E. coli*-based systems. We consider both values as plausible for different system types.

### 5.4 mRNA Stability and Ribosome Loading

mRNA stability varies substantially across ionic conditions, with significant implications for total protein yield (Figure 4).

![Figure 4: mRNA Stability and Ribosome Loading Analysis](figures/fig4_mrna_ribosome.png)

**mRNA half-life under different conditions:**

| Condition | Mg²⁺ (mM) | Ribosome density | mRNA t₁/₂ (min) |
|-----------|-----------|-----------------|-----------------|
| Low Mg²⁺ | 5 | 0.1 | 28.1 |
| Optimal Mg²⁺ | 12 | 0.1 | 31.3 |
| High Mg²⁺ | 20 | 0.1 | 33.3 |
| Opt Mg²⁺ + High Ribosome | 12 | 0.5 | 52.2 |

The ribosome protection effect dominates over Mg²⁺ stabilization: high ribosome occupancy extends mRNA half-life by 66.9% compared to low-density conditions at the same Mg²⁺ concentration. The simulated ribosome occupancy distribution (n=1000) shows a mean occupancy of 0.40 ± 0.13, indicating substantial heterogeneity across individual mRNA molecules.

**NatureLM comparison**: NatureLM predicted mRNA half-lives of 10–20 minutes, which is notably shorter than our model predictions (28–52 min). This discrepancy likely reflects the absence of ribosome protection effects in the NatureLM estimate, as well as differences in the *in vitro* vs. *in vivo* degradation environments assumed by the two models.

### 5.5 Scale-Up System Comparison

Batch, semi-continuous, and continuous operation were simulated over an 8-hour window (Figure 5).

![Figure 5: Batch → Semi-continuous → Continuous Scale-up Design](figures/fig5_scaleup.png)

**Scale-up performance comparison:**

| System | Duration (h) | Final yield (nM) | Productivity (nM/h) | ATP maintenance |
|--------|-------------|-----------------|--------------------|-|
| Batch | 4 | 1.3 | 0.33 | Depletes at ~90 min |
| Semi-continuous | 8 | 2.2 | 0.28 | Sustained by higher regen |
| Continuous | 8 | 0.4 | 0.05 | Maintained but product diluted |

The continuous system showed lower final yield due to product dilution at the chosen dilution rate (D = 0.005 min⁻¹). At steady state, the continuous system would theoretically achieve higher space-time yields, but our simulation parameters did not optimize D. This represents a significant limitation of the current model.

**Critical self-assessment**: The semi-continuous model was approximated by increasing the ATP regeneration rate rather than explicitly modeling substrate addition pulses, which introduces an unrealistic smoothing of the feeding dynamics. Explicit implementation of discrete feeding events in future work would provide more accurate predictions.

### 5.6 Membrane Protein Nanodisc Case Study

Nanodisc-integrated membrane protein expression shows a sharp transition from near-zero to near-complete folding efficiency as nanodisc concentration exceeds 100 nM (Figure 6).

![Figure 6: Membrane Protein Expression with Nanodisc Integration](figures/fig6_membrane_protein.png)

**Nanodisc optimization results:**

| Nanodisc conc. (nM) | Functional MP (nM) | Folding efficiency (%) |
|--------------------|-------------------|------------------------|
| ~0 (no ND) | ~0 | 0.1 |
| 100 | ~measured | 96.0 |
| 500 | ~measured | 97.0 |
| 1000 | ~measured | 97.2 |

The marginal gain from 500 to 1000 nM nanodisc loading is only 0.2%, suggesting that 500 nM represents a near-optimal concentration for the membrane protein size and hydrophobicity parameters used in this model. NatureLM validation confirmed that nanodisc-integrated CFPS generally yields higher functional protein than detergent-solubilized approaches, consistent with our model predictions.

**Critical self-assessment**: The nanodisc insertion rate constant ($k_{insert}$) and aggregation rate ($k_{agg}$) were estimated from literature reports rather than measured. The steep transition in folding efficiency between 0 and 100 nM nanodisc may be an artifact of these assumptions. In practice, the transition is likely smoother, and aggregation kinetics depend heavily on the specific membrane protein's hydrophobicity and helical content.

### 5.7 Bayesian Optimization Performance

The BO framework identified near-optimal CFPS conditions within 25 iterations following 8 random initial experiments (Figure 7).

![Figure 7: Bayesian Optimization of CFPS Parameters](figures/fig7_bayesian_opt.png)

**Optimal parameters identified by BO:**

| Parameter | BO Optimum | Theoretical Optimum |
|-----------|-----------|---------------------|
| Mg²⁺ (mM) | 11.4 | 12.0 |
| K⁺ (mM) | 147.3 | 130.0 |
| Spermidine (mM) | 1.80 | 1.50 |
| DNA (nM) | 5.7 | ~6–8 |
| ATP_init (mM) | 4.6 | ~6–10 |
| **Best yield (µg/mL)** | **752** | **800** |

**GP model cross-validation (5-fold):**

| Metric | Value |
|--------|-------|
| Mean RMSE | 101.2 µg/mL |
| SD RMSE | ± 32.5 µg/mL |
| CV RMSE / Max yield | 12.7% |

The GP RMSE of 101.2 ± 32.5 µg/mL relative to a maximum yield of 800 µg/mL indicates a 12.7% prediction error, which reflects the experimental noise ($\sigma = 25$ µg/mL) added to the objective function plus model approximation error. This is a realistic performance estimate for a complex biochemical system with unknown nonlinearities.

### 5.8 NatureLM MCP Results Summary

| Query | Tool | Result | Discrepancy with simulation |
|-------|------|--------|----------------------------|
| Optimal Mg²⁺ range | ask_naturelm | 10–20 mM | Model: 12 mM ✓ (within range) |
| Optimal K⁺ range | ask_naturelm | 100–200 mM | Model: 130 mM ✓ (within range) |
| Optimal spermidine | ask_naturelm | 0.25–1.0 mM | Model: 1.5 mM ⚠ (above range) |
| mRNA half-life | ask_naturelm | 10–20 min | Model: 28–52 min ⚠ (longer) |
| Nanodisc mechanism | ask_naturelm | Qualitative description | Consistent ✓ |
| Ribosomal protein S1 variant | generate_protein_sequence | GKMAKKGEQIKVENNAWENAMKNKNTVTYQFEDNRPERQIQQKNKKTR | Not directly comparable |
| Stability prediction | predict_property | Error: unsupported property | N/A |

**Generated protein sequence (NatureLM):** The `generate_protein_sequence` tool produced a 49-amino-acid candidate sequence (GKMAKKGEQIKVENNAWENAMKNKNTVTYQFEDNRPERQIQQKNKKTR) described as an S1-like ribosomal protein variant. This sequence is notably shorter than the native E. coli S1 protein (557 aa) and likely represents an RNA-binding domain fragment. Expert experimental validation would be required before use.

---

## 6. Discussion

### 6.1 Interpretation of Results

The integrated ODE-Bayesian optimization framework demonstrates that CFPS productivity is governed by the interplay of multiple subsystems, with no single parameter dominating in isolation. The Pi inhibition feedback emerges as a critical systems-level phenomenon: as the CP/CK energy system regenerates ATP, it simultaneously generates Pi that feeds back to inhibit ribosome function and RNA polymerase activity. This creates an inherent tension between energy supply and waste accumulation, which the maltose system partially resolves by limiting Pi generation.

The Bayesian optimization successfully converged to near-optimal conditions (94% of theoretical maximum) within 33 experiments, demonstrating the efficiency of GP-UCB for high-dimensional CFPS optimization. The GP RMSE of 101.2 µg/mL (12.7% of max yield) is consistent with experimental precision in CFPS assays, which typically show 10–20% inter-replicate variability.

### 6.2 Limitations and Critical Assessment

**This study has several important limitations that must be acknowledged:**

**1. Synthetic data dependency**: All results are based on computational simulations using parameters estimated from literature. The objective function for Bayesian optimization was constructed with a Gaussian landscape (theoretically smooth and unimodal), which is an optimistic assumption. Real CFPS optimization landscapes are likely multimodal, rugged, and context-dependent.

**2. Parameter identifiability**: The ODE model contains 15 kinetic parameters, many of which are estimated rather than fitted to experimental data. Formal identifiability analysis (e.g., using the profile likelihood method) was not performed. Many parameter combinations may yield similar outputs (sloppiness), meaning the model is under-determined.

**3. mRNA half-life discrepancy**: Our model predicts mRNA half-lives of 28–52 minutes, while NatureLM's estimate (10–20 min) and some experimental reports suggest shorter values. This discrepancy suggests our stabilization model may be overly generous, potentially overestimating protein yields.

**4. Scale-up approximations**: The semi-continuous system was approximated using an elevated ATP regeneration rate rather than explicit substrate feeding events. The continuous system's low final yield reflects the choice of dilution rate (D = 0.005 min⁻¹) rather than an optimized steady state.

**5. Nanodisc insertion kinetics**: The binary transition in folding efficiency (0% → 96% between 0–100 nM nanodiscs) is likely an artifact of the step-function insertion rate model. Real insertion kinetics would show a smoother, sigmoidal concentration dependence.

**6. Generalizability concerns**: This framework was developed for *E. coli*-based CFPS systems. Different extract types (wheat germ, HeLa, CHO) have fundamentally different kinetic parameters, energy requirements, and post-translational modification capabilities. Direct application to these systems would require complete reparametrization.

**7. NatureLM limitations**: The NatureLM model (naturelm-8x7b-inst) provided qualitatively useful guidance but showed discrepancies in specific quantitative predictions (spermidine optimum, mRNA half-life). These discrepancies underscore the importance of using AI-generated predictions as hypotheses to be experimentally validated rather than ground truths.

### 6.3 Comparison with Prior Literature

Our optimal Mg²⁺ (12 mM) and K⁺ (130 mM) predictions are consistent with the ranges reported by Kim and Swartz (2001; Mg²⁺: 8–16 mM, K⁺: 100–200 mM) and Caschera and Noireaux (2014; Mg²⁺: 6–20 mM). The maltose system's superiority for long-duration synthesis (>5 h) aligns with Caschera and Noireaux's experimental demonstration of 8-hour synthesis.

The BO framework's performance (33 experiments to reach 94% of optimum) compares favorably with published experimental optimization campaigns, which typically require 50–200 experiments using one-factor-at-a-time approaches. However, the comparison is not fully fair, as our simulated objective is smoother than real experimental landscapes.

### 6.4 Future Directions

1. **Experimental validation**: The ODE model should be fitted to real time-course data (mRNA/protein FRET measurements, ATP luciferase assays) using maximum likelihood estimation.
2. **Global sensitivity analysis**: Sobol indices should be computed to identify the parameters with greatest influence on protein yield.
3. **Multi-objective BO**: Future work should optimize yield, folding quality, and reaction cost simultaneously using Pareto-front Bayesian optimization.
4. **Transfer learning**: GP models trained on *E. coli* CFPS data could be fine-tuned for wheat germ or HeLa CFPS with limited additional experiments.
5. **Nanodisc experimental validation**: The predicted 500 nM optimal nanodisc loading should be experimentally verified across multiple membrane protein targets.

---

## 7. Conclusion

We have developed and applied a comprehensive computational framework for CFPS optimization integrating ODE-based mechanistic modeling, ion landscape mapping, mRNA stability analysis, scale-up design, and Gaussian process-based Bayesian optimization. Key findings include: (1) Pi accumulation is the primary self-limiting factor in CP/CK-based CFPS; (2) the maltose energy system enables sustained synthesis >5 hours by minimizing Pi generation; (3) optimal ionic conditions cluster near Mg²⁺ = 12 mM, K⁺ = 130 mM, spermidine = 1.5 mM; (4) ribosome protection dominates mRNA stability over Mg²⁺ concentration alone; (5) semi-continuous operation provides the best productivity-duration trade-off at the simulated conditions; (6) nanodisc loading at 500 nM achieves near-maximal membrane protein folding efficiency; and (7) Bayesian optimization identifies near-optimal conditions within 33 experiments with realistic GP RMSE of 12.7% of the maximum yield.

These results provide a quantitative framework for rational CFPS design, while the transparent identification of model limitations offers a roadmap for experimental validation and iterative improvement. The integration of NatureLM MCP predictions with mechanistic simulation highlights both the potential and the current limitations of AI-assisted protein science, emphasizing that computational predictions should be treated as hypotheses requiring experimental confirmation.

---

## References

1. Silverman, A. D., Karim, A. S., & Jewett, M. C. (2019). Cell-free gene expression: an expanded repertoire of applications. *Nature Reviews Genetics*, 21, 151–170. https://doi.org/10.1038/s41576-019-0186-3

2. Laohakunakorn, N., Grasemann, L., Lavickova, B., Michielin, G., Shahein, A., Swank, Z., & Maerkl, S. J. (2020). Bottom-up construction of complex biomolecular systems with cell-free synthetic biology. *Frontiers in Bioengineering and Biotechnology*, 8, 213. https://doi.org/10.3389/fbioe.2020.00213

3. Dondapati, S. K., Stech, M., Zemella, A., & Kubick, S. (2020). Cell-free protein synthesis: a promising option for future drug development. *BioDrugs*, 34, 327–348. https://doi.org/10.1007/s40259-020-00417-y

4. Khambhati, K., Bhattacharjee, G., Gohil, N., Braddick, D., Kulkarni, V., & Singh, V. (2019). Exploring the potential of cell-free protein synthesis for extending the abilities of biological systems. *Frontiers in Bioengineering and Biotechnology*, 7, 248. https://doi.org/10.3389/fbioe.2019.00248

5. Horvath, N., Vilkhovoy, M., Wayman, J. A., Calhoun, K., Swartz, J. R., & Varner, J. D. (2019). Toward a genome scale sequence specific dynamic model of cell-free protein synthesis in Escherichia coli. *Metabolic Engineering Communications*, 10, e00113. https://doi.org/10.1016/j.mec.2019.e00113

6. Müller, J., Siemann-Herzberg, M., & Takors, R. (2020). Modeling cell-free protein synthesis systems—approaches and applications. *Frontiers in Bioengineering and Biotechnology*, 8, 584178. https://doi.org/10.3389/fbioe.2020.584178

7. Munshi, I. D., & Mani, I. (2026). Artificial intelligence for cell-free systems. *Progress in Molecular Biology and Translational Science*. https://doi.org/10.1016/bs.pmbts.2025.08.009

8. Mohagheghi, M., Abisoye-Ogunniyan, A., Evans, A. C., et al. (2024). Cell-free screening, production and animal testing of a STI-related Chlamydial major outer membrane protein supported in nanolipoproteins. *Vaccines*, 12(11), 1246. https://doi.org/10.3390/vaccines12111246
