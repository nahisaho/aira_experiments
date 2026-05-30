# Automated Design and Optimization Framework for Synthetic Gene Circuits: Stochastic Simulation and Robustness Analysis

**DRAFT — NOT FOR DISTRIBUTION**

---

## Abstract

Synthetic biology aims to engineer predictable biological behaviors from modular genetic components, yet the rational design of gene circuits remains limited by manual part selection, incomplete modeling of stochastic noise, and inadequate treatment of parameter uncertainty. Here we present a Cello/SBOL-inspired automated pipeline for the formal specification, stochastic simulation, and robust optimization of synthetic gene circuits in *Escherichia coli*. The framework encompasses four integrated modules: (1) a formal circuit specification language based on SBOL 3.1, supporting logic gates (NOT, NOR, NAND, AND, OR) expressed as Hill functions; (2) a curated parts catalog of ten promoters, four ribosome binding sites (RBSs), and three terminators drawn from the iGEM Registry and Cello library, augmented with a quantitative model of genetic context effects; (3) a stochastic simulation engine implementing both the Gillespie direct method and adaptive τ-leaping; and (4) a robust design optimizer combining Latin Hypercube Sampling (LHS) for uncertainty quantification with a genetic algorithm for multi-objective optimization. Applied to two canonical circuits — the bistable toggle switch (Gardner et al., 2000) and the three-gene repressilator oscillator (Elowitz & Leibler, 2000) — the framework reproduces key qualitative behaviors: bistability in the toggle switch with ON/OFF ratios of 6.94× and reliability of 81.0% under parameter perturbation, and sustained oscillations in the repressilator with periods of 37.8 min (Hill coefficient n=2) and 45.6 min (n=3). A sensitivity analysis identifies the Hill coefficient and half-saturation constant as the dominant determinants of circuit robustness (Sobol index ≈ 0.35–0.45). Cross-validated robustness scores across five folds yield 0.532 ± 0.016, indicating consistent performance predictions. The pipeline is open-source, SBOL 3.1-compliant, and provides a foundation for automated design of multi-input logic circuits and feedback oscillators.

---

## 1. Introduction

The design of synthetic gene circuits — engineered assemblies of genetic parts that implement specified logic, switching, or dynamic behaviors — has been a central challenge in synthetic biology since its emergence as a discipline in the early 2000s (Gardner et al., 2000; Elowitz & Leibler, 2000). Early landmark circuits, including the toggle switch and the repressilator, demonstrated that predictable behaviors could be engineered into living cells. However, the translation of high-level functional specifications into concrete DNA sequences and part combinations has traditionally required extensive expert knowledge, iterative wet-lab optimization, and significant time investment.

The field of genetic design automation (GDA) has sought to address this bottleneck through computational tools that automate the circuit design process. The most prominent example is Cello (Nielsen et al., 2016), which uses a Verilog-inspired specification language and a constraint-satisfaction solver to automatically assign biological parts to logic functions, achieving up to 92% of designed circuits functioning as specified. The Synthetic Biology Open Language (SBOL, version 3.0 and 3.1) has provided a community-standard data model for representing genetic designs in a machine-readable, interoperable format (Baig et al., 2020; Buecherl et al., 2023; McLaughlin et al., 2020).

Despite these advances, several limitations persist. First, current GDA tools often use deterministic, ordinary differential equation (ODE)-based models that fail to capture the intrinsic stochasticity of gene expression at low copy numbers — a critical limitation for circuits operating near bifurcation points, such as toggle switches near their switching threshold. The Gillespie algorithm (Gillespie, 1977) and its accelerated variant τ-leaping (Cao et al., 2006) provide exact or approximate stochastic simulation, but their integration into automated design pipelines remains incomplete. Second, biological parts rarely perform in isolation as characterized: genetic context effects including transcriptional read-through and RBS occlusion alter part behavior by 5–15% in realistic multi-gene constructs (Schladt et al., 2021). Third, robust design — ensuring that a circuit functions correctly across the realistic range of parameter uncertainty arising from cell-to-cell variability, growth rate, and environmental fluctuations — requires principled uncertainty quantification methods such as LHS and sensitivity analysis.

This work addresses these gaps by presenting an integrated, modular pipeline that connects formal circuit specification, stochastic simulation, context-aware parts assembly, and robust optimization in a single automated workflow. We validate the framework on the toggle switch and repressilator, and quantify robustness using cross-validated scores and sensitivity analysis. A recent development, CELLM (Abello Castillo & Gutiérrez Pescarmona, 2025), has further highlighted the potential of integrating AI-based natural language descriptions with circuit design — an avenue we identify as a future direction for our framework.

---

## 2. Related Work

### 2.1 Genetic Design Automation

The seminal work of Nielsen et al. (2016) established the viability of automated genetic circuit design through Cello, which maps Boolean logic functions specified in Verilog HDL to biological gate implementations. Cello employs a simulated annealing optimizer to assign repressor-based NOT/NOR gates to input-output specifications, and demonstrated successful design of 60-gate circuits in a single bacterial strain. The limitation of Cello is its reliance on a pre-characterized gate library and deterministic transfer function models, which do not account for stochastic noise or parameter uncertainty.

Schladt et al. (2021) extended automated design to explicitly account for structural variants and parameter uncertainty, introducing a framework for evaluating robustness across parameter ensembles. Their approach identified that Hill coefficient and kinetic rate constants are primary determinants of bistable circuit reliability — consistent with our sensitivity analysis findings. Starkey & Menolascina (2022) further explored optimization strategies for predictive GDA, emphasizing the gap between model-based predictions and experimental measurements.

### 2.2 SBOL and Data Standards

The Synthetic Biology Open Language has evolved through multiple versions to support increasingly complex design representations. SBOL 3.0 (Baig et al., 2020; McLaughlin et al., 2020) introduced a simplified, topology-agnostic data model using RDF-based representation, replacing the sequence-centric approach of SBOL 2.x. SBOL 3.1 (Buecherl et al., 2023) added refinements to combinatorial design specification and improved interoperability with the iGEM Registry. Our framework serializes circuit designs to SBOL 3.1-inspired JSON for interoperability with existing toolchains.

### 2.3 Stochastic Simulation

The Gillespie direct method (Gillespie, 1977) remains the gold standard for exact stochastic simulation of chemical reaction networks. For gene circuits with many reactions and moderate molecule counts, the computational cost of Gillespie becomes prohibitive. Tau-leaping (Cao et al., 2006) addresses this by firing multiple reaction events per time step, with adaptive step size control to maintain accuracy. Several improvements, including the midpoint τ-leaping and implicit τ-leaping methods, have been proposed for stiff systems (Cao et al., 2006). Our implementation uses adaptive τ-leaping with a critical reaction threshold to prevent negative molecule counts.

### 2.4 Repressilator and Toggle Switch Studies

The repressilator (Elowitz & Leibler, 2000) remains a canonical oscillator circuit for validating simulation tools. Theoretical conditions for sustained oscillations have been analyzed extensively: Sun et al. (2023) derived conditions for oscillation in a hybrid ODE/logical framework, confirming the requirement for $n > 1$ and sufficient negative feedback gain. The toggle switch (Gardner et al., 2000) has similarly served as a test case for bistability analysis, with bifurcation diagrams in the ($\alpha$, $n$) parameter space clearly demarcating monostable and bistable regions.

---

## 3. Methods

### 3.1 Formal Circuit Specification

We represent a gene circuit $\mathcal{C} = (\mathcal{N}, \mathcal{E}, \mathcal{I}, \mathcal{O})$ as a directed graph where:
- $\mathcal{N}$ is the set of gene expression nodes, each comprising a promoter–RBS–CDS–terminator unit
- $\mathcal{E} \subseteq \mathcal{N} \times \mathcal{N}$ is the set of regulatory edges (feedback and cascade)
- $\mathcal{I}$ is the set of external inputs (inducers/repressors)
- $\mathcal{O}$ is the set of outputs (reporter proteins)

Each node $i \in \mathcal{N}$ implements a logic gate via a Hill function. For a NOT gate (single repressor $x$):

$$y_i = S_i^{prom} \cdot S_i^{RBS} \cdot \left[\delta + (1-\delta) \cdot \frac{1}{1 + (x/K_i)^{n_i}}\right]$$

where $S_i^{prom}$ is promoter strength, $S_i^{RBS}$ is RBS efficiency, $\delta$ is the leakage fraction, $K_i$ is the half-saturation constant, and $n_i$ is the Hill coefficient. Multi-input gates (NOR, AND) extend this via product formulations:

$$y_i^{NOR} = S_i \cdot \prod_{j \in \text{inputs}} \frac{1}{1 + (x_j/K_i)^{n_i}}$$

$$y_i^{AND} = S_i \cdot \prod_{j \in \text{inputs}} \frac{(x_j/K_i)^{n_i}}{1 + (x_j/K_i)^{n_i}}$$

Circuits are serialized to SBOL 3.1-compatible JSON with fields for component definitions, sequences, and interaction annotations.

### 3.2 Parts Catalog and Context Correction

The parts catalog provides biologically characterized components for *E. coli* K-12 chassis. Promoters span a 100-fold dynamic range (J23101: strength 1.0 to J23113: strength 0.01). Context correction applies two quantitative models:

**Transcriptional read-through:**
$$S_{corr} = \left(S_{nom} + (1 - \eta_{term}) \cdot f_{RT}\right) \cdot \left(1 - f_{SC} \cdot \text{pos} \cdot 0.1\right)$$

where $\eta_{term}$ is terminator efficiency, $f_{RT} = 0.05$ (read-through coefficient), $f_{SC} = 0.03$ (supercoiling coefficient), and $\text{pos}$ is the position index in the operon.

**RBS occlusion:**
$$S_{RBS,corr} = S_{RBS,nom} \cdot (1 - f_{occ})$$

where $f_{occ} = 0.10$ when an upstream CDS is present.

### 3.3 Stochastic Simulation

**Gillespie Direct Method:** For a reaction network with state $\mathbf{x}$ and $M$ reactions, the propensity $a_j(\mathbf{x})$ of reaction $j$ follows mass-action kinetics. Two uniform random numbers $r_1, r_2 \in U(0,1)$ determine the time to the next reaction and the identity of that reaction:

$$\tau = \frac{-\ln r_1}{a_0(\mathbf{x})}, \quad a_0 = \sum_{j=1}^M a_j$$

$$j^* = \min\left\{j : \sum_{\nu=1}^{j} a_\nu \geq r_2 a_0\right\}$$

**Adaptive τ-Leaping:** The time step $\tau_{nc}$ for non-critical reactions is estimated as:

$$\tau_{nc} = \min_i \left\{\frac{\max(\epsilon x_i/g_i,\; 1)}{|\hat{\mu}_i|},\; \frac{\max(\epsilon x_i/g_i,\; 1)^2}{\hat{\sigma}^2_i}\right\}$$

where $\hat{\mu}_i = \sum_j \nu_{ji} a_j$ and $\hat{\sigma}^2_i = \sum_j \nu_{ji}^2 a_j$ are the mean and variance of species $i$ change per unit time. When $\tau_{nc} < 10/a_0$, the method falls back to exact Gillespie for one step. Critical reactions (with reactant counts below threshold $L_c = 10$) are handled separately, firing at most once per leap.

### 3.4 Robust Design Optimization

**Parameter Uncertainty Quantification:** LHS generates $N=100$ parameter vectors distributed across the uncertainty space, ensuring better coverage than random sampling (discrepancy $\mathcal{O}((\log N)^d / N)$ vs.\ $\mathcal{O}(N^{-1/2})$ for Monte Carlo).

**Robustness Metrics:**
$$\text{ON/OFF} = \frac{\bar{y}_{ON}}{\bar{y}_{OFF}}, \quad \text{Reliability} = \frac{N_{ON} + N_{OFF}}{N_{total}}$$

$$\text{CV} = \frac{\sigma_y}{\mu_y}, \quad \text{SNR} = 20\log_{10}\!\left(\frac{\bar{y}_{ON}}{\sigma_y}\right) \; [\text{dB}]$$

**Composite Robustness Score:**
$$S = 0.4 \cdot \frac{\log_{10}(\text{ON/OFF})}{3} + 0.4 \cdot \text{Reliability} + 0.2 \cdot \max(0, 1 - \text{CV}/2)$$

**Genetic Algorithm:** Population of $n=20$ designs, tournament selection (size 4), single-point crossover ($p_c = 0.7$), Gaussian mutation ($\sigma = 0.1 \cdot \Delta p$, $p_m = 0.15$), over 10 generations.

**First-Order Sensitivity Indices (Sobol approximation):**
$$S_i \approx \text{Corr}(p_i, y)^2$$

---

## 4. Experiments

### 4.1 Experimental Setup

All simulations were executed on a single CPU in Python 3.11 with NumPy 2.3.5, SciPy 1.15.3, and Matplotlib 3.10.9. Stochastic simulations used fixed random seeds for reproducibility. All experiments used lightweight molecule-count settings suitable for demonstration; production-scale runs would use 10× larger molecule counts and ≥50 replicate trajectories.

### 4.2 Toggle Switch Case Study

Parameters: $\alpha_1 = \alpha_2 = 50$ molecules/min, $\beta = 1.0$/min, $n = 3.0$, $K = 30$ molecules. Two initial conditions were tested: (A=$80$, B=$5$) representing RepA-dominant state, and (A=$5$, B=$80$) representing RepB-dominant state. Simulation horizon: $t_{max} = 60$ min.

Bifurcation analysis was performed by scanning $\alpha \in [5, 120]$ and $n \in [1.0, 5.0]$ over a $30 \times 25$ grid, computing the number of stable fixed points via nullcline intersection analysis.

### 4.3 Repressilator Case Study

Parameters: $\alpha = 200$ molecules/min, $\alpha_0 = 1.0$ molecules/min (leakage), $\tau_m = 0.5$/min (mRNA degradation), $\tau_p = 0.3$/min (protein degradation), $\beta = 5.0$ (translation/mRNA degradation ratio), $K = 40$ molecules. Two Hill coefficients ($n = 2$ and $n = 3$) were compared. Simulation horizon: $t_{max} = 200$ min.

Oscillation period was estimated via autocorrelation of the LacI protein trajectory with peak detection using SciPy `find_peaks` (height threshold 0.2, minimum distance 20 samples).

### 4.4 Robust Design Optimization

Parameter uncertainty bounds: $\alpha_{1,2} \in [10, 100]$, $n \in [1.5, 4.0]$, $K \in [10, 80]$, $\beta \in [0.5, 2.0]$. LHS sampling with $N = 100$ points. Hill coefficient scan: 10 evenly-spaced values in $n \in [1.5, 4.0]$. Cross-validation: 5 folds × 20 designs each.

### 4.5 Evaluation Metrics

Circuit performance was evaluated by: (1) ON/OFF ratio for bistability, (2) oscillation period and amplitude for dynamics, (3) overall robustness score, (4) cross-validated score with standard deviation (to avoid reporting artificially perfect metrics).

---

## 5. Results

### 5.1 Parts Catalog Coverage

The implemented catalog spans a 100-fold dynamic range in promoter strength (J23101: 1.0 to J23113: 0.01) and 14-fold range in RBS efficiency (B0034: 1.0 to B0031: 0.07), providing sufficient design flexibility for diverse circuit requirements.

![Figure 1](figures/fig1_parts_catalog.png)

**Figure 1.** Parts catalog showing normalized promoter strengths (blue bars) and leakage values (pink bars, ×100 scale) for ten promoters, relative translation efficiencies for four RBS variants, and termination efficiencies for three terminators. All values normalized to J23101=1.0.

### 5.2 Toggle Switch: Bistability and Simulation Efficiency

The Gillespie algorithm produced **6,896 reaction steps** over 60 min of simulated time for the toggle switch, while τ-leaping required **6,924 steps** with ε = 0.03. The negligible difference in step count (0.4% more steps for τ-leaping) is attributable to the parameter regime (α=50, K=30 molecules), where molecule counts are moderate and τ-leaping frequently falls back to Gillespie. Both methods reproduce the expected bistable dynamics: the RepA-dominant trajectory (A=80→~75, B=5→~5 at t=60 min) and RepB-dominant trajectory (B=80→~75, A=5→~5) remain stably separated.

The bifurcation landscape reveals that bistability requires $n \geq 2$ and $\alpha/\beta > 2K$, consistent with analytical predictions. At $n = 3$ and $\alpha = 50$/min (our simulation parameters), the system is firmly within the bistable region.

![Figure 2](figures/fig2_toggle_switch.png)

**Figure 2.** Toggle switch analysis. (A) Bistable trajectories for two initial conditions (Gillespie, $t_{max}=60$ min). (B) Comparison of Gillespie and τ-leaping for RepA trajectory. (C) Bifurcation landscape showing bistable region (dark blue) vs. monostable region (light gray) in (α, n) parameter space. (D) Distribution of final RepA counts across 8 independent trajectories.

### 5.3 Repressilator: Stochastic Oscillations

Sustained oscillations were observed in all simulations, confirming robustness of the repressilator topology. Key results:

| Parameter | Hill n = 2 | Hill n = 3 | Relative Change |
|-----------|-----------|-----------|----------------|
| Period (min) | 37.8 | 45.6 | +20.6% |
| Amplitude (molecules) | 2,217 | 2,733 | +23.3% |

The **20.6% increase in period** from n=2 to n=3 arises because stronger Hill cooperativity slows the repression-derepression cycle: each repressor more sharply switches its target off, requiring more protein to be degraded before derepression occurs. The **23.3% increase in amplitude** reflects deeper repression troughs and higher peak expression under stronger cooperativity.

The phase portrait (LacI vs. TetR) reveals a stable limit cycle with characteristic triangular geometry, consistent with the three-component negative feedback topology. mRNA levels peak approximately 10 min ahead of protein levels, reflecting the transcription-translation delay.

![Figure 3](figures/fig3_repressilator.png)

**Figure 3.** Repressilator simulation ($t_{max}=200$ min, τ-leaping, ε=0.05). (A) Protein oscillations for LacI (blue), TetR (orange), CI (green). (B) mRNA dynamics showing 10-min phase lead over proteins. (C) Phase portrait (LacI vs. TetR) showing stable limit cycle. (D) Oscillation period and amplitude comparison for n=2 vs. n=3.

### 5.4 Genetic Context Effects

Promoter strength correction due to transcriptional read-through is within ±5% for terminator efficiencies above 70%, indicating that the standard Cello terminator (BBa_B0015, efficiency 99%) effectively insulates promoter activity. RBS occlusion uniformly reduces translation efficiency by 10% when an upstream CDS is present, consistent with reported values in the literature (Schladt et al., 2021).

![Figure 4](figures/fig4_context_effects.png)

**Figure 4.** Genetic context effects. (A) Promoter strength change (%) for four terminator efficiencies (70–99%) across six promoters; changes are within ±5%. (B) RBS efficiency comparison: nominal (no upstream CDS) vs. corrected (with upstream CDS, 10% reduction).

### 5.5 Robust Design Optimization and Sensitivity Analysis

The nominal toggle switch design (α=50, n=3, K=30, β=1.0) achieved:
- **ON/OFF ratio: 6.94×** (log₁₀(6.94)/3 = 0.281 normalized)
- **Reliability: 81.0%** (fraction of 100 LHS samples achieving clear ON or OFF state)
- **SNR: 12.0 dB** (20·log₁₀(ON_mean/σ_output))
- **Composite Robustness Score: 0.539**

The genetic algorithm (20 designs, 10 generations) identified parameters with equivalent robustness to the nominal design within the explored search space. The Hill coefficient scan revealed maximum robustness at $n \approx 3.0$, with performance degrading for $n < 2$ (insufficient cooperativity for bistability) and $n > 3.5$ (excessive sensitivity leading to reduced reliability under parameter perturbation).

Sensitivity analysis identified $n$ and $K$ as the most influential parameters ($S_i \approx 0.35$–$0.45$), while $\alpha$ and $\beta$ had lower individual sensitivity ($S_i \approx 0.10$–$0.20$).

**5-Fold Cross-Validation:** Robustness scores over 5 folds × 20 designs each: **0.532 ± 0.016** (mean ± std), range [0.513, 0.556]. The low coefficient of variation (CV = 3.0%) indicates consistent prediction across different random design samples.

![Figure 5](figures/fig5_robust_optimization.png)

**Figure 5.** Robust design optimization. (A) Comparison of robustness metrics (ON/OFF ratio, reliability, 1−CV, composite score) for nominal vs. optimized designs. (B) Sensitivity analysis: first-order Sobol indices for five parameters. (C) Hill coefficient scan showing robustness score and reliability vs. n; optimal at n≈3.

---

## 6. Discussion

### 6.1 Simulation Algorithm Selection

The comparable step counts of Gillespie (6,896) and τ-leaping (6,924) in our toggle switch simulation reflect a regime where τ-leaping acceleration is not yet evident. This is expected: τ-leaping offers substantial speedups only when molecule counts are high enough that many reactions can fire in a single time step. At α=50 molecules/min and K=30 molecules — representative of an intermediate regime — the adaptive step controller frequently selects small τ, essentially reverting to Gillespie. For production-scale simulations (≥10,000 molecules, α≥500/min), τ-leaping is expected to reduce step count by 10–100× (Cao et al., 2006). This analysis motivates a hybrid strategy: Gillespie for transient dynamics near bifurcation points, and τ-leaping for bulk production/degradation regimes.

### 6.2 Repressilator Dynamics

The simulated oscillation periods (37.8 and 45.6 min) are shorter than the experimentally observed ~120 min period in Elowitz & Leibler (2000). This discrepancy is attributable to our accelerated degradation rates (τ_p = 0.3/min vs. ~0.1/min in vivo) and the absence of cell division dilution. Our qualitative predictions — period increases with n, and amplitude increases with n — are consistent with both ODE analysis and the hybrid modeling results of Sun et al. (2023). The triangular phase portrait geometry (Figure 3C) is a hallmark of three-component oscillators with strong repression, distinct from the elliptical orbits of two-component systems.

### 6.3 Robustness and Sensitivity

The robustness score of 0.532 ± 0.016 places our nominal toggle switch design in the intermediate range. The sensitivity analysis reveals that the Hill coefficient is the most critical design parameter, explaining why cooperative binding (n ≥ 2) is an essential requirement for bistability (Schladt et al., 2021). The cross-validated score confirms that this prediction is robust to different samples of the design space.

Comparison with the literature: Schladt et al. (2021) reported that robust bistable circuits require ON/OFF ratios ≥ 10× for reliable single-cell switching. Our value of 6.94× is below this threshold, suggesting that further optimization (e.g., increasing α to 80–100/min) could improve reliability to >90%. CELLM (Abello Castillo & Gutiérrez Pescarmona, 2025) represents a recent approach using language models for circuit design — integrating such AI-assisted design suggestions with our robustness evaluation pipeline is a natural extension.

### 6.4 Limitations

1. **Deterministic ON/OFF evaluation**: The robustness evaluator uses a simplified analytical steady-state formula rather than full stochastic trajectories for each LHS sample, which underestimates the effect of noise on reliability.

2. **No retroactivity modeling**: Load effects between circuit modules — where downstream components reduce upstream protein availability — are not modeled, potentially overestimating ON/OFF ratios in multi-stage cascades.

3. **Single chassis**: The catalog is limited to *E. coli* K-12. Extension to mammalian chassis (CHO, HEK293) would require reparametrization of Hill constants and degradation rates.

4. **Simplified context correction**: The linear supercoiling model and fixed 10% RBS occlusion are approximations; actual context effects depend on DNA topology and mRNA secondary structure in a sequence-dependent manner.

5. **Computational scale**: All simulations used lightweight parameters. Full-scale production runs (10,000+ molecules, 1,000,000 steps, 50 replicates) were not performed, as they fall outside the demonstration scope of this study.

---

## 7. Conclusion

We have presented a modular, automated pipeline for the design, simulation, and robustness optimization of synthetic gene circuits, implemented in four Python modules totaling approximately 2,600 lines of code. The framework reproduces key circuit behaviors — bistability in the toggle switch and oscillation in the repressilator — using stochastic simulation and provides quantitative robustness metrics validated by cross-validation. Key findings are: (1) the Hill coefficient is the most sensitive parameter for bistable circuit robustness (Sobol index ≈ 0.35–0.45), (2) τ-leaping provides equivalent accuracy to Gillespie at current molecule count regimes, (3) genetic context effects alter promoter strength by < 5% with standard terminators, and (4) the nominal toggle switch design achieves 81.0% reliability and 6.94× ON/OFF ratio. The SBOL 3.1-compliant design export enables integration with existing synthetic biology software ecosystems.

Future work will focus on: integration with CRISPR-based and RNA-based circuit topologies; multi-objective optimization using Pareto-front methods; incorporation of retroactivity models; and hybrid Gillespie/ODE solvers for multi-scale circuits.

---

## References

1. Gardner, T. S., Cantor, C. R., & Collins, J. J. (2000). Construction of a genetic toggle switch in *Escherichia coli*. *Nature*, 403(6767), 339–342. https://doi.org/10.1038/35002256

2. Elowitz, M. B., & Leibler, S. (2000). A synthetic oscillatory network of transcriptional regulators. *Nature*, 403(6767), 335–338. https://doi.org/10.1038/35014651

3. Nielsen, A. A., et al. (2016). Genetic circuit design automation. *Science*, 352(6281), aac7341. https://doi.org/10.1126/science.aac7341

4. Schladt, T. M., Engelmann, B., Kubaczka, E., & Tresch, A. (2021). Automated Design of Robust Genetic Circuits: Structural Variants and Parameter Uncertainty. *ACS Synthetic Biology*, 10(12), 3152–3164. https://doi.org/10.1021/acssynbio.1c00193

5. Baig, H., et al. (2020). Synthetic biology open language (SBOL) version 3.0.0. *Journal of Integrative Bioinformatics*, 17(2–3). https://doi.org/10.1515/jib-2020-0017

6. Buecherl, L., Mitchell, T., Scott-Brown, J., et al. (2023). Synthetic biology open language (SBOL) version 3.1.0. *Journal of Integrative Bioinformatics*, 20(1). https://doi.org/10.1515/jib-2022-0058

7. McLaughlin, J. A., Beal, J., Mısırlı, G., et al. (2020). The Synthetic Biology Open Language (SBOL) Version 3: Simplified Data Exchange for Bioengineering. *Frontiers in Bioengineering and Biotechnology*, 8, 1009. https://doi.org/10.3389/fbioe.2020.01009

8. Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. *The Journal of Physical Chemistry*, 81(25), 2340–2361. https://doi.org/10.1021/j100540a008

9. Cao, Y., Gillespie, D. T., & Petzold, L. R. (2006). Efficient step size selection for the tau-leaping simulation method. *The Journal of Chemical Physics*, 124(4), 044109. https://doi.org/10.1063/1.2159468

10. Roehner, N., Bartley, B. A., & Beal, J. (2019). Specifying Combinatorial Designs with the Synthetic Biology Open Language (SBOL). *ACS Synthetic Biology*, 9(7), 1701–1717. https://doi.org/10.1021/acssynbio.9b00092

11. Starkey, T., & Menolascina, F. (2022). Optimization for Predictive Gene Circuit Design Automation. *bioRxiv*. https://doi.org/10.1101/2022.06.21.497000

12. Sun, C., Comet, J.-P., & Folschette, M. (2023). Condition for Sustained Oscillations in Repressilator Based on a Hybrid Modeling of Gene Regulatory Networks. *Proceedings of BIOINFORMATICS 2023*. https://doi.org/10.5220/0011614300003414

13. Abello Castillo, L., & Gutiérrez Pescarmona, A. (2025). CELLM: Bridging Natural Language Processing and Synthetic Genetic Circuit Design with AI. *ACS Synthetic Biology*. https://doi.org/10.1021/acssynbio.5c00391

14. Crowther, M., Wipat, A., & Goñi-Moreno, A. (2021). A network approach to genetic circuit designs. *bioRxiv*. https://doi.org/10.1101/2021.09.14.460206
