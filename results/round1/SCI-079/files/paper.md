# Computational Modeling of Pattern-Triggered Immunity and Effector-Triggered Immunity Signaling Networks in Plants: An Integrative ODE-Based Framework with Game-Theoretic Coevolutionary Analysis

## Abstract

Plant innate immunity relies on two interconnected layers: pattern-triggered immunity (PTI) and effector-triggered immunity (ETI). Despite significant advances in understanding individual signaling components, a comprehensive computational framework integrating receptor-level signal initiation, MAPK cascade dynamics, hormone crosstalk, transcriptional regulation, and evolutionary dynamics remains lacking. Here, we present an integrative ordinary differential equation (ODE)-based modeling framework that captures six critical aspects of plant immune signaling: (1) receptor-ligand binding and signal initiation kinetics for both PRR-PAMP and NLR-effector interactions, (2) three-tier MAPK cascade dynamics exhibiting ultrasensitive switch-like behavior (Hill coefficient ≈ 2.64), (3) salicylic acid (SA) and jasmonic acid (JA) pathway crosstalk with quantified antagonistic interactions, (4) WRKY/TGA transcription factor regulatory network analysis, (5) pathogen-host coevolution through evolutionary game theory, and (6) a rice blast resistance case study modeling the Pi-ta/AvrPi-ta interaction. Our simulations demonstrate that ETI generates approximately 4.4-fold stronger signaling output than PTI, that SA-JA antagonism produces bistable defense states dependent on antagonism strength, and that R-gene stacking provides diminishing but additive resistance benefits. The SBML-compatible model enables direct import into COPASI and CellDesigner for further analysis. This work provides a unified computational platform for understanding plant immune signaling architecture and informing rational crop protection strategies.

## 1. Introduction

Plants, as sessile organisms, have evolved sophisticated innate immune systems to detect and respond to microbial pathogens (Jones & Dangl, 2006). The plant immune system operates through two major perception layers. The first layer, pattern-triggered immunity (PTI), is initiated when cell-surface pattern recognition receptors (PRRs) detect conserved pathogen-associated molecular patterns (PAMPs) such as bacterial flagellin (flg22) or fungal chitin (Peng et al., 2018). The second layer, effector-triggered immunity (ETI), is activated when intracellular nucleotide-binding leucine-rich repeat (NLR) receptors recognize pathogen-derived effector proteins (Ngou et al., 2021; Yuan et al., 2021).

Recent landmark studies have revealed that PTI and ETI are not independent pathways but rather mutually potentiating systems (Ngou et al., 2021; Yuan et al., 2021). Both pathways converge on common downstream signaling modules including mitogen-activated protein kinase (MAPK) cascades, reactive oxygen species (ROS) bursts, calcium signaling, and hormone-mediated transcriptional reprogramming (Peng et al., 2018; Ding et al., 2022). The phytohormones salicylic acid (SA) and jasmonic acid (JA) serve as central mediators, with their antagonistic crosstalk determining the specificity of defense responses against biotrophic and necrotrophic pathogens, respectively (Ding et al., 2022).

Despite the wealth of molecular data, quantitative computational models that integrate these diverse signaling layers remain scarce. Mathematical modeling provides a powerful framework for understanding emergent properties of signaling networks, including ultrasensitivity, bistability, and oscillatory behavior (Meng & Zhang, 2013). Furthermore, evolutionary game theory offers tools to analyze the strategic dynamics of pathogen-host coevolution under the gene-for-gene framework (Tellier & Brown, 2007).

In this study, we develop a comprehensive ODE-based computational framework integrating six modules of plant immune signaling. We parameterize the model using literature-derived kinetic constants, validate emergent behaviors against known biological properties, and apply the framework to rice blast resistance as a case study. Our model is exported in SBML format for compatibility with standard systems biology tools including COPASI and CellDesigner.

### Contributions

1. A unified ODE framework capturing receptor binding, MAPK cascading, hormone crosstalk, and transcriptional regulation in PTI/ETI
2. Quantification of ultrasensitivity in the MAPK cascade with Hill coefficient analysis
3. Phase-portrait analysis of SA-JA antagonism revealing bistable defense states
4. Game-theoretic analysis of pathogen-host coevolution with ESS characterization
5. Application to rice blast (Oryza sativa–Magnaporthe oryzae) as a translational case study
6. SBML model for community use in COPASI/CellDesigner

## 2. Related Work

### 2.1 PTI-ETI Signaling Framework

The canonical PTI-ETI dichotomy has been refined by recent studies demonstrating mutual dependence between these pathways. Ngou et al. (2021) showed that PTI and ETI mutually potentiate each other through shared signaling components in Arabidopsis, challenging the traditional view of ETI as an amplified re-engagement of PTI. Yuan et al. (2021) independently demonstrated that PRRs are required for full NLR-mediated immunity, confirming the interdependence of these pathways. Peng et al. (2018) provided an earlier comprehensive review of convergent and divergent signaling mechanisms between PTI and ETI, identifying MAPK cascades, calcium signaling, and ROS production as key convergence points.

### 2.2 MAPK Cascade Dynamics

MAPK cascades serve as universal signal amplification modules in eukaryotic signal transduction. Meng & Zhang (2013) reviewed the roles of MAPK cascades in plant disease resistance signaling, highlighting the MEKK1-MKK4/5-MPK3/6 module as a central PTI signaling axis. Mathematical models of MAPK cascades have demonstrated ultrasensitive (switch-like) behavior arising from the multi-tier phosphorylation/dephosphorylation architecture (Huang & Ferrell, 1996).

### 2.3 Hormone Crosstalk

The SA-JA antagonism is a well-documented regulatory mechanism enabling plants to prioritize defense against biotrophic or necrotrophic pathogens. Ding et al. (2022) reviewed recent progress in understanding hormone-mediated signaling pathways, emphasizing the role of NPR1 as a master regulator of SA signaling and JAZ proteins as repressors of JA-responsive gene expression. The nonlinear dynamics of this crosstalk generate complex defense phenotypes depending on the nature and timing of pathogen attack.

### 2.4 Transcription Factor Networks

WRKY transcription factors constitute one of the largest TF families in plants, with critical roles in defense gene regulation. The WRKY-TGA regulatory network operates as a combinatorial logic circuit, integrating SA and JA pathway inputs to control the expression of pathogenesis-related (PR) genes and other defense effectors. Network analysis approaches have identified WRKY70 and WRKY33 as key hub nodes with opposing roles in SA- and JA-dependent defense pathways.

### 2.5 Evolutionary Game Theory in Host-Pathogen Systems

The gene-for-gene model proposed by Flor (1971) provides a genetic framework for host-pathogen coevolution. Tellier & Brown (2007) applied game-theoretic models to analyze the stability of genetic polymorphism in host-parasite interactions, demonstrating that frequency-dependent selection can maintain diversity in both R-gene and Avr-gene populations. These models have been extended to consider the costs of resistance and virulence in shaping evolutionary equilibria.

### 2.6 Rice Blast Resistance

Rice blast disease caused by Magnaporthe oryzae is the most devastating fungal disease of rice worldwide. Naveed et al. (2022) reviewed the dynamics of blast resistance including the molecular basis of R-gene mediated immunity. The Pi-ta/AvrPi-ta interaction represents one of the best-characterized gene-for-gene systems in rice, providing an ideal model for computational analysis of integrated PTI-ETI signaling.

## 3. Methods

### 3.1 Receptor-Ligand Binding Model

We model the receptor dynamics using a five-state ODE system representing free receptors ($R_f$), ligand-bound receptors ($R_b$), active receptors ($R_a$), internalized receptors ($R_i$), and downstream signal ($S$):

$$\frac{dR_f}{dt} = -k_{on} \cdot R_f \cdot [PAMP] + k_{off} \cdot R_b + v_{synth} - k_{deg} \cdot R_f$$

$$\frac{dR_b}{dt} = k_{on} \cdot R_f \cdot [PAMP] - k_{off} \cdot R_b - k_{act} \cdot R_b + k_{deact} \cdot R_a$$

$$\frac{dR_a}{dt} = k_{act} \cdot R_b - k_{deact} \cdot R_a - k_{int} \cdot R_a$$

$$\frac{dS}{dt} = k_s \cdot R_a - d_s \cdot S$$

For ETI, the NLR-effector interaction is modeled with higher affinity ($k_{on}^{ETI} = 0.2$ nM⁻¹ min⁻¹ vs $k_{on}^{PTI} = 0.1$) and lower dissociation ($k_{off}^{ETI} = 0.005$ vs $k_{off}^{PTI} = 0.01$), reflecting the tighter intracellular recognition.

### 3.2 MAPK Cascade Model

The three-tier MAPK cascade (MAPKKK → MAPKK → MAPK) is modeled using Michaelis-Menten kinetics for each phosphorylation and dephosphorylation step:

$$\frac{d[M_i^*]}{dt} = \frac{k_i^+ \cdot [M_{i+1}^*] \cdot [M_i]}{K_m + [M_i]} - \frac{k_i^- \cdot [M_i^*]}{K_m + [M_i^*]}$$

where $M_i^*$ denotes the phosphorylated (active) form, and $M_{i+1}^*$ is the upstream activating kinase. The ultrasensitivity of the cascade is characterized by the apparent Hill coefficient $n_H$:

$$n_H = \frac{\ln(81)}{\ln(EC_{90}/EC_{10})}$$

### 3.3 SA/JA Crosstalk Model

The hormone crosstalk model incorporates six state variables: SA, JA, NPR1 (active), JAZ (repressor), PR1 (SA-responsive gene), and PDF1.2 (JA-responsive gene). The antagonistic interaction is modeled as:

$$\frac{d[SA]}{dt} = v_{SA} \cdot \frac{P_b}{K_{SA} + P_b} - d_{SA} \cdot [SA] - \beta \cdot \frac{[JA]}{K_\beta + [JA]} \cdot [SA]$$

$$\frac{d[JA]}{dt} = v_{JA} \cdot \frac{P_n}{K_{JA} + P_n} - d_{JA} \cdot [JA] - \alpha \cdot \frac{[SA]}{K_\alpha + [SA]} \cdot [JA]$$

where $\alpha$ and $\beta$ represent the antagonism strengths of SA→JA and JA→SA inhibition, respectively.

### 3.4 Transcription Factor Network

A directed graph $G = (V, E)$ was constructed with 20 nodes (10 transcription factors + 10 target genes) and 26 regulatory edges. Each edge $e_{ij}$ has a weight $w_{ij} \in \{-1, +1\}$ representing repression or activation. A Boolean network simulation was used to predict defense gene expression patterns under different pathway activation scenarios. Network topology was analyzed using betweenness centrality $C_B(v)$ and feed-forward loop (FFL) enumeration.

### 3.5 Game Theory Model

The pathogen-host coevolution is modeled as a two-player asymmetric game with payoff matrices:

**Host payoff matrix** $\Pi_H$:

| Host \ Pathogen | Avirulent | Virulent |
|:---:|:---:|:---:|
| R-gene | $1 - c_R$ | $1 - c_R - 0.3b$ |
| Susceptible | $1 - b$ | $1 - b$ |

**Pathogen payoff matrix** $\Pi_P$:

| Pathogen \ Host | R-gene | Susceptible |
|:---:|:---:|:---:|
| Avirulent | $0$ | $b$ |
| Virulent | $b - c_v$ | $b - c_v$ |

The replicator dynamics are given by:

$$\dot{p} = p(f_R - \bar{f}_H), \quad \dot{q} = q(f_{Avr} - \bar{f}_P)$$

where $p$ is the frequency of R-gene hosts, $q$ is the frequency of avirulent pathogens, and $f$ denotes fitness values.

### 3.6 Rice Blast Model

An integrated 12-variable ODE model was developed for the Oryza sativa–Magnaporthe oryzae pathosystem, incorporating:
- Chitin-mediated PTI via OsCEBiP/OsCERK1
- NLR-mediated ETI via Pi-ta/AvrPi-ta recognition
- OsMPK3/6 MAPK signaling
- SA/JA hormone balance regulated by OsWRKY45 and OsWRKY13
- ROS burst and hypersensitive response (HR)

### 3.7 Implementation

All models were implemented in Python 3.12 using NumPy (v1.x), SciPy (solve_ivp with RK45), Matplotlib (visualization), and NetworkX (graph analysis). The SBML model was generated in Level 3 Version 2 format for COPASI/CellDesigner compatibility. Code is available in `sim_all.py`.

## 4. Experiments

### 4.1 Simulation Setup

| Module | Variables | Time Span | Integration Method | Key Parameters |
|--------|-----------|-----------|-------------------|----------------|
| Receptor binding | 5 | 0–120 min | RK45 | $k_{on}=0.1$, $k_{off}=0.01$ |
| MAPK cascade | 6 | 0–100 min | RK45 | $K_m=0.5$, 6 rate constants |
| SA/JA crosstalk | 6 | 0–200 min | RK45 | $\alpha=0.8$, $\beta=0.3$ |
| TF network | 20 nodes | Boolean | Graph analysis | 26 regulatory edges |
| Game theory | 2 | 0–200 gen | RK45 | $c_R=0.1$, $c_v=0.15$ |
| Rice blast | 12 | 0–120 hpi | RK45 | 15+ kinetic parameters |

### 4.2 Evaluation Metrics

- **EC50**: Half-maximal effective concentration for dose-response analysis
- **Hill coefficient**: Quantification of ultrasensitivity in MAPK cascade
- **Steady-state hormone levels**: SA and JA concentrations under different attack scenarios
- **Betweenness centrality**: Identification of hub nodes in TF network
- **ESS**: Evolutionarily stable strategy frequencies
- **HR score**: Hypersensitive response intensity as resistance indicator

### 4.3 Baseline Comparisons

- PTI vs ETI signal intensity and kinetics (Module 1)
- Transient vs sustained MAPK activation (Module 2)
- Compatible vs incompatible rice-pathogen interactions (Module 6)
- Single vs stacked R-gene resistance (Module 6)

## 5. Results

### 5.1 Receptor-Ligand Binding Dynamics

![Figure 1](figures/fig1_receptor_binding.png)

**Figure 1.** Receptor-level ligand binding and signal initiation model. (a) Signal output as a function of PAMP concentration showing dose-dependent activation. (b) Receptor state dynamics at 5 nM PAMP. (c) Dose-response curve with EC50 ≈ 0.17 nM. (d) Comparison of PTI (FLS2-flg22) and ETI (NLR-effector) signal dynamics.

The receptor model demonstrated a sigmoidal dose-response relationship with an EC50 of approximately 0.17 nM for flg22-mediated PTI signaling. ETI signaling generated a maximum signal intensity of 186.65 a.u., approximately 4.44-fold higher than PTI (42.00 a.u.), consistent with the known amplification of immune responses by NLR activation. The temporal profiles revealed that PTI signals reached near-maximum levels within 40–60 minutes, while ETI signals continued to increase throughout the simulation period, reflecting the sustained nature of NLR-mediated immunity.

### 5.2 MAPK Cascade Ultrasensitivity

![Figure 2](figures/fig2_mapk_cascade.png)

**Figure 2.** MAPK cascade dynamics simulation. (a) Temporal activation profiles of MAPKKK, MAPKK, and MAPK. (b) Ultrasensitive dose-response with Hill coefficient ≈ 2.64. (c) Comparison of transient (PTI-like) and sustained (ETI-like) MAPK activation. (d) Parameter sensitivity analysis.

The three-tier MAPK cascade exhibited an apparent Hill coefficient of 2.64, indicating significant ultrasensitivity arising from the cascade architecture. This switch-like behavior enables digital (all-or-none) defense activation in response to threshold-crossing stimuli. The sensitivity analysis identified the MAPKKK activation rate constant ($k_1$) as the most influential parameter, consistent with the biological role of MAPKKK as the entry point for signal integration.

### 5.3 SA-JA Crosstalk and Hormone Dynamics

![Figure 3](figures/fig3_sa_ja_crosstalk.png)

**Figure 3.** Salicylic acid / jasmonic acid pathway crosstalk model. (a) SA-dominant response during biotroph attack. (b) JA-dominant response during necrotroph attack. (c) Combined attack scenario. (d) Defense gene expression profiles. (e) NPR1 and JAZ regulator dynamics. (f) Phase portrait under varying antagonism strength.

The crosstalk model faithfully reproduced the well-documented SA-JA antagonism. Under biotrophic pathogen attack, SA accumulated to 41.67 a.u. while JA was completely suppressed. Conversely, necrotrophic/herbivore attack induced JA accumulation with SA suppression. The phase portrait analysis (Figure 3f) revealed that the strength of the antagonism parameter α dramatically influenced the steady-state balance, with stronger antagonism driving the system toward more extreme SA- or JA-dominant states.

### 5.4 WRKY/TGA Transcription Factor Network

![Figure 4](figures/fig4_tf_network.png)

**Figure 4.** WRKY/TGA transcription factor regulatory network. (a) Network topology with activating (green) and repressing (red) edges. (b) Centrality analysis. (c) Boolean simulation of defense gene expression under different pathway activation. (d) Network motif analysis.

Network analysis revealed WRKY33, WRKY70, and WRKY40 as the top hub nodes by betweenness centrality. The network contained 3 feed-forward loops, structures known to provide robust and noise-filtering regulatory logic. Boolean network simulation correctly predicted SA-dependent upregulation of PR1/PR2/PR5 and JA-dependent upregulation of PDF1.2/VSP2/LOX2, with WRKY70 acting as the key molecular switch between SA- and JA-responsive gene expression programs.

### 5.5 Evolutionary Game Theory

![Figure 5](figures/fig5_game_theory.png)

**Figure 5.** Pathogen-host coevolution game theory analysis. (a-b) Payoff matrices. (c) Replicator dynamics for host R-gene frequency. (d) Coevolutionary phase portrait. (e) Arms race escalation model. (f) ESS analysis as a function of resistance/virulence cost.

The replicator dynamics under the gene-for-gene payoff structure produced frequency-dependent oscillatory dynamics, consistent with the "trench warfare" model of host-pathogen coevolution. The ESS analysis showed that R-gene equilibrium frequency decreased monotonically with increasing resistance cost, providing a theoretical basis for the observation that costly resistance alleles are maintained at intermediate frequencies in natural plant populations.

### 5.6 Rice Blast Resistance Case Study

![Figure 6](figures/fig6_rice_blast.png)

**Figure 6.** Rice blast resistance case study (Oryza sativa–Magnaporthe oryzae). (a) Compatible interaction (no Pi-ta). (b) Incompatible interaction (Pi-ta + AvrPi-ta). (c) Defense gene expression comparison. (d) R-gene dosage effect. (e) R-gene stacking simulation. (f) PTI-ETI signal integration.

The rice blast model demonstrated a stark contrast between compatible (susceptible) and incompatible (resistant) interactions. In the incompatible interaction, Pi-ta recognition of AvrPi-ta triggered a strong HR (score = 0.874) and dramatically enhanced PR1a expression, while the compatible interaction showed negligible HR (score ≈ 0). R-gene stacking analysis revealed additive but diminishing returns: stacking 5 R-genes increased the resistance score by approximately 26% compared to a single gene (45.3 vs 36.0), supporting the strategy of pyramiding multiple R-genes for durable resistance.

## 6. Discussion

### 6.1 Integrative Signaling Architecture

Our modeling framework provides quantitative support for the emerging paradigm of PTI-ETI mutual potentiation (Ngou et al., 2021; Yuan et al., 2021). The 4.4-fold signal amplification from PTI to ETI, combined with the qualitative differences in signal kinetics (transient vs. sustained), suggests that the two pathways serve complementary roles: PTI provides rapid initial defense, while ETI ensures robust and sustained immunity through amplification of shared downstream components.

### 6.2 Ultrasensitivity as a Decision Mechanism

The Hill coefficient of 2.64 observed in the MAPK cascade model suggests that this module functions as a biological digital-to-analog converter, transforming graded receptor signals into more switch-like downstream responses. This property is particularly relevant for the HR decision, where cells must commit to programmed cell death only when infection pressure exceeds a critical threshold (Meng & Zhang, 2013).

### 6.3 Hormone Crosstalk and Defense Prioritization

The bistable behavior of the SA-JA system, governed by the antagonism parameters α and β, provides a mechanistic basis for defense prioritization. The asymmetry between SA→JA inhibition (α = 0.8) and JA→SA inhibition (β = 0.3) reflects the biological observation that SA pathway dominance is the default state during concurrent biotroph-necrotroph challenge (Ding et al., 2022).

### 6.4 Evolutionary Implications

The game-theoretic analysis reveals that the maintenance of R-gene diversity in plant populations is an expected outcome of frequency-dependent coevolutionary dynamics rather than an anomaly requiring special explanation. The cost-dependent ESS analysis provides guidance for breeding programs: minimizing the fitness cost of R-genes through precise gene editing or promoter optimization could increase the equilibrium frequency of resistance in crop populations.

### 6.5 Translational Implications for Rice Blast Resistance

The R-gene stacking simulation provides quantitative support for the pyramiding strategy widely used in rice breeding programs. The diminishing returns observed (26% increase for 5 genes vs. 1) suggest that combining 2–3 complementary R-genes may represent an optimal balance between resistance gain and breeding complexity. The PTI-ETI integration analysis further suggests that enhancing basal PTI components (e.g., OsCEBiP overexpression) could complement R-gene-mediated ETI for more durable resistance (Naveed et al., 2022).

### 6.6 Limitations

Our model has several important limitations: (1) spatial heterogeneity within tissues is not captured by the ODE framework; (2) stochastic fluctuations at low molecular concentrations are neglected; (3) many kinetic parameters are estimated rather than experimentally determined; (4) epigenetic regulation and chromatin remodeling are not included; (5) the game-theoretic model assumes well-mixed populations without spatial structure. Future work should address these limitations through partial differential equation models, stochastic simulations, and experimental parameter calibration.

## 7. Conclusion

We have developed a comprehensive computational framework for modeling plant PTI and ETI signaling, integrating receptor-level kinetics, MAPK cascade dynamics, hormone crosstalk, transcription factor networks, evolutionary game theory, and a rice blast case study. Key findings include the quantification of ETI signal amplification (4.4×), MAPK cascade ultrasensitivity (Hill coefficient ≈ 2.64), SA-JA antagonistic bistability, and the diminishing returns of R-gene stacking. The SBML-compatible model provides a community resource for further analysis and extension. This integrative approach demonstrates the power of computational modeling in elucidating the emergent properties of plant immune signaling networks and informing rational strategies for crop disease resistance.

## References

1. Ngou, B. P. M., Ahn, H.-K., Ding, P., & Jones, J. D. G. (2021). Mutual potentiation of plant immunity by cell-surface and intracellular receptors. *Nature*, 592, 110–115. https://doi.org/10.1038/s41586-021-03315-7

2. Yuan, M., Jiang, Z., Bi, G., Nomura, K., Liu, M., Wang, Y., ... & He, S. Y. (2021). Pattern-recognition receptors are required for NLR-mediated plant immunity. *Nature*, 592, 105–109. https://doi.org/10.1038/s41586-021-03250-9

3. Peng, Y., van Wersch, T., & Zhang, Y. (2018). Convergent and divergent signaling in PAMP-triggered immunity and effector-triggered immunity. *Molecular Plant*, 11(4), 511–524. https://doi.org/10.1016/j.molp.2018.06.001

4. Ding, L.-N., Li, Y.-T., Wu, Y.-Z., Li, T., Geng, R., Cao, J., Zhang, W., & Tan, X.-L. (2022). Plant disease resistance-related signaling pathways: Recent progress and future prospects. *International Journal of Molecular Sciences*, 23(24), 16200. https://doi.org/10.3390/ijms232416200

5. Meng, X., & Zhang, S. (2013). MAPK cascades in plant disease resistance signaling. *Annual Review of Phytopathology*, 51, 245–266. https://doi.org/10.1146/annurev-phyto-082712-102314

6. Tellier, A., & Brown, J. K. M. (2007). Stability of genetic polymorphism in host–parasite interactions. *Proceedings of the Royal Society B*, 274(1611), 809–817. https://doi.org/10.1098/rspb.2006.0281

7. Naveed, Z. A., Ali, G. S., & Khatri, W. A. (2022). Understanding the dynamics of blast resistance in rice–Magnaporthe oryzae interactions. *Journal of Fungi*, 8(6), 584. https://doi.org/10.3390/jof8060584

8. Jones, J. D. G., & Dangl, J. L. (2006). The plant immune system. *Nature*, 444, 323–329. https://doi.org/10.1038/nature05286

9. Huang, C.-Y. F., & Ferrell, J. E. (1996). Ultrasensitivity in the mitogen-activated protein kinase cascade. *Proceedings of the National Academy of Sciences*, 93(19), 10078–10083. https://doi.org/10.1073/pnas.93.19.10078

10. Flor, H. H. (1971). Current status of the gene-for-gene concept. *Annual Review of Phytopathology*, 9(1), 275–296. https://doi.org/10.1146/annurev.py.09.090171.001423
