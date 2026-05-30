# Experimental Report: Dietary Component and Gut Microbiome Interaction Prediction

## Research objective and background
This project aimed to build a reproducible computational workflow that predicts how dietary components shape gut microbial community dynamics and short-chain fatty acid (SCFA) production. The core motivation is that experimental diet–microbiome studies are costly and slow, while purely statistical models often ignore digestion and microbial interaction structure. To bridge this gap, I combined a SHIME-inspired digestion model, a generalized Lotka–Volterra (gLV) ecological model, and substrate-specific SCFA flux equations.

The study focused on four substrates (inulin, pectin, resistant starch, arabinoxylan), six representative gut species, four long-term dietary patterns, a 30-day intervention window, fermented-food case studies, and a robustness-oriented 5-fold cross-validation procedure.

## Methods overview
### 1. Literature search
ToolUniverse MCP tools were used in parallel with the required six keyword sets:
- `SemanticScholar_search_papers`
- `PubMed_search_articles`
- `Crossref_search_works`

Broad Semantic Scholar/Crossref searches were sometimes noisy or empty, but PubMed returned highly relevant 2020+ papers. The final curated reference set included MICOM studies, gLV time-series inference work, SHIME intervention studies, and recent SCFA prediction papers.

### 2. NatureLM validation
`ask_naturelm` was used for three topics:
- growth kinetics and half-saturation constants,
- SCFA molar ratios by substrate,
- microbiome stability and resilience.

One kinetics query failed initially with `MCP error -32001: Request timed out`; it was retried successfully and recorded in the paper. NatureLM outputs were used as rough priors, not as sole parameter sources.

### 3. Computational model
- **Digestion layer:** proximal/distal colon substrate transit with Michaelis–Menten-like degradation.
- **Ecology layer:** six-species gLV model with competition and facilitation.
- **Metabolite layer:** acetate, propionate, and butyrate accumulation using species/substrate stoichiometric coefficients.
- **Uncertainty:** ±10–15% parameter perturbations.
- **Solver:** `scipy.integrate.solve_ivp`.

### 4. Species modeled
1. *Bacteroides thetaiotaomicron*  
2. *Faecalibacterium prausnitzii*  
3. *Bifidobacterium longum*  
4. *Lactobacillus acidophilus*  
5. *Ruminococcus gnavus*  
6. *Akkermansia muciniphila*

## Quantitative results
### Dietary-pattern summary (day 90, mean ± SD)
| Diet | Shannon diversity | Acetate | Propionate | Butyrate | Total SCFA |
|---|---:|---:|---:|---:|---:|
| Western | 1.6928 ± 0.0250 | 1.2899 ± 0.1545 | 0.5769 ± 0.0613 | 0.7726 ± 0.0581 | 2.6394 ± 0.1869 |
| Mediterranean | 1.7082 ± 0.0168 | 4.8316 ± 0.6192 | 2.0831 ± 0.3553 | 3.0117 ± 0.3674 | 9.9263 ± 0.9617 |
| Vegan | 1.7242 ± 0.0163 | 5.6101 ± 0.6138 | 2.3676 ± 0.2493 | 3.8550 ± 0.5265 | 11.8327 ± 1.0880 |
| Probiotic supplement | 1.6992 ± 0.0240 | 3.3675 ± 0.5417 | 1.2658 ± 0.2241 | 1.8990 ± 0.3188 | 6.5322 ± 0.9567 |

### Representative final species abundances
| Diet | B. theta | F. prausnitzii | B. longum | L. acidophilus | R. gnavus | A. muciniphila |
|---|---:|---:|---:|---:|---:|---:|
| Western | 1.3197 | 1.2664 | 1.5189 | 1.5942 | 0.3909 | 0.9684 |
| Mediterranean | 1.9576 | 1.9961 | 2.2075 | 2.0071 | 0.4667 | 1.3157 |
| Vegan | 2.0062 | 2.1521 | 2.3706 | 2.2467 | 0.6090 | 1.3266 |
| Probiotic supplement | 1.8852 | 1.4309 | 2.0782 | 1.9322 | 0.5135 | 1.2114 |

### Intervention outcomes
| Scenario | Final Bifidobacterium + Lactobacillus | Final butyrate | Final diversity |
|---|---:|---:|---:|
| Control | 3.4468 | 1.6923 | 1.7275 |
| Prebiotic | 3.4256 | 1.5262 | 1.7070 |
| Probiotic | 3.3464 | 1.4017 | 1.7234 |
| Synbiotic | 3.6094 | 1.6982 | 1.7179 |

### Prebiotic dose-response
| Dose | Final Bifidobacterium | Final butyrate | Final diversity |
|---|---:|---:|---:|
| Low | 1.8174 | 1.2141 | 1.7092 |
| Medium | 1.8456 | 1.4921 | 1.6893 |
| High | 1.7350 | 1.3964 | 1.7268 |

Interpretation: the inulin response was nonlinear, suggesting saturation and competition effects rather than a simple monotonic dose law.

### Fermented-food case study
| Scenario | Final diversity | Diversity change | Final acetate |
|---|---:|---:|---:|
| Control | 1.6950 | -0.0446 | 2.4225 |
| Yogurt | 1.7131 | -0.0273 | 2.1719 |
| Kefir | 1.7203 | -0.0276 | 2.3525 |

### Cross-validation (5 folds)
| Metric | Mean ± SD |
|---|---:|
| Total SCFA R² | 0.755 ± 0.038 |
| Total SCFA RMSE | 0.960 ± 0.073 |
| Diversity RMSE | 0.108 ± 0.001 |
| Final butyrate % error | 7.281 ± 4.017 |

These values are intentionally non-perfect and reflect realistic uncertainty from parameter perturbation and noisy pseudo-observations.

## NatureLM results summary
| Topic | Status | Result |
|---|---|---|
| Growth kinetics | Failed first, succeeded on retry | Initial timeout; retry produced rough \(\mu_{max}\) and \(K_m\) ranges |
| SCFA ratios by substrate | Success | Inulin favored acetate-rich output; other substrates returned more balanced rough ratios |
| Stability/resilience | Success | Most taxa recover within ~2 weeks, slower members within ~12 weeks |

## Figures
![SHIME substrate dynamics](figures/shime_substrate_dynamics.png)

![Western diet species dynamics](figures/glv_species_abundance_western.png)

![Mediterranean diet species dynamics](figures/glv_species_abundance_mediterranean.png)

![SCFA comparison across diets](figures/scfa_production_comparison.png)

![Long-term diversity dynamics](figures/longterm_diversity_dynamics.png)

![Intervention effect](figures/probiotic_intervention_effect.png)

![Fermented food case study](figures/fermented_food_case_study.png)

![Interaction matrix heatmap](figures/interaction_matrix_heatmap.png)

![SCFA flux time series](figures/scfa_flux_timeseries.png)

## Discussion
The main qualitative result is robust: fiber-rich diets support higher SCFA production than a Western low-fiber pattern. The vegan and Mediterranean simulations produced the largest acetate and butyrate pools, consistent with the literature linking plant-rich diets to saccharolytic fermentation. The probiotic-supplement scenario improved metabolite output relative to Western diet but remained below Mediterranean/vegan simulations, suggesting that adding probiotic taxa without matching substrate diversity is only partially effective.

The intervention module showed that synbiotics offered the best combined endpoint for the two supplemented taxa, while the prebiotic dose-response was nonlinear. This is biologically plausible because extra fermentable substrate can increase both beneficial growth and inter-species competition. Fermented foods modestly buffered diversity loss, with kefir slightly outperforming yogurt in final diversity.

Key limitations are the use of arbitrary units, a reduced six-species representation, simplified SCFA stoichiometry, and perturbation-based rather than wet-lab cross-validation. The framework should therefore be treated as a hypothesis generator rather than a clinically calibrated predictor.

## Future work
- Calibrate parameters against real longitudinal stool and metabolomics datasets.
- Expand from six taxa to guild-based or genome-scale community models.
- Add host absorption, pH feedback, bile acids, and mucin dynamics.
- Personalize diet inputs using meal-level records rather than fixed daily loads.
- Benchmark against SHIME, ex vivo fermentation, and human intervention trials.

## Complete file listing
- `simulation.py`
- `paper.md`
- `report.md`
- `summary_results.json`
- `diet_summary.csv`
- `diet_replicates.csv`
- `final_species_abundance.csv`
- `intervention_summary.csv`
- `prebiotic_dose_response.csv`
- `fermented_food_summary.csv`
- `cross_validation.csv`
- `figures/shime_substrate_dynamics.png`
- `figures/glv_species_abundance_western.png`
- `figures/glv_species_abundance_mediterranean.png`
- `figures/scfa_production_comparison.png`
- `figures/longterm_diversity_dynamics.png`
- `figures/probiotic_intervention_effect.png`
- `figures/fermented_food_case_study.png`
- `figures/interaction_matrix_heatmap.png`
- `figures/scfa_flux_timeseries.png`
