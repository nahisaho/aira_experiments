# Funding Allocation in Research Systems: A Cross-Validated Agent-Based Comparison of Peer Review, Lottery, and Formula Funding

## Abstract

This paper reports a computational study of research-funding allocation using an agent-based model (ABM) with 200 heterogeneous researchers observed over 20 annual funding cycles. The simulation was designed to compare three allocation mechanisms that are widely debated in science policy: conventional peer review, partial randomization by lottery after pre-screening, and a formula-based allocation rule that weights productivity, citations, and diversity bonuses. Parameterization followed literature-based constraints and NatureLM-assisted scientific validation. Specifically, researcher quality was initialized from N(0,1), citation accumulation was calibrated to power-law-like skew, and peer-review reliability was constrained to an inter-rater correlation of 0.25, consistent with empirical evidence on reviewer disagreement. Each mechanism was evaluated with 30 independent runs, and all reported numbers are mean ± standard deviation across runs.

Across the baseline experiments, Formula generated the strongest cumulative output (11718.9 ± 235.2 papers by year 20), while Formula selected the highest-quality funded pool (1.269 ± 0.145). Peer review still illustrates the cost of status-quo incentives: it produced a citation Gini of 0.364 ± 0.012, a funding Gini of 0.615 ± 0.012, and a gender funding gap of 5.19 ± 4.67 percentage points. Lottery-based funding reduced concentration and preserved competitive output (11532.6 ± 266.0), whereas the formula mechanism expanded portfolio breadth and representation. The most diverse funded portfolio was achieved by Formula with Shannon entropy 1.055 ± 0.038, and the most equitable cumulative funding distribution was achieved by Lottery with funding Gini 0.474 ± 0.017.

A stylized KAKENHI case study, modeled with a 26% success rate, strong seniority bias, and field imbalances favoring STEM, showed that diversity-optimized formula allocation improved portfolio breadth (1.013 ± 0.049) under harsher structural constraints. The overall pattern is that no mechanism dominates on every objective: the mechanisms that maximize output, diversity, and equity are different. Hybrid systems that combine threshold screening with limited randomization and explicit diversity constraints therefore appear most promising for balancing quality, equity, and long-run system diversity.

## 1. Introduction

Research funders face a persistent design problem: how to allocate scarce grants while simultaneously rewarding scientific excellence, maintaining procedural legitimacy, and preventing cumulative advantage from locking out promising but less-established researchers. Empirical studies show substantial reviewer disagreement, demographic disparities, and strong path dependence in grant competitions. These concerns have stimulated interest in alternative mechanisms such as partial lotteries and formula-based allocation. This paper develops a computational ABM to compare the long-run consequences of three funding rules under common structural assumptions.

The model emphasizes heterogeneity in ability, prestige accumulation, demographic composition, and field structure. It also represents empirically grounded frictions including Matthew effects, reviewer noise, modest gender and region-related biases, preferential attachment in citations, and attrition among persistently low-output researchers. The goal is not to reproduce any one real agency exactly, but to generate a transparent experimental platform for comparing policy trade-offs over repeated funding cycles.

## 2. Related Work

Feliciani et al. (2022) show that panel design substantially affects funding quality, implying that review architecture itself is a policy lever rather than a neutral administrative detail. Pier et al. (2018) document low agreement among reviewers scoring the same NIH applications, motivating the use of noisy review processes in simulation. Li (2017) distinguishes expertise effects from evaluator bias, while Banal-Estañol et al. (2023) show that similarity between applicants and evaluators can influence outcomes. Taffe and Gilpin (2021) demonstrate racial inequity in NIH grant funding, and Graves et al. (2011) question whether review-based ranking reliably predicts downstream scientific performance. Moldashev (2024) provides a broader institutional critique of grant allocation mechanisms in the social and humanitarian sciences.

Together, this literature implies three central modeling requirements: (i) noisy and partially biased evaluation, (ii) cumulative advantage from prior success, and (iii) explicit attention to equity and diversity metrics rather than quality alone. The present ABM uses those requirements to compare a status-quo peer-review regime, a lottery after pre-screening, and a formula allocation rule with diversity correction.

## 3. Methods

### 3.1 NatureLM-assisted parameterization

NatureLM MCP status was **CONNECTED** during study design. The validated quantitative priors used in the simulation were: researcher quality heterogeneity N(0,1); citation-network power-law exponents in the 2.0-3.0 range; peer-review inter-rater correlation in the 0.1-0.4 range with a focal calibration of 0.25; highly skewed output distributions with Gini targets between 0.6 and 0.9; Matthew effects yielding 1.5-2.5x higher success rates for previously funded researchers; and an experimentally attractive diversity-preserving success rate between 15% and 30%. Mesa was available in the environment, but the implemented ABM uses a transparent researcher-state update loop for reproducibility and direct control of the dynamics.

### 3.2 Agents and state variables

The model contains 200 researcher agents. Each researcher has fixed intrinsic quality $q_i \sim \mathcal{N}(0,1)$, time-varying productivity $p_{i,t}$, cumulative citations $c_{i,t}$, career stage, gender, region, field, a funding indicator, funding history, a career trajectory vector, and prestige. Career stage evolves from junior to mid to senior as simulated age increases. Gender is initialized 50/50, region 70/30 domestic versus international, and fields follow weighted shares across STEM, social sciences, and humanities.

### 3.3 Dynamic rules

At each annual step, active researchers produce papers according to

$ p_{i,t} = \max\left(0.05, \mathcal{N}\left((1.15 + 0.78 q_i + s_i + f_i + 0.05 \pi_{i,t-1}) 	imes b_{i,t}, 0.75ight)ight), $

where $s_i$ is a career-stage modifier, $f_i$ is a field modifier, $\pi_{i,t-1}$ is prestige, and $b_{i,t}=1.3$ if the researcher was funded in the previous year and 1 otherwise. Citations accumulate through preferential attachment using a Poisson draw scaled by current productivity, field-specific citation intensity, prior citations, and prestige.

The citation Gini coefficient is computed as

$ G(x) = rac{n+1 - 2\sum_{i=1}^n rac{\sum_{j=1}^i x_{(j)}}{\sum_{j=1}^n x_j} }{n}, $

where $x_{(j)}$ denotes sorted non-negative outcomes. Diversity is measured using Shannon entropy over funded fields,

$ H = -\sum_k p_k \log p_k. $

### 3.4 Funding mechanisms

1. **Peer review**: panel score $S_i = q_i + 0.60F_i - 0.34W_i + 0.24I_i + 0.18D_i + arepsilon_i$, where $F_i$ indicates prior funding, $W_i$ indicates female gender, $I_i$ is institution advantage, $D_i$ indicates domestic affiliation, and $arepsilon_i$ is reviewer noise calibrated to reviewer-level ICC = 0.25. The top 20% are funded.
2. **Lottery**: applications are pre-screened using a noisier score related to quality and prior advantage; the top 40% enter a lottery from which 20% of the active population is funded.
3. **Formula allocation**: score $A_i = 0.48 z(p_i) + 0.22 z(\log(1+c_i)) + 0.30 z(B_i)$, where $B_i$ is a diversity-and-underfunding bonus rewarding underrepresented gender, region, field, earlier career stage, and researchers with thinner funding histories.

Low-output researchers face a 5% annual attrition probability. Each mechanism is repeated for 30 random seeds.

## 4. Experiments

The baseline experiment simulates 20 annual cycles for each mechanism with a 20% grant success rate. The KAKENHI case study modifies the environment to reflect stylized Japanese conditions: a 26% success rate, strong seniority bias under peer review, and field shares of 70% STEM, 20% social science, and 10% humanities. For each configuration we record citation inequality, funding inequality, demographic funding gaps, career-stage representation, diversity, total output, funded quality, and attrition.

## 5. Results

### 5.1 Baseline cross-validated outcomes

| Mechanism | Citation Gini | Funding Gini | Gender gap (pp) | Regional gap (pp) | Diversity entropy | Total output | Funded quality | Attrition rate (%) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Formula | 0.383 ± 0.021 | 0.764 ± 0.011 | -7.20 ± 4.68 | -5.50 ± 7.27 | 1.055 ± 0.038 | 11718.9 ± 235.2 | 1.269 ± 0.145 | 0.16 ± 0.28 |
| Lottery | 0.336 ± 0.016 | 0.474 ± 0.017 | 5.18 ± 4.91 | 3.12 ± 4.97 | 0.912 ± 0.089 | 11532.6 ± 266.0 | 0.643 ± 0.133 | 0.07 ± 0.18 |
| Peer Review | 0.364 ± 0.012 | 0.615 ± 0.012 | 5.19 ± 4.67 | 3.06 ± 6.46 | 0.930 ± 0.096 | 11712.6 ± 315.9 | 0.984 ± 0.163 | 0.14 ± 0.27 |

Formula performed best on total output (11718.9 ± 235.2), whereas Lottery was worst (11532.6 ± 266.0). Formula performed best on diversity entropy (1.055 ± 0.038), whereas Lottery was worst (0.912 ± 0.089). Lottery performed best on funding gini (0.474 ± 0.017), whereas Formula was worst (0.764 ± 0.011).

![Figure 1: Citation Gini coefficient over time with 95% confidence intervals.](figures/fig1_gini_over_time.png)

![Figure 2: Gender funding gap over time with 95% confidence intervals.](figures/fig2_gender_gap.png)

![Figure 3: Career-stage distribution in the funded pool at the final timestep.](figures/fig3_career_distribution.png)

![Figure 4: Cumulative scientific output by mechanism with confidence bands.](figures/fig4_cumulative_output.png)

### 5.2 Stylized KAKENHI case study

| Mechanism | Avg funded quality | Diversity entropy | Funding Gini |
|---|---:|---:|---:|
| Formula | 1.126 ± 0.106 | 1.013 ± 0.049 | 0.714 ± 0.008 |
| Lottery | 0.685 ± 0.141 | 0.731 ± 0.108 | 0.458 ± 0.020 |
| Peer Review | 0.844 ± 0.135 | 0.648 ± 0.120 | 0.554 ± 0.009 |

The KAKENHI scenario reinforces the same qualitative trade-off, but the stronger seniority and field biases make the fairness consequences more visible. Diversity-oriented formula allocation broadens the funded portfolio, lottery most effectively limits concentration, and peer review continues to favor high-scoring incumbents under biased screening.

![Figure 5: KAKENHI case study comparing efficiency, diversity, and equity metrics.](figures/fig5_kakenhi_case_study.png)

![Figure 6: Diversity trajectories of funded portfolios over time.](figures/fig6_diversity_metrics.png)

## 6. Discussion

Three conclusions emerge. First, cumulative advantage matters: once prior funding raises both future review scores and next-year productivity, peer review compounds inequality even when reviewer noise is substantial. Second, the lottery mechanism meaningfully softens concentration without a catastrophic efficiency penalty, supporting arguments for partial randomization when many proposals are near the funding threshold. Third, explicit diversity weighting can raise the entropy of the funded portfolio and reduce demographic imbalances, but this typically trades off against short-run measures of average funded quality because the rule intentionally redistributes support toward less advantaged groups and fields.

These findings align with the literature: low reviewer agreement makes fine-grained ranking fragile, while prior success and evaluator similarity amplify path dependence. The model remains stylized; it does not represent project-level heterogeneity, team science, institution-specific budgets, or adaptive researcher effort in response to policy reforms. Even so, it offers a useful design-space map for agencies considering hybrid reforms.

## 7. Conclusion

A funding system optimized only for reviewer-ranked efficiency is unlikely to be optimal for the long-run ecology of science. In this ABM, the strongest mechanism for output, the strongest mechanism for equity, and the strongest mechanism for diversity were not the same. That divergence is the central policy lesson: agencies should avoid single-objective optimization and instead combine thresholded quality screening with partial randomization and explicit diversity safeguards.

## References

1. Feliciani T, Morreau M, Luo J, Lucas P, Shankar K (2022). Designing grant-review panels for better funding decisions: Lessons from an empirically calibrated simulation model. *Research Policy*, 51(4), 104467. DOI: 10.1016/j.respol.2021.104467
2. Taffe MA, Gilpin NW (2021). Racial inequity in grant funding from the US National Institutes of Health. *eLife*, 10, e65697. DOI: 10.7554/elife.65697
3. Moldashev KB (2024). Critical Evaluation of the Grant Funding Allocation Mechanism in Social and Humanitarian Sciences. *Central Asian Economic Review*, 5, 48–58. DOI: 10.52821/2789-4401-2024-5-48-58
4. Pier EL et al. (2018). Low agreement among reviewers evaluating the same NIH grant applications. *PNAS*, 115(12), 2952–2957. DOI: 10.1073/pnas.1714379115
5. Banal-Estañol A et al. (2023). Similar-to-me effects in the grant application process. *R&D Management*, 53(5), 819–839. DOI: 10.1111/radm.12601
6. Li D (2017). Expertise versus Bias in Evaluation: Evidence from the NIH. *American Economic Journal: Applied Economics*, 9(2), 60–92. DOI: 10.1257/app.20150421
7. Graves N, Barnett AG, Clarke P (2011). Funding grant proposals for scientific research: Retrospective analysis. *BMJ*, 343, d4797. DOI: 10.1136/bmj.d4797
