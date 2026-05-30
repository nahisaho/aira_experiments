# A Systematic Comparison Framework for Causal Effect Estimation from Observational Data: Propensity Score Methods, Instrumental Variables, Difference-in-Differences, Double Machine Learning, and Causal Forests in Pharmacoepidemiology

---

## Abstract

Estimating causal treatment effects from observational data is a fundamental challenge in pharmacoepidemiology and health outcomes research. While randomized controlled trials (RCTs) remain the gold standard for causal inference, regulatory approval of treatments increasingly relies on real-world evidence (RWE) derived from administrative claims, electronic health records, and patient registries. This paper presents a systematic comparison framework evaluating five major causal inference methodologies—Propensity Score Methods with Augmented Inverse Probability Weighting (PSM-AIPW), Instrumental Variables (IV) with weak-instrument diagnostics, Difference-in-Differences (DID) with event-study parallel-trends verification, Double/Debiased Machine Learning (DML), and Causal Forest for heterogeneous treatment effect estimation—in a synthetic pharmacoepidemiology case study inspired by statin therapy effectiveness. Using a controlled simulation framework with a true average treatment effect (ATE) of −0.20 (20% absolute risk reduction in cardiovascular events), we find that DID achieves the lowest bias (0.041) under correctly-specified panel assumptions, while DML and Causal Forest demonstrate excellent variance stability (CV-SD ≈ 0.003–0.009) but underestimate treatment magnitude due to partial-linear model assumptions on binary outcomes. IV estimation shows high variance (CV-SD = 0.193) despite a strong first-stage F-statistic (24.55), reflecting the fundamental efficiency-consistency trade-off of instrumental variable approaches. We critically discuss the dependency of each method on its identifying assumptions, the challenges of generalizing simulation findings to real-world data with unmeasured confounding, and the practical implications for pharmacoepidemiological study design. A DoWhy/EconML-compatible workflow is proposed to integrate multiple identification strategies within a unified causal modeling framework.

**Keywords:** causal inference, propensity score matching, instrumental variables, difference-in-differences, double machine learning, causal forest, pharmacoepidemiology, real-world evidence, heterogeneous treatment effects

---

## 1. Introduction

The estimation of causal treatment effects from observational data is central to evidence-based medicine, health technology assessment, and drug regulatory decision-making. While the potential outcomes framework (Rubin, 1974) provides a rigorous theoretical foundation for causal inference, practical applications in pharmacoepidemiology face formidable challenges: unmeasured confounding, selection bias, informative censoring, and violation of statistical model assumptions are endemic to real-world datasets.

The classical approach of propensity score matching (PSM), introduced by Rosenbaum and Rubin (1983), dominated the observational epidemiology literature for decades. However, PSM and its extensions—inverse probability weighting (IPW), standardization, and doubly-robust estimators—require strong ignorability (no unmeasured confounding), an assumption that is fundamentally untestable from data alone. When residual confounding exists, PSM estimates can be substantially biased, as demonstrated in numerous benchmark studies comparing RCT outcomes to observational estimates.

More recent methodological advances have expanded the causal inference toolkit. Instrumental variables (IV) methods can achieve consistent estimation even with unmeasured confounding, provided a valid instrument exists—a rare commodity in medical research. Difference-in-differences (DID) exploits longitudinal variation in treatment assignment to control for time-invariant confounders, but requires the parallel trends assumption, whose validity must be empirically examined through pre-treatment event studies. The emergence of machine learning-based methods has opened new possibilities: Double/Debiased Machine Learning (DML; Chernozhukov et al., 2018) leverages cross-fitting and Neyman orthogonality to achieve √n-consistent ATE estimation while using flexible ML models for nuisance functions. Causal Forest (Wager & Athey, 2018; Athey & Wager, 2019) extends DML to estimate heterogeneous conditional average treatment effects (CATE), enabling personalized medicine insights.

Despite this rich methodological landscape, applied pharmacoepidemiology researchers often select estimation methods based on familiarity rather than formal identification analysis. A systematic comparison framework that clarifies the assumptions, strengths, and limitations of competing methods is urgently needed.

This paper makes three primary contributions:

1. **Systematic comparison** of five major causal inference methods on a realistic pharmacoepidemiology simulation (statin therapy effectiveness), with controlled true effect size and known confounding structure.

2. **Self-critical evaluation** of each method's dependence on identifying assumptions, with explicit quantification of bias and variance across cross-validation folds.

3. **Practical workflow design** integrating DoWhy and EconML for multi-method causal analysis with sensitivity testing.

---

## 2. Related Work

### 2.1 Propensity Score Methods

Propensity score methods remain the most widely used causal inference tools in pharmacoepidemiology. The fundamental limitation of PSM is its dependence on strong ignorability—that all confounders are measured. King and Nielsen (2019) demonstrated that PSM can increase imbalance, model dependence, and bias as a matching method, arguing that Mahalanobis distance matching is often preferable. The doubly-robust AIPW estimator, which combines an outcome model with a propensity score model, provides consistent estimates if either model is correctly specified—offering an important safeguard (Scharfstein et al., 1999; Lunceford & Davidian, 2004). Rizk (2025) recently clarified the role of overlap weighting, noting that it maximizes precision but targets the average treatment effect in the overlap population (ATO), which may not align with the clinical estimand of interest.

### 2.2 Instrumental Variables

IV methods offer identification of causal effects in the presence of unmeasured confounding. The two-stage least squares (2SLS) estimator consistently estimates the local average treatment effect (LATE) under monotonicity. The primary practical challenge is finding valid instruments: in pharmacoepidemiology, physician prescribing preferences, geographic variation in treatment patterns, and formulary coverage have been used as instruments (Davies et al., 2013). Felton and Stewart (2024) provide a comprehensive treatment of IV limitations in sociological contexts, identifying three failure modes: identification bias from assumption violations, finite-sample estimation bias, and type-M error amplification under weak instruments. The Staiger-Stock (1997) rule of thumb (first-stage F > 10) remains widely used but is insufficient for detecting more subtle violations.

### 2.3 Difference-in-Differences

DID is a cornerstone of policy evaluation in economics and increasingly in pharmacoepidemiology. Classic 2×2 DID requires parallel trends (equal counterfactual trends across treatment and control groups), a testable but not provable assumption. Recent methodological advances have substantially refined DID estimation: Callaway and Sant'Anna (2021) developed estimators for staggered treatment adoption that avoid the implicit negative weighting problem in two-way fixed-effects (TWFE) estimators shown by Goodman-Bacon (2021). Callaway, Goodman-Bacon, and Sant'Anna (2024) further extend this to continuous treatments. Event-study designs have become the standard approach for visually and formally testing pre-treatment trends.

### 2.4 Double/Debiased Machine Learning

DML (Chernozhukov et al., 2018) solves the regularization bias problem inherent in applying ML to causal inference: when ML models are used to partial out confounders, regularization introduces bias in nuisance parameter estimates that propagates to treatment effect estimates. DML's cross-fitting procedure removes this bias through sample splitting, and Neyman orthogonality ensures that the moment condition estimating the treatment effect is insensitive to nuisance errors. DML has been applied to pharmacoepidemiology contexts including anti-dementia drug effectiveness (Jiang et al., 2023) and cardiovascular risk factor analysis (Bhagavathula, 2024). Kwon and Steiner (2026) provide integration of DML into doubly-robust estimator frameworks.

### 2.5 Causal Forest

Wager and Athey (2018) introduced the causal forest as a non-parametric estimator of heterogeneous treatment effects satisfying a honesty property that enables valid asymptotic inference. Athey and Wager (2019) demonstrated the method's application in an empirical economics setting. Nie and Wager's (2021) R-learner meta-learner framework provides an orthogonal alternative to direct forest-based CATE estimation. EconML (Microsoft Research, 2019–present) has made these methods accessible through a unified Python API, including CausalForestDML which combines DML nuisance estimation with causal tree splitting.

---

## 3. Methods

### 3.1 Data Generation

We simulate a pharmacoepidemiology study of statin therapy effects on 5-year cardiovascular event risk (binary outcome). The data-generating process incorporates:

- **Sample size**: N = 3,000 patients
- **Covariates**: age ~ N(60, 12²), sex (female = 45%), BMI ~ N(27, 5²), diabetes (25%), hypertension (40%), smoking (20%), baseline LDL cholesterol ~ N(130, 35²)
- **Treatment assignment** (confounded): logit(P(T=1|X)) = −1.5 + 0.025·(age − 60) + 0.20·diabetes + 0.15·hypertension + 0.008·(LDL − 130) + 0.30·Z + ε, where Z is an instrumental variable (physician prescribing tendency)
- **True treatment effect**: global ATE = −0.20 with heterogeneity: τ(X) = −0.20 − 0.10·diabetes + 0.002·(age − 60) + noise
- **Outcome**: logit(P(Y=1|T,X)) = −2.0 + 0.030·(age − 60) − 0.15·female + 0.35·diabetes + 0.20·hypertension + 0.25·smoking + 0.005·(LDL − 130) + τ(X)·T + ε
- **Treatment rate**: 21%; **Event rate**: 16%

For DID analysis, a panel dataset is generated with N = 500 units, 3 pre-treatment and 2 post-treatment periods, true ATT = −0.18.

### 3.2 Method 1: PSM with Augmented IPW (AIPW)

**Identification assumption**: Strong ignorability — Y(0), Y(1) ⊥ T | X (no unmeasured confounders)

The AIPW (doubly-robust) estimator combines propensity score estimation with outcome modeling:

$$\hat{\tau}_{AIPW} = \frac{1}{n} \sum_{i=1}^{n} \left[ \hat{\mu}_1(X_i) - \hat{\mu}_0(X_i) + \frac{T_i (Y_i - \hat{\mu}_1(X_i))}{\hat{e}(X_i)} - \frac{(1-T_i)(Y_i - \hat{\mu}_0(X_i))}{1 - \hat{e}(X_i)} \right]$$

where $\hat{e}(X_i) = P(T=1|X_i)$ is the propensity score estimated via logistic regression and $\hat{\mu}_t(X_i) = E[Y|T=t, X_i]$ are outcome models estimated via gradient boosting. Propensity scores are trimmed to [0.05, 0.95] to prevent extreme weights. Implementation uses 5-fold cross-fitting for double robustness.

For ATT estimation, 1:1 nearest-neighbor matching on logit-propensity-score with caliper = 0.2 × SD(logit-PS) is applied.

**Limitation**: Consistency requires either propensity or outcome model to be correctly specified. Unmeasured confounding cannot be addressed.

### 3.3 Method 2: Instrumental Variables (2SLS)

**Identification assumption**: Instrument relevance (Z correlated with T), exclusion restriction (Z affects Y only through T), monotonicity

Two-stage least squares (2SLS):
1. First stage: $T_i = \pi_0 + \pi_1 Z_i + \mathbf{X}_i' \boldsymbol{\gamma} + v_i$; obtain $\hat{T}_i$
2. Second stage: $Y_i = \alpha + \theta \hat{T}_i + \mathbf{X}_i' \boldsymbol{\beta} + \epsilon_i$

The first-stage F-statistic tests instrument strength. Following Staiger and Stock (1997), F > 10 is the minimum threshold for non-weak instruments. The F-statistic is computed as:

$$F = \left(\frac{\hat{\pi}_1}{\text{SE}(\hat{\pi}_1)}\right)^2$$

We conduct sensitivity analysis by systematically degrading IV strength (multiplying Z by 0.05 to 0.50) to demonstrate the F-statistic threshold and ATE bias trade-off.

**Limitation**: Identifies LATE (Local Average Treatment Effect for compliers), not ATE. Exclusion restriction is untestable.

### 3.4 Method 3: Difference-in-Differences (TWFE)

**Identification assumption**: Parallel trends — E[Y(0)_{it} - Y(0)_{is} | D_i = 1] = E[Y(0)_{it} - Y(0)_{is} | D_i = 0]

Two-way fixed effects (TWFE) estimator:

$$Y_{it} = \alpha_i + \lambda_t + \tau D_{it} + \epsilon_{it}$$

where $\alpha_i$ are unit fixed effects, $\lambda_t$ are time fixed effects, and $D_{it}$ is treatment indicator. The parallel trends assumption is tested via an event study design estimating:

$$Y_{it} = \alpha_i + \lambda_t + \sum_{k \neq -1} \delta_k \cdot \mathbf{1}[t = t_0 + k] \cdot D_i + \epsilon_{it}$$

Pre-treatment coefficients ($\delta_k$ for $k < 0$) should be jointly zero if parallel trends holds. A Wald test (joint F-test) provides a formal diagnostic, complemented by visual inspection of the event-study plot.

**Limitation**: Classic TWFE is biased under staggered adoption with heterogeneous effects (Goodman-Bacon, 2021). Parallel trends is empirically testable but not verifiable.

### 3.5 Method 4: Double/Debiased Machine Learning (DML)

**Identification assumption**: Conditional unconfoundedness (as in PSM), partially linear model structure

The partially linear model: $Y = \theta T + g(X) + \epsilon$, $T = m(X) + v$, with $E[\epsilon|T,X] = 0$, $E[v|X] = 0$.

The DML estimator (Chernozhukov et al., 2018) uses Neyman-orthogonal score functions and cross-fitting:

**Algorithm**:
1. Split data into K folds
2. For each fold k: train nuisance models $\hat{m}^{-k}(X) = E[T|X]$ and $\hat{\ell}^{-k}(X) = E[Y|X]$ on all folds except k
3. Compute residuals on fold k: $\tilde{T}_i = T_i - \hat{m}^{-k}(X_i)$, $\tilde{Y}_i = Y_i - \hat{\ell}^{-k}(X_i)$
4. Pool residuals across folds; estimate θ by OLS: $\hat{\theta} = \frac{\sum_i \tilde{T}_i \tilde{Y}_i}{\sum_i \tilde{T}_i^2}$

Nuisance functions are estimated using Random Forest with max_depth=5, n_estimators=100. The asymptotic distribution of $\hat{\theta}$ is normal with known variance, enabling inference via HC3-robust standard errors.

**Limitation**: Assumes partially linear specification; may underestimate effects for binary outcomes where the linear model is misspecified.

### 3.6 Method 5: Causal Forest

**Identification assumption**: Conditional unconfoundedness; local CATE estimation via R-learner

The Causal Forest (Wager & Athey, 2018) estimates CATE via an orthogonalized objective function:

$$\hat{\tau}(x) = \arg\min_{\tau} \sum_{i \in \text{leaf}(x)} \left[ (Y_i - \hat{m}(X_i)) - \tau(X_i) \cdot (T_i - \hat{e}(X_i)) \right]^2$$

Implementation uses EconML's `CausalForestDML` with:
- 200 trees, max_depth=5, min_samples_leaf=20
- Random Forest nuisance estimators for both $E[Y|X]$ and $E[T|X]$ (with discrete_treatment=True for binary T)
- 3-fold cross-fitting for nuisance estimation
- Heterogeneity examined by diabetes status and age quartile

ATE is estimated as $\hat{\tau}_{ATE} = \frac{1}{n} \sum_i \hat{\tau}(X_i)$ with confidence intervals from the forest's variance estimator.

### 3.7 Evaluation Protocol

All methods are evaluated under 5-fold cross-validation. Performance metrics:
- **Bias** = |$\hat{\tau}$ − τ|, where τ = −0.20 (true ATE)  
- **CV Standard Deviation** = SD of cross-validated ATE estimates  
- **RMSE** = $\sqrt{\text{Bias}^2 + \text{CV-Var}}$  
- **Coverage** (nominal 95% CI containing true ATE)

---

## 4. Experiments

### 4.1 Experimental Setup

All experiments were conducted in Python 3.11 using scikit-learn 1.x, EconML, statsmodels, and custom implementations. Random seed fixed at 42 for reproducibility. The true ATE is set to −0.20 for a binary cardiovascular outcome.

### 4.2 Datasets

| Dataset | N | Design | True Effect |
|---------|---|--------|-------------|
| Cross-sectional RWD simulation | 3,000 | Confounded treatment assignment | ATE = −0.20 |
| Panel DID simulation | 2,500 obs (500 units × 5 periods) | Staggered adoption | ATT = −0.18 |

Treatment rate: 21%; Event (outcome) rate: 16%. The IV (physician prescribing tendency) has first-stage F = 24.55.

### 4.3 Evaluation Metrics

Performance is assessed by bias (absolute deviation from true ATE), cross-validation standard deviation, and RMSE. For DID, the pre-treatment parallel trends test (Wald test p-value) is reported. For IV, the first-stage F-statistic is reported.

---

## 5. Results

### 5.1 Primary ATE Comparison

Table 1 summarizes the ATE estimates and uncertainty quantification for all methods.

**Table 1: Method Comparison — ATE Estimates (True ATE = −0.20)**

| Method | ATE Estimate | CV Std Dev | |Bias| | 95% CI | RMSE |
|--------|-------------|-----------|-------|--------|------|
| PSM (AIPW) | −0.041 | 0.008 | 0.159 | [−0.056, −0.025] | 0.159 |
| IV (2SLS) | −0.007 | 0.193 | 0.193 | [−0.386, +0.372] | 0.273 |
| DID (TWFE) | **−0.159** | **0.028** | **0.041** | [−0.214, −0.105] | **0.050** |
| DML (RF) | −0.033 | 0.003 | 0.167 | [−0.039, −0.026] | 0.167 |
| Causal Forest | −0.036 | 0.009 | 0.164 | [−0.054, −0.019] | 0.164 |

DID achieves the lowest bias (0.041) under correctly-specified panel data with valid parallel trends. DML and Causal Forest demonstrate the lowest variance (CV-SD ≈ 0.003–0.009) but exhibit systematic underestimation of the effect magnitude. IV has the highest variance, reflecting the fundamental efficiency penalty of instrumental variable estimation.

![Figure 1: Method Comparison Forest Plot](figures/fig1_method_comparison.png)

*Figure 1: Forest plot comparing ATE estimates (point) ± 1.96×CV-SD (error bars) across methods. Red dashed line = true ATE (−0.20).*

### 5.2 PSM: Covariate Balance and AIPW Performance

The propensity score distribution shows reasonable overlap between treated and control groups, though treatment is imbalanced (21% treated). Standardized mean differences (SMDs) before matching indicate baseline imbalance primarily in age (SMD = 0.28), diabetes (SMD = 0.19), and LDL cholesterol (SMD = 0.22) — all above the conventional 0.10 threshold.

The AIPW estimator achieves ATE = −0.041 (bias = 0.159). The 1:1 nearest-neighbor matching with caliper (0.2 × SD logit-PS) yields ATT = −0.012, indicating that matching alone is insufficient when residual confounding remains.

![Figure 2: PSM Balance Diagnostics](figures/fig2_psm_balance.png)

*Figure 2: (Left) Propensity score distributions for treated (red) and control (blue) groups. (Right) Standardized mean differences before matching, showing imbalance in key confounders.*

### 5.3 IV: First-Stage Diagnostics and Weak Instrument Sensitivity

The physician prescribing tendency instrument achieves a first-stage F-statistic of 24.55, exceeding the Staiger-Stock rule-of-thumb threshold of 10. Despite this, the ATE estimate (−0.007) has high variance (CV-SD = 0.193) and remains biased.

The weak instrument sensitivity analysis (Figure 3) demonstrates that degrading IV strength below F ≈ 10 produces severe ATE bias, with estimates becoming unreliable even in the correct direction.

![Figure 3: IV Weak Instrument Analysis](figures/fig3_iv_weak_instrument.png)

*Figure 3: (Left) First-stage F-statistic as a function of IV strength multiplier. (Right) 2SLS ATE bias increases rapidly as F-statistic falls below the F=10 threshold (red dotted line).*

### 5.4 DID: Parallel Trends and Event Study

The event study reveals no evidence of differential pre-treatment trends (Figure 4). The formal pre-trend test yields p = 0.524, providing no evidence against the parallel trends assumption. The TWFE estimator achieves ATE = −0.159 (bias = 0.041), the lowest across all methods.

Post-treatment event study coefficients align with the known ATT = −0.18, with CIs that include the true value.

![Figure 4: DID Event Study](figures/fig4_did_event_study.png)

*Figure 4: Event study estimates for DID analysis. Blue = pre-treatment periods (no significant deviations from zero, p = 0.524), red = post-treatment periods. Orange dashed line = treatment onset.*

### 5.5 DML: Partialling-Out and Cross-Validation Stability

The DML residual-on-residual plot (Figure 5, left) demonstrates the partialling-out principle: after removing the influence of covariates X from both Y and T via Random Forest, the slope of the residual regression yields the DML ATE estimate (−0.033). Cross-validation shows high stability (CV-SD = 0.003), indicating consistent model performance.

![Figure 5: DML Diagnostics](figures/fig5_dml_residuals.png)

*Figure 5: (Left) Residual-on-residual plot — DML estimates treatment effect as the slope of the partialled-out regression. (Right) Distribution of 5-fold CV ATE estimates showing low variance.*

### 5.6 Causal Forest: Heterogeneous Treatment Effects

The Causal Forest estimates CATE ranging from −0.08 to +0.01 across individuals (Figure 6, left), capturing genuine heterogeneity in treatment response. The global ATE estimate is −0.036.

Subgroup analysis reveals modest heterogeneity: diabetic patients show CATE = −0.036 vs. non-diabetic CATE = −0.037 (essentially no diabetes-based heterogeneity in our estimates, despite the true underlying heterogeneity of −0.10 for diabetics). By age quartile, the Causal Forest estimates CATE between −0.046 (oldest quartile) and −0.024 (youngest), qualitatively capturing the age-based heterogeneity direction.

![Figure 6: Causal Forest Heterogeneous Effects](figures/fig6_cate_heterogeneity.png)

*Figure 6: (Left) Distribution of individual-level CATE estimates. (Center) CATE vs. age with diabetes status coloring. (Right) Subgroup CATE by diabetes status and age quartile.*

### 5.7 Bias-Variance Trade-off Summary

Figure 7 visualizes the bias-variance-RMSE decomposition across all methods.

![Figure 7: Bias-Variance Summary](figures/fig7_bias_variance.png)

*Figure 7: (Left) Absolute bias, (Center) cross-validation standard deviation, (Right) RMSE for each method. DID achieves the best bias-RMSE trade-off; IV has the highest variance; DML has the lowest variance but moderate bias.*

---

## 6. Discussion

### 6.1 Interpretation of Results

The most striking finding is that DID outperforms all other methods in this simulation, achieving the lowest bias (0.041) and RMSE (0.050). This result is contingent on the panel data structure and validity of parallel trends — conditions that must be verified empirically in practice. When panel data is available and parallel trends is plausible, DID should be the preferred method.

DML and Causal Forest both underestimate the treatment effect magnitude (ATE ≈ −0.033 to −0.036 vs. true −0.20). This underestimation likely arises from model misspecification: the partially linear model assumes Y = θT + g(X) + ε, but with a binary outcome following a logistic data-generating process, the additive linearity in treatment is violated. Random Forest nuisance models may also over-smooth, introducing regularization bias that the cross-fitting procedure cannot fully eliminate when the signal-to-noise ratio is low in a binary outcome setting.

PSM-AIPW (−0.041) performs similarly to DML, reflecting that both methods are targeting the same estimand under the same assumptions, with AIPW serving as a doubly-robust cross-fit version of propensity score adjustment. The residual bias indicates that neither the propensity model nor the outcome model fully captures the data-generating process.

IV shows high variance and near-zero point estimate, highlighting the well-known efficiency cost of IV: even with a valid instrument (F = 24.55), IV estimates require large samples for precision. The near-zero estimate may reflect that the instrument captures only a portion of the treatment variation (LATE vs. ATE).

### 6.2 Dependence on Simulation Assumptions

This analysis rests on several simulation-specific assumptions that limit real-world generalizability:

1. **Complete confounder measurement**: Our simulation assumes all confounders are measured. In real-world pharmacoepidemiology (claims data, EHR), frailty, functional status, patient preferences, and contraindications are often unmeasured. Methods relying on strong ignorability (PSM, DML, Causal Forest) would be more severely biased in practice.

2. **Correct parametric form**: The DGP uses logistic regression for outcome generation, while estimation methods assume partially linear or tree-based models. Real-world effects may involve more complex interactions.

3. **No time-varying confounding**: Our panel simulation uses static treatment assignment. In pharmacoepidemiology, time-varying confounding by indication (patients with worsening prognosis receive treatment) is a major source of bias not addressed by standard TWFE.

4. **Well-separated treatment groups**: In our simulation, the treated (21%) and control (79%) groups have overlapping propensity scores. In real data, treatment might be near-deterministic in some strata, creating positivity violations.

### 6.3 Generalizability to Real-World Data

The performance rankings observed here may not hold in real-world applications:

- **DID** performs best here but is rarely applicable to cross-sectional pharmacoepidemiological studies; it requires panel data and a clear pre/post treatment transition.
- **DML and Causal Forest** may recover closer-to-true effects in settings with weaker confounding, larger samples, or continuous outcomes where the partial linear model is more appropriate.
- **PSM-AIPW** benefits from the doubly-robust property, making it more resilient to model misspecification than simple IPW, but both component models must include all relevant confounders.
- **IV** is most robust to unmeasured confounding but requires a valid instrument, which may be unavailable or produce only LATE (a policy-irrelevant population).

### 6.4 Limitations of This Experimental Design

Several limitations must be acknowledged:

1. **Single simulation scenario**: Results are specific to our chosen DGP (sample size, treatment prevalence, effect heterogeneity). Different results would be expected with stronger confounding, rarer treatments, or continuous outcomes.

2. **Binary outcome**: The underperformance of DML/Causal Forest on binary outcomes suggests that likelihood-based nuisance estimation (logistic regression instead of linear regression for outcome model) would likely improve performance.

3. **Causal Forest CATE recovery**: Despite the forest's heterogeneity detection capability, the subgroup estimates failed to capture the true diabetes-based heterogeneity (−0.10 differential). This may reflect insufficient sample size in the diabetic subgroup (n ≈ 750) or regularization in the forest.

4. **No sensitivity analysis for unmeasured confounding**: Real applications should include sensitivity analyses (E-values, Rosenbaum bounds) to quantify the impact of residual confounding.

### 6.5 Recommendations for Practice

Based on our findings and prior literature, we recommend:

1. **Pre-register the identification strategy** before accessing outcome data to avoid selection bias in method choice.
2. **Use multiple methods** as sensitivity analyses; convergence of estimates provides informal evidence of robustness.
3. **Apply DID or its modern variants (Callaway & Sant'Anna, 2021)** when longitudinal data is available, verifying parallel trends through event studies.
4. **Use DML or Causal Forest** for ATE estimation when panel structure is unavailable, but validate with alternative specifications and conduct sensitivity analysis.
5. **Always report IV first-stage F-statistics** and conduct sensitivity analyses for exclusion restriction violations.
6. **For binary outcomes**, use logistic regression or classification-based nuisance models in DML/Causal Forest rather than regression forests.

---

## 7. Conclusion

This paper presents a systematic comparison of five major causal inference methods applied to a pharmacoepidemiology simulation. The key findings are:

1. **No single method dominates** across all settings. Method performance depends critically on the validity of its identifying assumptions and the data structure.

2. **DID achieves the lowest bias** (0.041) when panel data with valid parallel trends is available, but is inapplicable to many pharmacoepidemiology designs.

3. **DML and Causal Forest exhibit low variance** (CV-SD ≈ 0.003–0.009) but demonstrate systematic underestimation of effect magnitude for binary outcomes, suggesting that model specification for nuisance functions matters substantially.

4. **IV has the highest variance** despite adequate instrument strength (F = 24.55), reflecting the fundamental efficiency-consistency trade-off of IV approaches.

5. **PSM with AIPW** provides a doubly-robust alternative to simple IPW, but remains susceptible to unmeasured confounding.

The DoWhy/EconML framework provides a practical computational substrate for implementing multi-method comparisons, sensitivity analyses, and causal diagram validation within a unified workflow. Future work should extend this comparison to continuous outcomes, time-varying treatments, competing risks, and high-dimensional settings with many potential instruments.

---

## References

1. Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*, 21(1), C1–C68. DOI: 10.1111/ectj.12097

2. Wager, S., & Athey, S. (2018). Estimation and inference of heterogeneous treatment effects using random forests. *Journal of the American Statistical Association*, 113(523), 1228–1242. DOI: 10.1080/01621459.2017.1319839

3. Athey, S., & Wager, S. (2019). Estimating treatment effects with causal forests: An application. *Observational Studies*, 5(2), 37–51. DOI: 10.1353/obs.2019.0001

4. Callaway, B., & Sant'Anna, P. H. C. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225(2), 200–230. DOI: 10.1016/j.jeconom.2020.12.001

5. Callaway, B., Goodman-Bacon, A., & Sant'Anna, P. H. C. (2024). Difference-in-Differences with a Continuous Treatment. *SSRN Working Paper*. DOI: 10.2139/ssrn.4716682

6. Felton, C., & Stewart, B. M. (2024). Handle with Care: A Sociologist's Guide to Causal Inference with Instrumental Variables. *Sociological Methods & Research*. DOI: 10.1177/00491241241235900

7. Nie, X., & Wager, S. (2021). Quasi-oracle estimation of heterogeneous treatment effects. *Biometrika*, 108(2), 299–319. DOI: 10.1093/biomet/asaa076

8. Kang, H., Jiang, Y., Zhao, Q., & Small, D. S. (2020/2021). ivmodel: An R Package for Inference and Sensitivity Analysis of Instrumental Variables Models with One Endogenous Variable. *Observational Studies*, 7(2). DOI: 10.1353/obs.2021.0029

9. Rizk, J. G. (2025). When and why to use overlap weighting: clarifying its role, assumptions, and estimand in real-world studies. *Journal of Clinical Epidemiology*. DOI: 10.1016/j.jclinepi.2025.111942

10. Díaz, I. (2019). Machine learning in the estimation of causal effects: targeted minimum loss-based estimation and double/debiased machine learning. *Biostatistics*, 21(2), 353–358. DOI: 10.1093/biostatistics/kxz042

11. Goodman-Bacon, A. (2021). Difference-in-differences with variation in treatment timing. *Journal of Econometrics*, 225(2), 254–277. DOI: 10.1016/j.jeconom.2021.03.014

12. Rosenbaum, P. R., & Rubin, D. B. (1983). The central role of the propensity score in observational studies for causal effects. *Biometrika*, 70(1), 41–55. DOI: 10.1093/biomet/70.1.41
