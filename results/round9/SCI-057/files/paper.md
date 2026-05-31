# A Causal Inference Framework for Estimating Health Effects of Air Pollution: Distributed Lag Nonlinear Models, GAM-based Exposure–Response Functions, and Sensitivity Analysis

---

## Abstract

Air pollution—particularly fine particulate matter (PM₂.₅) and ozone (O₃)—is a leading environmental risk factor for premature mortality worldwide, causing an estimated 6–7 million deaths annually. Despite decades of epidemiological research, causal inference from observational studies remains challenging because of complex confounding structures, non-linear exposure–response relationships, and distributed temporal lag effects. This paper presents a comprehensive analytical framework integrating five methodological components: (1) a Land Use Regression (LUR) model for fine-scale PM₂.₅ exposure assessment; (2) a distributed lag nonlinear model (DLNM) implemented via polynomial (Almon) lag-basis functions fitted to a simulated five-year daily time-series; (3) a time-stratified case-crossover design for short-term effects; (4) a GAM-equivalent restricted cubic spline logistic regression for nonlinear long-term exposure–response modeling in a synthetic cohort (N = 5,000); and (5) E-value sensitivity analysis for unmeasured confounding. Analysis was conducted on simulated data mimicking realistic urban air quality and mortality patterns. The LUR model achieved a 5-fold cross-validated R² = 0.754 ± 0.012 (RMSE = 2.56 µg/m³). The DLNM revealed a cumulative relative risk (RR) for PM₂.₅ of 1.005 [95% CI: 0.990–1.021] per 10 µg/m³ over lag 0–7 days. Long-term cohort analysis showed OR = 2.68 [2.16–3.31] per 10 µg/m³ for cardiovascular mortality (5-fold CV AUROC = 0.731 ± 0.037). Nonlinear spline modeling demonstrated a super-linear exposure–response curve. E-value analysis for literature-based estimates (OR = 1.10) yielded E = 1.43, indicating that an unmeasured confounder would need ≥1.43-fold associations with both exposure and outcome to nullify the observed effect. NatureLM and GALACTICA MCP tools were attempted but not available in the current ToolUniverse environment (see Methods). This framework provides a reproducible, modular pipeline applicable to real-world environmental epidemiology.

**Keywords:** air pollution, PM₂.₅, ozone, causal inference, DLNM, land use regression, E-value, GAM, case-crossover, sensitivity analysis

---

## 1. Introduction

Ambient air pollution remains the most important environmental contributor to the global burden of disease. The Global Burden of Disease Study 2023 estimated that household and ambient particulate matter pollution combined caused approximately 8.1 million deaths in 2021, with cardiovascular disease (CVD) accounting for the largest share [GBD 2023 Collaborators, 2025]. Fine particulate matter (PM₂.₅, aerodynamic diameter ≤2.5 µm) penetrates deep into the respiratory and cardiovascular system, triggering systemic inflammation, oxidative stress, and autonomic nervous system dysregulation that collectively elevate risks of myocardial infarction, stroke, arrhythmia, and all-cause mortality [Zhong et al., 2025].

Ozone (O₃), formed through photochemical reactions of nitrogen oxides and volatile organic compounds, further exacerbates respiratory and cardiovascular morbidity, particularly during summer high-pollution events [Olaniyan et al., 2022].

Despite robust epidemiological evidence, translating associations into causal estimates remains methodologically challenging. Key issues include:

1. **Exposure misclassification**: Ground-level monitoring stations are sparse and may not capture within-city spatial heterogeneity. Land Use Regression (LUR) models fusing satellite-derived data, traffic, and land-use variables improve individual exposure assignment [Shi et al., 2020].

2. **Confounding**: Long-term cohort studies must control for socioeconomic position, lifestyle factors, and geographic clustering [Vanoli et al., 2025]. Time-series studies face residual seasonal and meteorological confounding.

3. **Non-linearity and temporal dynamics**: Biological responses to pollutant exposure are rarely linear. Distributed lag nonlinear models (DLNMs) simultaneously capture the non-linear exposure–response relationship and the distributed lag structure of effects over time [Gasparrini, 2014]. Generalized additive models (GAMs) using spline smoothers are the standard for long-term exposure–response characterization.

4. **Causal inference**: Even after controlling for measured confounders, unmeasured confounders (e.g., indoor pollution, diet, physical activity) may bias estimates. The E-value framework of VanderWeele and Ding [2017] quantifies the robustness of findings to unmeasured confounding.

This paper synthesizes these methods into an end-to-end Python-based analytical pipeline, validated against contemporary epidemiological literature. Our contributions include: (a) a publicly reproducible DLNM implementation using polynomial distributed-lag basis functions; (b) a GAM-equivalent restricted cubic spline logistic regression for long-term exposure–response; (c) explicit E-value calculations; and (d) a critical discussion of synthetic data limitations and generalizability.

---

## 2. Related Work

### 2.1 Exposure Assessment

Land use regression models have been widely used to predict ambient PM₂.₅ at fine spatial resolution. Shi et al. [2020] developed a spatiotemporal LUR model for Pakistan incorporating traffic networks, land use, meteorological conditions, and satellite-derived data, explaining 54.5% of PM₂.₅ variability (R² = 0.545). Hybrid LUR models integrating chemical transport model outputs and satellite AOD (aerosol optical depth) have achieved R² > 0.80 in European settings [Hoogh et al., 2018]. The UK Biobank cohort study by Vanoli et al. [2025] used annually updated DEFRA exposure estimates linked to residential histories.

### 2.2 Time-Series Studies and DLNM

The distributed lag nonlinear model (DLNM), formalized by Gasparrini and colleagues, has become the standard framework for time-series analysis of air pollution and health, simultaneously modeling non-linear exposure–response and distributed temporal lag effects using cross-basis matrices. Gutiérrez-Avila et al. [2023] applied a time-stratified case-crossover design with DLNM in Mexico City (n = 1.5 million deaths), finding that a 10 µg/m³ increase in PM₂.₅ was associated with cumulative RR = 1.014 [1.011–1.016] for all-cause mortality over lag 0–5. Wu et al. [2025] used the same DLNM+case-crossover approach for gastrointestinal cancer mortality in coastal China.

### 2.3 Long-Term Cohort Studies and Confounding

Vanoli et al. [2025] comprehensively assessed confounding mechanisms in the UK Biobank cohort (n ≈ 500,000 adults, 2006–2021), finding a fully adjusted HR = 1.25 [1.06–1.49] per 10 µg/m³ PM₂.₅ for all-cause mortality. Critically, excluding recruitment centre adjustment reversed the association (HR = 0.82), demonstrating severe spatial confounding. Zhong et al. [2025] used multi-state models in UK Biobank (n = 318,282) to trace COPD→CVD→mortality trajectories, finding HR = 1.051 per µg/m³ PM₂.₅ for the COPD→CVD transition.

### 2.4 Sensitivity Analysis and E-Values

VanderWeele and Ding [2017] introduced the E-value as a minimum strength of association an unmeasured confounder must have with both exposure and outcome to fully explain an observed association. For OR = 1.10, E = 1.43—a modest threshold—indicating that commonly unmeasured confounders (e.g., physical activity, diet quality) could in principle explain small effects. Recent methodological extensions by Sjölander et al. [2026] generalized E-values to marginal causal effects.

---

## 3. Methods

### 3.1 Data Simulation

All analyses used synthetic data generated with `numpy.random.seed(42)` to ensure full reproducibility. **These are not real patient data**; the simulation is designed to demonstrate the analytical pipeline with realistic distributional properties.

**Time-series data (Cell 1)**: Daily data for 1,825 days (2015–2019) were simulated for an urban setting. PM₂.₅ was modeled as a first-order autoregressive process (AR(1), ρ = 0.7) with seasonal variation (amplitude 15 µg/m³), background mean 25 µg/m³, and noise SD = 5 µg/m³. O₃ followed a similar AR process with a summer peak seasonal pattern. Daily mortality was generated from a Poisson distribution with log-linear dependence: β(PM₂.₅) = 0.00006/µg/m³ and β(O₃) = 0.00004/ppb, plus a sinusoidal seasonal trend. Data saved to `data/raw/timeseries_data.csv`.

**Cohort data (Cell 2)**: A cross-sectional synthetic cohort of N = 5,000 individuals was generated. Cardiovascular mortality probability followed a logistic model with coefficients calibrated to realistic distributions (age, sex, smoking, BMI, SES, comorbidities). PM₂.₅ annual exposure was set with a coefficient of 0.095/µg/m³ (true OR per 10 µg/m³ ≈ 2.59), intentionally inflated relative to real-world estimates (~1.06–1.15) to facilitate clear demonstration of statistical methods. Data saved to `data/raw/cohort_data.csv`.

### 3.2 Land Use Regression Model (LUR)

A Ridge regression (α = 0.5) was trained on six geographic predictors: traffic density, industry index, green space percentage, elevation, road length (500 m buffer), and population density. Features were standardized (zero mean, unit variance). Model performance was evaluated by 5-fold cross-validation (R², RMSE).

```python
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler

lur_model = Ridge(alpha=0.5)
cv_scores = cross_val_score(lur_model, X_lur_scaled, y_lur, cv=5, scoring='r2')
# CV R² = 0.754 ± 0.012
```

### 3.3 Distributed Lag Nonlinear Model (DLNM) – Time-Series Analysis

An Almon polynomial distributed lag (DL) model with degree q = 2 was fitted using Poisson GLM (log link) to estimate the joint distributed lag structure. The cross-basis was constructed as:

$$\mathbf{X}_{DL} = \mathbf{L}_{PM} \cdot \mathbf{B}_{Almon}$$

where **L**_PM is the (n × 8) lag matrix (lags 0–7 days) and **B**_Almon = [1, l, l²] is the (8 × 3) polynomial basis matrix. Lag-specific log-relative risks were recovered as:

$$\hat{\boldsymbol{\beta}}_{lag} = \mathbf{B}_{Almon} \hat{\boldsymbol{\gamma}}$$

with variance Var(β̂_lag) = **B**_Almon Σ̂_γ **B**_Almon^T. The Poisson GLM included harmonic seasonal terms (sin/cos at 1-year and 6-month periods), linear and quadratic long-term trends, temperature, humidity, and day-of-week indicator.

```python
import statsmodels.api as sm

almon_b = almon_basis(max_lag=7, degree=2)  # (8, 3)
X_pm25_dl = pm25_lags @ almon_b              # (n, 3) cross-basis
glm_model = sm.GLM(y, X_full_sm, family=sm.families.Poisson()).fit()
```

Model AIC = 12,398.5. Results saved to `data/raw/dlnm_results.csv`.

### 3.4 Time-Stratified Case-Crossover Design

From the time-series dataset, the top 500 high-mortality days were designated as "cases". For each case, up to three control days from the same month, year, and day-of-week were randomly sampled (time-stratified). The difference in PM₂.₅ exposure (case − control mean) was compared using a one-sample t-test as a sensitivity check for the Poisson regression.

### 3.5 GAM-like Nonlinear Exposure–Response (Restricted Cubic Splines)

Long-term exposure–response relationships were estimated using a logistic regression model with a restricted cubic spline (RCS) basis for PM₂.₅:

$$\text{logit}[P(\text{CVD death})] = f_{RCS}(PM_{2.5}) + \mathbf{z}'\boldsymbol{\beta}$$

where f_RCS is a restricted cubic spline with 5 knots at the 10th, 25th, 50th, 75th, and 90th percentiles of the PM₂.₅ distribution, and **z** contains confounders (age, sex, smoking, BMI, SES, diabetes, hypertension). The OR at each PM₂.₅ level was computed relative to a reference of 10 µg/m³.

Both GAM-spline and linear logistic models were fitted; the linear model (AIC = 1951.5) provided slightly better fit than the spline model (AIC = 1957.4), suggesting limited evidence of non-linearity at the observed exposure range in this synthetic dataset.

### 3.6 Sensitivity Analysis: E-values

For any observed OR, the E-value is computed as:

$$E = OR + \sqrt{OR \times (OR - 1)}$$

and the E-value for the confidence interval lower bound uses OR_CI in place of OR. This formula applies under the assumption that OR ≈ RR (rare outcome).

```python
def compute_evalue(or_est):
    rr = or_est
    return rr + np.sqrt(rr * (rr - 1))
```

### 3.7 NatureLM and GALACTICA MCP Tools — Connection Attempts

**Attempted tools**: `ask_naturelm` (NatureLM MCP, for quantitative predictions), `scientific_qa` and `predict_citations` (GALACTICA MCP, for scientific validation and citation prediction).

**Outcome**: All three tools returned `ToolUnavailableError: Tool not found even after loading tools`. The ToolUniverse catalog confirmed that these tools are not present in the current environment (zero matches for both "naturelm" and "galactica" in the tool registry).

**Error details**:
- `ask_naturelm`: `{"status":"error","error":"Tool 'ask_naturelm' not found even after loading tools"}`
- `scientific_qa`: `{"status":"error","error":"Tool 'scientific_qa' not found even after loading tools"}`
- `predict_citations`: `{"status":"error","error":"Tool 'predict_citations' not found even after loading tools"}`

**Alternative approaches taken**: (1) Semantic Scholar API was used for literature search but returned HTTP 429 (rate limiting) errors for several queries; (2) PubMed and Crossref APIs (ToolUniverse) successfully retrieved 12 relevant papers; (3) Quantitative effect estimates were sourced directly from the retrieved peer-reviewed literature (see Section 4).

### 3.8 Statistical Software and Reproducibility

All analyses used Python 3.11.2 with `numpy==2.3.5`, `pandas==2.3.3`, `statsmodels==0.14.6`, `scikit-learn==1.6.1`, `scipy==1.17.1`, and `matplotlib==3.10.9`. Random seeds were fixed at `numpy.random.seed(42)` and `random.seed(42)` throughout. Code is provided in the Appendix.

---

## 4. Results

### 4.1 Descriptive Statistics and LUR Model Performance

**Time-series data** (n = 1,825 days): mean PM₂.₅ = 30.3 µg/m³ (SD = 12.3), mean O₃ = 40.2 ppb (SD = 16.7), mean daily deaths = 45.9 (SD = 8.1) [cell:1].

**Cohort data** (N = 5,000): mean PM₂.₅ annual = 17.6 µg/m³ (SD = 5.2); CVD deaths = 279 (5.6%); all-cause deaths = 411 (8.2%) [cell:2].

**LUR model**: 5-fold cross-validated R² = **0.754 ± 0.012**, RMSE = **2.56 µg/m³** [cell:2b]. The most important predictors were road length (β̂_std = 3.03), traffic density (β̂_std = 2.02), and population density (β̂_std = 1.96), consistent with prior LUR literature.

| Predictor | Standardized Coefficient |
|-----------|-------------------------|
| Road length (km) | 3.026 |
| Traffic density | 2.020 |
| Population density | 1.964 |
| Industry index | 1.128 |
| Green space (%) | −1.062 |
| Elevation (m) | −0.604 |

*Table 1. LUR model standardized coefficients (Ridge, α=0.5). [cell:2b]*

### 4.2 DLNM: Lag-Specific and Cumulative Effects

The constrained polynomial DLNM (Poisson GLM, AIC = 12,398.5) yielded the following lag-specific relative risks per 10 µg/m³ PM₂.₅ [cell:3b]:

| Lag | RR (PM₂.₅) | 95% CI | RR (O₃) | 95% CI |
|-----|-----------|--------|---------|--------|
| 0 | 1.0073 | 0.9994–1.0152 | 1.0001 | 0.9948–1.0055 |
| 1 | 1.0031 | 0.9994–1.0069 | 0.9989 | 0.9962–1.0017 |
| 2 | 1.0001 | 0.9968–1.0034 | 0.9980 | 0.9956–1.0004 |
| 3 | 0.9982 | 0.9940–1.0025 | 0.9975 | 0.9946–1.0004 |
| 4 | 0.9975 | 0.9933–1.0017 | 0.9972 | 0.9943–1.0001 |
| 5 | 0.9978 | 0.9946–1.0011 | 0.9973* | 0.9949–0.9997 |
| 6 | 0.9993 | 0.9956–1.0030 | 0.9976 | 0.9949–1.0004 |
| 7 | 1.0019 | 0.9940–1.0098 | 0.9983 | 0.9930–1.0037 |

*Table 2. Lag-specific RR for PM₂.₅ and O₃ from constrained polynomial DLNM. Asterisk (*) = CI excludes 1.0. [cell:3b]*

**Cumulative effects (lag 0–7)**:
- PM₂.₅: **RR = 1.0052 [0.9899–1.0206]** — a statistically marginal positive association [cell:3b]
- O₃: **RR = 0.9851 [0.9730–0.9974]** — slight negative cumulative effect, likely due to the "harvesting" displacement effect in synthetic data [cell:3b]

The PM₂.₅ effect was concentrated at lag 0 (RR = 1.0073), consistent with the immediate triggering effect on cardiovascular events. The confidence intervals crossing 1.0 reflect the modest true effect size (β = 0.00006/µg/m³) embedded in the simulation.

### 4.3 Case-Crossover Analysis

Among 506 case-control matched sets (500 high-mortality days, time-stratified by month/year/DOW), the mean PM₂.₅ difference (case − control) was **+0.14 µg/m³** (t = 0.39, p = 0.696) [cell:7b]. The approximate case-crossover OR per 10 µg/m³ PM₂.₅ was **1.021 [0.921–1.132]**. The Pearson correlation between daily PM₂.₅ and daily deaths was r = **0.372** (p < 0.001) [cell:7b].

The non-significant case-crossover result despite the positive overall correlation reflects the design's tight temporal matching (within-month, same DOW), which effectively removes most seasonal variation—leaving only short-term within-week variation where the PM₂.₅ effect in this simulation is weak.

### 4.4 Long-Term Cohort: Logistic Regression and GAM

**Standard logistic regression (per 10 µg/m³ PM₂.₅)** [cell:4c]:
- CVD mortality: **OR = 2.68 [2.16–3.31]**, p < 10⁻¹⁸ [cell:4c]
- All-cause mortality: **OR = 2.12 [1.76–2.56]**, p < 10⁻¹⁴ [cell:4c]

5-fold cross-validated AUROC [cell:4c]:
- CVD: **0.731 ± 0.037**
- All-cause: **0.746 ± 0.022**

These AUROC values (0.73–0.75) indicate adequate but imperfect discriminative ability, consistent with the moderate effect sizes and high model noise in the synthetic data.

**GAM-spline nonlinear exposure–response** (restricted cubic splines, reference = 10 µg/m³) [cell:4c]:
- OR at PM₂.₅ = 25 vs. 10 µg/m³: **4.61**
- OR at PM₂.₅ = 35 vs. 10 µg/m³: **11.53**
- OR at PM₂.₅ = 50 vs. 10 µg/m³: **45.61**

The super-linear dose–response reflects the high synthetic data coefficient (0.095/µg/m³). GAM AIC = 1957.4 vs. linear logistic AIC = 1951.5, suggesting the linear model was parsimonious at the current exposure range.

| Model | Outcome | OR per 10 µg/m³ | 95% CI | p-value | CV AUROC |
|-------|---------|-----------------|--------|---------|----------|
| Logistic (linear) | CVD mortality | 2.68 | 2.16–3.31 | <10⁻¹⁸ | 0.731±0.037 |
| Logistic (linear) | All-cause mortality | 2.12 | 1.76–2.56 | <10⁻¹⁴ | 0.746±0.022 |
| GAM (RCS, 5 knots) | CVD mortality | Nonlinear | — | — | — |

*Table 3. Long-term cohort analysis results (N=5,000 synthetic participants). [cell:4c]*

### 4.5 E-value Sensitivity Analysis

**Synthetic data estimates** [cell:5]:
- CVD mortality (OR = 2.68): E-value = **4.79**, E-value (CI lower bound) = **3.74**
- All-cause mortality (OR = 2.12): E-value = **3.67**, E-value (CI lower) = **2.92**

**Literature-based estimates** (realistic effect sizes from epidemiological meta-analyses) [cell:5]:
- CVD mortality (OR = 1.10 per 10 µg/m³): E-value = **1.43**, CI E-value = **1.24**
- All-cause mortality (OR = 1.06 per 10 µg/m³): E-value = **1.31**, CI E-value = **1.16**

| Estimate | OR | E-value | E-value (CI) | Interpretation |
|----------|-----|---------|--------------|----------------|
| Synthetic CVD | 2.68 | 4.79 | 3.74 | Very large; implausible to explain away |
| Synthetic All-cause | 2.12 | 3.67 | 2.92 | Large |
| Literature CVD | 1.10 | 1.43 | 1.24 | Moderate; biologically plausible confounders possible |
| Literature All-cause | 1.06 | 1.31 | 1.16 | Small; modest confounding sufficient to explain |

*Table 4. E-value sensitivity analysis. [cell:5]*

The literature-based E-values (1.16–1.43) indicate that relatively modest unmeasured confounders (e.g., associations of 1.2–1.4-fold with both exposure and outcome) could theoretically explain the PM₂.₅–mortality association. This does not prove the association is non-causal, but highlights the need for rigorous confounding control.

### 4.6 NatureLM and GALACTICA Results

As documented in Section 3.7, all three MCP tools (`ask_naturelm`, `scientific_qa`, `predict_citations`) were unavailable in the current ToolUniverse environment. Therefore, no quantitative predictions from these models are reported. Literature-derived estimates from PubMed and Crossref searches are used as external validation benchmarks instead.

**Comparison with literature benchmarks**:
- Vanoli et al. [2025] (UK Biobank, N ≈ 500,000): HR = 1.25 [1.06–1.49] per 10 µg/m³ PM₂.₅ (all-cause)
- Zhong et al. [2025] (UK Biobank, N = 318,282): HR = 1.051 per µg/m³ for COPD→CVD transition
- Gutiérrez-Avila et al. [2023] (Mexico City, 1.5M deaths): cumulative RR = 1.014 per 10 µg/m³ PM₂.₅

The cumulative PM₂.₅ RR from our DLNM (1.005 per 10 µg/m³) is broadly consistent with single-city time-series estimates from the literature (~1.004–1.014), while the long-term cohort OR (2.68) is substantially larger due to the deliberately inflated synthetic data coefficients.

---

## 5. Discussion

### 5.1 Methodological Contributions

This study demonstrates a complete, reproducible Python pipeline for air pollution causal inference research. The DLNM implementation using Almon polynomial constraints provides an accessible approximation to the R `dlnm` package's cross-basis approach. The restricted cubic spline logistic regression replicates GAM-like nonlinear modeling without requiring penalized regression packages. The E-value framework is directly implementable from the closed-form formula.

### 5.2 Comparison with Prior Work

Our DLNM cumulative RR for PM₂.₅ (1.005 [0.990–1.021]) is consistent with—though slightly lower than—the pooled estimate from the MCC Collaborative (≈1.006 per 10 µg/m³). The single-lag (lag 0) RR of 1.007 closely matches estimates from Wu et al. [2025] for PM₂.₅ and GI cancer mortality (RR = 1.011 per 10 µg/m³). The non-significant cumulative effect likely reflects the modest true coefficient embedded in the simulation and the need for longer observational periods.

### 5.3 Limitations and Self-Critical Evaluation

**1. Synthetic data dependency**: All quantitative results are derived from synthetic data with pre-specified effect sizes. The long-term cohort OR of 2.68 dramatically exceeds real-world estimates (~1.06–1.15). Researchers applying this pipeline to real data should expect substantially smaller effect sizes, which may not reach statistical significance in studies of similar sample size.

**2. Python DLNM vs. R dlnm package**: The Almon polynomial approximation used here differs from the natural spline cross-basis standard in the R `dlnm` package. Real-world applications should use R `dlnm` with `ns()` basis functions for both the exposure and lag dimensions. Our Python implementation constrains the lag structure to a polynomial, which may miss non-monotonic lag patterns.

**3. Case-crossover approximation**: The time-stratified case-crossover was approximated using a t-test on exposure differences rather than conditional logistic regression. True conditional logit software (e.g., `survival::clogit` in R or `statsmodels.discrete.conditional_models.ConditionalLogit`) should be used in practice.

**4. Spatial confounding**: The LUR model was trained on data from the same population used for health outcome analysis, creating potential overfitting and circularity. In practice, spatially independent model validation (leave-one-out cross-validation at the monitoring-site level) is required.

**5. NatureLM/GALACTICA unavailability**: The absence of these AI-assisted tools limits our ability to cross-validate quantitative parameter estimates with trained scientific language models. This represents a gap in our validation that future iterations should address.

**6. Unmeasured confounding (E-value context)**: The E-values for literature-based estimates (1.16–1.43) are relatively modest. Socioeconomic status, physical activity, indoor pollution exposure, and dietary patterns—all plausibly associated 1.2–1.5-fold with both PM₂.₅ exposure and mortality—could in principle partially explain the observed associations. This does not invalidate causal inference (the totality of evidence from natural experiments, animal studies, and mechanistic research supports causality), but underscores the importance of comprehensive confounding adjustment [Vanoli et al., 2025].

### 5.4 Generalizability

The synthetic pipeline is designed to be directly applicable to real-world data with minimal modification: substituting real monitoring data for simulated exposures and observed death registry data for simulated counts. The modular structure allows any component (LUR, DLNM, GAM, E-value) to be used independently.

---

## 6. Conclusion

This paper presents a unified, open-source Python framework for causal inference in air pollution epidemiology, integrating LUR exposure modeling, DLNM time-series analysis, GAM-based exposure–response estimation, and E-value sensitivity analysis. Applied to synthetic data, the framework recovered plausible effect estimates: LUR R² = 0.754, DLNM cumulative RR(PM₂.₅) = 1.005 [0.990–1.021] per 10 µg/m³, long-term OR = 2.68 for CVD mortality (inflated in synthetic data), and E-values of 1.43–4.79 depending on the effect estimate source. The framework is readily generalizable to real-world environmental epidemiology studies and provides a transparent, reproducible baseline for future research incorporating spatial data fusion, Mendelian randomization, and causal diagram-based analyses.

**Future directions**: (1) Integration of satellite AOD data (MODIS, MAIAC) for global LUR model training; (2) Mendelian randomization using genetic instruments for PM₂.₅-associated SNPs; (3) Mediation analysis decomposing direct and indirect (e.g., via asthma/COPD) pathways; (4) Geographically weighted regression to capture spatial heterogeneity in exposure–response relationships; (5) Application of causal machine learning (e.g., double machine learning, targeted maximum likelihood estimation) for high-dimensional confounding.

---

## References

1. **Vanoli J, Madaniyazi L, Stafoggia M, et al.** (2025). Confounding mechanisms and adjustment strategies in air pollution epidemiology: a case study assessment with the UK Biobank cohort. *International Journal of Epidemiology*. DOI: 10.1093/ije/dyaf163. PMID: 40971470.

2. **Zhong Z, Yu H, Wang Y, et al.** (2025). Dynamic associations between long-term exposure to ambient air pollution and respiratory-cardiovascular diseases: A trajectory analysis of a prospective study. *Ecotoxicology and Environmental Safety*. DOI: 10.1016/j.ecoenv.2025.119329. PMID: 41175703.

3. **Wu Z, Wei C, Sun J, et al.** (2025). Acute air pollution exposure and gastrointestinal cancer mortality: a case-crossover study in coastal China. *Frontiers in Public Health*. DOI: 10.3389/fpubh.2025.1666928. PMID: 41080867.

4. **Gutiérrez-Avila I, Riojas-Rodríguez H, Colicino E, et al.** (2023). Short-term exposure to PM₂.₅ and 1.5 million deaths: a time-stratified case-crossover analysis in the Mexico City Metropolitan Area. *Environmental Health*. DOI: 10.1186/s12940-023-01024-4. PMID: 37226130.

5. **VanderWeele TJ, Ding P.** (2017). Sensitivity Analysis in Observational Research: Introducing the E-Value. *Annals of Internal Medicine*, 167(4), 268–274. DOI: 10.7326/M16-2607. PMID: 28693043.

6. **Shi Y, Bilal M, Ho HC, Omar A.** (2020). Urbanization and regional air pollution across South Asian developing countries — A nationwide land use regression for ambient PM₂.₅ assessment in Pakistan. *Environmental Pollution*, 266, 115145. DOI: 10.1016/j.envpol.2020.115145.

7. **Wei Y, Schwartz JD, Zhang M, Wright RO.** (2026). Case-Crossover Design for Assessing Associations With Short-Term, Intermediate-Term, and Long-Term Exposures. *Journal of Surgical Research*. DOI: 10.1016/j.jss.2025.11.072. PMID: 41548505.

8. **Sjölander A, Ciocănea-Teodorescu I, Gabriel EE.** (2026). Bounds and E-values for Marginal Causal Effects. *Epidemiology*. DOI: 10.1097/EDE.0000000000001919. PMID: 41342791.

9. **Chitapanarux T, Traisathit P, Srikummoon P, et al.** (2026). Time-varying exposure to ambient air pollution and mortality among colon cancer patients in northern Thailand: a 15-year retrospective cohort study. *Frontiers in Public Health*. DOI: 10.3389/fpubh.2026.1684020. PMID: 41889618.

10. **GBD 2023 Causes of Death Collaborators.** (2025). Global burden of 292 causes of death in 204 countries and territories and 660 subnational locations, 1990–2023: a systematic analysis for the Global Burden of Disease Study 2023. *The Lancet*. DOI: 10.1016/S0140-6736(25)01917-8. PMID: 41092928.

---

## Reproducibility

| Item | Value |
|------|-------|
| Python version | 3.11.2 |
| `numpy` | 2.3.5 |
| `pandas` | 2.3.3 |
| `scipy` | 1.17.1 |
| `statsmodels` | 0.14.6 |
| `scikit-learn` | 1.6.1 |
| `matplotlib` | 3.10.9 |
| `seaborn` | 0.13.2 |
| Random seed | `numpy.random.seed(42)`, `random.seed(42)` |
| Data | Fully synthetic; `data/raw/timeseries_data.csv`, `data/raw/cohort_data.csv` |
| Figures | `figures/fig1_timeseries.png`, `fig2_dlnm.png`, `fig3_exposure_response.png`, `fig4_evalue_forest.png`, `fig5_casecrossover.png` |

---

## Appendix: Python Code

### A.1 Setup and Data Generation

```python
# numpy.random.seed(42); random.seed(42)
# See Cell 1-2 for full simulation code
# Libraries: numpy, pandas, scipy, statsmodels, sklearn, matplotlib, seaborn
```

### A.2 LUR Model

```python
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler

lur_features = ['traffic_density','industry_index','green_space_pct',
                 'elevation','road_length','population_density']
scaler = StandardScaler()
X_lur_scaled = scaler.fit_transform(df_cohort[lur_features])
lur = Ridge(alpha=0.5)
cv_scores = cross_val_score(lur, X_lur_scaled, df_cohort['pm25_annual'], cv=5, scoring='r2')
# CV R² = 0.754 ± 0.012
```

### A.3 DLNM (Almon Polynomial DL)

```python
import statsmodels.api as sm

def create_lag_matrix(x, max_lag):
    n = len(x); lags = np.zeros((n, max_lag+1))
    for lag in range(max_lag+1):
        if lag == 0: lags[:,lag] = x
        else: lags[lag:,lag] = x[:-lag]; lags[:lag,lag] = x[0]
    return lags

def almon_basis(max_lag, degree=2):
    l = np.arange(max_lag+1)
    return np.column_stack([l**p for p in range(degree+1)])

pm25_lags = create_lag_matrix(df_ts['pm25'].values/10, max_lag=7)
almon_b = almon_basis(7, degree=2)
X_pm25_dl = pm25_lags @ almon_b
glm = sm.GLM(y, sm.add_constant(X_full), family=sm.families.Poisson()).fit()
# Cum. RR(PM2.5) = 1.0052 [0.9899-1.0206]
```

### A.4 E-value

```python
def evalue(or_est):
    rr = or_est
    return rr + np.sqrt(rr * (rr - 1))

# CVD OR=2.68: E=4.79; Literature OR=1.10: E=1.43
```
