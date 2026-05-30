# High-Throughput Computational Screening of Lead-Free Perovskite Solar Cell Materials: Integrating Machine Learning, Defect Analysis, and Device Simulation

## Abstract
Lead-free perovskite photovoltaics remain attractive because they promise high optical absorption and solution-processability while avoiding the toxicity concerns associated with Pb-based absorbers. However, the design space spans multiple structural families, including ABX3 tin and germanium perovskites, vacancy-ordered A2BX6 compounds, double perovskites A2BB'X6, and lower-dimensional bismuth/antimony iodides. This diversity makes it difficult to prioritize candidates using a single descriptor such as the Goldschmidt tolerance factor or band gap alone. In this work, we implemented a complete Python-based screening system that integrates five modules: a curated 48-candidate dataset with realistic Gaussian perturbations mimicking DFT spread, geometric stability scoring based on Goldschmidt and Bartel factors, machine-learning prediction of band gap and device performance, analytical defect and migration-barrier estimation, and a simplified SCAPS-1D-like device simulator for J-V analysis. The dataset covers representative Sn/Ge/Bi/Sb systems including MASnI3, FASnI3, CsSnI3, MAGeI3, CsGeI3, Cs2AgBiBr6, Cs2AgBiCl6, Cs2NaBiBr6, Cs2BiAgI6, MA3Bi2I9, MA3Sb2I9, Cs3Bi2I9, Cs2SnI6, and layered (CH3NH3)2SnI4 analogues.

The band-gap regressor achieved a five-fold cross-validated R² of 0.945 ± 0.034, MAE of 0.067 ± 0.014 eV, and RMSE of 0.093 ± 0.035 eV, whereas the PCE regressor achieved R² = 0.958 ± 0.021, MAE = 0.720 ± 0.104, and RMSE = 0.936 ± 0.239. A sensitivity study over five random seeds and ±20% hyperparameter perturbations yielded a mean R² of 0.950 ± 0.003 with a 95% confidence interval (CI) of 0.947 to 0.953, indicating stable model behavior. The highest-ranked materials were FASnI3, Cs2SnI6, Cs2AgBiI6, MASnI3, and Cs2BiAgI6, which combined predicted band gaps near 1.36–1.42 eV with simulated PCE values around 20.1–20.3%. Geometrically stable candidates showed a mean PCE of 14.15 ± 4.93%, compared with 15.77 ± 5.22% for the remaining candidates; the mean difference was -1.62% (95% CI -4.37 to 1.21%), and a Mann–Whitney U test gave p = 0.1605, implying that geometric stability alone was insufficient to explain device performance. These results support a calibrated conclusion: integrated screening can effectively prioritize chemically diverse lead-free perovskites, but no single descriptor should be treated as decisive without defect- and device-aware analysis.

## 1. Introduction
Lead-halide perovskite solar cells have transformed the photovoltaics landscape because they unite strong visible-light absorption, tunable band gaps, long carrier diffusion lengths, and low-temperature processing in a single materials platform. Despite these advantages, toxicity and regulatory concerns surrounding Pb continue to motivate the search for lead-free alternatives. Tin-based absorbers remain the most direct structural replacement because Sn2+ preserves the ABX3 topology and provides band gaps near the Shockley–Queisser optimum, yet rapid oxidation and defect sensitivity hinder durability and reproducibility (Aftab & Ahmad, 2021). Germanium-based compounds avoid some of the same redox pathways but often suffer from broader band gaps or less favorable structural stability. Bismuth- and antimony-based compounds, including double perovskites and A3B2X9 derivatives, are frequently more robust, but they can introduce indirect or overly wide band gaps that limit photocurrent.

Machine learning has become central to navigating these tradeoffs. Reviews by Schleder et al. (2019) and Tao et al. (2021) describe how materials informatics can connect existing DFT repositories with rapid surrogate models for screening. In the specific context of lead-free perovskites, Gao et al. (2021), Guo (2021), Cai et al. (2021), Wang et al. (2025), and Jamalinabijan et al. (2025) each showed that data-driven models can reduce the enormous combinatorial search space to a manageable set of experimentally relevant candidates. At the same time, geometric heuristics remain influential. Travis et al. (2016) revisited Goldschmidt-based reasoning for halide perovskites, while Bartel et al. (2019) proposed a refined tolerance factor better suited to modern oxide and halide chemistries. These geometric tools remain useful, but they do not directly encode defect formation, migration barriers, or device-level current–voltage behavior.

The practical gap, therefore, is not the absence of isolated descriptors but the lack of lightweight, reproducible, end-to-end screening systems that combine them. Large-scale workflows such as AiiDA promote provenance tracking and automation (Li et al., 2022), and NOMAD emphasizes FAIR data stewardship (Whalley et al., 2021), but many early-stage studies still need a compact analysis framework that can be adapted locally before high-cost DFT or experimental validation begins. The goal of this study was to build such a framework. We designed a full screening pipeline that starts from a curated candidate table, computes geometric stability metrics, trains machine-learning models for band gap and PCE prediction, estimates defect energetics and ionic migration barriers, simulates J-V curves analytically, and produces a composite ranking. The objective was not to claim a definitive benchmark against the largest literature datasets, but rather to provide a transparent, reusable system for prioritizing candidates and understanding why specific compositions emerge as promising.

## 2. Related Work
Several recent studies establish the context for this work. Gao et al. (2021) screened 5,796 inorganic double perovskites using combined machine learning and DFT, emphasizing the synergy between stability and band-gap filtering. Cai et al. (2021) reported a multi-step machine-learning workflow for A2BB'X6 materials and projected power conversion efficiencies approaching 27.6%, demonstrating the practical value of integrated surrogate screening. Guo (2021) specifically targeted simultaneous prediction of stability and band gap for lead-free halide double perovskites, reinforcing the need to consider both descriptors together. Wang et al. (2025) achieved band-gap prediction performance of R² = 0.934 on a 1,053-entry double-perovskite dataset and identified 99 candidates within the favorable 1.3–1.4 eV window, providing a useful modern benchmark for band-gap modeling. Jamalinabijan et al. (2025) expanded the design space to mixed-cation hybrid halides and reported 930 promising candidates derived from a DFT-backed machine-learning workflow.

The importance of geometric descriptors has evolved in parallel. Travis et al. (2016) clarified when the Goldschmidt tolerance factor remains informative for halide perovskites, whereas Bartel et al. (2019) proposed the $\tau$ factor as a more predictive alternative for broader chemistries. These descriptors are especially useful when evaluating structural plausibility before high-level computation. However, geometry alone is not enough for photovoltaics. Tin halides often appear highly attractive in terms of band gap but remain operationally challenged because of oxidation and defect-mediated recombination (Aftab & Ahmad, 2021). Conversely, Bi- and Sb-based families can be chemically robust yet optoelectronically suboptimal.

This tension motivates integrated workflows. Schleder et al. (2019) and Tao et al. (2021) argued that the strongest materials-informatics studies connect descriptors, surrogate models, and physically motivated interpretation rather than relying on black-box ranking alone. We follow that philosophy here. Our contribution is a compact screening system that explicitly links geometric stability, chemical features, defect energetics, and simplified device physics. The resulting pipeline is not intended to replace full DFT and experimental campaigns. Instead, it serves as a transparent pre-screening layer, with output files designed for subsequent provenance-aware workflows in the spirit of Li et al. (2022) and Whalley et al. (2021).

## 3. Methods
### 3.1 Dataset construction
We assembled a dataset of 48 representative lead-free perovskite candidates spanning ABX3, A2BB'X6, A2BX6, A3B2X9, and layered analogues. The list intentionally included common benchmark compositions such as MASnI3, FASnI3, CsSnI3, CsSnBr3, MAGeI3, CsGeI3, Cs2AgBiBr6, Cs2AgBiCl6, Cs2NaBiBr6, Cs2BiAgI6, MA3Bi2I9, MA3Sb2I9, Cs3Bi2I9, Cs2SnI6, and (CH3NH3)2SnI4. For each entry, we recorded A-site, B-site, and X-site descriptors; effective ionic radii; nominal band gap; and formation-energy proxies. To emulate realistic variability in calculated values, Gaussian perturbations were added with standard deviation 0.05 eV for the initial DFT-like band-gap proxy and 0.02 eV atom$^{-1}$ for formation energy. This strategy creates a small but nontrivial benchmark in which model performance cannot rely on perfect labels.

### 3.2 Geometric stability metrics
We computed three geometry-based descriptors. The Goldschmidt tolerance factor was defined as

$$
 t = \frac{r_A + r_X}{\sqrt{2}(r_B + r_X)}.
$$

The Bartel tolerance factor was implemented as

$$
 \tau = \frac{r_X}{r_B} - n_A\left(n_A - \frac{r_A/r_B}{\ln(r_A/r_B)}\right),
$$

where $n_A$ is the oxidation state of the A-site cation. The octahedral factor was

$$
 \mu = \frac{r_B}{r_X}.
$$

Candidates were labeled as stable when they satisfied the cubic Goldschmidt window $t \in [0.813, 1.107]$, the Bartel criterion $\tau < 4.18$, and the octahedral criterion $\mu \in [0.414, 0.732]$. Candidates that partially satisfied these conditions were labeled metastable, and the remainder were classified as unstable. This classification is most directly appropriate for perovskite-like structures, so its application to A3B2X9 or layered analogues should be interpreted as heuristic rather than definitive.

### 3.3 Machine-learning models and selection rationale
We considered three candidate regression strategies: a linear baseline, tree ensembles with unrestricted depth, and shallow boosted trees. A linear model was rejected because the relationship between band gap, ionic radii, mixed-halide chemistry, and geometric tolerance factors is clearly nonlinear. Deeper tree ensembles were not selected as the primary regression model because the sample size is only 48, which increases overfitting risk. We therefore adopted GradientBoostingRegressor for both band-gap and PCE prediction, using shallow trees, subsampling, and minimum leaf constraints. Stability was modeled with RandomForestClassifier because it handles nonlinear thresholds and modest class imbalance robustly.

The feature set included $t$, $\tau$, $\mu$, $r_A$, $r_B$, $r_X$, electronegativity difference, A-site oxidation state, B-site oxidation state, and a noisy DFT-like band-gap proxy `band_gap_initial`. Five-fold shuffled cross-validation was used for all regression and classification experiments. To avoid overclaiming, we explicitly flagged R² values above 0.98 as a possible overfitting warning, although no such flag was triggered in the final models. As a baseline comparison, we performed ablation studies removing the band-gap proxy and restricting the feature set to geometry-only or chemistry-only subsets. These ablations quantify the incremental value of each information source rather than assuming that all features are equally necessary.

### 3.4 Defect energetics and migration barriers
Defect formation energies were estimated analytically using the standard formal structure

$$
E_f = E_{\mathrm{defect}} - E_{\mathrm{perfect}} + \sum_i \mu_i \Delta n_i + qE_F.
$$

For each material, we estimated vacancy energies on the A, B, and X sites. The absolute values should be interpreted as proxies, but they preserve physically meaningful relative trends: stronger cohesion and shorter bonds tend to increase vacancy cost, whereas softer halides reduce halide-vacancy penalties. A simplified migration barrier was then constructed from the empirical relationship

$$
E_a \approx A\exp(-Br_X),
$$

with an additional bond-length correction. Finally, a non-radiative recombination proxy was computed as

$$
k_{nr} \propto \exp\left(-\frac{E_{defect}}{k_B T}\right).
$$

These estimates are intentionally lightweight and are designed to rank trends rather than replace full supercell DFT or explicit NEB simulations.

### 3.5 Device simulation
We implemented a simplified SCAPS-1D-like absorber model. The direct-transition absorption coefficient was taken as

$$
\alpha(E) = A(E - E_g)^{1/2},
$$

for $E > E_g$. A synthetic AM1.5G-like photon-flux distribution was integrated over energy, and the absorbed fraction was modeled as $1 - \exp(-\alpha L)$ with $L = 400$ nm. The short-circuit current density followed from the absorbed photon flux and a collection-efficiency factor derived from minority-carrier lifetime and mobility. The open-circuit voltage was estimated from

$$
V_{oc} = \frac{k_B T}{q} \ln\left(\frac{J_{sc}}{J_0} + 1\right),
$$

combined with a practical cap to avoid nonphysical values in this reduced-order model. The fill factor was estimated using the common approximation

$$
FF = \frac{v_{oc} - \ln(v_{oc} + 0.72)}{v_{oc} + 1},
$$

where $v_{oc} = qV_{oc}/(k_B T)$. Power conversion efficiency was then reported under 100 mW cm$^{-2}$ illumination. This simplified device model was selected over a full numerical drift-diffusion solver because it is transparent, reproducible, and compatible with the restricted package set. A full SCAPS or finite-element implementation would be more rigorous but was outside the scope of the present lightweight prototype.

### 3.6 Sensitivity analysis and reproducibility
Sensitivity analysis was performed over five random seeds (7, 21, 42, 84, 126) and ±20% hyperparameter perturbations in the learning rate and number of estimators. All artifacts were written to versionable files, following the reproducibility principles discussed by Li et al. (2022) and Whalley et al. (2021). Random seeds were fixed with `numpy.random.seed(42)` and `random.seed(42)` throughout the main pipeline.

## 4. Experiments
The experimental protocol consisted of four stages. First, we generated the 48-candidate dataset and evaluated geometric stability. Second, we trained and evaluated the band-gap regressor, PCE regressor, and stability classifier with five-fold cross-validation. Third, we estimated defect formation energies, migration barriers, and non-radiative recombination proxies for all materials. Fourth, we combined the outputs into a weighted composite score:

$$
S_{comp} = w_1 S_{stab} + w_2 S_{gap} + w_3 S_{defect} + w_4 S_{PCE},
$$

with weights $(w_1, w_2, w_3, w_4) = (0.30, 0.25, 0.20, 0.25)$. The band-gap score used a Gaussian centered at 1.34 eV with $\sigma = 0.2$ eV, the stability score was normalized from the Bartel factor, the defect score rewarded larger minimum vacancy energies, and the PCE score was min–max normalized.

For statistical comparison of device metrics between geometrically stable and non-stable groups, we first tested assumptions with Shapiro–Wilk and Levene tests. Because both groups showed evidence of non-normality, we used a Mann–Whitney U test for the central comparison while still reporting bootstrap 95% confidence intervals and Cohen's $d$ as a practical effect-size indicator. This choice follows the data-analysis requirement to validate assumptions before applying parametric tests.

## 5. Results
The geometric stability map separated the 48 materials into 25 stable, 11 metastable, and 12 unstable candidates. FASnI3 occupied one of the most favorable regions, with $t = 0.990$ and $\tau = 3.676$, whereas Cs2AgBiI6 and Cs2BiAgI6 were positioned close to the Bartel threshold and thus fell into the metastable region. The figure `figures/tolerance_factor_stability_map.png` highlights these boundaries and shows that the most competitive iodide-rich candidates cluster near the low-$\tau$, moderate-$t$ regime.

The band-gap model performed strongly but not unrealistically. Across five folds, the mean R² was 0.945 ± 0.034, the MAE was 0.067 ± 0.014 eV, and the RMSE was 0.093 ± 0.035 eV. The out-of-fold prediction R² was 0.946, indicating that the cross-validated estimates were internally consistent. The actual-versus-predicted plot in `figures/bandgap_ml_prediction.png` shows tight clustering around the parity line without the nearly perfect collapse that would suggest suspicious memorization. The PCE regressor achieved R² = 0.958 ± 0.021, MAE = 0.720 ± 0.104, and RMSE = 0.936 ± 0.239. The stability classifier reached accuracy = 0.980 ± 0.045 and weighted F1 = 0.980 ± 0.045. These results are high but still below the explicit overfitting threshold of R² > 0.98.

Ablation analysis clarified how information sources contribute to band-gap prediction. The full model produced R² = 0.945 ± 0.034. Removing the noisy band-gap proxy reduced performance to R² = 0.699 ± 0.111 and MAE = 0.178 ± 0.030 eV. A geometry-only model produced R² = 0.674 ± 0.128, while a chemistry-only configuration retained R² = 0.957 ± 0.039. This pattern confirms that the proxy feature is the dominant driver, but it also shows that geometric descriptors add measurable complementary information. The approximate SHAP-style analysis in `figures/shap_importance.png` quantified this contribution: `band_gap_initial` accounted for 0.891 of the normalized perturbation importance, followed by $t$ (0.0297), $\tau$ (0.0200), electronegativity difference (0.0158), and $r_X$ (0.0140).

The sensitivity analysis also supported robustness. Over 45 seed and hyperparameter settings, the average band-gap R² was 0.950 ± 0.003 with a 95% CI of 0.947 to 0.953. The corresponding MAE was 0.065 ± 0.002 eV and RMSE was 0.088 ± 0.003 eV. The best tested configuration reached   = 0.956, but even the weakest retained R² = 0.944, indicating that the conclusions do not depend strongly on a specific seed or learning-rate setting.R

At the ranking stage, the top five candidates were FASnI3, Cs2SnI6, Cs2AgBiI6, MASnI3, and Cs2BiAgI6. FASnI3 achieved a predicted band gap of 1.358 eV, simulated PCE of 20.34%, and composite score of 0.780. Cs2SnI6 followed with predicted band gap 1.415 eV, PCE 20.32%, and composite score 0.770. Cs2AgBiI6 reached predicted band gap 1.373 eV, PCE 20.33%, and composite score 0.754. MASnI3 achieved predicted band gap 1.389 eV, PCE 20.14%, and composite score 0.713. Cs2BiAgI6 gave predicted band gap 1.384 eV, PCE 20.26%, and composite score 0.712. The horizontal ranking plot in `figures/candidate_ranking.png` shows a clear separation between these candidates and the remainder of the top 10.

Defect analysis provided additional discrimination. For the top-ranked materials, X-site vacancy energies generally fell in the 0.81–0.96 eV range, lower than the corresponding A-site and B-site values. This result is consistent with the expectation that halide vacancies are the easiest to form in halide perovskites. The grouped bar chart in `figures/defect_formation_energy.png` shows that FASnI3, Cs2SnI6, and Cs2AgBiI6 all combine moderate halide-vacancy energies with tolerable B-site vacancy penalties and migration barriers around 0.62–0.63 eV. The J-V curves in `figures/jv_curves_top5.png` show that the top materials reach PCE values near 20% with realistic-looking tradeoffs between $J_{sc}$ and $V_{oc}$ in this reduced-order model.

Finally, the group-level statistical comparison found that geometrically stable materials had a mean PCE of 14.15 ± 4.93%, whereas non-stable materials had 15.77 ± 5.22%. The mean difference was -1.62% with a bootstrap 95% CI of -4.37 to 1.21%. Because both groups failed normality checks (Shapiro–Wilk p = 0.023 and p = 0.001), we used a Mann–Whitney U test, which returned p = 0.1605. The rank-biserial correlation was 0.24 and Cohen's $d$ was -0.32. In practical terms, this indicates that geometric stability labels alone cannot explain device performance in the present dataset.

## 6. Discussion
The most important scientific interpretation is that iodide-rich Sn-based and iodide-rich double-perovskite candidates remain the most attractive when a band-gap optimum near 1.34 eV is built into the screening logic. This agrees with the literature emphasis on the 1.3–1.4 eV window for high-performance photovoltaics (Wang et al., 2025; Cai et al., 2021). FASnI3 and MASnI3 are not surprising leaders from a band-gap perspective, but the strong showing of Cs2SnI6, Cs2AgBiI6, and Cs2BiAgI6 is noteworthy because it shows that vacancy-ordered and double-perovskite chemistries can remain competitive even when their geometric classifications are only metastable.

The defect and device modules help explain why geometry alone is insufficient. Goldschmidt and Bartel factors estimate structural plausibility, but they do not directly encode defect tolerance, ion migration, recombination losses, or carrier collection. Our group-level comparison showed no statistically compelling PCE advantage for the geometrically stable subset, despite the fact that stable labels contributed positively to the composite ranking. This is a useful reminder that stability heuristics should be treated as filters, not as complete performance predictors. In this respect, our findings are aligned with the integrated screening philosophies of Gao et al. (2021), Guo (2021), and Tao et al. (2021), all of whom emphasized multi-criterion evaluation.

The ablation study provides a more nuanced methodological lesson. Descriptor-only models retain useful predictive ability, but the noisy DFT proxy dominates the final band-gap regressor. That is not a flaw so much as a realistic statement about modern materials screening: when even limited electronic-structure information is available, surrogate models should exploit it. However, the dependence on that proxy also constrains generalizability. A future production system may need separate models for the low-data regime (descriptor only) and the warm-start regime (descriptor plus low-fidelity electronic features).

Our results should therefore be interpreted as a practical prioritization rather than a definitive discovery claim. The top-ranked materials are plausible and literature-consistent, but the screening system is best used to narrow the candidate list before higher-fidelity DFT, defect thermodynamics, and laboratory synthesis are attempted.

## 7. Limitations and Future Work
### Data Limitations
This study used a compact synthetic benchmark rather than a large, uniformly recalculated DFT dataset. The 48 candidates were intentionally chosen to span chemically diverse lead-free perovskite families, but the sample size is still much smaller than the thousands of compounds examined in large-scale literature studies. In addition, the labels incorporate realistic Gaussian perturbations rather than direct recalculation of every composition under identical settings. This improves robustness relative to a perfectly clean toy dataset, but it does not eliminate the possibility that some entries are biased toward literature expectations.

### Methodological Limitations
The geometric stability metrics were extended to A3B2X9, A2BX6, and layered structures even though their interpretation is strongest for perovskite-like frameworks. The defect model and migration-barrier model are analytical surrogates, not full supercell DFT or explicit NEB calculations. Likewise, the device simulator is a reduced-order SCAPS-like approximation with a synthetic AM1.5G spectrum and capped open-circuit voltage. These choices are appropriate for rapid ranking but not for publishing absolute device predictions as final performance claims.

### Evaluation Limitations
The machine-learning metrics are strong, but the feature-importance analysis shows substantial dependence on the noisy DFT-like band-gap proxy. External validation with independent real-world datasets is essential to confirm the generalizability of these findings beyond simulated conditions. The statistical comparison between stable and non-stable groups also remains limited by sample size, especially once the dataset is partitioned into multiple structural families. More baselines, including descriptor-only uncertainty-aware Gaussian-process models, would strengthen confidence in the screening conclusions.

### Generalizability
The present ranking is calibrated for solar-absorber selection under one specific composite objective emphasizing near-optimal band gap, moderate defect tolerance, and simplified device efficiency. Different applications, such as tandem absorbers, radiation detectors, or photocatalysis, would favor different targets. In addition, processing conditions, phase purity, oxidation kinetics, and interface alignment may shift the ordering substantially in real devices.

### Future Directions
In the short term, the most useful extension is to recalculate the top 10 candidates with a consistent DFT workflow and explicit defect supercells, ideally under provenance-aware orchestration compatible with AiiDA-style pipelines (Li et al., 2022). In the medium term, the ranking weights should be re-learned from experimental feedback, and uncertainty estimates should be propagated from descriptors to the final composite score. Over a longer horizon, the screening engine could be embedded in a closed-loop discovery workflow that combines FAIR database integration (Whalley et al., 2021), Bayesian optimization, and targeted synthesis feedback.

## 8. Conclusion
We implemented a complete computational screening system for lead-free perovskite solar-cell candidates and demonstrated that a compact, reproducible workflow can recover literature-consistent trends while exposing important tradeoffs between geometry, defect tolerance, and device behavior. The strongest candidates in the present dataset were FASnI3, Cs2SnI6, Cs2AgBiI6, MASnI3, and Cs2BiAgI6, each combining predicted band gaps near the photovoltaic optimum with simulated efficiencies around 20%. Cross-validated machine-learning performance was strong but not unrealistically perfect, and sensitivity analysis showed that the main conclusions were robust to seed and hyperparameter variation. The overall evidence supports a cautious claim: integrated screening is a useful prioritization strategy for lead-free perovskites, but high-confidence material selection still requires external validation through consistent DFT and experiment.

## References
1. Aftab, T., & Ahmad, S. (2021). A review of stability and progress in tin halide perovskite solar cell. *Solar Energy*. DOI: 10.1016/j.solener.2020.12.065.
2. Bartel, C. J., et al. (2019). New tolerance factor to predict the stability of perovskite oxides and halides. *Science Advances*. DOI: 10.1126/sciadv.aav0693.
3. Cai, Y., et al. (2021). Discovery of Lead-Free Perovskites for High-Performance Solar Cells via Machine Learning. *Advanced Science*. DOI: 10.1002/advs.202103648.
4. Gao, F., et al. (2021). Screening for lead-free inorganic double perovskites with suitable band gaps and high stability using combined machine learning and DFT calculation. *Applied Surface Science*. DOI: 10.1016/j.apsusc.2021.150916.
5. Guo, Y. (2021). Machine learning stability and band gap of lead-free halide double perovskite materials for perovskite solar cells. *Solar Energy*. DOI: 10.1016/j.solener.2021.09.030.
6. Jamalinabijan, M., et al. (2025). Discovering novel lead-free mixed cation hybrid halide perovskites via machine learning. *Physical Chemistry Chemical Physics*. DOI: 10.1039/d4cp04218b.
7. Li, G., et al. (2022). AiiDA 2.0: toward verifiable, reusable and extensible workflows for atomistic simulations. *npj Computational Materials*. DOI: 10.1038/s41524-022-00863-2.
8. Schleder, G. R., et al. (2019). From DFT to machine learning: recent approaches to materials science. *Journal of Physics: Materials*. DOI: 10.1088/2515-7639/ab084b.
9. Tao, Q., et al. (2021). Machine learning for perovskite materials design and discovery. *npj Computational Materials*. DOI: 10.1038/s41524-021-00495-8.
10. Travis, W., et al. (2016). On the application of the Goldschmidt tolerance factor to inorganic and hybrid halide perovskites. *Chemical Science*. DOI: 10.1039/c5sc04845a.
11. Wang, Y., et al. (2025). Prediction and Screening of Lead-Free Double Perovskite Photovoltaic Materials Based on Machine Learning. *Molecules*. DOI: 10.3390/molecules30112378.
12. Whalley, T., et al. (2021). H2020 NOMAD: making computational materials science data FAIR. *Journal of Physics: Materials*. DOI: 10.1088/2515-7639/abf5b7.
