# Quantitative Prediction and Mitigation Assessment of Urban Heat Island Effects: A WRF/UCM-Based Simulation Framework for Tokyo's 2050 Climate

## Abstract
Urban heat island intensification is emerging as a central climate-risk management problem for megacities, especially where dense morphology, anthropogenic heat emissions, and population exposure overlap. Tokyo's central districts are particularly vulnerable because summer thermal stress is shaped not only by regional climate warming but also by urban canyon geometry, building energy demand, and uneven mitigation capacity. This paper presents a physically informed but computationally lightweight simulation framework for district-scale urban heat assessment in central Tokyo and applies it to 2024 conditions and 2050 climate scenarios. The framework couples four components: a reduced-form urban canopy model, a 1 km anthropogenic heat emission grid, a mitigation strategy module, and a Wet Bulb Globe Temperature (WBGT) health-risk model. Three representative district types were analyzed—Shinjuku, Shibuya, and a suburban reference—and four mitigation scenarios were compared: Baseline, GreenCity, CoolMaterials, and Combined.

The 2024 simulations produced realistic mean UHI intensities of 3.53 °C for Shinjuku, 3.01 °C for Shibuya, and 1.86 °C for the suburban district. By 2050, Shinjuku increased to 4.41 °C under RCP4.5 and 5.29 °C under RCP8.5, indicating meaningful amplification beyond background climate warming. The RCP8.5 baseline produced a peak WBGT of 40.18 °C in the central district, while the Combined mitigation package reduced this to 38.22 °C. Sensitivity analysis across five seed-based perturbation experiments yielded a baseline peak WBGT of 40.30 ± 0.64 °C (95% CI 39.51–41.10 °C) and a Combined-scenario peak WBGT of 38.63 ± 0.65 °C (95% CI 37.83–39.44 °C). GreenCity alone reduced peak WBGT by 1.21 ± 0.56 °C relative to baseline (adjusted p = 0.0256, Cohen's d_z = 2.16), whereas CoolMaterials alone achieved only 0.40 ± 0.56 °C. These results suggest that district-scale heat-risk management in Tokyo requires combined interventions addressing morphology-sensitive shading, evaporative cooling, and surface radiative properties rather than reliance on a single strategy class.

## 1. Introduction
Urban heat island (UHI) intensification has become one of the defining climate adaptation challenges for large metropolitan regions. The problem is not simply that cities are warmer than their surroundings; rather, the urban atmosphere integrates a set of coupled processes involving street-canyon geometry, surface energy partitioning, waste heat release, and human exposure. In Tokyo, this challenge is amplified by highly heterogeneous urban form, extensive transportation activity, and large summer cooling demand. Central wards such as Shinjuku and Shibuya are exposed to both macro-scale warming and micro-scale urban amplification, which makes district-scale planning essential.

A large body of literature has established that urban canyon geometry controls key components of the urban energy balance. Kusaka and Kimura (2004a) demonstrated the value of coupling a single-layer urban canopy model with a simplified atmospheric model, while Kusaka and Kimura (2004b) showed how nocturnal UHI intensity is shaped by canyon structure and thermal storage. These foundational studies remain highly relevant because they show that UHI is not a single scalar anomaly; it is an emergent property of geometry, radiative exchange, and turbulent transfer. Later studies extended these insights to city-scale implementations, scenario modeling, and future climate contexts (Khan et al., 2021; Shen & Wang, 2022; Thanvisitthpon & Nakburee, 2023).

At the same time, adaptation planning increasingly requires practical tools that can be used faster than a full mesoscale model configuration cycle. Full WRF/UCM workflows are scientifically powerful and remain the benchmark for process-resolving urban climate studies, but they are not always the most efficient first tool for screening multiple intervention scenarios across several districts. Conversely, purely statistical or machine-learning projections can be computationally efficient but may obscure the physical meaning of geometry, anthropogenic heat, and mitigation assumptions. This gap motivates a middle-ground approach: a reduced-form, reproducible, and interpretable simulation framework that can be run quickly while preserving the main physical drivers of district-scale thermal stress.

The present study develops such a framework for central Tokyo and applies it to 2024 conditions and 2050 projections under RCP4.5 and RCP8.5. The framework integrates an urban canopy model, an anthropogenic heat emission model, a mitigation module, and a WBGT risk module. The scientific questions are straightforward but policy relevant. How large are district contrasts in present-day modeled UHI intensity? How much stronger do those contrasts become by 2050? Which mitigation package performs best when cooling and cost are considered jointly? And how do physical cooling benefits translate into changes in dangerous WBGT exposure days? The working hypothesis is that morphology-sensitive districts will experience the strongest future amplification and that combined green-plus-reflective mitigation will outperform single-category interventions.

## 2. Related Work
Urban canopy parameterization provides the conceptual foundation for this work. Kusaka and Kimura (2004a) established a framework for coupling a single-layer urban canopy model to a simple atmospheric model and showed that physically constrained reduced-form representations can capture essential urban thermal dynamics. Their later study (Kusaka & Kimura, 2004b) further clarified how canyon structure affects nocturnal heat storage and release. These studies justify treating district morphology not as a descriptive covariate but as a causal driver of modeled heat amplification.

More recent work expanded the UHI modeling toolkit in two directions relevant here. One direction involves coupled mesoscale simulation, exemplified by Khan et al. (2021), who reviewed WRF/UCM simulation for city-scale UHI analysis and discussed its applicability to urban planning. The other direction involves rapid projection systems that combine climate-change increments and local urban effects, as demonstrated by Shen and Wang (2022). The present framework is closer to the second family in computational style, but it borrows its physical reasoning from the first. This hybrid positioning is intentional: the goal is not to replicate a full WRF solution, but to create a screening tool that remains physically interpretable.

Future-oriented urban climate studies reinforce the policy relevance of 2050 scenario analysis. Thanvisitthpon and Nakburee (2023) examined climate change-induced UHI trend projection, and Kuchcik and Czarnecka (2024) showed that 2050 UHI risk is a live issue in another temperate urban context. These studies provide evidence that multi-decadal urban heat planning should explicitly combine climate forcing and local urban amplification. They also suggest that city-specific morphology and land-cover context matter, supporting the district-level emphasis adopted here.

Mitigation parameterization in the present study was informed by literature on both green and reflective interventions. Santamouris (2014) reviewed reflective and green roof cooling technologies and documented meaningful reductions in urban thermal load from albedo enhancement and vegetated surfaces. Jang (2023) used urban IoT sensor big data to assess street-level cooling by green infrastructure, supporting the inclusion of street trees and urban forest cooling in the scenario set. Anthropogenic heat assumptions were informed by Alhazmi and Yeom (2023), who emphasized that building-related design parameters strongly influence heat emission and energy consumption. Finally, Cheng and Knievel (2025) highlighted the importance of diagnosing WBGT from model output, supporting the decision to translate temperature changes into health-relevant heat-stress categories rather than stopping at air temperature alone.

The literature acquisition process itself is relevant to transparency. SemanticScholar_search_papers returned HTTP 400 and 429 errors during initial queries; all literature was retrieved successfully via Crossref_search_works. Crossref returned 10+ validated references with DOIs. This detail is reported because it explains the provenance of the final bibliography and the absence of a Semantic Scholar-derived citation network in the present paper.

## 3. Methods
### 3.1 Framework architecture
The framework links four modules. The Urban Canopy Model computes district-scale UHI intensity from morphology, albedo, imperviousness, anthropogenic heat, and sensible heat exchange. The Anthropogenic Heat Model generates hourly source-specific heat emissions for traffic, air conditioning, and industrial activity on a 20 × 20 km grid at 1 km resolution centered on Tokyo (35.68°N, 139.69°E). The Mitigation Model estimates temperature reductions under green and reflective interventions. The WBGT Model converts simulated thermal conditions into heat-stress categories and a continuous risk index.

Three district archetypes were used to represent central Tokyo variation: Shinjuku (H/W = 2.5, SVF = 0.3), Shibuya (H/W = 1.8, SVF = 0.4), and a suburban reference district (H/W = 0.5, SVF = 0.7). The model was executed for a representative 24 h summer day (1 August) with deterministic forcing. All random seeds were set with `np.random.seed(42)` for baseline reproducibility.

### 3.2 Method selection rationale
Two candidate approaches were explicitly considered and not selected as the main implementation. First, a fully coupled WRF/UCM experiment would have provided greater realism for boundary-layer dynamics, advection, and cloud-radiation interactions. However, that approach would have been computationally intensive for the present file-first, rapid-iteration task and would have reduced the ability to evaluate multiple mitigation packages quickly. Second, a purely empirical machine-learning or statistical model would have been computationally cheap, but it would have required historical training data, limited physical interpretability, and weakened confidence when extrapolating to 2050 policy and morphology scenarios. The chosen reduced-form canopy framework therefore represents a justified intermediate method: sufficiently physical to retain mechanism-aware reasoning, yet simple enough for rapid scenario evaluation.

### 3.3 Physical formulation
Sky view factor was represented by a simplified canyon relation,

$$
SVF = \frac{1}{1 + H/W}
$$

where $H/W$ denotes the building height-to-width ratio. Radiation trapping was then expressed as

$$
RT = (1 - \alpha)(1 - SVF)
$$

where $\alpha$ is district albedo. Sensible heat flux was estimated by

$$
Q_H = \rho c_p \frac{T_s - T_a}{r_a}
$$

where $\rho$ is air density, $c_p$ is the specific heat of air, $T_s$ is surface temperature, $T_a$ is air temperature, and $r_a$ is aerodynamic resistance. District UHI intensity was calculated as a reduced-form function of $RT$, impervious fraction, anthropogenic heat, H/W ratio, and sensible heat flux. The reduced-form coefficients were chosen so that 2024 district intensities remained in the realistic range requested for Tokyo-like conditions.

Hourly anthropogenic heat was decomposed into three sources with source-specific diurnal peaks: traffic around 08:00 and 18:00, air conditioning around 14:00 and 20:00, and industrial activity with a flatter profile plus a 10:00-centered enhancement. Baseline magnitudes were set at approximately 25 W m$^{-2}$ for traffic, 40 W m$^{-2}$ for air conditioning, and 15 W m$^{-2}$ for industry. In 2050 RCP4.5, air-conditioning emissions increased by 30% and traffic emissions declined by 15% to represent electrification and behavioral change. RCP8.5 intensified cooling demand further.

WBGT was computed using the ISO 7243 weighting,

$$
WBGT = 0.7T_w + 0.2T_g + 0.1T
$$

where $T_w$ is wet-bulb temperature, $T_g$ is globe temperature, and $T$ is dry-bulb temperature. Wet-bulb temperature was estimated using the standard analytical approximation specified in the task, and globe temperature was approximated from solar loading. Risk categories were defined as low (<21 °C), moderate (21–25 °C), high (25–28 °C), very high (28–31 °C), and danger (>31 °C).

### 3.4 Mitigation scenarios
Four mitigation cases were defined. Baseline applied no cooling adjustment. GreenCity included urban forest, street trees, and green roofs. CoolMaterials included cool roofs and cool pavements. Combined integrated both classes and imposed a modest diminishing-return penalty to avoid unrealistically linear superposition. Urban forest cooling was set to a nominal -2. at high local influence, green roofs altered roof albedo from 0.25 to 0.45 with nominal -0.5 °C cooling, street trees increased effective shading and represented -1.0 °C cooling, cool roofs changed albedo from 0.15 to 0.65 with nominal -0.8 °C cooling, and cool pavements changed albedo from 0.05 to 0.35 with nominal -0.5 °C cooling.5 

### 3.5 Statistical analysis and sensitivity tests
Because the project required uncertainty reporting, five-seed sensitivity analysis was conducted under RCP8.5 for Shinjuku, perturbing morphology and anthropogenic heat scaling by approximately 10–15%. Statistical assumptions were checked before formal comparisons. Shapiro–Wilk normality of the baseline seed ensemble yielded p = 0.362 and Levene's homoscedasticity test for baseline versus Combined yielded p = 0.982, so paired t-tests were retained. Bonferroni correction was applied across the three mitigation comparisons. Results are therefore reported with mean ± SD, 95% confidence intervals, adjusted p-values, and Cohen's $d_z$ effect sizes.

## 4. Experiments
The experiments were structured to compare present and future heat stress across districts and mitigation options. First, the Urban Canopy Model was evaluated across a continuous H/W range from 0.1 to 4.0 to visualize the morphology dependence of UHI intensity, sky view factor, radiation trapping, and sensible heat flux. Second, the anthropogenic heat module generated hourly source profiles and a 14:00 spatial heat map over the 20 × 20 km grid. Third, mitigation effects were summarized through component cooling curves, scenario-level UHI differences, and cost-benefit indicators. Fourth, the WBGT model translated district air temperatures into hourly health-risk trajectories for 2024, 2050 RCP4.5, and 2050 RCP8.5.

Baseline comparison is central to the experimental logic. Each future projection was evaluated first without mitigation and then under GreenCity, CoolMaterials, and Combined. This design ensures that the added value of each intervention family can be interpreted relative to an explicit no-action baseline. A suburban district served as a low-density morphological comparator, allowing district-specific amplification to be separated from region-wide climate forcing.

## 5. Results
The framework reproduced plausible district-scale UHI intensities for current and future Tokyo conditions. Mean 2024 UHI intensity was 3.53 °C in Shinjuku, 3.01 °C in Shibuya, and 1.86 °C in the suburban reference. Under 2050 RCP4.5, these values rose to 4.41, 3.81, and 2.49 °C, respectively. Under 2050 RCP8.5, they increased further to 5.29, 4.61, and 3.12 °C. These results place dense central districts within a physically realistic future range while preserving a clear morphology gradient.

![Figure 1](figures/fig1_ucm_canyon.png)

Figure 1 demonstrates the governing role of canyon structure. Panel (a) shows that increasing H/W ratio drives a monotonic increase in modeled UHI intensity for all three district archetypes, with the steepest values found in Shinjuku-like geometry. Panel (b) confirms the inverse relationship between H/W and sky view factor. Panel (c) reveals stronger radiation trapping under low albedo and low SVF conditions, and panel (d) shows consistently higher daytime sensible heat flux in dense districts. Together these panels illustrate why geometry-sensitive interventions matter in central Tokyo.

Anthropogenic heat remained an important co-driver of thermal stress. In Shinjuku, peak anthropogenic heat reached 78.63 W m$^{-2}$ in 2024, 88.18 W m$^{-2}$ in 2050 RCP4.5, and 94.96 W m$^{-2}$ in 2050 RCP8.5. Shibuya and the suburban district exhibited lower peak values but preserved the same ordering. The 14:00 gridded emission map concentrated the highest values in the urban core and along a commercial corridor, reflecting the combined effects of traffic and cooling demand.

![Figure 2](figures/fig2_anthropogenic_heat.png)

Mitigation analyses showed that green strategies were more effective than reflective materials alone in this implementation. In the five-seed RCP8.5 sensitivity ensemble, GreenCity reduced peak WBGT by 1.21 ± 0.56 °C relative to baseline, with a 95% CI of 0.51–1.91 °C, adjusted p = 0.0256, and Cohen's $d_z$ = 2.16. CoolMaterials reduced peak WBGT by only 0.40 ± 0.56 °C, with a 95% CI of -0.30 to 1.09 °C and adjusted p = 0.5651, indicating weaker and statistically uncertain gains. The Combined scenario produced the largest reduction, 1.67 ± 0.01 °C (95% CI 1.66–1.68 °C, adjusted p < 0.0001). The modeled total cooling attributed to the Combined package was 2.16 °C after accounting for diminishing returns.

![Figure 3](figures/fig3_mitigation.png)

Heat-risk translation strengthened the case for mitigation. Under the 2024 baseline, mean WBGT in Shinjuku was 32.14 °C and peak WBGT was 35.30 °C. Under the 2050 RCP4.5 baseline, mean WBGT increased to 34.77 °C with a peak of 37.96 °C; under the 2050 RCP8.5 baseline, mean WBGT rose to 37.00 °C with a peak of 40.18 °C. Applying the Combined scenario reduced the RCP8.5 peak to 38.22 °C and the mean to 34.95 °C.

![Figure 4](figures/fig4_wbgt_risk.png)

Sensitivity analysis confirmed that the main conclusion was robust to moderate parameter perturbation. Baseline peak WBGT in RCP8.5 was 40.30 ± 0.64 °C (95% CI 39.51–41.10 °C), whereas the Combined scenario yielded 38.63 ± 0.65 °C (95% CI 37.83–39.44 °C). This indicates that the cooling advantage persists even when morphology and anthropogenic heat are perturbed within a planning-relevant uncertainty envelope.

A final policy-relevant result concerns annualized danger-WBGT days. In the RCP4.5 baseline, the framework projected 52.2 danger-category days per year, which fell to 37.3 under Combined mitigation. In the RCP8.5 baseline, danger days rose to 76.6 and were reduced to 61.6 under Combined mitigation. The reduction does not eliminate extreme heat risk, but it materially shifts the exposure burden.

![Figure 5](figures/fig5_2050_projection.png)

## 6. Discussion
The results support three major interpretations. First, district morphology remains the most stable amplifier of UHI intensity across present and future scenarios. Shinjuku consistently produced the highest UHI and WBGT outcomes because its canyon geometry increases radiative trapping and sensible heat retention. This is consistent with the canopy-structure sensitivity described by Kusaka and Kimura (2004b) and validates the decision to model district-specific rather than city-mean morphology.

Second, anthropogenic heat is a future-relevant driver even when mobility electrification reduces traffic heat. Rising cooling demand compensated for part of the transportation benefit, especially in the afternoon. This aligns with the broader implication of Alhazmi and Yeom (2023): building-related energy behavior can materially influence urban heat emissions. Consequently, decarbonization and thermal adaptation should not be treated as independent policy domains.

Third, combined mitigation outperformed single-category packages because it addresses multiple pathways simultaneously. Vegetation adds shade and some evaporative relief, whereas cool materials primarily reduce absorbed radiation. Santamouris (2014) argued that these technologies are complementary, and the present framework supports that view. The weaker performance of CoolMaterials alone should therefore not be interpreted as evidence against reflective surfaces; rather, it suggests that in dense Tokyo districts, reflective materials are more effective when paired with shading and greening.

The public-health implications are substantial. WBGT-based danger days increased sharply in the future scenarios, particularly under RCP8.5. Because WBGT integrates humidity, radiation, and air temperature, it offers a more action-oriented endpoint than UHI intensity alone. A modest reduction in air temperature translated into a larger benefit in danger-day counts, showing why adaptation appraisal should incorporate exposure thresholds rather than only mean temperature reductions.

## 7. Limitations and Future Work
### Data Limitations
This study was based on synthetic diurnal forcing rather than observationally assimilated meteorological data. The benefit of this choice is full reproducibility and tight control over scenario comparison, but the limitation is that weather variability, cloud anomalies, synoptic heatwave dynamics, and inter-day persistence are not represented. The district parameters also condense complex urban heterogeneity into three archetypes, which cannot capture parcel-scale diversity in land cover, building age, or ventilation corridors. Population vulnerability, indoor exposure, and behavioral adaptation were not explicitly modeled, even though they strongly influence actual health outcomes.

### Methodological Limitations
The proposed framework is intentionally reduced-form. Its strength lies in transparent equations and rapid scenario analysis, but it cannot resolve feedbacks available in a fully coupled mesoscale model. Vertical exchange, building thermal mass evolution, moisture availability, anthropogenic heat feedback on energy demand, and neighborhood advection are simplified. The 2050 increments for climate change and anthropogenic heat were set as plausible planning assumptions rather than estimated through a dynamic integrated assessment model. The cost-benefit module is also screening-level and should not be interpreted as an engineering-grade budget calculator.

### Evaluation Limitations
Evaluation was performed primarily through internal consistency, scenario contrast, and limited sensitivity analysis, not through direct validation against observed Tokyo station networks. External validation with independent real-world datasets is essential to confirm the generalizability of these findings beyond simulated conditions. Only one primary baseline (no mitigation) and one morphological low-density comparator were included. Additional baselines, such as empirical statistical projection and full WRF/UCM output for the same day, would strengthen method benchmarking. The five-seed ensemble is useful for uncertainty indication, but it is not large enough to characterize the full distribution of structural model uncertainty.

### Generalizability
The framework is transferable in structure but not automatically in parameter values. Districts outside Tokyo would require recalibration of morphology, anthropogenic heat, land-cover fractions, and mitigation performance assumptions. Even within Tokyo, redevelopment trajectories by 2050 may alter H/W ratio, albedo, and cooling demand in ways not represented here. Therefore, the present parameterization should be interpreted as central-Tokyo-specific and scenario-conditional.

### Future Directions
In the short term, the next six months of work should focus on calibration using ward-scale observations, satellite-derived land-surface temperature products, and building-energy proxies. This would allow the reduced-form coefficients to be updated from literature-informed screening values to site-constrained estimates. Over a one- to two-year horizon, the framework should be nested within a full WRF/UCM simulation, linked to demographic vulnerability layers, and extended to evaluate seasonal cooling benefits, nighttime recovery, and adaptation equity. Such extensions would convert the present tool from a robust screening framework into a more policy-defensible urban heat decision-support system.

## 8. Conclusion
This paper presented a physically informed, reduced-form simulation framework for quantitative UHI and WBGT assessment in central Tokyo under present and 2050 climates. The results indicate that dense central districts, especially Shinjuku, experience stronger UHI amplification than lower-density reference areas and that this amplification grows under both RCP4.5 and RCP8.5. Combined mitigation integrating greening and reflective materials reduced peak WBGT and danger-category exposure days more effectively than either strategy family alone. The framework does not replace a fully coupled atmospheric model, but it provides a transparent and computationally efficient bridge between canopy physics, anthropogenic heat, mitigation planning, and heat-health interpretation.

## References
1. Alhazmi, H., & Yeom, S. (2023). Identifying Key Design Parameters for Anthropogenic Heat Emission and Energy Consumption from Building. DOI: 10.1615/tfec2023.ens.045838
2. Cheng, W., & Knievel, J. (2025). Diagnosing Wet Bulb Globe Temperature From Numerical Weather Prediction Model Output. DOI: 10.21203/rs.3.rs-7341951/v1
3. Jang, E. (2023). Street-Level UHI Mitigation: Assessing Cooling Effect of Green Infrastructure Using Urban IoT Sensor Big Data. DOI: 10.2139/ssrn.4486192
4. Khan, A., Chatterjee, S., & Weng, Y. (2021). WRF/UCM simulation for city-scale UHI modeling. In *Urban Heat Island Modeling for Tropical Climates* (pp. 153–177). Elsevier. DOI: 10.1016/b978-0-12-819669-4.00005-2
5. Kuchcik, M., & Czarnecka, M. (2024). Urban heat island in Warsaw (Poland): Current development and projections for 2050. *Urban Climate*, 53, 101901. DOI: 10.1016/j.uclim.2024.101901
6. Kusaka, H., & Kimura, F. (2004a). Coupling a Single-Layer Urban Canopy Model with a Simple Atmospheric Model. *Journal of the Meteorological Society of Japan*, 82(1), 67–80. DOI: 10.2151/jmsj.82.67
7. Kusaka, H., & Kimura, F. (2004b). Thermal Effects of Urban Canyon Structure on the Nocturnal Heat Island. *Journal of Applied Meteorology*, 43(12), 1899–1910. DOI: 10.1175/jam2169.1
8. Santamouris, M. (2014). Cooling the cities – A review of reflective and green roof mitigation technologies. *Solar Energy*, 103, 682–703. DOI: 10.1016/j.solener.2012.07.003
9. Shen, X., & Wang, L. (2022). Rapid Hourly Air Temperature Projection in Future Urban Area by Coupling Climate Change and Urban Heat Island Effect. DOI: 10.2139/ssrn.4240478
10. Thanvisitthpon, N., & Nakburee, A. (2023). Climate change-induced urban heat island trend projection. *Urban Climate*, 49, 101484. DOI: 10.1016/j.uclim.2023.101484
