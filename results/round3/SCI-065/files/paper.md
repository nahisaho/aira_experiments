# DRAFT — NOT FOR DISTRIBUTION

# An Integrated Simulation Framework for Brain Organoid Perfusion Bioreactor Design, Media Programming, and Monitoring

## Abstract

Brain organoid bioprocess development is constrained by competing requirements for low-shear hydrodynamics, adequate oxygen delivery, temporal media control, scalable manufacturing, and interpretable monitoring. We developed a reproducible Python-based simulation framework that integrates axisymmetric low-Reynolds-number computational fluid dynamics, spherical reaction-diffusion transport, shear-responsive maturation modeling, temporal media optimization, techno-economic scalability analysis, and noisy biomarker trajectory simulation. The framework was designed as an early-stage engineering surrogate for screening reactor configurations before wet-lab prototyping. Under the baseline perfusion setting, the predicted Reynolds number was 0.637, the chamber pressure drop was 6.01e-04 Pa, and the peak wall shear stress was 1.03e-04 Pa, indicating a highly laminar operating regime. Oxygen transport emerged as the dominant biological bottleneck: viability decreased from 0.988 for 0.5 mm-radius organoids to 0.257 for 2.5 mm-radius organoids, while the necrotic core increased to 2.259 mm. The mechanosensing surrogate predicted an optimal shear stress near 0.112 Pa and a day-90 marker score of 0.521 ± 0.079. Media optimization suggested that approximately 3.0-day exchanges maximize a combined maturation-viability objective (39.80), and biomarker simulation yielded a multivariate maturation index of 0.746 ± 0.028 at day 90. Collectively, these results support perfusion-enabled, sensor-informed brain organoid culture while identifying oxygen-demand management as the most critical challenge for larger constructs.

## Introduction

Brain organoids have become an important experimental platform because they recreate elements of self-organization, regional identity, and extended neurodevelopment that are difficult to reproduce in two-dimensional cultures. Lancaster et al. (2013) established cerebral organoids as a tractable model of human brain development and linked long-term culture success to dynamic bioreactor support. Qian et al. (2016) extended this idea by showing that mini-bioreactors enable region-specific brain organoid differentiation, thereby illustrating that reactor design is not merely a convenience but an active determinant of developmental outcome. These foundational studies motivated a growing body of organoid engineering research focused on hydrodynamics, mass transport, and protocol standardization.

Later reviews and methods papers expanded the engineering framing. Hofer and Lütolf (2021) described organoids as engineered multicellular systems whose quality depends on controllable microenvironmental parameters. Zhao et al. (2022) summarized the broader methodological landscape and emphasized the need for standardized production and characterization pipelines. Kim et al. (2020) positioned human organoids as increasingly important translational models but noted persistent reproducibility challenges. Yang et al. (2023) updated that translational perspective, arguing that scaling and quality control remain central barriers. These sources collectively suggest that organoid bioprocessing requires explicit coupling of biological protocol design with reactor engineering and process analytics.

Dynamic culture systems have already demonstrated practical benefits. Suong et al. (2021) showed that vertical-mixing bioreactors can alter brain organoid morphology, indicating that flow conditions can influence tissue architecture. Saglam-Metiner et al. (2023) reported that spatiotemporal dynamic culture enhanced cellular diversity and neuronal function in human cerebral organoids. Seiler et al. (2022) found that automated microfluidic culture reduced glycolytic stress in cortical organoids, directly linking controlled perfusion-like environments to metabolic outcomes. Licata et al. (2023) synthesized these and related developments in a dedicated review of bioreactor technologies for organoid culture, while Velasco et al. (2020) described how microtechnology can integrate sensing, flow control, and miniaturized culture. Zheng et al. (2021) further argued that organoid-on-chip systems and numerical simulation can be combined to rationalize reactor and microenvironment design.

Even with these advances, the field still lacks compact design frameworks that combine hydrodynamics, transport, maturation, media scheduling, economics, and monitoring within one reproducible workflow. Many studies focus on one of these domains in isolation. For example, brain organoid biology studies may emphasize markers and morphology without reactor optimization, whereas engineering reviews may discuss transport and device design without a coherent maturation endpoint. Decarli et al. (2021), although centered on spheroids more broadly, highlighted transferable principles regarding high-throughput production and transport-limited behavior. Birey et al. (2017) demonstrated that advanced neural phenotypes and circuit-like assembly can emerge in 3D forebrain spheroids, underscoring the importance of long-term maturation endpoints. Schafer et al. (2023) extended the maturation discussion into neuroimmune phenotyping, showing that late-stage organoid analysis increasingly depends on multiplexed measurements.

The objective of the present work was therefore to build a computationally efficient design environment for brain organoid perfusion bioreactors. We hypothesized that an integrated framework would reveal a coherent operating regime characterized by laminar flow, modest but nonzero mechanostimulation, phase-specific media programming, and multimodal monitoring. We also expected the framework to identify oxygen transport as a dominant limiting factor for larger organoids. Rather than pursuing maximum physical fidelity, we prioritized interpretability, reproducibility, and extensibility so that the resulting codebase could serve as an engineering scaffold for later calibration with experimental data.

## Methods / Proposed Approach

### Study design and modeling philosophy

We implemented the framework in Python using six modules covering hydrodynamics, reaction-diffusion transport, mechanosensing-driven maturation, media optimization, scalability analysis, and biomarker monitoring. The goal was to create a fast surrogate environment suitable for early-stage design reasoning. Two modeling philosophies were considered for each subsystem. The first was high-fidelity numerical or data-driven modeling, such as full transient CFD, agent-based developmental simulation, or machine-learning regression from large training datasets. The second was reduced-order mechanistic modeling with explicit assumptions and moderate computational cost. We selected the second option because the current task required a complete, understandable, and executable system without access to matched experimental calibration datasets. Simpler analytical and semi-empirical models were therefore more appropriate than deep learning or industrial-grade multiphysics solvers.

### CFD module

The perfusion bioreactor was represented as a 2D axisymmetric cylindrical chamber with radius 5 mm and length 20 mm containing a centered spherical organoid of radius 1.5 mm. Flow rate was set within the requested microperfusion range, and the baseline value used for the reported results was 5 × 10^-9 m^3/s. Water-like properties were assumed (ρ=1000 kg/m^3, μ=10^-3 Pa·s). Because the resulting flow regime is strongly laminar, we used a Stokes-flow approximation in which inertia is neglected and the axial velocity profile in each axial slice is obtained by finite differences.

The characteristic flow state was summarized with the Reynolds number,

$$Re = \frac{\rho v D}{\mu},$$

where $v$ is mean axial velocity and $D$ is hydraulic diameter. The radial momentum balance was discretized in cylindrical coordinates. For each axial slice, the effective fluid annulus changed with the local cross-section of the organoid obstruction. The no-slip condition was enforced at the chamber wall and at the organoid boundary. After solving for a unit pressure-gradient profile, the solution was rescaled to satisfy the prescribed volumetric flow rate. Pressure drop was computed by axial integration of the inferred pressure gradient.

Wall shear stress was estimated from the outer radial velocity gradient,

$$\tau_w = \mu\left.\frac{dv}{dr}\right|_{r=R_w},$$

which is adequate for screening bioreactor operating windows in low-Reynolds-number conditions. A full finite-volume Navier–Stokes solver was considered but rejected for this stage because the increased implementation burden would not materially improve reactor ranking under creeping-flow assumptions. Likewise, a lattice-Boltzmann approach was not adopted because the design task focused on interpretable metrics rather than benchmark-grade CFD. The module used a 50 × 30 axial-radial grid and generated publication-style contour and profile figures.

### Reaction-diffusion transport module

Organoids were modeled as spheres with radii of 0.5, 1.0, 1.5, 2.0, and 2.5 mm. Oxygen transport was described with a steady-state spherical reaction-diffusion equation with Michaelis–Menten consumption,

$$D_{eff}\left(\frac{1}{r^2}\frac{d}{dr}\left(r^2\frac{dC}{dr}\right)\right)-\frac{Q_{max}C}{K_m+C}=0,$$

using $D_{eff}=2.0\times10^{-9}$ m$^2$/s, $Q_{max}=0.01$ mol/(m$^3$·s), $K_m=0.001$ mol/m$^3$, and surface concentration $C_s=0.2$ mol/m$^3$. Symmetry at the center and a fixed surface boundary concentration were imposed, and the resulting boundary-value problem was solved numerically. Glucose transport used the same structural form but with $D_{gluc}=6.7\times10^{-10}$ m$^2$/s and $Q_{gluc}=5\times10^{-4}$ mol/(m$^3$·s).

A Thiele modulus provided a compact reaction-versus-diffusion summary,

$$\phi = R\sqrt{\frac{Q_{max}}{D_{eff}C_s}},$$

and necrotic core radius was defined by $C < 0.01$ mol/m$^3$. Viability fraction was computed from the fraction of organoid volume above the critical oxygen threshold. This approach was preferred to purely empirical size–viability regressions because it provides physically interpretable scaling with radius and consumption demand. A higher-fidelity alternative would be a spatially heterogeneous multiphase model with variable cell density and internal cavities; that option was not selected because it would exceed the scope of an initial design framework. Glucose results were retained because they provide a useful negative control: they showed that not every nutrient becomes limiting simultaneously.

### Shear-responsive maturation module

Mechanosensitive maturation was modeled using a non-monotonic shear-response curve that includes low-shear activation and high-shear penalty:

$$M(\tau)=b + A\left(1-e^{-\tau/\tau_{opt}}\right)e^{-\left(\tau/\tau_{max}\right)^\beta},$$

where $b$ is baseline maturation, $A$ is activation amplitude, $\tau_{opt}$ is the low-shear activation scale, $\tau_{max}$ is the high-shear penalty scale, and $\beta$ governs penalty steepness. Time-dependent maturation was then modulated by a sigmoidal developmental term,

$$S(t)=\frac{1}{1+e^{-k(t-t_{50})}},$$

with $k=0.08$ day^-1 and $t_{50}=45$ days. Marker-specific parameter sets were assigned to TBR1, CTIP2, SOX2, MAP2, SYN1, and MBP to represent different developmental sensitivities. This surrogate method was selected over machine learning because no matched longitudinal training data were provided. A simpler monotonic shear-benefit model was also rejected because it cannot represent overstimulation penalties that are biologically plausible in perfused culture.

### Media optimization, scalability, and monitoring

Temporal media programming was divided into five phases: embryoid body formation, neural induction, patterning, maturation, and long-term maintenance. Phase-specific growth-factor scales, lactate accumulation, glucose consumption, and pH drift were simulated over 90 days. The optimization objective was the time integral of maturation multiplied by viability. Bounded continuous optimization identified factor scales and exchange interval that maximize this objective. Alternative strategies such as discrete rule-based scheduling or reinforcement learning were considered, but gradient-based continuous optimization was preferred because it is simpler to reproduce and sufficient for the low-dimensional decision space.

Scalability analysis compared batch, perfusion, and continuous modes from 1 to 1000 organoids per batch. Each mode included oxygen supply-to-demand ratios, media consumption, cost-per-organoid estimates, batch variability, and technology readiness levels. For monitoring, noisy trajectories were simulated for lineage markers, maturation markers, metabolic markers, and structural features. Detection modalities included RT-qPCR, immunohistochemistry, live imaging, secretome analysis, and impedance spectroscopy. A multivariate maturation index (MMI) was derived using PCA-like weights from TBR1, CTIP2, MAP2, SYN1, and MBP.

### Statistical analysis and sensitivity analysis

Because the study used synthetic mechanistic outputs, statistics were treated as descriptive. Five stochastic observation seeds and ±15% hyperparameter perturbations were used to assess robustness. Mean ± SD and 95% confidence intervals were reported for core outputs, and Welch t-tests with Cohen's $d$ were used for pairwise comparisons between selected conditions. These results are documented in `results/statistical-summary.md` and `results/sensitivity-analysis.md`. No multiple-comparison-adjusted claims are made because the comparisons are illustrative engineering checks rather than confirmatory biological inference.

## Results

### Hydrodynamic regime and shear distribution

The baseline perfusion design remained in a firmly laminar regime, with Reynolds number 0.637. Peak wall shear stress was 1.03e-04 Pa and mean wall shear stress was 6.09e-05 Pa, while the total chamber pressure drop was 6.01e-04 Pa. These values indicate that the simulated reactor offers gentle convection rather than aggressive mechanical conditioning.

![Figure 1: Axisymmetric velocity field in the perfusion chamber showing flow acceleration around the organoid obstruction](figures/cfd_velocity_field.png)

![Figure 2: Wall shear stress distribution along the chamber axis under baseline perfusion conditions](figures/cfd_shear_stress.png) In the sensitivity sweep, scaling the flow rate by ±15% led to nearly proportional changes in both Reynolds number and maximum shear stress, confirming that the CFD module behaves consistently under expected operational adjustments.

### Oxygen limitation dominates size-dependent viability

Oxygen profiles degraded sharply with increasing radius. At 0.5 mm radius, the model predicted a necrotic core radius of 0.114 mm and viability fraction of 0.988. At 1.0 mm radius, viability declined to 0.600; at 1.5 mm, to 0.421; at 2.0 mm, to 0.324; and at 2.5 mm, to 0.257, with a necrotic core of 2.259 mm.

![Figure 3: Oxygen concentration profiles across organoid radius for five organoid sizes (0.5–2.5 mm)](figures/o2_profiles.png)

![Figure 4: Necrotic core radius and viability fraction as a function of organoid radius](figures/necrotic_core.png) By contrast, minimum glucose concentration remained above 4.24 mol/m^3 even at the largest radius, suggesting that oxygen, not glucose, is the main bottleneck in the current configuration. The seed-based synthetic replicate summary in `results/statistical-summary.md` showed large effect size differences between small and large organoid viability, reflecting a strong modeled transport penalty.

### Maturation is maximized in an intermediate shear regime

The maturation model predicted a non-monotonic optimum rather than monotonic benefit with increasing shear. The highest day-90 score occurred at approximately 0.112 Pa, where the overall maturation metric reached 0.603.

![Figure 5: Maturation score heatmap as a function of shear stress (x-axis) and culture day (y-axis)](figures/maturation_heatmap.png)

![Figure 6: Marker-specific maturation curves (TBR1, CTIP2, SOX2, MAP2, SYN1, MBP) vs. shear stress at day 90](figures/marker_maturation.png) Marker-specific analysis showed a day-90 mean of 0.521 ± 0.079 across TBR1, CTIP2, SOX2, MAP2, SYN1, and MBP. TBR1 and CTIP2 peaked at somewhat lower shear than SYN1 and MBP, indicating that later-stage synaptic and myelination-associated markers tolerate or benefit from slightly stronger mechanical stimulation. In the descriptive comparison between 0.01 Pa and the model-derived optimum, the difference in day-90 maturation remained large in standardized units, reinforcing that an intermediate target is preferable to either near-static or excessively stressed conditions.

### Media scheduling yields a broad optimum near three-day exchanges

The temporal media optimizer identified an optimal exchange interval of 3.0 days with integrated maturation-times-viability score 39.80.

![Figure 7: Temporal media composition profile across five culture phases (days 0–90)](figures/media_temporal.png)

![Figure 8: Optimization landscape showing integrated maturation×viability score vs. exchange interval](figures/media_optimization.png) Under this schedule, pH remained within the requested 7.2–7.6 range, lactate accumulation was moderated, and the maturation trajectory increased while viability stayed in an acceptable window. Hyperparameter sweeps across 2.4–3.6 day exchange intervals showed a shallow optimum around 3 days rather than a brittle single-point optimum. This is operationally useful because manufacturing processes benefit from schedules that are robust to modest deviations. The results suggest that carefully timed medium replacement can control metabolic stress without requiring extreme factor loading.

### Scale-up favors perfusion for quality, but economics remain transitional

The scalability module indicated that perfusion and continuous culture reduce modeled batch-to-batch variability more effectively than static batch processes as organoid number increases.

![Figure 9: Scalability curves comparing batch, perfusion, and continuous culture modes across 1–1000 organoids](figures/scalability_curves.png)

![Figure 10: Cost analysis and break-even comparison between culture modes](figures/cost_analysis.png) However, under the chosen stylized cost coefficients, break-even for perfusion relative to batch and for continuous relative to perfusion occurred only near the upper end of the explored range, approximately 1000 organoids. These findings imply that engineering sophistication improves process uniformity before it necessarily improves short-term cost. Technology readiness levels also followed an expected hierarchy: batch at TRL 8, perfusion at TRL 6, and continuous processing at TRL 4. Thus, the model supports perfusion as an intermediate translational path: more controllable than batch culture and more implementable than fully continuous production.

### Multimodal monitoring captures maturation progression without requiring exclusively destructive assays

The biomarker module generated noisy but realistic marker trajectories over 90 days. The multivariate maturation index increased from 0.224 at day 30 to 0.558 at day 60 and 0.746 ± 0.028 at day 90.

![Figure 11: Biomarker trajectories over 90 days with realistic noise (TBR1, CTIP2, SOX2, MAP2, SYN1, MBP, LGR) and multivariate maturation index](figures/biomarker_trajectories.png)

![Figure 12: Proposed monitoring strategy timeline showing destructive and non-destructive assay cadence](figures/monitoring_timeline.png) PCA-like weights were broadly balanced across TBR1, CTIP2, MAP2, SYN1, and MBP, preventing the index from collapsing onto a single readout. The proposed monitoring timeline combines weekly secretome analysis, frequent impedance readouts, weekly live imaging, and less frequent destructive RT-qPCR or IHC checkpoints. This strategy reflects the maturation complexity highlighted by Birey et al. (2017) and Schafer et al. (2023), while also remaining compatible with the process-control emphasis discussed by Velasco et al. (2020). The MMI also provides a practical decision variable for comparing bioreactor runs with different maturation and viability trade-offs.

## Discussion

The integrated results support the idea that brain organoid bioprocess optimization should be treated as a coupled transport–development–monitoring problem. In the present framework, low baseline perfusion is hydrodynamically safe but insufficient by itself to generate the mechanostimulatory optimum suggested by the maturation module. This creates an important engineering tension. If flow is increased solely to chase the predicted maturation optimum, oxygen delivery may improve, but chamber shear could eventually approach levels that suppress some markers. The appropriate strategy is therefore not indiscriminate flow escalation but coordinated adjustment of chamber geometry, residence time, and perhaps local mixing architecture.

The central role of oxygen limitation is consistent with the bioreactor literature. Licata et al. (2023) reviewed mass-transfer improvements as a principal advantage of organoid bioreactors. Seiler et al. (2022) showed that controlled microfluidic culture reduces glycolytic stress, indirectly supporting the importance of transport regulation. Saglam-Metiner et al. (2023) demonstrated that dynamic culture improves maturation and diversity, which our framework helps rationalize through improved resource delivery and waste clearance. The results also agree with the broader organoid engineering argument of Hofer and Lütolf (2021): environmental control is not ancillary to organoid biology; it is constitutive of organoid quality.

The cost and readiness results suggest a staged translation pathway. Batch culture remains attractive because it is mature, familiar, and operationally simple. Perfusion offers a compelling compromise by decreasing modeled variability while keeping technology readiness within reach. Continuous manufacturing, although appealing conceptually, still appears less mature and more equipment-intensive. This conclusion aligns with the broader translational caution expressed by Kim et al. (2020), Zhao et al. (2022), and Yang et al. (2023), each of whom stressed that scalable organoid production requires standardization and verification rather than device complexity alone.

The monitoring design is another practical contribution of the framework. Many organoid programs still rely heavily on terminal endpoint assays, but long-culture systems are better served by layered monitoring that can detect trajectory shifts early. Velasco et al. (2020) and Zheng et al. (2021) pointed toward this integration of microtechnology and computation, and the present biomarker module translates that idea into a specific monitoring cadence. RT-qPCR and IHC remain important for lineage confirmation and spatial validation, but non-destructive secretome, imaging, and impedance measurements can support weekly or subweekly control decisions. ⚠️ The exact cost–sensitivity matrix used here is a modeled heuristic and should be recalibrated against facility-specific workflows.

The framework also highlights how different literature strands can be combined. Lancaster et al. (2013) and Qian et al. (2016) established the developmental and reactor foundations. Birey et al. (2017) and Schafer et al. (2023) illustrate how late-stage phenotypes depend on extended culture quality. Decarli et al. (2021) contributed transferable scale-up principles from spheroid systems. Together, these studies support the need for integrated design logic rather than isolated protocol optimization. The present work does not claim to be a validated digital twin, but it does provide a structured scaffold for building one. Zheng et al. (2021) made a similar case for combining numerical simulation with organoid engineering, and the present study operationalizes that recommendation in a modular form.

## Limitations and Future Work

### Data Limitations

Only synthetic data were used in this study. The model parameters were selected from literature-informed ranges and engineering judgment rather than from a matched experimental calibration dataset. Consequently, the quantitative outputs are best interpreted as plausible screening estimates. No patient-derived organoids, lot-specific extracellular matrix measurements, dissolved oxygen traces, or facility-specific cost data were available. Different stem cell lines, regionalization protocols, and matrix formulations could materially change transport properties, mechanosensitivity, and biomarker timing.

### Methodological Limitations

The hydrodynamic model is a reduced-order axisymmetric approximation with steady flow, Newtonian fluid properties, and no transient pump pulsatility. The transport model assumes homogeneous tissue, constant diffusivity, and spherical symmetry. The maturation model is phenomenological rather than fit to matched longitudinal molecular data. The media and scalability modules also use stylized surrogate equations that intentionally trade realism for transparency and computational speed. These limitations mean the framework is better suited for comparative screening than for final engineering certification.

### Evaluation Limitations

The evaluation primarily assessed internal logic, numerical stability, parameter sensitivity, and consistency with published qualitative trends. It did not compare predictions directly with measured shear fields, oxygen microelectrode profiles, real media chemistry, or actual manufacturing records. External validation with independent real-world datasets is essential to confirm the generalizability of these findings beyond simulated conditions. The descriptive statistics provided here include mean ± SD, 95% confidence intervals, and effect sizes, but those quantities arise from synthetic replicate perturbations rather than independent biological cohorts.

### Generalizability

The conclusions are most applicable to nonvascularized brain organoids cultured under low-flow perfusion conditions similar to those assumed here. They may not transfer directly to assembloids, organoids with engineered vasculature, or organ-on-chip systems with highly structured microchannels. Likewise, the cost model is illustrative and will vary with labor, automation, and local procurement. Domain shift should therefore be expected across laboratories and organoid platforms.

### Future Directions

In the short term, approximately six months, the most important next step is calibration against experimentally measured oxygen, lactate, and biomarker time courses. A parallel effort should benchmark CFD outputs against particle-image velocimetry or tracer-based flow measurements in benchtop reactor prototypes. In the medium to long term, one to two years, the framework could be expanded to multi-organoid occupancy, adaptive control of perfusion and exchange intervals, and sensor-driven parameter updates that transform the present surrogate into a true process digital twin. Further integration with electrophysiology, calcium imaging, and single-cell omics would also strengthen the maturation module.

## Conclusion

We present a complete, executable simulation system for brain organoid bioreactor design and optimization. The framework integrates fluid mechanics, mass transport, developmental response modeling, temporal media control, scale-up economics, and monitoring strategy in a single reproducible environment. The major conclusion is not that one fixed reactor configuration is universally optimal, but that a consistent design logic emerges when these domains are considered together: use laminar perfusion, prioritize oxygen management as organoid size increases, target intermediate rather than maximal shear, maintain phase-aware media programming with roughly three-day exchanges, and monitor maturation with a hybrid destructive and non-destructive assay stack. These claims are calibrated to a synthetic modeling setting, but they provide a practical starting point for experimental refinement.

## References

1. Birey, F., Andersen, J., Makinson, C. D., Islam, S., Wei, W., Huber, N., ... & Paşca, S. P. (2017). Assembly of functionally integrated human forebrain spheroids. *Nature*, 545, 54–59. https://doi.org/10.1038/nature22330
2. Decarli, M. C., Amaral, R., Santos, D. P., Tofani, L. B., Katayama, E. T., Rezende, R. A., ... & Moraes, A. M. (2021). Cell spheroids as a versatile research platform: formation mechanisms, high throughput production, characterization and applications. *Biofabrication*, 13(3), 032002. https://doi.org/10.1088/1758-5090/abe6f2
3. Hofer, M., & Lütolf, M. P. (2021). Engineering organoids. *Nature Reviews Materials*, 6, 402–420. https://doi.org/10.1038/s41578-021-00279-y
4. Kim, J., Koo, B.-K., & Knoblich, J. A. (2020). Human organoids: model systems for human biology and medicine. *Nature Reviews Molecular Cell Biology*, 21, 571–584. https://doi.org/10.1038/s41580-020-0259-3
5. Lancaster, M. A., Renner, M., Martin, C.-A., Wenzel, D., Bicknell, L. S., Hurles, M. E., ... & Knoblich, J. A. (2013). Cerebral organoids model human brain development and microcephaly. *Nature*, 501, 373–379. https://doi.org/10.1038/nature12517
6. Licata, J. P., Schwab, K. H., Har-el, Y., Gerstenhaber, J. A., & Lelkes, P. I. (2023). Bioreactor Technologies for Enhanced Organoid Culture. *International Journal of Molecular Sciences*, 24(14), 11427. https://doi.org/10.3390/ijms241411427
7. Qian, X., Nguyen, H. N., Song, M. M., Hadiono, C., Ogden, S. C., Hammack, C., ... & Ming, G.-L. (2016). Brain-Region-Specific Organoids Using Mini-bioreactors for Modeling ZIKV Exposure. *Cell*, 165(5), 1238–1254. https://doi.org/10.1016/j.cell.2016.04.032
8. Saglam-Metiner, P., Devamoglu, U., Filiz, Y., Akbari, S., Beceren, G., Bagca, B. G., ... & Yeşil-Çeliktaş, Ö. (2023). Spatio-temporal dynamics enhance cellular diversity, neuronal function and further maturation of human cerebral organoids. *Communications Biology*, 6, 117. https://doi.org/10.1038/s42003-023-04547-1
9. Schafer, S. T., Mansour, A. A., Schlachetzki, J. C. M., Pena, M., Ghassemzadeh, S. S., Mitchell, L. E., ... & Gage, F. H. (2023). An in vivo neuroimmune organoid model to study human microglia phenotypes. *Cell*, 186(10), 2111–2126. https://doi.org/10.1016/j.cell.2023.04.022
10. Seiler, S. T., Mantalas, G. L., Selberg, J., Cordero, S., Torres, S., Baudin, P. V., ... & Teodorescu, M. (2022). Modular automated microfluidic cell culture platform reduces glycolytic stress in cerebral cortex organoids. *Scientific Reports*, 12, 15108. https://doi.org/10.1038/s41598-022-20096-9
11. Suong, D. N. A., Imamura, K., Inoue, I., Kabai, R., Sakamoto, S., Okumura, T., ... & Inoue, H. (2021). Induction of inverted morphology in brain organoids by vertical-mixing bioreactors. *Communications Biology*, 4, 1213. https://doi.org/10.1038/s42003-021-02719-5
12. Velasco, V., Shariati, S. A., & Esfandyarpour, R. (2020). Microtechnology-based methods for organoid models. *Microsystems & Nanoengineering*, 6, 76. https://doi.org/10.1038/s41378-020-00185-3
13. Yang, S., Hu, H.-J., Kung, H.-C., Zou, R.-Q., Dai, Y.-S., Hu, Y.-F., ... & Li, F.-Y. (2023). Organoids: The current status and biomedical applications. *MedComm*, 4(3), e274. https://doi.org/10.1002/mco2.274
14. Zhao, Z., Chen, X., Dowbaj, A. M., Sljukic, A., Bratlie, K. M., Lin, L., ... & Yu, H. (2022). Organoids. *Nature Reviews Methods Primers*, 2, 94. https://doi.org/10.1038/s43586-022-00174-y
15. Zheng, F., Xiao, Y., Liu, H., Fan, Y., & Dao, M. (2021). Patient-Specific Organoid and Organ-on-a-Chip: 3D Cell-Culture Meets 3D Printing and Numerical Simulation. *Advanced Biology*, 5(6), 2000024. https://doi.org/10.1002/adbi.202000024
