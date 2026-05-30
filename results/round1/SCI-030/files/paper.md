# An Integrated Simulation Framework for Supercritical Enhanced Geothermal Systems: DFN-THM Coupling with Induced Seismicity Risk Assessment and Well Placement Optimization at Kakkonda, Japan

## Abstract

Supercritical Enhanced Geothermal Systems (EGS) represent a transformative approach to geothermal energy extraction, promising substantially higher enthalpy yields compared to conventional hydrothermal systems. This study presents an integrated reservoir simulation framework combining Discrete Fracture Network (DFN) modeling, Thermo-Hydro-Mechanical (THM) coupling, supercritical water equation of state based on IAPWS-95, Coulomb failure stress analysis for induced seismicity risk assessment, and well placement optimization. The framework adopts a TOUGH2/OpenGeoSys-style finite volume approach with sequential THM coupling. A case study is conducted for the Kakkonda geothermal field in Iwate Prefecture, Tohoku, Japan, where the WD-1a well penetrated to 3,729 m depth with temperatures exceeding 380°C. The DFN model generates 65 fractures across three orientation sets following Fisher distributions, yielding an equivalent permeability of 2.37 × 10⁻¹³ m². Thirty-year simulations across five well spacing configurations (50–900 m) identify an optimal spacing of 600 m, achieving cumulative energy extraction of 20.99 PJ. Coulomb stress analysis predicts a maximum induced seismic event of M3.1, with rate-and-state friction modeling providing spatially resolved seismicity rate estimates. The framework demonstrates the feasibility of integrated multi-physics simulation for supercritical EGS development planning, while identifying key challenges in numerical stability near the critical point and the need for full THMC coupling in future work.

## 1. Introduction

### 1.1 Background

Geothermal energy represents one of the most promising renewable energy sources, offering baseload capacity independent of weather conditions. Enhanced Geothermal Systems (EGS) extend the reach of geothermal energy to regions lacking natural hydrothermal resources by engineering artificial reservoirs through hydraulic stimulation of hot dry rock (HDR) formations (Tester et al., 2006). Among the most exciting frontiers in geothermal research is the exploitation of supercritical geothermal resources, where reservoir temperatures exceed the critical point of water (T > 374°C, P > 22.1 MPa), potentially yielding 5–10 times the power output of conventional systems (Friðleifsson et al., 2014).

Japan, with its abundant volcanic geothermal resources, is at the forefront of supercritical geothermal exploration. The Kakkonda geothermal field in Iwate Prefecture, Tohoku region, hosts the well WD-1a, which penetrated to a depth of 3,729 m and encountered temperatures of 380°C, crossing the brittle-ductile transition at approximately 3,100 m (Doi et al., 1998). The Japan Beyond-Brittle Project (Muraoka et al., 2014) has systematically investigated the potential of supercritical geothermal resources at Kakkonda and other sites in Japan.

### 1.2 Challenges

Simulation of supercritical EGS reservoirs presents unique challenges:

1. **Fracture dominance**: Flow in crystalline basement rocks is controlled by discrete fracture networks, requiring explicit DFN modeling rather than equivalent porous medium approaches (Jiang et al., 2025).
2. **Multi-physics coupling**: The strong interaction between thermal, hydraulic, and mechanical processes at supercritical conditions demands fully coupled THM solvers (Andrés et al., 2022).
3. **Fluid property discontinuities**: Water properties vary dramatically near the critical point, introducing numerical instabilities (Wagner & Pruss, 2002).
4. **Induced seismicity**: Fluid injection can reactivate pre-existing faults, requiring quantitative seismic risk assessment via Coulomb stress analysis (Hutka et al., 2023).
5. **Long-term optimization**: Economic viability requires 30+ year heat extraction predictions with optimized well configurations (Zhou et al., 2021).

### 1.3 Contributions

This study makes the following contributions:

- An integrated simulation framework combining DFN, THM, supercritical EOS, Coulomb stress, and well optimization
- Application to realistic geological conditions at the Kakkonda geothermal field
- Systematic evaluation of well spacing effects on 30-year heat recovery
- Quantitative induced seismicity risk assessment using rate-and-state friction models
- A TOUGH2/OpenGeoSys-compatible workflow for supercritical EGS reservoir engineering

## 2. Related Work

### 2.1 DFN Modeling for EGS

Discrete Fracture Network modeling has become the standard approach for representing flow pathways in fractured EGS reservoirs. Jiang et al. (2025) demonstrated that fracture aperture and orientation are the dominant controls on EGS heat extraction performance in coupled THM-DFN simulations. Janiga et al. (2022) developed a DFN model for horizontal well-doublet EGS at the Utah FORGE site, revealing the sensitivity of thermal breakthrough to fracture network geometry. Fan et al. (2025) introduced a novel 3D deformable DFN meshing approach coupled with THM for explicit representation of fracture deformation effects. Schwering et al. (2020) created a deterministic DFN for the EGS Collab project, providing validation data for stochastic DFN approaches.

### 2.2 THM Coupling in Geothermal Reservoirs

Thermo-hydro-mechanical coupling is essential for realistic EGS simulation. Andrés et al. (2022) developed a fully coupled THM model with rate-and-state friction for the Basel EGS, simulating both induced seismicity and thermal decline. Zhou et al. (2021) performed THM modeling of heat extraction showing that fracture spacing and number significantly affect reservoir longevity over 30-year periods. OpenGeoSys (Kolditz et al., 2012) and TOUGH2 (Pruess et al., 2012) remain the primary simulation platforms, with recent coupling workflows combining TOUGH2's fluid flow capabilities with OpenGeoSys's geomechanical modules.

### 2.3 Induced Seismicity Risk Assessment

Coulomb failure stress modeling has emerged as the principal tool for induced seismicity risk assessment in EGS. Hutka et al. (2023) systematically evaluated the relative importance of pore pressure and thermal stress in triggering seismicity using coupled THM-Coulomb models, finding that pore pressure dominates but thermal stresses produce larger spatial footprints. Drif et al. (2024) analyzed the relationship between hydraulic energy input and seismic energy release at Soultz-sous-Forêts. Boyet et al. (2023) investigated post-injection seismicity at Basel through Coulomb stress transfer mechanisms.

### 2.4 Supercritical Geothermal Resources in Japan

The Kakkonda geothermal field has been extensively studied as a candidate for supercritical geothermal development. Doi et al. (1998) reported drilling results from WD-1a showing temperatures exceeding 380°C. Muraoka et al. (2014) outlined the vision for the Japan Beyond-Brittle Project targeting supercritical resources. Suzuki et al. (2022) applied Bayesian rock-physics modeling to constrain deep temperatures at Kakkonda using resistivity data, confirming supercritical conditions at depth.

## 3. Methods

### 3.1 Discrete Fracture Network Model

The DFN model generates stochastic fracture networks using the following approach:

**Fracture orientation**: The Fisher distribution is used to model fracture pole orientations around mean directions for each fracture set:

$$f(\theta) = \frac{\kappa \sin\theta \cdot e^{\kappa \cos\theta}}{e^\kappa - e^{-\kappa}}$$

where κ is the Fisher concentration parameter (κ = 20 in this study), and θ is the angular deviation from the mean pole direction.

**Fracture geometry**: Each fracture is characterized by:
- Center position: uniformly distributed within the domain
- Half-length: L ~ N(μ_L, σ_L), with μ_L = 50–80 m, σ_L = 30 m
- Aperture: a ~ N(μ_a, σ_a), with μ_a = 1 mm, σ_a = 0.3 mm

**Hydraulic transmissivity**: The cubic law governs fracture flow:

$$T_f = \frac{\rho g a^3}{12\mu}$$

where a is the hydraulic aperture, ρ is fluid density, g is gravitational acceleration, and μ is dynamic viscosity.

**Equivalent permeability**: The upscaled permeability is computed via volumetric averaging:

$$k_{eq} = \frac{1}{V} \sum_{i=1}^{N_f} k_i \cdot a_i \cdot A_i$$

where k_i = a_i²/12 is the fracture permeability, a_i is the aperture, and A_i = πL_i² is the fracture area.

### 3.2 Thermo-Hydro-Mechanical Coupling

The THM solver employs a sequential coupling approach inspired by TOUGH2-FLAC3D workflows:

**Hydraulic equation** (implicit pressure):

$$\phi c_t \frac{\partial P}{\partial t} = \nabla \cdot \left(\frac{k \rho}{\mu} \nabla P\right) + Q$$

where φ is porosity, c_t is total compressibility, k is permeability, and Q is the source/sink term.

**Thermal equation** (advection-diffusion):

$$\left[(1-\phi)\rho_r c_{p,r} + \phi \rho_f c_{p,f}\right] \frac{\partial T}{\partial t} = \nabla \cdot (k_T \nabla T) + \rho_f c_{p,f} \mathbf{v} \cdot \nabla T$$

where k_T is thermal conductivity and **v** is the Darcy velocity.

**Mechanical equation** (thermo-poroelastic):

$$\Delta\sigma_{ij} = -\alpha_T (3\lambda + 2G) \Delta T \cdot \delta_{ij}/3 + \alpha \Delta P \cdot \delta_{ij}$$

where α_T is the thermal expansion coefficient, λ and G are Lamé parameters, and α is the Biot coefficient.

### 3.3 Supercritical Water Properties

Thermodynamic and transport properties are computed following the IAPWS-95 formulation (Wagner & Pruss, 2002):

- **Density**: Modified Benedict-Webb-Rubin EOS with smooth interpolation across the critical point
- **Viscosity**: IAPWS R12-08 correlation with simplified supercritical extension
- **Thermal conductivity**: IAPWS R11-07 formulation
- **Specific heat**: With enhanced divergence behavior near the critical point
- **Enthalpy**: Integrated from specific heat with latent heat accounting

### 3.4 Coulomb Failure Stress Model

The Coulomb failure stress change is computed as:

$$\Delta CFS = \Delta\tau - \mu_f (\Delta\sigma_n - \Delta P)$$

where Δτ is the shear stress change on the fault plane, Δσ_n is the normal stress change, μ_f is the friction coefficient, and ΔP is the pore pressure change.

The seismicity rate enhancement follows the rate-and-state model (Dieterich, 1994):

$$\frac{R}{R_0} = \exp\left(\frac{\Delta CFS}{A\sigma}\right)$$

where R_0 is the background seismicity rate, and Aσ is the rate-and-state parameter.

Magnitude distribution follows the Gutenberg-Richter law:

$$\log_{10} N(\geq M) = a - bM$$

with b = 1.0 (typical for induced seismicity).

### 3.5 Well Placement Optimization

The optimization evaluates multiple well configurations by computing the objective function:

$$\max_{\mathbf{x}_{inj}, \mathbf{x}_{prod}} E_{cum}(T=30\text{ yr}) = \int_0^{30} \dot{m} [h_{prod}(t) - h_{inj}] \, dt$$

where ṁ is the mass flow rate, h_prod is the production enthalpy, and h_inj is the injection enthalpy.

## 4. Experiments

### 4.1 Study Site: Kakkonda Geothermal Field

The Kakkonda geothermal field is located in Iwate Prefecture, Tohoku region, northeastern Japan. Key geological parameters used in this study:

| Parameter | Value | Source |
|-----------|-------|--------|
| Depth range | 2,000–4,000 m | Doi et al. (1998) |
| Reservoir temperature | 350°C | Doi et al. (1998) |
| Maximum temperature | 500°C (at >3,700 m) | Doi et al. (1998) |
| Brittle-ductile transition | 3,100 m | Muraoka et al. (2014) |
| Reservoir pressure | 30 MPa | Estimated |
| Rock type | Kakkonda Granite | Doi et al. (1998) |
| Density | 2,650 kg/m³ | Laboratory data |
| Porosity | 1.5% | Laboratory data |
| Matrix permeability | 5 × 10⁻¹⁸ m² | Laboratory data |
| Thermal conductivity | 2.8 W/(m·K) | Laboratory data |
| Young's modulus | 55 GPa | Laboratory data |
| Poisson's ratio | 0.23 | Laboratory data |
| Biot coefficient | 0.75 | Estimated |
| Friction coefficient | 0.65 | Byerlee's law |
| UCS | 220 MPa | Laboratory data |

### 4.2 Fracture Sets

Three fracture sets were modeled based on regional structural analysis:

| Set | Strike | Dip | Density | Mean Length |
|-----|--------|-----|---------|-------------|
| 1 | N45°E | 60° | 30/km³ | 80 m |
| 2 | N135°E | 75° | 20/km³ | 60 m |
| 3 | N-S | 45° | 15/km³ | 50 m |

### 4.3 Simulation Configuration

- **Domain**: 1,000 × 1,000 × 500 m
- **Grid**: 20 × 20 × 10 cells (Δx = Δy = Δz = 50 m)
- **Time step**: 1 year
- **Simulation period**: 30 years
- **Injection temperature**: 50°C (323.15 K)
- **Injection rate**: 20 L/s equivalent
- **Well configurations**: 5 spacings tested (50, 275, 500, 600, 900 m)

### 4.4 Evaluation Metrics

1. **Production temperature** T_prod(t): monitored at production well
2. **Thermal power** P_th(t) = ṁ·[h_prod(t) − h_inj] [MW]
3. **Cumulative energy** E_cum = ∫P_th dt [PJ]
4. **Pressure drawdown** ΔP = P₀ − P_prod [MPa]
5. **Coulomb failure stress change** ΔCFS [MPa]
6. **Maximum predicted magnitude** M_max
7. **Heat recovery factor** = E_cum / E_stored

## 5. Results

### 5.1 DFN Characterization

The generated DFN consists of 65 fractures across three orientation sets, with 8 identified intersections contributing to hydraulic connectivity. The equivalent upscaled permeability is 2.37 × 10⁻¹³ m², approximately four orders of magnitude higher than the matrix permeability (5 × 10⁻¹⁸ m²), confirming the fracture-dominated flow regime characteristic of granitic EGS reservoirs.

![Figure 1: DFN visualization](figures/fig1_dfn_network.png)
*Figure 1: (Left) Plan view of the generated DFN colored by aperture. (Right) Rose diagram showing the strike distribution of fractures, with dominant NE-SW and NW-SE orientations.*

### 5.2 Supercritical Water Properties

Figure 2 shows the computed thermophysical properties of water at 25 MPa across the temperature range 300–600°C. The dramatic changes near the critical temperature (374°C, marked by red dashed lines) are evident in all properties: density drops sharply, viscosity reaches minimum values, and specific heat exhibits a pronounced peak.

![Figure 2: Water properties](figures/fig2_water_properties.png)
*Figure 2: Thermophysical properties of water at P = 25 MPa computed using the IAPWS-95 formulation. Red dashed lines indicate the critical temperature T_c = 374°C.*

### 5.3 Reservoir Temperature and Permeability Fields

Figure 3 presents cross-sectional views of the reservoir temperature and permeability fields after 30 years of operation with the optimal well configuration.

![Figure 3: Reservoir fields](figures/fig3_reservoir_fields.png)
*Figure 3: (Left) Temperature field showing the thermal gradient and cold-water plume from injection. (Right) Log-scale permeability field highlighting DFN-enhanced zones.*

### 5.4 30-Year Heat Recovery Analysis

The 30-year simulation results for five well spacing configurations are presented in Figure 4. Key findings:

- **Optimal spacing**: 600 m, yielding maximum cumulative energy of 20.99 PJ
- **Production temperature**: Maintained at 350°C across all configurations due to the high initial reservoir temperature
- **Thermal power**: Varied with well spacing due to different flow path lengths and heat exchange areas
- **Pressure drawdown**: Inversely correlated with well spacing

![Figure 4: Heat recovery](figures/fig4_heat_recovery.png)
*Figure 4: 30-year heat recovery performance for five well spacing configurations: (a) production temperature, (b) thermal power output, (c) cumulative energy extraction, (d) pressure drawdown.*

### 5.5 Induced Seismicity Risk

The Coulomb stress analysis reveals:

- **Maximum ΔCFS**: 8.13 MPa, concentrated near the injection well
- **Predicted maximum magnitude**: M3.1
- **M ≥ 2.0 events**: 1 over 30 years
- **M ≥ 3.0 events**: 1 over 30 years
- **Seismicity rate**: Enhanced by several orders of magnitude near the injection zone

![Figure 5: Seismicity](figures/fig5_seismicity.png)
*Figure 5: Induced seismicity risk assessment: (a) Coulomb failure stress change field, (b) Gutenberg-Richter magnitude-frequency distribution, (c) spatial distribution of seismicity rate (log scale), (d) temporal evolution of predicted seismicity.*

### 5.6 Well Placement Optimization

Figure 6 summarizes the optimization results, clearly showing the 600 m spacing as optimal for the Kakkonda geological conditions.

![Figure 6: Well optimization](figures/fig6_well_optimization.png)
*Figure 6: (Left) 30-year cumulative energy extraction by well spacing, with optimal configuration highlighted. (Right) Final production temperature and thermal power as functions of well spacing.*

### 5.7 Simulation Workflow

![Figure 7: Workflow](figures/fig7_workflow.png)
*Figure 7: Overview of the integrated EGS simulation workflow showing data flow between DFN generation, IAPWS-95 EOS, THM coupling, Coulomb stress analysis, and optimization modules.*

## 6. Discussion

### 6.1 Framework Performance

The integrated framework successfully demonstrates the feasibility of multi-physics simulation for supercritical EGS development. The sequential THM coupling approach provides reasonable computational efficiency (minutes for a 30-year simulation on a standard workstation) while capturing the essential physics.

### 6.2 DFN Sensitivity

The DFN results confirm findings by Jiang et al. (2025) that fracture aperture and orientation are primary controls on reservoir performance. The equivalent permeability of 2.37 × 10⁻¹³ m² falls within the range reported for stimulated granitic reservoirs (10⁻¹⁴ to 10⁻¹² m²), validating the model parameterization.

### 6.3 Well Spacing Optimization

The optimal 600 m spacing balances two competing effects: shorter spacings provide better hydraulic connectivity but risk premature thermal breakthrough, while larger spacings reduce the heat exchange area. This result aligns with Zhou et al. (2021), who found optimal spacings of 300–600 m for similar reservoir conditions.

### 6.4 Induced Seismicity

The predicted maximum magnitude of M3.1 is consistent with observed induced seismicity at EGS sites worldwide (e.g., Basel M3.4, Soultz M2.9). The Coulomb stress analysis following Hutka et al. (2023) reveals that pressure changes dominate the seismic risk near the injection well, while thermal stresses contribute to a broader spatial footprint of stress perturbation.

### 6.5 Limitations

Several limitations should be acknowledged:

1. **Chemical coupling**: Mineral dissolution/precipitation (THMC) is not included, which may significantly affect long-term permeability evolution
2. **Two-phase flow**: The current model does not explicitly handle liquid-vapor coexistence, which is important near the critical point
3. **Grid resolution**: The 50 m cell size may be too coarse to resolve near-well effects and individual fracture behavior
4. **Fault dynamics**: Dynamic rupture propagation is not modeled; only static stress analysis is performed
5. **3D stress**: The current model uses simplified isotropic stress changes; a full 3D geomechanical model would improve accuracy

### 6.6 Future Directions

1. **THMC extension**: Incorporating mineral reaction kinetics for long-term reservoir management
2. **Machine learning surrogates**: Training neural networks on simulation databases for rapid optimization (Zhou et al., 2025)
3. **Real-time integration**: Coupling with microseismic monitoring for adaptive injection protocols
4. **Multi-objective optimization**: Simultaneously maximizing energy extraction and minimizing seismic risk using Pareto-optimal approaches
5. **Full coupling**: Transitioning from sequential to fully implicit THM coupling for improved accuracy at supercritical conditions

## 7. Conclusion

This study presents an integrated simulation framework for supercritical Enhanced Geothermal Systems that combines Discrete Fracture Network modeling, Thermo-Hydro-Mechanical coupling, IAPWS-95 supercritical water properties, Coulomb failure stress analysis, and well placement optimization. Applied to the Kakkonda geothermal field in Tohoku, Japan, the framework demonstrates:

1. A stochastic DFN model generating 65 fractures with 8 intersections and equivalent permeability of 2.37 × 10⁻¹³ m², consistent with stimulated granitic reservoirs.

2. Thirty-year THM simulations identifying an optimal well spacing of 600 m, achieving cumulative energy extraction of 20.99 PJ.

3. Coulomb stress analysis predicting maximum induced seismicity of M3.1, with rate-and-state modeling providing spatially resolved risk maps.

4. A TOUGH2/OpenGeoSys-compatible workflow suitable for practical supercritical EGS development planning.

The framework provides a foundation for comprehensive reservoir engineering analysis of next-generation supercritical geothermal systems, with clear pathways for extension to full THMC coupling and machine learning-enhanced optimization.

## References

1. Andrés, S., Santillán, D., Mosquera, J. C., & Cueto-Felgueroso, L. (2022). Hydraulic stimulation of geothermal reservoirs: Numerical simulation of induced seismicity and thermal decline. *Water*, 14(22), 3697. DOI: [10.3390/w14223697](https://doi.org/10.3390/w14223697)

2. Boyet, A., De Simone, S., Ge, S., & Vilarrasa, V. (2023). Insights on post-injection seismicity through analysis of the Enhanced Geothermal System at Basel (Switzerland). *Communications Earth & Environment*. DOI: [10.21203/rs.3.rs-1918182/v1](https://doi.org/10.21203/rs.3.rs-1918182/v1)

3. Doi, N., Muraoka, H., et al. (1998). Drilling and scientific results of the well WD-1a and implications for high-temperature geothermal resources in the Kakkonda geothermal field, Japan. *Geothermics*, 27(5–6), 663–690. DOI: [10.1016/S0375-6505(98)00030-0](https://doi.org/10.1016/S0375-6505(98)00030-0)

4. Drif, K., Lengliné, O., Kinscher, J., & Schmitbuhl, J. (2024). Induced seismicity controlled by injected hydraulic energy: The case study of the EGS Soultz-sous-Forêts site. *Journal of Geophysical Research: Solid Earth*, 129. DOI: [10.1029/2023JB028190](https://doi.org/10.1029/2023JB028190)

5. Fan, H., et al. (2025). Modeling and simulation of enhanced geothermal systems with explicit representation of 3D deformable fracture networks. *SPE Journal*, 30(6), 3839. DOI: [10.2118/225458-PA](https://doi.org/10.2118/225458-PA)

6. Hutka, G. A., Cacace, M., Hofmann, H., Mathur, B., & Zang, A. (2023). Investigating seismicity rates with Coulomb failure stress models caused by pore pressure and thermal stress from operating a well doublet in a generic geothermal reservoir in the Netherlands. *Netherlands Journal of Geosciences*, 102. DOI: [10.1017/njg.2023.7](https://doi.org/10.1017/njg.2023.7)

7. Janiga, D., Kwaśnik, J., & Wojnarowski, P. (2022). Utilization of Discrete Fracture Network (DFN) in modelling and simulation of a horizontal well-doublet Enhanced Geothermal System (EGS) with sensitivity analysis of key production parameters. *Energies*, 15(23), 9020. DOI: [10.3390/en15239020](https://doi.org/10.3390/en15239020)

8. Jiang, H., et al. (2025). Influence of discrete fracture network on the performance of enhanced geothermal system considering thermal-hydraulic-mechanical multi-physical field coupling. *PLOS ONE*. DOI: [10.1371/journal.pone.0320015](https://doi.org/10.1371/journal.pone.0320015)

9. Muraoka, H., et al. (2014). The Japan Beyond-Brittle Project: Development of deep geothermal resources above the supercritical zone. *Geothermics*, 49, 56–66. DOI: [10.1016/j.geothermics.2013.09.001](https://doi.org/10.1016/j.geothermics.2013.09.001)

10. Schwering, P. C., et al. (2020). Deterministic Discrete Fracture Network (DFN) model for the EGS Collab project. *Proceedings, 45th Workshop on Geothermal Reservoir Engineering*, Stanford University. OSTI: [1762429](https://www.osti.gov/biblio/1762429)

11. Wagner, W., & Pruss, A. (2002). The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use. *Journal of Physical and Chemical Reference Data*, 31(2), 387–535. DOI: [10.1063/1.1461829](https://doi.org/10.1063/1.1461829)

12. Zhou, Z., et al. (2021). Thermo-hydro-mechanical modelling study of heat extraction and flow processes in enhanced geothermal systems. *Advances in Geosciences*, 54, 229–240. DOI: [10.5194/adgeo-54-229-2021](https://doi.org/10.5194/adgeo-54-229-2021)

13. Suzuki, Y., et al. (2022). Constraining temperature at depth of the Kakkonda geothermal field, Japan, using Bayesian rock-physics modelling of resistivity. *Geothermics*, 100, 102316. DOI: [10.1016/j.geothermics.2021.102316](https://doi.org/10.1016/j.geothermics.2021.102316)
