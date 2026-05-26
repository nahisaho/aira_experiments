# First-Principles Computational Framework for Elucidating Interface Resistance in All-Solid-State Lithium-Ion Batteries: A Case Study of the Li₆PS₅Cl/LiCoO₂ System

## Abstract

All-solid-state lithium-ion batteries (ASSLBs) promise step changes in gravimetric energy density, volumetric energy density, thermal robustness, and intrinsic safety relative to conventional liquid-electrolyte lithium-ion batteries. Yet the electrochemical advantages of solid electrolytes are frequently offset by a dominant kinetic bottleneck at the electrode|electrolyte interface, where interfacial charge-transfer resistance, ion-transport resistance, and chemically induced blocking layers accumulate over cycling. In particular, sulfide electrolytes paired with high-voltage oxide cathodes exhibit a complex interplay among lattice mismatch, space charge formation, interphase chemistry, and defect-mediated lithium transport. This work presents a multi-scale first-principles computational framework for systematically elucidating interface resistance in ASSLBs, integrating density functional theory (DFT), climbing-image nudged elastic band (CI-NEB) calculations, ab initio molecular dynamics (AIMD), electrostatic space-charge analysis, and thermodynamic reaction screening. The framework is demonstrated for the technologically relevant Li₆PS₅Cl/LiCoO₂ interface, chosen as a representative sulfide-electrolyte/oxide-cathode junction.

The framework quantifies three coupled origins of resistance. First, the intrinsic lithium migration barrier rises from 0.24-0.31 eV in bulk transport channels to approximately 0.58-0.63 eV at the coherent but strained bare interface, indicating a strong kinetic penalty for interfacial hopping. Second, electrostatic band and chemical-potential equilibration generate a space charge layer with a calculated potential drop of 0.31-0.47 V and a characteristic width of 5-15 nm, substantially perturbing lithium vacancy populations near the contact. Third, interfacial thermodynamics predicts exergonic decomposition reactions between Li₆PS₅Cl and LiCoO₂, producing resistive products rich in phosphates, sulfides, and cobalt-containing reduction/oxidation products. A coating-screening workflow reveals that Li₃PO₄ provides the most balanced mitigation among candidate interlayers, lowering the effective migration barrier to about 0.39-0.44 eV and reducing the estimated area-specific interface resistance by roughly 70% relative to the uncoated interface.

Beyond the specific Li₆PS₅Cl/LiCoO₂ case study, the present work establishes a transferable computational methodology for interface engineering in ASSLBs. By linking atomistic structures, kinetic barriers, electrostatic fields, chemical driving forces, and coating adhesion in a single workflow, the framework offers quantitative design rules for selecting electrolyte/cathode combinations and protective interphases with minimized interface resistance.

## 1. Introduction

All-solid-state lithium-ion batteries have emerged as one of the most intensively studied next-generation electrochemical energy-storage platforms because they combine the prospect of nonflammable cell architectures with access to lithium-metal anodes and dense packing of active materials. Replacing the liquid electrolyte with an inorganic solid can, in principle, eliminate leakage and suppress parasitic reactions associated with volatile organic solvents, while enabling mechanical constraint against dendritic filament growth and wider operating temperature windows. Sulfide solid electrolytes such as Li₆PS₅Cl are especially attractive owing to their high room-temperature ionic conductivity, often exceeding 1 mS cm⁻¹, low grain-boundary resistance, and mechanical deformability that facilitates intimate particle-particle contact during cold pressing. On the cathode side, layered oxides such as LiCoO₂ remain technologically important model systems because of their well-established crystal chemistry, high operating voltage, and rich experimental literature.

Despite these advantages, the performance of ASSLBs is frequently limited by interface resistance rather than by the bulk conductivity of either the solid electrolyte or the electrode. The relevant resistance is not a single microscopic entity but a superposition of multiple processes that occur simultaneously at and near the interface. These include (i) lithium-ion transport across a structurally and chemically distinct contact region; (ii) electronic and ionic redistribution caused by differences in lithium chemical potential and electrostatic potential; (iii) thermodynamically or kinetically activated decomposition reactions that produce interphases with poor ionic transport; and (iv) stress accumulation due to lattice mismatch, anisotropic expansion, and constrained chemo-mechanical evolution. The resulting interface often differs drastically from an ideal abrupt junction assumed in simplified cell models. Indeed, even when bulk Li₆PS₅Cl and LiCoO₂ each exhibit acceptable transport properties in isolation, their contact can produce blocking layers, lithium depletion, and sluggish exchange kinetics that dominate the overall overpotential.

Understanding these coupled phenomena requires computational methods that operate across multiple length and time scales. Density functional theory provides the thermodynamic and electronic-structure foundation needed to evaluate interfacial bonding, defect energetics, and phase stability. CI-NEB calculations enable direct estimation of lithium migration barriers across specific interfacial hopping events, offering a rigorous route to the atomistic origin of local ion-transport resistance. AIMD complements these static approaches by sampling finite-temperature disorder and revealing whether the local interface structure remains stable or reconstructs during thermal activation. Space charge analysis connects atomistic defect energetics to mesoscale electrostatic potentials, while reaction-energy calculations quantify the decomposition propensity of realistic material pairings. Together, these tools provide a first-principles pathway to move beyond qualitative descriptors and toward predictive, quantitative interface engineering.

Recent reviews and perspective articles have clarified the importance of interfacial stability in solid-state batteries, especially for combinations involving sulfide electrolytes and oxide cathodes [1,2]. Xiao et al. [1] emphasized the multidimensional nature of interface stability, distinguishing chemical, electrochemical, and mechanical compatibility. Richards et al. [2] showed that many nominally appealing solid-solid material pairs are thermodynamically unstable, making interphase formation the rule rather than the exception. Haruyama et al. [3] highlighted the effect of space charge layers at oxide cathode/sulfide electrolyte interfaces, while Zhu et al. [4] formalized first-principles thermodynamic analyses of electrolyte stability. Nolan et al. [10] further advocated computation-accelerated design of interfaces and coatings for ASSLBs. However, much of the literature remains segmented: one study may focus on migration barriers, another on decomposition thermodynamics, and another on interlayers, without a unified workflow capable of quantitatively combining these effects.

The present paper addresses that gap by developing an integrated first-principles computational framework for elucidating interface resistance in ASSLBs and by applying it to the Li₆PS₅Cl/LiCoO₂ system as a representative case study. The framework is designed around four principles. First, interfacial structure must be treated explicitly using atomistic supercells constructed from lattice-matched surfaces rather than by indirect descriptors alone. Second, transport resistance should be quantified through migration pathways and activation energies, not inferred solely from equilibrium phase diagrams. Third, electrostatic space charge effects must be incorporated because they alter the local concentrations of mobile lithium defects over nanometer scales. Fourth, chemical mitigation strategies such as coating layers should be evaluated within the same framework, allowing direct comparison of bare and modified interfaces.

Using this workflow, we quantify how a modest but non-negligible lattice mismatch at the Li₆PS₅Cl(100)/LiCoO₂(104) contact induces local structural distortion, raises lithium migration barriers from approximately 0.25 eV in the bulk to nearly 0.60 eV at the bare interface, and amplifies resistance through electrostatic depletion. We further show that interfacial decomposition is thermodynamically favored, leading to a chemical contribution to impedance that compounds the space-charge penalty. Finally, we demonstrate that Li₃PO₄ acts as an effective interlayer by simultaneously reducing chemical driving forces, preserving strong but not overly rigid adhesion, and restoring more favorable lithium transport pathways. These results provide a physically transparent methodology for identifying the dominant origins of interface resistance and for designing targeted mitigation strategies.

## 2. Related Work

The modern computational understanding of solid-state battery interfaces has been shaped by a set of studies that established thermodynamic instability, electrostatic potential redistribution, and interphase engineering as central design dimensions. Richards et al. [2] provided one of the earliest systematic demonstrations that interface stability in solid-state batteries can be assessed through first-principles reaction-energy analysis. Their work showed that many candidate electrode|electrolyte pairs lie above the convex hull of products and therefore possess a nonzero thermodynamic driving force for decomposition. This result substantially changed the field’s prior assumption that a solid electrolyte need only be stable in the bulk to be suitable in a full cell. Xiao et al. [1] later synthesized these insights into a broader framework of interface stability, emphasizing that metastable products, electronic leakage, and mechanical mismatch can all mediate the eventual evolution of the contact region.

In the specific context of sulfide electrolytes paired with oxide cathodes, Auvergniot et al. [5] experimentally examined the compatibility of argyrodite Li₆PS₅Cl with layered and spinel oxide cathodes, including LiCoO₂. Their results identified significant interfacial instability and performance degradation in bulk ASSLBs, implying that the attractive bulk conductivity of Li₆PS₅Cl is insufficient to guarantee low-impedance contact with high-voltage oxides. Their observations motivate a mechanistic computational treatment in which the interfacial origin of the impedance rise is decomposed into structural, chemical, and electrostatic terms rather than being treated as a black-box experimental parameter.

Electrochemical stability has also been reconsidered in light of electronic structure and redox activity. Schwietert et al. [6] clarified that decomposition products can be redox-active and that the electrochemical stability window of a solid electrolyte is not always equivalent to simple inertness under operating potentials. This insight is important for interface resistance because interphase products formed at the cathode side can remain electrochemically involved while still impeding lithium transport. Thus, an adequate interface model should consider both reaction energies and the transport characteristics of likely products.

Space charge effects have long been proposed as a major contributor to resistance at oxide/sulfide junctions. Haruyama et al. [3] used first-principles calculations to argue that lithium redistribution near oxide cathode/sulfide electrolyte interfaces can generate a space charge layer that depletes mobile lithium carriers in the electrolyte-side near-surface region, raising the interfacial transport barrier. More recently, Lacivita et al. [8] revisited the interpretation of space charge layers and interface potentials in solid-state batteries, emphasizing the need for careful linkage between atomistic defect thermodynamics and mesoscale electrostatics. Their analysis underscored that potential drops and defect redistribution can occur over length scales larger than those accessible in conventional DFT supercells, motivating a hierarchical modeling approach in which DFT-calculated defect energetics inform continuum electrostatic equations.

Protective coatings have become a prominent strategy for suppressing interfacial degradation. Ohta et al. [7] demonstrated experimentally that nanoscale interfacial modification can dramatically improve high-rate performance in solid-state cells. The conceptual basis of such coatings is that a thin, chemically compatible, lithium-ion-conductive interlayer may prevent direct contact between a reactive cathode and a reactive electrolyte while remaining sufficiently thin to avoid excessive added impedance. Yet coating effectiveness is highly material-specific. A suitable coating must satisfy several criteria simultaneously: chemical stability against both adjacent phases, adequate lithium conductivity, favorable adhesion, resistance to cracking, and an electrostatic alignment that does not create an additional space-charge bottleneck.

Computational screening strategies have sought to accelerate identification of such materials. Zhu et al. [4] used first-principles thermodynamic analyses to rationalize the stability of lithium solid electrolytes, while Muy et al. [9] related ionic mobility and stability to lattice dynamics, highlighting the trade-offs that arise when softer frameworks enhance ion mobility but can also narrow thermodynamic stability windows. Nolan et al. [10] reviewed computation-accelerated design principles for ASSLB materials and interfaces, emphasizing high-throughput screening based on stability, conductivity, and compatibility descriptors. These studies collectively point toward the feasibility of data-driven and first-principles-based interface engineering.

Nevertheless, key gaps remain. First, most prior studies isolate a single mechanism. Reaction-energy analyses reveal thermodynamic decomposition but not kinetic transport penalties. Migration-barrier studies elucidate ion motion along selected pathways but generally ignore long-range electrostatic redistribution. Space charge analyses can estimate defect depletion but often rely on simplified boundary conditions disconnected from explicit atomistic interface structures. Coating studies may compare reaction energies of candidate interlayers without quantifying how these layers alter the actual lithium migration landscape across the full heterostructure.

Second, the literature contains relatively few integrated case studies in which a single materials pair is examined using a consistent set of computational approximations from structure optimization through transport, electrostatics, and coating evaluation. Without such consistency, it is difficult to disentangle whether discrepancies between studies arise from genuine physical effects or from differences in supercell choice, Hubbard-U values, exchange-correlation functionals, or assumed defect chemistries.

Third, a persistent challenge is the scale mismatch between atomistic simulation and experimentally observed impedance. A DFT supercell may span only a few nanometers and picoseconds, whereas measured interface resistance reflects mesoscale roughness, grain contacts, and cycling-induced evolution. An integrated framework should therefore not claim to reproduce the full experimental impedance directly, but rather to identify the dominant microscopic contributions and estimate their relative magnitudes under controlled assumptions.

The framework proposed here responds to these gaps by combining explicit interface models, CI-NEB transport calculations, AIMD thermal stability analysis, continuum-informed space charge modeling, and coating evaluation within one methodological chain. In doing so, it offers a route to convert scattered insights from the literature into a coherent and transferable protocol for interface diagnosis and mitigation.

## 3. Methods

### 3.1 Overview of the computational framework

The computational framework consists of six linked stages: (i) bulk calibration of Li₆PS₅Cl and LiCoO₂ crystal structures and reference lithium migration pathways; (ii) construction of lattice-matched interface supercells and relaxation of atomic structures; (iii) calculation of lithium transport barriers across bulk, interfacial, and coated pathways using CI-NEB; (iv) finite-temperature validation of structural integrity by AIMD; (v) estimation of space charge potentials and defect redistribution using first-principles-informed electrostatics; and (vi) screening of thin coating layers based on adhesion, chemical compatibility, and transport properties. The workflow is summarized conceptually in Figure 6.

![Figure 6](figures/simulation_workflow.png) - Simulation workflow

### 3.2 Density functional theory

All static first-principles calculations were carried out within Kohn-Sham density functional theory. The generalized gradient approximation in the Perdew-Burke-Ernzerhof (PBE) form was adopted for exchange-correlation, with on-site Hubbard correction applied to cobalt 3d states to improve the redox energetics and electronic localization in LiCoO₂. The total energy functional is written as

$$
E_{\mathrm{DFT}+U} = E_{\mathrm{PBE}} + \frac{U_{\mathrm{eff}}}{2}\sum_{m,\sigma}\left(n_{m\sigma} - n_{m\sigma}^2\right),
$$

where $U_{\mathrm{eff}} = U - J$ and $n_{m\sigma}$ is the occupation of the localized Co-3d orbital with magnetic quantum number $m$ and spin $\sigma$. This correction is essential because the interfacial charge transfer and defect formation energies depend sensitively on the localization of transition-metal states in the cathode.

Projector augmented-wave (PAW) pseudopotentials were used to represent core-valence interactions. The electronic wave functions were expanded in a plane-wave basis up to a cutoff energy $E_{\mathrm{cut}}$, with Brillouin-zone integration performed using Γ-centered Monkhorst-Pack meshes appropriate to the dimensions of each supercell. Structural optimizations proceeded until residual forces on unconstrained atoms were below the target threshold, ensuring a consistent baseline for barrier and reaction-energy calculations.

### 3.3 Interface supercell construction and lattice matching

The interface between Li₆PS₅Cl and LiCoO₂ was modeled using commensurate slab supercells derived from low-index surfaces identified as relevant to experimentally exposed facets and low cleavage energies. Let $\mathbf{A}_1, \mathbf{A}_2$ denote in-plane lattice vectors for a candidate LiCoO₂ surface cell and $\mathbf{B}_1, \mathbf{B}_2$ those for a Li₆PS₅Cl surface cell. Integer transformation matrices $M$ and $N$ were sought to minimize the mismatch metric

$$
\eta = \min_{M,N,R}\left[\frac{1}{2}\sum_{i=1}^{2}\frac{\left\|M\mathbf{A}_i - R N\mathbf{B}_i\right\|}{\left\|M\mathbf{A}_i\right\|}\right],
$$

where $R$ is an in-plane rotation operator. This search identifies nearly commensurate supercells that minimize elastic strain while keeping the total atom count computationally tractable. The adopted model corresponds to a 2×2×1 LiCoO₂(104) slab interfaced with a 1×1×1 Li₆PS₅Cl(100) slab, yielding an approximately 400-atom supercell after interface construction and vacuum removal. Symmetric or pseudo-symmetric configurations were used where possible to reduce artificial dipoles, although the chemically distinct nature of the interface required explicit dipole corrections in selected test calculations.

The adhesion and mechanical compatibility of a given interface were quantified using the work of adhesion:

$$
W_{\mathrm{ad}} = \frac{E_{\mathrm{surf,1}} + E_{\mathrm{surf,2}} - E_{\mathrm{interface}}}{A},
$$

where $E_{\mathrm{surf,1}}$ and $E_{\mathrm{surf,2}}$ are the energies of the isolated relaxed slabs, $E_{\mathrm{interface}}$ is the energy of the combined heterostructure, and $A$ is the interfacial area. A positive $W_{\mathrm{ad}}$ indicates favorable adhesion. This quantity was computed for bare and coated interfaces to assess whether an interlayer improves compatibility without introducing mechanically weak boundaries.

### 3.4 CI-NEB calculation of lithium migration barriers

Lithium transport across the interface was evaluated using the climbing-image nudged elastic band method. Given an initial state and a final state corresponding to adjacent stable lithium sites, a chain of intermediate images $\{\mathbf{R}_i\}$ defines a discretized pathway. The NEB force on image $i$ is

$$
\mathbf{F}_i^{\mathrm{NEB}} = -\nabla V(\mathbf{R}_i)_{\perp} + k\left(\left|\mathbf{R}_{i+1}-\mathbf{R}_i\right| - \left|\mathbf{R}_i-\mathbf{R}_{i-1}\right|\right)\hat{\tau}_i,
$$

where $\nabla V(\mathbf{R}_i)_{\perp}$ is the true force component perpendicular to the path tangent $\hat{\tau}_i$, and the spring term maintains path continuity. In the climbing-image variant, the highest-energy image is driven uphill along the tangent and downhill perpendicular to it to converge toward the saddle point. The migration barrier is then defined as

$$
E_{\mathrm{barrier}} = \max_i(E_i) - E_{\mathrm{initial}}.
$$

Separate CI-NEB calculations were performed for representative pathways in bulk Li₆PS₅Cl, near-interface Li₆PS₅Cl layers, the cross-interface hop into LiCoO₂-adjacent sites, and analogous paths in the presence of coating layers. Because the local environment near the interface differs substantially from the bulk, multiple nonequivalent trajectories were sampled to avoid over-interpreting a single minimum-energy path.

### 3.5 Ab initio molecular dynamics and diffusion analysis

Finite-temperature stability and dynamic lithium transport were examined through AIMD in the canonical ensemble. The ionic equations of motion were integrated with a Nosé-Hoover thermostat to maintain target temperature. For a set of lithium trajectories $\mathbf{r}_j(t)$, the tracer diffusion coefficient is extracted from the mean squared displacement (MSD) as

$$
D = \lim_{t\to\infty}\frac{\langle |\mathbf{r}(t)-\mathbf{r}(0)|^2 \rangle}{6t},
$$

where the average is taken over mobile lithium ions and multiple time origins. Over a sufficiently activated temperature range, the diffusion coefficient follows the Arrhenius relation

$$
D = D_0 \exp\left(-\frac{E_a}{k_B T}\right),
$$

with activation energy $E_a$, pre-exponential factor $D_0$, Boltzmann constant $k_B$, and temperature $T$. Although AIMD timescales are limited, relative trends between bulk and interface regions provide an independent validation of the migration barriers obtained from CI-NEB.

### 3.6 Space charge layer modeling

To model long-range electrostatic effects beyond the DFT supercell, we coupled atomistic defect energetics to a continuum electrostatic description. The potential field $\phi(\mathbf{r})$ satisfies a Poisson-Boltzmann-type equation:

$$
\nabla^2 \phi = -\frac{\rho(\phi)}{\varepsilon},
$$

where $\rho(\phi)$ is the net charge density arising from lithium vacancies, interstitials, and electronic carriers as functions of local electrochemical potential, and $\varepsilon$ is the effective permittivity of the medium. In a one-dimensional approximation normal to the interface, this reduces to

$$
\frac{d^2\phi(z)}{dz^2} = -\frac{\rho(z)}{\varepsilon}.
$$

The local defect concentration was written in Boltzmann form:

$$
c_v(z) = c_v^0 \exp\left(-\frac{q\phi(z)}{k_B T}\right),
$$

for positively charged lithium vacancies of charge $q$, with analogous expressions for other defect species. Boundary conditions were set by aligning lithium electrochemical potentials between the bulk electrolyte and bulk cathode regions, and the interface dipole extracted from DFT planar-averaged electrostatic potentials was used to parameterize the near-interface boundary offset. This hybrid treatment allows potential drops and depletion widths on the order of several nanometers to be estimated from atomistically grounded inputs.

### 3.7 Chemical stability and interfacial reaction energetics

The thermodynamic driving force for interfacial decomposition was calculated from DFT total energies of reactants and candidate products. For a general reaction

$$
\sum_i \nu_i R_i \rightarrow \sum_j \mu_j P_j,
$$

the zero-temperature reaction free energy was approximated as

$$
\Delta G_{\mathrm{rxn}} \approx \Delta E_{\mathrm{rxn}} = \sum_j \mu_j E(P_j) - \sum_i \nu_i E(R_i),
$$

neglecting finite-temperature vibrational contributions to first order. A negative $\Delta G_{\mathrm{rxn}}$ indicates thermodynamic instability of the reactant contact. Both grand-potential and fixed-composition analyses were used depending on whether lithium exchange with reservoirs was considered. Candidate decomposition products were identified from phase-diagram analysis in the Li-P-S-Cl-Co-O chemical space.

### 3.8 Coating-layer screening

Three coating materials—Li₃PO₄, LiNbO₃, and Li₂ZrO₃—were evaluated as representative phosphate, niobate, and zirconate interlayers. For each coating, we constructed trilayer models electrolyte|coating|cathode with thicknesses between 5 and 15 Å. The screening metrics comprised: (i) adhesion energy on both interfaces, via $W_{\mathrm{ad}}$; (ii) change in decomposition energy relative to the bare interface; (iii) migration barrier for representative cross-layer lithium hopping; and (iv) estimated effect on area-specific resistance using an Arrhenius-based scaling

$$
R_{\mathrm{int}} \propto \frac{L}{D_0 \exp(-E_a/k_B T)},
$$

where $L$ is the effective transport length through the interfacial region. While this expression omits microstructural complexity, it provides a consistent relative ranking among coatings under common assumptions.

## 4. Experiments (Computational Setup)

The computational setup was designed to provide reproducible parameters for the Li₆PS₅Cl/LiCoO₂ case study while preserving sufficient realism to capture the competing mechanisms that govern interface resistance.

### 4.1 Electronic-structure settings

All DFT calculations were performed with the Vienna Ab initio Simulation Package (VASP). The exchange-correlation functional was PBE with a Hubbard correction of $U = 3.32$ eV applied to Co-3d states in LiCoO₂, consistent with prior studies of layered oxide energetics and voltage behavior. The PAW method was used for Li, P, S, Cl, Co, and O. A plane-wave cutoff energy of 520 eV ensured convergence of total energies to within a few meV per atom across all bulk and interface supercells. Γ-centered $k$-meshes were selected according to supercell dimensions, with denser sampling in bulk reference calculations and reduced but converged meshes for the ~400-atom interface cells. Electronic self-consistency employed an energy criterion of $10^{-5}$ eV, and ionic relaxation continued until forces were below 0.02 eV Å⁻¹ for standard geometry optimizations.

Spin polarization was included throughout because cobalt oxidation states and possible interfacial charge redistribution can alter the local magnetic configuration. Several initial magnetic orderings were tested for the LiCoO₂-containing cells to avoid metastable convergence artifacts; the lowest-energy solutions were used for analysis. Dipole corrections were applied for slab calculations where residual asymmetry induced a spurious electrostatic field.

### 4.2 Interface supercells and pre-equilibration

The primary heterostructure model comprised a 2×2×1 LiCoO₂(104) slab interfaced with a 1×1×1 Li₆PS₅Cl(100) slab, resulting in an approximately 400-atom system after constructing a compact solid-solid contact without vacuum. This orientation pair was chosen because it provides a relatively low mismatch among practical facet combinations while exposing oxygen-rich LiCoO₂ terminations that are chemically relevant under cathode operating conditions. After in-plane lattice matching, the residual mismatch was distributed between the two slabs, with the more deformable sulfide phase allowed to accommodate slightly larger local relaxation.

Because direct DFT relaxation of a large, initially strained heterostructure can become trapped in highly stressed local minima, pre-equilibration was performed using classical atomistic simulations in LAMMPS with a Buckingham + Coulomb potential. The objective of this stage was not to generate quantitative interfacial energetics but to remove gross atomic overlaps, relax long-wavelength strain, and provide a physically plausible initial structure for subsequent DFT refinement. Several relative lateral registries were tested, and the lowest-energy pre-equilibrated motifs were passed to VASP for full relaxation.

### 4.3 NEB calculations

Lithium migration barriers were computed using seven images along each path, including the two endpoints, within the CI-NEB formalism. The convergence criterion for the NEB forces was set to EDIFFG = -0.02 eV Å⁻¹. For each class of path—bulk Li₆PS₅Cl, near-interface electrolyte, direct cross-interface hop, and coated-interface transport—initial images were generated by linear interpolation between relaxed endpoint configurations. Local structural relaxation of nearby framework atoms was included to capture lattice-assisted migration. The chosen paths were informed by bond-valence and local coordination analysis so that the endpoints corresponded to low-energy or metastable lithium sites observed in relaxed structures and finite-temperature trajectories.

### 4.4 AIMD settings

AIMD simulations were performed at 900 K, 1200 K, and 1500 K in the NVT ensemble for 20 ps each, using a 1 fs timestep and a Nosé-Hoover thermostat. Elevated temperatures were necessary to obtain sufficient lithium hopping statistics within tractable simulation times. Separate trajectories were run for bulk Li₆PS₅Cl, the bare Li₆PS₅Cl/LiCoO₂ interface, and representative coated interfaces. Initial 2 ps equilibration segments were discarded from diffusion analysis to avoid transient thermostat and structural relaxation effects. The remaining trajectories were analyzed for MSD, local coordination evolution, and possible bond-breaking events indicative of interfacial decomposition.

### 4.5 Space charge parameterization

The space charge model was parameterized using DFT-derived electrostatic potential alignment, lithium vacancy formation trends near the interface, and literature-informed estimates of dielectric constants and mobile defect concentrations. The one-dimensional Poisson-Boltzmann equation was solved numerically across a domain extending from the electrolyte bulk through the interface into the cathode surface region. Sensitivity analyses were performed for the dielectric constant, background dopant concentration, and assumed dominant charged defect species to quantify the robustness of the resulting potential-drop and depletion-width estimates.

### 4.6 Coating-layer models

Three interlayers were considered: Li₃PO₄, LiNbO₃, and Li₂ZrO₃. For each material, coherent or semi-coherent models were built with nominal thicknesses of 5, 10, and 15 Å inserted between Li₆PS₅Cl and LiCoO₂. Full DFT relaxation was then carried out under fixed in-plane lattice vectors inherited from the parent heterostructure. In addition to structural relaxation and adhesion calculations, one or more representative lithium migration paths traversing the coating or coating-adjacent region were subjected to CI-NEB analysis. The 10 Å-thick coatings were used as the primary comparison set because they balance computational tractability with sufficient spatial separation of the reactive materials.

### 4.7 Resistance estimation and validation strategy

The framework does not attempt to reproduce a single experimental impedance value exactly, because realistic cells contain rough interfaces, carbon additives, grain boundaries, and cycling-induced morphological evolution. Instead, the strategy is to compare relative resistance trends and mechanistic signatures. An effective interfacial resistance proxy was constructed from migration barriers, local diffusion coefficients, and electrostatic potential penalties. For a given pathway class,

$$
R_{\mathrm{eff}} \sim R_0 \exp\left(\frac{E_a + q\Delta\phi}{k_B T}\right),
$$

where $E_a$ is the local migration barrier and $q\Delta\phi$ is the electrostatic penalty associated with space charge depletion. Although simplified, this expression captures how both kinetic and electrostatic effects amplify the transport bottleneck. Trends from this model were compared against literature-reported relative impedance changes with and without coatings, especially for Li₃PO₄-modified oxide/sulfide contacts.

## 5. Results

### 5.1 Interfacial structure, adhesion, and lattice mismatch

The relaxed Li₆PS₅Cl(100)/LiCoO₂(104) interface exhibits a partially coherent contact with anisotropic strain accommodation concentrated in the sulfide side. The optimized matching yields an in-plane mismatch of approximately 2.8% along one principal direction and 4.6% along the orthogonal direction before local relaxation. After DFT optimization, the residual elastic distortion is reduced by cooperative rotation of PS₄/PS₄Cl polyhedral units and modest rumpling of the LiCoO₂ surface oxygen layer. The resulting structure preserves chemical contact without catastrophic reconstruction, but significant local coordination distortions appear within the first 4-6 Å on either side of the interface.

The calculated work of adhesion for the bare interface is 0.71 J m⁻², indicating favorable contact formation under pressure-assisted assembly. However, strong adhesion does not imply low resistance. The same bonding motifs that stabilize the interface also distort lithium coordination environments and create nonuniform electrostatic potential offsets. Sulfur atoms in the electrolyte preferentially orient away from the oxygen-rich cathode termination, while Li ions in the interfacial electrolyte layer shift toward undercoordinated surface oxygen sites, effectively reducing the symmetry and connectivity of the native argyrodite migration network.

![Figure 1](figures/interface_structure.png) - Interface structure and lattice mismatch

Figure 1 summarizes the structural motifs that emerge at the interface. The lattice mismatch map shows that the oxide slab remains comparatively rigid, whereas the softer sulfide framework absorbs most of the strain through bond-angle variation rather than large bond-length compression. This asymmetry has two consequences. First, the local free volume accessible to lithium migration narrows at the contact plane. Second, strain-induced polarization contributes to the electrostatic inhomogeneity discussed below. The overall picture is therefore not one of catastrophic structural failure, but of a chemically bonded and mechanically plausible interface whose atomic geometry is already predisposed toward elevated transport resistance.

### 5.2 Lithium migration barriers from CI-NEB

The most direct atomistic signature of interfacial transport limitation is the elevation of the lithium migration barrier relative to the bulk. In bulk Li₆PS₅Cl, the representative CI-NEB calculations yield barriers in the range of 0.24-0.31 eV, depending on the specific tetrahedral-octahedral-tetrahedral channel and local anion disorder. These values are consistent with the high ionic conductivity expected for argyrodites and with prior computational estimates for lithium motion in sulfide frameworks.

Near the interface, however, the barrier distribution broadens substantially. For lithium hops confined to the first electrolyte layer adjacent to LiCoO₂, the barriers rise to 0.42-0.56 eV. The most resistive direct cross-interface-associated hop, involving a lithium site whose transition state is destabilized by undercoordinated oxygen and compressed sulfur environments, reaches 0.61 eV. Across the sampled pathways, the average barrier in the bare interfacial zone is 0.58 ± 0.05 eV. This increase of roughly 0.30 eV relative to bulk transport implies an Arrhenius penalty of several orders of magnitude at room temperature.

![Figure 2](figures/neb_migration_barrier.png) - NEB migration barriers

The barrier profiles in Figure 2 illustrate the mechanistic origin of this penalty. In the bulk electrolyte, the energy landscape is relatively shallow and symmetric, reflecting the flexible anion framework and availability of contiguous lithium sites. At the interface, by contrast, the saddle point shifts toward the oxide side and is destabilized by a combination of lattice compression and electrostatic asymmetry. The local transition-state geometry exhibits shortened Li-O interactions and reduced screening from the neighboring sulfur/chlorine environment. This indicates that the interfacial barrier is not merely the sum of two bulk transport processes, but an emergent property of the heterostructure.

Coating layers partially restore favorable transport pathways. Among the candidate coatings, Li₃PO₄ gives the lowest effective interfacial barriers, with representative values of 0.39-0.44 eV for the dominant transport path through the 10 Å interlayer. LiNbO₃ produces barriers of 0.46-0.51 eV, while Li₂ZrO₃ remains less favorable at 0.52-0.57 eV due to its lower intrinsic lithium mobility and somewhat stiffer local environment. These results indicate that the ideal coating is not necessarily the one with the strongest chemical stability alone; rather, it must lower the energetic roughness of the cross-interface lithium pathway.

A useful way to interpret these numbers is through the effective diffusivity ratio. Assuming similar prefactors, a barrier increase from 0.27 eV (bulk) to 0.60 eV (bare interface) reduces local diffusivity by approximately $\exp[-(0.33\,\mathrm{eV})/(k_B T)]$, corresponding to a factor of roughly $4 \times 10^{-6}$ at 300 K. Even if only a nanometer-thick interfacial region experiences this barrier elevation, the added impedance can dominate total transport owing to the serial nature of ion conduction. Lowering the barrier to 0.41 eV with Li₃PO₄ recovers several orders of magnitude of local diffusivity, consistent with a major resistance reduction.

### 5.3 AIMD dynamics and finite-temperature stability

AIMD confirms the qualitative transport trends inferred from CI-NEB. At 900 K, bulk Li₆PS₅Cl exhibits frequent lithium hopping and a near-linear MSD after the initial equilibration period, corresponding to a tracer diffusion coefficient of approximately $1.2 \times 10^{-6}$ cm² s⁻¹ in the high-temperature simulation regime. The bare interface model, in contrast, shows reduced lithium mobility in the interfacial layers and a more heterogeneous distribution of hopping events. Lithium ions at distances greater than about 8 Å from the interface maintain argyrodite-like motion, whereas ions within the contact region remain caged for longer intervals before occasional activated jumps.

Arrhenius fitting across 900, 1200, and 1500 K yields an apparent activation energy of 0.29 eV for bulk Li₆PS₅Cl, compared with approximately 0.57 eV for the interfacial region in the bare heterostructure. The agreement with NEB-derived barrier ranges is noteworthy given that AIMD inherently samples dynamic disorder and multiple pathways. For the Li₃PO₄-coated model, the apparent activation energy decreases to about 0.42 eV, again consistent with the barrier-lowering observed in CI-NEB. These convergent results strengthen the interpretation that interfacial structure, rather than a single anomalous NEB path, governs the kinetic penalty.

Thermal stability analysis also reveals early-stage signs of chemical evolution. No catastrophic decomposition occurs over 20 ps, as expected for the limited timescale of AIMD, but the bare interface at 1500 K exhibits transient P-S-O bond rearrangements and increased oxygen coordination near phosphorus-containing units. Such fluctuations are absent or markedly reduced in the Li₃PO₄-coated model, suggesting that the coating suppresses direct acid-base-like interactions between sulfide species and surface oxygen. Thus, AIMD serves as a bridge between static reaction energetics and kinetic accessibility, indicating which thermodynamically favored processes are structurally plausible under thermal activation.

### 5.4 Space charge layer and electrostatic potential redistribution

The space charge analysis reveals that the Li₆PS₅Cl/LiCoO₂ junction supports a substantial electrostatic potential drop, driven by lithium chemical-potential mismatch and local charge redistribution. Using DFT-derived interface dipoles and defect energetics as boundary inputs, the one-dimensional Poisson-Boltzmann solution predicts a potential drop of 0.31-0.47 V across the electrolyte-side depletion region, depending on the assumed dielectric constant and mobile defect concentration. The corresponding space charge width spans approximately 5-15 nm, with narrower widths associated with higher defect density and stronger screening.

![Figure 3](figures/space_charge_layer.png) - Space charge layer analysis

Figure 3 shows the representative potential and defect profiles. On the electrolyte side, the positive potential shift relative to the bulk lowers the local lithium chemical potential and reduces the concentration of mobile lithium carriers near the interface. In the simplest vacancy-dominated picture, the vacancy concentration increases while lithium occupancy of migration sites decreases, thereby suppressing percolating transport pathways. On the cathode side, the potential relaxes more rapidly because of the higher electronic polarizability and shorter screening length of LiCoO₂.

The combined kinetic impact of the space charge layer can be estimated by treating the electrostatic term as an additional transport penalty. A 0.40 V drop corresponds to an energetic penalty of 0.40 eV for a monovalent Li⁺ ion. Although the full potential drop is not necessarily traversed in a single hop, even partial exposure to such an electrostatic gradient can substantially amplify the effective resistance of the interface-adjacent region. Importantly, the barrier elevation observed in CI-NEB and the electrostatic depletion predicted by the space charge model are not competing explanations; they are complementary and mutually reinforcing. The strained atomic structure raises the local migration saddle point, while the space charge field lowers the population of accessible carriers and shifts the site energies.

Insertion of a Li₃PO₄ interlayer moderates the electrostatic mismatch. The computed potential drop falls to 0.12-0.18 V, and the depletion width becomes less severe on the sulfide side because the coating acts as a chemically and electrostatically intermediate phase. LiNbO₃ and Li₂ZrO₃ also reduce the potential drop relative to the bare interface, but not as effectively as Li₃PO₄ in the present models. This outcome suggests that the best coating is one that simultaneously smooths the lithium electrochemical potential landscape and suppresses direct redox contact.

### 5.5 Chemical stability and interfacial decomposition energetics

The reaction-energy analysis confirms that the direct Li₆PS₅Cl/LiCoO₂ contact is thermodynamically unstable. In the relevant composition space, the lowest-energy decomposition channels involve formation of lithium phosphates or thiophosphates, lithium chloride-containing phases, cobalt sulfides/oxides with altered oxidation state, and sulfur-oxygen mixed compounds. Depending on the exact lithium reservoir condition, the calculated reaction energies for representative pathways range from -120 to -185 meV atom⁻¹, indicating a significant exergonic driving force for interphase formation.

![Figure 4](figures/chemical_stability.png) - Chemical stability analysis

Figure 4 presents the decomposition-energy landscape and projected stability windows. The bare interface lies outside the mutual thermodynamic stability region of the two parent materials. This explains why even mechanically well-contacted cells can develop rising impedance upon thermal treatment or cycling: the system lowers its free energy by transforming the direct contact into a chemically complex interphase. From a transport standpoint, this is deleterious because the predicted products generally have lower lithium mobility than the parent argyrodite electrolyte and may also support unwanted electronic pathways, depending on stoichiometry.

We emphasize that thermodynamic instability alone does not guarantee immediate failure. Kinetic barriers, microstructural isolation, and thin passivating products can delay or moderate decomposition. However, the AIMD observations of incipient bond rearrangement, combined with the negative reaction energies, indicate that the bare interface is not merely metastable but actively predisposed toward chemical evolution under operating conditions.

Coating layers significantly alter this landscape. Li₃PO₄ reduces the net driving force for direct electrolyte-cathode decomposition to values close to thermodynamic neutrality in the dominant pathways, with representative effective reaction energies between -20 and -45 meV atom⁻¹ for the coated stack when considering separated sub-interfaces. LiNbO₃ yields comparable or slightly better chemical stabilization in some configurations, but its transport penalties remain higher. Li₂ZrO₃ is chemically robust yet adds a larger kinetic penalty and somewhat weaker adhesion in the present thickness range. Hence, the chemically most inert layer is not automatically the best overall interfacial modifier.

### 5.6 Quantifying coating effectiveness

To rank the candidate interlayers, we combined transport barriers, electrostatic potential drops, adhesion energies, and decomposition energetics into a comparative resistance metric. The bare interface is assigned a normalized resistance of 1.00. Under identical assumptions, the Li₃PO₄-coated interface falls to 0.30 ± 0.05, corresponding to an estimated reduction of about 70%. LiNbO₃ yields a normalized resistance of 0.43 ± 0.07, while Li₂ZrO₃ gives 0.58 ± 0.09. These values should be interpreted as relative rather than absolute, but the ranking is robust across reasonable variations in prefactors and defect concentrations.

![Figure 5](figures/coating_effectiveness.png) - Coating layer effectiveness

The superior performance of Li₃PO₄ arises from a favorable balance of properties. First, its wide band gap suppresses electronic leakage and reduces the chance of redox-assisted interfacial instability. Second, its phosphate chemistry is more compatible with oxygen-rich cathode surfaces and less prone to the soft-acid/soft-base interactions that destabilize sulfide|oxide contacts. Third, although Li₃PO₄ is not the fastest lithium conductor among all candidate materials, a 5-10 Å layer introduces only a moderate transport path length while substantially smoothing the local potential landscape. Fourth, its adhesion energies are high enough to maintain contact under moderate stress yet not so high as to create excessive structural distortion.

The thickness dependence is nontrivial. A 5 Å Li₃PO₄ layer already reduces chemical reactivity, but the transport and electrostatic benefits are less complete because the parent materials remain strongly coupled across the ultrathin interlayer. At 10 Å, the balance is optimal in our models. At 15 Å, the chemical protection remains excellent, but the added path length begins to offset some of the kinetic gains. This illustrates a general design principle: interlayers should be thick enough to decouple reactive surfaces but thin enough to avoid becoming the dominant ionic bottleneck.

### 5.7 Workflow validation and integrated case-study summary

The advantage of the proposed framework is that it converts multiple independent observables into a coherent mechanistic diagnosis. Figure 6 outlines the computational chain from structure selection to coating ranking, while Figure 7 consolidates the quantitative case-study outcomes.

![Figure 7](figures/case_study_summary.png) - Case study summary

The integrated view is as follows. The Li₆PS₅Cl/LiCoO₂ interface begins with moderate lattice mismatch and favorable adhesion, meaning that physical contact can be achieved experimentally. However, the same contact induces local coordination distortions that elevate lithium migration barriers from 0.24-0.31 eV in the bulk to about 0.58-0.63 eV at the bare interface. Concurrently, lithium chemical-potential mismatch establishes a 0.31-0.47 V space charge potential drop over roughly 5-15 nm, depleting mobile carriers near the contact. Finally, thermodynamic analysis shows that the direct interface is chemically unstable, with decomposition energies as favorable as -185 meV atom⁻¹. These three effects—kinetic, electrostatic, and chemical—are mutually reinforcing and together explain the high interface resistance.

The framework was also tested for internal consistency. NEB and AIMD produce comparable activation energies. Electrostatic trends correlate with the direction of DFT potential alignment. Coating rankings remain stable whether one emphasizes transport, chemistry, or combined resistance metrics. This methodological coherence is important because it suggests the framework can be transferred to other cathode|electrolyte combinations without relying on ad hoc interpretation.

## 6. Discussion

### 6.1 Comparison with experimental observations

Experimental studies of sulfide-electrolyte/oxide-cathode interfaces routinely report large interfacial impedances, especially for uncoated or insufficiently protected contacts involving LiCoO₂ and other layered oxides [5,7]. The present calculations do not seek to reproduce a single numeric impedance directly; nevertheless, the mechanisms and relative trends align closely with experimental observations. First, the prediction of thermodynamic instability for bare Li₆PS₅Cl/LiCoO₂ is consistent with reported interfacial degradation and capacity fade in bulk ASSLB architectures. Second, the calculated benefit of interlayers agrees with the longstanding experimental strategy of applying thin oxide or phosphate coatings to cathode particles. Third, the magnitude of the computed migration-barrier increase—roughly 0.3 eV relative to the bulk—is entirely sufficient to explain order-of-magnitude changes in local transport at room temperature.

A particularly important connection to experiment is that high interface resistance is rarely attributable to one factor alone. Cells with apparently good physical contact can still exhibit high impedance because chemical decomposition forms poorly conducting products. Conversely, chemically stable interfaces can remain resistive if the electrostatic alignment depletes mobile lithium carriers or if the coating itself is too thick. The present framework explicitly reproduces this multi-cause character and therefore provides a more realistic bridge between atomistic simulation and measured electrochemical behavior.

### 6.2 Origin of high interface resistance: coupled space charge and decomposition

The central conclusion of this case study is that interface resistance at the Li₆PS₅Cl/LiCoO₂ contact originates from coupled, not isolated, mechanisms. The common temptation is to identify a single dominant cause—either space charge or chemical decomposition or poor contact. Our results suggest that such a separation is artificial for this material pair.

The space charge layer alone imposes a large penalty by redistributing lithium and altering defect populations over nanometer scales. A potential drop of 0.3-0.5 V is enough to shift local site occupancies and reduce effective carrier density. However, space charge effects would be less damaging if the atomistic migration network remained highly conductive. Conversely, a structurally distorted interface with elevated barriers would still be partly tolerable if the electrostatic environment preserved abundant mobile lithium. The severe resistance emerges because both effects occur simultaneously.

Chemical decomposition then compounds the problem over longer times. The reaction-energy analysis shows that the system lowers its free energy by forming interphase products. Such products are expected to possess lower ionic conductivity than Li₆PS₅Cl and may also roughen the local potential landscape. Once formed, they can effectively widen the resistive region and lock in the transport penalties initially caused by the bare structural mismatch and space charge field. Thus, the initial contact condition sets the stage for long-term impedance growth.

This layered interpretation reconciles apparently conflicting viewpoints in the literature. Some studies emphasize space charge; others highlight chemical instability. In practice, space charge may dominate early-stage transport suppression in a freshly assembled interface, while chemical evolution determines long-term degradation and irreversible impedance growth. An integrated framework is necessary precisely because the time ordering and relative contribution of these mechanisms depend on the material pair and processing history.

### 6.3 Why Li₃PO₄ is effective

The strong performance of Li₃PO₄ in our screening can be rationalized by three complementary attributes: electronic insulation, chemical compatibility, and moderate ionic transport.

First, Li₃PO₄ has a wide band gap, which suppresses electron injection across the interlayer. This matters because electronic leakage into the solid electrolyte or into decomposition products can enable parasitic redox reactions and broaden the electrochemically active interphase. A wide-gap layer helps preserve the intended ionic selectivity of the junction.

Second, the phosphate chemistry of Li₃PO₄ is intrinsically more compatible with oxide cathode surfaces than a direct sulfide contact. Phosphate units are already oxygen-coordinated and therefore avoid the strong thermodynamic driving force associated with sulfide oxidation by high-valent oxide surfaces. At the same time, Li₃PO₄ is not so reactive toward Li₆PS₅Cl as to create a large new decomposition channel under the modeled conditions. This chemical moderation is reflected in the near-neutral reaction energies of the coated stack.

Third, Li₃PO₄ offers a transport compromise. Its ionic conductivity is lower than that of the best sulfide conductors, but when used as a nanoscale coating the relevant quantity is not bulk conductivity alone; it is the total resistance contribution of a very thin layer embedded within a strongly improved interface. By reducing the local migration barrier from ~0.60 eV to ~0.41 eV and diminishing the space charge penalty, Li₃PO₄ more than compensates for its own finite transport cost, yielding a net ~70% resistance reduction in the combined metric.

The comparison with LiNbO₃ and Li₂ZrO₃ is instructive. LiNbO₃ can be highly effective chemically and is indeed widely used experimentally, but in the present models its local transport penalty remains somewhat higher than that of Li₃PO₄. Li₂ZrO₃ is even more chemically robust yet too resistive and mechanically less optimal at the thicknesses studied. Therefore, coating selection cannot rely on chemical stability alone; one must optimize a multidimensional property set.

### 6.4 Limitations of the framework

Although comprehensive, the framework has several limitations that should inform interpretation.

First, DFT accuracy remains finite. The choice of exchange-correlation functional and Hubbard-U parameter affects redox energetics, defect levels, and reaction energies. While the adopted PBE+U setup is standard and physically justified for LiCoO₂, quantitative values such as decomposition energies may shift under alternative functionals or more advanced treatments.

Second, finite-size effects are unavoidable in explicit interface supercells. The modeled structure captures local coordination and strain but cannot represent the full disorder, roughness, grain boundaries, or mesoscale porosity of real composite cathodes. The computed barriers should therefore be viewed as local atomistic descriptors rather than literal cell-level transport coefficients.

Third, the AIMD temperatures are much higher than operating temperatures. Elevated-temperature simulations are necessary to sample sufficient lithium hops within tens of picoseconds, but extrapolation to room temperature assumes Arrhenius behavior and may miss rare-event mechanisms or subtle phase changes.

Fourth, the space charge treatment is hybrid rather than fully self-consistent across all scales. DFT provides boundary information, but the depletion widths are obtained from continuum electrostatics with assumed dielectric and defect parameters. This is a practical necessity because DFT cannot access 5-15 nm space charge regions directly, yet it introduces modeling dependence.

Fifth, the chemical stability analysis focuses on thermodynamic driving forces, not full reaction kinetics. Some exergonic decomposition channels may be kinetically inaccessible under realistic processing temperatures or short times. Conversely, metastable products not included in the phase-diagram set could form transiently and influence impedance.

Despite these limitations, the framework remains valuable because it integrates the dominant mechanisms in a transparent manner and yields robust relative trends. The main design conclusions—namely that the bare Li₆PS₅Cl/LiCoO₂ interface is both kinetically and chemically unfavorable, and that Li₃PO₄ is an effective mitigation layer—are unlikely to be reversed by reasonable methodological refinements.

### 6.5 Design guidelines for interface engineering

The present results suggest several general design rules for low-resistance interfaces in ASSLBs.

1. **Minimize direct sulfide|high-voltage oxide contact.** Even when mechanical contact is good, the combination of space charge and exergonic decomposition can create a severe transport bottleneck.

2. **Treat transport and chemistry simultaneously.** A chemically stable interface that blocks lithium transport is not useful, and a highly conductive contact that rapidly decomposes will fail during cycling. Candidate coatings should therefore be screened using both migration barriers and reaction energies.

3. **Use electrostatically intermediate layers.** Coatings that moderate lithium chemical-potential mismatch can reduce the potential drop and depletion width, thereby complementing purely chemical passivation.

4. **Optimize thickness, not just composition.** Interlayers that are too thin may not fully suppress reactive coupling, while layers that are too thick can become transport bottlenecks. For the present case, ~10 Å is near-optimal for Li₃PO₄.

5. **Incorporate adhesion and strain compatibility.** Mechanically weak or highly strained coatings may crack, delaminate, or create new defect-rich transport barriers. Positive but moderate adhesion is desirable.

6. **Link atomistic and mesoscale descriptors.** Explicit interface calculations reveal the local origin of resistance, but space charge and defect redistribution require a continuum extension. Future screening pipelines should preserve this multiscale philosophy.

These guidelines are broadly applicable beyond Li₆PS₅Cl/LiCoO₂. Similar principles should govern other sulfide|oxide and even oxide|oxide interfaces, although the relative weights of kinetic, electrostatic, and chemical contributions will differ by system.

## 7. Conclusion

This work has developed and demonstrated a first-principles computational framework for elucidating interface resistance in all-solid-state lithium-ion batteries. By integrating DFT structure optimization, CI-NEB migration-barrier analysis, AIMD finite-temperature dynamics, continuum-informed space charge modeling, chemical reaction thermodynamics, and coating-layer screening, the framework provides a unified route to diagnose and mitigate interfacial bottlenecks.

Applied to the Li₆PS₅Cl/LiCoO₂ system, the framework shows that the bare interface is problematic for three reinforcing reasons. First, local structural distortion caused by lattice mismatch and interfacial bonding raises lithium migration barriers from 0.24-0.31 eV in bulk Li₆PS₅Cl to approximately 0.58-0.63 eV at the interface. Second, lithium chemical-potential equilibration generates a space charge layer with a potential drop of 0.31-0.47 V and a characteristic width of 5-15 nm, depleting mobile carriers near the contact. Third, the direct interface is thermodynamically unstable, with decomposition energies reaching approximately -185 meV atom⁻¹ for representative pathways. These findings explain why interface resistance remains a primary bottleneck even when the bulk electrolyte is highly conductive.

Among the coatings examined, Li₃PO₄ provides the best overall mitigation. It lowers the effective transport barrier to about 0.39-0.44 eV, moderates the space charge potential, suppresses direct chemical decomposition, and reduces the estimated interface resistance by roughly 70% relative to the bare contact. The case study therefore highlights the importance of balancing chemical stability, electrostatic alignment, ionic transport, and adhesion when designing interlayers.

More broadly, the framework offers practical guidance for ASSLB development. It can be extended to alternative sulfide, oxide, halide, or polyanionic electrolytes; to Ni-rich layered cathodes and conversion cathodes; and to interface architectures involving gradient coatings or composite interphases. Future work should couple this methodology with machine-learning-accelerated materials screening, explicit defect-disorder sampling, and operando experimental validation through spectroscopy and microscopy. Such an integrated computational-experimental strategy will be essential for transforming solid-state battery interfaces from empirical bottlenecks into rationally engineered functional junctions.

## References

1. Xiao, Y., Wang, Y., Bo, S.-H., Kim, J. C., Miara, L. J., & Ceder, G. (2020). Understanding interface stability in solid-state batteries. *Nature Reviews Materials, 5*(2), 105-126. https://doi.org/10.1038/s41578-019-0157-5

2. Richards, W. D., Miara, L. J., Wang, Y., Kim, J. C., & Ceder, G. (2016). Interface stability in solid-state batteries. *Chemistry of Materials, 28*(1), 266-273. https://doi.org/10.1021/acs.chemmater.5b04082

3. Haruyama, J., Sodeyama, K., Han, L., Takada, K., & Tateyama, Y. (2014). Space-charge layer effect at interface between oxide cathode and sulfide electrolyte in all-solid-state lithium-ion battery. *Chemistry of Materials, 26*(14), 4248-4255. https://doi.org/10.1021/cm5016959

4. Zhu, Y., He, X., & Mo, Y. (2015). Origin of outstanding stability in the lithium solid electrolyte materials: insights from thermodynamic analyses based on first-principles calculations. *ACS Applied Materials & Interfaces, 7*(42), 23685-23693. https://doi.org/10.1021/acsami.5b07517

5. Auvergniot, S., Viallet, V., Seznec, V., Boulineau, A., Sougrati, M. T., Reynier, Y., ... & Masquelier, C. (2017). Interface stability of argyrodite Li6PS5Cl toward LiCoO2, LiNi1/3Co1/3Mn1/3O2, and LiMn2O4 in bulk all-solid-state batteries. *Chemistry of Materials, 29*(9), 3883-3890. https://doi.org/10.1021/acs.chemmater.6b04990

6. Schwietert, T. K., Arszelewska, V. A., Wang, C., Yu, C., Vasileiadis, A., de Klerk, N. J. J., ... & Wagemaker, M. (2020). Clarifying the relationship between redox activity and electrochemical stability in solid electrolytes. *Nature Materials, 19*(4), 428-435. https://doi.org/10.1038/s41563-019-0576-0

7. Ohta, N., Takada, K., Zhang, L., Ma, R., Osada, M., & Sasaki, T. (2006). Enhancement of the high-rate capability of solid-state lithium batteries by nanoscale interfacial modification. *Advanced Materials, 18*(17), 2226-2229. https://doi.org/10.1002/adma.200502604

8. Lacivita, V., Arber, A. S., Brunell, I. F., Cargiulo, B., & Ceder, G. (2024). Space charge layers and interface potentials in solid-state batteries. *Physical Review Materials, 8*(10), 105402. https://doi.org/10.1103/PhysRevMaterials.8.105402

9. Muy, S., Bachman, J. C., Giordano, L., Chang, H. H., Abernathy, D. L., Bansal, D., ... & Shao-Horn, Y. (2018). Tuning mobility and stability of lithium ion conductors based on lattice dynamics. *Energy & Environmental Science, 11*(4), 850-859. https://doi.org/10.1039/C7EE03364H

10. Nolan, A. M., Zhu, Y., He, X., Bai, Q., & Mo, Y. (2018). Computation-accelerated design of materials and interfaces for all-solid-state lithium-ion batteries. *Joule, 2*(10), 2016-2046. https://doi.org/10.1016/j.joule.2018.08.017
