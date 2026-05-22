# Preprocessing Log

- Dataset type: Synthetic spatial disease-risk grid
- Coordinate system: Cartesian grid coordinates (unit spacing)
- Grid size: 12 x 12
- Random seeds: numpy=42, random=42
- Spatial autocorrelation generation: Gaussian-smoothed latent field plus mild linear trend
- Disease process: Poisson cases simulated from expected counts multiplied by latent relative risk
- Analysis variable: log-transformed risk ratio using continuity correction
- Missing data: none
- Spatial weights: distance band approximating Queen adjacency with threshold sqrt(2)+0.05
- Local inference correction: Benjamini-Hochberg FDR
- Variogram assumption: isotropic semivariance model fitted by weighted least squares

## Malaria and dengue risk mapping preprocessing

- Timestamp: 2026-05-22T17:39:08.459451+00:00
- Dataset type: Synthetic monthly malaria and dengue surveillance panel
- Grid size: 14 x 14 administrative units (196 total)
- Observation horizon: 24 months
- Random seeds: numpy=20250314, random=20250314
- Spatial structure: 4-neighbor lattice with ICAR-compatible precision matrix for BYM modeling
- Fixed covariates: elevation, population density, urbanization index
- Time-varying covariates: temperature and precipitation
- Disease simulation model: Poisson counts with disease-specific seasonal signals and correlated latent spatial risk fields
- Modeling offset: population-based expected counts used for SMR and BYM Poisson mapping
- Standardization: log(population density) and all model covariates standardized to zero mean and unit variance before regression fitting
- Missing data: none introduced
- Output data tables: `data/synthetic_disease_data.csv`, `data/synthetic_disease_area_summary.csv`



## Ecological bias simulation (2026-05-22T17:39:55.934119+00:00)
- Data type: Synthetic individual-to-area ecological confounding scenario
- Random seeds: numpy=42, random=42
- Units: 7,576 individuals nested in 64 areas on a jittered 8 x 8 lattice
- Exposure model: Bernoulli treatment driven by individual confounder, area confounder, and spatial field
- Outcome model: Continuous Gaussian outcome with true treatment effect 1.0 plus individual and area-level confounding
- Aggregation: Individual outcomes and treatments summarized to area means/prevalences for ecological regression
- Overlap handling: propensity-score stratification trimmed observations outside the 0.05 to 0.95 propensity range

## Spatiotemporal spline simulation (2026-05-22T17:39:55.934119+00:00)
- Data type: Synthetic monthly disease counts across irregular spatial locations
- Random seeds: numpy=123, random=123
- Spatial units: 90 locations in the unit square
- Temporal span: 60 observed months plus 12 forecast months
- Spatial smoother: thin-plate radial basis with knot counts evaluated by grouped cross-validation
- Temporal smoother: cubic B-spline basis via scikit-learn `SplineTransformer`
- Outcome scale: log incidence using continuity-corrected counts divided by location-specific population
- Validation strategy: GroupKFold by location to reduce optimistic spatial leakage during knot selection
