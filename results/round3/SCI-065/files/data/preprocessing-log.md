# Preprocessing Log

No external biological dataset was ingested for this task. The simulation framework uses mechanistic and semi-empirical parameterizations directly encoded in the source modules. Preprocessing therefore consisted of parameter normalization, unit harmonization (m, mm, Pa, mol/m³, days), and reproducible seeding with `numpy.random.seed(42)`.

Derived quantities include radial and axial finite-difference grids, normalized media factor schedules, and standardized biomarker trajectories for PCA-derived maturation weighting. No raw subject-level data, patient identifiers, or wet-lab files were processed.
