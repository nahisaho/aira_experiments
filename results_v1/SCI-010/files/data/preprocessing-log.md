# Preprocessing Log

- Random seed fixed at 42 for numpy and random.
- All outputs were generated from deterministic simulation settings unless Monte Carlo sampling was explicitly requested.
- DAR sampling used binomial and Poisson models with support truncated to DAR 0-8.
- Linker kinetics were normalized to released fraction for cross-mechanism comparison.
- Reaction-diffusion PDE used a 1D finite-difference explicit solver with no-flux boundaries.
- Optimization landscape stored complete grid in `data/optimization_landscape.csv`; Pareto front stored in `results/optimization_results.csv`.
- PK/PD nominal dose assumed 6.4 mg/kg IV for a 70 kg patient.
- Monte Carlo sensitivity used Latin Hypercube Sampling over ±30% parameter ranges.
- Case study metrics were normalized only for visualization; raw values remain in `results/case_study_summary.csv`.
