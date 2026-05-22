# Preprocessing log

- Random seeds fixed to 9303 in all generated Python scripts.
- Medium optimization used bounded, phase-specific variables for bFGF, EGF, BDNF, retinoic acid, Matrigel, glucose, and oxygen tension.
- Scalability analysis converted reactor configuration assumptions into deterministic transport and cost metrics, then used CSV-safe quoting for fields containing commas.
- Biomarker monitoring data were generated on a 0.5 day grid with weekly offline markers, daily at-line assays, continuous online sensors, and an injected anomaly window at day 62-68.
- No raw external data were ingested; all datasets in this task are synthetic, processed outputs.

- Baseline transport boundary conditions were set to 0.20 mol m^-3 oxygen and 5.0 mol m^-3 glucose.
- Oxygen and glucose steady-state profiles were solved on synthetic radial grids; lactate was computed from glucose-consumption-derived source terms.
- Shear-response curves were sampled on a log-spaced 0.001-1.0 Pa grid, and the maturation heatmap was generated over a synthetic 0-60 day culture timeline.
