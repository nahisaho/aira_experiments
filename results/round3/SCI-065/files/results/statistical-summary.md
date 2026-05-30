# Statistical Summary

All statistics are descriptive because they arise from synthetic replicate observations and parameter sweeps generated from the mechanistic framework.

## Core metrics with uncertainty

| Metric | Mean ± SD | 95% CI |
|---|---|---|
| Viability at 2.5 mm | 0.257 ± 0.015 | [0.238, 0.276] |
| Day-90 maturation | 0.520 ± 0.020 | [0.496, 0.545] |
| Day-90 MMI | 0.744 ± 0.030 | [0.706, 0.781] |

## Pairwise comparisons

| Comparison | Group A mean ± SD | Group B mean ± SD | Welch t | p-value | Cohen d |
|---|---:|---:|---:|---:|---:|
| Viability: 0.5 mm vs 2.5 mm organoids | 0.986 ± 0.021 | 0.257 ± 0.015 | 63.66 | 2.86e-10 | 40.26 |
| Maturation: 0.01 Pa vs optimal shear at day 90 | 0.430 ± 0.015 | 0.597 ± 0.017 | -16.27 | 2.10e-07 | -10.29 |
| Cost at scales 100-1000: Batch vs Perfusion | 11.953 ± 4.278 | 12.462 ± 4.314 | -0.15 | 8.91e-01 | -0.12 |
| Cost at scales 100-1000: Perfusion vs Continuous | 12.462 ± 4.314 | 16.041 ± 6.770 | -0.77 | 4.91e-01 | -0.63 |

## Assumption checks

- Normality was not assumed for definitive inference because deterministic cores were combined with stochastic observation noise.
- Welch tests were used to reduce sensitivity to unequal variances.
- Because these are simulated replicates rather than independent biological cohorts, effect sizes and confidence intervals are more informative than nominal p-values.
