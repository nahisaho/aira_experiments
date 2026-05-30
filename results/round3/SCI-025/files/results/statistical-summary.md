# Statistical Summary

## Hydrolysis model comparison

| Model | Mean R² | SD | 95% CI |
|---|---:|---:|---:|
| SVR | 0.929 | 0.041 | [0.878, 0.979] |
| LinearRegression | 0.925 | 0.024 | [0.895, 0.955] |
| RandomForest | 0.877 | 0.036 | [0.832, 0.921] |
| GradientBoosting | 0.873 | 0.027 | [0.839, 0.907] |

## Pairwise comparison against best model

| Comparison | Mean ΔR² | 95% CI | p | Holm p | Effect size dz |
|---|---:|---:|---:|---:|---:|
| SVR vs LinearRegression | 0.004 | [-0.061, 0.069] | 0.8640 | 0.8640 | 0.082 |
| SVR vs RandomForest | 0.052 | [-0.018, 0.123] | 0.1086 | 0.2173 | 0.921 |
| SVR vs GradientBoosting | 0.056 | [-0.014, 0.126] | 0.0912 | 0.2737 | 0.990 |

Surface-zone PHA reached 100.000000% mass loss at day 365 with 95% CI [99.999997, 100.000000].
