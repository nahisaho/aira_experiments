# Sensitivity Analysis

## Seed variation

| Metric | Mean | SD | 95% CI |
|---|---:|---:|---:|
| best_hydrolysis_r2 | 0.916 | 0.016 | [0.896, 0.937] |
| pareto_count | 28.600 | 8.081 | [18.566, 38.634] |
| top_composition_score | 0.799 | 0.001 | [0.797, 0.801] |
| surface_pha_365_mean | 100.000 | 0.000 | [100.000, 100.000] |

## Hyperparameter perturbation (SVR C ±20%)

| Setting | Mean R² | SD | 95% CI |
|---|---:|---:|---:|
| SVR_C=0.4 | 0.927 | 0.041 | [0.876, 0.977] |
| SVR_C=0.5 | 0.929 | 0.041 | [0.878, 0.979] |
| SVR_C=0.6 | 0.929 | 0.040 | [0.879, 0.980] |

Hydrolysis performance remained stable across seeds and under ±20% variation in SVR regularization, indicating that the main ranking conclusions are not driven by a fragile configuration.