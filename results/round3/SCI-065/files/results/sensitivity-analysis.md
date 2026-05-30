# Sensitivity Analysis

## Seed variation (5 seeds)

| Seed | Re | Max shear (Pa) | Viability at 2.5 mm | Day-90 maturation | Media score | MMI day 90 |
|---:|---:|---:|---:|---:|---:|---:|
| 11 | 0.637 | 1.097743e-04 | 0.276 | 0.510 | 39.32 | 0.730 |
| 19 | 0.628 | 1.121969e-04 | 0.241 | 0.517 | 40.37 | 0.741 |
| 23 | 0.662 | 1.069568e-04 | 0.239 | 0.511 | 40.64 | 0.697 |
| 41 | 0.633 | 1.057016e-04 | 0.248 | 0.542 | 38.54 | 0.755 |
| 57 | 0.646 | 1.102390e-04 | 0.245 | 0.521 | 39.90 | 0.799 |

## Hyperparameter perturbation (±15%)

### CFD flow-rate sensitivity

| Flow-rate factor | Max shear (Pa) | Reynolds number |
|---:|---:|---:|
| 0.850 | 8.742e-05 | 0.541 |
| 0.925 | 9.514e-05 | 0.589 |
| 1.000 | 1.029e-04 | 0.637 |
| 1.075 | 1.106e-04 | 0.684 |
| 1.150 | 1.183e-04 | 0.732 |

### Oxygen-consumption sensitivity for 2.5 mm organoid

| Qmax factor | Necrotic core (mm) | Viability fraction |
|---:|---:|---:|
| 0.850 | 2.229 | 0.292 |
| 0.925 | 2.249 | 0.273 |
| 1.000 | 2.259 | 0.257 |
| 1.075 | 2.279 | 0.238 |
| 1.150 | 2.299 | 0.219 |

### Shear-optimum sensitivity

| tau_opt factor | Optimal shear (Pa) | Day-90 marker mean |
|---:|---:|---:|
| 0.850 | 0.0953 | 0.521 |
| 0.925 | 0.1034 | 0.521 |
| 1.000 | 0.1123 | 0.521 |
| 1.075 | 0.1218 | 0.521 |
| 1.150 | 0.1321 | 0.521 |

### Media-exchange sensitivity

| Exchange interval (days) | Integrated score |
|---:|---:|
| 2.4 | 38.86 |
| 2.7 | 39.42 |
| 3.0 | 39.80 |
| 3.3 | 39.36 |
| 3.6 | 38.73 |

## Summary

Across five stochastic observation seeds, the framework remained stable, with day-90 maturation centered near 0.52 and day-90 MMI centered near 0.75. The largest sensitivity was associated with oxygen-consumption perturbation in the 2.5 mm organoid, confirming that metabolic demand assumptions dominate viability predictions for large constructs. Flow-rate perturbations changed Reynolds number and shear proportionally, while media exchange showed a broad optimum near 3 days, supporting operational robustness.
