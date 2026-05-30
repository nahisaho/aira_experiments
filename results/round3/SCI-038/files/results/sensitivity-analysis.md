# Sensitivity Analysis

## Seed Variation
The mission-sequencing experiment was repeated for six independent catalog seeds.

| Seed | NN Δv (km/s) | GA+2opt Δv (km/s) | Improvement (%) |
| --- | --- | --- | --- |
| 7 | 13.075 | 12.933 | 1.08 |
| 11 | 21.303 | 21.521 | -1.02 |
| 19 | 10.832 | 10.098 | 6.78 |
| 23 | 10.872 | 10.872 | 0.00 |
| 31 | 19.962 | 20.570 | -3.05 |
| 42 | 12.963 | 12.848 | 0.89 |

Summary: nearest-neighbor averaged 14.834 ± 4.614 km/s, while GA+2opt averaged 14.807 ± 4.966 km/s. Mean improvement was 0.780 ± 3.304%, with 95% CI [-2.687, 4.248]. This indicates that the optimizer is reasonably stable, but the advantage of the evolutionary layer is not uniform across all random realizations.

## Hyperparameter Perturbation (±10% thrust)

| Thrust (N) | Final Δv (km/s) | Transfer time (h) | Fuel used (kg) | Final inclination (deg) |
| --- | --- | --- | --- | --- |
| 0.18 | 0.299 | 462.1 | 10.13 | 69.82 |
| 0.20 | 0.318 | 441.9 | 10.76 | 70.97 |
| 0.22 | 0.318 | 401.7 | 10.76 | 70.97 |

Summary: the low-thrust transfer remained well behaved under moderate thrust perturbations. Higher thrust reduced time of flight while preserving near-target inclination convergence.

## Data-Size Perturbation

| Catalog size | Top targets | GA+2opt Δv (km/s) |
| --- | --- | --- |
| 16 | 10 | 10.166 |
| 20 | 10 | 12.848 |
| 24 | 10 | 11.974 |

Summary: sequence cost changed with catalog size, but not monotonically. The dominant factor is the geometry of the highest-priority subset rather than the total background population alone.
