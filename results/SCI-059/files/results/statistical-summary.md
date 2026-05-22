# Statistical Summary

## Core metrics
- Total simulated cases: 352
- LGCP optimization success: True
- LGCP RMSE: 153.521
- KDE RMSE: 200.072
- LGCP MAE: 101.390
- KDE MAE: 119.595

## Spatial summary statistics
- First K-function value: 0.000955
- Final K-function value: 0.450903
- First Poisson reference K value: 0.000428
- Final Poisson reference K value: 0.384845
- Mean pair correlation: 1.385207

## Notes
No parametric hypothesis tests were required for this simulation task, so p-values, effect sizes, confidence intervals, and multiple-testing correction are not applicable here.

## Malaria and Dengue Risk Mapping Summary

Timestamp: 2026-05-22T17:39:08.459451+00:00

### Core metrics
- Study units: 196 administrative areas across 24 months
- Malaria total cases: 21,683
- Dengue total cases: 24,376
- Malaria posterior RR range: 0.259 to 4.576
- Dengue posterior RR range: 0.084 to 5.479
- High-risk overlap at exceedance > 0.80: 13 areas

### Posterior covariate effects (RR multipliers, 95% credible intervals)
- Malaria precipitation: 1.37 (1.12 to 1.68)
- Malaria elevation: 0.68 (0.42 to 1.10)
- Dengue urbanization: 1.46 (1.11 to 1.92)
- Dengue population density: 1.14 (0.82 to 1.58)

### Notes
- This workflow uses Bayesian posterior summaries and 95% credible intervals instead of frequentist p-values.
- Multiple-testing correction is not applicable because no family of frequentist hypothesis tests was performed.
- Spatial uncertainty is propagated through the Laplace-approximated posterior covariance matrix.



## Ecological bias metrics
- True treatment effect: 1.000
- Naive ecological estimate: 9.511 (95% CI 8.918 to 10.104)
- Stratified estimate: 1.117 (95% CI 1.041 to 1.202)
- Multilevel estimate: 1.004 (95% CI 0.934 to 1.074)
- Spatial basis estimate: 1.663 (95% CI 0.749 to 2.578)
- Propensity-score spatial matching estimate: 0.925 (95% CI 0.871 to 0.973)
- Absolute ecological bias: 8.511

## Spatiotemporal spline metrics
- Selected spatial knot count: 8
- Selected ridge alpha: 10.0
- Training RMSE: 0.1900
- Mean future log incidence: -2.5615
- Peak forecast month index: 71

## Notes
- Ecological correction methods report effect sizes with confidence intervals; bootstrap intervals were used for stratification and spatial matching.
- Spatiotemporal uncertainty bands are approximate model-based intervals derived from the penalized regression covariance.
- No multi-test adjustment was necessary because the deliverables emphasize estimation rather than a family of null-hypothesis tests.
