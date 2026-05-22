# Statistical Summary

## Assumption checks

- Reuse-count skewness: 0.813
- Shapiro-Wilk p-value on log(1 + reuse): 0.000001
- Recommended test family: nonparametric
- Rationale: Reuse counts are right-skewed count data, so rank-based inference is preferred.

## Documentation and reuse

- Spearman correlation: 0.318
- 95% CI: [0.272, 0.362]
- p-value: 1.584149e-36
- Pearson correlation (reference): 0.323, p = 1.208189e-37

## License effects on reuse

- Kruskal-Wallis H: 23.299
- p-value: 3.498106e-05
- Effect size (epsilon squared): 0.014

| Comparison | Rank-biserial r | Mean reuse diff. | 95% CI | Raw p | FDR-adjusted p |
|---|---:|---:|---:|---:|---:|
| CC-BY vs restricted | 0.149 | 1.204 | [0.634, 1.796] | 0.000473 | 0.000709 |
| CC0 vs restricted | 0.220 | 1.799 | [1.072, 2.532] | 0.000010 | 0.000030 |
| custom vs restricted | 0.073 | 0.626 | [-0.059, 1.279] | 0.114506 | 0.114506 |

## Interpretation notes

- Effect sizes and confidence intervals are reported alongside p-values.
- False-discovery-rate adjustment was applied to the three pairwise license comparisons.
- Practical significance is modest-to-moderate: open licenses improve reuse, but documentation and FAIR quality remain important co-drivers.
