# Power Analysis

This study is computational rather than wet-lab experimental, so the dominant source of uncertainty is variability across cross-validation folds rather than biological replicates. We therefore use fold-level effect estimation to justify a lightweight but informative design.

## Assumptions

- Primary endpoint: mean reciprocal rank (MRR)
- Minimal practically relevant effect size between embedding models: 0.08 MRR
- Expected fold-level standard deviation: 0.05 to 0.07 based on recent biomedical KGE benchmarking studies
- Significance level: α = 0.05
- Target power: 0.80

Using a paired standardized effect size of roughly d = 1.14 (0.08 / 0.07) suggests that 3 to 5 paired folds are sufficient for a preliminary computational comparison. Given the runtime constraint and synthetic-data setting, we use 3 folds and report means ± standard deviations, while treating inferential claims as exploratory rather than confirmatory.

## Interpretation

The design is adequately powered to distinguish clearly separated model families, but not to claim subtle superiority between methods whose fold-level MRR differs by less than ~0.04. Therefore the report emphasizes effect sizes, uncertainty bands, and qualitative ranking stability rather than binary significance claims.
