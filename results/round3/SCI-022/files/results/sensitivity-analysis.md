# Sensitivity Analysis

Seed variation and ±20% hyperparameter perturbations were evaluated for the bandgap regressor. Across 45 runs, the average cross-validated R² was 0.950 ± 0.003, with a 95% CI of 0.947 to 0.953 across seed-level means. The corresponding MAE was 0.065 ± 0.002 eV and RMSE was 0.088 ± 0.003 eV. The best configuration used seed 7, learning rate 0.04, and 100 estimators (R²=0.956), whereas the weakest tested configuration still retained R²=0.944.

## Selected runs

- seed=7, lr=0.04, estimators=100: R²=0.956 ± 0.043, MAE=0.061 ± 0.019 eV
- seed=7, lr=0.04, estimators=120: R²=0.956 ± 0.040, MAE=0.062 ± 0.017 eV
- seed=7, lr=0.04, estimators=140: R²=0.955 ± 0.039, MAE=0.063 ± 0.016 eV
- seed=84, lr=0.06, estimators=100: R²=0.953 ± 0.037, MAE=0.063 ± 0.016 eV
- seed=84, lr=0.06, estimators=120: R²=0.953 ± 0.036, MAE=0.064 ± 0.016 eV
- seed=84, lr=0.05, estimators=100: R²=0.953 ± 0.036, MAE=0.064 ± 0.015 eV
- seed=84, lr=0.05, estimators=120: R²=0.953 ± 0.035, MAE=0.064 ± 0.015 eV
- seed=84, lr=0.06, estimators=140: R²=0.952 ± 0.035, MAE=0.065 ± 0.015 eV
- seed=7, lr=0.06, estimators=100: R²=0.952 ± 0.038, MAE=0.065 ± 0.015 eV
- seed=126, lr=0.06, estimators=100: R²=0.952 ± 0.039, MAE=0.065 ± 0.015 eV
