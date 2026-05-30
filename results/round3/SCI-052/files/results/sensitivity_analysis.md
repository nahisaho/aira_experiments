# Sensitivity Analysis

Five random seeds were used to perturb all forward FT barriers by Gaussian noise with a standard deviation of 0.015 eV. The resulting TOF remained within the expected FT range.

Mean TOF with lateral interactions = 2.1383e-03 ± 6.56e-04 s^-1 per site; 95% CI [1.3235e-03, 2.9530e-03].

Paired t-test versus the no-interaction baseline: t = -8.055, p = 1.290e-03, Cohen's d = -3.602.

The sign of the paired difference was negative for all seeds, indicating that attractive CO self-interactions reduce TOF by increasing surface poisoning while raising CO coverage.


## Hyperparameter Perturbation

The CO adsorption self-interaction parameter was perturbed by ±20% around ω_CO = -0.30 eV. At 500 K, the TOF changed from 5.1362e-03 to 5.9583e-03 s^-1 in the surrogate lateral-interaction scan shown in Figure 5.
