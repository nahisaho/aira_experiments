# Statistical Summary

Assumptions: fold-level metrics are approximately symmetric because each fold contains a balanced mixture of SV classes; confidence intervals are reported from the t approximation across five folds. Multiple comparisons across methods were interpreted conservatively and effect sizes are emphasized over nominal p-values.

- split_read: F1 = 0.799 ± 0.032; 95% CI [0.771, 0.827]
- read_depth: F1 = 0.610 ± 0.030; 95% CI [0.583, 0.636]
- assembly: F1 = 0.769 ± 0.022; 95% CI [0.750, 0.789]
- hybrid: F1 = 0.784 ± 0.011; 95% CI [0.775, 0.793]

- Hybrid vs split_read: mean ΔF1 = -0.015, 95% CI [-0.042, 0.012], paired Cohen's d = -0.494.
- Hybrid vs read_depth: mean ΔF1 = 0.174, 95% CI [0.140, 0.209], paired Cohen's d = 4.416.
- Hybrid vs assembly: mean ΔF1 = 0.014, 95% CI [0.000, 0.029], paired Cohen's d = 0.882.
- Chromothripsis probability (synthetic positive sample): 0.925; ecDNA probability: 0.853; long-/short-read concordance F1: 0.571.