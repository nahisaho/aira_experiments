# Ablation Results

The composite NCC framework was benchmarked against reduced feature variants using logistic regression with five-fold stratified cross-validation.

- PCI only: macro-AUC = 0.690 ± 0.136 (95% CI ± 0.169)
- Phi + PCI: macro-AUC = 0.800 ± 0.106 (95% CI ± 0.131)
- Phi + PCI + GWT: macro-AUC = 0.932 ± 0.044 (95% CI ± 0.055)
- Full (Phi + PCI + GWT + Spectral Entropy): macro-AUC = 0.936 ± 0.045 (95% CI ± 0.055)
