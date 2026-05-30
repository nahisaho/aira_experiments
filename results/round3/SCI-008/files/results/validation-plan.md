# Validation Plan

Internal validation: 3-fold cross-validation with identical fold sizes and seed offsets (42, 43, 44).
Primary metrics: MRR and Hits@10.
Secondary metrics: Hits@1, Hits@3, AUROC.
Model selection criterion: best validation MRR within each fold.
Uncertainty reporting: mean ± standard deviation across folds.
Failure criterion: any model with MRR > 0.60 on this noisy synthetic graph is considered suspicious and triggers inspection for leakage.
