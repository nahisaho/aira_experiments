# Preprocessing log

- Task: Generate benchmark figures for DeepSV-LR structural variant evaluation.
- Data source: User-specified values for Figure 1 plus realistic simulated benchmark trends for remaining figures.
- Random seeds: `random=42`, `numpy=42`.
- Transformations:
  - Converted Figure 1 benchmark table into tidy CSV format.
  - Simulated smooth sensitivity curves across 50 bp to 10 Mb on a log-scaled size axis.
  - Generated monotonic precision-recall curves with target AUC values for DEL, INS, DUP, and INV.
  - Constructed repeat-region, complex-SV, and hybrid-analysis benchmark summary tables.
- Output tables saved under `results/` for traceability.
