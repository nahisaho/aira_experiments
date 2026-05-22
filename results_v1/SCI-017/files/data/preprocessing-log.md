# Preprocessing Log

Generated: 2026-05-22T13:29:42.205916Z

## Input Data
- Spike peptide: 142 amino acids (Wuhan reference)
- Human codon usage: Kazusa database (Homo sapiens)
- HLA frequencies: Global population survey

## Random Seeds
- numpy: 42
- All modules: seed=42

## Transformations
- Codon table: fraction-based weights
- CAI: log-ratio geometric mean
- Parker hydrophilicity: 7-mer sliding window
- LNP metrics: physics-inspired parametric model
