# Statistical Summary
Generated: 2026-05-22T14:20:23Z

## 1. Reaction Network
- Total species: 41
- Total reactions: 45
- Primary VOCs: 5
- Max generation: 2

## 2. Gas-Particle Partitioning
- Species analyzed: 34
- ELVOC (log C*=-3): pinic acid (Fpart≈1.0), norpinic acid (Fpart≈1.0)
- LVOC (log C*=-1): pinonic acid (Fpart≈0.99)
- SVOC (log C*=+1): pinaldehyde (Fpart≈0.53)
- IVOC (log C*=+3): methacrolein (Fpart≈0.001)

## 3. ML Rate Constant Model
- Training samples: 20
- R²: 0.9966
- RMSE: 0.0337 log units
- MAE: 0.0263 log units
- CV R² (5-fold): -10.1097 ± 17.5169
- Top features: BDE, delta_H_rxn, IP, n_double_bonds

## 4. Box Model Simulation
- alpha_pinene: final SOA = 6.664 μg/m³ (8h)
- beta_pinene: final SOA = 4.876 μg/m³ (8h)
- limonene: final SOA = 3.149 μg/m³ (8h)
- isoprene: final SOA = 18.042 μg/m³ (8h)
- toluene: final SOA = 10.364 μg/m³ (8h)

## 5. Sensitivity Analysis
### OAT (normalized sensitivity index):
- T: -51.3556
- VOC_ppb: 1.0833
- RH: -0.0002
- NOx_ppb: 0.0002
- JNO2: 0.0000
- O3_ppb: -0.0000

### Sobol First-Order Indices:
- T: S1 = 0.6696
- VOC_ppb: S1 = 0.1685
- RH: S1 = 0.0000
- NOx_ppb: S1 = 0.0000
- O3_ppb: S1 = 0.0000
- JNO2: S1 = 0.0000

## 6. SOA Yield Predictions
| VOC | Oxidant | Y_predicted | Y_literature | Uncertainty |
|-----|---------|-------------|--------------|-------------|
| alpha_pinene | OH | 0.185 | 0.300 | ±0.40 |
| alpha_pinene | O3 | 0.228 | 0.400 | ±0.39 |
| beta_pinene | OH | 0.122 | 0.150 | ±0.34 |
| beta_pinene | O3 | 0.123 | 0.130 | ±0.41 |
| limonene | OH | 0.258 | 0.390 | ±0.38 |
| limonene | O3 | 0.280 | 0.500 | ±0.49 |
| isoprene | OH | 0.028 | 0.030 | ±0.07 |
| isoprene | O3 | 0.005 | 0.010 | ±0.71 |
| toluene | OH | 0.187 | 0.280 | ±0.39 |
| toluene | O3 | 0.000 | N/A | ±0.00 |
