# Statistical Summary

## Hydrolysis model
- test_rmse: 0.00002
- test_r2: 0.93662
- cv_rmse_mean: 0.00002
- cv_r2_mean: 0.76632

## ML comparison
- GradientBoosting: RMSE=0.00001, R²=0.9459, 95% abs. error interval=[0.00000, 0.00002]
- RandomForest: RMSE=0.00001, R²=0.9003, 95% abs. error interval=[0.00000, 0.00004]
- LinearRegression: RMSE=0.00003, R²=0.5467, 95% abs. error interval=[0.00000, 0.00010]
- MLP: RMSE=0.14276, R²=-9521082.8196, 95% abs. error interval=[0.00320, 0.33073]

## Marine degradation summary
| polymer | final_weight_loss_pct | mean_k_total |
| --- | --- | --- |
| PBS | 4.11312 | 0.00011 |
| PHA | 8.14658 | 0.00023 |
| PLA | 5.03040 | 0.00014 |

## Top case-study variants
| polymer | variant | strategy_note | final_weight_loss_pct | apparent_k_h |
| --- | --- | --- | --- | --- |
| PHA | PHBV blend | blend + compatibilizer | 0.05546 | 0.00000 |
| PHA | PHBV 20% HV | HV lowers crystallinity | 0.04926 | 0.00000 |
| PLA | High-D stereodefect | higher D-unit content | 0.04671 | 0.00000 |
| PBS | Branched PBS | branching increases water uptake | 0.04452 | 0.00000 |
| PBS | Standard | baseline | 0.03649 | 0.00000 |

## Representative Pareto points
| source | crystallinity_pct | Mw | k_h | tensile_strength_MPa | modulus_GPa |
| --- | --- | --- | --- | --- | --- |
| grid_screen | 65.00000 | 300000.00000 | 0.00000 | 803.27622 | 366.10857 |

## Top copolymer candidates
| system | x1 | x2 | x3 | score | k_h_copol | strength_MPa | modulus_GPa |
| --- | --- | --- | --- | --- | --- | --- | --- |
| LA-GA-ε-CL | 0.00000 | 1.00000 | 0.00000 | 0.85000 | 0.01050 | 72.00000 | 2.80000 |
| LA-GA-ε-CL | 0.05000 | 0.95000 | 0.00000 | 0.83913 | 0.01047 | 70.92000 | 2.75100 |
| LA-GA-ε-CL | 0.10000 | 0.90000 | 0.00000 | 0.82812 | 0.01042 | 69.88000 | 2.70400 |
| LA-GA-ε-CL | 0.00000 | 0.95000 | 0.05000 | 0.81843 | 0.01020 | 69.62000 | 2.66350 |
| LA-GA-ε-CL | 0.15000 | 0.85000 | 0.00000 | 0.81697 | 0.01033 | 68.88000 | 2.65900 |

Note: Multi-model comparison involves 4 models; rankings were compared descriptively, and confidence intervals for absolute error are reported for transparency.