# DRAFT — NOT FOR DISTRIBUTION

# Statistical summary

This summary reports effect sizes and 95% confidence intervals for synthetic design analyses used to benchmark the requested deliverables. Bonferroni-adjusted p-values were applied across three comparisons.

| Comparison | Effect size | 95% CI | Adjusted p-value | Interpretation |
|---|---:|---:|---:|---|
| Optimized medium score vs. random feasible programs | Cohen's d = 7.55 | Δscore [60.29, 92.96] | 0.005988 | Optimization strongly outperformed naive random schedules. |
| Continuous vs. batch cost per organoid | Cohen's d = 10.52 | Cost saving [3.61, 6.17] USD/organoid | 0 | Continuous operation reduced modeled unit cost while supporting large yield gains. |
| Anomaly window vs. baseline anomaly score | Cohen's d = 30.64 | Δanomaly score [23.69, 24.78] | 5.731e-23 | Control-chart and anomaly-monitoring logic clearly separated deviating runs from baseline. |

## Additional deterministic comparisons
- Continuous scale-up factor vs. batch: 303.3x
- Continuous yield minus batch yield: 33011 organoids/batch (sensitivity interval [29230, 36879])
- Perfusion cost per organoid: $2.91
- Monitoring anomaly flags triggered at 184 of 241 sampled time points

## Assumptions and limitations
- These comparisons are synthetic and intended for design prioritization, not biological validation.
- Confidence intervals for scalability come from Monte Carlo perturbation of deterministic model parameters rather than replicate bioreactor runs.
- Online monitoring assumes stable baseline calibration during the first 20 culture days.


## Deterministic model summary: transport and shear-maturation

No hypothesis tests were applied to these new deliverables because both scripts are deterministic mechanistic/phenomenological models without replicate sampling. Practical design metrics are therefore emphasized instead of p-values.

| Deliverable | Metric | Value | Practical interpretation |
|---|---:|---:|---|
| Oxygen transport | Baseline critical radius | 0.799 mm | Center oxygen falls below 0.01 mol m^-3 once organoids exceed roughly 0.8 mm at baseline boundary oxygen. |
| Shear maturation | Best literature-window shear | 0.0452 Pa | Sits within the requested 0.01-0.1 Pa operating window and maximizes utility among sampled conditions. |
| Shear maturation | Composite maturation score at best shear | 0.902 | Indicates strong predicted neural maturation under moderate shear. |
| Shear maturation | Viability at best shear | 0.970 | Suggests minimal modeled damage in the optimal operating region. |

### Additional limitations
- Sensitivity outputs reflect model-form assumptions rather than experimental confidence intervals.
- Pareto optimality depends on the selected weights and damage definition used in the phenomenological model.
