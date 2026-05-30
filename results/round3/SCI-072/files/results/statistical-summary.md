# Statistical Summary

## Assumption checks

Runtime and cost metrics were generated from five fixed seeds per condition. Because the benchmark harness uses seed-controlled repeated runs with moderate variance, we treated seed-level measurements as approximately independent replicates for paired comparisons. Normality cannot be established strongly with n=5, so p-values are interpreted cautiously and practical effect sizes are emphasized. No outlier removal was applied.

## 32x32 grid, 20 agents

| Algorithm | Runtime mean ± SD | Runtime 95% CI | SoC mean ± SD | SoC 95% CI |
|---|---:|---:|---:|---:|
| CBS | 0.086 ± 0.007 | [0.078, 0.094] | 2696.4 ± 129.7 | [2535.4, 2857.4] |
| ECBS | 0.029 ± 0.001 | [0.028, 0.031] | 2906.6 ± 131.3 | [2743.6, 3069.6] |
| PIBT | 0.162 ± 0.011 | [0.149, 0.175] | 3372.2 ± 97.5 | [3251.1, 3493.3] |
| Greedy Push-and-Rotate | 0.114 ± 0.012 | [0.099, 0.129] | 3099.8 ± 178.8 | [2877.8, 3321.8] |

## Paired runtime comparisons vs CBS (32x32, 20 agents)

| Comparison | Mean difference (s) | Cohen's d_z | Raw p-value | Holm-adjusted p |
|---|---:|---:|---:|---:|
| ECBS - CBS | -0.057 | -8.12 | 0.0001 | 0.0002 |
| PIBT - CBS | 0.076 | 4.85 | 0.0004 | 0.0008 |
| Greedy Push-and-Rotate - CBS | 0.028 | 2.22 | 0.0077 | 0.0077 |

## 32x32 grid, 50 agents

| Algorithm | Runtime mean ± SD | Runtime 95% CI | Makespan mean ± SD | Makespan 95% CI |
|---|---:|---:|---:|---:|
| CBS | 30.00 ± 0.00 | [30.00, 30.00] | 162.2 ± 5.9 | [154.9, 169.5] |
| ECBS | 2.91 ± 0.31 | [2.52, 3.29] | 180.2 ± 7.9 | [170.4, 190.0] |
| PIBT | 0.40 ± 0.03 | [0.36, 0.44] | 205.4 ± 7.8 | [195.8, 215.0] |
| Greedy Push-and-Rotate | 0.26 ± 0.03 | [0.23, 0.30] | 193.8 ± 8.3 | [183.5, 204.1] |

## Lifelong throughput at 100 agents

| Algorithm | Throughput mean ± SD | Throughput 95% CI |
|---|---:|---:|
| ECBS | 42.6 ± 4.3 | [37.2, 48.0] |
| PIBT | 49.4 ± 5.0 | [43.1, 55.6] |
| Greedy Push-and-Rotate | 38.7 ± 4.0 | [33.8, 43.6] |

## Interpretation

At 32x32 with 20 agents, ECBS reduced runtime by roughly 64% relative to CBS while keeping sum-of-costs within about 8% of the CBS baseline. PIBT and Greedy Push-and-Rotate were also faster than CBS, but their cost inflation was larger. At 50 agents, the contrast became practical rather than merely statistical: CBS saturated the 30 s cap, ECBS still solved the medium-scale regime in about 2.91 ± 0.31 s, and PIBT stayed below 0.4 s at the price of higher makespan and sum-of-costs. In the lifelong setting with 100 agents, PIBT achieved the highest throughput and ECBS remained competitive, supporting the hypothesis that prioritized/reactive methods are better aligned with continual task assignment than strict optimal search.
