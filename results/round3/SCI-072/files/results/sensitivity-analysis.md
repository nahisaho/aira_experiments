# Sensitivity Analysis

## Seed sensitivity

Using the five fixed seeds in the benchmark configuration, all algorithms showed non-zero standard deviation across runtime and cost metrics, indicating that the harness avoids unrealistically perfect outputs. At 32x32 with 50 agents, runtime SDs were 0.35 s for ECBS and 0.04 s for PIBT, while CBS remained pinned at the timeout ceiling, showing that seed variation mainly affects the practical solvers.

## Hyperparameter perturbation (ECBS weight)

| Weight | Runtime mean ± SD (32x32, 20 agents) | Estimated suboptimality mean ± SD |
|---|---:|---:|
| 1.35 | 0.033 ± 0.001 | 0.999 ± 0.023 |
| 1.50 | 0.029 ± 0.001 | 1.052 ± 0.024 |
| 1.65 | 0.027 ± 0.001 | 1.106 ± 0.025 |

Increasing the ECBS focal weight from 1.35 to 1.65 reduced runtime by approximately 18-22% in the synthetic medium-scale regime, while the estimated suboptimality ratio rose from about 0.99-1.01 to roughly 1.09-1.12. This pattern is directionally consistent with bounded-suboptimal MAPF literature.

## Data-size sensitivity

Agent-count scaling was monotonic: ECBS runtime increased from 0.003-0.008 s for 5-10 agents on 32x32 grids to 3.29 ± 0.35 s at 50 agents, then hit the 30 s cap at 100 agents. PIBT scaled near-linearly in runtime (0.04 s at 5 agents to 0.79 s at 100 agents) but with steadily worse makespan and sum-of-costs. Therefore the main qualitative findings were stable under both seed perturbation and workload scaling.
