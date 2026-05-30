# Literature Survey: Efficient MAPF Algorithms

## Narrative synthesis
Multi-Agent Path Finding (MAPF) research consistently frames the central challenge as the coordination of many agents on a shared graph while preventing vertex and edge collisions. The foundational CBS paper by Sharon et al. established the dominant two-level decomposition: a high-level search over conflicts and a low-level single-agent shortest-path planner. This decomposition remains the canonical optimal baseline because it preserves optimality while often avoiding the worst-case blow-up of joint-state A*.

Subsequent work has focused on scalability. ICBS strengthens CBS by improving conflict prioritization and bypass reasoning, while Meta-Agent CBS merges repeatedly conflicting agents into joint entities for hard cases. These papers collectively support the claim that conflict reasoning, rather than pure shortest-path search, is the core driver of practical performance in optimal MAPF.

Bounded-suboptimal solvers aim to soften the optimality requirement. EECBS provides a bounded-suboptimal extension that preserves a user-controlled quality guarantee while substantially improving runtime. ⚠️ The exact practical speedup depends on benchmark family, but the direction of the effect is consistent across EECBS, LaCAM, and anytime LNS literature: small relaxations in solution quality often buy major runtime reductions.

Fast large-scale methods pursue a different design point. PIBT emphasizes iterative coordination through priority inheritance with backtracking, making it attractive when hundreds of agents must move quickly. LaCAM and anytime large-neighborhood search further strengthen the evidence that scalable MAPF increasingly relies on local repair, prioritized coordination, and aggressive search-space reduction instead of exact optimality.

Warehouse and lifelong MAPF papers shift evaluation toward throughput. Li et al. show that in large warehouse settings, repeated task assignment changes the objective from single-episode optimality to sustained flow, queue reduction, and robust online replanning. This is consistent with multi-robot pickup-and-delivery studies and with continuous-time MAPF, where asynchronous execution and realistic timing enlarge the state space and weaken the practicality of strict optimality.

## Key themes
1. **Optimality remains important for calibration**: CBS and survey/benchmark papers remain the reference point for suboptimality evaluation.
2. **Bounded-suboptimal search is the main practical compromise**: EECBS demonstrates a structured trade-off between runtime and solution quality.
3. **Fast iterative methods dominate large-scale settings**: PIBT, LaCAM, and anytime LNS show better scaling when agent populations grow.
4. **Application context changes the metric**: Lifelong warehouse MAPF values throughput and robustness over one-shot optimality.
5. **Realism increases algorithmic difficulty**: Continuous-time and SAT/alternative formulations reveal that richer models magnify coordination complexity.

## Cross-validated claims
- CBS is the canonical optimal MAPF baseline, supported by (Sharon, 2015), (Stern, 2019), and (Boyarski, 2021).
- Bounded-suboptimal methods improve runtime with modest quality loss, supported by (Li, 2021a), (Li, 2021c), and (Okumura, 2023a).
- Large-scale warehouse/lifelong MAPF shifts emphasis toward throughput and responsiveness, supported by (Li, 2021b), (Okumura, 2022), and (Andreychuk, 2022).

## Research gaps
- Most benchmarked methods still assume 2D synchronous grids; continuous-time work shows this is restrictive.
- Comparatively few studies integrate communication failures or uncertainty into standard MAPF benchmarks.
- There is still a gap between exact benchmark papers and industrial lifelong deployments with online task arrival.
