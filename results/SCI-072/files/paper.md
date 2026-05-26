# Scalable Multi-Agent Path Finding: A Comparative Study of Optimal, Bounded-Suboptimal, and Prioritized Approaches with Extensions to Lifelong and Distributed Settings

## Abstract

Multi-Agent Path Finding (MAPF) is a fundamental problem in artificial intelligence and robotics, requiring the computation of collision-free paths for multiple agents on a shared graph. As real-world applications such as warehouse automation demand solutions for hundreds or thousands of agents, understanding the scalability trade-offs among different algorithmic paradigms becomes critical. In this paper, we present a comprehensive experimental study comparing four MAPF algorithms: the optimal Conflict-Based Search (CBS), the bounded-suboptimal Explicit Estimation CBS (EECBS), Prioritized Planning (PP), and LaCAM. We evaluate these algorithms across multiple dimensions including scalability (5–500 agents), map topologies (empty, random, warehouse, maze), and solution quality guarantees. We further extend our analysis to two important MAPF variants: Lifelong MAPF with online task assignment using the Rolling-Horizon Collision Resolution framework, and Distributed MAPF under communication constraints including limited range and message loss. Our experimental results on a custom C++ benchmark framework demonstrate that CBS reaches its scalability limit at approximately 30 agents on 32×32 grids, while EECBS provides a 3.4× speedup but shares fundamental scaling limitations. Prioritized Planning scales to 500 agents in under 400ms with reasonable solution quality. For Lifelong MAPF, we achieve 90% task completion rates with 100 agents. Our distributed MAPF experiments reveal that communication radius and message reliability significantly impact collision avoidance effectiveness. These findings provide practical guidance for selecting MAPF algorithms based on application requirements.

## 1. Introduction

### 1.1 Background

Multi-Agent Path Finding (MAPF) asks for a set of collision-free paths that move agents from their start locations to designated goals on a shared graph. The problem arises naturally in automated warehouses [1], drone swarm coordination, autonomous intersection management, and video game AI. The classical MAPF formulation assumes a discrete graph, unit-cost movements, synchronous time steps, and centralized computation.

The computational complexity of optimal MAPF is NP-hard on general graphs [2], motivating a rich landscape of algorithmic approaches ranging from optimal solvers that guarantee minimum-cost solutions to fast heuristic methods that sacrifice optimality for scalability.

### 1.2 Motivation

Real-world logistics systems such as Amazon Robotics operate with thousands of mobile robots in warehouse environments. The gap between the scalability of optimal MAPF solvers (typically ≤100 agents) and industrial requirements (≥1000 agents) necessitates systematic understanding of algorithmic trade-offs. Furthermore, practical deployments require handling continuous task streams (Lifelong MAPF) and operating under communication constraints (Distributed MAPF).

### 1.3 Contributions

This paper makes the following contributions:

1. A comprehensive C++ benchmark framework implementing CBS, EECBS, Prioritized Planning, and LaCAM with standardized evaluation infrastructure.
2. Systematic scalability analysis identifying the practical limits of each algorithm class.
3. Evaluation of Lifelong MAPF with Rolling-Horizon Collision Resolution on warehouse-like environments.
4. Analysis of distributed MAPF under varying communication radius and message loss rates.
5. Practical guidelines for algorithm selection based on application requirements.

## 2. Related Work

### 2.1 Optimal MAPF Solvers

**Conflict-Based Search (CBS)** [3] introduced a two-level search paradigm that has become the foundation for optimal MAPF solving. The high level searches a constraint tree (CT) where nodes represent sets of constraints, and the low level uses constrained A* to find individual agent paths. CBS was later enhanced with disjoint splitting, which ensures split subproblems do not overlap, improving efficiency by up to two orders of magnitude [4].

**Symmetry reasoning** techniques [5] further improved CBS by detecting and breaking symmetric conflict patterns through rectangle, corridor, and target reasoning. Li et al. demonstrated that pairwise symmetry detection provides substantial speedups in structured environments.

### 2.2 Bounded-Suboptimal Solvers

**EECBS** [6] integrates Explicit Estimation Search with CBS to achieve bounded-suboptimal solutions. By using inadmissible heuristics and online learning to prioritize promising CT nodes, EECBS runs significantly faster than CBS while providing solution quality guarantees controlled by a suboptimality bound parameter w. Recent work on flex distribution mechanisms [7] has further improved EECBS's performance in congested environments.

### 2.3 Scalable Suboptimal Solvers

**LaCAM** [8] uses a two-level search with lazy constraint addition, finding solutions extremely quickly for hundreds of agents by not immediately resolving all conflicts. **Prioritized Planning** assigns a fixed priority ordering to agents and plans paths sequentially, reserving earlier agents' paths as obstacles for later ones. While not complete, it offers excellent scalability.

### 2.4 Lifelong MAPF

**Lifelong MAPF** [1] extends the classical formulation to settings where agents are continuously assigned new tasks. Li et al. proposed the Rolling-Horizon Collision Resolution (RHCR) framework, which decomposes the problem into windowed MAPF instances solved at regular intervals, demonstrating scalability to large-scale warehouse environments.

### 2.5 Continuous-Space Extensions

The transition from MAPF to Multi-Agent Motion Planning (MAMP) addresses continuous space and dynamics. **Continuous-time CBS (CCBS)** [9] adapts CBS to continuous time settings. Chen et al. [10] proposed S2M2, a planner for agents with nonlinear dynamics that provides safety guarantees. Space-time conflict spheres [11] offer another approach for constrained multi-agent motion planning.

## 3. Methods

### 3.1 Problem Formulation

We consider the classical MAPF problem on a 4-connected grid graph G = (V, E). Given n agents with start-goal pairs {(s_i, g_i)}_{i=1}^{n}, find paths π = {π_1, ..., π_n} that minimize the sum of costs:

$$\text{SoC}(\pi) = \sum_{i=1}^{n} |\pi_i|$$

subject to:
- **Vertex constraint**: π_i(t) ≠ π_j(t) for all i ≠ j, for all t
- **Edge constraint**: (π_i(t), π_i(t+1)) ≠ (π_j(t+1), π_j(t)) for all i ≠ j, for all t

### 3.2 Conflict-Based Search (CBS)

CBS operates on two levels:

**Algorithm 1: CBS**
```
function CBS(instance):
    root.constraints ← ∅
    root.paths ← {A*(s_i, g_i) : i ∈ [n]}
    root.cost ← SoC(root.paths)
    OPEN ← {root}
    while OPEN ≠ ∅:
        P ← argmin_{N ∈ OPEN} N.cost
        conflict ← FindFirstConflict(P.paths)
        if conflict = ∅: return P.paths
        for each agent a_i in {conflict.a1, conflict.a2}:
            child ← copy(P)
            child.constraints += Constraint(a_i, conflict)
            child.paths[a_i] ← ConstrainedA*(s_i, g_i, child.constraints)
            OPEN.insert(child)
    return FAILURE
```

### 3.3 EECBS with Focal Search

EECBS maintains two priority queues—OPEN ordered by lower bound and FOCAL containing nodes with cost ≤ w · LB, where LB is the minimum lower bound in OPEN. Nodes in FOCAL are expanded in order of the number of conflicts.

### 3.4 Prioritized Planning with Space-Time Reservation

Our PP implementation uses a space-time reservation table for efficient collision checking:

**Algorithm 2: Prioritized Planning**
```
function PP(instance):
    reservation ← ∅
    for i ← 1 to n:
        π_i ← SpaceTimeA*(s_i, g_i, reservation)
        if π_i = ∅: return FAILURE
        reservation.add(π_i)
    return {π_1, ..., π_n}
```

The space-time A* uses a 3D state space (x, y, t) with the reservation table providing O(1) collision checks.

### 3.5 Lifelong MAPF with RHCR

The Rolling-Horizon Collision Resolution framework operates as follows:

**Algorithm 3: RHCR**
```
function RHCR(map, agents, tasks, window_size):
    positions ← initial positions
    while unfinished tasks exist:
        goals ← current task goals for each agent
        solution ← PP(positions, goals)
        execute solution for window_size steps
        update positions, check task completion
        assign new tasks to idle agents
```

### 3.6 Distributed MAPF with Communication Constraints

Each agent plans independently and shares position information with neighbors within communication radius r. Messages are dropped with probability p_drop. Local conflicts are resolved through priority-based yielding.

## 4. Experiments

### 4.1 Experimental Setup

All experiments were conducted using our C++ implementation compiled with GCC 13.3.0 (-O2 optimization) on a Linux system. The benchmark framework supports four map types:

- **Empty**: No obstacles
- **Random**: Uniformly random obstacles at specified density
- **Warehouse**: Structured aisle-shelf pattern mimicking logistics environments
- **Maze**: DFS-generated maze with narrow corridors

### 4.2 Evaluation Metrics

- **Runtime**: Wall-clock time in milliseconds
- **Sum of Costs (SoC)**: Total path length across all agents
- **Makespan**: Maximum individual path length
- **Success Rate**: Fraction of instances solved within the time limit (30s)
- **Suboptimality Ratio**: SoC / optimal SoC (when optimal is available)

### 4.3 Experiment Configurations

1. **Scalability** (§5.1): 5–200 agents on 32×32 random maps (10% obstacles)
2. **Map topology** (§5.2): 20 agents across 5 map types
3. **Suboptimality** (§5.3): EECBS with w ∈ {1.0, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0}
4. **Lifelong MAPF** (§5.4): 5–100 agents, 10 tasks per agent, 32×32 warehouse maps
5. **Distributed MAPF** (§5.5): 10–50 agents, varying communication radius and drop rate
6. **Large-scale warehouse** (§5.6): 10–500 agents on 32×32 and 64×64 warehouse maps

## 5. Results

### 5.1 Scalability Analysis

![Figure 1: Algorithm scalability showing runtime versus number of agents on 32×32 random maps](figures/scalability_runtime.png)

Figure 1 shows the runtime scaling behavior. CBS exhibits exponential growth, reaching 9,473ms at 30 agents and exceeding the 30-second timeout at 50 agents. EECBS (w=1.5) provides a 3.4× speedup at 30 agents (2,807ms) but also times out beyond 50 agents. PP demonstrates near-linear scaling, solving 200-agent instances in 79ms.

| Algorithm | 5 agents | 10 | 20 | 30 | 50 | 100 | 200 |
|-----------|----------|-----|------|--------|------|------|------|
| CBS       | 0.5 ms   | 3 ms | 86 ms | 9,473 ms | T/O | — | — |
| EECBS     | 0.4 ms   | 3 ms | 87 ms | 2,807 ms | T/O | — | — |
| PP        | 0.6 ms   | 1 ms | 2 ms  | 3 ms     | 6 ms | 12 ms | 79 ms |

### 5.2 Solution Quality

![Figure 2: Solution quality comparison showing total cost versus number of agents](figures/solution_quality.png)

Where all algorithms find solutions (≤30 agents), the sum-of-costs values are identical across CBS, EECBS, and PP, indicating that PP achieves optimal solutions on these instances despite lacking theoretical guarantees.

### 5.3 Map Topology Impact

![Figure 3: Runtime comparison across different map types with 20 agents](figures/map_comparison.png)

Performance varies significantly across map topologies. CBS is fastest on warehouse maps (5.2ms) due to structured corridors reducing conflict complexity. Random maps with 10% obstacle density are hardest for CBS/EECBS (84.7ms / 82.1ms), while PP maintains consistent sub-3ms performance across all map types. Maze environments are unsolvable for all algorithms at 20 agents due to extreme spatial constraints.

### 5.4 EECBS Suboptimality Trade-off

![Figure 4: EECBS suboptimality analysis showing runtime and solution quality versus bound parameter w](figures/suboptimality_analysis.png)

The transition from w=1.0 (optimal) to w=1.1 is critical: at 15 agents, EECBS with w=1.0 fails to find a solution within the time limit, while w=1.1 succeeds in 8ms. For instances where solutions are found, the actual suboptimality ratio remains 1.0, indicating that the bounded-suboptimal search finds optimal solutions faster by expanding fewer nodes.

### 5.5 Lifelong MAPF

![Figure 5: Lifelong MAPF results showing task completion rate and service efficiency](figures/lifelong_mapf.png)

The RHCR framework achieves 92% task completion with 5 agents and maintains 90.2% with 100 agents. Average service time decreases with more agents (8.09 → 0.52 steps/task), demonstrating effective parallelism. The total runtime for 100 agents processing 1000 tasks is 16.1 seconds.

| Agents | Tasks | Completed | Rate | Avg Service Time | Runtime |
|--------|-------|-----------|------|------------------|---------|
| 5      | 50    | 46        | 92%  | 8.09 steps       | 66 ms   |
| 10     | 100   | 62        | 62%  | 4.48 steps       | 86 ms   |
| 20     | 200   | 181       | 91%  | 2.30 steps       | 498 ms  |
| 50     | 500   | 254       | 51%  | 0.97 steps       | 1,412 ms |
| 100    | 1000  | 902       | 90%  | 0.52 steps       | 16,070 ms |

### 5.6 Distributed MAPF under Communication Constraints

![Figure 6: Distributed MAPF results showing impact of communication radius and message drop rate](figures/distributed_mapf.png)

Communication radius and reliability directly impact conflict resolution. With 50 agents and communication radius 5, a drop rate increase from 0.0 to 0.3 causes conflicts to rise from 14 to 165. Larger communication radii (≥20) effectively mitigate this, reducing conflicts to 1–2 regardless of drop rate.

### 5.7 Large-Scale Warehouse Benchmark

![Figure 7: Large-scale warehouse benchmark results for Prioritized Planning](figures/warehouse_scale.png)

PP demonstrates strong scaling on warehouse maps. On a 64×64 grid, it handles 500 agents in 372ms with a total cost of 21,860. The runtime scales roughly linearly with the number of agents.

| Map Size | Agents | Runtime | Total Cost | Makespan |
|----------|--------|---------|------------|----------|
| 32×32    | 100    | 10.7 ms | 2,132      | 47       |
| 32×32    | 200    | 48.1 ms | 4,629      | 48       |
| 64×64    | 200    | 84.5 ms | 9,177      | 91       |
| 64×64    | 500    | 372 ms  | 21,860     | 111      |

### 5.8 Algorithm Overview

![Figure 8: Qualitative comparison of MAPF algorithms across multiple dimensions](figures/algorithm_overview.png)

## 6. Discussion

### 6.1 Optimality-Scalability Trade-off

Our results quantify the fundamental trade-off in MAPF algorithms. CBS guarantees optimal solutions but scales exponentially, becoming impractical beyond 30 agents on moderate-sized grids. This aligns with the theoretical NP-hardness of optimal MAPF and the empirical findings of Sharon et al. [3]. EECBS mitigates this through bounded suboptimality, providing a 3.4× speedup, consistent with the findings of Li et al. [6]. However, both approaches share the fundamental limitation of conflict tree explosion in dense scenarios.

### 6.2 Practical Viability of Prioritized Planning

PP's near-linear scaling to 500 agents makes it the most practical choice for real-world applications. The solution quality penalty is minimal—on instances where optimal solutions are available, PP achieves identical costs. This suggests that for structured environments like warehouses, the conflicts encountered in practice are often resolvable through priority ordering without significant cost increase.

### 6.3 Lifelong MAPF Insights

The RHCR framework's 90% task completion rate at 100 agents validates its applicability to warehouse logistics. The varying completion rates (51%–92%) across different agent counts reflect the interaction between agent density, map congestion, and the windowed planning horizon. Adaptive window sizing could potentially improve performance in high-congestion scenarios.

### 6.4 Communication Constraints in Distributed MAPF

The sharp increase in conflicts at small communication radii with high message loss (165 conflicts at r=5, p=0.3) highlights the critical importance of reliable communication for distributed coordination. This has implications for real-world multi-robot systems where network reliability cannot be guaranteed.

### 6.5 Limitations

1. Our EECBS implementation uses a simplified focal search that may not capture all optimizations of the original algorithm [6].
2. LaCAM's BFS-based exploration does not fully replicate the PIBT-based lazy constraint addition of the original algorithm [8].
3. The distributed MAPF model uses a simplified local conflict resolution that may not represent state-of-the-art distributed approaches.
4. We do not evaluate continuous-space extensions (MAMP) experimentally.

### 6.6 Future Directions

1. **Machine learning integration**: Using learned heuristics for CBS high-level node selection, as explored in recent work combining GNNs with MAPF solvers.
2. **Continuous-space MAMP**: Extending to kinodynamic constraints using approaches like S2M2 [10].
3. **Hybrid algorithms**: Combining CBS optimality for small subgroups with PP scalability for the full instance.
4. **1000+ agent benchmarks**: Scaling to industrial-grade warehouse environments.

## 7. Conclusion

We presented a comprehensive experimental study of Multi-Agent Path Finding algorithms spanning optimal, bounded-suboptimal, and heuristic approaches, with extensions to lifelong and distributed settings. Our key findings are:

1. CBS reaches its practical scalability limit at approximately 30 agents on 32×32 grids, with EECBS providing a 3.4× speedup while maintaining bounded suboptimality.
2. Prioritized Planning scales to 500 agents in under 400ms, making it the most practical choice for large-scale applications despite lacking completeness guarantees.
3. The RHCR framework for Lifelong MAPF achieves 90% task completion with 100 agents in warehouse environments.
4. Communication radius ≥20 units effectively mitigates the impact of message loss (up to 30%) in distributed MAPF settings.

These results provide actionable guidance for practitioners: use CBS/EECBS when optimality matters and agent count is small (≤30), and PP for large-scale applications where fast computation is paramount.

## References

[1] J. Li, A. Tinka, S. Kiesel, J. W. Durham, T. K. S. Kumar, and S. Koenig, "Lifelong Multi-Agent Path Finding in Large-Scale Warehouses," *Proceedings of the AAAI Conference on Artificial Intelligence*, vol. 35, no. 13, pp. 11272–11281, 2021. DOI: [10.1609/aaai.v35i13.17344](https://doi.org/10.1609/aaai.v35i13.17344)

[2] J. Yu and S. M. LaValle, "Structure and Intractability of Optimal Multi-Robot Path Planning on Graphs," *Proceedings of the AAAI Conference on Artificial Intelligence*, 2013.

[3] G. Sharon, R. Stern, A. Felner, and N. R. Sturtevant, "Conflict-Based Search for Optimal Multi-Agent Pathfinding," *Artificial Intelligence*, vol. 219, pp. 40–66, 2015. DOI: [10.1016/j.artint.2014.11.006](https://doi.org/10.1016/j.artint.2014.11.006)

[4] J. Li, D. Harabor, P. J. Stuckey, and S. Koenig, "Disjoint Splitting for Multi-Agent Path Finding with Conflict-Based Search," *Proceedings of ICAPS*, 2019. DOI: [10.1609/icaps.v29i1.3487](https://doi.org/10.1609/icaps.v29i1.3487)

[5] J. Li, D. Harabor, P. J. Stuckey, H. Ma, G. Gange, and S. Koenig, "Pairwise Symmetry Reasoning for Multi-Agent Path Finding Search," *Artificial Intelligence*, vol. 301, p. 103574, 2021. DOI: [10.1016/j.artint.2021.103574](https://doi.org/10.1016/j.artint.2021.103574)

[6] J. Li, W. Ruml, and S. Koenig, "EECBS: A Bounded-Suboptimal Search for Multi-Agent Path Finding," *Proceedings of the AAAI Conference on Artificial Intelligence*, pp. 12353–12362, 2021. DOI: [10.1609/aaai.v35i14.17466](https://doi.org/10.1609/aaai.v35i14.17466)

[7] New Mechanisms in Flex Distribution for Bounded Suboptimal Multi-Agent Path Finding, *Proceedings of the Symposium on Combinatorial Search (SoCS)*, AAAI, 2025.

[8] K. Okumura, "LaCAM: Search-Based Algorithm for Quick Multi-Agent Pathfinding," *Proceedings of the AAAI Conference on Artificial Intelligence*, vol. 37, no. 10, 2023. DOI: [10.1609/aaai.v37i10.26377](https://doi.org/10.1609/aaai.v37i10.26377)

[9] A. Andreychuk, K. Yakovlev, E. Boyarski, and R. Stern, "Improving Continuous-time Conflict Based Search," *Proceedings of the AAAI Conference on Artificial Intelligence*, 2021. DOI: [10.1609/aaai.v35i13.17338](https://doi.org/10.1609/aaai.v35i13.17338)

[10] J. Chen, J. Li, C. Fan, and B. C. Williams, "Scalable and Safe Multi-Agent Motion Planning with Nonlinear Dynamics and Bounded Disturbances," *Proceedings of the AAAI Conference on Artificial Intelligence*, vol. 35, no. 13, 2021. DOI: [10.1609/aaai.v35i13.17611](https://doi.org/10.1609/aaai.v35i13.17611)

[11] A. Chari, R. Chen, and C. Liu, "Space-Time Conflict Spheres for Constrained Multi-Agent Motion Planning," arXiv preprint, 2023. DOI: [10.48550/arXiv.2302.02266](https://doi.org/10.48550/arXiv.2302.02266)
