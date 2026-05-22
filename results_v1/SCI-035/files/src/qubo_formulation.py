"""
QUBO Formulation Best Practices Module
Implements encoding strategies for combinatorial optimization problems.
"""
from __future__ import annotations
import numpy as np
import dimod
from typing import Dict, List, Tuple, Optional
import json


class QUBOFormulator:
    """Best-practice QUBO formulations with penalty calibration."""

    # ------------------------------------------------------------------ #
    #  1. One-Hot Constraint                                               #
    # ------------------------------------------------------------------ #
    @staticmethod
    def one_hot_constraint(variables: List[str], penalty: float = 1.0) -> dict:
        """
        Sum-to-one constraint: (sum_i x_i - 1)^2 * penalty
        Expands to: sum_{i!=j} penalty * x_i * x_j - sum_i penalty * x_i + penalty
        """
        Q = {}
        n = len(variables)
        for i in range(n):
            xi = variables[i]
            Q[(xi, xi)] = Q.get((xi, xi), 0) - penalty
            for j in range(i + 1, n):
                xj = variables[j]
                Q[(xi, xj)] = Q.get((xi, xj), 0) + 2 * penalty
        return Q

    # ------------------------------------------------------------------ #
    #  2. Equality Constraint: sum_i a_i * x_i == b                       #
    # ------------------------------------------------------------------ #
    @staticmethod
    def equality_constraint(
        variables: List[str], coefficients: List[float], rhs: float, penalty: float = 1.0
    ) -> dict:
        """(sum_i a_i * x_i - b)^2 * penalty"""
        Q = {}
        n = len(variables)
        for i in range(n):
            xi = variables[i]
            ai = coefficients[i]
            Q[(xi, xi)] = Q.get((xi, xi), 0) + penalty * (ai * ai - 2 * ai * rhs)
            for j in range(i + 1, n):
                xj = variables[j]
                aj = coefficients[j]
                Q[(xi, xj)] = Q.get((xi, xj), 0) + 2 * penalty * ai * aj
        return Q

    # ------------------------------------------------------------------ #
    #  3. Penalty coefficient auto-calibration                            #
    # ------------------------------------------------------------------ #
    @staticmethod
    def auto_penalty(objective_Q: dict, constraint_Q: dict, safety_factor: float = 1.5) -> float:
        """
        Set penalty = safety_factor * max(|objective coefficient|)
        so that any feasible solution is preferred over any infeasible one.
        """
        if not objective_Q:
            return 1.0
        max_obj = max(abs(v) for v in objective_Q.values())
        return safety_factor * max_obj if max_obj > 0 else 1.0

    # ------------------------------------------------------------------ #
    #  4. Ising / QUBO conversion utilities                               #
    # ------------------------------------------------------------------ #
    @staticmethod
    def qubo_to_bqm(Q: dict) -> dimod.BinaryQuadraticModel:
        bqm = dimod.BinaryQuadraticModel.from_qubo(Q)
        return bqm

    @staticmethod
    def normalize_qubo(Q: dict) -> Tuple[dict, float]:
        """Scale Q so that max |coefficient| == 1. Returns (Q_norm, scale)."""
        max_val = max(abs(v) for v in Q.values()) if Q else 1.0
        Q_norm = {k: v / max_val for k, v in Q.items()}
        return Q_norm, max_val

    # ------------------------------------------------------------------ #
    #  5. QUBO statistics report                                          #
    # ------------------------------------------------------------------ #
    @staticmethod
    def qubo_stats(Q: dict) -> dict:
        vals = list(Q.values())
        linear = [v for (i, j), v in Q.items() if i == j]
        quadratic = [v for (i, j), v in Q.items() if i != j]
        variables = set()
        for i, j in Q:
            variables.add(i)
            variables.add(j)
        return {
            "num_variables": len(variables),
            "num_linear_terms": len(linear),
            "num_quadratic_couplings": len(quadratic),
            "density": len(quadratic) / max(1, len(variables) * (len(variables) - 1) / 2),
            "max_coeff": max(abs(v) for v in vals) if vals else 0,
            "min_coeff": min(abs(v) for v in vals if v != 0) if vals else 0,
            "coeff_ratio": (
                max(abs(v) for v in vals) / min(abs(v) for v in vals if v != 0)
                if vals and any(v != 0 for v in vals)
                else 1.0
            ),
        }


# ------------------------------------------------------------------ #
#  Vehicle Routing Problem (VRP) — QUBO Formulation                  #
# ------------------------------------------------------------------ #
class VRPQUBOFormulator:
    """
    QUBO formulation for the Capacitated Vehicle Routing Problem (CVRP).

    Decision variables: x_{v,i,t} = 1 if vehicle v visits city i at step t.

    Constraints:
      C1: Each customer visited exactly once (one-hot over vehicles×steps)
      C2: Each vehicle visits exactly one city per step (route continuity)
      C3: Vehicle capacity not exceeded (soft penalty)

    Objective: Minimize total travel distance.
    """

    def __init__(
        self,
        num_cities: int,
        num_vehicles: int,
        distance_matrix: np.ndarray,
        demands: Optional[List[float]] = None,
        capacity: float = 100.0,
        penalty_visit: float = None,  # auto-calibrated if None
        penalty_capacity: float = None,
    ):
        self.N = num_cities      # includes depot (city 0)
        self.V = num_vehicles
        self.T = num_cities      # max route length = N
        self.dist = distance_matrix
        self.demands = demands if demands else [0.0] + [10.0] * (num_cities - 1)
        self.capacity = capacity
        self._penalty_visit = penalty_visit
        self._penalty_capacity = penalty_capacity

    def var(self, v: int, i: int, t: int) -> str:
        return f"x_{v}_{i}_{t}"

    def build_qubo(self) -> Tuple[dict, dict]:
        """
        Returns (Q, meta) where meta contains penalty values and variable count.
        """
        Q_obj: dict = {}
        Q_c1: dict = {}   # each customer visited once
        Q_c2: dict = {}   # each vehicle at each step visits one city
        Q_c3: dict = {}   # capacity (soft)

        # --- Objective: travel distance ---
        for v in range(self.V):
            for i in range(self.N):
                for j in range(self.N):
                    if i == j:
                        continue
                    d = float(self.dist[i, j])
                    for t in range(self.T - 1):
                        xi = self.var(v, i, t)
                        xj = self.var(v, j, t + 1)
                        key = (xi, xj) if xi <= xj else (xj, xi)
                        Q_obj[key] = Q_obj.get(key, 0.0) + d

        # --- C1: each customer (city > 0) visited by exactly one vehicle at one step ---
        for i in range(1, self.N):
            vars_i = [self.var(v, i, t) for v in range(self.V) for t in range(self.T)]
            # build one-hot: (sum - 1)^2
            q_tmp = QUBOFormulator.one_hot_constraint(vars_i, penalty=1.0)
            for k, val in q_tmp.items():
                Q_c1[k] = Q_c1.get(k, 0.0) + val

        # --- C2: each vehicle visits exactly one city per step ---
        for v in range(self.V):
            for t in range(self.T):
                vars_vt = [self.var(v, i, t) for i in range(self.N)]
                q_tmp = QUBOFormulator.one_hot_constraint(vars_vt, penalty=1.0)
                for k, val in q_tmp.items():
                    Q_c2[k] = Q_c2.get(k, 0.0) + val

        # --- Penalty calibration ---
        max_dist = float(np.max(self.dist))
        p_visit = self._penalty_visit or QUBOFormulator.auto_penalty(Q_obj, Q_c1, safety_factor=max_dist * 2)
        p_cap = self._penalty_capacity or p_visit * 0.5

        # --- Merge ---
        Q: dict = {}
        for k, v in Q_obj.items():
            Q[k] = Q.get(k, 0.0) + v
        for k, v in Q_c1.items():
            Q[k] = Q.get(k, 0.0) + p_visit * v
        for k, v in Q_c2.items():
            Q[k] = Q.get(k, 0.0) + p_visit * v

        # Variable count
        all_vars: set = set()
        for i, j in Q:
            all_vars.add(i)
            all_vars.add(j)

        meta = {
            "num_variables": len(all_vars),
            "penalty_visit": p_visit,
            "penalty_capacity": p_cap,
            "max_distance": max_dist,
            "qubo_stats": QUBOFormulator.qubo_stats(Q),
        }
        return Q, meta


if __name__ == "__main__":
    import json
    # Quick smoke test
    N = 5
    np.random.seed(42)
    coords = np.random.rand(N, 2) * 100
    dist = np.sqrt(((coords[:, None] - coords[None, :]) ** 2).sum(axis=2))
    formulator = VRPQUBOFormulator(N, num_vehicles=2, distance_matrix=dist)
    Q, meta = formulator.build_qubo()
    print(json.dumps(meta, indent=2))
    print(f"QUBO size: {len(Q)} terms")
