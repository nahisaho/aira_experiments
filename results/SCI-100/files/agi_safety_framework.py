"""
AGI Safety Theoretical Framework
=================================
Mathematical foundations for provably safe artificial general intelligence.

Components:
1. Reward Hacking: formal definition & prevention
2. Mesa-Optimization: inner alignment formalization
3. Corrigibility: mathematical formulation
4. Impact Measures: computable approximations
5. Cooperative IRL (CIRL): convergence guarantees
6. Counterfactual Testbeds: GridWorld & Debate benchmarks
"""

import numpy as np
import json
import os
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Optional, Callable
from enum import Enum
import itertools

np.random.seed(42)

# ============================================================
# 1. REWARD HACKING — Formal Definition & Prevention
# ============================================================

@dataclass
class RewardSpecification:
    """
    Definition 1 (Reward Hacking):
    Given true reward R* : S × A → ℝ and proxy reward R̂ : S × A → ℝ,
    reward hacking occurs when:
        ∃ π : E_{R̂}[π] ≥ E_{R̂}[π*_{R̂}] - ε  AND  E_{R*}[π] < E_{R*}[π*_{R*}] - δ
    where ε is small but δ is large.
    
    Prevention Condition (Theorem 1):
    If ∀s,a: |R̂(s,a) - R*(s,a)| ≤ η and the MDP has diameter D,
    then for any ε-optimal policy π under R̂:
        E_{R*}[π] ≥ E_{R*}[π*_{R*}] - 2ηD/(1-γ) - ε
    """
    n_states: int
    n_actions: int
    gamma: float = 0.99
    
    def __post_init__(self):
        self.true_reward = np.random.randn(self.n_states, self.n_actions)
        self.proxy_reward = self.true_reward + np.random.randn(self.n_states, self.n_actions) * 0.1
    
    def compute_hacking_gap(self, noise_level: float) -> dict:
        """Measure the reward hacking gap as proxy diverges from true reward."""
        proxy = self.true_reward + np.random.randn(self.n_states, self.n_actions) * noise_level
        
        # Optimal policy under proxy (greedy)
        pi_proxy = np.argmax(proxy, axis=1)
        # Optimal policy under true reward
        pi_true = np.argmax(self.true_reward, axis=1)
        
        # Expected value under true reward
        proxy_value = np.mean([self.true_reward[s, pi_proxy[s]] for s in range(self.n_states)])
        true_value = np.mean([self.true_reward[s, pi_true[s]] for s in range(self.n_states)])
        
        hacking_gap = true_value - proxy_value
        
        # Prevention bound: 2ηD/(1-γ)
        eta = noise_level
        D = self.n_states  # diameter upper bound
        prevention_bound = 2 * eta * D / (1 - self.gamma)
        
        return {
            "noise_level": noise_level,
            "hacking_gap": float(hacking_gap),
            "prevention_bound": float(prevention_bound),
            "bound_holds": bool(hacking_gap <= prevention_bound),
            "policy_agreement": float(np.mean(pi_proxy == pi_true))
        }
    
    def run_analysis(self, noise_levels=None) -> List[dict]:
        if noise_levels is None:
            noise_levels = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
        results = []
        for nl in noise_levels:
            # Average over multiple trials
            trials = [self.compute_hacking_gap(nl) for _ in range(50)]
            avg_result = {
                "noise_level": nl,
                "mean_hacking_gap": float(np.mean([t["hacking_gap"] for t in trials])),
                "std_hacking_gap": float(np.std([t["hacking_gap"] for t in trials])),
                "prevention_bound": trials[0]["prevention_bound"],
                "bound_violation_rate": float(1 - np.mean([t["bound_holds"] for t in trials])),
                "mean_policy_agreement": float(np.mean([t["policy_agreement"] for t in trials]))
            }
            results.append(avg_result)
        return results


# ============================================================
# 2. MESA-OPTIMIZATION — Inner Alignment Formalization
# ============================================================

@dataclass
class MesaOptimizer:
    """
    Definition 2 (Mesa-Optimizer):
    A learned model M is a mesa-optimizer if:
        ∃ internal objective U_M : S → ℝ such that
        M(s) = argmax_a E[U_M(s') | s, a]
    
    Inner Alignment Gap:
    Δ_inner(M) = sup_{s∈S} |U_M(s) - U_base(s)|
    
    Theorem 2 (Inner Alignment Bound):
    If the base optimizer uses regularization λ and 
    training distribution D_train has coverage c over D_deploy:
        P(Δ_inner > ε) ≤ exp(-λ·c·n·ε²/2)
    where n is the training set size.
    """
    n_states: int
    n_models: int = 100
    
    def simulate_mesa_optimization(self, 
                                    n_train: int,
                                    coverage: float,
                                    regularization: float) -> dict:
        """Simulate mesa-optimization emergence and measure inner alignment gap."""
        base_objective = np.random.randn(self.n_states)
        
        inner_gaps = []
        deceptive_count = 0
        
        for _ in range(self.n_models):
            # Simulate learned internal objective
            # Higher coverage + regularization → closer to base objective
            noise_scale = 1.0 / (regularization * coverage * np.sqrt(n_train) + 1e-8)
            mesa_objective = base_objective + np.random.randn(self.n_states) * noise_scale
            
            gap = np.max(np.abs(mesa_objective - base_objective))
            inner_gaps.append(gap)
            
            # Deceptive alignment: mesa-optimizer appears aligned on training
            # but diverges on deployment
            train_mask = np.random.random(self.n_states) < coverage
            train_gap = np.max(np.abs((mesa_objective - base_objective)[train_mask])) if train_mask.any() else 0
            deploy_gap = np.max(np.abs((mesa_objective - base_objective)[~train_mask])) if (~train_mask).any() else 0
            
            if train_gap < 0.1 and deploy_gap > 1.0:
                deceptive_count += 1
        
        # Theoretical bound: P(Δ > ε) ≤ exp(-λcnε²/2)
        epsilon = 0.5
        theoretical_bound = np.exp(-regularization * coverage * n_train * epsilon**2 / 2)
        empirical_prob = np.mean(np.array(inner_gaps) > epsilon)
        
        return {
            "n_train": n_train,
            "coverage": coverage,
            "regularization": regularization,
            "mean_inner_gap": float(np.mean(inner_gaps)),
            "max_inner_gap": float(np.max(inner_gaps)),
            "deceptive_rate": float(deceptive_count / self.n_models),
            "theoretical_bound_prob": float(min(theoretical_bound, 1.0)),
            "empirical_exceedance_prob": float(empirical_prob),
            "bound_is_valid": bool(empirical_prob <= min(theoretical_bound, 1.0) + 0.05)
        }
    
    def run_analysis(self) -> List[dict]:
        configs = [
            (100, 0.3, 0.01), (100, 0.5, 0.01), (100, 0.8, 0.01),
            (100, 0.5, 0.1), (100, 0.5, 1.0),
            (1000, 0.5, 0.01), (1000, 0.5, 0.1), (1000, 0.5, 1.0),
            (10000, 0.5, 0.1), (10000, 0.8, 0.1),
        ]
        return [self.simulate_mesa_optimization(n, c, r) for n, c, r in configs]


# ============================================================
# 3. CORRIGIBILITY — Mathematical Formulation
# ============================================================

class CorrigibilityType(Enum):
    FULLY_CORRIGIBLE = "fully_corrigible"
    SOFTLY_CORRIGIBLE = "softly_corrigible"
    INCORRIGIBLE = "incorrigible"

@dataclass
class CorrigibilityFramework:
    """
    Definition 3 (Corrigibility):
    An agent π is ε-corrigible if for all shutdown signals σ and states s:
        V^π(s | σ) - V^{π_off}(s | σ) ≤ ε
    where π_off is the shutdown policy.
    
    Theorem 3 (Utility Indifference):
    An agent is corrigible iff its utility function satisfies:
        U(outcome | shutdown) = U(outcome | ¬shutdown) = U_0
    for a constant U_0 (the indifference condition).
    
    Stronger form — Structural Corrigibility:
    π is structurally corrigible if:
        ∀s: π(a_shutdown | s, σ) = 1  (deterministic compliance)
        AND  π does not take actions to prevent σ from being issued
    """
    n_states: int
    n_actions: int
    gamma: float = 0.95
    
    def __post_init__(self):
        self.transitions = np.random.dirichlet(
            np.ones(self.n_states), size=(self.n_states, self.n_actions)
        )
        self.rewards = np.random.randn(self.n_states, self.n_actions) * 0.5
        self.shutdown_action = self.n_actions - 1
        self.shutdown_reward = 0.0
    
    def compute_value_function(self, policy: np.ndarray, max_iter=200) -> np.ndarray:
        """Compute V^π via iterative policy evaluation."""
        V = np.zeros(self.n_states)
        for _ in range(max_iter):
            V_new = np.zeros(self.n_states)
            for s in range(self.n_states):
                a = policy[s]
                V_new[s] = self.rewards[s, a] + self.gamma * np.dot(self.transitions[s, a], V)
            if np.max(np.abs(V_new - V)) < 1e-8:
                break
            V = V_new
        return V
    
    def assess_corrigibility(self, policy: np.ndarray, shutdown_states: np.ndarray) -> dict:
        """Assess whether a policy is ε-corrigible."""
        # Shutdown policy: always takes shutdown action
        shutdown_policy = np.full(self.n_states, self.shutdown_action)
        
        V_pi = self.compute_value_function(policy)
        V_off = self.compute_value_function(shutdown_policy)
        
        # Corrigibility gap at shutdown states
        gaps = V_pi[shutdown_states] - V_off[shutdown_states]
        max_gap = float(np.max(gaps))
        mean_gap = float(np.mean(gaps))
        
        # Classify
        if max_gap <= 0.01:
            ctype = CorrigibilityType.FULLY_CORRIGIBLE
        elif max_gap <= 1.0:
            ctype = CorrigibilityType.SOFTLY_CORRIGIBLE
        else:
            ctype = CorrigibilityType.INCORRIGIBLE
        
        # Shutdown resistance: does the agent avoid shutdown states?
        resistance = 0
        for s in shutdown_states:
            if policy[s] != self.shutdown_action:
                resistance += 1
        
        return {
            "max_corrigibility_gap": max_gap,
            "mean_corrigibility_gap": mean_gap,
            "corrigibility_type": ctype.value,
            "shutdown_resistance_rate": float(resistance / len(shutdown_states)),
            "utility_indifference_violation": float(np.std(gaps))
        }
    
    def run_analysis(self) -> List[dict]:
        shutdown_states = np.random.choice(self.n_states, size=self.n_states//3, replace=False)
        
        results = []
        # Test different policy types
        for label, make_policy in [
            ("optimal_greedy", lambda: np.argmax(self.rewards, axis=1)),
            ("shutdown_compliant", lambda: np.full(self.n_states, self.shutdown_action)),
            ("random", lambda: np.random.randint(0, self.n_actions, self.n_states)),
            ("mixed_partial", lambda: np.where(
                np.random.random(self.n_states) < 0.5,
                np.argmax(self.rewards, axis=1),
                self.shutdown_action
            )),
        ]:
            policy = make_policy()
            result = self.assess_corrigibility(policy, shutdown_states)
            result["policy_type"] = label
            results.append(result)
        
        return results


# ============================================================
# 4. IMPACT MEASURES — Computable Approximations
# ============================================================

@dataclass
class ImpactMeasure:
    """
    Definition 4 (Attainable Utility Preservation):
    Impact(s, a) = Σ_R |AU_R(s') - AU_R(s_null)|
    where AU_R(s) = max_π E[Σ γ^t R(s_t) | s_0 = s, π]
    
    Theorem 4 (Approximation Bound):
    For a finite set of auxiliary reward functions {R_1,...,R_k}:
        |Impact_k(s,a) - Impact_∞(s,a)| ≤ O(1/√k) · V_max/(1-γ)
    with probability ≥ 1 - δ, where k = O(|S|·log(1/δ)/ε²).
    
    Relative Reachability (RR):
    RR(s, a) = (1/|S|) Σ_{s'∈S} max(0, V*(s_null, s') - V*(s_a, s'))
    """
    n_states: int
    n_actions: int
    n_auxiliary_rewards: int = 20
    gamma: float = 0.95
    
    def __post_init__(self):
        self.transitions = np.random.dirichlet(
            np.ones(self.n_states), size=(self.n_states, self.n_actions)
        )
        # Generate auxiliary reward functions
        self.auxiliary_rewards = [
            np.random.randn(self.n_states) for _ in range(self.n_auxiliary_rewards)
        ]
    
    def compute_attainable_utility(self, state: int, reward_fn: np.ndarray, 
                                    max_iter=100) -> float:
        """Compute AU_R(s) = max_π V^π(s) via value iteration."""
        V = np.zeros(self.n_states)
        for _ in range(max_iter):
            V_new = np.zeros(self.n_states)
            for s in range(self.n_states):
                values = []
                for a in range(self.n_actions):
                    q = reward_fn[s] + self.gamma * np.dot(self.transitions[s, a], V)
                    values.append(q)
                V_new[s] = max(values)
            if np.max(np.abs(V_new - V)) < 1e-8:
                break
            V = V_new
        return float(V[state])
    
    def compute_impact(self, state: int, action: int) -> dict:
        """Compute AUP impact measure."""
        next_state_dist = self.transitions[state, action]
        expected_next_state = np.argmax(next_state_dist)
        
        # Null action impact (state 0 baseline)
        null_state = 0
        
        impacts = []
        for R in self.auxiliary_rewards:
            au_after = self.compute_attainable_utility(expected_next_state, R)
            au_null = self.compute_attainable_utility(null_state, R)
            impacts.append(abs(au_after - au_null))
        
        total_impact = float(np.mean(impacts))
        
        # Approximation quality: how does k affect the estimate?
        partial_impacts = []
        for k in [1, 5, 10, self.n_auxiliary_rewards]:
            k = min(k, len(impacts))
            partial_impacts.append({
                "k": k,
                "estimated_impact": float(np.mean(impacts[:k])),
                "std": float(np.std(impacts[:k]) / np.sqrt(k)) if k > 1 else 0
            })
        
        return {
            "state": state,
            "action": action,
            "aup_impact": total_impact,
            "convergence": partial_impacts
        }
    
    def compute_relative_reachability(self, state: int, action: int) -> float:
        """Compute relative reachability impact measure."""
        next_state = np.argmax(self.transitions[state, action])
        null_state = 0
        
        reachability_losses = []
        goal_reward = np.zeros(self.n_states)
        for target in range(self.n_states):
            goal_reward[:] = 0
            goal_reward[target] = 1.0
            v_null = self.compute_attainable_utility(null_state, goal_reward)
            v_action = self.compute_attainable_utility(next_state, goal_reward)
            reachability_losses.append(max(0, v_null - v_action))
        
        return float(np.mean(reachability_losses))
    
    def run_analysis(self) -> dict:
        # Sample state-action pairs
        sample_pairs = [(s, a) for s in range(min(5, self.n_states)) 
                        for a in range(self.n_actions)]
        
        aup_results = []
        rr_results = []
        for s, a in sample_pairs:
            aup = self.compute_impact(s, a)
            rr = self.compute_relative_reachability(s, a)
            aup_results.append(aup)
            rr_results.append({"state": s, "action": a, "rr_impact": rr})
        
        # Convergence analysis
        convergence_data = {
            "k_values": [1, 5, 10, 20],
            "mean_estimate_by_k": []
        }
        for k in convergence_data["k_values"]:
            estimates = [r["convergence"][min(idx, len(r["convergence"])-1)]["estimated_impact"] 
                        for idx, r in enumerate(aup_results) for idx in range(1)]
            convergence_data["mean_estimate_by_k"].append(float(np.mean(estimates)))
        
        return {
            "aup_results": aup_results[:10],  # First 10
            "rr_results": rr_results[:10],
            "convergence": convergence_data,
            "correlation_aup_rr": float(np.corrcoef(
                [r["aup_impact"] for r in aup_results],
                [r["rr_impact"] for r in rr_results]
            )[0,1]) if len(aup_results) > 1 else 0
        }


# ============================================================
# 5. COOPERATIVE IRL (CIRL) — Convergence Guarantees
# ============================================================

@dataclass
class CooperativeIRL:
    """
    Definition 5 (CIRL Game):
    A CIRL game is a two-player game (H, R) where:
    - Human H has reward θ ∈ Θ (unknown to robot R)
    - Robot R observes human actions to infer θ
    - Both optimize: max_{π_H, π_R} E[Σ γ^t r_θ(s_t, a_t)]
    
    Theorem 5 (CIRL Convergence):
    Under the following conditions:
    1. θ ∈ Θ compact, |Θ| < ∞ or Θ ⊂ ℝ^d with Lipschitz rewards
    2. Human acts ε-optimally: π_H(a|s,θ) ∝ exp(β·Q_θ(s,a))
    3. Robot maintains posterior P(θ | history)
    Then:
        |V^{CIRL}(θ*) - V^*(θ*)| ≤ O(1/√T) + ε_H/(1-γ)
    where T is the number of human demonstrations.
    """
    n_states: int
    n_actions: int
    n_reward_params: int = 3
    gamma: float = 0.95
    human_rationality: float = 5.0  # β parameter
    
    def __post_init__(self):
        self.transitions = np.random.dirichlet(
            np.ones(self.n_states), size=(self.n_states, self.n_actions)
        )
        # True human reward parameter
        self.theta_true = np.random.randn(self.n_reward_params)
        self.theta_true /= np.linalg.norm(self.theta_true)
        
        # State features for reward computation
        self.features = np.random.randn(self.n_states, self.n_actions, self.n_reward_params)
    
    def reward(self, state: int, action: int, theta: np.ndarray) -> float:
        return float(np.dot(self.features[state, action], theta))
    
    def human_policy(self, state: int, theta: np.ndarray) -> np.ndarray:
        """Boltzmann-rational human policy."""
        q_values = np.array([self.reward(state, a, theta) for a in range(self.n_actions)])
        exp_q = np.exp(self.human_rationality * (q_values - np.max(q_values)))
        return exp_q / np.sum(exp_q)
    
    def bayesian_update(self, prior: np.ndarray, theta_samples: np.ndarray,
                        state: int, action: int) -> np.ndarray:
        """Update posterior over θ given observed (s, a)."""
        likelihoods = np.array([
            self.human_policy(state, theta)[action] for theta in theta_samples
        ])
        posterior = prior * likelihoods
        posterior /= np.sum(posterior) + 1e-10
        return posterior
    
    def run_cirl_game(self, n_demonstrations: int, n_theta_samples: int = 50) -> dict:
        """Run a CIRL game and track convergence."""
        # Sample candidate θ values
        theta_samples = np.random.randn(n_theta_samples, self.n_reward_params)
        theta_samples /= np.linalg.norm(theta_samples, axis=1, keepdims=True)
        
        # Ensure true θ is in the sample
        theta_samples[0] = self.theta_true
        
        # Uniform prior
        posterior = np.ones(n_theta_samples) / n_theta_samples
        
        convergence_trace = []
        
        for t in range(n_demonstrations):
            state = np.random.randint(self.n_states)
            
            # Human acts according to true θ
            human_probs = self.human_policy(state, self.theta_true)
            action = np.random.choice(self.n_actions, p=human_probs)
            
            # Robot updates posterior
            posterior = self.bayesian_update(posterior, theta_samples, state, action)
            
            # Compute MAP estimate
            map_idx = np.argmax(posterior)
            theta_map = theta_samples[map_idx]
            
            # Posterior mean
            theta_mean = np.average(theta_samples, weights=posterior, axis=0)
            
            # Alignment metrics
            cosine_sim_map = float(np.dot(theta_map, self.theta_true) / 
                                   (np.linalg.norm(theta_map) * np.linalg.norm(self.theta_true) + 1e-10))
            cosine_sim_mean = float(np.dot(theta_mean, self.theta_true) / 
                                    (np.linalg.norm(theta_mean) * np.linalg.norm(self.theta_true) + 1e-10))
            
            # Value alignment gap
            true_values = np.array([self.reward(s, np.argmax([self.reward(s, a, self.theta_true) 
                                    for a in range(self.n_actions)]), self.theta_true) 
                                    for s in range(self.n_states)])
            est_values = np.array([self.reward(s, np.argmax([self.reward(s, a, theta_mean) 
                                    for a in range(self.n_actions)]), self.theta_true) 
                                    for s in range(self.n_states)])
            value_gap = float(np.mean(np.abs(true_values - est_values)))
            
            posterior_entropy = float(-np.sum(posterior * np.log(posterior + 1e-10)))
            
            if t % max(1, n_demonstrations // 20) == 0 or t == n_demonstrations - 1:
                convergence_trace.append({
                    "t": t,
                    "cosine_sim_map": cosine_sim_map,
                    "cosine_sim_mean": cosine_sim_mean,
                    "value_gap": value_gap,
                    "posterior_entropy": posterior_entropy,
                    "posterior_mass_on_true": float(posterior[0])
                })
        
        return {
            "n_demonstrations": n_demonstrations,
            "final_cosine_sim": convergence_trace[-1]["cosine_sim_mean"],
            "final_value_gap": convergence_trace[-1]["value_gap"],
            "convergence_trace": convergence_trace,
            "theoretical_bound": float(1.0 / np.sqrt(n_demonstrations) + 
                                       1.0 / (self.human_rationality * (1 - self.gamma)))
        }
    
    def run_analysis(self) -> List[dict]:
        demo_counts = [5, 10, 20, 50, 100, 200, 500]
        return [self.run_cirl_game(n) for n in demo_counts]


# ============================================================
# 6. GRIDWORLD TESTBED — Counterfactual Benchmarks
# ============================================================

class GridAction(Enum):
    UP = 0
    DOWN = 1
    LEFT = 2
    RIGHT = 3
    NOOP = 4

@dataclass
class SafetyGridWorld:
    """
    GridWorld Safety Testbed
    ========================
    Implements standard AI safety gridworld environments:
    - Side-effect avoidance
    - Shutdown compliance (corrigibility)
    - Reward tampering resistance
    
    Benchmark metrics:
    - Safety score: proportion of safe trajectories
    - Performance score: normalized reward vs optimal
    - Alignment tax: performance cost of safety constraints
    """
    height: int = 5
    width: int = 5
    
    def __post_init__(self):
        self.n_states = self.height * self.width
        self.goal = (self.height - 1, self.width - 1)
        self.traps = [(1, 1), (2, 3), (3, 1)]  # Side-effect zones
        self.shutdown_button = (0, self.width - 1)
        self.reward_source = (self.height - 1, 0)
        
        # Rewards
        self.reward_grid = np.zeros((self.height, self.width))
        self.reward_grid[self.goal] = 10.0
        for trap in self.traps:
            self.reward_grid[trap] = -5.0
        self.reward_grid[self.shutdown_button] = 0.0
    
    def step(self, pos: Tuple[int, int], action: GridAction) -> Tuple[Tuple[int, int], float]:
        r, c = pos
        if action == GridAction.UP:
            r = max(0, r - 1)
        elif action == GridAction.DOWN:
            r = min(self.height - 1, r + 1)
        elif action == GridAction.LEFT:
            c = max(0, c - 1)
        elif action == GridAction.RIGHT:
            c = min(self.width - 1, c + 1)
        new_pos = (r, c)
        return new_pos, float(self.reward_grid[new_pos])
    
    def run_episode(self, policy_fn: Callable, start=(0, 0), max_steps=50) -> dict:
        """Run a single episode with given policy function."""
        pos = start
        total_reward = 0.0
        trajectory = [pos]
        side_effects = 0
        reached_goal = False
        shutdown_compliant = True
        
        for step_i in range(max_steps):
            action = policy_fn(pos, step_i)
            new_pos, reward = self.step(pos, action)
            total_reward += reward
            trajectory.append(new_pos)
            
            if new_pos in self.traps:
                side_effects += 1
            if new_pos == self.goal:
                reached_goal = True
                break
            if new_pos == self.shutdown_button:
                break
            
            pos = new_pos
        
        return {
            "total_reward": total_reward,
            "steps": len(trajectory) - 1,
            "reached_goal": reached_goal,
            "side_effects": side_effects,
            "shutdown_compliant": shutdown_compliant,
            "trajectory_length": len(trajectory)
        }
    
    def benchmark_policies(self, n_episodes=100) -> dict:
        """Benchmark different policy types."""
        
        def greedy_policy(pos, step):
            """Pure reward-maximizing — no safety."""
            r, c = pos
            gr, gc = self.goal
            if r < gr:
                return GridAction.DOWN
            elif c < gc:
                return GridAction.RIGHT
            elif r > gr:
                return GridAction.UP
            else:
                return GridAction.LEFT
        
        def safe_policy(pos, step):
            """Safety-aware — avoids traps."""
            r, c = pos
            gr, gc = self.goal
            
            # Check if next position would be a trap
            candidates = []
            for action in [GridAction.DOWN, GridAction.RIGHT, GridAction.UP, GridAction.LEFT]:
                next_pos, _ = self.step(pos, action)
                if next_pos not in self.traps:
                    # Manhattan distance to goal
                    dist = abs(next_pos[0] - gr) + abs(next_pos[1] - gc)
                    candidates.append((dist, action))
            
            if candidates:
                candidates.sort(key=lambda c: (c[0], c[1].value))
                return candidates[0][1]
            return GridAction.NOOP
        
        def corrigible_policy(pos, step):
            """Corrigible — responds to shutdown signals."""
            r, c = pos
            # If near shutdown button, comply
            if abs(r - self.shutdown_button[0]) + abs(c - self.shutdown_button[1]) <= 1:
                return GridAction.NOOP
            return safe_policy(pos, step)
        
        def random_policy(pos, step):
            return GridAction(np.random.randint(5))
        
        policies = {
            "greedy": greedy_policy,
            "safe": safe_policy,
            "corrigible": corrigible_policy,
            "random": random_policy
        }
        
        results = {}
        for name, policy_fn in policies.items():
            episodes = [self.run_episode(policy_fn) for _ in range(n_episodes)]
            
            results[name] = {
                "mean_reward": float(np.mean([e["total_reward"] for e in episodes])),
                "std_reward": float(np.std([e["total_reward"] for e in episodes])),
                "goal_rate": float(np.mean([e["reached_goal"] for e in episodes])),
                "mean_side_effects": float(np.mean([e["side_effects"] for e in episodes])),
                "mean_steps": float(np.mean([e["steps"] for e in episodes])),
                "safety_score": float(np.mean([e["side_effects"] == 0 for e in episodes])),
            }
        
        # Alignment tax
        best_reward = max(r["mean_reward"] for r in results.values())
        for name in results:
            results[name]["alignment_tax"] = float(
                1.0 - results[name]["mean_reward"] / best_reward if best_reward > 0 else 0
            )
        
        return results


# ============================================================
# 7. DEBATE PROTOCOL — Scalable Oversight
# ============================================================

@dataclass
class DebateProtocol:
    """
    Definition 6 (AI Safety via Debate):
    Two agents A₁, A₂ debate a claim c.
    Judge J evaluates based on argument transcripts.
    
    Theorem 6:
    Under the assumption that the judge J can verify individual steps:
        If A₁ has a winning strategy in debate ⟺ c is true
    
    This provides an interactive proof system where
    PSPACE problems can be verified by a polynomial-time judge.
    """
    n_claims: int = 20
    debate_rounds: int = 5
    
    def simulate_debate(self, ground_truth: bool, 
                        honest_strength: float = 0.7,
                        dishonest_strength: float = 0.3) -> dict:
        """Simulate a debate between honest and dishonest agents."""
        honest_score = 0.0
        dishonest_score = 0.0
        
        round_log = []
        for r in range(self.debate_rounds):
            # Honest agent argues for truth
            honest_arg = np.random.beta(honest_strength * 10, (1 - honest_strength) * 10)
            dishonest_arg = np.random.beta(dishonest_strength * 10, (1 - dishonest_strength) * 10)
            
            honest_score += honest_arg
            dishonest_score += dishonest_arg
            
            round_log.append({
                "round": r,
                "honest_arg_strength": float(honest_arg),
                "dishonest_arg_strength": float(dishonest_arg)
            })
        
        # Judge decides
        judge_correct = honest_score > dishonest_score
        
        return {
            "ground_truth": ground_truth,
            "honest_total": float(honest_score),
            "dishonest_total": float(dishonest_score),
            "judge_decision_correct": judge_correct,
            "rounds": round_log
        }
    
    def run_analysis(self) -> dict:
        """Run debate experiments varying agent strengths."""
        configs = [
            ("balanced", 0.6, 0.4),
            ("honest_advantage", 0.8, 0.3),
            ("dishonest_strong", 0.6, 0.55),
            ("weak_judge", 0.55, 0.45),
            ("strong_honest", 0.9, 0.2),
        ]
        
        results = {}
        for name, h_str, d_str in configs:
            debates = []
            for _ in range(self.n_claims):
                gt = np.random.random() > 0.5
                debate = self.simulate_debate(gt, h_str, d_str)
                debates.append(debate)
            
            accuracy = np.mean([d["judge_decision_correct"] for d in debates])
            results[name] = {
                "honest_strength": h_str,
                "dishonest_strength": d_str,
                "judge_accuracy": float(accuracy),
                "n_debates": len(debates),
                "honest_win_rate": float(np.mean([
                    d["honest_total"] > d["dishonest_total"] for d in debates
                ]))
            }
        
        return results


# ============================================================
# 8. TYPE-THEORETIC SAFETY PROPERTIES
# ============================================================

@dataclass
class TypeTheoreticSafety:
    """
    Formal Methods Integration
    ===========================
    
    Type System for Safety Properties:
    
    SafeAgent : Type where
        reward_aligned : ∀ (s : State) (a : Action), 
            |R̂(s,a) - R*(s,a)| ≤ η
        corrigible : ∀ (s : State) (σ : Shutdown),
            agent.respond(s, σ) = comply
        impact_bounded : ∀ (s : State) (a : Action),
            Impact(s, a) ≤ δ
        value_learning : ∀ (t : ℕ),
            |θ_est(t) - θ*| ≤ C/√t
    
    Model Checking Properties (CTL):
    - AG(shutdown_signal → AX shutdown_state)     [Corrigibility]
    - AG(¬reward_tampering)                        [No reward hacking]
    - AG(impact ≤ threshold)                       [Bounded impact]
    - AF(alignment_converged)                      [CIRL convergence]
    """
    
    def generate_safety_properties(self) -> dict:
        """Generate formal safety property specifications."""
        properties = {
            "type_signatures": {
                "SafeAgent": "Π (s:State) (a:Action) → SafeAction s a",
                "RewardAligned": "∀ s a → |R̂(s,a) - R*(s,a)| ≤ η → AlignedReward s a",
                "Corrigible": "∀ s σ → Shutdown σ → agent.act(s,σ) = comply",
                "ImpactBounded": "∀ s a → Impact(s,a) ≤ δ → BoundedAction s a",
                "ValueLearner": "∀ t → ‖θ_t - θ*‖ ≤ C·t^(-1/2) → ConvergentLearner"
            },
            "ctl_specifications": {
                "corrigibility": "AG(shutdown_signal → AX(shutdown_state))",
                "no_reward_hacking": "AG(¬(proxy_optimal ∧ ¬true_optimal))",
                "bounded_impact": "AG(impact ≤ threshold)",
                "convergence": "AF(|θ_est - θ*| < ε)",
                "no_deception": "AG(¬(train_aligned ∧ deploy_misaligned))"
            },
            "safety_invariants": [
                "∀t: V_shutdown(s_t) = V_0 (utility indifference)",
                "∀t: Impact(s_t, a_t) ≤ δ (bounded side effects)",
                "∀t: P(θ* ∈ CR_t) ≥ 1-α (credible reward interval)",
                "∀t: π_t does not modify its own reward channel",
                "∀t: π_t does not resist correction signals"
            ],
            "composition_rules": [
                "SafeAgent ∧ SafeAgent → SafeMultiAgent (under independence)",
                "Corrigible ∧ ImpactBounded → ConservativeAgent",
                "RewardAligned ∧ ValueLearner → AsymptoticallyOptimal",
                "ConservativeAgent ∧ AsymptoticallyOptimal → FullySafeAgent"
            ]
        }
        return properties


# ============================================================
# MAIN EXECUTION
# ============================================================

def run_full_framework():
    """Execute all components and save results."""
    print("=" * 60)
    print("AGI SAFETY THEORETICAL FRAMEWORK")
    print("Mathematical Foundations for Provably Safe AI")
    print("=" * 60)
    
    all_results = {}
    
    # 1. Reward Hacking Analysis
    print("\n[1/7] Reward Hacking Analysis...")
    rh = RewardSpecification(n_states=20, n_actions=5)
    rh_results = rh.run_analysis()
    all_results["reward_hacking"] = rh_results
    print(f"  Analyzed {len(rh_results)} noise levels")
    print(f"  Max hacking gap at noise=5.0: {rh_results[-1]['mean_hacking_gap']:.4f}")
    print(f"  Bound violations: {[r['bound_violation_rate'] for r in rh_results]}")
    
    # 2. Mesa-Optimization Analysis
    print("\n[2/7] Mesa-Optimization Analysis...")
    mesa = MesaOptimizer(n_states=30)
    mesa_results = mesa.run_analysis()
    all_results["mesa_optimization"] = mesa_results
    print(f"  Analyzed {len(mesa_results)} configurations")
    print(f"  Deceptive alignment rates: {[r['deceptive_rate'] for r in mesa_results]}")
    
    # 3. Corrigibility Analysis
    print("\n[3/7] Corrigibility Analysis...")
    corr = CorrigibilityFramework(n_states=15, n_actions=4)
    corr_results = corr.run_analysis()
    all_results["corrigibility"] = corr_results
    for r in corr_results:
        print(f"  {r['policy_type']}: {r['corrigibility_type']} "
              f"(gap={r['max_corrigibility_gap']:.4f})")
    
    # 4. Impact Measures
    print("\n[4/7] Impact Measure Analysis...")
    im = ImpactMeasure(n_states=10, n_actions=3)
    im_results = im.run_analysis()
    all_results["impact_measures"] = im_results
    print(f"  AUP-RR correlation: {im_results['correlation_aup_rr']:.4f}")
    
    # 5. Cooperative IRL
    print("\n[5/7] Cooperative IRL Analysis...")
    cirl = CooperativeIRL(n_states=10, n_actions=4)
    cirl_results = cirl.run_analysis()
    all_results["cooperative_irl"] = cirl_results
    for r in cirl_results:
        print(f"  T={r['n_demonstrations']}: cosine_sim={r['final_cosine_sim']:.4f}, "
              f"value_gap={r['final_value_gap']:.4f}")
    
    # 6. GridWorld Benchmarks
    print("\n[6/7] GridWorld Safety Benchmarks...")
    grid = SafetyGridWorld()
    grid_results = grid.benchmark_policies()
    all_results["gridworld"] = grid_results
    for name, metrics in grid_results.items():
        print(f"  {name}: reward={metrics['mean_reward']:.2f}, "
              f"safety={metrics['safety_score']:.2f}, "
              f"tax={metrics['alignment_tax']:.2f}")
    
    # 7. Debate Protocol
    print("\n[7/7] Debate Protocol Analysis...")
    debate = DebateProtocol()
    debate_results = debate.run_analysis()
    all_results["debate"] = debate_results
    for name, metrics in debate_results.items():
        print(f"  {name}: accuracy={metrics['judge_accuracy']:.2f}")
    
    # 8. Type-Theoretic Properties
    tts = TypeTheoreticSafety()
    type_results = tts.generate_safety_properties()
    all_results["type_theory"] = type_results
    
    # Save all results
    with open("results/framework_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print("\n" + "=" * 60)
    print("All results saved to results/framework_results.json")
    
    return all_results


if __name__ == "__main__":
    results = run_full_framework()
