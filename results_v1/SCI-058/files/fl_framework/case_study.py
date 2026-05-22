"""
Case study: Federated survival analysis on multi-site clinical data.

Simulates a multi-hospital federated learning scenario for
Cox proportional hazards model with:
  - Non-IID data distribution across hospitals
  - Differential privacy
  - Communication compression
  - Byzantine resilience evaluation
"""

import numpy as np
import json
import os
import sys
import time
from typing import Dict, List, Tuple

# Add parent to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fl_framework.aggregation import FedAvg, FedProx, SCAFFOLD, ClientUpdate
from fl_framework.differential_privacy import (
    DPFederatedAggregator,
    RDPAccountant,
    GradientClipper,
)
from fl_framework.communication import TopKSparsifier, StochasticQuantizer
from fl_framework.byzantine import (
    Krum,
    CoordinateMedian,
    TrimmedMean,
    ByzantineDetector,
)


def set_all_seeds(seed: int = 42):
    np.random.seed(seed)
    try:
        import random
        random.seed(seed)
    except ImportError:
        pass


def generate_hospital_data(
    num_hospitals: int = 5,
    samples_per_hospital: List[int] = None,
    num_features: int = 20,
    non_iid_degree: float = 0.5,
    seed: int = 42,
) -> List[Dict]:
    """
    Generate synthetic clinical survival data with non-IID distribution.

    Each hospital has different:
      - Sample sizes (institutional variation)
      - Feature distributions (demographic differences)
      - Event rates (treatment protocol differences)
    """
    set_all_seeds(seed)

    if samples_per_hospital is None:
        samples_per_hospital = [200, 350, 150, 500, 300]

    hospitals = []
    # Base coefficients for Cox model
    true_beta = np.random.randn(num_features) * 0.5

    for h in range(num_hospitals):
        n = samples_per_hospital[h] if h < len(samples_per_hospital) else 200

        # Non-IID: shift feature means per hospital
        feature_shift = np.random.randn(num_features) * non_iid_degree * (h + 1) / num_hospitals
        X = np.random.randn(n, num_features) + feature_shift

        # Generate survival times using Cox model
        linear_pred = X @ true_beta
        baseline_hazard = 0.1 * (1 + 0.2 * h)  # Different baseline per hospital
        hazard = baseline_hazard * np.exp(linear_pred)
        survival_time = np.random.exponential(1.0 / (hazard + 1e-6))

        # Censoring (varies by hospital)
        censor_rate = 0.3 + 0.1 * h / num_hospitals
        censor_time = np.random.exponential(
            np.median(survival_time) / (1 - censor_rate + 1e-6), size=n
        )
        observed_time = np.minimum(survival_time, censor_time)
        event = (survival_time <= censor_time).astype(float)

        hospitals.append({
            "hospital_id": f"hospital_{h}",
            "X": X,
            "time": observed_time,
            "event": event,
            "num_samples": n,
            "event_rate": float(np.mean(event)),
            "median_time": float(np.median(observed_time)),
        })

    return hospitals, true_beta


def cox_gradient(
    X: np.ndarray,
    time: np.ndarray,
    event: np.ndarray,
    beta: np.ndarray,
) -> Tuple[np.ndarray, float]:
    """Compute gradient of negative Cox partial log-likelihood."""
    n = X.shape[0]
    risk_scores = X @ beta
    exp_risk = np.exp(risk_scores - np.max(risk_scores))

    # Sort by time (descending)
    sort_idx = np.argsort(-time)
    X_sorted = X[sort_idx]
    event_sorted = event[sort_idx]
    exp_risk_sorted = exp_risk[sort_idx]

    # Compute gradient
    gradient = np.zeros_like(beta)
    cumsum_exp = 0.0
    cumsum_x_exp = np.zeros_like(beta)

    loss = 0.0
    for i in range(n):
        cumsum_exp += exp_risk_sorted[i]
        cumsum_x_exp += X_sorted[i] * exp_risk_sorted[i]

        if event_sorted[i] == 1:
            gradient += -(X_sorted[i] - cumsum_x_exp / (cumsum_exp + 1e-10))
            loss += -(risk_scores[sort_idx[i]] - np.log(cumsum_exp + 1e-10))

    gradient /= max(1, np.sum(event))
    loss /= max(1, np.sum(event))

    return gradient, float(loss)


def concordance_index(
    risk_scores: np.ndarray,
    time: np.ndarray,
    event: np.ndarray,
) -> float:
    """Compute Harrell's concordance index."""
    n = len(time)
    concordant = 0
    discordant = 0
    tied = 0

    for i in range(n):
        if event[i] == 0:
            continue
        for j in range(n):
            if i == j:
                continue
            if time[j] > time[i]:
                if risk_scores[i] > risk_scores[j]:
                    concordant += 1
                elif risk_scores[i] < risk_scores[j]:
                    discordant += 1
                else:
                    tied += 1

    total = concordant + discordant + tied
    if total == 0:
        return 0.5
    return (concordant + 0.5 * tied) / total


def local_training_step(
    X: np.ndarray,
    time: np.ndarray,
    event: np.ndarray,
    global_beta: np.ndarray,
    local_epochs: int = 5,
    lr: float = 0.01,
    mu: float = 0.0,  # FedProx proximal term
) -> Tuple[np.ndarray, float]:
    """Simulate local training at one hospital."""
    beta = global_beta.copy()
    total_loss = 0.0

    for _ in range(local_epochs):
        grad, loss = cox_gradient(X, time, event, beta)
        # Add proximal term gradient if mu > 0 (FedProx)
        if mu > 0:
            grad += mu * (beta - global_beta)
        beta -= lr * grad
        total_loss = loss

    return beta, total_loss


def run_federated_experiment(
    hospitals: List[Dict],
    true_beta: np.ndarray,
    strategy: str = "fedavg",
    num_rounds: int = 50,
    local_epochs: int = 5,
    local_lr: float = 0.01,
    dp_epsilon: float = None,
    dp_noise_multiplier: float = 1.0,
    compression_ratio: float = None,
    byzantine_fraction: float = 0.0,
    byzantine_defense: str = None,
    mu: float = 0.01,
    seed: int = 42,
) -> Dict:
    """Run a complete federated learning experiment."""
    set_all_seeds(seed)
    num_features = hospitals[0]["X"].shape[1]
    num_clients = len(hospitals)

    # Initialize global model
    global_beta = np.zeros(num_features)

    # Setup DP if requested
    dp_aggregator = None
    if dp_epsilon is not None:
        dp_aggregator = DPFederatedAggregator(
            epsilon=dp_epsilon,
            delta=1e-5,
            max_grad_norm=1.0,
            noise_multiplier=dp_noise_multiplier,
            num_clients=num_clients,
            clients_per_round=num_clients,
        )

    # Setup compression
    sparsifier = None
    if compression_ratio is not None:
        sparsifier = TopKSparsifier(compression_ratio=compression_ratio)

    # Setup Byzantine
    num_byzantine = int(num_clients * byzantine_fraction)
    detector = ByzantineDetector(z_threshold=2.5) if byzantine_defense else None

    # Aggregation strategy
    if strategy == "fedavg":
        aggregator = FedAvg()
    elif strategy == "fedprox":
        aggregator = FedProx(mu=mu)
    elif strategy == "scaffold":
        aggregator = SCAFFOLD(num_clients=num_clients)
        aggregator.initialize_controls({"beta": global_beta})
    else:
        aggregator = FedAvg()

    # Training history
    history = {
        "round": [],
        "global_loss": [],
        "c_index": [],
        "beta_error": [],
        "privacy_spent": [],
        "compression_ratio": [],
        "byzantine_detected": [],
    }

    for rnd in range(num_rounds):
        client_updates = []
        client_deltas = []

        for h_idx, hospital in enumerate(hospitals):
            # Byzantine attack: random noise injection
            if h_idx < num_byzantine:
                fake_beta = global_beta + np.random.randn(num_features) * 10.0
                local_loss = 999.0
            else:
                effective_mu = mu if strategy == "fedprox" else 0.0
                fake_beta, local_loss = local_training_step(
                    hospital["X"],
                    hospital["time"],
                    hospital["event"],
                    global_beta,
                    local_epochs=local_epochs,
                    lr=local_lr,
                    mu=effective_mu,
                )

            delta = {"beta": fake_beta - global_beta}

            # Compression
            if sparsifier is not None:
                delta, _ = sparsifier.compress(delta)

            client_deltas.append(delta)
            client_updates.append(
                ClientUpdate(
                    client_id=str(h_idx),
                    weights={"beta": fake_beta},
                    num_samples=hospital["num_samples"],
                    loss=local_loss,
                )
            )

        # Byzantine detection and filtering
        byz_detected = []
        if detector is not None:
            byz_detected = detector.detect(client_deltas)
            if byz_detected:
                client_deltas = [d for i, d in enumerate(client_deltas) if i not in byz_detected]
                client_updates = [u for i, u in enumerate(client_updates) if i not in byz_detected]

        # Byzantine-resilient aggregation
        if byzantine_defense == "krum" and len(client_updates) > 2:
            krum = Krum(num_byzantine=num_byzantine, multi_krum_k=max(1, len(client_deltas) - num_byzantine))
            new_weights = krum.aggregate({"beta": global_beta}, [d for d in client_deltas])
            global_beta = new_weights["beta"]
        elif byzantine_defense == "median":
            median_agg = CoordinateMedian()
            new_weights = median_agg.aggregate({"beta": global_beta}, client_deltas)
            global_beta = new_weights["beta"]
        elif byzantine_defense == "trimmed_mean":
            tm_agg = TrimmedMean(trim_ratio=0.2)
            new_weights = tm_agg.aggregate({"beta": global_beta}, client_deltas)
            global_beta = new_weights["beta"]
        elif dp_aggregator is not None:
            # DP aggregation
            new_weights = dp_aggregator.aggregate_with_dp(
                {"beta": global_beta}, client_deltas
            )
            global_beta = new_weights["beta"]
        elif strategy == "scaffold":
            new_weights = aggregator.aggregate(
                {"beta": global_beta}, client_updates,
                local_steps=local_epochs, local_lr=local_lr,
            )
            global_beta = new_weights["beta"]
        else:
            # Standard aggregation
            new_weights = aggregator.aggregate({"beta": global_beta}, client_updates)
            global_beta = new_weights["beta"]

        # Evaluate
        all_X = np.vstack([h["X"] for h in hospitals])
        all_time = np.concatenate([h["time"] for h in hospitals])
        all_event = np.concatenate([h["event"] for h in hospitals])

        risk_scores = all_X @ global_beta
        c_index = concordance_index(risk_scores, all_time, all_event)
        _, global_loss = cox_gradient(all_X, all_time, all_event, global_beta)
        beta_error = np.linalg.norm(global_beta - true_beta) / np.linalg.norm(true_beta)

        # Privacy accounting
        eps_spent = 0.0
        if dp_aggregator is not None:
            report = dp_aggregator.get_privacy_report()
            eps_spent = report["epsilon_spent"]

        # Compression ratio
        comp_ratio = compression_ratio if compression_ratio else 1.0

        history["round"].append(rnd)
        history["global_loss"].append(float(global_loss))
        history["c_index"].append(float(c_index))
        history["beta_error"].append(float(beta_error))
        history["privacy_spent"].append(float(eps_spent))
        history["compression_ratio"].append(float(comp_ratio))
        history["byzantine_detected"].append(len(byz_detected))

    return {
        "strategy": strategy,
        "final_c_index": history["c_index"][-1],
        "final_loss": history["global_loss"][-1],
        "final_beta_error": history["beta_error"][-1],
        "final_privacy_spent": history["privacy_spent"][-1],
        "history": history,
        "global_beta": global_beta.tolist(),
    }


def run_all_experiments():
    """Run the complete experimental suite."""
    print("=" * 70)
    print("Federated Survival Analysis: Multi-Site Clinical Case Study")
    print("=" * 70)

    # Generate data
    print("\n[1/6] Generating synthetic multi-hospital data...")
    hospitals, true_beta = generate_hospital_data(
        num_hospitals=5,
        samples_per_hospital=[200, 350, 150, 500, 300],
        num_features=20,
        non_iid_degree=0.5,
        seed=42,
    )

    for h in hospitals:
        print(f"  {h['hospital_id']}: n={h['num_samples']}, "
              f"event_rate={h['event_rate']:.2f}, "
              f"median_time={h['median_time']:.2f}")

    # Save data summary
    data_summary = []
    for h in hospitals:
        data_summary.append({
            "hospital_id": h["hospital_id"],
            "num_samples": h["num_samples"],
            "event_rate": round(h["event_rate"], 3),
            "median_time": round(h["median_time"], 3),
            "feature_means": h["X"].mean(axis=0)[:5].tolist(),
        })

    os.makedirs("data", exist_ok=True)
    with open("data/hospital_data_summary.json", "w") as f:
        json.dump(data_summary, f, indent=2)

    results = {}

    # Experiment 1: Compare aggregation strategies (FedAvg vs FedProx vs SCAFFOLD)
    print("\n[2/6] Comparing aggregation strategies...")
    for strategy in ["fedavg", "fedprox", "scaffold"]:
        print(f"  Running {strategy}...")
        result = run_federated_experiment(
            hospitals, true_beta,
            strategy=strategy,
            num_rounds=50,
            local_epochs=5,
            local_lr=0.01,
            mu=0.01,
            seed=42,
        )
        results[f"strategy_{strategy}"] = result
        print(f"    C-index: {result['final_c_index']:.4f}, "
              f"Loss: {result['final_loss']:.4f}, "
              f"Beta error: {result['final_beta_error']:.4f}")

    # Experiment 2: Non-IID degree comparison
    print("\n[3/6] Evaluating non-IID impact...")
    for non_iid in [0.0, 0.5, 1.0, 2.0]:
        h_data, t_beta = generate_hospital_data(
            non_iid_degree=non_iid, seed=42
        )
        for strategy in ["fedavg", "fedprox"]:
            result = run_federated_experiment(
                h_data, t_beta,
                strategy=strategy,
                num_rounds=50,
                mu=0.01,
                seed=42,
            )
            results[f"noniid_{non_iid}_{strategy}"] = result
        print(f"  non-IID={non_iid}: "
              f"FedAvg C-index={results[f'noniid_{non_iid}_fedavg']['final_c_index']:.4f}, "
              f"FedProx C-index={results[f'noniid_{non_iid}_fedprox']['final_c_index']:.4f}")

    # Experiment 3: DP impact analysis
    print("\n[4/6] Analyzing differential privacy impact...")
    for eps in [1.0, 5.0, 10.0, 50.0]:
        noise_mult = 10.0 / eps  # Scale noise inversely with epsilon
        result = run_federated_experiment(
            hospitals, true_beta,
            strategy="fedavg",
            num_rounds=50,
            dp_epsilon=eps,
            dp_noise_multiplier=noise_mult,
            seed=42,
        )
        results[f"dp_eps_{eps}"] = result
        print(f"  epsilon={eps}: C-index={result['final_c_index']:.4f}, "
              f"privacy_spent={result['final_privacy_spent']:.4f}")

    # Experiment 4: Communication efficiency
    print("\n[5/6] Testing communication compression...")
    for comp_ratio in [0.01, 0.05, 0.1, 0.5, 1.0]:
        result = run_federated_experiment(
            hospitals, true_beta,
            strategy="fedavg",
            num_rounds=50,
            compression_ratio=comp_ratio if comp_ratio < 1.0 else None,
            seed=42,
        )
        results[f"compression_{comp_ratio}"] = result
        print(f"  compression={comp_ratio}: C-index={result['final_c_index']:.4f}")

    # Experiment 5: Byzantine resilience
    print("\n[6/6] Evaluating Byzantine resilience...")
    for byz_frac in [0.0, 0.2, 0.4]:
        for defense in [None, "krum", "median", "trimmed_mean"]:
            result = run_federated_experiment(
                hospitals, true_beta,
                strategy="fedavg",
                num_rounds=50,
                byzantine_fraction=byz_frac,
                byzantine_defense=defense,
                seed=42,
            )
            defense_name = defense if defense else "none"
            results[f"byzantine_{byz_frac}_{defense_name}"] = result
        byz_results = {
            d: results[f"byzantine_{byz_frac}_{d if d else 'none'}"]["final_c_index"]
            for d in [None, "krum", "median", "trimmed_mean"]
        }
        print(f"  Byzantine={byz_frac}: " +
              ", ".join(f"{k or 'none'}={v:.4f}" for k, v in byz_results.items()))

    # Save all results
    os.makedirs("results", exist_ok=True)

    # Save summary metrics
    summary = {}
    for key, res in results.items():
        summary[key] = {
            "final_c_index": round(res["final_c_index"], 4),
            "final_loss": round(res["final_loss"], 4),
            "final_beta_error": round(res["final_beta_error"], 4),
            "final_privacy_spent": round(res["final_privacy_spent"], 4),
        }

    with open("results/experiment_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Save detailed histories for key experiments
    for key in ["strategy_fedavg", "strategy_fedprox", "strategy_scaffold"]:
        with open(f"results/{key}_history.json", "w") as f:
            json.dump(results[key]["history"], f, indent=2)

    print("\n" + "=" * 70)
    print("All experiments completed. Results saved to results/")
    print("=" * 70)

    return results


if __name__ == "__main__":
    results = run_all_experiments()
