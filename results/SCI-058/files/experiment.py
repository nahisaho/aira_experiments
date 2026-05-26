"""
Federated Learning Framework for Privacy-Preserving Medical Data Analysis
==========================================================================
Comprehensive experiment covering:
1. FedAvg convergence and improvements
2. Non-IID data handling (FedProx, SCAFFOLD)
3. Differential Privacy integration
4. Communication efficiency (gradient compression, knowledge distillation)
5. Byzantine fault tolerance
6. Multi-site survival analysis case study
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score
from scipy.special import expit as sigmoid
import os
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
FIGURES_DIR = "figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================
# Utility: Simple Logistic Regression Model (for FL simulation)
# ============================================================
class LogisticRegression:
    def __init__(self, n_features, lr=0.01, l2=0.0):
        self.w = np.zeros(n_features)
        self.b = 0.0
        self.lr = lr
        self.l2 = l2

    def predict_proba(self, X):
        z = X @ self.w + self.b
        return sigmoid(z)

    def predict(self, X):
        return (self.predict_proba(X) >= 0.5).astype(int)

    def get_params(self):
        return {'w': self.w.copy(), 'b': self.b}

    def set_params(self, params):
        self.w = params['w'].copy()
        self.b = params['b']

    def compute_loss(self, X, y):
        p = np.clip(self.predict_proba(X), 1e-7, 1 - 1e-7)
        loss = -np.mean(y * np.log(p) + (1 - y) * np.log(1 - p))
        loss += 0.5 * self.l2 * np.sum(self.w ** 2)
        return loss

    def compute_gradients(self, X, y):
        p = self.predict_proba(X)
        error = p - y
        grad_w = X.T @ error / len(y) + self.l2 * self.w
        grad_b = np.mean(error)
        return {'w': grad_w, 'b': grad_b}

    def step(self, grads):
        self.w -= self.lr * grads['w']
        self.b -= self.lr * grads['b']


# ============================================================
# Cox Proportional Hazards Model for Survival Analysis
# ============================================================
class CoxPH:
    def __init__(self, n_features, lr=0.001):
        self.beta = np.zeros(n_features)
        self.lr = lr

    def get_params(self):
        return {'beta': self.beta.copy()}

    def set_params(self, params):
        self.beta = params['beta'].copy()

    def negative_partial_log_likelihood(self, X, times, events):
        risk_scores = X @ self.beta
        # Sort by time (descending for risk set)
        order = np.argsort(-times)
        X_sorted = X[order]
        events_sorted = events[order]
        risk_sorted = risk_scores[order]

        log_risk = np.log(np.cumsum(np.exp(risk_sorted)) + 1e-10)
        nll = -np.mean(events_sorted * (risk_sorted - log_risk))
        return nll

    def compute_gradients(self, X, times, events):
        risk_scores = X @ self.beta
        order = np.argsort(-times)
        X_s = X[order]
        events_s = events[order]
        rs = risk_scores[order]

        exp_rs = np.exp(rs)
        cum_exp = np.cumsum(exp_rs)
        cum_weighted = np.cumsum((exp_rs[:, None] * X_s), axis=0)

        grad = np.zeros_like(self.beta)
        for i in range(len(events_s)):
            if events_s[i] == 1:
                grad += X_s[i] - cum_weighted[i] / (cum_exp[i] + 1e-10)
        grad = -grad / max(np.sum(events_s), 1)
        return {'beta': grad}

    def step(self, grads):
        self.beta -= self.lr * grads['beta']

    def concordance_index(self, X, times, events):
        risk = X @ self.beta
        concordant = 0
        discordant = 0
        for i in range(len(times)):
            if events[i] == 0:
                continue
            for j in range(len(times)):
                if times[j] > times[i]:
                    if risk[i] > risk[j]:
                        concordant += 1
                    elif risk[i] < risk[j]:
                        discordant += 1
        total = concordant + discordant
        return concordant / total if total > 0 else 0.5


# ============================================================
# Data Generation
# ============================================================
def generate_medical_data(n_samples=2000, n_features=20, n_clients=5, iid=True):
    """Generate synthetic medical classification data for FL."""
    X, y = make_classification(
        n_samples=n_samples, n_features=n_features,
        n_informative=15, n_redundant=3, n_clusters_per_class=3,
        flip_y=0.05, random_state=42
    )
    if iid:
        indices = np.random.permutation(n_samples)
        split = np.array_split(indices, n_clients)
    else:
        # Non-IID: sort by label, then distribute unevenly
        sorted_idx = np.argsort(y)
        # Create shards with label skew
        shard_size = n_samples // (n_clients * 2)
        shards = [sorted_idx[i*shard_size:(i+1)*shard_size] for i in range(n_clients * 2)]
        np.random.shuffle(shards)
        split = [np.concatenate(shards[i*2:(i+1)*2]) for i in range(n_clients)]

    clients = []
    for idx in split:
        X_train, X_test, y_train, y_test = train_test_split(
            X[idx], y[idx], test_size=0.2, random_state=42
        )
        clients.append({'X_train': X_train, 'y_train': y_train,
                        'X_test': X_test, 'y_test': y_test})
    return clients, n_features


def generate_survival_data(n_samples=300, n_features=10, n_clients=5):
    """Generate synthetic multi-site survival data."""
    clients = []
    for c in range(n_clients):
        np.random.seed(42 + c)
        X = np.random.randn(n_samples, n_features)
        # Site-specific effects
        true_beta = np.random.randn(n_features) * 0.5
        true_beta[:3] = [0.8, -0.6, 0.4]  # shared important features
        hazard = np.exp(X @ true_beta)
        times = np.random.exponential(1.0 / (hazard + 0.01))
        times = np.clip(times, 0.01, 10.0)
        # Random censoring
        censor_times = np.random.exponential(3.0, n_samples)
        events = (times <= censor_times).astype(float)
        times = np.minimum(times, censor_times)

        idx = np.arange(n_samples)
        train_idx, test_idx = train_test_split(idx, test_size=0.2, random_state=42)
        clients.append({
            'X_train': X[train_idx], 'times_train': times[train_idx],
            'events_train': events[train_idx],
            'X_test': X[test_idx], 'times_test': times[test_idx],
            'events_test': events[test_idx]
        })
    return clients, n_features


# ============================================================
# Federated Learning Algorithms
# ============================================================
def fedavg(clients, n_features, n_rounds=50, local_epochs=5, lr=0.05):
    """Standard Federated Averaging."""
    global_model = LogisticRegression(n_features, lr=lr)
    history = {'loss': [], 'accuracy': [], 'auc': []}

    for r in range(n_rounds):
        local_params = []
        weights = []
        for c in clients:
            model = LogisticRegression(n_features, lr=lr)
            model.set_params(global_model.get_params())
            for _ in range(local_epochs):
                grads = model.compute_gradients(c['X_train'], c['y_train'])
                model.step(grads)
            local_params.append(model.get_params())
            weights.append(len(c['X_train']))

        # Aggregate
        total_w = sum(weights)
        new_params = {
            'w': sum(p['w'] * w / total_w for p, w in zip(local_params, weights)),
            'b': sum(p['b'] * w / total_w for p, w in zip(local_params, weights))
        }
        global_model.set_params(new_params)

        # Evaluate
        all_X = np.vstack([c['X_test'] for c in clients])
        all_y = np.concatenate([c['y_test'] for c in clients])
        loss = global_model.compute_loss(all_X, all_y)
        acc = accuracy_score(all_y, global_model.predict(all_X))
        auc = roc_auc_score(all_y, global_model.predict_proba(all_X))
        history['loss'].append(loss)
        history['accuracy'].append(acc)
        history['auc'].append(auc)

    return global_model, history


def fedprox(clients, n_features, n_rounds=50, local_epochs=5, lr=0.05, mu=0.1):
    """FedProx with proximal term."""
    global_model = LogisticRegression(n_features, lr=lr)
    history = {'loss': [], 'accuracy': [], 'auc': []}

    for r in range(n_rounds):
        local_params = []
        weights = []
        global_p = global_model.get_params()
        for c in clients:
            model = LogisticRegression(n_features, lr=lr)
            model.set_params(global_p)
            for _ in range(local_epochs):
                grads = model.compute_gradients(c['X_train'], c['y_train'])
                # Proximal term
                grads['w'] += mu * (model.w - global_p['w'])
                grads['b'] += mu * (model.b - global_p['b'])
                model.step(grads)
            local_params.append(model.get_params())
            weights.append(len(c['X_train']))

        total_w = sum(weights)
        new_params = {
            'w': sum(p['w'] * w / total_w for p, w in zip(local_params, weights)),
            'b': sum(p['b'] * w / total_w for p, w in zip(local_params, weights))
        }
        global_model.set_params(new_params)

        all_X = np.vstack([c['X_test'] for c in clients])
        all_y = np.concatenate([c['y_test'] for c in clients])
        loss = global_model.compute_loss(all_X, all_y)
        acc = accuracy_score(all_y, global_model.predict(all_X))
        auc = roc_auc_score(all_y, global_model.predict_proba(all_X))
        history['loss'].append(loss)
        history['accuracy'].append(acc)
        history['auc'].append(auc)

    return global_model, history


def scaffold(clients, n_features, n_rounds=50, local_epochs=5, lr=0.05):
    """SCAFFOLD with control variates."""
    global_model = LogisticRegression(n_features, lr=lr)
    c_global = {'w': np.zeros(n_features), 'b': 0.0}
    c_locals = [{'w': np.zeros(n_features), 'b': 0.0} for _ in clients]
    history = {'loss': [], 'accuracy': [], 'auc': []}

    for r in range(n_rounds):
        local_params = []
        weights = []
        delta_c_list = []
        global_p = global_model.get_params()

        for i, c in enumerate(clients):
            model = LogisticRegression(n_features, lr=lr)
            model.set_params(global_p)
            for _ in range(local_epochs):
                grads = model.compute_gradients(c['X_train'], c['y_train'])
                # Variance reduction via control variates
                grads['w'] += c_global['w'] - c_locals[i]['w']
                grads['b'] += c_global['b'] - c_locals[i]['b']
                model.step(grads)

            local_params.append(model.get_params())
            weights.append(len(c['X_train']))

            # Update control variate
            new_ci_w = c_locals[i]['w'] - c_global['w'] + (global_p['w'] - model.get_params()['w']) / (local_epochs * lr)
            new_ci_b = c_locals[i]['b'] - c_global['b'] + (global_p['b'] - model.get_params()['b']) / (local_epochs * lr)
            delta_c_list.append({
                'w': new_ci_w - c_locals[i]['w'],
                'b': new_ci_b - c_locals[i]['b']
            })
            c_locals[i] = {'w': new_ci_w, 'b': new_ci_b}

        total_w = sum(weights)
        new_params = {
            'w': sum(p['w'] * w / total_w for p, w in zip(local_params, weights)),
            'b': sum(p['b'] * w / total_w for p, w in zip(local_params, weights))
        }
        global_model.set_params(new_params)

        # Update global control variate
        n_c = len(clients)
        c_global['w'] += sum(d['w'] for d in delta_c_list) / n_c
        c_global['b'] += sum(d['b'] for d in delta_c_list) / n_c

        all_X = np.vstack([c['X_test'] for c in clients])
        all_y = np.concatenate([c['y_test'] for c in clients])
        loss = global_model.compute_loss(all_X, all_y)
        acc = accuracy_score(all_y, global_model.predict(all_X))
        auc = roc_auc_score(all_y, global_model.predict_proba(all_X))
        history['loss'].append(loss)
        history['accuracy'].append(acc)
        history['auc'].append(auc)

    return global_model, history


def fedavg_dp(clients, n_features, n_rounds=50, local_epochs=5, lr=0.05,
              noise_multiplier=0.5, max_grad_norm=1.0, delta=1e-5):
    """FedAvg with Differential Privacy (DP-FedAvg)."""
    global_model = LogisticRegression(n_features, lr=lr)
    history = {'loss': [], 'accuracy': [], 'auc': [], 'epsilon': []}
    cumulative_epsilon = 0.0

    for r in range(n_rounds):
        local_params = []
        weights = []
        for c in clients:
            model = LogisticRegression(n_features, lr=lr)
            model.set_params(global_model.get_params())
            for _ in range(local_epochs):
                grads = model.compute_gradients(c['X_train'], c['y_train'])
                # Gradient clipping
                grad_norm = np.sqrt(np.sum(grads['w']**2) + grads['b']**2)
                clip_factor = min(1.0, max_grad_norm / (grad_norm + 1e-10))
                grads['w'] *= clip_factor
                grads['b'] *= clip_factor
                # Add Gaussian noise
                grads['w'] += np.random.normal(0, noise_multiplier * max_grad_norm / len(c['X_train']),
                                                size=grads['w'].shape)
                grads['b'] += np.random.normal(0, noise_multiplier * max_grad_norm / len(c['X_train']))
                model.step(grads)
            local_params.append(model.get_params())
            weights.append(len(c['X_train']))

        total_w = sum(weights)
        new_params = {
            'w': sum(p['w'] * w / total_w for p, w in zip(local_params, weights)),
            'b': sum(p['b'] * w / total_w for p, w in zip(local_params, weights))
        }
        global_model.set_params(new_params)

        # Privacy accounting (simplified Gaussian mechanism)
        sigma = noise_multiplier
        per_round_eps = np.sqrt(2 * np.log(1.25 / delta)) / sigma
        cumulative_epsilon += per_round_eps

        all_X = np.vstack([c['X_test'] for c in clients])
        all_y = np.concatenate([c['y_test'] for c in clients])
        loss = global_model.compute_loss(all_X, all_y)
        acc = accuracy_score(all_y, global_model.predict(all_X))
        auc = roc_auc_score(all_y, global_model.predict_proba(all_X))
        history['loss'].append(loss)
        history['accuracy'].append(acc)
        history['auc'].append(auc)
        history['epsilon'].append(cumulative_epsilon)

    return global_model, history


def fedavg_compressed(clients, n_features, n_rounds=50, local_epochs=5, lr=0.05,
                       compression_ratio=0.3):
    """FedAvg with gradient compression (Top-K sparsification)."""
    global_model = LogisticRegression(n_features, lr=lr)
    history = {'loss': [], 'accuracy': [], 'auc': [], 'comm_bytes': []}
    residuals = [{'w': np.zeros(n_features), 'b': 0.0} for _ in clients]
    total_bytes = 0

    for r in range(n_rounds):
        local_params = []
        weights = []
        global_p = global_model.get_params()

        for i, c in enumerate(clients):
            model = LogisticRegression(n_features, lr=lr)
            model.set_params(global_p)
            for _ in range(local_epochs):
                grads = model.compute_gradients(c['X_train'], c['y_train'])
                model.step(grads)

            # Compute update
            update_w = model.w - global_p['w'] + residuals[i]['w']
            update_b = model.b - global_p['b'] + residuals[i]['b']

            # Top-K sparsification
            k = max(1, int(n_features * compression_ratio))
            top_k_idx = np.argsort(np.abs(update_w))[-k:]
            sparse_w = np.zeros(n_features)
            sparse_w[top_k_idx] = update_w[top_k_idx]

            # Error feedback
            residuals[i]['w'] = update_w - sparse_w
            residuals[i]['b'] = 0.0

            compressed_params = {
                'w': global_p['w'] + sparse_w,
                'b': global_p['b'] + update_b
            }
            local_params.append(compressed_params)
            weights.append(len(c['X_train']))
            total_bytes += k * 8 + 8  # k float64 values + bias

        total_w = sum(weights)
        new_params = {
            'w': sum(p['w'] * w / total_w for p, w in zip(local_params, weights)),
            'b': sum(p['b'] * w / total_w for p, w in zip(local_params, weights))
        }
        global_model.set_params(new_params)

        all_X = np.vstack([c['X_test'] for c in clients])
        all_y = np.concatenate([c['y_test'] for c in clients])
        loss = global_model.compute_loss(all_X, all_y)
        acc = accuracy_score(all_y, global_model.predict(all_X))
        auc = roc_auc_score(all_y, global_model.predict_proba(all_X))
        history['loss'].append(loss)
        history['accuracy'].append(acc)
        history['auc'].append(auc)
        history['comm_bytes'].append(total_bytes)

    return global_model, history


def fedavg_byzantine(clients, n_features, n_rounds=50, local_epochs=5, lr=0.05,
                      n_byzantine=1, defense='krum'):
    """FedAvg with Byzantine attack and defense."""
    global_model = LogisticRegression(n_features, lr=lr)
    history = {'loss': [], 'accuracy': [], 'auc': []}

    for r in range(n_rounds):
        local_params = []
        weights = []
        global_p = global_model.get_params()

        for i, c in enumerate(clients):
            model = LogisticRegression(n_features, lr=lr)
            model.set_params(global_p)
            for _ in range(local_epochs):
                grads = model.compute_gradients(c['X_train'], c['y_train'])
                model.step(grads)
            params = model.get_params()

            # Byzantine attack: last n_byzantine clients send adversarial updates
            if i >= len(clients) - n_byzantine:
                params['w'] = -params['w'] + np.random.randn(n_features) * 5
                params['b'] = -params['b'] + np.random.randn() * 5

            local_params.append(params)
            weights.append(len(c['X_train']))

        # Defense mechanism
        if defense == 'krum':
            # Multi-Krum: select the update closest to others
            n = len(local_params)
            f = n_byzantine
            distances = np.zeros((n, n))
            for i in range(n):
                for j in range(n):
                    diff_w = local_params[i]['w'] - local_params[j]['w']
                    diff_b = local_params[i]['b'] - local_params[j]['b']
                    distances[i, j] = np.sum(diff_w**2) + diff_b**2

            scores = []
            for i in range(n):
                sorted_dists = np.sort(distances[i])
                scores.append(np.sum(sorted_dists[1:n - f]))

            # Select top-(n-f) closest updates
            selected = np.argsort(scores)[:n - f]
            sel_weights = [weights[i] for i in selected]
            total_w = sum(sel_weights)
            new_params = {
                'w': sum(local_params[i]['w'] * sel_weights[j] / total_w
                         for j, i in enumerate(selected)),
                'b': sum(local_params[i]['b'] * sel_weights[j] / total_w
                         for j, i in enumerate(selected))
            }

        elif defense == 'trimmed_mean':
            all_w = np.stack([p['w'] for p in local_params])
            all_b = np.array([p['b'] for p in local_params])
            trim = n_byzantine
            sorted_w = np.sort(all_w, axis=0)
            sorted_b = np.sort(all_b)
            new_params = {
                'w': np.mean(sorted_w[trim:-trim], axis=0) if trim > 0 else np.mean(sorted_w, axis=0),
                'b': np.mean(sorted_b[trim:-trim]) if trim > 0 else np.mean(sorted_b)
            }

        elif defense == 'median':
            all_w = np.stack([p['w'] for p in local_params])
            all_b = np.array([p['b'] for p in local_params])
            new_params = {
                'w': np.median(all_w, axis=0),
                'b': np.median(all_b)
            }
        else:
            # No defense (naive averaging)
            total_w = sum(weights)
            new_params = {
                'w': sum(p['w'] * w / total_w for p, w in zip(local_params, weights)),
                'b': sum(p['b'] * w / total_w for p, w in zip(local_params, weights))
            }

        global_model.set_params(new_params)

        all_X = np.vstack([c['X_test'] for c in clients])
        all_y = np.concatenate([c['y_test'] for c in clients])
        loss = global_model.compute_loss(all_X, all_y)
        acc = accuracy_score(all_y, global_model.predict(all_X))
        auc = roc_auc_score(all_y, global_model.predict_proba(all_X))
        history['loss'].append(loss)
        history['accuracy'].append(acc)
        history['auc'].append(auc)

    return global_model, history


def federated_survival(clients, n_features, n_rounds=50, local_epochs=5, lr=0.001,
                        algorithm='fedavg', mu=0.1):
    """Federated Cox PH survival analysis."""
    global_model = CoxPH(n_features, lr=lr)
    history = {'nll': [], 'c_index': []}

    for r in range(n_rounds):
        local_params = []
        weights = []
        global_p = global_model.get_params()

        for c in clients:
            model = CoxPH(n_features, lr=lr)
            model.set_params(global_p)
            for _ in range(local_epochs):
                grads = model.compute_gradients(c['X_train'], c['times_train'], c['events_train'])
                if algorithm == 'fedprox':
                    grads['beta'] += mu * (model.beta - global_p['beta'])
                model.step(grads)
            local_params.append(model.get_params())
            weights.append(len(c['X_train']))

        total_w = sum(weights)
        new_params = {
            'beta': sum(p['beta'] * w / total_w for p, w in zip(local_params, weights))
        }
        global_model.set_params(new_params)

        # Evaluate
        all_X = np.vstack([c['X_test'] for c in clients])
        all_times = np.concatenate([c['times_test'] for c in clients])
        all_events = np.concatenate([c['events_test'] for c in clients])
        nll = global_model.negative_partial_log_likelihood(all_X, all_times, all_events)
        ci = global_model.concordance_index(all_X, all_times, all_events)
        history['nll'].append(nll)
        history['c_index'].append(ci)

    return global_model, history


# ============================================================
# Experiments
# ============================================================
def run_all_experiments():
    results = {}
    print("=" * 60)
    print("FEDERATED LEARNING EXPERIMENTS")
    print("=" * 60)

    # --- Experiment 1: FedAvg Convergence (IID vs Non-IID) ---
    print("\n[Exp 1] FedAvg Convergence: IID vs Non-IID")
    clients_iid, nf = generate_medical_data(n_samples=2000, n_features=20, n_clients=5, iid=True)
    clients_noniid, _ = generate_medical_data(n_samples=2000, n_features=20, n_clients=5, iid=False)

    _, hist_iid = fedavg(clients_iid, nf, n_rounds=60, local_epochs=5, lr=0.05)
    _, hist_noniid = fedavg(clients_noniid, nf, n_rounds=60, local_epochs=5, lr=0.05)
    print(f"  IID final: acc={hist_iid['accuracy'][-1]:.4f}, AUC={hist_iid['auc'][-1]:.4f}")
    print(f"  Non-IID final: acc={hist_noniid['accuracy'][-1]:.4f}, AUC={hist_noniid['auc'][-1]:.4f}")
    results['exp1'] = {'iid': hist_iid, 'noniid': hist_noniid}

    # --- Experiment 2: Non-IID Methods Comparison ---
    print("\n[Exp 2] Non-IID Methods: FedAvg vs FedProx vs SCAFFOLD")
    _, hist_fedavg = fedavg(clients_noniid, nf, n_rounds=60, local_epochs=5, lr=0.05)
    _, hist_fedprox = fedprox(clients_noniid, nf, n_rounds=60, local_epochs=5, lr=0.05, mu=0.1)
    _, hist_scaffold = scaffold(clients_noniid, nf, n_rounds=60, local_epochs=5, lr=0.05)
    print(f"  FedAvg:    acc={hist_fedavg['accuracy'][-1]:.4f}, AUC={hist_fedavg['auc'][-1]:.4f}")
    print(f"  FedProx:   acc={hist_fedprox['accuracy'][-1]:.4f}, AUC={hist_fedprox['auc'][-1]:.4f}")
    print(f"  SCAFFOLD:  acc={hist_scaffold['accuracy'][-1]:.4f}, AUC={hist_scaffold['auc'][-1]:.4f}")
    results['exp2'] = {'fedavg': hist_fedavg, 'fedprox': hist_fedprox, 'scaffold': hist_scaffold}

    # --- Experiment 3: Differential Privacy ---
    print("\n[Exp 3] Differential Privacy Impact")
    noise_levels = [0.0, 0.3, 0.5, 1.0, 2.0]
    dp_results = {}
    for sigma in noise_levels:
        if sigma == 0.0:
            _, h = fedavg(clients_iid, nf, n_rounds=60, local_epochs=5, lr=0.05)
            h['epsilon'] = [0.0] * 60
        else:
            _, h = fedavg_dp(clients_iid, nf, n_rounds=60, local_epochs=5, lr=0.05,
                              noise_multiplier=sigma)
        dp_results[sigma] = h
        print(f"  σ={sigma}: acc={h['accuracy'][-1]:.4f}, AUC={h['auc'][-1]:.4f}, ε={h['epsilon'][-1]:.2f}")
    results['exp3'] = dp_results

    # --- Experiment 4: Communication Efficiency ---
    print("\n[Exp 4] Communication Efficiency (Gradient Compression)")
    comp_ratios = [1.0, 0.5, 0.3, 0.1]
    comp_results = {}
    for ratio in comp_ratios:
        if ratio == 1.0:
            _, h = fedavg(clients_iid, nf, n_rounds=60, local_epochs=5, lr=0.05)
            h['comm_bytes'] = list(range(1, 61))
            h['comm_bytes'] = [b * nf * 8 * len(clients_iid) for b in h['comm_bytes']]
        else:
            _, h = fedavg_compressed(clients_iid, nf, n_rounds=60, local_epochs=5, lr=0.05,
                                      compression_ratio=ratio)
        comp_results[ratio] = h
        print(f"  ratio={ratio}: acc={h['accuracy'][-1]:.4f}, total_bytes={h['comm_bytes'][-1]}")
    results['exp4'] = comp_results

    # --- Experiment 5: Byzantine Robustness ---
    print("\n[Exp 5] Byzantine Robustness")
    clients_byz, nf_byz = generate_medical_data(n_samples=2000, n_features=20, n_clients=7, iid=True)
    defenses = ['none', 'krum', 'trimmed_mean', 'median']
    byz_results = {}
    for defense in defenses:
        _, h = fedavg_byzantine(clients_byz, nf_byz, n_rounds=60, local_epochs=5, lr=0.05,
                                 n_byzantine=2, defense=defense)
        byz_results[defense] = h
        print(f"  {defense}: acc={h['accuracy'][-1]:.4f}, AUC={h['auc'][-1]:.4f}")
    results['exp5'] = byz_results

    # --- Experiment 6: Survival Analysis ---
    print("\n[Exp 6] Multi-Site Survival Analysis")
    surv_clients, surv_nf = generate_survival_data(n_samples=300, n_features=10, n_clients=5)
    _, hist_surv_avg = federated_survival(surv_clients, surv_nf, n_rounds=80, local_epochs=5,
                                           lr=0.001, algorithm='fedavg')
    _, hist_surv_prox = federated_survival(surv_clients, surv_nf, n_rounds=80, local_epochs=5,
                                            lr=0.001, algorithm='fedprox', mu=0.01)
    print(f"  FedAvg:  C-index={hist_surv_avg['c_index'][-1]:.4f}")
    print(f"  FedProx: C-index={hist_surv_prox['c_index'][-1]:.4f}")
    results['exp6'] = {'fedavg': hist_surv_avg, 'fedprox': hist_surv_prox}

    return results


# ============================================================
# Plotting
# ============================================================
def generate_figures(results):
    plt.rcParams.update({'font.size': 11, 'figure.dpi': 150})

    # Figure 1: FedAvg Convergence IID vs Non-IID
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, metric, label in zip(axes, ['loss', 'accuracy', 'auc'],
                                  ['Loss', 'Accuracy', 'AUC-ROC']):
        ax.plot(results['exp1']['iid'][metric], label='IID', linewidth=2)
        ax.plot(results['exp1']['noniid'][metric], label='Non-IID', linewidth=2, linestyle='--')
        ax.set_xlabel('Communication Round')
        ax.set_ylabel(label)
        ax.set_title(f'FedAvg {label}: IID vs Non-IID')
        ax.legend()
        ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig1_fedavg_convergence.png', bbox_inches='tight')
    plt.close()
    print("  Saved fig1_fedavg_convergence.png")

    # Figure 2: Non-IID Methods Comparison
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    methods = ['fedavg', 'fedprox', 'scaffold']
    labels = ['FedAvg', 'FedProx', 'SCAFFOLD']
    colors = ['#2196F3', '#FF5722', '#4CAF50']
    for ax, metric, ylabel in zip(axes, ['loss', 'accuracy', 'auc'],
                                   ['Loss', 'Accuracy', 'AUC-ROC']):
        for m, l, c in zip(methods, labels, colors):
            ax.plot(results['exp2'][m][metric], label=l, linewidth=2, color=c)
        ax.set_xlabel('Communication Round')
        ax.set_ylabel(ylabel)
        ax.set_title(f'Non-IID {ylabel} Comparison')
        ax.legend()
        ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig2_noniid_methods.png', bbox_inches='tight')
    plt.close()
    print("  Saved fig2_noniid_methods.png")

    # Figure 3: Differential Privacy Trade-off
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    dp = results['exp3']
    for sigma, h in dp.items():
        label = f'σ={sigma}' if sigma > 0 else 'No DP'
        axes[0].plot(h['accuracy'], label=label, linewidth=1.5)
        axes[1].plot(h['auc'], label=label, linewidth=1.5)
    axes[0].set_xlabel('Communication Round'); axes[0].set_ylabel('Accuracy')
    axes[0].set_title('Accuracy under DP'); axes[0].legend(fontsize=9)
    axes[0].grid(True, alpha=0.3)
    axes[1].set_xlabel('Communication Round'); axes[1].set_ylabel('AUC-ROC')
    axes[1].set_title('AUC under DP'); axes[1].legend(fontsize=9)
    axes[1].grid(True, alpha=0.3)

    # Privacy-utility trade-off bar chart
    sigmas = [s for s in dp.keys() if s > 0]
    final_accs = [dp[s]['accuracy'][-1] for s in sigmas]
    final_eps = [dp[s]['epsilon'][-1] for s in sigmas]
    x = np.arange(len(sigmas))
    ax2_twin = axes[2].twinx()
    bars1 = axes[2].bar(x - 0.2, final_accs, 0.35, label='Accuracy', color='#2196F3', alpha=0.8)
    bars2 = ax2_twin.bar(x + 0.2, final_eps, 0.35, label='ε (cumulative)', color='#FF5722', alpha=0.8)
    axes[2].set_xticks(x)
    axes[2].set_xticklabels([f'σ={s}' for s in sigmas])
    axes[2].set_ylabel('Final Accuracy')
    ax2_twin.set_ylabel('Cumulative ε')
    axes[2].set_title('Privacy-Utility Trade-off')
    axes[2].legend(loc='upper left', fontsize=9)
    ax2_twin.legend(loc='upper right', fontsize=9)
    axes[2].grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig3_differential_privacy.png', bbox_inches='tight')
    plt.close()
    print("  Saved fig3_differential_privacy.png")

    # Figure 4: Communication Efficiency
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    comp = results['exp4']
    for ratio, h in comp.items():
        label = f'Top-{int(ratio*100)}%' if ratio < 1.0 else 'Full'
        axes[0].plot(h['accuracy'], label=label, linewidth=1.5)
    axes[0].set_xlabel('Communication Round'); axes[0].set_ylabel('Accuracy')
    axes[0].set_title('Accuracy vs Compression Ratio'); axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    ratios = list(comp.keys())
    final_accs = [comp[r]['accuracy'][-1] for r in ratios]
    final_bytes = [comp[r]['comm_bytes'][-1] / 1e6 for r in ratios]
    axes[1].bar([f'{int(r*100)}%' for r in ratios], final_bytes, color=['#4CAF50','#2196F3','#FF9800','#F44336'])
    axes[1].set_xlabel('Compression Ratio (% params sent)')
    axes[1].set_ylabel('Total Communication (MB)')
    axes[1].set_title('Communication Cost Reduction')
    for i, (acc, b) in enumerate(zip(final_accs, final_bytes)):
        axes[1].text(i, b + 0.001, f'acc={acc:.3f}', ha='center', fontsize=9)
    axes[1].grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig4_communication_efficiency.png', bbox_inches='tight')
    plt.close()
    print("  Saved fig4_communication_efficiency.png")

    # Figure 5: Byzantine Robustness
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    byz = results['exp5']
    defense_labels = {'none': 'No Defense', 'krum': 'Multi-Krum',
                      'trimmed_mean': 'Trimmed Mean', 'median': 'Coordinate Median'}
    defense_colors = {'none': '#F44336', 'krum': '#2196F3',
                      'trimmed_mean': '#4CAF50', 'median': '#FF9800'}
    for d, h in byz.items():
        axes[0].plot(h['accuracy'], label=defense_labels[d], linewidth=2,
                     color=defense_colors[d])
        axes[1].plot(h['auc'], label=defense_labels[d], linewidth=2,
                     color=defense_colors[d])
    axes[0].set_xlabel('Communication Round'); axes[0].set_ylabel('Accuracy')
    axes[0].set_title('Byzantine Robustness: Accuracy'); axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[1].set_xlabel('Communication Round'); axes[1].set_ylabel('AUC-ROC')
    axes[1].set_title('Byzantine Robustness: AUC'); axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig5_byzantine_robustness.png', bbox_inches='tight')
    plt.close()
    print("  Saved fig5_byzantine_robustness.png")

    # Figure 6: Survival Analysis
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    surv = results['exp6']
    axes[0].plot(surv['fedavg']['c_index'], label='FedAvg', linewidth=2, color='#2196F3')
    axes[0].plot(surv['fedprox']['c_index'], label='FedProx', linewidth=2, color='#FF5722')
    axes[0].set_xlabel('Communication Round'); axes[0].set_ylabel('C-index')
    axes[0].set_title('Survival Analysis: Concordance Index'); axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(surv['fedavg']['nll'], label='FedAvg', linewidth=2, color='#2196F3')
    axes[1].plot(surv['fedprox']['nll'], label='FedProx', linewidth=2, color='#FF5722')
    axes[1].set_xlabel('Communication Round'); axes[1].set_ylabel('Neg. Partial Log-Likelihood')
    axes[1].set_title('Survival Analysis: Loss Convergence'); axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig6_survival_analysis.png', bbox_inches='tight')
    plt.close()
    print("  Saved fig6_survival_analysis.png")

    # Figure 7: Summary comparison table as heatmap
    fig, ax = plt.subplots(figsize=(10, 5))
    methods_all = ['FedAvg\n(IID)', 'FedAvg\n(Non-IID)', 'FedProx\n(Non-IID)',
                   'SCAFFOLD\n(Non-IID)', 'DP-FedAvg\n(σ=0.5)', 'Compressed\n(30%)',
                   'Byzantine\n(Krum)']
    metrics = ['Accuracy', 'AUC-ROC']
    data = np.array([
        [results['exp1']['iid']['accuracy'][-1], results['exp1']['iid']['auc'][-1]],
        [results['exp1']['noniid']['accuracy'][-1], results['exp1']['noniid']['auc'][-1]],
        [results['exp2']['fedprox']['accuracy'][-1], results['exp2']['fedprox']['auc'][-1]],
        [results['exp2']['scaffold']['accuracy'][-1], results['exp2']['scaffold']['auc'][-1]],
        [results['exp3'][0.5]['accuracy'][-1], results['exp3'][0.5]['auc'][-1]],
        [results['exp4'][0.3]['accuracy'][-1], results['exp4'][0.3]['auc'][-1]],
        [results['exp5']['krum']['accuracy'][-1], results['exp5']['krum']['auc'][-1]],
    ])
    im = ax.imshow(data, cmap='RdYlGn', aspect='auto', vmin=0.5, vmax=1.0)
    ax.set_xticks(np.arange(len(metrics)))
    ax.set_yticks(np.arange(len(methods_all)))
    ax.set_xticklabels(metrics)
    ax.set_yticklabels(methods_all)
    for i in range(len(methods_all)):
        for j in range(len(metrics)):
            ax.text(j, i, f'{data[i, j]:.3f}', ha='center', va='center', fontsize=11,
                    color='black', fontweight='bold')
    ax.set_title('Summary: Final Performance Across All Methods', fontsize=13)
    plt.colorbar(im, ax=ax, label='Score')
    plt.tight_layout()
    plt.savefig(f'{FIGURES_DIR}/fig7_summary_heatmap.png', bbox_inches='tight')
    plt.close()
    print("  Saved fig7_summary_heatmap.png")


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    results = run_all_experiments()
    print("\n" + "=" * 60)
    print("GENERATING FIGURES")
    print("=" * 60)
    generate_figures(results)

    # Save numerical results
    summary = {
        'exp1_iid_final_acc': results['exp1']['iid']['accuracy'][-1],
        'exp1_iid_final_auc': results['exp1']['iid']['auc'][-1],
        'exp1_noniid_final_acc': results['exp1']['noniid']['accuracy'][-1],
        'exp1_noniid_final_auc': results['exp1']['noniid']['auc'][-1],
        'exp2_fedavg_acc': results['exp2']['fedavg']['accuracy'][-1],
        'exp2_fedprox_acc': results['exp2']['fedprox']['accuracy'][-1],
        'exp2_scaffold_acc': results['exp2']['scaffold']['accuracy'][-1],
        'exp2_fedavg_auc': results['exp2']['fedavg']['auc'][-1],
        'exp2_fedprox_auc': results['exp2']['fedprox']['auc'][-1],
        'exp2_scaffold_auc': results['exp2']['scaffold']['auc'][-1],
        'exp3_dp_sigma0.5_acc': results['exp3'][0.5]['accuracy'][-1],
        'exp3_dp_sigma1.0_acc': results['exp3'][1.0]['accuracy'][-1],
        'exp3_dp_sigma0.5_eps': results['exp3'][0.5]['epsilon'][-1],
        'exp3_dp_sigma1.0_eps': results['exp3'][1.0]['epsilon'][-1],
        'exp4_full_bytes': results['exp4'][1.0]['comm_bytes'][-1],
        'exp4_30pct_bytes': results['exp4'][0.3]['comm_bytes'][-1],
        'exp4_30pct_acc': results['exp4'][0.3]['accuracy'][-1],
        'exp5_no_defense_acc': results['exp5']['none']['accuracy'][-1],
        'exp5_krum_acc': results['exp5']['krum']['accuracy'][-1],
        'exp5_trimmed_acc': results['exp5']['trimmed_mean']['accuracy'][-1],
        'exp5_median_acc': results['exp5']['median']['accuracy'][-1],
        'exp6_fedavg_cindex': results['exp6']['fedavg']['c_index'][-1],
        'exp6_fedprox_cindex': results['exp6']['fedprox']['c_index'][-1],
    }

    with open('results_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 60)
    print("ALL EXPERIMENTS COMPLETE")
    print("=" * 60)
    for k, v in summary.items():
        print(f"  {k}: {v:.4f}" if isinstance(v, float) else f"  {k}: {v}")
