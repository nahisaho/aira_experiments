"""
Multi-Objective Bayesian Optimization for HEA Composition Search
Uses Expected Hypervolume Improvement (EHVI) via botorch.
Objectives: maximize yield_strength, maximize elongation, maximize pitting_potential
Constraint: sum of compositions = 1, each fraction in [0.05, 0.50]
"""

import numpy as np
import pandas as pd
import torch
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings("ignore")

torch.manual_seed(42)
np.random.seed(42)


ELEMENTS_5 = ["Cr", "Mn", "Fe", "Co", "Ni"]   # Cantor alloy family
ELEMENTS_7 = ["Cr", "Mn", "Fe", "Co", "Ni", "Al", "Ti"]  # Extended


def _normalize_composition(x: np.ndarray, elements: List[str]) -> Dict[str, float]:
    """Normalize raw fractions to sum = 1."""
    x = np.clip(x, 0.01, 1.0)
    x = x / x.sum()
    return {el: float(v) for el, v in zip(elements, x)}


def sample_compositions(n: int, elements: List[str],
                        x_min: float = 0.05, x_max: float = 0.60,
                        seed: int = 42) -> List[Dict[str, float]]:
    """
    Dirichlet sampling of compositions with min/max fraction constraints.
    """
    rng = np.random.default_rng(seed)
    compositions = []
    n_el = len(elements)

    while len(compositions) < n:
        x = rng.dirichlet(np.ones(n_el))
        if x_min <= x.min() and x.max() <= x_max:
            compositions.append({el: float(v) for el, v in zip(elements, x)})

    return compositions


def latin_hypercube_compositions(n: int, elements: List[str],
                                  x_min: float = 0.05) -> List[Dict[str, float]]:
    """LHS sampling on the simplex via sequential rejection."""
    from scipy.stats.qmc import LatinHypercube
    n_el = len(elements)
    sampler = LatinHypercube(d=n_el, seed=42)
    raw = sampler.random(n * 5)  # oversample, filter
    comps = []
    for row in raw:
        row = np.clip(row, x_min, 1.0)
        row = row / row.sum()
        if row.min() >= x_min:
            comps.append({el: float(v) for el, v in zip(elements, row)})
        if len(comps) >= n:
            break
    # pad if needed
    while len(comps) < n:
        comps += sample_compositions(n - len(comps), elements, x_min)
    return comps[:n]


# -----------------------------------------------------------------------
# Surrogate-based multi-objective optimizer (Pareto front discovery)
# Uses GP surrogate predictions on large candidate pool — much faster
# than botorch EHVI for demonstration purposes.
# -----------------------------------------------------------------------
class MultiObjectiveBayesianOptimizer:
    """
    Surrogate-assisted multi-objective Bayesian optimization.
    Strategy: GP surrogate prediction on LHS candidate pool,
    Pareto front extraction, and upper-confidence-bound (UCB) acquisition
    for sequential batch selection.
    """

    def __init__(self,
                 elements: List[str] = ELEMENTS_5,
                 ref_point: List[float] = None,
                 x_min: float = 0.05,
                 x_max: float = 0.60):
        self.elements = elements
        self.x_min = x_min
        self.x_max = x_max
        self.ref_point = ref_point or [200.0, 5.0, -0.8]
        self._train_comps: List[Dict] = []
        self._train_obj: Optional[np.ndarray] = None
        self._surrogate = None

    def initialize(self, compositions: List[Dict[str, float]],
                   objectives: np.ndarray):
        """Seed with initial data. objectives shape: [n, 3]."""
        self._train_comps = list(compositions)
        self._train_obj = np.array(objectives)

    def _fit_surrogate(self):
        """Fit GP surrogates on current training set."""
        from src.hea_descriptors import descriptors_dataframe
        from src.surrogate_models import HEASurrogateModel, DESCRIPTOR_COLS, PROPERTY_COLS
        import pandas as pd

        desc_df = descriptors_dataframe(self._train_comps, T=1000.0)
        y_df = pd.DataFrame(self._train_obj, columns=PROPERTY_COLS)

        surrogate = HEASurrogateModel(features=DESCRIPTOR_COLS)
        surrogate.fit(desc_df, y_df, n_restarts=3)
        self._surrogate = surrogate
        self._desc_cols = DESCRIPTOR_COLS
        return surrogate

    def suggest_next(self, batch_size: int = 4,
                     n_candidates: int = 500) -> List[Dict[str, float]]:
        """
        Generate candidate pool → predict with GP UCB → return top-k by acquisition.
        Acquisition = UCB(mu + 2*sigma) balanced across objectives.
        """
        from src.hea_descriptors import descriptors_dataframe
        import pandas as pd

        surrogate = self._fit_surrogate()

        # Sample large candidate pool
        cands = sample_compositions(n_candidates, self.elements,
                                    x_min=self.x_min, x_max=self.x_max,
                                    seed=np.random.randint(1, 10000))
        cand_desc = descriptors_dataframe(cands, T=1000.0)
        preds = surrogate.predict(cand_desc, return_std=True)

        from src.surrogate_models import PROPERTY_COLS
        # UCB score per objective (normalized)
        ucb = np.zeros(len(cands))
        for prop, (mu, std) in preds.items():
            mu_n = (mu - mu.min()) / (mu.max() - mu.min() + 1e-9)
            std_n = (std - std.min()) / (std.max() - std.min() + 1e-9)
            ucb += (mu_n + 2 * std_n)

        # Select diverse top candidates
        idx_sorted = np.argsort(ucb)[::-1]
        selected = []
        for idx in idx_sorted:
            if len(selected) >= batch_size:
                break
            selected.append(cands[idx])
        return selected

    def update(self, new_compositions: List[Dict[str, float]],
               new_objectives: np.ndarray):
        """Append new observations."""
        self._train_comps.extend(new_compositions)
        self._train_obj = np.vstack([self._train_obj, new_objectives])

    def pareto_front(self) -> Tuple[np.ndarray, List[Dict[str, float]]]:
        """Extract current Pareto-optimal points."""
        Y = self._train_obj
        pareto_mask = np.ones(len(Y), dtype=bool)
        for i, yi in enumerate(Y):
            for j, yj in enumerate(Y):
                if i != j and np.all(yj >= yi) and np.any(yj > yi):
                    pareto_mask[i] = False
                    break

        pareto_Y = Y[pareto_mask]
        pareto_comps = [self._train_comps[i] for i in np.where(pareto_mask)[0]]
        return pareto_Y, pareto_comps

    def hypervolume(self) -> float:
        """Approximate hypervolume indicator using dominated area calculation."""
        Y = self._train_obj
        ref = np.array(self.ref_point)
        dominated = Y - ref
        dominated = dominated[np.all(dominated > 0, axis=1)]
        if len(dominated) == 0:
            return 0.0
        # Approximate 3-objective HV via Monte Carlo
        rng = np.random.default_rng(42)
        ub = dominated.max(axis=0) + ref
        samples = rng.uniform(ref, ub, size=(50000, 3))
        dominated_count = np.sum(
            np.any(np.all(samples[:, None, :] <= Y[None, :, :], axis=2), axis=1)
        )
        vol = np.prod(ub - ref)
        return float(vol * dominated_count / 50000)


# -----------------------------------------------------------------------
# Active learning selector (uncertainty-based query strategy)
# -----------------------------------------------------------------------
class ActiveLearningSelector:
    """
    Query-by-committee / uncertainty sampling for efficient experimental design.
    Selects candidates that maximize:
    (a) Predictive uncertainty (pure exploration)
    (b) Expected improvement (balanced exploitation/exploration)
    (c) Diversity (max-min distance from existing observations)
    """

    def __init__(self, strategy: str = "uncertainty"):
        assert strategy in ("uncertainty", "diversity", "hybrid")
        self.strategy = strategy

    def select(self,
               candidate_descriptors: pd.DataFrame,
               surrogate: "HEASurrogateModel",
               n_select: int = 5,
               existing_descriptors: pd.DataFrame = None) -> pd.DataFrame:
        """
        Returns top-n_select candidates from the candidate pool.
        """
        predictions = surrogate.predict(candidate_descriptors)

        if self.strategy == "uncertainty":
            total_std = np.zeros(len(candidate_descriptors))
            for prop, (mu, std) in predictions.items():
                total_std += std / (std.max() + 1e-10)
            scores = total_std

        elif self.strategy == "diversity":
            assert existing_descriptors is not None
            feat_cols = surrogate.features
            X_cand = candidate_descriptors[feat_cols].values
            X_exist = existing_descriptors[feat_cols].values
            from scipy.spatial.distance import cdist
            D = cdist(X_cand, X_exist, metric="euclidean")
            scores = D.min(axis=1)

        else:  # hybrid
            total_std = np.zeros(len(candidate_descriptors))
            for prop, (mu, std) in predictions.items():
                total_std += std / (std.max() + 1e-10)

            feat_cols = surrogate.features
            X_cand = candidate_descriptors[feat_cols].values
            X_exist = existing_descriptors[feat_cols].values if existing_descriptors is not None else X_cand
            from scipy.spatial.distance import cdist
            D = cdist(X_cand, X_exist, metric="euclidean")
            diversity = D.min(axis=1) / (D.min(axis=1).max() + 1e-10)

            scores = 0.5 * (total_std / total_std.max()) + 0.5 * diversity

        idx = np.argsort(scores)[::-1][:n_select]
        return candidate_descriptors.iloc[idx].copy(), scores[idx]



ELEMENTS_5 = ["Cr", "Mn", "Fe", "Co", "Ni"]   # Cantor alloy family
ELEMENTS_7 = ["Cr", "Mn", "Fe", "Co", "Ni", "Al", "Ti"]  # Extended


def _normalize_composition(x: np.ndarray, elements: List[str]) -> Dict[str, float]:
    """Normalize raw fractions to sum = 1."""
    x = np.clip(x, 0.01, 1.0)
    x = x / x.sum()
    return {el: float(v) for el, v in zip(elements, x)}


def sample_compositions(n: int, elements: List[str],
                        x_min: float = 0.05, x_max: float = 0.60,
                        seed: int = 42) -> List[Dict[str, float]]:
    """
    Dirichlet sampling of compositions with min/max fraction constraints.
    """
    rng = np.random.default_rng(seed)
    compositions = []
    n_el = len(elements)

    while len(compositions) < n:
        x = rng.dirichlet(np.ones(n_el))
        if x_min <= x.min() and x.max() <= x_max:
            compositions.append({el: float(v) for el, v in zip(elements, x)})

    return compositions


def latin_hypercube_compositions(n: int, elements: List[str],
                                  x_min: float = 0.05) -> List[Dict[str, float]]:
    """LHS sampling on the simplex via sequential rejection."""
    from scipy.stats.qmc import LatinHypercube
    n_el = len(elements)
    sampler = LatinHypercube(d=n_el, seed=42)
    raw = sampler.random(n * 5)  # oversample, filter
    comps = []
    for row in raw:
        row = np.clip(row, x_min, 1.0)
        row = row / row.sum()
        if row.min() >= x_min:
            comps.append({el: float(v) for el, v in zip(elements, row)})
        if len(comps) >= n:
            break
    # pad if needed
    while len(comps) < n:
        comps += sample_compositions(n - len(comps), elements, x_min)
    return comps[:n]


