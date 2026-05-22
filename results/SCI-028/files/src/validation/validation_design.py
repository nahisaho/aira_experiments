"""
Validation design for JET/KSTAR experimental data.
Covers:
1. Train/validation/test split strategy for disruption databases
2. Performance metrics (AUC, TPR at alarm threshold, FPR, Heidke Skill Score)
3. Bootstrap confidence intervals
4. Cross-device validation (leave-one-device-out)
5. Temporal validation (prevent future leakage)
6. Alarm latency vs. TPR Pareto analysis
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import stats
from scipy.interpolate import interp1d

logger = logging.getLogger("tokamak.validation")


# ─── Disruption database schema ───────────────────────────────────────────────

@dataclass
class ShotRecord:
    """Metadata for a single plasma discharge."""
    shot_id: int
    device: str           # "JET", "KSTAR", "ASDEX_UG"
    disrupted: bool
    disruption_time_s: Optional[float]   # None if not disrupted
    duration_s: float
    ip_max_ma: float
    betan_max: float
    cause: Optional[str]  # "VDE", "locked_mode", "NTM", "radiation_collapse", "other"


@dataclass
class ValidationSplit:
    """Train / validation / test shot indices."""
    train_shots: List[int]
    val_shots: List[int]
    test_shots: List[int]
    split_method: str
    device_split: Dict[str, List[int]]  # per-device


# ─── Split strategies ─────────────────────────────────────────────────────────

def temporal_split(
    shots: List[ShotRecord],
    train_frac: float = 0.70,
    val_frac: float = 0.15,
    # test_frac implied: 1 - train - val
) -> ValidationSplit:
    """
    Temporal split: train on earlier shots, test on later shots.
    Prevents look-ahead bias in time-series validation.

    JET has ~90,000 shots; use shots up to ~2015 for train,
    2015–2019 for validation, 2019+ for test.
    """
    shots_sorted = sorted(shots, key=lambda s: s.shot_id)
    n = len(shots_sorted)
    n_train = int(n * train_frac)
    n_val   = int(n * val_frac)

    train = [s.shot_id for s in shots_sorted[:n_train]]
    val   = [s.shot_id for s in shots_sorted[n_train:n_train + n_val]]
    test  = [s.shot_id for s in shots_sorted[n_train + n_val:]]

    device_split: Dict[str, List[int]] = {}
    for s in shots_sorted:
        device_split.setdefault(s.device, []).append(s.shot_id)

    return ValidationSplit(train, val, test, "temporal", device_split)


def leave_one_device_out_split(
    shots: List[ShotRecord],
    held_out_device: str = "KSTAR",
) -> ValidationSplit:
    """
    Leave-one-device-out split for cross-device generalisation validation.
    Train on all other devices, test on held_out_device.
    """
    train = [s.shot_id for s in shots if s.device != held_out_device]
    test  = [s.shot_id for s in shots if s.device == held_out_device]

    device_split: Dict[str, List[int]] = {}
    for s in shots:
        device_split.setdefault(s.device, []).append(s.shot_id)

    # Val = 15% of train
    n_val = int(len(train) * 0.15)
    rng = np.random.default_rng(42)
    rng.shuffle(train)
    val = train[:n_val]
    train = train[n_val:]

    return ValidationSplit(train, val, test, f"lodo_{held_out_device}", device_split)


# ─── Disruption-specific metrics ──────────────────────────────────────────────

@dataclass
class DisruptionMetrics:
    """Standard metrics for disruption alarm performance evaluation."""
    # Binary classification
    auc_roc: float
    auc_pr: float             # Precision-Recall AUC (imbalanced class aware)
    # At operating point (p_threshold)
    threshold: float
    tpr: float                # True positive rate (sensitivity)
    fpr: float                # False positive rate
    precision: float
    f1: float
    # Disruption-specific
    average_warning_time_ms: float   # Mean TTD when true alarm raised
    missed_disruptions: int
    false_alarms_per_hour: float
    # Heidke Skill Score (accounts for random chance)
    hss: float
    # Time-to-disruption error
    ttd_mae_ms: float
    ttd_rmse_ms: float
    # Bootstrap 95% CI
    auc_ci: Tuple[float, float] = (0.0, 0.0)
    tpr_ci: Tuple[float, float] = (0.0, 0.0)


def compute_disruption_metrics(
    y_true: np.ndarray,          # (N,) binary labels
    y_prob: np.ndarray,          # (N,) disruption probabilities
    ttd_true: np.ndarray,        # (N,) true TTD in ms (-1 if no disruption)
    ttd_pred: np.ndarray,        # (N,) predicted TTD in ms
    threshold: float = 0.5,
    shot_duration_h: float = 10.0,
) -> DisruptionMetrics:
    """
    Compute comprehensive disruption prediction metrics.
    """
    y_pred = (y_prob >= threshold).astype(int)

    # AUC-ROC
    auc_roc = _auc_roc(y_true, y_prob)
    auc_pr  = _auc_pr(y_true, y_prob)

    # Confusion matrix
    tp = int(np.sum((y_pred == 1) & (y_true == 1)))
    fp = int(np.sum((y_pred == 1) & (y_true == 0)))
    tn = int(np.sum((y_pred == 0) & (y_true == 0)))
    fn = int(np.sum((y_pred == 0) & (y_true == 1)))

    tpr       = tp / (tp + fn + 1e-9)
    fpr       = fp / (fp + tn + 1e-9)
    precision = tp / (tp + fp + 1e-9)
    recall    = tpr
    f1        = 2 * precision * recall / (precision + recall + 1e-9)

    # Heidke Skill Score
    n = len(y_true)
    expected_correct = (
        ((tp + fn) * (tp + fp) + (tn + fp) * (tn + fn)) / (n**2 + 1e-9)
    )
    hss = (tp + tn - n * expected_correct) / (n - n * expected_correct + 1e-9)

    # TTD metrics (only on true positives)
    disrupted_mask = (y_true == 1) & (y_pred == 1)
    if disrupted_mask.sum() > 0:
        ttd_err = ttd_pred[disrupted_mask] - ttd_true[disrupted_mask]
        ttd_mae  = float(np.mean(np.abs(ttd_err)))
        ttd_rmse = float(np.sqrt(np.mean(ttd_err**2)))
        avg_warning = float(np.mean(ttd_true[disrupted_mask]))
    else:
        ttd_mae = ttd_rmse = avg_warning = 0.0

    # False alarms per hour
    fa_per_hour = fp / (shot_duration_h + 1e-9)

    # Bootstrap CI (1000 iterations)
    auc_ci, tpr_ci = _bootstrap_ci(y_true, y_prob, threshold, n_boot=500)

    return DisruptionMetrics(
        auc_roc=auc_roc, auc_pr=auc_pr,
        threshold=threshold,
        tpr=tpr, fpr=fpr, precision=precision, f1=f1,
        average_warning_time_ms=avg_warning,
        missed_disruptions=fn,
        false_alarms_per_hour=fa_per_hour,
        hss=hss,
        ttd_mae_ms=ttd_mae, ttd_rmse_ms=ttd_rmse,
        auc_ci=auc_ci, tpr_ci=tpr_ci,
    )


def _auc_roc(y: np.ndarray, p: np.ndarray) -> float:
    """Mann-Whitney AUC-ROC."""
    pos = p[y == 1]
    neg = p[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return 0.5
    return float(np.mean(pos[:, None] > neg[None, :]))


def _auc_pr(y: np.ndarray, p: np.ndarray) -> float:
    """Precision-Recall AUC via trapezoidal integration."""
    thresholds = np.linspace(0, 1, 101)
    precisions, recalls = [], []
    for t in thresholds:
        yp = (p >= t).astype(int)
        tp = int(np.sum((yp == 1) & (y == 1)))
        fp = int(np.sum((yp == 1) & (y == 0)))
        fn = int(np.sum((yp == 0) & (y == 1)))
        precisions.append(tp / (tp + fp + 1e-9))
        recalls.append(tp / (tp + fn + 1e-9))
    return float(np.trapz(precisions[::-1], recalls[::-1]))


def _bootstrap_ci(
    y: np.ndarray,
    p: np.ndarray,
    threshold: float,
    n_boot: int = 500,
    alpha: float = 0.05,
) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    rng = np.random.default_rng(42)
    auc_boots, tpr_boots = [], []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yb, pb = y[idx], p[idx]
        auc_boots.append(_auc_roc(yb, pb))
        yp = (pb >= threshold).astype(int)
        tp = np.sum((yp == 1) & (yb == 1))
        fn = np.sum((yp == 0) & (yb == 1))
        tpr_boots.append(float(tp / (tp + fn + 1e-9)))
    al, ah = alpha / 2, 1 - alpha / 2
    return (
        (float(np.quantile(auc_boots, al)), float(np.quantile(auc_boots, ah))),
        (float(np.quantile(tpr_boots, al)), float(np.quantile(tpr_boots, ah))),
    )


# ─── Pareto analysis: warning time vs. TPR ───────────────────────────────────

def warning_time_tpr_pareto(
    shot_records: List[ShotRecord],
    predicted_probs: np.ndarray,   # (N_time,) probability time series
    timestamps: np.ndarray,        # (N_time,) time values [s]
    disruption_times: np.ndarray,  # (N_shots,) -1 if no disruption
    thresholds: Optional[np.ndarray] = None,
) -> Dict[str, np.ndarray]:
    """
    For each threshold, compute (mean_warning_time [ms], TPR).
    Returns the Pareto front for threshold selection.
    """
    if thresholds is None:
        thresholds = np.linspace(0.1, 0.99, 50)

    results = {"threshold": [], "tpr": [], "warning_time_ms": [], "fpr": []}

    for t in thresholds:
        alarms = predicted_probs >= t
        tpr_list, wt_list = [], []
        for i, rec in enumerate(shot_records):
            if not rec.disrupted:
                continue
            dt = disruption_times[i]
            # Find first alarm before disruption
            before_mask = (timestamps < dt) & alarms
            if before_mask.any():
                first_alarm_t = timestamps[before_mask][0]
                wt_list.append((dt - first_alarm_t) * 1000.0)  # → ms
                tpr_list.append(1.0)
            else:
                tpr_list.append(0.0)
        results["threshold"].append(float(t))
        results["tpr"].append(float(np.mean(tpr_list)) if tpr_list else 0.0)
        results["warning_time_ms"].append(float(np.mean(wt_list)) if wt_list else 0.0)

    return {k: np.array(v) for k, v in results.items()}


# ─── Cross-device validation summary ─────────────────────────────────────────

def cross_device_validation_summary(
    results_by_device: Dict[str, DisruptionMetrics],
) -> str:
    """Format a markdown table of cross-device validation results."""
    header = "| Device | AUC-ROC | TPR | FPR | HSS | Avg. Warning [ms] | FA/hr |\n"
    sep    = "|--------|---------|-----|-----|-----|-------------------|-------|\n"
    rows   = []
    for dev, m in results_by_device.items():
        rows.append(
            f"| {dev:<8} | {m.auc_roc:.3f}   | {m.tpr:.3f} | {m.fpr:.3f} | {m.hss:.3f} "
            f"| {m.average_warning_time_ms:>17.1f} | {m.false_alarms_per_hour:.2f}  |"
        )
    return header + sep + "\n".join(rows)


# ─── Validation protocol document ────────────────────────────────────────────

VALIDATION_PROTOCOL = """
# Validation Protocol: Tokamak Disruption Prediction AI

## 1. Dataset Requirements

### JET (primary training + test)
- Total shots: ~90,000 (1983–2022)
- Disruption rate: ~15–20% of high-performance shots
- Required signals: Mirnov array (32 coils), Thomson scattering, ECE, magnetics
- Recommended: JET Disruption Database (JDD) — contact EUROfusion
- Shot range for test: 90000–99000 (most recent campaigns)

### KSTAR (cross-device validation)
- Total shots: ~30,000 (2008–2023)
- Disruption rate: ~10% (higher H-mode fraction)
- Required signals: Same as JET (mapped via signal name registry)
- Data access: KSTAR data portal or IMAS MDSplus server

### ASDEX Upgrade (optional supplementary)
- Total shots: ~40,000
- Useful for: impurity-induced disruption scenarios

## 2. Preprocessing Checklist
- [ ] Synchronise signal timestamps to 0.1 ms grid
- [ ] Apply per-device normalisation (IP/IP_max, BT/BT_nom)
- [ ] Remove shots with incomplete data (<80% signal coverage)
- [ ] Label disruption time as last timestamp before dIp/dt < −0.5 MA/ms
- [ ] Assign disruption cause from Mirnov/radiation pattern (automated + manual audit)

## 3. Split Rules
- NEVER mix shot segments from the same discharge across train/test
- Use temporal split (shot ID ascending) to prevent future leakage
- Hold out entire JET campaigns for final test (e.g., C38, C40)
- KSTAR test set: all shots, zero-shot transfer evaluation

## 4. Evaluation Protocol
1. Train on JET train set (shots before 2018)
2. Validate hyperparameters on JET val set (2018–2020)
3. Report final metrics on JET test set (2020+)
4. Report zero-shot metrics on KSTAR (no fine-tuning)
5. Report few-shot metrics on KSTAR (10/50/100 shots fine-tuning)
6. Compare against:
   - APODIS baseline (Versace et al., 2010)
   - SVM baseline (Rattá et al., 2010)
   - LSTM baseline (Kates-Harbeck et al., 2019 — FRNN)

## 5. Minimum Acceptance Thresholds
| Metric              | Minimum | Target |
|---------------------|---------|--------|
| AUC-ROC (JET)       | 0.92    | 0.97   |
| TPR @ FPR=0.05      | 0.90    | 0.95   |
| Avg. Warning Time   | 30 ms   | 100 ms |
| KSTAR AUC (0-shot)  | 0.80    | 0.90   |
| Inference Latency   | <30 ms  | <20 ms |

## 6. Statistical Testing
- Compare models using DeLong test for AUC-ROC differences
- Bootstrap 95% CI (n=1000) for all point estimates
- McNemar's test for shot-level TPR/FPR comparison
- Bonferroni correction for multiple device comparisons (k=3)

## 7. Prospective Validation (future)
- Submit model to JET / KSTAR control room for shadow-mode validation
- Run in parallel with existing system for ≥100 disruptive shots
- Report prospective TPR/FPR without post-hoc threshold adjustment
"""


def save_validation_protocol(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "validation_protocol.md"
    path.write_text(VALIDATION_PROTOCOL)
    logger.info(f"Validation protocol saved → {path}")
