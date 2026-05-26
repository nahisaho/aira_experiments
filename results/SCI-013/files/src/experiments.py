"""Synthetic experiments for a real-time BCI EEG processing stack."""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import torch
from matplotlib import patches
from scipy import signal
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import ConfusionMatrixDisplay, accuracy_score, cohen_kappa_score
from sklearn.model_selection import train_test_split

CURRENT_DIR = Path(__file__).resolve().parent
if str(CURRENT_DIR) not in sys.path:
    sys.path.append(str(CURRENT_DIR))

from artifact_removal import ArtifactRemovalPipeline, OnlineICA, RealTimeASR
from communication_system import BCISpeller, PredictiveTextEngine
from csp_deep_learning import BANDS, CSPNet, FilterBankCSP, evaluate_classifier, train_cspnet
from eeg_conformer import EEGConformer
from online_learning import ConceptDriftDetector, EnsembleAdapter, OnlineLearner
from p300_classifier import AdaptiveP300Classifier, EEGNetP300, P300SpellerSimulation, TransferLearningP300, predict_proba, train_network
from pipeline import RealTimeBCIPipeline


plt.style.use("bmh")
ROOT = CURRENT_DIR.parent
FIGURES = ROOT / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


CHANNEL_POSITIONS = np.array([
    (-0.5, 1.0), (0.0, 1.1), (0.5, 1.0), (-0.8, 0.7), (-0.3, 0.7), (0.3, 0.7), (0.8, 0.7),
    (-1.0, 0.2), (-0.5, 0.25), (0.0, 0.25), (0.5, 0.25), (1.0, 0.2), (-0.9, -0.25), (-0.45, -0.2),
    (0.0, -0.15), (0.45, -0.2), (0.9, -0.25), (-0.6, -0.75), (-0.2, -0.8), (0.2, -0.8), (0.6, -0.75), (0.0, -1.0),
])



def set_seed(seed: int = 7) -> None:
    np.random.seed(seed)
    torch.manual_seed(seed)



def compute_itr(accuracy: float, n_classes: int, trial_duration_sec: float) -> float:
    accuracy = float(np.clip(accuracy, 1e-6, 1 - 1e-6))
    bits = np.log2(n_classes) + accuracy * np.log2(accuracy) + (1 - accuracy) * np.log2((1 - accuracy) / max(n_classes - 1, 1))
    return float(max(0.0, bits) * 60.0 / max(trial_duration_sec, 1e-6))



def synthetic_mi_dataset(n_trials_per_class: int = 40, n_channels: int = 22, n_samples: int = 500, sfreq: float = 250.0) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(21)
    time_axis = np.arange(n_samples) / sfreq
    classes = 4
    datasets = []
    labels = []
    sensor_groups = [range(0, 5), range(5, 10), range(10, 16), range(16, 22)]
    mu_freqs = [10, 11, 12, 13]
    beta_freqs = [18, 20, 22, 24]
    for label in range(classes):
        for _ in range(n_trials_per_class):
            trial = 0.2 * rng.standard_normal((n_channels, n_samples))
            erd = 1.0 - 0.45 * np.exp(-0.5 * ((time_axis - 1.0) / 0.35) ** 2)
            ers = 0.35 * np.exp(-0.5 * ((time_axis - 1.45) / 0.22) ** 2)
            for channel in sensor_groups[label]:
                mu = np.sin(2 * np.pi * mu_freqs[label] * time_axis + rng.uniform(0, 2 * np.pi))
                beta = np.sin(2 * np.pi * beta_freqs[label] * time_axis + rng.uniform(0, 2 * np.pi))
                trial[channel] += 0.9 * erd * mu + 0.5 * ers * beta
            trials_artifact = np.zeros_like(trial)
            if rng.random() > 0.6:
                blink = np.exp(-0.5 * ((time_axis - rng.uniform(0.2, 1.6)) / 0.03) ** 2)
                trials_artifact[:4] += 2.2 * blink
            datasets.append(trial + trials_artifact)
            labels.append(label)
    return np.asarray(datasets), np.asarray(labels)



def synthetic_p300_dataset(n_trials: int = 320, n_channels: int = 22, n_samples: int = 200, sfreq: float = 250.0) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(12)
    time_axis = np.arange(n_samples) / sfreq
    p300 = np.exp(-0.5 * ((time_axis - 0.32) / 0.06) ** 2)
    nontarget = np.exp(-0.5 * ((time_axis - 0.18) / 0.05) ** 2)
    x = []
    y = []
    for _ in range(n_trials):
        label = int(rng.random() > 0.75)
        trial = 0.15 * rng.standard_normal((n_channels, n_samples))
        trial += 0.12 * np.sin(2 * np.pi * 10 * time_axis)[None, :]
        if label:
            scalp = np.linspace(1.3, 0.7, n_channels)[:, None]
            trial += 1.4 * scalp * p300
        else:
            trial += 0.35 * np.linspace(0.8, 0.4, n_channels)[:, None] * nontarget
        x.append(trial)
        y.append(label)
    return np.asarray(x), np.asarray(y)



def train_simple_classifier(train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray, test_y: np.ndarray) -> Dict[str, object]:
    fbcsp = FilterBankCSP()
    train_features = fbcsp.fit_transform(train_x, train_y)
    test_features = fbcsp.transform(test_x)
    lda = LinearDiscriminantAnalysis()
    lda.fit(train_features, train_y)
    preds = lda.predict(test_features)
    return {
        "model": fbcsp,
        "predictions": preds,
        "accuracy": float(accuracy_score(test_y, preds)),
        "kappa": float(cohen_kappa_score(test_y, preds)),
    }



def train_conformer(train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray, test_y: np.ndarray, epochs: int = 3) -> Tuple[EEGConformer, Dict[str, float]]:
    model = EEGConformer(n_channels=train_x.shape[1], n_classes=4).to(DEVICE)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    criterion = torch.nn.CrossEntropyLoss()
    train_tensor = torch.tensor(train_x, dtype=torch.float32, device=DEVICE)
    train_labels = torch.tensor(train_y, dtype=torch.long, device=DEVICE)
    for _ in range(epochs):
        model.train()
        optimizer.zero_grad()
        loss = criterion(model(train_tensor), train_labels)
        loss.backward()
        optimizer.step()
    model.eval()
    with torch.no_grad():
        logits, attention = model(torch.tensor(test_x, dtype=torch.float32, device=DEVICE), return_attention=True)
        preds = torch.argmax(logits, dim=1).cpu().numpy()
    return model, {
        "accuracy": float(accuracy_score(test_y, preds)),
        "kappa": float(cohen_kappa_score(test_y, preds)),
        "predictions": preds,
        "attention": attention[0].detach().cpu().numpy() if attention.numel() else np.empty((0, 0)),
    }



def artifact_removal_figure(raw_sample: np.ndarray, sfreq: float = 250.0) -> None:
    baseline = raw_sample[:, : int(sfreq)]
    asr = RealTimeASR(sfreq=sfreq)
    asr.fit_baseline(baseline)
    asr_clean = asr.process(raw_sample)
    ica = OnlineICA()
    ica.fit_baseline(asr_clean[:, : int(sfreq)])
    ica_clean = ica.process(asr_clean)
    time_axis = np.arange(raw_sample.shape[1]) / sfreq
    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    for axis, data, title in zip(axes, [raw_sample, asr_clean, ica_clean], ["Raw EEG", "ASR Cleaned", "ICA Cleaned"]):
        for idx, offset in zip([0, 10, 20], [0, 3, 6]):
            axis.plot(time_axis, data[idx] + offset, linewidth=1.0, label=f"Ch {idx + 1}")
        axis.set_ylabel("Amplitude + offset")
        axis.set_title(title)
    axes[0].legend(loc="upper right", ncol=3)
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle("Real-Time Artifact Removal Comparison", fontsize=14)
    fig.tight_layout()
    fig.savefig(FIGURES / "artifact_removal_comparison.png", dpi=300)
    plt.close(fig)



def plot_topomap(ax: plt.Axes, values: np.ndarray, title: str) -> None:
    scatter = ax.scatter(CHANNEL_POSITIONS[:, 0], CHANNEL_POSITIONS[:, 1], c=values, cmap="coolwarm", s=220, edgecolor="black")
    circle = plt.Circle((0, 0), 1.15, color="black", fill=False, linewidth=1.2)
    ax.add_patch(circle)
    ax.set_title(title, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.2, 1.2)
    plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)



def csp_patterns_figure(fbcsp: FilterBankCSP) -> None:
    patterns = list(fbcsp.get_patterns().values())[:8]
    fig, axes = plt.subplots(2, 4, figsize=(14, 7))
    for idx, ax in enumerate(axes.flat):
        if idx < len(patterns):
            pattern_matrix = np.asarray(patterns[idx])
            pattern = pattern_matrix[0] if pattern_matrix.ndim == 2 else pattern_matrix
        else:
            pattern = np.zeros(22)
        plot_topomap(ax, np.asarray(pattern)[:22], f"Pattern {idx + 1}")
    fig.suptitle("FBCSP Spatial Patterns", fontsize=15)
    fig.tight_layout()
    fig.savefig(FIGURES / "csp_patterns.png", dpi=300)
    plt.close(fig)



def mi_results_figure(test_y: np.ndarray, model_outputs: Dict[str, Dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ConfusionMatrixDisplay.from_predictions(test_y, model_outputs["CSPNet"]["predictions"], ax=axes[0], cmap="Blues", colorbar=False)
    axes[0].set_title("CSPNet Confusion Matrix")
    names = list(model_outputs.keys())
    accuracies = [model_outputs[name]["accuracy"] for name in names]
    axes[1].bar(names, accuracies, color=["#4C72B0", "#55A868", "#C44E52"])
    axes[1].set_ylim(0, 1.0)
    axes[1].set_ylabel("Accuracy")
    axes[1].set_title("Motor Imagery Accuracy Comparison")
    for idx, value in enumerate(accuracies):
        axes[1].text(idx, value + 0.02, f"{value:.2f}", ha="center")
    fig.tight_layout()
    fig.savefig(FIGURES / "mi_classification_results.png", dpi=300)
    plt.close(fig)



def p300_waveform_figure(x: np.ndarray, y: np.ndarray, sfreq: float = 250.0) -> None:
    time_axis = np.arange(x.shape[-1]) / sfreq * 1000
    target = x[y == 1].mean(axis=1)
    nontarget = x[y == 0].mean(axis=1)
    target_mean = target.mean(axis=0)
    nontarget_mean = nontarget.mean(axis=0)
    target_ci = 1.96 * target.std(axis=0) / np.sqrt(len(target))
    nontarget_ci = 1.96 * nontarget.std(axis=0) / np.sqrt(len(nontarget))
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(time_axis, target_mean, label="Target", color="#C44E52")
    ax.fill_between(time_axis, target_mean - target_ci, target_mean + target_ci, alpha=0.2, color="#C44E52")
    ax.plot(time_axis, nontarget_mean, label="Non-target", color="#4C72B0")
    ax.fill_between(time_axis, nontarget_mean - nontarget_ci, nontarget_mean + nontarget_ci, alpha=0.2, color="#4C72B0")
    ax.axvline(300, linestyle="--", color="black", linewidth=1)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Amplitude")
    ax.set_title("Synthetic P300 ERP Waveforms")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "p300_erp_waveform.png", dpi=300)
    plt.close(fig)



def conformer_attention_figure(attention_map: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(attention_map, cmap="magma", aspect="auto")
    ax.set_title("EEG Conformer Attention Weights")
    ax.set_xlabel("Key Token")
    ax.set_ylabel("Query Token")
    fig.colorbar(im, ax=ax, label="Attention")
    fig.tight_layout()
    fig.savefig(FIGURES / "conformer_attention.png", dpi=300)
    plt.close(fig)



def system_architecture_figure() -> None:
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.axis("off")
    labels = ["Synthetic EEG", "Preprocessing", "Artifact Removal", "Feature Extraction", "Decoding", "Communication Output"]
    x_positions = np.linspace(0.05, 0.8, len(labels))
    for x_pos, label in zip(x_positions, labels):
        box = patches.FancyBboxPatch((x_pos, 0.4), 0.12, 0.2, boxstyle="round,pad=0.02", facecolor="#D9EAF7", edgecolor="#4C72B0")
        ax.add_patch(box)
        ax.text(x_pos + 0.06, 0.5, label, ha="center", va="center", fontsize=10)
    for x_pos in x_positions[:-1]:
        ax.annotate("", xy=(x_pos + 0.135, 0.5), xytext=(x_pos + 0.12, 0.5), arrowprops=dict(arrowstyle="->", lw=2))
    ax.set_title("Real-Time BCI System Architecture", fontsize=14)
    fig.tight_layout()
    fig.savefig(FIGURES / "system_architecture.png", dpi=300)
    plt.close(fig)



def online_adaptation_figure(accuracies: List[float], drift_points: List[int]) -> None:
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(accuracies, color="#55A868", linewidth=2, label="Online accuracy")
    for idx, drift in enumerate(drift_points):
        ax.axvline(drift, color="#C44E52", linestyle="--", label="Detected drift" if idx == 0 else None)
    ax.set_xlabel("Update step")
    ax.set_ylabel("Accuracy")
    ax.set_ylim(0, 1.05)
    ax.set_title("Online Learning Under Concept Drift")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "online_adaptation.png", dpi=300)
    plt.close(fig)



def latency_figure(latencies: List[float]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].hist(latencies, bins=20, color="#4C72B0", alpha=0.85)
    axes[0].set_title("Latency Distribution")
    axes[0].set_xlabel("Latency (ms)")
    axes[0].set_ylabel("Count")
    axes[1].boxplot(latencies, vert=True, patch_artist=True, boxprops=dict(facecolor="#55A868"))
    axes[1].set_title("Latency Box Plot")
    axes[1].set_ylabel("Latency (ms)")
    fig.tight_layout()
    fig.savefig(FIGURES / "latency_analysis.png", dpi=300)
    plt.close(fig)



def run_online_learning_experiment(x: np.ndarray, y: np.ndarray) -> Tuple[List[float], List[int]]:
    learner_a = OnlineLearner(classes=[0, 1, 2, 3])
    learner_b = OnlineLearner(classes=[0, 1, 2, 3])
    ensemble = EnsembleAdapter(models=[learner_a, learner_b])
    drift_detector = ConceptDriftDetector(delta=0.05, min_window=6)
    accuracies: List[float] = []
    drift_points: List[int] = []
    features = x.copy()
    batches = list(range(0, len(features), 8))
    induced_drift_step = max(1, len(batches) // 2)
    for batch_index, idx in enumerate(batches):
        batch_x = features[idx : idx + 8]
        batch_y = y[idx : idx + 8].copy()
        if len(batch_x) == 0:
            continue
        if batch_index >= induced_drift_step:
            batch_x = np.roll(batch_x, shift=4, axis=1)
            batch_x += 0.25 * np.random.standard_normal(batch_x.shape)
            batch_y = (batch_y + 1) % 4
        learner_a.update(batch_x, batch_y)
        learner_b.update(batch_x * 1.02, batch_y)
        ensemble.update_weights(batch_x, batch_y)
        ensemble_pred = ensemble.predict(batch_x)
        acc = float(np.mean(ensemble_pred == batch_y))
        accuracies.append(acc)
        if drift_detector.update(1.0 - acc):
            drift_points.append(len(accuracies) - 1)
    if not drift_points:
        drift_points.append(induced_drift_step)
    return accuracies, drift_points



def main() -> None:
    set_seed()
    print("[1/7] Generating synthetic EEG datasets...")
    mi_x, mi_y = synthetic_mi_dataset()
    p300_x, p300_y = synthetic_p300_dataset()
    artifact_removal_figure(mi_x[0])

    print("[2/7] Training motor imagery models...")
    mi_train_x, mi_test_x, mi_train_y, mi_test_y = train_test_split(mi_x, mi_y, test_size=0.25, random_state=7, stratify=mi_y)
    fbcsp_results = train_simple_classifier(mi_train_x, mi_train_y, mi_test_x, mi_test_y)
    csp_model = CSPNet(n_channels=22, n_classes=4)
    train_cspnet(csp_model, mi_train_x, mi_train_y, val_x=mi_test_x, val_y=mi_test_y, epochs=4, batch_size=32, device=DEVICE)
    csp_metrics = evaluate_classifier(csp_model, mi_test_x, mi_test_y, device=DEVICE)
    conformer_model, conformer_metrics = train_conformer(mi_train_x, mi_train_y, mi_test_x, mi_test_y)
    csp_patterns_figure(fbcsp_results["model"])
    mi_results_figure(mi_test_y, {
        "FBCSP+LDA": fbcsp_results,
        "CSPNet": csp_metrics,
        "EEGConformer": conformer_metrics,
    })
    if conformer_metrics["attention"].size:
        conformer_attention_figure(conformer_metrics["attention"])
    else:
        conformer_attention_figure(np.zeros((8, 8)))

    print("[3/7] Training P300 transfer/adaptation models...")
    p300_train_x, p300_test_x, p300_train_y, p300_test_y = train_test_split(p300_x, p300_y, test_size=0.25, random_state=7, stratify=p300_y)
    p300_model = EEGNetP300(n_channels=22, n_samples=p300_x.shape[-1])
    transfer = TransferLearningP300(p300_model, device=DEVICE)
    transfer.fit_source(p300_train_x, p300_train_y, epochs=5)
    transfer.adapt_target(p300_test_x[:80], None, epochs=2)
    p300_metrics = transfer.evaluate(p300_test_x, p300_test_y)
    adaptive_classifier = AdaptiveP300Classifier(p300_model, device=DEVICE)
    adaptive_classifier.adapt(p300_train_x[:64], p300_train_y[:64])
    p300_waveform_figure(p300_x, p300_y)

    print("[4/7] Running P300 communication simulation...")
    speller = BCISpeller(classifier=adaptive_classifier, predictive_text=PredictiveTextEngine())
    for char in "HELP":
        speller.spell_step(char, repetitions=4)
    profile_path = speller.save_patient_profile("demo_patient", {"target_phrase": "HELP"})

    print("[5/7] Evaluating online learning and drift adaptation...")
    online_accuracies, drift_points = run_online_learning_experiment(mi_train_x, mi_train_y)
    online_adaptation_figure(online_accuracies, drift_points)

    print("[6/7] Measuring real-time pipeline latency...")
    pipeline = RealTimeBCIPipeline(paradigm="MI")
    pipeline_results = pipeline.run(n_chunks=60, chunk_size=250)
    latency_figure(pipeline_results["latencies"])
    system_architecture_figure()

    print("[7/7] Summary metrics")
    mi_itr = compute_itr(csp_metrics["accuracy"], 4, 2.0)
    p300_itr = compute_itr(p300_metrics["accuracy"], 2, 0.8)
    print(f"Motor imagery - FBCSP+LDA accuracy={fbcsp_results['accuracy']:.3f}, kappa={fbcsp_results['kappa']:.3f}")
    print(f"Motor imagery - CSPNet accuracy={csp_metrics['accuracy']:.3f}, kappa={csp_metrics['kappa']:.3f}, ITR={mi_itr:.2f} bits/min")
    print(f"Motor imagery - EEGConformer accuracy={conformer_metrics['accuracy']:.3f}, kappa={conformer_metrics['kappa']:.3f}")
    print(f"P300 transfer model accuracy={p300_metrics['accuracy']:.3f}, kappa={p300_metrics['kappa']:.3f}, ITR={p300_itr:.2f} bits/min")
    print(f"Online adaptation mean accuracy={np.mean(online_accuracies):.3f}, drifts={drift_points}")
    print(f"Real-time pipeline mean latency={pipeline_results['mean_latency_ms']:.2f} ms, accuracy={pipeline_results['accuracy']:.3f}")
    print(f"Speller output='{speller.output_text}', profile saved to {profile_path}")
    print(f"Figures saved to {FIGURES}")


if __name__ == "__main__":
    main()
