"""
Visualization module for NCC Framework.
All figure text is in English (journal-ready).
Uses colorblind-friendly palettes (viridis, cividis, colorbrewer).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.patheffects as pe
from typing import List, Dict, Optional

# --- Colorblind-friendly palette (Wong 2011) ---
WONG = {
    "black":   "#000000",
    "orange":  "#E69F00",
    "sky":     "#56B4E9",
    "green":   "#009E73",
    "yellow":  "#F0E442",
    "blue":    "#0072B2",
    "vermillion": "#D55E00",
    "pink":    "#CC79A7",
}
PALETTE = list(WONG.values())

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def plot_phi_vs_consciousness(
    phi_values: np.ndarray,
    consciousness_levels: np.ndarray,
    save_path: str,
    phi_std: Optional[np.ndarray] = None,
) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    color = WONG["blue"]

    ax.plot(consciousness_levels, phi_values, "o-", color=color,
            linewidth=2, markersize=6, label="Φ (IIT)")
    if phi_std is not None:
        ax.fill_between(consciousness_levels,
                        phi_values - phi_std,
                        phi_values + phi_std,
                        alpha=0.2, color=color)

    ax.set_xlabel("Consciousness Level (0=deep anesthesia, 1=awake)")
    ax.set_ylabel("Integrated Information Φ (bits)")
    ax.set_title("IIT: Φ as a Function of Consciousness Level")
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Annotate regions
    ax.axvspan(0, 0.3, alpha=0.06, color="red", label="Anesthesia")
    ax.axvspan(0.3, 0.6, alpha=0.06, color="yellow")
    ax.axvspan(0.6, 1.0, alpha=0.06, color="green")
    ax.text(0.15, ax.get_ylim()[1] * 0.92, "Anesthesia", ha="center",
            fontsize=9, color="red", style="italic")
    ax.text(0.45, ax.get_ylim()[1] * 0.92, "Transition", ha="center",
            fontsize=9, color="#b07000", style="italic")
    ax.text(0.80, ax.get_ylim()[1] * 0.92, "Awake", ha="center",
            fontsize=9, color="green", style="italic")

    fig.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)


def plot_pci_spectrum(
    pci_results: list,
    save_path: str,
) -> None:
    levels = [r["consciousness_level"] for r in pci_results]
    pcis = [r["pci"] for r in pci_results]
    stds = [r["pci_std"] for r in pci_results]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # PCI vs level
    ax = axes[0]
    ax.errorbar(levels, pcis, yerr=stds, fmt="s-",
                color=WONG["orange"], linewidth=2, markersize=7,
                capsize=4, ecolor=WONG["vermillion"], label="PCI")
    ax.axhline(0.31, color="gray", linestyle="--", alpha=0.7,
               label="PCI threshold (0.31)")
    ax.set_xlabel("Consciousness Level")
    ax.set_ylabel("Perturbational Complexity Index (PCI)")
    ax.set_title("PCI vs. Consciousness Level")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Clinical reference bands
    clinical_refs = {
        "CTRL": (0.9, 1.0, WONG["green"]),
        "MCS+": (0.5, 0.65, WONG["sky"]),
        "MCS-": (0.3, 0.45, WONG["yellow"]),
        "VS":   (0.05, 0.25, WONG["vermillion"]),
    }
    ax2 = axes[1]
    for state, (lo, hi, col) in clinical_refs.items():
        mask = [(lo <= l <= hi) for l in levels]
        ls = [l for l, m in zip(levels, mask) if m]
        ps = [p for p, m in zip(pcis, mask) if m]
        ax2.scatter(ls, ps, color=col, s=80, label=state, zorder=5)

    ax2.plot(levels, pcis, "-", color="gray", alpha=0.4, linewidth=1)
    ax2.axhline(0.31, color="gray", linestyle="--", alpha=0.7)
    ax2.set_xlabel("Consciousness Level")
    ax2.set_ylabel("PCI")
    ax2.set_title("PCI with Clinical State Reference")
    ax2.legend(title="Clinical State", fontsize=9)
    ax2.grid(True, alpha=0.3)

    fig.suptitle("Perturbational Complexity Index (PCI) Simulation", fontweight="bold")
    fig.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)


def plot_gwt_metrics(
    gwt_results: list,
    consciousness_levels: np.ndarray,
    save_path: str,
) -> None:
    metrics_to_plot = [
        ("global_efficiency", "Global Efficiency", WONG["blue"]),
        ("clustering_coefficient", "Clustering Coefficient", WONG["orange"]),
        ("small_world_index", "Small World Index (σ)", WONG["green"]),
        ("ignition_index", "Ignition Index", WONG["pink"]),
        ("information_broadcast_capacity", "Broadcast Capacity", WONG["sky"]),
        ("gwt_index", "GWT Consciousness Index", WONG["vermillion"]),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    axes = axes.ravel()

    for ax, (key, label, color) in zip(axes, metrics_to_plot):
        values = [r.get(key, 0) for r in gwt_results]
        ax.plot(consciousness_levels, values, "o-", color=color,
                linewidth=2, markersize=5)
        ax.set_xlabel("Consciousness Level")
        ax.set_ylabel(label)
        ax.set_title(label)
        ax.grid(True, alpha=0.3)

        # Fit trend line
        try:
            coef = np.polyfit(consciousness_levels, values, 1)
            trend = np.poly1d(coef)
            ax.plot(consciousness_levels, trend(consciousness_levels),
                    "--", color="gray", alpha=0.5)
            r = np.corrcoef(consciousness_levels, values)[0, 1]
            ax.text(0.05, 0.93, f"r = {r:.2f}", transform=ax.transAxes,
                    fontsize=9, color="gray")
        except Exception:
            pass

    fig.suptitle("Global Workspace Theory Metrics vs. Consciousness",
                 fontweight="bold", fontsize=13)
    fig.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)


def plot_clinical_features(
    X: np.ndarray,
    y: np.ndarray,
    feature_names: List[str],
    class_names: Dict[int, str],
    save_path: str,
) -> None:
    """LDA scatter plot + feature importance for clinical classification."""
    from sklearn.preprocessing import StandardScaler
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    lda = LinearDiscriminantAnalysis(n_components=2)
    X_lda = lda.fit_transform(X_scaled, y)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # LDA scatter
    ax = axes[0]
    unique_classes = np.unique(y)
    colors = [PALETTE[i % len(PALETTE)] for i in range(len(unique_classes))]
    markers = ["o", "s", "^", "D", "v"]

    for cls, col, mk in zip(unique_classes, colors, markers):
        mask = y == cls
        ax.scatter(X_lda[mask, 0], X_lda[mask, 1],
                   c=col, marker=mk, s=60, alpha=0.8,
                   label=class_names.get(cls, str(cls)), edgecolors="white",
                   linewidths=0.5)
    ax.set_xlabel(f"LD1 ({100*lda.explained_variance_ratio_[0]:.1f}%)")
    ax.set_ylabel(f"LD2 ({100*lda.explained_variance_ratio_[1]:.1f}%)")
    ax.set_title("LDA Projection of Clinical EEG Features")
    ax.legend(title="Clinical State", framealpha=0.9)
    ax.grid(True, alpha=0.2)

    # Feature importance (via LDA coefficients norm)
    ax2 = axes[1]
    importance = np.abs(lda.coef_).mean(axis=0)
    sorted_idx = np.argsort(importance)[::-1][:10]
    top_names = [feature_names[i] for i in sorted_idx]
    top_vals = importance[sorted_idx]

    bars = ax2.barh(range(len(top_names)), top_vals,
                    color=[PALETTE[i % len(PALETTE)] for i in range(len(top_names))])
    ax2.set_yticks(range(len(top_names)))
    ax2.set_yticklabels(top_names)
    ax2.set_xlabel("Mean |LDA Coefficient|")
    ax2.set_title("Top 10 Discriminative Features")
    ax2.grid(True, alpha=0.2, axis="x")

    fig.suptitle("Clinical Consciousness Classification", fontweight="bold")
    fig.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)


def plot_multi_index_comparison(
    consciousness_levels: np.ndarray,
    phi_values: np.ndarray,
    pci_values: np.ndarray,
    gwt_indices: np.ndarray,
    save_path: str,
) -> None:
    """Combined plot of Φ, PCI, and GWT index vs consciousness level."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Normalize to [0, 1] for comparison
    def normalize(arr):
        lo, hi = arr.min(), arr.max()
        if hi == lo:
            return np.zeros_like(arr)
        return (arr - lo) / (hi - lo)

    phi_n = normalize(phi_values)
    pci_n = normalize(pci_values)
    gwt_n = normalize(gwt_indices)

    ax.plot(consciousness_levels, phi_n, "o-", color=WONG["blue"],
            label="IIT Φ (normalized)", linewidth=2, markersize=5)
    ax.plot(consciousness_levels, pci_n, "s-", color=WONG["orange"],
            label="PCI (normalized)", linewidth=2, markersize=5)
    ax.plot(consciousness_levels, gwt_n, "^-", color=WONG["green"],
            label="GWT Index (normalized)", linewidth=2, markersize=5)

    ax.set_xlabel("Consciousness Level (0=anesthesia, 1=awake)")
    ax.set_ylabel("Normalized Index Value")
    ax.set_title("Multi-Theory Consciousness Indices: Comparative Overview")
    ax.legend(framealpha=0.9)
    ax.grid(True, alpha=0.3)

    # Mark clinical states
    clinical_markers = {
        "VS": 0.15, "MCS-": 0.35, "MCS+": 0.55, "CTRL": 1.0
    }
    for label, level in clinical_markers.items():
        ax.axvline(level, color="gray", linestyle=":", alpha=0.5)
        ax.text(level + 0.01, 0.05, label, fontsize=9, color="gray", rotation=90)

    fig.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)


def plot_confusion_matrix(
    cm: np.ndarray,
    class_names: List[str],
    save_path: str,
    title: str = "Confusion Matrix",
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(cm, interpolation="nearest", cmap="Blues")
    plt.colorbar(im, ax=ax)

    ticks = range(len(class_names))
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(class_names, rotation=45, ha="right")
    ax.set_yticklabels(class_names)

    thresh = cm.max() / 2
    for i in range(len(class_names)):
        for j in range(len(class_names)):
            ax.text(j, i, format(cm[i, j], ".2f"),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black",
                    fontsize=10)

    ax.set_xlabel("Predicted Label")
    ax.set_ylabel("True Label")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)


def plot_eeg_examples(
    consciousness_levels: List[float],
    n_channels: int = 4,
    n_samples: int = 512,
    fs: float = 256.0,
    save_path: str = "figures/eeg_examples.png",
    seed: int = 42,
) -> None:
    from .utils import generate_anesthesia_data

    n_levels = len(consciousness_levels)
    fig, axes = plt.subplots(n_levels, n_channels, figsize=(16, 3 * n_levels),
                             sharex=True, sharey="row")

    t = np.arange(n_samples) / fs
    labels = {1.0: "Awake", 0.55: "MCS+", 0.35: "MCS-", 0.15: "VS", 0.3: "Light Anes."}

    for row, level in enumerate(consciousness_levels):
        data = generate_anesthesia_data(n_channels=n_channels, n_samples=n_samples,
                                        consciousness_level=level, seed=seed + row)
        state_label = labels.get(level, f"Level={level:.2f}")
        for ch in range(n_channels):
            ax = axes[row, ch] if n_levels > 1 else axes[ch]
            col = PALETTE[(row * 2) % len(PALETTE)]
            ax.plot(t, data[ch], color=col, linewidth=0.8)
            if ch == 0:
                ax.set_ylabel(f"{state_label}\nAmplitude (a.u.)", fontsize=9)
            if row == n_levels - 1:
                ax.set_xlabel("Time (s)")
            ax.set_title(f"Ch {ch+1}", fontsize=9) if row == 0 else None
            ax.grid(True, alpha=0.2)
            ax.tick_params(labelsize=8)

    fig.suptitle("Simulated EEG at Different Consciousness Levels",
                 fontweight="bold", fontsize=13)
    fig.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)


def plot_phi_heatmap(
    phi_matrix: np.ndarray,
    save_path: str,
    title: str = "Pairwise Integrated Information Matrix",
) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(phi_matrix, cmap="viridis", aspect="auto")
    plt.colorbar(im, ax=ax, label="Φ (bits)")
    n = phi_matrix.shape[0]
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels([f"Ch{i+1}" for i in range(n)])
    ax.set_yticklabels([f"Ch{i+1}" for i in range(n)])
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)
