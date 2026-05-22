from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def setup_style() -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "legend.fontsize": 10,
        }
    )


def save_figure(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def draw_box(ax, xy, width, height, text, color):
    patch = FancyBboxPatch(xy, width, height, boxstyle="round,pad=0.02,rounding_size=0.03", fc=color, ec="#2f2f2f", lw=1.5)
    ax.add_patch(patch)
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center", color="#1c1c1c", weight="bold")


def arrow(ax, start, end, text=None):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=16, lw=1.6, color="#4c4c4c"))
    if text:
        mid = ((start[0] + end[0]) / 2, (start[1] + end[1]) / 2)
        ax.text(mid[0], mid[1] + 0.03, text, ha="center", va="bottom", fontsize=10)


def figure_system_architecture() -> None:
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    colors = cm.cividis(np.linspace(0.18, 0.88, 6))
    draw_box(ax, (0.05, 0.68), 0.22, 0.14, "State Representation\n(mesh / particle / latent)", colors[0])
    draw_box(ax, (0.38, 0.68), 0.20, 0.14, "Planner\n(MPC / CEM)", colors[1])
    draw_box(ax, (0.70, 0.68), 0.22, 0.14, "Reactive Controller\n(visual feedback)", colors[2])
    draw_box(ax, (0.05, 0.30), 0.22, 0.14, "Physics Simulator\n(forward model)", colors[3])
    draw_box(ax, (0.38, 0.30), 0.20, 0.14, "Sim2Real Bridge\n(domain randomization)", colors[4])
    draw_box(ax, (0.70, 0.30), 0.22, 0.14, "Cloth Folding Task\n(application)", colors[5])
    arrow(ax, (0.27, 0.75), (0.38, 0.75), "state features")
    arrow(ax, (0.58, 0.75), (0.70, 0.75), "action sequence")
    arrow(ax, (0.27, 0.37), (0.38, 0.37), "sim traces")
    arrow(ax, (0.58, 0.37), (0.70, 0.37), "deployment")
    arrow(ax, (0.16, 0.68), (0.16, 0.44), "encoded state")
    arrow(ax, (0.48, 0.68), (0.48, 0.44), "planned rollouts")
    arrow(ax, (0.81, 0.68), (0.81, 0.44), "feedback actions")
    arrow(ax, (0.81, 0.44), (0.48, 0.68), "visual loop")
    ax.set_title("Cloth Folding System Architecture")
    save_figure(fig, "system_architecture")


def figure_state_representations() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))
    x = np.linspace(-1, 1, 8)
    y = np.linspace(-1, 1, 8)
    xx, yy = np.meshgrid(x, y)
    warp = 0.12 * np.sin(np.pi * xx) * np.cos(np.pi * yy)
    palette = cm.viridis(np.linspace(0.2, 0.85, 8))

    ax = axes[0]
    for i in range(xx.shape[0]):
        ax.plot(xx[i], yy[i] + warp[i], color=palette[i], lw=1.6)
        ax.plot(xx[:, i], yy[:, i] + warp[:, i], color=palette[i], lw=1.2, alpha=0.85)
    ax.set_title("Mesh Representation")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal")

    ax = axes[1]
    scatter = ax.scatter(xx.ravel(), (yy + warp).ravel(), c=(yy + warp).ravel(), cmap="cividis", s=42, edgecolor="black", linewidth=0.25)
    fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04, label="Height")
    ax.set_title("Particle Representation")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal")

    ax = axes[2]
    latent = np.outer(np.linspace(0.2, 1.0, 8), np.linspace(1.0, 0.25, 8))
    latent += 0.06 * np.sin(np.arange(64).reshape(8, 8) / 2.0)
    im = ax.imshow(latent, cmap="viridis", aspect="auto")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Activation")
    ax.set_title("Latent Representation")
    ax.set_xlabel("Latent Dimension")
    ax.set_ylabel("Channel")

    save_figure(fig, "state_representations")


def figure_planning_comparison() -> None:
    methods = ["MPC", "CEM", "MPPI", "Graph", "RL"]
    success = np.array([0.89, 0.84, 0.82, 0.76, 0.79])
    planning_time = np.array([0.42, 0.37, 0.46, 0.18, 0.06])
    smoothness = np.array([0.83, 0.71, 0.74, 0.61, 0.66])
    colors = cm.cividis(np.linspace(0.2, 0.85, len(methods)))
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6))
    for ax, values, title, ylabel in zip(
        axes,
        [success, planning_time, smoothness],
        ["Success Rate", "Planning Time", "Trajectory Smoothness"],
        ["Rate", "Seconds", "Score"],
    ):
        ax.bar(methods, values, color=colors)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x", rotation=20)
    axes[0].set_ylim(0, 1.0)
    axes[2].set_ylim(0, 1.0)
    save_figure(fig, "planning_comparison")


def figure_sim2real_pipeline() -> None:
    fig, ax = plt.subplots(figsize=(12, 4.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    colors = cm.viridis(np.linspace(0.15, 0.85, 5))
    labels = [
        "Simulation Assets",
        "Domain Randomization",
        "Learned Dynamics",
        "Reactive Policy",
        "Real Robot Deployment",
    ]
    xs = np.linspace(0.04, 0.8, len(labels))
    for idx, (x, label) in enumerate(zip(xs, labels)):
        draw_box(ax, (x, 0.38), 0.15, 0.22, label, colors[idx])
        if idx < len(labels) - 1:
            arrow(ax, (x + 0.15, 0.49), (xs[idx + 1], 0.49))
    ax.text(0.335, 0.24, "texture • lighting • material • sensor noise", ha="center", va="center", fontsize=10)
    ax.set_title("Sim-to-Real Transfer Pipeline")
    save_figure(fig, "sim2real_pipeline")


def cloth_polygon(stage: int) -> np.ndarray:
    base = np.array([[-1.0, -0.75], [1.0, -0.75], [1.0, 0.75], [-1.0, 0.75]])
    if stage == 0:
        return base
    if stage == 1:
        return np.array([[-1.0, -0.75], [0.0, -0.75], [0.35, 0.95], [-0.9, 0.85]])
    if stage == 2:
        return np.array([[-1.0, -0.75], [0.0, -0.75], [0.25, 1.05], [-0.8, 0.92]])
    if stage == 3:
        return np.array([[-1.0, -0.75], [0.0, -0.75], [-0.15, 0.76], [-1.0, 0.75]])
    if stage == 4:
        return np.array([[-1.0, -0.75], [0.0, -0.75], [-0.08, 0.75], [-0.92, 0.75]])
    return np.array([[-1.0, -0.75], [0.0, -0.75], [0.0, 0.75], [-1.0, 0.75]])


def figure_cloth_folding_sequence() -> None:
    labels = ["Initial", "Pick", "Lift", "Fold", "Release", "Final"]
    fig, axes = plt.subplots(1, 6, figsize=(15, 3.4))
    colors = cm.cividis(np.linspace(0.2, 0.85, 6))
    for idx, ax in enumerate(axes):
        poly = Polygon(cloth_polygon(idx), closed=True, facecolor=colors[idx], edgecolor="#333333", lw=1.5, alpha=0.95)
        ax.add_patch(poly)
        if idx in {1, 2}:
            ax.add_patch(Circle((0.05, 0.72), 0.08, facecolor="#2f2f2f", edgecolor="white", lw=0.8))
        ax.set_xlim(-1.15, 1.15)
        ax.set_ylim(-1.0, 1.1)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(labels[idx])
    save_figure(fig, "cloth_folding_sequence")


def figure_training_curves() -> None:
    episodes = np.arange(1, 121)
    rng = np.random.default_rng(7)
    loss = 0.6 * np.exp(-episodes / 28) + 0.03 + rng.normal(0.0, 0.01, size=episodes.size)
    reward = -1.8 + 1.7 * (1 - np.exp(-episodes / 35)) + rng.normal(0.0, 0.04, size=episodes.size)
    success = 0.25 + 0.68 * (1 - np.exp(-episodes / 42)) + rng.normal(0.0, 0.02, size=episodes.size)
    curves = [loss, reward, np.clip(success, 0, 1)]
    titles = ["Dynamics Model Loss", "RL Reward", "Planning Success Rate"]
    ylabels = ["Loss", "Reward", "Success Rate"]
    colors = [cm.viridis(0.2), cm.viridis(0.52), cm.viridis(0.82)]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    for ax, curve, title, ylabel, color in zip(axes, curves, titles, ylabels, colors):
        variance = 0.06 * np.maximum(np.abs(curve), 0.1)
        ax.plot(episodes, curve, color=color, lw=2.2)
        ax.fill_between(episodes, curve - variance, curve + variance, color=color, alpha=0.18)
        ax.set_title(title)
        ax.set_xlabel("Episode")
        ax.set_ylabel(ylabel)
    axes[2].set_ylim(0, 1.0)
    save_figure(fig, "training_curves")


def figure_domain_randomization_ablation() -> None:
    conditions = ["None", "Texture", "Lighting", "Material", "Dynamics", "Full"]
    values = np.array([0.52, 0.61, 0.64, 0.70, 0.74, 0.87])
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    bars = ax.bar(conditions, values, color=cm.cividis(np.linspace(0.2, 0.9, len(conditions))))
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Transfer Success Rate")
    ax.set_title("Domain Randomization Ablation")
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.02, f"{value:.2f}", ha="center", va="bottom")
    save_figure(fig, "domain_randomization_ablation")


def main() -> None:
    setup_style()
    figure_system_architecture()
    figure_state_representations()
    figure_planning_comparison()
    figure_sim2real_pipeline()
    figure_cloth_folding_sequence()
    figure_training_curves()
    figure_domain_randomization_ablation()


if __name__ == "__main__":
    main()
