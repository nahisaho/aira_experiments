"""Generate all figures for the federated learning case study."""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import json
import os

plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "figure.dpi": 300,
    "savefig.dpi": 300,
})

COLORS = plt.cm.viridis(np.linspace(0.15, 0.85, 6))
os.makedirs("figures", exist_ok=True)

# Load results
with open("results/experiment_summary.json") as f:
    summary = json.load(f)

for name in ["strategy_fedavg", "strategy_fedprox", "strategy_scaffold"]:
    with open(f"results/{name}_history.json") as f:
        globals()[f"hist_{name.split('_',1)[1]}"] = json.load(f)


# ── Figure 1: Convergence comparison ──
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

for i, (metric, ylabel) in enumerate([
    ("global_loss", "Negative Partial Log-Likelihood"),
    ("c_index", "Concordance Index (C-index)"),
    ("beta_error", "Relative Parameter Error"),
]):
    ax = axes[i]
    for j, (name, label, ls) in enumerate([
        ("fedavg", "FedAvg", "-"),
        ("fedprox", "FedProx", "--"),
        ("scaffold", "SCAFFOLD", "-."),
    ]):
        hist = globals()[f"hist_{name}"]
        ax.plot(hist["round"], hist[metric], label=label,
                color=COLORS[j], linewidth=2, linestyle=ls)
    ax.set_xlabel("Communication Round")
    ax.set_ylabel(ylabel)
    ax.legend(frameon=True, fancybox=True)
    ax.grid(True, alpha=0.3)

fig.suptitle("Convergence Comparison: FedAvg vs FedProx vs SCAFFOLD", fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig("figures/fig1_convergence_comparison.png")
plt.savefig("figures/fig1_convergence_comparison.svg")
plt.close()
print("Saved fig1_convergence_comparison")


# ── Figure 2: Non-IID impact ──
fig, ax = plt.subplots(figsize=(7, 5))

non_iid_levels = [0.0, 0.5, 1.0, 2.0]
fedavg_scores = [summary[f"noniid_{n}_fedavg"]["final_c_index"] for n in non_iid_levels]
fedprox_scores = [summary[f"noniid_{n}_fedprox"]["final_c_index"] for n in non_iid_levels]

x = np.arange(len(non_iid_levels))
width = 0.35
bars1 = ax.bar(x - width/2, fedavg_scores, width, label="FedAvg",
               color=COLORS[0], edgecolor="white", linewidth=0.5)
bars2 = ax.bar(x + width/2, fedprox_scores, width, label="FedProx",
               color=COLORS[2], edgecolor="white", linewidth=0.5)

ax.set_xlabel("Non-IID Degree (σ)")
ax.set_ylabel("Concordance Index")
ax.set_title("Impact of Data Heterogeneity on Model Performance")
ax.set_xticks(x)
ax.set_xticklabels([str(n) for n in non_iid_levels])
ax.legend()
ax.grid(axis="y", alpha=0.3)
ax.set_ylim(0.80, 0.90)

for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f"{height:.4f}", xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8)

plt.tight_layout()
plt.savefig("figures/fig2_noniid_impact.png")
plt.savefig("figures/fig2_noniid_impact.svg")
plt.close()
print("Saved fig2_noniid_impact")


# ── Figure 3: Privacy-utility tradeoff ──
fig, ax1 = plt.subplots(figsize=(7, 5))

epsilons = [1.0, 5.0, 10.0, 50.0]
c_indices = [summary[f"dp_eps_{e}"]["final_c_index"] for e in epsilons]
privacy_spent = [summary[f"dp_eps_{e}"]["final_privacy_spent"] for e in epsilons]

# No-DP baseline
baseline = summary["strategy_fedavg"]["final_c_index"]

ax1.plot(epsilons, c_indices, "o-", color=COLORS[0], linewidth=2,
         markersize=8, label="C-index (with DP)", zorder=3)
ax1.axhline(y=baseline, color=COLORS[4], linestyle="--", linewidth=1.5,
            label=f"No-DP Baseline ({baseline:.4f})")
ax1.set_xlabel("Privacy Budget (ε)")
ax1.set_ylabel("Concordance Index", color=COLORS[0])
ax1.tick_params(axis="y", labelcolor=COLORS[0])
ax1.set_xscale("log")
ax1.set_ylim(0.55, 0.90)

ax2 = ax1.twinx()
ax2.bar([str(e) for e in epsilons],
        [summary[f"dp_eps_{e}"]["final_privacy_spent"] for e in epsilons],
        alpha=0.3, color=COLORS[3], label="Privacy Spent")
ax2.set_ylabel("Actual Privacy Spent (ε)", color=COLORS[3])
ax2.tick_params(axis="y", labelcolor=COLORS[3])

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc="lower right")

ax1.set_title("Privacy-Utility Tradeoff: Differential Privacy Impact")
ax1.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("figures/fig3_privacy_utility_tradeoff.png")
plt.savefig("figures/fig3_privacy_utility_tradeoff.svg")
plt.close()
print("Saved fig3_privacy_utility_tradeoff")


# ── Figure 4: Communication compression ──
fig, ax = plt.subplots(figsize=(7, 5))

comp_ratios = [0.01, 0.05, 0.1, 0.5, 1.0]
comp_c = [summary[f"compression_{c}"]["final_c_index"] for c in comp_ratios]

ax.plot(comp_ratios, comp_c, "s-", color=COLORS[1], linewidth=2, markersize=8)
for i, (cr, ci) in enumerate(zip(comp_ratios, comp_c)):
    savings = (1 - cr) * 100 if cr < 1.0 else 0
    ax.annotate(f"C={ci:.4f}\n{savings:.0f}% saved",
                xy=(cr, ci), xytext=(10, -20 if i % 2 == 0 else 10),
                textcoords="offset points", fontsize=8,
                arrowprops=dict(arrowstyle="->", color="gray", lw=0.5))

ax.set_xlabel("Compression Ratio (fraction of parameters transmitted)")
ax.set_ylabel("Concordance Index")
ax.set_title("Communication Efficiency: Top-K Sparsification Impact")
ax.set_xscale("log")
ax.grid(True, alpha=0.3)
ax.set_ylim(0.80, 0.90)
plt.tight_layout()
plt.savefig("figures/fig4_communication_efficiency.png")
plt.savefig("figures/fig4_communication_efficiency.svg")
plt.close()
print("Saved fig4_communication_efficiency")


# ── Figure 5: Byzantine resilience ──
fig, ax = plt.subplots(figsize=(9, 5.5))

byz_fracs = [0.0, 0.2, 0.4]
defenses = ["none", "krum", "median", "trimmed_mean"]
defense_labels = ["No Defense", "Krum", "Coord. Median", "Trimmed Mean"]

x = np.arange(len(byz_fracs))
width = 0.2

for i, (defense, label) in enumerate(zip(defenses, defense_labels)):
    scores = [summary[f"byzantine_{bf}_{defense}"]["final_c_index"] for bf in byz_fracs]
    bars = ax.bar(x + i * width, scores, width, label=label,
                  color=COLORS[i], edgecolor="white", linewidth=0.5)
    for bar, score in zip(bars, scores):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f"{score:.3f}", ha="center", fontsize=7, rotation=45)

ax.set_xlabel("Byzantine Client Fraction")
ax.set_ylabel("Concordance Index")
ax.set_title("Byzantine Resilience: Defense Strategy Comparison")
ax.set_xticks(x + width * 1.5)
ax.set_xticklabels([f"{int(bf*100)}%" for bf in byz_fracs])
ax.legend(loc="upper right")
ax.grid(axis="y", alpha=0.3)
ax.set_ylim(0.3, 0.95)
plt.tight_layout()
plt.savefig("figures/fig5_byzantine_resilience.png")
plt.savefig("figures/fig5_byzantine_resilience.svg")
plt.close()
print("Saved fig5_byzantine_resilience")


# ── Figure 6: Architecture diagram (text-based) ──
fig, ax = plt.subplots(figsize=(12, 8))
ax.set_xlim(0, 12)
ax.set_ylim(0, 8)
ax.axis("off")

# Title
ax.text(6, 7.7, "Federated Learning Platform Architecture", ha="center",
        fontsize=16, fontweight="bold")

# Server box
server_rect = plt.Rectangle((1, 3.5), 10, 3.8, fill=True, facecolor="#e8f4fd",
                              edgecolor="#2196F3", linewidth=2, zorder=1)
ax.add_patch(server_rect)
ax.text(6, 7.0, "FL Server (Flower)", ha="center", fontsize=13, fontweight="bold",
        color="#1565C0")

# Server components
components = [
    (2.5, 5.8, "Strategy\nSelector"),
    (5.0, 5.8, "DP Guard\n(RDP Accountant)"),
    (7.5, 5.8, "Byzantine\nFilter"),
    (10.0, 5.8, "Comm.\nOptimizer"),
]
for x_pos, y_pos, label in components:
    rect = plt.Rectangle((x_pos - 0.9, y_pos - 0.5), 1.8, 1.0, fill=True,
                          facecolor="#bbdefb", edgecolor="#1976D2", linewidth=1.5,
                          zorder=2)
    ax.add_patch(rect)
    ax.text(x_pos, y_pos, label, ha="center", va="center", fontsize=8)

# Aggregation engine
agg_rect = plt.Rectangle((2.0, 3.8), 8, 1.0, fill=True,
                           facecolor="#c8e6c9", edgecolor="#388E3C", linewidth=1.5, zorder=2)
ax.add_patch(agg_rect)
ax.text(6, 4.3, "Aggregation Engine: FedAvg | FedProx | SCAFFOLD", ha="center",
        fontsize=10, fontweight="bold", color="#1B5E20")

# Arrows between components
for i in range(len(components) - 1):
    ax.annotate("", xy=(components[i+1][0] - 0.9, components[i+1][1]),
                xytext=(components[i][0] + 0.9, components[i][1]),
                arrowprops=dict(arrowstyle="->", color="#666", lw=1.5))

# Client boxes
clients = [
    (2.0, 1.0, "Hospital A\n(n=200)"),
    (4.5, 1.0, "Hospital B\n(n=350)"),
    (7.0, 1.0, "Hospital C\n(n=150)"),
    (9.5, 1.0, "Hospital D\n(n=500)"),
]
for x_pos, y_pos, label in clients:
    rect = plt.Rectangle((x_pos - 0.8, y_pos - 0.5), 1.6, 1.0, fill=True,
                          facecolor="#fff3e0", edgecolor="#F57C00", linewidth=1.5, zorder=2)
    ax.add_patch(rect)
    ax.text(x_pos, y_pos, label, ha="center", va="center", fontsize=8)
    # Arrows
    ax.annotate("", xy=(x_pos, 3.8), xytext=(x_pos, y_pos + 0.5),
                arrowprops=dict(arrowstyle="<->", color="#F57C00", lw=1.5))

# PySyft layer
pysyft_rect = plt.Rectangle((1, 2.2), 10, 0.6, fill=True, facecolor="#f3e5f5",
                              edgecolor="#7B1FA2", linewidth=1.5, linestyle="--", zorder=1)
ax.add_patch(pysyft_rect)
ax.text(6, 2.5, "PySyft: Privacy Budget Enforcement | Secure Data Access | Audit Logging",
        ha="center", fontsize=9, color="#4A148C")

# Legend
legend_items = [
    ("#e8f4fd", "Server Layer"),
    ("#c8e6c9", "Aggregation"),
    ("#fff3e0", "Client (Hospital)"),
    ("#f3e5f5", "Privacy Layer"),
]
for i, (color, label) in enumerate(legend_items):
    rect = plt.Rectangle((0.5, 0.1 + i * 0.3), 0.3, 0.2, facecolor=color,
                          edgecolor="#666", linewidth=0.5)
    ax.add_patch(rect)
    ax.text(1.0, 0.2 + i * 0.3, label, fontsize=8, va="center")

plt.tight_layout()
plt.savefig("figures/fig6_architecture.png")
plt.savefig("figures/fig6_architecture.svg")
plt.close()
print("Saved fig6_architecture")

print("\nAll figures generated successfully.")
