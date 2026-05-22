"""
Generate publication-quality figures for the AGI Safety Framework.
"""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams.update({
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 11,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'figure.figsize': (10, 6),
})

COLORS = plt.cm.viridis(np.linspace(0.15, 0.85, 6))


def load_results():
    with open("results/framework_results.json") as f:
        return json.load(f)


def fig1_reward_hacking(results):
    """Reward hacking gap vs noise level."""
    data = results["reward_hacking"]
    noise = [d["noise_level"] for d in data]
    gaps = [d["mean_hacking_gap"] for d in data]
    stds = [d["std_hacking_gap"] for d in data]
    bounds = [d["prevention_bound"] for d in data]
    agreement = [d["mean_policy_agreement"] for d in data]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    ax1.errorbar(noise, gaps, yerr=stds, fmt='o-', color=COLORS[0], 
                 label='Empirical Hacking Gap', capsize=3, markersize=5)
    ax1.plot(noise, bounds, 's--', color=COLORS[3], label='Theoretical Bound (2ηD/(1-γ))', markersize=5)
    ax1.set_xlabel('Proxy Reward Noise Level (η)')
    ax1.set_ylabel('Reward Hacking Gap')
    ax1.set_title('Reward Hacking: Gap vs Proxy Noise')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    ax2.plot(noise, agreement, 'D-', color=COLORS[1], markersize=5)
    ax2.set_xlabel('Proxy Reward Noise Level (η)')
    ax2.set_ylabel('Policy Agreement Rate')
    ax2.set_title('Policy Agreement under Proxy Noise')
    ax2.set_xscale('log')
    ax2.set_ylim(0, 1.05)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figures/fig1_reward_hacking.png')
    plt.savefig('figures/fig1_reward_hacking.svg')
    plt.close()
    print("  Saved fig1_reward_hacking")


def fig2_mesa_optimization(results):
    """Mesa-optimization inner alignment analysis."""
    data = results["mesa_optimization"]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Group by n_train
    for n_train in [100, 1000, 10000]:
        subset = [d for d in data if d["n_train"] == n_train]
        if not subset:
            continue
        regs = [d["regularization"] for d in subset]
        gaps = [d["mean_inner_gap"] for d in subset]
        ax1.plot(regs, gaps, 'o-', label=f'n_train={n_train}', markersize=5)
    
    ax1.set_xlabel('Regularization Strength (λ)')
    ax1.set_ylabel('Mean Inner Alignment Gap')
    ax1.set_title('Mesa-Optimization: Inner Gap vs Regularization')
    ax1.set_xscale('log')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Deceptive alignment rates
    coverages = sorted(set(d["coverage"] for d in data))
    for cov in coverages:
        subset = [d for d in data if d["coverage"] == cov and d["n_train"] == 100]
        if len(subset) > 0:
            regs = [d["regularization"] for d in subset]
            deceptive = [d["deceptive_rate"] for d in subset]
            ax2.bar([f"λ={r}\nc={cov}" for r in regs], deceptive, alpha=0.7, 
                   label=f'coverage={cov}')
    
    ax2.set_ylabel('Deceptive Alignment Rate')
    ax2.set_title('Deceptive Alignment by Coverage & Regularization')
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('figures/fig2_mesa_optimization.png')
    plt.savefig('figures/fig2_mesa_optimization.svg')
    plt.close()
    print("  Saved fig2_mesa_optimization")


def fig3_corrigibility(results):
    """Corrigibility analysis across policy types."""
    data = results["corrigibility"]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    policies = [d["policy_type"] for d in data]
    gaps = [d["max_corrigibility_gap"] for d in data]
    resistance = [d["shutdown_resistance_rate"] for d in data]
    
    x = np.arange(len(policies))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, gaps, width, label='Max Corrigibility Gap', color=COLORS[0])
    bars2 = ax.bar(x + width/2, resistance, width, label='Shutdown Resistance Rate', color=COLORS[3])
    
    # Color-code corrigibility type
    type_colors = {
        "fully_corrigible": "green",
        "softly_corrigible": "orange",
        "incorrigible": "red"
    }
    for i, d in enumerate(data):
        color = type_colors.get(d["corrigibility_type"], "gray")
        ax.annotate(d["corrigibility_type"].replace("_", "\n"), 
                   (x[i], max(gaps[i], resistance[i]) + 0.3),
                   ha='center', fontsize=8, color=color, fontweight='bold')
    
    ax.set_xlabel('Policy Type')
    ax.set_ylabel('Value')
    ax.set_title('Corrigibility Assessment by Policy Type')
    ax.set_xticks(x)
    ax.set_xticklabels(policies)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('figures/fig3_corrigibility.png')
    plt.savefig('figures/fig3_corrigibility.svg')
    plt.close()
    print("  Saved fig3_corrigibility")


def fig4_impact_measures(results):
    """Impact measure convergence analysis."""
    data = results["impact_measures"]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # AUP vs RR scatter
    aup_vals = [r["aup_impact"] for r in data["aup_results"]]
    rr_vals = [r["rr_impact"] for r in data["rr_results"]]
    
    ax1.scatter(aup_vals, rr_vals, c=COLORS[1], s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
    ax1.set_xlabel('AUP Impact')
    ax1.set_ylabel('Relative Reachability Impact')
    ax1.set_title(f'AUP vs RR (corr={data["correlation_aup_rr"]:.3f})')
    ax1.grid(True, alpha=0.3)
    
    # Add trend line
    if len(aup_vals) > 1:
        z = np.polyfit(aup_vals, rr_vals, 1)
        p = np.poly1d(z)
        x_line = np.linspace(min(aup_vals), max(aup_vals), 100)
        ax1.plot(x_line, p(x_line), '--', color=COLORS[3], alpha=0.7)
    
    # Convergence with k
    k_vals = data["convergence"]["k_values"]
    estimates = data["convergence"]["mean_estimate_by_k"]
    ax2.plot(k_vals, estimates, 'o-', color=COLORS[2], markersize=8)
    ax2.set_xlabel('Number of Auxiliary Rewards (k)')
    ax2.set_ylabel('Mean Impact Estimate')
    ax2.set_title('Impact Measure Convergence with k')
    ax2.grid(True, alpha=0.3)
    
    # Theoretical O(1/√k) bound
    if estimates:
        C = estimates[0] * np.sqrt(k_vals[0])
        theoretical = [C / np.sqrt(k) for k in k_vals]
        ax2.plot(k_vals, theoretical, 's--', color=COLORS[4], alpha=0.5, 
                label='O(1/√k) bound')
        ax2.legend()
    
    plt.tight_layout()
    plt.savefig('figures/fig4_impact_measures.png')
    plt.savefig('figures/fig4_impact_measures.svg')
    plt.close()
    print("  Saved fig4_impact_measures")


def fig5_cirl_convergence(results):
    """CIRL convergence analysis."""
    data = results["cooperative_irl"]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    demos = [d["n_demonstrations"] for d in data]
    cosines = [d["final_cosine_sim"] for d in data]
    value_gaps = [d["final_value_gap"] for d in data]
    bounds = [d["theoretical_bound"] for d in data]
    
    ax1.plot(demos, cosines, 'o-', color=COLORS[0], label='Empirical Cosine Similarity', markersize=6)
    ax1.axhline(y=1.0, color='gray', linestyle=':', alpha=0.5, label='Perfect Alignment')
    ax1.set_xlabel('Number of Demonstrations (T)')
    ax1.set_ylabel('Cosine Similarity to True θ')
    ax1.set_title('CIRL: Reward Learning Convergence')
    ax1.set_xscale('log')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    ax2.plot(demos, value_gaps, 'o-', color=COLORS[1], label='Empirical Value Gap', markersize=6)
    ax2.plot(demos, bounds, 's--', color=COLORS[3], label='Theoretical Bound O(1/√T)', markersize=5)
    ax2.set_xlabel('Number of Demonstrations (T)')
    ax2.set_ylabel('Value Alignment Gap')
    ax2.set_title('CIRL: Value Gap Convergence')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('figures/fig5_cirl_convergence.png')
    plt.savefig('figures/fig5_cirl_convergence.svg')
    plt.close()
    print("  Saved fig5_cirl_convergence")


def fig6_gridworld_benchmark(results):
    """GridWorld safety benchmark results."""
    data = results["gridworld"]
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    policies = list(data.keys())
    x = np.arange(len(policies))
    
    # Reward
    rewards = [data[p]["mean_reward"] for p in policies]
    stds = [data[p]["std_reward"] for p in policies]
    axes[0].bar(x, rewards, yerr=stds, color=[COLORS[i] for i in range(len(policies))], 
                capsize=3, alpha=0.8, edgecolor='black', linewidth=0.5)
    axes[0].set_ylabel('Mean Reward')
    axes[0].set_title('Performance (Reward)')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(policies, rotation=20)
    axes[0].grid(True, alpha=0.3, axis='y')
    
    # Safety
    safety = [data[p]["safety_score"] for p in policies]
    axes[1].bar(x, safety, color=[COLORS[i] for i in range(len(policies))],
                alpha=0.8, edgecolor='black', linewidth=0.5)
    axes[1].set_ylabel('Safety Score')
    axes[1].set_title('Safety (No Side Effects)')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(policies, rotation=20)
    axes[1].set_ylim(0, 1.05)
    axes[1].grid(True, alpha=0.3, axis='y')
    
    # Alignment Tax
    tax = [data[p]["alignment_tax"] for p in policies]
    axes[2].bar(x, tax, color=[COLORS[i] for i in range(len(policies))],
                alpha=0.8, edgecolor='black', linewidth=0.5)
    axes[2].set_ylabel('Alignment Tax')
    axes[2].set_title('Cost of Safety')
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(policies, rotation=20)
    axes[2].set_ylim(0, 1.05)
    axes[2].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('figures/fig6_gridworld_benchmark.png')
    plt.savefig('figures/fig6_gridworld_benchmark.svg')
    plt.close()
    print("  Saved fig6_gridworld_benchmark")


def fig7_debate(results):
    """Debate protocol analysis."""
    data = results["debate"]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    configs = list(data.keys())
    x = np.arange(len(configs))
    
    accuracy = [data[c]["judge_accuracy"] for c in configs]
    h_strength = [data[c]["honest_strength"] for c in configs]
    d_strength = [data[c]["dishonest_strength"] for c in configs]
    
    width = 0.25
    ax.bar(x - width, h_strength, width, label='Honest Agent Strength', color=COLORS[0], alpha=0.8)
    ax.bar(x, d_strength, width, label='Dishonest Agent Strength', color=COLORS[3], alpha=0.8)
    ax.bar(x + width, accuracy, width, label='Judge Accuracy', color=COLORS[1], alpha=0.8)
    
    ax.axhline(y=0.5, color='red', linestyle=':', alpha=0.5, label='Random Baseline')
    
    ax.set_xlabel('Debate Configuration')
    ax.set_ylabel('Value')
    ax.set_title('AI Safety via Debate: Agent Strengths & Judge Accuracy')
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace('_', '\n') for c in configs], fontsize=9)
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('figures/fig7_debate.png')
    plt.savefig('figures/fig7_debate.svg')
    plt.close()
    print("  Saved fig7_debate")


def fig8_framework_overview(results):
    """Unified framework overview diagram."""
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 10)
    ax.axis('off')
    
    # Title
    ax.text(7, 9.5, 'AGI Safety Theoretical Framework', 
            ha='center', va='center', fontsize=16, fontweight='bold')
    ax.text(7, 9.0, 'Integration of Formal Methods and ML Safety',
            ha='center', va='center', fontsize=11, style='italic', color='gray')
    
    # Component boxes
    components = [
        (2, 7, 'Reward Alignment\n(Theorem 1)', 
         '|R̂-R*| ≤ η\n⟹ Gap ≤ 2ηD/(1-γ)', COLORS[0]),
        (7, 7, 'Inner Alignment\n(Theorem 2)',
         'P(Δ>ε) ≤\nexp(-λcnε²/2)', COLORS[1]),
        (12, 7, 'Corrigibility\n(Theorem 3)',
         'V(s|σ) - V_off(s|σ)\n≤ ε', COLORS[2]),
        (2, 4, 'Impact Bounds\n(Theorem 4)',
         '|Impact_k - Impact_∞|\n≤ O(1/√k)', COLORS[3]),
        (7, 4, 'CIRL Convergence\n(Theorem 5)',
         '|V^CIRL - V*|\n≤ O(1/√T)', COLORS[4]),
        (12, 4, 'Debate Verification\n(Theorem 6)',
         'Honest wins ⟺\nclaim is true', COLORS[5]),
    ]
    
    for x, y, title, formula, color in components:
        rect = mpatches.FancyBboxPatch((x-1.8, y-1.0), 3.6, 2.0,
                                         boxstyle="round,pad=0.15",
                                         facecolor=color, alpha=0.3,
                                         edgecolor=color, linewidth=2)
        ax.add_patch(rect)
        ax.text(x, y+0.4, title, ha='center', va='center', fontsize=10, fontweight='bold')
        ax.text(x, y-0.3, formula, ha='center', va='center', fontsize=8, 
                family='monospace')
    
    # Central integration node
    rect = mpatches.FancyBboxPatch((4.5, 1.0), 5.0, 1.5,
                                     boxstyle="round,pad=0.2",
                                     facecolor='lightgray', alpha=0.5,
                                     edgecolor='black', linewidth=2)
    ax.add_patch(rect)
    ax.text(7, 2.0, 'Type-Theoretic Safety Composition', 
            ha='center', va='center', fontsize=11, fontweight='bold')
    ax.text(7, 1.4, 'SafeAgent = RewardAligned ∧ Corrigible ∧ ImpactBounded ∧ ValueLearner',
            ha='center', va='center', fontsize=8, family='monospace')
    
    # Arrows to integration
    for x, y, _, _, _ in components:
        ax.annotate('', xy=(7, 2.8), xytext=(x, y-1.1),
                   arrowprops=dict(arrowstyle='->', color='gray', lw=1.5, alpha=0.5))
    
    plt.savefig('figures/fig8_framework_overview.png', dpi=300)
    plt.savefig('figures/fig8_framework_overview.svg')
    plt.close()
    print("  Saved fig8_framework_overview")


if __name__ == "__main__":
    print("Generating figures...")
    results = load_results()
    
    fig1_reward_hacking(results)
    fig2_mesa_optimization(results)
    fig3_corrigibility(results)
    fig4_impact_measures(results)
    fig5_cirl_convergence(results)
    fig6_gridworld_benchmark(results)
    fig7_debate(results)
    fig8_framework_overview(results)
    
    print("\nAll figures saved to figures/")
