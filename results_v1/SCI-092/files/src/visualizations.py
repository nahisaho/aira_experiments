"""
Visualization module: Generate all figures for the social acceptance prediction model.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import json, os

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300, "font.size": 10,
    "axes.titlesize": 12, "axes.labelsize": 10,
    "figure.figsize": (10, 7)
})
COLORS = sns.color_palette("viridis", 5)
TECH_COLORS = {"gene_editing": "#440154", "AI": "#21918c", "nuclear_fusion": "#fde725"}


def fig1_meta_analysis_forest(meta_df, meta_results):
    """Forest plot of pooled acceptance rates by technology."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 6), sharey=False)
    techs = ["gene_editing", "AI", "nuclear_fusion"]
    titles = ["Gene Editing", "Artificial Intelligence", "Nuclear Fusion"]

    for ax, tech, title in zip(axes, techs, titles):
        sub = meta_df[meta_df["technology"] == tech].sort_values("acceptance_rate")
        y_pos = range(len(sub))
        ax.errorbar(sub["acceptance_rate"], y_pos, xerr=1.96*sub["se"],
                     fmt="o", color=TECH_COLORS[tech], markersize=4, capsize=2, linewidth=0.8)
        ax.set_yticks(list(y_pos))
        ax.set_yticklabels(sub["study_id"], fontsize=6)

        pooled = meta_results[tech]
        ax.axvline(pooled["pooled_effect"], color="red", linestyle="--", linewidth=1.5, alpha=0.8)
        ax.axvspan(pooled["ci_95"][0], pooled["ci_95"][1], alpha=0.15, color="red")
        ax.set_xlabel("Acceptance Rate")
        ax.set_title(f"{title}\nPooled={pooled['pooled_effect']:.3f} "
                     f"[{pooled['ci_95'][0]:.3f}, {pooled['ci_95'][1]:.3f}]\n"
                     f"I²={pooled['I2']:.1f}%")

    plt.tight_layout()
    plt.savefig("figures/fig1_meta_analysis_forest.png", bbox_inches="tight")
    plt.savefig("figures/fig1_meta_analysis_forest.svg", bbox_inches="tight")
    plt.close()


def fig2_sentiment_trends(sent_df, sent_results):
    """Temporal sentiment trends by technology."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel A: Hybrid score distribution by technology
    ax = axes[0]
    techs = ["gene_editing", "AI", "nuclear_fusion"]
    labels = ["Gene Editing", "AI", "Nuclear Fusion"]
    for tech, label in zip(techs, labels):
        sub = sent_df[sent_df["technology"] == tech]
        ax.hist(sub["hybrid_score"], bins=30, alpha=0.5, label=label, color=TECH_COLORS[tech], density=True)
    ax.set_xlabel("Hybrid Sentiment Score")
    ax.set_ylabel("Density")
    ax.set_title("(A) Sentiment Score Distribution")
    ax.legend()
    ax.axvline(0, color="gray", linestyle="--", alpha=0.5)

    # Panel B: Temporal trend
    ax = axes[1]
    for tech, label in zip(techs, labels):
        sub = sent_df[sent_df["technology"] == tech]
        yearly = sub.groupby("year")["hybrid_score"].agg(["mean", "sem"]).reset_index()
        ax.errorbar(yearly["year"], yearly["mean"], yerr=1.96*yearly["sem"],
                     marker="o", label=label, color=TECH_COLORS[tech], capsize=3)
    ax.set_xlabel("Year")
    ax.set_ylabel("Mean Hybrid Sentiment Score")
    ax.set_title("(B) Sentiment Temporal Trend")
    ax.legend()
    ax.axhline(0, color="gray", linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig("figures/fig2_sentiment_trends.png", bbox_inches="tight")
    plt.savefig("figures/fig2_sentiment_trends.svg", bbox_inches="tight")
    plt.close()


def fig3_psychometric_space(psych_df, psych_results):
    """Dread-Unknown risk space plot."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: Scatter in factor space
    ax = axes[0]
    techs = ["gene_editing", "AI", "nuclear_fusion"]
    labels = ["Gene Editing", "AI", "Nuclear Fusion"]
    for tech, label in zip(techs, labels):
        sub = psych_df[psych_df["technology"] == tech]
        ax.scatter(sub["dread_factor"], sub["unknown_factor"], alpha=0.15,
                   color=TECH_COLORS[tech], s=10)
        mean_d = sub["dread_factor"].mean()
        mean_u = sub["unknown_factor"].mean()
        ax.scatter(mean_d, mean_u, color=TECH_COLORS[tech], s=200, marker="*",
                   edgecolors="black", linewidth=1, zorder=5, label=label)
    ax.set_xlabel("Dread Factor")
    ax.set_ylabel("Unknown Factor")
    ax.set_title("(A) Psychometric Risk Space")
    ax.legend()
    ax.axhline(0, color="gray", linestyle="--", alpha=0.3)
    ax.axvline(0, color="gray", linestyle="--", alpha=0.3)

    # Panel B: Risk vs Acceptance
    ax = axes[1]
    for tech, label in zip(techs, labels):
        sub = psych_df[psych_df["technology"] == tech]
        ax.scatter(sub["overall_risk_perception"], sub["acceptance"], alpha=0.15,
                   color=TECH_COLORS[tech], s=10, label=label)
        # Regression line
        z = np.polyfit(sub["overall_risk_perception"], sub["acceptance"], 1)
        p = np.poly1d(z)
        x_line = np.linspace(1, 7, 50)
        ax.plot(x_line, p(x_line), color=TECH_COLORS[tech], linewidth=2, alpha=0.8)
    ax.set_xlabel("Overall Risk Perception")
    ax.set_ylabel("Acceptance")
    ax.set_title("(B) Risk Perception vs Acceptance")
    ax.legend()

    plt.tight_layout()
    plt.savefig("figures/fig3_psychometric_space.png", bbox_inches="tight")
    plt.savefig("figures/fig3_psychometric_space.svg", bbox_inches="tight")
    plt.close()


def fig4_framing_effects(framing_df, framing_results):
    """Framing condition comparison."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    techs = ["gene_editing", "AI", "nuclear_fusion"]
    titles = ["Gene Editing", "AI", "Nuclear Fusion"]
    frame_order = ["benefit_frame", "scientific_frame", "balanced_frame",
                   "ethical_frame", "risk_frame"]
    frame_labels = ["Benefit", "Scientific", "Balanced", "Ethical", "Risk"]

    for ax, tech, title in zip(axes, techs, titles):
        sub = framing_df[framing_df["technology"] == tech]
        means = [sub[sub["framing_condition"] == f]["acceptance"].mean() for f in frame_order]
        sems = [sub[sub["framing_condition"] == f]["acceptance"].sem() for f in frame_order]

        bars = ax.barh(range(len(frame_order)), means, xerr=[1.96*s for s in sems],
                       color=[plt.cm.viridis(i/4) for i in range(5)], capsize=3)
        ax.set_yticks(range(len(frame_order)))
        ax.set_yticklabels(frame_labels)
        ax.set_xlabel("Mean Acceptance Score")
        r = framing_results[tech]["anova"]
        ax.set_title(f"{title}\nF={r['F_statistic']:.1f}, η²={r['eta_squared']:.3f}")
        ax.axvline(sub["acceptance"].mean(), color="gray", linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig("figures/fig4_framing_effects.png", bbox_inches="tight")
    plt.savefig("figures/fig4_framing_effects.svg", bbox_inches="tight")
    plt.close()


def fig5_sem_path_diagram(sem_results):
    """SEM path diagram visualization."""
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.set_xlim(-0.5, 10.5)
    ax.set_ylim(-0.5, 8.5)
    ax.set_aspect("equal")
    ax.axis("off")

    # Node positions
    nodes = {
        "Literacy": (1, 7), "Media": (1, 1),
        "Trust": (4, 5.5), "Risk": (7, 7),
        "Benefit": (7, 1), "Acceptance": (10, 4)
    }

    # Draw nodes
    for name, (x, y) in nodes.items():
        circle = plt.Circle((x, y), 0.7, fill=True, facecolor="#e8e8e8",
                             edgecolor="black", linewidth=2)
        ax.add_patch(circle)
        ax.text(x, y, name, ha="center", va="center", fontsize=9, fontweight="bold")

    # Draw paths
    paths = sem_results["structural_paths"]
    for p in paths:
        if p["from"] in nodes and p["to"] in nodes:
            x1, y1 = nodes[p["from"]]
            x2, y2 = nodes[p["to"]]
            color = "#2ca02c" if p["estimate"] > 0 else "#d62728"
            width = abs(p["estimate"]) * 4
            sig = "***" if p["p_value"] < 0.001 else "**" if p["p_value"] < 0.01 else "*" if p["p_value"] < 0.05 else ""
            dx, dy = x2 - x1, y2 - y1
            dist = np.sqrt(dx**2 + dy**2)
            offset = 0.75
            ax.annotate("", xy=(x2 - dx/dist*offset, y2 - dy/dist*offset),
                        xytext=(x1 + dx/dist*offset, y1 + dy/dist*offset),
                        arrowprops=dict(arrowstyle="->", color=color, lw=max(width, 1)))
            mid_x, mid_y = (x1+x2)/2, (y1+y2)/2
            perp_x, perp_y = -dy/dist*0.3, dx/dist*0.3
            ax.text(mid_x + perp_x, mid_y + perp_y,
                    f"β={p['estimate']:.2f}{sig}",
                    fontsize=8, ha="center", color=color, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))

    # Fit indices
    fit = sem_results["fit_indices"]
    fit_text = f"Model Fit: CFI={fit.get('CFI', 'N/A')}, RMSEA={fit.get('RMSEA', 'N/A')}, GFI={fit.get('GFI', 'N/A')}"
    ax.text(5.5, -0.2, fit_text, ha="center", fontsize=10, style="italic")

    ax.set_title("Trust-Acceptance Structural Equation Model", fontsize=14, fontweight="bold", pad=20)
    plt.savefig("figures/fig5_sem_path_diagram.png", bbox_inches="tight")
    plt.savefig("figures/fig5_sem_path_diagram.svg", bbox_inches="tight")
    plt.close()


def fig6_japan_case_study(japan_df, japan_results):
    """Japan genome-edited food case study multi-panel figure."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))

    # Panel A: Acceptance distribution
    ax = axes[0, 0]
    ax.hist(japan_df["acceptance"], bins=30, color="#21918c", edgecolor="black", alpha=0.7)
    ax.axvline(japan_df["acceptance"].mean(), color="red", linestyle="--", linewidth=2,
               label=f'Mean={japan_df["acceptance"].mean():.2f}')
    ax.axvline(4, color="orange", linestyle=":", linewidth=2, label="Neutral (4.0)")
    ax.set_xlabel("Acceptance Score")
    ax.set_ylabel("Frequency")
    ax.set_title("(A) Acceptance Distribution")
    ax.legend()

    # Panel B: Regression coefficients
    ax = axes[0, 1]
    reg = japan_results["regression"]
    coefs = sorted(reg["coefficients"].items(), key=lambda x: x[1])
    names = [c[0].replace("_", " ").title() for c in coefs]
    values = [c[1] for c in coefs]
    colors = ["#d62728" if v < 0 else "#2ca02c" for v in values]
    ax.barh(range(len(names)), values, color=colors)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("Regression Coefficient (β)")
    ax.set_title(f"(B) Predictors of Acceptance (R²={reg['r_squared']:.3f})")
    ax.axvline(0, color="black", linewidth=0.5)

    # Panel C: Trust comparison
    ax = axes[1, 0]
    trust_vars = ["trust_government", "trust_scientists", "trust_corporations"]
    trust_labels = ["Government", "Scientists", "Corporations"]
    trust_means = [japan_df[v].mean() for v in trust_vars]
    trust_sems = [japan_df[v].sem() for v in trust_vars]
    bars = ax.bar(trust_labels, trust_means, yerr=[1.96*s for s in trust_sems],
                  color=["#440154", "#21918c", "#fde725"], capsize=5, edgecolor="black")
    ax.set_ylabel("Trust Score (1-7)")
    ax.set_title("(C) Institutional Trust Levels")
    ax.set_ylim(1, 7)
    ax.axhline(4, color="gray", linestyle="--", alpha=0.5, label="Neutral")

    # Panel D: Acceptance by age group
    ax = axes[1, 1]
    age_order = ["18-29", "30-44", "45-59", "60+"]
    age_data = [japan_df[japan_df["age_group"] == g]["acceptance"].values for g in age_order]
    bp = ax.boxplot(age_data, labels=age_order, patch_artist=True)
    for patch, color in zip(bp["boxes"], [plt.cm.viridis(i/3) for i in range(4)]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_xlabel("Age Group")
    ax.set_ylabel("Acceptance Score")
    ax.set_title("(D) Acceptance by Age Group")

    plt.suptitle("Genome-Edited Food Acceptance in Japan (n=600)", fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig("figures/fig6_japan_case_study.png", bbox_inches="tight")
    plt.savefig("figures/fig6_japan_case_study.svg", bbox_inches="tight")
    plt.close()


def fig7_integrated_heatmap(japan_df):
    """Correlation heatmap of key variables."""
    fig, ax = plt.subplots(figsize=(10, 8))
    key_vars = ["scientific_literacy", "media_exposure", "trust_composite",
                "knowledge_genome_editing", "perceived_risk", "perceived_benefit",
                "naturalness_concern", "acceptance", "purchase_intention"]
    labels = ["Sci. Literacy", "Media Exp.", "Trust", "Knowledge",
              "Risk Perc.", "Benefit Perc.", "Naturalness", "Acceptance", "Purchase Int."]
    corr = japan_df[key_vars].corr()
    mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
    sns.heatmap(corr, mask=mask, annot=True, fmt=".2f", cmap="RdBu_r",
                center=0, vmin=-1, vmax=1, square=True,
                xticklabels=labels, yticklabels=labels, ax=ax)
    ax.set_title("Correlation Matrix: Key Variables (Japan GE Food Study)", fontsize=12)
    plt.tight_layout()
    plt.savefig("figures/fig7_correlation_heatmap.png", bbox_inches="tight")
    plt.savefig("figures/fig7_correlation_heatmap.svg", bbox_inches="tight")
    plt.close()


def generate_all_figures(meta_df, meta_results, sent_df, sent_results,
                         psych_df, psych_results, framing_df, framing_results,
                         sem_results, japan_df, japan_results):
    """Generate all figures."""
    os.makedirs("figures", exist_ok=True)
    fig1_meta_analysis_forest(meta_df, meta_results)
    fig2_sentiment_trends(sent_df, sent_results)
    fig3_psychometric_space(psych_df, psych_results)
    fig4_framing_effects(framing_df, framing_results)
    fig5_sem_path_diagram(sem_results)
    fig6_japan_case_study(japan_df, japan_results)
    fig7_integrated_heatmap(japan_df)
    print("All 7 figures generated in figures/")
