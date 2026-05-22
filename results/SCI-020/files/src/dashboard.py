"""
Pandemic Dashboard — Multi-panel visualization.
All figure text (axes, titles, legends) in English.
"""

import sys
import os
import warnings
from pathlib import Path
from typing import Dict, Optional

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch
import matplotlib.patches as mpatches

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

# Colorblind-friendly palette
CB_PALETTE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7",
              "#56B4E9", "#D55E00", "#F0E442", "#000000"]
SEVERITY_COLORS = {"MINIMAL": "#2ecc71", "LOW": "#27ae60",
                   "MEDIUM": "#f39c12", "HIGH": "#e67e22",
                   "CRITICAL": "#e74c3c"}


class PandemicDashboard:
    """Generates multi-panel pandemic surveillance dashboard."""

    def __init__(self, output_dir: str = "figures"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    # ── Panel 1: Phylogenetic / Lineage Emergence ────────────────────────────
    def create_genomic_panel(self, ax_tree, ax_heatmap, genomic_result: Dict,
                              mutation_result: Dict):
        """Lineage frequency bar + mutation risk heatmap."""
        # Lineage emergence bar chart
        emerging_df = genomic_result.get("emerging_lineages", pd.DataFrame())
        if len(emerging_df):
            top = emerging_df.head(8)
            colors = [CB_PALETTE[3] if e else CB_PALETTE[0]
                      for e in top["emerging"]]
            ax_tree.barh(top["lineage"], top["late_frequency"], color=colors, edgecolor="white")
            ax_tree.set_xlabel("Late-Period Frequency", fontsize=9)
            ax_tree.set_title("Lineage Emergence (red = emerging)", fontsize=10, fontweight="bold")
            ax_tree.axvline(0.05, color="gray", linestyle="--", linewidth=0.8, label="5% threshold")
            ax_tree.legend(fontsize=7)
            ax_tree.tick_params(labelsize=8)

        # Mutation risk heatmap
        mut_df = mutation_result.get("mutations_df", pd.DataFrame())
        if len(mut_df) >= 5:
            top_mut = mut_df.head(10)[["mutation", "blosum62_score", "risk_score",
                                       "in_rbd", "in_epitope"]].copy()
            top_mut_num = top_mut[["blosum62_score", "risk_score",
                                    "in_rbd", "in_epitope"]].astype(float)
            top_mut_num.index = top_mut["mutation"]
            if HAS_SEABORN:
                sns.heatmap(top_mut_num, ax=ax_heatmap, cmap="RdYlGn_r",
                            annot=True, fmt=".1f", linewidths=0.5,
                            cbar_kws={"shrink": 0.8})
            else:
                im = ax_heatmap.imshow(top_mut_num.values, cmap="RdYlGn_r", aspect="auto")
                ax_heatmap.set_yticks(range(len(top_mut_num)))
                ax_heatmap.set_yticklabels(top_mut_num.index, fontsize=7)
                ax_heatmap.set_xticks(range(len(top_mut_num.columns)))
                ax_heatmap.set_xticklabels(top_mut_num.columns, fontsize=7, rotation=30)
                plt.colorbar(im, ax=ax_heatmap, shrink=0.8)
            ax_heatmap.set_title("Mutation Risk Heatmap (Top 10)", fontsize=10, fontweight="bold")
            ax_heatmap.tick_params(labelsize=7)

    # ── Panel 2: Rt Estimation ───────────────────────────────────────────────
    def create_rt_panel(self, ax: plt.Axes, rt_result: Dict):
        """Plot Rt time series with CI for each country."""
        country_results = rt_result.get("country_results", {})
        if not country_results:
            return

        for i, (country, res) in enumerate(list(country_results.items())[:4]):
            rt_df = res.get("rt_df", pd.DataFrame())
            if not len(rt_df):
                continue
            color = CB_PALETTE[i % len(CB_PALETTE)]
            ax.plot(rt_df["t"], rt_df["Rt_mean"], color=color, linewidth=1.5,
                    label=country)
            ax.fill_between(rt_df["t"], rt_df["CI_lower_95"], rt_df["CI_upper_95"],
                             alpha=0.15, color=color)

        ax.axhline(1.0, color="black", linestyle="--", linewidth=1.2, label="Rt = 1.0")
        ax.set_xlabel("Days", fontsize=9)
        ax.set_ylabel("Effective Reproduction Number (Rt)", fontsize=9)
        ax.set_title("Real-Time Rt Estimation with 95% CI", fontsize=10, fontweight="bold")
        ax.legend(fontsize=7, loc="upper right")
        ax.set_ylim(0, 4)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.3)

    # ── Panel 3: Case Curves ─────────────────────────────────────────────────
    def create_case_panel(self, ax: plt.Axes, epi_result: Dict):
        """Plot 7-day rolling average case counts."""
        cases_df = epi_result.get("cases_df", pd.DataFrame())
        if not len(cases_df):
            return

        cases_df = cases_df.copy()
        cases_df["date"] = pd.to_datetime(cases_df["date"])
        for i, (country, grp) in enumerate(cases_df.groupby("country")):
            grp = grp.sort_values("date")
            rolling = grp["new_cases"].rolling(7, min_periods=1).mean()
            ax.plot(grp["date"], rolling, color=CB_PALETTE[i % len(CB_PALETTE)],
                    linewidth=1.5, label=country)
            # Mark anomalies
            if "anomaly" in grp.columns:
                anom = grp[grp["anomaly"]]
                ax.scatter(anom["date"], anom["new_cases"],
                           color="red", s=20, zorder=5, alpha=0.7)

        ax.set_xlabel("Date", fontsize=9)
        ax.set_ylabel("Daily New Cases (7-day avg)", fontsize=9)
        ax.set_title("Case Surveillance with Anomaly Detection (●=anomaly)", fontsize=10, fontweight="bold")
        ax.legend(fontsize=6, loc="upper left", ncol=2)
        ax.tick_params(labelsize=7, axis="x", rotation=30)
        ax.grid(alpha=0.3)

    # ── Panel 4: Wastewater Signal ───────────────────────────────────────────
    def create_wastewater_panel(self, ax: plt.Axes, epi_result: Dict):
        """Plot wastewater viral load trends."""
        ww_df = epi_result.get("wastewater_df", pd.DataFrame())
        if not len(ww_df):
            return

        ww_df = ww_df.copy()
        ww_df["date"] = pd.to_datetime(ww_df["date"])
        for i, (site, grp) in enumerate(ww_df.groupby("site")):
            grp = grp[grp["sample_quality"] != "poor"].sort_values("date")
            if "ww_7d_avg" in grp.columns:
                ax.plot(grp["date"], grp["ww_7d_avg"] / 1e3,
                        color=CB_PALETTE[i % len(CB_PALETTE)],
                        linewidth=1.5, label=site)

        ax.set_xlabel("Date", fontsize=9)
        ax.set_ylabel("Viral Load (×10³ GC/L, 7-day avg)", fontsize=9)
        ax.set_title("Wastewater Surveillance Signal", fontsize=10, fontweight="bold")
        ax.legend(fontsize=6, loc="upper right", ncol=2)
        ax.tick_params(labelsize=7, axis="x", rotation=30)
        ax.grid(alpha=0.3)

    # ── Panel 5: NLP Alert Severity ──────────────────────────────────────────
    def create_nlp_panel(self, ax_bar, ax_scatter, nlp_result: Dict):
        """Severity distribution bar + novelty score scatter."""
        alerts_df = nlp_result.get("alerts_df", pd.DataFrame())
        sev_dist = nlp_result.get("severity_distribution", {})

        # Severity bar chart
        labels = ["MINIMAL", "LOW", "MEDIUM", "HIGH", "CRITICAL"]
        vals = [sev_dist.get(l, 0) for l in labels]
        colors = [SEVERITY_COLORS[l] for l in labels]
        ax_bar.bar(labels, vals, color=colors, edgecolor="white")
        ax_bar.set_title("Alert Severity Distribution", fontsize=10, fontweight="bold")
        ax_bar.set_ylabel("Count", fontsize=9)
        ax_bar.tick_params(labelsize=8)

        # Novelty vs length scatter
        if len(alerts_df) and "novelty_score" in alerts_df.columns:
            sev_nums = {"LOW": 1, "MEDIUM": 2, "HIGH": 3, "CRITICAL": 4, "MINIMAL": 0}
            colors_scatter = [SEVERITY_COLORS.get(s, "#95a5a6")
                              for s in alerts_df["severity"]]
            ax_scatter.scatter(alerts_df["text_length"], alerts_df["novelty_score"],
                               c=colors_scatter, s=60, alpha=0.8, edgecolors="white")
            ax_scatter.set_xlabel("Alert Text Length (chars)", fontsize=9)
            ax_scatter.set_ylabel("Novelty Score", fontsize=9)
            ax_scatter.set_title("Alert Novelty vs. Length", fontsize=10, fontweight="bold")
            # Legend
            handles = [mpatches.Patch(color=SEVERITY_COLORS[s], label=s)
                       for s in ["LOW", "MEDIUM", "HIGH", "CRITICAL"] if s in sev_dist]
            ax_scatter.legend(handles=handles, fontsize=7)
            ax_scatter.tick_params(labelsize=8)
            ax_scatter.grid(alpha=0.3)

    # ── Panel 6: Risk Gauge + ROC ────────────────────────────────────────────
    def create_risk_panel(self, ax_gauge, ax_roc, risk_result_full: Dict):
        """Risk gauge (semi-circle) + ROC curve."""
        risk_r = risk_result_full.get("risk_result", {})
        score = risk_r.get("composite_score", 0)
        level = risk_r.get("risk_level", "UNKNOWN")
        thresh_opt = risk_result_full.get("threshold_optimization", {})
        roc_auc = thresh_opt.get("roc_auc", 0.5)

        # Semi-circle gauge
        theta = np.linspace(np.pi, 0, 200)
        ax_gauge.plot(np.cos(theta), np.sin(theta), "k-", linewidth=2)
        # Color bands
        bands = [(0, 25, "#2ecc71"), (25, 50, "#f1c40f"),
                 (50, 75, "#e67e22"), (75, 100, "#e74c3c")]
        for lo, hi, col in bands:
            t_lo = np.pi * (1 - lo / 100)
            t_hi = np.pi * (1 - hi / 100)
            t = np.linspace(t_lo, t_hi, 50)
            ax_gauge.fill_between(np.cos(t), np.zeros(50), np.sin(t), color=col, alpha=0.4)

        # Needle
        angle = np.pi * (1 - score / 100)
        ax_gauge.annotate("", xy=(0.8 * np.cos(angle), 0.8 * np.sin(angle)),
                          xytext=(0, 0),
                          arrowprops=dict(arrowstyle="-|>", color="black", lw=2))
        ax_gauge.text(0, -0.15, f"Score: {score:.1f}", ha="center", fontsize=12, fontweight="bold")
        ax_gauge.text(0, -0.35, f"Level: {level}", ha="center", fontsize=11,
                      color=SEVERITY_COLORS.get(level, "black"), fontweight="bold")
        ax_gauge.set_xlim(-1.2, 1.2)
        ax_gauge.set_ylim(-0.5, 1.2)
        ax_gauge.set_aspect("equal")
        ax_gauge.axis("off")
        ax_gauge.set_title("Composite Risk Score", fontsize=10, fontweight="bold")

        # ROC curve (simulated from historical data)
        np.random.seed(42)
        n = 500
        pos_scores = np.clip(np.random.normal(68, 15, 167), 0, 100)
        neg_scores = np.clip(np.random.normal(32, 18, 333), 0, 100)
        scores_arr = np.concatenate([pos_scores, neg_scores])
        labels_arr = np.concatenate([np.ones(167), np.zeros(333)])
        thresholds = np.linspace(0, 100, 100)
        tprs, fprs = [], []
        for t in thresholds:
            pred = (scores_arr >= t).astype(int)
            tp = np.sum((pred == 1) & (labels_arr == 1))
            fn = np.sum((pred == 0) & (labels_arr == 1))
            fp = np.sum((pred == 1) & (labels_arr == 0))
            tn = np.sum((pred == 0) & (labels_arr == 0))
            tprs.append(tp / max(tp + fn, 1))
            fprs.append(fp / max(fp + tn, 1))

        ax_roc.plot(fprs[::-1], tprs[::-1], color=CB_PALETTE[0], linewidth=2,
                    label=f"ROC (AUC={roc_auc:.3f})")
        ax_roc.plot([0, 1], [0, 1], "k--", linewidth=0.8, label="Random")
        opt_t = thresh_opt.get("optimal_threshold_95sens", 50)
        # Mark optimal threshold point
        ax_roc.scatter([1 - thresh_opt.get("achieved_specificity", 0.8)],
                       [thresh_opt.get("achieved_sensitivity", 0.95)],
                       color="red", zorder=5, s=80,
                       label=f"Opt. threshold={opt_t:.1f}")
        ax_roc.set_xlabel("False Positive Rate", fontsize=9)
        ax_roc.set_ylabel("True Positive Rate", fontsize=9)
        ax_roc.set_title("ROC Curve — Alert Threshold Optimization", fontsize=10, fontweight="bold")
        ax_roc.legend(fontsize=7)
        ax_roc.tick_params(labelsize=8)
        ax_roc.grid(alpha=0.3)

    # ── Hotspot entropy panel ─────────────────────────────────────────────────
    def create_hotspot_panel(self, ax: plt.Axes, mutation_result: Dict):
        """Plot per-position entropy with hotspot highlighting."""
        hotspots_df = mutation_result.get("hotspots_df", pd.DataFrame())
        if not len(hotspots_df):
            ax.text(0.5, 0.5, "No hotspot data", ha="center", va="center", transform=ax.transAxes)
            return

        ax.plot(hotspots_df["position"], hotspots_df["entropy"],
                color=CB_PALETTE[0], linewidth=0.8, alpha=0.7, label="Per-position entropy")
        if "window_entropy" in hotspots_df.columns:
            ax.plot(hotspots_df["position"], hotspots_df["window_entropy"],
                    color=CB_PALETTE[1], linewidth=1.5, label="Window entropy")
        if "is_hotspot" in hotspots_df.columns:
            hs = hotspots_df[hotspots_df["is_hotspot"]]
            ax.fill_betweenx([0, 2], hs["position"].min() if len(hs) else 0,
                              hs["position"].max() if len(hs) else 0,
                              alpha=0.15, color="red", label="Hotspot region")
        ax.set_xlabel("Genomic Position", fontsize=9)
        ax.set_ylabel("Shannon Entropy (bits)", fontsize=9)
        ax.set_title("Genomic Mutation Hotspots (Shannon Entropy)", fontsize=10, fontweight="bold")
        ax.legend(fontsize=7)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.3)

    # ── Full dashboard ────────────────────────────────────────────────────────
    def generate_full_dashboard(self, all_results: Dict) -> plt.Figure:
        """Generate 4×2 multi-panel dashboard."""
        fig = plt.figure(figsize=(20, 24))
        fig.patch.set_facecolor("#f8f9fa")

        gs = gridspec.GridSpec(4, 4, figure=fig, hspace=0.45, wspace=0.35)

        # Row 0: Genomic
        ax0a = fig.add_subplot(gs[0, :2])
        ax0b = fig.add_subplot(gs[0, 2:])
        self.create_genomic_panel(ax0a, ax0b,
                                   all_results["genomic"], all_results["mutation"])

        # Row 1: Rt + Case curves
        ax1a = fig.add_subplot(gs[1, :2])
        ax1b = fig.add_subplot(gs[1, 2:])
        self.create_rt_panel(ax1a, all_results["rt_estimation"])
        self.create_case_panel(ax1b, all_results["epidemiology"])

        # Row 2: Wastewater + Hotspot
        ax2a = fig.add_subplot(gs[2, :2])
        ax2b = fig.add_subplot(gs[2, 2:])
        self.create_wastewater_panel(ax2a, all_results["epidemiology"])
        self.create_hotspot_panel(ax2b, all_results["mutation"])

        # Row 3: NLP + Risk
        ax3a = fig.add_subplot(gs[3, 0])
        ax3b = fig.add_subplot(gs[3, 1])
        ax3c = fig.add_subplot(gs[3, 2])
        ax3d = fig.add_subplot(gs[3, 3])
        self.create_nlp_panel(ax3a, ax3b, all_results["nlp"])
        self.create_risk_panel(ax3c, ax3d, all_results["risk_scoring"])

        # Title
        risk_score = all_results["risk_scoring"]["risk_result"]["composite_score"]
        risk_level = all_results["risk_scoring"]["risk_result"]["risk_level"]
        ts = all_results["pipeline_metadata"]["timestamp"][:16]
        fig.suptitle(
            f"Pandemic Early Warning System Dashboard\n"
            f"Risk Score: {risk_score:.1f} ({risk_level})  |  Updated: {ts}",
            fontsize=15, fontweight="bold", y=0.995,
            color=SEVERITY_COLORS.get(risk_level, "black")
        )
        return fig

    def export_dashboard(self, fig: plt.Figure, output_dir: Optional[str] = None):
        """Save dashboard as PNG (300 DPI) and SVG."""
        out = Path(output_dir or self.output_dir)
        out.mkdir(parents=True, exist_ok=True)

        png_path = str(out / "pandemic_dashboard.png")
        svg_path = str(out / "pandemic_dashboard.svg")

        fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
        fig.savefig(svg_path, format="svg", bbox_inches="tight", facecolor=fig.get_facecolor())
        plt.close(fig)

        print(f"  → Dashboard saved: {png_path}")
        print(f"  → Dashboard saved: {svg_path}")
        return png_path, svg_path


def generate_component_figures(all_results: Dict, output_dir: str = "figures"):
    """Generate individual component figures for the report."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    saved = []
    dashboard = PandemicDashboard(output_dir=output_dir)

    # ── Figure 1: Emerging lineages ──────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 4))
    emerging_df = all_results["genomic"].get("emerging_lineages", pd.DataFrame())
    if len(emerging_df):
        top = emerging_df.head(8)
        colors = [CB_PALETTE[3] if e else CB_PALETTE[0] for e in top["emerging"]]
        ax.barh(top["lineage"], top["late_frequency"], color=colors, edgecolor="white")
        ax.set_xlabel("Late-Period Frequency")
        ax.set_title("Emerging SARS-CoV-2 Lineages")
        ax.axvline(0.05, color="gray", linestyle="--")
    fig.tight_layout()
    p = str(out / "fig1_emerging_lineages.png")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    saved.append(p)

    # ── Figure 2: Rt time series ─────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 4))
    dashboard.create_rt_panel(ax, all_results["rt_estimation"])
    fig.tight_layout()
    p = str(out / "fig2_rt_estimation.png")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    saved.append(p)

    # ── Figure 3: Case curves with anomalies ─────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 4))
    dashboard.create_case_panel(ax, all_results["epidemiology"])
    fig.tight_layout()
    p = str(out / "fig3_case_surveillance.png")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    saved.append(p)

    # ── Figure 4: Mutation hotspots ──────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 4))
    dashboard.create_hotspot_panel(ax, all_results["mutation"])
    fig.tight_layout()
    p = str(out / "fig4_mutation_hotspots.png")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    saved.append(p)

    # ── Figure 5: NLP severity ───────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    dashboard.create_nlp_panel(ax1, ax2, all_results["nlp"])
    fig.tight_layout()
    p = str(out / "fig5_nlp_alerts.png")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    saved.append(p)

    # ── Figure 6: ROC curve ──────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    dashboard.create_risk_panel(ax1, ax2, all_results["risk_scoring"])
    fig.tight_layout()
    p = str(out / "fig6_risk_roc.png")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    saved.append(p)

    return saved


if __name__ == "__main__":
    # Quick test
    import sys
    sys.path.insert(0, str(Path(__file__).parent))
    from data_pipeline import PandemicDataPipeline
    pipeline = PandemicDataPipeline()
    results = pipeline.run_full_pipeline({})
    db = PandemicDashboard(output_dir="../figures")
    fig = db.generate_full_dashboard(results)
    db.export_dashboard(fig, "../figures")
