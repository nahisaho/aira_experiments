"""
07_visualization.py
-------------------
全脳コネクトーム解析 - 可視化
FC行列ヒートマップ / グラフ指標比較 / バイオマーカーAUC / 信頼性分布
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import json
import os

FIGURES_DIR = "figures"
RESULTS_DIR = "results"
DATA_DIR = "data"
os.makedirs(FIGURES_DIR, exist_ok=True)

NETWORK_LABELS = {
    "DMN": list(range(0, 15)),
    "FPN": list(range(15, 28)),
    "SMN": list(range(28, 42)),
    "VIS": list(range(42, 55)),
    "DAN": list(range(55, 65)),
    "SN":  list(range(65, 75)),
    "LIM": list(range(75, 90)),
}
NET_NAMES = list(NETWORK_LABELS.keys())
NET_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2"]

GROUP_COLORS = {"HC": "#4CAF50", "SCZ": "#F44336", "AD": "#FF9800"}


# ─────────────────────────────────────────────────────────────────────────────
def plot_fc_matrices():
    """Figure 1: 3群の平均FC行列ヒートマップ"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    groups = ["HC", "SCZ", "AD"]
    group_titles = ["Healthy Control (HC)", "Schizophrenia (SCZ)", "Alzheimer's Disease (AD)"]

    for ax, group, title in zip(axes, groups, group_titles):
        fp = f"{DATA_DIR}/FC_pearson_{group}.npy"
        if not os.path.exists(fp):
            ax.text(0.5, 0.5, "No data", ha="center", va="center")
            ax.set_title(title)
            continue
        fc = np.load(fp).mean(axis=0)
        im = ax.imshow(fc, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xlabel("ROI Index (AAL90)")
        ax.set_ylabel("ROI Index (AAL90)")
        # ネットワーク境界線
        boundaries = [0]
        for net in NET_NAMES:
            boundaries.append(boundaries[-1] + len(NETWORK_LABELS[net]))
        for b in boundaries[1:-1]:
            ax.axhline(b - 0.5, color="white", lw=0.7, alpha=0.8)
            ax.axvline(b - 0.5, color="white", lw=0.7, alpha=0.8)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Pearson r")

    plt.suptitle("Mean Functional Connectivity Matrices (AAL90)", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{FIGURES_DIR}/fig1_fc_matrices.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{FIGURES_DIR}/fig1_fc_matrices.pdf", bbox_inches="tight")
    plt.close()
    print("  → Figure 1 saved: FC matrices")


def plot_graph_metrics():
    """Figure 2: グラフ理論指標の群間比較"""
    fp = f"{RESULTS_DIR}/graph_metrics.json"
    if not os.path.exists(fp):
        print("  [skip] graph_metrics.json not found")
        return

    with open(fp) as f:
        gm = json.load(f)

    groups = [g for g in ["HC", "SCZ", "AD"] if g in gm]
    metrics_keys = [
        ("small_world", "sigma", "Small-world index (σ)"),
        ("small_world", "clustering_coeff", "Clustering Coefficient (C)"),
        ("small_world", "characteristic_path_length", "Characteristic Path Length (L)"),
        ("modularity", "modularity_Q", "Modularity (Q)"),
        ("hub_structure", "global_efficiency", "Global Efficiency (E_glob)"),
        ("hub_structure", "local_efficiency", "Local Efficiency (E_loc)"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.flatten()

    for ax, (outer_key, inner_key, label) in zip(axes, metrics_keys):
        vals = [gm[g][outer_key][inner_key] for g in groups if outer_key in gm[g]]
        bars = ax.bar(groups, vals, color=[GROUP_COLORS[g] for g in groups],
                      edgecolor="black", linewidth=0.8, width=0.5)
        ax.set_title(label, fontsize=10, fontweight="bold")
        ax.set_ylabel("Value")
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001,
                    f"{val:.3f}", ha="center", va="bottom", fontsize=8)
        # HC基準の有意差指標
        if len(vals) > 1:
            max_v = max(vals)
            for i, g in enumerate(groups[1:], 1):
                diff_pct = (vals[i] - vals[0]) / (abs(vals[0]) + 1e-9) * 100
                ax.text(i, vals[i] + max_v * 0.05,
                        f"Δ{diff_pct:+.1f}%", ha="center", fontsize=7, color="gray")

    plt.suptitle("Graph-Theoretic Metrics: HC vs SCZ vs AD", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{FIGURES_DIR}/fig2_graph_metrics.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{FIGURES_DIR}/fig2_graph_metrics.pdf", bbox_inches="tight")
    plt.close()
    print("  → Figure 2 saved: Graph metrics comparison")


def plot_biomarker_auc():
    """Figure 3: バイオマーカー分類器の AUC 比較"""
    fp = f"{RESULTS_DIR}/biomarker_results.json"
    if not os.path.exists(fp):
        print("  [skip] biomarker_results.json not found")
        return

    with open(fp) as f:
        br = json.load(f)

    diseases = list(br.keys())
    feat_sets = ["all_edges", "network_fc", "graph_metrics"]
    feat_labels = ["All Edges (4005D)", "Network FC (28D)", "Graph Metrics (180D)"]
    clf_names = ["SVM_RBF", "RandomForest", "LogisticRegression_L1", "GradientBoosting"]
    clf_markers = ["o", "s", "^", "D"]
    clf_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

    fig, axes = plt.subplots(1, len(diseases), figsize=(14, 6))
    if len(diseases) == 1:
        axes = [axes]

    for ax, dis in zip(axes, diseases):
        x = np.arange(len(feat_sets))
        for ci, (clf, marker, color) in enumerate(zip(clf_names, clf_markers, clf_colors)):
            aucs = []
            stds = []
            for feat in feat_sets:
                clf_res = br[dis]["classification"][feat]["classifiers"].get(clf, {})
                aucs.append(clf_res.get("auc_mean", 0.5))
                stds.append(clf_res.get("auc_std", 0))
            offset = (ci - 1.5) * 0.12
            ax.errorbar(x + offset, aucs, yerr=stds, fmt=marker, color=color,
                        label=clf, markersize=8, capsize=4, linewidth=1.5)

        ax.axhline(0.5, color="gray", linestyle="--", alpha=0.7, label="Chance (0.50)")
        ax.axhline(0.8, color="orange", linestyle=":", alpha=0.5, label="Good (0.80)")
        ax.set_xticks(x)
        ax.set_xticklabels(feat_labels, rotation=20, ha="right", fontsize=9)
        ax.set_ylim(0.4, 1.05)
        ax.set_ylabel("AUC (5-fold CV)")
        ax.set_title(f"HC vs {dis} Classification", fontsize=11, fontweight="bold")
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.suptitle("Biomarker Classification Performance", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{FIGURES_DIR}/fig3_biomarker_auc.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{FIGURES_DIR}/fig3_biomarker_auc.pdf", bbox_inches="tight")
    plt.close()
    print("  → Figure 3 saved: Biomarker AUC")


def plot_network_importance():
    """Figure 4: ネットワーク重要度（効果量ヒートマップ）"""
    fp = f"{RESULTS_DIR}/biomarker_results.json"
    if not os.path.exists(fp):
        return

    with open(fp) as f:
        br = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    diseases = ["SCZ", "AD"]

    for ax, dis in zip(axes, diseases):
        if dis not in br:
            continue
        n_nets = len(NET_NAMES)
        d_mat = np.zeros((n_nets, n_nets))
        for item in br[dis]["network_importance_top5"][:25]:
            ni_idx = NET_NAMES.index(item["network_i"]) if item["network_i"] in NET_NAMES else -1
            nj_idx = NET_NAMES.index(item["network_j"]) if item["network_j"] in NET_NAMES else -1
            if ni_idx >= 0 and nj_idx >= 0:
                d_mat[ni_idx, nj_idx] = item["cohen_d"]
                d_mat[nj_idx, ni_idx] = item["cohen_d"]

        # 全ネットワークペアの Cohen's d を fc_group_comparison から取得してマップ
        # ここではネットワーク重要度の top5 の情報のみ使用
        fc_path_hc = f"{DATA_DIR}/FC_pearson_HC.npy"
        fc_path_dis = f"{DATA_DIR}/FC_pearson_{dis}.npy"
        if os.path.exists(fc_path_hc) and os.path.exists(fc_path_dis):
            fc_hc = np.load(fc_path_hc)
            fc_dis_arr = np.load(fc_path_dis)
            for i, ni in enumerate(NET_NAMES):
                for j, nj in enumerate(NET_NAMES):
                    ri = np.array(NETWORK_LABELS[ni])
                    rj = np.array(NETWORK_LABELS[nj])
                    hc_v = fc_hc[:, ri[:, None], rj[None, :]].mean(axis=(1, 2))
                    dis_v = fc_dis_arr[:, ri[:, None], rj[None, :]].mean(axis=(1, 2))
                    d_mat[i, j] = abs(
                        (hc_v.mean() - dis_v.mean())
                        / np.sqrt((hc_v.std()**2 + dis_v.std()**2) / 2 + 1e-9)
                    )

        im = ax.imshow(d_mat, cmap="YlOrRd", vmin=0, vmax=1.5)
        ax.set_xticks(range(n_nets))
        ax.set_yticks(range(n_nets))
        ax.set_xticklabels(NET_NAMES, rotation=45, ha="right")
        ax.set_yticklabels(NET_NAMES)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="|Cohen's d|")
        ax.set_title(f"HC vs {dis}\nNetwork-level Effect Size", fontsize=11, fontweight="bold")
        # 値を表示
        for i in range(n_nets):
            for j in range(n_nets):
                ax.text(j, i, f"{d_mat[i,j]:.2f}", ha="center", va="center",
                        fontsize=7, color="black" if d_mat[i,j] < 1.0 else "white")

    plt.suptitle("Network-level FC Differences: Cohen's d", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{FIGURES_DIR}/fig4_network_importance.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{FIGURES_DIR}/fig4_network_importance.pdf", bbox_inches="tight")
    plt.close()
    print("  → Figure 4 saved: Network importance")


def plot_reliability():
    """Figure 5: ICC 分布 + Fingerprinting + 再現性戦略比較"""
    fp = f"{RESULTS_DIR}/reliability_results.json"
    if not os.path.exists(fp):
        return

    with open(fp) as f:
        rr = json.load(f)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # (a) ICC 分布
    ax = axes[0]
    cats = ["Poor\n(<0.5)", "Moderate\n(0.5-0.75)", "Good\n(0.75-0.9)", "Excellent\n(>0.9)"]
    icc2 = rr["icc2_1"]["icc_distribution"]
    counts = [icc2["poor_lt05"], icc2["moderate_05_075"],
              icc2["good_075_09"], icc2["excellent_gt09"]]
    colors_icc = ["#e74c3c", "#f39c12", "#27ae60", "#2980b9"]
    bars = ax.bar(cats, counts, color=colors_icc, edgecolor="black")
    ax.set_ylabel("Number of FC Edges")
    ax.set_title(f"ICC(2,1) Distribution\nMean ICC = {rr['icc2_1']['mean_icc']:.3f}",
                 fontsize=11, fontweight="bold")
    for bar, cnt in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                str(cnt), ha="center", va="bottom", fontsize=9)

    # (b) Fingerprinting
    ax = axes[1]
    fp_data = rr["fingerprinting"]
    net_accs = fp_data["network_fingerprint_accuracy"]
    nets = list(net_accs.keys())
    accs = list(net_accs.values())
    chance = fp_data["chance_level_pct"]
    overall = fp_data["identification_rate_pct"]
    bar_colors = [NET_COLORS[i % len(NET_COLORS)] for i in range(len(nets))]
    ax.barh(nets, accs, color=bar_colors, edgecolor="black")
    ax.axvline(chance, color="red", linestyle="--", label=f"Chance ({chance:.1f}%)")
    ax.axvline(overall, color="black", linestyle="-",
               linewidth=2, label=f"Overall ({overall:.1f}%)")
    ax.set_xlabel("Identification Accuracy (%)")
    ax.set_title("FC Fingerprinting by Network", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 100)

    # (c) 再現性戦略比較
    ax = axes[2]
    strategies = rr["reproducibility_strategies"]
    names = list(strategies.keys())
    names_short = ["Baseline\n(no GSR)", "With GSR", "Fisher-z", "Normalized"]
    icc_vals = [strategies[n]["mean_icc"] for n in names]
    pct_vals = [strategies[n]["pct_above_0_75"] for n in names]
    x = np.arange(len(names))
    ax2 = ax.twinx()
    bars1 = ax.bar(x - 0.2, icc_vals, width=0.35, color="#3498db",
                   alpha=0.8, label="Mean ICC", edgecolor="black")
    bars2 = ax2.bar(x + 0.2, pct_vals, width=0.35, color="#e67e22",
                    alpha=0.8, label="% ICC>0.75", edgecolor="black")
    ax.set_xticks(x)
    ax.set_xticklabels(names_short[:len(names)], fontsize=9)
    ax.set_ylabel("Mean ICC(2,1)", color="#3498db")
    ax2.set_ylabel("% Edges with ICC > 0.75", color="#e67e22")
    ax.set_title("Reproducibility Strategy Comparison", fontsize=11, fontweight="bold")
    ax.legend(loc="upper left", fontsize=8)
    ax2.legend(loc="upper right", fontsize=8)
    ax.set_ylim(0, 1.0)
    ax2.set_ylim(0, 100)

    plt.suptitle("Test-Retest Reliability Assessment", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{FIGURES_DIR}/fig5_reliability.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{FIGURES_DIR}/fig5_reliability.pdf", bbox_inches="tight")
    plt.close()
    print("  → Figure 5 saved: Reliability")


def plot_density_sweep():
    """Figure 6: 接続密度スイープ（閾値依存性）"""
    fp = f"{RESULTS_DIR}/graph_metrics.json"
    if not os.path.exists(fp):
        return

    with open(fp) as f:
        gm = json.load(f)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    groups = [g for g in ["HC", "SCZ", "AD"] if g in gm]

    metric_pairs = [
        ("clustering_coeff", "Clustering Coefficient (C)"),
        ("global_efficiency", "Global Efficiency (E_glob)"),
    ]

    for ax, (metric, ylabel) in zip(axes, metric_pairs):
        for group in groups:
            sweep = gm[group].get("density_sweep", [])
            if not sweep:
                continue
            densities = [s["density"] for s in sweep]
            vals = [s.get(metric) for s in sweep]
            vals = [v if v is not None else np.nan for v in vals]
            ax.plot(densities, vals, marker="o", label=group,
                    color=GROUP_COLORS[group], linewidth=2, markersize=6)
        ax.set_xlabel("Connection Density")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel, fontsize=11, fontweight="bold")
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.axvline(0.15, color="gray", linestyle="--", alpha=0.5, label="Default (0.15)")

    plt.suptitle("Graph Metric Stability Across Connection Densities", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{FIGURES_DIR}/fig6_density_sweep.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{FIGURES_DIR}/fig6_density_sweep.pdf", bbox_inches="tight")
    plt.close()
    print("  → Figure 6 saved: Density sweep")


def plot_pipeline_overview():
    """Figure 7: パイプライン概要図"""
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 7)
    ax.axis("off")

    steps = [
        (1,   5.5, "#3498DB", "Step 1\nPreprocessing\n(FSL/FNIRT/topup)"),
        (3.5, 5.5, "#9B59B6", "Step 2\nStructural\nConnectivity\n(iFOD2/SIFT2)"),
        (6.0, 5.5, "#27AE60", "Step 3\nFunctional\nConnectivity\n(Pearson/Partial/DFC)"),
        (8.5, 5.5, "#E67E22", "Step 4\nGraph Theory\n(σ/Q/Hubs)"),
        (11,  5.5, "#E74C3C", "Step 5\nBiomarkers\n(SVM/RF/LR)"),
        (7.0, 2.5, "#1ABC9C", "Step 6\nTest-Retest\nReliability\n(ICC/Fingerprint)"),
    ]

    # ノードを描画
    for x, y, color, label in steps:
        circle = plt.Circle((x, y), 0.9, color=color, zorder=3, alpha=0.9)
        ax.add_patch(circle)
        ax.text(x, y, label, ha="center", va="center", fontsize=7.5,
                fontweight="bold", color="white", zorder=4)

    # 矢印 (step 1→2→3→4→5)
    arrow_chain = [(1, 3.5), (3.5, 6.0), (6.0, 8.5), (8.5, 11.0)]
    for x1, x2 in arrow_chain:
        ax.annotate("", xy=(x2 - 0.95, 5.5), xytext=(x1 + 0.95, 5.5),
                    arrowprops=dict(arrowstyle="->", color="#333", lw=1.5))

    # Step 3→6 矢印
    ax.annotate("", xy=(6.8, 3.3), xytext=(6.2, 4.6),
                arrowprops=dict(arrowstyle="->", color="#333", lw=1.5))

    # 入力データ
    ax.text(1, 7.2, "fMRI (BOLD)", ha="center", fontsize=9, style="italic",
            bbox=dict(facecolor="#EBF5FB", edgecolor="#3498DB", boxstyle="round"))
    ax.text(3.5, 7.2, "dMRI (DWI)", ha="center", fontsize=9, style="italic",
            bbox=dict(facecolor="#F5EEF8", edgecolor="#9B59B6", boxstyle="round"))
    ax.annotate("", xy=(1, 6.5), xytext=(1, 7.0),
                arrowprops=dict(arrowstyle="->", color="#3498DB", lw=1.2))
    ax.annotate("", xy=(3.5, 6.5), xytext=(3.5, 7.0),
                arrowprops=dict(arrowstyle="->", color="#9B59B6", lw=1.2))

    # 出力
    ax.text(11, 7.2, "Report + Results", ha="center", fontsize=9, style="italic",
            bbox=dict(facecolor="#FDEDEC", edgecolor="#E74C3C", boxstyle="round"))
    ax.annotate("", xy=(11, 6.5), xytext=(11, 7.0),
                arrowprops=dict(arrowstyle="<-", color="#E74C3C", lw=1.2))

    ax.set_title("Whole-Brain Connectome Analysis Pipeline (FSL/FreeSurfer/NetworkX)",
                 fontsize=13, fontweight="bold", pad=10)
    plt.tight_layout()
    plt.savefig(f"{FIGURES_DIR}/fig7_pipeline_overview.png", dpi=150, bbox_inches="tight")
    plt.savefig(f"{FIGURES_DIR}/fig7_pipeline_overview.pdf", bbox_inches="tight")
    plt.close()
    print("  → Figure 7 saved: Pipeline overview")


def main():
    print("[07] 可視化中...")
    plot_pipeline_overview()
    plot_fc_matrices()
    plot_graph_metrics()
    plot_biomarker_auc()
    plot_network_importance()
    plot_reliability()
    plot_density_sweep()
    print("  → 全図表を figures/ に保存完了")


if __name__ == "__main__":
    main()
