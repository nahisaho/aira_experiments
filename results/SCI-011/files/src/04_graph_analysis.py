"""
04_graph_analysis.py
--------------------
全脳コネクトーム解析 - ステップ4: グラフ理論解析
スモールワールド性 / モジュール性 / ハブ構造 (NetworkX + 解析的計算)
"""

import numpy as np
import networkx as nx
import json
import os
from datetime import datetime
from itertools import combinations

RESULTS_DIR = "results"
LOGS_DIR = "logs"
DATA_DIR = "data"
np.random.seed(42)


# ─────────────────────────────────────────────────────────────────────────────
# 1.  コネクトームグラフの構築
# ─────────────────────────────────────────────────────────────────────────────
def fc_to_graph(
    fc_matrix: np.ndarray,
    threshold: float | None = None,
    density: float = 0.15,
    binarize: bool = False,
) -> nx.Graph:
    """
    FC 行列からグラフを構築。
    density 指定時: 上位 density% のエッジを保持 (プロポーショナル閾値)
    """
    n = fc_matrix.shape[0]
    fc = fc_matrix.copy()
    np.fill_diagonal(fc, 0)
    fc = np.abs(fc)  # 負の相関も接続として扱う

    if threshold is None:
        # プロポーショナル閾値
        n_edges_target = int(n * (n - 1) / 2 * density)
        upper = fc[np.triu_indices(n, k=1)]
        threshold = np.sort(upper)[::-1][n_edges_target]

    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if fc[i, j] >= threshold:
                weight = fc[i, j] if not binarize else 1.0
                G.add_edge(i, j, weight=weight)

    return G


# ─────────────────────────────────────────────────────────────────────────────
# 2.  スモールワールド指標
# ─────────────────────────────────────────────────────────────────────────────
def compute_small_world(
    G: nx.Graph, n_rand: int = 50, rng_seed: int = 42
) -> dict:
    """
    スモールワールド指標 σ = (C/C_rand) / (L/L_rand) を計算。
    σ > 1 → スモールワールドネットワーク

    C : クラスタリング係数 (ローカル効率の代理)
    L : 特性経路長
    """
    rng = np.random.default_rng(rng_seed)

    # 実際のグラフ指標
    if not nx.is_connected(G):
        largest_cc = max(nx.connected_components(G), key=len)
        G_conn = G.subgraph(largest_cc).copy()
    else:
        G_conn = G

    C_real = nx.average_clustering(G_conn)
    # 特性経路長 (重み逆数)
    try:
        L_real = nx.average_shortest_path_length(
            G_conn,
            weight=lambda u, v, d: 1.0 / d.get("weight", 1.0) if d.get("weight", 0) > 0 else 1e6,
        )
    except Exception:
        L_real = nx.average_shortest_path_length(G_conn)

    # ランダムグラフ (同密度・同次数列) との比較
    n_nodes = G_conn.number_of_nodes()
    n_edges = G_conn.number_of_edges()
    C_rands, L_rands = [], []

    for _ in range(n_rand):
        G_rand = nx.gnm_random_graph(n_nodes, n_edges, seed=int(rng.integers(1e6)))
        if not nx.is_connected(G_rand):
            largest_cc = max(nx.connected_components(G_rand), key=len)
            G_rand = G_rand.subgraph(largest_cc).copy()
        C_rands.append(nx.average_clustering(G_rand))
        try:
            L_rands.append(nx.average_shortest_path_length(G_rand))
        except Exception:
            L_rands.append(np.nan)

    C_rand = np.nanmean(C_rands)
    L_rand = np.nanmean(L_rands)

    gamma = C_real / (C_rand + 1e-9)  # 正規化クラスタリング係数
    lam = L_real / (L_rand + 1e-9)    # 正規化経路長
    sigma = gamma / (lam + 1e-9)      # スモールワールド指標

    return {
        "clustering_coeff": round(C_real, 4),
        "characteristic_path_length": round(L_real, 4),
        "clustering_coeff_rand": round(C_rand, 4),
        "path_length_rand": round(L_rand, 4),
        "gamma": round(gamma, 3),
        "lambda": round(lam, 3),
        "sigma": round(sigma, 3),
        "is_small_world": bool(sigma > 1),
    }


# ─────────────────────────────────────────────────────────────────────────────
# 3.  モジュール性解析
# ─────────────────────────────────────────────────────────────────────────────
def compute_modularity(G: nx.Graph) -> dict:
    """
    Louvain 法 (greedy modularity) によるコミュニティ検出。
    modularity Q: 0.3–0.7 → 中等度のモジュール構造
    """
    from networkx.algorithms.community import greedy_modularity_communities

    communities = list(greedy_modularity_communities(G, weight="weight"))
    Q = nx.algorithms.community.modularity(G, communities, weight="weight")

    # 各コミュニティのサイズ
    community_sizes = sorted([len(c) for c in communities], reverse=True)

    # 参加係数 (Participation Coefficient): 各ノードのコミュニティ間接続比率
    node_to_comm = {}
    for k, comm in enumerate(communities):
        for node in comm:
            node_to_comm[node] = k

    participation = {}
    for node in G.nodes():
        degree = G.degree(node, weight="weight")
        if degree == 0:
            participation[node] = 0.0
            continue
        comm_weights = {}
        for _, v, d in G.edges(node, data=True):
            c = node_to_comm.get(v, -1)
            comm_weights[c] = comm_weights.get(c, 0) + d.get("weight", 1.0)
        pc = 1 - sum((w / degree) ** 2 for w in comm_weights.values())
        participation[node] = round(float(pc), 4)

    return {
        "modularity_Q": round(float(Q), 4),
        "n_communities": len(communities),
        "community_sizes": community_sizes,
        "mean_participation_coeff": round(float(np.mean(list(participation.values()))), 4),
        "participation_coefficients": participation,
    }


# ─────────────────────────────────────────────────────────────────────────────
# 4.  ハブ構造解析
# ─────────────────────────────────────────────────────────────────────────────
def compute_hub_structure(G: nx.Graph, n_hubs: int = 12) -> dict:
    """
    ハブノードを次数中心性・媒介中心性・固有ベクトル中心性の複合スコアで同定。
    """
    n = G.number_of_nodes()

    # 各中心性指標
    degree_cent = nx.degree_centrality(G)
    betw_cent = nx.betweenness_centrality(G, weight="weight", normalized=True)
    try:
        eigen_cent = nx.eigenvector_centrality(G, weight="weight", max_iter=1000)
    except nx.PowerIterationFailedConvergence:
        eigen_cent = {node: 0.0 for node in G.nodes()}

    # 局所効率
    local_eff = nx.local_efficiency(G)

    # 複合ハブスコア (z スコア平均)
    def z_norm(d: dict) -> dict:
        vals = np.array(list(d.values()))
        mu, sigma = vals.mean(), vals.std() + 1e-9
        return {k: (v - mu) / sigma for k, v in d.items()}

    dc_z = z_norm(degree_cent)
    bc_z = z_norm(betw_cent)
    ec_z = z_norm(eigen_cent)

    hub_score = {
        node: (dc_z[node] + bc_z[node] + ec_z[node]) / 3
        for node in G.nodes()
    }

    # 上位ハブノード
    hub_nodes = sorted(hub_score, key=hub_score.get, reverse=True)[:n_hubs]

    # AAL90 ラベル（簡略化）
    aal90_labels = {
        0: "PreCG_L", 1: "PreCG_R", 2: "SFGdor_L", 3: "SFGdor_R",
        7: "MFG_L",  8: "MFG_R",  25: "PCC_L",   26: "PCC_R",
        30: "IPL_L", 31: "IPL_R", 67: "MTG_L",   68: "MTG_R",
    }

    hub_info = [
        {
            "node_id": int(h),
            "label": aal90_labels.get(h, f"ROI_{h}"),
            "hub_score": round(float(hub_score[h]), 4),
            "degree_centrality": round(float(degree_cent[h]), 4),
            "betweenness_centrality": round(float(betw_cent[h]), 4),
            "eigenvector_centrality": round(float(eigen_cent[h]), 4),
        }
        for h in hub_nodes
    ]

    return {
        "n_hub_nodes": n_hubs,
        "hub_nodes": hub_info,
        "global_efficiency": round(nx.global_efficiency(G), 4),
        "local_efficiency": round(local_eff, 4),
        "mean_degree": round(float(np.mean([d for _, d in G.degree()])), 2),
    }


# ─────────────────────────────────────────────────────────────────────────────
# 5.  マルチスケール解析（密度閾値スイープ）
# ─────────────────────────────────────────────────────────────────────────────
def density_threshold_sweep(
    fc_matrix: np.ndarray, densities: list[float] = None
) -> list[dict]:
    """
    複数の接続密度でグラフ指標を計算し，閾値依存性を評価。
    """
    if densities is None:
        densities = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]

    results = []
    for d in densities:
        G = fc_to_graph(fc_matrix, density=d)
        if not nx.is_connected(G):
            largest_cc = max(nx.connected_components(G), key=len)
            G = G.subgraph(largest_cc).copy()
        C = nx.average_clustering(G)
        try:
            L = nx.average_shortest_path_length(G)
        except Exception:
            L = np.nan
        results.append({
            "density": d,
            "n_edges": G.number_of_edges(),
            "clustering_coeff": round(C, 4),
            "path_length": round(float(L), 4) if not np.isnan(L) else None,
            "global_efficiency": round(nx.global_efficiency(G), 4),
            "local_efficiency": round(nx.local_efficiency(G), 4),
        })
    return results


# ─────────────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("[04] グラフ理論解析中...")
    ts = datetime.utcnow().isoformat()

    groups = ["HC", "SCZ", "AD"]
    graph_metrics = {}

    for group in groups:
        fc_path = f"{DATA_DIR}/FC_pearson_{group}.npy"
        if not os.path.exists(fc_path):
            print(f"  [{group}] FC データが見つかりません。スキップ。")
            continue

        fc_all = np.load(fc_path)
        fc_avg = fc_all.mean(axis=0)  # 群平均 FC

        print(f"  [{group}] グラフ構築中 (density=0.15) ...")
        G = fc_to_graph(fc_avg, density=0.15)

        sw = compute_small_world(G, n_rand=30)
        mod = compute_modularity(G)
        hub = compute_hub_structure(G)
        sweep = density_threshold_sweep(fc_avg)

        graph_metrics[group] = {
            "small_world": sw,
            "modularity": {k: v for k, v in mod.items() if k != "participation_coefficients"},
            "hub_structure": {k: v for k, v in hub.items() if k != "hub_nodes"},
            "top_hubs": hub["hub_nodes"][:5],
            "density_sweep": sweep,
        }

        print(f"    → σ={sw['sigma']:.3f} (SW={sw['is_small_world']}), "
              f"Q={mod['modularity_Q']:.3f}, "
              f"hubs_top={hub['hub_nodes'][0]['label']}")

    with open(f"{RESULTS_DIR}/graph_metrics.json", "w") as f:
        json.dump(graph_metrics, f, indent=2, ensure_ascii=False, default=str)

    log_entry = {
        "timestamp": ts,
        "phase": "graph_analysis",
        "event_type": "step_completed",
        "actor": "co-scientist",
        "skill_or_tool": "04_graph_analysis.py",
        "handoff_out": {
            g: {
                "sigma": graph_metrics[g]["small_world"]["sigma"],
                "Q": graph_metrics[g]["modularity"]["modularity_Q"],
            }
            for g in graph_metrics
        },
        "files_written": [f"{RESULTS_DIR}/graph_metrics.json"],
        "status": "ok",
    }
    with open(f"{LOGS_DIR}/process-log.jsonl", "a") as f:
        f.write(json.dumps(log_entry, ensure_ascii=False) + "\n")

    return graph_metrics


if __name__ == "__main__":
    main()
