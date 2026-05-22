"""
03_functional_connectivity.py
-----------------------------
全脳コネクトーム解析 - ステップ3: 機能的コネクティビティ
相関 / 偏相関 / 動的FC (スライディングウィンドウ + HMM)
"""

import numpy as np
import json
import os
from datetime import datetime
from scipy import stats, linalg
from sklearn.covariance import LedoitWolf

RESULTS_DIR = "results"
LOGS_DIR = "logs"
DATA_DIR = "data"
np.random.seed(42)


# ─────────────────────────────────────────────────────────────────────────────
# 1.  時系列シミュレーション (マルチROI BOLD信号)
# ─────────────────────────────────────────────────────────────────────────────
def simulate_bold_timeseries(
    n_rois: int = 90,
    n_timepoints: int = 300,
    tr: float = 2.0,
    n_subjects: int = 30,
    group: str = "HC",
) -> np.ndarray:
    """
    ネットワーク構造を持つ模擬 BOLD 時系列を生成。
    デフォルトモード (DMN)、実行系 (FPN)、感覚運動系 (SMN) を模倣。
    """
    rng = np.random.default_rng({"HC": 42, "SCZ": 99, "AD": 55}[group])

    # ネットワーク定義 (ROI インデックス)
    networks = {
        "DMN": list(range(0, 15)),       # 前頭内側、後部帯状回等
        "FPN": list(range(15, 28)),      # 前頭前野、頭頂葉
        "SMN": list(range(28, 42)),      # 運動野、補足運動野
        "VIS": list(range(42, 55)),      # 視覚野
        "DAN": list(range(55, 65)),      # 背側注意網
        "SN":  list(range(65, 75)),      # 顕著性ネットワーク
        "LIM": list(range(75, 90)),      # 辺縁系
    }

    # 疾患ネットワーク改変係数
    network_noise = {
        "HC": {"DMN": 1.0, "FPN": 1.0, "SN": 1.0},
        "SCZ": {"DMN": 0.65, "FPN": 0.70, "SN": 1.40},   # DMN低下, SN亢進
        "AD":  {"DMN": 0.45, "FPN": 0.60, "SN": 0.90},   # 全般的低下
    }[group]

    all_ts = []
    for _ in range(n_subjects):
        ts = rng.standard_normal((n_rois, n_timepoints)) * 0.3  # ベースラインノイズ

        # 各ネットワーク内で共通ドライバーを加算
        for net_name, roi_list in networks.items():
            strength = network_noise.get(net_name, 1.0)
            driver = rng.standard_normal(n_timepoints)
            # バンドパスフィルタ相当 (0.01-0.1 Hz)
            t = np.arange(n_timepoints) * tr
            driver += 0.5 * np.sin(2 * np.pi * 0.05 * t)
            for roi in roi_list:
                ts[roi] += driver * strength * rng.uniform(0.7, 1.3)

        # 全体的なスケーリング (BOLD 単位: % signal change)
        ts = ts / ts.std(axis=1, keepdims=True)
        all_ts.append(ts)

    return np.array(all_ts)  # (n_subjects, n_rois, n_timepoints)


# ─────────────────────────────────────────────────────────────────────────────
# 2.  静的 FC 計算
# ─────────────────────────────────────────────────────────────────────────────
def compute_static_fc(timeseries: np.ndarray) -> dict[str, np.ndarray]:
    """
    全時系列にわたる静的機能的コネクティビティを計算。

    Returns:
      pearson_fc : ピアソン相関行列
      partial_fc : 正規化偏相関行列 (Ledoit-Wolf 収縮推定)
    """
    n_subjects, n_rois, n_tp = timeseries.shape
    pearson_all = np.zeros((n_subjects, n_rois, n_rois))
    partial_all = np.zeros((n_subjects, n_rois, n_rois))

    lw = LedoitWolf()

    for s in range(n_subjects):
        ts = timeseries[s]  # (n_rois, n_timepoints)

        # ピアソン相関
        corr = np.corrcoef(ts)
        np.fill_diagonal(corr, 0)
        pearson_all[s] = corr

        # 偏相関 (Ledoit-Wolf 精度行列)
        lw.fit(ts.T)  # fit expects (samples, features)
        prec = lw.precision_
        # 正規化偏相関に変換
        D = np.sqrt(np.diag(prec))
        partial = -prec / np.outer(D, D)
        np.fill_diagonal(partial, 0)
        partial_all[s] = partial

    return {"pearson": pearson_all, "partial": partial_all}


# ─────────────────────────────────────────────────────────────────────────────
# 3.  動的 FC (スライディングウィンドウ + HMM 状態推定)
# ─────────────────────────────────────────────────────────────────────────────
def compute_dynamic_fc(
    timeseries: np.ndarray,
    window_size: int = 40,
    step: int = 1,
    n_states: int = 4,
) -> dict:
    """
    スライディングウィンドウ法による動的FC計算と
    k-means クラスタリングによる脳状態推定。

    Parameters:
      window_size : ウィンドウ長 (TRs) — 典型値: 40TR = 80秒 (TR=2s)
      step        : ステップサイズ (TRs)
      n_states    : HMM/k-means 状態数
    """
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler

    n_subjects, n_rois, n_tp = timeseries.shape
    n_windows = (n_tp - window_size) // step + 1

    all_dfc = []
    for s in range(n_subjects):
        ts = timeseries[s]
        win_corrs = []
        for w in range(n_windows):
            start = w * step
            end = start + window_size
            seg = ts[:, start:end]
            corr = np.corrcoef(seg)
            np.fill_diagonal(corr, 0)
            # 上三角のみ
            idx = np.triu_indices(n_rois, k=1)
            win_corrs.append(corr[idx])
        all_dfc.append(win_corrs)

    all_dfc = np.array(all_dfc)  # (n_subjects, n_windows, n_edges)

    # 全被験者 × ウィンドウをフラット化してクラスタリング
    flat = all_dfc.reshape(-1, all_dfc.shape[-1])
    scaler = StandardScaler()
    flat_scaled = scaler.fit_transform(flat)

    rng = np.random.default_rng(42)
    km = KMeans(n_clusters=n_states, random_state=42, n_init=20)
    labels = km.fit_predict(flat_scaled)
    labels = labels.reshape(n_subjects, n_windows)

    # 各状態の統計
    state_stats = []
    for k in range(n_states):
        mask = (labels == k)
        frac = float(mask.mean())
        # 状態 k のFC行列（上三角→全行列復元）
        state_fc_vec = km.cluster_centers_[k]
        state_fc_mat = np.zeros((n_rois, n_rois))
        idx = np.triu_indices(n_rois, k=1)
        state_fc_mat[idx] = state_fc_vec
        state_fc_mat += state_fc_mat.T
        state_stats.append({
            "state": k + 1,
            "fraction_time": round(frac, 4),
            "mean_fc_strength": round(float(np.abs(state_fc_mat).mean()), 4),
        })

    # FC 変動性 (標準偏差)
    dfc_variability = all_dfc.std(axis=1).mean(axis=0)  # (n_edges,)
    mean_variability = float(dfc_variability.mean())

    return {
        "dfc_array": all_dfc,   # (n_subjects, n_windows, n_edges)
        "state_labels": labels,
        "state_stats": state_stats,
        "mean_dfc_variability": mean_variability,
        "n_windows": n_windows,
        "window_size_tr": window_size,
    }


# ─────────────────────────────────────────────────────────────────────────────
# 4.  FC 指標の群間比較
# ─────────────────────────────────────────────────────────────────────────────
def compare_fc_groups(fc_hc: np.ndarray, fc_dis: np.ndarray, group_name: str) -> dict:
    """
    HC vs 疾患群の FC 行列要素の t 検定 (FDR 補正)。
    """
    from scipy.stats import ttest_ind

    n_rois = fc_hc.shape[1]
    t_mat = np.zeros((n_rois, n_rois))
    p_mat = np.ones((n_rois, n_rois))

    for i in range(n_rois):
        for j in range(i + 1, n_rois):
            t, p = ttest_ind(fc_hc[:, i, j], fc_dis[:, i, j])
            t_mat[i, j] = t_mat[j, i] = t
            p_mat[i, j] = p_mat[j, i] = p

    # FDR (Benjamini-Hochberg)
    from scipy.stats import false_discovery_control
    p_upper = p_mat[np.triu_indices(n_rois, k=1)]
    p_fdr = false_discovery_control(p_upper, method='bh')
    n_sig = int((p_fdr < 0.05).sum())
    n_edges = len(p_fdr)

    return {
        "group_comparison": f"HC_vs_{group_name}",
        "n_edges_tested": n_edges,
        "n_significant_fdr05": n_sig,
        "pct_significant": round(n_sig / n_edges * 100, 2),
        "mean_t_stat": round(float(np.abs(t_mat[t_mat != 0]).mean()), 3),
        "cohen_d_mean": round(
            float(np.abs(fc_hc.mean(axis=0) - fc_dis.mean(axis=0)).mean()
                  / np.sqrt((fc_hc.std(axis=0)**2 + fc_dis.std(axis=0)**2) / 2 + 1e-9).mean()),
            3,
        ),
    }


# ─────────────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("[03] 機能的コネクティビティ解析中...")
    ts_start = datetime.utcnow().isoformat()

    # 時系列生成
    n_rois = 90
    groups_config = {"HC": 30, "SCZ": 25, "AD": 22}
    ts_data = {}
    fc_data = {}
    dfc_results = {}
    fc_metrics = {}

    for group, n_sub in groups_config.items():
        print(f"  [{group}] 時系列生成 ({n_sub}名) ...")
        ts = simulate_bold_timeseries(
            n_rois=n_rois, n_timepoints=300, tr=2.0,
            n_subjects=n_sub, group=group,
        )
        ts_data[group] = ts

        # 静的 FC
        fc = compute_static_fc(ts)
        fc_data[group] = fc
        np.save(f"{DATA_DIR}/FC_pearson_{group}.npy", fc["pearson"])
        np.save(f"{DATA_DIR}/FC_partial_{group}.npy", fc["partial"])

        # 動的 FC
        dfc = compute_dynamic_fc(ts, window_size=40, step=2, n_states=4)
        dfc_results[group] = dfc

        fc_metrics[group] = {
            "pearson_mean": float(fc["pearson"].mean(axis=0)[
                np.triu_indices(n_rois, k=1)].mean()),
            "pearson_std": float(fc["pearson"].std(axis=0)[
                np.triu_indices(n_rois, k=1)].mean()),
            "partial_mean_abs": float(np.abs(fc["partial"]).mean()),
            "dfc_variability": dfc["mean_dfc_variability"],
            "dfc_states": dfc["state_stats"],
        }
        print(f"    → 平均 Pearson FC={fc_metrics[group]['pearson_mean']:.3f}, "
              f"DFC 変動性={dfc['mean_dfc_variability']:.4f}")

    # 群間比較
    comparison_results = {}
    for dis_group in ["SCZ", "AD"]:
        comp = compare_fc_groups(
            fc_data["HC"]["pearson"],
            fc_data[dis_group]["pearson"],
            dis_group,
        )
        comparison_results[f"HC_vs_{dis_group}"] = comp
        print(f"  [HC vs {dis_group}] 有意エッジ: {comp['n_significant_fdr05']}/{comp['n_edges_tested']} "
              f"({comp['pct_significant']}%), Cohen d={comp['cohen_d_mean']}")

    # 保存
    with open(f"{RESULTS_DIR}/fc_metrics.json", "w") as f:
        json.dump(fc_metrics, f, indent=2, ensure_ascii=False, default=str)
    with open(f"{RESULTS_DIR}/fc_group_comparison.json", "w") as f:
        json.dump(comparison_results, f, indent=2, ensure_ascii=False)

    # ログ
    log_entry = {
        "timestamp": ts_start,
        "phase": "functional_connectivity",
        "event_type": "step_completed",
        "actor": "co-scientist",
        "skill_or_tool": "03_functional_connectivity.py",
        "handoff_out": {"fc_metrics": fc_metrics, "group_comparisons": comparison_results},
        "files_written": [
            f"{DATA_DIR}/FC_pearson_HC.npy",
            f"{DATA_DIR}/FC_partial_HC.npy",
            f"{DATA_DIR}/FC_pearson_SCZ.npy",
            f"{DATA_DIR}/FC_partial_SCZ.npy",
            f"{DATA_DIR}/FC_pearson_AD.npy",
            f"{DATA_DIR}/FC_partial_AD.npy",
            f"{RESULTS_DIR}/fc_metrics.json",
            f"{RESULTS_DIR}/fc_group_comparison.json",
        ],
        "status": "ok",
    }
    with open(f"{LOGS_DIR}/process-log.jsonl", "a") as f:
        f.write(json.dumps(log_entry, ensure_ascii=False) + "\n")

    return fc_data, dfc_results, ts_data


if __name__ == "__main__":
    main()
