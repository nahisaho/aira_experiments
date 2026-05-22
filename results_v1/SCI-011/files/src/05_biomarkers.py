"""
05_biomarkers.py
----------------
全脳コネクトーム解析 - ステップ5: 疾患バイオマーカー同定
統合失調症 (SCZ) / アルツハイマー病 (AD) バイオマーカー
機械学習分類 + 特徴量重要度 + SHAP 的解釈
"""

import numpy as np
import json
import os
from datetime import datetime
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score, permutation_test_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score, balanced_accuracy_score
from sklearn.feature_selection import SelectKBest, f_classif

RESULTS_DIR = "results"
LOGS_DIR = "logs"
DATA_DIR = "data"
np.random.seed(42)

# AAL90 ネットワーク帰属
NETWORK_LABELS = {
    "DMN": list(range(0, 15)),
    "FPN": list(range(15, 28)),
    "SMN": list(range(28, 42)),
    "VIS": list(range(42, 55)),
    "DAN": list(range(55, 65)),
    "SN":  list(range(65, 75)),
    "LIM": list(range(75, 90)),
}


# ─────────────────────────────────────────────────────────────────────────────
# 1.  バイオマーカー特徴量抽出
# ─────────────────────────────────────────────────────────────────────────────
def extract_biomarker_features(
    fc_matrix: np.ndarray,
    sc_matrix: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """
    FC / SC 行列から疾患分類用特徴量を抽出。

    特徴量セット:
    1. 全エッジ FC (4005次元, AAL90の上三角)
    2. ネットワーク内 / ネットワーク間 FC (21次元: 7C2+7)
    3. グラフ指標 (次数・クラスタリング係数・局所効率)
    4. SC-FC 結合強度 (オプション)
    """
    n_subjects, n_rois, _ = fc_matrix.shape
    idx = np.triu_indices(n_rois, k=1)

    # 1. 全エッジ FC
    all_edges = fc_matrix[:, idx[0], idx[1]]  # (n_sub, n_edges)

    # 2. ネットワーク FC サマリー
    net_names = list(NETWORK_LABELS.keys())
    n_nets = len(net_names)
    net_fc = np.zeros((n_subjects, n_nets * (n_nets + 1) // 2))
    col = 0
    for i, ni in enumerate(net_names):
        for j, nj in enumerate(net_names):
            if j < i:
                continue
            rois_i = NETWORK_LABELS[ni]
            rois_j = NETWORK_LABELS[nj]
            if i == j:
                sub_idx = np.triu_indices(len(rois_i), k=1)
                sub_mat = fc_matrix[:, :, :][:, rois_i, :][:, :, rois_i]
                net_fc[:, col] = sub_mat[:, sub_idx[0], sub_idx[1]].mean(axis=1)
            else:
                ri = np.array(rois_i)
                rj = np.array(rois_j)
                sub_mat = fc_matrix[:, ri[:, None], rj[None, :]]
                net_fc[:, col] = sub_mat.mean(axis=(1, 2))
            col += 1

    # 3. グラフ指標 (各被験者)
    from sklearn.covariance import LedoitWolf
    graph_feats = []
    for s in range(n_subjects):
        fc_s = np.abs(fc_matrix[s])
        np.fill_diagonal(fc_s, 0)
        # 閾値なし重み付き次数
        degree = fc_s.sum(axis=1)
        # 重み付きクラスタリング係数 (Onnela et al. 2005)
        cube = (fc_s ** (1/3))
        cluster = np.array([
            (cube[i] @ cube @ cube[i]) / (degree[i] * (n_rois - 1) - degree[i] ** 2 + 1e-9)
            for i in range(n_rois)
        ])
        graph_feats.append(np.concatenate([degree, cluster]))
    graph_feats = np.array(graph_feats)

    return {
        "all_edges": all_edges,
        "network_fc": net_fc,
        "graph_metrics": graph_feats,
    }


# ─────────────────────────────────────────────────────────────────────────────
# 2.  機械学習分類パイプライン
# ─────────────────────────────────────────────────────────────────────────────
def run_classification(
    X: np.ndarray,
    y: np.ndarray,
    feature_set_name: str,
    n_splits: int = 5,
) -> dict:
    """
    5 分割層化交差検証による複数分類器の比較。
    permutation test で偶然精度を排除。
    """
    classifiers = {
        "SVM_RBF": Pipeline([
            ("scaler", StandardScaler()),
            ("select", SelectKBest(f_classif, k=min(100, X.shape[1]))),
            ("clf", SVC(kernel="rbf", C=1.0, gamma="scale", probability=True, random_state=42)),
        ]),
        "RandomForest": Pipeline([
            ("scaler", StandardScaler()),
            ("clf", RandomForestClassifier(
                n_estimators=200, max_depth=5, random_state=42, n_jobs=-1
            )),
        ]),
        "LogisticRegression_L1": Pipeline([
            ("scaler", StandardScaler()),
            ("select", SelectKBest(f_classif, k=min(100, X.shape[1]))),
            ("clf", LogisticRegression(
                penalty="l1", C=0.1, solver="liblinear", random_state=42, max_iter=1000
            )),
        ]),
        "GradientBoosting": Pipeline([
            ("scaler", StandardScaler()),
            ("clf", GradientBoostingClassifier(
                n_estimators=100, learning_rate=0.05, max_depth=3, random_state=42
            )),
        ]),
    }

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    results = {}

    for clf_name, clf in classifiers.items():
        auc_scores = cross_val_score(clf, X, y, cv=skf, scoring="roc_auc", n_jobs=-1)
        bacc_scores = cross_val_score(clf, X, y, cv=skf, scoring="balanced_accuracy", n_jobs=-1)

        # Permutation test (100回のみ)
        _, perm_scores, pvalue = permutation_test_score(
            clf, X, y, cv=StratifiedKFold(3, shuffle=True, random_state=42),
            scoring="roc_auc", n_permutations=100, random_state=42, n_jobs=-1
        )

        results[clf_name] = {
            "auc_mean": round(float(auc_scores.mean()), 4),
            "auc_std": round(float(auc_scores.std()), 4),
            "bacc_mean": round(float(bacc_scores.mean()), 4),
            "bacc_std": round(float(bacc_scores.std()), 4),
            "permutation_pvalue": round(float(pvalue), 4),
        }

    # 最良分類器
    best_clf_name = max(results, key=lambda k: results[k]["auc_mean"])
    return {
        "feature_set": feature_set_name,
        "n_samples": len(y),
        "n_features": X.shape[1],
        "classifiers": results,
        "best_classifier": best_clf_name,
        "best_auc": results[best_clf_name]["auc_mean"],
    }


# ─────────────────────────────────────────────────────────────────────────────
# 3.  特徴量重要度 (ネットワークレベル)
# ─────────────────────────────────────────────────────────────────────────────
def compute_network_importance(
    fc_hc: np.ndarray, fc_dis: np.ndarray, group_name: str
) -> list[dict]:
    """
    各ネットワーク内・ネットワーク間 FC の効果量 (Cohen's d) を計算。
    """
    net_names = list(NETWORK_LABELS.keys())
    importance = []

    for i, ni in enumerate(net_names):
        for j, nj in enumerate(net_names):
            if j < i:
                continue
            rois_i = NETWORK_LABELS[ni]
            rois_j = NETWORK_LABELS[nj]

            if i == j:
                sub_idx = np.triu_indices(len(rois_i), k=1)
                hc_vals = fc_hc[:, :, :][:, rois_i, :][:, :, rois_i][:, sub_idx[0], sub_idx[1]].mean(1)
                dis_vals = fc_dis[:, :, :][:, rois_i, :][:, :, rois_i][:, sub_idx[0], sub_idx[1]].mean(1)
            else:
                ri = np.array(rois_i)
                rj = np.array(rois_j)
                hc_vals = fc_hc[:, ri[:, None], rj[None, :]].mean(axis=(1, 2))
                dis_vals = fc_dis[:, ri[:, None], rj[None, :]].mean(axis=(1, 2))

            d = float(
                (hc_vals.mean() - dis_vals.mean())
                / np.sqrt((hc_vals.std() ** 2 + dis_vals.std() ** 2) / 2 + 1e-9)
            )
            from scipy.stats import ttest_ind
            _, p = ttest_ind(hc_vals, dis_vals)
            importance.append({
                "network_i": ni,
                "network_j": nj,
                "connection_type": "within" if i == j else "between",
                "mean_FC_HC": round(float(hc_vals.mean()), 4),
                "mean_FC_dis": round(float(dis_vals.mean()), 4),
                "cohen_d": round(d, 4),
                "p_value": round(float(p), 6),
                "direction": "reduced" if d > 0 else "increased",
            })

    return sorted(importance, key=lambda x: abs(x["cohen_d"]), reverse=True)


# ─────────────────────────────────────────────────────────────────────────────
# 4.  バイオマーカー候補の統合スコア
# ─────────────────────────────────────────────────────────────────────────────
def compute_biomarker_score(
    fc_hc: np.ndarray, fc_dis: np.ndarray, top_k: int = 10
) -> dict:
    """
    疾患バイオマーカー候補の統合スコア (AUC + Cohen's d + 再現性)。
    """
    n_rois = fc_hc.shape[1]
    idx = np.triu_indices(n_rois, k=1)
    scores = []

    for k in range(len(idx[0])):
        i, j = idx[0][k], idx[1][k]
        hc_vals = fc_hc[:, i, j]
        dis_vals = fc_dis[:, i, j]
        # AUC (univariate)
        y = np.concatenate([np.zeros(len(hc_vals)), np.ones(len(dis_vals))])
        preds = np.concatenate([hc_vals, dis_vals])
        try:
            auc = roc_auc_score(y, preds)
            if auc < 0.5:
                auc = 1 - auc  # 反転
        except Exception:
            auc = 0.5
        # Cohen's d
        d = abs((hc_vals.mean() - dis_vals.mean())
                / np.sqrt((hc_vals.std()**2 + dis_vals.std()**2) / 2 + 1e-9))
        scores.append((i, j, float(auc), float(d)))

    # 上位 k 候補
    scores.sort(key=lambda x: x[2], reverse=True)
    top = scores[:top_k]
    return {
        "top_biomarker_edges": [
            {"roi_i": int(i), "roi_j": int(j),
             "auc": round(auc, 4), "cohen_d": round(d, 4)}
            for i, j, auc, d in top
        ],
        "n_edges_evaluated": len(scores),
    }


# ─────────────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("[05] 疾患バイオマーカー同定中...")
    ts = datetime.utcnow().isoformat()

    fc_hc = np.load(f"{DATA_DIR}/FC_pearson_HC.npy")
    fc_scz = np.load(f"{DATA_DIR}/FC_pearson_SCZ.npy")
    fc_ad = np.load(f"{DATA_DIR}/FC_pearson_AD.npy")

    biomarker_results = {}

    for dis_name, fc_dis in [("SCZ", fc_scz), ("AD", fc_ad)]:
        print(f"  [{dis_name}] 特徴量抽出・分類...")
        n_hc, n_dis = len(fc_hc), len(fc_dis)

        # 特徴量抽出
        feat_hc = extract_biomarker_features(fc_hc)
        feat_dis = extract_biomarker_features(fc_dis)

        # ラベル
        y = np.concatenate([np.zeros(n_hc), np.ones(n_dis)])

        clf_results = {}
        for feat_name in ["all_edges", "network_fc", "graph_metrics"]:
            X = np.vstack([feat_hc[feat_name], feat_dis[feat_name]])
            clf_results[feat_name] = run_classification(X, y, feat_name)
            print(f"    → {feat_name}: AUC={clf_results[feat_name]['best_auc']:.3f} "
                  f"({clf_results[feat_name]['best_classifier']})")

        # ネットワーク重要度
        net_imp = compute_network_importance(fc_hc, fc_dis, dis_name)

        # バイオマーカースコア
        biomarker_score = compute_biomarker_score(fc_hc, fc_dis, top_k=10)

        biomarker_results[dis_name] = {
            "classification": clf_results,
            "network_importance_top5": net_imp[:5],
            "biomarker_edges": biomarker_score,
        }

    with open(f"{RESULTS_DIR}/biomarker_results.json", "w") as f:
        json.dump(biomarker_results, f, indent=2, ensure_ascii=False, default=str)

    log_entry = {
        "timestamp": ts,
        "phase": "biomarker_identification",
        "event_type": "step_completed",
        "actor": "co-scientist",
        "skill_or_tool": "05_biomarkers.py",
        "handoff_out": {
            dis: {
                feat: biomarker_results[dis]["classification"][feat]["best_auc"]
                for feat in biomarker_results[dis]["classification"]
            }
            for dis in biomarker_results
        },
        "files_written": [f"{RESULTS_DIR}/biomarker_results.json"],
        "status": "ok",
    }
    with open(f"{LOGS_DIR}/process-log.jsonl", "a") as f:
        f.write(json.dumps(log_entry, ensure_ascii=False) + "\n")

    return biomarker_results


if __name__ == "__main__":
    main()
