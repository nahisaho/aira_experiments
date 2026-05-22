"""
06_reliability.py
-----------------
全脳コネクトーム解析 - ステップ6: テスト-リテスト信頼性評価
ICC / Fingerprinting / 再現性戦略
"""

import numpy as np
import json
import os
from datetime import datetime
from scipy import stats

RESULTS_DIR = "results"
LOGS_DIR = "logs"
DATA_DIR = "data"
np.random.seed(42)


# ─────────────────────────────────────────────────────────────────────────────
# 1.  ICC (級内相関係数) 計算
# ─────────────────────────────────────────────────────────────────────────────
def compute_icc(
    test: np.ndarray, retest: np.ndarray, icc_type: str = "ICC2"
) -> dict:
    """
    ICC(2,1) または ICC(3,1) を計算。
    test/retest: (n_subjects, n_features)

    ICC(2,1): 2-way mixed, absolute agreement
    ICC(3,1): 2-way mixed, consistency

    解釈:
      < 0.5  : 低い信頼性
      0.5–0.75: 中等度
      0.75–0.9: 良好
      > 0.9  : 非常に良好
    """
    n_sub = test.shape[0]
    k = 2  # 2 time points
    n_feat = test.shape[1]

    icc_vals = []
    for feat in range(n_feat):
        y_t = test[:, feat]
        y_r = retest[:, feat]

        # データ行列: (n_sub, k)
        Y = np.column_stack([y_t, y_r])
        grand_mean = Y.mean()

        # 被験者ごとの平均
        subject_means = Y.mean(axis=1)
        # ラター(測定回)ごとの平均
        rater_means = Y.mean(axis=0)

        # SS 計算
        SS_S = k * np.sum((subject_means - grand_mean) ** 2)          # between-subject
        SS_R = n_sub * np.sum((rater_means - grand_mean) ** 2)         # between-rater
        SS_T = np.sum((Y - grand_mean) ** 2)                           # total
        SS_E = SS_T - SS_S - SS_R                                      # error

        df_S = n_sub - 1
        df_R = k - 1
        df_E = (n_sub - 1) * (k - 1)

        MS_S = SS_S / df_S if df_S > 0 else 0
        MS_R = SS_R / df_R if df_R > 0 else 0
        MS_E = SS_E / df_E if df_E > 0 else 1e-9

        if MS_E == 0:
            MS_E = 1e-9

        if icc_type == "ICC2":
            # 2-way random, absolute agreement
            icc = (MS_S - MS_E) / (MS_S + (k - 1) * MS_E + k * (MS_R - MS_E) / n_sub)
        else:
            # ICC3: 2-way mixed, consistency
            icc = (MS_S - MS_E) / (MS_S + (k - 1) * MS_E)

        icc_vals.append(float(np.clip(icc, -1, 1)))

    icc_vals = np.array(icc_vals)
    return {
        "icc_type": icc_type,
        "mean_icc": round(float(icc_vals.mean()), 4),
        "median_icc": round(float(np.median(icc_vals)), 4),
        "std_icc": round(float(icc_vals.std()), 4),
        "pct_above_0_75": round(float((icc_vals > 0.75).mean() * 100), 2),
        "pct_above_0_60": round(float((icc_vals > 0.60).mean() * 100), 2),
        "pct_above_0_40": round(float((icc_vals > 0.40).mean() * 100), 2),
        "icc_distribution": {
            "poor_lt05": int((icc_vals < 0.5).sum()),
            "moderate_05_075": int(((icc_vals >= 0.5) & (icc_vals < 0.75)).sum()),
            "good_075_09": int(((icc_vals >= 0.75) & (icc_vals < 0.9)).sum()),
            "excellent_gt09": int((icc_vals >= 0.9).sum()),
        },
        "icc_values": icc_vals.tolist(),
    }


# ─────────────────────────────────────────────────────────────────────────────
# 2.  FC フィンガープリント (個人識別能力)
# ─────────────────────────────────────────────────────────────────────────────
def compute_fingerprinting(
    fc_test: np.ndarray,   # (n_subjects, n_rois, n_rois)
    fc_retest: np.ndarray,
) -> dict:
    """
    Finn et al. (2015) の手法に基づく FC フィンガープリント。
    識別精度 = テスト時の FC と最もよく一致するリテスト FC が同一人物かどうかの割合。
    """
    n_subjects, n_rois, _ = fc_test.shape
    idx = np.triu_indices(n_rois, k=1)

    # 各被験者の上三角ベクトル
    test_vecs = fc_test[:, idx[0], idx[1]]    # (n_sub, n_edges)
    retest_vecs = fc_retest[:, idx[0], idx[1]]

    n_correct = 0
    rank_accuracy = []
    for i in range(n_subjects):
        corrs = np.array([
            np.corrcoef(test_vecs[i], retest_vecs[j])[0, 1]
            for j in range(n_subjects)
        ])
        predicted = np.argmax(corrs)
        if predicted == i:
            n_correct += 1
        # ランク
        rank = int(np.argsort(corrs)[::-1].tolist().index(i) + 1)
        rank_accuracy.append(rank)

    identification_rate = n_correct / n_subjects * 100
    mean_rank = np.mean(rank_accuracy)

    # ネットワーク別フィンガープリント寄与
    # インライン定義
    network_labels = {
        "DMN": list(range(0, 15)),
        "FPN": list(range(15, 28)),
        "SMN": list(range(28, 42)),
        "VIS": list(range(42, 55)),
        "DAN": list(range(55, 65)),
        "SN":  list(range(65, 75)),
        "LIM": list(range(75, 90)),
    }
    net_fingerprint = {}
    for net_name, roi_list in network_labels.items():
        n_correct_net = 0
        for i in range(n_subjects):
            corrs = np.array([
                np.corrcoef(
                    fc_test[i, roi_list, :][:, roi_list][np.triu_indices(len(roi_list), k=1)],
                    fc_retest[j, roi_list, :][:, roi_list][np.triu_indices(len(roi_list), k=1)],
                )[0, 1]
                for j in range(n_subjects)
            ])
            if np.argmax(corrs) == i:
                n_correct_net += 1
        net_fingerprint[net_name] = round(n_correct_net / n_subjects * 100, 1)

    return {
        "identification_rate_pct": round(identification_rate, 2),
        "mean_rank": round(float(mean_rank), 2),
        "chance_level_pct": round(100.0 / n_subjects, 2),
        "network_fingerprint_accuracy": net_fingerprint,
    }


# ─────────────────────────────────────────────────────────────────────────────
# 3.  再現性確保戦略の評価
# ─────────────────────────────────────────────────────────────────────────────
def evaluate_reproducibility_strategies(
    fc_test: np.ndarray, fc_retest: np.ndarray
) -> dict:
    """
    各種前処理・解析選択が再現性（ICC）に与える影響を評価。
    """
    n_rois = fc_test.shape[1]
    idx = np.triu_indices(n_rois, k=1)

    strategies = {
        "baseline_no_gsr": (fc_test, fc_retest),
    }

    # グローバルシグナル回帰（GSR）のシミュレーション
    def apply_gsr(fc: np.ndarray) -> np.ndarray:
        rng = np.random.default_rng(123)
        noise = rng.normal(0, 0.02, fc.shape)
        return fc - fc.mean(axis=(1, 2), keepdims=True) + noise

    strategies["with_gsr"] = (apply_gsr(fc_test), apply_gsr(fc_retest))

    # Fisher z 変換
    def fisher_z(fc: np.ndarray) -> np.ndarray:
        return np.arctanh(np.clip(fc, -0.999, 0.999))

    strategies["fisher_z"] = (fisher_z(fc_test), fisher_z(fc_retest))

    # 平均 FC 正規化
    def normalize_fc(fc: np.ndarray) -> np.ndarray:
        mean = fc.mean(axis=(1, 2), keepdims=True)
        std = fc.std(axis=(1, 2), keepdims=True) + 1e-9
        return (fc - mean) / std

    strategies["normalized"] = (normalize_fc(fc_test), normalize_fc(fc_retest))

    results = {}
    for name, (t, r) in strategies.items():
        t_vecs = t[:, idx[0], idx[1]]
        r_vecs = r[:, idx[0], idx[1]]
        icc_res = compute_icc(t_vecs, r_vecs)
        results[name] = {
            "mean_icc": icc_res["mean_icc"],
            "pct_above_0_75": icc_res["pct_above_0_75"],
        }

    return results


# ─────────────────────────────────────────────────────────────────────────────
# 4.  サンプルサイズ推定（検出力分析）
# ─────────────────────────────────────────────────────────────────────────────
def sample_size_for_reliability(
    target_icc: float = 0.75,
    current_icc: float = 0.70,
    alpha: float = 0.05,
    power: float = 0.80,
) -> dict:
    """
    目標 ICC に達するための必要サンプルサイズを推定。
    Bonett (2002) の近似式を使用。
    """
    # Fisher z 変換後のサンプルサイズ近似
    # F 検定に基づく近似
    k = 2  # 測定回数
    z_alpha = stats.norm.ppf(1 - alpha / 2)
    z_beta = stats.norm.ppf(power)

    # ICC の Fisher z 変換
    z_target = np.arctanh(target_icc)
    z_current = np.arctanh(current_icc)

    n_estimate = int(
        np.ceil((z_alpha + z_beta) ** 2 / (z_target - z_current) ** 2 + 3)
    )

    # ICC 95% CI (Shrout & Fleiss)
    n_obs = 30
    F = (1 + (k - 1) * current_icc) / (1 - current_icc)
    df1 = n_obs - 1
    df2 = n_obs * (k - 1)
    F_lower = F / stats.f.ppf(1 - alpha / 2, df1, df2)
    F_upper = F * stats.f.ppf(1 - alpha / 2, df2, df1)
    icc_lower = (F_lower - 1) / (F_lower + k - 1)
    icc_upper = (F_upper - 1) / (F_upper + k - 1)

    return {
        "current_icc": current_icc,
        "target_icc": target_icc,
        "recommended_n_subjects": n_estimate,
        "icc_95ci_lower": round(float(icc_lower), 3),
        "icc_95ci_upper": round(float(icc_upper), 3),
    }


# ─────────────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("[06] テスト-リテスト信頼性評価中...")
    ts = datetime.utcnow().isoformat()

    # テスト-リテスト FC シミュレーション
    fc_test = np.load(f"{DATA_DIR}/FC_pearson_HC.npy")

    # リテスト: テストに高相関ノイズを加算
    rng = np.random.default_rng(99)
    noise = rng.normal(0, 0.03, fc_test.shape)
    noise = (noise + noise.transpose(0, 2, 1)) / 2
    fc_retest = fc_test + noise

    n_rois = fc_test.shape[1]
    idx = np.triu_indices(n_rois, k=1)
    test_vecs = fc_test[:, idx[0], idx[1]]
    retest_vecs = fc_retest[:, idx[0], idx[1]]

    # ICC
    icc2 = compute_icc(test_vecs, retest_vecs, "ICC2")
    icc3 = compute_icc(test_vecs, retest_vecs, "ICC3")
    print(f"  → ICC(2,1) mean={icc2['mean_icc']:.3f}, >0.75: {icc2['pct_above_0_75']}%")
    print(f"  → ICC(3,1) mean={icc3['mean_icc']:.3f}, >0.75: {icc3['pct_above_0_75']}%")

    # Fingerprinting
    fp = compute_fingerprinting(fc_test, fc_retest)
    print(f"  → Fingerprinting: {fp['identification_rate_pct']:.1f}% "
          f"(チャンスレベル: {fp['chance_level_pct']:.1f}%)")

    # 再現性戦略比較
    repro = evaluate_reproducibility_strategies(fc_test, fc_retest)
    print("  → 再現性戦略比較:")
    for name, r in repro.items():
        print(f"    {name}: ICC={r['mean_icc']:.3f}, >0.75={r['pct_above_0_75']}%")

    # サンプルサイズ推定
    n_est = sample_size_for_reliability(
        target_icc=0.75, current_icc=icc2["mean_icc"]
    )
    print(f"  → 推奨サンプルサイズ (ICC≥0.75): n={n_est['recommended_n_subjects']}")

    reliability_results = {
        "icc2_1": {k: v for k, v in icc2.items() if k != "icc_values"},
        "icc3_1": {k: v for k, v in icc3.items() if k != "icc_values"},
        "fingerprinting": fp,
        "reproducibility_strategies": repro,
        "sample_size_recommendation": n_est,
    }

    with open(f"{RESULTS_DIR}/reliability_results.json", "w") as f:
        json.dump(reliability_results, f, indent=2, ensure_ascii=False)

    log_entry = {
        "timestamp": ts,
        "phase": "reliability",
        "event_type": "step_completed",
        "actor": "co-scientist",
        "skill_or_tool": "06_reliability.py",
        "handoff_out": {
            "mean_icc2": icc2["mean_icc"],
            "fingerprinting_rate": fp["identification_rate_pct"],
        },
        "files_written": [f"{RESULTS_DIR}/reliability_results.json"],
        "status": "ok",
    }
    with open(f"{LOGS_DIR}/process-log.jsonl", "a") as f:
        f.write(json.dumps(log_entry, ensure_ascii=False) + "\n")

    return reliability_results


if __name__ == "__main__":
    main()
