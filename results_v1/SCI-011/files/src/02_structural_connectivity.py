"""
02_structural_connectivity.py
-----------------------------
全脳コネクトーム解析 - ステップ2: 構造的コネクティビティ
確率的トラクトグラフィー（FSL BEDPOSTX / MRtrix3 iFOD2）+ 構造的接続行列の推定
"""

import numpy as np
import json
import os
from datetime import datetime

RESULTS_DIR = "results"
LOGS_DIR = "logs"
DATA_DIR = "data"
np.random.seed(42)


# ─────────────────────────────────────────────────────────────────────────────
# 1.  トラクトグラフィーコマンドジェネレータ
# ─────────────────────────────────────────────────────────────────────────────
class TractographyPipeline:
    """
    FSL BEDPOSTX + probtrackx2 / MRtrix3 iFOD2 コマンドを生成する。
    """

    def __init__(self, subject_id: str, atlas: str = "AAL90", base_dir: str = "/data"):
        self.subject_id = subject_id
        self.atlas = atlas
        self.base_dir = base_dir

    def bedpostx_commands(self) -> list[str]:
        """
        FSL BEDPOSTX: 拡散テンソル多成分モデルフィッティング。
        nfibres=3: 最大3本の繊維方向をモデル化
        model=2: bingham分布（デフォルトより精度高）
        """
        sub = self.subject_id
        bd = self.base_dir
        return [
            # BEDPOSTX セットアップ
            f"bedpostx_datacheck {bd}/{sub}/dwi/",
            f"bedpostx {bd}/{sub}/dwi/"
            f" -n 3 -w 1 -b 1000 --model=2"
            f" --rician --noard",
        ]

    def probtrackx_commands(self, n_samples: int = 5000) -> list[str]:
        """
        probtrackx2: 確率的トラクトグラフィー（全 ROI ペア）。
        - omatrix2: 接続確率行列モード
        - loopcheck: ループ検出
        - --pd: 距離補正
        """
        sub = self.subject_id
        bd = self.base_dir
        return [
            f"probtrackx2"
            f" -x {bd}/{sub}/atlas/{self.atlas}_rois.nii.gz"
            f" -l --onewaycondition"
            f" -c 0.2 -S 2000 --steplength=0.5"
            f" -P {n_samples}"
            f" --fibthresh=0.01"
            f" --distthresh=0.0"
            f" --sampvox=0.0"
            f" --forcedir --opd --pd"
            f" --loopcheck"
            f" --omatrix2"
            f" --target2={bd}/{sub}/atlas/{self.atlas}_rois.nii.gz"
            f" -s {bd}/{sub}/dwi.bedpostX/merged"
            f" -m {bd}/{sub}/dwi.bedpostX/nodif_brain_mask"
            f" --dir={bd}/{sub}/probtrackx/{self.atlas}"
            f" --stop={bd}/{sub}/anat/gmwmi.nii.gz",
        ]

    def mrtrix3_commands(self, n_streamlines: int = 10_000_000) -> list[str]:
        """
        MRtrix3 iFOD2: 繊維方向分布（FOD）ベースの確率的トラクトグラフィー。
        SIFT2 で接続強度を重み付け。
        """
        sub = self.subject_id
        bd = self.base_dir
        return [
            # 応答関数推定 (Dhollander アルゴリズム)
            f"dwi2response dhollander"
            f" {bd}/{sub}/dwi/dwi_biascorr.mif"
            f" {bd}/{sub}/dwi/wm_response.txt"
            f" {bd}/{sub}/dwi/gm_response.txt"
            f" {bd}/{sub}/dwi/csf_response.txt",
            # FOD 推定 (CSD)
            f"dwi2fod msmt_csd"
            f" {bd}/{sub}/dwi/dwi_biascorr.mif"
            f" {bd}/{sub}/dwi/wm_response.txt {bd}/{sub}/dwi/wm_fod.mif"
            f" {bd}/{sub}/dwi/gm_response.txt {bd}/{sub}/dwi/gm_fod.mif"
            f" {bd}/{sub}/dwi/csf_response.txt {bd}/{sub}/dwi/csf_fod.mif"
            f" -mask {bd}/{sub}/dwi/nodif_brain_mask.nii.gz",
            # FOD 正規化
            f"mtnormalise"
            f" {bd}/{sub}/dwi/wm_fod.mif {bd}/{sub}/dwi/wm_fod_norm.mif"
            f" {bd}/{sub}/dwi/gm_fod.mif {bd}/{sub}/dwi/gm_fod_norm.mif"
            f" {bd}/{sub}/dwi/csf_fod.mif {bd}/{sub}/dwi/csf_fod_norm.mif"
            f" -mask {bd}/{sub}/dwi/nodif_brain_mask.nii.gz",
            # GM/WM 境界マスク生成 (5tt)
            f"5ttgen fsl {bd}/{sub}/anat/T1w.nii.gz {bd}/{sub}/anat/5tt.mif",
            f"5tt2gmwmi {bd}/{sub}/anat/5tt.mif {bd}/{sub}/anat/gmwmi.mif",
            # iFOD2 トラクトグラフィー
            f"tckgen {bd}/{sub}/dwi/wm_fod_norm.mif"
            f" {bd}/{sub}/tractography/tracks_{n_streamlines//1000000}M.tck"
            f" -act {bd}/{sub}/anat/5tt.mif"
            f" -backtrack -seed_gmwmi {bd}/{sub}/anat/gmwmi.mif"
            f" -maxlength 250 -minlength 10"
            f" -select {n_streamlines}",
            # SIFT2 重み付け
            f"tcksift2"
            f" {bd}/{sub}/tractography/tracks_{n_streamlines//1000000}M.tck"
            f" {bd}/{sub}/dwi/wm_fod_norm.mif"
            f" {bd}/{sub}/tractography/sift2_weights.txt"
            f" -act {bd}/{sub}/anat/5tt.mif",
            # 構造的接続行列生成
            f"tck2connectome"
            f" {bd}/{sub}/tractography/tracks_{n_streamlines//1000000}M.tck"
            f" {bd}/{sub}/atlas/{self.atlas}_MNI.nii.gz"
            f" {bd}/{sub}/connectome/SC_matrix_{self.atlas}.csv"
            f" -tck_weights_in {bd}/{sub}/tractography/sift2_weights.txt"
            f" -symmetric -zero_diagonal"
            f" -scale_invnodevol",
        ]


# ─────────────────────────────────────────────────────────────────────────────
# 2.  構造的接続行列のシミュレーション
# ─────────────────────────────────────────────────────────────────────────────
def simulate_structural_connectome(
    n_rois: int = 90,
    n_subjects: int = 30,
    group: str = "HC",
) -> np.ndarray:
    """
    AAL90 アトラスを用いた構造的コネクトーム（SC）行列のシミュレーション。
    スモールワールド・コミュニティ構造を模倣。
    """
    rng = np.random.default_rng({"HC": 42, "SCZ": 7, "AD": 13}[group])

    # グループ平均 SC 行列を構築
    avg_sc = np.zeros((n_rois, n_rois))

    # 短距離結合（近隣 ROI）: 強い接続
    for i in range(n_rois):
        for j in range(max(0, i - 5), min(n_rois, i + 6)):
            if i != j:
                avg_sc[i, j] = rng.exponential(0.8)

    # 長距離ハブ結合（デフォルトモードネットワーク等）
    hub_nodes = [7, 8, 25, 26, 30, 31, 67, 68]  # 前頭葉・側頭葉・頭頂葉
    for h in hub_nodes:
        for other in hub_nodes:
            if h != other:
                avg_sc[h, other] += rng.exponential(1.2)

    # 半球内左右対称性
    for i in range(n_rois // 2):
        j = i + n_rois // 2
        avg_sc[i, j] += rng.exponential(0.5)
        avg_sc[j, i] = avg_sc[i, j]

    # 疾患群: 接続減弱
    if group == "SCZ":
        # 前頭葉-側頭葉接続の減弱 (-30%)
        for i in range(20, 30):
            for j in range(50, 65):
                avg_sc[i, j] *= rng.uniform(0.5, 0.75)
                avg_sc[j, i] = avg_sc[i, j]
    elif group == "AD":
        # デフォルトモードネットワーク接続の減弱 (-40%)
        for h in hub_nodes:
            for k in hub_nodes:
                if h != k:
                    avg_sc[h, k] *= rng.uniform(0.45, 0.65)

    # 対称化・対角ゼロ化
    avg_sc = (avg_sc + avg_sc.T) / 2
    np.fill_diagonal(avg_sc, 0)

    # 被験者ごとのノイズを加えてスタック
    all_sc = []
    for _ in range(n_subjects):
        noise = rng.normal(0, 0.05, (n_rois, n_rois))
        noise = (noise + noise.T) / 2
        sc_sub = np.clip(avg_sc + avg_sc * noise, 0, None)
        np.fill_diagonal(sc_sub, 0)
        all_sc.append(sc_sub)

    return np.array(all_sc)  # (n_subjects, n_rois, n_rois)


def compute_sc_metrics(sc_matrices: np.ndarray) -> dict:
    """
    構造的コネクトーム行列の群レベル指標を計算。
    """
    avg = sc_matrices.mean(axis=0)
    metrics = {
        "mean_connectivity_strength": float(avg[avg > 0].mean()),
        "density": float((avg > 0).sum() / (avg.shape[0] * (avg.shape[0] - 1))),
        "n_rois": avg.shape[0],
        "coefficient_of_variation": float(
            sc_matrices.std(axis=0)[avg > 0].mean() / avg[avg > 0].mean()
        ),
    }
    return metrics


# ─────────────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("[02] 構造的コネクティビティ解析中...")
    ts = datetime.utcnow().isoformat()

    # コマンド生成
    tpipe = TractographyPipeline("sub-01")
    commands = {
        "bedpostx": tpipe.bedpostx_commands(),
        "probtrackx2": tpipe.probtrackx_commands(n_samples=5000),
        "mrtrix3_iFOD2": tpipe.mrtrix3_commands(n_streamlines=10_000_000),
    }
    with open(f"{RESULTS_DIR}/tractography_commands.json", "w") as f:
        json.dump(commands, f, indent=2, ensure_ascii=False)

    # SC 行列シミュレーション（3群）
    sc_data = {}
    sc_metrics_all = {}
    for group in ["HC", "SCZ", "AD"]:
        n_sub = {"HC": 30, "SCZ": 25, "AD": 22}[group]
        sc = simulate_structural_connectome(n_rois=90, n_subjects=n_sub, group=group)
        sc_data[group] = sc
        sc_metrics_all[group] = compute_sc_metrics(sc)
        # 平均行列を保存
        np.save(f"{DATA_DIR}/SC_avg_{group}.npy", sc.mean(axis=0))
        print(f"  → {group}: {n_sub}名, 密度={sc_metrics_all[group]['density']:.3f}")

    with open(f"{RESULTS_DIR}/sc_metrics.json", "w") as f:
        json.dump(sc_metrics_all, f, indent=2, ensure_ascii=False)

    # ログ
    log_entry = {
        "timestamp": ts,
        "phase": "structural_connectivity",
        "event_type": "step_completed",
        "actor": "co-scientist",
        "skill_or_tool": "02_structural_connectivity.py",
        "handoff_out": {"sc_metrics": sc_metrics_all},
        "files_written": [
            f"{RESULTS_DIR}/tractography_commands.json",
            f"{RESULTS_DIR}/sc_metrics.json",
            f"{DATA_DIR}/SC_avg_HC.npy",
            f"{DATA_DIR}/SC_avg_SCZ.npy",
            f"{DATA_DIR}/SC_avg_AD.npy",
        ],
        "status": "ok",
    }
    with open(f"{LOGS_DIR}/process-log.jsonl", "a") as f:
        f.write(json.dumps(log_entry, ensure_ascii=False) + "\n")

    return sc_data


if __name__ == "__main__":
    main()
