"""
01_preprocessing.py
-------------------
全脳コネクトーム解析 - ステップ1: 前処理パイプライン
fMRI/dMRI 前処理シミュレーション（FSL/FreeSurfer コマンド生成 + NumPy デモ実装）
"""

import numpy as np
import json
import os
from datetime import datetime

# ── 出力ディレクトリ ──────────────────────────────────────────────────────────
RESULTS_DIR = "results"
LOGS_DIR = "logs"
DATA_DIR = "data"
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(LOGS_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

np.random.seed(42)


# ─────────────────────────────────────────────────────────────────────────────
# 1.  FSL/FreeSurfer コマンドジェネレータ
# ─────────────────────────────────────────────────────────────────────────────
class PreprocessingPipeline:
    """
    FSL / FreeSurfer ベースの前処理パイプラインコマンドを生成し、
    最適パラメータ選定ロジックを示す。
    """

    def __init__(self, subject_id: str, base_dir: str = "/data"):
        self.subject_id = subject_id
        self.base_dir = base_dir
        self.log: list[dict] = []

    # ── fMRI 前処理 ────────────────────────────────────────────────────────
    def fmri_preprocessing_commands(self, tr: float = 2.0) -> dict[str, list[str]]:
        """
        fMRI 前処理コマンドを段階別に返す。
        最適パラメータ:
          - MCFLIRT: cost=normcorr, smooth=1mm, refvol=middle
          - BET: f=0.3 (脳抽出閾値)
          - FLIRT→FNIRT: MNI152 2mm 標準化
          - SUSAN: FWHM=5mm スムージング
        """
        sub = self.subject_id
        bd = self.base_dir

        commands = {
            "1_motion_correction": [
                f"mcflirt -in {bd}/{sub}/func/bold.nii.gz"
                f" -out {bd}/{sub}/func/bold_mc"
                f" -cost normcorr -smooth 1.0 -meanvol -plots"
                f" -refvol middle",
                # 動き量パラメータを抽出
                f"fsl_motion_outliers -i {bd}/{sub}/func/bold_mc.nii.gz"
                f" -o {bd}/{sub}/func/confounds_motionoutliers.txt"
                f" --dvars --nomoco",
            ],
            "2_brain_extraction": [
                f"bet {bd}/{sub}/anat/T1w.nii.gz"
                f" {bd}/{sub}/anat/T1w_brain.nii.gz"
                f" -f 0.3 -R",
            ],
            "3_distortion_correction": [
                # B0 フィールドマップ使用 (topup)
                f"topup --imain={bd}/{sub}/fmap/AP_PA.nii.gz"
                f" --datain={bd}/{sub}/fmap/acqparams.txt"
                f" --config=b02b0.cnf"
                f" --out={bd}/{sub}/fmap/topup_results"
                f" --fout={bd}/{sub}/fmap/field_hz"
                f" --iout={bd}/{sub}/fmap/unwarped",
                f"applytopup --imain={bd}/{sub}/func/bold_mc.nii.gz"
                f" --inindex=1"
                f" --datain={bd}/{sub}/fmap/acqparams.txt"
                f" --topup={bd}/{sub}/fmap/topup_results"
                f" --out={bd}/{sub}/func/bold_mc_dc.nii.gz",
            ],
            "4_registration_MNI": [
                # T1 → MNI152 線形登録
                f"flirt -in {bd}/{sub}/anat/T1w_brain.nii.gz"
                f" -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz"
                f" -omat {bd}/{sub}/reg/T1_to_MNI_affine.mat"
                f" -cost corratio -dof 12",
                # 非線形登録 (FNIRT)
                f"fnirt --in={bd}/{sub}/anat/T1w.nii.gz"
                f" --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz"
                f" --aff={bd}/{sub}/reg/T1_to_MNI_affine.mat"
                f" --cout={bd}/{sub}/reg/T1_to_MNI_warp"
                f" --config=T1_2_MNI152_2mm"
                f" --warpres=10,10,10",
                # fMRI → MNI152 (BOLD→T1→MNI)
                f"flirt -in {bd}/{sub}/func/bold_mc_dc.nii.gz"
                f" -ref {bd}/{sub}/anat/T1w_brain.nii.gz"
                f" -omat {bd}/{sub}/reg/BOLD_to_T1.mat"
                f" -cost bbr -dof 6 -wmseg {bd}/{sub}/anat/wmseg.nii.gz",
                f"applywarp --in={bd}/{sub}/func/bold_mc_dc.nii.gz"
                f" --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz"
                f" --premat={bd}/{sub}/reg/BOLD_to_T1.mat"
                f" --warp={bd}/{sub}/reg/T1_to_MNI_warp.nii.gz"
                f" --out={bd}/{sub}/func/bold_MNI.nii.gz"
                f" --interp=spline",
            ],
            "5_nuisance_regression": [
                # CSF・WM・グローバルシグナル・運動パラメータの回帰
                f"fsl_regfilt -i {bd}/{sub}/func/bold_MNI.nii.gz"
                f" -d {bd}/{sub}/func/confounds.txt"
                f" -o {bd}/{sub}/func/bold_cleaned.nii.gz"
                f" -f '1,2,3,4,5,6,7,8,9'",
            ],
            "6_temporal_filtering": [
                # バンドパスフィルタ (0.01–0.1 Hz for resting-state)
                f"fslmaths {bd}/{sub}/func/bold_cleaned.nii.gz"
                f" -bptf {int(1/(2*tr*0.1))} {int(1/(2*tr*0.01))}"
                f" {bd}/{sub}/func/bold_final.nii.gz",
            ],
        }
        return commands

    # ── dMRI 前処理 ────────────────────────────────────────────────────────
    def dmri_preprocessing_commands(self) -> dict[str, list[str]]:
        """
        dMRI (DWI) 前処理コマンド。
        最適パラメータ:
          - eddy: repol=True (外れ値スライス置換), cnr_maps=True
          - b0 正規化: mean_s0_normalization
        """
        sub = self.subject_id
        bd = self.base_dir

        commands = {
            "1_denoising_gibbs": [
                # MP-PCA デノイジング (MRtrix3 dwidenoise)
                f"dwidenoise {bd}/{sub}/dwi/dwi.nii.gz"
                f" {bd}/{sub}/dwi/dwi_denoised.nii.gz"
                f" -noise {bd}/{sub}/dwi/noise_map.nii.gz",
                # Gibbs リンギング補正
                f"mrdegibbs {bd}/{sub}/dwi/dwi_denoised.nii.gz"
                f" {bd}/{sub}/dwi/dwi_degibbs.nii.gz",
            ],
            "2_eddy_correction": [
                # 渦電流・動き補正 (FSL eddy with topup)
                f"eddy_openmp --imain={bd}/{sub}/dwi/dwi_degibbs.nii.gz"
                f" --mask={bd}/{sub}/dwi/nodif_brain_mask.nii.gz"
                f" --index={bd}/{sub}/dwi/index.txt"
                f" --acqp={bd}/{sub}/dwi/acqparams.txt"
                f" --bvecs={bd}/{sub}/dwi/dwi.bvec"
                f" --bvals={bd}/{sub}/dwi/dwi.bval"
                f" --topup={bd}/{sub}/fmap/topup_results"
                f" --repol --cnr_maps"
                f" --out={bd}/{sub}/dwi/dwi_eddy",
            ],
            "3_bias_correction": [
                # N4 バイアス場補正 (ANTs)
                f"dwibiascorrect ants"
                f" {bd}/{sub}/dwi/dwi_eddy.nii.gz"
                f" {bd}/{sub}/dwi/dwi_biascorr.nii.gz"
                f" -fslgrad {bd}/{sub}/dwi/dwi_eddy.bvec {bd}/{sub}/dwi/dwi.bval",
            ],
            "4_registration_MNI": [
                # DWI → MNI152 (via T1)
                f"flirt -in {bd}/{sub}/dwi/nodif_brain.nii.gz"
                f" -ref {bd}/{sub}/anat/T1w_brain.nii.gz"
                f" -omat {bd}/{sub}/reg/DWI_to_T1.mat"
                f" -cost mutualinfo -dof 6",
            ],
        }
        return commands


# ─────────────────────────────────────────────────────────────────────────────
# 2.  前処理品質指標のシミュレーション
# ─────────────────────────────────────────────────────────────────────────────
def simulate_qc_metrics(n_subjects: int = 30) -> dict:
    """
    30 名分の前処理品質指標を生成し、最適パラメータ選定評価を行う。
    """
    rng = np.random.default_rng(42)

    # 動き量 (mm): 健常者の平均値を参考に設定
    mean_fd = rng.normal(0.15, 0.08, n_subjects).clip(0.05, 0.8)
    max_fd = mean_fd + rng.exponential(0.3, n_subjects)

    # DVARS: 時系列の RMS 変動
    dvars = rng.normal(1.2, 0.25, n_subjects).clip(0.6, 2.5)

    # 歪み補正後の幾何学的歪み (mm)
    distortion_before = rng.normal(3.2, 0.8, n_subjects).clip(1.0, 6.0)
    distortion_after = distortion_before * rng.uniform(0.05, 0.15, n_subjects)

    # 空間標準化精度 (Dice 係数, 脳マスク)
    dice_affine = rng.normal(0.92, 0.03, n_subjects).clip(0.80, 0.99)
    dice_nonlinear = rng.normal(0.96, 0.02, n_subjects).clip(0.88, 0.99)

    # tSNR (時間 SNR): 高いほど良い
    tsnr = rng.normal(55, 12, n_subjects).clip(20, 100)

    # 品質管理: FD > 0.5mm またはDVARS > 1.75 を除外推奨
    exclude_flag = (mean_fd > 0.5) | (dvars > 1.75)

    metrics = {
        "n_subjects": n_subjects,
        "n_excluded": int(exclude_flag.sum()),
        "mean_fd_mean": float(mean_fd.mean()),
        "mean_fd_std": float(mean_fd.std()),
        "max_fd_mean": float(max_fd.mean()),
        "dvars_mean": float(dvars.mean()),
        "distortion_reduction_pct": float(
            ((distortion_before - distortion_after) / distortion_before * 100).mean()
        ),
        "dice_affine_mean": float(dice_affine.mean()),
        "dice_nonlinear_mean": float(dice_nonlinear.mean()),
        "tsnr_mean": float(tsnr.mean()),
        "tsnr_std": float(tsnr.std()),
        "subjects_passing_qc": int((~exclude_flag).sum()),
    }

    # 個別被験者データを保存
    per_subject = {
        f"sub-{i+1:02d}": {
            "mean_fd": float(mean_fd[i]),
            "max_fd": float(max_fd[i]),
            "dvars": float(dvars[i]),
            "dice_nonlinear": float(dice_nonlinear[i]),
            "tsnr": float(tsnr[i]),
            "qc_pass": bool(~exclude_flag[i]),
        }
        for i in range(n_subjects)
    }

    return metrics, per_subject


# ─────────────────────────────────────────────────────────────────────────────
# 3.  最適パラメータ比較（グリッドサーチシミュレーション）
# ─────────────────────────────────────────────────────────────────────────────
def parameter_optimization_grid() -> dict:
    """
    fMRI 前処理パラメータのグリッドサーチ結果をシミュレーション。
    評価指標: tSNR, Dice, FD_corr (後続FCとの相関)
    """
    rng = np.random.default_rng(10)

    param_grid = {
        "smoothing_fwhm": [4, 5, 6, 8],
        "hpf_cutoff_hz": [0.008, 0.01, 0.012],
        "lpf_cutoff_hz": [0.08, 0.10, 0.12],
    }

    results = []
    for fwhm in param_grid["smoothing_fwhm"]:
        for hpf in param_grid["hpf_cutoff_hz"]:
            for lpf in param_grid["lpf_cutoff_hz"]:
                # 各パラメータでの性能をモデル化
                tsnr_score = 55 + (fwhm - 5) * 3 + rng.normal(0, 1)
                dice_score = 0.96 - abs(fwhm - 5) * 0.005 + rng.normal(0, 0.002)
                fc_reliability = 0.75 + (fwhm == 5) * 0.05 - abs(hpf - 0.01) * 10 + rng.normal(0, 0.02)
                results.append({
                    "smoothing_fwhm": fwhm,
                    "hpf_cutoff_hz": hpf,
                    "lpf_cutoff_hz": lpf,
                    "tsnr": round(tsnr_score, 2),
                    "dice": round(dice_score, 4),
                    "fc_reliability_icc": round(fc_reliability, 3),
                })

    # 最良パラメータ (fc_reliability で選定)
    best = max(results, key=lambda x: x["fc_reliability_icc"])
    return {"grid_results": results, "best_params": best}


# ─────────────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("[01] 前処理パイプライン実行中...")
    ts = datetime.utcnow().isoformat()

    # パイプラインコマンド生成
    pipeline = PreprocessingPipeline("sub-01")
    fmri_cmds = pipeline.fmri_preprocessing_commands(tr=2.0)
    dmri_cmds = pipeline.dmri_preprocessing_commands()

    # コマンドを保存
    cmd_output = {"fmri_commands": fmri_cmds, "dmri_commands": dmri_cmds}
    with open(f"{RESULTS_DIR}/preprocessing_commands.json", "w") as f:
        json.dump(cmd_output, f, indent=2, ensure_ascii=False)

    # QC 指標シミュレーション
    qc_summary, qc_per_subject = simulate_qc_metrics(n_subjects=30)
    with open(f"{RESULTS_DIR}/qc_summary.json", "w") as f:
        json.dump(qc_summary, f, indent=2, ensure_ascii=False)
    with open(f"{DATA_DIR}/qc_per_subject.json", "w") as f:
        json.dump(qc_per_subject, f, indent=2, ensure_ascii=False)

    # パラメータ最適化
    optim = parameter_optimization_grid()
    with open(f"{RESULTS_DIR}/parameter_optimization.json", "w") as f:
        json.dump(optim, f, indent=2, ensure_ascii=False)

    # ログ
    log_entry = {
        "timestamp": ts,
        "phase": "preprocessing",
        "event_type": "step_completed",
        "actor": "co-scientist",
        "skill_or_tool": "01_preprocessing.py",
        "handoff_out": {
            "qc_summary": qc_summary,
            "best_params": optim["best_params"],
        },
        "files_written": [
            f"{RESULTS_DIR}/preprocessing_commands.json",
            f"{RESULTS_DIR}/qc_summary.json",
            f"{DATA_DIR}/qc_per_subject.json",
            f"{RESULTS_DIR}/parameter_optimization.json",
        ],
        "status": "ok",
    }
    with open(f"{LOGS_DIR}/process-log.jsonl", "a") as f:
        f.write(json.dumps(log_entry, ensure_ascii=False) + "\n")

    print(f"  → QC: {qc_summary['subjects_passing_qc']}/{qc_summary['n_subjects']} 被験者合格")
    print(f"  → 最適パラメータ: {optim['best_params']}")
    print(f"  → 歪み補正効果: {qc_summary['distortion_reduction_pct']:.1f}% 低減")
    return qc_summary, optim


if __name__ == "__main__":
    main()
