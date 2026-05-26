Now I have sufficient literature. Let me proceed with the experiment implementation.Good. Now let me create the main experiment script.

Let me fix the error in the scenario comparison plot and the metrics issues.The issue is that `pred_temp` starts with an array (from `mean(axis=(-2,-1))` on seq_len frames) then appends scalars. Let me fix this and also the NaN pattern correlation for SLR.Now fix the NaN pattern correlation for sea level (constant spatial field):Now re-run:Experiment completed successfully. Now let me verify the figures were created, then create the report and paper.All 6 figures and metrics.json generated. Now let me read the metrics and create both report.md and paper.md.Now let me create both files in parallel:All files created. Let me verify the final workspace structure:実験が完了しました。以下のファイルを生成しました：

**成果物:**
- `report.md` — 実験レポート（手法・結果・考察・全図表埋め込み）
- `paper.md` — 学術論文形式（Abstract〜References、11文献引用、全6図埋め込み）
- `src/experiment.py` — 実験コード（データ生成・U-Net/ConvLSTM・物理制約学習・評価）

**生成図表（6点）:**
- `training_curves.png` — 学習曲線
- `spatial_predictions.png` — 空間場予測 vs 真値
- `scenario_comparison.png` — SSPシナリオ別温度比較
- `ensemble_uncertainty.png` — アンサンブル不確実性
- `metrics_comparison.png` — モデル間評価指標
- `physics_constraints.png` — 物理制約遵守評価

**主要知見:** ConvLSTMが気温予測で優位（RMSE 10.77K）、U-Netが降水パターン再現で優位（相関0.652）。物理制約は降水非負性を効果的に強制。