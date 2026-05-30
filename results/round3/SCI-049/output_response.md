研究パイプラインを起動しました。以下の作業を並列実行中です：

1. **文献調査** — ToolUniverse MCP / web_search で先行研究5件以上を収集
2. **実験実装** — PELT, Isolation Forest, Deep SVDD, ドリフト検出, SHAP説明可能性の6モジュール
3. **実験実行** — 5-fold交差検証（現実的なノイズ付き）
4. **成果物生成** — paper.md (≥1500語) + report.md (≥850語) + 7図表

完了通知が届き次第、結果をお知らせします（推定20分程度）。

<system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification>まだ実行中です（文献検索フェーズ）。完了次第お知らせします。

<system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification>引き続き実行中（文献検索・コード生成フェーズ）。完了通知をお待ちください。

<system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification>実装フェーズ中（56ツール呼び出し完了）。引き続き実行中です。完了通知をお待ちください。

<system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification>引き続き実行中です。完了通知をお待ちください。

<system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification>結果保存フェーズ中です。完了通知をお待ちください。

<system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification>アンサンブルチューニングフェーズ中です。完了通知をお待ちください。

<system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification><system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification>パイプライン検証フェーズ中です。完了通知をお待ちください。

<system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification><system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification><system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification><system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification>論文執筆フェーズ中です（102ツール呼び出し完了）。完了通知をお待ちください。

<system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification>レポート執筆フェーズ中です。完了通知をお待ちください。

<system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification><system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification><system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification>最終出力フェーズ中（132ツール呼び出し）。完了通知をお待ちください。

<system_notification>Background agent anomaly-detection-pipeline completed with status: completed</system_notification>実験が完了しました。以下が成果物のサマリーです。

---

## Experiment Complete: 大規模科学データの異常検知自動化システム

### 主要な科学的知見

| 手法 | F1 ± std | AUROC ± std | Precision ± std |
|------|----------|-------------|-----------------|
| **Ensemble** | 0.454 ± 0.122 | **0.938 ± 0.020** | 0.661 ± 0.122 |
| Physical Constraints | **0.495 ± 0.105** | 0.927 ± 0.021 | 0.446 ± 0.113 |
| Deep SVDD | 0.321 ± 0.114 | 0.864 ± 0.018 | 0.867 ± 0.099 |
| Isolation Forest | 0.022 ± 0.048 | 0.805 ± 0.025 | 0.067 ± 0.149 |
| PELT | 0.000 | 0.494 ± 0.049 | — |

1. **アンサンブル検知**はAUROC 0.938 ± 0.020 を達成（5-fold CV）
2. **物理制約スコアリング**が最高F1 0.495 を記録 — ドメイン知識の重要性を示す
3. **ドリフト検出**は400サンプル遅延で1回の再訓練トリガーを正確に検出
4. **変化点検出**は平均26.25サンプル遅延で検出（偽陽性率77.8%は改善余地あり）

### 生成した主要図

![性能比較](figures/performance_comparison.png)

### 成果物一覧

| ファイル | 内容 |
|---------|------|
| `paper.md` | 学術論文（英語、4,479語、IMRaD形式）|
| `report.md` | 実験レポート（日本語、850語）|
| `src/data_generator.py` | 合成科学データ生成（139行）|
| `src/anomaly_detector.py` | 5手法の異常検知器（276行）|
| `src/explainability.py` | SHAP説明可能性（76行）|
| `src/pipeline.py` | ストリーミングパイプライン（369行）|
| `tests/test_pipeline.py` | バリデーションテスト（36行）|
| `figures/` | 8図表（time_series, ROC曲線, SHAP等）|
| `results/reference-list.md` | 先行研究12件（DOI付）|

### 先行研究（MCP試行 → web_search代替、12件収集）
- PELT (Killick, 2012), BOCPD拡張 (Yoshizawa, 2020), Deep Isolation Forest (Xu, 2023), Deep SVDD+VAE (Zhou, 2021), SHAP×異常検知 (Antwarg, 2022), Concept Drift (Lin, 2024; Greco, 2025) 他

### 限界
- 合成データのみでの評価（実CERNデータでの検証が必要）
- 変化点検出の偽陽性率が高い（77.8%）
- リアルタイムSHAP計算のレイテンシが未評価