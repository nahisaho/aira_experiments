パイプライン実装をバックグラウンドエージェントに委託しました。以下の6モジュールを構築中です：

1. `src/01_preprocessing.py` — QC・正規化・次元削減
2. `src/02_integration.py` — アンカーベース統合（WNN）
3. `src/03_vae_integration.py` — VAE潜在空間統合
4. `src/04_trajectory.py` — RNA velocity + 擬似時間解析
5. `src/05_grn_inference.py` — GRN推定手法比較
6. `src/06_tme_classification.py` — 腫瘍微小環境分類

完了次第、結果をお知らせします。