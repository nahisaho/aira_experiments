Now I'll create the full pipeline. Let me build all components in parallel.Now I'll create all pipeline modules in parallel:Now let me create the main orchestrator script and the report:Now creating the comprehensive report:最後にログを更新し、ファイル構成を確認します。空間トランスクリプトミクス解析パイプラインの設計が完了しました。

## 生成物サマリ

| カテゴリ | ファイル |
|---------|---------|
| **設定** | `config.yaml`, `requirements.txt` |
| **オーケストレータ** | `run_pipeline.py` |
| **M0** データ前処理 | `pipeline/m00_data_loading.py` — QC・正規化・クラスタリング・空間グラフ |
| **M1** デコンボリューション | `pipeline/m01_deconvolution.py` — cell2location 2段階推論 |
| **M2** SVG 検出 | `pipeline/m02_spatial_patterns.py` — Moran's I + SpatialDE コンセンサス |
| **M3** 細胞間通信 | `pipeline/m03_communication.py` — LIANA 6手法統合 + Squidpy 空間共局在 |
| **M4** ニッチ同定 | `pipeline/m04_niche.py` — 近傍プロファイル + Leiden/KMeans |
| **M5** 3D 再構成 | `pipeline/m05_reconstruction_3d.py` — ICP アラインメント + z-stacking |
| **M6** 腫瘍免疫解析 | `pipeline/m06_tumor_immune.py` — 境界検出・浸潤勾配・疲弊スコア |
| **レポート** | `report.md` — 全手法・結果・考察・ファイル一覧 |

実行: `python run_pipeline.py --config config.yaml`（Visium データを `data/raw/` に配置後）