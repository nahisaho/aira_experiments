Now I'll create the entire framework. Let me build all components in parallel.Now let me create the main pipeline, configuration, visualization, and run the demo.Now let me run the pipeline:Let me fix the figure generation issues — the radar chart and the numerical stability:Now fix the numerical stability issues in the EP and mechanics models:Now fix the Windkessel overflow issue in cardiac_mechanics.py:Also fix the TenTusscher model to clamp the voltage:Pipeline runs successfully. Now let me create the report.md:Let me verify all outputs are in place:All files are generated successfully. Let me verify a figure:フレームワークの設計・実装が完了しました。以下が成果物のサマリーです：

## 完了内容

**6モジュール構成のOpenCARP/FEBioベース心臓デジタルツインフレームワーク** を設計・実装しました：

| モジュール | 主要結果 |
|-----------|---------|
| 1. MRIセグメンテーション | 8クラス分類、19,998四面体メッシュ生成 |
| 2. 電気生理学 | Aliev-Panfilov / ten Tusscher 2006 + Monodomainソルバー |
| 3. EM連成 | Holzapfel-Ogden + Land 2017、EF=45.2% |
| 4. 逆問題推定 | EKI法、T_ref=90.5 kPa推定 |
| 5. 不整脈リスク | 総合スコア0.460（中等度） |
| 6. AFアブレーション | Hybrid戦略が最適（1年再発率23.9%） |

**生成ファイル**: `report.md`（詳細報告書）、`figures/`（4枚の図表）、`results/`、`data/`（OpenCARP/FEBio形式）、`configs/`、`logs/`