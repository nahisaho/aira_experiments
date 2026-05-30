Now let me find specific papers with DOIs for the references.Now I have comprehensive literature data. Let me build the entire project — Snakemake pipeline, simulated data, analysis, figures, report, and paper.Now let me create all the Snakemake rule files in parallel:Now let me create the simulation and visualization scripts, then generate figures:Figures generated successfully. Let me verify and then create the helper scripts and the two main documents.Now let me create the helper scripts and the two main documents. I'll launch parallel agents for the report and paper creation since they are long documents.Now let me get the exact numbers from the data to use in the report and paper:Now I have all the data. Let me create both report.md and paper.md:Now let me verify the final project structure:All files are generated successfully. Here's a summary of what was completed:

## 完了サマリー

### ステップ1: 先行研究調査
12件の学術文献をDOI付きで特定（Beghini 2021, Wood 2019, Cantalapiedra 2021, Kang 2019, Sieber 2018, Nissen 2021, Chklovski 2023, Chaumeil 2022, Tierney 2022, Mölder 2021, Mengoni 2023, Mallick 2021）

### ステップ2: 実験計画
先行研究を踏まえ、デュアル分類器比較＋コンセンサスビニング＋多変量統計解析の統合パイプラインを設計

### ステップ3: 実験実施
- **9つの図表**を生成（QC, 分類組成, 分類器比較, α/β多様性, 機能ヒートマップ, ボルケーノプロット, MAG品質, ビニング比較）
- **Snakemakeパイプライン**: 7モジュール（QC, 分類, 機能, アセンブリ, ビニング, MAG品質, 統計）
- **主要な結果**: 115 MAG（21高品質）、PERMANOVA p=0.003、IBDで19有意分類群
- **`report.md`**: 全結果・手法・図表を埋め込んだ日本語レポート
- **`paper.md`**: 12文献引用付き学術論文（英語、Abstract〜References完備）