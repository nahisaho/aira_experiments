Now let me search for specific DOIs for key papers and then proceed with the experiment implementation.Now let me get a few more DOIs and then start the implementation.Good, I now have comprehensive literature data. Let me implement the experiment.Now let me create the main experiment script:

The ITER dataset is too small for the batch size. Let me fix this.

Experiments completed successfully. Now let me create the report.md and paper.md files.Verify all files exist:全ファイルが正常に生成されました。以下が実験の要約です：

---

**実験完了サマリー**

**先行研究**: 8件の主要論文を特定（全DOI付き）— Vega 2022 (Nature Physics), Churchill 2020, Rea 2020, Zheng 2023, Seo 2024 (Nature), Kim 2024, Sabbagh 2023, Conlin 2021

**モデル性能（JETテストセット）**:
| モデル | F1 | AUC | P99レイテンシ |
|--------|-----|-----|------------|
| LSTM | 0.989 | 1.000 | 4.24ms ✓ |
| CNN-LSTM | 0.986 | 1.000 | 10.37ms ✓ |
| PI-LSTM | 0.989 | 1.000 | 5.92ms ✓ |
| TM-MultiTask | 0.989 | 1.000 | 17.31ms ✓ |

**転移学習**: JET→KSTAR(F1=0.987) > KSTARネイティブ(F1=0.968)

**生成物**: `report.md`, `paper.md`, 10枚の図（`figures/`）, 3つのソースファイル（`src/`）, `results.json`