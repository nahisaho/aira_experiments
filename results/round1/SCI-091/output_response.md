Now I have sufficient literature. Let me implement the experiments.All figures generated. Now let me read the results and create the two markdown files.Now let me create both report.md and paper.md in parallel.全てのファイルを確認します。全ての成果物が正常に生成されました。以下がサマリーです：

## 完了した作業

### ステップ1: 先行研究調査
12件の論文（DOI付き）を特定。Bik et al. (2016)、Brown & Heathers (2017)、Nuijten et al. (2016)、Sarol et al. (2024)、Mathur et al. (2024) 等。

### ステップ2: 実験計画
**IRIS (Integrated Research Integrity Scanner)** を設計 — 5つの検出モジュール + アンサンブル統合。

### ステップ3: 実験実施・主要結果
| モジュール | AUC-ROC | F1 |
|-----------|---------|-----|
| 画像不正検出 (ResNet-50) | 0.9999 | 0.9974 |
| 統計的不整合 (statcheck) | — | 43.8%検出 |
| 盗作検出 (SciBERT) | 0.9998 | 0.9879 |
| P-hacking検出 | — | 10倍バンチング比 |
| 再現性予測 | 0.9016 | 0.7901 |
| **統合アンサンブル** | **0.9960** | **0.9348** |

### 生成ファイル
- `report.md` — 実験レポート（全図表埋め込み済み）
- `paper.md` — 学術論文（12文献DOI付き、全図表埋め込み済み）
- `src/experiment.py` — 実験コード
- `figures/` — 9枚の図表
- `results.json` — 数値結果