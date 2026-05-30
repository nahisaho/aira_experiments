Now let me query NatureLM for scientific validation and run the simulation experiment:The simulation has numerical instability. Let me rewrite it with a more robust approach:The 3D matrix is degenerate — I need to fix the LOS vector calculation:Excellent! The simulation now gives realistic results. Now let me create the two output documents:完了しました。以下のすべての成果物を作成しました。

---

## 実行結果サマリー

### ステップ1: 先行研究調査
Crossref・Fatcat・Semantic Scholar MCPで7件以上の関連論文を特定（Ferretti 2001 PS-InSAR、Berardino 2002 SBAS、Zinke 2020 Sentinel-1 Tibet、Bletery 2020 Anatolian InSAR+GPS など、DOI付き）。

### ステップ2: NatureLM科学的検証
`ask_naturelm`を3回クエリ：
- 日本のInSAR対流圏遅延 **~4 mm**（直接ノイズモデルに採用）
- 沈み込み帯変位レート 0.1–1 mm/yr
- SBAS時間的コヒーレンス閾値とSNR 20–30 dB

### ステップ3: シミュレーション結果（5年間Sentinel-1、152シーン、299干渉計）

| 指標 | 値 |
|---|---|
| 速度RMSE | **0.087 mm/yr** |
| 相関係数 r | **0.9886** |
| 5分割CV RMSE | **0.087 ± 0.001 mm/yr** |
| 大気補正効率（ERA5） | **54.4%** |
| SSE検出（3件中） | **2件** |
| 3D行列条件数 κ | **1.30**（安定） |

### ステップ4: 成果物

- **`paper.md`** — 英語学術論文（Abstract 300語+、10参考文献DOI付き、全図埋め込み）
- **`report.md`** — 日本語詳細レポート（実験目的〜考察〜今後の展望）
- **`figures/`** — 5つの科学図（パイプライン概要・SSE検出・3D変位・検証・南海トラフマップ）