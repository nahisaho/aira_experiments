I have gathered enough literature. Now let me implement the GNN weather prediction model and run the experiments.Now I'll create both output documents:すべてのファイルが正常に生成されました。以下に全ステップの成果をまとめます。

---

## 実験完了サマリー

### ステップ1: 先行研究調査
OpenAlex / Fatcat MCPツールで以下の10件の重要論文を特定しました：
- **GraphCast** (Lam et al., Science 2023)、**Pangu-Weather** (Bi et al., Nature 2023)、**FourCastNet** (Kurth et al., 2023)、**FuXi** (Chen et al., 2023)、**FengWu** (Chen et al., 2023) 等

### ステップ2: NatureLM検証
- 3回クエリ実施、接続成功したが定量的RMSEベンチマーク値の精度は不十分
- 得られた定性的知見（質量保存制約・エネルギー比）を物理メトリクス設計に活用

### ステップ3: 実験結果（PyTorch Geometric GNN, 65K params）

| 指標 | 6h | 24h | 120h |
|------|-----|-----|------|
| T RMSE (medium) | **1.90±0.07 K** | 3.58±0.39 K | 13.4±4.1 K |
| Z500 RMSE | **23.0±0.7 m²/s²** | 31.2±1.4 | 133±50 |
| Z500 Skill | **+0.20** | +0.08 | −0.81 |

### ステップ4: 成果物
- ������ **`paper.md`** — 英語学術論文（Abstract 300語+、References 10件、図8枚埋め込み）
- ������ **`report.md`** — 日本語実験レポート（全結果・考察・図8枚埋め込み）
- 📊 **`figures/`** — 8枚の可視化図（訓練曲線、RMSE比較、スキルスコア、垂直プロファイル、予測マップ、誤差分布、エネルギースペクトル、結果テーブル）