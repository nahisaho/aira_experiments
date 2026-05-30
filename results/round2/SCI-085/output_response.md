Now let me search for key foundational papers and run the analysis pipeline:Now I have enough literature context. Let me implement the analysis pipeline:Good. Now let me create the comprehensive analysis pipeline:All 8 figures generated. Now let me create the paper and report:All files are created and the pipeline ran successfully. Here's a complete summary:

---

## 完了サマリー

### ステップ1: 先行研究調査（ToolUniverse MCP）
- **PubMed** + **Semantic Scholar** を使って2020年以降の10件の主要論文を特定
- Norman et al. 2019 (Perturb-seq epistasis), Lotfollahi et al. 2023 (CPA), Barry et al. 2024 (GLM-EIV), Ishikawa et al. 2023 (RENGE), Zhang et al. 2026 (CRISPRi benchmarking) など

### ステップ2: NatureLM 科学的検証
- `ask_naturelm` を3回呼び出し、ガイド検出率(75%)、ドロップアウト率(50%)、DEG数(50-100/摂動) などの定量パラメータを取得し、シミュレーションに組み込んだ

### ステップ3: 実験実施（6モジュール）

| モジュール | 主要結果 |
|---|---|
| QC | 3,473/4,800細胞通過 (72.4%), ガイド検出率 74.6% |
| 差分発現+NMF | 平均 63.0 ± 45.2 DEGs/摂動, 8 遺伝子プログラム |
| 因果GRN | 10ノード, 26エッジ, 密度=0.289 |
| エピスタシス | 3/5 Buffering, 2/5 Additive (r=0.70–0.85) |
| CPA-style表現学習 | **CV R² = −0.129 ± 0.027** (線形モデルの汎化限界を反映) |
| 必須遺伝子 | RUNX1, FLI1, IRF1 ⚠️AUC=1.000は小サンプルアーティファクト |

### 生成ファイル
- **`paper.md`**: 学術論文形式（英語、7セクション、10参考文献、全図埋め込み）
- **`report.md`**: 実験全結果レポート（日本語、全図埋め込み）
- **`figures/`**: 8つのPNG図 + 2つのCSV結果ファイル