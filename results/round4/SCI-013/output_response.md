ToolUniverse の学術検索ツールが利用可能です。並行して複数のキーワードで検索を実行します。十分な文献情報が得られました。次に実験実装フェーズに移ります。CSPが完璧な精度を出すのはSNRが高すぎるため。データ生成をより現実的にして再実行します。実験完了。次に report.md と paper.md を作成します。すべての成果物が生成されました。以下が結果サマリーです。

---

## 実行結果サマリー

### ステップ1: 先行研究調査（MCP ToolUniverse）

**SemanticScholar**: HTTP 400/429エラー（rate limit）  
**OpenAlex API**: ✅ 成功 — 8件の関連論文を特定

| # | 論文 | 年 | 引用数 | DOI |
|---|------|-----|-------|-----|
| 1 | EEG-BCI Comprehensive Review (Rashid et al.) | 2020 | 479 | 10.3389/fnbot.2020.00025 |
| 2 | BCIs Survey on CI Approaches (Gu et al.) | 2021 | 376 | 10.1109/tcbb.2021.3052811 |
| 3 | P300 Transfer Learning XDAWN+Riemannian (Li et al.) | 2020 | 58 | 10.3390/app10051804 |
| 4 | STNN for P300 Detection (Zhang et al.) | 2021 | 12 | 10.1109/access.2021.3132024 |
| 5 | TCFormer for MI Decoding (Altaheri et al.) | 2025 | 6 | 10.1038/s41598-025-16219-7 |
| 6 | CSP-Net (Jiang et al.) | 2024 | 26 | 10.1016/j.knosys.2024.112668 |
| 7 | EEG is Better Left Alone (Delorme) | 2023 | 224 | 10.1038/s41598-023-27528-0 |
| 8 | Locked-in Spelling via Neurofeedback (Chaudhary et al.) | 2022 | 139 | 10.1038/s41467-022-28859-8 |

### ステップ3: 実験結果（5分割交差検証）

| モデル | Accuracy (mean±SD) | Macro F1 |
|-------|-------------------|---------|
| CSP+LDA | 1.000±0.000 ⚠️ | 1.000±0.000 |
| CSP+SVM | 1.000±0.000 ⚠️ | 1.000±0.000 |
| CSP+NN | 0.990±0.020 | 0.990±0.020 |
| **EEGNet** | **0.582±0.092** | 0.583±0.087 |
| **EEG Conformer** | **0.967±0.027** | 0.967±0.028 |

> ⚠️ CSP 系の完璧スコアは**合成データの構造的単純さ**（クラス別周波数成分が明確分離）に起因。実データでのCSP+LDA文献値は約72–78%。詳細は paper.md §6.1 参照。

### ステップ4: 成果物

- **`paper.md`** — 学術論文形式（Abstract 300語以上、8セクション、文献8件）
- **`report.md`** — 実験レポート（図10枚すべて埋め込み済み）
- **`figures/`** — 図10枚（パイプライン図・ERP・混同行列・転移学習比較など）