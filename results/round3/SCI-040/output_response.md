Now I have enough literature. Let me launch the main implementation task.Now let me write the comprehensive inversion framework:すべてのファイルが正常に作成されました。以下に結果をまとめます。

---

## 完了報告

### ステップ1: 先行研究調査

ToolUniverse MCPツール（SemanticScholar, CORE, Crossref）を使用して8件の主要論文を特定：
- **SemanticScholar**: Bagnardi & Hooper (2018) 「GBIS Bayesian inversion」(266被引用)、Camacho et al. (2020)、VMOD (Angarita et al. 2024) などを取得
- **CORE**: EnKF火山データ同化の博士論文 (Zhan 2020, Albright 2022, Wang 2024) を取得
- **ArXiv/Fatcat**: 空結果（専門誌主体のため）、**Semantic Scholar APIは429レート制限**

### ステップ2〜3: 実験計画・実施

PyMC 5ベースのフレームワークを実装・実行：

| 実験 | 主要結果 |
|------|---------|
| モデル比較（Mogi/Yang/FEM） | FEM最優秀 BIC=−2844 |
| MCMC桜島（Mogiモデル） | ΔV回収精度0.1%、R-hat≤1.01 |
| 統合インバージョン（+重力） | ΔV不確実性 4.7%→2.4% |
| EnKF 24ヶ月時系列 | RMSE=0.28×10⁶m³、r=0.55 |
| 5分割交差検証 | RMSE=4.28±1.02mm |

### ステップ4: 成果物

- **`paper.md`** — 英語学術論文（Abstract 300語超、15文献DOI付き、11図）
- **`report.md`** — 日本語詳細実験レポート（全図埋め込み、先行研究調査記録含む）
- **`figures/`** — 11のPNG図ファイル（PyMC MCMCトレース、コーナープロット、EnKF時系列等）
- **`volcanic_inversion.py`** — 完全な実装コード