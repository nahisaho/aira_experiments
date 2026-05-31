Now let me try NatureLM and GALACTICA MCP tools, and set up the Jupyter notebook:NatureLM and GALACTICA are not available in ToolUniverse. I'll document this. Now let me set up the Jupyter notebook and implement the SV detection pipeline:Now let me create the paper.md:Now let me also copy the figures to the workspace for the markdown image references:タスク完了です。以下にサマリーを示します。

---

## 完了サマリー

### ステップ1: 先行研究調査（ToolUniverse MCP使用）
PubMed・Semantic Scholar で **10件以上**の関連論文を特定：
- Sniffles2ベンチマーク (Cui et al. 2026)、Blackbird (Meleshko et al. 2025)、SqueezeCall (Zhu 2025)、ロングリード希少疾患診断 (Negi et al. 2025)、ハイブリッド解析 (Gambardella 2025, Hu et al. 2025)、ecDNA検出 (Li et al. 2026)など

### ステップ2: NatureLM/GALACTICA接続記録
⚠️ **両ツールともToolUniverseレジストリに未登録** → 接続失敗。Methods §3.4に詳細記録済み（ツール名・エラー内容・代替手段）

### ステップ3: Python実装・実行（Jupyter MCP）
- **2,000件の合成SVデータセット**生成（GIAB HG002統計に基づく14特徴量）
- 5戦略の比較実装（Split-Read / Read-Depth / Assembly / RF統合 / XGBoost Hybrid）
- XGBoost Hybrid: **AUROC=0.9383±0.0040, F1=0.9236±0.0051** [cell:5]
- 複雑なSV: クロモスリプシス AUROC=0.651、ecDNA AUROC=0.840 [cell:10]

### 生成ファイル
| ファイル | 内容 |
|---|---|
| `paper.md` | 学術論文（英語、400行、References 11件） |
| `report.md` | 実験レポート（日本語、270行、図表・数値引用） |
| `figures/*.png` | 5つの可視化図（パイプライン図・ROC・特徴量重要度・包括比較・データ探索） |