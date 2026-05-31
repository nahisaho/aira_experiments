NatureLM and GALACTICA are not available. I'll document this as required. Now let me search for literature and set up the Jupyter notebook.The notebook was created. Let me connect correctly and continue searching for more papers.The Jupyter server is for a different workspace. I'll run Python directly in bash for my workspace.Good, I have enough papers. Now let me run the main Python simulation analysis.I see a PV calculation bug. Let me fix it and re-run.All simulations completed successfully. Now let me also search for a few more papers and then write the final documents.All figures generated. Now let me create the paper.md and report.md.全ファイルが正常に生成されました。タスク完了です。

---

## 完了サマリー

### ステップ1: 先行研究調査
Semantic Scholar MCP (`SemanticScholar_search_papers`) で5つのクエリを実行し、**10件の論文**を特定：
- BIM-BPS変換 (Alexandrou 2025, Xu 2020, Nasyrov 2014)
- EnergyPlus自動化 (Osei-Owusu 2025, Westermann 2020)
- CFD換気 (Tai 2022, Li 2025)
- 昼光シミュレーション (Tong 2023, Mangkuto 2025, Abedini 2025)

### ステップ2: NatureLM/GALACTICA MCP
両ツール (`ask_naturelm`, `scientific_qa`) は **ToolUniverse未収録** → Methods セクションに試行記録・代替手段を明記

### ステップ3: Python実装・実行
`bim_simulation.py` を実行（全 Cell 0〜8 相当）。主要結果：
- **EUI** = 295.3 kWh/m²/yr (基準ケース), PVオフセット率 62% [cell:2]
- **ACH** = 10.38 (スタック効果), 40.02 at 3m/s [cell:3]
- **sDA₃₀₀** = 53.5%, **ASE₁₀₀₀** = 1.0% [cell:4]
- **RF R²** = 0.871±0.022, **GB R²** = 0.873±0.022 [cell:5]
- PV容量・窓U値が最重要設計変数（重要度0.61, 0.31）[cell:5]

### 成果物
| ファイル | 説明 |
|---------|------|
| **`paper.md`** | 学術論文（英語、Abstract 250語超、全セクション、DOI付き参考文献10件、[cell:N]引用） |
| **`report.md`** | 実験レポート（日本語、図表埋め込み） |
| `figures/fig1〜5.png` | 5件の図（エネルギー収支、CFD、昼光、ML、ダッシュボード） |
| `data/raw/` | 全生データ・pip_freeze |
| `bim_simulation.py` | 完全再現可能なPythonコード |