Now let me run literature searches in parallel and set up the Jupyter environment.NatureLM and GALACTICA are not available in ToolUniverse (I'll document this). Now let me build the full ABM simulation notebook.Now I have all the data. Let me create the paper.md and report.md files.The figures are in the Jupyter server's filesystem, not the workspace. Let me copy them over.すべてのファイルが揃いました。実験完了です。以下に成果をまとめます：

---

## 完了サマリー

### ステップ1: 先行研究調査（Crossref API使用）
Semantic Scholar（429レート制限）とCrossref APIを使って7件の先行研究を特定：
- Philipps (2021): 科研費抽選配分への科学者の意識調査
- Shaw (2022): 抽選型資金配分の系統的レビュー  
- Roshani et al. (2021): 研究資金と引用成果の関係分析
- Matveeva et al. (2026): 共著ネットワーク分析
- González-Salmón et al. (2026): ジェンダー・資金の不平等三角形
- Gundur & Kumar (2025): 計量書誌学におけるジェンダー不平等

### ステップ2: NatureLM / GALACTICA MCP接続試行
両ツールともToolUniverseに存在せず（0件マッチ）。エラー内容・代替手段をMethodsセクションに記録。

### ステップ3: Jupyter ABMシミュレーション実行
- 200エージェント × 20年 × 4メカニズムのシミュレーション完了
- RF AUROC = **0.811 ± 0.103** [cell:7]
- Kruskal-Wallis H=3.36, p=0.339（4メカニズム間に有意差なし）[cell:10]
- 科研費ジェンダーバイアス: -3.2pp効果 [cell:6]
- パレート最適点: h_weight=0.8, div_weight=0.2 [cell:8]

### 成果物
| ファイル | サイズ |
|---------|-------|
| `paper.md` | 25KB（学術論文形式） |
| `report.md` | 12KB（実験レポート） |
| `figures/abm_results_part1.png` | 388KB |
| `figures/abm_results_part2.png` | 185KB |
| `figures/abm_results_part3.png` | 156KB |
| `data/raw/researchers.csv` | 13KB |
| `research_funding_abm.ipynb` | Jupyter notebook |