NatureLM and GALACTICA MCPs are not available - will record this. Let me now search for papers in parallel.Good, I have sufficient literature. Let me now try NatureLM/GALACTICA tools and set up Jupyter simultaneously.NatureLM and GALACTICA are not available in ToolUniverse (confirmed: 0 results). Let me now proceed to implement the Python pipeline in Jupyter.All computations complete. Now let me write both output files.すべてのステップが完了しました。以下に成果をまとめます。

---

## 完了報告

### ステップ1: 先行研究調査（Crossref / Semantic Scholar）

5件以上の論文を特定（2020年以降）:
| # | 著者・年 | DOI | 主要知見 |
|---|---------|-----|---------|
| 1 | Pérez-Hernández & Benet (2022) | 10.1038/s43247-021-00337-x | アポフィスYarkovsky漂移 −2.923×10⁻³ AU/Myr 検出 |
| 2 | Greenberg et al. (2020) | 10.3847/1538-3881/ab62a3 | 247個NEAのYarkovsky漂移検出 |
| 3 | Liu et al. (2023) | 10.3847/1538-4357/accc81 | NEA軌道漂移の測光・位置天文融合測定 |
| 4 | Domínguez et al. (2023) | 10.1016/j.actaastro.2022.10.039 | 短警告時間小惑星偏向の動力学解析 |
| 5 | Zhao et al. (2025) | 10.1360/tb-2025-0041 | 2024 RW1リアルタイム衝突予測 |

### ステップ2: NatureLM / GALACTICA

- **NatureLM** (`ask_naturelm`): ToolUniverseに未登録（試行記録済み）
- **GALACTICA** (`scientific_qa`, `predict_citations`): ToolUniverseに未登録（試行記録済み）
- 代替：Crossref + Semantic Scholar による文献パラメータ検証

### ステップ3: Python実装・実行（Jupyter MCP）

12セルを実行、主要結果:
- **Monte Carlo**: 10,000クローン生成、MOID = 0.01325 AU `[cell:9]`
- **KS検定**: p = 0.9992（統計的妥当性確認）`[cell:9]`
- **衝突エネルギー**: アポフィス型 = 2,591 MT TNT, 爆風半径 54.9 km `[cell:5]`
- **DART偏向 Δv**: 0.211 ± 0.069 mm/s（10年警告、軌道偏差 0.021 R⊕）`[cell:6]`
- **ヤルコフスキー累積**: 100年で −3,740,000 km `[cell:9]`

### 成果物

| ファイル | サイズ |
|--------|-------|
| `paper.md` | 441行（英語学術論文） |
| `report.md` | 246行（日本語実験レポート） |
| `figures/fig1_neo_mc_uncertainty.png` | 491 KB（6パネル） |
| `figures/fig2_risk_deflection.png` | 225 KB（3パネル） |
| `figures/fig3_advanced_analysis.png` | 427 KB（3パネル） |
| `data/raw/*.csv, *.json` | 数値データ4ファイル |