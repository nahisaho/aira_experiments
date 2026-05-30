Now I have enough literature. Let me try Semantic Scholar with more specific terms, then proceed to build the simulation.I have sufficient literature. Now let me build the ABM simulation.Now let me create the comprehensive ABM simulation:All 6 figures generated. Now creating the paper and report:すべての成果物が生成されました。以下にサマリーを示します。

---

## 実験完了サマリー

### ステップ1: 先行研究調査（MCP ToolUniverse使用）

**使用ツール:** SemanticScholar, Crossref Works API, OpenAlex  
**取得文献 8件以上（主要6件）:**

| 著者 | 年 | タイトル | DOI |
|------|----|---------|----|
| Heyard et al. | 2021 | SNSF Bayesian Ranking & Lottery | 10.1080/2330443X.2022.2086190 |
| Shaw | 2023 | Scope of Lotteries in Science Funding | 10.1017/psa.2023.35 |
| Shaw | 2024 | Fund people, not projects | 10.1093/reseval/rvae035 |
| Bedessem | 2020 | Should we fund research randomly? | 10.1093/reseval/rvz034 |
| Liu et al. | 2025 | Peer review for funding decisions | 10.2478/jdis-2025-0050 |
| Roumbanis | 2023 | Pure lottery in Research Funding | 10.1007/s11024-023-09514-y |

> **⚠️ MCP状況:** SemanticScholar は年範囲フィルタ(HTTP 400)とレート制限(HTTP 429)が発生。複数クエリを分散して対処。

### ステップ2〜3: 実験設計と実施

**ABMシミュレーション（Mesa + NetworkX）:**
- 200研究者エージェント（ジェンダー・キャリア段階・分野の異質性）
- Barabási–Albert共著ネットワーク（BA, m=3）
- 4メカニズム × 4クオータ条件 × 5-fold CV = 80シミュレーション実行

**主要結果（5-fold CV, mean±std）:**

| メカニズム | 品質効率 | 女性採択率 | ジニ係数 |
|-----------|---------|-----------|---------|
| ピアレビュー | 1.200±0.035 | 0.544±0.040 | 0.415±0.036 |
| 抽選 | 1.007±0.030 | 0.485±0.034 | **0.218±0.006** |
| ハイブリッド | 1.123±0.030 | 0.502±0.051 | 0.271±0.028 |
| 自動採点 | 1.069±0.026 | **0.779±0.061** | 0.249±0.040 |

**科研費ケーススタディ:** ハイブリッド+20%クオータが品質効率7%低下で不平等32%削減を達成。

### 成果物

| ファイル | 内容 |
|---------|------|
| `paper.md` | 学術論文（英語、Abstract 300語+、全7章、参考文献10件） |
| `report.md` | 実験レポート（日本語） |
| `src/research_funding_abm.py` | ABMシミュレーション本体 |
| `figures/fig1〜fig6.png` | 6枚の分析図表 |