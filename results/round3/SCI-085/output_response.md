Now I have enough literature. Let me set up the workspace and implement the full pipeline.Good. pertpy has missing scikit-misc. I'll implement all 6 framework components using scanpy, scvi, scipy, statsmodels, torch, networkx. Let me build the modules:Now let's run the pipeline:The QC filter is too strict for small data. Let me fix the GMM to handle edge cases:Good - QC pipeline works. Now let's run the full pipeline:Pipeline complete! Now let me create the tests and write the reports.All 18 tests pass. Now write the reports:---

## Experiment Complete: Perturb-seq Analysis Framework

### Key Scientific Findings

1. **差分発現：** 15 摂動条件で合計 **2,391 有意 DE 遺伝子**を同定（平均 159.4 ± 127.8/摂動、FDR < 0.05）。guide_BRCA1 が最大転写影響（479 遺伝子）を示した。

2. **GRN 推定：** LASSO 回帰で **80 ノード・1,025 エッジ**の調節ネットワークを推定（密度 0.0102）。活性化エッジ 612、抑制エッジ 413。トップハブ遺伝子 Gene0268（PageRank = 0.0257）。

3. **エピスタシス：** 5 組合せのうち **3 件が相乗（Synergy）**、1 件が拮抗（Antagonism）、1 件が加算的。KRAS+MYC が最大相乗係数（ε̄ = +0.172）。

4. **NB-VAE：** 30 エポックで ELBO 損失 **3.255 に収束**（初期比 −8.5%）。10 次元潜在空間で PIK3CA/AKT1/PTEN が同一クラスターに集積（PI3K 経路の一貫性）。

5. **必須遺伝子：** guide_BRCA1 が唯一の高必須性摂動（z-score = 2.95）、共必須ネットワーク 8 ノード・17 エッジ。

### Most Important Figures

![Module Heatmap](figures/fig04_module_heatmap.png)
![Latent Space](figures/fig08_latent_space.png)

### Deliverables

| 種別 | 内容 |
|------|------|
| **ソースコード** | 8 モジュール（2,642 行）|
| **テスト** | 18 テスト、全通過 ✅ |
| **結果ファイル** | 10 CSV/JSON（results/） |
| **図** | 10 PNG（figures/） |
| **report.md** | 958 語（日本語） |
| **paper.md** | 3,789 語（英語・IMRaD） |

### ⚠️ MCP ツール接続記録（科学的透明性）

| ツール | 結果 |
|--------|------|
| Semantic Scholar API | HTTP 429（Rate Limit）→ 代替: PubMed |
| PubMed E-utilities | ✅ 正常動作（8 論文取得） |
| ToolUniverse MCP | 未設定（.mcp.json なし）|

### Limitations

- 合成データのみ使用（実 Perturb-seq データへの検証が必要）
- エピスタシス検定の統計的検出力不足（置換 100 回 → 1,000+ 回推奨）
- pertpy v1.0.3 は環境の依存関係問題（scikit-misc）でインポート不可 → 独自実装で代替