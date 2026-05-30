Now let me search for papers on TCR repertoire analysis in parallel:I now have sufficient literature. Let me implement the full pipeline and generate results:Now let me write the full analysis pipeline:Excellent! All results are in. Now let me create both output files:全ての成果物が作成されました。以下に結果をまとめます。

---

## 完了サマリー

### Step 1: 先行研究調査（ToolUniverse MCP使用）

**使用ツール:** SemanticScholar_search_papers（429エラーあり→部分成功）、Crossref_search_works（成功）、openalex_literature_search（全クエリ成功）

**特定した主要論文（10件以上、2020年以降）:**
- Krishna et al. 2020 — HLA多型とTCR多様性 (Immunity & Ageing)
- Sethna et al. 2020 — TCRレパトア生成・選択の集団変動 (PLoS Comp. Biol.)
- Song et al. 2021 (TRUST4) — RNA-seqからのTCR再構成 (Nature Methods)
- Zhang et al. 2021 (TCRAI) — TCR-抗原特異性予測 (Science Advances)
- Jiang et al. 2022 (TEINet) — TCR-epitope結合予測DL (AUROC=0.760)
- Luo et al. 2023 (NEO2IS) — ICB有効性統合予測 (Oncogene)
- Kim et al. 2023 (TSpred) — アンサンブルCNN+注意機構
- Lupyr et al. 2025 (TCRgrapher) — ALICE/TCRNET/GLIPH2/tcrdist3比較

### Step 2–3: 実験結果（80名コホート、145,015クローン）

| モジュール | 最良手法 | 性能指標 |
|---|---|---|
| TCR-Epitope結合予測 | Gradient Boosting | AUROC=**0.6458±0.0035**, F1=**0.685** |
| 免疫年齢推定 | Random Forest | MAE=**7.76年**, R²=**0.739** |
| ICB応答予測 | Random Forest | AUROC=**0.616±0.121** |
| 多様性-年齢相関 | Shannon entropy | r=**−0.72** (p<0.001) |

### Step 4: 成果物

- ������ **`paper.md`** — 学術論文形式（Abstract 250語+、6章、参考文献12件、図8枚）
- 📄 **`report.md`** — 実験全結果レポート（日本語、図8枚埋め込み）
- 🖼️ **`figures/`** — 8枚の図（多様性分布、年齢相関、ROC曲線、ヒートマップ等）
- 📊 **CSVファイル3件** — binding_results.csv, icb_results.csv, diversity_results.csv