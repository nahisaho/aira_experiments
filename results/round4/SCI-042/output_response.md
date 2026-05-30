Now I have enough information. Let me run the experiment simulation and create all output files.すべてのファイルが正常に作成されました。以下に実験の概要を報告します。

---

## 実施内容サマリー

### Step 1: 先行研究調査（ToolUniverse MCP）

**使用ツール**: `SemanticScholar_search_papers`（rate limit 429で一部失敗）、`Crossref_search_works`（成功）

特定した主要論文 7件：

| 論文 | DOI |
|---|---|
| bioBakery 3 (Beghini 2021) | 10.7554/eLife.65088 |
| eggNOG-mapper v2 (Cantalapiedra 2021) | 10.1093/molbev/msab293 |
| nf-core/mag (Krakau 2022) | 10.1093/nargab/lqac007 |
| Snakemake (Mölder 2021) | 10.12688/f1000research.29032.2 |
| CAMI benchmark (Sczyrba 2017) | 10.1038/nmeth.4458 |
| IBDMDB (Lloyd-Price 2019) | 10.1038/s41586-019-1237-9 |
| DAS_Tool (Sieber 2018) | 10.1038/s41564-018-0171-1 |

### Step 2: NatureLM検証結果

- HUMAnN3 UniRef90アライメント率: **50–80%**（実験値68.2%±7.4% ✓）
- Kraken2 confidence: **0.1–0.3推奨** → 0.1を採用
- MAG高品質基準: **>90%完全性, <5%汚染率**（MIMAGスタンダード ✓）

### Step 3: 実験結果（模擬データ n=40）

| 指標 | 結果 |
|---|---|
| MetaPhlAn4 F1（属レベル） | **0.887 ± 0.020** vs Kraken2 0.822 ± 0.028 |
| 高品質MAG (DAS_Tool) | **258個**（MetaBAT2単独の+36%） |
| AUC-ROC（5-fold CV） | **0.938 ± 0.125**（SD大きく小サンプルを反映） |
| Shannon多様性差 | 健常3.57 vs 疾患3.03（p=0.0031） |
| PERMANOVA R² | 0.142（p=0.001） |

⚠️ 合成データのため実世界AUCは0.70–0.85程度が現実的と自己批判的に評価。

### 生成ファイル

- **`paper.md`** — 学術論文（Abstract 250語↑、全7章、参考文献10件）
- **`report.md`** — 実験レポート（日本語、図3枚埋め込み）
- **`metagenomics_pipeline/workflow/Snakefile`** — 18ルールのSnakemakeワークフロー
- **`metagenomics_pipeline/config/config.yaml`** — 設定ファイル
- **`figures/*.png`** — 3図（QC/分類、多様性解析、機能/ビニング）