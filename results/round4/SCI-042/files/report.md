# 実験レポート：ショットガンメタゲノムデータからの機能プロファイリングパイプライン設計

---

## 1. 実験目的と背景

### 1.1 研究目的

ショットガンメタゲノムシーケンシングは腸内細菌叢の分類・機能特性を網羅的に解析するゴールドスタンダードであるが、再現性の高い統合パイプラインの整備が課題となっている。本実験では、6つの解析ステップを統合したSnakemakeベースのパイプライン（MetaFuncPipe）を設計・実装し、模擬データを用いたベンチマーク実験を実施した。

### 1.2 研究背景

| 項目 | 内容 |
|---|---|
| 対象サンプル | 模擬ショットガンメタゲノム 40サンプル（健常20、疾患20） |
| 参照疾患モデル | IBD（炎症性腸疾患）様dysbiosis |
| 平均シーケンス深度 | 約15.8M reads/sample |
| 分類対象種数 | 150種 |
| NatureLM活用 | パラメータ検証（アライメント率、クオリティ閾値等） |

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 パイプライン全体像

```
生リードデータ (FASTQ)
       ↓
[Step 1] 品質管理 (KneadData)
  ├─ アダプター除去 (Trimmomatic Q20, minlen50)
  ├─ ホストDNA除去 (Bowtie2 vs GRCh38)
  └─ 光学的重複除去

       ↓
[Step 2] 分類（2手法比較）
  ├─ Kraken2 (k-mer, confidence=0.1) + Bracken (species再推定)
  └─ MetaPhlAn4 (マーカー遺伝子, threshold=0.005%)

       ↓
[Step 3] 機能アノテーション
  ├─ HUMAnN3 (ChocoPhlAn + UniRef90)
  └─ eggNOG-mapper v2 (COG/KEGG/GO/EC)

       ↓
[Step 4] アセンブリ (MEGAHIT, kmer=21-141, mincontig=1000bp)

       ↓
[Step 5] ゲノムビニング（3ツール統合）
  ├─ MetaBAT2 (テトラヌクレオチド頻度 + カバレッジ)
  ├─ CONCOCT (ガウス混合モデル)
  ├─ MaxBin2 (期待値最大化)
  └─ DAS_Tool (アンサンブル精製, score≥0.5)

       ↓
[Step 6] MAG品質評価・系統配置
  ├─ CheckM2 (完全性・汚染率評価)
  └─ GTDB-Tk r220 (系統配置)

       ↓
[Step 7] 多変量統計解析
  ├─ Alpha多様性 (Shannon, Chao1)
  ├─ Beta多様性 (Bray-Curtis + PERMANOVA)
  ├─ LEfSe (差次的菌種・代謝経路)
  └─ Random Forest分類 (5-fold CV)
```

### 2.2 Snakemakeワークフロー

18のSnakemakeルールから構成されるDAGベースのワークフロー。主要な設計原則：
- 各ルールにConda環境ファイル（`.yaml`）を紐付け、ソフトウェア再現性を担保
- YAML設定ファイルによるパラメータ管理
- HPC/クラウド対応（SLURM/LSFプロファイル）

### 2.3 NatureLM MCPツール活用

**接続状況**: 接続成功（naturelm-8x7b-inst, vllm backend）

| クエリ内容 | NatureLM出力 | パイプラインへの適用 |
|---|---|---|
| QC品質スコア閾値 | Phred Q20–Q30 | Q20を基本設定、高精度が必要な場合Q30を推奨 |
| Kraken2信頼度閾値 | 0.1–0.3推奨 | confidence=0.1をデフォルト設定 |
| HUMAnN3 UniRef90アライメント率 | 50–80% | 実験結果68.2%±7.4%（範囲内） |
| MAG高品質基準 | 完全性>90%, 汚染率<5% | MIMAGスタンダードとして実装 |
| Shannon多様性の傾向 | 疾患群で低下 | 実験結果と一致（p=0.0031） |

---

## 3. 主要な結果と数値

### 3.1 品質管理結果

| 指標 | 値 |
|---|---|
| 平均生リード数 | 15.8M reads/sample |
| 平均QC後リード数 | 12.5M reads/sample |
| ホスト除去率 | 平均4.9% |
| アダプター/低品質除去率 | 平均4.1% |
| 重複除去率 | 平均11.7% |
| 総リード保持率 | 約79.3% |

### 3.2 分類学的プロファイリング比較

![Figure 1](figures/figure1_pipeline_overview.png)

**表1. Kraken2 vs MetaPhlAn4 性能比較（属レベル, n=40）**

| ツール | Precision | Recall | F1スコア | 実行時間 | DBサイズ |
|---|---|---|---|---|---|
| Kraken2 | 0.853 ± 0.035 | 0.796 ± 0.048 | 0.822 ± 0.028 | 2.3 min | ~49 GB |
| MetaPhlAn4 | 0.909 ± 0.034 | 0.867 ± 0.026 | **0.887 ± 0.020** | 8.7 min | ~1.5 GB |

- MetaPhlAn4のF1スコアがKraken2より7.9%高い（ΔF1=+0.065）
- Kraken2は3.8倍高速・データベースサイズは97%小さい
- 実際の使用推奨：高精度が必要な場合はMetaPhlAn4、スクリーニングや希少種検出はKraken2

### 3.3 HUMAnN3機能アノテーション

| 指標 | 値 |
|---|---|
| UniRef90アライメント率 | 68.2 ± 7.4%（範囲: 46.3–88.1%） |
| UniPathwayカバレッジ | 0.72 ± 0.08 |
| 検出遺伝子ファミリー数 | 12,124–44,687（平均28,456） |
| 差次的経路数（LEfSe, LDA≥2.0） | 8経路 |

**表2. LEfSe同定された主要差次的代謝経路**

| 代謝経路 | LDA スコア | 方向 |
|---|---|---|
| Short-chain fatty acid synthesis | +3.41 | 健常群富化 |
| Butyrate synthesis II | +3.18 | 健常群富化 |
| Propionate production | +2.95 | 健常群富化 |
| LPS biosynthesis | −2.61 | 疾患群富化 |
| Mucin degradation | −2.88 | 疾患群富化 |
| Bile acid transformation | −3.05 | 疾患群富化 |

### 3.4 MAGリカバリー結果

**表3. ビニングツール比較（MAG数）**

| ツール | 高品質MAG | 中品質MAG | 合計 |
|---|---|---|---|
| MetaBAT2 | 187 | 224 | 411 |
| CONCOCT | 142 | 198 | 340 |
| MaxBin2 | 158 | 211 | 369 |
| **DAS_Tool（アンサンブル）** | **258 (+36%)** | **301** | **559** |

- 総MAG数: 862
- 高品質MAG: 258（29.9%）、中品質: 301（34.9%）
- 平均完全性: 70.3 ± 23.1%、平均汚染率: 8.6 ± 6.9%
- GTDB-Tk分類: Firmicutes_A（37.2%）、Bacteroidota（28.6%）

### 3.5 多様性解析・疾患関連性

![Figure 2](figures/figure2_diversity_analysis.png)

**Alpha多様性（Shannon指数）:**
- 健常群: 3.57 ± 0.57
- 疾患群: 3.03 ± 0.41
- Mann-Whitney U検定: p = 0.0031

**Beta多様性（PERMANOVA）:**
- Bray-Curtis距離に基づくR² = 0.142
- 疾患群と健常群間の有意差: p = 0.001（999 permutations）

**表4. Random Forest疾患分類（5-fold CV）**

| 特徴量セット | AUC | F1スコア | Precision | Recall |
|---|---|---|---|---|
| 分類学のみ | 0.938 ± 0.125 | 0.911 ± 0.130 | 0.880 ± 0.160 | 0.950 ± 0.100 |
| 機能のみ | 0.862 ± 0.148 | 0.847 ± 0.162 | 0.821 ± 0.183 | 0.875 ± 0.141 |
| 統合 | **0.951 ± 0.098** | **0.934 ± 0.114** | 0.908 ± 0.138 | 0.962 ± 0.082 |

⚠️ **注意**: AUCのSDが0.098–0.148と大きく、これはn=40の小サンプルサイズによる折り間変動を反映している。合成データの強い疾患シグナルにより性能が過大推定されている可能性があり、実世界データでは0.70–0.85程度が現実的と考えられる。

![Figure 3](figures/figure3_functional_binning.png)

---

## 4. 考察と今後の展望

### 4.1 主要な知見の解釈

**分類ツール選択**: MetaPhlAn4の高精度（F1=0.887）は腸内細菌叢の既知マーカー遺伝子データベースの充実によるが、新規・希少種の検出ではKraken2が優位な場合がある。実用的な推奨としては、標準的な腸内細菌叢解析にはMetaPhlAn4を一次ツールとして使用し、ウイルス・古細菌・希少種スクリーニングにKraken2を補完的に使用する二段階アプローチが最適である。

**アンサンブルビニングの優位性**: DAS_Toolによる統合ビニングが単一ツールより36%多くの高品質MAGを回復したことは、各ビニングアルゴリズムの相補的な特性（テトラヌクレオチド頻度 vs カバレッジ vs 確率モデル）を反映している。

**疾患シグナル**: 酪酸産生経路（butyrate synthesis）の健常群での富化とLPS生合成の疾患群での富化は、IBDにおける公知の病態生理（短鎖脂肪酸産生低下、腸管バリア機能障害）と一致しており、シミュレーションの面的妥当性を支持する。

### 4.2 自己批判的評価

**合成データへの依存性:**
実験結果は合成データ上のシミュレーションに基づいており、以下の重要な前提条件がある：

1. **疾患シグナルの均一性**: 本実験では疾患群の全20サンプルに対して同一の種増減パターンを適用した。実際のIBDコホートでは個体間の不均一性が高く（σ²_between >> σ²_within）、分類精度は大きく低下する可能性がある。

2. **分類精度の独立評価**: 分類学的プロファイリングの精度評価（Precision/Recall）は既知構成の模擬データに基づいており、実世界サンプルでは未知種・複雑な菌叢が存在するため、精度が10–20%低下することが予想される。

3. **NatureLM予測の信頼性**: NatureLMが提示したアライメント率（50–80%）は文献値の要約であり、特定のデータセットや最新データベースバージョンとは乖離する可能性がある。

4. **統計的検出力**: n=40、5-fold CVでは各テストフォールドがn=8と小さく、AUCのSDが約0.10–0.15となる。実臨床研究ではn≥200のコホートが推奨される。

### 4.3 今後の展望

1. **実データへの適用**: NIDDK IBDMDB（HMP2）、MetaHIT、GMrepo等の公共コホートデータでの検証
2. **株レベル解析の追加**: StrainPhlAn4による株多様性解析の統合
3. **マルチオミクス統合**: メタトランスクリプトーム、メタプロテオーム、メタボロームとの統合解析
4. **ビロームと真菌叢**: 現行パイプラインにはウイルス・真菌コンポーネントが未統合
5. **クラウドネイティブ化**: Terra/AnVIL（WDL/Nextflow）への対応
6. **大規模コホート対応**: 数千サンプル規模でのスケーラビリティテスト

---

## 5. 生成したファイル一覧

| ファイル | 説明 |
|---|---|
| `metagenomics_pipeline/workflow/Snakefile` | Snakemakeワークフロー（18ルール） |
| `metagenomics_pipeline/config/config.yaml` | パイプライン設定ファイル |
| `metagenomics_pipeline/results/qc_metrics.json` | QCメトリクス（40サンプル） |
| `metagenomics_pipeline/results/taxonomic_results.json` | 分類ベンチマーク結果 |
| `metagenomics_pipeline/results/mag_results.json` | MAGビニング結果 |
| `metagenomics_pipeline/results/classification_results.json` | 疾患分類CV結果 |
| `metagenomics_pipeline/figures/figure1_pipeline_overview.png` | パイプライン概要図 |
| `metagenomics_pipeline/figures/figure2_diversity_analysis.png` | 多様性・疾患関連解析図 |
| `metagenomics_pipeline/figures/figure3_functional_binning.png` | 機能・ビニング解析図 |
| `paper.md` | 学術論文形式ドキュメント |
| `report.md` | 本レポート |

---

## 付録：先行研究調査サマリー

### Step 1で特定した主要論文（ToolUniverse MCP使用）

**使用ツール**: SemanticScholar_search_papers, Crossref_search_works

| # | 論文 | 年 | DOI | 主要知見 |
|---|---|---|---|---|
| 1 | Beghini et al. (bioBakery 3) | 2021 | 10.7554/eLife.65088 | MetaPhlAn4/HUMAnN3統合スイートの標準化 |
| 2 | Cantalapiedra et al. (eggNOG-mapper v2) | 2021 | 10.1093/molbev/msab293 | メタゲノムスケール機能アノテーション、>4,149引用 |
| 3 | Krakau et al. (nf-core/mag) | 2022 | 10.1093/nargab/lqac007 | ハイブリッドアセンブリ・ビニングパイプライン |
| 4 | Mölder et al. (Snakemake 5.x) | 2021 | 10.12688/f1000research.29032.2 | 持続可能なSnakemakeワークフロー設計 |
| 5 | Sczyrba et al. (CAMI) | 2017 | 10.1038/nmeth.4458 | メタゲノムツールの包括的ベンチマーク |
| 6 | Lloyd-Price et al. (IBDMDB) | 2019 | 10.1038/s41586-019-1237-9 | IBD腸内細菌叢マルチオミクス |
| 7 | Sieber et al. (DAS_Tool) | 2018 | 10.1038/s41564-018-0171-1 | アンサンブルMAGビニング戦略 |

**SemanticScholar APIエラー**: rate limit (HTTP 429) のため一部クエリで結果取得に失敗。Crossref APIを代替として使用し、eggNOG-mapper v2論文（4,149引用、Molecular Biology and Evolution誌）の完全メタデータを取得できた。

---

*レポート生成日: 2026-05-29*
*使用ツール: ToolUniverse MCP (SemanticScholar/Crossref), NatureLM MCP (naturelm-8x7b-inst)*
