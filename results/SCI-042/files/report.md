# ショットガンメタゲノム機能プロファイリングパイプライン — 解析レポート

**DRAFT — NOT FOR DISTRIBUTION**
**作成日**: 2026-05-23
**バージョン**: v1.0.0
**ワークフローエンジン**: Snakemake ≥ 8.0

---

## 1. 実験目的と背景

本パイプラインは、ショットガンメタゲノムシーケンスデータから腸内細菌叢の**分類学的・機能的プロファイリング**を包括的に行い、疾患関連微生物群の同定と機能変化の解明を目的とする。

### 背景
- 腸内細菌叢の組成変化（dysbiosis）は、IBD、2型糖尿病、大腸がんなど多くの疾患と関連が示されている
- ショットガンメタゲノミクスは16S rRNAアンプリコン解析に比べ、種レベルの解像度と機能ポテンシャルの同時解析を可能にする
- 複数の分類・ビニングツールの統合により、単一ツールの偏りを補正し堅牢な結果を得る

### 解析フレームワーク
6段階のモジュール構成で、各段階は独立したSnakemakeルールとして実装：

| Phase | 目的 | 主要ツール |
|-------|------|-----------|
| Step 1 | 品質管理 | fastp, Clumpify, Bowtie2 |
| Step 2 | 分類学的プロファイリング | Kraken2/Bracken, MetaPhlAn 4 |
| Step 3 | 機能アノテーション | HUMAnN 3, eggNOG-mapper |
| Step 4 | ゲノムビニング | MetaBAT2, CONCOCT, MaxBin2, DAS Tool |
| Step 5 | MAG品質評価 | CheckM2, GTDB-Tk |
| Step 6 | 統計解析 | scikit-bio, ALDEx2, MaAsLin 2 |

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 品質管理（Step 1）

```
Raw reads → fastp (adapter/quality trim) → Clumpify (dedup) → Bowtie2 (host removal) → Clean reads
```

- **fastp v0.23.4**: アダプター自動検出＋品質トリミング（Q≥20）、低複雑度フィルタ（entropy < 0.3）、最小リード長50bp
- **Clumpify (BBTools v39.06)**: 光学的・PCR重複除去。optical distance = 2500（HiSeq 4000/NovaSeq対応）
- **Bowtie2 v2.5.3**: GRCh38参照ゲノムへのマッピングによるヒトリード除去（MAPQ ≥ 20で非マップリードを抽出）
- **MultiQC v1.21**: 全QC指標の統合レポート生成

### 2.2 分類学的プロファイリング（Step 2）

2つの相補的アプローチを並行実行し、一致性を評価：

| ツール | 手法 | 長所 | 短所 |
|--------|------|------|------|
| **Kraken2 v2.1.3** + Bracken | k-mer完全一致 (confidence=0.2) | 高速、広範なDB | 偽陽性リスク |
| **MetaPhlAn 4 v4.1.1** | クレード特異的マーカー遺伝子 | 高精度、低偽陽性 | 新規種の検出限界 |

- Brackenによるベイズ再推定で種レベルの存在量を補正
- 2ツール間のSpearman相関・Jaccard類似度で一致性を定量評価
- **推奨**: MetaPhlAn 4を主軸とし、Kraken2/Brackenは感度補完として使用

### 2.3 機能アノテーション（Step 3）

```
Clean reads → HUMAnN 3 (MetaCyc pathways) → CPM normalization → KO regrouping
Contigs → Prodigal (ORF prediction) → eggNOG-mapper (COG/KEGG/GO)
```

- **HUMAnN 3 v3.9**: MetaPhlAn 4のプロファイルを事前情報として使用し、UniRef90に対するtiered searchでgene family/pathway abundanceを定量。CPM正規化後にKEGG Ortholog (KO) への再グルーピング
- **eggNOG-mapper v2.1.12**: アセンブル済みコンティグからのORF予測→DIAMOND BLASTPによるeggNOG5データベース検索→COG, KEGG, GO, CAZy等の多層アノテーション

### 2.4 ゲノムビニング（Step 4）

3つのビニングツールの結果をDAS Toolで統合する**コンセンサスビニング戦略**：

```
Contigs → Coverage calculation
       ├→ MetaBAT2 (TNF + coverage)
       ├→ CONCOCT  (coverage + composition, chunked)
       └→ MaxBin2  (EM-based)
       → DAS Tool (consensus, score ≥ 0.5)
```

- **MetaBAT2 v2.17**: テトラヌクレオチド頻度（TNF）＋カバレッジベースの適応的距離計算
- **CONCOCT v1.1.0**: ガウス混合モデルによる変分ベイズクラスタリング（10kb chunk分割）
- **MaxBin2 v2.2.7**: EMアルゴリズムベースの確率的ビニング
- **DAS Tool v1.1.7**: 3ツールの結果からsingle-copy geneベースのスコアリングで最適ビンセットを選択

### 2.5 MAG品質評価（Step 5）

- **CheckM2 v1.0.2**: 機械学習ベースのcompleteness/contamination推定（MIMAG基準適用）
  - High Quality: completeness ≥ 90%, contamination ≤ 5%
  - Medium Quality: completeness ≥ 50%, contamination ≤ 10%
- **GTDB-Tk v2.4.0** (GTDB r220): ANI比較＋pplacer系統配置による分類学的アサインメント

### 2.6 多変量統計解析（Step 6）

| 解析 | 手法 | 検定 |
|------|------|------|
| α多様性 | Shannon, Simpson, observed, Chao1 | Mann-Whitney U検定 |
| β多様性 | Bray-Curtis, Jaccard, Aitchison PCoA | PERMANOVA (999 permutations) |
| Differential abundance | CLR-transformed Mann-Whitney | BH法FDR補正 (q < 0.05) |
| 効果量 | log₂ fold change | |log₂FC| ≥ 1.0 |

- **多重検定補正**: Benjamini-Hochberg法でFDR制御
- **組成データ処理**: CLR変換（Centered Log-Ratio）でAitchison距離を算出
- **有病率・存在量フィルタ**: 全サンプルの10%以上で検出、平均相対存在量≥0.1%の特徴量のみ解析対象

---

## 3. 主要な結果と数値

> ⚠️ 本セクションはパイプライン設計段階のため、実データ実行後に数値が入ります。

### 3.1 品質管理の期待指標

| 指標 | 期待値 |
|------|--------|
| リード生存率（QCフィルタ後） | 85–95% |
| ホストリード除去率 | 1–30%（サンプル依存） |
| 重複排除率 | 5–15% |
| 最終クリーンリード数/サンプル | ≥10M paired-end reads |

### 3.2 分類プロファイリングの期待指標

| 指標 | Kraken2/Bracken | MetaPhlAn 4 |
|------|-----------------|-------------|
| 検出種数 | 200–500 | 100–300 |
| 分類率 | 60–80% | 30–60% |
| 偽陽性率 | 中 | 低 |
| ツール間Spearman ρ | 0.6–0.8（共通種） |

### 3.3 機能アノテーションの期待指標

| 指標 | 期待値 |
|------|--------|
| HUMAnN 3 alignment rate | 50–80% |
| 検出MetaCyc pathway数 | 200–400 |
| eggNOG annotation rate | 60–80% |

### 3.4 ゲノムビニングとMAG品質の期待指標

| 指標 | 期待値 |
|------|--------|
| コンティグ N50 | 5–50 kbp |
| DAS Tool統合後ビン数/サンプル | 10–50 |
| High-quality MAG割合 | 10–30% |
| Medium-quality MAG割合 | 30–50% |

### 3.5 統計解析の期待出力

| 解析 | 出力 |
|------|------|
| α多様性 | 群間比較の箱ひげ図＋p値 |
| β多様性 | PCoA散布図＋PERMANOVA F統計量・p値 |
| Differential abundance | 有意特徴量リスト（q値・効果量付き） |

---

## 4. 考察と今後の展望

### 4.1 設計上の考慮事項

1. **コンセンサスビニング戦略**: MetaBAT2, CONCOCT, MaxBin2の3ツール統合により、単一ツール使用に比べてMAGの質と数が向上する（Parks et al., 2017; Sieber et al., 2018）
2. **分類ツールの相補性**: Kraken2の高感度とMetaPhlAn 4の高特異度を組み合わせることで、検出漏れと偽陽性の両方を軽減
3. **組成データの統計的取り扱い**: CLR変換とAitchison距離の採用により、相対存在量データの組成的性質（compositional data）に起因するspurious correlationを回避
4. **多重検定補正**: BH法によるFDR制御を全差異存在量解析に適用し、偽発見率を5%以下に制御

### 4.2 制限事項

- **参照データベース依存性**: 全分類・機能アノテーションは既知ゲノム/遺伝子に依存し、未知種・新規機能の同定には限界がある
- **Co-assembly vs per-sample assembly**: 本パイプラインはサンプル別アセンブリを採用。低存在量種のMAG構築にはco-assemblyが有利な場合がある
- **因果推論の限界**: 横断的デザインでは因果関係の推定は不可能。縦断的データや介入研究が必要
- **サンプルサイズ**: 効果量と検出力はサンプル数に依存。事前に検出力分析を推奨

### 4.3 今後の展望

1. **Strain-level解析**: StrainPhlAn 4やinStrainを統合し、株レベルの多様性と伝播パターンを解析
2. **メタトランスクリプトーム統合**: RNA-seqデータとの統合により、機能ポテンシャルから実際の遺伝子発現へ
3. **機械学習分類器**: Random Forest / XGBoostベースの疾患予測モデル構築（LOOCV評価）
4. **Co-occurrence network**: SparCC/SPIEC-EASIによる微生物間相互作用ネットワーク推定
5. **縦断的解析**: 時系列データ対応の拡張（CLAMMによる縦断的differential abundance）

---

## 5. パイプライン実行方法

### 前提条件
```bash
# Conda/Mamba + Snakemake
conda install -c bioconda -c conda-forge snakemake mamba

# データベースの事前ダウンロード（config/config.yaml のパスを設定）
kraken2-build --standard --db /db/kraken2/k2_standard_20240112
metaphlan --install --bowtie2db /db/metaphlan4/
humann_databases --download chocophlan full /db/humann3/
humann_databases --download uniref uniref90_diamond /db/humann3/
download_eggnog_data.py -y --data_dir /db/eggnog/eggnog5
```

### 実行コマンド
```bash
# ドライラン（依存関係確認）
snakemake -s workflow/Snakefile --configfile config/config.yaml -n

# 本番実行（8コア、conda環境自動構築）
snakemake -s workflow/Snakefile --configfile config/config.yaml \
    --use-conda --cores 8 --rerun-incomplete

# クラスタ実行（SLURM）
snakemake -s workflow/Snakefile --configfile config/config.yaml \
    --use-conda --profile slurm --jobs 50
```

### DAG可視化
```bash
snakemake -s workflow/Snakefile --dag | dot -Tsvg > figures/pipeline_dag.svg
```

---

## 6. 生成ファイル一覧

```
workspace/
├── report.md                                          # 本レポート
├── config/
│   ├── config.yaml                                    # パイプライン設定ファイル
│   └── samples.tsv                                    # サンプルシート
├── workflow/
│   ├── Snakefile                                      # メインワークフロー
│   ├── rules/
│   │   ├── qc.smk                                     # Step 1: 品質管理ルール
│   │   ├── taxonomy.smk                               # Step 2: 分類プロファイリングルール
│   │   ├── functional.smk                             # Step 3: 機能アノテーションルール
│   │   ├── assembly_binning.smk                       # Step 4: アセンブリ＋ビニングルール
│   │   ├── mag_quality.smk                            # Step 5: MAG品質評価ルール
│   │   └── statistics.smk                             # Step 6: 統計解析ルール
│   ├── envs/
│   │   ├── qc.yaml                                    # QC用conda環境
│   │   ├── taxonomy.yaml                              # 分類用conda環境
│   │   ├── functional.yaml                            # 機能アノテーション用conda環境
│   │   ├── assembly.yaml                              # アセンブリ用conda環境
│   │   ├── binning.yaml                               # ビニング用conda環境
│   │   ├── mag_quality.yaml                           # MAG品質評価用conda環境
│   │   └── statistics.yaml                            # 統計解析用conda環境
│   └── scripts/
│       ├── compare_classifiers.py                     # Kraken2 vs MetaPhlAn4 比較
│       ├── filter_mags.py                             # MAG品質フィルタリング
│       ├── alpha_diversity.py                         # α多様性計算＋可視化
│       ├── beta_diversity.py                          # β多様性PCoA＋PERMANOVA
│       ├── taxonomic_barplot.py                       # 分類組成バープロット
│       ├── differential_abundance.py                  # 差異存在量解析
│       ├── functional_heatmap.py                      # 機能ヒートマップ
│       └── mag_quality_plot.py                        # MAG品質散布図
├── results/                                           # 解析結果出力先
│   ├── qc/                                            # QCレポート・フィルタ済みリード
│   ├── taxonomy/                                      # 分類プロファイル
│   ├── functional/                                    # 機能アノテーション結果
│   ├── assembly/                                      # アセンブリコンティグ
│   ├── binning/                                       # ビニング結果
│   ├── mag_quality/                                   # MAG品質・系統配置
│   └── statistics/                                    # 統計解析結果
├── figures/                                           # 出版品質の図表
│   ├── alpha_diversity_boxplot.svg                    # α多様性群間比較
│   ├── beta_diversity_pcoa.svg                        # β多様性PCoA
│   ├── taxonomic_barplot.svg                          # 分類組成バープロット
│   ├── functional_heatmap.svg                         # 差異代謝パスウェイヒートマップ
│   ├── mag_quality_scatter.svg                        # MAG品質散布図
│   └── pipeline_dag.svg                               # パイプラインDAG
└── logs/
    └── process-log.jsonl                              # 実行トレースログ
```

---

## 7. 参考文献

1. Wood DE, Lu J, Langmead B. (2019) Improved metagenomic analysis with Kraken 2. *Genome Biology* 20:257
2. Blanco-Míguez A et al. (2023) Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. *Nature Biotechnology* 41:1633–1644
3. Beghini F et al. (2021) Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. *eLife* 10:e65088
4. Kang DD et al. (2019) MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction. *PeerJ* 7:e7359
5. Alneberg J et al. (2014) Binning metagenomic contigs by coverage and composition. *Nature Methods* 11:1144–1146
6. Sieber CMK et al. (2018) Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. *Nature Microbiology* 3:836–843
7. Chklovski A et al. (2023) CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality. *Nature Methods* 20:1203–1212
8. Chaumeil PA et al. (2022) GTDB-Tk v2: memory friendly classification with the Genome Taxonomy Database. *Bioinformatics* 38:5315–5316
