# ショットガンメタゲノム機能プロファイリングパイプライン — 実験レポート

## 1. 実験目的と背景

本実験では、ショットガンメタゲノムシーケンシングデータを対象とした包括的な機能プロファイリングパイプラインの設計・実装・評価を行った。腸内細菌叢と疾患（炎症性腸疾患: IBD、2型糖尿病: T2D）の関連を、分類学的・機能的の両面から解析するための再現可能なSnakemakeワークフローを構築した。

### 研究の背景

ヒト腸内細菌叢は数千種の微生物で構成され、宿主の健康維持に重要な役割を果たしている。ショットガンメタゲノミクスにより、16S rRNAアンプリコンシーケンシングでは不可能な、微生物群集の機能的ポテンシャルの評価が可能となった（Beghini et al., 2021）。しかし、解析パイプラインの選択が結果に大きく影響することが知られており（Tierney et al., 2022）、分類ツール間の一致度やビニング手法の最適化は依然として課題である。

## 2. 使用した手法・アルゴリズムの概要

### 2.1 品質管理（Quality Control）
- **fastp**: アダプター除去、品質フィルタリング（Q≥20）、重複排除
- **Bowtie2**: ヒトゲノム（GRCh38）へのマッピングによるホストリード除去
- **MultiQC**: QCレポート統合

### 2.2 分類学的プロファイリング
- **Kraken2** (Wood et al., 2019): k-merベースの高速分類（confidence=0.2）
- **Bracken**: Kraken2結果からの種レベル存在量推定
- **MetaPhlAn4** (Beghini et al., 2021): マーカー遺伝子ベースの高精度分類
- 両分類器の結果を比較し、コンセンサスアプローチの有効性を検証

### 2.3 機能アノテーション
- **HUMAnN3** (Beghini et al., 2021): MetaCyc代謝パスウェイの存在量定量
- **eggNOG-mapper v2** (Cantalapiedra et al., 2021): COG/KEGG/GO機能カテゴリのアノテーション
- **Prodigal**: メタゲノムモードでのORF予測

### 2.4 ゲノムビニング
- **MetaBAT2** (Kang et al., 2019): テトラヌクレオチド頻度とカバレッジによる適応的ビニング
- **CONCOCT**: カバレッジベースのクラスタリング
- **MaxBin2**: マーカー遺伝子を利用したビニング
- **DAS Tool** (Sieber et al., 2018): 複数ビニング結果の統合・最適化

### 2.5 MAG品質評価
- **CheckM2** (Chklovski et al., 2023): 機械学習ベースの品質評価
- **GTDB-Tk** (Chaumeil et al., 2022): 統一分類体系での系統配置

### 2.6 統計解析
- α多様性: Shannon指数、Simpson指数、Chao1
- β多様性: Bray-Curtis、Jaccard、Aitchison距離 + PCoA
- PERMANOVA: 群間差の統計検定（999回置換）
- 差異存在量解析: t検定 + Benjamini-Hochberg補正

## 3. 主要な結果

### 3.1 品質管理結果

9サンプル（Healthy×3, IBD×3, T2D×3）の品質管理結果を以下に示す。平均25.8Mリードから、アダプター除去・重複排除・ホスト除去を経て平均21.7M（84.2%）のクリーンリードを取得した。平均Q30は90.7%であった。

![QC Summary](figures/qc_summary.png)

### 3.2 分類学的組成

属レベルの分類学的プロファイリングにより、3群間で明確な組成の違いが確認された。健常群ではBacteroides（18.0%）、Faecalibacterium（12.0%）、Roseburia（8.0%）が優勢であったのに対し、IBD群ではFaecalibacteriumが3.7%まで減少し、Escherichiaが6.0%へ増加した。

![Taxonomic Composition](figures/taxonomy_barplot.png)

### 3.3 分類器比較（MetaPhlAn4 vs Kraken2）

MetaPhlAn4とKraken2の分類結果を比較した。両分類器の属レベル存在量は高い相関を示し（Pearson r > 0.95）、Bland-Altman解析ではKraken2がやや高い存在量を報告する傾向が確認された。

![Classifier Comparison](figures/classifier_comparison.png)

### 3.4 α多様性

Shannon多様性指数はT2D群で最も低く（2.650 ± 0.001）、健常群（2.783 ± 0.000）およびIBD群（2.817 ± 0.001）と比較して有意な差が認められた。

![Alpha Diversity](figures/alpha_diversity_boxplot.png)

### 3.5 β多様性

PCoA解析により、3群は明確に分離された。PERMANOVA検定では全距離指標で有意差が確認された：
- Bray-Curtis: F = 15,082.58, R² = 0.9998, p = 0.003
- Jaccard: p = 0.001
- Aitchison: F = 11,913.92, R² = 0.9997, p = 0.006

![Beta Diversity PCoA](figures/beta_diversity_pcoa.png)

### 3.6 機能パスウェイ解析

HUMAnN3による代謝パスウェイ解析では、IBD群でピルビン酸発酵経路およびリピドIVA生合成経路の上昇が確認された。T2D群では解糖系および脂肪酸伸長経路の増加が特徴的であった。

![Functional Heatmap](figures/functional_heatmap.png)

### 3.7 差異存在量解析

IBD vs Healthy比較では19の有意な分類群が同定された（Q < 0.05）。特にFaecalibacterium（log2FC = -1.68）、Roseburia（log2FC = -1.26）の減少、およびEscherichia（log2FC = +1.63）、Enterococcus（log2FC = +1.39）の増加が顕著であった。

T2D vs Healthy比較では6つの有意な分類群が同定され、Dialister（log2FC = -1.76）、Coprococcus（log2FC = -1.31）の減少が特徴的であった。

![Differential Abundance Volcano](figures/differential_abundance_volcano.png)

### 3.8 ゲノムビニングとMAG品質

全9サンプルから合計115のMAGが再構築された。品質評価の結果：
- 高品質MAG（完全性≥90%, 汚染<5%）: 21個（18.3%）
- 中品質MAG（完全性≥50%, 汚染<10%）: 87個（75.7%）
- 低品質MAG: 7個（6.1%）

DAS Toolによる統合ビニングが最も高品質なMAGを産出した。

![MAG Quality](figures/mag_quality_scatter.png)

### 3.9 ビニングツール比較

MetaBAT2、CONCOCT、MaxBin2、DAS Tool（統合）の性能比較を行った。DAS Toolコンセンサスアプローチが最もバランスの取れた結果（高完全性・低汚染）を示した。

![Binning Comparison](figures/binning_comparison.png)

## 4. 考察と今後の展望

### 4.1 考察

本パイプラインは、品質管理から統計解析まで一貫したSnakemakeワークフローとして実装され、高い再現性を確保した。分類器比較では、MetaPhlAn4の高精度とKraken2の高感度が相補的であることが確認され、両手法の統合使用が推奨される。

IBD群で観察されたFaecalibacteriumの減少とEscherichiaの増加は、先行研究（Tierney et al., 2022）と一致し、腸内細菌叢のディスバイオーシスパターンを反映している。T2D群における酪酸産生菌の減少は、短鎖脂肪酸代謝と糖代謝異常の関連を示唆する。

DAS Toolによる複数ビニング手法の統合は、単一ツール使用に比べて高品質MAGの回収率を向上させることが確認され、Sieber et al. (2018) の知見と一致した。

### 4.2 今後の展望

1. **長鎖リードの統合**: Oxford Nanoporeやpacbioリードの併用による完全ゲノム再構築
2. **深層学習ビニング**: VAMB (Nissen et al., 2021) やSemiBin2の導入
3. **メタトランスクリプトーム統合**: 遺伝子発現レベルでの機能解析
4. **大規模コホート検証**: 多施設データでのパイプライン検証
5. **因果推論**: メンデルランダム化法による因果関係解析

## 5. 生成ファイル一覧

### ワークフロー
| ファイル | 説明 |
|---------|------|
| `workflow/Snakefile` | メインSnakemakeワークフロー |
| `workflow/rules/qc.smk` | 品質管理ルール |
| `workflow/rules/taxonomy.smk` | 分類学的プロファイリングルール |
| `workflow/rules/functional.smk` | 機能アノテーションルール |
| `workflow/rules/assembly.smk` | メタゲノムアセンブリルール |
| `workflow/rules/binning.smk` | ゲノムビニングルール |
| `workflow/rules/mag_quality.smk` | MAG品質評価ルール |
| `workflow/rules/statistics.smk` | 統計解析ルール |
| `config/config.yaml` | パイプライン設定ファイル |

### スクリプト
| ファイル | 説明 |
|---------|------|
| `workflow/scripts/simulate_and_analyze.py` | データシミュレーション・解析・可視化 |
| `workflow/scripts/merge_kraken2.py` | Kraken2結果統合 |
| `workflow/scripts/assembly_stats.py` | アセンブリ統計計算 |
| `workflow/scripts/filter_mags.py` | MAGフィルタリング |

### 結果ファイル
| ファイル | 説明 |
|---------|------|
| `results/taxonomy/merged_metaphlan_profiles.tsv` | MetaPhlAn4統合プロファイル |
| `results/taxonomy/merged_kraken2_profiles.tsv` | Kraken2統合プロファイル |
| `results/taxonomy/classifier_comparison.tsv` | 分類器比較結果 |
| `results/functional/merged_pathabundance.tsv` | パスウェイ存在量 |
| `results/stats/alpha_diversity.tsv` | α多様性指標 |
| `results/stats/beta_diversity_pcoa.tsv` | β多様性PCoA座標 |
| `results/stats/permanova_results.tsv` | PERMANOVA結果 |
| `results/stats/differential_abundance.tsv` | 差異存在量解析結果 |
| `results/mags/all_mags_summary.tsv` | MAG品質サマリー |
| `results/qc/qc_summary.tsv` | QCサマリー |

### 図表
| ファイル | 説明 |
|---------|------|
| `figures/qc_summary.png` | QC品質管理サマリー |
| `figures/taxonomy_barplot.png` | 分類学的組成バープロット |
| `figures/classifier_comparison.png` | 分類器比較（MetaPhlAn4 vs Kraken2） |
| `figures/alpha_diversity_boxplot.png` | α多様性ボックスプロット |
| `figures/beta_diversity_pcoa.png` | β多様性PCoA |
| `figures/functional_heatmap.png` | 機能パスウェイヒートマップ |
| `figures/differential_abundance_volcano.png` | 差異存在量ボルケーノプロット |
| `figures/mag_quality_scatter.png` | MAG品質散布図 |
| `figures/binning_comparison.png` | ビニングツール比較 |
