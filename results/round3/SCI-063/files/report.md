# 最小ゲノムの合理的設計と合成フレームワーク：統合バイオインフォマティクスパイプライン

**DRAFT — NOT FOR DISTRIBUTION**

---

## Abstract

#            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo .git .github .gitignore AGENTS.md data figures logs results src 0.0 → 15.1 per gene）、追加フィルタリングパスの必要性が示された。本パイプラインはオープンな計算フレームワークとして最小細胞研究への応用が期待される。"___Begin___Command_Done_Marker___$1）トランスポゾン変異導入（TN-seq）データに基づく機械学習モデルによる必須遺伝子予測（最良AUROC: 0.9875 ± 0.0055, F1: 0.9217 ± 0.0343）、（2）RSCU重みとGC含量制御を組み合わせたコドン最適化（CAI向上: +0.197; 0.420 → 0.617）、（3）複製方向バイアス解析（リーディング鎖遺伝子比率: 66.8%）、（4）オペロン構造に基づくゲノムリファクタリング戦略（推定圧縮量: 9,629 bp / 1.81%）、（5）3段階階層的 Gibson Assembly 設計（Tier-1: 107フラグメント, Tier-2: 11サブ染色体, Tier-3: 酵母TAR）。重要な知見として、コドン最適化後に反復配列数が増加するトレー} 

---

## 1. 実験目的と背景

### 1.1 研究背景

#.git .github .gitignore AGENTS.md data figures logs results src tests 
#            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             2016年、Hutchison ら（JCVI）は全ゲノム化学合成と設計・合成・試験の反復サイクルにより、473遺伝子・531 kbpの JCVI-syn3.0 を構築した（Hutchison et al., 2016）。この研究は合成生物学の革命的マイルストーンであったが、初期設計がトランスポゾンデータの不足により失敗したことは、データ駆動型設計手法の重要性を浮き彫りにした。JCVI-syn3.0 は 473 遺伝子のうち 149 .git .github .gitignore AGENTS.md data figures logs results Tests  src/__pycache__ src/codon_optimizer.py src/essential_gene_predictor.py src/genome_architect.py src /pipeline.py 


#
-転写干渉を最小化するための遺伝子配置と遺伝子方向のバイアス設計。第四に、機能的に重複した遺伝子を統合するリファクタリング。第五に、100 kbp 以上'REPORT_EOF' DNA を確実に組み立てる階層的アセンブリ設計。

### 1.2 先行研究調査（MCP ツール使用状況の記録）

ToolUniverse MCP経由でSemanticScholar_search_papersを試行したが、API error 400（パラメータ形式の不整合）が返された。代替としてPubMed_search_articles および Crossref_search_works ツールを成功裏に使用し、以下の文献を特定した。

| 文献 | 手法 | 主要知見 |
|------|------|---------|
| Hutchison et al. 2016 (Science) | 全ゲノム化学合成 + 反復設計試験 | JCVI-syn3.0: 473遺伝子, 531 kbp; 初期設計失敗はquasi-essential遺伝子の見落とし |
| Pelletier et al. 2022 (Trends Cell Biol) | JCVI-syn3A 細胞力学解析 | 最小細胞での細胞分裂は多遺伝子的基盤をもつ; 表面積/体積比と膜曲率が鍵 |
| Billmyre et al. 2025 (PLoS Biol) | TN-seq + Random Forest | 1,465 必須遺伝子を予測; ヒトオルソログなし302遺伝子は抗真菌薬ターゲット候補 |
| Levitan et al. 2020 (Curr Genet) | TN-seq + ML (3酵母種比較) | 挿入密度分布と挿入位置偏りがML予測精度に影響; Cross-species ortholog活用で精度向上 |
| Menuhin-Gruman et al. 2025 (Sci Adv) | AI予測遺伝子融合 (STABLES) | GOI-必須遺伝子融合で進化安定性向上; S. cerevisiaeで実験検証済み |
| Demissie et al. 2025 (JMB) | コドン最適化ツール比較 | 単一指標（CAI）では不十分; GC含量・mRNA二次構造・CPBの統合評価を推奨 |
| Simons 2021 (Studies Hist Phil Sci) | 技術科学的分析 | 最小ゲノム概念の認識論的意義; 機能未知遺伝子の多さは設計概念を複雑化 |
| Geng et al. 2026 (Sci Reports) | Tripleknock (深層学習) | E. coli三重遺伝子ノックアウト致死性予測F1=0.77; FBAの20倍高速 |

### 1.3 先行研究の課題・限界

1. **必須遺伝子予測の精度**: TN-seq データのみに依存した予測は、quasi-essential遺伝子（条件的に必須）を見落とす。複数の成長条件下でのデータ統合が必要。
2. **コドン最適化の多目標問題**: CAIの最大化だけでは不十分。反復配列の生成、GC含量バランス、mRNA二次構造との相互作用を同時に考慮する必要がある（Demissie et al. 2025）。
3. **機能未知遺伝子の扱い**: JCVI-syn3.0の149遺伝子（約31%）は機能未知であり、設計空間の探索を困難にしている。
4. **: 既存研究は各コンポーネントを独立に扱っており、設計フェーズ間の整合性（例: コドン最適化が反復配列に与える影響）を体系的に評価していない。**統合パイプラ'REPORT_EOF'

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 システム構成

4つの Python モジュールで構成される：

- **essential_gene_predictor.py**: TN-seq 特徴量シミュレーション + 3モデル交差検証
- **codon_optimizer.py**: RSCU-CAI 最適化 + 反復配列検出 + 安定性スコアリング
- **genome_architect.py**: オペロン設計 + 複製バイアス + リファクタリング + アセンブリ計画
- **pipeline.py**: 統合オーケストレーション + 可視化 + ログ出力

### 2.2 必須遺伝子予測の数理モデル

AUROC は Wilcoxon-Mann-Whitney 統計として：

$$\text{AUROC} = \int_0^1 \text{TPR}(t) \, d[\text{FPR}(t)] = P(\hat{y}_{pos} > \hat{y}_{neg})$$

CAI (Sharp & Li, 1987):

$$\text{CAI} = \exp\left(\frac{1}{L}\sum_{i=1}^{L} \ln w_i\right), \quad w_i = \frac{\text{RSCU}(c_i)}{\max_{c' \sim c_i} \text{RSCU}(c')}$$

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$GC補正付き）：

$$\text{score}(c) = \text{RSCU}(c) + \lambda \cdot (\text{GC}^* - \overline{\text{GC}}_{\text{current}}) \cdot \text{GC}(c)$$

 $\lambda = 0.3$, $\text{GC}^* = 0.33$（Mycoplasma mycoides 目標値）。

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }

$$S(g) = \text{CAI}(g) \times \max\!\left(0, 1 - \frac{|\text{GC}(g) - \text{GC}^*|}{0.2}\right) \times e^{-N_r(g)/5}$$

### 2.3 アセンブリ設計

3段階の階層的 Gibson Assembly を設計した：
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             } 
- **Tier 2**: Tier-1 × ~10 → Gibson Assembly による ~50 kbp サブ染色体
- **Tier 3**: Tier-2 サブ染色体 → 酵母 TAR (Transformation-Associated Recombination) による全ゲノム組み立て

---

## 3. 主要な結果と数値

### 3.1 必須遺伝子予測

![Figure 1: Model Comparison](figures/fig1_model_comparison.png)

**表1: 5分割交差検証メトリクス（平均 ± 標準偏差）**

| モデル | AUROC | F1 | Precision | Recall |
|--------|-------|-----|-----------|--------|
| Random Forest | 0.983 ± 0.005 | 0.907 ± 0.033 | 0.923 ± 0.023 | 0.892 ± 0.048 |
| Gradient Boosting | 0.981 ± 0.005 | 0.894 ± 0.028 | 0.906 ± 0.031 | 0.883 ± 0.041 |
| Logistic Regression | **0.988 ± 0.006** | **0.922 ± 0.034** | 0.907 ± 0.048 | **0.938 ± 0.028** |

3モデルすべてが AUROC > 0.98 を達成した。ロジスティック回帰が最高の AUROC（0.988 ± 0.006）および F1（0.922 ± 0.034）を示した。RF・GBが高い適合率（Precision > 0.9）を示した一方、LRが高い再現率（Recall 0.938）を示すことは、必須遺伝子の見落とし（偽陰性）を最小化する観点で重要である。

![Figure 2: Feature Importances](figures/fig2_feature_importances.png)

TN-seq 挿入密度（tn_insertion_density）が最重要特徴量であることが確認された。種間保存スコア、発現レベル、タンパク質間相互作用度数がこれに続いた。

### 3.2 コドン最適化

![Figure 3: Codon Optimization](figures/fig3_codon_optimization.png)

60遺伝子の平均 CAI は 0.420 から 0.617 へ有意に向上した（**+0.197, +46.9%**）。GC 含量は平均 0.465 から 0.342 へと目標値（0.33）に近づいた。

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$Demissie et al. 2025）が単一指標CAIの限界として指摘したトレードオフを定量的に実証するものである。追加の反復配列除去パス（例: 反復位置での同義コドンへの置換）が必要である。

### 3.3 ゲノム構造解析

![Figure 4: Genome Composition](figures/fig4_genome_composition.png)

syn3.0 模倣ゲノム（473遺伝子）の機能カテゴリ分布：翻訳関連 35.3%（167遺伝子） 30.7%（145遺伝子）、代謝 12.7%、ゲノム処理 10.1%、膜タンパク 9.1%、細胞分裂 2.1%。予測必須遺伝子数は 335（71% の遺伝子）。

**複製方向バイアス**: 遺伝子の 66.8% がリーディング鎖に配置（閾値 55% を大きく超える）。必須遺伝子の 65.7% がリーディング鎖に配置。これは複製-転写干渉を最小化し、ゲノム安定性を高める重要な特性である。

![Figure 6: Strand Bias](figures/fig6_strand_bias.png)

### 3.4 Assembly 設計

![Figure 5: Assembly Plan](figures/fig5_assembly_plan.png)

**表2: 3段階アセンブリ計画サマリー**

| 指標 | 値 |
|------|-----|
| Tier-1 合成ブロック数 | 107 |
| Tier-2 サブ染色体数 | 11 |
| 総合成 DNA 量 | 539,560 bp |
| 平均 Tier-1 サイズ | 4,962 bp |
| 平均 Tier-2 サイズ | 48,272 bp |
| アセンブリオーバーラップ | 40 bp |

### 3.5 JCVI-syn3.0 拡張ケーススタディ

![Figure 7: Case Study Dashboard](figures/fig7_case_study_dashboard.png)

.git .github .gitignore AGENTS.md data figures logs  9,629 bp（全ゲノムの 1.81%）の配列圧縮が可能と推定された。これは主に非必須遺伝子ペアの融合（プロモーター・RBS・遺伝子間領域の節約）による。Menuhin-Gruman et al. (2025) の STABLES 戦略をさらに活用することで、追加の圧縮が期待できる。Tests 

---

## 4. 考察と今後の展望

### 4.1 コドン最適化と反復配列のトレードオフ

# RSCU が高い特定コドン（例: AAAys (K), TTALeu など）を多用することで、同一コドンが連続・近接する構造が生まれるためと解釈される。この問題に対処するためには、（a）コドン選択時に同一コド'REPORT_EOF'
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$b）最適化後の反復スクリーニングと選択的コドン置換、（c）mRNA二次構造エネルギー（ΔG）の同時最小化、が必要である。

### 4.2 機械学習モデルの解釈可能性

#
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             **実データへの適用**: M. mycoides や Mycoplasma genitalium の実 TN-seq データへの適用}
2. **マルチオミクス統合**: トランスクリプトーム・プロテオームデータによる特徴量強化
3. **条件依存的必須性**: 異なる増殖条件下での必須遺伝子セットの変動の組み込み
4. **反復配列除去アルゴリズムの改良**: 2段階最適化（CAI最大化→反復最小化）の実装
5. **全ゲノムシミュレーション**: 提案ゲノム設計の in silico 増殖シミュレーション（FBA/ODEモデル）

---

## 5. 生成したファイル一覧

### ソースコード
| ファイル | 概要 | 行|
|---------|------|------|
| `src/essential_gene_predictor.py` | 必須遺伝子予測 ML モジュール | ~200 |
#| `src/codon_optimizer.py` | コドン最適'REPORT_EOF''REPORT_EOF'
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             } | ~250 |
| `src/genome_architect.py` | ゲノム設計・アセンブリ計画モジュール | ~280 |
| `src/pipeline.py` | 統合オーケストレーションスクリプト | ~360 |

### 図
| ファイル | 内容 |
|---------|------|
| `figures/fig1_model_comparison.png` | ML モデル比較（AUROC & F1, 5-CV） |
| `figures/fig2_feature_importances.png` | RF 特徴量重要度 |
| `figures/fig3_codon_optimization.png` | コドン最適化前後比較（CAI, GC, 反復） |
| `figures/fig4_genome_composition.png` | 機能カテゴリ分布 + オペロンサイズ分布 |
| `figures/fig5_assembly_plan.png` | 3段階階層的アセンブリ計画 |
| `figures/fig6_strand_bias.png` | 複製鎖バイアス解析 |
| `figures/fig7_case_study_dashboard.png` | JCVI-syn3.0 拡張ケーススタディ統合ダッシュボード |

### 結果ファイル
| ファイル | 内容 |
|---------|------|
| `results/model_cv_metrics.csv` | 交差検証メトリクス（3モデル） |
| `results/feature_importances.csv` | RF 特徴量重要度 |
| `results/essential_predictions.csv` | 各遺伝子の必須性予測確率 |
| `results/codon_optimization_results.csv` | コドン最適化前後の指標（60遺伝子） |
| `results/genome_gene_roster.csv` | シミュレーションゲノム遺伝子リスト（473遺伝子） |
| `results/refactoring_plan.csv` | リファクタリング推奨アクション |
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$Tier 1-3） |
| `results/pipeline_summary.json` | パイプライン全体サマリー指標 |

### ログ
| ファイル | 内容 |
|---------|------|
| `logs/process-log.jsonl` | 実行トレース（JSONL形式） |

---

## 参考文献

1. Hutchison CA 3rd, et al. (2016). Design and synthesis of a minimal bacterial genome. *Science*, 351(6280), aad6253. DOI: 10.1126/science.aad6253

2. Pelletier JF, Glass JI, Strychalski EA. (2022). Cellular mechanics during division of a genomically minimal cell. *Trends in Cell Biology*, 32(11), 900–909. DOI: 10.1016/j.tcb.2022.06.009

3. Billmyre RB, et al. (2025). Landscape of essential growth and fluconazole-resistance genes in Cryptococcus neoformans. *PLoS Biology*, 23(5), e3003184. DOI: 10.1371/journal.pbio.3003184

4. Levitan A, et al. (2020). Comparing the utility of in vivo transposon mutagenesis approaches in yeast species to infer gene essentiality. *Current Genetics*, 67, 49–65. DOI: 10.1007/s00294-020-01096-6

5. Menuhin-Gruman I, et al. (2025). AI-directed gene fusing prolongs the evolutionary half-life of synthetic gene circuits. *Science Advances*, 11(40), eadx0796. DOI: 10.1126/sciadv.adx0796

6. Demissie EA, et al. (2025). Comparative analysis of codon optimization tools: advancing toward a multi-criteria framework for synthetic gene design. *Journal of Microbiology and Biotechnology*, 35(4). DOI: 10.4014/jmb.2411.11066

7. Geng PX, et al. (2026). Tripleknock: predicting lethal effect of three-gene knockout in bacteria by deep learning. *Scientific Reports*, 16, 46272. DOI: 10.1038/s41598-026-46272-9

8. Segal ES, et al. (2018). Gene essentiality analyzed by in vivo transposon mutagenesis and machine learning in a stable haploid isolate of Candida albicans. *mBio*, 9(5), e02048-18. DOI: 10.1128/mBio.02048-18

9. Simons A. (2021). Synthetic biology as a technoscience: The case of minimal genomes and essential genes. *Studies in History and Philosophy of Science*, 85, 136–145. DOI: 10.1016/j.shpsa.2020.09.012

10. Sharp PM, Li WH. (1987). The codon Adaptation Index—a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281–1295. DOI: 10.1093/nar/15.3.1281

11. Gibson DG, et al. (2010). Creation of a bacterial cell controlled by a chemically synthesized genome. *Science*, 329(5987), 52–56. DOI: 10.1126/science.1190719

12. Cantore T, Gasperini D, Bevilacqua A. (2025). PRODE recovers essential and context-essential genes through neighborhood-informed scores. *Genome Biology*, 26, 77. DOI: 10.1186/s13059-025-03501-0
