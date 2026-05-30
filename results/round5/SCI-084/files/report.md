# 実験レポート: EpiTransMapper — RNAエピトランスクリプトーム解析パイプライン

---

## 1. 実験目的と背景

### 1.1 研究背景

RNA修飾（エピトランスクリプトーム）は、mRNA安定性・翻訳効率・RNA-タンパク質相互作用を制御する、遺伝子発現調節の新たな層を形成する。特にN6-メチルアデノシン（m6A）はmRNA内で最も豊富な内部修飾であり、転写物の約25%に存在し、METTl3/METTL14/WTAPライタータンパク質複合体によって動的に書き込まれ、FTO/ALKBH5イレーサーによって除去される。YTHDF1-3・YTHDC1-2・IGF2BP1-3などのリーダータンパク質が修飾を認識し、mRNA分解・翻訳促進・核外輸送などの機能的帰結をもたらす。

がんにおけるm6Aエピトランスクリプトームの異常は、急性骨髄性白血病（AML）、肺腺がん（LUAD）、肝細胞がん（HCC）、神経膠腫など20種類以上で報告されており、治療標的としての可能性が注目されている。

本実験では、MeRIP-seq・DART-seq・ナノポア直接RNA-seqデータを統合的に解析できるPythonパイプライン **EpiTransMapper** を設計・実装し、合成データを用いてその性能を評価した。

### 1.2 実験目的

1. MeRIP-seq・DART-seqのm6Aピーク検出アルゴリズムの実装と評価
2. ナノポアイオン電流シグナルに基づく機械学習m6A分類器の開発
3. 癌 vs. 正常における差分m6A修飾解析の実装
4. m6A修飾密度とmRNA安定性・翻訳効率の相関解析
5. ライター・リーダー・イレーサータンパク質発現とm6Aレベルの関連解析
6. 肺腺がんシミュレーションコホートにおけるm6Aエピトランスクリプトームケーススタディ

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 合成データ生成

実験には、現実的な生物学的ノイズを含む合成データを使用した:

| データセット | モデル | パラメータ |
|------------|--------|-----------|
| MeRIP-seq | 負の二項分布カウント | 5,000サイト、15% m6A、CV=0.35、3+3反復 |
| DART-seq | 二項編集率モデル | 背景編集率2%、m6A付近で高編集 |
| ナノポア | ガウスイオン電流シフト | 5 pAシフト（m6A: 75.5 pA vs 未修飾: 80.5 pA）|
| 差分修飾 | ベータ二項率モデル | 2,000遺伝子、12%差分、Δ=0.15-0.45、6反復 |
| WRE発現 | 対数正規発現モデル | 14遺伝子、30癌+30正常サンプル |

### 2.2 MeRIP-seqピーク検出

**アルゴリズム**: フィッシャーの正確確率検定 + ベンジャミニ・ホッホベルク（BH）FDR補正

各サイトに対して以下の2×2分割表を構築:

```
                | 対象サイト | その他全サイト
---------------------------------------------
IP ライブラリ  | X_ip       | N_ip - X_ip
Input ライブラリ| X_in       | N_in - X_in
```

統計的有意性（p_adj < 0.05）およびエンリッチメント（log₂(IP/Input) ≥ log₂(1.5)）の両条件を満たすサイトをm6Aピークとして呼び出す。

### 2.3 DART-seqピーク検出

**アルゴリズム**: 二項検定を用いた編集率のテスト

各サイトの C→U 編集率 $e_i = k_i/n_i$ を背景率 $e_0 = 0.02$ に対してテスト:

```
p_i = P(X ≥ k_i | X ~ Binomial(n_i, e_0))
```

BH補正後 p_adj < 0.05 かつ編集率 > 5% のサイトをピーク呼び出し。

### 2.4 ナノポア機械学習分類器

**特徴量（8次元）**:
- 平均イオン電流（pA）
- ドウェル時間（ms）
- 電流標準偏差
- 5-merコンテキスト特徴（5次元）

**モデル**:
- Random Forest（RF）: 100本の木、max_depth=8
- Gradient Boosting（GB）: 100 estimators、learning_rate=0.05
- Logistic Regression（LR）: C=1.0

評価: 層別5分割交差検証（AUC-ROC ± SD、F1スコア）

### 2.5 差分修飾解析

**アルゴリズム**: ウィルコクソン・マン・ホイットニー検定 + BH補正

遺伝子ごとに癌・正常の修飾率 $r_{g,c}^{(k)}$（修飾リード / 全リード）を比較:

```
H₀: r_g,cancer =ᵈ r_g,normal
```

有意基準: p_adj < 0.05（BH補正）かつ |Δ| > 0.05

### 2.6 機能アノテーション

- mRNA安定性: スピアマン相関（m6Aピーク密度 vs. mRNA半減期）
- 翻訳効率（TE）: スピアマン相関（m6Aピーク密度 vs. リボソームプロファイリングTE）

### 2.7 ライター/リーダー/イレーサー（WRE）解析

14種のWREタンパク質発現（癌 vs. 正常）:
- **ライター**: METTL3, METTL14, WTAP, METTL16
- **リーダー**: YTHDF1-3, YTHDC1-2, IGF2BP1-3
- **イレーサー**: FTO, ALKBH5

マン-ホイットニーU検定 + BH補正; 全体的m6Aレベルとのスピアマン相関も計算。

---

## 3. 主要な結果と数値

### 3.1 パイプライン全体像

![パイプラインアーキテクチャ](rna_modification_pipeline/figures/fig0_pipeline_overview.png)

*図0: EpiTransMapperのモジュール構成。5種の入力データが技術固有の処理モジュールを経て、多重修飾統合レイヤーで合流し、4つの解析モジュールへ分岐する。*

### 3.2 MeRIP-seq / DART-seq ピーク検出結果

![MeRIP-seq解析](rna_modification_pipeline/figures/fig1_merip_analysis.png)

*図1: MeRIP-seqピーク検出結果。(A) 真のm6Aサイトと背景のエンリッチメント分布。(B) 真陽性（緑）・偽陽性（橙）を示すボルケーノプロット。(C) 検出されたピークのゲノム領域分布（3'UTRが最多）。(D) IP vs. Inputカウント散布図。(E) DRACHモチーフ有無の比較。(F) FDR分布。*

**ピーク検出性能（5,000サイト、750真陽性）**:

| 手法 | TP | FP | FN | TN | 適合率 | 再現率 | F1 |
|------|----|----|----|----|--------|--------|-----|
| MeRIP-seq（Fisher検定+BH） | 482 | 76 | 268 | 4,174 | **0.864** | 0.643 | 0.737 |
| DART-seq（二項検定+BH） | 704 | 751 | 46 | 3,499 | 0.484 | **0.939** | 0.639 |

**主な知見**:
- MeRIP-seqは高適合率（0.864）を示す一方、再現率は0.643にとどまる（感度の低いピーク検出の典型例）
- DART-seqは高再現率（0.939）を示すが、偽陽性率が高い（適合率0.484）
- m6Aサイトの75%がDRACHモチーフを保有（背景: 30%）

### 3.3 ナノポア分類器性能

![ナノポア分類](rna_modification_pipeline/figures/fig2_nanopore_classification.png)

*図2: ナノポア直接RNA-seqによるm6A検出。(A) 修飾/未修飾サイトのイオン電流分布。(B) 3分類器のROC曲線（5分割CV）。(C) 交差検証性能比較（平均±SD）。(D) ランダムフォレスト特徴量重要度。*

**5分割交差検証結果（3,000サイト、600真m6A）**:

| 分類器 | AUC-ROC（平均±SD） | F1スコア（平均±SD） |
|--------|-------------------|-------------------|
| Random Forest | 0.753 ± 0.020 | 0.329 ± 0.047 |
| Gradient Boosting | 0.755 ± 0.016 | 0.360 ± 0.043 |
| Logistic Regression | **0.777 ± 0.018** | **0.371 ± 0.027** |

**特徴量重要度（RF）**: 平均イオン電流 > 5-merコンテキスト特徴 > ドウェル時間 > 電流標準偏差

### 3.4 差分修飾解析結果

![差分修飾解析](rna_modification_pipeline/figures/fig3_differential_modification.png)

*図3: 差分m6A修飾（癌 vs. 正常）。(A) MAプロット（緑: 真陽性、橙: 偽陽性）。(B) 差分修飾ボルケーノプロット。(C) FDR閾値に対する精度再現率曲線。(D) 効果量分布。*

**差分修飾解析性能（2,000遺伝子、240真の差分修飾遺伝子）**:

| 解析 | TP | FP | FN | TN | 適合率 | 再現率 | F1 |
|------|----|----|----|----|--------|--------|-----|
| Wilcoxon + BH（6反復） | 222 | 8 | 18 | 1,752 | **0.965** | **0.925** | **0.945** |

検出された差分修飾遺伝子のうち:
- 過剰修飾（Hyper-m6A）: ~120遺伝子（癌で増加）
- 低下修飾（Hypo-m6A）: ~102遺伝子（癌で減少）

### 3.5 機能アノテーション結果

| 機能関連 | スピアマン相関係数（ρ） | p値 | 解釈 |
|----------|------------------------|-----|------|
| m6A密度 vs. mRNA半減期 | **-0.244** | 5.5 × 10⁻¹⁵ | 高m6A → 加速分解（YTHDF2経路） |
| m6A密度 vs. 翻訳効率（TE） | **+0.186** | 3.1 × 10⁻⁹ | 高m6A → 翻訳促進（IGF2BP経路） |

### 3.6 WRE関連解析

**WRE遺伝子差分発現（14遺伝子、30癌 + 30正常）**:
- 有意な差分発現: 8/14遺伝子（57%）
- ライター全4遺伝子（METTL3, METTL14, WTAP, METTL16）が癌で有意高発現
- イレーサー2遺伝子（FTO, ALKBH5）も有意差あり

### 3.7 癌エピトランスクリプトームケーススタディ

![機能解析とWRE](rna_modification_pipeline/figures/fig4_functional_analysis.png)

*図4: 機能アノテーションとWRE解析。(A) m6A密度 vs. mRNA安定性散布図。(B) m6A密度 vs. 翻訳効率散布図。(C) WRE遺伝子発現ボルケーノプロット。(D) WRE-m6A相関棒グラフ。(E) DRACHモチーフシーケンスコンテキスト。(F) パイプライン性能サマリー。*

![癌ケーススタディ](rna_modification_pipeline/figures/fig5_cancer_case_study.png)

*図5: 模擬LUAD（肺腺がん）コホートにおけるm6Aエピトランスクリプトーム異常。(A) 癌における全体的m6Aレベルの有意な上昇（p < 10⁻⁶）。(B) 主要がん遺伝子転写産物のm6Aレベル比較。(C) 全体m6Aレベルによる生存層別化（KMカーブ類似）。(D) WRE発現ヒートマップ（Zスコア）。(E) m6A制御経路のエンリッチメントスコア。(F) 多重RNA修飾共存行列（m6A/m5C/Ψ）。*

**癌ケーススタディ主要結果**:

| 指標 | 癌（n=30） | 正常（n=30） | p値 |
|------|-----------|-------------|-----|
| 全体m6Aレベル（中央値） | 0.62 ± 0.11 | 0.45 ± 0.09 | < 10⁻⁶ |
| EGFR m6A修飾率 | 0.72 | 0.45 | — |
| KRAS m6A修飾率 | 0.81 | 0.50 | — |
| MYC m6A修飾率 | 0.68 | 0.42 | — |

**RNA修飾共存率**:
- m6A × m5C共存: 6.0%（m6Aサイト基準）
- m6A × Ψ共存: 3.7%
- m5C × Ψ共存: 6.4%

---

## 4. 考察と今後の展望

### 4.1 各手法の性能評価と自己批判

**MeRIP-seqピーク検出（F1=0.737）**:
このF1スコアは許容範囲内だが、実際のMeRIP-seqデータでは以下の要因でさらに低下が予想される:
- 抗体の非特異的結合（偽陽性率10-20%増加）
- RNAフラグメンテーションバイアス
- リピート配列・スプライスジャンクションでのマッピングアーチファクト
- 再現性のないピーク（5%未満の発現遺伝子で問題が顕著）

**DART-seqピーク検出（F1=0.639）**:
再現率は高い（0.939）が、偽陽性が多い（適合率0.484）。実際のDARTデータでは、APOBEC1融合タンパク質の発現レベルがオフターゲット編集に大きく影響し、本シミュレーションよりもノイズが高い可能性がある。

**ナノポア分類器（AUC=0.75-0.78）**:
m6AnetのAUC=0.86と比較して低い。差異の要因:
1. 浅い機械学習モデル（Neural networkの方が優れた特徴抽出が可能）
2. シングルリード特徴量（複数リードの集計による精度向上なし）
3. 保守的なシミュレーション設計（5 pAシフト、現実は文脈依存で変動する）
4. クラス不均衡（20%陽性）によるF1の低下

**差分修飾解析（F1=0.945）**:
一見優れた性能だが、この結果には重要な前提条件がある:
- 6反復/条件（実際の研究では3反復が一般的）
- バッチ効果なし（実際のマルチコホート研究では重大な問題）
- 大きな効果量（Δ=0.15-0.45）を仮定（実データでは小さな効果の遺伝子が多数）
実際の研究環境では、F1=0.70-0.85が現実的な期待値。

**機能アノテーション（ρ=-0.244、+0.186）**:
統計的有意性は高いが相関係数は中程度。m6Aの機能は高度にコンテキスト依存（3'UTR vs. CDSで効果が異なる）であり、単純線形モデルでは過度に単純化されている。YTHDF2によるmRNA分解とIGF2BP安定化の相反する効果が相殺され、弱い相関に見える可能性がある。

### 4.2 実世界への一般化可能性

本実験の合成データから実世界データへの一般化には以下の点を考慮する必要がある:

1. **サンプルサイズ**: 本研究の30癌/30正常は中規模。TCGAのような大規模コホート（>200サンプル）での検証が必要
2. **がんの不均一性**: 腫瘍内heterogeneityと間質細胞コンタミネーションがエピトランスクリプトーム信号を希釈する
3. **修飾化学量論**: 本シミュレーションは部位レベルの0/1ラベルを使用するが、実際にはサイトレベルの修飾化学量論（0-100%）が重要
4. **多重修飾相互作用**: m6A・m5C・Ψの機能的相互作用は未解明であり、独立して解析することでエピジェネティクス調節の複雑さが失われる

### 4.3 今後の展望

1. **シングルセルエピトランスクリプトーム**: scMeRIP-seqおよびナノポア単一細胞データへの対応
2. **深層学習統合**: ナノポアシグナル解析へのLSTM/Transformerアーキテクチャ導入
3. **多重修飾同時検出**: m6A/m5C/Ψの同時ML分類器開発
4. **臨床応用**: WRE阻害剤（STM2457/METTL3, CS1/CS2/FTO）の感受性予測モデル
5. **アレル特異的修飾**: SNPとの連鎖解析による機能的バリアント同定
6. **マルチコホート検証**: TCGA・ENCODE・GEO等の公開データを用いた外部検証

---

## 5. 生成したファイル一覧

### ソースコード

| ファイル | 内容 |
|---------|------|
| `rna_modification_pipeline/src/pipeline.py` | メインパイプライン（~600行） |

### 図表

| ファイル | 内容 |
|---------|------|
| `rna_modification_pipeline/figures/fig0_pipeline_overview.png` | パイプラインアーキテクチャ概要 |
| `rna_modification_pipeline/figures/fig1_merip_analysis.png` | MeRIP-seq解析結果（6パネル） |
| `rna_modification_pipeline/figures/fig2_nanopore_classification.png` | ナノポア分類器性能（4パネル） |
| `rna_modification_pipeline/figures/fig3_differential_modification.png` | 差分修飾解析（4パネル） |
| `rna_modification_pipeline/figures/fig4_functional_analysis.png` | 機能アノテーション・WRE（6パネル） |
| `rna_modification_pipeline/figures/fig5_cancer_case_study.png` | 癌ケーススタディ（6パネル） |

### レポート

| ファイル | 内容 |
|---------|------|
| `paper.md` | 学術論文形式（英語、Abstract/Introduction/Methods/Results/Discussion/Conclusion/References） |
| `report.md` | 本実験レポート（日本語、全結果・手法・考察） |

---

## 6. 参考文献

1. Dominissini D et al. (2012). Topology of the human and mouse m6A RNA methylomes revealed by m6A-seq. *Nature*, 485:201-206. DOI: 10.1038/nature11112

2. Meyer KD et al. (2012). Comprehensive analysis of mRNA methylation reveals enrichment in 3' UTRs and near stop codons. *Cell*, 149:1635-1646. DOI: 10.1016/j.cell.2012.05.003

3. Petri BJ et al. (2023). m6A readers, writers, erasers, and the m6A epitranscriptome in breast cancer. *J Mol Endocrinol*, 70(2). DOI: 10.1530/JME-22-0110

4. Wang W et al. (2022). FTO promotes Bortezomib resistance via m6A-dependent destabilization of SOD2 expression in multiple myeloma. *Cancer Gene Therapy*. DOI: 10.1038/s41417-022-00429-6

5. Meyer KD (2019). DART-seq: an antibody-free method for global m6A detection. *Nature Methods*, 16:1275-1280. DOI: 10.1038/s41592-019-0570-0

6. Hendra C et al. (2022). Detection of m6A from direct RNA sequencing using a multiple instance learning framework. *Nature Methods*, 19:1590-1598. DOI: 10.1038/s41592-022-01666-1

7. Pratanwanich PN et al. (2021). Identification of differential RNA modifications from nanopore direct RNA sequencing with xPore. *Nature Biotechnology*, 39:1394-1402. DOI: 10.1038/s41587-021-00949-w

8. Cheng G et al. (2026). Raw signal segmentation for estimating RNA modification from Nanopore direct RNA sequencing data. *eLife*, 14:e104618. DOI: 10.7554/elife.104618

9. Li X et al. (2025). NSUN2-mediated HCV RNA m5C Methylation Facilitates Viral RNA Stability and Replication. *Genomics, Proteomics & Bioinformatics*. DOI: 10.1093/gpbjnl/qzaf008

10. Zhang W, Pan T (2022). Pseudouridine RNA modification detection and quantification by RT-PCR. *Methods*, 203:1-4. DOI: 10.1016/j.ymeth.2021.05.010

11. Ge R et al. (2023). m6A-SAC-seq for quantitative whole transcriptome m6A profiling. *Nature Protocols*. DOI: 10.1038/s41596-023-00862-3

12. Zhong Z et al. (2023). Systematic comparison of tools used for m6A mapping from nanopore direct RNA sequencing. *Nature Communications*, 14:3714. DOI: 10.1038/s41467-023-37596-5

---

*EpiTransMapper v1.0 | 解析日: 2026年5月29日 | Python 3.11 | ランダムシード: 42*
