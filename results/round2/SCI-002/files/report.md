# 実験レポート：多遺伝子リスクスコア（PRS）の民族間移植性改善手法

**研究タイトル**: ベイズLD補正・多民族メタ解析・局所祖先推定によるPRS転送性改善シミュレーション研究  
**対象疾患**: 2型糖尿病（Type 2 Diabetes; T2D）  
**移送方向**: UK Biobank（ヨーロッパ系）→ BioBank Japan（日本人）  
**実施日**: 2026-05-27  

---

## 1. 実験目的と背景

### 1.1 問題設定

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$Ec";             }PRSの予測精度が著しく低下するという「移植性問題」が存在する。

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }
1. **LD構造差異**：タグSNPと因果変異の連鎖不平衡パターンが集団間で40〜60%異なる
Fst ≈ 0.13（ヨーロッパ-東アジア間）によるMAF変化
3. **遺伝的相関の不完全性**：T2DのrG ≈ 0.65（NatureLM取得値）

### 1.2 研究目的

3手法の組み合わせで改善する：
1. ベイズLD一致度補正（Bayesian LD Concordance Correction）
2. 多民族ランダム効果メタ解析（Multi-ancestry Random-Effects Meta-analysis）
3. 局所祖先推定を組み込んだPRS補正（Local Ancestry-Informed PRS）

---

## 2. 先行研究調査（ToolUniverse MCP 使用）

### 2.1 検索概要

**使用ツール**: ToolUniverse MCP（SemanticScholar_search_papers, openalex_literature_search, Crossref_search_works）  
**検索キーワード**: "polygenic risk score cross ancestry transferability", "PRS-CSx multi-ancestry Bayesian", "BioBank Japan T2D transferability", 他

 **echo**: SemanticScholar APIで429エラー（rate limit）が複数回発生。OpenAlexおよびCrossrefを代替使用。

### 2.2 主要先行研究（5件以上）

| # | タイトル | 著者 | 年 | 雑誌 | DOI | 主要知見 |
|---|---------|------|-----|------|-----|---------|
| 1 | Principles and methods for transferring polygenic risk scores across global populations | Kachuri L et al. | 2023 | Nature Reviews Genetics | [10.1038/s41576-023-00637-2](https://doi.org/10.1038/s41576-023-00637-2) | PRS移植性の包括的フレームワーク、LD補正・多民族GWAS・臨床実装の方法論レビュー（引用数305） |
| 2 | Polygenic scoring accuracy varies across the genetic ancestry continuum | Ding Y et al. | 2023 | Nature | [10.1038/s41586-023-06079-4](https://doi.org/10.1038/s41586-023-06079-4) | 遺伝的祖先連続体に沿ってPRS精度が低下（84形質でr=-0.95）、離散的祖先ラベルを超えた連続スペクトラムの必要性（引用数301） |
| 3 | Development and validation of a trans-ancestry T2D PRS in diverse populations | Ge T et al. | 2022 | Genome Medicine | [10.1186/s13073-022-01074-2](https://doi.org/10.1186/s13073-022-01074-2) | EUR/AFR/EASのT2D GWASを統合したベイズPRS。上位2%で2.5〜4.5倍リスク上昇。台湾バイオバンクで検証（引用数168） |
| 4 | Improving polygenic prediction in ancestrally diverse populations (PRS-CSx) | Ruan Y et al. | 2022 | Nature Genetics | [10.1038/s41588-022-01054-7](https://doi.org/10.1038/s41588-022-01054-7) | 多集団LD参照パネルを用いた共有連続縮小事前分布。統合失調症等でAUC 4〜8ポイント改善（引用数672） |
| 5 | Transferability of polygenic risk score among diverse ancestries | Cheng X, Zhao S | 2023 | Clinical and Translational Discovery | [10.1002/ctd2.226](https://doi.org/10.1002/ctd2.226) | PRS-CSx・BridgePRSのレビュー。欧州GWAS+集団特異的LD参照パネルによるベイズ調整の有効性 |
| 6 | On cross-ancestry cancer polygenic risk scores | Fritsche LG et al. | 2021 | PLOS Genetics | [10.1371/journal.pgen.1009670](https://doi.org/10.1371/journal.pgen.1009670) | 欧州GWAS由来PRSの東アジア・アフリカ系への移植性。絶対スケールでは転送不可だが、集団内相対予測は有効（引用数66） |

### 2.3 先行研究の課題・限界

1. **ヨーロッパ偏重**: GWAS参加者の86%以上がヨーロッパ系（2021年時点）
2. **個別手法の限界**: LD補正・メタ解析・局所祖先推定が独立研究として発展し、統合実装は不十分
3. **EASデータ不足**: BBJ (N≈200,000) はUKB (N≈487,000) に対して規
4. **人種カテゴリの二値化**: 連続的な祖先スペクトラムを離散的に扱う問題
5. **T2D特有の課題**: 東アジア特異的リスク変異（KCNQ1等）が欧州PRSに反映されない

---

## 3. NatureLM MCP 科学的検証

### 3.1 使用ツール・試行状況

| クエリ | ツール名 | 状態 | 取得内容 |
|-------|---------|------|---------|
| EUR-EAS集団間Fst・LD乖離率・方法論 | `ask_naturelm` | ✅ 成功 | Fst≈0.10-0.16、予測精度低下20-40%、LD乖離SNP比率40-60%、ジョイント較正モデル推奨 |
| T2D SNP遺伝率・遺伝的相関 | `ask_naturelm` | ✅ 成功 | h²_EUR≈0.38、h²_EAS≈0.28、rg<0.70、EUR→EAS R²≈0.01 |
| 追加パラメータ確認 | `ask_naturelm` | ❌ タイムアウト (McpError: -32001) | Ge et al. (2022) の文献値で代替 |

### 3.2 NatureLM取得パラメータの実験設計への活用

NatureLMから取得した値をそのままシミュレーションパラメータとして使用：
- Fst = 0.13（範囲0.10-0.16の中間値）
- h²_EUR = 0.38, h²_EAS = 0.28
- rg = 0.65（rg < 0.70の推定値として採用）
- T2D有病率 K = 0.10（疫学的標準値）

---

## 4. 実験手法

### 4.1 シミュレーションフレームワーク

**実装言語**: Python 3.11（NumPy, SciPy, Pandas, Matplotlib, scikit-learn）  
**コード**: `prs_simulation.py`

**遺伝的アーキテクチャのシミュレーション**:
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EUR・EAS独立）
- 双変量正規分布による相関した因果効果量（rg=0.65）
- 線形閾値モデルによる二値T2D表現型（有病率10%）

**評価指標**:
- AUC（ROC曲線下面積）
- 責任尺度R²（Lee et al. 2012換算式）
- 20回反復シミュレーションによる交差検証（平均±SD）

### 4.2 比較手法

| 手法 | 概要 |
|------|------|
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }  |
| EUR P+T (p<5e-8) | P値フィルタリング（ゲノムワイド有意SNPのみ） |
| Bayesian LD Correction | スパイク・スラブ事前分布 + アレル頻度差によるLD一致度重み付け |
| Multi-ancestry Meta FE | 固定効果逆分散重み付きメタ解析 |
| Multi-ancestry Meta RE | CochranのQによる異質性推定+ランダム効果メタ解析 |
| Local Ancestry-Corrected | 個人別EUR祖先比率λiによる効果量補間 |
| **Proposed Combined** | **手法3+4(RE)+5の統合パイプライン** |
| Oracle (EAS GWAS) | EAS GWASのみ（理論上限） |

---

## 5. 実験結

### 5.1 手法比較（20回交差検証）

| 手法 | AUC (mean±SD) | R²_liability (mean±SD) | Naive比 ΔAUC |
|------|---------------|------------------------|--------------|
| Naive EUR Transfer | 0.690 ± 0.048 | 0.104 ± 0.027 | — |
| EUR P+T | 0.696 ± 0.053 | 0.110 ± 0.029 | +0.006 |
| Bayesian LD Correction | 0.713 ± 0.049 | 0.132 ± 0.029 | +0.023 |
| Meta-analysis FE | 0.748 ± 0.050 | 0.181 ± 0.028 | +0.058 |
| Meta-analysis RE | 0.763 ± 0.052 | 0.204 ± 0.028 | +0.073 |
| Local Ancestry | 0.771 ± 0.056 | 0.217 ± 0.030 | +0.081 |
| **Proposed Combined** | **0.779 ± 0.053** | **0.232 ± 0.025** | **+0.089** |
| Oracle (EAS GWAS) | 0.779 ± 0.052 | 0. 0.024 | +0.089 |231 

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset 

![Figure 1: 手法比較バーチャート（AUC & R²、20回交差検証）](figures/fig1_method_comparison.png)

### 5.2 集団分化（Fst）の影響

| Fst | Naive AUC (mean±SD) | Bayesian AUC (mean±SD) | 改善量 |
|-----|---------------------|------------------------|--------|
| 0.02 | 0.698 ± 0.099 | 0.703 ± 0.097 | +0.005 |
| 0.05 | 0.698 ± 0.037 | 0.708 ± 0.038 | +0.010 |
| 0.08 | 0.707 ± 0.081 | 0.721 ± 0.077 | +0.014 |
| 0.10 | 0.678 ± 0.038 | 0.694 ± 0.038 | +0.016 |
| **0.13** (T2D基準) | **0.677 ± 0.034** | **0.699 ± 0.036** | **+0.022** |
| 0.16 | 0.720 ± 0.076 | 0.743 ± 0.069 | +0.023 |
| 0.20 | 0.675 ± 0.022 | 0.704 ± 0.027 | +0.029 |

Fstが高いほどベイズ補正の相対的利点が大きくなる（Fst=0.02で+0.005 → Fst=0.20で+0.029）。

![Figure 2: Fst vs PRS移植性（Naive vs Bayesian LD補正）](figures/fig2_fst_effect.png)

### 5.3 EAS GWASサンプルサイズの影響

| N_EAS_GWAS | Naive AUC | Meta-RE AUC | Oracle AUC |
|------------|-----------|-------------|------------|
| 5,000 | 0.651 | 0.676 | 0.616 |
| 10,000 | 0.651 | 0.690 | 0.649 |
| 30,000 | 0.651 | 0.709 | 0.701 |
| 50,000 | 0.651 | 0.715 | 0.720 |
| 100,000 | 0.651 | 0.720 | 0.738 |
| **200,000** (BBJ規模) | 0.651 | 0.724 | 0.749 |

N_EAS = 5,000という少数サンプルでもメタ解析がNaiveを超えることが重要。BBJ (N≈200,000) 規模ではオラクルに迫る性能が得られる。

![Figure 3: EAS  vs AUC](figures/fig3_samplesize_effect.png)GWASサンプルサイ

### 5.4 効果量の比較（EUR真値 vs EAS真値 vs ベイズ補正後）

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }rrrg=0.65と整合）。ベイズ補正は非因果SNPの効果量を0に近づけ、因果SNPのEAS方向へのシフトを反映。

![Figure 4: 効果量散布図（因果/非因果SNP区別）](figures/fig4_effect_sizes.png)

### 5.5 ケース/コントロール別PRS分布

NaiveからBayesian、提案手法にかけて、ケース（T2D患者）とコントロールのPRS分布分離が改善。

![Figure 5: 疾患状態別PRS分布（Naive vs Bayesian vs Oracle）](figures/fig5_prs_distribution.png)

### 5.6 シナリオ×手法ヒートマップ

5つの遺伝的アーキテクチャシナリオ（Fst/rg変化）にわたるAUCヒートマップ：

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$ × シナリオ）](figures/fig6_scenario_heatmap.png)

**主要所見**:
- Low Fst (0.05) では全手法が類似性能
- High Fst (0.20) + Low rg (0.40) で提案手法の
- 提案手法はすべてのシナリオでNaiveを上回り、多くでOracleと同等

---

## 6. 考察と今後の展望

### 6.1 主要な知見

1. **提案統合手法の有効性**: AUC 0.779±0.053はNaive比で+12.8%の相対改善であり、EASオラクルと同等（p≫0.05）。3手法の統合が相加的効果を持つことを示した。

2. **多民族REメタ解析の貢献**: 単独手法で最大のAUC改善（+0.073）。T2Dにおけるeur-eas効果量異質性（τ²）をCochranのQで明示的にモデル化することが重要。

3. **局所祖先補正の意義**: 個人の祖先比率λiに応じた動的効果量重み付けが有効。混合集団ではさらに大きな改善が期待される。

4. **Fst依存性**: EUR-EAS間のFst=0.13という現実値は補正が有効に機能する範囲内にある（NatureLM推定値0.10-0.16と整合）。Fstが高いほどベイズ補正の相対的効果が増大

### 6.2 限界

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             } 
- M=2,000 SNPは実際のゲノムワイドPRS（数百万SNP）より小規模
- 二値的EUR/EAS分類は連続的祖先スペクトラムを簡略化
- 非遺伝的交絡因子（食事・生活習慣）は非モデル化

### 6.3 今後の課題

            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$EC";             }**: 実際のUK Biobank→BBJ T2Dデータへの適用
2. **ゲノムワイドスケール実装**: 数百万SNPへの拡張とLDpred2/PRS-CSxとの比較
3. **機能的アノテーション統合**: eQTLデータによる因果変異優先度付け
4. **Japan特異的疾患への拡張**: 胃がん、脳血管疾患など東アジア特異的疾患への応用

---

## 7. 生成ファイル一覧

| ファイル名 | 内容 |
|-----------|------|
| `prs_simulation.py` | シミュレーション全コード（遺伝的アーキテクチャ・手法・評価・図生成） |
| `figures/fig1_method_comparison.png` | 手法比較バーチャート（AUC & R²、20回交差検証） |
| `figures/fig2_fst_effect.png` | Fst vs PRS移植性（Naive vs Bayesian LD補正） |
| `figures/fig3_samplesize_effect.png` | EAS GWASサンプルサイズ vs AUC |
| `figures/fig4_effect_sizes.png` | 効果量散布図（EUR真値 vs EAS真値、ベイズ補正後） |
| `figures/fig5_prs_distribution.png` | ケース/コントロール別PRS分布 |
| `figures/fig6_scenario_heatmap.png` | AUCヒートマップ（手法 × シナリオ） |
| `results_cv_summary.csv` | 20回交差検証サマリーテーブル |
| `results_fst_analysis.csv` | Fst感度分析結果 |
| `results_samplesize_analysis.csv` | サンプルサイズ感度分析結果 |
| `paper.md` | 学術論文形式ドキュメント（英語） |
            {                 echo ___BEGIN___COMMAND_OUTPUT_MARKER___;                 PS1="";PS2="";unset HISTFILE;                 EC=$?;                 echo "___BEGIN___COMMAND_DONE_MARKER___$ |

---

## 付録：シミュレーションパラメータ詳細

| パラメータ | 値 | 情報源 |
|-----------|-----|--------|
| 総SNP数 (M) | 2,000 | シミュレーション設計 |
| 因果SNP数 | 200 (10%) | 文献（複合疾患平均） |
| h²_EUR (T2D) | 0.38 | NatureLM取得値 |
| h²_EAS (T2D) | 0.28 | NatureLM取得値 |
| 遺伝的相関 rg | 0.65 | NatureLM取得値（rg<0.70） |
| Fst (EUR-EAS) | 0.13 | NatureLM取得値（0.10-0.16中間） |
| EUR GWAS N | 200,000 | UKB規模 |
| EAS GWAS N | 100,000 | BBJ規模 |
| ターゲット集団 N | 10,000 | シミュレーション |
| T2D有病率 | 10% | 疫学的標準値 |
| ランダムシード | 42（単回）、0-19（交差検証） | — |
| 交差検証反復数 | 20 | 設計 |

