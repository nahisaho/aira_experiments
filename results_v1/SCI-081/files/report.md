# がんプロテオゲノミクス統合解析パイプライン — 設計レポート

**DRAFT — NOT FOR DISTRIBUTION**

| 項目 | 内容 |
|------|------|
| 日付 | 2026-05-23 |
| 対象コホート | CPTAC 膵臓がん (PDAC): 腫瘍 140例 + 正常隣接組織 67例 |
| ゲノムビルド | GRCh38 / UniProt Human 2024_01 |
| 解析環境 | MaxQuant 2.4.13 + Perseus 2.1.2 + R 4.3.2 / Bioconductor 3.18 |

---

## 1. 実験目的と背景

### 1.1 研究目的

がんプロテオゲノミクスは、ゲノム・トランスクリプトーム・プロテオームを統合的に解析し、がんの分子機構をタンパク質レベルで解明するアプローチである。本パイプラインは以下の6つの解析モジュールを統合し、CPTAC膵管腺癌(PDAC)コホートへの適用を通じて、臨床的に実行可能な知見を導出することを目的とする。

### 1.2 背景

- **プロテオゲノミクスの意義**: ゲノム変異の約40%はタンパク質レベルでは検出されず（翻訳制御・タンパク質分解による）、マルチオミクス統合が不可欠
- **CPTAC**: Clinical Proteomic Tumor Analysis Consortium。NCI主導の大規模プロテオゲノミクスプロジェクトで、膵臓がんを含む複数がん種のデータを公開
- **PDAC**: 5年生存率 ~12%。分子サブタイプ（Basal-like / Classical）により予後・治療応答が異なるが、タンパク質レベルでの層別化は十分に確立されていない

---

## 2. 使用した手法・アルゴリズムの概要

### Module 1: Variant Peptide Search（ゲノム変異のプロテオーム検索反映）

| 要素 | 詳細 |
|------|------|
| **入力** | 体細胞変異 VCF、参照プロテオーム (UniProt)、MS/MS スペクトル |
| **手法** | customProDB (R/Bioconductor) で変異タンパク質配列を生成 → 参照DBと結合 → MaxQuant検索 |
| **FDR制御** | Target-decoy方式、PSM/Protein level ともに 1% FDR |
| **変異タイプ** | ミスセンス、フレームシフト、スプライスジャンクション、遺伝子融合 |
| **後処理** | Andromeda score ≥ 40、PEP < 0.01、COSMIC再発性アノテーション |
| **スクリプト** | `scripts/01_variant_peptide_search.R` |

**アルゴリズムの流れ:**
1. VCFファイルからPASS変異を抽出
2. Ensembl TxDb でコーディング変異をアノテーション
3. customProDB `OutputVarproseq()` で変異タンパク質FASTAを生成
4. 参照プロテオーム + 変異DB + コンタミナントDB を結合
5. MaxQuant XML を自動生成し検索実行
6. 変異ペプチド PSM を品質フィルタリング

### Module 2: RNA–Protein Discordance（翻訳制御推定）

| 要素 | 詳細 |
|------|------|
| **入力** | RNA-seq TPM マトリクス、MaxQuant LFQ 強度 (proteinGroups.txt) |
| **正規化** | メディアンセンタリング + log2変換 |
| **相関解析** | 遺伝子ごとの Spearman ρ（RNA vs タンパク質、サンプル横断） |
| **翻訳効率推定** | 回帰残差法: lm(Protein ~ RNA) の残差 = RNA で説明できないタンパク質量変動 |
| **乖離判定** | |Δz-score| > 1.5 を乖離遺伝子と判定 |
| **パスウェイ解析** | GO Biological Process 濃縮解析 (clusterProfiler) |
| **スクリプト** | `scripts/02_rna_protein_discordance.R` |

**期待される知見（PDAContext）:**
- 中央値 Spearman ρ ≈ 0.3–0.5（RNA-タンパク質相関は中程度）
- EMT関連遺伝子（VIM, CDH1）で強い翻訳制御
- リボソーム関連遺伝子は高い正の相関
- プロテアソーム・オートファジー経路は RNA > Protein（分解亢進）

### Module 3: Phosphoproteomics & Kinase Activity（リン酸化プロテオミクス）

| 要素 | 詳細 |
|------|------|
| **入力** | MaxQuant Phospho(STY)Sites.txt、proteinGroups.txt |
| **前処理** | Class I サイトフィルタ (局在化確率 ≥ 0.75)、タンパク質レベル正規化 |
| **欠損値補完** | MinProb 法（Perseus方式）: 左裾正規分布からランダムサンプリング (shift=1.8σ, width=0.3σ) |
| **差次解析** | limma moderated t-test (Tumor vs Normal)、adj.P < 0.05 & |log2FC| > 1.0 |
| **キナーゼ活性** | KSEA (Kinase-Substrate Enrichment Analysis): PhosphoSitePlus基質DBに基づく |
| **代替手法** | PhosR `kinaseSubstrateScore()` / `kinaseSubstratePred()` |
| **スクリプト** | `scripts/03_phosphoproteomics_kinase.R` |

**KSEA アルゴリズム:**
1. 各キナーゼの既知基質をPhosphoSitePlusから取得
2. 基質サイトのlog2FCを集約し、z-score = (mean_sub - global_mean) / (global_sd / √n)
3. 両側検定でp値計算、BH法で多重検定補正
4. 最低基質数 ≥ 3 のキナーゼのみ報告

### Module 4: Neoantigen Verification（ネオアンチゲン候補検証）

| 要素 | 詳細 |
|------|------|
| **入力** | 変異ペプチド (Module 1)、HLA型 (OptiType)、RNA発現量 |
| **結合予測** | NetMHCpan-4.1: MHC-I結合親和性予測 (8-11mer) |
| **フィルタ** | 結合親和性 < 500 nM or %Rank < 2.0% |
| **発現フィルタ** | RNA TPM ≥ 1.0 & MS検出あり（プロテオミクス検証） |
| **免疫原性** | PRIME スコア ≥ 0.5 |
| **スクリプト** | `scripts/04_neoantigen_verification.py` |

**多層フィルタリング戦略:**
```
全変異ペプチド → HLA結合予測 → 強結合体 → RNA発現確認 → MS検出確認 → 免疫原性スコア
  (~10,000)      (~2,000)      (~500)       (~200)         (~50)         (~20)
```

### Module 5: MOFA+ Multi-Omics Factor Analysis（マルチオミクス因子分解）

| 要素 | 詳細 |
|------|------|
| **入力** | 4 views: mRNA (5000遺伝子) + Proteome (5000) + Phosphoproteome (3000) + CNA |
| **モデル** | MOFA2 (R/Python): Group Factor Analysis の拡張 |
| **因子数** | 15（drop_factor_threshold = 0.01 で自動削減） |
| **収束** | slow モード、seed = 42 |
| **患者層別化** | k-means (k=3) on factor space |
| **生存解析** | Kaplan-Meier + log-rank検定 |
| **臨床関連** | Factor–clinical association (Kruskal-Wallis / Spearman) |
| **スクリプト** | `scripts/05_mofa_integration.R` |

**MOFA+ の統計的枠組み:**
- 確率的生成モデル: X_m = Z × W_m + ε_m （m = 各オミクスview）
- Z: 潜在因子行列 (N_samples × K_factors)
- W_m: view m の重み行列 (D_m × K_factors)
- ARD (Automatic Relevance Determination) で因子の各viewへの寄与度を推定

### Module 6: CPTAC PDAC Case Study

| 要素 | 詳細 |
|------|------|
| **コホート** | CPTAC PDAC: 140腫瘍 + 67正常隣接組織 |
| **プロテオミクス** | TMT-11plex, Orbitrap Fusion Lumos |
| **サブタイプ分類** | タンパク質ベースのBasal/Classical分類（マーカースコアリング） |
| **統合解析** | Module 1–5 の全出力を集約、サブタイプ別に比較 |
| **スクリプト** | `scripts/06_cptac_pdac_casestudy.R` |

---

## 3. 主要な結果と数値

### 3.1 パイプライン設計の主要パラメータ

| パラメータ | 値 | 根拠 |
|-----------|-----|------|
| PSM/Protein FDR | 1% | CPTAC標準 |
| Class I リン酸化サイト | 局在化確率 ≥ 0.75 | Olsen et al. (2006) |
| MinProb imputation shift | 1.8σ | Perseus デフォルト |
| RNA-Protein 乖離閾値 | |Δz| > 1.5 | Mertins et al. (2016) |
| HLA 結合親和性カットオフ | < 500 nM | NetMHCpan 推奨 |
| MOFA 因子数 | 15 (ARDで自動削減) | MOFA2 best practice |
| KSEA 最小基質数 | ≥ 3 | Casado et al. (2013) |

### 3.2 期待される結果概要（CPTAC PDAC文献ベース）

| 解析 | 期待される結果 |
|------|----------------|
| **変異ペプチド検索** | 数百のMS検出ミスセンス変異ペプチド。KRAS G12D/V, TP53 hotspot変異がプロテオームレベルで確認 |
| **RNA-Protein乖離** | 中央値 Spearman ρ ≈ 0.4。EMT・タンパク質分解経路で顕著な乖離。約500遺伝子が強い翻訳制御下 |
| **リン酸化/キナーゼ** | Basal型でSRC/FAK/EGFR活性化。Classical型でPKA/AMPK活性化。CDK-RB軸のリン酸化がBasal型で亢進 |
| **ネオアンチゲン** | KRAS/TP53変異由来のネオアンチゲン候補を同定。約20–50のMS検証済み候補 |
| **MOFA+ 層別化** | 3クラスター: Basal-immune / Classical / Mixed。OS差異あり (log-rank p < 0.01) |

### 3.3 分子サブタイプ特性（PDAC）

#### Basal-like サブタイプ
- **上方制御パスウェイ**: EMT、扁平上皮分化プログラム、NF-κB炎症シグナル
- **活性化キナーゼ**: SRC, FAK, AXL, EGFR, MET
- **マーカー**: KRT5, KRT17, S100A2, TP63
- **予後**: 不良（OS 中央値 ~18ヶ月）

#### Classical サブタイプ
- **上方制御パスウェイ**: 膵分泌プログラム、脂質代謝/PPAR、補体カスケード
- **活性化キナーゼ**: PKA, PKC, AMPK, CK2, CDK4/6
- **マーカー**: GATA6, HNF1A, TFF1, AGR2
- **予後**: 良好（OS 中央値 ~30ヶ月）

---

## 4. 考察と今後の展望

### 4.1 パイプライン設計の意義

本パイプラインは、がんプロテオゲノミクスの主要な解析軸を体系的にカバーする設計となっている。特に以下の点で既存手法を統合・拡張している：

1. **変異-プロテオーム連携**: customProDBによるサンプル特異的データベース構築は、標準プロテオームデータベースでは検出できない変異ペプチドの同定を可能にする
2. **翻訳制御の定量化**: 回帰残差法による翻訳効率推定は、RNA発現量のみでは説明できないプロテオーム変動を体系的に評価する
3. **キナーゼ活性の推定**: KSEAはリン酸化サイトの集約的変動からキナーゼ活性を推定し、薬物標的の優先順位付けに直結する
4. **ネオアンチゲンのMS検証**: 免疫療法標的の選定において、MS/MSによるタンパク質レベルでの確認は偽陽性を大幅に削減する
5. **MOFA+統合**: 4層オミクスデータの同時因子分解により、各オミクス層の独立・共有情報を分離し、より堅牢な患者層別化を実現する

### 4.2 技術的考慮事項

- **FDR制御**: 変異ペプチド検索では検索空間が拡大するため、厳格なFDR制御（PEP < 0.01 + Andromeda score ≥ 40）が必要
- **バッチ効果**: TMTプレックス間のバッチ効果はComBat/limma `removeBatchEffect()` で補正すべき
- **欠損値**: MinProb imputationはMAR (Missing At Random) を仮定しないため、プロテオミクスの検出限界以下の欠損に適している
- **複数検定**: KSEA, 差次リン酸化, GO enrichment の全てにBH法を適用

### 4.3 制限事項

- 本パイプラインはバルク組織解析を前提としており、腫瘍微小環境の細胞組成は反映しない
- ネオアンチゲン予測はMHC-Iのみ（MHC-II結合は今後の拡張対象）
- MOFA+の因子解釈は事後的であり、因子の生物学的意味付けには注意が必要

### 4.4 今後の展望

1. **シングルセルプロテオミクス (SCoPE2/plexDIA)** との統合により、腫瘍不均一性の解像度を向上
2. **空間プロテオミクス (MALDI-MSI)** との統合により、組織内局在情報を追加
3. **MHC-II ネオアンチゲン** 予測の追加（CD4+ T細胞応答の評価）
4. **薬剤感受性データ (GDSC/CCLE)** との統合によるバイオマーカー探索
5. **深層学習ベースのネオアンチゲン予測** (DeepImmuno, PRIME v2) への移行
6. **リアルタイムパイプライン化** (Nextflow/Snakemake) による自動化

---

## 5. 生成したファイル一覧

### 設定ファイル
| ファイル | 説明 |
|---------|------|
| `config/pipeline_config.yaml` | パイプライン全体の設定ファイル（全6モジュール + 計算環境定義） |

### 解析スクリプト
| ファイル | 説明 |
|---------|------|
| `scripts/00_run_pipeline.sh` | パイプラインマスター実行スクリプト（全モジュール順次実行） |
| `scripts/01_variant_peptide_search.R` | Module 1: 変異タンパク質DB構築 + MaxQuant検索設定 + 変異ペプチドフィルタリング |
| `scripts/02_rna_protein_discordance.R` | Module 2: RNA-タンパク質相関解析 + 翻訳効率推定 + 乖離遺伝子同定 |
| `scripts/03_phosphoproteomics_kinase.R` | Module 3: リン酸化サイト解析 + KSEA キナーゼ活性推定 + PhosR |
| `scripts/04_neoantigen_verification.py` | Module 4: ネオアンチゲン候補生成 + HLA結合予測 + MS検証 + 免疫原性スコア |
| `scripts/05_mofa_integration.R` | Module 5: MOFA+ 4層オミクス因子分解 + 患者クラスタリング + 生存解析 |
| `scripts/06_cptac_pdac_casestudy.R` | Module 6: CPTAC PDAC ケーススタディ（サブタイプ分類 + 統合サマリ） |

### 出力ファイル（パイプライン実行後に生成）
| ファイル | 説明 |
|---------|------|
| `results/variant_db/combined_search_db.fasta` | 変異+参照結合プロテオームDB |
| `results/variant_peptides_filtered.tsv` | MS検証済み変異ペプチド一覧 |
| `results/rna_protein_correlations.csv` | 遺伝子ごとのRNA-タンパク質相関 |
| `results/translational_efficiency_scores.csv` | 翻訳効率スコア |
| `results/discordant_genes.csv` | RNA-タンパク質乖離遺伝子 |
| `results/differential_phosphosites.csv` | 差次リン酸化サイト |
| `results/ksea_kinase_scores.csv` | KSEAキナーゼ活性スコア |
| `results/neoantigen_candidates.tsv` | ランク付きネオアンチゲン候補 |
| `results/mofa_model.hdf5` | MOFA+ 学習済みモデル |
| `results/mofa_patient_clusters.csv` | 患者クラスター割り当て |
| `results/mofa_variance_per_factor.csv` | 因子別分散説明率 |
| `results/mofa_factor_clinical_assoc.csv` | 因子–臨床変数関連 |
| `results/integrative_summary.csv` | 統合サマリテーブル |

### 図表（パイプライン実行後に生成）
| ファイル | 説明 |
|---------|------|
| `figures/pipeline_overview.txt` | パイプライン構成図（テキストベース） |
| `figures/rna_protein_correlation_hist.pdf` | RNA-タンパク質相関分布ヒストグラム |
| `figures/rna_protein_discordance_scatter.pdf` | RNA-タンパク質乖離散布図 |
| `figures/phospho_volcano.pdf` | 差次リン酸化ボルケーノプロット |
| `figures/ksea_barplot.pdf` | KSEAキナーゼ活性バープロット |
| `figures/mofa_variance_heatmap.pdf` | MOFA分散説明率ヒートマップ |
| `figures/mofa_factor_scatter.pdf` | MOFA因子空間患者散布図 |
| `figures/mofa_survival_km.pdf` | Kaplan-Meier生存曲線 |
| `figures/mofa_top_weights.pdf` | MOFA因子トップ重みプロット |
| `figures/pdac_subtype_scatter.pdf` | PDAC分子サブタイプ分類散布図 |

### ログ
| ファイル | 説明 |
|---------|------|
| `logs/process-log.jsonl` | パイプライン実行トレース（JSONL形式） |

---

## 6. 参考文献

1. Cao L, et al. "Proteogenomic characterization of pancreatic ductal adenocarcinoma." *Cell* (2021). CPTAC PDAC study.
2. Mertins P, et al. "Proteogenomics connects somatic mutations to signalling in breast cancer." *Nature* (2016).
3. Argelaguet R, et al. "MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data." *Genome Biology* (2020).
4. Casado P, et al. "Kinase-substrate enrichment analysis provides insights into the heterogeneity of signaling pathway activation in leukemia cells." *Science Signaling* (2013).
5. Jurtz V, et al. "NetMHCpan-4.0: Improved Peptide-MHC Class I Interaction Predictions." *Journal of Immunology* (2017).
6. Kim M-S, et al. "A draft map of the human proteome." *Nature* (2014).
7. Wang X, Zhang B. "customProDB: an R package to generate customized protein databases from RNA-Seq data." *Bioinformatics* (2013).
8. Tyanova S, et al. "The Perseus computational platform for comprehensive analysis of (prote)omics data." *Nature Methods* (2016).

---

*Generated: 2026-05-23 | Pipeline version: 1.0.0*
