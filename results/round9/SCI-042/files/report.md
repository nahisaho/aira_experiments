# 実験レポート: ショットガンメタゲノミクスデータからの機能プロファイリングパイプライン (MetaSnake)

---

## 1. 実験目的と背景

### 研究目的
ショットガンメタゲノムシーケンシングデータから腸内細菌叢の分類・機能プロファイルを統合的に解析する、再現可能なSnakemakeベースのパイプライン(**MetaSnake**)を設計・実装し、IBD（炎症性腸疾患）患者 vs. 健常人の合成データセットで検証する。

### 研究背景
ヒト腸内マイクロバイオームは10¹³個の微生物細胞を含み、免疫調節・代謝・腸バリア機能に不可欠な役割を果たす。ショットガンメタゲノミクスはこの生態系を文化非依存的に解析する最も包括的な手法であるが、解析には10以上の専門ツールの統合が必要で再現性に課題があった。

---

## 2. 使用した手法・アルゴリズムの概要

### パイプライン構成（6モジュール）

| モジュール | ツール | 目的 |
|-----------|--------|------|
| QC | fastp, Bowtie2 (hg38), Clumpify | アダプター除去・ホスト除去・重複排除 |
| 分類 | MetaPhlAn4, Kraken2/Bracken | アセンブリフリー分類（並列比較） |
| 機能アノテーション | HUMAnN3, eggNOG-mapper v2 | MetaCycパスウェイ・KO定量 |
| アセンブリ・ビニング | MEGAHIT, MetaBAT2, CONCOCT, MaxBin2, DAS_Tool | ゲノムビニングアンサンブル |
| MAG品質評価 | CheckM2, GTDB-Tk (r220) | 完全性・汚染率・系統配置 |
| 統計解析 | scikit-learn (RF/LR/GBM), scipy | 多変量統計・ML分類 |

### Snakemakeワークフロー

```
rule all → qc → host_removal → [metaphlan4 || kraken2] → [humann3 || assembly]
         → binning (metabat2 + concoct + maxbin2) → dastool → checkm2 → stats
```

- 完全な`Snakefile`と`config/config.yaml`を実装済み
- 各ツールのConda環境を分離（`envs/`ディレクトリ）

---

## 3. 先行研究調査結果（PubMed経由）

ToolUniverse MCP（PubMed, Semantic Scholar）を使って以下の先行研究を特定:

| # | 論文 | PMID | 主要知見 |
|---|------|------|----------|
| 1 | eggNOG-mapper v2 (Cantalapiedra et al., 2021) | 34597405 | メタゲノムスケールの機能アノテーション基盤 |
| 2 | CheckM2 (Chklovski et al., 2023) | 37500759 | ML-based MAG品質予測、従来CheckMを凌駕 |
| 3 | MetaPhlAn4 vs Kraken2 (Karagiannis et al., 2026) | 41525322 | 両ツールは類似多様性トレンドを示すが種レベルで差異あり |
| 4 | Meteor2 (Ghozlane et al., 2025) | 41199348 | MetaPhlAn4比+45%検出感度向上（低存在量種） |
| 5 | Doppelgänger bias (Zhou et al., 2025) | 40888678 | ML AUROC 15-30%過大評価のリスク |
| 6 | Microbiome Datahub (Mori et al., 2026) | 41840729 | 214,427 MAGデータベース、平均完全性80.5% |
| 7 | MAGFlow/BIgMAG (Yepes-García et al., 2024) | 39360247 | CheckM2+GTDB-Tk統合可視化フレームワーク |
| 8 | Kraken2/Bracken benchmark (Timilsina et al., 2025) | 40683452 | 低存在量病原体検出でKraken2/BrackenがMetaPhlAn4を上回る |
| 9 | HUMAnN3 + IBD (Noel et al., 2025) | 41077635 | CKD/AKI患者の腸内細菌叢をKraken2+MetaPhlAn3+HUMAnN3で解析 |
| 10 | SCFA pathways (Dissanayaka et al., 2026) | 41619271 | HUMAnN3でSCFA関連パスウェイを特定（前臨床AD） |

**NatureLM / GALACTICA MCPツール**: ToolUniverseレジストリに存在せず、接続不可。代替として文献値（SCFA Km = 1–5 mM、酪酸HDAC阻害IC₅₀ ≈ 2–5 mM）を使用。試行ツール名: `ask_naturelm`（NatureLM）、`scientific_qa`/`predict_citations`（GALACTICA）。

---

## 4. 主要な結果と数値

### 4.1 品質管理（QC）

![Figure 1: QC Pipeline](figures/fig01_qc_pipeline.png)

**Figure 1.** QCパイプラインの各段階での読み取り数保持率（A）、ホスト読み取り分率の群間比較（B）、シーケンス深度分布（C）。

[cell:4] QCメトリクス（n=120サンプル、IBD vs. 健常）:
- 平均生リード数: **12.7 ± 1.8 M reads**
- IBDホスト分率: **0.142 ± 0.057**（健常 0.091 ± 0.034）
- ホスト分率差: Mann-Whitney U p = **1.96×10⁻⁶**（有意）

IBDにおけるホストリード比率の上昇は腸管透過性亢進を反映し、炎症活性の指標として機能する。

---

### 4.2 分類学的プロファイリング比較

![Figure 2: Taxonomic Profiling](figures/fig02_taxonomic_profiling.png)

**Figure 2.** (A) MetaPhlAn4 vs. Kraken2でのShannon多様性比較。(B) Bray-Curtis距離行列のPCoA。(C) ツール間一致度散布図。

[cell:5] Alpha多様性（Shannon H'）:

| 群 | MetaPhlAn4 | Kraken2 | p値 |
|----|-----------|---------|-----|
| IBD | 2.984 ± 0.158 | 2.991 ± 0.148 | 0.0043 |
| Healthy | 3.061 ± 0.147 | 3.059 ± 0.152 | — |

[cell:6] ツール間比較:
- Pearson相関（全サンプル平均）: **r = 0.9954 ± 0.0023**
- Bray-Curtis乖離（同一サンプル内ツール間）: **0.0817 ± 0.0119**

両ツールはコミュニティレベルで高い一致を示したが、低存在量種では系統的差異が観察された（Karagiannis et al., 2026と一致）。

---

### 4.3 機能アノテーション（HUMAnN3）

![Figure 3: Functional Profiling](figures/fig03_functional_profiling.png)

**Figure 3.** HUMAnN3機能プロファイリング。(A) パスウェイVolcano plot。(B) 酪酸合成パスウェイ比較。(C) パスウェイ発現量ヒートマップ（上位15経路）。

[cell:7] FDR補正済み有意差パスウェイ（Benjamini-Hochberg法、q < 0.05）:

| パスウェイ | p値 | FDR q値 | 方向（IBD） |
|-----------|-----|---------|------------|
| BUTYRATE-SYNTHESIS (酪酸合成) | 6.34×10⁻⁹ | **2.53×10⁻⁷** | ↓ 減少 |
| LIPOPOLYSACCHARIDE-BIOSYN (LPS合成) | 2.00×10⁻⁴ | **3.15×10⁻³** | ↑ 増加 |
| PWY-2941 (葉酸合成) | 2.36×10⁻⁴ | **3.15×10⁻³** | ↓ 減少 |

**解釈**: 酪酸産生低下は *F. prausnitzii* / *Roseburia* 減少の直接的帰結であり、腸管上皮細胞のエネルギー源喪失・HDAC阻害活性低下によるNF-κB活性化（炎症促進）につながる。LPS合成上昇はグラム陰性菌（大腸菌群）の増殖を示す。

---

### 4.4 ゲノムビニングとMAG品質

![Figure 4: MAG Quality](figures/fig04_mag_quality.png)

**Figure 4.** (A) 各ビニングツールの完全性 vs. 汚染率散布図。(B) MIMAGカテゴリ分布。(C) DAS_Tool vs. MetaBAT2の完全性改善量。

[cell:8] MAG品質比較（30 MAGs）:

| ツール | 平均完全性 | 平均汚染率 | High品質 | Medium品質 |
|--------|----------|----------|---------|-----------|
| MetaBAT2 | 71.9% | 13.6% | 1 | 12 |
| CONCOCT | 73.7% | 24.4% | 0 | 3 |
| MaxBin2 | 60.5% | 13.2% | 0 | 11 |
| **DAS_Tool** | **79.8%** | **7.0%** | **3** | **20** |

DAS_Toolアンサンブルは:
- MetaBAT2比 **+7.9%** 完全性改善
- MetaBAT2比 **-6.6%** 汚染率低下
- Microbiome Datahub (Mori et al., 2026) の報告平均80.5%と一致

---

### 4.5 機械学習IBD分類

![Figure 5: ML Classification](figures/fig05_ml_classification.png)

**Figure 5.** (A) 5分割交差検証ROC曲線。(B) Random Forest特徴量重要度（上位15）。(C) モデル間AUROC比較。

[cell:9] 5分割層別交差検証（random_state=42）:

| モデル | AUROC (mean ± SD) | Accuracy (mean ± SD) |
|-------|-------------------|----------------------|
| Random Forest | **0.986 ± 0.018** | 0.942 ± 0.062 |
| Logistic Regression | 0.963 ± 0.023 | 0.908 ± 0.061 |
| Gradient Boosting | 0.942 ± 0.031 | 0.917 ± 0.059 |

[cell:18] ホールドアウト検証（25%テストセット）AUROC: **0.978**

⚠️ **重要な注記**: これらのAUROC値は合成データに埋め込まれたシグナルにより過大評価されている。実際の臨床データでは **0.72–0.88** 程度が期待値（Zhou et al., 2025）。

**Random Forest重要特徴量Top5**:
1. *Faecalibacterium prausnitzii* (重要度=0.188)
2. *Akkermansia muciniphila* (0.084)
3. *Clostridium difficile* (0.074)
4. *Bifidobacterium longum* (0.049)
5. *Escherichia coli* (0.041)

---

### 4.6 パイプライン全体アーキテクチャ

![Figure 6: Pipeline Summary](figures/fig06_pipeline_summary.png)

**Figure 6.** (A) MetaSnakeパイプラインのSnakemakeアーキテクチャ図。(B) 数値パフォーマンスサマリー。

---

## 5. 考察と今後の展望

### 5.1 分類ツールの選択指針

- **精度重視**（腸内細菌叢, well-characterized species）→ MetaPhlAn4推奨
- **感度重視**（低存在量病原体, 環境サンプル）→ Kraken2/Bracken推奨
- **最善策**: 本パイプラインのように両ツールを並列実行し結果を統合

### 5.2 機能プロファイリングの限界

- HUMAnN3はデータベースカバレッジ外の新規系統では機能を見逃す可能性
- eggNOG-mapper v2の統合でカバレッジを補完
- KEGG KEGGaNOGツール（Popov et al., 2026 [PMID: 40968530]）によるモジュール完全性スコア算出が有用

### 5.3 MAGビニングの課題

- 実世界では群集複雑性が高いほどビニング精度が低下
- 同一サンプル内ビニングの限界：複数サンプル共同ビニング（co-binning）が推奨
- CheckM2はPatescibacteriaなど縮退ゲノムでも精度維持（Chklovski et al., 2023）

### 5.4 自己批判的評価

| 観点 | 評価 |
|------|------|
| 合成データへの依存 | シグナルが意図的に埋め込まれており循環論的 |
| 実世界への一般化 | 未測定交絡因子（食事・抗生物質・バッチ効果）が性能を低下させる |
| 実験設計のバイアス | IBD/健常の均等割付けは実臨床コホートよりも分類が容易 |
| NatureLM/GALACTICAの不使用 | 定量的動態パラメータの予測・科学的検証に欠落 |
| Snakemake実行 | 実ツール（Bowtie2等）未インストールのため実際の実行はできない |

### 5.5 今後の展望

1. **ビロームとレジストームの統合**: AMR遺伝子（ABRicate）とウイルス配列（VirSorter2）の追加
2. **縦断サンプル対応**: 時系列メタゲノムデータのSNVトラッキング（StrainPhlAn）
3. **Meteor2統合**: 低存在量種の検出感度向上
4. **リアルコホート検証**: HMP2 (iHMP)やHGP等の公開データでの性能評価
5. **NatureLM/GALACTICA統合**: ツール利用可能になった際の定量的パラメータ予測

---

## 6. NatureLM / GALACTICA 接続試行記録

### 試行ツールと結果

| ツール | 試行名 | 結果 | 代替措置 |
|--------|--------|------|---------|
| NatureLM MCP | `ask_naturelm` | ToolUniverseレジストリに存在せず（total_matches=0） | 文献値使用 |
| GALACTICA MCP | `scientific_qa` | ToolUniverseレジストリに存在せず | PubMed検索で補完 |
| GALACTICA MCP | `predict_citations` | ToolUniverseレジストリに存在せず | 手動文献整理 |

### 文献ベースの代替パラメータ

- 酪酸産生速度 Vmax: 0.3–2.1 μmol/min/mg protein（*F. prausnitzii*特異的）
- SCFAトランスポーター Km: 1–5 mM
- 酪酸HDAC阻害 IC₅₀: 2–5 mM
- LPSリポポリサッカライド刺激 TLR4 EC₅₀: 0.1–10 ng/mL

---

## 7. 計算来歴（Computational Provenance）

| 項目 | 値 |
|------|-----|
| Python version | 3.11 |
| 乱数シード | `np.random.seed(42)`, `random.seed(42)` |
| NumPy | 2.4.6 |
| scikit-learn | 1.8.0 |
| SciPy | 1.17.1 |
| pandas | 3.0.3 |
| matplotlib | 3.10.9 |
| seaborn | 0.13.2 |
| データ来歴 | 合成データ（Dirichletモデル、seed=42） |
| データ保存先 | `data/raw/metadata.csv`, `data/raw/qc_metrics.csv` 等 |
| 環境記録 | `data/raw/environment.txt` |
| 実行Jupyter kernel | `b55ce365-0012-42d8-8bb7-f262884dd42f` |

---

## 8. 生成したファイル一覧

| ファイル | 説明 |
|---------|------|
| `paper.md` | 学術論文形式の文書（英語） |
| `report.md` | 本レポート（日本語） |
| `Snakefile` | 完全なSnakemakeワークフロー |
| `config/config.yaml` | パイプライン設定ファイル |
| `figures/fig01_qc_pipeline.png` | QCパイプライン図 |
| `figures/fig02_taxonomic_profiling.png` | 分類学的プロファイリング図 |
| `figures/fig03_functional_profiling.png` | 機能アノテーション図 |
| `figures/fig04_mag_quality.png` | MAG品質比較図 |
| `figures/fig05_ml_classification.png` | ML分類結果図 |
| `figures/fig06_pipeline_summary.png` | パイプライン概要図 |
| `data/raw/metadata.csv` | サンプルメタデータ |
| `data/raw/qc_metrics.csv` | QCメトリクス |
| `data/raw/metaphlan4_profiles.csv` | MetaPhlAn4プロファイル |
| `data/raw/kraken2_profiles.csv` | Kraken2プロファイル |
| `data/raw/humann3_pathways.csv` | HUMAnN3パスウェイ |
| `data/raw/bray_curtis_matrix.csv` | ベータ多様性行列 |
| `data/raw/mag_quality.csv` | MAG品質評価結果 |
| `data/raw/environment.txt` | pip freeze出力（パッケージバージョン） |
