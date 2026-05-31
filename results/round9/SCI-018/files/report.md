# 実験レポート：抗菌薬耐性（AMR）進化予測計算フレームワーク

**実験日:** 2026-05-31  
**フレームワーク:** AMR-EvoPredict v1.0  
**使用言語:** Python 3.11  

---

## 1. 実験目的と背景

### 目的

抗菌薬耐性（Antimicrobial Resistance: AMR）は現代医療における最も深刻な脅威のひとつであり、WHOは2050年までに年間1,000万人の死者が出ると予測している。本実験では、AMR進化を多スケールで予測する計算フレームワークを構築・実行し、以下の6つの研究課題を解決することを目的とした：

1. 全ゲノムシーケンス（WGS）からの耐性遺伝子（ARG）検出パイプライン
2. 耐性変異の適応度ランドスケープ（fitness landscape）構築
3. 進化経路予測（アクセシブルな変異パス列挙）
4. 水平遺伝子伝達（HGT）のネットワークモデリング
5. 抗菌薬使用量と耐性率の時空間動態モデル
6. 新規抗菌薬投与戦略（組み合わせ療法、サイクリング）の最適化

### 背景

AMRの進化は分子レベル（変異・選択）、ゲノムレベル（ARG・HGT）、集団レベル（疫学動態）の3つのスケールで同時進行する。既存研究は各スケールを独立に解析する傾向があり、統合的な予測フレームワークが不足している。本フレームワークはこのギャップを埋めることを目指した。

---

## 2. 使用した手法・アルゴリズムの概要

### Module 1: ARG検出（ランダムフォレスト）

- **データ:** 500細菌分離株 × 50ゲノム座位（ARG存在/非存在の二値行列）
- **手法:** ランダムフォレスト（n_estimators=100, random_state=42）
- **評価:** 5分割層化交差検証、AUROC
- **対象抗菌薬:** 8種類（Ampicillin, Ciprofloxacin, Tetracycline, Gentamicin, Cefotaxime, Meropenem, Azithromycin, Trimethoprim）

### Module 2: 適応度ランドスケープ

- **モデル:** 4サイト組み合わせ論的ランドスケープ（TEM-1 β-ラクタマーゼにインスパイア）
- **遺伝型数:** 2⁴ = 16
- **適応度関数:** 加法項 + 対形成エピスタシス + ガウスノイズ

$$f(\mathbf{g}) = 1 + \sum_i \alpha_i g_i + \frac{1}{2} \sum_{i \neq j} \beta_{ij} g_i g_j + \varepsilon, \quad \varepsilon \sim \mathcal{N}(0, 0.02^2)$$

### Module 3: 進化経路予測

- **対象:** 野生型（0000）→完全耐性（1111）の全4! = 24経路
- **アクセシブル判定:** 各ステップで適応度が単調非減少
- **経路確率:** Sella-Hirshモデルベースの経路確率

### Module 4: HGTネットワーク

- **ノード:** 30株（5種由来：E.coli, K.pneumoniae, P.aeruginosa, S.aureus, E.faecalis）
- **エッジ:** プラスミド伝達イベント（同種: p=0.18、異種: p=0.05）
- **コミュニティ検出:** 貪欲モジュラリティ最適化

### Module 5: 時空間動態（拡張SIRモデル）

$$\frac{dS}{dt} = -\lambda_S S - \lambda_R S, \quad \frac{dI_S}{dt} = \lambda_S S - \gamma I_S - \mu \phi I_S$$
$$\frac{dI_R}{dt} = \lambda_R S + \mu \phi I_S - \gamma I_R, \quad \frac{dR}{dt} = \gamma(I_S + I_R)$$

- 集団規模: N=10,000、シミュレーション期間: 365日
- パラメータ: β_S=0.35, β_R=0.30, γ=0.10, μ=0.40, φ=0.02, κ=0.15

### Module 6: 治療戦略最適化

3戦略比較：単剤療法（monotherapy）、サイクリング（cycling）、組み合わせ療法（combination）
サイクリング周期最適化: T ∈ {7, 14, 21, 30, 45, 60, 90} 日

### Module 7: 集団遺伝学（Wright-Fisher シミュレーション）

- 有効集団サイズ N_e=1000、200世代、50レプリケート
- 選択係数 s=0.05（抗生物質存在下）、フィットネスコスト 0.02（非存在下）
- 抗生物質存在確率 p_antibiotic=0.6

---

## 3. 先行研究調査結果（ToolUniverse MCP）

Semantic Scholar および PubMed を用いて以下の論文を特定した：

| # | タイトル | 著者 | 年 | DOI | 主要知見 |
|---|---------|------|-----|-----|---------|
| 1 | A convolutional neural network highlights mutations relevant to antimicrobial resistance in *M. tuberculosis* | Green et al. | 2022 | 10.1038/s41467-022-31236-0 | CNN による13薬剤AUC 82.6–99.5%、18の新規耐性関連部位を同定 |
| 2 | Drug resistance prediction in MTB based on hierarchical attentive neural network | Jiang et al. | 2022 | 10.1093/bib/bbac041 | HANN による全ゲノム変異を活用した耐性予測、AUC 97.90%（INH） |
| 3 | Accurate and rapid prediction of tuberculosis drug resistance using ML and CNN | Kuang et al. | 2022 | 10.1038/s41598-022-06449-4 | 10,575株の16ヶ国データ、1D-CNN が最高性能（F1 81.1–98.2%） |
| 4 | Environmental modulation of global epistasis in a drug resistance fitness landscape | Díaz-Colunga et al. | 2023 | 10.1038/s41467-023-43806-x | 薬物濃度がエピスタシスパターンを変調、進化軌跡の予測が困難に |
| 5 | Higher-order epistasis drives evolutionary unpredictability toward novel antibiotic resistance | Gaszek et al. | 2025 | 10.1101/2025.07.08.663783 | 55,296 TEM-1 変異体、8百万測定値、新規基質でより複雑なランドスケープ |
| 6 | A tale of two plasmids: contributions of plasmid associated phenotypes | De Silva et al. | 2022 | 10.1098/rspb.2022.0581 | プラスミド表現型が疫学的成功を決定する |
| 7 | Bacterial plasmid-associated and chromosomal proteins | Downing & Rahm | 2022 | 10.1038/s41598-022-20809-0 | プラスミド由来タンパク質は染色体由来より多くのPPI相互作用を持つ |
| 8 | Interspecies interaction alters trajectory of AMR evolution | Muzafar et al. | 2026 | 10.1093/ismejo/wrag014 | 種間競争が負のエピスタシスを増幅し、より高い耐性変異体を選択 |

**先行研究の課題・限界:**
- ML手法は主に*M. tuberculosis*に限定；他菌種への汎化性が未検証
- 適応度ランドスケープの実験的測定は少数サイトに限定
- HGT動態の定量的モデリングが不足
- WGS解析と疫学モデルの統合フレームワークが欠如

---

## 4. NatureLM / GALACTICA MCP 試行結果

| ツール | 試行ツール名 | エラー内容 | 対処法 |
|--------|------------|-----------|--------|
| NatureLM | `ask_naturelm` | ToolUniverse レジストリに登録なし（0マッチ） | 定量パラメータを文献値から直接取得 |
| GALACTICA | `scientific_qa`, `predict_citations` | ToolUniverse レジストリに登録なし（0マッチ） | 科学的検証を文献クロスリファレンスで代替 |

両MCPへの接続は現環境では不可能であった。実験の科学的透明性のため、この試行・失敗の記録を残す。代替手段として、定量パラメータ（結合自由エネルギー：-6.5 kcal/mol @β-ラクタム/β-ラクタマーゼ複合体; 耐性出現率：φ≈0.01-0.05/antibiotic-treated case）は文献値を使用した。

---

## 5. 主要な結果と数値

### 5.1 ARG検出パフォーマンス

| 抗菌薬 | AUROC平均 | 標準偏差（5-fold CV） |
|--------|-----------|----------------------|
| Ampicillin | 0.876 | ±0.028 |
| Ciprofloxacin | 0.936 | ±0.038 |
| Tetracycline | 0.861 | ±0.060 |
| Gentamicin | 0.950 | ±0.018 |
| Cefotaxime | 0.925 | ±0.028 |
| Meropenem | 0.950 | ±0.033 |
| Azithromycin | 0.932 | ±0.026 |
| Trimethoprim | 0.943 | ±0.033 |
| **平均** | **0.922** | **±0.032** |

![Figure 1: ARG Detection](figures/fig1_arg_detection.png)

### 5.2 適応度ランドスケープ

- 野生型適応度: 1.010
- 最大適応度（完全耐性1111）: 1.679
- 適応度ゲイン（WT→完全耐性）: **+0.669**
- アクセシブルパス: 24/24（**100%**）

![Figure 2: Fitness Landscape](figures/fig2_fitness_landscape.png)

### 5.3 進化経路予測

- **最確進化経路:** 0000 → 0010 → 0011 → 0111 → 1111（確率 0.116）
- 経路確率範囲: 0.020 – 0.116
- 全24経路がアクセシブル（フラットランドスケープ）

### 5.4 HGTネットワーク

- ノード数: 30株、エッジ数: 31転移イベント
- ネットワーク密度: 0.0356
- コミュニティ数: 14
- 最大ドナー株: 株27（外向き次数=6）
- 種内転移: 35.5%、種間転移: 64.5%

![Figure 3: HGT Network](figures/fig3_hgt_network.png)

### 5.5 時空間AMR動態

- 感受性菌ピーク感染者数: 3,300名（12日目）
- 耐性菌ピーク感染者数: 467名（16日目）
- R_eff（感受性株）：ピーク時 1.10
- μ=0.4→0.8時の耐性率上昇: ピーク比で約2.6倍

![Figure 4: AMR Dynamics and Optimization](figures/fig4_dynamics_optimization.png)

### 5.6 治療戦略最適化

| 戦略 | 365日後耐性割合 | 単剤療法比 |
|------|--------------|-----------|
| 単剤療法 | 0.1236 | — |
| サイクリング（30日周期） | 0.1102 | -10.8% |
| 組み合わせ療法 | 0.0813 | **-34.3%** |

**最適サイクリング周期: 45日** (最終耐性割合 = 0.1093)

| サイクリング周期 | 最終耐性割合 | 平均耐性割合 |
|---------------|------------|------------|
| 7日 | 0.1214 | 0.6011 |
| 14日 | 0.1176 | 0.5992 |
| 21日 | 0.1126 | 0.5959 |
| 30日 | 0.1102 | 0.5905 |
| **45日** | **0.1093** | 0.5808 |
| 60日 | 0.1101 | 0.5707 |
| 90日 | 0.1124 | 0.5499 |

### 5.7 集団遺伝学（Wright-Fisher）

- 固定確率: 0.12（50レプリケート中6回固定）
- 200世代後平均頻度: 0.118 ± 0.196
- t検定（vs 中立期待値0.5）: t = -13.669, **p < 0.0001**

![Figure 5: Population Genetics](figures/fig5_population_genetics.png)

---

## 6. 考察と今後の展望

### 自己批判的評価

**シミュレーションデータの前提依存:**
本研究の全結果は合成データに基づく。ARG-表現型関係が完全に既知という仮定の下でのAUROC=0.922は楽観的な推定値であり、実臨床データでは5–15%の低下が予想される。実世界ではラボ表現型判定に5–15%の誤差が存在し、ゲノムの特徴量も数千変数に及ぶ。

**適応度ランドスケープの単純化:**
100%のアクセシブルパスは、対形成エピスタシスのみを含むスムーズなモデルの産物であり、現実の多次エピスタシス（Gaszek et al., 2025）は景観をより粗くする。実際の耐性進化では、高次エピスタシスにより景観がより複雑になり、アクセシブルパスが大幅に減少する可能性がある。

**SIRモデルの限界:**
- 空間的不均一性の無視（地理的AMR拡散を過小評価）
- 宿主免疫多様性の無視
- 複数薬剤間の交差耐性の未考慮

**治療最適化の限界:**
組み合わせ療法の優位性（-34.3%）は薬物拮抗作用や耐性の二重獲得を考慮していない。実際の臨床設定では、相乗的・拮抗的薬物相互作用により結果が大きく変動する。

### 今後の展望

1. **実臨床データによる検証:** PATRICデータベースやResFinder等を用いた外部検証
2. **空間モデル拡張:** パッチモデル（地域間移動を考慮したSIRネットワーク）
3. **多薬剤耐性動態:** 複数ARGの共存・共進化をモデル化
4. **深層学習統合:** Transformer系モデルによるARG検出精度向上
5. **実験的ランドスケープとの統合:** TEM-1実測データ（Gaszek et al.）との照合

---

## 7. 計算来歴（Computational Provenance）

| 項目 | 詳細 |
|------|------|
| 乱数シード | `np.random.seed(42)`, `random.seed(42)` |
| Python バージョン | 3.11 |
| コードファイル | `amr_analysis.py` |
| データ出自 | 合成データ（シミュレーション、パラメータは文献値に基づく） |
| 保存先 | `data/raw/*.csv` |

---

## 8. 生成したファイル一覧

| ファイル | 種別 | 説明 |
|--------|------|------|
| `amr_analysis.py` | Python スクリプト | 全モジュールの実装コード |
| `figures/fig1_arg_detection.png` | 図 | ARG検出結果（ヒートマップ + AUROCバーチャート） |
| `figures/fig2_fitness_landscape.png` | 図 | 適応度ランドスケープ（バイオリン + 2D + パス散布図） |
| `figures/fig3_hgt_network.png` | 図 | HGTネットワーク（グラフ + 次数分布） |
| `figures/fig4_dynamics_optimization.png` | 図 | 時空間動態 + 治療戦略最適化 |
| `figures/fig5_population_genetics.png` | 図 | 集団遺伝学シミュレーション |
| `data/raw/genome_matrix.csv` | データ | 500株 × 50座位の遺伝子行列 |
| `data/raw/resistance_phenotypes.csv` | データ | 8抗菌薬の耐性表現型ラベル |
| `data/raw/auroc_results.csv` | データ | 5-fold CV AUROCスコア |
| `data/raw/fitness_landscape.csv` | データ | 16遺伝型の適応度値 |
| `data/raw/hgt_events.csv` | データ | HGT伝達イベントリスト |
| `data/raw/sir_dynamics.csv` | データ | SIRモデル時系列データ |
| `data/raw/cycling_optimization.csv` | データ | サイクリング周期最適化結果 |
| `paper.md` | 文書 | 学術論文形式レポート（英語） |
| `report.md` | 文書 | 実験レポート（日本語） |

---

## 参考文献

1. Green et al. (2022). Nature Communications. DOI: 10.1038/s41467-022-31236-0
2. Jiang et al. (2022). Briefings in Bioinformatics. DOI: 10.1093/bib/bbac041
3. Kuang et al. (2022). Scientific Reports. DOI: 10.1038/s41598-022-06449-4
4. Díaz-Colunga et al. (2023). Nature Communications. DOI: 10.1038/s41467-023-43806-x
5. Gaszek et al. (2025). bioRxiv. DOI: 10.1101/2025.07.08.663783
6. De Silva et al. (2022). Proc. R. Soc. B. DOI: 10.1098/rspb.2022.0581
7. Downing & Rahm (2022). Scientific Reports. DOI: 10.1038/s41598-022-20809-0
8. Muzafar et al. (2026). The ISME Journal. DOI: 10.1093/ismejo/wrag014
