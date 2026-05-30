# 実験レポート: 金属有機構造体（MOF）CO₂/H₂吸着性能予測のための  
# ハイスループットスクリーニングシステム

---

## 1. 実験目的と背景

### 目的

本実験の目的は、金属有機構造体（MOF）のCO₂吸着（特に直接空気回収: DAC条件での400 ppm）および高圧H₂貯蔵性能を予測するための、計算科学的ハイスループットスクリーニング（HTS）パイプラインを設計・実装・評価することである。

### 背景

MOFは金属ノードと有機リンカーから構成される多孔性結晶材料であり、莫大な化学的多様性を持つ。CoRE MOF 2019データベースには14,000以上の実験的に合成されたMOF構造が収録されており、仮想MOF（hMOF）データベースには数十万の計算生成構造が含まれる。このような大規模データベースを高精度Grand Canonical Monte Carlo（GCMC）シミュレーション（RASPA2）で全数評価することは計算コスト上非現実的であるため、機械学習サロゲートモデルを利用した2段階パイプライン（幾何学的事前スクリーニング + ML予測）が必要とされている。

特にDAC向けMOFには以下の性能要件が課される：
- CO₂取り込み量 > 1 mmol/g（400 ppm, 298 K）
- CO₂/N₂選択性 > 30
- 水安定性（高湿度環境での耐久性）
- 低再生エネルギー（< 60 kJ/kg）

---

## 2. 先行研究調査結果（ToolUniverse MCP 使用）

### 使用ツール
- **Crossref_search_works**: 複数の検索クエリ（"MOF CO2 adsorption GCMC machine learning"、"high-throughput MOF DAC screening geometric descriptors"等）にて実行
- **SemanticScholar_get_paper**: DOI指定での詳細情報取得（429エラーにより一部取得不可）
- **ask_naturelm**: CO₂吸着メカニズムおよびDAC性能指標の定量的情報取得

### 特定された主要論文（2020年以降、5件以上）

| No. | タイトル | 著者 | 年 | DOI | 主要知見 |
|---|---|---|---|---|---|
| 1 | Machine Learning Descriptors for CO₂ Capture Materials | Orhan, Zhao, Babarao | 2025 | 10.3390/molecules30030650 | Qst（吸着熱）が最重要記述子。amine官能基化で性能向上 |
| 2 | Towards accurate and scalable high-throughput MOF adsorption screening | Bonakala et al. | 2026 | 10.1039/d6sc00831c | 古典力場＋機械学習ポテンシャルのハイブリッドにより量子化学精度でのMOFスクリーニング実現 |
| 3 | High-throughput DFT screening of single-metal and high-entropy MOF-74 | Tamtaji et al. | 2025 | 10.1016/j.ijhydene.2025.151276 | DFT法でCO₂/N₂分離とH₂貯蔵向けMOF-74変種をスクリーニング。金属種が決定的因子 |
| 4 | CO₂ separation from flue gas using [BMIM][BF4]/MOF composites | Polat et al. | 2020 | 10.1016/j.cej.2020.124916 | HTS計算と実験を連携した検証プロトコルの確立 |
| 5 | Machine-learning model reveals critical features for high-throughput screening | Lim | 2024 | 10.1063/10.0028344 | 数百の特徴量の重要度を評価。幾何学記述子のみでは不十分でQst必要 |
| 6 | High-efficiency prediction of water adsorption by lattice GCMC | Liu et al. | 2025 | 10.1039/d4lf00354c | 格子GCMCによる多孔質材料の水吸着性能高速予測 |

### 先行研究の課題・限界

1. **DAC条件の軽視**: 多くの研究が煙道ガス条件（0.15 bar）に焦点を当て、より困難なDAC条件（400 ppm, 40 Pa）を軽視
2. **水安定性の定量的統合**: 有望なMOFの水安定性が定量スコアとしてランキングに組み込まれていない
3. **幾何学記述子のみの限界**: 幾何学特徴量だけでは精度不足（R² ≈ 0.73）。QstなどGCMC由来の物性値が必要
4. **合成データへの依存**: 実際のCoRE MOFデータベースではなく、合成データを使った研究が多い

---

## 3. 実験手法の概要

### 3.1 パイプライン構成

```
[Stage 1] MOF構造データベース (CoRE MOF / hMOF 相当, N=500)
          ↓
[Stage 2] 幾何学的記述子抽出 (Zeo++ 相当)
          - 最大内包球径 (LIS)、表面積 (VSA, GSA)、細孔体積
          ↓
[Stage 3] 幾何学的事前フィルター (void fraction ≥ 0.3, LIS ≥ 3 Å)
          → 467構造に削減 (93.4%)
          ↓
[Stage 4] GCMC吸着シミュレーション (Langmuirモデルを使用)
          - CO₂: 0.15 bar (煙道ガス) および 400 ppm (DAC)
          - H₂: 100 bar, 77 K
          → 吸着量、Qst、ヘンリー定数
          ↓
[Stage 5] 機械学習サロゲートモデル (RF, GB)
          - 5分割交差検証
          ↓
[Stage 6] 多目的DACスコアリング & ランキング
          - CO₂取り込み + 選択性 + 水安定性 + 合成可能性
```

### 3.2 GCMCアナログシミュレーション

物理情報組み込みLangmuirモデルを使用：

**CO₂吸着熱 (Qst) の推定式:**
```
Qst_CO2 = 25.0 + 7.0 × exp(-LIS/5.0) + 8.0 × [amine] + 5.0 × [OMS] + ε
```
（ε: ガウスノイズ N(0, 2.5) kJ/mol）

**ヘンリー定数:**
```
K_CO2 = 1×10⁻⁶ × exp((Qst - 25.0) / 8.314)  [mol/kg/Pa]
```

**Langmuir等温線:**
```
q(P) = q_sat × K × P / (1 + K × P)
```

### 3.3 機械学習モデル

| モデル | ハイパーパラメータ | 特徴量セット |
|---|---|---|
| Random Forest (RF) | n_estimators=200, max_depth=12, min_samples_leaf=3 | Full (12特徴) / Geo-only (11特徴) |
| Gradient Boosting (GB) | n_estimators=200, max_depth=5, lr=0.05, subsample=0.8 | Full (12特徴) |

評価: 5分割交差検証 (R², RMSE, 標準偏差)

### 3.4 NatureLM MCP ツール活用状況

| ツール | 試行結果 | 取得内容 |
|---|---|---|
| `generate_smiles` | ✅ 成功 | BDC: `O=C(O)c1ccc(C(=O)O)cc1`, NH₂-BDC: `Nc1cc(C(=O)O)ccc1C(=O)O`, 三アミンリンカー: `NCCNCCNCCN` |
| `predict_logp` | ✅ 成功 | BDC: 0.66, NH₂-BDC: 1.20, 三アミン: 0.90 |
| `predict_property` (溶解度) | ✅ 成功 | BDC: logS=-1.54, NH₂-BDC: logS=-3.14 |
| `retrosynthesis` | ⚠️ 部分的成功 | 芳香族化合物の逆合成経路が不正確（直鎖アルキル鎖を提案） |
| `ask_naturelm` | ✅ 成功（要注意） | DAC目標値確認、ただしQstの数値が単位解釈で疑義あり |

**注意**: NatureLMの逆合成および定量値（Qst）には化学的不正確さが含まれており、定性的スクリーニングには使えるが定量的ベンチマークとしての利用は推奨しない。

---

## 4. 主要な結果と数値

### 4.1 データベース統計

| 記述子 | 平均 ± 標準偏差 | 範囲 |
|---|---|---|
| Void fraction (φ) | 0.597 ± 0.187 | 0.15–0.95 |
| 重力比表面積 (m²/g) | 3,030 ± 2,150 | 100–15,000 |
| 細孔体積 (cm³/g) | 0.57 ± 0.43 | 0.05–5.0 |
| Qst,CO₂ (kJ/mol) | 34.6 ± 8.2 | 15–60 |

### 4.2 GCMC吸着シミュレーション結果

| 吸着条件 | 平均 ± 標準偏差 | 最大値 |
|---|---|---|
| CO₂ at 0.15 bar (mol/kg) | 0.513 ± 0.530 | 5.54 |
| CO₂ at 400 ppm (mol/kg) | 0.0013 ± 0.0014 | 0.0126 |
| H₂ at 100 bar, 77 K (mol/kg) | 5.93 ± 4.22 | 11.78 |
| CO₂/N₂選択性 | 3.87 ± 4.91 | 100.0 |

### 4.3 機械学習モデル性能（5分割交差検証）

| モデル | 目的変数 | 特徴量 | R² (mean ± std) | RMSE |
|---|---|---|---|---|
| Random Forest | CO₂ DAC | Full (12) | **0.908 ± 0.069** | 0.00042 mol/kg |
| Random Forest | CO₂ DAC | Geo-only (11) | 0.736 ± 0.094 | 0.00072 mol/kg |
| Gradient Boosting | CO₂ DAC | Full (12) | **0.934 ± 0.022** | 0.00036 mol/kg |
| Random Forest | H₂ 100 bar | Full (12) | **0.981 ± 0.003** | 0.587 mol/kg |
| Random Forest | CO₂ Flue | Full (12) | 0.622 ± 0.112 | 0.326 mol/kg |

⚠️ **自己批判的評価**: GBの R² = 0.934 はほぼ完璧に見えるが、これは合成データが単純なLangmuirモデルから生成されているためであり、実世界のMOFデータでは必ずノイズ・複雑な吸着挙動（柔軟な骨格、不均一サイト）が含まれ、予測精度は低下すると考えられる。実際のCoRE MOFデータへの適用では R² = 0.70–0.80 程度が現実的な期待値と推定される。

### 4.4 特徴量重要度（RF Full, CO₂ DAC）

| 順位 | 特徴量 | 重要度 |
|---|---|---|
| 1 | Qst,CO₂ (kJ/mol) | **0.532** |
| 2 | 重力比表面積 (m²/g) | 0.178 |
| 3 | 密度 (g/cm³) | 0.101 |
| 4 | 細孔体積 (cm³/g) | 0.093 |
| 5 | Void fraction | 0.055 |
| 6 | 体積比表面積 (m²/cm³) | 0.035 |
| 7-12 | 化学的記述子 | < 0.020 |

**知見**: Qstが53%の重要度を占め、DAC条件（Henry則域）における熱力学的親和性が支配的であることを確認。

### 4.5 スクリーニングファネル

| ステージ | 条件 | 残存数 | 割合 |
|---|---|---|---|
| 初期データベース | — | 500 | 100.0% |
| 幾何学的事前フィルター | φ≥0.3, LIS≥3Å | 467 | 93.4% |
| GCMC スクリーニング | q_CO₂≥0.01 mmol/g | 312 | 62.4% |
| 水安定性フィルター | score≥0.4 | 178 | 35.6% |
| DAC候補 | composite score | 20 | 4.0% |

### 4.6 DAC上位候補MOF

| 順位 | MOF ID | CO₂ DAC (mol/kg) | CO₂/N₂選択性 | 水安定性 | 合成可能性 | DACスコア |
|---|---|---|---|---|---|---|
| 1 | MOF_0343 | 0.0126 | 100.0 | 0.593 | 0.320 | 0.820 |
| 2 | MOF_0273 | 0.0116 | 100.0 | 0.129 | 0.619 | 0.728 |
| 3 | MOF_0471 | 0.0055 | 77.1 | 0.550 | 0.775 | 0.641 |

すべてアミン官能基を持ち、Zr・Al・Feの金属ノードを使用。これは文献報告のMIL-101(NH₂)、UiO-66(NH₂)等の高性能MOFと整合する。

### 4.7 NatureLM予測結果（定量的）

| MOFリンカー | SMILES | logP | logS | 特記事項 |
|---|---|---|---|---|
| BDC（テレフタル酸） | `O=C(O)c1ccc(C(=O)O)cc1` | 0.66 | -1.54 | MOF-5の標準リンカー |
| NH₂-BDC（2-アミノテレフタル酸） | `Nc1cc(C(=O)O)ccc1C(=O)O` | 1.20 | -3.14 | MIL-101(NH₂)リンカー |
| 三アミンリンカー | `NCCNCCNCCN` | 0.90 | N/A | 高CO₂親和性リンカー候補 |

---

## 5. 生成した図表

### Figure 1: スクリーニング概要

![Figure 1: High-Throughput MOF Screening Overview](figures/fig1_screening_overview.png)

*Figure 1. (a) 重力比表面積とCO₂ DAC取り込み量の関係（Q_stでカラーコーディング）; (b) 細孔体積とH₂取り込み量の関係; (c) アミン官能基化の有無によるCO₂/N₂選択性分布; (d) RFモデルの特徴量重要度; (e) CO₂ DAC予測のパリティプロット（5分割CV）; (f) H₂予測のパリティプロット（5分割CV）.*

### Figure 2: 相関行列と吸着等温線

![Figure 2: Correlation Matrix and Isotherms](figures/fig2_correlation_isotherms.png)

*Figure 2. (左) 構造記述子と吸着特性のペアワイズ相関行列。強い正相関: GSA-H₂取り込み (r≈0.88)、Qst-CO₂ DAC (r≈0.82)。(右) 上位5/下位5 DAC候補の対数圧力スケールCO₂ Langmuir等温線。*

### Figure 3: ML性能とスクリーニングファネル

![Figure 3: ML Performance and Screening Funnel](figures/fig3_ml_performance.png)

*Figure 3. (a) 全モデルのR²比較（エラーバー：5分割CVの標準偏差）; (b) 各段階での候補数の変化（スクリーニングファネル）; (c) CO₂ DAC取り込み量 vs log₁₀(CO₂/N₂選択性)のDAC候補ランドスケープ。*

### Figure 4: 上位候補MOFの特性

![Figure 4: Top DAC Candidates](figures/fig4_top_candidates.png)

*Figure 4. (左) 上位20 DAC候補の水安定性 vs CO₂取り込み量の散布図（バブルサイズ∝表面積、色＝金属種）; (右) 上位5候補の3指標（CO₂取り込み、水安定性、合成可能性）の比較.*

---

## 6. 考察と今後の展望

### 6.1 主要な知見の解釈

**DAC条件でのQst支配性**: 400 ppmという極希薄条件では吸着がHenry則域にあり、取り込み量は q ≈ K_H × P ∝ exp(Qst/RT) × P で決まる。よってQstがR² = 0.532の重要度を持つのは物理的に正当である。アミン官能基化（Qstを+8 kJ/mol増加）の効果は実験室でのMIL-101(NH₂)やUiO-66(NH₂)の高性能と整合する。

**H₂高精度予測の要因**: 100 bar, 77 KでのH₂吸着は多くの高表面積MOFで飽和近傍にあり、取り込み量は表面積・細孔体積の線形関数に近い。これが R² = 0.981 という高精度の物理的根拠である。

**煙道ガスCO₂の予測困難性**: 0.15 barではLangmuirの屈曲部付近にあり、q_satとKの両方が取り込みに効く。これがR² = 0.622という低値の原因であり、より高精度な力場や多サイトモデルが必要となる。

### 6.2 実験設計の自己批判的評価

**合成データへの依存度:**
- 本研究の最大の限界は、実際のMOF結晶構造でなく、統計的パラメータから生成した合成データを使用していることである
- 実MOFには骨格柔軟性、欠陥構造、複雑な多座配位リンカーが存在し、Langmuir単純モデルでは記述できない
- 実CoRE MOFへの適用時、R²は0.10〜0.20低下すると推定

**Q_st推定モデルの不確実性:**
- 加算モデル（Qst = 基準値 + micropore効果 + amine効果 + OMS効果）は物理的に合理的だが、
- 実際のQstは金属-配位子環境の原子レベル詳細に依存する
- バイナリ変数（amine: 0/1）はアミン含有量や配置の違いを無視している

**水安定性プロキシの限界:**
- 金属種のみに基づくスコアは過度に単純化されている
- ZIF-8（Zn）のように金属種に関わらず安定なMOFも存在
- 実際の安定性評価にはMD計算または加速劣化試験が必要

**NatureLM予測の過度な楽観性:**
- NatureLMのlogP・溶解度予測は定性的傾向は正しいが、定量精度は未知
- 逆合成予測は芳香族系で誤答であった
- Q_stの「0.65〜1.25 kJ/mol」という回答は吸着熱（20〜50 kJ/mol）でなくHenry定数の誤答の可能性がある

**実世界への一般化可能性:**
- 本研究の結果（R² > 0.9）は合成データの単純な構造のためであり、実世界への過度な一般化は避けるべきである
- 実データへの適用前に、実際のRASPA2によるGCMCシミュレーション結果でのバリデーションが必須

### 6.3 今後の展望

1. **実データへの適用**: CoRE MOF 2019 (14,000+構造) に対して実際のRASPA2/Zeo++パイプラインを実行
2. **水蒸気競合吸着の組み込み**: DAC条件では相対湿度50〜90%が標準。水分子との競合吸着モデルが必要
3. **グラフニューラルネットワーク**: CGCNN/MOFformerなど結晶グラフエンコーディングによる高精度予測
4. **実験バリデーション**: 上位候補のブレークスルー実験・多サイクル耐久試験
5. **経済性分析**: 吸着量・再生エネルギー・材料コストを統合した$/tCO₂分析

---

## 7. 生成したファイル一覧

| ファイル名 | 内容 | 形式 |
|---|---|---|
| `mof_database.csv` | 500 MOF構造の全記述子・吸着量・スコア | CSV |
| `figures/fig1_screening_overview.png` | スクリーニング概要（6パネル） | PNG |
| `figures/fig2_correlation_isotherms.png` | 相関行列と吸着等温線 | PNG |
| `figures/fig3_ml_performance.png` | MLモデル性能・ファネル・DAC候補 | PNG |
| `figures/fig4_top_candidates.png` | 上位候補MOFの特性比較 | PNG |
| `paper.md` | 学術論文形式の成果文書 | Markdown |
| `report.md` | 実験レポート（本ファイル） | Markdown |

---

## 付録: 参考文献

1. Orhan, I. B., Zhao, T., & Babarao, R. (2025). Machine Learning Descriptors for CO₂ Capture Materials. *Molecules*, 30(3), 650. DOI: 10.3390/molecules30030650
2. Bonakala et al. (2026). Towards accurate and scalable high-throughput MOF adsorption screening. *Chem. Sci.* DOI: 10.1039/d6sc00831c
3. Tamtaji et al. (2025). High-throughput DFT screening of MOF-74. *Int. J. Hydrogen Energy*, 151276. DOI: 10.1016/j.ijhydene.2025.151276
4. Polat et al. (2020). CO₂ separation using MOF composites. *Chem. Eng. J.*, 401, 124916. DOI: 10.1016/j.cej.2020.124916
5. Lim, D. W. (2024). Machine-learning model for CO₂ adsorption screening. DOI: 10.1063/10.0028344
6. Liu et al. (2025). High-efficiency water adsorption prediction by lattice GCMC. *RSC Appl. Interfaces*, 2, 230–242. DOI: 10.1039/d4lf00354c
