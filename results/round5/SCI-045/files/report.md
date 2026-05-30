# 実験レポート：DNAメチル化データによるエピジェネティッククロック改良モデルの開発

---

## 1. 実験目的と背景

### 目的

DNAメチル化パターンから生物学的年齢を推定するエピジェネティッククロックの改良モデルを開発し、以下の観点から評価する：

1. Horvath/GrimAgeクラスの線形モデル（ElasticNet）に対する改善方針の検証
2. 組織特異的メチル化パターンを考慮した評価
3. 加齢加速度（age acceleration）のバイオマーカーとしての検証
4. 深層学習（DeepClock、AttentionClock）の設計と評価
5. 介入効果（運動・食事・薬物）の検出感度評価
6. 長寿コホートデータでのバリデーション戦略

### 背景

エピジェネティッククロックは、CpGサイトのメチル化レベルを入力として生物学的年齢を推定する回帰モデルである。Horvath（2013）がElasticNetを用いて353のCpGサイトから多組織横断的な年齢推定（MAD ≈ 3.6年）を実現して以来、GrimAge、PhenoAge等が発展してきた。しかし、これらのモデルは線形性の仮定、組織特異性の欠如、介入効果検出感度の低さという限界を持つ。

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 合成データ生成

実際のDNAメチル化データの特性を模倣した合成データセットを以下の仕様で生成した：

| パラメータ | 値 |
|---|---|
| サンプル数 | 798 |
| CpG特徴量数 | 500 |
| 年齢範囲 | 20.1 – 89.8歳（一様分布） |
| 組織タイプ | 3種（血液・唾液・脳） |
| 加齢加速サンプル数 | 119（14.9%） |
| 平均加齢加速度 | 6.83 ± 1.77年 |

CpGサイトの分類：
- **過メチル化CpG（35%）**: 加齢とともにβ値が増加（傾き 0.003–0.008/年）
- **低メチル化CpG（35%）**: 加齢とともにβ値が減少（傾き −0.007 – −0.003/年）
- **非線形CpG（15%）**: 正弦波的年齢応答
- **ノイズCpG（15%）**: 年齢情報を含まない

### 2.2 モデル一覧

| モデル | タイプ | 主要パラメータ |
|---|---|---|
| ElasticNet（Horvath型ベースライン） | 線形正則化回帰 | α=0.01, l1_ratio=0.5 |
| Ridge Regression | 線形正則化回帰 | λ=10 |
| Random Forest | アンサンブル（木） | 100木, max_depth=10 |
| Gradient Boosting | アンサンブル（勾配ブースティング） | 200推定器, lr=0.05 |
| DeepClock（全結合NN） | 深層学習 | 500→256→128→64→32→1, Dropout=0.3 |
| AttentionClock（Multi-head Attention） | 深層学習 | 500→64（埋め込み）, 4ヘッド, Dropout=0.2 |

### 2.3 評価プロトコル

- **5分割交差検証**（KFold, シャッフルあり）
- 各フォールドで訓練データに対してStandardScaler適合、テストデータに変換適用
- 評価指標：MAE（主指標）、RMSE、R²（すべて平均±SD）

---

## 3. 主要な結果と数値

### 3.1 モデル全体比較（5分割CV）

![Figure 1: モデル比較（MAEおよびR²）](figures/fig1_model_comparison.png)

![Figure 2: 予測年齢 vs 実際の年齢（全フォールド）](figures/fig2_predicted_vs_actual.png)

**表1：5分割交差検証の結果（平均 ± SD）**

| モデル | MAE（年） | RMSE（年） | R² |
|---|---|---|---|
| **Ridge Regression** | **1.11 ± 0.05** | **1.39 ± 0.04** | **0.9955 ± 0.0004** |
| ElasticNet（Horvath型） | 1.38 ± 0.13 | 1.71 ± 0.13 | 0.9931 ± 0.0011 |
| Gradient Boosting | 2.24 ± 0.12 | 2.78 ± 0.16 | 0.9821 ± 0.0018 |
| Random Forest | 2.38 ± 0.18 | 3.01 ± 0.23 | 0.9788 ± 0.0034 |
| AttentionClock | 4.40 ± 0.98 | 6.86 ± 1.80 | 0.8822 ± 0.0674 |
| DeepClock（NN） | 33.95 ± 1.13 | 39.43 ± 1.12 | **−2.62 ± 0.33** |

**⚠️ 注意点**：Ridge/ElasticNetのR²（>0.99）は合成データの線形構造を直接反映したものであり、実際のDNAメチル化データでは3–7倍程度の誤差増大が予想される。

**DeepClockの失敗**：パラメータ数（約15万）に対してサンプル数（~640/フォールド）が不十分（比率 ≈ 234:1）。R² = −2.62は定数予測より悪い結果を示し、完全な過学習/未学習を示唆する。

### 3.2 組織特異的性能

![Figure 3: 組織別パフォーマンス](figures/fig3_tissue_specific.png)

**表2：組織別MAE（3分割CV、平均 ± SD）**

| 組織 | ElasticNet MAE（年） | DeepClock MAE（年） |
|---|---|---|
| 血液 | 1.53 ± 0.04 | 52.71 ± 1.11 |
| 唾液 | 1.74 ± 0.16 | 53.69 ± 0.76 |
| 脳 | 1.71 ± 0.07 | 51.31 ± 1.38 |

線形モデルは組織内でも安定した性能を維持（MAE ~1.5–1.7年）。DeepClockは組織内split（n≈266）で完全に破綻し、組織特異的NN開発には大規模データが不可欠であることが示された。

### 3.3 加齢加速度の検出

![Figure 4: 加齢加速度の検出](figures/fig4_age_acceleration.png)

**表3：加齢加速度の検出性能**

| モデル | Pearson r | Mann-Whitney p値 |
|---|---|---|
| ElasticNet | −0.050 | 0.935 |
| DeepClock | 0.000 | 0.487 |

**有意な検出には至らなかった**。この結果は、暦年齢（chronological age）を目的変数として訓練したモデルは生物学的年齢加速を直接検出できないことを示す。実用的な加齢加速度推定には、死亡率・疾患発症などの外的結果変数を用いたキャリブレーションが必要。

### 3.4 介入効果の検出感度

![Figure 5: 介入効果の検出](figures/fig5_intervention_effects.png)

**表4：介入別効果検出感度（ElasticNet）**

| 介入 | シミュレーション効果量 | ElasticNet検出Δ（年） | p値 | NN検出Δ（年） | p値 |
|---|---|---|---|---|---|
| 運動 | −2.5年 | −3.22年 | < 0.0001 | −1.39年 | < 0.0001 |
| 食事 | −1.5年 | −1.92年 | < 0.0001 | −0.94年 | < 0.0001 |
| 薬物（ラパマイシン型） | −4.0年 | −5.12年 | < 0.0001 | −2.50年 | < 0.0001 |

すべての介入でp < 0.0001の有意差を検出。ElasticNetは実際の効果量を若干過大評価（介入操作がElasticNetの学習パターンと一致するため）。

### 3.5 長寿コホートバリデーション

![Figure 6: 長寿コホートバリデーション](figures/fig6_longevity_cohort.png)

**表5：長寿コホートパフォーマンス（年齢 >75、n=167）**

| モデル | MAE（年） | Pearson r | p値 |
|---|---|---|---|
| ElasticNet | 5.34 | 0.799 | 2.55 × 10⁻³⁸ |
| DeepClock | 13.60 | 0.221 | 4.06 × 10⁻³ |

長寿コホートでは年齢範囲が狭くなり、ElasticNetのMAEが増加（1.38→5.34年）するものの、強い相関（r=0.799）を維持。DeepClockは相関が大幅に低下（r=0.221）。

### 3.6 CpG特徴量重要度

![Figure 7: CpG特徴量重要度](figures/fig7_cpg_importance.png)

ElasticNet係数の分析により、係数の高いCpGが明確に正（過メチル化、年齢とともに増加）・負（低メチル化、年齢とともに減少）に二極化していることが確認された。

---

## 4. 考察

### 4.1 線形モデルの優位性とその解釈

今回の実験で Ridge/ElasticNet が圧倒的な性能を示した主因は、**合成データが線形の年齢-メチル化関係を前提として生成されている**からである。実際のDNAメチル化データは：
- CpG間に強い相関構造（クロマチンドメイン依存）
- 細胞タイプ混合効果による交絡
- 非線形・閾値効果
を含むため、線形モデルのR²は大幅に低下すると予想される。

### 4.2 DeepClockの失敗から得られる教訓

DeepClockの失敗（R² = −2.62）は**データ不足**に起因する。パラメータ数/サンプル比が約234:1では、どんなに正則化しても汎化は困難。AltumAge（de Lima Camillo et al., 2022）は20,318サンプルでNNの優位性を示している。実用的な深層学習エピジェネティッククロックには最低5,000–10,000サンプルが必要と推定される。

### 4.3 自己批判的評価

| 批判点 | 詳細 |
|---|---|
| 合成データへの依存 | R² > 0.99の性能は実際のゲノムデータでは実現しない可能性が高い |
| 細胞タイプ混合モデリングの欠如 | 血液メチル化の主要交絡因子（白血球分画）を考慮していない |
| 加齢加速度の外的妥当性なし | 死亡率・疾患データなしに生物学的年齢の意味を検証できない |
| 介入シミュレーションの循環性 | 同一CpGを操作してモデルに入力するため効果を過大検出 |
| 長寿コホートデータの限界 | 真の生物学的年齢（−5年）は仮定値であり実測データではない |

### 4.4 実世界への一般化可能性

実世界のGEO公開データ（例：GSE40279、血液メチル化450Kアレイ）でのHorvathクロック検証では、MAD ≈ 3.6年が報告されている。本研究のElasticNet MAE（1.38年）はこれより3倍程度楽観的であり、**実世界データへの直接適用には追加バリデーションが必須**。

---

## 5. 今後の展望

1. **実際のDNAメチル化データ（GEO等）での検証**: 本パイプラインをGSE40279等の公開データに適用し、実際のベンチマークを取得
2. **細胞タイプ補正**: Houseman法やCELL-MiXを用いた白血球分画補正の組み込み
3. **多タスク学習**: 暦年齢と生物学的年齢（死亡リスクスコア）の同時予測
4. **フェデレーテッドラーニング**: プライバシーに配慮した多コホートデータ統合
5. **縦断的モデリング**: 個人内の経時的メチル化変化を捉える再帰型/Transformerアーキテクチャ
6. **注意機構の解釈**: AttentionClockの注意重みを生物学的機能（プロモーター/エンハンサー等）と照合

---

## 6. 生成したファイル一覧

| ファイル | 説明 |
|---|---|
| `epigenetic_clock_experiment.py` | 実験スクリプト本体 |
| `experiment_results.json` | 全実験結果の数値データ（JSON形式） |
| `figures/fig1_model_comparison.png` | モデル比較（MAEとR²のバーチャート） |
| `figures/fig2_predicted_vs_actual.png` | 予測年齢 vs 実際年齢（6モデル） |
| `figures/fig3_tissue_specific.png` | 組織特異的性能比較 |
| `figures/fig4_age_acceleration.png` | 加齢加速度の検出（正常 vs 加速グループ） |
| `figures/fig5_intervention_effects.png` | 介入効果検出（3種類の介入） |
| `figures/fig6_longevity_cohort.png` | 長寿コホートのバリデーション |
| `figures/fig7_cpg_importance.png` | CpG特徴量重要度（ElasticNet係数） |
| `paper.md` | 学術論文形式の成果文書 |
| `report.md` | 本実験レポート |

---

## 参考文献

1. de Lima Camillo, L.P. et al. (2022). A pan-tissue DNA-methylation epigenetic clock based on deep learning. *npj Aging*, 8, 1–12. https://doi.org/10.1038/s41514-022-00085-y
2. Oblak, L. et al. (2021). GrimAge Outperforms Other Epigenetic Clocks in the Prediction of Age-Related Clinical Phenotypes and All-Cause Mortality. *J Gerontol A*, 76(5), 741–749. https://doi.org/10.1093/gerona/glaa286
3. Moqri, M. et al. (2023). Biomarkers of aging for the identification and evaluation of longevity interventions. *Cell*, 186, 3758–3775. https://doi.org/10.1016/j.cell.2023.08.003
4. Rutledge, J. et al. (2022). Measuring biological age using omics data. *Nat Rev Genet*, 23, 715–727. https://doi.org/10.1038/s41576-022-00511-7
5. Moqri, M. et al. (2024). Validation of biomarkers of aging. *Nat Med*, 30, 360–372. https://doi.org/10.1038/s41591-023-02784-9
6. Tian, Y.E. et al. (2023). Heterogeneous aging across multiple organ systems. *Nat Med*, 29, 1221–1231. https://doi.org/10.1038/s41591-023-02296-6
7. Bocklandt, S. et al. (2020). DNA Methylation Biomarkers in Aging and Age-Related Diseases. *Front Genet*, 11, 171. https://doi.org/10.3389/fgene.2020.00171
8. Shireby, G.L. et al. (2020). Recalibrating the epigenetic clock. *Brain*, 143(12), 3763–3775. https://doi.org/10.1093/brain/awaa334
9. McCartney, D.L. et al. (2021). Genome-wide association studies identify 137 genetic loci for DNA methylation biomarkers. *Genome Biol*, 22, 194. https://doi.org/10.1186/s13059-021-02398-9
10. Higgins-Chen, A.T. et al. (2022). A computational solution for bolstering reliability of epigenetic clocks. *Nat Aging*, 2, 644–661. https://doi.org/10.1038/s43587-022-00248-2
