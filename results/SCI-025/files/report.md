# 生分解性ポリマーの分子設計フレームワーク：実験レポート

## 1. 実験目的と背景

本研究では、環境中で制御的に分解される生分解性ポリマーの分子設計フレームワークを開発した。具体的には、加水分解速度予測モデル、機械的性質と分解性のトレードオフ最適化、酵素分解のMichaelis-Mentenモデリング、海洋環境分解シミュレーション、コンビナトリアル共重合体設計、およびPLA/PHA/PBSの改質設計ケーススタディを統合的に実施した。さらに、分子記述子と機械学習を用いた構造-分解性関係モデルを構築した。

### 背景

プラスチック汚染は地球規模の環境問題であり、生分解性ポリマーの開発が急務である。しかし、分解性と機械的性質のトレードオフ、環境条件（温度、pH、微生物叢）による分解挙動の変動など、設計上の課題が多い。先行研究では、機械学習を用いたポリマー物性予測（Zhao et al., 2023）、高スループット実験による生分解性ポリエステルの発見（Fransen et al., 2023）、酵素分解のMichaelis-Mentenモデリングなどが報告されているが、これらを統合した設計フレームワークは未だ確立されていない。

## 2. 使用した手法・アルゴリズムの概要

### 2.1 加水分解速度予測モデル

Arrhenius型の速度式に主鎖結合種、結晶度、分子量の依存性を組み込んだ半経験的モデルを構築：

$$k_h = A \cdot f_{bond} \cdot \exp\left(-\frac{E_a}{RT}\right) \cdot (1 - X_c)^{\alpha} \cdot \left(\frac{M_w}{M_{ref}}\right)^{-\beta} \cdot g(pH)$$

### 2.2 機械的性質モデル

結晶度、分子量、架橋密度の関数として引張強度と弾性率を推定し、分解速度とのPareto最適化を実施。

### 2.3 Michaelis-Menten酵素分解モデル

酵素失活を考慮したODE系により、Proteinase K、PHA depolymerase、Lipase、Cutinaseによる分解をシミュレーション：

$$\frac{d[S]}{dt} = -\frac{V_{max} \cdot [E] \cdot [S]}{K_m + [S]}, \quad \frac{d[E]}{dt} = -k_d \cdot [E]$$

### 2.4 海洋環境シミュレーション

ポリマー→オリゴマー→モノマーの逐次分解とMonod型微生物成長を結合した4変数ODE系。

### 2.5 機械学習モデル

9種の分子記述子（log Mw、結晶度、結合因子、親水性、架橋密度、分岐度、Tg、表面積、多孔性）を特徴量とし、Random ForestおよびGradient Boostingで分解速度を予測。

## 3. 主要な結果と数値

### 3.1 加水分解速度の構造依存性

6種の主鎖結合について加水分解速度を比較した。オルトエステル結合が最も速く、アミド結合が最も遅い。結晶度の増加により分解速度は非線形に低下し、分子量の増加もk_hを減少させる。

![Figure 1: 加水分解速度モデル](figures/fig1_hydrolysis_rate.png)

### 3.2 機械的性質-分解性トレードオフ

500サンプルのランダム探索により、引張強度と分解速度のPareto front（赤星印）を同定した。高結晶度は強度を高めるが分解性を低下させる。

![Figure 2: トレードオフ最適化](figures/fig2_tradeoff.png)

### 3.3 Michaelis-Menten酵素分解

4種の酵素について分解動態を解析。PHA depolymeraseが最も高い初期速度を示し、Lipase (PBS)は緩やかだが持続的な分解パターンを示した。

![Figure 3: Michaelis-Menten酵素分解](figures/fig3_michaelis_menten.png)

### 3.4 海洋環境分解シミュレーション

4条件下（熱帯表層、温帯表層、深海、酸性化海洋）での分解挙動：

| 環境条件 | 半減期（日） |
|---|---|
| 熱帯表層 (28°C, pH 8.1) | 102.6 |
| 温帯表層 (15°C, pH 8.1) | 396.2 |
| 深海 (4°C, pH 7.8) | 1095.0 |
| 酸性化海洋 (15°C, pH 7.6) | 653.2 |

![Figure 4: 海洋分解シミュレーション](figures/fig4_marine_degradation.png)

![Figure 5: 海洋環境別半減期](figures/fig5_marine_halflife.png)

### 3.5 コンビナトリアル共重合体設計

6種モノマーのDirichlet分布からの2000組成のライブラリを生成し、分解速度-強度-コストの3次元設計空間を可視化した。

![Figure 6: コンビナトリアル設計](figures/fig6_combinatorial.png)

### 3.6 PLA/PHA/PBS改質ケーススタディ

各ポリマーの改質効果（ベースに対する相対値）：

**PLA系:**
- PLA + 10% PEG: 強度 0.75x, 分解 1.59x
- PLA + Nanoclay: 強度 1.35x, 分解 0.48x
- PLA-co-GA (90:10): 強度 0.61x, 分解 2.00x
- Stereocomplex PLA: 強度 1.48x, 分解 0.37x

**PHA系:**
- P(HB-co-HV) 80:20: 強度 0.71x, 分解 1.97x
- PHA + Chain Extender: 強度 1.07x, 分解 0.76x

**PBS系:**
- PBS-co-BF (70:30): 強度 0.57x, 分解 2.28x（最大分解加速）

![Figure 7: ケーススタディ](figures/fig7_case_studies.png)

### 3.7 機械学習モデル性能

| モデル | R² | RMSE | MAE | CV-R²（5-fold） |
|---|---|---|---|---|
| Random Forest | 0.9424 | 0.1833 | 0.1402 | 0.9386 ± 0.0100 |
| Gradient Boosting | 0.9670 | 0.1387 | 0.1070 | 0.9650 ± 0.0028 |

Gradient Boostingが最高性能を示した（R² = 0.967）。

![Figure 8: 機械学習モデル結果](figures/fig8_ml_model.png)

![Figure 9: 残差分析](figures/fig9_residuals.png)

特徴量重要度（Random Forest）上位5：
1. crosslink_density: 0.6029
2. bond_factor: 0.2185
3. crystallinity: 0.0733
4. log_Mw: 0.0587
5. porosity: 0.0142

![Figure 10: 相関マトリクス](figures/fig10_correlation.png)

## 4. 考察と今後の展望

### 考察

- 加水分解モデルは主鎖結合種を最重要因子として正しく捉えており、結晶度・分子量の寄与も先行研究と整合する。
- 機械的性質と分解性のPareto frontは、実用設計における意思決定の指針となる。
- 海洋シミュレーションでは、深海条件で半減期が約3年と極めて長く、環境条件が分解挙動を決定的に支配することが示された。
- MLモデルにおいて架橋密度が最も重要な特徴量となった。これは架橋が水分子や酵素のアクセスを物理的に阻害するためである。
- Gradient BoostingがRandom Forestを上回った要因として、逐次的な誤差補正による非線形関係の精密な捕捉が挙げられる。

### 今後の展望

1. 実験データによるモデル校正とバリデーション
2. Graph Neural Networkの導入による分子構造直接入力モデル
3. 実環境サンプリングによる微生物叢データの統合
4. 多目的最適化（NSGA-II等）による自動設計
5. ライフサイクルアセスメント（LCA）との統合

## 5. 生成したファイル一覧

| ファイル名 | 説明 |
|---|---|
| `experiment.py` | 全実験のPythonスクリプト |
| `ml_results.csv` | MLモデル性能結果 |
| `case_study_results.csv` | ケーススタディ数値結果 |
| `figures/fig1_hydrolysis_rate.png` | 加水分解速度モデル |
| `figures/fig2_tradeoff.png` | トレードオフ最適化 |
| `figures/fig3_michaelis_menten.png` | Michaelis-Menten分解 |
| `figures/fig4_marine_degradation.png` | 海洋分解シミュレーション |
| `figures/fig5_marine_halflife.png` | 海洋半減期比較 |
| `figures/fig6_combinatorial.png` | コンビナトリアル設計 |
| `figures/fig7_case_studies.png` | PLA/PHA/PBSケーススタディ |
| `figures/fig8_ml_model.png` | MLモデル結果 |
| `figures/fig9_residuals.png` | 残差分析 |
| `figures/fig10_correlation.png` | 構造-物性相関 |
| `report.md` | 本レポート |
| `paper.md` | 学術論文 |
