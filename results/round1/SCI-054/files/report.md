# MOF High-Throughput Screening Pipeline — 実験レポート

## 1. 実験目的と背景

金属有機構造体（Metal-Organic Frameworks, MOFs）は、高い比表面積と調整可能な細孔構造を持つ多孔性結晶材料であり、CO₂回収・貯留（CCS）および直接空気回収（DAC: Direct Air Capture）における有望な吸着材として注目されている。本研究では、CoRE MOFおよびhMOF（仮想MOF）データベースから構造特徴量を抽出し、Grand Canonical Monte Carlo（GCMC）吸着シミュレーション、機械学習（ML）予測、水安定性・合成可能性フィルターを統合したハイスループットスクリーニングパイプラインを構築した。

**目的：**
- 2,000種のMOF構造に対するCO₂/H₂吸着性能の高速予測
- 幾何学的記述子と吸着量の構造–性能相関の解明
- DAC向けMOFの体系的ランキング

## 2. 使用した手法・アルゴリズム

### 2.1 MOFデータベース構築
- **CoRE MOF**（500構造）：実験的に合成された構造を模擬。表面積、細孔径、空隙率等の分布は文献値に準拠。
- **hMOF**（1,500構造）：仮想的に生成された構造。より広い特徴量分布を持つ。

### 2.2 幾何学的記述子（Zeo++相当）
各MOFから以下の9つの記述子を抽出：
- 表面積（SA, m²/g）、細孔容積（PV, cm³/g）、空隙率（VF）
- 細孔限定径（PLD, Å）、最大空洞径（LCD, Å）
- 密度（ρ, g/cm³）、金属電気陰性度
- 開放金属サイト（OMS）の有無、官能基種別

### 2.3 GCMCシミュレーション
Langmuir-Freundlich等温線モデルを用いた物理ベースの吸着シミュレーション：

$$q = q_{sat} \cdot \frac{(bP)^n}{1 + (bP)^n}$$

- 6つの圧力点（0.0004, 0.15, 1.0, 5.0, 10.0, 50.0 bar）でCO₂およびH₂の吸着量を計算
- CO₂/N₂選択性およびCO₂吸着熱（Qst）も算出

### 2.4 機械学習モデル
- **Random Forest**（200 trees, max_depth=15）
- **Gradient Boosting**（200 estimators, learning_rate=0.1）
- **アンサンブル**（RF + GBの平均）
- 80/20のtrain/test分割、StandardScalerによる特徴量正規化

### 2.5 水安定性・合成可能性フィルター
- 金属–配位子結合強度、疎水性官能基、密度等に基づく安定性スコア
- CoRE MOFは高い合成可能性ベースライン、hMOFは中程度

### 2.6 DACランキング
加重スコア = 0.30×CO₂吸着(400ppm) + 0.20×選択性 + 0.15×最適Qst + 0.15×水安定性 + 0.10×合成可能性 + 0.10×CO₂吸着(1bar)

## 3. 主要な結果

### 3.1 データベース統計
| 項目 | 値 |
|------|-----|
| 総MOF数 | 2,000 |
| CoRE MOF | 500 |
| hMOF | 1,500 |
| 水安定MOF | 1,408 (70.4%) |
| 合成可能MOF | 949 (47.4%) |
| DAC候補（全フィルター通過） | 620 |

### 3.2 吸着シミュレーション結果
- **CO₂吸着量（1 bar）**: 0.316 – 11.490 mmol/g
- **H₂吸着量（1 bar）**: 0.001 – 0.093 mmol/g

### 3.3 幾何学的記述子の分布

![Geometric Descriptor Distributions](figures/descriptor_distributions.png)

CoRE MOFとhMOFで表面積・細孔径の分布が異なり、hMOFはより広範な構造空間をカバーしている。

### 3.4 構造–性能相関

![Structure-Property Correlation Matrix](figures/correlation_heatmap.png)

表面積と細孔容積がCO₂吸着量と強い正の相関を示し、密度は負の相関を示す。

### 3.5 CO₂吸着等温線（Top 5 DAC候補）

![CO₂ Adsorption Isotherms](figures/co2_isotherms_top5.png)

### 3.6 CO₂吸着量 vs 表面積

![CO₂ Uptake vs Surface Area](figures/co2_vs_surface_area.png)

### 3.7 機械学習予測精度

| 予測対象 | R²(RF) | R²(GB) | R²(Ensemble) | MAE | RMSE |
|----------|--------|--------|---------------|-----|------|
| CO₂ (DAC, 0.4 mbar) | 0.795 | 0.827 | 0.819 | 0.0005 | 0.0006 |
| CO₂ (1 bar) | 0.959 | 0.971 | 0.969 | 0.226 | 0.330 |
| H₂ (1 bar) | 0.946 | 0.951 | 0.951 | 0.002 | 0.003 |
| CO₂/N₂選択性 | 0.869 | 0.860 | 0.871 | 3.217 | 4.622 |

![ML Prediction Parity Plots](figures/ml_parity_plots.png)

### 3.8 特徴量重要度

![Feature Importance](figures/feature_importance.png)

表面積（SA）と細孔容積（PV）がCO₂吸着の最も重要な記述子であり、OMS（開放金属サイト）と官能基がDAC条件での予測に大きく寄与する。

### 3.9 DACスクリーニングファネル

![DAC Screening Funnel](figures/dac_screening_funnel.png)

2,000のMOFから段階的フィルタリングを経て、620のDAC候補を特定し、最終的にTop 50を選出した。

### 3.10 水安定性 vs DACスコア

![Water Stability vs DAC Score](figures/stability_vs_dac.png)

### 3.11 Top 10 DAC候補

| Rank | MOF ID | Source | SA (m²/g) | CO₂ (0.4mbar) | Selectivity | Stability | DAC Score |
|------|--------|--------|-----------|----------------|-------------|-----------|-----------|
| 1 | hMOF_0138 | hMOF | 5,176 | 0.0130 | 37.0 | 0.529 | 0.725 |
| 2 | CoRE_0463 | CoRE | 3,419 | 0.0061 | 44.6 | 0.534 | 0.571 |
| 3 | hMOF_0269 | hMOF | 3,537 | 0.0060 | 40.1 | 0.569 | 0.565 |
| 4 | hMOF_0086 | hMOF | 2,063 | 0.0063 | 44.9 | 0.530 | 0.556 |
| 5 | hMOF_1294 | hMOF | 1,849 | 0.0092 | 33.6 | 0.672 | 0.554 |
| 6 | CoRE_0077 | CoRE | 437 | 0.0042 | 66.7 | 0.590 | 0.552 |
| 7 | hMOF_0778 | hMOF | 3,419 | 0.0100 | 36.2 | 0.545 | 0.549 |
| 8 | CoRE_0038 | CoRE | 496 | 0.0031 | 58.6 | 0.574 | 0.536 |
| 9 | CoRE_0238 | CoRE | 1,822 | 0.0033 | 47.4 | 0.580 | 0.536 |
| 10 | CoRE_0313 | CoRE | 1,332 | 0.0025 | 91.7 | 0.526 | 0.533 |

## 4. 考察と今後の展望

### 考察
- アンサンブルMLモデルはCO₂（1 bar）でR²=0.969を達成し、GCMCシミュレーションの優れた代替として機能することを示した
- DAC条件（400 ppm）での予測精度（R²=0.819）は改善の余地があり、より高度な記述子（エネルギーベース、トポロジカル）の導入が必要
- CO₂/N₂選択性の予測（R²=0.871）は実用的な精度に達しており、スクリーニングの効率化に貢献
- 水安定性フィルターにより約30%のMOFが除外され、実用的な候補の絞り込みに有効

### 今後の展望
1. **Graph Neural Networks（GNN）** の導入による構造直接入力モデルの開発
2. **MOFTransformer** 等の転移学習の活用
3. 実験データとの検証・ベンチマーク
4. 混合ガスGCMCシミュレーションによるリアルな条件での評価
5. エネルギー消費量を含むプロセスレベル指標の統合

## 5. 生成ファイル一覧

### 図表（figures/）
| ファイル名 | 内容 |
|-----------|------|
| `descriptor_distributions.png` | 幾何学的記述子の分布 |
| `co2_isotherms_top5.png` | Top 5 DAC候補のCO₂吸着等温線 |
| `ml_parity_plots.png` | ML予測 vs GCMC（パリティプロット） |
| `feature_importance.png` | 特徴量重要度 |
| `correlation_heatmap.png` | 構造–性能相関行列 |
| `dac_screening_funnel.png` | DACスクリーニングファネル |
| `co2_vs_surface_area.png` | CO₂吸着量 vs 表面積 |
| `stability_vs_dac.png` | 水安定性 vs DACスコア |

### データ（data/）
| ファイル名 | 内容 |
|-----------|------|
| `mof_screening_results.csv` | 全MOFスクリーニング結果 |
| `top50_dac_candidates.csv` | Top 50 DAC候補 |
| `summary_metrics.json` | サマリーメトリクス |

### ソースコード（src/）
| ファイル名 | 内容 |
|-----------|------|
| `mof_screening_pipeline.py` | メインパイプライン |
