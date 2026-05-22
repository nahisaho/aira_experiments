# DRAFT — NOT FOR DISTRIBUTION

# 脳オルガノイド大量培養のためのバイオリアクター設計と最適化

**作成日**: 2026-05-23  
**バージョン**: 1.0  
**ステータス**: シミュレーション設計・解析完了

---

## 目次

1. [実験目的と背景](#1-実験目的と背景)
2. [使用した手法・アルゴリズムの概要](#2-使用した手法アルゴリズムの概要)
3. [主要な結果と数値](#3-主要な結果と数値)
4. [考察と今後の展望](#4-考察と今後の展望)
5. [生成したファイル一覧](#5-生成したファイル一覧)

---

## 1. 実験目的と背景

### 1.1 研究目的

脳オルガノイドの大量培養に適した灌流型バイオリアクターの設計最適化を、計算流体力学（CFD）シミュレーションおよび数理モデリングを通じて実施する。具体的には以下の6つの設計課題に対してシミュレーションフレームワークを構築した：

1. 灌流型バイオリアクターの流体力学解析（CFD）
2. 酸素・栄養素輸送の反応-拡散モデリング
3. せん断応力と組織成熟度の関係モデリング
4. 培地組成の時間プログラム最適化
5. スケーラビリティ設計（バッチ→灌流→連続）
6. バイオマーカーモニタリング戦略

### 1.2 背景

脳オルガノイドは、ヒト脳の発生・疾患モデルとして創薬スクリーニングや基礎研究に不可欠である。しかし、従来の静置培養では以下の制約がある：

- **酸素・栄養素の拡散限界**: オルガノイド中心部の壊死（直径 > 1 mm で顕著）
- **スケーラビリティの欠如**: 手動操作に依存した低スループット
- **成熟度のばらつき**: 環境条件の不均一性による品質の変動

灌流型バイオリアクターはこれらの課題を解決する有力なアプローチであるが、設計パラメータ（流速、せん断応力、培地供給プロファイル）の最適化には体系的なシミュレーション手法が必要である。

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 CFD 流体力学シミュレーション

#### バイオリアクター形状

| パラメータ | 値 |
|---|---|
| 容器直径 | 80 mm |
| 容器高さ | 120 mm |
| 入出口ポート径 | 6 mm |
| オルガノイドバスケット直径 | 60 mm |
| バスケット高さ | 80 mm |
| バスケット多孔質度 | 0.4 |
| 透過率 | 1×10⁻¹⁰ m² |

#### 計算手法

- **2D 軸対称 Navier-Stokes 方程式**（定常解）を有限差分法で離散化
- **Darcy-Forchheimer モデル**でオルガノイドバスケットの多孔質領域をモデリング
- 培養液の物性：密度 1007 kg/m³、粘度 0.001 Pa·s
- 流量条件：0.5, 1.0, 2.0, 5.0 mL/min の4条件でパラメトリックスイープ
- OpenFOAM 互換のケースファイル（blockMeshDict, controlDict 等）を自動生成
- COMSOL Multiphysics 用のセットアップガイドも作成

### 2.2 酸素・栄養素輸送モデリング

**反応-拡散方程式**: ∂C/∂t = D∇²C − R(C)

球座標系での定常・非定常解を数値的に求解した。

| 物質 | 拡散係数 D (m²/s) | Vmax (mol/m³·s) | Km (mol/m³) |
|---|---|---|---|
| 酸素 (O₂) | 2.5×10⁻⁹ | 5×10⁻³ | 4.6×10⁻³ |
| グルコース | 6.7×10⁻¹⁰ | 1.2×10⁻² | 0.5 |

- **Michaelis-Menten 型消費速度式**で酸素・グルコースの消費をモデル化
- オルガノイド半径 0.5〜2.0 mm の範囲で濃度プロファイルを計算
- **臨界半径**（中心O₂濃度が 0.01 mol/m³ 以下となる半径）を外部O₂濃度・Vmaxの関数として算出

### 2.3 せん断応力-成熟度関係モデリング

- せん断応力 0.001〜1.0 Pa の範囲で成熟度指標をモデル化
- **神経マーカー発現**（MAP2, TUJ1, GFAP）: シグモイド応答関数
- **細胞生存率**: 高せん断域で減少（閾値 ~0.5 Pa）
- **多目的最適化**: 成熟度最大化 × 細胞損傷最小化の Pareto 最適解を算出

### 2.4 培地最適化

- 4段階の分化プロトコル（Neural induction → Patterning → Cortical differentiation → Maturation）
- 7成分（bFGF, EGF, BDNF, レチノイン酸, Matrigel, グルコース, O₂張力）の時間プロファイルを最適化
- **ベイズ最適化**による複合成熟スコアの最大化
- コスト制約を含む実現可能性評価

### 2.5 スケーラビリティ解析

- 3モード（バッチ 50 mL / 灌流 500 mL / 連続 5 L）の比較設計
- **無次元解析**: Reynolds 数, Damköhler 数, Peclet 数, 単位体積あたり動力

### 2.6 バイオマーカーモニタリング

- 15種バイオマーカーの時系列プロファイル生成
- **オンライン計測**（pH, DO, グルコース, 乳酸）: 連続センサー
- **Shewhart/CUSUM 管理図**による工程管理
- 異常検知アルゴリズムによるプロセス偏差検出

---

## 3. 主要な結果と数値

### 3.1 CFD シミュレーション結果

4つの流量条件における主要な流動特性を以下に示す：

| 流量 (mL/min) | 最大流速 (m/s) | 圧力損失 (Pa) | 壁面最大せん断 (Pa) | バスケット平均せん断 (Pa) |
|---|---|---|---|---|
| 0.5 | 1.38×10⁻⁵ | 0.00293 | 2.13×10⁻⁵ | 6.81×10⁻⁸ |
| 1.0 | 2.76×10⁻⁵ | 0.00587 | 4.27×10⁻⁵ | 1.36×10⁻⁷ |
| 2.0 | 5.52×10⁻⁵ | 0.01174 | 8.53×10⁻⁵ | 2.72×10⁻⁷ |
| 5.0 | 1.38×10⁻⁴ | 0.02935 | 2.13×10⁻⁴ | 6.81×10⁻⁷ |

**主要所見**: バスケット内のせん断応力はすべての条件で 10⁻⁷〜10⁻⁶ Pa と極めて低く、オルガノイドへの機械的損傷リスクは最小限である。圧力損失と流速は流量に線形比例し、層流条件が確認された。

![Velocity Field](figures/velocity_field.png)
![Pressure Field](figures/pressure_field.png)
![Shear Stress Distribution](figures/shear_stress_distribution.png)

### 3.2 酸素・栄養素輸送解析

#### 臨界半径

| 外部 O₂ (mol/m³) | Vmax 倍率 | 臨界半径 (mm) |
|---|---|---|
| 0.10 | 1.0× | 0.563 |
| 0.15 | 1.0× | 0.692 |
| **0.20** | **1.0×** | **0.799** |
| 0.25 | 1.0× | 0.891 |
| 0.30 | 1.0× | 0.975 |

**基準条件**（外部O₂ = 0.2 mol/m³, 基準Vmax）での臨界半径は **0.799 mm** であり、半径 0.8 mm 以上のオルガノイドでは中心部の低酸素壊死が予想される。

- 外部O₂濃度を 0.3 mol/m³ に上げれば臨界半径は ~0.975 mm に拡大
- 細胞密度（Vmax）を0.6倍に下げれば臨界半径は ~1.03 mm まで拡大可能

![Oxygen Radial Profile](figures/oxygen_radial_profile.png)
![Glucose Radial Profile](figures/glucose_radial_profile.png)
![Critical Radius Analysis](figures/critical_radius_analysis.png)
![Oxygen Time Evolution](figures/oxygen_time_evolution.png)

### 3.3 せん断応力-成熟度関係

パレート最適化の結果、以下の最適運転条件を同定した：

| 指標 | 値 |
|---|---|
| **最適せん断応力** | **0.0452 Pa** |
| 複合成熟スコア | 0.902 |
| 細胞生存率 | 0.970（97.0%）|
| 文献推奨範囲 | 0.01〜0.1 Pa |

最適せん断応力 0.0452 Pa は文献で報告されている至適範囲（0.01〜0.1 Pa）の中央付近に位置し、高い成熟スコア（0.902）と生存率（97.0%）を両立する。

![Shear-Maturation Response](figures/shear_maturation_response.png)
![Pareto Frontier](figures/pareto_frontier.png)
![Maturation Heatmap](figures/maturation_heatmap.png)

### 3.4 培地最適化プロファイル

| フェーズ | 期間 | bFGF | EGF | BDNF | RA (µM) | Matrigel (%) | Glucose (mM) | O₂ (%) | コスト ($/L) |
|---|---|---|---|---|---|---|---|---|---|
| Neural induction | Day 0–6 | 12 | 0 | 0 | 0 | 0 | 17 | 5.5 | 9.11 |
| Neural patterning | Day 6–25 | 24 | 12 | 5 | 0.12 | 1.2 | 16 | 8.0 | 65.11 |
| Cortical diff. | Day 25–50 | 8 | 4 | 28 | 0.55 | 0.8 | 11 | 12 | 50.03 |
| Maturation | Day 50–90 | 0 | 0 | 42 | 0.20 | 0 | 8 | 18 | 11.10 |

ベイズ最適化の結果、最適化プログラムはランダムプログラムに対して Cohen's d = 7.55（95% CI: Δscore [60.29, 92.96]）と有意に優位であった。

![Medium Temporal Program](figures/medium_temporal_program.png)
![Optimization Convergence](figures/optimization_convergence.png)
![Medium Cost Analysis](figures/medium_cost_analysis.png)

### 3.5 スケーラビリティ比較

| モード | 容量 | 収量/バッチ | 培地消費 (L/day) | コスト/個 ($) | サイズ均一性 (CV%) | 生存率 (%) |
|---|---|---|---|---|---|---|
| バッチ | 50 mL | 109 | 0.020 | 6.44 | 19.0 | 83 |
| 灌流 | 500 mL | 2,552 | 0.240 | 2.91 | 11.5 | 91 |
| 連続 | 5 L | 33,120 | 1.800 | 1.56 | 8.5 | 94 |

**スケールアップ因子**: 連続モードはバッチに対して **303倍** の生産性向上を達成し、コストは **75.8% 削減**（$6.44 → $1.56/個）された。

#### 無次元数解析

| 無次元数 | バッチ | 灌流 | 連続 |
|---|---|---|---|
| Reynolds | 18 | 3,000 | 22,320 |
| Damköhler | 0.141 | 0.082 | 0.046 |
| Peclet | 56.4 | 1,070 | 1,203 |

連続モードの Reynolds 数（22,320）は乱流域に入るため、適切なインペラ設計と流路制御が必要である。

![Scalability Comparison](figures/scalability_comparison.png)
![Dimensionless Analysis](figures/dimensionless_analysis.png)
![Cost Scaling](figures/cost_scaling.png)

### 3.6 バイオマーカーモニタリング

#### モニタリング階層

| 分類 | 計測項目 | 頻度 | 手法 |
|---|---|---|---|
| オンライン | pH, DO, グルコース, 乳酸 | 連続（30分間隔） | インラインセンサー |
| アットライン | LDH, サイトカインパネル | 毎日 | サンプリング→即時分析 |
| オフライン | OCT4, PAX6, MAP2, SYN1 等 | 毎週 | qPCR, 免疫染色, 電気生理 |

#### 品質管理

- Shewhart 管理図による即時逸脱検出
- CUSUM 管理図によるドリフト検出
- 多変量異常スコア > 4.5 でプロセス保留・確認

![Biomarker Timecourse](figures/biomarker_timecourse.png)
![Control Charts](figures/control_charts.png)
![Monitoring Decision Tree](figures/monitoring_decision_tree.png)

---

## 4. 考察と今後の展望

### 4.1 設計上の主要知見

1. **せん断応力の安全域**: CFD解析により、0.5〜5.0 mL/min のすべての灌流条件でバスケット内せん断応力は 10⁻⁷〜10⁻⁶ Pa であり、最適せん断（0.0452 Pa）の達成にはバスケット構造の透過率調整もしくは直接灌流方式への変更が必要である。

2. **酸素輸送がサイズ制約**: 臨界半径 0.799 mm は、直径約 1.6 mm 以上のオルガノイドで中心壊死が発生することを示す。実用的には (a) 外部O₂分圧の増加、(b) 灌流による対流輸送の強化、(c) 血管網の誘導 のいずれかが必要。

3. **スケールアップの実現性**: 連続灌流モード（5L）はバッチに対して 303 倍の生産性・76% のコスト削減を達成するが、Reynolds 数 22,320 は乱流遷移域であり、低せん断インペラ（錨型・ピッチドブレード）の採用と適切なバッフル設計が重要である。

4. **培地コスト**: Patterning フェーズ（$65.11/L）が最もコストが高く、bFGF・EGFの使用量最適化が全体コスト削減の鍵となる。

### 4.2 モデルの限界

- CFD シミュレーションは 2D 軸対称の簡略モデルであり、3D 効果（旋回流、バスケット支持構造の影響）は評価していない
- せん断応力-成熟度モデルは文献値に基づく現象論的モデルであり、実験的検証が必要
- 培地最適化は合成データに基づくデモンストレーションであり、実際の細胞応答データでの再最適化が前提
- スケーラビリティ比較の信頼区間はモンテカルロ摂動に基づく推定値であり、実機実験データではない

### 4.3 COMSOL/OpenFOAM 連携

- **OpenFOAM**: ケースファイル一式を `results/openfoam_cases/` に生成済み。blockMeshDict, controlDict, fvSchemes, fvSolution, 境界条件ファイルを含む
- **COMSOL**: 詳細なセットアップガイドを `results/comsol_setup.md` に記載。Laminar Flow (spf) + Transport of Diluted Species (tds) の2物理モジュール構成
- 両ソルバーのクロスバリデーション（同一条件での速度場・圧力場の比較）を推奨

### 4.4 今後の展望

1. **3D CFD 解析**: OpenFOAM/COMSOL による 3D フルジオメトリ解析の実施
2. **実験検証**: PIV（粒子画像速度測定）による流速場の実験的検証
3. **マルチフィジックス連成**: 流体-構造-物質輸送の完全連成解析
4. **機械学習統合**: 実培養データを用いた培地最適化モデルの更新
5. **血管誘導モデル**: 内皮細胞共培養による血管網形成のモデリング
6. **GMP 準拠設計**: 臨床応用に向けた GMP 準拠バイオリアクターへの展開
7. **デジタルツイン**: リアルタイムモニタリングデータとシミュレーションを連携したデジタルツインシステムの構築

---

## 5. 生成したファイル一覧

### シミュレーションコード

| ファイル | 内容 |
|---|---|
| `data/bioreactor_geometry.py` | バイオリアクター形状のパラメトリック生成スクリプト |
| `results/cfd_openfoam_setup.py` | OpenFOAM ケースファイル自動生成スクリプト |
| `results/cfd_simulation.py` | 2D 軸対称 CFD シミュレーション（有限差分法） |
| `results/oxygen_transport.py` | 酸素・栄養素の反応-拡散方程式ソルバー |
| `results/shear_maturation_model.py` | せん断応力-成熟度関係のモデリング・最適化 |
| `results/medium_optimization.py` | 培地組成のベイズ最適化 |
| `results/scalability_design.py` | スケーラビリティ比較・無次元解析 |
| `results/biomarker_monitoring.py` | バイオマーカーモニタリング・管理図 |

### データファイル

| ファイル | 内容 |
|---|---|
| `data/geometry_params.json` | バイオリアクター形状パラメータ（OpenFOAM メッシュ情報含む） |
| `data/biomarker_timecourse.csv` | 合成バイオマーカー時系列データ |
| `results/velocity_field.csv` | 速度場データ |
| `results/shear_stress.csv` | せん断応力分布データ |
| `results/cfd_summary.csv` | CFD 解析サマリー |
| `results/oxygen_profiles.csv` | 酸素濃度プロファイル |
| `results/nutrient_profiles.csv` | グルコース濃度プロファイル |
| `results/critical_radius_scan.csv` | 臨界半径パラメトリックスキャン結果 |
| `results/shear_maturation_data.csv` | せん断-成熟度応答データ |
| `results/pareto_frontier.csv` | パレート最適フロンティアデータ |
| `results/optimized_medium_profiles.csv` | 最適化培地プロファイル |
| `results/scalability_comparison.csv` | スケーラビリティ比較データ |

### 図表（すべて 300 DPI, colorblind-friendly パレット）

| ファイル | 内容 |
|---|---|
| `figures/velocity_field.png` | 速度場コンター図 |
| `figures/pressure_field.png` | 圧力場分布図 |
| `figures/shear_stress_distribution.png` | せん断応力分布図 |
| `figures/oxygen_radial_profile.png` | 酸素濃度の半径方向プロファイル |
| `figures/glucose_radial_profile.png` | グルコース濃度プロファイル |
| `figures/critical_radius_analysis.png` | 臨界半径の感度解析 |
| `figures/oxygen_time_evolution.png` | 酸素濃度の時間発展 |
| `figures/shear_maturation_response.png` | せん断-成熟度多パネル応答図 |
| `figures/pareto_frontier.png` | パレート最適フロンティア |
| `figures/maturation_heatmap.png` | 成熟スコアヒートマップ |
| `figures/medium_temporal_program.png` | 培地組成時間プログラム |
| `figures/optimization_convergence.png` | ベイズ最適化収束曲線 |
| `figures/medium_cost_analysis.png` | 培地コスト-品質トレードオフ |
| `figures/scalability_comparison.png` | スケーラビリティ比較棒グラフ |
| `figures/dimensionless_analysis.png` | 無次元数比較 |
| `figures/cost_scaling.png` | コストスケーリング曲線 |
| `figures/biomarker_timecourse.png` | バイオマーカー時系列プロファイル |
| `figures/control_charts.png` | Shewhart 管理図 |
| `figures/monitoring_decision_tree.png` | モニタリング判定フローチャート |

### ドキュメント

| ファイル | 内容 |
|---|---|
| `report.md` | 本レポート |
| `results/comsol_setup.md` | COMSOL Multiphysics セットアップガイド |
| `results/monitoring_protocol.md` | バイオマーカーモニタリングプロトコル |
| `results/statistical-summary.md` | 統計解析サマリー |
| `results/cfd_summary.md` | CFD 解析結果サマリー |
| `data/preprocessing-log.md` | データ前処理ログ |

### OpenFOAM ケースファイル

| ディレクトリ | 内容 |
|---|---|
| `results/openfoam_cases/` | 全流量条件のケースディレクトリ一式 |

### ログ

| ファイル | 内容 |
|---|---|
| `logs/process-log.jsonl` | 実行トレースログ |
| `logs/learnings-log.jsonl` | 学習記録 |

---

*本レポートはシミュレーション設計段階のものであり、実験的検証前の計算結果に基づいています。実機での検証と最適化の反復が推奨されます。*
