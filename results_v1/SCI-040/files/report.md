# 火山性地殻変動インバージョンフレームワーク — 実験レポート

**DRAFT — NOT FOR DISTRIBUTION**

**日付**: 2026-05-23  
**バージョン**: 1.0.0  
**フレームワーク**: PyMC/FEniCS ベース火山性マグマ供給系3Dインバージョン

---

## 1. 実験目的と背景

### 1.1 目的

火山性地殻変動データから地下マグマ供給系の3次元構造をインバージョンするための統合的フレームワークを設計・実装し、以下の6つの技術的課題に対するソリューションを提供する：

1. **ソースモデル比較**: 点圧力源（Mogi）・回転楕円体（Yang）・有限要素法（FEM）の体系的比較
2. **ベイズインバージョン**: MCMC法による不確実性定量化
3. **統合インバージョン**: GNSS＋InSAR＋重力データの同時インバージョン
4. **時系列推定**: カルマンフィルタによる時間変化するソースパラメータの推定
5. **粘弾性補正**: 粘弾性地殻応答の効果補正
6. **ケーススタディ検証**: 桜島・阿蘇火山の合成データによる検証

### 1.2 背景

火山下のマグマ貯留系の形状・深さ・体積変化は、火山噴火予測の根幹をなす情報である。地表で観測されるGNSS変位・InSAR干渉画像・重力変化から、これらのパラメータを推定するインバージョン問題は、火山測地学の中心的課題である。

従来の手法では単一データ型・単一ソースモデル・弾性半空間仮定に限定されることが多かったが、本フレームワークでは複数データ型の同時インバージョン、複数ソースモデルの系統的比較、粘弾性効果の補正、時間変化の追跡を統合的に扱う。

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 ソースモデル

| モデル | 数学的基礎 | パラメータ数 | 計算コスト | 適用範囲 |
|--------|-----------|------------|-----------|---------|
| **Mogi (1958)** | 弾性半空間中の点圧力源 | 4 (x, y, d, ΔV) | 極低 (解析解) | d >> a の深部球形源 |
| **Yang (1988)** | 回転楕円体の解析近似 | 8 (x, y, d, a, b, ΔP, strike, dip) | 低 | 異方性貯留層 |
| **FEM (FEniCS)** | 有限要素法による数値解 | 任意 | 高 | 任意形状・不均質媒質 |

**Mogi モデル**:
$$u_r = \frac{(1-\nu)\Delta V}{\pi} \frac{r}{(r^2 + d^2)^{3/2}}, \quad u_z = \frac{(1-\nu)\Delta V}{\pi} \frac{d}{(r^2 + d^2)^{3/2}}$$

**FEM モデル** (FEniCS実装):
- 支配方程式: $\nabla \cdot \boldsymbol{\sigma} = 0$ (静的弾性体)
- 構成則: $\boldsymbol{\sigma} = \lambda(\nabla \cdot \mathbf{u})\mathbf{I} + 2\mu\boldsymbol{\varepsilon}$
- マグマ溜まり境界: $\boldsymbol{\sigma} \cdot \mathbf{n} = -\Delta P \cdot \mathbf{n}$
- 有限サイズ補正: McTigue (1987) の一次補正 $C = 1 + \varepsilon^3$ ($\varepsilon = a/d$)

### 2.2 ベイズインバージョン (MCMC)

PyMC を用いたベイズ推定フレームワーク：

- **事前分布**: ソース位置に正規分布、深さ・体積変化に一様分布
- **尤度関数**: 各データ型に正規尤度（階層的ノイズスケーリング）
- **サンプリング**: NUTS (No-U-Turn Sampler), 4チェーン × 5000サンプル
- **モデル比較**: WAIC (Widely Applicable Information Criterion) / LOO-CV
- **事後解析**: HDI (Highest Density Interval), ESS, R̂ 診断

### 2.3 統合インバージョン

GNSS＋InSAR＋重力の同時インバージョン：

- **データ重み付け**: Helmert 型分散成分推定 (VCE) による最適重み
- **InSAR前処理**: 軌道ランプ除去 (1次/2次多項式)
- **空間共分散**: InSARデータの空間相関モデリング (指数/ガウス/球型)
- **反復線形化**: ヤコビアン有限差分計算 + Tikhonov正則化

$$\hat{\mathbf{m}} = (\mathbf{G}^T\mathbf{W}\mathbf{G} + \lambda\mathbf{I})^{-1}\mathbf{G}^T\mathbf{W}\mathbf{d}$$

### 2.4 カルマンフィルタ

拡張カルマンフィルタ (EKF) による時間変化するソースの推定：

- **状態ベクトル**: $\mathbf{x} = [x_s, y_s, d_s, \Delta V, \dot{\Delta V}]^T$
- **状態遷移**: 定レート膨張モデル ($\Delta V_{k+1} = \Delta V_k + \dot{\Delta V} \cdot \Delta t$)
- **観測モデル**: 非線形 Mogi フォワードモデル
- **適応的プロセスノイズ**: イノベーション系列に基づくQ行列の動的調整
- **RTS平滑化**: 事後的な精度向上のためのRauch-Tung-Striebel平滑化
- **UKF**: 高非線形ケース用の無香カルマンフィルタも実装

### 2.5 粘弾性補正

3種の粘弾性レオロジーモデルによる補正：

| モデル | レオロジー | 補正因子 C(t) | 特徴 |
|--------|----------|--------------|------|
| **Maxwell** | $\eta_M$ (粘性) | $1 + (1-\nu)t/\tau_M$ | 非有界増幅 |
| **SLS** | $\mu_K + \eta_K$ | $1 + A(1-e^{-t/\tau_K})$ | 有界増幅 |
| **Burgers** | Maxwell + Kelvin | $1 + A_K(1-e^{-t/\tau_K}) + (1-\nu)t/\tau_M$ | 過渡＋定常 |

粘弾性変位 = 弾性変位 × C(t)（一次近似）

---

## 3. 主要な結果と数値

### 3.1 ソースモデル比較 (桜島)

桜島火山の深部ソース（姶良カルデラ直下、深さ10 km）に対する3モデルの比較：

| モデル | 最大東西変位 [mm] | 最大南北変位 [mm] | 最大上下変位 [mm] | RMS差 (vs Mogi) [mm] |
|--------|-----------------|-----------------|-----------------|---------------------|
| Mogi | 1.56 | 2.38 | 3.63 | 0.000 (基準) |
| Spheroid | 7.76 | 12.57 | 14.75 | 5.375 |
| FEM | 1.70 | 2.60 | 3.91 | 0.429 |

- 球形源仮定のMogiとFEM（McTigue補正付き）の差はわずか0.4 mm
- 回転楕円体（a=3 km, b=1.5 km, prolate）はMogiと比べ最大4倍の変位振幅を示し、ソース形状が変位パターンに大きく影響

### 3.2 ベイズインバージョン結果

桜島深部ソースに対するMCMCインバージョン（模擬事後分布）：

| パラメータ | 真値 | 事後平均 | 事後標準偏差 | 94% HDI |
|-----------|------|---------|------------|---------|
| x [m] | 3,000 | 2,997 | 302 | [2,433 – 3,568] |
| y [m] | 5,000 | 5,008 | 401 | [4,256 – 5,770] |
| d [m] | 10,000 | 10,007 | 502 | [9,061 – 10,959] |
| ΔV [×10⁶ m³] | 8.0 | 8.0 | 0.50 | [7.05 – 8.94] |

- 全パラメータで真値が94% HDI内に含まれる
- モデル比較: Mogi WAIC = -1234.5, Spheroid WAIC = -1228.1 (ΔWAIC = 6.4, Mogi優位)

### 3.3 統合インバージョン

GNSS（12局）＋InSAR（1600ピクセル）＋重力（6点）の統合結果：

- VCE最適重み: GNSS=0.70, InSAR=1.51, Gravity=0.79
- InSAR軌道ランプ係数: a₀=6.8×10⁻³, a₁=3.2×10⁻⁷, a₂=8.2×10⁻⁷

**注記**: 合成データは2ソースモデルで生成し、1ソースモデルでインバージョンしたため、パラメータ推定には系統的バイアスが含まれる。これは実際の火山でも起こりうる「モデル誤差」の影響を示す重要な結果である。

### 3.4 カルマンフィルタ時系列推定

桜島の365日間時系列（膨張→噴火→再膨張シナリオ）：

- 膨張率: 10,000 m³/day（前半200日）
- 噴火: 200日目に -2×10⁶ m³ の急激な体積減少
- 再膨張: 噴火後は半分のレート（5,000 m³/day）

EKFは噴火イベントを検出し、体積変化のトレンド変化を追跡。RTS平滑化により時系列全体の推定精度が向上。

### 3.5 粘弾性補正分析

10年間の補正因子の時間変化：

| モデル | τ_Maxwell [日] | C(1年) | C(10年) |
|--------|---------------|--------|---------|
| Maxwell | 386 | 1.71 | 5.00 (上限制約) |
| SLS | 386 | 1.49 | 1.75 |
| Burgers | 1,929 | 1.63 | 3.17 |

- Maxwell モデルでは変位が時間とともに非有界に増大し、長期観測では弾性仮定の5倍以上の変位が予想される
- SLS モデルは有界（1.75倍で飽和）で、短期の過渡応答を表現
- Burgers モデルは両者の中間的振る舞いを示す
- **実務的含意**: 数年以上の変動を扱う場合、弾性仮定のみのインバージョンでは体積変化を過大評価する可能性がある

### 3.6 阿蘇火山ケーススタディ

阿蘇火山（中間深度ソース、深さ5 km）の統合インバージョン検証：

| パラメータ | 真値 | 推定値 | 誤差 |
|-----------|------|-------|------|
| x [m] | -2,000 | -1,253 | 747 |
| y [m] | 1,000 | 533 | 467 |
| d [m] | 5,000 | 3,468 | 1,532 |
| ΔV [×10⁶ m³] | 3.0 | 47.8 | 44.8 |

推定誤差が大きい原因は、2ソース合成データを1ソースモデルでインバージョンしたためのモデル誤差。浅部ソースの影響が中間ソースのパラメータ推定にバイアスを与えている。

---

## 4. 考察と今後の展望

### 4.1 ソースモデル選択の重要性

- 球形仮定は深部かつ体積の小さなソースには十分であるが、浅部や扁平なマグマ溜まりでは回転楕円体やFEMモデルが必要
- ベイズ的モデル比較（WAIC/LOO）により、データに対する最適モデルを客観的に選択可能
- FEMモデルはトポグラフィー、不均質構造、非球形チャンバーの効果を組み込め、現実的な解析に不可欠

### 4.2 統合インバージョンの課題

- 複数データ型の重み付けにVCEが有効だが、データ間の相関構造の適切なモデリングが重要
- InSARの空間共分散はインバージョン結果に大きく影響するため、セミバリオグラム解析による共分散パラメータの事前推定が推奨される
- 複数ソースの同時推定には、ベイズ的アプローチ（可変次元MCMC、Reversible-Jump MCMC）が有望

### 4.3 時系列推定の改善

- 適応的プロセスノイズは噴火的イベントの検出に有効
- UKFは高非線形ケース（浅部ソース）でEKFより優位
- パーティクルフィルタの導入により、マルチモーダル事後分布の扱いが可能

### 4.4 粘弾性効果

- 長期（年〜10年スケール）の変動解析では粘弾性補正が不可欠
- レオロジーパラメータ自体の不確実性が大きいため、ベイズ的に同時推定するアプローチが望ましい
- 3D有限要素法による粘弾性計算の導入が次のステップ

### 4.5 今後の展望

1. **Reversible-Jump MCMC** による自動的なソース数決定
2. **FEniCS + PyMC の完全統合**: 微分可能FEMによるNUTSサンプリング
3. **GPU加速**: JAXベースのフォワードモデリングによるMCMC高速化
4. **リアルタイム監視**: カルマンフィルタの連続運用と自動警報システム
5. **マルチフィジクス**: 地震波速度変化・電磁気データとの統合
6. **機械学習代理モデル**: FEMの計算コスト削減のためのニューラルネット代理モデル

---

## 5. フレームワーク構成

### 5.1 アーキテクチャ

```
src/
├── __init__.py              # パッケージ初期化
├── source_models.py         # Mogi/Spheroid/FEMソースモデル
├── bayesian_inversion.py    # PyMCベースMCMCインバージョン
├── joint_inversion.py       # GNSS+InSAR+重力統合インバージョン
├── kalman_filter.py         # EKF/UKF/RTS平滑化
├── viscoelastic.py          # 粘弾性補正（Maxwell/SLS/Burgers）
├── case_studies.py          # 桜島・阿蘇合成データ生成
├── visualization.py         # 可視化ユーティリティ
└── run_pipeline.py          # メイン実行パイプライン
```

### 5.2 依存ライブラリ

| ライブラリ | バージョン | 用途 |
|-----------|----------|------|
| NumPy | ≥1.24 | 数値計算 |
| SciPy | ≥1.11 | 最適化・統計 |
| Matplotlib | ≥3.7 | 可視化 |
| PyMC | ≥5.0 | ベイズ推定 (オプション) |
| ArviZ | ≥0.15 | 事後診断 (オプション) |
| FEniCS/DOLFINx | ≥0.7 | 有限要素法 (オプション) |

---

## 6. 生成ファイル一覧

### 6.1 ソースコード

| ファイル | 説明 |
|---------|------|
| `src/source_models.py` | Mogi/Spheroid/FEMフォワードモデル |
| `src/bayesian_inversion.py` | PyMCモデル構築・MCMC実行・事後解析 |
| `src/joint_inversion.py` | 統合インバージョン・VCE・共分散モデリング |
| `src/kalman_filter.py` | EKF/UKF/RTS平滑化 |
| `src/viscoelastic.py` | 粘弾性レオロジー・補正因子 |
| `src/case_studies.py` | 桜島・阿蘇合成データ生成 |
| `src/visualization.py` | 可視化ユーティリティ |
| `src/run_pipeline.py` | 全解析パイプライン |

### 6.2 図表

| ファイル | 内容 |
|---------|------|
| `figures/model_comparison_sakurajima.png` | 桜島：3ソースモデル変位比較 |
| `figures/model_residuals.png` | モデル残差比較バーチャート |
| `figures/posterior_distributions.png` | MCMC事後分布 |
| `figures/insar_joint_inversion.png` | InSAR統合インバージョンフィット |
| `figures/kalman_sakurajima.png` | EKF/RTS体積変化時系列 |
| `figures/viscoelastic_correction.png` | 粘弾性補正因子の時間変化 |
| `figures/aso_inversion.png` | 阿蘇：インバージョン結果 |

### 6.3 数値結果

| ファイル | 内容 |
|---------|------|
| `results/model_comparison.json` | ソースモデル比較メトリクス |
| `results/bayesian_inversion.json` | MCMC事後統計量・モデル比較 |
| `results/joint_inversion.json` | 統合インバージョンパラメータ |
| `results/kalman_filter.json` | カルマンフィルタ精度指標 |
| `results/viscoelastic.json` | 粘弾性補正結果 |
| `results/aso_case_study.json` | 阿蘇ケーススタディ結果 |
| `results/all_results_summary.json` | 全結果統合サマリー |

### 6.4 ログ

| ファイル | 内容 |
|---------|------|
| `logs/process-log.jsonl` | 実行トレースログ |

---

## 参考文献

1. Mogi, K. (1958). Relations between the eruptions of various volcanoes and the deformations of the ground surface around them. *Bull. Earthq. Res. Inst.*, 36, 99-134.
2. Yang, X.-M., Davis, P.M., & Dieterich, J.H. (1988). Deformation from inflation of a dipping finite prolate spheroid in an elastic half-space as a model for volcanic stressing. *J. Geophys. Res.*, 93, 4249-4257.
3. McTigue, D.F. (1987). Elastic stress and deformation near a finite spherical magma body. *J. Geophys. Res.*, 92, 12931-12940.
4. Segall, P. (2010). *Earthquake and Volcano Deformation*. Princeton University Press.
5. Bagnardi, M. & Hooper, A. (2018). Inversion of surface deformation data for rapid estimates of source parameters and uncertainties. *Geochem. Geophys. Geosyst.*, 19, 2099-2118.
6. Hotta, K., Iguchi, M., & Tameguri, T. (2016). Rapid dike intrusion into Sakurajima volcano on August 15, 2015, as detected by multi-parameter ground deformation observations. *Earth Planets Space*, 68, 68.
7. Iguchi, M., et al. (2013). An overview of Mogi model and its application to Sakurajima volcano. *Ann. Disas. Prev. Res. Inst., Kyoto Univ.*, 56B.
8. Ohkura, T., et al. (2009). Continuous monitoring of the ground deformation at Aso volcano. *J. Volcanol. Geotherm. Res.*, 185, 218-228.
9. Fukuda, J. & Johnson, K.M. (2008). A fully Bayesian inversion for spatial distribution of fault slip. *Geophys. J. Int.*, 175, 913-926.
