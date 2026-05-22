# 高次元パラメータ空間でのベイズ最適化フレームワーク — 実験レポート

**日付**: 2026-05-23  
**ステータス**: DRAFT — NOT FOR DISTRIBUTION  
**フレームワーク**: BOTorch 0.17.2 / GPyTorch 1.15.2 / Ax 1.2.4  
**乱数シード**: 42（再現性保証）

---

## 1. 実験目的と背景

本プロジェクトでは、高次元パラメータ空間における実験計画を効率化するベイズ最適化（Bayesian Optimization; BO）フレームワークを設計・実装した。従来のグリッドサーチやランダムサーチでは、次元数が増加するにつれて評価回数が指数関数的に増大する「次元の呪い」問題が深刻化する。ベイズ最適化は、ガウス過程（GP）による代理モデルと獲得関数を組み合わせることで、少ない評価回数で最適解に到達するデータ効率の高い手法である。

本フレームワークは以下の6つの要素技術を統合的に実装・評価した：

1. **GPカーネル選択と超パラメータ最適化** — 問題構造に適応したカーネル選択
2. **獲得関数（EI, UCB, KG）の比較** — 問題依存の選択基準の確立
3. **バッチ最適化** — 並列実験提案による実験効率化
4. **多目的ベイズ最適化（EHVI）** — 複数目的関数の同時最適化
5. **高次元次元削減統合（REMBO, HeSBO）** — 20変数以上への対応
6. **化学反応条件最適化** — 実応用ケーススタディ

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 ガウス過程回帰とカーネル選択

5種類のカーネルを交差検証（5-fold CV）により比較した：

| カーネル | 数学的表現 | 特性 |
|---------|-----------|------|
| RBF (SE) | $k(x,x') = \sigma^2 \exp(-\|x-x'\|^2 / 2l^2)$ | 無限回微分可能、滑らかな関数向け |
| Matérn 5/2 | $k = \sigma^2(1 + \sqrt{5}r + \frac{5}{3}r^2)\exp(-\sqrt{5}r)$ | 2回微分可能、物理シミュレーション向け |
| Matérn 3/2 | $k = \sigma^2(1 + \sqrt{3}r)\exp(-\sqrt{3}r)$ | 1回微分可能、やや粗い関数向け |
| Rational Quadratic | $k = \sigma^2(1 + r^2/2\alpha l^2)^{-\alpha}$ | マルチスケール構造に適応 |
| RBF + Periodic | 加法カーネル | 周期成分を含む関数向け |

全カーネルで **ARD（Automatic Relevance Determination）** を使用し、各次元ごとの長さスケールを独立に学習した。超パラメータの最適化には周辺尤度最大化（Type-II MLE）を使用した。

### 2.2 獲得関数

| 獲得関数 | アルゴリズム | 計算コスト |
|---------|------------|-----------|
| **EI** (Expected Improvement) | $\alpha(x) = \mathbb{E}[\max(f(x) - f^*, 0)]$ | 低（解析解） |
| **UCB** (Upper Confidence Bound) | $\alpha(x) = \mu(x) + \beta\sigma(x)$ | 低 |
| **KG** (Knowledge Gradient) | 次ステップのベイズ最適値の期待改善量 | 高（Fantasy必要） |

問題依存の選択ヒューリスティックを実装（ノイズ水準、予算、次元数、バッチサイズに基づくスコアリング）。

### 2.3 バッチ最適化

並列実験提案のため3手法を実装：

- **q-EI**: Monte Carlo近似による複数点同時最適化
- **q-UCB**: UCBのバッチ拡張（Sobol QMCサンプリング）
- **Kriging Believer**: GP平均予測を仮想観測として逐次的にバッチ構築

### 2.4 多目的ベイズ最適化

**q-EHVI（Expected Hypervolume Improvement）** を実装。非劣解集合のハイパーボリューム指標を最大化する。FastNondominatedPartitioning による効率的な分割を使用。

### 2.5 高次元次元削減

| 手法 | アプローチ | 特性 |
|------|-----------|------|
| **REMBO** | ランダム線形射影 $x = Az$ | 低次元部分空間での最適化 |
| **HeSBO** | ハッシュベース次元削減 | メモリ効率が高い |
| **Vanilla BO** | 全次元でのGP | ベースライン |
| **Random Search** | ランダムサンプリング | 比較基準 |

### 2.6 化学反応最適化

6パラメータの反応条件空間（温度、圧力、触媒量、溶媒比、滞留時間、pH）において、収率と選択性を最適化した。反応モデルにはArrhenius型温度依存性、平衡論的圧力効果、触媒の飽和挙動などの物理化学的知見を組み込んだ。

---

## 3. 主要な結果と数値

### 3.1 カーネル比較結果

| カーネル | 平均NLPD | 標準偏差 | 順位 |
|---------|---------|---------|------|
| **Matérn 5/2** | **0.8085** | **0.2322** | **1** |
| Matérn 3/2 | 0.8202 | 0.1852 | 2 |
| RBF | 0.8358 | 0.3382 | 3 |
| RQ | 0.8368 | 0.3394 | 4 |
| RBF+Periodic | 0.9933 | 0.1895 | 5 |

**Matérn 5/2** が最良のCV性能を達成。ARD長さスケール解析により、次元3（$l=0.203$）と次元5（$l=0.309$）が最も感度が高い（短い長さスケール＝高感度）パラメータであることが判明した。次元4（$l=68.307$）はほぼ不活性であった。

![Kernel Comparison](figures/kernel_comparison.png)
![Lengthscales](figures/lengthscales.png)

### 3.2 獲得関数比較結果（Hartmann-6）

| 獲得関数 | 最終最良値（平均±SD） | 3試行 |
|---------|---------------------|-------|
| **UCB** | **3.0605 ± 0.1419** | 探索と活用のバランス |
| EI | 2.8191 ± 0.0247 | 安定だがやや保守的 |

UCBが平均的により高い最良値を達成したが、分散も大きい。EIはより安定した収束を示した（Hartmann-6の最適値: 3.3224）。

**問題依存選択の推奨例**:
- 低ノイズ（σ=0.01）、中予算（n=50）→ **EI** 推奨
- 高ノイズ（σ=0.3）、小予算（n=15）、バッチ（q=4）→ **KG** 推奨
- 中ノイズ（σ=0.1）、高次元（d=25）→ **UCB** 推奨

![Acquisition Comparison](figures/acquisition_comparison.png)

### 3.3 バッチ最適化結果（batch_size=4, Hartmann-6）

| 手法 | 最終最良値 | 計算時間 |
|------|-----------|---------|
| **q-EI** | **3.1236** | **20.9s** |
| q-UCB | 3.1041 | 44.6s |
| Kriging Believer | 2.7885 | 83.2s |

q-EIが最良の性能・効率バランスを達成。Kriging Believerは逐次的なGP再構築のため計算コストが高く、性能も劣った。

![Batch Comparison](figures/batch_comparison.png)

### 3.4 多目的ベイズ最適化結果

- **最終ハイパーボリューム**: 3.96
- **パレート最適解数**: 18個（15初期 + 25反復から）
- ハイパーボリュームは反復とともに単調増加し、安定した収束を確認

![Pareto Front](figures/pareto_front.png)

### 3.5 高次元BO結果（D=25, 有効次元=6）

| 手法 | 最終最良値 | 最適値比（%） |
|------|-----------|-------------|
| **Vanilla BO** | **3.2721** | **98.5%** |
| REMBO | 2.9780 | 89.6% |
| HeSBO | 2.5949 | 78.1% |
| Random Search | 1.2936 | 38.9% |

意外な結果として、Vanilla BOがREMBOを上回った。これは以下の理由による：
- 有効次元が6次元と相対的に低いため、25次元でもGPが適応可能
- REMBOのランダム射影行列が最適部分空間を完全に捕捉できていない可能性
- n_init=15, n_iter=40（計55評価）は25次元としては少ないが、有効次元の低さが補償

ランダムサーチとの差は全手法で顕著であり、BO手法の有効性を確認した。

![High-Dim Comparison](figures/highdim_comparison.png)

### 3.6 化学反応条件最適化結果

#### 単目的最適化（収率最大化）

- **最良収率**: 72.17%
- **必要実験回数**: 55回（15初期 + 40反復）

**最適反応条件**:

| パラメータ | 最適値 | 単位 | 物理的解釈 |
|-----------|--------|------|-----------|
| 温度 | 140.18 | °C | Arrhenius最適温度近傍に収束 |
| 圧力 | 50.0 | bar | 上限に到達（平衡シフト） |
| 触媒量 | 10.0 | mol% | 上限に到達（飽和前） |
| 溶媒比 | 0.62 | v/v | 最適混合比に収束 |
| 滞留時間 | 90.11 | min | 十分な反応時間 |
| pH | 7.04 | - | 中性付近の最適窓 |

#### 多目的最適化（収率 × 選択性）

- **ハイパーボリューム**: 4392.65
- **パレート最適解**: 14個
- 収率と選択性のトレードオフを明確に可視化

![Chemical Convergence](figures/chemical_convergence.png)
![Chemical Pareto](figures/chemical_pareto.png)

---

## 4. 考察と今後の展望

### 4.1 主要な知見

1. **カーネル選択の重要性**: Matérn 5/2が多くの実応用で推奨される。ARD長さスケールは特徴量重要度の代理指標として活用でき、実験パラメータの感度分析に直結する。

2. **獲得関数の使い分け**: 問題特性に応じた選択が性能に大きく影響する。低ノイズ・十分な予算ではEI、探索重視やβチューニング可能な場合はUCB、高ノイズ・バッチ設定ではKGが推奨される。

3. **バッチ最適化**: q-EIが計算効率・性能のバランスで最優秀。実験室での並列実験提案に直接適用可能。

4. **高次元対応**: 有効次元が低い場合、標準BOでも比較的良好な性能を示す。REMBOは射影行列の品質に依存するため、複数回の独立実行と結果の併合が推奨される。

5. **化学反応最適化**: 55回の実験で72%の収率を達成。物理化学的に妥当な条件（温度≈140°C、中性pH）に収束しており、BOの実用性を実証した。

### 4.2 限界と注意事項

- **合成データによる評価**: 実実験データでの検証が必要
- **GPのスケーラビリティ**: O(n³)の計算量制約。n > 500では近似GP（SVGP、VariationalGP）の導入が必要
- **KG未計測**: 計算コストが高く、本実験ではEI/UCBのみを定量比較。KGの優位性はバッチ+高ノイズ設定で特に発揮される
- **SAASBO未実装**: 完全ベイズ推論（NUTS）によるスパースGPは高次元で有望だが、計算時間の制約から今回は除外

### 4.3 今後の展望

1. **Transfer Learning**: 類似反応系からの事前知識の転移学習によるwarm-start
2. **制約付き最適化**: 安全性・コスト制約を考慮した制約付きBO
3. **SAASBO統合**: 高次元でのスパース軸整列カーネルによる自動次元選択
4. **実験自動化連携**: ラボオートメーション（OT-2等）との統合による完全自動BO-実験ループ
5. **qLogEHVI**: 数値安定性を改善したlog変換版EHVIの採用（BOTorch推奨）

---

## 5. 生成ファイル一覧

### コードモジュール

| ファイル | 説明 |
|---------|------|
| `bayesopt_framework/__init__.py` | パッケージ初期化 |
| `bayesopt_framework/kernel_selection.py` | GPカーネル選択・CV比較 |
| `bayesopt_framework/acquisition_functions.py` | 獲得関数比較・問題依存選択 |
| `bayesopt_framework/batch_optimization.py` | バッチBO（q-EI, q-UCB, Kriging Believer） |
| `bayesopt_framework/multi_objective.py` | 多目的BO（q-EHVI） |
| `bayesopt_framework/high_dimensional.py` | 高次元BO（REMBO, HeSBO） |
| `bayesopt_framework/chemical_optimization.py` | 化学反応最適化ケーススタディ |
| `bayesopt_framework/visualization.py` | 可視化モジュール |
| `run_experiments.py` | 全実験実行スクリプト |

### 結果ファイル

| ファイル | 説明 |
|---------|------|
| `results/all_results.json` | 全実験の数値結果（JSON） |

### 図表

| ファイル | 説明 |
|---------|------|
| `figures/kernel_comparison.png` | カーネルCV比較（棒グラフ） |
| `figures/lengthscales.png` | ARD長さスケール（特徴量重要度） |
| `figures/acquisition_comparison.png` | 獲得関数収束曲線 |
| `figures/batch_comparison.png` | バッチ手法比較 |
| `figures/pareto_front.png` | 多目的最適化パレートフロント＋HV収束 |
| `figures/highdim_comparison.png` | 高次元BO手法比較 |
| `figures/chemical_convergence.png` | 化学反応収率の収束曲線 |
| `figures/chemical_pareto.png` | 収率-選択性パレートフロント |

### ログ

| ファイル | 説明 |
|---------|------|
| `logs/process-log.jsonl` | 実行トレース（JSONL） |

---

## 参考文献

1. Balandat, M. et al. (2020). BoTorch: A Framework for Efficient Monte-Carlo Bayesian Optimization. NeurIPS.
2. Snoek, J. et al. (2012). Practical Bayesian Optimization of Machine Learning Algorithms. NeurIPS.
3. Wang, Z. et al. (2016). Bayesian Optimization in a Billion Dimensions via Random Embeddings. JAIR.
4. Daulton, S. et al. (2020). Differentiable Expected Hypervolume Improvement for Parallel Multi-Objective BO. NeurIPS.
5. Eriksson, D. & Jankowiak, M. (2021). High-Dimensional Bayesian Optimization with Sparse Axis-Aligned Subspaces. UAI.
6. Ament, S. et al. (2023). Unexpected Improvements to Expected Improvement for BO. arXiv:2310.20708.
