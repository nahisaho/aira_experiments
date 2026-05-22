# 大規模科学データの品質管理と異常検知 — 実験レポート

**DRAFT — NOT FOR DISTRIBUTION**
**作成日**: 2026-05-23
**バージョン**: 1.0.0

---

## 1. 実験目的と背景

### 目的
CERN/LIGO規模の大規模科学実験データに対応可能な、自動品質管理・異常検知パイプラインを設計・実装・検証する。特に以下の6つの要件を満たすシステムの構築を目的とする：

1. **時系列変化点検出**（PELT/BOCPD）
2. **多変量外れ値検出**（Isolation Forest/Deep SVDD）
3. **物理的制約を組み込んだ異常スコアリング**
4. **ドリフト検出**と自動モデル再訓練トリガー
5. **説明可能な異常検知**（異常原因の自動特定）
6. **ストリーミング対応アーキテクチャ**（CERN/LIGO型）

### 背景
高エネルギー物理学実験（LHC）や重力波検出（LIGO）では、毎秒数千万イベントが生成される。このデータストリームからリアルタイムに異常を検出し、物理的に妥当な信号を選別するパイプラインが不可欠である。本実験では、合成粒子物理データを用いてパイプライン全体の有効性を検証した。

---

## 2. 使用した手法・アルゴリズム

### 2.1 変化点検出（Module 1）

| アルゴリズム | タイプ | 計算量 | 用途 |
|---|---|---|---|
| **PELT** (Pruned Exact Linear Time) | オフラインバッチ | O(n) | 蓄積データの後方解析 |
| **BOCPD** (Bayesian Online Changepoint Detection) | オンラインストリーミング | O(n) per update | リアルタイム変化点検出 |

- PELT: RBFカーネルを用いたセグメント分割。ペナルティ項は BIC (log(n)·σ²) で自動設定
- BOCPD: Student-t予測分布を用いたベイズ更新。ハザードレート λ=1/250

### 2.2 多変量外れ値検出（Module 2）

| 手法 | 特徴 | contamination |
|---|---|---|
| **Isolation Forest** | アンサンブル木ベース、高次元に強い | 0.05 |
| **Deep SVDD** | ニューラルネット（ReLU AE proxy）、一クラス分類 | ν=0.05 |

- Isolation Forest: 200本の木、置換特徴量重要度を同時計算
- Deep SVDD: 8→4次元エンコーダ、50エポック学習、中心からの距離でスコアリング

### 2.3 物理制約付きスコアリング（Module 3）

統計的異常スコアと物理法則違反スコアを重み付き結合：

```
combined = w_stat × normalize(stat_score) + w_phys × normalize(phys_score)
```

実装した物理制約ライブラリ：
- エネルギー保存則 (E_in ≈ E_out)
- 運動量保存則 (|p_total| ≈ 0)
- 値域制約 (pT < 500 GeV, |η| < 5)
- 正定値制約 (質量 > 0, エネルギー > 0)
- 因果律制約 (Δr ≤ c·Δt)

### 2.4 ドリフト検出（Module 4）

| 手法 | アプローチ | 感度 |
|---|---|---|
| **ADWIN** | 適応的ウィンドウ分割 | δ=0.01 (高感度) |
| **Page-Hinkley** | 累積和統計量 | threshold=30 |

- **RetrainingTrigger**: 性能劣化（>5%低下）、累積ドリフト（3回連続検出）、急激な性能崩壊で再訓練を発火

### 2.5 説明可能な異常検知（Module 5）

1. **置換特徴量重要度**: 各特徴量をシャッフルしスコア変化を測定
2. **局所説明**: 正常分布からのz-scoreで各異常の主要寄与特徴を特定
3. **決定ルール抽出**: サロゲート決定木（max_depth=4）で解釈可能なIF-THENルールを生成
4. **根本原因分析**: 異常サンプルをz-score支配特徴量でクラスタリング

### 2.6 ストリーミングパイプライン（Module 6）

4層階層アーキテクチャ（CERN LHCトリガーシステム準拠）：

| 層 | レイテンシ | データ削減 | 実装 |
|---|---|---|---|
| **L0: Hardware Trigger** | <1μs | 40MHz → 100kHz (400×) | FPGA/ASIC |
| **L1: Online Filter** | <10ms | 100kHz → 1kHz (100×) | z-score, EWMA |
| **L2: Nearline Analysis** | <1s | 1kHz → 10Hz (100×) | Isolation Forest, PELT |
| **L3: Offline Deep** | min〜hours | 保存イベント | Deep SVDD, SHAP |

総データ削減率: **約10⁶倍**

---

## 3. 主要な結果と数値

### 3.1 変化点検出

| 指標 | PELT | BOCPD |
|---|---|---|
| 真の変化点 | 4 (t=500, 900, 1200, 1600) | 4 |
| 検出数 | **4** | 0 (閾値0.3) |
| 検出位置 | [500, 900, 1200, 1600] | — |
| 完全一致率 | **100%** | 0% |

**考察**: PELTはバッチ処理で4つの変化点を完全に検出。BOCPDは閾値 0.3 では検出できず、パラメータチューニングが必要。

![図1: 変化点検出](figures/fig1_changepoint_detection.png)

### 3.2 多変量外れ値検出

| 指標 | Isolation Forest | Deep SVDD |
|---|---|---|
| 検出数 | 250 | 390 |
| Precision | **0.672** | 0.415 |
| Recall | **0.672** | 0.648 |
| F1 score | **0.672** | 0.506 |
| 異常率 | 5.0% | 7.8% |

**特徴量重要度** (Isolation Forest):
| 特徴量 | 重要度 |
|---|---|
| pseudorapidity | 0.149 |
| transverse_momentum | 0.133 |
| track_isolation | 0.127 |
| azimuthal_angle | 0.126 |

![図2: 多変量外れ値検出](figures/fig2_multivariate_outliers.png)

### 3.3 物理制約付きスコアリング

| 制約 | 違反数 | 違反率 | 最大違反量 |
|---|---|---|---|
| energy_positive | 0 | 0.0% | 0.0 |
| pT_range (<500 GeV) | 1 | 0.02% | 118.4 |
| mass_positive | 0 | 0.0% | 0.0 |
| eta_range (\|η\|<5) | 63 | 1.26% | 3.84 |

- 統計スコアと物理スコアの相関: 制約違反は統計的異常と部分的に一致
- 物理制約の導入により、物理的に不可能なイベントを優先的にフラグ付け

![図3: 物理制約スコアリング](figures/fig3_physics_constraints.png)

### 3.4 ドリフト検出

| 手法 | 検出ドリフト数 |
|---|---|
| ADWIN | 136 (高感度) |
| Page-Hinkley | 3 (保守的) |
| 再訓練トリガー | **8回** 発火 |

- ADWINは高感度設定（δ=0.01）のため多数の微小ドリフトを検出
- Page-Hinkleyは保守的に3回の主要ドリフトのみを検出
- 再訓練トリガーは性能閾値とドリフト累積の両方を考慮し、8回の再訓練を推奨

![図4: ドリフト検出](figures/fig4_drift_detection.png)

### 3.5 説明可能な異常検知

**根本原因分析**:
| 主因特徴量 | 異常の割合 | 平均z-score | 方向 |
|---|---|---|---|
| transverse_momentum | 23.2% | 8.74 | 高 |
| track_isolation | 21.2% | 7.89 | 高 |
| missing_ET | 21.2% | 10.38 | 高 |
| energy | 18.8% | 10.45 | 高 |
| invariant_mass | 15.6% | 3.18 | 高 |

- サロゲート決定木から **9個の解釈可能なルール**を抽出
- 異常の主要因は transverse_momentum と missing_ET の高値（BSM信号の特徴と一致）

![図5: 説明可能な異常検知](figures/fig5_explainable_anomalies.png)

### 3.6 ストリーミングパイプライン

| 指標 | 値 |
|---|---|
| 処理イベント数 | 2,000 |
| 検出異常数 | 117 |
| 異常率 | 5.9% |
| スループット | **37,449 events/sec** |
| Critical アラート | 37 |
| Warning アラート | 20 |
| Info アラート | 60 |

![図6: アーキテクチャ設計](figures/fig6_architecture.png)

---

## 4. 考察と今後の展望

### 考察

1. **Isolation Forestの優位性**: F1=0.672でDeep SVDD (F1=0.506) を上回った。高次元データにおいてアンサンブル木ベースの手法が安定した性能を示す
2. **PELTの精度**: 4つの真の変化点を100%検出。バッチ処理が可能な場合はPELTが最も信頼性が高い
3. **物理制約の重要性**: η範囲制約で1.26%の物理的に非現実的なイベントを特定。統計的手法のみでは見逃されるアーティファクトを捕捉可能
4. **説明可能性**: 根本原因分析により、検出された異常の物理的意味（高pT・高missing ETはBSM信号の特徴）が自動的に解釈可能
5. **スループット**: 単一プロセスで約37,000 events/secを達成。分散処理により実運用レベル（>MHz）にスケール可能

### 制限事項

- Deep SVDDはPyTorchなしの軽量実装のため、深層モデルの本来の性能を反映していない
- BOCPDのハイパーパラメータ（ハザードレート、閾値）はデータ依存であり、自動チューニングが必要
- ストリーミングパイプラインのスループットは単一スレッドの結果であり、実運用では分散フレームワーク（Apache Flink等）が必要
- 合成データでの検証であり、実際のCERN/LIGOデータでの追加検証が必要

### 今後の展望

1. **PyTorch/JAXベースのDeep SVDD**: 本格的な深層ネットワークによる一クラス分類の実装
2. **連合学習**: 複数の検出器/実験施設にまたがるプライバシー保護型異常検知
3. **強化学習ベースの閾値最適化**: 動的な異常閾値の自動調整
4. **GNN (Graph Neural Networks)**: 検出器ジオメトリを活用した構造的異常検知
5. **Apache Flink/Kafka統合**: 本番環境向けストリーミング基盤の構築
6. **A/Bテストフレームワーク**: 新旧異常検知モデルのオンライン比較評価

---

## 5. 生成ファイル一覧

### ソースコード (`src/`)

| ファイル | 説明 |
|---|---|
| `src/__init__.py` | パッケージ初期化 |
| `src/changepoint_detection.py` | PELT/BOCPD変化点検出 |
| `src/multivariate_outlier.py` | Isolation Forest/Deep SVDD外れ値検出 |
| `src/physics_constraints.py` | 物理制約付きスコアリング |
| `src/drift_detection.py` | ADWIN/Page-Hinkleyドリフト検出 |
| `src/explainable_anomaly.py` | 説明可能な異常検知 |
| `src/streaming_pipeline.py` | ストリーミングパイプライン・アーキテクチャ設計 |

### 実験スクリプト

| ファイル | 説明 |
|---|---|
| `run_experiment.py` | 全実験の実行スクリプト |
| `generate_figures.py` | 全図表の生成スクリプト |

### データ (`data/`)

| ファイル | 説明 |
|---|---|
| `data/synthetic_events.csv` | 合成粒子物理イベントデータ (5,000件×8特徴量) |
| `data/timeseries.csv` | 変化点検出用時系列データ (2,000点) |

### 結果 (`results/`)

| ファイル | 説明 |
|---|---|
| `results/experiment_results.json` | 全実験の数値結果 |

### 図表 (`figures/`)

| ファイル | 説明 |
|---|---|
| `figures/fig1_changepoint_detection.png` | PELT/BOCPD変化点検出結果 |
| `figures/fig2_multivariate_outliers.png` | 多変量外れ値検出（6パネル） |
| `figures/fig3_physics_constraints.png` | 物理制約スコアリング |
| `figures/fig4_drift_detection.png` | ドリフト検出・再訓練トリガー |
| `figures/fig5_explainable_anomalies.png` | 説明可能な異常検知 |
| `figures/fig6_architecture.png` | CERN/LIGO型アーキテクチャ設計図 |

### ログ (`logs/`)

| ファイル | 説明 |
|---|---|
| `logs/process-log.jsonl` | 実行トレースログ |

---

*本レポートは自動生成されたドラフトです。実データへの適用前に、ドメイン専門家によるレビューを推奨します。*
