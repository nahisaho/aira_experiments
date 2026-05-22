# 都市交通ミクロシミュレーションとリアルタイム制御最適化の統合システム設計

**DRAFT — NOT FOR DISTRIBUTION**

**日付**: 2026-05-23  
**著者**: Co-Scientist  
**システム**: SUMO/Flow/RLlib ベース統合交通シミュレーション最適化フレームワーク

---

## 1. 実験目的と背景

### 1.1 目的

東京都心部（丸の内・日本橋エリア、3km×3km）を対象に、ミクロ交通シミュレーションとマルチエージェント強化学習（MARL）による信号制御最適化を統合したフレームワークを設計・実装する。以下の6要素を統一的に扱うシステムアーキテクチャを構築することを目的とする：

1. 車両挙動モデル（IDM/MOBIL）のパラメータ化とSUMO統合
2. 交差点信号制御のMAPPO（Multi-Agent PPO）による最適化
3. マルチモーダル交通（車・バス・自転車・歩行者）の統合
4. プローブデータを活用したリアルタイム交通需要推定
5. 事故・工事時の動的リルーティング
6. 東京都心3km四方のケーススタディ

### 1.2 背景

都市交通の効率化は、渋滞緩和、CO₂排出削減、公共交通の信頼性向上において重要課題である。従来の固定時間制御や単純なアクチュエーテッド制御では、時々刻々変化する交通需要への適応が困難である。近年、マルチエージェント強化学習（MARL）がスケーラブルな適応型信号制御として注目されているが、需要推定やインシデント対応との統合は十分に研究されていない。

本研究では、SUMOミクロシミュレータをベースに、Flow制御フレームワークとRay/RLlibを統合し、エンドツーエンドの交通最適化パイプラインを構築する。

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 システムアーキテクチャ

![System Architecture](figures/fig1_system_architecture.png)

システムは以下の6モジュールから構成される：

| モジュール | 技術 | ファイル |
|---|---|---|
| 車両挙動モデル | IDM + MOBIL | `src/models/idm_model.py` |
| 信号制御 | MAPPO (CTDE) | `src/agents/marl_signal_control.py` |
| 需要推定 | 拡張カルマンフィルタ | `src/models/demand_estimation.py` |
| 動的リルーティング | A* + Yen's K-shortest | `src/models/dynamic_routing.py` |
| SUMO連携 | TraCI / Flow | `src/network/sumo_environment.py` |
| 統合オーケストレータ | パイプライン制御 | `src/main_orchestrator.py` |

### 2.2 車両挙動モデル（IDM + MOBIL）

**Intelligent Driver Model (IDM)** を縦方向の追従行動モデルとして採用。加速度は以下の式で計算される：

```
dv/dt = a × [1 - (v/v₀)^δ - (s*(v,Δv)/s)²]
```

ここで `s*(v,Δv) = s₀ + vT + vΔv/(2√(ab))` は希望車間距離である。

| パラメータ | 乗用車 | バス | 自転車 |
|---|---|---|---|
| v₀ (希望速度, m/s) | 13.89 (50km/h) | 11.11 (40km/h) | 4.17 (15km/h) |
| T (安全車頭時間, s) | 1.5 | 2.0 | 1.0 |
| a (最大加速度, m/s²) | 1.4 | 1.0 | 1.2 |
| b (快適減速度, m/s²) | 2.0 | 1.5 | 2.5 |
| s₀ (最小車間, m) | 2.0 | 3.0 | 1.0 |

![IDM Acceleration](figures/fig5_idm_acceleration.png)

車線変更には **MOBIL** (Minimizing Overall Braking Induced by Lane changes) モデルを使用。安全基準と礼譲パラメータに基づく車線変更判断を行う。

**キャリブレーション**: Levenberg-Marquardt法による非線形最小二乗法で、実測軌跡データからパラメータ推定を行う機能を実装。

### 2.3 交差点信号制御（MAPPO）

**Multi-Agent Proximal Policy Optimization (MAPPO)** による分散型信号制御を設計。

**アーキテクチャ**: Centralized Training with Decentralized Execution (CTDE)
- 各交差点が独立エージェント（48エージェント）
- パラメータ共有による効率的学習
- 近傍2ホップの観測情報を利用

**観測空間** (52次元):
- 各アプローチの待ち行列長（4次元）
- 各アプローチの待ち時間（4次元）
- 現在フェーズ・経過時間（2次元）
- 近傍交差点情報（16次元）
- 時刻・需要推定値（2次元）
- バス接近フラグ（4次元）
- 予備次元（20次元）

**行動空間**: 離散4フェーズ選択（NS直進、NS左折、EW直進、EW左折）

**報酬関数** (複合報酬):
```
r = -0.4×Δ待ち時間 - 0.3×Δ待ち行列 + 0.2×スループット + 0.1×バス優先ボーナス
```

**Transit Signal Priority (TSP)**: バス接近時にグリーン延長（最大10秒）または早期グリーン（最大5秒）を適用。

**ハイパーパラメータ**:

| パラメータ | 値 |
|---|---|
| 学習率 | 3×10⁻⁴ |
| 割引率 γ | 0.99 |
| GAE λ | 0.95 |
| クリップ範囲 | 0.2 |
| エントロピー係数 | 0.01 |
| SGD反復数 | 10 |
| バッチサイズ | 4,096 |
| ワーカー数 | 8 |

### 2.4 交通需要推定（拡張カルマンフィルタ）

プローブ車両データ（浸透率15%）を活用したリアルタイムOD需要推定。

**状態空間モデル**:
- 状態ベクトル: ODフロー（25ゾーン × 25ゾーン = 625次元）
- 観測: リンク交通量（100リンク）
- 配分行列 H: 経路選択モデルから導出

**データ融合** (Multi-Source Fusion):
- ループ検知器（重み: 0.5）: 高精度・固定地点
- プローブ車両（重み: 0.3）: 広域カバレッジ
- Bluetooth/WiFiセンサ（重み: 0.2）: 旅行時間推定

![Demand Profile](figures/fig3_demand_profile.png)

### 2.5 動的リルーティング

**インシデント検出**: 速度異常検出アルゴリズム
- しきい値: 自由流速度の30%以下
- 確認時間: 120秒の持続後に確定

**リルーティングアルゴリズム**:
- A*（リアルタイム旅行時間コスト）
- Yen's K-shortest pathsによる経路多様化（K=3）
- 遵守率70%の下での車両配分

### 2.6 ネットワーク構成

![Network Topology](figures/fig4_network_topology.png)

東京都心部を8×6の格子ネットワーク（48交差点）でモデル化：
- 東西ブロック長: 500m（主要道路3車線）
- 南北ブロック長: 375m（補助道路2車線）
- 制限速度: 50 km/h
- バス路線: 12路線（各路線8停留所、5分間隔）

---

## 3. 主要な結果と数値

### 3.1 合成評価結果

3エピソード（各2時間、朝7:00-9:00ピーク含む）の合成評価を実施。

| 指標 | Episode 1 | Episode 2 | Episode 3 | 平均 ± 標準偏差 |
|---|---|---|---|---|
| 平均速度 (km/h) | 22.3 | 22.3 | 22.3 | 22.3 ± 1.0 |
| 平均遅延 (s) | 13.6 | 13.5 | 13.6 | 13.6 ± 1.7 |
| 平均待ち行列長 (台) | 4.5 | 4.5 | 4.5 | 4.5 ± 0.5 |
| スループット (台/区間) | 72.2 | 72.1 | 71.6 | 72.0 ± 10.6 |
| バス遅延 (s) | 8.9 | 8.9 | 8.9 | 8.9 ± 1.7 |
| CO₂排出量 (g/ステップ) | 460.3 | 460.4 | 460.5 | 460.4 ± 42.5 |

![Performance Metrics](figures/fig2_performance_metrics.png)

### 3.2 インシデント対応テスト

Episode 2でシミュレーション時刻2400秒にlink_3_2_E上で事故を発生させ（容量80%削減）、時刻3600秒に解消するシナリオをテスト。

- インシデント検出 → 容量を360 veh/hに低減
- A*による代替経路計算が正常動作
- K-shortest paths（K=3）による経路分散を確認
- インシデント解消後の状態復旧を確認

### 3.3 コンポーネント性能

| コンポーネント | パラメータ数 | 計算時間/ステップ |
|---|---|---|
| IDMモデル | 5パラメータ/車種 | < 1 μs |
| MAPPOネットワーク | 48エージェント × 観測52次元 | ~10 ms (推論) |
| カルマンフィルタ | 625状態 × 100観測 | ~50 ms |
| A*ルーティング | 48ノード × 82リンク | < 1 ms |

---

## 4. 考察と今後の展望

### 4.1 設計の特徴

1. **モジュラー設計**: 各コンポーネントが独立して動作・テスト可能であり、SUMOなしでも合成データによる検証が可能。

2. **スケーラビリティ**: MAPPOのパラメータ共有により、交差点数の増加に対してポリシーパラメータは一定。RLlibの分散学習（8ワーカー）により学習を高速化。

3. **マルチモーダル対応**: IDMパラメータのモード別設定、TSPによるバス優先、自転車専用レーン対応を統合。

4. **リアルタイム性**: カルマンフィルタによる5分間隔の需要更新、60秒間隔のリルーティング更新により、変動する交通状況に適応。

### 4.2 現在の制約

- **合成データ評価**: 本設計検証は合成データに基づくため、実環境での性能は別途検証が必要。
- **歩行者モデル**: Social Force Modelの完全実装は今後の課題（現在はIDMベースの近似）。
- **通信遅延**: エッジコンピューティング環境でのレイテンシは未考慮。
- **天候・イベント**: 気象条件やイベントによる需要変動は Historical Profile の拡張で対応予定。

### 4.3 今後の展望

1. **SUMO実環境統合**: TraCIインターフェースによるSUMOとの完全結合、実測データでのIDMキャリブレーション。
2. **分散学習の本格化**: Ray Clusterでの大規模MAPPO学習（5,000エピソード）、GPU活用。
3. **実データ統合**: 東京メトロオープンデータ、タクシープローブデータとの連携。
4. **Connected Vehicle対応**: V2X通信を想定したプローブ浸透率の段階的向上シナリオ。
5. **公平性評価**: 各交通モード間のサービス水準の公平性指標の導入。
6. **Transfer Learning**: 学習済みポリシーの他都市への転移学習の検証。

---

## 5. 生成したファイル一覧

### ソースコード

| ファイル | 内容 |
|---|---|
| `src/models/idm_model.py` | IDM + MOBIL 車両挙動モデル、キャリブレーション、SUMO vType XML生成 |
| `src/agents/marl_signal_control.py` | MAPPO信号制御エージェント、報酬関数、RLlib設定、東京グリッド生成 |
| `src/models/demand_estimation.py` | 拡張カルマンフィルタOD推定、マルチソース融合、時間帯需要プロファイル |
| `src/models/dynamic_routing.py` | A*ルーティング、Yen's K-shortest paths、インシデント検出・対応 |
| `src/network/sumo_environment.py` | SUMOネットワーク生成、フロー定義、メトリクス収集 |
| `src/main_orchestrator.py` | 統合パイプラインオーケストレータ |
| `src/utils/visualize.py` | 可視化スクリプト（5種類の図を生成） |

### 設定ファイル

| ファイル | 内容 |
|---|---|
| `configs/simulation_config.yaml` | 全体設定（IDMパラメータ、MARL設定、需要推定、リルーティング） |
| `requirements.txt` | Python依存パッケージ |
| `README.md` | プロジェクト概要・使用方法 |

### 結果ファイル

| ファイル | 内容 |
|---|---|
| `results/evaluation_summary.json` | 3エピソードの評価サマリー |
| `results/metrics_history.json` | 360レコードの時系列メトリクス |
| `results/mappo_config.json` | RLlib MAPPO設定 |
| `results/network/tokyo_nodes.nod.xml` | SUMOノード定義 |
| `results/network/tokyo_edges.edg.xml` | SUMOエッジ定義 |
| `results/network/vtypes.add.xml` | 車両タイプ定義 |

### 図表

| ファイル | 内容 |
|---|---|
| `figures/fig1_system_architecture.png/svg` | システムアーキテクチャ図 |
| `figures/fig2_performance_metrics.png/svg` | 交通性能指標の時系列 |
| `figures/fig3_demand_profile.png/svg` | 東京都心交通需要プロファイル |
| `figures/fig4_network_topology.png/svg` | ネットワークトポロジー |
| `figures/fig5_idm_acceleration.png/svg` | IDM加速度関数 |

### ログ

| ファイル | 内容 |
|---|---|
| `logs/process-log.jsonl` | 実行トレースログ |

---

## 参考文献

1. Treiber, M., Hennecke, A., & Helbing, D. (2000). Congested traffic states in empirical observations and microscopic simulations. *Physical Review E*, 62(2), 1805.
2. Kesting, A., Treiber, M., & Helbing, D. (2007). General lane-changing model MOBIL for car-following models. *Transportation Research Record*, 1999(1), 86-94.
3. Yu, C., et al. (2022). The surprising effectiveness of PPO in cooperative multi-agent games. *NeurIPS*.
4. Wei, H., et al. (2019). PressLight: Learning max pressure control to coordinate traffic signals in arterial network. *KDD*.
5. Cascetta, E. (2009). *Transportation Systems Analysis*. Springer.
6. Lopez, P. A., et al. (2018). Microscopic traffic simulation using SUMO. *IEEE ITSC*.
7. Wu, C., et al. (2021). Flow: A modular learning framework for mixed autonomy traffic. *IEEE T-RO*.
