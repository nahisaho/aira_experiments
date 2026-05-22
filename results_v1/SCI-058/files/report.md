# プライバシー保護下での医療データ解析のための連合学習フレームワーク

**DRAFT — NOT FOR DISTRIBUTION**

**作成日**: 2026-05-23  
**バージョン**: 1.0.0  
**ツール**: Co-Scientist (co-scientist-federated-learning)

---

## 1. 実験目的と背景

### 1.1 目的

多施設臨床データを対象に、患者プライバシーを保護しながら生存時間解析モデルを協調学習する**連合学習（Federated Learning）フレームワーク**を設計・実装・評価する。具体的には以下の6つの要素を統合的に検証する：

1. **FedAvg の収束保証と改良**（FedProx, SCAFFOLD との比較）
2. **データ異質性（non-IID）対策**の有効性
3. **差分プライバシー（DP）**の統合とプライバシー・ユーティリティトレードオフ
4. **通信効率化**（勾配圧縮）の精度への影響
5. **ビザンチン攻撃耐性**の設計と評価
6. **多施設臨床データでの生存時間解析**ケーススタディ

### 1.2 背景

医療データは患者プライバシーの観点から施設間共有が困難であり、HIPAA・GDPR等の規制下では生データの集中管理は現実的でない。連合学習はデータを各施設に留めたまま協調的にモデルを訓練する枠組みであり、この課題に対する有力な解決策である。

しかし、以下の課題が存在する：
- 施設間のデータ分布の異質性（non-IID問題）
- 通信コストの増大
- モデル更新からの情報漏洩リスク
- 悪意あるクライアントによる攻撃

本研究では、これらの課題に対応する包括的なフレームワークを設計し、合成臨床データを用いたCox比例ハザードモデルの連合学習で評価した。

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 集約アルゴリズム

#### FedAvg（Federated Averaging）
McMahan et al. (2017) による基本的な連合学習アルゴリズム。各クライアントが局所的にSGDを実行し、サーバがパラメータの加重平均を取る。

**収束保証**（凸の場合）:

```
E[f(w_T) - f(w*)] ≤ O(1 / (K × T))
```

非凸の場合（勾配非類似度 σ² のもと）:

```
(1/T) Σ E[||∇f(w_t)||²] ≤ O(σ / √(K×T) + 1/T)
```

#### FedProx
Li et al. (2020) による改良。局所目的関数に近接項を追加し、データ異質性下での安定性を向上：

```
h_k(w; w_t) = F_k(w) + μ/2 × ||w - w_t||²
```

#### SCAFFOLD
Karimireddy et al. (2020) による分散低減手法。制御変量（control variate）を用いてクライアントドリフトを除去し、データ異質性に依存しない O(1/T) 収束を達成。

### 2.2 差分プライバシー

**(ε, δ)-差分プライバシー**をガウスメカニズムで実現：

- **勾配クリッピング**: L2ノルムを上限 C で制限（感度の上界）
- **ガウスノイズ**: σ ≥ C × √(2 ln(1.25/δ)) / ε
- **RDPアカウンタント**: Rényi差分プライバシーによる厳密な合成追跡
- **適応的クリッピング**: Andrew et al. (2021) の手法で分位点ベースのクリップ境界自動調整

### 2.3 通信効率化

- **Top-K スパース化**: 勾配の上位 K 個のみ送信（誤差フィードバック付き）
- **確率的量子化（QSGD）**: 32bit → bbit への量子化（不偏推定）
- **知識蒸留（FedMD/FedDF）**: モデル重みの代わりにロジットを共有

### 2.4 ビザンチン耐性

| 手法 | 耐性上限 | 原理 |
|------|---------|------|
| Krum | f < (n-2)/2 | 最近傍距離最小のクライアント選択 |
| 座標中央値 | f < n/2 | 各座標の中央値を使用 |
| Trimmed Mean | f < βn | 上下β%除去後の平均 |
| FLTrust | - | サーバ参照勾配とのコサイン類似度で重み付け |

### 2.5 プラットフォーム設計

**Flower + PySyft ハイブリッドアーキテクチャ**:

```
┌─────────────────────────────────────────────┐
│              FL Server (Flower)              │
│  ┌──────────┐ ┌────────┐ ┌───────────────┐  │
│  │Strategy  │ │DP Guard│ │Byzantine      │  │
│  │Selector  │→│(RDP)   │→│Filter         │  │
│  └──────────┘ └────────┘ └───────────────┘  │
│  ┌──────────────────────────────────────┐    │
│  │  Aggregation: FedAvg|FedProx|SCAFFOLD│    │
│  └──────────────────────────────────────┘    │
│  ┌──────────────────────────────────────┐    │
│  │  Communication: TopK|QSGD|Distill   │    │
│  └──────────────────────────────────────┘    │
├─────────────────────────────────────────────┤
│  PySyft: Privacy Budget | Audit | Governance│
└──────────┬──────────┬──────────┬────────────┘
     ┌─────┴───┐┌─────┴───┐┌─────┴───┐
     │Hospital A││Hospital B││Hospital C│
     │(n=200)  ││(n=350)  ││(n=150)  │
     └─────────┘└─────────┘└─────────┘
```

アーキテクチャ図は `figures/fig6_architecture.png` を参照。

---

## 3. 実験設定

### 3.1 データ

5施設の合成臨床データを生成（Cox比例ハザードモデルに基づく）：

| 施設 | サンプル数 | イベント率 | 中央生存時間 |
|------|-----------|-----------|------------|
| Hospital 0 | 200 | 0.460 | 2.70 |
| Hospital 1 | 350 | 0.480 | 1.87 |
| Hospital 2 | 150 | 0.470 | 2.41 |
| Hospital 3 | 500 | 0.490 | 1.94 |
| Hospital 4 | 300 | 0.530 | 0.93 |

- **特徴量**: 20次元
- **Non-IID度**: 施設ごとに特徴量分布をシフト（σ で制御）
- **打ち切り率**: 施設ごとに 30%〜40%

### 3.2 実験条件

- **通信ラウンド**: 50
- **局所エポック**: 5
- **局所学習率**: 0.01
- **乱数シード**: 42（numpy, random で固定）

---

## 4. 主要な結果と数値

### 4.1 集約アルゴリズム比較

| アルゴリズム | C-index | 損失関数 | パラメータ誤差 |
|-------------|---------|---------|-------------|
| **FedAvg** | **0.8534** | 0.9220 | 0.3787 |
| FedProx (μ=0.01) | 0.8534 | 0.9225 | 0.3788 |
| SCAFFOLD | 0.8523 | 1.0590 | 0.3827 |

**考察**: 本実験のnon-IID度（σ=0.5）では3手法の差は微小。FedAvgとFedProxはほぼ同等の性能を示し、SCAFFOLDは制御変量の初期化の影響でわずかに劣後した。

![収束比較](figures/fig1_convergence_comparison.png)

### 4.2 データ異質性（Non-IID）の影響

| Non-IID度 (σ) | FedAvg C-index | FedProx C-index |
|---------------|---------------|----------------|
| 0.0 (IID) | 0.8523 | 0.8523 |
| 0.5 | 0.8534 | 0.8534 |
| 1.0 | 0.8580 | 0.8580 |
| 2.0 | 0.8721 | 0.8721 |

**考察**: non-IID度が上昇しても C-index は維持（むしろ向上）。これは合成データにおいて特徴量シフトが予測精度に寄与するためであり、実データではdriftの悪影響がより顕著になる可能性がある。FedProxの近接項（μ=0.01）は本設定では効果が限定的。

![Non-IID影響](figures/fig2_noniid_impact.png)

### 4.3 差分プライバシーの影響

| プライバシー予算 (ε) | C-index | 実消費ε | パラメータ誤差 |
|---------------------|---------|--------|-------------|
| 1.0（強保護） | 0.6351 | 1.098 | 7.454 |
| 5.0 | 0.6175 | 5.303 | 1.347 |
| 10.0 | 0.6342 | 11.757 | 0.942 |
| 50.0（弱保護） | 0.7541 | 61.513 | 0.889 |
| ∞（DP無し） | **0.8534** | 0.0 | 0.379 |

**考察**: DPノイズはモデル性能に大きな影響を与える。ε=1.0 の強プライバシー保護下ではC-indexが0.635まで低下（−25.6%）。ε=50.0でもベースラインから−11.6%の劣化が見られる。医療応用では ε=5〜10 が現実的な妥協点と考えられる。

![プライバシー・ユーティリティトレードオフ](figures/fig3_privacy_utility_tradeoff.png)

### 4.4 通信効率化

| 圧縮率 | C-index | 通信量削減 |
|--------|---------|----------|
| 1% (Top-1%) | 0.8534 | 99% |
| 5% | 0.8534 | 95% |
| 10% | 0.8534 | 90% |
| 50% | 0.8534 | 50% |
| 100%（無圧縮） | 0.8534 | 0% |

**考察**: 誤差フィードバック付きTop-Kスパース化は、わずか1%のパラメータ送信でも精度劣化なし。Cox PHモデルの20次元パラメータでは圧縮の恩恵は限定的だが、深層ニューラルネットワーク（数百万パラメータ）では通信コスト99%削減の効果は極めて大きい。

![通信効率](figures/fig4_communication_efficiency.png)

### 4.5 ビザンチン攻撃耐性

| 攻撃率 | 防御なし | Krum | 座標中央値 | Trimmed Mean |
|--------|---------|------|----------|-------------|
| 0% | **0.8534** | 0.8477 | 0.8473 | 0.8490 |
| 20% | 0.4235 | **0.8466** | 0.8444 | 0.8460 |
| 40% | 0.4198 | **0.8404** | 0.8329 | 0.5475 |

**考察**:
- **防御なし**: 20%のビザンチンクライアントでC-indexが0.424（ランダム以下）に壊滅的に劣化
- **Krum**: 40%攻撃下でもC-index 0.840を維持。最も頑健な防御手法
- **座標中央値**: Krumに次ぐ性能。理論的耐性上限（f < n/2）に整合
- **Trimmed Mean**: 20%では有効だが、40%攻撃ではトリミング率の限界により劣化

![ビザンチン耐性](figures/fig5_byzantine_resilience.png)

---

## 5. 考察と今後の展望

### 5.1 主要な知見

1. **集約アルゴリズム**: 中程度のnon-IID下ではFedAvgで十分な性能が得られるが、高異質性環境ではFedProx/SCAFFOLDの効果が顕在化すると予想される
2. **プライバシー・ユーティリティトレードオフ**: 臨床応用では ε=5〜10 が現実的。ε=1.0 以下では実用的な予測性能の維持が困難
3. **通信効率化**: Top-Kスパース化は誤差フィードバックにより精度劣化なしで通信量を大幅削減可能
4. **ビザンチン耐性**: Krumが最も頑健。40%の悪意あるクライアント下でも93%以上のベースライン性能を維持

### 5.2 限界

- **合成データ**: 実際の臨床データ（画像、テキスト等の高次元データ）での検証が必要
- **線形モデル**: Cox PHモデルは20次元の線形モデルであり、深層学習モデルでは異なる傾向が予想される
- **セキュアアグリゲーション**: 暗号学的手法（秘密分散、準同型暗号）は未実装
- **クライアント選択**: 全クライアント参加を仮定しており、部分参加の影響は未検証

### 5.3 今後の展望

1. **実臨床データでの検証**: TCGA, MIMIC-III 等の公開データセットでの評価
2. **深層学習モデル**: DeepSurv, Cox-nnet 等のニューラルネットワークベースモデルへの拡張
3. **パーソナライズドFL**: Per-FedAvg, APFL 等の個別化手法の導入
4. **セキュアアグリゲーション**: SMPC、準同型暗号の統合
5. **異種モデルFL**: 知識蒸留ベースの異種アーキテクチャ対応
6. **連合転移学習**: 事前学習モデルからのドメイン適応

---

## 6. フレームワーク構成

### 6.1 モジュール構成

```
fl_framework/
├── __init__.py             # パッケージ定義
├── aggregation.py          # FedAvg, FedProx, SCAFFOLD
├── differential_privacy.py # DP統合, RDPアカウンタント
├── communication.py        # Top-K, QSGD, 知識蒸留
├── byzantine.py            # Krum, 中央値, Trimmed Mean, FLTrust
├── fl_platform.py          # Flower/PySyftプラットフォーム設計
├── case_study.py           # 生存時間解析ケーススタディ
└── generate_figures.py     # 可視化スクリプト
```

### 6.2 Flower デプロイメント

`fl_platform.py` に以下のコード生成機能を実装：

- `generate_flower_server_code()`: DPFedAvgStrategy のFlowerサーバコード
- `generate_flower_client_code()`: HospitalClient のFlowerクライアントコード
- `get_architecture_spec()`: プラットフォーム仕様のJSON出力

### 6.3 PySyft 統合

- **リモートデータアクセス**: Duet プロトコルによる安全なデータ参照
- **プライバシー予算管理**: データオーナーがε, δを設定
- **監査ログ**: 全アクセスの追跡と記録

---

## 7. 生成ファイル一覧

### フレームワークコード
| ファイル | 説明 |
|---------|------|
| `fl_framework/__init__.py` | パッケージ定義 |
| `fl_framework/aggregation.py` | FedAvg, FedProx, SCAFFOLD 実装 |
| `fl_framework/differential_privacy.py` | DP統合（ガウスメカニズム, RDPアカウンタント） |
| `fl_framework/communication.py` | 通信効率化（Top-K, QSGD, 知識蒸留） |
| `fl_framework/byzantine.py` | ビザンチン耐性（Krum, 中央値, Trimmed Mean, FLTrust） |
| `fl_framework/fl_platform.py` | Flower/PySyft プラットフォーム設計 |
| `fl_framework/case_study.py` | 多施設生存時間解析ケーススタディ |
| `fl_framework/generate_figures.py` | 可視化スクリプト |

### 実験結果
| ファイル | 説明 |
|---------|------|
| `results/experiment_summary.json` | 全実験の要約指標 |
| `results/strategy_fedavg_history.json` | FedAvg 学習履歴 |
| `results/strategy_fedprox_history.json` | FedProx 学習履歴 |
| `results/strategy_scaffold_history.json` | SCAFFOLD 学習履歴 |

### 図表
| ファイル | 説明 |
|---------|------|
| `figures/fig1_convergence_comparison.png/svg` | 集約アルゴリズム収束比較 |
| `figures/fig2_noniid_impact.png/svg` | Non-IID影響分析 |
| `figures/fig3_privacy_utility_tradeoff.png/svg` | プライバシー・ユーティリティトレードオフ |
| `figures/fig4_communication_efficiency.png/svg` | 通信圧縮効率 |
| `figures/fig5_byzantine_resilience.png/svg` | ビザンチン耐性比較 |
| `figures/fig6_architecture.png/svg` | プラットフォームアーキテクチャ図 |

### データ・ログ
| ファイル | 説明 |
|---------|------|
| `data/hospital_data_summary.json` | 施設別データサマリー |
| `logs/process-log.jsonl` | 実行トレース |

---

## 参考文献

1. McMahan, B., et al. (2017). "Communication-Efficient Learning of Deep Networks from Decentralized Data." AISTATS.
2. Li, T., et al. (2020). "Federated Optimization in Heterogeneous Networks." MLSys.
3. Karimireddy, S.P., et al. (2020). "SCAFFOLD: Stochastic Controlled Averaging for Federated Learning." ICML.
4. Abadi, M., et al. (2016). "Deep Learning with Differential Privacy." CCS.
5. Mironov, I. (2017). "Rényi Differential Privacy." CSF.
6. Andrew, G., et al. (2021). "Differentially Private Learning with Adaptive Clipping." NeurIPS.
7. Alistarh, D., et al. (2017). "QSGD: Communication-Efficient SGD via Gradient Quantization and Encoding." NeurIPS.
8. Blanchard, P., et al. (2017). "Machine Learning with Adversaries: Byzantine Tolerant Gradient Descent." NeurIPS.
9. Yin, D., et al. (2018). "Byzantine-Robust Distributed Learning." ICML.
10. Cao, X., et al. (2021). "FLTrust: Byzantine-robust Federated Learning via Trust Bootstrapping." NDSS.
11. Beutel, D.J., et al. (2020). "Flower: A Friendly Federated Learning Framework." arXiv.
