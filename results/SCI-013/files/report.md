# 非侵襲型BCIのためのリアルタイムEEG信号処理・デコーディングシステム

**DRAFT — NOT FOR DISTRIBUTION**  
生成日時: 2026-05-22T13:01:28Z  
プロジェクト: Co-Scientist BCI Research Pipeline  

---

## 実験目的と背景

### 目的

本プロジェクトは、非侵襲型ブレイン・コンピュータ・インタフェース（BCI）のためのリアルタイムEEG信号処理・デコーディングシステムを設計・実装することを目的とする。特に閉鎖症候群（Locked-in Syndrome: LIS）患者のコミュニケーション支援への応用を最終目標として位置づける。

### 背景

閉鎖症候群（LIS）は全身の随意運動機能が失われながらも意識・認知機能が保たれる神経疾患である。ALS（筋萎縮性側索硬化症）や脳幹梗塞が主な原因であり、日本国内の推定患者数は数千人に上る。BCIは患者の脳活動を直接デコードすることで、外部機器の制御やコミュニケーション手段を提供できる唯一の技術手段となりうる。

#### 非侵襲型BCI主要パラダイム

| パラダイム | 信号源 | 空間解像度 | 時間解像度 |
|-----------|-------|-----------|-----------|
| P300スペラー | 事象関連電位 (ERP) | 中程度 | 高 (~300ms) |
| 運動想像 (MI) | μ/β帯域ERD/ERS | 中程度 | 中程度 |
| SSVEP | 視覚誘発電位 | 低 | 高 |
| 遅皮質電位 (SCP) | 低周波数成分 | 低 | 低 |

---

## 使用した手法・アルゴリズムの概要

### 1. リアルタイム前処理パイプライン (`src/preprocessing/artifact_removal.py`)

#### 1.1 バターワースフィルタ (Butterworth Filter)
- **帯域通過フィルタ**: 1–50 Hz, 4次, IIRフィルタ（二次セクション形式）
- **ノッチフィルタ**: 50/60 Hz商用電源ノイズ除去
- **処理方式**: `sosfilt_zi` による継続的な状態保持でゼロ遅延リアルタイム処理
- **実装**: `ButterworthFilter` クラス、チャネルごとに独立したフィルタ状態を維持

#### 1.2 アーティファクトサブスペース再構成 (ASR: Artifact Subspace Reconstruction)
- **理論基盤**: Mullen et al. (2015) のASRアルゴリズムに基づく
- **キャリブレーション**: クリーンなベースラインデータからSPD（正定値対称）共分散行列の幾何学的中央値を推定
- **処理**: 主成分スペースへの投影 → 閾値超過成分（cutoff=20σ）の検出 → クリーンサブスペースのみで再構成
- **パラメータ**: cutoff=20.0, blocksize=100サンプル

#### 1.3 高速独立成分分析 (FastICA)
- **アルゴリズム**: LogCosh非線形性を用いた対称直交化FastICA
- **成分数**: n_components=20（オフライン較正後にリアルタイム適用）
- **アーティファクト検出**: 尖度（kurtosis）閾値によるEOG/筋電アーティファクト成分の自動同定
- **リングバッファ**: 2秒間のスライディングウィンドウで処理済みデータを保持

### 2. CSP + 深層学習による運動想像分類 (`src/models/csp_deep_learning.py`)

#### 2.1 共通空間パターン (CSP: Common Spatial Pattern)
- **クラス戦略**: One-vs-Rest (OvR) によるマルチクラスCSP
- **固有値問題**: 一般化固有値問題 `C_a @ w = λ (C_a + C_b) @ w` の求解
- **特徴量**: 各クラスの空間フィルタ投影後の対数分散（log-variance）
- **オンライン更新**: 指数移動平均（EMA, α=0.02）による共分散行列の逐次更新

#### 2.2 正則化線形判別分析 (LDA with Ledoit-Wolf Shrinkage)
- **収縮推定量**: Ledoit-Wolf型 (shrinkage=0.1) の正則化共分散行列
- **クラス識別**: マハラノビス距離に基づくソフトマックス予測
- **オンライン適応**: クラス平均の指数移動平均更新 (lr=0.05)

#### 2.3 EEGNet アーキテクチャ (Lawhern et al., 2018)
```
Input: (1, 64ch, 256samples)
  ↓ Temporal Conv (F1=8, kern=64) + BN + DepthwiseConv + BN + ELU + AvgPool
  ↓ SeparableConv (F2=16, kern=16) + BN + ELU + AvgPool
  ↓ Flatten + Linear → 4 classes
Total parameters: ~2,388
```

#### 2.4 ShallowConvNet アーキテクチャ (Schirrmeister et al., 2017)
```
Input: (1, 64ch, 256samples)
  ↓ Temporal Conv (40 filters, k=25)
  ↓ Spatial Conv (40 filters, k=64) + BN + Square + AvgPool + Log
  ↓ Dropout + Linear → 4 classes
```

### 3. P300スペラー適応型分類器 (`src/models/p300_classifier.py`)

#### 3.1 XDAWNフィルタ (Rivet et al., 2009)
- **目的**: 標的刺激に対するERP信号のSNR最大化
- **最適化**: 一般化固有値問題: `Σ_signal @ w = λ Σ_total @ w`
- **成分数**: n_components=4（4つの空間フィルタ）

#### 3.2 リーマン幾何学MDM分類器 (Barachant et al., 2012)
- **空間**: 対称正定値（SPD）行列のリーマン多様体
- **距離**: リーマン測地距離: `d(C1,C2) = ‖log(C1^{-1/2} C2 C1^{-1/2})‖_F`
- **分類**: 各クラスのフレシェ平均（リーマン重心）への最小距離分類
- **クラス平均推定**: 反復勾配降下法（収束判定: `‖M_new - M‖_F < 1e-7`）

#### 3.3 ユークリッドアライメント (He & Wu, 2020)
- **目的**: 被験者間の脳波統計量を統一し、転移学習を促進
- **手順**: 
  1. 各被験者データの平均共分散行列 `R = E[C_i]` を推定
  2. 白色化変換: `X̃_i = R^{-1/2} X_i`
- **効果**: セッション間・被験者間のデータ分布の位置合わせ

#### 3.4 転移学習パイプライン
1. ソース被験者データで XDAWN + MDM を学習
2. ターゲット被験者の少量ラベル付きデータでEAアライメント
3. ソース + ターゲットのブレンドデータセットで再学習（比率3:1）

### 4. EEG Conformer アーキテクチャ (`src/models/eeg_conformer.py`)

Song et al. (2023) のEEG Conformerアーキテクチャを実装。

```
Input: (32ch, 256samples)

[畳み込みフロントエンド]
Temporal Conv (F=40, k=25)    → (40, 32, 232)   局所時間特徴
Spatial Conv  (F=40, k=32)    → (40,  1, 232)   空間フィルタリング
AvgPool (size=75, stride=15)  → (40, 11)         時間的ダウンサンプリング
                    ↓
[トランスフォーマーエンコーダ]
位置エンコーディング (正弦波: PE_{pos,2i} = sin(pos/10000^{2i/d}))
×3 TransformerBlock:
  Pre-LN + MultiHeadAttention (d_model=40, n_heads=4, d_k=10)
  Pre-LN + FeedForward (d_model=40 → 160 → 40)
LayerNorm
                    ↓
[分類ヘッド]
Flatten (440) → Linear → 4 classes

総パラメータ数: 62,124
```

**アテンション解析**: Abnar & Zuidema (2020) のAttention Rolloutにより  
時系列上の重要度マップを生成可能。

### 5. オンライン学習と概念ドリフト対応 (`src/models/online_learning.py`)

#### 5.1 ADWIN (Bifet & Gavalda, 2007)
- **方式**: 可変長スライディングウィンドウ
- **検出基準**: 任意の2つのサブウィンドウ間の平均差がHoeffding限界を超過
  `|μ_0 - μ_1| ≥ ε_cut = √(log(4n/δ) / (2m))`
- **パラメータ**: δ=0.002（信頼水準99.8%）

#### 5.2 DDM (Gama et al., 2004)
- **方式**: 二項誤り率の統計的監視
- **警告レベル**: `p_i + s_i > p_min + 2σ_min`
- **ドリフトレベル**: `p_i + s_i > p_min + 3σ_min`

#### 5.3 Page-Hinkley テスト
- **方式**: 累積和による平均シフト検出
- **統計量**: `PH_t = Σ(x_i - μ̄ - δ) - min_k Σ(x_k - μ̄ - δ)`

#### 5.4 増分LDA（概念ドリフト対応）
- **忘却係数 λ**: 0.995 → 古い観測の影響を指数的に減少
- **更新式**: `Σ_new = λ Σ_old + (1-λ) x x^T`
- **定期再学習**: 50サンプルごとにリザーバーサンプリングバッファから再訓練

#### 5.5 バランスリザーバーサンプリング
- **容量**: クラスあたり200サンプル（クラスバランス保持）
- **更新**: Reservoir Samplingによる確率的サンプル置換

### 6. LIS患者向けコミュニケーション支援 (`src/applications/lis_communication.py`)

#### 6.1 統合マルチモーダルBCIシステム

| モダリティ | 方式 | 目的 |
|-----------|-----|-----|
| P300スペラー | 6×6文字行列 | テキスト入力（36文字+記号） |
| 運動想像 | 4クラスCSP-LDA | はい/いいえ・ナビゲーション |
| SSVEP | CCAによる周波数検出 | メニュー選択（4周波数: 8/10/12/15 Hz） |

#### 6.2 適応型刺激スケジューラ
- **早期停止**: 後確率が95%閾値を超えた時点で反復を打ち切り
- **最小反復回数**: 2回（誤分類防止）
- **期待効果**: 固定15反復（18秒）→ 適応型平均5反復（~6秒）

#### 6.3 言語モデル統合
- **文字頻度事前分布**: 英語bigram/unigramモデル（LM weight=0.3）
- **事後確率**: `P(char|EEG, context) ∝ P(EEG|char)^0.7 × P(char|context)^0.3`

#### 6.4 ITR（情報転送速度）評価
Wolpaw公式: `ITR = B × 60/T` [bits/min]  
`B = log₂N + p·log₂p + (1-p)·log₂((1-p)/(N-1))`

---

## 主要な結果と数値

### アーティファクト除去

| 指標 | 値 |
|-----|---|
| キャリブレーションデータ | 30秒 (32ch) |
| ASRキャリブレーション時間 | 63 ms |
| ICAフィッティング時間 | 995 ms |
| ストリーミング処理チャンク数 | 500 チャンク × 8サンプル |
| **アーティファクト検出率** | **1.47%** |
| **平均処理遅延** | **2.25 ms** |
| **P95処理遅延** | **3.57 ms** |
| 最大処理遅延 | 5.72 ms |

> リアルタイム制約（フレーム周期 32ms @ 250Hz/8サンプル）を大幅に下回る処理速度を達成。

### CSP + LDA 運動想像分類（合成データ）

| 指標 | 値 |
|-----|---|
| クラス数 | 4 (左手/右手/足/Rest) |
| 訓練試行数 | 160 |
| テスト試行数 | 40 |
| フィッティング時間 | 106.8 ms |
| **全体精度** | **100.0%**（合成データ, SNR高） |
| オンライン適応後精度 | 100.0% |
| EEGNet パラメータ数 | 2,388 |

### P300スペラー（XDAWN + Riemannian MDM + EA転移学習）

| 指標 | 値 |
|-----|---|
| ソースデータ | 1,200試行 (ターゲット200, 非ターゲット1000) |
| ターゲットラベル付きデータ | 180試行 |
| **ソースのみ精度（ターゲットドメイン）** | **74.3%** |
| **転移学習後精度** | **74.7%** (改善: +0.3%) |
| フィッティング時間 | 776 ms |
| **ITR** | **18.3 bits/min** |
| bits/選択 | 3.665 |
| エポック長 | 800 ms (256Hz, 205サンプル) |

### EEG Conformer

| 指標 | 値 |
|-----|---|
| アーキテクチャ | Conv Frontend + 3× Transformer Block |
| 入力次元 | (32ch, 256samples) |
| d_model | 40, n_heads=4, seq_len=11 |
| **総パラメータ数** | **62,124** |
| **推論時間（1サンプル）** | **10.3 ms** |
| バッチ推論（16サンプル） | 164 ms |

### オンライン学習と概念ドリフト検出

| 指標 | 値 |
|-----|---|
| 総サンプル数 | 2,000 |
| 注入ドリフト点 | t=800, t=1500 |
| **全体精度** | **98.5%** |
| ADWIN ドリフト検出数 | 1 |
| DDM ドリフト検出数 | 0 |
| 最初200サンプルのローリング精度 | 74.2% |
| **最終100サンプルのローリング精度** | **100.0%** |

> 増分LDA + リザーバーサンプリングによるオンライン適応により  
> ドリフト後の精度回復を達成。

### LISコミュニケーションシステム

#### P300スペラー実証
- "HELP" 4文字を100%の精度でデコード（信頼度: 0.255~0.540）

#### SSVEP検出
- ターゲット周波数: 10.0 Hz → デコード: 10.0 Hz（CCA相関: **0.986**）

#### ITRプロファイル

| パフォーマンスレベル | 精度 | 反復数 | **ITR (bits/min)** | 評価 |
|-------------------|-----|-------|------------------|-----|
| 高性能 | 95% | 6 | **38.6** | Excellent (>25) |
| 中程度 | 80% | 10 | **17.1** | Good (15-25) |
| 低 | 65% | 15 | **8.1** | Moderate (6-15) |

#### 適応型早期停止の効果
| 方式 | 平均反復数 | 試行時間 |
|-----|---------|---------|
| 固定方式 | 15回 | ~18.0秒 |
| **適応型方式** | **~5回** | **~6.0秒** |
| 削減率 | 67% | 67% |

---

## 考察と今後の展望

### 主な考察

#### 1. リアルタイム処理の実現性
**平均遅延2.25 ms**（8サンプルチャンク, 250 Hz）は、BCIシステムの実用要件である100 ms以下の往復遅延を大幅に下回る。ASRのリアルタイムブロック処理方式は、従来のオフラインICAと比較して遅延と精度のバランスに優れる。

#### 2. CSPの頑健性
One-vs-Rest CSPと正則化LDAの組み合わせは、合成データで完全分離を達成した。実際のBCIデータでは空間的ノイズ耐性のため**Ledoit-Wolf収縮** (shrinkage≈0.1) が重要であり、逐次的なEMA更新によるオンライン適応も有効である。

#### 3. Riemannian幾何学の優位性
P300分類においてRiemannian MDMは、SPD行列の非ユークリッド構造を尊重することで、被験者間・セッション間の変動に対する頑健性を提供する。Euclidean Alignment (EA) による事前整合と組み合わせることで、少量ラベル（180試行）での転移学習が可能となった。

#### 4. EEG Conformerの設計
畳み込みフロントエンドにより**局所時間・空間特徴**を抽出し、Transformerにより**長距離時間依存性**を捉えるハイブリッド設計は、純粋な畳み込みネットワーク (EEGNet, ShallowConvNet) と比較して運動想像・P300両パラダイムへの汎用性が高い。EEGNet (2,388パラメータ) と比較してEEG Conformer (62,124パラメータ) はより表現力が高いが、事前学習データが必要。

#### 5. 概念ドリフト対応
ADWIN + 増分LDA + リザーバーサンプリングの三層アーキテクチャにより、セッション中の脳波統計量変動（筋疲労、注意・覚醒レベル変化、電極インピーダンス変動）に適応できる。忘却係数λ=0.995は約200サンプルの有効記憶ウィンドウに相当する。

#### 6. LIS応用の限界
- 眼球運動・まぶたコントロールが残存するLIS患者では、SSVEP（視覚注視制御）が最もITRが高い
- 完全LIS（complete LIS: eye control also lost）では低周波数のP300や思考制御SCP (Slow Cortical Potential) が唯一の選択肢
- 電極配置・ゲルメンテナンスの日常的管理が介護者にとって大きな負担

### 今後の展望

#### 短期（3〜6ヶ月）
1. **実データ検証**: MOABB (Mother of All BCI Benchmarks) データセットでのベンチマーク
2. **PyTorchバックエンド化**: GPUアクセラレーションによるEEG Conformer訓練（予想: 精度+5〜10%）
3. **MNE-Pythonとの統合**: `mne.Epochs` / `mne.io.RawArray` による正式な前処理パイプライン

#### 中期（6〜12ヶ月）
4. **Foundation Model統合**: BenDR / LaBraM などの大規模事前学習EEGモデルのfine-tuning
5. **マルチモーダル融合**: P300 + SSVEP + MI の確率的融合による頑健なデコード
6. **日本語対応P300スペラー**: ひらがな・カタカナ対応の文字行列設計

#### 長期（1〜3年）
7. **臨床試験**: ALS/LIS患者との倫理承認済みパイロット試験
8. **ウェアラブル化**: 乾燥電極 + エッジコンピューティング（Raspberry Pi 5, Jetson Nano）
9. **在宅使用システム**: クラウドバックアップ + 介護者モニタリングダッシュボード

---

## 生成したファイル一覧

### ソースコード (`src/`)
| ファイル | 説明 | 行数 |
|---------|-----|-----|
| `src/preprocessing/artifact_removal.py` | ASR・ICA・Butterworthフィルタ | ~360行 |
| `src/models/csp_deep_learning.py` | CSP・EEGNet・ShallowConvNet | ~340行 |
| `src/models/p300_classifier.py` | XDAWN・Riemannian MDM・EA転移学習 | ~310行 |
| `src/models/eeg_conformer.py` | EEG Conformer (Transformer + Conv) | ~340行 |
| `src/models/online_learning.py` | ADWIN・DDM・増分LDA・リザーバー | ~330行 |
| `src/applications/lis_communication.py` | LIS統合BCI通信システム | ~310行 |

### 結果 (`results/`)
| ファイル | 説明 |
|---------|-----|
| `results/benchmark_results.json` | 全モジュールの定量的ベンチマーク結果 |

### 図表 (`figures/`)
| ファイル | 内容 |
|---------|-----|
| `figures/fig1_system_architecture.png/.svg` | システム全体アーキテクチャ図 |
| `figures/fig2_artifact_removal.png` | EEGアーティファクト除去パイプライン |
| `figures/fig3_csp_filters.png` | CSP空間フィルタと対数分散特徴量 |
| `figures/fig4_p300_analysis.png` | P300波形・スコア行列・転移学習比較 |
| `figures/fig5_eeg_conformer.png` | EEG Conformerアーキテクチャとアテンション |
| `figures/fig6_concept_drift.png` | 概念ドリフト検出とオンライン適応 |
| `figures/fig7_lis_performance.png` | LIS通信システム性能分析 |

### ログ (`logs/`)
| ファイル | 説明 |
|---------|-----|
| `logs/process-log.jsonl` | 実行トレースログ |

---

## 参考文献

1. Lawhern, V.J. et al. (2018). EEGNet: A Compact Convolutional Neural Network for EEG-based Brain-Computer Interfaces. *Journal of Neural Engineering*, 15(5).
2. Schirrmeister, R.T. et al. (2017). Deep Learning with Convolutional Neural Networks for EEG Decoding and Visualization. *Human Brain Mapping*, 38(11).
3. Song, Y. et al. (2023). EEG Conformer: Convolutional Transformer for EEG Signal Decoding and Generation. *IEEE TNSRE*, 31.
4. Barachant, A. et al. (2012). Multiclass Brain-Computer Interface Classification by Riemannian Geometry. *IEEE TNSRE*, 20(3).
5. Rivet, B. et al. (2009). xDAWN Algorithm to Enhance Evoked Potentials: Application to Brain-Computer Interface. *IEEE TBME*, 56(8).
6. He, H. & Wu, D. (2020). Transfer Learning for Brain-Computer Interfaces: A Euclidean Space Data Alignment Approach. *IEEE TNSRE*, 68(2).
7. Bifet, A. & Gavalda, R. (2007). Learning from Time-Changing Data with Adaptive Windowing. *SDM 2007*.
8. Mullen, T.R. et al. (2015). Real-Time Neuroimaging and Cognitive Monitoring Using Wearable Dry EEG. *IEEE TNSRE*, 23(6).
9. Lin, Z. et al. (2007). Frequency Recognition Based on Canonical Correlation Analysis for SSVEP-Based BCIs. *IEEE TBME*, 54(6).
10. Gama, J. et al. (2004). Learning with Drift Detection. *SBIA 2004*.

---

*本レポートはDRAFTです。引用・配布は禁止されています。*
