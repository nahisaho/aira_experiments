# NeuroSense mHealth Framework: 実験レポート

## スマートフォンセンサーデータからの神経変性疾患早期バイオマーカー検出

---

## 1. 実験目的と背景

本研究では、スマートフォンに搭載された各種センサーから取得可能なデータを活用し、神経変性疾患（パーキンソン病、ALS、認知機能低下）の早期バイオマーカーを検出するモバイルヘルス（mHealth）データ解析フレームワーク「**NeuroSense**」を設計・実装した。

### 研究の動機

- 神経変性疾患は早期発見が治療予後を大きく左右する
- 従来の臨床評価は病院受診が必要であり、頻度と客観性に限界がある
- スマートフォンは加速度センサー、ジャイロスコープ、マイク、タッチスクリーンなど多様なセンサーを搭載しており、日常的なデジタルバイオマーカー収集が可能

### 対象疾患と検出モダリティ

| モダリティ | 対象疾患 | センサー |
|---|---|---|
| 歩行パターン | パーキンソン病 | 加速度計・ジャイロスコープ |
| 音声特徴量 | ALS（筋萎縮性側索硬化症） | マイク |
| タッチ操作パターン | 認知機能低下 | タッチスクリーン |

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 システムアーキテクチャ

![NeuroSense System Architecture](figures/system_architecture.png)

### 2.2 特徴量抽出

**歩行特徴量（18特徴量）：**
- 加速度統計量：平均、標準偏差、範囲、RMS、歪度、尖度
- 歩行規則性：ステップ規則性、ストライド規則性（自己相関ベース）
- 横方向非対称性
- ジャイロスコープ統計量：平均、標準偏差、範囲、RMS
- 周波数領域：支配的周波数、スペクトルエントロピー
- ジャーク（加速度の微分）：平均、標準偏差

**音声特徴量（17特徴量）：**
- 基本周波数（F0）
- ジッター（周波数擾乱）
- シマー（振幅擾乱）
- 調波対雑音比（HNR）
- MFCC 13係数

**タッチ特徴量（8特徴量）：**
- 反応時間、タップ精度、スワイプ速度
- ダブルタップ間隔変動、タイピング速度
- エラー率、圧力変動、トレイル追跡時間

### 2.3 分類アルゴリズム

5種類の機械学習モデルを比較評価：
1. **ロジスティック回帰（LR）**
2. **ランダムフォレスト（RF）** — 100本のデシジョンツリー
3. **勾配ブースティング（GB）** — 100ステージ
4. **SVM（RBFカーネル）**
5. **多層パーセプトロン（MLP）** — 隠れ層 [64, 32]

### 2.4 マルチモーダル融合

4つの融合戦略を比較：
- **平均融合**：各モダリティの予測確率の算術平均
- **重み付き融合**：歩行(0.45), 音声(0.30), タッチ(0.25)
- **メタ学習器（LR）**：ロジスティック回帰によるスタッキング
- **メタ学習器（GB）**：勾配ブースティングによるスタッキング

### 2.5 変化点検出アルゴリズム

3つの変化点検出手法を実装・比較：
- **CUSUM**：累積和法
- **PELT**：Pruned Exact Linear Time法
- **ベイジアンオンラインCPD**：Adams & MacKay (2007) のアルゴリズム

---

## 3. 主要な結果

### 3.1 歩行ベースパーキンソン病スクリーニング

![Gait Model Comparison and ROC Curves](figures/gait_model_comparison.png)

| モデル | Accuracy | Precision | Recall | F1 | AUC-ROC |
|---|---|---|---|---|---|
| Logistic Regression | 0.995 | 0.990 | 1.000 | 0.995 | 1.000 |
| Random Forest | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Gradient Boosting | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| SVM (RBF) | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| MLP | 0.995 | 1.000 | 0.990 | 0.995 | 1.000 |

**特徴量重要度：**

![Gait Feature Importance](figures/gait_feature_importance.png)

**混同行列（最良モデル）：**

![Confusion Matrix](figures/gait_confusion_matrix.png)

### 3.2 音声ベースALS進行モニタリング

![Voice Analysis Results](figures/voice_analysis.png)

| モデル | Accuracy | Precision | Recall | F1 | AUC-ROC |
|---|---|---|---|---|---|
| Logistic Regression | 0.963 | 0.953 | 0.972 | 0.962 | 0.993 |
| Random Forest | 0.965 | 0.965 | 0.965 | 0.965 | 0.990 |
| Gradient Boosting | 0.963 | 0.960 | 0.967 | 0.963 | 0.991 |
| SVM (RBF) | 0.963 | 0.956 | 0.969 | 0.962 | 0.992 |
| MLP | 0.957 | 0.953 | 0.961 | 0.957 | 0.991 |

**経時的音声特徴量変化：**

![Voice Feature Progression](figures/voice_progression.png)

ALS群ではセッション経過に伴い、ジッターとシマーが上昇し、基本周波数（F0）が低下する傾向が明確に確認された。

### 3.3 タッチスクリーン認知機能低下検出

![Touch Analysis Results](figures/touch_analysis.png)

| モデル | Accuracy | Precision | Recall | F1 | AUC-ROC |
|---|---|---|---|---|---|
| Logistic Regression | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Random Forest | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Gradient Boosting | 0.996 | 0.993 | 1.000 | 0.996 | 0.998 |
| SVM (RBF) | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| MLP | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

![Touch ROC Curves](figures/touch_roc_curves.png)

### 3.4 マルチモーダル融合

![Multimodal Fusion Results](figures/multimodal_fusion.png)

| 融合戦略 | Accuracy | F1 | AUC-ROC |
|---|---|---|---|
| Average | 1.000 | 1.000 | 1.000 |
| Weighted Average | 1.000 | 1.000 | 1.000 |
| Meta-Learner (LR) | 1.000 | 1.000 | 1.000 |
| Meta-Learner (GB) | 1.000 | 1.000 | 1.000 |

### 3.5 縦断データの変化点検出

![Change Point Detection](figures/change_point_detection.png)

![CPD Method Comparison](figures/cpd_method_comparison.png)

| 手法 | Precision | Recall | F1 |
|---|---|---|---|
| CUSUM | 0.000 | 0.000 | 0.000 |
| PELT | 0.167 | 0.167 | 0.167 |
| Bayesian | 0.097 | 0.333 | 0.150 |

変化点検出は最も困難なタスクであり、パラメータ調整の余地が大きい。ベイジアン手法が最も高いRecallを示した。

### 3.6 臨床エンドポイントとの相関

![Clinical Validation](figures/clinical_validation.png)

デジタルバイオマーカーと臨床スコアの相関分析により、歩行スコア、音声スコア、タッチスコアのいずれも臨床評価との有意な相関を示した。特にコンポジットスコア（多モーダル融合）は単一モダリティよりも高い相関を達成した。

---

## 4. 考察と今後の展望

### 4.1 主要な知見

1. **歩行解析**は最も高い弁別性能を示し、パーキンソン病スクリーニングに有効
2. **音声解析**はALSの縦断的モニタリングに適しており、特にジッターとシマーが有力なバイオマーカー
3. **タッチスクリーン解析**は認知機能低下の検出において高精度
4. **マルチモーダル融合**は個々のモダリティよりもロバストな予測を提供
5. **変化点検出**は現在最も改善の余地がある領域

### 4.2 限界

- シミュレーションデータを使用しているため、実世界データでの検証が必要
- 対象間変動や環境要因の影響が十分にモデル化されていない
- 変化点検出のパラメータ最適化が不十分
- プライバシーとデータセキュリティの観点が未検討

### 4.3 今後の方向性

- 実臨床データを用いた検証（PhysioNet、mPowerデータセットなど）
- ディープラーニング（CNN-LSTM、Transformer）の導入
- フェデレーテッドラーニングによるプライバシー保護
- リアルタイム推論のためのエッジ展開最適化
- 大規模前向き臨床試験の設計

---

## 5. 生成したファイル一覧

### ソースコード
| ファイル | 説明 |
|---|---|
| `src/data_generation.py` | 合成センサーデータ生成 |
| `src/models.py` | 機械学習モデル（分類・融合） |
| `src/change_point_detection.py` | 変化点検出アルゴリズム |
| `src/run_experiments.py` | 実験実行・可視化メインスクリプト |

### 生成データ
| ファイル | 説明 |
|---|---|
| `data/gait_features.csv` | 歩行特徴量データ（200名） |
| `data/voice_features.csv` | 音声特徴量データ（150名×10セッション） |
| `data/touch_features.csv` | タッチ操作データ（180名×8セッション） |
| `data/longitudinal_data.csv` | 縦断マルチモーダルデータ（100名×52週） |

### 図表
| ファイル | 説明 |
|---|---|
| `figures/system_architecture.png` | システムアーキテクチャ図 |
| `figures/gait_model_comparison.png` | 歩行モデル比較・ROC曲線 |
| `figures/gait_feature_importance.png` | 歩行特徴量重要度 |
| `figures/gait_confusion_matrix.png` | 混同行列 |
| `figures/voice_analysis.png` | 音声特徴量分布・モデル比較 |
| `figures/voice_progression.png` | 音声特徴量の経時変化 |
| `figures/touch_analysis.png` | タッチ操作特徴量の群間比較 |
| `figures/touch_roc_curves.png` | 認知機能低下検出ROC曲線 |
| `figures/change_point_detection.png` | 変化点検出結果 |
| `figures/cpd_method_comparison.png` | CPD手法比較 |
| `figures/multimodal_fusion.png` | マルチモーダル融合結果 |
| `figures/clinical_validation.png` | 臨床バリデーション |

### ドキュメント
| ファイル | 説明 |
|---|---|
| `report.md` | 本レポート |
| `paper.md` | 学術論文形式の文書 |
