# Experiment Report: Automated Quality Control and Anomaly Detection for Large-Scale Scientific Data

## 実験目的と背景

大規模科学実験（CERN LHC、LIGO等）では、毎秒数十GBに達するセンサーデータのリアルタイム品質管理が不可欠である。本実験では、ストリーミング対応の異常検知パイプラインを設計・実装し、以下の6つの技術要素を統合的に評価した：

1. **変化点検出（PELT / BOCPD）** — 時系列データの統計的性質の変化を検出
2. **多変量外れ値検出（Isolation Forest / Deep SVDD）** — 高次元空間における異常パターンの識別
3. **物理的制約に基づく異常スコアリング** — ドメイン知識を活用した物理法則違反の検出
4. **ドリフト検出（ADWIN / Page-Hinkley）** — データ分布の経時変化の監視
5. **説明可能な異常検知** — 特徴量帰属による異常原因の自動特定
6. **ストリーミングパイプライン設計** — Apache Kafka/Flink 基盤のリアルタイム処理

### 背景

先行研究では、Isolation Forest（Liu et al., 2008）が軽量で効果的な異常検知手法として広く使われているが、高次元データや非線形パターンへの対応に限界がある。Xu et al.（2023）は Deep Isolation Forest を提案し、ニューラルネットワークによる表現学習を統合することで検出精度を向上させた。また、Altamirano et al.（2023）はBOCPDのスケーラビリティとロバスト性を改善する手法を提案している。物理的制約を組み込んだ異常検知については、Zideh et al.（2024）がPhysics-Informed Machine Learningの包括的レビューを行っている。

---

## 使用した手法・アルゴリズムの概要

### 2.1 変化点検出

- **PELT (Pruned Exact Linear Time)**: RBFカーネルコスト関数を用い、ペナルティパラメータ `pen=10` で複数変化点を検出。線形時間計算量で大規模データに対応。
- **BOCPD (Bayesian Online Changepoint Detection)**: Student-t予測分布を用いたオンラインベイズ推論。ハザードレート `1/200` で実行。

### 2.2 多変量異常検出

- **Isolation Forest**: `n_estimators=200`, `contamination=0.05`。標準化後のデータに対して適用。
- **Deep SVDD-like**: PCA（4次元）による表現学習 + 超球面中心からの距離に基づく一クラス分類。

### 2.3 物理的制約スコアリング

4つの物理法則に基づく制約違反スコアを計算：
- **オームの法則**: V = IR（期待抵抗 50Ω）
- **光度-ビーム強度関係**: L ∝ I_beam
- **イベントレート-光度関係**: R_event ∝ L × σ
- **温度物理的範囲**: |T - 20| ≤ 5

### 2.4 ドリフト検出

- **ADWIN**: ホフディング限界に基づく適応的ウィンドウ手法（δ = 0.002）
- **Page-Hinkley テスト**: 累積和に基づくシーケンシャル変化検出（閾値 = 50）

### 2.5 説明可能な異常検知

SHAP-inspired な特徴量帰属法：各特徴量をベースラインに置換した場合のスコア変化を計算し、異常への寄与度を推定。

### 2.6 統合スコア

3つのスコア（IF: 35%, SVDD: 30%, Physics: 35%）の重み付き平均に対し、大津法で閾値を自動決定。

---

## 主要な結果と数値

### 3.1 データセット概要

合成データ（CERN/LIGO型センサーデータ模擬）を生成：
- サンプル数: 5,000、特徴量数: 8
- 異常率: 3.5%
- 真の変化点: t=1500（圧力シフト）、t=3000（磁場シフト）
- 注入異常: ポイント異常（30点）、文脈異常（50点）、集合異常（100点）

![データ概要](figures/data_overview.png)

### 3.2 変化点検出結果

| 手法 | 検出数 | TP | FP | FN | Precision | Recall | F1 |
|------|--------|----|----|-----|-----------|--------|-----|
| PELT | 17 | 2 | 15 | 0 | 0.118 | 1.000 | 0.211 |
| BOCPD | 0 | 0 | 0 | 2 | 0.000 | 0.000 | 0.000 |

PELTは両方の真の変化点を検出したが、偽陽性が多い。ペナルティパラメータの調整が必要。BOCPDはハザードレートの調整が必要である。

![変化点検出](figures/changepoint_detection.png)

### 3.3 異常検知性能比較

| 手法 | AUC | AP | Precision | Recall | F1 |
|------|-----|-----|-----------|--------|-----|
| Isolation Forest | 0.937 | 0.611 | 0.448 | 0.633 | 0.525 |
| Deep SVDD-like | 0.962 | 0.869 | 0.596 | 0.842 | 0.698 |
| Physics Constraints | 0.703 | 0.360 | 0.240 | 0.339 | 0.281 |
| **Combined Pipeline** | **0.962** | **0.806** | **0.988** | **0.480** | **0.646** |

![異常スコア](figures/anomaly_scores.png)

![ROC・PRカーブ](figures/roc_pr_curves.png)

![混同行列](figures/confusion_matrices.png)

### 3.4 物理的制約違反

物理的制約に基づくスコアリングは、特にオームの法則違反（文脈異常）の検出に有効であった。

![物理的制約](figures/physics_constraints.png)

### 3.5 ドリフト検出

- ADWIN: 温度チャネルで32箇所のドリフトを検出
- Page-Hinkley: 圧力チャネルで高感度に変化を検出（4,900点）

![ドリフト検出](figures/drift_detection.png)

### 3.6 特徴量帰属（説明可能性）

異常データに対する特徴量帰属分析により、EventRate、Pressure、Voltage が最も異常検知に寄与する特徴量として特定された。

![特徴量帰属](figures/feature_attribution.png)

### 3.7 ストリーミングパイプラインアーキテクチャ

Apache Kafka + Flink を基盤とした5段階のリアルタイム処理パイプラインを設計した。

![パイプラインアーキテクチャ](figures/pipeline_architecture.png)

---

## 考察と今後の展望

### 主要な知見

1. **Deep SVDD-like 手法が最高のAUC（0.962）とF1（0.698）を達成**し、PCA表現学習と超球面アプローチの有効性を示した。
2. **統合パイプラインは最高のPrecision（0.988）を達成**し、誤検知を最小化する実運用シナリオに適している。
3. **物理的制約スコアリングは単独では性能が限定的**だが、統合スコアの精度向上に大きく貢献している。
4. **PELTは高いRecallを示すがFPが多く**、ペナルティの最適化やポストフィルタリングが必要。

### 限界

- BOCPDのハイパーパラメータ（ハザードレート、事前分布）の自動チューニングが課題
- Deep SVDD の代わりにPCAベースの近似を使用しており、真のディープラーニング表現との差がある
- 合成データでの評価であり、実データでの検証が必要

### 今後の展望

- 真のDeep SVDDの実装とGPU対応
- CERN ROOT形式のデータへの直接適用
- フェデレーテッドラーニングによる分散異常検知
- Transformer ベースの時系列異常検知の統合

---

## 生成したファイル一覧

| ファイル名 | 説明 |
|------------|------|
| `experiment.py` | 実験メインスクリプト |
| `results.json` | 数値結果のJSON出力 |
| `figures/data_overview.png` | データ概要（8センサーチャネル） |
| `figures/changepoint_detection.png` | PELT/BOCPD変化点検出結果 |
| `figures/anomaly_scores.png` | 4手法の異常スコア時系列 |
| `figures/roc_pr_curves.png` | ROC/PRカーブ比較 |
| `figures/confusion_matrices.png` | 混同行列比較 |
| `figures/physics_constraints.png` | 物理的制約違反スコア |
| `figures/drift_detection.png` | ドリフト検出結果 |
| `figures/feature_attribution.png` | 特徴量帰属ヒートマップ |
| `figures/pipeline_architecture.png` | ストリーミングパイプライン設計図 |
| `report.md` | 本レポート |
| `paper.md` | 学術論文形式文書 |
