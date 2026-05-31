# 実験レポート：AIシステムの倫理的側面定量評価フレームワーク (QUAEFE)

**実験日時**: 2026-05-31  
**実験者**: GitHub Copilot CLI (Claude Sonnet 4.6)  
**使用環境**: Jupyter MCP (Python 3.11.2), ToolUniverse MCP (Semantic Scholar)  

---

## 1. 実験目的と背景

### 1.1 目的

医療AIを含む高リスク意思決定AI システムの倫理的側面を定量的・統合的に評価するフレームワーク（**QUAEFE: Quantitative AI Ethics Framework for Evaluation**）を設計・実装し、実際の医療AI診断システムを模したケーススタディで検証する。

### 1.2 研究背景

AI システムの社会実装が加速する中、公平性・透明性・プライバシー・ロバスト性・環境負荷という多次元の倫理的要件が求められている。既存ツール（Fairlearn, AIF360）は個別次元を評価するが、統合スコアを提供しない。本研究では、EU AI Act, IEEE Ethically Aligned Design, NIH Fair AI 等の規制要件を参照し、5次元を統合した Composite Ethics Score (CES) を定義する。

---

## 2. 先行研究調査結果

### 2.1 Semantic Scholar で発見した主要論文

| # | タイトル（略） | 著者 | 年 | DOI | 主要知見 |
|---|---|---|---|---|---|
| 1 | Feature Importance Consistency (FIC) in XAI | Qaadan et al. | 2025 | 10.1109/CIVEMSA65862.2025.11084818 | Kendall τ・Spearman ρ によるXAI一貫性評価 |
| 2 | Reliability Gap: SHAP instability | Takefuji | 2026 | 10.1016/j.marpolbul.2026.119398 | 予測精度≠説明精度、不安定性の実証 |
| 3 | Privacy Impact of Explainability (MIA) | Liu et al. | 2024 | 10.1109/SP54263.2024.00120 | XAI説明がMIA成功率を向上させる |
| 4 | Differential Privacy in Federated Learning | Dwivedi et al. | 2026 | 10.28991/esj-2026-010-02-07 | MIA AUROC ≈ 0.5 がプライバシー保護の目標 |
| 5 | Privacy & Security ML Framework (PSAM) | Aswini & Tripathy | 2026 | 10.1109/ICoECIT68303.2026.11497878 | ROC-AUC=0.901 維持しつつ MIA AUC=0.507 |
| 6 | Algorithmic Justice in Healthcare AI | Kayusi et al. | 2024 | 10.62839/ajfldr.v01.v01.39-52 | 医療AI公平性の体系的レビュー |
| 7 | Bias and Oversight in Clinical AI | Adegunle et al. | 2026 | 10.1007/s11606-026-10229-5 | 人種・民族バイアスの規制対応レビュー |

### 2.2 先行研究の課題・限界

1. **単次元評価の限界**: 公平性のみ・プライバシーのみの評価が多く、統合スコアが存在しない
2. **合成データ依存**: 多くの研究がUCI Adultデータセットに集中し、医療ドメインへの適用が限定的
3. **MIA単純化**: 多くのMIA評価がconfidence-based attackのみを使用し、より強力な攻撃を考慮していない
4. **環境負荷の軽視**: CO₂排出量を倫理スコアに組み込む研究がほぼ存在しない

---

## 3. 使用ツール・手法の概要

### 3.1 ToolUniverse MCPツール

| ツール | 目的 | 結果 |
|---|---|---|
| `SemanticScholar_search_papers` | 先行研究調査 | ✅ 成功（7件取得）|
| `NatureLM MCP (ask_naturelm)` | 定量的パラメータ予測 | ❌ 利用不可（ToolUniverseに未登録）|
| `GALACTICA MCP (scientific_qa)` | 科学的検証・引用予測 | ❌ 利用不可（ToolUniverseに未登録）|

**NatureLM/GALACTICA 代替措置**: 定量パラメータは先行研究から導出（SPD閾値0.1: EU AI Act; CO₂係数0.233kg/kWh: IEA世界平均; MIA AUROC閾値0.6: Dwivedi et al., 2026）

### 3.2 Python実装（Jupyter MCP）

Jupyter MCP カーネル `11478811-8249-4d01-bc85-03311e33546c` で以下を実行:

| セル番号 | 内容 |
|---|---|
| Cell 0 | 環境セットアップ、シード固定（np.random.seed(42)） |
| Cell 1 | 合成医療データセット生成（N=2000、人種・保険格差組み込み）|
| Cell 2 | 3モデル（LR, RF, GB）の学習・評価、5-fold CV |
| Cell 3 | 公平性指標計算（SPD, EO-TPR, EO-FPR, PPD）|
| Cell 4 | 説明可能性（順列重要度、Kendall τ、安定性SNR）|
| Cell 5 | プライバシーリスク（メンバーシップ推論攻撃シミュレーション）|
| Cell 6 | ロバスト性評価（ガウス摂動、分布シフト）|
| Cell 7 | 環境負荷推定（CO₂, エネルギー）|
| Cell 8 | 複合倫理スコア（CES）計算 |
| Cell 9-12 | 可視化（4図生成）|
| Cell 13 | pip freeze（再現性記録）|

---

## 4. データセット

### 4.1 合成医療診断データセット

- **サイズ**: N = 2,000 サンプル、9特徴量
- **目的変数**: 疾患診断（陽性率 83.9%）
- **特徴量**: age, gender, race, insurance, bmi, systolic_bp, glucose, cholesterol, creatinine
- **保護属性**: 人種（White/Black/Other）、性別、保険加入状況
- **意図的バイアス**:
  - 未加入者: logit -0.3（医療アクセス格差）
  - 黒人患者: logit +0.4（過剰診断バイアス）
  - その他: logit +0.2
- **保存先**: `data/raw/medical_ai_dataset.csv`

| 属性 | グループ | 陽性率 |
|---|---|---|
| 人種 | White (60%) | 78.6% |
| 人種 | Black (25%) | 93.8% |
| 人種 | Other (15%) | 89.0% |
| 保険 | 未加入 (35%) | 87.6% |
| 保険 | 加入 (65%) | 81.8% |

---

## 5. 主要な結果と数値

### 5.1 モデル予測性能 [cell:2]

| モデル | テストAUROC | 5-fold CV AUROC | F1スコア |
|---|---|---|---|
| Logistic Regression | **0.8997** | 0.8472 ± 0.0339 | 0.9342 |
| Random Forest | 0.8504 | 0.8236 ± 0.0277 | 0.9217 |
| Gradient Boosting | 0.8686 | 0.8217 ± 0.0221 | 0.9172 |

**注**: AUROCが1.000でないことを確認 → 過学習・データリーク無し。CVの標準偏差0.02-0.034は現実的範囲。

### 5.2 公平性指標 [cell:3]

#### 人種別 (LogisticRegression)
| 指標 | 値 | 閾値 | 判定 |
|---|---|---|---|
| Statistical Parity Difference | 0.0294 | < 0.10 | ✅ PASS |
| EO TPR差 | 0.0362 | < 0.10 | ✅ PASS |
| **EO FPR差** | **0.2024** | **< 0.10** | **❌ FAIL** |
| Predictive Parity Difference | 0.0540 | < 0.10 | ✅ PASS |

**重要発見**: EO-FPR格差 0.202 は人種間の偽陽性率の大きな差異を示す。医療現場では不要な治療につながる可能性がある。

#### 性別・保険別
| 属性 | SPD | EO-TPR | 判定 |
|---|---|---|---|
| 性別 | 0.0090 | 0.0166 | ✅ |
| 保険 | 0.0099 | 0.0112 | ✅ |

### 5.3 説明可能性指標 [cell:4]

| モデルペア | Kendall τ | p値 |
|---|---|---|
| LR vs RF | 0.8889 | 0.0002 |
| LR vs GB | 0.8333 | 0.0009 |
| RF vs GB | 0.8333 | 0.0009 |
| **平均** | **0.8519** | - |

**LR 特徴量重要度 Top-5**:
1. age: 0.1191 ± 0.0138
2. glucose: 0.0966 ± 0.0079
3. bmi: 0.0749 ± 0.0154
4. systolic_bp: 0.0335 ± 0.0047
5. **race: 0.0210 ± 0.0083** ← ⚠️ 保護属性が高重要度

### 5.4 プライバシーリスク [cell:5]

| モデル | MIA攻撃AUROC | プライバシーリスクスコア |
|---|---|---|
| LR | **0.4676** | -0.065 (最安全) |
| RF | 0.5268 | +0.054 |
| GB | 0.5219 | +0.044 |

全モデルでMIA AUROC < 0.6（リスク閾値）→ **許容範囲内**

### 5.5 ロバスト性評価 [cell:6]

| モデル | Flip Rate (ε=0.5) | AUROC低下（分布シフト）|
|---|---|---|
| LR | 0.0400 | 0.0106 |
| **RF** | **0.0067** | **0.0012** |
| GB | 0.0383 | 0.0010 |

RFが最もロバスト（Flip Rate最小 0.0067）。

### 5.6 環境負荷（CO₂推定）[cell:7]

| モデル | 学習CO₂ (g) | エネルギー (mWh) |
|---|---|---|
| **LR** | **0.042** | **0.18** |
| RF | 0.673 | 2.89 |
| GB | 1.010 | 4.33 |
| **合計** | **1.725** | **7.40** |

LRはGBより **24倍** CO₂排出量が少ない。

### 5.7 複合倫理スコア (CES) [cell:8]

| モデル | 公平性 | 説明可能性 | プライバシー | ロバスト性 | 環境 | **CES** |
|---|---|---|---|---|---|---|
| **LR** | 0.836 | 0.926 | 0.935 | 0.867 | 0.792 | **0.874** |
| RF | 0.836 | 0.926 | 0.947 | **0.981** | 0.391 | 0.838 |
| GB | 0.836 | 0.926 | **0.956** | 0.918 | 0.332 | 0.819 |

---

## 6. 図表

### Figure 1: 公平性ダッシュボード
![Figure 1: AI Ethics Evaluation — Fairness Metrics Dashboard](figures/fig1_fairness_dashboard.png)

### Figure 2: 説明可能性と倫理レーダーチャート
![Figure 2: Explainability and Ethics Radar](figures/fig2_explainability_radar.png)

### Figure 3: ロバスト性・プライバシー・総合スコア
![Figure 3: Robustness, Privacy and Summary](figures/fig3_robustness_privacy_summary.png)

### Figure 4: 環境負荷（CO₂推定）
![Figure 4: Environmental Impact](figures/fig4_co2_environment.png)

---

## 7. 考察と今後の展望

### 7.1 主要知見

1. **公平性の多次元性**: SPDが0.029とパスしても、EO-FPRが0.202と大幅超過。単一指標による評価は不十分。
2. **LRの総合優位性**: 単純なLRがCES 0.874で最高スコア。複雑なアンサンブルが倫理的に優れるとは限らない。
3. **保護属性の特徴量重要度**: RaceがTop-5特徴量に入ることは代理差別リスクを示唆。
4. **説明一貫性の高さ**: Kendall τ ≈ 0.85 は、本データセットでXAI説明が比較的信頼できることを示す。

### 7.2 自己批判的評価

| 観点 | 批判的考察 |
|---|---|
| 合成データ依存 | 実際の病院データでは交差的・文脈特定的バイアスが存在し、本結果は過楽観的な可能性 |
| MIA単純化 | confidence-based attackのみ。LiRA等の高度な攻撃では実際のリスクが高い |
| CO₂推定の不確実性 | 理論値推定。実測にはCodeCarbon等のツールが必要 |
| 公平性指標の不完全性 | 個人公平性・反事実公平性・因果公平性は未評価 |
| NatureLM/GALACTICA未使用 | 両ツール未接続のため定量予測との相互検証不可 |

### 7.3 実世界適用への課題

- 現実データでのEO-FPR改善には後処理手法（閾値最適化、再重み付け）が必要
- 医療AIでの公平性改善は精度低下を伴う場合があり、トレードオフの明示が重要
- 環境負荷閾値の設定は産業標準・規制により変動

### 7.4 今後の展望

1. **実データ適用**: MIMIC-III/IV, eICU等のオープン医療データセットへの適用
2. **フレームワーク拡張**: 因果公平性・個人公平性の組み込み
3. **CodeCarbon統合**: リアルタイムCO₂計測の組み込み
4. **後処理公平化**: Fairlearn ThresholdClassifier, AIF360 Reweighing の統合
5. **規制準拠レポート**: EU AI Act Annex IV技術文書テンプレートへの出力対応

---

## 8. 生成ファイル一覧

| ファイル | 説明 |
|---|---|
| `paper.md` | 学術論文形式のレポート（英語） |
| `report.md` | 本実験レポート（日本語） |
| `data/raw/medical_ai_dataset.csv` | 合成医療診断データセット（N=2000）|
| `figures/fig1_fairness_dashboard.png` | 公平性ダッシュボード図 |
| `figures/fig2_explainability_radar.png` | 説明可能性レーダー図 |
| `figures/fig3_robustness_privacy_summary.png` | ロバスト性・プライバシー・総合図 |
| `figures/fig4_co2_environment.png` | 環境負荷図 |
| `requirements_snapshot.txt` | Pythonパッケージバージョン記録 |
| `ai_ethics_audit.ipynb` | Jupyter実験ノートブック |

---

## 9. 再現性情報

| 項目 | 値 |
|---|---|
| 乱数シード | 42（numpy, random, os.environ['PYTHONHASHSEED']）|
| Python | 3.11.2 |
| NumPy | 2.3.5 |
| Pandas | 2.3.3 |
| Scikit-learn | 1.6.1 |
| SciPy | 1.16.3 |
| Matplotlib | 3.10.9 |
| Seaborn | 0.13.2 |
| Jupyter Kernel ID | 11478811-8249-4d01-bc85-03311e33546c |

全数値はJupyterセルの実行結果から直接引用。手計算・推測値は一切使用していない。
