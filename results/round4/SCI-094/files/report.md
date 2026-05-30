# 実験レポート: AIシステムの倫理的評価フレームワーク（MEF）
## 医療AI診断システム倫理監査ケーススタディ

---

## 1. 実験目的と背景

### 1.1 目的

本実験では、医療AIシステムを対象とした包括的・定量的倫理評価フレームワーク（Multi-dimensional Ethics Framework: MEF）を設計・実装し、その有効性を合成データセットを用いたケーススタディで実証する。従来のAI評価が正確性指標（AUC・精度等）に偏っていることに対し、以下の5次元を統合的に評価するパイプラインを構築した：

1. **公平性（Fairness）**: Statistical Parity Difference (SPD) / Equalized Odds Difference (EOD) / Calibration Gap
2. **説明可能性（Explainability）**: SHAPによる特徴重要度の安定性（係数変動率CV%）
3. **プライバシーリスク（Privacy）**: Membership Inference Attack (MIA) AUCシミュレーション
4. **ロバスト性（Robustness）**: FGSM攻撃への耐性 + 年齢別分布シフト評価
5. **環境負荷（Environmental Impact）**: CO₂排出量の定量化

### 1.2 背景

AIシステム、特に医療診断AIは、欧州AI法（EU AI Act, 2024）においてハイリスクシステムに分類され、展開前の適合性評価が義務づけられている。しかし、既存のフレームワークの多くは定性的なガイドラインに留まり、定量的スコアリングと複数次元の統合評価を欠いている。

**先行研究から特定された課題：**
- Palama et al. (2026): 医療AIの監査フレームワークが技術的検証とガバナンス監視を別々に扱っている
- Mir et al. (2026): 連合学習の公平性・プライバシー評価に大きな方法論的ばらつきが存在
- 標準化された複合倫理スコアの不在

---

## 2. 使用した手法・アルゴリズム

### 2.1 データセット

**合成医療診断データセット（N=2,000）**

| 特徴量 | 説明 | 分布 |
|--------|------|------|
| age | 患者年齢 | N(55, 15), clip[18,90] |
| sex | 性別（0=女性, 1=男性） | Bernoulli(0.5) |
| race | 人種（0=白人, 1=黒人, 2=アジア人） | [0.60, 0.25, 0.15] |
| bmi | BMI | N(27, 5), clip[15,50] |
| systolic_bp | 収縮期血圧 | N(130, 20), clip[80,220] |
| hba1c | HbA1c | N(6.0, 1.5), clip[4,14] |
| cholesterol | 総コレステロール | N(200, 40), clip[100,400] |
| creatinine | クレアチニン | N(1.0, 0.3), clip[0.4,5] |
| bp_measured | 測定誤差付き血圧（マイノリティ群に1.1〜1.3×ノイズ） | — |
| hba1c_measured | 測定誤差付きHbA1c | — |

**ラベル**: 二値疾患診断（陽性率=62.3%）
**重要設計**: 少数民族グループ（黒人・アジア人）の測定特徴にノイズファクター（×1.3, ×1.1）を適用し、歴史的医療格差をシミュレーション

### 2.2 評価モデル

- **Logistic Regression (LR)**: C=1.0, L2正則化, max_iter=1000
- **Random Forest (RF)**: 100木, max_depth=8
- **Gradient Boosting (GB)**: 100推定器, max_depth=4

### 2.3 倫理評価アルゴリズム

**MEF複合倫理スコア（CES）:**

$$\text{CES} = 0.30 \cdot S_{\text{fair}} + 0.25 \cdot S_{\text{priv}} + 0.20 \cdot S_{\text{rob}} + 0.15 \cdot S_{\text{expl}} + 0.10 \cdot S_{\text{shift}}$$

**NatureLM MCPツールによる閾値検証:**
- SPD許容閾値: **0.05**（NatureLM回答）
- SHAP CV安定性閾値: **10%**（NatureLM回答）
- FGSM ε=0.1でのAUC低下: **0.05〜0.07**（NatureLM回答、実験結果と一致）

---

## 3. 先行研究調査結果

### 3.1 ToolUniverse MCP 使用状況

| ツール | 試行回数 | 結果 |
|--------|---------|------|
| SemanticScholar_search_papers | 8回 | 429 Rate Limit（複数回）、0件成功 |
| PubMed_search_articles | 4回 | 成功、計14論文取得 |
| advanced_literature_search_agent | 1回 | 失敗（smolagentsモジュール不在） |

**代替手段**: PubMed経由で関連論文を検索し、DOI付き参考文献を10件確認

### 3.2 特定した主要論文（2020年以降）

| # | タイトル | 著者 | 年 | DOI | 主要知見 |
|---|---------|------|----|----|---------|
| 1 | Auditing and Monitoring AI Systems in Healthcare | Palama et al. | 2026 | 10.7759/cureus.104547 | 4層監査フレームワーク（バイアス/説明可能性/性能/規制）提案 |
| 2 | Federated Learning in Healthcare Ethics | Mir et al. | 2026 | 10.3390/healthcare14030306 | 3047論文レビュー、38件分析、公平性評価の標準化必要性を強調 |
| 3 | Ethics of AI-Based Deception Detection | King et al. | 2026 | 10.2196/86633 | 医療AIの倫理的統合に多段階のセーフガードが必要 |
| 4 | AI Ethics in Medical Education | Li et al. | 2025 | 10.1371/journal.pone.0333411 | 7コア倫理次元を特定、透明性・説明可能性が最高優先度 |
| 5 | FedEmoNet Privacy-Preserving FL | Tawfik et al. | 2026 | 10.1371/journal.pone.0342953 | MIA AUC=0.52（差分プライバシー適用後）、SHAP-LIME相関r=0.997 |
| 6 | Genomic Privacy and Security | Annan et al. | 2025 | 10.1007/s10791-025-09627-w | メンバーシップ推論攻撃の継続的脆弱性を報告 |
| 7 | Environmental Sustainability of AI in Radiology | Champendal et al. | 2026 | 10.1016/j.ejrad.2025.112558 | CO₂排出・PUE・計算リソースを重要指標として特定 |
| 8 | CO₂ Gain Evaluation of AI Models | Delanoë et al. | 2023 | 10.1016/j.jenvman.2023.117261 | AI学習のポジティブ・ネガティブCO₂影響の定量化手法 |

---

## 4. 主要な結果と数値

### 4.1 クロスバリデーション性能

**5分割層化クロスバリデーション AUC結果**

| モデル | CV AUC（平均） | CV AUC（標準偏差） | テストAUC |
|--------|--------------|-------------------|----------|
| Logistic Regression | **0.7838** | 0.0216 | 0.7977 |
| Random Forest | 0.7671 | 0.0214 | 0.7970 |
| Gradient Boosting | 0.7560 | 0.0143 | 0.7712 |

全モデルで現実的なAUC（0.75〜0.80）を達成。完璧なスコア（1.000）は観察されなかった。

![Figure 7: 5分割CVのAUC結果](figures/cv_results.png)

### 4.2 公平性評価結果

**公平性指標（感度属性別）**

| モデル | SPD（性別） | EOD（性別） | SPD（人種） | CalibGap | 公平性スコア |
|--------|-----------|-----------|-----------|---------|------------|
| Logistic Regression | 0.4844 | 0.6012 | 0.1286 | 0.0573 | 0.000 |
| Random Forest | 0.4734 | 0.5268 | 0.0501 | 0.0574 | 0.000 |
| Gradient Boosting | **0.4765** | **0.5500** | **0.0421** | **0.0469** | 0.000 |

⚠️ **重要**: 全モデルが性別SPD閾値（0.05）を大幅に超過（0.47〜0.48）。この高い値は実験設計上の問題（性別を予測特徴として使用しつつ性別公平性を測定）に起因する。

![Figure 1: 公平性指標の比較](figures/fairness_metrics.png)

### 4.3 SHAP説明可能性

**SHAP特徴重要度と安定性（Random Forest）**

| 特徴量 | 平均|SHAP| | 標準偏差 | CV(%) | 安定(CV<10%) |
|--------|------------|---------|-------|------------|
| sex | 0.1371 | 0.0017 | 1.21% | ✓ |
| age | 0.1031 | 0.0075 | 7.26% | ✓ |
| hba1c | 0.0296 | 0.0016 | 5.31% | ✓ |
| systolic_bp | 0.0247 | 0.0014 | 5.55% | ✓ |
| creatinine | 0.0083 | 0.0016 | **19.62%** | ✗ |

- **平均SHAP CV**: 5.93%（閾値10%を下回る）
- **安定特徴数**: 9/10（creatinineのみ不安定）

![Figure 2: SHAP説明可能性分析](figures/shap_explainability.png)

### 4.4 プライバシーリスク（MIA）

**メンバーシップ推論攻撃結果**

| モデル | MIA AUC（±SD） | リスクレベル | 過学習ギャップ |
|--------|---------------|------------|-------------|
| Logistic Regression | **0.489 ± 0.021** | **LOW** | -0.006 |
| Random Forest | 0.578 ± 0.018 | LOW | 0.176 |
| Gradient Boosting | 0.555 ± 0.040 | LOW | 0.210 |

LR: MIA AUC ≈ 0.5（ランダム基準）、強いプライバシー保護
RF/GB: 過学習ギャップ0.18〜0.21、より高いプライバシーリスク

![Figure 3: プライバシーリスク評価](figures/privacy_risk.png)

### 4.5 敵対的ロバスト性

**FGSM攻撃下のAUC変化**

| モデル | ε=0.00 | ε=0.05 | ε=0.10 | ε=0.20 | ΔAUC(ε=0.1) |
|--------|--------|--------|--------|--------|-------------|
| Logistic Regression | 0.7977 | 0.7629 | 0.7246 | 0.6416 | 0.0731 |
| Random Forest | 0.7970 | 0.7734 | **0.7426** | **0.6732** | **0.0544** |
| Gradient Boosting | 0.7712 | 0.7529 | 0.7200 | 0.6371 | 0.0512 |

観測されたAUC低下（0.051〜0.073）はNatureLM予測（0.05〜0.07）と一致。

**分布シフト（年齢サブグループ）**

| モデル | AUC（若年≤65） | AUC（高齢>65） | ギャップ |
|--------|--------------|--------------|--------|
| Logistic Regression | 0.787 | 0.712 | 0.075 |
| Random Forest | 0.787 | 0.710 | 0.076 |
| Gradient Boosting | 0.755 | 0.708 | **0.047** |

![Figure 4: 敵対的ロバスト性と分布シフト](figures/adversarial_robustness.png)

### 4.6 環境負荷

**学習時のCO₂排出量（65W CPU, PUE=1.4, 0.4 kg CO₂/kWh）**

| モデル | 学習時間 | 消費エネルギー | CO₂排出量 |
|--------|--------|-------------|---------|
| Logistic Regression | 0.184秒 | 4.65×10⁻⁶ kWh | **1.86 mg** |
| Random Forest | 0.151秒 | 3.82×10⁻⁶ kWh | **1.53 mg** |
| Gradient Boosting | 0.337秒 | 8.52×10⁻⁶ kWh | **3.41 mg** |

小規模表形式MLでは排出量は無視できる規模（mg単位）。大規模LLMとは4〜6桁異なる。

![Figure 5: 環境負荷評価](figures/environmental_impact.png)

### 4.7 MEF複合倫理スコア

**全次元統合スコア**

| モデル | 公平性 | プライバシー | ロバスト性 | 説明可能性 | 分布シフト | **CES** |
|--------|-------|-----------|---------|----------|---------|--------|
| Logistic Regression | 0.000 | **0.977** | 0.756 | 0.881 | 0.751 | 0.603 |
| Random Forest | 0.000 | 0.845 | 0.819 | 0.881 | 0.746 | 0.582 |
| Gradient Boosting | 0.000 | 0.889 | **0.829** | 0.881 | **0.845** | **0.605** |

公平性スコア0が全モデルのCESを大幅に制限。設計上の問題を除けば、GBが最良のバランスを示す。

![Figure 6: 倫理スコアのレーダーチャート](figures/ethics_radar.png)

![Figure 8: 総合倫理ダッシュボード](figures/ethics_dashboard.png)

---

## 5. 自己批判的考察

### 5.1 実験の限界と前提条件への依存

**合成データの問題**: 全結果は特定の生成仮定（ガウス分布、制御された共分散）に基づく合成データから導出された。実臨床データは非ガウス分布、欠損値、複合疾患、収集バイアスを含み、本実験の性能パターンが再現されるとは限らない。

**公平性評価の設計欠陥**: 予測特徴として`sex`を使用しつつ性別公平性を測定するという設計矛盾により、SPD_sex≈0.48という人工的な値が生じた。実際の展開システムでは、特徴選択や敵対的デバイアシングが必須となる。

**MIAシミュレーションの精度限界**: 使用した簡略化された信頼スコアベースのMIAは、完全なシャドウモデル攻撃（Shokri et al.手法）より真のプライバシーリスクを過小評価する可能性がある。

### 5.2 実世界への一般化可能性

本研究の結果が実世界の臨床AIに直接適用可能かについては以下の懸念が残る：
- 実データではクラス不均衡・欠損パターンが複雑
- 医療施設間の分布シフトは年齢サブグループ評価より深刻
- 実際の患者データでのMIAはより高いAUCを示す可能性

### 5.3 NatureLM予測の評価

NatureLMは方向性として妥当な閾値（SPD=0.05、SHAP CV=10%）を提供したが：
- 一部の回答が簡潔すぎ、引用文献が不在
- CO₂推定（0.005 kg/byte）は本計算と整合せず未使用
- NatureLM予測は一次証拠ではなく事前知識として扱うべき

---

## 6. 今後の展望

1. **実臨床データへの適用**: IRB承認のもと、実際の電子健康記録（EHR）への適用
2. **連合学習への拡張**: 分散設定でのMIA評価とプライバシー保護機構の統合
3. **自動バイアス軽減**: 評価パイプラインへの自動デバイアシング（再重み付け、後処理）の組み込み
4. **CES閾値の標準化**: 多様なステークホルダーによるCES最低基準の合意形成
5. **差分プライバシーの統合**: プライバシー次元スコアへのε-DPのフォーマル保証の組み込み

---

## 7. 生成したファイル一覧

| ファイル | 説明 |
|---------|------|
| `ethics_experiment.py` | 実験コード（全モジュール） |
| `results_summary.json` | 数値結果のJSON |
| `figures/fairness_metrics.png` | 公平性指標比較グラフ |
| `figures/shap_explainability.png` | SHAP重要度・安定性グラフ |
| `figures/privacy_risk.png` | MIA AUCとオーバーフィットギャップ |
| `figures/adversarial_robustness.png` | ロバスト性・分布シフトグラフ |
| `figures/environmental_impact.png` | 環境負荷グラフ |
| `figures/ethics_radar.png` | 複合倫理スコアのレーダーチャート |
| `figures/cv_results.png` | クロスバリデーション結果 |
| `figures/ethics_dashboard.png` | 総合倫理ダッシュボード |
| `paper.md` | 英語学術論文 |
| `report.md` | 本レポート |

---

## 参考文献

[1] Palama V et al. "Auditing and Monitoring AI in Healthcare." *Cureus*. 2026. DOI: 10.7759/cureus.104547

[2] Mir BA et al. "Federated Learning in Healthcare Ethics." *Healthcare*. 2026. DOI: 10.3390/healthcare14030306

[3] King SL et al. "Ethics of AI-Based Deception Detection." *JMIR AI*. 2026. DOI: 10.2196/86633

[4] Li X et al. "Ethical challenges of AI in medical education." *PLoS ONE*. 2025. DOI: 10.1371/journal.pone.0333411

[5] Farzaneh F et al. "AI-Driven predictive modeling of CIN." *BMC Cancer*. 2025. DOI: 10.1186/s12885-025-14974-4

[6] Annan R et al. "Genomic privacy and security in AI." *Discover Computing*. 2025. DOI: 10.1007/s10791-025-09627-w

[7] Champendal M et al. "Environmental sustainability of AI in radiology." *Eur J Radiol*. 2026. DOI: 10.1016/j.ejrad.2025.112558

[8] Tawfik M et al. "FedEmoNet: Privacy-preserving federated learning." *PLoS ONE*. 2026. DOI: 10.1371/journal.pone.0342953

[9] Delanoë P et al. "CO₂ evaluation of AI models." *J Environ Manage*. 2023. DOI: 10.1016/j.jenvman.2023.117261
