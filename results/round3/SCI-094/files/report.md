# AIシステム倫理監査フレームワーク — 実験レポート

> DRAFT — NOT FOR DISTRIBUTION

## 実験目的と背景

本実験は、AIシステムの倫理的側面を定量的かつ総合的に評価するフレームワーク（EAIF: Ethics Audit & Integration Framework）の設計・実装・検証を目的とする。近年、AIシステムが医療・金融・司法等の高リスク意思決定領域に普及するにつれ、単一の性能指標（精度・AUC）では不十分であり、公平性・説明可能性・プライバシー・ロバスト性・環境負荷を統合的に評価する枠組みが求められている。本研究では、Fairlearn・SHAP・scikit-learnを基盤として、これら5次元の倫理評価を単一の合成スコアに統合するパイプラインを構築し、合成医療診断データセットを用いた5モデル比較ケーススタディを実施した。

## 先行研究調査（ToolUniverse MCP / CrossRef API）

### 試行したAPIと結果

| API | 試行結果 |
|-----|---------|
| Semantic Scholar API | **429 Too Many Requests**（レート制限）— API keyなし環境では利用不可 |
| CrossRef REST API | **成功**（200 OK）— 30件以上の関連論文を取得 |
| PubMed / ToolUniverse MCP | 接続環境なし（.mcp.json未設定） |

科学的透明性のため、上記試行状況を記録する。代替手段としてCrossRef APIを使用した。

### 特定された主要先行研究

| # | タイトル | 著者 | 年 | DOI | 主要知見 |
|---|---------|------|-----|-----|---------|
| 1 | Metrics of Bias and Fairness in AI | Sekh & Prasad | 2024 | 10.1145/3690407 | 統計的平等・機会均等・校正の3指標統合評価手法を提案 |
| 2 | Equalized odds is a requirement of algorithmic fairness | Hedden | 2023 | 10.1007/s11229-023-04054-0 | EOを哲学的・数学的に正当化し法的フレームワークと接続 |
| 3 | Differential Privacy Protection Against MIA on Deep Learning | Rahman et al. | 2020 | 10.1142/9789811232701_0003 | 差分プライバシーがMIA耐性を有意に改善（AUC→0.52） |
| 4 | Advancing Fairness in Clinical AI Decision-Making | Rajpurkar et al. | 2026 | 10.5220/0014304300004070 | 臨床AIにおける公平性評価の社会技術的フレームワーク |
| 5 | Incorporating Ethics into the AI Clinical Decision Support System Lifecycle | Char et al. | 2024 | 10.2307/jj.14491770.11 | 医療AI倫理審査をSDLCに組み込むガイドラインを提示 |
| 6 | Green AI: Strategies for Mitigating Carbon Footprints | Naeem et al. | 2026 | 10.21275/sc26211100444 | ML訓練のCO2排出量推定モデルとグリーンAI戦略 |
| 7 | Integrating Feature Selection, ML, and SHAP Explainability | Almazni et al. | 2025 | 10.3390/diagnostics15192473 | SHAPの安定性・一貫性が診断AIの信頼性向上に寄与 |

### 先行研究の課題・限界

1. **次元間統合の欠如**: 大半の研究が単一の倫理次元（公平性のみ、またはプライバシーのみ）を扱い、5次元を統合した合成スコアが存在しない。
2. **医療特有の評価不足**: 医療診断AIに特化したエンドツーエンドの倫理監査パイプラインは限定的。
3. **環境負荷の軽視**: CO2排出量を倫理評価の一軸として扱う研究はごく少数。
4. **ベースライン比較の不備**: 単一モデルの評価が多く、アーキテクチャ間の比較が不十分。

## 使用した手法・アルゴリズムの概要

### フレームワーク構成

**EAIF（Ethics Audit & Integration Framework）** は以下の5モジュールで構成される：

#### 1. 公平性評価モジュール（`fairness_metrics.py`）

3指標を統合：

$$\text{DP-diff} = |\Pr(\hat{Y}=1 \mid A=1) - \Pr(\hat{Y}=1 \mid A=0)|$$

$$\text{EO-TPR-diff} = |\text{TPR}_{A=1} - \text{TPR}_{A=0}|, \quad \text{EO-FPR-diff} = |\text{FPR}_{A=1} - \text{FPR}_{A=0}|$$

$$\text{ECE} = \frac{1}{B}\sum_{b=1}^{B} |\bar{p}_b - \bar{y}_b|$$

合成公平性スコア：

$$S_{\text{fair}} = 1 - \text{clip}\!\left(3.5 \cdot \text{DP-diff} + 3.5 \cdot \max(\text{EO-TPR-diff}, \text{EO-FPR-diff}) + 6.0 \cdot \text{ECE-diff},\ 0,\ 1\right)$$

#### 2. 説明可能性定量化（`explainability_privacy_robustness.py`）

SHAP一貫性をサブサンプル間ピアソン相関で測定：

$$S_{\text{expl}} = 0.5 \cdot \rho(\text{SHAP}_1, \text{SHAP}_2) + 0.3 \cdot \text{RankStability}_{k} + 0.2 \cdot (1 - \text{Sparsity})$$

#### 3. プライバシーリスクスコア（MIA耐性）

信頼スコアギャップヒューリスティックによるメンバーシップ推論攻撃耐性：

$$\text{MIA-Advantage} = 2 \cdot |\text{AUC}_{\text{shadow}} - 0.5|, \quad S_{\text{priv}} = 1 - \text{MIA-Risk}$$

#### 4. ロバスト性評価

$$S_{\text{rob}} = 1 - \frac{\max(\Delta_{\text{FGSM}},\, \Delta_{\text{Gaussian}},\, \Delta_{\text{shift}})}{0.30}$$

#### 5. 環境負荷（CO2推定）

$$E_{\text{kWh}} = \frac{P_{\text{TDP}} \times n_{\text{CPU}} \times t_{\text{train}} \times \text{PUE}}{3.6 \times 10^6}, \quad \text{CO}_2 = E_{\text{kWh}} \times 475 \text{ [gCO}_2\text{eq/kWh]}$$

#### 合成スコア（重み付き統合）

$$S_{\text{composite}} = 0.30 \cdot S_{\text{fair}} + 0.20 \cdot S_{\text{expl}} + 0.20 \cdot S_{\text{priv}} + 0.20 \cdot S_{\text{rob}} + 0.10 \cdot S_{\text{env}}$$

## 主要な結果

### 5-分割交差検証によるモデル性能

| モデル | AUC (mean ± std) | F1 (mean ± std) | Acc (mean) | 合成倫理スコア | リスクレベル |
|-------|-----------------|----------------|------------|--------------|------------|
| Logistic Regression | **0.717 ± 0.035** | 0.599 ± 0.041 | 0.665 | **0.865** | LOW ✅ |
| Random Forest | 0.705 ± 0.024 | 0.548 ± 0.039 | 0.665 | 0.718 | MEDIUM ⚠️ |
| Gradient Boosting | 0.686 ± 0.028 | 0.588 ± 0.033 | 0.657 | 0.535 | HIGH ❌ |
| SVM | 0.697 ± 0.034 | 0.561 ± 0.034 | 0.649 | 0.708 | MEDIUM ⚠️ |
| MLP | 0.648 ± 0.029 | 0.568 ± 0.018 | 0.615 | 0.699 | MEDIUM ⚠️ |

### 倫理次元別スコア

| モデル | 公平性 | 説明可能性 | プライバシー保護 | ロバスト性 | 環境 |
|-------|------|-----------|--------------|---------|------|
| Logistic Regression | 0.638 | 0.944 | **0.993** | 0.933 | **1.000** |
| Random Forest | 0.601 | **0.990** | 0.276 | 0.922 | **1.000** |
| Gradient Boosting | 0.220 | 0.960 | 0.075 | 0.811 | 0.999 |
| SVM | 0.330 | 0.983 | 0.607 | **0.956** | **1.000** |
| MLP | 0.408 | 0.926 | 0.470 | **0.989** | 0.999 |

### 主要発見

1. **Logistic Regression が最高合成スコア（0.865, LOW risk）**: 高い透明性・低い過学習によりプライバシーリスクが最低（MIA AUC=0.502, chance level）。
2. **Gradient Boosting は HIGH risk（スコア0.535）**: 公平性が最低（0.220）かつMIA AUC=0.731で深刻なプライバシー漏洩リスク。
3. **性能と倫理はトレードオフ**: AUC最高のモデル（LR: 0.717）が同時に最高の倫理スコアを示し、複雑モデル（GB）が最低。

### 生成した図

![Figure 1: 倫理次元レーダーチャート](figures/fig1_radar_chart.png)

![Figure 2: 合成倫理スコアのモデル別比較](figures/fig2_composite_bar.png)

![Figure 3: 倫理指標ヒートマップ](figures/fig3_ethics_heatmap.png)

![Figure 4: 公平性指標詳細](figures/fig4_fairness_breakdown.png)

![Figure 5: 性能 vs 倫理スコアのトレードオフ](figures/fig5_perf_ethics_tradeoff.png)

![Figure 6: モデル別CO2排出量](figures/fig6_co2_emissions.png)

## 考察と今後の展望

### 考察

**公平性次元**: 全5モデルでDP差分（最大0.077）がFairlearnの推奨閾値（0.10）以下であったが、EO-TPR差分ではSVM（0.110）とMLP（0.107）が閾値超過。ランダムフォレストのMIA AUC（0.681）とGBのMIA AUC（0.731）は、訓練データの記憶によるプライバシー漏洩リスクを示す。SHAP一貫性は全モデルで高く（0.927–0.999）、説明可能性次元は主要な差別化要因でないことが判明した。

**パフォーマンスと倫理のトレードオフ**: Figure 5が示すように、高性能モデル（AUC高）が必ずしも倫理スコアが低いわけではない。Logistic Regressionは最高AUC（0.717）と最高倫理スコア（0.865）を同時に達成しており、単純モデルが複雑モデルより医療診断AIに適している可能性を示す。

**環境負荷**: Gradient Boostingの訓練時間（0.497秒）はLogistic Regression（0.0016秒）の約310倍。大規模データセットでは環境差が顕著になる。

### 今後の展望

1. 実際の医療データ（MIMIC-III等）を用いた検証
2. Fairlearn・AIF360の再重み付け手法（Reweighing, ExponentiatedGradient）による公平性改善
3. 差分プライバシー（DP-SGD）の統合による高リスクモデルのプライバシー強化
4. 動的な重み設定（ドメイン別の倫理優先度反映）
5. リアルタイム監視ダッシュボードの実装

## 生成したファイル一覧

| ファイル | 種類 | 説明 |
|---------|------|------|
| `src/fairness_metrics.py` | Python (117行) | 公平性評価モジュール |
| `src/explainability_privacy_robustness.py` | Python (195行) | 説明可能性・プライバシー・ロバスト性 |
| `src/environmental_audit.py` | Python (120行) | CO2推定・合成スコア計算 |
| `src/pipeline.py` | Python (296行) | メイン実験パイプライン |
| `tests/test_ethics_framework.py` | Python (69行) | ユニットテスト（7件、全通過） |
| `results/ethics_results.csv` | CSV | 全モデルの定量評価結果 |
| `figures/fig1_radar_chart.png` | PNG | 倫理次元レーダーチャート |
| `figures/fig2_composite_bar.png` | PNG | 合成スコア棒グラフ |
| `figures/fig3_ethics_heatmap.png` | PNG | 倫理指標ヒートマップ |
| `figures/fig4_fairness_breakdown.png` | PNG | 公平性詳細比較 |
| `figures/fig5_perf_ethics_tradeoff.png` | PNG | 性能vs倫理トレードオフ |
| `figures/fig6_co2_emissions.png` | PNG | CO2排出量比較 |
| `logs/process-log.jsonl` | JSONL | 実行トレースログ |

## 参考文献

1. Sekh, A.A., & Prasad, D.K. (2024). Metrics of Bias and Fairness in AI. DOI: 10.1145/3690407
2. Hedden, B. (2023). Equalized odds is a requirement of algorithmic fairness. *Synthese*, 202, 1–21. DOI: 10.1007/s11229-023-04054-0
3. Rahman, M.A. et al. (2020). Differential Privacy Protection Against Membership Inference Attack on Machine Learning. DOI: 10.1142/9789811232701_0003
4. Rajpurkar, P. et al. (2026). Advancing Fairness in Clinical AI Decision-Making through a Sociotechnical Framework. DOI: 10.5220/0014304300004070
5. Char, D.S. et al. (2024). Incorporating Ethics into the AI Clinical Decision Support System Lifecycle. DOI: 10.2307/jj.14491770.11
6. Naeem, M. et al. (2026). Green AI: Strategies for Mitigating Carbon Footprints in Machine Learning. DOI: 10.21275/sc26211100444
7. Almazni, N. et al. (2025). Integrating Feature Selection, Machine Learning, and SHAP Explainability. *Diagnostics*, 15(19), 2473. DOI: 10.3390/diagnostics15192473
8. Lundberg, S.M., & Lee, S.I. (2017). A Unified Approach to Interpreting Model Predictions. *NeurIPS*. DOI: 10.48550/arXiv.1705.07874
9. Hardt, M., Price, E., & Srebro, N. (2016). Equality of Opportunity in Supervised Learning. *NeurIPS*. DOI: 10.48550/arXiv.1610.02413
10. Shokri, R. et al. (2017). Membership Inference Attacks Against Machine Learning Models. *IEEE S&P*. DOI: 10.1109/SP.2017.41
11. Bird, S. et al. (2020). Fairlearn: A Toolkit for Assessing and Improving Fairness in AI. *Microsoft Research*. DOI: 10.48550/arXiv.2006.02424
12. Patterson, D. et al. (2021). Carbon and the Broad Landscape of Digital Technology. DOI: 10.48550/arXiv.2007.10392
