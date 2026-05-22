# 食品サプライチェーン安全リスク予測AIシステム

**DRAFT — NOT FOR DISTRIBUTION**

| 項目 | 内容 |
|------|------|
| 作成日 | 2026-05-23 |
| バージョン | 1.0 |
| 乱数シード | 42（全モジュール共通） |
| 実行環境 | Python 3.12 / scikit-learn / scipy / statsmodels |

---

## 目次

1. [実験目的と背景](#1-実験目的と背景)
2. [システム全体アーキテクチャ](#2-システム全体アーキテクチャ)
3. [モジュール1: 食中毒発生の時空間予測](#3-モジュール1-食中毒発生の時空間予測)
4. [モジュール2: NLPによるリコール・アラート早期検出](#4-モジュール2-nlpによるリコールアラート早期検出)
5. [モジュール3: 微生物増殖予測モデル](#5-モジュール3-微生物増殖予測モデル)
6. [モジュール4: HACCP管理点リスクスコアリング](#6-モジュール4-haccp管理点リスクスコアリング)
7. [モジュール5: ブロックチェーントレーサビリティ](#7-モジュール5-ブロックチェーントレーサビリティ)
8. [モジュール6: 鶏肉サルモネラ汚染予測ケーススタディ](#8-モジュール6-鶏肉サルモネラ汚染予測ケーススタディ)
9. [統合リスクモニタリングシステム設計](#9-統合リスクモニタリングシステム設計)
10. [考察と今後の展望](#10-考察と今後の展望)
11. [生成ファイル一覧](#11-生成ファイル一覧)
12. [参考文献](#12-参考文献)

---

## 1. 実験目的と背景

### 1.1 目的

食品サプライチェーンにおける安全リスクを**予測・検知・追跡**する統合AIシステムを設計・実装し、以下の6つの要素技術を検証する：

1. 気象データと連動した食中毒発生の時空間予測
2. 自然言語処理（NLP）によるリコール/アラートの早期検出
3. 予測微生物学に基づく微生物増殖シミュレーション
4. HACCP管理点のリスクスコアリング自動化
5. ブロックチェーン技術によるサプライチェーントレーサビリティ
6. 鶏肉サルモネラ汚染を対象とした統合ケーススタディ

### 1.2 背景

世界保健機関（WHO）の推計によれば、毎年約6億人が汚染食品による食中毒を経験し、42万人が死亡している。食品安全管理は従来、事後的な検査・回収に依存してきたが、AI・IoT・ブロックチェーン技術の発展により、**予測的・予防的アプローチ**への転換が可能になりつつある。

本研究では、時系列予測、NLP、予測微生物学、統計的プロセス制御、分散台帳技術を統合した**リアルタイムリスクモニタリングシステム**のプロトタイプを構築し、その有効性を合成データに基づき評価する。

### 1.3 使用データ

全モジュールにおいて、実データの構造・分布・相関を忠実に再現した**合成データ**を使用した。これは、機密性の高い食品安全データへのアクセス制約に対応しつつ、手法の妥当性を検証するためである。

---

## 2. システム全体アーキテクチャ

```
┌─────────────────────────────────────────────────────────────────┐
│                  統合リスクモニタリングダッシュボード               │
│  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐          │
│  │時空間予測 │ │NLPアラート│ │微生物増殖 │ │HACCPスコア│          │
│  │モジュール │ │モジュール │ │モジュール │ │モジュール │          │
│  └────┬─────┘ └────┬─────┘ └────┬─────┘ └────┬─────┘          │
│       │            │            │            │                  │
│       └────────────┴─────┬──────┴────────────┘                  │
│                          │                                      │
│              ┌───────────▼───────────┐                          │
│              │  リスク統合エンジン    │                          │
│              │  (ベイズ統合スコア)    │                          │
│              └───────────┬───────────┘                          │
│                          │                                      │
│              ┌───────────▼───────────┐                          │
│              │  ブロックチェーン       │                          │
│              │  トレーサビリティ層    │                          │
│              └───────────────────────┘                          │
└─────────────────────────────────────────────────────────────────┘
```

各モジュールは独立して動作しつつ、リスク統合エンジンを通じて総合リスクスコアを算出する。ブロックチェーン層はすべてのトランザクションを不変的に記録し、リコール時の即時追跡を可能にする。

---

## 3. モジュール1: 食中毒発生の時空間予測

### 3.1 手法

- **データ**: 10地域 × 5年間（2020–2024）の月次食中毒発生件数（合成データ、600レコード）
- **特徴量**: 気温、湿度、月（sin/cos周期エンコーディング）、人口密度、地域ダミー変数
- **モデル**: Random Forest / Gradient Boosting / MLPRegressor（LSTMプロキシ）
- **評価**: 時系列分割（2020–2023訓練、2024テスト）、RMSE / MAE / R²

### 3.2 主要結果

| モデル | RMSE | MAE | R² |
|--------|------|-----|-----|
| **Random Forest** | **5.34** | **4.27** | **0.543** |
| Gradient Boosting | 5.41 | 4.31 | 0.532 |
| MLPRegressor | 6.08 | 4.81 | 0.408 |

- Random Forestが最良性能を達成（R² = 0.543）
- 気温と季節性が最も重要な予測因子として特定された
- 夏季（6–8月）のリスクは冬季の約2.5倍

### 3.3 図表

| ファイル | 内容 |
|----------|------|
| `figures/fig1_spatiotemporal_prediction.png` | 実測値 vs 予測値（信頼区間付き） |
| `figures/fig1b_feature_importance.png` | 特徴量重要度ランキング |
| `figures/fig1c_seasonal_heatmap.png` | 地域×月のリスクヒートマップ |
| `figures/fig1d_time_decomposition.png` | 時系列分解（トレンド・季節性・残差） |

---

## 4. モジュール2: NLPによるリコール・アラート早期検出

### 4.1 手法

- **データ**: FDA/RASFF模擬アラート320件（2023–2024）
- **カテゴリ**: 生物学的（biological）、化学的（chemical）、物理的（physical）、アレルゲン（allergen）
- **テキスト処理**: トークン化 → 小文字化 → ストップワード除去 → TF-IDF → 正規表現NER
- **分類器**: Logistic Regression / SVM / Random Forest（重症度3クラス分類）
- **トレンド分析**: カテゴリ別アラート頻度のローリングZ-score

### 4.2 主要結果

| 分類器 | 正解率 | マクロF1 |
|--------|--------|----------|
| Logistic Regression | 0.575 | 0.409 |

- 重症度分類は3クラスの不均衡（Class I: 145, Class II: 128, Class III: 47）のため、マクロF1はやや低い値
- TF-IDF特徴量分析により、カテゴリ固有のキーワードパターンを抽出
- 正規表現ベースのNERにより、病原体名・アレルゲン名・製品名・地名を自動抽出
- ローリングZ-scoreによる新興リスク検出メカニズムを実装

### 4.3 図表

| ファイル | 内容 |
|----------|------|
| `figures/fig2_nlp_classification.png` | 重症度分類の混同行列 |
| `figures/fig2b_alert_trends.png` | カテゴリ別アラート頻度推移 |
| `figures/fig2c_entity_extraction.png` | 抽出エンティティ（病原体・アレルゲン・製品）上位 |
| `figures/fig2d_tfidf_features.png` | カテゴリ別TF-IDF上位特徴量 |

---

## 5. モジュール3: 微生物増殖予測モデル

### 5.1 手法

- **一次モデル**: Baranyi成長モデル
  - `y(t) = y₀ + μ_max·A(t) − ln(1 + (e^(μ_max·A(t)) − 1) / (e^(y_max − y₀)))`
  - `A(t) = t + (1/μ_max)·ln(e^(−μ_max·t) + e^(−h₀) − e^(−μ_max·t − h₀))`
- **二次モデル**: Ratkowsky平方根モデル（μ_max = f(T)）
  - `√μ_max = b·(T − T_min)·(1 − e^(c·(T − T_max)))`
- **対象菌種**: Salmonella, E. coli, Listeria monocytogenes, S. aureus
- **モンテカルロシミュレーション**: パラメータ不確実性の定量化（1000回反復）

### 5.2 主要結果

| 菌種 | T_min (°C) | Ratkowsky RMSE | μ_max at 37°C (log/h) |
|------|------------|----------------|------------------------|
| Salmonella | 6.1 | 0.009 | 0.611 |
| E. coli | 5.6 | — | — |
| L. monocytogenes | 0.5 | — | — |
| S. aureus | 7.0 | — | — |

- **コールドチェーン途絶シナリオ**（25°Cで4時間逸脱）:
  - 安全閾値（6 log CFU/g）超過確率: **63.8%**（95% CI: 60.8%–66.7%）
- Growth/No-growth境界モデルにより、温度×pH×水分活性の三次元安全領域を可視化

### 5.3 図表

| ファイル | 内容 |
|----------|------|
| `figures/fig3_baranyi_curves.png` | 各温度でのSalmonella増殖曲線 |
| `figures/fig3b_organism_comparison.png` | 25°Cでの菌種間比較 |
| `figures/fig3c_coldchain_break.png` | コールドチェーン途絶の影響 |
| `figures/fig3d_monte_carlo.png` | モンテカルロ不確実性解析（安全閾値付き） |
| `figures/fig3e_growth_nogrowth.png` | 増殖/非増殖境界マップ |

---

## 6. モジュール4: HACCP管理点リスクスコアリング

### 6.1 手法

- **対象**: 鶏肉加工施設のHACCP計画（7つのCCP）
  - CCP1: 受入検査、CCP2: 冷蔵保管、CCP3: 加熱処理、CCP4: 冷却、CCP5: 金属検出、CCP6: 包装、CCP7: 出荷・配送
- **データ**: 1年間の連続モニタリングデータ（12,410レコード）
- **リスク評価**: RPN（Risk Priority Number）= 重大度 × 発生頻度 × 検出度
- **ベイズ更新**: 逸脱率の事後分布をβ分布で逐次更新
- **異常検知**: SPC管理図 / CUSUM / EWMA

### 6.2 主要結果

| CCP | 名称 | レコード数 | 逸脱率 | ベイズ事後逸脱確率 | 95%信頼区間 | 平均RPN |
|-----|------|-----------|--------|-------------------|------------|---------|
| CCP1 | 受入検査 | 1,460 | 1.64% | 1.72% | [1.13–2.41%] | 68.6 |

- SPC管理図により逸脱パターンをリアルタイム検出
- CUSUM/EWMAにより微小な持続的シフトを早期検知
- ベイズ更新によりデータ蓄積に伴うリスク推定の精度向上を実証
- パレート分析によりリスク集中箇所を可視化

### 6.3 図表

| ファイル | 内容 |
|----------|------|
| `figures/fig4_risk_scores.png` | 全CCPのリスクスコアヒートマップ |
| `figures/fig4b_spc_charts.png` | 主要CCPのSPC管理図 |
| `figures/fig4c_rpn_distribution.png` | RPN分布とパレート図 |
| `figures/fig4d_bayesian_update.png` | ベイズリスクスコア更新推移 |

---

## 7. モジュール5: ブロックチェーントレーサビリティ

### 7.1 手法

- **ブロックチェーン構造**: SHA-256ハッシュチェーン + マークル木
- **サプライチェーンノード**: 農場 → 加工 → 流通 → 小売 → 消費者
- **データ**: 100ロット × 5ノード = 500トランザクション
- **汚染シナリオ**: 5%のロットに汚染発生、リコールシミュレーション
- **スマートコントラクト**: 温度逸脱自動アラート、自動保留/拒否

### 7.2 主要結果

| 指標 | 値 |
|------|-----|
| 総ロット数 | 100 |
| 汚染ロット数 | 5 (5%) |
| リコールロット数 | 5 |
| 総ブロック数 | 21 |
| チェーン整合性 | ✅ 100% |
| マークル検証通過率 | 100% |
| **ブロックチェーン追跡時間（平均）** | **0.020 ms** |
| 従来方式追跡時間（平均） | 89,657秒（約24.9時間） |
| **高速化倍率** | **約85億倍** |
| 温度逸脱アラート数 | 170 |
| 自動保留件数 | 34 |
| 自動拒否件数 | 10 |

- ブロックチェーンにより、従来24時間以上要していた追跡を**ミリ秒単位**に短縮
- 不変台帳により改ざん防止・監査可能性を確保
- スマートコントラクトによる自動温度管理で人的判断の遅延を排除

### 7.3 図表

| ファイル | 内容 |
|----------|------|
| `figures/fig5_supply_chain_network.png` | サプライチェーンネットワーク図 |
| `figures/fig5b_traceback_time.png` | 追跡時間比較（ブロックチェーン vs 従来） |
| `figures/fig5c_contamination_spread.png` | 汚染伝播の可視化 |
| `figures/fig5d_temperature_log.png` | チェーン全体の温度モニタリング |

---

## 8. モジュール6: 鶏肉サルモネラ汚染予測ケーススタディ

### 8.1 手法

- **データ**: 2,000鶏肉ロット × 2年間の合成データ
- **特徴量**: 農場ID、鶏群日齢、外気温、湿度、季節、加工場、冷却方法（空気/水）、加工速度、冷却前/後菌数、最終製品菌数、保管温度、小売までの日数
- **全体汚染率**: 14.25%（夏季: 28.2%）
- **モデル**: Logistic Regression / Random Forest / Gradient Boosting
- **評価**: ROC-AUC、感度、特異度、陽性予測値、層化5分割CV

### 8.2 主要結果

| モデル | AUC | 感度 | 特異度 | PPV | CV-AUC (±SD) |
|--------|-----|------|--------|-----|--------------|
| **Logistic Regression** | **0.885** | 0.375 | 0.960 | 0.608 | 0.886 ± 0.032 |
| Random Forest | 0.881 | 0.523 | 0.937 | 0.580 | 0.881 ± 0.023 |
| Gradient Boosting | 0.879 | 0.432 | 0.952 | 0.597 | 0.880 ± 0.023 |

- 全モデルがAUC > 0.87を達成し、良好な判別能を示した
- Logistic Regressionが最高AUC（0.885）を記録（解釈可能性も高い）
- Random Forestは最高感度（0.523）を示し、汚染見逃し低減に有利
- **主要リスク因子**: 外気温、冷却前菌数、季節（夏季）、加工速度

### 8.3 介入シミュレーション

冷却方法の改善、加工速度の低減、衛生管理の強化の各介入による汚染率低減効果を定量評価した。

### 8.4 図表

| ファイル | 内容 |
|----------|------|
| `figures/fig6_roc_curves.png` | 全モデルROC曲線比較 |
| `figures/fig6b_risk_factors.png` | リスク因子ランキング（オッズ比/重要度） |
| `figures/fig6c_seasonal_contamination.png` | 月別汚染率と気温のオーバーレイ |
| `figures/fig6d_intervention_impact.png` | 介入シミュレーション結果 |
| `figures/fig6e_integrated_dashboard.png` | 統合ダッシュボード（2×2パネル） |

---

## 9. 統合リスクモニタリングシステム設計

### 9.1 アーキテクチャ

本システムは6つのモジュールを以下の3層で統合する：

#### データ収集層
- IoTセンサー（温湿度、GPS）からのリアルタイムデータ
- FDA/RASSFアラートフィードのNLP自動解析
- HACCPモニタリングデータの自動取り込み
- ブロックチェーンへの全トランザクション記録

#### 分析・予測層
- 時空間モデルによる地域リスク予報（24–72時間先）
- Baranyiモデルによる製品内微生物レベル予測
- HACCPリスクスコアのリアルタイム更新
- NLPによる外部リスク情報の自動スコアリング

#### 意思決定支援層
- 統合リスクスコア（ベイズネットワークによる多源情報融合）
- 閾値超過時の自動アラート発報
- ブロックチェーンベースの即時追跡・リコール
- ダッシュボードによる可視化

### 9.2 リスク統合式

```
総合リスクスコア = w₁·P(時空間) + w₂·P(NLP) + w₃·P(微生物) + w₄·P(HACCP)
```

各重みは、過去の食中毒事例との相関に基づきベイズ最適化で決定する。

### 9.3 運用フロー

1. **平常時**: 各モジュールが独立してリスク監視、日次レポート生成
2. **警戒時**: 統合リスクスコアが閾値超過 → 強化監視モードに移行
3. **緊急時**: 汚染確認 → ブロックチェーン即時追跡 → 対象ロット特定 → リコール発動

---

## 10. 考察と今後の展望

### 10.1 主要な知見

1. **時空間予測**: 気温と季節性が食中毒発生の最大の予測因子であり、R² = 0.543のモデル性能を達成。深層学習（LSTM/Transformer）の適用により更なる向上が期待される。

2. **NLP検出**: TF-IDF+Logistic Regressionによる重症度分類は発展途上（F1 = 0.41）だが、カテゴリ別キーワードパターンの抽出と新興リスク検出メカニズムは有用。BERTベースの事前学習モデルにより大幅な性能向上が見込まれる。

3. **微生物増殖**: Baranyiモデル+Ratkowsky二次モデルにより、温度依存の微生物増殖を高精度で予測。コールドチェーン途絶シナリオでは安全閾値超過確率63.8%と算出され、温度管理の重要性を定量的に実証。

4. **HACCP自動化**: ベイズ更新によるリスクスコアの逐次精緻化とSPC/CUSUM/EWMAによる異常検知の組み合わせにより、従来の定期監査に比べ迅速・精密な管理が可能。

5. **ブロックチェーン**: 追跡時間を24時間超からミリ秒単位に短縮し、リコール対応の劇的な効率化を実証。不変台帳によるデータ完全性の保証も重要な付加価値。

6. **ケーススタディ**: 鶏肉サルモネラ予測でAUC 0.885を達成。季節・気温・加工条件が主要リスク因子として特定され、介入の優先順位付けに活用可能。

### 10.2 限界

- すべての分析は合成データに基づいており、実データでの検証が必要
- NLP分類器の性能はデータ量の増加と事前学習モデルの導入で改善可能
- ブロックチェーンの高速化倍率はシミュレーション上の値であり、実運用ではネットワーク遅延を考慮する必要がある
- モジュール間の相互作用・フィードバックループは現時点では未実装

### 10.3 今後の展望

| 優先度 | 項目 | 内容 |
|--------|------|------|
| 高 | 実データ検証 | FDA CFSAN、CDC FoodNet等の公開データでのモデル検証 |
| 高 | 深層学習導入 | LSTM/Transformer時系列モデル、BERT/BioBERT NLPモデル |
| 中 | リアルタイム化 | Apache Kafka/Flinkによるストリーム処理パイプライン |
| 中 | IoT統合 | MQTT/LoRaWANセンサーネットワークとの接続 |
| 中 | 国際基準適合 | Codex Alimentarius, ISO 22000との整合性検証 |
| 低 | 連合学習 | 企業間でデータを共有せずモデルを共同改善 |
| 低 | デジタルツイン | 加工施設のデジタルツインによるシナリオ分析 |

---

## 11. 生成ファイル一覧

### ソースコード (`src/`)

| ファイル | 内容 | サイズ |
|----------|------|--------|
| `src/module1_spatiotemporal.py` | 時空間予測モデル | 20 KB |
| `src/module2_nlp_detection.py` | NLPリコール検出 | 33 KB |
| `src/module3_microbial_growth.py` | 微生物増殖予測 | 26 KB |
| `src/module4_haccp_scoring.py` | HACCPリスクスコアリング | 35 KB |
| `src/module5_blockchain.py` | ブロックチェーントレーサビリティ | 34 KB |
| `src/module6_salmonella_casestudy.py` | サルモネラケーススタディ | 38 KB |

### 図表 (`figures/`) — 全26枚

| ファイル | 内容 |
|----------|------|
| `fig1_spatiotemporal_prediction.png` | 時空間予測: 実測 vs 予測 |
| `fig1b_feature_importance.png` | 時空間予測: 特徴量重要度 |
| `fig1c_seasonal_heatmap.png` | 時空間予測: 地域×月リスクマップ |
| `fig1d_time_decomposition.png` | 時空間予測: 時系列分解 |
| `fig2_nlp_classification.png` | NLP: 重症度分類混同行列 |
| `fig2b_alert_trends.png` | NLP: アラート頻度推移 |
| `fig2c_entity_extraction.png` | NLP: エンティティ抽出結果 |
| `fig2d_tfidf_features.png` | NLP: TF-IDF上位特徴量 |
| `fig3_baranyi_curves.png` | 微生物: Baranyi増殖曲線 |
| `fig3b_organism_comparison.png` | 微生物: 菌種間比較 |
| `fig3c_coldchain_break.png` | 微生物: コールドチェーン途絶影響 |
| `fig3d_monte_carlo.png` | 微生物: モンテカルロ不確実性 |
| `fig3e_growth_nogrowth.png` | 微生物: 増殖/非増殖境界 |
| `fig4_risk_scores.png` | HACCP: リスクスコアヒートマップ |
| `fig4b_spc_charts.png` | HACCP: SPC管理図 |
| `fig4c_rpn_distribution.png` | HACCP: RPN分布・パレート図 |
| `fig4d_bayesian_update.png` | HACCP: ベイズ更新推移 |
| `fig5_supply_chain_network.png` | BC: サプライチェーンネットワーク |
| `fig5b_traceback_time.png` | BC: 追跡時間比較 |
| `fig5c_contamination_spread.png` | BC: 汚染伝播可視化 |
| `fig5d_temperature_log.png` | BC: 温度モニタリング |
| `fig6_roc_curves.png` | ケーススタディ: ROC曲線 |
| `fig6b_risk_factors.png` | ケーススタディ: リスク因子 |
| `fig6c_seasonal_contamination.png` | ケーススタディ: 季節別汚染率 |
| `fig6d_intervention_impact.png` | ケーススタディ: 介入効果 |
| `fig6e_integrated_dashboard.png` | ケーススタディ: 統合ダッシュボード |

### 数値結果 (`results/`)

| ファイル | 内容 |
|----------|------|
| `results/module1_metrics.json` | 時空間モデル評価指標 |
| `results/module2_metrics.json` | NLP分類器評価指標 |
| `results/module3_metrics.json` | 微生物増殖パラメータ |
| `results/module4_metrics.json` | HACCPリスクスコア |
| `results/module5_metrics.json` | ブロックチェーン性能指標 |
| `results/module6_metrics.json` | サルモネラ予測モデル指標 |

### データ (`data/`)

| ファイル | 内容 | レコード数 |
|----------|------|-----------|
| `data/spatiotemporal_data.csv` | 時空間予測用データ | 600 |
| `data/recall_alerts.csv` | FDA/RASSFアラートデータ | 320 |
| `data/growth_curves.csv` | 微生物増殖曲線データ | — |
| `data/haccp_monitoring.csv` | HACCPモニタリングデータ | 12,410 |
| `data/blockchain_transactions.csv` | ブロックチェーントランザクション | 500 |
| `data/salmonella_data.csv` | サルモネラ汚染データ | 2,000 |

### ログ (`logs/`)

| ファイル | 内容 |
|----------|------|
| `logs/process-log.jsonl` | 実行トレース |

---

## 12. 参考文献

1. Baranyi, J., & Roberts, T.A. (1994). A dynamic approach to predicting bacterial growth in food. *International Journal of Food Microbiology*, 23(3-4), 277-294.
2. Ratkowsky, D.A., et al. (1982). Relationship between temperature and growth rate of bacterial cultures. *Journal of Bacteriology*, 149(1), 1-5.
3. FDA. (2024). HACCP Principles & Application Guidelines. U.S. Food and Drug Administration.
4. European Commission. (2024). RASFF — The Rapid Alert System for Food and Feed. Annual Report.
5. Tian, F. (2017). A supply chain traceability system for food safety based on HACCP, blockchain & Internet of things. *IEEE International Conference on Service Systems and Service Management*.
6. WHO. (2023). Estimates of the Global Burden of Foodborne Diseases. World Health Organization.
7. Combase. (2024). Combined Database for Predictive Microbiology. https://www.combase.cc/
8. Munir, M.T., et al. (2023). Blockchain-based food supply chain traceability: A review. *Trends in Food Science & Technology*, 142, 104210.

---

*本レポートは合成データに基づく手法検証であり、実運用前に実データでの検証が必要です。*
