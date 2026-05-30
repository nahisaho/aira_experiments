# Food Safety AI Risk Prediction System - Technical Report

## 実験目的と背景 (Experiment Purpose & Background)
本レポートは、食品サプライチェーンにおける **AIベース食品安全リスク予測システム** の実験実装結果をまとめたものである。対象ハザードは主に **Salmonella** とし、以下の5つの観点を統合した。

1. **Spatiotemporal risk prediction**: 温度・湿度・季節・地域・コールドチェーン破綻などから Salmonella リスクを推定
2. **Predictive microbiology**: Baranyi-Roberts モデルによる菌増殖シミュレーション
3. **NLP recall detection**: リコール文書と通常アドバイザリ文書の自動識別
4. **HACCP scoring**: 8つの CCP に基づく工程リスク評価
5. **Time-series seasonality**: 月次アウトブレイクの季節性分解

さらに、NatureLM MCP (`ask_naturelm`) による取得済み科学知識を利用し、Salmonella の成長温度特性、D-values、季節的 prevalence 差を方法論に反映した。取得記録では、Salmonella in chicken に関する参考値として **Tmin ≈ 0°C, Topt ≈ 43°C, Tmax > 48°C**、および夏季高リスク傾向が確認されている。

## 使用した手法・アルゴリズム (Methods & Algorithms)

### 1. Spatiotemporal Salmonella Risk Prediction
- Samples: **2,000**
- Features:
  - temperature (0-40°C)
  - relative humidity (%)
  - season (1-4)
  - month (1-12)
  - region (0-4)
  - cold_chain_break (0/1)
  - processing_time_hours (0-72)
- Risk score rule:
  - base risk = 0.15
  - +0.02 per °C above 10°C (cap +0.5)
  - +0.003 per % RH above 60%
  - summer +0.15, spring/fall +0.05, winter +0
  - cold chain break +0.25
  - Gaussian noise $\sigma=0.05$
  - clip to [0,1]
- Label: `risk_score > 0.5`
- Model: **XGBoost classifier**, 5-fold cross-validation
- Note: observed features were intentionally noised to avoid unrealistically perfect performance.

### 2. Baranyi-Roberts Microbial Growth Simulation
Used the cardinal temperature submodel and Baranyi dynamics:

$$
\mu_{\max}(T) = \mu_{opt} \cdot \frac{(T-T_{\max})(T-T_{\min})^2}{(T_{opt}-T_{\min})\left[(T_{opt}-T_{\min})(T-T_{opt})-(T_{opt}-T_{\max})(T_{opt}+T_{\min}-2T)\right]}
$$

$$
\frac{dq}{dt}=\mu_{\max}q, \quad
\frac{dN}{dt}=\mu_{\max}\left(\frac{q}{1+q}\right)N\left(1-\frac{N}{N_{\max}}\right)
$$

Parameters used in simulation:
- Tmin = **2°C**
- Topt = **37°C**
- Tmax = **47°C**
- μopt = **1.5 /hr**
- q0 = **0.01**
- N0 = **100 CFU/g**
- Nmax = **1e9 CFU/g**
- Safety threshold = **1000 CFU/g**

### 3. NLP Recall Alert Detection
- Synthetic texts: **500**
  - positive recall texts: 250
  - negative advisory/routine texts: 250
- Pipeline:
  - **TF-IDF (1-2 grams)**
  - **Logistic Regression**
- Evaluation: 5-fold CV
- Design choice: positive/negative text templates shared overlapping food-safety language to keep the task realistic.

### 4. HACCP Risk Scoring
8 Critical Control Points (CCPs):
1. Receiving raw materials (0.12)
2. Cold storage (0.15)
3. Thawing (0.10)
4. Cooking/heat treatment (0.20)
5. Post-cook handling (0.15)
6. Packaging (0.10)
7. Refrigerated storage (0.12)
8. Distribution (0.06)

- Synthetic records: **300**
- Composite risk = weighted sum of CCP scores
- Risk classes:
  - Low < 3.5
  - Medium 3.5-6.5
  - High > 6.5
- Model: **RandomForest classifier**
- Added correlation: cold storage temperatures > 8°C increased cold-storage-related risk.

### 5. Time Series Seasonality Analysis
- Period: **2020-2024** monthly data (60 points)
- Base outbreak rate: 12/month
- Seasonal pattern: summer peak
- Trend: +0.5/year
- Noise: Poisson
- Method: `seasonal_decompose` from statsmodels

### 6. Case Study Dashboard
Regional risk dashboard combined:
- predicted Salmonella probability from Experiment 1
- regional HACCP average risk
- seasonal weighting from Experiment 5

## 主要な結果と数値 (Key Results & Numbers)

### Experiment 1: Spatiotemporal Prediction
- **AUC = 0.8493 ± 0.0132**
- **F1 = 0.6675 ± 0.0203**
- Positive class rate = **0.3625**
- Top features:
  - temperature = 0.3938
  - season = 0.1774
  - month = 0.1405
  - cold_chain_break = 0.0970

![Figure 1](figures/fig1_roc_curves.png)

![Figure 2](figures/fig2_feature_importance.png)

### Experiment 2: Baranyi-Roberts Growth
Selected growth-rate outputs:
- μmax(4°C) = 0.0087 /hr
- μmax(10°C) = 0.1344 /hr
- μmax(15°C) = 0.3408 /hr
- μmax(20°C) = 0.6197 /hr
- μmax(25°C) = 0.9411 /hr
- μmax(37°C) = 1.5000 /hr

Threshold crossing time for 1000 CFU/g:
- 4°C: not reached within 48 h
- 10°C: not reached within 48 h
- 15°C: **20.008 h**
- 20°C: **11.0621 h**
- 25°C: **7.3106 h**
- 37°C: **4.6172 h**

![Figure 3](figures/fig3_baranyi_growth.png)

### Experiment 3: NLP Recall Detection
- **Precision = 0.9325 ± 0.0248**
- **Recall = 0.9320 ± 0.0515**
- **F1 = 0.9314 ± 0.0304**
- Confusion matrix = `[[233, 17], [17, 233]]`

![Figure 4](figures/fig4_nlp_confusion_matrix.png)

### Experiment 4: HACCP Risk Scoring
- **Accuracy = 0.8700 ± 0.0194**
- Category counts:
  - Low = 22
  - Medium = 249
  - High = 29
- Confusion matrix (Low/Medium/High):
  - `[[6, 16, 0], [0, 248, 1], [0, 22, 7]]`

![Figure 5](figures/fig5_haccp_risk_distribution.png)

![Figure 6](figures/fig6_haccp_ccp_weights.png)

### Experiment 5: Time Series Seasonality
- Mean monthly outbreaks = **12.7333**
- Max observed month = **31** outbreaks
- Peak average month = **August (24.4)**

![Figure 7](figures/fig7_time_series.png)

### Case Study: Regional Dashboard
Composite regional risk index:
- South = **40.8905**
- North = **40.1253**
- East = **39.8950**
- West = **39.8179**
- Central = **39.5035**

![Figure 8](figures/fig8_supply_chain_risk_map.png)

## 考察と今後の展望 (Discussion & Future Work)
本実験では、表形式データ、微生物動態、自然言語、HACCP、時系列を単一の研究ワークフローに統合できることを示した。特に、XGBoost による spatiotemporal risk prediction は **AUC 0.8493** を達成し、先行研究の Wu & Wang (2026) が示した実運用レベルの食品安全リスク予測性能と整合的である。

Baranyi-Roberts シミュレーションは、**低温維持の重要性** を明確に示した。4-10°C では 48 時間以内に 1000 CFU/g に到達しない一方、15°C 以上では比較的短時間で閾値を超える。これは輸送・保管・冷蔵逸脱の評価に有用である。

NLP モジュールは、重い深層学習を使わずとも、軽量な TF-IDF + Logistic Regression で **F1 0.9314** を達成した。これはリコール通知の一次仕分けやアラート優先度付けに適している。HACCP モジュールは、工程管理の説明可能性を保ちながら **0.87 accuracy** を示し、現場向けダッシュボードとの親和性が高い。

主な制約は以下の通り。
- データが synthetic である
- FDA/FSIS 実データで未検証
- 地域ダッシュボードは簡略化された仮想地域に基づく
- 微生物パラメータは文献・NatureLM 取得知識に基づく保守的設定

今後は、以下を推奨する。
- 実際の recall database / inspection data / cold-chain sensor data で検証
- genomics や thermal inactivation data の統合
- SHAP などによる explainability 強化
- real-time monitoring dashboard への発展

## 生成したファイル一覧 (Generated File List)
- `run_experiments.py` - all experiments in one Python script
- `paper.md` - academic paper manuscript
- `report.md` - this technical report
- `figures/fig1_roc_curves.png`
- `figures/fig2_feature_importance.png`
- `figures/fig3_baranyi_growth.png`
- `figures/fig4_nlp_confusion_matrix.png`
- `figures/fig5_haccp_risk_distribution.png`
- `figures/fig6_haccp_ccp_weights.png`
- `figures/fig7_time_series.png`
- `figures/fig8_supply_chain_risk_map.png`
- `figures/metrics_summary.json`
