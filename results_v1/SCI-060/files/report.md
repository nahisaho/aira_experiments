# リアルワールドデータからのワクチン有効性（VE）推定：方法論フレームワーク

**DRAFT — NOT FOR DISTRIBUTION**  
**作成日**: 2026-05-23  
**バージョン**: 1.0  
**解析環境**: R 4.x（survival, gnm, mgcv, MatchIt, WeightIt, cobalt）

---

## 目次

1. [実験目的と背景](#1-実験目的と背景)
2. [方法論概要](#2-方法論概要)
3. [コンポーネント詳細](#3-コンポーネント詳細)
4. [主要な結果と数値](#4-主要な結果と数値)
5. [考察と今後の展望](#5-考察と今後の展望)
6. [生成ファイル一覧](#6-生成ファイル一覧)
7. [参考文献](#7-参考文献)

---

## 1. 実験目的と背景

### 目的

リアルワールドデータ（RWD）からワクチン有効性（Vaccine Effectiveness; VE）を推定するための包括的な統計的フレームワークを設計する。本フレームワークは以下の6つの方法論的課題に対応する：

1. **Test-Negative Design（TND）** の統計的性質と仮定の形式的検証
2. **経時的ワクチン効果減衰（waning）** の柔軟な推定モデル
3. **変異株特異的VE** 推定のための統計手法
4. **健康バイアス（healthy vaccinee bias）** の定量化と補正
5. **ブースター接種** の追加効果の因果推定
6. **mRNAワクチンの入院予防効果** 評価のケーススタディ

### 背景

COVID-19パンデミックにおいて、ランダム化比較試験（RCT）で示されたワクチンの効能（efficacy）が、実社会でどの程度維持されるかを評価することが不可欠となった。リアルワールドデータによるVE推定は、以下の理由から重要である：

- **変異株の出現**: 新規変異株に対するVEはRCTの対象外
- **効果の減衰**: 時間経過に伴うVEの低下はRCTの追跡期間では十分に評価できない
- **ブースター接種**: 追加接種の効果は市販後データでの評価が必要
- **特定集団**: 高齢者、基礎疾患者、医療従事者などのサブグループ解析

Test-Negative Design（TND）は、観察研究においてVEを推定するための最も広く使用されるデザインであり、医療受診バイアスをコントロールする利点を持つ（Jackson & Nelson, 2013; Fukushima & Hirota, 2017）。

---

## 2. 方法論概要

### 解析パイプライン全体構成

```
[01_simulate_data.R]  合成データ生成 (N=50,000)
        ↓
[02_tnd_analysis.R]   TND基本解析 + 仮定検証
        ↓
[03_waning_model.R]   VE減衰モデリング (4手法)
        ↓
[04_variant_specific_ve.R]  変異株別VE推定
        ↓
[05_healthy_vaccinee_bias.R]  バイアス補正 (IPTW/PS matching)
        ↓
[06_booster_causal.R]  ブースター因果推定 (AIPW)
        ↓
[07_hospitalization_case_study.R]  入院予防VEケーススタディ
```

### 使用パッケージと役割

| パッケージ | 用途 |
|-----------|------|
| `survival` | 条件付きロジスティック回帰 (`clogit`)、Cox比例ハザードモデル |
| `gnm` | 高次元層変数の効率的な条件付きロジスティック回帰 (`eliminate`) |
| `mgcv` | GAMによる柔軟なワニング曲線推定 |
| `splines` | 自然スプラインによる連続的VE変化 |
| `MatchIt` | 傾向スコアマッチング |
| `WeightIt` | 逆確率重み付け（IPTW） |
| `cobalt` | 共変量バランス診断（Love plot） |
| `sandwich` | ロバスト分散推定（Huber-White） |
| `EValue` | 未測定交絡に対するE-value |

---

## 3. コンポーネント詳細

### 3.1 Test-Negative Design (TND)

#### 統計的定式化

TNDでは、症状を有して検査を受けた者全員を対象とし、標的病原体（SARS-CoV-2）検査陽性者を**ケース**、陰性者を**コントロール**とする。

**モデル式（条件付きロジスティック回帰）**:

```
logit(P(Y=1 | X, S=s)) = α_s + β₁·Vax + β₂·Age + β₃·Sex + β₄·Comorbidity + β₅·PriorInf
```

- `Y`: 検査陽性（1=ケース, 0=コントロール）
- `Vax`: ワクチン接種状態
- `α_s`: カレンダー週 `s` の層別切片（`strata()` or `eliminate`）
- **VE = 1 − exp(β₁) = 1 − OR**

#### 実装方法

```r
# Method 1: survival::clogit
m_clogit <- clogit(
  tnd_case ~ vaccinated + age + sex + comorbidity + 
    prior_infection + strata(calendar_week),
  data = tnd_data)

# Method 2: gnm (大規模データ向け)
m_gnm <- gnm(
  tnd_case ~ vaccinated + age + sex + comorbidity + prior_infection,
  family = binomial(link = "logit"),
  eliminate = factor(calendar_week),
  data = tnd_data)
```

`gnm` の `eliminate` 引数は、多数の層（例: 52週 × 地域）がある場合に `clogit` より計算効率が大幅に改善される。

#### TND の仮定と検証方法

| 仮定 | 内容 | 検証方法 |
|------|------|---------|
| **A1**: 医療受診行動の独立性 | ワクチン接種状態と医療受診確率が独立 | ケース・コントロール間の共変量バランスチェック |
| **A2**: ワクチンの非標的病原体への無効果 | ワクチンが非標的病原体の感染確率に影響しない | **偽薬テスト** (falsification test) |
| **A3**: 検査感度・特異度の独立性 | 検査性能がワクチン接種状態に依存しない | 感度分析（検査種別で層別） |
| **A4**: 交絡のない比較 | 残余交絡がない | E-value計算、negative control outcome |

**偽薬テスト（Falsification Test）**:

```r
# 非標的病原体に対するVE → 0%が期待値
df_controls <- tnd_data %>% filter(target_positive == 0)
m_falsi <- glm(nontarget_positive ~ vaccinated + covariates,
  family = binomial, data = df_controls)
# 期待: OR ≈ 1.0 (VE ≈ 0%)
```

### 3.2 ワクチン効果減衰 (Waning) の推定

#### 4つのモデリングアプローチ

**Method 1: 区間別定数モデル (Piecewise)**

接種後の時間を離散区間に分割し、各区間のVEを独立に推定する。

```r
m_piece <- clogit(
  tnd_case ~ time_since_vax_cat + covariates + strata(calendar_week),
  data = df_vax)
```

時間区間: 0-13日, 14-59日, 60-119日, 120-179日, 180日以上

**Method 2: 自然スプラインモデル (Continuous)**

```r
m_spline <- glm(
  tnd_case ~ ns(days_since_vax, df = 4) + covariates,
  family = binomial, data = df_vax)
```

自由度4の自然スプラインにより、接種後の時間に対する滑らかなVE変化を推定。

**Method 3: GAMモデル (Data-Driven Smoothing)**

```r
m_gam <- gam(
  tnd_case ~ s(days_since_vax, bs = "cr", k = 10) + 
    s(age, bs = "cr", k = 5) + covariates,
  family = binomial, data = df_vax, method = "REML")
```

Cubic regression spline + REML による自動平滑化パラメータ選択。

**Method 4: 指数減衰モデル (Parametric Half-Life)**

```r
# VE(t) = VE₀ × exp(−λt)
# → logit(P) = β₀ + λ·t + ...
m_exp <- glm(tnd_case ~ days_since_vax + covariates,
  family = binomial, data = df_vax)
half_life <- log(2) / abs(coef(m_exp)["days_since_vax"])
```

**各手法の比較**:

| 手法 | 利点 | 欠点 | 推奨場面 |
|------|------|------|---------|
| Piecewise | 解釈が容易、CIが明確 | 区間境界が恣意的 | 政策報告 |
| Natural Spline | 滑らか、過学習を抑制 | dfの選択が必要 | 論文発表 |
| GAM | 自動平滑化、柔軟 | ブラックボックス的 | 探索的分析 |
| 指数減衰 | 半減期が計算可能 | 形状の仮定が必要 | メカニズム解釈 |

### 3.3 変異株特異的VE推定

#### 層別解析

各変異株の流行期間でデータを分割し、変異株別にTND解析を実施。

```r
for (v in variants) {
  df_v <- tnd_data %>% filter(variant == v)
  m_v <- glm(tnd_case ~ vaccinated + covariates,
    family = binomial, data = df_v)
  VE_v <- 1 - exp(coef(m_v)["vaccinated"])
}
```

#### 交互作用モデル

ワクチン効果が変異株によって異なるかを尤度比検定で評価:

```r
m_full <- glm(tnd_case ~ vaccinated * variant_f + covariates,
  family = binomial, data = tnd_data)
m_null <- glm(tnd_case ~ vaccinated + variant_f + covariates,
  family = binomial, data = tnd_data)
anova(m_null, m_full, test = "LRT")
```

#### 変異株特性との関連 (Meta-Regression)

変異株の免疫逃避スコアやRBD変異数とVEの関連を重み付き線形回帰で評価。

### 3.4 健康バイアス (Healthy Vaccinee Bias) の補正

#### バイアスのメカニズム

```
健康意識が高い → ワクチン接種 + 感染予防行動 → 見かけ上のVE過大評価
```

#### 3つの補正方法

**Method A: 逆確率重み付け (IPTW)**

```r
W <- weightit(vaccinated ~ age + sex + comorbidity + ...,
  data = tnd_data, method = "ps", estimand = "ATE")
m_iptw <- glm(tnd_case ~ vaccinated + covariates,
  family = binomial, data = tnd_data, weights = W$weights)
# + Huber-White robust SE
coeftest(m_iptw, vcov = vcovHC(m_iptw, type = "HC0"))
```

安定化重み（stabilized weights）を使用し、1/99パーセンタイルで切断。

**Method B: 傾向スコアマッチング**

```r
m_match <- matchit(vaccinated ~ covariates,
  data = tnd_data, method = "nearest", ratio = 1, caliper = 0.1)
matched_data <- match.data(m_match)
m_matched <- clogit(tnd_case ~ vaccinated + strata(subclass),
  data = matched_data)
```

`cobalt::love.plot()` で共変量バランスを可視化。

**Method C: Negative Control Outcome**

非標的病原体に対する見かけのVEをバイアスの指標として使用し、補正:

```r
OR_adjusted = OR_naive / OR_negative_control
VE_adjusted = 1 - OR_adjusted
```

**Method D: E-value (未測定交絡の評価)**

観察されたVEをnullにするために必要な未測定交絡の最小強度を計算。E-valueが大きいほど、結果がロバスト。

### 3.5 ブースター接種の因果推定

#### Target Trial Emulation フレームワーク

| 要素 | 内容 |
|------|------|
| 対象 | 初回接種シリーズ（2回）完了者 |
| 介入 | ブースター接種 vs 非ブースター |
| 割り付け | 観察データ + 逆確率重み付け |
| 追跡開始 | ブースター適格日 |
| アウトカム | COVID-19検査陽性（TND） |
| 因果対比 | Per-protocol効果 |

#### Augmented IPW (AIPW: 二重ロバスト推定)

```r
# 傾向スコアモデル (treatment model)
ps <- glm(booster ~ covariates, family = binomial)$fitted

# アウトカムモデル (outcome model)
m_out <- glm(tnd_case ~ booster + covariates, family = binomial)
mu1 <- predict(m_out, newdata = transform(df, booster = 1), type = "response")
mu0 <- predict(m_out, newdata = transform(df, booster = 0), type = "response")

# AIPW推定量
E[Y(1)] = mean(mu1 + booster/ps * (Y - mu1))
E[Y(0)] = mean(mu0 + (1-booster)/(1-ps) * (Y - mu0))
rVE = 1 - E[Y(1)] / E[Y(0)]
```

**二重ロバスト性**: 傾向スコアモデルまたはアウトカムモデルのいずれか一方が正しく特定されていれば、一致推定量を与える。

信頼区間はノンパラメトリックブートストラップ（B=500）で推定。

### 3.6 入院予防VEケーススタディ

#### 解析設計

- **対象**: 全検査受診者のうち入院した者（COVID入院 vs 非COVID入院のTND）
- **主要解析**: 全体VE、接種回数別VE、変異株別VE、年齢群別VE
- **ワニング**: GAMによる入院VE減衰曲線
- **感度分析**: 
  - 接種後14日未満の除外
  - 65歳以上に限定
  - 基礎疾患者に限定

---

## 4. 主要な結果と数値

### 4.1 シミュレーションデータの概要

| パラメータ | 値 |
|-----------|-----|
| 総対象者数 | 50,000 |
| ワクチン接種率 | 約58% |
| ケース（検査陽性） | 約12,500 (25%) |
| コントロール（検査陰性） | 約37,500 (75%) |
| 変異株分布 | Alpha/Delta/BA.1/BA.5 (各25%) |

### 4.2 TND基本解析結果

| モデル | 全体VE (%) | 95% CI |
|--------|-----------|--------|
| `clogit` (survival) | 推定: 45-55% | ― |
| `gnm` (eliminable) | 推定: 45-55% | ― |

**偽薬テスト**: 非標的病原体に対するVE ≈ 0%（仮定A2の支持）

### 4.3 Waning推定結果

| 接種後期間 | 推定VE (%) |
|-----------|-----------|
| 14-59日 | 最大（参照） |
| 60-119日 | 中程度の減衰 |
| 120-179日 | 顕著な減衰 |
| 180日以上 | 大幅な減衰 |

**推定半減期**: データ依存（指数減衰モデルの λ から計算）

### 4.4 変異株特異的VE

| 変異株 | 設定真値VE (%) | 推定VE (%) |
|--------|--------------|-----------|
| Alpha | 88 | ~85-90 |
| Delta | 82 | ~78-85 |
| Omicron BA.1 | 55 | ~50-60 |
| Omicron BA.5 | 45 | ~40-50 |

交互作用検定（ワクチン × 変異株）: p < 0.001（VEは変異株間で有意に異なる）

### 4.5 バイアス補正後VE

| 手法 | VE (%) |
|------|--------|
| 未調整（Naive） | やや過大 |
| IPTW | 補正後やや低下 |
| PS Matching | 補正後やや低下 |
| Negative Control Adjusted | 補正後やや低下 |

### 4.6 ブースター効果

| 推定法 | 相対VE (Booster vs Dose2) (%) | 95% CI |
|--------|------------------------------|--------|
| IPTW | 推定: 15-30 | ― |
| AIPW (二重ロバスト) | 推定: 15-30 | ブートストラップCI |

### 4.7 入院予防VE

| 解析 | VE (%) |
|------|--------|
| 全体（入院予防） | 推定: 70-90 |
| Dose 1 | 推定: 40-60 |
| Dose 2 | 推定: 70-85 |
| Booster | 推定: 85-95 |

> **注**: 上記の数値は合成データのシミュレーション設定に基づく期待値範囲です。実際の推定値は `Rscript run_all.R` を実行して得られます。

---

## 5. 考察と今後の展望

### 5.1 方法論的考察

#### TND の強みと限界

**強み**:
- 医療受診バイアスのコントロール（ケース・コントロール双方が検査を受けた者）
- 比較的実施が容易で、サーベイランスデータに適用可能
- カレンダー時間の交絡を層別で制御

**限界**:
- ワクチンが医療受診行動自体に影響する場合（例: 接種者が安心して受診を遅らせる）
- 検査アクセスの格差がある場合
- 繰り返し検査の扱い（個人内相関）

#### Waning推定の方法論的課題

- **時間起点の定義**: 2回目接種日 vs ブースター日
- **免疫カレンダー時間の交絡**: waning と変異株交代が同時進行
- **Depletion-of-susceptibles bias**: 感染者が追跡から脱落することによるVEの見かけの変化

#### 因果推定の課題

- **時間変動交絡**: ワクチン接種の意思決定と感染リスクが同時に変化
- **正値性の仮定**: 極端な傾向スコアを持つ個人の存在
- **Immortal time bias**: ブースター接種までの生存が必要

### 5.2 今後の展望

1. **Marginal Structural Models (MSMs)**: 時間変動交絡に対応した因果推定
2. **g-computation**: ブースター接種タイミングの最適化
3. **クラスタリング**: 地域・施設レベルのランダム効果を導入
4. **ゲノムデータ統合**: 全ゲノム配列データと連結し、系統ごとのVEを推定
5. **Hybrid immunity**: 自然感染 + ワクチンの複合免疫効果のモデリング
6. **相関免疫マーカー (Correlates of Protection)**: 中和抗体価とVEの関連性モデル
7. **Multi-state model**: 感染→入院→重症→死亡の段階的VE推定

### 5.3 実データ適用時の推奨事項

- **データ品質**: 接種記録と検査記録のリンケージ精度を確認
- **対象期間**: 変異株の流行期間とカレンダー時間を明確に定義
- **サンプルサイズ**: サブグループ解析（変異株 × 時間区間 × 年齢群）に十分な検体数を確保
- **報告**: STROBE-VE ガイドラインに準拠
- **再現性**: 乱数シードの固定、パッケージバージョンの記録

---

## 6. 生成ファイル一覧

### R スクリプト (`R/`)

| ファイル | 内容 | 行数 |
|---------|------|------|
| `R/00_setup.R` | パッケージ管理・環境設定 | ~25行 |
| `R/01_simulate_data.R` | 合成データ生成（N=50,000） | ~90行 |
| `R/02_tnd_analysis.R` | TND基本解析・偽薬テスト・ロバストSE | ~80行 |
| `R/03_waning_model.R` | 4手法のwaning推定（Piecewise/Spline/GAM/指数） | ~100行 |
| `R/04_variant_specific_ve.R` | 変異株別VE・交互作用検定 | ~80行 |
| `R/05_healthy_vaccinee_bias.R` | IPTW/PS Matching/NCO/E-value | ~180行 |
| `R/06_booster_causal.R` | Target Trial Emulation/IPTW/AIPW | ~140行 |
| `R/07_hospitalization_case_study.R` | 入院予防VE・サブグループ・感度分析 | ~180行 |
| `run_all.R` | 全パイプライン実行 | ~30行 |

### 図表 (`figures/`) — R実行時に生成

| ファイル | 内容 |
|---------|------|
| `figures/waning_curve.{png,svg}` | VE減衰曲線（GAM + Natural Spline） |
| `figures/ve_piecewise.png` | 時間区間別VE棒グラフ |
| `figures/ve_by_variant.{png,svg}` | 変異株別VEフォレストプロット |
| `figures/ve_vs_immune_escape.png` | VE vs 免疫逃避スコア散布図 |
| `figures/propensity_score_dist.png` | 傾向スコア分布 |
| `figures/bias_correction_comparison.png` | バイアス補正手法比較 |
| `figures/love_plot_balance.png` | 共変量バランスLoveプロット |
| `figures/booster_effectiveness.png` | ブースター効果推定 |
| `figures/hospitalization_case_study.{png,svg}` | 入院VE 4パネル図 |

### 結果テーブル (`results/`) — R実行時に生成

| ファイル | 内容 |
|---------|------|
| `results/tnd_ve_estimates.csv` | TND全体VE推定値 |
| `results/falsification_test.csv` | 偽薬テスト結果 |
| `results/waning_piecewise.csv` | 区間別VE |
| `results/waning_continuous.csv` | 連続waning曲線データ |
| `results/waning_halflife.csv` | 半減期推定 |
| `results/ve_by_variant.csv` | 変異株別VE |
| `results/bias_correction_comparison.csv` | バイアス補正比較 |
| `results/booster_ve_estimates.csv` | ブースターVE |
| `results/booster_ve_by_variant.csv` | 変異株別ブースター効果 |
| `results/hosp_ve_by_dose.csv` | 接種回数別入院VE |
| `results/hosp_ve_by_variant.csv` | 変異株別入院VE |
| `results/hosp_ve_by_age.csv` | 年齢群別入院VE |
| `results/hosp_sensitivity.csv` | 感度分析結果 |
| `results/hosp_waning_curve.csv` | 入院VE減衰曲線 |

### その他

| ファイル | 内容 |
|---------|------|
| `data/tnd_simulated.{rds,csv}` | シミュレーションデータ |
| `logs/process-log.jsonl` | 実行ログ |
| `report.md` | 本レポート |

---

## 7. 参考文献

1. Jackson ML, Nelson JC. The test-negative design for estimating influenza vaccine effectiveness. *Vaccine*. 2013;31(17):2165-2168.
2. Fukushima W, Hirota Y. Basic principles of test-negative design in evaluating influenza vaccine effectiveness. *Vaccine*. 2017;35(36):4796-4800.
3. Tartof SY, Slezak JM, Fischer H, et al. Effectiveness of mRNA BNT162b2 COVID-19 vaccine up to 6 months in a large integrated health system in the USA. *Lancet*. 2021;398(10309):1407-1416.
4. Andrews N, Tessier E, Stowe J, et al. Duration of protection against mild and severe disease by COVID-19 vaccines. *N Engl J Med*. 2022;386(4):340-350.
5. Lewnard JA, Patel MM, Jewell NP, et al. Theoretical framework for retrospective studies of the effectiveness of SARS-CoV-2 vaccines. *Epidemiology*. 2021;32(4):508-517.
6. Hernán MA, Robins JM. Using big data to emulate a target trial when a randomized trial is not available. *Am J Epidemiol*. 2016;183(8):758-764.
7. Robins JM, Rotnitzky A, Zhao LP. Estimation of regression coefficients when some regressors are not always observed. *J Am Stat Assoc*. 1994;89(427):846-866.
8. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. *Ann Intern Med*. 2017;167(4):268-274.
9. Sullivan SG, Tchetgen Tchetgen EJ, Cowling BJ. Theoretical basis of the test-negative study design for assessment of influenza vaccine effectiveness. *Am J Epidemiol*. 2016;184(5):345-353.
10. Thompson MG, Natarajan K, Irving SA, et al. Effectiveness of a third dose of mRNA vaccines against COVID-19–associated emergency department and urgent care encounters and hospitalizations. *MMWR*. 2022;71(4):139-145.

---

*本レポートは合成データに基づく方法論フレームワークの設計文書です。実データへの適用時には、倫理審査委員会の承認を得た上で、データ品質の検証と適切な感度分析を実施してください。*
